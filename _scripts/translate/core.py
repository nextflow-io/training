"""Core translation orchestration.

Provides the main translation functions used by both the CLI and CI commands.
"""

from __future__ import annotations

import asyncio
import hashlib
import time
from datetime import datetime

import anthropic

from .api import call_claude_async
from .chunking import Chunk, diff_chunks, reassemble_chunks, split_into_chunks
from .config import MAX_VERIFY_RETRIES, make_console
from .git_utils import (
    get_file_at_commit,
    get_file_diff,
    get_prompt_diff,
    get_translation_baseline,
)
from .models import FileLogEntry, TranslationFile, TranslationLog
from .paths import get_languages
from .postprocess import post_process_file
from .progress import TranslationProgress, progress_logger
from .prompts import build_chunk_translation_prompt, build_translation_prompt
from .verify import verify_translation_async

# Minimum number of ## chunks to attempt chunk-based translation
_MIN_CHUNKS_FOR_CHUNKING = 3
# If more than this fraction of chunks changed, fall back to full-file
_MAX_CHANGED_CHUNK_RATIO = 0.5


async def _translate_chunks_async(
    tf: TranslationFile,
    en_content: str,
    existing_content: str,
    baseline: str,
    client: anthropic.AsyncAnthropic,
    entry: FileLogEntry,
) -> str | None:
    """Attempt chunk-based translation.  Returns reassembled text or None on fallback.

    Returns None if chunking is not applicable or fails, signalling the caller
    to fall back to full-file translation.
    """
    console = make_console()
    filename = str(tf.relative_path)
    langs = get_languages()

    # Get old English content at baseline
    old_en_content = get_file_at_commit(tf.en_path, baseline)
    if old_en_content is None:
        return None

    # Split into chunks
    old_en_chunks = split_into_chunks(old_en_content)
    new_en_chunks = split_into_chunks(en_content)
    trans_chunks = split_into_chunks(existing_content)

    # Check feasibility
    if len(new_en_chunks) < _MIN_CHUNKS_FOR_CHUNKING:
        return None

    # Verify translation chunk count matches old English (they should have
    # identical heading structure from previous verification)
    if len(trans_chunks) != len(old_en_chunks):
        console.print(
            f"  [{filename}] Chunk count mismatch: translation has "
            f"{len(trans_chunks)} vs old English {len(old_en_chunks)}, "
            f"falling back to full-file"
        )
        return None

    # Diff old vs new English chunks
    cdiff = diff_chunks(old_en_chunks, new_en_chunks)

    # Check if too many chunks changed
    num_changed = len(cdiff.modified) + len(cdiff.added)
    if num_changed == 0:
        # No chunks changed — file may have changed only in whitespace or
        # non-heading structural elements. Fall back to full-file.
        return None
    if num_changed > len(new_en_chunks) * _MAX_CHANGED_CHUNK_RATIO:
        console.print(
            f"  [{filename}] {num_changed}/{len(new_en_chunks)} chunks changed, "
            f"falling back to full-file"
        )
        return None

    trans_by_old_idx = {c.index: c for c in trans_chunks}
    result_chunks: list[Chunk] = []

    chunk_plan: dict[int, tuple[str, int | None]] = {}
    for old_idx, new_idx in cdiff.unchanged:
        chunk_plan[new_idx] = ("unchanged", old_idx)
    for old_idx, new_idx in cdiff.modified:
        chunk_plan[new_idx] = ("modified", old_idx)
    for new_idx in cdiff.added:
        chunk_plan[new_idx] = ("added", None)

    entry.chunked = True
    entry.total_chunks = len(new_en_chunks)
    entry.unchanged_chunks = len(cdiff.unchanged)

    console.print(
        f"  [{filename}] Chunk mode: {len(cdiff.unchanged)} unchanged, "
        f"{len(cdiff.modified)} modified, {len(cdiff.added)} added, "
        f"{len(cdiff.removed)} removed"
    )

    chunks_to_translate: list[tuple[int, str, str | None]] = []

    for ni in range(len(new_en_chunks)):
        plan = chunk_plan.get(ni)
        if plan is None:
            # Shouldn't happen, but fall back to translating it
            chunks_to_translate.append((ni, new_en_chunks[ni].content, None))
            continue

        kind, old_idx = plan
        if kind == "unchanged":
            result_chunks.append(
                Chunk(
                    heading=new_en_chunks[ni].heading,
                    content=trans_by_old_idx[old_idx].content,
                    index=ni,
                )
            )
        elif kind == "modified":
            existing_chunk = trans_by_old_idx[old_idx].content
            chunks_to_translate.append((ni, new_en_chunks[ni].content, existing_chunk))
        elif kind == "added":
            chunks_to_translate.append((ni, new_en_chunks[ni].content, None))

    entry.translated_chunks = len(chunks_to_translate)

    translated: dict[int, str] = {}
    for ni, en_chunk, existing_chunk in chunks_to_translate:
        prev_heading = new_en_chunks[ni - 1].heading if ni > 0 else None
        next_heading = (
            new_en_chunks[ni + 1].heading if ni + 1 < len(new_en_chunks) else None
        )

        prompt = build_chunk_translation_prompt(
            tf.language,
            langs[tf.language],
            en_chunk,
            existing_chunk,
            prev_heading,
            next_heading,
        )
        result = await call_claude_async(prompt, f"{filename}[chunk {ni}]", client)
        translated[ni] = result.text.strip()

        entry.input_tokens += result.input_tokens
        entry.output_tokens += result.output_tokens
        entry.model = result.model
        entry.stop_reason = result.stop_reason

    for ni, en_chunk, _ in chunks_to_translate:
        result_chunks.append(
            Chunk(
                heading=new_en_chunks[ni].heading,
                content=translated[ni],
                index=ni,
            )
        )

    result_chunks.sort(key=lambda c: c.index)
    return reassemble_chunks(result_chunks)


async def translate_file_async(
    tf: TranslationFile,
    semaphore: asyncio.Semaphore,
    progress: TranslationProgress,
    client: anthropic.AsyncAnthropic,
    log: TranslationLog | None = None,
) -> None:
    """Translate a single file with verification and retry.

    On the first attempt, tries chunk-based translation if applicable.
    Falls back to full-file translation if chunking is not feasible or
    if verification fails.
    """
    console = make_console()
    filename = str(tf.relative_path)
    entry = FileLogEntry(filename=filename, started_at=datetime.now().isoformat())
    start_time = time.time()

    async with semaphore:
        await progress.start_one(filename)
        try:
            langs = get_languages()
            en_content = tf.en_path.read_text(encoding="utf-8")
            original_existing = (
                tf.lang_path.read_text(encoding="utf-8") if tf.exists else None
            )

            entry.input_lines = en_content.count("\n") + 1
            entry.input_hash = hashlib.md5(en_content.encode()).hexdigest()[:12]

            # Resolve baseline once (lru_cached, but avoids duplicate try/except)
            try:
                baseline = get_translation_baseline()
            except Exception:
                baseline = None

            for attempt in range(1 + MAX_VERIFY_RETRIES):
                output_content = None

                # First attempt: try chunk-based translation
                if (
                    attempt == 0
                    and not tf.force_full
                    and not tf.prompt_changed
                    and original_existing is not None
                    and baseline
                ):
                    output_content = await _translate_chunks_async(
                        tf, en_content, original_existing, baseline, client, entry
                    )

                # Full-file translation (initial or fallback)
                if output_content is None:
                    if tf.force_full or attempt > 0:
                        use_existing = None
                        en_diff = None
                        prompt_diff_text = None
                    else:
                        use_existing = original_existing
                        en_diff = None
                        prompt_diff_text = None
                        if use_existing and baseline:
                            if tf.prompt_changed:
                                prompt_diff_text = get_prompt_diff(
                                    tf.language, baseline
                                )
                            else:
                                en_diff = get_file_diff(tf.en_path, baseline)

                    prompt = build_translation_prompt(
                        tf.language,
                        langs[tf.language],
                        en_content,
                        use_existing,
                        en_diff,
                        prompt_diff_text,
                    )
                    result = await call_claude_async(prompt, filename, client)
                    output_content = f"{result.text}\n"

                    entry.model = result.model
                    entry.input_tokens += result.input_tokens
                    entry.output_tokens += result.output_tokens
                    entry.stop_reason = result.stop_reason
                    entry.continuations = result.continuations

                tf.lang_path.parent.mkdir(parents=True, exist_ok=True)
                tf.lang_path.write_text(
                    output_content if output_content.endswith("\n") else output_content + "\n",
                    encoding="utf-8",
                    newline="\n",
                )
                post_process_file(tf.lang_path, tf.language)

                issues = await verify_translation_async(
                    tf.en_path, tf.lang_path, client
                )

                if not issues:
                    break

                if attempt < MAX_VERIFY_RETRIES:
                    console.print(
                        f"[yellow]Verification failed for {filename} "
                        f"(attempt {attempt + 1}/{1 + MAX_VERIFY_RETRIES}), "
                        f"retrying...[/yellow]"
                    )
                    for issue in issues:
                        console.print(f"  [dim]{issue}[/dim]")
                else:
                    entry.error = (
                        f"Verification failed after {MAX_VERIFY_RETRIES} retries: "
                        + "; ".join(issues)
                    )
                    console.print(
                        f"[red]Verification failed for {filename} "
                        f"after {MAX_VERIFY_RETRIES} retries:[/red]"
                    )
                    for issue in issues:
                        console.print(f"  [red]{issue}[/red]")

            # Compute hash and changed flag AFTER post-processing
            final_content = tf.lang_path.read_text(encoding="utf-8")
            entry.output_lines = final_content.count("\n")
            entry.output_hash = hashlib.md5(final_content.encode()).hexdigest()[:12]
            entry.changed = (
                original_existing is None or original_existing != final_content
            )

            await progress.finish_one(filename, success=True)
        except Exception as e:
            entry.error = str(e)
            await progress.finish_one(filename, success=False)
            raise
        finally:
            entry.finished_at = datetime.now().isoformat()
            entry.duration_s = time.time() - start_time
            if log:
                await log.add_entry(entry)


async def translate_all(
    files: list[TranslationFile],
    parallel: int,
    log: TranslationLog | None = None,
) -> None:
    """Translate all files in parallel with progress logging."""
    console = make_console()

    # Pre-compute line counts (used for sorting and progress tracking)
    file_lines = {
        str(f.relative_path): f.en_path.read_text(encoding="utf-8").count("\n")
        for f in files
    }

    # Sort largest files first to avoid "long pole" problem
    files.sort(key=lambda f: file_lines[str(f.relative_path)], reverse=True)

    semaphore = asyncio.Semaphore(parallel)
    progress = TranslationProgress(files, file_lines)
    client = anthropic.AsyncAnthropic()

    tasks = [translate_file_async(tf, semaphore, progress, client, log) for tf in files]

    # Run translations with progress logger
    logger_task = asyncio.create_task(progress_logger(progress))
    try:
        results = await asyncio.gather(*tasks, return_exceptions=True)
    finally:
        logger_task.cancel()
        try:
            await logger_task
        except asyncio.CancelledError:
            pass

    if log:
        log.write()
        console.print(f"[dim]Log written to: {log.path}[/dim]")

    console.print(f"[bold green]Translation complete:[/bold green] {progress.status()}")

    # Raise first error if any failed
    errors = [r for r in results if isinstance(r, Exception)]
    if errors:
        for e in errors:
            console.print(f"[red]Error: {e}[/red]")
        raise errors[0]
