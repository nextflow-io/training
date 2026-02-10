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
from .config import MAX_VERIFY_RETRIES, make_console
from .git_utils import get_file_diff, get_translation_baseline
from .models import FileLogEntry, TranslationFile, TranslationLog
from .paths import get_languages
from .postprocess import post_process_file
from .progress import TranslationProgress, progress_logger
from .prompts import build_translation_prompt
from .verify import verify_translation_async


async def translate_file_async(
    tf: TranslationFile,
    semaphore: asyncio.Semaphore,
    progress: TranslationProgress,
    client: anthropic.AsyncAnthropic,
    log: TranslationLog | None = None,
) -> None:
    """Translate a single file with verification and retry.

    After translation and post-processing, verifies the output for structural
    and semantic correctness. If verification fails, retries up to
    MAX_VERIFY_RETRIES times with a full re-translation (no incremental).
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

            for attempt in range(1 + MAX_VERIFY_RETRIES):
                # On first attempt, use incremental if available (unless force_full).
                # On retries, force full re-translation without existing reference.
                if tf.force_full or attempt > 0:
                    use_existing = None
                    en_diff = None
                else:
                    use_existing = original_existing
                    en_diff = None
                    if use_existing and not tf.prompt_changed:
                        baseline = get_translation_baseline()
                        if baseline:
                            en_diff = get_file_diff(tf.en_path, baseline)

                prompt = build_translation_prompt(
                    tf.language, langs[tf.language], en_content, use_existing, en_diff
                )
                result = await call_claude_async(prompt, filename, client)
                output_content = f"{result.text}\n"

                # Log API response metadata (overwritten on retry)
                entry.model = result.model
                entry.input_tokens += result.input_tokens
                entry.output_tokens += result.output_tokens
                entry.stop_reason = result.stop_reason
                entry.continuations = result.continuations

                tf.lang_path.parent.mkdir(parents=True, exist_ok=True)
                tf.lang_path.write_text(output_content, encoding="utf-8", newline="\n")
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
