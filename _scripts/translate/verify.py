"""Translation verification — structural and semantic checks."""

from __future__ import annotations

import asyncio
import json
import re
from pathlib import Path

import anthropic

from .api import call_claude_async
from .config import (
    DEFAULT_PARALLEL,
    EN_DOCS,
    TranslationError,
    check_api_key,
    make_console,
)
from .models import TranslationFile
from .paths import en_to_lang_path, iter_en_docs
from .postprocess import (
    extract_code_blocks,
    extract_headers,
    iter_lines_outside_code,
)

# ---------------------------------------------------------------------------
# Semantic verification prompt
# ---------------------------------------------------------------------------

VERIFY_PROMPT = """You are verifying the quality of a translated document.
Given the first and last meaningful sentences from an English source document and its translation, check:
1. Does the translated first sentence correspond semantically to the English first sentence?
2. Is the translated first sentence complete (not truncated mid-sentence)?
3. Does the translated last sentence correspond semantically to the English last sentence?
4. Is the translated last sentence complete (not truncated mid-sentence)?

Reply ONLY with a JSON object (no markdown fencing):
{{"first_match": true/false, "first_complete": true/false, "last_match": true/false, "last_complete": true/false}}

English first sentence:
{en_first}

Translated first sentence:
{trans_first}

English last sentence:
{en_last}

Translated last sentence:
{trans_last}
"""

# Lines to skip when extracting meaningful content
_SKIP_LINE_RE = re.compile(
    r"^("
    r"#{1,6}\s"  # headings
    r"|(!{3}|\?{3}\+?)\s"  # admonitions
    r"|===\s"  # tabs
    r"|</?[a-zA-Z]"  # HTML tags
    r"|!\["  # images
    r"|\[.*\]:\s"  # link definitions
    r"|[-*_]{3,}\s*$"  # horizontal rules
    r")"
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _extract_meaningful_lines(lines: list[str]) -> list[str]:
    """Extract meaningful text lines, skipping frontmatter, code, headings, etc.

    Reuses ``iter_lines_outside_code`` from postprocess to avoid duplicating
    code-fence tracking logic.
    """
    meaningful: list[str] = []

    # Handle frontmatter: skip lines between opening --- and closing ---
    start_idx = 0
    if lines and lines[0].strip() == "---":
        for i, line in enumerate(lines[1:], 1):
            if line.strip() == "---":
                start_idx = i + 1
                break
        else:
            # Unclosed frontmatter — skip all
            return []

    remaining = lines[start_idx:]

    for _line_no, line in iter_lines_outside_code(remaining):
        stripped = line.strip()
        if not stripped or _SKIP_LINE_RE.match(stripped):
            continue
        meaningful.append(stripped)

    return meaningful


# ---------------------------------------------------------------------------
# Structural verification
# ---------------------------------------------------------------------------


def verify_translation_structural(en_path: Path, lang_path: Path) -> list[str]:
    """Check structural integrity of a translation. Returns list of issues (empty = OK)."""
    issues: list[str] = []

    if not lang_path.exists():
        issues.append("Translation file missing")
        return issues

    en_content = en_path.read_text(encoding="utf-8")
    trans_content = lang_path.read_text(encoding="utf-8")
    en_lines = en_content.splitlines()
    trans_lines = trans_content.splitlines()

    # Check for LLM wrapping entire output in a code fence
    first_line = trans_lines[0].strip() if trans_lines else ""
    if re.match(r"^```\w*$", first_line):
        issues.append(
            f"Translation appears wrapped in a code fence (starts with '{first_line}')"
        )

    # Header count and level agreement
    en_headers = extract_headers(en_lines)
    trans_headers = extract_headers(trans_lines)
    if len(en_headers) != len(trans_headers):
        issues.append(
            f"Header count mismatch: {len(trans_headers)} vs {len(en_headers)} in English"
        )
    else:
        for i, (eh, th) in enumerate(zip(en_headers, trans_headers)):
            if eh.level != th.level:
                issues.append(
                    f"Header level mismatch at header {i + 1}: "
                    f"{th.level} vs {eh.level} in English"
                )
                break

    # Code block count
    en_blocks = extract_code_blocks(en_lines)
    trans_blocks = extract_code_blocks(trans_lines)
    if len(en_blocks) != len(trans_blocks):
        issues.append(
            f"Code block count mismatch: {len(trans_blocks)} vs {len(en_blocks)} in English"
        )

    # Frontmatter presence
    en_has_fm = en_lines and en_lines[0].strip() == "---"
    trans_has_fm = trans_lines and trans_lines[0].strip() == "---"
    if en_has_fm and not trans_has_fm:
        issues.append("Missing frontmatter (English has it, translation does not)")

    # Admonition count (allow 30% variance)
    admonition_re = re.compile(r"^(!{3}|\?{3}\+?)\s+(\w[\w-]*)")
    en_admonitions = sum(1 for ln in en_lines if admonition_re.match(ln.strip()))
    trans_admonitions = sum(1 for ln in trans_lines if admonition_re.match(ln.strip()))
    if (
        en_admonitions
        and abs(en_admonitions - trans_admonitions) > en_admonitions * 0.3
    ):
        issues.append(
            f"Admonition count differs significantly: {trans_admonitions} vs {en_admonitions} in English"
        )

    # Admonition keywords should be English (not translated)
    valid_keywords = {
        "note",
        "tip",
        "warning",
        "info",
        "abstract",
        "success",
        "failure",
        "danger",
        "bug",
        "example",
        "quote",
        "question",
        "hint",
        "exercise",
        "solution",
        "full-code",
    }
    for ln in trans_lines:
        m = admonition_re.match(ln.strip())
        if m:
            keyword = m.group(2).lower()
            if keyword not in valid_keywords:
                issues.append(
                    f"Admonition keyword appears translated: '{m.group(2)}' "
                    f"(expected English keyword)"
                )
                break  # report once

    # Tab label consistency — all Before/After labels in a file should match
    tab_re = re.compile(r'^===\s+"([^"]+)"')
    tab_labels = set()
    for ln in trans_lines:
        m = tab_re.match(ln.strip())
        if m:
            tab_labels.add(m.group(1))
    # If we have translated Before/After, there should be at most 2 unique
    # labels for the pair (plus any custom labels). Check for known-bad
    # patterns: more than 4 unique tab labels often signals inconsistency.
    en_tab_labels = set()
    for ln in en_lines:
        m = tab_re.match(ln.strip())
        if m:
            en_tab_labels.add(m.group(1))
    if en_tab_labels and len(tab_labels) > len(en_tab_labels) + 2:
        issues.append(
            f"Tab labels may be inconsistent: {len(tab_labels)} unique labels "
            f"vs {len(en_tab_labels)} in English"
        )

    # Translation notice presence (skip for English source files)
    if lang_path != en_path:
        has_notice = "ai-translation-notice" in trans_content
        if not has_notice:
            issues.append("Missing AI translation notice span")

    return issues


# ---------------------------------------------------------------------------
# Semantic verification
# ---------------------------------------------------------------------------


async def verify_translation_semantic_async(
    en_path: Path,
    lang_path: Path,
    client: anthropic.AsyncAnthropic,
) -> list[str]:
    """Check first/last meaningful sentences match semantically. Returns issues."""
    if not lang_path.exists():
        return ["Translation file missing"]

    en_lines = en_path.read_text(encoding="utf-8").splitlines()
    trans_lines = lang_path.read_text(encoding="utf-8").splitlines()

    en_meaningful = _extract_meaningful_lines(en_lines)
    trans_meaningful = _extract_meaningful_lines(trans_lines)
    if not en_meaningful or not trans_meaningful:
        return []

    en_first, en_last = en_meaningful[0], en_meaningful[-1]
    trans_first, trans_last = trans_meaningful[0], trans_meaningful[-1]

    prompt = VERIFY_PROMPT.format(
        en_first=en_first,
        trans_first=trans_first,
        en_last=en_last,
        trans_last=trans_last,
    )

    result = await call_claude_async(prompt, f"verify:{lang_path.name}", client)

    # Parse JSON response
    text = result.text.strip()
    if text.startswith("```"):
        text = re.sub(r"^```\w*\n?", "", text)
        text = re.sub(r"\n?```$", "", text)

    try:
        data = json.loads(text)
    except json.JSONDecodeError as e:
        raise TranslationError(
            f"Semantic verification returned invalid JSON for {lang_path.name}: {e}\nResponse: {text!r}"
        )

    issues: list[str] = []
    if not data.get("first_match"):
        issues.append(f"First sentence does not match English meaning: {trans_first!r}")
    if not data.get("first_complete"):
        issues.append(f"First sentence appears truncated: {trans_first!r}")
    if not data.get("last_match"):
        issues.append(f"Last sentence does not match English meaning: {trans_last!r}")
    if not data.get("last_complete"):
        issues.append(f"Last sentence appears truncated: {trans_last!r}")
    return issues


# ---------------------------------------------------------------------------
# Combined verification
# ---------------------------------------------------------------------------


async def verify_translation_async(
    en_path: Path,
    lang_path: Path,
    client: anthropic.AsyncAnthropic,
) -> list[str]:
    """Full verification: structural + semantic. Returns list of issues."""
    issues = verify_translation_structural(en_path, lang_path)
    if lang_path.exists():
        semantic_issues = await verify_translation_semantic_async(
            en_path, lang_path, client
        )
        issues.extend(semantic_issues)
    return issues


# ---------------------------------------------------------------------------
# Broken-file scanner
# ---------------------------------------------------------------------------


async def get_broken_files(
    lang: str, parallel: int = DEFAULT_PARALLEL
) -> list[TranslationFile]:
    """Scan all translations for a language, returning those with issues.

    Runs full verification (structural + semantic) on every file in parallel.
    """
    console = make_console(stderr=True)
    check_api_key()
    client = anthropic.AsyncAnthropic()
    semaphore = asyncio.Semaphore(parallel)

    en_files = iter_en_docs()

    async def check_one(en_path: Path) -> TranslationFile | None:
        lang_path = en_to_lang_path(en_path, lang)
        async with semaphore:
            issues = await verify_translation_async(en_path, lang_path, client)
        if issues:
            rel = en_path.relative_to(EN_DOCS)
            console.print(f"  [red]Broken:[/red] {rel}")
            for issue in issues:
                console.print(f"    {issue}")
            return TranslationFile(
                en_path=en_path,
                lang_path=lang_path,
                language=lang,
                force_full=True,
            )
        return None

    console.print(f"[cyan]Scanning {len(en_files)} files for {lang}...[/cyan]")
    results = await asyncio.gather(
        *[check_one(p) for p in en_files], return_exceptions=True
    )
    broken: list[TranslationFile] = []
    for r in results:
        if isinstance(r, Exception):
            console.print(f"  [red]Error during scan: {r}[/red]")
        elif r is not None:
            broken.append(r)

    if broken:
        console.print(f"[yellow]Found {len(broken)} broken files for {lang}[/yellow]")
    else:
        console.print(f"[green]All files OK for {lang}[/green]")

    return broken
