"""Translation verification — structural and semantic checks."""

from __future__ import annotations

import json
import re
from collections import Counter
from pathlib import Path

import anthropic

from .api import call_claude_async
from .config import (
    DEFAULT_PARALLEL,
    EN_DOCS,
    TranslationError,
    make_console,
)
from .models import TranslationFile
from .paths import en_to_lang_path, iter_en_docs
from .postprocess import (
    FENCE_LABEL_RE,
    IMAGE_RE,
    MARKDOWN_LINK_RE,
    TAG_RE,
    URL_RE,
    parse_frontmatter,
    extract_code_blocks,
    extract_directives,
    extract_headers,
    extract_imports,
    iter_lines_outside_code,
    iter_tags,
    mask_inline_code,
    strip_translation_notice,
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
    r"|:{3,}"  # asides
    r"|import\s"  # code imports
    r"|\{/\*"  # comments
    r"|</?[a-zA-Z]"  # HTML/JSX tags
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

_ESCAPES_RE = re.compile(r"\\[{}<]|&lt;|\{/\*.*?\*/\}")
_HEADING_ID_RE = re.compile(r"\s*\{#[^}]+\}\s*$")


def _syntax_counts(lines: list[str]) -> Counter:
    """Counts of the MDX syntax that must match English: tags, asides, bare braces and `<`."""
    if lines and lines[0].strip() == "---":
        end = next((i for i in range(1, len(lines)) if lines[i].strip() == "---"), 0)
        lines = lines[end + 1 :]
    counts: Counter = Counter()
    for _, m in iter_tags(lines):
        counts[f"<{m['close']}{m['name']}{m['self']}>"] += 1
    for _, m in extract_directives(lines):
        counts[f":::{m.group(3) or ' (close)'}"] += 1
    for _, line in iter_lines_outside_code(lines):
        rest = _ESCAPES_RE.sub("", TAG_RE.sub("", mask_inline_code(line)))
        rest = _HEADING_ID_RE.sub("", rest) if line.startswith("#") else rest
        for ch in "{}<":
            counts[f"bare {ch}"] += rest.count(ch)
    return counts


def _link_targets(lines: list[str]) -> list[str]:
    """Link and ``href`` targets without their ``#anchor``, which translations re-point."""
    targets = []
    for _, line in iter_lines_outside_code(lines):
        targets += [m["url"].split("#")[0] for m in MARKDOWN_LINK_RE.finditer(line)]
        targets += [m[1].split("#")[0] for m in re.finditer(r'\bhref="([^"]*)"', line)]
    return targets


def _asset_targets(lines: list[str], page_dir: Path) -> list[Path]:
    """Files that relative images and ``?raw`` imports resolve to."""
    paths = [m["src"] for _, l in iter_lines_outside_code(lines) for m in IMAGE_RE.finditer(l)]
    paths += [m["path"].split("?")[0] for _, m in extract_imports(lines)]
    return [(page_dir / p).resolve() for p in paths if not URL_RE.match(p)]


def _key_paths(data: dict, prefix: str = "") -> set[str]:
    keys = set()
    for k, v in data.items():
        keys.add(prefix + str(k))
        if isinstance(v, dict):
            keys |= _key_paths(v, f"{prefix}{k}.")
    return keys


def verify_translation_structural(en_path: Path, lang_path: Path) -> list[str]:
    """Check structural integrity of a translation. Returns list of issues (empty = OK)."""
    issues: list[str] = []

    if not lang_path.exists():
        issues.append("Translation file missing")
        return issues

    en_content = en_path.read_text(encoding="utf-8")
    trans_content = lang_path.read_text(encoding="utf-8")
    en_lines = en_content.splitlines()
    trans_lines = strip_translation_notice(trans_content.splitlines())

    # Check for LLM wrapping entire output in a code fence
    first_line = trans_lines[0].strip() if trans_lines else ""
    if re.match(r"^```\w*$", first_line):
        issues.append(
            f"Translation appears wrapped in a code fence (starts with '{first_line}')"
        )

    # Frontmatter: present, valid YAML, same keys
    en_fm, _ = parse_frontmatter(en_content)
    trans_fm, _ = parse_frontmatter(trans_content)
    if en_fm is not None:
        if trans_fm is None:
            issues.append("Missing or invalid frontmatter")
        elif _key_paths(en_fm) != _key_paths(trans_fm):
            diff = sorted(_key_paths(en_fm) ^ _key_paths(trans_fm))
            issues.append(f"Frontmatter keys differ from English: {', '.join(diff)}")

    # Header count, levels and pinned ids
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
        for eh, th in zip(en_headers, trans_headers):
            # A missing id is tolerated (post-processing adds it); a changed one is not.
            if th.anchor and th.anchor.strip() != eh.anchor.strip():
                issues.append(
                    f"Heading id changed: {th.anchor.strip()} vs {eh.anchor.strip() or 'none'} in English"
                )

    # Code blocks: same count and info strings (step/label values are translated)
    en_blocks = extract_code_blocks(en_lines)
    trans_blocks = extract_code_blocks(trans_lines)
    if len(en_blocks) != len(trans_blocks):
        issues.append(
            f"Code block count mismatch: {len(trans_blocks)} vs {len(en_blocks)} in English"
        )
    else:
        for eb, tb in zip(en_blocks, trans_blocks):
            info = lambda b: FENCE_LABEL_RE.sub(r'\1=""', b.lines[0].strip())  # noqa: E731
            if info(eb) != info(tb):
                issues.append(
                    f"Code block info string changed at line {tb.start_line + 1}: "
                    f"{tb.lines[0].strip()!r} vs {eb.lines[0].strip()!r} in English"
                )

    # MDX syntax: JSX/HTML tags, ::: asides, unescaped { } <
    en_counts = _syntax_counts(en_lines)
    trans_counts = _syntax_counts(trans_lines)
    for key in sorted(en_counts.keys() | trans_counts.keys()):
        if en_counts[key] != trans_counts[key]:
            issues.append(
                f"MDX syntax count differs for {key}: "
                f"{trans_counts[key]} vs {en_counts[key]} in English"
            )

    # Link targets, images and code imports
    if _link_targets(en_lines) != _link_targets(trans_lines):
        issues.append("Link targets differ from English")
    if _asset_targets(en_lines, en_path.parent) != _asset_targets(
        trans_lines, lang_path.parent
    ):
        issues.append("Image or import paths do not point at the English files")

    # Translation notice presence (span on pages, aside on the homepage)
    if lang_path != en_path and not (
        "ai-translation-notice" in trans_content or "TRANSLATING.md" in trans_content
    ):
        issues.append("Missing AI translation notice")

    # Truncation check: translation should be roughly comparable in length
    if len(trans_lines) < len(en_lines) * 0.7:
        issues.append(
            f"Translation may be truncated: {len(trans_lines)} lines "
            f"vs {len(en_lines)} in English"
        )

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
    """Full verification: structural checks only.

    Semantic verification (LLM-based first/last sentence check) was removed
    from the hot path because it produced false positives that triggered
    expensive retries.  Truncation is now caught by the line-count ratio
    check in ``verify_translation_structural``, combined with the API
    continuation loop in ``api.py`` and header/code-block count checks.

    The ``verify_translation_semantic_async`` function is still available
    for manual or debug use.
    """
    return verify_translation_structural(en_path, lang_path)


# ---------------------------------------------------------------------------
# Broken-file scanner
# ---------------------------------------------------------------------------


async def get_broken_files(
    lang: str, parallel: int = DEFAULT_PARALLEL
) -> list[TranslationFile]:
    """Scan all translations for a language, returning those with issues.

    Runs structural verification on every file.  The ``parallel`` parameter
    is kept for API compatibility but is no longer used since verification
    is now purely structural (no LLM calls).
    """
    console = make_console(stderr=True)
    en_files = iter_en_docs()
    broken: list[TranslationFile] = []

    console.print(f"[cyan]Scanning {len(en_files)} files for {lang}...[/cyan]")
    for en_path in en_files:
        lang_path = en_to_lang_path(en_path, lang)
        issues = verify_translation_structural(en_path, lang_path)
        if issues:
            rel = en_path.relative_to(EN_DOCS)
            console.print(f"  [red]Broken:[/red] {rel}")
            for issue in issues:
                console.print(f"    {issue}")
            broken.append(
                TranslationFile(
                    en_path=en_path,
                    lang_path=lang_path,
                    language=lang,
                    force_full=True,
                )
            )

    if broken:
        console.print(f"[yellow]Found {len(broken)} broken files for {lang}[/yellow]")
    else:
        console.print(f"[green]All files OK for {lang}[/green]")

    return broken
