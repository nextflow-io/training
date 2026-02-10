"""Post-processing: fix structural elements in translated markdown.

After Claude translates a file, this module ensures that headers, links,
HTML links, code blocks, admonitions, tab labels, frontmatter, and the
AI translation notice match the English source structurally. Code blocks
are restored from source with only translated comments preserved.
Glossary-controlled elements (admonitions, tab labels, notices) are
enforced deterministically from per-language YAML glossary files.
"""

from __future__ import annotations

import re
from functools import lru_cache
from pathlib import Path

import yaml

from .config import (
    DOCS_ROOT,
    HASH_COMMENT_LANGS,
    MIXED_COMMENT_LANGS,
    REPO_ROOT,
    SLASH_COMMENT_LANGS,
    StructureMismatchError,
    make_console,
)
from .models import CodeBlock, Header, HtmlLink, Link

# ---------------------------------------------------------------------------
# Compiled regex patterns
# ---------------------------------------------------------------------------

HEADER_RE = re.compile(r"^(#{1,6})\s+(.+?)(\s*\{#[^}]+\})?\s*$")
MARKDOWN_LINK_RE = re.compile(
    r"(?<!!)\[(?P<text>[^\]]*)\]\((?P<url>[^)\s]+)(?:\s+\"(?P<title>[^\"]*)\")?\)(?:\{(?P<attrs>[^}]*)\})?"
)
HTML_LINK_RE = re.compile(r"<a\s+([^>]*)>(.*?)</a>", re.DOTALL)
CODE_FENCE_RE = re.compile(r"^(`{3,})([\w-]*)")
HASH_COMMENT_RE = re.compile(r"^(.*?)((?:\s+#|^#)\s.*)$")
SLASH_COMMENT_RE = re.compile(r"^(.*?)((?:\s+//|^//)\s.*)$")
BLOCK_COMMENT_START_RE = re.compile(r"^\s*/\*")
BLOCK_COMMENT_END_RE = re.compile(r"\*/\s*$")


# ---------------------------------------------------------------------------
# Iterators
# ---------------------------------------------------------------------------


def iter_lines_outside_code(lines: list[str]):
    """Yield ``(line_no, line)`` for lines outside code blocks.

    Public so that the verification module can reuse it instead of
    reimplementing code-fence tracking.
    """
    in_code = False
    fence = ""
    for i, line in enumerate(lines):
        stripped = line.lstrip()
        if not in_code:
            if m := CODE_FENCE_RE.match(stripped):
                in_code, fence = True, m.group(1)
                continue
        elif stripped.startswith(fence) and stripped.strip() == fence:
            in_code, fence = False, ""
            continue
        if not in_code:
            yield i, line


# ---------------------------------------------------------------------------
# Extractors
# ---------------------------------------------------------------------------


def extract_headers(lines: list[str]) -> list[Header]:
    """Extract headers outside code blocks."""
    return [
        Header(i, m.group(1), m.group(2).strip(), m.group(3) or "")
        for i, line in iter_lines_outside_code(lines)
        if (m := HEADER_RE.match(line))
    ]


def extract_links(lines: list[str]) -> list[Link]:
    """Extract markdown links outside code blocks."""
    links: list[Link] = []
    for i, line in iter_lines_outside_code(lines):
        for m in MARKDOWN_LINK_RE.finditer(line):
            links.append(
                Link(
                    i,
                    m.start(),
                    m.end(),
                    m.group("text"),
                    m.group("url"),
                    m.group("title"),
                    m.group("attrs"),
                )
            )
    return links


def extract_html_links(lines: list[str]) -> list[HtmlLink]:
    """Extract HTML ``<a>`` links outside code blocks."""
    links: list[HtmlLink] = []
    for i, line in iter_lines_outside_code(lines):
        for m in HTML_LINK_RE.finditer(line):
            links.append(HtmlLink(i, m.start(), m.end(), m.group(1), m.group(2)))
    return links


def extract_code_blocks(lines: list[str]) -> list[CodeBlock]:
    """Extract code blocks with their content lines."""
    blocks: list[CodeBlock] = []
    in_code, fence, lang, start = False, "", "", 0
    block_lines: list[str] = []

    for i, line in enumerate(lines):
        stripped = line.lstrip()
        if not in_code:
            if m := CODE_FENCE_RE.match(stripped):
                in_code, fence, lang, start = True, m.group(1), m.group(2), i
                block_lines = [line]
        else:
            block_lines.append(line)
            if stripped.startswith(fence) and stripped.strip() == fence:
                blocks.append(CodeBlock(start, i, lang, block_lines))
                in_code, fence, block_lines = False, "", []
    return blocks


# ---------------------------------------------------------------------------
# Fixers
# ---------------------------------------------------------------------------


def _extract_comment(line: str, lang: str) -> str | None:
    """Extract comment portion from a code line based on language."""
    if lang in HASH_COMMENT_LANGS or lang in MIXED_COMMENT_LANGS:
        if m := HASH_COMMENT_RE.match(line):
            return m.group(2)
    if lang in SLASH_COMMENT_LANGS or lang in MIXED_COMMENT_LANGS:
        if m := SLASH_COMMENT_RE.match(line):
            return m.group(2)
    return None


def fix_headers(trans: list[str], source_headers: list[Header]) -> list[str]:
    """Replace header anchors with those from source."""
    trans_headers = extract_headers(trans)
    if len(trans_headers) != len(source_headers):
        raise StructureMismatchError(
            f"Header count: {len(trans_headers)} vs {len(source_headers)}"
        )

    result = trans.copy()
    for th, sh in zip(trans_headers, source_headers):
        if th.level != sh.level:
            raise StructureMismatchError(f"Header level mismatch at line {th.line_no}")
        result[th.line_no] = f"{th.level} {th.title}{sh.anchor}"
    return result


def fix_links(trans: list[str], source_links: list[Link]) -> list[str]:
    """Replace link URLs/attrs with source values, keeping translated text."""
    trans_links = extract_links(trans)
    if len(trans_links) != len(source_links):
        raise StructureMismatchError(
            f"Link count: {len(trans_links)} vs {len(source_links)}"
        )

    result = trans.copy()
    for tl, sl in reversed(list(zip(trans_links, source_links))):
        new = f"[{tl.text}]({sl.url}"
        if tl.title:
            new += f' "{tl.title}"'
        new += ")"
        if sl.attrs:
            new += f"{{{sl.attrs}}}"
        result[tl.line_no] = (
            result[tl.line_no][: tl.start] + new + result[tl.line_no][tl.end :]
        )
    return result


def fix_html_links(trans: list[str], source_links: list[HtmlLink]) -> list[str]:
    """Replace HTML link attributes with source, keeping translated text."""
    trans_links = extract_html_links(trans)
    if len(trans_links) != len(source_links):
        raise StructureMismatchError(
            f"HTML link count: {len(trans_links)} vs {len(source_links)}"
        )

    result = trans.copy()
    for tl, sl in reversed(list(zip(trans_links, source_links))):
        new = f"<a {sl.attrs}>{tl.text}</a>"
        result[tl.line_no] = (
            result[tl.line_no][: tl.start] + new + result[tl.line_no][tl.end :]
        )
    return result


def fix_code_block(trans_block: CodeBlock, src_block: CodeBlock) -> list[str]:
    """Fix a single code block: use source code, preserve translated comments."""
    if trans_block.lang != src_block.lang:
        raise StructureMismatchError(
            f"Code lang mismatch at line {trans_block.start_line}"
        )
    if len(trans_block.lines) != len(src_block.lines):
        raise StructureMismatchError(
            f"Code lines mismatch at line {trans_block.start_line}"
        )

    lang = trans_block.lang.lower()

    # Mermaid diagrams are structural — use source entirely
    if lang == "mermaid":
        return src_block.lines.copy()

    result: list[str] = []
    in_block_comment = False

    for trans_line, src_line in zip(trans_block.lines, src_block.lines):
        # Fence lines: use source
        if trans_line.lstrip().startswith("```"):
            result.append(src_line)
            continue

        # Block comments in mixed languages: keep translation
        if lang in MIXED_COMMENT_LANGS:
            if BLOCK_COMMENT_START_RE.match(trans_line):
                in_block_comment = True
            if in_block_comment:
                result.append(trans_line)
                if BLOCK_COMMENT_END_RE.search(trans_line):
                    in_block_comment = False
                continue

        # Line comments: use source code with translated comment
        trans_comment = _extract_comment(trans_line, lang)
        if trans_comment:
            src_comment = _extract_comment(src_line, lang)
            if src_comment:
                result.append(src_line.replace(src_comment, trans_comment))
                continue

        # Default: use source line
        result.append(src_line)

    return result


def fix_code_blocks(trans: list[str], source_blocks: list[CodeBlock]) -> list[str]:
    """Fix all code blocks in translation."""
    trans_blocks = extract_code_blocks(trans)
    if len(trans_blocks) != len(source_blocks):
        raise StructureMismatchError(
            f"Code block count: {len(trans_blocks)} vs {len(source_blocks)}"
        )

    result = trans.copy()
    for tb, sb in zip(trans_blocks, source_blocks):
        fixed = fix_code_block(tb, sb)
        for i, line in enumerate(fixed):
            result[tb.start_line + i] = line
    return result


# ---------------------------------------------------------------------------
# Glossary loading and glossary-driven fixers
# ---------------------------------------------------------------------------

TAB_LABEL_RE = re.compile(r'^(\s*===\s*)"(.*)"(.*)$')
ADMONITION_RE = re.compile(r"^(!{3}|\?{3}\+?)\s+(\w[\w-]*)\s*(.*)")
NOTICE_SPAN_RE = re.compile(r'class="ai-translation-notice"')


@lru_cache
def load_glossary(lang: str) -> dict:
    """Load the per-language glossary YAML file."""
    glossary_path = DOCS_ROOT / lang / "glossary.yml"
    if not glossary_path.exists():
        return {}
    return yaml.safe_load(glossary_path.read_text(encoding="utf-8")) or {}


def fix_tab_labels(trans: list[str], source: list[str], lang: str) -> list[str]:
    """Replace tab labels with canonical glossary translations.

    Ensures consistent Before/After labels for MkDocs tab synchronization.
    Uses source lines to identify which English label each tab corresponds to,
    then applies the canonical translation. This handles variant translations
    (e.g., Korean using three different Before/After phrasings).
    Labels not in the glossary are left as-is (LLM translation kept).
    """
    glossary = load_glossary(lang)
    tab_labels = glossary.get("tab_labels", {})
    if not tab_labels:
        return trans

    # Collect source tab labels in order (outside code blocks)
    source_labels: list[str] = []
    for _, line in iter_lines_outside_code(source):
        m = TAB_LABEL_RE.match(line)
        if m:
            source_labels.append(m.group(2))

    # Collect translation tab label positions (outside code blocks)
    trans_positions: list[int] = []
    for i, line in iter_lines_outside_code(trans):
        if TAB_LABEL_RE.match(line):
            trans_positions.append(i)

    result = trans.copy()

    # Paired mode: match source and translation tab labels by position
    if len(source_labels) == len(trans_positions):
        for src_label, trans_idx in zip(source_labels, trans_positions):
            if src_label in tab_labels:
                m = TAB_LABEL_RE.match(result[trans_idx])
                prefix, suffix = m.group(1), m.group(3)
                result[trans_idx] = f'{prefix}"{tab_labels[src_label]}"{suffix}'
    else:
        # Fallback: replace any English labels found in translation
        for i in trans_positions:
            m = TAB_LABEL_RE.match(result[i])
            label = m.group(2)
            if label in tab_labels:
                prefix, suffix = m.group(1), m.group(3)
                result[i] = f'{prefix}"{tab_labels[label]}"{suffix}'

    return result


def fix_admonitions(trans: list[str], source: list[str], lang: str) -> list[str]:
    """Fix admonition keywords and titles using glossary.

    - Keywords are preserved from English source (prevents translated keywords)
    - Bare admonitions get glossary titles (e.g., !!! tip -> !!! tip "Tipp")
    - Standard titled admonitions get glossary translations
    - Custom one-off titles are left as-is (LLM translation kept)
    """
    glossary = load_glossary(lang)
    bare_titles = glossary.get("admonition_titles", {})
    title_glossary = glossary.get("admonition_title_glossary", {})
    if not bare_titles and not title_glossary:
        return trans

    # Extract admonition lines from source for pairing (outside code blocks)
    src_admns = []
    for i, line in iter_lines_outside_code(source):
        m = ADMONITION_RE.match(line.strip())
        if m:
            src_admns.append((i, m.group(1), m.group(2), m.group(3).strip()))

    trans_admns = []
    for i, line in iter_lines_outside_code(trans):
        m = ADMONITION_RE.match(line.strip())
        if m:
            trans_admns.append((i, m.group(1), m.group(2), m.group(3).strip()))

    result = trans.copy()

    if len(src_admns) == len(trans_admns):
        # Paired mode: fix using source as reference
        for (_, s_marker, s_kw, s_title), (t_idx, t_marker, _, _) in zip(
            src_admns, trans_admns
        ):
            indent = len(result[t_idx]) - len(result[t_idx].lstrip())
            pad = result[t_idx][:indent]
            # Use source keyword (preserves English + casing)
            kw = s_kw
            # Determine title
            if s_title:
                # Source has a title — strip quotes to get text
                en_title = s_title.strip('"')
                if en_title in title_glossary:
                    title = f' "{title_glossary[en_title]}"'
                else:
                    # Custom title: keep whatever the LLM produced
                    m2 = ADMONITION_RE.match(result[t_idx].strip())
                    title = (
                        f" {m2.group(3).strip()}" if m2 and m2.group(3).strip() else ""
                    )
            else:
                # Source is bare: apply glossary default title if available
                kw_lower = kw.lower()
                if kw_lower in bare_titles:
                    title = f' "{bare_titles[kw_lower]}"'
                else:
                    title = ""
            result[t_idx] = f"{pad}{s_marker} {kw}{title}"
    else:
        # Unpaired fallback: fix only recognizable standard titles
        # and bare admonitions using glossary, without source pairing
        for t_idx, t_marker, t_kw, t_title in trans_admns:
            indent = len(result[t_idx]) - len(result[t_idx].lstrip())
            pad = result[t_idx][:indent]
            if t_title:
                en_title = t_title.strip('"')
                # Check if this title matches a known English standard title
                if en_title in title_glossary:
                    result[t_idx] = (
                        f'{pad}{t_marker} {t_kw} "{title_glossary[en_title]}"'
                    )
                # Also check if this is already a translated title that
                # matches a glossary value — leave as-is in that case
            else:
                # Bare admonition: add glossary title
                kw_lower = t_kw.lower()
                if kw_lower in bare_titles:
                    result[t_idx] = f'{pad}{t_marker} {t_kw} "{bare_titles[kw_lower]}"'
    return result


def _split_frontmatter(text: str) -> tuple[str | None, str]:
    """Extract the raw frontmatter string and body from text.

    Returns (raw_fm_text_between_dashes, body) or (None, text).
    Handles code-fence-wrapped frontmatter (```yaml / ```markdown).
    """
    lines = text.split("\n")

    # Detect code-fence-wrapped frontmatter
    if lines and lines[0].strip().startswith("```"):
        fence_end = None
        for i in range(1, len(lines)):
            if lines[i].strip() == "```":
                fence_end = i
                break
        if fence_end is not None:
            inner = "\n".join(lines[1:fence_end])
            rest = "\n".join(lines[fence_end + 1 :])
            inner_raw, inner_body = _split_frontmatter(inner)
            if inner_raw is not None:
                return inner_raw, (inner_body + "\n" + rest).lstrip("\n")
        return None, text

    if not lines or lines[0].strip() != "---":
        return None, text

    for i in range(1, len(lines)):
        if lines[i].strip() == "---":
            return "\n".join(lines[1:i]), "\n".join(lines[i + 1 :])

    return None, text


def _parse_frontmatter(text: str) -> tuple[dict | None, str]:
    """Split text into (frontmatter_dict, body).

    Returns (None, text) if no valid frontmatter found.
    """
    raw, body = _split_frontmatter(text)
    if raw is None:
        return None, text
    try:
        fm = yaml.safe_load(raw)
    except yaml.YAMLError:
        return None, text
    if not isinstance(fm, dict):
        return None, text
    return fm, body


# Keys whose values should be taken from the translation (translatable).
_TRANSLATABLE_FM_KEYS = {"title", "description"}
_TRANSLATABLE_FM_NESTED = {
    # parent_key -> set of child keys that are translatable
    "additional_information": {"learning_objectives", "audience_prerequisites"},
}


def _replace_fm_values(en_raw: str, en_fm: dict, trans_fm: dict | None) -> str:
    """Replace translatable values in the English raw frontmatter string.

    Uses string-level replacement to preserve exact formatting.
    """
    if trans_fm is None:
        return en_raw

    result = en_raw

    # Replace top-level translatable scalar keys
    for key in _TRANSLATABLE_FM_KEYS:
        if key in en_fm and key in trans_fm and en_fm[key] != trans_fm[key]:
            en_val = en_fm[key]
            trans_val = trans_fm[key]
            # Replace the value after "key: " — handle both quoted and unquoted
            for pattern in [
                re.compile(
                    rf"^({key}:\s*)" + re.escape(str(en_val)) + r"$", re.MULTILINE
                ),
                re.compile(
                    rf'^({key}:\s*)"' + re.escape(str(en_val)) + r'"$',
                    re.MULTILINE,
                ),
                re.compile(
                    rf"^({key}:\s*)'" + re.escape(str(en_val)) + r"'$",
                    re.MULTILINE,
                ),
            ]:
                m = pattern.search(result)
                if m:
                    # Preserve the original quoting style
                    prefix = m.group(1)
                    original = m.group(0)
                    if original.endswith('"'):
                        replacement = f'{prefix}"{trans_val}"'
                    elif original.endswith("'"):
                        replacement = f"{prefix}'{trans_val}'"
                    else:
                        replacement = f"{prefix}{trans_val}"
                    result = result[: m.start()] + replacement + result[m.end() :]
                    break

    # Replace translatable nested list values
    for parent_key, child_keys in _TRANSLATABLE_FM_NESTED.items():
        en_parent = en_fm.get(parent_key)
        trans_parent = trans_fm.get(parent_key) if trans_fm else None
        if not isinstance(en_parent, dict) or not isinstance(trans_parent, dict):
            continue
        for child_key in child_keys:
            en_list = en_parent.get(child_key)
            trans_list = trans_parent.get(child_key)
            if (
                not isinstance(en_list, list)
                or not isinstance(trans_list, list)
                or len(en_list) != len(trans_list)
            ):
                continue
            for en_item, trans_item in zip(en_list, trans_list):
                if en_item != trans_item and isinstance(en_item, str):
                    # Replace anchored to YAML list item context
                    # Handle optional YAML quoting (", ', or unquoted)
                    escaped = re.escape(en_item)
                    pattern = re.compile(
                        rf'^(\s*-\s*)["\']?{escaped}["\']?$',
                        re.MULTILINE,
                    )
                    m = pattern.search(result)
                    if m:
                        # Preserve original quoting style
                        original_line = m.group(0)
                        prefix = m.group(1)
                        after_prefix = original_line[len(prefix) :]
                        if after_prefix.startswith('"'):
                            replacement = f'{prefix}"{trans_item}"'
                        elif after_prefix.startswith("'"):
                            replacement = f"{prefix}'{trans_item}'"
                        else:
                            replacement = f"{prefix}{trans_item}"
                        result = result[: m.start()] + replacement + result[m.end() :]

    return result


def fix_frontmatter(text: str, source_text: str, lang: str) -> str:
    """Ensure frontmatter matches English source structure.

    Uses string-level replacement to preserve exact formatting from
    the English source, substituting only translatable values.
    If English has no frontmatter, any spurious frontmatter is stripped.
    If translation is missing frontmatter, English structure is restored.
    """
    en_fm, _ = _parse_frontmatter(source_text)
    if en_fm is None:
        # English has no frontmatter — strip any the LLM may have added
        _, body = _parse_frontmatter(text)
        return body.lstrip("\n")

    en_raw, _ = _split_frontmatter(source_text)
    trans_fm, body = _parse_frontmatter(text)
    fm_text = _replace_fm_values(en_raw, en_fm, trans_fm)

    # Ensure blank line between closing --- and body content
    if body and not body.startswith("\n"):
        body = "\n" + body
    return f"---\n{fm_text}\n---\n{body}"


def _build_notice_span(lang: str) -> str | None:
    """Build the canonical translation notice span for a language."""
    glossary = load_glossary(lang)
    notice = glossary.get("translation_notice", {})
    text = notice.get("text")
    if not text:
        return None
    return (
        f'<span class="ai-translation-notice">'
        f":material-information-outline:{{ .ai-translation-notice-icon }} "
        f"{text}</span>"
    )


def _build_notice_admonition(lang: str) -> str | None:
    """Build the canonical homepage notice admonition for a language."""
    glossary = load_glossary(lang)
    notice = glossary.get("translation_notice", {})
    homepage = notice.get("homepage", {})
    title = homepage.get("title")
    body = homepage.get("body")
    if not title or not body:
        return None
    lines = [f'!!! note "{title}"', ""]
    for bline in body.strip().split("\n"):
        lines.append(f"    {bline}")
    return "\n".join(lines)


def ensure_translation_notice(text: str, lang: str, is_homepage: bool = False) -> str:
    """Ensure the AI translation notice is present exactly once.

    For regular pages: inserts a <span> after the first H1 heading.
    For the homepage: inserts a !!! note admonition before the first H2.
    Idempotent: replaces malformed notices, never duplicates.
    """
    lines = text.split("\n")

    # Check for existing notice and its location
    existing_idx = None
    for i, line in enumerate(lines):
        if NOTICE_SPAN_RE.search(line):
            existing_idx = i
            break

    if is_homepage:
        notice = _build_notice_admonition(lang)
        if not notice:
            return text

        # Check for existing admonition notice on homepage
        existing_admn_idx = None
        for i, line in enumerate(lines):
            if NOTICE_SPAN_RE.search(line):
                # Span on homepage — remove it, we'll use admonition
                existing_admn_idx = i
                break
            if "ai-translation" in line.lower() or (
                line.strip().startswith("!!! note")
                and i + 2 < len(lines)
                and "TRANSLATING.md" in "".join(lines[i : i + 5])
            ):
                existing_admn_idx = i
                break

        if existing_admn_idx is not None:
            # Remove existing notice block (admonition is ~5 lines)
            end = existing_admn_idx + 1
            while end < len(lines) and (
                lines[end].startswith("    ") or lines[end].strip() == ""
            ):
                end += 1
            # Remove blank lines before the notice too
            start = existing_admn_idx
            while start > 0 and lines[start - 1].strip() == "":
                start -= 1
            lines = lines[:start] + lines[end:]

        # Find insertion point: before first ## heading
        insert_idx = None
        for i, line in enumerate(lines):
            if line.startswith("## "):
                insert_idx = i
                break

        if insert_idx is not None:
            lines = lines[:insert_idx] + [notice, ""] + lines[insert_idx:]
            # Ensure blank line before notice
            if insert_idx > 0 and lines[insert_idx - 1].strip() != "":
                lines.insert(insert_idx, "")

        return "\n".join(lines)

    # Regular page: span format
    notice_span = _build_notice_span(lang)
    if not notice_span:
        return text

    if existing_idx is not None:
        # Replace existing (possibly malformed) notice with canonical version
        lines[existing_idx] = notice_span
        return "\n".join(lines)

    # No existing notice: find first H1 heading and insert after it
    h1_idx = None
    for i, line in enumerate(lines):
        if line.startswith("# ") and not line.startswith("## "):
            h1_idx = i
            break

    if h1_idx is None:
        return text  # No H1 found, skip

    # Insert: blank line, notice, blank line after H1
    insert_at = h1_idx + 1
    # Skip any existing blank line after heading
    if insert_at < len(lines) and lines[insert_at].strip() == "":
        insert_at += 1
    lines = lines[: h1_idx + 1] + ["", notice_span, ""] + lines[insert_at:]
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Combined translation fixer and file-level post-processing
# ---------------------------------------------------------------------------


def fix_translation(translation: list[str], source: list[str], lang: str) -> list[str]:
    """Apply all structural fixes to a translation against English source.

    Paired fixes (headers, links, code blocks) require matching structure
    between source and translation and may raise StructureMismatchError.
    Glossary-based fixes (admonitions, tab labels) are applied separately
    and never fail.
    """
    result = translation.copy()
    # Paired fixes — may raise StructureMismatchError
    result = fix_headers(result, extract_headers(source))
    result = fix_links(result, extract_links(source))
    result = fix_html_links(result, extract_html_links(source))
    result = fix_code_blocks(result, extract_code_blocks(source))
    # Glossary-based fixes — always safe
    result = fix_admonitions(result, source, lang)
    result = fix_tab_labels(result, source, lang)
    return result


# ---------------------------------------------------------------------------
# File-level post-processing (convenience wrapper with I/O)
# ---------------------------------------------------------------------------


def post_process_file(lang_path: Path, lang: str) -> bool:
    """Post-process a translated file in-place. Returns True if changes were made.

    Reads the translation and its English source, runs structural fixes
    (headers, links, code blocks, admonitions, tab labels), then applies
    frontmatter merging and translation notice injection.
    """
    from .paths import lang_to_en_path

    en_path = lang_to_en_path(lang_path, lang)
    if not en_path.exists():
        return False

    original = lang_path.read_text(encoding="utf-8")
    source_text = en_path.read_text(encoding="utf-8")

    source_lines = source_text.splitlines()
    try:
        fixed = fix_translation(original.splitlines(), source_lines, lang)
    except StructureMismatchError as e:
        console = make_console()
        console.print(
            f"[yellow]Post-process paired fixes skipped for "
            f"{lang_path.relative_to(REPO_ROOT)}: {e}[/yellow]"
        )
        # Paired fixes failed, but still apply glossary-based fixes
        fixed = original.splitlines()
        fixed = fix_admonitions(fixed, source_lines, lang)
        fixed = fix_tab_labels(fixed, source_lines, lang)

    fixed_text = "\n".join(fixed) + "\n"

    # Frontmatter: merge structural keys from English, keep translated values
    fixed_text = fix_frontmatter(fixed_text, source_text, lang)

    # Translation notice: ensure present exactly once
    is_homepage = lang_path.name == "index.md" and lang_path.parent.name == "docs"
    fixed_text = ensure_translation_notice(fixed_text, lang, is_homepage=is_homepage)

    # Ensure file ends with single newline
    fixed_text = fixed_text.rstrip("\n") + "\n"

    if fixed_text != original:
        lang_path.write_text(fixed_text, encoding="utf-8", newline="\n")
        return True
    return False
