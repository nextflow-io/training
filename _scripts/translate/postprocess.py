"""Post-processing: fix structural elements in translated MDX.

After Claude translates a file, this module makes headings, links, code
blocks, JSX tags, asides, frontmatter, imports, image paths and the AI
translation notice match the English source structurally. Code blocks are
restored from source with only translated comments, ``step``/``label``
values and ``[!annotate]`` text preserved. Glossary-controlled elements
(aside and ``<summary>`` titles, Before/After labels, notices) are enforced
deterministically from per-language YAML glossary files.
"""

from __future__ import annotations

import copy
import html
import os
import re
import unicodedata
from functools import lru_cache
from pathlib import Path

import yaml

from .config import (
    DOCS_ROOT,
    HASH_COMMENT_LANGS,
    LANGUAGES_DIR,
    MIXED_COMMENT_LANGS,
    REPO_ROOT,
    SLASH_COMMENT_LANGS,
    StructureMismatchError,
    make_console,
)
from .models import CodeBlock, Header, Link

# ---------------------------------------------------------------------------
# Compiled regex patterns
# ---------------------------------------------------------------------------

HEADER_RE = re.compile(r"^(#{1,6})\s+(.+?)(\s*\{#[^}]+\})?\s*$")
MARKDOWN_LINK_RE = re.compile(
    r"(?<!!)\[(?P<text>[^\]]*)\]\((?P<url>[^)\s]+)(?:\s+\"(?P<title>[^\"]*)\")?\)(?:\{(?P<attrs>[^}]*)\})?"
)
IMAGE_RE = re.compile(r"!\[[^\]]*\]\((?P<src>[^)\s]+)")
CODE_FENCE_RE = re.compile(r"^(`{3,})([\w-]*)")
HASH_COMMENT_RE = re.compile(r"^(.*?)((?:\s+#|^#)\s.*)$")
SLASH_COMMENT_RE = re.compile(r"^(.*?)((?:\s+//|^//)\s.*)$")
BLOCK_COMMENT_START_RE = re.compile(r"^\s*/\*")
BLOCK_COMMENT_END_RE = re.compile(r"\*/\s*$")

INLINE_CODE_RE = re.compile(r"(`+)(.+?)(?<!`)\1(?!`)")
# An HTML/JSX tag on one line; attribute values may hold `>` inside quotes or {expressions}.
TAG_RE = re.compile(
    r"<(?P<close>/?)(?P<name>[A-Za-z][\w.-]*)(?P<attrs>(?:[^<>\"'{}]|\"[^\"]*\"|'[^']*'|\{[^{}]*\})*?)(?P<self>/?)>"
)
ATTR_RE = re.compile(r'([\w-]+)="([^"]*)"')
DIRECTIVE_RE = re.compile(r"^(\s*)(:{3,})([\w-]*)(?:\[(.*)\])?\s*$")
SUMMARY_RE = re.compile(r"<summary>(.*?)</summary>")
IMPORT_RE = re.compile(r"^import (?P<name>\w+) from '(?P<path>[^']+)';\s*$")
FENCE_LABEL_RE = re.compile(r'\b(step|label)="([^"]*)"')
COMMENT_DIRECTIVE_RE = re.compile(r"\[!\w[^\]]*\]")
OUTPUT_MARKER_RE = re.compile(r"^\s*# \[!output\b")
ANCHOR_LINK_RE = re.compile(
    r'(?:\]\(|href=")(?P<path>[^)\s"#:]*)#(?P<hash>[^)\s"]+)'
)
URL_RE = re.compile(r"^([a-z][a-z\d+.-]*:|/|#)", re.I)

# JSX attributes whose values are prose; every other attribute is copied from English.
TRANSLATABLE_ATTRS = {
    "TabItem": {"label"},
    "Card": {"title"},
    "LinkCard": {"title", "description"},
}


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


def mask_inline_code(line: str) -> str:
    """Blank out inline code spans, keeping positions, so tag patterns skip them."""
    return INLINE_CODE_RE.sub(lambda m: " " * len(m.group(0)), line)


def iter_tags(lines: list[str]):
    """Yield ``(line_no, match)`` for HTML/JSX tags outside code."""
    for i, line in iter_lines_outside_code(lines):
        for m in TAG_RE.finditer(mask_inline_code(line)):
            yield i, m


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


def extract_directives(lines: list[str]) -> list[tuple[int, re.Match]]:
    """``:::name[title]`` openers and ``:::`` closers outside code."""
    return [
        (i, m)
        for i, line in iter_lines_outside_code(lines)
        if (m := DIRECTIVE_RE.match(line))
    ]


def extract_imports(lines: list[str]) -> list[tuple[int, re.Match]]:
    """ESM ``import x from '...';`` lines outside code."""
    return [
        (i, m) for i, line in iter_lines_outside_code(lines) if (m := IMPORT_RE.match(line))
    ]


# ---------------------------------------------------------------------------
# Paired fixers (translation vs English source)
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


def _fix_fence_line(src_line: str, trans_line: str, tab_labels: dict) -> str:
    """Source fence line, with the translation's ``step``/``label`` values."""
    trans_values = dict(FENCE_LABEL_RE.findall(trans_line))

    def value(m: re.Match) -> str:
        key, en = m.group(1), m.group(2)
        return f'{key}="{tab_labels.get(en) or trans_values.get(key, en)}"'

    return FENCE_LABEL_RE.sub(value, src_line)


def fix_code_block(
    trans_block: CodeBlock, src_block: CodeBlock, tab_labels: dict | None = None
) -> list[str]:
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

    result = [_fix_fence_line(src_block.lines[0], trans_block.lines[0], tab_labels or {})]
    in_block_comment = False
    in_output = False

    for trans_line, src_line in zip(trans_block.lines[1:], src_block.lines[1:]):
        # Closing fence and terminal output (from `# [!output]` to the end): use source
        in_output = in_output or bool(OUTPUT_MARKER_RE.match(src_line))
        if in_output or trans_line.lstrip().startswith("```"):
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

        # Line comments: use source code with translated comment, as long as
        # the comment keeps its `[!annotate]`-style markers
        trans_comment = _extract_comment(trans_line, lang)
        if trans_comment:
            src_comment = _extract_comment(src_line, lang)
            if src_comment and COMMENT_DIRECTIVE_RE.findall(
                src_comment
            ) == COMMENT_DIRECTIVE_RE.findall(trans_comment):
                result.append(src_line.replace(src_comment, trans_comment))
                continue

        # Default: use source line
        result.append(src_line)

    return result


def fix_code_blocks(
    trans: list[str], source_blocks: list[CodeBlock], tab_labels: dict | None = None
) -> list[str]:
    """Fix all code blocks in translation."""
    trans_blocks = extract_code_blocks(trans)
    if len(trans_blocks) != len(source_blocks):
        raise StructureMismatchError(
            f"Code block count: {len(trans_blocks)} vs {len(source_blocks)}"
        )

    result = trans.copy()
    for tb, sb in zip(trans_blocks, source_blocks):
        fixed = fix_code_block(tb, sb, tab_labels)
        for i, line in enumerate(fixed):
            result[tb.start_line + i] = line
    return result


def fix_tags(trans: list[str], source: list[str], tab_labels: dict | None = None) -> list[str]:
    """Copy every HTML/JSX tag from source, keeping translated prose attributes."""
    src_tags = list(iter_tags(source))
    trans_tags = list(iter_tags(trans))
    shape = lambda tags: [(m["close"], m["name"], m["self"]) for _, m in tags]  # noqa: E731
    if shape(src_tags) != shape(trans_tags):
        raise StructureMismatchError(
            f"Tags differ: {len(trans_tags)} vs {len(src_tags)} in English"
        )

    labels = tab_labels or {}
    result = trans.copy()
    for (si, sm), (ti, tm) in reversed(list(zip(src_tags, trans_tags))):
        new = source[si][sm.start() : sm.end()]
        keep = TRANSLATABLE_ATTRS.get(sm["name"], set())
        if keep:
            trans_values = dict(ATTR_RE.findall(tm["attrs"]))

            def value(a: re.Match) -> str:
                if a[1] not in keep:
                    return a[0]
                if a[1] == "label" and a[2] in labels:
                    return f'{a[1]}="{labels[a[2]]}"'
                return f'{a[1]}="{trans_values.get(a[1], a[2])}"'

            new = ATTR_RE.sub(value, new)
        result[ti] = result[ti][: tm.start()] + new + result[ti][tm.end() :]
    return result


def fix_directives(trans: list[str], source: list[str], titles: dict | None = None) -> list[str]:
    """Restore ``:::`` aside names and nesting from source; glossary titles win."""
    src = extract_directives(source)
    tr = extract_directives(trans)
    if len(src) != len(tr):
        raise StructureMismatchError(f"Aside count: {len(tr)} vs {len(src)}")

    titles = titles or {}
    result = trans.copy()
    for (_, sm), (ti, tm) in zip(src, tr):
        indent, colons, name, en_title = sm.groups()
        title = titles.get(en_title) or tm.group(4) or en_title
        result[ti] = f"{indent}{colons}{name}" + (f"[{title}]" if name and title else "")
    return result


def fix_summaries(trans: list[str], source: list[str], titles: dict) -> list[str]:
    """Apply glossary titles to ``<summary>`` text whose English is in the glossary."""
    find = lambda lines: [  # noqa: E731
        (i, m) for i, line in iter_lines_outside_code(lines) for m in SUMMARY_RE.finditer(line)
    ]
    src, tr = find(source), find(trans)
    if len(src) != len(tr):
        raise StructureMismatchError(f"<summary> count: {len(tr)} vs {len(src)}")

    result = trans.copy()
    for (_, sm), (ti, tm) in reversed(list(zip(src, tr))):
        if sm.group(1) in titles:
            line = result[ti]
            result[ti] = line[: tm.start(1)] + titles[sm.group(1)] + line[tm.end(1) :]
    return result


# ---------------------------------------------------------------------------
# Asset paths and anchors (translations live one directory deeper)
# ---------------------------------------------------------------------------


def rebase(path: str, rel: str, lang: str) -> str:
    """Re-point a path written relative to English page ``rel`` from its ``<lang>/`` copy."""
    if URL_RE.match(path):
        return path
    file, sep, query = path.partition("?")
    target = os.path.normpath(os.path.join(os.path.dirname(rel), file))
    new = os.path.relpath(target, os.path.join(lang, os.path.dirname(rel)))
    new = new.replace(os.sep, "/")
    return (new if new.startswith(".") else f"./{new}") + sep + query


def fix_images(trans: list[str], rel: str, lang: str, root: Path = DOCS_ROOT) -> list[str]:
    """Point relative images at the English tree, where images live."""
    lang_dir = (root / lang / rel).parent
    en_dir = (root / rel).parent
    result = trans.copy()
    for i, line in iter_lines_outside_code(trans):

        def sub(m: re.Match) -> str:
            src = m["src"]
            if URL_RE.match(src) or (lang_dir / src).exists() or not (en_dir / src).exists():
                return m[0]
            return m[0][: m.start("src") - m.start()] + rebase(src, rel, lang)

        result[i] = IMAGE_RE.sub(sub, line)
    return result


def _github_slug(text: str) -> str:
    """github-slugger: keeps letters, marks, digits, `_` and `-`; spaces become `-`."""

    def keep(c: str) -> bool:
        cat = unicodedata.category(c)
        return c in " -" or cat[0] in "LM" or cat in ("Nd", "Nl", "Pc")

    return "".join(c for c in text.lower() if keep(c)).replace(" ", "-")


def _heading_text(md: str) -> str:
    """Text of inline Markdown as the heading slugger sees it."""
    spans: list[str] = []

    def stash(m: re.Match) -> str:
        code = re.sub(r"^([^`]*\S)\{:[\w#+-][\w#+.-]*\}$", r"\1", m[2])
        spans.append(code)
        return f"\0{len(spans) - 1}\0"

    s = INLINE_CODE_RE.sub(stash, md)
    s = re.sub(r"!\[[^\]]*\]\([^)]*\)", "", s)
    s = re.sub(r"\[([^\]]*)\]\([^)]*\)", r"\1", s)
    s = re.sub(r"\{/\*.*?\*/\}", "", s)
    s = re.sub(r"</?[A-Za-z][^>]*>", "", s)
    s = re.sub(r"\*\*|__|~~|\*|(?<!\w)_|_(?!\w)", "", s)
    s = re.sub(r"\\([!-/:-@\[-`{-~])", r"\1", s)
    s = html.unescape(s)
    return re.sub(r"\0(\d+)\0", lambda m: spans[int(m[1])], s)


def heading_ids(lines: list[str]) -> list[tuple[int, str]]:
    """``(level, id)`` of each heading as Starlight renders it."""
    seen: dict[str, int] = {}
    out = []
    for h in extract_headers(lines):
        if h.anchor:
            out.append((len(h.level), h.anchor.strip()[2:-1]))
            continue
        base = slug = _github_slug(_heading_text(h.title))
        n = seen.get(base, 0)
        while slug in seen:
            n += 1
            slug = f"{base}-{n}"
        if slug != base:
            seen[base] = n
        seen[slug] = 0
        out.append((len(h.level), slug))
    return out


def fix_anchors(
    trans: list[str], source: list[str], rel: str, lang: str, root: Path = DOCS_ROOT
) -> list[str]:
    """Point each ``#anchor`` at the id its heading has in this locale's target page.

    Translated headings get ids from translated text, so an anchor is mapped
    through the English page: same heading position, when both pages have the
    same heading levels.
    """
    cache: dict[str, tuple[list, list]] = {}

    def targets(path: str) -> tuple[list, list] | None:
        if not path:
            return heading_ids(trans), heading_ids(source)
        if not path.endswith(".mdx"):
            return None
        target = os.path.normpath(os.path.join(os.path.dirname(rel), path))
        if target not in cache:
            en, local = root / target, root / lang / target
            if not en.is_file():
                return None
            en_ids = heading_ids(en.read_text(encoding="utf-8").splitlines())
            # A missing translation falls back to the English page and its ids.
            local_ids = (
                heading_ids(local.read_text(encoding="utf-8").splitlines())
                if local.is_file()
                else en_ids
            )
            cache[target] = (local_ids, en_ids)
        return cache[target]

    def find(lines):
        return [
            (i, m)
            for i, line in iter_lines_outside_code(lines)
            for m in ANCHOR_LINK_RE.finditer(mask_inline_code(line))
        ]

    src_links, trans_links = find(source), find(trans)
    paired = [m["path"] for _, m in src_links] == [m["path"] for _, m in trans_links]
    result = trans.copy()
    for k, (i, m) in reversed(list(enumerate(trans_links))):
        found = targets(m["path"])
        if not found:
            continue
        local, en = found
        if m["hash"] in [id_ for _, id_ in local]:
            continue
        key = src_links[k][1]["hash"] if paired else m["hash"]
        en_ids = [id_ for _, id_ in en]
        if key not in en_ids or [lv for lv, _ in local] != [lv for lv, _ in en]:
            continue
        new = local[en_ids.index(key)][1]
        result[i] = result[i][: m.start("hash")] + new + result[i][m.end("hash") :]
    return result


def fix_anchors_in_language(lang: str, root: Path = DOCS_ROOT) -> int:
    """Re-point anchors in every page of a locale. Returns the number of files changed.

    Re-translating a page can change its heading ids, which breaks anchors
    pointing at it from untouched pages (and the build's link check).
    """
    changed = 0
    for path in sorted((root / lang).rglob("*.mdx")):
        rel = path.relative_to(root / lang).as_posix()
        en = root / rel
        if not en.is_file():
            continue
        text = path.read_text(encoding="utf-8")
        lines = text.split("\n")
        fixed = fix_anchors(lines, en.read_text(encoding="utf-8").split("\n"), rel, lang, root)
        if fixed != lines:
            path.write_text("\n".join(fixed), encoding="utf-8", newline="\n")
            changed += 1
    return changed


# ---------------------------------------------------------------------------
# Glossary loading and the translation notice
# ---------------------------------------------------------------------------

NOTICE_SPAN_RE = re.compile(r'class="ai-translation-notice"')

# Matches the full notice span line (used for stripping LLM-generated notices)
_NOTICE_SPAN_LINE_RE = re.compile(
    r'^\s*<span\s+class="ai-translation-notice".*?</span>\s*$'
)


def strip_translation_notice(lines: list[str]) -> list[str]:
    """Remove any AI translation notice the LLM may have inserted.

    Handles both the span (regular pages) and the ``:::note`` aside
    (homepage) so that ``ensure_translation_notice`` can add the canonical
    version later. Returns a new list with the notice lines removed.
    """
    result: list[str] = []
    i = 0
    while i < len(lines):
        if NOTICE_SPAN_RE.search(lines[i]):
            # Also skip one surrounding blank line to avoid double-blanks:
            # the notice is normally padded by a blank line on each side.
            i += 1
            if i < len(lines) and lines[i].strip() == "":
                i += 1
            continue
        if lines[i].startswith(":::note"):
            end = next(
                (j for j in range(i + 1, len(lines)) if lines[j].strip() == ":::"), None
            )
            if end is not None and "TRANSLATING.md" in "\n".join(lines[i:end]):
                i = end + 1
                if i < len(lines) and lines[i].strip() == "":
                    i += 1
                continue
        result.append(lines[i])
        i += 1
    return result


@lru_cache
def load_glossary(lang: str) -> dict:
    """Load the per-language glossary YAML file."""
    glossary_path = LANGUAGES_DIR / lang / "glossary.yml"
    if not glossary_path.exists():
        return {}
    return yaml.safe_load(glossary_path.read_text(encoding="utf-8")) or {}


def title_glossary(lang: str) -> dict:
    """English aside / ``<summary>`` title -> canonical translation."""
    glossary = load_glossary(lang)
    titles = {
        kw.capitalize(): title
        for kw, title in (glossary.get("admonition_titles") or {}).items()
    }
    titles.update(glossary.get("admonition_title_glossary") or {})
    return titles


def _build_notice_span(lang: str) -> str | None:
    """Build the canonical translation notice span for a language."""
    text = load_glossary(lang).get("translation_notice", {}).get("text")
    return f'<span class="ai-translation-notice">{text}</span>' if text else None


def _build_notice_admonition(lang: str) -> str | None:
    """Build the canonical homepage notice aside for a language."""
    homepage = load_glossary(lang).get("translation_notice", {}).get("homepage", {})
    title = homepage.get("title")
    body = homepage.get("body")
    if not title or not body:
        return None
    return f":::note[{title}]\n{body.strip()}\n:::"


def ensure_translation_notice(text: str, lang: str, is_homepage: bool = False) -> str:
    """Ensure the AI translation notice is present exactly once.

    Regular pages get a ``<span>`` after the frontmatter and imports; the
    homepage gets a ``:::note`` aside before its first ``##`` heading.
    Idempotent: replaces existing notices, never duplicates.
    """
    notice = _build_notice_admonition(lang) if is_homepage else _build_notice_span(lang)
    if not notice:
        return text

    lines = strip_translation_notice(text.split("\n"))
    if is_homepage:
        idx = next(
            (i for i, line in iter_lines_outside_code(lines) if line.startswith("## ")),
            len(lines),
        )
    else:
        idx = 0
        if lines and lines[0].strip() == "---":
            idx = next(
                (i + 1 for i in range(1, len(lines)) if lines[i].strip() == "---"), 0
            )
        while idx < len(lines) and (lines[idx].strip() == "" or IMPORT_RE.match(lines[idx])):
            idx += 1

    before, after = lines[:idx], lines[idx:]
    while before and before[-1].strip() == "":
        before.pop()
    while after and after[0].strip() == "":
        after.pop(0)
    return "\n".join(before + ([""] if before else []) + [notice, ""] + after)


# ---------------------------------------------------------------------------
# Frontmatter and imports
# ---------------------------------------------------------------------------


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


def parse_frontmatter(text: str) -> tuple[dict | None, str]:
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


# Dotted key paths whose values are taken from the translation.
TRANSLATABLE_FM_KEYS = {
    "title",
    "description",
    "sidebar.label",
    "additional_information.learning_objectives",
    "additional_information.audience_prerequisites",
}


def _merge_frontmatter(en: dict, trans: dict | None, prefix: str = "") -> dict:
    """English frontmatter with the translated values of translatable keys."""
    trans = trans if isinstance(trans, dict) else {}
    out = {}
    for key, value in en.items():
        path, t = f"{prefix}{key}", trans.get(key)
        if isinstance(value, dict):
            out[key] = _merge_frontmatter(value, t, f"{path}.")
        elif (
            path in TRANSLATABLE_FM_KEYS
            and type(t) is type(value)
            and (not isinstance(value, list) or len(t) == len(value))
        ):
            out[key] = copy.deepcopy(t)
        else:
            out[key] = copy.deepcopy(value)
    return out


def dump_frontmatter(fm: dict) -> str:
    """YAML in the style the Starlight converter wrote, so untouched values diff cleanly."""
    return yaml.safe_dump(
        fm, sort_keys=False, allow_unicode=True, width=10_000, default_flow_style=False
    )


def fix_frontmatter(text: str, source_text: str, lang: str) -> str:
    """Ensure frontmatter matches English source structure.

    Keys, order and untranslatable values come from English; translatable
    values (``TRANSLATABLE_FM_KEYS``) from the translation. If English has
    no frontmatter, any spurious frontmatter is stripped.
    """
    en_fm, _ = parse_frontmatter(source_text)
    if en_fm is None:
        # English has no frontmatter — strip any the LLM may have added
        _, body = parse_frontmatter(text)
        body = body.lstrip("\n")
        # parse_frontmatter only recognizes a *closed* "---...---" block.
        # An unclosed leading "---" (e.g. the LLM emitted a stray/empty
        # frontmatter opener) isn't real frontmatter but would otherwise
        # survive as spurious content, so strip it here too.
        first_line, _, rest = body.partition("\n")
        if first_line.strip() == "---":
            body = rest.lstrip("\n")
        return body

    trans_fm, body = parse_frontmatter(text)
    if trans_fm is None:
        # Unparseable or missing: keep everything after a closed block as body
        _, body = _split_frontmatter(text)
    fm = _merge_frontmatter(en_fm, trans_fm)
    body = body.lstrip("\n")
    return f"---\n{dump_frontmatter(fm)}---\n\n{body}"


def fix_imports(text: str, source_text: str, rel: str, lang: str) -> str:
    """Use the English ``import`` lines, re-pointed from the translation's directory."""
    en_imports = [
        f"import {m['name']} from '{rebase(m['path'], rel, lang)}';"
        for _, m in extract_imports(source_text.split("\n"))
    ]
    raw, body = _split_frontmatter(text)
    body_lines = body.split("\n")
    drop = {i for i, _ in extract_imports(body_lines)}
    if not en_imports and not drop:
        return text
    body = "\n".join(line for i, line in enumerate(body_lines) if i not in drop).lstrip("\n")
    head = f"---\n{raw}\n---\n" if raw is not None else ""
    if en_imports:
        head += "\n" + "\n".join(en_imports) + "\n"
    return f"{head}\n{body}"


# ---------------------------------------------------------------------------
# Combined translation fixer and file-level post-processing
# ---------------------------------------------------------------------------


def fix_translation(
    translation: list[str], source: list[str], lang: str
) -> tuple[list[str], list[str]]:
    """Apply all structural fixes to a translation against English source.

    Each paired fix runs on its own: one that finds a structure mismatch is
    skipped and reported, and the others still apply.

    Returns:
        (fixed lines, mismatch messages)
    """
    tab_labels = load_glossary(lang).get("tab_labels") or {}
    titles = title_glossary(lang)
    fixers = [
        lambda t: fix_headers(t, extract_headers(source)),
        lambda t: fix_links(t, extract_links(source)),
        lambda t: fix_code_blocks(t, extract_code_blocks(source), tab_labels),
        lambda t: fix_tags(t, source, tab_labels),
        lambda t: fix_directives(t, source, titles),
        lambda t: fix_summaries(t, source, titles),
    ]
    result, problems = translation.copy(), []
    for fix in fixers:
        try:
            result = fix(result)
        except StructureMismatchError as e:
            problems.append(str(e))
    return result, problems


def post_process_text(
    text: str, source_text: str, rel: str, lang: str, root: Path = DOCS_ROOT
) -> tuple[str, list[str]]:
    """Post-process translated text of English page ``rel``.

    Runs the structural fixes, then re-points images and anchors, merges
    frontmatter and imports and adds the translation notice.

    Returns:
        (fixed text, mismatch messages from skipped fixes)
    """
    source_lines = source_text.splitlines()

    # Strip any notice the LLM inserted before the paired fixes: its link
    # does not exist in the English source and would skew the link count.
    lines = strip_translation_notice(text.splitlines())

    lines, problems = fix_translation(lines, source_lines, lang)
    lines = fix_images(lines, rel, lang, root)
    lines = fix_anchors(lines, source_lines, rel, lang, root)

    fixed_text = "\n".join(lines) + "\n"
    fixed_text = fix_frontmatter(fixed_text, source_text, lang)
    fixed_text = fix_imports(fixed_text, source_text, rel, lang)
    fixed_text = ensure_translation_notice(fixed_text, lang, is_homepage=rel == "index.mdx")
    return fixed_text.rstrip("\n") + "\n", problems


def post_process_file(lang_path: Path, lang: str) -> bool:
    """Post-process a translated file in-place. Returns True if changes were made."""
    from .paths import lang_to_en_path

    en_path = lang_to_en_path(lang_path, lang)
    if not en_path.exists():
        return False

    original = lang_path.read_text(encoding="utf-8")
    rel = en_path.relative_to(DOCS_ROOT).as_posix()
    fixed_text, problems = post_process_text(
        original, en_path.read_text(encoding="utf-8"), rel, lang
    )
    if problems:
        make_console().print(
            f"[yellow]Post-process fixes skipped for "
            f"{lang_path.relative_to(REPO_ROOT)}: {'; '.join(problems)}[/yellow]"
        )

    if fixed_text != original:
        lang_path.write_text(fixed_text, encoding="utf-8", newline="\n")
        return True
    return False
