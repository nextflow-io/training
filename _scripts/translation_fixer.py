#!/usr/bin/env python3
# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "typer>=0.12.0",
#     "rich>=13.0.0",
# ]
# ///
"""
Post-process translations to fix common LLM mistakes.

Fixes:
- Header anchors/permalinks (preserve from English)
- Markdown links (preserve URLs, keep translated text)
- HTML links (preserve attributes, keep translated text)
- Code blocks (preserve code, keep translated comments)

Usage:
    uv run translation_fixer.py fix pt
    uv run translation_fixer.py fix-file docs/pt/docs/index.md
"""

from __future__ import annotations

import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import NamedTuple

import typer
from rich.console import Console

# =============================================================================
# Configuration
# =============================================================================

REPO_ROOT = Path(__file__).parent.parent
DOCS_ROOT = REPO_ROOT / "docs"
EN_DOCS = DOCS_ROOT / "en" / "docs"

# =============================================================================
# Exceptions
# =============================================================================


class StructureMismatchError(Exception):
    """Translation structure doesn't match source."""


# =============================================================================
# Regex patterns
# =============================================================================

# Header with optional permalink: ## Title {#anchor}
HEADER_RE = re.compile(r"^(#{1,6})\s+(.+?)(\s*\{#[^}]+\})?\s*$")

# Markdown link: [text](url "title"){attrs}
MARKDOWN_LINK_RE = re.compile(
    r"(?<!!)"  # not preceded by ! (not an image)
    r"\[(?P<text>[^\]]*)\]"
    r"\((?P<url>[^)\s]+)"
    r'(?:\s+"(?P<title>[^"]*)")?'
    r"\)"
    r"(?:\{(?P<attrs>[^}]*)\})?"
)

# HTML link: <a href="...">text</a>
HTML_LINK_RE = re.compile(r"<a\s+([^>]*)>(.*?)</a>", re.DOTALL)
HTML_ATTR_RE = re.compile(r'(\w+)=(["\'])([^"\']*)\2')

# Code block fence
CODE_FENCE_RE = re.compile(r"^(`{3,4})([\w-]*)")

# Comment patterns for different languages
HASH_COMMENT_RE = re.compile(r"^(.*?)((?:\s+#|^#)\s.*)$")
SLASH_COMMENT_RE = re.compile(r"^(.*?)((?:\s+//|^//)\s.*)$")
BLOCK_COMMENT_START = re.compile(r"^\s*/\*")
BLOCK_COMMENT_END = re.compile(r"\*/\s*$")


# =============================================================================
# Data structures
# =============================================================================


class Header(NamedTuple):
    line_no: int  # 0-based
    level: str  # "##"
    title: str
    anchor: str  # "{#anchor}" or ""


class Link(NamedTuple):
    line_no: int
    start: int  # position in line
    end: int
    text: str
    url: str
    title: str | None
    attrs: str | None
    full_match: str


class HtmlLink(NamedTuple):
    line_no: int
    start: int
    end: int
    attrs: str
    text: str
    full_match: str


@dataclass
class CodeBlock:
    start_line: int  # 0-based, inclusive
    end_line: int  # 0-based, inclusive
    lang: str
    lines: list[str]


# =============================================================================
# Extraction functions
# =============================================================================


def extract_headers(lines: list[str]) -> list[Header]:
    """Extract headers, skipping those inside code blocks."""
    headers = []
    in_code = False
    fence_pattern = ""

    for i, line in enumerate(lines):
        stripped = line.lstrip()

        # Track code blocks
        if not in_code:
            if m := CODE_FENCE_RE.match(stripped):
                in_code = True
                fence_pattern = m.group(1)
                continue
        else:
            if stripped.startswith(fence_pattern) and stripped.strip() == fence_pattern:
                in_code = False
                fence_pattern = ""
            continue

        # Parse header
        if m := HEADER_RE.match(line):
            headers.append(
                Header(
                    line_no=i,
                    level=m.group(1),
                    title=m.group(2).strip(),
                    anchor=m.group(3) or "",
                )
            )

    return headers


def extract_links(lines: list[str]) -> list[Link]:
    """Extract markdown links, skipping those inside code blocks."""
    links = []
    in_code = False
    fence_pattern = ""

    for i, line in enumerate(lines):
        stripped = line.lstrip()

        # Track code blocks
        if not in_code:
            if m := CODE_FENCE_RE.match(stripped):
                in_code = True
                fence_pattern = m.group(1)
                continue
        else:
            if stripped.startswith(fence_pattern) and stripped.strip() == fence_pattern:
                in_code = False
                fence_pattern = ""
            continue

        for m in MARKDOWN_LINK_RE.finditer(line):
            links.append(
                Link(
                    line_no=i,
                    start=m.start(),
                    end=m.end(),
                    text=m.group("text"),
                    url=m.group("url"),
                    title=m.group("title"),
                    attrs=m.group("attrs"),
                    full_match=m.group(0),
                )
            )
    return links


def extract_html_links(lines: list[str]) -> list[HtmlLink]:
    """Extract HTML <a> links, skipping those inside code blocks."""
    links = []
    in_code = False
    fence_pattern = ""

    for i, line in enumerate(lines):
        stripped = line.lstrip()

        # Track code blocks
        if not in_code:
            if m := CODE_FENCE_RE.match(stripped):
                in_code = True
                fence_pattern = m.group(1)
                continue
        else:
            if stripped.startswith(fence_pattern) and stripped.strip() == fence_pattern:
                in_code = False
                fence_pattern = ""
            continue

        for m in HTML_LINK_RE.finditer(line):
            links.append(
                HtmlLink(
                    line_no=i,
                    start=m.start(),
                    end=m.end(),
                    attrs=m.group(1),
                    text=m.group(2),
                    full_match=m.group(0),
                )
            )
    return links


def extract_code_blocks(lines: list[str]) -> list[CodeBlock]:
    """Extract code blocks with their content."""
    blocks = []
    in_code = False
    fence_pattern = ""
    start_line = 0
    lang = ""
    block_lines: list[str] = []

    for i, line in enumerate(lines):
        stripped = line.lstrip()

        if not in_code:
            if m := CODE_FENCE_RE.match(stripped):
                in_code = True
                fence_pattern = m.group(1)
                lang = m.group(2)
                start_line = i
                block_lines = [line]
        else:
            block_lines.append(line)
            if stripped.startswith(fence_pattern) and stripped.strip() == fence_pattern:
                blocks.append(
                    CodeBlock(
                        start_line=start_line,
                        end_line=i,
                        lang=lang,
                        lines=block_lines,
                    )
                )
                in_code = False
                fence_pattern = ""
                block_lines = []

    return blocks


# =============================================================================
# Replacement functions
# =============================================================================


def fix_headers(trans_lines: list[str], source_headers: list[Header]) -> list[str]:
    """Replace header anchors with those from source."""
    trans_headers = extract_headers(trans_lines)

    if len(trans_headers) != len(source_headers):
        raise StructureMismatchError(
            f"Header count mismatch: translation has {len(trans_headers)}, "
            f"source has {len(source_headers)}"
        )

    result = trans_lines.copy()
    for trans_h, src_h in zip(trans_headers, source_headers):
        if trans_h.level != src_h.level:
            raise StructureMismatchError(
                f"Header level mismatch at line {trans_h.line_no + 1}: "
                f"'{trans_h.level}' vs '{src_h.level}'"
            )
        # Reconstruct header with source anchor
        new_line = f"{trans_h.level} {trans_h.title}"
        if src_h.anchor:
            new_line += src_h.anchor
        result[trans_h.line_no] = new_line

    return result


def fix_links(trans_lines: list[str], source_links: list[Link]) -> list[str]:
    """Replace link URLs/attrs with source, keeping translated text."""
    trans_links = extract_links(trans_lines)

    if len(trans_links) != len(source_links):
        raise StructureMismatchError(
            f"Link count mismatch: translation has {len(trans_links)}, "
            f"source has {len(source_links)}"
        )

    result = trans_lines.copy()

    # Process in reverse order to preserve positions
    for trans_l, src_l in reversed(list(zip(trans_links, source_links))):
        # Build new link with source URL but translated text
        new_link = f"[{trans_l.text}]({src_l.url}"
        if trans_l.title:
            new_link += f' "{trans_l.title}"'
        new_link += ")"
        if src_l.attrs:
            new_link += f"{{{src_l.attrs}}}"

        line = result[trans_l.line_no]
        result[trans_l.line_no] = line[: trans_l.start] + new_link + line[trans_l.end :]

    return result


def fix_html_links(trans_lines: list[str], source_links: list[HtmlLink]) -> list[str]:
    """Replace HTML link attributes with source, keeping translated text."""
    trans_links = extract_html_links(trans_lines)

    if len(trans_links) != len(source_links):
        raise StructureMismatchError(
            f"HTML link count mismatch: translation has {len(trans_links)}, "
            f"source has {len(source_links)}"
        )

    result = trans_lines.copy()

    for trans_l, src_l in reversed(list(zip(trans_links, source_links))):
        new_link = f"<a {src_l.attrs}>{trans_l.text}</a>"
        line = result[trans_l.line_no]
        result[trans_l.line_no] = line[: trans_l.start] + new_link + line[trans_l.end :]

    return result


def fix_code_block(trans_block: CodeBlock, source_block: CodeBlock) -> list[str]:
    """Fix a code block: use source code, keep translated comments."""
    if trans_block.lang != source_block.lang:
        raise StructureMismatchError(
            f"Code block language mismatch at line {trans_block.start_line + 1}: "
            f"'{trans_block.lang}' vs '{source_block.lang}'"
        )

    if len(trans_block.lines) != len(source_block.lines):
        raise StructureMismatchError(
            f"Code block line count mismatch at line {trans_block.start_line + 1}: "
            f"{len(trans_block.lines)} vs {len(source_block.lines)}"
        )

    lang = trans_block.lang.lower()

    # Skip mermaid diagrams - they're structural
    if lang == "mermaid":
        return source_block.lines.copy()

    # Determine comment style based on language
    hash_langs = {"python", "py", "sh", "bash", "dockerfile", "yaml", "yml", "toml"}
    slash_langs = {"console", "json"}
    mixed_langs = {"nextflow", "groovy", "nf", "java", "kotlin"}

    result = []
    in_block_comment = False

    for trans_line, src_line in zip(trans_block.lines, source_block.lines):
        # First and last lines are fences - use source
        if trans_line.lstrip().startswith("```"):
            result.append(src_line)
            continue

        # Handle block comments for Groovy-like languages
        if lang in mixed_langs:
            if BLOCK_COMMENT_START.match(trans_line):
                in_block_comment = True
                result.append(trans_line)  # Keep translated block comment
                continue
            if in_block_comment:
                if BLOCK_COMMENT_END.search(trans_line):
                    in_block_comment = False
                result.append(trans_line)  # Keep translated block comment
                continue

        # Extract and preserve translated comments
        trans_comment = None
        if lang in hash_langs or lang in mixed_langs:
            if m := HASH_COMMENT_RE.match(trans_line):
                trans_comment = m.group(2)
        elif lang in slash_langs:
            if m := SLASH_COMMENT_RE.match(trans_line):
                trans_comment = m.group(2)

        # Use source line, but replace comment with translated version if present
        if trans_comment:
            # Check if source has a comment we can replace
            src_comment = None
            if lang in hash_langs or lang in mixed_langs:
                if m := HASH_COMMENT_RE.match(src_line):
                    src_comment = m.group(2)
            elif lang in slash_langs:
                if m := SLASH_COMMENT_RE.match(src_line):
                    src_comment = m.group(2)

            if src_comment:
                result.append(src_line.replace(src_comment, trans_comment))
            else:
                result.append(src_line)
        else:
            result.append(src_line)

    return result


def fix_code_blocks(
    trans_lines: list[str], source_blocks: list[CodeBlock]
) -> list[str]:
    """Fix all code blocks in translation."""
    trans_blocks = extract_code_blocks(trans_lines)

    if len(trans_blocks) != len(source_blocks):
        raise StructureMismatchError(
            f"Code block count mismatch: translation has {len(trans_blocks)}, "
            f"source has {len(source_blocks)}"
        )

    result = trans_lines.copy()

    for trans_b, src_b in zip(trans_blocks, source_blocks):
        fixed = fix_code_block(trans_b, src_b)
        for i, line in enumerate(fixed):
            result[trans_b.start_line + i] = line

    return result


# =============================================================================
# Main fix function
# =============================================================================


def fix_translation(translation: list[str], source: list[str]) -> list[str]:
    """
    Fix a translation against its English source.

    Preserves URLs, anchors, and code from source while keeping translated text.
    """
    result = translation.copy()

    # Fix headers (preserve anchors)
    source_headers = extract_headers(source)
    result = fix_headers(result, source_headers)

    # Fix markdown links (preserve URLs)
    source_links = extract_links(source)
    result = fix_links(result, source_links)

    # Fix HTML links (preserve attributes)
    source_html = extract_html_links(source)
    result = fix_html_links(result, source_html)

    # Fix code blocks (preserve code, keep translated comments)
    source_blocks = extract_code_blocks(source)
    result = fix_code_blocks(result, source_blocks)

    return result


# =============================================================================
# CLI
# =============================================================================

app = typer.Typer(help="Post-process translations to fix common LLM mistakes")
console = Console()


def lang_to_en_path(lang_path: Path, lang: str) -> Path:
    """Convert language path to English path."""
    rel = lang_path.relative_to(DOCS_ROOT / lang / "docs")
    return EN_DOCS / rel


def fix_file(path: Path, lang: str) -> bool:
    """Fix a single file. Returns True if changes were made."""
    en_path = lang_to_en_path(path, lang)
    if not en_path.exists():
        raise FileNotFoundError(f"English source not found: {en_path}")

    original = path.read_text(encoding="utf-8")
    source = en_path.read_text(encoding="utf-8")

    fixed = fix_translation(original.splitlines(), source.splitlines())
    fixed_text = "\n".join(fixed) + "\n"

    if fixed_text != original:
        path.write_text(fixed_text, encoding="utf-8", newline="\n")
        return True
    return False


@app.command()
def fix(lang: str = typer.Argument(..., help="Language code")):
    """Fix all translations for a language."""
    lang_docs = DOCS_ROOT / lang / "docs"
    if not lang_docs.exists():
        raise typer.Exit(f"Language docs not found: {lang_docs}")

    files = list(lang_docs.rglob("*.md"))
    console.print(f"Processing {len(files)} files for {lang}...")

    fixed = 0
    for path in files:
        try:
            if fix_file(path, lang):
                console.print(f"  [green]Fixed:[/green] {path.name}")
                fixed += 1
            else:
                console.print(f"  [dim]OK:[/dim] {path.name}")
        except (StructureMismatchError, FileNotFoundError) as e:
            console.print(f"  [yellow]Skip:[/yellow] {path.name}: {e}")

    console.print(f"\n[green]Fixed {fixed} files[/green]")


@app.command("fix-file")
def fix_file_cmd(paths: list[Path] = typer.Argument(..., help="Files to fix")):
    """Fix specific translation files."""
    for path in paths:
        if not path.exists():
            raise typer.Exit(f"File not found: {path}")

        # Determine language from path
        try:
            lang = path.relative_to(DOCS_ROOT).parts[0]
        except ValueError:
            raise typer.Exit(f"Path not in docs/: {path}")

        if lang == "en":
            console.print(f"[yellow]Skipping English file:[/yellow] {path}")
            continue

        try:
            if fix_file(path, lang):
                console.print(f"[green]Fixed:[/green] {path}")
            else:
                console.print(f"[dim]OK:[/dim] {path}")
        except (StructureMismatchError, FileNotFoundError) as e:
            console.print(f"[yellow]Skip:[/yellow] {path}: {e}")


if __name__ == "__main__":
    if not DOCS_ROOT.exists():
        sys.exit("Error: Must run from repository root")
    app()
