#!/usr/bin/env python3
# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "anthropic>=0.40.0",
#     "gitpython>=3.1.0",
#     "typer>=0.12.0",
#     "pyyaml>=6.0",
#     "rich>=13.0.0",
#     "pygithub>=2.0.0",
# ]
# ///
"""
LLM-powered translation CLI for Nextflow training docs.

Usage:
    uv run translate.py sync pt --dry-run
    uv run translate.py sync pt
    uv run translate.py translate docs/en/docs/index.md --lang pt
"""

from __future__ import annotations

import json
import os
import re
import secrets
import subprocess
import sys
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import NamedTuple

import time

import anthropic
import git
import typer
import yaml
from rich.console import Console


# =============================================================================
# Configuration
# =============================================================================

REPO_ROOT = Path(__file__).parent.parent
DOCS_ROOT = REPO_ROOT / "docs"
EN_DOCS = DOCS_ROOT / "en" / "docs"
SCRIPTS_DIR = REPO_ROOT / "_scripts"

MODEL = "claude-sonnet-4-5"
MAX_TOKENS = 16384
REQUEST_TIMEOUT = 600.0  # 10 minute timeout per API request

# Retry configuration for transient API errors
MAX_RETRIES = 8
BASE_DELAY = 2.0  # seconds, doubles each retry: 2, 4, 8, 16, 32, 64, 128, 256

PRIORITY_DIRS = ["hello_nextflow", "hello_nf-core", "nf4_science", "envsetup"]

# Comment styles by language
HASH_COMMENT_LANGS = {"python", "py", "sh", "bash", "dockerfile", "yaml", "yml", "toml"}
SLASH_COMMENT_LANGS = {"console"}  # Note: json has no comments
MIXED_COMMENT_LANGS = {"nextflow", "groovy", "nf", "java", "kotlin"}  # # and /* */


# =============================================================================
# Exceptions
# =============================================================================


class TranslationError(Exception):
    """Base exception for translation errors."""


class ConfigError(TranslationError):
    """Configuration or setup error."""


class StructureMismatchError(TranslationError):
    """Translation structure doesn't match source."""


# =============================================================================
# Data structures
# =============================================================================


@dataclass
class TranslationFile:
    """A file that needs translation."""

    en_path: Path
    lang_path: Path
    language: str

    @property
    def exists(self) -> bool:
        return self.lang_path.exists()

    @property
    def relative_path(self) -> Path:
        return self.en_path.relative_to(EN_DOCS)


class Header(NamedTuple):
    line_no: int
    level: str
    title: str
    anchor: str


class Link(NamedTuple):
    line_no: int
    start: int
    end: int
    text: str
    url: str
    title: str | None
    attrs: str | None


class HtmlLink(NamedTuple):
    line_no: int
    start: int
    end: int
    attrs: str
    text: str


@dataclass
class CodeBlock:
    start_line: int
    end_line: int
    lang: str
    lines: list[str]


# =============================================================================
# Path utilities
# =============================================================================


@lru_cache
def get_languages() -> dict[str, str]:
    """Load language code -> name mapping."""
    path = DOCS_ROOT / "language_names.yml"
    if not path.exists():
        raise ConfigError(f"Language names file not found: {path}")
    return yaml.safe_load(path.read_text(encoding="utf-8"))


def get_translation_languages() -> list[str]:
    """Get all translation language codes (excludes English)."""
    return [
        d.name
        for d in DOCS_ROOT.iterdir()
        if d.is_dir() and (d / "mkdocs.yml").exists() and d.name != "en"
    ]


def en_to_lang_path(en_path: Path, lang: str) -> Path:
    """Convert English doc path to equivalent path in target language."""
    return DOCS_ROOT / lang / "docs" / en_path.relative_to(EN_DOCS)


def lang_to_en_path(lang_path: Path, lang: str) -> Path:
    """Convert language doc path to equivalent English path."""
    return EN_DOCS / lang_path.relative_to(DOCS_ROOT / lang / "docs")


def iter_en_docs() -> list[Path]:
    """List all English docs in priority order."""
    paths: list[Path] = []
    seen: set[Path] = set()

    def add(p: Path) -> None:
        if p not in seen:
            paths.append(p)
            seen.add(p)

    # Root files first
    for p in sorted(EN_DOCS.glob("*.md")):
        add(p)

    # Priority directories
    for dir_name in PRIORITY_DIRS:
        dir_path = EN_DOCS / dir_name
        if dir_path.exists():
            for p in sorted(dir_path.rglob("*.md")):
                add(p)

    # Remaining
    for p in sorted(EN_DOCS.rglob("*.md")):
        add(p)

    return paths


# =============================================================================
# Prompt loading
# =============================================================================


@lru_cache
def get_general_prompt() -> str:
    """Load the general translation prompt."""
    path = SCRIPTS_DIR / "general-llm-prompt.md"
    if not path.exists():
        raise ConfigError(f"General prompt not found: {path}")
    return path.read_text(encoding="utf-8")


@lru_cache
def get_lang_prompt(lang: str) -> str:
    """Load language-specific prompt."""
    path = DOCS_ROOT / lang / "llm-prompt.md"
    if not path.exists():
        raise ConfigError(f"Language prompt not found: {path}")
    return path.read_text(encoding="utf-8")


def build_translation_prompt(
    lang: str, lang_name: str, en_content: str, existing: str | None = None
) -> str:
    """Build the full prompt for translation."""
    parts = [get_general_prompt(), get_lang_prompt(lang)]

    if existing:
        parts.append(
            "## Existing Translation\n"
            "Update minimally: add new content, remove deleted content, "
            "fix guideline violations, preserve correct lines exactly.\n\n"
            f"Previous translation:\n%%%\n{existing}%%%"
        )

    parts.append(
        f"## Task\nTranslate to {lang} ({lang_name}).\n\n"
        f"Original content:\n%%%\n{en_content}%%%"
    )

    return "\n\n".join(parts)


# =============================================================================
# Git utilities
# =============================================================================


@lru_cache
def _get_repo() -> git.Repo:
    return git.Repo(REPO_ROOT)


def get_file_commit_time(path: Path) -> int | None:
    """Get timestamp of last commit affecting a file."""
    try:
        commits = list(_get_repo().iter_commits(paths=str(path), max_count=1))
        return commits[0].committed_date if commits else None
    except git.GitCommandError:
        return None


def get_outdated_files(lang: str) -> list[TranslationFile]:
    """Find translations older than their English source."""
    outdated = []
    for en_path in iter_en_docs():
        lang_path = en_to_lang_path(en_path, lang)
        if not lang_path.exists():
            continue

        en_time = get_file_commit_time(en_path)
        lang_time = get_file_commit_time(lang_path)

        if en_time and lang_time and lang_time < en_time:
            outdated.append(TranslationFile(en_path, lang_path, lang))

    return outdated


def get_missing_files(lang: str) -> list[TranslationFile]:
    """Find English files without translations."""
    return [
        TranslationFile(en_path, en_to_lang_path(en_path, lang), lang)
        for en_path in iter_en_docs()
        if not en_to_lang_path(en_path, lang).exists()
    ]


def get_orphaned_files(lang: str) -> list[Path]:
    """Find translation files without English source."""
    lang_docs = DOCS_ROOT / lang / "docs"
    if not lang_docs.exists():
        return []
    return [p for p in lang_docs.rglob("*.md") if not lang_to_en_path(p, lang).exists()]


# =============================================================================
# Translation API
# =============================================================================


def check_api_key() -> None:
    """Verify ANTHROPIC_API_KEY is set."""
    if not os.environ.get("ANTHROPIC_API_KEY"):
        raise ConfigError(
            "ANTHROPIC_API_KEY not set. Get key at https://console.anthropic.com/settings/keys"
        )


def call_claude(prompt: str, console: Console | None = None) -> str:
    """Call Claude API with retry logic and exponential backoff.

    Retries on transient errors (connection, timeout, rate limit, server errors)
    with exponential backoff: 2s, 4s, 8s, 16s, 32s, 64s, 128s, 256s (~8.5 min total).
    """
    client = anthropic.Anthropic()

    for attempt in range(MAX_RETRIES + 1):
        try:
            message = client.messages.create(
                model=MODEL,
                max_tokens=MAX_TOKENS,
                timeout=REQUEST_TIMEOUT,
                messages=[{"role": "user", "content": prompt}],
            )
            return message.content[0].text.strip()
        except (
            anthropic.APIConnectionError,
            anthropic.APITimeoutError,
            anthropic.RateLimitError,
            anthropic.InternalServerError,
        ) as e:
            if attempt == MAX_RETRIES:
                raise
            delay = BASE_DELAY * (2**attempt)
            if console:
                console.print(
                    f"    [yellow]Retry {attempt + 1}/{MAX_RETRIES} "
                    f"after {delay:.0f}s: {type(e).__name__}[/yellow]"
                )
            time.sleep(delay)

    raise RuntimeError("Unreachable")


def translate_file(tf: TranslationFile, console: Console) -> None:
    """Translate a single file."""
    langs = get_languages()
    if tf.language not in langs:
        raise ConfigError(f"Unknown language: {tf.language}")

    en_content = tf.en_path.read_text(encoding="utf-8")
    existing = tf.lang_path.read_text(encoding="utf-8") if tf.exists else None

    action = "[yellow]Updating[/yellow]" if existing else "[green]Translating[/green]"
    console.print(f"  {action} [magenta]{tf.relative_path}[/magenta]")

    prompt = build_translation_prompt(
        tf.language, langs[tf.language], en_content, existing
    )
    result = call_claude(prompt, console)

    tf.lang_path.parent.mkdir(parents=True, exist_ok=True)
    tf.lang_path.write_text(f"{result}\n", encoding="utf-8", newline="\n")


# =============================================================================
# Post-processing
# =============================================================================

HEADER_RE = re.compile(r"^(#{1,6})\s+(.+?)(\s*\{#[^}]+\})?\s*$")
MARKDOWN_LINK_RE = re.compile(
    r"(?<!!)\[(?P<text>[^\]]*)\]\((?P<url>[^)\s]+)(?:\s+\"(?P<title>[^\"]*)\")?\)(?:\{(?P<attrs>[^}]*)\})?"
)
HTML_LINK_RE = re.compile(r"<a\s+([^>]*)>(.*?)</a>", re.DOTALL)
CODE_FENCE_RE = re.compile(r"^(`{3,4})([\w-]*)")
HASH_COMMENT_RE = re.compile(r"^(.*?)((?:\s+#|^#)\s.*)$")
SLASH_COMMENT_RE = re.compile(r"^(.*?)((?:\s+//|^//)\s.*)$")
BLOCK_COMMENT_START_RE = re.compile(r"^\s*/\*")
BLOCK_COMMENT_END_RE = re.compile(r"\*/\s*$")


def _iter_lines_outside_code(lines: list[str]):
    """Yield (line_no, line) for lines outside code blocks."""
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


def extract_headers(lines: list[str]) -> list[Header]:
    """Extract headers outside code blocks."""
    return [
        Header(i, m.group(1), m.group(2).strip(), m.group(3) or "")
        for i, line in _iter_lines_outside_code(lines)
        if (m := HEADER_RE.match(line))
    ]


def extract_links(lines: list[str]) -> list[Link]:
    """Extract markdown links outside code blocks."""
    links = []
    for i, line in _iter_lines_outside_code(lines):
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
    """Extract HTML <a> links outside code blocks."""
    links = []
    for i, line in _iter_lines_outside_code(lines):
        for m in HTML_LINK_RE.finditer(line):
            links.append(HtmlLink(i, m.start(), m.end(), m.group(1), m.group(2)))
    return links


def extract_code_blocks(lines: list[str]) -> list[CodeBlock]:
    """Extract code blocks with content."""
    blocks = []
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


def _extract_comment(line: str, lang: str) -> str | None:
    """Extract comment from a line based on language."""
    if lang in HASH_COMMENT_LANGS or lang in MIXED_COMMENT_LANGS:
        if m := HASH_COMMENT_RE.match(line):
            return m.group(2)
    if lang in SLASH_COMMENT_LANGS:
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
    """Replace link URLs/attrs with source, keeping translated text."""
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
    """Fix code block: use source code, preserve translated comments."""
    if trans_block.lang != src_block.lang:
        raise StructureMismatchError(
            f"Code lang mismatch at line {trans_block.start_line}"
        )
    if len(trans_block.lines) != len(src_block.lines):
        raise StructureMismatchError(
            f"Code lines mismatch at line {trans_block.start_line}"
        )

    lang = trans_block.lang.lower()

    # Mermaid diagrams are structural, use source entirely
    if lang == "mermaid":
        return src_block.lines.copy()

    result = []
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


def fix_translation(translation: list[str], source: list[str]) -> list[str]:
    """Fix translation against English source."""
    result = translation.copy()
    result = fix_headers(result, extract_headers(source))
    result = fix_links(result, extract_links(source))
    result = fix_html_links(result, extract_html_links(source))
    result = fix_code_blocks(result, extract_code_blocks(source))
    return result


def post_process_file(lang_path: Path, lang: str) -> bool:
    """Post-process a translation. Returns True if changes were made."""
    en_path = lang_to_en_path(lang_path, lang)
    if not en_path.exists():
        return False

    original = lang_path.read_text(encoding="utf-8")
    source = en_path.read_text(encoding="utf-8")

    try:
        fixed = fix_translation(original.splitlines(), source.splitlines())
    except StructureMismatchError:
        return False

    fixed_text = "\n".join(fixed) + "\n"
    if fixed_text != original:
        lang_path.write_text(fixed_text, encoding="utf-8", newline="\n")
        return True
    return False


# =============================================================================
# CLI
# =============================================================================

app = typer.Typer(help="Translation CLI for Nextflow training docs")


def _resolve_en_path(path: Path) -> Path:
    """Resolve path to English doc file."""
    # Absolute path
    if path.is_absolute():
        if path.exists():
            return path
        raise ConfigError(f"File not found: {path}")

    # Relative to EN_DOCS
    if (EN_DOCS / path).exists():
        return EN_DOCS / path

    # Relative to cwd
    if path.exists():
        return path.resolve()

    raise ConfigError(f"File not found: {path}")


@app.command()
def translate(
    path: Path = typer.Argument(..., help="English file to translate"),
    lang: str = typer.Option(..., "--lang", "-l", help="Target language code"),
):
    """Translate a single file."""
    check_api_key()
    console = Console(force_terminal=True if os.getenv("GITHUB_ACTIONS") else None)

    if lang == "en":
        raise typer.Exit("Cannot translate to English")

    en_path = _resolve_en_path(path)
    tf = TranslationFile(en_path, en_to_lang_path(en_path, lang), lang)

    translate_file(tf, console)
    post_process_file(tf.lang_path, lang)
    console.print(f"[green]Done:[/green] {tf.lang_path}")


@app.command()
def sync(
    lang: str = typer.Argument(..., help="Language code"),
    include: str | None = typer.Option(None, "--include", "-i", help="Filter pattern"),
    dry_run: bool = typer.Option(False, "--dry-run", "-n", help="Preview only"),
):
    """Sync translations: update outdated, add missing, remove orphaned."""
    console = Console(force_terminal=True if os.getenv("GITHUB_ACTIONS") else None)

    # Gather work
    orphaned = get_orphaned_files(lang)
    outdated = get_outdated_files(lang)
    missing = get_missing_files(lang)

    if include:
        outdated = [tf for tf in outdated if include in str(tf.en_path)]
        missing = [tf for tf in missing if include in str(tf.en_path)]

    # Dry run: just print
    if dry_run:
        if orphaned:
            console.print(f"[red]Would remove {len(orphaned)} orphaned:[/red]")
            for p in orphaned:
                console.print(f"  {p.relative_to(REPO_ROOT)}")
        if outdated:
            console.print(f"[yellow]Would update {len(outdated)} outdated:[/yellow]")
            for tf in outdated:
                console.print(f"  {tf.relative_path}")
        if missing:
            console.print(f"[blue]Would add {len(missing)} missing:[/blue]")
            for tf in missing:
                console.print(f"  {tf.relative_path}")
        if not (orphaned or outdated or missing):
            console.print(f"[green]Nothing to do for {lang}[/green]")
        return

    check_api_key()

    # Remove orphaned
    if orphaned:
        console.print(f"[bold]Removing {len(orphaned)} orphaned...[/bold]")
        for p in orphaned:
            p.unlink()
            console.print(f"  [red]Removed:[/red] {p.relative_to(REPO_ROOT)}")

    # Update outdated
    if outdated:
        console.print(f"[bold]Updating {len(outdated)} outdated...[/bold]")
        for i, tf in enumerate(outdated, 1):
            console.print(f"[{i}/{len(outdated)}]", end="", style="blue")
            translate_file(tf, console)
            post_process_file(tf.lang_path, lang)

    # Add missing
    if missing:
        console.print(f"[bold]Adding {len(missing)} missing...[/bold]")
        for i, tf in enumerate(missing, 1):
            console.print(f"[{i}/{len(missing)}]", end="")
            translate_file(tf, console)
            post_process_file(tf.lang_path, lang)

    console.print("[green]:white_check_mark: Sync complete[/green]")


# =============================================================================
# CI Commands
# =============================================================================


@app.command("ci-detect")
def ci_detect(language: str | None = typer.Option(None, "--language")):
    """Detect languages needing sync (GitHub Actions output)."""
    console = Console(
        stderr=True, force_terminal=True if os.getenv("GITHUB_ACTIONS") else None
    )
    all_langs = get_translation_languages()

    if language:
        if language not in all_langs:
            raise typer.Exit(f"Unknown language: {language}")
        langs = [language]
    else:
        langs = all_langs

    # Check which have work and collect details
    need_sync = []
    for lang in langs:
        missing = get_missing_files(lang)
        outdated = get_outdated_files(lang)
        orphaned = get_orphaned_files(lang)

        if missing or outdated or orphaned:
            need_sync.append(lang)
            console.print(f"[bold cyan]{lang}[/bold cyan]:")
            if outdated:
                console.print(f"  [yellow]Outdated:[/yellow] {len(outdated)}")
                for f in outdated:
                    console.print(f"    {f.relative_path}")
            if missing:
                console.print(f"  [green]Missing:[/green] {len(missing)}")
                for f in missing:
                    console.print(f"    {f.relative_path}")
            if orphaned:
                console.print(f"  [red]Orphaned:[/red] {len(orphaned)}")
                for f in orphaned:
                    console.print(f"    {f.relative_path}")

    print(f"languages={json.dumps(need_sync)}")
    print(f"has_work={'true' if need_sync else 'false'}")


@app.command("ci-run")
def ci_run(
    languages: str = typer.Option("[]", "--languages"),
    dry_run: bool = typer.Option(False, "--dry-run"),
):
    """Run sync for multiple languages (CI)."""
    console = Console(force_terminal=True if os.getenv("GITHUB_ACTIONS") else None)
    for lang in json.loads(languages):
        console.rule(f"[bold]{lang}[/bold]")
        sync(lang, include=None, dry_run=dry_run)


@app.command("ci-pr")
def ci_pr(
    github_token: str = typer.Option(..., envvar="GITHUB_TOKEN"),
    github_repository: str = typer.Option(..., envvar="GITHUB_REPOSITORY"),
):
    """Create PR with translation changes (CI)."""
    from github import Github

    repo = _get_repo()
    if not repo.is_dirty(untracked_files=True):
        print("No changes to commit")
        return

    subprocess.run(
        ["git", "config", "user.name", "github-actions[bot]"], check=True, cwd=REPO_ROOT
    )
    subprocess.run(
        ["git", "config", "user.email", "github-actions[bot]@users.noreply.github.com"],
        check=True,
        cwd=REPO_ROOT,
    )

    branch = f"translate-{secrets.token_hex(4)}"
    subprocess.run(["git", "checkout", "-b", branch], check=True, cwd=REPO_ROOT)
    subprocess.run(["git", "add", "docs/"], check=True, cwd=REPO_ROOT)
    subprocess.run(
        ["git", "commit", "-m", "Update translations"], check=True, cwd=REPO_ROOT
    )
    subprocess.run(["git", "push", "origin", branch], check=True, cwd=REPO_ROOT)

    gh = Github(github_token)
    pr = gh.get_repo(github_repository).create_pull(
        title="Update translations",
        body="Automated translation update.\n\nTo improve translations, edit prompt files:\n"
        "- `_scripts/general-llm-prompt.md`\n- `docs/*/llm-prompt.md`",
        base="master",
        head=branch,
    )
    print(f"Created PR: {pr.html_url}")


if __name__ == "__main__":
    if not DOCS_ROOT.exists():
        sys.exit("Error: Must run from repository root")
    app()
