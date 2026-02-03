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

# Translation priority order
PRIORITY_DIRS = ["hello_nextflow", "hello_nf-core", "nf4_science", "envsetup"]


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
# Core data structures
# =============================================================================


@dataclass
class TranslationFile:
    """Represents a file that needs translation."""

    en_path: Path
    lang_path: Path
    language: str

    @property
    def exists(self) -> bool:
        return self.lang_path.exists()

    @property
    def relative_path(self) -> Path:
        return self.en_path.relative_to(EN_DOCS)


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
    """Get all languages except English."""
    return [
        d.name
        for d in DOCS_ROOT.iterdir()
        if d.is_dir() and (d / "mkdocs.yml").exists() and d.name != "en"
    ]


def en_to_lang_path(en_path: Path, lang: str) -> Path:
    """Convert English path to language-specific path."""
    rel = en_path.relative_to(EN_DOCS)
    return DOCS_ROOT / lang / "docs" / rel


def lang_to_en_path(lang_path: Path, lang: str) -> Path:
    """Convert language path to English path."""
    rel = lang_path.relative_to(DOCS_ROOT / lang / "docs")
    return EN_DOCS / rel


def resolve_path(path: Path) -> Path:
    """Resolve a path that may be relative to various locations."""
    for base in [Path.cwd(), REPO_ROOT, EN_DOCS]:
        candidate = base / path if not path.is_absolute() else path
        if candidate.exists():
            return candidate.resolve()
    raise ConfigError(f"File not found: {path}")


def iter_en_docs() -> list[Path]:
    """Iterate English docs in priority order."""
    paths: list[Path] = []

    # Root files first
    paths.extend(sorted(EN_DOCS.glob("*.md")))

    # Priority directories
    for dir_name in PRIORITY_DIRS:
        dir_path = EN_DOCS / dir_name
        if dir_path.exists():
            paths.extend(sorted(dir_path.rglob("*.md")))

    # Remaining directories
    seen = {p for p in paths}
    for path in sorted(EN_DOCS.rglob("*.md")):
        if path not in seen:
            paths.append(path)

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
    lang: str, lang_name: str, en_content: str, existing_translation: str | None = None
) -> str:
    """Build the full prompt for translation."""
    parts = [get_general_prompt(), get_lang_prompt(lang)]

    if existing_translation:
        parts.extend(
            [
                "## Existing Translation",
                "Update the existing translation minimally:",
                "- Add new content from English source",
                "- Remove deleted content",
                "- Update changed content",
                "- Fix guideline violations",
                "- Preserve correct lines exactly (minimize diff)",
                "",
                f"Previous translation:\n%%%\n{existing_translation}%%%",
            ]
        )

    parts.extend(
        [
            "## Task",
            f"Translate to {lang} ({lang_name}).",
            f"Original content:\n%%%\n{en_content}%%%",
        ]
    )

    return "\n\n".join(parts)


# =============================================================================
# Git utilities
# =============================================================================


def get_file_commit_time(repo: git.Repo, path: Path) -> int | None:
    """Get the timestamp of the last commit affecting a file."""
    try:
        commits = list(repo.iter_commits(paths=str(path), max_count=1))
        return commits[0].committed_date if commits else None
    except git.GitCommandError:
        return None


def get_outdated_files(lang: str) -> list[TranslationFile]:
    """Find translations older than their English source."""
    repo = git.Repo(REPO_ROOT)
    outdated = []

    for en_path in iter_en_docs():
        tf = TranslationFile(
            en_path=en_path,
            lang_path=en_to_lang_path(en_path, lang),
            language=lang,
        )
        if not tf.exists:
            continue

        en_time = get_file_commit_time(repo, en_path)
        lang_time = get_file_commit_time(repo, tf.lang_path)

        if en_time and lang_time and lang_time < en_time:
            outdated.append(tf)

    return outdated


def get_missing_files(lang: str) -> list[TranslationFile]:
    """Find English files without translations."""
    return [
        TranslationFile(
            en_path=en_path,
            lang_path=en_to_lang_path(en_path, lang),
            language=lang,
        )
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
            "ANTHROPIC_API_KEY not set. "
            "Get your key at https://console.anthropic.com/settings/keys"
        )


def call_claude(prompt: str) -> str:
    """Call Claude API and return the response text."""
    client = anthropic.Anthropic()
    message = client.messages.create(
        model=MODEL,
        max_tokens=MAX_TOKENS,
        messages=[{"role": "user", "content": prompt}],
    )
    return message.content[0].text.strip()


def translate_file(tf: TranslationFile, console: Console) -> None:
    """Translate a single file."""
    langs = get_languages()
    if tf.language not in langs:
        raise ConfigError(f"Unknown language: {tf.language}")

    lang_name = langs[tf.language]
    en_content = tf.en_path.read_text(encoding="utf-8")
    existing = tf.lang_path.read_text(encoding="utf-8") if tf.exists else None

    action = "Updating" if existing else "Translating"
    console.print(f"  {action} [cyan]{tf.relative_path}[/cyan]")

    prompt = build_translation_prompt(tf.language, lang_name, en_content, existing)
    result = call_claude(prompt)

    tf.lang_path.parent.mkdir(parents=True, exist_ok=True)
    tf.lang_path.write_text(f"{result}\n", encoding="utf-8", newline="\n")


# =============================================================================
# Post-processing (fix common LLM mistakes)
# =============================================================================

# Regex patterns
HEADER_RE = re.compile(r"^(#{1,6})\s+(.+?)(\s*\{#[^}]+\})?\s*$")
MARKDOWN_LINK_RE = re.compile(
    r"(?<!!)"  # not preceded by ! (not an image)
    r"\[(?P<text>[^\]]*)\]"
    r"\((?P<url>[^)\s]+)"
    r'(?:\s+"(?P<title>[^"]*)")?'
    r"\)"
    r"(?:\{(?P<attrs>[^}]*)\})?"
)
HTML_LINK_RE = re.compile(r"<a\s+([^>]*)>(.*?)</a>", re.DOTALL)
CODE_FENCE_RE = re.compile(r"^(`{3,4})([\w-]*)")
HASH_COMMENT_RE = re.compile(r"^(.*?)((?:\s+#|^#)\s.*)$")
SLASH_COMMENT_RE = re.compile(r"^(.*?)((?:\s+//|^//)\s.*)$")
BLOCK_COMMENT_START = re.compile(r"^\s*/\*")
BLOCK_COMMENT_END = re.compile(r"\*/\s*$")


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


def _in_code_block(lines: list[str]):
    """Generator that yields (line_no, line, in_code) for each line."""
    in_code = False
    fence = ""
    for i, line in enumerate(lines):
        stripped = line.lstrip()
        if not in_code:
            if m := CODE_FENCE_RE.match(stripped):
                in_code = True
                fence = m.group(1)
                yield i, line, True
                continue
        else:
            if stripped.startswith(fence) and stripped.strip() == fence:
                in_code = False
                fence = ""
                yield i, line, True
                continue
        yield i, line, in_code


def extract_headers(lines: list[str]) -> list[Header]:
    """Extract headers, skipping those inside code blocks."""
    headers = []
    for i, line, in_code in _in_code_block(lines):
        if in_code:
            continue
        if m := HEADER_RE.match(line):
            headers.append(Header(i, m.group(1), m.group(2).strip(), m.group(3) or ""))
    return headers


def extract_links(lines: list[str]) -> list[Link]:
    """Extract markdown links, skipping those inside code blocks."""
    links = []
    for i, line, in_code in _in_code_block(lines):
        if in_code:
            continue
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
    """Extract HTML <a> links, skipping those inside code blocks."""
    links = []
    for i, line, in_code in _in_code_block(lines):
        if in_code:
            continue
        for m in HTML_LINK_RE.finditer(line):
            links.append(HtmlLink(i, m.start(), m.end(), m.group(1), m.group(2)))
    return links


def extract_code_blocks(lines: list[str]) -> list[CodeBlock]:
    """Extract code blocks with their content."""
    blocks = []
    in_code = False
    fence = ""
    start = 0
    lang = ""
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
    """Fix a code block: use source code, keep translated comments."""
    if trans_block.lang != src_block.lang:
        raise StructureMismatchError(
            f"Code lang mismatch at line {trans_block.start_line}"
        )
    if len(trans_block.lines) != len(src_block.lines):
        raise StructureMismatchError(
            f"Code line count mismatch at line {trans_block.start_line}"
        )

    lang = trans_block.lang.lower()
    if lang == "mermaid":
        return src_block.lines.copy()

    hash_langs = {"python", "py", "sh", "bash", "dockerfile", "yaml", "yml", "toml"}
    slash_langs = {"console", "json"}
    mixed_langs = {"nextflow", "groovy", "nf", "java", "kotlin"}

    result = []
    in_block_comment = False

    for trans_line, src_line in zip(trans_block.lines, src_block.lines):
        if trans_line.lstrip().startswith("```"):
            result.append(src_line)
            continue

        if lang in mixed_langs:
            if BLOCK_COMMENT_START.match(trans_line):
                in_block_comment = True
                result.append(trans_line)
                continue
            if in_block_comment:
                if BLOCK_COMMENT_END.search(trans_line):
                    in_block_comment = False
                result.append(trans_line)
                continue

        trans_comment = None
        if lang in hash_langs or lang in mixed_langs:
            if m := HASH_COMMENT_RE.match(trans_line):
                trans_comment = m.group(2)
        elif lang in slash_langs:
            if m := SLASH_COMMENT_RE.match(trans_line):
                trans_comment = m.group(2)

        if trans_comment:
            src_comment = None
            if lang in hash_langs or lang in mixed_langs:
                if m := HASH_COMMENT_RE.match(src_line):
                    src_comment = m.group(2)
            elif lang in slash_langs:
                if m := SLASH_COMMENT_RE.match(src_line):
                    src_comment = m.group(2)
            result.append(
                src_line.replace(src_comment, trans_comment)
                if src_comment
                else src_line
            )
        else:
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
    """Fix a translation against its English source."""
    result = translation.copy()
    result = fix_headers(result, extract_headers(source))
    result = fix_links(result, extract_links(source))
    result = fix_html_links(result, extract_html_links(source))
    result = fix_code_blocks(result, extract_code_blocks(source))
    return result


def post_process_file(lang_path: Path, lang: str) -> bool:
    """Post-process a translation to fix common LLM mistakes."""
    en_path = lang_to_en_path(lang_path, lang)
    if not en_path.exists():
        return False

    original = lang_path.read_text(encoding="utf-8")
    en_content = en_path.read_text(encoding="utf-8")

    try:
        fixed = fix_translation(original.splitlines(), en_content.splitlines())
    except StructureMismatchError:
        return False

    fixed_text = "\n".join(fixed) + "\n"
    if fixed_text != original:
        lang_path.write_text(fixed_text, encoding="utf-8", newline="\n")
        return True
    return False


def post_process_language(lang: str, console: Console) -> int:
    """Post-process all files for a language. Returns count of fixed files."""
    lang_docs = DOCS_ROOT / lang / "docs"
    if not lang_docs.exists():
        return 0

    fixed = 0
    for path in lang_docs.rglob("*.md"):
        if post_process_file(path, lang):
            console.print(f"  Fixed: [cyan]{path.name}[/cyan]")
            fixed += 1
    return fixed


def get_languages_needing_sync(language: str | None = None) -> list[str]:
    """Get languages that have work to do."""
    all_langs = get_translation_languages()
    if language:
        if language not in all_langs:
            raise ConfigError(f"Unknown language: {language}. Available: {all_langs}")
        langs = [language]
    else:
        langs = all_langs

    return [
        lang
        for lang in langs
        if get_missing_files(lang)
        or get_outdated_files(lang)
        or get_orphaned_files(lang)
    ]


# =============================================================================
# CLI Commands
# =============================================================================

app = typer.Typer(help="Translation tools for Nextflow training docs")
console = Console()


@app.command()
def translate(
    path: Path = typer.Argument(..., help="English file to translate"),
    lang: str = typer.Option(..., "--lang", "-l", help="Target language code"),
):
    """Translate a single file."""
    check_api_key()

    if lang == "en":
        raise typer.Exit("Cannot translate to English (source language)")

    en_path = resolve_path(path)
    tf = TranslationFile(en_path, en_to_lang_path(en_path, lang), lang)

    translate_file(tf, console)
    post_process_file(tf.lang_path, lang)
    console.print(f"[green]Done:[/green] {tf.lang_path}")


@app.command()
def sync(
    lang: str = typer.Argument(..., help="Language code"),
    include: str | None = typer.Option(None, "--include", "-i", help="Filter pattern"),
    dry_run: bool = typer.Option(
        False, "--dry-run", "-n", help="Show what would be done"
    ),
):
    """Sync translations: update outdated, add missing, remove orphaned."""
    orphaned = get_orphaned_files(lang)
    outdated = get_outdated_files(lang)
    missing = get_missing_files(lang)

    if include:
        outdated = [tf for tf in outdated if include in str(tf.en_path)]
        missing = [tf for tf in missing if include in str(tf.en_path)]

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

    if orphaned:
        console.print(f"[bold]Removing {len(orphaned)} orphaned files...[/bold]")
        for p in orphaned:
            p.unlink()
            console.print(f"  [red]Removed:[/red] {p.relative_to(REPO_ROOT)}")

    if outdated:
        console.print(f"[bold]Updating {len(outdated)} outdated translations...[/bold]")
        for i, tf in enumerate(outdated, 1):
            console.print(f"[{i}/{len(outdated)}]", end="")
            translate_file(tf, console)
            post_process_file(tf.lang_path, lang)

    if missing:
        console.print(f"[bold]Adding {len(missing)} missing translations...[/bold]")
        for i, tf in enumerate(missing, 1):
            console.print(f"[{i}/{len(missing)}]", end="")
            translate_file(tf, console)
            post_process_file(tf.lang_path, lang)

    console.print("[green]Sync complete[/green]")


@app.command("post-process")
def post_process_cmd(
    lang: str = typer.Argument(..., help="Language code"),
    files: list[Path] | None = typer.Option(
        None, "--file", "-f", help="Specific files"
    ),
):
    """Run post-processing on translations to fix common issues."""
    if files:
        for path in files:
            if not path.exists():
                raise typer.Exit(f"File not found: {path}")
            if post_process_file(path, lang):
                console.print(f"[green]Fixed:[/green] {path}")
            else:
                console.print(f"[dim]OK:[/dim] {path}")
    else:
        fixed = post_process_language(lang, console)
        console.print(f"[green]Fixed {fixed} files[/green]")


# =============================================================================
# CI Commands
# =============================================================================


@app.command("ci-detect")
def ci_detect(language: str | None = typer.Option(None, "--language")):
    """Detect languages needing sync (outputs GitHub Actions format)."""
    langs = get_languages_needing_sync(language)
    print(f"languages={json.dumps(langs)}")
    print(f"has_work={'true' if langs else 'false'}")


@app.command("ci-run")
def ci_run(
    languages: str = typer.Option("[]", "--languages"),
    dry_run: bool = typer.Option(False, "--dry-run"),
):
    """Run translation sync for CI."""
    for lang in json.loads(languages):
        console.rule(f"[bold]{lang}[/bold]")
        sync(lang, include=None, dry_run=dry_run)


@app.command("ci-pr")
def ci_pr(
    github_token: str = typer.Option(..., envvar="GITHUB_TOKEN"),
    github_repository: str = typer.Option(..., envvar="GITHUB_REPOSITORY"),
):
    """Create PR with translation changes (for CI)."""
    from github import Github

    repo = git.Repo(REPO_ROOT)
    if not repo.is_dirty(untracked_files=True):
        console.print("[yellow]No changes to commit[/yellow]")
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
    gh_repo = gh.get_repo(github_repository)
    pr = gh_repo.create_pull(
        title="Update translations",
        body="""Update translations

This PR was created automatically using LLM translation.

To improve translations, edit the prompt files rather than individual translations:
- General: `_scripts/general-llm-prompt.md`
- Language-specific: `docs/*/llm-prompt.md`
""",
        base="master",
        head=branch,
    )
    console.print(f"[green]Created PR:[/green] {pr.html_url}")


if __name__ == "__main__":
    if not DOCS_ROOT.exists():
        sys.exit("Error: Must run from repository root")
    app()
