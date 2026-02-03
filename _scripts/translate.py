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
    uv run translate.py list-missing pt
    uv run translate.py translate docs/en/docs/index.md --lang pt
    uv run translate.py sync pt
"""

from __future__ import annotations

import json
import os
import secrets
import subprocess
import sys
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import TYPE_CHECKING

import anthropic
import git
import typer
import yaml
from rich.console import Console

if TYPE_CHECKING:
    from github import Github

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
            f"## Task",
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


def get_changed_files(before_sha: str, after_sha: str) -> list[str]:
    """Get files changed between two commits."""
    repo = git.Repo(REPO_ROOT)
    try:
        diff = repo.git.diff("--name-only", before_sha, after_sha)
        return diff.split("\n") if diff else []
    except git.GitCommandError:
        return []


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
# Post-processing
# =============================================================================


def post_process_file(lang_path: Path, lang: str) -> bool:
    """
    Post-process a translation to fix common LLM mistakes.

    Fixes: links, header anchors, code blocks.
    Returns True if changes were made.
    """
    from translation_fixer import StructureMismatchError, fix_translation

    en_path = lang_to_en_path(lang_path, lang)
    if not en_path.exists():
        return False

    original = lang_path.read_text(encoding="utf-8")
    en_content = en_path.read_text(encoding="utf-8")

    try:
        fixed = fix_translation(
            translation=original.splitlines(),
            source=en_content.splitlines(),
        )
    except StructureMismatchError:
        # Structure mismatch is not fatal - just skip post-processing
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
    """Get languages that have work to do (missing, outdated, or orphaned files)."""
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
    tf = TranslationFile(
        en_path=en_path,
        lang_path=en_to_lang_path(en_path, lang),
        language=lang,
    )

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
    orphaned_files = get_orphaned_files(lang)
    outdated_files = get_outdated_files(lang)
    missing_files = get_missing_files(lang)

    if include:
        outdated_files = [tf for tf in outdated_files if include in str(tf.en_path)]
        missing_files = [tf for tf in missing_files if include in str(tf.en_path)]

    if dry_run:
        if orphaned_files:
            console.print(f"[red]Would remove {len(orphaned_files)} orphaned:[/red]")
            for p in orphaned_files:
                console.print(f"  {p.relative_to(REPO_ROOT)}")
        if outdated_files:
            console.print(
                f"[yellow]Would update {len(outdated_files)} outdated:[/yellow]"
            )
            for tf in outdated_files:
                console.print(f"  {tf.relative_path}")
        if missing_files:
            console.print(f"[blue]Would add {len(missing_files)} missing:[/blue]")
            for tf in missing_files:
                console.print(f"  {tf.relative_path}")
        if not (orphaned_files or outdated_files or missing_files):
            console.print(f"[green]Nothing to do for {lang}[/green]")
        return

    check_api_key()

    # Remove orphaned
    if orphaned_files:
        console.print(f"[bold]Removing {len(orphaned_files)} orphaned files...[/bold]")
        for p in orphaned_files:
            p.unlink()
            console.print(f"  [red]Removed:[/red] {p.relative_to(REPO_ROOT)}")

    # Update outdated
    if outdated_files:
        console.print(
            f"[bold]Updating {len(outdated_files)} outdated translations...[/bold]"
        )
        for i, tf in enumerate(outdated_files, 1):
            console.print(f"[{i}/{len(outdated_files)}]", end="")
            translate_file(tf, console)
            post_process_file(tf.lang_path, lang)

    # Add missing
    if missing_files:
        console.print(
            f"[bold]Adding {len(missing_files)} missing translations...[/bold]"
        )
        for i, tf in enumerate(missing_files, 1):
            console.print(f"[{i}/{len(missing_files)}]", end="")
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
def ci_detect(
    language: str | None = typer.Option(None, "--language"),
):
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
    langs: list[str] = json.loads(languages)

    for lang in langs:
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

    # Configure git
    subprocess.run(
        ["git", "config", "user.name", "github-actions[bot]"],
        check=True,
        cwd=REPO_ROOT,
    )
    subprocess.run(
        ["git", "config", "user.email", "github-actions[bot]@users.noreply.github.com"],
        check=True,
        cwd=REPO_ROOT,
    )

    # Create branch and commit
    branch = f"translate-{secrets.token_hex(4)}"
    subprocess.run(["git", "checkout", "-b", branch], check=True, cwd=REPO_ROOT)
    subprocess.run(["git", "add", "docs/"], check=True, cwd=REPO_ROOT)

    msg = "Update translations"
    subprocess.run(["git", "commit", "-m", msg], check=True, cwd=REPO_ROOT)
    subprocess.run(["git", "push", "origin", branch], check=True, cwd=REPO_ROOT)

    # Create PR
    g = Github(github_token)
    gh_repo = g.get_repo(github_repository)
    body = """Update translations

This PR was created automatically using LLM translation.

To improve translations, edit the prompt files rather than individual translations:
- General: `_scripts/general-llm-prompt.md`
- Language-specific: `docs/*/llm-prompt.md`
"""

    pr = gh_repo.create_pull(title=msg, body=body, base="master", head=branch)
    console.print(f"[green]Created PR:[/green] {pr.html_url}")


if __name__ == "__main__":
    if not DOCS_ROOT.exists():
        sys.exit("Error: Must run from repository root")
    app()
