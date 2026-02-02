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
LLM-powered translation using Anthropic Claude.

Usage:
    uv run translate.py list-missing pt
    uv run translate.py translate-page --language pt --en-path docs/en/docs/index.md
    uv run translate.py update-outdated --language pt
"""

import json
import secrets
import subprocess
import sys
from functools import lru_cache
from pathlib import Path

import anthropic
import git
import typer
import yaml
from github import Github
from rich.console import Console
from rich.progress import Progress, SpinnerColumn, TextColumn

app = typer.Typer(help="Translation tools for Nextflow training docs")
console = Console()

# Must be run from repo root
REPO_ROOT = Path(__file__).parent.parent
DOCS_PATH = REPO_ROOT / "docs"
EN_DOCS_PATH = DOCS_PATH / "en" / "docs"
SCRIPTS_PATH = REPO_ROOT / "_scripts"

# Model to use for translations
# Using alias to automatically get the latest Sonnet version
# Sonnet is ~40% cheaper than Opus and recommended for translation tasks
MODEL = "claude-sonnet-4-5"

# Priority directories for translation (most important first)
PRIORITY_DIRS = ["hello_nextflow", "hello_nf-core", "nf4_science", "envsetup"]


@lru_cache
def get_general_prompt() -> str:
    """Load the general translation prompt."""
    prompt_path = SCRIPTS_PATH / "general-llm-prompt.md"
    return prompt_path.read_text(encoding="utf-8")


@lru_cache
def get_langs() -> dict[str, str]:
    """Load language names from YAML file."""
    return yaml.safe_load(
        (DOCS_PATH / "language_names.yml").read_text(encoding="utf-8")
    )


def get_lang_prompt(lang: str) -> str:
    """Load language-specific prompt."""
    prompt_path = DOCS_PATH / lang / "llm-prompt.md"
    if not prompt_path.exists():
        console.print(f"[red]Error:[/red] Prompt file not found: {prompt_path}")
        raise typer.Exit(code=1)
    return prompt_path.read_text(encoding="utf-8")


def to_lang_path(en_path: Path, lang: str) -> Path:
    """Convert an English docs path to a language-specific path."""
    rel_path = en_path.relative_to(EN_DOCS_PATH)
    return DOCS_PATH / lang / "docs" / rel_path


def to_en_path(lang_path: Path, lang: str) -> Path:
    """Convert a language-specific path to an English docs path."""
    rel_path = lang_path.relative_to(DOCS_PATH / lang / "docs")
    return EN_DOCS_PATH / rel_path


def iter_en_docs() -> list[Path]:
    """
    Iterate English markdown files in priority order.

    Returns files in this order:
    1. Root-level .md files
    2. Priority directories (hello_nextflow, etc.)
    3. Remaining directories
    """
    paths: list[Path] = []

    # Root level files first
    paths.extend(sorted(EN_DOCS_PATH.glob("*.md")))

    # Priority directories
    for dir_name in PRIORITY_DIRS:
        dir_path = EN_DOCS_PATH / dir_name
        if dir_path.exists():
            paths.extend(sorted(dir_path.rglob("*.md")))

    # Remaining directories (skip already processed)
    priority_prefixes = tuple(str(EN_DOCS_PATH / d) for d in PRIORITY_DIRS)
    for path in sorted(EN_DOCS_PATH.rglob("*.md")):
        if path.parent == EN_DOCS_PATH:
            continue  # Already added root files
        if str(path).startswith(priority_prefixes):
            continue  # Already added priority dirs
        paths.append(path)

    return paths


def check_api_key() -> None:
    """Check that ANTHROPIC_API_KEY is set and provide helpful error if not."""
    import os

    if not os.environ.get("ANTHROPIC_API_KEY"):
        console.print(
            "[red]Error:[/red] ANTHROPIC_API_KEY environment variable not set"
        )
        console.print()
        console.print("[dim]To fix this, set the environment variable:[/dim]")
        console.print("  export ANTHROPIC_API_KEY='your-api-key-here'")
        console.print()
        console.print("[dim]Get your API key from:[/dim]")
        console.print("  https://console.anthropic.com/settings/keys")
        raise typer.Exit(code=1)


def resolve_en_path(en_path: Path) -> Path:
    """
    Resolve an English docs path to an absolute path.

    Accepts paths in several formats:
    - Absolute path: /full/path/to/docs/en/docs/file.md
    - Relative to repo root: docs/en/docs/file.md
    - Relative to EN_DOCS_PATH: hello_nextflow/index.md
    """
    # If already absolute and exists, use it
    if en_path.is_absolute() and en_path.exists():
        return en_path

    # Try as relative to current directory
    if en_path.exists():
        return en_path.resolve()

    # Try as relative to repo root
    repo_relative = REPO_ROOT / en_path
    if repo_relative.exists():
        return repo_relative

    # Try as relative to EN_DOCS_PATH (e.g., "hello_nextflow/index.md")
    en_docs_relative = EN_DOCS_PATH / en_path
    if en_docs_relative.exists():
        return en_docs_relative

    # Return original for error message
    return en_path


@app.command()
def translate_page(
    language: str = typer.Option(..., "--language", "-l", envvar="LANGUAGE"),
    en_path: Path = typer.Option(..., "--en-path", "-p", envvar="EN_PATH"),
):
    """Translate a single page from English to another language."""
    check_api_key()

    if language == "en":
        console.print("[red]Error:[/red] Cannot translate to English (source language)")
        raise typer.Exit(code=1)

    # Resolve the path
    en_path = resolve_en_path(en_path)

    if not en_path.exists():
        console.print(f"[red]Error:[/red] Source file not found: {en_path}")
        raise typer.Exit(code=1)

    langs = get_langs()
    if language not in langs:
        console.print(f"[red]Error:[/red] Unknown language: {language}")
        raise typer.Exit(code=1)

    language_name = langs[language]
    out_path = to_lang_path(en_path, language)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    original_content = en_path.read_text(encoding="utf-8")
    old_translation = (
        out_path.read_text(encoding="utf-8") if out_path.exists() else None
    )

    console.print(f"Translating [cyan]{en_path.name}[/cyan] to {language_name}")

    # Build prompt
    prompt_parts = [get_general_prompt(), get_lang_prompt(language)]

    if old_translation:
        console.print("  [dim]Updating existing translation[/dim]")
        prompt_parts.extend(
            [
                "There is an existing translation that may be outdated. Update only where necessary:",
                "- Add new parts from the English source",
                "- Remove parts deleted from the English source",
                "- Update changed parts",
                "- Fix any violations of current instructions",
                "- Otherwise preserve the translation LINE-BY-LINE, AS-IS",
                "",
                "Do NOT rephrase correct lines or change formatting unless required.",
                "Minimize diffs for easier human review.",
                "",
                "Previous translation:",
                f"%%%\n{old_translation}%%%",
            ]
        )

    prompt_parts.extend(
        [
            f"Translate to {language} ({language_name}).",
            "Original content:",
            f"%%%\n{original_content}%%%",
        ]
    )

    prompt = "\n\n".join(prompt_parts)

    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        console=console,
    ) as progress:
        progress.add_task("Calling Claude API...", total=None)
        client = anthropic.Anthropic()
        message = client.messages.create(
            model=MODEL,
            max_tokens=16384,
            messages=[{"role": "user", "content": prompt}],
        )

    out_content = f"{message.content[0].text.strip()}\n"
    out_path.write_text(out_content, encoding="utf-8", newline="\n")
    console.print(f"  [green]Saved:[/green] {out_path}")


@app.command()
def list_missing(language: str):
    """List English files that haven't been translated."""
    missing = []
    for en_path in iter_en_docs():
        lang_path = to_lang_path(en_path, language)
        if not lang_path.exists():
            missing.append(en_path)

    if not missing:
        console.print(f"[green]All files translated for {language}[/green]")
        return

    console.print(f"[yellow]Missing translations for {language}:[/yellow]")
    for path in missing:
        console.print(f"  {path.relative_to(REPO_ROOT)}")
    console.print(f"\nTotal: {len(missing)} files")


@app.command()
def list_outdated(language: str):
    """List translations older than their English source (by git commit time)."""
    repo = git.Repo(REPO_ROOT)
    outdated = []

    for en_path in iter_en_docs():
        lang_path = to_lang_path(en_path, language)
        if not lang_path.exists():
            continue

        en_commits = list(repo.iter_commits(paths=str(en_path), max_count=1))
        lang_commits = list(repo.iter_commits(paths=str(lang_path), max_count=1))

        if not en_commits or not lang_commits:
            continue

        if lang_commits[0].committed_datetime < en_commits[0].committed_datetime:
            outdated.append(en_path)

    if not outdated:
        console.print(f"[green]All translations up to date for {language}[/green]")
        return

    console.print(f"[yellow]Outdated translations for {language}:[/yellow]")
    for path in outdated:
        console.print(f"  {path.relative_to(REPO_ROOT)}")
    console.print(f"\nTotal: {len(outdated)} files")


@app.command()
def list_removable(language: str):
    """List translation files without English source (orphaned)."""
    lang_docs = DOCS_PATH / language / "docs"
    if not lang_docs.exists():
        console.print(f"[red]Error:[/red] Language docs not found: {lang_docs}")
        raise typer.Exit(code=1)

    removable = []
    for lang_path in lang_docs.rglob("*.md"):
        en_path = to_en_path(lang_path, language)
        if not en_path.exists():
            removable.append(lang_path)

    if not removable:
        console.print(f"[green]No orphaned files for {language}[/green]")
        return

    console.print(f"[yellow]Orphaned files for {language}:[/yellow]")
    for path in removable:
        console.print(f"  {path.relative_to(REPO_ROOT)}")
    console.print(f"\nTotal: {len(removable)} files")


@app.command()
def update_outdated(
    language: str = typer.Option(..., "--language", "-l", envvar="LANGUAGE"),
    verify: bool = typer.Option(
        False,
        "--verify",
        "-v",
        help="Run verification pass after translation to ensure compliance with prompt guidelines",
    ),
    verify_iterations: int = typer.Option(
        2,
        "--verify-iterations",
        help="Number of verification iterations (only used with --verify)",
    ),
):
    """Update all outdated translations for a language."""
    repo = git.Repo(REPO_ROOT)
    outdated = []

    for en_path in iter_en_docs():
        lang_path = to_lang_path(en_path, language)
        if not lang_path.exists():
            continue

        en_commits = list(repo.iter_commits(paths=str(en_path), max_count=1))
        lang_commits = list(repo.iter_commits(paths=str(lang_path), max_count=1))

        if en_commits and lang_commits:
            if lang_commits[0].committed_datetime < en_commits[0].committed_datetime:
                outdated.append(en_path)

    if not outdated:
        console.print(f"[green]All translations up to date for {language}[/green]")
        return

    console.print(f"Updating {len(outdated)} outdated translations...")
    translated_paths = []
    for i, en_path in enumerate(outdated, 1):
        console.print(f"\n[{i}/{len(outdated)}] {en_path.name}")
        translate_page(language=language, en_path=en_path)
        translated_paths.append(to_lang_path(en_path, language))

    # Run verification pass if requested
    if verify and translated_paths:
        console.print(
            f"\n[bold]Running verification pass on {len(translated_paths)} translated files...[/bold]"
        )
        for iteration in range(1, verify_iterations + 1):
            console.print(
                f"\n[cyan]Verification pass {iteration}/{verify_iterations}[/cyan]"
            )
            changes = 0
            for lang_path in translated_paths:
                if fix_single_file(lang_path, language):
                    console.print(f"  [yellow]Fixed:[/yellow] {lang_path.name}")
                    changes += 1
            if changes == 0:
                console.print("[green]All translations verified![/green]")
                break
            console.print(f"  [dim]Fixed {changes} files in this pass[/dim]")


@app.command()
def add_missing(
    language: str = typer.Option(..., "--language", "-l", envvar="LANGUAGE"),
    include: str = typer.Option(
        None,
        "--include",
        "-i",
        help="Only translate files matching this pattern (e.g., 'hello_nextflow')",
    ),
    verify: bool = typer.Option(
        False,
        "--verify",
        "-v",
        help="Run verification pass after translation to ensure compliance with prompt guidelines",
    ),
    verify_iterations: int = typer.Option(
        2,
        "--verify-iterations",
        help="Number of verification iterations (only used with --verify)",
    ),
):
    """Translate all missing files for a language."""
    missing = []
    for en_path in iter_en_docs():
        # Filter by include pattern if specified
        if include and include not in str(en_path):
            continue
        lang_path = to_lang_path(en_path, language)
        if not lang_path.exists():
            missing.append(en_path)

    if not missing:
        if include:
            console.print(
                f"[green]All matching files already translated for {language}[/green]"
            )
        else:
            console.print(f"[green]All files already translated for {language}[/green]")
        return

    console.print(f"Translating {len(missing)} missing files...")
    translated_paths = []
    for i, en_path in enumerate(missing, 1):
        console.print(f"\n[{i}/{len(missing)}] {en_path.name}")
        translate_page(language=language, en_path=en_path)
        translated_paths.append(to_lang_path(en_path, language))

    # Run verification pass if requested
    if verify and translated_paths:
        console.print(
            f"\n[bold]Running verification pass on {len(translated_paths)} translated files...[/bold]"
        )
        for iteration in range(1, verify_iterations + 1):
            console.print(
                f"\n[cyan]Verification pass {iteration}/{verify_iterations}[/cyan]"
            )
            changes = 0
            for lang_path in translated_paths:
                if fix_single_file(lang_path, language):
                    console.print(f"  [yellow]Fixed:[/yellow] {lang_path.name}")
                    changes += 1
            if changes == 0:
                console.print("[green]All translations verified![/green]")
                break
            console.print(f"  [dim]Fixed {changes} files in this pass[/dim]")


@app.command()
def remove_removable(
    language: str = typer.Option(..., "--language", "-l", envvar="LANGUAGE"),
):
    """Remove orphaned translation files (no English source)."""
    lang_docs = DOCS_PATH / language / "docs"
    removed = 0

    for lang_path in lang_docs.rglob("*.md"):
        en_path = to_en_path(lang_path, language)
        if not en_path.exists():
            lang_path.unlink()
            console.print(f"[red]Removed:[/red] {lang_path.relative_to(REPO_ROOT)}")
            removed += 1

    console.print(f"\nRemoved {removed} orphaned files")


def fix_single_file(
    lang_path: Path,
    language: str,
    run_post_fixer: bool = True,
) -> bool:
    """
    Fix a single translation file to comply with current prompt guidelines.

    Returns True if changes were made, False if file was already compliant.
    """
    from translation_fixer import fix_one_page

    if not lang_path.exists():
        console.print(f"[red]Error:[/red] File not found: {lang_path}")
        return False

    # Find English source
    en_path = to_en_path(lang_path, language)
    if not en_path.exists():
        console.print(
            f"[yellow]Warning:[/yellow] No English source for {lang_path.name}, skipping"
        )
        return False

    original_content = lang_path.read_text(encoding="utf-8")
    en_content = en_path.read_text(encoding="utf-8")

    # Build fix prompt
    prompt_parts = [
        "# Translation Fix Task",
        "",
        "You are reviewing an existing translation of Nextflow training materials.",
        "",
        "## Your Task",
        "Check the translation against the guidelines below and fix any violations.",
        "Make MINIMAL changes - only fix what violates the guidelines.",
        "If the translation is already correct, output it unchanged.",
        "",
        "## Guidelines",
        get_general_prompt(),
        get_lang_prompt(language),
        "",
        "## English Source (for reference)",
        f"%%%\n{en_content}%%%",
        "",
        "## Current Translation to Fix",
        f"%%%\n{original_content}%%%",
        "",
        "## Output",
        "Output ONLY the corrected markdown file. No explanations or commentary.",
    ]

    prompt = "\n\n".join(prompt_parts)

    client = anthropic.Anthropic()
    message = client.messages.create(
        model=MODEL,
        max_tokens=16384,
        messages=[{"role": "user", "content": prompt}],
    )

    fixed_content = f"{message.content[0].text.strip()}\n"

    # Check if content changed
    if fixed_content == original_content:
        return False

    # Write fixed content
    lang_path.write_text(fixed_content, encoding="utf-8", newline="\n")

    # Run post-processing fixer (restores code blocks, links, etc.)
    if run_post_fixer:
        try:
            fix_one_page(lang_path)
        except ValueError as e:
            console.print(f"[yellow]Warning:[/yellow] Post-fixer issue: {e}")

    return True


@app.command()
def fix_translations(
    language: str = typer.Option(..., "--language", "-l", help="Language code to fix"),
    max_iterations: int = typer.Option(
        8, "--max-iterations", "-m", help="Maximum fix iterations"
    ),
    files: list[Path] = typer.Option(
        None,
        "--files",
        "-f",
        help="Specific files to fix (relative to docs/{lang}/docs/)",
    ),
    dry_run: bool = typer.Option(
        False, "--dry-run", help="Show what would be fixed without making changes"
    ),
):
    """
    Fix existing translations to comply with current prompt guidelines.

    Uses LLM to review and fix translations, with multiple passes until
    no more changes are needed or max iterations reached.
    """
    langs = get_langs()
    if language not in langs:
        console.print(f"[red]Error:[/red] Unknown language: {language}")
        raise typer.Exit(code=1)

    lang_docs = DOCS_PATH / language / "docs"
    if not lang_docs.exists():
        console.print(f"[red]Error:[/red] Language docs not found: {lang_docs}")
        raise typer.Exit(code=1)

    # Get list of files to process
    if files:
        # Resolve specified files
        paths_to_fix = []
        for f in files:
            # Try as relative to lang docs
            full_path = lang_docs / f
            if full_path.exists():
                paths_to_fix.append(full_path)
            elif f.exists():
                paths_to_fix.append(f.resolve())
            else:
                console.print(f"[red]Error:[/red] File not found: {f}")
                raise typer.Exit(code=1)
    else:
        # Get all translation files
        paths_to_fix = sorted(lang_docs.rglob("*.md"))

    if not paths_to_fix:
        console.print(f"[yellow]No translation files found for {language}[/yellow]")
        return

    console.print(
        f"[bold]Fixing {len(paths_to_fix)} translation files for {langs[language]}[/bold]"
    )

    if dry_run:
        console.print("[yellow]DRY RUN - no changes will be made[/yellow]")
        for p in paths_to_fix:
            console.print(f"  Would check: {p.relative_to(REPO_ROOT)}")
        return

    # Check API key only when actually running (not dry-run)
    check_api_key()

    total_changes = 0
    for iteration in range(1, max_iterations + 1):
        console.print(f"\n[bold cyan]Pass {iteration}/{max_iterations}[/bold cyan]")

        changes_this_pass = 0
        for i, lang_path in enumerate(paths_to_fix, 1):
            console.print(f"  [{i}/{len(paths_to_fix)}] {lang_path.name}", end=" ")

            with Progress(
                SpinnerColumn(),
                TextColumn("[progress.description]{task.description}"),
                console=console,
                transient=True,
            ) as progress:
                progress.add_task("checking...", total=None)
                changed = fix_single_file(lang_path, language)

            if changed:
                console.print("[yellow]fixed[/yellow]")
                changes_this_pass += 1
            else:
                console.print("[green]OK[/green]")

        total_changes += changes_this_pass
        console.print(
            f"  [dim]Pass {iteration} complete: {changes_this_pass} files changed[/dim]"
        )

        if changes_this_pass == 0:
            console.print(
                f"\n[green]All translations now comply with guidelines![/green]"
            )
            console.print(f"Total changes across all passes: {total_changes}")
            return

    console.print(f"\n[yellow]Max iterations ({max_iterations}) reached.[/yellow]")
    console.print(f"Total changes: {total_changes}")
    console.print("Some issues may remain - consider running again or manual review.")


@app.command()
def make_pr(
    language: str = typer.Option(None, "--language", "-l", envvar="LANGUAGE"),
    command: str = typer.Option(None, "--command", "-c", envvar="COMMAND"),
    github_token: str = typer.Option(..., envvar="GITHUB_TOKEN"),
    github_repository: str = typer.Option(..., envvar="GITHUB_REPOSITORY"),
):
    """Create a PR with translation changes (for CI use)."""
    repo = git.Repo(REPO_ROOT)

    if not repo.is_dirty(untracked_files=True):
        console.print("[yellow]No changes to commit[/yellow]")
        return

    # Configure git
    subprocess.run(
        ["git", "config", "user.name", "github-actions[bot]"], check=True, cwd=REPO_ROOT
    )
    subprocess.run(
        ["git", "config", "user.email", "github-actions[bot]@users.noreply.github.com"],
        check=True,
        cwd=REPO_ROOT,
    )

    # Create branch
    branch_name = (
        f"translate-{language or 'all'}-{command or 'update'}-{secrets.token_hex(4)}"
    )
    subprocess.run(["git", "checkout", "-b", branch_name], check=True, cwd=REPO_ROOT)

    # Commit
    subprocess.run(["git", "add", "docs/"], check=True, cwd=REPO_ROOT)
    message = "Update translations"
    if language:
        message += f" for {language}"
    if command:
        message += f" ({command})"
    subprocess.run(["git", "commit", "-m", message], check=True, cwd=REPO_ROOT)

    # Push
    subprocess.run(["git", "push", "origin", branch_name], check=True, cwd=REPO_ROOT)

    # Create PR
    g = Github(github_token)
    gh_repo = g.get_repo(github_repository)
    body = f"""{message}

This PR was created automatically using LLM translation.

Prompt file: https://github.com/{github_repository}/blob/master/docs/{language}/llm-prompt.md

To improve translations, edit the prompt file rather than individual translations."""

    pr = gh_repo.create_pull(title=message, body=body, base="master", head=branch_name)
    console.print(f"[green]Created PR:[/green] {pr.html_url}")


if __name__ == "__main__":
    if not DOCS_PATH.exists():
        console.print("[red]Error:[/red] Must run from repository root")
        sys.exit(1)
    app()
