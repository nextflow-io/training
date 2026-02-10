"""CLI commands for the translation tool.

Provides the typer application with all user-facing commands:
translate, sync, ci-detect, ci-delete-orphans, ci-run, ci-pr.
"""

from __future__ import annotations

import asyncio
import json
import secrets
import subprocess
from pathlib import Path

import typer
from rich.console import Console

from .config import (
    DOCS_ROOT,
    EN_DOCS,
    REPO_ROOT,
    ConfigError,
    DEFAULT_PARALLEL,
    check_api_key,
    make_console,
)
from .core import translate_all
from .git_utils import (
    delete_orphans,
    gather_work,
    get_missing_files,
    get_orphaned_files,
    prompt_changed_since,
    resolve_baseline,
)
from .verify import get_broken_files
from .models import TranslationFile, TranslationLog
from .paths import en_to_lang_path, get_translation_languages


app = typer.Typer(help="Translation CLI for Nextflow training docs")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _resolve_en_path(path: Path) -> Path:
    """Resolve *path* to an English doc file."""
    if path.is_absolute():
        if path.exists():
            return path
        raise ConfigError(f"File not found: {path}")

    if (EN_DOCS / path).exists():
        return EN_DOCS / path

    if path.exists():
        return path.resolve()

    raise ConfigError(f"File not found: {path}")


# ---------------------------------------------------------------------------
# sync sub-modes
# ---------------------------------------------------------------------------


def _sync_fix(
    lang: str,
    include: str | None,
    parallel: int,
    log_file: Path | None,
    console: Console,
) -> None:
    """Fix mode: scan all translations for issues and re-translate broken ones."""
    console.print(
        f"[bold magenta]Fix mode:[/bold magenta] scanning all {lang} translations for issues..."
    )

    check_api_key()

    orphaned = get_orphaned_files(lang)
    missing = get_missing_files(lang)
    broken = asyncio.run(get_broken_files(lang, parallel))

    if include:
        missing = [tf for tf in missing if include in str(tf.en_path)]
        broken = [tf for tf in broken if include in str(tf.en_path)]

    delete_orphans(orphaned, console)

    # Deduplicate broken + missing (missing files may appear in both)
    all_files = list({str(tf.en_path): tf for tf in broken + missing}.values())

    if all_files:
        console.print(
            f"[bold]Re-translating {len(all_files)} files "
            f"({len(broken)} broken + {len(missing)} missing, deduplicated) "
            f"with parallel={parallel}[/bold]"
        )
        log = (
            TranslationLog(log_file, "fix-mode", lang, parallel, len(all_files))
            if log_file
            else None
        )
        asyncio.run(translate_all(all_files, parallel, log))
    else:
        console.print(f"[green]No broken or missing files found for {lang}[/green]")

    console.print("[green]✓ Fix complete[/green]")


def _sync_normal(
    lang: str,
    include: str | None,
    since: str | None,
    parallel: int,
    log_file: Path | None,
    console: Console,
) -> None:
    """Normal mode: baseline-based sync."""
    check_api_key()
    baseline = resolve_baseline(since, console)

    orphaned, missing, outdated = gather_work(lang, baseline)

    prompts_changed = baseline and prompt_changed_since(lang, baseline)
    if prompts_changed:
        console.print(
            f"[magenta]Prompt changed:[/magenta] all {lang} translations will be updated"
        )

    if include:
        outdated = [tf for tf in outdated if include in str(tf.en_path)]
        missing = [tf for tf in missing if include in str(tf.en_path)]

    delete_orphans(orphaned, console)

    all_files = outdated + missing
    if all_files:
        console.print(
            f"Translating {len(all_files)} files with parallel={parallel}. "
            f"Largest files first."
        )
        log = (
            TranslationLog(log_file, baseline, lang, parallel, len(all_files))
            if log_file
            else None
        )
        asyncio.run(translate_all(all_files, parallel, log))

    console.print("[green]✓ Sync complete[/green]")


# ---------------------------------------------------------------------------
# Commands
# ---------------------------------------------------------------------------


@app.command()
def translate(
    path: Path = typer.Argument(..., help="English file to translate"),
    lang: str = typer.Option(..., "--lang", "-l", help="Target language code"),
) -> None:
    """Translate a single file."""
    check_api_key()

    valid_langs = get_translation_languages()
    if lang not in valid_langs:
        console = make_console(stderr=True)
        console.print(
            f"[red]Unknown language: {lang}. Valid: {', '.join(sorted(valid_langs))}[/red]"
        )
        raise typer.Exit(code=1)

    en_path = _resolve_en_path(path)
    tf = TranslationFile(en_path, en_to_lang_path(en_path, lang), lang)

    # Use the async path for consistency (includes verification)
    asyncio.run(_translate_single(tf))

    console = make_console()
    console.print(f"[green]Done:[/green] {tf.lang_path}")


async def _translate_single(tf: TranslationFile) -> None:
    """Translate a single file through the unified async path."""
    import anthropic as _anthropic

    from .core import translate_file_async
    from .progress import TranslationProgress

    semaphore = asyncio.Semaphore(1)
    file_lines = {str(tf.relative_path): tf.en_path.read_text().count("\n")}
    progress = TranslationProgress([tf], file_lines)
    client = _anthropic.AsyncAnthropic()

    await translate_file_async(tf, semaphore, progress, client)


@app.command()
def sync(
    lang: str = typer.Argument(..., help="Language code"),
    include: str | None = typer.Option(None, "--include", "-i", help="Filter pattern"),
    since: str | None = typer.Option(
        None,
        "--since",
        help="Compare since this commit (default: auto-detect from last translation PR)",
    ),
    fix: bool = typer.Option(
        False,
        "--fix",
        help="Fix mode: scan all translations for issues and re-translate broken files",
    ),
    parallel: int = typer.Option(
        DEFAULT_PARALLEL,
        "--parallel",
        "-p",
        help="Max concurrent translations (default: 50)",
    ),
    log_file: Path | None = typer.Option(
        None,
        "--log",
        help="Write detailed JSON log to file",
    ),
) -> None:
    """Sync translations: update outdated, add missing, remove orphaned."""
    console = make_console()

    if fix:
        _sync_fix(lang, include, parallel, log_file, console)
    else:
        _sync_normal(lang, include, since, parallel, log_file, console)


# ---------------------------------------------------------------------------
# CI Commands
# ---------------------------------------------------------------------------


@app.command("ci-detect")
def ci_detect(
    language: str | None = typer.Option(None, "--language"),
    since: str | None = typer.Option(
        None,
        "--since",
        help="Compare since this commit (default: auto-detect from last translation PR)",
    ),
    fix: bool = typer.Option(
        False,
        "--fix",
        help="Fix mode: detect languages with broken translations instead of outdated ones",
    ),
) -> None:
    """Detect languages needing sync (GitHub Actions output)."""
    console = make_console(stderr=True)
    all_langs = get_translation_languages()

    if language:
        if language not in all_langs:
            console.print(
                f"[red]Unknown language: {language}. Valid: {', '.join(sorted(all_langs))}[/red]"
            )
            raise typer.Exit(code=1)
        langs = [language]
    else:
        langs = all_langs

    if fix:
        need_sync = _ci_detect_fix(langs, console)
    else:
        need_sync = _ci_detect_normal(langs, since, console)

    print(f"languages={json.dumps(need_sync)}")
    print(f"has_work={'true' if need_sync else 'false'}")


def _ci_detect_fix(langs: list[str], console: Console) -> list[str]:
    """Fix mode: scan all translations for structural/semantic issues."""
    console.print(
        "[bold magenta]Fix mode:[/bold magenta] scanning for broken translations..."
    )
    check_api_key()

    async def _detect_broken():
        need_sync = []
        for lang in langs:
            orphaned = get_orphaned_files(lang)
            missing = get_missing_files(lang)
            broken = await get_broken_files(lang)

            if orphaned or missing or broken:
                need_sync.append(lang)
                console.print(f"[bold cyan]{lang}[/bold cyan]:")
                if broken:
                    console.print(f"  [yellow]Broken:[/yellow] {len(broken)}")
                    for tf in broken:
                        console.print(f"    {tf.relative_path}")
                if missing:
                    console.print(f"  [green]Missing:[/green] {len(missing)}")
                    for f in missing:
                        console.print(f"    {f.relative_path}")
                if orphaned:
                    console.print(f"  [red]Orphaned:[/red] {len(orphaned)}")
                    for f in orphaned:
                        console.print(f"    {f.relative_to(DOCS_ROOT)}")
        return need_sync

    return asyncio.run(_detect_broken())


def _ci_detect_normal(
    langs: list[str], since: str | None, console: Console
) -> list[str]:
    """Normal mode: baseline-based detection."""
    baseline = resolve_baseline(since, console)

    need_sync = []
    for lang in langs:
        orphaned, missing, outdated = gather_work(lang, baseline)
        prompts_changed = baseline and prompt_changed_since(lang, baseline)

        if missing or outdated or orphaned:
            need_sync.append(lang)
            console.print(f"[bold cyan]{lang}[/bold cyan]:")
            if prompts_changed:
                console.print(
                    f"  [magenta]Prompt changed:[/magenta] all translations will be updated"
                )
            if outdated:
                console.print(f"  [yellow]Outdated:[/yellow] {len(outdated)}")
                if not prompts_changed:
                    for f in outdated:
                        console.print(f"    {f.relative_path}")
            if missing:
                console.print(f"  [green]Missing:[/green] {len(missing)}")
                for f in missing:
                    console.print(f"    {f.relative_path}")
            if orphaned:
                console.print(f"  [red]Orphaned:[/red] {len(orphaned)}")
                for f in orphaned:
                    console.print(f"    {f.relative_to(DOCS_ROOT)}")

    return need_sync


@app.command("ci-delete-orphans")
def ci_delete_orphans() -> None:
    """Delete orphaned translation files for all languages (CI).

    Runs in the merge job to delete translation files that no longer
    have a corresponding English source file. Files are deleted and staged
    for commit using git rm.
    """
    console = make_console(stderr=True)
    all_langs = get_translation_languages()
    total_deleted = 0

    for lang in all_langs:
        orphaned = get_orphaned_files(lang)
        if orphaned:
            console.print(
                f"[bold cyan]{lang}[/bold cyan]: deleting {len(orphaned)} orphaned files"
            )
            delete_orphans(orphaned, console, use_git_rm=True)
            total_deleted += len(orphaned)

    if total_deleted:
        console.print(f"\n[bold]Total deleted:[/bold] {total_deleted} orphaned files")
    else:
        console.print("[dim]No orphaned files found[/dim]")


@app.command("ci-run")
def ci_run(
    languages: str = typer.Option("[]", "--languages"),
    fix: bool = typer.Option(False, "--fix", help="Fix mode"),
) -> None:
    """Run sync for multiple languages (CI)."""
    console = make_console()
    lang_list = json.loads(languages)
    if not isinstance(lang_list, list):
        raise ConfigError(
            f"--languages must be a JSON array, got {type(lang_list).__name__}"
        )
    for lang in lang_list:
        console.rule(f"[bold]{lang}[/bold]")
        if fix:
            _sync_fix(
                lang,
                include=None,
                parallel=DEFAULT_PARALLEL,
                log_file=None,
                console=console,
            )
        else:
            _sync_normal(
                lang,
                include=None,
                since=None,
                parallel=DEFAULT_PARALLEL,
                log_file=None,
                console=console,
            )


@app.command("ci-pr")
def ci_pr(
    github_token: str = typer.Option(..., envvar="GITHUB_TOKEN"),
    github_repository: str = typer.Option(..., envvar="GITHUB_REPOSITORY"),
    run_url: str = typer.Option(
        "", envvar="GITHUB_RUN_URL", help="GitHub Actions run URL"
    ),
    trigger: str = typer.Option(
        "", help="What triggered this run (commit subject or 'manual')"
    ),
    commit_sha: str = typer.Option(
        "", envvar="GITHUB_SHA", help="Git commit SHA that triggered this run"
    ),
    commit_message: str = typer.Option(
        "",
        envvar="TRIGGER_COMMIT_MESSAGE",
        help="Full commit message that triggered this run",
    ),
) -> None:
    """Create PR with translation changes (CI)."""
    from github import Github

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

    # Stage docs/ and check if anything was actually added
    subprocess.run(["git", "add", "docs/"], check=True, cwd=REPO_ROOT)
    result = subprocess.run(
        ["git", "diff", "--cached", "--quiet"],
        cwd=REPO_ROOT,
    )
    if result.returncode == 0:
        print("No changes to commit")
        return

    branch = f"translate-{secrets.token_hex(4)}"
    subprocess.run(["git", "checkout", "-b", branch], check=True, cwd=REPO_ROOT)
    subprocess.run(
        ["git", "commit", "-m", "Update translations"], check=True, cwd=REPO_ROOT
    )
    subprocess.run(["git", "push", "origin", branch], check=True, cwd=REPO_ROOT)

    # Build PR title
    if trigger:
        trigger_short = trigger[:150] + "..." if len(trigger) > 150 else trigger
        title = f"Update translations ({trigger_short})"
    else:
        title = "Update translations"

    # Build PR body
    if not commit_sha:
        raise ConfigError(
            "GITHUB_SHA not set. Required to record translation baseline in PR body."
        )
    trigger_line = f"commit {commit_sha[:10]}"
    commit_quote = f"\n\n> {commit_message}" if commit_message else ""

    body = (
        f"Automated translation update, generated by workflow run {run_url}\n\n"
        f"To improve translations, see "
        f"[TRANSLATING.md](https://github.com/{github_repository}/blob/master/TRANSLATING.md) "
        f"or edit prompt files:\n"
        f"- `_scripts/general-llm-prompt.md`\n"
        f"- `docs/*/llm-prompt.md`\n\n"
        f"Update triggered by: {trigger_line}{commit_quote}"
    )

    gh = Github(github_token)
    pr = gh.get_repo(github_repository).create_pull(
        title=title,
        body=body,
        base="master",
        head=branch,
    )
    print(f"Created PR: {pr.html_url}")
