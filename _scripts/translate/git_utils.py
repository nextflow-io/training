"""Git operations and file-change detection."""

from __future__ import annotations

import os
import re
from functools import lru_cache
from pathlib import Path

import git

from .config import DOCS_ROOT, REPO_ROOT, SCRIPTS_DIR, ConfigError
from .models import TranslationFile
from .paths import en_to_lang_path, iter_en_docs, lang_to_en_path


@lru_cache
def _get_repo() -> git.Repo:
    return git.Repo(REPO_ROOT)


def file_changed_since(path: Path, since_commit: str) -> bool:
    """Check if file has any commits after the given commit."""
    commits = list(
        _get_repo().iter_commits(
            rev=f"{since_commit}..HEAD", paths=str(path), max_count=1
        )
    )
    return len(commits) > 0


def get_file_diff(path: Path, since_commit: str) -> str | None:
    """Get the git diff for a file since a specific commit.

    Returns the unified diff string, or None if no changes or
    diff is too large (>20% of current file lines).
    """
    repo = _get_repo()
    rel_path = path.relative_to(REPO_ROOT)

    diff = repo.git.diff(f"{since_commit}..HEAD", "--", str(rel_path))
    if not diff:
        return None

    # Check if diff is too large (>20% of file lines)
    file_lines = path.read_text(encoding="utf-8").count("\n") + 1
    diff_changes = sum(
        1
        for line in diff.split("\n")
        if line.startswith(("+", "-")) and not line.startswith(("+++", "---"))
    )
    if diff_changes > file_lines * 0.2:
        return None

    return diff


@lru_cache
def get_translation_baseline(
    github_token: str | None = None,
    github_repository: str | None = None,
) -> str:
    """Get the commit SHA that the last translation was based on.

    Finds the most recent merged translation PR and extracts the
    trigger commit from its body.

    Raises:
        ConfigError: If no translation PR found or commit cannot be extracted
    """
    from github import Auth, Github

    token = github_token or os.environ.get("GITHUB_TOKEN")
    repo_name = github_repository or os.environ.get("GITHUB_REPOSITORY")

    if not token:
        raise ConfigError(
            "GITHUB_TOKEN not set. Required to find translation baseline. "
            "Use --since to specify a commit manually."
        )
    if not repo_name:
        raise ConfigError(
            "GITHUB_REPOSITORY not set. Required to find translation baseline. "
            "Use --since to specify a commit manually."
        )

    gh = Github(auth=Auth.Token(token))
    repo = gh.get_repo(repo_name)

    prs = repo.get_pulls(state="closed", sort="updated", direction="desc")
    commit_pattern = re.compile(r"commit:?\s*([a-f0-9]{7,40})")

    for pr in prs:
        if not pr.merged:
            continue
        if "Update translations" not in pr.title:
            continue
        if pr.body:
            match = commit_pattern.search(pr.body)
            if match:
                return match.group(1)

    raise ConfigError(
        "Could not find baseline commit from previous translation PRs. "
        "Use --since to specify a commit manually."
    )


def resolve_baseline(since: str | None, console=None) -> str:
    """Resolve the baseline commit SHA, from explicit value or auto-detection.

    Args:
        since: Explicit commit SHA, or None to auto-detect from last PR.
        console: Optional Rich console for status messages.

    Returns:
        The baseline commit SHA.

    Raises:
        ConfigError: If auto-detection fails.
    """
    if since:
        if console:
            console.print(f"[cyan]Using baseline commit (manual):[/cyan] {since}")
        return since

    baseline = get_translation_baseline()
    if console:
        console.print(f"[cyan]Using baseline commit (from last PR):[/cyan] {baseline}")
    return baseline


def prompt_changed_since(lang: str, baseline: str) -> bool:
    """Check if any prompt files affecting this language changed since baseline."""
    general_prompt = SCRIPTS_DIR / "general-llm-prompt.md"
    lang_prompt = DOCS_ROOT / lang / "llm-prompt.md"

    return file_changed_since(general_prompt, baseline) or file_changed_since(
        lang_prompt, baseline
    )


def get_outdated_files(lang: str, baseline: str | None = None) -> list[TranslationFile]:
    """Find translations needing updates since the baseline commit.

    If prompt files changed since baseline, ALL existing translations
    for that language are considered outdated.
    """
    if not baseline:
        return []

    outdated = []
    prompts_changed = prompt_changed_since(lang, baseline)

    for en_path in iter_en_docs():
        lang_path = en_to_lang_path(en_path, lang)
        if not lang_path.exists():
            continue

        if prompts_changed:
            outdated.append(
                TranslationFile(en_path, lang_path, lang, prompt_changed=True)
            )
        elif file_changed_since(en_path, baseline):
            outdated.append(
                TranslationFile(en_path, lang_path, lang, prompt_changed=False)
            )

    return outdated


def get_missing_files(lang: str) -> list[TranslationFile]:
    """Find English files without translations."""
    return [
        TranslationFile(en_path, lang_path, lang)
        for en_path in iter_en_docs()
        if not (lang_path := en_to_lang_path(en_path, lang)).exists()
    ]


def get_orphaned_files(lang: str) -> list[Path]:
    """Find translation files without English source."""
    lang_docs = DOCS_ROOT / lang / "docs"
    if not lang_docs.exists():
        return []
    return [p for p in lang_docs.rglob("*.md") if not lang_to_en_path(p, lang).exists()]


def gather_work(
    lang: str, baseline: str | None = None
) -> tuple[list[Path], list[TranslationFile], list[TranslationFile]]:
    """Gather all work items for a language.

    Returns:
        (orphaned, missing, outdated) tuple.
    """
    orphaned = get_orphaned_files(lang)
    missing = get_missing_files(lang)
    outdated = get_outdated_files(lang, baseline=baseline)
    return orphaned, missing, outdated


def delete_orphans(orphaned: list[Path], console=None, use_git_rm: bool = False) -> int:
    """Delete orphaned translation files.

    Args:
        orphaned: List of file paths to delete.
        console: Optional Rich console for status messages.
        use_git_rm: If True, use `git rm` (for CI). Otherwise use unlink.

    Returns:
        Number of files deleted.
    """
    import subprocess

    if not orphaned:
        return 0

    if console:
        console.print(f"[bold]Removing {len(orphaned)} orphaned...[/bold]")

    for p in orphaned:
        if use_git_rm:
            subprocess.run(["git", "rm", "-f", str(p)], check=True, cwd=REPO_ROOT)
        else:
            p.unlink()
        if console:
            console.print(f"  [red]Removed:[/red] {p.relative_to(REPO_ROOT)}")

    return len(orphaned)
