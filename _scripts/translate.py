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
    uv run translate.py sync pt --parallel 20
    uv run translate.py translate docs/en/docs/index.md --lang pt
"""

from __future__ import annotations

import asyncio
import json
import os
import random
import re
import secrets
import subprocess
import sys
import time
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import NamedTuple

import anthropic
import git
import typer
import yaml
from rich.console import Console

# Module-level console for output (force color in CI)
console = Console(force_terminal=True if os.getenv("GITHUB_ACTIONS") else None)


# =============================================================================
# Configuration
# =============================================================================

REPO_ROOT = Path(__file__).parent.parent
DOCS_ROOT = REPO_ROOT / "docs"
EN_DOCS = DOCS_ROOT / "en" / "docs"
SCRIPTS_DIR = REPO_ROOT / "_scripts"

MODEL = "claude-sonnet-4-5"
MAX_TOKENS = 32768  # Large enough for biggest docs (~60KB source)
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


def file_changed_since(path: Path, since_commit: str) -> bool:
    """Check if file has any commits after the given commit.

    Args:
        path: File path to check
        since_commit: Commit SHA to compare against

    Returns:
        True if the file has commits newer than since_commit
    """
    try:
        # git log since_commit..HEAD -- path
        commits = list(
            _get_repo().iter_commits(
                rev=f"{since_commit}..HEAD", paths=str(path), max_count=1
            )
        )
        return len(commits) > 0
    except git.GitCommandError:
        return False


def get_translation_baseline(
    github_token: str | None = None,
    github_repository: str | None = None,
) -> str:
    """Get the commit SHA that the last translation was based on.

    Finds the most recent merged translation PR and extracts the
    trigger commit from its body.

    Args:
        github_token: GitHub token for API access (uses GITHUB_TOKEN env if not provided)
        github_repository: Repository in owner/repo format (uses GITHUB_REPOSITORY env if not provided)

    Returns:
        The commit SHA that triggered the last translation

    Raises:
        ConfigError: If no translation PR found or commit cannot be extracted
    """
    from github import Github

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

    from github import Auth

    gh = Github(auth=Auth.Token(token))
    repo = gh.get_repo(repo_name)

    # Find the most recent merged translation PR
    # Search for PRs with "Update translations" in title, merged, sorted by updated desc
    prs = repo.get_pulls(state="closed", sort="updated", direction="desc")

    commit_pattern = re.compile(r"commit:?\s*([a-f0-9]{7,40})")

    for pr in prs:
        if not pr.merged:
            continue
        if "Update translations" not in pr.title:
            continue

        # Extract commit SHA from PR body
        if pr.body:
            match = commit_pattern.search(pr.body)
            if match:
                return match.group(1)

        # If we found a translation PR but couldn't extract commit, keep looking
        # (older PRs might have different format)

    raise ConfigError(
        "Could not find baseline commit from previous translation PRs. "
        "Use --since to specify a commit manually."
    )


def prompt_changed_since(lang: str, baseline: str) -> bool:
    """Check if any prompt files affecting this language changed since baseline.

    Returns True if either:
    - The general prompt (_scripts/general-llm-prompt.md) changed
    - The language-specific prompt (docs/{lang}/llm-prompt.md) changed
    """
    general_prompt = SCRIPTS_DIR / "general-llm-prompt.md"
    lang_prompt = DOCS_ROOT / lang / "llm-prompt.md"

    if file_changed_since(general_prompt, baseline):
        return True
    if file_changed_since(lang_prompt, baseline):
        return True
    return False


def get_outdated_files(lang: str, baseline: str | None = None) -> list[TranslationFile]:
    """Find translations needing updates.

    Args:
        lang: Target language code
        baseline: Commit SHA to compare against. Finds files that changed
            in English since this commit. If None, returns empty list.

    Returns:
        List of TranslationFile objects needing updates

    Note:
        If prompt files (general or language-specific) changed since baseline,
        ALL existing translations for that language are considered outdated.
    """
    if not baseline:
        return []

    outdated = []

    # Check if prompts changed - if so, all existing translations are outdated
    prompts_changed = prompt_changed_since(lang, baseline)

    for en_path in iter_en_docs():
        lang_path = en_to_lang_path(en_path, lang)
        if not lang_path.exists():
            continue

        if prompts_changed:
            # Prompt changed: all existing translations need re-translation
            outdated.append(TranslationFile(en_path, lang_path, lang))
        elif file_changed_since(en_path, baseline):
            # Check if English file changed since baseline
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


MAX_CONTINUATIONS = 5  # Max continuation requests for very large files


def _call_claude_once(
    client: anthropic.Anthropic,
    messages: list[dict],
    console: Console | None = None,
) -> anthropic.types.Message:
    """Make a single Claude API call with retry logic for transient errors.

    Returns the raw Message object so caller can check stop_reason.
    """
    for attempt in range(MAX_RETRIES + 1):
        try:
            return client.messages.create(
                model=MODEL,
                max_tokens=MAX_TOKENS,
                timeout=REQUEST_TIMEOUT,
                messages=messages,
            )
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


def call_claude(prompt: str, console: Console | None = None) -> str:
    """Call Claude API with automatic continuation for large responses.

    If the response hits max_tokens, automatically continues the conversation
    to get the complete output. Includes retry logic with exponential backoff
    for transient errors.
    """
    client = anthropic.Anthropic()
    messages: list[dict] = [{"role": "user", "content": prompt}]

    # Initial request
    message = _call_claude_once(client, messages, console)
    result_parts = [message.content[0].text]

    # Handle continuations if response was truncated
    continuation = 0
    while message.stop_reason == "max_tokens":
        continuation += 1
        if continuation > MAX_CONTINUATIONS:
            raise TranslationError(
                f"Response still incomplete after {MAX_CONTINUATIONS} continuations. "
                f"File may be too large to translate."
            )
        if console:
            console.print(
                f"    [yellow]Response truncated, requesting continuation "
                f"{continuation}/{MAX_CONTINUATIONS}...[/yellow]"
            )

        # Append assistant response and continuation request to conversation history
        messages.append({"role": "assistant", "content": message.content[0].text})
        messages.append(
            {"role": "user", "content": "Continue exactly where you left off."}
        )

        message = _call_claude_once(client, messages, console)
        result_parts.append(message.content[0].text)

    return "".join(result_parts).strip()


# =============================================================================
# Async Translation API (for parallel execution)
# =============================================================================


async def _call_claude_once_async(
    client: anthropic.AsyncAnthropic,
    messages: list[dict],
    file_label: str = "",
) -> anthropic.types.Message:
    """Make a single async Claude API call with jittered exponential backoff.

    Args:
        client: Async Anthropic client
        messages: Message history
        file_label: Optional label for log messages (e.g., file path)
    """
    for attempt in range(MAX_RETRIES + 1):
        try:
            return await client.messages.create(
                model=MODEL,
                max_tokens=MAX_TOKENS,
                timeout=REQUEST_TIMEOUT,
                messages=messages,
            )
        except (
            anthropic.APIConnectionError,
            anthropic.APITimeoutError,
            anthropic.RateLimitError,
            anthropic.InternalServerError,
        ) as e:
            if attempt == MAX_RETRIES:
                raise
            # Jittered exponential backoff to avoid thundering herd
            base_delay = BASE_DELAY * (2**attempt)
            jitter = random.uniform(0, base_delay * 0.5)
            delay = base_delay + jitter
            print(
                f"  [{file_label}] Retry {attempt + 1}/{MAX_RETRIES} "
                f"after {delay:.1f}s: {type(e).__name__}"
            )
            await asyncio.sleep(delay)

    raise RuntimeError("Unreachable")


async def call_claude_async(
    prompt: str, file_label: str = "", client: anthropic.AsyncAnthropic | None = None
) -> str:
    """Async version of call_claude with continuation support.

    Args:
        prompt: The translation prompt
        file_label: Optional label for log messages
        client: Optional shared client (creates one if not provided)
    """
    client = client or anthropic.AsyncAnthropic()
    messages: list[dict] = [{"role": "user", "content": prompt}]

    # Initial request
    message = await _call_claude_once_async(client, messages, file_label)
    result_parts = [message.content[0].text]

    # Handle continuations if response was truncated
    continuation = 0
    while message.stop_reason == "max_tokens":
        continuation += 1
        if continuation > MAX_CONTINUATIONS:
            raise TranslationError(
                f"Response still incomplete after {MAX_CONTINUATIONS} continuations. "
                f"File may be too large to translate."
            )
        print(
            f"  [{file_label}] Response truncated, continuation "
            f"{continuation}/{MAX_CONTINUATIONS}..."
        )

        messages.append({"role": "assistant", "content": message.content[0].text})
        messages.append(
            {"role": "user", "content": "Continue exactly where you left off."}
        )

        message = await _call_claude_once_async(client, messages, file_label)
        result_parts.append(message.content[0].text)

    return "".join(result_parts).strip()


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


class TranslationProgress:
    """Track translation progress across parallel tasks."""

    def __init__(self, files: list["TranslationFile"], file_lines: dict[str, int]):
        self.total = len(files)
        self.queued: list[str] = [str(f.relative_path) for f in files]
        self.working: list[str] = []
        self.complete = 0
        self.failed = 0
        self._file_lines = file_lines
        self.lines_total = sum(file_lines.values())
        self.lines_complete = 0
        self._lock = asyncio.Lock()
        self._start_time = time.time()

    async def start_one(self, filename: str):
        async with self._lock:
            self.queued.remove(filename)
            self.working.append(filename)

    async def finish_one(self, filename: str, success: bool = True):
        async with self._lock:
            self.working.remove(filename)
            if success:
                self.complete += 1
                self.lines_complete += self._file_lines[filename]
            else:
                self.failed += 1

    def status(self) -> str:
        """Return a status line with aligned columns, using Rich colors."""
        elapsed = int(time.time() - self._start_time)
        mins, secs = divmod(elapsed, 60)
        lines_remaining = self.lines_total - self.lines_complete
        w = len(str(self.total))  # width for padding
        line = (
            f"[dim][{mins:02d}:{secs:02d}][/dim] "
            f"[green]{self.complete:>{w}}/{self.total}[/green] files complete, "
            f"[cyan]{_format_lines(lines_remaining):>5}[/cyan] lines remaining, "
            f"[yellow]{len(self.working):>{w}}[/yellow] files underway"
        )
        if self.failed:
            line += f", [red]{self.failed} failed[/red]"
        return line


def _format_lines(n: int) -> str:
    """Format line count with k suffix for thousands."""
    if n >= 1000:
        return f"{n / 1000:.1f}k"
    return str(n)


async def _progress_logger(progress: TranslationProgress, interval: float = 10.0):
    """Print progress status every interval seconds, starting immediately."""
    console.print(
        f"[dim][00:00][/dim] Starting: [bold]{progress.total}[/bold] files, [bold]{_format_lines(progress.lines_total)}[/bold] lines"
    )
    while progress.complete + progress.failed < progress.total:
        await asyncio.sleep(interval)
        console.print(progress.status())


async def translate_file_async(
    tf: TranslationFile,
    semaphore: asyncio.Semaphore,
    progress: TranslationProgress,
    client: anthropic.AsyncAnthropic,
) -> None:
    """Translate a single file (async version for parallel execution)."""
    filename = str(tf.relative_path)
    async with semaphore:
        await progress.start_one(filename)
        try:
            langs = get_languages()
            en_content = tf.en_path.read_text(encoding="utf-8")
            existing = tf.lang_path.read_text(encoding="utf-8") if tf.exists else None

            prompt = build_translation_prompt(
                tf.language, langs[tf.language], en_content, existing
            )
            result = await call_claude_async(prompt, filename, client)

            tf.lang_path.parent.mkdir(parents=True, exist_ok=True)
            tf.lang_path.write_text(f"{result}\n", encoding="utf-8", newline="\n")
            post_process_file(tf.lang_path, tf.language)
            await progress.finish_one(filename, success=True)
        except Exception:
            await progress.finish_one(filename, success=False)
            raise


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

    if lang == "en":
        raise typer.Exit("Cannot translate to English")

    en_path = _resolve_en_path(path)
    tf = TranslationFile(en_path, en_to_lang_path(en_path, lang), lang)

    translate_file(tf, console)
    post_process_file(tf.lang_path, lang)
    console.print(f"[green]Done:[/green] {tf.lang_path}")


DEFAULT_PARALLEL = 50


async def _translate_all(files: list[TranslationFile], parallel: int) -> None:
    """Translate all files in parallel with progress logging."""
    # Pre-compute line counts once (used for sorting and progress tracking)
    file_lines = {
        str(f.relative_path): f.en_path.read_text().count("\n") for f in files
    }

    # Sort largest files first to avoid "long pole" problem
    files.sort(key=lambda f: file_lines[str(f.relative_path)], reverse=True)

    semaphore = asyncio.Semaphore(parallel)
    progress = TranslationProgress(files, file_lines)
    client = anthropic.AsyncAnthropic()

    # Create translation tasks
    tasks = [translate_file_async(tf, semaphore, progress, client) for tf in files]

    # Run translations with progress logger
    logger_task = asyncio.create_task(_progress_logger(progress))
    try:
        results = await asyncio.gather(*tasks, return_exceptions=True)
    finally:
        logger_task.cancel()
        try:
            await logger_task
        except asyncio.CancelledError:
            pass

    # Print final status
    console.print(f"[bold green]Translation complete:[/bold green] {progress.status()}")

    # Raise if any failed
    errors = [r for r in results if isinstance(r, Exception)]
    if errors:
        for e in errors:
            console.print(f"[red]Error: {e}[/red]")
        raise errors[0]


@app.command()
def sync(
    lang: str = typer.Argument(..., help="Language code"),
    include: str | None = typer.Option(None, "--include", "-i", help="Filter pattern"),
    since: str | None = typer.Option(
        None,
        "--since",
        help="Compare since this commit (default: auto-detect from last translation PR)",
    ),
    dry_run: bool = typer.Option(False, "--dry-run", "-n", help="Preview only"),
    parallel: int = typer.Option(
        DEFAULT_PARALLEL,
        "--parallel",
        "-p",
        help="Max concurrent translations (default: 50)",
    ),
):
    """Sync translations: update outdated, add missing, remove orphaned."""
    console = Console(force_terminal=True if os.getenv("GITHUB_ACTIONS") else None)

    # Determine baseline commit
    if since:
        baseline = since
        console.print(f"[cyan]Using baseline commit (manual):[/cyan] {baseline}")
    else:
        try:
            baseline = get_translation_baseline()
            console.print(
                f"[cyan]Using baseline commit (from last PR):[/cyan] {baseline}"
            )
        except ConfigError as e:
            console.print(f"[red]Error:[/red] {e}")
            raise typer.Exit(1)

    # Gather work
    orphaned = get_orphaned_files(lang)
    outdated = get_outdated_files(lang, baseline=baseline)
    missing = get_missing_files(lang)

    # Check if prompt changed (for informational output)
    prompts_changed = baseline and prompt_changed_since(lang, baseline)
    if prompts_changed:
        console.print(
            f"[magenta]Prompt changed:[/magenta] all {lang} translations will be updated"
        )

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
            if not prompts_changed:
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

    # Run translations in parallel
    all_files = outdated + missing
    if all_files:
        print(
            f"Translating {len(all_files)} files with parallel={parallel}. Largest files first."
        )
        asyncio.run(_translate_all(all_files, parallel))

    console.print("[green]✓ Sync complete[/green]")


# =============================================================================
# CI Commands
# =============================================================================


@app.command("ci-detect")
def ci_detect(
    language: str | None = typer.Option(None, "--language"),
    since: str | None = typer.Option(
        None,
        "--since",
        help="Compare since this commit (default: auto-detect from last translation PR)",
    ),
):
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

    # Determine baseline commit
    if since:
        baseline = since
        console.print(f"[cyan]Using baseline commit (manual):[/cyan] {baseline}")
    else:
        try:
            baseline = get_translation_baseline()
            console.print(
                f"[cyan]Using baseline commit (from last PR):[/cyan] {baseline}"
            )
        except ConfigError as e:
            console.print(f"[red]Error:[/red] {e}")
            raise typer.Exit(1)

    # Check which have work and collect details
    need_sync = []
    for lang in langs:
        missing = get_missing_files(lang)
        outdated = get_outdated_files(lang, baseline=baseline)
        orphaned = get_orphaned_files(lang)

        # Check if prompt changed (for informational output)
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
                    # Only list individual files if not all files are outdated due to prompt change
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
):
    """Create PR with translation changes (CI)."""
    from github import Github

    subprocess.run(
        ["git", "config", "user.name", "github-actions[bot]"], check=True, cwd=REPO_ROOT
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
        # Nothing staged - exit cleanly
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
    # Always include commit SHA so future runs can detect baseline
    if not commit_sha:
        raise ConfigError(
            "GITHUB_SHA not set. Required to record translation baseline in PR body."
        )
    trigger_line = f"commit {commit_sha[:10]}"
    commit_quote = f"\n\n> {commit_message}" if commit_message else ""

    body = f"""Automated translation update, generated by workflow run {run_url}

To improve translations, see [TRANSLATING.md](https://github.com/{github_repository}/blob/master/TRANSLATING.md) or edit prompt files:
- `_scripts/general-llm-prompt.md`
- `docs/*/llm-prompt.md`

Update triggered by: {trigger_line}{commit_quote}"""

    gh = Github(github_token)
    pr = gh.get_repo(github_repository).create_pull(
        title=title,
        body=body,
        base="master",
        head=branch,
    )
    print(f"Created PR: {pr.html_url}")


if __name__ == "__main__":
    if not DOCS_ROOT.exists():
        sys.exit("Error: Must run from repository root")
    app()
