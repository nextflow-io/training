"""Data structures used across the translation pipeline."""

from __future__ import annotations

import asyncio
import json
from dataclasses import asdict, dataclass
from datetime import datetime
from pathlib import Path
from typing import NamedTuple

from .config import EN_DOCS


# ---------------------------------------------------------------------------
# Translation file descriptor
# ---------------------------------------------------------------------------


@dataclass
class TranslationFile:
    """A file that needs translation."""

    en_path: Path
    lang_path: Path
    language: str
    prompt_changed: bool = False  # True if update triggered by prompt change
    force_full: bool = False  # True to force full re-translation (skip incremental)

    @property
    def exists(self) -> bool:
        return self.lang_path.exists()

    @property
    def relative_path(self) -> Path:
        return self.en_path.relative_to(EN_DOCS)


# ---------------------------------------------------------------------------
# Structural elements (used by post-processing and verification)
# ---------------------------------------------------------------------------


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


# ---------------------------------------------------------------------------
# API result
# ---------------------------------------------------------------------------


@dataclass
class TranslationResult:
    """Result from Claude API including metadata."""

    text: str
    model: str
    input_tokens: int
    output_tokens: int
    stop_reason: str
    continuations: int


# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------


@dataclass
class FileLogEntry:
    """Log entry for a single file translation."""

    filename: str
    started_at: str
    finished_at: str = ""
    duration_s: float = 0.0
    input_lines: int = 0
    input_hash: str = ""
    output_lines: int = 0
    output_hash: str = ""
    changed: bool = False
    error: str = ""
    # API response metadata
    model: str = ""
    input_tokens: int = 0
    output_tokens: int = 0
    stop_reason: str = ""
    continuations: int = 0


class TranslationLog:
    """Verbose log for debugging translation runs."""

    def __init__(
        self,
        log_path: Path,
        baseline: str,
        language: str,
        parallel: int,
        total_files: int,
    ):
        self.path = log_path
        self.baseline = baseline
        self.language = language
        self.parallel = parallel
        self.total_files = total_files
        self.entries: list[FileLogEntry] = []
        self._lock = asyncio.Lock()
        self.started_at = datetime.now().isoformat()

    async def add_entry(self, entry: FileLogEntry) -> None:
        async with self._lock:
            self.entries.append(entry)

    def write(self) -> None:
        entries = sorted(self.entries, key=lambda x: x.started_at)
        for e in entries:
            e.duration_s = round(e.duration_s, 2)
        data = {
            "started_at": self.started_at,
            "finished_at": datetime.now().isoformat(),
            "baseline": self.baseline,
            "language": self.language,
            "parallel": self.parallel,
            "total_files": self.total_files,
            "files_changed": sum(1 for e in entries if e.changed),
            "files_unchanged": sum(1 for e in entries if not e.changed and not e.error),
            "files_failed": sum(1 for e in entries if e.error),
            "total_input_tokens": sum(e.input_tokens for e in entries),
            "total_output_tokens": sum(e.output_tokens for e in entries),
            "files": [asdict(e) for e in entries],
        }
        self.path.write_text(json.dumps(data, indent=2), encoding="utf-8")
