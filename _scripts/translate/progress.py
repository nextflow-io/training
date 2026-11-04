"""Progress tracking for parallel translation runs."""

from __future__ import annotations

import asyncio
import time

from .config import make_console
from .models import TranslationFile


def _format_lines(n: int) -> str:
    """Format line count with k suffix for thousands."""
    if n >= 1000:
        return f"{n / 1000:.1f}k"
    return str(n)


class TranslationProgress:
    """Track translation progress across parallel tasks."""

    def __init__(self, files: list[TranslationFile], file_lines: dict[str, int]):
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

    async def start_one(self, filename: str) -> None:
        async with self._lock:
            self.queued.remove(filename)
            self.working.append(filename)

    async def finish_one(self, filename: str, success: bool = True) -> None:
        async with self._lock:
            self.working.remove(filename)
            if success:
                self.complete += 1
                self.lines_complete += self._file_lines[filename]
            else:
                self.failed += 1

    def status(self) -> str:
        """Return a Rich-formatted status line."""
        elapsed = int(time.time() - self._start_time)
        mins, secs = divmod(elapsed, 60)
        lines_remaining = self.lines_total - self.lines_complete
        w = len(str(self.total))
        line = (
            f"[dim][{mins:02d}:{secs:02d}][/dim] "
            f"[green]{self.complete:>{w}}/{self.total}[/green] files complete, "
            f"[cyan]{_format_lines(lines_remaining):>5}[/cyan] lines remaining, "
            f"[yellow]{len(self.working):>{w}}[/yellow] files underway"
        )
        if self.failed:
            line += f", [red]{self.failed} failed[/red]"
        return line


async def progress_logger(
    progress: TranslationProgress, interval: float = 10.0
) -> None:
    """Print progress status at regular intervals."""
    console = make_console()
    console.print(
        f"[dim][00:00][/dim] Starting: [bold]{progress.total}[/bold] files, "
        f"[bold]{_format_lines(progress.lines_total)}[/bold] lines"
    )
    while progress.complete + progress.failed < progress.total:
        await asyncio.sleep(interval)
        console.print(progress.status())
