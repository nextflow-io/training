"""Heading-based chunking for incremental translation.

Splits markdown files at ``## ``-level headings and compares old vs new
English chunks to determine which sections actually changed.  Unchanged
sections reuse the existing translation, avoiding expensive LLM calls.
"""

from __future__ import annotations

import hashlib
import re
from dataclasses import dataclass, field
from difflib import SequenceMatcher

from .postprocess import CODE_FENCE_RE


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------


@dataclass
class Chunk:
    """A section of a markdown file delimited by ``## `` headings."""

    heading: str
    """The heading line itself (e.g. ``## 2. Publish outputs``), or empty
    string for the preamble (chunk 0)."""

    content: str
    """Full text of the chunk **including** the heading line."""

    index: int
    """Position of this chunk in the file (0-based)."""

    content_hash: str = field(default="", repr=False)
    """MD5 hex digest of *content* (set automatically by ``split_into_chunks``)."""


@dataclass
class ChunkDiff:
    """Result of comparing old and new chunk lists."""

    unchanged: list[tuple[int, int]]
    """(old_idx, new_idx) pairs for chunks with identical content."""

    modified: list[tuple[int, int]]
    """(old_idx, new_idx) pairs for chunks whose content changed."""

    added: list[int]
    """new_idx values with no match in the old list."""

    removed: list[int]
    """old_idx values with no match in the new list."""


# ---------------------------------------------------------------------------
# Splitting
# ---------------------------------------------------------------------------


def _content_hash(text: str) -> str:
    """Hash content with trailing whitespace stripped for stable comparison.

    Trailing blank lines shift between chunks when sections are added or
    removed, so we strip them before hashing to avoid false "modified"
    classifications.
    """
    return hashlib.md5(text.rstrip().encode()).hexdigest()[:12]


def split_into_chunks(text: str) -> list[Chunk]:
    """Split markdown *text* into chunks at ``## `` headings.

    * Chunk 0 is the "preamble" — everything before the first ``## ``.
    * Each subsequent chunk starts at a ``## `` line (outside code blocks)
      and includes all content up to (but not including) the next ``## ``.
    * ``## `` lines inside fenced code blocks are ignored.
    """
    lines = text.split("\n")
    chunks: list[Chunk] = []
    current_lines: list[str] = []
    current_heading = ""
    in_code = False
    fence = ""

    def _flush():
        content = "\n".join(current_lines)
        chunks.append(
            Chunk(
                heading=current_heading,
                content=content,
                index=len(chunks),
                content_hash=_content_hash(content),
            )
        )

    for line in lines:
        stripped = line.lstrip()

        # Track code fences
        if not in_code:
            if m := CODE_FENCE_RE.match(stripped):
                in_code, fence = True, m.group(1)
        elif stripped.startswith(fence) and stripped.strip() == fence:
            in_code, fence = False, ""

        # Start a new chunk at each ## heading (outside code blocks)
        if not in_code and line.startswith("## "):
            _flush()
            current_lines = [line]
            current_heading = line
        else:
            current_lines.append(line)

    _flush()
    return chunks


def reassemble_chunks(chunks: list[Chunk]) -> str:
    """Reassemble chunks into a single document string."""
    return "\n".join(c.content for c in chunks)


# ---------------------------------------------------------------------------
# Diffing
# ---------------------------------------------------------------------------


def _normalize_heading(heading: str) -> str:
    """Normalize heading text for comparison (strip whitespace, anchors)."""
    # Remove trailing anchor IDs like { #some-id }
    h = re.sub(r"\s*\{[^}]*\}\s*$", "", heading)
    return h.strip()


def diff_chunks(old_chunks: list[Chunk], new_chunks: list[Chunk]) -> ChunkDiff:
    """Compare old and new chunk lists to identify changes.

    Matching strategy (applied in order):
    1. **Heading match**: exact heading text (after normalization).
    2. **Content hash match**: same content, different heading (detects
       renames/reorders without content changes).
    3. **Fuzzy heading match**: heading similarity >0.8 via SequenceMatcher
       (detects heading renames with content changes).
    """
    unchanged: list[tuple[int, int]] = []
    modified: list[tuple[int, int]] = []

    # Track which chunks have been matched
    matched_old: set[int] = set()
    matched_new: set[int] = set()

    # Build lookup maps
    old_by_heading: dict[str, list[int]] = {}
    for c in old_chunks:
        key = _normalize_heading(c.heading)
        old_by_heading.setdefault(key, []).append(c.index)

    old_by_hash: dict[str, list[int]] = {}
    for c in old_chunks:
        old_by_hash.setdefault(c.content_hash, []).append(c.index)

    # Pass 1: exact heading match
    for nc in new_chunks:
        key = _normalize_heading(nc.heading)
        candidates = old_by_heading.get(key, [])
        for oi in candidates:
            if oi not in matched_old:
                matched_old.add(oi)
                matched_new.add(nc.index)
                if old_chunks[oi].content_hash == nc.content_hash:
                    unchanged.append((oi, nc.index))
                else:
                    modified.append((oi, nc.index))
                break

    # Pass 2: content hash match (for reordered/renamed chunks)
    for nc in new_chunks:
        if nc.index in matched_new:
            continue
        candidates = old_by_hash.get(nc.content_hash, [])
        for oi in candidates:
            if oi not in matched_old:
                matched_old.add(oi)
                matched_new.add(nc.index)
                unchanged.append((oi, nc.index))
                break

    # Pass 3: fuzzy heading match
    unmatched_old = [c for c in old_chunks if c.index not in matched_old]
    unmatched_new = [c for c in new_chunks if c.index not in matched_new]

    for nc in list(unmatched_new):
        best_score = 0.0
        best_oc = None
        nc_heading = _normalize_heading(nc.heading)
        for oc in unmatched_old:
            oc_heading = _normalize_heading(oc.heading)
            score = SequenceMatcher(None, oc_heading, nc_heading).ratio()
            if score > best_score:
                best_score = score
                best_oc = oc
        if best_oc is not None and best_score > 0.8:
            matched_old.add(best_oc.index)
            matched_new.add(nc.index)
            unmatched_old.remove(best_oc)
            unmatched_new.remove(nc)
            if best_oc.content_hash == nc.content_hash:
                unchanged.append((best_oc.index, nc.index))
            else:
                modified.append((best_oc.index, nc.index))

    # Remaining unmatched
    added = [nc.index for nc in unmatched_new]
    removed = [oc.index for oc in unmatched_old]

    return ChunkDiff(
        unchanged=unchanged,
        modified=modified,
        added=added,
        removed=removed,
    )
