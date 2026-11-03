"""Shared model and helpers for the MkDocs -> Starlight converter.

A rule is a module in `rules/` with `apply(page: Page) -> None` that edits
`page.meta` and `page.body` in place. Rules run in the order of `rules.RULES`.
Every rule must be deterministic and must leave text it does not own untouched.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

import yaml

LOCALES = ["en", "ca", "de", "es", "fr", "hi", "it", "ko", "pl", "pt", "tr"]

FENCE_RE = re.compile(r"^(?P<indent>[ \t]*)(?P<fence>`{3,}|~{3,})(?P<info>.*)$")


@dataclass
class Run:
    """State shared by every page in one converter run."""

    src_root: Path  # MkDocs checkout (worktree of master)
    out_root: Path  # docs/src/content/docs
    # Published path under out_root -> source file. Rules add images and other assets here.
    assets: dict[str, Path] = field(default_factory=dict)
    # Published path under out_root -> callable(bytes) -> bytes, applied when copying.
    asset_transforms: dict[str, object] = field(default_factory=dict)
    warnings: list[str] = field(default_factory=list)

    def docs_dir(self, lang: str) -> Path:
        return self.src_root / "docs" / lang / "docs"


@dataclass
class Page:
    run: Run
    lang: str
    rel: str  # source path relative to the locale's docs dir, e.g. "hello_nextflow/index.md"
    meta: dict
    body: str
    imports: list[str] = field(default_factory=list)  # explicit MDX import lines, rare

    @property
    def src(self) -> Path:
        return self.run.docs_dir(self.lang) / self.rel

    @property
    def out_rel(self) -> str:
        """Output path under the content root: English at the root, other locales in `<lang>/`."""
        rel = str(Path(self.rel).with_suffix(".mdx"))
        return rel if self.lang == "en" else f"{self.lang}/{rel}"

    def warn(self, msg: str) -> None:
        self.run.warnings.append(f"{self.lang}/{self.rel}: {msg}")

    def render(self) -> str:
        out = "---\n" + dump_yaml(self.meta) + "---\n"
        if self.imports:
            out += "\n" + "\n".join(sorted(set(self.imports))) + "\n"
        body = self.body.strip("\n")
        return out + ("\n" + body + "\n" if body else "")


def dump_yaml(meta: dict) -> str:
    return yaml.safe_dump(meta, sort_keys=False, allow_unicode=True, width=10_000, default_flow_style=False)


def split_frontmatter(text: str) -> tuple[dict, str]:
    if text.startswith("---\n"):
        end = text.find("\n---\n", 4)
        if end != -1:
            return (yaml.safe_load(text[4:end]) or {}), text[end + 5 :]
    return {}, text


# --- line helpers -------------------------------------------------------------


def indent_of(line: str) -> int:
    return len(line) - len(line.lstrip(" "))


def dedent(lines: list[str], n: int) -> list[str]:
    return [l[n:] if l[:n].strip() == "" else l.lstrip(" ") for l in lines]


def indent(lines: list[str], n: int) -> list[str]:
    pad = " " * n
    return [pad + l if l.strip() else "" for l in lines]


@dataclass
class Fence:
    start: int  # line index of the opening fence
    end: int  # line index of the closing fence (inclusive); len(lines)-1 if unclosed
    indent: int
    marker: str  # e.g. "```"
    info: str  # everything after the marker, stripped: "groovy title=\"x.nf\" ..."

    @property
    def lang(self) -> str:
        return self.info.split(" ", 1)[0] if self.info else ""

    @property
    def meta(self) -> str:
        return self.info.split(" ", 1)[1].strip() if " " in self.info else ""


def fences(lines: list[str]) -> list[Fence]:
    """Every fenced code block, at any indentation. Fences inside fences are not reported."""
    out, i = [], 0
    while i < len(lines):
        m = FENCE_RE.match(lines[i])
        if not m:
            i += 1
            continue
        marker = m["fence"]
        j = i + 1
        while j < len(lines):
            s = lines[j].strip()
            if s.startswith(marker[0] * len(marker)) and s.strip(marker[0]) == "":
                break
            j += 1
        out.append(Fence(i, min(j, len(lines) - 1), len(m["indent"]), marker, m["info"].strip()))
        i = j + 1
    return out


def fenced_lines(lines: list[str]) -> set[int]:
    """Indexes of every line that belongs to a fenced code block, fences included."""
    s: set[int] = set()
    for f in fences(lines):
        s.update(range(f.start, f.end + 1))
    return s


@dataclass
class Container:
    """An MkDocs indented block: a header line, then a body indented 4 more spaces."""

    start: int
    end: int  # last body line (inclusive), trailing blank lines excluded
    indent: int
    match: re.Match
    body: list[str]  # dedented to column 0


def containers(lines: list[str], header: re.Pattern) -> list[Container]:
    """Top-most containers whose header line matches `header` (matched on the stripped line)."""
    in_fence = fenced_lines(lines)
    out, i = [], 0
    while i < len(lines):
        m = header.match(lines[i].strip()) if i not in in_fence else None
        if not m:
            i += 1
            continue
        ind = indent_of(lines[i])
        j, last = i + 1, i
        while j < len(lines):
            l = lines[j]
            if l.strip() and indent_of(l) < ind + 4:
                break
            if l.strip():
                last = j
            j += 1
        body = dedent(lines[i + 1 : last + 1], ind + 4)
        while body and not body[0].strip():
            body.pop(0)
        out.append(Container(i, last, ind, m, body))
        i = last + 1
    return out


def replace_lines(lines: list[str], start: int, end: int, new: list[str]) -> list[str]:
    return lines[:start] + new + lines[end + 1 :]


def apply_containers(body: str, header: re.Pattern, fn) -> str:
    """Replace each top-most container with fn(container) -> list of lines at column 0.

    The result is re-indented to the container's indent. Runs to a fixed point,
    so containers nested in a replaced body are handled too.
    """
    lines = body.split("\n")
    while True:
        found = containers(lines, header)
        if not found:
            return "\n".join(lines)
        for c in reversed(found):
            lines = replace_lines(lines, c.start, c.end, indent(fn(c), c.indent))
