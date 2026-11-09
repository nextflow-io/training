"""Links: page links point at the .mdx files, and #anchors use Starlight's heading ids.

MkDocs slugged headings with python-markdown (ASCII only, runs of `-` and spaces
collapse, repeats get `_1`); Starlight uses github-slugger on the rendered text.
English headings whose two slugs differ get a `{#mkdocs-slug}` heading attribute,
so every English heading keeps its MkDocs id and inbound links keep working.
Translated headings keep their github-slugger ids (accepted drift).

Translated pages often link to the English anchor of a translated heading, which
never worked; such anchors map to the heading with the same number ("2.3."), or
else at the same position when both pages have the same heading levels.
"""

import html
import os
import re
import sys
import unicodedata
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

from core import Page, Run, fenced_lines, split_frontmatter

LINK = re.compile(r"(?P<pre>\]\(|href=\")(?P<path>(?![a-z][a-z\d+.-]*:|/)[^)\s\"#]*)(?P<hash>#[^)\s\"]*)?(?=[)\"])")
DEAD_LINK = re.compile(r"(?<!!)\[(?P<text>[^\[\]]*)\]\((?P<path>(?![a-z][a-z\d+.-]*:|/|#|<)[^){}\s#]+)(?:#[^)\s]*)?\)")
HEADING = re.compile(r"^[ \t]*(?P<hashes>#{1,6})[ \t]+(?P<text>.*?)(?:[ \t]+#+)?[ \t]*$")
NUMBER = re.compile(r"^(\d+(?:\.\d+)*)\.?\s")
ATTR_LIST = re.compile(r"\s*\{[^}]*\}$")
CODE = re.compile(r"(`+)(.+?)(?<!`)\1(?!`)")
INLINE_LANG = re.compile(r"^([^`]*\S)\{:[\w#+-][\w#+.-]*\}$")
SHORTCODE = re.compile(r":[a-z0-9_+-]+:(\{[^}]*\})?")
MKDOCS_REPEAT = re.compile(r"^(.*)_(\d+)$")


@dataclass
class Heading:
    line: int
    level: int
    text: str  # inline Markdown
    number: str | None


@dataclass
class Target:
    level: int
    number: str | None
    mkdocs: str  # id the heading had in MkDocs
    site: str  # id it has in Starlight
    gh: str  # github-slugger id, which some links were written for
    pinned: bool


def plain(md: str, code, mkdocs: bool) -> str:
    """Text content of inline Markdown, as the heading slugger sees it."""
    spans: list[str] = []

    def stash(m):
        spans.append(code(m[2]))
        return f"\0{len(spans) - 1}\0"

    s = CODE.sub(stash, md)
    if mkdocs:
        s = SHORTCODE.sub("", s)
    s = re.sub(r"!\[[^\]]*\]\([^)]*\)", "", s)
    s = re.sub(r"\[([^\]]*)\]\([^)]*\)", r"\1", s)
    s = re.sub(r"\{/\*.*?\*/\}", "", s)
    s = re.sub(r"</?[A-Za-z][^>]*>", "", s)
    s = re.sub(r"\*\*|__|~~|\*|(?<!\w)_|_(?!\w)", "", s)
    s = re.sub(r"\\([!-/:-@\[-`{-~])", r"\1", s)
    s = html.unescape(s)
    return re.sub(r"\0(\d+)\0", lambda m: spans[int(m[1])], s)


def mkdocs_text(md: str) -> str:
    return plain(ATTR_LIST.sub("", md), lambda c: re.sub(r"^#!\w+ ", "", c).strip(), mkdocs=True)


def site_text(md: str) -> str:
    return plain(md, lambda c: (m[1] if (m := INLINE_LANG.match(c)) else c), mkdocs=False)


def mkdocs_slug(text: str) -> str:
    text = unicodedata.normalize("NFKD", text).encode("ascii", "ignore").decode("ascii")
    text = re.sub(r"[^\w\s-]", "", text).strip().lower()
    return re.sub(r"[-\s]+", "-", text)


def github_slug(text: str) -> str:
    """github-slugger: keeps letters, marks, digits, `_` and `-`; spaces become `-`."""

    def keep(c: str) -> bool:
        cat = unicodedata.category(c)
        return c in " -" or cat[0] in "LM" or cat in ("Nd", "Nl", "Pc")

    return "".join(c for c in text.lower() if keep(c)).replace(" ", "-")


def mkdocs_ids(texts: list[str]) -> list[str]:
    ids: set[str] = set()
    out = []
    for t in texts:
        id_ = mkdocs_slug(t)
        while id_ in ids or not id_:
            m = MKDOCS_REPEAT.match(id_)
            id_ = f"{m[1]}_{int(m[2]) + 1}" if m else f"{id_}_1"
        ids.add(id_)
        out.append(id_)
    return out


class Slugger:
    """github-slugger's repeat counting. `commit=False` only peeks."""

    def __init__(self):
        self.seen: dict[str, int] = {}

    def slug(self, text: str, commit: bool = True) -> str:
        base = result = github_slug(text)
        n = self.seen.get(base, 0)
        while result in self.seen:
            n += 1
            result = f"{base}-{n}"
        if commit:
            if result != base:
                self.seen[base] = n
            self.seen[result] = 0
        return result


def headings(body: str) -> list[Heading]:
    lines = body.split("\n")
    in_fence = fenced_lines(lines)
    out = []
    for i, line in enumerate(lines):
        if i not in in_fence and (m := HEADING.match(line)):
            n = NUMBER.match(m["text"])
            out.append(Heading(i, len(m["hashes"]), m["text"], n[1] if n else None))
    return out


def mkdocs_ids_of(source: list[Heading], converted: list[Heading]) -> list[str] | None:
    """MkDocs id of each converted heading, or None when conversion changed the headings."""
    old = mkdocs_ids([mkdocs_text(h.text) for h in source])
    # The frontmatter rule drops a leading H1; every other heading survives, in order.
    if len(source) == len(converted) + 1 and source[0].level == 1:
        source, old = source[1:], old[1:]
    return old if [h.level for h in source] == [h.level for h in converted] else None


def targets(lang: str, old: list[str] | None, converted: list[Heading]) -> list[Target]:
    """Old and new id of each converted heading. English headings are pinned to their MkDocs id."""
    slugger, out = Slugger(), []
    for h, md in zip(converted, old or [None] * len(converted)):
        text = site_text(h.text)
        gh = slugger.slug(text, commit=False)
        pin = lang == "en" and md is not None and gh != md
        site = md if pin else slugger.slug(text)
        out.append(Target(h.level, h.number, md or site, site, gh, pin))
    return out


@lru_cache(maxsize=None)
def source_headings(path: Path) -> list[Heading]:
    return headings(split_frontmatter(path.read_text(encoding="utf-8"))[1]) if path.is_file() else []


@lru_cache(maxsize=None)
def page_targets(src_root: Path, lang: str, rel: str) -> list[Target]:
    """Headings of another page as Starlight will see them: converted by the rules before this one."""
    from rules import RULES

    run = Run(src_root=src_root, out_root=Path("/nonexistent"))
    path = run.docs_dir(lang) / rel
    if not path.is_file():
        return []
    meta, body = split_frontmatter(path.read_text(encoding="utf-8"))
    page = Page(run=run, lang=lang, rel=rel, meta=meta, body=body)
    for rule in RULES[: RULES.index(sys.modules[__name__])]:
        rule.apply(page)
    converted = headings(page.body)
    return targets(lang, mkdocs_ids_of(source_headings(path), converted), converted)


def resolve_anchor(target: list[Target], english: list[Target], anchor: str) -> str | None:
    """The new id for `anchor`, matched on the target page, else on the English page by number or position."""
    for t in target:
        if anchor in (t.mkdocs, t.site, t.gh):
            return t.site
    i = next((i for i, t in enumerate(english) if anchor in (t.mkdocs, t.site, t.gh)), None)
    if i is None or english is target:
        return None
    if english[i].number:
        same = [t for t in target if t.number == english[i].number]
        return same[0].site if len(same) == 1 else None
    if [t.level for t in english] == [t.level for t in target]:
        return target[i].site
    return None


def unlink_missing(page, line: str) -> str:
    """A relative link to a file that exists in no locale keeps its text and loses the link."""

    def sub(m):
        rel = os.path.normpath(os.path.join(os.path.dirname(page.rel), m["path"]))
        if any((page.run.docs_dir(l) / rel).exists() for l in (page.lang, "en")):
            return m[0]
        page.warn(f"link {m['path']} matches no file; kept the text only")
        return m["text"]

    return DEAD_LINK.sub(sub, line)


def apply(page) -> None:
    lines = page.body.split("\n")
    converted = headings(page.body)
    # Fixtures have no source file: their body is the source.
    old = mkdocs_ids_of(source_headings(page.src) if page.src.is_file() else converted, converted)
    if old is None:
        page.warn("headings changed in conversion; MkDocs ids not kept")
    own = targets(page.lang, old, converted)
    for h, t in zip(converted, own):
        if t.pinned:
            lines[h.line] = lines[h.line].rstrip() + f" {{#{t.site}}}"
    in_fence = fenced_lines(lines)

    def sub(m):
        path, hash_ = m["path"], m["hash"] or ""
        if (not path and not hash_) or (path and not path.endswith(".md")):
            return m[0]
        rel = os.path.normpath(os.path.join(os.path.dirname(page.rel), path)) if path else page.rel
        lang = page.lang if not path else next((l for l in (page.lang, "en") if (page.run.docs_dir(l) / rel).is_file()), None)
        if lang is None:
            page.warn(f"link {path} matches no page")
            return m[0]
        if len(hash_) > 1:
            target = own if not path else page_targets(page.run.src_root, lang, rel)
            english = target if lang == "en" else page_targets(page.run.src_root, "en", rel)
            new = resolve_anchor(target, english, hash_[1:])
            if new is None:
                page.warn(f"anchor {path}{hash_} matches no heading")
            else:
                hash_ = "#" + new
        return f"{m['pre']}{path[:-3] + '.mdx' if path else ''}{hash_}"

    page.body = "\n".join(l if i in in_fence else LINK.sub(sub, unlink_missing(page, l)) for i, l in enumerate(lines))
