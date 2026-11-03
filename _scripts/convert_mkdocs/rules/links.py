"""Links: page links point at the .mdx files, and #anchors use Starlight's heading slugs.

MkDocs slugged headings with python-markdown (ASCII only, runs of `-` collapsed);
Starlight uses github-slugger. Translated pages often link to the English anchor
of a translated heading, which never worked; the heading number ("2.3.") is the
same in every locale, so such anchors map to the heading with that number.
"""

import os
import re
import unicodedata
from functools import lru_cache
from pathlib import Path

from core import fenced_lines

LINK = re.compile(r"(?P<pre>\]\(|href=\")(?P<path>(?![a-z][a-z\d+.-]*:|/)[^)\s\"#]*)(?P<hash>#[^)\s\"]*)?(?=[)\"])")
HEADING = re.compile(r"^(#{2,6})\s+(?P<text>.+?)\s*$")
NUMBER = re.compile(r"^(\d+(?:\.\d+)*)\.?\s")


def heading_text(md: str) -> str:
    md = re.sub(r"`#!\w+ ([^`]*)`", r"\1", md)
    md = re.sub(r"\{:\w+\}`", "`", md)
    md = re.sub(r"!?\[([^\]]*)\]\([^)]*\)", r"\1", md)
    md = re.sub(r"<[^>]+>", "", md)
    return re.sub(r"[`*]|(?<!\w)_|_(?!\w)", "", md).strip()


def mkdocs_slug(text: str) -> str:
    text = unicodedata.normalize("NFKD", text).encode("ascii", "ignore").decode("ascii")
    text = re.sub(r"[^\w\s-]", "", text).strip().lower()
    return re.sub(r"[-\s]+", "-", text)


def github_slug(text: str) -> str:
    keep = (c for c in text.lower() if c.isalnum() or c in "_- " or unicodedata.category(c).startswith("M"))
    return "".join(keep).replace(" ", "-")


def unique(slugs: list[str], sep: str) -> list[str]:
    seen: dict[str, int] = {}
    out = []
    for s in slugs:
        n = seen.get(s, 0)
        seen[s] = n + 1
        out.append(s if n == 0 else f"{s}{sep}{n}")
    return out


@lru_cache(maxsize=None)
def headings(path: Path) -> list[tuple[str | None, str, str]]:
    """(number, mkdocs slug, github slug) for each heading of a source page."""
    if not path.is_file():
        return []
    lines = path.read_text(encoding="utf-8").split("\n")
    in_fence = fenced_lines(lines)
    texts = [heading_text(m["text"]) for i, l in enumerate(lines) if i not in in_fence and (m := HEADING.match(l))]
    numbers = [(m[1] if (m := NUMBER.match(t)) else None) for t in texts]
    return list(zip(numbers, unique([mkdocs_slug(t) for t in texts], "_"), unique([github_slug(t) for t in texts], "-")))


def resolve_anchor(page, target_rel: str, anchor: str) -> str | None:
    target = headings(page.run.docs_dir(page.lang) / target_rel) or headings(page.run.docs_dir("en") / target_rel)
    for _, md, gh in target:
        if anchor in (md, gh):
            return gh
    if page.lang != "en":
        for number, md, gh in headings(page.run.docs_dir("en") / target_rel):
            if anchor in (md, gh) and number:
                return next((g for n, _, g in target if n == number), None)
    return None


def apply(page) -> None:
    lines = page.body.split("\n")
    in_fence = fenced_lines(lines)

    def sub(m):
        path, hash_ = m["path"], m["hash"]
        if not path and not hash_:
            return m[0]
        if path and not path.endswith(".md"):
            return m[0]
        target_rel = os.path.normpath(os.path.join(os.path.dirname(page.rel), path)) if path else page.rel
        exists = any((page.run.docs_dir(l) / target_rel).is_file() for l in (page.lang, "en"))
        new_path = path[:-3] + ".mdx" if path and exists else path
        new_hash = hash_
        if hash_ and len(hash_) > 1 and exists:
            gh = resolve_anchor(page, target_rel, hash_[1:])
            if gh is None:
                page.warn(f"anchor {path}{hash_} matches no heading")
            else:
                new_hash = "#" + gh
        return f"{m['pre']}{new_path}{new_hash or ''}"

    page.body = "\n".join(l if i in in_fence else LINK.sub(sub, l) for i, l in enumerate(lines))
