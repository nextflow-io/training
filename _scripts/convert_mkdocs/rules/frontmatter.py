"""A1. Title from the H1, `hide`, and a generated `description`."""

import re

from core import fenced_lines

H1 = re.compile(r"^# (?P<text>.+?)\s*(\{[^}]*\})?\s*$")
SENTENCE_END = re.compile(r"(?<=[.!?。।])\s+")
SKIP_BLOCK = re.compile(r"^\s*(<|/// |!!! |\?\?\?|=== |#|```|~~~|\||- |\* |\d+\. |!\[|---|\[.*\]\(.*\)\{)")


def plain(text: str) -> str:
    text = re.sub(r"\{[^}]*\}", "", text)  # attr lists
    text = re.sub(r":[a-z0-9_+-]+(?:-[a-z0-9_-]+)*:", "", text)  # icon and emoji shortcodes
    text = re.sub(r"!?\[([^\]]*)\]\([^)]*\)", r"\1", text)  # links and images
    text = re.sub(r"`#!\w+ ([^`]*)`", r"\1", text)
    text = re.sub(r"(\*\*|__|`|\*|(?<!\w)_(?=\S)|(?<=\S)_(?!\w))", "", text)
    text = re.sub(r"<[^>]+>", "", text)
    return re.sub(r"\s+", " ", text).strip()


def first_sentence(lines: list[str]) -> str | None:
    in_fence = fenced_lines(lines)
    para: list[str] = []
    container_indent = None
    in_caption = False
    for i, line in enumerate(lines):
        if i in in_fence:
            continue
        if line.startswith("/// "):
            in_caption = True
            continue
        if in_caption:
            in_caption = line.strip() != "///"
            continue
        if container_indent is not None:
            # skip the indented body of an admonition, tab or caption
            if not line.strip() or line.startswith(" " * 4):
                continue
            container_indent = None
        if not line.strip():
            if para:
                break
            continue
        if not para and SKIP_BLOCK.match(line):
            if re.match(r"^\s*(!!! |\?\?\?|=== )", line):
                container_indent = 4
            continue
        if not para and line.startswith("    "):
            continue
        para.append(line.strip())
    if not para:
        return None
    text = plain(" ".join(para))
    return SENTENCE_END.split(text, 1)[0] if text else None


def apply(page) -> None:
    lines = page.body.split("\n")
    in_fence = fenced_lines(lines)
    title = None
    for i, line in enumerate(lines):
        if i in in_fence or not line.strip():
            continue
        m = H1.match(line)
        if m:
            title = plain(m["text"])
            del lines[i]
            if i < len(lines) and not lines[i].strip():
                del lines[i]
        break
    if title is None:
        page.warn("no H1 on the first line")
        title = page.meta.get("title") or page.rel
    page.body = "\n".join(lines)

    meta = dict(page.meta)
    old_title = meta.pop("title", None)
    out = {"title": title}
    description = meta.pop("description", None) or first_sentence(lines)
    if description:
        out["description"] = description
    if old_title and old_title != title:
        out["sidebar"] = {"label": old_title}
    hide = meta.pop("hide", None) or []
    if "toc" in hide:
        out["tableOfContents"] = False
    for h in hide:
        if h != "toc":
            page.warn(f"hide: {h} has no Starlight equivalent")
    out.update(meta)
    page.meta = out
