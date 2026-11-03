"""A1. Title from the H1, `hide`, and a `description` for every page."""

import html
import re
from pathlib import Path

import yaml

from core import fenced_lines, split_frontmatter

HEADING = re.compile(r"^(?P<level>#{1,6}) (?P<text>.+?)(?:\s+#+)?\s*$")
# Translation pipeline output that leaked into the body: LLM chatter, then a second frontmatter block, then the H1.
STRAY_FRONTMATTER = re.compile(r"^---\n(?P<yaml>.*?)\n---\n(?=\s*^# )", re.S | re.M)
SENTENCE = re.compile(r"(?<=[.!?])\s+(?=\S)|(?<=[।。])\s*")
# A period that does not end a sentence: "e.g.", "z. B.", Turkish and German ordinals ("1.", "1-4.").
NOT_AN_END = re.compile(r"(?:^|\W)(?:\w|\d+(?:-\d+)?|etc|vs|ex|ej|np|bzw|usw|örn|vb|ecc|cf)\.$", re.I)
SKIP_BLOCK = re.compile(
    r"^(<|/// |!!!|\?\?\?|=== |#|\||[-*+] |\d+[.)] |!\[|\[!\[|\[\^|\*\[|>|-{3,}\s*$|\*{3,}\s*$|_{3,}\s*$|--8<--|\[[^\]]*\]\([^)]*\)\{)"
)
CODE = re.compile(r"(`+)(?P<code>.+?)\1")
DESCRIPTION_MAX = 160


def plain(text: str) -> str:
    """Markdown and HTML inline markup to plain text."""
    stashed: list[str] = []

    def stash(literal: str) -> str:
        stashed.append(literal)
        return f"\0{len(stashed) - 1}\0"

    text = CODE.sub(lambda m: stash(re.sub(r"^#!\w+ ", "", m["code"].strip())), text)
    text = re.sub(r"\\([\\`*_{}\[\]()#+\-.!<>|])", lambda m: stash(m[1]), text)
    text = re.sub(r"\{:?\s*(?:[.#]|\w[\w-]*=)[^}]*\}", "", text)  # attr lists
    text = re.sub(r"!\[[^\]]*\](?:\([^)]*\)|\[[^\]]*\])", "", text)
    text = re.sub(r"\[\^[^\]]+\]", "", text)
    text = re.sub(r"\[([^\]]*)\](?:\([^)]*\)|\[[^\]]*\])", r"\1", text)
    text = re.sub(r"<((?:https?|mailto):[^>\s]+)>", r"\1", text)
    text = re.sub(r"</?[A-Za-z][^>]*>", "", text)
    text = re.sub(r"(?<![\w:]):[a-z][a-z0-9_+-]*:(?![\w:])", "", text)  # icon and emoji shortcodes
    text = re.sub(r"(\*\*|__|~~|==|\^\^)(?=\S)(.+?)(?<=\S)\1", r"\2", text)
    text = re.sub(r"\*(?=\S)(.+?)(?<=\S)\*", r"\1", text)
    text = re.sub(r"(?<!\w)_(?=\S)(.+?)(?<=\S)_(?!\w)", r"\1", text)
    text = re.sub(r"\0(\d+)\0", lambda m: stashed[int(m[1])], text)
    return re.sub(r"\s+", " ", html.unescape(text)).strip()


def intro_paragraph(lines: list[str]) -> str | None:
    """The first paragraph between the H1 and the first H2, or None."""
    fenced = fenced_lines(lines)
    para: list[str] = []
    skipping, closer = False, None
    for i, line in enumerate(lines):
        if closer:
            if closer(line):
                closer = None
            continue
        if i in fenced or not line.strip():
            if para:
                break
            skipping = False
            continue
        if para:
            if SKIP_BLOCK.match(line):
                break
            para.append(line.strip())
            continue
        if line.startswith("## "):
            return None
        if skipping:
            continue
        if line.startswith("/// "):
            closer = lambda l: l.strip() == "///"  # noqa: E731
        elif line.startswith("<!--") and "-->" not in line:
            closer = lambda l: "-->" in l  # noqa: E731
        elif line[0].isspace() or SKIP_BLOCK.match(line):
            skipping = True
        else:
            para.append(line.strip())
    return " ".join(para) if para else None


def first_sentence(text: str) -> str:
    parts = SENTENCE.split(text)
    sentence = parts[0]
    for part in parts[1:]:
        if not NOT_AN_END.search(sentence.rsplit(" ", 1)[-1]):
            break
        sentence += " " + part
    return sentence


def summarise(text: str) -> str:
    sentence = first_sentence(plain(text))
    if len(sentence) <= DESCRIPTION_MAX:
        return sentence
    return sentence[: DESCRIPTION_MAX - 1].rsplit(" ", 1)[0].rstrip(" ,;:-–—(") + "…"


def h1_of(path: Path) -> str | None:
    if not path.is_file():
        return None
    _, body = split_frontmatter(path.read_text(encoding="utf-8"))
    m = next((HEADING.match(l) for l in body.split("\n") if l.startswith("# ")), None)
    return plain(m["text"]) if m else None


def fallback_description(page, title: str) -> str:
    """`Course: Page title`, from the nearest index.md above the page, in the page's language."""
    root = page.run.docs_dir(page.lang)
    parent = Path(page.rel).parent
    if Path(page.rel).name == "index.md":
        parent = parent.parent
    while parent != Path("."):
        course = h1_of(root / parent / "index.md")
        if course:
            return f"{course}: {title}"
        parent = parent.parent
    return title


def english_description(page) -> str | None:
    path = page.run.docs_dir("en") / page.rel
    return split_frontmatter(path.read_text(encoding="utf-8"))[0].get("description") if path.is_file() else None


def leaked_frontmatter(body: str) -> tuple[dict, str] | None:
    if next((l for l in body.split("\n") if l.strip()), "").startswith("# "):
        return None
    m = STRAY_FRONTMATTER.search(body)
    try:
        meta = yaml.safe_load(m["yaml"]) if m else None
    except yaml.YAMLError:
        return None
    return (meta, body[m.end() :]) if isinstance(meta, dict) else None


def apply(page) -> None:
    meta, body = dict(page.meta), page.body
    if leaked := leaked_frontmatter(body):
        page.warn("dropped text and a second frontmatter block before the H1 (translation pipeline output)")
        meta.update(leaked[0])
        body = leaked[1]

    lines = body.split("\n")
    fenced = fenced_lines(lines)
    title = None
    for i, line in enumerate(lines):
        m = HEADING.match(line) if i not in fenced else None
        if not m:
            continue
        if m["level"] == "#":
            title = plain(m["text"])
            del lines[i]
            if i < len(lines) and not lines[i].strip():
                del lines[i]
        break
    old_title = meta.pop("title", None)
    if title is None:
        page.warn("no H1")
        title = old_title or Path(page.rel).stem.replace("_", " ").capitalize()
    page.body = "\n".join(lines)

    out = {"title": title}
    description = str(meta.pop("description", "") or "").strip()
    if description and page.lang != "en" and description == english_description(page):
        page.warn("description: is untranslated, generated one from the page instead")
        description = ""
    if not description:
        intro = intro_paragraph(lines)
        description = (summarise(intro) if intro else "") or fallback_description(page, title)
    out["description"] = description
    if old_title and old_title != title:
        out["sidebar"] = {"label": old_title}
    hide = meta.pop("hide", None) or []
    if "toc" in hide:
        out["tableOfContents"] = False
    if "footer" in hide:
        out["prev"] = False
        out["next"] = False
    for h in hide:
        if h not in ("toc", "footer"):
            page.warn(f"hide: {h} has no Starlight equivalent")
    out.update(meta)
    page.meta = out
