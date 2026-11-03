"""A8. Components: quiz, index page marker, grid cards, buttons, icons, keys."""

import re

from core import fenced_lines, indent, indent_of

GRID = re.compile(r'^<div class="grid cards[^"]*" markdown>$')
ITEM = re.compile(r"^[-*+]\s+(?P<head>\S.*)$")
RULE = re.compile(r"^(-{3,}|\*{3,}|_{3,})$")
BOLD = re.compile(r"^(\*\*|__)(?P<text>.*)\1$")
BUTTON = re.compile(r"\[(?P<text>[^\]]+)\]\((?P<href>[^)\s]+)\)\{\s*\.md-button(?P<classes>[^}]*)\}")
ICON = re.compile(r":(?P<name>(?:material|fontawesome|octicons|simple)-[a-z0-9-]+):(?:\{[^}]*\})?(?P<space> ?)")
ARROW = re.compile(r":material-arrow-(?P<dir>left|right|up|down):")
CLASSES = re.compile(r"\s*\{\s*\.[\w-]+(?:\s+\.[\w-]+)*\s*\}")
CODE = re.compile(r"(`+).*?(?<!`)\1(?!`)")
KEYS = re.compile(r"\+\+(?P<keys>[a-z0-9]+(?:\+[a-z0-9]+)*)\+\+", re.I)
CLEARFIX = re.compile(r"^(<!-- Clearfix for float -->\n)?<div style=\"content: ''; clear: both; display: table;\"></div>\n?", re.M)
MARKER = "<!-- additional_information -->"

# Every icon and emoji shortcode in the corpus. Icons not listed here are dropped.
ICONS = {
    "material-run-fast": "🏃",
    "material-bookmark": "🔖",
    "material-file-edit-outline": "📝",
    "material-subtitles": "CC",
}
EMOJI = {":green_book:": "📗"}


def button(m) -> str:
    text = m["text"]
    attrs = f'href="{m["href"]}" variant="{"primary" if "md-button--primary" in m["classes"] else "secondary"}"'
    if arrow := ARROW.search(text):
        attrs += f' icon="{arrow["dir"]}-arrow"'
        if not text[: arrow.start()].strip():
            attrs += ' iconPlacement="start"'
    return f"<LinkButton {attrs}>{ICON.sub('', text).strip()}</LinkButton>"


def icon(m) -> str:
    char = ICONS.get(m["name"])
    return char + m["space"] if char else ""


def keys(m) -> str:
    return "+".join(f"<kbd>{k.capitalize() if len(k) > 1 else k.upper()}</kbd>" for k in m["keys"].split("+"))


def prose(text: str) -> str:
    text = ICON.sub(icon, BUTTON.sub(button, text))
    text = CLASSES.sub("", text)
    for code, char in EMOJI.items():
        text = text.replace(code, char)
    return KEYS.sub(keys, text)


def outside_code(line: str) -> str:
    out, i = [], 0
    for m in CODE.finditer(line):
        out += [prose(line[i : m.start()]), m[0]]
        i = m.end()
    return "".join(out) + prose(line[i:])


def card(head: str, body: list[str]) -> list[str]:
    title = BOLD.sub(r"\g<text>", ICON.sub("", head).strip()).replace('"', "&quot;")
    n = min((indent_of(l) for l in body if l.strip()), default=0)
    body = [l[n:] for l in body]
    while body and not body[0].strip():
        body.pop(0)
    if body and RULE.match(body[0].strip()):
        body.pop(0)
    text = "\n".join(body).strip("\n")
    return [f'<Card title="{title}">', "", *([text, ""] if text else []), "</Card>", ""]


def grids(page, lines: list[str]) -> list[str]:
    """`<div class="grid cards" markdown>` with one card per list item becomes <CardGrid> with <Card>s."""
    in_fence = fenced_lines(lines)
    out, i = [], 0
    while i < len(lines):
        if i in in_fence or not GRID.match(lines[i].strip()):
            out.append(lines[i])
            i += 1
            continue
        ind = indent_of(lines[i])
        end = next((j for j in range(i + 1, len(lines)) if lines[j].strip() == "</div>" and indent_of(lines[j]) == ind), None)
        items: list[list[str]] = []
        for l in lines[i + 1 : end]:
            m = ITEM.match(l[ind:]) if l.strip() and indent_of(l) == ind else None
            if m:
                items.append([m["head"]])
            elif items:
                items[-1].append(l[ind:])
        if end is None or not items:
            page.warn("grid cards without a closing </div> or list items")
            out.append(lines[i])
            i += 1
            continue
        grid = ["<CardGrid>", ""]
        for head, *body in items:
            grid += card(head, body)
        out += indent(grid + ["</CardGrid>"], ind)
        i = end + 1
    return out


def block(out: list[str], tag: str, after: str) -> None:
    """Append a JSX tag on its own line, with blank lines to separate it from Markdown."""
    if out and out[-1].strip():
        out.append("")
    out.append(tag)
    if after.strip():
        out.append("")


def apply(page) -> None:
    lines = grids(page, CLEARFIX.sub("", page.body).split("\n"))
    in_fence = fenced_lines(lines)
    index_page = page.meta.get("page_type") == "index_page"
    if index_page and MARKER not in (l.strip() for l in lines):
        page.warn(f"index_page without the {MARKER} marker")
    out: list[str] = []
    quizzes = 0
    for i, line in enumerate(lines):
        after = lines[i + 1] if i + 1 < len(lines) else ""
        pad = line[: indent_of(line)]
        s = line.strip()
        if i in in_fence:
            out.append(line)
        elif s == "<quiz>":
            # Untitled quizzes otherwise share one id, and answers persist per id.
            quizzes += 1
            block(out, f'{pad}<Quiz id="quiz-{quizzes}">', after)
        elif s == "</quiz>":
            block(out, f"{pad}</Quiz>", after)
        elif index_page and s == MARKER:
            block(out, f"{pad}<CourseInfo frontmatter={{frontmatter}} />", after)
        else:
            out.append(outside_code(line))
    page.body = "\n".join(out)
