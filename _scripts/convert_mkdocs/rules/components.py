"""A8. Components: quiz, index page marker, grid cards, buttons, icons, keys."""

import re

from core import fenced_lines

BUTTON = re.compile(r"\[(?P<text>[^\]]+)\]\((?P<href>[^)\s]+)\)\{\s*\.md-button(?P<classes>[^}]*)\}")
ICON = re.compile(r":(?:material|fontawesome|octicons|simple)-[a-z0-9-]+:(?:\{[^}]*\})? ?")
EMOJI = {":green_book:": "📗", ":white_check_mark:": "✅", ":warning:": "⚠️", ":tada:": "🎉"}
KEYS = re.compile(r"\+\+(?P<keys>[a-z0-9]+(?:\+[a-z0-9]+)*)\+\+", re.I)
CLEARFIX = re.compile(r"^(<!-- Clearfix for float -->\n)?<div style=\"content: ''; clear: both; display: table;\"></div>\n?", re.M)


def button(m) -> str:
    text = m["text"]
    icon = ""
    if ":material-arrow-right:" in text:
        icon = ' icon="right-arrow"'
    elif ":material-arrow-left:" in text:
        icon = ' icon="left-arrow" iconPlacement="start"'
    text = ICON.sub("", text).strip()
    variant = "primary" if "md-button--primary" in m["classes"] else "secondary"
    return f'<LinkButton href="{m["href"]}" variant="{variant}"{icon}>{text}</LinkButton>'


def keys(m) -> str:
    return "+".join(f"<kbd>{k.capitalize() if len(k) > 1 else k.upper()}</kbd>" for k in m["keys"].split("+"))


def prose(line: str) -> str:
    line = BUTTON.sub(button, line)
    line = ICON.sub("", line)
    for code, char in EMOJI.items():
        line = line.replace(code, char)
    return KEYS.sub(keys, line)


def apply(page) -> None:
    body = CLEARFIX.sub("", page.body)
    if page.meta.get("page_type") == "index_page":
        if "<!-- additional_information -->" not in body:
            page.warn("index_page without the <!-- additional_information --> marker")
        body = body.replace("<!-- additional_information -->", "<CourseInfo frontmatter={frontmatter} />")
    lines = body.split("\n")
    in_fence = fenced_lines(lines)
    quizzes = 0
    for i, line in enumerate(lines):
        if i in in_fence:
            continue
        s = line.strip()
        if s == "<quiz>":
            # Untitled quizzes share one id unless given one; progress is keyed per page.
            quizzes += 1
            lines[i] = line.replace("<quiz>", f'<Quiz id="quiz-{quizzes}">\n')
        elif s == "</quiz>":
            lines[i] = line.replace("</quiz>", "\n</Quiz>")
        else:
            lines[i] = prose(line)
    page.body = "\n".join(lines)
