"""A2. `!!!` admonitions become `:::type[Title]` directives; `???` collapsibles add `{collapsed}` (`???+`: `{open}`).

Headers are parsed as Python-Markdown and pymdownx do: optional space after the
marker, lower-cased classes (type first, then modifiers such as `inline end`),
and a title that runs to the last quote on the line.
"""

import re

from core import apply_containers, fences

HEADER = re.compile(
    r'^(?P<mark>!!!|\?\?\?\+?) ?(?P<classes>[\w-]+(?: +[\w-]+)*?)(?: +"(?P<title>.*?)")?'
    r'(?P<attrs>(?: +[\w-]+="[^"]*")*) *$'
)
ATTR = re.compile(r'([\w-]+)="([^"]*)"')
DIRECTIVE = re.compile(r"^\s*(:{3,})")

# Types docs/src/markdown/admonitions.mjs renders; keep the two lists in step. Others become notes.
TYPES = {
    "note", "abstract", "info", "tip", "success", "question", "warning", "caution", "failure", "danger", "bug",
    "example", "quote", "exercise", "solution", "full-code", "optional", "courses", "catalog",
}
# An untitled aside (`!!! tip ""`) keeps Starlight's aside markup in one of its four variants.
VARIANT = {"note": "note", "info": "note", "tip": "tip", "warning": "caution", "danger": "danger"}


def convert(c, page) -> list[str]:
    body = apply_containers("\n".join(c.body), HEADER, lambda d: convert(d, page)).split("\n")
    classes = c.match["classes"].lower().split()
    kind, title = classes[0], c.match["title"]
    if title and title.count("[") != title.count("]"):
        page.warn(f"admonition title {title!r} has unbalanced brackets")
    # Attributes on a header are authoring slips: `hl_lines` was meant for the code block inside.
    for key, value in ATTR.findall(c.match["attrs"]):
        f = fences(body)[0] if key == "hl_lines" and fences(body) else None
        if f:
            body[f.start] += " {" + ",".join(value.split()) + "}"
        page.warn(f"{key}=\"{value}\" on admonition header " + ("moved to its code block" if f else "dropped"))

    collapsible = c.match["mark"].startswith("?")
    if title == "" and not collapsible:
        # A directive always renders a title, so use Starlight's aside markup without one.
        cls = " ".join(["starlight-aside", f"starlight-aside--{VARIANT.get(kind, 'note')}", *classes[1:]])
        return [f'<div class="{cls}">', "", *body, "", "</div>"]
    if kind not in TYPES:
        kind, title = "note", title or kind.capitalize()
    attrs = []
    if collapsible:
        attrs = ["open" if c.match["mark"].endswith("+") else "collapsed", *(f".{x}" for x in classes[1:])]
    # An outer directive needs more colons than any directive inside it.
    colons = ":" * max([3] + [len(m[1]) + 1 for l in body if (m := DIRECTIVE.match(l))])
    opener = f"{colons}{kind}" + (f"[{title}]" if title else "") + (f"{{{' '.join(attrs)}}}" if attrs else "")
    return [opener, *body, colons]


def apply(page) -> None:
    page.body = apply_containers(page.body, HEADER, lambda c: convert(c, page))
