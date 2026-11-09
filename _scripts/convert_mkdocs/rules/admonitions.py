"""A2. `!!!` admonitions become aside directives; `???` collapsibles become <details>.

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

# An untitled admonition of these types takes Starlight's (translated) default title.
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

    if c.match["mark"].startswith("?"):
        opener = f'<details class="{" ".join(classes)}"' + (" open" if c.match["mark"].endswith("+") else "") + ">"
        summary = kind.capitalize() if title is None else title
        return [opener, f"<summary>{summary}</summary>", "", *body, "", "</details>"]

    variant = VARIANT.get(kind, "note")
    if title == "":
        # A directive always renders a title, so use Starlight's aside markup without one.
        cls = " ".join(["starlight-aside", f"starlight-aside--{variant}", *classes[1:]])
        return [f'<div class="{cls}">', "", *body, "", "</div>"]
    if title is None and kind not in VARIANT:
        title = kind.capitalize()
    # An outer directive needs more colons than any directive inside it.
    colons = ":" * max([3] + [len(m[1]) + 1 for l in body if (m := DIRECTIVE.match(l))])
    return [f"{colons}{variant}" + (f"[{title}]" if title else ""), *body, colons]


def apply(page) -> None:
    page.body = apply_containers(page.body, HEADER, lambda c: convert(c, page))
