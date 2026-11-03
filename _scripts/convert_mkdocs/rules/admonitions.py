"""A2. `!!!` admonitions become aside directives; `???` collapsibles become <details>."""

import re

from core import apply_containers, containers

ADMONITION = re.compile(r'^!!! (?P<type>[\w-]+)(?: "(?P<title>[^"]*)")?$')
DETAILS = re.compile(r'^\?\?\?(?P<open>\+)? (?P<type>[\w-]+)(?: "(?P<title>[^"]*)")?$')

ASIDE = {
    **dict.fromkeys(["tip", "hint", "success", "check", "done"], "tip"),
    **dict.fromkeys(["warning", "caution", "attention"], "caution"),
    **dict.fromkeys(["danger", "error", "failure", "fail", "missing", "bug"], "danger"),
}


def depth(body: list[str]) -> int:
    """How many levels of admonitions are nested inside this body."""
    return max((1 + depth(c.body) for c in containers(body, ADMONITION)), default=0)


def aside(c) -> list[str]:
    kind = ASIDE.get(c.match["type"], "note")
    title = c.match["title"]
    colons = ":" * (3 + depth(c.body))
    label = f"[{title}]" if title else ""
    return [f"{colons}{kind}{label}", *c.body, colons]


def details(c) -> list[str]:
    kind = c.match["type"]
    title = c.match["title"]
    if title is None:
        title = kind.replace("-", " ").capitalize()
    opener = f'<details class="{kind}"' + (" open" if c.match["open"] else "") + ">"
    summary = f"<summary>{title}</summary>" if title else ""
    return [opener, *([summary] if summary else []), "", *c.body, "", "</details>"]


def apply(page) -> None:
    body = apply_containers(page.body, DETAILS, details)
    page.body = apply_containers(body, ADMONITION, aside)
