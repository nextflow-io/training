"""A3. Raw HTML `<details>` blocks become collapsible directives.

    <details class="tip" open>        :::tip[Title]{open}
    <summary>Title</summary>    →     ...
    ...                               :::
    </details>

A `<details>` with no class is a contents box, so it becomes `abstract`.
"""

import re

OPEN = re.compile(r'^(?P<indent>[ \t]*)<details(?: class="(?P<cls>[^"]*)")?(?P<open> open)?>\s*$')
SUMMARY = re.compile(r"^\s*<summary>(?P<title>.*)</summary>\s*$")
CLOSE = re.compile(r"^\s*</details>\s*$")
FENCE = re.compile(r"^\s*(`{3,}|~{3,})")
DIRECTIVE = re.compile(r"^\s*(:{3,})")


def rewrite(text: str, warn=lambda msg: None) -> str:
    lines: list[str | None] = text.split("\n")
    stack: list[int] = []
    fence = ""
    for i, line in enumerate(lines):
        if fence:
            m = FENCE.match(line)
            if m and m[1][0] == fence[0] and len(m[1]) >= len(fence) and not line.strip()[len(m[1]) :]:
                fence = ""
            continue
        if m := FENCE.match(line):
            fence = m[1]
        elif OPEN.match(line) and i + 1 < len(lines) and SUMMARY.match(lines[i + 1] or ""):
            stack.append(i)
        elif CLOSE.match(line) and stack:
            start = stack.pop()
            m = OPEN.match(lines[start])
            classes = (m["cls"] or "abstract").split()
            title = SUMMARY.match(lines[start + 1])["title"].strip()
            if title.count("[") != title.count("]"):
                warn(f"collapsible title {title!r} has unbalanced brackets")
            attrs = " ".join(["open" if m["open"] else "collapsed", *(f".{c}" for c in classes[1:])])
            # An outer directive needs more colons than any directive inside it.
            inner = [len(d[1]) for l in lines[start + 2 : i] if l is not None and (d := DIRECTIVE.match(l))]
            colons = ":" * max([3] + [n + 1 for n in inner])
            lines[start] = f"{m['indent']}{colons}{classes[0]}[{title}]{{{attrs}}}"
            lines[start + 1] = None
            lines[i] = f"{m['indent']}{colons}"
    if stack:
        warn(f"{len(stack)} <details> block(s) never closed")
    return "\n".join(l for l in lines if l is not None)


def apply(page) -> None:
    page.body = rewrite(page.body, page.warn)
