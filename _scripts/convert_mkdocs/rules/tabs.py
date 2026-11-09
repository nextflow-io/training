"""A3. Content tabs.

- An After/Before pair, each tab one code block, becomes a <CodeWalkthrough>:
  Before first, then After, labels kept from the source (so translations keep theirs).
- Other tab groups where every tab is one code block become `:::code-tabs`.
- Anything else becomes Starlight <Tabs>. Tabs nested in a tab are converted in a later pass.
"""

import re

from core import containers, fences, indent, replace_lines

TAB = re.compile(r'^===\+? "(?P<label>[^"]*)"$')

# Tab labels in every locale of the corpus. Any locale accepts any of them.
AFTER = {"After", "Després", "Danach", "Después", "Après", "बाद में", "Dopo", "후", "Po", "Depois", "Sonra"}
BEFORE = {"Before", "Abans", "Vorher", "Antes", "Avant", "पहले", "Prima", "전", "Przed", "Önce"}


def groups(lines: list[str]):
    """Runs of tab containers at one indent, separated only by blank lines."""
    found = containers(lines, TAB)
    group: list = []
    for c in found:
        if group:
            prev = group[-1]
            gap = lines[prev.end + 1 : c.start]
            if c.indent == prev.indent and all(not l.strip() for l in gap):
                group.append(c)
                continue
            yield group
        group = [c]
    if group:
        yield group


def single_fence(body: list[str]):
    """The one fence that is the whole tab body, or None."""
    fs = fences(body)
    nonblank = [i for i, l in enumerate(body) if l.strip()]
    if len(fs) == 1 and nonblank and nonblank[0] == fs[0].start and nonblank[-1] == fs[0].end:
        return fs[0]
    return None


def fence_with(body: list[str], f, attr: str) -> list[str]:
    block = body[f.start : f.end + 1]
    block[0] = f"{block[0]} {attr}"
    return block


def pair(labels: list[str]) -> tuple[int, int] | None:
    """Indexes (before, after) if the labels are one After and one Before, in either order."""
    if len(labels) != 2:
        return None
    if labels[0] in AFTER and labels[1] in BEFORE:
        return 1, 0
    if labels[0] in BEFORE and labels[1] in AFTER:
        return 0, 1
    return None


def convert(group, page) -> list[str]:
    labels = [c.match["label"] for c in group]
    single = [single_fence(c.body) for c in group]
    steps = pair(labels)
    # code-tabs needs two blocks or more; a lone tab falls through to <Tabs>.
    if all(single) and len(group) > 1:
        if steps:
            out = ["<CodeWalkthrough>", ""]
            for i in steps:
                out += [*fence_with(group[i].body, single[i], f'step="{labels[i]}"'), ""]
            return out + ["</CodeWalkthrough>"]
        if set(labels) & (AFTER | BEFORE):
            page.warn(f"tab group {labels} looks like After/Before but is not a pair")
        out = [":::code-tabs"]
        for c, f, label in zip(group, single, labels):
            out += fence_with(c.body, f, f'label="{label}"')
        return out + [":::"]
    if steps:
        page.warn(f"After/Before tabs {labels} hold prose or several code blocks: kept as <Tabs>")
    out = ["<Tabs>", ""]
    for c, label in zip(group, labels):
        out += [f'<TabItem label="{label}">', "", *c.body, "", "</TabItem>", ""]
    return out + ["</Tabs>"]


def apply(page) -> None:
    lines = page.body.split("\n")
    while True:
        gs = list(groups(lines))
        if not gs:
            break
        for g in reversed(gs):
            lines = replace_lines(lines, g[0].start, g[-1].end, indent(convert(g, page), g[0].indent))
    page.body = "\n".join(lines)
