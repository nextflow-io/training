"""A5. A shell command block followed by an output box (`??? success|failure|warning "…"`)
becomes one block with scripted output: the command, then `# [!output]` and the output lines.

The box is matched by structure, not by its title, because the title is translated.
Boxes that do not fit stay as they are and become <details> in the admonitions rule:
a box not right after a shell block, a box with no or several code blocks, and `???+`,
which asks for the output to show without a click.
"""

import re

from core import containers, dedent, fences, indent, indent_of

# Shell languages whose comments start with `#` in starlight-codeblocks.
SHELL = {"bash", "sh", "shell", "console", "zsh"}
OUTPUT_BOX = re.compile(r'^\?\?\? (?P<type>success|failure|warning) ".*"$')
HL_LINES = re.compile(r'hl_lines="(?P<v>[^"]*)"')
WHOLE_COMMENT = re.compile(r"^[+-]?\s*#")
LONG = 60


def line_set(spec: str) -> set[int]:
    out: set[int] = set()
    for part in spec.split():
        a, _, b = part.partition("-")
        out.update(range(int(a), int(b or a) + 1))
    return out


def escape(line: str) -> str:
    """`[!` after a `#` would be read as a directive; `[\\!` shows as `[!` there."""
    h = line.find("#")
    return line if h == -1 else line[:h] + line[h:].replace("[!", "[\\!")


def scripted(raw: list[str], hl: set[int], stderr: bool) -> list[str]:
    """`# [!output]` and the output lines, with `hl` (1-based) highlighted."""
    # MkDocs (Pygments stripnl) drops leading and trailing blank lines before it counts hl_lines.
    keep = [n for n, l in enumerate(raw) if l.strip()]
    if not keep:
        return []
    raw = raw[keep[0] : keep[-1] + 1]
    err = next((n for n, l in enumerate(raw) if l.lstrip().startswith("ERROR ~")), None) if stderr else None
    out = [] if err == 0 else ["# [!output]"]
    for n, line in enumerate(raw):
        if n == err:
            while out and not out[-1].strip():
                out.pop()
            out.append("# [!output stderr]")
        line = escape(line)
        # A blank line holding only a directive is dropped and passes the directive to the next line.
        if n + 1 in hl and line.strip():
            # A line that starts with `#` is all comment to the plugin: a second `#` would show.
            line += " [!code highlight]" if WHOLE_COMMENT.match(line) else " # [!code highlight]"
        out.append(line)
    return out


def apply(page) -> None:
    lines = page.body.split("\n")
    for cmd in reversed(fences(lines)):
        if cmd.lang not in SHELL:
            continue
        j = cmd.end + 1
        while j < len(lines) and not lines[j].strip():
            j += 1
        if j >= len(lines) or indent_of(lines[j]) != cmd.indent or not OUTPUT_BOX.match(lines[j].strip()):
            continue
        box = containers(lines[j:], OUTPUT_BOX)[0]
        found = fences(box.body)
        if len(found) != 1:
            page.warn(f"output box at line {j + 1} holds {len(found)} code blocks; left as <details>")
            continue
        of = found[0]
        raw = dedent(box.body[of.start + 1 : of.end], of.indent)
        hl = HL_LINES.search(of.info)
        out = scripted(raw, line_set(hl["v"]) if hl else set(), box.match["type"] != "success")
        # Prose in the box comments on the output, so it follows the block.
        prose: list[str] = []
        for part in (box.body[: of.start], box.body[of.end + 1 :]):
            text = "\n".join(part).strip()
            if text:
                prose += ([""] if prose else []) + text.split("\n")

        marker = of.marker if len(of.marker) > len(cmd.marker) else cmd.marker
        pad = " " * cmd.indent
        opener = lines[cmd.start].replace(cmd.marker, marker, 1)
        if sum(not l.startswith("# [!output") for l in out) > LONG:
            opener += " runnable.outputDelay=0"
        new = [opener, *lines[cmd.start + 1 : cmd.end], *indent(out, cmd.indent), pad + marker]
        if prose:
            new += ["", *indent(prose, cmd.indent)]
        lines = lines[: cmd.start] + new + lines[j + box.end + 1 :]
    page.body = "\n".join(lines)
