"""A5. A shell command block followed by a success/failure output box becomes one
block with scripted output: the command, then `# [!output]` and the output lines.

The box is matched by structure, not by its title, because the title is translated.
"""

import re

from core import containers, fences, indent_of

SHELL = {"bash", "sh", "shell", "console", "zsh"}
OUTPUT_BOX = re.compile(r'^\?\?\?\+? (?P<type>success|failure) ".*"$')
HL_LINES = re.compile(r'\s*hl_lines="(?P<v>[^"]*)"')
LONG = 60


def line_set(spec: str) -> set[int]:
    out: set[int] = set()
    for part in spec.split():
        a, _, b = part.partition("-")
        out.update(range(int(a), int(b or a) + 1))
    return out


def apply(page) -> None:
    lines = page.body.split("\n")
    # Work bottom-up so earlier indexes stay valid.
    for cmd in reversed(fences(lines)):
        if cmd.lang not in SHELL:
            continue
        j = cmd.end + 1
        while j < len(lines) and not lines[j].strip():
            j += 1
        if j >= len(lines) or indent_of(lines[j]) != cmd.indent or not OUTPUT_BOX.match(lines[j].strip()):
            continue
        box = containers(lines[j:], OUTPUT_BOX)[0]
        nonblank = [i for i, l in enumerate(box.body) if l.strip()]
        out_fences = fences(box.body)
        if len(out_fences) != 1 or (nonblank[0], nonblank[-1]) != (out_fences[0].start, out_fences[0].end):
            page.warn(f"output box at line {j + 1} holds more than one code block; left as <details>")
            continue
        of = out_fences[0]
        output = [l[of.indent :] for l in box.body[of.start + 1 : of.end]]
        while output and not output[-1].strip():
            output.pop()
        while output and not output[0].strip():
            output.pop(0)

        hl = HL_LINES.search(of.info)
        marked = line_set(hl["v"]) if hl else set()
        failure = box.match["type"] == "failure"
        new_out: list[str] = []
        stderr = False
        for n, l in enumerate(output, 1):
            if failure and not stderr and l.lstrip().startswith("ERROR ~"):
                new_out.append("# [!output stderr]")
                stderr = True
            new_out.append(f"{l} # [!code highlight]" if n in marked else l)
        if not output:
            block = []  # the command prints nothing: drop the empty box
        else:
            block = ["# [!output]", *new_out]

        pad = " " * cmd.indent
        opener = lines[cmd.start]
        if len(output) > LONG:
            opener += " runnable.outputDelay=0"
        new = [opener, *lines[cmd.start + 1 : cmd.end], *[pad + l if l else l for l in block], lines[cmd.end]]
        lines = lines[: cmd.start] + new + lines[j + box.end + 1 :]
    page.body = "\n".join(lines)
