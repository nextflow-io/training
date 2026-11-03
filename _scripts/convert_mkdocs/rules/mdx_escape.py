"""A6. Make prose safe for MDX.

Runs last, so it sees the JSX that earlier rules wrote and must leave it alone:
known tags and `{/* */}` comments pass through; `{`, `}` and a `<` that does not
open a known tag are escaped. Fenced code and inline code are never touched.
"""

import re
from core import fenced_lines

HTML_TAGS = set(
    "a abbr b br button code del details div em figcaption figure hr i iframe img input ins kbd label li mark "
    "ol p path picture pre s section small source span strong sub summary sup svg table tbody td th thead tr u ul video".split()
)
TOKEN = re.compile(
    r"(?P<code>`+)"
    r"|(?P<comment>\{/\*.*?\*/\})"
    r"|(?P<autolink><(?P<url>https?://[^>\s]+)>)"
    r"|(?P<tag></?(?P<name>[A-Za-z][\w.:-]*)(?:\s(?:[^<>\"']|\"[^\"]*\"|'[^']*')*?)?/?>)"
    r"|(?P<special>[{}<])"
)
ESM_START = re.compile(r"^(\s*)(import|export)(\s)")


def is_tag(name: str) -> bool:
    return name.lower() in HTML_TAGS or name[0].isupper()


def escape_line(line: str) -> str:
    out, i = [], 0
    while True:
        m = TOKEN.search(line, i)
        if not m:
            out.append(line[i:])
            return "".join(out)
        out.append(line[i : m.start()])
        if m["code"]:
            run = m["code"]
            close = re.compile(rf"(?<!`){run}(?!`)").search(line, m.end())
            end = close.end() if close else m.end()
            out.append(line[m.start() : end])
            i = end
            continue
        if m["autolink"]:
            out.append(f"[{m['url']}]({m['url']})")
        elif m["tag"]:
            if is_tag(m["name"]):
                out.append(m["tag"])
            else:
                out.append("&lt;" + escape_line(m["tag"][1:]))
        elif m["special"] == "<":
            out.append("&lt;")
        elif m["special"]:
            out.append("\\" + m["special"])
        else:
            out.append(m[0])
        i = m.end()


def comments(lines: list[str], in_fence: set[int]) -> list[str]:
    """HTML comments become MDX comments. A comment over several lines is joined onto one."""
    out: list[str] = []
    i = 0
    while i < len(lines):
        line = lines[i]
        if i in in_fence or "<!--" not in line:
            out.append(line)
            i += 1
            continue
        while "-->" not in line[line.index("<!--") :] and i + 1 < len(lines):
            i += 1
            line += " " + lines[i].strip()
        out.append(re.sub(r"<!--\s*(.*?)\s*-->", lambda m: "{/* " + m[1].replace("*/", "* /") + " */}", line))
        i += 1
    return out



def apply(page) -> None:
    lines = page.body.split("\n")
    lines = comments(lines, fenced_lines(lines))
    in_fence = fenced_lines(lines)
    for i, line in enumerate(lines):
        if i in in_fence:
            continue
        line = escape_line(line)
        lines[i] = ESM_START.sub(lambda m: f"{m[1]}&#{ord(m[2][0])};{m[2][1:]}{m[3]}", line)
    page.body = "\n".join(lines)
