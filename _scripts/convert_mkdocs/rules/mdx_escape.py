"""A6. Make prose safe for MDX.

Runs last, so it sees the JSX that earlier rules wrote and must leave it alone:
known tags and `{/* */}` comments pass through; `{`, `}` and a `<` that does not
open a known tag are escaped. Fenced code and inline code are never touched.
"""

import re
from core import fenced_lines

HTML_TAGS = set(
    "a abbr b br button code del details div em figcaption figure hr i iframe img input ins kbd label li mark "
    "ol p path picture pre s script section small source span strong sub summary sup svg table tbody td th thead tr u ul video".split()
)
VOID_TAGS = {"area", "br", "col", "embed", "hr", "img", "input", "source", "track", "wbr"}
TOKEN = re.compile(
    r"(?P<esc>\\[!-/:-@\[-`{-~])"
    r"|(?P<code>`+)"
    r"|(?P<dest>\]\(<[^<>\n]*>)"
    r"|(?P<comment>\{/\*.*?\*/\})"
    r"|<!--(?P<html_comment>.*?)-->"
    r"|(?P<autolink><(?P<url>(?:https?|ftp)://[^>\s]+|mailto:[^>\s]+)>)"
    r"|(?P<tag></?(?P<name>[A-Za-z][\w.:-]*)(?P<attrs>(?:\s(?:[^<>\"'{}]|\"[^\"]*\"|'[^']*'|\{[^{}]*\})*?)?)(?P<self>/?)>)"
    r"|(?P<special>[{}<])"
)
CODE_SPAN = re.compile(r"(`+).*?(?<!`)\1(?!`)")
HEADING_ID = re.compile(r"^(?P<head>\s*#{1,6}\s.*?)(?P<id>\s\{#[\w-]+\})\s*$")
ESM_START = re.compile(r"^(\s*)(import|export)(\s)")
ATTR_VALUE = re.compile(r"\"[^\"]*\"|'[^']*'|\{[^{}]*\}|(?P<name>\s[\w:-]+=)(?P<bare>[^\s\"'{>/][^\s>]*)")


def is_tag(name: str) -> bool:
    return name.lower() in HTML_TAGS or name[0].isupper()


def tag(m: re.Match) -> str:
    """A known tag as valid JSX: attribute values quoted, void elements self-closed."""
    name, attrs = m["name"], ATTR_VALUE.sub(lambda v: f'{v["name"]}"{v["bare"]}"' if v["bare"] else v[0], m["attrs"])
    if m[0].startswith("</"):
        return f"</{name}{attrs}>"
    if m["self"] or name.lower() in VOID_TAGS:
        return f"<{name}{attrs.rstrip()} />"
    return f"<{name}{attrs}>"


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
        if m["html_comment"] is not None:
            out.append("{/* " + m["html_comment"].strip().replace("*/", "* /") + " */}")
        elif m["autolink"]:
            url = m["url"]
            out.append(f"[{url}](<{url}>)" if re.search(r"[()]", url) else f"[{url}]({url})")
        elif m["tag"]:
            out.append(tag(m) if is_tag(m["name"]) else "&lt;" + escape_line(m["tag"][1:]))
        elif m["special"] == "<":
            out.append("&lt;")
        elif m["special"]:
            out.append("\\" + m["special"])
        else:
            out.append(m[0])
        i = m.end()


def join_comments(lines: list[str], in_fence: set[int]) -> list[str]:
    """An HTML comment over several lines is joined onto its first line."""
    out: list[str] = []
    i = 0
    while i < len(lines):
        line = lines[i]
        prose = CODE_SPAN.sub("", line)
        if i not in in_fence and "<!--" in prose:
            while "-->" not in prose[prose.rindex("<!--") :] and i + 1 < len(lines):
                i += 1
                line += " " + lines[i].strip()
                prose = CODE_SPAN.sub("", line)
        out.append(line)
        i += 1
    return out


def apply(page) -> None:
    lines = page.body.split("\n")
    lines = join_comments(lines, fenced_lines(lines))
    in_fence = fenced_lines(lines)
    for i, line in enumerate(lines):
        if i in in_fence:
            continue
        # A trailing `{#id}` on a heading is a heading attribute (written by the links rule), not an expression.
        m = HEADING_ID.match(line)
        line = escape_line(m["head"]) + m["id"] if m else escape_line(line)
        lines[i] = ESM_START.sub(lambda m: f"{m[1]}&#{ord(m[2][0])};{m[2][1:]}{m[3]}", line)
    page.body = "\n".join(lines)
