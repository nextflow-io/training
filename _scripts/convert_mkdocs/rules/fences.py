"""A4. Fence lines (language, hl_lines, linenums), code annotations and inline highlight."""

import re

from core import fenced_lines, fences, indent_of

ATTR = re.compile(r'(?P<key>[\w.-]+)="(?P<value>[^"]*)"')
ANNOTATION = re.compile(r"\((?P<n>\d+)\)!")
LIST_ITEM = re.compile(r"^(?P<n>\d+)\.\s+(?P<text>.*)$")
INLINE = re.compile(r"`#!(?P<lang>\w+) (?P<code>[^`]+)`")
GROOVY_TITLE = re.compile(r"\.(groovy|gradle|nf\.test|java)$|^build\.gradle|\.test$", re.I)
GROOVY_CODE = re.compile(r"^\s*(package |import |@\w+|class \w+|(public |protected |private )?(static )?(void|String|def) \w+\()", re.M)
# Pages about the Groovy that sits under Nextflow: untitled snippets there are Groovy.
GROOVY_PAGES = ("side_quests/plugin_development/",)


def nextflow_or_groovy(title: str | None, code: str, rel: str) -> str:
    if title:
        return "groovy" if GROOVY_TITLE.search(title.strip()) else "nextflow"
    if rel.startswith(GROOVY_PAGES) or GROOVY_CODE.search(code):
        return "groovy"
    return "nextflow"


def is_config_title(title: str | None) -> bool:
    return bool(title) and "config" in title.lower() and not title.strip().endswith(".config")


def ranges(spec: str) -> str:
    return ",".join(spec.split())


def fence_line(page, info: str, code: str) -> str:
    lang, _, rest = info.partition(" ")
    attrs = dict((m["key"], m["value"]) for m in ATTR.finditer(rest))
    leftover = ATTR.sub("", rest).strip()
    title = attrs.get("title")
    if lang == "groovy":
        lang = nextflow_or_groovy(title, code, page.rel)
    parts = [lang]
    if "hl_lines" in attrs:
        parts.append("{" + ranges(attrs["hl_lines"]) + "}")
    if title is not None:
        parts.append(f'title="{title}"')
    if lang == "nextflow" and is_config_title(title):
        parts.append('nextflow="config"')
    if "linenums" in attrs:
        parts.append("showLineNumbers")
        if attrs["linenums"] != "1":
            parts.append(f"startLineNumber={attrs['linenums']}")
    for k, v in attrs.items():
        if k not in ("hl_lines", "title", "linenums"):
            parts.append(f'{k}="{v}"')
    if leftover:
        parts.append(leftover)
    return " ".join(parts)


def inline(text: str, rel: str) -> str:
    def sub(m):
        lang = m["lang"]
        if lang == "groovy":
            lang = "groovy" if rel.startswith(GROOVY_PAGES) else "nextflow"
        return f"`{m['code']}{{:{lang}}}`"

    return INLINE.sub(sub, text)


def annotation_list(lines: list[str], start: int, ind: int) -> tuple[dict[int, str], int]:
    """The numbered list after a fence: {n: text} and the index after the list."""
    j = start
    while j < len(lines) and not lines[j].strip():
        j += 1
    notes: dict[int, str] = {}
    while j < len(lines):
        line = lines[j]
        m = LIST_ITEM.match(line[ind:]) if indent_of(line) == ind else None
        if m:
            notes[int(m["n"])] = m["text"].strip()
            j += 1
            continue
        if line.strip() and indent_of(line) > ind and notes:
            notes[max(notes)] += " " + line.strip()
            j += 1
            continue
        if not line.strip() and notes:
            k = j
            while k < len(lines) and not lines[k].strip():
                k += 1
            if k < len(lines) and indent_of(lines[k]) == ind and LIST_ITEM.match(lines[k][ind:]):
                j = k
                continue
        break
    return notes, j


def annotate(text: str) -> str:
    return re.sub(r"`#!\w+ ([^`]*)`", r"`\1`", text)


def apply(page) -> None:
    lines = page.body.split("\n")
    for f in reversed(fences(lines)):
        code_lines = lines[f.start + 1 : f.end]
        code = "\n".join(code_lines)
        if any(ANNOTATION.search(l) for l in code_lines):
            notes, after = annotation_list(lines, f.end + 1, f.indent)
            if notes:
                for k, l in enumerate(code_lines):
                    code_lines[k] = ANNOTATION.sub(lambda m: f"[!annotate] {annotate(notes.get(int(m['n']), ''))}".rstrip(), l)
                lines = lines[: f.start + 1] + code_lines + [lines[f.end]] + [""] + lines[after:]
            else:
                page.warn(f"annotation markers with no list after the block at line {f.start + 1}")
        pad = " " * f.indent
        lines[f.start] = pad + f.marker + fence_line(page, f.info, code)

    in_fence = fenced_lines(lines)
    page.body = "\n".join(l if i in in_fence else inline(l, page.rel) for i, l in enumerate(lines))
