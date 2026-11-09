"""A4. Fence lines (language, hl_lines, linenums), code annotations and inline highlight.

`groovy` blocks are Nextflow unless they are real Groovy (plugin sources, Gradle, nf-test).
A file name in the title decides; titles are translated, so blocks without one go by content.
"""

import re

from core import fenced_lines, fences, indent_of

ATTR = re.compile(r'(?P<key>[\w.-]+)="(?P<value>[^"]*)"')
ANNOTATION = re.compile(r"\((?P<n>\d+)\)!")
LIST_ITEM = re.compile(r"^(?P<n>\d+)\.\s+(?P<text>.*)$")
INLINE = re.compile(r"`#!(?P<lang>\w+) (?P<code>[^`]+)`")
TITLE_EXT = re.compile(r"\.(nf|config|groovy|gradle|test|java)\b(?!\.)")
GROOVY_EXT = {"groovy", "gradle", "test", "java"}
GROOVY_CODE = re.compile(
    r"^\s*(package [\w.]+|import [\w.]+|@[A-Z]\w*|(public|protected|private|static|final|abstract) |class \w+"
    r"|def '|nextflow_(process|workflow|pipeline|function) \{|nextflowPlugin \{|dependencies \{)",
    re.M,
)
SCOPES = (
    "process|profiles|plugins|docker|singularity|apptainer|podman|charliecloud|conda|spack|executor"
    "|manifest|env|aws|azure|google|k8s|wave|fusion|tower|timeline|report|trace|dag|notification|mail|lineage"
)
# Fence languages in the corpus that Shiki doesn't know: MkDocs used the file name as the language.
LANGUAGES = {"Dockerfile": "dockerfile", "conda.yml": "yaml"}
# Not `params.x = ...`: old-style scripts set params that way too.
CONFIG_START = re.compile(rf"^((withName|withLabel)\s*:|(params|{SCOPES})\s*\{{|({SCOPES})\.\w|ext\.\w|(outputDir|workDir|includeConfig)\b)")


def classify(title: str | None, code: str) -> tuple[str, bool]:
    """(language, is a config block) for a `groovy` block."""
    m = TITLE_EXT.search(title or "")
    if m:
        return ("groovy", False) if m[1] in GROOVY_EXT else ("nextflow", m[1] == "config")
    if GROOVY_CODE.search(code):
        return "groovy", False
    first = next((s for s in (l.strip() for l in code.split("\n")) if s and not s.startswith(("//", "/*", "*"))), "")
    return "nextflow", bool(CONFIG_START.match(first))


def fence_line(info: str, code: str) -> str:
    lang, _, rest = info.partition(" ")
    lang = LANGUAGES.get(lang, lang)
    attrs = dict((m["key"], m["value"]) for m in ATTR.finditer(rest))
    leftover = ATTR.sub("", rest).strip()
    title = attrs.get("title")
    config = False
    if lang == "groovy":
        lang, config = classify(title, code)
    parts = [lang]
    if "hl_lines" in attrs:
        parts.append("{" + ",".join(attrs["hl_lines"].split()) + "}")
    if title is not None:
        parts.append(f'title="{title}"')
    if config and not (title or "").strip().endswith(".config"):
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


def inline(text: str) -> str:
    def sub(m):
        lang = m["lang"]
        if lang == "groovy":
            lang = classify(None, m["code"])[0]
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
        markers = sorted({int(m["n"]) for l in code_lines for m in ANNOTATION.finditer(l)})
        if markers:
            notes, after = annotation_list(lines, f.end + 1, f.indent)
            if missing := [n for n in markers if n not in notes]:
                page.warn(f"annotation markers {missing} have no list item after the block at line {f.start + 1}")
            if notes:
                for k, l in enumerate(code_lines):
                    code_lines[k] = ANNOTATION.sub(lambda m: f"[!annotate] {annotate(notes.get(int(m['n']), ''))}".rstrip(), l)
                lines = lines[: f.start + 1] + code_lines + [lines[f.end]] + lines[after:]
        lines[f.start] = " " * f.indent + f.marker + fence_line(f.info, "\n".join(code_lines))

    in_fence = fenced_lines(lines)
    page.body = "\n".join(l if i in in_fence else inline(l) for i, l in enumerate(lines))
