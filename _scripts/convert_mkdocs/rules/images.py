"""A7. Images, Excalidraw SVG includes, Mermaid, code includes, figures and captions.

Images live only in the English tree. Every relative image a page uses is
published at its English path, and translated pages point there.
"""

import base64
import hashlib
import os
import re
from pathlib import Path

from core import dedent, fenced_lines, fences

from .fences import fence_line

CONTENT_ROOT = "docs/src/content/docs"  # repo-relative; `?raw` imports climb out of it to the repo root
MERMAID_CACHE = Path(__file__).resolve().parents[1] / "mermaid"

EXCALIDRAW = re.compile(r'<figure class="excalidraw">\s*--8<--\s*"(?P<path>[^"]+)"\s*</figure>', re.S)
INCLUDE = re.compile(r'^--8<--\s*"(?P<path>[^"]+)"\s*$')
MD_IMAGE = re.compile(r'!\[(?P<alt>[^\]]*)\]\((?P<src>[^)\s]+)(?P<title>\s+"[^"]*")?\)(?:\{(?P<attrs>[^}\n]*)\})?')
ATTR = re.compile(r'(?P<key>[\w-]+)=(?:"(?P<q>[^"]*)"|(?P<v>[^\s}]+))')
THEME_ONLY = {"only-light": "dark:sl-hidden", "only-dark": "light:sl-hidden"}
VIDEO = re.compile(r'<div\b[^>]*>\s*(?P<iframe><iframe\b[^>]*\bsrc="https://www\.youtube\.com/[^"]*"[^>]*>\s*</iframe>)\s*</div>')
INLINE_CODE = re.compile(r"(`+).*?(?<!`)\1(?!`)")
CAPTION = re.compile(r"^/// caption\s*$")
FIRST_RECT = re.compile(rb"<rect\b[^>]*?(?:/>|></rect>)")
BG_RECT = re.compile(rb'<rect (?:x="0" y="0" )?width="(?P<w>[\d.]+)" height="(?P<h>[\d.]+)" fill="#ffffff"\s*(?:/>|></rect>)')
VIEWBOX = re.compile(rb'<svg\b[^>]*\bviewBox="0 0 (?P<w>[\d.]+) (?P<h>[\d.]+)"')
FONT_FACE = re.compile(rb"@font-face\s*\{[^}]*\}")
FACE_FAMILY = re.compile(rb"font-family:\s*[\"']?(?P<f>[^\"';]+)")
FACE_URL = re.compile(rb"url\(\s*[\"']?(?P<u>[^\"')]+)")
USED_FAMILY = re.compile(rb'font-family="(?P<f>[^",]+)')
FONT_DIR = "assets/fonts/excalidraw"  # in the English tree; the fonts MkDocs served for inline Excalidraw SVGs
# No local Excalifont: Virgil, the older Excalidraw hand-drawn font, stands in for it.
FONT_FILES = {b"Virgil": "Virgil.woff2", b"Cascadia": "Cascadia.woff2", b"Excalifont": "Virgil.woff2"}


def strip_background(svg: bytes) -> bytes:
    """Drop the white full-canvas rect that Excalidraw draws first, so dark mode can invert the drawing."""
    rect, box = FIRST_RECT.search(svg), VIEWBOX.search(svg)
    bg = rect and box and BG_RECT.fullmatch(rect[0])
    if bg and all(abs(float(bg[k]) - float(box[k])) < 0.01 for k in "wh"):
        return svg[: rect.start()] + svg[rect.end() :]
    return svg


def embed_fonts(svg: bytes, font_dir: Path) -> bytes:
    """An <img> loads no external fonts: embed the ones the drawing uses and drop the other external faces."""
    used = {f.strip() for f in USED_FAMILY.findall(svg)}

    def face(m):
        family, url = FACE_FAMILY.search(m[0]), FACE_URL.search(m[0])
        if not family or not url or url["u"].startswith(b"data:"):
            return m[0]
        name = family["f"].strip()
        if name not in used or name not in FONT_FILES:
            return b""
        data = base64.b64encode((font_dir / FONT_FILES[name]).read_bytes())
        return b"@font-face { font-family: " + name + b"; src: url(data:font/woff2;base64," + data + b"); }"

    return FONT_FACE.sub(face, svg)


def excalidraw_name(path: str) -> str:
    return path if path.endswith(".excalidraw.svg") else path[: -len(".svg")] + ".excalidraw.svg"


def published(page, dst: str) -> str:
    rel = os.path.relpath(dst, os.path.dirname(page.out_rel) or ".")
    return rel if rel.startswith(".") else f"./{rel}"


def markdown_image(page, alt: str, dst: str, src: Path, transform=None, title: str = "") -> str:
    """Register the asset and return the image. A missing file is dropped, since it would fail the build."""
    if not src.is_file():
        page.warn(f"missing image dropped: {dst}")
        return ""
    page.run.assets[dst] = src
    if transform:
        page.run.asset_transforms[dst] = transform
    return f"![{alt}]({published(page, dst)}{title})"


def en_rel(page, src: str) -> str | None:
    """The English-tree path of a relative image; None for URLs and absolute paths."""
    if re.match(r"^([a-z][a-z\d+.-]*:|/|#)", src, re.I):
        return None
    return os.path.normpath(os.path.join(os.path.dirname(page.rel), src)).replace(os.sep, "/")


def image(page, m) -> str:
    src, _, frag = m["src"].partition("#")
    rel = en_rel(page, src)
    if not rel:
        out = f"![{m['alt']}]({m['src']}{m['title'] or ''})"
    elif not (out := markdown_image(page, m["alt"], rel, page.run.docs_dir("en") / rel, title=m["title"] or "")):
        return ""
    classes = [THEME_ONLY[frag]] if rel and frag in THEME_ONLY else []
    style = ""
    for a in ATTR.finditer(m["attrs"] or ""):
        key, value = a["key"], a["v"] or a["q"]
        if key == "align" and value in ("left", "right"):
            classes.append(f"align-{value}")
        elif key == "width":
            style = f"display: inline-block; width: {value}{'px' if value.isdigit() else ''}"
        else:
            page.warn(f"image attribute dropped: {a[0]}")
    if ATTR.sub("", m["attrs"] or "").strip():
        page.warn(f"image attribute dropped: {{{m['attrs']}}}")
    if classes or style:
        out = "<span" + (f' class="{" ".join(classes)}"' if classes else "") + (f' style="{style}"' if style else "") + f">{out}</span>"
    return out


def mermaid_key(code: list[str]) -> str:
    return hashlib.sha256(("\n".join(code) + "\n").encode()).hexdigest()[:16]


def mermaid(page, code: list[str]) -> list[str] | None:
    key = mermaid_key(code)
    svg = MERMAID_CACHE / f"{key}.svg"
    if not svg.is_file():
        page.warn(f"no pre-rendered Mermaid SVG {key}.svg; run _scripts/convert_mkdocs/mermaid/render.py")
        return None
    return [markdown_image(page, "", f"{os.path.dirname(page.rel)}/img/mermaid-{key}.svg", svg)]


def include(page, info: str, path: str) -> list[str]:
    """A fence that only includes a repo file becomes <Code> fed by a `?raw` import of that file."""
    src = page.run.src_root / path
    if not src.is_file():
        page.warn(f"missing code include {path}")
    lang, _, meta = fence_line(info, src.read_text(encoding="utf-8") if src.is_file() else "").partition(" ")
    title = re.search(r'\btitle="([^"]*)"', meta)
    meta = re.sub(r'\s*\btitle="[^"]*"', "", meta).strip()
    name = re.sub(r"\W+", "_", path).strip("_")
    rel = os.path.relpath(path, os.path.join(CONTENT_ROOT, os.path.dirname(page.out_rel))).replace(os.sep, "/")
    page.imports.append(f"import {name} from '{rel}?raw';")
    q = "'" if '"' in meta else '"'
    props = f'code={{{name}}} lang="{lang}"' + (f' title="{title[1]}"' if title else "") + (f" meta={q}{meta}{q}" if meta else "")
    return [f"<Code {props} />"]


def blocks(page, lines: list[str]) -> list[str]:
    for f in reversed(fences(lines)):
        code = dedent(lines[f.start + 1 : f.end], f.indent)
        new = None
        if f.lang == "mermaid":
            new = mermaid(page, code)
        elif len(code) == 1 and (m := INCLUDE.match(code[0].strip())):
            new = include(page, f.info, m["path"])
        if new:
            lines[f.start : f.end + 1] = [" " * f.indent + l for l in new]
    return lines


def excalidraw(page, m) -> str:
    path, prefix = m["path"], "docs/en/docs/"
    if not path.startswith(prefix):
        page.warn(f"Excalidraw include outside docs/en/docs: {path}")
        return m[0]
    rel = path[len(prefix) :]
    src = page.run.docs_dir("en") / rel
    font_dir = page.run.docs_dir("en") / FONT_DIR

    def transform(svg: bytes) -> bytes:
        return embed_fonts(strip_background(svg), font_dir)

    if page.lang == "en" and src.is_file():
        svg = src.read_bytes()
        if strip_background(svg) == svg:
            page.warn(f"no background rect to strip in {rel}")
        if re.search(rb"url\(\s*[\"']?https?:", transform(svg)):
            page.warn(f"external URL left in {rel}")
    return markdown_image(page, "", excalidraw_name(rel), src, transform)


def video(m) -> str:
    iframe = re.sub(r'\s+data-ruffle-polyfilled="[^"]*"', "", m["iframe"]).replace('allowfullscreen=""', "allowfullscreen")
    return f'<div class="video-wrapper">\n  {iframe}\n</div>'


def outside_code(line: str, fn) -> str:
    out, i = [], 0
    for m in INLINE_CODE.finditer(line):
        out += [fn(line[i : m.start()]), m[0]]
        i = m.end()
    return "".join(out) + fn(line[i:])


def apply(page) -> None:
    lines = captions(page.body.split("\n"))
    in_fence = fenced_lines(lines)
    for i, line in enumerate(lines):
        if i not in in_fence:
            lines[i] = outside_code(line, lambda s: MD_IMAGE.sub(lambda m: image(page, m), s))
    lines = blocks(page, lines)

    body = EXCALIDRAW.sub(lambda m: excalidraw(page, m), "\n".join(lines))
    lines = VIDEO.sub(video, body).split("\n")
    for line in lines:
        if "--8<--" in line:
            page.warn(f"unhandled snippet include: {line.strip()}")
        if re.search(r"<img\b|<figure\s", line):
            page.warn(f"raw HTML image left as is: {line.strip()[:80]}")

    page.body = "\n".join(lines)


def captions(lines: list[str]) -> list[str]:
    """`/// caption` after a block wraps that block and the caption in <figure>."""
    out: list[str] = []
    i = 0
    while i < len(lines):
        if not CAPTION.match(lines[i]):
            out.append(lines[i])
            i += 1
            continue
        j = i + 1
        while j < len(lines) and lines[j].strip() != "///":
            j += 1
        caption = "\n".join(lines[i + 1 : j]).strip("\n").split("\n")
        # the captioned block is the previous paragraph or HTML block
        while out and not out[-1].strip():
            out.pop()
        k = len(out)
        while k > 0 and out[k - 1].strip():
            k -= 1
        block = out[k:]
        del out[k:]
        out += ["<figure>", *block, "<figcaption>", "", *caption, "", "</figcaption>", "</figure>"]
        i = j + 1
    return out
