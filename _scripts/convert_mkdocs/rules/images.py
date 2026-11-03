"""A7. Images, Excalidraw SVG includes, figures and captions.

Images live only in the English tree. Every relative image a page uses is
published at its English path, and translated pages point there.
"""

import os
import re
from pathlib import Path

from core import fenced_lines

EXCALIDRAW = re.compile(r'<figure class="excalidraw">\s*--8<--\s*"(?P<path>[^"]+)"\s*</figure>', re.S)
MD_IMAGE = re.compile(r"(!\[[^\]]*\]\()(?P<src>[^)\s]+)((?:\s+\"[^\"]*\")?\))")
HTML_IMG = re.compile(r'(<img\b[^>]*?\bsrc=")(?P<src>[^"]+)(")')
CAPTION = re.compile(r"^/// caption\s*$")
BG_RECT = re.compile(rb'<rect x="0" y="0" [^>]*fill="#ffffff"[^>]*?(/>|></rect>)', re.I)


def strip_background(svg: bytes) -> bytes:
    """Drop the first full-canvas white rect that Excalidraw exports as the background."""
    return BG_RECT.sub(b"", svg, count=1)


def excalidraw_name(path: str) -> str:
    return path if path.endswith(".excalidraw.svg") else path[: -len(".svg")] + ".excalidraw.svg"


def publish(page, en_rel: str, transform=None) -> str:
    """Register an English asset and return its path relative to the page's output file."""
    src = page.run.docs_dir("en") / en_rel
    dst = en_rel
    if transform:
        dst = excalidraw_name(en_rel)
        page.run.asset_transforms[dst] = transform
    if not src.is_file():
        page.warn(f"missing image {en_rel}")
    page.run.assets[dst] = src
    rel = os.path.relpath(dst, os.path.dirname(page.out_rel) or ".")
    return rel if rel.startswith(".") else f"./{rel}"


def relative_src(page, src: str) -> str | None:
    if re.match(r"^([a-z][a-z\d+.-]*:|/|#|data:)", src, re.I):
        return None
    return os.path.normpath(os.path.join(os.path.dirname(page.rel), src)).replace(os.sep, "/")


def apply(page) -> None:
    def rewrite(m):
        rel = relative_src(page, m["src"])
        return m[0] if rel is None else m[1] + publish(page, rel) + m[3]

    lines = page.body.split("\n")
    in_fence = fenced_lines(lines)
    for i, line in enumerate(lines):
        if i not in in_fence:
            lines[i] = HTML_IMG.sub(rewrite, MD_IMAGE.sub(rewrite, line))

    def excalidraw(m):
        path = m["path"]
        prefix = "docs/en/docs/"
        if not path.startswith(prefix):
            page.warn(f"Excalidraw include outside docs/en/docs: {path}")
            return m[0]
        return f"![]({publish(page, path[len(prefix):], strip_background)})"

    lines = EXCALIDRAW.sub(excalidraw, "\n".join(lines)).split("\n")
    in_fence = fenced_lines(lines)
    for i, line in enumerate(lines):
        if "--8<--" in line and i not in in_fence:
            page.warn(f"unhandled snippet include: {line.strip()}")

    page.body = "\n".join(captions(lines))


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
        caption = lines[i + 1 : j]
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
