#!/usr/bin/env python
"""
Check that Excalidraw-generated SVGs still embed their editable scene data.

An SVG exported from Excalidraw without "Embed scene" has the Excalidraw header
but no payload, so it cannot be reopened for editing.

Usage: check_excalidraw.py [FILES...]   (defaults to every SVG under docs/en)
"""

import sys
from pathlib import Path

HEADER = "svg-source:excalidraw"
PAYLOAD = "payload-type:application/vnd.excalidraw+json"


def main(args: list[str]) -> int:
    files = [Path(a) for a in args] or sorted(Path("docs/en").rglob("*.svg"))
    bad = []
    for f in files:
        if f.suffix != ".svg" or not f.is_file():
            continue
        text = f.read_text(errors="ignore")
        if (HEADER in text or f.name.endswith(".excalidraw.svg")) and PAYLOAD not in text:
            bad.append(f)
    for f in bad:
        print(f"{f}: Excalidraw SVG has no embedded scene data")
    if bad:
        print(
            f"\n{len(bad)} file(s) not editable. Re-export from Excalidraw with 'Embed scene' ticked."
        )
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
