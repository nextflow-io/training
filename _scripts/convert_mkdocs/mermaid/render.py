# /// script
# requires-python = ">=3.11"
# dependencies = ["pyyaml>=6"]
# ///
"""Pre-render every Mermaid diagram in the MkDocs source to <key>.svg in this folder.

    uv run _scripts/convert_mkdocs/mermaid/render.py [--src ../mkdocs-m0]

The converter stays offline: it only reads the SVGs cached here, keyed by a hash of the
diagram source. Run this when it warns about a missing SVG, then commit the new files.
Needs Node and a Chrome that Puppeteer can find (set PUPPETEER_EXECUTABLE_PATH if needed).
"""

import argparse
import re
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE.parent))

from core import LOCALES, dedent, fences  # noqa: E402
from rules.images import mermaid_key  # noqa: E402

MMDC = "@mermaid-js/mermaid-cli@12.2.0"


def fixed_size(svg: str) -> str:
    """Mermaid sizes the SVG to its container (width="100%"); an <img> needs the drawing's own size."""
    w, h = re.search(r'viewBox="[-\d.]+ [-\d.]+ ([\d.]+) ([\d.]+)"', svg).groups()
    return svg.replace('width="100%"', f'width="{w}" height="{h}"', 1)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--src", type=Path, default=HERE.parents[2].parent / "mkdocs-m0")
    args = ap.parse_args()
    for lang in LOCALES:
        for md in sorted((args.src / "docs" / lang / "docs").rglob("*.md")):
            lines = md.read_text(encoding="utf-8").split("\n")
            for f in fences(lines):
                if f.lang != "mermaid":
                    continue
                code = dedent(lines[f.start + 1 : f.end], f.indent)
                key = mermaid_key(code)
                if (HERE / f"{key}.svg").is_file():
                    continue
                (HERE / f"{key}.mmd").write_text("\n".join(code) + "\n", encoding="utf-8")
                print(f"rendering {key} from {md.relative_to(args.src)}", file=sys.stderr)
                subprocess.run(
                    ["npx", "-y", MMDC, "-q", "-b", "white", "-i", f"{key}.mmd", "-o", f"{key}.svg"],
                    cwd=HERE,
                    check=True,
                )
                svg = HERE / f"{key}.svg"
                svg.write_text(fixed_size(svg.read_text(encoding="utf-8")), encoding="utf-8")
    return 0


if __name__ == "__main__":
    sys.exit(main())
