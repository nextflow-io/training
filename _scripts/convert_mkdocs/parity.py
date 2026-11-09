# /// script
# requires-python = ">=3.11"
# ///
"""URL parity: every page in the live /latest/ sitemaps must exist in a build.

    uv run _scripts/convert_mkdocs/parity.py docs/dist

dist is one build (any BASE); its root stands for /latest/. A page passes if
dist/<path>/index.html exists (Astro writes redirect pages there too).
"""

import re
import sys
import urllib.request
from pathlib import Path

SITE = "https://training.nextflow.io/"
LOCALES = ["", "ca/", "de/", "es/", "fr/", "hi/", "it/", "ko/", "pl/", "pt/", "tr/"]


def sitemap(locale: str) -> list[str]:
    with urllib.request.urlopen(f"{SITE}latest/{locale}sitemap.xml") as r:
        xml = r.read().decode()
    return [loc.removeprefix(SITE) for loc in re.findall(r"<loc>([^<]+)</loc>", xml)]


def main() -> int:
    dist = Path(sys.argv[1])
    missing, total = [], 0
    for locale in LOCALES:
        for path in sitemap(locale):
            total += 1
            if not (dist / path / "index.html").is_file():
                missing.append(path)
    for m in missing:
        print(f"missing: /{m}")
    print(f"{total - len(missing)}/{total} sitemap URLs present")
    return 1 if missing else 0


if __name__ == "__main__":
    sys.exit(main())
