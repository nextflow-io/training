# /// script
# requires-python = ">=3.11"
# dependencies = ["pyyaml>=6"]
# ///
"""Convert the MkDocs training content to the Starlight content tree.

    uv run _scripts/convert_mkdocs/convert.py --src ../mkdocs-m0 --out docs/src/content/docs
    uv run _scripts/convert_mkdocs/convert.py --src ../mkdocs-m0 --out ... --lang en,pt --only hello_nextflow/index.md
    uv run _scripts/convert_mkdocs/convert.py --src ../mkdocs-m0 --renames   # print `git mv` pairs for commit 1

The output is deterministic: running twice gives the same files.
"""

from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from core import LOCALES, Page, Run, split_frontmatter  # noqa: E402
from generate import GENERATORS  # noqa: E402
from rules import RULES  # noqa: E402

# Not pages: MkDocs skips README.md when index.md exists, and assets/ holds theme files.
SKIP = {"README.md"}
SKIP_DIRS = {"assets"}


def source_pages(run: Run, langs: list[str], only: list[str]) -> list[tuple[str, str]]:
    out = []
    for lang in langs:
        root = run.docs_dir(lang)
        for p in sorted(root.rglob("*.md")):
            rel = p.relative_to(root).as_posix()
            if p.name in SKIP or SKIP_DIRS & set(Path(rel).parts[:-1]):
                continue
            if only and rel not in only:
                continue
            out.append((lang, rel))
    return out


def convert_page(run: Run, lang: str, rel: str) -> Page:
    text = (run.docs_dir(lang) / rel).read_text(encoding="utf-8")
    meta, body = split_frontmatter(text)
    page = Page(run=run, lang=lang, rel=rel, meta=meta, body=body)
    for rule in RULES:
        rule.apply(page)
    return page


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--src", required=True, type=Path, help="MkDocs checkout (worktree of master)")
    ap.add_argument("--out", type=Path, default=Path("docs/src/content/docs"))
    ap.add_argument("--app", type=Path, default=Path("docs"), help="Astro app root for generated config")
    ap.add_argument("--no-generate", action="store_true", help="skip sidebar, redirects and i18n generators")
    ap.add_argument("--lang", default=",".join(LOCALES))
    ap.add_argument("--only", default="", help="comma-separated page paths relative to docs/<lang>/docs")
    ap.add_argument("--clean", action="store_true", help="delete --out first")
    ap.add_argument("--renames", action="store_true", help="print `src dst` pairs for the rename-only commit")
    args = ap.parse_args()

    langs = args.lang.split(",")
    only = [o for o in args.only.split(",") if o]
    run = Run(src_root=args.src.resolve(), out_root=args.out.resolve())

    pages = [convert_page(run, lang, rel) for lang, rel in source_pages(run, langs, only)]

    if args.renames:
        for p in pages:
            print(f"docs/{p.lang}/docs/{p.rel}\tdocs/src/content/docs/{p.out_rel}")
        for dst, src in sorted(run.assets.items()):
            # Assets made by the converter (pre-rendered Mermaid) are new files, added with its output.
            if src.is_relative_to(run.src_root):
                print(f"{src.relative_to(run.src_root).as_posix()}\tdocs/src/content/docs/{dst}")
        return 0

    if args.clean and run.out_root.exists():
        shutil.rmtree(run.out_root)
    for p in pages:
        out = run.out_root / p.out_rel
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(p.render(), encoding="utf-8")
    for dst, src in sorted(run.assets.items()):
        out = run.out_root / dst
        out.parent.mkdir(parents=True, exist_ok=True)
        data = src.read_bytes()
        transform = run.asset_transforms.get(dst)
        out.write_bytes(transform(data) if transform else data)

    if not args.no_generate:
        for gen in GENERATORS:
            # Always every locale: the sidebar and redirects describe the whole site.
            for rel, content in sorted(gen.generate(run, LOCALES).items()):
                out = args.app / rel
                out.parent.mkdir(parents=True, exist_ok=True)
                out.write_text(content, encoding="utf-8")

    for w in run.warnings:
        print(f"warning: {w}", file=sys.stderr)
    print(f"{len(pages)} pages, {len(run.assets)} assets, {len(run.warnings)} warnings", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
