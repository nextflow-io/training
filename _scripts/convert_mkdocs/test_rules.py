# /// script
# requires-python = ">=3.11"
# dependencies = ["pyyaml>=6"]
# ///
"""Run converter rule fixtures: MkDocs snippet in, MDX snippet out.

    uv run _scripts/convert_mkdocs/test_rules.py                 # every rule
    uv run _scripts/convert_mkdocs/test_rules.py fences tabs     # some rules
    uv run _scripts/convert_mkdocs/test_rules.py --update fences # rewrite expected output

Fixtures live in tests/<rule>/<case>.md (input) and tests/<case>.mdx (expected).
Only that one rule runs on the input. If the expected file starts with `---`,
the whole rendered page is compared; otherwise only the body.

An optional first line sets the page the snippet pretends to come from:
    <!-- fixture: lang=pt rel=hello_nextflow/01_hello_world.md -->
Rules that look at other files read them from --src (default ../mkdocs-m0).
"""

from __future__ import annotations

import argparse
import difflib
import importlib
import re
import sys
from pathlib import Path

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))

from core import Page, Run, split_frontmatter  # noqa: E402

HEADER = re.compile(r"^<!-- fixture:(?P<opts>.*?)-->\n")


def run_case(rule, src: Path, inp: Path) -> str:
    text = inp.read_text(encoding="utf-8")
    opts = {"lang": "en", "rel": "fixture/page.md"}
    m = HEADER.match(text)
    if m:
        opts.update(kv.split("=", 1) for kv in m["opts"].split())
        text = text[m.end() :]
    meta, body = split_frontmatter(text)
    page = Page(run=Run(src_root=src, out_root=Path("/nonexistent")), lang=opts["lang"], rel=opts["rel"], meta=meta, body=body)
    rule.apply(page)
    expected = inp.with_suffix(".mdx")
    full = expected.is_file() and expected.read_text(encoding="utf-8").startswith("---")
    out = page.render() if full else page.body.strip("\n") + "\n"
    for w in page.run.warnings:
        out += f"<!-- warning: {w.split(': ', 1)[1]} -->\n"
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("rules", nargs="*")
    ap.add_argument("--update", action="store_true")
    ap.add_argument("--src", type=Path, default=HERE.parents[2] / "mkdocs-m0")
    args = ap.parse_args()

    dirs = sorted(d for d in (HERE / "tests").iterdir() if d.is_dir() and (not args.rules or d.name in args.rules))
    failed = passed = 0
    for d in dirs:
        rule = importlib.import_module(f"rules.{d.name}")
        for inp in sorted(d.glob("*.md")):
            got = run_case(rule, args.src.resolve(), inp)
            exp_path = inp.with_suffix(".mdx")
            if args.update:
                exp_path.write_text(got, encoding="utf-8")
                continue
            exp = exp_path.read_text(encoding="utf-8") if exp_path.is_file() else ""
            if got == exp:
                passed += 1
                continue
            failed += 1
            print(f"FAIL {d.name}/{inp.name}")
            sys.stdout.writelines(difflib.unified_diff(exp.splitlines(True), got.splitlines(True), "expected", "got"))
    print(f"{passed} passed, {failed} failed")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
