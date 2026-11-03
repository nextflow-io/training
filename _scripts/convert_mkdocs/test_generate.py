# /// script
# requires-python = ">=3.11"
# dependencies = ["pyyaml>=6"]
# ///
"""Check the sidebar, redirects and i18n generators against the MkDocs checkout.

    uv run _scripts/convert_mkdocs/test_generate.py [--src ../mkdocs-m0]
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))

from core import LOCALES, Run  # noqa: E402
from generate import GENERATORS, i18n, redirects, sidebar, source  # noqa: E402


def js_export(text: str, name: str, base: str = "/"):
    """Evaluate a generated .mjs module with node and return one export."""
    with tempfile.NamedTemporaryFile("w", suffix=".mjs", delete=False) as f:
        f.write(text)
    script = f"import('{f.name}').then(m => console.log(JSON.stringify(m.{name})))"
    out = subprocess.run(["node", "-e", script], capture_output=True, text=True, check=True, env={**os.environ, "BASE": base})
    os.unlink(f.name)
    return json.loads(out.stdout)


def find(items: list, label: str) -> dict:
    return next(i for i in items if isinstance(i, dict) and i.get("label") == label)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--src", type=Path, default=HERE.parents[2] / "mkdocs-m0")
    run = Run(src_root=ap.parse_args().src.resolve(), out_root=Path("/nonexistent"))

    assert source.slug("index.md") == ""
    assert source.slug("Foo/index.md") == "foo"
    assert source.slug("nf4_science/rnaseq/02_single-sample.md") == "nf4_science/rnaseq/02_single-sample"
    assert redirects.url("a/index.md#x", "de") == "/de/a/#x"
    assert redirects.url("https://example.com/a/", "de") == "https://example.com/a/"

    first = [g.generate(run, LOCALES) for g in GENERATORS]
    assert first == [g.generate(run, LOCALES) for g in GENERATORS], "output is not deterministic"

    bar = js_export(sidebar.generate(run, LOCALES)["src/sidebar.mjs"], "sidebar")
    assert bar[0] == "index"
    assert [g["label"] for g in bar[1:]] == ["Users", "Developers", "Setup & Help"]
    users = find(bar, "Users")
    assert users["translations"]["de"] == "Nutzer*innen" and len(users["translations"]) == 10
    run_group = find(users["items"], "Nextflow Run")
    assert run_group["collapsed"] is True and "translations" not in run_group
    overview = run_group["items"][0]
    assert (overview["label"], overview["slug"], overview["translations"]["fr"]) == ("Overview", "nextflow_run", "Aperçu")
    quests = find(find(bar, "Developers")["items"], "Side Quests")
    assert quests["items"][1] == "side_quests/orientation"
    assert [i["label"] for i in quests["items"][2:]] == [
        "Developer Tools & Tricks",
        "Deep Dives into Dataflow",
        "Modular Architecture in Action",
        "The Nextflow Extended Universe",
    ]
    plugins = find(quests["items"][-1]["items"], "Plugin Development")
    assert plugins["items"][0] == "side_quests/plugin_development"
    setup = find(bar, "Setup & Help")
    assert find(setup["items"], "Training Environment")["translations"]["de"] == "Schulungsumgebung"
    assert setup["translations"]["de"] == "Einrichtung & Hilfe"

    red_js = redirects.generate(run, LOCALES)["src/redirects.mjs"]
    red = js_export(red_js, "redirects")
    assert len(red) % len(LOCALES) == 0 and len(red) // len(LOCALES) > 40
    assert red["/training_collections/"] == "/side_quests/"
    assert red["/pt/training_collections/architects_toolkit_1/"] == "/pt/side_quests/#suggested-path-the-architects-toolkit-i"
    assert red["/ko/nextflow_run/01_basics/"] == "https://training.nextflow.io/3.6.1/nextflow_run/01_basics/"
    based = js_export(red_js, "redirects", base="/latest/")
    assert based["/pt/training_collections/"] == "/latest/pt/side_quests/"
    assert based["/ko/nextflow_run/01_basics/"] == red["/ko/nextflow_run/01_basics/"]

    strings = {k: json.loads(v) for k, v in i18n.generate(run, LOCALES).items()}
    assert sorted(strings) == sorted(f"src/content/i18n/{l}.json" for l in LOCALES)
    en, de = strings["src/content/i18n/en.json"], strings["src/content/i18n/de.json"]
    assert list(en) == list(i18n.UI_STRINGS) + list(i18n.CONSENT)
    assert de["consent.title"] == "Cookie-Zustimmung" and "<a href=" in de["consent.description"]
    assert "[Environment options](/envsetup/)" in en["courseInfo.defaultTechnicalRequirements"]

    print("generate: all checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
