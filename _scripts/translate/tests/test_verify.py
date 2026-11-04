"""Tests for structural verification of MDX translations."""

from __future__ import annotations

from pathlib import Path

import pytest

from translate.postprocess import _build_notice_span
from translate.verify import verify_translation_structural

FIXTURES = Path(__file__).parent / "fixtures"


@pytest.fixture
def pages(tmp_path):
    """Write the English fixture and its good German translation; return a writer for the latter."""
    en = tmp_path / "hello" / "page.mdx"
    de = tmp_path / "de" / "hello" / "page.mdx"
    en.parent.mkdir(parents=True)
    de.parent.mkdir(parents=True)
    en.write_text((FIXTURES / "page.en.mdx").read_text())
    good = (FIXTURES / "page.expected.mdx").read_text().replace("{notice}", _build_notice_span("de"))

    def write(old: str = "", new: str = "") -> list[str]:
        de.write_text(good.replace(old, new, 1) if old else good)
        return verify_translation_structural(en, de)

    return write


def test_good_translation_passes(pages):
    assert pages() == []


@pytest.mark.parametrize(
    "old, new, issue",
    [
        ("</details>\n", "", "MDX syntax count differs for </details>"),
        (":::tip[", ":::consejo[", "MDX syntax count differs for :::consejo"),
        ("\\{a: 1\\}", "{a: 1}", "MDX syntax count differs for bare {"),
        ("a &lt; b", "a < b", "MDX syntax count differs for bare <"),
        ('title="Terminal"', 'title="Fenster"', "Code block info string changed"),
        ("```\n\n<details", "\n<details", "Code block count mismatch"),
        ("{#2-change-the-workflow}", "{#2-aendern}", "Heading id changed"),
        ("## 1. Ausführen", "### 1. Ausführen", "Header level mismatch"),
        ("tableOfContents: false\n", "", "Frontmatter keys differ from English: tableOfContents"),
        ("02_next.mdx", "02_weiter.mdx", "Link targets differ from English"),
        ("../../hello/img/diagram.svg", "img/diagram.svg", "Image or import paths"),
        ("../../../../../../hello-nextflow", "../../../../../hello-nextflow", "Image or import paths"),
        ("---\ntitle", "```markdown\n---\ntitle", "wrapped in a code fence"),
    ],
)
def test_broken_translation_flagged(pages, old, new, issue):
    issues = pages(old, new)
    assert any(i.startswith(issue) or issue in i for i in issues), issues


def test_missing_notice_flagged(pages):
    assert "Missing AI translation notice" in pages(_build_notice_span("de") + "\n\n", "")


def test_dropped_heading_id_tolerated(pages):
    """Post-processing restores a missing id, so verification only rejects a changed one."""
    assert pages(" {#2-change-the-workflow}", "") == []
