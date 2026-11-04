"""Tests for postprocess edge cases."""

from __future__ import annotations

from pathlib import Path

import yaml

from translate.postprocess import (
    _build_notice_span,
    ensure_translation_notice,
    fix_anchors,
    fix_anchors_in_language,
    fix_frontmatter,
    heading_ids,
    post_process_text,
    rebase,
    strip_translation_notice,
)

FIXTURES = Path(__file__).parent / "fixtures"


def make_tree(root: Path) -> None:
    """English page hello/page.mdx with an image and a sibling page, and a de sibling."""
    (root / "hello" / "img").mkdir(parents=True)
    (root / "hello" / "img" / "diagram.svg").write_text("<svg/>")
    (root / "hello" / "page.mdx").write_text((FIXTURES / "page.en.mdx").read_text())
    (root / "hello" / "02_next.mdx").write_text("## 1. Start\n")
    (root / "de" / "hello").mkdir(parents=True)
    (root / "de" / "hello" / "02_next.mdx").write_text("## 1. Anfangen\n")


def test_post_process_fixture(tmp_path):
    """A German LLM output with typical MDX mistakes is repaired from the English source."""
    make_tree(tmp_path)
    en = (FIXTURES / "page.en.mdx").read_text()
    llm = (FIXTURES / "page.llm.mdx").read_text()
    expected = (FIXTURES / "page.expected.mdx").read_text()
    expected = expected.replace("{notice}", _build_notice_span("de"))

    result, problems = post_process_text(llm, en, "hello/page.mdx", "de", root=tmp_path)

    assert problems == []
    assert result == expected
    # Idempotent
    assert post_process_text(result, en, "hello/page.mdx", "de", root=tmp_path)[0] == result


def test_post_process_mismatch_keeps_other_fixes(tmp_path):
    """A missing code block skips the code fix but still restores tags and asides."""
    make_tree(tmp_path)
    en = (FIXTURES / "page.en.mdx").read_text()
    llm = (FIXTURES / "page.llm.mdx").read_text()
    llm = llm.replace("<Code code={hello_config} lang=\"nextflow\" title=\"nextflow.config\" />\n", "")
    llm = llm.replace('```bash title="Terminal-Fenster"', "```bash")

    result, problems = post_process_text(llm, en, "hello/page.mdx", "de", root=tmp_path)

    assert any(p.startswith("Tags differ") for p in problems)
    assert '<details class="success">' not in result  # tag fix skipped as a whole
    assert ":::tip[Übung]" in result
    assert '```bash title="Terminal"' in result


def test_rebase():
    assert rebase("img/x.svg", "hello/page.mdx", "de") == "../../hello/img/x.svg"
    assert rebase("../assets/x.png", "a/b/page.mdx", "pt") == "../../../a/assets/x.png"
    assert rebase("https://x.org/a.png", "hello/page.mdx", "de") == "https://x.org/a.png"
    assert rebase("../../x.nf?raw", "index.mdx", "de") == "../../../x.nf?raw"


def test_heading_ids_match_starlight():
    lines = [
        "## 1. Ausführen",
        "### Takeaway",
        "### Takeaway",
        "## Use `x{:nextflow}` and [links](a.mdx)",
        "## 2. Pinned {#old-id}",
    ]
    assert [i for _, i in heading_ids(lines)] == [
        "1-ausführen",
        "takeaway",
        "takeaway-1",
        "use-x-and-links",
        "old-id",
    ]


def test_fix_anchors_maps_through_english_headings(tmp_path):
    en = ["## 1. Run it", "", "See [this](#1-run-it) and [pinned](#keep)", "## 2. Two {#keep}"]
    trans = ["## 1. Ausführen", "", "Siehe [dies](#1-run-it) und [fest](#keep)", "## 2. Zwei {#keep}"]
    result = fix_anchors(trans, en, "page.mdx", "de", root=tmp_path)
    assert result[2] == "Siehe [dies](#1-ausführen) und [fest](#keep)"


def test_fix_anchors_in_language_repairs_inbound_links(tmp_path):
    """A re-translated page changed its heading ids; pages linking to it are re-pointed."""
    (tmp_path / "de").mkdir()
    (tmp_path / "a.mdx").write_text("## Intro\n\n[b](b.mdx#1-start)\n")
    (tmp_path / "b.mdx").write_text("## 1. Start\n")
    (tmp_path / "de" / "a.mdx").write_text("## Einleitung\n\n[b](b.mdx#1-anfang)\n")
    (tmp_path / "de" / "b.mdx").write_text("## 1. Beginnen\n")

    assert fix_anchors_in_language("de", root=tmp_path) == 1
    assert "[b](b.mdx#1-beginnen)" in (tmp_path / "de" / "a.mdx").read_text()
    assert fix_anchors_in_language("de", root=tmp_path) == 0


def test_fix_frontmatter_merges_translated_values():
    source = (
        "---\ntitle: 'Part 1: Setup'\ndescription: Set up.\nsidebar:\n  label: Setup\n"
        "tableOfContents: false\nadditional_information:\n  learning_objectives:\n  - Run it\n---\n\nBody\n"
    )
    translation = (
        "---\ntitle: \"Partie 1 : l'environnement\"\ndescription: Installer.\nsidebar:\n  label: Installation\n"
        "tableOfContents: true\nextra: x\nadditional_information:\n  learning_objectives:\n  - L'exécuter\n---\nCorps\n"
    )
    result = fix_frontmatter(translation, source, "fr")
    fm = yaml.safe_load(result.split("---\n")[1])
    assert fm == {
        "title": "Partie 1 : l'environnement",
        "description": "Installer.",
        "sidebar": {"label": "Installation"},
        "tableOfContents": False,
        "additional_information": {"learning_objectives": ["L'exécuter"]},
    }
    assert result.endswith("---\n\nCorps\n")


def test_notice_inserted_after_frontmatter_and_imports():
    text = "---\ntitle: Test\n---\n\nimport x from './x.nf?raw';\n\nSome content here."
    lines = ensure_translation_notice(text, "pt").split("\n")
    assert lines[5] == ""
    assert lines[6].startswith('<span class="ai-translation-notice">')
    assert lines[7:] == ["", "Some content here."]


def test_notice_inserted_without_frontmatter():
    """Files without frontmatter should still get the notice."""
    text = "Some content without any heading."
    result = ensure_translation_notice(text, "pt")
    assert result.startswith('<span class="ai-translation-notice">')


def test_homepage_notice_before_first_h2():
    text = "---\ntitle: Home\n---\n\n<CardGrid>\n</CardGrid>\n\n## Catalog\n"
    lines = ensure_translation_notice(text, "de", is_homepage=True).split("\n")
    start = lines.index(next(l for l in lines if l.startswith(":::note[")))
    end = lines.index(":::", start)
    assert lines[start - 1] == "" and lines[end + 1] == "" and lines[end + 2] == "## Catalog"
    assert "TRANSLATING.md" in "\n".join(lines[start:end])


def test_notice_not_duplicated():
    """Calling ensure_translation_notice twice should not duplicate."""
    text = "---\ntitle: T\n---\n\nContent."
    result1 = ensure_translation_notice(text, "pt")
    assert ensure_translation_notice(result1, "pt") == result1
    homepage = ensure_translation_notice(text + "\n\n## H\n", "pt", is_homepage=True)
    assert ensure_translation_notice(homepage, "pt", is_homepage=True) == homepage


def test_strip_then_reinsert_notice_no_double_blank():
    """Stripping an existing well-formed notice and re-inserting the
    canonical one should not leave a double blank line behind.

    Regression test for the bug behind PR #987's pre-commit failure.
    """
    text = (
        "---\ntitle: T\n---\n"
        "\n"
        '<span class="ai-translation-notice">notice</span>\n'
        "\n"
        "## Next heading\n"
    )
    stripped = "\n".join(strip_translation_notice(text.split("\n")))
    result = ensure_translation_notice(stripped, "pt")
    lines = result.split("\n")
    notice_idx = next(i for i, l in enumerate(lines) if "ai-translation-notice" in l)
    # Exactly one blank line on each side of the notice
    assert lines[notice_idx - 1] == "" and lines[notice_idx - 2] == "---"
    assert lines[notice_idx + 1].strip() == ""
    assert lines[notice_idx + 2].strip() != ""


def test_fix_frontmatter_strips_unclosed_stray_dashes():
    """An unclosed leading '---' the LLM emitted when the English source has
    no frontmatter should be stripped, not left as stray content.

    Regression test for docs/tr/docs/config_exec/next_steps.md in PR #987.
    """
    source = "# Course summary\n\nSome content.\n"
    translation = "---\n\n# Resume du cours\n\nContenu.\n"
    result = fix_frontmatter(translation, source, "fr")
    assert not result.lstrip().startswith("---")
    assert result.lstrip().startswith("# Resume du cours")
