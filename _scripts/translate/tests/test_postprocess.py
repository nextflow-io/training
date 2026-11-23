"""Tests for postprocess edge cases."""

from __future__ import annotations

from translate.postprocess import (
    ensure_translation_notice,
    fix_frontmatter,
    strip_translation_notice,
)


def test_notice_inserted_after_h1():
    text = "# My Title\n\nSome content here."
    result = ensure_translation_notice(text, "pt")
    assert "ai-translation-notice" in result
    lines = result.split("\n")
    # Notice should be after the H1 line
    h1_idx = next(i for i, l in enumerate(lines) if l.startswith("# "))
    notice_idx = next(i for i, l in enumerate(lines) if "ai-translation-notice" in l)
    assert notice_idx > h1_idx


def test_notice_inserted_without_h1():
    """Files without H1 should still get the notice."""
    text = "Some content without any heading."
    result = ensure_translation_notice(text, "pt")
    assert "ai-translation-notice" in result


def test_notice_inserted_after_frontmatter_without_h1():
    """Files with frontmatter but no H1 should get notice after frontmatter."""
    text = "---\ntitle: Test\nhide:\n  - toc\n---\n\nSome content here."
    result = ensure_translation_notice(text, "pt")
    assert "ai-translation-notice" in result
    lines = result.split("\n")
    # Find the closing --- of frontmatter
    fm_close = None
    in_fm = False
    for i, line in enumerate(lines):
        if line.strip() == "---":
            if not in_fm:
                in_fm = True
            else:
                fm_close = i
                break
    assert fm_close is not None
    # Notice should be after frontmatter
    notice_idx = next(i for i, l in enumerate(lines) if "ai-translation-notice" in l)
    assert notice_idx > fm_close


def test_notice_not_duplicated():
    """Calling ensure_translation_notice twice should not duplicate."""
    text = "# My Title\n\nContent."
    result1 = ensure_translation_notice(text, "pt")
    result2 = ensure_translation_notice(result1, "pt")
    # The span contains "ai-translation-notice" in both class and icon class,
    # so count the span opening tag instead
    assert result2.count('<span class="ai-translation-notice"') == 1


def test_strip_then_reinsert_notice_no_double_blank():
    """Stripping an existing well-formed notice and re-inserting the
    canonical one should not leave a double blank line behind.

    Regression test for the bug behind PR #987's pre-commit failure:
    strip_translation_notice removed the notice line but not the blank
    line(s) around it, and the leftover blank survived re-insertion.
    """
    text = (
        "# My Title\n"
        "\n"
        '<span class="ai-translation-notice">notice</span>\n'
        "\n"
        "## Next heading\n"
    )
    stripped = "\n".join(strip_translation_notice(text.split("\n")))
    result = ensure_translation_notice(stripped, "pt")
    lines = result.split("\n")
    notice_idx = next(i for i, l in enumerate(lines) if "ai-translation-notice" in l)
    # Exactly one blank line between the notice and the following heading
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
