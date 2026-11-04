"""Tests for heading-based chunking."""

from __future__ import annotations

from translate.chunking import Chunk, ChunkDiff, diff_chunks, reassemble_chunks, split_into_chunks


# ---------------------------------------------------------------------------
# split_into_chunks
# ---------------------------------------------------------------------------


def test_split_basic():
    text = "# Title\n\nIntro text.\n\n## Section 1\n\nContent 1.\n\n## Section 2\n\nContent 2."
    chunks = split_into_chunks(text)
    assert len(chunks) == 3
    assert chunks[0].heading == ""  # preamble
    assert "# Title" in chunks[0].content
    assert chunks[1].heading == "## Section 1"
    assert "Content 1." in chunks[1].content
    assert chunks[2].heading == "## Section 2"
    assert "Content 2." in chunks[2].content


def test_split_no_h2():
    text = "# Title\n\nJust a single section with no ## headings."
    chunks = split_into_chunks(text)
    assert len(chunks) == 1
    assert chunks[0].heading == ""


def test_split_frontmatter():
    text = "---\ntitle: Test\n---\n\n# Title\n\n## Section 1\n\nContent."
    chunks = split_into_chunks(text)
    assert len(chunks) == 2
    assert "---" in chunks[0].content
    assert "title: Test" in chunks[0].content
    assert chunks[1].heading == "## Section 1"


def test_split_code_block_with_h2():
    """## inside code blocks should not start a new chunk."""
    text = (
        "# Title\n\n"
        "## Section 1\n\n"
        "```markdown\n"
        "## This is not a heading\n"
        "```\n\n"
        "## Section 2\n\n"
        "Content."
    )
    chunks = split_into_chunks(text)
    assert len(chunks) == 3
    assert chunks[1].heading == "## Section 1"
    assert "## This is not a heading" in chunks[1].content
    assert chunks[2].heading == "## Section 2"


def test_split_nested_headings():
    """### and #### should not start new chunks."""
    text = (
        "## Section 1\n\n"
        "### Subsection 1.1\n\nContent.\n\n"
        "#### Sub-subsection\n\nMore.\n\n"
        "## Section 2\n\nContent 2."
    )
    chunks = split_into_chunks(text)
    assert len(chunks) == 3  # preamble (empty), section 1, section 2
    assert "### Subsection 1.1" in chunks[1].content
    assert "#### Sub-subsection" in chunks[1].content


def test_split_preserves_content_hash():
    text = "## Section 1\n\nContent."
    chunks = split_into_chunks(text)
    assert chunks[1].content_hash != ""
    # Same content should produce same hash
    chunks2 = split_into_chunks(text)
    assert chunks[1].content_hash == chunks2[1].content_hash


# ---------------------------------------------------------------------------
# reassemble_chunks
# ---------------------------------------------------------------------------


def test_reassemble_roundtrip():
    text = "# Title\n\nIntro.\n\n## Section 1\n\nContent 1.\n\n## Section 2\n\nContent 2."
    chunks = split_into_chunks(text)
    result = reassemble_chunks(chunks)
    assert result == text


def test_reassemble_preserves_blank_lines():
    text = "# Title\n\n\n## Section 1\n\nLine 1.\n\nLine 2.\n\n## Section 2\n\nEnd."
    chunks = split_into_chunks(text)
    result = reassemble_chunks(chunks)
    assert result == text


# ---------------------------------------------------------------------------
# diff_chunks
# ---------------------------------------------------------------------------


def test_diff_no_changes():
    text = "# Title\n\n## Sec 1\n\nContent.\n\n## Sec 2\n\nMore."
    old = split_into_chunks(text)
    new = split_into_chunks(text)
    diff = diff_chunks(old, new)
    assert len(diff.unchanged) == 3
    assert diff.modified == []
    assert diff.added == []
    assert diff.removed == []


def test_diff_single_edit():
    old_text = "# Title\n\n## Sec 1\n\nOld content.\n\n## Sec 2\n\nMore."
    new_text = "# Title\n\n## Sec 1\n\nNew content.\n\n## Sec 2\n\nMore."
    old = split_into_chunks(old_text)
    new = split_into_chunks(new_text)
    diff = diff_chunks(old, new)
    assert len(diff.unchanged) == 2  # preamble + sec 2
    assert len(diff.modified) == 1  # sec 1
    assert diff.modified[0] == (1, 1)
    assert diff.added == []
    assert diff.removed == []


def test_diff_added_section():
    old_text = "# Title\n\n## Sec 1\n\nContent."
    new_text = "# Title\n\n## Sec 1\n\nContent.\n\n## Sec 2\n\nNew section."
    old = split_into_chunks(old_text)
    new = split_into_chunks(new_text)
    diff = diff_chunks(old, new)
    assert len(diff.unchanged) == 2  # preamble + sec 1
    assert len(diff.added) == 1
    assert diff.added[0] == 2  # new sec 2 at index 2
    assert diff.removed == []


def test_diff_removed_section():
    old_text = "# Title\n\n## Sec 1\n\nContent.\n\n## Sec 2\n\nMore."
    new_text = "# Title\n\n## Sec 1\n\nContent."
    old = split_into_chunks(old_text)
    new = split_into_chunks(new_text)
    diff = diff_chunks(old, new)
    assert len(diff.unchanged) == 2  # preamble + sec 1
    assert len(diff.removed) == 1
    assert diff.removed[0] == 2  # old sec 2
    assert diff.added == []


def test_diff_reorder():
    """Swapping two sections should detect them as unchanged (via content hash)."""
    old_text = "# Title\n\n## Sec A\n\nContent A.\n\n## Sec B\n\nContent B."
    new_text = "# Title\n\n## Sec B\n\nContent B.\n\n## Sec A\n\nContent A."
    old = split_into_chunks(old_text)
    new = split_into_chunks(new_text)
    diff = diff_chunks(old, new)
    # Both sections should match (headings match, and content hashes match for swapped ones)
    assert len(diff.unchanged) == 3  # preamble + both sections (matched by heading)
    assert diff.modified == []
    assert diff.added == []
    assert diff.removed == []


def test_diff_heading_rename():
    """Renaming a heading with similar text should match via fuzzy matching."""
    old_text = "# Title\n\n## 1. Hello World\n\nContent."
    new_text = "# Title\n\n## 1. Hello World Example\n\nContent."
    old = split_into_chunks(old_text)
    new = split_into_chunks(new_text)
    diff = diff_chunks(old, new)
    # Should match via fuzzy heading (>0.8 similarity)
    # The section is "modified" because the heading text changed (content hash differs)
    assert len(diff.unchanged) == 1  # preamble
    assert len(diff.modified) == 1  # the section (heading renamed)
    assert diff.added == []
    assert diff.removed == []


def test_diff_heading_rename_with_content_change():
    old_text = "# Title\n\n## 1. Hello World\n\nOld content."
    new_text = "# Title\n\n## 1. Hello World Example\n\nNew content."
    old = split_into_chunks(old_text)
    new = split_into_chunks(new_text)
    diff = diff_chunks(old, new)
    assert len(diff.unchanged) == 1  # preamble
    assert len(diff.modified) == 1  # the section (heading fuzzy matched, content changed)
