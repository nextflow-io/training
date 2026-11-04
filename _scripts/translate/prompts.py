"""Prompt loading and building for translation requests."""

from __future__ import annotations

from functools import lru_cache

from .config import DOCS_ROOT, SCRIPTS_DIR, ConfigError


@lru_cache
def get_general_prompt() -> str:
    """Load the general translation prompt."""
    path = SCRIPTS_DIR / "general-llm-prompt.md"
    if not path.exists():
        raise ConfigError(f"General prompt not found: {path}")
    return path.read_text(encoding="utf-8")


@lru_cache
def get_lang_prompt(lang: str) -> str:
    """Load language-specific prompt."""
    path = DOCS_ROOT / lang / "llm-prompt.md"
    if not path.exists():
        raise ConfigError(f"Language prompt not found: {path}")
    return path.read_text(encoding="utf-8")


def build_translation_prompt(
    lang: str,
    lang_name: str,
    en_content: str,
    existing: str | None = None,
    en_diff: str | None = None,
    prompt_diff: str | None = None,
) -> str:
    """Build the full prompt for translation.

    Args:
        lang: Language code (e.g., 'pt', 'es')
        lang_name: Human-readable language name (e.g., 'português')
        en_content: Full English source content
        existing: Existing translation content, if updating
        en_diff: Git diff of English changes, for incremental updates
        prompt_diff: Git diff of prompt changes, for targeted prompt updates
    """
    parts = [get_general_prompt(), get_lang_prompt(lang)]

    if existing is not None:
        if prompt_diff:
            # Prompt-change update mode: show what changed in the guidelines
            parts.append(
                "## Prompt Update Mode\n\n"
                "The translation guidelines have been updated. "
                "The following diff shows exactly what changed:\n\n"
                f"```diff\n{prompt_diff}\n```\n\n"
                "**Instructions:**\n\n"
                "1. Review the existing translation below\n"
                "2. Update ONLY the parts affected by these guideline changes\n"
                "3. Preserve ALL other content exactly as-is, character for character\n\n"
                f"## Existing Translation\n\n%%%\n{existing}%%%"
            )
        elif en_diff:
            # Diff-aware incremental update mode
            parts.append(
                "## Incremental Update Mode\n\n"
                "The English source has been updated. "
                "The following diff shows exactly what changed:\n\n"
                f"```diff\n{en_diff}\n```\n\n"
                "**Instructions:**\n\n"
                "1. Locate the corresponding section(s) in your existing translation\n"
                "2. Update ONLY those specific sections to reflect the English changes\n"
                "3. Preserve ALL other content exactly as-is, character for character\n"
                "4. Do NOT rephrase, improve, or modify any sections unrelated to the diff\n\n"
                "**Exception:** If you find major issues beyond the diff (e.g., truncated "
                "content, missing sections, or significant errors), you may fix those as well. "
                "The goal is to avoid unnecessary rephrasing of translations that are already "
                "correct.\n\n"
                f"## Existing Translation\n\n%%%\n{existing}%%%"
            )
        else:
            # Full re-translation with existing as reference
            # (prompt changed, diff too large, or diff unavailable)
            parts.append(
                "## Existing Translation\n"
                "Update minimally: add new content, remove deleted content, "
                "fix guideline violations, preserve correct lines exactly.\n\n"
                f"Previous translation:\n%%%\n{existing}%%%"
            )

    parts.append(
        f"## Task\nTranslate to {lang} ({lang_name}).\n\n"
        f"Original content:\n%%%\n{en_content}%%%"
    )

    return "\n\n".join(parts)


def build_chunk_translation_prompt(
    lang: str,
    lang_name: str,
    en_chunk: str,
    existing_chunk: str | None = None,
    prev_heading: str | None = None,
    next_heading: str | None = None,
) -> str:
    """Build a prompt for translating a single chunk of a document.

    Args:
        lang: Language code (e.g., 'pt', 'es')
        lang_name: Human-readable language name
        en_chunk: The English chunk content to translate
        existing_chunk: Existing translation of this chunk, if updating
        prev_heading: Heading of the previous section (for context)
        next_heading: Heading of the next section (for context)
    """
    parts = [get_general_prompt(), get_lang_prompt(lang)]

    context_lines = []
    if prev_heading:
        context_lines.append(f"Previous section: {prev_heading}")
    if next_heading:
        context_lines.append(f"Next section: {next_heading}")
    context = "\n".join(context_lines) if context_lines else "This is a standalone section."

    if existing_chunk is not None:
        parts.append(
            "## Section Update Mode\n\n"
            "You are updating a single section of a larger translated document.\n"
            "Update the existing translation to match the new English content.\n"
            "Preserve correct translations exactly; only change what differs.\n\n"
            f"### Context\n{context}\n\n"
            f"### Existing Translation of This Section\n%%%\n{existing_chunk}%%%"
        )
    else:
        parts.append(
            "## Section Translation Mode\n\n"
            "You are translating a single section of a larger document.\n\n"
            f"### Context\n{context}"
        )

    parts.append(
        f"## Task\nTranslate this section to {lang} ({lang_name}).\n"
        "Output ONLY the translated section, nothing else.\n\n"
        f"Section to translate:\n%%%\n{en_chunk}%%%"
    )

    return "\n\n".join(parts)
