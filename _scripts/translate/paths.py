"""Path utilities and language discovery."""

from __future__ import annotations

from functools import lru_cache
from pathlib import Path

import yaml

from .config import DOCS_ROOT, EN_DOCS, PRIORITY_DIRS, ConfigError


@lru_cache
def get_languages() -> dict[str, str]:
    """Load language code -> name mapping from language_names.yml."""
    path = DOCS_ROOT / "language_names.yml"
    if not path.exists():
        raise ConfigError(f"Language names file not found: {path}")
    data = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise ConfigError(f"Expected dict in {path}, got {type(data).__name__}")
    return data


def get_translation_languages() -> list[str]:
    """Get all translation language codes (excludes English)."""
    return sorted(
        d.name
        for d in DOCS_ROOT.iterdir()
        if d.is_dir() and (d / "mkdocs.yml").exists() and d.name != "en"
    )


def en_to_lang_path(en_path: Path, lang: str) -> Path:
    """Convert English doc path to equivalent path in target language."""
    return DOCS_ROOT / lang / "docs" / en_path.relative_to(EN_DOCS)


def lang_to_en_path(lang_path: Path, lang: str) -> Path:
    """Convert language doc path to equivalent English path."""
    return EN_DOCS / lang_path.relative_to(DOCS_ROOT / lang / "docs")


def iter_en_docs() -> list[Path]:
    """List all English docs in priority order.

    Priority directories are listed first, then remaining files.
    Within each group, files are sorted alphabetically.
    """
    paths: list[Path] = []
    seen: set[Path] = set()

    def add(p: Path) -> None:
        if p not in seen:
            paths.append(p)
            seen.add(p)

    # Root files first
    for p in sorted(EN_DOCS.glob("*.md")):
        add(p)

    # Priority directories
    for dir_name in PRIORITY_DIRS:
        dir_path = EN_DOCS / dir_name
        if dir_path.exists():
            for p in sorted(dir_path.rglob("*.md")):
                add(p)

    # Remaining
    for p in sorted(EN_DOCS.rglob("*.md")):
        add(p)

    return paths
