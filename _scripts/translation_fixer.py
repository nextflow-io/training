#!/usr/bin/env python3
# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "typer>=0.12.0",
#     "rich>=13.0.0",
# ]
# ///
"""
Post-process translations to fix common LLM mistakes.

Usage:
    uv run translation_fixer.py fix-all pt
    uv run translation_fixer.py fix-pages docs/pt/docs/index.md
"""

import sys
from pathlib import Path

import typer
from rich.console import Console

from doc_parsing_utils import check_translation

app = typer.Typer(help="Fix common issues in translated documents")
console = Console()

# Must be run from repo root
REPO_ROOT = Path(__file__).parent.parent
DOCS_PATH = REPO_ROOT / "docs"
EN_DOCS_PATH = DOCS_PATH / "en" / "docs"


def fix_one_page(path: Path) -> bool:
    """
    Fix a translated document by comparing it to the English version.

    Returns True if successful. Raises on error (fail loud).
    """
    # Determine language from path
    try:
        lang_code = path.relative_to(DOCS_PATH).parts[0]
    except ValueError:
        console.print(f"[red]Error:[/red] Path not in docs/: {path}")
        raise typer.Exit(code=1)

    if lang_code == "en":
        console.print(f"[yellow]Skipping English file:[/yellow] {path}")
        return True

    # Find English source
    rel_path = path.relative_to(DOCS_PATH / lang_code / "docs")
    en_path = EN_DOCS_PATH / rel_path

    if not en_path.exists():
        console.print(f"[red]Error:[/red] English source not found: {en_path}")
        raise typer.Exit(code=1)

    doc_lines = path.read_text(encoding="utf-8").splitlines()
    en_doc_lines = en_path.read_text(encoding="utf-8").splitlines()

    # This will raise ValueError if structure doesn't match - that's intentional
    fixed_lines = check_translation(
        doc_lines=doc_lines,
        en_doc_lines=en_doc_lines,
        auto_fix=True,
        path=str(path),
    )

    # Write back
    fixed_lines.append("")  # Ensure trailing newline
    path.write_text("\n".join(fixed_lines), encoding="utf-8")
    return True


@app.command()
def fix_all(language: str):
    """Fix all translation files for a language."""
    lang_docs = DOCS_PATH / language / "docs"

    if not lang_docs.exists():
        console.print(f"[red]Error:[/red] Language docs not found: {lang_docs}")
        raise typer.Exit(code=1)

    files = list(lang_docs.rglob("*.md"))
    console.print(f"Fixing {len(files)} files for {language}...")

    for path in files:
        console.print(f"  {path.name}", end=" ")
        fix_one_page(path)
        console.print("[green]OK[/green]")

    console.print(f"\n[green]Fixed {len(files)} files[/green]")


@app.command()
def fix_pages(paths: list[Path]):
    """Fix specific translation files."""
    for path in paths:
        if not path.exists():
            console.print(f"[red]Error:[/red] File not found: {path}")
            raise typer.Exit(code=1)

        console.print(f"Fixing {path}...", end=" ")
        fix_one_page(path)
        console.print("[green]OK[/green]")


if __name__ == "__main__":
    if not DOCS_PATH.exists():
        console.print("[red]Error:[/red] Must run from repository root")
        sys.exit(1)
    app()
