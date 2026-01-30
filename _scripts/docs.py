#!/usr/bin/env python3
# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "typer>=0.12.0",
#     "pyyaml>=6.0",
#     "rich>=13.0.0",
#     # MkDocs and plugins
#     "mkdocs",
#     "mkdocs-material",
#     "mkdocs-static-i18n",
#     "pymdown-extensions>=10.12",
#     "pillow",
#     "cairosvg",
#     "mkdocs-enumerate-headings-plugin>=0.6.0",
#     "mkdocs-quiz>=1.5.0",
#     "mike",
# ]
# ///
"""
Build and language management commands for Nextflow training docs.

Usage:
    uv run docs.py build-lang en
    uv run docs.py build-all
    uv run docs.py live pt
"""

import json
import os
import shutil
import subprocess
import sys
from http.server import HTTPServer, SimpleHTTPRequestHandler
from multiprocessing import Pool
from pathlib import Path

import typer
import yaml
from rich.console import Console

app = typer.Typer(help="Nextflow Training Documentation Tools")
console = Console()

# Must be run from repo root
REPO_ROOT = Path(__file__).parent.parent
DOCS_PATH = REPO_ROOT / "docs"
EN_DOCS_PATH = DOCS_PATH / "en"
SITE_PATH = REPO_ROOT / "site"
BUILD_SITE_PATH = REPO_ROOT / "site_build"


def get_supported_langs() -> set[str]:
    """
    Load supported languages from language_names.yml.

    This is the single source of truth for which languages are supported.
    """
    lang_file = DOCS_PATH / "language_names.yml"
    if not lang_file.exists():
        console.print(f"[red]Error:[/red] Language file not found: {lang_file}")
        raise typer.Exit(code=1)
    langs = yaml.safe_load(lang_file.read_text(encoding="utf-8"))
    return set(langs.keys())


def get_lang_paths() -> list[Path]:
    """Get all language directories."""
    return sorted(DOCS_PATH.iterdir())


@app.command()
def new_lang(lang: str):
    """Generate a new docs translation directory for a language."""
    lang = lang.lower()
    new_path = DOCS_PATH / lang

    if new_path.exists():
        console.print(f"[red]Error:[/red] Language already exists: {lang}")
        raise typer.Exit(code=1)

    new_path.mkdir()
    (new_path / "mkdocs.yml").write_text(
        f"INHERIT: ../en/mkdocs.yml\ntheme:\n  language: {lang}\n",
        encoding="utf-8",
    )
    (new_path / "llm-prompt.md").write_text(
        f"# Translation Rules for {lang}\n\nTODO: Add language-specific translation rules here.\n",
        encoding="utf-8",
    )
    (new_path / "docs").mkdir()

    console.print(f"[green]Created:[/green] {new_path}")
    console.print("\nNext steps:")
    console.print(f"  1. Edit {new_path}/llm-prompt.md to add translation rules")
    console.print(f"  2. Add '{lang}: <native name>' to docs/language_names.yml")
    console.print(f"  3. Add language to extra.alternate in docs/en/mkdocs.yml")


@app.command()
def build_lang(lang: str):
    """Build the docs for a language."""
    lang = lang.lower()
    lang_path = DOCS_PATH / lang

    if not lang_path.is_dir():
        console.print(f"[red]Error:[/red] Language directory not found: {lang_path}")
        raise typer.Exit(code=1)

    console.print(f"Building docs for: [cyan]{lang}[/cyan]")
    build_site_dist_path = BUILD_SITE_PATH / lang

    if lang == "en":
        dist_path = SITE_PATH
    else:
        dist_path = SITE_PATH / lang
        shutil.rmtree(dist_path, ignore_errors=True)

    shutil.rmtree(build_site_dist_path, ignore_errors=True)
    subprocess.run(
        ["mkdocs", "build", "--site-dir", str(build_site_dist_path)],
        cwd=lang_path,
        check=True,
    )
    shutil.copytree(build_site_dist_path, dist_path, dirs_exist_ok=True)

    console.print(f"[green]Built:[/green] {lang}")


@app.command()
def build_all():
    """Build docs for all supported languages."""
    shutil.rmtree(SITE_PATH, ignore_errors=True)
    supported = get_supported_langs()
    langs = [p.name for p in get_lang_paths() if p.is_dir() and p.name in supported]

    console.print(f"Building {len(langs)} languages: {', '.join(langs)}")

    cpu_count = os.cpu_count() or 1
    pool_size = min(cpu_count * 2, len(langs))

    with Pool(pool_size) as p:
        p.map(build_lang, langs)

    console.print("[green]All builds complete[/green]")


@app.command()
def serve():
    """Serve the built site (run build-all first)."""
    if not SITE_PATH.exists():
        console.print("[red]Error:[/red] Site not built. Run 'build-all' first.")
        raise typer.Exit(code=1)

    os.chdir(SITE_PATH)
    server = HTTPServer(("", 8008), SimpleHTTPRequestHandler)
    console.print("Serving at: [link]http://127.0.0.1:8008[/link]")
    server.serve_forever()


@app.command()
def live(lang: str = "en", dirty: bool = False):
    """Serve docs with live reload for a specific language."""
    lang = lang.lower()
    lang_path = DOCS_PATH / lang

    if not lang_path.is_dir():
        console.print(f"[red]Error:[/red] Language not found: {lang}")
        raise typer.Exit(code=1)

    args = ["mkdocs", "serve", "--dev-addr", "127.0.0.1:8008"]
    if dirty:
        args.append("--dirty")

    console.print(f"Serving [cyan]{lang}[/cyan] at [link]http://127.0.0.1:8008[/link]")
    subprocess.run(args, cwd=lang_path, check=True)


@app.command()
def langs_json():
    """Output supported languages as JSON (for CI)."""
    supported = get_supported_langs()
    langs = [p.name for p in get_lang_paths() if p.is_dir() and p.name in supported]
    print(json.dumps(sorted(langs)))


@app.command()
def check_config(lang: str):
    """Verify a language's mkdocs.yml inherits correctly."""
    lang = lang.lower()
    config_path = DOCS_PATH / lang / "mkdocs.yml"

    if not config_path.exists():
        console.print(f"[red]Error:[/red] Config not found: {config_path}")
        raise typer.Exit(code=1)

    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))

    if config.get("INHERIT") != "../en/mkdocs.yml":
        console.print(
            f"[red]Error:[/red] {config_path} must have 'INHERIT: ../en/mkdocs.yml'"
        )
        raise typer.Exit(code=1)

    console.print(f"[green]OK:[/green] {config_path}")


if __name__ == "__main__":
    # Ensure we're running from repo root
    if not DOCS_PATH.exists():
        console.print("[red]Error:[/red] Must run from repository root")
        sys.exit(1)
    app()
