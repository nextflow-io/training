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
#     "pymdown-extensions>=10.12",
#     "mkdocs-enumerate-headings-plugin>=0.6.0",
#     "mkdocs-quiz>=1.5.4",
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
    return set(get_language_names().keys())


def get_language_names() -> dict[str, str]:
    """
    Load language code to native name mapping from language_names.yml.

    This is the single source of truth for which languages are supported.
    """
    lang_file = DOCS_PATH / "language_names.yml"
    if not lang_file.exists():
        console.print(f"[red]Error:[/red] Language file not found: {lang_file}")
        raise typer.Exit(code=1)
    return yaml.safe_load(lang_file.read_text(encoding="utf-8"))


def get_lang_paths() -> list[Path]:
    """Get all language directories."""
    return sorted(DOCS_PATH.iterdir())


MKDOCS_TEMPLATE = """\
INHERIT: ../en/mkdocs.yml
theme:
  language: {lang}
  custom_dir: ../en/overrides
extra:
  consent:
    title: "Cookie consent"
    description: >-
      We use cookies to recognize your repeated visits and preferences, as well
      as to measure the effectiveness of our documentation and whether users
      find what they're searching for. With your consent, you're helping us to
      make our training materials better.
      Find out more on
      <a href="https://seqera.io/privacy-policy/#cookies" target="_blank" rel="noopener">how we use cookies</a>.
    cookies:
      posthog: "PostHog Analytics"
"""


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
        MKDOCS_TEMPLATE.format(lang=lang),
        encoding="utf-8",
    )
    (new_path / "llm-prompt.md").write_text(
        f"# Translation Rules for {lang}\n\nTODO: Add language-specific translation rules here.\n",
        encoding="utf-8",
    )

    # Copy ui-strings.yml from English as a starting point
    en_ui_strings = EN_DOCS_PATH / "ui-strings.yml"
    if en_ui_strings.exists():
        ui_strings_content = en_ui_strings.read_text(encoding="utf-8")
        # Add TODO comment at the top
        ui_strings_content = (
            f"# TODO: Translate these UI strings to {lang}\n" + ui_strings_content
        )
        (new_path / "ui-strings.yml").write_text(ui_strings_content, encoding="utf-8")

    (new_path / "docs").mkdir()

    console.print(f"[green]Created:[/green] {new_path}")
    console.print("\nNext steps:")
    console.print(f"  1. Edit {new_path}/llm-prompt.md to add translation rules")
    console.print(f"  2. Translate {new_path}/ui-strings.yml")
    console.print(f"  3. Translate extra.consent in {new_path}/mkdocs.yml")
    console.print(f"  4. Add '{lang}: <native name>' to docs/language_names.yml")
    console.print(
        "  5. Run 'uv run docs.py sync-language-picker' to update language picker"
    )


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
    result = subprocess.run(
        ["mkdocs", "build", "--site-dir", str(build_site_dist_path)],
        cwd=lang_path,
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        console.print(f"[red]Error building {lang}:[/red]")
        console.print(result.stdout)
        console.print(result.stderr)
        raise subprocess.CalledProcessError(result.returncode, result.args)
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


@app.command()
def sync_language_picker():
    """Update the language picker in docs/en/mkdocs.yml from language_names.yml."""
    import re

    mkdocs_path = EN_DOCS_PATH / "mkdocs.yml"
    if not mkdocs_path.exists():
        console.print(f"[red]Error:[/red] mkdocs.yml not found: {mkdocs_path}")
        raise typer.Exit(code=1)

    lang_names = get_language_names()
    content = mkdocs_path.read_text(encoding="utf-8")

    # Build the new alternate section
    # English first (primary language), then others alphabetically
    alternate_lines = ["  alternate:"]
    lang_codes = ["en"] + sorted(k for k in lang_names.keys() if k != "en")
    for lang_code in lang_codes:
        name = lang_names[lang_code]
        # English goes to root, others to /{lang}/
        link = "/" if lang_code == "en" else f"/{lang_code}/"
        alternate_lines.append(f"    - name: {name}")
        alternate_lines.append(f"      link: {link}")
        alternate_lines.append(f"      lang: {lang_code}")
    new_alternate = "\n".join(alternate_lines)

    # Pattern to match the existing alternate section
    # Matches from "  alternate:" to the next unindented or less-indented line
    pattern = r"(  alternate:\n(?:    .*\n)+)"

    if not re.search(pattern, content):
        console.print(
            "[red]Error:[/red] Could not find alternate section in mkdocs.yml"
        )
        raise typer.Exit(code=1)

    new_content = re.sub(pattern, new_alternate + "\n", content)

    if new_content == content:
        console.print("[green]Language picker already up to date[/green]")
        return

    mkdocs_path.write_text(new_content, encoding="utf-8")
    console.print(f"[green]Updated:[/green] {mkdocs_path}")
    console.print(f"Languages: {', '.join(sorted(lang_names.keys()))}")


if __name__ == "__main__":
    # Ensure we're running from repo root
    if not DOCS_PATH.exists():
        console.print("[red]Error:[/red] Must run from repository root")
        sys.exit(1)
    app()
