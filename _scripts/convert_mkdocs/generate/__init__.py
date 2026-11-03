"""Generators for config the converter derives from MkDocs (stream A9).

Each module has `generate(run: Run, langs: list[str]) -> dict[str, str]`, mapping a
path relative to the Astro app (docs/) to file content. The driver writes them.
"""

from . import i18n, redirects, sidebar

GENERATORS = [sidebar, redirects, i18n]
