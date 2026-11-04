"""src/sidebar.mjs: the Starlight sidebar from the MkDocs nav, with labels for every locale.

- `nav_group_labels` become the top-level groups, each starting at the section it labels.
- `nav_section_separators` become nested groups, each starting at the page it precedes.
- A section is labelled by `nav_title_overrides`, else by its index page title.
- Labels from config go in `translations`. Plain pages are bare slugs: Starlight labels them with the locale's page title.
"""

from __future__ import annotations

import json

from core import LOCALES

from . import source


def slug(path: str) -> str:
    return source.slug(path) or "index"


def nav_path(item) -> str | None:
    if isinstance(item, str):
        return item
    value = next(iter(item.values()))
    return value if isinstance(value, str) else None


def generate(run, langs):
    langs = [l for l in LOCALES if l == "en" or l in langs]
    extras = {lang: source.extra(run, lang) for lang in langs}

    def setting(lang: str, name: str, key: str | None) -> str | None:
        return extras[lang].get(name, {}).get(key)

    def labelled(by_lang: dict, **rest) -> dict:
        en = by_lang["en"]
        tr = {l: v for l, v in by_lang.items() if l != "en" and v and v != en}
        return {"label": en, **({"translations": tr} if tr else {}), **rest}

    def section_title(lang: str, key: str, children: list) -> str:
        first = nav_path(children[0]) if children else None
        if not (first and first.endswith("index.md")):
            return key
        return setting(lang, "nav_title_overrides", first) or source.page_title(run, lang, first) or key

    def entry(item):
        if isinstance(item, str):
            return slug(item)
        ((title, value),) = item.items()
        if isinstance(value, str):
            return labelled(
                {l: setting(l, "nav_child_title_overrides", value) or title for l in langs}, slug=slug(value)
            )
        titles = {l: section_title(l, title, value) for l in langs}
        return labelled(titles, collapsed=True, items=grouped(value, separator))

    def separator(item) -> dict | None:
        path = nav_path(item)
        if setting("en", "nav_section_separators", path):
            return {l: setting(l, "nav_section_separators", path) for l in langs}

    def group_label(item) -> dict | None:
        if isinstance(item, str) or nav_path(item):
            return None
        ((key, children),) = item.items()
        if setting("en", "nav_group_labels", section_title("en", key, children)):
            return {l: setting(l, "nav_group_labels", section_title(l, key, children)) for l in langs}

    def grouped(items: list, starts_group) -> list:
        """Entries in nav order. An item that starts a group goes in it, with the items after it."""
        out = target = []
        for item in items:
            if labels := starts_group(item):
                out.append(labelled(labels, items=[]))
                target = out[-1]["items"]
            target.append(entry(item))
        return out

    sidebar = grouped(source.config(run, "en")["nav"], group_label)
    return {"src/sidebar.mjs": "export const sidebar = " + json.dumps(sidebar, indent=2, ensure_ascii=False) + ";\n"}
