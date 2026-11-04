"""Read MkDocs config and page titles from the source checkout."""

from __future__ import annotations

import re
from pathlib import Path

import yaml

from core import fenced_lines, split_frontmatter
from rules.frontmatter import H1, plain


class _Loader(yaml.SafeLoader):
    pass


# mkdocs.yml names Python objects for Markdown extensions; the generators never read them.
_Loader.add_multi_constructor("tag:yaml.org,2002:python/", lambda loader, suffix, node: None)


def load_yaml(path: Path) -> dict:
    return yaml.load(path.read_text(encoding="utf-8"), Loader=_Loader) or {}


def config(run, lang: str) -> dict:
    """The locale's own mkdocs.yml, without what it INHERITs from English."""
    return load_yaml(run.src_root / "docs" / lang / "mkdocs.yml")


def extra(run, lang: str) -> dict:
    """`extra` as MkDocs saw it: INHERIT deep-merges the locale's dicts over English."""
    out = dict(config(run, "en").get("extra") or {})
    if lang != "en":
        for k, v in (config(run, lang).get("extra") or {}).items():
            out[k] = {**out[k], **v} if isinstance(v, dict) and isinstance(out.get(k), dict) else v
    return out


def slug(path: str) -> str:
    """Starlight's slug for a docs path: lowercase, no extension, `index` is its directory."""
    return re.sub(r"(^|/)index$", "", path.removesuffix(".md").lower())


def page_title(run, lang: str, rel: str) -> str | None:
    """The page's MkDocs title (frontmatter `title`, else the first H1), from English if untranslated."""
    path = run.docs_dir(lang) / rel
    if not path.is_file():
        path = run.docs_dir("en") / rel
    if not path.is_file():
        return None
    meta, body = split_frontmatter(path.read_text(encoding="utf-8"))
    if meta.get("title"):
        return str(meta["title"])
    lines = body.split("\n")
    fenced = fenced_lines(lines)
    return next((plain(m["text"]) for i, l in enumerate(lines) if i not in fenced and (m := H1.match(l))), None)
