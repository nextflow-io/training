"""src/redirects.mjs: the mkdocs-redirects `redirect_maps`, for the root and every locale prefix."""

from __future__ import annotations

import json

from core import LOCALES

from . import source


def url(path: str, lang: str) -> str:
    if "://" in path:
        return path
    page, _, anchor = path.partition("#")
    parts = [p for p in ("" if lang == "en" else lang, source.slug(page)) if p]
    return "/" + "".join(f"{p}/" for p in parts) + (f"#{anchor}" if anchor else "")


def generate(run, langs):
    plugins = source.config(run, "en")["plugins"]
    maps = next(p["redirects"]["redirect_maps"] for p in plugins if isinstance(p, dict) and "redirects" in p)
    redirects = {url(old, lang): url(new, lang) for lang in LOCALES if lang in langs for old, new in maps.items()}
    body = json.dumps(redirects, indent=4)[:-1] + "  }"
    return {"src/redirects.mjs": JS.replace("REDIRECTS", body)}


JS = """\
// Astro puts redirect pages under `base` but leaves their destinations as written.
const base = (process.env.BASE ?? '/').replace(/\\/$/, '');

export const redirects = Object.fromEntries(
  Object.entries(REDIRECTS).map(([from, to]) => [from, to.startsWith('/') ? base + to : to]),
);
"""
