"""src/content/i18n/<lang>.json from each locale's ui-strings.yml and its own `extra.consent`.

Only keys the locale sets are written: Starlight falls back to English for the rest.
"""

from __future__ import annotations

import json

from core import LOCALES

from . import source

# i18n key -> path in ui-strings.yml
UI_STRINGS = {
    "courseInfo.courseSummary": ("index_page", "course_summary"),
    "courseInfo.additionalInformation": ("index_page", "additional_information"),
    "courseInfo.technicalRequirements": ("index_page", "technical_requirements"),
    "courseInfo.learningObjectives": ("index_page", "learning_objectives"),
    "courseInfo.audiencePrerequisites": ("index_page", "audience_prerequisites"),
    "courseInfo.courseVideos": ("index_page", "course_videos"),
    "courseInfo.defaultTechnicalRequirements": ("defaults", "technical_requirements"),
    "courseInfo.defaultVideos": ("defaults", "videos"),
    "courseInfo.viewPlaylist": ("defaults", "view_playlist"),
}
# i18n key -> path in mkdocs.yml `extra.consent`
CONSENT = {
    "consent.title": ("title",),
    "consent.description": ("description",),
    "consent.posthog": ("cookies", "posthog"),
}


def lookup(data: dict, path: tuple) -> str | None:
    for key in path:
        data = data.get(key) if isinstance(data, dict) else None
    return data


def generate(run, langs):
    out = {}
    for lang in (l for l in LOCALES if l in langs):
        ui_path = run.src_root / "docs" / lang / "ui-strings.yml"
        ui = source.load_yaml(ui_path) if ui_path.is_file() else {}
        consent = (source.config(run, lang).get("extra") or {}).get("consent") or {}
        strings = {k: lookup(ui, p) for k, p in UI_STRINGS.items()}
        strings |= {k: lookup(consent, p) for k, p in CONSENT.items()}
        strings = {k: v for k, v in strings.items() if v}
        if strings:
            out[f"src/content/i18n/{lang}.json"] = json.dumps(strings, indent=2, ensure_ascii=False) + "\n"
    return out
