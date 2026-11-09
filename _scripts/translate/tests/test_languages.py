"""Tests for page paths, the translation baseline and adding a language."""

from __future__ import annotations

import json
from types import SimpleNamespace

import github
import pytest
import yaml
from typer.testing import CliRunner

from translate import cli, git_utils, paths
from translate.config import DOCS_ROOT, EN_DOCS


def test_page_paths():
    en = EN_DOCS / "hello_nextflow" / "01_hello_world.mdx"
    de = DOCS_ROOT / "de" / "hello_nextflow" / "01_hello_world.mdx"
    assert paths.en_to_lang_path(en, "de") == de
    assert paths.lang_to_en_path(de, "de") == en
    assert paths.is_en_doc(en)
    assert not paths.is_en_doc(de)
    assert not paths.is_en_doc(EN_DOCS / "nf4_science" / "_template" / "_README.mdx")


def test_iter_en_docs_skips_locales():
    docs = paths.iter_en_docs()
    assert docs and all(paths.is_en_doc(p) for p in docs)
    assert len(docs) == len(set(docs))


@pytest.mark.parametrize(
    "pr_commit, floor, ancestor, expected",
    [
        ("aaaaaaa", None, None, "aaaaaaa"),  # no floor: last translation PR
        ("aaaaaaa", "bbbbbbb", True, "bbbbbbb"),  # PR older than the floor
        ("ccccccc", "bbbbbbb", False, "ccccccc"),  # a newer translation PR wins
        (None, "bbbbbbb", None, "bbbbbbb"),  # no translation PR yet
    ],
)
def test_baseline_floor(monkeypatch, pr_commit, floor, ancestor, expected):
    prs = [SimpleNamespace(merged=True, title="Update translations", body=f"commit {pr_commit}")]
    repo = SimpleNamespace(get_pulls=lambda **_: prs if pr_commit else [])
    monkeypatch.setattr(github, "Github", lambda **_: SimpleNamespace(get_repo=lambda _: repo))
    monkeypatch.setattr(git_utils, "get_translate_config", lambda: {"baseline": floor})
    monkeypatch.setattr(
        git_utils, "_get_repo", lambda: SimpleNamespace(is_ancestor=lambda a, b: ancestor)
    )
    git_utils.get_translation_baseline.cache_clear()
    try:
        assert git_utils.get_translation_baseline("token", "org/repo") == expected
    finally:
        git_utils.get_translation_baseline.cache_clear()


def test_new_lang(monkeypatch, tmp_path):
    names = tmp_path / "language_names.yml"
    names.write_text("de: Deutsch\nen: English\n")
    i18n = tmp_path / "i18n"
    i18n.mkdir()
    (i18n / "en.json").write_text('{"a": "b"}')
    monkeypatch.setattr(paths, "LANGUAGE_NAMES", names)
    monkeypatch.setattr(cli, "LANGUAGE_NAMES", names)
    monkeypatch.setattr(cli, "LANGUAGES_DIR", tmp_path / "languages")
    monkeypatch.setattr(cli, "I18N_DIR", i18n)
    monkeypatch.setattr(cli, "REPO_ROOT", tmp_path)
    paths.get_languages.cache_clear()
    try:
        result = CliRunner().invoke(cli.app, ["new-lang", "JA", "--name", "日本語"])
        assert result.exit_code == 0, result.output
        assert "docs/src/locales.mjs" in result.output
        assert yaml.safe_load(names.read_text()) == {"de": "Deutsch", "en": "English", "ja": "日本語"}
        assert names.read_text() == "de: Deutsch\nen: English\nja: 日本語\n"
        assert json.loads((i18n / "ja.json").read_text()) == {"a": "b"}
        glossary = yaml.safe_load((tmp_path / "languages" / "ja" / "glossary.yml").read_text())
        assert "TRANSLATING.md" in glossary["translation_notice"]["text"]
        assert (tmp_path / "languages" / "ja" / "llm-prompt.md").exists()

        again = CliRunner().invoke(cli.app, ["new-lang", "ja", "--name", "x"])
        assert again.exit_code == 1
    finally:
        paths.get_languages.cache_clear()
