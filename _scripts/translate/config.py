"""Configuration constants, exceptions, and shared utilities."""

from __future__ import annotations

import os
from pathlib import Path

from rich.console import Console

# Paths (resolved relative to this file: _scripts/translate/config.py)
REPO_ROOT = Path(__file__).parent.parent.parent
DOCS_ROOT = REPO_ROOT / "docs"
EN_DOCS = DOCS_ROOT / "en" / "docs"
SCRIPTS_DIR = REPO_ROOT / "_scripts"

# Claude API settings
MODEL = "claude-sonnet-4-5"
MAX_TOKENS = 32768  # Large enough for biggest docs (~60KB source)
REQUEST_TIMEOUT = 600.0  # 10 minute timeout per API request

# Retry configuration for transient API errors
MAX_RETRIES = 8
BASE_DELAY = 2.0  # seconds, doubles each retry: 2, 4, 8, 16, 32, 64, 128, 256

# Translation settings
MAX_CONTINUATIONS = 5  # Max continuation requests for very large files
MAX_VERIFY_RETRIES = 2  # Re-translation attempts after verification failure
DEFAULT_PARALLEL = 50
PRIORITY_DIRS = ["hello_nextflow", "hello_nf-core", "nf4_science", "envsetup"]

# Comment styles by language (used for code block post-processing)
HASH_COMMENT_LANGS = {"python", "py", "sh", "bash", "dockerfile", "yaml", "yml", "toml"}
SLASH_COMMENT_LANGS = {"console"}  # Note: json has no comments
MIXED_COMMENT_LANGS = {"nextflow", "groovy", "nf", "java", "kotlin"}  # //, # and /* */


# ---------------------------------------------------------------------------
# Exceptions
# ---------------------------------------------------------------------------


class TranslationError(Exception):
    """Base exception for translation errors."""


class ConfigError(TranslationError):
    """Configuration or setup error."""


class StructureMismatchError(TranslationError):
    """Translation structure doesn't match source."""


# ---------------------------------------------------------------------------
# Console factory
# ---------------------------------------------------------------------------


def make_console(*, stderr: bool = False) -> Console:
    """Create a Rich console with CI-friendly settings.

    Centralizes the force_terminal logic so it isn't repeated everywhere.
    """
    return Console(
        stderr=stderr,
        force_terminal=True if os.getenv("GITHUB_ACTIONS") else None,
    )


def check_api_key() -> None:
    """Verify ANTHROPIC_API_KEY is set."""
    if not os.getenv("ANTHROPIC_API_KEY"):
        raise ConfigError(
            "ANTHROPIC_API_KEY not set. Get key at https://console.anthropic.com/settings/keys"
        )
