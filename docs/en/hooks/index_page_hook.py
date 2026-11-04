"""
MkDocs hook for processing index pages with custom frontmatter.

This hook transforms pages with `template: index_page` frontmatter,
generating Material for MkDocs grid cards from structured frontmatter data.

UI strings are loaded from ui-strings.yml in the language directory,
falling back to English if translations are not available.
"""

import re
from functools import lru_cache
from pathlib import Path

import yaml

# Path to docs directory (parent of hooks directory's parent)
DOCS_PATH = Path(__file__).parent.parent.parent


@lru_cache
def load_ui_strings(lang: str) -> dict:
    """
    Load UI strings for a language, falling back to English.

    Args:
        lang: Language code (e.g., 'en', 'ko', 'de')

    Returns:
        Dictionary with UI strings
    """
    # Try language-specific file first
    lang_file = DOCS_PATH / lang / "ui-strings.yml"
    if lang_file.exists():
        strings = yaml.safe_load(lang_file.read_text(encoding="utf-8"))
        if strings:
            return strings

    # Fall back to English
    en_file = DOCS_PATH / "en" / "ui-strings.yml"
    if en_file.exists():
        return yaml.safe_load(en_file.read_text(encoding="utf-8"))

    # Ultimate fallback - hardcoded English defaults
    return {
        "index_page": {
            "course_summary": "Course summary",
            "additional_information": "Additional information",
            "technical_requirements": "Technical requirements",
            "learning_objectives": "Learning objectives",
            "audience_prerequisites": "Audience & prerequisites",
            "course_videos": "Course videos",
        },
        "defaults": {
            "technical_requirements": (
                "You will need a GitHub account OR a local installation of Nextflow. "
                "See [Environment options](../envsetup/index.md) for more details."
            ),
            "videos": (
                "Videos are available for each chapter, featuring an instructor working "
                "through the exercises. The video for each part of the course is embedded "
                "at the top of the corresponding page."
            ),
        },
    }


def get_lang_from_config(config) -> str:
    """Extract language code from MkDocs config."""
    # Try theme.language first
    lang = config.get("theme", {}).get("language", "en")
    return lang if lang else "en"


def on_page_markdown(markdown, page, config, files):
    """
    Process pages with template: index_page frontmatter.

    Extracts summary content and generates grid cards with additional information
    admonitions from frontmatter data.
    """
    # Skip pages without the index_page page_type
    if page.meta.get("page_type") != "index_page":
        return None

    # Load UI strings for this language
    lang = get_lang_from_config(config)
    ui = load_ui_strings(lang)
    labels = ui.get("index_page", {})
    defaults = ui.get("defaults", {})

    # Extract H1 heading
    h1_match = re.search(r"^# (.+)$", markdown, re.MULTILINE)
    if not h1_match:
        raise ValueError(
            f"Page '{page.file.src_path}' uses template: index_page but has no H1 heading. "
            "Add a line starting with '# ' for the page title."
        )

    h1_title = h1_match.group(1)
    h1_end = h1_match.end()

    # Find the additional_information marker
    marker = "<!-- additional_information -->"
    marker_pos = markdown.find(marker)
    if marker_pos == -1:
        raise ValueError(
            f"Page '{page.file.src_path}' uses template: index_page but is missing "
            f"the '{marker}' marker. Add this marker after your summary content."
        )

    # Extract summary (content between H1 and marker)
    summary_content = markdown[h1_end:marker_pos].strip()

    # Extract content after the marker
    rest_content = markdown[marker_pos + len(marker) :].strip()

    # Generate additional information admonitions
    additional_info = page.meta.get("additional_information", {})
    admonitions = generate_admonitions(
        additional_info, page.file.src_path, labels, defaults
    )

    # Generate badge HTML if index_type is specified
    badge_html = ""
    index_type = page.meta.get("index_type")
    if index_type:
        badge_html = f'<span class="index-type-badge">{index_type}</span>\n\n'

    # Get translated labels
    course_summary = labels.get("course_summary", "Course summary")
    additional_information = labels.get(
        "additional_information", "Additional information"
    )

    # Build the grid cards structure
    grid_content = f"""# {h1_title}{badge_html}

<div class="grid cards" markdown>

-   :material-book-open-variant:{{ .lg .middle }} __{course_summary}__

    ---

{indent_content(summary_content, 4)}

-   :material-information-outline:{{ .lg .middle }} __{additional_information}__

    ---

{indent_content(admonitions, 4)}

</div>

{rest_content}"""

    return grid_content


def generate_admonitions(additional_info, src_path, labels, defaults):
    """Generate the admonition markdown from additional_information frontmatter."""
    if not additional_info:
        return ""

    admonitions = []

    # Get translated labels with English fallbacks
    tech_req_label = labels.get("technical_requirements", "Technical requirements")
    learning_obj_label = labels.get("learning_objectives", "Learning objectives")
    audience_label = labels.get("audience_prerequisites", "Audience & prerequisites")
    videos_label = labels.get("course_videos", "Course videos")

    # Get default content with fallbacks
    default_tech_req = defaults.get(
        "technical_requirements",
        "You will need a GitHub account OR a local installation of Nextflow. "
        "See [Environment options](../envsetup/index.md) for more details.",
    )
    default_videos = defaults.get(
        "videos",
        "Videos are available for each chapter, featuring an instructor working "
        "through the exercises. The video for each part of the course is embedded "
        "at the top of the corresponding page.",
    )

    # Technical requirements
    tech_req = additional_info.get("technical_requirements")
    if tech_req:
        if tech_req is True:
            content = default_tech_req
        elif isinstance(tech_req, str):
            content = tech_req
        else:
            raise ValueError(
                f"Page '{src_path}': technical_requirements must be true or a string"
            )
        admonitions.append(
            f'??? terminal "{tech_req_label}"\n{indent_content(content, 4)}'
        )

    # Learning objectives (must be a list)
    learning_obj = additional_info.get("learning_objectives")
    if learning_obj:
        if learning_obj is True:
            raise ValueError(
                f"Page '{src_path}': learning_objectives must be a list of strings, not true"
            )
        if not isinstance(learning_obj, list):
            raise ValueError(
                f"Page '{src_path}': learning_objectives must be a list of strings"
            )
        items = "\n".join(f"- {item}" for item in learning_obj)
        admonitions.append(
            f'??? learning "{learning_obj_label}"\n{indent_content(items, 4)}'
        )

    # Audience & prerequisites (must be a list)
    audience = additional_info.get("audience_prerequisites")
    if audience:
        if audience is True:
            raise ValueError(
                f"Page '{src_path}': audience_prerequisites must be a list of strings, not true"
            )
        if not isinstance(audience, list):
            raise ValueError(
                f"Page '{src_path}': audience_prerequisites must be a list of strings"
            )
        items = "\n".join(f"- {item}" for item in audience)
        admonitions.append(f'??? people "{audience_label}"\n{indent_content(items, 4)}')

    # Videos (videos and videos_playlist are mutually exclusive)
    videos = additional_info.get("videos")
    videos_playlist = additional_info.get("videos_playlist")

    if videos and videos_playlist:
        raise ValueError(
            f"Page '{src_path}': 'videos' and 'videos_playlist' are mutually exclusive. "
            "Use 'videos' for custom text or 'videos_playlist' for default text with a playlist link."
        )

    if videos_playlist:
        view_playlist = defaults.get("view_playlist", "View the playlist on YouTube")
        content = f'{default_videos}\n\n[{view_playlist}]({videos_playlist}){{:target="_blank"}}'
        admonitions.append(f'??? videos "{videos_label}"\n{indent_content(content, 4)}')
    elif videos:
        admonitions.append(f'??? videos "{videos_label}"\n{indent_content(videos, 4)}')

    return "\n\n".join(admonitions)


def indent_content(content, spaces):
    """Indent all lines of content by the specified number of spaces."""
    indent = " " * spaces
    lines = content.split("\n")
    return "\n".join(indent + line if line else line for line in lines)
