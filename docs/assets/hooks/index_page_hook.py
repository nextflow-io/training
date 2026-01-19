"""
MkDocs hook for processing index pages with custom frontmatter.

This hook transforms pages with `template: index_page` frontmatter,
generating Material for MkDocs grid cards from structured frontmatter data.
"""

import re

# Default content for various sections
DEFAULT_CONTENT = {
    "technical_requirements": (
        "You will need a GitHub account OR a local installation of Nextflow. "
        "See [Environment options](../envsetup/index.md) for more details."
    ),
    "videos": (
        "Videos are available for each chapter, featuring an instructor working "
        "through the exercises. The video for each part of the course is embedded "
        "at the top of the corresponding page."
    ),
}


def on_page_markdown(markdown, page, config, files):
    """
    Process pages with template: index_page frontmatter.

    Extracts summary content and generates grid cards with additional information
    admonitions from frontmatter data.
    """
    # Skip pages without the index_page page_type
    if page.meta.get("page_type") != "index_page":
        return None

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
    admonitions = generate_admonitions(additional_info, page.file.src_path)

    # Generate badge HTML if index_type is specified
    badge_html = ""
    index_type = page.meta.get("index_type")
    if index_type:
        badge_html = f'<span class="index-type-badge">{index_type}</span>\n\n'

    # Build the grid cards structure
    grid_content = f"""# {h1_title}{badge_html}

<div class="grid cards" markdown>

-   :material-book-open-variant:{{ .lg .middle }} __Course summary__

    ---

{indent_content(summary_content, 4)}

-   :material-information-outline:{{ .lg .middle }} __Additional information__

    ---

{indent_content(admonitions, 4)}

</div>

{rest_content}"""

    return grid_content


def generate_admonitions(additional_info, src_path):
    """Generate the admonition markdown from additional_information frontmatter."""
    if not additional_info:
        return ""

    admonitions = []

    # Technical requirements
    tech_req = additional_info.get("technical_requirements")
    if tech_req:
        if tech_req is True:
            content = DEFAULT_CONTENT["technical_requirements"]
        elif isinstance(tech_req, str):
            content = tech_req
        else:
            raise ValueError(
                f"Page '{src_path}': technical_requirements must be true or a string"
            )
        admonitions.append(
            f'??? terminal "Technical requirements"\n{indent_content(content, 4)}'
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
            f'??? learning "Learning objectives"\n{indent_content(items, 4)}'
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
        admonitions.append(
            f'??? people "Audience & prerequisites"\n{indent_content(items, 4)}'
        )

    # Videos (videos and videos_playlist are mutually exclusive)
    videos = additional_info.get("videos")
    videos_playlist = additional_info.get("videos_playlist")

    if videos and videos_playlist:
        raise ValueError(
            f"Page '{src_path}': 'videos' and 'videos_playlist' are mutually exclusive. "
            "Use 'videos' for custom text or 'videos_playlist' for default text with a playlist link."
        )

    if videos_playlist:
        content = f'{DEFAULT_CONTENT["videos"]}\n\n[View the playlist on YouTube]({videos_playlist}){{:target="_blank"}}'
        admonitions.append(f'??? videos "Course videos"\n{indent_content(content, 4)}')
    elif videos:
        admonitions.append(f'??? videos "Course videos"\n{indent_content(videos, 4)}')

    return "\n\n".join(admonitions)


def indent_content(content, spaces):
    """Indent all lines of content by the specified number of spaces."""
    indent = " " * spaces
    lines = content.split("\n")
    return "\n".join(indent + line if line else line for line in lines)
