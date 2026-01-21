# Contributing

Table of contents:

- [Contributing](#contributing)
  - [Contribution model](#contribution-model)
  - [Installation](#installation)
    - [Docker](#docker)
    - [Python](#python)
    - [Social cards](#social-cards)
  - [Figures \& diagrams](#figures--diagrams)
  - [Content style and formatting](#content-style-and-formatting)
    - [Formatting / linting](#formatting--linting)
    - [Nextflow linting](#nextflow-linting)
    - [Headings CI tests](#headings-ci-tests)
    - [Admonitions](#admonitions)
    - [Index page template](#index-page-template)
  - [Known limitations](#known-limitations)
    - [Code annotations](#code-annotations)
    - [Word highlighting](#word-highlighting)
  - [TODO / FIXME](#todo--fixme)

## Contribution model

The typical workflow for contributing to training content is as follows:

1. Make a _fork_ of the GitHub repository to your own account
2. Develop locally (see below) and make your changes
3. Commit and push to your forked repository
4. Open a pull-request against the main repo, which can be reviewed and merged

The training website is built using [Material for MkDocs](https://squidfunk.github.io/mkdocs-material/).
It's a static site generator designed for documentation websites which is fast and lightweight and comes with a lot of nice features.
We use the open-source version of the tool (not any of the "insiders" features, currently).

To make changes, you should run the website locally so that you can preview changes.

## Installation

See the [mkdocs material docs](https://squidfunk.github.io/mkdocs-material/getting-started/) for full installation instructions.
A short version for this training site is below.

### Docker

If you are used to using Docker and don't want to mess around with Python, you can run the following command to preview the site:

```bash
docker run --rm -it -p 8000:8000 -v ${PWD}:/docs ghcr.io/nextflow-io/training-mkdocs:latest
```

This uses a custom image with all required mkdocs plugins.
You should get some output that looks like this:

```console
INFO     -  Documentation built in 27.56 seconds
INFO     -  [21:52:17] Watching paths for changes: 'docs', 'mkdocs.yml'
INFO     -  [21:52:17] Serving on http://0.0.0.0:8000/
```

Visit <http://0.0.0.0:8000/> in your web browser to view the site.
Pages will automatically refresh when you save changes in your editor.

### Python

If you have a recent version of Python installed, then local installation should be as simple as:

```bash
pip install -r requirements.txt
```

Once installed, you can view the site by running:

```bash
mkdocs serve
```

The log output will show a URL, probably <http://127.0.0.1:8000/> - open this in your browser to view the site.
Pages will automatically refresh when you save changes in your editor.

### Social cards

If you're having trouble with the social sharing card images, set the environment variable `CARDS` to `false`:

```bash
CARDS=false mkdocs serve
```

```bash
docker run --rm -it -p 8000:8000 -e 'CARDS=false' -v ${PWD}:/docs ghcr.io/nextflow-io/training-mkdocs:latest
```

## Announcement banner

If there is an announcement banner, you can enable and customise it using the following config in `mkdocs.yml`:

```yaml
extra:
  # Announcement banner for upcoming training
  announcement:
    active: false
    date_text: March 5-6, 2024
    register_url: https://nf-co.re/events/2024/training-foundational-march
```

If you need more customisation, edit `docs/assets/overrides/main.html`

## Figures & diagrams

Graphics should be drawn using [Excalidraw](https://excalidraw.com/).
Please use the [VSCode extension](https://marketplace.visualstudio.com/items?itemName=pomdtr.excalidraw-editor) and give files a `.excalidraw.svg` filename suffix.
Files will continue to be editable by others using this method.

Excalidraw SVGs should be embedded as follows:

<!-- prettier-ignore-start -->
```html
<figure class="excalidraw">
--8<-- "docs/basic_training/img/channel-files.excalidraw.svg"
</figure>
```
<!-- prettier-ignore-end -->

> Note: The file path is from the root of the repo, not the markdown file!

We inline the SVG into the content like this to make remotely loaded fonts work, as well as dark-mode compatibility.

## Content style and formatting

All training content must be written as markdown.

### Formatting / linting

Please make sure that you have Prettier installed and working locally: <https://prettier.io/> (ideally via the VSCode plugin or similar, formatting on save).

There is a GitHub action that checks pull-requests for valid formatting.

For paragraphs containing multiple sentences, please start each sentence on a new line.
Doing that produces cleaner diffs for future code reviews.
The sentences will still be rendered as a continuous paragraph.
For an example, have a look at the code for this paragraph.

### Nextflow linting

This repository uses `nextflow lint` to check Nextflow scripts for errors and warnings.
The linter runs automatically in CI and posts results as a PR comment.

To run locally:

```bash
nextflow lint .
```

### Headings CI tests

This repository includes a Python tool to validate markdown heading numbering consistency across training materials.

The `check_headings.py` script ensures:

- Sequential numbering at each level (1., 1.1., 1.2., etc.)
- Trailing periods after heading numbers
- Heading levels match numbering depth (## for 1., ### for 1.1.)

The easiest way to run it is [with `uv`](https://docs.astral.sh/uv/), which handles dependencies for you automatically:

```bash
# Check files for issues
uv run .github/check_headings.py docs/**/*.md
```

```bash
# Auto-fix detected issues
uv run .github/check_headings.py --fix docs/**/*.md
```

Otherwise, run `pip install typer rich` then `python .github/check_headings.py`.

The script runs automatically in CI on markdown file changes via GitHub Actions,
and will cause a CI failure if any incorrect headings are found.

### Admonitions

We use admonitions extensively to make certain pieces of content stand out.
Please see the [official docs](https://squidfunk.github.io/mkdocs-material/reference/admonitions/) for an explanation.

- Note that we have two custom admonitions: `exercise` and `result` (alias `solution`).
- `!!!` does a regular admonition, `???` makes it collapsed (click to expand).
- Indentation is important! Make sure you check the rendered site, as it's easy to make a mistake.

### Index page template

Course and module index pages can use a special template system that generates Material for MkDocs grid cards from structured frontmatter data.
This keeps the index pages consistent and makes them easier to maintain.

To use the template, add `page_type: index_page` to your frontmatter and include the `<!-- additional_information -->` marker in your content.

#### Basic structure

```markdown
---
title: Course Title
hide:
  - toc
page_type: index_page
index_type: course
additional_information:
  technical_requirements: true
  learning_objectives:
    - First objective
    - Second objective
  audience_prerequisites:
    - "**Audience:** Description of target audience"
    - "**Skills:** Required skills"
  videos_playlist: https://www.youtube.com/playlist?list=...
---

# Course Title

Summary paragraph describing the course.
This content appears in the "Course summary" card.

<!-- additional_information -->

## Rest of page

Content after the marker appears below the grid cards.
```

#### Frontmatter fields

| Field                    | Type           | Required | Description                                                   |
| ------------------------ | -------------- | -------- | ------------------------------------------------------------- |
| `page_type`              | `"index_page"` | Yes      | Enables the index page template                               |
| `index_type`             | string         | No       | Badge label displayed in top-right (e.g., "course", "module") |
| `additional_information` | object         | No       | Container for the collapsible admonitions                     |

#### Additional information fields

All fields within `additional_information` are optional:

| Field                    | Type             | Description                                                                                       |
| ------------------------ | ---------------- | ------------------------------------------------------------------------------------------------- |
| `technical_requirements` | `true` or string | If `true`, uses default text about GitHub/local installation. If a string, uses that custom text. |
| `learning_objectives`    | list of strings  | Rendered as a bulleted list. Must be a list, not `true`.                                          |
| `audience_prerequisites` | list of strings  | Rendered as a bulleted list. Supports markdown formatting. Must be a list, not `true`.            |
| `videos_playlist`        | URL string       | Uses default video description text plus a link to the playlist.                                  |
| `videos`                 | string           | Custom video description text (no link). Mutually exclusive with `videos_playlist`.               |

#### Default content

When `technical_requirements: true` is set:

> You will need a GitHub account OR a local installation of Nextflow. See [Environment options](../envsetup/index.md) for more details.

When `videos_playlist` is set, the following text precedes the link:

> Videos are available for each chapter, featuring an instructor working through the exercises. The video for each part of the course is embedded at the top of the corresponding page.

#### Requirements

- The page must have an H1 heading (`# Title`)
- The page must include the `<!-- additional_information -->` marker
- `learning_objectives` and `audience_prerequisites` must be lists (not `true`)
- `videos` and `videos_playlist` are mutually exclusive

The build will fail with a descriptive error if these requirements are not met.

## Known limitations

There are a couple of known limitations that I haven't figured out how to get around yet

### Code annotations

Mkdocs Material uses `// code comments` to anchor the annotations. That's great, until you want an annotation in the middle of a large multi-line string (say, like a `script` block).

There will hopefully be a way to add annotations at arbitrary line numbers in the future.
See the footnote on the mkdocs material docs:

> Code annotations [are] currently not compatible with [..] languages that do not have comments in their grammar. However, we're actively working on supporting alternate ways of defining code annotations, allowing to always place code annotations at the end of lines.

See [this GitHub discussions thread](https://github.com/squidfunk/mkdocs-material/discussions/3832#discussioncomment-4871068) for updates.

Note also that annotations cannot be added to script blocks generated by importing an external file.

### Word highlighting

Code blocks can have lines highlighted with `hl_lines` in the code block header. However, specific words / characters can not have additional focus (as in a GitHub diff, for example).

## TODO / FIXME

Please use `TODO` comments when you see something that needs coming back to.
I recommend the [Todo Tree VSCode extension](https://marketplace.visualstudio.com/items?itemName=Gruntfuggly.todo-tree) to find these comments easily.

A list of key ones also included here:

- Remove plugin install from Phil's GitHub fork in `requirements.txt` and `.github/mkdocs.Dockerfile` when [this PR](https://github.com/timvink/mkdocs-enumerate-headings-plugin/pull/33) is merged
