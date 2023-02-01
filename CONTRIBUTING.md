# Contributing

## Contribution model

The typical workflow for contributing to training content is as follows:

1. Make a _fork_ of the GitHub repository to your own account
2. Develop locally (see below) and make your changes
3. Commit and push to your forked repository
4. Open a pull-request against the main repo, which can be reviewed and merged

## Material for MkDocs

The training website is built using [Material for MkDocs](https://squidfunk.github.io/mkdocs-material/).
It's a static site generator designed for documentation websites which is fast and lightweight and comes with a lot of nice features.
We use the open-source version of the tool (not any of the "insiders" features, currently).

To make changes, you should run the website locally so that you can preview changes.

### Installation

If you have a recent version of Python installed, then local installation should be as simple as:

```bash
pip install -r requirements.txt
```

See the main mkdocs material docs for installation instructions: <https://squidfunk.github.io/mkdocs-material/getting-started/>

Note that we use the non-standard [mkdocs-enumerate-headings-plugin](https://github.com/timvink/mkdocs-enumerate-headings-plugin) plugin, which also needs to be installed.

### Social cards

If you're having trouble with the social sharing card images, set the environment variable `CARDS` to `false`:

```bash
CARDS=false mkdocs serve
```

or

```bash
export CARDS=false
mkdocs serve
```

## Content style and formatting

All training content must be written as markdown.

### Formatting / linting

Please make sure that you have Prettier installed and working locally: <https://prettier.io/> (ideally via the VSCode plugin or similar, formatting on save).

There is a GitHub action that checks pull-requests for valid formatting.

### Admonitions

We use admonitions extensively to make certain pieces of content stand out.
Please see the [official docs](https://squidfunk.github.io/mkdocs-material/reference/admonitions/) for an explanation.

-   Note that we have two custom admonitions: `exercise` and `result` (alias `solution`).
-   `!!!` does a regular admonition, `???` makes it collapsed (click to expand).
-   Intendation is important! Make sure you check the rendered site, as it's easy to make a mistake.
