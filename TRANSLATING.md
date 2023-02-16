# 🌐 <abbr title="translation">Translation</abbr> Guide

Table of contents:
- [Contribution model](#contribution-model)
- [Preview](#preview)
    - [Docker](#docker)
    - [Python](#python)
- [Translation Glossaries](#translation-glossaries)

## Contribution model

The typical workflow for contributing with translation is as follows:

1. Make a _fork_ of the GitHub repository to your own account
2. Work locally (see below) and make your changes
    - Check if the language you want to translate to is already enabled in the `mkdocs.yml` file. If it isn't do it according to [this](https://github.com/nextflow-io/training/pull/163/files#diff-98d0f806abc9af24e6a7c545d3d77e8f9ad57643e27211d7a7b896113e420ed2) example.
    - If you want to improve an already existing translation, the file already exists and the language is already set up. Simply open the file and work on it. Otherwise, create the new file with the following pattern: If the original file in English is `filename.md`, you will create in the same folder a new file named `filename.language_code.md`, where `language_code` is the language code that you can find [here](https://en.wikipedia.org/wiki/ISO_639-1), for the language you wish to translate to. Pay attention to the language code, as it has to be the same that is specified in the `mkdocs.yml` file for that language.
3. Commit and push to your forked repository
4. Open a pull-request against the main repo, which can be reviewed and merged

The training website is built using [Material for MkDocs](https://squidfunk.github.io/mkdocs-material/).
It's a static site generator designed for documentation websites which is fast and lightweight and comes with a lot of nice features.
We use the open-source version of the tool (not any of the "insiders" features, currently).

To make changes, you should run the website locally so that you can preview changes.

## Preview

See the [mkdocs material docs](https://squidfunk.github.io/mkdocs-material/getting-started/) for full installation instructions in order to see the preview of the training material with your changes/translations. A short version for this training site is below.

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

Even though it's required that you read this document in order to contribute with translations, it's really important that you also check the [CONTRIBUTING.md](https://github.com/nextflow-io/training/blob/master/CONTRIBUTING.md) file.

## Translation Glossaries
In order to keep consistency in translations, every language should have a translation glossary where common technical terms and terms related to Nextflow have an official translation to be followed by future translators. Ideally, these links should point to an online spreadsheet where anyone can comment and make suggestions, but not edit.

* [Portuguese](https://docs.google.com/spreadsheets/d/1HUa3BO2kwukhX4EXQ-1blXeP5iueUdM23OwDRpfarDg/edit?usp=sharing)

