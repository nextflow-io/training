---
title: Nextflow versions
description:
  Understanding and managing the evolution of Nextflow's syntax versions
  - toc
  - footer
---

# Nextflow versions

TODO [Blurb about Nextflow versions disambiguating DSL1 / DSL2 and strict syntax v2 parser stuff.]

TODO [tip about trainings being versioned and how to switch between them?]

## Current supported version & requirements

TODO [Latest materials are based on Nextflow 25.10.2 and use typed inputs at the workflow level, which requires use of the V2 syntax parser] [A few exceptions: nf-core courses, any course that runs an nf-core pipeline, and deprecated materials. Noted in course instructions.]

If you plan to use the environment we provide through [Github Codespaces](envsetup/01_setup.md) or [local devcontainers](envsetup/03_devcontainer.md), you don't need to do anything unless specifically noted in the course instructions.

However, if you are planning to work through the trainings in your own environment ([Manual install](envsetup/02_local.md)), you will need to make sure to use Nextflow version 25.10.2 or later with the v2 syntax parser enabled.

### Checking and setting the Nextflow version

You can check what version of Nextflow is installed on your system using the command `nextflow --version`.

TODO [refer to Nextflow docs for version management instructions]

### Enabling the v2 syntax parser

To **enable** the v2 syntax parser for your current session, run the following command in your terminal:

```bash
export NXF_SYNTAX_PARSER=v2
```

To make this permanent (pending v2 becoming the default in Nextflow 26.04), add the export command to your shell profile (`~/.bashrc`, `~/.zshrc`, etc.):

```bash
echo 'export NXF_SYNTAX_PARSER=v2' >> ~/.bashrc
source ~/.bashrc
```

!!! note "Temporary requirement"

    The `NXF_SYNTAX_PARSER=v2` environment variable is a temporary requirement.
    From Nextflow 26.04 onward, the v2 parser will become the default and this setting will no longer be needed.

### Disabling the v2 syntax parser

To **disable** the v2 syntax parser for your current session, run the following command in your terminal:

```bash
export NXF_SYNTAX_PARSER=v1
```

### Migrating your code

TODO: [refer to the documentation for full details]

TODO: [Provide some summary guidance on using the linter with the `-format` option.]

---

<div markdown class="homepage_logos">

![Seqera](assets/img/seqera_logo.png#only-light)

![Seqera](assets/img/seqera_logo_dark.png#only-dark)

</div>
