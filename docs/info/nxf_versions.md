---
title: Nextflow versions
description:
  Understanding and managing the evolution of Nextflow's syntax versions
  - toc
  - footer
---

## Current supported Nextflow syntax version & requirements

As of version 3.0 of the training portal, all our training courses are based on the 25.10.2 version release of Nextflow unless otherwise specified on the course index page (except deprecated or otherwise archived materials which may not include a version notice).

Because the courses now use typed inputs at the workflow level as well as workflow-level output directives, they require the use of the V2 syntax parser.
If you plan to use the environment we provide through [Github Codespaces](envsetup/01_setup.md) or [local devcontainers](envsetup/03_devcontainer.md), you don't need to do anything unless specifically noted in the course instructions.
However, if you are planning to work through the trainings in your own environment ([Manual install](envsetup/02_local.md)), you will need to make sure to use Nextflow version 25.10.2 or later with the v2 syntax parser enabled.

## Older versions of the training materials

Our training materials have been versioned since February 2025.

You can access older version of the training materials that work with versions of Nextflow **before 25.10.2** via the dropdown menu item at the top of every page that shows the numbered version of the training materials.
When you select an older version of the training materials, links to the training environment will automatically specify the corresponding version of the environment.

## Other information about Nextflow syntax versions

<!-- TODO: Blurb about Nextflow versions disambiguating DSL1 / DSL2 and strict syntax v2 parser stuff. -->

### Checking and setting the Nextflow version

You can check what version of Nextflow is installed on your system using the command `nextflow --version`.

For more information about how to update your version of Nextflow, please see the reference documentation on [Updating Nextflow](https://www.nextflow.io/docs/latest/updating-nextflow.html).

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

Note that the `NXF_SYNTAX_PARSER=v2` environment variable is a temporary requirement.
From Nextflow 26.04 onward, the v2 parser will become the default and this setting will no longer be needed.

### Disabling the v2 syntax parser

To **disable** the v2 syntax parser for your current session, run the following command in your terminal:

```bash
export NXF_SYNTAX_PARSER=v1
```

<!-- Will it be possible to disable it in versions after 26.04? -->

### Migrating existing code

For guidance regarding the migration of existing code to comply with more recent versions of Nextflow, please see the [Migration Notes](https://www.nextflow.io/docs/latest/migrations/index.html) in the reference documentation.

These two articles are particularly helpful for migrating to the most recent release:

- [Migrating to workflow outputs](https://www.nextflow.io/docs/latest/tutorials/workflow-outputs.html)
- [Migrating to static types](https://www.nextflow.io/docs/latest/tutorials/static-types.html)

Both of these features are covered as part of the beginner training starting in version 3.0 of the training materials.

Depending on the generation of Nextflow code you intend to migrate, you may be able to get most of it done by the Nextflow linter using the `nextflow lint -format` command.
See the CLI reference for [`lint`](https://www.nextflow.io/docs/latest/reference/cli.html#lint) for more details.

We hope this will be helpful.
Don't hesitate to reach out on Slack or on the forum if you need help.

---

<div markdown class="homepage_logos">

![Seqera](assets/img/seqera_logo.png#only-light)

![Seqera](assets/img/seqera_logo_dark.png#only-dark)

</div>
