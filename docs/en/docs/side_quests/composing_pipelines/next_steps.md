# Next Steps

Congratulations on completing the **Composing Pipelines** training course.

---

## 1. Top 3 ways to continue your composition journey

Here are three recommendations for what to do next.

### 1.1. Study a real-world composed pipeline

The greetings and report pipelines in this course are deliberately small.
To see the same patterns at full scale, look at [rnaseq-diffabundance-meta](https://github.com/pinin4fjords/rnaseq-diffabundance-meta), which composes nf-core/rnaseq and nf-core/differentialabundance into a single run.
The quantification outputs of one pipeline feed straight into the other, in one DAG, without an intermediate samplesheet or a second launch.

It also answers practical questions that this course only touches on:

- **Bringing pipelines into the project**: the component pipelines are vendored, copied into `pipelines/nf-core/`, with a `pipelines.json` file recording the repository, branch and commit of each, in the same shape as an nf-core `modules.json`.
- **Configuration**: since an included pipeline's `nextflow.config` isn't loaded, each component pipeline keeps its defaults and process configuration in separate files that both its own `nextflow.config` and the composed pipeline include.
- **Tool arguments**: arguments are passed in as process inputs rather than read from `params` in config closures, so options like `--rnaseq.extra_star_align_args` reach the right tool when the pipeline is included.

Treat it as a worked example of the patterns rather than a template to copy as-is.

### 1.2. Read the documentation on pipeline composition

The Nextflow documentation describes the full syntax and its limits.
Start with the [pipeline composition](https://www.nextflow.io/docs/latest/workflow-typed.html#pipeline-composition) section of the typed workflows page, then read [static typing](https://www.nextflow.io/docs/latest/static-typing.html) for the type system it builds on.

### 1.3. Make one of your own pipelines composable

Pick a pipeline you maintain and work through the checklist from Part 6.
Even if nobody includes it yet, typed interfaces, declared inputs and a clean `params {}` and `output {}` contract make it easier to test, to document and to reuse later.
If your pipeline isn't typed yet, the [migration guide](https://www.nextflow.io/docs/latest/tutorials/static-types.html) shows how to convert it step by step.

---

## 2. Get help from the community

- [Nextflow Slack](https://www.nextflow.io/slack-invite.html): Ask questions and share your work
- [Community forum](https://community.seqera.io/): Discuss ideas and get advice
- [GitHub discussions](https://github.com/nextflow-io/nextflow/discussions): Technical questions and feature requests

Feedback from people trying pipeline composition on real pipelines is especially valuable.

---

## 3. Continue your Nextflow training

If you haven't already, check out our other training courses:

- **[Types and Records](../types_and_records/index.md)**: Static types and records in depth
- **[Build with nf-core](../../nfcore_build/index.md)**: nf-core pipelines and best practices
- **[Side Quests](../index.md)**: Deep dives into specific topics
