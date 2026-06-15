# Next Steps

Congratulations on completing the **Building with Seqera AI** training course.

As you put CoScientist to work on your own pipelines, keep doing what you practised here: let the agent move fast, and stay the reviewer who checks its output before relying on it.

---

## 1. What else CoScientist can do.

This side quest followed one end-to-end flow, but CoScientist does considerably more.

- **Build, plan, and goal modes**: make changes directly, plan read-only before acting, or let goal mode keep working until an objective is met.
- **Built-in slash commands**: ready-made workflows in the CLI such as `/debug-last-run-on-seqera`, `/nextflow-config`, and `/nextflow-schema`.
- **Session resumption**: continue your most recent session from the CLI with `seqera ai -c`.
- **Command approval**: choose how much the agent runs automatically (`basic`, `default`, or `full`) so you stay in control of what executes.
- **Code intelligence**: language-server support for Nextflow, Python, and R that flags errors as code is written.
- **Script conversion**: convert R and Python scripts and Jupyter notebooks toward Nextflow.
- **Wave containers**: build containers from `conda` or `pip` packages without writing a Dockerfile.
- **nf-core modules**: discover and reuse from over 1,000 standardized nf-core modules.

See the [CoScientist documentation](https://docs.seqera.io/platform-cloud/co-scientist/) for the full feature set.

!!! info "Cloud and Enterprise"

    This training runs on Seqera Cloud, which uses Seqera's hosted model and meters usage through credits.
    On Seqera Enterprise (currently AWS only), CoScientist can run Claude inference inside your own AWS account through Amazon Bedrock, so prompts and data stay within that account.

---

## 2. Keep building with CoScientist.

Here are three ways to continue developing with CoScientist.

### 2.1. Apply CoScientist to your own pipelines.

Try CoScientist against one of your own Nextflow pipelines or a real nf-core pipeline.
Working on a codebase you know well is the fastest way to see where CoScientist adds value and where you need to steer it.

### 2.2. Use the CLI in your development environment.

The `seqera ai` CLI brings the same agent into your local terminal.
Install it and run it alongside your usual development tools so CoScientist can read and edit files directly in your working directory.

### 2.3. Read the documentation.

The [CoScientist documentation](https://docs.seqera.io/platform-cloud/co-scientist/) covers capabilities, configuration, and usage patterns in detail.
Refer to it when you want to understand what the agent can and cannot do, or when you need to configure it for a specific environment.

---

## 3. Get help from the community.

- [Nextflow Slack](https://www.nextflow.io/slack-invite.html): Ask questions and share your work
- [Community forum](https://community.seqera.io/): Discuss ideas and get advice

---

## 4. Continue your Nextflow training.

If you haven't already, check out our other training courses:

- **[Hello Nextflow](../../hello_nextflow/index.md)**: Foundational Nextflow concepts
- **[Testing with nf-test](../nf_test/index.md)**: Writing and running tests for workflows
- **[Troubleshooting Workflows](../debugging/index.md)**: Identifying and fixing common workflow errors
- **[Side Quests](../index.md)**: Deep dives into specific topics
