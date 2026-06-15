# Develop on the CLI

The user journey in the previous lessons all happened in the browser chat.
Now you become a developer and move to the `seqera ai` command-line agent to change the pipeline itself.
You fork `rnaseq-nf`, add a real processing step, and run the result on your own machine.

---

## 1. Why move to the CLI.

The web chat is built for the user journey; developing a pipeline is better suited to your own terminal:

- Compute: the CLI uses your machine, which has more resources than the in-browser environment.
- Docker: it is available locally, so containerized processes and pipeline runs work.
- Local access: the agent works against your own files, repos, editor, and personal credentials and config.
- Iteration: you edit, run, and see results in one terminal, without round-tripping through the web chat.

## 2. Install and start the CLI.

Install the CLI (it needs Node.js 18 or later) and sign in:

```bash
npm install -g seqera
seqera login
```

`seqera login` opens a browser to authenticate against your Seqera account.
Then start an interactive session:

```bash
seqera ai
```

## 3. Fork rnaseq-nf.

Ask the agent to fork the pipeline into your own GitHub account and clone it locally so you can work on it:

```text
Fork nextflow-io/rnaseq-nf into my GitHub account and clone it locally so I can work on it.
```

The agent uses git and the GitHub CLI in your terminal to create the fork and bring it down.

!!! note "Checkpoint"

    A local clone of your fork of `rnaseq-nf` exists on your machine.

<!-- TODO: verify the agent's GitHub fork+clone flow on the CLI (gh auth) against the product -->

## 4. Add a trimming step.

Ask CoScientist to add a `fastp` read-trimming process before quantification, reusing an existing nf-core module rather than writing one from scratch:

```text
Add a fastp read-trimming step to the pipeline before quantification, reusing the nf-core fastp module, and wire its trimmed reads into the existing QUANT process.
```

This exercises CoScientist's nf-core module discovery: it finds a tested module and integrates it instead of generating a new process by hand.

??? example "What CoScientist typically does"

    It finds the nf-core/fastp module, adds the process to the pipeline, and updates the workflow to pass the trimmed reads into `QUANT` downstream.
    The exact wording will differ from run to run.

<!-- TODO: verify nf-core fastp module integration into rnaseq-nf works as described -->

## 5. Run it.

Run the modified pipeline locally to confirm the new step works end to end:

```text
Run the pipeline locally with the docker profile and tell me if it completes.
```

`rnaseq-nf` has no `test` profile; its default parameters point at the bundled `ggal_gut` test data, so the command is `nextflow run . -profile docker`.
If the first run fails, which a newly added process often does, ask the agent to read the error and fix the wiring.

!!! note "Checkpoint"

    The pipeline runs to completion with the `fastp` step included.

### Takeaway

You moved from the web chat to the `seqera ai` CLI, forked `rnaseq-nf` into your own account, and used CoScientist to add a real processing step.
The step reuses a tested nf-core module rather than new code you have to maintain.

### What's next?

In the next lesson, [test the change, add CI, and open a PR](04_test_ci_and_contribute.md).
