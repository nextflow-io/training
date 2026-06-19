# Develop with CoScientist

Adding a step to a pipeline normally means finding a tool, writing a process, wiring its inputs and outputs, and getting the change reviewed.
CoScientist can do that legwork for you, and the developer journey starts here, still in the web chat.
You connect your GitHub account, fork `rnaseq-nf`, and have CoScientist add a new process by reusing an nf-core module and open a pull request, letting it work until the change is done.

---

## 1. Connect your GitHub account.

Give CoScientist access to your GitHub account by configuring a GitHub token, so it can commit changes and open pull requests on your behalf.
After adding the token, start a new session so CoScientist picks it up; an existing session keeps using the old token.

!!! note "Checkpoint"

    CoScientist reports that it can access your GitHub account.

## 2. Fork rnaseq-nf.

Fork the pipeline into your own GitHub account so you have a copy to change.
On GitHub, open [nextflow-io/rnaseq-nf](https://github.com/nextflow-io/rnaseq-nf) and select **Fork**.

!!! note "Checkpoint"

    The repository `your-user/rnaseq-nf` exists on your GitHub account.

## 3. Add a trimming step and open a pull request.

Ask CoScientist to add a `fastp` read-trimming step before quantification, reusing the nf-core fastp module, and to open a pull request with the change:

```text
In my fork, add a fastp read-trimming step before quantification, reusing the nf-core fastp module, wire its trimmed reads into the QUANT process, and open a pull request.
```

Reusing the nf-core module keeps the work small: CoScientist pulls in a tested process and mainly has to get the input and output wiring right, rather than reimplementing fastp.
This exercises its nf-core module discovery.

??? example "What CoScientist typically does"

    It finds the nf-core/fastp module, adds it under `modules/nf-core/fastp/`, wires the trimmed reads into the workflow, commits to a branch on your fork, and opens a pull request titled something like "add fastp trimming step before quantification".
    The exact wording will differ from run to run.

## 4. Let the agent work until it is done.

Adding a process rarely lands in one attempt.
Let CoScientist iterate: when something is wrong it adjusts the wiring and tries again, until the change is complete and the pull request is open.
Keep it on track with the habits from the previous lesson: check what it changed, and redirect it if it goes too far.

<!-- TODO: verify whether the web interface offers a goal or iterate-until-done mode, and how it surfaces -->

!!! note "Checkpoint"

    A pull request with the `fastp` step is open on your fork.

### Takeaway

Working entirely in the web chat, you connected GitHub, forked `rnaseq-nf`, and had CoScientist add a real processing step by reusing an nf-core module and open a pull request.

### What's next?

In the next lesson, [move to the CLI](04_move_to_the_cli.md) to test the change and see why the CLI suits development work.
