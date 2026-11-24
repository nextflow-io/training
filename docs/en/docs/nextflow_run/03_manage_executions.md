# Part 3: Manage workflow executions

As you run and re-run pipelines, you accumulate execution history and old `work/` directories.
In [Part 1](./01_run_nextflow.md#21-run-the-workflow-with-resume) you already used `-resume` to skip work that was already done.
Here you'll learn how to inspect the history of past runs with [`nextflow log`](https://nextflow.io/docs/latest/reference/cli.html#log) and delete old work directories you no longer need with [`nextflow clean`](https://nextflow.io/docs/latest/reference/cli.html#clean).

This part is optional.
If time is short, skip ahead to [Part 4](./04_remote_repositories.md), which does not depend on it.

---

## 1. Inspect the log of past executions

Whether you're developing a pipeline or running it in production, at some point you'll need to look up information about past runs.

### 1.1. The history file

Every time you launch a Nextflow workflow, a line gets written to a log file called `history`, under a hidden directory called `.nextflow` in the current working directory.

??? abstract "File contents"

    ```txt title=".nextflow/history" linenums="1"
    2026-11-23 23:37:27	2.5s	soggy_heyrovsky	OK	f5f2664c7a42831ad27cde8269282db8	43d52f9c-245a-4b8c-a22a-e972440278fe	nextflow run main.nf --input data/greetings.csv --character turkey
    2026-11-23 23:37:31	2.3s	cheesy_ptolemy	OK	f5f2664c7a42831ad27cde8269282db8	b19e5c7e-693b-411e-924b-827e0cea9423	nextflow run main.nf --input data/greetings.csv --character tux
    2026-11-23 23:37:35	2.6s	silly_mercator	OK	f5f2664c7a42831ad27cde8269282db8	aeb89dce-8133-4a2f-b622-ef67fb9309e1	nextflow run main.nf --input data/greetings.csv --character stegosaurus
    2026-11-23 23:37:39	1.9s	silly_galileo	OK	f5f2664c7a42831ad27cde8269282db8	aeb89dce-8133-4a2f-b622-ef67fb9309e1	nextflow run main.nf --input data/greetings.csv --character stegosaurus -resume
    ```

Each line gives you the timestamp, duration, run name, status, revision ID, session ID, and full command line for a run launched from this directory.

Look at the last two lines: they're two separate invocations (one plain, one with `-resume`) of the exact same command, and they share the same session ID.
The session ID only changes when you launch a genuinely new run; using `-resume` keeps it, which is how Nextflow knows which cache to reuse.

### 1.2. Use `nextflow log` for a friendlier view

Reading the raw history file works, but `nextflow log` formats the same information with a header:

```bash
nextflow log
```

??? success "Command output"

    ```console linenums="1"
    TIMESTAMP          	DURATION	RUN NAME       	STATUS	REVISION ID	SESSION ID                          	COMMAND
    2026-11-23 23:37:27	2.5s    	soggy_heyrovsky	OK    	f5f2664c7a 	43d52f9c-245a-4b8c-a22a-e972440278fe	nextflow run main.nf --input data/greetings.csv --character turkey
    2026-11-23 23:37:31	2.3s    	cheesy_ptolemy 	OK    	f5f2664c7a 	b19e5c7e-693b-411e-924b-827e0cea9423	nextflow run main.nf --input data/greetings.csv --character tux
    2026-11-23 23:37:35	2.6s    	silly_mercator 	OK    	f5f2664c7a 	aeb89dce-8133-4a2f-b622-ef67fb9309e1	nextflow run main.nf --input data/greetings.csv --character stegosaurus
    2026-11-23 23:37:39	1.9s    	silly_galileo  	OK    	f5f2664c7a 	aeb89dce-8133-4a2f-b622-ef67fb9309e1	nextflow run main.nf --input data/greetings.csv --character stegosaurus -resume
    ```

Nextflow groups the caching information it uses for `-resume` under `.nextflow/cache`, keyed by session ID.
That's why looking up the right run name or session ID here is the first step whenever you need to investigate or clean up a past execution.

### Takeaway

You know where Nextflow records the history of past runs, and how to inspect it with `nextflow log`.

### What's next?

Learn how to remove old work directories you no longer need.

---

## 2. Delete older work directories

Every run leaves its task directories behind under `work/`, even after you've copied the outputs you care about to `results/`.
Run enough pipelines during development and those subdirectories add up, so Nextflow provides `nextflow clean` to remove the ones you no longer need.

### 2.1. Determine deletion criteria

`nextflow clean` supports several ways to select what to remove; see the [reference documentation](https://www.nextflow.io/docs/latest/reference/cli.html#clean) for the full list.
Here you'll delete everything from runs before a given run, using its run name.

Look up the most recent run you want to keep using `nextflow log`; in the [example from 1.2](#12-use-nextflow-log-for-a-friendlier-view) that's `silly_mercator`, the last plain run before the `-resume` one.
The run name is the machine-generated two-part string shown in the `Launching (...)` console line, or in the `RUN NAME` column of `nextflow log`.

### 2.2. Do a dry run

Add `-n` first to check what a given command would delete without actually deleting anything:

```bash
nextflow clean -before silly_mercator -n
```

??? success "Command output"

    ```console
    Would remove /workspaces/training/nextflow-run/work/7c/5801ce6e01159a03b56278824151e5
    Would remove /workspaces/training/nextflow-run/work/82/5a06fd4cb701738f88022546aaedea
    Would remove /workspaces/training/nextflow-run/work/2f/130379895837b7a95fe652adc335a3
    Would remove /workspaces/training/nextflow-run/work/3a/b9f6c79dd7bd225e63a2024f10ea63
    Would remove /workspaces/training/nextflow-run/work/7a/b2502a27862d8b859b46f86fe110e1
    Would remove /workspaces/training/nextflow-run/work/99/b0f0a37ab91cd109e4f40f7eeccc88
    Would remove /workspaces/training/nextflow-run/work/23/2ba361a071b20dba45fd6b271247f1
    Would remove /workspaces/training/nextflow-run/work/d2/c79a535f01076389275c2b66574c8f
    Would remove /workspaces/training/nextflow-run/work/bf/da869e3a9d07e755c86287f451bfdb
    Would remove /workspaces/training/nextflow-run/work/be/8f9ba84dab8e42e4b8e6c7e619fe03
    Would remove /workspaces/training/nextflow-run/work/4b/6553fb233c0ae93f0827c0bcdd0e88
    Would remove /workspaces/training/nextflow-run/work/a4/e13b1ad066c5ab77a6731e68b3cdaa
    Would remove /workspaces/training/nextflow-run/work/9f/b630d59ea3eef8c1b712e34fbf2d96
    Would remove /workspaces/training/nextflow-run/work/a0/2a5226725af0509c4b00f800364853
    Would remove /workspaces/training/nextflow-run/work/20/6212cc5a86a58e100935966c9f74eb
    Would remove /workspaces/training/nextflow-run/work/8a/0707e0d57c39f1e0ae62bb2c6cc034
    ```

That's 16 task directories: the 8 tasks from the `turkey` run plus the 8 from the `tux` run, exactly as many as you'd expect for two full runs of this four-process pipeline.
The `silly_mercator` run itself, and the cached tasks the `-resume` run reused from it, are left alone.

Your output will list different directory names, and how many lines you get depends on how many runs you've done; if you don't see any lines, either the run name doesn't match one in your log, or there's nothing to delete before it.

### 2.3. Proceed with deletion

Once the dry run looks right, re-run the same command with `-f` instead of `-n`:

```bash
nextflow clean -before silly_mercator -f
```

??? success "Command output"

    ```console
    Removed /workspaces/training/nextflow-run/work/7c/5801ce6e01159a03b56278824151e5
    Removed /workspaces/training/nextflow-run/work/82/5a06fd4cb701738f88022546aaedea
    Removed /workspaces/training/nextflow-run/work/2f/130379895837b7a95fe652adc335a3
    Removed /workspaces/training/nextflow-run/work/3a/b9f6c79dd7bd225e63a2024f10ea63
    Removed /workspaces/training/nextflow-run/work/7a/b2502a27862d8b859b46f86fe110e1
    Removed /workspaces/training/nextflow-run/work/99/b0f0a37ab91cd109e4f40f7eeccc88
    Removed /workspaces/training/nextflow-run/work/23/2ba361a071b20dba45fd6b271247f1
    Removed /workspaces/training/nextflow-run/work/d2/c79a535f01076389275c2b66574c8f
    Removed /workspaces/training/nextflow-run/work/bf/da869e3a9d07e755c86287f451bfdb
    Removed /workspaces/training/nextflow-run/work/be/8f9ba84dab8e42e4b8e6c7e619fe03
    Removed /workspaces/training/nextflow-run/work/4b/6553fb233c0ae93f0827c0bcdd0e88
    Removed /workspaces/training/nextflow-run/work/a4/e13b1ad066c5ab77a6731e68b3cdaa
    Removed /workspaces/training/nextflow-run/work/9f/b630d59ea3eef8c1b712e34fbf2d96
    Removed /workspaces/training/nextflow-run/work/a0/2a5226725af0509c4b00f800364853
    Removed /workspaces/training/nextflow-run/work/20/6212cc5a86a58e100935966c9f74eb
    Removed /workspaces/training/nextflow-run/work/8a/0707e0d57c39f1e0ae62bb2c6cc034
    ```

`nextflow clean` empties the task directories but leaves the two-character parent directories (like `7c/`) in place.

!!! warning

    Deleting work directories from past runs removes them from Nextflow's cache and deletes any outputs stored only there.
    That breaks Nextflow's ability to resume execution without re-running the corresponding processes, so only clean up runs you're confident you won't need to resume from.
    This is also why it's worth publishing anything you care about to `results/` with `mode 'copy'` rather than relying on the `work/` directory or a `symlink` publish mode.

    `nextflow clean -f` also removes the cleaned runs from the execution history, without telling you: afterwards they no longer appear in `nextflow log`.
    Note down any run names or session IDs you still need before you clean.

### Takeaway

You know how to remove old work directories with `nextflow clean`, and why doing so trades away the ability to resume from those runs.

### What's next?

Learn how to run pipelines directly from remote repositories such as GitHub in [Part 4](./04_remote_repositories.md).

---

## Summary

In this part you learned to:

- Inspect the history of past runs with `nextflow log`
- Remove old work directories with `nextflow clean`, and understand the resume trade-off that comes with it
