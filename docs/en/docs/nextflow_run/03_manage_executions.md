# Part 3: Manage workflow executions

As you run and re-run pipelines, you accumulate execution history and old `work/` directories.
In [Part 1](./01_run_nextflow.md#23-use-resume-to-skip-completed-work) you already used `-resume` to skip work that was already done.
Here you'll learn how to inspect the history of past runs with [`nextflow log`](https://nextflow.io/docs/latest/reference/cli.html#log), and how to delete old work directories you no longer need with [`nextflow clean`](https://nextflow.io/docs/latest/reference/cli.html#clean).

---

## 1. Inspect the log of past executions

Whether you're developing a pipeline or running it in production, at some point you'll need to look up information about past runs.

### 1.1. The history file

Every time you launch a Nextflow workflow, a line gets written to a log file called `history`, under a hidden directory called `.nextflow` in the current working directory.

??? abstract "File contents"

    ```txt title=".nextflow/history" linenums="1"
    2026-09-12 20:22:37	3s	scruffy_wescoff	OK	c3c85dec78d428e60f0168172f4a7b10	f20fdf93-d052-4559-8e55-779b2d1d26ad	nextflow run main.nf --input data/greetings.csv --character turkey
    2026-09-12 20:22:42	3.2s	sleepy_einstein	OK	c3c85dec78d428e60f0168172f4a7b10	3ccd353e-d7a6-4363-9fd3-3145ca537bda	nextflow run main.nf --input data/greetings.csv --character tux
    2026-09-12 20:22:48	1.9s	golden_golick	OK	c3c85dec78d428e60f0168172f4a7b10	073dea98-65b6-4ece-8bf9-f74c1f5813e9	nextflow run main.nf --input data/greetings.csv --character stegosaurus
    2026-09-12 20:22:55	1.1s	distraught_heyrovsky	OK	c3c85dec78d428e60f0168172f4a7b10	073dea98-65b6-4ece-8bf9-f74c1f5813e9	nextflow run main.nf --input data/greetings.csv --character stegosaurus -resume
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
    TIMESTAMP          	DURATION	RUN NAME            	STATUS	REVISION ID	SESSION ID                          	COMMAND
    2026-09-12 20:22:37	3s      	scruffy_wescoff     	OK    	c3c85dec78 	f20fdf93-d052-4559-8e55-779b2d1d26ad	nextflow run main.nf --input data/greetings.csv --character turkey
    2026-09-12 20:22:42	3.2s    	sleepy_einstein     	OK    	c3c85dec78 	3ccd353e-d7a6-4363-9fd3-3145ca537bda	nextflow run main.nf --input data/greetings.csv --character tux
    2026-09-12 20:22:48	1.9s    	golden_golick       	OK    	c3c85dec78 	073dea98-65b6-4ece-8bf9-f74c1f5813e9	nextflow run main.nf --input data/greetings.csv --character stegosaurus
    2026-09-12 20:22:55	1.1s    	distraught_heyrovsky	OK    	c3c85dec78 	073dea98-65b6-4ece-8bf9-f74c1f5813e9	nextflow run main.nf --input data/greetings.csv --character stegosaurus -resume
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

Look up the most recent run you want to keep using `nextflow log`; in the example above that's `golden_golick`, the last plain run before the `-resume` one.
The run name is the machine-generated two-part string shown in the `Launching (...)` console line, or in the `RUN NAME` column of `nextflow log`.

### 2.2. Do a dry run

Add `-n` first to check what a given command would delete without actually deleting anything:

```bash
nextflow clean -before golden_golick -n
```

??? success "Command output"

    ```console
    Would remove /workspaces/training/nextflow-run/work/95/303e4371ccfd87fcc86eb0fd5bae83
    Would remove /workspaces/training/nextflow-run/work/54/65b396bba6cf13e920e4a2d7c92e06
    Would remove /workspaces/training/nextflow-run/work/1d/663cba3fe358fdef2bce4371541660
    Would remove /workspaces/training/nextflow-run/work/7b/b8eb1ddc338e8f8362231f6474cf77
    Would remove /workspaces/training/nextflow-run/work/cd/0e51e58df1f07a15fd66184fea8eac
    Would remove /workspaces/training/nextflow-run/work/26/e360181f20ea2eb1eb59fdb211cb63
    Would remove /workspaces/training/nextflow-run/work/e3/6606b42049a5c795eab3377248aef3
    Would remove /workspaces/training/nextflow-run/work/7e/2b4474a7c7a1fd715dac963d18f5d1
    Would remove /workspaces/training/nextflow-run/work/de/20a5c34be23b0a96f990a1452352f5
    Would remove /workspaces/training/nextflow-run/work/66/ff4cccb589594c18316670eaca9f94
    Would remove /workspaces/training/nextflow-run/work/39/b13463fee8af8f9f4215feea42d22e
    Would remove /workspaces/training/nextflow-run/work/0e/1ec827e23bf79fa99a18f3dc928ec3
    Would remove /workspaces/training/nextflow-run/work/b8/3058ea8a0a2e9104abafc5716ed0fe
    Would remove /workspaces/training/nextflow-run/work/d0/7a9c634b5516e56f9523c27d0a6ec5
    Would remove /workspaces/training/nextflow-run/work/ba/7ccecf594af6b10ccf50e953bbe777
    Would remove /workspaces/training/nextflow-run/work/c2/bcd477ecb624ed38bd26e32d57c8af
    ```

That's 16 task directories: the 8 tasks from the `turkey` run plus the 8 from the `tux` run, exactly as many as you'd expect for two full runs of this four-process pipeline.
The `golden_golick` run itself, and the cached tasks the `-resume` run reused from it, are left alone.

Your output will list different directory names, and how many lines you get depends on how many runs you've done; if you don't see any lines, either the run name doesn't match one in your log, or there's nothing to delete before it.

### 2.3. Proceed with deletion

Once the dry run looks right, re-run the same command with `-f` instead of `-n`:

```bash
nextflow clean -before golden_golick -f
```

??? success "Command output"

    ```console
    Removed /workspaces/training/nextflow-run/work/95/303e4371ccfd87fcc86eb0fd5bae83
    Removed /workspaces/training/nextflow-run/work/54/65b396bba6cf13e920e4a2d7c92e06
    Removed /workspaces/training/nextflow-run/work/1d/663cba3fe358fdef2bce4371541660
    Removed /workspaces/training/nextflow-run/work/7b/b8eb1ddc338e8f8362231f6474cf77
    Removed /workspaces/training/nextflow-run/work/cd/0e51e58df1f07a15fd66184fea8eac
    Removed /workspaces/training/nextflow-run/work/26/e360181f20ea2eb1eb59fdb211cb63
    Removed /workspaces/training/nextflow-run/work/e3/6606b42049a5c795eab3377248aef3
    Removed /workspaces/training/nextflow-run/work/7e/2b4474a7c7a1fd715dac963d18f5d1
    Removed /workspaces/training/nextflow-run/work/de/20a5c34be23b0a96f990a1452352f5
    Removed /workspaces/training/nextflow-run/work/66/ff4cccb589594c18316670eaca9f94
    Removed /workspaces/training/nextflow-run/work/39/b13463fee8af8f9f4215feea42d22e
    Removed /workspaces/training/nextflow-run/work/0e/1ec827e23bf79fa99a18f3dc928ec3
    Removed /workspaces/training/nextflow-run/work/b8/3058ea8a0a2e9104abafc5716ed0fe
    Removed /workspaces/training/nextflow-run/work/d0/7a9c634b5516e56f9523c27d0a6ec5
    Removed /workspaces/training/nextflow-run/work/ba/7ccecf594af6b10ccf50e953bbe777
    Removed /workspaces/training/nextflow-run/work/c2/bcd477ecb624ed38bd26e32d57c8af
    ```

`nextflow clean` empties the task directories but leaves the two-character parent directories (like `95/`) in place.

!!! warning

    Deleting work directories from past runs removes them from Nextflow's cache and deletes any outputs stored only there.
    That breaks Nextflow's ability to resume execution without re-running the corresponding processes, so only clean up runs you're confident you won't need to resume from.
    This is also why it's worth publishing anything you care about to `results/` with `mode 'copy'` rather than relying on the `work/` directory or a `symlink` publish mode.

### Takeaway

You know how to remove old work directories with `nextflow clean`, and why doing so trades away the ability to resume from those runs.

### What's next?

Learn how to run pipelines directly from remote repositories such as GitHub in [Part 4](./04_remote_repositories.md).

---

## Summary

In this part you learned to:

- Inspect the history of past runs with `nextflow log`
- Remove old work directories with `nextflow clean`, and understand the resume trade-off that comes with it
