# Part 2: The Nextflow debugging toolkit

When an error message doesn't immediately tell you what's wrong, you reach for the toolkit.

This lesson introduces six techniques, applied in order to the same small pipeline:

1. **Work-directory forensics** when a process fails.
2. **`-preview`** to validate workflow logic before running anything.
3. **`debug true`** to stream process output as it runs.
4. **`-stub-run`** to iterate on workflow logic without running the real commands.
5. **`-dump-hashes`** to find out why `-resume` re-ran something you thought was cached.
6. **A systematic methodology** that ties the tools together.

Each section uses the same pipeline so you can see how the tools fit together in a real development workflow.

---

## 0. Get started

Move into the working directory if you aren't already there:

```bash
cd side-quests/debugging
```

The pipeline we'll use throughout this lesson is `sample_processing.nf`.
It reads a sample manifest, counts lines in each sample's gzipped FASTQ, and writes a small report:

```groovy title="sample_processing.nf" linenums="1"
#!/usr/bin/env nextflow

params.input = 'data/sample_data.csv'
params.outdir = 'results'

process COUNT_LINES {

    container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
    publishDir params.outdir, mode: 'copy'

    input:
    tuple val(sample_id), path(fastq)

    output:
    path "${sample_id}.count.txt"

    script:
    """
    lines=\$(zcat ${fastq} | wc -l)
    cowpy "${sample_id} has \${lines} lines" > ${sample_id}_count.txt
    """

    stub:
    """
    echo "${sample_id} has 0 lines" > ${sample_id}_count.txt
    """
}

process REPORT {

    publishDir params.outdir, mode: 'copy'

    input:
    path count_files

    output:
    path 'report.txt'

    script:
    """
    cat ${count_files} > report.txt
    """
}

workflow {

    samples_ch = channel
        .fromPath(params.input)
        .splitCsv(header: true)
        .map { row -> [row.sample_id, file(row.fastq_path)] }

    counts_ch = COUNT_LINES(samples_ch)

    REPORT(counts_ch.collect())
}
```

The pipeline ships with one intentional bug, which we'll fix in section 1.
The rest of the sections build on the fixed pipeline.

If at any point you want to start over, run:

```bash
git checkout sample_processing.nf
```

to restore the original file.

---

## 1. Work-directory forensics

When a process fails, Nextflow creates a work directory containing everything that ran: the command, its output, its error stream, and its exit code.
This is the first place to look when an error message doesn't tell you enough on its own.

### 1.1. Run the pipeline

```bash
nextflow run sample_processing.nf -profile docker
```

??? failure "Command output"

    ```console
    N E X T F L O W   ~  version 25.10.4

    Launching `sample_processing.nf` [chaotic_meninsky] DSL2 - revision: 4a8a16d3a2

    executor >  local (5)
    [34/6278a1] COUNT_LINES (3) | 0 of 5 ✘
    [-        ] REPORT          -
    ERROR ~ Error executing process > 'COUNT_LINES (5)'

    Caused by:
      Missing output file(s) `sample_005.count.txt` expected by process `COUNT_LINES (5)`


    Command executed:

      lines=$(zcat sample_005.fastq.gz | wc -l)
      cowpy "sample_005 has ${lines} lines" > sample_005_count.txt

    Command exit status:
      0

    Command output:
      (empty)

    Work dir:
      /workspaces/training/side-quests/debugging/work/18/8e48f60e7b68419a62909ca8c5bd3d

    Tip: you can try to figure out what's wrong by changing to the process work dir and showing the script file named `.command.sh`

     -- Check '.nextflow.log' file for details
    ```

The process exited cleanly (`Command exit status: 0`) but Nextflow couldn't find the output file it was promised.
The error already hints at the cause: the script writes `sample_005_count.txt`, but the process declares `sample_005.count.txt`.
Let's confirm that by looking inside the work directory.

### 1.2. Walk the work directory

The error message gives you the path. Copy it from your own terminal output (yours will have a different hash from ours).

Every failed task has the same set of hidden files. Let's go through them in turn.

#### 1.2.1. `.command.sh` — the executed command

This is the exact script Nextflow ran, after variable substitution:

```bash
cat work/18/8e48f60e7b68419a62909ca8c5bd3d/.command.sh
```

```console title="Output"
#!/bin/bash -ue
lines=$(zcat sample_005.fastq.gz | wc -l)
cowpy "sample_005 has ${lines} lines" > sample_005_count.txt
```

Notice that Nextflow has already substituted `${sample_id}` for `sample_005` and the redirect target is `sample_005_count.txt`.
The script ran exactly as written.

#### 1.2.2. `.command.err` and `.command.out` — what the command printed

```bash
cat work/18/8e48f60e7b68419a62909ca8c5bd3d/.command.err
```

```console title="Output"
(empty)
```

```bash
cat work/18/8e48f60e7b68419a62909ca8c5bd3d/.command.out
```

```console title="Output"
(empty)
```

Nothing on either stream, because `cowpy` wrote its output to the redirected file.

#### 1.2.3. `.exitcode` — how the command exited

```bash
cat work/18/8e48f60e7b68419a62909ca8c5bd3d/.exitcode
```

```console title="Output"
0
```

Zero, as we already knew from the error.
This confirms the failure is on Nextflow's side, not the command's: the command succeeded, but it produced a file with a name Nextflow wasn't expecting.

Common exit codes to recognise:

- **0**: success (and yet here we are — usually means an output mismatch)
- **127**: command not found (software not installed in the container)
- **137**: killed by the scheduler (memory or time limit exceeded)

#### 1.2.4. `ls` the directory — what files actually exist

```bash
ls work/18/8e48f60e7b68419a62909ca8c5bd3d/
```

```console title="Output"
sample_005.fastq.gz  sample_005_count.txt
```

There it is: the file `sample_005_count.txt` exists, but the process expects `sample_005.count.txt`.

### 1.3. Fix the bug

The output declaration uses a `.` (dot), the script uses an `_` (underscore). Make them match.
Fix the script and the stub block:

=== "After"

    ```groovy title="sample_processing.nf" hl_lines="11 16 21" linenums="6"
    process COUNT_LINES {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
        publishDir params.outdir, mode: 'copy'

        input:
        tuple val(sample_id), path(fastq)

        output:
        path "${sample_id}.count.txt"

        script:
        """
        lines=\$(zcat ${fastq} | wc -l)
        cowpy "${sample_id} has \${lines} lines" > ${sample_id}.count.txt
        """

        stub:
        """
        echo "${sample_id} has 0 lines" > ${sample_id}.count.txt
        """
    }
    ```

=== "Before"

    ```groovy title="sample_processing.nf" hl_lines="11 16 21" linenums="6"
    process COUNT_LINES {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
        publishDir params.outdir, mode: 'copy'

        input:
        tuple val(sample_id), path(fastq)

        output:
        path "${sample_id}.count.txt"

        script:
        """
        lines=\$(zcat ${fastq} | wc -l)
        cowpy "${sample_id} has \${lines} lines" > ${sample_id}_count.txt
        """

        stub:
        """
        echo "${sample_id} has 0 lines" > ${sample_id}_count.txt
        """
    }
    ```

### 1.4. Re-run and confirm

```bash
nextflow run sample_processing.nf -profile docker -resume
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 25.10.4

    Launching `sample_processing.nf` [friendly_tuckerman] DSL2 - revision: 62a54f6852

    executor >  local (6)
    [65/82c6f9] COUNT_LINES (2) | 5 of 5 ✔
    [2d/d0078c] REPORT          | 1 of 1 ✔
    ```

Note we used `-resume`. With this single-fault example resume doesn't save you anything, but in real pipelines it lets you retry just the failed steps after a fix.

### Takeaway

The work directory contains the complete record of what the process did.
`.command.sh` shows the executed command, `.command.err` and `.command.out` show its streams, `.exitcode` reveals how it finished, and a plain `ls` shows what files exist.
For any process failure: read the error, find the work directory, walk through these files until the cause is obvious.

### What's next?

Learn how to validate workflow logic before running anything, with `-preview`.

---

## 2. Validate workflow logic with `-preview`

Before running an expensive workflow, you can ask Nextflow to _compile_ it and report what _would_ execute, without running any tasks.
This catches syntax errors, missing files, and bad workflow structure in seconds rather than minutes.

### 2.1. Preview the working pipeline

```bash
nextflow run sample_processing.nf -preview
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 25.10.4

    Launching `sample_processing.nf` [ecstatic_curran] DSL2 - revision: 62a54f6852

    [-        ] COUNT_LINES -
    [-        ] REPORT      -
    ```

The dashes mean "this process would run, but nothing was executed".
The script parsed cleanly, the workflow DAG is well-formed, and the two processes show up as expected.

### 2.2. Catch a syntax error before running

Now intentionally break the file.
Open `sample_processing.nf` and delete the closing brace at the end of the workflow block (the `}` on line 55):

=== "Broken"

    ```groovy title="sample_processing.nf" hl_lines="10" linenums="45"
    workflow {

        samples_ch = channel
            .fromPath(params.input)
            .splitCsv(header: true)
            .map { row -> [row.sample_id, file(row.fastq_path)] }

        counts_ch = COUNT_LINES(samples_ch)

        REPORT(counts_ch.collect())
    ```

=== "Working"

    ```groovy title="sample_processing.nf" hl_lines="10" linenums="45"
    workflow {

        samples_ch = channel
            .fromPath(params.input)
            .splitCsv(header: true)
            .map { row -> [row.sample_id, file(row.fastq_path)] }

        counts_ch = COUNT_LINES(samples_ch)

        REPORT(counts_ch.collect())
    }
    ```

Try `-preview` again:

```bash
nextflow run sample_processing.nf -preview
```

??? failure "Command output"

    ```console
    N E X T F L O W   ~  version 25.10.4

    Launching `sample_processing.nf` [...] DSL2 - revision: ...

    ERROR ~ Script compilation error
    - file : /workspaces/training/side-quests/debugging/sample_processing.nf
    - cause: Unexpected input: '{' @ line 45, column 10.
       workflow {
                ^

    1 error
    ```

The error is identical to what you'd see if you ran the workflow normally, but it took milliseconds and didn't pull a container or run any tasks.

Restore the brace before moving on:

```bash
git checkout sample_processing.nf
```

This also re-introduces the output-mismatch bug from section 1, so re-apply your earlier fix.

### Takeaway

`-preview` parses your workflow and shows the planned execution without running anything.
Use it as a fast sanity check before launching expensive pipelines or after refactoring workflow logic.

### What's next?

Learn how to see what a process is actually receiving as input at runtime, with `debug true`.

---

## 3. Stream process output with `debug true`

The work directory is a great forensic tool _after_ a failure, but sometimes you want to see what a process is doing _while_ it's running.
The `debug true` directive streams `stdout` and `stderr` to your terminal as the process runs.
Combined with a strategic `echo`, this gives you the same kind of visibility you'd get from `print` statements during local development.

### 3.1. Add `debug true` and an `echo`

Edit `sample_processing.nf` and add two lines to the COUNT_LINES process:

=== "After"

    ```groovy title="sample_processing.nf" hl_lines="5 18" linenums="6"
    process COUNT_LINES {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
        publishDir params.outdir, mode: 'copy'
        debug true

        input:
        tuple val(sample_id), path(fastq)

        output:
        path "${sample_id}.count.txt"

        script:
        """
        lines=\$(zcat ${fastq} | wc -l)
        cowpy "${sample_id} has \${lines} lines" > ${sample_id}.count.txt
        echo "DEBUG: processed ${sample_id} (\${lines} lines)"
        """
    ```

=== "Before"

    ```groovy title="sample_processing.nf" linenums="6"
    process COUNT_LINES {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
        publishDir params.outdir, mode: 'copy'

        input:
        tuple val(sample_id), path(fastq)

        output:
        path "${sample_id}.count.txt"

        script:
        """
        lines=\$(zcat ${fastq} | wc -l)
        cowpy "${sample_id} has \${lines} lines" > ${sample_id}.count.txt
        """
    ```

The cache from section 1 will be invalidated because the script body changed (more on that in section 5), so the tasks will re-run.

### 3.2. Run the pipeline

```bash
nextflow run sample_processing.nf -profile docker -resume
```

??? success "Command output"

    ```console
    executor >  local (6)
    [a3/5fe041] COUNT_LINES (3) | 5 of 5 ✔
    [bd/16c4b2] REPORT          | 1 of 1 ✔
    DEBUG: processed sample_001 (12 lines)

    DEBUG: processed sample_002 (12 lines)

    DEBUG: processed sample_003 (12 lines)

    DEBUG: processed sample_004 (12 lines)

    DEBUG: processed sample_005 (12 lines)
    ```

The `DEBUG:` lines come from inside the running process.
You can use this to print the value of any variable Nextflow has substituted, or any value computed by your script, at the moment the process actually sees it.
This is invaluable when the issue is "the process is doing something, but not what I expected".

### 3.3. Reset

`debug true` is a development tool. Remove it (and the `echo`) before continuing:

```bash
git checkout sample_processing.nf
```

And re-apply the section 1 fix once more.

### Takeaway

`debug true` plus an `echo` inside the script gives you live visibility into what a process actually receives and computes.
Reach for it when a process completes successfully but produces something you didn't expect.

### What's next?

Learn how to test workflow logic when running the real commands is slow or impossible, with `-stub-run`.

---

## 4. Iterate on logic with `-stub-run`

Sometimes the real work a process does is slow, expensive, or depends on software you don't have available locally.
The `stub:` directive lets you declare a "fake" command that produces files of the right shape without doing the real computation.
Run the workflow with `-stub-run` and Nextflow uses the stubs instead of the real scripts.

This is invaluable when you're iterating on downstream logic and don't want to wait for the upstream heavy lifting every time.

### 4.1. Look at the stub directive

The COUNT_LINES process in our pipeline already has one:

```groovy title="sample_processing.nf" hl_lines="1-4" linenums="23"
    stub:
    """
    echo "${sample_id} has 0 lines" > ${sample_id}.count.txt
    """
```

The stub produces a file with the same name as the real output, but with placeholder content and no `cowpy` dependency.

### 4.2. Run with `-stub-run` (no Docker required)

Notice that you don't need `-profile docker`. Stubs run on the host, so they bypass the container declaration entirely:

```bash
nextflow run sample_processing.nf -stub-run
```

??? success "Command output"

    ```console
    executor >  local (6)
    [af/649b28] COUNT_LINES (1) | 5 of 5 ✔
    [ff/8fddbe] REPORT          | 1 of 1 ✔
    ```

The pipeline runs end-to-end in seconds, producing the same output shape as the real run.
You can now iterate on REPORT (or any downstream change) without waiting on `cowpy`.

### Takeaway

`-stub-run` is the fastest way to validate a workflow change end-to-end.
Treat stubs as part of process design: every process you write should ship with a stub that produces the right output files.

### What's next?

Learn how to investigate why `-resume` didn't hit the cache when you expected it to, with `-dump-hashes`.

---

## 5. Debug cache invalidation with `-dump-hashes`

You make a small change to a pipeline, run it with `-resume`, and watch every task re-run.
Why?
Nextflow decides whether to reuse a cached task by hashing its inputs: script body, container, input values, and so on.
If any of those change, the hash changes and the cache misses.

The `-dump-hashes` flag writes every component of every task hash to `.nextflow.log`, so you can see exactly which input changed.

We'll run three experiments on the working pipeline, each making one small change and looking at the hash output to see what happened.

### 5.1. Establish a baseline

Make sure your pipeline is the working version from section 1 (output declarations and script both using `.count.txt`).
Run the pipeline once with `-dump-hashes` to populate the cache and the log:

```bash
nextflow run sample_processing.nf -profile docker -dump-hashes
```

Look at the cache hash entries that were written to `.nextflow.log`:

```bash
grep "cache hash" .nextflow.log | head
```

```console title="Output"
... [COUNT_LINES (4)] cache hash: 766733e20df3908928ff9e5f7cd2ba11; mode: STANDARD; entries:
... [COUNT_LINES (2)] cache hash: 5ba7de0fa92670ce4e414e3468df82ac; mode: STANDARD; entries:
... [COUNT_LINES (1)] cache hash: 1eef4a425bc5dc2ae4e06059a08404e4; mode: STANDARD; entries:
... [COUNT_LINES (5)] cache hash: 1d05640ac563902c089065adb578c616; mode: STANDARD; entries:
... [COUNT_LINES (3)] cache hash: 60025c26c50d23e7af16d9ca35b1b664; mode: STANDARD; entries:
```

Each task has a single hash, but it's computed from many components.
Look at the full entry for one task:

```bash
grep -A 18 "COUNT_LINES (1)" .nextflow.log
```

```console title="Output"
[COUNT_LINES (1)] cache hash: 1eef4a425bc5dc2ae4e06059a08404e4; mode: STANDARD; entries:
  ...UUID... [java.util.UUID] 326d1c67-b353-45b9-afda-38018e57a635
  ...      [java.lang.String] COUNT_LINES
  87b5860a [java.lang.String]     """
    lines=\$(zcat ${fastq} | wc -l)
    cowpy "${sample_id} has \${lines} lines" > ${sample_id}.count.txt
    """
  1b554354 [java.lang.String] community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273
  ...      [java.lang.String] sample_id
  e245add3 [java.lang.String] sample_001
  ...      [java.lang.String] fastq
  23747ba3 [nextflow.util.ArrayBag] [FileHolder(sourceObj:.../sample_001.fastq.gz, ...)]
  ...      [java.lang.String] $
  ...      [java.lang.Boolean] true
```

The components include the session UUID, the process name, the full script body, the container, each input variable name and value, and a few framework constants.
Note the script body in particular: its hash (`87b5860a...`) is computed from the exact text shown.
Any change to that text - including comments and whitespace - will produce a different hash.

!!! tip "JSON output"

    For longer pipelines, `-dump-hashes json` writes the entries in JSON format which is easier to diff with `jq` or a script.

### 5.2. Experiment 1 — a "harmless" comment

Add a single-line comment inside the COUNT_LINES script:

=== "After"

    ```groovy title="sample_processing.nf" hl_lines="3" linenums="17"
        script:
        """
        # Count lines in the gzipped fastq
        lines=\$(zcat ${fastq} | wc -l)
        cowpy "${sample_id} has \${lines} lines" > ${sample_id}.count.txt
        """
    ```

=== "Before"

    ```groovy title="sample_processing.nf" linenums="17"
        script:
        """
        lines=\$(zcat ${fastq} | wc -l)
        cowpy "${sample_id} has \${lines} lines" > ${sample_id}.count.txt
        """
    ```

Re-run with `-resume`:

```bash
nextflow run sample_processing.nf -profile docker -resume -dump-hashes
```

??? "Command output"

    ```console
    executor >  local (6)
    [da/470bba] COUNT_LINES (1) | 5 of 5 ✔
    [43/4659bf] REPORT          | 1 of 1 ✔
    ```

Despite using `-resume`, every COUNT_LINES task ran again. REPORT ran too, because its inputs depend on COUNT_LINES outputs.

Look at the new hash for COUNT_LINES (1):

```bash
grep -A 18 "COUNT_LINES (1)" .nextflow.log
```

```console title="Output"
[COUNT_LINES (1)] cache hash: 1a1e1e398f69e0ed797c7859c6c31b90; mode: STANDARD; entries:
  ...
  179bc8d7 [java.lang.String]     """
    # Count lines in the gzipped fastq
    lines=\$(zcat ${fastq} | wc -l)
    cowpy "${sample_id} has \${lines} lines" > ${sample_id}.count.txt
    """
  ...
```

The overall cache hash changed (`1eef4a42...` → `1a1e1e39...`).
The script hash changed too (`87b5860a...` → `179bc8d7...`) - because Nextflow hashes the script body literally, including the comment.

**Lesson:** any change inside `script:`, `stub:`, or `shell:` blocks invalidates the cache, even comments and whitespace.

Remove the comment before the next experiment.

### 5.3. Experiment 2 — a resource directive

Add a `memory` directive to COUNT_LINES:

=== "After"

    ```groovy title="sample_processing.nf" hl_lines="5" linenums="6"
    process COUNT_LINES {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
        publishDir params.outdir, mode: 'copy'
        memory '2.GB'
    ```

=== "Before"

    ```groovy title="sample_processing.nf" linenums="6"
    process COUNT_LINES {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
        publishDir params.outdir, mode: 'copy'
    ```

Re-run with `-resume`:

```bash
nextflow run sample_processing.nf -profile docker -resume
```

??? success "Command output"

    ```console
    [da/470bba] COUNT_LINES (1) | 5 of 5, cached: 5 ✔
    [43/4659bf] REPORT          | 1 of 1, cached: 1 ✔
    ```

Everything cached. The hashes are unchanged because resource directives like `memory`, `cpus`, and `time` are not part of the cache key.

**Lesson:** tuning resources never busts the cache.
That makes resource adjustment cheap, but it also means you can't force a re-run by changing memory or CPUs - you'd need to touch something that _is_ hashed (such as the script body) or pass `-resume <session_id>` from an earlier session.

Remove the `memory` directive before the next experiment.

### 5.4. Experiment 3 — a channel transformation

Change the workflow's `.map { }` to upper-case the sample ID:

=== "After"

    ```groovy title="sample_processing.nf" hl_lines="4" linenums="45"
    workflow {

        samples_ch = channel
            .fromPath(params.input)
            .splitCsv(header: true)
            .map { row -> [row.sample_id.toUpperCase(), file(row.fastq_path)] }

        counts_ch = COUNT_LINES(samples_ch)

        REPORT(counts_ch.collect())
    }
    ```

=== "Before"

    ```groovy title="sample_processing.nf" hl_lines="4" linenums="45"
    workflow {

        samples_ch = channel
            .fromPath(params.input)
            .splitCsv(header: true)
            .map { row -> [row.sample_id, file(row.fastq_path)] }

        counts_ch = COUNT_LINES(samples_ch)

        REPORT(counts_ch.collect())
    }
    ```

Re-run with `-resume`:

```bash
nextflow run sample_processing.nf -profile docker -resume -dump-hashes
```

??? "Command output"

    ```console
    executor >  local (6)
    [06/01db08] COUNT_LINES (1) | 5 of 5 ✔
    [d7/4c994c] REPORT          | 1 of 1 ✔
    ```

Every task ran again - including REPORT, which we didn't touch.
Look at one of the new COUNT_LINES hash entries:

```console title="Excerpt"
... [java.lang.String] sample_id
... [java.lang.String] SAMPLE_001       <-- was 'sample_001' before
... [java.lang.String] fastq
... [nextflow.util.ArrayBag] [FileHolder(...sample_001.fastq.gz...)]
```

The input value for `sample_id` is now `SAMPLE_001` rather than `sample_001`, so the COUNT_LINES task hash changed.
But why did REPORT re-run?
Because REPORT's input is the _files_ COUNT_LINES produced - and those files now live in different work directories (the new COUNT_LINES tasks).
A change upstream cascaded into a cache miss downstream, even though REPORT's own script and directives were untouched.

**Lesson:** when you trace a cache miss, look upstream too.
A change to a `.map { }` or input file can quietly invalidate every downstream process whose inputs depend on the changed task's outputs.

Restore the original `.map`:

```bash
git checkout sample_processing.nf
```

(and re-apply your section 1 fix once more).

### Takeaway

`-dump-hashes` shows the full input set Nextflow used to compute each task's cache key.
When a `-resume` re-runs more than you expected, find the affected task in `.nextflow.log` and check which component changed.
Comments, whitespace and any text inside the script block bust the cache. Resource directives don't. Upstream input changes cascade downstream.

### What's next?

Tie all the tools together into a systematic methodology.

---

## 6. A systematic debugging approach

With each technique covered individually, here is how they combine into a workflow you can apply to any pipeline failure.

### 6.1. The four-phase method

**Phase 1 — Parse first (seconds).**
Run `nextflow run workflow.nf -preview`.
Catches syntax errors, missing process definitions, and bad workflow structure before anything expensive runs.

**Phase 2 — Read the error (minutes).**
For runtime failures, the Nextflow error message names the failing process and includes the work directory.
Decide whether the error is structural (channel shape, missing software, wrong configuration) or process-internal (the command itself failed).

**Phase 3 — Investigate (minutes to hours).**
For process-internal failures, walk the work directory: `.command.sh`, `.command.err`, `.command.out`, `.exitcode`, `ls`.
For structural issues, use `.view()` on the relevant channel, or add `debug true` plus an `echo` inside a process to see what's actually flowing through.
If running the real commands is slow, switch to `-stub-run` while you iterate.

**Phase 4 — Fix and verify.**
Make the smallest change that addresses the root cause.
Re-run with `-resume`.
If the fix worked but tasks you didn't change also re-ran, use `-dump-hashes` to find out why.

### 6.2. A debugging profile

You can bake several of these tools into a profile so they're a single flag away:

```groovy title="nextflow.config (debug profile)"
profiles {
    debug {
        process {
            debug = true
            cleanup = false
            maxForks = 1
        }
    }
}
```

Then run with `-profile debug` whenever you need maximum visibility.
The `cleanup = false` keeps work directories around for inspection, and `maxForks = 1` makes parallel output easier to follow.

### Takeaway

Reach for tools in order of cost: `-preview` is free, work-directory inspection is the first step for any process failure, `debug true` and `.view()` give you channel and runtime visibility, `-stub-run` lets you iterate fast, and `-dump-hashes` is your last resort for cache mysteries.

### What's next?

Apply the toolkit to an unfamiliar pipeline.

---

## 7. Practical exercise

`buggy_workflow.nf` is a different pipeline that contains several intentional bugs covering all the categories from Part 1 and this lesson.
Use the four-phase method to fix it.

!!! exercise

    Run the workflow and start the debugging loop:

    ```bash
    nextflow run buggy_workflow.nf
    ```

    ??? failure "Command output"

        ```console
        N E X T F L O W   ~  version 25.10.4

        Launching `buggy_workflow.nf` [wise_ramanujan] DSL2 - revision: d51a8e83fd

        ERROR ~ Range [11, 12) out of bounds for length 11

         -- Check '.nextflow.log' file for details
        ```

        This cryptic error indicates a parsing problem in the `params{}` block. Apply Phase 1 (preview) first.

    **Suggested approach:**

    1. **Phase 1 — Parse first.** Use `-preview` to identify syntax issues. Fix them.
    2. **Phase 2 — Read each runtime error.** Decide whether the cause is structural or process-internal.
    3. **Phase 3 — Investigate.** For process failures, walk the work directory. For channel-shape problems, use `.view()` or `debug true`.
    4. **Phase 4 — Fix and verify.** After each fix, re-run with `-resume`.

    Stop when the workflow runs to completion. There are roughly 9 or 10 bugs depending on how you count.

    ??? solution

        The bugs in `buggy_workflow.nf`, in the order you'll typically encounter them:

        **Bug 1 — Trailing comma in output declaration**

        ```groovy linenums="21"
        output:
            path "${sample_id}_result.txt",  // remove the trailing comma
        ```

        **Bug 2 — Missing closing brace on `processFiles`**

        Add the missing `}` after the script block of `processFiles`.

        **Bug 3 — Variable name mismatch**

        ```groovy linenums="26"
        echo "Processing: ${sample}"     // should be ${sample_id}
        cat ${input_file} > ${sample}_result.txt  // should be ${sample_id}
        ```

        **Bug 4 — Undefined channel reference**

        ```groovy linenums="87"
        heavy_ch = heavyProcess(sample_ids)  // sample_ids doesn't exist; use input_ch
        ```

        **Bug 5 — Wrong channel shape for `processFiles`**

        ```groovy linenums="83"
        .map { row -> row.sample_id }  // processFiles expects a tuple
        // fix:
        .map { row -> [row.sample_id, file(row.fastq_path)] }
        ```

        **Bug 6 — `heavyProcess` now receives two-element tuples it doesn't want**

        ```groovy linenums="87"
        heavy_ch = heavyProcess(input_ch.map { it[0] })
        ```

        **Bug 7 — Unescaped Bash variable**

        ```groovy linenums="48"
        echo "Heavy computation \${i} for ${sample_id}"
        ```

        **Bug 8 — Unrealistic time limit**

        ```groovy linenums="36"
        time '100 s'   // not '1 ms'
        ```

        **Bug 9 — Output filename mismatch in `heavyProcess`**

        ```groovy linenums="49"
        done > ${sample_id}_heavy.txt   // not ${sample_id}.txt
        ```

        **Bug 10 — `handleFiles` is reading from `pwd` instead of an upstream channel**

        ```groovy linenums="88"
        file_ch = handleFiles(heavy_ch)
        ```

        Once all are fixed, the workflow runs to completion.

### Takeaway

The four-phase method scales: start with cheap checks (`-preview`), let error messages drive your next move, drop into forensic detail when needed, and verify each fix with `-resume`.

---

## Summary

You learned a sequence of debugging tools, each applied to the same small pipeline.

| Tool              | Best for                                                    |
| ----------------- | ----------------------------------------------------------- |
| Work directory    | Any process failure - the complete record of what ran       |
| `-preview`        | Catching syntax and structural issues without running tasks |
| `debug true`      | Seeing what a process actually receives at runtime          |
| `-stub-run`       | Iterating on workflow logic without real commands           |
| `-dump-hashes`    | Diagnosing unexpected cache invalidation                    |
| Four-phase method | A repeatable order in which to apply the tools              |

The single biggest debugging skill is matching the tool to the problem.
Most failures don't need every tool - they need the right one.

---

## What's next?

Return to the [menu of Side Quests](../index.md) or click the button in the bottom right of the page to move on to the next topic in the list.
