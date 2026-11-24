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

The pipeline used throughout this lesson is `sample_processing.nf`.
It reads a sample manifest, counts lines in each sample's gzipped FASTQ, and writes a small report:

```groovy title="sample_processing.nf" linenums="1"
#!/usr/bin/env nextflow

process COUNT_LINES {

    container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

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

    main:
    samples_ch = channel
        .fromPath(params.input)
        .splitCsv(header: true)
        .map { row -> [row.sample_id, file(row.fastq_path)] }

    counts_ch = COUNT_LINES(samples_ch)
    report_ch = REPORT(counts_ch.collect())

    publish:
    report = report_ch
}

output {
    report {
        path '.'
    }
}
```

The input path comes from `params.input`, which is set in the side quest's `nextflow.config`.
The `publish:` section and `output` block publish the report to the `results/` directory set in that config.

The pipeline ships with one intentional bug, which you'll fix in section 1.
The later sections build on the fixed pipeline, and each one tells you how to undo its changes when you're done.

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
     N E X T F L O W   ~  version 26.10.7

    Launching `sample_processing.nf` [goofy_yonath] revision: d9d7021aa3

    executor >  local (5)
    [a5/e9b813] COUNT_LINES (1) | 0 of 5
    [-        ] REPORT          -
    WARN: Killing running tasks (4)
    ERROR ~ Error executing process > 'COUNT_LINES (1)'

    Caused by:
      Missing output file(s) `sample_001.count.txt` expected by process `COUNT_LINES (1)`


    Command executed:

      lines=$(zcat sample_001.fastq.gz | wc -l)
      cowpy "sample_001 has ${lines} lines" > sample_001_count.txt

    Command exit status:
      0

    Command output:
      (empty)

    Work dir:
      /workspaces/training/side-quests/debugging/work/a5/e9b8130f972624eddbabdf1cfd8555

    Container:
      community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273

    Tip: you can replicate the issue by changing to the process work dir and entering the command `bash .command.run`

     -- Check '.nextflow.log' file for details
    ```

The process exited cleanly (`Command exit status: 0`) but Nextflow couldn't find the output file it was promised.
The error already hints at the cause: the script writes `sample_001_count.txt`, but the process declares `sample_001.count.txt`.
The work directory confirms it.

### 1.2. Walk the work directory

The error message gives you the path.
Copy it from your own terminal output; yours will have a different hash.

Every failed task has the same set of hidden files.
This section goes through them in turn.

#### 1.2.1. `.command.sh`: the executed command

This is the exact script Nextflow ran, after variable substitution:

```bash
cat work/a5/e9b8130f972624eddbabdf1cfd8555/.command.sh
```

```console title="Output"
#!/bin/bash -ue
lines=$(zcat sample_001.fastq.gz | wc -l)
cowpy "sample_001 has ${lines} lines" > sample_001_count.txt
```

Nextflow has already substituted `sample_001` for `${sample_id}`, and the redirect target is `sample_001_count.txt`.
The script ran exactly as written.

#### 1.2.2. `.command.err` and `.command.out`: what the command printed

```bash
cat work/a5/e9b8130f972624eddbabdf1cfd8555/.command.err
```

```console title="Output"
(empty)
```

```bash
cat work/a5/e9b8130f972624eddbabdf1cfd8555/.command.out
```

```console title="Output"
(empty)
```

Nothing on either stream, because `cowpy` wrote its output to the redirected file.

#### 1.2.3. `.exitcode`: how the command exited

```bash
cat work/a5/e9b8130f972624eddbabdf1cfd8555/.exitcode
```

```console title="Output"
0
```

Zero, as the error message already reported.
The failure is on Nextflow's side, not the command's: the command succeeded, but it produced a file with a name Nextflow wasn't expecting.

Common exit codes to recognise:

- **0**: success; if the task still failed, the cause is usually an output mismatch
- **127**: command not found (software not installed in the container)
- **137**: process killed, usually because it exceeded its memory limit

#### 1.2.4. `ls` the directory: what files actually exist

```bash
ls work/a5/e9b8130f972624eddbabdf1cfd8555/
```

```console title="Output"
sample_001.fastq.gz  sample_001_count.txt
```

The file `sample_001_count.txt` exists, but the process expects `sample_001.count.txt`.

### 1.3. Fix the bug

The output declaration uses a `.` (dot) and the script uses an `_` (underscore).
Make them match by fixing the script and the stub block:

=== "After"

    ```groovy title="sample_processing.nf" hl_lines="14 19" linenums="3"
    process COUNT_LINES {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

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

    ```groovy title="sample_processing.nf" hl_lines="14 19" linenums="3"
    process COUNT_LINES {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

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
     N E X T F L O W   ~  version 26.10.7

    Launching `sample_processing.nf` [serene_kilby] revision: f195ff9bb4

    executor >  local (6)
    [5b/e1b51a] COUNT_LINES (2) | 5 of 5 ✔
    [bf/148f61] REPORT          | 1 of 1 ✔

    Outputs:

      /workspaces/training/side-quests/debugging/results

      report: report.txt
    ```

This run uses `-resume`.
With this example every task failed, so nothing is reused, but in real pipelines `-resume` lets you retry just the failed steps after a fix.

### Takeaway

The work directory contains the complete record of what the process did.
`.command.sh` shows the executed command, `.command.err` and `.command.out` show its streams, `.exitcode` reveals how it finished, and a plain `ls` shows what files exist.
For any process failure: read the error, find the work directory, and walk through these files until the cause is obvious.

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
     N E X T F L O W   ~  version 26.10.7

    Launching `sample_processing.nf` [crazy_magritte] revision: f195ff9bb4

    [-        ] COUNT_LINES -
    [-        ] REPORT      -
    ```

The dashes mean "this process would run, but nothing was executed".
The script parsed cleanly, the workflow is well-formed, and the two processes show up as expected.

### 2.2. Catch a syntax error before running

Now break the file on purpose.
Open `sample_processing.nf` and delete the closing brace at the end of the workflow block (the `}` on line 52):

=== "Broken"

    ```groovy title="sample_processing.nf" hl_lines="14" linenums="39"
    workflow {

        main:
        samples_ch = channel
            .fromPath(params.input)
            .splitCsv(header: true)
            .map { row -> [row.sample_id, file(row.fastq_path)] }

        counts_ch = COUNT_LINES(samples_ch)
        report_ch = REPORT(counts_ch.collect())

        publish:
        report = report_ch
        // missing closing brace
    ```

=== "Working"

    ```groovy title="sample_processing.nf" hl_lines="14" linenums="39"
    workflow {

        main:
        samples_ch = channel
            .fromPath(params.input)
            .splitCsv(header: true)
            .map { row -> [row.sample_id, file(row.fastq_path)] }

        counts_ch = COUNT_LINES(samples_ch)
        report_ch = REPORT(counts_ch.collect())

        publish:
        report = report_ch
    }
    ```

Run `-preview` again:

```bash
nextflow run sample_processing.nf -preview
```

??? failure "Command output"

    ```console
     N E X T F L O W   ~  version 26.10.7

    Launching `sample_processing.nf` [cranky_allen] revision: 42fe8ff05d

    Error sample_processing.nf:58:1: Unexpected end of file -- check for an unclosed brace, bracket, parenthesis, or string


    ERROR ~ Script compilation failed

     -- Check '.nextflow.log' file for details
    ```

The error points at line 58, the end of the file, not at line 52: the parser only discovers the brace is missing when it runs out of input.
The error is identical to what you'd see if you ran the workflow normally, but it took seconds and didn't pull a container or run any tasks.

Add the closing brace back before moving on.

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

Edit `sample_processing.nf` and add two lines to the `COUNT_LINES` process:

=== "After"

    ```groovy title="sample_processing.nf" hl_lines="4 16" linenums="3"
    process COUNT_LINES {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
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

    ```groovy title="sample_processing.nf" linenums="3"
    process COUNT_LINES {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

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

The script body changed, so the cached results from section 1 can't be reused and the tasks will re-run (section 5 explains why).

### 3.2. Run the pipeline

```bash
nextflow run sample_processing.nf -profile docker -resume
```

??? success "Command output"

    ```console
     N E X T F L O W   ~  version 26.10.7

    Launching `sample_processing.nf` [admiring_maxwell] revision: eb63926ea6

    executor >  local (6)
    [86/8bee05] COUNT_LINES (1) | 5 of 5 ✔
    [1c/927b1f] REPORT          | 1 of 1 ✔
    DEBUG: processed sample_003 (12 lines)
    DEBUG: processed sample_004 (12 lines)
    DEBUG: processed sample_002 (12 lines)
    DEBUG: processed sample_001 (12 lines)
    DEBUG: processed sample_005 (12 lines)

    Outputs:

      /workspaces/training/side-quests/debugging/results

      report: report.txt
    ```

The `DEBUG:` lines come from inside the running process.
You can use this to print the value of any variable Nextflow has substituted, or any value computed by your script, at the moment the process actually sees it.
Reach for it when a process does something, but not what you expected.

### 3.3. Reset

`debug true` is a development tool.
Remove the `debug true` line and the `echo` line before continuing.

### Takeaway

`debug true` plus an `echo` inside the script gives you live visibility into what a process actually receives and computes.
Use it when a process completes successfully but produces something you didn't expect.

### What's next?

Learn how to test workflow logic when running the real commands is slow or impossible, with `-stub-run`.

---

## 4. Iterate on logic with `-stub-run`

Sometimes the real work a process does is slow, expensive, or depends on software you don't have available locally.
The `stub:` section lets you declare a placeholder command that produces files of the right shape without doing the real computation.
Run the workflow with `-stub-run` and Nextflow uses the stubs instead of the real scripts.

This helps when you're iterating on downstream logic and don't want to wait for the upstream heavy lifting every time.

### 4.1. Look at the stub

The `COUNT_LINES` process in the pipeline already has one:

```groovy title="sample_processing.nf" linenums="19"
    stub:
    """
    echo "${sample_id} has 0 lines" > ${sample_id}.count.txt
    """
```

The stub produces a file with the same name as the real output, but with placeholder content and no `cowpy` dependency.

### 4.2. Run with `-stub-run`

This time, leave out `-profile docker`.
Without it, Docker isn't enabled, and the stub's `echo` command doesn't need any software that isn't already installed:

```bash
nextflow run sample_processing.nf -stub-run
```

??? success "Command output"

    ```console
     N E X T F L O W   ~  version 26.10.7

    Launching `sample_processing.nf` [cranky_banach] revision: f195ff9bb4

    executor >  local (6)
    [e1/8383ee] COUNT_LINES (3) | 5 of 5 ✔
    [5b/16ab71] REPORT          | 1 of 1 ✔

    Outputs:

      /workspaces/training/side-quests/debugging/results

      report: report.txt
    ```

The pipeline runs end-to-end in seconds, producing the same output shape as the real run.
You can now iterate on `REPORT` (or any downstream change) without waiting on `cowpy`.

### Takeaway

`-stub-run` is the fastest way to validate a workflow change end-to-end.
Treat stubs as part of process design: every process you write should ship with a stub that produces the right output files.

### What's next?

Learn how to investigate why `-resume` didn't hit the cache when you expected it to, with `-dump-hashes`.

---

## 5. Debug cache invalidation with `-dump-hashes`

You make a small change to a pipeline, run it with `-resume`, and every task re-runs.
Nextflow decides whether to reuse a cached task by hashing its inputs: script body, container, input values, and so on.
If any of those change, the hash changes and the cache misses.

The `-dump-hashes` flag reports every component of every task hash, so you can see exactly which input changed.
The entries print to the console during the run and are also written to `.nextflow.log`, which is easier to search.
The command outputs in this section leave out the console entries.

The three experiments below each make one small change and use the hash output to show what happened.

### 5.1. Establish a baseline

Make sure your pipeline is the working version from section 1, with the output declaration and the script both using `.count.txt`.
Run the pipeline once with `-dump-hashes` to populate the cache and the log:

```bash
nextflow run sample_processing.nf -profile docker -dump-hashes
```

Find the cache hash entries in `.nextflow.log`:

```bash
grep "cache hash" .nextflow.log
```

```console title="Output (shortened)"
... [COUNT_LINES (4)] cache hash: 10ecd4117bfed742f548bb92e8185040; mode: STANDARD; entries:
... [COUNT_LINES (2)] cache hash: 369cd022bcbb34e326758f19c35b0e8b; mode: STANDARD; entries:
... [COUNT_LINES (5)] cache hash: 8829bd0bde9129e43cd68d31ee5930e3; mode: STANDARD; entries:
... [COUNT_LINES (1)] cache hash: 3947cc02a00a6a485f7c52692effade3; mode: STANDARD; entries:
... [COUNT_LINES (3)] cache hash: ec790c5e3d54f8a4ed175c1536a99856; mode: STANDARD; entries:
... [REPORT] cache hash: b359be7dbd9338a5186fb2dda896da14; mode: STANDARD; entries:
```

Each task has a single hash, computed from many components.
Look at the full entry for one task:

```bash
grep -A 15 "COUNT_LINES (1)\] cache hash" .nextflow.log
```

```console title="Output (shortened)"
... [COUNT_LINES (1)] cache hash: 3947cc02a00a6a485f7c52692effade3; mode: STANDARD; entries:
  6b3923b3035efad5c457ce378889777d [java.util.UUID] c6ea828a-f131-443b-a5cd-75def6189f19
  7984c034f9fa679dc4c66aa6caa67e30 [java.lang.String] COUNT_LINES
  87b5860ad73d6368e4e4c7e4386387de [java.lang.String]     """
    lines=\$(zcat ${fastq} | wc -l)
    cowpy "${sample_id} has \${lines} lines" > ${sample_id}.count.txt
    """

  1b554354b5bbbfb2ae34d99f6f78914b [java.lang.String] community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273
  7ec7cbd71ff757f5fcdbaa760c9ce6de [java.lang.String] sample_id
  e245add314c044224bc7d8be9df83b10 [java.lang.String] sample_001
  12d4e14db590163dc1217292bb444f41 [java.lang.String] fastq
  b52295e79996a373c2bd17dcd09880cf [nextflow.util.ArrayBag] [FileHolder(sourceObj:/workspaces/training/side-quests/debugging/data/sample_001.fastq.gz, ...)]
  4f9d4b0d22865056c37fb6d9c2a04a67 [java.lang.String] $
  16fe7483905cce7a85670e43e4678877 [java.lang.Boolean] true
```

The components include the session ID, the process name, the full script body, the container, each input variable name and value, and a few framework constants.
The script body's hash (`87b5860a...`) is computed from the exact text shown, so any change to that text, including comments and whitespace, produces a different hash.

!!! tip "JSON output"

    For longer pipelines, `-dump-hashes json` writes the entries in JSON format, which is easier to compare with `jq` or a script.

### 5.2. Experiment 1: a "harmless" comment

Add a single-line comment inside the `COUNT_LINES` script:

=== "After"

    ```groovy title="sample_processing.nf" hl_lines="3" linenums="13"
        script:
        """
        # Count lines in the gzipped fastq
        lines=\$(zcat ${fastq} | wc -l)
        cowpy "${sample_id} has \${lines} lines" > ${sample_id}.count.txt
        """
    ```

=== "Before"

    ```groovy title="sample_processing.nf" linenums="13"
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

??? success "Command output"

    ```console
     N E X T F L O W   ~  version 26.10.7

    Launching `sample_processing.nf` [loquacious_edison] revision: de89b0603d

    executor >  local (6)
    [f3/d619eb] COUNT_LINES (3) | 5 of 5 ✔
    [2e/1786ec] REPORT          | 1 of 1 ✔

    Outputs:

      /workspaces/training/side-quests/debugging/results

      report: report.txt
    ```

Despite `-resume`, every `COUNT_LINES` task ran again.
`REPORT` ran too, because its inputs depend on `COUNT_LINES` outputs.

Look at the new hash entry for `COUNT_LINES (1)`:

```bash
grep -A 15 "COUNT_LINES (1)\] cache hash" .nextflow.log
```

```console title="Output (excerpt)"
... [COUNT_LINES (1)] cache hash: cfdcf25abce68338b2aaa86a87a6d410; mode: STANDARD; entries:
  ...
  179bc8d788ad65e3f789c86209f14d55 [java.lang.String]     """
    # Count lines in the gzipped fastq
    lines=\$(zcat ${fastq} | wc -l)
    cowpy "${sample_id} has \${lines} lines" > ${sample_id}.count.txt
    """
  ...
```

The overall cache hash changed (`3947cc02...` → `cfdcf25a...`).
So did the script hash (`87b5860a...` → `179bc8d7...`), because Nextflow hashes the script body literally, including the comment.

Any change inside a `script:`, `stub:`, or `shell:` block invalidates the cache, even comments and whitespace.

Remove the comment before the next experiment.

### 5.3. Experiment 2: a resource directive

Add a `memory` directive to `COUNT_LINES`:

=== "After"

    ```groovy title="sample_processing.nf" hl_lines="4" linenums="3"
    process COUNT_LINES {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
        memory '2.GB'
    ```

=== "Before"

    ```groovy title="sample_processing.nf" linenums="3"
    process COUNT_LINES {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
    ```

Re-run with `-resume`:

```bash
nextflow run sample_processing.nf -profile docker -resume -dump-hashes
```

??? success "Command output"

    ```console
     N E X T F L O W   ~  version 26.10.7

    Launching `sample_processing.nf` [hungry_bhabha] revision: 6c225249a5

    [00/723b86] COUNT_LINES (4) | 5 of 5, cached: 5 ✔
    [9e/02243a] REPORT          | 1 of 1, cached: 1 ✔

    Outputs:

      /workspaces/training/side-quests/debugging/results

      report: report.txt
    ```

Everything was cached.
The task hashes are unchanged because resource directives such as `memory`, `cpus`, and `time` are not part of the cache key.

Tuning resources never invalidates the cache, which makes resource adjustments cheap to test.
Failed tasks are never cached, so after raising the `memory` of a task that ran out of memory, `-resume` re-runs that task and anything downstream of it, and reuses everything else.

Remove the `memory` directive before the next experiment.

### 5.4. Experiment 3: a channel transformation

Change the workflow's `.map { }` to upper-case the sample ID:

=== "After"

    ```groovy title="sample_processing.nf" hl_lines="7" linenums="39"
    workflow {

        main:
        samples_ch = channel
            .fromPath(params.input)
            .splitCsv(header: true)
            .map { row -> [row.sample_id.toUpperCase(), file(row.fastq_path)] }

        counts_ch = COUNT_LINES(samples_ch)
        report_ch = REPORT(counts_ch.collect())

        publish:
        report = report_ch
    }
    ```

=== "Before"

    ```groovy title="sample_processing.nf" hl_lines="7" linenums="39"
    workflow {

        main:
        samples_ch = channel
            .fromPath(params.input)
            .splitCsv(header: true)
            .map { row -> [row.sample_id, file(row.fastq_path)] }

        counts_ch = COUNT_LINES(samples_ch)
        report_ch = REPORT(counts_ch.collect())

        publish:
        report = report_ch
    }
    ```

Re-run with `-resume`:

```bash
nextflow run sample_processing.nf -profile docker -resume -dump-hashes
```

??? success "Command output"

    ```console
     N E X T F L O W   ~  version 26.10.7

    Launching `sample_processing.nf` [special_hamilton] revision: 97c7b66d09

    executor >  local (6)
    [e4/b597cd] COUNT_LINES (2) | 5 of 5 ✔
    [6c/751e5e] REPORT          | 1 of 1 ✔

    Outputs:

      /workspaces/training/side-quests/debugging/results

      report: report.txt
    ```

Every task ran again, including `REPORT`, which you didn't touch.
Look at the input entries in one of the new `COUNT_LINES` hash entries:

```bash
grep -A 15 "COUNT_LINES (1)\] cache hash" .nextflow.log
```

```console title="Output (excerpt)"
... [COUNT_LINES (1)] cache hash: 814e19c4b67ffc429ff5f11b8788a86c; mode: STANDARD; entries:
  ...
  7ec7cbd71ff757f5fcdbaa760c9ce6de [java.lang.String] sample_id
  af673c57ffefa779e4d736c2b4617348 [java.lang.String] SAMPLE_001
  12d4e14db590163dc1217292bb444f41 [java.lang.String] fastq
  b52295e79996a373c2bd17dcd09880cf [nextflow.util.ArrayBag] [FileHolder(sourceObj:/workspaces/training/side-quests/debugging/data/sample_001.fastq.gz, ...)]
  ...
```

The input value for `sample_id` is now `SAMPLE_001` rather than `sample_001`, so the `COUNT_LINES` task hash changed.
`REPORT` re-ran because its input is the set of files `COUNT_LINES` produced, and those are now different files: they have new names (`SAMPLE_001.count.txt`) and live in new work directories.
A change upstream cascaded into a cache miss downstream, even though `REPORT`'s own script and directives were untouched.

When you trace a cache miss, look upstream too.
A change to a `.map { }` or an input file can invalidate every downstream process whose inputs depend on the changed task's outputs.

Restore the original `.map { }` before moving on.

### Takeaway

`-dump-hashes` shows the full input set Nextflow used to compute each task's cache key.
When `-resume` re-runs more than you expected, find the affected task in `.nextflow.log` and check which component changed.
Any text change inside the script block invalidates the cache, including comments and whitespace.
Resource directives don't.
Upstream input changes cascade downstream.

### What's next?

Tie all the tools together into a systematic methodology.

---

## 6. A systematic debugging approach

With each technique covered individually, here is how they combine into a workflow you can apply to any pipeline failure.

### 6.1. The four-phase method

**Phase 1: Parse first (seconds).**
Run `nextflow run workflow.nf -preview`.
This catches syntax errors, missing process definitions, and bad workflow structure before anything expensive runs.

**Phase 2: Read the error (minutes).**
For runtime failures, the Nextflow error message names the failing process and includes the work directory.
Decide whether the error is structural (channel shape, missing software, wrong configuration) or process-internal (the command itself failed).

**Phase 3: Investigate (minutes to hours).**
For process-internal failures, walk the work directory: `.command.sh`, `.command.err`, `.command.out`, `.exitcode`, `ls`.
For structural issues, use `.view()` on the relevant channel, or add `debug true` plus an `echo` inside a process to see what's actually flowing through.
If running the real commands is slow, switch to `-stub-run` while you iterate.

**Phase 4: Fix and verify.**
Make the smallest change that addresses the root cause.
Re-run with `-resume`.
If the fix worked but tasks you didn't change also re-ran, use `-dump-hashes` to find out why.

### 6.2. A debugging profile

You can bundle several of these tools into a profile so they're a single flag away:

```groovy title="nextflow.config (debug profile)" linenums="1"
profiles {
    debug {
        process {
            debug = true
            maxForks = 1
        }
    }
}
```

Then run with `-profile debug` whenever you need maximum visibility.
`debug = true` streams every process's output to your terminal, and `maxForks = 1` runs one task at a time so that output is easier to follow.

### Takeaway

Reach for tools in order of cost: `-preview` is nearly free, work-directory inspection is the first step for any process failure, `debug true` and `.view()` give you channel and runtime visibility, `-stub-run` lets you iterate fast, and `-dump-hashes` answers cache questions.

### What's next?

Apply the toolkit to an unfamiliar pipeline.

---

## 7. Practical exercise

`buggy_workflow.nf` is a different pipeline that contains several intentional bugs covering the categories from Part 1 and this lesson.
Use the four-phase method to fix it.

!!! exercise

    Run the workflow and start the debugging loop:

    ```bash
    nextflow run buggy_workflow.nf
    ```

    ??? failure "Command output"

        ```console
         N E X T F L O W   ~  version 26.10.7

        Launching `buggy_workflow.nf` [nauseous_tuckerman] revision: 6a296ff695

        Error buggy_workflow.nf:25:12: Unexpected input: '\n'
        │  25 |     script:
        ╰     |            ^

        ERROR ~ Script compilation failed

         -- Check '.nextflow.log' file for details
        ```

        The parser points at line 25 (`script:`), but the real culprit is just above it: the trailing comma after the `output:` declaration on line 23 leaves the parser expecting another output, so it fails when it reaches `script:`. This is the first of several syntax errors to work through.

    **Suggested approach:**

    1. **Phase 1: Parse first.** Use `-preview` to identify syntax issues. Fix them.
    2. **Phase 2: Read each runtime error.** Decide whether the cause is structural or process-internal.
    3. **Phase 3: Investigate.** For process failures, walk the work directory. For channel-shape problems, use `.view()` or `debug true`.
    4. **Phase 4: Fix and verify.** After each fix, re-run with `-resume`.

    Stop when the workflow runs to completion. There are 10 bugs in total.

    ??? solution
        The `buggy_workflow.nf` contains 10 distinct errors covering all major debugging categories. Here's a systematic breakdown of each error and how to fix it, in the order you actually encounter them on Nextflow 26.10. The compiler resolves the workflow in two passes: first it parses the syntax, then it statically checks that every variable is defined. So you clear the syntax errors first, then a batch of undefined-variable errors, before the workflow runs at all and the runtime errors begin.

        Let's start with those syntax errors:

        **Error 1: Syntax Error - Trailing Comma**
        ```groovy linenums="21"
        output:
            path "${sample_id}_result.txt",  // ERROR: Trailing comma
        ```
        **Fix:** Remove the trailing comma
        ```groovy linenums="21"
        output:
            path "${sample_id}_result.txt"
        ```

        With the comma gone, the parser runs to the end of the file looking for the brace that should close `processFiles` and reports `Unexpected end of file`.

        **Error 2: Syntax Error - Missing Closing Brace**
        ```groovy linenums="24"
        script:
        """
        echo "Processing: ${sample}"
        cat ${input_file} > ${sample}_result.txt
        """
        // ERROR: Missing closing brace for processFiles process
        ```
        **Fix:** Add the missing closing brace
        ```groovy linenums="29"
        """
        echo "Processing: ${sample_id}"
        cat ${input_file} > ${sample_id}_result.txt
        """
        }  // Add missing closing brace
        ```

        Now the syntax parses, so the static type checker runs. It reports every undefined variable at once, before the workflow runs:

        ```console
        Error buggy_workflow.nf:86:29: `sample_ids` is not defined
        Error buggy_workflow.nf:27:25: `sample` is not defined
        Error buggy_workflow.nf:28:27: `sample` is not defined
        Error buggy_workflow.nf:49:33: `i` is not defined
        ```

        These four lines correspond to three distinct bugs, Errors 3, 4 and 5 below. The last of them, `i`, is a Bash variable that the type checker can't tell apart from a Nextflow variable, so it surfaces here at compile time rather than as a runtime failure. Fix all three before re-running.

        **Error 3: Variable Name Error**
        ```groovy linenums="26"
        echo "Processing: ${sample}"     // ERROR: should be sample_id
        cat ${input_file} > ${sample}_result.txt  // ERROR: should be sample_id
        ```
        **Fix:** Use the correct input variable name
        ```groovy linenums="26"
        echo "Processing: ${sample_id}"
        cat ${input_file} > ${sample_id}_result.txt
        ```

        **Error 4: Undefined Variable Error**
        ```groovy linenums="87"
        heavy_ch = heavyProcess(sample_ids)  // ERROR: sample_ids undefined
        ```
        **Fix:** Use the correct channel
        ```groovy linenums="87"
        heavy_ch = heavyProcess(input_ch)
        ```

        **Error 5: Bash Variable Escaping Error**
        ```groovy linenums="48"
        echo "Heavy computation $i for ${sample_id}"  // ERROR: $i looks like an undefined Nextflow variable
        ```
        **Fix:** Escape the bash variable so Nextflow leaves it for the shell
        ```groovy linenums="48"
        echo "Heavy computation \${i} for ${sample_id}"
        ```

        With those resolved the workflow compiles and starts to run. The first runtime error comes from `processFiles`, which expects a tuple but is being fed a bare value: `Input tuple does not match tuple declaration in process 'processFiles' -- offending value: sample_003`.

        **Error 6: Channel Structure Error - Wrong Map Output**
        ```groovy linenums="83"
        .map { row -> row.sample_id }  // ERROR: processFiles expects tuple
        ```
        **Fix:** Return the tuple structure that processFiles expects
        ```groovy linenums="83"
        .map { row -> [row.sample_id, file(row.fastq_path)] }
        ```

        That fixes `processFiles`, but `input_ch` now emits a two-element tuple, and `heavyProcess` is still being handed the whole tuple where it expects a single value. The tuple gets rendered into the script as `[sample_005, /path/sample_005.fastq.gz]`, which breaks the Bash command with a syntax error and an exit status of 2.
        Depending on timing, Nextflow may instead report the time limit of Error 8 first; the `Command executed` section still shows the tuple rendered into the script.

        **Error 7: Bad channel structure for heavyProcess**
        ```groovy linenums="87"
        heavy_ch = heavyProcess(input_ch)  // ERROR: input_ch now emits a 2-element tuple; heavyProcess needs only the first element
        ```
        **Fix:** Pass just the sample IDs
        ```groovy linenums="87"
        heavy_ch = heavyProcess(input_ch.map { sample_id, fastq -> sample_id })
        ```

        Now `heavyProcess` runs, but hits its time limit: `Process exceeded running time limit (1ms)`, which points to its `time` directive:

        **Error 8: Resource Configuration Error**
        ```groovy linenums="36"
        time '1 ms'  // ERROR: Unrealistic time limit
        ```
        **Fix:** Increase to a realistic time limit
        ```groovy linenums="36"
        time '100 s'
        ```

        Next we have a `Missing output file(s)` error to resolve, because the script writes `${sample_id}.txt` but the output declaration expects `${sample_id}_heavy.txt`:

        **Error 9: Output File Name Mismatch**
        ```groovy linenums="49"
        done > ${sample_id}.txt  // ERROR: Wrong filename, should match output declaration
        ```
        **Fix:** Match the output declaration
        ```groovy linenums="49"
        done > ${sample_id}_heavy.txt
        ```

        The workflow now completes without an error, but the `files` output is empty: `handleFiles` never ran. Its input channel, `channel.fromPath("*.txt")`, matches no files in the launch directory, so the process is simply skipped rather than failing loudly.

        **Error 10: Wrong Channel Source**
        ```groovy linenums="88"
        file_ch = channel.fromPath("*.txt") // Error: attempting to take input from the pwd rather than a process
        handleFiles(file_ch)
        ```
        **Fix:** Take the output from the previous process
        ```groovy linenums="88"
        file_ch = handleFiles(heavy_ch)
        ```

        With that, the whole workflow runs end to end and all three outputs are populated.

        **Complete Corrected Workflow:**
        ```groovy linenums="1"
        #!/usr/bin/env nextflow

        /*
        * Buggy workflow for debugging exercises
        * This workflow contains several intentional bugs for learning purposes
        */

        params{
            // Parameters with missing validation
            input: Path = 'data/sample_data.csv'
            output: String = 'results'
        }

        /*
        * Process with input/output mismatch
        */
        process processFiles {

            input:
                tuple val(sample_id), path(input_file)

            output:
                path "${sample_id}_result.txt"

            script:
            """
            echo "Processing: ${sample_id}"
            cat ${input_file} > ${sample_id}_result.txt
            """
        }

        /*
        * Process with resource issues
        */
        process heavyProcess {

            time '100 s'

            input:
                val sample_id

            output:
                path "${sample_id}_heavy.txt"

            script:
            """
            # Simulate heavy computation
            for i in {1..10000}; do
                echo "Heavy computation \$i for ${sample_id}"
            done > ${sample_id}_heavy.txt
            """
        }

        /*
        * Process with file handling issues
        */
        process handleFiles {

            input:
                path input_file

            output:
                path "processed_${input_file}"

            script:
            """
            if [ -f "${input_file}" ]; then
                cp ${input_file} processed_${input_file}
            fi
            """
        }

        /*
        * Main workflow with channel issues
        */
        workflow {
            main:
            // Channel with incorrect usage
            input_ch = channel
                .fromPath(params.input)
                .splitCsv(header: true)
                .map { row -> [row.sample_id, file(row.fastq_path)] }

            processed_ch = processFiles(input_ch)

            heavy_ch = heavyProcess(input_ch.map { sample_id, fastq -> sample_id })

            file_ch = handleFiles(heavy_ch)

            publish:
            processed = processed_ch
            heavy = heavy_ch
            files = file_ch
        }

        output {
            processed {
                path 'processed'
            }
            heavy {
                path 'heavy'
            }
            files {
                path 'files'
            }
        }
        ```

**Error Categories Covered:**

- **Syntax errors**: Missing braces, trailing commas, undefined variables
- **Channel structure errors**: Wrong data shapes, undefined channels
- **Process errors**: Output file mismatches, variable escaping
- **Resource errors**: Unrealistic time limits

**Key Debugging Lessons:**

1. **Read error messages carefully** - they often point directly to the problem
2. **Use systematic approaches** - fix one error at a time and test with `-resume`
3. **Understand data flow** - channel structure errors are often the most subtle
4. **Check work directories** - when processes fail, the logs tell you exactly what went wrong

### Takeaway

The four-phase method scales: start with cheap checks (`-preview`), let error messages drive your next move, drop into forensic detail when needed, and verify each fix with `-resume`.

---

You've reached the end of Part 2.
For a recap and quick-reference cheat sheet, continue to the summary.

[Continue to the summary :material-arrow-right:](summary.md){ .md-button .md-button--primary }
