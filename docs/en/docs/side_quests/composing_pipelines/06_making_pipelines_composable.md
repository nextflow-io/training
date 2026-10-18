<!-- TODO(26.10): recapture console output on 26.10 release -->

# Part 6: Making Pipelines Composable

At the end of Part 5, you could set every parameter of both pipelines, from the command line or from config, and you'd brought the report pipeline's process configuration into the composed run.
One problem remained: the report pipeline works on its own, but when your pipeline includes it, its banner disappears and its title reads `null`, even though you set `report.title` in config.

Nobody made an obvious mistake.
The report pipeline was developed and tested as a standalone pipeline, and it works that way.
The problem only shows up when another pipeline includes it, because a few habits that are harmless in a standalone pipeline depend on things that change when the pipeline is included.

In this part, you'll find those habits in the report pipeline, fix them, check that the pipeline works both ways, and finish with a checklist for writing pipelines that other pipelines can include.

!!! tip "Starting from here?"

    If you're joining at this part, copy the solution from Part 5 into your working directory, then run the greetings and report pipelines on their own, and the composed pipeline, once each:

    ```bash
    cd side-quests/composing_pipelines
    cp -r ../solutions/composing_pipelines/5/* .
    nextflow -C pipelines/greetings/nextflow.config run ./pipelines/greetings --names data/names.csv -output-dir results_greetings
    nextflow -C pipelines/report/nextflow.config run ./pipelines/report --input greetings_samplesheet.csv -output-dir results_report
    nextflow run main.nf
    ```

    The report pipeline in the Part 5 solution still contains the problem this part fixes.

### Learning goals

By the end of this part, you'll be able to:

- Explain what `params` refers to inside an included pipeline, and why reading it in a process breaks
- Explain the difference between `projectDir` and `moduleDir` when a pipeline is included
- Fix a pipeline so that it behaves the same on its own and when included
- Apply a checklist of practices for writing composable pipelines

---

## 1. Find the cause

Start from the symptoms, and narrow them down to the code responsible.

### 1.1. Compare the two reports

In Part 4, you ran the report pipeline on its own and published its report to `results_report/`.
Look at it next to the report from your composed pipeline:

```bash
cat results_report/report/report.txt
```

??? abstract "File contents"

    ```console title="results_report/report/report.txt"
    ==============================
        T E X T   R E P O R T
    ==============================
    Text report

    Alice: 36 chars
    Bob: 34 chars
    Charlie: 38 chars
    ```

```bash
cat results/report/report.txt
```

??? abstract "File contents"

    ```console title="results/report/report.txt"
    null

    Alice: 36 chars
    Bob: 34 chars
    Charlie: 38 chars
    ```

The counts are identical, so `COUNT_TEXT` behaves the same in both runs.
Both differences are at the top of the report, which `SUMMARIZE` writes: the banner is missing, and the title is `null`.

### 1.2. Lint the report pipeline

`nextflow lint` checks for some habits that cause trouble when a pipeline is included.
Lint the report pipeline:

```bash
nextflow lint pipelines/report
```

??? success "Command output"

    ```console
    Linting Nextflow code..
    Linting: pipelines/report/conf/modules.config
    Linting: pipelines/report/types.nf
    Linting: pipelines/report/modules/count_text.nf
    Linting: pipelines/report/modules/summarize.nf
    Linting: pipelines/report/nextflow.config
    Linting: pipelines/report/main.nf
    Warn  pipelines/report/modules/summarize.nf:15:14: The use of `projectDir` in a process is discouraged -- input files should be provided as process inputs
    │  15 |     banner=${projectDir}/assets/banner.txt
    ╰     |              ^^^^^^^^^^


    Nextflow linting complete!
     ⚠️  1 file had 1 warning
     ✅ 6 files had no errors
    ```

The linter flags `projectDir` in the `SUMMARIZE` process, the same process that writes the top of the report.

### 1.3. Read the summary process

Open `pipelines/report/modules/summarize.nf`:

```groovy title="pipelines/report/modules/summarize.nf" linenums="6" hl_lines="10 14"
process SUMMARIZE {
    input:
    counts: Bag<Path>

    output:
    file('report.txt')

    script:
    """
    banner=${projectDir}/assets/banner.txt
    if [ -f \$banner ]; then
        cat \$banner > report.txt
    fi
    echo "${params.title}" >> report.txt
    echo "" >> report.txt
    sort --parallel=${task.cpus} ${counts.join(' ')} >> report.txt
    """
}
```

Two lines reach outside the process for something it didn't receive as an input.

**`projectDir`** is the directory of the top-level pipeline, the one you launched.
When you run the report pipeline on its own, that's `pipelines/report/`, and the banner is at `pipelines/report/assets/banner.txt`.
When you run your composed pipeline, `projectDir` is the directory of the top-level `main.nf`, which here is also your working directory, and there's no `assets/banner.txt` there.
The banner is optional branding, so the script checks for it with `if [ -f ... ]` and leaves it out when the file is missing.
In the composed run the file isn't found, so the banner is skipped without any error.

**`params.title`** reads the run's parameters, not the report pipeline's.
When the report pipeline runs on its own, they're the same thing, so it works.
When it's included, only its entry workflow receives the report's own parameters, from the record passed in the call.
Everywhere else, including inside its processes, `params` means the top-level pipeline's parameters: `greetings` and `report`.
There's no `params.title` there, so Nextflow warns about an undefined parameter, and the script prints `null`.

Compare this with `metric`, which you changed successfully in Part 5.
The entry workflow reads `params.metric` and passes it to `COUNT_TEXT` as an input, so `COUNT_TEXT` gets the right value however the pipeline is run.

### Takeaway

Inside an included pipeline, `params` means the report's own parameters only in the entry workflow, and `projectDir` always points at the top-level pipeline.
A process that reads either one works on its own and breaks when included.

### What's next?

Make `SUMMARIZE` take the title and the banner as inputs, like `COUNT_TEXT` does with `metric`.

---

## 2. Fix the report pipeline

The fix is the same for both problems: everything the process needs arrives as an input, and only the entry workflow reads parameters and locates files.

### 2.1. Pass the title and the banner as inputs

Update `pipelines/report/modules/summarize.nf`:

=== "After"

    ```groovy title="pipelines/report/modules/summarize.nf" linenums="6" hl_lines="4-5 12-13"
    process SUMMARIZE {
        input:
        counts: Bag<Path>
        title: String
        banner: Path

        output:
        file('report.txt')

        script:
        """
        cat ${banner} > report.txt
        echo "${title}" >> report.txt
        echo "" >> report.txt
        sort --parallel=${task.cpus} ${counts.join(' ')} >> report.txt
        """
    }
    ```

=== "Before"

    ```groovy title="pipelines/report/modules/summarize.nf" linenums="6" hl_lines="10-14"
    process SUMMARIZE {
        input:
        counts: Bag<Path>

        output:
        file('report.txt')

        script:
        """
        banner=${projectDir}/assets/banner.txt
        if [ -f \$banner ]; then
            cat \$banner > report.txt
        fi
        echo "${params.title}" >> report.txt
        echo "" >> report.txt
        sort --parallel=${task.cpus} ${counts.join(' ')} >> report.txt
        """
    }
    ```

`title` is now a `String` input and `banner` a `Path` input.
Because `banner` is a `Path`, Nextflow stages the file into the task directory, so the script reads it directly.
A staged input also works on remote executors, such as cloud batch services, where a task can't read an absolute path on the launch machine like `#!groovy "${projectDir}/assets/banner.txt"`.
The existence check is gone too: the entry workflow now always supplies a banner, so a missing file fails the task instead of being skipped.

### 2.2. Supply them from the entry workflow

Now update the entry workflow in `pipelines/report/main.nf` to provide the two new inputs:

=== "After"

    ```groovy title="pipelines/report/main.nf" linenums="20" hl_lines="4-5"
    workflow {
        main:
        counts = COUNT_TEXT(params.input, params.metric)
        banner = file("${moduleDir}/assets/banner.txt")
        summary = SUMMARIZE(counts.collect(), params.title, banner)
    ```

=== "Before"

    ```groovy title="pipelines/report/main.nf" linenums="20" hl_lines="4"
    workflow {
        main:
        counts = COUNT_TEXT(params.input, params.metric)
        summary = SUMMARIZE(counts.collect())
    ```

The entry workflow is the one place where `params` always means the report's own parameters, so it's the right place to read `params.title`.

**`moduleDir`** is the directory of the script file that uses it, here `pipelines/report/`, no matter which pipeline was launched.
That makes it the right way for a pipeline to find files it ships with, such as the banner.

### 2.3. Lint the project

Lint both pipelines and the top-level pipeline:

```bash
nextflow lint main.nf pipelines
```

??? success "Command output"

    ```console
    Linting Nextflow code..
    Linting: main.nf
    Linting: pipelines/greetings/types.nf
    Linting: pipelines/greetings/modules/reverse_text.nf
    Linting: pipelines/greetings/modules/validate_name.nf
    Linting: pipelines/greetings/modules/say_hello.nf
    Linting: pipelines/greetings/modules/timestamp_greeting.nf
    Linting: pipelines/greetings/modules/say_hello_upper.nf
    Linting: pipelines/greetings/nextflow.config
    Linting: pipelines/greetings/workflows/transform.nf
    Linting: pipelines/greetings/workflows/greeting.nf
    Linting: pipelines/greetings/main.nf
    Linting: pipelines/report/conf/modules.config
    Linting: pipelines/report/types.nf
    Linting: pipelines/report/modules/count_text.nf
    Linting: pipelines/report/modules/summarize.nf
    Linting: pipelines/report/nextflow.config
    Linting: pipelines/report/main.nf
    Nextflow linting complete!
     ✅ 17 files had no errors
    ```

The `projectDir` warning is gone, and the whole project is lint-clean.

### Takeaway

Pass values into processes as inputs, read `params` only in the entry workflow, and locate a pipeline's own files with `moduleDir`.
The linter catches `projectDir` in a process, but not every use of `params` outside the entry workflow, so keep both habits in mind as you write.

### What's next?

Check that the fix works in the composed pipeline, and that it didn't break the report pipeline on its own.

---

## 3. Check both ways of running it

A composable pipeline has to work in both situations, so test both.

### 3.1. Run the composed pipeline

```bash
nextflow run main.nf -resume
```

??? success "Command output"

    ```console
     N E X T F L O W   ~  version 26.09.1-edge

    Launching `main.nf` [nasty_hilbert] revision: b973113816

    executor >  local (1)
    [29/4af209] GRE…_NAME (validating Charlie) | 3 of 3, cached: 3 ✔
    [1d/9d5962] GRE…Y_HELLO (greeting Charlie) | 3 of 3, cached: 3 ✔
    [38/ce3c7d] GRE…ing timestamp to greeting) | 3 of 3, cached: 3 ✔
    [e9/e29598] GRE…imestamped_Bob-output.txt) | 3 of 3, cached: 3 ✔
    [ef/66132d] GRE…tamped_Charlie-output.txt) | 3 of 3, cached: 3 ✔
    [aa/f1bff1] REPORT:COUNT_TEXT (Bob)        | 3 of 3, cached: 3 ✔
    [98/920867] REPORT:SUMMARIZE               | 1 of 1 ✔

    Outputs:

      /workspaces/training/side-quests/composing_pipelines/results

      ...

      summary: report/report.txt
    ```

The greetings and the counts came from the cache, and only `REPORT:SUMMARIZE` ran again, since its script changed.
There's no warning about an undefined parameter this time.

```bash
cat results/report/report.txt
```

??? abstract "File contents"

    ```console title="results/report/report.txt"
    ==============================
        T E X T   R E P O R T
    ==============================
    Greetings report

    Alice: 36 chars
    Bob: 34 chars
    Charlie: 38 chars
    ```

The banner is back, and the title is the one you set in Part 5, `Greetings report`.

### 3.2. Run the report pipeline on its own

Run the report pipeline on its own, on the samplesheet from Part 4:

```bash
nextflow -C pipelines/report/nextflow.config run ./pipelines/report --input greetings_samplesheet.csv -output-dir results_report
```

??? success "Command output"

    ```console
     N E X T F L O W   ~  version 26.09.1-edge

    Launching `./pipelines/report/main.nf` [nostalgic_meitner] revision: eb8e2580c6

    executor >  local (4)
    [37/2bc8ea] COUNT_TEXT (Alice) | 3 of 3 ✔
    [fa/cca921] SUMMARIZE          | 1 of 1 ✔

    Outputs:

      /workspaces/training/side-quests/composing_pipelines/results_report

      summary: report/report.txt
    ```

```bash
cat results_report/report/report.txt
```

??? abstract "File contents"

    ```console title="results_report/report/report.txt"
    ==============================
        T E X T   R E P O R T
    ==============================
    Text report

    Alice: 36 chars
    Bob: 34 chars
    Charlie: 38 chars
    ```

On its own, the report has its banner and its default title, `Text report`, just as before the fix.

The fix lives in your vendored copy of the report pipeline.
In a real project, you'd also upstream the fix to the report pipeline, so that the next version you vendor already includes it and nobody else who includes the pipeline hits the same problem.

### Takeaway

The fixed report pipeline produces the same report on its own as before, and the right report when included.
Test a pipeline both ways whenever you change something that affects how it's included.

### What's next?

Turn what you've learned across Parts 3 to 6 into a checklist.

---

## 4. A checklist for composable pipelines

When a pipeline is included, Nextflow brings in its code: its `main.nf` and everything that file includes.
Anything the pipeline relies on outside that code, such as its `nextflow.config`, its project directory or the run's parameters, belongs to the top-level pipeline instead.
The items below keep a pipeline working when that happens, and keep its interface easy to call.

### 4.1. What this course covered

1. **Declare every parameter in the `params {}` block, with a type.**
   The top-level pipeline imports them as a type and passes them in the call (Part 3, section 2.1; Part 5, section 1).
   A parameter set only in your `nextflow.config` doesn't reach your pipeline when it's included, because that file isn't loaded (Part 5, section 4.1).

2. **Type the `output {}` block.**
   When your pipeline is included, its outputs are emitted to the caller with these types, and `nextflow lint` checks the caller's code against them (Part 3, sections 2.2 and 3.3).
   The top-level `main.nf` in this course leaves its own outputs untyped because nothing includes it; type them as soon as something might.

3. **Read `params` only in the entry workflow and the `output {}` block.**
   Everywhere else, `params` means the top-level pipeline's parameters, so pass values to processes and workflows as inputs (sections 1.3 and 2 of this part).

4. **Locate your pipeline's own files with `moduleDir`, not `projectDir`.**
   `projectDir` is the top-level pipeline's directory (section 2.2 of this part).

5. **Keep process config in a separate file, and include it from your `nextflow.config`.**
   Your `nextflow.config` isn't loaded when the pipeline is included, but the separate file can be, by the top-level config (Part 5, sections 4.2 and 4.3).

6. **Mind selector collisions in shared process config.**
   `withName: 'SUMMARIZE'` matches both `SUMMARIZE` and `REPORT:SUMMARIZE`, so one file works on its own and when included.
   Once the top-level config includes that file, the selector can also match a same-named process from another included pipeline.
   When names could collide, the top-level config shouldn't include that file, and should set the same values with the alias instead, such as `withName: 'REPORT:SUMMARIZE'` (Part 5, section 4.3).

### 4.2. Also consider

The report pipeline didn't need these, but other pipelines often do:

- **Publish only through the `output {}` block, not with `publishDir`.**
  `publishDir` writes to a path fixed in the process, relative to the launch directory, outside the run's output directory and the top-level `output {}` block.
  Moving that publishing or turning it off means overriding `publishDir` for each process in the top-level config.
  Publishing through the pipeline's `output {}` block avoids this.
- **Don't rely on `bin/` or `lib/`.**
  Like `projectDir`, these directories are resolved for the top-level pipeline, so an included pipeline's scripts and classes there aren't found.
  A task calling one of its `bin/` scripts fails with "command not found", and code using one of its `lib/` classes reports it as not defined.
- **Set software dependencies (`container`, `conda`) and tool arguments in the process, in the shared process config from item 5, or pass them as inputs.**
  They don't belong in your `nextflow.config`, which isn't loaded when the pipeline is included.

None of these rules is absolute.
A top-level pipeline can always work around an included pipeline that breaks one, by redeclaring a parameter, copying a config setting or recreating a file.
Following them means it doesn't have to, so including your pipeline takes one `include` line and one call.

Run through the checklist on any pipeline you'd like others to be able to include, and run `nextflow lint` on it, which flags some of these habits for you.

### Takeaway

A composable pipeline depends only on its declared interface and on files it can find relative to its own code.
The checklist turns that into concrete habits, most of which make a pipeline cleaner even if nobody ever includes it.

---

## Takeaway

In this part, you made the report pipeline composable:

- **Diagnosing**: comparing standalone and composed output, and linting, pointed to `SUMMARIZE`
- **`params` outside the entry workflow**: means the top-level pipeline's parameters when included, so pass values as inputs
- **`projectDir` vs `moduleDir`**: `projectDir` belongs to the top-level pipeline, `moduleDir` to the file that uses it
- **Testing both ways**: a composable pipeline gives the right result on its own and when included
- **A checklist**: the practices you applied in this course, plus a few more, that let other pipelines include yours without workarounds

---

## What's next?

You've followed one project from processes, to workflows, to pipelines, to pipelines composed with other pipelines.
Head to the summary for a recap of the key patterns at each level.

[Continue to Summary :material-arrow-right:](summary.md){ .md-button .md-button--primary }
