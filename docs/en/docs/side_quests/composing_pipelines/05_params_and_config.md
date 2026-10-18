<!-- TODO(26.10): recapture console output on 26.10 release -->

# Part 5: Parameters and Configuration

At the end of Part 4, one command ran both pipelines: the greetings went straight from the greetings pipeline into the report pipeline, in a single run with a single `-resume`.

But look at how parameters reach the two pipelines.
The top-level pipeline redeclares the greetings pipeline's one parameter, `names`, just to pass it through.
The report pipeline's `metric` and `title` can't be set at all, so they always take their defaults.
Redeclaring every parameter by hand might be manageable for two small pipelines, but a real pipeline can have dozens of parameters, and every one you forget is one your users can't set.

Configuration has a similar gap.
Each pipeline ships its own `nextflow.config`, and it's not obvious which of those settings apply when the pipeline is included by another one.

In this part, you'll give users control over every parameter of both pipelines, from the command line and from config, and find out what happens to each pipeline's own configuration.

!!! tip "Starting from here?"

    If you're joining at this part, copy the solution from Part 4 into your working directory, run the greetings pipeline on its own once so that the files listed in `greetings_samplesheet.csv` exist, then run the composed pipeline once to build the cache that the next sections reuse:

    ```bash
    cd side-quests/composing_pipelines
    cp -r ../solutions/composing_pipelines/4/* .
    nextflow -C pipelines/greetings/nextflow.config run ./pipelines/greetings --names data/names.csv -output-dir results_greetings
    nextflow run main.nf --names data/names.csv
    ```

### Learning goals

By the end of this part, you'll be able to:

- Import an included pipeline's `params {}` block as a record type
- Set an included pipeline's parameters from the command line and from config
- Explain where an included pipeline's defaults come from, and when missing parameters are reported
- Override a parameter from the calling code, and recognize when a user's value is silently ignored
- Explain why an included pipeline's `nextflow.config` doesn't apply, and reuse its process configuration with `includeConfig`
- Target an included pipeline's processes with `withName` selectors

---

## 1. Import the parameters of included pipelines

Instead of redeclaring the parameters of each included pipeline, you can import a pipeline's whole `params {}` block as a type, and declare one parameter per pipeline.

### 1.1. Import each pipeline's `params {}` block

Update `main.nf`:

=== "After"

    ```groovy title="main.nf" linenums="1" hl_lines="5-6 9-10 15 19"
    #!/usr/bin/env nextflow

    nextflow.enable.types = true

    include { params as GreetingsParams ; workflow as GREETINGS } from './pipelines/greetings'
    include { params as ReportParams ; workflow as REPORT } from './pipelines/report'

    params {
        greetings: GreetingsParams
        report: ReportParams
    }

    workflow {
        main:
        greetings = GREETINGS(params.greetings)

        // Reshape each greeting into the sample record the report expects
        samples = greetings.reversed.map { greeting -> record(id: greeting.name, file: greeting.file) }
        summary = REPORT(params.report + record(input: samples))
    ```

=== "Before"

    ```groovy title="main.nf" linenums="1" hl_lines="5-7 10 15 19"
    #!/usr/bin/env nextflow

    nextflow.enable.types = true

    include { workflow as GREETINGS } from './pipelines/greetings'
    include { workflow as REPORT } from './pipelines/report'
    include { Person } from './pipelines/greetings/types'

    params {
        names: Channel<Person>
    }

    workflow {
        main:
        greetings = GREETINGS(record(names: params.names))

        // Reshape each greeting into the sample record the report expects
        samples = greetings.reversed.map { greeting -> record(id: greeting.name, file: greeting.file) }
        summary = REPORT(record(input: samples))
    ```

Here's what changed:

- **`#!groovy include { params as GreetingsParams ; ... }`** imports the greetings pipeline's `params {}` block as a record type called `GreetingsParams`, with one field per parameter.
  The same `include` also brings in the pipeline itself.
- **`greetings: GreetingsParams`** and **`report: ReportParams`** declare one parameter per included pipeline.
  Each pipeline's parameters live in their own namespace, so two pipelines can both have a parameter called `input` without clashing.
- **`GREETINGS(params.greetings)`** passes the greetings parameters through unchanged.
- **`#!groovy params.report + record(input: samples)`** takes whatever report parameters the user set and adds the `input` channel.
  With `+`, the record on the right wins, so `input` always comes from the greetings, as in Part 4.

The `Person` include is gone, since the top-level pipeline no longer declares `names` itself.

### 1.2. Run without parameters

Run the pipeline without passing any parameters:

```bash
nextflow run main.nf -resume
```

??? failure "Command output"

    ```console
     N E X T F L O W   ~  version 26.09.1-edge

    Launching `main.nf` [loquacious_avogadro] revision: b973113816

    Parameter `names` of pipeline `GREETINGS` is required but no value was provided

     -- Check script 'main.nf' at line: 15 or see '.nextflow.log' file for more details
    ```

The imported types are **partial**: every field is optional, so Nextflow doesn't check them when the top-level pipeline launches.
The greetings pipeline checks its own parameters when it's called, so the error names the pipeline (`GREETINGS`), the missing parameter (`names`) and the line of the call.

### 1.3. Pass parameters on the command line

Each pipeline's parameters are nested under the name you declared, so you set them with a dotted name.
Pass the greetings samplesheet, and set one of the report's parameters, which you couldn't do at all in Part 4:

```bash
nextflow run main.nf -resume --greetings.names data/names.csv --report.metric words
```

??? success "Command output"

    ```console
     N E X T F L O W   ~  version 26.09.1-edge

    Launching `main.nf` [furious_fermi] revision: b973113816

    executor >  local (4)
    [bd/4e703e] GRE…TE_NAME (validating Alice) | 3 of 3, cached: 3 ✔
    [3c/8b6b53] GRE…W:SAY_HELLO (greeting Bob) | 3 of 3, cached: 3 ✔
    [2d/cfb439] GRE…ing timestamp to greeting) | 3 of 3, cached: 3 ✔
    [48/210fff] GRE…tamped_Charlie-output.txt) | 3 of 3, cached: 3 ✔
    [3a/2d138c] GRE…tamped_Charlie-output.txt) | 3 of 3, cached: 3 ✔
    [69/b00198] REPORT:COUNT_TEXT (Charlie)    | 3 of 3 ✔
    [1b/e17976] REPORT:SUMMARIZE               | 1 of 1 ✔

    Outputs:

      /workspaces/training/side-quests/composing_pipelines/results

      timestamped:
        - {name: Bob, file: timestamped/timestamped_Bob-output.txt}
        - {name: Alice, file: timestamped/timestamped_Alice-output.txt}
        - {name: Charlie, file: timestamped/timestamped_Charlie-output.txt}

      reversed:
        - {name: Bob, file: reversed/REVERSED-UPPER-timestamped_Bob-output.txt}
        - {name: Charlie, file: reversed/REVERSED-UPPER-timestamped_Charlie-output.txt}
        - {name: Alice, file: reversed/REVERSED-UPPER-timestamped_Alice-output.txt}

      summary: report/report.txt

    WARN: Access to undefined parameter `title` -- Initialise it to a default value eg. `params.title = some_value`
    ```

`--greetings.names` reached the greetings pipeline as its `names` parameter, and was loaded from the samplesheet as a channel of `Person` records, exactly as `--names` was when you ran the pipeline on its own.
The greetings came from the cache, since restructuring how parameters reach the pipelines didn't change any of their tasks.
Only the report tasks ran again, because `metric` changed.

```bash
cat results/report/report.txt
```

??? abstract "File contents"

    ```console title="results/report/report.txt"
    null

    Alice: 4 words
    Bob: 4 words
    Charlie: 4 words
    ```

The report now counts words instead of characters.
The title still reads `null`, which you'll come back to shortly.

### Takeaway

`include { params as XParams }` imports an included pipeline's parameters as a partial record type, so the top-level pipeline can declare them all with one line per pipeline.
Users set them as `--<name>.<param>`, the included pipeline applies its own defaults, and it reports any required parameter that's still missing when it's called.

### What's next?

Typing `--greetings.names data/names.csv` on every run gets old quickly.
Put the values you always use into config instead.

---

## 2. Set parameters in config

Parameters for included pipelines can go in the top-level pipeline's `nextflow.config`, nested in the same way as on the command line.

### 2.1. Add nested parameters to `nextflow.config`

Set the greetings samplesheet, and give the report a more fitting title:

=== "After"

    ```groovy title="nextflow.config" linenums="1" hl_lines="4-11"
    outputDir = 'results'
    workflow.output.mode = 'copy'

    params {
        greetings {
            names = 'data/names.csv'
        }
        report {
            title = 'Greetings report'
        }
    }
    ```

=== "Before"

    ```groovy title="nextflow.config" linenums="1"
    outputDir = 'results'
    workflow.output.mode = 'copy'
    ```

### 2.2. Run with config parameters

```bash
nextflow run main.nf -resume
```

??? success "Command output"

    ```console
     N E X T F L O W   ~  version 26.09.1-edge

    Launching `main.nf` [big_faggin] revision: b973113816

    [bd/4e703e] GRE…TE_NAME (validating Alice) | 3 of 3, cached: 3 ✔
    [e3/4a2cdb] GRE…SAY_HELLO (greeting Alice) | 3 of 3, cached: 3 ✔
    [2d/cfb439] GRE…ing timestamp to greeting) | 3 of 3, cached: 3 ✔
    [48/210fff] GRE…tamped_Charlie-output.txt) | 3 of 3, cached: 3 ✔
    [48/5400c9] GRE…imestamped_Bob-output.txt) | 3 of 3, cached: 3 ✔
    [48/5c389b] REPORT:COUNT_TEXT (Charlie)    | 3 of 3, cached: 3 ✔
    [bb/cc8fb4] REPORT:SUMMARIZE               | 1 of 1, cached: 1 ✔

    Outputs:

      /workspaces/training/side-quests/composing_pipelines/results

      ...

      summary: report/report.txt

    WARN: Access to undefined parameter `title` -- Initialise it to a default value eg. `params.title = some_value`
    ```

The pipeline ran without any command-line parameters.
Nothing set `report.metric` this time, so the report pipeline applied its own default, `'chars'`, when it was called, and the report's tasks came from the cache of the earlier character counts.
Now check the title:

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

The title is still `null`, and the same warning appeared.
`report.metric` came from the command line and `report.title` from config, but both are nested under `report` the same way.
Config parameters do reach included pipelines: `greetings.names` is set only in config, and the greetings pipeline found it.
Yet whatever prints the title isn't seeing `report.title`, which points to the report pipeline itself rather than to how you passed the title.
Part 6 tracks it down.

### Takeaway

Parameters for included pipelines can be set in config, nested under the same names as on the command line.

### What's next?

You've set parameters from the command line, from config and from the calling code.
When more than one of them sets the same parameter, find out which one wins.

---

## 3. Know which value wins

In `main.nf`, the report's `input` comes from the calling code: `params.report + record(input: samples)`.
A user can also set `--report.input` on the command line.

### 3.1. Try to override the report's input

Make a samplesheet with only Alice in it:

```bash
head -n 2 greetings_samplesheet.csv > one_sample.csv
cat one_sample.csv
```

??? abstract "File contents"

    ```console title="one_sample.csv"
    id,file
    Alice,results_greetings/reversed/REVERSED-UPPER-timestamped_Alice-output.txt
    ```

Pass it as the report's input:

```bash
nextflow run main.nf -resume --report.input one_sample.csv
```

??? success "Command output"

    ```console
     N E X T F L O W   ~  version 26.09.1-edge

    Launching `main.nf` [pensive_plateau] revision: b973113816

    [a0/748675] GRE…DATE_NAME (validating Bob) | 3 of 3, cached: 3 ✔
    [3c/8b6b53] GRE…W:SAY_HELLO (greeting Bob) | 3 of 3, cached: 3 ✔
    [2d/cfb439] GRE…ing timestamp to greeting) | 3 of 3, cached: 3 ✔
    [5b/0b7665] GRE…imestamped_Bob-output.txt) | 3 of 3, cached: 3 ✔
    [3a/2d138c] GRE…tamped_Charlie-output.txt) | 3 of 3, cached: 3 ✔
    [02/637ce5] REPORT:COUNT_TEXT (Alice)      | 3 of 3, cached: 3 ✔
    [bb/cc8fb4] REPORT:SUMMARIZE               | 1 of 1, cached: 1 ✔

    Outputs:

      /workspaces/training/side-quests/composing_pipelines/results

      ...

      summary: report/report.txt

    WARN: Access to undefined parameter `title` -- Initialise it to a default value eg. `params.title = some_value`
    ```

Every task was cached, and the report still covers all three greetings.
Nextflow loaded `one_sample.csv`, then `+ record(input: samples)` replaced it with the greetings channel, without any message.

!!! warning "Overridden parameters are dropped silently"

    When the calling code sets a parameter with `params.x + record(...)`, any value the user gave for that parameter is discarded without a warning.
    Nextflow still loads the user's value at launch, so a samplesheet that doesn't exist or doesn't match the record type fails the run, even though it would never have been used.

    If your pipeline overrides a parameter of an included pipeline, say so in its documentation, so that users don't expect `--report.input` to do anything.

### Takeaway

A value set by the calling code with `+` wins over anything the user sets on the command line or in config.
The user's value is dropped silently, so document the parameters your pipeline overrides.

### What's next?

Parameters cross pipeline boundaries in a predictable way.
Configuration is a different story.

---

## 4. Configure included pipelines

The report pipeline's config sets two CPUs for its summary step, in `pipelines/report/nextflow.config`:

```groovy title="pipelines/report/nextflow.config" linenums="1"
outputDir = 'results'
workflow.output.mode = 'copy'

process {
    withName: 'SUMMARIZE' {
        cpus = 2
    }
}
```

### 4.1. Check what the report's processes got

`nextflow log` can show the resources each task of the last run was given.
Filter it to the report's tasks:

```bash
nextflow log last -f name,cpus,memory | grep REPORT
```

??? success "Command output"

    ```console
    REPORT:COUNT_TEXT (Bob)	1	-
    REPORT:COUNT_TEXT (Charlie)	1	-
    REPORT:COUNT_TEXT (Alice)	1	-
    REPORT:SUMMARIZE	1	-
    ```

`REPORT:SUMMARIZE` ran with one CPU, not two.

When you include a pipeline, Nextflow includes its code, meaning its `main.nf` and everything that file includes.
It doesn't load the included pipeline's `nextflow.config`.
Only the top-level pipeline's config applies to the run, so every setting in the report's config was ignored, without a warning.

Some of those settings are global, like `outputDir`, and the top-level pipeline sets its own anyway.
Process settings are different: they were chosen for the report pipeline's processes, and you'd usually want to keep them.

### 4.2. Move the report's process config into its own file

The pattern for sharing process config is to keep it in a separate file that both the pipeline's own `nextflow.config` and a top-level config can load.

Create `pipelines/report/conf/modules.config`:

```bash
mkdir -p pipelines/report/conf
```

```groovy title="pipelines/report/conf/modules.config" linenums="1"
process {
    withName: 'SUMMARIZE' {
        cpus = 2
    }
}
```

Then replace the `process` block in `pipelines/report/nextflow.config` with an `includeConfig` of the new file:

=== "After"

    ```groovy title="pipelines/report/nextflow.config" linenums="1" hl_lines="4"
    outputDir = 'results'
    workflow.output.mode = 'copy'

    includeConfig 'conf/modules.config'
    ```

=== "Before"

    ```groovy title="pipelines/report/nextflow.config" linenums="1" hl_lines="4-8"
    outputDir = 'results'
    workflow.output.mode = 'copy'

    process {
        withName: 'SUMMARIZE' {
            cpus = 2
        }
    }
    ```

When the report pipeline runs on its own, nothing changes: its `nextflow.config` loads the same settings, just from a different file.

### 4.3. Include the process config from the top-level config

In your top-level `nextflow.config`, include the report's process config, and add a setting of your own for the report's processes:

=== "After"

    ```groovy title="nextflow.config" linenums="1" hl_lines="13-19"
    outputDir = 'results'
    workflow.output.mode = 'copy'

    params {
        greetings {
            names = 'data/names.csv'
        }
        report {
            title = 'Greetings report'
        }
    }

    includeConfig 'pipelines/report/conf/modules.config'

    process {
        withName: 'REPORT:.*' {
            memory = 1.GB
        }
    }
    ```

=== "Before"

    ```groovy title="nextflow.config" linenums="1"
    outputDir = 'results'
    workflow.output.mode = 'copy'

    params {
        greetings {
            names = 'data/names.csv'
        }
        report {
            title = 'Greetings report'
        }
    }
    ```

Two kinds of selector are at work here:

- **`withName: 'SUMMARIZE'`**, in the included file, matches the report's process both when the report runs on its own (as `SUMMARIZE`) and when it's included (as `REPORT:SUMMARIZE`).
  That's what lets one file serve both situations.
  A simple name matches every process with that name in the composed run, though, including a same-named process from another included pipeline.
- **`withName: 'REPORT:.*'`**, in your config, matches every process whose name starts with `REPORT:`, that is, every process of the included report pipeline, and nothing in the greetings pipeline.
  Use the alias in a selector to adjust an included pipeline's settings from the top-level config, without editing the vendored copy.

If another included pipeline also had a `SUMMARIZE` process, adding `withName: 'REPORT:SUMMARIZE'` to your config wouldn't be enough.
It would override the report's process, but the included file's `withName: 'SUMMARIZE'` would still match the other pipeline's process.
When names could collide, don't include that file: write its settings in the top-level config with the alias instead, such as `withName: 'REPORT:SUMMARIZE'`.
This course has no such collision, so including the file works.

### 4.4. Check the result

Run the pipeline without `-resume`.
Changing a task's resources doesn't invalidate its cache entry, so a resumed run would reuse the old tasks and their old settings.

```bash
nextflow run main.nf
```

??? success "Command output"

    ```console
     N E X T F L O W   ~  version 26.09.1-edge

    Launching `main.nf` [ridiculous_descartes] revision: b973113816

    executor >  local (19)
    [0c/cc5ab4] GRE…TE_NAME (validating Alice) | 3 of 3 ✔
    [49/8fd953] GRE…SAY_HELLO (greeting Alice) | 3 of 3 ✔
    [65/cd1bb4] GRE…ing timestamp to greeting) | 3 of 3 ✔
    [e9/e29598] GRE…imestamped_Bob-output.txt) | 3 of 3 ✔
    [3a/dd7b4d] GRE…imestamped_Bob-output.txt) | 3 of 3 ✔
    [aa/f1bff1] REPORT:COUNT_TEXT (Bob)        | 3 of 3 ✔
    [ed/2ef2ee] REPORT:SUMMARIZE               | 1 of 1 ✔

    Outputs:

      /workspaces/training/side-quests/composing_pipelines/results

      ...

      summary: report/report.txt

    WARN: Access to undefined parameter `title` -- Initialise it to a default value eg. `params.title = some_value`
    ```

```bash
nextflow log last -f name,cpus,memory
```

??? success "Command output"

    ```console
    GREETINGS:GREETING_WORKFLOW:VALIDATE_NAME (validating Bob)	1	-
    GREETINGS:GREETING_WORKFLOW:VALIDATE_NAME (validating Charlie)	1	-
    GREETINGS:GREETING_WORKFLOW:VALIDATE_NAME (validating Alice)	1	-
    GREETINGS:GREETING_WORKFLOW:SAY_HELLO (greeting Bob)	1	-
    GREETINGS:GREETING_WORKFLOW:SAY_HELLO (greeting Charlie)	1	-
    GREETINGS:GREETING_WORKFLOW:SAY_HELLO (greeting Alice)	1	-
    GREETINGS:GREETING_WORKFLOW:TIMESTAMP_GREETING (adding timestamp to greeting)	1	-
    GREETINGS:GREETING_WORKFLOW:TIMESTAMP_GREETING (adding timestamp to greeting)	1	-
    GREETINGS:GREETING_WORKFLOW:TIMESTAMP_GREETING (adding timestamp to greeting)	1	-
    GREETINGS:TRANSFORM_WORKFLOW:SAY_HELLO_UPPER (converting timestamped_Alice-output.txt)	1	-
    GREETINGS:TRANSFORM_WORKFLOW:SAY_HELLO_UPPER (converting timestamped_Charlie-output.txt)	1	-
    GREETINGS:TRANSFORM_WORKFLOW:SAY_HELLO_UPPER (converting timestamped_Bob-output.txt)	1	-
    GREETINGS:TRANSFORM_WORKFLOW:REVERSE_TEXT (reversing UPPER-timestamped_Alice-output.txt)	1	-
    GREETINGS:TRANSFORM_WORKFLOW:REVERSE_TEXT (reversing UPPER-timestamped_Charlie-output.txt)	1	-
    GREETINGS:TRANSFORM_WORKFLOW:REVERSE_TEXT (reversing UPPER-timestamped_Bob-output.txt)	1	-
    REPORT:COUNT_TEXT (Alice)	1	1 GB
    REPORT:COUNT_TEXT (Charlie)	1	1 GB
    REPORT:COUNT_TEXT (Bob)	1	1 GB
    REPORT:SUMMARIZE	2	1 GB
    ```

`REPORT:SUMMARIZE` now gets its two CPUs from the report's own process config, and every report process gets the memory setting from your `withName: 'REPORT:.*'` selector.
The greetings processes are untouched.

!!! note "Global settings belong to the top-level pipeline"

    Settings that apply to the whole run, such as `outputDir`, executors, profiles, reports and the `plugins` block, can't come from an included pipeline.
    If an included pipeline relies on a plugin or a profile, declare it again in the top-level config.

### Takeaway

An included pipeline's `nextflow.config` isn't loaded, so its settings silently don't apply.
Keep process config in a separate file, load it with `includeConfig` from both the pipeline's own config and the top-level config, and use the include alias in `withName` selectors to target one pipeline's processes.

---

## Takeaway

In this part, you took control of how parameters and configuration reach included pipelines:

- **Imported parameters**: `include { params as XParams }` declares all of a pipeline's parameters with one line, as a partial record type
- **Setting parameters**: `--<name>.<param>` on the command line, or a nested `params` block in config, while the included pipeline applies its own defaults
- **Overrides**: `params.x + record(...)` in the calling code wins, and silently drops any value the user gave
- **Config isolation**: an included pipeline's `nextflow.config` isn't loaded, so share process config through `includeConfig`
- **Scoped selectors**: `withName: 'REPORT:.*'` targets one included pipeline's processes, and `withName: 'REPORT:SUMMARIZE'` avoids matching a same-named process in another pipeline

One problem is still open: the report's title.
`report.title` is set, `report.metric` works, and yet the report prints `null`.

---

## What's next?

The report pipeline works perfectly on its own and misbehaves when included.
In Part 6, you'll find out why, fix it, and turn what you've learned into a checklist for writing pipelines that work both on their own and inside another pipeline.

[Continue to Part 6 :material-arrow-right:](06_making_pipelines_composable.md){ .md-button .md-button--primary }
