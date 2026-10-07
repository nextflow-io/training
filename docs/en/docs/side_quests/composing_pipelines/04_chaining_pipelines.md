<!-- TODO(26.10): recapture console output on 26.10 release -->

# Part 4: Chaining Pipelines

At the end of Part 3, you had two pipelines that each run on their own.
Your greetings pipeline produces reversed greetings, and a new top-level pipeline includes it and publishes some of its outputs.
The report pipeline, developed separately, turns a samplesheet of text files into a combined report.

Now you want the report on your greetings.
In this part, you'll first connect the two pipelines the traditional way, with two separate runs and a samplesheet you write by hand, so you can see exactly what that costs.
Then you'll include the report pipeline next to the greetings pipeline and pass the greetings straight into it, so the whole chain runs as one pipeline.

!!! tip "Starting from here?"

    If you're joining at this part, copy the solution from Part 3 into your working directory, and run the greetings pipeline on its own once so that `results_greetings/` exists:

    ```bash
    cd side-quests/composing_pipelines
    cp -r ../solutions/composing_pipelines/3/* .
    nextflow -C pipelines/greetings/nextflow.config run ./pipelines/greetings --names data/names.csv -output-dir results_greetings
    ```

### Learning goals

By the end of this part, you'll be able to:

- Vendor a pipeline into your project so that you can include it
- Chain two pipelines through a samplesheet, and explain the costs of doing so, including a samplesheet that goes stale when the inputs change
- Connect one included pipeline's output channel to another included pipeline's input
- Adapt records between two pipelines whose types don't match exactly
- Run and resume two composed pipelines as a single run, so that a change to the input reruns only the affected tasks

---

## 1. Bring the report pipeline into your project

You can only include a pipeline from a local path, so the first step is to put a copy of the report pipeline inside your project.

### 1.1. Vendor the report pipeline

Copy the report pipeline next to your greetings pipeline:

```bash
cp -r extras/report pipelines/report
```

```bash
tree pipelines/report
```

??? abstract "Directory contents"

    ```console
    pipelines/report
    ├── assets
    │   └── banner.txt
    ├── main.nf
    ├── modules
    │   ├── count_text.nf
    │   └── summarize.nf
    ├── nextflow.config
    └── types.nf

    2 directories, 6 files
    ```

Copying a dependency's code into your own project is called **vendoring**.
Each vendored pipeline brings its own `modules/`, so two pipelines can depend on different versions of the same module without interfering with each other.

### Takeaway

The report pipeline now sits in `pipelines/report/`, a complete copy with its own modules, config and assets.

### What's next?

Before composing the two pipelines, connect them the way you would without composition.

---

## 2. Chain the pipelines by hand

Without composition, the only way to feed one pipeline's results into another is through files: run the first pipeline, find its outputs, describe them in a samplesheet, and run the second pipeline on that samplesheet.

### 2.1. Find the greetings pipeline's outputs

In Part 3, you ran the greetings pipeline on its own with `-output-dir results_greetings`.
Its reversed greetings are still there:

```bash
ls results_greetings/reversed
```

??? success "Command output"

    ```console
    REVERSED-UPPER-timestamped_Alice-output.txt
    REVERSED-UPPER-timestamped_Bob-output.txt
    REVERSED-UPPER-timestamped_Charlie-output.txt
    ```

### 2.2. Write a samplesheet for the report

The report pipeline wants a samplesheet with an `id` column and a `file` column.
Create `greetings_samplesheet.csv` in your working directory, with one row per greeting:

```csv title="greetings_samplesheet.csv" linenums="1"
id,file
Alice,results_greetings/reversed/REVERSED-UPPER-timestamped_Alice-output.txt
Bob,results_greetings/reversed/REVERSED-UPPER-timestamped_Bob-output.txt
Charlie,results_greetings/reversed/REVERSED-UPPER-timestamped_Charlie-output.txt
```

You had to copy each file name exactly, and you had to know where the greetings pipeline published its outputs.

### 2.3. Run the report pipeline

Run the report pipeline on its own, with its own config, on that samplesheet:

```bash
nextflow -C pipelines/report/nextflow.config run ./pipelines/report --input greetings_samplesheet.csv -output-dir results_report
```

??? success "Command output"

    ```console
     N E X T F L O W   ~  version 26.09.1-edge

    Launching `./pipelines/report/main.nf` [soggy_shirley] revision: 57103de971

    executor >  local (4)
    [ad/dc59a5] COUNT_TEXT (Charlie) | 3 of 3 ✔
    [d2/3fb0ff] SUMMARIZE            | 1 of 1 ✔

    Outputs:

      /workspaces/training/side-quests/composing_pipelines/results_report

      summary: report/report.txt
    ```

Look at the report:

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

The report opens with a banner and its default title, then lists the number of characters in each greeting.

### 2.4. Add a name

A new name, Diana, joins the input.
`data/more_names.csv` is the same samplesheet with one more row:

```bash
cat data/more_names.csv
```

??? abstract "File contents"

    ```csv title="data/more_names.csv"
    name
    Alice
    Bob
    Charlie
    Diana
    ```

Run the greetings pipeline on it, with `-resume` so that the three existing names can come from the cache:

```bash
nextflow -C pipelines/greetings/nextflow.config run ./pipelines/greetings --names data/more_names.csv -output-dir results_greetings -resume
```

??? success "Command output"

    ```console
     N E X T F L O W   ~  version 26.09.1-edge

    Launching `./pipelines/greetings/main.nf` [adoring_ochoa] revision: 7c8686ea3f

    executor >  local (20)
    [25/b4e134] GRE…TE_NAME (validating Alice) | 4 of 4 ✔
    [87/0c17f7] GRE…SAY_HELLO (greeting Alice) | 4 of 4 ✔
    [92/2765b3] GRE…ing timestamp to greeting) | 4 of 4 ✔
    [bc/f7e731] TRA…estamped_Diana-output.txt) | 4 of 4 ✔
    [49/46c487] TRA…estamped_Diana-output.txt) | 4 of 4 ✔

    Outputs:

      /workspaces/training/side-quests/composing_pipelines/results_greetings

      greetings:
        - {name: Bob, file: greetings/Bob-output.txt}
        - {name: Charlie, file: greetings/Charlie-output.txt}
        - {name: Diana, file: greetings/Diana-output.txt}
        - {name: Alice, file: greetings/Alice-output.txt}

      timestamped:
        - {name: Charlie, file: timestamped/timestamped_Charlie-output.txt}
        - {name: Alice, file: timestamped/timestamped_Alice-output.txt}
        - {name: Bob, file: timestamped/timestamped_Bob-output.txt}
        - {name: Diana, file: timestamped/timestamped_Diana-output.txt}

      upper:
        - {name: Alice, file: upper/UPPER-timestamped_Alice-output.txt}
        - {name: Charlie, file: upper/UPPER-timestamped_Charlie-output.txt}
        - {name: Bob, file: upper/UPPER-timestamped_Bob-output.txt}
        - {name: Diana, file: upper/UPPER-timestamped_Diana-output.txt}

      reversed:
        - {name: Alice, file: reversed/REVERSED-UPPER-timestamped_Alice-output.txt}
        - {name: Charlie, file: reversed/REVERSED-UPPER-timestamped_Charlie-output.txt}
        - {name: Bob, file: reversed/REVERSED-UPPER-timestamped_Bob-output.txt}
        - {name: Diana, file: reversed/REVERSED-UPPER-timestamped_Diana-output.txt}
    ```

All 20 tasks ran, not just Diana's five.
`-resume` restores the cache of the most recent run in this directory, whichever pipeline that was, and the last run was the report.
With two separate pipelines, keeping each one's `-resume` on its own history is up to you.

Diana's reversed greeting is now in `results_greetings/reversed/`.
Run the report again:

```bash
nextflow -C pipelines/report/nextflow.config run ./pipelines/report --input greetings_samplesheet.csv -output-dir results_report -resume
```

??? success "Command output"

    ```console
     N E X T F L O W   ~  version 26.09.1-edge

    Launching `./pipelines/report/main.nf` [insane_curry] revision: 57103de971

    executor >  local (4)
    [d3/8a28bb] COUNT_TEXT (Charlie) | 3 of 3 ✔
    [7b/ba1112] SUMMARIZE            | 1 of 1 ✔

    Outputs:

      /workspaces/training/side-quests/composing_pipelines/results_report

      summary: report/report.txt
    ```

The report's tasks ran again, even though its inputs are the same three greetings.
The greetings run has just rewritten those files in `results_greetings/`, and the report pipeline can't tell a rewritten file from a new one.
Look at the report:

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

Both runs succeeded, and Diana isn't in the report.
The samplesheet still lists three files, so the report pipeline never saw hers.
Nothing warned you: to include her, you'd have to notice and edit `greetings_samplesheet.csv` by hand again, with the exact path to her output.

This is the same problem as Part 1, one level up.
In Part 1, `transform.nf` picked up whatever `greeting.nf` had left in `results/`.
Here, the report pipeline reads its input from paths you copied into a samplesheet.
Either way, the second step depends on files an earlier run left behind, not on what the first step produced in this run.

Leave `greetings_samplesheet.csv` as it is; later parts use it with three names.

### 2.5. Count the costs

Chaining by hand gets the job done, and in a production setting you'd automate it: a script or scheduler that launches the first pipeline, polls until it finishes, checks that its outputs exist, writes the samplesheet, and launches the second pipeline.
Automated or not, the chain has the same weaknesses:

- **Two runs**: each pipeline has its own run, its own work directory history and its own `-resume`.
  Nothing tracks the chain as a whole: the report only sees files, so whenever the greetings pipeline republishes them, every report task runs again.
- **Waiting**: the report can't start until the greetings pipeline has finished completely, even though each greeting is ready long before the last one.
- **Path coupling**: the samplesheet depends on the greetings pipeline's output directory and file names, and on someone keeping it up to date.
  When the inputs change, it goes stale without an error, as it did with Diana.
  If the greetings pipeline renames an output or changes its layout, the chain breaks, and nothing tells you until the report pipeline fails or reads the wrong files.

There's also something worth keeping: each pipeline runs on its own, and you can re-run the report on old greetings results without touching the greetings pipeline.
Composition keeps that, as you'll see.

### Takeaway

Chaining pipelines through files works, but it splits one analysis into two runs, makes the second wait for the first, and ties the second to the first pipeline's output layout.
When the inputs change, nothing updates the samplesheet, so the report silently leaves out the new data.

### What's next?

Include the report pipeline next to the greetings pipeline, and connect them with a channel instead of a samplesheet.

---

## 3. Compose the two pipelines

In the composed version, the greetings go straight from one pipeline into the other, as a channel.
No samplesheet, no paths, no second launch.

### 3.1. Include the report pipeline and connect it

Update `main.nf` to include the report pipeline and pass it the reversed greetings:

=== "After"

    ```groovy title="main.nf" linenums="1" hl_lines="6 17-19 24 34-36"
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

        publish:
        timestamped = greetings.timestamped
        reversed = greetings.reversed
        summary = summary
    }

    output {
        timestamped {
            path 'timestamped'
        }
        reversed {
            path 'reversed'
        }
        summary {
            path 'report'
        }
    }
    ```

=== "Before"

    ```groovy title="main.nf" linenums="1"
    #!/usr/bin/env nextflow

    nextflow.enable.types = true

    include { workflow as GREETINGS } from './pipelines/greetings'
    include { Person } from './pipelines/greetings/types'

    params {
        names: Channel<Person>
    }

    workflow {
        main:
        greetings = GREETINGS(record(names: params.names))

        publish:
        timestamped = greetings.timestamped
        reversed = greetings.reversed
    }

    output {
        timestamped {
            path 'timestamped'
        }
        reversed {
            path 'reversed'
        }
    }
    ```

Three things changed:

- **`#!groovy include { workflow as REPORT }`** includes the report pipeline, the same way you included the greetings pipeline.
- **`#!groovy greetings.reversed.map { ... }`** reshapes each `Greeting` into the record the report expects: `name` becomes `id`, and `file` stays `file`.
  This one line replaces the samplesheet you wrote by hand.
  It names fields, not paths, so it doesn't care where or under what file names anything is published.
- **`REPORT(record(input: samples))`** passes the channel as the report's `input` parameter.
  When the report pipeline runs on its own, `input` is loaded from a samplesheet; here it receives a live channel instead.
  The report pipeline doesn't need to know the difference.

The report pipeline has a single output, so the call returns it directly, and the top-level pipeline publishes it as `summary`.

### 3.2. Run the composed pipeline

```bash
nextflow run main.nf --names data/names.csv
```

??? success "Command output"

    ```console
     N E X T F L O W   ~  version 26.09.1-edge

    Launching `main.nf` [curious_easley] revision: d3220442cd

    executor >  local (19)
    [bd/4e703e] GRE…TE_NAME (validating Alice) | 3 of 3 ✔
    [3c/8b6b53] GRE…W:SAY_HELLO (greeting Bob) | 3 of 3 ✔
    [2d/cfb439] GRE…ing timestamp to greeting) | 3 of 3 ✔
    [26/ec242e] GRE…estamped_Alice-output.txt) | 3 of 3 ✔
    [e2/26ecf7] GRE…estamped_Alice-output.txt) | 3 of 3 ✔
    [02/637ce5] REPORT:COUNT_TEXT (Alice)      | 3 of 3 ✔
    [bb/cc8fb4] REPORT:SUMMARIZE               | 1 of 1 ✔

    Outputs:

      /workspaces/training/side-quests/composing_pipelines/results

      timestamped:
        - {name: Bob, file: timestamped/timestamped_Bob-output.txt}
        - {name: Charlie, file: timestamped/timestamped_Charlie-output.txt}
        - {name: Alice, file: timestamped/timestamped_Alice-output.txt}

      reversed:
        - {name: Bob, file: reversed/REVERSED-UPPER-timestamped_Bob-output.txt}
        - {name: Charlie, file: reversed/REVERSED-UPPER-timestamped_Charlie-output.txt}
        - {name: Alice, file: reversed/REVERSED-UPPER-timestamped_Alice-output.txt}

      summary: report/report.txt

    WARN: Access to undefined parameter `title` -- Initialise it to a default value eg. `params.title = some_value`
    ```

One command ran all 19 tasks from both pipelines in a single run.
The report's processes appear as `REPORT:COUNT_TEXT` and `REPORT:SUMMARIZE`, scoped under their alias just like the greetings pipeline's processes.

Because the two pipelines share one dataflow graph, nothing waits for a pipeline to finish.
`REPORT:COUNT_TEXT` can start on a greeting as soon as `REVERSE_TEXT` emits it, while other greetings are still being processed.
`REPORT:SUMMARIZE` is different: the report pipeline collects every count before summarizing, so it runs once, after the last count.
That fan-in is part of the report pipeline's own design, and composition keeps it.

### 3.3. Check the report

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

The counts match the report from section 2.3, so the data made it through.
But compare the top of the file with the report you produced by hand: the banner is gone, and the title reads `null` instead of `Text report`.
The warning at the end of the run points the same way: something in the report pipeline read a parameter called `title` that didn't exist.

The report pipeline works on its own, and the composed run didn't fail, so this is easy to miss.
Note it for now; you'll track it down in Part 6, once you've seen how parameters and configuration reach an included pipeline.

### Takeaway

Including both pipelines and connecting them with a channel replaces the samplesheet, the second launch and the wait.
A `map` that reshapes records is all the adapter you need when two pipelines' types differ.
Composition can also surface behavior a pipeline never showed on its own, like the report's missing banner and title.

### What's next?

Make the change that left Diana out of the report, this time in the composed pipeline.

---

## 4. Add a name to the composed pipeline

In section 2.4, adding Diana reran both pipelines in full and still left her out of the report.
The composed pipeline is a single run with a single cache, and no samplesheet in between.

### 4.1. Run on the new samplesheet with `-resume`

```bash
nextflow run main.nf --names data/more_names.csv -resume
```

??? success "Command output"

    ```console
     N E X T F L O W   ~  version 26.09.1-edge

    Launching `main.nf` [chaotic_venter] revision: d3220442cd

    executor >  local (7)
    [b1/f58667] GRE…TE_NAME (validating Diana) | 4 of 4, cached: 3 ✔
    [c4/1fbe04] GRE…SAY_HELLO (greeting Diana) | 4 of 4, cached: 3 ✔
    [76/d5b86b] GRE…ing timestamp to greeting) | 4 of 4, cached: 3 ✔
    [9c/446012] GRE…estamped_Diana-output.txt) | 4 of 4, cached: 3 ✔
    [a5/2ba01e] GRE…estamped_Diana-output.txt) | 4 of 4, cached: 3 ✔
    [47/1d53d9] REPORT:COUNT_TEXT (Diana)      | 4 of 4, cached: 3 ✔
    [14/60b114] REPORT:SUMMARIZE               | 1 of 1 ✔

    Outputs:

      /workspaces/training/side-quests/composing_pipelines/results

      timestamped:
        - {name: Charlie, file: timestamped/timestamped_Charlie-output.txt}
        - {name: Bob, file: timestamped/timestamped_Bob-output.txt}
        - {name: Alice, file: timestamped/timestamped_Alice-output.txt}
        - {name: Diana, file: timestamped/timestamped_Diana-output.txt}

      reversed:
        - {name: Alice, file: reversed/REVERSED-UPPER-timestamped_Alice-output.txt}
        - {name: Bob, file: reversed/REVERSED-UPPER-timestamped_Bob-output.txt}
        - {name: Charlie, file: reversed/REVERSED-UPPER-timestamped_Charlie-output.txt}
        - {name: Diana, file: reversed/REVERSED-UPPER-timestamped_Diana-output.txt}

      summary: report/report.txt

    WARN: Access to undefined parameter `title` -- Initialise it to a default value eg. `params.title = some_value`
    ```

Seven tasks ran: Diana's five greetings tasks, her `REPORT:COUNT_TEXT` task, and `REPORT:SUMMARIZE`, which runs again because the counts it collects now include hers.
Everything for Alice, Bob and Charlie, in both pipelines, came from the cache.

```bash
cat results/report/report.txt
```

??? abstract "File contents"

    ```console title="results/report/report.txt"
    null

    Alice: 36 chars
    Bob: 34 chars
    Charlie: 38 chars
    Diana: 36 chars
    ```

Diana is in the report, and you didn't edit anything by hand.
Her greeting reached the report pipeline through the channel, so there's no samplesheet to go stale and no output path to get wrong.

Composition didn't take anything away, either.
Both pipelines are unchanged, so you can still run either one on its own with `nextflow -C`, as you did in section 2.

### Takeaway

A composed pipeline is one run with one cache, so a change to the input reruns only the affected tasks, across both pipelines.
Each included pipeline still runs on its own.

---

## Takeaway

In this part, you connected two independently developed pipelines:

- **Vendoring**: a copy of the report pipeline in `pipelines/report/` makes it includable
- **Chaining by hand**: a samplesheet built from published paths connects two runs, at the cost of waiting, path coupling, separate resumes and a samplesheet that goes stale when the inputs change
- **Composing**: one included pipeline's output channel feeds another's `Channel` parameter, with a `map` to adapt the records
- **One run**: both pipelines appear in one run, with one `-resume` that reruns only what a change affects, and each still runs standalone

You also saw a first hint of trouble: the report's title and banner went missing when it was included.

---

## What's next?

So far you've passed the greetings pipeline a single parameter, by redeclaring it in the top-level pipeline, and left every report parameter at its default.
In Part 5, you'll give users control over every parameter of both pipelines, from the command line and from config, and find out which parts of each pipeline's configuration apply when it's included.

[Continue to Part 5 :material-arrow-right:](05_params_and_config.md){ .md-button .md-button--primary }
