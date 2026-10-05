<!-- TODO(26.10): recapture console output on 26.10 release -->

# Part 3: Pipelines as Workflows

At the end of Part 2, your greetings pipeline was a complete, typed unit.
`GREETING_WORKFLOW` and `TRANSFORM_WORKFLOW` have typed `take:` and `emit:` blocks, and the entry workflow in `main.nf` wires them together and publishes the results.

Now you want to build on those results.
A **report** pipeline, developed independently of yours, takes a samplesheet of text files, counts what's in each one and writes a combined report.
It runs as a standalone pipeline, and you want to run it on your reversed greetings.

Nextflow can include a whole pipeline and call it like a workflow, so that two pipelines run as one.
Before that's possible, your greetings pipeline needs an interface at the pipeline level, the same way `GREETING_WORKFLOW` needed `take:` and `emit:` in Part 1.
In this part, you'll give it that interface, run it on its own, and call it from a new pipeline.

!!! tip "Starting from here?"

    If you're joining at this part, copy the solution from Part 2 into your working directory to use as your starting point:

    ```bash
    cd side-quests/composing_pipelines
    cp -r ../solutions/composing_pipelines/2/* .
    ```

### Learning goals

By the end of this part, you'll be able to:

- Explain how a pipeline's `params {}` and `output {}` blocks map onto a workflow's `take:` and `emit:`
- Load a `Channel` parameter from a samplesheet
- Run a pipeline from a subdirectory with `nextflow -C`
- Include a pipeline with `include { workflow as NAME }` and call it with a record of parameters
- Explain why an included pipeline's outputs are emitted to the caller rather than published

---

## 1. Give the pipeline its own directory

To include a pipeline, you need a pipeline to include and a separate pipeline to include it from.
Right now your greetings pipeline occupies the top level of your working directory, which is where the new top-level pipeline, the one that includes the others, will go.

### 1.1. Move the pipeline into `pipelines/greetings/`

Move everything that belongs to the greetings pipeline into its own directory, and clear out the results from earlier parts:

```bash
mkdir -p pipelines/greetings
mv main.nf nextflow.config types.nf modules workflows pipelines/greetings/
rm -rf results
```

```bash
tree pipelines
```

??? abstract "Directory contents"

    ```console
    pipelines
    └── greetings
        ├── main.nf
        ├── modules
        │   ├── reverse_text.nf
        │   ├── say_hello.nf
        │   ├── say_hello_upper.nf
        │   ├── timestamp_greeting.nf
        │   └── validate_name.nf
        ├── nextflow.config
        ├── types.nf
        └── workflows
            ├── greeting.nf
            └── transform.nf

    3 directories, 10 files
    ```

Nothing inside the pipeline changed.
Its includes are relative (`./workflows/greeting`, `../modules/say_hello`), so the pipeline works the same from its new location.

Your working directory also has a `data/` directory, which holds samplesheets of names, and an `extras/` directory, which holds the report pipeline.
You'll use both later in the course.

### 1.2. Compare the pipeline's interface with a workflow's

A named workflow has an interface: `take:` declares what goes in and `emit:` declares what comes out.
A pipeline has an interface too, at a larger scale:

| Interface | Named workflow             | Pipeline                                 |
| --------- | -------------------------- | ---------------------------------------- |
| Inputs    | `take:` block              | `params {}` block                        |
| Outputs   | `emit:` block              | `publish:` section and `output {}` block |
| Called by | `GREETING_WORKFLOW(names)` | `GREETINGS(record(names: ...))`          |

Now look at `pipelines/greetings/main.nf` with that table in mind:

```groovy title="pipelines/greetings/main.nf" linenums="8"
workflow {
    main:
    names = channel.of('Alice', 'Bob', 'Charlie')
```

The output half of the interface is already there: the `publish:` section and the `output {}` block from Part 2.
The input half is missing.
The names are hard-coded, so nobody, neither a user on the command line nor another pipeline, can give the pipeline different names.

### Takeaway

A pipeline's `params {}` block plays the role of `take:`, and its `publish:` section with the `output {}` block plays the role of `emit:`.
Your greetings pipeline has outputs but no inputs yet.

### What's next?

Give the pipeline a `params {}` block, so that it can receive its names from outside.

---

## 2. Make the pipeline includable

A pipeline that takes its inputs from a `params {}` block and returns its outputs through an `output {}` block can be run on its own and can be included by another pipeline.
You'll add the missing input, tighten up the output block, and check that the pipeline still runs on its own.

### 2.1. Add a `params {}` block

The names should arrive as a channel, because `GREETING_WORKFLOW` takes a `Channel<String>`.
A `Channel` parameter is loaded from a samplesheet, one record per row, as section 3 of the [Types and Records](../types_and_records/index.md) side quest showed.

Look at the samplesheet provided in `data/names.csv`:

```csv title="data/names.csv"
name
Alice
Bob
Charlie
```

Each row has a single `name` column.
Add a record type that matches it to `pipelines/greetings/types.nf`:

=== "After"

    ```groovy title="pipelines/greetings/types.nf" linenums="1" hl_lines="5-7"
    #!/usr/bin/env nextflow

    nextflow.enable.types = true

    record Person {
        name: String
    }

    record Greeting {
        name: String
        file: Path
    }
    ```

=== "Before"

    ```groovy title="pipelines/greetings/types.nf" linenums="1"
    #!/usr/bin/env nextflow

    nextflow.enable.types = true

    record Greeting {
        name: String
        file: Path
    }
    ```

Then, in `pipelines/greetings/main.nf`, include the new type, declare the `params {}` block, and read the names from it instead of hard-coding them:

=== "After"

    ```groovy title="pipelines/greetings/main.nf" linenums="1" hl_lines="7 9-11 15"
    #!/usr/bin/env nextflow

    nextflow.enable.types = true

    include { GREETING_WORKFLOW } from './workflows/greeting'
    include { TRANSFORM_WORKFLOW } from './workflows/transform'
    include { Greeting ; Person } from './types'

    params {
        names: Channel<Person>
    }

    workflow {
        main:
        names = params.names.map { person -> person.name }
    ```

=== "Before"

    ```groovy title="pipelines/greetings/main.nf" linenums="1" hl_lines="10"
    #!/usr/bin/env nextflow

    nextflow.enable.types = true

    include { GREETING_WORKFLOW } from './workflows/greeting'
    include { TRANSFORM_WORKFLOW } from './workflows/transform'

    workflow {
        main:
        names = channel.of('Alice', 'Bob', 'Charlie')
    ```

`names: Channel<Person>` has no default, so it's required: the pipeline won't start without it.
The entry workflow turns each `Person` record into a plain name with `map`, so `GREETING_WORKFLOW` and everything below it stay exactly as they were.
The `Greeting` type is included too, because the next step uses it.

!!! note "Why not `Channel<String>`?"

    A `Channel<String>` parameter is rejected with an error saying that the element type should be `Map`, `Record` or a record type.

### 2.2. Type the output block

The `output {}` block is now the outward-facing half of your pipeline's interface.
In Part 2 you typed `emit:` so that callers knew what they'd get back; do the same for the outputs:

=== "After"

    ```groovy title="pipelines/greetings/main.nf" linenums="30" hl_lines="2 5 8 11"
    output {
        greetings: Channel<Greeting> {
            path 'greetings'
        }
        timestamped: Channel<Greeting> {
            path 'timestamped'
        }
        upper: Channel<Greeting> {
            path 'upper'
        }
        reversed: Channel<Greeting> {
            path 'reversed'
        }
    }
    ```

=== "Before"

    ```groovy title="pipelines/greetings/main.nf" linenums="30" hl_lines="2 5 8 11"
    output {
        greetings {
            path 'greetings'
        }
        timestamped {
            path 'timestamped'
        }
        upper {
            path 'upper'
        }
        reversed {
            path 'reversed'
        }
    }
    ```

When another pipeline includes this one, these types are what it sees.
With them, `nextflow lint` can check the caller's code, for example reporting `Unrecognized property 'nme' for type Greeting` if the caller misspells a field.

Lint the pipeline:

```bash
nextflow lint pipelines/greetings
```

??? success "Command output"

    ```console
    Linting Nextflow code..
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
    Nextflow linting complete!
     ✅ 10 files had no errors
    ```

### 2.3. Run the pipeline on its own

Being includable doesn't take anything away: the greetings pipeline is still a pipeline you can run directly.
Run it from your working directory, passing the samplesheet as `--names`:

```bash
nextflow -C pipelines/greetings/nextflow.config run ./pipelines/greetings --names data/names.csv -output-dir results_greetings
```

??? success "Command output"

    ```console
     N E X T F L O W   ~  version 26.09.1-edge

    Launching `./pipelines/greetings/main.nf` [kickass_curie] revision: 7c8686ea3f

    executor >  local (15)
    [c0/b81049] GRE…DATE_NAME (validating Bob) | 3 of 3 ✔
    [3e/935045] GRE…Y_HELLO (greeting Charlie) | 3 of 3 ✔
    [80/f65bcc] GRE…ing timestamp to greeting) | 3 of 3 ✔
    [71/0d1c3d] TRA…tamped_Charlie-output.txt) | 3 of 3 ✔
    [da/c4773b] TRA…tamped_Charlie-output.txt) | 3 of 3 ✔

    Outputs:

      /workspaces/training/side-quests/composing_pipelines/results_greetings

      greetings:
        - {name: Bob, file: greetings/Bob-output.txt}
        - {name: Charlie, file: greetings/Charlie-output.txt}
        - {name: Alice, file: greetings/Alice-output.txt}

      timestamped:
        - {name: Bob, file: timestamped/timestamped_Bob-output.txt}
        - {name: Alice, file: timestamped/timestamped_Alice-output.txt}
        - {name: Charlie, file: timestamped/timestamped_Charlie-output.txt}

      upper:
        - {name: Bob, file: upper/UPPER-timestamped_Bob-output.txt}
        - {name: Alice, file: upper/UPPER-timestamped_Alice-output.txt}
        - {name: Charlie, file: upper/UPPER-timestamped_Charlie-output.txt}

      reversed:
        - {name: Bob, file: reversed/REVERSED-UPPER-timestamped_Bob-output.txt}
        - {name: Alice, file: reversed/REVERSED-UPPER-timestamped_Alice-output.txt}
        - {name: Charlie, file: reversed/REVERSED-UPPER-timestamped_Charlie-output.txt}
    ```

Nextflow loaded `data/names.csv` into a channel of `Person` records, and the pipeline produced the same greetings as before.

The command has three parts worth a closer look:

- **`-C pipelines/greetings/nextflow.config`**: by default, Nextflow loads the `nextflow.config` next to the pipeline script and also the `nextflow.config` in your launch directory, and merges them.
  In a moment, your launch directory gets a config of its own, for the top-level pipeline.
  `-C` tells Nextflow to use only the file you name.
  It's an option of `nextflow` itself, so it goes before `run`.
- **`./pipelines/greetings`**: a directory is enough; Nextflow runs the `main.nf` inside it.
- **`-output-dir results_greetings`**: the pipeline's config publishes to `results/`, which the top-level pipeline will use.
  This keeps the standalone results apart.

### Takeaway

Your greetings pipeline now has a complete interface: a `params {}` block for its input and a typed `output {}` block for its outputs.
A `Channel` parameter is loaded from a samplesheet when you run the pipeline directly, and you can still run the pipeline on its own with `nextflow -C`.

### What's next?

Write the top-level pipeline, which includes the greetings pipeline and calls it.

---

## 3. Include the pipeline in a new pipeline

The **top-level pipeline** is the one you launch: it lives at the top of your working directory and includes the other pipelines.
For now, it does only one thing: it calls the greetings pipeline.

### 3.1. Write the new `main.nf`

Create a new `main.nf` in your working directory:

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

Compare it with the old `main.nf` that called the two workflows:

- **`include { workflow as GREETINGS }`** includes the whole pipeline, meaning its `params {}` block, entry workflow and `output {}` block.
  The `workflow` keyword selects the unnamed entry workflow, and the alias gives it a name you can call.
- **`GREETINGS(record(names: params.names))`** calls it with a single record that has one field per parameter.
  A parameter with a default can be left out of the record; `names` has no default, so it must be there.
- **`greetings.timestamped`** and **`greetings.reversed`** are two of the outputs from the greetings pipeline's `output {}` block.
  As with a workflow that has several `emit:` outputs, the call returns them as fields of a record.
- **`#!groovy params { names: Channel<Person> }`** gives the top-level pipeline its own input, which it passes straight through.
  The `Person` type is included from the greetings pipeline so the two match.

The top-level pipeline only publishes two of the four outputs.
That's its choice to make, as you'll see when you run it.

### 3.2. Add a config for the top-level pipeline

The greetings pipeline's `nextflow.config` moved into `pipelines/greetings/`, so the top-level pipeline needs its own.
Create `nextflow.config` in your working directory:

```groovy title="nextflow.config" linenums="1"
outputDir = 'results'
workflow.output.mode = 'copy'
```

### 3.3. Run the top-level pipeline

```bash
nextflow run main.nf --names data/names.csv
```

??? success "Command output"

    ```console
     N E X T F L O W   ~  version 26.09.1-edge

    Launching `main.nf` [marvelous_mandelbrot] revision: e0f2fe97d8

    executor >  local (15)
    [b6/ea7c5d] GRE…TE_NAME (validating Alice) | 3 of 3 ✔
    [76/68c242] GRE…Y_HELLO (greeting Charlie) | 3 of 3 ✔
    [67/c2ad32] GRE…ing timestamp to greeting) | 3 of 3 ✔
    [28/bccb7e] GRE…estamped_Alice-output.txt) | 3 of 3 ✔
    [90/d2b0e0] GRE…estamped_Alice-output.txt) | 3 of 3 ✔

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
    ```

The same 15 tasks ran, but the console truncates their names.
Ask Nextflow for the full names of the tasks in the last run:

```bash
nextflow log last -f name
```

??? success "Command output"

    ```console
    GREETINGS:GREETING_WORKFLOW:VALIDATE_NAME (validating Bob)
    GREETINGS:GREETING_WORKFLOW:VALIDATE_NAME (validating Alice)
    GREETINGS:GREETING_WORKFLOW:VALIDATE_NAME (validating Charlie)
    GREETINGS:GREETING_WORKFLOW:SAY_HELLO (greeting Charlie)
    GREETINGS:GREETING_WORKFLOW:SAY_HELLO (greeting Alice)
    GREETINGS:GREETING_WORKFLOW:SAY_HELLO (greeting Bob)
    GREETINGS:GREETING_WORKFLOW:TIMESTAMP_GREETING (adding timestamp to greeting)
    GREETINGS:GREETING_WORKFLOW:TIMESTAMP_GREETING (adding timestamp to greeting)
    GREETINGS:GREETING_WORKFLOW:TIMESTAMP_GREETING (adding timestamp to greeting)
    GREETINGS:TRANSFORM_WORKFLOW:SAY_HELLO_UPPER (converting timestamped_Bob-output.txt)
    GREETINGS:TRANSFORM_WORKFLOW:SAY_HELLO_UPPER (converting timestamped_Charlie-output.txt)
    GREETINGS:TRANSFORM_WORKFLOW:SAY_HELLO_UPPER (converting timestamped_Alice-output.txt)
    GREETINGS:TRANSFORM_WORKFLOW:REVERSE_TEXT (reversing UPPER-timestamped_Bob-output.txt)
    GREETINGS:TRANSFORM_WORKFLOW:REVERSE_TEXT (reversing UPPER-timestamped_Charlie-output.txt)
    GREETINGS:TRANSFORM_WORKFLOW:REVERSE_TEXT (reversing UPPER-timestamped_Alice-output.txt)
    ```

Each name now shows all three levels of composition: the pipeline (`GREETINGS`), the workflow (`GREETING_WORKFLOW`) and the process (`VALIDATE_NAME`).
The pipeline's alias works like a workflow name: it scopes every process inside it.

Now look at what was published:

```bash
tree results
```

??? abstract "Directory contents"

    ```console
    results
    ├── reversed
    │   ├── REVERSED-UPPER-timestamped_Alice-output.txt
    │   ├── REVERSED-UPPER-timestamped_Bob-output.txt
    │   └── REVERSED-UPPER-timestamped_Charlie-output.txt
    └── timestamped
        ├── timestamped_Alice-output.txt
        ├── timestamped_Bob-output.txt
        └── timestamped_Charlie-output.txt

    2 directories, 6 files
    ```

There are no `greetings/` or `upper/` directories, even though the greetings pipeline's `output {}` block declares them.
When a pipeline is included, its outputs are **emitted** to the caller instead of being published, exactly like a workflow's `emit:`.
Only the top-level pipeline publishes anything, through its own `output {}` block, so it decides which results land in `results/` and where.

### Takeaway

`include { workflow as GREETINGS }` brings in a whole pipeline, which you call with a record of its parameters.
Its processes are scoped under the alias, and its outputs come back to you as channels rather than being published.
The pipeline itself didn't change, and it still runs on its own.

### What's next?

Look at the report pipeline that will consume the greetings.

---

## 4. Meet the report pipeline

The report pipeline is in `extras/report/`.
It was developed independently of your greetings pipeline, and it follows the same conventions: a typed `params {}` block for its inputs, and an `output {}` block for its outputs.

### 4.1. Read its interface

Open `extras/report/main.nf`:

```groovy title="extras/report/main.nf" linenums="1"
#!/usr/bin/env nextflow

nextflow.enable.types = true

include { COUNT_TEXT } from './modules/count_text'
include { SUMMARIZE } from './modules/summarize'
include { Sample } from './types'

params {
    // Samplesheet with one text file per sample (columns: id, file)
    input: Channel<Sample>

    // What to count in each file: 'chars' or 'words'
    metric: String = 'chars'

    // Title printed at the top of the report
    title: String = 'Text report'
}

workflow {
    main:
    counts = COUNT_TEXT(params.input, params.metric)
    summary = SUMMARIZE(counts.collect())

    publish:
    summary = summary
}

output {
    summary: Path {
        path 'report'
    }
}
```

You can tell how to use this pipeline from the `params {}` and `output {}` blocks alone, without opening its modules:

- **`input: Channel<Sample>`** is required: a samplesheet of text files, loaded into `Sample` records.
  `extras/report/types.nf` declares `Sample` with two fields, `id` and `file`.
- **`metric`** and **`title`** have defaults, so you can leave them out.
- **`summary: Path`** is a single file, the combined report.
  `counts.collect()` gathers every per-sample count before `SUMMARIZE` runs, so the report covers all samples at once.

Your greetings are `Greeting` records with a `name` and a `file`.
The report wants `Sample` records with an `id` and a `file`.
They're close, but not the same, which is normal for two pipelines developed separately.

### Takeaway

A pipeline with a typed `params {}` block and `output {}` block documents its own interface.
The report pipeline needs a channel of `Sample` records and returns a single report file.

---

## Takeaway

In this part, you turned your greetings pipeline into something another pipeline can include:

- **A pipeline-level interface**: the `params {}` block plays the role of `take:`, and the `output {}` block plays the role of `emit:`
- **`Channel` parameters**: loaded from a samplesheet into records when the pipeline runs on its own
- **Standalone runs**: `nextflow -C <config> run <dir>` runs a pipeline from a subdirectory using only its own config
- **Including a pipeline**: `include { workflow as GREETINGS }`, called with a record of parameters
- **Emitted, not published**: an included pipeline's outputs come back to the caller, which decides what to publish

---

## What's next?

You have two pipelines, each runnable on its own, and you know how to include one in another.
In Part 4, you'll connect them: first the traditional way, with two runs and a samplesheet in between, and then by including both in one pipeline and passing the greetings straight into the report.

[Continue to Part 4 :material-arrow-right:](04_chaining_pipelines.md){ .md-button .md-button--primary }
