<!-- TODO(26.10): recapture console output on 26.10 release -->

# Part 1: Workflows of Workflows

You maintain a small greetings pipeline with two stages: one creates timestamped greetings for a list of names, and the other turns those greetings into uppercase and reversed versions.
Right now each stage is a standalone script.
To chain them, you run the first script, and the second one picks up whatever files the first one left in `results/`.
That only works as long as both scripts agree on where those files live and what they're called, and nothing in either script says so.

Nextflow lets you build a pipeline from named workflows, each with a declared set of inputs and outputs, and pass channels from one workflow to the next in a single run.
In this part, you'll turn each stage into a named workflow and compose them into one pipeline, replacing the hand-off through `results/` with a channel.
That pipeline carries through the rest of the course: in Part 2 you'll give its workflows typed interfaces, and in Parts 3-6 the whole pipeline becomes a building block that another pipeline includes.

### Learning goals

By the end of this part, you'll be able to:

- Turn a standalone workflow script into a named workflow with `take:` and `emit:` blocks
- Explain why a named workflow can't run on its own, and what the entry workflow is for
- Include named workflows in an entry workflow and pass one workflow's outputs to another
- Replace a hand-off through published files with a channel between workflows

---

## 0. Get started

If you haven't yet done so, open the training codespace as described on the [course overview](index.md) page.

#### Move into the project directory

Move into the directory where the files for this tutorial are located.

```bash
cd side-quests/composing_pipelines
```

You can set VSCode to focus on this directory:

```bash
code .
```

The editor opens with the project directory in focus.

#### Review the materials

You'll find a `modules` directory with process definitions, a `workflows` directory with two pre-written workflow scripts, and a `main.nf` file that you will progressively update:

```console title="Directory contents"
├── data/                        # Used from Part 3
├── extras/                      # Used from Part 3
├── main.nf
├── modules/
│   ├── reverse_text.nf          # Reverses text content
│   ├── say_hello.nf             # Creates a greeting (from Hello Nextflow)
│   ├── say_hello_upper.nf       # Converts to uppercase (from Hello Nextflow)
│   ├── timestamp_greeting.nf    # Adds timestamps to greetings
│   └── validate_name.nf         # Validates input names
├── nextflow.config
└── workflows/
    ├── greeting.nf              # Standalone greeting workflow (to be made composable)
    └── transform.nf             # Standalone transform workflow (to be made composable)
```

The `modules/` directory contains the individual process definitions, and the `workflows/` directory contains the two pre-written workflow scripts you will work with in this part.

#### Review the assignment

Your challenge is to turn the two workflow scripts into named workflows, then compose them in `main.nf`:

- A `GREETING_WORKFLOW` that validates names, creates greetings, and adds timestamps
- A `TRANSFORM_WORKFLOW` that converts text to uppercase and reverses it

The finished pipeline chains both workflows into a single data flow:

1. **Validate**: check that each name is well-formed
2. **Greet**: generate a greeting for each valid name
3. **Timestamp**: record when each greeting was created
4. **Uppercase**: convert the timestamped greeting to uppercase
5. **Reverse**: reverse the uppercased text

The first three steps belong to `GREETING_WORKFLOW`, and the last two belong to `TRANSFORM_WORKFLOW`.
Building them as separate, composable workflows lets you develop and test each stage independently before wiring them together.

#### Readiness checklist

Think you're ready to dive in?

- [ ] I understand the goal of this course and its prerequisites
- [ ] My codespace is up and running
- [ ] I've set my working directory appropriately
- [ ] I understand the assignment

If you can check all the boxes, you're good to go.

---

## 1. Add the greeting workflow to the pipeline

The greeting workflow validates names and generates timestamped greetings.

### 1.1. Review and run the greeting workflow

Open `workflows/greeting.nf` and take a look at the code:

```groovy title="workflows/greeting.nf" linenums="1"
#!/usr/bin/env nextflow

include { VALIDATE_NAME } from '../modules/validate_name'
include { SAY_HELLO } from '../modules/say_hello'
include { TIMESTAMP_GREETING } from '../modules/timestamp_greeting'

workflow {
    main:
    names_ch = channel.of('Alice', 'Bob', 'Charlie')

    // Chain processes: validate -> create greeting -> add timestamp
    validated_ch = VALIDATE_NAME(names_ch)
    greetings_ch = SAY_HELLO(validated_ch)
    timestamped_ch = TIMESTAMP_GREETING(greetings_ch)

    publish:
    greetings = greetings_ch
    timestamped = timestamped_ch
}

output {
    greetings {
    }
    timestamped {
    }
}
```

This is a complete, self-contained workflow with the same structure you saw in the 'Hello Nextflow' tutorial.
It hardcodes the input names, chains three processes, and publishes two outputs.

Run it to verify everything works:

```bash
nextflow run workflows/greeting.nf
```

??? success "Command output"

    ```console
     N E X T F L O W   ~  version 26.09.1-edge

    Launching `workflows/greeting.nf` [kickass_maxwell] revision: 04c0a33f79

    executor >  local (9)
    [c3/97e273] VALIDATE_NAME (validating Bob) | 3 of 3 ✔
    [86/7b5f97] SAY_HELLO (greeting Bob)       | 3 of 3 ✔
    [18/eefc42] TIM…ing timestamp to greeting) | 3 of 3 ✔

    Outputs:

      /workspaces/training/side-quests/composing_pipelines/results

      greetings:
        - Alice-output.txt
        - Charlie-output.txt
        - Bob-output.txt

      timestamped:
        - timestamped_Alice-output.txt
        - timestamped_Bob-output.txt
        - timestamped_Charlie-output.txt
    ```

To make it composable with other workflows, a few things need to change.

### 1.2. Make the workflow composable

To make a workflow composable, three things need to change: the workflow gets a name, inputs move to a `take:` block, and outputs move to an `emit:` block.
The `emit:` block replaces the standalone `publish:` and `output {}` blocks, which belong in the entry workflow instead.

The following sections cover these changes one by one.

#### 1.2.1. Name the workflow

Give the workflow a name so it can be imported from a parent workflow.

=== "After"

    ```groovy title="workflows/greeting.nf" linenums="7" hl_lines="1"
    workflow GREETING_WORKFLOW {
    ```

=== "Before"

    ```groovy title="workflows/greeting.nf" linenums="7" hl_lines="1"
    workflow {
    ```

With a name, the workflow can be imported into other scripts.

#### 1.2.2. Declare inputs with `take:`

Replace the hardcoded channel declaration with a `take:` block that declares what inputs the workflow expects.
The `take:` block goes before `main:`, and the `names_ch = channel.of(...)` line is removed.

=== "After"

    ```groovy title="workflows/greeting.nf" linenums="7" hl_lines="2 3"
    workflow GREETING_WORKFLOW {
        take:
        names_ch // Input channel with names

        main:
        // Chain processes: validate -> create greeting -> add timestamp
        validated_ch = VALIDATE_NAME(names_ch)
        greetings_ch = SAY_HELLO(validated_ch)
        timestamped_ch = TIMESTAMP_GREETING(greetings_ch)
    ```

=== "Before"

    ```groovy title="workflows/greeting.nf" linenums="7" hl_lines="3"
    workflow GREETING_WORKFLOW {
        main:
        names_ch = channel.of('Alice', 'Bob', 'Charlie')

        // Chain processes: validate -> create greeting -> add timestamp
        validated_ch = VALIDATE_NAME(names_ch)
        greetings_ch = SAY_HELLO(validated_ch)
        timestamped_ch = TIMESTAMP_GREETING(greetings_ch)
    ```

The `take:` block declares the channel by name only.
The parent workflow defines what goes into it.

#### 1.2.3. Declare outputs with `emit:`

Remove the `publish:` section and the `output {}` block, and add an `emit:` block that names the outputs.

=== "After"

    ```groovy title="workflows/greeting.nf" linenums="16" hl_lines="2 3 4"

        emit:
        greetings = greetings_ch // Original greetings
        timestamped = timestamped_ch // Timestamped greetings
    }
    ```

=== "Before"

    ```groovy title="workflows/greeting.nf" linenums="16" hl_lines="2 3 4 7 8 9 10 11 12"

        publish:
        greetings = greetings_ch
        timestamped = timestamped_ch
    }

    output {
        greetings {
        }
        timestamped {
        }
    }
    ```

The `emit:` block exposes named outputs that parent workflows can access via `GREETING_WORKFLOW.out.greetings` and `GREETING_WORKFLOW.out.timestamped`.

#### 1.2.4. Verify the result and test it

After all three changes, the complete file should look like this:

```groovy title="workflows/greeting.nf" linenums="1" hl_lines="7 8 9 17 18 19"
#!/usr/bin/env nextflow

include { VALIDATE_NAME } from '../modules/validate_name'
include { SAY_HELLO } from '../modules/say_hello'
include { TIMESTAMP_GREETING } from '../modules/timestamp_greeting'

workflow GREETING_WORKFLOW {
    take:
    names_ch // Input channel with names

    main:
    // Chain processes: validate -> create greeting -> add timestamp
    validated_ch = VALIDATE_NAME(names_ch)
    greetings_ch = SAY_HELLO(validated_ch)
    timestamped_ch = TIMESTAMP_GREETING(greetings_ch)

    emit:
    greetings = greetings_ch // Original greetings
    timestamped = timestamped_ch // Timestamped greetings
}
```

Now try running it directly:

```bash
nextflow run workflows/greeting.nf
```

??? failure "Command output"

    ```console
     N E X T F L O W   ~  version 26.09.1-edge

    Launching `workflows/greeting.nf` [curious_poisson] revision: 89ac88c6c2

    No entry workflow specified -- script must define an entry workflow, a single process or named workflow, or be a code snippet
    ```

This introduces a key concept: the **entry workflow**.
Nextflow uses an unnamed `workflow {}` block as the entry point when you run a script directly.
`GREETING_WORKFLOW` is named, so Nextflow doesn't know how to run it on its own.

That's intentional.
Composable workflows are designed to be called from an entry workflow, not run directly.
The solution is an entry workflow in `main.nf` that imports and calls `GREETING_WORKFLOW`.

### 1.3. Update and test the main workflow

Now update the main workflow to call the greeting workflow.

#### 1.3.1. Include the greeting workflow and call it

Add the `include` statement, update the workflow body to call `GREETING_WORKFLOW` and replace the `channel.empty()` placeholder in `publish:`:

=== "After"

    ```groovy title="main.nf" linenums="1" hl_lines="3 9 10 13"
    #!/usr/bin/env nextflow

    include { GREETING_WORKFLOW } from './workflows/greeting'

    workflow {
        main:
        names = channel.of('Alice', 'Bob', 'Charlie')

        // Run the greeting workflow
        GREETING_WORKFLOW(names)

        publish:
        greetings = GREETING_WORKFLOW.out.greetings
    }
    ```

=== "Before"

    ```groovy title="main.nf" linenums="1" hl_lines="8"
    #!/usr/bin/env nextflow

    workflow {
        main:
        names = channel.of('Alice', 'Bob', 'Charlie')

        publish:
        greetings = channel.empty()
    }
    ```

The entry workflow stays un-named so that Nextflow will use it as the pipeline entry point.

#### 1.3.2. Update the output block

Add a `path` directive to route published greetings into a `greetings/` subdirectory:

=== "After"

    ```groovy title="main.nf" linenums="16" hl_lines="3"
    output {
        greetings {
            path 'greetings'
        }
    }
    ```

=== "Before"

    ```groovy title="main.nf" linenums="16"
    output {
        greetings {
        }
    }
    ```

#### 1.3.3. Run the workflow

Run the workflow to test that it works:

```bash
nextflow run main.nf
```

??? success "Command output"

    ```console
     N E X T F L O W   ~  version 26.09.1-edge

    Launching `main.nf` [gloomy_legentil] revision: 6b0e98dd05

    executor >  local (9)
    [db/3e1397] GRE…_NAME (validating Charlie) | 3 of 3 ✔
    [be/2fc96a] GRE…W:SAY_HELLO (greeting Bob) | 3 of 3 ✔
    [ae/711cf6] GRE…ing timestamp to greeting) | 3 of 3 ✔

    Outputs:

      /workspaces/training/side-quests/composing_pipelines/results

      greetings:
        - greetings/Alice-output.txt
        - greetings/Charlie-output.txt
        - greetings/Bob-output.txt
    ```

??? abstract "New content added under `results/`"

    ```console
    results/
    └── greetings
        ├── Alice-output.txt
        ├── Bob-output.txt
        └── Charlie-output.txt
    ```

    The six files from the standalone run in section 1.1 (`Alice-output.txt`, `timestamped_Alice-output.txt`, and so on) are also in `results/`.
    Leave them there: section 2.1 reads the timestamped ones.

??? abstract "File contents"

    ```console title="results/greetings/Alice-output.txt"
    Hello, Alice!
    ```

The greeting files are published to `results/greetings/`.
The main workflow calls `GREETING_WORKFLOW` and wires its output directly to the `publish:` section.

`GREETING_WORKFLOW.out.timestamped` isn't wired to `publish:` here.
Starting in section 2, that channel becomes the input to `TRANSFORM_WORKFLOW` instead of a published output in its own right.
Contrast this with the standalone run in section 1.1, where `timestamped` had nothing else to feed into, so it was published directly.

### Takeaway

A named workflow declares its inputs in a `take:` block and its outputs in an `emit:` block, and it only runs when an entry workflow calls it.
The entry workflow in `main.nf` includes `GREETING_WORKFLOW`, passes it a channel of names, and reads its outputs through `GREETING_WORKFLOW.out`.

### What's next?

`GREETING_WORKFLOW` emits the timestamped greetings, but nothing in the pipeline consumes them yet.
The transform stage does, and at the moment it gets them from `results/`.

---

## 2. Add the transformation workflow to the pipeline

The transform workflow applies text transformations to the timestamped greetings.
Before making it composable, look at how it finds its input.

### 2.1. Review and run the workflow

Open `workflows/transform.nf` and take a look at the code:

```groovy title="workflows/transform.nf" linenums="1"
#!/usr/bin/env nextflow

include { SAY_HELLO_UPPER } from '../modules/say_hello_upper'
include { REVERSE_TEXT } from '../modules/reverse_text'

workflow {
    main:
    input_ch = channel.fromPath('results/timestamped_*.txt')

    // Apply transformations in sequence
    upper_ch = SAY_HELLO_UPPER(input_ch)
    reversed_ch = REVERSE_TEXT(upper_ch)

    publish:
    upper = upper_ch
    reversed = reversed_ch
}

output {
    upper {
    }
    reversed {
    }
}
```

This standalone workflow reads timestamped greeting files from the `results/` directory produced by `greeting.nf`, converts them to uppercase, then reverses the text.

The first line of `main:` is the only connection between the two stages: `#!groovy channel.fromPath('results/timestamped_*.txt')`.
It only finds anything if `greeting.nf` has already run, published to `results/`, and named its files `timestamped_*.txt`.
None of that is written down in `greeting.nf`, so a change there can break `transform.nf` without either script showing it.

Run it to verify it works with the greeting results from section 1.1:

```bash
nextflow run workflows/transform.nf
```

??? success "Command output"

    ```console
     N E X T F L O W   ~  version 26.09.1-edge

    Launching `workflows/transform.nf` [sick_snyder] revision: 4e35790d32

    executor >  local (6)
    [c8/123b29] SAY…imestamped_Bob-output.txt) | 3 of 3 ✔
    [74/8da573] REV…estamped_Alice-output.txt) | 3 of 3 ✔

    Outputs:

      /workspaces/training/side-quests/composing_pipelines/results

      upper:
        - UPPER-timestamped_Charlie-output.txt
        - UPPER-timestamped_Bob-output.txt
        - UPPER-timestamped_Alice-output.txt

      reversed:
        - REVERSED-UPPER-timestamped_Charlie-output.txt
        - REVERSED-UPPER-timestamped_Bob-output.txt
        - REVERSED-UPPER-timestamped_Alice-output.txt
    ```

It works because the files from section 1.1 are still in `results/`.
Without them, for example in a fresh directory or after someone renames the outputs of `greeting.nf`, the glob matches nothing.
The run then completes without running a single task and publishes nothing, with no error.

Composing the two stages removes that dependency: `TRANSFORM_WORKFLOW` will receive the timestamped greetings from `GREETING_WORKFLOW` as a channel.

### 2.2. Make it composable

Apply the same three changes as in section 1.2: name the workflow, replace the `channel.fromPath(...)` input with `take:`, and replace `publish:`/`output {}` with `emit:`.

The finished file should look like this:

```groovy title="workflows/transform.nf" linenums="1" hl_lines="6 7 8 15 16 17"
#!/usr/bin/env nextflow

include { SAY_HELLO_UPPER } from '../modules/say_hello_upper'
include { REVERSE_TEXT } from '../modules/reverse_text'

workflow TRANSFORM_WORKFLOW {
    take:
    input_ch // Input channel with greetings

    main:
    // Apply transformations in sequence
    upper_ch = SAY_HELLO_UPPER(input_ch)
    reversed_ch = REVERSE_TEXT(upper_ch)

    emit:
    upper = upper_ch // Uppercase greetings
    reversed = reversed_ch // Reversed uppercase greetings
}
```

The transform workflow is now composable and ready to be imported into the main workflow.

### 2.3. Update and test the main workflow

Now update the main workflow to call the transformation workflow.

#### 2.3.1. Include the transformation workflow and call it

Add the `include` statement, a call to `TRANSFORM_WORKFLOW` chained on the timestamped greetings, and the two new `publish:` entries:

=== "After"

    ```groovy title="main.nf" linenums="1" hl_lines="4 13 14 18 19"
    #!/usr/bin/env nextflow

    include { GREETING_WORKFLOW } from './workflows/greeting'
    include { TRANSFORM_WORKFLOW } from './workflows/transform'

    workflow {
        main:
        names = channel.of('Alice', 'Bob', 'Charlie')

        // Run the greeting workflow
        GREETING_WORKFLOW(names)

        // Run the transform workflow
        TRANSFORM_WORKFLOW(GREETING_WORKFLOW.out.timestamped)

        publish:
        greetings = GREETING_WORKFLOW.out.greetings
        upper = TRANSFORM_WORKFLOW.out.upper
        reversed = TRANSFORM_WORKFLOW.out.reversed
    }
    ```

=== "Before"

    ```groovy title="main.nf" linenums="1"
    #!/usr/bin/env nextflow

    include { GREETING_WORKFLOW } from './workflows/greeting'

    workflow {
        main:
        names = channel.of('Alice', 'Bob', 'Charlie')

        // Run the greeting workflow
        GREETING_WORKFLOW(names)

        publish:
        greetings = GREETING_WORKFLOW.out.greetings
    }
    ```

`TRANSFORM_WORKFLOW` takes `GREETING_WORKFLOW.out.timestamped` as its input, so the timestamped greetings go straight from one workflow to the other.

#### 2.3.2. Update the output block

Add `upper` and `reversed` entries to the `output {}` block, each with a `path` directive for its subdirectory:

=== "After"

    ```groovy title="main.nf" linenums="22" hl_lines="5 6 7 8 9 10"
    output {
        greetings {
            path 'greetings'
        }
        upper {
            path 'upper'
        }
        reversed {
            path 'reversed'
        }
    }
    ```

=== "Before"

    ```groovy title="main.nf" linenums="22"
    output {
        greetings {
            path 'greetings'
        }
    }
    ```

Each output gets its own subdirectory of `results/`.

#### 2.3.3. Run the complete pipeline

Run the pipeline to test that it all works:

```bash
nextflow run main.nf
```

??? success "Command output"

    ```console
     N E X T F L O W   ~  version 26.09.1-edge

    Launching `main.nf` [loquacious_stallman] revision: 2a0e54fecd

    executor >  local (15)
    [54/4d7e12] GRE…_NAME (validating Charlie) | 3 of 3 ✔
    [21/d16068] GRE…W:SAY_HELLO (greeting Bob) | 3 of 3 ✔
    [a3/4b46ff] GRE…ing timestamp to greeting) | 3 of 3 ✔
    [0a/9af50d] TRA…estamped_Alice-output.txt) | 3 of 3 ✔
    [73/992c26] TRA…estamped_Alice-output.txt) | 3 of 3 ✔

    Outputs:

      /workspaces/training/side-quests/composing_pipelines/results

      greetings:
        - greetings/Charlie-output.txt
        - greetings/Alice-output.txt
        - greetings/Bob-output.txt

      upper:
        - upper/UPPER-timestamped_Bob-output.txt
        - upper/UPPER-timestamped_Charlie-output.txt
        - upper/UPPER-timestamped_Alice-output.txt

      reversed:
        - reversed/REVERSED-UPPER-timestamped_Bob-output.txt
        - reversed/REVERSED-UPPER-timestamped_Charlie-output.txt
        - reversed/REVERSED-UPPER-timestamped_Alice-output.txt
    ```

??? abstract "New content added under `results/`"

    ```console
    results/
    ├── greetings
    │   ├── Alice-output.txt
    │   ├── Bob-output.txt
    │   └── Charlie-output.txt
    ├── reversed
    │   ├── REVERSED-UPPER-timestamped_Alice-output.txt
    │   ├── REVERSED-UPPER-timestamped_Bob-output.txt
    │   └── REVERSED-UPPER-timestamped_Charlie-output.txt
    └── upper
        ├── UPPER-timestamped_Alice-output.txt
        ├── UPPER-timestamped_Bob-output.txt
        └── UPPER-timestamped_Charlie-output.txt
    ```

    The files at the top level of `results/` come from the standalone runs in sections 1.1 and 2.1.
    This run doesn't read or write them.

??? abstract "File contents"

    ```console title="results/reversed/REVERSED-UPPER-timestamped_Alice-output.txt"
    !ECILA ,OLLEH ]43:21:71 50-01-6202[
    ```

The pipeline is working end-to-end: the greeting has been uppercased and reversed.

`TRANSFORM_WORKFLOW` doesn't read anything from `results/` in this run.
It gets the timestamped greetings from `GREETING_WORKFLOW` as a channel, so it no longer depends on `greeting.nf` having run first, on the output directory, or on the file names.

### Takeaway

`TRANSFORM_WORKFLOW` takes its input as a channel argument, so the two stages are connected in `main.nf` and run as one pipeline.
The entry workflow decides what gets published and where.

---

## Takeaway

In this part, you composed your greetings pipeline from two named workflows:

- **Named workflows**: `take:` declares a workflow's inputs and `emit:` declares its outputs, so it can be included and called by name
- **The entry workflow**: an unnamed `workflow {}` block is where a run starts, so a script that only defines named workflows can't run on its own
- **Wiring workflows together**: the entry workflow passes `GREETING_WORKFLOW.out.timestamped` to `TRANSFORM_WORKFLOW` as a channel
- **No hand-off through files**: the transform stage no longer depends on the greeting stage's output directory or file names

---

## What's next?

The two workflows are connected, but the connection is informal: nothing in the code says what `GREETING_WORKFLOW` emits or what `TRANSFORM_WORKFLOW` expects.
In Part 2, you'll write that contract down with static types and records, so that `nextflow lint` can check it.

[Continue to Part 2 :material-arrow-right:](02_typed_interfaces.md){ .md-button .md-button--primary }
