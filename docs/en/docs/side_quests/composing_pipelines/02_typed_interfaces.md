<!-- TODO(26.10): recapture console output on 26.10 release -->

# Part 2: Typed Interfaces

At the end of Part 1, your greetings pipeline was built from two named workflows, and one line in `main.nf` connected them:

```groovy title="main.nf" linenums="14"
    TRANSFORM_WORKFLOW(GREETING_WORKFLOW.out.timestamped)
```

That line is a contract between the two workflows: whatever `GREETING_WORKFLOW` emits as `timestamped` has to be what `TRANSFORM_WORKFLOW` expects.
Neither workflow writes it down.
`TRANSFORM_WORKFLOW` describes its input with a comment, `// Input channel with greetings`, and the only way to find out what `GREETING_WORKFLOW` emits is to read its processes.
If one side changes, nothing tells the other side.

The [Types and Records](../types_and_records/index.md) side quest showed how types describe a workflow's interface to the entry workflow that calls it.
This part applies the same tools to the connection between two workflows in one pipeline.
You'll change one side of the connection and watch the other side break, then write the contract down on both sides so that `nextflow lint` checks the connection itself.

There's a second reason to type the whole pipeline: in Part 3, another pipeline will include yours, and only a fully typed pipeline can be included.

!!! tip "Starting from here?"

    If you're joining at this part, copy the solution from Part 1 into your working directory to use as your starting point:

    ```bash
    cd side-quests/composing_pipelines
    cp -r ../solutions/composing_pipelines/1/* .
    ```

### Learning goals

By the end of this part, you'll be able to:

- Explain what the connection between two workflows depends on, and what happens when one side changes
- Convert a pipeline's modules to typed processes around a shared record type
- Declare typed `emit:` and `take:` blocks on both sides of a connection between workflows
- Use `nextflow lint` to find a broken connection at the call that makes it
- Recognize what type checking doesn't catch

### Prerequisites

This part builds on the [Types and Records](../types_and_records/index.md) side quest.
It uses typed processes, `record` types, `Channel<T>` and workflow call results without re-explaining them, so take that side quest first if any of these are new to you.

---

## 1. Change one side of the connection

Typing a pipeline starts with its processes, and a natural way to work is one workflow at a time.
Start with the modules behind `GREETING_WORKFLOW`, and run the pipeline once they're done.

### 1.1. Define a shared `Greeting` record type

Every step after `SAY_HELLO` passes around a greeting file, and each file belongs to a name.
At the moment that name only survives in the file name: `TIMESTAMP_GREETING` recovers it with `greeting_file.baseName`, and every later step relies on the naming convention.

A record carries the name and the file together, as named fields.
Both workflows will use it, so declare it once, in its own file.
Create `types.nf` in the project directory:

```groovy title="types.nf" linenums="1"
#!/usr/bin/env nextflow

nextflow.enable.types = true

record Greeting {
    name: String
    file: Path
}
```

Any script that needs to talk about greetings includes `Greeting` from `types.nf` by name.

### 1.2. Type `TIMESTAMP_GREETING`

`TIMESTAMP_GREETING` takes a greeting and returns a timestamped one, so both its input and its output become a `Greeting`.
Replace the contents of `modules/timestamp_greeting.nf`:

=== "After"

    ```groovy title="modules/timestamp_greeting.nf" linenums="1" hl_lines="1 3 9 12 16"
    nextflow.enable.types = true

    include { Greeting } from '../types'

    process TIMESTAMP_GREETING {
        tag "adding timestamp to greeting"

        input:
        greeting: Greeting

        output:
        record(name: greeting.name, file: file("timestamped_${greeting.file.baseName}.txt"))

        script:
        """
        echo "[\$(date '+%Y-%m-%d %H:%M:%S')] \$(cat ${greeting.file})" > timestamped_${greeting.file.baseName}.txt
        """
    }
    ```

=== "Before"

    ```groovy title="modules/timestamp_greeting.nf" linenums="1" hl_lines="5 8 11 13"
    process TIMESTAMP_GREETING {
        tag "adding timestamp to greeting"

        input:
        path greeting_file

        output:
        path 'timestamped_*.txt'

        script:
        def base_name = greeting_file.baseName
        """
        echo "[\$(date '+%Y-%m-%d %H:%M:%S')] \$(cat ${greeting_file})" > timestamped_${base_name}.txt
        """
    }
    ```

The `include` makes the `Greeting` type available in this script.
Because `file` is declared as a `Path`, Nextflow stages it into the task directory like a `path` input, and the script reads it as `greeting.file`.
The output passes `greeting.name` through unchanged, so later steps don't have to reconstruct it from a file name.

### 1.3. Type the other greeting modules

Apply the same treatment to the other two modules that `GREETING_WORKFLOW` uses, leaving their `script:` blocks as they are:

- `modules/validate_name.nf`: enable static typing, take `name: String`, and return `name`.
- `modules/say_hello.nf`: enable static typing, take `name: String`, and return `#!groovy record(name: name, file: file("${name}-output.txt"))`.

`SAY_HELLO` returns a record with exactly the fields of `Greeting`, so it can be used wherever a `Greeting` is expected without including the type.
The finished files are in the solution directory, if you'd rather copy them:

```bash
cp ../solutions/composing_pipelines/2/modules/{validate_name,say_hello}.nf modules/
```

### 1.4. Run the pipeline

The workflows and `main.nf` haven't changed, and typed modules can be called from untyped scripts, so run the pipeline to check that the greeting side still works:

```bash
nextflow run main.nf
```

??? failure "Command output"

    ```console
     N E X T F L O W   ~  version 26.09.1-edge

    Launching `main.nf` [dreamy_watson] revision: 2a0e54fecd

    executor >  local (9)
    [9c/e184da] GRE…_NAME (validating Charlie) | 3 of 3 ✔
    [e8/11cb47] GRE…SAY_HELLO (greeting Alice) | 3 of 3 ✔
    [a6/a39cfa] GRE…ing timestamp to greeting) | 3 of 3 ✔
    [-        ] TRA…M_WORKFLOW:SAY_HELLO_UPPER -
    [-        ] TRA…FORM_WORKFLOW:REVERSE_TEXT -
    ERROR ~ Error executing process > 'TRANSFORM_WORKFLOW:SAY_HELLO_UPPER (3)'

    Caused by:
      Not a valid path value type: nextflow.util.RecordMap ([name:Charlie, file:/workspaces/training/side-quests/composing_pipelines/work/a6/a39cfa47ca144245fff4bb91f06925/timestamped_Charlie-output.txt])

    Tip: when you have fixed the problem you can continue the execution adding the option `-resume` to the run command line

     -- Check '.nextflow.log' file for details
    ```

All three greeting processes succeed.
The failure comes from `TRANSFORM_WORKFLOW`, which you didn't touch: its first process, `SAY_HELLO_UPPER`, received a record where it expected a file.

Changing `TIMESTAMP_GREETING` changed what `GREETING_WORKFLOW` emits as `timestamped`, from files to records, and `TRANSFORM_WORKFLOW` still expects files.
Nothing at the call in `main.nf` records either side of that, so the run started anyway, and the error is reported inside the other workflow, in terms of its internals.
In a larger pipeline, the producing workflow might be edited months later by someone who never reads the consuming one, and the failure might come hours into a run.

### Takeaway

The connection between two workflows depends on what one emits and what the other expects, and an untyped interface records neither.
Changing what `GREETING_WORKFLOW` emits broke `TRANSFORM_WORKFLOW` at run time, with an error that points at a process rather than at the connection.

### What's next?

Converting the transformation modules would fix this run, but the connection would still be unrecorded.
The next section converts them and types both sides of the connection.

---

## 2. Type both sides of the connection

The connection has two sides: the `emit:` block of `GREETING_WORKFLOW` and the `take:` block of `TRANSFORM_WORKFLOW`.
Both need typed modules behind them, then a declared type, and the entry workflow that wires them together becomes typed too.

### 2.1. Type the transformation modules

`SAY_HELLO_UPPER` and `REVERSE_TEXT` follow the same pattern as `TIMESTAMP_GREETING`: include `Greeting`, take `greeting: Greeting`, and return a record with the same `name` and a new file.
In each script, use `greeting.file` for the input file and `greeting.file.name` where the file name is part of a string.
The output file names stay the same.
To copy the finished files instead:

```bash
cp ../solutions/composing_pipelines/2/modules/{say_hello_upper,reverse_text}.nf modules/
```

### 2.2. Type `GREETING_WORKFLOW`

Open `workflows/greeting.nf`.
Enable static typing, include the `Greeting` type, and add types to the `take:` and `emit:` declarations:

=== "After"

    ```groovy title="workflows/greeting.nf" linenums="1" hl_lines="3 8 12 16 21 22"
    #!/usr/bin/env nextflow

    nextflow.enable.types = true

    include { VALIDATE_NAME } from '../modules/validate_name'
    include { SAY_HELLO } from '../modules/say_hello'
    include { TIMESTAMP_GREETING } from '../modules/timestamp_greeting'
    include { Greeting } from '../types'

    workflow GREETING_WORKFLOW {
        take:
        names: Channel<String>

        main:
        // Chain processes: validate -> create greeting -> add timestamp
        validated_ch = VALIDATE_NAME(names)
        greetings_ch = SAY_HELLO(validated_ch)
        timestamped_ch = TIMESTAMP_GREETING(greetings_ch)

        emit:
        greetings: Channel<Greeting> = greetings_ch
        timestamped: Channel<Greeting> = timestamped_ch
    }
    ```

=== "Before"

    ```groovy title="workflows/greeting.nf" linenums="1" hl_lines="9 13 18 19"
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

The `emit:` block now states what this side of the connection promises: two channels of `Greeting` records.
The comments are gone because the types say the same thing, and Nextflow can check types.
The input is renamed from `names_ch` to `names`, since the type already says it's a channel.

### 2.3. Type `TRANSFORM_WORKFLOW`

Apply the same changes to `workflows/transform.nf`: enable static typing, include `Greeting`, declare the input as `greetings: Channel<Greeting>`, and type both emits as `Channel<Greeting>`.

??? abstract "File contents"

    ```groovy title="workflows/transform.nf" linenums="1" hl_lines="3 7 11 15 19 20"
    #!/usr/bin/env nextflow

    nextflow.enable.types = true

    include { SAY_HELLO_UPPER } from '../modules/say_hello_upper'
    include { REVERSE_TEXT } from '../modules/reverse_text'
    include { Greeting } from '../types'

    workflow TRANSFORM_WORKFLOW {
        take:
        greetings: Channel<Greeting>

        main:
        // Apply transformations in sequence
        upper_ch = SAY_HELLO_UPPER(greetings)
        reversed_ch = REVERSE_TEXT(upper_ch)

        emit:
        upper: Channel<Greeting> = upper_ch
        reversed: Channel<Greeting> = reversed_ch
    }
    ```

The `take:` block is the other side of the contract: `TRANSFORM_WORKFLOW` expects a channel of `Greeting` records, which is exactly what `GREETING_WORKFLOW` declares for `timestamped`.

### 2.4. Update the caller

A typed script reads a workflow's outputs from the result of the call rather than through `.out`, as in the [Types and Records](../types_and_records/index.md) side quest.
Open `main.nf`, enable static typing, assign each workflow call to a variable, and read the outputs from those variables.
The timestamped greetings are now useful records in their own right, so publish them too:

=== "After"

    ```groovy title="main.nf" linenums="1" hl_lines="3 13 16 19 20 21 22"
    #!/usr/bin/env nextflow

    nextflow.enable.types = true

    include { GREETING_WORKFLOW } from './workflows/greeting'
    include { TRANSFORM_WORKFLOW } from './workflows/transform'

    workflow {
        main:
        names = channel.of('Alice', 'Bob', 'Charlie')

        // Run the greeting workflow
        greeting = GREETING_WORKFLOW(names)

        // Run the transform workflow
        transform = TRANSFORM_WORKFLOW(greeting.timestamped)

        publish:
        greetings = greeting.greetings
        timestamped = greeting.timestamped
        upper = transform.upper
        reversed = transform.reversed
    }
    ```

=== "Before"

    ```groovy title="main.nf" linenums="1" hl_lines="11 14 17 18 19"
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

`greeting.timestamped` carries the type declared in the `emit:` block of `GREETING_WORKFLOW`, so line 16 is where the two declared types meet.

Then add a `timestamped` entry to the `output {}` block:

=== "After"

    ```groovy title="main.nf" linenums="25" hl_lines="5 6 7"
    output {
        greetings {
            path 'greetings'
        }
        timestamped {
            path 'timestamped'
        }
    ```

=== "Before"

    ```groovy title="main.nf" linenums="25"
    output {
        greetings {
            path 'greetings'
        }
    ```

### 2.5. Run the typed pipeline

Run the pipeline:

```bash
nextflow run main.nf
```

??? success "Command output"

    ```console
     N E X T F L O W   ~  version 26.09.1-edge

    Launching `main.nf` [pensive_bose] revision: 936591edc4

    executor >  local (15)
    [0e/8caee4] GRE…TE_NAME (validating Alice) | 3 of 3 ✔
    [39/76fc5c] GRE…W:SAY_HELLO (greeting Bob) | 3 of 3 ✔
    [92/8e2aa9] GRE…ing timestamp to greeting) | 3 of 3 ✔
    [3a/8b3930] TRA…imestamped_Bob-output.txt) | 3 of 3 ✔
    [64/32f6ed] TRA…imestamped_Bob-output.txt) | 3 of 3 ✔

    Outputs:

      /workspaces/training/side-quests/composing_pipelines/results

      greetings:
        - {name: Bob, file: greetings/Bob-output.txt}
        - {name: Alice, file: greetings/Alice-output.txt}
        - {name: Charlie, file: greetings/Charlie-output.txt}

      timestamped:
        - {name: Charlie, file: timestamped/timestamped_Charlie-output.txt}
        - {name: Alice, file: timestamped/timestamped_Alice-output.txt}
        - {name: Bob, file: timestamped/timestamped_Bob-output.txt}

      upper:
        - {name: Charlie, file: upper/UPPER-timestamped_Charlie-output.txt}
        - {name: Alice, file: upper/UPPER-timestamped_Alice-output.txt}
        - {name: Bob, file: upper/UPPER-timestamped_Bob-output.txt}

      reversed:
        - {name: Charlie, file: reversed/REVERSED-UPPER-timestamped_Charlie-output.txt}
        - {name: Alice, file: reversed/REVERSED-UPPER-timestamped_Alice-output.txt}
        - {name: Bob, file: reversed/REVERSED-UPPER-timestamped_Bob-output.txt}
    ```

The published outputs are now records rather than bare files.
The `output {}` block only needs a `path` for each output: Nextflow publishes every file it finds in each record (here, the `file` field) to that directory.
The results are the same as at the end of Part 1, plus the timestamped greetings, and every script in the pipeline is now typed.

### Takeaway

Both sides of the connection now declare their types: `GREETING_WORKFLOW` emits `timestamped` as a `Channel<Greeting>`, and `TRANSFORM_WORKFLOW` takes a `Channel<Greeting>`.
The entry workflow wires them together through typed call results instead of `.out`.

### What's next?

Time to find out what the declared types change when one side of the connection changes again.

---

## 3. Let the contract catch the change

In section 1, a change to one workflow broke the other at run time.
This section repeats that kind of change with the connection typed, then looks at a mistake the types can't see.

### 3.1. Change what `GREETING_WORKFLOW` emits

Suppose someone decides that downstream steps only need the timestamped files, and changes `GREETING_WORKFLOW` to emit those instead of the records:

=== "After"

    ```groovy title="workflows/greeting.nf" linenums="20" hl_lines="3"
        emit:
        greetings: Channel<Greeting> = greetings_ch
        timestamped: Channel<Path> = timestamped_ch.map { g -> g.file }
    ```

=== "Before"

    ```groovy title="workflows/greeting.nf" linenums="20" hl_lines="3"
        emit:
        greetings: Channel<Greeting> = greetings_ch
        timestamped: Channel<Greeting> = timestamped_ch
    ```

`greeting.nf` is consistent on its own, and nothing in `main.nf` or `transform.nf` changed.
Check the pipeline from its entry point:

```bash
nextflow lint main.nf
```

??? failure "Command output"

    ```console
    Linting Nextflow code..
    Linting: main.nf
    Error main.nf:16:36: Argument with type Channel<Path> is not compatible with parameter of type Channel<Greeting>
    │  16 |     transform = TRANSFORM_WORKFLOW(greeting.timestamped)
    ╰     |                                    ^^^^^^^^^^^^^^^^^^^^

    Nextflow linting complete!
     ❌ 1 file had 1 error
     ✅ 8 files had no errors
    ```

The change was made in `greeting.nf`, but the error is reported in `main.nf`, at the call where the two workflows meet, before anything runs.
It states the broken contract: `GREETING_WORKFLOW` now emits a `Channel<Path>`, and `TRANSFORM_WORKFLOW` declares a `Channel<Greeting>`.
Compare that with section 1.4, where the same kind of change surfaced as a failed task inside `TRANSFORM_WORKFLOW`.

`nextflow lint` is the check to rely on here.
`nextflow run` performs the same type check but only reports a warning and carries on, as covered in the [Types and Records](../types_and_records/index.md) side quest.

Put the line back to `timestamped: Channel<Greeting> = timestamped_ch`.

### 3.2. Swap two outputs of the same type

Types check the shape of what crosses a connection, not its meaning.
In `main.nf`, pass the original greetings to `TRANSFORM_WORKFLOW` instead of the timestamped ones:

=== "After"

    ```groovy title="main.nf" linenums="15" hl_lines="2"
        // Run the transform workflow
        transform = TRANSFORM_WORKFLOW(greeting.greetings)
    ```

=== "Before"

    ```groovy title="main.nf" linenums="15" hl_lines="2"
        // Run the transform workflow
        transform = TRANSFORM_WORKFLOW(greeting.timestamped)
    ```

```bash
nextflow lint main.nf
```

??? success "Command output"

    ```console
    Linting Nextflow code..
    Linting: main.nf
    Nextflow linting complete!
     ✅ 9 files had no errors
    ```

`greeting.greetings` and `greeting.timestamped` are both `Channel<Greeting>`, so the call is valid as far as the types go.
The pipeline would run without complaint and publish uppercased and reversed greetings without timestamps, as `UPPER-Alice-output.txt` and so on.
Choosing the right output is still up to whoever writes the call, and clear output names are what help them.

Put the call back to `transform = TRANSFORM_WORKFLOW(greeting.timestamped)` and lint the pipeline once more:

```bash
nextflow lint main.nf
```

??? success "Command output"

    ```console
    Linting Nextflow code..
    Linting: main.nf
    Nextflow linting complete!
     ✅ 9 files had no errors
    ```

Only `main.nf` is listed, but `nextflow lint` follows each `include`, so the nine files are `main.nf`, `types.nf`, the five modules and both workflows.

### Takeaway

With both sides of the connection typed, changing what one workflow emits produces a `nextflow lint` error at the call where the workflows meet, before anything runs.
Type checking can't tell apart two outputs of the same type, so the wiring still needs a human eye.

---

## Takeaway

In this part, you turned the connection between your two workflows into an explicit, typed contract:

- **The connection depends on both sides**: changing what one workflow emits can break another workflow you didn't touch
- **A shared record type**: `Greeting`, declared once in `types.nf`, carries a greeting's name and file through every step
- **Typed connections**: `emit:` on one side and `take:` on the other declare the same `Channel<Greeting>`
- **Errors at the call**: `nextflow lint` reports a mismatched connection at the line that wires the workflows together
- **The limit**: types check shape, not meaning, so two outputs of the same type can still be swapped

---

## What's next?

Your greetings pipeline is now fully typed, from the processes up to the entry workflow.

Every interface you've typed so far is internal to the pipeline.
The pipeline as a whole has an interface too, but it's still informal: the names are hard-coded in `main.nf`, and other code can only use the results by reading files from `results/`.

In Part 3, you'll give the pipeline a typed `params {}` block for its inputs and use its `output {}` block as its outputs, so that another pipeline can include your whole pipeline and call it like a workflow.
That step relies on the work in this part: only a fully typed pipeline can be included by another one.

[Continue to Part 3 :material-arrow-right:](03_pipelines_as_workflows.md){ .md-button .md-button--primary }
