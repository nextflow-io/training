# Summary

You have completed the Composing Pipelines training.
This page retells the story of the course, recaps the key patterns at each level of composition, and points you to further reading.

---

## The story so far

You started with a greetings pipeline whose two stages were standalone scripts, the second reading the first one's published files.
By the end, the whole greetings pipeline was itself a building block, included by another pipeline and connected to a separately developed report pipeline in one run.
Each part moved one step up the ladder from modules to workflows, from workflows to a pipeline, and from a pipeline to a pipeline of pipelines, and each step was driven by a problem the previous level couldn't solve.
Throughout, each pipeline kept running on its own as well.

### Part 1: Workflows of Workflows

The greeting and transform stages were standalone scripts, chained by hand: `transform.nf` picked up whatever `greeting.nf` had left in `results/`, so it depended on the other script's output directory and file names.
You gave each stage a name, a `take:` block and an `emit:` block, and called both from an entry workflow in `main.nf`.
The transform stage now receives the timestamped greetings as a channel instead, and you saw why a named workflow can't run on its own.

### Part 2: Typed Interfaces

When you typed the greeting side's modules, `GREETING_WORKFLOW` started emitting records, and `TRANSFORM_WORKFLOW` failed at run time inside a process you hadn't touched.
You wrote the contract down on both sides of the connection: a shared `Greeting` record type, typed `emit:` and `take:` blocks, and a typed caller.
After that, changing what one workflow emits produced a `nextflow lint` error at the call where the two workflows meet, before anything ran.
You also saw the limit: two outputs of the same type can be swapped without any error, because types check shape, not meaning.
The result was a fully typed pipeline, which is what Part 3 needs to make it includable.

### Part 3: Pipelines as Workflows

The pipeline's own interface was still informal, with hard-coded inputs and results that others could only use by reading files after the run.
You gave the greetings pipeline a typed `params {}` block and a complete `output {}` block, which play the roles of `take:` and `emit:` for a whole pipeline, and ran it on its own.
Then you included it in a new top-level pipeline and called it like a workflow.

### Part 4: Chaining Pipelines

A separately developed report pipeline summarizes text files, and you wanted it to run on the greetings pipeline's results.
You first connected the two the traditional way, with two runs and a samplesheet built from published paths.
When Diana joined the input, both pipelines reran every task and the report still left her out, because the samplesheet had gone stale: the same path coupling as Part 1, one level up.
Then you composed them, feeding the greetings pipeline's outputs straight into the report pipeline.
The same change, with `-resume`, then reran only Diana's tasks and the summary, across both pipelines, and put her in the report.

### Part 5: Parameters and Configuration

Once pipelines are nested, parameters and configuration have to cross pipeline boundaries.
You passed parameters to the included pipelines from the command line, from config and from the calling code, and found out that an included pipeline's own `nextflow.config` is not loaded, so the top-level config has to bring in the configuration it needs.
You scoped process settings to one included pipeline with its alias, which also keeps a simple-name selector from matching a same-named process in another pipeline.

### Part 6: Making Pipelines Composable

The report pipeline worked perfectly on its own but misbehaved as soon as another pipeline included it.
You traced the problem to a process that read `params` and `projectDir`, both of which belong to the top-level pipeline, fixed it with explicit inputs and `moduleDir`, and finished with a checklist: the practices you applied in this course, plus a few more to consider.

---

## Key patterns

### Within a pipeline

1.  **Named workflows with an interface**: `take:` declares the inputs, `emit:` declares the outputs, and an unnamed entry workflow calls them.

    ```groovy
    workflow GREETING_WORKFLOW {
        take:
        names: Channel<String>

        main:
        validated_ch = VALIDATE_NAME(names)
        greetings_ch = SAY_HELLO(validated_ch)
        timestamped_ch = TIMESTAMP_GREETING(greetings_ch)

        emit:
        greetings: Channel<Greeting> = greetings_ch
        timestamped: Channel<Greeting> = timestamped_ch
    }
    ```

2.  **Shared record types at the boundary**: declare a record type once and include it by name wherever it's used.

    ```groovy
    record Greeting {
        name: String
        file: Path
    }
    ```

    ```groovy
    include { Greeting } from '../types'
    ```

3.  **A typed caller**: capture each workflow call's result in a variable and read its outputs as fields, so the connection between the two workflows is visible to `nextflow lint`.

    ```groovy
    greeting = GREETING_WORKFLOW(names)
    transform = TRANSFORM_WORKFLOW(greeting.timestamped)
    ```

4.  **Check the contract**: run `nextflow lint` on typed code, since type errors found during `nextflow run` are only reported as warnings.
    Linting the entry point follows every `include`.

    ```bash
    nextflow lint main.nf
    ```

### Between pipelines

5.  **A pipeline's interface**: the `params {}` block plays the role of `take:`, and `publish:` with an `output {}` block plays the role of `emit:`.
    A pipeline written this way runs on its own and can also be included by another pipeline.

    ```groovy
    params {
        names: Channel<Person>
    }
    ```

    ```bash
    nextflow -C pipelines/greetings/nextflow.config run ./pipelines/greetings --names data/names.csv
    ```

6.  **Including a pipeline**: an included pipeline is called like a named workflow, with a record of parameters as its input.
    Its outputs are emitted on the call result, not published, so the caller re-publishes what it wants in its own `output {}` block.

    ```groovy
    include { workflow as GREETINGS } from './pipelines/greetings'

    greetings = GREETINGS(record(names: params.names))
    ```

7.  **Chaining pipelines**: passing one pipeline's output channel into another pipeline's parameters connects them in a single dataflow, so they run and resume as one pipeline, with no hand-written samplesheet or published paths in between.
    A change to the input reruns only the tasks it affects, in both pipelines.

    ```groovy
    samples = greetings.reversed.map { greeting -> record(id: greeting.name, file: greeting.file) }
    summary = REPORT(params.report + record(input: samples))
    ```

8.  **Parameters and configuration across boundaries**: import an included pipeline's parameters as a record type, so they can be set from the command line (`--report.metric words`) or nested config, with the included pipeline's defaults applied.
    Values supplied by the caller with `+ record(...)` win, and a user-supplied value for the same parameter is silently dropped.
    The included pipeline's own `nextflow.config` is not loaded, so the caller brings in the configuration it needs.

    ```groovy
    include { params as ReportParams ; workflow as REPORT } from './pipelines/report'

    params {
        report: ReportParams
    }
    ```

    ```groovy
    includeConfig 'pipelines/report/conf/modules.config'

    process {
        withName: 'REPORT:.*' {
            memory = 1.GB
        }
    }
    ```

    When a simple process name in the included file could match processes in more than one included pipeline, don't include that file; set its values in the top-level config with the alias instead, such as `withName: 'REPORT:SUMMARIZE'`.

9.  **Composable habits**: a pipeline that passes everything its processes need as declared inputs, and locates its own files with `moduleDir` rather than `projectDir`, keeps working when another pipeline includes it.
    The checklist at the end of Part 6 collects these habits.

    ```groovy
    banner = file("${moduleDir}/assets/banner.txt")
    summary = SUMMARIZE(counts.collect(), params.title, banner)
    ```

---

## Additional resources

- [Workflows](https://www.nextflow.io/docs/latest/workflow.html): named workflows, `take:`, `emit:` and entry workflows
- [Typed workflows](https://www.nextflow.io/docs/latest/workflow-typed.html): typed `take:` and `emit:`, and pipeline composition
- [Typed processes](https://www.nextflow.io/docs/latest/process-typed.html): typed inputs and outputs, including records
- [Static typing](https://www.nextflow.io/docs/latest/static-typing.html): records and the `Channel` and `Value` types
- [Migrating to static typing](https://www.nextflow.io/docs/latest/tutorials/static-types.html): converting existing code step by step
- [Workflow outputs](https://www.nextflow.io/docs/latest/workflow-outputs.html): `publish:` and the `output {}` block

---

## What's next?

Continue to the [Next Steps](next_steps.md) for ideas on where to take what you've learned.

[Continue to Next Steps :material-arrow-right:](next_steps.md){ .md-button .md-button--primary }
