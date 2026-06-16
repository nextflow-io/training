# Data Lineage

Scientific workflows produce chains of computations — each output depends on specific inputs, parameters, and tool versions.
Nextflow's data lineage feature captures this provenance automatically, giving you a queryable record of every workflow run, task execution, and output file.

!!! warning "Experimental feature"

    Data lineage is an experimental feature introduced in Nextflow 25.04.
    The feature is under active development and its behavior may change in future releases.
    Consult the [official documentation](https://docs.seqera.io/nextflow/tutorials/data-lineage) for the latest details.

### Learning goals

By the end of this side quest, you'll be able to:

- Enable lineage tracking in a Nextflow workflow
- Explore lineage records using the CLI (`list`, `view`, `find`, `render`, `diff`)
- Understand the three lineage record types: `WorkflowRun`, `TaskRun`, and `FileOutput`
- Understand why the `output {}` block is required for file-level provenance tracking
- Compare two workflow runs to identify exactly what changed between them
- Use `channel.fromLineage()` to create a channel from lineage records and feed tracked outputs into downstream workflows

### Prerequisites

Before taking on this side quest, you should:

- Have completed the [Hello Nextflow](../../hello_nextflow/index.md) tutorial or equivalent beginner's course.
- Be comfortable with basic Nextflow concepts: processes, channels, and configuration.

---

## 0. Get started

#### Open the training codespace

If you haven't yet done so, make sure to open the training environment as described in the [Environment Setup](../../envsetup/index.md).

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Move into the project directory

Move into the directory where the files for this tutorial are located.

```bash
cd side-quests/lineage
```

You can set VSCode to focus on this directory:

```bash
code .
```

#### Review the materials

You'll find a workflow file, a configuration file, and a `data` directory with a CSV input file.

??? abstract "Directory contents"

    ```console
    .
    ├── data
    │   └── greetings.csv
    ├── main.nf
    └── nextflow.config
    ```

The workflow in `main.nf` has two processes.
`SAYHELLO` writes each greeting to a text file, and `CONVERTTOUPPER` converts that file to uppercase.
This is the same pipeline logic used in [Hello Nextflow](../../hello_nextflow/index.md).

The `greetings.csv` file contains four greetings:

```console title="data/greetings.csv"
Hello
Bonjour
Hola
Ciao
```

#### Review the assignment

Your challenge is to add lineage tracking to the workflow, explore the records it produces, and use the lineage CLI to query and compare runs.

For each section:

1. **Apply the code changes** to the workflow or configuration file as shown
2. **Run the workflow** and observe the output
3. **Use the lineage CLI** to inspect the records produced
4. **Complete the exercises** before reading ahead

#### Readiness checklist

Think you're ready to dive in?

- [ ] I understand the goal of this course and its prerequisites
- [ ] My codespace is up and running
- [ ] I've set my working directory appropriately
- [ ] I understand the assignment

If you can check all the boxes, you're good to go.

---

## 1. Enable lineage tracking

Open `main.nf` to examine the starting workflow:

```groovy title="main.nf" linenums="1"
#!/usr/bin/env nextflow

process SAYHELLO {
    input:
    val greeting

    output:
    path "${greeting}-output.txt"

    script:
    """
    echo '${greeting}' > '${greeting}-output.txt'
    """
}

process CONVERTTOUPPER {
    input:
    path input_file

    output:
    path "UPPER-${input_file}"

    script:
    """
    cat '${input_file}' | tr '[a-z]' '[A-Z]' > 'UPPER-${input_file}'
    """
}

params {
    input: Path = 'data/greetings.csv'
}

workflow {
    main:
    greeting_ch = channel.fromPath(params.input)
        .splitCsv()
        .map { line -> line[0] }

    SAYHELLO(greeting_ch)
    CONVERTTOUPPER(SAYHELLO.out)

    publish:
    uppercased = CONVERTTOUPPER.out
}

output {
    uppercased {
        path 'results'
        mode 'copy'
    }
}
```

The workflow reads greetings from a CSV file, writes each one to a text file with `SAYHELLO`, converts it to uppercase with `CONVERTTOUPPER`, and publishes the results to a `results/` directory.
This is the same pipeline used in [Hello Nextflow](../../hello_nextflow/index.md), extended with an `output {}` block to publish results.

Run it once to verify everything works before adding lineage tracking:

```bash
nextflow run main.nf
```

??? success "Command output"

    ```console
     N E X T F L O W   ~  version 25.10.4

    Launching `main.nf` [sharp_bohr] DSL2 - revision: abc123def4

    executor >  local (8)
    [1a/2b3c4d] SAYHELLO (4)       [100%] 4 of 4 ✔
    [5e/6f7a8b] CONVERTTOUPPER (4) [100%] 4 of 4 ✔
    ```

### 1.1. Add lineage configuration

Add `lineage.enabled = true` to `nextflow.config`:

=== "After"

    ```groovy title="nextflow.config" linenums="1" hl_lines="3"
    docker.enabled = true

    lineage.enabled = true
    ```

=== "Before"

    ```groovy title="nextflow.config" linenums="1"
    docker.enabled = true
    ```

!!! tip "Changing the lineage store location"

    By default, lineage data is stored in a `.lineage` directory within your working directory.
    To store it elsewhere, add:

    ```groovy
    lineage.store.location = '/path/to/your/lineage/store'
    ```

Run the workflow with lineage tracking enabled:

```bash
nextflow run main.nf
```

??? success "Command output"

    ```console
     N E X T F L O W   ~  version 25.10.4

    Launching `main.nf` [happy_dijkstra] DSL2 - revision: abc123def4

    executor >  local (8)
    [1a/2b3c4d] SAYHELLO (4)       [100%] 4 of 4 ✔
    [5e/6f7a8b] CONVERTTOUPPER (4) [100%] 4 of 4 ✔
    ```

The workflow output looks the same as before, but Nextflow has now written provenance records to the `.lineage` directory.

??? abstract "Directory contents — `.lineage/`"

    ```console
    .lineage/
    ├── .history/
    │   └── a1b2c3d4e5f6789012345678          ← one entry per workflow run
    ├── a1b2c3d4e5f6789012345678/
    │   └── .data.json                         ← WorkflowRun
    ├── b2c3d4e5f6789012345678a1/
    │   └── .data.json                         ← TaskRun (SAYHELLO/Hello)
    ├── b3c4d5e6f7890123456789a2/
    │   └── .data.json                         ← TaskRun (SAYHELLO/Bonjour)
    ├── b4c5d6e7f8901234567890a3/
    │   └── .data.json                         ← TaskRun (SAYHELLO/Hola)
    ├── b5c6d7e8f9012345678901a4/
    │   └── .data.json                         ← TaskRun (SAYHELLO/Ciao)
    ├── c3d4e5f6789012345678a1b2/
    │   └── .data.json                         ← TaskRun (CONVERTTOUPPER/Hello)
    ├── c4d5e6f7890123456789a2b3/
    │   └── .data.json                         ← TaskRun (CONVERTTOUPPER/Bonjour)
    ├── c5d6e7f8901234567890a3b4/
    │   └── .data.json                         ← TaskRun (CONVERTTOUPPER/Hola)
    └── c6d7e8f9012345678901a4b5/
        └── .data.json                         ← TaskRun (CONVERTTOUPPER/Ciao)
    ```

    Each record is a `.data.json` file stored in a directory named after its LID hash.
    The `.history/` directory contains one file per workflow run, linking run names to LIDs.

??? abstract "Contents of a record file — `WorkflowRun`"

    ```json title=".lineage/a1b2c3d4e5f6789012345678/.data.json"
    {
      "version": "lineage/v1beta1",
      "kind": "WorkflowRun",
      "spec": {
        "workflow": {
          "scriptFiles": [
            {
              "path": "/workspace/side-quests/lineage/main.nf",
              "checksum": {
                "value": "78910abc",
                "algorithm": "nextflow",
                "mode": "standard"
              }
            }
          ]
        },
        "sessionId": "4f02559e-9ebd-41d8-8ee2-a8d1e4f09c67",
        "name": "happy_dijkstra",
        "params": [
          {
            "type": "Path",
            "name": "input",
            "value": "data/greetings.csv"
          }
        ],
        "config": {
          "docker": { "enabled": true },
          "lineage": { "enabled": true }
        },
        "metadata": {
          "start": "2026-06-15T10:00:00Z",
          "end": "2026-06-15T10:00:45Z",
          "nextflow": {
            "version": "25.10.4",
            "enable": { "dsl": 2 }
          }
        }
      }
    }
    ```

    Every record follows the same envelope: `version`, `kind`, and `spec`.
    The `kind` field identifies the record type (`WorkflowRun`, `TaskRun`, `FileOutput`, etc.).

### 1.2. Lineage record types

Nextflow records provenance using three primary record types:

| Type          | Stores                                                             |
| ------------- | ------------------------------------------------------------------ |
| `WorkflowRun` | Overall pipeline execution — run name, parameters, start/end times |
| `TaskRun`     | Individual task execution — inputs, script, container, timing      |
| `FileOutput`  | A published output file — path, checksum, size, producing task     |

Each record is identified by a unique **Lineage ID (LID)** formatted as `lid://[hash]`.
The LID is the entry point for all lineage queries.

### Takeaway

In this section, you've learned:

- **How to enable lineage tracking:** Set `lineage.enabled = true` in `nextflow.config`.
- **Record types:** `WorkflowRun` (the overall pipeline execution), `TaskRun` (each individual task), and `FileOutput` (each published output file).
- **Where records are stored:** In `.lineage` within your working directory by default and how to modify the configuration.

---

## 2. Explore lineage records with the CLI

Nextflow provides a `lineage` subcommand for querying the lineage store directly from the terminal.
It gives you five operations: `list` to enumerate runs, `view` to inspect a record by LID, `find` to search across records, `render` to generate a graph visualization, and `diff` to compare two records.

### 2.1. List workflow runs

To see all workflow runs in the lineage store, run:

```bash
nextflow lineage list
```

??? success "Command output"

    ```console
    RUN NAME         SESSION ID                             LID                              TIMESTAMP
    happy_dijkstra   4f02559e-9ebd-41d8-8ee2-a8d1e4f09c67  lid://a1b2c3d4e5f6789012345678   2026-06-15 10:00:00
    ```

Copy the LID for your run: you'll need it in the following steps.

### 2.2. View a WorkflowRun record

Pass the LID to `nextflow lineage view` to inspect the full record for the workflow run:

```bash
nextflow lineage view lid://a1b2c3d4e5f6789012345678
```

??? success "Command output"

    ```json
    {
      "type": "WorkflowRun",
      "lid": "lid://a1b2c3d4e5f6789012345678",
      "name": "happy_dijkstra",
      "params": {
        "input": "data/greetings.csv"
      },
      "start": "2026-06-15T10:00:00Z",
      "end": "2026-06-15T10:00:45Z"
    }
    ```

This record captures the run-level context: which parameters were used and when the pipeline ran.

### 2.3. Find records by field

The `find` subcommand searches the lineage store by record type or field value, without needing to know a LID upfront.
It is useful when you want to locate all records of a given type — for example, all `TaskRun` records across every workflow run, or all `FileOutput` records matching a particular path.
Use it to list all `TaskRun` records:

```bash
nextflow lineage find --type TaskRun
```

??? success "Command output"

    ```console
    lid://b2c3d4e5f6789012345678a1   SAYHELLO        2026-06-15 10:00:10
    lid://b3c4d5e6f7890123456789a2   SAYHELLO        2026-06-15 10:00:11
    lid://b4c5d6e7f8901234567890a3   SAYHELLO        2026-06-15 10:00:12
    lid://b5c6d7e8f9012345678901a4   SAYHELLO        2026-06-15 10:00:13
    lid://c3d4e5f6789012345678a1b2   CONVERTTOUPPER  2026-06-15 10:00:20
    lid://c4d5e6f7890123456789a2b3   CONVERTTOUPPER  2026-06-15 10:00:22
    lid://c5d6e7f8901234567890a3b4   CONVERTTOUPPER  2026-06-15 10:00:24
    lid://c6d7e8f9012345678901a4b5   CONVERTTOUPPER  2026-06-15 10:00:26
    ```

Copy the LID of any `SAYHELLO` task — you'll use it in the next step.

### 2.4. View a TaskRun record

Pass a TaskRun LID to `view` to inspect the full record:

```bash
nextflow lineage view lid://b2c3d4e5f6789012345678a1
```

??? success "Command output"

    ```json
    {
      "type": "TaskRun",
      "lid": "lid://b2c3d4e5f6789012345678a1",
      "name": "SAYHELLO",
      "inputs": [
        {
          "name": "greeting",
          "value": "Hello"
        }
      ],
      "workflowRun": "lid://a1b2c3d4e5f6789012345678",
      "start": "2026-06-15T10:00:10Z",
      "end": "2026-06-15T10:00:15Z"
    }
    ```

This record shows exactly which input value was used for this specific task execution and which workflow run it belongs to.

### 2.5. Render the lineage graph

The `render` subcommand generates an HTML visualization of the lineage for a given run:

```bash
nextflow lineage render lid://a1b2c3d4e5f6789012345678 lineage.html
```

Open the generated HTML file in your browser.
It shows a directed graph connecting `WorkflowRun` → `TaskRun` → `FileOutput` for each published file.
Section 3 explores what this graph looks like when the `output {}` block is absent.

### Takeaway

In this section, you've learned how to use:

- **`nextflow lineage list`:** View all workflow runs in the store with their LIDs.
- **`nextflow lineage find`:** Query records by type or field without knowing the LID in advance.
- **`nextflow lineage view [LID]`:** Inspect the full JSON record for any `WorkflowRun` or `TaskRun`.
- **`nextflow lineage render [LID]`:** Generate a graph visualization of a run's provenance.

---

## 3. File outputs and lineage

The `output {}` block is enables the creation `FileOutput` records in the lineage store.
To understand why it matters, try removing it temporarily and observe what changes.

### 3.1. Remove the output block

Remove the `publish:` section and `output {}` block from `main.nf`:

=== "After"

    ```groovy title="main.nf" linenums="33"
    workflow {
        main:
        greeting_ch = channel.fromPath(params.input)
            .splitCsv()
            .map { line -> line[0] }

        SAYHELLO(greeting_ch)
        CONVERTTOUPPER(SAYHELLO.out)
    }
    ```

=== "Before"

    ```groovy title="main.nf" linenums="33" hl_lines="10 11 14 15 16 17 18 19"
    workflow {
        main:
        greeting_ch = channel.fromPath(params.input)
            .splitCsv()
            .map { line -> line[0] }

        SAYHELLO(greeting_ch)
        CONVERTTOUPPER(SAYHELLO.out)

        publish:
        uppercased = CONVERTTOUPPER.out
    }

    output {
        uppercased {
            path 'results'
            mode 'copy'
        }
    }
    ```

Run the workflow:

```bash
nextflow run main.nf
```

??? success "Command output"

    ```console
     N E X T F L O W   ~  version 25.10.4

    Launching `main.nf` [tender_turing] DSL2 - revision: abc123def4

    executor >  local (8)
    [1a/2b3c4d] SAYHELLO (4)       [100%] 4 of 4, cached: 4 ✔
    [5e/6f7a8b] CONVERTTOUPPER (4) [100%] 4 of 4, cached: 4 ✔
    ```

Now try to view the published outputs of this run using the `#output` suffix:

```bash
nextflow lineage view lid://a1b2c3d4e5f6789012345678#output
```

??? failure "Command output"

    ```console
    No output records found for lid://a1b2c3d4e5f6789012345678
    ```

The lineage store has no knowledge of the files produced — only that the tasks ran.
Render the graph for this run and you'll see it ends at the `TaskRun` nodes with no file nodes attached.

### Takeaway

In this section, you've seen:

- **Without `output {}` block:** The lineage store records `WorkflowRun` and `TaskRun` entries, but no `FileOutput` records. The provenance trail stops at the task level.
- **With `output {}` block:** Each published file becomes a `FileOutput` record with a checksum, size, and a direct link back to the task that produced it.

Before continuing, restore the `publish:` section and `output {}` block to `main.nf`.

---

## 4. Comparing runs with diff

The `diff` subcommand compares two lineage records to show exactly what changed between them.
This is most useful when you re-run a workflow after modifying inputs or parameters.

### 4.1. Modify the input data

Edit `data/greetings.csv` to replace `Hello` with `Hi`:

=== "After"

    ```console title="data/greetings.csv" hl_lines="1"
    Hi
    Bonjour
    Hola
    Ciao
    ```

=== "Before"

    ```console title="data/greetings.csv" hl_lines="1"
    Hello
    Bonjour
    Hola
    Ciao
    ```

### 4.2. Run the workflow again

```bash
nextflow run main.nf
```

??? success "Command output"

    ```console
     N E X T F L O W   ~  version 25.10.4

    Launching `main.nf` [modest_euler] DSL2 - revision: abc123def4

    executor >  local (8)
    [2b/3c4d5e] SAYHELLO (4)       [100%] 4 of 4, cached: 3 ✔
    [6f/7a8b9c] CONVERTTOUPPER (4) [100%] 4 of 4, cached: 3 ✔
    ```

Three tasks were cached (the greetings that didn't change) and one task in each process ran fresh for `Hi`.

List the runs to find the new LID:

```bash
nextflow lineage list
```

??? success "Command output"

    ```console
    RUN NAME         LID                              TIMESTAMP
    modest_euler     lid://a3b4c5d6e7f8901234567890   2026-06-15 10:10:00
    tender_turing    lid://a9b8c7d6e5f4321098765432   2026-06-15 10:05:00
    happy_dijkstra   lid://a1b2c3d4e5f6789012345678   2026-06-15 10:00:00
    ```

### 4.3. Compare two task runs

Use `nextflow lineage find` to locate the `SAYHELLO` task LIDs from the first run and the third run, then compare them:

```bash
nextflow lineage diff lid://b2c3d4e5f6789012345678a1 lid://b4c5d6e7f8901234567890a3
```

??? success "Command output"

    ```console
    TODO add real output
    ```

The `diff` identifies exactly what changed between the two task executions.
For tasks with many inputs or complex parameters, `diff` quickly pinpoints which values differed — without manually comparing JSON records.

### Takeaway

In this section, you've seen:

- **`nextflow lineage diff [LID-A] [LID-B]`:** Compare any two lineage records and surface exactly what changed.
- **Practical use:** When debugging a pipeline or verifying reproducibility, `diff` confirms whether two runs processed identical inputs — and if not, what differed.

---

## Summary

Data lineage in Nextflow gives you a queryable record of every workflow run, task execution, and published file.

### Key commands

| Command                                 | Purpose                               |
| --------------------------------------- | ------------------------------------- |
| `nextflow lineage list`                 | List all workflow runs in the store   |
| `nextflow lineage view [LID]`           | Inspect any record as JSON            |
| `nextflow lineage view [LID]#output`    | View published file records for a run |
| `nextflow lineage find`                 | Query records by type or field        |
| `nextflow lineage render [LID]`         | Generate an HTML lineage graph        |
| `nextflow lineage diff [LID-A] [LID-B]` | Compare two records                   |

### Key patterns

1. **Enable tracking** in `nextflow.config`:

   ```groovy
   lineage.enabled = true
   ```

2. **List all workflow runs** in the lineage store:

   ```bash
   nextflow lineage list
   ```

3. **Inspect any record** by LID:

   ```bash
   nextflow lineage view <LID>
   ```

4. **View published file records** for a run using the `#output` suffix:

   ```bash
   nextflow lineage view <LID>#output
   ```

5. **Query records by type** without knowing the LID:

   ```bash
   nextflow lineage find --type TaskRun
   ```

6. **Render a lineage graph** as HTML:

   ```bash
   nextflow lineage render <LID>
   ```

7. **Compare any two records** to surface exactly what changed between them:

   ```bash
   nextflow lineage diff <LID-A> <LID-B>
   ```

### Additional resources

- [Nextflow data lineage documentation](https://docs.seqera.io/nextflow/tutorials/data-lineage)

---

## What's next?

Return to the [menu of Side Quests](../index.md) or click the button in the bottom right of the page to move on to the next topic in the list.
