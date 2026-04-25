# Metadata and Meta Maps

In any scientific analysis, we rarely work with just the raw data files.
Each file comes with its own additional information: what it is, where it came from, and what makes it special.
This extra information is what we call metadata.

Metadata is data describing other data.
Metadata tracks important details about files and experimental conditions, and helps tailor analyses to each dataset's unique characteristics.

Think of it like a library catalog: while books contain the actual content (raw data), the catalog cards provide essential information about each book—when it was published, who wrote it, where to find it (metadata).
In Nextflow pipelines, metadata can be used to:

- Track file-specific information throughout the workflow
- Configure processes based on file characteristics
- Group related files for joint analysis

### Learning goals

In this side quest, we'll explore how to handle metadata in workflows.
Starting with a simple datasheet (often called a samplesheet in bioinformatics) containing basic file information, you'll learn how to:

- Read and parse file metadata from CSV files
- Understand why the "meta map + data file" interface is a widely used convention
- Add new metadata fields during workflow execution
- Use metadata to customize process behavior and organize outputs

These skills will help you build more robust and flexible pipelines that can handle complex file relationships and processing requirements.

### Prerequisites

Before taking on this side quest, you should:

- Have completed the [Hello Nextflow](../hello_nextflow/README.md) tutorial or equivalent beginner's course.
- Be comfortable using basic Nextflow concepts and mechanisms (processes, channels, operators)

---

## 0. Get started

#### Open the training codespace

If you haven't yet done so, make sure to open the training environment as described in the [Environment Setup](../envsetup/index.md).

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Move into the project directory

Move into the directory where the files for this tutorial are located.

```bash
cd side-quests/metadata
```

You can set VSCode to focus on this directory:

```bash
code .
```

The editor opens with the project directory in focus.

#### Review the materials

You'll find a main workflow file and a `data` directory containing a datasheet and a handful of data files.

??? abstract "Directory contents"

    ```console
    .
    ├── data
    │   ├── bonjour.txt
    │   ├── ciao.txt
    │   ├── guten_tag.txt
    │   ├── hallo.txt
    │   ├── hello.txt
    │   ├── hola.txt
    │   ├── salut.txt
    │   └── datasheet.csv
    ├── main.nf
    └── nextflow.config
    ```

The workflow in the `main.nf` file is a stub that you will gradually expand into a fully functioning workflow.

The datasheet lists the paths to the data files and some associated metadata, organized in 3 columns:

- `id`: self-explanatory, an ID given to the file
- `character`: a character name, that we will use later to draw different creatures
- `data`: paths to `.txt` files that contain greetings in different languages

```console title="datasheet.csv"
id,character,recording
sampleA,squirrel,/workspaces/training/side-quests/metadata/data/bonjour.txt
sampleB,tux,/workspaces/training/side-quests/metadata/data/guten_tag.txt
sampleC,sheep,/workspaces/training/side-quests/metadata/data/hallo.txt
sampleD,turkey,/workspaces/training/side-quests/metadata/data/hello.txt
sampleE,stegosaurus,/workspaces/training/side-quests/metadata/data/hola.txt
sampleF,moose,/workspaces/training/side-quests/metadata/data/salut.txt
sampleG,turtle,/workspaces/training/side-quests/metadata/data/ciao.txt
```

Each data file contains some greeting text in one of five languages (fr: French, de: German, es: Spanish, it: Italian, en: English).

We're going to use a tool called [`COWPY`](https://github.com/jeffbuttars/cowpy) to generate ASCII art of each character speaking its recorded greeting.

??? info "What does `COWPY` do?"

    `COWPY` is a command-line tool that generates ASCII art to display arbitrary text inputs in a fun way.
    It is a Python implementation of the classic [cowsay](https://en.wikipedia.org/wiki/Cowsay) tool by Tony Monroe.

    ```console
    cowpy "Hello Nextflow"
    ```

    ```console
    ______________________________________________________
    < Hello Nextflow >
    ------------------------------------------------------
        \   ^__^
          \  (oo)\_______
            (__)\       )\/\
              ||----w |
              ||     ||
    ```

    Optionally, you can select a character (or 'cowacter') to use instead of the default cow.

    ```console
    cowpy "Hello Nextflow" -c tux
    ```

    ```console
    __________________
    < Hello Nextflow >
    ------------------
      \
        \
            .--.
          |o_o |
          |:_/ |
          //   \ \
        (|     | )
        /'\_   _/`\
        \___)=(___/
    ```

In addition, we'll use a language analysis tool called `langid` to identify what language each character speaks and organize the outputs of the pipeline accordingly.

#### Review the assignment

Your challenge is to write a Nextflow workflow that will:

1. **Generate ASCII art** of each character
2. **Organize** outputs by language group (Germanic vs Romance languages)

This represents a typical workflow pattern where file-specific metadata drives processing decisions; exactly the kind of problem that metadata maps solve elegantly.

#### Readiness checklist

Think you're ready to dive in?

- [ ] I understand the goal of this course and its prerequisites
- [ ] My codespace is up and running
- [ ] I've set my working directory appropriately
- [ ] I understand the assignment

If you can check all the boxes, you're good to go.

---

## 1. Basic options for loading and using metadata

Open the `main.nf` workflow file to examine the workflow stub we're giving you as a starting point.

```groovy title="main.nf" linenums="1"
#!/usr/bin/env nextflow

workflow  {
    main:
    ch_datasheet = channel.fromPath("./data/datasheet.csv")
        .splitCsv(header: true)
        .view()

    publish:
    cowpy_art = channel.empty()
}

output {
    cowpy_art {
    }
}
```

The [`splitCsv`](https://www.nextflow.io/docs/latest/reference/operator.html#splitcsv) operator reads each row in the file as a channel element.
This is the same approach we use to load CSV data in Hello Nextflow, our beginner course.
Have a look at [this section](../../hello_nextflow/02_hello_channels/#4-read-input-values-from-a-csv-file) if you need a reminder of how that works.

With `header: true`, the first row is treated as column headers, so each element becomes a map of key-value pairs keyed by column name.

Note that since we're not running any processes on the data yet, the `publish` and `output` blocks are just stubs.

### 1.1. Run the workflow

Run the workflow to see how the channel contents are structured once everything is loaded:

```bash
nextflow run main.nf
```

??? success "Command output"

    ```console
     N E X T F L O W   ~  version 25.10.4

    Launching `main.nf` [exotic_albattani] DSL2 - revision: c0d03cec83

    [id:sampleA, character:squirrel, recording:/workspaces/training/side-quests/metadata/data/bonjour.txt]
    [id:sampleB, character:tux, recording:/workspaces/training/side-quests/metadata/data/guten_tag.txt]
    [id:sampleC, character:sheep, recording:/workspaces/training/side-quests/metadata/data/hallo.txt]
    [id:sampleD, character:turkey, recording:/workspaces/training/side-quests/metadata/data/hello.txt]
    [id:sampleE, character:stegosaurus, recording:/workspaces/training/side-quests/metadata/data/hola.txt]
    [id:sampleF, character:moose, recording:/workspaces/training/side-quests/metadata/data/salut.txt]
    [id:sampleG, character:turtle, recording:/workspaces/training/side-quests/metadata/data/ciao.txt]
    ```

As you can see, the operator has constructed a map of key-value pairs for each row in the CSV file, with the column headers as keys for the corresponding values.

Each map entry corresponds to a column in our datasheet:

- `id`
- `character`
- `recording`

This makes it easy to access specific fields from each row.
For example, we could access the file ID with `id` or the txt file path with `recording`.

??? info "(Optional) More about Groovy maps"

    In Groovy, the programming language that Nextflow is built on, a map is a key-value data structure similar to dictionaries in Python, objects in JavaScript, or hashes in Ruby.

    Here's a runnable script that shows how you can define a map and access its contents in practice:

    ```groovy title="examples/map_demo.nf"
    #!/usr/bin/env nextflow

    // Create a simple map
    def my_map = [id:'sampleA', character:'squirrel']

    // Print the whole map
    println "map: ${my_map}"

    // Access individual values using dot notation
    println "id: ${my_map.id}"
    println "character: ${my_map.character}"
    ```

    Even though it doesn't have a proper `workflow` block, Nextflow can run this as if it were a workflow:

    ```bash
    nextflow run examples/map_demo.nf
    ```

    And here's what you can expect to see in the output:

    ```console title="Output"
     N E X T F L O W   ~  version 25.10.4

    Launching `map_demo.nf` [cheesy_plateau] DSL2 - revision: fae5b8496e

    map: [id:sampleA, character:squirrel]
    id: sampleA
    character: squirrel
    ```

### 1.2. Pick out a specific field with `map`

We're going to use the `map` operator to iterate over each element in a channel and pick out just the character field, which we can access by name using dot notation.

#### 1.2.1. Add the map operation

To access the `character` column, add the `map` operation before the `.view()` operation as follows:

=== "After"

    ```groovy title="main.nf" linenums="5" hl_lines="3-5"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                row.character
            }
            .view()
    ```

=== "Before"

    ```groovy title="main.nf" linenums="5"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .view()
    ```

This way of accessing a specific field is explained in more detail in [this section](../../hello_nextflow/02_hello_channels/#43-use-the-map-operator-to-extract-the-greetings) of Hello Nextflow, if you need a reminder of how that works.

#### 1.2.2. Run the workflow

Run the workflow to verify that you can view the extracted character names.

```bash
nextflow run main.nf
```

??? success "Command output"

    ```console
     N E X T F L O W   ~  version 25.10.4

    Launching `main.nf` [exotic_albattani] DSL2 - revision: c0d03cec83

    squirrel
    tux
    sheep
    turkey
    stegosaurus
    moose
    turtle
    ```

This shows that we're able to access the values from the `character` column for each row.

Now let's do something with this data: use the `character` and `recording` fields together to generate ASCII art using `COWPY`.

### 1.3. Emit sub-channels with `multiMap`

We provide you with a pre-written `COWPY` process module, so first you need to examine the input requirements of the process.

You can open the file to see what the process looks like:

```groovy title="modules/cowpy.nf" linenums="1" hl_lines="7 8"
// Generate ASCII art with cowpy
process COWPY {

    container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

    input:
    path input_file
    val character

    output:
    path "cowpy-${input_file}"

    script:
    """
    cat ${input_file} | cowpy -c ${character} > cowpy-${input_file}
    """
}
```

As you can see, the process takes two separate inputs: a recording file and a character name.
Importantly, we have values for both, but they are currently bundled inside each element in the channel.

One way to extract multiple fields into separate channels is the [`multiMap`](https://www.nextflow.io/docs/latest/reference/operator.html#multimap) operator, which splits one channel into multiple named sub-channels in a single operation.

#### 1.3.1. Add the multiMap operation

Replace the `map` operation with `multiMap`:

=== "After"

    ```groovy title="main.nf" linenums="5" hl_lines="3-6"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .multiMap { row ->
                file: row.recording
                character: row.character
            }
    ```

=== "Before"

    ```groovy title="main.nf" linenums="5" hl_lines="3-6"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                row.character
            }
            .view()
    ```

The `multiMap` block defines two named sub-channels (`file` and `character`) from each row, which we can access as `ch_datasheet.file` and `ch_datasheet.character`.

#### 1.3.2. Call COWPY on the sub-channels

Now, include the `COWPY` process and give it each sub-channel as a separate argument:

=== "After"

    ```groovy title="main.nf" linenums="1" hl_lines="3 14"
    #!/usr/bin/env nextflow

    include { COWPY } from './modules/cowpy.nf'

    workflow {
        main:
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .multiMap { row ->
                file: row.recording
                character: row.character
            }

        COWPY(ch_datasheet.file, ch_datasheet.character)

        publish:
        cowpy_art = channel.empty()
    }

    output {
        cowpy_art {
        }
    }
    ```

=== "Before"

    ```groovy title="main.nf" linenums="1"
    #!/usr/bin/env nextflow

    workflow {
        main:
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .multiMap { row ->
                file: row.recording
                character: row.character
            }

        publish:
        cowpy_art = channel.empty()
    }

    output {
        cowpy_art {
        }
    }
    ```

This allows us to pass the two fields separately as `COWPY` requires.

#### 1.3.3. Set up the output publishing

Finally, add the output of `COWPY` to the `publish:` block:

=== "After"

    ```groovy title="main.nf" linenums="14" hl_lines="4"
        COWPY(ch_datasheet.file, ch_datasheet.character)

        publish:
        cowpy_art = COWPY.out
    ```

=== "Before"

    ```groovy title="main.nf" linenums="14" hl_lines="4"
        COWPY(ch_datasheet.file, ch_datasheet.character)

        publish:
        cowpy_art = channel.empty()
    ```

This will allow us to easily view the outputs produced by the workflow.

#### 1.3.4. Run the workflow

Run the workflow to check that `COWPY` runs on the inputs we provided:

```bash
nextflow run main.nf
```

??? success "Command output"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [clever_dijkstra] DSL2 - revision: a1b2c3d4e5

    executor >  local (7)
    [3a/f1c290] COWPY (7) [100%] 7 of 7 ✔
    ```

As you can see, `COWPY` ran on each file using the correct character for each one.

??? abstract "Results directory contents"

    ```console
    results/
    ├── cowpy-bonjour.txt
    ├── cowpy-ciao.txt
    ├── cowpy-guten_tag.txt
    ├── cowpy-hallo.txt
    ├── cowpy-hello.txt
    ├── cowpy-hola.txt
    └── cowpy-salut.txt
    ```

??? example "Content of results/cowpy-guten_tag.txt"

    ```console
    $ cat results/cowpy-guten_tag.txt
     _____________________________
    / Guten Tag, wie geht es dir? \
    \ Auf Wiedersehen, bis morgen /
     -----------------------------
       \
        \
            .--.
           |o_o |
           |:_/ |
          //   \ \
         (|     | )
        /'\_   _/`\
        \___)=(___/
    ```

This approach works, but has a limitation: we had to split the channel into two separate sub-channels.
If we wanted to pass more fields to the process, we'd need to split them out into more sub-channels.
That could get annoying and messy.

Good news: there is a simpler way to do this.

### 1.4. Group everything as a single input to the process

Rather than splitting the fields into separate channels, we can update the process to receive all inputs as a single tuple, which simplifies the call to the process.

#### 1.4.1. Update the COWPY process

Update `COWPY` to accept a tuple corresponding to the three elements in each row:

=== "After"

    ```groovy title="modules/cowpy.nf" linenums="1" hl_lines="7 10 14"
    // Generate ASCII art with cowpy
    process COWPY {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

        input:
        tuple val(id), val(character), path(recording)

        output:
        path "cowpy-${recording}"

        script:
        """
        cat ${recording} | cowpy -c ${character} > cowpy-${recording}
        """
    }
    ```

=== "Before"

    ```groovy title="modules/cowpy.nf" linenums="1" hl_lines="7-8 11 15"
    // Generate ASCII art with cowpy
    process COWPY {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

        input:
        path input_file
        val character

        output:
        path "cowpy-${input_file}"

        script:
        """
        cat ${input_file} | cowpy -c ${character} > cowpy-${input_file}
        """
    }
    ```

Now the process takes just one input containing all the things we might want to give it.

#### 1.4.2. Use `map()` to create the input tuple

We still need to use a mapping operation to enumerate the elements we want to pass in the tuple to the process:

=== "After"

    ```groovy title="main.nf" linenums="5" hl_lines="3-5"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [row.id, row.character, row.recording]
            }
    ```

=== "Before"

    ```groovy title="main.nf" linenums="5" hl_lines="3-6"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .multiMap { row ->
                file: row.recording
                character: row.character
            }
    ```

You might wonder why we can't just pass the entire Groovy map coming from `splitCsv` as is.
It's because we need to tell Nextflow explicitly that the recording file needs to be handled as a path (i.e. it needs to be staged properly).
That happens at the level of `COWPY`'s input interface, where the `recording` element is explicitly designated as a `path`.

#### 1.4.3. Update the call to the process

Finally, let's replace the two separate inputs in the process call with the single tuple we just created:

=== "After"

    ```groovy title="main.nf" linenums="5" hl_lines="7"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [row.id, row.character, row.recording]
            }

        COWPY(ch_datasheet)
    ```

=== "Before"

    ```groovy title="main.nf" linenums="5" hl_lines="7"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [row.id, row.character, row.recording]
            }

        COWPY(ch_datasheet.file, ch_datasheet.character)
    ```

This simplifies the call to the process a little.

#### 1.4.4. Run the workflow

Run the workflow to verify that `COWPY` can still process the data correctly:

```bash
nextflow run main.nf
```

??? success "Command output"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [pedantic_lovelace] DSL2 - revision: b2c3d4e5f6

    executor >  local (7)
    [5e/2a1b34] COWPY (7) [100%] 7 of 7 ✔
    ```

The output is the same seven `cowpy-*.txt` files as before, now produced with a simpler call to `COWPY`.

??? abstract "Results directory contents"

    ```console
    results/
    ├── cowpy-bonjour.txt
    ├── cowpy-ciao.txt
    ├── cowpy-guten_tag.txt
    ├── cowpy-hallo.txt
    ├── cowpy-hello.txt
    ├── cowpy-hola.txt
    └── cowpy-salut.txt
    ```

??? example "Content of results/cowpy-guten_tag.txt"

    ```console
    $ cat results/cowpy-guten_tag.txt
     _____________________________
    / Guten Tag, wie geht es dir? \
    \ Auf Wiedersehen, bis morgen /
     -----------------------------
       \
        \
            .--.
           |o_o |
           |:_/ |
          //   \ \
         (|     | )
        /'\_   _/`\
        \___)=(___/
    ```

This is a slight improvement over the `multiMap` approach.
But we still had to unpack the original Groovy map to create the input tuple, and there's a tight coupling between the process and the datasheet: the `COWPY` input definition now references the column names `id`, `character`, and `recording` directly.

```groovy
input:
tuple val(id), val(character), path(recording)
```

If a collaborator uses a differently structured datasheet — with additional columns, or columns in a different order — this process won't work without modification.
This makes the process fragile, because its input structure is tied to the exact composition of the datasheet.

To solve this, we need a way to pass all the metadata as a bundle without hard-coding its exact structure into the process interface.

### 1.5. Use a meta map + file interface

The solution is to separate two distinct concerns in the channel: the **metadata about a sample**, and the **data file** itself.
By bundling all metadata into a single map — the "meta map" — we get a consistent two-element tuple regardless of how many metadata columns the datasheet contains:

```groovy title="Syntax example"
input:
tuple val(meta), path(file)
```

Adding or removing columns from the datasheet changes what's inside `meta`, but the tuple shape `[meta, file]` stays constant.
Processes that accept this structure don't need to know or care how many metadata fields exist.

#### 1.5.1. Re-organize the tuple content into a meta map

Let's restructure the `map` operation to produce a `[meta, file]` tuple:

=== "After"

    ```groovy title="main.nf" linenums="5" hl_lines="4 6 8 11"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [[id: row.id, character: row.character], row.recording]
            }
            .view()

        // COWPY(ch_datasheet)  // Will update in the next step

        publish:
        cowpy_art = channel.empty() // COWPY.out
    ```

=== "Before"

    ```groovy title="main.nf" linenums="5" hl_lines="4 7"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [row.id, row.character, row.recording]
            }

        COWPY(ch_datasheet)

        publish:
        cowpy_art = COWPY.out
    ```

You'll notice that we also added a `view()` statement, commented out the `COWPY` call and replaced `COWPY.out` with `channel.empty()` because the process input definition doesn't match the new structure yet.

#### 1.5.2. Run the workflow to inspect the re-organized content

Run the workflow to see the new channel shape:

```bash
nextflow run main.nf
```

??? success "Command output"

    ```console title="View meta map"
     N E X T F L O W   ~  version 25.10.4

    Launching `main.nf` [lethal_booth] DSL2 - revision: 0d8f844c07

    [[id:sampleA, character:squirrel], /workspaces/training/side-quests/metadata/data/bonjour.txt]
    [[id:sampleB, character:tux], /workspaces/training/side-quests/metadata/data/guten_tag.txt]
    [[id:sampleC, character:sheep], /workspaces/training/side-quests/metadata/data/hallo.txt]
    [[id:sampleD, character:turkey], /workspaces/training/side-quests/metadata/data/hello.txt]
    [[id:sampleE, character:stegosaurus], /workspaces/training/side-quests/metadata/data/hola.txt]
    [[id:sampleF, character:moose], /workspaces/training/side-quests/metadata/data/salut.txt]
    [[id:sampleG, character:turtle], /workspaces/training/side-quests/metadata/data/ciao.txt]
    ```

Each element in the channel is now a two-element tuple: the meta map first, the file second.

```console title="Example element structure"
[
  [id:sampleA, character:squirrel],
  /workspaces/training/side-quests/metadata/data/bonjour.txt
]
```

If we later add a `language` column to the datasheet, it will become available as `meta.language` without requiring any changes to the process input definition.

#### 1.5.3. Update the `COWPY` process to use the meta map

Update `COWPY` to accept the `[meta, file]` tuple structure:

=== "After"

    ```groovy title="modules/cowpy.nf" linenums="1" hl_lines="7 10 14"
    // Generate ASCII art with cowpy
    process COWPY {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

        input:
        tuple val(meta), path(input_file)

        output:
        path "cowpy-${input_file}"

        script:
        """
        cat ${input_file} | cowpy -c ${meta.character} > cowpy-${input_file}
        """
    }
    ```

=== "Before"

    ```groovy title="modules/cowpy.nf" linenums="1" hl_lines="7 10 14"
    // Generate ASCII art with cowpy
    process COWPY {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

        input:
        tuple val(id), val(character), path(recording)

        output:
        path "cowpy-${recording}"

        script:
        """
        cat ${recording} | cowpy -c ${character} > cowpy-${recording}
        """
    }
    ```

Inside the script block, `meta.character` accesses the `character` field from the meta map.
Any field in the meta map is accessible the same way.

#### 1.5.4. Update the process call

Restore the `COWPY` call and connect its output for publishing:

=== "After"

    ```groovy title="main.nf" linenums="5" hl_lines="7 10"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [[id: row.id, character: row.character], row.recording]
            }

        COWPY(ch_datasheet)

        publish:
        cowpy_art = COWPY.out
    ```

=== "Before"

    ```groovy title="main.nf" linenums="5" hl_lines="6 8 11"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [[id: row.id, character: row.character], row.recording]
            }
            .view()

        // COWPY(ch_datasheet)  // Will update in the next step

        publish:
        cowpy_art = channel.empty() // COWPY.out
    ```

We've also restored the output publishing.

#### 1.5.5. Run the workflow

Run the workflow to check that it all works:

```bash
nextflow run main.nf
```

??? success "Command output"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [wise_sammet] DSL2 - revision: 99797b1e92

    executor >  local (7)
    [5d/dffd4e] COWPY (7) [100%] 7 of 7 ✔
    ```

The results directory now contains the ASCII art files.

??? abstract "Directory contents"

    ```console
    results/
    ├── cowpy-bonjour.txt
    ├── cowpy-ciao.txt
    ├── cowpy-guten_tag.txt
    ├── cowpy-hallo.txt
    ├── cowpy-hello.txt
    ├── cowpy-hola.txt
    └── cowpy-salut.txt
    ```

??? example "Content of results/cowpy-guten_tag.txt"

    ```console
    $ cat results/cowpy-guten_tag.txt
     _____________________________
    / Guten Tag, wie geht es dir? \
    \ Auf Wiedersehen, bis morgen /
     -----------------------------
       \
        \
            .--.
           |o_o |
           |:_/ |
          //   \ \
         (|     | )
        /'\_   _/`\
        \___)=(___/
    ```

The process now receives all metadata as a bundle via `meta`, uses what it needs (`meta.character`), and ignores the rest.

This is the standard interface used by all [nf-core](https://nf-co.re/) modules.
The `tuple val(meta), path(file)` pattern appears consistently throughout the nf-core module library, which is why workflows that adopt this convention can swap in nf-core modules with minimal friction.

### Takeaway

In this section, you've learned:

- **How to read in datasheets:** Using `splitCsv` to parse CSV files with header information
- **Why the meta map convention exists:** Separating metadata from data files into `[meta, file]` tuples keeps the channel structure stable as the datasheet evolves
- **How to use meta map fields inside a process:** Any field in the meta map is accessible via dot notation in the script block

---

## 2. Additional metadata manipulations

Now that the meta map interface is in place, we can enrich it as data flows through the pipeline.

We're going to use a tool called [`langid`](https://github.com/saffsd/langid.py) to identify the language in each recording file.
Given a snippet of text, it outputs a language prediction and a probability score to `stdout`.

### 2.1. Add a language identification step

We provide a pre-written process module called `IDENTIFY_LANGUAGE` that wraps the `langid` tool.

Open the module file to examine its code:

```groovy title="modules/langid.nf" linenums="1" hl_lines="7 10"
// Use langid to predict the language of each input file
process IDENTIFY_LANGUAGE {

    container 'community.wave.seqera.io/library/pip_langid:b2269f456a5629ff'

    input:
    tuple val(meta), path(file)

    output:
    tuple val(meta), path(file), stdout

    script:
    """
    langid < ${file} -l en,de,fr,es,it | sed -E "s/.*\\('([a-z]+)'.*/\\1/" | tr -d '\\n'
    """
}
```

The input definition uses the same `tuple val(meta), path(file)` structure we just built in Section 1, so `ch_datasheet` can feed directly into this process without any adaptation.

The output adds `stdout` as a third element: this captures the language prediction that `langid` prints to the console.
The `sed` command strips the probability score and trailing newline, leaving just the two-letter language code.

#### 2.1.1. Add a call to `IDENTIFY_LANGUAGE`

Include the `IDENTIFY_LANGUAGE` process module and call it on the datasheet channel:

=== "After"

    ```groovy title="main.nf" linenums="1" hl_lines="4 14-16"
    #!/usr/bin/env nextflow

    include { COWPY } from './modules/cowpy.nf'
    include { IDENTIFY_LANGUAGE } from './modules/langid.nf'

    workflow {
        main:
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [[id: row.id, character: row.character], row.recording]
            }

        // Run langid to identify the language of each greeting
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out.view()

        COWPY(ch_datasheet)
    ```

=== "Before"

    ```groovy title="main.nf" linenums="1"
    #!/usr/bin/env nextflow

    include { COWPY } from './modules/cowpy.nf'

    workflow {
        main:
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [[id: row.id, character: row.character], row.recording]
            }

        COWPY(ch_datasheet)
    ```

The main output of this process is just a string, so there are no output files to publish.
Instead, we use `IDENTIFY_LANGUAGE.out.view()` to view the results of the operation.

#### 2.1.2. Run the workflow

Run the workflow to produce the language identification, using `-resume` to avoid re-running the `COWPY` tasks:

```bash
nextflow run main.nf -resume
```

??? success "Command output"

    ```console
     N E X T F L O W   ~  version 25.10.4

    Launching `main.nf` [voluminous_mcnulty] DSL2 - revision: f9bcfebabb

    executor >  local (14)
    [5d/dffd4e] COWPY (7)             [100%] 7 of 7, cached: 7 ✔
    [4e/f722fe] IDENTIFY_LANGUAGE (7) [100%] 7 of 7 ✔
    [[id:sampleA, character:squirrel], /workspaces/training/side-quests/metadata/work/eb/f7148ebdd898fbe1136bec6a714acb/bonjour.txt, fr]
    [[id:sampleB, character:tux], /workspaces/training/side-quests/metadata/work/16/71d72410952c22cd0086d9bca03680/guten_tag.txt, de]
    [[id:sampleD, character:turkey], /workspaces/training/side-quests/metadata/work/c4/b7562adddc1cc0b7d414ec45d436eb/hello.txt, en]
    [[id:sampleC, character:sheep], /workspaces/training/side-quests/metadata/work/ea/04f5d979429e4455e14b9242fb3b45/hallo.txt, de]
    [[id:sampleF, character:moose], /workspaces/training/side-quests/metadata/work/5a/6c2b84bf8fadb98e28e216426be079/salut.txt, fr]
    [[id:sampleE, character:stegosaurus], /workspaces/training/side-quests/metadata/work/af/ee7c69bcab891c40d0529305f6b9e7/hola.txt, es]
    [[id:sampleG, character:turtle], /workspaces/training/side-quests/metadata/work/4e/f722fe47271ba7ebcd69afa42964ca/ciao.txt, it]
    ```

We now have a language prediction for each file in the dataset.

Note that the output tuple is composed of `[meta, file, lang_id]`, meaning the meta map and file are carried through alongside the new result.

!!! note

    This pattern of keeping the meta map associated with results makes it easier to join results across channels later on.
    You can't rely on the order of items in channels to associate data correctly.
    You must use keys instead.
    Meta maps provide an ideal structure for this purpose.

    This use case is explored in detail in the [Splitting & Grouping](../splitting_and_grouping/) side quest.

### 2.2. Augment metadata with process outputs

The language prediction is itself metadata about the data in the file.
Rather than keeping it as a separate element, let's fold it back into the meta map.

#### 2.2.1. Create a new and expanded meta map

We can create a new meta map to replace the original one using the Groovy `+` operator:

=== "After"

    ```groovy title="main.nf" linenums="14" hl_lines="3-7"
        // Run langid to identify the language of each greeting
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out
            .map { meta, file, lang_id ->
                [meta + [lang: lang_id], file]
            }
            .view()
    ```

=== "Before"

    ```groovy title="main.nf" linenums="14" hl_lines="3"
        // Run langid to identify the language of each greeting
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out.view()
    ```

The heart of this operation is `#!groovy meta + [lang: lang_id]`.

That code essentially creates a temporary map with a single key-value pair containing the (`[lang: lang_id]`), then uses the Groovy `+` operator to combine it with the original `meta` map containing the pre-existing metadata, producing a new and expanded meta map.

For a more detailed explanation, see the box below.

??? info "Creation of the new meta map using the `+` operator"

    **First, you need to know that we can merge the contents of two maps using the Groovy operator `+`.**

    Let's say we have the following maps:

    ```groovy
    map1 = [id: 'sampleA', character: 'squirrel']
    map2 = [lang: 'fr']
    ```

    We can merge them like this:

    ```groovy
    new_map = map1 + map2
    ```

    The contents of `new_map` will be:

    ```groovy
    [id: 'sampleA', character: 'squirrel', lang: 'fr']
    ```

    Great!

    **But what if you need to add a field that's not already part of a map?**

    Let's say you start again from `map1`, but the language prediction is not in its own map (there is no `map2`).
    Instead, it is held in a variable called `lang_id`, and you know you want to store its value (`'fr'`) with the key `lang`.

    You can actually do the following:

    ```groovy
    new_map = map1 + [lang: lang_id]
    ```

    Here, `[lang: lang_id]` creates a new unnamed map on the fly, and `map1 + ` merges `map1` with the new unnamed map, producing the same `new_map` contents as before.

    Neat, right?

    **Now let's transpose that into the context of a Nextflow `channel.map()` operation.**

    The code becomes:

    ```groovy
    .map { map1, lang_id ->
        map1 + [lang: lang_id]
    }
    ```

    This does the following:

    - `map1, lang_id ->` takes the two items in the tuple
    - `map1 + [lang: lang_id]` creates the new map as detailed above

    The output is a single unnamed map with the same contents as `new_map` in our example above.
    So we've effectively transformed:

    ```groovy
    [id: 'sampleA', character: 'squirrel'], 'fr'
    ```

    into:

    ```groovy
    [id: 'sampleA', character: 'squirrel', lang: 'fr']
    ```

    Hopefully you can see that if we change `map1` to `meta`, that's basically all we need in order to add the language predication to our meta map in our workflow.

    Except for one thing!

    In the case of our workflow, **we also need to account for the presence of the `file` object in the tuple**, which is composed of `meta, file, lang_id`.

    So the code here would become:

    ```groovy
    .map { meta, file, lang_id ->
        [meta + [lang: lang_id], file]
    }
    ```

    If you're having a hard time following why the `file` seems to be moving around in the `map` operation, imagine that instead of `[meta + [lang: lang_id], file]`, that line reads `[new_map, file]`.
    This should make it more clear that we're simply leaving the `file` in its original place in second position in the tuple. We've just taken the `new_info` value and folded it into the map that's in first position.

    **And this brings us back to the `tuple val(meta), path(file)` channel structure!**

#### 2.2.2. Run the workflow

Once you're confident you understand what the code is doing, run the workflow to see if it works:

```bash
nextflow run main.nf -resume
```

??? success "Command output"

    ```console
     N E X T F L O W   ~  version 25.10.4

    Launching `main.nf` [cheeky_fermat] DSL2 - revision: d096281ee4

    [5d/dffd4e] COWPY (7)             [100%] 7 of 7, cached: 7 ✔
    [4e/f722fe] IDENTIFY_LANGUAGE (7) [100%] 7 of 7, cached: 7 ✔
    [[id:sampleA, character:squirrel, lang:fr], /workspaces/training/side-quests/metadata/work/eb/f7148ebdd898fbe1136bec6a714acb/bonjour.txt]
    [[id:sampleB, character:tux, lang:de], /workspaces/training/side-quests/metadata/work/16/71d72410952c22cd0086d9bca03680/guten_tag.txt]
    [[id:sampleC, character:sheep, lang:de], /workspaces/training/side-quests/metadata/work/ea/04f5d979429e4455e14b9242fb3b45/hallo.txt]
    [[id:sampleD, character:turkey, lang:en], /workspaces/training/side-quests/metadata/work/c4/b7562adddc1cc0b7d414ec45d436eb/hello.txt]
    [[id:sampleF, character:moose, lang:fr], /workspaces/training/side-quests/metadata/work/5a/6c2b84bf8fadb98e28e216426be079/salut.txt]
    [[id:sampleE, character:stegosaurus, lang:es], /workspaces/training/side-quests/metadata/work/af/ee7c69bcab891c40d0529305f6b9e7/hola.txt]
    [[id:sampleG, character:turtle, lang:it], /workspaces/training/side-quests/metadata/work/4e/f722fe47271ba7ebcd69afa42964ca/ciao.txt]
    ```

Yep, that checks out!
We've neatly reorganized the output of the process from `meta, file, lang_id` so that `lang_id` is now one of the keys in the meta map, and the channel's tuples fit the `meta, file` model once again.

!!! tip "Removing keys from a meta map"

    You can remove a key from a meta map using the Groovy [`subMap`](https://docs.groovy-lang.org/latest/html/groovy-jdk/java/util/Map.html#subMap(java.util.Collection)) method, which returns a new map containing only the keys you specify:

    ```groovy
    meta.subMap(['id', 'character'])  // returns a map with only 'id' and 'character'
    ```

    This is useful when a downstream process or module does not need all the fields that have accumulated in the meta map.

### 2.3. Assign a language group using conditionals

With the language prediction in the meta map, we can derive further metadata from it.
The languages in our dataset fall into two families: Germanic (English, German) and Romance (French, Spanish, Italian).
Adding a `lang_group` field will make that classification available downstream.

#### 2.3.1. Add a `map` operation with the conditional logic

We're going to use a second `map` operation with conditional logic to assign the language family:

```groovy
.map { meta, file ->

    // conditional logic defining lang_group goes here

    [meta + [lang_group: lang_group], file]
}
```

Here's the logic to apply:

- Start with `lang_group = 'unknown'` as the default.
- If `meta.lang` is `'de'` or `'en'`, set `lang_group` to `'germanic'`.
- Else if `meta.lang` is in `['fr', 'es', 'it']`, set `lang_group` to `'romance'`.

!!! tip

    You can access the value of `lang` within the map operation with `meta.lang`.

Make the following changes to the workflow:

=== "After"

    ```groovy title="main.nf" linenums="14" hl_lines="7-19 21"
        // Run langid to identify the language of each greeting
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out
            .map { meta, file, lang_id ->
                [meta + [lang: lang_id], file]
            }
            .map { meta, file ->

                def lang_group = "unknown"
                if (meta.lang.equals("de") || meta.lang.equals('en')) {
                    lang_group = "germanic"
                }
                else if (meta.lang in ["fr", "es", "it"]) {
                    lang_group = "romance"
                }

                [meta + [lang_group: lang_group], file]
            }
            .set { ch_languages }

        ch_languages.view()
    ```

=== "Before"

    ```groovy title="main.nf" linenums="14" hl_lines="7"
        // Run langid to identify the language of each greeting
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out
            .map { meta, file, lang_id ->
                [meta + [lang: lang_id], file]
            }
            .view()
    ```

Key points:

- `def lang_group = "unknown"` initializes the variable with a safe default.
- The `if / else if` structure handles the two language families; anything else stays `'unknown'`.
- `#!groovy .set { ch_languages }` gives the resulting channel a name for use in the next step.

<!-- TODO (future) Add note/links to relevant docs in additional resources section -->

#### 2.3.2. Run the workflow:

Run the workflow to verify that it works:

```bash
nextflow run main.nf -resume
```

??? success "Command output"

    ```console
     N E X T F L O W   ~  version 25.10.4

    Launching `main.nf` [wise_almeida] DSL2 - revision: 46778c3cd0

    [5d/dffd4e] COWPY (7)             [100%] 7 of 7, cached: 7 ✔
    [da/652cc6] IDENTIFY_LANGUAGE (7) [100%] 7 of 7, cached: 7 ✔
    [[id:sampleA, character:squirrel, lang:fr, lang_group:romance], /workspaces/training/side-quests/metadata/data/bonjour.txt]
    [[id:sampleB, character:tux, lang:de, lang_group:germanic], /workspaces/training/side-quests/metadata/data/guten_tag.txt]
    [[id:sampleC, character:sheep, lang:de, lang_group:germanic], /workspaces/training/side-quests/metadata/data/hallo.txt]
    [[id:sampleD, character:turkey, lang:en, lang_group:germanic], /workspaces/training/side-quests/metadata/data/hello.txt]
    [[id:sampleE, character:stegosaurus, lang:es, lang_group:romance], /workspaces/training/side-quests/metadata/data/hola.txt]
    [[id:sampleF, character:moose, lang:fr, lang_group:romance], /workspaces/training/side-quests/metadata/data/salut.txt]
    [[id:sampleG, character:turtle, lang:it, lang_group:romance], /workspaces/training/side-quests/metadata/data/ciao.txt]
    ```

The meta map now carries four fields: `id`, `character`, `lang`, and `lang_group`.
The channel structure is still `[meta, file]`.

### 2.4. Use metadata to name and organize outputs

With `lang` and `lang_group` now available in the meta map, we can use them to give the output files meaningful names and organize them into subdirectories by language family.

This requires three changes: updating the `COWPY` process to rename its output and include `meta` in what it emits, wiring `ch_languages` into `COWPY`, and updating the output block to specify the subdirectory path.

#### 2.4.1. Update the `COWPY` process

Rename the output file using the language code from the meta map, and add `meta` to the output so the output block can access `lang_group` for subdirectory routing:

=== "After"

    ```groovy title="modules/cowpy.nf" linenums="9" hl_lines="2 6"
        output:
        tuple val(meta), path("${meta.lang}-${input_file}")

        script:
        """
        cat ${input_file} | cowpy -c ${meta.character} > ${meta.lang}-${input_file}
        """
    ```

=== "Before"

    ```groovy title="modules/cowpy.nf" linenums="9" hl_lines="2 6"
        output:
        path "cowpy-${input_file}"

        script:
        """
        cat ${input_file} | cowpy -c ${meta.character} > cowpy-${input_file}
        """
    ```

This shows how we can take advantage of other metadata fields to customize the behavior of a process, without having to modify the input definition at all.

#### 2.4.2. Wire `ch_languages` into `COWPY` and update the output block

Replace `COWPY(ch_datasheet)` with `COWPY(ch_languages)` and update the `output {}` block to route each file into its language group subdirectory:

=== "After"

    ```groovy title="main.nf" linenums="36" hl_lines="1 7-9"
        COWPY(ch_languages)

        publish:
        cowpy_art = COWPY.out
    }

    output {
        cowpy_art {
            path { meta, file -> meta.lang_group }
        }
    }
    ```

=== "Before"

    ```groovy title="main.nf" linenums="36" hl_lines="1 7-9"
        COWPY(ch_datasheet)

        publish:
        cowpy_art = COWPY.out
    }

    output {
        cowpy_art {
        }
    }
    ```

Also remove the temporary `ch_languages.view()` line.

#### 2.4.3. Run the full pipeline

Delete the previous results and run the full pipeline:

```bash
rm -r results
nextflow run main.nf
```

??? success "Command output"

    ```console
     N E X T F L O W   ~  version 25.10.4

    Launching `main.nf` [suspicious_crick] DSL2 - revision: 25541014c5

    executor >  local (14)
    [5d/dffd4e] IDENTIFY_LANGUAGE (7) [100%] 7 of 7 ✔
    [e7/317c18] COWPY (7)             [100%] 7 of 7 ✔
    ```

The results directory is now organized by language family, with each file named after its detected language:

```console title="Results directory contents"
results/
├── germanic
│   ├── de-guten_tag.txt
│   ├── de-hallo.txt
│   └── en-hello.txt
└── romance
    ├── es-hola.txt
    ├── fr-bonjour.txt
    ├── fr-salut.txt
    └── it-ciao.txt
```

The `path` closure in the `output {}` block receives each `[meta, file]` tuple and returns `meta.lang_group` as the subdirectory name.
The file name itself comes from what the process outputs (`#!groovy "${meta.lang}-${input_file}"`).
Both pieces of metadata (language code and language group) come from the enriched meta map built up in this section.

### Takeaway

In this section, you've learned:

- **How to augment the meta map with process outputs:** Adding new keys with `#!groovy meta + [key: value]` keeps the `[meta, file]` channel structure intact while enriching the metadata.
- **How to derive metadata from metadata:** Conditional logic inside a `map` operation can compute new fields from existing ones.
- **How to use metadata for output organization:** The `path` closure in the `output {}` block can read from the meta map to route files into subdirectories.

---

## 3. Robustness considerations

When metadata values drive process behavior, missing or incomplete data can cause problems that are difficult to diagnose.
Here's what to expect and how to handle it.

### 3.1. What happens when a required metadata field is missing

The `character` value is required for the `COWPY` process to produce a valid result.
The failure mode depends on whether the column exists in the datasheet but is empty, or is absent entirely.

#### 3.1.1. The column exists but a value is empty

Suppose one entry in the datasheet has a blank `character` field:

```csv title="datasheet.csv" linenums="1" hl_lines="2"
id,character,recording
sampleA,,/workspaces/training/side-quests/metadata/data/bonjour.txt
sampleB,tux,/workspaces/training/side-quests/metadata/data/guten_tag.txt
...
```

The `character` key is created for all entries when the datasheet is parsed, but `meta.character` for `sampleA` will be an empty string.
When Nextflow substitutes `#!groovy ${meta.character}` into the command, the `COWPY` tool receives an empty argument for `-c` and fails:

??? failure "Command output"

    ```console hl_lines="8 11 16 28"
     N E X T F L O W   ~  version 25.10.4

    Launching `main.nf` [marvelous_hirsch] DSL2 - revision: 0dfeee3cc1

    executor >  local (9)
    [c1/c5dd4f] process > IDENTIFY_LANGUAGE (7) [ 85%] 6 of 7
    [d3/b7c415] process > COWPY (2)             [  0%] 0 of 6
    ERROR ~ Error executing process > 'COWPY (1)'

    Caused by:
      Process `COWPY (1)` terminated with an error exit status (2)


    Command executed:

      cat bonjour.txt | cowpy -c  > fr-bonjour.txt

    Command exit status:
      2

    Command output:
      (empty)

    Command error:
      usage: cowpy [-h] [-l] [-L] [-t] [-u] [-e EYES] [-c COWACTER] [-E] [-r] [-x]
                  [-C]
                  [msg ...]
      cowpy: error: argument -c/--cowacter: expected one argument

    Work dir:
      /workspaces/training/side-quests/metadata/work/ca/9d49796612a54dec5ed466063c809b

    Container:
      community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273

    Tip: you can try to figure out what's wrong by changing to the process work dir and showing the script file named `.command.sh`

    -- Check '.nextflow.log' file for details
    ```

The error message (`expected one argument`) points to the empty `-c` flag.
Checking the work directory's `.command.sh` file confirms the command was run with an empty value.

#### 3.1.2. The column does not exist in the datasheet

If the `character` column is absent entirely:

```csv title="datasheet.csv" linenums="1"
id,recording
sampleA,/workspaces/training/side-quests/metadata/data/bonjour.txt
sampleB,/workspaces/training/side-quests/metadata/data/guten_tag.txt
...
```

The `character` key is never created in the meta map.
When the process script evaluates `#!groovy ${meta.character}`, the missing key returns `null`, and Nextflow literally substitutes the string `null` into the command:

??? failure "Command output"

    ```console hl_lines="8 11 16"
     N E X T F L O W   ~  version 25.10.4

    Launching `main.nf` [jovial_bohr] DSL2 - revision: eaaf375827

    executor >  local (9)
    [0d/ada9db] process > IDENTIFY_LANGUAGE (5) [ 85%] 6 of 7
    [06/28065f] process > COWPY (2)             [  0%] 0 of 6
    ERROR ~ Error executing process > 'COWPY (2)'

    Caused by:
      Process `COWPY (2)` terminated with an error exit status (1)


    Command executed:

      cat guten_tag.txt | cowpy -c null > de-guten_tag.txt

    Command exit status:
      1

    Command output:
      (empty)

    Command error:
      Traceback (most recent call last):
        File "/opt/conda/bin/cowpy", line 10, in <module>
          sys.exit(main())
                  ~~~~^^
        File "/opt/conda/lib/python3.13/site-packages/cowpy/cow.py", line 1215, in main
          print(cow(eyes=args.eyes,
                ~~~^^^^^^^^^^^^^^^^
                tongue=args.tongue,
                ^^^^^^^^^^^^^^^^^^^
                thoughts=args.thoughts
                ^^^^^^^^^^^^^^^^^^^^^^
                    ).milk(msg)
                    ^
      TypeError: 'str' object is not callable

    Work dir:
      /workspaces/training/side-quests/metadata/work/06/28065f7d9fd7d22bba084aa941b6d6

    Container:
      community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273

    Tip: you can replicate the issue by changing to the process work dir and entering the command `bash .command.run`

    -- Check '.nextflow.log' file for details
    ```

The `cowpy -c null` in the executed command is the diagnostic clue.

### 3.2. Strategies for handling missing metadata

There are two complementary approaches to make workflows more robust against missing metadata.

**1. Input validation**

The most reliable solution is to validate the datasheet before any processing begins, so problems are caught early with a clear error message rather than surfacing as a cryptic process failure mid-run.
The [Hello nf-core](../hello_nf-core/05_input_validation.md) training covers how to add input validation using the nf-schema plugin. <!-- TODO (future) pending a proper Validation side quest -->

**2. Explicit process inputs for required values**

If you want the process interface itself to communicate that a particular value is mandatory, consider extracting it from the meta map as an explicit input:

=== "Process definition"

    ```groovy title="modules/cowpy.nf" linenums="6"
    input:
    tuple val(meta), val(character), path(input_file)
    ```

=== "Workflow call"

    ```groovy title="main.nf"
    COWPY(ch_languages.map { meta, file -> [meta, meta.character, file] })
    ```

This approach makes `character` a visible, required part of the process contract.
Anyone reading the module can immediately see that a character value must be provided.
If the field is absent, the workflow fails clearly at the channel level before the process even runs.

This highlights a useful design principle:

**Use the meta map for optional or descriptive information; extract required values as explicit inputs.**

The meta map keeps channel structures clean and stable, but for values that are genuinely required by a process, surfacing them as named inputs improves clarity and makes the module easier to use correctly in other contexts.

### Takeaway

In this section, you've seen:

- **How missing metadata manifests:** An empty field produces an empty argument; an absent field produces `null` substituted literally into the command.
- **Two complementary strategies:** Input validation to catch problems early, and explicit process inputs to communicate requirements clearly.

---

## Summary

In this side quest, you've explored how to effectively work with metadata in Nextflow workflows.

The "meta map + data file" tuple pattern is a core convention in Nextflow, offering several advantages over passing metadata as individual values:

- The channel structure stays stable as the datasheet evolves
- Process behavior can be customized per sample without hard-coding field names
- Metadata is available throughout the pipeline for naming, grouping, and organizing outputs
- Modules written to this interface are interchangeable, including nf-core modules

### Key patterns

1.  **Reading and structuring metadata:** Parse a CSV datasheet and create a meta map.

    ```groovy
    channel.fromPath('datasheet.csv')
        .splitCsv(header: true)
        .map { row ->
            [ [id: row.id, character: row.character], row.recording ]
        }
    ```

2.  **Expanding metadata during the workflow:** Add new keys from process outputs or derived logic.

    ```groovy
    // From a process output
    .map { meta, file, lang ->
        [ meta + [lang: lang], file ]
    }

    // From conditional logic
    .map { meta, file ->
        def lang_group = "unknown"
        if (meta.lang in ["de", "en"]) { lang_group = "germanic" }
        else if (meta.lang in ["fr", "es", "it"]) { lang_group = "romance" }
        [ meta + [lang_group: lang_group], file ]
    }
    ```

3.  **Using metadata inside a process:** Access any field via dot notation in the script block.

    ```groovy
    cat ${input_file} | cowpy -c ${meta.character} > ${meta.lang}-${input_file}
    ```

4.  **Organizing outputs by metadata value:** Use the `path` closure in the `output {}` block.

    ```groovy
    output {
        cowpy_art {
            path { meta, file -> meta.lang_group }
        }
    }
    ```

### Additional resources

- [map operator](https://www.nextflow.io/docs/latest/operator.html#map)
- [multiMap operator](https://www.nextflow.io/docs/latest/reference/operator.html#multimap)
- [stdout output qualifier](https://www.nextflow.io/docs/latest/process.html#outputs)

---

## What's next?

Return to the [menu of Side Quests](../) or click the button in the bottom right of the page to move on to the next topic in the list.
