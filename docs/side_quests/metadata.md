# Metadata in workflows

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
- Create and manipulate metadata maps
- Add new metadata fields during workflow execution
- Use metadata to customize process behavior

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

Let's move into the directory where the files for this tutorial are located.

```bash
cd side-quests/metadata
```

You can set VSCode to focus on this directory:

```bash
code .
```

#### Review the materials

You'll find a main workflow file and a `data` directory containing a datasheet and a handful of data files.

```console title="Directory contents"
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

The datasheet list the paths to the data files and some associated metadata, organized in 3 columns:

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

We will also provide you with a containerized language analysis tool called `langid`.

#### Review the assignment

Your challenge is to write a Nextflow workflow that will:

1. **Identify** the language in each file automatically
2. **Group** files by language family (Germanic vs Romance languages)
3. **Customize** the processing for each file based on its language and metadata
4. **Organize** outputs by language group

This represents a typical workflow pattern where file-specific metadata drives processing decisions; exactly the kind of problem that metadata maps solve elegantly.

#### Readiness checklist

Think you're ready to dive in?

- [ ] I understand the goal of this course and its prerequisites
- [ ] My codespace is up and running
- [ ] I've set my working directory appropriately
- [ ] I understand the assignment

If you can check all the boxes, you're good to go.

---

## 1. Load metadata from a datasheet

Open the `main.nf` workflow file to examine the workflow stub we're giving you as a starting point.

```groovy title="main.nf" linenums="1"
#!/usr/bin/env nextflow

workflow  {

    ch_datasheet = channel.fromPath("./data/datasheet.csv")

}
```

You can see we've set up a basic channel factory to load the example datasheet as a file, but that won't yet read in the contents of the file.
Let's start by adding that.

### 1.1. Read in contents with `splitCsv`

We need to choose an operator that will parse the file contents appropriately with minimal effort on our part.
Since our datasheet is in CSV format, this is a job for the [`splitCsv`](https://www.nextflow.io/docs/latest/reference/operator.html#splitcsv) operator, which loads each row in the file as an element in the channel.

Make the following changes to add a `splitCsv()` operation to the channel construction code, plus a `view()` operation to check that the contents of the file are getting loaded into the channel correctly.

=== "After"

    ```groovy title="main.nf" linenums="3" hl_lines="6-7"
    workflow  {

        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .view()

    }
    ```

=== "Before"

    ```groovy title="main.nf" linenums="3"
    workflow  {

        ch_datasheet = channel.fromPath("./data/datasheet.csv")

    }
    ```

Note that we're using the `header: true` option to tell Nextflow to read the first row of the CSV file as the header row.

Let's what comes out of that, shall we?
Run the workflow:

```bash
nextflow run main.nf
```

??? example "Output"

    ```console
    N E X T F L O W   ~  version 24.10.4

    Launching `main.nf` [exotic_albattani] DSL2 - revision: c0d03cec83

    [id:sampleA, character:squirrel, recording:/workspaces/training/side-quests/metadata/data/bonjour.txt]
    [id:sampleB, character:tux, recording:/workspaces/training/side-quests/metadata/data/guten_tag.txt]
    [id:sampleC, character:sheep, recording:/workspaces/training/side-quests/metadata/data/hallo.txt]
    [id:sampleD, character:turkey, recording:/workspaces/training/side-quests/metadata/data/hello.txt]
    [id:sampleE, character:stegosaurus, recording:/workspaces/training/side-quests/metadata/data/hola.txt]
    [id:sampleF, character:moose, recording:/workspaces/training/side-quests/metadata/data/salut.txt]
    [id:sampleG, character:turtle, recording:/workspaces/training/side-quests/metadata/data/ciao.txt]
    ```

We can see that operator has constructed a map of key-value pairs for each row in the CSV file, with the column headers as keys for the corresponding values.

Each map entry corresponds to a column in our datasheet:

- `id`
- `character`
- `recording`

This is great! It makes it easy to access specific fields from each file.
For example, we could access the file ID with `id` or the txt file path with `recording`.

### 1.2. Pick out specific fields with `map`

Let's say we want to access the `character` column from the datasheet and print it.
We can use the Nextflow `map` operator to iterate over each item in our channel and specifically pick out the `character` entry from the map object.

Make the following edits to the workflow:

=== "After"

    ```groovy title="main.nf" linenums="3" hl_lines="5-7"
    workflow  {

        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map{ row ->
                row.character
            }
            .view()

    }
    ```

=== "Before"

    ```groovy title="main.nf" linenums="3"
    workflow  {

        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .view()

    }
    ```

Now run the workflow again:

```bash
nextflow run main.nf
```

??? example "Output"

    ```console
    N E X T F L O W   ~  version 24.10.4

    Launching `main.nf` [exotic_albattani] DSL2 - revision: c0d03cec83

    squirrel
    tux
    sheep
    turkey
    stegosaurus
    moose
    turtle
    ```

Success! We've taken advantage of the map structure derived from our datasheet to access the values from individual columns for each row.

<details>
  <summary>More about maps</summary>

In Groovy, the programming language that Nextflow is built on, a map is a key-value data structure similar to dictionaries in Python, objects in JavaScript, or hashes in Ruby.

For example:

```groovy title="Groovy map"
def my_map = [id:'sampleA', character:'squirrel']
println my_map.id  // Prints: sampleA
```

And here's a runnable script that applies this in practice:

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
Nextflow 25.10.0 is available - Please consider updating your version to it

N E X T F L O W   ~  version 25.04.3

Launching `map_demo.nf` [cheesy_plateau] DSL2 - revision: fae5b8496e

map: [id:sampleA, character:squirrel]
id: sampleA
character: squirrel
```

</details>

Now that we've successfully read in the datasheet and have access to the data in each row, we can begin implementing our pipeline logic.

### 1.3. Organize the metadata into a 'meta map'

In the current state of the workflow, the input files (under the `recording` key) and associated metadata (`id`, `character`) are all on the same footing, like they're all in one big bag.
The practical consequence is that every process that consumes this channel would need to be configured with this structure in mind:

```groovy title="Syntax example"
    input:
    tuple val(id), val(character), file(recording)
```

That's fine as long as the number of columns in the datasheet doesn't change.
However, if you add even just one column to the datasheet, the shape of the channel will no longer match what the process expects, and the workflow will produce errors.
It also makes the process hard to share with others who might have slightly different input data, and you might end up having to hard-code variables into the process that aren't needed by the script block.

To avoid this problem, we need to find a way of keeping the channel structure consistent irrespective of how many columns that datasheet contains.

We can do that by collecting all the metadata into a item within the tuple, which we'll call the metadata map, or more simply 'meta map'.

Make the following edits to the `map` operation:

=== "After"

    ```groovy title="main.nf" linenums="5" hl_lines="8"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                    [ [id: row.id, character: row.character], row.recording ]
            }
            .view()
    ```

=== "Before"

    ```groovy title="main.nf" linenums="5" hl_lines="8"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map{ row ->
                row.character
            }
            .view()
    ```

We've restructured our channel elements into a tuple consisting of two elements, the meta map and the corresponding file object.

Let's run the workflow:

```bash title="View meta map"
nextflow run main.nf
```

??? example "Output"

    ```console title="View meta map"
    N E X T F L O W   ~  version 24.10.4

    Launching `main.nf` [lethal_booth] DSL2 - revision: 0d8f844c07

    [[id:sampleA, character:squirrel], /workspaces/training/side-quests/metadata/data/bonjour.txt]
    [[id:sampleB, character:tux], /workspaces/training/side-quests/metadata/data/guten_tag.txt]
    [[id:sampleC, character:sheep], /workspaces/training/side-quests/metadata/data/hallo.txt]
    [[id:sampleD, character:turkey], /workspaces/training/side-quests/metadata/data/hello.txt]
    [[id:sampleE, character:stegosaurus], /workspaces/training/side-quests/metadata/data/hola.txt]
    [[id:sampleF, character:moose], /workspaces/training/side-quests/metadata/data/salut.txt]
    [[id:sampleG, character:turtle], /workspaces/training/side-quests/metadata/data/ciao.txt]
    ```

Now, each element in the channel contains the metadata map first (_e.g._ `[id:sampleA, character:squirrel]`) and the corresponding file object second (_e.g._ `/workspaces/training/side-quests/metadata/data/bonjour.txt`).

As a result, adding more columns in the datasheet will make more metadata available in the `meta` map, but won't change the channel shape.
This enables us to write processes that consume the channel without having to hard-code the metadata items into the input specification:

```groovy title="Syntax example"
    input:
    tuple val(meta), file(recording)
```

This is a widely used convention for organizing metadata in Nextflow workflows.

### Takeaway

In this section, you've learned:

- **Why metadata is important:** Keeping metadata with your data preserves important file information throughout the workflow.
- **How to read in a datasheets:** Using `splitCsv` to read CSV files with header information and transform rows into structured data
- **How to create a meta map:** Separating metadata from file data using the tuple structure `[ [id:value, ...], file ]`

---

## 2. Manipulating metadata

Now that we have our metadata loaded, let's do something with it!

We're going to use a tool called [`langid`](https://github.com/saffsd/langid.py) to identify the language contained in each creature's recording file.
The tool comes pre-trained on a set of languages, and given a snippet of text, it will output a language prediction and an associated probability score, both to `stdout`.

### 2.1. Import the process and examine the code

We provide you with a pre-written process module called `IDENTIFY_LANGUAGE` that wraps the `langid` tool, so you just need to add an include statement before the workflow block.

Make the following edit to the workflow:

=== "After"

    ```groovy title="main.nf" linenums="1" hl_lines="3"
    #!/usr/bin/env nextflow

    include { IDENTIFY_LANGUAGE } from './modules/langid.nf'

    workflow {
    ```

=== "Before"

    ```groovy title="main.nf" linenums="1"
    #!/usr/bin/env nextflow

    workflow {
    ```

You can open the module file to examine its code:

```groovy title="modules/langid.nf" linenums="1" hl_lines="11 14"
#!/usr/bin/env nextflow

/*
* Use langid to predict the language of each input file
*/
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

As you can see, the input definition uses the same `tuple val(meta), path(file)` structure that we just applied to our input channel.

Note that in the `script` section, we use `sed` to remove the probability score, clean up the string by removing newline characters, and return only the language prediction.
Since the output is printed directly to the console, we use Nextflow’s [`stdout` output qualifier](https://www.nextflow.io/docs/latest/process.html#outputs) to capture the string as output.

The output definition is structured as a tuple with a similar structure to the input's, except it also contains `stdout` as a third element.
This `tuple val(meta), path(file), <output>` pattern keeps the metadata associated with both the input data and the outputs as it flows through the pipeline.

### 2.2. Identify the language of each greeting

We still need to add a call to the `IDENTIFY_LANGUAGE` process, so let's do that now.
Make the following edits to the workflow:

=== "After"

    ```groovy title="main.nf" linenums="8" hl_lines="7-9"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [[id: row.id, character: row.character], row.recording]
            }

        // Run langid to identify the language of each greeting
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out.view()
    ```

=== "Before"

    ```groovy title="main.nf" linenums="8" hl_lines="6"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [[id: row.id, character: row.character], row.recording]
            }
            .view()
    ```

And now run the workflow:

```bash
nextflow run main.nf
```

??? example "Output"

    ```console
    N E X T F L O W   ~  version 24.10.4

    Launching `main.nf` [voluminous_mcnulty] DSL2 - revision: f9bcfebabb

    executor >  local (7)
    [4e/f722fe] IDENTIFY_LANGUAGE (7) [100%] 7 of 7 ✔
    [[id:sampleA, character:squirrel], /workspaces/training/side-quests/metadata/work/eb/f7148ebdd898fbe1136bec6a714acb/bonjour.txt, fr]
    [[id:sampleB, character:tux], /workspaces/training/side-quests/metadata/work/16/71d72410952c22cd0086d9bca03680/guten_tag.txt, de]
    [[id:sampleD, character:turkey], /workspaces/training/side-quests/metadata/work/c4/b7562adddc1cc0b7d414ec45d436eb/hello.txt, en]
    [[id:sampleC, character:sheep], /workspaces/training/side-quests/metadata/work/ea/04f5d979429e4455e14b9242fb3b45/hallo.txt, de]
    [[id:sampleF, character:moose], /workspaces/training/side-quests/metadata/work/5a/6c2b84bf8fadb98e28e216426be079/salut.txt, fr]
    [[id:sampleE, character:stegosaurus], /workspaces/training/side-quests/metadata/work/af/ee7c69bcab891c40d0529305f6b9e7/hola.txt, es]
    [[id:sampleG, character:turtle], /workspaces/training/side-quests/metadata/work/4e/f722fe47271ba7ebcd69afa42964ca/ciao.txt, it]
    ```

Excellent! For each of our creatures, we now have a prediction for what language the file contains.

As noted earlier, we've also included the input file and the meta map in the output, which means that both remain associated with the new information we've just produced.
This will prove useful in the next step.

!!! note

    More generally, this pattern of keeping the meta map associated with results makes it easier to associate related results that share the same identifiers.

    As you will have learned already, you can't rely on the order of items in channels to match results across them.
    Instead, you must use keys to associate data correctly, and meta maps provide an ideal structure for this purpose.

    We explore this use case in detail in the [Splitting & Grouping](./splitting_and_grouping.md) side quest.

### 2.3. Augment metadata with process outputs

Given that the results we just produced are in themselves a form of metadata about the contents of the files, we're going to want to add them to our meta map.

However, we don't want to modify the existing meta map in place.
From a technical standpoint, it is _possible_ to do that, but it's unsafe.
So instead, we're going to create a new meta map containing the contents of the existing meta map plus a new key-value pair holding the new information.

We do that using the [`map`](https://www.nextflow.io/docs/latest/operator.html#map) operator.

<!-- TODO Explain first how the map merge stuff works

```groovy title="Syntax example"
.map { meta, new_info ->
    [meta + [new_key: new_info]]
}
```

-->

In the case of our workflow, we also need to account for the presence of the `file` object in the tuple.

Make the following edits to the workflow:

=== "After"

    ```groovy title="main.nf" linenums="28" hl_lines="2-6"
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out
            .map { meta, file, lang ->
                [meta + [lang: lang], file]
            }
            .view()
    }
    ```

=== "Before"

    ```groovy title="main.nf" linenums="28" hl_lines="2"
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out.view()
    ```

The `map` operator takes each channel element and processes it to create a modified version.
Inside the closure `{ meta, file, lang -> ... }`, we then take the existing `meta` map, create a new map `[lang:lang]`, and merge both together using `+`.

The `+` operator in Groovy merges two maps together.
So if our original `meta` was `[id:sampleA, character:squirrel]`, then `meta + [lang:'fr']` creates a new map: `[id:sampleA, character:squirrel, lang:fr]`.

!!! Note

    The `+` notation with maps creates an entirely new map object, which is what we want.
    If we'd done something like `meta.lang = lang` we'd have been modifying the original object, which can lead to unpredictable effects.

```bash title="View new meta map key"
nextflow run main.nf -resume
```

```console title="View new meta map key"

 N E X T F L O W   ~  version 24.10.4

Launching `main.nf` [cheeky_fermat] DSL2 - revision: d096281ee4

[4e/f722fe] IDENTIFY_LANGUAGE (7) [100%] 7 of 7, cached: 7 ✔
[[id:sampleA, character:squirrel, lang:fr], /workspaces/training/side-quests/metadata/work/eb/f7148ebdd898fbe1136bec6a714acb/bonjour.txt]
[[id:sampleB, character:tux, lang:de], /workspaces/training/side-quests/metadata/work/16/71d72410952c22cd0086d9bca03680/guten_tag.txt]
[[id:sampleC, character:sheep, lang:de], /workspaces/training/side-quests/metadata/work/ea/04f5d979429e4455e14b9242fb3b45/hallo.txt]
[[id:sampleD, character:turkey, lang:en], /workspaces/training/side-quests/metadata/work/c4/b7562adddc1cc0b7d414ec45d436eb/hello.txt]
[[id:sampleF, character:moose, lang:fr], /workspaces/training/side-quests/metadata/work/5a/6c2b84bf8fadb98e28e216426be079/salut.txt]
[[id:sampleE, character:stegosaurus, lang:es], /workspaces/training/side-quests/metadata/work/af/ee7c69bcab891c40d0529305f6b9e7/hola.txt]
[[id:sampleG, character:turtle, lang:it], /workspaces/training/side-quests/metadata/work/4e/f722fe47271ba7ebcd69afa42964ca/ciao.txt]
```

Nice, we expanded our meta map with new information we gathered in the pipeline.
After running the language prediction, each element in the output channel looks like this:

```console
[meta, file, lang]  // e.g. [[id:sampleA, character:squirrel], bonjour.txt, fr]
```

<!-- TODO Should we also show how to remove a key using subMap?! -->

### 2.4. Assign a language group using conditionals

Now that we have our language predictions, let's use the information to assign them into new groups. In our example data, we have provided data sets that belong either to `germanic` (either English or German) or `romance` (French, Spanish, Italian) languages.

We can add another `map` operator to assign either group (add this below your last map operation).

=== "After"

    ```groovy title="main.nf" linenums="31" hl_lines="4-15"
            .map { meta, file, lang ->
                [meta + [lang: lang], file]
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
            .view()
    ```

=== "Before"

    ```groovy title="main.nf" linenums="31"
            .map { meta, file, lang ->
                [meta + [lang: lang], file]
            }
            .view()
    ```

Let's rerun it

```bash title="View language groups"
nextflow run main.nf -resume
```

```console title="View language groups"
 N E X T F L O W   ~  version 24.10.4

Launching `main.nf` [wise_almeida] DSL2 - revision: 46778c3cd0

[da/652cc6] IDENTIFY_LANGUAGE (7) [100%] 7 of 7, cached: 7 ✔
[[id:sampleA, character:squirrel, lang:fr, lang_group:romance], /workspaces/training/side-quests/metadata/data/bonjour.txt]
[[id:sampleB, character:tux, lang:de, lang_group:germanic], /workspaces/training/side-quests/metadata/data/guten_tag.txt]
[[id:sampleC, character:sheep, lang:de, lang_group:germanic], /workspaces/training/side-quests/metadata/data/hallo.txt]
[[id:sampleD, character:turkey, lang:en, lang_group:germanic], /workspaces/training/side-quests/metadata/data/hello.txt]
[[id:sampleE, character:stegosaurus, lang:es, lang_group:romance], /workspaces/training/side-quests/metadata/data/hola.txt]
[[id:sampleF, character:moose, lang:fr, lang_group:romance], /workspaces/training/side-quests/metadata/data/salut.txt]
[[id:sampleG, character:turtle, lang:it, lang_group:romance], /workspaces/training/side-quests/metadata/data/ciao.txt]
```

Let's understand how this transformation works.
The `map` operator here again takes a closure that processes each element in the channel.
Inside the closure, we're using an if-clause to create a new language group classification.

Here's what's happening step by step:

- Create a new field `lang_group`: We set a default to `unknown`
- Extract existing metadata: We access `meta.lang` (the language we predicted earlier) from the existing meta map
- Apply conditional logic: We use an if-clause to determine the language group based on the language: Is `meta.lang` either `de` or `en`, we re-assign `lang_group` to `germanic`, if `fr`, `es`, or `it`, then we re-assign to `romance`
- Merge with existing metadata: We use `meta + [lang_group:lang_group]` in the same way as before to combine the existing meta map with our new field

The resulting channel elements maintain their `[meta, file]` structure, but the meta map now includes this new classification.

### Takeaway

In this section, you've learned how to :

- **Apply input metadata to output channels**: Copying metadata in this way allows us to associate results later on based on metadata content.
- **Create custom keys**: You created two new keys in your meta map, merging them with `meta + [new_key:value]` into the existing meta map.One based on a computed value from a process, and one based on a condition you set in the `map` operator.

These allow you to associate new and existing metadata with files as you progress through your pipeline.

---

## 3. Use meta map information in a process

Let's make some fun characters say the phrases from the files our channel.
In the [hello-nextflow training](../hello_nextflow/05_hello_containers.md), you already encountered the `cowpy` package, a python implementation of a tool called `cowsay` that generates ASCII art to display arbitrary text inputs in a fun way. We will re-use a process from there.

Copy in the process before your workflow block:

=== "After"

    ```groovy title="main.nf" linenums="20"
    /*
     * Generate ASCII art with cowpy
    */
    process COWPY {

        publishDir "results/", mode: 'copy'

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

        input:
        tuple val(meta), path(input_file)

        output:
        tuple val(meta), path("cowpy-${input_file}")

        script:
        """
        cat $input_file | cowpy -c "stegosaurus" > cowpy-${input_file}
        """

    }

    workflow{}
    ```

=== "Before"

    ```groovy title="main.nf" linenums="20"
    workflow{
    ```

### 3.1. Use meta map information in the process definition

Let's run our files through `COWPY` and remove our `view` statement:

=== "After"

    ```groovy title="main.nf" linenums="59" hl_lines="3"
                [ meta + [lang_group: lang_group], file ]
            }
        COWPY(ch_languages)
    ```

=== "Before"

    ```groovy title="main.nf" linenums="59" hl_lines="3"
                [meta + [lang_group: lang_group], file]
            }
            .view()
    ```

The process definition as provided would direct results to the `results` folder, but let's make a tweak to be a little smarter.
Given we have been trying to figure out what languages our samples were in, let's group the samples by language in the output directory.

Earlier, we added the predicted language to the `meta` map.
We can access this `key` in the process and use it in the `publishDir` directive:

=== "After"

    ```groovy title="main.nf" linenums="23" hl_lines="3"
    process COWPY {

        publishDir "results/${meta.lang_group}", mode: 'copy'

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
    ```

=== "Before"

    ```groovy title="main.nf" linenums="23" hl_lines="3"
    process COWPY {

        publishDir "results/", mode: 'copy'

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
    ```

Let's run this:

```bash title="Use cowpy"
nextflow run main.nf -resume
```

You should now see a new folder called `results`:

```console title="Results folder"
results/
├── germanic
│   ├── cowpy-guten_tag.txt
│   ├── cowpy-hallo.txt
│   └── cowpy-hello.txt
└── romance
    ├── cowpy-bonjour.txt
    ├── cowpy-ciao.txt
    ├── cowpy-hola.txt
    └── cowpy-salut.txt
```

Success! All our phrases are correctly sorted and we now see which of them correspond to which language.

Let's take a look at `cowpy-salut.txt`:

```console title="cowpy-salut.txt"
 ____________________
/ Salut, ça va?      \
\ à plus dans le bus /
 --------------------
\                             .       .
 \                           / `.   .' "
  \                  .---.  <    > <    >  .---.
   \                 |    \  \ - ~ ~ - /  /    |
         _____          ..-~             ~-..-~
        |     |   \~~~\.'                    `./~~~/
       ---------   \__/                        \__/
      .'  O    \     /               /       \  "
     (_____,    `._.'               |         }  \/~~~/
      `----.          /       }     |        /    \__/
            `-.      |       /      |       /      `. ,~~|
                ~-.__|      /_ - ~ ^|      /- _      `..-'
                     |     /        |     /     ~-.     `-. _  _  _
                     |_____|        |_____|         ~ - . _ _ _ _ _>
```

Look through the other files.
All phrases should be spoken by the fashionable stegosaurus.

How did this work? The `publishDir` directive is evaluated at runtime when the process executes.
Each process task gets its own meta map from the input tuple When the directive is evaluated, `${meta.lang_group}` is replaced with the actual group language value for that dataset creating the dynamic paths like `results/romance`.

### 3.2. Customize the character

In our datasheet, we have another column: `character`.
To tailor the tool parameters per file, we can also access information from the `meta` map in the script section.
This is really useful in cases where a tool should have different parameters for each file.

Let's customize the characters by changing the `cowpy` command:

=== "After"

    ```groovy title="main.nf" linenums="37" hl_lines="3"
        script:
        """
        cat $input_file | cowpy -c ${meta.character} > cowpy-${input_file}
        """
    ```

=== "Before"

    ```groovy title="main.nf" linenums="37"
        script:
        """
        cat $input_file | cowpy -c "stegosaurus" > cowpy-${input_file}
        """
    ```

Let's run this:

```bash title="Use cowpy"
nextflow run main.nf -resume
```

and take another look at our french phrase:

```console title="romance/cowpy-salut.txt"
 ____________________
/ Salut, ça va?      \
\ à plus dans le bus /
 --------------------
  \
   \   \_\_    _/_/
    \      \__/
           (oo)\_______
           (__)\       )\/\
               ||----w |
               ||     ||
```

This approach differs from using pipeline parameters (`params`), which generally apply the same configuration to all files in your workflow.
By leveraging metadata applied to each item in a channel, you can fine-tune process behavior on a per-file basis.

#### 3.2.1. Exploiting metadata at the workflow level

In the example above, by using a property of the meta map in the script block, we introduce a hard requirement on the properties that must be present.
Anyone running with a sample sheet that did not contain the `character` property would encounter an error.
The process `input:` only says that the `meta` map is required, so someone trying to use this process in another workflow might not notice immediately that the `character` property was required.

A better approach is to make the required metadata property an explicit input rather than accessing it from within the meta map. This makes the requirement clear and provides better error messages.
Here's how to refactor the COWPY process:

=== "After"

    ```groovy title="main.nf" linenums="30" hl_lines="2"
        input:
        tuple val(meta), val(character), path(input_file)
    ```

=== "Before"

    ```groovy title="main.nf" linenums="30" hl_lines="2"
        input:
        tuple val(meta), path(input_file)
    ```

Then, at the workflow level, we extract the `character` property from the metadata and pass it as a separate input:

=== "After"

    ```groovy title="main.nf" linenums="61" hl_lines="1"
        COWPY(ch_languages.map{meta, file -> [meta, meta.character, file]})
    ```

=== "Before"

    ```groovy title="main.nf" linenums="61" hl_lines="1"
        COWPY(ch_languages)
    ```

**Why this approach is better:**

1. **Clear Requirements**: The process input explicitly shows that `character` is required
2. **Better Error Messages**: If `character` is missing, Nextflow will fail at the input stage with a clear error
3. **Self-Documenting**: Anyone using this process can immediately see what inputs are needed
4. **Reusable**: The process is easier to redeploy in other contexts

This highlights an important design principle:
Use the meta map for optional, descriptive information, but extract required values as explicit inputs.
The meta map is excellent for keeping channel structures clean and preventing arbitrary channel structures, but for mandatory elements that are directly referenced in a process, extracting them as explicit inputs creates more robust and maintainable code.

### Takeaway

In this section, you've learned how to:

- **Tweak directives using meta values**: Using meta map values in `publishDir` directives to create dynamic output paths based on the file's metadata

- **Tweak the script section based on meta values**: Customizing tool parameters per file using meta information in the `script` section

---

## Summary

In this side quest, you've explored how to effectively work with metadata in Nextflow workflows.

This pattern of keeping metadata explicit and attached to the data is a core best practice in Nextflow, offering several advantages over hardcoding file information:

- File metadata stays associated with files throughout the workflow
- Process behavior can be customized per file
- Output organization can reflect file metadata
- File information can be expanded during pipeline execution

Applying this pattern in your own work will enable you to build robust, maintainable bioinformatics workflows.

### Key patterns

1.  **Reading and Structuring Metadata:** Reading CSV files and creating organized metadata maps that stay associated with your data files.

    ```groovy
    channel.fromPath('datasheet.csv')
      .splitCsv(header: true)
      .map { row ->
          [ [id:row.id, character:row.character], row.recording ]
      }
    ```

2.  **Expanding Metadata During Workflow** Adding new information to your metadata as your pipeline progresses by adding process outputs and deriving values through conditional logic.

    - Adding new keys based on process output

    ```groovy
    .map { meta, file, lang ->
      [ meta + [lang:lang], file ]
    }
    ```

    - Adding new keys using a conditional clause

    ```groovy
    .map{ meta, file ->
        if ( meta.lang.equals("de") || meta.lang.equals('en') ){
            lang_group = "germanic"
        } else if ( meta.lang in ["fr", "es", "it"] ) {
            lang_group = "romance"
        } else {
            lang_group = "unknown"
        }
    }
    ```

3.  **Customizing Process Behavior:** Using metadata to adapt how processes handle different files.

    - Using meta values in Process Directives

    ```groovy
    publishDir "results/${meta.lang_group}", mode: 'copy'
    ```

    - Adapting tool parameters for individual files

    ```groovy
    cat $input_file | cowpy -c ${meta.character} > cowpy-${input_file}
    ```

### Additional resources

- [map](https://www.nextflow.io/docs/latest/operator.html#map)
- [stdout](https://www.nextflow.io/docs/latest/process.html#outputs)

---

## What's next?

Return to the [menu of Side Quests](./index.md) or click the button in the bottom right of the page to move on to the next topic in the list.
