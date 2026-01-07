# Metadata and meta maps

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

    ```groovy title="main.nf" linenums="3" hl_lines="4-5"
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

??? success "Command output"

    ```console
     N E X T F L O W   ~  version 25.10.2

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

??? info "(Optional) More about maps"

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
     N E X T F L O W   ~  version 25.10.2

    Launching `map_demo.nf` [cheesy_plateau] DSL2 - revision: fae5b8496e

    map: [id:sampleA, character:squirrel]
    id: sampleA
    character: squirrel
    ```

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

??? success "Command output"

    ```console
     N E X T F L O W   ~  version 25.10.2

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

Now that we've successfully read in the datasheet and have access to the data in each row, we can begin implementing our pipeline logic.

### 1.3. Organize the metadata into a 'meta map'

In the current state of the workflow, the input files (under the `recording` key) and associated metadata (`id`, `character`) are all on the same footing, like they're all in one big bag.
The practical consequence is that every process that consumes this channel would need to be configured with this structure in mind:

```groovy
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

    ```groovy title="main.nf" linenums="5" hl_lines="4"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                    [ [id: row.id, character: row.character], row.recording ]
            }
            .view()
    ```

=== "Before"

    ```groovy title="main.nf" linenums="5" hl_lines="4"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map{ row ->
                row.character
            }
            .view()
    ```

We've restructured our channel elements into a tuple consisting of two elements, the meta map and the corresponding file object.

Let's run the workflow:

```bash
nextflow run main.nf
```

??? success "Command output"

    ```console title="View meta map"
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [lethal_booth] DSL2 - revision: 0d8f844c07

    [[id:sampleA, character:squirrel], /workspaces/training/side-quests/metadata/data/bonjour.txt]
    [[id:sampleB, character:tux], /workspaces/training/side-quests/metadata/data/guten_tag.txt]
    [[id:sampleC, character:sheep], /workspaces/training/side-quests/metadata/data/hallo.txt]
    [[id:sampleD, character:turkey], /workspaces/training/side-quests/metadata/data/hello.txt]
    [[id:sampleE, character:stegosaurus], /workspaces/training/side-quests/metadata/data/hola.txt]
    [[id:sampleF, character:moose], /workspaces/training/side-quests/metadata/data/salut.txt]
    [[id:sampleG, character:turtle], /workspaces/training/side-quests/metadata/data/ciao.txt]
    ```

Now, each element in the channel contains the metadata map first and the corresponding file object second:

```console title="Example output structure"
[
  [id:sampleA, character:squirrel],
  /workspaces/training/side-quests/metadata/data/bonjour.txt
]
```

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

```groovy title="modules/langid.nf" linenums="1" hl_lines="9 12"
#!/usr/bin/env nextflow

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

As you can see, the input definition uses the same `tuple val(meta), path(file)` structure that we just applied to our input channel.

The output definition is structured as a tuple with a similar structure to the input's, except it also contains `stdout` as a third element.
This `tuple val(meta), path(file), <output>` pattern keeps the metadata associated with both the input data and the outputs as it flows through the pipeline.

Note that we're using Nextflow’s [`stdout`](https://www.nextflow.io/docs/latest/process.html#outputs) output qualifier here because the tool prints its output directly to the console rather than writing a file; and we use `sed` in the command line to remove the probability score, clean up the string by removing newline characters, and return only the language prediction.

### 2.2. Add a call to `IDENTIFY_LANGUAGE`

Now that the process is available to the workflow, we can add a call to the `IDENTIFY_LANGUAGE` process to run it on the data channel.

Make the following edits to the workflow:

=== "After"

    ```groovy title="main.nf" linenums="7" hl_lines="7-9"
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

    ```groovy title="main.nf" linenums="7" hl_lines="6"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [[id: row.id, character: row.character], row.recording]
            }
            .view()
    ```

Note that we've removed the original `.view()` operation in the channel construction.

We can now run the workflow:

```bash
nextflow run main.nf
```

??? success "Command output"

    ```console
     N E X T F L O W   ~  version 25.10.2

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

Excellent! We now have a prediction for what language each character speaks.

And as noted earlier, we've also included the input file and the meta map in the output, which means that both remain associated with the new information we've just produced.
This will prove useful in the next step.

!!! note

    More generally, this pattern of keeping the meta map associated with results makes it easier to associate related results that share the same identifiers.

    As you will have learned already, you can't rely on the order of items in channels to match results across them.
    Instead, you must use keys to associate data correctly, and meta maps provide an ideal structure for this purpose.

    We explore this use case in detail in the [Splitting & Grouping](./splitting_and_grouping.md) side quest.

### 2.3. Augment metadata with process outputs

Given that the results we just produced are in themselves a form of metadata about the contents of the files, it would be useful to add them to our meta map.

However, we don't want to modify the existing meta map in place.
From a technical standpoint, it is _possible_ to do that, but it's unsafe.

So instead, we'll create a new meta map containing the contents of the existing meta map plus a new `lang: lang_id` key-value pair holding the new information, using the `+` operator (a Groovy feature).
And we'll combine this with a [`map`](https://www.nextflow.io/docs/latest/operator.html#map) operation to replace the old map with the new one.

Here are the edits you need to make to the workflow:

=== "After"

    ```groovy title="main.nf" linenums="13" hl_lines="3-7"
        // Run langid to identify the language of each greeting
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out
            .map { meta, file, lang_id ->
                [meta + [lang: lang_id], file]
            }
            .view()
    ```

=== "Before"

    ```groovy title="main.nf" linenums="13" hl_lines="3"
        // Run langid to identify the language of each greeting
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out.view()
    ```

If you're not yet familiar with the `+` operator, or if this seems confusing, take a few minutes to go through the detailed explanation below.

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
    new_map = [map1 + [lang: lang_id]]
    ```

    Here, `[lang: new_info]` creates a new unnamed map on the fly, and `map1 + ` merges `map1` with the new unnamed map, producing the same `new_map` contents as before.

    Neat, right?

    **Now let's transpose that into the context of a Nextflow `channel.map()` operation.**

    The code becomes:

    ```groovy
    .map { map1, lang_id ->
        [map1 + [lang: lang_id]]
    }
    ```

    This does the following:

    - `map1, lang_id ->` takes the two items in the tuple
    - `[map1 + [lang: lang_id]]` creates the new map as detailed above

    The output is a single unnamed map with the same contents as `new_map` in our example above.
    So we've effectively transformed:

    ```groovy
    [id: 'sampleA', character: 'squirrel'], 'it'
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

Once you're confident you understand what this code is doing, run the workflow to see if it worked:

```bash
nextflow run main.nf -resume
```

??? success "Command output"

    ```console
     N E X T F L O W   ~  version 25.10.2

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

Yep, that checks out!
We've neatly reorganized the output of the process from `meta, file, lang_id` so that `lang_id` is now one of the keys in the meta map, and the channel's tuples fit the `meta, file` model once again.

<!-- TODO (future) Should we also show how to remove a key using subMap?! Or note where to find that. -->

### 2.4. Assign a language group using conditionals

Now that we have our language predictions, let's use the information to assign some new groupings.

In our example data, the languages used by our characters can be grouped into germanic languages (English, German) and romance language (French, Spanish, Italian).
It might be useful to have that classification readily available somewhere later in the pipeline, so let's add that information in the meta map.

And, good news, this is yet another case that lends itself perfectly to using the `map` operator!

Specifically, we're going to define a variable called `lang_group`, use some simple conditional logic to determine what value to assign to the `lang_group` for each piece of data.

The general syntax is going to look like this:

```groovy
.map { meta, file ->

    // conditional logic defining lang_group goes here

    [meta + [lang_group: lang_group], file]
}
```

You can see this is very similar to the on-the-fly map merging operation we used in the previous step.
We just need to write the conditional statements.

Here's the conditional logic we want to apply:

- Define a variable called `lang_group` with default value `'unknown'`.
- If `lang` is either German (`'de'`) or English (`'en'`), change `lang_group` to `germanic`.
- Else if `lang` is included in a list containing French (`'fr'`), Spanish (`'es'`) and Italian (`'it'`), change `lang_group` to `romance`.

Try taking a stab at writing it yourself if you already know how to write conditional statements in Nextflow.

!!! tip

    You can access the value of `lang` within the map operation with `meta.lang`.

You should end up making the following changes to the workflow:

=== "After"

    ```groovy title="main.nf" linenums="13" hl_lines="7-19"
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
            .view()
    ```

=== "Before"

    ```groovy title="main.nf" linenums="13" hl_lines="7"
        // Run langid to identify the language of each greeting
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out
            .map { meta, file, lang_id ->
                [meta + [lang: lang_id], file]
            }
            .view()
    ```

Here are the key points:

- We use `def lang_group = "unknown"` to create the `lang_group` variable with default value set to `unknown`.
- We use an `if {} else if {}` structure for the conditional logic, with alternative `.equals()` tests for the two germanic languages, and a test for existence in a list for the three romance languages.
- We use the `meta + [lang_group:lang_group]` merge operation as previously to generate the updated meta map.

<!-- TODO (future) Add note/links to relevant docs in additional resources section -->

Once that all makes sense, run the workflow again to see the result:

```bash
nextflow run main.nf -resume
```

??? success "Command output"

    ```console
     N E X T F L O W   ~  version 25.10.2

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

As you can see, the channel elements maintain their `[meta, file]` structure, but the meta map now includes this new classification.

### Takeaway

In this section, you've learned how to :

- **Apply input metadata to output channels**: Copying metadata in this way allows us to associate results later on based on metadata content.
- **Create custom keys**: You created two new keys in your meta map, merging them with `meta + [new_key:value]` into the existing meta map.One based on a computed value from a process, and one based on a condition you set in the `map` operator.

These allow you to associate new and existing metadata with files as you progress through your pipeline.
Even if you're not using metadata as part of a process, keeping the meta map associated with the data like this makes it easy to keep all the relevant information together.

---

## 3. Using meta map information in a process

Now that you know how to create and update the meta map, we can get to the really fun bit: actually using the metadata in a process.

More specifically, we're going to add a second step to our workflow to draw each animal as ASCII art and make it say the recorded text in a speech bubble.
We're going to do this using a tool called [`cowpy`](https://github.com/jeffbuttars/cowpy).

??? info "What does `cowpy` do?"

    `cowpy` is a command-line tool that generates ASCII art to display arbitrary text inputs in a fun way.
    It is a python implementation of the classic [cowsay](https://en.wikipedia.org/wiki/Cowsay) tool by Tony Monroe.

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

If you worked through the Hello Nextflow course, you've already seen this tool in action.
If not, don't worry; we'll cover everything you need to know as we go.

### 3.1. Import the process and examine the code

We provide you with a pre-written process module called `COWPY` that wraps the `cowpy` tool, so you just need to add an include statement before the workflow block.

Make the following edit to the workflow:

=== "After"

    ```groovy title="main.nf" linenums="1" hl_lines="4"
    #!/usr/bin/env nextflow

    include { IDENTIFY_LANGUAGE } from './modules/langid.nf'
    include { COWPY } from './modules/cowpy.nf'

    workflow {
    ```

=== "Before"

    ```groovy title="main.nf" linenums="1"
    #!/usr/bin/env nextflow

    include { IDENTIFY_LANGUAGE } from './modules/langid.nf'

    workflow {
    ```

You can open the module file to examine its code:

```groovy title="modules/cowpy.nf" linenums="1"
#!/usr/bin/env nextflow

// Generate ASCII art with cowpy
process COWPY {

    publishDir "results/", mode: 'copy'

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

As you can see, this process is currently designed to take an input file (containing the text to be displayed) and a value specifying the character that should be drawn in ASCII art, usually provided at the workflow level by a command-line parameter.

### 3.2. Pass a meta map field as an input

When we used the `cowpy` tool in the Hello Nextflow course, we used a command-line parameter to determine what character to use to draw the final image.
That made sense, because we were only generating one image per run of the pipeline.

However, in this tutorial, we want to generate an appropriate image for each subject that we're processing, so using a command-line parameter would be too limiting.

Good news: we have a `character` column in our datasheet and therefore, in our meta map.
Let's use that to set the character that the process should use for each entry.

To that end, we'll need to do three things:

1. Give a name to the output channel coming out of the previous process so we can operate on it more conveniently.
2. Determine how to access the information of interest
3. Add a call to the second process and feed in the information appropriately.

Let's get started.

#### 3.2.1. Name the previous output channel

We applied the previous manipulations directly on the output channel of the first process, `IDENTIFY_LANGUAGE.out`.
In order to feed the contents of that channel to the next process (and do so in a way that is clear and easy to read) we want to give it its own name, `ch_languages`.

We can do that using the [`set`](https://www.nextflow.io/docs/latest/reference/operator.html#set) operator.

In the main workflow, replace the `.view()` operator with `.set { ch_languages }`, and add a line testing that we can refer to the channel by name.

=== "After"

    ```groovy title="main.nf" linenums="14" hl_lines="19 21 22"
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

        // Temporary: peek into ch_languages
        ch_languages.view()
    ```

=== "Before"

    ```groovy title="main.nf" linenums="14" hl_lines="19"
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
            .view()
    ```

Let's run this:

```bash
nextflow run main.nf
```

??? success "Command output"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `./main.nf` [friendly_austin] DSL2 - revision: 3dbe460fd6

    [36/cca6a7] IDENTIFY_LANGUAGE (7) | 7 of 7 ✔
    [[id:sampleB, character:tux, lang:de, lang_group:germanic], /workspaces/training/side-quests/metadata/work/e2/6db2402d83cf72081bcd2d11784714/guten_tag.txt]
    [[id:sampleA, character:squirrel, lang:fr, lang_group:romance], /workspaces/training/side-quests/metadata/work/6c/114c818317d169457d6e7336d5d55b/bonjour.txt]
    [[id:sampleC, character:sheep, lang:de, lang_group:germanic], /workspaces/training/side-quests/metadata/work/55/68c69c5efb527f3604ddb3daab8057/hallo.txt]
    [[id:sampleD, character:turkey, lang:en, lang_group:germanic], /workspaces/training/side-quests/metadata/work/2a/4752055ccb5d1370b0ef9da41d3993/hello.txt]
    [[id:sampleE, character:stegosaurus, lang:es, lang_group:romance], /workspaces/training/side-quests/metadata/work/f4/fcd3186dc666d5d239ffa6c37d125d/hola.txt]
    [[id:sampleF, character:moose, lang:fr, lang_group:romance], /workspaces/training/side-quests/metadata/work/c3/3b2627f733f278a7088332a5806108/salut.txt]
    [[id:sampleG, character:turtle, lang:it, lang_group:romance], /workspaces/training/side-quests/metadata/work/36/cca6a7dbfa26ac24f9329787a32e9d/ciao.txt]
    ```

This confirms we can now refer to the channel by name.

#### 3.2.2. Access the file and character metadata

We know from looking at the module code that the `COWPY` process expects to be given a text file and a `character` value.
To write the call to the `COWPY` process call, we just need to know how to extract the corresponding file object and metadata from each element in the channel.

As is often the case, the simplest way to do that is to use a `map` operation.

Our channel contains tuples structured as `[meta, file]`, so we can access the `file` object directly, and we can access the `character` value stored inside the meta map by referring to it as `meta.character`.

In the main workflow, make the following code changes:

=== "After"

    ```groovy title="main.nf" linenums="34"
        // Temporary: access the file and character
        ch_languages.map { meta, file -> file }.view { file -> "File: " + file }
        ch_languages.map { meta, file -> meta.character }.view { character -> "Character: " + character }
    ```

=== "Before"

    ```groovy title="main.nf" linenums="34"
        // Temporary: peek into ch_languages
        ch_languages.view()
    ```

Note that we're using closures (such as `{ file -> "File: " + file }`) to make the output of the `.view` operations more readable.

Let's run this:

```bash
nextflow run main.nf -resume
```

??? success "Command output"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `./main.nf` [cheesy_cantor] DSL2 - revision: 15af9c1ec7

    [43/05df08] IDENTIFY_LANGUAGE (7) [100%] 7 of 7, cached: 7 ✔
    Character: squirrel
    File: /workspaces/training/side-quests/metadata/work/8d/4b9498bbccb7a74f04e41877cdc3e5/bonjour.txt
    File: /workspaces/training/side-quests/metadata/work/d3/604274985406e40d79021dea658e60/guten_tag.txt
    Character: tux
    Character: turkey
    File: /workspaces/training/side-quests/metadata/work/d4/fafcc9415b61d2b0fea872e6a05e8a/hello.txt
    File: /workspaces/training/side-quests/metadata/work/02/468ac9efb27f636715e8144b37e9a7/hallo.txt
    Character: sheep
    Character: moose
    Character: stegosaurus
    File: /workspaces/training/side-quests/metadata/work/d4/61a7e1188b4f2742bc72004e226eca/salut.txt
    File: /workspaces/training/side-quests/metadata/work/ae/68364be238c11149c588bf6fc858b1/hola.txt
    File: /workspaces/training/side-quests/metadata/work/43/05df081af5d879ab52e5828fa0357e/ciao.txt
    Character: turtle
    ```

_The file paths and character values may come out in a different order in your output._

This confirms we're able to access the file and the character for each element in the channel.

#### 3.2.3. Call the `COWPY` process

Now let's put it all together and actually call the `COWPY` process on the `ch_languages` channel.

In the main workflow, make the following code changes:

=== "After"

    ```groovy title="main.nf" linenums="34"
        // Run cowpy to generate ASCII art
        COWPY(
            ch_languages.map { meta, file -> file },
            ch_languages.map { meta, file -> meta.character }
        )
    ```

=== "Before"

    ```groovy title="main.nf" linenums="34"
        // Temporary: access the file and character
        ch_languages.map { meta, file -> [file, meta.character] }
            .view()
    ```

You see we simply copy the two map operations (minus the `.view()` statements) as the inputs to the process call.
Just make sure you don't forget the comma between them!

It's a bit clunky, but we'll see how to make that better in the next section.

Let's run this:

```bash
nextflow run main.nf -resume
```

??? success "Command output"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [suspicious_crick] DSL2 - revision: 25541014c5

    executor >  local (7)
    [43/05df08] IDENTIFY_LANGUAGE (7) [100%] 7 of 7, cached: 7 ✔
    [e7/317c18] COWPY (6)             [100%] 7 of 7 ✔
    ```

If you look in the results directory, you should see the individual files containing the ASCII art of each greeting spoken by the corresponding character.

??? abstract "Directory and example file contents"

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

    ```text title="results/cowpy-bonjour.txt"
     _________________
    / Bonjour         \
    \ Salut, à demain /
    -----------------
      \
        \
                      _ _
          | \__/|  .~    ~.
          /oo `./      .'
          {o__,   \    {
            / .  . )    \
            `-` '-' \    }
          .(   _(   )_.'
          '---.~_ _ _|
    ```

This shows we were able to use the information in the meta map to parameterize the command in the second step of the pipeline.

However, as noted above, some of the code involved was a bit clunky, since we had to unpack meta data while still in the context of the workflow body.
That approach works fine for using a small number of fields from the meta map, but would scale poorly if we wanted to use a lot more.

There is another operator called `multiMap()` that allows us to streamline this a little bit, but even then it's not ideal.

??? info "(Optional) Alternative version with `multiMap()`"

    In case you're wondering, we couldn't just write a single `map()` operation that outputs both the `file` and the `character`, because that would return them as a tuple.
    We had to write two separate `map()` operations in order to feed the `file` and `character` elements to the process separately.

    Technically there is another way to do this through a single mapping operation, using the `multiMap()` operator, which is capable of emitting multiple channels.
    For example, you could replace the call to `COWPY` above with the following code:

    === "After"

        ```groovy title="main.nf" linenums="34"
            // Run cowpy to generate ASCII art
            COWPY(
                ch_languages.multiMap { meta, file ->
                    file: file
                    character: meta.character
                }
            )
        ```

    === "Before"

        ```groovy title="main.nf" linenums="34"
            // Run cowpy to generate ASCII art
            COWPY(
                ch_languages.map { meta, file -> file },
                ch_languages.map { meta, file -> meta.character }
            )
        ```

    This produces exactly the same result.

In either case, it's awkward that we have to do some unpacking at the workflow level.

It would be better if we could feed the entire meta map into the process and pick what we needed once there.

### 3.3. Pass and use the entire meta map

The point of the meta map is after all to pass all the metadata together as a bundle.
The only reason we couldn't do that above is that the process is not set up to accept a meta map.
But since we control the process code, we can change that.

Let's modify the `COWPY` process to accept the `[meta, file]` tuple structure that we used in the first process so we can streamline the workflow.

To that end, we'll need to do three things:

1. Modify the `COWPY` process module's input definitions
2. Update the process command to use the meta map
3. Update the process call in the workflow body

Ready? Let's go!

#### 3.3.1. Modify the `COWPY` module input

Make the following edits to the `cowpy.nf` module file:

=== "After"

    ```groovy title="cowpy.nf" linenums="10" hl_lines="2"
    input:
    tuple val(meta), path(input_file)
    ```

=== "Before"

    ```groovy title="cowpy.nf" linenums="10" hl_lines="2-3"
    input:
    path(input_file)
    val character
    ```

This enables us to use the `[meta, file]` tuple structure we covered earlier in the tutorial.

Note that we did not update the process output definition to output the meta map, in order to keep the tutorial streamlined, but feel free do that yourself as an exercise following the model of the `IDENTIFY_LANGUAGE` process.

#### 3.3.2. Update the command to use the meta map field

The entire meta map is now available inside the process, so we can refer to the information it contains directly from inside the command block.

Make the following edits to the `cowpy.nf` module file:

=== "After"

    ```groovy title="cowpy.nf" linenums="16" hl_lines="3"
    script:
    """
    cat ${input_file} | cowpy -c ${meta.character} > cowpy-${input_file}
    """
    ```

=== "Before"

    ```groovy title="cowpy.nf" linenums="16" hl_lines="3"
    script:
    """
    cat ${input_file} | cowpy -c ${character} > cowpy-${input_file}
    """
    ```

We've replaced the reference to the `character` value previously passed as a standalone input with the value held in the meta map, which we refer to using `meta.character`.

Now let's update the process call accordingly.

#### 3.3.3. Update the process call and run it

The process now expects its input to use the `[meta, file]` tuple structure, which is what the previous process outputs, so we can simply pass the whole `ch_languages` channel to the `COWPY` process.

Make the following edits to the main workflow:

=== "After"

    ```groovy title="main.nf" linenums="34" hl_lines="2"
    // Run cowpy to generate ASCII art
    COWPY(ch_languages)
    ```

=== "Before"

    ```groovy title="main.nf" linenums="34" hl_lines="3-4"
    // Run cowpy to generate ASCII art
    COWPY(
        ch_languages.map { meta, file -> file },
        ch_languages.map { meta, file -> meta.character }
    )
    ```

That simplifies the call significantly!

Let's delete the results of the previous execution and run it:

```bash
rm -r results
nextflow run main.nf
```

??? success "Command output"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [wise_sammet] DSL2 - revision: 99797b1e92

    executor >  local (14)
    [5d/dffd4e] process > IDENTIFY_LANGUAGE (7) [100%] 7 of 7 ✔
    [25/9243df] process > COWPY (7)             [100%] 7 of 7 ✔
    ```

If you look in the results directory, you should see the same outputs as previously, _i.e._ individual files containing the ASCII art of each greeting spoken by the corresponding character.

??? abstract "Directory contents"

    ```console
    ./results/
    ├── cowpy-bonjour.txt
    ├── cowpy-ciao.txt
    ├── cowpy-guten_tag.txt
    ├── cowpy-hallo.txt
    ├── cowpy-hello.txt
    ├── cowpy-hola.txt
    └── cowpy-salut.txt
    ```

So this produces the same results as before with simpler code.

Of course, this assumes you are able to modify the process code.
In some cases, you may have to rely on existing processes that you're not at liberty to modify, which limits your options.
The good news, if you're planning to use modules from the [nf-core](https://nf-co.re/) project, is that nf-core modules are all set up to use the `[meta, file]` tuple structure as a standard.

### 3.4. Troubleshooting missing required inputs

The `character` value is required for the `COWPY` process to run successfully.
If we do not set a default value for it in a configuration file, we MUST provide a value for it in the datasheet.

**What happens if we do not?**
It depends on what the input datasheet contains and which version of the workflow we're running.

#### 3.4.1. The character column exists but is empty

Let's say we delete the character value for one of the entries in our datasheet to simulate a data collection error:

```csv title="datasheet.csv" linenums="1" hl_lines="2"
id,character,recording
sampleA,,/workspaces/training/side-quests/metadata/data/bonjour.txt
sampleB,tux,/workspaces/training/side-quests/metadata/data/guten_tag.txt
sampleC,sheep,/workspaces/training/side-quests/metadata/data/hallo.txt
sampleD,turkey,/workspaces/training/side-quests/metadata/data/hello.txt
sampleE,stegosaurus,/workspaces/training/side-quests/metadata/data/hola.txt
sampleF,moose,/workspaces/training/side-quests/metadata/data/salut.txt
sampleG,turtle,/workspaces/training/side-quests/metadata/data/ciao.txt
```

For either version of the workflow we've used above, the `character` key will be created for all entries when the datasheet is read in, but for `sampleA` the value will be an empty string.

This will cause an error.

??? failure "Command output"

    ```console hl_lines="8 11 16 28"
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [marvelous_hirsch] DSL2 - revision: 0dfeee3cc1

    executor >  local (9)
    [c1/c5dd4f] process > IDENTIFY_LANGUAGE (7) [ 85%] 6 of 7
    [d3/b7c415] process > COWPY (2)             [  0%] 0 of 6
    ERROR ~ Error executing process > 'COWPY (1)'

    Caused by:
      Process `COWPY (1)` terminated with an error exit status (2)


    Command executed:

      cat bonjour.txt | cowpy -c  > cowpy-bonjour.txt

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

When Nextflow runs the `cowpy` command line for that sample, `${meta.character}` is filled with an empty string in the `cowpy` command line, so the `cowpy` tool throws an error saying no value was provided for the `-c` argument.

#### 3.4.2. The character column does not exist in the datasheet

Now let's say we delete the `character` column entirely from our datasheet:

```csv title="datasheet.csv" linenums="1"
id,recording
sampleA,/workspaces/training/side-quests/metadata/data/bonjour.txt
sampleB,/workspaces/training/side-quests/metadata/data/guten_tag.txt
sampleC,/workspaces/training/side-quests/metadata/data/hallo.txt
sampleD,/workspaces/training/side-quests/metadata/data/hello.txt
sampleE,/workspaces/training/side-quests/metadata/data/hola.txt
sampleF,/workspaces/training/side-quests/metadata/data/salut.txt
sampleG,/workspaces/training/side-quests/metadata/data/ciao.txt
```

In this case the `character` key will not be created at all when the datasheet is read in.

##### 3.4.2.1. Value accessed at the workflow level

If we're using the version of the code we wrote in section 3.2, Nextflow will attempt to access the `character` key in the meta map BEFORE calling the `COWPY` process.

It will not find any elements that match the instruction, so it will not run `COWPY` at all.

??? success "Command output"

    ```console hl_lines="7"
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [desperate_montalcini] DSL2 - revision: 0dfeee3cc1

    executor >  local (7)
    [1a/df2544] process > IDENTIFY_LANGUAGE (7) [100%] 7 of 7 ✔
    [-        ] process > COWPY                 -
    ```

As far as Nextflow is concerned, this workflow ran successfully!
However, none of the outputs we want will be produced.

##### 3.4.2.2. Value accessed at the process level

If we're using the version in section 3.3, Nextflow will pass the entire meta map to the `COWPY` process and attempt to run the command.

This will cause an error, but a different one compared to the first case.

??? failure "Command output"

    ```console hl_lines="8 11 16"
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [jovial_bohr] DSL2 - revision: eaaf375827

    executor >  local (9)
    [0d/ada9db] process > IDENTIFY_LANGUAGE (5) [ 85%] 6 of 7
    [06/28065f] process > COWPY (2)             [  0%] 0 of 6
    ERROR ~ Error executing process > 'COWPY (2)'

    Caused by:
      Process `COWPY (2)` terminated with an error exit status (1)


    Command executed:

      cat guten_tag.txt | cowpy -c null > cowpy-guten_tag.txt

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

This happens because `meta.character` does not exist, so our attempt to access it returns `null`. As a result, Nextflow literally plugs in `null` into the command-line, which is of course not recognized by the `cowpy` tool.

#### 3.4.3. Solutions

Aside from supplying a default value as part of the workflow configuration, there are two things we can do to handle this more robustly:

1. Implement input validation to your workflow to ensure that the datasheet contains all the required information. You can find an [introduction to input validation](../hello_nf-core/05_input_validation.md) in the Hello nf-core training course. <!-- TODO (future) pending a proper Validation side quest -->

2. If you want to make sure anyone who uses your process module can immediately identify required inputs, you can also make the required metadata property an explicit input.

Here's an example of how that would work.

First, at the process level, update the input definition as follows:

=== "After"

    ```groovy title="cowpy.nf" linenums="12" hl_lines="2"
        input:
        tuple val(meta), val(character), path(input_file)
    ```

=== "Before"

    ```groovy title="cowpy.nf" linenums="12" hl_lines="2"
        input:
        tuple val(meta), path(input_file)
    ```

Then, at the workflow level, use a mapping operation to extract the `character` property from the metadata and make it an explicit component of the input tuple:

=== "After"

    ```groovy title="main.nf" linenums="37" hl_lines="1"
        COWPY(ch_languages.map{meta, file -> [meta, meta.character, file]})
    ```

=== "Before"

    ```groovy title="main.nf" linenums="37" hl_lines="1"
        COWPY(ch_languages)
    ```

This approach has the advantage of showing explicitly that `character` is required, and makes the process easier to redeploy in other contexts.

This highlights an important design principle:

**Use the meta map for optional, descriptive information, but extract required values as explicit inputs.**

The meta map is excellent for keeping channel structures clean and preventing arbitrary channel structures, but for mandatory elements that are directly referenced in a process, extracting them as explicit inputs creates more robust and maintainable code.

### Takeaway

In this section, you've learned how to utilize metadata to customize the execution of a process, accessing it either at the workflow level or at the process level.

---

## Supplemental exercise

If you'd like to practice using meta map information from inside a process, try using other pieces of information from the meta map such as `lang` and `lang_group` to customize how the outputs are named and/or organized.

For example, try to modify the code to produce this result:

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

<!-- TODO (future) Provide worked out solution -->
<!-- the renaming should use the meta inside the process -->
<!-- the output org should use the meta in the workflow outputs -->

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

3.  **Customizing Process Behavior:** Using metadata inside the process.

    ```groovy
    cat $input_file | cowpy -c ${meta.character} > cowpy-${input_file}
    ```

### Additional resources

- [map](https://www.nextflow.io/docs/latest/operator.html#map)
- [stdout](https://www.nextflow.io/docs/latest/process.html#outputs)

---

## What's next?

Return to the [menu of Side Quests](./index.md) or click the button in the bottom right of the page to move on to the next topic in the list.
