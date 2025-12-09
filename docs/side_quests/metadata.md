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

```groovy
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

```groovy
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
So instead, we're going to create a new meta map containing the contents of the existing meta map plus a new key-value pair holding the new information, and use the [`map`](https://www.nextflow.io/docs/latest/operator.html#map) operator to replace the old map with the new one.

This is going to take only a very small amount of code, but it's going to have a lot packed into it, so let's break it down in stages.

First, you need to know that we can merge the contents of two maps using the Groovy operator `+`.
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
[id: 'sampleA', character: 'squirrel', 'fr']
```

Great! But now let's say the starting point is a little different: you still have `map1` but the language prediction is not in its own map.
It's held in a variable called `lang_id`, and you know you want to store its value (`'fr'`) with the key `lang`.

You can actually do the following:

```groovy
new_map = [map1 + [lang: lang_id]]
```

Here, `[lang: new_info]` creates a new unnamed map on the fly, and `map1 + ` merges `map1` with the new unnamed map, producing the same `new_map` contents as before.

Neat, right?

Now let's transpose that into the context of a Nextflow `channel.map()` operation.
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
So we've effectively transformed `[id: 'sampleA', character: 'squirrel'], 'it'` into `[id: 'sampleA', character: 'squirrel', lang: 'fr']`.

Hopefully you can see that if we change `map1` to `meta`, that's basically all we need in order to add the language predication to our meta map in our workflow.

Except for one thing!

In the case of our workflow, we also need to account for the presence of the `file` object in the tuple, which contains `meta, file, lang_id`, so the code here would become:

```groovy
.map { meta, file, lang_id ->
    [meta + [lang: lang_id], file]
}
```

If you're having a hard time following why the `file` seems to be moving around, imagine that instead of `[meta + [lang: lang_id], file]`, that line reads `[new_map, file]`.
This should make it more clear that we're simply leaving the `file` in its original place in second position in the tuple; we've just taken the `new_info` value and folded it into the map that's in first position.

And this brings us back to the `tuple val(meta), path(file)` channel structure!

With that, let's make the following edits to the workflow:

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

See, that's not a lot of code, but it helps to walk through the steps before you plug it in.

Let's run the workflow to see if it worked:

```bash
nextflow run main.nf -resume
```

??? example "Output"

    ```console

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

Your goal is to make the following changes to the workflow:

=== "After"

    ```groovy title="main.nf" linenums="14" hl_lines="7-19"
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

    ```groovy title="main.nf" linenums="14" hl_lines="7"
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

??? example "Output"

    ```console
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

As you can see, the channel elements maintain their `[meta, file]` structure, but the meta map now includes this new classification.

### Takeaway

In this section, you've learned how to :

- **Apply input metadata to output channels**: Copying metadata in this way allows us to associate results later on based on metadata content.
- **Create custom keys**: You created two new keys in your meta map, merging them with `meta + [new_key:value]` into the existing meta map.One based on a computed value from a process, and one based on a condition you set in the `map` operator.

These allow you to associate new and existing metadata with files as you progress through your pipeline.
Even if you're not using metadata as part of a process, keeping the meta map associated with the data like this makes it easy to keep all the relevant information together.

---

## 3. Use meta map information in a process

Now that you know how to create and update the meta map, we can get to the really fun bit: actually using the metadata in a process.

More specifically, we're going to add a second step to our workflow to draw each animal as ASCII art and make it say the recorded text in a speech bubble.
We're going to do this using a tool called [`cowpy`](https://github.com/jeffbuttars/cowpy).

<details>
  <summary>What does `cowpy` do?</summary>

`cowpy` is a command-line tool that generates ASCII art to display arbitrary text inputs in a fun way.
It is a python implementation of the classic [`cowsay`](https://en.wikipedia.org/wiki/Cowsay) tool by Tony Monroe.

```console
> cowpy "Hello Nextflow"
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
> cowpy "Hello Nextflow" -c tux
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

</details>

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

```groovy title="modules/cowpy.nf" linenums="1" hl_lines=""
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

<!-- TODO Need to continue adapting the instructions from here after the reordering of steps / slight change in code and logic. -->

### 3.2. Use meta map information to set the character

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

### 3.2. Use meta map information to organize results

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
