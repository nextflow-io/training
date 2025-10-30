# Metadata in workflows

In any scientific analysis, we rarely work with just the raw data files. Each file comes with its own additional information: what it is, where it came from, and what makes it special. This extra information is what we call metadata.

Metadata is data describing other data.
Metadata tracks important details about files and experimental conditions, and helps tailor analyses to each dataset's unique characteristics.

Think of it like a library catalog: while books contain the actual content (raw data), the catalog cards provide essential information about each book—when it was published, who wrote it, where to find it (metadata). In Nextflow pipelines, metadata can be used to:

- Track file-specific information throughout the workflow
- Configure processes based on file characteristics
- Group related files for joint analysis

We'll explore how to handle metadata in workflows. Starting with a simple datasheet containing basic file information, you'll learn how to:

- Read and parse file metadata from CSV files
- Create and manipulate metadata maps
- Add new metadata fields during workflow execution
- Use metadata to customize process behavior

These skills will help you build more robust and flexible pipelines that can handle complex file relationships and processing requirements.

In this tutorial, we'll tackle a common data processing scenario: handling multiple files where each file needs different processing based on its characteristics.

We have greeting files in different languages (French, German, Spanish, Italian, English), but we don't know which language each file contains. Our challenge is to:

1. **Identify** the language in each file automatically
2. **Group** files by language family (Germanic vs Romance languages)
3. **Customize** the processing for each file based on its language and metadata
4. **Organize** outputs by language group

This represents a typical workflow pattern where file-specific metadata drives processing decisions - exactly the kind of problem that metadata maps solve elegantly.

Let's dive in!

## 0. Warmup

### 0.1. Prerequisites

Before taking on this side quest you should:

- Complete the [Hello Nextflow](../hello_nextflow/README.md) tutorial
- Understand basic Nextflow concepts (processes, channels, operators)

### 0.2. Starting Point

Let's move into the project directory.

```bash title="Navigate to the project directory"
cd side-quests/metadata
```

You can set VSCode to focus on this directory:

```bash title="Open VSCode in current directory"
code .
```

You'll find a `data` directory containing a samplesheet and a main workflow file.

```console title="Directory contents"
> tree
.
├── data
│   ├── bonjour.txt
│   ├── ciao.txt
│   ├── guten_tag.txt
│   ├── hallo.txt
│   ├── hello.txt
│   ├── hola.txt
│   ├── salut.txt
│   └── samplesheet.csv
├── main.nf
└── nextflow.config
```

The datasheet contains information about different files and some associated data that we will use in this exercise to tailor our analysis to each file. Each data file contains greetings in different languages, but we don't know what language they are in. The datasheet has 3 columns:

- `id`: self-explanatory, an ID given to the file
- `character`: a character name, that we will use later to draw different creatures
- `data`: paths to `.txt` files that contain phrases in different languages

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

## 1. Read in datasheet

### 1.1. Read in datasheet with splitCsv

Let's start by reading in the datasheet with `splitCsv`. In the main workflow file, you'll see that we've already started the workflow:

```groovy title="main.nf" linenums="1"
workflow  {

    ch_samplesheet = channel.fromPath("./data/samplesheet.csv")

}
```

!!! note

    Throughout this tutorial, we'll use the `ch_` prefix for all channel variables to clearly indicate they are Nextflow channels.

=== "After"

    ```groovy title="main.nf" linenums="3" hl_lines="2-3"
        ch_samplesheet = channel.fromPath("./data/samplesheet.csv")
            .splitCsv(header: true)
            .view()
    ```

=== "Before"

    ```groovy title="main.nf" linenums="3"
        ch_samplesheet = channel.fromPath("./data/samplesheet.csv")
    ```

We can use the [`splitCsv` operator](https://www.nextflow.io/docs/latest/operator.html#splitcsv) to split the datasheet into a channel of maps, where each map represents a row from the CSV file.

A map is a key-value data structure similar to dictionaries in Python, objects in JavaScript, or hashes in Ruby. For example:

```groovy
// Groovy map
def my_map = [id:'sampleA', character:'squirrel']
println my_map.id  // Prints: sampleA
```

!!! example "Try it yourself"

      You can run this example to see how maps look like with:

      ```bash title="Run map demo example"
      nextflow run examples/map_demo.nf
      ```

The `header: true` option tells Nextflow to use the first row of the CSV file as the header row, which will be used as keys for the values. Let's see what Nextflow can see after reading with `splitCsv`. Run the pipeline with the `view()` operator we added above:

Run the pipeline:

```bash title="Read the datasheet"
nextflow run main.nf
```

```console title="Read datasheet with splitCsv"
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

We can see that each row from the CSV file has been converted into a map with keys matching the header row.
Each map entry corresponds to a column in our datasheet:

- `id`
- `character`
- `recording`

This format makes it easy to access specific fields from each file. For example, we could access the file ID with `id` or the txt file path with `recording`. The output above shows each row from the CSV file converted into a map with keys matching the header row.

Let's access a specific column, the `character` column, that we read in from the datasheet, and print it. We can use the Nextflow `map` operator to iterate over each item in our channel and return specific entry of our map object:

=== "After"

    ```groovy title="main.nf" linenums="4" hl_lines="2-5"
            .splitCsv(header: true)
            .map{ row ->
                row.character
            }
            .view()
    ```

=== "Before"

    ```groovy title="main.nf" linenums="4" hl_lines="2-3"
            .splitCsv(header: true)
            .view()
    ```

and rerun:

```bash title="Print the creatures"
nextflow run main.nf
```

```console title="Print the creatures"
squirrel
tux
sheep
turkey
stegosaurus
moose
turtle
```

Success, we can use the map structures derived from our datasheet to access the values from individual columns for each row.

Now that we've successfully read in the datasheet and have access to the data in each row, we can begin implementing our pipeline logic.

### 1.2. Separate metadata and data

In the datasheet, we have both the input files and data about the input files (`id`, `character`), the metadata.
As we progress through the workflow, we generate more metadata about each file.
So far, every column from the datasheet has become a separate item in the channel items we derived using `splitCsv()`. Every process that consumes this channel would need to be configured with this structure in mind:

```groovy
    input:
    tuple val(id), val(character), file(recording)
```

This means that as soon as somebody added an extra column to the datasheet, the workflow would start producing errors, because the shape of the channel would no longer match what the process expected. It also makes the process hard to share with others who might have slightly different input data, and we end up hard-coding variables into the process that aren't needed by the script block.

To avoid this, we need to find a way of keeping the channel structure consistent however many columns that datasheet contains, and we can do that by keeping all the metadata in a single part of the channel tuples we call simply the `meta map`.

Let's use this and separate our metadata from the file path. We'll use the `map` operator to restructure our channel elements into a tuple consisting of the meta map and file:

=== "After"

    ```groovy title="main.nf" linenums="5" hl_lines="2"
        .map { row ->
                [ [id: row.id, character: row.character], row.recording ]
        }
    ```

=== "Before"

    ```groovy title="main.nf" linenums="5" hl_lines="2"
        .map { row ->
                row.character
        }
    ```

Let's run it:

```bash title="View meta map"
nextflow run main.nf
```

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

Now, each item in the channel has only two items, the metadata map (e.g. `[id:sampleA, character:squirrel]`) and the data file described by that metadata (e.g. `/workspaces/training/side-quests/metadata/data/bonjour.txt`).

Now we can write processes to consume the channel without hard-coding the metadata items into the input specification:

```groovy
    input:
    tuple val(meta), file(recording)
```

Additional columns in the datasheet will make more metadata available in the `meta` map, but won't change the channel shape.

### Takeaway

In this section, you've learned:

- **Why metadata is important**: Keeping metadata with your data preserves important file information throughout the workflow.
- **How to read in a datasheets**: Using `splitCsv` to read CSV files with header information and transform rows into structured data
- **How to create a meta map**: Separating metadata from file data using tuple structure `[ [id:value, ...], file ]`

---

## 2. Manipulating metadata

### 2.1. Copying input metadata to output channels

Now we want to process our files with unidentified languages. Let's add a process definition before the `workflow` that can identify the language in each file:

=== "After"

    ```groovy title="main.nf" linenums="1"
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

    workflow {
    ```

=== "Before"

    ```groovy title="main.nf" linenums="1"
    workflow  {
    ```

The tool [langid](https://github.com/saffsd/langid.py) is used for language identification. It comes pre-trained on a set of languages. For a given phrase, it outputs a language prediction and a probability score for each guess to the console. In the `script` section, we use sed to remove the probability score, clean up the string by removing newline characters, and return only the language prediction. Since the output is printed directly to the console, we use Nextflow’s [`stdout` output qualifier](https://www.nextflow.io/docs/latest/process.html#outputs) to capture and pass the string as output.

Let's include the process, then run, and view it:

=== "After"

    ```groovy title="main.nf" linenums="25" hl_lines="4-5"
                [ [id: row.id, character: row.character], row.recording ]
            }

        ch_prediction = IDENTIFY_LANGUAGE(ch_samplesheet)
        ch_prediction.view()
    ```

=== "Before"

    ```groovy title="main.nf" linenums="25" hl_lines="3"
                [ [id:row.id, character:row.character], row.recording ]
            }
            .view()
    ```

```bash title="Identify languages"
nextflow run main.nf
```

```console title="Identify languages"
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

Neat, for each of our files, we now have a language predicted.

!!! note

    "de" stands for "deutsch", the German word for "german"

You may have noticed something else: we kept the metadata of our files and associated it with our new piece of information. We achieved this by adding the `meta` map to the output tuple in the process:

```groovy title="main.nf" linenums="12"
output:
tuple val(meta), path(file), stdout
```

This is a useful way to ensure the metadata stays connected with any new information that is generated.

!!! note

    Another compelling reason to use meta maps in this way is that they make it easier to associate related results that share the same identifiers. As you learned in "Hello Nextflow", you can't rely on the order of items in channels to match results across them. Instead, you must use keys to associate data correctly - and meta maps provide an ideal structure for this purpose. We explore this use case in detail in [Splitting & Grouping](./splitting_and_grouping.md).

### 2.2. Using process outputs to augment metadata

Given that this is more data about the files, let's add it to our meta map. We can use the [`map` operator](https://www.nextflow.io/docs/latest/operator.html#map) again to create a new key `lang` and set the value to the predicted language:

=== "After"

    ```groovy title="main.nf" linenums="28" hl_lines="2-6"
        ch_prediction = IDENTIFY_LANGUAGE(ch_samplesheet)
        ch_languages = ch_prediction
            .map { meta, file, lang ->
                [meta + [lang: lang], file]
            }
            .view()
    }
    ```

=== "Before"

    ```groovy title="main.nf" linenums="28" hl_lines="2"
        ch_prediction = IDENTIFY_LANGUAGE(ch_samplesheet)
        ch_prediction.view()
    ```

The `map` operator takes each channel element and processes it to create a modified version. Inside the closure `{ meta, file, lang -> ... }`, we then take the existing `meta` map, create a new map `[lang:lang]`, and merge both together using `+`.

The `+` operator in Groovy merges two maps together. So if our original `meta` was `[id:sampleA, character:squirrel]`, then `meta + [lang:'fr']` creates a new map: `[id:sampleA, character:squirrel, lang:fr]`.

!!! Note

    The `+` notation with maps creates an entirely new map object, which is what we want. If we'd done something like `meta.lang = lang` we'd have been modifying the original object, which can lead to unpredictable effects.

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

Nice, we expanded our meta map with new information we gathered in the pipeline. After running the language prediction, each element in the output channel looks like this:

```console
[meta, file, lang]  // e.g. [[id:sampleA, character:squirrel], bonjour.txt, fr]
```

<!-- TODO Should we also show how to remove a key using subMap?! -->

### 2.3. Assign a language group using conditionals

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

Let's understand how this transformation works. The `map` operator here again takes a closure that processes each element in the channel. Inside the closure, we're using an if-clause to create a new language group classification.

Here's what's happening step by step:

- Create a new field `lang_group`: We set a default to `unknown`
- Extract existing metadata: We access `meta.lang` (the language we predicted earlier) from the existing meta map
- Apply conditional logic: We use an if-clause to determine the language group based on the language: Is `meta.lang` either `de` or `en`, we re-assign `lang_group` to `germanic`, if `fr`, `es`, or `it`, then we re-assign to `romance`
- Merge with existing metadata: We use `meta + [lang_group:lang_group]` in the same way as before to combine the existing meta map with our new field

The resulting channel elements maintain their `[meta, file]` structure, but the meta map now includes this new classification.

### Takeaway

In this section, you've learned how to :

- **Apply input metadata to output channels**: Copying metadata in this way allows us to associate results later on based on metadata content.
- **Create custom keys**: You created two new keys in your meta map, merging them with `meta + [new_key:value]` into the existing meta map. One based on a computed value from a process, and one based on a condition you set in the `map` operator.

These allow you to associate new and existing metadata with files as you progress through your pipeline.

---

## 3. Use meta map information in a process

Let's make some fun characters say the phrases from the files our channel. In the [hello-nextflow training](../hello_nextflow/05_hello_containers.md), you already encountered the `cowpy` package, a python implementation of a tool called `cowsay` that generates ASCII art to display arbitrary text inputs in a fun way. We will re-use a process from there.

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

The process definition as provided would direct results to the `results` folder, but let's make a tweak to be a little smarter. Given we have been trying to figure out what languages our samples were in, let's group the samples by language in the output directory.
Earlier, we added the predicted language to the `meta` map. We can access this `key` in the process and use it in the `publishDir` directive:

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

Look through the other files. All phrases should be spoken by the fashionable stegosaurus.

How did this work? The `publishDir` directive is evaluated at runtime when the process executes. Each process task gets its own meta map from the input tuple When the directive is evaluated, `${meta.lang_group}` is replaced with the actual group language value for that dataset creating the dynamic paths like `results/romance`.

### 3.2. Customize the character

In our datasheet, we have another column: `character`. To tailor the tool parameters per file, we can also access information from the `meta` map in the script section. This is really useful in cases where a tool should have different parameters for each file.

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

This approach differs from using pipeline parameters (`params`), which generally apply the same configuration to all files in your workflow. By leveraging metadata applied to each item in a channel, you can fine-tune process behavior on a per-file basis.

#### 3.2.1. Exploiting metadata at the workflow level

In the example above, by using a property of the meta map in the script block, we introduce a hard requirement on the properties that must be present. Anyone running with a sample sheet that did not contain the `character` property would encounter an error. The process `input:` only says that the `meta` map is required, so someone trying to use this process in another workflow might not notice immediately that the `character` property was required.

A better approach is to make the required metadata property an explicit input rather than accessing it from within the meta map. This makes the requirement clear and provides better error messages. Here's how to refactor the COWPY process:

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

This highlights an important design principle: Use the meta map for optional, descriptive information, but extract required values as explicit inputs. The meta map is excellent for keeping channel structures clean and preventing arbitrary channel structures, but for mandatory elements that are directly referenced in a process, extracting them as explicit inputs creates more robust and maintainable code.

### Takeaway

In this section, you've learned how to:

- **Tweak directives using meta values**: Using meta map values in `publishDir` directives to create dynamic output paths based on the file's metadata

- **Tweak the script section based on meta values**: Customizing tool parameters per file using meta information in the `script` section

---

## Summary

In this side quest, you've explored how to effectively work with metadata in Nextflow workflows. This pattern of keeping metadata explicit and attached to the data is a core best practice in Nextflow that enables building robust, maintainable bioinformatics workflows.

Here's what you've learned:

1. **Reading and Structuring Metadata**: Reading CSV files and creating organized metadata maps that stay associated with your data files

2. **Expanding Metadata During Workflow**: Adding new information to your metadata as your pipeline progresses by adding process outputs and deriving values through conditional logic

3. **Customizing Process Behavior**: Using metadata to adapt how processes handle different files

This approach offers several advantages over hardcoding file information:

- File metadata stays associated with files throughout the workflow
- Process behavior can be customized per file
- Output organization can reflect file metadata
- File information can be expanded during pipeline execution

### Key Concepts

- **Reading Datasheets & creating meta maps**

  ```nextflow
  channel.fromPath('samplesheet.csv')
    .splitCsv(header: true)
    .map { row ->
        [ [id:row.id, character:row.character], row.recording ]
    }
  ```

- **Adding new keys to the meta map**

  1.  based on process output:

```nextflow
.map { meta, file, lang ->
  [ meta + [lang:lang], file ]
}
```

2.  and using a conditional clause

```nextflow
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

- **Using meta values in Process Directives**

  ```nextflow
  publishDir "results/${meta.lang_group}", mode: 'copy'
  ```

- **Adapting tool parameters for individual files**

  ```nextflow
  cat $input_file | cowpy -c ${meta.character} > cowpy-${input_file}
  ```

## Resources

- [map](https://www.nextflow.io/docs/latest/operator.html#map)
- [stdout](https://www.nextflow.io/docs/latest/process.html#outputs)
