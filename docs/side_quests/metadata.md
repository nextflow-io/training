# Working with sample-specific data

In any scientific analysis, we rarely work with just the raw data files. Each sample comes with its own additional information: what it is, where it came from, and what makes it special. This extra information is what we call sample-specific data.

Instead of saying “sample-specific data” (which is a mouthful), it is often called "metadata" instead: data describing other data.
Metadata tracks important details about samples and experimental conditions, and helps tailor analyses to each dataset’s unique characteristics.

Think of it like a library catalog: while books contain the actual content (raw data), the catalog cards provide essential information about each book—when it was published, who wrote it, where to find it (metadata). In Nextflow pipelines, metadata can be used to:

- Track sample-specific information throughout the workflow
- Configure processes based on sample characteristics
- Group related samples for joint analysis

We'll explore how to handle metadata in workflows. Starting with a simple samplesheet containing basic sample information, you'll learn how to:

- Read and parse sample metadata from CSV files
- Create and manipulate metadata maps
- Add new metadata fields during workflow execution
- Use metadata to customize process behavior

These skills will help you build more robust and flexible pipelines that can handle complex sample relationships and processing requirements.

Let's dive in!

## 0. Warmup

### 0.1. Prerequisites

Before taking on this side quest you should:

- Complete the [Hello Nextflow](../hello_nextflow/README.md) tutorial
- Understand basic Nextflow concepts (processes, channels, operators)

### 0.2. Starting Point

Let's move into the project directory.

```bash
cd side-quests/metadata
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

The samplesheet contains information about different samples and some associated data that we will use in this exercise to tailor our analysis to each sample. Each data files contains greetings in different languages, but we don't know what language they are in. In particular, the samplesheet has 3 columns:

- `id`: self-explanatory, an ID given to the sample
- `character`: a character name, that we will use later to draw different creatures
- `data`: paths to `.txt` files that contain phrases in different languages

```console title="samplesheet.csv"
id,character,recording
sampleA,squirrel,/workspaces/training/side-quests/metadata/data/bonjour.txt
sampleB,tux,/workspaces/training/side-quests/metadata/data/guten_tag.txt
sampleC,sheep,/workspaces/training/side-quests/metadata/data/hallo.txt
sampleD,turkey,/workspaces/training/side-quests/metadata/data/hello.txt
sampleE,stegosaurus,/workspaces/training/side-quests/metadata/data/hola.txt
sampleF,moose,/workspaces/training/side-quests/metadata/data/salut.txt
sampleG,turtle,/workspaces/training/side-quests/metadata/data/ciao.txt
```

## 1. Read in samplesheet

### 1.1. Read in samplesheet with splitCsv

Let's start by reading in the samplesheet with `splitCsv`. In the main workflow file, you'll see that we've already started the workflow:

```groovy title="main.nf" linenums="1"
workflow  {

    ch_samplesheet = Channel.fromPath("./data/samplesheet.csv")

}
```

!!! note

    Throughout this tutorial, we'll use the `ch_` prefix for all channel variables to clearly indicate they are Nextflow channels.

=== "After"

    ```groovy title="main.nf" linenums="3" hl_lines="2-3"
    ch_samplesheet = Channel.fromPath("./data/samplesheet.csv")
                            .splitCsv(header: true)
                            .view()
    ```

=== "Before"

    ```groovy title="main.nf" linenums="3"
    ch_samplesheet = Channel.fromPath("./data/samplesheet.csv")
    ```

We can use the [`splitCsv` operator](https://www.nextflow.io/docs/latest/operator.html#splitcsv) to split the samplesheet into a channel of maps, where each map represents a row from the CSV file.

A map is a key-value data structure similar to dictionaries in Python, objects in JavaScript, or hashes in Ruby. For example:

```groovy
// Groovy map
def my_map = [id:'sampleA', character:'squirrel']
println my_map.id  // Prints: sampleA
```

```python
# Python equivalent dictionary
my_map = {'id': 'sampleA', 'character': 'squirrel'}
print(my_map['id'])  # Prints: sampleA
```

The `header: true` option tells Nextflow to use the first row of the CSV file as the header row, which will be used as keys for the values. Let's see what Nextflow can see after reading with `splitCsv`. To do this, we can use the `view` operator.

Run the pipeline:

```bash title="Read the samplesheet"
nextflow run main.nf
```

```console title="Read samplesheet with splitCsv"
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
Each map entry corresponds to a column in our samplesheet:

- `id`
- `character`
- `recording`

This format makes it easy to access specific fields from each sample. For example, we could access the sample ID with `id` or the txt file path with `recording`. The output above shows each row from the CSV file converted into a map with keys matching the header row.

Let's access a specific column, the `character` column, that we read in from the samplesheet, and print it. We can use the Nextflow `map` operator to iterate over each `row` in our samplesheet and return specific entry of our map object:

=== "After"

    ```groovy title="main.nf" linenums="4" hl_lines="2-4"
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

Success, we can access individual entries or entire rows from our samplesheet.

Now that we've successfully read in the samplesheet and have access to the data in each row, we can begin implementing our pipeline logic.

### 1.2. Separate meta data and data

In the samplesheet, we have both the input files and data about the input files (`id`, `character`), the meta data.
As we progress through the workflow, we generate more meta data about each sample.
If we keep all fields as separate Channel elements, it quickly becomes messy: every time we add a new field, we’d have to update every downstream process to expect a different number of inputs, making the code brittle and hard to maintain.

By grouping the metadata into its own key-value map, and keeping the file path as a distinct element, we make our workflow much more robust and flexible. We can add or remove metadata at any stage without having to rewrite process inputs. This approach also keeps process code cleaner and makes it easier to share and reuse code, both within a workflow and across different workflows.

Let's use this and separate our metadata from the file path. We'll use the `map` operator to restructure our channel elements into a tuple consisting of the meta map and file:

=== "After"

  ```groovy title="main.nf" linenums="5" hl_lines="2"
  .map { row ->
    [ [id:row.id, character:row.character], row.recording ]
  }
  .view()
  ```

=== "Before"

    ```groovy title="main.nf" linenums="5"
    .map{ row ->
      row.character
    }
    ```

Let's run it:

```bash title="View meta map"
nextflow run main.nf
````

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

We have successfully separated our values into their own map to separate it from the file. Each of our channel elements now has the shape `[ meta,  file ]`. Each process that now uses the data can be of similar input shape even if we keep adding or removing fields.

### Takeaway

In this section, you've learned:

- **Why metadata is important**: Keeping metadata with your data preserves important sample information throughout the workflow.
- **How to read in a samplesheet**: Using `splitCsv` to read CSV files with header information and transform rows into structured data
- **How to create a meta map**: Separating metadata from file data using tuple structure `[ [id:value, ...], file ]`

---

## 2. Create new meta map keys

### 2.1. Passing the meta map through a process

Now we want to process our samples with unidentified languages. Let's add a process definition before the `workflow` that can identify the language in each file:

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
    ```

=== "Before"

    ```groovy title="main.nf" linenums="1"
    workflow  {
    ```

The tool [langid](https://github.com/saffsd/langid.py) is used for language identification. It comes pre-trained on a set of languages. For a given phrase, it outputs a language prediction and a probability score for each guess to the console. In the `script` section, we remove the probability score, clean up the string by removing newline characters, and return only the language prediction. Since the output is printed directly to the console, we use Nextflow’s [`stdout` output qualifier](https://www.nextflow.io/docs/latest/process.html#outputs) to capture and pass the string as output.

Let's include the process, then run, and view it:

=== "After"

    ```groovy title="main.nf" linenums="25" hl_lines="3-5"
                          [ [id:row.id, character:row.character], row.recording ]
                      }

          ch_prediction = IDENTIFY_LANGUAGE(ch_samplesheet)
          ch_prediction.view()
    ```

=== "Before"

    ```groovy title="main.nf" linenums="25"
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

Neat, for each of our samples, we now have a language predicted.

!!! note

    "de" stands for "deutsch", the German word for "german"

You may have noticed something else: we kept the meta data of our samples and associated it with our new piece of information. We achieved this by adding the `meta` map to the output tuple in the process:

```groovy title="main.nf" linenums="12"
output:
      tuple val(meta), path(file), stdout
```

This is a useful way to ensure the sample-specific information stays connected with any new information that is generated.

### 2.2. Add the language prediction to the meta map

Given that this is more data about the files, let's add it to our meta map. We can use the [`map` operator](https://www.nextflow.io/docs/latest/operator.html#map) again to create a new key `lang` and set the value to the predicted language:

=== "After"

    ```groovy title="main.nf" linenums="28" hl_lines="2-5"
      ch_prediction = IDENTIFY_LANGUAGE(ch_samplesheet)
      ch_languages = ch_prediction.map { meta, file, lang ->
                                      [ meta + [lang:lang], file ]
                                  }
                                  .view()

    }
    ```

=== "Before"

    ```groovy title="main.nf" linenums="28"
      ch_prediction = IDENTIFY_LANGUAGE(ch_samplesheet)
    ```

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

Nice, we expanded our meta map with new information we gathered in the pipeline. Let's take a look at what happened here:

After running the language prediction, each element in the output channel looks like this:

```console
[meta, file, lang]  // e.g. [[id:sampleA, character:squirrel], bonjour.txt, fr]
```

The `map` operator takes each channel element and processes it to create a modified version. Inside the closure `{ meta, file, lang -> ... }`, we then take the existing `meta` map, create a new map `[lang:lang]`, and merge both together using `+`.

The `+` operator in Groovy merges two maps together. So if our original `meta` was `[id:sampleA, character:squirrel]`, then `meta + [lang:'fr']` creates a new map: `[id:sampleA, character:squirrel, lang:fr]`.

<!-- TODO Should we also show how to remove a key using subMap?! -->

### 2.3. Assign a language group using conditionals

Alright, now that we have our language predictions, let's use the information to assign them into new groups. In our example data, we have provided data sets that belong either to `germanic` (either English or German) or `romance` (French, Spanish, Italian) languages.

We can use the `map` operator to assign either group.

=== "After"

    ```groovy title="main.nf" linenums="31" hl_lines="3-11"
    }
    .map { meta, file ->

        if ( meta.lang.equals("de") || meta.lang.equals('en') ){
            lang_group = "germanic"
        } else if ( meta.lang in ["fr", "es", "it"] ) {
            lang_group = "romance"
        } else {
            lang_group = "unknown"
        }

        [ meta + [lang_group:lang_group], file ]
    }
    .view()
    ```

=== "Before"

    ```groovy title="main.nf" linenums="31"
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

- Create a new field `lang_group`: We set a default to `romance`
- Extract existing metadata: We access `meta.lang` (the language we predicted earlier) from the existing meta map
- Apply conditional logic: We use an if-clause to determine the language group based on the language: Is `meta.lang` either `de` or `en`, then we re-assign `lang_group` to `germanic`
- Merge with existing metadata: We use `meta + [lang_group:lang_group]` to combine the existing meta map with our new field

The resulting channel elements maintain their `[meta, file]` structure, but the meta map now includes this new classification.

### Takeaway

In this section, you've learned how to :

- **Create custom keys**: You created two new keys in your meta map, merging them with `meta + [new_key:value]` into the existing meta map. One based on a computed value from a process, and one based on a condition you set in the `map` operator.

These allow you to associated new and existing meta data with files as you progress through your pipeline.

---

## 3. Customize a process with meta map

Let's let some fun characters say the phrases that we have passed in. In the [hello-nextflow training](../hello_nextflow/05_hello_containers.md), you already encountered the `cowpy` package, a python implementation of a tool called `cowsay` that generates ASCII art to display arbitrary text inputs in a fun way. We will re-use a process from there.

Copy in the process before your workflow block:

=== "After"

    ```groovy title="main.nf" linenums="20"
    /*
     * Generate ASCII art with cowpy
    */
    process COWPY {

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
    ```

=== "Before"

    ```groovy title="main.nf" linenums="20"
    workflow{
      ...
    }
    ```

### 3.1. Add a custom publishing location

Let's run our samples through `COWPY` and remove our `view` statement:

=== "After"

    ```groovy title="main.nf" linenums="59" hl_lines="3"
                                      [ meta + [lang_group:lang_group], file ]
                                  }
    COWPY(ch_languages)
    ```

=== "Before"

    ```groovy title="main.nf" linenums="59"
                                      [ meta + [lang_group:lang_group], file ]
                                  }
                                  .view()
    ```

We are still missing a publishing location. Given we have been trying to figure out what languages our samples were in, let's group the samples by language in the output directory. Earlier, we added the predicted language to the `meta` map. We can access this `key` in the process and use it in the `publishDir` directive:

=== "After"

    ```groovy title="main.nf" linenums="23" hl_lines="3"
    process COWPY {

        publishDir "results/${meta.lang_group}", mode: 'copy'

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
    ```

=== "Before"

    ```groovy title="main.nf" linenums="23"
    process COWPY {

      container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
    ```

Let's run this:

```bash title="Use cowpy"
nextflow run main.nf
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

How did this work? The `publishDir` directive is evaluated at runtime when the process executes. Each process task gets its own meta map from the input tuple When the directive is evaluated, `${meta.lang_group}` is replaced with the actual group language value for that sample creating the dynamic paths like `results/romance`.

### 3.2. Customize the character

In our samplesheet, we have another column: `character`. To tailor the tool parameters per sample, we can also access information from the `meta` map in the script section. This is really useful in cases were a tool should have different parameters for each sample.

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

This is a subtle difference to other parameters that we have set in the pipelines in previous trainings. A parameter that is passed as part of the `params` object is generally applied to all samples. When a more surgical approache is necessary, using the sample specific information is a good alternative.

### Takeaway

In this section, you've learned how to:

- **Tweak directives using meta values**: Using meta map values in `publishDir` directives to create dynamic output paths based on sample properties

- **Tweak the script section based on meta values**: Customizing tool parameters per sample using meta information in the `script` section

---

## Summary

In this side quest, you've explored how to effectively work with metadata in Nextflow workflows. This pattern of keeping metadata explicit and attached to the data is a core best practice in Nextflow that enables building robust, maintainable bioinformatics workflows.

Here's what you've learned:

1. **Reading and Structuring Metadata**: Reading CSV files and creating organized metadata maps that stay associated with your data files

2. **Expanding Metadata During Workflow**: Adding new information to your metadata as your pipeline progresses by adding process outputs and deriving values through conditional logic

3. **Customizing Process Behavior**: Using metadata to adapt how processes handle different samples

This approach offers several advantages over hardcoding sample information:

- Sample metadata stays associated with files throughout the workflow
- Process behavior can be customized per sample
- Output organization can reflect sample properties
- Sample information can be expanded during pipeline execution

### Key Concepts

- **Reading Samplesheets & creating meta maps**

  ```nextflow
  Channel.fromPath('samplesheet.csv')
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

- **Adapting tool parameters for individual samples**

  ```nextflow
  cat $input_file | cowpy -c ${meta.character} > cowpy-${input_file}
  ```

## Resources

- [map](https://www.nextflow.io/docs/latest/operator.html#map)
- [stdout](https://www.nextflow.io/docs/latest/process.html#outputs)
