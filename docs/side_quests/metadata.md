# Metadata

Metadata is crucial information that describes and gives context to your data. In workflows, it helps track important details about samples, experimental conditions, and processing parameters. Metadata is sample-specific and helps tailor analyses to each dataset's unique characteristics.

Think of it like a library catalog: while books contain the actual content (raw data), the catalog cards provide essential information about each book - when it was published, who wrote it, where to find it (metadata). In Nextflow pipelines, metadata helps us:

- Track sample-specific information throughout the workflow
- Configure processes based on sample characteristics
- Group related samples for joint analysis

In this side quest, we'll explore how to handle metadata effectively in Nextflow workflows. Starting with a simple samplesheet containing basic sample information, you'll learn how to:

- Read and parse sample metadata from CSV files
- Create and manipulate metadata maps
- Add new metadata fields during workflow execution
- Use metadata to customize process behavior

These skills will help you build more robust and flexible pipelines that can handle complex sample relationships and processing requirements.

Let's dive in and see how metadata can make our workflows smarter and more maintainable!

## 0. Warmup

### 0.1 Prerequisites

Before taking on this side quest you should:

- Complete the [Hello Nextflow](../hello_nextflow/README.md) tutorial
- Understand basic Nextflow concepts (processes, channels, operators)

### 0.2 Starting Point

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

The samplesheet contains information about different samples and some associated data that we will use in this exercise to tailor our analysis to each sample. In particular, the samplesheet has 3 columns:

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

Let's start by reading in the samplesheet with `splitCsv`. In the main workflow file, you'll see that we've already started the workflow.

```groovy title="main.nf" linenums="1"s
workflow  {

    ch_samplesheet = Channel.fromPath("./data/samplesheet.csv")
                            .splitCsv(header: true)

}
```

!!! note

    Throughout this tutorial, we'll use the `ch_` prefix for all channel variables to clearly indicate they are Nextflow channels.

=== "After"

    ```groovy title="main.nf" linenums="2" hl_lines="2-3"
    ch_samplesheet = Channel.fromPath("./data/samplesheet.csv")
                            .splitCsv(header: true)
                            .view()
    ```

=== "Before"

    ```groovy title="main.nf" linenums="2"
    ch_samplesheet = Channel.fromPath("./data/samplesheet.csv")
    ```

We can use the [`splitCsv` operator](https://www.nextflow.io/docs/latest/operator.html#splitcsv) to split the samplesheet into a channel of maps, where each map represents a row from the CSV file.

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

We can see that each row from the CSV file has been converted into a map with keys matching the header row. A map is a key-value data structure similar to dictionaries in Python, objects in JavaScript, or hashes in Ruby.

Each map contains:

- `id`: an ID given to the sample
- `character`: a character name, that we will use later to draw different creatures
- `data`: paths to `.txt` files that contain phrases in different languages

This format makes it easy to access specific fields from each sample. For example, we could access the sample ID with `id` or the txt file path with `data`. The output above shows each row from the CSV file converted into a map with keys matching the header row. Now that we've successfully read in the samplesheet and have access to the data in each row, we can begin implementing our pipeline logic.

### 1.2 Separate meta data and data

In the samplesheet, we have both the input files and data about the input files (`id`, `character`), the meta data. As we progress through the workflow, we generate more meta data about each sample. To avoid having to keep track of how many fields we have at any point in time and making the input of our processes more robust, we can combine the meta information into its own key-value paired map:

=== "After"

    ```groovy title="main.nf" linenums="3" hl_lines="5-8"
    ch_samplesheet = Channel.fromPath("./data/samplesheet.csv")
                            .splitCsv(header: true)
                            .map { row ->
                              [ [id:row.id, character:row.character], row.recording ]
                            }
                            .view()
    ```

=== "Before"

    ```groovy title="main.nf" linenums="3"
    ch_samplesheet = Channel.fromPath("./data/samplesheet.csv")
                            .splitCsv(header: true)
                            .view()
    ```

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

We have successfully separate our meta data into its own map to keep it next to the file data. Each of our channel elements now has the shape `[ metamap,  file ]`.

### Takeaway

In this section, you've learned:

- **Reading in a samplesheet**: How to read in a samplesheet with `splitCsv`
- **Creating a meta map**: Moving columns with meta information into a separate data structure and keep it next to the input data

---

## 2. Create new meta map keys

### 2.1 Passing the meta map through a process

Now we want to process our samples. These samples are language samples, but we don't know what language they are in. Let's add a process definition before the `workflow` that can identify the language in each file:

=== "After"

```groovy title="main.nf" linenums="1" hl_lines="1-19"
/*
 * Use langid to predict the language of each input file
 */
process IDENTIFY_LANGUAGE {

    container 'community.wave.seqera.io/library/pip_langid:b2269f456a5629ff'

    input:
    tuple val(meta), path(greeting)

    output:
    tuple val(meta), stdout

    script:
    """
    langid < ${greeting} -l en,de,fr,es,it | sed -E "s/.*\\('([a-z]+)'.*/\\1/" | tr -d '\\n'
    """
}

```

=== "Before"

    ```groovy title="main.nf" linenums="1"
    workflow  {
        ch_samplesheet = Channel.fromPath("./data/samplesheet.csv")
                    .splitCsv(header: true)
                    .map { row ->
                        [ [id:row.id, character:row.character], row.recording ]
                    }
                    .view()
    }
    ```

The tool [langid](https://github.com/saffsd/langid.py) is a language identification tool. It is pre-trained on a set of languages. For a given phrase, it prints a language guess and a probability score for each guess to the console. In the `script` section, we are removing the probability score, clean up the string by removing a newline character and return the language guess. Since it is printed directly to the console, we are using Nextflow's [`stdout` output qualifier](https://www.nextflow.io/docs/latest/process.html#outputs), passing the string on as output.

Let's include the process, run, and view it:

=== "After"

```groovy title="main.nf" linenums="20" hl_lines="27-29"
workflow  {

      ch_samplesheet = Channel.fromPath("./data/samplesheet.csv")
                  .splitCsv(header: true)
                  .map { row ->
                      [ [id:row.id, character:row.character], row.recording ]
                  }

      ch_prediction = IDENTIFY_LANGUAGE(ch_samplesheet)
      ch_prediction.view()

  }

```

=== "Before"

    ```groovy title="main.nf" linenums="20"
    workflow  {

        ch_samplesheet = Channel.fromPath("./data/samplesheet.csv")
                    .splitCsv(header: true)
                    .map { row ->
                        [ [id:row.id, character:row.character], row.recording ]
                    }
                    .view()
    }
    ```

```bash title="Identify languages"
nextflow run main.nf
```

```console title="Identify languages"
 N E X T F L O W   ~  version 24.10.4

Launching `main.nf` [voluminous_mcnulty] DSL2 - revision: f9bcfebabb

executor >  local (7)
[2c/888abb] IDENTIFY_LANGUAGE (7) [100%] 7 of 7 ✔
[[id:sampleA, character:squirrel], fr]
[[id:sampleB, character:tux], de]
[[id:sampleC, character:sheep], de]
[[id:sampleD, character:turkey], en]
[[id:sampleE, character:stegosaurus], es]
[[id:sampleF, character:moose], fr]
[[id:sampleG, character:turtle], it]
```

Neat, for each of our samples, we now have a language predicted. You may have noticed something else: we kept the meta data of our samples and associated it with our new piece of information. We achieved this by adding the `meta` map to the output tuple in the process:

```groovy title="main.nf" linenums="12"
output:
      tuple val(meta), stdout
```

This is a useful tool to ensure the sample-specific meta information stays connected with any new information that is generated.

### 2.2 Associate the language prediction with the input file

At the moment, our sample files and their language prediction are separated in two different channels: `ch_samplesheet` and `ch_predictions`. But both channels have the same meta information associated with the interesting data points. We can use the meta map to combine our channels back together.

Nextflow includes many methods for combining channels, but in this case the most appropriate operator is [`join`](https://www.nextflow.io/docs/latest/operator.html#join). If you are familiar with SQL, it acts like the `JOIN` operation, where we specify the key to join on and the type of join to perform.

If we check the [`join`](https://www.nextflow.io/docs/latest/operator.html#join) documentation, we can see that it joins two channels based on a defined item, by default the first item in each tuple. If you don't have the console output still available, let's run the pipeline to check our data structures. If you removed the `view()` operator from the `ch_samplesheet` add it back in:

=== "After"

```groovy title="main.nf" linenums="20" hl_lines="27"
workflow  {

      ch_samplesheet = Channel.fromPath("./data/samplesheet.csv")
                  .splitCsv(header: true)
                  .map { row ->
                      [ [id:row.id, character:row.character], row.recording ]
                  }
                  .view()

      ch_prediction = IDENTIFY_LANGUAGE(ch_samplesheet)
      ch_prediction.view()

  }

```

=== "Before"

    ```groovy title="main.nf" linenums="20"
    workflow  {

        ch_samplesheet = Channel.fromPath("./data/samplesheet.csv")
                    .splitCsv(header: true)
                    .map { row ->
                        [ [id:row.id, character:row.character], row.recording ]
                    }

        ch_prediction = IDENTIFY_LANGUAGE(ch_samplesheet)
        ch_prediction.view()

    }
    ```

```bash title="View samplesheet and prediction channel content"
nextflow run main.nf
```

```console title="View samplesheet and prediction channel content"
 N E X T F L O W   ~  version 24.10.4

Launching `main.nf` [trusting_blackwell] DSL2 - revision: de90745ea4

executor >  local (7)
[d6/0f2efd] IDENTIFY_LANGUAGE (7) [100%] 7 of 7 ✔
[[id:sampleA, character:squirrel], /workspaces/training/side-quests/metadata/data/bonjour.txt]
[[id:sampleB, character:tux], /workspaces/training/side-quests/metadata/data/guten_tag.txt]
[[id:sampleC, character:sheep], /workspaces/training/side-quests/metadata/data/hallo.txt]
[[id:sampleD, character:turkey], /workspaces/training/side-quests/metadata/data/hello.txt]
[[id:sampleE, character:stegosaurus], /workspaces/training/side-quests/metadata/data/hola.txt]
[[id:sampleF, character:moose], /workspaces/training/side-quests/metadata/data/salut.txt]
[[id:sampleG, character:turtle], /workspaces/training/side-quests/metadata/data/ciao.txt]
[[id:sampleB, character:tux], de]
[[id:sampleA, character:squirrel], fr]
[[id:sampleC, character:sheep], de]
[[id:sampleD, character:turkey], en]
[[id:sampleE, character:stegosaurus], es]
[[id:sampleF, character:moose], fr]
[[id:sampleG, character:turtle], it]
```

We can see that the meta map is the first element in each map and the map is the same for both channels. We can simply use the `join` operator to combine the two channels:

=== "After"

```groovy title="main.nf" linenums="20" hl_lines="27"
workflow  {

  ch_samplesheet = Channel.fromPath("./data/samplesheet.csv")
                  .splitCsv(header: true)
                  .map { row ->
                      [ [id:row.id, character:row.character], row.recording ]
                  }

  ch_prediction = IDENTIFY_LANGUAGE(ch_samplesheet)

  ch_languages = ch_samplesheet.join(ch_prediction)
                               .view()
}

```

=== "Before"

    ```groovy title="main.nf" linenums="20"
    workflow  {

        ch_samplesheet = Channel.fromPath("./data/samplesheet.csv")
                    .splitCsv(header: true)
                    .map { row ->
                        [ [id:row.id, character:row.character], row.recording ]
                    }
                    .view()

        ch_prediction = IDENTIFY_LANGUAGE(ch_samplesheet)
        ch_prediction.view()

    }
    ```

```bash title="View joined channel"
nextflow run main.nf
```

```console title="View joined channel"
[[id:sampleA, character:squirrel], /workspaces/training/side-quests/metadata/data/bonjour.txt, fr]
[[id:sampleB, character:tux], /workspaces/training/side-quests/metadata/data/guten_tag.txt, de]
[[id:sampleC, character:sheep], /workspaces/training/side-quests/metadata/data/hallo.txt, de]
[[id:sampleD, character:turkey], /workspaces/training/side-quests/metadata/data/hello.txt, en]
[[id:sampleE, character:stegosaurus], /workspaces/training/side-quests/metadata/data/hola.txt, es]
[[id:sampleF, character:moose], /workspaces/training/side-quests/metadata/data/salut.txt, fr]
[[id:sampleG, character:turtle], /workspaces/training/side-quests/metadata/data/ciao.txt, it]
```

It is becoming a bit hard to see, but if you look all the way on the right side, you can see that now each of our language predictions is associated with our input files.

!!! warning

    The `join` operator will discard any un-matched tuples. In this example, we made sure all samples were matched for tumor and normal but if this is not true you must use the parameter `remainder: true` to keep the unmatched tuples. Check the [documentation](https://www.nextflow.io/docs/latest/operator.html#join) for more details.

### 2.3 Add the language prediction to the meta map

Given that this is more data about the files, let's add it to our meta map. We can use the [`map` operator](https://www.nextflow.io/docs/latest/operator.html#map) to create a new key `lang` and set the value to prediction:

=== "After"

```groovy title="main.nf" linenums="20" hl_lines="31-33"
workflow  {

  ch_samplesheet = Channel.fromPath("./data/samplesheet.csv")
                  .splitCsv(header: true)
                  .map { row ->
                      [ [id:row.id, character:row.character], row.recording ]
                  }

  ch_prediction = IDENTIFY_LANGUAGE(ch_samplesheet)

  ch_languages = ch_samplesheet.join(ch_prediction)
                              .map { meta, file, lang ->
                                  [ meta + [lang:lang], file ]
                              }
                              .view()

}
```

=== "Before"

    ```groovy title="main.nf" linenums="20"
    workflow  {

      ch_samplesheet = Channel.fromPath("./data/samplesheet.csv")
                      .splitCsv(header: true)
                      .map { row ->
                          [ [id:row.id, character:row.character], row.recording ]
                      }

      ch_prediction = IDENTIFY_LANGUAGE(ch_samplesheet)

      ch_languages = ch_samplesheet.join(ch_prediction)
                                   .view()
    }
    ```

```bash title="View new meta map key"
nextflow run main.nf -resume
```

```console title="View new meta map key"

 N E X T F L O W   ~  version 24.10.4

Launching `main.nf` [cheeky_fermat] DSL2 - revision: d096281ee4

[da/652cc6] IDENTIFY_LANGUAGE (7) [100%] 7 of 7, cached: 7 ✔
[[id:sampleA, character:squirrel, lang:fr], /workspaces/training/side-quests/metadata/data/bonjour.txt]
[[id:sampleB, character:tux, lang:de], /workspaces/training/side-quests/metadata/data/guten_tag.txt]
[[id:sampleC, character:sheep, lang:de], /workspaces/training/side-quests/metadata/data/hallo.txt]
[[id:sampleD, character:turkey, lang:en], /workspaces/training/side-quests/metadata/data/hello.txt]
[[id:sampleE, character:stegosaurus, lang:es], /workspaces/training/side-quests/metadata/data/hola.txt]
[[id:sampleF, character:moose, lang:fr], /workspaces/training/side-quests/metadata/data/salut.txt]
[[id:sampleG, character:turtle, lang:it], /workspaces/training/side-quests/metadata/data/ciao.txt]
```

Nice, we expanded our meta map with new information we gathered in the pipeline.

### 2.4 Assign a language group using a ternary operator

Alright, now that we have our language predictions, let's use the information to assign them into new groups. In our example data, we have provided data sets that belong either to `germanic` (either English or German) or `romanic` (French, Spanish, Italian) languages.

We can use the `map` operator and an [ternary operator](https://groovy-lang.org/operators.html#_ternary_operator) to assign either group. The ternary operator, is a short cut to an if/else clause. It says:

```
variable = <condition> ? 'Value' : 'Default'
```

and is the same as:

```
if (<condition>){
  variable = 'Value'
} else {
  variable = 'Default'
}
```

=== "After"

```groovy title="main.nf" linenums="20" hl_lines="34-37"
workflow  {

  ch_samplesheet = Channel.fromPath("./data/samplesheet.csv")
                  .splitCsv(header: true)
                  .map { row ->
                      [ [id:row.id, character:row.character], row.recording ]
                  }

  ch_prediction = IDENTIFY_LANGUAGE(ch_samplesheet)

  ch_languages = ch_samplesheet.join(ch_prediction)
                              .map { meta, file, lang ->
                                  [ meta + [lang:lang], file ]
                              }
                              .map{ meta, file ->
                                  def lang_group = (meta.lang.equals('de') || meta.lang.equals('en')) ? 'germanic' : 'romanic'
                                  [ meta + [lang_group:lang_group], file ]
                              }
                              .view()

}

```

=== "Before"

    ```groovy title="main.nf" linenums="20"
    workflow  {

      ch_samplesheet = Channel.fromPath("./data/samplesheet.csv")
                      .splitCsv(header: true)
                      .map { row ->
                          [ [id:row.id, character:row.character], row.recording ]
                      }

      ch_prediction = IDENTIFY_LANGUAGE(ch_samplesheet)

      ch_languages = ch_samplesheet.join(ch_prediction)
                                  .map { meta, file, lang ->
                                      [ meta + [lang:lang], file ]
                                  }
                                  .view()

    }
    ```

Let's rerun it

```bash title="View language groups"
nextflow run main.nf -resume
```

```console title="View language groups"
 N E X T F L O W   ~  version 24.10.4

Launching `main.nf` [wise_almeida] DSL2 - revision: 46778c3cd0

[da/652cc6] IDENTIFY_LANGUAGE (7) [100%] 7 of 7, cached: 7 ✔
[[id:sampleA, character:squirrel, lang:fr, lang_group:romanic], /workspaces/training/side-quests/metadata/data/bonjour.txt]
[[id:sampleB, character:tux, lang:de, lang_group:germanic], /workspaces/training/side-quests/metadata/data/guten_tag.txt]
[[id:sampleC, character:sheep, lang:de, lang_group:germanic], /workspaces/training/side-quests/metadata/data/hallo.txt]
[[id:sampleD, character:turkey, lang:en, lang_group:germanic], /workspaces/training/side-quests/metadata/data/hello.txt]
[[id:sampleE, character:stegosaurus, lang:es, lang_group:romanic], /workspaces/training/side-quests/metadata/data/hola.txt]
[[id:sampleF, character:moose, lang:fr, lang_group:romanic], /workspaces/training/side-quests/metadata/data/salut.txt]
[[id:sampleG, character:turtle, lang:it, lang_group:romanic], /workspaces/training/side-quests/metadata/data/ciao.txt]
```

### Takeaway

In this section, you've learned:

- **Merging on meta maps**: You used `join` to combine two channels based on their meta maps
- **Creating custom keys**: You created two new keys in your meta map. One based on a computed value from a process, and one based on a condition you set in the `map` operator.

Both of these allow you to associated new and existing meta data with files as you progress through your pipeline.

---

## 3. Filter data based on meta map values

We can use the [`filter` operator](https://www.nextflow.io/docs/latest/operator.html#filter) to filter the data based on a condition. Let's say we only want to process romanic language samples firther. We can do this by filtering the data based on the `lang_group` field. Let's create a new channel that only contains romanic languages and `view` it:

=== "After"

```groovy title="main.nf" linenums="20" hl_lines="38-42"
workflow  {

  ch_samplesheet = Channel.fromPath("./data/samplesheet.csv")
                  .splitCsv(header: true)
                  .map { row ->
                      [ [id:row.id, character:row.character], row.recording ]
                  }

  ch_prediction = IDENTIFY_LANGUAGE(ch_samplesheet)

  ch_languages = ch_samplesheet.join(ch_prediction)
                              .map { meta, file, lang ->
                                  [ meta + [lang:lang], file ]
                              }
                              .map{ meta, file ->
                                  def lang_group = (meta.lang.equals('de') || meta.lang.equals('en')) ? 'germanic' : 'romanic'
                                  [ meta + [lang_group:lang_group], file ]
                              }

romanic_languages = ch_language_groups.filter { meta, file ->
                                        meta.lang_group == 'romanic'
                                      }
                                      .view()
}

```

=== "Before"

    ```groovy title="main.nf" linenums="20"
    workflow  {

      ch_samplesheet = Channel.fromPath("./data/samplesheet.csv")
                      .splitCsv(header: true)
                      .map { row ->
                          [ [id:row.id, character:row.character], row.recording ]
                      }

      ch_prediction = IDENTIFY_LANGUAGE(ch_samplesheet)

      ch_languages = ch_samplesheet.join(ch_prediction)
                                  .map { meta, file, lang ->
                                      [ meta + [lang:lang], file ]
                                  }
                                  .map{ meta, file ->
                                      def lang_group = (meta.lang.equals('de') || meta.lang.equals('en')) ? 'germanic' : 'romanic'
                                      [ meta + [lang_group:lang_group], file ]
                                  }
                                  .view()

    }
    ```

Let's rerun it

```bash title="View romanic samples"
nextflow run main.nf -resume
```

```console title="View romanic samples"
 N E X T F L O W   ~  version 24.10.4

Launching `main.nf` [drunk_brattain] DSL2 - revision: 453fdd4e91

[da/652cc6] IDENTIFY_LANGUAGE (7) [100%] 7 of 7, cached: 7 ✔
[[id:sampleA, character:squirrel, lang:fr, lang_group:romanic], /workspaces/training/side-quests/metadata/data/bonjour.txt]
[[id:sampleE, character:stegosaurus, lang:es, lang_group:romanic], /workspaces/training/side-quests/metadata/data/hola.txt]
[[id:sampleF, character:moose, lang:fr, lang_group:romanic], /workspaces/training/side-quests/metadata/data/salut.txt]
[[id:sampleG, character:turtle, lang:it, lang_group:romanic], /workspaces/training/side-quests/metadata/data/ciao.txt]
```

We have successfully filtered the data to only include romanic samples. Let's recap how this works. The `filter` operator takes a closure that is applied to each element in the channel. If the closure returns `true`, the element is included in the output channel. If the closure returns `false`, the element is excluded from the output channel.

In this case, we want to keep only the samples where `meta.lang_group == 'romanic'`. In the closure, we first know that our channel elements are all of shape `[meta, file]` and we can then access the individual keys of the meta map. We then check if `meta.lang_group` is equal to `'romanic'`. If it is, the sample is included in the output channel. If it is not, the sample is excluded from the output channel.

```groovy title="main.nf" linenums="4"
.filter { meta,file -> meta.lang_group == 'romanic' }
```

### Takeaway

In this section, you've learned:

- How to filter data with `filter`

we now have only the romanic language samples left and can process those further. Next we want to make characters say the phrases.

---

## 4. Customize a process with meta map

Let's let the characters say the phrases that we have passed in. In the [hello-nextflow training](../hello_nextflow/05_hello_containers.md), you already encountered the `cowpy` package, a python implementation of a tool called `cowsay` that generates ASCII art to display arbitrary text inputs in a fun way.

We will re-use a process from there.

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

### 4.1 Add a custom publishing location

Let's run our romanic languages through `COWPY` and remove our `view` statement:

=== "After"

```groovy title="main.nf" linenums="40" hl_lines="61-63"
workflow  {

  ch_samplesheet = Channel.fromPath("./data/samplesheet.csv")
                  .splitCsv(header: true)
                  .map { row ->
                      [ [id:row.id, character:row.character], row.recording ]
                  }

  ch_prediction = IDENTIFY_LANGUAGE(ch_samplesheet)

  ch_languages = ch_samplesheet.join(ch_prediction)
                              .map { meta, file, lang ->
                                  [ meta + [lang:lang], file ]
                              }
                              .map{ meta, file ->
                                  def lang_group = (meta.lang.equals('de') || meta.lang.equals('en')) ? 'germanic' : 'romanic'
                                  [ meta + [lang_group:lang_group], file ]
                              }

  romanic_languages = ch_languages.filter { meta, file ->
                                      meta.lang_group == 'romanic'
                                  }

  COWPY(romanic_languages)

}
```

=== "Before"

    ```groovy title="main.nf" linenums="20"
    workflow  {

    ch_samplesheet = Channel.fromPath("./data/samplesheet.csv")
                    .splitCsv(header: true)
                    .map { row ->
                        [ [id:row.id, character:row.character], row.recording ]
                    }

    ch_prediction = IDENTIFY_LANGUAGE(ch_samplesheet)

    ch_languages = ch_samplesheet.join(ch_prediction)
                                .map { meta, file, lang ->
                                    [ meta + [lang:lang], file ]
                                }
                                .map{ meta, file ->
                                    def lang_group = (meta.lang.equals('de') || meta.lang.equals('en')) ? 'germanic' : 'romanic'
                                    [ meta + [lang_group:lang_group], file ]
                                }

    romanic_languages = ch_languages.filter { meta, file ->
                                        meta.lang_group == 'romanic'
                                    }.view()
    }
    ```

We are still missing a publishing location. Given we have been trying to figure out what languages our samples were in, let's group the samples by language in the output directory. Earlier, we added the predicted language to the `meta` map. We can access this `key` in the process and use it in the `publishDir` directive:

=== "After"

```groovy title="main.nf" linenums="23" hl_lines="25"
/*
 * Generate ASCII art with cowpy
*/
process COWPY {

    publishDir "results/${meta.lang}", mode: 'copy'

    container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

```

=== "Before"

    ```groovy title="main.nf" linenums="23"
    process COWPY {

      container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

    ```

Let's run this:

```bash title="Use cowpy"
nextflow run main.nf -resume
```

You should now see a new folder called `results`:

```console title="Results folder"
results/
├── es
│   └── cowpy-hola.txt
├── fr
│   ├── cowpy-bonjour.txt
│   └── cowpy-salut.txt
└── it
    └── cowpy-ciao.txt
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

### 4.2 Customize the character

In our samplesheet, we have another column: `character`. To tailor the tool parameters per sample, we can also access information from the `meta` map in the script section. This is really useful in cases were a tool should have different parameters for each sample.

Let's customize the characters by changing the `cowpy` command:

=== "After"

```groovy title="main.nf" linenums="35" hl_lines="37"
script:
"""
cat $input_file | cowpy -c ${meta.character} > cowpy-${input_file}
"""
```

=== "Before"

    ```groovy title="main.nf" linenums="23"
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

```console title="fr/cowpy-salut.txt"
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

In this section, you've learned:

- **Tweaking directives using meta values**
- **Tweaking script section based on meta values**

---

## Summary
