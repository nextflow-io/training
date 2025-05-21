# Metadata

introduction:

what is meta data
why is it important
sample specific, not something that is the same for all samples

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

  ``` groovy title="main.nf" linenums="1" hl_lines="1-19"
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

  ``` groovy title="main.nf" linenums="20" hl_lines="27-29"
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

  ``` groovy title="main.nf" linenums="20" hl_lines="27"
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

  ``` groovy title="main.nf" linenums="20" hl_lines="27"
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

  ``` groovy title="main.nf" linenums="20" hl_lines="31-33"
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

  ``` groovy title="main.nf" linenums="20" hl_lines="34-37"
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

## 3. Filter data with certain values in the meta map


Repeat a grouping example, using this new value as grouping or branching decider (Objective: use the meta map to decide on workflow paths)
Run a module with the groups

### Takeaway

In this section, you've learned:

- **Extracting an arbitray value to group or filter on**

---

## 4. Publishing location based on meta map value

tweak the publishing directory based on a field in the meta map (Objective: use the meta map in the module)

### Takeaway

In this section, you've learned:

- **Tweaking directives using meta values**

---

## 5. Tweak tool arguments based on meta map value

tweak the publishing directory based on a field in the meta map (Objective: use the meta map in the module)

### Takeaway

In this section, you've learned:

- **Tweaking script section based on meta values**

---

## Summary
