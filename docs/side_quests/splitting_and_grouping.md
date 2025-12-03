# Splitting and Grouping

Nextflow provides powerful tools for working with data flexibly. A key capability is splitting data into different streams and then grouping related items back together. This is especially valuable in bioinformatics workflows where you need to process different types of samples separately before combining results for analysis.

Think of it like sorting mail: you separate letters by destination, process each pile differently, then recombine items going to the same person. Nextflow uses special operators to accomplish this with scientific data. This approach is also commonly known as the **scatter/gather** pattern in distributed computing and bioinformatics workflows.

Nextflow's channel system is at the heart of this flexibility. Channels connect different parts of your workflow, allowing data to flow through your analysis. You can create multiple channels from a single data source, process each channel differently, and then merge channels back together when needed. This approach lets you design workflows that naturally mirror the branching and converging paths of complex bioinformatics analyses.

### Learning goals

In this side quest, you'll learn to split and group data using Nextflow's channel operators.
We'll start with a CSV file containing sample information and associated data files, then manipulate and reorganize this data.

By the end of this side quest, you'll be able to separate and combine data streams effectively, using the following techniques:

- Read data from files using `splitCsv`
- Filter and transform data with `filter` and `map`
- Combine related data using `join` and `groupTuple`
- Create data combinations with `combine` for parallel processing
- Optimize data structure using `subMap` and deduplication strategies
- Build reusable functions with named closures to help you manipulate channel structures

These skills will help you build workflows that can handle multiple input files and different types of data efficiently, while maintaining clean, maintainable code structure.

### Prerequisites

Before taking on this side quest, you should:

- Have completed the [Hello Nextflow](../hello_nextflow/README.md) tutorial or equivalent beginner's course.
- Be comfortable using basic Nextflow concepts and mechanisms (processes, channels, operators, working with files, meta data)

**Optional:** We recommend completing the [Metadata in workflows](./metadata.md) side quest first.
That covers the fundamentals of reading CSV files with `splitCsv` and creating meta maps, which we'll use heavily here.

---

## 0. Get started

#### Open the training codespace

If you haven't yet done so, make sure to open the training environment as described in the [Environment Setup](../envsetup/index.md).

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Move into the project directory

Let's move into the directory where the files for this tutorial are located.

```bash
cd side-quests/splitting_and_grouping
```

You can set VSCode to focus on this directory:

```bash
code .
```

#### Review the materials

You'll find a main workflow file and a `data` directory containing a samplesheet named `samplesheet.csv`.

```console title="Directory contents"
.
├── data
│   └── samplesheet.csv
└── main.nf
```

The samplesheet contains information about samples from different patients, including the patient ID, sample repeat number, type (normal or tumor), and paths to hypothetical data files (which don't actually exist, but we will pretend they do).

```console title="samplesheet.csv"
id,repeat,type,bam
patientA,1,normal,patientA_rep1_normal.bam
patientA,1,tumor,patientA_rep1_tumor.bam
patientA,2,normal,patientA_rep2_normal.bam
patientA,2,tumor,patientA_rep2_tumor.bam
patientB,1,normal,patientB_rep1_normal.bam
patientB,1,tumor,patientB_rep1_tumor.bam
patientC,1,normal,patientC_rep1_normal.bam
patientC,1,tumor,patientC_rep1_tumor.bam
```

This samplesheet lists eight samples from three patients (A, B, C).

For each patient, we have samples that are of type `tumor` (typically originating from tumor biopsies) or `normal` (taken from healthy tissue or blood).
If you're not familiar with cancer analysis, just know that this corresponds to an experimental model that uses paired tumor/normal samples to perform contrastive analyses.

For patient A specifically, we have two sets of technical replicates (repeats).

!!! note

    Don't worry if you're not familiar with this experimental design, it's not critical for understanding this tutorial.

#### Review the assignment

Your challenge is to write a Nextflow workflow that will group and split the samples based on the associated metadata.

<!-- TODO: give a bit more details, similar to how it's done in the Metadata side quest -->

#### Readiness checklist

Think you're ready to dive in?

- [ ] I understand the goal of this course and its prerequisites
- [ ] My codespace is up and running
- [ ] I've set my working directory appropriately
- [ ] I understand the assignment

If you can check all the boxes, you're good to go.

---

## 1. Read in sample data

### 1.1. Read in sample data with splitCsv and create meta maps

Let's start by reading in the sample data with `splitCsv` and organizing it into the meta map pattern. In the `main.nf`, you'll see that we've already started the workflow.

```groovy title="main.nf" linenums="1" hl_lines="2"
workflow {
    ch_samplesheet = channel.fromPath("./data/samplesheet.csv")
}
```

!!! note

    Throughout this tutorial, we'll use the `ch_` prefix for all channel variables to clearly indicate they are Nextflow channels.

If you completed the [Metadata in workflows](./metadata.md) side quest, you'll recognize this pattern. We'll use `splitCsv` to read the CSV and immediately structure the data with a meta map to separate metadata from file paths.

!!!Note

    We'll encounter two different concepts called `map` in this training:

    - **Data structure**: The Groovy map (equivalent to dictionaries/hashes in other languages) that stores key-value pairs
    - **Channel operator**: The `.map()` operator that transforms items in a channel

    We'll clarify which one we mean in context, but this distinction is important to understand when working with Nextflow.

Apply these changes to `main.nf`:

=== "After"

    ```groovy title="main.nf" linenums="2" hl_lines="2-6"
        ch_samples = channel.fromPath("./data/samplesheet.csv")
            .splitCsv(header: true)
            .map{ row ->
              [[id:row.id, repeat:row.repeat, type:row.type], row.bam]
            }
            .view()
    ```

=== "Before"

    ```groovy title="main.nf" linenums="2" hl_lines="1"
        ch_samplesheet = channel.fromPath("./data/samplesheet.csv")
    ```

This combines the `splitCsv` operation (reading the CSV with headers) and the `map` operation (structuring data as `[meta, file]` tuples) in one step. Apply that change and run the pipeline:

```bash title="Verify the data structure"
nextflow run main.nf
```

```console title="Sample data with meta maps"
 N E X T F L O W   ~  version 25.04.3

Launching `main.nf` [deadly_mercator] DSL2 - revision: bd6b0224e9

[[id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam]
[[id:patientA, repeat:1, type:tumor], patientA_rep1_tumor.bam]
[[id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam]
[[id:patientA, repeat:2, type:tumor], patientA_rep2_tumor.bam]
[[id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam]
[[id:patientB, repeat:1, type:tumor], patientB_rep1_tumor.bam]
[[id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam]
[[id:patientC, repeat:1, type:tumor], patientC_rep1_tumor.bam]
```

We now have a channel where each item is a `[meta, file]` tuple - metadata separated from file paths. This structure allows us to split and group our workload based on metadata fields.

---

## 2. Filter and transform data

### 2.1. Filter data with `filter`

We can use the [`filter` operator](https://www.nextflow.io/docs/latest/operator.html#filter) to filter the data based on a condition. Let's say we only want to process normal samples. We can do this by filtering the data based on the `type` field. Let's insert this before the `view` operator.

=== "After"

    ```groovy title="main.nf" linenums="2" hl_lines="6"
        ch_samples = channel.fromPath("./data/samplesheet.csv")
            .splitCsv(header: true)
            .map{ row ->
              [[id:row.id, repeat:row.repeat, type:row.type], row.bam]
            }
            .filter { meta, file -> meta.type == 'normal' }
            .view()
    ```

=== "Before"

    ```groovy title="main.nf" linenums="2"
        ch_samples = channel.fromPath("./data/samplesheet.csv")
            .splitCsv(header: true)
            .map{ row ->
              [[id:row.id, repeat:row.repeat, type:row.type], row.bam]
            }
            .view()
    ```

Run the workflow again to see the filtered result:

```bash title="Test the filter operation"
nextflow run main.nf
```

```console title="Filtered normal samples"
 N E X T F L O W   ~  version 25.04.3

Launching `main.nf` [admiring_brown] DSL2 - revision: 194d61704d

[[id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam]
[[id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam]
[[id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam]
[[id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam]
```

We have successfully filtered the data to only include normal samples. Let's recap how this works.

The `filter` operator takes a closure that is applied to each element in the channel. If the closure returns `true`, the element is included; if it returns `false`, the element is excluded.

In our case, we want to keep only samples where `meta.type == 'normal'`. The closure uses the tuple `meta,file` to refer to each sample, accesses the sample type with `meta.type`, and checks if it equals `'normal'`.

This is accomplished with the single closure we introduced above:

```groovy title="main.nf" linenums="7"
    .filter { meta, file -> meta.type == 'normal' }
```

### 2.2. Create separate filtered channels

Currently we're applying the filter to the channel created directly from the CSV, but we want to filter this in more ways than one, so let's re-write the logic to create a separate filtered channel for normal samples:

=== "After"

    ```groovy title="main.nf" linenums="2" hl_lines="6 8"
        ch_samples = channel.fromPath("./data/samplesheet.csv")
            .splitCsv(header: true)
            .map{ row ->
                [[id:row.id, repeat:row.repeat, type:row.type], row.bam]
            }
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
        ch_normal_samples
            .view()
    ```

=== "Before"

    ```groovy title="main.nf" linenums="2"
        ch_samples = channel.fromPath("./data/samplesheet.csv")
            .splitCsv(header: true)
            .map{ row ->
              [[id:row.id, repeat:row.repeat, type:row.type], row.bam]
            }
            .filter { meta, file -> meta.type == 'normal' }
            .view()
    ```

Once again, run the pipeline to see the results:

```bash title="Test separate channel creation"
nextflow run main.nf
```

```console title="Filtered normal samples"
 N E X T F L O W   ~  version 25.04.3

Launching `main.nf` [trusting_poisson] DSL2 - revision: 639186ee74

[[id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam]
[[id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam]
[[id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam]
[[id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam]
```

We've successfully filtered the data and created a separate channel for normal samples. Let's create a filtered channel for the tumor samples as well:

=== "After"

    ```groovy title="main.nf" linenums="7" hl_lines="3-8"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
        ch_normal_samples
            .view{'Normal sample: ' + it}
        ch_tumor_samples
            .view{'Tumor sample: ' + it}
    ```

=== "Before"

    ```groovy title="main.nf" linenums="7" hl_lines="3 4"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
        ch_normal_samples
            .view()
    ```

```bash title="Test filtering both sample types"
nextflow run main.nf
```

```console title="Normal and tumor samples"
 N E X T F L O W   ~  version 25.04.3

Launching `main.nf` [maniac_boltzmann] DSL2 - revision: 3636b6576b

Tumor sample: [[id:patientA, repeat:1, type:tumor], patientA_rep1_tumor.bam]
Tumor sample: [[id:patientA, repeat:2, type:tumor], patientA_rep2_tumor.bam]
Normal sample: [[id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam]
Normal sample: [[id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam]
Normal sample: [[id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam]
Normal sample: [[id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam]
Tumor sample: [[id:patientB, repeat:1, type:tumor], patientB_rep1_tumor.bam]
Tumor sample: [[id:patientC, repeat:1, type:tumor], patientC_rep1_tumor.bam]
```

We've separated out the normal and tumor samples into two different channels, and used a closure supplied to `view()` to label them differently in the output: `ch_tumor_samples.view{'Tumor sample: ' + it}`.

### Takeaway

In this section, you've learned:

- **Filtering data**: How to filter data with `filter`
- **Splitting data**: How to split data into different channels based on a condition
- **Viewing data**: How to use `view` to print the data and label output from different channels

We've now separated out the normal and tumor samples into two different channels. Next, we'll join the normal and tumor samples on the `id` field.

---

## 3. Joining channels by identifiers

In the previous section, we separated out the normal and tumor samples into two different channels. These could be processed independently using specific processes or workflows based on their type. But what happens when we want to compare the normal and tumor samples from the same patient? At this point, we need to join them back together making sure to match the samples based on their `id` field.

Nextflow includes many methods for combining channels, but in this case the most appropriate operator is [`join`](https://www.nextflow.io/docs/latest/operator.html#join). If you are familiar with SQL, it acts like the `JOIN` operation, where we specify the key to join on and the type of join to perform.

### 3.1. Use `map` and `join` to combine based on patient ID

If we check the [`join`](https://www.nextflow.io/docs/latest/operator.html#join) documentation, we can see that by default it joins two channels based on the first item in each tuple. If you don't have the console output still available, let's run the pipeline to check our data structure and see how we need to modify it to join on the `id` field.

```bash title="Check current data structure"
nextflow run main.nf
```

```console title="Normal and tumor samples"
 N E X T F L O W   ~  version 25.04.3

Launching `main.nf` [maniac_boltzmann] DSL2 - revision: 3636b6576b

Tumor sample: [[id:patientA, repeat:1, type:tumor], patientA_rep1_tumor.bam]
Tumor sample: [[id:patientA, repeat:2, type:tumor], patientA_rep2_tumor.bam]
Normal sample: [[id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam]
Normal sample: [[id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam]
Normal sample: [[id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam]
Normal sample: [[id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam]
Tumor sample: [[id:patientB, repeat:1, type:tumor], patientB_rep1_tumor.bam]
Tumor sample: [[id:patientC, repeat:1, type:tumor], patientC_rep1_tumor.bam]
```

We can see that the `id` field is the first element in each meta map. For `join` to work, we should isolate the `id` field in each tuple. After that, we can simply use the `join` operator to combine the two channels.

To isolate the `id` field, we can use the [`map` operator](https://www.nextflow.io/docs/latest/operator.html#map) to create a new tuple with the `id` field as the first element.

=== "After"

    ```groovy title="main.nf" linenums="7" hl_lines="3 6"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
            .map { meta, file -> [meta.id, meta, file] }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
            .map { meta, file -> [meta.id, meta, file] }
        ch_normal_samples
            .view{'Normal sample: ' + it}
        ch_tumor_samples
            .view{'Tumor sample: ' + it}
    ```

=== "Before"

    ```groovy title="main.nf" linenums="7"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
        ch_normal_samples
            .view{'Normal sample: ' + it}
        ch_tumor_samples
            .view{'Tumor sample: ' + it}
    ```

```bash title="Test the map transformation"
nextflow run main.nf
```

```console title="Samples with ID as first element"
 N E X T F L O W   ~  version 25.04.3

Launching `main.nf` [mad_lagrange] DSL2 - revision: 9940b3f23d

Tumor sample: [patientA, [id:patientA, repeat:1, type:tumor], patientA_rep1_tumor.bam]
Tumor sample: [patientA, [id:patientA, repeat:2, type:tumor], patientA_rep2_tumor.bam]
Normal sample: [patientA, [id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam]
Normal sample: [patientA, [id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam]
Tumor sample: [patientB, [id:patientB, repeat:1, type:tumor], patientB_rep1_tumor.bam]
Tumor sample: [patientC, [id:patientC, repeat:1, type:tumor], patientC_rep1_tumor.bam]
Normal sample: [patientB, [id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam]
Normal sample: [patientC, [id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam]
```

It might be subtle, but you should be able to see the first element in each tuple is the `id` field. Now we can use the `join` operator to combine the two channels based on the `id` field.

Once again, we will use `view` to print the joined outputs.

=== "After"

    ```groovy title="main.nf" linenums="7" hl_lines="7-9"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
            .map { meta, file -> [meta.id, meta, file] }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
            .map { meta, file -> [meta.id, meta, file] }
        ch_joined_samples = ch_normal_samples
            .join(ch_tumor_samples)
        ch_joined_samples.view()
    ```

=== "Before"

    ```groovy title="main.nf" linenums="7" hl_lines="7-10"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
            .map { meta, file -> [meta.id, meta, file] }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
            .map { meta, file -> [meta.id, meta, file] }
        ch_normal_samples
            .view{'Normal sample: ' + it}
        ch_tumor_samples
            .view{'Tumor sample: ' + it}
    ```

```bash title="Test the join operation"
nextflow run main.nf
```

```console title="Joined normal and tumor samples"
 N E X T F L O W   ~  version 25.04.3

Launching `main.nf` [soggy_wiles] DSL2 - revision: 3bc1979889

[patientA, [id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam, [id:patientA, repeat:1, type:tumor], patientA_rep1_tumor.bam]
[patientA, [id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam, [id:patientA, repeat:2, type:tumor], patientA_rep2_tumor.bam]
[patientB, [id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam, [id:patientB, repeat:1, type:tumor], patientB_rep1_tumor.bam]
[patientC, [id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam, [id:patientC, repeat:1, type:tumor], patientC_rep1_tumor.bam]
```

It's a little hard to tell because it's so wide, but you should be able to see the samples have been joined by the `id` field. Each tuple now has the format:

- `id`: The sample ID
- `normal_meta_map`: The normal sample meta data including type, replicate and path to bam file
- `normal_sample_file`: The normal sample file
- `tumor_meta_map`: The tumor sample meta data including type, replicate and path to bam file
- `tumor_sample`: The tumor sample including type, replicate and path to bam file

!!! warning

    The `join` operator will discard any un-matched tuples. In this example, we made sure all samples were matched for tumor and normal but if this is not true you must use the parameter `remainder: true` to keep the unmatched tuples. Check the [documentation](https://www.nextflow.io/docs/latest/operator.html#join) for more details.

### Takeaway

In this section, you've learned:

- How to use `map` to isolate a field in a tuple
- How to use `join` to combine tuples based on the first field

With this knowledge, we can successfully combine channels based on a shared field. Next, we'll consider the situation where you want to join on multiple fields.

### 3.2. Join on multiple fields

We have 2 replicates for sampleA, but only 1 for sampleB and sampleC. In this case we were able to join them effectively by using the `id` field, but what would happen if they were out of sync? We could mix up the normal and tumor samples from different replicates!

To avoid this, we can join on multiple fields. There are actually multiple ways to achieve this but we are going to focus on creating a new joining key which includes both the sample `id` and `replicate` number.

Let's start by creating a new joining key. We can do this in the same way as before by using the [`map` operator](https://www.nextflow.io/docs/latest/operator.html#map) to create a new tuple with the `id` and `repeat` fields as the first element.

=== "After"

    ```groovy title="main.nf" linenums="7" hl_lines="3 6"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
            .map { meta, file -> [[meta.id, meta.repeat], meta, file] }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
            .map { meta, file -> [[meta.id, meta.repeat], meta, file] }
    ```

=== "Before"

    ```groovy title="main.nf" linenums="7" hl_lines="3 6"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
            .map { meta, file -> [meta.id, meta, file] }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
            .map { meta, file -> [meta.id, meta, file] }
    ```

Now we should see the join is occurring but using both the `id` and `repeat` fields. Run the workflow:

```bash title="Test multi-field joining"
nextflow run main.nf
```

```console title="Samples joined on multiple fields"
 N E X T F L O W   ~  version 25.04.3

Launching `main.nf` [prickly_wing] DSL2 - revision: 3bebf22dee

[[patientA, 1], [id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam, [id:patientA, repeat:1, type:tumor], patientA_rep1_tumor.bam]
[[patientA, 2], [id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam, [id:patientA, repeat:2, type:tumor], patientA_rep2_tumor.bam]
[[patientB, 1], [id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam, [id:patientB, repeat:1, type:tumor], patientB_rep1_tumor.bam]
[[patientC, 1], [id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam, [id:patientC, repeat:1, type:tumor], patientC_rep1_tumor.bam]
```

Note how we have a tuple of two elements (`id` and `repeat` fields) as the first element of each joined result. This demonstrates how complex items can be used as a joining key, enabling fairly intricate matching between samples from the same conditions.

If you want to explore more ways to join on different keys, check out the [join operator documentation](https://www.nextflow.io/docs/latest/operator.html#join) for additional options and examples.

### 3.3. Use subMap to create a new joining key

The previous approach loses the field names from our joining key - the `id` and `repeat` fields become just a list of values. To retain the field names for later access, we can use the [`subMap` method](<https://docs.groovy-lang.org/latest/html/groovy-jdk/java/util/Map.html#subMap(java.util.Collection)>).

The `subMap` method extracts only the specified key-value pairs from a map. Here we'll extract just the `id` and `repeat` fields to create our joining key.

=== "After"

    ```groovy title="main.nf" linenums="7" hl_lines="3 6"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
            .map { meta, file -> [meta.subMap(['id', 'repeat']), meta, file] }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
            .map { meta, file -> [meta.subMap(['id', 'repeat']), meta, file] }
    ```

=== "Before"

    ```groovy title="main.nf" linenums="7" hl_lines="3 6"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
            .map { meta, file -> [[meta.id, meta.repeat], meta, file] }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
            .map { meta, file -> [[meta.id, meta.repeat], meta, file] }
    ```

```bash title="Test subMap joining keys"
nextflow run main.nf
```

```console title="Samples with subMap joining keys"
 N E X T F L O W   ~  version 25.04.3

Launching `main.nf` [reverent_wing] DSL2 - revision: 847016c3b7

[[id:patientA, repeat:1], [id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam, [id:patientA, repeat:1, type:tumor], patientA_rep1_tumor.bam]
[[id:patientA, repeat:2], [id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam, [id:patientA, repeat:2, type:tumor], patientA_rep2_tumor.bam]
[[id:patientB, repeat:1], [id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam, [id:patientB, repeat:1, type:tumor], patientB_rep1_tumor.bam]
[[id:patientC, repeat:1], [id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam, [id:patientC, repeat:1, type:tumor], patientC_rep1_tumor.bam]
```

Now we have a new joining key that not only includes the `id` and `repeat` fields but also retains the field names so we can access them later by name, e.g. `meta.id` and `meta.repeat`.

### 3.4. Use a named closure in map

To avoid duplication and reduce errors, we can use a named closure. A named closure allows us to create a reusable function that we can call in multiple places.

To do so, first we define the closure as a new variable:

=== "After"

    ```groovy title="main.nf" linenums="2" hl_lines="7"
        ch_samples = channel.fromPath("./data/samplesheet.csv")
            .splitCsv(header: true)
            .map{ row ->
                [[id:row.id, repeat:row.repeat, type:row.type], row.bam]
            }

        getSampleIdAndReplicate = { meta, bam -> [ meta.subMap(['id', 'repeat']), meta, file(bam) ] }

        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
    ```

=== "Before"

    ```groovy title="main.nf" linenums="2"
        ch_samples = channel.fromPath("./data/samplesheet.csv")
            .splitCsv(header: true)
            .map{ row ->
                [[id:row.id, repeat:row.repeat, type:row.type], row.bam]
            }
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
    ```

We've defined the map transformation as a named variable that we can reuse. Note that we also convert the file path to a Path object using `file()` so that any process receiving this channel can handle the file correctly (for more information see [Working with files](./working_with_files.md)).

Let's implement the closure in our workflow:

=== "After"

    ```groovy title="main.nf" linenums="10" hl_lines="3 6"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
             .map ( getSampleIdAndReplicate )
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
             .map ( getSampleIdAndReplicate )

    ```

=== "Before"

    ```groovy title="main.nf" linenums="10" hl_lines="3 6"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
            .map { meta, file -> [meta.subMap(['id', 'repeat']), meta, file] }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
            .map { meta, file -> [meta.subMap(['id', 'repeat']), meta, file] }
    ```

!!! note

    The `map` operator has switched from using `{ }` to using `( )` to pass the closure as an argument. This is because the `map` operator expects a closure as an argument and `{ }` is used to define an anonymous closure. When calling a named closure, use the `( )` syntax.

Just run the workflow once more to check everything is still working:

```bash title="Test the named closure"
nextflow run main.nf
```

```console title="Samples using named closure"
 N E X T F L O W   ~  version 25.04.3

Launching `main.nf` [angry_meninsky] DSL2 - revision: 2edc226b1d

[[id:patientA, repeat:1], [id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam, [id:patientA, repeat:1, type:tumor], patientA_rep1_tumor.bam]
[[id:patientA, repeat:2], [id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam, [id:patientA, repeat:2, type:tumor], patientA_rep2_tumor.bam]
[[id:patientB, repeat:1], [id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam, [id:patientB, repeat:1, type:tumor], patientB_rep1_tumor.bam]
[[id:patientC, repeat:1], [id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam, [id:patientC, repeat:1, type:tumor], patientC_rep1_tumor.bam]
```

Using a named closure allows us to reuse the same transformation in multiple places, reducing the risk of errors and making the code more readable and maintainable.

### 3.5. Reduce duplication of data

We have a lot of duplicated data in our workflow. Each item in the joined samples repeats the `id` and `repeat` fields. Since this information is already available in the grouping key, we can avoid this redundancy. As a reminder, our current data structure looks like this:

```groovy
[
  [
    "id": "sampleC",
    "repeat": "1",
  ],
  [
    "id": "sampleC",
    "repeat": "1",
    "type": "normal",
  ],
  "sampleC_rep1_normal.bam"
  [
    "id": "sampleC",
    "repeat": "1",
    "type": "tumor",
  ],
  "sampleC_rep1_tumor.bam"
]
```

Since the `id` and `repeat` fields are available in the grouping key, let's remove them from the rest of each channel item to avoid duplication. We can do this by using the `subMap` method to create a new map with only the `type` field. This approach allows us to maintain all necessary information while eliminating redundancy in our data structure.

=== "After"

    ```groovy title="main.nf" linenums="8" hl_lines="1"
        getSampleIdAndReplicate = { meta, bam -> [ meta.subMap(['id', 'repeat']), meta.subMap(['type']), file(bam) ] }
    ```

=== "Before"

    ```groovy title="main.nf" linenums="8" hl_lines="1"
        getSampleIdAndReplicate = { meta, bam -> [ meta.subMap(['id', 'repeat']), meta, file(bam) ] }
    ```

Now the closure returns a tuple where the first element contains the `id` and `repeat` fields, and the second element contains only the `type` field. This eliminates redundancy by storing the `id` and `repeat` information once in the grouping key, while maintaining all necessary information.

Run the workflow to see what this looks like:

```bash title="Test data deduplication"
nextflow run main.nf
```

```console title="Deduplicated sample data"
[[id:patientA, repeat:1], [type:normal], /workspaces/training/side-quests/splitting_and_grouping/patientA_rep1_normal.bam, [type:tumor], /workspaces/training/side-quests/splitting_and_grouping/patientA_rep1_tumor.bam]
[[id:patientA, repeat:2], [type:normal], /workspaces/training/side-quests/splitting_and_grouping/patientA_rep2_normal.bam, [type:tumor], /workspaces/training/side-quests/splitting_and_grouping/patientA_rep2_tumor.bam]
[[id:patientB, repeat:1], [type:normal], /workspaces/training/side-quests/splitting_and_grouping/patientB_rep1_normal.bam, [type:tumor], /workspaces/training/side-quests/splitting_and_grouping/patientB_rep1_tumor.bam]
[[id:patientC, repeat:1], [type:normal], /workspaces/training/side-quests/splitting_and_grouping/patientC_rep1_normal.bam, [type:tumor], /workspaces/training/side-quests/splitting_and_grouping/patientC_rep1_tumor.bam]
```

We can see we only state the `id` and `repeat` fields once in the grouping key and we have the `type` field in the sample data. We haven't lost any information but we managed to make our channel contents more succinct.

### 3.6. Remove redundant information

We removed duplicated information above, but we still have some other redundant information in our channels.

In the beginning, we separated the normal and tumor samples using `filter`, then joined them based on `id` and `repeat` keys. The `join` operator preserves the order in which tuples are merged, so in our case, with normal samples on the left side and tumor samples on the right, the resulting channel maintains this structure: `id, <normal elements>, <tumor elements>`.

Since we know the position of each element in our channel, we can simplify the structure further by dropping the `[type:normal]` and `[type:tumor]` metadata.

=== "After"

    ```groovy title="main.nf" linenums="8" hl_lines="1"
        getSampleIdAndReplicate = { meta, file -> [ meta.subMap(['id', 'repeat']), file ] }
    ```

=== "Before"

    ```groovy title="main.nf" linenums="8" hl_lines="1"
        getSampleIdAndReplicate = { meta, file -> [ meta.subMap(['id', 'repeat']), meta.subMap(['type']), file ] }
    ```

Run again to see the result:

```bash title="Test streamlined data structure"
nextflow run main.nf
```

```console title="Streamlined sample data"
 N E X T F L O W   ~  version 25.04.3

Launching `main.nf` [confident_leavitt] DSL2 - revision: a2303895bd

[[id:patientA, repeat:1], patientA_rep1_normal.bam, patientA_rep1_tumor.bam]
[[id:patientA, repeat:2], patientA_rep2_normal.bam, patientA_rep2_tumor.bam]
[[id:patientB, repeat:1], patientB_rep1_normal.bam, patientB_rep1_tumor.bam]
[[id:patientC, repeat:1], patientC_rep1_normal.bam, patientC_rep1_tumor.bam]
```

### Takeaway

In this section, you've learned:

- **Manipulating Tuples**: How to use `map` to isolate a field in a tuple
- **Joining Tuples**: How to use `join` to combine tuples based on the first field
- **Creating Joining Keys**: How to use `subMap` to create a new joining key
- **Named Closures**: How to use a named closure in map
- **Multiple Field Joining**: How to join on multiple fields for more precise matching
- **Data Structure Optimization**: How to streamline channel structure by removing redundant information

You now have a workflow that can split a samplesheet, filter the normal and tumor samples, join them together by sample ID and replicate number, then print the results.

This is a common pattern in bioinformatics workflows where you need to match up samples or other types of data after processing independently, so it is a useful skill. Next, we will look at repeating a sample multiple times.

## 4. Spread samples intervals

A key pattern in bioinformatics workflows is distributing analysis across genomic regions. For instance, variant calling can be parallelized by dividing the genome into intervals (like chromosomes or smaller regions). This parallelization strategy significantly improves pipeline efficiency by distributing computational load across multiple cores or nodes, reducing overall execution time.

In the following section, we'll demonstrate how to distribute our sample data across multiple genomic intervals. We'll pair each sample with every interval, allowing parallel processing of different genomic regions. This will multiply our dataset size by the number of intervals, creating multiple independent analysis units that can be brought back together later.

### 4.1. Spread samples over intervals using `combine`

Let's start by creating a channel of intervals. To keep life simple, we will just use 3 intervals we will manually define. In a real workflow, you could read these in from a file input or even create a channel with lots of interval files.

=== "After"

    ```groovy title="main.nf" linenums="17" hl_lines="3"
            .join(ch_tumor_samples)
        ch_intervals = channel.of('chr1', 'chr2', 'chr3')
    ```

=== "Before"

    ```groovy title="main.nf" linenums="17" hl_lines="2"
            .join(ch_tumor_samples)
        ch_joined_samples.view()
    ```

Now remember, we want to repeat each sample for each interval. This is sometimes referred to as the Cartesian product of the samples and intervals. We can achieve this by using the [`combine` operator](https://www.nextflow.io/docs/latest/operator.html#combine). This will take every item from channel 1 and repeat it for each item in channel 2. Let's add a combine operator to our workflow:

=== "After"

    ```groovy title="main.nf" linenums="18" hl_lines="3-5"
        ch_intervals = channel.of('chr1', 'chr2', 'chr3')

        ch_combined_samples = ch_joined_samples
            .combine(ch_intervals)
            .view()
    ```

=== "Before"

    ```groovy title="main.nf" linenums="18"
        ch_intervals = channel.of('chr1', 'chr2', 'chr3')
    ```

Now let's run it and see what happens:

```bash title="Test the combine operation"
nextflow run main.nf
```

```console title="Samples combined with intervals"
 N E X T F L O W   ~  version 25.04.3

Launching `main.nf` [mighty_tesla] DSL2 - revision: ae013ab70b

[[id:patientA, repeat:1], patientA_rep1_normal.bam, patientA_rep1_tumor.bam, chr1]
[[id:patientA, repeat:1], patientA_rep1_normal.bam, patientA_rep1_tumor.bam, chr2]
[[id:patientA, repeat:1], patientA_rep1_normal.bam, patientA_rep1_tumor.bam, chr3]
[[id:patientA, repeat:2], patientA_rep2_normal.bam, patientA_rep2_tumor.bam, chr1]
[[id:patientA, repeat:2], patientA_rep2_normal.bam, patientA_rep2_tumor.bam, chr2]
[[id:patientA, repeat:2], patientA_rep2_normal.bam, patientA_rep2_tumor.bam, chr3]
[[id:patientB, repeat:1], patientB_rep1_normal.bam, patientB_rep1_tumor.bam, chr1]
[[id:patientB, repeat:1], patientB_rep1_normal.bam, patientB_rep1_tumor.bam, chr2]
[[id:patientB, repeat:1], patientB_rep1_normal.bam, patientB_rep1_tumor.bam, chr3]
[[id:patientC, repeat:1], patientC_rep1_normal.bam, patientC_rep1_tumor.bam, chr1]
[[id:patientC, repeat:1], patientC_rep1_normal.bam, patientC_rep1_tumor.bam, chr2]
[[id:patientC, repeat:1], patientC_rep1_normal.bam, patientC_rep1_tumor.bam, chr3]
```

Success! We have repeated every sample for every single interval in our 3 interval list. We've effectively tripled the number of items in our channel. It's a little hard to read though, so in the next section we will tidy it up.

### 4.2. Organise the channel

We can use the `map` operator to tidy and refactor our sample data so it's easier to understand. Let's move the intervals string to the joining map at the first element.

=== "After"

    ```groovy title="main.nf" linenums="20" hl_lines="3-9"
        ch_combined_samples = ch_joined_samples
            .combine(ch_intervals)
            .map { grouping_key, normal, tumor, interval ->
                [
                    grouping_key + [interval: interval],
                    normal,
                    tumor
                ]
            }
            .view()
    ```

=== "Before"

    ```groovy title="main.nf" linenums="20"
        ch_combined_samples = ch_joined_samples
            .combine(ch_intervals)
            .view()
    ```

Let's break down what this map operation does step by step.

First, we use named parameters to make the code more readable. By using the names `grouping_key`, `normal`, `tumor` and `interval`, we can refer to the elements in the tuple by name instead of by index:

```groovy
        .map { grouping_key, normal, tumor, interval ->
```

Next, we combine the `grouping_key` with the `interval` field. The `grouping_key` is a map containing `id` and `repeat` fields. We create a new map with the `interval` and merge them using Groovy's map addition (`+`):

```groovy
                grouping_key + [interval: interval],
```

Finally, we return this as a tuple with three elements: the combined metadata map, the normal sample file, and the tumor sample file:

```groovy
            [
                grouping_key + [interval: interval],
                normal,
                tumor
            ]
```

Let's run it again and check the channel contents:

```bash title="Test the reorganized structure"
nextflow run main.nf
```

```console title="Samples combined with intervals"
 N E X T F L O W   ~  version 25.04.3

Launching `main.nf` [sad_hawking] DSL2 - revision: 1f6f6250cd

[[id:patientA, interval:chr1], patientA_rep1_normal.bam, patientA_rep1_tumor.bam]
[[id:patientA, interval:chr2], patientA_rep1_normal.bam, patientA_rep1_tumor.bam]
[[id:patientA, interval:chr3], patientA_rep1_normal.bam, patientA_rep1_tumor.bam]
[[id:patientA, interval:chr1], patientA_rep2_normal.bam, patientA_rep2_tumor.bam]
[[id:patientA, interval:chr2], patientA_rep2_normal.bam, patientA_rep2_tumor.bam]
[[id:patientA, interval:chr3], patientA_rep2_normal.bam, patientA_rep2_tumor.bam]
[[id:patientB, interval:chr1], patientB_rep1_normal.bam, patientB_rep1_tumor.bam]
[[id:patientB, interval:chr2], patientB_rep1_normal.bam, patientB_rep1_tumor.bam]
[[id:patientB, interval:chr3], patientB_rep1_normal.bam, patientB_rep1_tumor.bam]
[[id:patientC, interval:chr1], patientC_rep1_normal.bam, patientC_rep1_tumor.bam]
[[id:patientC, interval:chr2], patientC_rep1_normal.bam, patientC_rep1_tumor.bam]
[[id:patientC, interval:chr3], patientC_rep1_normal.bam, patientC_rep1_tumor.bam]
```

Using `map` to coerce your data into the correct structure can be tricky, but it's crucial for effective data manipulation.

We now have every sample repeated across all genomic intervals, creating multiple independent analysis units that can be processed in parallel. But what if we want to bring related samples back together? In the next section, we'll learn how to group samples that share common attributes.

### Takeaway

In this section, you've learned:

- **Spreading samples over intervals**: How to use `combine` to repeat samples over intervals
- **Creating Cartesian products**: How to generate all combinations of samples and intervals
- **Organizing channel structure**: How to use `map` to restructure data for better readability
- **Parallel processing preparation**: How to set up data for distributed analysis

## 5. Aggregating samples using `groupTuple`

In the previous sections, we learned how to split data from an input file and filter by specific fields (in our case normal and tumor samples). But this only covers a single type of joining. What if we want to group samples by a specific attribute? For example, instead of joining matched normal-tumor pairs, we might want to process all samples from "sampleA" together regardless of their type. This pattern is common in bioinformatics workflows where you may want to process related samples separately for efficiency reasons before comparing or combining the results at the end.

Nextflow includes built in methods to do this, the main one we will look at is `groupTuple`.

Let's start by grouping all of our samples that have the same `id` and `interval` fields, this would be typical of an analysis where we wanted to group technical replicates but keep meaningfully different samples separated.

To do this, we should separate out our grouping variables so we can use them in isolation.

The first step is similar to what we did in the previous section. We must isolate our grouping variable as the first element of the tuple. Remember, our first element is currently a map of `id`, `repeat` and `interval` fields:

```groovy title="main.nf" linenums="1"
{
  "id": "sampleA",
  "repeat": "1",
  "interval": "chr1"
}
```

We can reuse the `subMap` method from before to isolate our `id` and `interval` fields from the map. Like before, we will use `map` operator to apply the `subMap` method to the first element of the tuple for each sample.

=== "After"

    ```groovy title="main.nf" linenums="20" hl_lines="11-19"
        ch_combined_samples = ch_joined_samples
            .combine(ch_intervals)
            .map { grouping_key, normal, tumor, interval ->
                [
                    grouping_key + [interval: interval],
                    normal,
                    tumor
                ]
            }

        ch_grouped_samples = ch_combined_samples
            .map { grouping_key, normal, tumor ->
                [
                    grouping_key.subMap('id', 'interval'),
                    normal,
                    tumor
                ]
              }
              .view()
    ```

=== "Before"

    ```groovy title="main.nf" linenums="20" hl_lines="10"
        ch_combined_samples = ch_joined_samples
            .combine(ch_intervals)
            .map { grouping_key, normal, tumor, interval ->
                [
                    grouping_key + [interval: interval],
                    normal,
                    tumor
                ]
            }
            .view()
    ```

Let's run it again and check the channel contents:

```bash title="Test grouping key isolation"
nextflow run main.nf
```

```console title="Samples prepared for grouping"
 N E X T F L O W   ~  version 25.04.3

Launching `main.nf` [hopeful_brenner] DSL2 - revision: 7f4f7fea76

[[id:patientA, interval:chr1], patientA_rep1_normal.bam, patientA_rep1_tumor.bam]
[[id:patientA, interval:chr2], patientA_rep1_normal.bam, patientA_rep1_tumor.bam]
[[id:patientA, interval:chr3], patientA_rep1_normal.bam, patientA_rep1_tumor.bam]
[[id:patientA, interval:chr1], patientA_rep2_normal.bam, patientA_rep2_tumor.bam]
[[id:patientA, interval:chr2], patientA_rep2_normal.bam, patientA_rep2_tumor.bam]
[[id:patientA, interval:chr3], patientA_rep2_normal.bam, patientA_rep2_tumor.bam]
[[id:patientB, interval:chr1], patientB_rep1_normal.bam, patientB_rep1_tumor.bam]
[[id:patientB, interval:chr2], patientB_rep1_normal.bam, patientB_rep1_tumor.bam]
[[id:patientB, interval:chr3], patientB_rep1_normal.bam, patientB_rep1_tumor.bam]
[[id:patientC, interval:chr1], patientC_rep1_normal.bam, patientC_rep1_tumor.bam]
[[id:patientC, interval:chr2], patientC_rep1_normal.bam, patientC_rep1_tumor.bam]
[[id:patientC, interval:chr3], patientC_rep1_normal.bam, patientC_rep1_tumor.bam]
```

We can see that we have successfully isolated the `id` and `interval` fields, but not grouped the samples yet.

!!! note

    We are discarding the `replicate` field here. This is because we don't need it for further downstream processing. After completing this tutorial, see if you can include it without affecting the later grouping!

Let's now group the samples by this new grouping element, using the [`groupTuple` operator](https://www.nextflow.io/docs/latest/operator.html#grouptuple).

=== "After"

    ```groovy title="main.nf" linenums="30" hl_lines="9"
        ch_grouped_samples = ch_combined_samples
            .map { grouping_key, normal, tumor ->
                [
                    grouping_key.subMap('id', 'interval'),
                    normal,
                    tumor
                ]
              }
              .groupTuple()
              .view()
    ```

=== "Before"

    ```groovy title="main.nf" linenums="30"
        ch_grouped_samples = ch_combined_samples
            .map { grouping_key, normal, tumor ->
                [
                    grouping_key.subMap('id', 'interval'),
                    normal,
                    tumor
                ]
              }
              .view()
    ```

That's all there is to it! We just added a single line of code. Let's see what happens when we run it:

```bash title="Test the groupTuple operation"
nextflow run main.nf
```

```console title="Grouped samples by ID and interval"
 N E X T F L O W   ~  version 25.04.3

Launching `main.nf` [friendly_jang] DSL2 - revision: a1bee1c55d

[[id:patientA, interval:chr1], [patientA_rep1_normal.bam, patientA_rep2_normal.bam], [patientA_rep1_tumor.bam, patientA_rep2_tumor.bam]]
[[id:patientA, interval:chr2], [patientA_rep1_normal.bam, patientA_rep2_normal.bam], [patientA_rep1_tumor.bam, patientA_rep2_tumor.bam]]
[[id:patientA, interval:chr3], [patientA_rep1_normal.bam, patientA_rep2_normal.bam], [patientA_rep1_tumor.bam, patientA_rep2_tumor.bam]]
[[id:patientB, interval:chr1], [patientB_rep1_normal.bam], [patientB_rep1_tumor.bam]]
[[id:patientB, interval:chr2], [patientB_rep1_normal.bam], [patientB_rep1_tumor.bam]]
[[id:patientB, interval:chr3], [patientB_rep1_normal.bam], [patientB_rep1_tumor.bam]]
[[id:patientC, interval:chr1], [patientC_rep1_normal.bam], [patientC_rep1_tumor.bam]]
[[id:patientC, interval:chr2], [patientC_rep1_normal.bam], [patientC_rep1_tumor.bam]]
[[id:patientC, interval:chr3], [patientC_rep1_normal.bam], [patientC_rep1_tumor.bam]]
```

Note our data has changed structure and within each channel element the files now contained in tuples like `[patientA_rep1_normal.bam, patientA_rep2_normal.bam]`. This is because when we use `groupTuple`, Nextflow combines the single files for each sample of a group. This is important to remember when trying to handle the data downstream.

!!! note

    [`transpose`](https://www.nextflow.io/docs/latest/reference/operator.html#transpose) is the opposite of groupTuple. It unpacks the items in a channel and flattens them. Try and add `transpose` and undo the grouping we performed above!

### Takeaway

In this section, you've learned:

- **Grouping related samples**: How to use `groupTuple` to aggregate samples by common attributes
- **Isolating grouping keys**: How to use `subMap` to extract specific fields for grouping
- **Handling grouped data structures**: How to work with the nested structure created by `groupTuple`
- **Technical replicate handling**: How to group samples that share the same experimental conditions

---

## Summary

In this side quest, you've learned how to split and group data using channels.

By modifying the data as it flows through the pipeline, you can construct a scalable pipeline without using loops or while statements, offering several advantages over more traditional approaches:

- We can scale to as many or as few inputs as we want with no additional code
- We focus on handling the flow of data through the pipeline, instead of iteration
- We can be as complex or simple as required
- The pipeline becomes more declarative, focusing on what should happen rather than how it should happen
- Nextflow will optimize execution for us by running independent operations in parallel

Mastering these channel operations will enable you to build flexible, scalable pipelines that handle complex data relationships without resorting to loops or iterative programming, allowing Nextflow to optimize execution and parallelize independent operations automatically.

### Key patterns

1.  **Creating structured input data:** Starting from a CSV file with meta maps (building on patterns from [Metadata in workflows](./metadata.md))

    ```groovy
    ch_samples = channel.fromPath("./data/samplesheet.csv")
        .splitCsv(header: true)
        .map{ row ->
          [[id:row.id, repeat:row.repeat, type:row.type], row.bam]
        }
    ```

2.  **Splitting data into separate channels:** We used `filter` to divide data into independent streams based on the `type` field

    ```groovy
    channel.filter { it.type == 'tumor' }
    ```

3.  **Joining matched samples:** We used `join` to recombine related samples based on `id` and `repeat` fields

    - Join two channels by key (first element of tuple)

    ```groovy
    tumor_ch.join(normal_ch)
    ```

    - Extract joining key and join by this value

    ```groovy
    tumor_ch.map { meta, file -> [meta.id, meta, file] }
        .join(
          normal_ch.map { meta, file -> [meta.id, meta, file] }
        )
    ```

    - Join on multiple fields using subMap

    ```groovy
    tumor_ch.map { meta, file -> [meta.subMap(['id', 'repeat']), meta, file] }
        .join(
          normal_ch.map { meta, file -> [meta.subMap(['id', 'repeat']), meta, file] }
        )
    ```

4.  **Distributing across intervals:** We used `combine` to create Cartesian products of samples with genomic intervals for parallel processing.

    ```groovy
    samples_ch.combine(intervals_ch)
    ```

5.  **Aggregating by grouping keys:** We used `groupTuple` to group by the first element in each tuple, thereby collecting samples sharing `id` and `interval` fields and merging technical replicates.

    ```groovy
    channel.groupTuple()
    ```

6.  **Optimizing the data structure:** We used `subMap` to extract specific fields and created a named closure for making transformations reusable.

    - Extract specific fields from a map

    ```groovy
    meta.subMap(['id', 'repeat'])
    ```

    - Use named closure for reusable transformations

    ```groovy
    getSampleIdAndReplicate = { meta, file -> [meta.subMap(['id', 'repeat']), file] }
    channel.map(getSampleIdAndReplicate)
    ```

### Additional resources

- [filter](https://www.nextflow.io/docs/latest/operator.html#filter)
- [map](https://www.nextflow.io/docs/latest/operator.html#map)
- [join](https://www.nextflow.io/docs/latest/operator.html#join)
- [groupTuple](https://www.nextflow.io/docs/latest/operator.html#grouptuple)
- [combine](https://www.nextflow.io/docs/latest/operator.html#combine)

---

## What's next?

Return to the [menu of Side Quests](./index.md) or click the button in the bottom right of the page to move on to the next topic in the list.
