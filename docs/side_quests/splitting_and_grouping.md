# Splitting and Grouping

Nextflow helps you work with your data in flexible ways. One of the most useful things you can do is split your data into different streams and then group related items back together. This capability is particularly valuable in bioinformatics workflows where you often need to process sub-groups of input data separately before combining results for comparison or joint analysis.

Think of it like sorting mail: you might first separate letters by their destination, process each pile differently, and then recombine items going to the same person. In Nextflow, we use special operators to do this with our scientific data.

Nextflow's channel system is at the heart of this flexibility. Channels connect different parts of your workflow, allowing data to flow through your analysis. You can create multiple channels from a single data source, process each channel differently, and then merge channels back together when needed. This approach lets you design workflows that naturally mirror the branching and converging paths of complex bioinformatics analyses.

In this side quest, we'll explore how to split and group data using Nextflow's powerful channel operators. We'll start with a CSV file containing information about different samples and their associated data, which we'll read and manipulate. By the end of this side quest, you'll be able to separate and combine data streams effectively, making your workflows more efficient and easier to understand.

You will:

- Read data from files using `splitCsv`
- Filter and transform data with `filter` and `map`
- Combine related data using `join` and `groupTuple`

These skills will help you build workflows that can handle multiple input files and different types of data efficiently.

---

## 0. Warmup

### 0.1. Prerequisites

Before taking on this side quest you should:

- Complete the [Hello Nextflow](../hello_nextflow/README.md) tutorial
- Understand basic Nextflow concepts (processes, channels, operators, working with files, meta data)
- Basic understanding of using sample-specific data (meta data): [Working with metadata](./metadata.md)

### 0.2. Starting Point

Let's move into the project directory.

```bash
cd side-quests/splitting_and_grouping
```

You'll find a `data` directory containing a samplesheet and a main workflow file.

```console title="Directory contents"
> tree
.
├── data
│   └── samplesheet.csv
└── main.nf
```

samplesheet.csv contains information about samples taken from different patients and their associated data. In particular, it contains information about the patient's ID, the sample's repeat number and type (normal or tumor), and the paths to the BAM files (which don't actually exist, but we will pretend they do).

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

Note there are 8 samples in total from 3 patients (patientA has 2 repeats), 4 normal and 4 tumor.

We're going to read in samplesheet.csv, then group and split the patients based on their data.

---

## 1. Read in sample data

### 1.1. Read in sample data with splitCsv

Let's start by reading in the sample data with `splitCsv`. In the main workflow file, you'll see that we've already started the workflow.

```groovy title="main.nf" linenums="1"
workflow {
    ch_samplesheet = Channel.fromPath("./data/data.csv")
}
```

!!! note

    Throughout this tutorial, we'll use the `ch_` prefix for all channel variables to clearly indicate they are Nextflow channels.

=== "After"

    ```groovy title="main.nf" linenums="2" hl_lines="1-3"
      ch_samples = Channel.fromPath("./data/data.csv")
          .splitCsv(header: true)
          .view()
    ```

=== "Before"

    ```groovy title="main.nf" linenums="2" hl_lines="1"
        ch_samplesheet = Channel.fromPath("./data/data.csv")
    ```

We can use the [`splitCsv` operator](https://www.nextflow.io/docs/latest/operator.html#splitcsv) to split the data into a channel of maps, where each map represents a row from the CSV file.

The `header: true` option tells Nextflow to use the first row of the CSV file as the header row, which will be used as keys for the values. Let's see what Nextflow can see after reading with splitCsv. To do this, we can use the `view` operator.

Run the pipeline:

```bash title="Read the data"
nextflow run main.nf
```

```console title="Read data with splitCsv"
 N E X T F L O W   ~  version 24.10.5

Launching `main.nf` [deadly_mercator] DSL2 - revision: bd6b0224e9

[id:patientA, repeat:1, type:normal, bam:patientA_rep1_normal.bam]
[id:patientA, repeat:1, type:tumor, bam:patientA_rep1_tumor.bam]
[id:patientA, repeat:2, type:normal, bam:patientB_rep1_normal.bam]
[id:patientA, repeat:2, type:tumor, bam:patientB_rep1_tumor.bam]
[id:patientB, repeat:1, type:normal, bam:patientC_rep1_normal.bam]
[id:patientB, repeat:1, type:tumor, bam:patientC_rep1_tumor.bam]
[id:patientC, repeat:1, type:normal, bam:patientD_rep1_normal.bam]
[id:patientC, repeat:1, type:tumor, bam:patientD_rep1_tumor.bam]
```

Each row from the CSV file has been converted into a map with keys matching the header row.
Each map contains:

- `id`: The patient identifier (patientA, patientB, patientC)
- `repeat`: The replicate number (1 or 2)
- `type`: The sample type (normal or tumor)
- `bam`: Path to the BAM file

This format makes it easy to access specific fields from each sample via their keys in the map. We can access the BAM file path with the `bam` key, but also any of the 'metadata' fields that describe the file via `id`, `repeat`, `type`.

!!!Note

    For a more extensive introduction on working with metadatadata, you can work through the training [Working with metadata](./metadata.md)

Let's separate the metadata from the files:

=== "After"

    ```groovy title="main.nf" linenums="2" hl_lines="3-5"
        ch_samples = Channel.fromPath("./data/samplesheet.csv")
            .splitCsv(header: true)
            .map{ row ->
              [[id:row.id, repeat:row.repeat, type:row.type], row.bam]
            }
            .view()
    ```

=== "Before"

    ```groovy title="main.nf" linenums="2"
        ch_samples = Channel.fromPath("./data/samplesheet.csv")
            .splitCsv(header: true)
            .view()
    ```

and re-run the pipeline:

```bash title="Convert the sample information into a map"
nextflow run main.nf
```

```console title="Convert the sample information into a map"
 N E X T F L O W   ~  version 24.10.5

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

We separated the sample meta data from the files into a map.
We now have a channel of maps and files, each representing a row from the input sample sheet, which we will use in this training to split and group our workload.

### Takeaway

In this section, you've learned:

- **Reading in a data sheet**: How to read in data sheet with `splitCsv`
- **Combining patient-specific information**: Using groovy maps to hold information about a patient
- **Viewing data**: How to use `view` to print the data

---

## 2. Filter and transform data

### 2.1. Filter data with `filter`

We can use the [`filter` operator](https://www.nextflow.io/docs/latest/operator.html#filter) to filter the data based on a condition. Let's say we only want to process normal samples. We can do this by filtering the data based on the `type` field. Let's insert this before the `view` operator.

=== "After"

    ```groovy title="main.nf" linenums="2" hl_lines="6"
        ch_samples = Channel.fromPath("./data/samplesheet.csv")
            .splitCsv(header: true)
            .map{ row ->
              [[id:row.id, repeat:row.repeat, type:row.type], row.bam]
            }
            .filter { meta, file -> meta.type == 'normal' }
            .view()
    ```

=== "Before"

    ```groovy title="main.nf" linenums="2"
        ch_samples = Channel.fromPath("./data/samplesheet.csv")
            .splitCsv(header: true)
            .map{ row ->
              [[id:row.id, repeat:row.repeat, type:row.type], row.bam]
            }
            .view()
    ```

```bash title="View normal samples"
nextflow run main.nf
```

```console title="View normal samples"
 N E X T F L O W   ~  version 24.10.5

Launching `main.nf` [admiring_brown] DSL2 - revision: 194d61704d

[[id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam]
[[id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam]
[[id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam]
[[id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam]
```

We have successfully filtered the data to only include normal samples. Let's recap how this works. The `filter` operator takes a closure that is applied to each element in the channel. If the closure returns `true`, the element is included in the output channel. If the closure returns `false`, the element is excluded from the output channel.

In this case, we want to keep only the samples where `meta.type == 'normal'`. In the closure, we use the tuple `meta,file` to refer to each sample in the channel. We can then access the type of the particular sample with `meta.type` and check if it is equal to `'normal'`. If it is, the sample is included in the output channel. If it is not, the sample is excluded from the output channel.

```groovy title="main.nf" linenums="4"
    .filter { meta, file -> meta.type == 'normal' }
```

### 2.2. Filter to just the tumor samples

Currently we're applying the filter to the channel created directly from the CSV, but we want to filter this in more ways than one, so let's re-write the logic to create a separate filtered channel for normal samples:

=== "After"

    ```groovy title="main.nf" linenums="2" hl_lines="6 7"
        ch_samples = Channel.fromPath("./data/samplesheet.csv")
            .splitCsv(header: true)
            .map{ row ->
                [[id:row.id, repeat:row.repeat, type:row.type], row.bam]
            }
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
        ch_normal_samples.view()
    ```

=== "Before"

    ```groovy title="main.nf" linenums="2"
        ch_samples = Channel.fromPath("./data/samplesheet.csv")
            .splitCsv(header: true)
            .filter { meta, file -> meta.type == 'normal' }
            .view()
    ```

Once again, run the pipeline to see the results:

```bash title="View normal samples"
nextflow run main.nf
```

```console title="View normal samples"
 N E X T F L O W   ~  version 24.10.5

Launching `main.nf` [trusting_poisson] DSL2 - revision: 639186ee74

[[id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam]
[[id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam]
[[id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam]
[[id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam]
```

Success! We have filtered the data to only include normal samples. Note that we can use view and save the new channel. If we wanted, we still have access to the tumor samples within the `ch_samples` channel. Since we managed it for the normal samples, let's do it for the tumor samples as well:

=== "After"

    ```groovy title="main.nf" linenums="8" hl_lines="8-11"
        ch_samples = Channel.fromPath("./data/samplesheet.csv")
            .splitCsv(header: true)
            .map{ row ->
                [[id:row.id, repeat:row.repeat, type:row.type], row.bam]
            }
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
        ch_normal_samples.view{'Normal sample: ' + it}
        ch_tumor_samples.view{'Tumour sample: ' + it}
    ```

=== "Before"

    ```groovy title="main.nf" linenums="2"
        ch_samples = Channel.fromPath("./data/samplesheet.csv")
            .splitCsv(header: true)
            .map{ row ->
                [[id:row.id, repeat:row.repeat, type:row.type], row.bam]
            }
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
        ch_normal_samples.view()
    ```

```bash title="View tumor samples"
nextflow run main.nf
```

```console title="View tumor samples"
 N E X T F L O W   ~  version 24.10.5

Launching `main.nf` [big_bernard] DSL2 - revision: 897c9e44cc

Tumour sample: [[id:patientA, repeat:1, type:tumor], patientA_rep1_tumor.bam]
Normal sample: [[id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam]
Tumour sample: [[id:patientA, repeat:2, type:tumor], patientA_rep2_tumor.bam]
Normal sample: [[id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam]
Normal sample: [[id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam]
Normal sample: [[id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam]
Tumour sample: [[id:patientB, repeat:1, type:tumor], patientB_rep1_tumor.bam]
Tumour sample: [[id:patientC, repeat:1, type:tumor], patientC_rep1_tumor.bam]
```

We've managed to separate out the normal and tumor samples into two different channels, and note the closure we supplied to `view()` to label them differently: `ch_tumor_samples.view{'Tumour sample: ' + it}`.

### Takeaway

In this section, you've learned:

- **Filtering data**: How to filter data with `filter`
- **Splitting data**: How to split data into different channels based on a condition
- **Viewing data**: How to use `view` to print the data and label output from different channels

We've now separated out the normal and tumor samples into two different channels. Next, we'll join the normal and tumor samples on the `id` field.

---

## 3. Join on patient ID

In the previous section, we separated out the normal and tumor samples into two different channels. These could be processed independently using specific processes or workflows based on their type. But what happens when we want to compare the normal and tumor samples from the same patient? At this point, we need to join them back together making sure to match the samples based on their `id` field.

Nextflow includes many methods for combining channels, but in this case the most appropriate operator is [`join`](https://www.nextflow.io/docs/latest/operator.html#join). If you are familiar with SQL, it acts like the `JOIN` operation, where we specify the key to join on and the type of join to perform.

### 3.1. Use `map` and `join` to combine based on patient ID

If we check the [`join`](https://www.nextflow.io/docs/latest/operator.html#join) documentation, we can see that by default it joins two channels based on the first item in each tuple. If you don't have the console output still available, let's run the pipeline to check our data structure and see how we need to modify it to join on the `id` field.

```bash title="View normal and tumor samples"
nextflow run main.nf
```

```console title="View normal and tumor samples"
 N E X T F L O W   ~  version 24.10.5

Launching `main.nf` [loving_bardeen] DSL2 - revision: 012d38e59f

Tumour sample: [[id:patientA, repeat:1, type:tumor], patientA_rep1_tumor.bam]
Normal sample: [[id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam]
Tumour sample: [[id:patientA, repeat:2, type:tumor], patientA_rep2_tumor.bam]
Normal sample: [[id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam]
Normal sample: [[id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam]
Normal sample: [[id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam]
Tumour sample: [[id:patientB, repeat:1, type:tumor], patientB_rep1_tumor.bam]
Tumour sample: [[id:patientC, repeat:1, type:tumor], patientC_rep1_tumor.bam]
```

We can see that the `id` field is the first element in each meta map. For `join` to work, we should isolate the `id` field in each tuple. After that, we can simply use the `join` operator to combine the two channels.

To isolate the `id` field, we can use the [`map` operator](https://www.nextflow.io/docs/latest/operator.html#map) to create a new tuple with the `id` field as the first element.

=== "After"

    ```groovy title="main.nf" linenums="8" hl_lines="3 6"
        ch_normal_samples = ch_samplesheet
            .filter { meta, file -> meta.type == 'normal' }
            .map { meta, file -> [meta.id, meta, file] }
        ch_tumor_samples = ch_samplesheet
            .filter { meta, file -> meta.type == 'tumor' }
            .map { meta, file -> [meta.id, meta, file] }
        ch_normal_samples.view{'Normal sample: ' + it}
        ch_tumor_samples.view{'Tumour sample: ' + it}
    ```

=== "Before"

    ```groovy title="main.nf" linenums="8"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
        ch_normal_samples.view{'Normal sample: ' + it}
        ch_tumor_samples.view{'Tumour sample: ' + it}
    ```

```bash title="View normal and tumor samples with ID as element 0"
nextflow run main.nf
```

```console title="View normal and tumor samples with ID as element 0"
 N E X T F L O W   ~  version 24.10.5

Launching `main.nf` [dreamy_sax] DSL2 - revision: 882ae9add4

Tumour sample: [patientA, [id:patientA, repeat:1, type:tumor], patientA_rep1_tumor.bam]
Tumour sample: [patientA, [id:patientA, repeat:2, type:tumor], patientA_rep2_tumor.bam]
Tumour sample: [patientB, [id:patientB, repeat:1, type:tumor], patientB_rep1_tumor.bam]
Tumour sample: [patientC, [id:patientC, repeat:1, type:tumor], patientC_rep1_tumor.bam]
'Normal sample: [patientA, [id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam]
'Normal sample: [patientA, [id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam]
'Normal sample: [patientB, [id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam]
'Normal sample: [patientC, [id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam]
```

It might be subtle, but you should be able to see the first element in each tuple is the `id` field. Now we can use the `join` operator to combine the two channels based on the `id` field.

Once again, we will use `view` to print the joined outputs.

=== "After"
    ```groovy title="main.nf" linenums="11" hl_lines="7-9"
        ch_normal_samples = ch_samplesheet
            .filter { meta, file -> meta.type == 'normal' }
            .map { meta, file -> [meta.id, meta, file] }
        ch_tumor_samples = ch_samplesheet
            .filter { meta, file -> meta.type == 'tumor' }
            .map { meta, file -> [meta.id, meta, file] }
        ch_joined_samples = ch_normal_samples
            .join(ch_tumor_patients)
        ch_joined_samples.view()
    ```

=== "Before"

    ```groovy title="main.nf" linenums="2" hl_lines="3 6"
        ch_normal_samples = ch_samplesheet
            .filter { meta, file -> meta.type == 'normal' }
            .map { meta, file -> [meta.id, meta, file] }
        ch_tumor_samples = ch_samplesheet
            .filter { meta, file -> meta.type == 'tumor' }
            .map { meta, file -> [meta.id, meta, file] }
        ch_normal_samples.view{'Normal sample: ' + it}
        ch_tumor_samples.view{'Tumour sample: ' + it}
    ```

```bash title="View normal and tumor samples"
nextflow run main.nf
```

```console title="View joined normal and tumor samples"
 N E X T F L O W   ~  version 24.10.5

Launching `main.nf` [elegant_waddington] DSL2 - revision: c552f22069

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

We have 2 replicates for sampleA, but only 1 for sampleB and sampleC. In this case we were able to join them effectively by using the `id` field, but what would happen if they were out of sync? We could mix up the normal and tumor samples from different replicates! This could be disastrous!

To avoid this, we can join on multiple fields. There are actually multiple ways to achieve this but we are going to focus on creating a new joining key which includes both the sample `id` and `replicate` number.

Let's start by creating a new joining key. We can do this in the same way as before by using the [`map` operator](https://www.nextflow.io/docs/latest/operator.html#map) to create a new tuple with the `id` and `repeat` fields as the first element.

=== "After"

    ```groovy title="main.nf" linenums="8" hl_lines="3 6"
        ch_normal_samples = ch_samplesheet
            .filter { meta, file -> meta.type == 'normal' }
            .map { meta, file -> [[meta.id, meta.repeat], meta, file] }
        ch_tumor_samples = ch_samplesheet
            .filter { meta, file -> meta.type == 'tumor' }
            .map { meta, file -> [[meta.id, meta.repeat], meta, file] }
    ```

=== "Before"

    ```groovy title="main.nf" linenums="8"
        ch_normal_samples = ch_samplesheet
            .filter { meta, file -> meta.type == 'normal' }
            .map { meta, file -> [meta.id, meta, file] }
        ch_tumor_samples = ch_samplesheet
            .filter { meta, file -> meta.type == 'tumor' }
            .map { meta, file -> [meta.id, meta, file] }
    ```

Now we should see the join is occurring but using both the `id` and `repeat` fields.

```bash title="View normal and tumor samples"
nextflow run main.nf
```

```console title="View normal and tumor samples"
 N E X T F L O W   ~  version 24.10.5

Launching `main.nf` [prickly_wing] DSL2 - revision: 3bebf22dee

[[patientA, 1], [id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam, [id:patientA, repeat:1, type:tumor], patientA_rep1_tumor.bam]
[[patientA, 2], [id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam, [id:patientA, repeat:2, type:tumor], patientA_rep2_tumor.bam]
[[patientB, 1], [id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam, [id:patientB, repeat:1, type:tumor], patientB_rep1_tumor.bam]
[[patientC, 1], [id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam, [id:patientC, repeat:1, type:tumor], patientC_rep1_tumor.bam]
```

Note how we have a tuple of two elements (`id` and `repeat` fields) as the first element of each joined result. This demonstrates how complex items can be used as a joining key, enabling fairly intricate matching between samples from the same conditions.

If you want to explore more ways to join on different keys, check out the [join operator documentation](https://www.nextflow.io/docs/latest/operator.html#join) for additional options and examples.

### 3.3. Use subMap to create a new joining key

We have an issue from the above example. We have lost the field names from the original joining key, i.e. the `id` and `repeat` fields are just a list of two values. If we want to retain the field names so we can access them later by name we can use the [`subMap` method](<https://docs.groovy-lang.org/latest/html/groovy-jdk/java/util/Map.html#subMap(java.util.Collection)>).

The `subMap` method takes a map and returns a new map with only the key-value pairs specified in the argument. In this case we want to specify the `id` and `repeat` fields.

=== "After"

    ```groovy title="main.nf" linenums="4" hl_lines="3 6"
        ch_normal_samples = ch_samplesheet
            .filter { meta, file -> meta.type == 'normal' }
            .map { meta, file -> [meta.subMap(['id', 'repeat'], meta, file] }
        ch_tumor_samples = ch_samplesheet
            .filter { meta, file -> meta.type == 'tumor' }
            .map { meta, file -> [meta.subMap(['id', 'repeat'], meta, file] }
    ```

=== "Before"

    ```groovy title="main.nf" linenums="4" hl_lines="3 6"
        ch_normal_samples = ch_samplesheet
            .filter { meta, file -> meta.type == 'normal' }
            .map { meta, file -> [[meta.id, meta.repeat], meta, file] }
        ch_tumor_samples = ch_samplesheet
            .filter { meta, file -> meta.type == 'tumor' }
            .map { meta, file -> [[meta.id, meta.repeat], meta, file] }
    ```

```bash title="View normal and tumor samples"
nextflow run main.nf
```

```console title="View normal and tumor samples"
 N E X T F L O W   ~  version 24.10.5

Launching `main.nf` [curious_hopper] DSL2 - revision: 90283e523d

[[id:patientA, repeat:1], [id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam, [id:patientA, repeat:1, type:tumor], patientA_rep1_tumor.bam]
[[id:patientA, repeat:2], [id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam, [id:patientA, repeat:2, type:tumor], patientA_rep2_tumor.bam]
[[id:patientB, repeat:1], [id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam, [id:patientB, repeat:1, type:tumor], patientB_rep1_tumor.bam]
[[id:patientC, repeat:1], [id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam, [id:patientC, repeat:1, type:tumor], patientC_rep1_tumor.bam]
```

Now we have a new joining key that not only includes the `id` and `repeat` fields but also retains the field names so we can access them later by name, e.g. `meta.id` and `meta.repeat`.

### 3.4. Use a named closure in map

Since we are re-using the same map in multiple places, we run the risk of introducing errors if we accidentally change the map in one place but not the other. To avoid this, we can use a named closure in the map. A named closure allows us to make a reusable function we can call later within a map.

To do so, first we define the closure as a new variable:

=== "After"

    ```groovy title="main.nf" linenums="2" hl_lines="7"
        ch_samples = Channel.fromPath("./data/samplesheet.csv")
            .splitCsv(header: true)
            .map{ row ->
                [[id:row.id, repeat:row.repeat, type:row.type], row.bam]
            }

        getSampleIdAndReplicate = { meta, bam -> [ meta.subMap(['id', 'repeat']), meta, file(bam) ] }

        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
    ```

=== "Before"

    ```groovy title="main.nf" linenums="2" hl_lines="6"
        ch_samples = Channel.fromPath("./data/samplesheet.csv")
            .splitCsv(header: true)
            .map{ row ->
                [[id:row.id, repeat:row.repeat, type:row.type], row.bam]
            }
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
    ```

We have taken the map we used previously and defined it as a named variable we can call later. We can also use this to convert the file path to a Path object using `file()` so that any process we passed the channel to could handle the file correctly (for more information see [Working with files](./working_with_files.md)).

Let's implement the closure in our workflow:

=== "After"

    ```groovy title="main.nf" linenums="7" hl_lines="3 6"
        ch_normal_samples = ch_samplesheet
            .filter { meta, file -> meta.type == 'normal' }
             .map ( getSampleIdAndReplicate )
        ch_tumor_samples = ch_samplesheet
            .filter { meta, file -> meta.type == 'tumor' }
             .map ( getSampleIdAndReplicate )

    ```

=== "Before"

    ```groovy title="main.nf" linenums="3 6"
        ch_normal_samples = ch_samplesheet
            .filter { meta, file -> meta.type == 'normal' }
            .map { meta, file -> [meta.subMap(['id', 'repeat'], meta, file] }
        ch_tumor_samples = ch_samplesheet
            .filter { meta, file -> meta.type == 'tumor' }
            .map { meta, file -> [meta.subMap(['id', 'repeat'], meta, file] }
    ```

!!! note

    The `map` operator has switched from using `{ }` to using `( )` to pass the closure as an argument. This is because the `map` operator expects a closure as an argument and `{ }` is used to define an anonymous closure. When calling a named closure, use the `( )` syntax.

```bash title="View normal and tumor samples"
nextflow run main.nf
```

```console title="View normal and tumor samples"
 N E X T F L O W   ~  version 24.10.5

Launching `main.nf` [angry_meninsky] DSL2 - revision: 2edc226b1d

[[id:patientA, repeat:1], [id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam, [id:patientA, repeat:1, type:tumor], patientA_rep1_tumor.bam]
[[id:patientA, repeat:2], [id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam, [id:patientA, repeat:2, type:tumor], patientA_rep2_tumor.bam]
[[id:patientB, repeat:1], [id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam, [id:patientB, repeat:1, type:tumor], patientB_rep1_tumor.bam]
[[id:patientC, repeat:1], [id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam, [id:patientC, repeat:1, type:tumor], patientC_rep1_tumor.bam]
```

Using a named closure in the map allows us to reuse the same map in multiple places which reduces our risk of introducing errors. It also makes the code more readable and easier to maintain.

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

Since the `id` and `repeat` fields are available in the grouping key, let's remove them from the patient data to avoid duplication. We can do this by using the `subMap` method to create a new map with only the `type` field. This approach allows us to maintain all necessary information while eliminating redundancy in our data structure.

=== "After"

    ```groovy title="main.nf" linenums="5" hl_lines="3-5"
    getSampleIdAndReplicate = { meta, bam -> [ meta.subMap(['id', 'repeat']), meta.subMap(['type']), file(bam) ] }
    ```

=== "Before"

    ```groovy title="main.nf" linenums="5"
        getSampleIdAndReplicate = { meta, bam -> [ meta.subMap(['id', 'repeat']), meta, file(bam) ] }
    ```

Now, when the closure returns the tuple, the first element is the `id` and `repeat` fields and the second element is the `type` field. We have effectively removed the `id` and `repeat` fields from the sample data and uniquely store them in the grouping key. This approach eliminates redundancy while maintaining all necessary information.

```bash title="View deduplicated data"
nextflow run main.nf
```

```console title="View deduplicated data"
 N E X T F L O W   ~  version 24.10.5

Launching `main.nf` [trusting_pike] DSL2 - revision: 09d3c7a81b

[[id:patientA, repeat:1], [type:normal], patientA_rep1_normal.bam, [type:tumor], patientA_rep1_tumor.bam]
[[id:patientA, repeat:2], [type:normal], patientA_rep2_normal.bam, [type:tumor], patientA_rep2_tumor.bam]
[[id:patientB, repeat:1], [type:normal], patientB_rep1_normal.bam, [type:tumor], patientB_rep1_tumor.bam]
[[id:patientC, repeat:1], [type:normal], patientC_rep1_normal.bam, [type:tumor], patientC_rep1_tumor.bam]
```

We can see we only state the `id` and `repeat` fields once in the grouping key and we have the `type` field in the sample data. We haven't lost any information but we managed to make our channel contents more succinct.

### 3.6. Remove redundant information

We removed duplicated information above, but we have some other redundant information in our channels.
In the beginning, we separated the normal and tumor samples using `filter`.
We then joined them based on `id` and `repeat` keys.
The `join` operator preserves the order in which a tuple is merged. In this example, we are using the normal on the left side and tumor on the right with the `id` as the join key.
The resulting channel preserves this order with `id, <elements normal>, <elements tumor>`.
Therefore, we know where in our current channel each element is.
We can simplify our channel structure further by dropping `[type:normal]` and `[type:tumor]`.

=== "After"

    ```groovy title="main.nf" linenums="5" hl_lines="3-4"
    getSampleIdAndReplicate = { meta, file -> [ meta.subMap(['id', 'repeat']), file ] }

    ```

=== "Before"

    ```groovy title="main.nf" linenums="5"
    getSampleIdAndReplicate = { meta, file -> [ meta.subMap(['id', 'repeat']), meta.subMap(['type']), file ] }

    ```

```bash title="Remove redundant information"
nextflow run main.nf
```

```console title="Remove redundant information"
 N E X T F L O W   ~  version 24.10.5

Launching `main.nf` [trusting_pike] DSL2 - revision: 09d3c7a81b

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
- **Deduplicating Data**: How to remove duplicate data from the channel

You now have a workflow that can split a samplesheet, filter the normal and tumor samples, join them together by sample ID and replicate number, then print the results.

This is a common pattern in bioinformatics workflows where you need to match up samples or other types of data after processing independently, so it is a useful skill. Next, we will look at repeating a sample multiple times.

## 4. Spread patients over intervals

A key pattern in bioinformatics workflows is distributing analysis across genomic regions. For instance, variant calling can be parallelized by dividing the genome into intervals (like chromosomes or smaller regions). This parallelization strategy significantly improves pipeline efficiency by distributing computational load across multiple cores or nodes, reducing overall execution time.

In the following section, we'll demonstrate how to distribute our sample data across multiple genomic intervals. We'll pair each sample with every interval, allowing parallel processing of different genomic regions. This will multiply our dataset size by the number of intervals, creating multiple independent analysis units that can be brought back together later.

### 4.1. Spread samples over intervals using `combine`

Let's start by creating a channel of intervals. To keep life simple, we will just use 3 intervals we will manually define. In a real workflow, you could read these in from a file input or even create a channel with lots of interval files.

=== "After"

    ```groovy title="main.nf" linenums="24" hl_lines="3"
            .join(ch_tumor_patients)
        ch_intervals = Channel.of('chr1', 'chr2', 'chr3')
    ```

=== "Before"

    ```groovy title="main.nf" linenums="24" hl_lines="2"
            .join(ch_tumor_patients)
        ch_joined_samples.view()
    ```

Now remember, we want to repeat each sample for each interval. This is sometimes referred to as the Cartesian product of the samples and intervals. We can achieve this by using the [`combine` operator](https://www.nextflow.io/docs/latest/operator.html#combine). This will take every item from channel 1 and repeat it for each item in channel 2. Let's add a combine operator to our workflow:

=== "After"

    ```groovy title="main.nf" linenums="26" hl_lines="3-4"
        ch_intervals = Channel.of('chr1', 'chr2', 'chr3')

        ch_combined_samples = ch_joined_samples
            .combine(ch_intervals)
            .view()
    ```

=== "Before"

    ```groovy title="main.nf" linenums="26"
        ch_intervals = Channel.of('chr1', 'chr2', 'chr3')
    ```

Now let's run it and see what happens:

```bash title="View combined samples"
nextflow run main.nf
```

```console title="View combined samples"
 N E X T F L O W   ~  version 24.10.5

Launching `main.nf` [soggy_fourier] DSL2 - revision: fa8f5edb22

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

    ```groovy title="main.nf" linenums="25" hl_lines="3-9"
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

    ```groovy title="main.nf" linenums="28"
        ch_combined_samples = ch_joined_samples
            .combine(ch_intervals)
            .view()
    ```

Wait? What did we do here? Let's go over it piece by piece.

First, we use a map operator to iterate over every item in the channel. By using the names `grouping_key` ,`normal`, `tumor` and `interval`, we can refer to the elements in the tuple by name instead of by index. This makes the code more readable and easier to understand.

```groovy
        .map { grouping_key, normal, tumor, interval ->
```

Next, create a new map by combining the `grouping_key` with the `interval` field. Remember, the `grouping_key` is the first element of the tuple, which is a map of `id` and `repeat` fields. The `interval` is just a string, but we make it into a new map with the key `interval` and value the string. By 'adding' them (`+`), Groovy will merge them together to produce the union of the two maps.

```groovy
                grouping_key + [interval: interval],
```

Finally, we return all of this as one tuple of the 3 elements, the new map, the normal sample data and the tumor sample data.

```groovy
            [
                grouping_key + [interval: interval],
                normal,
                tumor
            ]
```

Let's run it again and check the channel contents:

```bash title="View combined samples"
nextflow run main.nf
```

```console title="View combined samples"
 N E X T F L O W   ~  version 24.10.5

Launching `main.nf` [sad_hawking] DSL2 - revision: 1f6f6250cd

[[id:patientA, repeat:1, interval:chr1], patientA_rep1_normal.bam, patientA_rep1_tumor.bam]
[[id:patientA, repeat:1, interval:chr2], patientA_rep1_normal.bam, patientA_rep1_tumor.bam]
[[id:patientA, repeat:1, interval:chr3], patientA_rep1_normal.bam, patientA_rep1_tumor.bam]
[[id:patientA, repeat:2, interval:chr2], patientA_rep2_normal.bam, patientA_rep2_tumor.bam]
[[id:patientA, repeat:2, interval:chr1], patientA_rep2_normal.bam, patientA_rep2_tumor.bam]
[[id:patientA, repeat:2, interval:chr3], patientA_rep2_normal.bam, patientA_rep2_tumor.bam]
[[id:patientB, repeat:1, interval:chr1], patientB_rep1_normal.bam, patientB_rep1_tumor.bam]
[[id:patientB, repeat:1, interval:chr2], patientB_rep1_normal.bam, patientB_rep1_tumor.bam]
[[id:patientB, repeat:1, interval:chr3], patientB_rep1_normal.bam, patientB_rep1_tumor.bam]
[[id:patientC, repeat:1, interval:chr1], patientC_rep1_normal.bam, patientC_rep1_tumor.bam]
[[id:patientC, repeat:1, interval:chr2], patientC_rep1_normal.bam, patientC_rep1_tumor.bam]
[[id:patientC, repeat:1, interval:chr3], patientC_rep1_normal.bam, patientC_rep1_tumor.bam]
```

Using `map` to coerce your data into the correct structure can be tricky, but it's crucial to correctly splitting and grouping effectively.

### Takeaway

In this section, you've learned:

- **Spreading samples over intervals**: How to use `combine` to repeat samples over intervals

## 5. Aggregating samples using `groupTuple`

In the previous section, we learned how to split data from an input file and filter by specific fields (in our case normal and tumor samples). But this only covers a single type of joining. What if we want to group samples by a specific attribute? For example, instead of joining matched normal-tumor pairs, we might want to process all samples from "sampleA" together regardless of their type. This pattern is common in bioinformatics workflows where you may want to process related samples separately for efficiency reasons before comparing or combining the results at the end.

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

    ```groovy title="main.nf" linenums="28" hl_lines="11-20"
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

        ch_grouped_samples = ch_combined_patients
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

    ```groovy title="main.nf" linenums="25"
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

```bash title="View grouped samples"
nextflow run main.nf
```

```console title="View grouped samples"
 N E X T F L O W   ~  version 24.10.5

Launching `main.nf` [loving_escher] DSL2 - revision: 3adccba898

[[id:sampleA, interval:chr1], sampleA_rep1_normal.bam, sampleA_rep1_tumor.bam]
[[id:sampleA, interval:chr2], sampleA_rep1_normal.bam, sampleA_rep1_tumor.bam]
[[id:sampleA, interval:chr3], sampleA_rep1_normal.bam, sampleA_rep1_tumor.bam]
[[id:sampleA, interval:chr1], sampleA_rep2_normal.bam, sampleA_rep2_tumor.bam]
[[id:sampleA, interval:chr2], sampleA_rep2_normal.bam, sampleA_rep2_tumor.bam]
[[id:sampleA, interval:chr3], sampleA_rep2_normal.bam, sampleA_rep2_tumor.bam]
[[id:sampleB, interval:chr1], sampleB_rep1_normal.bam, sampleB_rep1_tumor.bam]
[[id:sampleB, interval:chr2], sampleB_rep1_normal.bam, sampleB_rep1_tumor.bam]
[[id:sampleB, interval:chr3], sampleB_rep1_normal.bam, sampleB_rep1_tumor.bam]
[[id:sampleC, interval:chr1], sampleC_rep1_normal.bam, sampleC_rep1_tumor.bam]
[[id:sampleC, interval:chr2], sampleC_rep1_normal.bam, sampleC_rep1_tumor.bam]
[[id:sampleC, interval:chr3], sampleC_rep1_normal.bam, sampleC_rep1_tumor.bam]
```

We can see that we have successfully isolated the `id` and `interval` fields, but not grouped the samples yet.

!!! note

    We are discarding the `replicate` field here. This is because we don't need it for further downstream processing. After completing this tutorial, see if you can include it without affecting the later grouping!

Let's now group the samples by this new grouping element, using the [`groupTuple` operator](https://www.nextflow.io/docs/latest/operator.html#grouptuple).

=== "After"

    ```groovy title="main.nf" linenums="39" hl_lines="8"
        ch_grouped_samples = ch_combined_patients
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

    ```groovy title="main.nf" linenums="39"
        ch_grouped_samples = ch_combined_patients
            .map { grouping_key, normal, tumor ->
                [
                    grouping_key.subMap('id', 'interval'),
                    normal,
                    tumor
                ]
            }
            .view()
    ```

Simple, huh? We just added a single line of code. Let's see what happens when we run it:

```bash title="View grouped samples"
nextflow run main.nf
```

```console title="View grouped samples"
 N E X T F L O W   ~  version 24.10.5

Launching `main.nf` [festering_almeida] DSL2 - revision: 78988949e3

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

Note our data has changed structure. What was previously a list of tuples is now a list of lists of tuples. This is because when we use `groupTuple`, Nextflow creates a new list for each group. This is important to remember when trying to handle the data downstream.

!!! note

    [`transpose`](https://www.nextflow.io/docs/latest/reference/operator.html#transpose) is the opposite of groupTuple. It unpacks the items in a channel and flattens them. Try and add `transpose` and undo the grouping we performed above!

## Summary

In this side quest, you've learned how to split and group data using channels. By modifying the data as it flows through the pipeline, you can construct a pipeline that handles as many items as possible with no loops or while statements. It gracefully scales to large numbers of items. Here's what we achieved:

1. **Read in samplesheet with splitCsv**: We read in a CSV file with sample data and viewed the contents.

2. **Use filter (and/or map) to manipulate into 2 separate channels**: We used `filter` to split the data into two channels based on the `type` field.

3. **Join on ID**: We used `join` to join the two channels on the `id` field.

4. **Use groupTuple to group up samples by ID**: We used `groupTuple` to group the samples by the `id` field.

5. **Combine by intervals**: We used `combine` to combine the two channels on the `interval` field.

6. **Group after intervals**: We used `groupTuple` to group the samples by the `interval` field.

This approach offers several advantages over writing a pipeline as more standard code, such as using for and while loops:

- We can scale to as many or as few inputs as we want with no additional code
- We focus on handling the flow of data through the pipeline, instead of iteration
- We can be as complex or simple as required
- The pipeline becomes more declarative, focusing on what should happen rather than how it should happen
- Nextflow will optimize execution for us by running independent operations in parallel

By mastering these channel operations, you can build flexible, scalable pipelines that handle complex data relationships without resorting to loops or iterative programming. This declarative approach allows Nextflow to optimize execution and parallelize independent operations automatically.

### Key Concepts

- **Reading data sheets**

  ```nextflow
  // Read CSV with header
  Channel.fromPath('samplesheet.csv')
      .splitCsv(header: true)
  ```

- **Filtering**

  ```nextflow
  // Filter channel based on condition
  channel.filter { it.type == 'tumor' }
  ```

- **Joining Channels**

  ```nextflow
  // Join two channels by key
  tumor_ch.join(normal_ch)

  // Extract a key and join by this value
  tumor_ch.map { [it.patient_id, it] }
      .join(
         normal_ch.map { [it.patient_id, it] }
       )
  ```

- **Grouping Data**

  ```nextflow
  // Group by the first element in each tuple
  channel.groupTuple()
  ```

- **Combining Channels**

  ```nextflow
  // Combine with Cartesian product
  samples_ch.combine(intervals_ch)
  ```

## Resources

- [filter](https://www.nextflow.io/docs/latest/operator.html#filter)
- [map](https://www.nextflow.io/docs/latest/operator.html#map)
- [join](https://www.nextflow.io/docs/latest/operator.html#join)
- [groupTuple](https://www.nextflow.io/docs/latest/operator.html#grouptuple)
- [combine](https://www.nextflow.io/docs/latest/operator.html#combine)
