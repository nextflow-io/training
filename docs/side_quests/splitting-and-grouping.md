# Splitting and Grouping

Nextflow helps you work with your data in flexible ways. One of the most useful things you can do is split your data into different streams and then group related items back together.

Think of it like sorting mail: you might first separate letters by their destination, process each pile differently, and then recombine items going to the same person. In Nextflow, we use special operators to do this with our scientific data.

Nextflow's channel system is at the heart of this flexibility. Channels act as pipelines that connect different parts of your workflow, allowing data to flow through your analysis. You can create multiple channels from a single data source, process each channel differently, and then merge channels back together when needed. This approach lets you design workflows that naturally mirror the branching and converging paths of complex bioinformatics analyses.

In this side quest, we'll explore how to split and group data using Nextflow's powerful channel operators. We'll start with a samplesheet containing information about different samples and their associated data. By the end of this side quest, you'll be able to manipulate and combine data streams effectively, making your workflows more efficient and easier to understand.

- Read data from files using `splitCsv`
- Filter and transform data with `filter` and `map`
- Combine related data using `join` and `groupTuple`

These skills will help you build workflows that can handle multiple samples and different types of data efficiently.

---

## 0. Warmup

### 0.1 Prerequisites

Before taking on this side quest you should:

- Complete the [Hello Nextflow](../hello_nextflow/README.md) tutorial
- Understand basic Nextflow concepts (processes, channels, operators)

### 0.2 Starting Point

Let's move into the project directory.

```bash
cd side-quests/splitting-and-grouping
```

You'll find a `data` directory containing a samplesheet and a main workflow file.

```console title="Directory contents"
> tree
.
├── data
│   └── samplesheet.csv
└── main.nf
```

The samplesheet contains information about different samples and their associated data. In particular, it contains information about the sample's ID, repeat number, type (normal or tumor), and the paths to the fastq files.

```console title="samplesheet.csv"
id,repeat,type,fastq1,fastq2
sampleA,1,normal,sampleA_rep1_normal_R1.fastq.gz,sampleA_rep1_normal_R2.fastq.gz
sampleA,1,tumor,sampleA_rep1_tumor_R1.fastq.gz,sampleA_rep1_tumor_R2.fastq.gz
sampleA,2,normal,sampleA_rep2_normal_R1.fastq.gz,sampleA_rep2_normal_R2.fastq.gz
sampleA,2,tumor,sampleA_rep2_tumor_R1.fastq.gz,sampleA_rep2_tumor_R2.fastq.gz
sampleB,1,normal,sampleB_rep1_normal_R1.fastq.gz,sampleB_rep1_normal_R2.fastq.gz
sampleB,1,tumor,sampleB_rep1_tumor_R1.fastq.gz,sampleB_rep1_tumor_R2.fastq.gz
sampleC,1,normal,sampleC_rep1_normal_R1.fastq.gz,sampleC_rep1_normal_R2.fastq.gz
sampleC,1,tumor,sampleC_rep1_tumor_R1.fastq.gz,sampleC_rep1_tumor_R2.fastq.gz
```

Note there are 8 samples in total, 4 normal and 4 tumor. sampleA has 2 repeats, while sampleB and sampleC only have 1.

We're going to read in this samplesheet, then group and split the samples based on their data.

---

## 1. Read in samplesheet

### 1.1. Read in samplesheet with splitCsv

Let's start by reading in the samplesheet with `splitCsv`. In the main workflow file, you'll see that we've already started the workflow.

```groovy title="main.nf" linenums="1"
workflow {
    ch_samplesheet = Channel.fromPath("./data/samplesheet.csv")
}
```

!!! note
Throughout this tutorial, we'll use the `ch_` prefix for all channel variables to clearly indicate they are Nextflow channels.

We can use the [`splitCsv` operator](https://www.nextflow.io/docs/latest/operator.html#splitcsv) to split the samplesheet into a channel of maps, where each map represents a row from the CSV file.

_Before:_

```groovy title="main.nf" linenums="1"
workflow {
    ch_samplesheet = Channel.fromPath("./data/samplesheet.csv")
                        .splitCsv(header: true)
}
```

The `header: true` option tells Nextflow to use the first row of the CSV file as the header row, which will be used as keys for the values. Let's see what Nextflow can see after reading with splitCsv. To do this, we can use the `view` operator.

_After:_

```groovy title="main.nf" linenums="1"
workflow {
    ch_samplesheet = Channel.fromPath("./data/samplesheet.csv")
                        .splitCsv(header: true)
                        .view()
}
```

```bash title="Read the samplesheet"
nextflow run main.nf
```

```console title="Read samplesheet with splitCsv"
 N E X T F L O W   ~  version 24.10.5

Launching `main.nf` [berserk_cray] DSL2 - revision: 8f31622c03

[id:sampleA, repeat:1, type:normal, fastq1:sampleA_rep1_normal_R1.fastq.gz, fastq2:sampleA_rep1_normal_R2.fastq.gz]
[id:sampleA, repeat:1, type:tumor, fastq1:sampleA_rep1_tumor_R1.fastq.gz, fastq2:sampleA_rep1_tumor_R2.fastq.gz]
[id:sampleA, repeat:2, type:normal, fastq1:sampleA_rep2_normal_R1.fastq.gz, fastq2:sampleA_rep2_normal_R2.fastq.gz]
[id:sampleA, repeat:2, type:tumor, fastq1:sampleA_rep2_tumor_R1.fastq.gz, fastq2:sampleA_rep2_tumor_R2.fastq.gz]
[id:sampleB, repeat:1, type:normal, fastq1:sampleB_rep1_normal_R1.fastq.gz, fastq2:sampleB_rep1_normal_R2.fastq.gz]
[id:sampleB, repeat:1, type:tumor, fastq1:sampleB_rep1_tumor_R1.fastq.gz, fastq2:sampleB_rep1_tumor_R2.fastq.gz]
[id:sampleC, repeat:1, type:normal, fastq1:sampleC_rep1_normal_R1.fastq.gz, fastq2:sampleC_rep1_normal_R2.fastq.gz]
[id:sampleC, repeat:1, type:tumor, fastq1:sampleC_rep1_tumor_R1.fastq.gz, fastq2:sampleC_rep1_tumor_R2.fastq.gz]
```

We can see that each row from the CSV file has been converted into a map with keys matching the header row. A map is a key-value data structure similar to dictionaries in Python, objects in JavaScript, or hashes in Ruby.

Each map contains:

- `id`: The sample identifier (sampleA, sampleB, sampleC)
- `repeat`: The replicate number (1 or 2)
- `type`: The sample type (normal or tumor)
- `fastq1`: Path to the first FASTQ file
- `fastq2`: Path to the second FASTQ file

This format makes it easy to access specific fields from each sample. For example, we could access the sample ID with `row.id` or the FASTQ paths with `row.fastq1` and `row.fastq2`.

This means we have successfully read in the samplesheet and have access to the data in each row. We can start to implement this in our pipeline.

### 1.2. Use dump to pretty print the data

For a prettier output format, we can use the [`dump` operator](https://www.nextflow.io/docs/latest/operator.html#dump) instead of `view`:

```groovy title="main.nf" linenums="1"
workflow {
    ch_samplesheet = Channel.fromPath("./data/samplesheet.csv")
                        .splitCsv(header: true)
                        .dump(tag: 'samples', pretty: true)
}
```

```bash title="Read the samplesheet"
nextflow run main.nf
```

```console title="Read samplesheet with dump"
 N E X T F L O W   ~  version 24.10.5

Launching `./main.nf` [grave_stone] DSL2 - revision: b2bafa8755
```

Wait?! Where is our output? `dump` is a special operator that prints the data to the console only when specifically enabled. That is what the `tag` parameter is for. Let's enable it:

```bash title="Enable dump"
nextflow run main.nf -dump-channels samples
```

```console title="Read samplesheet with dump"
 N E X T F L O W   ~  version 24.10.5

Launching `main.nf` [wise_kirch] DSL2 - revision: 7f194f2473

[DUMP: samples] {
    "id": "sampleA",
    "repeat": "1",
    "type": "normal",
    "fastq1": "sampleA_rep1_normal_R1.fastq.gz",
    "fastq2": "sampleA_rep1_normal_R2.fastq.gz"
}
[DUMP: samples] {
    "id": "sampleA",
    "repeat": "1",
    "type": "tumor",
    "fastq1": "sampleA_rep1_tumor_R1.fastq.gz",
    "fastq2": "sampleA_rep1_tumor_R2.fastq.gz"
}
[DUMP: samples] {
    "id": "sampleA",
    "repeat": "2",
    "type": "normal",
    "fastq1": "sampleA_rep2_normal_R1.fastq.gz",
    "fastq2": "sampleA_rep2_normal_R2.fastq.gz"
}
[DUMP: samples] {
    "id": "sampleA",
    "repeat": "2",
    "type": "tumor",
    "fastq1": "sampleA_rep2_tumor_R1.fastq.gz",
    "fastq2": "sampleA_rep2_tumor_R2.fastq.gz"
}
[DUMP: samples] {
    "id": "sampleB",
    "repeat": "1",
    "type": "normal",
    "fastq1": "sampleB_rep1_normal_R1.fastq.gz",
    "fastq2": "sampleB_rep1_normal_R2.fastq.gz"
}
[DUMP: samples] {
    "id": "sampleB",
    "repeat": "1",
    "type": "tumor",
    "fastq1": "sampleB_rep1_tumor_R1.fastq.gz",
    "fastq2": "sampleB_rep1_tumor_R2.fastq.gz"
}
[DUMP: samples] {
    "id": "sampleC",
    "repeat": "1",
    "type": "normal",
    "fastq1": "sampleC_rep1_normal_R1.fastq.gz",
    "fastq2": "sampleC_rep1_normal_R2.fastq.gz"
}
[DUMP: samples] {
    "id": "sampleC",
    "repeat": "1",
    "type": "tumor",
    "fastq1": "sampleC_rep1_tumor_R1.fastq.gz",
    "fastq2": "sampleC_rep1_tumor_R2.fastq.gz"
}
```

This is a long output, but we can see that each row from the CSV file has been converted into a map with keys matching the header row. It's more clear to read, at the cost of being too much content for a small terminal. If you want it to be more concise, you can remove the `pretty: true` parameter and the console output will be similar to `view`.

!!! note
If the output is too tall for your terminal but you have a very wide terminal, you can remove `pretty: true` from the `dump` operator to make it more concise.

Both dump and view are useful for debugging and we will continue to use them throughout this side quest. Feel free to intersperse them if you need additional clarification at any step.

### Takeaway

In this section, you've learned:

- **Reading in a samplesheet**: How to read in a samplesheet with `splitCsv`
- **Viewing data**: How to use `view` to print the data
- **Dumping data**: How to use `dump` to pretty print the data

We now have a channel of maps, each representing a row from the samplesheet. Next, we'll transform this data into a format suitable for our pipeline by extracting metadata and organizing the file paths.

---

## 2. Filter and transform data

### 2.1. Filter data with `filter`

We can use the [`filter` operator](https://www.nextflow.io/docs/latest/operator.html#filter) to filter the data based on a condition. Let's say we only want to process normal samples. We can do this by filtering the data based on the `type` field. Let's insert this before the `dump` operator.

_Before:_

```groovy title="main.nf" linenums="1"
workflow {
    ch_samplesheet = Channel.fromPath("./data/samplesheet.csv")
                        .splitCsv(header: true)
                        .dump(tag: 'samples', pretty: true)
}
```

_After:_

```groovy title="main.nf" linenums="1"
workflow {
    ch_samplesheet = Channel.fromPath("./data/samplesheet.csv")
                        .splitCsv(header: true)
                        .filter { sample -> sample.type == 'normal' }
                        .dump(tag: 'samples')
}
```

!!! note
We drop the `pretty: true` parameter from `dump` because it makes it easier to see the difference

```bash title="View normal samples"
nextflow run main.nf
```

```console title="View normal samples"
 N E X T F L O W   ~  version 24.10.5

Launching `main.nf` [stupefied_pike] DSL2 - revision: 8761d1b103

[DUMP: samples] ['id':'sampleA', 'repeat':'1', 'type':'normal', 'fastq1':'sampleA_rep1_normal_R1.fastq.gz', 'fastq2':'sampleA_rep1_normal_R2.fastq.gz']
[DUMP: samples] ['id':'sampleA', 'repeat':'2', 'type':'normal', 'fastq1':'sampleA_rep2_normal_R1.fastq.gz', 'fastq2':'sampleA_rep2_normal_R2.fastq.gz']
[DUMP: samples] ['id':'sampleB', 'repeat':'1', 'type':'normal', 'fastq1':'sampleB_rep1_normal_R1.fastq.gz', 'fastq2':'sampleB_rep1_normal_R2.fastq.gz']
[DUMP: samples] ['id':'sampleC', 'repeat':'1', 'type':'normal', 'fastq1':'sampleC_rep1_normal_R1.fastq.gz', 'fastq2':'sampleC_rep1_normal_R2.fastq.gz']
```

We have successfully filtered the data to only include normal samples. Let's recap how this works. The `filter` operator takes a closure that is applied to each element in the channel. If the closure returns `true`, the element is included in the output channel. If the closure returns `false`, the element is excluded from the output channel.

In this case, we want to keep only the samples where `sample.type == 'normal'`. In the closure, we use the variable name `sample` to refer to each element in the channel, which then checks if `sample.type` is equal to `'normal'`. If it is, the sample is included in the output channel. If it is not, the sample is excluded from the output channel.

```groovy title="main.nf" linenums="4"
.filter { sample -> sample.type == 'normal' }
```

### 2.2. Save results of filter to a new channel

While useful, we are discarding the tumor samples. Instead, let's rewrite our pipeline to save all the samples to one channel called `samplesheet`, then filter that channel to just the normal samples and save the results to a new channel called `normal_samples`.

_Before:_

```groovy title="main.nf" linenums="1"
workflow {
    ch_samplesheet = Channel.fromPath("./data/samplesheet.csv")
                        .splitCsv(header: true)
                        .filter { sample -> sample.type == 'normal' }
                        .dump(tag: 'samples')
}
```

_After:_

```groovy title="main.nf" linenums="1"
workflow {
    ch_samplesheet = Channel.fromPath("./data/samplesheet.csv")
                        .splitCsv(header: true)
    ch_normal_samples = ch_samplesheet
                        .filter { sample -> sample.type == 'normal' }
                        .view()
}
```

Once again, run the pipeline to see the results:

```bash title="View normal samples"
nextflow run main.nf
```

```console title="View normal samples"
 N E X T F L O W   ~  version 24.10.5

Launching `main.nf` [lonely_miescher] DSL2 - revision: 7e26f19fd3

[id:sampleA, repeat:1, type:normal, fastq1:sampleA_rep1_normal_R1.fastq.gz, fastq2:sampleA_rep1_normal_R2.fastq.gz]
[id:sampleA, repeat:2, type:normal, fastq1:sampleA_rep2_normal_R1.fastq.gz, fastq2:sampleA_rep2_normal_R2.fastq.gz]
[id:sampleB, repeat:1, type:normal, fastq1:sampleB_rep1_normal_R1.fastq.gz, fastq2:sampleB_rep1_normal_R2.fastq.gz]
[id:sampleC, repeat:1, type:normal, fastq1:sampleC_rep1_normal_R1.fastq.gz, fastq2:sampleC_rep1_normal_R2.fastq.gz]
```

Success! We have filtered the data to only include normal samples. If we wanted, we still have access to the tumor samples within the `samplesheet` channel. Since we managed it for the normal samples, let's do it for the tumor samples as well:

_Before:_

```groovy title="main.nf" linenums="1"
workflow {
    ch_samplesheet = Channel.fromPath("./data/samplesheet.csv")
                         .splitCsv(header: true)
    ch_normal_samples = ch_samplesheet
                         .filter { sample -> sample.type == 'normal' }
                         .view()
}
```

_After:_

```groovy title="main.nf" linenums="1"
workflow {
    ch_samplesheet = Channel.fromPath("./data/samplesheet.csv")
                        .splitCsv(header: true)
    ch_normal_samples = ch_samplesheet
                        .filter { sample -> sample.type == 'normal' }
                        .view()
    ch_tumor_samples = ch_samplesheet
                        .filter { sample -> sample.type == 'tumor' }
                        .view()
}
```

```bash title="View tumor samples"
nextflow run main.nf
```

```console title="View tumor samples"
 N E X T F L O W   ~  version 24.10.5

Launching `main.nf` [focused_kirch] DSL2 - revision: 87d6672658

[id:sampleA, repeat:1, type:normal, fastq1:sampleA_rep1_normal_R1.fastq.gz, fastq2:sampleA_rep1_normal_R2.fastq.gz]
[id:sampleA, repeat:1, type:tumor, fastq1:sampleA_rep1_tumor_R1.fastq.gz, fastq2:sampleA_rep1_tumor_R2.fastq.gz]
[id:sampleA, repeat:2, type:normal, fastq1:sampleA_rep2_normal_R1.fastq.gz, fastq2:sampleA_rep2_normal_R2.fastq.gz]
[id:sampleA, repeat:2, type:tumor, fastq1:sampleA_rep2_tumor_R1.fastq.gz, fastq2:sampleA_rep2_tumor_R2.fastq.gz]
[id:sampleB, repeat:1, type:normal, fastq1:sampleB_rep1_normal_R1.fastq.gz, fastq2:sampleB_rep1_normal_R2.fastq.gz]
[id:sampleB, repeat:1, type:tumor, fastq1:sampleB_rep1_tumor_R1.fastq.gz, fastq2:sampleB_rep1_tumor_R2.fastq.gz]
[id:sampleC, repeat:1, type:normal, fastq1:sampleC_rep1_normal_R1.fastq.gz, fastq2:sampleC_rep1_normal_R2.fastq.gz]
[id:sampleC, repeat:1, type:tumor, fastq1:sampleC_rep1_tumor_R1.fastq.gz, fastq2:sampleC_rep1_tumor_R2.fastq.gz]
```

We've managed to separate out the normal and tumor samples into two different channels but they're mixed up when we `view` them in the console! Here's where dump could be useful, because it can label the different channels with a tag.

_Before:_

```groovy title="main.nf" linenums="1"
workflow {
    ch_samplesheet = Channel.fromPath("./data/samplesheet.csv")
                        .splitCsv(header: true)
    ch_normal_samples = ch_samplesheet
                        .filter { sample -> sample.type == 'normal' }
                        .view()
    ch_tumor_samples = ch_samplesheet
                        .filter { sample -> sample.type == 'tumor' }
                        .view()
}
```

_After:_

```groovy title="main.nf" linenums="1"
workflow {
    ch_samplesheet = Channel.fromPath("./data/samplesheet.csv")
                        .splitCsv(header: true)
    ch_normal_samples = ch_samplesheet
                        .filter { sample -> sample.type == 'normal' }
                        .dump(tag: 'normal')
    ch_tumor_samples = ch_samplesheet
                        .filter { sample -> sample.type == "tumor" }
                        .dump(tag: 'tumor')
}
```

```bash title="View normal and tumor samples"
nextflow run main.nf -dump-channels normal,tumor
```

```console title="View normal and tumor samples"
 N E X T F L O W   ~  version 24.10.5

Launching `main.nf` [spontaneous_jones] DSL2 - revision: 0e794240ef

[DUMP: tumor] ['id':'sampleA', 'repeat':'1', 'type':'tumor', 'fastq1':'sampleA_rep1_tumor_R1.fastq.gz', 'fastq2':'sampleA_rep1_tumor_R2.fastq.gz']
[DUMP: normal] ['id':'sampleA', 'repeat':'1', 'type':'normal', 'fastq1':'sampleA_rep1_normal_R1.fastq.gz', 'fastq2':'sampleA_rep1_normal_R2.fastq.gz']
[DUMP: tumor] ['id':'sampleA', 'repeat':'2', 'type':'tumor', 'fastq1':'sampleA_rep2_tumor_R1.fastq.gz', 'fastq2':'sampleA_rep2_tumor_R2.fastq.gz']
[DUMP: normal] ['id':'sampleA', 'repeat':'2', 'type':'normal', 'fastq1':'sampleA_rep2_normal_R1.fastq.gz', 'fastq2':'sampleA_rep2_normal_R2.fastq.gz']
[DUMP: tumor] ['id':'sampleB', 'repeat':'1', 'type':'tumor', 'fastq1':'sampleB_rep1_tumor_R1.fastq.gz', 'fastq2':'sampleB_rep1_tumor_R2.fastq.gz']
[DUMP: normal] ['id':'sampleB', 'repeat':'1', 'type':'normal', 'fastq1':'sampleB_rep1_normal_R1.fastq.gz', 'fastq2':'sampleB_rep1_normal_R2.fastq.gz']
[DUMP: tumor] ['id':'sampleC', 'repeat':'1', 'type':'tumor', 'fastq1':'sampleC_rep1_tumor_R1.fastq.gz', 'fastq2':'sampleC_rep1_tumor_R2.fastq.gz']
[DUMP: normal] ['id':'sampleC', 'repeat':'1', 'type':'normal', 'fastq1':'sampleC_rep1_normal_R1.fastq.gz', 'fastq2':'sampleC_rep1_normal_R2.fastq.gz']
```

Note how the `normal` and `tumor` tags are used to label the different channels. This is useful for debugging and for understanding the data flow in our pipeline.

### Takeaway

In this section, you've learned:

- **Filtering data**: How to filter data with `filter`
- **Splitting data**: How to split data into different channels based on a condition
- **Dumping data**: How to use `dump` to label and print the data

We've now separated out the normal and tumor samples into two different channels. Next, we'll join the normal and tumor samples on the `id` field.

---

## 3. Join on ID

In the previous section, we separated out the normal and tumor samples into two different channels. These could be processed independently using specific processes or workflows based on their type. But what happens when we want to compare the normal and tumor samples from the same patient? At this point, we need to join them back together making sure to match the samples based on their `id` field.

Nextflow includes many methods for combing channels, but in this case the most appropriate operator is [`join`](https://www.nextflow.io/docs/latest/operator.html#join). This acts like a SQL `JOIN` operation, where we specify the key to join on and the type of join to perform.

### 3.1. Use `map` and `join` to combine based on sample ID

If we check the [`join`](https://www.nextflow.io/docs/latest/operator.html#join) documentation, we can see that it joins two channels based on the first item in each tuple. Let's run the pipeline to check our data structure and see how we need to modify it to join on the `id` field.

```bash title="View normal and tumor samples"
nextflow run main.nf -dump-channels normal,tumor
```

```console title="View normal and tumor samples"
 N E X T F L O W   ~  version 24.10.5

Launching `main.nf` [spontaneous_jones] DSL2 - revision: 0e794240ef

[DUMP: tumor] ['id':'sampleA', 'repeat':'1', 'type':'tumor', 'fastq1':'sampleA_rep1_tumor_R1.fastq.gz', 'fastq2':'sampleA_rep1_tumor_R2.fastq.gz']
[DUMP: normal] ['id':'sampleA', 'repeat':'1', 'type':'normal', 'fastq1':'sampleA_rep1_normal_R1.fastq.gz', 'fastq2':'sampleA_rep1_normal_R2.fastq.gz']
[DUMP: tumor] ['id':'sampleA', 'repeat':'2', 'type':'tumor', 'fastq1':'sampleA_rep2_tumor_R1.fastq.gz', 'fastq2':'sampleA_rep2_tumor_R2.fastq.gz']
[DUMP: normal] ['id':'sampleA', 'repeat':'2', 'type':'normal', 'fastq1':'sampleA_rep2_normal_R1.fastq.gz', 'fastq2':'sampleA_rep2_normal_R2.fastq.gz']
[DUMP: tumor] ['id':'sampleB', 'repeat':'1', 'type':'tumor', 'fastq1':'sampleB_rep1_tumor_R1.fastq.gz', 'fastq2':'sampleB_rep1_tumor_R2.fastq.gz']
[DUMP: normal] ['id':'sampleB', 'repeat':'1', 'type':'normal', 'fastq1':'sampleB_rep1_normal_R1.fastq.gz', 'fastq2':'sampleB_rep1_normal_R2.fastq.gz']
[DUMP: tumor] ['id':'sampleC', 'repeat':'1', 'type':'tumor', 'fastq1':'sampleC_rep1_tumor_R1.fastq.gz', 'fastq2':'sampleC_rep1_tumor_R2.fastq.gz']
[DUMP: normal] ['id':'sampleC', 'repeat':'1', 'type':'normal', 'fastq1':'sampleC_rep1_normal_R1.fastq.gz', 'fastq2':'sampleC_rep1_normal_R2.fastq.gz']
```

We can see that the `id` field is the first element in each map. For `join` to work, we should isolate the `id` field in each tuple. After that, we can simply use the `join` operator to combine the two channels.

To isolate the `id` field, we can use the [`map` operator](https://www.nextflow.io/docs/latest/operator.html#map) to create a new tuple with the `id` field as the first element.

_Before:_

```groovy title="main.nf" linenums="1"
workflow {
    ch_samplesheet = Channel.fromPath("./data/samplesheet.csv")
                        .splitCsv(header: true)
    ch_normal_samples = ch_samplesheet
                        .filter { sample -> sample.type == 'normal' }
                        .dump(tag: 'normal')
    ch_tumor_samples = ch_samplesheet
                        .filter { sample -> sample.type == "tumor" }
                        .dump(tag: 'tumor')
}
```

_After:_

```groovy title="main.nf" linenums="1"
workflow {
    ch_samplesheet = Channel.fromPath("./data/samplesheet.csv")
                        .splitCsv(header: true)
    ch_normal_samples = ch_samplesheet
                        .filter { sample -> sample.type == 'normal' }
                        .map { sample -> [sample.id, sample] }
                        .dump(tag: 'normal')
    ch_tumor_samples = ch_samplesheet
                        .filter { sample -> sample.type == "tumor" }
                        .map { sample -> [sample.id, sample] }
                        .dump(tag: 'tumor')
}
```

```bash title="View normal and tumor samples with ID as element 0"
nextflow run main.nf -dump-channels normal,tumor
```

```console title="View normal and tumor samples with ID as element 0"
 N E X T F L O W   ~  version 24.10.5

Launching `main.nf` [sick_jones] DSL2 - revision: 9b183fbc7c

[DUMP: tumor] ['sampleA', ['id':'sampleA', 'repeat':'1', 'type':'tumor', 'fastq1':'sampleA_rep1_tumor_R1.fastq.gz', 'fastq2':'sampleA_rep1_tumor_R2.fastq.gz']]
[DUMP: normal] ['sampleA', ['id':'sampleA', 'repeat':'1', 'type':'normal', 'fastq1':'sampleA_rep1_normal_R1.fastq.gz', 'fastq2':'sampleA_rep1_normal_R2.fastq.gz']]
[DUMP: tumor] ['sampleA', ['id':'sampleA', 'repeat':'2', 'type':'tumor', 'fastq1':'sampleA_rep2_tumor_R1.fastq.gz', 'fastq2':'sampleA_rep2_tumor_R2.fastq.gz']]
[DUMP: normal] ['sampleA', ['id':'sampleA', 'repeat':'2', 'type':'normal', 'fastq1':'sampleA_rep2_normal_R1.fastq.gz', 'fastq2':'sampleA_rep2_normal_R2.fastq.gz']]
[DUMP: tumor] ['sampleB', ['id':'sampleB', 'repeat':'1', 'type':'tumor', 'fastq1':'sampleB_rep1_tumor_R1.fastq.gz', 'fastq2':'sampleB_rep1_tumor_R2.fastq.gz']]
[DUMP: normal] ['sampleB', ['id':'sampleB', 'repeat':'1', 'type':'normal', 'fastq1':'sampleB_rep1_normal_R1.fastq.gz', 'fastq2':'sampleB_rep1_normal_R2.fastq.gz']]
[DUMP: tumor] ['sampleC', ['id':'sampleC', 'repeat':'1', 'type':'tumor', 'fastq1':'sampleC_rep1_tumor_R1.fastq.gz', 'fastq2':'sampleC_rep1_tumor_R2.fastq.gz']]
[DUMP: normal] ['sampleC', ['id':'sampleC', 'repeat':'1', 'type':'normal', 'fastq1':'sampleC_rep1_normal_R1.fastq.gz', 'fastq2':'sampleC_rep1_normal_R2.fastq.gz']]
```

It might be subtle, but you should be able to see the first element in each tuple is the `id` field. Now we can use the `join` operator to combine the two channels based on the `id` field.

Once again, we will use `dump` to selectively print the joined outputs.

_Before:_

```groovy title="main.nf" linenums="1"
workflow {
    samplesheet = Channel.fromPath("./data/samplesheet.csv")
                        .splitCsv(header: true)
    ch_normal_samples = ch_samplesheet
                        .filter { sample -> sample.type == 'normal' }
                        .map { sample -> [sample.id, sample] }
                        .dump(tag: 'normal')
    ch_tumor_samples = ch_samplesheet
                        .filter { sample -> sample.type == "tumor" }
                        .map { sample -> [sample.id, sample] }
                        .dump(tag: 'tumor')
}
```

_After:_

```groovy title="main.nf" linenums="1"
workflow {
    samplesheet = Channel.fromPath("./data/samplesheet.csv")
                        .splitCsv(header: true)
    ch_normal_samples = ch_samplesheet
                        .filter { sample -> sample.type == 'normal' }
                        .map { sample -> [sample.id, sample] }
                        .dump(tag: 'normal')
    ch_tumor_samples = ch_samplesheet
                        .filter { sample -> sample.type == "tumor" }
                        .map { sample -> [sample.id, sample] }
                        .dump(tag: 'tumor')
    joined_samples = ch_normal_samples
                        .join(ch_tumor_samples)
                        .dump(tag: 'joined')
}
```

```bash title="View normal and tumor samples"
nextflow run main.nf -dump-channels joined
```

```console title="View joined normal and tumor samples"
 N E X T F L O W   ~  version 24.10.5

Launching `main.nf` [thirsty_poitras] DSL2 - revision: 95a2b8902b

[DUMP: joined] ['sampleA', ['id':'sampleA', 'repeat':'1', 'type':'normal', 'fastq1':'sampleA_rep1_normal_R1.fastq.gz', 'fastq2':'sampleA_rep1_normal_R2.fastq.gz'], ['id':'sampleA', 'repeat':'1', 'type':'tumor', 'fastq1':'sampleA_rep1_tumor_R1.fastq.gz', 'fastq2':'sampleA_rep1_tumor_R2.fastq.gz']]
[DUMP: joined] ['sampleA', ['id':'sampleA', 'repeat':'2', 'type':'normal', 'fastq1':'sampleA_rep2_normal_R1.fastq.gz', 'fastq2':'sampleA_rep2_normal_R2.fastq.gz'], ['id':'sampleA', 'repeat':'2', 'type':'tumor', 'fastq1':'sampleA_rep2_tumor_R1.fastq.gz', 'fastq2':'sampleA_rep2_tumor_R2.fastq.gz']]
[DUMP: joined] ['sampleB', ['id':'sampleB', 'repeat':'1', 'type':'normal', 'fastq1':'sampleB_rep1_normal_R1.fastq.gz', 'fastq2':'sampleB_rep1_normal_R2.fastq.gz'], ['id':'sampleB', 'repeat':'1', 'type':'tumor', 'fastq1':'sampleB_rep1_tumor_R1.fastq.gz', 'fastq2':'sampleB_rep1_tumor_R2.fastq.gz']]
[DUMP: joined] ['sampleC', ['id':'sampleC', 'repeat':'1', 'type':'normal', 'fastq1':'sampleC_rep1_normal_R1.fastq.gz', 'fastq2':'sampleC_rep1_normal_R2.fastq.gz'], ['id':'sampleC', 'repeat':'1', 'type':'tumor', 'fastq1':'sampleC_rep1_tumor_R1.fastq.gz', 'fastq2':'sampleC_rep1_tumor_R2.fastq.gz']]
```

It's a little hard to tell because it's so wide, but you should be able to see the samples have been joined by the `id` field. Each tuple now has the format:

- `id`: The sample ID
- `normal_sample`: The normal sample including type, replicate and path to fastq files
- `tumor_sample`: The tumor sample including type, replicate and path to fastq files

If you want you can use the `pretty` parameter of `dump` to make it easier to read:

_After:_

```groovy title="main.nf" linenums="1"
workflow {
    ch_samplesheet = Channel.fromPath("./data/samplesheet.csv")
                        .splitCsv(header: true)
    ch_normal_samples = ch_samplesheet
                        .filter { sample -> sample.type == 'normal' }
                        .map { sample -> [sample.id, sample] }
                        .dump(tag: 'normal')
    ch_tumor_samples = ch_samplesheet
                        .filter { sample -> sample.type == "tumor" }
                        .map { sample -> [sample.id, sample] }
                        .dump(tag: 'tumor')
    joined_samples = ch_normal_samples
                        .join(ch_tumor_samples)
                        .dump(tag: 'joined', pretty: true)
}
```

```bash title="View normal and tumor samples"
nextflow run main.nf -dump-channels joined
```

```console title="View normal and tumor samples"
 N E X T F L O W   ~  version 24.10.5

Launching `main.nf` [tender_feynman] DSL2 - revision: 3505c6a732

[DUMP: joined] [
    "sampleA",
    {
        "id": "sampleA",
        "repeat": "1",
        "type": "normal",
        "fastq1": "sampleA_rep1_normal_R1.fastq.gz",
        "fastq2": "sampleA_rep1_normal_R2.fastq.gz"
    },
    {
        "id": "sampleA",
        "repeat": "1",
        "type": "tumor",
        "fastq1": "sampleA_rep1_tumor_R1.fastq.gz",
        "fastq2": "sampleA_rep1_tumor_R2.fastq.gz"
    }
]
[DUMP: joined] [
    "sampleA",
    {
        "id": "sampleA",
        "repeat": "2",
        "type": "normal",
        "fastq1": "sampleA_rep2_normal_R1.fastq.gz",
        "fastq2": "sampleA_rep2_normal_R2.fastq.gz"
    },
    {
        "id": "sampleA",
        "repeat": "2",
        "type": "tumor",
        "fastq1": "sampleA_rep2_tumor_R1.fastq.gz",
        "fastq2": "sampleA_rep2_tumor_R2.fastq.gz"
    }
]
[DUMP: joined] [
    "sampleB",
    {
        "id": "sampleB",
        "repeat": "1",
        "type": "normal",
        "fastq1": "sampleB_rep1_normal_R1.fastq.gz",
        "fastq2": "sampleB_rep1_normal_R2.fastq.gz"
    },
    {
        "id": "sampleB",
        "repeat": "1",
        "type": "tumor",
        "fastq1": "sampleB_rep1_tumor_R1.fastq.gz",
        "fastq2": "sampleB_rep1_tumor_R2.fastq.gz"
    }
]
[DUMP: joined] [
    "sampleC",
    {
        "id": "sampleC",
        "repeat": "1",
        "type": "normal",
        "fastq1": "sampleC_rep1_normal_R1.fastq.gz",
        "fastq2": "sampleC_rep1_normal_R2.fastq.gz"
    },
    {
        "id": "sampleC",
        "repeat": "1",
        "type": "tumor",
        "fastq1": "sampleC_rep1_tumor_R1.fastq.gz",
        "fastq2": "sampleC_rep1_tumor_R2.fastq.gz"
    }
]
```

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

_Before:_

```groovy title="main.nf" linenums="1"
workflow {
    ch_samplesheet = Channel.fromPath("./data/samplesheet.csv")
                        .splitCsv(header: true)
    ch_normal_samples = ch_samplesheet
                        .filter { sample -> sample.type == 'normal' }
                        .map { sample -> [sample.id, sample] }
                        .dump(tag: 'normal')
    ch_tumor_samples = ch_samplesheet
                        .filter { sample -> sample.type == "tumor" }
                        .map { sample -> [sample.id, sample] }
                        .dump(tag: 'tumor')
    ch_joined_samples = ch_normal_samples
                        .join(ch_tumor_samples)
                        .dump(tag: 'joined', pretty: true)
}
```

_After:_

```groovy title="main.nf" linenums="1"
workflow {
    ch_samplesheet = Channel.fromPath("./data/samplesheet.csv")
                        .splitCsv(header: true)
    ch_normal_samples = ch_samplesheet
                        .filter { sample -> sample.type == 'normal' }
                        .map { sample -> [
                                [sample.id, sample.repeat],
                                sample
                            ]
                        }
                        .dump(tag: 'normal')
    ch_tumor_samples = ch_samplesheet
                        .filter { sample -> sample.type == "tumor" }
                        .map { sample -> [
                                [sample.id, sample.repeat],
                                sample
                            ]
                        }
                        .dump(tag: 'tumor')
    ch_joined_samples = ch_normal_samples
                        .join(ch_tumor_samples)
                        .dump(tag: 'joined', pretty: true)
}
```

Now we should see the join is occurring but using both the `id` and `repeat` fields.

```bash title="View normal and tumor samples"
nextflow run main.nf -dump-channels joined
```

```console title="View normal and tumor samples"
 N E X T F L O W   ~  version 24.10.5

Launching `main.nf` [cranky_lorenz] DSL2 - revision: 2be25de1df

[DUMP: joined] [
    [
        "sampleA",
        "1"
    ],
    {
        "id": "sampleA",
        "repeat": "1",
        "type": "normal",
        "fastq1": "sampleA_rep1_normal_R1.fastq.gz",
        "fastq2": "sampleA_rep1_normal_R2.fastq.gz"
    },
    {
        "id": "sampleA",
        "repeat": "1",
        "type": "tumor",
        "fastq1": "sampleA_rep1_tumor_R1.fastq.gz",
        "fastq2": "sampleA_rep1_tumor_R2.fastq.gz"
    }
]
[DUMP: joined] [
    [
        "sampleA",
        "2"
    ],
    {
        "id": "sampleA",
        "repeat": "2",
        "type": "normal",
        "fastq1": "sampleA_rep2_normal_R1.fastq.gz",
        "fastq2": "sampleA_rep2_normal_R2.fastq.gz"
    },
    {
        "id": "sampleA",
        "repeat": "2",
        "type": "tumor",
        "fastq1": "sampleA_rep2_tumor_R1.fastq.gz",
        "fastq2": "sampleA_rep2_tumor_R2.fastq.gz"
    }
]
[DUMP: joined] [
    [
        "sampleB",
        "1"
    ],
    {
        "id": "sampleB",
        "repeat": "1",
        "type": "normal",
        "fastq1": "sampleB_rep1_normal_R1.fastq.gz",
        "fastq2": "sampleB_rep1_normal_R2.fastq.gz"
    },
    {
        "id": "sampleB",
        "repeat": "1",
        "type": "tumor",
        "fastq1": "sampleB_rep1_tumor_R1.fastq.gz",
        "fastq2": "sampleB_rep1_tumor_R2.fastq.gz"
    }
]
[DUMP: joined] [
    [
        "sampleC",
        "1"
    ],
    {
        "id": "sampleC",
        "repeat": "1",
        "type": "normal",
        "fastq1": "sampleC_rep1_normal_R1.fastq.gz",
        "fastq2": "sampleC_rep1_normal_R2.fastq.gz"
    },
    {
        "id": "sampleC",
        "repeat": "1",
        "type": "tumor",
        "fastq1": "sampleC_rep1_tumor_R1.fastq.gz",
        "fastq2": "sampleC_rep1_tumor_R2.fastq.gz"
    }
]
```

Note how we have a tuple of two elements (`id` and `repeat` fields) as the first element of each joined result. This demonstrates how complex items can be used as a joining key, enabling fairly intricate matching between samples from the same conditions.

### 3.3. Use subMap to create a new joining key

We have an issue from the above example. We have lost the field names from the original joining key, i.e. the `id` and `repeat` fields are just a list of two values. If we want to retain the field names so we can access them later by name we can use the [`subMap` method](<https://docs.groovy-lang.org/latest/html/groovy-jdk/java/util/Map.html#subMap(java.util.Collection)>).

The `subMap` method takes a map and returns a new map with only the key-value pairs specified in the argument. In this case we want to specify the `id` and `repeat` fields.

_Before:_

```groovy title="main.nf" linenums="1"
workflow {
    ch_samplesheet = Channel.fromPath("./data/samplesheet.csv")
                        .splitCsv(header: true)
    ch_normal_samples = ch_samplesheet
                        .filter { sample -> sample.type == 'normal' }
                        .map { sample -> [
                                [sample.id, sample.repeat],
                                sample
                            ]
                        }
                        .dump(tag: 'normal')
    ch_tumor_samples = ch_samplesheet
                        .filter { sample -> sample.type == "tumor" }
                        .map { sample -> [
                                [sample.id, sample.repeat],
                                sample
                            ]
                        }
                        .dump(tag: 'tumor')
    ch_joined_samples = ch_normal_samples
                        .join(ch_tumor_samples)
                        .dump(tag: 'joined', pretty: true)
}
```

_After:_

```groovy title="main.nf" linenums="1"
workflow {
    ch_samplesheet = Channel.fromPath("./data/samplesheet.csv")
                        .splitCsv(header: true)
    ch_normal_samples = ch_samplesheet
                        .filter { sample -> sample.type == 'normal' }
                        .map { sample -> [
                                sample.subMap(['id', 'repeat']),
                                sample
                            ]
                        }
                        .dump(tag: 'normal')
    ch_tumor_samples = ch_samplesheet
                        .filter { sample -> sample.type == "tumor" }
                        .map { sample -> [
                                sample.subMap(['id', 'repeat']),
                                sample
                            ]
                        }
                        .dump(tag: 'tumor')
    ch_joined_samples = ch_normal_samples
                        .join(ch_tumor_samples)
                        .dump(tag: 'joined', pretty: true)
}
```

```bash title="View normal and tumor samples"
nextflow run main.nf -dump-channels joined
```

```console title="View normal and tumor samples"
 N E X T F L O W   ~  version 24.10.5

Launching `main.nf` [insane_gautier] DSL2 - revision: bf5b9a6d37

[DUMP: joined] [
    {
        "id": "sampleA",
        "repeat": "1"
    },
    {
        "id": "sampleA",
        "repeat": "1",
        "type": "normal",
        "fastq1": "sampleA_rep1_normal_R1.fastq.gz",
        "fastq2": "sampleA_rep1_normal_R2.fastq.gz"
    },
    {
        "id": "sampleA",
        "repeat": "1",
        "type": "tumor",
        "fastq1": "sampleA_rep1_tumor_R1.fastq.gz",
        "fastq2": "sampleA_rep1_tumor_R2.fastq.gz"
    }
]
[DUMP: joined] [
    {
        "id": "sampleA",
        "repeat": "2"
    },
    {
        "id": "sampleA",
        "repeat": "2",
        "type": "normal",
        "fastq1": "sampleA_rep2_normal_R1.fastq.gz",
        "fastq2": "sampleA_rep2_normal_R2.fastq.gz"
    },
    {
        "id": "sampleA",
        "repeat": "2",
        "type": "tumor",
        "fastq1": "sampleA_rep2_tumor_R1.fastq.gz",
        "fastq2": "sampleA_rep2_tumor_R2.fastq.gz"
    }
]
[DUMP: joined] [
    {
        "id": "sampleB",
        "repeat": "1"
    },
    {
        "id": "sampleB",
        "repeat": "1",
        "type": "normal",
        "fastq1": "sampleB_rep1_normal_R1.fastq.gz",
        "fastq2": "sampleB_rep1_normal_R2.fastq.gz"
    },
    {
        "id": "sampleB",
        "repeat": "1",
        "type": "tumor",
        "fastq1": "sampleB_rep1_tumor_R1.fastq.gz",
        "fastq2": "sampleB_rep1_tumor_R2.fastq.gz"
    }
]
[DUMP: joined] [
    {
        "id": "sampleC",
        "repeat": "1"
    },
    {
        "id": "sampleC",
        "repeat": "1",
        "type": "normal",
        "fastq1": "sampleC_rep1_normal_R1.fastq.gz",
        "fastq2": "sampleC_rep1_normal_R2.fastq.gz"
    },
    {
        "id": "sampleC",
        "repeat": "1",
        "type": "tumor",
        "fastq1": "sampleC_rep1_tumor_R1.fastq.gz",
        "fastq2": "sampleC_rep1_tumor_R2.fastq.gz"
    }
]
```

Now we have a new joining key that not only includes the `id` and `repeat` fields but also retains the field names so we can access them later by name, e.g. `sample.id` and `sample.repeat`.

### 3.4. Use a named closure in map

Since we are re-using the same map in multiple places, we run the risk of introducing errors if we accidentally change the map in one place but not the other. To avoid this, we can use a named closure in the map. A named closure allows us to make a reusable function we can call later within a map.

To do so, first we define the closure as a new variable:

_Before:_

```groovy title="main.nf" linenums="1"
workflow {
    samplesheet = Channel.fromPath("./data/samplesheet.csv")
                        .splitCsv(header: true)
```

_After:_

```groovy title="main.nf" linenums="1"
workflow {
    getSampleIdAndReplicate = { sample -> [ sample.subMap(['id', 'repeat']), sample ] }
    samplesheet = Channel.fromPath("./data/samplesheet.csv")
                        .splitCsv(header: true)
```

We have taken the map we used previously and defined it as a named variable we can call later. Let's implement it in our workflow:

_Before:_

```groovy title="main.nf" linenums="5"
    ch_normal_samples = ch_samplesheet
                        .filter { sample -> sample.type == 'normal' }
                        .map { sample -> [
                            sample.subMap(['id', 'repeat']),
                            sample
                          ]
                        }
                        .dump(tag: 'normal')
    ch_tumor_samples = ch_samplesheet
                        .filter { sample -> sample.type == "tumor" }
                        .map { sample -> [
                            sample.subMap(['id', 'repeat']),
                            sample
                          ]
                        }
                        .dump(tag: 'tumor')
```

_After:_

```groovy title="main.nf" linenums="5"
    ch_normal_samples = ch_samplesheet
                        .filter { sample -> sample.type == 'normal' }
                        .map ( getSampleIdAndReplicate )
                        .dump(tag: 'normal')
    ch_tumor_samples = ch_samplesheet
                        .filter { sample -> sample.type == "tumor" }
                        .map ( getSampleIdAndReplicate )
                        .dump(tag: 'tumor')
```

!!! note
The `map` operator has switched from using `{ }` to using `( )` to pass the closure as an argument. This is because the `map` operator expects a closure as an argument and `{ }` is used to define an anonymous closure. When calling a named closure, use the `( )` syntax.

```bash title="View normal and tumor samples"
nextflow run main.nf -dump-channels joined
```

```console title="View normal and tumor samples"
 N E X T F L O W   ~  version 24.10.5

Launching `main.nf` [trusting_boltzmann] DSL2 - revision: 0b1cd77e3b

[DUMP: joined] [
    {
        "id": "sampleA",
        "repeat": "1"
    },
    {
        "id": "sampleA",
        "repeat": "1",
        "type": "normal",
        "fastq1": "sampleA_rep1_normal_R1.fastq.gz",
        "fastq2": "sampleA_rep1_normal_R2.fastq.gz"
    },
    {
        "id": "sampleA",
        "repeat": "1",
        "type": "tumor",
        "fastq1": "sampleA_rep1_tumor_R1.fastq.gz",
        "fastq2": "sampleA_rep1_tumor_R2.fastq.gz"
    }
]
[DUMP: joined] [
    {
        "id": "sampleA",
        "repeat": "2"
    },
    {
        "id": "sampleA",
        "repeat": "2",
        "type": "normal",
        "fastq1": "sampleA_rep2_normal_R1.fastq.gz",
        "fastq2": "sampleA_rep2_normal_R2.fastq.gz"
    },
    {
        "id": "sampleA",
        "repeat": "2",
        "type": "tumor",
        "fastq1": "sampleA_rep2_tumor_R1.fastq.gz",
        "fastq2": "sampleA_rep2_tumor_R2.fastq.gz"
    }
]
[DUMP: joined] [
    {
        "id": "sampleB",
        "repeat": "1"
    },
    {
        "id": "sampleB",
        "repeat": "1",
        "type": "normal",
        "fastq1": "sampleB_rep1_normal_R1.fastq.gz",
        "fastq2": "sampleB_rep1_normal_R2.fastq.gz"
    },
    {
        "id": "sampleB",
        "repeat": "1",
        "type": "tumor",
        "fastq1": "sampleB_rep1_tumor_R1.fastq.gz",
        "fastq2": "sampleB_rep1_tumor_R2.fastq.gz"
    }
]
[DUMP: joined] [
    {
        "id": "sampleC",
        "repeat": "1"
    },
    {
        "id": "sampleC",
        "repeat": "1",
        "type": "normal",
        "fastq1": "sampleC_rep1_normal_R1.fastq.gz",
        "fastq2": "sampleC_rep1_normal_R2.fastq.gz"
    },
    {
        "id": "sampleC",
        "repeat": "1",
        "type": "tumor",
        "fastq1": "sampleC_rep1_tumor_R1.fastq.gz",
        "fastq2": "sampleC_rep1_tumor_R2.fastq.gz"
    }
]
```

Using a named closure in the map allows us to reuse the same map in multiple places which reduces our risk of introducing errors. It also makes the code more readable and easier to maintain.

### Takeaway

In this section, you've learned:

- **Manipulating Tuples**: How to use `map` to isolate a field in a tuple
- **Joining Tuples**: How to use `join` to combine tuples based on the first field
- **Creating Joining Keys**: How to use `subMap` to create a new joining key
- **Named Closures**: How to use a named closure in map

You now have a workflow that can split a samplesheet, filter the normal and tumor samples, join them together by sample ID and replicate number, then dump the results.

This is a common pattern in bioinformatics workflows where you need to match up samples after processing independently, so it is a useful skill. Next, we will look at repeating a sample multiple times.

## 4. Spread samples over intervals

Spreading samples over different conditions is a common pattern in bioinformatics workflows. For example, it is used to spread variant calling over a range of intervals. This can help distribute work across multiple cores or nodes and make the pipelines more efficient and be turned around faster.

In the next section, we will demonstrate how to take our existing samples and repeat each one for every interval. In this way, we will have a single sample for each input interval. We will also multiply our number of samples by the number of intervals, so get ready for a busy terminal!

### 4.1. Spread samples over intervals using `combine`

Let's start by creating a channel of intervals. To keep life simple, we will just use 3 intervals we will manually define. In a real workflow, you could read these in from a file input or even create a channel with lots of interval files.

_Before:_

```groovy title="main.nf" linenums="15"
                        .dump(tag: 'joined', pretty: true)
}
```

_After:_

```groovy title="main.nf" linenums="24"
                        .dump(tag: 'joined', pretty: true)
    ch_intervals = Channel.of('chr1', 'chr2', 'chr3')
                    .dump(tag: "intervals")
}
```

Now remember, we want to repeat each sample for each interval. This is sometimes referred to as the Cartesian product of the samples and intervals. We can achieve this by using the [`combine` operator](https://www.nextflow.io/docs/latest/operator.html#combine). This will take every item from channel 1 and repeat it for each item in channel 2. Let's add a combine operator to our workflow:

_Before:_

```groovy title="main.nf" linenums="26"
    ch_intervals = Channel.of('chr1', 'chr2', 'chr3')
                    .dump(tag: "intervals")
}
```

_After:_

```groovy title="main.nf" linenums="26"
    ch_intervals = Channel.of('chr1', 'chr2', 'chr3')
                    .dump(tag: "intervals")

    ch_combined_samples = ch_joined_samples.combine(ch_intervals)
                        .dump(tag: 'combined')
}
```

Now let's run it and see what happens:

```bash title="View combined samples"
nextflow run main.nf -dump-channels combined
```

```console title="View combined samples"
 N E X T F L O W   ~  version 24.10.5

Launching `main.nf` [extravagant_maxwell] DSL2 - revision: 459bde3584

[DUMP: combined] [['id':'sampleA', 'repeat':'1'], ['id':'sampleA', 'repeat':'1', 'type':'normal', 'fastq1':'sampleA_rep1_normal_R1.fastq.gz', 'fastq2':'sampleA_rep1_normal_R2.fastq.gz'], ['id':'sampleA', 'repeat':'1', 'type':'tumor', 'fastq1':'sampleA_rep1_tumor_R1.fastq.gz', 'fastq2':'sampleA_rep1_tumor_R2.fastq.gz'], 'chr1']
[DUMP: combined] [['id':'sampleA', 'repeat':'1'], ['id':'sampleA', 'repeat':'1', 'type':'normal', 'fastq1':'sampleA_rep1_normal_R1.fastq.gz', 'fastq2':'sampleA_rep1_normal_R2.fastq.gz'], ['id':'sampleA', 'repeat':'1', 'type':'tumor', 'fastq1':'sampleA_rep1_tumor_R1.fastq.gz', 'fastq2':'sampleA_rep1_tumor_R2.fastq.gz'], 'chr2']
[DUMP: combined] [['id':'sampleA', 'repeat':'1'], ['id':'sampleA', 'repeat':'1', 'type':'normal', 'fastq1':'sampleA_rep1_normal_R1.fastq.gz', 'fastq2':'sampleA_rep1_normal_R2.fastq.gz'], ['id':'sampleA', 'repeat':'1', 'type':'tumor', 'fastq1':'sampleA_rep1_tumor_R1.fastq.gz', 'fastq2':'sampleA_rep1_tumor_R2.fastq.gz'], 'chr3']
[DUMP: combined] [['id':'sampleA', 'repeat':'2'], ['id':'sampleA', 'repeat':'2', 'type':'normal', 'fastq1':'sampleA_rep2_normal_R1.fastq.gz', 'fastq2':'sampleA_rep2_normal_R2.fastq.gz'], ['id':'sampleA', 'repeat':'2', 'type':'tumor', 'fastq1':'sampleA_rep2_tumor_R1.fastq.gz', 'fastq2':'sampleA_rep2_tumor_R2.fastq.gz'], 'chr1']
[DUMP: combined] [['id':'sampleA', 'repeat':'2'], ['id':'sampleA', 'repeat':'2', 'type':'normal', 'fastq1':'sampleA_rep2_normal_R1.fastq.gz', 'fastq2':'sampleA_rep2_normal_R2.fastq.gz'], ['id':'sampleA', 'repeat':'2', 'type':'tumor', 'fastq1':'sampleA_rep2_tumor_R1.fastq.gz', 'fastq2':'sampleA_rep2_tumor_R2.fastq.gz'], 'chr2']
[DUMP: combined] [['id':'sampleA', 'repeat':'2'], ['id':'sampleA', 'repeat':'2', 'type':'normal', 'fastq1':'sampleA_rep2_normal_R1.fastq.gz', 'fastq2':'sampleA_rep2_normal_R2.fastq.gz'], ['id':'sampleA', 'repeat':'2', 'type':'tumor', 'fastq1':'sampleA_rep2_tumor_R1.fastq.gz', 'fastq2':'sampleA_rep2_tumor_R2.fastq.gz'], 'chr3']
[DUMP: combined] [['id':'sampleB', 'repeat':'1'], ['id':'sampleB', 'repeat':'1', 'type':'normal', 'fastq1':'sampleB_rep1_normal_R1.fastq.gz', 'fastq2':'sampleB_rep1_normal_R2.fastq.gz'], ['id':'sampleB', 'repeat':'1', 'type':'tumor', 'fastq1':'sampleB_rep1_tumor_R1.fastq.gz', 'fastq2':'sampleB_rep1_tumor_R2.fastq.gz'], 'chr1']
[DUMP: combined] [['id':'sampleB', 'repeat':'1'], ['id':'sampleB', 'repeat':'1', 'type':'normal', 'fastq1':'sampleB_rep1_normal_R1.fastq.gz', 'fastq2':'sampleB_rep1_normal_R2.fastq.gz'], ['id':'sampleB', 'repeat':'1', 'type':'tumor', 'fastq1':'sampleB_rep1_tumor_R1.fastq.gz', 'fastq2':'sampleB_rep1_tumor_R2.fastq.gz'], 'chr2']
[DUMP: combined] [['id':'sampleB', 'repeat':'1'], ['id':'sampleB', 'repeat':'1', 'type':'normal', 'fastq1':'sampleB_rep1_normal_R1.fastq.gz', 'fastq2':'sampleB_rep1_normal_R2.fastq.gz'], ['id':'sampleB', 'repeat':'1', 'type':'tumor', 'fastq1':'sampleB_rep1_tumor_R1.fastq.gz', 'fastq2':'sampleB_rep1_tumor_R2.fastq.gz'], 'chr3']
[DUMP: combined] [['id':'sampleC', 'repeat':'1'], ['id':'sampleC', 'repeat':'1', 'type':'normal', 'fastq1':'sampleC_rep1_normal_R1.fastq.gz', 'fastq2':'sampleC_rep1_normal_R2.fastq.gz'], ['id':'sampleC', 'repeat':'1', 'type':'tumor', 'fastq1':'sampleC_rep1_tumor_R1.fastq.gz', 'fastq2':'sampleC_rep1_tumor_R2.fastq.gz'], 'chr1']
[DUMP: combined] [['id':'sampleC', 'repeat':'1'], ['id':'sampleC', 'repeat':'1', 'type':'normal', 'fastq1':'sampleC_rep1_normal_R1.fastq.gz', 'fastq2':'sampleC_rep1_normal_R2.fastq.gz'], ['id':'sampleC', 'repeat':'1', 'type':'tumor', 'fastq1':'sampleC_rep1_tumor_R1.fastq.gz', 'fastq2':'sampleC_rep1_tumor_R2.fastq.gz'], 'chr2']
[DUMP: combined] [['id':'sampleC', 'repeat':'1'], ['id':'sampleC', 'repeat':'1', 'type':'normal', 'fastq1':'sampleC_rep1_normal_R1.fastq.gz', 'fastq2':'sampleC_rep1_normal_R2.fastq.gz'], ['id':'sampleC', 'repeat':'1', 'type':'tumor', 'fastq1':'sampleC_rep1_tumor_R1.fastq.gz', 'fastq2':'sampleC_rep1_tumor_R2.fastq.gz'], 'chr3']
```

Success! We have repeated every sample for every single interval in our 3 interval list. We've effectively tripled the number of items in our channel. It's a little hard to read though, so in the next section we will tidy it up.

### 4.2. Organise the channel

We can use the `map` operator to tidy and refactor our sample data so it's easier to understand. Let's move the intervals string to the joining map at the first element.

_Before:_

```groovy title="main.nf" linenums="19"
    ch_combined_samples = joined_samples.combine(ch_intervals)
                        .dump(tag: 'combined')
}
```

_After:_

```groovy title="main.nf" linenums="19"
    ch_combined_samples = ch_joined_samples.combine(ch_intervals)
                        .map { grouping_key, normal, tumor, interval ->
                            [
                                grouping_key + [interval: interval],
                                normal,
                                tumor
                            ]

                        }
                        .dump(tag: 'combined')
}
```

Wait? What did we do here? Let's go over it piece by piece.

First, we use a map operator to iterate over every item in the channel. By using the names `grouping_key`, `normal`, `tumor` and `interval`, we can refer to the elements in the tuple by name instead of by index. This makes the code more readable and easier to understand.

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
nextflow run main.nf -dump-channels combined
```

```console title="View combined samples"
 N E X T F L O W   ~  version 24.10.5

Launching `main.nf` [focused_curie] DSL2 - revision: 9953685fec

[DUMP: combined] [['id':'sampleA', 'repeat':'1', 'interval':'chr1'], ['id':'sampleA', 'repeat':'1', 'type':'normal', 'fastq1':'sampleA_rep1_normal_R1.fastq.gz', 'fastq2':'sampleA_rep1_normal_R2.fastq.gz'], ['id':'sampleA', 'repeat':'1', 'type':'tumor', 'fastq1':'sampleA_rep1_tumor_R1.fastq.gz', 'fastq2':'sampleA_rep1_tumor_R2.fastq.gz']]
[DUMP: combined] [['id':'sampleA', 'repeat':'1', 'interval':'chr2'], ['id':'sampleA', 'repeat':'1', 'type':'normal', 'fastq1':'sampleA_rep1_normal_R1.fastq.gz', 'fastq2':'sampleA_rep1_normal_R2.fastq.gz'], ['id':'sampleA', 'repeat':'1', 'type':'tumor', 'fastq1':'sampleA_rep1_tumor_R1.fastq.gz', 'fastq2':'sampleA_rep1_tumor_R2.fastq.gz']]
[DUMP: combined] [['id':'sampleA', 'repeat':'1', 'interval':'chr3'], ['id':'sampleA', 'repeat':'1', 'type':'normal', 'fastq1':'sampleA_rep1_normal_R1.fastq.gz', 'fastq2':'sampleA_rep1_normal_R2.fastq.gz'], ['id':'sampleA', 'repeat':'1', 'type':'tumor', 'fastq1':'sampleA_rep1_tumor_R1.fastq.gz', 'fastq2':'sampleA_rep1_tumor_R2.fastq.gz']]
[DUMP: combined] [['id':'sampleA', 'repeat':'2', 'interval':'chr1'], ['id':'sampleA', 'repeat':'2', 'type':'normal', 'fastq1':'sampleA_rep2_normal_R1.fastq.gz', 'fastq2':'sampleA_rep2_normal_R2.fastq.gz'], ['id':'sampleA', 'repeat':'2', 'type':'tumor', 'fastq1':'sampleA_rep2_tumor_R1.fastq.gz', 'fastq2':'sampleA_rep2_tumor_R2.fastq.gz']]
[DUMP: combined] [['id':'sampleA', 'repeat':'2', 'interval':'chr2'], ['id':'sampleA', 'repeat':'2', 'type':'normal', 'fastq1':'sampleA_rep2_normal_R1.fastq.gz', 'fastq2':'sampleA_rep2_normal_R2.fastq.gz'], ['id':'sampleA', 'repeat':'2', 'type':'tumor', 'fastq1':'sampleA_rep2_tumor_R1.fastq.gz', 'fastq2':'sampleA_rep2_tumor_R2.fastq.gz']]
[DUMP: combined] [['id':'sampleA', 'repeat':'2', 'interval':'chr3'], ['id':'sampleA', 'repeat':'2', 'type':'normal', 'fastq1':'sampleA_rep2_normal_R1.fastq.gz', 'fastq2':'sampleA_rep2_normal_R2.fastq.gz'], ['id':'sampleA', 'repeat':'2', 'type':'tumor', 'fastq1':'sampleA_rep2_tumor_R1.fastq.gz', 'fastq2':'sampleA_rep2_tumor_R2.fastq.gz']]
[DUMP: combined] [['id':'sampleB', 'repeat':'1', 'interval':'chr1'], ['id':'sampleB', 'repeat':'1', 'type':'normal', 'fastq1':'sampleB_rep1_normal_R1.fastq.gz', 'fastq2':'sampleB_rep1_normal_R2.fastq.gz'], ['id':'sampleB', 'repeat':'1', 'type':'tumor', 'fastq1':'sampleB_rep1_tumor_R1.fastq.gz', 'fastq2':'sampleB_rep1_tumor_R2.fastq.gz']]
[DUMP: combined] [['id':'sampleB', 'repeat':'1', 'interval':'chr2'], ['id':'sampleB', 'repeat':'1', 'type':'normal', 'fastq1':'sampleB_rep1_normal_R1.fastq.gz', 'fastq2':'sampleB_rep1_normal_R2.fastq.gz'], ['id':'sampleB', 'repeat':'1', 'type':'tumor', 'fastq1':'sampleB_rep1_tumor_R1.fastq.gz', 'fastq2':'sampleB_rep1_tumor_R2.fastq.gz']]
[DUMP: combined] [['id':'sampleB', 'repeat':'1', 'interval':'chr3'], ['id':'sampleB', 'repeat':'1', 'type':'normal', 'fastq1':'sampleB_rep1_normal_R1.fastq.gz', 'fastq2':'sampleB_rep1_normal_R2.fastq.gz'], ['id':'sampleB', 'repeat':'1', 'type':'tumor', 'fastq1':'sampleB_rep1_tumor_R1.fastq.gz', 'fastq2':'sampleB_rep1_tumor_R2.fastq.gz']]
[DUMP: combined] [['id':'sampleC', 'repeat':'1', 'interval':'chr1'], ['id':'sampleC', 'repeat':'1', 'type':'normal', 'fastq1':'sampleC_rep1_normal_R1.fastq.gz', 'fastq2':'sampleC_rep1_normal_R2.fastq.gz'], ['id':'sampleC', 'repeat':'1', 'type':'tumor', 'fastq1':'sampleC_rep1_tumor_R1.fastq.gz', 'fastq2':'sampleC_rep1_tumor_R2.fastq.gz']]
[DUMP: combined] [['id':'sampleC', 'repeat':'1', 'interval':'chr2'], ['id':'sampleC', 'repeat':'1', 'type':'normal', 'fastq1':'sampleC_rep1_normal_R1.fastq.gz', 'fastq2':'sampleC_rep1_normal_R2.fastq.gz'], ['id':'sampleC', 'repeat':'1', 'type':'tumor', 'fastq1':'sampleC_rep1_tumor_R1.fastq.gz', 'fastq2':'sampleC_rep1_tumor_R2.fastq.gz']]
[DUMP: combined] [['id':'sampleC', 'repeat':'1', 'interval':'chr3'], ['id':'sampleC', 'repeat':'1', 'type':'normal', 'fastq1':'sampleC_rep1_normal_R1.fastq.gz', 'fastq2':'sampleC_rep1_normal_R2.fastq.gz'], ['id':'sampleC', 'repeat':'1', 'type':'tumor', 'fastq1':'sampleC_rep1_tumor_R1.fastq.gz', 'fastq2':'sampleC_rep1_tumor_R2.fastq.gz']]
```

Using `map` to coerce your data into the correct structure can be tricky, but it's crucial to correctly splitting and grouping effectively.

### Takeaway

In this section, you've learned:

- **Spreading samples over intervals**: How to use `combine` to repeat samples over intervals

### 5. Aggregating samples

In the previous section, we learned how to split a samplesheet and filter the normal and tumor samples. But this only covers a single type of joining. What if we want to group samples by a specific attribute? For example, instead of joining matched normal-tumor pairs, we might want to process all samples from "sampleA" together regardless of their type. This pattern is common in bioinformatics workflows where you may want to process related samples separately for efficiency reasons before comparing or combining the results at the end.

Nextflow includes built in methods to do this, the main one we will look at is `groupTuple`.

### 5.1. Grouping samples using `groupTuple`

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

_Before:_

```groovy title="main.nf" linenums="19"
    ch_combined_samples = ch_joined_samples.combine(ch_intervals)
                        .map { grouping_key, normal, tumor, interval ->
                            [
                                grouping_key + [interval: interval],
                                normal,
                                tumor
                            ]

                        }
                        .dump(tag: 'combined')
}
```

_After:_

```groovy title="main.nf" linenums="19"
    ch_combined_samples = ch_joined_samples.combine(ch_intervals)
                        .map { grouping_key, normal, tumor, interval ->
                            [
                                grouping_key + [interval: interval],
                                normal,
                                tumor
                            ]

                        }
                        .dump(tag: 'combined')

    ch_grouped_samples = ch_combined_samples.map { grouping_key, normal, tumor ->
                            [
                                grouping_key.subMap('id', 'interval'),
                                normal,
                                tumor
                            ]

                        }
                      .dump(tag: 'grouped')
}
```

Let's run it again and check the channel contents:

```bash title="View grouped samples"
nextflow run main.nf -dump-channels grouped
```

```console title="View grouped samples"
 N E X T F L O W   ~  version 24.10.5

Launching `main.nf` [fabulous_baekeland] DSL2 - revision: 5d2d687351

[DUMP: grouped] [['id':'sampleA', 'interval':'chr1'], ['id':'sampleA', 'repeat':'1', 'type':'normal', 'fastq1':'sampleA_rep1_normal_R1.fastq.gz', 'fastq2':'sampleA_rep1_normal_R2.fastq.gz'], ['id':'sampleA', 'repeat':'1', 'type':'tumor', 'fastq1':'sampleA_rep1_tumor_R1.fastq.gz', 'fastq2':'sampleA_rep1_tumor_R2.fastq.gz']]
[DUMP: grouped] [['id':'sampleA', 'interval':'chr2'], ['id':'sampleA', 'repeat':'1', 'type':'normal', 'fastq1':'sampleA_rep1_normal_R1.fastq.gz', 'fastq2':'sampleA_rep1_normal_R2.fastq.gz'], ['id':'sampleA', 'repeat':'1', 'type':'tumor', 'fastq1':'sampleA_rep1_tumor_R1.fastq.gz', 'fastq2':'sampleA_rep1_tumor_R2.fastq.gz']]
[DUMP: grouped] [['id':'sampleA', 'interval':'chr3'], ['id':'sampleA', 'repeat':'1', 'type':'normal', 'fastq1':'sampleA_rep1_normal_R1.fastq.gz', 'fastq2':'sampleA_rep1_normal_R2.fastq.gz'], ['id':'sampleA', 'repeat':'1', 'type':'tumor', 'fastq1':'sampleA_rep1_tumor_R1.fastq.gz', 'fastq2':'sampleA_rep1_tumor_R2.fastq.gz']]
[DUMP: grouped] [['id':'sampleA', 'interval':'chr1'], ['id':'sampleA', 'repeat':'2', 'type':'normal', 'fastq1':'sampleA_rep2_normal_R1.fastq.gz', 'fastq2':'sampleA_rep2_normal_R2.fastq.gz'], ['id':'sampleA', 'repeat':'2', 'type':'tumor', 'fastq1':'sampleA_rep2_tumor_R1.fastq.gz', 'fastq2':'sampleA_rep2_tumor_R2.fastq.gz']]
[DUMP: grouped] [['id':'sampleA', 'interval':'chr2'], ['id':'sampleA', 'repeat':'2', 'type':'normal', 'fastq1':'sampleA_rep2_normal_R1.fastq.gz', 'fastq2':'sampleA_rep2_normal_R2.fastq.gz'], ['id':'sampleA', 'repeat':'2', 'type':'tumor', 'fastq1':'sampleA_rep2_tumor_R1.fastq.gz', 'fastq2':'sampleA_rep2_tumor_R2.fastq.gz']]
[DUMP: grouped] [['id':'sampleA', 'interval':'chr3'], ['id':'sampleA', 'repeat':'2', 'type':'normal', 'fastq1':'sampleA_rep2_normal_R1.fastq.gz', 'fastq2':'sampleA_rep2_normal_R2.fastq.gz'], ['id':'sampleA', 'repeat':'2', 'type':'tumor', 'fastq1':'sampleA_rep2_tumor_R1.fastq.gz', 'fastq2':'sampleA_rep2_tumor_R2.fastq.gz']]
[DUMP: grouped] [['id':'sampleB', 'interval':'chr1'], ['id':'sampleB', 'repeat':'1', 'type':'normal', 'fastq1':'sampleB_rep1_normal_R1.fastq.gz', 'fastq2':'sampleB_rep1_normal_R2.fastq.gz'], ['id':'sampleB', 'repeat':'1', 'type':'tumor', 'fastq1':'sampleB_rep1_tumor_R1.fastq.gz', 'fastq2':'sampleB_rep1_tumor_R2.fastq.gz']]
[DUMP: grouped] [['id':'sampleB', 'interval':'chr2'], ['id':'sampleB', 'repeat':'1', 'type':'normal', 'fastq1':'sampleB_rep1_normal_R1.fastq.gz', 'fastq2':'sampleB_rep1_normal_R2.fastq.gz'], ['id':'sampleB', 'repeat':'1', 'type':'tumor', 'fastq1':'sampleB_rep1_tumor_R1.fastq.gz', 'fastq2':'sampleB_rep1_tumor_R2.fastq.gz']]
[DUMP: grouped] [['id':'sampleB', 'interval':'chr3'], ['id':'sampleB', 'repeat':'1', 'type':'normal', 'fastq1':'sampleB_rep1_normal_R1.fastq.gz', 'fastq2':'sampleB_rep1_normal_R2.fastq.gz'], ['id':'sampleB', 'repeat':'1', 'type':'tumor', 'fastq1':'sampleB_rep1_tumor_R1.fastq.gz', 'fastq2':'sampleB_rep1_tumor_R2.fastq.gz']]
[DUMP: grouped] [['id':'sampleC', 'interval':'chr1'], ['id':'sampleC', 'repeat':'1', 'type':'normal', 'fastq1':'sampleC_rep1_normal_R1.fastq.gz', 'fastq2':'sampleC_rep1_normal_R2.fastq.gz'], ['id':'sampleC', 'repeat':'1', 'type':'tumor', 'fastq1':'sampleC_rep1_tumor_R1.fastq.gz', 'fastq2':'sampleC_rep1_tumor_R2.fastq.gz']]
[DUMP: grouped] [['id':'sampleC', 'interval':'chr2'], ['id':'sampleC', 'repeat':'1', 'type':'normal', 'fastq1':'sampleC_rep1_normal_R1.fastq.gz', 'fastq2':'sampleC_rep1_normal_R2.fastq.gz'], ['id':'sampleC', 'repeat':'1', 'type':'tumor', 'fastq1':'sampleC_rep1_tumor_R1.fastq.gz', 'fastq2':'sampleC_rep1_tumor_R2.fastq.gz']]
[DUMP: grouped] [['id':'sampleC', 'interval':'chr3'], ['id':'sampleC', 'repeat':'1', 'type':'normal', 'fastq1':'sampleC_rep1_normal_R1.fastq.gz', 'fastq2':'sampleC_rep1_normal_R2.fastq.gz'], ['id':'sampleC', 'repeat':'1', 'type':'tumor', 'fastq1':'sampleC_rep1_tumor_R1.fastq.gz', 'fastq2':'sampleC_rep1_tumor_R2.fastq.gz']]
```

We can see that we have successfully isolated the `id` and `interval` fields, but not grouped the samples yet.

Let's now group the samples by this new grouping element, using the [`groupTuple` operator](https://www.nextflow.io/docs/latest/operator.html#grouptuple).

_Before:_

```groovy title="main.nf" linenums="30"
    ch_grouped_samples = ch_combined_samples.map { grouping_key, normal, tumor ->
                            [
                                grouping_key.subMap('id', 'interval'),
                                grouping_key,
                                normal,
                                tumor
                            ]

                        }
                        .dump(tag: 'grouped')
}
```

_After:_

```groovy title="main.nf" linenums="29"
    ch_grouped_samples = ch_combined_samples.map { grouping_key, normal, tumor ->
                            [
                                grouping_key.subMap('id', 'interval'),
                                normal,
                                tumor
                            ]

                        }
                        .groupTuple()
                        .dump(tag: 'grouped')
}
```

Simple, huh? We just added a single line of code. Let's see what happens when we run it:

```bash title="View grouped samples"
nextflow run main.nf -dump-channels grouped
```

```console title="View grouped samples"
 N E X T F L O W   ~  version 24.10.5

Launching `main.nf` [reverent_nightingale] DSL2 - revision: 72c6664d6f

[DUMP: grouped] [['id':'sampleA', 'interval':'chr1'], [['id':'sampleA', 'repeat':'1', 'type':'normal', 'fastq1':'sampleA_rep1_normal_R1.fastq.gz', 'fastq2':'sampleA_rep1_normal_R2.fastq.gz'], ['id':'sampleA', 'repeat':'2', 'type':'normal', 'fastq1':'sampleA_rep2_normal_R1.fastq.gz', 'fastq2':'sampleA_rep2_normal_R2.fastq.gz']], [['id':'sampleA', 'repeat':'1', 'type':'tumor', 'fastq1':'sampleA_rep1_tumor_R1.fastq.gz', 'fastq2':'sampleA_rep1_tumor_R2.fastq.gz'], ['id':'sampleA', 'repeat':'2', 'type':'tumor', 'fastq1':'sampleA_rep2_tumor_R1.fastq.gz', 'fastq2':'sampleA_rep2_tumor_R2.fastq.gz']]]
[DUMP: grouped] [['id':'sampleA', 'interval':'chr2'], [['id':'sampleA', 'repeat':'1', 'type':'normal', 'fastq1':'sampleA_rep1_normal_R1.fastq.gz', 'fastq2':'sampleA_rep1_normal_R2.fastq.gz'], ['id':'sampleA', 'repeat':'2', 'type':'normal', 'fastq1':'sampleA_rep2_normal_R1.fastq.gz', 'fastq2':'sampleA_rep2_normal_R2.fastq.gz']], [['id':'sampleA', 'repeat':'1', 'type':'tumor', 'fastq1':'sampleA_rep1_tumor_R1.fastq.gz', 'fastq2':'sampleA_rep1_tumor_R2.fastq.gz'], ['id':'sampleA', 'repeat':'2', 'type':'tumor', 'fastq1':'sampleA_rep2_tumor_R1.fastq.gz', 'fastq2':'sampleA_rep2_tumor_R2.fastq.gz']]]
[DUMP: grouped] [['id':'sampleA', 'interval':'chr3'], [['id':'sampleA', 'repeat':'1', 'type':'normal', 'fastq1':'sampleA_rep1_normal_R1.fastq.gz', 'fastq2':'sampleA_rep1_normal_R2.fastq.gz'], ['id':'sampleA', 'repeat':'2', 'type':'normal', 'fastq1':'sampleA_rep2_normal_R1.fastq.gz', 'fastq2':'sampleA_rep2_normal_R2.fastq.gz']], [['id':'sampleA', 'repeat':'1', 'type':'tumor', 'fastq1':'sampleA_rep1_tumor_R1.fastq.gz', 'fastq2':'sampleA_rep1_tumor_R2.fastq.gz'], ['id':'sampleA', 'repeat':'2', 'type':'tumor', 'fastq1':'sampleA_rep2_tumor_R1.fastq.gz', 'fastq2':'sampleA_rep2_tumor_R2.fastq.gz']]]
[DUMP: grouped] [['id':'sampleB', 'interval':'chr1'], [['id':'sampleB', 'repeat':'1', 'type':'normal', 'fastq1':'sampleB_rep1_normal_R1.fastq.gz', 'fastq2':'sampleB_rep1_normal_R2.fastq.gz']], [['id':'sampleB', 'repeat':'1', 'type':'tumor', 'fastq1':'sampleB_rep1_tumor_R1.fastq.gz', 'fastq2':'sampleB_rep1_tumor_R2.fastq.gz']]]
[DUMP: grouped] [['id':'sampleB', 'interval':'chr2'], [['id':'sampleB', 'repeat':'1', 'type':'normal', 'fastq1':'sampleB_rep1_normal_R1.fastq.gz', 'fastq2':'sampleB_rep1_normal_R2.fastq.gz']], [['id':'sampleB', 'repeat':'1', 'type':'tumor', 'fastq1':'sampleB_rep1_tumor_R1.fastq.gz', 'fastq2':'sampleB_rep1_tumor_R2.fastq.gz']]]
[DUMP: grouped] [['id':'sampleB', 'interval':'chr3'], [['id':'sampleB', 'repeat':'1', 'type':'normal', 'fastq1':'sampleB_rep1_normal_R1.fastq.gz', 'fastq2':'sampleB_rep1_normal_R2.fastq.gz']], [['id':'sampleB', 'repeat':'1', 'type':'tumor', 'fastq1':'sampleB_rep1_tumor_R1.fastq.gz', 'fastq2':'sampleB_rep1_tumor_R2.fastq.gz']]]
[DUMP: grouped] [['id':'sampleC', 'interval':'chr1'], [['id':'sampleC', 'repeat':'1', 'type':'normal', 'fastq1':'sampleC_rep1_normal_R1.fastq.gz', 'fastq2':'sampleC_rep1_normal_R2.fastq.gz']], [['id':'sampleC', 'repeat':'1', 'type':'tumor', 'fastq1':'sampleC_rep1_tumor_R1.fastq.gz', 'fastq2':'sampleC_rep1_tumor_R2.fastq.gz']]]
[DUMP: grouped] [['id':'sampleC', 'interval':'chr2'], [['id':'sampleC', 'repeat':'1', 'type':'normal', 'fastq1':'sampleC_rep1_normal_R1.fastq.gz', 'fastq2':'sampleC_rep1_normal_R2.fastq.gz']], [['id':'sampleC', 'repeat':'1', 'type':'tumor', 'fastq1':'sampleC_rep1_tumor_R1.fastq.gz', 'fastq2':'sampleC_rep1_tumor_R2.fastq.gz']]]
[DUMP: grouped] [['id':'sampleC', 'interval':'chr3'], [['id':'sampleC', 'repeat':'1', 'type':'normal', 'fastq1':'sampleC_rep1_normal_R1.fastq.gz', 'fastq2':'sampleC_rep1_normal_R2.fastq.gz']], [['id':'sampleC', 'repeat':'1', 'type':'tumor', 'fastq1':'sampleC_rep1_tumor_R1.fastq.gz', 'fastq2':'sampleC_rep1_tumor_R2.fastq.gz']]]
```

It's a little awkward to read! If you're having trouble visualizing it, you can use the `pretty` flag of `dump` to make it easier to read:

_Before:_

```groovy title="main.nf" linenums="40"
                    .dump(tag: 'grouped')
}
```

_After:_

```groovy title="main.nf" linenums="40"
                    .dump(tag: 'grouped', pretty: true)
}
```

Note, we only include the first sample to keep this concise!

```console title="View grouped samples"
 N E X T F L O W   ~  version 24.10.5

Launching `main.nf` [dreamy_lichterman] DSL2 - revision: 953a5dd264

[DUMP: grouped] [
    {
        "id": "sampleA",
        "interval": "chr1"
    },
    [
        {
            "id": "sampleA",
            "repeat": "1",
            "type": "normal",
            "fastq1": "sampleA_rep1_normal_R1.fastq.gz",
            "fastq2": "sampleA_rep1_normal_R2.fastq.gz"
        },
        {
            "id": "sampleA",
            "repeat": "2",
            "type": "normal",
            "fastq1": "sampleA_rep2_normal_R1.fastq.gz",
            "fastq2": "sampleA_rep2_normal_R2.fastq.gz"
        }
    ],
    [
        {
            "id": "sampleA",
            "repeat": "1",
            "type": "tumor",
            "fastq1": "sampleA_rep1_tumor_R1.fastq.gz",
            "fastq2": "sampleA_rep1_tumor_R2.fastq.gz"
        },
        {
            "id": "sampleA",
            "repeat": "2",
            "type": "tumor",
            "fastq1": "sampleA_rep2_tumor_R1.fastq.gz",
            "fastq2": "sampleA_rep2_tumor_R2.fastq.gz"
        }
    ]
]
...
```

Note our data has changed structure. What was previously a list of tuples is now a list of lists of tuples. This is because when we use `groupTuple`, Nextflow creates a new list for each group. This is important to remember when trying to handle the data downstream.

It's possible to use a simpler data structure than this, by separating our the sample information from the sequencing data. We generally refer to this as a `metamap`, but this will be covered in a later side quest. For now, you should just understand that we can group up samples using the `groupTuple` operator and that the data structure will change as a result.

!!! note
[`transpose`](https://www.nextflow.io/docs/latest/reference/operator.html#transpose) is the opposite of groupTuple. It unpacks the items in a channel and flattens them. Try and add `transpose` and undo the grouping we performed above!

# 5.2. Reduce duplication of data

We have a lot of duplicated data in our workflow. Each item in the grouped sample repeats the `id` and `interval` fields. Since this information is available in the metamap, let's just save it once. As a reminder, our data is structured like this:

```groovy
[
    {
        "id": "sampleC",
        "interval": "chr3"
    },
    [
        {
            "id": "sampleC",
            "repeat": "1",
            "type": "normal",
            "fastq1": "sampleC_rep1_normal_R1.fastq.gz",
            "fastq2": "sampleC_rep1_normal_R2.fastq.gz"
        }
    ],
    [
        {
            "id": "sampleC",
            "repeat": "1",
            "type": "tumor",
            "fastq1": "sampleC_rep1_tumor_R1.fastq.gz",
            "fastq2": "sampleC_rep1_tumor_R2.fastq.gz"
        }
    ]
]
```

We could parse the data after grouping to remove the duplication, but this requires us to handle all of the outputs. Instead, we can parse the data before grouping, which will mean the results are never included in the first place.

In the same `map` operator where we isolate the `id` and `interval` fields, we can also grab the `fastq1` and `fastq2` fields for our sample data and _not_ include the `id` and `interval` fields.

_Before:_

```groovy title="main.nf" linenums="30"
    ch_grouped_samples = ch_combined_samples.map { grouping_key, normal, tumor ->
                            [
                                grouping_key.subMap('id', 'interval'),
                                normal,
                                tumor
                            ]

                        }
                        .groupTuple()
                        .dump(tag: 'grouped', pretty: true)
}
```

_After:_

```groovy title="main.nf" linenums="30"
    ch_grouped_samples = ch_combined_samples.map { grouping_key, normal, tumor ->
                            [
                                grouping_key.subMap('id', 'interval'),
                                normal.subMap("fastq1", "fastq2"),
                                tumor.subMap("fastq1", "fastq2")
                            ]

                        }
                        .groupTuple()
                        .dump(tag: 'grouped', pretty: true)
}
```

```bash title="View grouped samples"
nextflow run main.nf -dump-channels grouped
```

```console title="View grouped samples"
 N E X T F L O W   ~  version 24.10.5

Launching `main.nf` [modest_stallman] DSL2 - revision: 5be827a6e8

[DUMP: grouped] [
    {
        "id": "sampleA",
        "interval": "chr1"
    },
    [
        {
            "fastq1": "sampleA_rep1_normal_R1.fastq.gz",
            "fastq2": "sampleA_rep1_normal_R2.fastq.gz"
        },
        {
            "fastq1": "sampleA_rep2_normal_R1.fastq.gz",
            "fastq2": "sampleA_rep2_normal_R2.fastq.gz"
        }
    ],
    [
        {
            "fastq1": "sampleA_rep1_tumor_R1.fastq.gz",
            "fastq2": "sampleA_rep1_tumor_R2.fastq.gz"
        },
        {
            "fastq1": "sampleA_rep2_tumor_R1.fastq.gz",
            "fastq2": "sampleA_rep2_tumor_R2.fastq.gz"
        }
    ]
]
...
```

Now we have a much cleaner output. We can see that the `id` and `interval` fields are only included once, and the `fastq1` and `fastq2` fields are included in the sample data

### Takeaway

In this section, you've learned:

- **Grouping samples**: How to use `groupTuple` to group samples by a specific attribute

You now have a workflow that can split a samplesheet, filter the normal and tumor samples, join them together by sample ID and replicate number, then group them by `id`.

## Summary

In this side quest, you've learned how to split and group data using channels. By modifying the data as it flows through the pipeline, you can construct a pipeline that handles as many samples as possible with no loops or while statements. It gracefully scales to large numbers of samples. Here's what we achieved:

1. **Read in samplesheet with splitCsv**

- Samplesheet details here
- Show with view, then show with dump (is prettier!)

2. **Use filter (and/or map) to manipulate into 2 separate channels**

- Use named closure in map here?
- Show that elements can be in two channels by filtering twice

3. **Join on ID**

- Show that elements can be in two channels by filtering twice

4. **Use groupTuple to group up samples by ID**

- Show that elements can be in two channels by filtering twice

5. **Combine by intervals**

- Show that elements can be in two channels by filtering twice

6. **Group after intervals**

- Show that elements can be in two channels by filtering twice

This approach offers several advantages over writing a pipeline as more standard code, such as using for and while loops:

- We can scale to as many or as few samples as we want with no additional code
- We focus on handling the flow of data through the pipeline, instead of iterating over samples
- We can be as complex or simple as required
- The pipeline becomes more declarative, focusing on what should happen rather than how it should happen
- Nextflow will optimize execution for us by running independent operations in parallel

By mastering these channel operations, you can build flexible, scalable pipelines that handle complex data relationships without resorting to loops or iterative programming. This declarative approach allows Nextflow to optimize execution and parallelize independent operations automatically.

### Key Concepts

1. **Reading Samplesheets**

   ```nextflow
   // Read CSV with header
   Channel.fromPath('samplesheet.csv')
       .splitCsv(header: true)
   ```

2. **Filtering**

   ```nextflow
   // Filter channel based on condition
   channel.filter { it.type == 'tumor' }
   ```

3. **Joining Channels**

   ```nextflow
   // Join two channels by key
   tumor_ch.join(normal_ch)

   // Extract a key and join by this value
   tumor_ch.map { [it.patient_id, it] }
       .join(
          normal_ch.map { [it.patient_id, it] }
        )
   ```

4. **Grouping Data**

   ```nextflow
   // Group by the first element in each tuple
   channel.groupTuple()
   ```

5. **Combining Channels**

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
