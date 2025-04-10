# Splitting and Grouping

Nextflow helps you work with your data in flexible ways. One of the most useful things you can do is split your data into different streams and then group related items back together. This capability is particularly valuable in bioinformatics workflows where you often need to process different sample types separately before combining results for comparison or joint analysis.

Think of it like sorting mail: you might first separate letters by their destination, process each pile differently, and then recombine items going to the same person. In Nextflow, we use special operators to do this with our scientific data.

Nextflow's channel system is at the heart of this flexibility. Channels connect different parts of your workflow, allowing data to flow through your analysis. You can create multiple channels from a single data source, process each channel differently, and then merge channels back together when needed. This approach lets you design workflows that naturally mirror the branching and converging paths of complex bioinformatics analyses.

In this side quest, we'll explore how to split and group data using Nextflow's powerful channel operators. We'll start with a samplesheet containing information about different samples and their associated data. By the end of this side quest, you'll be able to manipulate and combine data streams effectively, making your workflows more efficient and easier to understand.

You will:

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

The samplesheet contains information about different samples and their associated data. In particular, it contains information about the sample's ID, repeat number, type (normal or tumor), and the paths to the BAM files (which don't actually exist, but we will pretend they do).

```console title="samplesheet.csv"
id,repeat,type,bam
sampleA,1,normal,sampleA_r1_normal.bam
sampleA,1,tumor,sampleA_rep1_tumor.bam
sampleB,1,normal,sampleB_rep1_normal.bam
sampleB,1,tumor,sampleB_rep1_tumor.bam
sampleC,1,normal,sampleC_rep1_normal.bam
sampleC,1,tumor,sampleC_rep1_tumor.bam
sampleD,1,normal,sampleD_rep1_normal.bam
sampleD,1,tumor,sampleD_rep1_tumor.bam
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

Launching `main.nf` [deadly_mercator] DSL2 - revision: bd6b0224e9

[id:sampleA, repeat:1, type:normal, bam:sampleA_rep1_normal.bam]
[id:sampleA, repeat:1, type:tumor, bam:sampleA_rep1_tumor.bam]
[id:sampleA, repeat:2, type:normal, bam:sampleB_rep1_normal.bam]
[id:sampleA, repeat:2, type:tumor, bam:sampleB_rep1_tumor.bam]
[id:sampleB, repeat:1, type:normal, bam:sampleC_rep1_normal.bam]
[id:sampleB, repeat:1, type:tumor, bam:sampleC_rep1_tumor.bam]
[id:sampleC, repeat:1, type:normal, bam:sampleD_rep1_normal.bam]
[id:sampleC, repeat:1, type:tumor, bam:sampleD_rep1_tumor.bam]
```

We can see that each row from the CSV file has been converted into a map with keys matching the header row. A map is a key-value data structure similar to dictionaries in Python, objects in JavaScript, or hashes in Ruby.

Each map contains:

- `id`: The sample identifier (sampleA, sampleB, sampleC)
- `repeat`: The replicate number (1 or 2)
- `type`: The sample type (normal or tumor)
- `bam`: Path to the BAM file

This format makes it easy to access specific fields from each sample. For example, we could access the sample ID with `sample.id` or the BAM file path with `sample.bam`. The output above shows each row from the CSV file converted into a map with keys matching the header row. Now that we've successfully read in the samplesheet and have access to the data in each row, we can begin implementing our pipeline logic.

### Takeaway

In this section, you've learned:

- **Reading in a samplesheet**: How to read in a samplesheet with `splitCsv`
- **Viewing data**: How to use `view` to print the data

We now have a channel of maps, each representing a row from the samplesheet. Next, we'll transform this data into a format suitable for our pipeline by extracting metadata and organizing the file paths.

---

## 2. Filter and transform data

### 2.1. Filter data with `filter`

We can use the [`filter` operator](https://www.nextflow.io/docs/latest/operator.html#filter) to filter the data based on a condition. Let's say we only want to process normal samples. We can do this by filtering the data based on the `type` field. Let's insert this before the `view` operator.

_Before:_

```groovy title="main.nf" linenums="1"
workflow {
    ch_samplesheet = Channel.fromPath("./data/samplesheet.csv")
                        .splitCsv(header: true)
                        .view()
}
```

_After:_

```groovy title="main.nf" linenums="1"
workflow {
    ch_samplesheet = Channel.fromPath("./data/samplesheet.csv")
                        .splitCsv(header: true)
                        .filter { sample -> sample.type == 'normal' }
                        .view()
}
```

```bash title="View normal samples"
nextflow run main.nf
```

```console title="View normal samples"
 N E X T F L O W   ~  version 24.10.5

Launching `main.nf` [admiring_brown] DSL2 - revision: 194d61704d

[id:sampleA, repeat:1, type:normal, bam:sampleA_rep1_normal.bam]
[id:sampleA, repeat:2, type:normal, bam:sampleB_rep1_normal.bam]
[id:sampleB, repeat:1, type:normal, bam:sampleC_rep1_normal.bam]
[id:sampleC, repeat:1, type:normal, bam:sampleD_rep1_normal.bam]
```

We have successfully filtered the data to only include normal samples. Let's recap how this works. The `filter` operator takes a closure that is applied to each element in the channel. If the closure returns `true`, the element is included in the output channel. If the closure returns `false`, the element is excluded from the output channel.

In this case, we want to keep only the samples where `sample.type == 'normal'`. In the closure, we use the variable name `sample` to refer to each element in the channel, which then checks if `sample.type` is equal to `'normal'`. If it is, the sample is included in the output channel. If it is not, the sample is excluded from the output channel.

```groovy title="main.nf" linenums="4"
.filter { sample -> sample.type == 'normal' }
```

### 2.2. Filter to just the tumor samples

While useful, we are discarding the tumor samples. Instead, let's rewrite our pipeline to save all the samples to one channel called `ch_samplesheet`, then filter that channel to just the normal samples and save the results to a new channel called `ch_normal_samples`.

_Before:_

```groovy title="main.nf" linenums="1"
workflow {
    ch_samplesheet = Channel.fromPath("./data/samplesheet.csv")
                        .splitCsv(header: true)
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
    ch_normal_samples.view()
}
```

Once again, run the pipeline to see the results:

```bash title="View normal samples"
nextflow run main.nf
```

```console title="View normal samples"
 N E X T F L O W   ~  version 24.10.5

Launching `main.nf` [trusting_poisson] DSL2 - revision: 639186ee74

[id:sampleA, repeat:1, type:normal, bam:sampleA_rep1_normal.bam]
[id:sampleA, repeat:2, type:normal, bam:sampleB_rep1_normal.bam]
[id:sampleB, repeat:1, type:normal, bam:sampleC_rep1_normal.bam]
[id:sampleC, repeat:1, type:normal, bam:sampleD_rep1_normal.bam]
```

Success! We have filtered the data to only include normal samples. Note that we can use view and save the new channel. If we wanted, we still have access to the tumor samples within the `ch_samplesheet` channel. Since we managed it for the normal samples, let's do it for the tumor samples as well:

_Before:_

```groovy title="main.nf" linenums="1"
workflow {
    ch_samplesheet = Channel.fromPath("./data/samplesheet.csv")
                        .splitCsv(header: true)
    ch_normal_samples = ch_samplesheet
                        .filter { sample -> sample.type == 'normal' }
    ch_normal_samples.view()
}
```

_After:_

```groovy title="main.nf" linenums="1"
workflow {
    ch_samplesheet = Channel.fromPath("./data/samplesheet.csv")
                        .splitCsv(header: true)
    ch_normal_samples = ch_samplesheet
                        .filter { sample -> sample.type == 'normal' }
    ch_tumor_samples = ch_samplesheet
                        .filter { sample -> sample.type == 'tumor' }
    ch_normal_samples.view()
    ch_tumor_samples.view()
}
```

```bash title="View tumor samples"
nextflow run main.nf
```

```console title="View tumor samples"
 N E X T F L O W   ~  version 24.10.5

Launching `main.nf` [big_bernard] DSL2 - revision: 897c9e44cc

[id:sampleA, repeat:1, type:normal, bam:sampleA_rep1_normal.bam]
[id:sampleA, repeat:1, type:tumor, bam:sampleA_rep1_tumor.bam]
[id:sampleA, repeat:2, type:normal, bam:sampleB_rep1_normal.bam]
[id:sampleA, repeat:2, type:tumor, bam:sampleB_rep1_tumor.bam]
[id:sampleB, repeat:1, type:normal, bam:sampleC_rep1_normal.bam]
[id:sampleB, repeat:1, type:tumor, bam:sampleC_rep1_tumor.bam]
[id:sampleC, repeat:1, type:normal, bam:sampleD_rep1_normal.bam]
[id:sampleC, repeat:1, type:tumor, bam:sampleD_rep1_tumor.bam]
```

We've managed to separate out the normal and tumor samples into two different channels but they're mixed up when we `view` them in the console! If we want, we can remove one of the `view` operators to see the data in each channel separately. Let's remove the `view` operator for the normal samples:

_Before:_

```groovy title="main.nf" linenums="1"
workflow {
    ch_samplesheet = Channel.fromPath("./data/samplesheet.csv")
                        .splitCsv(header: true)
    ch_normal_samples = ch_samplesheet
                        .filter { sample -> sample.type == 'normal' }
    ch_tumor_samples = ch_samplesheet
                        .filter { sample -> sample.type == 'tumor' }
    ch_normal_samples.view()
    ch_tumor_samples.view()
}
```

_After:_

```groovy title="main.nf" linenums="1"
workflow {
    ch_samplesheet = Channel.fromPath("./data/samplesheet.csv")
                        .splitCsv(header: true)
    ch_normal_samples = ch_samplesheet
                        .filter { sample -> sample.type == 'normal' }
    ch_tumor_samples = ch_samplesheet
                        .filter { sample -> sample.type == 'tumor' }
    ch_tumor_samples.view()
}
```

```bash title="View normal and tumor samples"
nextflow run main.nf
```

```console title="View normal and tumor samples"
 N E X T F L O W   ~  version 24.10.5

Launching `main.nf` [loving_bardeen] DSL2 - revision: 012d38e59f

[id:sampleA, repeat:1, type:tumor, bam:sampleA_rep1_tumor.bam]
[id:sampleA, repeat:2, type:tumor, bam:sampleB_rep1_tumor.bam]
[id:sampleB, repeat:1, type:tumor, bam:sampleC_rep1_tumor.bam]
[id:sampleC, repeat:1, type:tumor, bam:sampleD_rep1_tumor.bam]
```

Note how we can only see the tumor samples in the output. This is because we removed the `view` operator for the normal samples.

### Takeaway

In this section, you've learned:

- **Filtering data**: How to filter data with `filter`
- **Splitting data**: How to split data into different channels based on a condition
- **Viewing data**: How to use `view` to print the data

We've now separated out the normal and tumor samples into two different channels. Next, we'll join the normal and tumor samples on the `id` field.

---

## 3. Join on sample ID

In the previous section, we separated out the normal and tumor samples into two different channels. These could be processed independently using specific processes or workflows based on their type. But what happens when we want to compare the normal and tumor samples from the same patient? At this point, we need to join them back together making sure to match the samples based on their `id` field.

Nextflow includes many methods for combining channels, but in this case the most appropriate operator is [`join`](https://www.nextflow.io/docs/latest/operator.html#join). This acts like a SQL `JOIN` operation, where we specify the key to join on and the type of join to perform.

### 3.1. Use `map` and `join` to combine based on sample ID

If we check the [`join`](https://www.nextflow.io/docs/latest/operator.html#join) documentation, we can see that it joins two channels based on the first item in each tuple. If you don't have the console output still available, let's run the pipeline to check our data structure and see how we need to modify it to join on the `id` field.

```bash title="View normal and tumor samples"
nextflow run main.nf
```

```console title="View normal and tumor samples"
 N E X T F L O W   ~  version 24.10.5

Launching `main.nf` [loving_bardeen] DSL2 - revision: 012d38e59f

[id:sampleA, repeat:1, type:tumor, bam:sampleA_rep1_tumor.bam]
[id:sampleA, repeat:2, type:tumor, bam:sampleB_rep1_tumor.bam]
[id:sampleB, repeat:1, type:tumor, bam:sampleC_rep1_tumor.bam]
[id:sampleC, repeat:1, type:tumor, bam:sampleD_rep1_tumor.bam]
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
    ch_tumor_samples = ch_samplesheet
                        .filter { sample -> sample.type == 'tumor' }
    ch_tumor_samples.view()
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
    ch_tumor_samples = ch_samplesheet
                        .filter { sample -> sample.type == 'tumor' }
                        .map { sample -> [sample.id, sample] }
    ch_normal_samples.view()
    ch_tumor_samples.view()
}
```

```bash title="View normal and tumor samples with ID as element 0"
nextflow run main.nf
```

```console title="View normal and tumor samples with ID as element 0"
 N E X T F L O W   ~  version 24.10.5

Launching `main.nf` [dreamy_sax] DSL2 - revision: 882ae9add4

[sampleA, [id:sampleA, repeat:1, type:normal, bam:sampleA_rep1_normal.bam]]
[sampleA, [id:sampleA, repeat:1, type:tumor, bam:sampleA_rep1_tumor.bam]]
[sampleA, [id:sampleA, repeat:2, type:normal, bam:sampleB_rep1_normal.bam]]
[sampleA, [id:sampleA, repeat:2, type:tumor, bam:sampleB_rep1_tumor.bam]]
[sampleB, [id:sampleB, repeat:1, type:normal, bam:sampleC_rep1_normal.bam]]
[sampleB, [id:sampleB, repeat:1, type:tumor, bam:sampleC_rep1_tumor.bam]]
[sampleC, [id:sampleC, repeat:1, type:normal, bam:sampleD_rep1_normal.bam]]
[sampleC, [id:sampleC, repeat:1, type:tumor, bam:sampleD_rep1_tumor.bam]]
```

It might be subtle, but you should be able to see the first element in each tuple is the `id` field. Now we can use the `join` operator to combine the two channels based on the `id` field.

Once again, we will use `view` to print the joined outputs.

_Before:_

```groovy title="main.nf" linenums="1"
workflow {
    ch_samplesheet = Channel.fromPath("./data/samplesheet.csv")
                        .splitCsv(header: true)
    ch_normal_samples = ch_samplesheet
                        .filter { sample -> sample.type == 'normal' }
                        .map { sample -> [sample.id, sample] }
    ch_tumor_samples = ch_samplesheet
                        .filter { sample -> sample.type == 'tumor' }
                        .map { sample -> [sample.id, sample] }
    ch_normal_samples.view()
    ch_tumor_samples.view()
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
    ch_tumor_samples = ch_samplesheet
                        .filter { sample -> sample.type == 'tumor' }
                        .map { sample -> [sample.id, sample] }
    ch_joined_samples = ch_normal_samples
                        .join(ch_tumor_samples)
    ch_joined_samples.view()
}
```

```bash title="View normal and tumor samples"
nextflow run main.nf
```

```console title="View joined normal and tumor samples"
 N E X T F L O W   ~  version 24.10.5

Launching `main.nf` [elegant_waddington] DSL2 - revision: c552f22069

[sampleA, [id:sampleA, repeat:1, type:normal, bam:sampleA_rep1_normal.bam], [id:sampleA, repeat:1, type:tumor, bam:sampleA_rep1_tumor.bam]]
[sampleA, [id:sampleA, repeat:2, type:normal, bam:sampleB_rep1_normal.bam], [id:sampleA, repeat:2, type:tumor, bam:sampleB_rep1_tumor.bam]]
[sampleB, [id:sampleB, repeat:1, type:normal, bam:sampleC_rep1_normal.bam], [id:sampleB, repeat:1, type:tumor, bam:sampleC_rep1_tumor.bam]]
[sampleC, [id:sampleC, repeat:1, type:normal, bam:sampleD_rep1_normal.bam], [id:sampleC, repeat:1, type:tumor, bam:sampleD_rep1_tumor.bam]]
```

It's a little hard to tell because it's so wide, but you should be able to see the samples have been joined by the `id` field. Each tuple now has the format:

- `id`: The sample ID
- `normal_sample`: The normal sample including type, replicate and path to bam file
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

_Before:_

```groovy title="main.nf" linenums="1"
workflow {
    ch_samplesheet = Channel.fromPath("./data/samplesheet.csv")
                        .splitCsv(header: true)
    ch_normal_samples = ch_samplesheet
                        .filter { sample -> sample.type == 'normal' }
                        .map { sample -> [sample.id, sample] }
    ch_tumor_samples = ch_samplesheet
                        .filter { sample -> sample.type == 'tumor' }
                        .map { sample -> [sample.id, sample] }
    ch_joined_samples = ch_normal_samples
                        .join(ch_tumor_samples)
    ch_joined_samples.view()
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
    ch_tumor_samples = ch_samplesheet
                        .filter { sample -> sample.type == 'tumor' }
                        .map { sample -> [
                                [sample.id, sample.repeat],
                                sample
                            ]
                        }
    ch_joined_samples = ch_normal_samples
                        .join(ch_tumor_samples)
    ch_joined_samples.view()
}
```

Now we should see the join is occurring but using both the `id` and `repeat` fields.

```bash title="View normal and tumor samples"
nextflow run main.nf
```

```console title="View normal and tumor samples"
 N E X T F L O W   ~  version 24.10.5

Launching `main.nf` [prickly_wing] DSL2 - revision: 3bebf22dee

[[sampleA, 1], [id:sampleA, repeat:1, type:normal, bam:sampleA_rep1_normal.bam], [id:sampleA, repeat:1, type:tumor, bam:sampleA_rep1_tumor.bam]]
[[sampleA, 2], [id:sampleA, repeat:2, type:normal, bam:sampleB_rep1_normal.bam], [id:sampleA, repeat:2, type:tumor, bam:sampleB_rep1_tumor.bam]]
[[sampleB, 1], [id:sampleB, repeat:1, type:normal, bam:sampleC_rep1_normal.bam], [id:sampleB, repeat:1, type:tumor, bam:sampleC_rep1_tumor.bam]]
[[sampleC, 1], [id:sampleC, repeat:1, type:normal, bam:sampleD_rep1_normal.bam], [id:sampleC, repeat:1, type:tumor, bam:sampleD_rep1_tumor.bam]]
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
    ch_tumor_samples = ch_samplesheet
                        .filter { sample -> sample.type == 'tumor' }
                        .map { sample -> [
                                [sample.id, sample.repeat],
                                sample
                            ]
                        }
    ch_joined_samples = ch_normal_samples
                        .join(ch_tumor_samples)
    ch_joined_samples.view()
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
    ch_tumor_samples = ch_samplesheet
                        .filter { sample -> sample.type == 'tumor' }
                        .map { sample -> [
                                sample.subMap(['id', 'repeat']),
                                sample
                            ]
                        }
    ch_joined_samples = ch_normal_samples
                        .join(ch_tumor_samples)
    ch_joined_samples.view()
}
```

```bash title="View normal and tumor samples"
nextflow run main.nf
```

```console title="View normal and tumor samples"
 N E X T F L O W   ~  version 24.10.5

Launching `main.nf` [curious_hopper] DSL2 - revision: 90283e523d

[[id:sampleA, repeat:1], [id:sampleA, repeat:1, type:normal, bam:sampleA_rep1_normal.bam], [id:sampleA, repeat:1, type:tumor, bam:sampleA_rep1_tumor.bam]]
[[id:sampleA, repeat:2], [id:sampleA, repeat:2, type:normal, bam:sampleB_rep1_normal.bam], [id:sampleA, repeat:2, type:tumor, bam:sampleB_rep1_tumor.bam]]
[[id:sampleB, repeat:1], [id:sampleB, repeat:1, type:normal, bam:sampleC_rep1_normal.bam], [id:sampleB, repeat:1, type:tumor, bam:sampleC_rep1_tumor.bam]]
[[id:sampleC, repeat:1], [id:sampleC, repeat:1, type:normal, bam:sampleD_rep1_normal.bam], [id:sampleC, repeat:1, type:tumor, bam:sampleD_rep1_tumor.bam]]
```

Now we have a new joining key that not only includes the `id` and `repeat` fields but also retains the field names so we can access them later by name, e.g. `sample.id` and `sample.repeat`.

### 3.4. Use a named closure in map

Since we are re-using the same map in multiple places, we run the risk of introducing errors if we accidentally change the map in one place but not the other. To avoid this, we can use a named closure in the map. A named closure allows us to make a reusable function we can call later within a map.

To do so, first we define the closure as a new variable:

_Before:_

```groovy title="main.nf" linenums="1"
workflow {
    ch_samplesheet = Channel.fromPath("./data/samplesheet.csv")
                        .splitCsv(header: true)
```

_After:_

```groovy title="main.nf" linenums="1"
workflow {
    getSampleIdAndReplicate = { sample -> [ sample.subMap(['id', 'repeat']), sample ] }
    ch_samplesheet = Channel.fromPath("./data/samplesheet.csv")
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
    ch_tumor_samples = ch_samplesheet
                        .filter { sample -> sample.type == "tumor" }
                        .map { sample -> [
                            sample.subMap(['id', 'repeat']),
                            sample
                          ]
                        }
```

_After:_

```groovy title="main.nf" linenums="5"
    ch_normal_samples = ch_samplesheet
                        .filter { sample -> sample.type == 'normal' }
                        .map ( getSampleIdAndReplicate )
    ch_tumor_samples = ch_samplesheet
                        .filter { sample -> sample.type == "tumor" }
                        .map ( getSampleIdAndReplicate )
```

!!! note

    The `map` operator has switched from using `{ }` to using `( )` to pass the closure as an argument. This is because the `map` operator expects a closure as an argument and `{ }` is used to define an anonymous closure. When calling a named closure, use the `( )` syntax.

```bash title="View normal and tumor samples"
nextflow run main.nf
```

```console title="View normal and tumor samples"
 N E X T F L O W   ~  version 24.10.5

Launching `main.nf` [angry_meninsky] DSL2 - revision: 2edc226b1d

[[id:sampleA, repeat:1], [id:sampleA, repeat:1, type:normal, bam:sampleA_rep1_normal.bam], [id:sampleA, repeat:1, type:tumor, bam:sampleA_rep1_tumor.bam]]
[[id:sampleA, repeat:2], [id:sampleA, repeat:2, type:normal, bam:sampleB_rep1_normal.bam], [id:sampleA, repeat:2, type:tumor, bam:sampleB_rep1_tumor.bam]]
[[id:sampleB, repeat:1], [id:sampleB, repeat:1, type:normal, bam:sampleC_rep1_normal.bam], [id:sampleB, repeat:1, type:tumor, bam:sampleC_rep1_tumor.bam]]
[[id:sampleC, repeat:1], [id:sampleC, repeat:1, type:normal, bam:sampleD_rep1_normal.bam], [id:sampleC, repeat:1, type:tumor, bam:sampleD_rep1_tumor.bam]]
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
    "bam": "sampleC_rep1_normal.bam"
  ],
  [
    "id": "sampleC",
    "repeat": "1",
    "type": "tumor",
    "bam": "sampleC_rep1_tumor.bam"
  ]
]
```

Since the `id` and `repeat` fields are available in the grouping key, let's remove them from the sample data to avoid duplication. We can do this by using the `subMap` method to create a new map with only the `type` and `bam` fields. This approach allows us to maintain all necessary information while eliminating redundancy in our data structure.

_Before:_

```groovy title="main.nf" linenums="15"
workflow {
    getSampleIdAndReplicate = { sample -> [ sample.subMap(['id', 'repeat']), sample ] }
    ch_samplesheet = Channel.fromPath("./data/samplesheet.csv")
                        .splitCsv(header: true)
    ch_normal_samples = ch_samplesheet
                        .filter { sample -> sample.type == 'normal' }
                        .map ( getSampleIdAndReplicate )
    ch_tumor_samples = ch_samplesheet
                        .filter { sample -> sample.type == "tumor" }
                        .map ( getSampleIdAndReplicate )
    ch_joined_samples = ch_normal_samples
                        .join(ch_tumor_samples)
    ch_joined_samples.view()
}
```

_After:_

```groovy title="main.nf" linenums="15"
workflow {
    getSampleIdAndReplicate = { sample ->
                                  [
                                    sample.subMap(['id', 'repeat']),
                                    sample.subMap(['type', 'bam'])
                                  ]
                              }
    ch_samplesheet = Channel.fromPath("./data/samplesheet.csv")
                        .splitCsv(header: true)
    ch_normal_samples = ch_samplesheet
                        .filter { sample -> sample.type == 'normal' }
                        .map ( getSampleIdAndReplicate )
    ch_tumor_samples = ch_samplesheet
                        .filter { sample -> sample.type == "tumor" }
                        .map ( getSampleIdAndReplicate )
    ch_joined_samples = ch_normal_samples
                        .join(ch_tumor_samples)
    ch_joined_samples.view()
}
```

Now, when the closure returns the tuple, the first element is the `id` and `repeat` fields and the second element is the `type` and `bam` fields. We have effectively removed the `id` and `repeat` fields from the sample data and uniquely store them in the grouping key. This approach eliminates redundancy while maintaining all necessary information.

```bash title="View deduplicated data"
nextflow run main.nf
```

```console title="View deduplicated data"
 N E X T F L O W   ~  version 24.10.5

Launching `main.nf` [trusting_pike] DSL2 - revision: 09d3c7a81b

[[id:sampleA, repeat:1], [type:normal, bam:sampleA_rep1_normal.bam], [type:tumor, bam:sampleA_rep1_tumor.bam]]
[[id:sampleA, repeat:2], [type:normal, bam:sampleB_rep1_normal.bam], [type:tumor, bam:sampleB_rep1_tumor.bam]]
[[id:sampleB, repeat:1], [type:normal, bam:sampleC_rep1_normal.bam], [type:tumor, bam:sampleC_rep1_tumor.bam]]
[[id:sampleC, repeat:1], [type:normal, bam:sampleD_rep1_normal.bam], [type:tumor, bam:sampleD_rep1_tumor.bam]]
```

We can see we only state the `id` and `repeat` fields once in the grouping key and we have the `type` and `bam` fields in the sample data. We haven't lost any information but we managed to make our channel contents more succinct.

### Takeaway

In this section, you've learned:

- **Manipulating Tuples**: How to use `map` to isolate a field in a tuple
- **Joining Tuples**: How to use `join` to combine tuples based on the first field
- **Creating Joining Keys**: How to use `subMap` to create a new joining key
- **Named Closures**: How to use a named closure in map
- **Deduplicating Data**: How to remove duplicate data from the channel
  You now have a workflow that can split a samplesheet, filter the normal and tumor samples, join them together by sample ID and replicate number, then print the results.

This is a common pattern in bioinformatics workflows where you need to match up samples after processing independently, so it is a useful skill. Next, we will look at repeating a sample multiple times.

## 4. Spread samples over intervals

A key pattern in bioinformatics workflows is distributing analysis across genomic regions. For instance, variant calling can be parallelized by dividing the genome into intervals (like chromosomes or smaller regions). This parallelization strategy significantly improves pipeline efficiency by distributing computational load across multiple cores or nodes, reducing overall execution time.

In the following section, we'll demonstrate how to distribute our sample data across multiple genomic intervals. We'll pair each sample with every interval, allowing parallel processing of different genomic regions. This will multiply our dataset size by the number of intervals, creating multiple independent analysis units that can be brought back together later.

### 4.1. Spread samples over intervals using `combine`

Let's start by creating a channel of intervals. To keep life simple, we will just use 3 intervals we will manually define. In a real workflow, you could read these in from a file input or even create a channel with lots of interval files.

_Before:_

```groovy title="main.nf" linenums="15"
    ch_joined_samples.view()
}
```

_After:_

```groovy title="main.nf" linenums="15"
    ch_intervals = Channel.of('chr1', 'chr2', 'chr3')
}
```

Now remember, we want to repeat each sample for each interval. This is sometimes referred to as the Cartesian product of the samples and intervals. We can achieve this by using the [`combine` operator](https://www.nextflow.io/docs/latest/operator.html#combine). This will take every item from channel 1 and repeat it for each item in channel 2. Let's add a combine operator to our workflow:

_Before:_

```groovy title="main.nf" linenums="15"
    ch_intervals = Channel.of('chr1', 'chr2', 'chr3')
}
```

_After:_

```groovy title="main.nf" linenums="15"
    ch_intervals = Channel.of('chr1', 'chr2', 'chr3')

    ch_combined_samples = ch_joined_samples.combine(ch_intervals)
                        .view()
}
```

Now let's run it and see what happens:

```bash title="View combined samples"
nextflow run main.nf
```

```console title="View combined samples"
 N E X T F L O W   ~  version 24.10.5

Launching `main.nf` [soggy_fourier] DSL2 - revision: fa8f5edb22

[[id:sampleA, repeat:1], [type:normal, bam:sampleA_rep1_normal.bam], [type:tumor, bam:sampleA_rep1_tumor.bam], chr1]
[[id:sampleA, repeat:1], [type:normal, bam:sampleA_rep1_normal.bam], [type:tumor, bam:sampleA_rep1_tumor.bam], chr2]
[[id:sampleA, repeat:1], [type:normal, bam:sampleA_rep1_normal.bam], [type:tumor, bam:sampleA_rep1_tumor.bam], chr3]
[[id:sampleA, repeat:2], [type:normal, bam:sampleB_rep1_normal.bam], [type:tumor, bam:sampleB_rep1_tumor.bam], chr1]
[[id:sampleA, repeat:2], [type:normal, bam:sampleB_rep1_normal.bam], [type:tumor, bam:sampleB_rep1_tumor.bam], chr2]
[[id:sampleA, repeat:2], [type:normal, bam:sampleB_rep1_normal.bam], [type:tumor, bam:sampleB_rep1_tumor.bam], chr3]
[[id:sampleB, repeat:1], [type:normal, bam:sampleC_rep1_normal.bam], [type:tumor, bam:sampleC_rep1_tumor.bam], chr1]
[[id:sampleB, repeat:1], [type:normal, bam:sampleC_rep1_normal.bam], [type:tumor, bam:sampleC_rep1_tumor.bam], chr2]
[[id:sampleB, repeat:1], [type:normal, bam:sampleC_rep1_normal.bam], [type:tumor, bam:sampleC_rep1_tumor.bam], chr3]
[[id:sampleC, repeat:1], [type:normal, bam:sampleD_rep1_normal.bam], [type:tumor, bam:sampleD_rep1_tumor.bam], chr1]
[[id:sampleC, repeat:1], [type:normal, bam:sampleD_rep1_normal.bam], [type:tumor, bam:sampleD_rep1_tumor.bam], chr2]
[[id:sampleC, repeat:1], [type:normal, bam:sampleD_rep1_normal.bam], [type:tumor, bam:sampleD_rep1_tumor.bam], chr3]
```

Success! We have repeated every sample for every single interval in our 3 interval list. We've effectively tripled the number of items in our channel. It's a little hard to read though, so in the next section we will tidy it up.

### 4.2. Organise the channel

We can use the `map` operator to tidy and refactor our sample data so it's easier to understand. Let's move the intervals string to the joining map at the first element.

_Before:_

```groovy title="main.nf" linenums="19"
    ch_combined_samples = ch_joined_samples.combine(ch_intervals)
                        .view()
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
                        .view()
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
nextflow run main.nf
```

```console title="View combined samples"
 N E X T F L O W   ~  version 24.10.5

Launching `main.nf` [sad_hawking] DSL2 - revision: 1f6f6250cd

[[id:sampleA, repeat:1, interval:chr1], [type:normal, bam:sampleA_rep1_normal.bam], [type:tumor, bam:sampleA_rep1_tumor.bam]]
[[id:sampleA, repeat:1, interval:chr2], [type:normal, bam:sampleA_rep1_normal.bam], [type:tumor, bam:sampleA_rep1_tumor.bam]]
[[id:sampleA, repeat:1, interval:chr3], [type:normal, bam:sampleA_rep1_normal.bam], [type:tumor, bam:sampleA_rep1_tumor.bam]]
[[id:sampleA, repeat:2, interval:chr1], [type:normal, bam:sampleB_rep1_normal.bam], [type:tumor, bam:sampleB_rep1_tumor.bam]]
[[id:sampleA, repeat:2, interval:chr2], [type:normal, bam:sampleB_rep1_normal.bam], [type:tumor, bam:sampleB_rep1_tumor.bam]]
[[id:sampleA, repeat:2, interval:chr3], [type:normal, bam:sampleB_rep1_normal.bam], [type:tumor, bam:sampleB_rep1_tumor.bam]]
[[id:sampleB, repeat:1, interval:chr1], [type:normal, bam:sampleC_rep1_normal.bam], [type:tumor, bam:sampleC_rep1_tumor.bam]]
[[id:sampleB, repeat:1, interval:chr2], [type:normal, bam:sampleC_rep1_normal.bam], [type:tumor, bam:sampleC_rep1_tumor.bam]]
[[id:sampleB, repeat:1, interval:chr3], [type:normal, bam:sampleC_rep1_normal.bam], [type:tumor, bam:sampleC_rep1_tumor.bam]]
[[id:sampleC, repeat:1, interval:chr1], [type:normal, bam:sampleD_rep1_normal.bam], [type:tumor, bam:sampleD_rep1_tumor.bam]]
[[id:sampleC, repeat:1, interval:chr2], [type:normal, bam:sampleD_rep1_normal.bam], [type:tumor, bam:sampleD_rep1_tumor.bam]]
[[id:sampleC, repeat:1, interval:chr3], [type:normal, bam:sampleD_rep1_normal.bam], [type:tumor, bam:sampleD_rep1_tumor.bam]]
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
                        .view()
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

    ch_grouped_samples = ch_combined_samples.map { grouping_key, normal, tumor ->
                            [
                                grouping_key.subMap('id', 'interval'),
                                normal,
                                tumor
                            ]

                        }
                        .view()
}
```

Let's run it again and check the channel contents:

```bash title="View grouped samples"
nextflow run main.nf
```

```console title="View grouped samples"
 N E X T F L O W   ~  version 24.10.5

Launching `main.nf` [loving_escher] DSL2 - revision: 3adccba898

[[id:sampleA, interval:chr1], [type:normal, bam:sampleA_rep1_normal.bam], [type:tumor, bam:sampleA_rep1_tumor.bam]]
[[id:sampleA, interval:chr2], [type:normal, bam:sampleA_rep1_normal.bam], [type:tumor, bam:sampleA_rep1_tumor.bam]]
[[id:sampleA, interval:chr3], [type:normal, bam:sampleA_rep1_normal.bam], [type:tumor, bam:sampleA_rep1_tumor.bam]]
[[id:sampleA, interval:chr1], [type:normal, bam:sampleB_rep1_normal.bam], [type:tumor, bam:sampleB_rep1_tumor.bam]]
[[id:sampleA, interval:chr2], [type:normal, bam:sampleB_rep1_normal.bam], [type:tumor, bam:sampleB_rep1_tumor.bam]]
[[id:sampleA, interval:chr3], [type:normal, bam:sampleB_rep1_normal.bam], [type:tumor, bam:sampleB_rep1_tumor.bam]]
[[id:sampleB, interval:chr1], [type:normal, bam:sampleC_rep1_normal.bam], [type:tumor, bam:sampleC_rep1_tumor.bam]]
[[id:sampleB, interval:chr2], [type:normal, bam:sampleC_rep1_normal.bam], [type:tumor, bam:sampleC_rep1_tumor.bam]]
[[id:sampleB, interval:chr3], [type:normal, bam:sampleC_rep1_normal.bam], [type:tumor, bam:sampleC_rep1_tumor.bam]]
[[id:sampleC, interval:chr1], [type:normal, bam:sampleD_rep1_normal.bam], [type:tumor, bam:sampleD_rep1_tumor.bam]]
[[id:sampleC, interval:chr2], [type:normal, bam:sampleD_rep1_normal.bam], [type:tumor, bam:sampleD_rep1_tumor.bam]]
[[id:sampleC, interval:chr3], [type:normal, bam:sampleD_rep1_normal.bam], [type:tumor, bam:sampleD_rep1_tumor.bam]]
```

We can see that we have successfully isolated the `id` and `interval` fields, but not grouped the samples yet.

!!! note

    We are discarding the `replicate` field here. This is because we don't need it for further downstream processing. After completing this tutorial, see if you can include it without affecting the later grouping!

Let's now group the samples by this new grouping element, using the [`groupTuple` operator](https://www.nextflow.io/docs/latest/operator.html#grouptuple).

_Before:_

```groovy title="main.nf" linenums="30"
    ch_grouped_samples = ch_combined_samples.map { grouping_key, normal, tumor ->
                            [
                                grouping_key.subMap('id', 'interval'),
                                normal,
                                tumor
                            ]

                        }
                        .view()
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
                        .view()
}
```

Simple, huh? We just added a single line of code. Let's see what happens when we run it:

```bash title="View grouped samples"
nextflow run main.nf
```

```console title="View grouped samples"
 N E X T F L O W   ~  version 24.10.5

Launching `main.nf` [festering_almeida] DSL2 - revision: 78988949e3

[[id:sampleA, interval:chr1], [[type:normal, bam:sampleA_rep1_normal.bam], [type:normal, bam:sampleB_rep1_normal.bam]], [[type:tumor, bam:sampleA_rep1_tumor.bam], [type:tumor, bam:sampleB_rep1_tumor.bam]]]
[[id:sampleA, interval:chr2], [[type:normal, bam:sampleA_rep1_normal.bam], [type:normal, bam:sampleB_rep1_normal.bam]], [[type:tumor, bam:sampleA_rep1_tumor.bam], [type:tumor, bam:sampleB_rep1_tumor.bam]]]
[[id:sampleA, interval:chr3], [[type:normal, bam:sampleA_rep1_normal.bam], [type:normal, bam:sampleB_rep1_normal.bam]], [[type:tumor, bam:sampleA_rep1_tumor.bam], [type:tumor, bam:sampleB_rep1_tumor.bam]]]
[[id:sampleB, interval:chr1], [[type:normal, bam:sampleC_rep1_normal.bam]], [[type:tumor, bam:sampleC_rep1_tumor.bam]]]
[[id:sampleB, interval:chr2], [[type:normal, bam:sampleC_rep1_normal.bam]], [[type:tumor, bam:sampleC_rep1_tumor.bam]]]
[[id:sampleB, interval:chr3], [[type:normal, bam:sampleC_rep1_normal.bam]], [[type:tumor, bam:sampleC_rep1_tumor.bam]]]
[[id:sampleC, interval:chr1], [[type:normal, bam:sampleD_rep1_normal.bam]], [[type:tumor, bam:sampleD_rep1_tumor.bam]]]
[[id:sampleC, interval:chr2], [[type:normal, bam:sampleD_rep1_normal.bam]], [[type:tumor, bam:sampleD_rep1_tumor.bam]]]
[[id:sampleC, interval:chr3], [[type:normal, bam:sampleD_rep1_normal.bam]], [[type:tumor, bam:sampleD_rep1_tumor.bam]]]
```

Note our data has changed structure. What was previously a list of tuples is now a list of lists of tuples. This is because when we use `groupTuple`, Nextflow creates a new list for each group. This is important to remember when trying to handle the data downstream.

It's possible to use a simpler data structure than this, by separating our the sample information from the sequencing data. We generally refer to this as a `metamap`, but this will be covered in a later side quest. For now, you should just understand that we can group up samples using the `groupTuple` operator and that the data structure will change as a result.

!!! note

    [`transpose`](https://www.nextflow.io/docs/latest/reference/operator.html#transpose) is the opposite of groupTuple. It unpacks the items in a channel and flattens them. Try and add `transpose` and undo the grouping we performed above!

# 5.2. Reorganise the data

Let's consider the inputs to a typical Nextflow process. Generally, inputs can be in the form of values or files. In this example, we have a set of values for sample information (`id` and `interval`) and a set of files for sequencing data (`normal` and `tumor`). The `input` block of a process might look like this:

```groovy title="main.nf"
    input:
      tuple val(sampleInfo), path(normalBam), path(tumorBam)
```

Currently, our data structure isn't optimal for this. We have a tuple where the first element is a map and the second and third elements are lists of maps. Here is the current data structure:

```groovy
[
    [id:sampleA, interval:chr1],
    [[type:normal, bam:sampleA_rep1_normal.bam], [type:normal, bam:sampleB_rep1_normal.bam]],
    [[type:tumor, bam:sampleA_rep1_tumor.bam], [type:tumor, bam:sampleB_rep1_tumor.bam]]
],
[
    [id:sampleA, interval:chr1],
    [[type:normal, bam:sampleA_rep1_normal.bam]]
    [[type:tumor, bam:sampleA_rep1_tumor.bam]]
]
```

What we need is to flatten the data structure so that the second and third elements are lists of BAM files - essentially removing the `type` field and simplifying the structure.

Once again, we'll employ the `map` operator to manipulate our data structure as it flows through the channel. This time, we'll take the lists of data associated with each BAM file and specifically extract the `bam` field, dropping the `type` field in the process.

To operate on a list of maps, we can use the [`collect` method](<https://docs.groovy-lang.org/2.4.3/html/groovy-jdk/java/util/Collection.html#collect(groovy.lang.Closure)>), a built-in method in Groovy which iterates over each element in the list and applies the closure to it. You can see examples in the Nextflow documentation [here](https://www.nextflow.io/docs/latest/script.html#closures). In our case, we'll extract the `bam` field from each map and return a list of just BAM file paths. Our closure will look like this:

```groovy
{ bam_data -> bam_data.bam }
```

Note this is conceptually similar but distinct to the Nextflow [`collect` operator](https://www.nextflow.io/docs/latest/reference/operator.html#collect). This method takes a closure that will be applied to each element in the list. In this case, we want to extract the `bam` field from each map.

Let's append our map to the end of our pipeline and show the resulting data structure:

_Before:_

```groovy title="main.nf" linenums="38"
                        .groupTuple()
                        .view()
}
```

_After:_

```groovy title="main.nf" linenums="38"
                        .groupTuple()
                        .map { sample_info, normal, tumor ->
                            [
                                sample_info,
                                normal.collect { bam_data -> bam_data.bam },
                                tumor.collect { bam_data -> bam_data.bam }
                            ]
                        }
                        .view()
}
```

```bash title="View flattened samples"
nextflow run main.nf
```

```console title=""

 N E X T F L O W   ~  version 24.10.5

Launching `main.nf` [compassionate_moriondo] DSL2 - revision: 71a8c35ed9

[[id:sampleA, interval:chr1], [sampleA_rep1_normal.bam, sampleB_rep1_normal.bam], [sampleA_rep1_tumor.bam, sampleB_rep1_tumor.bam]]
[[id:sampleA, interval:chr2], [sampleA_rep1_normal.bam, sampleB_rep1_normal.bam], [sampleA_rep1_tumor.bam, sampleB_rep1_tumor.bam]]
[[id:sampleA, interval:chr3], [sampleA_rep1_normal.bam, sampleB_rep1_normal.bam], [sampleA_rep1_tumor.bam, sampleB_rep1_tumor.bam]]
[[id:sampleB, interval:chr1], [sampleC_rep1_normal.bam], [sampleC_rep1_tumor.bam]]
[[id:sampleB, interval:chr2], [sampleC_rep1_normal.bam], [sampleC_rep1_tumor.bam]]
[[id:sampleB, interval:chr3], [sampleC_rep1_normal.bam], [sampleC_rep1_tumor.bam]]
[[id:sampleC, interval:chr1], [sampleD_rep1_normal.bam], [sampleD_rep1_tumor.bam]]
[[id:sampleC, interval:chr2], [sampleD_rep1_normal.bam], [sampleD_rep1_tumor.bam]]
[[id:sampleC, interval:chr3], [sampleD_rep1_normal.bam], [sampleD_rep1_tumor.bam]]
```

Note how the channel is now structured as a 3-part tuple:

- `sample_info` is a map of sample information
- `normal` is a list of BAM file paths
- `tumor` is a list of BAM file paths

`groupTuple` is a powerful operator but can generate complex data structures. It's important to understand how the data structure changes as it flows through the pipeline so you can manipulate it as needed. Using a `map` at the end of a pipeline helps refine the output into a structure that fits our processes pipeline.

## 5.3. Simplify the data

One issue we have faced in this pipeline is that we have a moderately complicated data structure which we have had to coerce throughout the pipeline. What if we could simplify it at the start? Then we would only handle the relevant fields in the pipeline and avoid the need for the final `map` operator.

If we parse the data right at the start of our pipeline to _only_ include the `bam` field, we can avoid passing the `type` field through the pipeline which makes our entire pipeline cleaner while retaining the same functionality:

_Before:_

```groovy title="main.nf" linenums="1"
workflow {
    getSampleIdAndReplicate = { sample ->
                                  [
                                    sample.subMap(['id', 'repeat']),
                                    sample.subMap(['type', 'bam'])
                                  ]
                              }
```

_After:_

```groovy title="main.nf" linenums="1"
workflow {
    getSampleIdAndReplicate = { sample ->
                                  [
                                    sample.subMap(['id', 'repeat']),
                                    sample.bam
                                  ]
                              }
```

A reminder, this will select only the BAM files once we have separated the channels into normal and tumor. We are losing the `type` field, but we know which samples are normal and tumor because they have been filtered and the channel should only contain one type per sample. Once we have done this we can remove the `map` operator from the end of the pipeline:

_Before:_

```groovy title="main.nf" linenums="38"
                        .groupTuple()
                        .map { sample_info, normal, tumor ->
                            [
                                sample_info,
                                normal.collect { bam_data -> bam_data.bam },
                                tumor.collect { bam_data -> bam_data.bam }
                            ]
                        }
                        .view()
}
```

_After:_

```groovy title="main.nf" linenums="38"
                        .groupTuple()
                        .view()
}
```

Sometimes parsing data earlier in the pipeline is the right choice to avoid complicated code.

```bash title="View flattened samples"
nextflow run main.nf
```

```console title="View flattened samples"
 N E X T F L O W   ~  version 24.10.5

Launching `main.nf` [reverent_angela] DSL2 - revision: 656a31b305

[[id:sampleA, interval:chr1], [sampleA_rep1_normal.bam, sampleB_rep1_normal.bam], [sampleA_rep1_tumor.bam, sampleB_rep1_tumor.bam]]
[[id:sampleA, interval:chr2], [sampleA_rep1_normal.bam, sampleB_rep1_normal.bam], [sampleA_rep1_tumor.bam, sampleB_rep1_tumor.bam]]
[[id:sampleA, interval:chr3], [sampleA_rep1_normal.bam, sampleB_rep1_normal.bam], [sampleA_rep1_tumor.bam, sampleB_rep1_tumor.bam]]
[[id:sampleB, interval:chr1], [sampleC_rep1_normal.bam], [sampleC_rep1_tumor.bam]]
[[id:sampleB, interval:chr2], [sampleC_rep1_normal.bam], [sampleC_rep1_tumor.bam]]
[[id:sampleB, interval:chr3], [sampleC_rep1_normal.bam], [sampleC_rep1_tumor.bam]]
[[id:sampleC, interval:chr1], [sampleD_rep1_normal.bam], [sampleD_rep1_tumor.bam]]
[[id:sampleC, interval:chr2], [sampleD_rep1_normal.bam], [sampleD_rep1_tumor.bam]]
[[id:sampleC, interval:chr3], [sampleD_rep1_normal.bam], [sampleD_rep1_tumor.bam]]
```

### Takeaway

In this section, you've learned:

- **Grouping samples**: How to use `groupTuple` to group samples by a specific attribute
- **Flattening data structure**: How to use `map` to flatten the data structure
  You now have a workflow that can split a samplesheet, filter the normal and tumor samples, join them together by sample ID and replicate number, then group them by `id`.

## Summary

In this side quest, you've learned how to split and group data using channels. By modifying the data as it flows through the pipeline, you can construct a pipeline that handles as many samples as possible with no loops or while statements. It gracefully scales to large numbers of samples. Here's what we achieved:

1. **Read in samplesheet with splitCsv**

- Samplesheet details here
- Show with view, then show with view

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
