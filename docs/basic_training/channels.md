---
description: Basic Nextflow Training Workshop
---

# Channels

Channels are a key data structure of Nextflow that allows the implementation of reactive-functional oriented computational workflows based on the [Dataflow](https://en.wikipedia.org/wiki/Dataflow_programming) programming paradigm.

They are used to logically connect tasks to each other or to implement functional style data transformations.

<figure class="excalidraw">
--8<-- "docs/basic_training/img/channel-files.excalidraw.svg"
</figure>

## Channel types

Nextflow distinguishes two different kinds of channels: **queue** channels and **value** channels.

### Queue channel

A _queue_ channel is an _asynchronous_ unidirectional _FIFO_ queue that connects two processes or operators.

-   _asynchronous_ means that operations are non-blocking.
-   _unidirectional_ means that data flows from a producer to a consumer.
-   _FIFO_ means that the data is guaranteed to be delivered in the same order as it is produced. First In, First Out.

A queue channel is implicitly created by process output definitions or using channel factories such as [Channel.of](https://www.nextflow.io/docs/latest/channel.html#of) or [Channel.fromPath](https://www.nextflow.io/docs/latest/channel.html#frompath).

Try the following snippets:

!!! info ""

    Click the :material-plus-circle: icons in the code for explanations.

```groovy linenums="1"
ch = Channel.of(1, 2, 3)
println(ch) // (1)!
ch.view() // (2)!
```

1. Use the built-in print line function `println` to print the `ch` channel
2. Apply the `view` channel operator to the `ch` channel prints each item emitted by the channels

!!! exercise

    Try to execute this snippet. You can do that by creating a new `.nf` file or by editing an already existing `.nf` file.

    ```groovy linenums="1"
    ch = Channel.of(1, 2, 3)
    ch.view()
    ```

### Value channels

A **value** channel (a.k.a. singleton channel) by definition is bound to a single value and it can be read unlimited times without consuming its contents. A `value` channel is created using the [value](https://www.nextflow.io/docs/latest/channel.html#value) channel factory or by operators returning a single value, such as [first](https://www.nextflow.io/docs/latest/operator.html#first), [last](https://www.nextflow.io/docs/latest/operator.html#last), [collect](https://www.nextflow.io/docs/latest/operator.html#operator-collect), [count](https://www.nextflow.io/docs/latest/operator.html#operator-count), [min](https://www.nextflow.io/docs/latest/operator.html#operator-min), [max](https://www.nextflow.io/docs/latest/operator.html#operator-max), [reduce](https://www.nextflow.io/docs/latest/operator.html#operator-reduce), and [sum](https://www.nextflow.io/docs/latest/operator.html#operator-sum).

To better understand the difference between value and queue channels, save the snippet below as `example.nf`.

```groovy linenums="1" title="example.nf" linenums="1"
ch1 = Channel.of(1, 2, 3)
ch2 = Channel.of(1)

process SUM {
    input:
    val x
    val y

    output:
    stdout

    script:
    """
    echo \$(($x+$y))
    """
}

workflow {
    SUM(ch1, ch2).view()
}
```

When you run the script, it prints only 2, as you can see below:

```console
2
```

To understand why, we can inspect the queue channel and running Nextflow with DSL1 gives us a more explicit comprehension of what is behind the curtains.

```groovy linenums="1"
ch1 = Channel.of(1)
println ch1
```

```console
$ nextflow run example.nf -dsl1
...
DataflowQueue(queue=[DataflowVariable(value=1), DataflowVariable(value=groovyx.gpars.dataflow.operator.PoisonPill@34be065a)])
```

We have the value 1 as the single element of our queue channel and a poison pill, which will tell the process that there’s nothing left to be consumed. That’s why we only have one output for the example above, which is 2. Let’s inspect a value channel now.

```groovy linenums="1"
ch1 = Channel.value(1)
println ch1
```

```console
$ nextflow run example.nf -dsl1
...
DataflowVariable(value=1)
```

There is no poison pill, and that’s why we get a different output with the code below, where `ch2` is turned into a value channel through the `first` operator.

```groovy linenums="1"
ch1 = Channel.of(1, 2, 3)
ch2 = Channel.of(1)

process SUM {
    input:
    val x
    val y

    output:
    stdout

    script:
    """
    echo \$(($x+$y))
    """
}

workflow {
    SUM(ch1, ch2.first()).view()
}
```

```console title="Output"
4

3

2
```

Besides, in many situations, Nextflow will implicitly convert variables to value channels when they are used in a process invocation. For example, when you invoke a process with a pipeline parameter (`params.example`) which has a string value, it is automatically cast into a value channel.

## Channel factories

These are Nextflow commands for creating channels that have implicit expected inputs and functions.

### `value()`

The `value` channel factory is used to create a _value_ channel. An optional not `null` argument can be specified to bind the channel to a specific value. For example:

```groovy linenums="1"
ch1 = Channel.value() // (1)!
ch2 = Channel.value('Hello there') // (2)!
ch3 = Channel.value([1, 2, 3, 4, 5]) // (3)!
```

1. Creates an _empty_ value channel
2. Creates a value channel and binds a string to it
3. Creates a value channel and binds a list object to it that will be emitted as a sole emission

### `of()`

The factory `Channel.of` allows the creation of a queue channel with the values specified as arguments.

```groovy linenums="1"
ch = Channel.of(1, 3, 5, 7)
ch.view { "value: $it" }
```

The first line in this example creates a variable `ch` which holds a channel object. This channel emits the values specified as a parameter in the `of` channel factory. Thus the second line will print the following:

```console
value: 1
value: 3
value: 5
value: 7
```

The `Channel.of` channel factory works in a similar manner to `Channel.from` (which is now [deprecated](https://www.nextflow.io/docs/latest/channel.html#of)), fixing some inconsistent behaviors of the latter and providing better handling when specifying a range of values. For example, the following works with a range from 1 to 23:

```groovy linenums="1"
Channel
    .of(1..23, 'X', 'Y')
    .view()
```

### `fromList()`

The `Channel.fromList` channel factory creates a channel emitting the elements provided by a list object specified as an argument:

```groovy linenums="1"
list = ['hello', 'world']

Channel
    .fromList(list)
    .view()
```

### `fromPath()`

The `fromPath` channel factory creates a queue channel emitting one or more files matching the specified glob pattern.

```groovy linenums="1"
Channel.fromPath('./data/meta/*.csv')
```

This example creates a channel and emits as many items as there are files with a `csv` extension in the `./data/meta` folder. Each element is a file object implementing the [Path](https://docs.oracle.com/javase/8/docs/api/java/nio/file/Paths.html) interface.

!!! tip

    Two asterisks, i.e. `**`, works like `*` but cross directory boundaries. This syntax is generally used for matching complete paths. Curly brackets specify a collection of sub-patterns.

| Name          | Description                                                                                                                                |
| ------------- | ------------------------------------------------------------------------------------------------------------------------------------------ |
| glob          | When `true` interprets characters `*`, `?`, `[]` and `{}` as glob wildcards, otherwise handles them as normal characters (default: `true`) |
| type          | Type of path returned, either `file`, `dir` or `any` (default: `file`)                                                                     |
| hidden        | When `true` includes hidden files in the resulting paths (default: `false`)                                                                |
| maxDepth      | Maximum number of directory levels to visit (default: `no limit`)                                                                          |
| followLinks   | When `true` symbolic links are followed during directory tree traversal, otherwise they are managed as files (default: `true`)             |
| relative      | When `true` return paths are relative to the top-most common directory (default: `false`)                                                  |
| checkIfExists | When `true` throws an exception when the specified path does not exist in the file system (default: `false`)                               |

Learn more about the glob patterns syntax at [this link](https://docs.oracle.com/javase/tutorial/essential/io/fileOps.html#glob).

!!! exercise

    Use the `Channel.fromPath` channel factory to create a channel emitting all files with the suffix `.fq` in the `data/ggal/` directory and any subdirectory, in addition to hidden files. Then print the file names.

    ??? solution

        ```groovy linenums="1"
        Channel.fromPath('./data/ggal/**.fq', hidden: true)
            .view()
        ```

### `fromFilePairs()`

The `fromFilePairs` channel factory creates a channel emitting the file pairs matching a glob pattern provided by the user. The matching files are emitted as tuples, in which the first element is the grouping key of the matching pair and the second element is the list of files (sorted in lexicographical order).

```groovy linenums="1"
Channel
    .fromFilePairs('./data/ggal/*_{1,2}.fq')
    .view()
```

It will produce an output similar to the following:

```groovy
[liver, [/workspace/gitpod/nf-training/data/ggal/liver_1.fq, /workspace/gitpod/nf-training/data/ggal/liver_2.fq]]
[gut, [/workspace/gitpod/nf-training/data/ggal/gut_1.fq, /workspace/gitpod/nf-training/data/ggal/gut_2.fq]]
[lung, [/workspace/gitpod/nf-training/data/ggal/lung_1.fq, /workspace/gitpod/nf-training/data/ggal/lung_2.fq]]
```

!!! warning

    The glob pattern _must_ contain at least a star wildcard character (`*`).

| Name          | Description                                                                                                                    |
| ------------- | ------------------------------------------------------------------------------------------------------------------------------ |
| type          | Type of paths returned, either `file`, `dir` or `any` (default: `file`)                                                        |
| hidden        | When `true` includes hidden files in the resulting paths (default: `false`)                                                    |
| maxDepth      | Maximum number of directory levels to visit (default: `no limit`)                                                              |
| followLinks   | When `true` symbolic links are followed during directory tree traversal, otherwise they are managed as files (default: `true`) |
| size          | Defines the number of files each emitted item is expected to hold (default: `2`). Set to `-1` for any                          |
| flat          | When `true` the matching files are produced as sole elements in the emitted tuples (default: `false`)                          |
| checkIfExists | When `true`, it throws an exception of the specified path that does not exist in the file system (default: `false`)            |

!!! exercise

    Use the `fromFilePairs` channel factory to create a channel emitting all pairs of fastq read in the `data/ggal/` directory and print them. Then use the `flat: true` option and compare the output with the previous execution.

    ??? solution

        Use the following, with or without `flat: true`:

        ```groovy linenums="1"
        Channel.fromFilePairs('./data/ggal/*_{1,2}.fq', flat: true)
            .view()
        ```

        Then check the square brackets around the file names, to see the difference with `flat`.

### `fromSRA()`

The `Channel.fromSRA` channel factory makes it possible to query the [NCBI SRA](https://www.ncbi.nlm.nih.gov/sra) archive and returns a channel emitting the FASTQ files matching the specified selection criteria.

The query can be project ID(s) or accession number(s) supported by the [NCBI ESearch API](https://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.ESearch).

!!! info

    This function now requires an API key you can only get by logging into your NCBI account.

??? example "Instructions for NCBI login and key acquisition"

    1. Go to: <https://www.ncbi.nlm.nih.gov/>
    2. Click the top right "Log in" button to sign into NCBI. Follow their instructions.
    3. Once into your account, click the button at the top right, usually your ID.
    4. Go to Account settings
    5. Scroll down to the API Key Management section.
    6. Click on "Create an API Key".
    7. The page will refresh and the key will be displayed where the button was. Copy your key.

For example, the following snippet will print the contents of an NCBI project ID:

```groovy linenums="1"
params.ncbi_api_key = '<Your API key here>'

Channel
    .fromSRA(['SRP073307'], apiKey: params.ncbi_api_key)
    .view()
```

!!! info ""

    :material-lightbulb: Replace `<Your API key here>` with your API key.

This should print:

```groovy
[SRR3383346, [/vol1/fastq/SRR338/006/SRR3383346/SRR3383346_1.fastq.gz, /vol1/fastq/SRR338/006/SRR3383346/SRR3383346_2.fastq.gz]]
[SRR3383347, [/vol1/fastq/SRR338/007/SRR3383347/SRR3383347_1.fastq.gz, /vol1/fastq/SRR338/007/SRR3383347/SRR3383347_2.fastq.gz]]
[SRR3383344, [/vol1/fastq/SRR338/004/SRR3383344/SRR3383344_1.fastq.gz, /vol1/fastq/SRR338/004/SRR3383344/SRR3383344_2.fastq.gz]]
[SRR3383345, [/vol1/fastq/SRR338/005/SRR3383345/SRR3383345_1.fastq.gz, /vol1/fastq/SRR338/005/SRR3383345/SRR3383345_2.fastq.gz]]
// (remaining omitted)
```

Multiple accession IDs can be specified using a list object:

```groovy linenums="1"
ids = ['ERR908507', 'ERR908506', 'ERR908505']
Channel
    .fromSRA(ids, apiKey: params.ncbi_api_key)
    .view()
```

```groovy
[ERR908507, [/vol1/fastq/ERR908/ERR908507/ERR908507_1.fastq.gz, /vol1/fastq/ERR908/ERR908507/ERR908507_2.fastq.gz]]
[ERR908506, [/vol1/fastq/ERR908/ERR908506/ERR908506_1.fastq.gz, /vol1/fastq/ERR908/ERR908506/ERR908506_2.fastq.gz]]
[ERR908505, [/vol1/fastq/ERR908/ERR908505/ERR908505_1.fastq.gz, /vol1/fastq/ERR908/ERR908505/ERR908505_2.fastq.gz]]
```

!!! info

    Read pairs are implicitly managed and are returned as a list of files.

It’s straightforward to use this channel as an input using the usual Nextflow syntax. The code below creates a channel containing two samples from a public SRA study and runs FASTQC on the resulting files. See:

```groovy linenums="1"
params.ncbi_api_key = '<Your API key here>'

params.accession = ['ERR908507', 'ERR908506']

process FASTQC {
    input:
    tuple val(sample_id), path(reads_file)

    output:
    path("fastqc_${sample_id}_logs")

    script:
    """
    mkdir fastqc_${sample_id}_logs
    fastqc -o fastqc_${sample_id}_logs -f fastq -q ${reads_file}
    """
}

workflow {
    reads = Channel.fromSRA(params.accession, apiKey: params.ncbi_api_key)
    FASTQC(reads)
}
```

If you want to run the pipeline above and do not have fastqc installed in your machine, don’t forget what you learned in the previous section. Run this pipeline with `-with-docker biocontainers/fastqc:v0.11.5`, for example.

### Text files

The `splitText` operator allows you to split multi-line strings or text file items, emitted by a source channel into chunks containing n lines, which will be emitted by the resulting channel. See:

```groovy linenums="1"
Channel
    .fromPath('data/meta/random.txt') // (1)!
    .splitText() // (2)!
    .view() // (3)!
```

1. Instructs Nextflow to make a channel from the path `data/meta/random.txt`
2. The `splitText` operator splits each item into chunks of one line by default.
3. View contents of the channel.

You can define the number of lines in each chunk by using the parameter `by`, as shown in the following example:

```groovy linenums="1"
Channel
    .fromPath('data/meta/random.txt')
    .splitText(by: 2)
    .subscribe {
        print it;
        print "--- end of the chunk ---\n"
    }
```

!!! info

    The `subscribe` operator permits execution of user defined functions each time a new value is emitted by the source channel.

An optional closure can be specified in order to transform the text chunks produced by the operator. The following example shows how to split text files into chunks of 10 lines and transform them into capital letters:

```groovy linenums="1"
Channel
    .fromPath('data/meta/random.txt')
    .splitText(by: 10) { it.toUpperCase() }
    .view()
```

You can also make counts for each line:

```groovy linenums="1"
count = 0

Channel
    .fromPath('data/meta/random.txt')
    .splitText()
    .view { "${count++}: ${it.toUpperCase().trim()}" }
```

Finally, you can also use the operator on plain files (outside of the channel context):

```groovy linenums="1"
def f = file('data/meta/random.txt')
def lines = f.splitText()
def count = 0
for (String row : lines) {
    log.info "${count++} ${row.toUpperCase()}"
}
```

### Comma separate values (.csv)

The `splitCsv` operator allows you to parse text items emitted by a channel, that are CSV formatted.

It then splits them into records or groups them as a list of records with a specified length.

In the simplest case, just apply the `splitCsv` operator to a channel emitting a CSV formatted text files or text entries. For example, to view only the first and fourth columns:

```groovy linenums="1"
Channel
    .fromPath("data/meta/patients_1.csv")
    .splitCsv()
    // row is a list object
    .view { row -> "${row[0]},${row[3]}" }
```

When the CSV begins with a header line defining the column names, you can specify the parameter `header: true` which allows you to reference each value by its column name, as shown in the following example:

```groovy linenums="1"
Channel
    .fromPath("data/meta/patients_1.csv")
    .splitCsv(header: true)
    // row is a list object
    .view { row -> "${row.patient_id},${row.num_samples}" }
```

Alternatively, you can provide custom header names by specifying a list of strings in the header parameter as shown below:

```groovy linenums="1"
Channel
    .fromPath("data/meta/patients_1.csv")
    .splitCsv(header: ['col1', 'col2', 'col3', 'col4', 'col5'])
    // row is a list object
    .view { row -> "${row.col1},${row.col4}" }
```

You can also process multiple CSV files at the same time:

```groovy linenums="1"
Channel
    .fromPath("data/meta/patients_*.csv") // <-- just use a pattern
    .splitCsv(header: true)
    .view { row -> "${row.patient_id}\t${row.num_samples}" }
```

!!! tip

    Notice that you can change the output format simply by adding a different delimiter.

Finally, you can also operate on CSV files outside the channel context:

```groovy linenums="1"
def f = file('data/meta/patients_1.csv')
def lines = f.splitCsv()
for (List row : lines) {
    log.info "${row[0]} -- ${row[2]}"
}
```

!!! exercise

    Try inputting fastq reads into the RNA-Seq workflow from earlier using `.splitCsv`.

    ??? solution

        Add a CSV text file containing the following, as an example input with the name "fastq.csv":

        ```csv
        gut,/workspace/gitpod/nf-training/data/ggal/gut_1.fq,/workspace/gitpod/nf-training/data/ggal/gut_2.fq
        ```

        Then replace the input channel for the reads in `script7.nf`. Changing the following lines:

        ```groovy linenums="1"
        Channel
            .fromFilePairs(params.reads, checkIfExists: true)
            .set { read_pairs_ch }
        ```

        To a splitCsv channel factory input:

        ```groovy linenums="1" hl_lines="2 3 4"
        Channel
            .fromPath("fastq.csv")
            .splitCsv()
            .view { row -> "${row[0]},${row[1]},${row[2]}" }
            .set { read_pairs_ch }
        ```

        Finally, change the cardinality of the processes that use the input data. For example, for the quantification process, change it from:

        ```groovy linenums="1"
        process QUANTIFICATION {
            tag "$sample_id"

            input:
            path salmon_index
            tuple val(sample_id), path(reads)

            output:
            path sample_id, emit: quant_ch

            script:
            """
            salmon quant --threads $task.cpus --libType=U -i $salmon_index -1 ${reads[0]} -2 ${reads[1]} -o $sample_id
            """
        }
        ```

        To:

        ```groovy linenums="1" hl_lines="6 13"
        process QUANTIFICATION {
            tag "$sample_id"

            input:
            path salmon_index
            tuple val(sample_id), path(reads1), path(reads2)

            output:
            path sample_id, emit: quant_ch

            script:
            """
            salmon quant --threads $task.cpus --libType=U -i $salmon_index -1 ${reads1} -2 ${reads2} -o $sample_id
            """
        }
        ```

        Repeat the above for the fastqc step.

        ```groovy linenums="1"  hl_lines="5 13"
        process FASTQC {
            tag "FASTQC on $sample_id"

            input:
            tuple val(sample_id), path(reads1), path(reads2)

            output:
            path "fastqc_${sample_id}_logs"

            script:
            """
            mkdir fastqc_${sample_id}_logs
            fastqc -o fastqc_${sample_id}_logs -f fastq -q ${reads1} ${reads2}
            """
        }
        ```

        Now the workflow should run from a CSV file.

### Tab separated values (.tsv)

Parsing TSV files works in a similar way, simply add the `sep: '\t'` option in the `splitCsv` context:

```groovy linenums="1"
Channel
    .fromPath("data/meta/regions.tsv", checkIfExists: true)
    // use `sep` option to parse TAB separated files
    .splitCsv(sep: '\t')
    .view()
```

!!! exercise

    Try using the tab separation technique on the file `data/meta/regions.tsv`, but print just the first column, and remove the header.


    ??? solution

        ```groovy linenums="1"
        Channel
            .fromPath("data/meta/regions.tsv", checkIfExists: true)
            // use `sep` option to parse TAB separated files
            .splitCsv(sep: '\t', header: true)
            // row is a list object
            .view { row -> "${row.patient_id}" }
        ```

## More complex file formats

### JSON

We can also easily parse the JSON file format using the following groovy schema:

```groovy linenums="1"
import groovy.json.JsonSlurper

def f = file('data/meta/regions.json')
def records = new JsonSlurper().parse(f)


for (def entry : records) {
    log.info "$entry.patient_id -- $entry.feature"
}
```

!!! warning

    When using an older JSON version, you may need to replace `parse(f)` with `parseText(f.text)`

### YAML

This can also be used as a way to parse YAML files:

```groovy linenums="1"
import org.yaml.snakeyaml.Yaml

def f = file('data/meta/regions.yml')
def records = new Yaml().load(f)


for (def entry : records) {
    log.info "$entry.patient_id -- $entry.feature"
}
```

### Storage of parsers into modules

The best way to store parser scripts is to keep them in a Nextflow module file.

See the following Nextflow script:

```groovy linenums="1"
include { parseJsonFile } from './modules/parsers.nf'

process FOO {
    input:
    tuple val(patient_id), val(feature)

    output:
    stdout

    script:
    """
    echo $patient_id has $feature as feature
    """
}

workflow {
    Channel.fromPath('data/meta/regions*.json')
        | flatMap { parseJsonFile(it) }
        | map { record -> [record.patient_id, record.feature] }
        | unique
        | FOO
        | view
}
```

For this script to work, a module file called `parsers.nf` needs to be created and stored in a modules folder in the current directory.

The `parsers.nf` file should contain the `parseJsonFile` function. For example:

```groovy linenums="1"
import groovy.json.JsonSlurper

def parseJsonFile(json_file) {
    def f = file('data/meta/regions.json')
    def records = new JsonSlurper().parse(f)
    return records
}
```

Nextflow will use this as a custom function within the workflow scope.

!!! tip

    You will learn more about module files later in the [Modularization section](../modules/) of this tutorial.
