---
title: Operators
description: Fundamentals Nextflow Training Workshop
---

# Operators

Nextflow operators are methods that allow you to manipulate channels. Every operator, with the exception of `set` and `subscribe`, produces one or more new channels, allowing you to chain operators to fit your needs.

There are seven main groups of operators are described in greater detail within the Nextflow Reference Documentation, linked below:

1. [Filtering operators](https://www.nextflow.io/docs/latest/operator.html#filtering-operators)
2. [Transforming operators](https://www.nextflow.io/docs/latest/operator.html#transforming-operators)
3. [Splitting operators](https://www.nextflow.io/docs/latest/operator.html#splitting-operators)
4. [Combining operators](https://www.nextflow.io/docs/latest/operator.html#combining-operators)
5. [Forking operators](https://www.nextflow.io/docs/latest/operator.html#forking-operators)
6. [Maths operators](https://www.nextflow.io/docs/latest/operator.html#maths-operators)
7. [Other operators](https://www.nextflow.io/docs/latest/operator.html#other-operators)

## Basic example

The `map` operator applies a function of your choosing to every item emitted by a channel, and returns the items so obtained as a new channel. The function applied is called the mapping function and is expressed with a closure as shown in the example below:

!!! info ""

    Click the :material-plus-circle: icons in the code for explanations.

```groovy linenums="1" title="snippet.nf"
nums = channel.of(1, 2, 3, 4) // (1)!
square = nums.map { it -> it * it } // (2)!
square.view() // (3)!
```

1. Creates a queue channel emitting four values
2. Creates a new channel, transforming each number into its square
3. Prints the channel content

<figure class="excalidraw">
--8<-- "docs/basic_training/img/channel-map.excalidraw.svg"
</figure>

Operators can also be chained to implement custom behaviors, so the previous snippet can also be written as:

```groovy linenums="1" title="snippet.nf"
channel
    .of(1, 2, 3, 4)
    .map { it -> it * it }
    .view()
```

!!! cboard-list-2 "Summary"

    In this step you have learned:

    1. The basic features of an operator

## Commonly used operators

Here you will explore some of the most commonly used operators.

### `view()`

The `view` operator prints the items emitted by a channel to the console standard output, appending a _new line_ character to each item. For example:

```groovy linenums="1" title="snippet.nf"
channel
    .of('foo', 'bar', 'baz')
    .view()
```

```console title="Output"
foo
bar
baz
```

An optional _closure_ parameter can be specified to customize how items are printed. For example:

```groovy linenums="1" title="snippet.nf"
channel
    .of('foo', 'bar', 'baz')
    .view { "- $it" }
```

```console title="Output"
- foo
- bar
- baz
```

### `map()`

The `map` operator applies a function of your choosing to every item emitted by a channel and returns the items obtained as a new channel. The function applied is called the _mapping_ function and is expressed with a _closure_. In the example below the groovy `reverse` method has been used to reverse the order of the characters in each string emitted by the channel.

```groovy linenums="1" title="snippet.nf"
channel
    .of('hello', 'world')
    .map { it -> it.reverse() }
    .view()
```

A `map` can associate a generic _tuple_ to each element and can contain any data. In the example below the groovy `size` method is used to return the length of each string emitted by the channel.

```groovy linenums="1" title="snippet.nf"
channel
    .of('hello', 'world')
    .map { word -> [word, word.size()] }
    .view()
```

```console title="Output"
[hello, 5]
[world, 5]
```

!!! question "Exercise"

    Use `fromPath` to create a channel emitting the _fastq_ files matching the pattern `data/ggal/*.fq`, then use `map` to return a pair containing the file name and the file path. Finally, use `view` to print the resulting channel.

    !!! hint

        You can use the `name` method to get the file name.

    ??? solution

        Here is one possible solution:

        ```groovy linenums="1" title="snippet.nf"
        channel
            .fromPath('data/ggal/*.fq')
            .map { file -> [file.name, file] }
            .view()
        ```

        Your output should look like this:

        ```console title="Output"
        [gut_1.fq, /workspaces/training/nf-training/data/ggal/gut_1.fq]
        [gut_2.fq, /workspaces/training/nf-training/data/ggal/gut_2.fq]
        [liver_1.fq, /workspaces/training/nf-training/data/ggal/liver_1.fq]
        [liver_2.fq, /workspaces/training/nf-training/data/ggal/liver_2.fq]
        [lung_1.fq, /workspaces/training/nf-training/data/ggal/lung_1.fq]
        [lung_2.fq, /workspaces/training/nf-training/data/ggal/lung_2.fq]
        ```

### `mix()`

The `mix` operator combines the items emitted by two (or more) channels.

```groovy linenums="1" title="snippet.nf"
my_channel_1 = channel.of(1, 2, 3)
my_channel_2 = channel.of('a', 'b')
my_channel_3 = channel.of('z')

my_channel_1
    .mix(my_channel_2, my_channel_3)
    .view()
```

It prints a single channel containing all the items emitted by the three channels:

```console title="Output"
1
2
a
3
b
z
```

!!! warning

    The items in the resulting channel have the same order as in the respective original channels. However, there is no guarantee that the elements of the second channel are appended after the elements of the first. Indeed, in the example above, the element `a` has been printed before `3`.

### `flatten()`

The `flatten` operator transforms a channel in such a way that every _tuple_ is flattened so that each entry is emitted as a sole element by the resulting channel.

```groovy linenums="1" title="snippet.nf"
foo = [1, 2, 3]
bar = [4, 5, 6]

channel
    .of(foo, bar)
    .flatten()
    .view()
```

```console title="Output"
1
2
3
4
5
6
```

### `collect()`

The `collect` operator collects all of the items emitted by a channel in a list and returns the object as a sole emission.

```groovy linenums="1" title="snippet.nf"
channel
    .of(1, 2, 3, 4)
    .collect()
    .view()
```

```console title="Output"
[1, 2, 3, 4]
```

!!! info

    The result of the `collect` operator is a **value** channel.

### `groupTuple()`

The `groupTuple` operator collects tuples (or lists) of values emitted by the source channel, grouping the elements that share the same key. Finally, it emits a new tuple object for each distinct key collected.

```groovy linenums="1" title="snippet.nf"
channel
    .of([1, 'A'], [1, 'B'], [2, 'C'], [3, 'B'], [1, 'C'], [2, 'A'], [3, 'D'])
    .groupTuple()
    .view()
```

```console title="Output"
[1, [A, B, C]]
[2, [C, A]]
[3, [B, D]]
```

This operator is especially useful to process a group together with all the elements that share a common property or grouping key.

!!! question "Exercise"

    Use `fromPath` to create a channel emitting all of the files in the folder `data/meta/`, then use a `map` to associate the `baseName` method to each file. Finally, group all files that have the same common prefix.

    ??? solution

        ```groovy linenums="1" title="snippet.nf"
        channel
            .fromPath('data/meta/*')
            .map { file -> tuple(file.baseName, file) }
            .groupTuple()
            .view()
        ```

        ```console title="Output"
        [patients_1, [/workspaces/training/nf-training/data/meta/patients_1.csv]]
        [patients_2, [/workspaces/training/nf-training/data/meta/patients_2.csv]]
        [random, [/workspaces/training/nf-training/data/meta/random.txt]]
        [regions, [/workspaces/training/nf-training/data/meta/regions.json, /workspaces/training/nf-training/data/meta/regions.tsv, /workspaces/training/nf-training/data/meta/regions.yml]]
        [regions2, [/workspaces/training/nf-training/data/meta/regions2.json]]
        ```

### `join()`

The `join` operator creates a channel that joins together the items emitted by two channels with a matching key. The key is defined, by default, as the first element in each item emitted.

```groovy linenums="1" title="snippet.nf"
left = channel.of(['X', 1], ['Y', 2], ['Z', 3], ['P', 7])
right = channel.of(['Z', 6], ['Y', 5], ['X', 4])
left.join(right).view()
```

```console title="Output"
[Z, 3, 6]
[Y, 2, 5]
[X, 1, 4]
```

!!! note

    Notice _P_ is missing in the final result.

### `branch()`

The `branch` operator allows you to forward the items emitted by a source channel to one or more output channels.

The selection criterion is defined by specifying a closure that provides one or more boolean expressions, each of which is identified by a unique label. For the first expression that evaluates to a true value, the item is bound to a named channel as the label identifier.

```groovy linenums="1" title="snippet.nf"
channel
    .of(1, 2, 3, 40, 50)
    .branch {
        small: it < 10
        large: it > 10
    }
    .set { result }

result.small.view { "$it is small" }
result.large.view { "$it is large" }
```

```console title="Output"
1 is small
40 is large
2 is small
3 is small
50 is large
```

!!! info

    The `branch` operator returns a multi-channel object (i.e., a variable that holds more than one channel object).

!!! note

    In the above example, what would happen to a value of 10? To deal with this, you can also use `>=`.

!!! cboard-list-2 "Summary"

    In this step you have learned:

    1. How to use the `view` operator to print the content of a channel
    2. How to use the `map` operator to transform the content of a channel
    3. How to use the `mix` operator to combine the content of two or more channels
    4. How to use the `flatten` operator to flatten the content of a channel
    5. How to use the `collect` operator to collect the content of a channel
    6. How to use the `groupTuple` operator to group the content of a channel
    7. How to use the `join` operator to join the content of two channels
    8. How to use the `branch` operator to split the content of a channel

## Text files

### `splitText()`

The `splitText` operator allows you to split multi-line strings or text file items, emitted by a source channel into chunks containing _n_ lines, which will be emitted by the resulting channel.

```groovy linenums="1" title="snippet.nf"
channel
    .fromPath('data/meta/random.txt') // (1)!
    .splitText() // (2)!
    .view() // (3)!
```

1. Instructs Nextflow to make a channel from the path `data/meta/random.txt`
2. The `splitText` operator splits each item into chunks of one line by default.
3. View contents of the channel.

```console title="Output"
Lorem Ipsum is simply dummy text of the printing and typesetting industry.
Lorem Ipsum has been the industry's standard dummy text ever since the 1500s,
when an unknown printer took a galley of type and scrambled it to make a type specimen book.
It has survived not only five centuries, but also the leap into electronic typesetting,
...
```

You can define the number of lines in each chunk by using the parameter `by`, as shown in the following example:

```groovy linenums="1" title="snippet.nf"
channel
    .fromPath('data/meta/random.txt')
    .splitText(by: 2)
    .view()
```

```console title="Output"
Lorem Ipsum is simply dummy text of the printing and typesetting industry.
Lorem Ipsum has been the industry's standard dummy text ever since the 1500s,

when an unknown printer took a galley of type and scrambled it to make a type specimen book.
It has survived not only five centuries, but also the leap into electronic typesetting,
...
```

An optional closure can also be specified in order to transform the text chunks produced by the operator. The following example shows how to split text files into chunks of 2 lines and transform them into capital letters:

```groovy linenums="1" title="snippet.nf"
channel
    .fromPath('data/meta/random.txt')
    .splitText(by: 2) { it.toUpperCase() }
    .view()
```

```console title="Output"
LOREM IPSUM IS SIMPLY DUMMY TEXT OF THE PRINTING AND TYPESETTING INDUSTRY.
LOREM IPSUM HAS BEEN THE INDUSTRY'S STANDARD DUMMY TEXT EVER SINCE THE 1500S,

WHEN AN UNKNOWN PRINTER TOOK A GALLEY OF TYPE AND SCRAMBLED IT TO MAKE A TYPE SPECIMEN BOOK.
IT HAS SURVIVED NOT ONLY FIVE CENTURIES, BUT ALSO THE LEAP INTO ELECTRONIC TYPESETTING,
...
```

### `splitCsv()`

The `splitCsv` operator allows you to parse text items emitted by a channel, that are CSV formatted.

It then splits them into records or groups them as a list of records with a specified length.

In the simplest case, just apply the `splitCsv` operator to a channel emitting a CSV formatted text file or text entries. For example, to view only the first and fourth columns:

```groovy linenums="1" title="snippet.nf"
channel
    .fromPath("data/meta/patients_1.csv")
    .splitCsv()
    .view { row -> "${row[0]}, ${row[3]}" }
```

```console title="Output"
patient_id, num_samples
ATX-TBL-001-GB-02-117, 3
ATX-TBL-001-GB-01-110, 3
ATX-TBL-001-GB-03-101, 3
ATX-TBL-001-GB-04-201, 3
ATX-TBL-001-GB-02-120, 3
ATX-TBL-001-GB-04-102, 3
ATX-TBL-001-GB-03-104, 3
ATX-TBL-001-GB-03-103, 3
```

When the CSV begins with a header line defining the column names, you can specify the parameter `header: true` which allows you to reference each value by its column name, as shown in the following example:

```groovy linenums="1" title="snippet.nf"
channel
    .fromPath("data/meta/patients_1.csv")
    .splitCsv(header: true)
    // row is a list object
    .view { row -> "${row.patient_id}, ${row.num_samples}" }
```

Alternatively, you can provide custom header names by specifying a list of strings in the header parameter as shown below:

```groovy linenums="1" title="snippet.nf"
channel
    .fromPath("data/meta/patients_1.csv")
    .splitCsv(header: ['col1', 'col2', 'col3', 'col4', 'col5'])
    .view { row -> "${row.col1}, ${row.col4}" }
```

```console title="Output"
patient_id, num_samples
ATX-TBL-001-GB-02-117, 3
ATX-TBL-001-GB-01-110, 3
ATX-TBL-001-GB-03-101, 3
ATX-TBL-001-GB-04-201, 3
ATX-TBL-001-GB-02-120, 3
ATX-TBL-001-GB-04-102, 3
ATX-TBL-001-GB-03-104, 3
ATX-TBL-001-GB-03-103, 3
```

You can also process multiple CSV files at the same time:

```groovy linenums="1" title="snippet.nf"
channel
    .fromPath("data/meta/patients_*.csv") // <-- just use a pattern
    .splitCsv(header: true)
    .view { row -> "${row.patient_id}\t${row.num_samples}" }
```

```console title="Output"
ATX-TBL-001-GB-02-117   3
ATX-TBL-001-GB-01-110   3
ATX-TBL-001-GB-03-101   3
ATX-TBL-001-GB-04-201   3
ATX-TBL-001-GB-02-120   3
ATX-TBL-001-GB-04-102   3
ATX-TBL-001-GB-03-104   3
ATX-TBL-001-GB-03-103   3
ATX-TBL-001-GB-01-111   2
ATX-TBL-001-GB-01-112   3
ATX-TBL-001-GB-04-202   3
ATX-TBL-001-GB-02-124   3
ATX-TBL-001-GB-02-107   3
ATX-TBL-001-GB-01-105   3
ATX-TBL-001-GB-02-108   3
ATX-TBL-001-GB-01-113   3
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

!!! question "Exercise"

    Create a CSV file and use it as input for `script7.nf`, part of the [Simple RNA-Seq workflow tutorial](https://training.nextflow.io/basic_training/rnaseq_pipeline/).

    ??? solution

        Add a CSV text file containing the following, as an example input with the name "fastq.csv":

        ```csv title="fastq.csv"
        gut,/workspaces/training/nf-training/data/ggal/gut_1.fq,/workspaces/training/nf-training/data/ggal/gut_2.fq
        ```

        Then replace the input channel for the reads in `script7.nf`. Changing the following lines:

        ```groovy linenums="1"
        channel
            .fromFilePairs(params.reads, checkIfExists: true)
            .set { read_pairs_ch }
        ```

        To a splitCsv channel factory input:

        ```groovy linenums="1" title="script7.nf"
        channel
            .fromPath("fastq.csv")
            .splitCsv()
            .view { row -> "${row[0]}, ${row[1]}, ${row[2]}" }
            .set { read_pairs_ch }
        ```

        Finally, change the cardinality of the processes that use the input data:

        ```groovy linenums="1" title="script7.nf"
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

        ```groovy linenums="1" hl_lines="5 13" title="script7.nf"
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

Parsing TSV files works in a similar way. Simply add the `sep: '\t'` option in the `splitCsv` context:

```groovy linenums="1" title="snippet.nf"
channel
    .fromPath("data/meta/regions.tsv", checkIfExists: true)
    // use `sep` option to parse TAB separated files
    .splitCsv(sep: '\t')
    .view()
```

!!! question "Exercise"

    Use the tab separation technique on the file `data/meta/regions.tsv`, but print just the first column, and remove the header.

    ??? solution

        ```groovy linenums="1" title="snippet.nf"
        channel
            .fromPath("data/meta/regions.tsv", checkIfExists: true)
            // use `sep` option to parse TAB separated files
            .splitCsv(sep: '\t', header: true)
            // row is a list object
            .view { row -> "${row.patient_id}" }
        ```

### `splitJson()`

You can parse the JSON file format using the `splitJson` channel operator.

The `splitJson` operator supports JSON arrays:

```groovy linenums="1" title="snippet.nf"
channel
    .of('["Sunday", "Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday"]')
    .splitJson()
    .view()
```

```console title="Output"
Sunday
Monday
Tuesday
Wednesday
Thursday
Friday
Saturday
```

As well as JSON arrays in objects:

```groovy linenums="1" title="snippet.nf"
channel
    .of('{"player": {"name": "Bob", "height": 180, "champion": false}}')
    .splitJson()
    .view()
```

```console title="Output"
[value:[name:Bob, height:180, champion:false], key:player]
```

And even a JSON array of JSON objects:

```groovy linenums="1" title="snippet.nf"
channel
    .of('[{"name": "Bob", "height": 180, "champion": false}, \
          {"name": "Alice", "height": 170, "champion": false}]')
    .splitJson()
    .view()
```

```console title="Output"
[name:Bob, height:180, champion:false]
[name:Alice, height:170, champion:false]
```

You can also parse JSON files directly:

```json title="file.json"
[
  { "name": "Bob", "height": 180, "champion": false },
  { "name": "Alice", "height": 170, "champion": false }
]
```

```groovy linenums="1" title="snippet.nf"
channel
    .fromPath('file.json')
    .splitJson()
    .view()
```

```console title="Output"
[name:Bob, height:180, champion:false]
[name:Alice, height:170, champion:false]
```

!!! cboard-list-2 "Summary"

    In this step you have learned:

    1. How to use the `splitText` operator to split text files of various formats
    2. How to use the `splitJson` operator to split JSON files of various formats

## More resources

Check the [operators documentation](https://www.nextflow.io/docs/latest/operator.html) on Nextflow web site.
