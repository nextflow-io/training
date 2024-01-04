---
description: Basic Nextflow Training Workshop
---

# Operators

Operators are methods that allow you to connect channels, transform values emitted by a channel, or apply some user-provided rules.

There are seven main groups of operators are described in greater detail within the Nextflow Reference Documentation, linked below:

1. [Filtering operators](https://www.nextflow.io/docs/latest/operator.html#filtering-operators)
2. [Transforming operators](https://www.nextflow.io/docs/latest/operator.html#transforming-operators)
3. [Splitting operators](https://www.nextflow.io/docs/latest/operator.html#splitting-operators)
4. [Combining operators](https://www.nextflow.io/docs/latest/operator.html#combining-operators)
5. [Forking operators](https://www.nextflow.io/docs/latest/operator.html#forking-operators)
6. [Maths operators](https://www.nextflow.io/docs/latest/operator.html#maths-operators)
7. [Other operators](https://www.nextflow.io/docs/latest/operator.html#other-operators)

## Basic example

!!! info ""

    Click the :material-plus-circle: icons in the code for explanations.

```groovy linenums="1"
nums = Channel.of(1, 2, 3, 4) // (1)!
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

```groovy linenums="1"
Channel
    .of(1, 2, 3, 4)
    .map { it -> it * it }
    .view()
```

## Basic operators

Here we explore some of the most commonly used operators.

### `view()`

The `view` operator prints the items emitted by a channel to the console standard output, appending a _new line_ character to each item. For example:

```groovy linenums="1"
Channel
    .of('foo', 'bar', 'baz')
    .view()
```

```console title="Output"
foo
bar
baz
```

An optional _closure_ parameter can be specified to customize how items are printed. For example:

```groovy linenums="1"
Channel
    .of('foo', 'bar', 'baz')
    .view { "- $it" }
```

```console title="Output"
- foo
- bar
- baz
```

### `map()`

The `map` operator applies a function of your choosing to every item emitted by a channel and returns the items obtained as a new channel. The function applied is called the _mapping_ function and is expressed with a _closure_ as shown in the example below:

```groovy linenums="1"
Channel
    .of('hello', 'world')
    .map { it -> it.reverse() }
    .view()
```

A `map` can associate a generic _tuple_ to each element and can contain any data.

```groovy linenums="1"
Channel
    .of('hello', 'world')
    .map { word -> [word, word.size()] }
    .view { word, len -> "$word contains $len letters" }
```

!!! question "Exercise"

    Use `fromPath` to create a channel emitting the _fastq_ files matching the pattern `data/ggal/*.fq`, then use `map` to return a pair containing the file name and the path itself, and finally, use `view` to print the resulting channel.

    ??? solution

        ```groovy linenums="1"
        Channel
            .fromPath('data/ggal/*.fq')
            .map { file -> [file.name, file] }
            .view { name, file -> "> $name : $file" }
        ```

### `mix()`

The `mix` operator combines the items emitted by two (or more) channels into a single channel.

```groovy linenums="1"
my_channel_1 = Channel.of(1, 2, 3)
my_channel_2 = Channel.of('a', 'b')
my_channel_3 = Channel.of('z')

my_channel_1
    .mix(my_channel_2, my_channel_3)
    .view()
```

```console title="Output"
1
2
a
3
b
z
```

!!! warning

    The items in the resulting channel have the same order as in the respective original channels. However, there is no guarantee that the element of the second channel are appended after the elements of the first. Indeed, in the example above, the element `a` has been printed before `3`.

### `flatten()`

The `flatten` operator transforms a channel in such a way that every _tuple_ is flattened so that each entry is emitted as a sole element by the resulting channel.

```groovy linenums="1"
foo = [1, 2, 3]
bar = [4, 5, 6]

Channel
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

```groovy linenums="1"
Channel
    .of(1, 2, 3, 4)
    .collect()
    .view()
```

It prints a single value:

```console title="Output"
[1, 2, 3, 4]
```

!!! info

    The result of the `collect` operator is a **value** channel.

### `groupTuple()`

The `groupTuple` operator collects tuples (or lists) of values emitted by the source channel, grouping the elements that share the same key. Finally, it emits a new tuple object for each distinct key collected.

Try the following example:

```groovy linenums="1"
Channel
    .of([1, 'A'], [1, 'B'], [2, 'C'], [3, 'B'], [1, 'C'], [2, 'A'], [3, 'D'])
    .groupTuple()
    .view()
```

```console title="Output"
[1, [A, B, C]]
[2, [C, A]]
[3, [B, D]]
```

This operator is useful to process a group together with all the elements that share a common property or grouping key.

!!! question "Exercise"

    Use `fromPath` to create a channel emitting all of the files in the folder `data/meta/`, then use a `map` to associate the `baseName` prefix to each file. Finally, group all files that have the same common prefix.

    ??? solution

        ```groovy linenums="1"
        Channel
            .fromPath('data/meta/*')
            .map { file -> tuple(file.baseName, file) }
            .groupTuple()
            .view { baseName, file -> "> $baseName : $file" }
        ```

### `join()`

The `join` operator creates a channel that joins together the items emitted by two channels with a matching key. The key is defined, by default, as the first element in each item emitted.

```groovy linenums="1"
left = Channel.of(['X', 1], ['Y', 2], ['Z', 3], ['P', 7])
right = Channel.of(['Z', 6], ['Y', 5], ['X', 4])
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

The selection criterion is defined by specifying a closure that provides one or more boolean expressions, each of which is identified by a unique label. For the first expression that evaluates to a true value, the item is bound to a named channel as the label identifier. For example:

```groovy linenums="1"
Channel
    .of(1, 2, 3, 40, 50)
    .branch {
        small: it < 10
        large: it > 10
    }
    .set { result }

result.small.view { "$it is small" }
result.large.view { "$it is large" }
```

!!! info

    The `branch` operator returns a multi-channel object (i.e., a variable that holds more than one channel object).

!!! note

    In the above example, what would happen to a value of 10? To deal with this, you can also use `>=`.

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
    .view { row -> "${row[0]}, ${row[3]}" }
```

When the CSV begins with a header line defining the column names, you can specify the parameter `header: true` which allows you to reference each value by its column name, as shown in the following example:

```groovy linenums="1"
Channel
    .fromPath("data/meta/patients_1.csv")
    .splitCsv(header: true)
    // row is a list object
    .view { row -> "${row.patient_id}, ${row.num_samples}" }
```

Alternatively, you can provide custom header names by specifying a list of strings in the header parameter as shown below:

```groovy linenums="1"
Channel
    .fromPath("data/meta/patients_1.csv")
    .splitCsv(header: ['col1', 'col2', 'col3', 'col4', 'col5'])
    // row is a list object
    .view { row -> "${row.col1}, ${row.col4}" }
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

!!! question "Exercise"

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
            .view { row -> "${row[0]}, ${row[1]}, ${row[2]}" }
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

!!! question "Exercise"

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

We can also easily parse the JSON file format using the `splitJson` channel operator.

The `splitJson` operator supports JSON arrays:

=== "Source code"

    ```groovy linenums="1"
    Channel
        .of('["Sunday", "Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday"]')
        .splitJson()
        .view { "Item: ${it}" }
    ```

=== "Output"

    ```console
    Item: Sunday
    Item: Monday
    Item: Tuesday
    Item: Wednesday
    Item: Thursday
    Item: Friday
    Item: Saturday
    ```

JSON objects:

=== "Source code"

    ```groovy linenums="1"
    Channel
        .of('{"player": {"name": "Bob", "height": 180, "champion": false}}')
        .splitJson()
        .view { "Item: ${it}" }
    ```

=== "Output"

    ```console
    Item: [key:player, value:[name:Bob, height:180, champion:false]]
    ```

And even a JSON array of JSON objects!

=== "Source code"

    ```groovy linenums="1"
    Channel
        .of('[{"name": "Bob", "height": 180, "champion": false}, \
            {"name": "Alice", "height": 170, "champion": false}]')
        .splitJson()
        .view { "Item: ${it}" }
    ```

=== "Output"

    ```console
    Item: [name:Bob, height:180, champion:false]
    Item: [name:Alice, height:170, champion:false]
    ```

Files containing JSON content can also be parsed:

=== "Source code"

    ```groovy linenums="1"
    Channel
        .fromPath('file.json')
        .splitJson()
        .view { "Item: ${it}" }
    ```

=== "file.json"

    ```json
    [
      { "name": "Bob", "height": 180, "champion": false },
      { "name": "Alice", "height": 170, "champion": false }
    ]
    ```

=== "Output"

    ```console
    Item: [name:Bob, height:180, champion:false]
    Item: [name:Alice, height:170, champion:false]
    ```

### YAML

This can also be used as a way to parse YAML files:

=== "Source code"

    ```groovy linenums="1"
    import org.yaml.snakeyaml.Yaml

    def f = file('data/meta/regions.yml')
    def records = new Yaml().load(f)


    for (def entry : records) {
        log.info "$entry.patient_id -- $entry.feature"
    }
    ```

=== "data/meta/regions.yml"

    ```yaml
    --8<-- "nf-training/data/meta/regions.yml"
    ```

=== "Output"

    ```console
    ATX-TBL-001-GB-01-105 -- pass_vafqc_flag
    ATX-TBL-001-GB-01-105 -- pass_stripy_flag
    ATX-TBL-001-GB-01-105 -- pass_manual_flag
    ATX-TBL-001-GB-01-105 -- other_region_selection_flag
    ATX-TBL-001-GB-01-105 -- ace_information_gained
    ATX-TBL-001-GB-01-105 -- concordance_flag
    ATX-TBL-001-GB-01-105 -- pass_vafqc_flag
    ATX-TBL-001-GB-01-105 -- pass_stripy_flag
    ATX-TBL-001-GB-01-105 -- pass_manual_flag
    ATX-TBL-001-GB-01-105 -- other_region_selection_flag
    ATX-TBL-001-GB-01-105 -- ace_information_gained
    ATX-TBL-001-GB-01-105 -- concordance_flag
    ATX-TBL-001-GB-01-105 -- pass_vafqc_flag
    ATX-TBL-001-GB-01-105 -- pass_stripy_flag
    ```

### Storage of parsers into modules

The best way to store parser scripts is to keep them in a Nextflow module file.

Let's say we don't have a JSON channel operator, but we create a function instead. The `parsers.nf` file should contain the `parseJsonFile` function. See the contente below:

=== "Source code"

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
        Channel
            .fromPath('data/meta/regions*.json')
            | flatMap { parseJsonFile(it) }
            | map { record -> [record.patient_id, record.feature] }
            | unique
            | FOO
            | view
    }
    ```

=== "./modules/parsers.nf"

    ```groovy linenums="1"
    import groovy.json.JsonSlurper

    def parseJsonFile(json_file) {
        def f = file(json_file)
        def records = new JsonSlurper().parse(f)
        return records
    }
    ```

=== "Output"

    ```console
    ATX-TBL-001-GB-01-105 has pass_stripy_flag as feature

    ATX-TBL-001-GB-01-105 has ace_information_gained as feature

    ATX-TBL-001-GB-01-105 has concordance_flag as feature

    ATX-TBL-001-GB-01-105 has pass_vafqc_flag as feature

    ATX-TBL-001-GB-01-105 has pass_manual_flag as feature

    ATX-TBL-001-GB-01-105 has other_region_selection_flag as feature
    ```

Nextflow will use this as a custom function within the workflow scope.

!!! tip

    You will learn more about module files later in the [Modularization section](../modules/) of this tutorial.

## More resources

Check the [operators documentation](https://www.nextflow.io/docs/latest/operator.html) on Nextflow web site.
