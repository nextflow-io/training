# Operator Tour

In this chapter, we take a curated tour of the Nextflow operators. Commonly used and well understood operators are not covered here - only those that we've seen could use more attention or those where the usage could be more elaborate. These set of operators have been chosen to illustrate tangential concepts and Nextflow features.

## `map`

### Basics

Map is certainly the most commonly used of the operators covered here. It's a way to supply a closure through which each element in the channel is passed. The return value of the closure is emitted as an element in a new output channel. A canonical example is a closure that multiplies two numbers:

```groovy linenums="1"
workflow {
    channel.of( 1, 2, 3, 4, 5 )
        .map { num -> num * num }
        .view()
}
```

The code above is available in a starter `main.nf` file available at `advanced/operators/main.nf`. It is recommended to open and edit this file to follow along with the examples given in the rest of this chapter. The workflow can be executed with:

```bash
cd operators
nextflow run .
```

The code creates a workflow that emits the numbers 1 through 5 into a channel, applies the `map` operator to square each number, and then prints the results. This demonstrates how `map` transforms each element in a channel by applying a closure, producing a new channel with the squared values.

The value passed into `map` is assigned to a variable using the `->` syntax (e.g., `num` in `{ num -> ... }`). This variable is local to the closure and used for any transformations inside it.

Groovy is an optionally typed language, and it is possible to specify the type of the argument passed to the closure.

```groovy linenums="1" hl_lines="3"
workflow {
    channel.of( 1, 2, 3, 4, 5 )
        .map { Integer num -> num * num }
        .view()
}
```

### Named Closures

If you find yourself re-using the same closure multiple times in your pipeline, the closure can be named and referenced:

```groovy linenums="1" hl_lines="2 6"
workflow {
    def squareIt = { Integer num -> num * num }

    channel.of( 1, 2, 3, 4, 5 )
        .map( squareIt )
        .view()
}
```

If you have these re-usable closures defined, you can compose them together.

```groovy linenums="1" hl_lines="2-3 6"
workflow {
    def squareIt = { num -> num * num }
    def addTwo = { num -> num + 2 }

    channel.of( 1, 2, 3, 4, 5 )
        .map( squareIt >> addTwo )
        .view()
}
```

```console title="Output"
N E X T F L O W  ~  version 23.04.1
Launching `./main.nf` [focused_borg] DSL2 - revision: f3c3e751fe
3
6
11
18
27
```

The above is the same as writing:

```groovy linenums="1" hl_lines="6-7"
workflow {
    def squareIt = { num -> num * num }
    def addTwo = { num -> num + 2 }

    channel.of( 1, 2, 3, 4, 5 )
        .map( squareIt )
        .map( addTwo )
        .view()
}
```

For those inclined towards functional programming, you'll be happy to know that closures can be curried:

```groovy linenums="1" hl_lines="2-3 6"
workflow {
    def timesN = { multiplier, num -> num * multiplier }
    def timesTen = timesN.curry(10)

    channel.of( 1, 2, 3, 4, 5 )
        .map( timesTen )
        .view()
}
```

## `view`

In addition to the argument-less usage of `view` as shown above, this operator can also take a closure to customize the stdout message. We can create a closure to print the value of the elements in a channel as well as their type, for example:

```groovy linenums="1" hl_lines="7"
workflow {
    def timesN = { multiplier, num -> num * multiplier }
    def timesTen = timesN.curry(10)

    channel.of( 1, 2, 3, 4, 5 )
        .map( timesTen )
        .view { value -> "Found '$value' (${value.getClass()})"}
}
```

!!! note "Most closures will remain anonymous"

    In many cases, it is simply cleaner to keep the closure anonymous, defined inline. Giving closures a name is only recommended when you find yourself defining the same or similar closures repeatedly in a given workflow.

## `splitCsv`

A common Nextflow pattern is for a simple samplesheet to be passed as primary input into a workflow. We'll see some more complicated ways to manage these inputs later on in the workshop, but the `splitCsv` ([docs](https://www.nextflow.io/docs/latest/operator.html#splitcsv)) is an excellent tool to have in a pinch. This operator will parse a csv/tsv and return a channel where each item is a row in the csv/tsv:

```groovy linenums="1" hl_lines="2"
workflow {
    channel.fromPath("data/samplesheet.csv")
        .splitCsv( header: true )
        .view()
}
```

!!! exercise

    From the directory `advanced/operators`, use the `splitCsv` and `map` operators to read the file `data/samplesheet.csv` and return a channel that would be suitable input to the process below. Feel free to consult the [splitCsv documentation](https://www.nextflow.io/docs/latest/operator.html#splitcsv) for tips.

    ```groovy linenums="1"
    process FastQC {
        input:
        tuple val(id), path(fastqs)
        // ... rest of the process
    ```

    ??? solution
        Specifying the `header` argument in the `splitCsv` operator, we have convenient named access to csv elements. The closure returns a list of two elements where the second element a list of paths.

        ```groovy linenums="1" hl_lines="4-6"
        workflow {
            channel.fromPath("data/samplesheet.csv")
                .splitCsv( header: true )
                .map { row ->
                    [row.id, [file(row.fastq1), file(row.fastq2)]]
                }
                .view()
        }
        ```

        !!! warning "Convert Strings to Paths"

            The fastq paths are simple strings in the context of a csv row. In order to pass them as paths to a Nextflow process, they need to be converted into objects that adjere to the `Path` interface. This is accomplished by wrapping them in `file`.

        In the sample above, we've lost an important piece of metadata - the tumor/normal classification, choosing only the sample id as the first element in the output list.

        In the next chapter, we'll discuss the "meta map" pattern in more detail, but we can preview that here.

        ```groovy linenums="1" hl_lines="4-7"
        workflow {
            channel.fromPath("data/samplesheet.csv")
                .splitCsv( header: true )
                .map { row ->
                    def metaMap = [id: row.id, type: row.type, repeat: row.repeat]
                    [metaMap, [file(row.fastq1), file(row.fastq2)]]
                }
                .view()
        }
        ```

        The construction of this map is very repetitive, and in the next chapter, we'll discuss some Groovy methods available on the `Map` class that can make this pattern more concise and less error-prone.

## `multiMap`

The `multiMap` ([documentation](https://www.nextflow.io/docs/latest/operator.html#multimap)) operator is a way of taking a single input channel and emitting into **multiple channels for each input element**.

Let's assume we've been given a samplesheet that has tumor/normal pairs bundled together on the same row. View the example samplesheet with:

```bash
cd operators
cat data/samplesheet.ugly.csv
```

Using the `splitCsv` operator would give us one entry that would contain all four fastq files. Let's consider that we wanted to split these fastqs into separate channels for tumor and normal. In other words, for every row in the samplesheet, we would like to emit an entry into two new channels. To do this, we can use the `multiMap` operator:

```groovy linenums="1" hl_lines="4-11"
workflow {
    channel.fromPath("data/samplesheet.ugly.csv")
        .splitCsv( header: true )
        .multiMap { row ->
            tumor:
                def tumor_meta = [id: row.id, type:'tumor', repeat:row.repeat]
                [tumor_meta, file(row.tumor_fastq_1), file(row.tumor_fastq_2)]
            normal:
                def normal_meta = [id: row.id, type:'normal', repeat:row.repeat]
                [normal_meta, file(row.normal_fastq_1), file(row.normal_fastq_2)]
        }
        .set { samples }

    samples.tumor.view { sample -> "Tumor: $sample"}
    samples.normal.view { sample -> "Normal: $sample"}
}
```

!!! tip "multiMapCriteria"

    The closure supplied to `multiMap` needs to return multiple channels, so using named closures as described in the `map` section above will not work. Fortunately, Nextflow provides the convenience `multiMapCriteria` method to allow you to define named `multiMap` closures should you need them. See the [`multiMap` documentation](https://www.nextflow.io/docs/latest/operator.html#multimap) for more info.

## `branch`

The `branch` operator ([documentation](https://www.nextflow.io/docs/latest/operator.html#branch)) is a way of taking a single input channel and emitting a new element into one (and only one) of a selection of output channels.

In the example above, the `multiMap` operator was necessary because we were supplied with a samplesheet that combined two pairs of fastq per row and we wanted to turn each row into new elements in multiple channels. If we were to use the neater samplesheet that had tumor/normal pairs on separate rows, we could use the `branch` operator to achieve the same result as we are routing each input element into a single output channel.

```groovy linenums="1" hl_lines="5-8"
workflow {
    channel.fromPath("data/samplesheet.csv")
        .splitCsv( header: true )
        .map { row -> [[id: row.id, repeat: row.repeat, type: row.type], [file(row.fastq1), file(row.fastq2)]] }
        .branch { meta, _reads ->
            tumor: meta.type == "tumor"
            normal: meta.type == "normal"
        }
        .set { samples }

    samples.tumor.view { sample -> "Tumor: $sample"}
    samples.normal.view { sample -> "Normal: $sample"}
}
```

An element is only emitted to the first channel were the test condition is met. If an element does not meet any of the tests, it is not emitted to any of the output channels. You can 'catch' any such samples by specifying `true` as a condition. If we knew that all samples would be either tumor or normal and no third 'type', we could write

```groovy linenums="5" hl_lines="2-3"
.branch { meta, reads ->
    tumor: meta.type == "tumor"
    normal: true
}
```

We may want to emit a slightly different element than the one passed as input. The `branch` operator can (optionally) return a _new_ element to a channel. For example, to add an extra key in the meta map of the tumor samples, we add a new line under the condition and return our new element. In this example, we modify the first element of the `List` to be a new list that is the result of merging the existing meta map with a new map containing a single key:

```groovy linenums="5" hl_lines="2-4"
.branch { meta, reads ->
    tumor: meta.type == "tumor"
        return [meta + [newKey: 'myValue'], reads]
    normal: true
}
```

!!! exercise

    How would you modify the element returned in the `tumor` channel to have the key:value pair `type:'abnormal'` instead of `type:'tumor'`?

    ??? solution

        There are many ways to accomplish this, but the map merging pattern introduced above can also be used to safely and concisely rename values in a map.

        ```groovy linenums="5"
        .branch { meta, reads ->
            tumor: meta.type == "tumor"
                return [meta + [type: 'abnormal'], reads]
            normal: true
        }
        ```

        !!! note "Merging maps is safe"

            Using the `+` operator to merge two or more Maps returns a _new_ Map. There are rare edge cases where modification of map rather than returning a new map can affect other channels. We discuss this further in the next chapter, but just be aware that this `+` operator is safer and often more convenient than modifying the `meta` object directly.

            See the Groovy [Map documentation](https://docs.groovy-lang.org/latest/html/groovy-jdk/java/util/Map.html#plus(java.util.Map)) for details.

### Multi-channel Objects

Certain Nextflow operators, such as `multiMap` and `branch`, return special objects containing multiple named channels. These multi-channel objects allow you to split data into separate streams while maintaining the ability to reference each stream independently.

```groovy linenums="1" hl_lines="3-6"
workflow {
    numbers = channel.of( 1, 2, 3, 4, 5 )
        .multiMap { num ->
            small: num
            large: num * 10
        }
    numbers.small.view { num -> "Small: $num"}
    numbers.large.view { num -> "Large: $num"}
}
```

This creates two channels accessible through the `numbers` object: `small` containing the original values (1, 2, 3, 4, 5) and `large` containing the values multiplied by 10 (10, 20, 30, 40, 50).

When a process requires multiple input channels, Nextflow automatically synchronizes the values from these channels. Each channel provides values as separate input parameters, and Nextflow pairs them together for each process execution:

```groovy linenums="1"
process MultiInput {
    debug true
    input:
    val(smallNum)
    val(bigNum)

    script:
    "echo -n small is $smallNum and big is $bigNum"
}

workflow {
```

The following will be kept synchronous, allowing you to supply multiple channel inputs to a process and keeping the order. This is the only place where order is guaranteed in Nextflow!

```groovy linenums="11" hl_lines="8"
workflow {
    numbers = channel.of( 1, 2, 3, 4, 5 )
        .multiMap { num ->
            small: num
            large: num * 10
        }

    MultiInput(numbers.small, numbers.large)
}
```

This workflow produces synchronized output:

```console title="Multi-channel Input Output"
small is 1 and big is 10
small is 2 and big is 20
small is 3 and big is 30
small is 4 and big is 40
small is 5 and big is 50
```

Multi-channel objects, created with `multiMap` or `branch`, let you split data streams while keeping named references. They allow related data to be processed in sync—Nextflow guarantees ordering only when multiple channels are used as process inputs.

## `groupTuple`

A common operation is to group elements from a _single_ channel where those elements share a common key. Take this example samplesheet as an example:

```groovy linenums="1" hl_lines="6"
workflow {
    channel.fromPath("data/samplesheet.csv")
        .splitCsv(header: true)
        .map { row ->
            def meta = [id: row.id, type: row.type]
            [meta, row.repeat, [row.fastq1, row.fastq2]]
        }
        .view()
}
```

We see that there are multiple rows where the first element in the item emitted by the channel is the Map `[id:sampleA, type:normal]` and items in the channel where the first element is the Map `[id:sampleA, type:tumor]`.

The `groupTuple` operator allows us to combine elements that share a common key:

```groovy linenums="1" hl_lines="8"
workflow {
    channel.fromPath("data/samplesheet.csv")
        .splitCsv(header: true)
        .map { row ->
            def meta = [id: row.id, type: row.type]
            [meta, row.repeat, [row.fastq1, row.fastq2]]
        }
        .groupTuple()
        .view()
}
```

## `transpose`

The transpose operator is often misunderstood. It can be thought of as the inverse of the `groupTuple` operator. Give the following workflow, the `groupTuple` and `transpose` operators cancel each other out. Removing lines 8 and 9 returns the same result.

Given a workflow that returns one element per sample, where we have grouped the samplesheet lines on a meta containing only id and type:

```groovy linenums="1"
workflow {
    channel.fromPath("data/samplesheet.csv")
        .splitCsv(header: true)
        .map { row ->
            def meta = [id: row.id, type: row.type]
            [meta, row.repeat, [row.fastq1, row.fastq2]]
        }
        .groupTuple()
        .view()
}
```

```console title="Output"
N E X T F L O W  ~  version 23.04.1
Launching `./main.nf` [spontaneous_rutherford] DSL2 - revision: 7dc1cc0039
[[id:sampleA, type:normal], [1, 2], [[data/reads/sampleA_rep1_normal_R1.fastq.gz, data/reads/sampleA_rep1_normal_R2.fastq.gz], [data/reads/sampleA_rep2_normal_R1.fastq.gz, data/reads/sampleA_rep2_normal_R2.fastq.gz]]]
[[id:sampleA, type:tumor], [1, 2], [[data/reads/sampleA_rep1_tumor_R1.fastq.gz, data/reads/sampleA_rep1_tumor_R2.fastq.gz], [data/reads/sampleA_rep2_tumor_R1.fastq.gz, data/reads/sampleA_rep2_tumor_R2.fastq.gz]]]
[[id:sampleB, type:normal], [1], [[data/reads/sampleB_rep1_normal_R1.fastq.gz, data/reads/sampleB_rep1_normal_R2.fastq.gz]]]
[[id:sampleB, type:tumor], [1], [[data/reads/sampleB_rep1_tumor_R1.fastq.gz, data/reads/sampleB_rep1_tumor_R2.fastq.gz]]]
[[id:sampleC, type:normal], [1], [[data/reads/sampleC_rep1_normal_R1.fastq.gz, data/reads/sampleC_rep1_normal_R2.fastq.gz]]]
[[id:sampleC, type:tumor], [1], [[data/reads/sampleC_rep1_tumor_R1.fastq.gz, data/reads/sampleC_rep1_tumor_R2.fastq.gz]]]
```

If we add in a `transpose`, each repeat number is matched back to the appropriate list of reads:

```groovy linenums="1" hl_lines="9"
workflow {
    channel.fromPath("data/samplesheet.csv")
        .splitCsv(header: true)
        .map { row ->
            def meta = [id: row.id, type: row.type]
            [meta, row.repeat, [row.fastq1, row.fastq2]]
        }
        .groupTuple()
        .transpose()
        .view()
}
```

```console title="Output"
N E X T F L O W  ~  version 23.04.1
Launching `./main.nf` [elegant_rutherford] DSL2 - revision: 2c5476b133
[[id:sampleA, type:normal], 1, [data/reads/sampleA_rep1_normal_R1.fastq.gz, data/reads/sampleA_rep1_normal_R2.fastq.gz]]
[[id:sampleA, type:normal], 2, [data/reads/sampleA_rep2_normal_R1.fastq.gz, data/reads/sampleA_rep2_normal_R2.fastq.gz]]
[[id:sampleA, type:tumor], 1, [data/reads/sampleA_rep1_tumor_R1.fastq.gz, data/reads/sampleA_rep1_tumor_R2.fastq.gz]]
[[id:sampleA, type:tumor], 2, [data/reads/sampleA_rep2_tumor_R1.fastq.gz, data/reads/sampleA_rep2_tumor_R2.fastq.gz]]
[[id:sampleB, type:normal], 1, [data/reads/sampleB_rep1_normal_R1.fastq.gz, data/reads/sampleB_rep1_normal_R2.fastq.gz]]
[[id:sampleB, type:tumor], 1, [data/reads/sampleB_rep1_tumor_R1.fastq.gz, data/reads/sampleB_rep1_tumor_R2.fastq.gz]]
[[id:sampleC, type:normal], 1, [data/reads/sampleC_rep1_normal_R1.fastq.gz, data/reads/sampleC_rep1_normal_R2.fastq.gz]]
[[id:sampleC, type:tumor], 1, [data/reads/sampleC_rep1_tumor_R1.fastq.gz, data/reads/sampleC_rep1_tumor_R2.fastq.gz]]
```

## `flatMap`

As the name suggests, the `flatMap` operator allows you to modify the elements in a channel and then flatten the resulting collection. This is useful if you need to "expand" elements in a channel an incoming element can turn into zero or more elements in the output channel. For example:

```groovy linenums="1" hl_lines="5"
workflow {
    numbers = channel.of(1, 2)

    numbers
        .flatMap { n -> [ n, n*10, n*100 ] }
        .view()
}
```

The input channel has two elements. For each element in the input channel, we return a List of length three. The List is flattened and each element in our returned list is emitted independently into the output channel:

```output
1
10
100
2
20
200
```

!!! exercise

    The `flatten` operation only "unfolds" one layer from the returned collection. Given this information, what do you expect the following workflow to return?

    ```
    workflow {
        numbers = channel.of(1, 2)

        numbers
            .flatMap { n -> [ n, [n*10, n*100] ] }
            .view()
    }
    ```

!!! exercise

    Let's say we have some collection of data for two samples:

    ```bash
    mkdir -p data/datfiles/sample{1,2}
    touch data/datfiles/sample1/data.{1,2,3,4,5,6,7}.dat
    touch data/datfiles/sample2/data.{1,2,3,4}.dat
    tree data/datfiles
    ```

    You would like to process these datfiles in batches - up to three at a time. The catch is that your batch process requires that the samples are processed independently. You have begun your workflow and grouped the samples together:

    ```groovy linenums="1" hl_lines="3-4"
    workflow {
        channel.fromPath("data/datfiles/sample*/*.dat", checkIfExists: true)
            .map { myfile -> [myfile.getParent().name, myfile] }
            .groupTuple()
            .view()
    }
    ```

    Which returns a channel with two elements corresponding to each sample:

    ```
    [sample2, [sample2/data.1.dat, sample2/data.3.dat, sample2/data.2.dat, sample2/data.4.dat]]
    [sample1, [sample1/data.1.dat, sample1/data.3.dat, sample1/data.2.dat, sample1/data.6.dat, sample1/data.5.dat, sample1/data.4.dat]]
    ```

    You would like to turn this channel into something that has at most three datfiles in each element. Something like:

    ```
    [sample1, [sample1/data.1.dat, sample1/data.3.dat, sample1/data.2.dat]]
    [sample1, [sample1/data.6.dat, sample1/data.5.dat, sample1/data.4.dat]]
    [sample2, [sample2/data.1.dat, sample2/data.3.dat, sample2/data.2.dat]]
    [sample2, [sample2/data.4.dat]]
    ```

    This is a challenging problem that pushes slightly beyond what we have covered.

    !!! tip
        You may find it helpful to use Groovy's `collate()` method ([docs](http://docs.groovy-lang.org/docs/groovy-2.3.5/html/groovy-jdk/java/util/List.html#collate(int)), [tutorial](https://blog.mrhaki.com/2012/04/groovy-goodness-collate-list-into-sub.html)). This method segments a list into sub-lists of specified size.

        You may also need the `collect()` method, which will perform an operation on every element in the collection (array, tuple, etc) and return a new collection of the outputs. This is similar to the `map` operator for but works on collections instead of channels ([docs](http://docs.groovy-lang.org/2.4.3/html/groovy-jdk/java/util/Collection.html#collect(groovy.lang.Closure)), [tutorial](https://www.baeldung.com/groovy-lists#Collecting)).

    ??? solution
        This is a sensible approach. For each item in the channel, the `flatmap` will collate the files into groups of 3 using `.collate(3)`, then collect them into a single tuple using the id as the first value and the collection of 3 files as the second value.

        ```groovy linenums="1" hl_lines="6-11"
        workflow {
            channel.fromPath("data/datfiles/sample*/*.dat", checkIfExists: true)
                .map { myfile -> [myfile.getParent().name, myfile] }
                .groupTuple()
                .flatMap { id, files ->
                    files
                        .collate(3)
                        .collect { chunk -> [ id, chunk ]
                    }
                }
                .view()
        }
        ```

## `collectFile`

The `collectFile` operator allows you to write one or more new files based on the contents of a channel.

### Writing strings to a single file

At its most basic, this operator writes the contents of the elements of a channel directly into a file. If the objecting being passed into the `collectFile` operator is a string, the strings are written to the collected file:

```groovy linenums="1"
workflow {
    characters = channel.of(
        ['name': 'Jake', 'title': 'Detective'],
        ['name': 'Rosa', 'title': 'Detective'],
        ['name': 'Terry', 'title': 'Sergeant'],
        ['name': 'Amy', 'title': 'Detective'],
        ['name': 'Charles', 'title': 'Detective'],
        ['name': 'Gina', 'title': 'Administrator'],
        ['name': 'Raymond', 'title': 'Captain'],
        ['name': 'Michael', 'title': 'Detective'],
        ['name': 'Norm', 'title': 'Detective']
    )

    characters
        .map { character -> character.name }
        .collectFile()
        .view()
}
```

The operator returns a channel containing a new file `collect-file.data`:

```console title="Output"
N E X T F L O W  ~  version 23.10.1
Launching `./main.nf` [distracted_ritchie] DSL2 - revision: 1e5061406d
work/tmp/57/8052c566f657c59679f07958031231/collect-file.data
```

If we view the contents of this file, we see the strings written (without delimiter) to a single file:

```console title="collect-file.data"
RaymondTerryNormRosaCharlesJakeGinaMichaelAmy
```

We can supply arguments `name` and `newLine` to the `collectFile` operator to return a file with a more informative name and newlines separating each entry:

```groovy linenums="14" hl_lines="3"
characters
    .map { it.name }
    .collectFile(name: 'people.txt', newLine: true)
    .view()
```

```console title="people.txt"
Raymond
Terry
Norm
Rosa
Charles
Jake
Gina
Michael
Amy
```

### Collecting to a new location

By default, the collected file is written into the work directory, which makes it suitable for input into a downstream process. If the collected file is an output of the workflow instead of an intermediate, it can be written to a directory of your choosing using the `storeDir` argument:

```groovy linenums="14" hl_lines="3"
characters
    .map { it.name }
    .collectFile(name: 'characters.txt', newLine: true, storeDir: 'results')
    .view()
```

### Collecting file contents

If the contents of the input channel is a file, its _contents_ are appended to the collected file. Here we make a small process `WriteBio` that generates a CSV from the Map object supplied as input:

!!! note "Groovy in the script block"

    In the example below, we include a line of groovy to define a variable `article` which is used in the interpolated script string. This is a convenient way to avoid crowding the final string block with too much logic.

    This line includes two Groovy syntax features:

    1. The [ternary operator](https://docs.groovy-lang.org/latest/html/documentation/core-operators.html#_ternary_operator) - a terse if/else block
    2. The [find operator](https://docs.groovy-lang.org/latest/html/documentation/core-operators.html#_find_operator) `=~`

```groovy linenums="1"
process WriteBio {
    input: val(character)
    output: path('bio.txt')
    script:
    def article = character.title.toLowerCase() =~ ~/^[aeiou]/ ? 'an' : 'a'
    """
    echo ${character.name} is ${article} ${character.title} > bio.txt
    """
}

workflow {
    characters = channel.of(
        ['name': 'Jake', 'title': 'Detective'],
        ['name': 'Rosa', 'title': 'Detective'],
        ['name': 'Terry', 'title': 'Sergeant'],
        ['name': 'Amy', 'title': 'Detective'],
        ['name': 'Charles', 'title': 'Detective'],
        ['name': 'Gina', 'title': 'Administrator'],
        ['name': 'Raymond', 'title': 'Captain'],
        ['name': 'Michael', 'title': 'Detective'],
        ['name': 'Norm', 'title': 'Detective']
    )

    WriteBio(characters)
      .collectFile()
      .view()
}
```

... to produce the collected file

```text title="bio.txt"
Charles is a Detective
Amy is a Detective
Norm is a Detective
Jake is a Detective
Rosa is a Detective
Raymond is a Captain
Gina is an Administrator
Terry is a Sergeant
Michael is a Detective
```

### Collecting into multiple files

Instead of writing all entries to a single file, you can direct entries from the input channel to different files by supplying a closure to the `collectFile` operator. The closure _must_ return a `List` of two entries where the two elements in the `List` are

1. the name of the file into which the data should be written, and
2. the data to write.

For example:

```groovy linenums="12" hl_lines="7"
characters
    .collectFile(newLine: true, storeDir: 'results') { character ->
      def filename = "${character.title}s.txt"
      def article = character.title.toLowerCase() =~ ~/^[aeiou]/ ? 'an' : 'a'
      def text = "${character.name} is ${article} ${character.title}"
        [filename, text]
    }
    .view()
```

The `collectFile` operator now returns a channel containing four files where we have the outputs grouped by the filename we specified.

```output
Launching `./main.nf` [marvelous_legentil] DSL2 - revision: a571cfd449
results/Administrators.txt
results/Captains.txt
results/Detectives.txt
results/Sergeants.txt
```

### Dealing with headers

Let's say we have a process that returns a CSV. In this example, we're going to create very small two-line CVSs, but the same approach is applicable to larger files as well.

```groovy linenums="1"
process WriteBio {
    input: val(character)
    output: path('bio.csv')
    script:
    """
    echo "precinct,name,title" > bio.csv
    echo 99th,${character.name},${character.title} >> bio.csv
    """
}
```

If we run this with the same workflow as before:

```groovy linenums="24"
WriteBio(characters)
    .collectFile(name: 'characters.csv', storeDir: 'results')
    .view()
```

... the CSVs are simply concatenated with the header included each time:

```text title="results/characters.csv"
precinct,name,title
99th,Jake,Detective
precinct,name,title
99th,Terry,Sergeant
precinct,name,title
99th,Rosa,Detective
...
```

To keep the header from only the first entry, we can use the `keepHeader` argument to `collectFile`:

```groovy linenums="24" hl_lines="2"
WriteBio(characters)
    .collectFile(name: 'characters.csv', storeDir: 'results', keepHeader: true)
    .view()
```

!!! exercise

    If we modify the `WriteBio` to also emit the `character` Map into the output channel:

    ```groovy linenums="1"
    process WriteBio {
        input: val(character)
        output: tuple val(character), path('bio.csv')
        script:
        """
        echo "precinct,name,title" > bio.csv
        echo 99th,${character.name},${character.title} >> bio.csv
        """
    }
    ```

    How might we produce a `results` directory that has one csv for each character title/rank where the csv includes the appropriate header?

    ```output
    results
    ├── Administrators.csv
    ├── Captains.csv
    ├── Detectives.csv
    └── Sergeants.csv
    ```

    ??? solution
        A good solution would be to pass a closure to the `collectFile` operator. The closure will return the filename and the file in a List:

        ```groovy linenums="24" hl_lines="2-4"
        WriteBio(characters)
            .collectFile(storeDir: 'results', keepHeader: true) { character, file ->
                ["${character.title}s.csv", file]
            }
            .view()
        ```

        Another viable option would be to `map` over the channel before `collectFile`:

        ```groovy linenums="24" hl_lines="2"
        WriteBio(characters)
            .map { character, file -> ["${character.title}s.csv", file] }
            .collectFile(storeDir: 'results', keepHeader: true)
            .view()
        ```
