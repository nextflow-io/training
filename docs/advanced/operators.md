# Operator Tour

In this chapter, we take a curated tour of the Nextflow operators. Commonly used and well understood operators are not covered here - only those that we've seen could use more attention or those where the usage could be more elaborate.

These modest set of operators have been chosen to simultaneously demo tangential concepts and Nextflow features.

## `map`

### Basics

Map is certainly the most commonly used of the operators covered here. It's a way to supply a closure through which each element in the channel is passed. The return value of the closure is emitted as an element in a new output channel. A canonical example is a closure that multiplies two numbers:

```groovy linenums="1"
workflow {
    Channel.of( 1, 2, 3, 4, 5 )
    | map { it * it }
    | view
}
```

By default, the element being passed to the closure is given the default name `it`. The variable can be named by using the `->` notation:

```groovy linenums="1"
workflow {
    Channel.of( 1, 2, 3, 4, 5 )
    | map { num -> num * num }
    | view
}
```

Groovy is an optionally typed language, and it is possible to specify the type of the argument passed to the closure.

```groovy linenums="1"
workflow {
    Channel.of( 1, 2, 3, 4, 5 )
    | map { Integer num -> num * num }
    | view
}
```

### Named Closures

If you find yourself re-using the same closure multiple times in your pipeline, the closure can be named and referenced:

```groovy linenums="1"
def squareIt = { Integer num -> num * num }

workflow {
    Channel.of( 1, 2, 3, 4, 5 )
    | map( squareIt )
    | view
}
```

If you have these re-usable closures defined, you can compose them together.

```groovy linenums="1"
def squareIt = { it * it }
def addTwo = { it + 2 }

workflow {
    Channel.of( 1, 2, 3, 4, 5 )
    | map( squareIt >> addTwo )
    | view
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

```groovy linenums="1"
def squareIt = { it * it }
def addTwo = { it + 2 }

workflow {
    Channel.of( 1, 2, 3, 4, 5 )
    | map( squareIt )
    | map( addTwo )
    | view
}
```

For those inclined towards functional programming, you'll be happy to know that closures can be curried:

```groovy linenums="1"
def timesN = { multiplier, it -> it * multiplier }
def timesTen = timesN.curry(10)

workflow {
    Channel.of( 1, 2, 3, 4, 5 )
    | map( timesTen )
    | view
}
```

## `view`

In addition to the argument-less usage of `view` as shown above, this operator can also take a closure to customize the stdout message. We can create a closure to print the value of the elements in a channel as well as their type, for example:

```groovy linenums="1"
def timesN = { multiplier, it -> it * multiplier }
def timesTen = timesN.curry(10)
def prettyPrint = { "Found '$it' (${it.getClass()})"}

workflow {
    Channel.of( 1, 2, 3, 4, 5 )
    | map( timesTen )
    | view( prettyPrint )
}
```

!!! note "Most closures will remain anonymous"

    In many cases, it is simply cleaner to keep the closure anonymous, defined inline. Giving closures a name is only recommended when you find yourself defining the same or similar closures repeatedly in a given workflow.

## `splitCsv`

It is common that a samplesheet is passed as input into a Nextflow workflow. We'll see some more complicated ways to manage these inputs later on in the workshop, but the `splitCsv` is an excellent tool to have in a pinch.

```groovy linenums="1"
workflow {
    Channel.fromPath("data/samplesheet.csv")
    | splitCsv( header: true )
    | view
}
```

!!! exercise

    From the directory `chapter_01_operators`, use the `splitCsv` and `map` operators to create a channel that would be suitable input to the

    ```groovy linenums="1"
    process FastQC {
        input:
        tuple val(id), path(fastqs)
        //
    ```

    ??? solution
        Specifying the `header` argument in the `splitCsv` operator, we have convenient named access to csv elements. The closure returns a list of two elements where the second element a list of paths.

        ```groovy linenums="1"
        workflow {
            Channel.fromPath("data/samplesheet.csv")
            | splitCsv( header: true )
            | map { row ->
                [row.id, [file(row.fastq1), file(row.fastq2)]]
            }
            | view
        }
        ```

        !!! warning "Convert Strings to Paths"

            The fastq paths are simple strings in the context of a csv row. In order to pass them as paths to a Nextflow process, they need to be converted into objects that adjere to the `Path` interface. This is accomplished by wrapping them in `file`.

        In the sample above, we've lost an important piece of metadata - the tumor/normal classification, choosing only the sample id as the first element in the output list.

        In the next chapter, we'll discuss the "meta map" pattern in more detail, but we can preview that here.

        ```groovy linenums="1"
        workflow {
            Channel.fromPath("data/samplesheet.csv")
            | splitCsv( header: true )
            | map { row ->
                metaMap = [id: row.id, type: row.type, repeat: row.repeat]
                [metaMap, [file(row.fastq1), file(row.fastq2)]]
            }
            | view
        }
        ```

        The construction of this map is very repetitive, and in the next chapter, we'll discuss some Groovy methods available on the `Map` class that can make this pattern more concise and less error-prone.

## `multiMap`

The `multiMap` operator is a way of creating multiple channels from a single source.

Let's assume we've been given a samplesheet that has tumor/normal pairs bundled together on the same row.

```bash
cd chapter_01_operators
cat data/samplesheet.ugly.csv
```

Using the `splitCsv` operator would give us one entry that would contain all four fastq files. Let's consider that we wanted to split these fastqs into separate channels for tumor and normal, we could use `multiMap`:

```groovy linenums="1"
workflow {
    Channel.fromPath("data/samplesheet.ugly.csv")
    | splitCsv( header: true )
    | multiMap { row ->
        tumor:
            metamap = [id: row.id, type:'tumor', repeat:row.repeat]
            [metamap, file(row.tumor_fastq_1), file(row.tumor_fastq_2)]
        normal:
            metamap = [id: row.id, type:'normal', repeat:row.repeat]
            [metamap, file(row.normal_fastq_1), file(row.normal_fastq_2)]
    }
    | set { samples }

    samples.tumor | view { "Tumor: $it"}
    samples.normal | view { "Normal: $it"}
}
```

!!! tip "multiMapCriteria"

    The closure supplied to `multiMap` needs to return multiple channels, so using named closures as described in the `map` section above will not work. Fortunately, Nextflow provides the convenience `multiMapCriteria` method to allow you to define named `multiMap` closures should you need them. See the [`multiMap` documentation](https://www.nextflow.io/docs/latest/operator.html#multimap) for more info.

## `branch`

In the example above, the `multiMap` operator was necessary because we were supplied with a samplesheet that combined two pairs of fastq per row. If we were to use the neater samplesheet, we could use the `branch` operator to achieve the same result.

```groovy linenums="1"
workflow {
    Channel.fromPath("data/samplesheet.csv")
    | splitCsv( header: true )
    | map { row -> [[id: row.id, repeat: row.repeat, type: row.type], [file(row.fastq1), file(row.fastq2)]] }
    | branch { meta, reads ->
        tumor: meta.type == "tumor"
        normal: meta.type == "normal"
    }
    | set { samples }

    samples.tumor | view { "Tumor: $it"}
    samples.normal | view { "Normal: $it"}
}
```

An element is only emitted to the first channel were the test condition is met. If an element does not meet any of the tests, it is not emitted to any of the output channels. You can 'catch' any such samples by specifying `true` as a condition. If we knew that all samples would be either tumor or normal and no third 'type', we could write

```groovy linenums="1"
branch { meta, reads ->
    tumor: meta.type == "tumor"
    normal: true
}
```

We can optionally return a new element to one or more of the output channels. For example, to add an extra key in the meta map of the tumor samples, we add a new line under the condition and return our new element. In this example, we modify the first element of the `List` to be a new list that is the result of merging the existing meta map with a new map containing a single key:

```groovy linenums="1"
branch { meta, reads ->
    tumor: meta.type == "tumor"
        return [meta + [newKey: 'myValue'], reads]
    normal: true
}
```

!!! exercise

    How would you modify the element returned in the `tumor` channel to have the key:value pair `type:'abnormal'` instead of `type:'tumor'`?

    ??? solution

        There are many ways to accomplish this, but the map merging pattern introduced above can also be used to safely and concisely rename values in a map.

        ```groovy linenums="1"
        branch { meta, reads ->
            tumor: meta.type == "tumor"
                return [meta + [type: 'abnormal'], reads]
            normal: true
        }
        ```

        !!! note "Merging maps is safe"

            Using the `+` operator to merge two or more Maps returns a _new_ Map. There are rare edge cases where modification of map rather than returning a new map can affect other channels. We discuss this further in the next chapter, but just be aware that this `+` operator is safer and often more convenient than modifying the `meta` object directly.

            See the Groovy [Map documentation](https://docs.groovy-lang.org/latest/html/groovy-jdk/java/util/Map.html#plus(java.util.Map)) for details.

### Multi-channel Objects

Some Nextflow operators return objects that contain _multiple_ channels. The `multiMap` and `branch` operators are excellent examples. In most instances, the output is assigned to a variable and then addressed by name:

```groovy linenums="1"
numbers = Channel.from(1,2,3,4,5)
| multiMap {
    small: it
    large: it * 10
}
numbers.small | view { num -> "Small: $num"}
numbers.large | view { num -> "Large: $num"}
```

or by using `set`:

```groovy linenums="1"
Channel.from(1,2,3,4,5)
| multiMap {
    small: it
    large: it * 10
}
| set { numbers }

numbers.small | view { num -> "Small: $num"}
numbers.large | view { num -> "Large: $num"}
```

Given a process that takes multiple channels

```groovy linenums="1"
process MultiInput {
    debug true
    input:
    val(smallNum)
    val(bigNum)

    "echo -n small is $smallNum and big is $bigNum"
}
```

You can either provide the channels individually:

```groovy linenums="1"
Channel.from(1,2,3,4,5)
| multiMap {
    small: it
    large: it * 10
}
| set { numbers }

MultiInput(numbers.small, numbers.large)
```

or you can provide the multichannel as a single input:

```groovy linenums="1"
Channel.from(1,2,3,4,5)
| multiMap {
    small: it
    large: it * 10
}
| set { numbers }

MultiInput(numbers)
```

This also means you can skip the `set` operator for the cleanest solution:

```groovy linenums="1"
Channel.from(1,2,3,4,5)
| multiMap {
    small: it
    large: it * 10
}
| MultiInput
```

If you have processes that output multiple channels and input multiple channels and the cardinality matches, they can be chained together in the same manner.

## `transpose`

The transpose operator is often misunderstood. It can be thought of as the inverse of the `groupTuple` operator. Give the following workflow, the `groupTuple` and `transpose` operators cancel each other out. Removing lines 8 and 9 returns the same result.

Given a workflow that returns one element per sample, where we have grouped the samplesheet lines on a meta containing only id and type:

```groovy linenums="1"
workflow {
    Channel.fromPath("data/samplesheet.csv")
    | splitCsv(header: true)
    | map { row ->
        meta = [id: row.id, type: row.type]
        [meta, row.repeat, [row.fastq1, row.fastq2]]
    }
    | groupTuple
    | view
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

```groovy linenums="1"
workflow {
    Channel.fromPath("data/samplesheet.csv")
    | splitCsv(header: true)
    | map { row ->
        meta = [id: row.id, type: row.type]
        [meta, row.repeat, [row.fastq1, row.fastq2]]
    }
    | groupTuple
    | transpose
    | view
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
