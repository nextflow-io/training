# Groovy Imports

There exists in Groovy a wealth of helper classes that can be imported into Nextflow scripts. In this chapter, we create a very small Workflow using the FastP tool to investigate importing the Groovy JSONSlurper class.

First, let's move into the chapter 4 directory:

```bash
cd groovy
```

Let's assume that we would like to pull in a samplesheet, parse the entries and run them through the FastP tool. So far, we have been concerned with local files, but Nextflow will handle remote files transparently:

```groovy linenums="3"
params.input = "https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq/samplesheet/v3.10/samplesheet_test.csv"

workflow {

    channel.fromPath(params.input)
        .splitCsv(header: true)
        .view()
}
```

Let's write a small closure to parse each row into the now-familiar map + files shape. We might start by constructing the meta-map:

```groovy linenums="5" hl_lines="5-8"
workflow {

    samples = channel.fromPath(params.input)
        .splitCsv(header: true)
        .map { row ->
            def meta = row.subMap('sample', 'strandedness')
            meta
        }
        .view()
}
```

... but this precludes the possibility of adding additional columns to the samplesheet. We might to ensure the parsing will capture any extra metadata columns should they be added. Instead, let's partition the column names into those that begin with "fastq" and those that don't. Within the map closure, let's add an additional line to partition the column names:

```groovy linenums="10"
def (readKeys, metaKeys) = row.keySet().split { key -> key =~ /^fastq/ }
```

!!! note "New methods"

    We've introduced a new keySet method here. This is a method on Java's LinkedHashMap class ([docs here](https://docs.oracle.com/javase/8/docs/api/java/util/LinkedHashMap.html#keySet--))

    We're also using the `.split()` method, which divides collection based on the return value of the closure. The mrhaki blog [provides a succinct summary](https://blog.mrhaki.com/2009/12/groovy-goodness-splitting-with-closures.html).

From here, let's add another line collect the values of the read keys into a list of file objects:

```groovy linenums="11"
def reads = row.subMap(readKeys).values().collect { value -> file(value) }
```

... but we run into an error:

```groovy
Argument of `file` function cannot be empty
```

If we have a closer look at the samplesheet, we notice that not all rows have two read pairs. Let's add a condition to the collect method to only include the values that are not empty:

```groovy linenums="11"
def reads = row.subMap(readKeys).values()
    .findAll { value -> value != "" } // Single-end reads will have an empty string
    .collect { path -> file(path) }
```

Now we need to construct the meta map. Let's have a quick look at the FASTP module that I've already pre-defined:

```groovy linenums="1"
process FASTP {
    container 'quay.io/biocontainers/fastp:0.23.2--h79da9fb_0'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('*.fastp.fastq.gz') , optional:true, emit: reads
    tuple val(meta), path('*.json')           , emit: json

    script:
    def prefix = task.ext.prefix ?: meta.id
    if (meta.single_end) {
        // SNIP
    } else {
        // SNIP
    }
```

I can see that we require two extra keys, `id` and `single_end`:

```groovy linenums="14" hl_lines="1-3"
def meta = row.subMap(metaKeys)
meta = meta + [ id: meta.sample, single_end: reads.size == 1 ]
[meta, reads]
```

This is now able to be passed through to our FASTP process:

```groovy linenums="5" hl_lines="15 17"
workflow {

    samples = channel.fromPath(params.input)
        .splitCsv(header: true)
        .map { row ->
            def (readKeys, metaKeys) = row.keySet().split { key -> key =~ /^fastq/ }
            def reads = row.subMap(readKeys).values()
                .findAll { value -> value != "" } // Single-end reads will have an empty string
                .collect { path -> file(path) }
            def meta = row.subMap(metaKeys)
            meta = meta + [ id: meta.sample, single_end: reads.size == 1 ]
            [meta, reads]
        }

    FASTP(samples)

    FASTP.out.json.view()
}
```

Let's assume that we want to pull some information out of these JSON files. To make our lives a little more convenient, let's "publish" these json files so that they are more convenient. We're going to discuss configuration more completely in a later chapter, but that's no reason not to dabble a bit here.

We'd like to add a `publishDir` directive to our FASTP process.

```groovy linenums="3"
process {
    withName: 'FASTP' {
        publishDir = [
            path: { "results/fastp/json" },
            saveAs: { filename -> filename.endsWith('.json') ? filename : null },
        ]
    }
}
```

!!! note "Groovy Tip: Elvis Operator"

    This pattern of returning something if it is true and `somethingElse` if not:

    ```groovy linenums="1"
    somethingThatMightBeFalsey ? somethingThatMightBeFalsey : somethingElse
    ```

    has a shortcut in Groovy - the "Elvis" operator:

    ```groovy linenums="1"
    somethingThatMightBeFalsey ?: somethingElse
    ```

This enables us to iterate quickly to test out our JSON parsing without waiting on the FASTP caching to calculate on these slow virtual machines.

```bash
nextflow run . -resume
```

Let's consider the possibility that we'd like to capture some of these metrics so that they can be used downstream. First, we'll have a quick peek at the [Groovy docs](https://groovy-lang.org/documentation.html) and I see that I need to use `JsonSlurper`.

Now let's create a second entrypoint to quickly pass these JSON files through some tests:

!!! note "Entrypoint developing"

    Using a second Entrypoint allows us to do quick debugging or development using a small section of the workflow without disturbing the main flow.

```groovy linenums="5"
workflow Jsontest {
    channel.fromPath("results/fastp/json/*.json")
        .view()
}
```

which we run with

```bash
nextflow run . -resume -entry Jsontest
```

Let's create a small function inside the workflow to take the JSON path and pull out some basic metrics:

```groovy linenums="5"
def getFilteringResult(json_file) {
    return new groovy.json.JsonSlurper().parseText(json_file.text)
}

workflow Jsontest {
    channel.fromPath("results/fastp/json/*.json")
        .view()
}
```

The `fastpResult` returned from the `parseText` method is a large Map - a class which we're already familiar with. Modify the `getFilteringResult` function to return just the `after_filtering` section of the report.

In the interest of brevity, here is the solution to return just the `after_filtering` section of the report:

```groovy linenums="5"
def getFilteringResult(json_file) {
    return new groovy.json.JsonSlurper().parseText(json_file.text)
        ?.summary
        ?.after_filtering
}
```

!!! note

    `?.` is new notation is a null-safe access operator. The `?.summary` will access the summary property if the property exists.

We can then join this new map back to the original reads using the `join` operator:

```groovy linenums="31"
    FASTP.out.json
        .map { meta, json -> [meta, getFilteringResult(json)] }
        .join( FASTP.out.reads )
        .view()
}
```

!!! exercise

    Can you amend this pipeline to create two channels that filter the reads to exclude any samples where the Q30 rate is less than 93.5%?

    ??? solution

        ```groovy linenums="31"
            reads = FASTP.out.json
                .map { meta, json -> [meta, getFilteringResult(json)] }
                .join( FASTP.out.reads )
                .map { meta, fastpMap, reads -> [meta + fastpMap, reads] }
                .branch { meta, reads ->
                    pass: meta.q30_rate >= 0.935
                    fail: true
                }

            reads.fail.view { meta, _reads -> "Failed: ${meta.id}" }
            reads.pass.view { meta, _reads -> "Passed: ${meta.id}" }
        }
        ```
