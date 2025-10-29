# Grouping and Splitting

## Grouping using subMap

Now that we have a channel that conforms to the `tuple val(meta) path(<something>)` pattern, we can investigate splitting and grouping patterns.

We'll start with a simple main.nf in the `grouping` directory

```bash
cd grouping
```

```groovy linenums="1"
workflow {
    channel.fromPath("data/samplesheet.csv")
        .splitCsv( header:true )
        .map { row ->
            def meta = [id:row.id, repeat:row.repeat, type:row.type]
            [
                meta,
                [
                    file(row.fastq1, checkIfExists: true),
                    file(row.fastq2, checkIfExists: true)
                ]
            ]
        }
        .view()
}
```

The first change we're going to make is to correct some repetitive code that we've seen quite a lot already in this workshop. The construction of the meta map from this row stutters quite a lot. We can make use of the [`subMap`](<https://docs.groovy-lang.org/latest/html/groovy-jdk/java/util/Map.html#subMap(java.util.Collection)>) method available Maps to quickly return a new map constructed from the subset of an existing map:

```groovy linenums="1" hl_lines="5"
workflow {
    channel.fromPath("data/samplesheet.csv")
        .splitCsv( header:true )
        .map { row ->
            def meta = row.subMap('id', 'repeat', 'type')
            [
                meta,
                [
                    file(row.fastq1, checkIfExists: true),
                    file(row.fastq2, checkIfExists: true)
                ]
            ]
        }
        .view()
}
```

!!! note "Complete meta map safety"

    The `subMap` method will take a collection of keys and construct a _new_ map with just the keys listed in the collection. This method, in combination with the `plus` or `+` method for combining maps and resetting values should allow all contraction, expansion and modification of maps safely.

!!! exercise

    Can you extend our workflow in an _unsafe_ manner? Use the `set` operator to name the channel in our workflow above, and then `map` (the operator) over that without modification. In a separate map operation, try modifying the meta map in a way that is reflected in the first map.

    Note that we're trying to do the _wrong_ thing in this example to clarify what the correct approach might be.

    ??? solution

        To ensure that the modification of the map happens first, we introduce a `sleep` into the first map operation. This `sleep` emulates a long-running Nextflow process.

        ```groovy linenums="1"
        workflow {
            channel.fromPath("data/samplesheet.csv")
                .splitCsv( header:true )
                .map { row ->
                    def meta = row.subMap('id', 'repeat', 'type')
                    [
                        meta,
                        [
                            file(row.fastq1, checkIfExists: true),
                            file(row.fastq2, checkIfExists: true)
                        ]
                    ]
                }
                .set { samples }


            samples
                .map { element -> sleep 10; element }
                .view { meta, reads -> "Should be unmodified: $meta" }

            samples
                .map { meta, reads ->
                    meta.type = meta.type == "tumor" ? "abnormal" : "normal"
                    [meta, reads]
                }
                .view { meta, reads -> "Should be modified: $meta" }
        }
        ```

!!! exercise

    How would you fix the example above to use the safe operators `plus` and `subMap` to ensure that the original map remains unmodified?

    ??? solution

        ```groovy linenums="1"
        workflow {
            channel.fromPath("data/samplesheet.csv")
                .splitCsv( header:true )
                .map { row ->
                    def meta = row.subMap('id', 'repeat', 'type')
                    [
                        meta,
                        [
                            file(row.fastq1, checkIfExists: true),
                            file(row.fastq2, checkIfExists: true)
                        ]
                    ]
                }
                .set { samples }


            samples
                .map { element -> sleep 10; element }
                .view { meta, reads -> "Should be unmodified: $meta" }

            samples
                .map { meta, reads ->
                    def newmap = [type: meta.type == "tumor" ? "abnormal" : "normal"]
                    [meta + newmap, reads]
                }
                .view { meta, reads -> "Should be modified: $meta" }
        }
        ```

## Passing maps through processes

Let's construct a dummy read mapping process. This is not a bioinformatics workshop, so we can 'cheat' in the interests of time.

```groovy linenums="1" hl_lines="3 29"
process MapReads {
    input:
    tuple val(meta), path(reads)
    path(genome)

    output:
    tuple val(meta), path("*.bam")

    script:
    "touch out.bam"
}

workflow {
    reference = channel.fromPath("data/genome.fasta").first()

    samples = channel.fromPath("data/samplesheet.csv")
        .splitCsv( header:true )
        .map { row ->
            def meta = row.subMap('id', 'repeat', 'type')
            [
              meta,
              [
                file(row.fastq1, checkIfExists: true),
                file(row.fastq2, checkIfExists: true)
              ]
            ]
        }

    mapped_reads = MapReads( samples, reference )
    mapped_reads.view()
}
```

Let's consider that we might now want to merge the repeats. We'll need to group bams that share the `id` and `type` attributes.

```groovy linenums="29" hl_lines="2"
mapped_reads = MapReads( samples, reference )
    .map { meta, bam -> [meta.subMap('id', 'type'), bam]}
    .groupTuple()
mapped_reads.view()
```

This is easy enough, but the `groupTuple` operator has to wait until all items are emitted from the incoming queue before it is able to reassemble the output queue. If even one read mapping job takes a long time, the processing of all other samples is held up. We need a way of signalling to Nextflow how many items are in a given group so that items can be emitted as early as possible.

By default, the `groupTuple` operator groups on the first item in the element, which at the moment is a `Map`. We can turn this map into a special class using the `groupKey` method, which takes our grouping object as a first parameter and the number of expected elements in the second parameter.

```groovy linenums="29" hl_lines="3"
mapped_reads = MapReads( samples, reference )
    .map { meta, bam ->
        def key = groupKey(meta.subMap('id', 'type'), NUMBER_OF_ITEMS_IN_GROUP)
        [key, bam]
    }
    .groupTuple()
mapped_reads.view()
```

!!! exercise

    How might we modify the upstream channels to the number of repeats into the metamap?

    ??? solution

        ```groovy linenums="13" hl_lines="16-20"
        workflow {
            reference = channel.fromPath("data/genome.fasta").first()

            samples = channel.fromPath("data/samplesheet.csv")
                .splitCsv( header:true )
                .map { row ->
                    def meta = row.subMap('id', 'repeat', 'type')
                    [
                      meta,
                      [
                        file(row.fastq1, checkIfExists: true),
                        file(row.fastq2, checkIfExists: true)
                      ]
                    ]
                }
                .map { meta, reads -> [meta.subMap('id', 'type'), meta.repeat, reads] }
                .groupTuple()
                .map { meta, repeats, reads -> [meta + [repeatcount:repeats.size()], repeats, reads] }
                .transpose()
                .map { meta, repeat, reads -> [meta + [repeat:repeat], reads]}

            mapped_reads = MapReads( samples, reference )
                .map { meta, bam ->
                    def key = groupKey(meta.subMap('id', 'type'), meta.repeatcount)
                    [key, bam]
                }
                .groupTuple()
            mapped_reads.view()
        }
        ```

Now that we have our repeats together in an output channel, we can combine them using "advanced bioinformatics":

```groovy linenums="13"
process CombineBams {
    input:
    tuple val(meta), path("input/in_*_.bam")

    output:
    tuple val(meta), path("combined.bam")

    script:
    "cat input/*.bam > combined.bam"
}
```

In our workflow:

```groovy linenums="45" hl_lines="3"
mapped_reads = MapReads( samples, reference )
    .map { meta, bam ->
        def key = groupKey(meta.subMap('id', 'type'), meta.repeatcount)
        [key, bam]
    }
    .groupTuple()

CombineBams(mapped_reads)
  .view()
```

## Fanning out over intervals

The previous exercise demonstrated the fan-in approach using `groupTuple` and `groupKey`, but we might want to fan out our processes. An example might be computing over some intervals - genotyping over intervals, for example.

We can take an existing bed file, for example and turn it into a channel of Maps.

```groovy linenums="26"
intervals = channel.fromPath("data/intervals.bed")
    .splitCsv(header: ['chr', 'start', 'stop', 'name'], sep: '\t')
    .collectFile { entry -> ["${entry.name}.bed", entry*.value.join("\t")] }
    .view()
```

Given a dummy genotyping process:

```groovy linenums="24"
process GenotypeOnInterval {
    input:
    tuple val(meta), path(bam), path(bed)

    output:
    tuple val(meta), path("genotyped.vcf")

    script:
    "cat $bam $bed > genotyped.vcf"
}
```

We can use the `combine` operator to emit a new channel where each combined bam is attached to each bed file. These can then be piped into the genotyping process:

```groovy linenums="60" hl_lines="9"
mapped_reads = MapReads( samples, reference )
    .map { meta, bam ->
        def key = groupKey(meta.subMap('id', 'type'), meta.repeatcount)
        [key, bam]
    }
    .groupTuple()

combined_bams = CombineBams(mapped_reads)
    .combine( intervals )

genotyped_bams = GenotypeOnInterval(combined_bams)
    .view()
```

Finally, we can combine these genotyped bams back using `groupTuple` and another bam merge process. We construct our "merge" process that will combine the bam files from multiple intervals:

```groovy linenums="35"
process MergeGenotyped {
    input:
    tuple val(meta), path("input/in_*_.vcf")

    output:
    tuple val(meta), path("merged.genotyped.vcf")

    script:
    "cat input/*.vcf > merged.genotyped.vcf"
}
```

We might be tempted to pipe the output of `GenotypeOnInterval` directly into groupTuple, but the `meta` object we are passing down is still the `groupKey` we created earlier:

```groovy linenums="71" hl_lines="13"
mapped_reads = MapReads( samples, reference )
    .map { meta, bam ->
        def key = groupKey(meta.subMap('id', 'type'), meta.repeatcount)
        [key, bam]
    }
    .groupTuple()

combined_bams = CombineBams(mapped_reads)
    .map { meta, bam -> [meta.subMap('id', 'type'), bam] }
    .combine( intervals )

genotyped_bams = GenotypeOnInterval(combined_bams)
    .view { meta, bamfile -> "Meta is of ${meta.getClass()}" }
```

To ensure that grouping is performed only on the relevant elements, we can unwrap the `groupKey` to return the underlying `Map` using the `getGroupTarget()` method available on groupKeys. This allows the `groupTuple` operator to group by just the keys present in the map, similar to how `subMap` works. This approach ensures that downstream grouping and merging steps operate on the intended sample attributes.

```groovy linenums="71" hl_lines="13"
mapped_reads = MapReads( samples, reference )
    .map { meta, bam ->
        def key = groupKey(meta.subMap('id', 'type'), meta.repeatcount)
        [key, bam]
    }
    .groupTuple()

combined_bams = CombineBams(mapped_reads)
    .map { meta, bam -> [meta.subMap('id', 'type'), bam] }
    .combine( intervals )

genotyped_bams = GenotypeOnInterval(combined_bams)
    .map { groupKey, bamfile -> [groupKey.getGroupTarget(), bamfile] }
    .groupTuple()

merged_bams = MergeGenotyped(genotyped_bams)
merged_bams.view()
```

## Publishing the bams

This will return us six bam files - a tumor and normal pair for each of the three samples:

```output title="Final channel output"
[[id:sampleB, type:normal], merged.genotyped.vcf]
[[id:sampleB, type:tumor], merged.genotyped.vcf]
[[id:sampleA, type:normal], merged.genotyped.vcf]
[[id:sampleA, type:tumor], merged.genotyped.vcf]
[[id:sampleC, type:normal], merged.genotyped.vcf]
[[id:sampleC, type:tumor], merged.genotyped.vcf]
```

If we would like to save the output of our `MergeGenotyped` process, we can "publish" the outputs of a process using the `publishDir` directive. Try modifying the `MergeGenotyped` process to include the directive:

```groovy linenums="35" hl_lines="2"
process MergeGenotyped {
    publishDir 'results/genotyped'

    input:
    tuple val(meta), path("input/in_*_.vcf")

    output:
    tuple val(meta), path("merged.genotyped.vcf")

    script:
    "cat input/*.vcf > merged.genotyped.vcf"
}
```

This will publish all of the files in the `output` block of this process to the `results/genotyped` directory.

!!! tip "Workflow outputs"

    As of Nextflow 24.04 you can also manage result publication at the **workflow** level using the new `output { }` block. This allows you to publish files by pushing them from a channel instead of a process, similar to the channel operations we have been exploring in this workshop.

    See the [Nextflow documentation](https://www.nextflow.io/docs/latest/workflow.html#workflow-outputs) for a full description.

    ```groovy
    nextflow.preview.output = true

    workflow {
        main:
        merged_bams = MergeGenotyped(genotyped_bams)

        publish:
        results = merged_bams
    }

    output {
        results {
            path 'genotyped'
        }
    }
    ```

!!! exercise

    Inspect the contents of the `results` directory. Does this match what you were expecting? What is missing here?

    Can you modify the `MergeGenotyped` process to ensure we are capturing all of the expected output files?

    ??? solution

        One solution might be to modify the `script` block to ensure that each file has a unique name:

        ```groovy linenums="35" hl_lines="11"
        process MergeGenotyped {
            publishDir 'results/genotyped'

            input:
            tuple val(meta), path("input/in_*_.vcf")

            output:
            tuple val(meta), path("*.vcf")

            script:
            "cat input/*.vcf > ${meta.id}.${meta.type}.genotyped.vcf"
        }
        ```

        Another option might be to use the `saveAs` argument to the `publishDir` directive:

        ```groovy linenums="35" hl_lines="2"
        process MergeGenotyped {
            publishDir 'results/genotyped', saveAs: { "${meta.id}.${meta.type}.genotyped.vcf" }

            input:
            tuple val(meta), path("input/in_*_.vcf")

            output:
            tuple val(meta), path("merged.genotyped.vcf")

            script:
            "cat input/*.vcf > merged.genotyped.vcf"
        }
        ```

        If you are using the new output syntax described above, you can edit the filename using the `path` directive:

        ```groovy
        output {
            results {
                path { meta, vcf ->
                    "genotyped/${meta.id}.${meta.type}.genotyped.vcf"
                }
            }
        }
        ```
