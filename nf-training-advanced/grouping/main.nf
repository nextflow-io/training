process MapReads {
    input:
    tuple val(meta), path(reads)
    path(genome)

    output:
    tuple val(meta), path("*.bam")

    script:
    "touch out.bam"
}

process CombineBams {
    input:
    tuple val(meta), path("input/in_*_.bam")

    output:
    tuple val(meta), path("combined.bam")

    script:
    "cat input/*.bam > combined.bam"
}

process GenotypeOnInterval {
    input:
    tuple val(meta), path(bam), path(bed)

    output:
    tuple val(meta), path("genotyped.vcf")

    script:
    "cat $bam $bed > genotyped.vcf"
}

process MergeGenotyped {
    input:
    tuple val(meta), path("input/in_*_.vcf")

    output:
    tuple val(meta), path("merged.genotyped.vcf")

    script:
    "cat input/*.vcf > merged.genotyped.vcf"
}

workflow {
    reference = Channel.fromPath("data/genome.fasta").first()
    intervals = Channel.fromPath("data/intervals.bed")
        .splitCsv(header: ['chr', 'start', 'stop', 'name'], sep: '\t')
        .collectFile { entry -> ["${entry.name}.bed", entry*.value.join("\t")] }
        .view()

    samples = Channel.fromPath("data/samplesheet.csv")
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

    combined_bams = CombineBams(mapped_reads)
        .map { meta, bam -> [meta.subMap('id', 'type'), bam] }
        .combine( intervals )

    genotyped_bams = GenotypeOnInterval(combined_bams)
        .map { groupKey, bamfile -> [groupKey as Map, bamfile] }
        .groupTuple()

    merged_bams = MergeGenotyped(genotyped_bams)
    merged_bams.view()
}
