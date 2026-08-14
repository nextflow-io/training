#!/usr/bin/env nextflow

/*
 * RNA-seq analysis pipeline (starter).
 *
 * You will fill in the TODO sections through the tutorial.
 */

include { FASTQC } from './modules/fastqc'
include { FASTP } from './modules/fastp'
include { UNTAR; SALMON_QUANT } from './modules/salmon'

params.samples = '../data/samples.csv'
params.salmon_index = 'https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq/reference/salmon.tar.gz'

workflow {

    main:
    // Each sample becomes [id, [read1, read2]]
    ch_samples = channel
        .fromPath(params.samples)
        .splitCsv(header: true)
        .map { row -> [row.sample, [file(row.fastq_1), file(row.fastq_2)]] }

    ch_salmon_index = channel.fromPath(params.salmon_index)
    UNTAR(ch_salmon_index)

    // TODO: Call FASTQC with ch_samples

    // TODO: Call FASTP with ch_samples

    // TODO: Call SALMON_QUANT with FASTP output and the index
    // Hint: Use FASTP.out.reads and UNTAR.out.index.first()

    publish:
    // TODO: Publish the outputs you want kept after the run.
    // Hint: each line names a published channel, e.g.
    //     fastqc_html = FASTQC.out.html
}

output {
    // TODO: Configure output paths for each published channel.
    // Hint: each block names a subdirectory, e.g.
    //     fastqc_html { path 'fastqc' }
}
