#!/usr/bin/env nextflow

/*
 * RNA-seq Analysis Pipeline - Nextflow Version
 */

// Include processes
include { FASTQC } from './modules/fastqc'
include { FASTP } from './modules/fastp'
include { UNTAR; SALMON_QUANT } from './modules/salmon'

// Pipeline parameters
params.samples = '../data/samples.csv'
params.outdir = '../results'
params.salmon_index = 'https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq/reference/salmon.tar.gz'

// Main workflow
workflow {
    // Create channel from sample sheet
    ch_samples = Channel
        .fromPath(params.samples)
        .splitCsv(header: true)
        .map { row ->
            def meta = [id: row.sample]
            def reads = [file(row.fastq_1), file(row.fastq_2)]
            return [meta, reads]
        }

    // Run processes
    FASTQC(ch_samples)
    FASTP(ch_samples)

    // Salmon needs the index AND fastp's output
    ch_salmon_index = Channel.fromPath(params.salmon_index)
    UNTAR(ch_salmon_index)

    SALMON_QUANT(FASTP.out.reads, UNTAR.out.index.first())
}
