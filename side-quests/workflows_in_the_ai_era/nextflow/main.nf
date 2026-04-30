#!/usr/bin/env nextflow

/*
 * RNA-seq Analysis Pipeline - Nextflow Version
 */

// Include processes from modules
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
    // Each sample becomes [meta, [read1, read2]]
    ch_samples = Channel
        .fromPath(params.samples)
        .splitCsv(header: true)
        .map { row ->
            def meta = [id: row.sample]
            def reads = [file(row.fastq_1), file(row.fastq_2)]
            return [meta, reads]
        }

    // Prepare salmon index
    ch_salmon_index = Channel.fromPath(params.salmon_index)
    UNTAR(ch_salmon_index)

    // TODO: Call FASTQC process with ch_samples

    // TODO: Call FASTP process with ch_samples

    // TODO: Call SALMON_QUANT with FASTP output and the index
    // Hint: Use FASTP.out.reads for the trimmed reads
    // Hint: Use UNTAR.out.index.first() for the index
}
