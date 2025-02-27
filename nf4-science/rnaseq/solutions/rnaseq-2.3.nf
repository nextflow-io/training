#!/usr/bin/env nextflow

// Module INCLUDE statements
include { FASTQC } from './modules/fastqc'
include { TRIM_GALORE } from './modules/trim_galore'
include { HISAT2_ALIGN } from './modules/hisat2_align'

/*
 * Pipeline parameters
 */
params.hisat2_index = "path/to/hisat2/index"
params.splice_sites = "path/to/splice_sites.txt"

// Primary input
params.reads = "path/to/reads/*.fastq.gz"

workflow {
    // Create input channel from a file path
    read_ch = channel.fromPath(params.reads)

    // Initial quality control
    FASTQC(read_ch)

    // Adapter trimming and post-trimming QC
    TRIM_GALORE(read_ch)

    // Alignment to a reference genome
    HISAT2_ALIGN(TRIM_GALORE.out.trimmed_reads, params.hisat2_index, params.splice_sites)

}
