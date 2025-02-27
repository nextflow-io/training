#!/usr/bin/env nextflow

// Module INCLUDE statements
include { FASTQC } from './modules/fastqc'
include { TRIM_GALORE } from './modules/trim_galore'
include { HISAT2_ALIGN } from './modules/hisat2_align'
include { MULTIQC } from './modules/multiqc'

/*
 * Pipeline parameters
 */
params.hisat2_index = "path/to/hisat2/index"
params.splice_sites = "path/to/splice_sites.txt"
params.report_id = "test_batch"

// Primary input
params.input_csv = "path/to/reads/*.csv"

workflow {
    // Create input channel from the contents of a CSV file
    read_ch = Channel.fromPath(params.input_csv)
        .splitCsv(header:true)
        .map { row -> file(row.fastq_path) }

    /// Initial quality control
    FASTQC(read_ch)

    // Adapter trimming and post-trimming QC
    TRIM_GALORE(read_ch)

    // Alignment to a reference genome
    HISAT2_ALIGN(TRIM_GALORE.out.trimmed_reads, params.hisat2_index, params.splice_sites)

    // Comprehensive QC report generation
    MULTIQC(
        FASTQC.out.zip.mix(
        FASTQC.out.html,
        TRIM_GALORE.out.trimming_reports,
        TRIM_GALORE.out.fastqc_reports,
        HISAT2_ALIGN.out.log
        ).collect(),
        params.report_id
    )
}
