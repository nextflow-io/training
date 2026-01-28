#!/usr/bin/env nextflow

// Module INCLUDE statements
include { FASTQC } from './modules/fastqc.nf'
include { TRIM_GALORE } from './modules/trim_galore.nf'
include { HISAT2_ALIGN } from './modules/hisat2_align.nf'
include { MULTIQC } from './modules/multiqc.nf'

/*
 * Pipeline parameters
 */
params {
    // Primary input
    input_csv: Path = "data/single-end.csv"

    // Reference genome archive
    hisat2_index_zip: Path = "data/genome_index.tar.gz"

    // Report ID
    report_id: String = "all_single-end"
}

workflow {
    // Create input channel from the contents of a CSV file
    read_ch = channel.fromPath(params.input_csv)
        .splitCsv(header: true)
        .map { row -> file(row.fastq_path) }

    /// Initial quality control
    FASTQC(read_ch)

    // Adapter trimming and post-trimming QC
    TRIM_GALORE(read_ch)

    // Alignment to a reference genome
    HISAT2_ALIGN(TRIM_GALORE.out.trimmed_reads, file(params.hisat2_index_zip))

    // Comprehensive QC report generation
    MULTIQC(
        FASTQC.out.zip.mix(
            FASTQC.out.html,
            TRIM_GALORE.out.trimming_reports,
            TRIM_GALORE.out.fastqc_reports,
            HISAT2_ALIGN.out.log,
        ).collect(),
        params.report_id,
    )
}
