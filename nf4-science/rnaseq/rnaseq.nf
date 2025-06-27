#!/usr/bin/env nextflow

// Module INCLUDE statements
include { FASTQC } from './modules/fastqc.nf'
include { TRIM_GALORE } from './modules/trim_galore.nf'
include { HISAT2_ALIGN } from './modules/hisat2_align.nf'
include { MULTIQC } from './modules/multiqc.nf'

/*
 * Pipeline parameters
 */
params.hisat2_index_zip = "data/genome_index.tar.gz"
params.report_id = "all_single-end"

// Primary input
params.input_csv = "data/single-end.csv"
// params.reads = "data/reads/ENCSR000COQ1_1.fastq.gz"


workflow {

    // Create input channel from a file path
    read_ch = channel.fromPath(params.input_csv)
        .splitCsv(header:true)   // operator to split CSV into rows
        .map { row -> file(row.fastq_path) } // operator to grab the fastq_path column data

    // read_ch = channel.fromPath(params.reads)

    // Initial quality control
    FASTQC(read_ch)

    // Adapter trimming and post-trimming QC
    TRIM_GALORE(read_ch)

    // Alignment to a reference genome
    HISAT2_ALIGN(TRIM_GALORE.out.trimmed_reads, file (params.hisat2_index_zip))

    // Comprehensive QC report generation
    MULTIQC(
        FASTQC.out.zip.mix(    // mix combines output from multiple channels https://nextflow.io/docs/latest/reference/operator.html#mix
        FASTQC.out.html,
        TRIM_GALORE.out.trimming_reports,
        TRIM_GALORE.out.fastqc_reports,
        HISAT2_ALIGN.out.log
        ).collect(),
        params.report_id
    )

}
