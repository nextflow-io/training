#!/usr/bin/env nextflow

// Module INCLUDE statements
include { FASTQC } from './modules/fastqc_pe.nf'
include { TRIM_GALORE } from './modules/trim_galore_pe.nf'
include { HISAT2_ALIGN } from './modules/hisat2_align_pe.nf'
include { MULTIQC } from './modules/multiqc.nf'

/*
 * Pipeline parameters
 */
params.hisat2_index_zip = "data/genome_index.tar.gz"
params.report_id = "all_paired-end"

// Primary input
params.input_csv = "data/paired-end1.csv"
// params.input_csv = "data/single-end.csv"
// params.reads = "data/reads/ENCSR000COQ1_1.fastq.gz"




workflow {

    // Create input channel from a file path - https://nextflow.io/docs/latest/reference/channel.html#frompath
    read_ch = channel.fromPath(params.input_csv)
        .splitCsv(header:true)   // operator to split CSV into rows
        .map { row -> [file(row.fastq_1), file(row.fastq_2)] } // operator to grab column data https://nextflow.io/docs/latest/reference/operator.html#map
    //    .view { fastq1f, fastq2f -> "Extracted Fastq1 and Fastq2 files are $fastq1f and $fastq2f" }

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
        TRIM_GALORE.out.fastqc_reports_1,
        TRIM_GALORE.out.fastqc_reports_2,
        HISAT2_ALIGN.out.log
        ).collect(),
        params.report_id
    )

}
