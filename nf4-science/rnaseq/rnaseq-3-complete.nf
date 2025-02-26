#!/usr/bin/env nextflow

include { FASTQC as FASTQC_RAW } from './modules/fastqc'
include { FASTQC as FASTQC_TRIMMED } from './modules/fastqc'
include { TRIM_GALORE } from './modules/trim_galore'
include { HISAT2_ALIGN } from './modules/hisat2_align'
include { MULTIQC as MULTIQC_QC } from './modules/multiqc'
include { MULTIQC as MULTIQC_ALIGN } from './modules/multiqc'
include { MULTIQC as MULTIQC_FINAL } from './modules/multiqc'

workflow {
    // Read the CSV file and create a channel
    read_ch = Channel.fromPath(params.input_csv)
        .splitCsv(header:true)
        .map { row -> tuple(row.sample_id, file(row.fastq_path)) }

    // Part 1: QC and Trimming
    FASTQC_RAW(read_ch)
    TRIM_GALORE(read_ch)
    FASTQC_TRIMMED(TRIM_GALORE.out.trimmed_reads)

    MULTIQC_QC(
        FASTQC_RAW.out.zip.mix(
        FASTQC_RAW.out.html,
        FASTQC_TRIMMED.out.zip,
        FASTQC_TRIMMED.out.html,
        TRIM_GALORE.out.trimming_reports
        ).collect(),
        'multiqc_qc'
    )

    // Part 2: Alignment
    HISAT2_ALIGN(TRIM_GALORE.out.trimmed_reads, params.hisat2_index, params.splice_sites)
    MULTIQC_ALIGN(HISAT2_ALIGN.out.log.collect(), 'multiqc_align')

    // Part 3: Final MultiQC Report
    MULTIQC_FINAL(
        MULTIQC_QC.out.report.mix(
        MULTIQC_QC.out.data,
        MULTIQC_ALIGN.out.report,
        MULTIQC_ALIGN.out.data,
        HISAT2_ALIGN.out.bam
        ).collect(),
        'multiqc_final'
    )
}
