#!/usr/bin/env nextflow

include { FASTQC as FASTQC_RAW } from './modules/fastqc'
include { FASTQC as FASTQC_TRIMMED } from './modules/fastqc'
include { TRIM_GALORE } from './modules/trim_galore'
include { HISAT2_ALIGN } from './modules/hisat2_align'
include { MULTIQC as MULTIQC_QC } from './modules/multiqc'
include { MULTIQC as MULTIQC_ALIGN } from './modules/multiqc'

workflow {
    read_ch = channel.fromPath(params.reads)

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
}
