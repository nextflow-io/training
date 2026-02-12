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
    input: Path

    // Reference genome archive
    hisat2_index_zip: Path

    // Report ID
    report_id: String
}

workflow {

    main:
    // Create input channel from the contents of a CSV file
    read_ch = channel.fromPath(params.input)
        .splitCsv(header: true)
        .map { row -> file(row.fastq_path) }

    // Initial quality control
    FASTQC(read_ch)

    // Adapter trimming and post-trimming QC
    TRIM_GALORE(read_ch)

    // Alignment to a reference genome
    HISAT2_ALIGN(TRIM_GALORE.out.trimmed_reads, file(params.hisat2_index_zip))

    // Comprehensive QC report generation
    multiqc_files_ch = channel.empty().mix(
        FASTQC.out.zip,
        FASTQC.out.html,
        TRIM_GALORE.out.trimming_reports,
        TRIM_GALORE.out.fastqc_reports,
        HISAT2_ALIGN.out.log,
    )
    multiqc_files_list = multiqc_files_ch.collect()
    MULTIQC(multiqc_files_list, params.report_id)

    publish:
    fastqc_zip = FASTQC.out.zip
    fastqc_html = FASTQC.out.html
    trimmed_reads = TRIM_GALORE.out.trimmed_reads
    trimming_reports = TRIM_GALORE.out.trimming_reports
    trimming_fastqc = TRIM_GALORE.out.fastqc_reports
    bam = HISAT2_ALIGN.out.bam
    align_log = HISAT2_ALIGN.out.log
    multiqc_report = MULTIQC.out.report
    multiqc_data = MULTIQC.out.data
}

output {
    fastqc_zip {
        path 'fastqc'
    }
    fastqc_html {
        path 'fastqc'
    }
    trimmed_reads {
        path 'trimming'
    }
    trimming_reports {
        path 'trimming'
    }
    trimming_fastqc {
        path 'trimming'
    }
    bam {
        path 'align'
    }
    align_log {
        path 'align'
    }
    multiqc_report {
        path 'multiqc'
    }
    multiqc_data {
        path 'multiqc'
    }
}
