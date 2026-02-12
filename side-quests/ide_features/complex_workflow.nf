#!/usr/bin/env nextflow

/*
 * Complex Workflow Example for IDE Navigation Training
 */

include { FASTQC } from './modules/fastqc'
include { STAR_ALIGN } from './modules/star'
include { MULTIQC } from './modules/utils'

process TRIM_GALORE {
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_trimmed_*.fq.gz")

    script:
    """
    trim_galore --paired ${reads}
    """
}

process FEATURECOUNTS {
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(bam), path(gtf)

    output:
    tuple val(sample_id), path("${sample_id}_counts.txt")

    script:
    """
    featureCounts -p -t exon -g gene_id -a ${gtf} -o ${sample_id}_counts.txt ${bam}
    """
}

workflow RNASEQ_PIPELINE {
    take:
    reads_ch
    reference_ch
    gtf_ch

    main:
    // Quality control
    FASTQC(reads_ch)

    // Trimming
    TRIM_GALORE(reads_ch)

    // Alignment
    STAR_ALIGN(TRIM_GALORE.out, reference_ch)

    // Quantification
    count_input = STAR_ALIGN.out.combine(gtf_ch)
    FEATURECOUNTS(count_input)

    // Aggregate QC
    qc_files = FASTQC.out.mix(TRIM_GALORE.out)
    MULTIQC(qc_files.collect())

    emit:
    fastqc = FASTQC.out
    trimmed = TRIM_GALORE.out
    alignments = STAR_ALIGN.out
    counts = FEATURECOUNTS.out
    qc_report = MULTIQC.out
}

workflow {
    main:
    // Input channels
    ch_reads = channel.fromPath(params.input)
        .splitCsv(header: true)
        .map { row -> [row.sample_id, [file(row.fastq_1), file(row.fastq_2)]] }

    ch_reference = channel.fromPath(params.reference)
    ch_gtf = channel.fromPath(params.gtf)

    // Execute pipeline
    RNASEQ_PIPELINE(ch_reads, ch_reference, ch_gtf)

    // Output summary
    RNASEQ_PIPELINE.out.counts
        .collectFile(name: 'expression_matrix.txt', sort: true)
        .view { "Expression matrix created: ${it}" }

    publish:
    fastqc = RNASEQ_PIPELINE.out.fastqc
    trimmed = RNASEQ_PIPELINE.out.trimmed
    alignments = RNASEQ_PIPELINE.out.alignments
    counts = RNASEQ_PIPELINE.out.counts
    qc_report = RNASEQ_PIPELINE.out.qc_report
}

output {
    directory params.output
    mode 'copy'

    fastqc {
        path 'fastqc'
    }
    trimmed {
        path 'trimmed'
    }
    alignments {
        path 'alignments'
    }
    counts {
        path 'counts'
    }
    qc_report {
        path 'multiqc'
    }
}
