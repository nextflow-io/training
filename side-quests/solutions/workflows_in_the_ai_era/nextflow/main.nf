#!/usr/bin/env nextflow

/*
 * RNA-seq analysis pipeline.
 *
 * Quality control, adapter trimming, and transcript quantification
 * for paired-end RNA-seq samples listed in a CSV sample sheet.
 */

include { FASTQC } from './modules/fastqc'
include { FASTP } from './modules/fastp'
include { UNTAR; SALMON_QUANT } from './modules/salmon'

params.samples = '../data/samples.csv'
params.salmon_index = 'https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq/reference/salmon.tar.gz'

workflow {

    main:
    ch_samples = channel
        .fromPath(params.samples)
        .splitCsv(header: true)
        .map { row -> [row.sample, [file(row.fastq_1), file(row.fastq_2)]] }

    ch_salmon_index = channel.fromPath(params.salmon_index)
    UNTAR(ch_salmon_index)

    FASTQC(ch_samples)
    FASTP(ch_samples)
    SALMON_QUANT(FASTP.out.reads, UNTAR.out.index.first())

    publish:
    fastqc_html = FASTQC.out.html
    fastqc_zip = FASTQC.out.zip
    fastp_reads = FASTP.out.reads
    fastp_json = FASTP.out.json
    fastp_html = FASTP.out.html
    salmon_quant = SALMON_QUANT.out.results
}

output {
    fastqc_html {
        path 'fastqc'
    }
    fastqc_zip {
        path 'fastqc'
    }
    fastp_reads {
        path 'fastp'
    }
    fastp_json {
        path 'fastp'
    }
    fastp_html {
        path 'fastp'
    }
    salmon_quant {
        path 'salmon'
    }
}
