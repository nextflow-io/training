#!/usr/bin/env nextflow

include { FASTP } from './modules/local/fastp/main.nf'

workflow {
    params.input = "https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq/samplesheet/v3.10/samplesheet_test.csv"

    Channel.fromPath(params.input)
    | splitCsv(header: true)
    | view
}
