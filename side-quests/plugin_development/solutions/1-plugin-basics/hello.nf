#!/usr/bin/env nextflow

include { randomString } from 'plugin/nf-hello'
include { samplesheetToList } from 'plugin/nf-schema'

params.input = 'greetings.csv'

process SAY_HELLO {
    input:
        val greeting
    output:
        stdout
    script:
    """
    echo '$greeting'
    """
}

workflow {
    def samplesheet_list = samplesheetToList(params.input, 'greetings_schema.json')
    greeting_ch = Channel.fromList(samplesheet_list)
        .map { row -> "${row[0]}_${randomString(8)}" }
    SAY_HELLO(greeting_ch)
    SAY_HELLO.out.view { result -> "Output: ${result.trim()}" }
}
