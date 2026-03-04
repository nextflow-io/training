#!/usr/bin/env nextflow

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
    greeting_ch = Channel.fromList(samplesheetToList(params.input, 'greetings_schema.json'))
        .map { row -> row[0] }
    SAY_HELLO(greeting_ch)
    SAY_HELLO.out.view { result -> "Output: ${result.trim()}" }
}
