#!/usr/bin/env nextflow

/*
 * Simple greeting pipeline
 * Will be enhanced to use custom plugin functions
 */

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
    greeting_ch = channel.fromPath(params.input)
                        .splitCsv(header: true)
                        .map { row -> row.greeting }

    SAY_HELLO(greeting_ch)

    SAY_HELLO.out.view { result -> "Output: ${result.trim()}" }
}
