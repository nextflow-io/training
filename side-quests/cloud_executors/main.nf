#!/usr/bin/env nextflow

/*
 * Simple greeting pipeline for demonstrating executor configuration
 */

params.input = 'greetings.csv'

process SAY_HELLO {

    input:
        val greeting

    output:
        path "${greeting}-output.txt"

    script:
    """
    echo '$greeting' > '${greeting}-output.txt'
    """
}

process CONVERT_TO_UPPER {

    input:
        path input_file

    output:
        path "UPPER-${input_file}"

    script:
    """
    cat '$input_file' | tr '[a-z]' '[A-Z]' > 'UPPER-${input_file}'
    """
}

workflow {

    greeting_ch = channel.fromPath(params.input)
                        .splitCsv(header: true)
                        .map { row -> row.greeting }

    SAY_HELLO(greeting_ch)
    CONVERT_TO_UPPER(SAY_HELLO.out)
}
