#!/usr/bin/env nextflow

// Import custom functions from our plugin
include { decorateGreeting } from 'plugin/nf-greeting'
include { shoutAll } from 'plugin/nf-greeting'

params.input = 'greetings.csv'

process SAY_HELLO {
    input:
        val greeting

    output:
        stdout

    script:
    def decorated = decorateGreeting(greeting)
    """
    echo '$decorated'
    """
}

workflow {
    greeting_ch = channel.fromPath(params.input)
                        .splitCsv(header: true)
                        .map { row -> row.greeting }

    // Demonstrate using the shoutAll operator
    greeting_ch
        .shoutAll()
        .view { shouted -> "SHOUTED: $shouted" }

    SAY_HELLO(greeting_ch).view{ result -> "Decorated with custom prefix: ${result.trim()}" }
}
