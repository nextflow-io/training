#!/usr/bin/env nextflow

include { samplesheetToList } from 'plugin/nf-schema'
// Import custom functions from our plugin
include { reverseGreeting } from 'plugin/nf-greeting'
include { decorateGreeting } from 'plugin/nf-greeting'

params.input = 'greetings.csv'

process SAY_HELLO {
    input:
        val greeting
    output:
        stdout
    script:
    // Use our custom plugin function to decorate the greeting
    def decorated = decorateGreeting(greeting)
    """
    echo '$decorated'
    """
}

workflow {
    greeting_ch = Channel.fromList(samplesheetToList(params.input, 'greetings_schema.json'))
        .map { row -> row[0] }

    // Demonstrate using reverseGreeting function
    greeting_ch
        .map { greeting -> reverseGreeting(greeting) }
        .view { reversed -> "Reversed: $reversed" }

    SAY_HELLO(greeting_ch)
    SAY_HELLO.out.view { result -> "Decorated: ${result.trim()}" }
}
