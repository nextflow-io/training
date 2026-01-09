#!/usr/bin/env nextflow

// Import custom functions from our plugin
include { reverseGreeting; decorateGreeting } from 'plugin/nf-greeting'

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
    greeting_ch = channel.fromPath(params.input)
                        .splitCsv(header: true)
                        .map { row -> row.greeting }

    // Demonstrate using reverseGreeting function
    greeting_ch
        .map { greeting -> reverseGreeting(greeting) }
        .view { reversed -> "Reversed: $reversed" }

    SAY_HELLO(greeting_ch)
    SAY_HELLO.out.view { result -> "Decorated: ${result.trim()}" }
}
