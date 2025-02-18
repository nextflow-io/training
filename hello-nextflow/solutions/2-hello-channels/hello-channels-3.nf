#!/usr/bin/env nextflow

/*
 * Use echo to print 'Hello World!' to a file
 */
process sayHello {

    publishDir 'results', mode: 'copy'

    input:
        val greeting

    output:
        path "${greeting}-output.txt"

    script:
    """
    echo '$greeting' > '$greeting-output.txt'
    """
}

/*
 * Pipeline parameters
 */
params.greeting = 'Holà mundo'

workflow {

    greetings_array = ['Hello','Bonjour','Holà']

    // create a channel for inputs
    greeting_ch = Channel.of(greetings_array)
                    .view { "Before flatten: $it" }
                    .flatten()
                    .view { "After flatten: $it" }

    // emit a greeting
    sayHello(greeting_ch)
}
