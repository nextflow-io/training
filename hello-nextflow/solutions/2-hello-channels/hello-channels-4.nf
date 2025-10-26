#!/usr/bin/env nextflow

/*
 * Use echo to print 'Hello World!' to a file
 */
process sayHello {

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
params.greeting = 'greetings.csv'

workflow {

    main:
    greetings_array = ['Hello','Bonjour','Holà']

    // create a channel for inputs from a CSV file
    greeting_ch = Channel.fromPath(params.greeting)
                        .view { "Before splitCsv: $it" }
                        .splitCsv()
                        .view { "After splitCsv: $it" }
                        .map { line -> line[0] }
                        .view { "After map: $it" }

    // emit a greeting
    ch_output = sayHello(greeting_ch)

    publish:
    greetings = ch_output
}

output {
    greetings {
        path '.'
    }
}
