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
    // create a channel for inputs from a CSV file
    greeting_ch = Channel.fromPath(params.greeting)
                        .view { csv -> "Before splitCsv: $csv" }
                        .splitCsv()
                        .view { csv -> "After splitCsv: $csv" }
                        .map { line -> line[0] }
                        .view { csv -> "After map: $csv" }

    // emit a greeting
    sayHello(greeting_ch)

    publish:
    greetings = sayHello.out
}

output {
    greetings {
        path '.'
    }
}
