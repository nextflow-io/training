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
params.greeting = 'Holà mundo!'

workflow {

    main:
    // create a channel for inputs
    greeting_ch = Channel.of('Hello','Bonjour','Holà')

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
