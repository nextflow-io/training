#!/usr/bin/env nextflow

/*
 * Use echo to print 'Hello World!' to a file
 */
process sayHello {

    input:
        val greeting

    output:
        path 'output.txt'

    script:
    """
    echo '$greeting' > output.txt
    """
}

/*
 * Pipeline parameters
 */
params.greeting = 'Holà mundo!'

workflow {

    // create a channel for inputs
    greeting_ch = Channel.of('Hello Channels!')

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
