#!/usr/bin/env nextflow

/*
 * Use echo to print 'Hello World!' to a file
 */
process sayHello {
    
    publishDir 'results', mode: 'copy'

    input:
    val greeting

    output:
        path "$greeting-output.txt"

    script:
    """
    echo '$greeting' > '$greeting-output.txt'
    """
}

/* * Passing default parameters to the workflow */
params.greeting = 'Hello World!'

workflow {
    greetings_array = [params.greeting, 'Hello, Nextflow!', 'Greetings, NF!']
    greeting_ch = channel.of(greetings_array)
                         .view{ greeting -> "Before flatten: $greeting" }
                         .flatten()
                         .view{ greeting -> "After flatten: $greeting" }

    // emit a greeting
    // nf created the channels automatically
    sayHello(greeting_ch)
}