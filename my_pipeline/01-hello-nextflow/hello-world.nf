#!/usr/bin/env nextflow

/*
 * Use echo to print 'Hello World!' to a file
 */
process sayHello {
    
    publishDir 'results', mode: 'copy'

    input:
    val greeting

    output:
        path 'output.txt'

    script:
    """
    echo '$greeting' > output.txt
    """
}

/* * Passing default parameters to the workflow */
params.greeting = 'Hello World!'

workflow {

    // emit a greeting
    sayHello(params.greeting)
}