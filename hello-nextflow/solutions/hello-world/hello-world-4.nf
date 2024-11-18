#!/usr/bin/env nextflow

/*
 * Use echo to print 'Hello World!' to standard out
 */
process sayHello {

    publishDir 'results', mode: 'copy'

    output:
        path 'output.txt'

    script:
    """
    echo 'Hello World!' > output.txt
    """
}

workflow {

    // emit a greeting
    sayHello()
}
