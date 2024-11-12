#!/usr/bin/env nextflow

/*
 * Use echo to print 'Hello World!' to standard out
 */
process sayHello {

    output:
        stdout

    script:
    """
    echo 'Hello World!'
    """
}

workflow {

    // emit a greeting
    sayHello()
}
