#!/usr/bin/env nextflow
workflow {

    // emit a greeting
    sayHello()
}


/*
 * Use echo to print 'Hello World!' to a file
 */
process sayHello {
    output:
    path 'output.txt'

    script:
    """
    echo 'Hello World!' > output.txt
    """
}
