#!/usr/bin/env nextflow

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

workflow {

    main:
    // emit a greeting
    sayHello()

    publish:
    greetings = sayHello.out
}

output {
    greetings {
        path '.'
    }
}
