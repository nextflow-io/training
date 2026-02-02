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
    first_output = sayHello.out
}

output {
    first_output {
        path 'hello_world'
        mode 'copy'
    }
}
