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
    greeting = sayHello()

    publish:
    greetings = greeting
}

output {
    greetings {
        path 'hello'
        mode 'copy'
    }
}
