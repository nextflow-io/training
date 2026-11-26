#!/usr/bin/env nextflow

/*
 * Use echo to print a greeting to a file
 */
process sayHello {

    input:
    val name

    output:
    path 'output.txt'

    script:
    """
    echo 'Hello, ${name}!' > output.txt
    """
}

/*
 * Pipeline parameters
 */
params {
    who: String = 'Friend'
}

workflow {

    main:
    // emit a greeting
    greeting = sayHello(params.who)

    publish:
    greetings = greeting
}

output {
    greetings {
        path 'hello'
        mode 'copy'
    }
}
