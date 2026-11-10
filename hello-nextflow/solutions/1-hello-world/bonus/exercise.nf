#!/usr/bin/env nextflow

/*
 * Use echo to print a greeting to a file
 */
process sayHello {

    input:
    val name

    output:
    path 'greeting.txt'

    script:
    """
    echo 'Hello, ${name}!' > greeting.txt
    """
}

/*
 * Pipeline parameters
 */
params {
    name: String = 'Friend'
}

workflow {

    main:
    // emit a greeting
    sayHello(params.name)

    publish:
    greetings = sayHello.out
}

output {
    greetings {
        path 'hello'
        mode 'copy'
    }
}
