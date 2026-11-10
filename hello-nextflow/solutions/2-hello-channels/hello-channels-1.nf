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
    name: String = 'Friend'
}

workflow {

    main:
    // create a channel for inputs
    name_ch = channel.of('Mateo')
                        .view()
    // emit a greeting
    sayHello(name_ch)

    publish:
    greetings = sayHello.out
}

output {
    greetings {
        path 'hello'
        mode 'copy'
    }
}
