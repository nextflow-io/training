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
    name_ch = channel.of('Ziggy')
                        .view()
    // emit a greeting
    greetings_ch = sayHello(name_ch)

    publish:
    greetings = greetings_ch
}

output {
    greetings {
        path 'hello'
        mode 'copy'
    }
}
