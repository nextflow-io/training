#!/usr/bin/env nextflow

/*
 * Use echo to print 'Hello World!' to a file
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
    name: String = 'World'
}

workflow {

    main:
    // create a channel for inputs
    name_ch = channel.of('Mateo')
                        .view()
    // emit a greeting
    sayHello(name_ch)

    publish:
    first_output = sayHello.out
}

output {
    first_output {
        path '.'
    }
}
