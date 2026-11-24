#!/usr/bin/env nextflow

/*
 * Use echo to print a greeting to a file
 */
process sayHello {

    tag "${name}"

    input:
    val name

    output:
    path "${name}.txt"

    script:
    """
    echo 'Hello, ${name}!' > '${name}.txt'
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
    name_ch = channel.of('Mateo','Raphael','Celine')
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
