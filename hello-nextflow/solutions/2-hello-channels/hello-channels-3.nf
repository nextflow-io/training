#!/usr/bin/env nextflow

/*
 * Use echo to print 'Hello World!' to a file
 */
process sayHello {

    input:
        val greeting

    output:
        path "${greeting}-output.txt"

    script:
    """
    echo '$greeting' > '$greeting-output.txt'
    """
}

/*
 * Pipeline parameters
 */
params.greeting = 'Holà mundo'

workflow {

    main:
    greetings_array = ['Hello','Bonjour','Holà']

    // create a channel for inputs
    greeting_ch = Channel.of(greetings_array)
                    .view { greeting -> "Before flatten: $greeting" }
                    .flatten()
                    .view { greeting -> "After flatten: $greeting" }

    // emit a greeting
    sayHello(greeting_ch)

    publish:
    greetings = sayHello.out
}

output {
    greetings {
        path '.'
    }
}
