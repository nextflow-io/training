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
    echo '${greeting}' > '${greeting}-output.txt'
    """
}

/*
 * Pipeline parameters
 */
params {
    input: String = 'Holà mundo!'
}

workflow {

    main:
    // declare an array of input greetings
    greetings_array = ['Hello','Bonjour','Holà']
    // create a channel for inputs
    greeting_ch = channel.of(greetings_array)
                        .view { greeting -> "Before flatten: $greeting" }
                        .flatten()
                        .view { greeting -> "After flatten: $greeting" }
    // emit a greeting
    sayHello(greeting_ch)

    publish:
    first_output = sayHello.out
}

output {
    first_output {
        path 'hello_channels'
        mode 'copy'
    }
}
