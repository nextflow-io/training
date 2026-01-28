#!/usr/bin/env nextflow

/*
 * Use echo to print 'Hello World!' to a file
 */
process sayHello {

    input:
    val greeting

    output:
    path 'output.txt'

    script:
    """
    echo '${greeting}' > output.txt
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
    // emit a greeting
    sayHello(params.input)

    publish:
    first_output = sayHello.out
}

output {
    first_output {
        path 'hello_world'
        mode 'copy'
    }
}
