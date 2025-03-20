#!/usr/bin/env nextflow
/*
 * Pipeline parameters
 */
params.greeting = 'Holà mundo!'

workflow {

    // emit a greeting
    sayHello(params.greeting)
}


/*
 * Use echo to print 'Hello World!' to a file
 */
process sayHello {

    publishDir 'results', mode: 'copy'

    input:
    val greeting

    output:
    path 'output.txt'

    script:
    """
    echo '${greeting}' > output.txt
    """
}
