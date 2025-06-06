#!/usr/bin/env nextflow
/*
 * Pipeline parameters
 */
params.greeting = 'greetings.csv'

workflow {

    greetings_array = ['Hello', 'Bonjour', 'Holà']

    // create a channel for inputs from a CSV file
    greeting_ch = Channel
        .fromPath(params.greeting)
        .view { "Before splitCsv: ${it}" }
        .splitCsv()
        .view { "After splitCsv: ${it}" }
        .map { line -> line[0] }
        .view { "After map: ${it}" }

    // emit a greeting
    sayHello(greeting_ch)
}


/*
 * Use echo to print 'Hello World!' to a file
 */
process sayHello {

    publishDir 'results', mode: 'copy'

    input:
    val greeting

    output:
    path "${greeting}-output.txt"

    script:
    """
    echo '${greeting}' > '${greeting}-output.txt'
    """
}
