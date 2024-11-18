#!/usr/bin/env nextflow

/*
 * Pipeline parameters
 */

params.input_file = "data/greetings.csv"
process sayHello {
    publishDir 'results', mode: 'copy'

    input:
    val greeting

    output:
    path "${greeting}-output.txt"

    script:
    """
    echo '$greeting' > '$greeting-output.txt'
    """
}

process convertToUpper {
    publishDir 'results', mode: 'copy'

    input:
    path input_file

    output:
    path "UPPER-${input_file}"

    script:
    """
    cat '$input_file' | tr '[a-z]' '[A-Z]' > "UPPER-${input_file}"
    """
}

workflow {
    // Create a channel for inputs from a CSV file
    greeting_ch = Channel.fromPath(params.input_file)
                         .splitCsv()
                         .flatten()


    // Emit a greeting
    sayHello(greeting_ch)

    // Convert the greeting to uppercase
    convertToUpper(sayHello.out)
}

