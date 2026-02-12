#!/usr/bin/env nextflow

/*
 * Pipeline parameters
 */
params.input_file = "greetings.csv"

/*
 * Use echo to print 'Hello World!' to standard out
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
 * Use a text replace utility to convert the greeting to uppercase
 */
process convertToUpper {

    input:
    path input_file

    output:
    path "UPPER-${input_file}"

    script:
    """
    cat '$input_file' | tr '[a-z]' '[A-Z]' > UPPER-${input_file}
    """
}

workflow {
    main:
    // create a channel for inputs from a CSV file
    greeting_ch = channel.fromPath(params.input_file).splitCsv().flatten()

    // emit a greeting
    sayHello(greeting_ch)

    // convert the greeting to uppercase
    convertToUpper(sayHello.out)

    publish:
    greetings = sayHello.out
    upper_greetings = convertToUpper.out
}

output {
    directory 'results'
    mode 'copy'
}
