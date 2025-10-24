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
 * Use a text replacement tool to convert the greeting to uppercase
 */
process convertToUpper {

    input:
        path input_file

    output:
        path "UPPER-${input_file}"

    script:
    """
    cat '$input_file' | tr '[a-z]' '[A-Z]' > 'UPPER-${input_file}'
    """
}

/*
 * Collect uppercase greetings into a single output file
 */
process collectGreetings {

    input:
        path input_files

    output:
        path "COLLECTED-output.txt"

    script:
    """
    cat ${input_files} > 'COLLECTED-output.txt'
    """
}

/*
 * Pipeline parameters
 */
params.greeting = 'greetings.csv'

workflow {

    // create a channel for inputs from a CSV file
    greeting_ch = Channel.fromPath(params.greeting)
                        .splitCsv()
                        .map { line -> line[0] }

    // emit a greeting
    ch_hello = sayHello(greeting_ch)

    // convert the greeting to uppercase
    ch_upper = convertToUpper(ch_hello)

    // collect all the greetings into one file
    ch_collected = collectGreetings(ch_upper.collect())

    // optional view statements
    ch_upper.view { "Before collect: $it" }
    ch_upper.collect().view { "After collect: $it" }

    publish:
    greetings = ch_hello
    uppercase = ch_upper
    collected = ch_collected
}

output {
    greetings {
        path '.'
    }
    uppercase {
        path '.'
    }
    collected {
        path '.'
    }
}
