#!/usr/bin/env nextflow

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
        val batch_name

    output:
        path "COLLECTED-${batch_name}-output.txt" , emit: outfile
        val count_greetings , emit: count

    script:
        count_greetings = input_files.size()
    """
    cat ${input_files} > 'COLLECTED-${batch_name}-output.txt'
    """
}

/*
 * Pipeline parameters
 */
params.greeting = 'greetings.csv'
params.batch = 'test-batch'

// Include modules
include { sayHello } from './modules/sayHello.nf'

workflow {

    main:
    // create a channel for inputs from a CSV file
    greeting_ch = Channel.fromPath(params.greeting)
                        .splitCsv()
                        .map { line -> line[0] }

    // emit a greeting
    ch_hello = sayHello(greeting_ch)

    // convert the greeting to uppercase
    ch_upper = convertToUpper(ch_hello)

    // collect all the greetings into one file
    ch_collected = collectGreetings(ch_upper.collect(), params.batch)

    // emit a message about the size of the batch
    ch_collected.count.view { "There were $it greetings in this batch" }

    publish:
    greetings = ch_hello
    uppercase = ch_upper
    collected = ch_collected.outfile
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
