#!/usr/bin/env nextflow

/*
 * Pipeline parameters
 */
params.greeting = 'greetings.csv'
params.batch = 'test-batch'

// Include modules
include { sayHello } from './modules/sayHello.nf'
include { convertToUpper } from './modules/convertToUpper.nf'
include { collectGreetings } from './modules/collectGreetings.nf'

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
