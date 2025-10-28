#!/usr/bin/env nextflow

/*
 * Pipeline parameters
 */
params.input = 'greetings.csv'

// Include modules
include { sayHello } from './modules/sayHello.nf'
include { convertToUpper } from './modules/convertToUpper.nf'
include { collectGreetings } from './modules/collectGreetings.nf'

workflow {

    // create a channel for inputs from a CSV file
    greeting_ch = channel.fromPath(params.input)
                        .splitCsv()
                        .map { line -> line[0] }

    // emit a greeting
    sayHello(greeting_ch)

    // convert the greeting to uppercase
    convertToUpper(sayHello.out)

    // collect all the greetings into one file
    collectGreetings(convertToUpper.out.collect())

    // emit a message about the size of the batch
    sayHello.out.count().view { "There were $it greetings in this batch" }
}
