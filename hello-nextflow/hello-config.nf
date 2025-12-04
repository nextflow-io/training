#!/usr/bin/env nextflow

// Include modules
include { sayHello } from './modules/sayHello.nf'
include { convertToUpper } from './modules/convertToUpper.nf'
include { collectGreetings } from './modules/collectGreetings.nf'
include { cowpy } from './modules/cowpy.nf'

/*
 * Pipeline parameters
 */
params {
    greeting: Path = 'greetings.csv'
    batch: String = 'test-batch'
    character: String = 'turkey'
}

workflow {

    // create a channel for inputs from a CSV file
    greeting_ch = channel.fromPath(params.greeting)
                        .splitCsv()
                        .map { line -> line[0] }

    // emit a greeting
    sayHello(greeting_ch)

    // convert the greeting to uppercase
    convertToUpper(sayHello.out)

    // collect all the greetings into one file
    collectGreetings(convertToUpper.out.collect(), params.batch)

    // emit a message about the size of the batch
    collectGreetings.out.count.view { num_greetings -> "There were $num_greetings greetings in this batch" }

    // generate ASCII art of the greetings with cowpy
    cowpy(collectGreetings.out.outfile, params.character)
}
