#!/usr/bin/env nextflow

// import processes from modules
include { sayHello } from './modules/sayHello.nf' 
include { convertToUpper } from './modules/convertToUpper.nf'
include { collectGreetings } from './modules/collectGreetings.nf'

/*
 * Pipeline parameters
 */
params.greeting = 'greetings.csv'
params.batch_name = 'test-batch'

workflow {

    // create a channel for inputs from a CSV file
    greeting_ch = Channel.fromPath(params.greeting)
                        .splitCsv()
                        .map { line -> line[0] }

    // emit a greeting
    sayHello(greeting_ch)

    convertToUpper(sayHello.out)

    collectGreetings(convertToUpper.out.collect(), params.batch_name)
    // convertToUpper.out.view { greetings -> "before collect: $greetings"}
    // convertToUpper.out.collect().view { greetings -> "after collect: $greetings" }

    collectGreetings.out.count.view { num_greetings -> "there were $num_greetings greetings" }
 
}
