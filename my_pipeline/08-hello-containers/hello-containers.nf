#!/usr/bin/env nextflow

include { sayHello } from './modules/sayHello.nf'
include { convertToUppercase } from './modules/convertToUppercase.nf'
include { collectGreetings } from './modules/collectGreetings.nf'


/* * Pipeline parameters */
params.greeting = 'greetings.csv'
params.batch_name = 'test-batch'

workflow {
    // create a channel from the CSV file
    greeting_ch = channel.fromPath(params.greeting)
                         .view{ csv -> "Before splitCsv: $csv" }
                         .splitCsv()
                         .view{ csv -> "After splitCsv: $csv" }
                         .map{ item -> item[0] }
                         .view{ csv -> "After map: $csv" }

    // emit a greeting
    // nf created the channels automatically
    sayHello(greeting_ch)

    // Convert the output of sayHello to uppercase
    convertToUppercase(sayHello.out)

    // Collect all greetings into a single file
    collectGreetings(sayHello.out.collect(), params.batch_name)
    convertToUppercase.out.view{ greetings -> "Before collectGreetings: $greetings" }
    convertToUppercase.out.collect()
        .view{ greetings -> "After collectGreetings: $greetings" }

    collectGreetings.out.count.view{ num_greetings -> "There were $num_greetings greetings!" }
}