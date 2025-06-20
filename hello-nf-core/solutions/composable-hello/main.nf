#!/usr/bin/env nextflow

// import the workflow code from the hello.nf file
include { HELLO } from './hello.nf'

// declare input parameter
params.greeting = 'greetings.csv'

workflow {
    // create a channel for inputs from a CSV file
    greeting_ch = Channel.fromPath(params.greeting)
                        .splitCsv()
                        .map { line -> line[0] }

    // call the imported workflow on the channel of greetings
    HELLO(greeting_ch)

    // view the outputs emitted by the workflow
    HELLO.out.view { "Outputs: $it" }
}
