#!/usr/bin/env nextflow

params.greeting = 'Hello world!'
greeting_ch = Channel.from(params.greeting)

include { SPLITLETTERS as SPLITLETTERS_one } from './modules.nf'
include { SPLITLETTERS as SPLITLETTERS_two } from './modules.nf'

include { CONVERTTOUPPER as CONVERTTOUPPER_one } from './modules.nf'
include { CONVERTTOUPPER as CONVERTTOUPPER_two } from './modules.nf'

workflow {
    letters_ch1 = SPLITLETTERS_one(greeting_ch)
    results_ch1 = CONVERTTOUPPER_one(letters_ch1.flatten())
    results_ch1.view{ it }

    letters_ch2 = SPLITLETTERS_two(greeting_ch)
    results_ch2 = CONVERTTOUPPER_two(letters_ch2.flatten())
    results_ch2.view{ it }
}