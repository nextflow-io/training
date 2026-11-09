#!/usr/bin/env nextflow

nextflow.enable.types = true

include { SAY_HELLO_UPPER } from '../modules/say_hello_upper'
include { REVERSE_TEXT } from '../modules/reverse_text'
include { Greeting } from '../types'

workflow TRANSFORM_WORKFLOW {
    take:
    greetings: Channel<Greeting>

    main:
    // Apply transformations in sequence
    upper_ch = SAY_HELLO_UPPER(greetings)
    reversed_ch = REVERSE_TEXT(upper_ch)

    emit:
    upper: Channel<Greeting> = upper_ch
    reversed: Channel<Greeting> = reversed_ch
}
