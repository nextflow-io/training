#!/usr/bin/env nextflow

nextflow.enable.types = true

include { VALIDATE_NAME } from '../modules/validate_name'
include { SAY_HELLO } from '../modules/say_hello'
include { TIMESTAMP_GREETING } from '../modules/timestamp_greeting'
include { Greeting } from '../types'

workflow GREETING_WORKFLOW {
    take:
    names: Channel<String>

    main:
    // Chain processes: validate -> create greeting -> add timestamp
    validated_ch = VALIDATE_NAME(names)
    greetings_ch = SAY_HELLO(validated_ch)
    timestamped_ch = TIMESTAMP_GREETING(greetings_ch)

    emit:
    greetings: Channel<Greeting> = greetings_ch
    timestamped: Channel<Greeting> = timestamped_ch
}
