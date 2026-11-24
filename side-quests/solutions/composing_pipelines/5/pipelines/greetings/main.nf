#!/usr/bin/env nextflow

nextflow.enable.types = true

include { GREETING_WORKFLOW } from './workflows/greeting'
include { TRANSFORM_WORKFLOW } from './workflows/transform'
include { Greeting ; Person } from './types'

params {
    names: Channel<Person>
}

workflow {
    main:
    names = params.names.map { person -> person.name }

    // Run the greeting workflow
    greeting = GREETING_WORKFLOW(names)

    // Run the transform workflow
    transform = TRANSFORM_WORKFLOW(greeting.timestamped)

    publish:
    greetings = greeting.greetings
    timestamped = greeting.timestamped
    upper = transform.upper
    reversed = transform.reversed
}

output {
    greetings: Channel<Greeting> {
        path 'greetings'
    }
    timestamped: Channel<Greeting> {
        path 'timestamped'
    }
    upper: Channel<Greeting> {
        path 'upper'
    }
    reversed: Channel<Greeting> {
        path 'reversed'
    }
}
