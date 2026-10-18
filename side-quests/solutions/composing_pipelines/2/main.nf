#!/usr/bin/env nextflow

nextflow.enable.types = true

include { GREETING_WORKFLOW } from './workflows/greeting'
include { TRANSFORM_WORKFLOW } from './workflows/transform'

workflow {
    main:
    names = channel.of('Alice', 'Bob', 'Charlie')

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
    greetings {
        path 'greetings'
    }
    timestamped {
        path 'timestamped'
    }
    upper {
        path 'upper'
    }
    reversed {
        path 'reversed'
    }
}
