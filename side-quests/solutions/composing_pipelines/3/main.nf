#!/usr/bin/env nextflow

nextflow.enable.types = true

include { workflow as GREETINGS } from './pipelines/greetings'
include { Person } from './pipelines/greetings/types'

params {
    names: Channel<Person>
}

workflow {
    main:
    greetings = GREETINGS(record(names: params.names))

    publish:
    timestamped = greetings.timestamped
    reversed = greetings.reversed
}

output {
    timestamped {
        path 'timestamped'
    }
    reversed {
        path 'reversed'
    }
}
