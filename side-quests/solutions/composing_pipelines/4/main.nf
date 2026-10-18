#!/usr/bin/env nextflow

nextflow.enable.types = true

include { workflow as GREETINGS } from './pipelines/greetings'
include { workflow as REPORT } from './pipelines/report'
include { Person } from './pipelines/greetings/types'

params {
    names: Channel<Person>
}

workflow {
    main:
    greetings = GREETINGS(record(names: params.names))

    // Reshape each greeting into the sample record the report expects
    samples = greetings.reversed.map { greeting -> record(id: greeting.name, file: greeting.file) }
    summary = REPORT(record(input: samples))

    publish:
    timestamped = greetings.timestamped
    reversed = greetings.reversed
    summary = summary
}

output {
    timestamped {
        path 'timestamped'
    }
    reversed {
        path 'reversed'
    }
    summary {
        path 'report'
    }
}
