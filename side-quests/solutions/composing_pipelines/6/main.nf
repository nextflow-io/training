#!/usr/bin/env nextflow

nextflow.enable.types = true

include { params as GreetingsParams ; workflow as GREETINGS } from './pipelines/greetings'
include { params as ReportParams ; workflow as REPORT } from './pipelines/report'

params {
    greetings: GreetingsParams
    report: ReportParams
}

workflow {
    main:
    greetings = GREETINGS(params.greetings)

    // Reshape each greeting into the sample record the report expects
    samples = greetings.reversed.map { greeting -> record(id: greeting.name, file: greeting.file) }
    summary = REPORT(params.report + record(input: samples))

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
