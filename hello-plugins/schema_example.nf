#!/usr/bin/env nextflow

params.input = 'greetings.csv'

workflow {
    Channel.fromPath(params.input)
        .splitCsv(header: true)
        .map { row -> row.greeting }
        .view { greeting -> "Greeting: $greeting" }
}
