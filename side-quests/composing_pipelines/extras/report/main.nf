#!/usr/bin/env nextflow

nextflow.enable.types = true

include { COUNT_TEXT } from './modules/count_text'
include { SUMMARIZE } from './modules/summarize'
include { Sample } from './types'

params {
    // Samplesheet with one text file per sample (columns: id, file)
    input: Channel<Sample>

    // What to count in each file: 'chars' or 'words'
    metric: String = 'chars'

    // Title printed at the top of the report
    title: String = 'Text report'
}

workflow {
    main:
    counts = COUNT_TEXT(params.input, params.metric)
    summary = SUMMARIZE(counts.collect())

    publish:
    summary = summary
}

output {
    summary: Path {
        path 'report'
    }
}
