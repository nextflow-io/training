#!/usr/bin/env nextflow

// Enable workflow output definition syntax (required for Nextflow < 25.10)
nextflow.preview.output = true

/*
 * Pipeline parameters
 */
params.input = 'greetings.csv'

// Include modules
include { SAY_HELLO } from './modules/greetings.nf'
include { CONVERT_TO_UPPER } from './modules/greetings.nf'

workflow {
    main:
    // Create a channel from the CSV file with metadata
    greeting_ch = channel.fromPath(params.input)
                        .splitCsv(header: true)
                        .map { row -> [[id: row.greeting, language: row.language], row.greeting] }

    // Create greeting files
    SAY_HELLO(greeting_ch)

    // Convert to uppercase
    CONVERT_TO_UPPER(SAY_HELLO.out)

    publish:
    greetings = SAY_HELLO.out
    uppercase = CONVERT_TO_UPPER.out
}

/*
 * Output block defines how published outputs are organized
 */
output {
    greetings {
        mode 'copy'
        path { meta, file -> "greetings/${meta.language}" }
        index {
            path 'greetings/index.json'
        }
    }

    uppercase {
        mode 'copy'
        path { meta, file -> "uppercase/${meta.language}" }
        index {
            path 'uppercase/index.json'
        }
    }
}
