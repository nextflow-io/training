#!/usr/bin/env nextflow

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
                        .map { row -> [row.greeting, row.language] }

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
        path { greeting, language, file -> "greetings/${language}" }
        index {
            path 'greetings/index.csv'
        }
    }

    uppercase {
        mode 'copy'
        path { greeting, language, file -> "uppercase/${language}" }
        index {
            path 'uppercase/index.csv'
        }
    }
}
