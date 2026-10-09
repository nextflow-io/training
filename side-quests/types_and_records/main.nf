#!/usr/bin/env nextflow

include { SAY_HELLO } from './modules/say_hello.nf'
include { SHOUT } from './modules/shout.nf'
include { COLLECT_GREETINGS } from './modules/collect_greetings.nf'

params {
    input: Path = 'data/people.csv'
    batch: String = 'batch'
}

workflow {
    main:
    // Read the CSV into a channel of meta maps, one per person
    people = channel.fromPath(params.input)
        .splitCsv(header: true)
        .map { row ->
            [id: row.id, name: row.name, greeting: row.greeting]
        }

    greetings = SAY_HELLO(people)
    shouted = SHOUT(greetings)
    collected = COLLECT_GREETINGS(shouted.map { _meta, file -> file }.collect(), params.batch)

    publish:
    greetings = greetings
    shouted = shouted
    collected = collected
}

output {
    greetings {
        path 'greetings'
    }
    shouted {
        path 'shouted'
    }
    collected {
    }
}
