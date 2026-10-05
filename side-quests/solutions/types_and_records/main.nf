#!/usr/bin/env nextflow

nextflow.enable.types = true

include { SAY_HELLO } from './modules/say_hello.nf'
include { SHOUT } from './modules/shout.nf'
include { COLLECT_GREETINGS } from './modules/collect_greetings.nf'

params {
    input: Path = 'data/people.csv'
    batch: String = 'batch'
}

record Person {
    id: String
    name: String
    language: String
    greeting: String
}

workflow GREET {
    take:
    people: Channel<Person>
    batch: Value<String>

    main:
    greetings = SAY_HELLO(people)
    shouted = SHOUT(greetings)
    collected = COLLECT_GREETINGS(shouted.map { r -> r.shouted }.collect(), batch)

    emit:
    greetings: Channel<Record> = greetings
    shouted: Channel<Record> = shouted
    collected: Value<Path> = collected
}

workflow {
    main:
    // Read the CSV into a channel of records, one per person
    people = channel.of(params.input)
        .flatMap { csv -> csv.splitCsv(header: true) }
        .map { row ->
            record(id: row.id, name: row.name, language: row.language, greeting: row.greeting)
        }
        .map { person ->
            person + record(name: person.name.toLowerCase().capitalize())
        }

    res = GREET(people, channel.value(params.batch))

    publish:
    greetings = res.greetings
    shouted = res.shouted
    collected = res.collected
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
