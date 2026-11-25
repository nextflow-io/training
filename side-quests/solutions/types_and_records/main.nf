#!/usr/bin/env nextflow

nextflow.enable.types = true

include { SAY_HELLO } from './modules/say_hello.nf'
include { SHOUT } from './modules/shout.nf'
include { COLLECT_GREETINGS } from './modules/collect_greetings.nf'
include { Person ; Greeting ; Shouted } from './types.nf'

params {
    input: Channel<Person>
    batch: String = 'batch'
}

workflow GREET {
    take:
    people: Channel<Person>
    batch: Value<String>

    main:
    greetings = SAY_HELLO(people)
    shouted = SHOUT(greetings)
    collected = COLLECT_GREETINGS(shouted.map { s -> s.shouted }.collect(), batch)

    emit:
    greetings: Channel<Greeting> = greetings
    shouted: Channel<Shouted> = shouted
    collected: Value<Path> = collected
}

workflow {
    main:
    res = GREET(params.input, channel.value(params.batch))

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
