#!/usr/bin/env nextflow

nextflow.enable.types = true

record Person {
    name: String
}

record Greeting {
    name: String
    file: Path
}
