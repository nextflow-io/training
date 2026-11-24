#!/usr/bin/env nextflow

nextflow.enable.types = true

record Greeting {
    name: String
    file: Path
}
