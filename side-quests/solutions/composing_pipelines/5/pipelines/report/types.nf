#!/usr/bin/env nextflow

nextflow.enable.types = true

record Sample {
    id: String
    file: Path
}
