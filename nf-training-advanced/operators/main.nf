#!/usr/bin/env nextflow

workflow {
    Channel.of( 1, 2, 3, 4, 5 )
        .map { num -> num * num }
        .view()
}
