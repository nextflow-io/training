#!/usr/bin/env nextflow

/*
 * Example: Local function vs plugin function
 * Edit this file to replace the local function with the nf-hello plugin
 */

// Local function - defined in this file
def sayHello(name) {
    return "Hello, ${name}!"
}

workflow {
    Channel.of('Alice', 'Bob', 'Carol')
        | map { name -> sayHello(name) }
        | view
}
