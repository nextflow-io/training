#!/usr/bin/env nextflow

def counts = ["China": 1, "India": 2, "USA": 3]

def result = 0
counts.keySet().each { v ->
    result += counts[v]
}

println result