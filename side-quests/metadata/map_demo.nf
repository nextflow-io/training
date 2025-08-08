#!/usr/bin/env nextflow

// Create a simple map
def my_map = [id:'sampleA', character:'squirrel']

// Print the whole map
println "map: ${my_map}"

// Access individual values using dot notation
println "id: ${my_map.id}"
println "character: ${my_map.character}"
