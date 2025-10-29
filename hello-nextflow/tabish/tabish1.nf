#!/usr/bin/env nextflow

/*
 * Me playing around with stuff
 * Part 1
 */

process printGreeting {
    publishDir 'results', mode: 'copy'

    input:
    val greeting

    output:
        path 'tabish1.txt'

    script:
    """
    echo '$greeting' > 'tabish1.txt'
    """
}

process printName {

    publishDir 'results', mode: 'copy'

    input:
    val name

    output:
        path 'tabish1.txt'

    script:
    """
    echo '$name' > tabish1.txt
    """

}

params.greeting = "Good Morning"
params.name = "Tabish"

workflow {
    printGreeting(params.greeting)
    printName(params.name)
}
