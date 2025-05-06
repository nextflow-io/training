#!/usr/bin/env nextflow

process collectGreetings {

    publishDir 'results', mode: 'copy'

    input:  // not order of parameters matters - affects input order
    path input_files
    val batch_name

    output: // outputs output in oder/positional manner - adding emit - allows to emit output under emit name
    path "COLLECTED-${batch_name}-output.txt", emit: outFile
    val count_greetings, emit: count
    // val count_greetings = input_files.size()  // this will not work - it is not a valid syntax

    script:  // script-block triple quotes is a script that will be executed in bash - can add code above block - will be executed but not in bash
    count_greetings = input_files.size()   // size() operator counts input size 
    """
    cat $input_files > COLLECTED-${batch_name}-output.txt
    """
}