#!/usr/bin/env nextflow
process collectGreetings {
    
    publishDir 'results', mode: 'copy'

    input:
    path input_file
    val batch_name 

    output:
        path "COLLECTED-${batch_name}-output.txt", emit: outfile
        val count_greetings, emit: count

    script:
    count_greetings=input_file.size()
    """
    cat $input_file > "COLLECTED-${batch_name}-output.txt"
    """
}
