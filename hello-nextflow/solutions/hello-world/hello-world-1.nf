#!/usr/bin/env nextflow

process sayHello {

    output:
        stdout

    """
    echo 'Hello World!'
    """
}

workflow {
    sayHello()
}
