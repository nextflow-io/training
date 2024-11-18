#!/usr/bin/env nextflow

process sayHello {

    output:
        stdout

    script:
    """
    echo 'Hello World!'
    """
}

workflow {
    sayHello()
}
