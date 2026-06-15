#!/usr/bin/env nextflow

process SAYHELLO {
    input:
    val greeting

    output:
    path "${greeting}-output.txt"

    script:
    """
    echo '${greeting}' > '${greeting}-output.txt'
    """
}

process CONVERTTOUPPER {
    input:
    path input_file

    output:
    path "UPPER-${input_file}"

    script:
    """
    cat '${input_file}' | tr '[a-z]' '[A-Z]' > 'UPPER-${input_file}'
    """
}

params {
    input: Path = 'data/greetings.csv'
}

workflow {
    main:
    greeting_ch = channel.fromPath(params.input)
        .splitCsv()
        .map { line -> line[0] }

    SAYHELLO(greeting_ch)
    CONVERTTOUPPER(SAYHELLO.out)

    publish:
    uppercased = CONVERTTOUPPER.out
}

output {
    uppercased {
        path 'results'
        mode 'copy'
    }
}
