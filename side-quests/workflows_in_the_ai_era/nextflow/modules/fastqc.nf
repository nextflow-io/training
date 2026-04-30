process FASTQC {
    tag "$id"
    container 'quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0'

    input:
    // TODO: Define input - a tuple with sample id and the read files
    ???

    output:
    // TODO: Define outputs - HTML reports and ZIP files
    ???

    script:
    // TODO: Add the fastqc command
    """
    ???
    """
}
