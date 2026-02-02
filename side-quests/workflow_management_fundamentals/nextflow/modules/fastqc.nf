process FASTQC {
    tag "$meta.id"
    container 'quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0'
    publishDir "${params.outdir}/fastqc", mode: 'copy'

    input:
    // TODO: Define input - a tuple with sample metadata and read files
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
