/*
 * FASTQC - Quality control for sequencing reads
 */
process FASTQC {
    tag "$id"
    container 'quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0'

    input:
    tuple val(id), path(reads)

    output:
    tuple val(id), path("*.html"), emit: html
    tuple val(id), path("*.zip"), emit: zip

    script:
    """
    fastqc --quiet --threads 2 ${reads}
    """
}
