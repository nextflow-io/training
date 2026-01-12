/*
 * FASTQC - Quality control for sequencing reads
 *
 * Each process declares its own container - no manual installation needed.
 * The exact version (0.12.1) is locked, ensuring reproducible results.
 */
process FASTQC {
    tag "$meta.id"
    container 'biocontainers/fastqc:0.12.1--hdfd78af_0'
    publishDir "${params.outdir}/fastqc", mode: 'copy'

    cpus 2
    memory '2.GB'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.html"), emit: html
    tuple val(meta), path("*.zip"), emit: zip

    script:
    """
    fastqc --quiet --threads $task.cpus ${reads}
    """
}
