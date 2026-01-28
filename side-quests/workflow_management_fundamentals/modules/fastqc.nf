process FASTQC {
    tag "$meta.id"
    container 'biocontainers/fastqc:0.12.1'
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
    fastqc -q -t $task.cpus ${reads}
    """
}
