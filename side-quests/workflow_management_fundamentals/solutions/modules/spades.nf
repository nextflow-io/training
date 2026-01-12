process SPADES {
    tag "$meta.id"
    container 'biocontainers/spades:3.15.5'
    publishDir "${params.outdir}/assemblies", mode: 'copy'

    cpus 8
    memory '16.GB'
    time '6.h'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}/contigs.fasta"), emit: assembly

    script:
    """
    spades.py \\
        -1 ${reads[0]} \\
        -2 ${reads[1]} \\
        -o ${meta.id} \\
        --threads $task.cpus \\
        --memory ${task.memory.toGiga()}
    """
}
