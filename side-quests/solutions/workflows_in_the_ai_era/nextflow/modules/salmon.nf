/*
 * UNTAR - Extract tarball archive
 *
 * Helper process to extract the pre-built salmon index.
 */
process UNTAR {
    container 'ubuntu:22.04'

    input:
    path archive

    output:
    path "salmon", emit: index

    script:
    """
    tar -xzf $archive
    """
}

/*
 * SALMON - Fast transcript quantification
 *
 * Salmon uses pseudo-alignment for rapid quantification.
 * We use a pre-built index to keep the tutorial fast.
 */
process SALMON_QUANT {
    tag "$meta.id"
    container 'quay.io/biocontainers/salmon:1.10.3--h6dccd9a_2'
    publishDir "${params.outdir}/salmon", mode: 'copy'

    cpus 4
    memory '8.GB'

    input:
    tuple val(meta), path(reads)
    path index

    output:
    tuple val(meta), path("${meta.id}"), emit: results

    script:
    """
    salmon quant \\
        --index $index \\
        --libType A \\
        --mates1 ${reads[0]} \\
        --mates2 ${reads[1]} \\
        --output ${meta.id} \\
        --threads $task.cpus
    """
}
