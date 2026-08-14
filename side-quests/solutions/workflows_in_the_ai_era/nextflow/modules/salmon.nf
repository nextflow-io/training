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
    tar -xzf ${archive}
    """
}

/*
 * SALMON_QUANT - Fast transcript quantification
 */
process SALMON_QUANT {
    tag "$id"
    container 'quay.io/biocontainers/salmon:1.10.3--h6dccd9a_2'

    cpus 4
    memory '8.GB'

    input:
    tuple val(id), path(reads)
    path index

    output:
    tuple val(id), path("${id}"), emit: results

    script:
    """
    salmon quant \\
        --index ${index} \\
        --libType A \\
        --mates1 ${reads[0]} \\
        --mates2 ${reads[1]} \\
        --output ${id} \\
        --threads ${task.cpus}
    """
}
