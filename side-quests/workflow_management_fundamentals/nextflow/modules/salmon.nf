// Helper process to extract the pre-built salmon index
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

process SALMON_QUANT {
    tag "$meta.id"
    container 'quay.io/biocontainers/salmon:1.10.3--h6dccd9a_2'
    publishDir "${params.outdir}/salmon", mode: 'copy'

    cpus 4
    memory '8.GB'

    input:
    // TODO: Define inputs - trimmed reads AND the salmon index
    // Hint: This process needs two inputs (one tuple, one path)
    ???

    output:
    // TODO: Define output - the quantification results directory
    ???

    script:
    // TODO: Add the salmon quant command
    // Hint: Use $index for the index path
    // Hint: Use ${reads[0]} and ${reads[1]} for input reads
    // Hint: Use $task.cpus for thread count
    """
    ???
    """
}
