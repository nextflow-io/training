process QUAST {
    tag "$meta.id"
    container 'biocontainers/quast:5.2.0'
    publishDir "${params.outdir}/quast", mode: 'copy'

    cpus 2
    memory '4.GB'

    input:
    tuple val(meta), path(assembly)

    output:
    tuple val(meta), path("${meta.id}/*"), emit: report

    script:
    """
    quast.py \\
        $assembly \\
        -o ${meta.id} \\
        --threads $task.cpus
    """
}
