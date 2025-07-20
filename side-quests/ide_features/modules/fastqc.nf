process FASTQC {
    tag "${sample_id}"
    publishDir "${params.output_dir}/fastqc", mode: 'copy'

    // Container directive for reproducibility
    container 'biocontainers/fastqc:v0.11.9_cv8'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("*.html"), emit: html
    tuple val(sample_id), path("*.zip"), emit: zip

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    fastqc \\
        ${args} \\
        --threads ${task.cpus} \\
        ${reads}
    """
}
