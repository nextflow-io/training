process FASTP {
    tag "$meta.id"
    container 'biocontainers/fastp:0.23.4'
    publishDir "${params.outdir}/trimmed", mode: 'copy'

    cpus 4
    memory '8.GB'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_R{1,2}.fastq.gz"), emit: reads
    path "*.{json,html}", emit: reports

    script:
    def prefix = meta.id
    """
    fastp \\
        -i ${reads[0]} \\
        -I ${reads[1]} \\
        -o ${prefix}_R1.fastq.gz \\
        -O ${prefix}_R2.fastq.gz \\
        --json ${prefix}.json \\
        --html ${prefix}.html \\
        --thread $task.cpus
    """
}
