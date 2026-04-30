/*
 * FASTP - Adapter trimming and quality filtering
 */
process FASTP {
    tag "$id"
    container 'quay.io/biocontainers/fastp:0.23.4--hadf994f_2'

    cpus 4
    memory '4.GB'

    input:
    tuple val(id), path(reads)

    output:
    tuple val(id), path("*_trimmed*.fastq.gz"), emit: reads
    tuple val(id), path("*.json"), emit: json
    tuple val(id), path("*.html"), emit: html

    script:
    """
    fastp \\
        --in1 ${reads[0]} \\
        --in2 ${reads[1]} \\
        --out1 ${id}_trimmed_R1.fastq.gz \\
        --out2 ${id}_trimmed_R2.fastq.gz \\
        --json ${id}.fastp.json \\
        --html ${id}.fastp.html \\
        --thread ${task.cpus}
    """
}
