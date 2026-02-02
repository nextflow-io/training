/*
 * FASTP - Adapter trimming and quality filtering
 *
 * Different container than FastQC - each process has its own
 * isolated software environment. No version conflicts!
 */
process FASTP {
    tag "$meta.id"
    container 'quay.io/biocontainers/fastp:0.23.4--hadf994f_2'
    publishDir "${params.outdir}/fastp", mode: 'copy'

    cpus 4
    memory '4.GB'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_trimmed*.fastq.gz"), emit: reads
    tuple val(meta), path("*.json"), emit: json
    tuple val(meta), path("*.html"), emit: html

    script:
    def prefix = meta.id
    """
    fastp \\
        --in1 ${reads[0]} \\
        --in2 ${reads[1]} \\
        --out1 ${prefix}_trimmed_R1.fastq.gz \\
        --out2 ${prefix}_trimmed_R2.fastq.gz \\
        --json ${prefix}.fastp.json \\
        --html ${prefix}.fastp.html \\
        --thread $task.cpus
    """
}
