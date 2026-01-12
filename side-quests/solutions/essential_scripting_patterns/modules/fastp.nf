process FASTP {
    container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

    cpus { meta.depth > 40000000 ? 4 : 2 }
    memory { 1.GB * task.attempt }
    errorStrategy 'retry'
    maxRetries 2

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_trimmed*.fastq.gz"), emit: reads
    path "*.{json,html}"                             , emit: reports

    script:
    // Simple single-end vs paired-end detection
    def is_single = reads instanceof List ? reads.size() == 1 : true

    if (is_single) {
        def input_file = reads instanceof List ? reads[0] : reads
        """
        fastp \\
            --in1 ${input_file} \\
            --out1 ${meta.id}_trimmed.fastq.gz \\
            --json ${meta.id}.fastp.json \\
            --html ${meta.id}.fastp.html \\
            --thread $task.cpus
        """
    } else {
        """
        fastp \\
            --in1 ${reads[0]} \\
            --in2 ${reads[1]} \\
            --out1 ${meta.id}_trimmed_R1.fastq.gz \\
            --out2 ${meta.id}_trimmed_R2.fastq.gz \\
            --json ${meta.id}.fastp.json \\
            --html ${meta.id}.fastp.html \\
            --thread $task.cpus
        """
    }
}
