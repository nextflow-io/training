process TRIMGALORE {
    container 'quay.io/biocontainers/trim-galore:0.6.10--hdfd78af_0'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_trimmed*.fq"), emit: reads
    path "*_trimming_report.txt"                   , emit: reports

    script:
    // Simple single-end vs paired-end detection
    def is_single = reads instanceof List ? reads.size() == 1 : true

    if (is_single) {
        def input_file = reads instanceof List ? reads[0] : reads
        """
        trim_galore \\
            --cores $task.cpus \\
            ${input_file}

        # Rename output to match expected pattern
        mv *_trimmed.fq ${meta.id}_trimmed.fq
        """
    } else {
        """
        trim_galore \\
            --paired \\
            --cores $task.cpus \\
            ${reads[0]} ${reads[1]}

        # Rename outputs to match expected pattern
        mv *_val_1.fq ${meta.id}_trimmed_R1.fq
        mv *_val_2.fq ${meta.id}_trimmed_R2.fq
        """
    }
}
