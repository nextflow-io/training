process FASTP {
    tag "$meta.id"
    container 'quay.io/biocontainers/fastp:0.23.4--hadf994f_2'
    publishDir "${params.outdir}/fastp", mode: 'copy'

    input:
    // TODO: Define input - same pattern as FASTQC
    ???

    output:
    // TODO: Define outputs - trimmed reads, JSON report, HTML report
    // Important: The trimmed reads output needs 'emit: reads' for downstream processes
    ???

    script:
    // TODO: Add the fastp command
    // Hint: Use ${reads[0]} and ${reads[1]} for input files
    // Hint: Use ${meta.id} for output file prefixes
    """
    ???
    """
}
