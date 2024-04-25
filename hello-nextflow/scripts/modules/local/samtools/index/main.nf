/*
 * Generate BAM index file
 */
process SAMTOOLS_INDEX {

    container 'quay.io/biocontainers/samtools:1.19.2--h50ea8bc_1' 

    input:
        tuple val(id), path(input_bam)

    output:
        tuple val(id), path(input_bam), path("${input_bam}.bai")

    """
    samtools index '$input_bam'

    """
}
