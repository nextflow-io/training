/*
 * Pipeline parameters
 */

// Primary input
params.bams = "${workflow.projectDir}/../data/bam/reads_mother.bam"


/*
 * Generate BAM index file
 */
process SAMTOOLS_INDEX {

    container 'community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464'
    conda "bioconda::samtools=1.19.2"

    input:
        path input_bam

    output:
        path "${input_bam}.bai"

    """
    samtools index '$input_bam'

    """
}

workflow {

    // Create input channel
    reads_ch = Channel.fromPath(params.bams, checkIfExists: true)

    // Create index file for input BAM file
    SAMTOOLS_INDEX(reads_ch)
}
