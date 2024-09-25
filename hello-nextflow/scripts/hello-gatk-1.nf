/*
 * Pipeline parameters
 */

// Primary input
params.reads_bam = "${projectDir}/data/bam/reads_mother.bam"

/*
 * Generate BAM index file
 */
process SAMTOOLS_INDEX {

    container 'community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464'

    publishDir 'results', mode: 'copy'

    input:
        path input_bam

    output:
        path "${input_bam}.bai"

    """
    samtools index '$input_bam'
    """
}


workflow {

    // Create input channel (single file via CLI parameter)
    reads_ch = Channel.fromPath(params.reads_bam)

    // Create index file for input BAM file
    SAMTOOLS_INDEX(reads_ch)
    
}