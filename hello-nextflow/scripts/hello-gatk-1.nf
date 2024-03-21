/*
 * Pipeline parameters
 */

// Execution environment setup
params.baseDir = "/workspace/gitpod/hello-nextflow" 
$baseDir = params.baseDir

// Primary input
params.reads_bam = "${baseDir}/data/bam/reads_mother.bam"

/*
 * Generate BAM index file
 */
process SAMTOOLS_INDEX {

    container 'quay.io/biocontainers/samtools:1.19.2--h50ea8bc_1' 

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
    reads_ch = Channel.from(params.reads_bam)

    // Create index file for input BAM file
    SAMTOOLS_INDEX(reads_ch)
}