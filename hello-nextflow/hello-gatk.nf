
params.baseDir = "/workspace/gitpod/nf-training" 
$baseDir = params.baseDir

params.reads_bam = "${baseDir}/data/gatk/bam/reads_mother.bam"

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

    reads_ch = Channel.from(params.reads_bam)

    SAMTOOLS_INDEX(reads_ch)
}