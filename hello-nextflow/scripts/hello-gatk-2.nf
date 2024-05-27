/*
 * Pipeline parameters
 */

// Execution environment setup
params.projectDir = "/workspace/gitpod/hello-nextflow" 
$projectDir = params.projectDir

// Primary input
params.reads_bam = "${projectDir}/data/bam/reads_mother.bam"

// Accessory files
params.genome_reference = "${projectDir}/data/ref/ref.fasta"
params.genome_reference_index = "${projectDir}/data/ref/ref.fasta.fai"
params.genome_reference_dict = "${projectDir}/data/ref/ref.dict"
params.calling_intervals = "${projectDir}/data/intervals.list"

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

/*
 * Call variants with GATK HapolotypeCaller in GVCF mode
 */
process GATK_HAPLOTYPECALLER {

    container "docker.io/broadinstitute/gatk:4.5.0.0"

    input:
        path input_bam
        path input_bam_index
        path ref_fasta
        path ref_index
        path ref_dict
        path interval_list

    output:
        path "${input_bam}.g.vcf"
        path "${input_bam}.g.vcf.idx"

    """
    gatk HaplotypeCaller \
        -R ${ref_fasta} \
        -I ${input_bam} \
        -O ${input_bam}.g.vcf \
        -L ${interval_list} \
        -ERC GVCF
    """
}

workflow {

    // Create input channel
    reads_ch = Channel.of(params.reads_bam)

    // Create index file for input BAM file
    SAMTOOLS_INDEX(reads_ch)

    // Call variants from the indexed BAM file
    GATK_HAPLOTYPECALLER(
        reads_ch,
        SAMTOOLS_INDEX.out,
        params.genome_reference,
        params.genome_reference_index,
        params.genome_reference_dict,
        params.calling_intervals
    )
}
