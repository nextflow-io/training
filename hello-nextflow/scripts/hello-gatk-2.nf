/*
 * Pipeline parameters
 */


// Primary input
params.reads_bam = "${projectDir}/../data/bam/reads_mother.bam"

// Accessory files
params.reference = "${workflow.projectDir}/../data/ref/ref.fasta"
params.reference_index = "${workflow.projectDir}/../data/ref/ref.fasta.fai"
params.reference_dict = "${workflow.projectDir}/../data/ref/ref.dict"
params.calling_intervals = "${workflow.projectDir}/../data/ref/intervals.bed"

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

/*
 * Call variants with GATK HapolotypeCaller in GVCF mode
 */
process GATK_HAPLOTYPECALLER {

    container "community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867"
    conda "bioconda::gatk4=4.5.0.0"

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
    bam_ch = Channel.of(params.reads_bam)

    // Create reference channels using the fromPath channel factory
    // The collect converts from a queue channel to a value channel
    // See https://www.nextflow.io/docs/latest/channel.html#channel-types for details
    ref_ch               = Channel.fromPath(params.reference, checkIfExists: true).collect()
    ref_index_ch         = Channel.fromPath(params.reference_index, checkIfExists: true).collect()
    ref_dict_ch          = Channel.fromPath(params.reference_dict, checkIfExists: true).collect()
    calling_intervals_ch = Channel.fromPath(params.calling_intervals, checkIfExists: true).collect()


    // Create index file for input BAM file
    SAMTOOLS_INDEX(bam_ch)

    // Call variants from the indexed BAM file
    GATK_HAPLOTYPECALLER(
        bam_ch,
        SAMTOOLS_INDEX.out,
        ref_ch,
        ref_index_ch,
        ref_dict_ch,
        calling_intervals_ch
    )
}
