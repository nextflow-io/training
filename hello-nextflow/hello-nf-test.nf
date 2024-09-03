// Include modules
include { SAMTOOLS_INDEX } from './modules/local/samtools/index/main.nf'
include { GATK_HAPLOTYPECALLER } from './modules/local/gatk/haplotypecaller/main.nf'
include { GATK_JOINTGENOTYPING } from './modules/local/gatk/jointgenotyping/main.nf'

workflow {

    // Create input channel from BAM files
    // We convert it to a tuple with the file name and the file path
    // See https://www.nextflow.io/docs/latest/script.html#getting-file-attributes
    bam_ch = Channel.fromPath(params.reads_bam, checkIfExists: true)

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
        SAMTOOLS_INDEX.out,
        ref_ch,
        ref_index_ch,
        ref_dict_ch,
        calling_intervals_ch
    )

    all_vcfs = GATK_HAPLOTYPECALLER.out[0].collect()
    all_tbis = GATK_HAPLOTYPECALLER.out[1].collect()

    // Consolidate GVCFs and apply joint genotyping analysis
    GATK_JOINTGENOTYPING(
        all_vcfs,
        all_tbis,
        params.cohort_name,
        ref_ch,
        ref_index_ch,
        ref_dict_ch,
        calling_intervals_ch
    )
}
