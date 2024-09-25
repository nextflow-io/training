// Include modules
include { SAMTOOLS_INDEX       } from './modules/local/samtools/index/main.nf'
include { GATK_HAPLOTYPECALLER } from './modules/local/gatk/haplotypecaller/main.nf'
include { GATK_JOINTGENOTYPING } from './modules/local/gatk/jointgenotyping/main.nf'

workflow {

    // Create input channel from BAM files
    bam_ch = Channel.fromPath(params.reads_bam, checkIfExists: true)


    // Reference objects
    ref_file               = file(params.reference)
    ref_index_file         = file(params.reference_index)
    ref_dict_file          = file(params.reference_dict)
    calling_intervals_file = file(params.calling_intervals)

    // Create index file for input BAM file
    SAMTOOLS_INDEX(bam_ch)

    // Call variants from the indexed BAM file
    GATK_HAPLOTYPECALLER(
        SAMTOOLS_INDEX.out,
        ref_file,
        ref_index_file,
        ref_dict_file,
        calling_intervals_file
    )

    all_vcfs = GATK_HAPLOTYPECALLER.out[0].collect()
    all_tbis = GATK_HAPLOTYPECALLER.out[1].collect()

    // Consolidate GVCFs and apply joint genotyping analysis
    GATK_JOINTGENOTYPING(
        all_vcfs,
        all_tbis,
        params.cohort_name,
        ref_file,
        ref_index_file,
        ref_dict_file,
        calling_intervals_file
    )
}
