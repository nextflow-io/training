#!/usr/bin/env nextflow

/*
 * Pipeline parameters
 */
params {
    // Primary input (file of input files, one per line)
    reads_bam: Path = "${projectDir}/data/sample_bams.txt"

    // Output directory
    outdir: String = "results_genomics"

    // Accessory files
    reference: Path = "${projectDir}/data/ref/ref.fasta"
    reference_index: Path = "${projectDir}/data/ref/ref.fasta.fai"
    reference_dict: Path = "${projectDir}/data/ref/ref.dict"
    intervals: Path = "${projectDir}/data/ref/intervals.bed"

    // Base name for final output file
    cohort_name: String = "family_trio"
}

include { SAMTOOLS_INDEX } from './modules/samtools/index/main.nf'
include { GATK_HAPLOTYPECALLER } from './modules/gatk/haplotypecaller/main.nf'
include { GATK_JOINTGENOTYPING } from './modules/gatk/jointgenotyping/main.nf'


workflow {

    // Create input channel from a text file listing input file paths
    reads_ch = channel.fromPath(params.reads_bam).splitText()

    // Load the file paths for the accessory files (reference and intervals)
    ref_file        = file(params.reference)
    ref_index_file  = file(params.reference_index)
    ref_dict_file   = file(params.reference_dict)
    intervals_file  = file(params.intervals)

    // Create index file for input BAM file
    SAMTOOLS_INDEX(reads_ch)

    // Call variants from the indexed BAM file
    GATK_HAPLOTYPECALLER(
        SAMTOOLS_INDEX.out,
        ref_file,
        ref_index_file,
        ref_dict_file,
        intervals_file
    )

    // Collect variant calling outputs across samples
    all_gvcfs_ch = GATK_HAPLOTYPECALLER.out.vcf.collect()
    all_idxs_ch = GATK_HAPLOTYPECALLER.out.idx.collect()

    // Combine GVCFs into a GenomicsDB data store and apply joint genotyping
    GATK_JOINTGENOTYPING(
        all_gvcfs_ch,
        all_idxs_ch,
        intervals_file,
        params.cohort_name,
        ref_file,
        ref_index_file,
        ref_dict_file
    )
}
