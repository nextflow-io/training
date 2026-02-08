#!/usr/bin/env nextflow

/*
 * Pipeline parameters
 */
params {
    // Primary input (file of input files, one per line)
    reads_bam: Path = "${projectDir}/data/sample_bams.txt"

    // Accessory files
    reference: Path = "${projectDir}/data/ref/ref.fasta"
    reference_index: Path = "${projectDir}/data/ref/ref.fasta.fai"
    reference_dict: Path = "${projectDir}/data/ref/ref.dict"
    intervals: Path = "${projectDir}/data/ref/intervals.bed"
}

// Include modules
include { SAMTOOLS_INDEX } from './modules/samtools_index.nf'
include { GATK_HAPLOTYPECALLER } from './modules/gatk_haplotypecaller.nf'

workflow {

    main:
    // Create input channel from a CSV file listing input file paths
    reads_ch = Channel.fromPath(params.reads_bam)
            .splitCsv()
            .map { line -> file(line[0]) }

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

    publish:
    indexed_bam = SAMTOOLS_INDEX.out
    vcf = GATK_HAPLOTYPECALLER.out.vcf
    vcf_idx = GATK_HAPLOTYPECALLER.out.idx
}

output {
    indexed_bam {
        path 'bam'
    }
    vcf {
        path 'vcf'
    }
    vcf_idx {
        path 'vcf'
    }
}
