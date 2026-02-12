#!/usr/bin/env nextflow

// Module INCLUDE statements
include { SAMTOOLS_INDEX } from './modules/samtools_index.nf'
include { GATK_HAPLOTYPECALLER } from './modules/gatk_haplotypecaller.nf'

/*
 * Pipeline parameters
 */
params {
    // Primary input (file of input files, one per line)
    input: Path

    // Accessory files
    reference: Path
    reference_index: Path
    reference_dict: Path
    intervals: Path
}

workflow {

    main:
    // Create input channel from a CSV file listing input file paths
    reads_ch = channel.fromPath(params.input)
            .splitCsv(header: true)
            .map { row -> file(row.reads_bam) }

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
