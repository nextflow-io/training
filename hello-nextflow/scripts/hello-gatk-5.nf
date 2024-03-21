/*
 * Pipeline parameters
 */

// Execution environment setup
params.baseDir = "/workspace/gitpod/hello-nextflow" 
$baseDir = params.baseDir

// Primary input (samplesheet in CSV format with ID and file path, one sample per line)
params.reads_bam = "${baseDir}/data/samplesheet.csv"

// Accessory files
params.genome_reference = "${baseDir}/data/ref/ref.fasta"
params.genome_reference_index = "${baseDir}/data/ref/ref.fasta.fai"
params.genome_reference_dict = "${baseDir}/data/ref/ref.dict"
params.calling_intervals = "${baseDir}/data/intervals.list"

/*
 * Generate BAM index file
 */
process SAMTOOLS_INDEX {

    container 'quay.io/biocontainers/samtools:1.19.2--h50ea8bc_1' 

    input:
        tuple val(id), path(input_bam)

    output:
        tuple path(input_bam), path("${input_bam}.bai")

    """
    samtools index '$input_bam'

    """
}

/*
 * Call variants with GATK HapolotypeCaller in GVCF mode
 */
process GATK_HAPLOTYPECALLER {

    container "broadinstitute/gatk:4.5.0.0"

    input:
        tuple path(input_bam), path(input_bam_index)
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

    // Create input channel from samplesheet in CSV format
    reads_ch = Channel.fromPath(params.reads_bam)
                        .splitCsv(header: true)
                        .map{row -> [row.id, file(row.reads_bam)]}

    // Create index file for input BAM file
    SAMTOOLS_INDEX(reads_ch)

    // Call variants from the indexed BAM file
    GATK_HAPLOTYPECALLER(
        SAMTOOLS_INDEX.out,
        params.genome_reference,
        params.genome_reference_index,
        params.genome_reference_dict,
        params.calling_intervals
    )
}
