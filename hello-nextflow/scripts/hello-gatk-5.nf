
params.baseDir = "/workspace/gitpod/nf-training" 
$baseDir = params.baseDir

params.reads_bam = "${baseDir}/data/gatk/samplesheet.csv"

params.genome_reference = "${baseDir}/data/gatk/ref/ref.fasta"
params.genome_reference_index = "${baseDir}/data/gatk/ref/ref.fasta.fai"
params.genome_reference_dict = "${baseDir}/data/gatk/ref/ref.dict"
params.calling_intervals = "${baseDir}/data/gatk/intervals.list"

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

    reads_ch = Channel.fromPath(params.reads_bam)
                        .splitCsv(header: true)
                        .map{row -> [row.id, file(row.reads_bam)]}

    SAMTOOLS_INDEX(reads_ch)

    GATK_HAPLOTYPECALLER(
        SAMTOOLS_INDEX.out,
        params.genome_reference,
        params.genome_reference_index,
        params.genome_reference_dict,
        params.calling_intervals
    )
}
