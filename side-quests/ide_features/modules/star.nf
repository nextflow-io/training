process STAR_ALIGN {
    tag "${sample_id}"
    publishDir "${params.output}/alignments", mode: 'copy'
    container 'biocontainers/star:2.7.10a_cv1'
    cpus 8
    memory '32.GB'

    input:
    tuple val(sample_id), path(reads)
    path reference_dir

    output:
    tuple val(sample_id), path("${sample_id}_Aligned.sortedByCoord.out.bam")

    script:
    """
    STAR \\
        --runThreadN ${task.cpus} \\
        --genomeDir ${reference_dir} \\
        --readFilesIn ${reads} \\
        --readFilesCommand zcat \\
        --outFileNamePrefix ${sample_id}_ \\
        --outSAMtype BAM SortedByCoordinate \\
        --outSAMunmapped Within \\
        --outSAMattributes Standard
    """

    stub:
    """
    touch ${sample_id}_Aligned.sortedByCoord.out.bam
    """
}
