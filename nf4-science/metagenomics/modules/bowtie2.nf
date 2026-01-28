process BOWTIE2 {
    tag "${sample_id}"
    publishDir "$params.outdir/${sample_id}", pattern: "*.sam", mode:'copy'
    container "community.wave.seqera.io/library/bowtie2:2.5.4--d51920539234bea7"

    input:
    tuple val(sample_id), path(reads)
    path bowtie2_index

    output:
    tuple val("${sample_id}"), path("${sample_id}.1"), path("${sample_id}.2"), path("${sample_id}.sam")

    script:
    """
    export BOWTIE2_INDEXES=/workspaces/training/nf4-science/metagenomics/data/yeast
    bowtie2 -x $bowtie2_index -1 ${reads[0]} -2 ${reads[1]} -p 2 -S ${sample_id}.sam --un-conc-gz ${sample_id}
    """
}
