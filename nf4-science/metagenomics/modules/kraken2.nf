process KRAKEN2 {
    tag "${sample_id}"
    publishDir "$params.outdir/${sample_id}", mode:'copy'
    container "community.wave.seqera.io/library/kraken2:2.14--83aa57048e304f01"

    input:
    tuple val(sample_id), path(reads_1), path(reads_2), path(sam)
    path kraken2_db

    output:
    tuple val("${sample_id}"), path("${sample_id}.k2report"), path("${sample_id}.kraken2")

    script:
    """
    kraken2 --db $kraken2_db --threads 2 \
    --report ${sample_id}.k2report \
    --report-minimizer-data \
    --minimum-hit-groups 2 \
    --gzip-compressed \
    --paired \
    ${reads_1} ${reads_2} > ${sample_id}.kraken2
    """
}
