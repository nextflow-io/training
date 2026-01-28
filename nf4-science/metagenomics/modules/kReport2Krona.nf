process K_REPORT_TO_KRONA {
    tag "${sample_id}"
    publishDir "$params.outdir/${sample_id}", mode:'copy'
    container "community.wave.seqera.io/library/krakentools:1.2--db94e0b19cfa397b"

    input:
    tuple val(sample_id), path(b_report), path(bracken)

    output:
    tuple val("${sample_id}"), path("${sample_id}.b.krona.txt")

    script:
    """
    kreport2krona.py -r ${b_report} \
    -o ${sample_id}.b.krona.txt \
    --no-intermediate-ranks
    """
}
