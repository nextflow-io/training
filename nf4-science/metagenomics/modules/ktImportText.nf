process KT_IMPORT_TEXT {
    tag "${sample_id}"
    publishDir "$params.outdir/${sample_id}", mode:'copy'
    container "community.wave.seqera.io/library/krona:2.8.1--2f750080982f027e"

    input:
    tuple val(sample_id), path(krona_txt)

    output:
    path "${sample_id}.krona.html"

    script:
    """
    ktImportText ${krona_txt} \
    -o ${sample_id}.krona.html
    """
}
