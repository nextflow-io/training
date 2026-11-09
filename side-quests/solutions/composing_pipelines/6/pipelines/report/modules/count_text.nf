nextflow.enable.types = true

include { Sample } from '../types'

/*
 * Count the characters or words in a sample's text file
 */
process COUNT_TEXT {
    tag "${sample.id}"

    input:
    sample: Sample
    metric: String

    output:
    file("${sample.id}.count")

    script:
    def flag = metric == 'words' ? '-w' : '-m'
    """
    echo "${sample.id}: \$(wc ${flag} < ${sample.file}) ${metric}" > ${sample.id}.count
    """
}
