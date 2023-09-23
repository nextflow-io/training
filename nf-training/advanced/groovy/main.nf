include { FASTP } from './modules/local/fastp/main.nf'

workflow {
    params.input = "https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq/samplesheet/v3.4/samplesheet_test.csv"

    Channel.fromPath(params.input)
    | splitCsv(header: true)
    | view
}