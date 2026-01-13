include { FASTP } from './modules/fastp.nf'
include { TRIMGALORE } from './modules/trimgalore.nf'
include { GENERATE_REPORT } from './modules/generate_report.nf'

def validateInputs() {
    // Check input parameter is provided
    if (!params.input) {
        error("Input CSV file path not provided. Please specify --input <file.csv>")
    }

    // Check CSV file exists
    if (!file(params.input).exists()) {
        error("Input CSV file not found: ${params.input}")
    }
}

def separateMetadata(row) {
    def sample_meta = [
        id: row.sample_id.toLowerCase(),
        organism: row.organism,
        tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
        depth: row.sequencing_depth.toInteger(),
        quality: row.quality_score.toDouble()
    ]
    def run_id = row.run_id?.toUpperCase() ?: 'UNSPECIFIED'
    sample_meta.run = run_id
    def fastq_path = file(row.file_path)

    def m = (fastq_path.name =~ /^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$/)
    def file_meta = m ? [
        sample_num: m[0][2].toInteger(),
        lane: m[0][3],
        read: m[0][4],
        chunk: m[0][5]
    ] : [:]

    def priority = sample_meta.quality > 40 ? 'high' : 'normal'

    // Validate data makes sense
    if (sample_meta.depth < 30000000) {
        log.warn "Low sequencing depth for ${sample_meta.id}: ${sample_meta.depth}"
    }

    return tuple(sample_meta + file_meta + [priority: priority], fastq_path)
}

workflow {
    validateInputs()

    ch_samples = channel.fromPath(params.input)
        .splitCsv(header: true)
        .map{ row -> separateMetadata(row) }

    // Filter out invalid or low-quality samples
    ch_valid_samples = ch_samples
        .filter { meta, reads ->
            meta.id && meta.organism && meta.depth >= 25000000
        }

    trim_branches = ch_valid_samples
        .branch { meta, reads ->
            fastp: meta.organism == 'human' && meta.depth >= 30000000
            trimgalore: true
        }

    ch_fastp = FASTP(trim_branches.fastp)
    ch_trimgalore = TRIMGALORE(trim_branches.trimgalore)
    GENERATE_REPORT(ch_samples)

    workflow.onComplete = {
        println ""
        println "Pipeline execution summary:"
        println "=========================="
        println "Completed at: ${workflow.complete}"
        println "Duration    : ${workflow.duration}"
        println "Success     : ${workflow.success}"
        println "workDir     : ${workflow.workDir}"
        println "exit status : ${workflow.exitStatus}"
        println ""

        if (workflow.success) {
            println "✅ Pipeline completed successfully!"
        } else {
            println "❌ Pipeline failed!"
            println "Error: ${workflow.errorMessage}"
        }
    }
}
