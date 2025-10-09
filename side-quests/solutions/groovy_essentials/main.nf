include { FASTP } from './modules/fastp.nf'
include { TRIMGALORE } from './modules/trimgalore.nf'
include { GENERATE_REPORT } from './modules/generate_report.nf'

def validateInputs() {
    // Check CSV file exists
    if (!file(params.input ?: './data/samples.csv').exists()) {
        error("Input CSV file not found: ${params.input ?: './data/samples.csv'}")
    }

    // Warn if output directory already exists
    if (file(params.outdir ?: 'results').exists()) {
        log.warn "Output directory already exists: ${params.outdir ?: 'results'}"
    }

    // Check for required genome parameter
    if (params.run_gatk && !params.genome) {
        error("Genome reference required when running GATK. Please provide --genome")
    }
}

def separateMetadata(row) {
    def sample_meta = [
        id: row.sample_id.toLowerCase(),
        organism: row.organism?.toLowerCase() ?: 'unknown',
        tissue: row.tissue_type?.replaceAll('_', ' ')?.toLowerCase() ?: 'unknown',
        depth: row.sequencing_depth.toInteger(),
        quality: row.quality_score?.toDouble() ?: 0.0
    ]
    def fastq_path = file(row.file_path)

    def m = (fastq_path.name =~ /^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$/)
    def file_meta = m ? [
        sample_num: m[0][2].toInteger(),
        lane: m[0][3],
        read: m[0][4],
        chunk: m[0][5]
    ] : [:]

    def priority = sample_meta.quality > 40 ? 'high' : 'normal'
    return [sample_meta + file_meta + [priority: priority], fastq_path]
}

workflow {
    validateInputs()

    ch_samples = Channel.fromPath("./data/samples.csv")
        .splitCsv(header: true)
        .map(separateMetadata)
        .filter { meta, reads ->
            meta.organism != 'unknown' && (meta.quality ?: 0) > 0
        }

    GENERATE_REPORT(ch_samples)

    trim_branches = ch_samples
        .branch { meta, reads ->
            fastp: meta.organism == 'human' && meta.depth >= 30000000
            trimgalore: true
        }

    ch_fastp = FASTP(trim_branches.fastp)
    ch_trimgalore = TRIMGALORE(trim_branches.trimgalore)
}
