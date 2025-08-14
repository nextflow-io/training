#!/usr/bin/env nextflow

/*
 * Groovy Essentials Demo Workflow
 *
 * This workflow demonstrates essential Groovy concepts in Nextflow contexts
 * Follow along with docs/side_quests/groovy_essentials.md for detailed explanations
 */

// Basic pipeline parameters
params.input = "./data/samples.csv"
params.quality_threshold_min = 25
params.quality_threshold_high = 40

//=============================================================================
// SECTION 1: Nextflow vs Groovy Boundaries
//=============================================================================

workflow {
    println "=== Groovy Essentials Demo ==="

    // 1.1: Basic sample processing (Nextflow + Groovy)
    ch_samples = Channel.fromPath(params.input)
        .splitCsv(header: true)
        .map { row ->
            // Groovy: Map operations and string manipulation
            def sample_meta = [
                id: row.sample_id.toLowerCase(),
                organism: row.organism,
                tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                depth: row.sequencing_depth.toInteger(),
                quality: row.quality_score.toDouble()
            ]

            // Groovy: Conditional logic and string interpolation
            def priority = sample_meta.quality > params.quality_threshold_high ? 'high' : 'normal'

            // Nextflow: Return tuple for channel
            return [sample_meta + [priority: priority], file(row.file_path)]
        }

    // 1.2: Demonstrate Groovy vs Nextflow collect
    demonstrateCollectDifference()

    // 2: String processing and pattern matching
    demonstrateStringProcessing()

    // 3: Strategy selection and conditional logic
    ch_enriched_samples = ch_samples.map { meta, file_path ->
        def strategy = selectAnalysisStrategy(meta)
        return [meta + strategy, file_path]
    }

        // 4: Validation and error handling
    ch_validated_samples = ch_enriched_samples.map { meta, file_path ->
        def errors = validateSample(meta)
        if (errors) {
            log.warn "Sample ${meta.id} has validation issues: ${errors.join(', ')}"
        }
        return [meta, file_path]
    }

    // 5: Essential Groovy operators demonstration
    demonstrateGroovyOperators()

    // 6: Advanced closures demonstration
    demonstrateAdvancedClosures()

    // 7: Collection operations demonstration
    demonstrateCollectionOperations()

    // Display final processed samples
    ch_validated_samples.view { meta, file_path ->
        "Final: ${meta.id} (${meta.organism}, ${meta.approach}, priority: ${meta.priority})"
    }
}

//=============================================================================
// SECTION 1.2: Collect Confusion Demo
//=============================================================================

def demonstrateCollectDifference() {
    println "\n=== Groovy vs Nextflow Collect ==="

    def sample_ids = ['sample_001', 'sample_002', 'sample_003']

    // Groovy collect: transform each element
    def formatted_ids = sample_ids.collect { id ->
        id.toUpperCase().replace('SAMPLE_', 'SPECIMEN_')
    }
    println "Groovy collect result: ${formatted_ids}"

    // Nextflow collect: gather channel elements
    ch_collected = Channel.of('sample_001', 'sample_002', 'sample_003')
        .collect()
    ch_collected.view { "Nextflow collect result: ${it}" }
}

//=============================================================================
// SECTION 2: String Processing and Pattern Matching
//=============================================================================

def demonstrateStringProcessing() {
    println "\n=== String Processing Demo ==="

    // Pattern matching for sample file names
    def sample_files = [
        'Human_Liver_001.fastq',
        'mouse_brain_002.fastq',
        'SRR12345678.fastq'
    ]

    // Simple pattern to extract organism and tissue
    def pattern = ~/^(\w+)_(\w+)_(\d+)\.fastq$/

    sample_files.each { filename ->
        def matcher = filename =~ pattern
        if (matcher) {
            println "${filename} -> Organism: ${matcher[0][1]}, Tissue: ${matcher[0][2]}, ID: ${matcher[0][3]}"
        } else {
            println "${filename} -> No standard pattern match"
        }
    }
}

// Function to extract sample metadata from filename
def parseSampleName(String filename) {
    def pattern = ~/^(\w+)_(\w+)_(\d+)\.fastq$/
    def matcher = filename =~ pattern

    if (matcher) {
        return [
            organism: matcher[0][1].toLowerCase(),
            tissue: matcher[0][2].toLowerCase(),
            sample_id: matcher[0][3],
            valid: true
        ]
    } else {
        return [
            filename: filename,
            valid: false
        ]
    }
}

//=============================================================================
// SECTION 3: Conditional Logic and Strategy Selection
//=============================================================================

def selectAnalysisStrategy(Map sample_meta) {
    def strategy = [:]

    // Sequencing depth determines processing approach
    if (sample_meta.depth < 10_000_000) {
        strategy.approach = 'low_depth'
        strategy.processes = ['quality_check', 'simple_alignment']
        strategy.sensitivity = 'high'
    } else if (sample_meta.depth < 50_000_000) {
        strategy.approach = 'standard'
        strategy.processes = ['quality_check', 'bwa_alignment', 'variant_calling']
        strategy.sensitivity = 'medium'
    } else {
        strategy.approach = 'high_depth'
        strategy.processes = ['quality_check', 'bwa_alignment', 'variant_calling', 'structural_variants']
        strategy.sensitivity = 'low'
    }

    // Organism-specific adjustments
    switch(sample_meta.organism) {
        case 'human':
            strategy.reference = 'GRCh38'
            strategy.known_variants = 'dbSNP'
            break
        case 'mouse':
            strategy.reference = 'GRCm39'
            strategy.known_variants = 'mgp_variants'
            break
        default:
            strategy.reference = 'custom'
            strategy.known_variants = null
    }

    // Quality-based modifications
    if (sample_meta.quality < 30) {
        strategy.extra_qc = true
        strategy.processes = ['extensive_qc'] + strategy.processes
    }

    return strategy
}

//=============================================================================
// SECTION 4: Error Handling and Validation
//=============================================================================

def validateSample(Map sample) {
    def errors = []

    // Check required fields
    if (!sample.id) errors << "Missing sample ID"
    if (!sample.organism) errors << "Missing organism"
    if (!sample.quality) errors << "Missing quality score"

    // Validate data types and ranges
    if (sample.quality && (sample.quality < 0 || sample.quality > 50)) {
        errors << "Quality score out of range (0-50)"
    }

    if (sample.depth && sample.depth < 1000) {
        errors << "Sequencing depth too low (<1000)"
    }

    return errors
}

def processSample(Map sample) {
    try {
        def errors = validateSample(sample)
        if (errors) {
            throw new RuntimeException("Validation failed: ${errors.join(', ')}")
        }

        // Simulate processing
        println "Processing sample: ${sample.id}"
        return [success: true, sample: sample]

    } catch (Exception e) {
        println "Error processing ${sample.id}: ${e.message}"
        return [
            success: false,
            sample: sample,
            error: e.message,
            // Partial result for recovery
            partial_result: [id: sample.id, status: 'failed']
        ]
    }
}

//=============================================================================
// SECTION 5: Essential Groovy Operators and Patterns
//=============================================================================

def demonstrateGroovyOperators() {
    println "\n=== Essential Groovy Operators Demo ==="

    // 5.1: Safe navigation and Elvis operators for robust data handling
    def simulateQcMetrics = { hasData ->
        if (hasData) {
            return [
                summary: [
                    before_filtering: [total_reads: 25_000_000],
                    after_filtering: [q30_rate: 95.2]
                ],
                adapter_cutting: [adapter_trimmed_reads: 2_500_000]
            ]
        } else {
            return [:]  // Empty or incomplete data
        }
    }

    println "Safe navigation with QC data:"
    def goodQc = simulateQcMetrics(true)
    def badQc = simulateQcMetrics(false)

    println "  Good QC total reads: ${goodQc?.summary?.before_filtering?.total_reads ?: 0}"
    println "  Bad QC total reads: ${badQc?.summary?.before_filtering?.total_reads ?: 0}"
    println "  Good QC with Elvis: ${goodQc?.summary?.after_filtering?.q30_rate ?: 'No data'}"

    // 5.2: Groovy Truth in workflow context
    println "Groovy Truth for workflow validation:"
    def samples = [
        [id: 'sample1', files: ['file1.fastq', 'file2.fastq']],
        [id: 'sample2', files: []],
        [id: 'sample3', files: null],
        [id: 'sample4', organism: 'human', files: ['data.fastq']]
    ]

    samples.each { sample ->
        // Groovy Truth simplifies validation
        def status = sample.files ? "Ready (${sample.files.size()} files)" : "No files"
        def organism = sample.organism ?: 'unknown'
        println "  ${sample.id}: ${status}, organism: ${organism}"
    }

    // 5.3: Slashy strings for filename parsing in workflows
    println "Filename parsing for pipeline routing:"
    def sampleFiles = [
        'Human_Liver_001_R1.fastq',
        'Mouse_Brain_002_R2.fastq',
        'SRR12345678.fastq',
        'invalid_file.txt'
    ]

    // Use slashy strings to avoid escaping
    def humanPattern = /^(Human)_(\w+)_(\d+)_R([12])\.fastq$/
    def mousePattern = /^(Mouse)_(\w+)_(\d+)_R([12])\.fastq$/
    def sraPattern = /^(SRR\d+)\.fastq$/

    sampleFiles.each { filename ->
        def route = 'unknown'
        if (filename =~ humanPattern) route = 'human_pipeline'
        else if (filename =~ mousePattern) route = 'mouse_pipeline'
        else if (filename =~ sraPattern) route = 'sra_pipeline'

        println "  ${filename} -> ${route}"
    }

    // 5.4: Numeric literals in bioinformatics context
    def evaluateSequencingDepth = { depth ->
        def category = 'unknown'
        def recommendation = 'unknown'

        if (depth < 1_000_000) {
            category = 'very_low'
            recommendation = 'Consider resequencing'
        } else if (depth < 10_000_000) {
            category = 'low'
            recommendation = 'Adequate for basic analysis'
        } else if (depth < 50_000_000) {
            category = 'good'
            recommendation = 'Suitable for variant calling'
        } else {
            category = 'excellent'
            recommendation = 'Perfect for comprehensive analysis'
        }

        return [category: category, recommendation: recommendation]
    }

    println "Depth analysis with readable numbers:"
    [500_000, 5_000_000, 25_000_000, 75_000_000].each { depth ->
        def analysis = evaluateSequencingDepth(depth)
        def formatted = depth.toString().replaceAll(/(\d)(?=(\d{3})+$)/, '$1,')
        println "  ${formatted} reads -> ${analysis.category}: ${analysis.recommendation}"
    }
}

//=============================================================================
// SECTION 6: Advanced Closures and Functional Programming
//=============================================================================

def demonstrateAdvancedClosures() {
    println "\n=== Advanced Closures Demo ==="

        // 6.1: Named closures for reusability
    def extractSampleInfo = { row ->
        [
            id: row.sample_id?.toLowerCase() ?: 'unknown',
            organism: row.organism ?: 'unknown',
            quality: (row.quality_score as Double) ?: 0.0
        ]
    }

    def addPriority = { meta ->
        meta + [priority: meta.quality > 40 ? 'high' : 'normal']
    }

    println "Named closures example:"
    def testRow = [sample_id: 'TEST_001', organism: 'human', quality_score: '42.5']
    def processed = addPriority(extractSampleInfo(testRow))
    println "  Processed: ${processed}"

    // 6.2: Function composition with >>
    def normalizeId = { meta ->
        meta + [id: meta.id.toLowerCase().replaceAll(/[^a-z0-9_]/, '_')]
    }

    def addQualityCategory = { meta ->
        def category = meta.quality > 40 ? 'excellent' :
                      meta.quality > 30 ? 'good' : 'acceptable'
        meta + [quality_category: category]
    }

    def enrichSample = normalizeId >> addQualityCategory
    def testSample = [id: 'Sample-001', quality: 42.0]
    def enriched = enrichSample(testSample)
    println "Function composition: ${enriched}"

    // 6.3: Currying example
    def qualityFilter = { threshold, meta -> meta.quality >= threshold }
    def highQualityFilter = qualityFilter.curry(40)
    def standardQualityFilter = qualityFilter.curry(30)

    println "Currying example:"
    println "  High quality filter (40+): ${highQualityFilter([quality: 42])}"
    println "  Standard quality filter (30+): ${standardQualityFilter([quality: 35])}"

    // 6.4: Closure with scope access
    def stats = [total: 0, high_quality: 0]
    def collectStats = { meta ->
        stats.total++
        if (meta.quality > 40) stats.high_quality++
        return meta
    }

    // Process some test data
    [[quality: 45], [quality: 30], [quality: 42]].each(collectStats)
    println "Statistics collected: ${stats}"
}

//=============================================================================
// SECTION 7: Collection Operations
//=============================================================================

def demonstrateCollectionOperations() {
    println "\n=== Collection Operations Demo ==="

    // Sample data with mixed quality and organisms
    def samples = [
        [id: 'sample_001', organism: 'human', quality: 42, files: ['data1.txt', 'data2.txt']],
        [id: 'sample_002', organism: 'mouse', quality: 28, files: ['data3.txt']],
        [id: 'sample_003', organism: 'human', quality: 35, files: ['data4.txt', 'data5.txt', 'data6.txt']],
        [id: 'sample_004', organism: 'rat', quality: 45, files: ['data7.txt']],
        [id: 'sample_005', organism: 'human', quality: 30, files: ['data8.txt', 'data9.txt']]
    ]

    // findAll - filter collections based on conditions
    def high_quality_samples = samples.findAll { it.quality > 40 }
    println "High quality samples: ${high_quality_samples.collect { it.id }.join(', ')}"

    // groupBy - group samples by organism
    def samples_by_organism = samples.groupBy { it.organism }
    println "Grouping by organism:"
    samples_by_organism.each { organism, sample_list ->
        println "  ${organism}: ${sample_list.size()} samples"
    }

    // unique - get unique organisms
    def organisms = samples.collect { it.organism }.unique()
    println "Unique organisms: ${organisms.join(', ')}"

    // flatten - flatten nested file lists
    def all_files = samples.collect { it.files }.flatten()
    println "All files: ${all_files.take(5).join(', ')}... (${all_files.size()} total)"

    // sort - sort samples by quality
    def sorted_by_quality = samples.sort { it.quality }
    println "Quality range: ${sorted_by_quality.first().quality} to ${sorted_by_quality.last().quality}"

    // count - count items matching condition
    def human_samples = samples.count { it.organism == 'human' }
    println "Human samples: ${human_samples} out of ${samples.size()}"

    // any/every - check conditions across collection
    def has_high_quality = samples.any { it.quality > 40 }
    def all_have_files = samples.every { it.files.size() > 0 }
    println "Has high quality samples: ${has_high_quality}"
    println "All samples have files: ${all_have_files}"

    // Spread operator demonstration
    demonstrateSpreadOperator()
}

def demonstrateSpreadOperator() {
    println "\n=== Spread Operator Demo ==="

    def file_paths = [
        '/data/sample1.fastq',
        '/data/sample2.fastq',
        '/results/output1.bam',
        '/results/output2.bam'
    ]

    // Convert to file objects
    def files = file_paths.collect { file(it) }

    // Using spread operator - equivalent to files.collect { it.getName() }
    def filenames = files*.getName()
    println "Filenames: ${filenames.join(', ')}"

    // Get all parent directories
    def parent_dirs = files*.getParent()*.getName()
    println "Parent directories: ${parent_dirs.unique().join(', ')}"

    // Get all extensions
    def extensions = files*.getExtension().unique()
    println "File types: ${extensions.join(', ')}"
}

//=============================================================================
// SAMPLE PROCESSES (for demonstration)
//=============================================================================

process QUALITY_FILTER {
    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}_filtered.fastq")

    script:
    // Groovy logic to determine parameters based on metadata
    def quality_threshold = meta.organism == 'human' ? 30 :
                           meta.organism == 'mouse' ? 28 : 25
    def min_length = meta.priority == 'high' ? 75 : 50

    // Conditional script sections
    def extra_qc = meta.priority == 'high' ? '--strict-quality' : ''

    """
    echo "Processing ${meta.id} (${meta.organism}, priority: ${meta.priority})"

    # Dynamic quality filtering based on sample characteristics
    fastp \\
        --in1 ${reads} \\
        --out1 ${meta.id}_filtered.fastq \\
        --qualified_quality_phred ${quality_threshold} \\
        --length_required ${min_length} \\
        ${extra_qc}

    echo "Applied quality threshold: ${quality_threshold}"
    echo "Applied length threshold: ${min_length}"
    """

    stub:
    """
    echo "STUB: Processing ${meta.id} (${meta.organism}, priority: ${meta.priority})"
    touch ${meta.id}_filtered.fastq
    """
}

process JOINT_ANALYSIS {
    input:
    path sample_files  // This will be a list of files
    path reference

    output:
    path "joint_results.txt"

    script:
    // Transform file list into command arguments
    def file_args = sample_files.collect { file -> "--input ${file}" }.join(' ')
    def sample_names = sample_files.collect { file ->
        file.baseName.replaceAll(/\..*$/, '')
    }.join(',')

    """
    echo "Processing ${sample_files.size()} samples"
    echo "Sample names: ${sample_names}"

    # Use the transformed arguments in the actual command
    analysis_tool \\
        ${file_args} \\
        --reference ${reference} \\
        --output joint_results.txt \\
        --samples ${sample_names}
    """

    stub:
    """
    echo "STUB: Processing ${sample_files.size()} samples"
    echo "STUB: Sample names: ${sample_names}"
    touch joint_results.txt
    """
}
