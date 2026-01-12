/*
 * FastQC - Quality control for sequencing reads
 *
 * Note the 'container' directive - this ensures:
 * - Exact version 0.12.1 is always used
 * - No manual installation required
 * - Same results on any machine
 */
process FASTQC {
    tag "$meta.id"
    container 'biocontainers/fastqc:0.12.1'
    publishDir "${params.outdir}/fastqc", mode: 'copy'

    // Declare resource requirements
    // Nextflow uses these to schedule jobs intelligently
    cpus 2
    memory '2.GB'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.html"), emit: html
    tuple val(meta), path("*.zip"), emit: zip

    script:
    // Mock implementation for demonstration
    // Real implementation would be: fastqc -q -t $task.cpus ${reads}
    """
    echo "Running FastQC on ${meta.id}..."
    sleep 2  # Simulate processing time

    for read in ${reads}; do
        BASENAME=\$(basename \$read .fastq.gz)
        echo "<html><body><h1>FastQC Report: \${BASENAME}</h1><p>Quality: PASS</p></body></html>" > \${BASENAME}_fastqc.html
        echo "FastQC data" > \${BASENAME}_data.txt
        zip -q \${BASENAME}_fastqc.zip \${BASENAME}_data.txt
        rm \${BASENAME}_data.txt
    done

    echo "FastQC complete for ${meta.id}"
    """
}
