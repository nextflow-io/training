/*
 * fastp - Adapter trimming and quality filtering
 *
 * Different container than FastQC - each process has its own
 * isolated software environment. No version conflicts!
 */
process FASTP {
    tag "$meta.id"
    container 'biocontainers/fastp:0.23.4'
    publishDir "${params.outdir}/trimmed", mode: 'copy'

    cpus 4
    memory '8.GB'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}_trimmed_R{1,2}.fastq.gz"), emit: reads
    path "*.{json,html}", emit: reports

    script:
    def prefix = meta.id
    // Mock implementation for demonstration
    // Real: fastp -i ${reads[0]} -I ${reads[1]} -o ${prefix}_trimmed_R1.fastq.gz ...
    """
    echo "Running fastp on ${meta.id}..."
    sleep 3  # Simulate processing time

    # Pass through reads (mock trimming)
    zcat ${reads[0]} | gzip > ${prefix}_trimmed_R1.fastq.gz
    zcat ${reads[1]} | gzip > ${prefix}_trimmed_R2.fastq.gz

    # Create mock reports
    echo '{"summary": {"before_filtering": {"total_reads": 1000000}, "after_filtering": {"total_reads": 950000}}}' > ${prefix}.json
    echo "<html><body><h1>fastp Report: ${prefix}</h1><p>Reads passed filter: 95%</p></body></html>" > ${prefix}.html

    echo "fastp complete for ${meta.id}"
    """
}
