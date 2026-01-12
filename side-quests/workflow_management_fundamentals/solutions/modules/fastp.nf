process FASTP {
    tag "$meta.id"
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
    // Mock fastp - creates trimmed reads and reports
    // In production, you would use: container 'biocontainers/fastp:0.23.4'
    """
    echo "Running fastp on ${meta.id}..."
    sleep 2  # Simulate processing time

    # "Trim" reads by creating new output files
    # In reality, fastp filters and trims - here we just pass through
    zcat ${reads[0]} | gzip > ${prefix}_trimmed_R1.fastq.gz
    zcat ${reads[1]} | gzip > ${prefix}_trimmed_R2.fastq.gz

    # Create mock JSON report with static values
    cat > ${prefix}.json << 'ENDJSON'
{
    "summary": {
        "before_filtering": {
            "total_reads": 3,
            "total_bases": 120,
            "q20_rate": 0.95,
            "q30_rate": 0.90
        },
        "after_filtering": {
            "total_reads": 3,
            "total_bases": 114,
            "q20_rate": 0.98,
            "q30_rate": 0.95
        }
    },
    "filtering_result": {
        "passed_filter_reads": 3,
        "low_quality_reads": 0,
        "too_short_reads": 0
    }
}
ENDJSON

    # Create mock HTML report
    cat > ${prefix}.html << 'ENDHTML'
<!DOCTYPE html>
<html>
<head><title>fastp Report</title></head>
<body>
<h1>fastp Report</h1>
<h2>Summary</h2>
<table>
<tr><th>Metric</th><th>Before</th><th>After</th></tr>
<tr><td>Total Reads</td><td>3</td><td>3</td></tr>
<tr><td>Q20 Rate</td><td>95%</td><td>98%</td></tr>
<tr><td>Q30 Rate</td><td>90%</td><td>95%</td></tr>
</table>
<p><strong>Filtering: PASS</strong></p>
</body>
</html>
ENDHTML

    echo "fastp complete for ${meta.id}"
    """
}
