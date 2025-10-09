process GENERATE_REPORT {

    publishDir 'results/reports', mode: 'copy'

    input:
    tuple val(meta), path(reads)

    output:
    path "${meta.id}_report.txt"

    script:
    def report_type = meta.priority == 'high' ? 'PRIORITY' : 'STANDARD'
    """
    echo "=== ${report_type} SAMPLE REPORT ===" > ${meta.id}_report.txt
    echo "Processing ${reads}" >> ${meta.id}_report.txt
    echo "Sample: ${meta.id}" >> ${meta.id}_report.txt
    echo "Quality: ${meta.quality}" >> ${meta.id}_report.txt
    echo "Priority: ${meta.priority}" >> ${meta.id}_report.txt
    echo "---" >> ${meta.id}_report.txt
    echo "Processed by: \${USER}" >> ${meta.id}_report.txt
    echo "Hostname: \$(hostname)" >> ${meta.id}_report.txt
    echo "Date: \$(date)" >> ${meta.id}_report.txt
    """
}
