process GENERATE_REPORT {

    publishDir 'results/reports', mode: 'copy'

    input:
    tuple val(meta), path(reads)

    output:
    path "${meta.id}_report.txt"

    script:
    """
    echo "Processing ${reads}" > ${meta.id}_report.txt
    echo "Sample: ${meta.id}" >> ${meta.id}_report.txt
    echo "Processed by: \${USER}" >> ${meta.id}_report.txt
    echo "Hostname: \$(hostname)" >> ${meta.id}_report.txt
    echo "Date: \$(date)" >> ${meta.id}_report.txt
    """
}
