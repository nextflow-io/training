process QUAST {
    tag "$meta.id"
    publishDir "${params.outdir}/quast", mode: 'copy'

    cpus 2
    memory '4.GB'

    input:
    tuple val(meta), path(assembly)

    output:
    tuple val(meta), path("${meta.id}/*"), emit: report

    script:
    // Mock QUAST - creates placeholder quality reports
    // In production, you would use: container 'biocontainers/quast:5.2.0'
    """
    echo "Running QUAST on ${meta.id}..."
    sleep 1  # Simulate processing time

    mkdir -p ${meta.id}

    # Count contigs and calculate mock stats
    CONTIG_COUNT=\$(grep -c '^>' ${assembly} || echo "0")
    TOTAL_LENGTH=\$(grep -v '^>' ${assembly} | tr -d '\\n' | wc -c)

    # Create mock report TSV
    cat > ${meta.id}/report.tsv << EOF
Assembly\t${meta.id}
# contigs\t\${CONTIG_COUNT}
Total length\t\${TOTAL_LENGTH}
Largest contig\t300
GC (%)\t50.00
N50\t250
N90\t180
L50\t2
L90\t3
EOF

    # Create mock HTML report
    cat > ${meta.id}/report.html << EOF
<!DOCTYPE html>
<html>
<head><title>QUAST Report: ${meta.id}</title></head>
<body>
<h1>QUAST Report for ${meta.id}</h1>
<h2>Assembly Statistics</h2>
<table border="1">
<tr><th>Metric</th><th>Value</th></tr>
<tr><td># contigs</td><td>\${CONTIG_COUNT}</td></tr>
<tr><td>Total length</td><td>\${TOTAL_LENGTH}</td></tr>
<tr><td>Largest contig</td><td>300</td></tr>
<tr><td>GC (%)</td><td>50.00</td></tr>
<tr><td>N50</td><td>250</td></tr>
</table>
<p><strong>Quality: GOOD</strong></p>
</body>
</html>
EOF

    echo "QUAST complete for ${meta.id}"
    """
}
