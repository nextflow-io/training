/*
 * QUAST - Assembly quality assessment
 */
process QUAST {
    tag "$meta.id"
    container 'biocontainers/quast:5.2.0'
    publishDir "${params.outdir}/quast", mode: 'copy'

    cpus 2
    memory '4.GB'

    input:
    tuple val(meta), path(assembly)

    output:
    tuple val(meta), path("${meta.id}/*"), emit: report

    script:
    // Mock implementation for demonstration
    // Real: quast.py $assembly -o ${meta.id} --threads $task.cpus
    """
    echo "Running QUAST on ${meta.id}..."
    sleep 1  # Simulate processing time

    mkdir -p ${meta.id}

    # Calculate mock stats from the assembly
    CONTIG_COUNT=\$(grep -c '^>' ${assembly} || echo "0")
    TOTAL_LENGTH=\$(grep -v '^>' ${assembly} | tr -d '\\n' | wc -c)

    cat > ${meta.id}/report.tsv << EOF
Assembly	${meta.id}
# contigs	\${CONTIG_COUNT}
Total length	\${TOTAL_LENGTH}
Largest contig	50000
GC (%)	50.00
N50	35000
EOF

    cat > ${meta.id}/report.html << EOF
<!DOCTYPE html>
<html>
<head><title>QUAST Report: ${meta.id}</title></head>
<body>
<h1>QUAST Report: ${meta.id}</h1>
<table border="1">
<tr><th>Metric</th><th>Value</th></tr>
<tr><td># contigs</td><td>\${CONTIG_COUNT}</td></tr>
<tr><td>Total length</td><td>\${TOTAL_LENGTH}</td></tr>
<tr><td>N50</td><td>35000</td></tr>
</table>
<p><strong>Quality: GOOD</strong></p>
</body>
</html>
EOF

    echo "QUAST complete for ${meta.id}"
    """
}
