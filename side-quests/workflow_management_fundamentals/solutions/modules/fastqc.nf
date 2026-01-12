process FASTQC {
    tag "$meta.id"
    publishDir "${params.outdir}/fastqc", mode: 'copy'

    cpus 2
    memory '2.GB'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.html"), emit: html
    tuple val(meta), path("*.zip"), emit: zip

    script:
    // Mock FastQC - creates placeholder reports
    // In production, you would use: container 'biocontainers/fastqc:0.12.1'
    """
    echo "Running FastQC on ${meta.id}..."
    sleep 1  # Simulate processing time

    # Count reads and create mock report
    READ_COUNT=\$(zcat ${reads[0]} | grep -c '^@' || echo "0")

    for read in ${reads}; do
        BASENAME=\$(basename \$read .fastq.gz)

        # Create mock HTML report
        cat > \${BASENAME}_fastqc.html << EOF
<!DOCTYPE html>
<html>
<head><title>FastQC Report: \${BASENAME}</title></head>
<body>
<h1>FastQC Report for \${BASENAME}</h1>
<p>Sample: ${meta.id}</p>
<p>Organism: ${meta.organism}</p>
<p>Total Sequences: \${READ_COUNT}</p>
<p>Sequence Length: 40 bp</p>
<p>%GC: 50</p>
<p><strong>Status: PASS</strong></p>
</body>
</html>
EOF

        # Create mock ZIP (just a placeholder)
        echo "FastQC data for \${BASENAME}" > \${BASENAME}_fastqc_data.txt
        zip -q \${BASENAME}_fastqc.zip \${BASENAME}_fastqc_data.txt
        rm \${BASENAME}_fastqc_data.txt
    done

    echo "FastQC complete for ${meta.id}"
    """
}
