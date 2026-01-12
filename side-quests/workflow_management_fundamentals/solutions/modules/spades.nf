process SPADES {
    tag "$meta.id"
    publishDir "${params.outdir}/assemblies", mode: 'copy'

    cpus 8
    memory '16.GB'
    time '6.h'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}/contigs.fasta"), emit: assembly

    script:
    // Mock SPAdes - creates a placeholder assembly
    // In production, you would use: container 'biocontainers/spades:3.15.5'
    """
    echo "Running SPAdes assembly on ${meta.id}..."
    sleep 3  # Simulate processing time (real SPAdes takes hours)

    mkdir -p ${meta.id}

    # Create mock assembly with a few "contigs"
    cat > ${meta.id}/contigs.fasta << EOF
>NODE_1_length_5000_cov_50.0
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
>NODE_2_length_3000_cov_45.5
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
>NODE_3_length_2000_cov_40.0
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
EOF

    echo "SPAdes assembly complete for ${meta.id}"
    echo "Generated 3 contigs"
    """
}
