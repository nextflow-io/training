/*
 * SPAdes - Genome assembly
 *
 * Note the higher resource requirements - Nextflow uses these
 * to schedule jobs. Won't run 10 SPAdes jobs if you only have
 * enough memory for 2.
 */
process SPADES {
    tag "$meta.id"
    container 'biocontainers/spades:3.15.5'
    publishDir "${params.outdir}/assemblies", mode: 'copy'

    // SPAdes needs more resources than FastQC
    // Nextflow schedules accordingly
    cpus 8
    memory '16.GB'
    time '6.h'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}/contigs.fasta"), emit: assembly

    script:
    // Mock implementation for demonstration
    // Real: spades.py -1 ${reads[0]} -2 ${reads[1]} -o ${meta.id} --threads $task.cpus
    """
    echo "Running SPAdes assembly on ${meta.id}..."
    sleep 5  # Simulate processing time (real SPAdes takes hours!)

    mkdir -p ${meta.id}
    cat > ${meta.id}/contigs.fasta << 'FASTA'
>NODE_1_length_50000_cov_45.5
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
>NODE_2_length_35000_cov_42.1
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
>NODE_3_length_28000_cov_38.7
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
FASTA

    echo "SPAdes assembly complete for ${meta.id}"
    echo "Generated 3 contigs"
    """
}
