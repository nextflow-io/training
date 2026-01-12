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
    """
    echo "Running FastQC on ${meta.id}..."
    sleep 1

    for read in ${reads}; do
        BASENAME=\$(basename \$read .fastq.gz)
        echo "<html><body><h1>FastQC: \${BASENAME}</h1></body></html>" > \${BASENAME}_fastqc.html
        echo "data" > \${BASENAME}_fastqc_data.txt
        zip -q \${BASENAME}_fastqc.zip \${BASENAME}_fastqc_data.txt
        rm \${BASENAME}_fastqc_data.txt
    done
    """
}
