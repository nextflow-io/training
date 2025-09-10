process COUNT_LINES {
    debug true

    input:
    val fastq_file

    script:
    """
    set -o pipefail
    echo "Processing file: $fastq_file"
    gzip -dc $fastq_file | wc -l
    """
}

workflow {
    myFile = 'data/patientA_rep1_normal_R1_001.fastq.gz'
    COUNT_LINES(myFile)
}
