// This will cause problems!
myFile = 'data/patientA_rep1_normal_R1_001.fastq.gz'

process COUNT_LINES {
    debug true

    input:
    val fastq_file

    script:
    """
    echo "Processing file: $fastq_file"
    gzip -dc $fastq_file | wc -l
    """
}

workflow {
    COUNT_LINES(myFile)
}
