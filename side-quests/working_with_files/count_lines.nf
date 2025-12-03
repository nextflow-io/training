process COUNT_LINES {
    debug true

    input:
    val input_file

    script:
    """
    set -o pipefail
    echo "Processing file: $input_file"
    gzip -dc $input_file | wc -l
    """
}

workflow {
    myFile = 'data/patientA_rep1_normal_R1_001.fastq.gz'
    COUNT_LINES(myFile)
}
