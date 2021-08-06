process example {
    script:
    """
    grep '^@S' data/ggal/gut_1.fq > Fastq_headers
    cat Fastq_headers | cut -d '.' -f2 > Uniq_Fastq_Names
    
    """
}
