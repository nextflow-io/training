process HISAT2_ALIGN {
    input:
    path read
    path index
    path splice_sites

    output:
    path "${read.simpleName}.bam", emit: bam
    path "${read.simpleName}.hisat2.log", emit: log

    script:
    def index_base = index[0].toString() - ~/.\d.ht2/
    """
    hisat2 -x $index_base -U $read --known-splicesite-infile $splice_sites --new-summary --summary-file ${read.simpleName}.hisat2.log | \
    samtools view -bS - > ${read.simpleName}.bam
    """
}
