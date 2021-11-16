process makeBams {
    publishDir "/some/directory/bam_files", mode: 'copy'

    input:
    file index from index_ch
    tuple val(name), file(reads) from reads_ch

    output:
    tuple val(name), file ('*.bam') into star_aligned

    """
    echo STAR --genomeDir $index --readFilesIn $reads
    """
}