reads = Channel.fromPath( 'data/ggal/*.fq' )

process foo {
    input:
    file 'sample.fastq' from reads
    script:
    """
    echo your_command --reads sample.fastq
    """
}