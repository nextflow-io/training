reads = Channel.fromPath( 'data/ggal/*.fq' )

process foo {
    input:
    file sample from reads
    script:
    """
    echo your_command --reads $sample
    """
}
