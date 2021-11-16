reads = Channel.fromPath( 'data/ggal/*.fq' )

process foo {
    input:
    file sample from reads.collect()
    script:
    """
    echo your_command --reads $sample
    """
}