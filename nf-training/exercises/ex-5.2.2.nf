reads = Channel.fromPath( 'data/ggal/*_1.fq' )

process foo {
    input:
    file all_samples from reads.collect()
    script:
    """
    cat $all_samples | head -n 20
    """
}
