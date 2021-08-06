nextflow.enable.dsl=2

reads = Channel.fromPath( 'data/ggal/*.fq' )

process foo {
    input:
    file sample
    script:
    """
    ls -lsh $sample
    """
}

workflow{

   foo(reads.collect())

}
