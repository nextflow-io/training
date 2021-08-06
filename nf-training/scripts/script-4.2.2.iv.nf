nextflow.enable.dsl=2

params.genome = "$baseDir/data/ggal/transcriptome.fa"

process foo {
    input:
    path genome
    script:
    """
    ls -lsh $genome
    """
}

workflow{

   foo(params.genome)

}
