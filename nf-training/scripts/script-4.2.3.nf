nextflow.enable.dsl=2

params.genome = "$baseDir/data/ggal/transcriptome.fa"

process foo {
    input:
    path genome
    script:
    """
    echo your_command --reads $genome
    """
}

workflow {
	
    foo (params.genome)

}