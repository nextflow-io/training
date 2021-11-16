params.genome = "$baseDir/data/ggal/transcriptome.fa"

process foo {
    input:
    path genome from params.genome
    script:
    """
    echo your_command --reads $genome
    """
}