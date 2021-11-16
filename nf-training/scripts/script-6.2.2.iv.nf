params.genome = 'data/ggal/transcriptome.fa'

process foo {
    input:
    file genome from params.genome
    script:
    """
    echo your_command --reads $genome
    """
}