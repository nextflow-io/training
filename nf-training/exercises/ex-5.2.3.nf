params.reads = 'data/ggal/*_1.fq'
params.transcript = 'data/ggal/transcriptome.fa'

transcript = file(params.transcript)
reads = Channel.fromPath(params.reads)

process foo {
    echo true
    input:
    file fasta from transcript 
    file sample from reads 
    script:
    """
    echo mapper_command_here --ref $fasta --reads $sample
    """
}
