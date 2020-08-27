params.reads = 'data/ggal/*_1.fq'
params.transcript = 'data/ggal/transcriptome.fa'

transcript = file(params.transcript)
reads = Channel.fromPath(params.reads)
alingners = ['salmon','kallisto']

process foo {
    echo true
    input:
    each aligner from alingners
    file fasta from transcript 
    file sample from reads
    script:
    """
    echo $aligner --ref $fasta --reads $sample
    """
}
