params.outdir = 'my-results'
params.prot = 'data/prots/*.tfa'
proteins = Channel.fromPath(params.prot)


process blastSeq {
    publishDir "$params.outdir/bam_files", mode: 'copy'

    input:
    file fasta from proteins

    output:
    file ('*.txt') into blast_ch

    """
    echo blastp $fasta > ${fasta}_result.txt
    """
}

blast_ch.view()