workflow {
    channel.fromFilePairs("data/reads/*/*_R{1,2}.fastq.gz")
        .view()
}
