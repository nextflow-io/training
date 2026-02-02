workflow {
    ch_samples = Channel.fromPath("./data/samples.csv")
        .splitCsv(header: true)
        .view()
}
