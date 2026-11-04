workflow {
    main:
    ch_samples = channel.fromPath("./data/samples.csv")
        .splitCsv(header: true)
        .view()

    publish:
    reports = channel.empty()
}

output {
    reports {
        path 'reports'
    }
}
