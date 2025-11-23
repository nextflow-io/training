workflow {
    channel.fromPath("data/samplesheet.csv")
        .splitCsv( header:true )
        .map { row ->
            def meta = [id:row.id, repeat:row.repeat, type:row.type]
            [
                meta,
                [
                    file(row.fastq1, checkIfExists: true),
                    file(row.fastq2, checkIfExists: true)
                ]
            ]
        }
        .view()
}
