workflow {
    getSampleIdAndReplicate = { sample -> [ sample.subMap(['id', 'repeat']), sample ] }
    samplesheet = Channel.fromPath("./data/samplesheet.csv")
                        .splitCsv(header: true)
    normal_samples = samplesheet
                        .filter { sample -> sample.type == 'normal' }
                        .map ( getSampleIdAndReplicate )
                        .dump(tag: 'normal')
    tumor_samples = samplesheet
                        .filter { sample -> sample.type == "tumor" }
                        .map ( getSampleIdAndReplicate )
                        .dump(tag: 'tumor')
    joined_samples = normal_samples
                        .join(tumor_samples)
                        .dump(tag: 'joined', pretty: true)
    grouped_samples = joined_samples.map { samples, normal, tumor ->
                        [
                            samples.subMap('id'),
                            normal,
                            tumor
                        ]
                    }
                    .groupTuple()
                    .dump(tag: 'grouped', pretty: true)

    intervals = Channel.of('chr1', 'chr2', 'chr3')
                    .dump(tag: "intervals")

    grouped_samples.combine(intervals)
                        .dump(tag: 'combined')
}
