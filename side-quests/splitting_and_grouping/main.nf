workflow {
    getSampleIdAndReplicate = { sample -> [ sample.subMap(['id', 'repeat']), sample ] }
    ch_samplesheet = Channel.fromPath("./data/samplesheet.csv")
                        .splitCsv(header: true)
    ch_normal_samples = ch_samplesheet
                        .filter { sample -> sample.type == 'normal' }
                        .map ( getSampleIdAndReplicate )
    ch_tumor_samples = ch_samplesheet
                        .filter { sample -> sample.type == "tumor" }
                        .map ( getSampleIdAndReplicate )
    ch_joined_samples = ch_normal_samples
                        .join(ch_tumor_samples)
    ch_intervals = Channel.of('chr1', 'chr2', 'chr3')

    ch_combined_samples = ch_joined_samples.combine(ch_intervals)
                        .map { grouping_key, normal, tumor, interval ->
                            [
                                grouping_key + [interval: interval],
                                normal,
                                tumor
                            ]

                        }

    ch_grouped_samples = ch_combined_samples.map { grouping_key, normal, tumor ->
                            [
                                grouping_key.subMap('id', 'interval'),
                                normal.subMap("bam"),
                                tumor.subMap("bam")
                            ]

                        }
                        .groupTuple()
                        .view()
}
