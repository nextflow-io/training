workflow {
    ch_samples = channel.fromPath("./data/samplesheet.csv")
        .splitCsv(header: true)
        .map{ row ->
            [[id:row.id, repeat:row.repeat, type:row.type], row.bam]
        }

    getSampleIdAndReplicate = { meta, file -> [ meta.subMap(['id', 'repeat']), file ] }

    ch_normal_samples = ch_samples
        .filter { meta, file -> meta.type == 'normal' }
        .map ( getSampleIdAndReplicate )
    ch_tumor_samples = ch_samples
        .filter { meta, file -> meta.type == 'tumor' }
        .map ( getSampleIdAndReplicate )
    ch_joined_samples = ch_normal_samples
        .join(ch_tumor_samples)
    ch_intervals = channel.of('chr1', 'chr2', 'chr3')

    ch_combined_samples = ch_joined_samples
        .combine(ch_intervals)
        .map { grouping_key, normal, tumor, interval ->
            [
                grouping_key + [interval: interval],
                normal,
                tumor
            ]
        }

    ch_grouped_samples = ch_combined_samples
        .map { grouping_key, normal, tumor ->
            [
                grouping_key.subMap('id', 'interval'),
                normal,
                tumor
            ]
            }
            .groupTuple()
            .view()

}
