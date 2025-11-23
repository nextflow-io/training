process ANALYZE_READS {
    tag "${meta.id}"

    publishDir "results/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(fastqs)

    output:
    tuple val(meta.id), path("${meta.id}_stats.txt")

    script:
    """
    echo "Sample metadata: ${meta.id}" > ${meta.id}_stats.txt
    echo "Replicate: ${meta.replicate}" >> ${meta.id}_stats.txt
    echo "Type: ${meta.type}" >> ${meta.id}_stats.txt
    echo "Read 1: ${fastqs[0]}" >> ${meta.id}_stats.txt
    echo "Read 2: ${fastqs[1]}" >> ${meta.id}_stats.txt
    echo "Read 1 size: \$(gunzip -dc ${fastqs[0]} | wc -l | awk '{print \$1/4}') reads" >> ${meta.id}_stats.txt
    echo "Read 2 size: \$(gunzip -dc ${fastqs[1]} | wc -l | awk '{print \$1/4}') reads" >> ${meta.id}_stats.txt
    """
}

workflow {
    // Create a file object from a string path

    ch_fastq = channel.fromFilePairs('data/sampleA_rep1_normal_R{1,2}_001.fastq.gz')
    ch_samples = ch_fastq.map { id, fastqs ->
        def (sample, replicate, type, readNum) = id.tokenize('_')
        [
            [
                id: sample,
                replicate: replicate.replace('rep', ''),
                type: type,
                readNum: readNum,
            ],
            fastqs
        ]
    }
    ANALYZE_READS(ch_samples)
}
