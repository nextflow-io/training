process ANALYZE_READS {
    tag "${meta.id}"

    publishDir "results/${meta.type}/${meta.id}/${meta.replicate}", mode: 'copy'

    input:
    tuple val(meta), path(files)

    output:
    tuple val(meta.id), path("${meta.id}_stats.txt")

    script:
    """
    echo "Sample metadata: ${meta.id}" > ${meta.id}_stats.txt
    echo "Replicate: ${meta.replicate}" >> ${meta.id}_stats.txt
    echo "Type: ${meta.type}" >> ${meta.id}_stats.txt
    echo "Read 1: ${files[0]}" >> ${meta.id}_stats.txt
    echo "Read 2: ${files[1]}" >> ${meta.id}_stats.txt
    echo "Read 1 size: \$(gunzip -dc ${files[0]} | wc -l | awk '{print \$1/4}') reads" >> ${meta.id}_stats.txt
    echo "Read 2 size: \$(gunzip -dc ${files[1]} | wc -l | awk '{print \$1/4}') reads" >> ${meta.id}_stats.txt
    """
}

workflow {
    // Create a file object from a string path

    ch_files = channel.fromFilePairs('data/*_R{1,2}_001.fastq.gz')
    ch_samples = ch_files.map { id, files ->
        def (sample, replicate, type) = id.tokenize('_')
        [
            [
                id: sample,
                replicate: replicate.replace('rep', ''),
                type: type
            ],
            files
        ]
    }
    ANALYZE_READS(ch_samples)
}
