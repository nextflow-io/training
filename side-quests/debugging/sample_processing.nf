#!/usr/bin/env nextflow

process COUNT_LINES {

    container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

    input:
    tuple val(sample_id), path(fastq)

    output:
    path "${sample_id}.count.txt"

    script:
    """
    lines=\$(zcat ${fastq} | wc -l)
    cowpy "${sample_id} has \${lines} lines" > ${sample_id}_count.txt
    """

    stub:
    """
    echo "${sample_id} has 0 lines" > ${sample_id}_count.txt
    """
}

process REPORT {

    input:
    path count_files

    output:
    path 'report.txt'

    script:
    """
    cat ${count_files} > report.txt
    """
}

workflow {

    main:
    samples_ch = channel
        .fromPath(params.input)
        .splitCsv(header: true)
        .map { row -> [row.sample_id, file(row.fastq_path)] }

    counts_ch = COUNT_LINES(samples_ch)
    report_ch = REPORT(counts_ch.collect())

    publish:
    report = report_ch
}

output {
    report {
        path '.'
    }
}
