nextflow.enable.types = true

/*
 * Combine the per-sample counts into a single report
 */
process SUMMARIZE {
    input:
    counts: Bag<Path>
    title: String
    banner: Path

    output:
    file('report.txt')

    script:
    """
    cat ${banner} > report.txt
    echo "${title}" >> report.txt
    echo "" >> report.txt
    sort --parallel=${task.cpus} ${counts.join(' ')} >> report.txt
    """
}
