nextflow.enable.types = true

/*
 * Combine the per-sample counts into a single report
 */
process SUMMARIZE {
    input:
    counts: Bag<Path>

    output:
    file('report.txt')

    script:
    """
    banner=${projectDir}/assets/banner.txt
    if [ -f \$banner ]; then
        cat \$banner > report.txt
    fi
    echo "${params.title}" >> report.txt
    echo "" >> report.txt
    sort --parallel=${task.cpus} ${counts.join(' ')} >> report.txt
    """
}
