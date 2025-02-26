process TRIM_GALORE {
    input:
    path read

    output:
    path "${read.simpleName}_trimmed.fq.gz", emit: trimmed_reads
    path "${read.simpleName}_trimming_report.txt", emit: trimming_reports

    script:
    """
    trim_galore --gzip -o . $read
    """
}
