process MULTIQC {
    publishDir "${params.output}/multiqc", mode: 'copy'
    container 'biocontainers/multiqc:1.13_cv1'
    cpus 2
    memory '4.GB'

    input:
    path qc_files

    output:
    path "multiqc_report.html"
    path "multiqc_data/"

    script:
    """
    multiqc \\
        --force \\
        --title "Pipeline QC Report" \\
        --filename multiqc_report.html \\
        .
    """

    stub:
    """
    mkdir multiqc_data
    touch multiqc_report.html
    """
}

process CUSTOM_DUMPSOFTWAREVERSIONS {
    publishDir "${params.output}/pipeline_info", mode: 'copy'

    input:
    path versions

    output:
    path "software_versions.yml"
    path "software_versions_mqc.yml"

    script:
    """
    cat ${versions} > software_versions.yml
    echo 'id: "software_versions"' > software_versions_mqc.yml
    echo 'section_name: "Software Versions"' >> software_versions_mqc.yml
    """
}
