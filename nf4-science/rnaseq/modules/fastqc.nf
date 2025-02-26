process FASTQC {
    input:
    path read

    output:
    path "${read.simpleName}_fastqc.zip", emit: zip
    path "${read.simpleName}_fastqc.html", emit: html

    script:
    """
    fastqc $read
    """
}
