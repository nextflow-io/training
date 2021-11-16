params.ncbi_api_key = '<Your API key here>'

params.accession = ['ERR908507', 'ERR908506']
reads = Channel.fromSRA(params.accession, apiKey: params.ncbi_api_key)

process fastqc {
    input:
    tuple sample_id, file(reads_file) from reads

    output:
    file("fastqc_${sample_id}_logs") into fastqc_ch

    script:
    """
    mkdir fastqc_${sample_id}_logs
    fastqc -o fastqc_${sample_id}_logs -f fastq -q ${reads_file}
    """
}