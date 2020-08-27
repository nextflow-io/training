params.outdir = 'results'

process INDEX {
    tag "$transcriptome.simpleName"

    input:
    path transcriptome

    output:
    path 'index'

    script:
    """
    salmon index --threads $task.cpus -t $transcriptome -i index
    """
}


process QUANT {
    tag "$pair_id"

    input:
    path index
    tuple pair_id, path(reads)

    output:
    path pair_id

    script:
    """
    salmon quant --threads $task.cpus --libType=U -i $index -1 ${reads[0]} -2 ${reads[1]} -o $pair_id
    """
}

process FASTQC {
    tag "FASTQC on $sample_id"
    publishDir params.outdir

    input:
    tuple sample_id, path(reads)

    output:
    path "fastqc_${sample_id}_logs"


    script:
    """
    mkdir fastqc_${sample_id}_logs
    fastqc -o fastqc_${sample_id}_logs -f fastq -q ${reads}
    """
}


process MULTIQC {
    publishDir params.outdir, mode:'copy'
    errorStrategy 'ignore'
    input:
    path 'data*/*'
    path config

    output:
    path 'multiqc_report.html'

    script:
    """
    cp $config/* .
    echo "custom_logo: \$PWD/logo.png" >> multiqc_config.yaml
    multiqc -v .
    """
}