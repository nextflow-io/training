nextflow.enable.dsl=2

reads_ch = Channel.fromFilePairs('data/ggal/*_{1,2}.fq')
params.transcriptome = "data/ggal/transcriptome.fa"
index = Channel.fromPath(params.transcriptome)

process makeBams {
    publishDir "$baseDir/bam_files", mode: 'copy'
    echo true

    input:
    file index
    tuple val(name), file(reads)

    output:
    tuple val(name), file ('*.bam')

    """
    echo STAR --genomeDir $index --readFilesIn $reads
    """
}


workflow {
    
    makeBams(index,reads_ch)

}
