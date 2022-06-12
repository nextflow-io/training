/* 
 * pipeline input parameters 
 */
params.reads = "$projectDir/data/ggal/gut_{1,2}.fq"
params.transcriptome_file = "$projectDir/data/ggal/transcriptome.fa"
params.multiqc = "$projectDir/multiqc"
params.outdir = "results"

log.info """\
         R N A S E Q - N F   P I P E L I N E    
         ===================================
         transcriptome: ${params.transcriptome_file}
         reads        : ${params.reads}
         outdir       : ${params.outdir}
         """
         .stripIndent()

 
/* 
 * define the `index` process that creates a binary index 
 * given the transcriptome file
 */
process index {
    
    input:
    path transcriptome
     
    output:
    path 'salmon_index'

    script: 
    """
    salmon index --threads $task.cpus -t $transcriptome -i salmon_index
    """
}


Channel 
    .fromFilePairs( params.reads, checkIfExists: true )
    .set { read_pairs_ch } 

process quantification {
     
    input:
    path salmon_index
    tuple val(sample_id), path(reads)
 
    output:
    path "$sample_id"
 
    script:
    """
    salmon quant --threads $task.cpus --libType=U -i $salmon_index -1 ${reads[0]} -2 ${reads[1]} -o $sample_id
    """
}

process fastqc {
    tag "FASTQC on $sample_id"

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

workflow {

    index_ch = index(Channel.from(params.transcriptome))

    Channel
    .fromFilePairs( params.reads, checkIfExists: true )
    .set { read_pairs_ch } 

    quant_ch = quantification(index_ch, read_pairs_ch)

    fastqc_ch = fastqc(read_pairs_ch)   

}
