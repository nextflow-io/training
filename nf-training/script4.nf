nextflow.enable.dsl=2

/* 
 * pipeline input parameters 
 */
params.reads = "$baseDir/data/ggal/gut_{1,2}.fq"
params.transcriptome = "$baseDir/data/ggal/transcriptome.fa"
params.multiqc = "$baseDir/multiqc"
params.outdir = "results"

log.info """\
         R N A S E Q - N F   P I P E L I N E    
         ===================================
         transcriptome: ${params.transcriptome}
         reads        : ${params.reads}
         outdir       : ${params.outdir}
         """
         .stripIndent()

 
/* 
 * define the `index` process that create a binary index 
 * given the transcriptome file
 */
process index {

    input:
    path transcriptome

    output:
    path 'index'

    script:
    """
    salmon index --threads $task.cpus -t $transcriptome -i index
    """
}

process quantification {
     
    input:
    path index 
    tuple val(sample_id), path(reads)
 
    output:
    path "$sample_id"
 
    script:
    """
    salmon quant --threads $task.cpus --libType=U -i $index -1 ${reads[0]} -2 ${reads[1]} -o $sample_id
    """
}



workflow {

    index_ch = index(Channel.from(params.transcriptome))

    Channel
    .fromFilePairs( params.reads, checkIfExists: true )
    .set { read_pairs_ch } 

    quant_ch = quantification(index_ch, read_pairs_ch)

}
