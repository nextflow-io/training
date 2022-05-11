/* 
 * pipeline input parameters 
 */
params.reads = "$projectDir/data/ggal/*_{1,2}.fq"
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
 * define the `index` process that create a binary index 
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

workflow {

    index_ch = index(Channel.from(params.transcriptome))

}