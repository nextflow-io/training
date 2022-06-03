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
 * define the `index` process that create a binary index 
 * given the transcriptome file
 */
process index {
    
    input:
    path transcriptome from params.transcriptome_file
     
    output:
    path 'salmon_index' into index_ch

    script:       
    """
    salmon index --threads $task.cpus -t $transcriptome -i salmon_index
    """
}


Channel 
    .fromFilePairs( params.reads, checkIfExists: true )
    .into { read_pairs_ch; read_pairs2_ch } 

process quantification {
    tag "$sample_id"
         
    input:
    path salmon_index from index_ch
    tuple val(sample_id), path(reads) from read_pairs_ch
 
    output:
    path sample_id into quant_ch
 
    script:
    """
    salmon quant --threads $task.cpus --libType=U -i $salmon_index -1 ${reads[0]} -2 ${reads[1]} -o $sample_id
    """
}

process fastqc {
    tag "FASTQC on $sample_id"

    input:
    tuple val(sample_id), path(reads) from read_pairs2_ch

    output:
    path "fastqc_${sample_id}_logs" into fastqc_ch


    script:
    """
    mkdir fastqc_${sample_id}_logs
    fastqc -o fastqc_${sample_id}_logs -f fastq -q ${reads}
    """  
}  
 

process multiqc {
    publishDir params.outdir, mode:'copy'
       
    input:
    path '*' from quant_ch.mix(fastqc_ch).collect()
    
    output:
    path 'multiqc_report.html'
     
    script:
    """
    multiqc . 
    """
} 


workflow.onComplete { 
	log.info ( workflow.success ? "\nDone! Open the following report in your browser --> $params.outdir/multiqc_report.html\n" : "Oops .. something went wrong" )
}
