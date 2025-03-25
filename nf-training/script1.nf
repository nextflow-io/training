#!/usr/bin/env nextflow

params.reads = "$projectDir/data/ggal/gut_{1,2}.fq"
params.transcriptome_file = "$projectDir/data/ggal/transcriptome.fa"
params.multiqc = "$projectDir/multiqc"
params.outdir = "$projectDir/results"


log.info """\
         R N A S E Q - N F   P I P E L I N E
         ===================================
         transcriptome : ${params.transcriptome_file}
         reads         : ${params.reads}
         multiqc       : ${params.multiqc}
         outdir        : ${params.outdir}
         """
         .stripIndent()

/*
 * define o processo INDEX que cria um índice binário
 * dado um arquivo de transcriptoma
 */

process INDEX {
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
    index_ch = INDEX(params.transcriptome_file)
}
