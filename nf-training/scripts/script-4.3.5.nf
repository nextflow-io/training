nextflow.enable.dsl=2

reads_ch = Channel.fromFilePairs('data/ggal/*_{1,2}.fq')

process foo {
  input:
    tuple val(sample_id), file(sample_id)
  output:
    tuple val(sample_id), file('sample.bam')
  script:
  """
    echo your_command_here --reads $sample_id > sample.bam
  """
}


workflow{

    bam_ch = foo(reads_ch)
    bam_ch.view()

}
