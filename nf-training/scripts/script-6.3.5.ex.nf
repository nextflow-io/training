reads_ch = Channel.fromFilePairs('data/ggal/*_{1,2}.fq')

process foo {
  input:
    tuple val(sample_id), file(sample_files) from reads_ch
  output:
    tuple val(sample_id), file("${sample_id}.bam") into bam_ch
  script:
  """
    echo your_command_here --reads $sample_id > ${sample_id}.bam
  """
}