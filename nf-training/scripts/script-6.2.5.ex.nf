params.reads = "$baseDir/data/ggal/*_1.fq"
params.transcriptome_file = "$baseDir/data/ggal/transcriptome.fa"
methods= ['salmon', 'kallisto']

Channel
    .fromPath( params.reads )
    .set { read_ch }

process command {
  tag "Run_command"

  input:
  path reads from read_ch
  path transcriptome from params.transcriptome_file
  each mode from methods

  output:
  path result into concat_ch

  script:
  """
  echo $mode $reads $transcriptome > result
  """
  }

concat_ch
    .view { "To run : ${it.text}" }