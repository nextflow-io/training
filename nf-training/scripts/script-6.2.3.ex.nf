params.reads = "$baseDir/data/ggal/*_1.fq"

Channel
    .fromPath( params.reads )
    .set { read_ch }

process concat {
  tag "Concat all files"

  input:
  path '*' from read_ch.collect()

  output:
  path 'top_10_lines' into concat_ch

  script:
  """
  cat * > concatenated.txt
  head -n 20 concatenated.txt > top_10_lines
  """
  }

concat_ch.view()
}

concat_ch.view{ it }