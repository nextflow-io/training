params.compress = 'gzip'
params.file2compress = "$projectDir/data/ggal/transcriptome.fa"

process foo {

  input:
  path file

  script:
  if( params.compress == 'gzip' )
    """
    gzip -c $file > ${file}.gz
    """
  else if( params.compress == 'bzip2' )
    """
    bzip2 -c $file > ${file}.bz2
    """
  else
    throw new IllegalArgumentException("Unknown aligner $params.compress")
}   

workflow{
  foo(params.file2compress)
}
