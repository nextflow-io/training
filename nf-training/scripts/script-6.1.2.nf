params.compress = 'gzip'
params.file2compress = "$baseDir/data/ggal/transcriptome.fa"

process foo {

  input:
  path file from params.file2compress

  script:
  if( params.compress == 'gzip' )
    """
    gzip -c $file > ${file}.fa.gz
    """
  else if( params.compress == 'bzip2' )
    """
    bzip2 -c $file > ${file}.fa.bz2
    """
  else
    throw new IllegalArgumentException("Unknown aligner $params.compress")
}