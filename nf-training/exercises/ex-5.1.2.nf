params.aligner = 'kallisto'

process foo {
  script:
  alignCommand(params.aligner)
}


def alignCommand( aligner ) {
  if( aligner == 'kallisto' ) 
    """
    echo kallisto --reads /some/data.fastq
    """
  else if( aligner == 'salmon' ) 
    """
    echo salmon --reads /some/data.fastq
    """
  else 
    throw new IllegalArgumentException("Unknown aligner $aligner")
}
