species = ['cat','dog', 'sloth']
sequences = ['AGATAG','ATGCTCT', 'ATCCCAA']

process align {
  input:
  val x from species
  val seq from sequences

  output:
  file "${x}.aln" into genomes

  """
  echo align -in $seq > ${x}.aln
  """
}