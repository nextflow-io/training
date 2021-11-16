params.dbtype = 'nr'
params.prot = 'data/prots/*.tfa'
proteins = Channel.fromPath(params.prot)

process find {
  input:
  file fasta from proteins
  val type from params.dbtype

  when:
  fasta.name =~ /^BB11.*/ && type == 'nr'

  script:
  """
  echo blastp -query $fasta -db nr
  """
}