nextflow.enable.dsl=2

params.dbtype = 'nr'
params.prot = 'data/prots/*.tfa'
proteins = Channel.fromPath(params.prot)

process find {
  input:
  file fasta
  val type

  when:
  fasta.name =~ /^BB11.*/ && type == 'nr'

  script:
  """
  echo blastp -query $fasta -db nr
  """
}

workflow{

    result = find(proteins, params.dbtype)

}
