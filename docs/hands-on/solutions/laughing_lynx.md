<div class="formalpara-title">

**Solution**

</div>

``` nextflow
process '6C_ASE_knownSNPs' {
  tag "$sampleId"
  publishDir "$params.results/$sampleId"

  input:
    path genome from params.genome
    path index from genome_index_ch
    path dict from genome_dict_ch
    tuple val(sampleId), path(vcf), path(bam), path(bai) from grouped_vcf_bam_bai_ch

  output:
    path "ASE.tsv"

  script:
  """
  echo "${bam.join('\n')}" > bam.list

  java -jar $GATK -R ${genome} \
                  -T ASEReadCounter \
                  -o ASE.tsv \
                  -I bam.list \
                  -sites ${vcf}
  """
}
```
