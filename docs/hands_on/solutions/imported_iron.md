<div class="formalpara-title">

**Solution**

</div>

``` nextflow
process '5_rnaseq_call_variants' {
  tag "$sampleId" 

  input:
    path genome from params.genome 
    path index from genome_index_ch 
    path dict from genome_dict_ch 
    tuple val(sampleId), path(bam), path(bai) from final_output_ch.groupTuple() 

  output:
    tuple val(sampleId), path('final.vcf') into vcf_files 

  script:
  """
  echo "${bam.join('\n')}" > bam.list

  # Variant calling
  java -jar $GATK -T HaplotypeCaller \
                  -R $genome -I bam.list \
                  -dontUseSoftClippedBases \
                  -stand_call_conf 20.0 \
                  -o output.gatk.vcf.gz

  # Variant filtering
  java -jar $GATK -T VariantFiltration \
                  -R $genome -V output.gatk.vcf.gz \
                  -window 35 -cluster 3 \
                  -filterName FS -filter "FS > 30.0" \
                  -filterName QD -filter "QD < 2.0" \
                  -o final.vcf 
  """
}
```

- [`tag`](https://www.nextflow.io/docs/latest/process.html#tag) line with the using the sample id as the tag.

- the genome fasta file.

- the genome index from the `genome_index_ch` channel created in the process `1A_prepare_genome_samtools`.

- the genome dictionary from the `genome_dict_ch` channel created in the process `1B_prepare_genome_picard`.

- the sets grouped by sampleID from the `final_output_ch` channel created in the process `4_rnaseq_gatk_recalibrate`.

- the set containing the sample ID and final VCF file.

- the line specifing the name resulting final vcf file.
