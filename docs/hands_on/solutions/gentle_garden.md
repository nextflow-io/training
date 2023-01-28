<div class="formalpara-title">

**Solution**

</div>

``` nextflow
process '3_rnaseq_gatk_splitNcigar' {
  tag "$replicateId" 

  input:
    path genome from params.genome  
    path index from genome_index_ch  
    path genome_dict from genome_dict_ch  
    tuple val(replicateId), path(bam), path(bai) from aligned_bam_ch  

  output:
    tuple val(replicateId), path('split.bam'), path('split.bai') into splitted_bam_ch  

  script:
  """
  # SplitNCigarReads and reassign mapping qualities
  java -jar $GATK -T SplitNCigarReads \
                  -R $genome -I $bam \//
                  -o split.bam \//
                  -rf ReassignOneMappingQuality \
                  -RMQF 255 -RMQT 60 \
                  -U ALLOW_N_CIGAR_READS \
                  --fix_misencoded_quality_scores

  """
}
```

- [`tag`](https://www.nextflow.io/docs/latest/process.html#tag) line with the using the replicate id as the tag.

- the genome fasta file

- the genome index from the `genome_index_ch` channel created in the process `1A_prepare_genome_samtools`

- the genome dictionary from the `genome_dict_ch` channel created in the process `1B_prepare_genome_picard`

- the set containing the aligned reads from the `aligned_bam_ch` channel created in the process `2 _rnaseq_mapping_star`

- a set containing the sample id, the split bam file and the split bam index

- specifies the input file names `$genome` and `$bam` to GATK

- specifies the output file names to GATK
