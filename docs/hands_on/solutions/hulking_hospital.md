<div class="formalpara-title">

**Solution**

</div>

``` nextflow
process '4_rnaseq_gatk_recalibrate' {
  tag "$replicateId"

  input:
    path genome from params.genome 
    path index from genome_index_ch 
    path dict from genome_dict_ch 
    tuple val(replicateId), path(bam), path(bai) from splitted_bam_ch 
    tuple path(prepared_variants_file), path(prepared_variants_file_index) from prepared_vcf_ch 

  output:
    tuple val(sampleId), path("${replicateId}.final.uniq.bam"), path("${replicateId}.final.uniq.bam.bai") into (final_output_ch, bam_for_ASE_ch) 

  script:
  sampleId = replicateId.replaceAll(/[12]$/,'')
  """
  # Indel Realignment and Base Recalibration
  java -jar $GATK -T BaseRecalibrator \
                  --default_platform illumina \
                  -cov ReadGroupCovariate \
                  -cov QualityScoreCovariate \
                  -cov CycleCovariate \
                  -knownSites ${prepared_variants_file} \
                  -cov ContextCovariate \
                  -R ${genome} -I ${bam} \
                  --downsampling_type NONE \
                  -nct ${task.cpus} \
                  -o final.rnaseq.grp

  java -jar $GATK -T PrintReads \
                  -R ${genome} -I ${bam} \
                  -BQSR final.rnaseq.grp \
                  -nct ${task.cpus} \
                  -o final.bam

  # Select only unique alignments, no multimaps
  (samtools view -H final.bam; samtools view final.bam| grep -w 'NH:i:1') \
  |samtools view -Sb -  > ${replicateId}.final.uniq.bam 

  # Index BAM files
  samtools index ${replicateId}.final.uniq.bam
  """
}
```

- the genome fasta file.

- the genome index from the `genome_index_ch` channel created in the process `1A_prepare_genome_samtools`.

- the genome dictionary from the `genome_dict_ch` channel created in the process `1B_prepare_genome_picard`.

- the set containing the split reads from the `splitted_bam_ch` channel created in the process `3_rnaseq_gatk_splitNcigar`.

- the set containing the filtered/recoded VCF file and the tab index (TBI) file from the `prepared_vcf_ch` channel created in the process `1D_prepare_vcf_file`.

- the set containing the replicate id, the unique bam file and the unique bam index file which goes into two channels.

- line specifying the filename of the output bam file
