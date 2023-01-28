<div class="formalpara-title">

**Solution**

</div>

``` nextflow
process '2_rnaseq_mapping_star' {

  input:
    path genome from params.genome 
    path genomeDir from genome_dir_ch 
    tuple val(replicateId), path(reads) from reads_ch 

  output:
    tuple val(replicateId), path('Aligned.sortedByCoord.out.bam'), path('Aligned.sortedByCoord.out.bam.bai') into aligned_bam_ch 

  script:
  """
  # ngs-nf-dev Align reads to genome
  STAR --genomeDir $genomeDir \
       --readFilesIn $reads \
       --runThreadN ${task.cpus} \
       --readFilesCommand zcat \
       --outFilterType BySJout \
       --alignSJoverhangMin 8 \
       --alignSJDBoverhangMin 1 \
       --outFilterMismatchNmax 999

  # 2nd pass (improve alignmets using table of splice junctions and create a new index)
  mkdir genomeDir
  STAR --runMode genomeGenerate \
       --genomeDir genomeDir \
       --genomeFastaFiles $genome \
       --sjdbFileChrStartEnd SJ.out.tab \
       --sjdbOverhang 75 \
       --runThreadN ${task.cpus}

  # Final read alignments
  STAR --genomeDir genomeDir \
       --readFilesIn $reads \
       --runThreadN ${task.cpus} \
       --readFilesCommand zcat \
       --outFilterType BySJout \
       --alignSJoverhangMin 8 \
       --alignSJDBoverhangMin 1 \
       --outFilterMismatchNmax 999 \
       --outSAMtype BAM SortedByCoordinate \
       --outSAMattrRGline ID:$replicateId LB:library PL:illumina PU:machine SM:GM12878

  # Index the BAM file
  samtools index Aligned.sortedByCoord.out.bam 
  """
}
```

- the genome fasta file.

- the STAR genome index directory from the `genome_dir_ch` channel created in the process `1C_prepare_star_genome_index`.

- set containing replicate ID and pairs of reads.

- set containing the replicate ID, resulting bam file and bam index.

- line specifying the name of the resulting bam file which is indexed with samtools to create a bam index file (`.bai`).
