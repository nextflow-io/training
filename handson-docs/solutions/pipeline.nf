params.genome     = "$baseDir/data/genome.fa" 
params.variants   = "$baseDir/data/known_variants.vcf.gz"
params.blacklist  = "$baseDir/data/blacklist.bed"
params.reads      = "$baseDir/data/reads/ENCSR000COQ1_{1,2}.fastq.gz" 
params.results    = "results" 
params.gatk       = "/opt/broad/GenomeAnalysisTK.jar"


reads_ch        =  Channel.fromFilePairs(params.reads) 
GATK            =  params.gatk 

process '1A_prepare_genome_samtools' {

  input:
    path genome from params.genome 

  output:
    path "${genome}.fai" into genome_index_ch

  script:
  """
  samtools faidx ${genome}
  """
}

/*
 * Process 1B: Create a FASTA genome sequence dictionary with Picard for GATK
 */

process '1B_prepare_genome_picard' {

  input:
    path genome from params.genome 

  output:
    path "${genome.baseName}.dict" into genome_dict_ch 

  script:
  """
  PICARD=`which picard.jar`
  java -jar \$PICARD CreateSequenceDictionary R= $genome O= ${genome.baseName}.dict
  """
}

process '1C_prepare_star_genome_index' {

  input:
    path genome from params.genome 

  output:
    path 'genome_dir' into genome_dir_ch 

  script:
  """
  mkdir genome_dir 

  STAR --runMode genomeGenerate \
       --genomeDir genome_dir \
       --genomeFastaFiles ${genome} \
       --runThreadN ${task.cpus}
  """
}

process '1D_prepare_vcf_file' {

  input:
    path variantsFile from params.variants 
    path blacklisted from params.blacklist 

  output:
    tuple path("${variantsFile.baseName}.filtered.recode.vcf.gz"), \
          path("${variantsFile.baseName}.filtered.recode.vcf.gz.tbi") into prepared_vcf_ch 

  script:
  """
  vcftools --gzvcf $variantsFile -c \
           --exclude-bed ${blacklisted} \
           --recode | bgzip -c \
           > ${variantsFile.baseName}.filtered.recode.vcf.gz 

  tabix ${variantsFile.baseName}.filtered.recode.vcf.gz 
  """
}

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
                  -R $genome -I $bam \
                  -o split.bam \
                  -rf ReassignOneMappingQuality \
                  -RMQF 255 -RMQT 60 \
                  -U ALLOW_N_CIGAR_READS \
                  --fix_misencoded_quality_scores

  """
}

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

process '6A_post_process_vcf' {
  tag "$sampleId"
  publishDir "$params.results/$sampleId"

  input:
    tuple val(sampleId), path('final.vcf') from vcf_files
    tuple path('filtered.recode.vcf.gz'), path('filtered.recode.vcf.gz.tbi') from prepared_vcf_ch
  output:
    tuple val(sampleId), path('final.vcf'), path('commonSNPs.diff.sites_in_files') into vcf_and_snps_ch

  script:
  '''
  grep -v '#' final.vcf | awk '$7~/PASS/' |perl -ne 'chomp($_); ($dp)=$_=~/DP\\=(\\d+)\\;/; if($dp>=8){print $_."\\n"};' > result.DP8.vcf

  vcftools --vcf result.DP8.vcf --gzdiff filtered.recode.vcf.gz  --diff-site --out commonSNPs
  '''
}

process '6B_prepare_vcf_for_ase' {
  tag "$sampleId"
  publishDir "$params.results/$sampleId"

  input:
    tuple val(sampleId), path('final.vcf'), path('commonSNPs.diff.sites_in_files') from vcf_and_snps_ch
  output:
    tuple val(sampleId), path('known_snps.vcf') into vcf_for_ASE
    path 'AF.histogram.pdf' into gghist_pdfs

  script:
  '''
  awk 'BEGIN{OFS="\t"} $4~/B/{print $1,$2,$3}' commonSNPs.diff.sites_in_files  > test.bed

  vcftools --vcf final.vcf --bed test.bed --recode --keep-INFO-all --stdout > known_snps.vcf

  grep -v '#'  known_snps.vcf | awk -F '\\t' '{print $10}' \
               |awk -F ':' '{print $2}'|perl -ne 'chomp($_); \
               @v=split(/\\,/,$_); if($v[0]!=0 ||$v[1] !=0)\
               {print  $v[1]/($v[1]+$v[0])."\\n"; }' |awk '$1!=1' \
               >AF.4R

  gghist.R -i AF.4R -o AF.histogram.pdf
  '''
}

bam_for_ASE_ch
  .groupTuple()
  .phase(vcf_for_ASE)
  .map{ left, right ->
    def sampleId = left[0]
    def bam = left[1]
    def bai = left[2]
    def vcf = right[1]
    tuple(sampleId, vcf, bam, bai)
  }
  .set { grouped_vcf_bam_bai_ch }


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

