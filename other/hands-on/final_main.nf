#!/usr/bin/env nextflow

/*
 * Define the default parameters
 */

params.genome     = "${projectDir}/data/genome.fa"
params.variants   = "${projectDir}/data/known_variants.vcf.gz"
params.blacklist  = "${projectDir}/data/blacklist.bed"
params.reads      = "${projectDir}/data/reads/ENCSR000C*_{1,2}.fastq.gz"
params.results    = "results"

/*
 * Process 1A: Create a FASTA genome index with samtools
 */

process prepare_genome_samtools {
    container 'quay.io/biocontainers/samtools:1.3.1--h0cf4675_11'

    input:
    path genome

    output:
    path "${genome}.fai"

    script:
    """
    samtools faidx ${genome}
    """
}

/*
 * Process 1B: Create a FASTA genome sequence dictionary with Picard for GATK
 */

process prepare_genome_picard {
    container 'quay.io/biocontainers/picard:1.141--hdfd78af_6'

    input:
    path genome

    output:
    path "${genome.baseName}.dict"

    script:
    """
    picard CreateSequenceDictionary R= ${genome} O= ${genome.baseName}.dict
    """
}

/*
* Process 1C: Create the genome index file for STAR
*/

process prepare_star_genome_index {
    container 'quay.io/biocontainers/star:2.7.10b--h6b7c446_1'

    input:
    path genome

    output:
    path 'genome_dir'

    script:
    """
    mkdir -p genome_dir

    STAR --runMode genomeGenerate \
        --genomeDir genome_dir \
        --genomeFastaFiles ${genome} \
        --runThreadN ${task.cpus}
    """
}

/*
 * Process 1D: Create a file containing the filtered and recoded set of variants
 */

process prepare_vcf_file {
    container 'quay.io/biocontainers/mulled-v2-b9358559e3ae3b9d7d8dbf1f401ae1fcaf757de3:ac05763cf181a5070c2fdb9bb5461f8d08f7b93b-0'

    input:
    path variantsFile
    path blacklisted

    output:
    tuple path("${variantsFile.baseName}.filtered.recode.vcf.gz"),
        path("${variantsFile.baseName}.filtered.recode.vcf.gz.tbi")

    script:
    """
    vcftools --gzvcf ${variantsFile} -c \
            --exclude-bed ${blacklisted} \
            --recode | bgzip -c \
            > ${variantsFile.baseName}.filtered.recode.vcf.gz

    tabix ${variantsFile.baseName}.filtered.recode.vcf.gz
    """
}

/*
 * Process 2: Align RNA-Seq reads to the genome with STAR
 */

process rnaseq_mapping_star {
    container 'quay.io/biocontainers/mulled-v2-52f8f283e3c401243cee4ee45f80122fbf6df3bb:e3bc54570927dc255f0e580cba1789b64690d611-0'

    input:
    path genome
    path genomeDir
    tuple val(replicateId), path(reads)

    output:
    tuple val(replicateId),
        path('Aligned.sortedByCoord.out.bam'),
        path('Aligned.sortedByCoord.out.bam.bai')

    script:
    """
    # ngs-nf-dev Align reads to genome
    STAR --genomeDir ${genomeDir} \
        --readFilesIn ${reads} \
        --runThreadN ${task.cpus} \
        --readFilesCommand zcat \
        --outFilterType BySJout \
        --alignSJoverhangMin 8 \
        --alignSJDBoverhangMin 1 \
        --outFilterMismatchNmax 999

    # 2nd pass (improve alignments using table of splice
    # junctions and create a new index)
    mkdir -p genomeDir
    STAR --runMode genomeGenerate \
        --genomeDir genomeDir \
        --genomeFastaFiles ${genome} \
        --sjdbFileChrStartEnd SJ.out.tab \
        --sjdbOverhang 75 \
        --runThreadN ${task.cpus}

    # Final read alignments
    STAR --genomeDir genomeDir \
        --readFilesIn ${reads} \
        --runThreadN ${task.cpus} \
        --readFilesCommand zcat \
        --outFilterType BySJout \
        --alignSJoverhangMin 8 \
        --alignSJDBoverhangMin 1 \
        --outFilterMismatchNmax 999 \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMattrRGline ID:${replicateId} LB:library PL:illumina \
                            PU:machine SM:GM12878

    # Index the BAM file
    samtools index Aligned.sortedByCoord.out.bam
    """
}

/*
 * Process 3: GATK Split on N
 */

process rnaseq_gatk_splitNcigar {
    container 'quay.io/broadinstitute/gotc-prod-gatk:1.0.0-4.1.8.0-1626439571'
    tag "${replicateId}"

    input:
    path genome
    path index
    path genome_dict
    tuple val(replicateId),
        path(bam),
        path(bai)

    output:
    tuple val(replicateId), path('split.bam'), path('split.bai')

    script:
    """
    # SplitNCigarReads and reassign mapping qualities
    java -jar /usr/gitc/GATK35.jar -T SplitNCigarReads \
                                -R ${genome} -I ${bam} \
                                -o split.bam \
                                -rf ReassignOneMappingQuality \
                                -RMQF 255 -RMQT 60 \
                                -U ALLOW_N_CIGAR_READS \
                                --fix_misencoded_quality_scores

    """
}

/*
 * Process 4: GATK Recalibrate
 */

process rnaseq_gatk_recalibrate {
    container 'quay.io/biocontainers/mulled-v2-aa1d7bddaee5eb6c4cbab18f8a072e3ea7ec3969:f963c36fd770e89d267eeaa27cad95c1c3dbe660-0'
    tag "${replicateId}"

    input:
    path genome
    path index
    path dict
    tuple val(replicateId), path(bam), path(bai)
    tuple path(prepared_variants_file), path(prepared_variants_file_index)

    output:
    tuple val(sampleId),
        path("${replicateId}.final.uniq.bam"),
        path("${replicateId}.final.uniq.bam.bai")

    script:
    sampleId = replicateId.replaceAll(/[12]$/,'')
    """
    # Indel Realignment and Base Recalibration
    gatk3 -T BaseRecalibrator \
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

    gatk3 -T PrintReads \
        -R ${genome} -I ${bam} \
        -BQSR final.rnaseq.grp \
        -nct ${task.cpus} \
        -o final.bam

    # Select only unique alignments, no multimaps
    (samtools view -H final.bam; samtools view final.bam | \
    grep -w 'NH:i:1') | samtools view -Sb -  > ${replicateId}.final.uniq.bam

    # Index BAM files
    samtools index ${replicateId}.final.uniq.bam
    """
}


/*
* Process 5: GATK Variant Calling
*/

process rnaseq_call_variants {
    container 'quay.io/broadinstitute/gotc-prod-gatk:1.0.0-4.1.8.0-1626439571'
    tag "${sampleId}"

    input:
    path genome
    path index
    path dict
    tuple val(sampleId), path(bam), path(bai)

    output:
    tuple val(sampleId), path('final.vcf')

    script:
    """
    echo "${bam.join('\n')}" > bam.list

    # Variant calling
    java -jar /usr/gitc/GATK35.jar -T HaplotypeCaller \
                    -R ${genome} -I bam.list \
                    -dontUseSoftClippedBases \
                    -stand_call_conf 20.0 \
                    -o output.gatk.vcf.gz

    # Variant filtering
    java -jar /usr/gitc/GATK35.jar -T VariantFiltration \
                                -R ${genome} -V output.gatk.vcf.gz \
                                -window 35 -cluster 3 \
                                -filterName FS -filter "FS > 30.0" \
                                -filterName QD -filter "QD < 2.0" \
                                -o final.vcf
    """
}

/*
 * Process 6: ASE & RNA Editing
 */

process post_process_vcf {
    container 'quay.io/biocontainers/mulled-v2-b9358559e3ae3b9d7d8dbf1f401ae1fcaf757de3:ac05763cf181a5070c2fdb9bb5461f8d08f7b93b-0'
    tag "${sampleId}"
    publishDir "${params.results}/${sampleId}"

    input:
    tuple val(sampleId), path('final.vcf')
    tuple path('filtered.recode.vcf.gz'),
        path('filtered.recode.vcf.gz.tbi')

    output:
    tuple val(sampleId),
        path('final.vcf'),
        path('commonSNPs.diff.sites_in_files')

    script:
    '''
    grep -v '#' final.vcf | awk '$7~/PASS/' | perl -ne 'chomp($_); \
                            ($dp)=$_=~/DP\\=(\\d+)\\;/; \
                            if($dp>=8){print $_."\\n"};' > result.DP8.vcf

    vcftools --vcf result.DP8.vcf --gzdiff filtered.recode.vcf.gz  --diff-site \
            --out commonSNPs
    '''
}

process prepare_vcf_for_ase {
    container 'cbcrg/callings-with-gatk:latest'
    tag "${sampleId}"
    publishDir "${params.results}/${sampleId}"

    input:
    tuple val(sampleId),
        path('final.vcf'),
        path('commonSNPs.diff.sites_in_files')

    output:
    tuple val(sampleId), path('known_snps.vcf'), emit: vcf_for_ASE
    path 'AF.histogram.pdf'                    , emit: gghist_pdfs

    script:
    '''
    awk 'BEGIN{OFS="\t"} $4~/B/{print $1,$2,$3}' commonSNPs.diff.sites_in_files  > test.bed

    vcftools --vcf final.vcf --bed test.bed --recode --keep-INFO-all \
            --stdout > known_snps.vcf

    grep -v '#' known_snps.vcf | awk -F '\\t' '{print $10}' \
                | awk -F ':' '{print $2}' | perl -ne 'chomp($_); \
                @v=split(/\\,/,$_); if($v[0]!=0 ||$v[1] !=0)\
                {print  $v[1]/($v[1]+$v[0])."\\n"; }' | awk '$1!=1' \
                > AF.4R

    gghist.R -i AF.4R -o AF.histogram.pdf
    '''
}

/*
 * Process 7: Allele-Specific Expression analysis with GATK ASEReadCounter
 */

process ASE_knownSNPs {
    container 'quay.io/broadinstitute/gotc-prod-gatk:1.0.0-4.1.8.0-1626439571'
    tag "${sampleId}"
    publishDir "${params.results}/${sampleId}"

    input:
    path genome
    path index
    path dict
    tuple val(sampleId), path(vcf), path(bam), path(bai)

    output:
    path "ASE.tsv"

    script:
    """
    echo "${bam.join('\n')}" > bam.list

    java -jar /usr/gitc/GATK35.jar -R ${genome} \
                                -T ASEReadCounter \
                                -o ASE.tsv \
                                -I bam.list \
                                -sites ${vcf}
    """
}


workflow {
    reads_ch =  channel.fromFilePairs(params.reads)

    prepare_genome_samtools(params.genome)
    prepare_genome_picard(params.genome)
    prepare_star_genome_index(params.genome)
    prepare_vcf_file(params.variants, params.blacklist)

    rnaseq_mapping_star(params.genome, prepare_star_genome_index.out, reads_ch)

    rnaseq_gatk_splitNcigar(params.genome,
                            prepare_genome_samtools.out,
                            prepare_genome_picard.out,
                            rnaseq_mapping_star.out)

    rnaseq_gatk_recalibrate(params.genome,
                        prepare_genome_samtools.out,
                        prepare_genome_picard.out,
                        rnaseq_gatk_splitNcigar.out,
                        prepare_vcf_file.out)

    // New channel to aggregate bam from different replicates into sample level.
    rnaseq_gatk_recalibrate.out
        | groupTuple
        | set { recalibrated_samples }

    rnaseq_call_variants(params.genome,
                        prepare_genome_samtools.out,
                        prepare_genome_picard.out,
                        recalibrated_samples)

    post_process_vcf(rnaseq_call_variants.out,
                        prepare_vcf_file.out)

    prepare_vcf_for_ase(post_process_vcf.out)

    recalibrated_samples
        .join(prepare_vcf_for_ase.out.vcf_for_ASE)
        .map { meta, bams, bais, vcf -> [meta, vcf, bams, bais] }
        .set { grouped_vcf_bam_bai_ch }

    ASE_knownSNPs(params.genome,
                prepare_genome_samtools.out,
                prepare_genome_picard.out,
                grouped_vcf_bam_bai_ch)
}
