The tutorial in this hands-on session shows how to implement a Variant Calling analysis pipeline for RNA-seq data based on GATK best practices and using Nextflow as the pipeline framework.

# Data Description

The input data used to test the pipeline implementation is described below. For the purpose of this project, only a subset of the original data is used for most of the data types.

Genome assembly  
`genome.fa`

The human genome assembly <span class="crg">hg19 (GRCh37)</span> from [GenBank](https://www.ncbi.nlm.nih.gov/assembly/GCA_000001405.1), chromosome 22 only.

RNA-seq reads  
`ENCSR000COQ[12]_[12].fastq.gz`

The RNA-seq data comes from the human <span class="crg">GM12878</span> cell line from whole cell, cytosol and nucleous extraction (see table below).

The libraries are <span class="crg">stranded PE76 Illumina GAIIx</span> RNA-Seq from <span class="crg">rRNA-depleted Poly-A+</span> long RNA (`> 200` nucleotides in size).

Only reads mapped to the [22q11^](http://genome-euro.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr22%3A14700001-25900000&hgsid=221945779_QucOFSFGagd1cn9uVki0TFjrxSBU) locus of the human genome (`chr22:16000000-18000000`) are used.

<table>
<colgroup>
<col style="width: 25%" />
<col style="width: 25%" />
<col style="width: 25%" />
<col style="width: 25%" />
</colgroup>
<tbody>
<tr class="odd">
<td style="text-align: left;"><p>ENCODE ID</p></td>
<td style="text-align: left;"><p>Cellular fraction</p></td>
<td style="text-align: left;"><p>replicate ID</p></td>
<td style="text-align: left;"><p>file names</p></td>
</tr>
<tr class="even">
<td style="text-align: left;"><p><a href="https://www.encodeproject.org/experiments/ENCSR000COQ/">ENCSR000COQ</a></p></td>
<td style="text-align: left;"><p>Whole Cell</p></td>
<td style="text-align: left;"><p>1</p></td>
<td style="text-align: left;"><pre><code>ENCSR000COQ1_1.fastq.gz
ENCSR000COQ1_2.fastq.gz</code></pre></td>
</tr>
<tr class="odd">
<td style="text-align: left;"><p>2</p></td>
<td style="text-align: left;"><pre><code>ENCSR000COQ2_1.fastq.gz
ENCSR000COQ2_2.fastq.gz</code></pre></td>
<td></td>
<td></td>
</tr>
<tr class="even">
<td style="text-align: left;"><p><a href="https://www.encodeproject.org/experiments/ENCSR000CPO/">ENCSR000CPO</a></p></td>
<td style="text-align: left;"><p>Nuclear</p></td>
<td style="text-align: left;"><p>1</p></td>
<td style="text-align: left;"><pre><code>ENCSR000CPO1_1.fastq.gz
ENCSR000CPO1_2.fastq.gz</code></pre></td>
</tr>
<tr class="odd">
<td style="text-align: left;"><p>2</p></td>
<td style="text-align: left;"><pre><code>ENCSR000CPO2_1.fastq.gz
ENCSR000CPO2_2.fastq.gz</code></pre></td>
<td></td>
<td></td>
</tr>
<tr class="even">
<td style="text-align: left;"><p><a href="https://www.encodeproject.org/experiments/ENCSR000COR/">ENCSR000COR</a></p></td>
<td style="text-align: left;"><p>Cytosolic</p></td>
<td style="text-align: left;"><p>1</p></td>
<td style="text-align: left;"><pre><code>ENCSR000COR1_1.fastq.gz
ENCSR000COR1_1.fastq.gz</code></pre></td>
</tr>
<tr class="odd">
<td style="text-align: left;"><p>2</p></td>
<td style="text-align: left;"><pre><code>ENCSR000COR2_1.fastq.gz
ENCSR000COR2_1.fastq.gz</code></pre></td>
<td></td>
<td></td>
</tr>
</tbody>
</table>

"Known" variants  
`known_variants.vcf.gz`

Known variants come from high confident variant calls for <span class="crg">GM12878</span> from the [Illumina Platinum Genomes](https://www.illumina.com/platinumgenomes.html) project. These variant calls were obtained by taking into account pedigree information and the concordance of calls across different methods.

We’re using the subset from chromosome 22 only.

Blacklisted regions  
`blacklist.bed`

Blacklisted regions are regions of the genomes with anomalous coverage. We use regions for the <span class="crg">hg19</span> assembly, taken from the [ENCODE project portal](https://www.encodeproject.org/annotations/ENCSR636HFF/). These regions were identified with DNAse and ChiP-seq samples over ~60 human tissues/cell types, and had a very high ratio of multi-mapping to unique-mapping reads and high variance in mappability.

# Workflow Description

The aim of the pipeline is to process raw RNA-seq data (in FASTQ format) and obtain the list of small variants, SNVs (SNPs and INDELs) for the downstream analysis. The pipeline is based on the [GATK best practices for variant calling with RNAseq data](https://software.broadinstitute.org/gatk/guide/article?id=3891) and includes all major steps. In addition the pipeline includes SNVs postprocessing and quantification for allele specific expression.

Samples processing is done **independently** for **each replicate**. This includes mapping of the reads, splitting at the CIGAR, reassigning mapping qualities and recalibrating base qualities.

Variant calling is done **simultaneously** on bam files from **all replicates**. This allows to improve coverage of genomic regions and obtain more reliable results.

## Software manuals

Documentation for all software used in the workflow can be found at the following links:

- [samtools](http://www.htslib.org/doc/samtools.html)

- [picard `CreateSequenceDictionary`](https://broadinstitute.github.io/picard/command-line-overview.html#CreateSequenceDictionary)

- [STAR](http://labshare.cshl.edu/shares/gingeraslab/www-data/dobin/STAR/STAR.posix/doc/STARmanual.pdf)

- [vcftools](https://vcftools.github.io/man_latest.html)

- [GATK tools](https://software.broadinstitute.org/gatk/gatkdocs/3.6-0/index)

  - [`SplitNCigarReads`](https://software.broadinstitute.org/gatk/gatkdocs/3.6-0/org_broadinstitute_gatk_tools_walkers_rnaseq_SplitNCigarReads.php)

  - [`BaseRecalibrator`](https://software.broadinstitute.org/gatk/gatkdocs/3.6-0/org_broadinstitute_gatk_tools_walkers_bqsr_BaseRecalibrator.php)

  - [`PrintReads`](https://software.broadinstitute.org/gatk/gatkdocs/3.6-0/org_broadinstitute_gatk_tools_walkers_readutils_PrintReads.php)

  - [`HaplotypeCaller`](https://software.broadinstitute.org/gatk/gatkdocs/3.6-0/org_broadinstitute_gatk_tools_walkers_haplotypecaller_HaplotypeCaller.php)

  - [`VariantFiltration`](https://software.broadinstitute.org/gatk/gatkdocs/3.6-0/org_broadinstitute_gatk_tools_walkers_filters_VariantFiltration.php)

  - [`ASEReadCounter`](https://software.broadinstitute.org/gatk/gatkdocs/3.6-0/org_broadinstitute_gatk_tools_walkers_rnaseq_ASEReadCounter.php)

## Pipeline steps

In order to get a general idea of the workflow, all the composing steps, together with the corresponding commands, are explained in the next sections.

### Preparing data

This step prepares input files for the analysis. Genome indexes are created and variants overlapping blacklisted regions are filtered out.

Genome indices with `samtools` and `picard` are produced first. They will be needed for GATK commands such as `Split'N'Trim`:

    samtools faidx genome.fa
    java -jar picard.jar CreateSequenceDictionary R= genome.fa O= genome.dict

Genome index for `STAR`, needed for RNA-seq reads mappings, is created next. Index files are written to the folder `genome_dir` :

    STAR --runMode genomeGenerate \
         --genomeDir genome_dir \
         --genomeFastaFiles genome.fa \
         --runThreadN 4

Variants overlapping blacklisted regions are then filtered in order to reduce false positive calls \[optional\]:

    vcftools --gzvcf known_variants.vcf.gz -c \
             --exclude-bed blacklist.bed \
             --recode | bgzip -c \
             > known_variants.filtered.recode.vcf.gz

### Mapping RNA-seq reads to the reference

To align RNA-seq reads to the genome we’re using STAR 2-pass approach. The first alignment creates a table with splice-junctions that is used to guide final alignments. The alignments at both steps are done with default parameters.

Additional fields with the read groups, libraries and sample information are added into the final bam file at the second mapping step. As a result we do not need to run Picard processing step from GATK best practices.

STAR 1-pass:

     STAR --genomeDir genome_dir \
          --readFilesIn ENCSR000COQ1_1.fastq.gz ENCSR000COQ1_2.fastq.gz \
          --runThreadN 4 \
          --readFilesCommand zcat \
          --outFilterType BySJout \
          --alignSJoverhangMin 8 \
          --alignSJDBoverhangMin 1 \
          --outFilterMismatchNmax 999

Create new genome index using splice-junction table:

    STAR --runMode genomeGenerate \
         --genomeDir genome_dir \
         --genomeFastaFiles genome.fa \
         --sjdbFileChrStartEnd SJ.out.tab \
         --sjdbOverhang 75 \
         --runThreadN 4

STAR 2-pass, final alignments:

      STAR --genomeDir genome_dir \
           --readFilesIn ENCSR000COQ1_1.fastq.gz ENCSR000COQ1_2.fastq.gz \
           --runThreadN 4 \
           --readFilesCommand zcat \
           --outFilterType BySJout \
           --alignSJoverhangMin 8 \
           --alignSJDBoverhangMin 1 \
           --outFilterMismatchNmax 999 \
           --outSAMtype BAM SortedByCoordinate \
           --outSAMattrRGline ID:ENCSR000COQ1 LB:library PL:illumina PU:machine SM:GM12878

Index the resulting bam file:

    samtools index final_alignments.bam

### Split’N'Trim and reassign mapping qualities

The RNA-seq reads overlapping exon-intron junctions can produce false positive variants due to inaccurate splicing. To solve this problem the GATK team recommend to hard-clip any sequence that overlap intronic regions and developed a speciall tool for this purpose: `SplitNCigarReads`. The tool identifies Ns in the CIGAR string of the alignment and split reads at this position so that few new reads are created.

At this step we also reassign mapping qualities to the alignments. This is important because STAR assign the value `255` (high quality) to “unknown” mappings that are meaningless to GATK and to variant calling in general.

This step is done with recommended parameters from the GATK best practices.

    java -jar GenomeAnalysisTK.jar -T SplitNCigarReads \
                      -R genome.fa -I final_alignments.bam \
                      -o split.bam \
                      -rf ReassignOneMappingQuality \
                      -RMQF 255 -RMQT 60 \
                      -U ALLOW_N_CIGAR_READS \
                      --fix_misencoded_quality_scores

### Base Recalibration

The proposed worflow does not include an indel re-alignment step, which is an optional step in the GATK best practices. We excluded that since it is quite time-intensive and does not really improve variant calling.

We instead include a base re-calibration step. This step allows to remove possible systematic errors introduced by the sequencing machine during the assignment of read qualities. To do this, the list of known variants is used as a training set to the machine learning algorithm that models possible errors. Base quality scores are then adjusted based on the obtained results.

    java -jar GenomeAnalysisTK.jar -T BaseRecalibrator \
                      --default_platform illumina \
                      -cov ReadGroupCovariate \
                      -cov QualityScoreCovariate \
                      -cov CycleCovariate \
                      -knownSites known_variants.filtered.recode.vcf.gz\
                      -cov ContextCovariate \
                      -R genome.fa -I split.bam \
                      --downsampling_type NONE \
                      -nct 4 \
                      -o final.rnaseq.grp

      java -jar GenomeAnalysisTK.jar -T PrintReads \
                      -R genome.fa -I split.bam \
                      -BQSR final.rnaseq.grp \
                      -nct 4 \
                      -o final.bam

### Variant Calling and Variant filtering

The variant calling is done on the uniquely aligned reads only in order to reduce the number of false positive variants called:

    (samtools view -H final.bam; samtools view final.bam| grep -w 'NH:i:1') \
      | samtools view -Sb -  > final.uniq.bam

    samtools index final.uniq.bam

For variant calling we’re using the GATK tool `HaplotypeCaller` with default parameters:

    ls final.uniq.bam  > bam.list
    java -jar GenomeAnalysisTK.jar -T HaplotypeCaller \
                      -R genome.fa -I bam.list \
                      -dontUseSoftClippedBases \
                      -stand_call_conf 20.0 \
                      -o output.gatk.vcf.gz

Variant filtering is done as recommended in the GATK best practices:

- keep clusters of at least 3 SNPs that are within a window of 35 bases between them

- estimate strand bias using Fisher’s Exact Test with values \> 30.0 (Phred-scaled p-value)

- use variant call confidence score `QualByDepth` (QD) with values \< 2.0. The QD is the QUAL score normalized by allele depth (AD) for a variant.

<!-- -->

     java -jar GenomeAnalysisTK.jar -T VariantFiltration \
                      -R genome.fa -V output.gatk.vcf.gz \
                      -window 35 -cluster 3 \
                      -filterName FS -filter "FS > 30.0" \
                      -filterName QD -filter "QD < 2.0" \
                      -o final.vcf

### Variant Post-processing

For downstream analysis we’re considering only sites that pass all filters and are covered with at least 8 reads:

    grep -v '#' final.vcf | awk '$7~/PASS/' \
    | perl -ne 'chomp($_); ($dp)=$_=~/DP\\=(\\d+)\\;/; if($dp>=8){print $_."\\n"};' > result.DP8.vcf

Filtered RNA-seq variants are compared with those obtained from DNA sequencing (from Illumina platinum genome project). Variants that are common to these two datasets are "known" SNVs. The ones present only in the RNA-seq cohort only are "novel".

<div class="note">

**Known SNVs** will be used for **allele specific expression** analysis.

**Novel variants** will be used to detect **RNA-editing events**.

</div>

We compare two variants files to detect common and different sites:

    vcftools --vcf result.DP8.vcf --gzdiff known_SNVs.filtered.recode.vcf.gz --diff-site --out commonSNPs

Here we select sites present in both files ("known" SNVs only):

    awk 'BEGIN{OFS="\t"} $4~/B/{print $1,$2,$3}' commonSNPs.diff.sites_in_files  > test.bed

    vcftools --vcf final.vcf --bed test.bed --recode --keep-INFO-all --stdout > known_snps.vcf

Plot a histogram with allele frequency distribution for "known" SNVs:

    grep -v '#'  known_snps.vcf | awk -F '\\t' '{print $10}' \
                   |awk -F ':' '{print $2}'|perl -ne 'chomp($_); \
                   @v=split(/\\,/,$_); if($v[0]!=0 ||$v[1] !=0)\
                   {print  $v[1]/($v[1]+$v[0])."\\n"; }' |awk '$1!=1' \
                   >AF.4R

    gghist.R -i AF.4R -o AF.histogram.pdf

Calculate read counts for each "known" SNVs per allele for allele specific expression analysis:

    java -jar GenomeAnalysisTK.jar -R genome.fa \
                   -T ASEReadCounter \
                   -o ASE.tsv \
                   -I bam.list \
                   -sites known_snps.vcf

# Environment Setup

## Pipeline data

All the files needed for the hands-on activity are stored in the directory shown below:

``` cmd
tree $HOME/environment/hands-on
```

    $HOME/environment/hands-on
    ├── bin
    │   └── gghist.R
    ├── data
    │   ├── blacklist.bed
    │   ├── genome.fa
    │   ├── known_variants.vcf.gz
    │   └── reads
    │       ├── ENCSR000COQ1_1.fastq.gz
    │       ├── ENCSR000COQ1_2.fastq.gz
    │       ├── ENCSR000COQ2_1.fastq.gz
    │       ├── ENCSR000COQ2_2.fastq.gz
    │       ├── ENCSR000COR1_1.fastq.gz
    │       ├── ENCSR000COR1_2.fastq.gz
    │       ├── ENCSR000COR2_1.fastq.gz
    │       ├── ENCSR000COR2_2.fastq.gz
    │       ├── ENCSR000CPO1_1.fastq.gz
    │       ├── ENCSR000CPO1_2.fastq.gz
    │       ├── ENCSR000CPO2_1.fastq.gz
    │       └── ENCSR000CPO2_2.fastq.gz
    ├── nextflow.config
    └── README.md

    3 directories, 19 files

## Pulling the Docker image

Nextflow can pull Docker images at runtime, but let’s just download it manually to see how Docker works:

``` cmd
docker pull cbcrg/callings-with-gatk:latest
```

You should see the progress of the download:

    sha256:93910bf77bc197cb790eca776e42950bc8eff117bdc6e67157295e09c98fc381: Pulling from cbcrg/callings-with-gatk
    915665fee719: Downloading [=============================================>     ] 47.08 MB/51.36 MB
    f332de2321e6: Downloading [===========>                                       ] 41.96 MB/187.8 MB
    1577a6dd9e43: Downloading [===============================>                   ] 46.72 MB/73.45 MB
    7059d9bb5245: Waiting
    71863f70269f: Waiting
    ce2a2879246d: Waiting
    e38ba5d5f9fb: Waiting
    90158da87bb2: Waiting

and the following message when the pull is completed:

    Digest: sha256:93910bf77bc197cb790eca776e42950bc8eff117bdc6e67157295e09c98fc381
    Status: Downloaded newer image for cbcrg/callings-with-gatk:latest

## Script permission

Make sure the following R script has execute permissions:

``` cmd
chmod +x $HOME/environment/hands-on/bin/gghist.R
```

# Pipeline Implementation

## Data preparation

A first step in any pipeline is to prepare the input data. You will find all the data required to run the pipeline in the folder `data` within the `$HOME/environment/hands-on` repository directory.

There are four data inputs that we will use in this tutorial:

1.  **Genome File** (`data/genome.fa`)

    - Human chromosome 22 in FASTA file format

2.  **Read Files** (`data/reads/`)

    - Sample ENCSR000COQ1: 76bp paired-end reads (`ENCSR000COQ1_1.fq.gz` and `ENCSR000COQ1_2.fq.gz`).

3.  **Variants File** (`data/known_variants.vcf.gz`)

    - Known variants, gzipped as a Variant Calling File (VCF) format.

4.  **Blacklist File** (`data/blacklist.bed`)

    - Genomic locations which are known to produce artifacts and spurious variants in Browser Extensible Data (BED) format.

**\***

## Input parameters

We can begin writing the pipeline by creating and editing a text file called `main.nf` from the `$HOME/nf-course/hands-on` repository directory with your favourite text editor. In this example we are using `nano`:

``` cmd
cd $HOME/nf-course/hands-on
nano main.nf
```

Edit this file to specify the input files as script parameters. Using this notation allows you to override them by specifying different values when launching the pipeline execution.

``` nextflow
/*
 * Define the default parameters 
 */

params.genome     = "$baseDir/data/genome.fa" 
params.variants   = "$baseDir/data/known_variants.vcf.gz"
params.blacklist  = "$baseDir/data/blacklist.bed"
params.reads      = "$baseDir/data/reads/ENCSR000COQ1_{1,2}.fastq.gz" 
params.results    = "results" 
params.gatk       = "/opt/broad/GenomeAnalysisTK.jar" 
```

<div class="tip">

You can copy the above text by using the kbd:\[Cmd+C\] keys, then move in the terminal window, open `nano` and paste the above text by using the kbd:\[Cmd+V\] keys shortcut.

</div>

- The `/\*`, `*` and `*/` specify comment lines which are ignored by Nextflow.

- The `baseDir` variable represents the main script path location.

- The `reads` parameter uses a glob pattern to specify the forward (`ENCSR000COQ1_1.fq.gz`) and reverse (`ENCSR000COQ1_2.fq.gz`) reads are pairs of the same sample.

- The `results` parameter is used to specify a directory called `results`.

- The `gatk` parameter specifies the location of the GATK jar file.

Once you have the default parameters in the `main.nf` file, you can save and run the main script for the first time.

<div class="tip">

With `nano` you can save and close the file with kbd:\[Ctrl+O\], then kbd:\[Enter\], followed by kbd:\[Ctrl+X\].

</div>

To run the main script use the following command:

``` cmd
nextflow run main.nf
```

You should see the script execute, print Nextflow version and pipeline revision and then exit.

    N E X T F L O W  ~  version 20.10.0
    Launching `main.nf` [lethal_faggin] - revision: 4c9a5c830c

**\*** === Problem \#1 Great, now we need to define a [channel](https://www.nextflow.io/docs/latest/channel.html) variable to handle the read-pair files. To do that open the `main.nf` file and copy the lines below at the end of the file.

<div class="tip">

In `nano` you can move to the end of the file using kbd:\[Ctrl+W\] and then kbd:\[Ctrl+V\].

</div>

This time you must fill the `BLANK` space with the correct function and parameter.

``` nextflow
/*
 *  Parse the input parameters
 */

reads_ch        = BLANK
GATK            = params.gatk
```

<div class="tip">

Use the [fromFilePairs](https://www.nextflow.io/docs/latest/channel.html#fromfilepairs) channel factory method. The second one, declares a variable named `GATK` specifying the path of the GATK application file.

</div>

Once you think you have data organised, you can again run the pipeline. However this time, we can use the the `-resume` flag.

``` cmd
nextflow run main.nf -resume
```

<div class="tip">

See [here](https://www.nextflow.io/docs/latest/getstarted.html?highlight=resume#modify-and-resume) for more details about using the `resume` option.

</div>

[Solution](#solutions/anxious_advice.adoc#)

**\***

## Process 1A

Now we have our inputs set up we can move onto the processes. In our first process we will create a genome index using [samtools](http://www.htslib.org/).

You should implement a process having the following structure:

Name  
1A_prepare_genome_samtools

Command  
create a genome index for the genome fasta with samtools

Input  
the genome fasta file

Output  
the samtools genome index file

## Problem \#2

Copy the code below and paste it at the end of `main.nf`.

Your aim is to replace `BLANK` placeholder with the the correct variable name of the genome file that you have defined in previous problem.

``` nextflow
/*
 * Process 1A: Create a FASTA genome index with samtools
 */

process '1A_prepare_genome_samtools' { 

  input:
    path genome from BLANK 

  output:
    path "${genome}.fai" into genome_index_ch 

  script:
  """
  samtools faidx ${genome} 
  """
}
```

In plain english, the process could be written as:

- A **process** called 1A_prepare_genome_samtools

- takes as **input** the genome file from `BLANK`

- and creates as **output** a genome index file which goes into channel `genome_index_ch`

- **script**: using samtools create the genome index from the genome file

Now when we run the pipeline, we see that the process 1A is submitted:

``` cmd
nextflow run main.nf -resume
```

    N E X T F L O W  ~  version 20.10.0
    Launching `main.nf` [cranky_bose] - revision: d1df5b7267
    executor >  local (1)
    [cd/47f882] process > 1A_prepare_genome_samtools [100%] 1 of 1 ✔

[Solution](#solutions/busy_building.adoc#)

**\*** == Process 1B

Our first process created the genome index for GATK using samtools. For the next process we must do something very similar, this time creating a genome sequence dictionary using [Picard](https://broadinstitute.github.io/picard/).

You should implement a process having the following structure:

Name  
1B_prepare_genome_picard

Command  
create a genome dictionary for the genome fasta with Picard tools

Input  
the genome fasta file

Output  
the genome dictionary file

## Problem \#3

Fill in the `BLANK` words for both the input and output sections.

Copy the code below and paste it at the end of `main.nf`.

Your aim is to insert the correct input name from into the input step (written as `BLANK`) of the process and run the pipeline.

<div class="tip">

You can choose any channel output name that makes sense to you.

</div>

``` nextflow
/*
 * Process 1B: Create a FASTA genome sequence dictionary with Picard for GATK
 */

process '1B_prepare_genome_picard' {

  input:
    path genome BLANK BLANK

  output:
    path "${genome.baseName}.dict" BLANK BLANK

  script:
  """
  PICARD=`which picard.jar`
  java -jar \$PICARD CreateSequenceDictionary R= $genome O= ${genome.baseName}.dict
  """
}
```

<div class="note">

`.baseName` returns the filename without the file suffix. If `"${genome}"` is `human.fa`, then `"${genome.baseName}.dict"` would be `human.dict`.

</div>

[Solution](#solutions/cooing_clock.adoc#)

**\***

## Process 1C

Next we must create a genome index for the [STAR](https://github.com/alexdobin/STAR) mapping software.

You should implement a process having the following structure:

Name  
1C_prepare_star_genome_index

Command  
create a STAR genome index for the genome fasta

Input  
the genome fasta file

Output  
a directory containing the STAR genome index

## Problem \#4

This is a similar exercise as problem 3, except this time both `input` and `output` lines have been left `BLANK` and must be completed.

``` nextflow
/*
 * Process 1C: Create the genome index file for STAR
 */

process '1C_prepare_star_genome_index' {

  input:
      BLANK_LINE

  output:
      BLANK_LINE

  script:
  """
  mkdir genome_dir

  STAR --runMode genomeGenerate \
       --genomeDir genome_dir \
       --genomeFastaFiles ${genome} \
       --runThreadN ${task.cpus}
  """
}
```

<div class="tip">

The output of the STAR genomeGenerate command is specified here as `genome_dir`.

</div>

[Solution](#solutions/discreet_direction.adoc#)

**\***

## Process 1D

Next on to something a little more tricky. The next process takes two inputs: the variants file and the blacklist file.

It should output a channel named `prepared_vcf_ch` which emitting a tuple of two files.

<div class="note">

In Nextflow, tuples can be defined in the input or output using the [`tuple`](https://www.nextflow.io/docs/latest/process.html#input-of-type-tuple) qualifier.

</div>

You should implement a process having the following structure:

Name  
1D_prepare_vcf_file

Command  
create a filtered and recoded set of variants

Input  
the variants file  
the blacklisted regions file

Output  
a tuple containing the filtered/recoded VCF file and the tab index (TBI) file.

## Problem \#5

You must fill in the two `BLANK_LINES` in the input and the two `BLANK` output files.

``` nextflow
/*
 * Process 1D: Create a file containing the filtered and recoded set of variants
 */

process '1D_prepare_vcf_file' {

  input:
      BLANK_LINE
      BLANK_LINE

  output:
      tuple BLANK, BLANK into prepared_vcf_ch

  script:
  """
  vcftools --gzvcf $variantsFile -c \//
           --exclude-bed ${blacklisted} \//
           --recode | bgzip -c \
           > ${variantsFile.baseName}.filtered.recode.vcf.gz 

  tabix ${variantsFile.baseName}.filtered.recode.vcf.gz 
  """
}
```

- The input variable for the variants file

- The input variable for the blacklist file

- The first of the two output files

- Generates the second output file named `"${variantsFile.baseName}.filtered.recode.vcf.gz.tbi"`

Try run the pipeline from the project directory with:

``` cmd
nextflow run main.nf -resume
```

[Solution](#solutions/expensive_ear.adoc#)

Congratulations! Part 1 is now complete.

**\***

We have all the data prepared and into channels ready for the more serious steps

## Process 2

In this process, for each sample, we align the reads to our genome using the STAR index we created previously.

You should implement a process having the following structure:

Name  
2_rnaseq_mapping_star

Command  
mapping of the RNA-Seq reads using STAR

Input  
the genome fasta file  
the STAR genome index  
a tuple containing the replicate id and paired read files

Output  
a tuple containing replicate id, aligned bam file & aligned bam file index

## Problem \#6

Copy the code below and paste it at the end of `main.nf`.

You must fill in the three `BLANK_LINE` lines in the input and the one `BLANK_LINE` line in the output.

``` nextflow
/*
 * Process 2: Align RNA-Seq reads to the genome with STAR
 */

process '2_rnaseq_mapping_star' {

  input:
      BLANK_LINE
      BLANK_LINE
      BLANK_LINE

  output:
      BLANK_LINE

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

<div class="tip">

The final command produces an bam index which is the full filename with an additional `.bai` suffix.

</div>

[Solution](#solutions/fat_floor.adoc#)

**\***

The next step is a filtering step using GATK. For each sample, we split all the reads that contain N characters in their [CIGAR](http://genome.sph.umich.edu/wiki/SAM#What_is_a_CIGAR.3F) string.

## Process 3

The process creates k+1 new reads (where k is the number of N cigar elements) that correspond to the segments of the original read beside/between the splicing events represented by the Ns in the original CIGAR.

You should implement a process having the following structure:

Name  
3_rnaseq_gatk_splitNcigar

Command  
split reads on Ns in CIGAR string using GATK

Input  
the genome fasta file  
the genome index made with samtools  
the genome dictionary made with picard  
a tuple containing replicate id, aligned bam file and aligned bam file index from the STAR mapping

Output  
a tuple containing the replicate id, the split bam file and the split bam index file

## Problem \#7

Copy the code below and paste it at the end of `main.nf`.

You must fill in the four `BLANK_LINE` lines in the input and the one `BLANK_LINE` line in the output.

<div class="caution">

There is an optional [`tag`](https://www.nextflow.io/docs/latest/process.html#tag) line added to the start of this process. The [`tag`](https://www.nextflow.io/docs/latest/process.html#tag) line allows you to assign a name to a specific task (single execution of a process). This is particularly useful when there are many samples/replicates which pass through the same process.

</div>

``` nextflow
process '3_rnaseq_gatk_splitNcigar' {
  tag OPTIONAL_BLANK

  input:
      BLANK_LINE
      BLANK_LINE
      BLANK_LINE
      BLANK_LINE

  output:
      BLANK_LINE

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
```

<div class="tip">

The GATK command above automatically creates a bam index (`.bai`) of the `split.bam` output file

</div>

<div class="tip">

A `tag` line would also be useful in [???](#Process 2)

</div>

[Solution](#solutions/gentle_garden.adoc#)

**\***

Next we perform a Base Quality Score Recalibration step using GATK.

## Process 4

This step uses GATK to detect systematic errors in the base quality scores, select unique alignments and then index the resulting bam file with samtools. You can find details of the specific GATK BaseRecalibrator parameters [here](https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_bqsr_BaseRecalibrator.php).

You should implement a process having the following structure:

Name  
4_rnaseq_gatk_recalibrate

Command  
recalibrate reads from each replicate using GATK

Input  
the genome fasta file  
the genome index made with samtools  
the genome dictionary made with picard  
a tuple containing replicate id, aligned bam file and aligned bam file index from process 3  
a tuple containing the filtered/recoded VCF file and the tab index (TBI) file from process 1D  

Output  
a tuple containing the sample id, the unique bam file and the unique bam index file

## Problem \#8

Copy the code below and paste it at the end of `main.nf`.

You must fill in the five `BLANK_LINE` lines in the input and the one `BLANK` in the output line.

``` nextflow
process '4_rnaseq_gatk_recalibrate' {
  tag "$replicateId"

  input:
      BLANK_LINE
      BLANK_LINE
      BLANK_LINE
      BLANK_LINE
      BLANK_LINE

  output:
      BLANK into (final_output_ch, bam_for_ASE_ch) 

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

- The files resulting from this process will be used in two downstream processes. If a process is executed more than once, and the downstream channel is used by more than one process, we must duplicate the channel. We can do this using the `into` operator with parenthesis in the output section. See [here](https://www.nextflow.io/docs/latest/operator.html#into) for more information on using `into`.

- The unique bam file

- The index of the unique bam file (bam file name + `.bai`)

[Solution](#solutions/hulking_hospital.adoc#)

**\***

Now we are ready to perform the variant calling with GATK.

## Process 5

This steps call variants with GATK HaplotypeCaller. You can find details of the specific GATK HaplotypeCaller parameters [here](https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_gatk_tools_walkers_haplotypecaller_HaplotypeCaller.php).

You should implement a process having the following structure:

Name  
5_rnaseq_call_variants

Command  
variant calling of each sample using GATK

Input  
the genome fasta file  
the genome index made with samtools  
the genome dictionary made with picard  
a tuple containing replicate id, aligned bam file and aligned bam file index from process 4

Output  
a tuple containing the sample id the resulting variant calling file (vcf)

## Problem \#9

In this problem we will introduce the use of a channel operator in the input section. The [groupTuple](https://www.nextflow.io/docs/latest/operator.html#grouptuple) operator groups together the tuples emitted by a channel which share a common key.

<div class="caution">

Note that in process 4, we used the sampleID (not replicateID) as the first element of the tuple in the output. Now we combine the replicates by grouping them on the sample ID. It follows from this that process 4 is run one time per replicate and process 5 is run one time per sample.

</div>

Fill in the `BLANK_LINE` lines and `BLANK` words as before.

``` nextflow
process '5_rnaseq_call_variants' {
  tag BLANK

  input:
      BLANK_LINE
      BLANK_LINE
      BLANK_LINE
      BLANK from BLANK.groupTuple()

  output:
      BLANK_LINE

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

[Solution](#solutions/imported_iron.adoc#)

**\***

## Processes 6A and 6B

In the final steps we will create processes for Allele-Specific Expression and RNA Editing Analysis.

We must process the VCF result to prepare variants file for allele specific expression (ASE) analysis. We will implement both processes together.

You should implement two processes having the following structure:

Name  
6A_post_process_vcf

Command  
post-process the variant calling file (vcf) of each sample

Input  
tuple containing the sample ID and vcf file  
a tuple containing the filtered/recoded VCF file and the tab index (TBI) file from process 1D  

Output  
a tuple containing the sample id, the variant calling file (vcf) and a file containing common SNPs

<!-- -->

Name  
6B_prepare_vcf_for_ase

Command  
prepare the VCF for allele specific expression (ASE) and generate a figure in R.

Input  
a tuple containing the sample id, the variant calling file (vcf) and a file containing common SNPs

Output  
a tuple containing the sample ID and known SNPs in the sample for ASE  
a figure of the SNPs generated in R as a PDF file

## Problem \#10

Here we introduce the `publishDir` directive. This allows us to specifiy a location for the outputs of the process. See [here](https://www.nextflow.io/docs/latest/process.html#publishdir) for more details.

You must have the output of process 6A become the input of process 6B.

``` nextflow
process '6A_post_process_vcf' {
  tag BLANK
  publishDir "$params.results/$sampleId" 

  input:
      BLANK_LINE
      BLANK_LINE

  output:
      BLANK_LINE

  script:
  '''
  grep -v '#' final.vcf | awk '$7~/PASS/' |perl -ne 'chomp($_); ($dp)=$_=~/DP\\=(\\d+)\\;/; if($dp>=8){print $_."\\n"};' > result.DP8.vcf

  vcftools --vcf result.DP8.vcf --gzdiff filtered.recode.vcf.gz  --diff-site --out commonSNPs
  '''
}


process '6B_prepare_vcf_for_ase' {
  tag BLANK
  publishDir BLANK

  input:
      BLANK_LINE
  output:
      BLANK_LINE
      BLANK_LINE

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
```

- here the output location is specified as a combination of a pipeline parameter and a process input variable

[Solution](#solutions/jumping_jack.adoc#)

**\*** The final step is the GATK ASEReadCounter.

## Problem \#11

We have seen the basics of using processes in Nextflow. Yet one of the features of Nextflow is the operations that can be performed on channels outside of processes. See [here](https://www.nextflow.io/docs/latest/operator.html) for details on the specific operators.

Before we perform the GATK ASEReadCounter process, we must group the data for allele-specific expression. To do this we must combine channels.

The `bam_for_ASE_ch` channel emites tuples having the following structure, holding the final BAM/BAI files:

``` bash
< sample_id, file_bam, file_bai >
```

The `vcf_for_ASE` channel emits tuples having the following structure:

``` bash
< sample_id, output.vcf >
```

In the first operation, the BAMs are grouped together by sample id.

Next, this resulting channel is merged with the VCFs (vcf_for_ASE) having the same sample id.

We must take the merged channel and creates a channel named `grouped_vcf_bam_bai_ch` emitting the following tuples:

``` bash
< sample_id, file_vcf, List[file_bam], List[file_bai] >
```

Your aim is to fill in the `BLANKS` below.

``` nextflow
bam_for_ASE_ch
  .BLANK                            
  .phase(vcf_for_ASE)               
  .map{ left, right ->              
    def sampleId = left[0]          
    def bam = left[1]               
    def bai = left[2]               
    def vcf = right [1]             
    tuple(BLANK, vcf, BLANK, BLANK) 
  }
  .set { grouped_vcf_bam_bai_ch }   
```

- an operator that groups tuples that contain a common first element.

- the phase operator synchronizes the values emitted by two other channels. See [here](https://www.nextflow.io/docs/latest/operator.html?phase#phase) for more details

- the map operator can apply any function to every item on a channel. In this case we take our tuple from the phase operation, define the seperate elements and create a new tuple.

- define `sampleId` to be the first element of left.

- define bam to be the second element of left.

- define bai to be the third element of left.

- define vcf to be the first element of right.

- create a new tuple made of four elements

- rename the resulting as `grouped_vcf_bam_bai_ch`

<div class="caution">

`left` and `right` above are arbitary names. From the phase operator documentation, we see that phase returns pairs of items. So here `left` originates from contents of the `bam_for_ASE_ch` channel and `right` originates from the contents of `vcf_for_ASE` channel.

</div>

[Solution](#solutions/kind_koala.adoc#)

**\***

## Process 6C

Now we are ready for the final process.

You should implement a process having the following structure:

Name  
6C_ASE_knownSNPs

Command  
calculate allele counts at a set of positions with GATK tools

Input  
genome fasta file  
genome index file from samtools  
genome dictionary file  
the \`grouped_vcf_bam_bai_ch\`channel

Output  
the allele specific expression file (`ASE.tsv`)

## Problem \#12

You should construct the process and run the pipeline in its entirety.

``` nextflow
  echo "${bam.join('\n')}" > bam.list

  java -jar $GATK -R ${genome} \
                  -T ASEReadCounter \
                  -o ASE.tsv \
                  -I bam.list \
                  -sites ${vcf}
```

[Solution](#solutions/laughing_lynx.adoc#)

Congratulations! If you made it this far you now have all the basics to create your own Nextflow workflows.

**\***

# Results overview

For each processed sample the pipeline stores results into a folder named after the sample identifier. These folders are created in the directory specified as a parameter in `params.results`.

Result files for this workshop can be found in the folder `results` within the current folder. There you should see a directory called `ENCSR000COQ/` containing the following files:

Variant calls  
`final.vcf`

This file contains all somatic variants (SNVs) called from RNAseq data. You will see variants that pass all filters, with the `PASS` keyword in the <span class="red">7th</span> field of the vcf file (`filter status`), and also those that did not pass one or more filters.

`commonSNPs.diff.sites_in_files`

Tab-separated file with comparison between variants obtained from RNAseq and "known" variants from DNA.

The file is sorted by genomic position and contains 8 fields:

<table>
<colgroup>
<col style="width: 33%" />
<col style="width: 33%" />
<col style="width: 33%" />
</colgroup>
<tbody>
<tr class="odd">
<td style="text-align: left;"><p><span class="red">1</span></p></td>
<td style="text-align: center;"><pre><code> CHROM   </code></pre></td>
<td style="text-align: left;"><p>chromosome name;</p></td>
</tr>
<tr class="even">
<td style="text-align: left;"><p><span class="red">2</span></p></td>
<td style="text-align: center;"><pre><code> POS1    </code></pre></td>
<td style="text-align: left;"><p>position of the SNV in file #1 (RNAseq data);</p></td>
</tr>
<tr class="odd">
<td style="text-align: left;"><p><span class="red">3</span></p></td>
<td style="text-align: center;"><pre><code> POS2    </code></pre></td>
<td style="text-align: left;"><p>position of SNV in file #2 (DNA "known" variants);</p></td>
</tr>
<tr class="even">
<td style="text-align: left;"><p><span class="red">4</span></p></td>
<td style="text-align: center;"><pre><code> IN_FILE </code></pre></td>
<td style="text-align: left;"><p>flag whether SNV is present in the file #1 <em>1</em>, in the file #2 <em>2</em>, or in both files <em>B</em>;</p></td>
</tr>
<tr class="odd">
<td style="text-align: left;"><p><span class="red">5</span></p></td>
<td style="text-align: center;"><pre><code> REF1    </code></pre></td>
<td style="text-align: left;"><p>reference sequence in the file 1;</p></td>
</tr>
<tr class="even">
<td style="text-align: left;"><p><span class="red">6</span></p></td>
<td style="text-align: center;"><pre><code> REF2    </code></pre></td>
<td style="text-align: left;"><p>reference sequence in the file 2;</p></td>
</tr>
<tr class="odd">
<td style="text-align: left;"><p><span class="red">7</span></p></td>
<td style="text-align: center;"><pre><code> ALT1    </code></pre></td>
<td style="text-align: left;"><p>alternative sequence in the file 1;</p></td>
</tr>
<tr class="even">
<td style="text-align: left;"><p><span class="red">8</span></p></td>
<td style="text-align: center;"><pre><code> ALT2    </code></pre></td>
<td style="text-align: left;"><p>alternative sequence in the file 2</p></td>
</tr>
</tbody>
</table>

`known_snps.vcf`

Variants that are common to RNAseq and "known" variants from DNA.

Allele specific expression quantification  
`ASE.tsv`

Tab-separated file with allele counts at common SNVs positions (only SNVs from the file `known_snps.vcf`)

The file is sorted by coordinates and contains 13 fields:

<table>
<colgroup>
<col style="width: 33%" />
<col style="width: 33%" />
<col style="width: 33%" />
</colgroup>
<tbody>
<tr class="odd">
<td style="text-align: left;"><p><span class="red">1</span></p></td>
<td style="text-align: center;"><pre><code> contig        </code></pre></td>
<td style="text-align: left;"><p>contig, scaffold or chromosome name of the variant</p></td>
</tr>
<tr class="even">
<td style="text-align: left;"><p><span class="red">2</span></p></td>
<td style="text-align: center;"><pre><code> position      </code></pre></td>
<td style="text-align: left;"><p>position of the variant</p></td>
</tr>
<tr class="odd">
<td style="text-align: left;"><p><span class="red">3</span></p></td>
<td style="text-align: center;"><pre><code> variant ID    </code></pre></td>
<td style="text-align: left;"><p>variant ID in the dbSNP</p></td>
</tr>
<tr class="even">
<td style="text-align: left;"><p><span class="red">4</span></p></td>
<td style="text-align: center;"><pre><code> refAllele     </code></pre></td>
<td style="text-align: left;"><p>reference allele sequence</p></td>
</tr>
<tr class="odd">
<td style="text-align: left;"><p><span class="red">5</span></p></td>
<td style="text-align: center;"><pre><code> altAllele     </code></pre></td>
<td style="text-align: left;"><p>alternate allele sequence</p></td>
</tr>
<tr class="even">
<td style="text-align: left;"><p><span class="red">6</span></p></td>
<td style="text-align: center;"><pre><code> refCount      </code></pre></td>
<td style="text-align: left;"><p>number of reads that support the reference allele</p></td>
</tr>
<tr class="odd">
<td style="text-align: left;"><p><span class="red">7</span></p></td>
<td style="text-align: center;"><pre><code> altCount      </code></pre></td>
<td style="text-align: left;"><p>number of reads that support the alternate allele</p></td>
</tr>
<tr class="even">
<td style="text-align: left;"><p><span class="red">8</span></p></td>
<td style="text-align: center;"><pre><code> totalCount    </code></pre></td>
<td style="text-align: left;"><p>total number of reads at the site that support both reference and alternate allele and any other alleles present at the site</p></td>
</tr>
<tr class="odd">
<td style="text-align: left;"><p><span class="red">9</span></p></td>
<td style="text-align: center;"><pre><code> lowMAPQDepth  </code></pre></td>
<td style="text-align: left;"><p>number of reads that have low mapping quality</p></td>
</tr>
<tr class="even">
<td style="text-align: left;"><p><span class="red">10</span></p></td>
<td style="text-align: center;"><pre><code> lowBaseQDepth </code></pre></td>
<td style="text-align: left;"><p>number of reads that have low base quality</p></td>
</tr>
<tr class="odd">
<td style="text-align: left;"><p><span class="red">11</span></p></td>
<td style="text-align: center;"><pre><code> rawDepth      </code></pre></td>
<td style="text-align: left;"><p>total number of reads at the site that support both reference and alternate allele and any other alleles present at the site</p></td>
</tr>
<tr class="even">
<td style="text-align: left;"><p><span class="red">12</span></p></td>
<td style="text-align: center;"><pre><code> otherBases    </code></pre></td>
<td style="text-align: left;"><p>number of reads that support bases other than reference and alternate bases</p></td>
</tr>
<tr class="odd">
<td style="text-align: left;"><p><span class="red">13</span></p></td>
<td style="text-align: center;"><pre><code> improperPairs </code></pre></td>
<td style="text-align: left;"><p>number of reads that have malformed pairs</p></td>
</tr>
</tbody>
</table>

Allele frequency histogram  
`AF.histogram.pdf`

This file contains a histogram plot of allele frequency for SNVs common to RNA-seq and "known" variants from DNA.

## Bonus step

Until now the pipeline has been executed using just a single sample (`ENCSR000COQ1`).

Now we can re-execute the pipeline specifying a large set of samples by using the command shown below:

``` cmd
nextflow run main.nf -resume --reads 'data/reads/ENCSR000C*_{1,2}.fastq.gz'
```

It will print an output similar to the one below:

    N E X T F L O W  ~  version 20.10.0
    Launching `main.nf` [hungry_wing] - revision: a6359031a1
    executor >  local (27)
    [cd/47f882] process > 1A_prepare_genome_samtools               [100%] 1 of 1, cached: 1 ✔
    [5f/216ba8] process > 1B_prepare_genome_picard                 [100%] 1 of 1, cached: 1 ✔
    [76/5fdc20] process > 1C_prepare_star_genome_index             [100%] 1 of 1, cached: 1 ✔
    [19/f8842c] process > 1D_prepare_vcf_file                      [100%] 1 of 1, cached: 1 ✔
    [f1/d66ba8] process > 2_rnaseq_mapping_star (6)                [100%] 6 of 6, cached: 1 ✔
    [74/c0f3a3] process > 3_rnaseq_gatk_splitNcigar (ENCSR000CPO2) [100%] 6 of 6, cached: 1 ✔
    [b6/59d9f7] process > 4_rnaseq_gatk_recalibrate (ENCSR000CPO2) [100%] 6 of 6, cached: 1 ✔
    [22/4a07fa] process > 5_rnaseq_call_variants (ENCSR000CPO)     [100%] 3 of 3 ✔
    [1a/c68bfe] process > 6A_post_process_vcf (ENCSR000CPO)        [100%] 3 of 3 ✔
    [dc/e58d02] process > 6B_prepare_vcf_for_ase (ENCSR000CPO)     [100%] 3 of 3 ✔
    [2a/0e4e7b] process > 6C_ASE_knownSNPs (ENCSR000CPO)           [100%] 3 of 3 ✔

You can notice that this time the pipeline spawns the execution of more tasks because three samples have been provided instead of one.

This shows the ability of Nextflow to implicitly handle multiple parallel task executions depending on the specified pipeline input dataset.

A fully functional version of this pipeline is available at the following GitHub repository: [CalliNGS-NF](https://github.com/CRG-CNAG/CalliNGS-NF).
