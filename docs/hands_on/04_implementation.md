# Pipeline Implementation

## Data preparation

A first step in any pipeline is to prepare the input data. You will find all the data required to run the pipeline in the folder `data` within the `/workspace/gitpod/hands-on` repository directory.

There are four data inputs that we will use in this tutorial:

1. **Genome File** (`data/genome.fa`)
    - Human chromosome 22 in FASTA file format
2. **Read Files** (`data/reads/`)
    - Sample ENCSR000COQ1: 76bp paired-end reads (`ENCSR000COQ1_1.fq.gz` and `ENCSR000COQ1_2.fq.gz`).
3. **Variants File** (`data/known_variants.vcf.gz`)
    - Known variants, gzipped as a Variant Calling File (VCF) format.
4. **Blacklist File** (`data/blacklist.bed`)
    - Genomic locations which are known to produce artifacts and spurious variants in Browser Extensible Data (BED) format.

## Input parameters

We can begin writing the pipeline by creating and editing a text file called `main.nf` from the `/workspace/gitpod/hands-on` repository directory with your favourite text editor. In this example we are using `code`:

```bash
cd /workspace/gitpod/hands-on
code main.nf
```

Edit this file to specify the input files as script parameters. Using this notation allows you to override them by specifying different values when launching the pipeline execution.

!!! info

    Click the :material-plus-circle: icons in the code for explanations.

```groovy linenums="1"
/*
 * Define the default parameters (1)
 */

params.genome     = "$baseDir/data/genome.fa" // (2)!
params.variants   = "$baseDir/data/known_variants.vcf.gz"
params.blacklist  = "$baseDir/data/blacklist.bed"
params.reads      = "$baseDir/data/reads/ENCSR000COQ1_{1,2}.fastq.gz" // (3)!
params.results    = "results" // (4)!
params.gatk       = "/opt/broad/GenomeAnalysisTK.jar" // (5)!
```

1. The `/*`, `*` and `*/` specify comment lines which are ignored by Nextflow.
2. The `baseDir` variable represents the main script path location.
3. The `reads` parameter uses a glob pattern to specify the forward (`ENCSR000COQ1_1.fq.gz`) and reverse (`ENCSR000COQ1_2.fq.gz`) reads are pairs of the same sample.
4. The `results` parameter is used to specify a directory called `results`.
5. The `gatk` parameter specifies the location of the GATK jar file. This is a path within the container.

!!! tip

    You can copy the above text (:material-content-copy: top right or ++cmd+c++), then move in the terminal window, open `code` and paste using the keyboard ++cmd+v++ shortcut.

Once you have the default parameters in the `main.nf` file, you can save and run the main script for the first time.

!!! tip

    With `code` you can save and close the file with ++ctrl+o++, then ++enter++, followed by ++ctrl+x++.

To run the main script use the following command:

```bash
nextflow run main.nf
```

You should see the script execute, print Nextflow version and pipeline revision and then exit.

```console
N E X T F L O W  ~  version 23.04.1
Launching `main.nf` [lethal_faggin] - revision: 4c9a5c830c
```

!!! exercise "Problem #1"

    Great, now we need to define a [channel](https://www.nextflow.io/docs/latest/channel.html) variable to handle the read-pair files. To do that open the `main.nf` file and copy the lines below at the end of the file.

    !!! tip

        In `code` you can move to the end of the file using ++ctrl+w++ and then ++ctrl+v++.

    This time you must fill the `BLANK` space with the correct function and parameter.

    ```groovy linenums="1" hl_lines="5"
    /*
     *  Parse the input parameters
     */
    workflow {
        reads_ch = BLANK
    }
    ```

    !!! tip

        Use the [fromFilePairs](https://www.nextflow.io/docs/latest/channel.html#fromfilepairs) channel factory method. The second one, declares a variable named `GATK` specifying the path of the GATK application file.

    Once you think you have data organised, you can again run the pipeline. However this time, we can use the the `-resume` flag.

    ```bash
    nextflow run main.nf -resume
    ```

    !!! tip

        See [here](https://www.nextflow.io/docs/latest/getstarted.html?highlight=resume#modify-and-resume) for more details about using the `resume` option.

    ??? solution


        ```groovy linenums="1" hl_lines="2"
        workflow {
            reads_ch = Channel.fromFilePairs(params.reads) // (1)!
        }
        ```

        1. Create a channel using [fromFilePairs()](https://www.nextflow.io/docs/latest/channel.html#fromfilepairs).

## Process 1A: Create a FASTA genome index

Now we have our inputs set up we can move onto the processes. In our first process we will create a genome index using [samtools](http://www.htslib.org/).

You should implement a process having the following structure:

-   **Name**: `prepare_genome_samtools`
-   **Command**: create a genome index for the genome fasta with samtools
-   **Input**: the genome fasta file
-   **Output**: the samtools genome index file

!!! exercise "Problem #2"

    Copy the code below and paste it at the end of `main.nf`.

    Your aim is to replace the `BLANK` placeholder with the the correct process call.

    ```groovy linenums="1" hl_lines="20"
    /*
     * Process 1A: Create a FASTA genome index with samtools
     */

    process prepare_genome_samtools {
        input:
        path genome

        output:
        path "${genome}.fai"

        script:
        """
        samtools faidx ${genome}
        """
    }

    workflow {
        reads_ch = Channel.fromFilePairs(params.reads)
        BLANK
    }
    ```

    In plain english, the process could be written as:

    -   A **process** called `prepare_genome_samtools`
    -   takes as **input** the genome file
    -   and creates as **output** a genome index file
    -   **script**: using samtools create the genome index from the genome file

    Now when we run the pipeline, we see that the process 1A is submitted:

    ```bash
    nextflow run main.nf -resume
    ```
    ```console
    N E X T F L O W  ~  version 23.04.1
    Launching `main.nf` [cranky_bose] - revision: d1df5b7267
    executor >  local (1)
    [cd/47f882] process > 1A_prepare_genome_samtools [100%] 1 of 1 ✔
    ```

    ??? solution

        ```groovy linenums="1" hl_lines="20"
        /*
         * Process 1A: Create a FASTA genome index with samtools
         */

        process prepare_genome_samtools {
            input:
            path genome

            output:
            path "${genome}.fai"

            script:
            """
            samtools faidx ${genome}
            """
        }

        workflow {
            reads_ch = Channel.fromFilePairs(params.reads)
            prepare_genome_samtools(params.genome) // (1)!
        }
        ```

        1. The solution is to provide **`params.genome`** as input to the `prepare_genome_samtools` process.

        !!! warning

            `params.genome` is just a regular variable, not a channel, but when passed as input to a process, it's automatically converted into a value channel.

## Process 1B: Create a FASTA genome sequence dictionary with Picard for GATK

Our first process created the genome index for GATK using samtools. For the next process we must do something very similar, this time creating a genome sequence dictionary using [Picard](https://broadinstitute.github.io/picard/).

You should implement a process having the following structure:

-   **Name**: `1B_prepare_genome_picard`
-   **Command**: create a genome dictionary for the genome fasta with Picard tools
-   **Input**: the genome fasta file
-   **Output**: the genome dictionary file

!!! exercise "Problem #3"

    Fill in the `BLANK` words for both the input and output sections.

    Copy the code below and paste it at the end of `main.nf`.

    Your aim is to insert the correct input name from into the input step (written as `BLANK`) of the process and run the pipeline.

    !!! info

        You can choose any channel output name that makes sense to you.

    ```groovy linenums="1" hl_lines="8 11"
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

    !!! note

        `.baseName` returns the filename without the file suffix. If `"${genome}"` is `human.fa`, then `"${genome.baseName}.dict"` would be `human.dict`.

    ??? solution

        ```groovy linenums="1" hl_lines="8 11"
        /*
         * Process 1B: Create a FASTA genome sequence dictionary with Picard for GATK
         */

        process '1B_prepare_genome_picard' {

            input:
            path genome from params.genome // (1)!

            output:
            path "${genome.baseName}.dict" into genome_dict_ch // (2)!

            script:
            """
            PICARD=`which picard.jar`
            java -jar \$PICARD CreateSequenceDictionary R= $genome O= ${genome.baseName}.dict
            """
        }
        ```

        1. Take as input the `genome` file from the `params.genome` parameter
        2. Give as output the file `${genome.baseName}.dict` and adds it to the channel `genome_dict_ch`

## Process 1C: Create STAR genome index file

Next we must create a genome index for the [STAR](https://github.com/alexdobin/STAR) mapping software.

You should implement a process having the following structure:

-   **Name**: 1C_prepare_star_genome_index
-   **Command**: create a STAR genome index for the genome fasta
-   **Input**: the genome fasta file
-   **Output**: a directory containing the STAR genome index

!!! exercise "Problem #4"

    This is a similar exercise as problem 3, except this time both `input` and `output` lines have been left `BLANK` and must be completed.

    ```groovy linenums="1" hl_lines="8 11"
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

    !!! info

        The output of the STAR genomeGenerate command is specified here as `genome_dir`.

    ??? solution

        ```groovy linenums="1" hl_lines="8 11"
        /*
        * Process 1C: Create the genome index file for STAR
        */

        process '1C_prepare_star_genome_index' {

            input:
            path genome from params.genome // (1)!

            output:
            path 'genome_dir' into genome_dir_ch // (2)!

            script: // (3)!
            """
            mkdir genome_dir

            STAR --runMode genomeGenerate \
            --genomeDir genome_dir \
            --genomeFastaFiles ${genome} \
            --runThreadN ${task.cpus}
            """
        }
        ```

        1. Take as input the `genome` file from the `params.genome` parameter.
        2. The `output` is a `file`\* called `genome_dir` and is added `into` a channel called `genome_dir_ch`. You can call the channel whatever you wish.
        3. Creates the output directory that will contain the resulting STAR genome index.

        !!! note

            The file in this case is a directory however it makes no difference.

## Process 1D: Filtered and recoded set of variants

Next on to something a little more tricky. The next process takes two inputs: the variants file and the blacklist file.

It should output a channel named `prepared_vcf_ch` which emitting a tuple of two files.

!!! info

    In Nextflow, tuples can be defined in the input or output using the [`tuple`](https://www.nextflow.io/docs/latest/process.html#input-of-type-tuple) qualifier.

You should implement a process having the following structure:

-   **Name**
    -   `1D_prepare_vcf_file`
-   **Command**
    -   create a filtered and recoded set of variants
-   **Input**
    -   the variants file
    -   the blacklisted regions file
-   **Output**
    -   a tuple containing the filtered/recoded VCF file and the tab index (TBI) file.

!!! exercise "Problem #5"

    You must fill in the two `BLANK_LINES` in the input and the two `BLANK` output files.

    ```groovy linenums="1" hl_lines="8-9 12"
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
        vcftools --gzvcf $variantsFile -c \ #
                 --exclude-bed ${blacklisted} \
                 --recode | bgzip -c \
                 > ${variantsFile.baseName}.filtered.recode.vcf.gz

        tabix ${variantsFile.baseName}.filtered.recode.vcf.gz
        """
    }
    ```


    Broken down, here is what the script is doing:

    ```bash
    vcftools --gzvcf $variantsFile -c \ # (1)!
             --exclude-bed ${blacklisted} \ # (2)!
             --recode | bgzip -c \
             > ${variantsFile.baseName}.filtered.recode.vcf.gz # (3)!

    tabix ${variantsFile.baseName}.filtered.recode.vcf.gz # (4)!
    ```

    1.   The input variable for the variants file
    2.   The input variable for the blacklist file
    3.   The first of the two output files
    4.   Generates the second output file named `"${variantsFile.baseName}.filtered.recode.vcf.gz.tbi"`

    Try run the pipeline from the project directory with:

    ```bash
    nextflow run main.nf -resume
    ```

    ??? solution

        ```groovy linenums="1" hl_lines="8-9 12-13"
        /*
         * Process 1D: Create a file containing the filtered and recoded set of variants
         */

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
        ```

        - Take as input the variants file, assigning the name `${variantsFile}`.
        - Take as input the blacklisted file, assigning the name `${blacklisted}`.
        - Out a tuple (or set) of two files into the `prepared_vcf_ch` channel.
        - Defines the name of the first output file.
        - Generates the second output file (with `.tbi` suffix).

Congratulations! Part 1 is now complete.

We have all the data prepared and into channels ready for the more serious steps

## Process 2: STAR Mapping

In this process, for each sample, we align the reads to our genome using the STAR index we created previously.

You should implement a process having the following structure:

-   **Name**
    -   `2_rnaseq_mapping_star`
-   **Command**
    -   mapping of the RNA-Seq reads using STAR
-   **Input**
    -   the genome fasta file
    -   the STAR genome index
    -   a tuple containing the replicate id and paired read files
-   **Output**
    -   a tuple containing replicate id, aligned bam file & aligned bam file index

!!! Exercise "Problem #6"

    Copy the code below and paste it at the end of `main.nf`.

    You must fill in the three `BLANK_LINE` lines in the input and the one `BLANK_LINE` line in the output.

    ```groovy linenums="1" hl_lines="8-10 13"
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

        # 2nd pass (improve alignments using table of splice junctions and create a new index)
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

    !!! info

        The final command produces an bam index which is the full filename with an additional `.bai` suffix.

    ??? solution

        ```groovy linenums="1" hl_lines="8-10 13"
        /*
         * Process 2: Align RNA-Seq reads to the genome with STAR
         */

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

            # 2nd pass (improve alignments using table of splice junctions and create a new index)
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

The next step is a filtering step using GATK. For each sample, we split all the reads that contain N characters in their [CIGAR](http://genome.sph.umich.edu/wiki/SAM#What_is_a_CIGAR.3F) string.

## Process 3: GATK Split on N

The process creates `k+1` new reads (where `k` is the number of `N` cigar elements) that correspond to the segments of the original read beside/between the splicing events represented by the `N`s in the original CIGAR.

You should implement a process having the following structure:

-   **Name**
    -   3_rnaseq_gatk_splitNcigar
-   **Command**
    -   split reads on Ns in CIGAR string using GATK
-   **Input**
    -   the genome fasta file
    -   the genome index made with samtools
    -   the genome dictionary made with picard
    -   a tuple containing replicate id, aligned bam file and aligned bam file index from the STAR mapping
-   **Output**
    -   a tuple containing the replicate id, the split bam file and the split bam index file

!!! exercise "Problem #7"

    Copy the code below and paste it at the end of `main.nf`.

    You must fill in the four `BLANK_LINE` lines in the input and the one `BLANK_LINE` line in the output.

    !!! warning

        There is an optional [`tag`](https://www.nextflow.io/docs/latest/process.html#tag) line added to the start of this process. The [`tag`](https://www.nextflow.io/docs/latest/process.html#tag) line allows you to assign a name to a specific task (single execution of a process). This is particularly useful when there are many samples/replicates which pass through the same process.


    ```groovy linenums="1" hl_lines="2 5-8 11"
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

    !!! info

        The GATK command above automatically creates a bam index (`.bai`) of the `split.bam` output file

    !!! example

        A `tag` line would also be useful in [Process 2](#process-2)

    ??? solution


        ```groovy linenums="1" hl_lines="2 5-8 11"
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

Next we perform a Base Quality Score Recalibration step using GATK.

## Process 4: GATK Recalibrate

This step uses GATK to detect systematic errors in the base quality scores, select unique alignments and then index the resulting bam file with samtools. You can find details of the specific GATK BaseRecalibrator parameters [here](https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_bqsr_BaseRecalibrator.php).

You should implement a process having the following structure:

-   **Name**
    -   4_rnaseq_gatk_recalibrate
-   **Command**
    -   recalibrate reads from each replicate using GATK
-   **Input**
    -   the genome fasta file
    -   the genome index made with samtools
    -   the genome dictionary made with picard
    -   a tuple containing replicate id, aligned bam file and aligned bam file index from process 3
    -   a tuple containing the filtered/recoded VCF file and the tab index (TBI) file from process 1D
-   **Output**
    -   a tuple containing the sample id, the unique bam file and the unique bam index file

!!! exercise "Problem #8"

    Copy the code below and paste it at the end of `main.nf`.

    You must fill in the five `BLANK_LINE` lines in the input and the one `BLANK` in the output line.

    ```groovy linenums="1" hl_lines="5-9 12"
    process '4_rnaseq_gatk_recalibrate' {
        tag "$replicateId"

        input:
        BLANK_LINE
        BLANK_LINE
        BLANK_LINE
        BLANK_LINE
        BLANK_LINE

        output:
        BLANK into (final_output_ch, bam_for_ASE_ch) // (1)!

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

    1. The files resulting from this process will be used in two downstream processes. If a process is executed more than once, and the downstream channel is used by more than one process, we must duplicate the channel. We can do this using the `into` operator with parenthesis in the output section. See [here](https://www.nextflow.io/docs/latest/operator.html#into) for more information on using `into`.

    -   The unique bam file
    -   The index of the unique bam file (bam file name + `.bai`)

    ??? solution


        ```groovy linenums="1" hl_lines="5-9 12"
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

Now we are ready to perform the variant calling with GATK.

## Process 5: GATK Variant Calling

This steps call variants with GATK HaplotypeCaller. You can find details of the specific GATK HaplotypeCaller parameters [here](https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_gatk_tools_walkers_haplotypecaller_HaplotypeCaller.php).

You should implement a process having the following structure:

-   **Name**
    -   5_rnaseq_call_variants
-   **Command**
    -   variant calling of each sample using GATK
-   **Input**
    -   the genome fasta file
    -   the genome index made with samtools
    -   the genome dictionary made with picard
    -   a tuple containing replicate id, aligned bam file and aligned bam file index from process 4
-   **Output**
    -   a tuple containing the sample id the resulting variant calling file (vcf)

!!! exercise "Problem #9"

    In this problem we will introduce the use of a channel operator in the input section. The [groupTuple](https://www.nextflow.io/docs/latest/operator.html#grouptuple) operator groups together the tuples emitted by a channel which share a common key.

    !!! warning

        Note that in process 4, we used the sampleID (not replicateID) as the first element of the tuple in the output. Now we combine the replicates by grouping them on the sample ID. It follows from this that process 4 is run one time per replicate and process 5 is run one time per sample.

    Fill in the `BLANK_LINE` lines and `BLANK` words as before.

    ```groovy linenums="1" hl_lines="2 5-8 11"
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

    ??? solution

        ```groovy linenums="1" hl_lines="2 5-8 11"
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

        -   [`tag`](https://www.nextflow.io/docs/latest/process.html#tag) line with the using the sample id as the tag.
        -   the genome fasta file.
        -   the genome index from the `genome_index_ch` channel created in the process `1A_prepare_genome_samtools`.
        -   the genome dictionary from the `genome_dict_ch` channel created in the process `1B_prepare_genome_picard`.
        -   the sets grouped by sampleID from the `final_output_ch` channel created in the process `4_rnaseq_gatk_recalibrate`.
        -   the set containing the sample ID and final VCF file.
        -   the line specifying the name resulting final vcf file.

## Processes 6A and 6B: ASE & RNA Editing

In the final steps we will create processes for Allele-Specific Expression and RNA Editing Analysis.

We must process the VCF result to prepare variants file for allele specific expression (ASE) analysis. We will implement both processes together.

You should implement two processes having the following structure:

-   _1st process_
    -   **Name**
        -   6A_post_process_vcf
    -   **Command**
        -   post-process the variant calling file (vcf) of each sample
    -   **Input**
        -   tuple containing the sample ID and vcf file
        -   a tuple containing the filtered/recoded VCF file and the tab index (TBI) file from process 1D
    -   **Output**
        -   a tuple containing the sample id, the variant calling file (vcf) and a file containing common SNPs
-   _2nd process_
    -   **Name**
        -   6B_prepare_vcf_for_ase
    -   **Command**
        -   prepare the VCF for allele specific expression (ASE) and generate a figure in R.
    -   **Input**
        -   a tuple containing the sample id, the variant calling file (vcf) and a file containing common SNPs
    -   **Output**
        -   a tuple containing the sample ID and known SNPs in the sample for ASE
        -   a figure of the SNPs generated in R as a PDF file

!!! exercise "Problem #10"

    Here we introduce the `publishDir` directive. This allows us to specify a location for the outputs of the process. See [here](https://www.nextflow.io/docs/latest/process.html#publishdir) for more details.

    You must have the output of process 6A become the input of process 6B.

    ```groovy linenums="1" hl_lines="2 6-7 10 21-22 25 28-29"
    process '6A_post_process_vcf' {
        tag BLANK
        publishDir "$params.results/$sampleId" // (1)!

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

    1. here the output location is specified as a combination of a pipeline parameter and a process input variable

    ??? solution


        ```groovy linenums="1" hl_lines="2 6-7 10 21-22 25 28-29"
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
        ```

The final step is the GATK ASEReadCounter.

!!! exercise "Problem #11"

    We have seen the basics of using processes in Nextflow. Yet one of the features of Nextflow is the operations that can be performed on channels outside of processes. See [here](https://www.nextflow.io/docs/latest/operator.html) for details on the specific operators.

    Before we perform the GATK ASEReadCounter process, we must group the data for allele-specific expression. To do this we must combine channels.

    The `bam_for_ASE_ch` channel emites tuples having the following structure, holding the final BAM/BAI files:

    ```bash
    < sample_id, file_bam, file_bai >
    ```

    The `vcf_for_ASE` channel emits tuples having the following structure:

    ```bash
    < sample_id, output.vcf >
    ```

    In the first operation, the BAMs are grouped together by sample id.

    Next, this resulting channel is merged with the VCFs (vcf_for_ASE) having the same sample id.

    We must take the merged channel and creates a channel named `grouped_vcf_bam_bai_ch` emitting the following tuples:

    ```bash
    < sample_id, file_vcf, List[file_bam], List[file_bai] >
    ```

    Your aim is to fill in the `BLANKS` below.

    ```groovy linenums="1" hl_lines="2 9"
    bam_for_ASE_ch
        .BLANK // (1)!
        .phase(vcf_for_ASE) // (2)!
        .map{ left, right -> // (3)!
            def sampleId = left[0] // (4)!
            def bam = left[1] // (5)!
            def bai = left[2] // (6)!
            def vcf = right [1] // (7)!
            tuple(BLANK, vcf, BLANK, BLANK) // (8)!
        }
        .set { grouped_vcf_bam_bai_ch } // (9)!
    ```

    1.   an operator that groups tuples that contain a common first element.
    2.   the phase operator synchronizes the values emitted by two other channels. See [here](https://www.nextflow.io/docs/latest/operator.html?phase#phase) for more details
    3.   the map operator can apply any function to every item on a channel. In this case we take our tuple from the phase operation, define the separate elements and create a new tuple.
    4.   define `sampleId` to be the first element of left.
    5.   define bam to be the second element of left.
    6.   define bai to be the third element of left.
    7.   define vcf to be the first element of right.
    8.   create a new tuple made of four elements
    9.   rename the resulting as `grouped_vcf_bam_bai_ch`

    !!! note

        `left` and `right` above are arbitrary names. From the phase operator documentation, we see that phase returns pairs of items. So here `left` originates from contents of the `bam_for_ASE_ch` channel and `right` originates from the contents of `vcf_for_ASE` channel.


    ??? solution

        ```groovy linenums="1" hl_lines="2 9"
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
        ```

## Process 6C: Allele-Specific Expression analysis with GATK ASEReadCounter

Now we are ready for the final process.

You should implement a process having the following structure:

-   **Name**
    -   6C_ASE_knownSNPs
-   **Command**
    -   calculate allele counts at a set of positions with GATK tools
-   **Input**
    -   genome fasta file
    -   genome index file from samtools
    -   genome dictionary file
    -   the `grouped_vcf_bam_bai_ch`channel
-   **Output**
    -   the allele specific expression file (`ASE.tsv`)

!!! exercise "Problem #12"

    You should construct the process and run the pipeline in its entirety.

    ```bash linenums="1"
    echo "${bam.join('\n')}" > bam.list

    java -jar $GATK -R ${genome} \
                    -T ASEReadCounter \
                    -o ASE.tsv \
                    -I bam.list \
                    -sites ${vcf}
    ```

    ??? solution

        ```groovy linenums="1"
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

Congratulations! If you made it this far you now have all the basics to create your own Nextflow workflows.

## Results overview

For each processed sample the pipeline stores results into a folder named after the sample identifier. These folders are created in the directory specified as a parameter in `params.results`.

Result files for this workshop can be found in the folder `results` within the current folder. There you should see a directory called `ENCSR000COQ/` containing the following files:

### Variant calls

`final.vcf`

This file contains all somatic variants (SNVs) called from RNAseq data. You will see variants that pass all filters, with the `PASS` keyword in the <span class="red">7th</span> field of the vcf file (`filter status`), and also those that did not pass one or more filters.

`commonSNPs.diff.sites_in_files`

Tab-separated file with comparison between variants obtained from RNAseq and "known" variants from DNA.

The file is sorted by genomic position and contains 8 fields:

|     |           |                                                                                           |
| --- | --------- | ----------------------------------------------------------------------------------------- |
| 1   | `CHROM`   | chromosome name;                                                                          |
| 2   | `POS1`    | position of the SNV in file #1 (RNAseq data);                                             |
| 3   | `POS2`    | position of SNV in file #2 (DNA "known" variants);                                        |
| 4   | `IN_FILE` | flag whether SNV is present in the file #1 _1_, in the file #2 _2_, or in both files _B_; |
| 5   | `REF1`    | reference sequence in the file 1;                                                         |
| 6   | `REF2`    | reference sequence in the file 2;                                                         |
| 7   | `ALT1`    | alternative sequence in the file 1;                                                       |
| 8   | `ALT2`    | alternative sequence in the file 2                                                        |

`known_snps.vcf`

Variants that are common to RNAseq and "known" variants from DNA.

### Allele specific expression quantification

`ASE.tsv`

Tab-separated file with allele counts at common SNVs positions (only SNVs from the file `known_snps.vcf`)

The file is sorted by coordinates and contains 13 fields:

|     |                 |                                                                                                                              |
| --- | --------------- | ---------------------------------------------------------------------------------------------------------------------------- |
| 1   | `contig`        | contig, scaffold or chromosome name of the variant                                                                           |
| 2   | `position`      | position of the variant                                                                                                      |
| 3   | `variant ID`    | variant ID in the dbSNP                                                                                                      |
| 4   | `refAllele`     | reference allele sequence                                                                                                    |
| 5   | `altAllele`     | alternate allele sequence                                                                                                    |
| 6   | `refCount`      | number of reads that support the reference allele                                                                            |
| 7   | `altCount`      | number of reads that support the alternate allele                                                                            |
| 8   | `totalCount`    | total number of reads at the site that support both reference and alternate allele and any other alleles present at the site |
| 9   | `lowMAPQDepth`  | number of reads that have low mapping quality                                                                                |
| 10  | `lowBaseQDepth` | number of reads that have low base quality                                                                                   |
| 11  | `rawDepth`      | total number of reads at the site that support both reference and alternate allele and any other alleles present at the site |
| 12  | `otherBases`    | number of reads that support bases other than reference and alternate bases                                                  |
| 13  | `improperPairs` | number of reads that have malformed pairs                                                                                    |

### Allele frequency histogram

`AF.histogram.pdf`

This file contains a histogram plot of allele frequency for SNVs common to RNA-seq and "known" variants from DNA.

## Bonus step

Until now the pipeline has been executed using just a single sample (`ENCSR000COQ1`).

Now we can re-execute the pipeline specifying a large set of samples by using the command shown below:

```bash
nextflow run main.nf -resume --reads 'data/reads/ENCSR000C*_{1,2}.fastq.gz'
```

Or run the final version of the Nextflow pipeline that is already prepared for you:

```bash
nextflow run final_main.nf -resume
```

It will print an output similar to the one below:

```console
N E X T F L O W  ~  version 23.04.1
Launching `main.nf` [nice_kirch] DSL2 - revision: 45de8f41e4
executor >  local (43)
[8f/0010ac] process > prepare_genome_samtools                [100%] 1 of 1 ✔
[7e/9404fb] process > prepare_genome_picard                  [100%] 1 of 1 ✔
[c1/52f4c2] process > prepare_star_genome_index              [100%] 1 of 1 ✔
[d7/edc7d4] process > prepare_vcf_file                       [100%] 1 of 1 ✔
[c6/469237] process > rnaseq_mapping_star (5)                [100%] 6 of 6 ✔
[c4/c689ec] process > rnaseq_gatk_splitNcigar (ENCSR000CPO1) [100%] 6 of 6 ✔
[8f/3b28cd] process > rnaseq_gatk_recalibrate (ENCSR000CPO2) [100%] 6 of 6 ✔
[8e/3305f0] process > rnaseq_call_variants (ENCSR000CPO)     [100%] 6 of 6 ✔
[fd/02e430] process > post_process_vcf (ENCSR000CPO)         [100%] 6 of 6 ✔
[11/9af6f7] process > prepare_vcf_for_ase (ENCSR000CPO)      [100%] 6 of 6 ✔
[33/6f0d90] process > ASE_knownSNPs (ENCSR000CPO)            [100%] 3 of 3 ✔
```

You can notice that this time the pipeline spawns the execution of more tasks because three samples have been provided instead of one.

This shows the ability of Nextflow to implicitly handle multiple parallel task executions depending on the specified pipeline input dataset.

A fully functional version of this pipeline is available at the following GitHub repository: [CalliNGS-NF](https://github.com/CRG-CNAG/CalliNGS-NF).
