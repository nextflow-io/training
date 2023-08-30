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
```

1. The `/*`, `*` and `*/` specify comment lines which are ignored by Nextflow.
2. The `baseDir` variable represents the main script path location.
3. The `reads` parameter uses a glob pattern to specify the forward (`ENCSR000COQ1_1.fq.gz`) and reverse (`ENCSR000COQ1_2.fq.gz`) reads (paired-end) of a sample.
4. The `results` parameter is used to specify a directory called `results`.

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
Launching `main.nf` [elated_davinci] DSL2 - revision: 5187dd3166
```

!!! exercise "Problem #1"

    Great, now we need to define a [channel](https://www.nextflow.io/docs/latest/channel.html) variable to handle the read-pair files. To do that open the `main.nf` file and copy the lines below at the end of the file.

    !!! tip

        In `code` you can move to the end of the file using ++ctrl+w++ and then ++ctrl+v++.

    This time you must fill the `BLANK` space with a channel factory that will create a channel out of the `params.reads` information.

    ```groovy linenums="1" hl_lines="2"
    workflow {
        reads_ch = BLANK
    }
    ```

    !!! tip

        Use the [fromFilePairs](https://www.nextflow.io/docs/latest/channel.html#fromfilepairs) channel factory.

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

        1. Creates a channel using the [fromFilePairs()](https://www.nextflow.io/docs/latest/channel.html#fromfilepairs) channel factory.

## Process 1A: Create a FASTA genome index

Now we have our inputs set up we can move onto the processes. In our first process we will create a genome index using [samtools](http://www.htslib.org/).

The first process has the following structure:

-   **Name**: `prepare_genome_samtools`
-   **Command**: create a genome index for the genome fasta with samtools
-   **Input**: the genome fasta file
-   **Output**: the samtools genome index file

!!! exercise "Problem #2"

    Copy the code below and paste it at the end of `main.nf`, removing the previous workflow block. Be careful not to accidently have multiple workflow blocks.

    Your aim is to replace the `BLANK` placeholder with the the correct process call.

    ```groovy linenums="1" hl_lines="23"
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

    Now when we run the pipeline, we see that the process `prepare_genome_samtools` is submitted:

    ```bash
    nextflow run main.nf -resume
    ```
    ```console
    N E X T F L O W  ~  version 23.04.1
    Launching `main.nf` [cranky_bose] - revision: d1df5b7267
    executor >  local (1)
    [cd/47f882] process > prepare_genome_samtools [100%] 1 of 1 ✔
    ```

    ??? solution

        ```groovy linenums="1" hl_lines="23"
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

The next process should have the following structure:

-   **Name**: `prepare_genome_picard`
-   **Command**: create a genome dictionary for the genome fasta with Picard tools
-   **Input**: the genome fasta file
-   **Output**: the genome dictionary file

!!! exercise "Problem #3"

    Your aim is to replace the `BLANK` placeholder with the the correct process call.

    Copy the code below and paste it at the end of `main.nf`, removing the previous workflow block.

    !!! info

        You can choose any channel output name that makes sense to you.

    ```groovy linenums="1" hl_lines="24"
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
        picard CreateSequenceDictionary R= $genome O= ${genome.baseName}.dict
        """
    }

    workflow {
        reads_ch = Channel.fromFilePairs(params.reads)

        prepare_genome_samtools(params.genome)
        BLANK
    }
    ```

    !!! note

        `.baseName` returns the filename without the file suffix. If `"${genome}"` is `human.fa`, then `"${genome.baseName}.dict"` would be `human.dict`.

    ??? solution

        ```groovy linenums="1" hl_lines="24"
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
            picard CreateSequenceDictionary R= $genome O= ${genome.baseName}.dict
            """
        }

        workflow {
            reads_ch = Channel.fromFilePairs(params.reads)

            prepare_genome_samtools(params.genome)
            prepare_genome_picard(params.genome) // (1)!
        }
        ```

        1. The solution is to provide **`params.genome`** as input to the `prepare_genome_picard` process.

## Process 1C: Create STAR genome index file

Next we must create a genome index for the [STAR](https://github.com/alexdobin/STAR) mapping software.

The next process has the following structure:

-   **Name**: `prepare_star_genome_index`
-   **Command**: create a STAR genome index for the genome fasta
-   **Input**: the genome fasta file
-   **Output**: a directory containing the STAR genome index

!!! exercise "Problem #4"

    This is a similar exercise as problem 3.

    ```groovy linenums="1" hl_lines="30"
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
        mkdir genome_dir

        STAR --runMode genomeGenerate \
             --genomeDir genome_dir \
             --genomeFastaFiles ${genome} \
             --runThreadN ${task.cpus}
        """
    }

    workflow {
        reads_ch = Channel.fromFilePairs(params.reads)

        prepare_genome_samtools(params.genome)
        prepare_genome_picard(params.genome)
        BLANK
    }
    ```

    !!! info

        The output of the STAR genomeGenerate command is specified here as `genome_dir`.

    ??? solution

        ```groovy linenums="1" hl_lines="30"
        /*
        * Process 1C: Create the genome index file for STAR
        */

        process prepare_star_genome_index {
            container 'quay.io/biocontainers/star:2.7.10b--h6b7c446_1'

            input:
            path genome // (1)!

            output:
            path 'genome_dir' // (2)!

            script: // (3)!
            """
            mkdir genome_dir

            STAR --runMode genomeGenerate \
                 --genomeDir genome_dir \
                 --genomeFastaFiles ${genome} \
                 --runThreadN ${task.cpus}
            """
        }

        workflow {
            reads_ch = Channel.fromFilePairs(params.reads)

            prepare_genome_samtools(params.genome)
            prepare_genome_picard(params.genome)
            prepare_star_genome_index(params.genome)
        }
        ```

        1. Take as input the `genome` file from the `params.genome` parameter.
        2. The `output` is a `path` called `genome_dir`
        3. Creates the output directory that will contain the resulting STAR genome index.

        !!! note

            The path in this case is a directory however it makes no difference, it could be a text file, for example.

## Process 1D: Filtered and recoded set of variants

Next on to something a little more tricky. The next process takes two inputs: the variants file and the blacklist file.

!!! info

    In Nextflow, tuples can be defined in the input or output using the [`tuple`](https://www.nextflow.io/docs/latest/process.html#input-type-tuple) qualifier.

The next process has the following structure:

-   **Name**: `prepare_vcf_file`
-   **Command**: create a filtered and recoded set of variants
-   **Input**:
    -   the variants file
    -   the blacklisted regions file
-   **Output**: a tuple containing the filtered/recoded VCF file and the tab index (TBI) file.

!!! exercise "Problem #5"

    You must fill in the `BLANK`.

    ```groovy linenums="1" hl_lines="33"
    /*
     * Process 1D: Create a file containing the filtered and recoded set of variants
     */

    process prepare_vcf_file {
        container 'quay.io/biocontainers/mulled-v2-b9358559e3ae3b9d7d8dbf1f401ae1fcaf757de3:ac05763cf181a5070c2fdb9bb5461f8d08f7b93b-0'

        input:
        path variantsFile
        path blacklisted

        output:
        tuple path("${variantsFile.baseName}.filtered.recode.vcf.gz"), \
              path("${variantsFile.baseName}.filtered.recode.vcf.gz.tbi")

        script:
        """
        vcftools --gzvcf $variantsFile -c \ #
                 --exclude-bed ${blacklisted} \
                 --recode | bgzip -c \
                 > ${variantsFile.baseName}.filtered.recode.vcf.gz

        tabix ${variantsFile.baseName}.filtered.recode.vcf.gz
        """
    }

    workflow {
        reads_ch = Channel.fromFilePairs(params.reads)

        prepare_genome_samtools(params.genome)
        prepare_genome_picard(params.genome)
        prepare_star_genome_index(params.genome)
        BLANK
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

    1.   The `$variantsFile` variable contains the path to the file with the known variants
    2.   The `$blacklisted` variable contains the path to the file with the genomic locations which are known to produce artifacts and spurious variants
    3.   The `>` symbol is used to redirect the output to the file specified after it
    4.   `tabix` is used here to create the second output that we want to consider from this process

    Try run the pipeline from the project directory with:

    ```bash
    nextflow run main.nf -resume
    ```

    ??? solution

        ```groovy linenums="1" hl_lines="33"
        /*
         * Process 1D: Create a file containing the filtered and recoded set of variants
         */

        process prepare_vcf_file {
            container 'quay.io/biocontainers/mulled-v2-b9358559e3ae3b9d7d8dbf1f401ae1fcaf757de3:ac05763cf181a5070c2fdb9bb5461f8d08f7b93b-0'

            input:
            path variantsFile
            path blacklisted

            output:
            tuple path("${variantsFile.baseName}.filtered.recode.vcf.gz"), \
                  path("${variantsFile.baseName}.filtered.recode.vcf.gz.tbi")

            script:
            """
            vcftools --gzvcf $variantsFile -c \
                     --exclude-bed ${blacklisted} \
                     --recode | bgzip -c \
                     > ${variantsFile.baseName}.filtered.recode.vcf.gz

            tabix ${variantsFile.baseName}.filtered.recode.vcf.gz
            """
        }

        workflow {
            reads_ch = Channel.fromFilePairs(params.reads)

            prepare_genome_samtools(params.genome)
            prepare_genome_picard(params.genome)
            prepare_star_genome_index(params.genome)
            prepare_vcf_file(params.variants, params.blacklist)
        }
        ```

        - Take as input the variants file, assigning the name `${variantsFile}`.
        - Take as input the blacklisted file, assigning the name `${blacklisted}`.
        - Out a tuple of two files
        - Defines the name of the first output file.
        - Generates the second output file (with `.tbi` suffix).

Congratulations! Part 1 is now complete.

We have all the data prepared and into channels ready for the more serious steps.

## Process 2: STAR Mapping

In this process, for each sample, we align the reads to our genome using the STAR index we created previously.

The process has the following structure:

-   **Name**: `rnaseq_mapping_star`
-   **Command**: mapping of the RNA-Seq reads using STAR
-   **Input**:
    -   the genome fasta file
    -   the STAR genome index
    -   a tuple containing the replicate id and paired read files
-   **Output**: a tuple containing replicate id, aligned bam file & aligned bam file index

!!! Exercise "Problem #6"

    Copy the code below and paste it at the end of `main.nf`, removing the previous workflow block.

    You must fill the `BLANK` space with the correct process call and inputs.

    ```groovy linenums="1" hl_lines="62"
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
        tuple val(replicateId), path('Aligned.sortedByCoord.out.bam'), path('Aligned.sortedByCoord.out.bam.bai')

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

    workflow {
        reads_ch = Channel.fromFilePairs(params.reads)

        prepare_genome_samtools(params.genome)
        prepare_genome_picard(params.genome)
        prepare_star_genome_index(params.genome)
        prepare_vcf_file(params.variants, params.blacklist)

        BLANK
    }
    ```

    !!! info

        The final command produces an bam index which is the full filename with an additional `.bai` suffix.

    ??? solution

        ```groovy linenums="1" hl_lines="62-64"
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
            tuple val(replicateId), path('Aligned.sortedByCoord.out.bam'), path('Aligned.sortedByCoord.out.bam.bai')

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

        workflow {
            reads_ch = Channel.fromFilePairs(params.reads)

            prepare_genome_samtools(params.genome)
            prepare_genome_picard(params.genome)
            prepare_star_genome_index(params.genome)
            prepare_vcf_file(params.variants, params.blacklist)

            rnaseq_mapping_star(params.genome,
                                prepare_star_genome_index.out,
                                reads_ch)
        }
        ```

The next step is a filtering step using GATK. For each sample, we split all the reads that contain N characters in their [CIGAR](http://genome.sph.umich.edu/wiki/SAM#What_is_a_CIGAR.3F) string.

## Process 3: GATK Split on N

The process creates `k+1` new reads (where `k` is the number of `N` cigar elements) that correspond to the segments of the original read beside/between the splicing events represented by the `N`s in the original CIGAR.

The next process has the following structure:

-   **Name**: `rnaseq_gatk_splitNcigar`
-   **Command**: split reads on Ns in CIGAR string using GATK
-   **Input**:
    -   the genome fasta file
    -   the genome index made with samtools
    -   the genome dictionary made with picard
    -   a tuple containing replicate id, aligned bam file and aligned bam file index from the STAR mapping
-   **Output**: a tuple containing the replicate id, the split bam file and the split bam index file

!!! exercise "Problem #7"

    Copy the code below and paste it at the end of `main.nf`, removing the previous workflow block.

    You must fill the `BLANK` space with the correct process call and inputs.

    !!! warning

        There is an optional [`tag`](https://www.nextflow.io/docs/latest/process.html#tag) line added to the start of this process. The [`tag`](https://www.nextflow.io/docs/latest/process.html#tag) line allows you to assign a name to a specific task (single instance of a process). This is particularly useful when there are many samples/replicates which pass through the same process.


    ```groovy linenums="1" hl_lines="43"
    /*
     * Process 3: GATK Split on N
     */

    process rnaseq_gatk_splitNcigar {
        container 'quay.io/broadinstitute/gotc-prod-gatk:1.0.0-4.1.8.0-1626439571'
        tag "$replicateId"

        input:
        path genome
        path index
        path genome_dict
        tuple val(replicateId), path(bam), path(bai)

        output:
        tuple val(replicateId), path('split.bam'), path('split.bai')

        script:
        """
        # SplitNCigarReads and reassign mapping qualities
        java -jar /usr/gitc/GATK35.jar -T SplitNCigarReads \
                                       -R $genome -I $bam \
                                       -o split.bam \
                                       -rf ReassignOneMappingQuality \
                                       -RMQF 255 -RMQT 60 \
                                       -U ALLOW_N_CIGAR_READS \
                                       --fix_misencoded_quality_scores
        """
    }

    workflow {
        reads_ch = Channel.fromFilePairs(params.reads)

        prepare_genome_samtools(params.genome)
        prepare_genome_picard(params.genome)
        prepare_star_genome_index(params.genome)
        prepare_vcf_file(params.variants, params.blacklist)

        rnaseq_mapping_star(params.genome,
                            prepare_star_genome_index.out,
                            reads_ch)

        BLANK
    }
    ```

    !!! info

        The GATK command above automatically creates a bam index (`.bai`) of the `split.bam` output file

    !!! example

        A `tag` line would also be useful in [Process 2](#process-2-star-mapping)

    ??? solution


        ```groovy linenums="1" hl_lines="44-47"
        /*
         * Process 3: GATK Split on N
         */

        process rnaseq_gatk_splitNcigar {
            container 'quay.io/broadinstitute/gotc-prod-gatk:1.0.0-4.1.8.0-1626439571'
            tag "$replicateId"

            input:
            path genome
            path index
            path genome_dict
            tuple val(replicateId), path(bam), path(bai)

            output:
            tuple val(replicateId), path('split.bam'), path('split.bai')

            script:
            """
            # SplitNCigarReads and reassign mapping qualities
            java -jar /usr/gitc/GATK35.jar -T SplitNCigarReads \
                                           -R $genome -I $bam \
                                           -o split.bam \
                                           -rf ReassignOneMappingQuality \
                                           -RMQF 255 -RMQT 60 \
                                           -U ALLOW_N_CIGAR_READS \
                                           --fix_misencoded_quality_scores

            """
        }

        workflow {
            reads_ch = Channel.fromFilePairs(params.reads)

            prepare_genome_samtools(params.genome)
            prepare_genome_picard(params.genome)
            prepare_star_genome_index(params.genome)
            prepare_vcf_file(params.variants, params.blacklist)

            rnaseq_mapping_star(params.genome,
                                prepare_star_genome_index.out,
                                reads_ch)

            rnaseq_gatk_splitNcigar(params.genome,
                                    prepare_genome_samtools.out,
                                    prepare_genome_picard.out,
                                    rnaseq_mapping_star.out)
        }
        ```

        - [`tag`](https://www.nextflow.io/docs/latest/process.html#tag) line using the replicate id as the tag.
        - the genome fasta file
        - the genome index in the output channel from the `prepare_genome_samtools` process
        - the genome dictionary in the output channel from the `prepare_genome_picard` process
        - the set containing the aligned reads in the output channel from the `rnaseq_mapping_star` process
        - a set containing the sample id, the split bam file and the split bam index
        - specifies the input file names `$genome` and `$bam` to GATK
        - specifies the output file names to GATK

Next we perform a Base Quality Score Recalibration step using GATK.

## Process 4: GATK Recalibrate

This step uses GATK to detect systematic errors in the base quality scores, select unique alignments and then index the resulting bam file with samtools. You can find details of the specific GATK BaseRecalibrator parameters [here](https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_bqsr_BaseRecalibrator.php).

The next process has the following structure:

-   **Name**: `rnaseq_gatk_recalibrate`
-   **Command**: recalibrate reads from each replicate using GATK
-   **Input**
    -   the genome fasta file
    -   the genome index made with samtools
    -   the genome dictionary made with picard
    -   a tuple containing replicate id, aligned bam file and aligned bam file index from process 3
    -   a tuple containing the filtered/recoded VCF file and the tab index (TBI) file from process 1D
-   **Output**: a tuple containing the sample id, the unique bam file and the unique bam index file

!!! exercise "Problem #8"

    Copy the code below and paste it at the end of `main.nf`, removing the previous workflow block.

    Your aim is to replace the `BLANK` placeholder with the the correct process call.

    ```groovy linenums="1" hl_lines="67"
    /*
     * Process 4: GATK Recalibrate
     */

    process rnaseq_gatk_recalibrate {
        container 'quay.io/biocontainers/mulled-v2-aa1d7bddaee5eb6c4cbab18f8a072e3ea7ec3969:f963c36fd770e89d267eeaa27cad95c1c3dbe660-0'
        tag "$replicateId"

        input:
        path genome
        path index
        path dict
        tuple val(replicateId), path(bam), path(bai)
        tuple path(prepared_variants_file), path(prepared_variants_file_index)

        output:
        tuple val(sampleId), path("${replicateId}.final.uniq.bam"), path("${replicateId}.final.uniq.bam.bai")

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
        (samtools view -H final.bam; samtools view final.bam| grep -w 'NH:i:1') \
        |samtools view -Sb -  > ${replicateId}.final.uniq.bam

        # Index BAM files
        samtools index ${replicateId}.final.uniq.bam
        """
    }

    workflow {
        reads_ch = Channel.fromFilePairs(params.reads)

        prepare_genome_samtools(params.genome)
        prepare_genome_picard(params.genome)
        prepare_star_genome_index(params.genome)
        prepare_vcf_file(params.variants, params.blacklist)

        rnaseq_mapping_star(params.genome,
                            prepare_star_genome_index.out,
                            reads_ch)

        rnaseq_gatk_splitNcigar(params.genome,
                                prepare_genome_samtools.out,
                                prepare_genome_picard.out,
                                rnaseq_mapping_star.out)

        BLANK
    }
    ```

    -   The unique bam file
    -   The index of the unique bam file (bam file name + `.bai`)

    ??? solution


        ```groovy linenums="1" hl_lines="67-71"
        /*
         * Process 4: GATK Recalibrate
         */

        process rnaseq_gatk_recalibrate {
            container 'quay.io/biocontainers/mulled-v2-aa1d7bddaee5eb6c4cbab18f8a072e3ea7ec3969:f963c36fd770e89d267eeaa27cad95c1c3dbe660-0'
            tag "$replicateId"

            input:
            path genome
            path index
            path dict
            tuple val(replicateId), path(bam), path(bai)
            tuple path(prepared_variants_file), path(prepared_variants_file_index)

            output:
            tuple val(sampleId), path("${replicateId}.final.uniq.bam"), path("${replicateId}.final.uniq.bam.bai")

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
            (samtools view -H final.bam; samtools view final.bam| grep -w 'NH:i:1') \
            |samtools view -Sb -  > ${replicateId}.final.uniq.bam

            # Index BAM files
            samtools index ${replicateId}.final.uniq.bam
            """
        }

        workflow {
            reads_ch = Channel.fromFilePairs(params.reads)

            prepare_genome_samtools(params.genome)
            prepare_genome_picard(params.genome)
            prepare_star_genome_index(params.genome)
            prepare_vcf_file(params.variants, params.blacklist)

            rnaseq_mapping_star(params.genome,
                                prepare_star_genome_index.out,
                                reads_ch)

            rnaseq_gatk_splitNcigar(params.genome,
                                    prepare_genome_samtools.out,
                                    prepare_genome_picard.out,
                                    rnaseq_mapping_star.out)

            rnaseq_gatk_recalibrate(params.genome,
                                    prepare_genome_samtools.out,
                                    prepare_genome_picard.out,
                                    rnaseq_gatk_splitNcigar.out,
                                    prepare_vcf_file.out)
        }
        ```

        - the genome fasta file.
        - the genome index in the output channel from the `prepare_genome_samtools` process.
        - the genome dictionary in the output channel from the `prepare_genome_picard` process.
        - the set containing the split reads in the output channel from the `rnaseq_gatk_splitNcigar` process.
        - the set containing the filtered/recoded VCF file and the tab index (TBI) file in the output channel from the `prepare_vcf_file` process.
        - the set containing the replicate id, the unique bam file and the unique bam index file which goes into two channels.
        - line specifying the filename of the output bam file

Now we are ready to perform the variant calling with GATK.

## Process 5: GATK Variant Calling

This steps call variants with GATK HaplotypeCaller. You can find details of the specific GATK HaplotypeCaller parameters [here](https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_gatk_tools_walkers_haplotypecaller_HaplotypeCaller.php).

The next process has the following structure:

-   **Name**: `rnaseq_call_variants`
-   **Command**: variant calling of each sample using GATK
-   **Input**:
    -   the genome fasta file
    -   the genome index made with samtools
    -   the genome dictionary made with picard
    -   a tuple containing replicate id, aligned bam file and aligned bam file index from process 4
-   **Output**: a tuple containing the sample id the resulting variant calling file (vcf)

!!! exercise "Problem #9"

    Copy the code below and paste it at the end of `main.nf`, removing the previous workflow block. Be careful not to accidently have multiple workflow blocks.

    !!! warning

        Note that in process 4, we used the sampleID (not replicateID) as the first element of the tuple in the output. Now we combine the replicates by grouping them on the sample ID. It follows from this that process 4 is run one time per replicate and process 5 is run one time per sample.

    Your aim is to replace the `BLANK` placeholder with the the correct process call.

    ```groovy linenums="1" hl_lines="62"
    /*
     * Process 5: GATK Variant Calling
     */

    process rnaseq_call_variants {
        container 'quay.io/broadinstitute/gotc-prod-gatk:1.0.0-4.1.8.0-1626439571'
        tag "$sampleId"

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
                                       -R $genome -I bam.list \
                                       -dontUseSoftClippedBases \
                                       -stand_call_conf 20.0 \
                                       -o output.gatk.vcf.gz

        # Variant filtering
        java -jar /usr/gitc/GATK35.jar -T VariantFiltration \
                                       -R $genome -V output.gatk.vcf.gz \
                                       -window 35 -cluster 3 \
                                       -filterName FS -filter "FS > 30.0" \
                                       -filterName QD -filter "QD < 2.0" \
                                       -o final.vcf
        """
    }

    workflow {
        reads_ch = Channel.fromFilePairs(params.reads)

        prepare_genome_samtools(params.genome)
        prepare_genome_picard(params.genome)
        prepare_star_genome_index(params.genome)
        prepare_vcf_file(params.variants, params.blacklist)

        rnaseq_mapping_star(params.genome,
                            prepare_star_genome_index.out,
                            reads_ch)

        rnaseq_gatk_splitNcigar(params.genome,
                                prepare_genome_samtools.out,
                                prepare_genome_picard.out,
                                rnaseq_mapping_star.out)

        rnaseq_gatk_recalibrate(params.genome,
                                prepare_genome_samtools.out,
                                prepare_genome_picard.out,
                                rnaseq_gatk_splitNcigar.out,
                                prepare_vcf_file.out)

        BLANK
    }
    ```

    ??? solution

        ```groovy linenums="1" hl_lines="61-64"
        /*
         * Process 5: GATK Variant Calling
         */

        process rnaseq_call_variants {
            tag "$sampleId"

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
                                           -R $genome -I bam.list \
                                           -dontUseSoftClippedBases \
                                           -stand_call_conf 20.0 \
                                           -o output.gatk.vcf.gz

            # Variant filtering
            java -jar /usr/gitc/GATK35.jar -T VariantFiltration \
                                           -R $genome -V output.gatk.vcf.gz \
                                           -window 35 -cluster 3 \
                                           -filterName FS -filter "FS > 30.0" \
                                           -filterName QD -filter "QD < 2.0" \
                                           -o final.vcf
            """
        }

        workflow {
            reads_ch = Channel.fromFilePairs(params.reads)

            prepare_genome_samtools(params.genome)
            prepare_genome_picard(params.genome)
            prepare_star_genome_index(params.genome)
            prepare_vcf_file(params.variants, params.blacklist)

            rnaseq_mapping_star(params.genome,
                                prepare_star_genome_index.out,
                                reads_ch)

            rnaseq_gatk_splitNcigar(params.genome,
                                    prepare_genome_samtools.out,
                                    prepare_genome_picard.out,
                                    rnaseq_mapping_star.out)

            rnaseq_gatk_recalibrate(params.genome,
                                    prepare_genome_samtools.out,
                                    prepare_genome_picard.out,
                                    rnaseq_gatk_splitNcigar.out,
                                    prepare_vcf_file.out)

            rnaseq_call_variants(params.genome,
                                 prepare_genome_samtools.out,
                                 prepare_genome_picard.out,
                                 rnaseq_gatk_recalibrate.out)
        }
        ```

        -   [`tag`](https://www.nextflow.io/docs/latest/process.html#tag) line with the using the sample id as the tag.
        -   the genome fasta file.
        -   the genome index in the output channel from the `prepare_genome_samtools` process.
        -   the genome dictionary in the output channel from the `prepare_genome_picard` process.
        -   the sets grouped by sampleID in the output channel from the `rnaseq_gatk_recalibrate` process.
        -   the set containing the sample ID and final VCF file.
        -   the line specifying the name resulting final VCF file.

## Processes 6A and 6B: ASE & RNA Editing

In the final steps we will create processes for Allele-Specific Expression and RNA Editing Analysis.

We must process the VCF result to prepare variants file for allele specific expression (ASE) analysis. We will implement both processes together.

You should implement two processes having the following structure:

-   _1st process_
    -   **Name**: `post_process_vcf`
    -   **Command**: post-process the variant calling file (vcf) of each sample
    -   **Input**:
        -   tuple containing the sample ID and vcf file
        -   a tuple containing the filtered/recoded VCF file and the tab index (TBI) file from process 1D
    -   **Output**: a tuple containing the sample id, the variant calling file (vcf) and a file containing common SNPs
-   _2nd process_
    -   **Name**: `prepare_vcf_for_ase`
    -   **Command**: prepare the VCF for allele specific expression (ASE) and generate a figure in R.
    -   **Input**: a tuple containing the sample id, the variant calling file (vcf) and a file containing common SNPs
    -   **Output**:
        -   a tuple containing the sample ID and known SNPs in the sample for ASE
        -   a figure of the SNPs generated in R as a PDF file

!!! exercise "Problem #10"

    Here we introduce the `publishDir` directive. This allows us to specify a location for the outputs of the process. See [here](https://www.nextflow.io/docs/latest/process.html#publishdir) for more details.

    You must have the output of process 6A become the input of process 6B.

    ```groovy linenums="1" hl_lines="80"
    /*
     * Processes 6: ASE & RNA Editing
     */

    process post_process_vcf {
        container 'quay.io/biocontainers/mulled-v2-b9358559e3ae3b9d7d8dbf1f401ae1fcaf757de3:ac05763cf181a5070c2fdb9bb5461f8d08f7b93b-0'
        tag "$sampleId"
        publishDir "$params.results/$sampleId" // (1)!

        input:
        tuple val(sampleId), path('final.vcf')
        tuple path('filtered.recode.vcf.gz'), path('filtered.recode.vcf.gz.tbi')

        output:
        tuple val(sampleId), path('final.vcf'), path('commonSNPs.diff.sites_in_files')

        script:
        '''
        grep -v '#' final.vcf | awk '$7~/PASS/' |perl -ne 'chomp($_); ($dp)=$_=~/DP\\=(\\d+)\\;/; if($dp>=8){print $_."\\n"};' > result.DP8.vcf

        vcftools --vcf result.DP8.vcf --gzdiff filtered.recode.vcf.gz  --diff-site --out commonSNPs
        '''
    }

    process prepare_vcf_for_ase {
        tag "$sampleId"
        publishDir "$params.results/$sampleId"

        input:
        tuple val(sampleId), path('final.vcf'), path('commonSNPs.diff.sites_in_files')

        output:
        tuple val(sampleId), path('known_snps.vcf'), emit: vcf_for_ASE
        path 'AF.histogram.pdf'                    , emit: gghist_pdfs

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

    workflow {
        reads_ch = Channel.fromFilePairs(params.reads)

        prepare_genome_samtools(params.genome)
        prepare_genome_picard(params.genome)
        prepare_star_genome_index(params.genome)
        prepare_vcf_file(params.variants, params.blacklist)

        rnaseq_mapping_star(params.genome,
                            prepare_star_genome_index.out,
                            reads_ch)

        rnaseq_gatk_splitNcigar(params.genome,
                                prepare_genome_samtools.out,
                                prepare_genome_picard.out,
                                rnaseq_mapping_star.out)

        rnaseq_gatk_recalibrate(params.genome,
                                prepare_genome_samtools.out,
                                prepare_genome_picard.out,
                                rnaseq_gatk_splitNcigar.out,
                                prepare_vcf_file.out)

        rnaseq_call_variants(params.genome,
                             prepare_genome_samtools.out,
                             prepare_genome_picard.out,
                             rnaseq_gatk_recalibrate.out)

        BLANK
    }
    ```

    1. The path to the `publishDir` process directive consists of variables that will be evaluated before saving the files over there

    ??? solution


        ```groovy linenums="1" hl_lines="81-83"
        /*
         * Processes 6: ASE & RNA Editing
         */

        process post_process_vcf {
            container 'quay.io/biocontainers/mulled-v2-b9358559e3ae3b9d7d8dbf1f401ae1fcaf757de3:ac05763cf181a5070c2fdb9bb5461f8d08f7b93b-0'
            tag "$sampleId"
            publishDir "$params.results/$sampleId"

            input:
            tuple val(sampleId), path('final.vcf')
            tuple path('filtered.recode.vcf.gz'), path('filtered.recode.vcf.gz.tbi')

            output:
            tuple val(sampleId), path('final.vcf'), path('commonSNPs.diff.sites_in_files')

            script:
            '''
            grep -v '#' final.vcf | awk '$7~/PASS/' |perl -ne 'chomp($_); ($dp)=$_=~/DP\\=(\\d+)\\;/; if($dp>=8){print $_."\\n"};' > result.DP8.vcf

            vcftools --vcf result.DP8.vcf --gzdiff filtered.recode.vcf.gz  --diff-site --out commonSNPs
            '''
        }

        process prepare_vcf_for_ase {
            container 'cbcrg/callings-with-gatk:latest'
            tag "$sampleId"
            publishDir "$params.results/$sampleId"

            input:
            tuple val(sampleId), path('final.vcf'), path('commonSNPs.diff.sites_in_files')

            output:
            tuple val(sampleId), path('known_snps.vcf'), emit: vcf_for_ASE
            path 'AF.histogram.pdf'                    , emit: gghist_pdfs

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

        workflow {
            reads_ch = Channel.fromFilePairs(params.reads)

            prepare_genome_samtools(params.genome)
            prepare_genome_picard(params.genome)
            prepare_star_genome_index(params.genome)
            prepare_vcf_file(params.variants, params.blacklist)

            rnaseq_mapping_star(params.genome,
                                prepare_star_genome_index.out,
                                reads_ch)

            rnaseq_gatk_splitNcigar(params.genome,
                                    prepare_genome_samtools.out,
                                    prepare_genome_picard.out,
                                    rnaseq_mapping_star.out)

            rnaseq_gatk_recalibrate(params.genome,
                                    prepare_genome_samtools.out,
                                    prepare_genome_picard.out,
                                    rnaseq_gatk_splitNcigar.out,
                                    prepare_vcf_file.out)

            rnaseq_call_variants(params.genome,
                                 prepare_genome_samtools.out,
                                 prepare_genome_picard.out,
                                 rnaseq_gatk_recalibrate.out)

            post_process_vcf(rnaseq_call_variants.out,
                             prepare_vcf_file.out)
            prepare_vcf_for_ase(post_process_vcf.out)
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

    Next, this resulting channel is merged with the VCFs having the same sample id.

    We must take the merged channel and creates a channel named `grouped_vcf_bam_bai_ch` emitting the following tuples:

    ```bash
    < sample_id, file_vcf, List[file_bam], List[file_bai] >
    ```

    Your aim is to fill in the `BLANKS` below.

    ```groovy linenums="1" hl_lines="3-6"
    rnaseq_gatk_recalibrate
        .out
        .BLANK // (1)!
        .BLANK // (2)!
        .map { BLANK } // (3)!
        .set { BLANK } // (4)!
    ```

    1.   an operator that groups tuples that contain a common first element.
    2.   an operator that joins two channels taking a key into consideration. See [here](https://www.nextflow.io/docs/latest/operator.html?join#join) for more details
    3.   the map operator can apply any function to every item on a channel. In this case we take our tuple from the previous setp, define the separate elements and create a new tuple.
    4.   rename the resulting as `grouped_vcf_bam_bai_ch`

    ??? solution

        ```groovy linenums="1" hl_lines="35-38"
        workflow {
            reads_ch = Channel.fromFilePairs(params.reads)

            prepare_genome_samtools(params.genome)
            prepare_genome_picard(params.genome)
            prepare_star_genome_index(params.genome)
            prepare_vcf_file(params.variants, params.blacklist)

            rnaseq_mapping_star(params.genome,
                                prepare_star_genome_index.out,
                                reads_ch)

            rnaseq_gatk_splitNcigar(params.genome,
                                    prepare_genome_samtools.out,
                                    prepare_genome_picard.out,
                                    rnaseq_mapping_star.out)

            rnaseq_gatk_recalibrate(params.genome,
                                    prepare_genome_samtools.out,
                                    prepare_genome_picard.out,
                                    rnaseq_gatk_splitNcigar.out,
                                    prepare_vcf_file.out)

            rnaseq_call_variants(params.genome,
                                 prepare_genome_samtools.out,
                                 prepare_genome_picard.out,
                                 rnaseq_gatk_recalibrate.out)

            post_process_vcf(rnaseq_call_variants.out,
                             prepare_vcf_file.out)
            prepare_vcf_for_ase(post_process_vcf.out)

            rnaseq_gatk_recalibrate
                .out
                .groupTuple()
                .join(prepare_vcf_for_ase.out.vcf_for_ASE)
                .map { meta, bams, bais, vcf -> [meta, vcf, bams, bais] }
                .set { grouped_vcf_bam_bai_ch }
        }
        ```

## Process 7: Allele-Specific Expression analysis with GATK ASEReadCounter

Now we are ready for the final process.

The next process has the following structure:

-   **Name**: `ASE_knownSNPs`
-   **Command**: calculate allele counts at a set of positions with GATK tools
-   **Input**:
    -   genome fasta file
    -   genome index file from samtools
    -   genome dictionary file
    -   the `grouped_vcf_bam_bai_ch` channel
-   **Output**: the allele specific expression file (`ASE.tsv`)

!!! exercise "Problem #12"

    You should construct the process from scratch, add the process call and inputs to the workflow block and run the pipeline in its entirety.

    ```bash linenums="1"
    echo "${bam.join('\n')}" > bam.list

    java -jar /usr/gitc/GATK35.jar -R ${genome} \
                                   -T ASEReadCounter \
                                   -o ASE.tsv \
                                   -I bam.list \
                                   -sites ${vcf}
    ```

    ??? solution

        ```groovy linenums="1" hl_lines="5-29 70-73"
        /*
         * Processes 7: Allele-Specific Expression analysis with GATK ASEReadCounter
         */

        process ASE_knownSNPs {
            container 'quay.io/broadinstitute/gotc-prod-gatk:1.0.0-4.1.8.0-1626439571'
            tag "$sampleId"
            publishDir "$params.results/$sampleId"

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
            reads_ch = Channel.fromFilePairs(params.reads)

            prepare_genome_samtools(params.genome)
            prepare_genome_picard(params.genome)
            prepare_star_genome_index(params.genome)
            prepare_vcf_file(params.variants, params.blacklist)

            rnaseq_mapping_star(params.genome,
                                prepare_star_genome_index.out,
                                reads_ch)

            rnaseq_gatk_splitNcigar(params.genome,
                                prepare_genome_samtools.out,
                                prepare_genome_picard.out,
                                rnaseq_mapping_star.out)

            rnaseq_gatk_recalibrate(params.genome,
                                    prepare_genome_samtools.out,
                                    prepare_genome_picard.out,
                                    rnaseq_gatk_splitNcigar.out,
                                    prepare_vcf_file.out)

            rnaseq_call_variants(params.genome,
                                 prepare_genome_samtools.out,
                                 prepare_genome_picard.out,
                                 rnaseq_gatk_recalibrate.out)

            post_process_vcf(rnaseq_call_variants.out,
                             prepare_vcf_file.out)
            prepare_vcf_for_ase(post_process_vcf.out)

            rnaseq_gatk_recalibrate
                .out
                .groupTuple()
                .join(prepare_vcf_for_ase.out.vcf_for_ASE)
                .map { meta, bams, bais, vcf -> [meta, vcf, bams, bais] }
                .set { grouped_vcf_bam_bai_ch }

            ASE_knownSNPs(params.genome,
                          prepare_genome_samtools.out,
                          prepare_genome_picard.out,
                          grouped_vcf_bam_bai_ch)
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
