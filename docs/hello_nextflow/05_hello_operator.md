# Part 4: Hello Operator

In Part 3, you built a pipeline that was completely linear and processed each sample's data independently of the others.
However, in real pipelines, you may need to combine data from multiple samples, or combine different kinds of data.
Here we show you how to use channels and channel operators to implement a pipeline with more interesting plumbing.

Specifically, we show you how to implement joint variant calling with GATK, building on the pipeline from Part 2.

!!! note

    Don't worry if you're not familiar with GATK or genomics in general. We'll summarize the necessary concepts as we go, and the workflow implementation principles we demonstrate here apply broadly to any use case that follows a similar pattern.

### Method overview

The GATK variant calling method we used in Part 3 simply generated variant calls per sample.
That's fine if you only want to look at the variants from each sample in isolation, but that yields limited information.
It's often more interesting to look at variant calls differ across multiple samples, and to do so, GATK offers an alternative method called joint variant calling, which we demonstrate here.

Joint variant calling involves generating a special kind of variant output called GVCF (for Genomic VCF) for each sample, then combining the GVCF data from all the samples and finally, running a 'joint genotyping' statistical analysis.

![Joint analysis](img/joint-calling.png)

What's special about a sample's GVCF is that it contains records summarizing sequence data statistics about all positions in the targeted area of the genome, not just the positions where the program found evidence of variation.
This is critical for the joint genotyping calculation ([further reading](https://gatk.broadinstitute.org/hc/en-us/articles/360035890431-The-logic-of-joint-calling-for-germline-short-variants)).

The GVCF is produced by GATK HaplotypeCaller, the same tool we used in Part 3, with an additional parameter (`-ERC GVCF`).
Combining the GVCFs is done with GATK GenomicsDBImport, which combines the per-sample calls into a data store (analogous to a database), then the actual 'joint genotyping' analysis is done with GATK GenotypeGVCFs.

So to recap, we're going to develop a workflow that does the following:

<figure class="excalidraw">
--8<-- "docs/hello_nextflow/img/hello-gatk-2.svg"
</figure>

1. Generate an index file for each BAM input file using Samtools
2. Run the GATK HaplotypeCaller on each BAM input file to generate a GVCF of per-sample genomic variant calls
3. Collect all the GVCFs and combine them into a GenomicsDB data store
4. Run joint genotyping on the combined GVCF data store to produce a cohort-level VCF

### Dataset

-   **A reference genome** consisting of a small region of the human chromosome 20 (from hg19/b37) and its accessory files (index and sequence dictionary).
-   **Three whole genome sequencing samples** corresponding to a family trio (mother, father and son), which have been subset to a small portion on chromosome 20 to keep the file sizes small.
    The sequencing data is in [BAM](https://samtools.github.io/hts-specs/SAMv1.pdf) (Binary Alignment Map) format, _i.e._ genome sequencing reads that have already been mapped to the reference genome.
-   **A list of genomic intervals**, _i.e._ coordinates on the genome where our samples have data suitable for calling variants, provided in BED format.

---

## 0. Warmup: Run Samtools and GATK directly

Just like previously, we want to try out the commands manually before we attempt to wrap them in a workflow.

!!! note

     Make sure you're in the correct working directory:
     `cd /workspace/gitpod/hello-nextflow`

### 0.1. Index a BAM input file with Samtools

This first step is the same as in Part 3: Hello-Science, so it should feel very familiar, but this time we need to do it for all three samples.

!!! note

    We've technically already generated index files for the three samples through our pipeline, so we could go fish those out of the results directory. However, it's cleaner to just redo this manually, and it'll only take a minute.

#### 0.1.1. Spin up the Samtools container interactively

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464
```

#### 0.1.2. Run the indexing command for the three samples

```bash
samtools index /data/bam/reads_mother.bam
samtools index /data/bam/reads_father.bam
samtools index /data/bam/reads_son.bam
```

Just like previously, this should produce the index files in the same directory as the corresponding BAM files.

```console title="Directory contents"
data/bam/
├── reads_father.bam
├── reads_father.bam.bai
├── reads_mother.bam
├── reads_mother.bam.bai
├── reads_son.bam
└── reads_son.bam.bai
```

Now that we have index files for all three samples, we can proceed to generating the GVCFs for each of them.

#### 0.1.3. Exit the Samtools container

```bash
exit
```

### 0.2. Call variants with GATK HaplotypeCaller in GVCF mode

This second step is very similar to what we did Part 3: Hello-Science, but we are now going to run GATK in 'GVCF mode'.

#### 0.2.1. Spin up the GATK container interactively

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867
```

#### 0.2.2. Run the variant calling command with the GVCF option

In order to produce a genomic VCF (GVCF), we add the `-ERC GVCF` option to the base command, which switches on the HaplotypeCaller's GVCF mode.

We also change the file extension for the output file from `.vcf` to `.g.vcf`.
This is technically not a requirement, but it is a strongly recommended convention.

```bash
gatk HaplotypeCaller \
        -R /data/ref/ref.fasta \
        -I /data/bam/reads_mother.bam \
        -O reads_mother.g.vcf \
        -L /data/ref/intervals.bed \
        -ERC GVCF
```

This creates the GVCF output file `reads_mother.g.vcf` in the current working directory in the container.

If you `cat` it to view the contents, you'll see it's much longer than the equivalent VCF we generated in Part 3. You can't even scroll up to the start of the file, and most of the lines look quite different from what we saw in the VCF in Part 3.

```console title="Output" linenums="1674"
20_10037292_10066351    14714   .       T       <NON_REF>       .       .       END=14718       GT:DP:GQ:MIN_DP:PL       0/0:37:99:37:0,99,1192
20_10037292_10066351    14719   .       T       <NON_REF>       .       .       END=14719       GT:DP:GQ:MIN_DP:PL       0/0:36:82:36:0,82,1087
20_10037292_10066351    14720   .       T       <NON_REF>       .       .       END=14737       GT:DP:GQ:MIN_DP:PL       0/0:42:99:37:0,100,1160
```

These represent non-variant regions where the variant caller found no evidence of variation, so it captured some statistics describing its level of confidence in the absence of variation. This makes it possible to distinguish between two very different case figures: (1) there is good quality data showing that the sample is homozygous-reference, and (2) there is not enough good data available to make a determination either way.

In a GVCF, there are typically lots of such non-variant lines, with a smaller number of variant records sprinkled among them. Try running `head -176` on the GVCF to load in just the first 176 lines of the file to find an actual variant call.

```console title="Output" linenums="174"
20_10037292_10066351    3479    .       T       <NON_REF>       .       .       END=3479        GT:DP:GQ:MIN_DP:PL       0/0:34:36:34:0,36,906
20_10037292_10066351    3480    .       C       CT,<NON_REF>    503.03  .       DP=23;ExcessHet=0.0000;MLEAC=2,0;MLEAF=1.00,0.00;RAW_MQandDP=82800,23    GT:AD:DP:GQ:PL:SB       1/1:0,18,0:18:54:517,54,0,517,54,517:0,0,7,11
20_10037292_10066351    3481    .       T       <NON_REF>       .       .       END=3481        GT:DP:GQ:MIN_DP:PL       0/0:21:51:21:0,51,765
```

The second line shows the first variant record in the file, which corresponds to the first variant in the VCF file we looked at it Part 3.

Just like the original VCF was, the output GVCF file is also accompanied by an index file, called `reads_mother.g.vcf.idx`.

#### 0.2.3. Repeat the process on the other two samples

In order to test the joint genotyping step, we need GVCFs for all three samples, so let's generate those manually now.

```bash
gatk HaplotypeCaller \
        -R /data/ref/ref.fasta \
        -I /data/bam/reads_father.bam \
        -O reads_father.g.vcf \
        -L /data/ref/intervals.bed \
        -ERC GVCF
```

```bash
gatk HaplotypeCaller \
        -R /data/ref/ref.fasta \
        -I /data/bam/reads_son.bam \
        -O reads_son.g.vcf \
        -L /data/ref/intervals.bed \
        -ERC GVCF
```

Once this completes, you should have three files ending in `.g.vcf` in your work directory (one per sample) and their respective index files ending in `.g.vcf.idx`.

### 0.3. Run joint genotyping

Now that we have all the GVCFs, we can finally try out the joint genotyping approach to generating variant calls for a cohort of samples.
As a reminder, it's a two-step method that consists of combining the data from all the GVCFs into a data store, then running the joint genotyping analysis proper to generate the final VCF of joint-called variants.

#### 0.3.1. Combine all the per-sample GVCFs

This first step uses another GATK tool, called GenomicsDBImport, to combine the data from all the GVCFs into a GenomicsDB data store.

```bash
gatk GenomicsDBImport \
    -V reads_mother.g.vcf \
    -V reads_father.g.vcf \
    -V reads_son.g.vcf \
    -L /data/ref/intervals.bed \
    --genomicsdb-workspace-path family_trio_gdb
```

The output of this step is effectively a directory containing a set of further nested directories holding the combined variant data in the form of multiple different files.
You can poke around it but you'll quickly see this data store format is not intended to be read directly by humans.

!!! note

    GATK includes tools that make it possible to inspect and extract variant call data from the data store as needed.

#### 0.3.2. Run the joint genotyping analysis proper

This second step uses yet another GATK tool, called GenotypeGVCFs, to recalculate variant statistics and individual genotypes in light of the data available across all samples in the cohort.

```bash
gatk GenotypeGVCFs \
    -R /data/ref/ref.fasta \
    -V gendb://family_trio_gdb \
    -O family_trio.vcf
```

This creates the cohort-level VCF output file `family_trio.vcf` in the current working directory in the container.
It's another reasonably small file so you can `cat` this file to view its contents, and scroll up to find the first few variant lines.

````console title="family_trio.vcf" linenums="40"
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  reads_father    reads_mother    reads_son
20_10037292_10066351    3480    .       C       CT      1625.89 .       AC=5;AF=0.833;AN=6;BaseQRankSum=0.220;DP=85;ExcessHet=0.0000;FS=2.476;MLEAC=5;MLEAF=0.833;MQ=60.00;MQRankSum=0.00;QD=21.68;ReadPosRankSum=-1.147e+00;SOR=0.487    GT:AD:DP:GQ:PL  0/1:15,16:31:99:367,0,375       1/1:0,18:18:54:517,54,0 1/1:0,26:26:78:756,78,0
20_10037292_10066351    3520    .       AT      A       1678.89 .       AC=5;AF=0.833;AN=6;BaseQRankSum=1.03;DP=80;ExcessHet=0.0000;FS=2.290;MLEAC=5;MLEAF=0.833;MQ=60.00;MQRankSum=0.00;QD=22.39;ReadPosRankSum=0.701;SOR=0.730 GT:AD:DP:GQ:PL   0/1:18,13:31:99:296,0,424       1/1:0,18:18:54:623,54,0 1/1:0,26:26:78:774,78,0
20_10037292_10066351    3529    .       T       A       154.29  .       AC=1;AF=0.167;AN=6;BaseQRankSum=-5.440e-01;DP=104;ExcessHet=0.0000;FS=1.871;MLEAC=1;MLEAF=0.167;MQ=60.00;MQRankSum=0.00;QD=7.71;ReadPosRankSum=-1.158e+00;SOR=1.034       GT:AD:DP:GQ:PL  0/0:44,0:44:99:0,112,1347       0/1:12,8:20:99:163,0,328        0/0:39,0:39:99:0,105,1194

This looks more like the original VCF we generated in Part 3, except this time we have genotype-level information for all three samples. The last three columns in the file are the genotype blocks for the samples, listed in alphabetical order.

If we look at the genotypes called for our test family trio for the very first variant, we see that the father is heterozygous-variant (`0/1`), and the mother and son are both homozygous-variant (`1/1`).

That is ultimately the information we're looking to extract from the dataset! So let's go wrap all this into a Nextflow workflow so we can do this at scale.

#### 0.3.3. Exit the GATK container

```bash
exit
````

### Takeaway

You know how to run the individual commands in the terminal to verify that they will produce the information you want.

### What's next?

Wrap these commands into an actual pipeline.

---

## 1. Modify the per-sample variant calling step to produce a GVCF

The good news is that we don't need to start all over, since we already wrote a workflow that does some of this work in Part 3.
However, that pipeline produces VCF files, whereas now we want GVCF files in order to do the joint genotyping.
So we need to start by switching on the GVCF variant calling mode and updating the output file extension.

!!! note

    For convenience, we are going to work with a fresh copy of the GATK workflow as it stands at the end of Part 3, but under a different name: `hello-operator.nf`.

### 1.1. Tell HaplotypeCaller to emit a GVCF and update the output file path

Let's open the `hello-operator.nf` file in the code editor.
It should look very familiar, but feel free to run it if you want to satisfy yourself that it runs as expected.

We're going to start by making two changes:

-   Add the `-ERC GVCF` parameter to the GATK HaplotypeCaller command;
-   Update the output file path to use the corresponding `.g.vcf` extension, as per GATK convention.

_Before:_

```groovy title="hello-operator.nf"
    """
    gatk HaplotypeCaller \
        -R ${ref_fasta} \
        -I ${input_bam} \
        -O ${input_bam}.vcf \
        -L ${interval_list}
    """
```

_After:_

```groovy title="hello-operator.nf"
    """
    gatk HaplotypeCaller \
        -R ${ref_fasta} \
        -I ${input_bam} \
        -O ${input_bam}.g.vcf \
        -L ${interval_list} \
        -ERC GVCF
    """
```

That is all it takes to switch HaplotypeCaller to generating GVCFs instead of VCFs.

### 1.2. Run the pipeline to verify that you can generate GVCFs

The Nextflow execution command is the same as before, save for the workflow filename itself.
Make sure to update that appropriately.

```bash
nextflow run hello-operator.nf
```

Should produce something like:

_TODO: COPY OUTPUT_

```console title="Output"

```

If you open the file and scroll through it, you can see that GATK HaplotypeCaller produced a GVCF file, which contains additional information compared to the VCF file.

### Takeaway

Okay, this one was minimal in terms of Nextflow learning...
But hey, now you know the difference between a regular VCF and a Genomic VCF!

### What's next?

Learn to collect the contents of a channel and pass them in as a single input.

---

## 2. Collect and combine the GVCF data across all samples

We now need to collect all the GVCF files together and bundle them together as a single input to the next process, to feed into GATK GenomicsDBImport.

###

### 2.1. Write a process called GATK_GENOMICSDB

```groovy title="hello-operator.nf"
/*
 * Combine GVCFs into GenomicsDB datastore
 */
process GATK_GENOMICSDB {

    container "community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867"

    publishDir 'results', mode: 'copy'

    input:
        path gvcfs
        path idxs
        val cohort_name

    output:
        path "${cohort_name}_gdb"

    """
    gatk GenomicsDBImport \
        -V ${gvcfs} \
        --genomicsdb-workspace-path ${cohort_name}_gdb
    """
}
```

### 2.2. Add default value for the cohort name parameter up top

We need to provide an arbitrary name for the cohort.
Later in the series you'll learn how to use the available metadata for this, but for now we just declare a CLI parameter using `params` and give it a default value for convenience.

```groovy title="hello-operator.nf"
// Base name for final output file
params.cohort_name = "family_trio"
```

### 2.3. Gather the outputs of GATK_HAPLOTYPECALLER across samples using `collect()`

We collect the VCFs and their index files (`.idx`) separately in order to list only the VCFs in the command we're going to construct.
Since we'll give all of those files together to the joint genotyping process, we don't have to worry about the order of files like we did in Part 2.

Add this to the `workflow` body, right after the call to GATK_HAPLOTYPECALLER:

```groovy title="hello-operator.nf"
// Collect variant calling outputs across samples
all_gvcfs = GATK_HAPLOTYPECALLER.out[0].collect()
all_idxs = GATK_HAPLOTYPECALLER.out[1].collect()
```

!!! tip

    You can view the contents of the channel after performing this collect operation using `.view()`

### 2.4. Add a call to the workflow block to run GATK_GENOMICSDB

```groovy title="hello-operator.nf"
// Combine GVCFs into a GenomicsDB datastore
GATK_GENOMICSDB(
    all_gvcfs,
    all_idxs,
    params.cohort_name
)
```

### 2.5. Run the workflow

```bash
nextflow run hello-operator.nf -resume
```

Oh no! The pipeline produces an error. When we dig into the console output, we can see the command executed isn't correct:

```bash
Command executed:

  gatk GenomicsDBImport -V reads_mother.bam.g.vcf reads_father.bam.g.vcf reads_son.bam.g.vcf --genomicsdb-workspace-path family_trio_gdb
```

Can you spot the error? We gave `gatk GenomicsDBImport` multiple VCF files for a single `-V` argument, but the tool expects a separate `-V` argument for each GVCF file.

As a reminder, this was the command we ran in the container:

```bash
gatk GenomicsDBImport \
    -V /data/bam/reads_mother.vcf \
    -V /data/bam/reads_father.vcf \
    -V /data/bam/reads_son.vcf \
    --genomicsdb-workspace-path family_trio_gdb
```

### 2.6. Construct a command line with a separate `-V` argument for each input GVCF

We add the reserved keyword `script:` to declare the start of the script section explicitly, and we use some string manipulations to construct the `${gvcfs_line}`:

_Before:_

```groovy title="hello-operator.nf"
    """
    gatk GenomicsDBImport \
        -V ${gvcfs} \
        --genomicsdb-workspace-path ${cohort_name}_gdb
    """
```

_After:_

```groovy title="hello-operator.nf"
    script:
    def gvcfs_line = gvcfs.collect { "-V ${it}" }.join(' ')
    """
    gatk GenomicsDBImport \
        ${gvcfs_line} \
        --genomicsdb-workspace-path ${cohort_name}_gdb
    """
```

In this new code, `gvcfs.collect { "-V ${it}" }.join(' ')` takes each GVCF path (referenced by `${it}`), adds the `-V ` prefix then joins the resulting strings into a single string formatted correctly for the GenomicsDBImport tool.

### 2.7. Run the workflow to verify that it generates the GenomicsDB output as expected

```bash
nextflow run hello-operator.nf -resume
```

Now we see the additional process show up in the log output (showing the compact view):

_TODO: Update console output_

_TODO: check the output_

### Takeaway

Now you know how to collect outputs from a channel and bundle them as a single input to another process.
You also know how to construct a command line to provide the file names to the tool with the appropriate syntax.

### What's next?

Add a second command to the same process.

---

## 3. Run the joint genotyping step as part of the same process

Now that we have the combined genomic variant calls, we can run the joint genotyping tool, which will produce the final output that we actually care about: the VCF of cohort-level variant calls.

### 3.1. Rename the process from GATK_GENOMICSDB to GATK_JOINTGENOTYPING

Since the process will be running more than one tool, we change its name to refer to the overall operation rather than a single tool name.

_Before:_

```groovy title="hello-operator.nf"
/*
 * Combine GVCFs into GenomicsDB datastore
 */
process GATK_GENOMICSDB {
```

_After:_

```groovy title="hello-operator.nf"
/*
 * Combine GVCFs into GenomicsDB datastore and run joint genotyping to produce cohort-level calls
 */
process GATK_JOINTGENOTYPING {
```

### 3.2. Add the joint genotyping command to the GATK_JOINTGENOTYPING process

Simply add the second command after the first one inside the script section.

_Before:_

```groovy title="hello-operator.nf"
    """
    gatk GenomicsDBImport \
        -V ${gvcfs} \
        --genomicsdb-workspace-path ${cohort_name}_gdb
    """
```

_After:_

```groovy title="hello-operator.nf"
    """
    gatk GenomicsDBImport \
        ${gvcfs_line} \
        --genomicsdb-workspace-path ${cohort_name}_gdb

    gatk GenotypeGVCFs \
        -R ${ref_fasta} \
        -V gendb://${cohort_name}_gdb \
        -O ${cohort_name}.joint.vcf
    """
```

### 3.3. Add the reference genome files to the GATK_JOINTGENOTYPING process input definitions

The second command requires the reference genome files, so we need to add those to the process inputs.

_Before:_

```groovy title="hello-operator.nf"
input:
    path gvcfs
    path idxs
    val cohort_name
```

_After:_

```groovy title="hello-operator.nf"
input:
    path gvcfs
    path idxs
    val cohort_name
    path ref_fasta
    path ref_index
    path ref_dict
```

### 3.4. Update the process output definition to emit the VCF of cohort-level variant calls

We don't really care to save the GenomicsDB datastore; the output we're actually interested in is the VCF produced by the joint genotyping command.

_Before:_

````groovy title="hello-operator.nf"
output:
    path "${cohort_name}_gdb"
_After:_

```groovy title="hello-operator.nf"
output:
    path "${cohort_name}.joint.vcf"
    path "${cohort_name}.joint.vcf.idx"
````

### 3.5. Update the process call from GATK_GENOMICSDB to GATK_JOINTGENOTYPING

Let's rename the process call from GATK_GENOMICSDB to GATK_JOINTGENOTYPING.

And while we're at it, let's add the reference genome files as inputs, since we need to provide them to the joint genotyping tool.

_Before:_

```groovy title="hello-operator.nf"
// Combine GVCFs into a GenomicsDB datastore
GATK_GENOMICSDB(
    all_gvcfs,
    all_idxs,
    params.cohort_name
)
```

_After:_

```groovy title="hello-operator.nf"
// Combine GVCFs into a GenomicsDB datastore and apply joint genotyping
GATK_JOINTGENOTYPING(
    all_gvcfs,
    all_idxs,
    params.cohort_name,
    ref_file,
    ref_index_file,
    ref_dict_file
)
```

### 3.6. Run the workflow

Finally, we can run the modified workflow!

```bash
nextflow run hello-operator.nf -resume
```

The output should look like this:

```console title="Output"
N E X T F L O W  ~  version 24.02.0-edge
Launching `hello-operator.nf` [nauseous_thompson] DSL2 - revision: b346a53aae
executor >  local (7)
[d1/43979a] process > SAMTOOLS_INDEX (2)       [100%] 3 of 3 ✔
[20/247592] process > GATK_HAPLOTYPECALLER (3) [100%] 3 of 3 ✔
[14/7145b6] process > GATK_JOINTGENOTYPING (1)      [100%] 1 of 1 ✔
```

You can find the final output file, `family_trio.joint.vcf`, in the work directory for the last process.
Click on it to open it and you'll see 40 lines of metadata header followed by just under 30 jointly genotyped variant records (meaning at least one of the family members has a variant genotype at each genomic position listed).

!!! tip

    Keep in mind the data files covered only a tiny portion of chromosome 20; the real size of a variant callset would be counted in millions of variants. That's why we use only tiny subsets of data for training purposes!

### Takeaway

You know how to modify how channel contents are grouped using the collect() operator, in order to generate a cohort-level VCF of variant calls.

### What's next?

Celebrate your success and take an extra long break! This was tough and you deserve it.

In the next training, you'll learn how to manage metadata and use a more sophisticated samplesheet for passing inputs.

**Good luck!**
