# Part 3: Hello Channels

In Part 2, you built a pipeline that was completely linear and processed each sample's data independently of the others. However, in real pipelines, you may need to combine data from multiple samples, or combine different kinds of data. Here we show you how to use channels and channel operators to implement a pipeline with more interesting plumbing.

Specifically, we show you how to implement joint variant calling with GATK, building on the pipeline from Part 2.

!!! note

    Don't worry if you're not familiar with GATK or genomics in general. We'll summarize the necessary concepts as we go, and the workflow implementation principles we demonstrate here apply broadly to any use case that follows a similar pattern.

### Method overview

The GATK variant calling method we used in Part 2 simply generated variant calls per sample. That's fine if you only want to look at the variants from each sample in isolation, but that only yields limited information. It's often more interesting to look at variant calls across multiple samples, and to do so, GATK offers an alternative method called joint variant calling, which we demonstrate here.

Joint variant calling involves generating a special kind of variant output called GVCF (for Genomic VCF) for each sample, then combining the GVCF data from all the samples and finally, running a 'joint genotyping' analysis to produce the final cohort-level VCF containing variant calls calculated based on the information from all samples.

![Joint analysis](img/joint-calling.png)

What's special about a sample's GVCF is that it contains records summarizing sequence data statistics about all positions in the targeted area of the genome, not just the positions where the program found evidence of variation. This is critical for the joint genotyping calculation. _TODO: ADD LINK FOR MORE READING_

The GVCF is produced by GATK HaplotypeCaller, the same tool we used in Part 2, with an additional parameter (`-ERC GVCF`). Combining the GVCFs is done with GATK GenomicsDBImport, which combines the per-sample calls into a data store (analogous to a database), then the actual 'joint genotyping' analysis is done with GATK GenotypeGVCFs.

So to recap, we're going to develop a workflow that does the following:

_ADD joint genotyping FLOWCHART_
<figure class="excalidraw">
--8<-- "docs/hello_nextflow/img/haplotype-caller.excalidraw.svg"
</figure>

1. Generate an index file for each BAM input file using Samtools
2. Run the GATK HaplotypeCaller on each BAM input file to generate a GVCF of per-sample genomic variant calls
3. Collect all the GVCFs and run joint genotyping on them to produce a cohort-level VCF.

### Dataset

-   **A reference genome** consisting of a small region of the human chromosome 20 (from hg19/b37) and its accessory files (index and sequence dictionary).
-   **Three whole genome sequencing samples** corresponding to a family trio (mother, father and son), which have been subset to a small portion on chromosome 20 to keep the file sizes small. The sequencing data is in [BAM](https://samtools.github.io/hts-specs/SAMv1.pdf) (Binary Alignment Map) format, i.e. genome sequencing reads that have already been mapped to the reference genome.
-   **A list of genomic intervals**, i.e. coordinates on the genome where our samples have data suitable for calling variants, provided in BED format.

---

## 0. Warmup: Run Samtools and GATK directly

Just like previously, we want to try out the commands manually before we attempt to wrap them in a workflow.

!!! note

     Make sure you're in the correct working directory:
     `cd /workspace/gitpod/hello-nextflow`

### 0.1. Index a BAM input file with Samtools

This first step is the same as in Part 2: Hello-GATK so you can skip it if you've already done that in this session. 

#### 0.1.1. Pull the samtools container

```bash
docker pull community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464
```

#### 0.1.2. Spin up the container interactively

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464
```

#### 0.1.3. Run the indexing command

```bash
samtools index /data/bam/reads_mother.bam
```

#### 0.1.4. Check that the BAM index has been produced

```bash
ls /data/bam/
```

This should show:

```console title="Output"
reads_father.bam      reads_mother.bam      reads_mother.bam.bai  reads_son.bam
```

Where `reads_mother.bam.bai` has been created as an index to `reads_mother.bam`.

#### 0.1.5. Exit the container

```bash
exit
```

### 0.2. Call variants with GATK HaplotypeCaller in GVCF mode

#### 0.2.1. Pull the GATK container

```bash
docker pull community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867
```

#### 0.2.2. Spin up the container interactively

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867
```

#### 0.2.3. Run the variant calling command with the GVCF option

This is the same base command as in Part 2, except this time we add the `-ERC GVCF` option, which switches on the HaplotypeCaller's GVCF mode to produce genomic VCFs.

```bash
gatk HaplotypeCaller \
        -R /data/ref/ref.fasta \
        -I /data/bam/reads_mother.bam \
        -O reads_mother.g.vcf \
        -L /data/ref/intervals.bed \
        -ERC GVCF
```

#### 0.2.4. Check the contents of the output file

```bash
cat reads_mother.g.vcf
```

_POINT OUT WHAT'S DIFFERENT COMPARED TO THE VCF IN PART 2_

#### 0.2.5. Repeat the process on the other two samples

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

### 0.3. Run joint genotyping 

#### 0.3.1. Combine all the per-sample GVCFs

```bash
_INSERT COMMAND_
```

#### 0.3.2. Run the joint genotyper

```bash
_INSERT COMMAND_
```

#### 0.3.3. Exit the container

```bash
exit
```

### Takeaway

You know how to _DO THE THING_.

### What's next?

_DO THE NEXT THING_.

---

## 1. Modify the per-sample variant calling step to produce GVCFs

#### 1.1. STEP NAME

```groovy title="hello-channels.nf"
RELEVANT CODE
```

#### 1.2. NEXT STEP

```groovy title="hello-channels.nf"
RELEVANT CODE
```

#### 1.X. Run it to verify that you can generate GVCFs

```bash
nextflow run hello-channels.nf
```

Should produce something like:

```console title="Output"
COPY OUTPUT
```

_ADD OUTPUT CHECK_

### Takeaway

You know how to _DO THE THING_.

### What's next?

_DO THE NEXT THING_.

---

## 2. Add joint genotyping step _UPDATE BASED ON ADAM'S PROPOSAL_



One slight complication is that these tools require the use of individually specified VCF files, and the syntax of the GenomicsDBImport tool looks like this:

```bash title="hello-gatk.nf"
gatk GenomicsDBImport \
    -V sample1.vcf.gz \
    -V sample2.vcf.gz \
    -V sample3.vcf.gz \
    ...
```

So to perform joint genotyping, we will need to collect all VCF files together in a single process and construct a command line for GenomicsDBImport.

#### 5.1. Write a process called GATK_JOINTGENOTYPING that wraps GenomicsDBImport and GenotypeGVCFs

```groovy title="hello-gatk.nf"
/*
 * Consolidate GVCFs and apply joint genotyping analysis
 */
process GATK_JOINTGENOTYPING {

    container "community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867"

    publishDir 'results', mode: 'copy'

    input:
        path vcfs
        path idxs
        val cohort_name
        path ref_fasta
        path ref_index
        path ref_dict
        path interval_list

    output:
        path "${cohort_name}.joint.vcf"
        path "${cohort_name}.joint.vcf.idx"

    script:
    """
    gatk GenomicsDBImport \
        -V ${vcfs} \
        --genomicsdb-workspace-path ${cohort_name}_gdb \
        -L ${interval_list}

    gatk GenotypeGVCFs \
        -R ${ref_fasta} \
        -V gendb://${cohort_name}_gdb \
        -O ${cohort_name}.joint.vcf \
        -L ${interval_list}
    """
}
```

#### 5.2. Add default value for the cohort name parameter up top

```groovy title="hello-gatk.nf"
// Base name for final output file
params.cohort_name = "family_trio"
```

#### 5.3. Gather the outputs of GATK_HAPLOTYPECALLER across samples using `collect()`

We collect the VCFs and their indices separately in order to list only the VCFs in the command we're going to construct. Since we'll give all of those files together to the joint genotyping process, we don't have to worry about the order of files like we did earlier.

Add this after the call to GATK_HAPLOTYPECALLER:

```groovy title="hello-gatk.nf"
// Collect variant calling outputs across samples
all_vcfs = GATK_HAPLOTYPECALLER.out[0].collect()
all_tbis = GATK_HAPLOTYPECALLER.out[1].collect()
```

!!! tip

    You can view the contents of the channel after performing this collect operation using `.view()`

#### 5.4. Add a call to the workflow block to run GATK_JOINTGENOTYPING

```groovy title="hello-gatk.nf"
// Consolidate GVCFs and apply joint genotyping analysis
GATK_JOINTGENOTYPING(
    all_vcfs,
    all_tbis,
    params.cohort_name,
    ref_file,
    ref_index_file,
    ref_dict_file,
    intervals_file
)
```

#### 5.5. Run the workflow

```bash
nextflow run hello-gatk.nf -resume
```

Oh no! The pipeline produces an error. When we dig into the console output, we can see the command executed isn't correct:

```bash
Command executed:

  gatk GenomicsDBImport -V reads_mother.bam.g.vcf reads_father.bam.g.vcf reads_son.bam.g.vcf --genomicsdb-workspace-path family_trio_gdb -L intervals.bed

  gatk GenotypeGVCFs -R ref.fasta -V gendb://family_trio_gdb -O family_trio.joint.vcf -L intervals.bed
```

Can you spot the error? We gave `gatk GenomicsDBImport` multiple VCF files for a single `-V` argument, but the tool expects a separate `-V` argument for each VCF file.

#### 5.6. Construct a command line with a separate `-V` argument for each input VCF

We use some string manipulations to repeat the VCFs with the argument `-V`, and replace `-V ${vcfs}` with the resulting `${vcfs_line}`:

_Before:_

```groovy title="hello-gatk.nf"
    script:
    """
    gatk GenomicsDBImport \
        -V ${vcfs} \
        --genomicsdb-workspace-path ${cohort_name}_gdb \
        -L ${interval_list}
```

_After:_

```groovy title="hello-gatk.nf"
    script:
    def vcfs_line = vcfs.collect { "-V ${it}" }.join(' ')
    """
    gatk GenomicsDBImport \
        ${vcfs_line} \
        --genomicsdb-workspace-path ${cohort_name}_gdb \
        -L ${interval_list}
```

#### 5.7. Run the workflow to verify that it generates the final VCF output as expected

```bash
nextflow run hello-gatk.nf -resume
```

Now we see the additional process show up in the log output (showing the compact view):

```console title="Output"
N E X T F L O W  ~  version 24.02.0-edge
Launching `hello-gatk.nf` [nauseous_thompson] DSL2 - revision: b346a53aae
executor >  local (7)
[d1/43979a] process > SAMTOOLS_INDEX (2)       [100%] 3 of 3 ✔
[20/247592] process > GATK_HAPLOTYPECALLER (3) [100%] 3 of 3 ✔
[14/7145b6] process > GATK_JOINTGENOTYPING (1) [100%] 1 of 1 ✔
```

You can find the final output file, `family_trio.joint.vcf`, in the work directory for the last process. Click on it to open it and you'll see 40 lines of metadata header followed by just under 30 jointly genotyped variant records (meaning at least one of the family members has a variant genotype at each genomic position listed).

!!! tip

    Keep in mind the data files covered only a tiny portion of chromosome 20; the real size of a variant callset would be counted in millions of variants. That's why we use only tiny subsets of data for training purposes!

### Takeaway

You know how to make a joint variant calling workflow that outputs a cohort VCF.

### What's next?

Celebrate your success and take an extra long break! This was tough and you deserve it.

In future trainings, you'll learn more sophisticated methods for managing inputs and outputs, as well as modularizing your code, testing it, and configuring execution.

**Good luck!**
