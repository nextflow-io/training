# Part 2: Hello GATK

The [GATK](https://gatk.broadinstitute.org/) (Genome Analysis Toolkit) is a widely used software package developed by the Broad Institute to analyze high-throughput sequencing data. We're going to use GATK and a related tool, [Samtools](https://www.htslib.org/), in a very basic pipeline that identifies genomic variants through a method called **variant calling**.

![GATK pipeline](img/gatk-pipeline.png)

!!! note

    Don't worry if you're not familiar with GATK or genomics in general. We'll summarize the necessary concepts as we go, and the workflow implementation principles we demonstrate here apply broadly to any command line tool that takes in some input files and produce some output files.

A full variant calling pipeline typically involves a lot of steps. For simplicity, we are only going to look at the core variant calling steps.

### Method overview

1. Generate an index file for each BAM input file using Samtools
2. Run the GATK HaplotypeCaller on each BAM input file to generate per-sample variant calls in GVCF (Genomic Variant Call Format)

<figure class="excalidraw">
--8<-- "docs/hello_nextflow/img/haplotype-caller.excalidraw.svg"
</figure>

### Dataset

-   **A reference genome** consisting of a small region of the human chromosome 20 (from hg19/b37) and its accessory files (index and sequence dictionary).
-   **Three whole genome sequencing samples** corresponding to a family trio (mother, father and son), which have been subset to a small portion on chromosome 20 to keep the file sizes small. The sequencing data is in [BAM](https://samtools.github.io/hts-specs/SAMv1.pdf) (Binary Alignment Map) format, i.e. genome sequencing reads that have already been mapped to the reference genome.
-   **A list of genomic intervals**, i.e. coordinates on the genome where our samples have data suitable for calling variants, provided in BED format.

---

## 0. Warmup: Run Samtools and GATK directly

Just like in the Hello World example, we want to try out the commands manually before we attempt to wrap them in a workflow. The difference here is that we're going to use Docker containers to obtain and run the tools.

!!! note

     Make sure you're in the correct working directory:
     `cd /workspace/gitpod/hello-nextflow`

### 0.1. Index a BAM input file with Samtools

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

### 0.2. Call variants with GATK HaplotypeCaller

#### 0.2.1. Pull the GATK container

```bash
docker pull community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867
```

#### 0.2.2. Spin up the container interactively

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867
```

#### 0.2.3. Run the variant calling command

```bash
gatk HaplotypeCaller \
        -R /data/ref/ref.fasta \
        -I /data/bam/reads_mother.bam \
        -O reads_mother.vcf \
        -L /data/ref/intervals.bed
```

#### 0.2.4. Check the contents of the output file

The output file `reads_mother.vcf` is a small test file, so you can `cat` it or click on it to open it and view the contents.

```bash
cat reads_mother.vcf
```

If you scroll through, you'll find a header composed of many lines of metadata, followed by a list of variant calls, one per line.

#### 0.2.5. Exit the container

```bash
exit
```

### Takeaway

You know how to _DO THE THING_.

### What's next?

_DO THE NEXT THING_.

---

## 1. Write a single-stage workflow that runs Samtools index on a BAM file

#### 1.1. Define the indexing process

```groovy title="hello-gatk.nf"
/*
 * Generate BAM index file
 */
process SAMTOOLS_INDEX {

    container 'community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464'

    publishDir 'results', mode: 'copy'

    input:
        path input_bam

    output:
        path "${input_bam}.bai"

    """
    samtools index '$input_bam'
    """
}
```

#### 1.2. Add parameter declarations up top

```groovy title="hello-gatk.nf"
/*
 * Pipeline parameters
 */

// Primary input
params.reads_bam = "${projectDir}/data/bam/reads_mother.bam"
```

#### 1.3. Add workflow block to run SAMTOOLS_INDEX

```groovy title="hello-gatk.nf"
workflow {

    // Create input channel (single file via CLI parameter)
    reads_ch = Channel.fromPath(params.reads_bam)

    // Create index file for input BAM file
    SAMTOOLS_INDEX(reads_ch)
}
```

#### 1.4. Run it to verify you can run the indexing step

```bash
nextflow run hello-gatk.nf
```

Should produce something like:

```console title="Output"
N E X T F L O W  ~  version 24.02.0-edge
Launching `hello-gatk.nf` [compassionate_cray] DSL2 - revision: 9b97744397
executor >  local (1)
[bf/072bd7] process > SAMTOOLS_INDEX (1) [100%] 1 of 1 ✔
```

### Takeaway

You know how to wrap a real bioinformatics tool in a single-step Nextflow workflow.

### What's next?

Add a second step that consumes the output of the first.

---

## 2. Add a second step that runs GATK HaplotypeCaller on the indexed BAM file

#### 2.1. Define the variant calling process

```groovy title="hello-gatk.nf"
/*
 * Call variants with GATK HaplotypeCaller
 */
process GATK_HAPLOTYPECALLER {

    container "community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867"

    publishDir 'results', mode: 'copy'

    input:
        path input_bam
        path input_bam_index
        path ref_fasta
        path ref_index
        path ref_dict
        path interval_list

    output:
        path "${input_bam}.vcf"
        path "${input_bam}.vcf.idx"

    """
    gatk HaplotypeCaller \
        -R ${ref_fasta} \
        -I ${input_bam} \
        -O ${input_bam}.vcf \
        -L ${interval_list}
    """
}
```

#### 2.2. Add definitions for accessory inputs

```groovy title="hello-gatk.nf"
// Accessory files
params.reference        = "${projectDir}/data/ref/ref.fasta"
params.reference_index  = "${projectDir}/data/ref/ref.fasta.fai"
params.reference_dict   = "${projectDir}/data/ref/ref.dict"
params.intervals        = "${projectDir}/data/ref/intervals.bed"
```

#### 2.3. Make a value channel for each of the accessory files

Add this to the workflow block (after the `reads_ch` creation):

```groovy title="hello-gatk.nf"
// Create channels for the accessory files (reference and intervals)
ref_file        = file(params.reference)
ref_index_file  = file(params.reference_index)
ref_dict_file   = file(params.reference_dict)
intervals_file  = file(params.intervals)
```

This will load each of the accessory files in its own single-element value channel.

#### 2.4. Add a call to the workflow block to run GATK_HAPLOTYPECALLER

```groovy title="hello-gatk.nf"
// Call variants from the indexed BAM file
GATK_HAPLOTYPECALLER(
    reads_ch,
    SAMTOOLS_INDEX.out,
    ref_file,
    ref_index_file,
    ref_dict_file,
    intervals_file
)
```

#### 2.4. Run the workflow to verify that the variant calling step works

```bash
nextflow run hello-gatk.nf -resume
```

Now we see the two processes being run:

```console title="Output"
N E X T F L O W  ~  version 24.02.0-edge
Launching `hello-gatk.nf` [lethal_keller] DSL2 - revision: 30a64b9325
executor >  local (2)
[97/0f85bf] process > SAMTOOLS_INDEX (1)       [100%] 1 of 1 ✔
[2d/43c247] process > GATK_HAPLOTYPECALLER (1) [100%] 1 of 1 ✔
```

You'll find the output file `reads_mother.bam.vcf` in the results directory. If you open it, you should see the same contents as in the file you generated by running the GATK command directly in the container.

### Takeaway

You know how to make a very basic two-step variant calling workflow.

### What's next?

Make the workflow handle multiple samples in bulk.

---

## 3. Adapt the workflow to run on a batch of samples

#### 3.1. Turn the input parameter declaration into a list of the three samples

_Before:_

```groovy title="hello-gatk.nf"
// Primary input
params.reads_bam = "${projectDir}/data/bam/reads_mother.bam"
```

_After:_

```groovy title="hello-gatk.nf"
// Primary input (list of three samples)
params.reads_bam = [
    "${projectDir}/data/bam/reads_mother.bam",
    "${projectDir}/data/bam/reads_father.bam",
    "${projectDir}/data/bam/reads_son.bam"
]
```

#### 3.2. Run the workflow to verify that it runs on all three samples

```bash
nextflow run hello-gatk.nf -resume
```

Uh-oh! It will probably fail with an error like this:

```console title="Output"
ERROR ~ Error executing process > 'GATK_HAPLOTYPECALLER (2)'

Caused by:
  Process `GATK_HAPLOTYPECALLER (2)` terminated with an error exit status (2)

Command executed:

  gatk HaplotypeCaller         -R ref.fasta         -I reads_father.bam         -O reads_father.bam.vcf         -L intervals.bed         -ERC GVCF

Command exit status:
  2

Command error:
```

And buried in the GATK command error output, there will be a line like this:

```
A USER ERROR has occurred: Traversal by intervals was requested but some input files are not indexed.
```

This is because the script as written so far is not safe for running on multiple samples, because the order of items in the output channel is not guaranteed to match the order of items in the original input channel. This causes the wrong files to be paired up in the second step. So we need to make sure the BAM files and their index files travel together through the channels.

#### 3.3. Change the output of the SAMTOOLS_INDEX process into a tuple that keeps the input file and its index together

_Before:_

```groovy title="hello-gatk.nf"
output:
    path "${input_bam}.bai"
```

_After:_

```groovy title="hello-gatk.nf"
output:
    tuple path(input_bam), path("${input_bam}.bai")
```

#### 3.4. Change the input to the GATK_HAPLOTYPECALLER process to be a tuple

_Before:_

```groovy title="hello-gatk.nf"
input:
    path input_bam
    path input_bam_index
```

_After:_

```groovy title="hello-gatk.nf"
input:
    tuple path(input_bam), path(input_bam_index)
```

#### 3.5. Update the call to GATK_HAPLOTYPECALLER in the workflow block

_Before:_

```groovy title="hello-gatk.nf"
GATK_HAPLOTYPECALLER(
    reads_ch,
    SAMTOOLS_INDEX.out,
```

_After:_

```groovy title="hello-gatk.nf"
GATK_HAPLOTYPECALLER(
    SAMTOOLS_INDEX.out,
```

#### 3.6. Run the workflow to verify it works correctly on all three samples now

```bash
nextflow run hello-gatk.nf -ansi-log false
```

This time everything should run correctly:

```console title="Output"
N E X T F L O W  ~  version 24.02.0-edge
Launching `hello-gatk.nf` [adoring_hopper] DSL2 - revision: 8cad21ea51
[e0/bbd6ef] Submitted process > SAMTOOLS_INDEX (3)
[71/d26b2c] Submitted process > SAMTOOLS_INDEX (2)
[e6/6cad6d] Submitted process > SAMTOOLS_INDEX (1)
[26/73dac1] Submitted process > GATK_HAPLOTYPECALLER (1)
[23/12ed10] Submitted process > GATK_HAPLOTYPECALLER (2)
[be/c4a067] Submitted process > GATK_HAPLOTYPECALLER (3)
```

### Takeaway

You know how to make a variant calling workflow run on multiple samples (independently).

### What's next?

Make it easier to handle samples in bulk.

---

## 4. Make it nicer to run on arbitrary samples by using a list of files as input

#### 4.1. Create a text file listing the input paths

```csv title="sample_bams.txt"
/workspace/gitpod/hello-nextflow/data/bam/reads_mother.bam
/workspace/gitpod/hello-nextflow/data/bam/reads_father.bam
/workspace/gitpod/hello-nextflow/data/bam/reads_son.bam
```

#### 4.2. Update the parameter default

_Before:_

```groovy title="hello-gatk.nf"
// Primary input
params.reads_bam = [
    "${projectDir}/data/bam/reads_mother.bam",
    "${projectDir}/data/bam/reads_father.bam",
    "${projectDir}/data/bam/reads_son.bam"
]
```

_After:_

```groovy title="hello-gatk.nf"
// Primary input (list of input files, one per line)
params.reads_bam = "${projectDir}/data/sample_bams.txt"
```

#### 4.3. Update the channel factory to read lines from a file

_Before:_

```groovy title="hello-gatk.nf"
// Create input channel (single file via CLI parameter)
reads_ch = Channel.fromPath(params.reads_bam)
```

_After:_

```groovy title="hello-gatk.nf"
// Create input channel from list of input files in plain text
reads_ch = Channel.fromPath(params.reads_bam).splitText()
```

#### 4.4. Run the workflow to verify that it works correctly

```bash
nextflow run hello-gatk.nf -resume -ansi-log false
```

This should produce essentially the same result as before:

```console title="Output"
N E X T F L O W  ~  version 24.02.0-edge
Launching `hello-gatk.nf` [backstabbing_raman] DSL2 - revision: 5378632b71
[56/5f8548] Cached process > SAMTOOLS_INDEX (1)
[69/7a0ce1] Cached process > SAMTOOLS_INDEX (3)
[f1/6ae50e] Cached process > SAMTOOLS_INDEX (2)
[49/abf910] Cached process > GATK_HAPLOTYPECALLER (1)
[03/4407cd] Cached process > GATK_HAPLOTYPECALLER (2)
[af/2ea71a] Cached process > GATK_HAPLOTYPECALLER (3)
```

### Takeaway

You know how to make a multi-step variant calling workflow handle a file containing input samples.

### What's next?

Celebrate your success and take an extra long break! 

In the next training module, you'll learn how to _DO MORE THINGS_.

**Good luck!**
