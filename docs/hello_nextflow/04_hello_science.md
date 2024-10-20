# Part 3: Hello Science

In Part 1, you learned how to use the basic building blocks of Nextflow to assemble a simple pipeline capable of processing some text and parallelizing execution if there were multiple inputs. Now, we show you how to use the same components and principles to build a pipeline that does something a bit more interesting, and hopefully a bit more relatable to your work.

Specifically, we show you how to implement a simple variant calling pipeline with [GATK](https://gatk.broadinstitute.org/) (Genome Analysis Toolkit), a widely used software package for analyzing high-throughput sequencing data.

!!! note

    Don't worry if you're not familiar with GATK or genomics in general. We'll summarize the necessary concepts as we go, and the workflow implementation principles we demonstrate here apply broadly to any command line tool that takes in some input files and produce some output files.

### Method overview

Variant calling is a genomic analysis method that aims to identify variations in a genome sequence relative to a reference genome. Here we are going to use tools and methods designed for calling short variants, _i.e._ SNPs and indels.

A full variant calling pipeline typically involves a lot of steps, including mapping to the reference and variant filtering and prioritization. For simplicity, we are going to focus on the core variant calling step, which takes as its main input a file of short-read sequencing data in BAM format (Binary-compressed version of SAM, for Sequence Alignment Map), as well as a reference genome and a list of genomic intervals to analyze.

![GATK pipeline](img/gatk-pipeline.png)

For this exercise, we provide you with three samples in BAM format (see Dataset below). However, GATK requires an index file for each BAM file, which we did not provide (on purpose), so the workflow will have to create one as a preliminary step.

!!! note

    Index files are a common feature of bioinformatics file formats; they contain information about the structure of the main file that allows tools like GATK to access a subset of the data without having to read through the whole file. This is important because of how big these files can get.

So to recap, we're going to develop a workflow that does the following:

1. Generate an index file for each BAM input file using [Samtools](https://www.htslib.org/)
2. Run the GATK HaplotypeCaller on each BAM input file to generate per-sample variant calls in VCF (Variant Call Format)

<figure class="excalidraw">
--8<-- "docs/hello_nextflow/img/haplotype-caller.excalidraw.svg"
</figure>

### Dataset

-   **A reference genome** consisting of a small region of the human chromosome 20 (from hg19/b37) and its accessory files (index and sequence dictionary).
-   **Three whole genome sequencing samples** corresponding to a family trio (mother, father and son), which have been subset to a small portion on chromosome 20 to keep the file sizes small. The sequencing data is in [BAM](https://samtools.github.io/hts-specs/SAMv1.pdf) (Binary Alignment Map) format, i.e. genome sequencing reads that have already been mapped to the reference genome.
-   **A list of genomic intervals**, i.e. coordinates on the genome where our samples have data suitable for calling variants, provided in BED format.

---

## 0. Warmup: Test the Samtools and GATK commands interactively

Just like in the Hello World example, we want to try out the commands manually before we attempt to wrap them in a workflow. The tools we need (Samtools and GATK) are not installed in the Gitpod environment, but that's not a problem since you learned how to work with containers in Part 2 of this training series (Hello Containers).

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

The [Samtools documentation](https://www.htslib.org/doc/samtools-index.html) gives us the command line to run to index a BAM file.

We only need to provide the input file; the tool will automatically generate a name for the output by appending `.bai` to the input filename.

```bash
samtools index /data/bam/reads_mother.bam
```

This should complete immediately, and you should now see a file called `reads_mother.bam.bai` in the same directory as the original BAM input file.

#### 0.1.4. Exit the container

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

The [GATK documentation](https://gatk.broadinstitute.org/hc/en-us/articles/21905025322523-HaplotypeCaller) gives us the command line to run to perform variant calling on a BAM file.

We need to provide the BAM input file (`-I`) as well as the reference genome (`-R`), a name for the output file (`-O`) and a list of genomic intervals to analyze (`-L`).

However, we don't need to specify the path to the index file; the tool will automatically look for it in the same directory, based on the established naming and co-location convention. The same applies to the reference genome's accessory files (index and sequence dictionary files).

```bash
gatk HaplotypeCaller \
        -R /data/ref/ref.fasta \
        -I /data/bam/reads_mother.bam \
        -O reads_mother.vcf \
        -L /data/ref/intervals.bed
```

The output file `reads_mother.vcf` is a small test file, so you can `cat` it or click on it to open it and view the contents. If you scroll through, you'll find a header composed of many lines of metadata, followed by a list of variant calls, one per line.

TODO _SHOW A FEW LINES_

#### 0.2.4. Exit the container

```bash
exit
```

### Takeaway

You know how to test the Samtools indexing and GATK variant calling commands in their respective containers.

### What's next?

Learn how to wrap those same commands into a two-step workflow that uses containers to execute the work.

---

## 1. Write a single-stage workflow that runs Samtools index on a BAM file

We provide you with a workflow file, `hello-gatk.nf`, that outlines the main parts of the workflow. It's not functional; its purpose is just to serve as a skeleton that you'll use to write the actual workflow.

### 1.1. Define the indexing process

Let's start by writing a process, which we'll call `SAMTOOLS_INDEX`, describing the indexing operation.

```groovy title="hello-gatk.nf"
/*
 * Generate BAM index file
 */
process SAMTOOLS_INDEX {

    container 'community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464'

    publishDir 'results_genomics', mode: 'symlink'

    input:
        path input_bam

    output:
        path "${input_bam}.bai"

    """
    samtools index '$input_bam'
    """
}
```

You should recognize all the pieces from what you learned in Part 1 & Part 2 of this training series; the only notable change is that this time we're using `mode: symlink` for the `publishDir` directive.

!!! note
Even though the data files we're using here are very small, in genomics they can get very large, so we should get into the habit of using symbolic links rather than making actual copies of these files, unless there's a compelling reason to do so.

This process is going to require passing in a filepath via the `input_bam` input, so let's set that up next.

### 1.2. Add an input parameter declaration

At the top of the file, under the `Pipeline parameters` section, we declare a CLI parameter called `reads_bam` and give it a default value. That way, we can be lazy and not specify the input when we type the command to launch the pipeline (for development purposes).

```groovy title="hello-gatk.nf"
/*
 * Pipeline parameters
 */

// Primary input
params.reads_bam = "${projectDir}/data/bam/reads_mother.bam"
```

Now we have a process ready, and a parameter to give it an input to run on, so let's wire those things up together.

### 1.3. Add workflow block to run SAMTOOLS_INDEX

In the `workflow` block, we need to set up a **channel** to feed the input to the `SAMTOOLS_INDEX` process; then we can call the process itself to run on the contents of that channel.

```groovy title="hello-gatk.nf"
workflow {

    // Create input channel (single file via CLI parameter)
    reads_ch = Channel.fromPath(params.reads_bam)

    // Create index file for input BAM file
    SAMTOOLS_INDEX(reads_ch)
}
```

You'll notice we're using the same `.fromPath` channel constructor as we used at the end of Part 1 (Hello World) of this training series. Indeed, we're doing something very similar; the difference is that this time we're telling Nextflow to load the filepath itself into the channel as an input element, rather than reading in its contents.

### 1.4. Run the workflow to verify that the indexing step works

Let's run the workflow! As a reminder, we don't need to specify an input in the command line because we set up a default value for the input when we declared the input parameter.

```bash
nextflow run hello-gatk.nf
```

The command should produce something like this:

```console title="Output"
N E X T F L O W  ~  version 24.02.0-edge
Launching `hello-gatk.nf` [compassionate_cray] DSL2 - revision: 9b97744397
executor >  local (1)
[bf/072bd7] process > SAMTOOLS_INDEX (1) [100%] 1 of 1 ✔
```

You can check that the index file has been generated correctly by looking in the work directory or in the directory set up with `publishDir`.

### Takeaway

You know how to wrap a real bioinformatics tool in a single-step Nextflow workflow and have it run using a container.

### What's next?

Add a second step that consumes the output of the first.

---

## 2. Add a second process to run GATK HaplotypeCaller on the indexed BAM file

Now that we have an index for our input file, we can move on to setting up the variant calling step, which is the interesting part of the workflow.

### 2.1. Define the variant calling process

Let's write a process, which we'll call `GATK_HAPLOTYPECALLER`, describing the variant calling operation.

```groovy title="hello-gatk.nf"
/*
 * Call variants with GATK HaplotypeCaller
 */
process GATK_HAPLOTYPECALLER {

    container "community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867"

    publishDir 'results_genomics', mode: 'symlink'

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

This command takes quite a few more inputs, because GATK needs more information to perform the analysis compared to a simple indexing job. But you'll note that there are even more inputs defined in the inputs block than are listed in the GATK command. Why is that?

!!! note
The GATK knows to look for the BAM index file and the reference genome's accessory files because it is aware of the conventions surrounding those files. However, Nextflow is designed to be domain-agnostic and doesn't know anything about bioinformatics file format requirements. So we need to tell it explicitly that it has to stage those files in the working directory at runtime; otherwise it won't do it, and GATK will (correctly) throw an error about the index files being missing.

    Similarly, we have to list the output VCF's index file (the `"${input_bam}.vcf.idx"` file) explicitly so that Nextflow will know to keep track of that file in case it's needed in subsequent steps.

### 2.2. Add definitions for accessory inputs

Since our new process expects a handful of additional files to be provided, we set up some CLI parameters for them under the `Pipeline parameters` section, along with some default values (same reasons as before).

```groovy title="hello-gatk.nf"
// Accessory files
params.reference        = "${projectDir}/data/ref/ref.fasta"
params.reference_index  = "${projectDir}/data/ref/ref.fasta.fai"
params.reference_dict   = "${projectDir}/data/ref/ref.dict"
params.intervals        = "${projectDir}/data/ref/intervals.bed"
```

### 2.3. Create variables to hold the accessory file paths

Unlike the main data inputs, which must be fed to processes through channels, the accessory files can be handled a bit more simply: we can use the `file()` function to create variables to hold those file paths.

Add this to the workflow block (after the `reads_ch` creation):

```groovy title="hello-gatk.nf"
// Load the file paths for the accessory files (reference and intervals)
ref_file        = file(params.reference)
ref_index_file  = file(params.reference_index)
ref_dict_file   = file(params.reference_dict)
intervals_file  = file(params.intervals)
```

This will make the accessory file paths available for input to processes.

### 2.4. Add a call to the workflow block to run GATK_HAPLOTYPECALLER

Now that we've got our second process set up and all the inputs and accessory files are ready and available, we can add a call to the `GATK_HAPLOTYPECALLER` process in the workflow body.

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

You should recognize the `SAMTOOLS_INDEX.out` syntax from Part 1 of this training series; we are telling Nextflow to take the channel output by `SAMTOOLS_INDEX` and plugging that into the `GATK_HAPLOTYPECALLER` process.

!!! note
You'll notice that the inputs are provided in the exact same order in the call to the process as they are listed in the process' input block. In Nextflow, inputs are positional, meaning you _must_ follow the same order; and of course there have to be the same number of elements.

### 2.5. Run the workflow to verify that the variant calling step works

Let's run the expanded workflow with `-resume` so that we don't have to run the indexing step again.

```bash
nextflow run hello-gatk.nf -resume
```

Now if we look at the console output, we see the two processes being run:

```console title="Output"
N E X T F L O W  ~  version 24.02.0-edge
Launching `hello-gatk.nf` [lethal_keller] DSL2 - revision: 30a64b9325
executor >  local (2)
[97/0f85bf] process > SAMTOOLS_INDEX (1)       [100%] 1 of 1 ✔
[2d/43c247] process > GATK_HAPLOTYPECALLER (1) [100%] 1 of 1 ✔
```

You'll find the output file `reads_mother.bam.vcf` in the results directory. If you open it, you should see the same contents as in the file you generated by running the GATK command directly in the container.

### Takeaway

You know how to make a very basic two-step workflow that does real analysis work and is capable of dealing with bioinformatics idiosyncrasies like the accessory files.

### What's next?

Make the workflow handle multiple samples in bulk.

---

## 3. Adapt the workflow to run on a batch of samples

It's all well and good to have a workflow that can automate processing on a single sample, but what if you have 1000 samples? Do you need to write a bash script that loops through all your samples?

No, thank goodness! Just make a minor tweak to the code and Nextflow will handle that for you too.

### 3.1. Turn the input parameter declaration into a list of the three samples

Let's turn that default file path in the input BAM file declaration into a list of file paths, up under the `Pipeline parameters` section.

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

And that's actually all you need to do, because the channel constructor we use in the workflow body (`.fromPath`) is just as happy to accept multiple file paths to load into the input channel as it was to load a single one.

!!! note
Normally, you wouldn't want to hardcode the list of samples into your workflow file, but we're doing that here to keep things simple. We'll go over more elegant ways of handling inputs later in this training series.

### 3.2. Run the workflow to verify that it runs on all three samples

Let's try running the workflow now that the plumbing is set up to run on all three test samples.

```bash
nextflow run hello-gatk.nf -resume
```

Funny thing: this might work, OR it might fail with an error like this:

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

Well, that's weird, considering we explicitly indexed the BAM files in the first step of the workflow. Could there be something wrong with the plumbing? Let's check the work directories for the relevant calls.

TODO _Add snippet of output showing mismatch between main file and index_

What the heck? Nextflow has staged an index file in this process call's work directory, but it's the wrong one. How could this have happened?

!!! note
When you call a Nextflow process on a channel containing multiple elements, Nextflow will try to parallelize execution as much as possible. The consequence is that the corresponding outputs may be collected in a different order than the original inputs were fed in.

As currently written, our workflow script assumes that the index files will come out of the indexing step listed in the same mother/father/son order as the inputs were given. But that is not guaranteed to be the case, which is why sometimes (though not always) the wrong files get paired up in the second step.

To fix this, we need to make sure the BAM files and their index files travel together through the channels.

### 3.3. Change the output of the SAMTOOLS_INDEX process into a tuple that keeps the input file and its index together

The simplest way to ensure a BAM file and its index stay closely associated is to package them together into a tuple coming out of the index task.

!!! note
A **tuple** is a finite, ordered list of elements that is commonly used for returning multiple values from a function.

First, let's change the output of the `SAMTOOLS_INDEX` process to include the BAM file in its output declaration.

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

This way, each index file will be tightly coupled with its original BAM file, and the overall output of the indexing step will be a single channel containing pairs of files.

### 3.4. Change the input to the GATK_HAPLOTYPECALLER process to be a tuple

Since we've changed the 'shape' of the output of the first process in the workflow, we need to update the input definition of the second process to match.

Specifically, where we previously declared two separate input paths in the input block of the `GATK_HAPLOTYPECALLER` process, we now declare a single input matching the structure of the tuple emitted by `SAMTOOLS_INDEX`.

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

Of course, since we've now changed the shape of the inputs that `GATK_HAPLOTYPECALLER` expects, we need to update the process call accordingly in the workflow body.

### 3.5. Update the call to GATK_HAPLOTYPECALLER in the workflow block

We no longer need to provide the original `reads_ch` to the `GATK_HAPLOTYPECALLER` process, since the BAM file is now bundled (in the form of a symlink) into the channel output by `SAMTOOLS_INDEX`.

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

### 3.6. Run the workflow to verify it works correctly on all three samples now

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

### 4.1. Create a text file listing the input paths

```csv title="sample_bams.txt"
/workspace/gitpod/hello-nextflow/data/bam/reads_mother.bam
/workspace/gitpod/hello-nextflow/data/bam/reads_father.bam
/workspace/gitpod/hello-nextflow/data/bam/reads_son.bam
```

### 4.2. Update the parameter default

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

### 4.3. Update the channel factory to read lines from a file

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

### 4.4. Run the workflow to verify that it works correctly

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

You know how to make a multi-step workflow handle a file containing input samples.

### What's next?

Celebrate your success and take an extra long break!

In the next training module, you'll learn how to use channel operators to develop pipelines with more interesting plumbing.

**Good luck!**
