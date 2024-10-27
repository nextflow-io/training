# Part 6: Hello Modules

This section covers how to organize your workflow code to make development and maintenance of your pipeline more efficient and sustainable.
Specifically, we are going to demonstrate how to use **modules**.

In Nextflow, a **module** is a single process definition that is encapsulated by itself in a standalone code file.
To use a module in a workflow, you just add a single-line import statement to your workflow code file; then you can integrate the process into the workflow the same way you normally would.

Putting processes into individual modules makes it possible to reuse process definitions in multiple workflows without producing multiple copies of the code.
This makes the code more shareable, flexible and maintainable.

!!!note

    It is also possible to encapsulate a section of a workflow as a 'subworkflow' that can be imported into a larger pipeline, but that is outside the scope of this training.

---

## 0. Warmup

When we started developing our workflow, we put everything in one single code file.
In Part 5, we started turning our one-script workflow into a proper pipeline project by teasing out parameter settings and building out the configuration file(s).
We also moved to the standard Nextflow convention of naming the workflow file `main.nf`.

Now it's time to tackle **modularizing** our code.

We're going to pick up where we left off in Part 5 (Hello Config), this time working inside the project directory called `projectM` (for Modules).

### 0.1. Explore the projectM directory¶

Let's move into the project directory.
If you're continuing on directly from Part 5, you'll need to move up one directory first.

```bash
cd projectM
```

The `projectM` directory has the same content and structure that you're expected to end up with on completion of Part 5.

```console title="Directory contents"
projectC/
├── demo-params.json
├── intermediates
├── main.nf
└── nextflow.config
```

For a detailed description of these files, see the warmup in Part 5.

### 0.2. Create a symbolic link to the data¶

Run this command from inside the projectC directory:

```bash
ln -s ../data data
```

This creates a symbolic link called data pointing to the data directory, which allows us to avoid having to change anything to how the file paths are set up.

### 0.3 Run the workflow using the appropriate profiles

Now that everything is in place, we should be able to run the workflow successfully using the profiles we set up in Part 5.

```bash
nextflow run main.nf -profile my_laptop,demo
```

This should run successfully:

```console title="Output"
**TODO: rerun & add updated output**
```

There will now be a `work` directory and a `results_genomics` directory inside your project directory.

---

**TODO: update numbering and code snippets**

## 2. Create a module for the `SAMTOOLS_INDEX` process

By convention, a process module should be written to a standalone file named `main.nf` and stored in a directory structure named after the tool (and optionally, toolkit) that the process invokes. For example, the module we create for the `SAMTOOLS_INDEX` process will live under `samtools/index/`.

Additionally, local modules are grouped under the directory structure `modules/local/`. Modules can also belong to remote repositories; we will not cover those in this training.

### 2.1 Create a folder to house the local module code for the `SAMTOOLS_INDEX` process

```bash
mkdir -p modules/local/samtools/index
```

### 2.2 Create a file stub for the `SAMTOOLS_INDEX` process module

```bash
touch modules/local/samtools/index/main.nf
```

This creates an empty file named `main.nf` under the appropriate directory structure.

### 2.3 Move the `SAMTOOLS_INDEX` process code to the module file

```groovy
/*
 * Generate BAM index file
 */
process SAMTOOLS_INDEX {

    container 'quay.io/biocontainers/samtools:1.19.2--h50ea8bc_1'

    input:
        tuple val(id), path(input_bam)

    output:
        tuple val(id), path(input_bam), path("${input_bam}.bai")

    """
    samtools index '$input_bam'
    """
}
```

### 2.4 Add an import declaration before the workflow block

_Before:_

```groovy
workflow {
```

_After:_

```groovy
// Include modules
include { SAMTOOLS_INDEX } from './modules/local/samtools/index/main.nf'

workflow {
```

### 2.5 Run the workflow to verify that it does the same thing as before

```bash
nextflow run hello-modules.nf
```

### Takeaway

You know how to extract a process into a local module.

### What's next?

Practice making more modules.

---

## 3. Repeat procedure for the remaining processes

### 3.1 Create folders to house the local module code

```bash
mkdir -p modules/local/gatk/haplotypecaller
mkdir -p modules/local/gatk/jointgenotyping
```

Here we are grouping the two modules under `gatk/` because the corresponding processes both invoke tools that belong to the GATK toolkit.

### 3.2 Create file stubs for the process modules

```bash
touch modules/local/gatk/haplotypecaller/main.nf
touch modules/local/gatk/jointgenotyping/main.nf
```

### 3.3 Move the process code to the module files

Move this code to `modules/local/gatk/haplotypecaller/main.nf`:

```groovy
/*
 * Call variants with GATK HaplotypeCaller in GVCF mode
 */
process GATK_HAPLOTYPECALLER {

    container "docker.io/broadinstitute/gatk:4.5.0.0"

    input:
        tuple val(id), path(input_bam), path(input_bam_index)
        path ref_fasta
        path ref_index
        path ref_dict
        path interval_list

    output:
        tuple val(id), path("${input_bam}.g.vcf"), path("${input_bam}.g.vcf.idx")

    """
    gatk HaplotypeCaller \
        -R ${ref_fasta} \
        -I ${input_bam} \
        -O ${input_bam}.g.vcf \
        -L ${interval_list} \
        -ERC GVCF
    """
}
```

And move this code to `modules/local/gatk/jointgenotyping/main.nf`:

```groovy
/*
 * Consolidate GVCFs and apply joint genotyping analysis
 */
process GATK_JOINTGENOTYPING {

    container "docker.io/broadinstitute/gatk:4.5.0.0"

    input:
        path(sample_map)
        val(cohort_name)
        path ref_fasta
        path ref_index
        path ref_dict
        path interval_list

    output:
        path "${cohort_name}.joint.vcf"
        path "${cohort_name}.joint.vcf.idx"

    """
    gatk GenomicsDBImport \
        --sample-name-map ${sample_map} \
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

### 3.4 Add import declarations before the workflow block

_Before:_

```groovy
// Include modules
include { SAMTOOLS_INDEX } from './modules/local/samtools/index/main.nf'

workflow {
```

_After:_

```groovy
// Include modules
include { SAMTOOLS_INDEX } from './modules/local/samtools/index/main.nf'
include { GATK_HAPLOTYPECALLER } from './modules/local/gatk/haplotypecaller/main.nf'
include { GATK_JOINTGENOTYPING } from './modules/local/gatk/jointgenotyping/main.nf'

workflow {
```

### 3.5 Run the workflow to verify that it does the same thing as before

```bash
nextflow run hello-modules.nf
```

This should still produce the same output. Congratulations, you've done all this work and absolutely nothing has changed to how the pipeline works! But now your code is more modular, and if you decide to write another pipeline that calls on one of those processes, you just need to type one short import statement to use the relevant module. This is better than just copy-pasting the code, because if later you decide to improve the module, all your pipelines will inherit the improvement.

### Takeaway

You know how to modularize multiple processes in a workflow.

### What's next?

Learn to add tests to your pipeline using the nf-test framework.
