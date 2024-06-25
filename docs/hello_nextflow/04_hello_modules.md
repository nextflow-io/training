# Part 3: Hello Modules

In Nextflow, a module is a way to encapsulate a single process by itself in a standalone code file. To use a module in a workflow, you just add a single-line import statement to your workflow code file; then you can integrate the process into the workflow the same way you normally would.

Putting processes into individual modules makes it possible to reuse process descriptions in multiple workflows instead of just replicating the code. This makes the code more shareable, flexible and maintainable.

Note: It is also possible to encapsulate a section of a workflow as a 'subworkflow' that can be imported into a larger pipeline, but that is outside the scope of this training.

---

## 0. Warmup

We start from a base workflow called `hello-modules.nf`, which corresponds to the workflow we produced in Part 2: Hello
GATK (equivalent to `scripts/hello-gatk-6.nf`).

Note: This is a basic variant calling pipeline consisting of three processes. You can find a complete description of the pipeline in the previous section of this training.

This workflow relies on reference files that are provided in compressed form in the Gitpod environment. If you completed the previous parts of the training course, then you already have everything you need in the working directory. However, if you're picking this up here, you need to run the following command to expand the reference files:

```bash
cd /workspace/gitpod/hello-nextflow
tar -zxvf data/ref.tar.gz -C data/
```

### 0.1 Run the workflow to verify that it produces the expected outputs

```bash
nextflow run hello-modules.nf
```

The pipeline takes in three BAM files, each one containing sequencing data for one of three samples from a human family trio (mother, father and son), and outputs a VCF file containing variant calls. For more details, see the previous section of this training.

---

## 1. Move parameter declarations from the workflow file to a config file

### 1.1 Move the parameter definitions into `nextflow.config`

```groovy
/*
 * Pipeline parameters
 */

// Execution environment setup
params.projectDir = "/workspace/gitpod/hello-nextflow"
$projectDir = params.projectDir

// Primary input (samplesheet in CSV format with ID and file path, one sample per line)
params.reads_bam = "${projectDir}/data/samplesheet.csv"

// Accessory files
params.genome_reference = "${projectDir}/data/ref/ref.fasta"
params.genome_reference_index = "${projectDir}/data/ref/ref.fasta.fai"
params.genome_reference_dict = "${projectDir}/data/ref/ref.dict"
params.calling_intervals = "${projectDir}/data/intervals.list"

// Base name for final output file
params.cohort_name = "family_trio"
```

### 1.2 Run the workflow to verify that it does the same thing as before

```bash
nextflow run hello-modules.nf
```

You should see the same output as earlier.

Note: Workflow parameters and other elements of configuration can be specified in several different places. When the same element is defined in more than one place, the conflict is resolved according to the order of precedence detailed in [this documentation](https://nextflow.io/docs/latest/config.html#configuration-file).

### Takeaway

You know how to move parameter definitions to a configuration file.

### What's next?

Learn how to extract processes into local modules for code reuse.

---

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
