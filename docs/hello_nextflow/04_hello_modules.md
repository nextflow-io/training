# Part 3: Hello Modules

[short blurb about modularizing a pipeline]
[note you can also do subworkflows but that's out of scope of this training]

---

## 0. Warmup

We start from a base workflow called `hello-modules.nf`, which corresponds to the workflow we produced in Part 2: Hello 
GATK (equivalent to `scripts/hello-gatk-6.nf`). 

[TODO: add brief description of what the pipeline does]

This workflow relies on reference file that are provided in compressed form in the Gitpod environment. If you completed the previous parts of the training course, then you already have everything you need in the working directory. However, if you're picking this up here, you need to run the following command to expand the reference files:

```bash
tar -zxvf data/ref.tar.gz -C data/
```

### 0.1 Run the workflow to verify that it produces the expected outputs 

```bash
nextflow run hello-modules.nf
```

[TODO: briefly describe outputs; refer to Part 2 for full details]

---

## 1. Move parameter declarations from the workflow file to a config file

### 1.1 Move the parameter definitions into `nextflow.config`

```groovy
/*
 * Pipeline parameters
 */

// Execution environment setup
params.baseDir = "/workspace/gitpod/hello-nextflow"
$baseDir = params.baseDir

// Primary input (samplesheet in CSV format with ID and file path, one sample per line)
params.reads_bam = "${baseDir}/data/samplesheet.csv"

// Base name for final output file
params.cohort_name = "family_trio"

// Accessory files
params.genome_reference = "${baseDir}/data/ref/ref.fasta"
params.genome_reference_index = "${baseDir}/data/ref/ref.fasta.fai"
params.genome_reference_dict = "${baseDir}/data/ref/ref.dict"
params.calling_intervals = "${baseDir}/data/intervals.list"
```

### 1.2 Run the workflow to verify that it does the same thing as before

```bash
nextflow run hello-modules.nf
```

[TODO: Note where else config elements can be found and mention the order of precedence]


### Takeaway

You know how to move parameter definitions to a configuration file.

### What's next?

Learn how to extract processes into local modules for code reuse.

---

## 2. Create a module for the `SAMTOOLS_INDEX` process

[TODO: Briefly explain that for local modules, there should be a directory for each process, containing a file called `main.nf`; folder structure follows the toolkit/tool naming convention]

[TODO: Note non-local, nf-core cases exist but no details; covered later]

### 2.1 Create a folder to house the local module code for the `SAMTOOLS_INDEX` process

```bash
mkdir -p modules/local/samtools/index
```

### 2.2 Create a file stub for the `SAMTOOLS_INDEX` process module

```bash
touch modules/local/samtools/index/main.nf
```

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

### 2.4 Add an import declaration before the workflow 

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

    container "broadinstitute/gatk:4.5.0.0"

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

    container "broadinstitute/gatk:4.5.0.0"

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

### 3.4 Add import declarations before the workflow 

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

### Takeaway

You know how to modularize an entire workflow.

### What's next?

Learn to add tests to your pipeline using the nf-test framework.

