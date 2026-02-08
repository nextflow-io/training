# Part 4: Run nf-core/sarek

In Parts 2 and 3, you built a variant calling pipeline from scratch using GATK tools.
This gave you a solid understanding of how Nextflow processes work and how to control data flow.

In this part, we introduce [**nf-core/sarek**](https://nf-co.re/sarek), a production-ready variant calling pipeline from the nf-core community.
Sarek implements a more comprehensive version of the variant calling workflow you built, with additional features, quality control, and support for multiple variant callers.

## 1. Understanding nf-core pipelines

Before we run Sarek, here is some background on what nf-core is and why it matters.

### 1.1. What is nf-core?

[nf-core](https://nf-co.re/) is a community-driven collection of high-quality Nextflow pipelines.
All nf-core pipelines follow the same structure and conventions, which means once you learn to run one, you can run any of them.

Key features of nf-core pipelines:

- **Standardized structure**: All pipelines have consistent parameter names and usage patterns
- **Built-in test data**: Every pipeline includes test profiles for quick validation
- **Comprehensive documentation**: Detailed usage instructions and parameter descriptions
- **Quality control**: Automated QC reports using MultiQC
- **Container support**: Pre-built containers for reproducibility

!!! tip "Want to learn more about nf-core?"

    For an in-depth introduction to nf-core pipeline development, check out the [Hello nf-core](../../hello_nf-core/index.md) training course.
    It covers how to create and customize nf-core pipelines from scratch.

### 1.2. The Sarek pipeline

![nf-core/sarek pipeline](https://raw.githubusercontent.com/nf-core/sarek/master/docs/images/sarek_workflow.png)

The [nf-core/sarek](https://nf-co.re/sarek) pipeline is designed for detecting germline and somatic variants from whole genome or targeted sequencing data.
It supports humans, mice, and other species with reference genomes.

The pipeline processes data through several stages:

1. **Preprocessing**: Quality control, read mapping, duplicate marking, and base quality score recalibration (BQSR)
2. **Variant calling**: Multiple tools available including GATK HaplotypeCaller, DeepVariant, Strelka, FreeBayes, and more
3. **Annotation**: Variant annotation using SnpEff or Ensembl VEP
4. **Quality control**: Comprehensive QC reports via MultiQC

Key outputs include:

- Aligned and processed BAM/CRAM files
- Variant calls in VCF format
- Annotated variants
- MultiQC quality control report

---

## 2. Run Sarek with the built-in test data

nf-core pipelines include built-in test profiles that use small datasets to validate the pipeline works correctly.
This is the recommended way to verify your environment is set up properly.

!!! note

    Make sure you're in the correct working directory:
    `cd /workspaces/training/nf4-science/genomics`

### 2.1. Understanding the test profile

[TODO: brief intro text]

#### 2.1.1. Container requirements

Like most nf-core pipelines, Sarek requires containerization (Docker, Singularity, or similar) because it uses many specialized bioinformatics tools.

When you run an nf-core pipeline, each process specifies which container image provides its required software.
However, Nextflow only uses these containers if you tell it to via a profile or configuration.

#### 2.1.2. Test profile configuration

The test profile configures Sarek with:

- A small subset of sequencing data (to run quickly)
- Minimal resource requirements
- A reduced reference genome
- Default variant caller (Strelka)

You can inspect what the test profile sets by looking at the pipeline's configuration:

```bash
nextflow config nf-core/sarek -r 3.5.1 -profile test
```

[TODO: add a collapsible 'success' admonition showing the command output]

This shows all parameters that the test profile configures, including the input samplesheet URL and reference genome settings.

### 2.2. Run the test profile

The simplest way to run Sarek is using its built-in test profile, which includes a small test dataset.
To that end, we're going to give the following arguments to the `nextflow run` command:

- `nf-core/sarek`: The pipeline to run (Nextflow downloads it from GitHub automatically)
- `-r 3.5.1`: The pipeline version (pinning versions ensures reproducibility)
- `-profile test,docker`: Use the test profile AND enable Docker containers
- `--outdir results_sarek`: Where to save the results

Putting it all together, we get this:

```bash
nextflow run nf-core/sarek -r 3.5.1 -profile test,docker --outdir results_sarek
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 24.10.0

    Launching `https://github.com/nf-core/sarek` [friendly_torvalds] DSL2 - revision: 47a59ed14e [3.5.1]

    ------------------------------------------------------
                                            ,--./,-.
            ___     __   __   __   ___     /,-._.--~'
      |\ | |__  __ /  ` /  \ |__) |__         }  {
      | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                            `._,._,'
      nf-core/sarek 3.5.1
    ------------------------------------------------------
    ...
    executor >  local (X)
    [ab/123456] NFCORE_SAREK:SAREK:FASTQC                    [100%] 2 of 2 ✔
    [cd/789012] NFCORE_SAREK:SAREK:FASTP                     [100%] 2 of 2 ✔
    [ef/345678] NFCORE_SAREK:SAREK:BWAMEM2_MEM               [100%] 2 of 2 ✔
    ...
    -[nf-core/sarek] Pipeline completed successfully-
    ```

Note that the first time you run Sarek, Nextflow will download the pipeline code and Docker containers.
This can take several minutes depending on your internet connection.

### 2.3. Explore the outputs

When the pipeline completes, you'll find results organized in the output directory.

#### 2.3.1. Output directory structure

List the contents of the results directory:

```bash
ls results_sarek/
```

You should see directories like:

```console
results_sarek/
├── csv/
├── multiqc/
├── pipeline_info/
├── preprocessing/
├── reports/
└── variant_calling/
```

Key directories:

- **csv/**: Samplesheets for resuming from different pipeline stages
- **preprocessing/**: Aligned and processed BAM files, recalibration tables
- **variant_calling/**: VCF files with variant calls
- **reports/**: Quality control reports from various tools
- **multiqc/**: Aggregated quality control report
- **pipeline_info/**: Execution reports and logs

#### 2.3.2. Examine the MultiQC report

The MultiQC report aggregates quality metrics from all pipeline steps.
Open it in your browser or use VS Code's preview feature:

```bash
ls results_sarek/multiqc/
```

The report includes:

- FastQC metrics for raw reads
- Alignment statistics
- Duplicate rates
- Coverage metrics
- Variant calling statistics

#### 2.3.3. Find the variant calls

The variant calls are in the `variant_calling/` directory, organized by variant caller:

```bash
ls results_sarek/variant_calling/
```

For the test profile (which uses Strelka by default), you'll find VCF files containing the detected variants.

#### 2.3.4. Pipeline execution reports

Check the pipeline_info directory for execution details:

```bash
ls results_sarek/pipeline_info/
```

This includes:

- **execution_report.html**: Resource usage and timeline
- **execution_trace.txt**: Detailed task metrics
- **pipeline_dag.html**: Workflow structure diagram

### Takeaway

nf-core pipelines produce well-organized outputs with comprehensive quality control reports.

### What's next?

Learn how to configure Sarek for different use cases.

---

## 3. Configure Sarek for your data

The test profile is useful for validation, but running Sarek on your own data requires additional configuration.

### 3.1. Key parameters

Sarek has many configurable parameters. Here are the most commonly used:

| Parameter     | Description                 | Example                                        |
| ------------- | --------------------------- | ---------------------------------------------- |
| `--input`     | Path to input samplesheet   | `samplesheet.csv`                              |
| `--outdir`    | Output directory            | `results/`                                     |
| `--genome`    | Reference genome            | `GATK.GRCh38`, `GRCh37`                        |
| `--tools`     | Variant callers to use      | `haplotypecaller,strelka`                      |
| `--step`      | Pipeline starting point     | `mapping`, `markduplicates`, `variant_calling` |
| `--intervals` | Target regions BED file     | `targets.bed`                                  |
| `--wes`       | Whole exome sequencing mode | `true`                                         |

### 3.2. Input samplesheet format

Sarek expects a CSV samplesheet with specific columns.
For starting from FASTQ files:

```csv
patient,sample,lane,fastq_1,fastq_2
patient1,sample1,lane1,sample1_R1.fastq.gz,sample1_R2.fastq.gz
patient1,sample2,lane1,sample2_R1.fastq.gz,sample2_R2.fastq.gz
```

For starting from BAM files (using `--step markduplicates` or later):

```csv
patient,sample,bam,bai
patient1,sample1,sample1.bam,sample1.bam.bai
```

### 3.3. Using GATK joint calling

Remember the joint calling workflow you built in Part 3?
Sarek can do joint calling with GATK HaplotypeCaller:

```bash
nextflow run nf-core/sarek -r 3.5.1 -profile docker \
  --input samplesheet.csv \
  --outdir results \
  --genome GATK.GRCh38 \
  --tools haplotypecaller \
  --joint_germline
```

The `--joint_germline` flag enables the same joint genotyping approach you implemented manually: generating GVCFs per sample, combining them, and running GenotypeGVCFs.

??? info "Other variant callers"

    Sarek supports multiple variant callers beyond GATK HaplotypeCaller.
    Use the `--tools` parameter to specify which ones to run:

    ```bash
    nextflow run nf-core/sarek -r 3.5.1 -profile docker \
      --input samplesheet.csv \
      --outdir results \
      --genome GATK.GRCh38 \
      --tools haplotypecaller,deepvariant,strelka
    ```

    - For germline variant calling (inherited variants), use tools like `haplotypecaller`, `deepvariant`, `freebayes`, or `strelka`
    - For somatic variant calling (tumor mutations), use tools like `mutect2`, `strelka` (with tumor/normal pairs), or `manta`

### Takeaway

Sarek is highly configurable and supports many use cases through its parameters and profiles.

### What's next?

Congratulations on completing the genomics training course!

---

## 4. Summary

In this part, you learned how to run nf-core/sarek, a production-ready variant calling pipeline.

??? info "Troubleshooting common issues"

    **Container download failures**

    If container downloads fail, try running again. Network issues are often transient.
    For persistent issues, consider pre-pulling containers:

    ```bash
    docker pull nfcore/sarek:3.5.1
    ```

    **Out of memory errors**

    Sarek's default resource allocations may not suit all environments.
    Create a custom config to adjust resources:

    ```groovy title="custom.config"
    process {
        withName: 'BWAMEM2_MEM' {
            memory = 32.GB
        }
        withName: 'GATK4_APPLYBQSR' {
            memory = 16.GB
        }
    }
    ```

    Then include it in your run:

    ```bash
    nextflow run nf-core/sarek -r 3.5.1 -profile test,docker \
      -c custom.config \
      --outdir results_sarek
    ```

    **Debugging failed tasks**

    When a task fails, Nextflow shows the work directory hash.
    Navigate there to inspect logs:

    ```bash
    # Example: if hash is [ab/123456]
    ls work/ab/123456*/
    cat work/ab/123456*/.command.err
    cat work/ab/123456*/.command.log
    ```

??? info "Getting help"

    nf-core provides excellent support resources:

    - **Documentation**: [nf-co.re/sarek/docs](https://nf-co.re/sarek/docs)
    - **Parameter reference**: [nf-co.re/sarek/parameters](https://nf-co.re/sarek/parameters)
    - **Slack community**: [nf-co.re/join](https://nf-co.re/join)
    - **GitHub issues**: [github.com/nf-core/sarek/issues](https://github.com/nf-core/sarek/issues)

Consider exploring:

- Other [nf-core pipelines](https://nf-co.re/pipelines) for your research needs
- The [Hello nf-core](../../hello_nf-core/index.md) course for learning to build nf-core-compatible pipelines
- Advanced Nextflow topics in the [Side Quests](../../side_quests/index.md)
