# Part 3: Joint calling on a cohort

In Part 2, you built a per-sample variant calling pipeline that processed each sample's data independently.
Now we're going to extend it to implement joint variant calling, as covered in [Part 1](01_method.md).

## Assignment

In this part of the course, we're going to extend the workflow to do the following:

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/genomics/img/hello-gatk-2.svg"
</figure>

1. Generate an index file for each BAM input file using Samtools
2. Run the GATK HaplotypeCaller on each BAM input file to generate a GVCF of per-sample genomic variant calls
3. Collect all the GVCFs and combine them into a GenomicsDB data store
4. Run joint genotyping on the combined GVCF data store to produce a cohort-level VCF

This part builds directly on the workflow produced by Part 2.

??? info "How to begin from this section"

    This section of the course assumes you have completed [Part 2: Per-sample variant calling](./02_per_sample_variant_calling.md) and have a working `genomics.nf` pipeline.

    If you did not complete Part 2 or want to start fresh for this part, you can use the Part 2 solution as your starting point.
    Run these commands from inside the `nf4-science/genomics/` directory:

    ```bash
    cp solutions/part2/genomics-2.nf genomics.nf
    cp solutions/part2/nextflow.config .
    cp solutions/part2/modules/* modules/
    ```

    This gives you a complete per-sample variant calling workflow.
    You can test that it runs successfully by running the following command:

    ```bash
    nextflow run genomics.nf -profile test
    ```

## Lesson plan

We've broken this down into two steps:

1. **Modify the per-sample variant calling step to produce a GVCF.**
   This covers updating process commands and outputs.
2. **Add a joint genotyping step that combines and genotypes the per-sample GVCFs.**
   This introduces the `collect()` operator, Groovy closures for command-line construction, and multi-command processes.

!!! note

     Make sure you're in the correct working directory:
     `cd /workspaces/training/nf4-science/genomics`

---

## 1. Modify the per-sample variant calling step to produce a GVCF

The pipeline from Part 2 produces VCF files, but joint calling requires GVCF files.
We need to switch on the GVCF variant calling mode and update the output file extension.

Recall the GVCF variant calling command from [Part 1](01_method.md):

```bash
gatk HaplotypeCaller \
        -R /data/ref/ref.fasta \
        -I /data/bam/reads_mother.bam \
        -O reads_mother.g.vcf \
        -L /data/ref/intervals.bed \
        -ERC GVCF
```

Compared to the base HaplotypeCaller command we wrapped in Part 2, the differences are the `-ERC GVCF` parameter and the `.g.vcf` output extension.

### 1.1. Tell HaplotypeCaller to emit a GVCF and update the output extension

Open the `modules/gatk_haplotypecaller.nf` module file to make two changes:

- Add the `-ERC GVCF` parameter to the GATK HaplotypeCaller command;
- Update the output file path to use the corresponding `.g.vcf` extension, as per GATK convention.

Make sure you add a backslash (`\`) at the end of the previous line when you add `-ERC GVCF`.

=== "After"

    ```groovy title="modules/gatk_haplotypecaller.nf" linenums="22" hl_lines="5 7"
        """
        gatk HaplotypeCaller \
            -R ${ref_fasta} \
            -I ${input_bam} \
            -O ${input_bam}.g.vcf \
            -L ${interval_list} \
            -ERC GVCF
        """
    ```

=== "Before"

    ```groovy title="modules/gatk_haplotypecaller.nf" linenums="22" hl_lines="5"
        """
        gatk HaplotypeCaller \
            -R ${ref_fasta} \
            -I ${input_bam} \
            -O ${input_bam}.vcf \
            -L ${interval_list}
        """
    ```

We also need to update the output block to match the new file extension.
Since we changed the command output from `.vcf` to `.g.vcf`, the process `output:` block must reflect the same change.

### 1.2. Update the output file extension in the process outputs block

=== "After"

    ```groovy title="modules/gatk_haplotypecaller.nf" linenums="17" hl_lines="2 3"
        output:
        path "${input_bam}.g.vcf"     , emit: vcf
        path "${input_bam}.g.vcf.idx" , emit: idx
    ```

=== "Before"

    ```groovy title="modules/gatk_haplotypecaller.nf" linenums="17" hl_lines="2 3"
        output:
        path "${input_bam}.vcf"     , emit: vcf
        path "${input_bam}.vcf.idx" , emit: idx
    ```

We also need to update the workflow's publish and output configuration to reflect the new GVCF outputs.

### 1.3. Update the publish targets for the new GVCF outputs

Since we're now producing GVCFs instead of VCFs, we should update the workflow's `publish:` section to use more descriptive names.
We'll also organize the GVCF files into their own subdirectory for clarity.

=== "After"

    ```groovy title="genomics.nf" linenums="47" hl_lines="3 4"
        publish:
        indexed_bam = SAMTOOLS_INDEX.out
        gvcf = GATK_HAPLOTYPECALLER.out.vcf
        gvcf_idx = GATK_HAPLOTYPECALLER.out.idx
    ```

=== "Before"

    ```groovy title="genomics.nf" linenums="47"
        publish:
        indexed_bam = SAMTOOLS_INDEX.out
        vcf = GATK_HAPLOTYPECALLER.out.vcf
        vcf_idx = GATK_HAPLOTYPECALLER.out.idx
    ```

Now update the output block to match.

### 1.4. Update the output block for the new directory structure

We also need to update the `output` block to put the GVCF files in a `gvcf` subdirectory.

=== "After"

    ```groovy title="genomics.nf" linenums="53" hl_lines="3 5 6 8 9"
    output {
        indexed_bam {
            path 'indexed_bam'
        }
        gvcf {
            path 'gvcf'
        }
        gvcf_idx {
            path 'gvcf'
        }
    }
    ```

=== "Before"

    ```groovy title="genomics.nf" linenums="53"
    output {
        indexed_bam {
            path 'bam'
        }
        vcf {
            path 'vcf'
        }
        vcf_idx {
            path 'vcf'
        }
    }
    ```

With the module, publish targets, and output block all updated, we can test the changes.

### 1.5. Run the pipeline

Run the workflow to verify the changes work.

```bash
nextflow run genomics.nf
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics.nf` [nostalgic_franklin] DSL2 - revision: f2c0a93c6a

    executor >  local (6)
    [cc/fbc705] SAMTOOLS_INDEX (3)       | 3 of 3 ✔
    [27/0d7eb9] GATK_HAPLOTYPECALLER (2) | 3 of 3 ✔
    ```

The Nextflow output looks the same as before, but the `.g.vcf` files and their index files are now organized in subdirectories.

??? abstract "Directory contents (symlinks shortened)"

    ```console
    results/
    ├── gvcf/
    │   ├── reads_father.bam.g.vcf -> */27/0d7eb9*/reads_father.bam.g.vcf
    │   ├── reads_father.bam.g.vcf.idx -> */27/0d7eb9*/reads_father.bam.g.vcf.idx
    │   ├── reads_mother.bam.g.vcf -> */e4/4ed55e*/reads_mother.bam.g.vcf
    │   ├── reads_mother.bam.g.vcf.idx -> */e4/4ed55e*/reads_mother.bam.g.vcf.idx
    │   ├── reads_son.bam.g.vcf -> */08/e95962*/reads_son.bam.g.vcf
    │   └── reads_son.bam.g.vcf.idx -> */08/e95962*/reads_son.bam.g.vcf.idx
    └── indexed_bam/
        ├── reads_father.bam -> */9a/c7a873*/reads_father.bam
        ├── reads_father.bam.bai -> */9a/c7a873*/reads_father.bam.bai
        ├── reads_mother.bam -> */f1/8d8486*/reads_mother.bam
        ├── reads_mother.bam.bai -> */f1/8d8486*/reads_mother.bam.bai
        ├── reads_son.bam -> */cc/fbc705*/reads_son.bam
        └── reads_son.bam.bai -> */cc/fbc705*/reads_son.bam.bai
    ```

If you open one of the GVCF files and scroll through it, you can verify that GATK HaplotypeCaller produced GVCF files as requested.

### Takeaway

When you change a tool command's output filename, the process `output:` block and publish/output configuration must be updated to match.

### What's next?

Learn to collect the contents of a channel and pass them on to the next process as a single input.

---

## 2. Add a joint genotyping step

We now need to collect the per-sample GVCFs, combine them into a GenomicsDB data store, and run joint genotyping to produce a cohort-level VCF.
As covered in [Part 1](01_method.md), this is a two-tool operation: GenomicsDBImport combines the GVCFs, then GenotypeGVCFs produces the final variant calls.
We'll wrap both tools in a single process called `GATK_JOINTGENOTYPING`.

Recall the two commands from [Part 1](01_method.md):

```bash
gatk GenomicsDBImport \
    -V reads_mother.g.vcf \
    -V reads_father.g.vcf \
    -V reads_son.g.vcf \
    -L /data/ref/intervals.bed \
    --genomicsdb-workspace-path family_trio_gdb
```

```bash
gatk GenotypeGVCFs \
    -R /data/ref/ref.fasta \
    -V gendb://family_trio_gdb \
    -O family_trio.vcf
```

The first command takes the per-sample GVCFs and an intervals file, and produces a GenomicsDB data store.
The second takes that data store, a reference genome, and produces the final cohort-level VCF.
The container URI is the same as for HaplotypeCaller: `community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867`.

### 2.1. Set up the inputs

The joint genotyping process needs two kinds of inputs that we don't have yet: an arbitrary cohort name, and the collected GVCF outputs from all samples bundled together.

#### 2.1.1. Add a `cohort_name` parameter

We need to provide an arbitrary name for the cohort.
Later in the training series you'll learn how to use sample metadata for this sort of thing, but for now we just declare a CLI parameter using `params` and give it a default value for convenience.

=== "After"

    ```groovy title="genomics.nf" linenums="14" hl_lines="3-4"
        intervals: Path = "${projectDir}/data/ref/intervals.bed"

        // Base name for final output file
        cohort_name: String = "family_trio"
    }
    ```

=== "Before"

    ```groovy title="genomics.nf" linenums="14"
        intervals: Path = "${projectDir}/data/ref/intervals.bed"
    }
    ```

#### 2.1.2. Gather the HaplotypeCaller outputs across samples

If we were to plug the output channel from `GATK_HAPLOTYPECALLER` directly into the new process, Nextflow would call the process on each sample GVCF separately.
We want to bundle all three GVCFs (and their index files) so that Nextflow hands all of them together to a single process call.

We can do that using the `collect()` channel operator.
Add the following lines to the `workflow` body, right after the call to GATK_HAPLOTYPECALLER:

=== "After"

    ```groovy title="genomics.nf" hl_lines="4-6"
            intervals_file
        )

        // Collect variant calling outputs across samples
        all_gvcfs_ch = GATK_HAPLOTYPECALLER.out.vcf.collect()
        all_idxs_ch = GATK_HAPLOTYPECALLER.out.idx.collect()
    ```

=== "Before"

    ```groovy title="genomics.nf"
            intervals_file
        )
    ```

Breaking this down:

1. We take the output channel from `GATK_HAPLOTYPECALLER` using the `.out` property.
2. Because we named the outputs using `emit:` in section 1, we can select the GVCFs with `.vcf` and the index files with `.idx`. Without named outputs, we would have to use `.out[0]` and `.out[1]`.
3. The `collect()` operator bundles all files into a single element, so `all_gvcfs_ch` contains all three GVCFs together, and `all_idxs_ch` contains all three index files together.

We can collect the GVCFs and their index files separately (as opposed to keeping them together in tuples) because Nextflow will stage all input files together for execution, so the index files will be present alongside the GVCFs.

!!! tip

    You can use the `view()` operator to inspect the contents of channels before and after applying channel operators.

### 2.2. Write the joint genotyping process and call it in the workflow

Following the same pattern we used in Part 2, we'll write the process definition in a module file, import it into the workflow, and call it on the inputs we just prepared.

#### 2.2.1. Construct a string to give each GVCF a `-V` argument

Before we start filling in the process definition, there's one thing to work out.
The GenomicsDBImport command expects a separate `-V` argument for each GVCF file, like this:

```bash
gatk GenomicsDBImport \
    -V reads_mother.bam.g.vcf \
    -V reads_father.bam.g.vcf \
    -V reads_son.bam.g.vcf \
    ...
```

If we were to write `-V ${all_gvcfs_ch}`, Nextflow would simply concatenate the filenames and that part of the command would look like this:

```groovy
-V reads_mother.bam.g.vcf reads_father.bam.g.vcf reads_son.bam.g.vcf
```

But we need the string to look like this:

```groovy
-V reads_mother.bam.g.vcf -V reads_father.bam.g.vcf -V reads_son.bam.g.vcf
```

Importantly, we need to construct this string dynamically from whatever files are in the collected channel.
Nextflow (via Groovy) provides a concise way to do this:

```groovy
def gvcfs_line = all_gvcfs.collect { gvcf -> "-V ${gvcf}" }.join(' ')
```

Breaking this down:

1. `all_gvcfs.collect { gvcf -> "-V ${gvcf}" }` iterates over each file path and prepends `-V ` to it, producing `["-V A.g.vcf", "-V B.g.vcf", "-V C.g.vcf"]`.
2. `.join(' ')` concatenates them with spaces: `"-V A.g.vcf -V B.g.vcf -V C.g.vcf"`.
3. The result is assigned to a local variable `gvcfs_line` (defined with `def`), which we can interpolate into the command template.

This line goes inside the process's `script:` block, before the command template.
You can place arbitrary Groovy code between `script:` and the opening `"""` of the command template.

Then you'll be able to refer to that whole string as `gvcfs_line` in the `script:` block of the process.

#### 2.2.2. Fill in the module for the joint genotyping process

Now we can tackle writing the full process.

Open `modules/gatk_jointgenotyping.nf` and examine the outline of the process definition.

Go ahead and fill in the process definition using the information provided above, then check your work against the solution in the "After" tab below.

=== "Before"

    ```groovy title="modules/gatk_jointgenotyping.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
     * Combine GVCFs into GenomicsDB datastore and run joint genotyping to produce cohort-level calls
     */
    process GATK_JOINTGENOTYPING {

        container

        input:

        output:

        script:
        """

        """
    }
    ```

=== "After"

    ```groovy title="modules/gatk_jointgenotyping.nf" linenums="1" hl_lines="8 11-17 20-21 24-25 29-33"
    #!/usr/bin/env nextflow

    /*
     * Combine GVCFs into GenomicsDB datastore and run joint genotyping to produce cohort-level calls
     */
    process GATK_JOINTGENOTYPING {

        container "community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867"

        input:
        path all_gvcfs
        path all_idxs
        path interval_list
        val cohort_name
        path ref_fasta
        path ref_index
        path ref_dict

        output:
        path "${cohort_name}.joint.vcf"     , emit: vcf
        path "${cohort_name}.joint.vcf.idx" , emit: idx

        script:
        def gvcfs_line = all_gvcfs.collect { gvcf -> "-V ${gvcf}" }.join(' ')
        """
        gatk GenomicsDBImport \
            ${gvcfs_line} \
            -L ${interval_list} \
            --genomicsdb-workspace-path ${cohort_name}_gdb

        gatk GenotypeGVCFs \
            -R ${ref_fasta} \
            -V gendb://${cohort_name}_gdb \
            -L ${interval_list} \
            -O ${cohort_name}.joint.vcf
        """
    }
    ```

There are several things worth calling out here.

As previously, several inputs are listed even though the commands don't reference them directly: `all_idxs`, `ref_index`, and `ref_dict`.
Listing them ensures Nextflow stages these files in the working directory alongside the files that do appear in the commands, which GATK expects to find based on naming conventions.

The `gvcfs_line` variable uses the Groovy closure described above to construct the `-V` arguments for GenomicsDBImport.

This process runs two commands in serial, just as you would in the terminal.
GenomicsDBImport combines the per-sample GVCFs into a data store, then GenotypeGVCFs reads that data store and produces the final cohort-level VCF.
The GenomicsDB data store (`${cohort_name}_gdb`) is an intermediate artifact used only within the process; it doesn't appear in the output block.

Once you've completed this, the process is ready to use.
To use it in the workflow, you'll need to import the module and add a process call.

#### 2.2.3. Import the module

Add the import statement to `genomics.nf`, below the existing import statements:

=== "After"

    ```groovy title="genomics.nf" linenums="21" hl_lines="3"
    include { SAMTOOLS_INDEX } from './modules/samtools_index.nf'
    include { GATK_HAPLOTYPECALLER } from './modules/gatk_haplotypecaller.nf'
    include { GATK_JOINTGENOTYPING } from './modules/gatk_jointgenotyping.nf'
    ```

=== "Before"

    ```groovy title="genomics.nf" linenums="21"
    include { SAMTOOLS_INDEX } from './modules/samtools_index.nf'
    include { GATK_HAPLOTYPECALLER } from './modules/gatk_haplotypecaller.nf'
    ```

The process is now available in the workflow scope.

#### 2.2.4. Add the process call

Add the call to `GATK_JOINTGENOTYPING` in the workflow body, after the `collect()` lines:

=== "After"

    ```groovy title="genomics.nf" hl_lines="3-12"
        all_idxs_ch = GATK_HAPLOTYPECALLER.out.idx.collect()

        // Combine GVCFs into a GenomicsDB data store and apply joint genotyping
        GATK_JOINTGENOTYPING(
            all_gvcfs_ch,
            all_idxs_ch,
            intervals_file,
            params.cohort_name,
            ref_file,
            ref_index_file,
            ref_dict_file
        )
    ```

=== "Before"

    ```groovy title="genomics.nf"
        all_idxs_ch = GATK_HAPLOTYPECALLER.out.idx.collect()
    ```

The process is now fully wired up.
Next, we configure how the outputs are published.

### 2.3. Configure the output handling

We need to publish the joint VCF outputs.
Add publish targets and output block entries for the joint genotyping results.

#### 2.3.1. Add publish targets for the joint VCF

Add the joint VCF and its index to the `publish:` section of the workflow:

=== "After"

    ```groovy title="genomics.nf" hl_lines="5-6"
        publish:
        indexed_bam = SAMTOOLS_INDEX.out
        gvcf = GATK_HAPLOTYPECALLER.out.vcf
        gvcf_idx = GATK_HAPLOTYPECALLER.out.idx
        joint_vcf = GATK_JOINTGENOTYPING.out.vcf
        joint_vcf_idx = GATK_JOINTGENOTYPING.out.idx
    ```

=== "Before"

    ```groovy title="genomics.nf"
        publish:
        indexed_bam = SAMTOOLS_INDEX.out
        gvcf = GATK_HAPLOTYPECALLER.out.vcf
        gvcf_idx = GATK_HAPLOTYPECALLER.out.idx
    ```

Now update the output block to match.

#### 2.3.2. Add output block entries for the joint VCF

Add entries for the joint VCF files.
We'll put them at the root of the results directory since this is the final output.

=== "After"

    ```groovy title="genomics.nf" hl_lines="11-16"
    output {
        indexed_bam {
            path 'indexed_bam'
        }
        gvcf {
            path 'gvcf'
        }
        gvcf_idx {
            path 'gvcf'
        }
        joint_vcf {
            path '.'
        }
        joint_vcf_idx {
            path '.'
        }
    }
    ```

=== "Before"

    ```groovy title="genomics.nf"
    output {
        indexed_bam {
            path 'indexed_bam'
        }
        gvcf {
            path 'gvcf'
        }
        gvcf_idx {
            path 'gvcf'
        }
    }
    ```

With the process, publish targets, and output block all in place, we can test the complete workflow.

### 2.4. Run the workflow

Run the workflow to verify everything works.

```bash
nextflow run genomics.nf -resume
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics.nf` [crazy_marconi] DSL2 - revision: 5da9afc841

    executor >  local (1)
    [9a/c7a873] SAMTOOLS_INDEX (2)       | 3 of 3, cached: 3 ✔
    [e4/4ed55e] GATK_HAPLOTYPECALLER (2) | 3 of 3, cached: 3 ✔
    [a6/7cc8ed] GATK_JOINTGENOTYPING     | 1 of 1 ✔
    ```

The first two steps are cached from the previous run, and the new `GATK_JOINTGENOTYPING` step runs once on the collected inputs from all three samples.
The final output file, `family_trio.joint.vcf` (and its index), are in the results directory.

??? abstract "Directory contents (symlinks shortened)"

    ```console
    results/
    ├── family_trio.joint.vcf -> */a6/7cc8ed*/family_trio.joint.vcf
    ├── family_trio.joint.vcf.idx -> */a6/7cc8ed*/family_trio.joint.vcf.idx
    ├── gvcf/
    │   ├── reads_father.bam.g.vcf -> */27/0d7eb9*/reads_father.bam.g.vcf
    │   ├── reads_father.bam.g.vcf.idx -> */27/0d7eb9*/reads_father.bam.g.vcf.idx
    │   ├── reads_mother.bam.g.vcf -> */e4/4ed55e*/reads_mother.bam.g.vcf
    │   ├── reads_mother.bam.g.vcf.idx -> */e4/4ed55e*/reads_mother.bam.g.vcf.idx
    │   ├── reads_son.bam.g.vcf -> */08/e95962*/reads_son.bam.g.vcf
    │   └── reads_son.bam.g.vcf.idx -> */08/e95962*/reads_son.bam.g.vcf.idx
    └── indexed_bam/
        ├── reads_father.bam -> */9a/c7a873*/reads_father.bam
        ├── reads_father.bam.bai -> */9a/c7a873*/reads_father.bam.bai
        ├── reads_mother.bam -> */f1/8d8486*/reads_mother.bam
        ├── reads_mother.bam.bai -> */f1/8d8486*/reads_mother.bam.bai
        ├── reads_son.bam -> */cc/fbc705*/reads_son.bam
        └── reads_son.bam.bai -> */cc/fbc705*/reads_son.bam.bai
    ```

If you open the joint VCF file, you can verify that the workflow produced the expected variant calls.

```console title="family_trio.joint.vcf" linenums="40"
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	reads_father	reads_mother	reads_son
20_10037292_10066351	3480	.	C	CT	1625.89	.	AC=5;AF=0.833;AN=6;BaseQRankSum=0.220;DP=85;ExcessHet=0.0000;FS=2.476;MLEAC=5;MLEAF=0.833;MQ=60.00;MQRankSum=0.00;QD=21.68;ReadPosRankSum=-1.147e+00;SOR=0.487	GT:AD:DP:GQ:PL	0/1:15,16:31:99:367,0,375	1/1:0,18:18:54:517,54,0	1/1:0,26:26:78:756,78,0
20_10037292_10066351	3520	.	AT	A	1678.89	.	AC=5;AF=0.833;AN=6;BaseQRankSum=1.03;DP=80;ExcessHet=0.0000;FS=2.290;MLEAC=5;MLEAF=0.833;MQ=60.00;MQRankSum=0.00;QD=22.39;ReadPosRankSum=0.701;SOR=0.730	GT:AD:DP:GQ:PL	0/1:18,13:31:99:296,0,424	1/1:0,18:18:54:623,54,0	1/1:0,26:26:78:774,78,0
20_10037292_10066351	3529	.	T	A	154.29	.	AC=1;AF=0.167;AN=6;BaseQRankSum=-5.440e-01;DP=104;ExcessHet=0.0000;FS=1.871;MLEAC=1;MLEAF=0.167;MQ=60.00;MQRankSum=0.00;QD=7.71;ReadPosRankSum=-1.158e+00;SOR=1.034	GT:AD:DP:GQ:PL	0/0:44,0:44:99:0,112,1347	0/1:12,8:20:99:163,0,328	0/0:39,0:39:99:0,105,1194
```

You now have an automated, fully reproducible joint variant calling workflow!

!!! note

    Keep in mind the data files we gave you cover only a tiny portion of chromosome 20.
    The real size of a variant callset would be counted in millions of variants.
    That's why we use only tiny subsets of data for training purposes!

### Takeaway

You know how to collect outputs from a channel and bundle them as a single input to another process.
You also know how to construct a command line using Groovy closures, and how to run multiple commands in a single process.

### What's next?

Give yourself a big pat on the back! You have completed the Nextflow for Genomics course.

Head on to the final [course summary](./next_steps.md) to review what you learned and find out what comes next.
