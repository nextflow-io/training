# Part 2: Per-sample variant calling

In Part 1, you tested the Samtools and GATK commands manually in their respective containers.
Now we're going to wrap those same commands into a Nextflow workflow.

## Assignment

In this part of the course, we're going to develop a workflow that does the following:

1. Generate an index file for each BAM input file using [Samtools](https://www.htslib.org/)
2. Run the GATK HaplotypeCaller on each BAM input file to generate per-sample variant calls in VCF (Variant Call Format)

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/genomics/img/hello-gatk-1.svg"
</figure>

[TODO: sentence about how this replicates what we did in section X of Part 1.]

As a starting point, we provide you with a workflow file, `genomics-1.nf`, that outlines the main parts of the workflow, as well as two module files, samtools_index.nf and gatk_haplotypecaller.nf, that outline the structure of the modules.
These files are not functional; their purpose is just to serve as scaffolds for you to fill in with the interesting parts of the code.

## Lesson plan

In order to make the development process more educational, we've broken this down into four distinct steps:

1. **Write a single-stage workflow that runs Samtools index on a BAM file.**
   [TODO: single-sentence summary]
2. **Add a second process to run GATK HaplotypeCaller on the indexed BAM file.**
   [TODO: single-sentence summary]
3. **Adapt the workflow to run on a batch of samples.**
   [TODO: single-sentence summary]
4. **Make the workflow accept a text file containing a batch of input files.**
   [TODO: single-sentence summary]
5. **Adapt the input handling to use an nf-core style samplesheet.**
   [TODO: single-sentence summary]

[TODO: sentence about how this will allow us to focus on a specific aspect of workflow development.]

---

## 1. Write a single-stage workflow that runs Samtools index on a BAM file

[TODO: sentence about what this focuses on, ie getting the basics in place: loading a BAM file and doing something to it.]

### 1.1. Set up the inputs

#### 1.1.1. Add an input parameter declaration

In the main workflow file `genomics-1.nf`, under the `Pipeline parameters` section, declare a CLI parameter called `reads_bam`.

[TODO: change `reads_bam` to `input`]

=== "After"

    ```groovy title="genomics-1.nf" linenums="5" hl_lines="4-7"
    /*
     * Pipeline parameters
     */
    params {
        // Primary input
        reads_bam: Path
    }
    ```

=== "Before"

    ```groovy title="genomics-1.nf" linenums="5"
    /*
     * Pipeline parameters
     */

    // Primary input
    ```

[TODO: sentence about how we'll set up a default value later in this section]

#### 1.1.2. Create a test profile with default values in `nextflow.config`

[TODO: sentence about purpose of test profile, with link to Part 6 of Hello Nextflow]

[TODO: add instructions for creating a test profile with a value for the reads_bam input parameter, pointing to the 'mother' bam file and using ${projectDir} for the path]

!!! note

    `${projectDir}` is a built-in Nextflow variable that points to the directory where the current Nextflow workflow script (`genomics-1.nf`) is located.

    This makes it easy to reference files, data directories, and other resources included in the workflow repository without hardcoding absolute paths.

#### 1.1.2. Set up the input channel

[TODO: instructions for setting up the input channel]

We're using the same `.fromPath` channel factory as in [Hello Channels](../../hello_nextflow/02_hello_channels.md).
The difference is that we're telling Nextflow to just load the file path itself into the channel as an input element, rather than reading in its contents.

[TODO: transition sentence saying now we need a process to run indexing on the input]

### 1.2. Set up the indexing step

#### 1.2.1. Fill in the module for the indexing process

Open `modules/samtools_index.nf` and examine the outline of the process definition.
[TODO: sentence about how you should recognize all the main elements if you worked through the beginner training. If not, consider reading through that for a refresher before continuing.]

```groovy title="modules/samtools_index.nf" linenums="1"
#!/usr/bin/env nextflow

/*
 * Generate BAM index file
 */
process SAMTOOLS_INDEX {

    container

    input:

    output:

    script:
    """

    """
}
```

[TODO: instructions to recall the commands we ran in Part 1 and fill out the process based on that]

[TODO: recap the command we used in Part 1]

[TODO: sentence about how you can try to fill out the process definition by yourself if you're feeling confident. If not, take a peek at the solution below]

[TODO: put the before/after below inside an admonition with an appropriate type]

=== "After"

    ```groovy title="modules/samtools_index.nf" linenums="1" hl_lines="8 11 14 18"
    #!/usr/bin/env nextflow

    /*
     * Generate BAM index file
     */
    process SAMTOOLS_INDEX {

        container 'community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464'

        input:
        path input_bam

        output:
        path "${input_bam}.bai"

        script:
        """
        samtools index '$input_bam'
        """
    }
    ```

=== "Before"

    ```groovy title="modules/samtools_index.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
     * Generate BAM index file
     */
    process SAMTOOLS_INDEX {

        container

        input:

        output:

        script:
        """

        """
    }
    ```

[TODO: sentence about how the process is now ready, so we need to make it available in the workflow]

#### 1.2.2. Include the module

In `genomics-1.nf`, import the module:

=== "After"

    ```groovy title="genomics-1.nf" linenums="3" hl_lines="2"
    // Module INCLUDE statements
    include { SAMTOOLS_INDEX } from './modules/samtools_index.nf'
    ```

=== "Before"

    ```groovy title="genomics-1.nf" linenums="3"
    // Module INCLUDE statements
    ```

[TODO: sentence about how the process is now available, we can call it]

#### 1.2.3. Call the indexing process on the input

[TODO: instructions to add a call in the workflow block under `main:`]

=== "After"

    ```groovy title="genomics-1.nf" linenums="14" hl_lines="4-5 7-8 11 15-17"
    workflow {

        main:
        // Create input channel (single file via CLI parameter)
        reads_ch = channel.fromPath(params.reads_bam)

        // Create index file for input BAM file
        SAMTOOLS_INDEX(reads_ch)

        publish:
        bam_index = SAMTOOLS_INDEX.out
    }

    output {
        bam_index {
            path 'bam'
        }
    }
    ```

=== "Before"

    ```groovy title="genomics-1.nf" linenums="14"
    workflow {

        main:
        // Create input channel

        // Call processes

        publish:
        // Declare outputs to publish
    }

    output {
        // Configure publish targets
    }
    ```



    By default, Nextflow publishes output files as symbolic links, which avoids unnecessary duplication.
    Even though the data files we're using here are very small, in genomics they can get very large.
    Be aware that symlinks will break when you clean up your `work` directory, so for production workflows you may want to override the default publish mode to `'copy'`.

### 1.3. Set up the output handling

#### 1.3.1. Add an output declaration under `publish:`

[TODO: add instructions to add an output declaration in the workflow block]

#### 1.3.2. Specify how the output should be handled in the `output {}` block

[TODO: add instructions to add the output handling in the output block]

### 1.4. Run the workflow to verify that the indexing step works

As a reminder, we don't need to specify an input in the command line because we set up a default value for the input when we declared the input parameter.

```bash
nextflow run genomics-1.nf
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics-1.nf` [reverent_sinoussi] DSL2 - revision: 41d43ad7fe

    executor >  local (1)
    [2a/e69536] SAMTOOLS_INDEX (1) | 1 of 1 ✔
    ```

You can check that the index file has been generated correctly by looking in the work directory or in the results directory.

??? abstract "Work directory contents"

    ```console
    work/2a/e695367b2f60df09cf826b07192dc3
    ├── reads_mother.bam -> /workspaces/training/nf4-science/genomics/data/bam/reads_mother.bam
    └── reads_mother.bam.bai
    ```

??? abstract "Results directory contents"

    ```console
    results/
    └── bam/
        └── reads_mother.bam.bai -> ...
    ```

There it is!

### Takeaway

You know how to create a module containing a process and use it in a Nextflow workflow.
[TODO: expand this with a bit more detail]

### What's next?

Add a second step that consumes the output of the first.
[TODO: expand this with a bit more detail]

---

## 2. Add a second process to run GATK HaplotypeCaller on the indexed BAM file

Now that we have an index for our input file, we can move on to setting up the variant calling step, which is the interesting part of the workflow.

[TODO: adapt this section to match the progression used in the previous section]

### 2.1. Fill in the module for the variant calling process

Open the skeleton file `modules/gatk_haplotypecaller.nf` and fill in the `GATK_HAPLOTYPECALLER` process:

=== "After"

    ```groovy title="modules/gatk_haplotypecaller.nf" linenums="1" hl_lines="8 11-16 19-20 24-28"
    #!/usr/bin/env nextflow

    /*
     * Call variants with GATK HaplotypeCaller
     */
    process GATK_HAPLOTYPECALLER {

        container "community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867"

        input:
        path input_bam
        path input_bam_index
        path ref_fasta
        path ref_index
        path ref_dict
        path interval_list

        output:
        path "${input_bam}.vcf"     , emit: vcf
        path "${input_bam}.vcf.idx" , emit: idx

        script:
        """
        gatk HaplotypeCaller \
            -R ${ref_fasta} \
            -I ${input_bam} \
            -O ${input_bam}.vcf \
            -L ${interval_list}
        """
    }
    ```

=== "Before"

    ```groovy title="modules/gatk_haplotypecaller.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
     * Call variants with GATK HaplotypeCaller
     */
    process GATK_HAPLOTYPECALLER {

        container

        input:

        output:

        script:
        """

        """
    }
    ```

You'll notice that we've introduced some new syntax here (`emit:`) to uniquely name each of our output channels, and the reasons for this will become clear soon.

This command takes quite a few more inputs, because GATK needs more information to perform the analysis compared to a simple indexing job.
But you'll note that there are even more inputs defined in the inputs block than are listed in the GATK command. Why is that?

!!! note

    The GATK knows to look for the BAM index file and the reference genome's accessory files because it is aware of the conventions surrounding those files.
    However, Nextflow is designed to be domain-agnostic and doesn't know anything about bioinformatics file format requirements.

We need to tell Nextflow explicitly that it has to stage those files in the working directory at runtime; otherwise it won't do it, and GATK will (correctly) throw an error about the index files being missing.

Similarly, we have to list the output VCF's index file (the `"${input_bam}.vcf.idx"` file) explicitly so that Nextflow will know to keep track of that file in case it's needed in subsequent steps.

### 2.2. Add definitions for accessory inputs

Since our new process expects a handful of additional files to be provided, add some CLI parameters for them in `genomics-1.nf` under the `Pipeline parameters` section, along with some default values:

=== "After"

    ```groovy title="genomics-1.nf" linenums="9" hl_lines="5-9"
    params {
        // Primary input
        reads_bam: Path = "${projectDir}/data/bam/reads_mother.bam"

        // Accessory files
        reference: Path = "${projectDir}/data/ref/ref.fasta"
        reference_index: Path = "${projectDir}/data/ref/ref.fasta.fai"
        reference_dict: Path = "${projectDir}/data/ref/ref.dict"
        intervals: Path = "${projectDir}/data/ref/intervals.bed"
    }
    ```

=== "Before"

    ```groovy title="genomics-1.nf" linenums="9"
    params {
        // Primary input
        reads_bam: Path = "${projectDir}/data/bam/reads_mother.bam"
    }
    ```

### 2.3. Import the new module and create variables for accessory files

Update `genomics-1.nf` to import the new module:

=== "After"

    ```groovy title="genomics-1.nf" linenums="3" hl_lines="3"
    // Module INCLUDE statements
    include { SAMTOOLS_INDEX } from './modules/samtools_index.nf'
    include { GATK_HAPLOTYPECALLER } from './modules/gatk_haplotypecaller.nf'
    ```

=== "Before"

    ```groovy title="genomics-1.nf" linenums="3"
    // Module INCLUDE statements
    include { SAMTOOLS_INDEX } from './modules/samtools_index.nf'
    ```

Then add variables for the accessory file paths inside the workflow block:

=== "After"

    ```groovy title="genomics-1.nf" linenums="21" hl_lines="7-11"
    workflow {

        main:
        // Create input channel (single file via CLI parameter)
        reads_ch = channel.fromPath(params.reads_bam)

        // Load the file paths for the accessory files (reference and intervals)
        ref_file        = file(params.reference)
        ref_index_file  = file(params.reference_index)
        ref_dict_file   = file(params.reference_dict)
        intervals_file  = file(params.intervals)

        // Create index file for input BAM file
        SAMTOOLS_INDEX(reads_ch)
    ```

=== "Before"

    ```groovy title="genomics-1.nf" linenums="21"
    workflow {

        main:
        // Create input channel (single file via CLI parameter)
        reads_ch = channel.fromPath(params.reads_bam)

        // Create index file for input BAM file
        SAMTOOLS_INDEX(reads_ch)
    ```

### 2.4. Add a call to run GATK_HAPLOTYPECALLER

Now add the process call in the workflow body, under `main:`:

=== "After"

    ```groovy title="genomics-1.nf" linenums="33" hl_lines="4-12"
        // Create index file for input BAM file
        SAMTOOLS_INDEX(reads_ch)

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

=== "Before"

    ```groovy title="genomics-1.nf" linenums="33"
        // Create index file for input BAM file
        SAMTOOLS_INDEX(reads_ch)
    ```

You should recognize the `*.out` syntax from the Hello Nextflow training series; we are telling Nextflow to take the channel output by `SAMTOOLS_INDEX` and plugging that into the `GATK_HAPLOTYPECALLER` process call.

!!! note

    You'll notice that the inputs are provided in the exact same order in the call to the process as they are listed in the input block of the process.
    In Nextflow, inputs are positional, meaning you _must_ follow the same order; and of course there have to be the same number of elements.

### 2.5. Update the publish and output blocks

Add the GATK HaplotypeCaller outputs to the `publish:` section and the `output {}` block:

=== "After"

    ```groovy title="genomics-1.nf" linenums="45" hl_lines="3-4 11-16"
        publish:
        bam_index = SAMTOOLS_INDEX.out
        vcf = GATK_HAPLOTYPECALLER.out.vcf
        vcf_idx = GATK_HAPLOTYPECALLER.out.idx
    }

    output {
        bam_index {
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

=== "Before"

    ```groovy title="genomics-1.nf" linenums="45"
        publish:
        bam_index = SAMTOOLS_INDEX.out
    }

    output {
        bam_index {
            path 'bam'
        }
    }
    ```

The VCF and its index are published as separate targets that both go into the `vcf/` subdirectory.

### 2.6. Run the workflow to verify that the variant calling step works

Run the expanded workflow with `-resume` so that we don't have to run the indexing step again.

```bash
nextflow run genomics-1.nf -resume
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics-1.nf` [grave_volta] DSL2 - revision: 4790abc96a

    executor >  local (1)
    [2a/e69536] SAMTOOLS_INDEX (1)       | 1 of 1, cached: 1 ✔
    [53/e18e98] GATK_HAPLOTYPECALLER (1) | 1 of 1 ✔
    ```

Now if we look at the console output, we see the two processes listed.

The first process was skipped thanks to the caching, as expected, whereas the second process was run since it's brand new.

You'll find the output files in the results directory (as symbolic links to the work directory).

??? abstract "Directory contents"

    ```console
    results/
    ├── bam/
    │   └── reads_mother.bam.bai -> ...
    └── vcf/
        ├── reads_mother.bam.vcf -> ...
        └── reads_mother.bam.vcf.idx -> ...
    ```

If you open the VCF file, you should see the same contents as in the file you generated by running the GATK command directly in the container.

```console title="reads_mother.bam.vcf" linenums="26"
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	reads_mother
20_10037292_10066351	3480	.	C	CT	503.03	.	AC=2;AF=1.00;AN=2;DP=23;ExcessHet=0.0000;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=27.95;SOR=1.179	GT:AD:DP:GQ:PL	1/1:0,18:18:54:517,54,0
20_10037292_10066351	3520	.	AT	A	609.03	.	AC=2;AF=1.00;AN=2;DP=18;ExcessHet=0.0000;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=33.83;SOR=0.693	GT:AD:DP:GQ:PL	1/1:0,18:18:54:623,54,0
20_10037292_10066351	3529	.	T	A	155.64	.	AC=1;AF=0.500;AN=2;BaseQRankSum=-0.544;DP=21;ExcessHet=0.0000;FS=1.871;MLEAC=1;MLEAF=0.500;MQ=60.00;MQRankSum=0.000;QD=7.78;ReadPosRankSum=-1.158;SOR=1.034	GT:AD:DP:GQ:PL	0/1:12,8:20:99:163,0,328
```

This is the output we care about generating for each sample in our study.

### Takeaway

You know how to make a two-step modular workflow that does real analysis work and is capable of dealing with genomics file format idiosyncrasies like the accessory files.

### What's next?

Make the workflow handle multiple samples in bulk.

---

## 3. Adapt the workflow to run on a batch of samples

It's all well and good to have a workflow that can automate processing on a single sample, but what if you have 1000 samples?
Do you need to write a bash script that loops through all your samples?

No, thank goodness! Just make a minor tweak to the code and Nextflow will handle that for you too.

### 3.1. Turn the input parameter declaration into an array listing the three samples

Turn that default file path in the input BAM file declaration into an array listing file paths for our three test samples, up under the `Pipeline parameters` section.

=== "After"

    ```groovy title="genomics-1.nf" linenums="10" hl_lines="1-6"
        // Primary input (array of three samples)
        reads_bam = [
            "${projectDir}/data/bam/reads_mother.bam",
            "${projectDir}/data/bam/reads_father.bam",
            "${projectDir}/data/bam/reads_son.bam"
        ]
    ```

=== "Before"

    ```groovy title="genomics-1.nf" linenums="10"
        // Primary input
        reads_bam: Path = "${projectDir}/data/bam/reads_mother.bam"
    ```

!!! note

    When using typed parameter declarations (like `reads_bam: Path`), you cannot assign an array value.
    For arrays, omit the type annotation.

And that's actually all we need to do, because the channel factory we use in the workflow body (`.fromPath`) is just as happy to accept multiple file paths to load into the input channel as it was to load a single one.

!!! note

    Normally, you wouldn't want to hardcode the list of samples into your workflow file, but we're doing that here to keep things simple.
    We'll present more elegant ways for handling inputs later in this training series.

### 3.2. Run the workflow to verify that it runs on all three samples

Try running the workflow now that the plumbing is set up to run on all three test samples.

```bash
nextflow run genomics-1.nf -resume
```

Funny thing: this _might work_, OR it _might fail_. For example, here's a run that succeeded:

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics-1.nf` [peaceful_yalow] DSL2 - revision: a256d113ad

    executor >  local (6)
    [4f/7071b0] SAMTOOLS_INDEX (3)       | 3 of 3, cached: 1 ✔
    [7a/89bc43] GATK_HAPLOTYPECALLER (2) | 3 of 3, cached: 1 ✔
    ```

If your workflow run succeeded, run it again until you get an error like this:

??? failure "Command output"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics-1.nf` [loving_pasteur] DSL2 - revision: d2a8e63076

    executor >  local (4)
    [01/eea165] SAMTOOLS_INDEX (2)       | 3 of 3, cached: 1 ✔
    [a5/fa9fd0] GATK_HAPLOTYPECALLER (3) | 1 of 3, cached: 1
    ERROR ~ Error executing process > 'GATK_HAPLOTYPECALLER (2)'

    Caused by:
      Process `GATK_HAPLOTYPECALLER (2)` terminated with an error exit status (2)

    Command executed:

      gatk HaplotypeCaller         -R ref.fasta         -I reads_father.bam         -O reads_father.bam.vcf         -L intervals.bed

    Command exit status:
      2

    Command error:
      ...
      A USER ERROR has occurred: Traversal by intervals was requested but some input files are not indexed.
      ...
    ```

If you look at the GATK command error output, there will be a line like this:

```console
A USER ERROR has occurred: Traversal by intervals was requested but some input files are not indexed.
```

Well, that's weird, considering we explicitly indexed the BAM files in the first step of the workflow. Could there be something wrong with the plumbing?

#### 3.2.1. Check the work directories for the relevant calls

Take a look inside the work directory for the failed `GATK_HAPLOTYPECALLER` process call listed in the console output.

??? abstract "Directory contents"

    ```console
    work/a5/fa9fd0994b6beede5fb9ea073596c2
    ├── intervals.bed -> /workspaces/training/nf4-science/genomics/data/ref/intervals.bed
    ├── reads_father.bam.bai -> /workspaces/training/nf4-science/genomics/work/01/eea16597bd6e810fb4cf89e60f8c2d/reads_father.bam.bai
    ├── reads_son.bam -> /workspaces/training/nf4-science/genomics/data/bam/reads_son.bam
    ├── reads_son.bam.vcf
    ├── reads_son.bam.vcf.idx
    ├── ref.dict -> /workspaces/training/nf4-science/genomics/data/ref/ref.dict
    ├── ref.fasta -> /workspaces/training/nf4-science/genomics/data/ref/ref.fasta
    └── ref.fasta.fai -> /workspaces/training/nf4-science/genomics/data/ref/ref.fasta.fai
    ```

Pay particular attention to the names of the BAM file and the BAM index that are listed in this directory: `reads_son.bam` and `reads_father.bam.bai`.

What the heck? Nextflow has staged an index file in this process call's work directory, but it's the wrong one. How could this have happened?

#### 3.2.2. Use the [view() operator](https://www.nextflow.io/docs/latest/reference/operator.html#view) to inspect channel contents

Add these two lines in the workflow body before the `GATK_HAPLOTYPECALLER` process call:

=== "After"

    ```groovy title="genomics-1.nf" hl_lines="3-5"
        SAMTOOLS_INDEX(reads_ch)

        // temporary diagnostics
        reads_ch.view()
        SAMTOOLS_INDEX.out.view()

        // Call variants from the indexed BAM file
        GATK_HAPLOTYPECALLER(
    ```

=== "Before"

    ```groovy title="genomics-1.nf"
        SAMTOOLS_INDEX(reads_ch)

        // Call variants from the indexed BAM file
        GATK_HAPLOTYPECALLER(
    ```

Then run the workflow command again.

```bash
nextflow run genomics-1.nf
```

Once again, this may succeed or fail. Here's what the output of the two `.view()` calls looks like for a failed run:

```console
/workspaces/training/nf4-science/genomics/data/bam/reads_mother.bam
/workspaces/training/nf4-science/genomics/data/bam/reads_father.bam
/workspaces/training/nf4-science/genomics/data/bam/reads_son.bam
/workspaces/training/nf4-science/genomics/work/9c/53492e3518447b75363e1cd951be4b/reads_father.bam.bai
/workspaces/training/nf4-science/genomics/work/cc/37894fffdf6cc84c3b0b47f9b536b7/reads_son.bam.bai
/workspaces/training/nf4-science/genomics/work/4d/dff681a3d137ba7d9866e3d9307bd0/reads_mother.bam.bai
```

The first three lines correspond to the input channel and the second, to the output channel.
You can see that the BAM files and index files for the three samples are not listed in the same order!

!!! note

    When you call a Nextflow process on a channel containing multiple elements, Nextflow will try to parallelize execution as much as possible, and will collect outputs in whatever order they become available.
    The consequence is that the corresponding outputs may be collected in a different order than the original inputs were fed in.

As currently written, our workflow script assumes that the index files will come out of the indexing step listed in the same mother/father/son order as the inputs were given.
But that is not guaranteed to be the case, which is why sometimes (though not always) the wrong files get paired up in the second step.

To fix this, we need to make sure the BAM files and their index files travel together through the channels.

!!! tip

    The `view()` statements in the workflow code don't do anything, so it's not a problem to leave them in.
    However they will clutter up your console output, so we recommend removing them when you're done troubleshooting the issue.

### 3.3. Change the output of the SAMTOOLS_INDEX module to a tuple

The simplest way to ensure a BAM file and its index stay closely associated is to package them together into a tuple coming out of the index task.

!!! note

    A **tuple** is a finite, ordered list of elements that is commonly used for returning multiple values from a function. Tuples are particularly useful for passing multiple inputs or outputs between processes while preserving their association and order.

Update the output in `modules/samtools_index.nf` to include the BAM file:

=== "After"

    ```groovy title="modules/samtools_index.nf" linenums="14" hl_lines="2"
        output:
        tuple path(input_bam), path("${input_bam}.bai")
    ```

=== "Before"

    ```groovy title="modules/samtools_index.nf" linenums="14"
        output:
        path "${input_bam}.bai"
    ```

This way, each index file will be tightly coupled with its original BAM file, and the overall output of the indexing step will be a single channel containing pairs of files.

### 3.4. Change the input of the GATK_HAPLOTYPECALLER module to accept a tuple

Since we've changed the 'shape' of the output of the first process, we need to update the input definition of the second process to match.

Update `modules/gatk_haplotypecaller.nf`:

=== "After"

    ```groovy title="modules/gatk_haplotypecaller.nf" linenums="11" hl_lines="2"
        input:
        tuple path(input_bam), path(input_bam_index)
    ```

=== "Before"

    ```groovy title="modules/gatk_haplotypecaller.nf" linenums="11"
        input:
        path input_bam
        path input_bam_index
    ```

### 3.5. Update the call to GATK_HAPLOTYPECALLER in the workflow

We no longer need to provide the original `reads_ch` to the `GATK_HAPLOTYPECALLER` process, since the BAM file is now bundled into the channel output by `SAMTOOLS_INDEX`.

Update the call in `genomics-1.nf`:

=== "After"

    ```groovy title="genomics-1.nf" linenums="42" hl_lines="2"
        GATK_HAPLOTYPECALLER(
            SAMTOOLS_INDEX.out,
    ```

=== "Before"

    ```groovy title="genomics-1.nf" linenums="42"
        GATK_HAPLOTYPECALLER(
            reads_ch,
            SAMTOOLS_INDEX.out,
    ```

### 3.6. Update the publish target for the indexed BAM output

Since the SAMTOOLS_INDEX output is now a tuple containing both the BAM file and its index, rename the publish target from `bam_index` to `indexed_bam` to better reflect its contents:

=== "After"

    ```groovy title="genomics-1.nf" linenums="46" hl_lines="2 8"
        publish:
        indexed_bam = SAMTOOLS_INDEX.out
        vcf = GATK_HAPLOTYPECALLER.out.vcf
        vcf_idx = GATK_HAPLOTYPECALLER.out.idx
    }

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

=== "Before"

    ```groovy title="genomics-1.nf" linenums="46"
        publish:
        bam_index = SAMTOOLS_INDEX.out
        vcf = GATK_HAPLOTYPECALLER.out.vcf
        vcf_idx = GATK_HAPLOTYPECALLER.out.idx
    }

    output {
        bam_index {
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

### 3.7. Run the workflow to verify it works correctly on all three samples every time

Run the workflow again to make sure this will work reliably going forward.

```bash
nextflow run genomics-1.nf
```

This time (and every time) everything should run correctly:

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics-1.nf` [special_goldstine] DSL2 - revision: 4cbbf6ea3e

    executor >  local (6)
    [d6/10c2c4] SAMTOOLS_INDEX (1)       | 3 of 3 ✔
    [88/1783aa] GATK_HAPLOTYPECALLER (2) | 3 of 3 ✔
    ```

The results directory now contains both BAM and BAI files for each sample (from the tuple), along with the VCF outputs:

??? abstract "Results directory contents"

    ```console
    results/
    ├── bam/
    │   ├── reads_father.bam -> ...
    │   ├── reads_father.bam.bai -> ...
    │   ├── reads_mother.bam -> ...
    │   ├── reads_mother.bam.bai -> ...
    │   ├── reads_son.bam -> ...
    │   └── reads_son.bam.bai -> ...
    └── vcf/
        ├── reads_father.bam.vcf -> ...
        ├── reads_father.bam.vcf.idx -> ...
        ├── reads_mother.bam.vcf -> ...
        ├── reads_mother.bam.vcf.idx -> ...
        ├── reads_son.bam.vcf -> ...
        └── reads_son.bam.vcf.idx -> ...
    ```

### Takeaway

You know how to make your workflow run on multiple samples (independently).

### What's next?

Make it easier to handle samples in bulk.

---

## 4. Make the workflow accept a text file containing a batch of input files

A very common way to provide multiple data input files to a workflow is to do it with a text file containing the file paths.
It can be as simple as a text file listing one file path per line and nothing else, or the file can contain additional metadata, in which case it's often called a samplesheet.

Here we are going to show you how to do the simple case.

### 4.1. Examine the provided text file listing the input file paths

We already made a text file listing the input file paths, called `sample_bams.txt`, which you can find in the `data/` directory.

```txt title="sample_bams.txt"
/workspaces/training/nf4-science/genomics/data/bam/reads_mother.bam
/workspaces/training/nf4-science/genomics/data/bam/reads_father.bam
/workspaces/training/nf4-science/genomics/data/bam/reads_son.bam
```

As you can see, we listed one file path per line, and they are absolute paths.

!!! note

    The files we are using here are just on your GitHub Codespaces's local filesystem, but we could also point to files in cloud storage.

### 4.2. Update the parameter default

Switch the default value for our `reads_bam` input parameter to point to the `sample_bams.txt` file.

=== "After"

    ```groovy title="genomics-1.nf" linenums="10" hl_lines="1-2"
        // Primary input (file of input files, one per line)
        reads_bam: Path = "${projectDir}/data/sample_bams.txt"
    ```

=== "Before"

    ```groovy title="genomics-1.nf" linenums="10"
    // Primary input (array of three samples)
        reads_bam = [
            "${projectDir}/data/bam/reads_mother.bam",
            "${projectDir}/data/bam/reads_father.bam",
            "${projectDir}/data/bam/reads_son.bam"
        ]
    ```

This way we can continue to be lazy, but the list of files no longer lives in the workflow code itself, which is a big step in the right direction.

### 4.3. Update the channel factory to read lines from a file

Currently, our input channel factory treats any files we give it as the data inputs we want to feed to the indexing process.
Since we're now giving it a file that lists input file paths, we need to change its behavior to parse the file and treat the file paths it contains as the data inputs.

Fortunately we can do that very simply, just by adding the [`.splitText()` operator](https://www.nextflow.io/docs/latest/reference/operator.html#operator-splittext) to the channel construction step.

=== "After"

    ```groovy title="genomics-1.nf" linenums="24" hl_lines="1-2"
        // Create input channel from a text file listing input file paths
        reads_ch = channel.fromPath(params.reads_bam).splitText()
    ```

=== "Before"

    ```groovy title="genomics-1.nf" linenums="24"
        // Create input channel (single file via CLI parameter)
        reads_ch = channel.fromPath(params.reads_bam)
    ```

!!! tip

    This is another great opportunity to use the `.view()` operator to look at what the channel contents look like before and after applying an operator.

### 4.4. Run the workflow to verify that it works correctly

Run the workflow one more time. This should produce the same result as before, right?

```bash
nextflow run genomics-1.nf -resume
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics-1.nf` [sick_albattani] DSL2 - revision: 46d84642f6

    [18/23b4bb] SAMTOOLS_INDEX (1)       | 3 of 3, cached: 3 ✔
    [12/f727bb] GATK_HAPLOTYPECALLER (3) | 3 of 3, cached: 3 ✔
    ```

Yes! In fact, Nextflow correctly detects that the process calls are exactly the same, and doesn't even bother re-running everything, since we were running with `-resume`.

And that's it! Our simple variant calling workflow has all the basic features we wanted.

### Takeaway

You know how to make a multi-step modular workflow to index a BAM file and apply per-sample variant calling using GATK.

More generally, you've learned how to use essential Nextflow components and logic to build a simple genomics pipeline that does real work, taking into account the idiosyncrasies of genomics file formats and tool requirements.

### What's next?

Celebrate your success and take an extra long break!

In the next part of this course, you'll learn how to use a few additional Nextflow features (including more channel operators) to apply joint variant calling to the data.
