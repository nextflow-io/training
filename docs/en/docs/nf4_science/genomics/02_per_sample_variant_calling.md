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

This replicates the steps from Part 1, where you ran these commands manually in their containers.

As a starting point, we provide you with a workflow file, `genomics.nf`, that outlines the main parts of the workflow, as well as two module files, samtools_index.nf and gatk_haplotypecaller.nf, that outline the structure of the modules.
These files are not functional; their purpose is just to serve as scaffolds for you to fill in with the interesting parts of the code.

## Lesson plan

In order to make the development process more educational, we've broken this down into four steps:

1. **Write a single-stage workflow that runs Samtools index on a BAM file.**
   This covers creating a module, importing it, and calling it in a workflow.
2. **Add a second process to run GATK HaplotypeCaller on the indexed BAM file.**
   This introduces chaining process outputs to inputs and handling accessory files.
3. **Adapt the workflow to run on a batch of samples.**
   This covers parallel execution and introduces tuples to keep associated files together.
4. **Make the workflow accept a text file containing a batch of input files.**
   This demonstrates a common pattern for providing inputs in bulk.

Each step focuses on a specific aspect of workflow development.

---

## 1. Write a single-stage workflow that runs Samtools index on a BAM file

This first step focuses on the basics: loading a BAM file and generating an index for it.

Recall the `samtools index` command from [Part 1](01_method.md):

```bash
samtools index '<input_bam>'
```

The command takes a BAM file as input and produces a `.bai` index file alongside it.
The container URI was `community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464`.

We're going to take this information and wrap it in Nextflow in three stages:

1. Set up the input
2. Write the indexing process and call it in the workflow
3. Configure the output handling

### 1.1. Set up the input

We need to declare an input parameter, create a test profile to provide a convenient default value, and create an input channel.

#### 1.1.1. Add an input parameter declaration

In the main workflow file `genomics.nf`, under the `Pipeline parameters` section, declare a CLI parameter called `reads_bam`.

=== "After"

    ```groovy title="genomics.nf" linenums="5" hl_lines="4-7"
    /*
     * Pipeline parameters
     */
    params {
        // Primary input
        reads_bam: Path
    }
    ```

=== "Before"

    ```groovy title="genomics.nf" linenums="5"
    /*
     * Pipeline parameters
     */

    // Primary input
    ```

That sets up the CLI parameter, but we don't want to type out the file path every time we run the workflow during development.
There are multiple options for providing a default value; here we use a test profile.

#### 1.1.2. Create a test profile with a default value in `nextflow.config`

A test profile provides convenient default values for trying out a workflow without specifying inputs on the command line.
This is a common convention in the Nextflow ecosystem (see [Hello Config](../../hello_nextflow/06_hello_config.md) for more detail).

Add a `profiles` block to `nextflow.config` with a `test` profile that sets the `reads_bam` parameter to one of the test BAM files.

=== "After"

    ```groovy title="nextflow.config" linenums="1" hl_lines="3-7"
    docker.enabled = true

    profiles {
        test {
            params.reads_bam = "${projectDir}/data/bam/reads_mother.bam"
        }
    }
    ```

=== "Before"

    ```groovy title="nextflow.config" linenums="1"
    docker.enabled = true
    ```

Here, we're using `${projectDir}`, a built-in Nextflow variable that points to the directory where the workflow script is located.
This makes it easy to reference data files and other resources without hardcoding absolute paths.

#### 1.1.3. Set up the input channel

In the workflow block, create an input channel from the parameter value using the `.fromPath` channel factory (as used in [Hello Channels](../../hello_nextflow/02_hello_channels.md)).

=== "After"

    ```groovy title="genomics.nf" linenums="13" hl_lines="4-5"
    workflow {

        main:
        // Create input channel (single file via CLI parameter)
        reads_ch = channel.fromPath(params.reads_bam)
    ```

=== "Before"

    ```groovy title="genomics.nf" linenums="13"
    workflow {

        main:
        // Create input channel
    ```

Now we need to create the process to run the indexing on this input.

### 1.2. Write the indexing process and call it in the workflow

We need to write the process definition in the module file, import it into the workflow using an include statement, and call it on the input.

#### 1.2.1. Fill in the module for the indexing process

Open `modules/samtools_index.nf` and examine the outline of the process definition.
You should recognize the main structural elements; if not, consider reading through [Hello Nextflow](../../hello_nextflow/01_hello_world.md) for a refresher.

Go ahead and fill in the process definition by yourself using the information provided above, then check your work against the solution in the "After" tab below.

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

Once you've completed this, the process is complete.
To use it in the workflow, you'll need to import the module and add a process call.

#### 1.2.2. Include the module

In `genomics.nf`, add an `include` statement to make the process available to the workflow:

=== "After"

    ```groovy title="genomics.nf" linenums="3" hl_lines="2"
    // Module INCLUDE statements
    include { SAMTOOLS_INDEX } from './modules/samtools_index.nf'
    ```

=== "Before"

    ```groovy title="genomics.nf" linenums="3"
    // Module INCLUDE statements
    ```

The process is now available in the workflow scope.

#### 1.2.3. Call the indexing process on the input

Now, let's add a call to `SAMTOOLS_INDEX` in the workflow block, passing the input channel as an argument.

=== "After"

    ```groovy title="genomics.nf" linenums="14" hl_lines="7-8"
    workflow {

        main:
        // Create input channel (single file via CLI parameter)
        reads_ch = channel.fromPath(params.reads_bam)

        // Create index file for input BAM file
        SAMTOOLS_INDEX(reads_ch)
    ```

=== "Before"

    ```groovy title="genomics.nf" linenums="14"
    workflow {

        main:
        // Create input channel (single file via CLI parameter)
        reads_ch = channel.fromPath(params.reads_bam)

        // Call processes
    ```

The workflow now loads the input and runs the indexing process on it.
Next, we need to configure how the output is published.

### 1.3. Configure the output handling

We need to declare which process outputs to publish and specify where they should go.

#### 1.3.1. Declare an output in the `publish:` section

The `publish:` section inside the workflow block declares which process outputs should be published.
Assign the output of `SAMTOOLS_INDEX` to a named target called `bam_index`.

=== "After"

    ```groovy title="genomics.nf" linenums="22" hl_lines="2"
        publish:
        bam_index = SAMTOOLS_INDEX.out
    }
    ```

=== "Before"

    ```groovy title="genomics.nf" linenums="22"
        publish:
        // Declare outputs to publish
    }
    ```

Now we need to tell Nextflow where to put the published output.

#### 1.3.2. Configure the output target in the `output {}` block

The `output {}` block sits outside the workflow and specifies where each named target is published.
Let's add a target for `bam_index` that publishes into a `bam/` subdirectory.

=== "After"

    ```groovy title="genomics.nf" linenums="26" hl_lines="2-4"
    output {
        bam_index {
            path 'bam'
        }
    }
    ```

=== "Before"

    ```groovy title="genomics.nf" linenums="26"
    output {
        // Configure publish targets
    }
    ```

!!! note

    By default, Nextflow publishes output files as symbolic links, which avoids unnecessary duplication.
    Even though the data files we're using here are very small, in genomics they can get very large.
    Symlinks will break when you clean up your `work` directory, so for production workflows you may want to override the default publish mode to `'copy'`.

### 1.4. Run the workflow

At this point, we have a one-step indexing workflow that should be fully functional. Let's test that it works!

We can run it with `-profile test` to use the default value set up in the test profile and avoid having to write the path on the command line.

```bash
nextflow run genomics.nf -profile test
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics.nf` [reverent_sinoussi] DSL2 - revision: 41d43ad7fe

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

You know how to create a module containing a process, import it into a workflow, call it with an input channel, and publish the results.

### What's next?

Add a second step that takes the output of the indexing process and uses it to run variant calling.

---

## 2. Add a second process to run GATK HaplotypeCaller on the indexed BAM file

Now that we have an index for our input file, we can move on to setting up the variant calling step.

Recall the `gatk HaplotypeCaller` command from [Part 1](01_method.md):

```bash
gatk HaplotypeCaller \
        -R /data/ref/ref.fasta \
        -I /data/bam/reads_mother.bam \
        -O reads_mother.vcf \
        -L /data/ref/intervals.bed
```

The command takes a BAM file (`-I`), a reference genome (`-R`), and an intervals file (`-L`), and produces a VCF file (`-O`) along with its index.
The tool also expects the BAM index, reference index, and reference dictionary to be co-located with their respective files.
The container URI was `community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867`.

We follow the same three stages as before:

1. Set up the inputs
2. Write the variant calling process and call it in the workflow
3. Configure the output handling

### 2.1. Set up the inputs

The variant calling step requires several additional input files.
We need to declare parameters for them, add default values to the test profile, and create variables to load them.

#### 2.1.1. Add parameter declarations for accessory inputs

Since our new process expects a handful of additional files to be provided, add parameter declarations for them in `genomics.nf` under the `Pipeline parameters` section:

=== "After"

    ```groovy title="genomics.nf" linenums="9" hl_lines="5-9"
    params {
        // Primary input
        reads_bam: Path

        // Accessory files
        reference: Path
        reference_index: Path
        reference_dict: Path
        intervals: Path
    }
    ```

=== "Before"

    ```groovy title="genomics.nf" linenums="9"
    params {
        // Primary input
        reads_bam: Path
    }
    ```

As before, we provide default values through the test profile rather than inline.

#### 2.1.2. Add accessory file defaults to the test profile

Just as we did for `reads_bam` in section 1.1.2, add default values for the accessory files to the test profile in `nextflow.config`:

=== "After"

    ```groovy title="nextflow.config" linenums="4" hl_lines="3-6"
    test {
        params.reads_bam = "${projectDir}/data/bam/reads_mother.bam"
        params.reference = "${projectDir}/data/ref/ref.fasta"
        params.reference_index = "${projectDir}/data/ref/ref.fasta.fai"
        params.reference_dict = "${projectDir}/data/ref/ref.dict"
        params.intervals = "${projectDir}/data/ref/intervals.bed"
    }
    ```

=== "Before"

    ```groovy title="nextflow.config" linenums="4"
    test {
        params.reads_bam = "${projectDir}/data/bam/reads_mother.bam"
    }
    ```

Now we need to create variables that load these file paths for use in the workflow.

#### 2.1.3. Create variables for the accessory files

Add variables for the accessory file paths inside the workflow block:

=== "After"

    ```groovy title="genomics.nf" linenums="21" hl_lines="7-11"
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

    ```groovy title="genomics.nf" linenums="21"
    workflow {

        main:
        // Create input channel (single file via CLI parameter)
        reads_ch = channel.fromPath(params.reads_bam)

        // Create index file for input BAM file
        SAMTOOLS_INDEX(reads_ch)
    ```

The `file()` syntax tells Nextflow explicitly to handle these inputs as file paths.
You can learn more about this in the Side Quest [Working with files](../../side_quests/working_with_files.md).

### 2.2. Write the variant calling process and call it in the workflow

We need to write the process definition in the module file, import it into the workflow using an include statement, and call it on the input reads plus the output of the indexing step and the accessory files.

#### 2.2.1. Fill in the module for the variant calling process

Open `modules/gatk_haplotypecaller.nf` and examine the outline of the process definition.

Go ahead and fill in the process definition by yourself using the information provided above, then check your work against the solution in the "After" tab below.

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

You'll notice that this process has more inputs than the GATK command itself requires.
The GATK knows to look for the BAM index file and the reference genome's accessory files based on naming conventions, but Nextflow is domain-agnostic and doesn't know about these conventions.
We need to list them explicitly so that Nextflow stages them in the working directory at runtime; otherwise GATK will throw an error about missing files.

Similarly, we list the output VCF's index file (`"${input_bam}.vcf.idx"`) explicitly so that Nextflow keeps track of it for subsequent steps.
We use the `emit:` syntax to assign a name to each output channel, which will become useful when we wire the outputs into the publish block.

Once you've completed this, the process is complete.
To use it in the workflow, you'll need to import the module and add a process call.

#### 2.2.2. Import the new module

Update `genomics.nf` to import the new module:

=== "After"

    ```groovy title="genomics.nf" linenums="3" hl_lines="3"
    // Module INCLUDE statements
    include { SAMTOOLS_INDEX } from './modules/samtools_index.nf'
    include { GATK_HAPLOTYPECALLER } from './modules/gatk_haplotypecaller.nf'
    ```

=== "Before"

    ```groovy title="genomics.nf" linenums="3"
    // Module INCLUDE statements
    include { SAMTOOLS_INDEX } from './modules/samtools_index.nf'
    ```

The process is now available in the workflow scope.

#### 2.2.3. Add the process call

Add the process call in the workflow body, under `main:`:

=== "After"

    ```groovy title="genomics.nf" linenums="33" hl_lines="4-12"
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

    ```groovy title="genomics.nf" linenums="33"
        // Create index file for input BAM file
        SAMTOOLS_INDEX(reads_ch)
    ```

You should recognize the `*.out` syntax from the Hello Nextflow training series; we are telling Nextflow to take the channel output by `SAMTOOLS_INDEX` and plugging that into the `GATK_HAPLOTYPECALLER` process call.

!!! note

    Notice that the inputs are provided in the exact same order in the call to the process as they are listed in the input block of the process.
    In Nextflow, inputs are positional, meaning you _must_ follow the same order; and of course there have to be the same number of elements.

### 2.3. Configure the output handling

We need to add the new outputs to the publish declaration and configure where they go.

#### 2.3.1. Add publish targets for the variant calling outputs

Add the VCF and index outputs to the `publish:` section:

=== "After"

    ```groovy title="genomics.nf" linenums="45" hl_lines="3-4"
        publish:
        bam_index = SAMTOOLS_INDEX.out
        vcf = GATK_HAPLOTYPECALLER.out.vcf
        vcf_idx = GATK_HAPLOTYPECALLER.out.idx
    }
    ```

=== "Before"

    ```groovy title="genomics.nf" linenums="45"
        publish:
        bam_index = SAMTOOLS_INDEX.out
    }
    ```

Now we need to tell Nextflow where to put the new outputs.

#### 2.3.2. Configure the new output targets

Add entries for the `vcf` and `vcf_idx` targets in the `output {}` block, publishing both into a `vcf/` subdirectory:

=== "After"

    ```groovy title="genomics.nf" linenums="51" hl_lines="5-10"
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

    ```groovy title="genomics.nf" linenums="49"
    output {
        bam_index {
            path 'bam'
        }
    }
    ```

The VCF and its index are published as separate targets that both go into the `vcf/` subdirectory.

### 2.4. Run the workflow

Run the expanded workflow, adding `-resume` this time so that we don't have to run the indexing step again.

```bash
nextflow run genomics.nf -profile test -resume
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics.nf` [grave_volta] DSL2 - revision: 4790abc96a

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

??? abstract "File contents"

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

### 3.1. Update the input to list three samples

To run on multiple samples, update the test profile to provide an array of file paths instead of a single one.
This is a quick way to test multi-sample execution; in the next step we'll switch to a more scalable approach using a file of inputs.

First, comment out the type annotation in the parameter declaration, since arrays cannot use typed declarations:

=== "After"

    ```groovy title="genomics.nf" linenums="10" hl_lines="1-2"
        // Primary input (array of three samples)
        reads_bam //: Path
    ```

=== "Before"

    ```groovy title="genomics.nf" linenums="10"
        // Primary input
        reads_bam: Path
    ```

Then update the test profile to list all three samples:

=== "After"

    ```groovy title="nextflow.config" linenums="4" hl_lines="2-6"
    test {
        params.reads_bam = [
            "${projectDir}/data/bam/reads_mother.bam",
            "${projectDir}/data/bam/reads_father.bam",
            "${projectDir}/data/bam/reads_son.bam"
        ]
        params.reference = "${projectDir}/data/ref/ref.fasta"
        params.reference_index = "${projectDir}/data/ref/ref.fasta.fai"
        params.reference_dict = "${projectDir}/data/ref/ref.dict"
        params.intervals = "${projectDir}/data/ref/intervals.bed"
    }
    ```

=== "Before"

    ```groovy title="nextflow.config" linenums="4"
    test {
        params.reads_bam = "${projectDir}/data/bam/reads_mother.bam"
        params.reference = "${projectDir}/data/ref/ref.fasta"
        params.reference_index = "${projectDir}/data/ref/ref.fasta.fai"
        params.reference_dict = "${projectDir}/data/ref/ref.dict"
        params.intervals = "${projectDir}/data/ref/intervals.bed"
    }
    ```

The channel factory in the workflow body (`.fromPath`) accepts multiple file paths just as well as a single one, so no other changes are needed.

### 3.2. Run the workflow

Try running the workflow now that the plumbing is set up to run on all three test samples.

```bash
nextflow run genomics.nf -profile test -resume
```

Funny thing: this _might work_, OR it _might fail_. For example, here's a run that succeeded:

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics.nf` [peaceful_yalow] DSL2 - revision: a256d113ad

    executor >  local (6)
    [4f/7071b0] SAMTOOLS_INDEX (3)       | 3 of 3, cached: 1 ✔
    [7a/89bc43] GATK_HAPLOTYPECALLER (2) | 3 of 3, cached: 1 ✔
    ```

If your workflow run succeeded, run it again until you get an error like this:

??? failure "Command output"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics.nf` [loving_pasteur] DSL2 - revision: d2a8e63076

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

### 3.3. Troubleshoot the problem

We'll inspect the work directories and use the `view()` operator to figure out what went wrong.

#### 3.3.1. Check the work directories for the relevant calls

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

#### 3.3.2. Use the [view() operator](https://www.nextflow.io/docs/latest/reference/operator.html#view) to inspect channel contents

Add these two lines in the workflow body before the `GATK_HAPLOTYPECALLER` process call to view the contents of the channel:

=== "After"

    ```groovy title="genomics.nf" hl_lines="3-5"
        SAMTOOLS_INDEX(reads_ch)

        // temporary diagnostics
        reads_ch.view()
        SAMTOOLS_INDEX.out.view()

        // Call variants from the indexed BAM file
        GATK_HAPLOTYPECALLER(
    ```

=== "Before"

    ```groovy title="genomics.nf"
        SAMTOOLS_INDEX(reads_ch)

        // Call variants from the indexed BAM file
        GATK_HAPLOTYPECALLER(
    ```

Then run the workflow command again.

```bash
nextflow run genomics.nf -profile test
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

### 3.4. Update the workflow to handle the index files correctly

The fix is to bundle each BAM file with its index into a tuple, then update the downstream process and workflow plumbing to match.

#### 3.4.1. Change the output of the SAMTOOLS_INDEX module to a tuple

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

#### 3.4.2. Change the input of the GATK_HAPLOTYPECALLER module to accept a tuple

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

Now we need to update the workflow to reflect the new tuple structure in the process call and the publish targets.

#### 3.4.3. Update the call to GATK_HAPLOTYPECALLER in the workflow

We no longer need to provide the original `reads_ch` to the `GATK_HAPLOTYPECALLER` process, since the BAM file is now bundled into the channel output by `SAMTOOLS_INDEX`.

Update the call in `genomics.nf`:

=== "After"

    ```groovy title="genomics.nf" linenums="42" hl_lines="2"
        GATK_HAPLOTYPECALLER(
            SAMTOOLS_INDEX.out,
    ```

=== "Before"

    ```groovy title="genomics.nf" linenums="42"
        GATK_HAPLOTYPECALLER(
            reads_ch,
            SAMTOOLS_INDEX.out,
    ```

Finally, we need to update the publish targets to reflect the new output structure.

#### 3.4.4. Update the publish target for the indexed BAM output

Since the SAMTOOLS_INDEX output is now a tuple containing both the BAM file and its index, rename the publish target from `bam_index` to `indexed_bam` to better reflect its contents:

=== "After"

    ```groovy title="genomics.nf" linenums="46" hl_lines="2 8"
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

    ```groovy title="genomics.nf" linenums="46"
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

With these changes, the BAM and its index are guaranteed to travel together, so the pairing will always be correct.

### 3.5. Run the corrected workflow

Run the workflow again to make sure this will work reliably going forward.

```bash
nextflow run genomics.nf -profile test
```

This time (and every time) everything should run correctly:

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics.nf` [special_goldstine] DSL2 - revision: 4cbbf6ea3e

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

By bundling associated files into tuples, we ensured the correct files always travel together through the workflow.
The workflow now processes any number of samples reliably, but listing them individually in the config is not very scalable.
In the next step, we'll switch to reading inputs from a file.

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
    If you are not using the provided Codespaces environment, you may need to adapt the file paths to match your local setup.

### 4.2. Update the parameter and test profile

Switch the `reads_bam` parameter to point to the `sample_bams.txt` file instead of listing individual samples.

Restore the type annotation in the params block (since it's a single path again):

=== "After"

    ```groovy title="genomics.nf" linenums="10" hl_lines="1-2"
        // Primary input (file of input files, one per line)
        reads_bam: Path
    ```

=== "Before"

    ```groovy title="genomics.nf" linenums="10"
        // Primary input (array of three samples)
        reads_bam
    ```

Then update the test profile to point to the text file:

=== "After"

    ```groovy title="nextflow.config" linenums="4" hl_lines="2"
    test {
        params.reads_bam = "${projectDir}/data/sample_bams.txt"
        params.reference = "${projectDir}/data/ref/ref.fasta"
        params.reference_index = "${projectDir}/data/ref/ref.fasta.fai"
        params.reference_dict = "${projectDir}/data/ref/ref.dict"
        params.intervals = "${projectDir}/data/ref/intervals.bed"
    }
    ```

=== "Before"

    ```groovy title="nextflow.config" linenums="4"
    test {
        params.reads_bam = [
            "${projectDir}/data/bam/reads_mother.bam",
            "${projectDir}/data/bam/reads_father.bam",
            "${projectDir}/data/bam/reads_son.bam"
        ]
        params.reference = "${projectDir}/data/ref/ref.fasta"
        params.reference_index = "${projectDir}/data/ref/ref.fasta.fai"
        params.reference_dict = "${projectDir}/data/ref/ref.dict"
        params.intervals = "${projectDir}/data/ref/intervals.bed"
    }
    ```

The list of files no longer lives in the code at all, which is a big step in the right direction.

### 4.3. Update the channel factory to read lines from a file

Currently, our input channel factory treats any files we give it as the data inputs we want to feed to the indexing process.
Since we're now giving it a file that lists input file paths, we need to change its behavior to parse the file and treat the file paths it contains as the data inputs.

We can do this using the same pattern we used in [Part 2 of Hello Nextflow](../../hello_nextflow/02_hello_channels.md#42-use-the-splitcsv-operator-to-parse-the-file): applying the [`splitCsv()`](https://nextflow.io/docs/latest/reference/operator.html#splitcsv) operator to parse the file, then a `map` operation to select the first field of each line.

=== "After"

    ```groovy title="genomics.nf" linenums="24" hl_lines="1-4"
        // Create input channel from a CSV file listing input file paths
        reads_ch = Channel.fromPath(params.reads_bam)
                .splitCsv()
                .map { line -> file(line[0]) }
    ```

=== "Before"

    ```groovy title="genomics.nf" linenums="24"
        // Create input channel (single file via CLI parameter)
        reads_ch = channel.fromPath(params.reads_bam)
    ```

Technically we could do this more simply using the [`.splitText()`](https://www.nextflow.io/docs/latest/reference/operator.html#operator-splittext) operator, since our input file currently only contains file paths.
However, by using the more versatile `splitCsv` operator (supplemented by `map`), we can future-proof our workflow in case we decide to add metadata to the file containing file paths.

!!! tip

    If you're not confident you understand what the operators are doing here, this is another great opportunity to use the `.view()` operator to look at what the channel contents look like before and after applying them.

### 4.4. Run the workflow

Run the workflow one more time. This should produce the same result as before, right?

```bash
nextflow run genomics.nf -profile test -resume
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics.nf` [sick_albattani] DSL2 - revision: 46d84642f6

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

In the next part of this course, you'll learn how to transform this simple per-sample variant calling workflow to apply joint variant calling to the data.
