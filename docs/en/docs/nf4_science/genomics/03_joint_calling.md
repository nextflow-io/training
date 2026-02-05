# Part 3: Joint calling on a cohort

In Part 2, you built a variant calling pipeline that was completely linear and processed each sample's data independently of the others.
However, in a real genomics use case, you'll typically need to look at the variant calls of multiple samples together.

In this part, we show you how to use channels and channel operators to implement joint variant calling with GATK, building on the pipeline from Part 2.

### Method overview

The GATK variant calling method we used in Part 2 simply generated variant calls per sample.
That's fine if you only want to look at the variants from each sample in isolation, but that yields limited information.
It's often more interesting to look at how variant calls differ across multiple samples, and to do so, GATK offers an alternative method called joint variant calling, which we demonstrate here.

Joint variant calling involves generating a special kind of variant output called GVCF (for Genomic VCF) for each sample, then combining the GVCF data from all the samples and finally, running a 'joint genotyping' statistical analysis.

![Joint analysis](img/joint-calling.png)

What's special about a sample's GVCF is that it contains records summarizing sequence data statistics about all positions in the targeted area of the genome, not just the positions where the program found evidence of variation.
This is critical for the joint genotyping calculation ([further reading](https://gatk.broadinstitute.org/hc/en-us/articles/360035890431-The-logic-of-joint-calling-for-germline-short-variants)).

The GVCF is produced by GATK HaplotypeCaller, the same tool we tested in Part 1, with an additional parameter (`-ERC GVCF`).
Combining the GVCFs is done with GATK GenomicsDBImport, which combines the per-sample calls into a data store (analogous to a database), then the actual 'joint genotyping' analysis is done with GATK GenotypeGVCFs.

### Workflow

So to recap, in this part of the course, we're going to develop a workflow that does the following:

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/genomics/img/hello-gatk-2.svg"
</figure>

1. Generate an index file for each BAM input file using Samtools
2. Run the GATK HaplotypeCaller on each BAM input file to generate a GVCF of per-sample genomic variant calls
3. Collect all the GVCFs and combine them into a GenomicsDB data store
4. Run joint genotyping on the combined GVCF data store to produce a cohort-level VCF

We'll apply this to the same dataset as in Parts 1 and 2.

!!! note

     Make sure you're in the correct working directory:
     `cd /workspaces/training/nf4-science/genomics`

---

## 1. Modify the per-sample variant calling step to produce a GVCF

The good news is that we don't need to start all over, since we already wrote a workflow that does some of this work in Part 2.
However, that pipeline produces VCF files, whereas now we want GVCF files in order to do the joint genotyping.
So we need to start by switching on the GVCF variant calling mode and updating the output file extension.

!!! note

    For convenience, we are going to work with a fresh copy of the GATK workflow as it stands at the end of Part 2, but under a different name: `genomics-2.nf`.

### 1.1. Tell HaplotypeCaller to emit a GVCF and update the output extension

Let's open the `genomics-2.nf` file in the code editor.
It should look very familiar, but feel free to run it if you want to satisfy yourself that it runs as expected.

We're going to start by making two changes:

- Add the `-ERC GVCF` parameter to the GATK HaplotypeCaller command;
- Update the output file path to use the corresponding `.g.vcf` extension, as per GATK convention.

Make sure you add a backslash (`\`) at the end of the previous line when you add `-ERC GVCF`.

=== "After"

    ```groovy title="genomics-2.nf" linenums="56" hl_lines="4 6"
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

    ```groovy title="genomics-2.nf" linenums="56" hl_lines="4"
        """
        gatk HaplotypeCaller \
            -R ${ref_fasta} \
            -I ${input_bam} \
            -O ${input_bam}.vcf \
            -L ${interval_list}
        """
    ```

And that's all it takes to switch HaplotypeCaller to generating GVCFs instead of VCFs, right?

### 1.2. Run the pipeline to verify that you can generate GVCFs

The Nextflow execution command is the same as before, save for the workflow filename itself.
Make sure to update that appropriately.

```bash
nextflow run genomics-2.nf
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics-2.nf` [crazy_venter] DSL2 - revision: a2d6f6f09f

    executor >  local (6)
    [f1/8d8486] SAMTOOLS_INDEX (1)       | 3 of 3 ✔
    [72/3249ca] GATK_HAPLOTYPECALLER (3) | 0 of 3
    ERROR ~ Error executing process > 'GATK_HAPLOTYPECALLER (2)'

    Caused by:
      Missing output file(s) `reads_son.bam.vcf` expected by process `GATK_HAPLOTYPECALLER (2)`

    Command executed:

      gatk HaplotypeCaller         -R ref.fasta         -I reads_son.bam         -O reads_son.bam.g.vcf         -L intervals.bed         -ERC GVCF
    ```

And the output is... all red! Oh no.

The command that was executed is correct, so we were right that that was enough to change the GATK tool's behavior.
But look at that line about the missing output file. Notice anything?

That's right, we forgot to tell Nextflow to expect a new file name. Oops.

### 1.3. Update the output file extension in the process outputs block too

Because it's not enough to just change the file extension in the tool command itself, you also have to tell Nextflow that the expected output filename has changed.

=== "After"

    ```groovy title="genomics-2.nf" linenums="50" hl_lines="2 3"
        output:
        path "${input_bam}.g.vcf"     , emit: vcf
        path "${input_bam}.g.vcf.idx" , emit: idx
    ```

=== "Before"

    ```groovy title="genomics-2.nf" linenums="50" hl_lines="2 3"
        output:
        path "${input_bam}.vcf"     , emit: vcf
        path "${input_bam}.vcf.idx" , emit: idx
    ```

### 1.4. Update the publish targets for the new GVCF outputs

Since we're now producing GVCFs instead of VCFs, we should update the workflow's `publish:` section to use more descriptive names.
We'll also organize the GVCF files into their own subdirectory for clarity.

=== "After"

    ```groovy title="genomics-2.nf" linenums="88" hl_lines="3 4"
        publish:
        indexed_bam = SAMTOOLS_INDEX.out
        gvcf = GATK_HAPLOTYPECALLER.out.vcf
        gvcf_idx = GATK_HAPLOTYPECALLER.out.idx
    ```

=== "Before"

    ```groovy title="genomics-2.nf" linenums="88"
        publish:
        indexed_bam = SAMTOOLS_INDEX.out
        vcf = GATK_HAPLOTYPECALLER.out.vcf
        vcf_idx = GATK_HAPLOTYPECALLER.out.idx
    ```

### 1.5. Update the output block for the new directory structure

We also need to update the `output` block to put the GVCF files in a `gvcf` subdirectory.

=== "After"

    ```groovy title="genomics-2.nf" linenums="94" hl_lines="3 5 6 8 9"
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

    ```groovy title="genomics-2.nf" linenums="94"
    output {
        indexed_bam {
            path '.'
        }
        vcf {
            path '.'
        }
        vcf_idx {
            path '.'
        }
    }
    ```

### 1.6. Run the pipeline again

Let's run it with `-resume` this time.

```bash
nextflow run genomics-2.nf -resume
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics-2.nf` [nostalgic_franklin] DSL2 - revision: f2c0a93c6a

    executor >  local (3)
    [cc/fbc705] SAMTOOLS_INDEX (3)       | 3 of 3, cached: 3 ✔
    [27/0d7eb9] GATK_HAPLOTYPECALLER (2) | 3 of 3 ✔
    ```

This time it works.

The Nextflow output itself doesn't look any different (compared to a successful run in normal VCF mode), but now we can find the `.g.vcf` files and their respective index files, for all three samples, organized in subdirectories.

??? abstract "Directory contents (symlinks shortened)"

    ```console
    results_genomics/
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

Okay, this one was minimal in terms of Nextflow learning...
But it was a nice opportunity to reiterate the importance of the process output block!

### What's next?

Learn to collect the contents of a channel and pass them on to the next process as a single input.

---

## 2. Collect and combine the GVCF data across all samples

We now need to combine the data from all the per-sample GVCFs into a form that supports the joint genotyping analysis we want to do.

### 2.1. Create a module for the process that will combine the GVCFs

As a reminder of what we covered earlier in Part 1, combining the GVCFs is a job for the GATK tool GenomicsDBImport, which will produce a data store in the so-called GenomicsDB format.

Following the same pattern we used in Part 2, we'll create this as a module file.

```bash
mkdir -p modules/gatk/genomicsdb
```

Now create the module file:

```groovy title="modules/gatk/genomicsdb/main.nf"
/*
 * Combine GVCFs into GenomicsDB datastore
 */
process GATK_GENOMICSDB {

    container "community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867"

    input:
    path all_gvcfs
    path all_idxs
    path interval_list
    val cohort_name

    output:
    path "${cohort_name}_gdb"

    script:
    """
    gatk GenomicsDBImport \
        -V ${all_gvcfs} \
        -L ${interval_list} \
        --genomicsdb-workspace-path ${cohort_name}_gdb
    """
}
```

What do you think, looks reasonable?

### 2.2. Import the module into the workflow

Add the import statement to `genomics-2.nf`, below the existing import statements:

```groovy title="genomics-2.nf" linenums="10"
include { GATK_GENOMICSDB } from './modules/gatk/genomicsdb/main.nf'
```

Now we can wire it up and see what happens.

### 2.3. Add a `cohort_name` parameter with a default value

We need to provide an arbitrary name for the cohort.
Later in the training series you'll learn how to use sample metadata for this sort of thing, but for now we just declare a CLI parameter using `params` and give it a default value for convenience.

```groovy title="genomics-2.nf" linenums="16"
    // Base name for final output file
    cohort_name: String = "family_trio"
```

### 2.4. Gather the outputs of GATK_HAPLOTYPECALLER across samples

If we were to just plug the output channel from the `GATK_HAPLOTYPECALLER` process as is, Nextflow would call the process on each sample GVCF separately.
However, we want to bundle all three GVCFs (and their index files) in such a way that Nextflow hands all of them together to a single process call.

Good news: we can do that using the `collect()` channel operator. Let's add the following lines to the `workflow` body, right after the call to GATK_HAPLOTYPECALLER:

```groovy title="genomics-2.nf" linenums="118"
// Collect variant calling outputs across samples
all_gvcfs_ch = GATK_HAPLOTYPECALLER.out.vcf.collect()
all_idxs_ch = GATK_HAPLOTYPECALLER.out.idx.collect()
```

Does that seem a bit complicated? Let's break this down and translate it into plain language.

1. We're taking the output channel from the `GATK_HAPLOTYPECALLER` process, referred to using the `.out` property.
2. Each 'element' coming out of the channel is a pair of files: the GVCF and its index file, in that order because that's the order they're listed in the process output block. Conveniently, because in the last session we named the outputs of this process (using `emit:`), we can pick out the GVCFs on one hand by adding `.vcf` and the index files on the other by adding `.idx` after the `.out` property. If we had not named those outputs, we would have had to refer to them by `.out[0]` and `.out[1]`, respectively.
3. We append the `collect()` channel operator to bundle all the GVCF files together into a single element in a new channel called `all_gvcfs_ch`, and do the same with the index files to form the new channel called `all_idxs_ch`.

!!! tip

    If you're having a hard time envisioning exactly what is happening here, remember that you can use the `view()` operator to inspect the contents of channels before and after applying channel operators.

The resulting `all_gvcfs_ch` and `all_idxs_ch` channels are what we're going to plug into the `GATK_GENOMICSDB` process we just wrote.

!!! note

    In case you were wondering, we collect the GVCFs and their index files separately because the GATK GenomicsDBImport command only wants to see the GVCF file paths. Fortunately, since Nextflow will stage all the files together for execution, we don't have to worry about the order of files like we did for BAMs and their index in Part 1.

### 2.5. Add a call to the workflow block to run GATK_GENOMICSDB

We've got a process, and we've got input channels. We just need to add the process call.

```groovy title="genomics-2.nf" linenums="122"
    // Combine GVCFs into a GenomicsDB datastore
    GATK_GENOMICSDB(
        all_gvcfs_ch,
        all_idxs_ch,
        intervals_file,
        params.cohort_name
    )
```

Ok, everything is wired up.

### 2.6. Run the workflow

Let's see if this works.

```bash
nextflow run genomics-2.nf -resume
```

??? failure "Command output"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics-2.nf` [disturbed_bell] DSL2 - revision: 57942246cc

    executor >  local (1)
    [f1/8d8486] SAMTOOLS_INDEX (1)       | 3 of 3, cached: 3 ✔
    [e4/4ed55e] GATK_HAPLOTYPECALLER (2) | 3 of 3, cached: 3 ✔
    [51/d350ea] GATK_GENOMICSDB          | 0 of 1
    ERROR ~ Error executing process > 'GATK_GENOMICSDB'

    Caused by:
      Process `GATK_GENOMICSDB` terminated with an error exit status (1)

    Command executed:

      gatk GenomicsDBImport         -V reads_son.bam.g.vcf reads_father.bam.g.vcf reads_mother.bam.g.vcf         -L intervals.bed         --genomicsdb-workspace-path family_trio_gdb
    ```

It runs fairly quickly, since we're running with `-resume`, but it fails!

Ah. On the bright side, we see that Nextflow has picked up the `GATK_GENOMICSDB` process, and specifically called it just once.
That suggests that the `collect()` approach worked, to a point.
But, and it's a big one, the process call failed.

When we dig into the console output above, we can see the command executed isn't correct.

Can you spot the error?
Look at this bit: `-V reads_father.bam.g.vcf reads_son.bam.g.vcf reads_mother.bam.g.vcf`

We gave `gatk GenomicsDBImport` multiple GVCF files for a single `-V` argument, but the tool expects a separate `-V` argument for each GVCF file.

As a reminder, this was the command we ran in the container:

```bash
gatk GenomicsDBImport \
    -V reads_mother.g.vcf \
    -V reads_father.g.vcf \
    -V reads_son.g.vcf \
    -L /data/ref/intervals.bed \
    --genomicsdb-workspace-path family_trio_gdb
```

So that means we need to somehow transform our bundle of GVCF files into a properly formatted command string.

### 2.7. Construct a command line with a separate `-V` argument for each input GVCF

This is where Nextflow being based on Groovy comes in handy, because it's going to allow us to use some fairly straightforward string manipulations to construct the necessary command string.

Specifically, using this syntax: `all_gvcfs.collect { gvcf -> "-V ${gvcf}" }.join(' ')`

Once again, let's break it down into its components.

1. First, we take the contents of the `all_gvcfs` input channel and apply `.collect()` on it (just like earlier).
2. That allows us to pass each individual GVCF file path in the bundle to the **closure**, `{ gvcf -> "-V ${gvcf}" }`, where `gvcf` refers to that GVCF file path.
   The closure is a mini-function that we use to prepend `-V ` to the file path, in the form of `"-V ${gvcf}"`.
3. Then we use `.join(' ')` to concatenate all three strings with a single space as separator.

With a concrete example, it looks like this:

1. We have three files:

   `[A.ext, B.ext, C.ext]`

2. The closure modifies each one to create the strings:

   `"-V A.ext", "-V B.ext", "-V C.ext"`

3. The `.join(' ')` operation generates the final string:

   `"-V A.ext -V B.ext -V C.ext"`

Once we have that string, we can assign it to a local variable, `gvcfs_line`, defined with the `def` keyword:

`def gvcfs_line = all_gvcfs.collect { gvcf -> "-V ${gvcf}" }.join(' ')`

Ok, so we have our string manipulation thingy. Where do we put it?

We want this to go inside the process definition somewhere, because we want to do it _after_ we've channeled the GVCF file paths into the process.
That is because Nextflow must see them as file paths in order to stage the files themselves correctly for execution.

But _where_ in the process can we add this?

Fun fact: you can add arbitrary code after `script:` and before the `"""` !

Great, let's add our string manipulation line there then, and update the `gatk GenomicsDBImport` command to use the concatenated string it produces.

=== "After"

    ```groovy title="modules/gatk/genomicsdb/main.nf" linenums="17"  hl_lines="2 5"
        script:
        def gvcfs_line = all_gvcfs.collect { gvcf -> "-V ${gvcf}" }.join(' ')
        """
        gatk GenomicsDBImport \
            ${gvcfs_line} \
            -L ${interval_list} \
            --genomicsdb-workspace-path ${cohort_name}_gdb
        """
    ```

=== "Before"

    ```groovy title="modules/gatk/genomicsdb/main.nf" linenums="17"  hl_lines="4"
        script:
        """
        gatk GenomicsDBImport \
            -V ${all_gvcfs} \
            -L ${interval_list} \
            --genomicsdb-workspace-path ${cohort_name}_gdb
        """
    ```

That should be all that's needed to provide the inputs to `gatk GenomicsDBImport` correctly.

!!! tip

    When you update the `gatk GenomicsDBImport` command, make sure to remove the `-V ` prefix when you swap in the `${gvcfs_line}` variable.

### 2.8. Run the workflow to verify that it generates the GenomicsDB output as expected

Alright, let's see if that addressed the issue.

```bash
nextflow run genomics-2.nf -resume
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics-2.nf` [peaceful_gates] DSL2 - revision: ca0bf847ed

    executor >  local (1)
    [cc/fbc705] SAMTOOLS_INDEX (3)       | 3 of 3, cached: 3 ✔
    [27/0d7eb9] GATK_HAPLOTYPECALLER (2) | 3 of 3, cached: 3 ✔
    [76/d13861] GATK_GENOMICSDB          | 1 of 1 ✔
    ```

Aha! It seems to be working now.

The first two steps were successfully skipped, and the third step worked like a charm this time.
The GenomicsDB data store is created in the work directory but not published to results, since it's just an intermediate format that we'll use for joint genotyping.

By the way, we didn't have to do anything special to handle the output being a directory instead of a single file.

### Takeaway

Now you know how to collect outputs from a channel and bundle them as a single input to another process.
You also know how to construct a command line to provide inputs to a given tool with the appropriate syntax.

### What's next?

Learn how to add a second command to the same process.

---

## 3. Run the joint genotyping step as part of the same process

Now that we have the combined genomic variant calls, we can run the joint genotyping tool, which will produce the final output that we actually care about: the VCF of cohort-level variant calls.

For logistical reasons, we decide to include the joint genotyping inside the same process.

### 3.1. Rename the module from GATK_GENOMICSDB to GATK_JOINTGENOTYPING

Since the process will be running more than one tool, we change its name to refer to the overall operation rather than a single tool name.
This means renaming both the module directory and the process inside it.

First, rename the module directory:

```bash
mv modules/gatk/genomicsdb modules/gatk/jointgenotyping
```

Then update the process name in the module file:

=== "After"

    ```groovy title="modules/gatk/jointgenotyping/main.nf"
    /*
     * Combine GVCFs into GenomicsDB datastore and run joint genotyping to produce cohort-level calls
     */
    process GATK_JOINTGENOTYPING {
    ```

=== "Before"

    ```groovy title="modules/gatk/jointgenotyping/main.nf"
    /*
     * Combine GVCFs into GenomicsDB datastore
     */
    process GATK_GENOMICSDB {
    ```

Finally, update the import statement in `genomics-2.nf`:

=== "After"

    ```groovy title="genomics-2.nf" linenums="10"
    include { GATK_JOINTGENOTYPING } from './modules/gatk/jointgenotyping/main.nf'
    ```

=== "Before"

    ```groovy title="genomics-2.nf" linenums="10"
    include { GATK_GENOMICSDB } from './modules/gatk/genomicsdb/main.nf'
    ```

Remember to keep your process names as descriptive as possible, to maximize readability for your colleagues —and your future self!

### 3.2. Add the joint genotyping command to the GATK_JOINTGENOTYPING module

Simply add the second command after the first one inside the script section of the module file.

=== "After"

    ```groovy title="modules/gatk/jointgenotyping/main.nf" linenums="19"  hl_lines="6-10"
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
    ```

=== "Before"

    ```groovy title="modules/gatk/jointgenotyping/main.nf" linenums="19"
        """
        gatk GenomicsDBImport \
            ${gvcfs_line} \
            -L ${interval_list} \
            --genomicsdb-workspace-path ${cohort_name}_gdb
        """
    ```

The two commands will be run in serial, in the same way that they would if we were to run them manually in the terminal.

### 3.3. Add the reference genome files to the GATK_JOINTGENOTYPING module input definitions

The second command requires the reference genome files, so we need to add those to the module's input block.

=== "After"

    ```groovy title="modules/gatk/jointgenotyping/main.nf" linenums="8"  hl_lines="6-8"
        input:
        path all_gvcfs
        path all_idxs
        path interval_list
        val cohort_name
        path ref_fasta
        path ref_index
        path ref_dict
    ```

=== "Before"

    ```groovy title="modules/gatk/jointgenotyping/main.nf" linenums="8"
        input:
        path all_gvcfs
        path all_idxs
        path interval_list
        val cohort_name
    ```

It may seem annoying to type these out, but remember, you only type them once, and then you can run the workflow a million times. Worth it?

### 3.4. Update the module output definition to emit the VCF of cohort-level variant calls

We don't really care about saving the GenomicsDB datastore, which is just an intermediate format that only exists for logistical reasons, so we can just remove it from the output block if we want.

The output we're actually interested in is the VCF produced by the joint genotyping command.

=== "After"

    ```groovy title="modules/gatk/jointgenotyping/main.nf" linenums="17" hl_lines="2 3"
        output:
        path "${cohort_name}.joint.vcf"     , emit: vcf
        path "${cohort_name}.joint.vcf.idx" , emit: idx
    ```

=== "Before"

    ```groovy title="modules/gatk/jointgenotyping/main.nf" linenums="17"
        output:
        path "${cohort_name}_gdb"
    ```

We're almost done!

### 3.5. Update the process call to use the renamed module

We already updated the import statement in section 3.1.
Now we need to rename the process call in the workflow body and add the reference genome files as inputs, since we need to provide them to the joint genotyping tool.

=== "After"

    ```groovy title="genomics-2.nf" linenums="126"
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

    ```groovy title="genomics-2.nf" linenums="126"
    // Combine GVCFs into a GenomicsDB data store
    GATK_GENOMICSDB(
        all_gvcfs_ch,
        all_idxs_ch,
        intervals_file,
        params.cohort_name
    )
    ```

Now the process is completely wired up.

### 3.6. Add the joint VCF to the publish section

We need to publish the joint VCF outputs from the new process.
Add these lines to the `publish:` section of the workflow:

```groovy title="genomics-2.nf" linenums="145"
    joint_vcf = GATK_JOINTGENOTYPING.out.vcf
    joint_vcf_idx = GATK_JOINTGENOTYPING.out.idx
```

### 3.7. Add the joint VCF targets to the output block

Finally, add output targets for the joint VCF files.
We'll put them at the root of the results directory since this is the final output.

```groovy title="genomics-2.nf" linenums="157"
    joint_vcf {
        path '.'
    }
    joint_vcf_idx {
        path '.'
    }
```

Now everything should be completely wired up.

### 3.8. Run the workflow

Finally, we can run the modified workflow...

```bash
nextflow run genomics-2.nf -resume
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics-2.nf` [crazy_marconi] DSL2 - revision: 5da9afc841

    executor >  local (1)
    [9a/c7a873] SAMTOOLS_INDEX (2)       | 3 of 3, cached: 3 ✔
    [e4/4ed55e] GATK_HAPLOTYPECALLER (2) | 3 of 3, cached: 3 ✔
    [a6/7cc8ed] GATK_JOINTGENOTYPING     | 1 of 1 ✔
    ```

And it works!

You'll find the final output file, `family_trio.joint.vcf` (and its file index), in the results directory.

??? abstract "Directory contents (symlinks shortened)"

    ```console
    results_genomics/
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

If you're the skeptical type, you can click on the joint VCF file to open it and verify that the workflow has generated the same variant calls that you obtained by running the tools manually at the start of this section.

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

You know how to use some common operators as well as Groovy closures to control the flow of data in your workflow.

### What's next?

Celebrate your success and take a well-deserved break.

In the next part of this course, you'll learn how to add automated testing to verify that your pipeline produces correct outputs.
