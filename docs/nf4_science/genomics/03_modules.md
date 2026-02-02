# Part 3: Moving code into modules

In the first part of this course, you built a variant calling pipeline that was completely linear and processed each sample's data independently of the others.

In the second part, we showed you how to use channels and channel operators to implement joint variant calling with GATK, building on the pipeline from Part 1.

In this part, we'll show you how to convert the code in that workflow into modules. To follow this part of the training, you should have completed Part 1 and Part 2, as well as [Hello Modules](../../../hello_nextflow/hello_modules.md), which covers the basics of modules.

---

## 0. Warmup

When we started developing our workflow, we put everything in one single code file.
Now it's time to tackle **modularizing** our code, _i.e._ extracting the process definitions into modules.

We're going to start with the same workflow as in Part 2, which we've provided for you in the file `genomics-3.nf`.

!!! note

     Make sure you're in the correct working directory:
     `cd /workspaces/training/nf4-science/genomics`

Run the workflow to verify the starting point:

```bash
nextflow run genomics-3.nf -resume
```

```console title="Output"
 N E X T F L O W   ~  version 25.10.2

Launching `genomics-3.nf` [serene_borg] DSL2 - revision: 0cbebb67a1

executor >  local (7)
[6f/83ee72] SAMTOOLS_INDEX (3)       | 3 of 3 ✔
[53/b9d342] GATK_HAPLOTYPECALLER (1) | 3 of 3 ✔
[0c/fa6d15] GATK_JOINTGENOTYPING     | 1 of 1 ✔
```

There will now be a `work` directory and a `results_genomics` directory inside your project directory.

### Takeaway

You're ready to start modularizing your workflow.

### What's next?

Move the Genomics workflow's processes into modules.

---

## 1. Move processes into modules

As you learned in [Hello Modules](../../../hello_nextflow/hello_modules.md), you can create a module simply by copying the process definition into its own file, in any directory, and you can name that file anything you want.

For reasons that will become clear later (in particular when we come to testing), in this training we'll follow the convention of naming the file `main.nf`, and placing it in a directory structure named after the tool kit and the command.

### 1.1. Create a module for the `SAMTOOLS_INDEX` process

In the case of the `SAMTOOLS_INDEX` process, 'samtools' is the toolkit and 'index' is the command. So, we'll create a directory structure `modules/samtools/index` and put the `SAMTOOLS_INDEX` process definition in the `main.nf` file inside that directory.

```bash
mkdir -p modules/samtools/index
touch modules/samtools/index/main.nf
```

Open the `main.nf` file and copy the `SAMTOOLS_INDEX` process definition into it.

```groovy title="modules/samtools/index/main.nf" linenums="1"
/*
 * Generate BAM index file
 */
process SAMTOOLS_INDEX {

    container 'community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464'

    input:
    path input_bam

    output:
    tuple path(input_bam), path("${input_bam}.bai")

    script:
    """
    samtools index '$input_bam'
    """
}
```

Then, remove the `SAMTOOLS_INDEX` process definition from `genomics-3.nf`, and add an import declaration for the module before the next process definition, like this:

=== "After"

    ```groovy title="genomics-3.nf" linenums="1" hl_lines="1 2"
    // Include modules
    include { SAMTOOLS_INDEX } from './modules/samtools/index/main.nf'

    /*
     * Call variants with GATK HaplotypeCaller
     */
    process GATK_HAPLOTYPECALLER {
    ```

=== "Before"

    ```groovy title="genomics-3.nf" linenums="1" hl_lines="1"
    /*
     * Call variants with GATK HaplotypeCaller
     */
    process GATK_HAPLOTYPECALLER {
    ```

You can now run the workflow again, and it should still work the same way as before. If you supply the `-resume` flag, no new tasks should even need to be run:

```bash
nextflow run genomics-3.nf -resume
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `genomics-3.nf` [sleepy_snyder] DSL2 - revision: aa68d06c43

    [0f/71b55e] SAMTOOLS_INDEX (1)       | 3 of 3, cached: 3 ✔
    [f1/18971b] GATK_HAPLOTYPECALLER (3) | 3 of 3, cached: 3 ✔
    [0c/fa6d15] GATK_JOINTGENOTYPING     | 1 of 1, cached: 1 ✔
    ```

### 1.2. Create modules for the `GATK_HAPLOTYPECALLER` and `GATK_JOINTGENOTYPING` processes

Repeat the same steps for the remaining processes.
For each process:

1. Create the directory structure (`modules/gatk/haplotypecaller/` and `modules/gatk/jointgenotyping/`)
2. Create a `main.nf` file containing the process definition
3. Remove the process definition from `genomics-3.nf`
4. Add an import declaration for the module

Once you're done, check that your modules directory structure is correct by running:

```bash
tree modules/
```

??? abstract "Directory contents"

    ```console
    modules/
    ├── gatk
    │   ├── haplotypecaller
    │   │   └── main.nf
    │   └── jointgenotyping
    │       └── main.nf
    └── samtools
        └── index
            └── main.nf

    5 directories, 3 files
    ```

You should also have something like this in the main workflow file, after the parameters section:

```
include { SAMTOOLS_INDEX } from './modules/samtools/index/main.nf'
include { GATK_HAPLOTYPECALLER } from './modules/gatk/haplotypecaller/main.nf'
include { GATK_JOINTGENOTYPING } from './modules/gatk/jointgenotyping/main.nf'

workflow {
```

### Takeaway

You've practiced modularizing a workflow, with the genomics workflow as an example.

### What's next?

Test the modularized workflow.

---

## 2. Test the modularized workflow

Run the modularized workflow to verify everything still works.

```bash
nextflow run genomics-3.nf -resume
```

```console title="Output"
 N E X T F L O W   ~  version 25.10.2

Launching `genomics-3.nf` [astonishing_venter] DSL2 - revision: ca27264c13

[6f/83ee72] SAMTOOLS_INDEX (3)       | 3 of 3, cached: 3 ✔
[53/b9d342] GATK_HAPLOTYPECALLER (3) | 3 of 3, cached: 3 ✔
[0c/fa6d15] GATK_JOINTGENOTYPING     | 1 of 1, cached: 1 ✔
```

Everything still works, including the resumability of the pipeline.
The results continue to be published to the `results_genomics` directory.

```console title="Directory contents"
results_genomics/
├── family_trio.joint.vcf
├── family_trio.joint.vcf.idx
├── gvcf
│   ├── reads_father.bam.g.vcf
│   ├── reads_father.bam.g.vcf.idx
│   ├── reads_mother.bam.g.vcf
│   ├── reads_mother.bam.g.vcf.idx
│   ├── reads_son.bam.g.vcf
│   └── reads_son.bam.g.vcf.idx
└── indexed_bam
    ├── reads_father.bam
    ├── reads_father.bam.bai
    ├── reads_mother.bam
    ├── reads_mother.bam.bai
    ├── reads_son.bam
    └── reads_son.bam.bai
```

### Takeaway

You've modularized a workflow and verified it still works the same way as before.

### What's next?

Review what you've learned and look ahead to testing.

---

## 3. Summary

You've modularized the workflow, and nothing has changed to how the pipeline works.
This is intentional: you've restructured the code without impacting its function.

The modules contain only the process logic, making them clean and reusable.
The main script controls what gets published and where, while modules remain focused on their computational task.

You've laid a foundation for things that will make your code easier to maintain.
For example, you can now add tests to your pipeline using the nf-test framework.
This is what we'll look at in the next part of this course.
