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

Let's try running that now.

```bash
nextflow run genomics-3.nf -resume
```

And it works!

```console title="Output"
 N E X T F L O W   ~  version 24.10.0

Launching `genomics-3.nf` [gloomy_poincare] DSL2 - revision: 43203316e0

executor >  local (7)
[18/89dfa4] SAMTOOLS_INDEX (1)       | 3 of 3 ✔
[30/b2522b] GATK_HAPLOTYPECALLER (2) | 3 of 3 ✔
[a8/d2c189] GATK_JOINTGENOTYPING     | 1 of 1 ✔
```

Like previously, there will now be a `work` directory and a `results_genomics` directory inside your project directory.

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

Open the `main.nf` file and copy the `SAMTOOLS_INDEX` process definition into it, so you end up with something like this:

```groovy title="modules/samtools/index/main.nf" linenums="1"
#!/usr/bin/env nextflow

/*
 * Generate BAM index file
 */
process SAMTOOLS_INDEX {

    container 'community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464'

    publishDir params.outdir, mode: 'symlink'

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

_Before:_

```groovy title="tests/main.nf.test" linenums="1" hl_lines="1"
/*
 * Call variants with GATK HaplotypeCaller
 */
process GATK_HAPLOTYPECALLER {
```

_After:_

```groovy title="genomics-3.nf" linenums="1" hl_lines="1 2"
// Include modules
include { SAMTOOLS_INDEX } from './modules/samtools/index/main.nf'

/*
 * Call variants with GATK HaplotypeCaller
 */
process GATK_HAPLOTYPECALLER {
```

You can now run the workflow again, and it should still work the same way as before. If you supply the `-resume` flag, no new tasks should even need to be run:

```bash
nextflow run genomics-3.nf -resume
```

```console title="Re-used Output after moving SAMTOOLS_INDEX to a module"
 N E X T F L O W   ~  version 24.10.0

Launching `genomics-3.nf` [ridiculous_jones] DSL2 - revision: c5a13e17a1

[cf/289c2d] SAMTOOLS_INDEX (2)       | 3 of 3, cached: 3 ✔
[30/b2522b] GATK_HAPLOTYPECALLER (1) | 3 of 3, cached: 3 ✔
[a8/d2c189] GATK_JOINTGENOTYPING     | 1 of 1, cached: 1 ✔
```

### 1.2. Create a modules for the `GATK_HAPLOTYPECALLER` and `GATK_JOINTGENOTYPING` processes

Repeat the same steps for the remaining processes. You'll need to create a directory for each process, and then create a `main.nf` file inside that directory, removing the process definition from the workflow's `main.nf` file and adding an import declaration for the module. Once you're done, check that your modules directory structure is correct by running:

```bash
tree modules/
```

```console title="Directory structure"
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

Test the modularised workflow.

---

## 2. Test the modularised workflow

Let's try running that now.

```bash
nextflow run genomics-3.nf -resume
```

And it works!

```console title="Output"
 N E X T F L O W   ~  version 24.10.0

Launching `genomics-3.nf` [gloomy_poincare] DSL2 - revision: 43203316e0

executor >  local (7)
[18/89dfa4] SAMTOOLS_INDEX (1)       | 3 of 3 ✔
[30/b2522b] GATK_HAPLOTYPECALLER (2) | 3 of 3 ✔
[a8/d2c189] GATK_JOINTGENOTYPING     | 1 of 1 ✔
```

Yep, everything still works, including the resumability of the pipeline.

### Takeaway

You've practiced modularizing a workflow, and you've seen that it still works the same way as before.

---

## 3. Summary

So, once again (assuming you followed [Hello Modules](../../../hello_nextflow/hello_modules.md)), you've done all this work and absolutely nothing has changed to how the pipeline works! This is a good thing, because it means that you've modularised your workflow without impacting its function. Importantly, you've laid a foundation for doing things that will make your code more modular and easier to maintain- for example, you can now add tests to your pipeline using the nf-test framework. This is what we'll be looking at in the next part of this course.
