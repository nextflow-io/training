# Part 5: Hello Config

This section will explore how to configure Nextflow pipelines using configuration files, profiles, process directives, executors, and parameter files.
Configuration management is an essential aspect of Nextflow pipeline development, allowing you to customize the behavior of your pipeline, adapt it to different environments, and optimize resource usage.
By understanding and effectively utilizing these configuration options, you can enhance the flexibility, scalability, and performance of your pipelines.

TODO (refine)-> Learn to modify the execution of a pipeline _without altering a single line of code_. This highlights the power of Nextflow's configuration; enabling you to control how the pipeline runs without changing what it runs. Use this flexibility to adapt your pipeline to run in any environment.

## 0. Warmup: Moving to a formal project structure

So far we've been working with a very loose structure, with just one workflow code file and a tiny configuration file that we've mostly ignored, because we were very focused on learning how to implement the workflow itself.
However, we're now moving into the phase of this training series that is more focused on code development and maintenance practices.

As part of that, we're going to adopt a formal project structure. We're going to work inside a dedicated project directory called `projectC` (C for configuration), and we've renamed the workflow file `main.nf` to match the recommended Nextflow convention.

```console title="Directory contents"
projectC/
├── data -> ../data
├── main.nf
└── nextflow.config
```

The `main.nf` file is a copy of the workflow produced by completing Part 4 of this training course.

The `nextflow.config` file is a copy of the original `nextflow.config` file from the `hello-nextflow` directory (which we've ben using so far). It is the default configuration file that is expected by Nextflow.

```console title="nextflow.config" linenums="1"
docker.fixOwnership = true
docker.enabled = true
```

The `docker.fixOwnership = true` line is not really interesting.
It's a workaround for an issue that sometimes occur with containerized tools that set the wrong permissions on the files they write (which is the case with GATK GenomicsDBImport in our workflow).

The `docker.enabled = true` line is what we care about here.
It specifies that Nextflow should use Docker containers to execute process calls.
We're going to be playing with that shortly.

We've also created a symbolic link called `data` pointing to the data directory, to avoid having to change anything to how the file paths are set up.
Later we'll cover a better way of handling this, but this works for now.

You can test that it runs properly:

```bash
cd projectC
nextflow run main.nf
```

This should run successfully:

```console title="Output"
Nextflow 24.09.2-edge is available - Please consider updating your version to it

 N E X T F L O W   ~  version 24.02.0-edge

 ┃ Launching `main.nf` [tender_brahmagupta] DSL2 - revision: 848ff2f9b5

executor >  local (7)
[fb/f755b1] SAMTOOLS_INDEX (1)       [100%] 3 of 3 ✔
[d8/467767] GATK_HAPLOTYPECALLER (1) [100%] 3 of 3 ✔
[ee/2c7855] GATK_JOINTGENOTYPING     [100%] 1 of 1 ✔
```

There will now be a `work` directory and a `results_genomics` directory inside your current `projectC` directory.

### Takeaway

You know what are the two most important files in a Nextflow project: `main.nf` and its `nextflow.config`.

### What's next?

Learn how to modify basic configuration properties to adapt to your compute environment's requirements.

---

## 1. Configuration basics

In the very first part of this training course (Part 1: Hello World) we just used locally installed software in our workflow. Then from Part 2 onward, we've been using Docker containers.

Now, let's pretend we're working on an HPC cluster and the admin doesn't allow the use of Docker for security reasons.

### 1.1 Disable Docker in the config file

First, we have to switch the value of `docker.enabled` to false.

_Before:_

```console title="nextflow.config" linenums="1"
docker.fixOwnership = true
docker.enabled = true
```

_After:_

```console title="nextflow.config" linenums="1"
docker.fixOwnership = true
docker.enabled = false
```

Let's see what happens if we run that.

### 1.2 Run the workflow without Docker

We are now launching the `main.nf` workflow from inside the `projectC` directory.

```bash
nextflow run main.nf
```

As expected, the run fails with an error message that looks like this:

```console title="Output"
 N E X T F L O W   ~  version 24.02.0-edge

 ┃ Launching `projectC/main.nf` [silly_ramanujan] DSL2 - revision: 9129bc4618

executor >  local (3)
[93/4417d0] SAMTOOLS_INDEX (1)   [  0%] 0 of 3
[-        ] GATK_HAPLOTYPECALLER -
[-        ] GATK_JOINTGENOTYPING -
ERROR ~ Error executing process > 'SAMTOOLS_INDEX (2)'

Caused by:
  Process `SAMTOOLS_INDEX (2)` terminated with an error exit status (127)

Command executed:

  samtools index 'reads_father.bam'

Command exit status:
  127

Command output:
  (empty)

Command error:
  .command.sh: line 2: samtools: command not found
```

Command not found? Of course, we don't have Samtools installed in our environment, and we can no longer use the Docker container. What to do?

!!!note

    Nextflow supports multiple other container technologies such as including Singularity (which is more widely used on HPC), and software package managers such as Conda.

Let's try using Conda environments for our workflow.

### 1.3 Add a conda environment to the Samtools process definition

We know that Bioconda provides a Samtools environment, so we just need to retrieve its URI and add it to the process definition using the `conda` directive.

_Before:_

```console title="main.nf" linenums="22"
process SAMTOOLS_INDEX {

    container 'community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464'

    publishDir 'results_genomics', mode: 'symlink'
```

_After:_

```console title="main.nf" linenums="22"
process SAMTOOLS_INDEX {

    container "community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464"
    conda "bioconda::samtools=1.20"

    publishDir 'results_genomics', mode: 'symlink'
```

Make sure to _add_ the conda directive.
We're not _replacing_ the docker directive, just adding an alternative option.

### 1.4 Enable Conda in the configuration file

We need to add a line enabling the conda directive.
And while we're add it, let's put a blank line before thos two to emphasize the logical grouping.

_Before:_

```groovy title="nextflow.config" linenums="1"
docker.fixOwnership = true
docker.enabled = false
```

_After:_

```groovy title="nextflow.config" linenums="1"
docker.fixOwnership = true

docker.enabled = false
conda.enabled = true
```

This should allow Nextflow to use the Conda environment we added for Samtools.

### 1.5 Run it to see if it works

Let's try it out.

```bash
nextflow run main.nf
```

This will take a bit longer than usual the first time, and you might see the console output stay 'stuck' at this stage for a minute or so:

```console title="Output"
 N E X T F L O W   ~  version 24.02.0-edge

 ┃ Launching `main.nf` [extravagant_thompson] DSL2 - revision: 848ff2f9b5

[-        ] SAMTOOLS_INDEX       -
[-        ] GATK_HAPLOTYPECALLER -
[-        ] GATK_JOINTGENOTYPING -
Creating env using conda: bioconda::gatk4=4.5.0.0 [cache /workspace/gitpod/hello-nextflow/projectC/work/conda/env-6684ea23d69ceb1742019ff36904f612]
```

That's because Nextflow has to retrieve and spin up the Conda environment, which takes a bit of work behind the scenes.
Then after a moment it'll spit out the rest of the output. It may look a little messy:

```console title="Output"
  N E X T F L O W   ~  version 24.02.0-edge

 ┃ Launching `main.nf` [big_engelbart] DSL2 - revision: 63efa30da7

[-        ] SAMTOOLS_INDEX       -
executor >  local (3)
executor >  local (6)
[6e/78985e] SAMTOOLS_INDEX (1)       [100%] 3 of 3 ✔
[eb/70a150] GATK_HAPLOTYPECALLER (3) [  0%] 0 of 3
[-        ] GATK_JOINTGENOTYPING     -
Creating env using conda: bioconda::samtools=1.20 [cache /workspace/gitpod/hello-nextflow/projectC/work/conda/env-1f9e4747cd58cb43e7ca4da34bb66eee]
ERROR ~ Error executing process > 'GATK_HAPLOTYPECALLER (1)'

Caused by:
  Process `GATK_HAPLOTYPECALLER (1)` terminated with an error exit status (127)

Command executed:

  gatk HaplotypeCaller         -R ref.fasta         -I reads_son.bam         -O reads_son.bam.g.vcf         -L intervals.bed         -ERC GVCF

Command exit status:
  127

Command output:
  (empty)

Command error:
  .command.sh: line 2: gatk: command not found
```

But this is great!
We see that the Samtools process calls were executed successfully.
The overall run just fails because we haven't yet added a Conda environment for the GATK processes, so let's do that now.

### 1.6 Add conda environment to GATK processes

There are two GATK processes that need to be updated, GATK_HAPLOTYPECALLER and GATK_JOINTGENOTYPING.

#### 1.6.1. Update GATK_HAPLOTYPECALLER

_Before:_

```console title="main.nf" linenums="43"
process GATK_HAPLOTYPECALLER {

    container "community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867"

    publishDir 'results_genomics', mode: 'symlink'
```

_After:_

```console title="main.nf" linenums="43"
process GATK_HAPLOTYPECALLER {

    container "community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867"
    conda "bioconda::gatk4=4.5.0.0"

    publishDir 'results_genomics', mode: 'symlink'
```

#### 1.6.2. Update GATK_JOINTGENOTYPING

_Before:_

```console title="main.nf" linenums="74"
process GATK_JOINTGENOTYPING {

    container "community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867"

    publishDir 'results_genomics', mode: 'symlink'
```

_After:_

```console title="main.nf" linenums="74"
process GATK_JOINTGENOTYPING {

    container "community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867"
    conda "bioconda::gatk4=4.5.0.0"

    publishDir 'results_genomics', mode: 'symlink'
```

Once both processes are updated, we can try it out again.

### 1.7 Run the workflow again

This time it should all work.

```bash
nextflow run main.nf
```

And so it does!

```console title="Output"
 N E X T F L O W   ~  version 24.02.0-edge

 ┃ Launching `main.nf` [silly_goldstine] DSL2 - revision: a60f9fd6af

executor >  local (7)
[23/b59106] SAMTOOLS_INDEX (1)       [100%] 3 of 3 ✔
[da/e1bf1d] GATK_HAPLOTYPECALLER (1) [100%] 3 of 3 ✔
[2e/e6ffca] GATK_JOINTGENOTYPING     [100%] 1 of 1 ✔
```

This means we're all set to run with Conda environments if need.

### Takeaway

You know how to switch software packaging systems using a configuration file.

### What's next?

Learn how to use profiles to make selecting an option easier.

---

## 2. Profiles

Profiles are a way to customize the behavior by selecting an option at runtime, to avoid having to edit a file every time you want to do something differently.

### 2.1. Create a profile

We're still working inside the `nextflow.config` file, but we're going to restructure how we write the configuration.

_Before:_

```groovy title="nextflow.config" linenums="1"
docker.fixOwnership = true

docker.enabled = false
conda.enabled = true
```

_After:_

```groovy title="nextflow.config" linenums="1"
docker.fixOwnership = true

profiles {
    docker {
        docker.enabled = true
    }
    conda {
        conda.enabled = true
    }
}
```

This

### 2.2. Run the pipeline with a profile

```bash
nextflow run seqeralabs/nf-hello-gatk -r main -profile docker
```

or

```bash
nextflow run seqeralabs/nf-hello-gatk -r main -profile conda
```

As demonstrated above, by creating and using profiles, we've enhanced our pipeline's flexibility and ease of use.
We can now run our pipeline with Docker or Conda using a single command line argument by specifying the appropriate profile (`-profile docker` or `-profile conda`).
This method of configuration management improves the portability and maintainability of our Nextflow pipeline, enabling us to accommodate various execution scenarios easily.

### Takeaway

You know how to use profiles to customize the configuration of your pipeline.

### What's next?

Learn how to change the executor used by Nextflow.

---

## 3. Executors

Until now, we have been running our pipeline with the local executor.
This runs each step on the same machine that Nextflow is running on.
However, for large genomics pipelines, you will want to use a distributed executor.
Nextflow supports several different distributed executors, including:

-   HPC (SLURM, PBS, SGE)
-   AWS Batch
-   Google Batch
-   Azure Batch
-   Kubernetes

The executor is subject to process directive called `executor`. By default it is set to `local`, so the following configuration is implied:

```groovy title="Built-in configuration"
process {
    executor = 'local'
}
```

### 3.1. Switch to a different executor

!!! note

    This is for demonstration purposes but will not work since we don't have access to an external executor!

TODO a common HPC executor is slurm

```groovy title="nextflow.config" linenums="1"
process {
    executor = 'slurm'
}
```

Assumes slurm is installed, connection to a cluster; generates an sbatch command

### 3.2. Run to generate the job submission script

Run it so we can see what nextflow does

TODO run command

Generates an error which is expected since we don't have slurm

```console
Cannot run program "sbatch"
```

TODO but if we check inside the `.command.run` file created in the work directory, we can see that Nextflow has created a script to submit the job to Slurm.

```bash title=".command.run" linenums="1"
#!/bin/bash
#SBATCH -J nf-SAMTOOLS_INDEX_(1)
#SBATCH -o /home/gitpod/work/34/850fe31af0eb62a0eb1643ed77b84f/.command.log
#SBATCH --no-requeue
#SBATCH --signal B:USR2@30
NXF_CHDIR=/home/gitpod/work/34/850fe31af0eb62a0eb1643ed77b84f
### ---
### name: 'SAMTOOLS_INDEX (1)'
### container: 'community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464'
### outputs:
### - 'reads_father.bam'
### - 'reads_father.bam.bai'
### ...
```

Can control other options using other process directives such as `clusterOptions`, `cpus`, `memory`, `queue`, and `time`, these would also be included in the `.command.run` file and directly passed to the Slurm execution.
TODO: demonstrate this later in module.

If using a different executor, values translated into equivalent options.

This is how Nextflow creates the commands required to correctly submit a job to the sbatch cluster via a single configuration change.

### 3.3. Set the executor in a profile

Let's combine profiles with executors. Make the following changes to your configuration file:

Remove the following lines:

```groovy title="nextflow.config" linenums="17"
process {
    executor = 'slurm'
}
```

Now create profiles for local and slurm executors

Before:

```groovy title="nextflow.config" linenums="1"
profiles {
    docker {
        docker.enabled = true
        conda.enabled = false
    }
    conda {
        docker.enabled = false
        conda.enabled = true
    }
}
```

After:

```groovy title="nextflow.config"
profiles {
    docker {
        docker.enabled = true
        conda.enabled = false
    }
    conda {
        docker.enabled = false
        conda.enabled = true
    }
    local {
        process.executor = 'local'
    }
    slurm {
        process.executor = 'slurm'
    }
}
```

### 3.4. Run with a combination of profiles

Now run the pipeline using two profiles, `docker` and `local`:

```bash
nextflow run seqeralabs/nf-hello-gatk -profile docker,local
```

We have returned to the original configuration of using Docker containers with local execution.
However, now we can use profiles to switch to a different software packaging system (conda) or a different executor (slurm) with a single command-line option. Eg could do `-profile conda,slurm` if running on an institutional HPC.

### Takeaway

You now know how to change the executor and combine that with other environment settings.

### What's next?

Learn how to change process resource allocations.

---

## 4. Allocate compute resources with process directives

In a previous training module, we used process directives to modify the behavior of a process when we added the `publishDir` directive to export files from the working directory.
Let's look into directives in more detail.

### 4.1. Set default process resources

By default, Nextflow will use a single CPU and 2GB of memory for each process.
We can modify this behavior by setting the `cpu` and `memory` directives in the `process` block.
Add the following to the end of your `nextflow.config` file:

```groovy title="nextflow.config" linenums="11"
process {
    cpus = 4
    memory = 2.GB
}
```

Run the pipeline again with the modified configuration:

```bash
nextflow run seqeralabs/nf-hello-gatk -r main -profile docker
```

You may not notice any difference in terms of how quickly this runs, especially if you have a powerful machine.
But if you have a lightweight machine and allocate high numbers, you might see process calls getting queued behind each other.
This is because Nextflow will ensure we aren't using more CPUs than are available.

!!! tip

    You can check the number of CPUs allocated to a given process by looking at the `.command.run` log in its work directory.
    There will be a function called `nxf_launch()` that includes the command `docker run -i --cpu-shares 1024`, where `--cpu-shares` is the number of CPUs given to the process multiplied by 1024.

TODO transition: how do you know what to give/ when to tweak? first step is finding what the tools are actually using. good news: can generate a report of resource utilization

### 4.2. Generate a report on the pipeline's resource utilization

TODO run this

nextflow run hello-config.nf -with-report report-config-1.html

TODO: output is an html file; open in browser

TODO: screenshot?

TODO: based on findings, might want to tweak resources for specific process

### 4.3. Adjust resource allocations for a specific process

We can allocate the resources for a specific process using the `withName` directive.
Add the following to the end of your `nextflow.config` file:

```groovy title="nextflow.config" linenums="11"
process {
    withName: 'GATK_HAPLOTYPECALLER' {
        cpus = 8
        memory = 4.GB
    }
}
```

TODO This should do X

### 4.4. Run again with the modified configuration

TODO Run the pipeline again with the modified configuration and get report

```bash
nextflow
```

Now, the default settings apply to everything but the GATK HaplotypeCaller process gets special treatment.

TODO: open the report and look at the differences.

This is useful when your processes have different resource requirements; you can right-size your resources for each process.

### 4.5. Add a custom conf file to a profile

TODO: If e.g. your cluster has more generous resources, can set up a custom conf file specifying higher process allocations. Make a new profile with include statement

TODO code (+ make a conf file)

### 4.6. Run with profile conda,slurm,beefy

TODO: scenario is you have access to heavyweight nodes on HPC

TODO run with profile conda,slurm,beefy

TODO: check inside the `.command.run` file created in the work directory, look at values assigned

TODO: wrap up; this is how we can have alternate confs to override settings in various contexts

### Takeaway

You know how to allocate process resources, tweak based on utilization, and adapt based on compute environment.

### What's next?

Configuring pipeline parameters (on the scientific side)

---

## 5. Configure pipeline parameters

TODO introduce the idea of also modifying algorithm settings without modifying workflow code

### 5.1. Move params to nextflow.conf

TODO literally just copy over

TODO show code

### 5.2. Run and show is all same

TODO run with docker profile and -resume

all is same

### 5.3. Change to params{} syntax

TODO explain is nicer

### 5.4. Run and show is all still same

TODO run with docker profile and -resume

all is still same

### 5.5. Run with a params-file to override defaults

TODO explain want to have presets just like for infra conf

TODO show example params (still same files for now but in future could point to cloud files maybe?)

TODO run with params file

### 5.6. Remove default files

TODO now that we have a better way of providing input files, would be better to remove values from defaults

TODO remove file paths but leave a sensible cohort name default (any other non file values?)

TODO That's much cleaner

### 5.7. Create demo profile

TODO explain want a runnable demo example without having to point to params-file

TODO create demo profile with include statement loading demo params file

### 5.8. Run with demo profile

TODO Now we can be lazy again

TODO run with profile docker,demo

TODO explain you can have a bunch of other configuration files. Point to [order of precedence](https://www.nextflow.io/docs/latest/config.html#configuration-file) in docs.

Will focus on nextflow.config (and directly related) in this module.

### Takeaway

You know how to manage parameter defaults, inputs and set up profiles.

### What's next?

Celebrate and relax. Then move on to learning how to modularize the code for reuse.
