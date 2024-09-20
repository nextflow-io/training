# Part 5: Hello Config

In this section, we will explore how to configure Nextflow pipelines using configuration files, profiles, process directives, and executors. Configuration management is an essential aspect of Nextflow pipeline development, allowing you to customize the behavior of your pipeline, adapt it to different environments, and optimize resource usage. By understanding and utilizing these configuration options effectively, you can enhance the flexibility, scalability, and performance of your pipelines.

## 1. Check and modify configuration

### 1.1. Run nf-hello-gatk with default settings

```bash
nextflow run seqeralabs/nf-hello-gatk
```

When you run the pipeline with the default settings using the command above, the following happens:

1. Nextflow downloads the pipeline from the GitHub repository `seqeralabs/nf-hello-gatk`.
2. It then executes the pipeline using the default configuration.
3. The pipeline will likely use Docker containers to run the required tools (Samtools and GATK).
4. It processes the input BAM files, creating index files and performing variant calling.
5. The results are generated in the `results` directory by default.
6. Nextflow also creates a `work` directory containing intermediate files and logs.
7. Upon completion, Nextflow displays a summary of the run, including any errors or warnings.

Now let's move on to seeing how this was configured and set up.

### 1.2. Check configuration

Open the `nextflow.config` file and inspect the contents:

```bash
code nextflow.config
```

The contents should look like this:

```groovy title="nextflow.config" linenums="1"
docker.enabled = true
```

This means that the pipeline will use Docker containers to run the required tools.

### 1.3. Modify configuration

Let's modify the configuration to use Conda instead of Docker.

_Before:_

```groovy title="nextflow.config" linenums="1"
docker.enabled = true
```

_After:_

```groovy title="nextflow.config" linenums="1"
docker.enabled = false
conda.enabled = true
```

Now let's run the pipeline again with the modified configuration:

```bash title="nextflow.config"
nextflow run seqeralabs/nf-hello-gatk
```

This time, the pipeline will use Conda environments to run the required tools.

## 2. Profiles

Profiles are a way to customize the behavior of Nextflow pipelines by selection, rather than setting them permanently.

### 2.1. Create a profile

_Before:_

```groovy title="nextflow.config" linenums="1"
docker.enabled = false
conda.enabled = true
```

_After:_

```groovy title="nextflow.config" linenums="1"
profiles {
    'docker' {
        docker.enabled = true
    }
    'conda' {
        conda.enabled = true
    }
}
```

### 2.2. Run the pipeline with a profile

```bash
nextflow run seqeralabs/nf-hello-gatk -profile docker
```

or

```bash
nextflow run seqeralabs/nf-hello-gatk -profile conda
```

By creating and using profiles as demonstrated above, we've enhanced our pipeline's flexibility and ease of use. Profiles allow us to switch between Docker and Conda environments with a simple command-line option, without modifying the main configuration file. This approach streamlines the process of adapting our pipeline to different execution environments or tool preferences. We can now run our pipeline with either Docker or Conda by simply specifying the appropriate profile (`-profile docker` or `-profile conda`), making it more versatile and user-friendly. This method of configuration management improves the portability and maintainability of our Nextflow pipeline, enabling us to easily accommodate various execution scenarios.

## 3. Process directives and resources

### 3.1. Process directives

Previously, we have seen the use of process directives to modify the behaviour of a process when we added the `publishDir` directive to export files from the working directory. Let's look into directives in more detail.

### 3.1.1 Set process resources

By default, Nextflow will use a single CPU and 2GB of memory for each process. We can modify this behaviour by setting the `cpu` and `memory` directives in the `process` block. Add the following to the end of your `nextflow.config` file:

```groovy title="nextflow.config" linenums="11"
process {
    cpus = 8
    memory = 4GB
}
```

Run the pipeline again with the modified configuration:

```bash
nextflow run seqeralabs/nf-hello-gatk -profile docker
```

You shouldn't see any difference, however you might notice the three processes get bottlenecked behind each other. This is because Nextflow will make sure we aren't using more than CPUs than are available.

### 3.1.2 Modify process resources for a specific process

We can also modify the resources for a specific process using the `withName` directive. Add the following to the end of your `nextflow.config` file:

```groovy title="nextflow.config" linenums="11"
process {
    withName: 'GATK_HAPLOTYPECALLER' {
        cpus = 8
        memory = 4.GB
    }
}
```

Run the pipeline again with the modified configuration:

```bash
nextflow run seqeralabs/nf-hello-gatk -profile docker
```

Now the settings are only applied to the GATK HaplotypeCaller process. This is useful when your processes have different resource requirements, so you can right-size your resources for each process.

## 4. Executor

### 4.1. Local executor

Up until now, we have been running our pipeline with the local executor. This runs each step on the same machine that Nextflow is running on. However, for large genomics pipelines, you will want to use a distributed executor. Nextflow supports a number of different distributed executors, including:

-   HPC (SLURM, PBS, SGE)
-   AWS Batch
-   Google Batch
-   Azure Batch
-   Kubernetes

We can modify the executor used by nextflow using the `executor` process directive. Until now, we have been using the `local` executor (default). The following configuration is implied:

```groovy title="nextflow.config" linenums="18"
process {
    executor = 'local'
}
```

### 4.2. Other executors

!!! note

    This is a demonstration and designed to go wrong!

If we wish to change executor, we could simply set this to one of the values in the documentation:

```groovy title="nextflow.config" linenums="18"
process {
    executor = 'slurm'
}
```

However, if we add this to our config and run the pipeline we will that includes this error:

```console
Cannot run program "sbatch"
```

Nextflow has interpreted that we wish to submit to a Slurm cluster, which requires te use of the command `sbatch`. However, because we are not _actually_ using Slurm the `sbatch` command is not installed therefore there is an error.

If we check inside the `.command.run` file created in the work directory, we can see that Nextflow has created a script to submit the job to Slurm.

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

If our process had more directives, such as `clusterOptions`, `cpus`, `memory`, `queue`, and `time`, these would also be included in the `.command.run` file and directly passed to the Slurm execution. They would also be translated to the equivalent options for other executors. This is how Nextflow creates the commands required to correctly submit a job to the sbatch cluster via a single configuration change.

### 4.3. Using Executors in Profiles

Let's combine `profiles` with `executors`. Add the following to your configuration file:

before:

```groovy title="nextflow.config" linenums="1"
profiles {
    'docker' {
        docker.enabled = true
        conda.enabled = false
    }
    'conda' {
        docker.enabled = false
        conda.enabled = true
    }
}
```

After:

```groovy title="nextflow.config"
profiles {
    'docker' {
        docker.enabled = true
        conda.enabled = false
    }
    'conda' {
        docker.enabled = false
        conda.enabled = true
    }
    'local' {
        process.executor = 'local'
    }
    'slurm' {
        process.executor = 'slurm'
    }
}
```

Now run the pipeline using two profiles, `docker` and `local`:

```bash
nextflow run seqeralabs/nf-hello-gatk -profile docker,local
```

Note, we have returned back to the original configuration of using Docker containers with a local execution, but we can now use profiles to switch to a different software packaging system (conda) or a different executor (slurm) with a single command-line option.
