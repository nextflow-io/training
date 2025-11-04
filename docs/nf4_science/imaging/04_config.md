# Part 4: Configuration

In Parts 1-3, we learned how to run Nextflow, run an nf-core pipeline, and manage inputs with parameter files and samplesheets.
Now we'll explore how to configure pipelines for different computing environments using **configuration files** and **profiles**.

## Learning objectives

By the end of this part, you'll be able to:

- Understand how Nextflow resolves configuration from multiple sources
- Use nf-core built-in profiles for containers and testing
- Create custom profiles for different computing environments
- Customize resource requests using process labels
- Manage resource limits in constrained environments
- Inspect resolved configuration with `nextflow config`

---

## 1. Understanding Nextflow configuration

### 1.1. What is a configuration file?

Nextflow uses configuration files to separate **workflow logic** (what to do) from **execution settings** (how and where to do it).

Configuration files control:

- Container engines (Docker, Singularity, Conda)
- Compute resources (CPUs, memory, time)
- Execution platforms (local, HPC, cloud)
- Pipeline parameters

### 1.2. Configuration precedence

Nextflow loads configuration from multiple sources, with later sources overriding earlier ones:

1. **Pipeline config**: `nextflow.config` in the pipeline repository
2. **Directory config**: `nextflow.config` in your present working directory
3. **User config**: `~/.nextflow/config`
4. **Command-line**: Parameters and options passed directly

This layered approach lets you keep defaults in the pipeline, override with user-specific settings, and make quick adjustments on the command line.

### 1.3. Our current configuration

Let's look at the configuration we've been using:

```groovy title="nextflow.config"
docker.enabled = true
process {
    resourceLimits = [
        cpus: 2,
        memory: '7.GB',
    ]
}

```

Let's comment out or change back the `docker.enabled = true` line from Part 2, and figure out how we can achieve the same result using a profile in molkart instead.

---

## 2. Using profiles

### 2.1. What are profiles?

Profiles are named sets of configuration that can be activated with the `-profile` flag via the `nextflow run` command.
They make it easy to switch between different compute scenarios without editing config files.

All nf-core pipelines come with a number of default profiles we can take use.

### 2.2. Inspecting built-in profiles

Let's inspect them in the `molkart/nextflow.config` file associated with the pipeline codebase:

```bash
code molkart/nextflow.config
```

Search for the `profiles` block:

```groovy title="molkart/nextflow.config (excerpt)"
profiles {
    docker {
        docker.enabled          = true
        singularity.enabled     = false
        conda.enabled           = false
    }
    singularity {
        singularity.enabled     = true
        singularity.autoMounts  = true
        docker.enabled          = false
    }
    conda {
        conda.enabled           = true
        docker.enabled          = false
        conda.channels          = ['conda-forge', 'bioconda']
    }
}
```

Common container profiles:

- `docker`: Use Docker containers (most common for local development)
- `singularity`: Use Singularity/Apptainer (common on HPC)
- `conda`: Use Conda environments
- `apptainer`: Use Apptainer containers

### 2.3. Re-running with profiles instead of nextflow.config

Now that we've disabled the docker configuration in our local `nextflow.config` file and understand profiles, let's re-run the pipeline using the `-profile` flag.

Previously in Part 3, we created a `params.yaml` file with our custom parameters.
We can now combine that with the built-in Docker profile:

```bash
nextflow run ./molkart \
  -profile docker \
  -params-file params.yaml \
  -resume
```

Let's break down what each flag does:

- `-profile docker`: Activates the Docker profile from molkart's `nextflow.config`, which sets `docker.enabled = true`
- `-params-file params.yaml`: Loads all pipeline parameters from our YAML file
- `-resume`: Reuses cached results from previous runs

Because we're using `-resume`, Nextflow will check if anything changed since the last run.
If the parameters, inputs, and code are the same, all tasks will be retrieved from cache and the pipeline will complete almost instantly.

```console title="Output (excerpt)"
executor >  local (12)
...
[1a/2b3c4d] NFCORE_MOLKART:MOLKART:MINDAGAP_MINDAGAP (mem_only)   [100%] 2 of 2, cached: 2 ✔
[5e/6f7g8h] NFCORE_MOLKART:MOLKART:CLAHE (mem_only)               [100%] 2 of 2, cached: 2 ✔
...
-[nf-core/molkart] Pipeline completed successfully-
```

Notice all processes show `cached: 2` or `cached: 1` - nothing was re-executed!

### 2.4. Test profiles

Test profiles provide quick ways to specify default input parameters adn datafiles to let you verify the pipeline works.
nf-core pipelines will always include at least two test profiles:

- `test`: Small dataset with fast parameters for quick testing
- `test_full`: More comprehensive test with larger data

Let's take a closer look at the `test` profile in molkart which is included using the `includeConfig` directive:

```groovy title="molkart/nextflow.config (excerpt)"
profiles {
  ...
    test      { includeConfig 'conf/test.config'      }
}
```

This means whenever we run the pipeline with `-profile test`, Nextflow will load the configuration from `conf/test.config`.

```groovy title="molkart/conf/test.config (excerpt)"
params {
    config_profile_name        = 'Test profile'
    config_profile_description = 'Minimal test dataset to check pipeline function'

    input = 'https://raw.githubusercontent.com/nf-core/test-datasets/molkart/test_data/samplesheets/samplesheet_membrane.csv'
    mindagap_tilesize = 90
    mindagap_boxsize = 7
    mindagap_loopnum = 100
    clahe_pyramid_tile = 368
    segmentation_method = "mesmer,cellpose,stardist"
}

process {
    resourceLimits = [
        cpus: 4,
        memory: '15.GB',
        time: '1.h'
    ]
}
```

Notice that this profile contains the same parameters we used in our `params.yaml` file earlier.

You can activate multiple profiles by separating them with commas.
Let's use that to test our pipeline without needing our params file:

```bash
nextflow run ./molkart -profile docker,test --outdir results -resume
```

This combines:

- `docker`: Enable Docker containers
- `test`: Use test dataset and parameters

Profiles are applied left to right, so later profiles override earlier ones if they set the same values.

### Takeaway

nf-core pipelines come with built-in profiles for containers, testing, and special environments.
You can combine multiple profiles to build the configuration you need.

### What's next?

Learn how to create your own custom profiles for different computing environments.

---

## 3. Creating custom profiles

### 3.1. Create profiles for switching between local development and execution on HPC

Let's create custom profiles for two scenarios:

1. Local development with Docker
2. University HPC with Slurm scheduler and Singularity

Add the following to your `nextflow.config`:

```groovy title="nextflow.config"
profiles {
    local_dev {
        docker.enabled = true
        process.executor = 'local'
    }

    hpc_cluster {
        singularity.enabled = true
        process.executor = 'slurm'
        process.queue = 'standard_queue'
        singularity.cacheDir = '/shared/containers'
    }
}
```

Now you can switch between environments easily:

```bash
# For local development
nextflow run ./molkart -profile local_dev --input data/samplesheet.csv --outdir results

# For HPC (when available)
nextflow run ./molkart -profile hpc_cluster --input data/samplesheet.csv --outdir results
```

!!! Note

    We can't test the HPC profile in this training environment since we don't have access to a Slurm scheduler.
    But this shows how you would configure it for real-world use.

### 3.2. Use `nextflow config` to inspect configuration

The `nextflow config` command shows the fully resolved configuration without running the pipeline.

View the default configuration:

```bash
nextflow config ./molkart
```

View configuration with a specific profile:

```bash
nextflow config -profile local_dev ./molkart
```

This is extremely useful for:

- Debugging configuration issues
- Understanding what values will actually be used
- Checking how multiple profiles interact

### Takeaway

Custom profiles let you switch between different computing environments with a single command-line flag.
Use `nextflow config` to inspect the resolved configuration before running.

### What's next?

Learn how to customize resource requests for individual processes using nf-core's process label system.

---

## 4. Customizing resource requests

### 4.1. Understanding process labels in nf-core pipelines

For simplicity, nf-core pipelines use [**process labels**](https://www.nextflow.io/docs/latest/reference/process.html#process-label) to standardize resource allocation across all pipelines.
Each process is tagged with a label like `process_low`, `process_medium`, or `process_high` to describe low, medium, or high compute resource requirements, respectively.
These labels get converted into specific resource requests in one of the configuration files located in the `conf/` directory of the pipeline.

```groovy title="molkart/conf/base.config (excerpt)"
process {
    cpus   = { 1      * task.attempt }
    memory = { 6.GB   * task.attempt }
    time   = { 4.h    * task.attempt }

    errorStrategy = { task.exitStatus in ((130..145) + 104) ? 'retry' : 'finish' }
    maxRetries    = 1
    maxErrors     = '-1'

    withLabel:process_single {
        cpus   = { 1                   }
        memory = { 6.GB * task.attempt }
        time   = { 4.h  * task.attempt }
    }
    withLabel:process_low {
        cpus   = { 2     * task.attempt }
        memory = { 12.GB * task.attempt }
        time   = { 4.h   * task.attempt }
    }
    withLabel:process_medium {
        cpus   = { 6     * task.attempt }
        memory = { 36.GB * task.attempt }
        time   = { 8.h   * task.attempt }
    }
    withLabel:process_high {
        cpus   = { 12    * task.attempt }
        memory = { 72.GB * task.attempt }
        time   = { 16.h  * task.attempt }
    }
}
```

Notice the `task.attempt` multiplier - this allows subsequent task retries to request more resources, if the pipeline is set with `process.maxRetries > 1`.

### 4.2. Overriding resources for specific processes

For fine-grained control, target individual processes by name:

```groovy title="nextflow.config"
process {
    withName: 'NFCORE_MOLKART:MOLKART:CELLPOSE' {
        cpus   = 16
        memory = 32.GB
    }
}
```

If we try to run this pipeline with the above override, the `CELLPOSE` process will request 16 CPUs and 32 GB memory instead of the default defined by its label.
This will cause the pipeline to fail in our current environment since we don't have that much RAM available.
We'll learn how to prevent these types of failures in the next section.

!!! Tip

    To find process names, look at the pipeline execution output or check `.nextflow.log`.
    Process names follow the pattern `WORKFLOW:SUBWORKFLOW:PROCESS`.

### Takeaway

nf-core pipelines use process labels to standardize resource allocation.
You can override resources by label (affects multiple processes) or by name (affects one specific process).

### What's next?

Learn how to manage resource limits in constrained environments like GitHub Codespaces.

---

## 5. Managing resources in constrained environments

### 5.1. The resource limits problem

If we tried to run molkart with a process requesting 16 CPUs and 32 GB memory (as shown in section 4.2), it would fail in our current environment because we don't have that many resources available.
In a cluster environment with larger nodes, such requests would be submitted to the scheduler.

In constrained environments like GitHub Codespaces, without limits, Nextflow would refuse to run processes that exceed available resources.

### 5.2. Setting resource limits

The `resourceLimits` directive caps resource requests at specified values:

```groovy title="nextflow.config"
process {
    resourceLimits = [ cpus: 2, memory: 7.GB ]
}
```

This tells Nextflow: "If any process requests more than 2 CPUs or 7 GB memory, cap it at these limits instead."

### 5.3. Adding resource limits to custom profiles

Update your custom profiles to include appropriate limits:

```groovy title="nextflow.config"
profiles {
    local_dev {
        docker.enabled = true
        process.executor = 'local'
        process.resourceLimits = [
            cpus: 2,
            memory: 7.GB
        ]
    }

    hpc_cluster {
        singularity.enabled = true
        process.executor = 'slurm'
        process.queue = 'batch'
        process.resourceLimits = [
            cpus: 32,
            memory: 128.GB,
            time: 24.h
        ]
    }
}
```

!!! Warning

    Setting resource limits too low may cause processes to fail or run slowly.
    The pipeline may need to use less memory-intensive algorithms or process data in smaller chunks.

### Takeaway

Use `resourceLimits` to run pipelines in resource-constrained environments by capping process resource requests.
Different profiles can have different limits appropriate for their environment.

### What's next?

You've completed the core Nextflow for Bioimaging training!

---

## Conclusion

You now understand how to configure Nextflow pipelines for different computing environments.

Key skills you've learned:

- **Configuration precedence**: How Nextflow resolves settings from multiple sources
- **nf-core profiles**: Using built-in profiles for containers, testing, and utilities
- **Custom profiles**: Creating your own profiles for different environments
- **Process labels**: Understanding and overriding resource requests by label
- **Resource limits**: Managing constrained environments with `resourceLimits`
- **Configuration inspection**: Using `nextflow config` to debug and verify settings

These configuration skills are transferable to any Nextflow pipeline and will help you run workflows efficiently across local machines, HPC clusters, and cloud platforms.

### What's next?

Congratulations on completing the Nextflow for Bioimaging course!

Next steps:

- Fill out the course survey to provide feedback
- Check out [Hello Nextflow](../hello_nextflow/index.md) to learn more about developing workflows
- Explore [Hello nf-core](../hello_nf-core/index.md) to dive deeper into nf-core tooling
- Browse other courses in the [training collections](../training_collections/index.md)
