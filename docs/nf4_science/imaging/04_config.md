# Part 4: Configuration

In Parts 1-3, we learned how to run molkart, manage inputs with parameter files and samplesheets, and train custom segmentation models.
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

```bash
cat nextflow.config
```

```groovy title="nextflow.config"
docker {
    enabled          = true
    runOptions       = '-u $(id -u):$(id -g) --entrypoint ""'
}

process {
    resourceLimits = [ cpus: 2, memory: 7.GB]
}
```

This minimal config:

- Enables Docker containers
- Sets Docker run options to avoid permission issues
- Caps resource requests at 2 CPUs and 7 GB (for GitHub Codespaces)

### Takeaway

Configuration files separate execution settings from workflow code, allowing you to adapt pipelines to different environments without modifying the workflow itself.

### What's next?

Learn about the built-in profiles that come with every nf-core pipeline.

---

## 2. nf-core built-in profiles

### 2.1. What are profiles?

Profiles are named sets of configuration that can be activated with the `-profile` flag via the `nextflow run` command.
They make it easy to switch between different compute scenarios without editing config files.

All nf-core pipelines come with standard profiles built in.

### 2.2. Container engine profiles

nf-core pipelines support multiple container engines:

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

### 2.3. Test profiles

Test profiles provide quick ways to verify the pipeline works:

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

Test profiles:

- `test`: Small dataset with fast parameters for quick testing
- `test_full`: More comprehensive test with larger data

### 2.4. Utility profiles

Additional profiles for specific use cases:

- `debug`: Enable detailed debugging information
- `wave`: Use Wave containers for on-demand container building
- `gitpod`: Optimized for Gitpod cloud development environment

### 2.5. Combining profiles

You can activate multiple profiles by separating them with commas:

```bash
nextflow run nf-core/molkart -profile docker,test --outdir results
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
        process.queue = 'batch'
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
nextflow config
```

View configuration with a specific profile:

```bash
nextflow config -profile local_dev
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

nf-core pipelines use **process labels** to standardize resource allocation across all pipelines.
Each process is tagged with a label like `process_low`, `process_medium`, or `process_high` to describe low, medium, or high compute resource requirements, respectively.

Here are the standard labels from molkart:

```groovy title="molkart/conf/base.config (excerpt)"
process {
    cpus   = { 1    * task.attempt }
    memory = { 6.GB * task.attempt }
    time   = { 4.h  * task.attempt }

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

Notice the `task.attempt` multiplier - this allows automatic retries with increased resources if a process fails.

### 4.2. Overriding resources by label

You can override all processes with a specific label:

```groovy title="nextflow.config"
process {
    withLabel:process_medium {
        cpus   = 4
        memory = 16.GB
    }
}
```

This changes resource requests for all `process_medium` processes at once.

### 4.3. Overriding resources for specific processes

For fine-grained control, target individual processes by name:

```groovy title="nextflow.config"
process {
    withName: 'NFCORE_MOLKART:MOLKART:CELLPOSE' {
        cpus   = 8
        memory = 32.GB
    }
}
```

!!! Tip

    To find process names, look at the pipeline execution output or check `.nextflow.log`.
    Process names follow the pattern `WORKFLOW:SUBWORKFLOW:PROCESS`.

### 4.4. Combining label and name selectors

You can use both approaches together:

```groovy title="nextflow.config"
process {
    // Set defaults for all high-resource processes
    withLabel:process_high {
        cpus   = 8
        memory = 32.GB
    }

    // Override specific process
    withName: 'NFCORE_MOLKART:MOLKART:CELLPOSE' {
        cpus   = 12
        memory = 64.GB
    }
}
```

### Takeaway

nf-core pipelines use process labels to standardize resource allocation.
You can override resources by label (affects multiple processes) or by name (affects one specific process).

### What's next?

Learn how to manage resource limits in constrained environments like GitHub Codespaces.

---

## 5. Managing resources in constrained environments

### 5.1. The resource limits problem

When we first ran molkart in GitHub Codespaces, we encountered a problem:

Some molkart processes request 12 GB or more of memory, but Codespaces only has ~8 GB available.

Without limits, Nextflow would refuse to run processes that exceed available resources.

### 5.2. Setting resource limits

The `resourceLimits` directive caps resource requests at specified values:

```groovy title="nextflow.config"
process {
    resourceLimits = [ cpus: 2, memory: 7.GB ]
}
```

This tells Nextflow: "If any process requests more than 2 CPUs or 7 GB memory, cap it at these limits."

This is why the test profile includes resourceLimits:

```groovy title="molkart/conf/test.config"
process {
    resourceLimits = [
        cpus: 4,
        memory: '15.GB',
        time: '1.h'
    ]
}
```

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
