# Part 2: Configure the pipeline

In [Part 1](./01_run_nextflow.md), you ran a complete multi-step pipeline that processes multiple inputs in parallel using containers.
Now we're going to look at how to configure pipeline behavior using `nextflow.config`: first by examining the configuration file we already gave you, then by exploring a couple of other ways to supply configuration, and finally by controlling how and where outputs get published.

---

## 1. Examine the main configuration file

Nextflow automatically picks up `nextflow.config` from the working directory and applies its settings to every run.

We provide you with a configuration file that covers four areas: software packaging, process settings, pipeline parameters, and execution profiles.

??? full-code "nextflow.config"

    ```groovy title="nextflow.config" linenums="1"
    /*
     * Software packaging
     */
    docker.enabled = true

    /*
     * Process settings
     */
    process {
        cpus = 1
        memory = 1.GB
    }

    /*
     * Pipeline parameters
     */
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }

    /*
     * Profiles
     */
    profiles {
        test {
            params.input = 'data/greetings.csv'
            params.batch = 'test'
            params.character = 'tux'
        }
        conda {
            docker.enabled = false
            conda.enabled = true
        }
    }
    ```

Let's go through each one, then put profiles to use by running the pipeline with one.

!!! note

    This config covers local execution on a single machine.
    Nextflow also supports HPC schedulers (SLURM, PBS, LSF) and cloud executors (AWS Batch, Google Cloud Batch, Azure Batch), all configured through the same `nextflow.config` mechanism.
    See [Part 1: Adapt to your compute environment](../nextflow_config/01_packaging_execution_resources.md) in the [Nextflow Config](../nextflow_config/index.md) course for a full walkthrough of these options.

### 1.1. Software packaging

Software packaging is how Nextflow supplies the actual tools your processes need, whether that's a container image, a Conda environment, or something else.

```groovy title="nextflow.config" linenums="1"
/*
 * Software packaging
 */
docker.enabled = true
```

This line enables Docker for every process.
Any process that declares a `container` directive runs inside the specified image.

### 1.2. Process settings

Remember that a process is a single step in your pipeline, like `sayHello` or `cowpy`.
Nextflow lets you configure a number of things about how each one actually runs: how much CPU and memory it gets, which container or Conda environment it uses, and more.

```groovy title="nextflow.config" linenums="6"
/*
 * Process settings
 */
process {
    cpus = 1
    memory = 1.GB
}
```

This caps every process at a single CPU and 1 GB of memory.

Nextflow also lets you set different values for individual named processes or groups of processes; you'll learn how in [Part 1: Adapt to your compute environment](../nextflow_config/01_packaging_execution_resources.md#32-set-resource-allocations-for-a-specific-process) of the [Nextflow Config](../nextflow_config/index.md) course.

### 1.3. Pipeline parameters

Parameters are the pipeline's command-line inputs, the same `--input`, `--batch` and `--character` flags you've already been setting directly on the command line.
Setting defaults for them here means you don't have to type them out every time, though as you'll see later in this part, there are a couple of other ways to supply them too.

```groovy title="nextflow.config" linenums="14"
/*
 * Pipeline parameters
 */
params {
    input = 'data/greetings.csv'
    batch = 'batch'
    character = 'turkey'
}
```

These defaults kick in whenever a parameter isn't supplied on the command line, so running `nextflow run main.nf` with no flags still works.

### 1.4. Profiles

Profiles let you bundle up a set of settings under a single name, so you can switch between whole configurations with one flag instead of changing values by hand every time.

```groovy title="nextflow.config" linenums="23"
/*
 * Profiles
 */
profiles {
    test {
        params.input = 'data/greetings.csv'
        params.batch = 'test'
        params.character = 'tux'
    }
    conda {
        docker.enabled = false
        conda.enabled = true
    }
}
```

The `test` profile overrides three parameters to run the pipeline with a small, well-defined input set; every nf-core pipeline ships with one of these for quick validation, and it's a convention worth following in your own pipelines too.

The `conda` profile switches software packaging from Docker to Conda.

You activate a profile by passing `-profile <name>` on the command line.

Let's put the `test` profile to use.

```bash
nextflow run main.nf -profile test
```

??? success "Command output"

    ```console hl_lines="6 7 8 9"
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [reverent_heisenberg] DSL2 - revision: ce74f81996

    executor >  local (8)
    [3d/8a12c7] sayHello (3)       | 3 of 3 ✔
    [e7/e0934f] convertToUpper (2) | 3 of 3 ✔
    [1d/616569] collectGreetings   | 1 of 1 ✔
    [44/7d46cf] cowpy              | 1 of 1 ✔
    ```

The pipeline runs with `batch = 'test'` and `character = 'tux'`.
Check `results/test/`: the batch name is now part of the directory path itself, and the ASCII art features the tux penguin instead of a turkey.

!!! note

    You can activate several profiles at once, and use `nextflow config -profile <name>,<name>` to see the fully resolved result before running anything.
    Combining profiles, and how Nextflow resolves conflicts between them, is covered in depth in [Part 2: Use profiles to switch configurations](../nextflow_config/02_profiles.md) of the [Nextflow Config](../nextflow_config/index.md) course.

### Takeaway

You know what the most common elements of a `nextflow.config` file do, and how to activate a profile.

### What's next?

Learn a couple of other ways to supply configuration values without modifying the main `nextflow.config` file, useful for configuring individual runs and for sharing an exact set of settings with someone else.

---

## 2. Provide configuration via supplemental files

Setting defaults in `nextflow.config` works well for values that rarely change.
Nextflow also gives you two more targeted mechanisms: a run-specific configuration file for adapting execution to a particular environment, and a parameter file for sharing an exact set of input values with a collaborator.

### 2.1. Use a run-specific configuration file

Say you're moving the pipeline to a machine that doesn't have Docker, and you want to give every process more room to work with.
Create a new configuration file with just the overrides you need:

```groovy title="custom.config" linenums="1"
process {
    cpus = 2
    memory = 2.GB
}

docker.enabled = false
conda.enabled = true
```

Pass it alongside your main pipeline with `-c`:

```bash
nextflow run main.nf -c custom.config
```

??? success "Command output"

    ```console hl_lines="6 7 8 9"
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [exotic_cray] DSL2 - revision: ce74f81996

    executor >  local (8)
    [77/e02315] sayHello (1)       | 3 of 3 ✔
    [a6/ccf44b] convertToUpper (3) | 3 of 3 ✔
    [56/fd1296] collectGreetings   | 1 of 1 ✔
    [2c/205a94] cowpy              | 1 of 1 ✔
    ```

Nextflow merges `custom.config` on top of the pipeline's own `nextflow.config`, so every process now gets 2 CPUs and 2 GB of memory instead of the defaults, and runs through Conda instead of Docker.
`cowpy` is the only process with a Conda package declared alongside its container, so it's the one you'll see Nextflow actually build an environment for:

```console
Creating env using conda: conda-forge::cowpy==1.1.5 [cache /path/to/work/conda/env-898314d566668b6587ad714ae06b8520]
```

A small file that only overrides resource allocation and packaging, without touching pipeline parameters, is exactly the pattern nf-core pipelines expect from institutional configs.
Browse the [nf-core/configs](https://github.com/nf-core/configs) repository for real-world examples.

That gives you a disposable way to adapt a pipeline to a new environment without touching your normal configuration.

### 2.2. Use a parameter file

Say instead you need to share an exact set of run parameters with a collaborator, or record them for a publication.

Nextflow allows you to supply [parameter files](https://nextflow.io/docs/latest/config.html#parameter-file) in YAML or JSON format, which are a simpler way to distribute an exact, reproducible set of values.

A parameter file called `test-params.yaml` is already provided in your working directory:

```yaml title="test-params.yaml" linenums="1"
input: "data/greetings.csv"
batch: "yaml"
character: "stegosaurus"
```

The syntax uses colons (`:`) instead of the equal signs (`=`) used in `nextflow.config`, since this file is plain YAML rather than Groovy.

!!! info

    A JSON version, `test-params.json`, is also provided. Feel free to try it on your own; the syntax for passing it is identical.

Pass the file with `-params-file`:

```bash
nextflow run main.nf -params-file test-params.yaml
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [sharp_faraday] DSL2 - revision: ce74f81996

    executor >  local (8)
    [1c/9ff63e] sayHello (1)       | 3 of 3 ✔
    [3b/bb5691] convertToUpper (2) | 3 of 3 ✔
    [cd/2c1f6e] collectGreetings   | 1 of 1 ✔
    [89/c333bc] cowpy              | 1 of 1 ✔
    ```

??? abstract "File contents"

    ```console title="results/yaml/cowpy-COLLECTED-yaml-output.txt"
     _________
    / BONJOUR \
    | HOLA    |
    \ HELLO   /
     ---------
    \                             .       .
     \                           / `.   .' "
      \                  .---.  <    > <    >  .---.
       \                 |    \  \ - ~ ~ - /  /    |
             _____          ..-~             ~-..-~
            |     |   \~~~\.'                    `./~~~/
           ---------   \__/                        \__/
          .'  O    \     /               /       \  "
         (_____,    `._.'               |         }  \/~~~/
          `----.          /       }     |        /    \__/
                `-.      |       /      |       /      `. ,~~|
                    ~-.__|      /_ - ~ ^|      /- _      `..-'
                        |     /        |     /     ~-.     `-. _  _  _
                        |_____|        |_____|         ~ - . _ _ _ _ _>
    ```

A parameter file is especially valuable once a pipeline has more than a handful of parameters: it lets you supply them all at once, without a sprawling command line or any change to the workflow script, and it's easy to distribute alongside your results.

### Takeaway

You know two more ways to supply configuration: a run-specific configuration file for adapting execution to a new environment, and a parameter file for sharing exact, reproducible input values.

### What's next?

Learn how to control how and where your pipeline's outputs get published.

---

## 3. Manage pipeline outputs

A pipeline author decides how outputs are organized in code, but you don't need to touch that code to control where they end up or how they get there.
Nextflow gives you config-level ways to do that instead: set a base output directory, and choose whether files get copied or symlinked.

### 3.1. Customize the output directory

By default, Nextflow publishes outputs under `results/`.
Point it elsewhere with `-output-dir` (or its short form, `-o`):

```bash
nextflow run main.nf -output-dir outputs
```

??? success "Command output"

    ```console hl_lines="6 7 8 9"
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [serene_kimura] DSL2 - revision: ce74f81996

    executor >  local (8)
    [31/df5c15] sayHello (3)       | 3 of 3 ✔
    [08/8bb2e5] convertToUpper (1) | 3 of 3 ✔
    [e5/5814da] collectGreetings   | 1 of 1 ✔
    [cf/8ab8c6] cowpy              | 1 of 1 ✔
    ```

??? abstract "Directory contents"

    ```console
    outputs/batch
    ├── batch-report.txt
    ├── cowpy-COLLECTED-batch-output.txt
    └── intermediates
        ├── Bonjour-output.txt
        ├── COLLECTED-batch-output.txt
        ├── Hello-output.txt
        ├── Hola-output.txt
        ├── UPPER-Bonjour-output.txt
        ├── UPPER-Hello-output.txt
        └── UPPER-Hola-output.txt
    ```

The outputs now land under `outputs/batch/` instead of the built-in `results/batch/` default.
The pipeline's own code still decides the structure within that base directory, like the `batch/` and `intermediates/` subdirectories; `-output-dir` only controls where that structure starts.

`-output-dir` is really just a command-line shortcut for the `outputDir` configuration option, so it can go anywhere configuration can: directly in `nextflow.config`, inside a profile, or in a `-c` overlay file like the one you used earlier in this part.
For example, this snippet shows the same setting placed directly in `nextflow.config` instead of passed on the command line:

```groovy title="nextflow.config"
outputDir = 'outputs'
```

See [Configuration file](https://nextflow.io/docs/latest/config.html) in the Nextflow reference for the full list of places a configuration option like this can live.

### 3.2. Choose how outputs get published

By default, Nextflow publishes outputs as symlinks that point to the locations of the outputs under `work/`, not real copies:

```console
$ ls -l results/batch/intermediates/Hello-output.txt
lrwxr-xr-x  ...  Hello-output.txt -> /workspaces/training/nextflow-run/work/b7/b8c4d1.../Hello-output.txt
```

Pipeline authors can set the 'publish mode' to either `'copy'` or `'move'` for each individual process in the workflow code.
They typically do this for the final outputs of the pipeline, while leaving the default `'symlink'` behavior set for intermediate files that can be deleted once the full pipeline has been run.

That avoids duplicating data on disk, but it means you can't delete the task directories under `work/` without breaking the link, losing the ability to use `-resume`.
If you want all output files to be properly copied instead, set [`workflow.output.mode`](https://nextflow.io/docs/latest/reference/config.html#workflow) to `'copy'` in your pipeline configuration. (Unlike `-output-dir`, there's no command-line flag for this; it's config-only.)

Try setting it in `nextflow.config`:

=== "After"

    ```groovy title="nextflow.config" hl_lines="7"
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }

    workflow.output.mode = 'copy'
    ```

=== "Before"

    ```groovy title="nextflow.config"
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }
    ```

Then run the pipeline, changing the batch name so that you can see the difference in the outputs:

```bash
nextflow run main.nf --batch withmode
```

??? success "Command output"

    ```console hl_lines="6 7 8 9"
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [angry_noether] DSL2 - revision: ce74f81996

    executor >  local (8)
    [41/2b478d] sayHello (3)       | 3 of 3 ✔
    [bf/dd2840] convertToUpper (2) | 3 of 3 ✔
    [ea/364e97] collectGreetings   | 1 of 1 ✔
    [10/76fe7b] cowpy              | 1 of 1 ✔
    ```

Have a look at one of the output files like before:

```console
$ ls -l results/withmode/intermediates/Hello-output.txt
-rw-r--r--  ...  Hello-output.txt
```

Now it's a real, independent file that will stay available even if `work/` gets cleaned up.

!!! warning

    The `workflow.output.mode` setting only fills in a default for outputs that don't already have a mode set in the pipeline code.
    It cannot override a mode the author hardcoded, no matter what you set it to.

### Takeaway

You know how to customize the base output directory and choose between copied and symlinked outputs, both without touching the pipeline's code.

### What's next?

Head on to [Part 3](./03_manage_executions.md), where you'll learn how to inspect the history of past runs, generate execution reports, and clean up old work directories.

---

## Summary

In this part you learned to:

- Configure pipeline behavior using `nextflow.config` and profiles
- Supply configuration via a run-specific configuration file or a parameter file
- Customize the output directory and choose between copied and symlinked outputs
