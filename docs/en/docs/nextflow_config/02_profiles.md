# Part 2: Use profiles to switch configurations

In the last part, you accumulated a few configuration options: software packaging and resource allocations.
In practice, you'll often want to switch between whole sets of these options depending on where you're running, for example a laptop for development and an HPC cluster for production.

Nextflow lets you set up any number of [profiles](https://nextflow.io/docs/latest/config.html#profiles) describing different configurations, and select one (or several) at runtime with a single flag.

You've already used one: the `test` profile from [Nextflow Run](../nextflow_run/index.md) overrides the input parameters to a small, well-defined set.
Now you'll create your own infrastructure profiles and combine them with it.

---

## 1. Create profiles for different environments

### 1.1. Set up the profiles

Add two profiles to `nextflow.config`: one for running on a regular laptop with Docker, and one for a university HPC cluster with a Slurm scheduler and Conda.

=== "After"

    ```groovy title="nextflow.config" linenums="35" hl_lines="10-19"
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
        my_laptop {
            process.executor = 'local'
            docker.enabled = true
        }
        univ_hpc {
            process.executor = 'slurm'
            conda.enabled = true
            process.resourceLimits = [
                memory: 750.GB,
                cpus: 200,
                time: 30.d
            ]
        }
    }
    ```

=== "Before"

    ```groovy title="nextflow.config" linenums="35"
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

The `univ_hpc` profile also sets resource limits, since that's typically required on shared HPC infrastructure.

### 1.2. Run the workflow with a profile

Select a profile at runtime with `-profile`.

```bash
nextflow run main.nf -profile my_laptop
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [cheeky_goldstine] DSL2 - revision: c3c85dec78

    executor >  local (8)
    [da/0111e3] sayHello (1)       | 3 of 3 ✔
    [db/e19e47] convertToUpper (2) | 3 of 3 ✔
    [71/ad2bde] collectGreetings   | 1 of 1 ✔
    [85/afb958] cowpy              | 1 of 1 ✔
    ```

!!! warning

    The `univ_hpc` profile won't run in the training environment, since there's no Slurm scheduler available.

If you find other settings that always belong together, add them to the corresponding profile.
You can also create additional profiles to group any other combination you need.

### Takeaway

You know how to define profiles that bundle infrastructure-specific configuration, and select one at runtime with `-profile`.

### What's next?

Learn how to combine multiple profiles in a single run.

---

## 2. Combine profiles

Profiles aren't mutually exclusive.
You can activate several at once with `-profile <profile1>,<profile2>`.

### 2.1. Run with multiple profiles

Combine `my_laptop` with the `test` profile you already know from Nextflow Run.

```bash
nextflow run main.nf -profile my_laptop,test
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [sleepy_minsky] DSL2 - revision: c3c85dec78

    executor >  local (8)
    [b2/103b0d] sayHello (3)       | 3 of 3 ✔
    [00/6a7642] convertToUpper (2) | 3 of 3 ✔
    [b1/3df765] collectGreetings   | 1 of 1 ✔
    [d1/4e5ae1] cowpy              | 1 of 1 ✔
    ```

The individual file names correctly pick up `batch = 'test'` from the `test` profile (`COLLECTED-test-output.txt`, and so on).

If you combine profiles that set the same option, Nextflow resolves the conflict using whichever value it reads last, that is, whichever comes later in the file.
If the conflicting settings come from different configuration sources entirely, the standard [order of precedence](https://www.nextflow.io/docs/latest/config.html) applies.

### Takeaway

You know how to combine multiple profiles in a single run, and how Nextflow resolves conflicts when more than one profile sets the same option.

### What's next?

Learn how to inspect the fully resolved configuration before you run anything.

---

## 3. Inspect the resolved configuration

You already used `nextflow config -profile test` in [Nextflow Run](../nextflow_run/02_configure_pipeline.md) to check what a single profile resolves to.
That command becomes especially useful once you're combining multiple profiles: as you just saw, when two profiles set the same option, it can be tricky to work out by hand which value actually wins.
The `nextflow config` command resolves all of that for you, without running the pipeline.

### 3.1. Resolve the default configuration

```bash
nextflow config
```

??? success "Command output"

    ```groovy
    docker {
       enabled = true
    }

    process {
       memory = '1 GB'
       withName:cowpy {
          memory = '2 GB'
          cpus = 2
       }
    }

    params {
       input = 'data/greetings.csv'
       batch = 'batch'
       character = 'turkey'
    }
    ```

This is exactly what would apply if you ran the pipeline with no extra flags.

### 3.2. Resolve the configuration with profiles activated

Add the same profiles you'd use for an actual run.

```bash
nextflow config -profile my_laptop,test
```

??? success "Command output"

    ```groovy
    docker {
       enabled = true
    }

    process {
       memory = '1 GB'
       withName:cowpy {
          memory = '2 GB'
          cpus = 2
       }
       executor = 'local'
    }

    params {
       input = 'data/greetings.csv'
       batch = 'test'
       character = 'tux'
    }
    ```

Comparing the two confirms what changed: `params.batch`, `params.character`, and `process.executor` all reflect the `my_laptop,test` profiles.
This gets especially valuable for pipelines with many layers of configuration, where working out the resolved settings by hand would be tedious and error-prone.

### Takeaway

You know how to use `nextflow config` to inspect the fully resolved configuration for any combination of profiles, before running anything.

### What's next?

You've covered the essentials of configuring Nextflow pipelines.
See [Course summary](next_steps.md) for where to go from here.

---

## Summary

In this part you learned to:

- Define profiles that bundle infrastructure-specific configuration
- Combine multiple profiles in a single run, and understand how conflicts between them resolve
- Use `nextflow config` to inspect the fully resolved configuration
