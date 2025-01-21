# Part 6: Hello Config

This section will explore how to set up and manage the configuration of your Nextflow pipeline so that you'll be able to customize its behavior, adapt it to different environments, and optimize resource usage _without altering a single line of the workflow code itself_.

There are multiple ways to do this; here we are going to use the simplest and most common configuration file mechanism, the `nextflow.config` file.
Whenever there is a file named `nextflow.config` in the current directory, Nextflow will automatically load configuration from it.

!!!note

    Anything you put into the `nextflow.config` can be overridden at runtime by providing the relevant process directives or parameters and values on the command line, or by importing another configuration file, according to the order of precedence described [here](https://www.nextflow.io/docs/latest/config.html).

In this part of the training, we're going to use the `nextflow.config` file to demonstrate essential components of Nextflow configuration such as process directives, executors, profiles, and parameter files.
By learning to utilize these configuration options effectively, you can enhance the flexibility, scalability, and performance of your pipelines.

---

## 0. Warmup: Run the Hello Config workflow

[TODO] [usual run workflow to check everything works]

Verify that the initial workflow runs properly:

```bash
nextflow run hello-config.nf
```

This should run successfully:

```console title="Output"
Nextflow 24.09.2-edge is available - Please consider updating your version to it

 N E X T F L O W   ~  version 24.10.0

 ┃ Launching `main.nf` [tender_brahmagupta] DSL2 - revision: 848ff2f9b5

[TODO] [UPDATE]
```

If everything works, you're ready to learn how to modify basic configuration properties to adapt to your compute environment's requirements.

---

## 1. Determine what software packaging technology to use

The first step toward adapting your workflow configuration to your compute environment is specifying where the software packages that will get run in each step are going to be coming from.
Are they already installed in the local compute environment? Do we need to retrieve images and run them via a container system? Or do we need to retrieve Conda packages and build a local Conda environment?

In the very first part of this training course (Parts 1-4) we just used locally installed software in our workflow.
Then in Part 5, we introduced Docker containers, using the `-with-docker` command-line argument.

Now let's look at how we can configure Nextflow to use Docker or other container systems without having to specify that every time, using a `nextflow.config` file.

### 1.1. Enable Docker in the config file

There is a `nextflow.config` file in the current directory but it's a stub; there's nothing in it.

Let's add the line `docker.enabled = true` to the file.

```console title="nextflow.config" linenums="1"
docker.enabled = true
```

This instruction specifies that Nextflow should use Docker to run process calls that specify a Docker container image.

### 1.2. Run the workflow without the Docker CLI argument

```bash
nextflow run hello-config.nf
```

This should produce the following output:

```console title="Output"
TODO add updated output
```

This shows how you can get Nextflow to use Docker for any processes that specify a container with stating so everytime on the command line.

### 1.3. Disable Docker and enable Conda in the config file

Now, let's pretend we're working on an HPC cluster and the admin doesn't allow the use of Docker for security reasons.

Fortunately for us, Nextflow supports multiple other container technologies such as including Singularity (which is more widely used on HPC), and software package managers such as Conda.

We can change our configuration file to use Conda instead of Docker.
To do so, we switch the value of `docker.enabled` to `false`, and add a directive enabling the use of Conda:

_Before:_

```groovy title="nextflow.config" linenums="1"
docker.enabled = true
```

_After:_

```groovy title="nextflow.config" linenums="1"
docker.enabled = false
conda.enabled = true
```

This should allow Nextflow to create and utilize Conda environments for processes that have Conda packages specified. Which means we now need to add one to our `cowSay` process!

### 1.4. Specify a Conda package in the process definition

We've already retrieved the URI for a Conda package containing the `cowsay` tool:

`bioconda::cowsay==6.1`

!!! note

    There are a few different ways to get the URI for a given conda package.
    We recommend using the [Seqera Containers](https://seqera.io/containers/) search query, which will give you a URI that you can copy and paste, even if you're not planning to create a container from it.

Now we add the URI to the `cowSay` process definition using the `conda` directive:

_Before:_

```console title="hello-config.nf" linenums="22"
process cowSay {

    container 'community.wave.seqera.io/library/pip_cowsay:131d6a1b707a8e65'

    publishDir 'results', mode: 'copy'
```

_After:_

```console title="hello-config.nf" linenums="22"
process cowSay {

    container 'community.wave.seqera.io/library/pip_cowsay:131d6a1b707a8e65'
    conda 'bioconda::cowsay==6.1'

    publishDir 'results', mode: 'copy'
```

To be clear, we're not _replacing_ the `docker` directive, we're _adding_ an alternative option.

### 1.5. Run the workflow to verify that it can use Conda

Let's try it out.

```bash
nextflow run hello-config.nf
```

This may take a bit longer than usual the first time, and you might see the console output stay 'stuck' at this stage for a minute or so:

```console title="Output"
 N E X T F L O W   ~  version 24.10.0

 ┃ Launching `hello-config.nf` [extravagant_thompson] DSL2 - revision: 848ff2f9b5

TODO [UPDATE]
```

That's because Nextflow has to retrieve the Conda packages and create the environment, which takes a bit of work behind the scenes.
The good news is that you don't need to deal with any of it yourself!

After a few moments, it should spit out some more output, and eventually complete without error.

```console title="Output"
 N E X T F L O W   ~  version 24.10.0

 ┃ Launching `hello-config.nf` [silly_goldstine] DSL2 - revision: a60f9fd6af

[UPDATE]
```

From our standpoint, it looks like it works exactly the same as running with Docker, even though on the backend the mechanics are a bit different.

This means we're all set to run with Conda environments if needed.

!!!note

    Since these directives are assigned per process, it is possible 'mix and match', _i.e._ configure some of the processes in your workflow to run with Docker and others with Conda, for example, if the compute infrastructure you are using supports both.
    In that case, you would enable both Docker and Conda in your configuration file.
    If both are available for a given process, Nextflow will prioritize containers.

    And as noted earlier, Nextflow supports multiple other software packaging and container technologies, so you are not limited to just those two.

### Takeaway

You know how to configure which software package each process should use, and how to switch between technologies.

### What's next?

Learn how to change the executor used by Nextflow to actually do the work.

---

## 2. Determine what executor(s) should be used to do the work

Until now, we have been running our pipeline with the local executor.
This runs each step on the same machine that Nextflow is running on.
However, for large workloads, you will typically want to use a distributed executor such as an HPC or cloud.
Nextflow supports several different distributed executors, including:

-   HPC (SLURM, PBS, SGE)
-   AWS Batch
-   Google Batch
-   Azure Batch
-   Kubernetes
-   GA4GH TES

The executor is subject to a process directive called `executor`. By default it is set to `local`, so the following configuration is implied:

```groovy title="Built-in configuration"
process {
    executor = 'local'
}
```

Let's look at what it would take to using a Slurm scheduler, assuming we had a connection to a cluster and Slurm was installed appropriately.

!!! warning

    What follows is for demonstration purposes but **will not execute the work** since we don't have access to an external executor.

### 2.1. Set up a Slurm executor

Add the following lines to the `nextflow.config` file:

```groovy title="nextflow.config" linenums="12"
process {
    executor = 'slurm'
}
```

And... that's it! As noted before, this does assume that Slurm itself is already set up for you, but this is really all Nextflow itself needs to know.

Basically we are telling Nextflow to generate a Slurm submission script and submit it using an `sbatch` command.

### 2.2. Launch the workflow to generate the job submission script

TODO: THIS WAS CONFUSING — EXPLAIN BETTER OR CUT

Let's try running this; even though we know it won't execute in the training environment, we'll be able to see what the submission script looks like.

```bash
nextflow run main.nf -profile conda_on
```

As expected, this fails with a fairly unambiguous error:

```console title="Output"
nextflow
 N E X T F L O W   ~  version 24.10.0

 ┃ Launching `hello-config.nf` [grave_gauss] DSL2 - revision: 66cd7c255a

[UPDATE]
Caused by:
  java.io.IOException: Cannot run program "sbatch" (in directory "/workspace/gitpod/hello-nextflow/hello-config/work/eb/2962ce167b3025a41ece6ce6d7efc2"): error=2, No such file or directory

Command executed:

  sbatch .command.run
```

However, it did produce what we are looking for: the `.command.run` file that Nextflow tried to submit to Slurm via the `sbatch` command.

Let's take a look inside. <!-- **TODO: UPDATE NEXTFLOW VERSION SO WE CAN HAVE THIS SWEET OUTPUT** -->

```bash title=".command.run" linenums="1"
#!/bin/bash
TODO: UPDATE
```

This shows the job submission details that Nextflow is trying to hand over to Slurm.

You can try using any of the other supported executors in the same way. Nextflow will translate the values submitted to the executor into the appropriate equivalent instructions.

!!! warning

    Before continuing the training, either delete the line `executor = 'slurm'` or change it to `executor = 'local'` in your configuration file.
    You will not be able to run subsequent commands otherwise.

### Takeaway

You now know how to change the executor to use different kinds of computing infrastructure.

### What's next?

Learn how to control the resources allocated for executing processes.

---

## 3. Allocate compute resources with process directives

Most high-performance computing platforms allow (and sometimes require) that you specify certain resource allocation parameters such as number of CPUs and memory.

By default, Nextflow will use a single CPU and 2GB of memory for each process.
The corresponding process directives are called `cpus` and `memory`, so the following configuration is implied:

```groovy title="Built-in configuration"
process {
    cpus = 1
    memory = 2.GB
}
```

You can modify these values, either for all processes or for specific named processes, using additional process directives in your configuration file.
Nextflow will translate them into the appropriate instructions for the chosen executor.

But how do you know what values to use?

### 3.1. Run the workflow to generate a resource utilization report

If you don't know up front how much CPU and memory your processes are likely to need, you can do some resource profiling, meaning you run the workflow with some default allocations, record how much each process used, and from there, estimate how to adjust the base allocations.

Conveniently, Nextflow includes built-in tools for doing this, and will happily generate a report for you on request.

To do so, add `-with-report <filename>.html` to your command line.

```bash
nextflow run hello-config.nf -with-report report-config-1.html
```

The report is an html file, which you can download and open in your browser. You can also right click it in the file explorer on the left and click on `Show preview` in order to view it in the training environment.

Take a few minutes to look through the report and see if you can identify some opportunities for adjusting resources.
Make sure to click on the tabs that show the utilization results as a percentage of what was allocated.
There is some [documentation](https://www.nextflow.io/docs/latest/reports.html) describing all the available features.

<!-- TODO: insert images -->

### 3.2. Set resource allocations for all processes

The profiling shows that the processes in our training workflow are very lightweight, so let's reduce the default memory allocation to 1GB per process.

Add the following to your `nextflow.config` file:

```groovy title="nextflow.config" linenums="20"
process {
    memory = 1.GB
}
```

### 3.3. Set resource allocations for an individual process

At the same time, we're going to pretend that the `cowSay` process requires more resources than the others, just so we can demonstrate how to adjust allocations for an individual process.

_Before:_

```groovy title="nextflow.config" linenums="1"
process {
    memory = 1.GB
}
```

_After:_

```groovy title="nextflow.config" linenums="1"
process {
    memory = 1.GB
    withName: 'cowSay' {
        memory = 2.GB
        cpus = 2
    }
}
```

With this configuration, all processes will request 1GB of memory and a single CPU (the implied default), except the `cowSay` process, which will request 2GB and 2 CPUs.

!!! note

    If you have a machine with few CPUs and you allocate a high number per process, you might see process calls getting queued behind each other.
    This is because Nextflow ensures we don't request more CPUs than are available.

### 3.4. Run the workflow with the modified configuration

Let's try that out, supplying a different filename for the profiling report so we can compare performance before and after the configuration changes.

```bash
nextflow run hello-config.nf -with-report report-config-2.html
```

You will probably not notice any real difference since this is such a small workload, but this is the approach you would use to analyze the performance and resource requirements of a real-world workflow.

It is very useful when your processes have different resource requirements. It empowers you to right-size the resource allocations you set up for each process based on actual data, not guesswork.

!!!note

    This is just a tiny taster of what you can do to optimize your use of resources.
    Nextflow itself has some really neat [dynamic retry logic](https://training.nextflow.io/basic_training/debugging/#dynamic-resources-allocation) built in to retry jobs that fail due to resource limitations.
    Additionally, the Seqera Platform offers AI-driven tooling for optimizing your resource allocations automatically as well.

    We'll cover both of those approaches in an upcoming part of this training course.

### 3.5. Add resource limits

Depending on what computing executor and compute infrastructure you're using, there may be some constraints on what you can (or must) allocate.
For example, your cluster may require you to stay within certain limits.

You can use the `resourceLimits` directive to set the relevant limitations. The syntax looks like this when it's by itself in a process block:

```groovy title="Syntax example"
process {
    resourceLimits = [
        memory: 750.GB,
        cpus: 200,
        time: 30.d
    ]
}
```

Nextflow will translate these values into the appropriate instructions depending on the executor that you specified.

We're not going to run this, since we don't have access to relevant infrastructure in the training environment.
However, if you were to try running the workflow with resource allocations that exceed these limits, then look up the `sbatch` command in the `.command.run` script file, you would see that the requests that actually get sent to the executor are capped at the values specified by `resourceLimits`.

!!!note

    The nf-core project has compiled a [collection of configuration files](https://nf-co.re/configs/) shared by various institutions around the world, covering a wide range of HPC and cloud executors.

    Those shared configs are valuable both for people who work there and can therefore just utilize their institution's configuration out of the box, and as a model for people who are looking to develop a configuration for their own infrastructure.

### Takeaway

You know how to generate a profiling report to assess resource utilization and how to modify resource allocations for all processes and/or for individual processes, as well as set resource limitations for running on HPC.

### What's next?

Learn to use a parameter file to store workflow parameters.

---

## 4. Use a parameter file to store workflow parameters

So far we've been looking at configuration from the technical point of view of the compute infrastructure.
Now let's consider another aspect of workflow configuration that is very important for reproducibility: the configuration of the workflow parameters.

Currently, our workflow is set up to accept several parameter values via the command-line, with default values set in the workflow script itself.
This is fine for a simple workflow with very few parameters that need to be set for a given run.
However, many real-world workflows will have many more parameters that may be run-specific, and putting all of them in the commend line would be tedious and error-prone.

Nextflow allows us to specify parameters via a parameter file in JSON format, which makes it very convenient to manage and distribute alternative sets of default values, for example, as well as run-specific parameter values.

We provide an example parameter file in the current directory, called `test-params.json`:

```json title="test-params.json" linenums="1"
{
    "greet": "Dobrý den",
    "batch": "Trio",
    "character": "pig"
}
```

This parameter file contains a key-value pair for each of the inputs our workflow expects.

### 4.1. Run the workflow using a parameter file

To run the workflow with this parameter file, simply add `-params-file <filename>` to the base command.

```bash
nextflow run hello-config.nf -params-file test-params.json
```

It works! And as expected, this produces the same outputs as previously.

```console title="Output"
TODO: UPDATE OUTPUT
```

This may seem like overkill when you only have a few parameter to specify, but some pipelines expect dozens of parameters.
In those cases, using a parameter file will allow us to provide parameter values at runtime without having to type massive command lines and without modifying the workflow script.

### Takeaway

You know how to manage parameter defaults and override them at runtime using a parameter file.

### What's next?

Learn how to use profiles to conveniently switch between alternative configurations.

---

## 5. Use profiles to select preset configurations

You may want to switch between alternative settings depending on what computing infrastructure you're using. For example, you might want to develop and run small-scale tests locally on your laptop, then run full-scale workloads on HPC or cloud.

Nextflow lets you set up profiles that describe different configurations, which you can then select at runtime using a command-line argument, rather than having to modify the configuration file itself.

### 5.1. Create profiles for switching between local development and execution on HPC

Let's set up two alternative profiles; one for running small scale loads on a regular computer, where we'll use Docker containers, and one for running on a university HPC with a Slurm scheduler, where we'll use Conda packages.

Add the following to your `nextflow.config` file:

```groovy title="nextflow.config" linenums="3"
profiles {
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

You see that for the university HPC, we're also specifying resource limitations.

### 5.2. Run the workflow with a profile

To specify a profile in our Nextflow command line, we use the `-profile` argument.

Let's try running the workflow with the `my_laptop` configuration.

```bash
nextflow run hello-config.nf -profile my_laptop
```

This should produce the following output:

```
TODO UPDATE
```

As you can see, this allows us to toggle between configurations very conveniently at runtime.

!!! warning

    The `univ_hpc` profile will not run properly in the training environment since we do not have access to a Slurm scheduler.

If in the future we find other elements of configuration that are always co-occurring with these, we can simply add them to the corresponding profile(s).
We can also create additional profiles if there are other elements of configuration that we want to group together.

### 5.3. Create a test profile

Profiles are not only for infrastructure configuration.
We can also use them to set default values for workflow parameters, to make it easier for others to try out the workflow without having to gather appropriate input values themselves.
This is intended as an alternative to using a parameter file.

The syntax for expressing default values is the same as when writing them into the workflow file itself, except we wrap them in a block named `test`:

```groovy title="Syntax example"
    test {
        params.<parameter1>
        params.<parameter2>
        ...
    }
```

If we add a test profile for our workflow, the `profiles` block becomes:

```groovy title="nextflow.config" linenums="3"
profiles {
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
    test { TODO UPDATE THIS
        params.greetings    = "greetings.csv"
        params.thing        = "thing.txt"
    }
}
```

Just like for technical configuration profiles, you can set up multiple different profiles specifying parameters under any arbitrary name you like.

### 5.4. Run the workflow locally with the test profile

Conveniently, profiles are not mutually exclusive, so we can specify multiple profiles in our command line using the following syntax `-profile <profile1>,<profile2>` (for any number of profiles).

!!! note

    If you combine profiles that set values for the same elements of configuration and are described in the same configuration file, Nextflow will resolve the conflict by using whichever value it read in last (_i.e._ whatever comes later in the file).
    If the conflicting settings are set in different configuration sources, the default [order of precedence](https://www.nextflow.io/docs/latest/config.html) applies.

Let's try adding the test profile to our previous command:

```bash
nextflow run hello-config.nf -profile my_laptop,test
```

This should produce the following:

```console title="Output"
TODO: UPDATE
```

This means that as long as we distribute any test data files with the workflow code, anyone can quickly try out the workflow without having to supply their own inputs via the command line or a parameter file.

!!! note

    We can even point to URLs for larger files that are stored externally.
    Nextflow will download them automatically as long as there is an open connection.

### Takeaway

You know how to use profiles to select a preset configuration at runtime with minimal hassle. More generally, you know how to configure your workflow executions to suit different compute platforms and enhance the reproducibility of your analyses.

### What's next?

Celebrate and give yourself a big pat on the back! You have completed your very first Nextflow developer course.
Then check out the training portal homepage for more training content that may be of interest.
