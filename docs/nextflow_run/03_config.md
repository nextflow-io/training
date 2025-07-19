# Part 3: Configuration

This section will explore how to manage the configuration of a Nextflow pipeline in order to customize its behavior, adapt it to different environments, and optimize resource usage _without altering a single line of the workflow code itself_.

There are multiple ways to do this; here we are going to use the simplest and most common configuration file mechanism, the `nextflow.config` file.
As noted previously, whenever there is a file named `nextflow.config` in the current directory, Nextflow will automatically load configuration from it.

!!! Tip

    Anything you put into the `nextflow.config` can be overridden at runtime by providing the relevant process directives or parameters and values on the command line, or by importing another configuration file, according to the order of precedence described [here](https://www.nextflow.io/docs/latest/config.html).

In this part of the training, we're going to use the `nextflow.config` file to demonstrate essential components of Nextflow configuration such as process directives, executors, profiles, and parameter files.
By learning to utilize these configuration options effectively, you can enhance the flexibility, scalability, and performance of your pipelines.

To exercise these elements of configuration, we're going to be running a fresh copy of the workflow we last ran at the end of Part 2 of this training course, renamed `3-main.nf`.

---

## 1. Determine what software packaging technology to use

The first step toward adapting your workflow configuration to your compute environment is specifying where the software packages that will get run in each step are going to be coming from.
Are they already installed in the local compute environment?
Do we need to retrieve images and run them via a container system?
Or do we need to retrieve Conda packages and build a local Conda environment?

For most of this training course so far, we just used locally installed software in our workflow.
Then in the last section of Part 2, we introduced Docker containers and the `nextflow.config` file, which we used to enable the use of Docker containers.

Now let's see how we can configure an alternative software packaging option via the `nextflow.config` file.

### 1.1. Disable Docker and enable Conda in the config file

Let's pretend we're working on an HPC cluster and the admin doesn't allow the use of Docker for security reasons.

Fortunately for us, Nextflow supports multiple other container technologies such as including Singularity/Apptainer (which is more widely used on HPC), and software package managers such as Conda.

We can change our configuration file to use Conda instead of Docker.
To do so, we switch the value of `docker.enabled` to `false`, and add a directive enabling the use of Conda:

=== "After"

    ```groovy title="nextflow.config" linenums="1" hl_lines="1-2"
    docker.enabled = false
    conda.enabled = true
    ```

=== "Before"

    ```groovy title="nextflow.config" linenums="1"

    docker.enabled = true
    ```

This will allow Nextflow to create and utilize Conda environments for processes that have Conda packages specified.
Which means we now need to add one of those to the `cowpy` process definition!

### 1.2. Specify a Conda package in the process definition

We've already retrieved the URI for a Conda package containing the `cowpy` tool: `conda-forge::cowpy==1.1.5`

!!! Tip

    There are a few different ways to get the URI for a given conda package.
    We recommend using the [Seqera Containers](https://seqera.io/containers/) search query, which will give you a URI that you can copy and paste, even if you're not planning to create a container from it.

Now we add the URI to the `cowpy` process definition using the `conda` directive:

=== "After"

    ```console title="modules/cowpy.nf" linenums="4" hl_lines="4"
    process cowpy {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
        conda 'conda-forge::cowpy==1.1.5'

        publishDir 'results', mode: 'copy'
    ```

=== "Before"

    ```console title="modules/cowpy.nf" linenums="4"
    process cowpy {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

        publishDir 'results', mode: 'copy'
    ```

To be clear, we're not _replacing_ the `container` directive, we're _adding_ an alternative option.

### 1.3. Run the workflow to verify that it can use Conda

Let's try it out.

```bash
nextflow run 3-main.nf --input greetings.csv --character turkey
```

This should work without error.

<details>
  <summary>Command output</summary>

```console title="Output"
 N E X T F L O W   ~  version 25.04.3

Launching `3-main.nf` [trusting_lovelace] DSL2 - revision: 028a841db1

executor >  local (8)
[ee/4ca1f2] sayHello (3)       | 3 of 3 ✔
[20/2596a7] convertToUpper (1) | 3 of 3 ✔
[b3/e15de5] collectGreetings   | 1 of 1 ✔
[c5/af5f88] cowpy              | 1 of 1 ✔
```

</details>

Behind the scenes, Nextflow has retrieved the Conda packages and created the environment, which normally takes a bit of work; so it's nice that we don't have to do any of that ourselves!

!!! Tip

    This runs quickly because the `cowpy` package is quite small, but if you're working with large packages, it may take a bit longer than usual the first time, and you might see the console output stay 'stuck' for a minute or so before completing.
    This is normal and is due to the extra work Nextflow does the first time you use a new package.

From our standpoint, it looks like it works exactly the same as running with Docker, even though on the backend the mechanics are a bit different.

This means we're all set to run with Conda environments if needed.

!!! Tip

    Since these directives are assigned per process, it is possible 'mix and match', _i.e._ configure some of the processes in your workflow to run with Docker and others with Conda, for example, if the compute infrastructure you are using supports both.
    In that case, you would enable both Docker and Conda in your configuration file.
    If both are available for a given process, Nextflow will prioritize containers.

    And as noted earlier, Nextflow supports multiple other software packaging and container technologies, so you are not limited to just those two.

### Takeaway

You know how to configure which software package each process should use, and how to switch between technologies.

### What's next?

Learn how to specify what executor Nextflow should use to actually do the work.

---

## 2. Specify what executor should be used to do the work

Until now, we have been running our pipeline with the local executor.
This executes each task on the machine that Nextflow is running on.
When Nextflow begins, it looks at the available CPUs and memory.
If the resources of the tasks ready to run exceed the available resources, Nextflow will hold the last tasks back from execution until one or more of the earlier tasks have finished, freeing up the necessary resources.

For very large workloads, you may discover that your local machine is a bottleneck, either because you have a single task that requires more resources than you have available, or because you have so many tasks that waiting for a single machine to run them would take too long.
The local executor is convenient and efficient, but is limited to that single machine.
Nextflow supports [many different execution backends](https://www.nextflow.io/docs/latest/executor.html), including HPC schedulers (Slurm, LSF, SGE, PBS, Moab, OAR, Bridge, HTCondor and others) as well as cloud execution backends such (AWS Batch, Google Cloud Batch, Azure Batch, Kubernetes and more).

### 2.1. Targeting a different backend

The choice of executor is set by a process directive called `executor`.
By default it is set to `local`, so the following configuration is implied:

```groovy title="Built-in configuration"
process {
    executor = 'local'
}
```

To set the executor to target a different backend, simply specify the executor you want using similar syntax as described above for resource allocations (see [documentation](https://www.nextflow.io/docs/latest/executor.html) for all options).

```groovy title="nextflow.config"
process {
    executor = 'slurm'
}
```

!!! Warning

    We can't actually test this in the training environment because it's not set up to connect to an HPC.

### 2.2. Dealing with backend-specific syntax for execution parameters

Most high-performance computing platforms allow (and sometimes require) that you specify certain parameters such as resource allocation requests and limitations (for e.g. number of CPUs and memory) and name of the job queue to use.

Unfortunately, each of these systems uses different technologies, syntaxes and configurations for defining how a job should be defined and submitted to the relevant scheduler.

For example, a job requiring 8 CPUs and 4GB of RAM to be executed on the queue "my-science-work" needs to be expressed in the following ways depending on the backend:

<details>
  <summary>Config for SLURM / submit using `sbatch`</summary>

```bash
#SBATCH -o /path/to/my/task/directory/my-task-1.log
#SBATCH --no-requeue
#SBATCH -c 8
#SBATCH --mem 4096M
#SBATCH -p my-science-work
```

</details>

<details>
  <summary>Config for PBS / submit using `qsub`</summary>

```bash
#PBS -o /path/to/my/task/directory/my-task-1.log
#PBS -j oe
#PBS -q my-science-work
#PBS -l nodes=1:ppn=5
#PBS -l mem=4gb
```

</details>

<details>
  <summary>Config for SGE / submit using `qsub`</summary>

```bash
#$ -o /path/to/my/task/directory/my-task-1.log
#$ -j y
#$ -terse
#$ -notify
#$ -q my-science-work
#$ -l slots=5
#$ -l h_rss=4096M,mem_free=4096M
```

</details>

Fortunately, Nextflow simplifies all of this.
It provides a standardized syntax so that you can specify the relevant properties such as `cpus`, `memory` and `queue` (see documentation for other properties) just once.
Then, at runtime, Nextflow will use those settings to generate the appropriate backend-specific scripts based on the executor setting.

We'll cover that standardized syntax in the next section.

### Takeaway

You now know how to change the executor to use different kinds of computing infrastructure.

### What's next?

Learn how to evaluate and express resource allocations and limitations in Nextflow.

---

## 3. Allocate compute resources with process directives

As noted above, high-performance computing system generally allow or require you to specify request allocations and set limitations for compute resources such as the number of CPUs and memory to use.

By default, Nextflow will use a single CPU and 2GB of memory for each process.
The corresponding process directives are called `cpus` and `memory`, so the following configuration is implied:

```groovy title="Built-in configuration" linenums="1"
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
nextflow run 3-main.nf -with-report report-config-1.html
```

The report is an html file, which you can download and open in your browser. You can also right click it in the file explorer on the left and click on `Show preview` in order to view it in the training environment.

<!-- TODO: insert images -->

Take a few minutes to look through the report and see if you can identify some opportunities for adjusting resources.
Make sure to click on the tabs that show the utilization results as a percentage of what was allocated.
There is some [documentation](https://www.nextflow.io/docs/latest/reports.html) describing all the available features.

### 3.2. Set resource allocations for all processes

The profiling shows that the processes in our training workflow are very lightweight, so let's reduce the default memory allocation to 1GB per process.

Add the following to your `nextflow.config` file:

```groovy title="nextflow.config" linenums="4"
process {
    memory = 1.GB
}
```

### 3.3. Set resource allocations for an individual process

At the same time, we're going to pretend that the `cowpy` process requires more resources than the others, just so we can demonstrate how to adjust allocations for an individual process.

=== "After"

    ```groovy title="nextflow.config" linenums="4" hl_lines="3-6"
    process {
        memory = 1.GB
        withName: 'cowpy' {
            memory = 2.GB
            cpus = 2
        }
    }
    ```

=== "Before"

    ```groovy title="nextflow.config" linenums="14"
    process {
        memory = 1.GB
    }
    ```

With this configuration, all processes will request 1GB of memory and a single CPU (the implied default), except the `cowpy` process, which will request 2GB and 2 CPUs.

!!! Tip

    If you have a machine with few CPUs and you allocate a high number per process, you might see process calls getting queued behind each other.
    This is because Nextflow ensures we don't request more CPUs than are available.

You could then run the workflow again, supplying a different filename for the profiling report, and compare performance before and after the configuration changes.
You may not notice any real difference since this is such a small workload, but this is the approach you would use to analyze the performance and resource requirements of a real-world workflow.

It is very useful when your processes have different resource requirements. It empowers you to right-size the resource allocations you set up for each process based on actual data, not guesswork.

!!! Tip

    This is just a tiny taster of what you can do to optimize your use of resources.
    Nextflow itself has some really neat [dynamic retry logic](https://training.nextflow.io/basic_training/debugging/#dynamic-resources-allocation) built in to retry jobs that fail due to resource limitations.
    Additionally, the Seqera Platform offers AI-driven tooling for optimizing your resource allocations automatically.

### 3.4. Add resource limits

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

As previously, we can't demonstrate this in action since we don't have access to relevant infrastructure in the training environment.
However, if you were to set the executor to `slurm`, try running the workflow with resource allocations that exceed these limits, then look up the `sbatch` command in the `.command.run` script file (which will be generated even though the run is doomed to fail), you would see that the requests that would get sent to the executor are capped at the values specified by `resourceLimits`.

!!! Tip

    The nf-core project has compiled a [collection of configuration files](https://nf-co.re/configs/) shared by various institutions around the world, covering a wide range of HPC and cloud executors.

    Those shared configs are valuable both for people who work there and can therefore just utilize their institution's configuration out of the box, and as a model for people who are looking to develop a configuration for their own infrastructure.

### Takeaway

You know how to generate a profiling report to assess resource utilization and how to modify resource allocations for all processes and/or for individual processes, as well as set resource limitations for running on HPC.

### What's next?

Learn how to manage workflow parameters.

---

## 4. Manage workflow parameters

So far we've been looking at configuration from the technical point of view of the compute infrastructure.
Now let's consider another aspect of workflow configuration that is very important for reproducibility: the configuration of the workflow parameters.

Currently, our workflow is set up to accept a couple of parameter values via the command-line.
This is fine for a simple workflow with very few parameters that need to be set for a given run.
However, many real-world workflows will have many more parameters that may be run-specific, and putting all of them in the command line would be tedious and error-prone.

### 4.1. Specify default parameter values

It is possible to specify default values in the workflow script itself; for example you may see something like this in the main body of the workflow:

```groovy title="Syntax example"
params.input = 'greetings.csv'
params.character = 'turkey'
```

The same syntax can also be used to store parameter defaults in the `nextflow.config` file.
Let's try that out.

Open the `nextflow.config` file and add the following lines to it:

```groovy title="nextflow.config" linenums="1"
/*
 * Pipeline parameters
 */
params.input = 'greetings.csv'
params.character = 'turkey'
```

<details>
  <summary>Code (full file)</summary>

```groovy title="nextflow.config" linenums="1"
docker.enabled = false
conda.enabled = true

process {
    memory = 1.GB
    withName: 'cowpy' {
        memory = 2.GB
        cpus = 2
    }
}

/*
 * Pipeline parameters
 */
params.input = 'greetings.csv'
params.character = 'turkey'
```

</details>

Now you can run the workflow without specifying the parameters on the command line.

```bash
nextflow run 3-main.nf
```

This will produce the same output, but is more convenient to type, especially when the workflow requires multiple parameters.

<details>
  <summary>Command output</summary>

```console linenums="1"
 N E X T F L O W   ~  version 25.04.3

Launching `3-main.nf` [wise_mahavira] DSL2 - revision: 356df0818d

executor >  local (8)
[2e/d12fcb] sayHello (2)       [100%] 3 of 3 ✔
[a0/5799b6] convertToUpper (3) [100%] 3 of 3 ✔
[db/d3bbb6] collectGreetings   [100%] 1 of 1 ✔
[a9/f75d13] cowpy              [100%] 1 of 1 ✔
```

</details>

The final output file should contain the turkey character saying the greetings.

<details>
  <summary>File contents</summary>

```console title="results/cowpy-COLLECTED-output.txt"
 _________
/ HELLO   \
| BONJOUR |
\ HOLà    /
 ---------
  \                                  ,+*^^*+___+++_
   \                           ,*^^^^              )
    \                       _+*                     ^**+_
     \                    +^       _ _++*+_+++_,         )
              _+^^*+_    (     ,+*^ ^          \+_        )
             {       )  (    ,(    ,_+--+--,      ^)      ^\
            { (\@)    } f   ,(  ,+-^ __*_*_  ^^\_   ^\       )
           {:;-/    (_+*-+^^^^^+*+*<_ _++_)_    )    )      /
          ( /  (    (        ,___    ^*+_+* )   <    <      \
           U _/     )    *--<  ) ^\-----++__)   )    )       )
            (      )  _(^)^^))  )  )\^^^^^))^*+/    /       /
          (      /  (_))_^)) )  )  ))^^^^^))^^^)__/     +^^
         (     ,/    (^))^))  )  ) ))^^^^^^^))^^)       _)
          *+__+*       (_))^)  ) ) ))^^^^^^))^^^^^)____*^
          \             \_)^)_)) ))^^^^^^^^^^))^^^^)
           (_             ^\__^^^^^^^^^^^^))^^^^^^^)
             ^\___            ^\__^^^^^^))^^^^^^^^)\\
                  ^^^^^\uuu/^^\uuu/^^^^\^\^\^\^\^\^\^\
                     ___) >____) >___   ^\_\_\_\_\_\_\)
                    ^^^//\\_^^//\\_^       ^(\_\_\_\)
                      ^^^ ^^ ^^^ ^

```

</details>

You can override those defaults by providing parameter values on the command line, or by providing them through another source of configuration information.

### 4.2. Override defaults with a run-specific config file

You may want to override those defaults without having to either specify parameters on the command line, or modify the original script file.

A clean way to do this is to create a new `nextflow.config` file in a run-specific working directory.

Let's start by creating a new directory:

```bash
mkdir -p tux-run
cd tux-run
```

Then, create a blank configuration file in that directory:

```bash
touch nextflow.config
```

Now open the new file and add the parameters you want to customize:

```groovy title="tux-run/nextflow.config" linenums="1"
params.input = '../greetings.csv'
params.character = 'tux'
```

Note that the path to the input file must reflect the directory structure.

We can now run our pipeline from within our new working directory:

```bash
nextflow run ../3-main.nf
```

This will create a new set of directories under `tux-run/` including `tux-run/work/` and `tux-run/results/`.

<details>
  <summary>Command output</summary>

```console linenums="1"
 N E X T F L O W   ~  version 25.04.3

Launching `../3-main.nf` [trusting_escher] DSL2 - revision: 356df0818d

executor >  local (8)
[59/b66913] sayHello (2)       [100%] 3 of 3 ✔
[ad/f06364] convertToUpper (3) [100%] 3 of 3 ✔
[10/714895] collectGreetings   [100%] 1 of 1 ✔
[88/3ece98] cowpy              [100%] 1 of 1 ✔
```

</details>

In this run, Nextflow combines the `nextflow.config` in our current directory with the `nextflow.config` in the root directory of the pipeline, and thereby overrides the default character (turkey) with the tux character.

The final output file should contain the tux character saying the greetings.

<details>
  <summary>File contents</summary>

```console title="results/cowpy-COLLECTED-output.txt"
 _________
/ HELLO   \
| BONJOUR |
\ HOLà    /
 ---------
   \
    \
        .--.
       |o_o |
       |:_/ |
      //   \ \
     (|     | )
    /'\_   _/`\
    \___)=(___/

```

</details>

That's it!

Make sure to change back to the previous directory before moving to the next section.

```bash
cd ..
```

Now let's look at another useful way to set parameter values.

### 4.3. Specify parameters using a parameter file

Nextflow also allows us to specify parameters via a parameter file in either YAML or JSON format.
This makes it very convenient to manage and distribute alternative sets of default values, for example, as well as run-specific parameter values.

We provide an example YAML parameter file in the current directory, called `test-params.yaml`, which contains a key-value pair for each of the inputs our workflow expects.

<details>
  <summary>File contents</summary>

```yaml title="test-params.yaml" linenums="1"
input: "greetings.csv"
character: "stegosaurus"
```

</details>

To run the workflow with this parameter file, simply add `-params-file <filename>` to the base command.

```bash
nextflow run 3-main.nf -params-file test-params.yaml
```

This should run without error.

<details>
  <summary>Command output</summary>

```console title="Output"
 N E X T F L O W   ~  version 25.04.3

Launching `3-main.nf` [disturbed_sammet] DSL2 - revision: ede9037d02

executor >  local (8)
[f0/35723c] sayHello (2)       | 3 of 3 ✔
[40/3efd1a] convertToUpper (3) | 3 of 3 ✔
[17/e97d32] collectGreetings   | 1 of 1 ✔
[98/c6b57b] cowpy              | 1 of 1 ✔
```

</details>

The final output file should contain the stegosaurus character saying the greetings.

<details>
  <summary>File contents</summary>

```console title="results/cowpy-COLLECTED-output.txt"
_________
/ HELLO   \
| HOLà    |
\ BONJOUR /
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

</details>

Using a parameter file may seem like overkill when you only have a few parameters to specify, but some pipelines expect dozens of parameters.
In those cases, using a parameter file will allow us to provide parameter values at runtime without having to type massive command lines and without modifying the workflow script.
It also makes it easier to distribute sets of parameters to collaborators.

### Takeaway

You know how to manage parameter defaults and override them at runtime using command-line arguments or a parameter file.
There are a few more options but these are the ones you are most likely to encounter.

### What's next?

Learn how to bring it all together by using profiles to switch between alternative configurations more conveniently.

---

## 5. Use profiles to select preset configurations

You may want to switch between alternative settings depending on what computing infrastructure you're using. For example, you might want to develop and run small-scale tests locally on your laptop, then run full-scale workloads on HPC or cloud.

This applies to workflow parameters too: you may have different sets of reference files or groups of settings that you want to swap out depending on the data you're analyzing (e.g. mouse vs human data etc).

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
nextflow run 3-main.nf -profile my_laptop
```

This should run without error and produce the same results as previously.

<details>
  <summary>Command output</summary>

```console
 N E X T F L O W   ~  version 25.04.3

Launching `3-main.nf` [gigantic_brazil] DSL2 - revision: ede9037d02

executor >  local (8)
[58/da9437] sayHello (3)       | 3 of 3 ✔
[35/9cbe77] convertToUpper (2) | 3 of 3 ✔
[67/857d05] collectGreetings   | 1 of 1 ✔
[37/7b51b5] cowpy              | 1 of 1 ✔
```

</details>

As you can see, this allows us to toggle between configurations very conveniently at runtime.

!!! Warning

    The `univ_hpc` profile will not run properly in the training environment since we do not have access to a Slurm scheduler.

If in the future we find other elements of configuration that are always co-occurring with these, we can simply add them to the corresponding profile(s).
We can also create additional profiles if there are other elements of configuration that we want to group together.

### 5.3. Create a test profile

As noted above, profiles are not only for infrastructure configuration.
We can also use them to swap out sets of default values for workflow parameters, or to make it easier for ourselves and for others to try out the workflow without having to gather appropriate input values themselves.

Let's take the example of creating a test profile to make it easy to test the workflow with minimal effort.

The syntax for expressing default values is the same as when writing them into the workflow file itself, except we wrap them in a block named `test`:

```groovy title="Syntax example"
    test {
        params.<parameter1>
        params.<parameter2>
        ...
    }
```

If we add a test profile for our workflow, the `profiles` block becomes:

```groovy title="nextflow.config" linenums="4"
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
    test {
        params.input = 'greetings.csv'
        params.character = 'turtle'
    }
}
```

Just like for technical configuration profiles, you can set up multiple different profiles specifying workflow parameters under any arbitrary name you like.

### 5.4. Run the workflow locally with the test profile

Conveniently, profiles are not mutually exclusive, so we can specify multiple profiles in our command line using the following syntax `-profile <profile1>,<profile2>` (for any number of profiles).

!!! Tip

    If you combine profiles that set values for the same elements of configuration and are described in the same configuration file, Nextflow will resolve the conflict by using whichever value it read in last (_i.e._ whatever comes later in the file).

Let's try adding the test profile to our previous command:

```bash
nextflow run 3-main.nf -profile my_laptop,test
```

This should run without error.

<details>
  <summary>Command output</summary>

```console
 N E X T F L O W   ~  version 25.04.3

Launching `3-main.nf` [gigantic_brazil] DSL2 - revision: ede9037d02

executor >  local (8)
[58/da9437] sayHello (3)       | 3 of 3 ✔
[35/9cbe77] convertToUpper (2) | 3 of 3 ✔
[67/857d05] collectGreetings   | 1 of 1 ✔
[37/7b51b5] cowpy              | 1 of 1 ✔
```

</details>

The final output file should contain the turtle character saying the greetings.

<details>
  <summary>File contents</summary>

```console title="results/cowpy-COLLECTED-output.txt"
 _________
/ BONJOUR \
| HOLà    |
\ HELLO   /
 ---------
    \                                  ___-------___
     \                             _-~~             ~~-_
      \                         _-~                    /~-_
             /^\__/^\         /~  \                   /    \
           /|  O|| O|        /      \_______________/        \
          | |___||__|      /       /                \          \
          |          \    /      /                    \          \
          |   (_______) /______/                        \_________ \
          |         / /         \                      /            \
           \         \^\\         \                  /               \     /
             \         ||           \______________/      _-_       //\__//
               \       ||------_-~~-_ ------------- \ --/~   ~\    || __/
                 ~-----||====/~     |==================|       |/~~~~~
                  (_(__/  ./     /                    \_\      \.
                         (_(___/                         \_____)_)
```

</details>

This means that as long as we distribute any test data files with the workflow code, anyone can quickly try out the workflow without having to supply their own inputs via the command line or a parameter file.

!!! Tip

    You can even point to URLs for larger files that are stored externally.
    Nextflow will download them automatically as long as there is an open connection.

### 5.5. Use `nextflow config` to see the resolved configuration

As noted above, sometimes the same parameter can be set to different values in profiles that you want to combine.
And more generally, there are numerous places where elements of configuration can be stored, and sometimes the same properties can be set to different values in different places.

Nextflow applies a set [order of precedence](https://www.nextflow.io/docs/latest/config.html) to resolve any conflicts, but that can be tricky to determine yourself.
And even if nothing is conflicting, it can be tedious to look up all the possible places where things could be configured.

Fortunately, Nextflow includes a convenient utility tool called `config` that can automate that whole process for you.

The `config` tool will explore all the contents in your current working directory, hoover up any configuration files, and produce the fully resolved configuration that Nextflow would use to run the workflow.
This allows you to find out what settings will be used without having to launch anything.

#### 5.5.1. Resolve the default configuration

Run this command to resolve the configuration that would be applied by default.

```bash
nextflow config
```

<details>
  <summary>Command output</summary>

```groovy
docker {
   enabled = false
}

conda {
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
   input = 'greetings.csv'
   character = 'turkey'
}
```

</details>

#### 5.5.2. Resolve the configuration with specific settings activated

If you provide command-line parameters, e.g. enabling one or more profiles or loading a parameter file, the command will additionally take those into account.

```bash
nextflow config -profile my_laptop,test
```

<details>
  <summary>Command output</summary>

```groovy
docker {
   enabled = true
}

conda {
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
   input = 'greetings.csv'
   character = 'turtle'
}
```

</details>

This gets especially useful for complex projects that involve multiple layers of configuration.

### Takeaway

You know how to use profiles to select a preset configuration at runtime with minimal hassle.
More generally, you know how to configure your workflow executions to suit different compute platforms and enhance the reproducibility of your analyses.

### What's next?

Give yourself a big pat on the back!
You know everything you need to know to get started running and managing Nextflow pipelines.

That concludes this course, but if you're eager to keep learning, we have two main recommendations:

- If you want to dig deeper into developing your own pipelines, have a look at [Hello Nextflow](../hello_nextflow/index.md), a course for beginners that covers the same general progression as this one but goes into much more detail about channels and operators.
- If you would like to continue learning how to run Nextflow pipelines without going deeper into the code, have a look at the first part of [Hello nf-core](../hello_nf-core/index.md), which introduces the tooling for finding and running pipelines from the hugely popular [nf-core](https://nf-co.re/) project.

Have fun!
