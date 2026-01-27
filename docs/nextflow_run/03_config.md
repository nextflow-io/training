# Part 3: Run configuration

This section will explore how to manage the configuration of a Nextflow pipeline in order to customize its behavior, adapt it to different environments, and optimize resource usage _without altering a single line of the workflow code itself_.

There are multiple ways to do this, which can be used in combination and are interpreted according to the order of precedence described [here](https://www.nextflow.io/docs/latest/config.html).

In this part of the course, we are going to show you the simplest and most common configuration file mechanism, the `nextflow.config` file, which you already encountered in the section on containers in Part 2.

We'll go over essential components of Nextflow configuration such as process directives, executors, profiles, and parameter files.
By learning to utilize these configuration options effectively, you can take full advantage of the flexibility, scalability, and performance of Nextflow pipelines.

To exercise these elements of configuration, we're going to be running a fresh copy of the workflow we last ran at the end of Part 2 of this training course, renamed `3-main.nf`.

If you're not familiar with the Hello pipeline or you could use a reminder, see [this info page](../info/hello_pipeline.md).

---

## 1. Manage workflow input parameters

We're going to start with an aspect of configuration that is simply an extension of what we've been working with so far: the management of input parameters.

Currently, our workflow is set up to accept several parameter values via the command-line, declared in a `params` block in the workflow script itself.
One has a default value set as part of its declaration.

However, you might want to set defaults for all of them, or override the existing default without having to either specify parameters on the command line, or modify the original script file.

There are multiple ways of doing that; we're going to show you three basic ways that are very commonly used.

### 1.1. Set up values in `nextflow.config`

This is the simplest approach, though it's possibly the least flexible since the main `nextflow.config` file is not something you want to be editing for every run.
But it does have the advantage of separating the concerns of _declaring_ the parameters in the workflow (which definitely belongs there) versus supplying _default values_, which are more at home in a configuration file.

Let's do this in two steps.

#### 1.1.1. Create a `params` block in the configuration file

Make the following code changes in the `nextflow.config` file:

=== "After"

    ```groovy title="nextflow.config" linenums="1" hl_lines="3-10"
    docker.enabled = true

    /*
    * Pipeline parameters
    */
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }
    ```

=== "Before"

    ```groovy title="nextflow.config" linenums="1"
    docker.enabled = true
    ```

Note that we didn't simply copy the `params` block from the workflow to the configuration file.
For the `batch` parameter that had a default value declared already, the syntax is a little different.
In the workflow file, that's a typed declaration.
In the configuration, those are value assignments.

Technically, this is sufficient for overriding the default values still specified in the workflow file.
You could modify the default value for `batch` and run the workflow to satisfy yourself that the value set in the configuration file overrides the one set in the workflow file.

But in the spirit of moving configuration completely to the configuration file, let's remove that default value from the workflow file entirely.

#### 1.1.2. Remove the default value for `batch` in the workflow file

Make the following code change to the `3-main.nf` workflow file:

=== "After"

    ```groovy title="3-main.nf" linenums="9" hl_lines="6"
    /*
    * Pipeline parameters
    */
    params {
        input: Path
        batch: String
        character: String
    }
    ```

=== "Before"

    ```groovy title="3-main.nf" linenums="9" hl_lines="6"
    /*
    * Pipeline parameters
    */
    params {
        input: Path
        batch: String = 'batch'
        character: String
    }
    ```

Now the workflow file itself does not set any default values for these parameters.

#### 1.1.3. Run the pipeline

Let's test that it works correctly without specifying any parameters in the command line.

```bash
nextflow run 3-main.nf
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `3-main.nf` [disturbed_einstein] DSL2 - revision: ede9037d02

    executor >  local (8)
    [f0/35723c] sayHello (2)       | 3 of 3 ✔
    [40/3efd1a] convertToUpper (3) | 3 of 3 ✔
    [17/e97d32] collectGreetings   | 1 of 1 ✔
    [98/c6b57b] cowpy              | 1 of 1 ✔
    ```

This still produces the same output as previously.

The final ASCII art output is in the `results/3-main/` directory, under the name `cowpy-COLLECTED-batch-output.txt`, same as before.

??? abstract "File contents"

    ```console title="results/3-main/cowpy-COLLECTED-batch-output.txt"
    _________
    / HOLà    \
    | HELLO   |
    \ BONJOUR /
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

Functionally, this move has changed nothing, but conceptually it's a little cleaner to have the default values set in the configuration file.

### 1.2. Use a run-specific configuration file

That's great, but sometimes you might want to run some temporary experiments with different default values without messing with the main configuration file.
You can do that by creating a new `nextflow.config` file in a subdirectory that you'll use as working directory for your experiments.

#### 1.2.1. Create the working directory with a blank configuration

Let's start by creating a new directory and moving into it:

```bash
mkdir -p tux-run
cd tux-run
```

Then, create a blank configuration file in that directory:

```bash
touch nextflow.config
```

This produces an empty file.

#### 1.2.2. Set up the experimental configuration

Now open the new file and add the parameters you want to customize:

```groovy title="tux-run/nextflow.config" linenums="1"
params {
    input = '../data/greetings.csv'
    batch = 'experiment'
    character = 'tux'
}
```

Note that the path to the input file must reflect the directory structure.

#### 1.2.3. Run the pipeline

We can now run our pipeline from within our new working directory.
Make sure to adapt the path accordingly!

```bash
nextflow run ../3-main.nf
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `../3-main.nf` [trusting_escher] DSL2 - revision: 356df0818d

    executor >  local (8)
    [59/b66913] sayHello (2)       [100%] 3 of 3 ✔
    [ad/f06364] convertToUpper (3) [100%] 3 of 3 ✔
    [10/714895] collectGreetings   [100%] 1 of 1 ✔
    [88/3ece98] cowpy              [100%] 1 of 1 ✔
    ```

This will create a new set of directories under `tux-run/` including `tux-run/work/` and `tux-run/results/`.

In this run, Nextflow combines the `nextflow.config` in our current directory with the `nextflow.config` in the root directory of the pipeline, and thereby overrides the default character (turkey) with the tux character.

The final output file should contain the tux character saying the greetings.

??? abstract "File contents"

    ```console title="tux-run/results/3-main/cowpy-COLLECTED-experiment-output.txt"
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

That's it; now you have a space for experimenting without modifying your 'normal' configuration.

!!! warning

    Make sure to change back to the previous directory before moving to the next section!

    ```bash
    cd ..
    ```

Now let's look at another useful way to set parameter values.

### 1.3. Use a parameter file

The subdirectory approach works great for experimenting, but it does involve a bit of setup and requires that you adapt paths accordingly.
There's a simpler approach for when you want to run your pipeline with a specific set of values, or enable someone else to do it with minimal effort.

Nextflow allows us to specify parameters via a parameter file in either YAML or JSON format, which makes it very convenient to manage and distribute alternative sets of default values, for example, as well as run-specific parameter values.

#### 1.3.1. Examine the example parameter file

To demonstrate this, we provide an example parameter file in the current directory, called `test-params.yaml`:

```yaml title="test-params.yaml" linenums="1"
input: "data/greetings.csv"
batch: "yaml"
character: "stegosaurus"
```

This parameter file contains a key-value pair for each of the inputs that we want to specify.
Note the use of colons (`:`) instead of equal signs (`=`) if you compare the syntax to the configuration file.
The config file is written in Groovy, whereas the parameter file is written in YAML.

!!! info

    We also provide a JSON version of the parameter file as an example but we're not going to run with it here.
    Feel free to try that one on your own.

#### 1.3.2. Run the pipeline

To run the workflow with this parameter file, simply add `-params-file <filename>` to the base command.

```bash
nextflow run 3-main.nf -params-file test-params.yaml
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `3-main.nf` [disturbed_sammet] DSL2 - revision: ede9037d02

    executor >  local (8)
    [f0/35723c] sayHello (2)       | 3 of 3 ✔
    [40/3efd1a] convertToUpper (3) | 3 of 3 ✔
    [17/e97d32] collectGreetings   | 1 of 1 ✔
    [98/c6b57b] cowpy              | 1 of 1 ✔
    ```

The final output file should contain the stegosaurus character saying the greetings.

??? abstract "File contents"

    ```console title="results/3-main/cowpy-COLLECTED-yaml-output.txt"
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

Using a parameter file may seem like overkill when you only have a few parameters to specify, but some pipelines expect dozens of parameters.
In those cases, using a parameter file will allow us to provide parameter values at runtime without having to type massive command lines and without modifying the workflow script.

It also makes it easier to distribute sets of parameters to collaborators, or as supporting information for a publication, for example.
This makes your work more reproducible by others.

### Takeaway

You know how to take advantage of key configuration options for managing workflow inputs.

### What's next?

Learn how to manage where and how your workflow outputs get published.

---

## 2. Manage workflow outputs

The workflow we inherited uses paths for workflow-level output declarations, which isn't terribly flexible and involves a lot of repetition.

Let's look at a few common ways you might configure this to be more flexible.

### 2.1. Customize the `outputDir` directory name

Each version of the workflow we've run so far has published its outputs to a different subdirectory hardcoded into the output definitions.

Let's change that to use a user-configurable parameter.
We could create a whole new parameter for this, but let's use the `batch` parameter since it's right there.

#### 2.1.1. Set a value for `outputDir` in the configuration file

The path Nextflow uses for publishing outputs is controlled by the `outputDir` option.
To change the path for all outputs, you can set a value for this option in the `nextflow.config` configuration file.

Add the following code to the `nextflow.config` file:

=== "After"

    ```groovy title="nextflow.config" linenums="9" hl_lines="10-13"
    /*
    * Pipeline parameters
    */
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }

    /*
    * Output settings
    */
    outputDir = "results/${params.batch}"
    ```

=== "Before"

    ```groovy title="nextflow.config" linenums="9"
    /*
    * Pipeline parameters
    */
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }
    ```

This will replace the built-in default path, `results/`, with `results/` plus the value of the `batch` parameter as subdirectory.
You could also change the `results` part if you wanted.

For a temporary change, you could set this option from the command-line using the `-output-dir` parameter in your command (but then you couldn't use the `batch` parameter value).

#### 2.1.2. Remove the repeated part of the hardcoded path

We still have a subdirectory hardcoded in the output options, so let's get rid of that now.

Make the following code changes in the workflow file:

=== "After"

    ```groovy title="3-main.nf" linenums="42" hl_lines="3 7 11 15 19"
    output {
        first_output {
            path 'intermediates'
            mode 'copy'
        }
        uppercased {
            path 'intermediates'
            mode 'copy'
        }
        collected {
            path 'intermediates'
            mode 'copy'
        }
        batch_report {
            path ''
            mode 'copy'
        }
        cowpy_art {
            path ''
            mode 'copy'
        }
    }
    ```

=== "Before"

    ```groovy title="3-main.nf" linenums="42" hl_lines="3 7 11 15 19"
    output {
        first_output {
            path '3-main/intermediates'
            mode 'copy'
        }
        uppercased {
            path '3-main/intermediates'
            mode 'copy'
        }
        collected {
            path '3-main/intermediates'
            mode 'copy'
        }
        batch_report {
            path '3-main'
            mode 'copy'
        }
        cowpy_art {
            path '3-main'
            mode 'copy'
        }
    }
    ```

We could also have just added `${params.batch}` to each path instead of modifying the `outputDir` default, but this is more concise.

#### 2.1.3. Run the pipeline

Let's test that it works correctly, setting the batch name to `outdir` from the command line.

```bash
nextflow run 3-main.nf --batch outdir
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `3-main.nf` [disturbed_einstein] DSL2 - revision: ede9037d02

    executor >  local (8)
    [f0/35723c] sayHello (2)       | 3 of 3 ✔
    [40/3efd1a] convertToUpper (3) | 3 of 3 ✔
    [17/e97d32] collectGreetings   | 1 of 1 ✔
    [98/c6b57b] cowpy              | 1 of 1 ✔
    ```

This still produces the same output as previously, except this time we find our outputs under `results/outdir/`.

??? abstract "Directory contents"

    ```console
    results/outdir/
    ├── cowpy-COLLECTED-outdir-output.txt
    ├── intermediates
    │   ├── Bonjour-output.txt
    │   ├── COLLECTED-outdir-output.txt
    │   ├── Hello-output.txt
    │   ├── Holà-output.txt
    │   ├── UPPER-Bonjour-output.txt
    │   ├── UPPER-Hello-output.txt
    │   └── UPPER-Holà-output.txt
    └── outdir-report.txt
    ```

You can combine this approach with custom path definitions to construct any directory hierarchy you like.

### 2.2. Organize outputs by process

One popular way to organize outputs further is to do it by process, _i.e._ create subdirectories for each process run in the pipeline.

#### 2.2.1. Replace the output paths by a reference to process names

All you need to do is reference the name of the process as `<task>.name` in the output path declaration.

Make the following changes in the workflow file:

=== "After"

    ```groovy title="3-main.nf" linenums="42" hl_lines="3 7 11 15 19"
    output {
        first_output {
            path { sayHello.name }
            mode 'copy'
        }
        uppercased {
            path { convertToUpper.name }
            mode 'copy'
        }
        collected {
            path { collectGreetings.name }
            mode 'copy'
        }
        batch_report {
            path { collectGreetings.name }
            mode 'copy'
        }
        cowpy_art {
            path { cowpy.name }
            mode 'copy'
        }
    }
    ```

=== "Before"

    ```groovy title="3-main.nf" linenums="42" hl_lines="3 7 11 15 19"
    output {
        first_output {
            path 'intermediates'
            mode 'copy'
        }
        uppercased {
            path 'intermediates'
            mode 'copy'
        }
        collected {
            path 'intermediates'
            mode 'copy'
        }
        batch_report {
            path ''
            mode 'copy'
        }
        cowpy_art {
            path ''
            mode 'copy'
        }
    }
    ```

This removes the remaining hardcoded elements from the output path configuration.

#### 2.2.2. Run the pipeline

Let's test that it works correctly, setting the batch name to `pnames` from the command line.

```bash
nextflow run 3-main.nf --batch pnames
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `3-main.nf` [jovial_mcclintock] DSL2 - revision: ede9037d02

    executor >  local (8)
    [f0/35723c] sayHello (2)       | 3 of 3 ✔
    [40/3efd1a] convertToUpper (3) | 3 of 3 ✔
    [17/e97d32] collectGreetings   | 1 of 1 ✔
    [98/c6b57b] cowpy              | 1 of 1 ✔
    ```

This still produces the same output as previously, except this time we find our outputs under `results/pnames/`, and they are grouped by process.

??? abstract "Directory contents"

    ```console
    results/pnames/
    ├── collectGreetings
    │   ├── COLLECTED-pnames-output.txt
    │   └── pnames-report.txt
    ├── convertToUpper
    │   ├── UPPER-Bonjour-output.txt
    │   ├── UPPER-Hello-output.txt
    │   └── UPPER-Holà-output.txt
    ├── cowpy
    │   └── cowpy-COLLECTED-pnames-output.txt
    └── sayHello
        ├── Bonjour-output.txt
        ├── Hello-output.txt
        └── Holà-output.txt
    ```

Note that here we've erased the distinction between `intermediates` versus final outputs being at the top level.
You could of course mix and match these approaches, for example by setting the first output's path as `intermediates/${sayHello.name}`

### 2.3. Set the publish mode at the workflow level

Finally, in the spirit of reducing the amount of repetitive code, we can replace the per-output `mode` declarations with a single line in the configuration.

#### 2.3.1. Add `workflow.output.mode` to the configuration file

Add the following code to the `nextflow.config` file:

=== "After"

    ```groovy title="nextflow.config" linenums="2" hl_lines="5"
    /*
    * Output settings
    */
    outputDir = "results/${params.batch}"
    workflow.output.mode = 'copy'
    ```

=== "Before"

    ```groovy title="nextflow.config" linenums="12"
    /*
    * Output settings
    */
    outputDir = "results/${params.batch}"
    ```

Just like the `outputDir` option, giving `workflow.output.mode` a value in the configuration file would be sufficient to override what is set in the workflow file, but let's remove the unnecessary code anyway.

#### 2.3.2. Remove output mode from the workflow file

Make the following changes in the workflow file:

=== "After"

    ```groovy title="3-main.nf" linenums="42"
    output {
        first_output {
            path { sayHello.name }
        }
        uppercased {
            path { convertToUpper.name }
        }
        collected {
            path { collectGreetings.name }
        }
        batch_report {
            path { collectGreetings.name }
        }
        cowpy_art {
            path { cowpy.name }
        }
    }
    ```

=== "Before"

    ```groovy title="3-main.nf" linenums="42" hl_lines="3 7 11 15 19"
    output {
        first_output {
            path { sayHello.name }
            mode 'copy'
        }
        uppercased {
            path { convertToUpper.name }
            mode 'copy'
        }
        collected {
            path { collectGreetings.name }
            mode 'copy'
        }
        batch_report {
            path { collectGreetings.name }
            mode 'copy'
        }
        cowpy_art {
            path { cowpy.name }
            mode 'copy'
        }
    }
    ```

That's more concise, isn't it?

#### 2.3.3. Run the pipeline

Let's test that it works correctly, setting the batch name to `outmode` from the command line.

```bash
nextflow run 3-main.nf --batch outmode
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `3-main.nf` [rowdy_sagan] DSL2 - revision: ede9037d02

    executor >  local (8)
    [f0/35723c] sayHello (2)       | 3 of 3 ✔
    [40/3efd1a] convertToUpper (3) | 3 of 3 ✔
    [17/e97d32] collectGreetings   | 1 of 1 ✔
    [98/c6b57b] cowpy              | 1 of 1 ✔
    ```

This still produces the same output as previously, except this time we find our outputs under `results/outmode/`.
They are still all proper copies, not symlinks.

??? abstract "Directory contents"

    ```console
    results/outmode/
    ├── collectGreetings
    │   ├── COLLECTED-outmode-output.txt
    │   └── outmode-report.txt
    ├── convertToUpper
    │   ├── UPPER-Bonjour-output.txt
    │   ├── UPPER-Hello-output.txt
    │   └── UPPER-Holà-output.txt
    ├── cowpy
    │   └── cowpy-COLLECTED-outmode-output.txt
    └── sayHello
        ├── Bonjour-output.txt
        ├── Hello-output.txt
        └── Holà-output.txt
    ```

The main reason you might still want to use the per-output way of setting mode is if you want to mix and match within the same workflow, _i.e._ have some outputs be copied and some be symlinked.

There are plenty of other options that you can customize in this way, but hopefully this gives you a sense of the range of options and how to utilize them effectively to suit your preferences.

### Takeaway

You know how to control the naming and structure of the directories where your outputs are published, as well as the workflow output publishing mode.

### What's next?

Learn how to adapt your workflow configuration to your compute environment, starting with the software packaging technology.

---

## 3. Select a software packaging technology

So far we've been looking at configuration elements that control how inputs go in and where inputs come out. Now it's time to focus more specifically on adapting your workflow configuration to your compute environment.

The first step on that path is specifying where the software packages that will get run in each step are going to be coming from.
Are they already installed in the local compute environment?
Do we need to retrieve images and run them via a container system?
Or do we need to retrieve Conda packages and build a local Conda environment?

In the very first part of this training course (Parts 1-4) we just used locally installed software in our workflow.
Then in Part 5, we introduced Docker containers and the `nextflow.config` file, which we used to enable the use of Docker containers.

Now let's see how we can configure an alternative software packaging option via the `nextflow.config` file.

### 3.1. Disable Docker and enable Conda in the config file

Let's pretend we're working on an HPC cluster and the admin doesn't allow the use of Docker for security reasons.
Fortunately for us, Nextflow supports multiple other container technologies such as including Singularity (which is more widely used on HPC), and software package managers such as Conda.

We can change our configuration file to use Conda instead of Docker.
To do so, let's switch the value of `docker.enabled` to `false`, and add a directive enabling the use of Conda:

=== "After"

    ```groovy title="nextflow.config" linenums="1" hl_lines="1-2"
    docker.enabled = false
    conda.enabled = true
    ```

=== "Before"

    ```groovy title="nextflow.config" linenums="1" hl_lines="1"
    docker.enabled = true
    ```

This will allow Nextflow to create and utilize Conda environments for processes that have Conda packages specified.
Which means we now need to add one of those to our `cowpy` process!

### 3.2. Specify a Conda package in the process definition

We've already retrieved the URI for a Conda package containing the `cowpy` tool: `conda-forge::cowpy==1.1.5`

Now we add the URI to the `cowpy` process definition using the `conda` directive:

=== "After"

    ```groovy title="modules/cowpy.nf" linenums="4" hl_lines="4"
    process cowpy {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
        conda 'conda-forge::cowpy==1.1.5'

        input:
    ```

=== "Before"

    ```groovy title="modules/cowpy.nf" linenums="4"
    process cowpy {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

        input:
    ```

To be clear, we're not _replacing_ the `docker` directive, we're _adding_ an alternative option.

!!! tip

    There are a few different ways to get the URI for a given conda package.
    We recommend using the [Seqera Containers](https://seqera.io/containers/) search query, which will give you a URI that you can copy and paste, even if you're not planning to create a container from it.

### 3.3. Run the workflow to verify that it can use Conda

Let's try it out.

```bash
nextflow run 3-main.nf --batch conda
```

??? success "Command output"

    ```console title="Output"
    N E X T F L O W   ~  version 25.10.2

    Launching `3-main.nf` [trusting_lovelace] DSL2 - revision: 028a841db1

    executor >  local (8)
    [ee/4ca1f2] sayHello (3)       | 3 of 3 ✔
    [20/2596a7] convertToUpper (1) | 3 of 3 ✔
    [b3/e15de5] collectGreetings   | 1 of 1 ✔
    [c5/af5f88] cowpy              | 1 of 1 ✔
    ```

This should work without issue and produce the same outputs as previously under `results/conda`.

Behind the scenes, Nextflow has retrieved the Conda packages and created the environment, which normally takes a bit of work; so it's nice that we don't have to do any of that ourselves!

!!! info

    This runs quickly because the `cowpy` package is quite small, but if you're working with large packages, it may take a bit longer than usual the first time, and you might see the console output stay 'stuck' for a minute or so before completing.
    This is normal and is due to the extra work Nextflow does the first time you use a new package.

From our standpoint, it looks like it works exactly the same as running with Docker, even though on the backend the mechanics are a bit different.

This means we're all set to run with Conda environments if needed.

??? info "Mixing and matching Docker and Conda"

    Since these directives are assigned per process, it is possible 'mix and match', _i.e._ configure some of the processes in your workflow to run with Docker and others with Conda, for example, if the compute infrastructure you are using supports both.
    In that case, you would enable both Docker and Conda in your configuration file.
    If both are available for a given process, Nextflow will prioritize containers.

    And as noted earlier, Nextflow supports multiple other software packaging and container technologies, so you are not limited to just those two.

### Takeaway

You know how to configure which software package each process should use, and how to switch between technologies.

### What's next?

Learn how to change the execution platform used by Nextflow to actually do the work.

---

## 4. Select an execution platform

Until now, we've been running our pipeline with the local executor.
This executes each task on the machine that Nextflow is running on.
When Nextflow begins, it looks at the available CPUs and memory.
If the resources of the tasks ready to run exceed the available resources, Nextflow will hold the last tasks back from execution until one or more of the earlier tasks have finished, freeing up the necessary resources.

The local executor is convenient and efficient, but it is limited to that single machine. For very large workloads, you may discover that your local machine is a bottleneck, either because you have a single task that requires more resources than you have available, or because you have so many tasks that waiting for a single machine to run them would take too long.

Nextflow supports [many different execution backends](https://www.nextflow.io/docs/latest/executor.html), including HPC schedulers (Slurm, LSF, SGE, PBS, Moab, OAR, Bridge, HTCondor and others) as well as cloud execution backends such (AWS Batch, Google Cloud Batch, Azure Batch, Kubernetes and more).

### 4.1. Targeting a different backend

The choice of executor is set by a process directive called `executor`.
By default it is set to `local`, so the following configuration is implied:

```groovy title="Built-in configuration"
process {
    executor = 'local'
}
```

To set the executor to target a different backend, you would simply specify the executor you want using similar syntax as described above for resource allocations (see [documentation](https://www.nextflow.io/docs/latest/executor.html) for all options).

```groovy title="nextflow.config"
process {
    executor = 'slurm'
}
```

!!! warning

    We can't actually test this in the training environment because it's not set up to connect to an HPC.

### 4.2. Dealing with backend-specific syntax for execution parameters

Most high-performance computing platforms allow (and sometimes require) that you specify certain parameters such as resource allocation requests and limitations (for e.g. number of CPUs and memory) and name of the job queue to use.

Unfortunately, each of these systems uses different technologies, syntaxes and configurations for defining how a job should be defined and submitted to the relevant scheduler.

??? abstract "Examples"

    For example, the same job requiring 8 CPUs and 4GB of RAM to be executed on the queue "my-science-work" needs to be expressed in different following ways depending on the backend.

    ```bash title="Config for SLURM / submit using sbatch"
    #SBATCH -o /path/to/my/task/directory/my-task-1.log
    #SBATCH --no-requeue
    #SBATCH -c 8
    #SBATCH --mem 4096M
    #SBATCH -p my-science-work
    ```

    ```bash title="Config for PBS / submit using qsub"
    #PBS -o /path/to/my/task/directory/my-task-1.log
    #PBS -j oe
    #PBS -q my-science-work
    #PBS -l nodes=1:ppn=5
    #PBS -l mem=4gb
    ```

    ```bash title="Config for SGE / submit using qsub"
    #$ -o /path/to/my/task/directory/my-task-1.log
    #$ -j y
    #$ -terse
    #$ -notify
    #$ -q my-science-work
    #$ -l slots=5
    #$ -l h_rss=4096M,mem_free=4096M
    ```

Fortunately, Nextflow simplifies all of this.
It provides a standardized syntax so that you can specify the relevant properties such as `cpus`, `memory` and `queue` (see documentation for other properties) just once.
Then, at runtime, Nextflow will use those settings to generate the appropriate backend-specific scripts based on the executor setting.

We'll cover that standardized syntax in the next section.

### Takeaway

You now know how to change the executor to use different kinds of computing infrastructure.

### What's next?

Learn how to evaluate and express resource allocations and limitations in Nextflow.

---

## 5. Control compute resource allocations

Most high-performance computing platforms allow (and sometimes require) that you specify certain resource allocation parameters such as number of CPUs and memory.

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

### 5.1. Run the workflow to generate a resource utilization report

If you don't know up front how much CPU and memory your processes are likely to need, you can do some resource profiling, meaning you run the workflow with some default allocations, record how much each process used, and from there, estimate how to adjust the base allocations.

Conveniently, Nextflow includes built-in tools for doing this, and will happily generate a report for you on request.

To do so, add `-with-report <filename>.html` to your command line.

```bash
nextflow run 3-main.nf -with-report report-config-1.html
```

The report is an html file, which you can download and open in your browser. You can also right click it in the file explorer on the left and click on `Show preview` in order to view it in the training environment.

Take a few minutes to look through the report and see if you can identify some opportunities for adjusting resources.
Make sure to click on the tabs that show the utilization results as a percentage of what was allocated.
There is some [documentation](https://www.nextflow.io/docs/latest/reports.html) describing all the available features.

### 5.2. Set resource allocations for all processes

The profiling shows that the processes in our training workflow are very lightweight, so let's reduce the default memory allocation to 1GB per process.

Add the following to your `nextflow.config` file, before the pipeline parameters section:

```groovy title="nextflow.config" linenums="4"
/*
* Process settings
*/
process {
    memory = 1.GB
}
```

That will help reduce the amount of compute we consume.

### 5.3. Set resource allocations for a specific process

At the same time, we're going to pretend that the `cowpy` process requires more resources than the others, just so we can demonstrate how to adjust allocations for an individual process.

=== "After"

    ```groovy title="nextflow.config" linenums="4" hl_lines="6-9"
    /*
    * Process settings
    */
    process {
        memory = 1.GB
        withName: 'cowpy' {
            memory = 2.GB
            cpus = 2
        }
    }
    ```

=== "Before"

    ```groovy title="nextflow.config" linenums="4"
    /*
    * Process settings
    */
    process {
        memory = 1.GB
    }
    ```

With this configuration, all processes will request 1GB of memory and a single CPU (the implied default), except the `cowpy` process, which will request 2GB and 2 CPUs.

!!! info

    If you have a machine with few CPUs and you allocate a high number per process, you might see process calls getting queued behind each other.
    This is because Nextflow ensures we don't request more CPUs than are available.

### 5.4. Run the workflow with the updated configuration

Let's try that out, supplying a different filename for the profiling report so we can compare performance before and after the configuration changes.

```bash
nextflow run 3-main.nf -with-report report-config-2.html
```

You will probably not notice any real difference since this is such a small workload, but this is the approach you would use to analyze the performance and resource requirements of a real-world workflow.

It is very useful when your processes have different resource requirements. It empowers you to right-size the resource allocations you set up for each process based on actual data, not guesswork.

!!! tip

    This is just a tiny taster of what you can do to optimize your use of resources.
    Nextflow itself has some really neat [dynamic retry logic](https://training.nextflow.io/basic_training/debugging/#dynamic-resources-allocation) built in to retry jobs that fail due to resource limitations.
    Additionally, the Seqera Platform offers AI-driven tooling for optimizing your resource allocations automatically as well.

### 5.5. Add resource limits

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

??? info "Institutional reference configurations"

    The nf-core project has compiled a [collection of configuration files](https://nf-co.re/configs/) shared by various institutions around the world, covering a wide range of HPC and cloud executors.

    Those shared configs are valuable both for people who work there and can therefore just utilize their institution's configuration out of the box, and as a model for people who are looking to develop a configuration for their own infrastructure.

### Takeaway

You know how to generate a profiling report to assess resource utilization and how to modify resource allocations for all processes and/or for individual processes, as well as set resource limitations for running on HPC.

### What's next?

Learn how to set up preset configuration profiles and switch between them at runtime.

---

## 6. Use profiles to switch between preset configurations

We've shown you a number of ways that you can customize your pipeline configuration depending on the project you're working on or the compute environment you're using.

You may want to switch between alternative settings depending on what computing infrastructure you're using. For example, you might want to develop and run small-scale tests locally on your laptop, then run full-scale workloads on HPC or cloud.

Nextflow lets you set up any number of profiles that describe different configurations, which you can then select at runtime using a command-line argument, rather than having to modify the configuration file itself.

### 6.1. Create profiles for switching between local development and execution on HPC

Let's set up two alternative profiles; one for running small scale loads on a regular computer, where we'll use Docker containers, and one for running on a university HPC with a Slurm scheduler, where we'll use Conda packages.

#### 6.1.1. Set up the profiles

Add the following to your `nextflow.config` file, after the pipeline parameters section but before the output settings:

```groovy title="nextflow.config" linenums="24"
/*
* Profiles
*/
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

#### 6.1.2. Run the workflow with a profile

To specify a profile in our Nextflow command line, we use the `-profile` argument.

Let's try running the workflow with the `my_laptop` configuration.

```bash
nextflow run 3-main.nf -profile my_laptop
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `3-main.nf` [gigantic_brazil] DSL2 - revision: ede9037d02

    executor >  local (8)
    [58/da9437] sayHello (3)       | 3 of 3 ✔
    [35/9cbe77] convertToUpper (2) | 3 of 3 ✔
    [67/857d05] collectGreetings   | 1 of 1 ✔
    [37/7b51b5] cowpy              | 1 of 1 ✔
    ```

As you can see, this allows us to toggle between configurations very conveniently at runtime.

!!! warning

    The `univ_hpc` profile will not run properly in the training environment since we do not have access to a Slurm scheduler.

If in the future we find other elements of configuration that are always co-occurring with these, we can simply add them to the corresponding profile(s).
We can also create additional profiles if there are other elements of configuration that we want to group together.

### 6.2. Create a profile of test parameters

Profiles are not only for infrastructure configuration.
We can also use them to set default values for workflow parameters, to make it easier for others to try out the workflow without having to gather appropriate input values themselves.
You can consider this an alternative to using a parameter file.

#### 6.2.1. Set up the profile

The syntax for expressing default values in this context looks like this, for a profile that we name `test`:

```groovy title="Syntax example"
    test {
        params.<parameter1>
        params.<parameter2>
        ...
    }
```

If we add a test profile for our workflow, the `profiles` block becomes:

```groovy title="nextflow.config" linenums="24"
/*
* Profiles
*/
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
        params.input = 'data/greetings.csv'
        params.batch = 'test'
        params.character = 'dragonandcow'
    }
}
```

Just like for technical configuration profiles, you can set up multiple different profiles specifying parameters under any arbitrary name you like.

#### 6.2.2. Run the workflow locally with the test profile

Conveniently, profiles are not mutually exclusive, so we can specify multiple profiles in our command line using the following syntax `-profile <profile1>,<profile2>` (for any number of profiles).

If you combine profiles that set values for the same elements of configuration and are described in the same configuration file, Nextflow will resolve the conflict by using whichever value it read in last (_i.e._ whatever comes later in the file).
If the conflicting settings are set in different configuration sources, the default [order of precedence](https://www.nextflow.io/docs/latest/config.html) applies.

Let's try adding the test profile to our previous command:

```bash
nextflow run 3-main.nf -profile my_laptop,test
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `3-main.nf` [jovial_coulomb] DSL2 - revision: 46a6763141

    executor >  local (8)
    [9b/687cdc] sayHello (2)       | 3 of 3 ✔
    [ca/552187] convertToUpper (3) | 3 of 3 ✔
    [e8/83e306] collectGreetings   | 1 of 1 ✔
    [fd/e84fa9] cowpy              | 1 of 1 ✔
    ```

This will use Docker where possible and produce outputs under `results/test`, and this time the character is the comedic duo `dragonandcow`.

??? abstract "File contents"

    ```console title="results/test/"
     _________
    / HOLà    \
    | HELLO   |
    \ BONJOUR /
    ---------
                \                    ^    /^
                  \                  / \  // \
                  \   |\___/|      /   \//  .\
                    \  /O  O  \__  /    //  | \ \           *----*
                      /     /  \/_/    //   |  \  \          \   |
                      \@___\@`    \/_   //    |   \   \         \/\ \
                    0/0/|       \/_ //     |    \    \         \ \
                0/0/0/0/|        \///      |     \     \       | |
              0/0/0/0/0/_|_ /   (  //       |      \     _\     |  /
          0/0/0/0/0/0/`/,_ _ _/  ) ; -.    |    _ _\.-~       /   /
                      ,-}        _      *-.|.-~-.           .~    ~
      \     \__/        `/\      /                 ~-. _ .-~      /
      \____(oo)           *.   }            {                   /
      (    (--)          .----~-.\        \-`                 .~
      //__\\  \__ Ack!   ///.----..<        \             _ -~
      //    \\               ///-._ _ _ _ _ _ _{^ - - - - ~
    ```

This means that as long as we distribute any test data files with the workflow code, anyone can quickly try out the workflow without having to supply their own inputs via the command line or a parameter file.

!!! tip

    We can point to URLs for larger files that are stored externally.
    Nextflow will download them automatically as long as there is an open connection.

    For more details, see the Side Quest [Working with Files](../side_quests/working_with_files.md)

### 6.3. Use `nextflow config` to see the resolved configuration

As noted above, sometimes the same parameter can be set to different values in profiles that you want to combine.
And more generally, there are numerous places where elements of configuration can be stored, and sometimes the same properties can be set to different values in different places.

Nextflow applies a set [order of precedence](https://www.nextflow.io/docs/latest/config.html) to resolve any conflicts, but that can be tricky to determine yourself.
And even if nothing is conflicting, it can be tedious to look up all the possible places where things could be configured.

Fortunately, Nextflow includes a convenient utility tool called `config` that can automate that whole process for you.

The `config` tool will explore all the contents in your current working directory, hoover up any configuration files, and produce the fully resolved configuration that Nextflow would use to run the workflow.
This allows you to find out what settings will be used without having to launch anything.

#### 6.3.1. Resolve the default configuration

Run this command to resolve the configuration that would be applied by default.

```bash
nextflow config
```

??? success "Command output"

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
      batch = 'batch'
      character = 'turkey'
    }
    ```

This shows you the base configuration you get if you don't specify anything extra in the command line.

#### 6.3.2. Resolve the configuration with specific settings activated

If you provide command-line parameters, e.g. enabling one or more profiles or loading a parameter file, the command will additionally take those into account.

```bash
nextflow config -profile my_laptop,test
```

??? success "Command output"

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
      batch = 'test'
      character = 'dragonandcow'
    }
    ```

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

---

## Quiz

<quiz>
When parameter values are set in both the workflow file and nextflow.config, which takes precedence?
- [ ] The workflow file value
- [x] The configuration file value
- [ ] The first value encountered
- [ ] It causes an error
</quiz>

<quiz>
What is the syntax difference between setting a parameter default in a workflow file vs. a config file?
- [ ] They use the same syntax
- [x] Workflow uses typed declaration (param: Type = value), config uses assignment (param = value)
- [ ] Config uses typed declaration, workflow uses assignment
- [ ] Only config files can set default values
</quiz>

<quiz>
How do you specify a parameter file when running a workflow?
- [ ] --params params.yaml
- [ ] -config params.yaml
- [x] -params-file params.yaml
- [ ] --input-params params.yaml
</quiz>

<quiz>
What does the `outputDir` configuration option control?
- [ ] The location of the work directory
- [x] The base path where workflow outputs are published
- [ ] The directory for log files
- [ ] The location of module files
</quiz>

<quiz>
How do you reference a process name dynamically in output path configuration?
- [ ] ${processName}
- [ ] process.name
- [x] { processName.name }
- [ ] @processName
</quiz>

<quiz>
If both Docker and Conda are enabled and a process has both directives, which is prioritized?
- [x] Docker (containers)
- [ ] Conda
- [ ] The first one defined in the process
- [ ] It causes an error
</quiz>

<quiz>
What is the default executor in Nextflow?
- [x] local
- [ ] slurm
- [ ] kubernetes
- [ ] aws
</quiz>

<quiz>
What command generates a resource utilization report?
- [ ] nextflow run workflow.nf -with-metrics
- [ ] nextflow run workflow.nf -with-stats
- [x] nextflow run workflow.nf -with-report report.html
- [ ] nextflow run workflow.nf -profile report
</quiz>

<quiz>
How do you set resource requirements for a specific process named 'cowpy' in the config file?
- [ ] cowpy.memory = '2.GB'
- [ ] process.cowpy.memory = '2.GB'
- [x] process { withName: 'cowpy' { memory = '2.GB' } }
- [ ] resources.cowpy.memory = '2.GB'
</quiz>

<quiz>
What does the `resourceLimits` directive do?
- [ ] Sets minimum resource requirements
- [ ] Allocates resources to processes
- [x] Caps the maximum resources that can be requested
- [ ] Monitors resource usage in real-time
</quiz>

<quiz>
How do you specify multiple profiles in a single command?
- [ ] -profile profile1 -profile profile2
- [ ] -profiles profile1,profile2
- [x] -profile profile1,profile2
- [ ] --profile profile1 --profile profile2
</quiz>

<quiz>
What command shows the fully resolved configuration that Nextflow would use?
- [ ] nextflow show-config
- [ ] nextflow settings
- [x] nextflow config
- [ ] nextflow resolve
</quiz>

<quiz>
What can profiles be used for? (Select all that apply)
- [x] Defining infrastructure-specific settings (executors, containers)
- [x] Setting resource limits for different environments
- [x] Providing test parameters for easy workflow testing
- [ ] Defining new processes
</quiz>
