# Part 1: Manage inputs and outputs

In [Nextflow Run](../nextflow_run/index.md), you saw that a pipeline's default parameter values can be set in `nextflow.config` instead of the workflow file, and that `-with-report` generates an execution summary.
This part builds on that foundation: you'll learn two more ways to supply parameter values, then take full control of where and how your pipeline's outputs get published.

---

## 1. Manage workflow input parameters

Setting parameter defaults in `nextflow.config` works well for values that rarely change.
For values that change often, per run or per collaborator, two other mechanisms are more convenient: a run-specific configuration file, and a parameter file.

### 1.1. Use a run-specific configuration file

Nextflow automatically picks up any `nextflow.config` file present in the current working directory, in addition to the one in the pipeline's own directory.
You can take advantage of this to create a dedicated subdirectory for experimenting with alternative settings, without touching your main configuration.

Create a new directory and an empty configuration file inside it.

```bash
mkdir -p tux-run
touch tux-run/nextflow.config
```

Add the parameters you want to override.

```groovy title="tux-run/nextflow.config" linenums="1"
params {
    input = '../data/greetings.csv'
    batch = 'experiment'
    character = 'tux'
}
```

The input path is relative to the new subdirectory, so it has to point back up one level.

Run the pipeline from inside `tux-run/`, pointing at the pipeline script in the parent directory.

```bash
cd tux-run
nextflow run ../main.nf
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `../main.nf` [trusting_escher] DSL2 - revision: 356df0818d

    executor >  local (8)
    [51/d6a7ea] sayHello (3)       | 3 of 3 ✔
    [af/5e1684] convertToUpper (1) | 3 of 3 ✔
    [8f/33bdde] collectGreetings   | 1 of 1 ✔
    [59/05e97e] cowpy              | 1 of 1 ✔

    Outputs:

      /workspaces/training/nextflow-config/tux-run/results

      first_output:
        - full_pipeline/intermediates/Bonjour-output.txt
        - full_pipeline/intermediates/Hello-output.txt
        - full_pipeline/intermediates/Hola-output.txt

      uppercased:
        - full_pipeline/intermediates/UPPER-Bonjour-output.txt
        - full_pipeline/intermediates/UPPER-Hello-output.txt
        - full_pipeline/intermediates/UPPER-Hola-output.txt

      collected: full_pipeline/intermediates/COLLECTED-experiment-output.txt

      batch_report: full_pipeline/experiment-report.txt

      cowpy_art: full_pipeline/cowpy-COLLECTED-experiment-output.txt
    ```

Nextflow merges the `nextflow.config` in your current directory with the one in the pipeline's directory, so the tux character overrides the turkey default.

??? abstract "File contents"

    ```console title="tux-run/results/full_pipeline/cowpy-COLLECTED-experiment-output.txt"
     _________
    / BONJOUR \
    | HOLA    |
    \ HELLO   /
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

This gives you a disposable space for experimenting without touching your normal configuration.

!!! warning

    Change back to the parent directory before continuing.

    ```bash
    cd ..
    ```

### 1.2. Use a parameter file

For sharing an exact set of parameter values with a collaborator, or recording them for a publication, Nextflow supports [parameter files](https://nextflow.io/docs/latest/config.html#parameter-file) in YAML or JSON format.

A parameter file called `test-params.yaml` is already provided in your working directory.

```yaml title="test-params.yaml" linenums="1"
input: "data/greetings.csv"
batch: "yaml"
character: "stegosaurus"
```

The syntax uses colons (`:`) instead of the equal signs (`=`) used in `nextflow.config`, since this file is plain YAML rather than Groovy.

!!! info

    A JSON version, `test-params.json`, is also provided. Feel free to try it on your own; the syntax for passing it is identical.

Pass the file with `-params-file`.

```bash
nextflow run main.nf -params-file test-params.yaml
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [disturbed_sammet] DSL2 - revision: 6e18cd130e

    executor >  local (8)
    [0c/0e490e] sayHello (2)       | 3 of 3 ✔
    [58/755545] convertToUpper (2) | 3 of 3 ✔
    [47/1c1484] collectGreetings   | 1 of 1 ✔
    [65/15d8dd] cowpy              | 1 of 1 ✔

    Outputs:

      /workspaces/training/nextflow-config/results

      first_output:
        - full_pipeline/intermediates/Bonjour-output.txt
        - full_pipeline/intermediates/Hello-output.txt
        - full_pipeline/intermediates/Hola-output.txt

      uppercased:
        - full_pipeline/intermediates/UPPER-Bonjour-output.txt
        - full_pipeline/intermediates/UPPER-Hola-output.txt
        - full_pipeline/intermediates/UPPER-Hello-output.txt

      collected: full_pipeline/intermediates/COLLECTED-yaml-output.txt

      batch_report: full_pipeline/yaml-report.txt

      cowpy_art: full_pipeline/cowpy-COLLECTED-yaml-output.txt
    ```

??? abstract "File contents"

    ```console title="results/full_pipeline/cowpy-COLLECTED-yaml-output.txt"
     _________
    / HOLA    \
    | BONJOUR |
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

A parameter file becomes especially valuable once a pipeline has more than a handful of parameters: it lets you supply them all at once, without a sprawling command line or any change to the workflow script, and it's easy to distribute alongside your results.

### Takeaway

You know three ways to manage workflow input parameters: defaults in `nextflow.config`, a run-specific configuration file for disposable experiments, and a parameter file for sharing exact, reproducible parameter sets.

### What's next?

Learn how to control where and how your workflow outputs get published.

---

## 2. Manage workflow outputs

So far, the pipeline's `output` block has hardcoded a `full_pipeline` subdirectory into every output path.
That's not very flexible, and it repeats itself five times.
Let's make it more dynamic.

### 2.1. Customize the output directory

The path Nextflow uses as the base for publishing outputs is controlled by the `outputDir` option.
Set it in `nextflow.config`, using the `batch` parameter to keep runs separate.

=== "After"

    ```groovy title="nextflow.config" linenums="20" hl_lines="7-10"
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }

    /*
     * Output settings
     */
    outputDir = "results_config/${params.batch}"
    ```

=== "Before"

    ```groovy title="nextflow.config" linenums="20"
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }
    ```

Now remove the hardcoded `full_pipeline` prefix from the workflow's output paths, since `outputDir` already covers it.

=== "After"

    ```groovy title="main.nf" linenums="37" hl_lines="3 7 11 15 19"
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

    ```groovy title="main.nf" linenums="37" hl_lines="3 7 11 15 19"
    output {
        first_output {
            path 'full_pipeline/intermediates'
            mode 'copy'
        }
        uppercased {
            path 'full_pipeline/intermediates'
            mode 'copy'
        }
        collected {
            path 'full_pipeline/intermediates'
            mode 'copy'
        }
        batch_report {
            path 'full_pipeline'
            mode 'copy'
        }
        cowpy_art {
            path 'full_pipeline'
            mode 'copy'
        }
    }
    ```

Run it, setting the batch name to `outdir` so it's easy to spot in the output path.

```bash
nextflow run main.nf --batch outdir
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [amazing_church] DSL2 - revision: 6e18cd130e

    executor >  local (8)
    [21/591b73] sayHello (2)       | 3 of 3 ✔
    [9c/4e98d3] convertToUpper (2) | 3 of 3 ✔
    [f2/5c5f24] collectGreetings   | 1 of 1 ✔
    [16/d4f61f] cowpy              | 1 of 1 ✔

    Outputs:

      /workspaces/training/nextflow-config/results_config/outdir

      first_output:
        - intermediates/Hello-output.txt
        - intermediates/Hola-output.txt
        - intermediates/Bonjour-output.txt

      uppercased:
        - intermediates/UPPER-Bonjour-output.txt
        - intermediates/UPPER-Hello-output.txt
        - intermediates/UPPER-Hola-output.txt

      collected: intermediates/COLLECTED-outdir-output.txt

      batch_report: outdir-report.txt

      cowpy_art: cowpy-COLLECTED-outdir-output.txt
    ```

The outputs now land under `results_config/outdir/` instead of the built-in `results/` default.

??? abstract "Directory contents"

    ```console
    results_config/outdir
    ├── cowpy-COLLECTED-outdir-output.txt
    ├── intermediates
    │   ├── Bonjour-output.txt
    │   ├── COLLECTED-outdir-output.txt
    │   ├── Hello-output.txt
    │   ├── Hola-output.txt
    │   ├── UPPER-Bonjour-output.txt
    │   ├── UPPER-Hello-output.txt
    │   └── UPPER-Hola-output.txt
    └── outdir-report.txt
    ```

You can combine `outputDir` with custom path definitions to build any directory hierarchy you like.

### 2.2. Organize outputs by process

A common way to organize outputs further is by process: a subdirectory per process that ran.
Reference the process name as `<process>.name` in each output path.

=== "After"

    ```groovy title="main.nf" linenums="37" hl_lines="3 7 11 15 19"
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

    ```groovy title="main.nf" linenums="37" hl_lines="3 7 11 15 19"
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

This removes the last hardcoded strings from the output path configuration.

```bash
nextflow run main.nf --batch pnames
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [jovial_mcclintock] DSL2 - revision: 6e18cd130e

    executor >  local (8)
    [f3/2a731e] sayHello (1)       | 3 of 3 ✔
    [16/fe06f6] convertToUpper (3) | 3 of 3 ✔
    [1d/2d23bd] collectGreetings   | 1 of 1 ✔
    [20/f466ca] cowpy              | 1 of 1 ✔

    Outputs:

      /workspaces/training/nextflow-config/results_config/pnames

      first_output:
        - sayHello/Bonjour-output.txt
        - sayHello/Hola-output.txt
        - sayHello/Hello-output.txt

      uppercased:
        - convertToUpper/UPPER-Hello-output.txt
        - convertToUpper/UPPER-Hola-output.txt
        - convertToUpper/UPPER-Bonjour-output.txt

      collected: collectGreetings/COLLECTED-pnames-output.txt

      batch_report: collectGreetings/pnames-report.txt

      cowpy_art: cowpy/cowpy-COLLECTED-pnames-output.txt
    ```

The outputs are now grouped by process under `results_config/pnames/`.

??? abstract "Directory contents"

    ```console
    results_config/pnames
    ├── collectGreetings
    │   ├── COLLECTED-pnames-output.txt
    │   └── pnames-report.txt
    ├── convertToUpper
    │   ├── UPPER-Bonjour-output.txt
    │   ├── UPPER-Hello-output.txt
    │   └── UPPER-Hola-output.txt
    ├── cowpy
    │   └── cowpy-COLLECTED-pnames-output.txt
    └── sayHello
        ├── Bonjour-output.txt
        ├── Hello-output.txt
        └── Hola-output.txt
    ```

!!! note

    You can mix and match these approaches, and even combine multiple variables, for example `#!groovy "${params.batch}/intermediates/${sayHello.name}"`.

### 2.3. Set the publish mode for the whole workflow

Finally, the repeated `mode 'copy'` line in every output block can be replaced with a single setting in the configuration file.

=== "After"

    ```groovy title="nextflow.config" linenums="26" hl_lines="5"
    /*
     * Output settings
     */
    outputDir = "results_config/${params.batch}"
    workflow.output.mode = 'copy'
    ```

=== "Before"

    ```groovy title="nextflow.config" linenums="26"
    /*
     * Output settings
     */
    outputDir = "results_config/${params.batch}"
    ```

=== "After"

    ```groovy title="main.nf" linenums="37"
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

    ```groovy title="main.nf" linenums="37" hl_lines="3 7 11 15 19"
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

```bash
nextflow run main.nf --batch outmode
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [rowdy_sagan] DSL2 - revision: 6e18cd130e

    executor >  local (8)
    [e5/4ab319] sayHello (1)       | 3 of 3 ✔
    [18/a06d83] convertToUpper (2) | 3 of 3 ✔
    [01/60bc8f] collectGreetings   | 1 of 1 ✔
    [87/77dedd] cowpy              | 1 of 1 ✔

    Outputs:

      /workspaces/training/nextflow-config/results_config/outmode

      first_output:
        - sayHello/Hola-output.txt
        - sayHello/Hello-output.txt
        - sayHello/Bonjour-output.txt

      uppercased:
        - convertToUpper/UPPER-Hola-output.txt
        - convertToUpper/UPPER-Bonjour-output.txt
        - convertToUpper/UPPER-Hello-output.txt

      collected: collectGreetings/COLLECTED-outmode-output.txt

      batch_report: collectGreetings/outmode-report.txt

      cowpy_art: cowpy/cowpy-COLLECTED-outmode-output.txt
    ```

The outputs land in the same structure as before, still as real copies rather than symlinks.
The main reason to keep the per-output `mode` setting instead is if you want to mix and match within the same workflow, copying some outputs and symlinking others.

### Takeaway

You know how to control the base output directory, organize outputs by process, and set the publish mode once for the whole workflow instead of repeating it for every output.

### What's next?

Head on to [Part 2](./02_packaging_execution_resources.md), where you'll learn how to adapt your configuration to different compute environments.

---

## Summary

In this part you learned to:

- Manage input parameters with a run-specific configuration file and a parameter file
- Customize the base output directory with `outputDir`
- Organize outputs by process
- Set the publish mode once for the whole workflow with `workflow.output.mode`
