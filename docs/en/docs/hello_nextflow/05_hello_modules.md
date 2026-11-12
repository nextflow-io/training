# Part 5: Hello Modules

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/43Ot-f0iOME?si=y8lAedhEHWaTV4zd&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } See [the whole playlist](https://youtube.com/playlist?list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&si=eF7cLR62goy-lc6n) on the Nextflow YouTube channel.

:green_book: The video transcript is available [here](./transcripts/04_hello_modules.md).
///

This section covers how to organize your workflow code to make development and maintenance of your pipeline more efficient and sustainable.
Specifically, we are going to demonstrate how to use [**modules**](https://nextflow.io/docs/latest/modules.html).

In Nextflow, a **module** is a standalone code file, often encapsulating a single process definition.
To use a module in a workflow, you just add a single-line `include` statement to your workflow code file; then you can integrate the process into the workflow the same way you normally would.
That makes it possible to reuse process definitions in multiple workflows without producing multiple copies of the code.

When we started developing our workflow, we wrote everything in one single code file.
Now we're going to give the project the standard layout of a Nextflow pipeline: the workflow file becomes `main.nf`, the processes go into individual modules, and the parameter values go into the configuration file.

This will make our code more shareable, flexible and maintainable.

This part reorganizes the code of the pipeline you built in Parts 1 to 4, without changing what it does:

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-course-code.svg"
</figure>

??? info "How to begin from this section"

    This section of the course assumes you have completed Parts 1-4 of the [Hello Nextflow](./index.md) course and have a complete working pipeline.

    If you're starting the course from this point, you'll need to copy the `nextflow.config` file over from the solutions, since it enables Docker for the `cowpy` process:

    ```bash
    cp solutions/4-hello-containers/nextflow.config .
    ```

---

## 0. Warmup: Run `hello-modules.nf`

We're going to use the workflow script `hello-modules.nf` as a starting point.
It is equivalent to the script produced by working through Part 4 of this training course.

Just to make sure everything is working, run the script once before making any changes:

```bash
nextflow run hello-modules.nf
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `hello-modules.nf` [suspicious_mahavira] revision: a94f0f7888

    executor >  local (8)
    [4c/5eb228] sayHello (Celine)                  | 3 of 3 ✔
    [c1/c55d1e] translateGreeting (es-Raphael.txt) | 3 of 3 ✔
    [06/430547] collectGreetings (batch)           | 1 of 1 ✔
    [50/bbafd2] cowpy (turkey)                     | 1 of 1 ✔

    Outputs:

      /workspaces/training/hello-nextflow/results

      greetings:
        - batch/intermediates/Mateo.txt
        - batch/intermediates/Raphael.txt
        - batch/intermediates/Celine.txt

      translated:
        - batch/intermediates/es-Mateo.txt
        - batch/intermediates/es-Celine.txt
        - batch/intermediates/es-Raphael.txt

      collected: batch/intermediates/batch-es-greetings.txt

      batch_report: batch/batch-es-report.txt

      cowpy_art: batch/cowpy-batch-es-greetings.txt
    ```

As previously, you will find the output files under `results/batch/`.

??? abstract "Directory contents"

    ```console
    results/batch
    ├── batch-es-report.txt
    ├── cowpy-batch-es-greetings.txt
    └── intermediates
        ├── batch-es-greetings.txt
        ├── Celine.txt
        ├── es-Celine.txt
        ├── es-Mateo.txt
        ├── es-Raphael.txt
        ├── Mateo.txt
        └── Raphael.txt
    ```

If that worked for you, you're ready to organize your pipeline code.

---

## 1. Organize the code into modules

The first step is to give the code the standard file layout: an entry workflow file called `main.nf`, and one module file per process in a `modules/` directory.

### 1.1. Rename the workflow file to `main.nf`

By convention, the workflow file that holds the entry workflow of a Nextflow project is called `main.nf`.
It is the entry script of the finished Hello pipeline in Nextflow Run, and it is the file Nextflow looks for when you run a pipeline from a repository without naming a script.

Rename `hello-modules.nf` to `main.nf`:

```bash
mv hello-modules.nf main.nf
```

From now on, run the workflow with `nextflow run main.nf`.
Run it now with the `-resume` flag:

```bash
nextflow run main.nf -resume
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [awesome_marconi] revision: a94f0f7888

    [4c/5eb228] sayHello (Celine)                 | 3 of 3, cached: 3 ✔
    [63/eb68b8] translateGreeting (es-Celine.txt) | 3 of 3, cached: 3 ✔
    [06/430547] collectGreetings (batch)          | 1 of 1, cached: 1 ✔
    [50/bbafd2] cowpy (turkey)                    | 1 of 1, cached: 1 ✔

    Outputs:

      /workspaces/training/hello-nextflow/results

      greetings:
        - batch/intermediates/Raphael.txt
        - batch/intermediates/Mateo.txt
        - batch/intermediates/Celine.txt

      translated:
        - batch/intermediates/es-Mateo.txt
        - batch/intermediates/es-Raphael.txt
        - batch/intermediates/es-Celine.txt

      collected: batch/intermediates/batch-es-greetings.txt

      batch_report: batch/batch-es-report.txt

      cowpy_art: batch/cowpy-batch-es-greetings.txt
    ```

Every task is retrieved from the cache: the name of the workflow file doesn't change the work to be done.

### 1.2. Create a directory to store modules

It is best practice to store your modules in a specific directory.
You can call that directory anything you want, but the convention is to call it `modules/`.

```bash
mkdir modules
```

### 1.3. Create a module for `sayHello()`

In its simplest form, turning an existing process into a module is little more than a copy-paste operation.
We're going to create a file stub for the module, copy the relevant code over then delete it from the main workflow file.

Then all we'll need to do is add an `include` statement so that Nextflow will know to pull in the relevant code at runtime.

#### 1.3.1. Create a file stub for the new module

Create an empty file for the module called `sayHello.nf`.

```bash
touch modules/sayHello.nf
```

This gives us a place to put the process code.

#### 1.3.2. Move the `sayHello` process code to the module file

Copy the whole process definition over from the workflow file to the module file.

```groovy title="modules/sayHello.nf" linenums="1"
/*
 * Use echo to print a greeting to a file
 */
process sayHello {

    tag "${name}"

    input:
    val name

    output:
    path "${name}.txt"

    script:
    """
    echo 'Hello, ${name}!' > '${name}.txt'
    """
}
```

Once that is done, delete the process definition from the workflow file.

#### 1.3.3. Add an include declaration before the `params` block

The syntax for including a process from a module is fairly straightforward:

```groovy title="Syntax: include declaration"
include { <PROCESS_NAME> } from '<path_to_module>'
```

Insert that above the `params` block and fill it out appropriately.

=== "After"

    ```groovy title="main.nf" linenums="90" hl_lines="1 2"
    // Include modules
    include { sayHello } from './modules/sayHello.nf'

    /*
     * Pipeline parameters
     */
    params {
        input: Path = 'data/names.csv'
        language: String = 'es'
        batch: String = 'batch'
        character: String = 'turkey'
    }
    ```

=== "Before"

    ```groovy title="main.nf" linenums="90"
    /*
     * Pipeline parameters
     */
    params {
        input: Path = 'data/names.csv'
        language: String = 'es'
        batch: String = 'batch'
        character: String = 'turkey'
    }
    ```

You see we've filled in the process name, `sayHello`, and the path to the file containing the module code, `./modules/sayHello.nf`.

Once all four processes are in modules, the workflow file will include each module file this way:

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/modules.svg"
</figure>

#### 1.3.4. Run the workflow

We're running the workflow with essentially the same code and inputs as before, so run it with the `-resume` flag and see what happens.

```bash
nextflow run main.nf -resume
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [nasty_heisenberg] revision: 89ffc7e3d2

    [4c/5eb228] sayHello (Celine)                  | 3 of 3, cached: 3 ✔
    [c1/c55d1e] translateGreeting (es-Raphael.txt) | 3 of 3, cached: 3 ✔
    [06/430547] collectGreetings (batch)           | 1 of 1, cached: 1 ✔
    [50/bbafd2] cowpy (turkey)                     | 1 of 1, cached: 1 ✔

    Outputs:

      /workspaces/training/hello-nextflow/results

      greetings:
        - batch/intermediates/Mateo.txt
        - batch/intermediates/Celine.txt
        - batch/intermediates/Raphael.txt

      translated:
        - batch/intermediates/es-Mateo.txt
        - batch/intermediates/es-Celine.txt
        - batch/intermediates/es-Raphael.txt

      collected: batch/intermediates/batch-es-greetings.txt

      batch_report: batch/batch-es-report.txt

      cowpy_art: batch/cowpy-batch-es-greetings.txt
    ```

This should run very quickly because everything is cached.
Feel free to check the published outputs.

Nextflow recognized that it's still all the same work to be done, even if the code is split up into multiple files.

### 1.4. Modularize the remaining processes

The `translateGreeting`, `collectGreetings` and `cowpy` processes follow exactly the same pattern as `sayHello`.

#### 1.4.1. Move the process code into module files

Create an empty file for each of the three modules:

```bash
touch modules/translateGreeting.nf modules/collectGreetings.nf modules/cowpy.nf
```

Copy each process definition over from the workflow file to its module file:

```groovy title="modules/translateGreeting.nf" linenums="1"
/*
 * Use Python to translate the greeting into another language
 */
process translateGreeting {

    tag "${language}-${input_file}"

    input:
    path input_file
    val language

    output:
    path "${language}-${input_file}"

    script:
    """
    #!/usr/bin/env python3
    import sys

    translations = {
        'es': 'Hola',        # Spanish
        'fr': 'Bonjour',     # French
        'it': 'Ciao',        # Italian
        'nl': 'Hallo',       # Dutch
        'ar': 'Marhaban',    # Arabic
        'ja': 'Konnichiwa',  # Japanese
        'en': 'Hello',       # English
    }

    if '${language}' not in translations:
        sys.exit('Unsupported language code: ${language}')

    with open('${input_file}') as greeting_file:
        greeting = greeting_file.read()

    with open('${language}-${input_file}', 'w') as translated_file:
        translated_file.write(greeting.replace('Hello', translations['${language}']))
    """
}
```

```groovy title="modules/collectGreetings.nf" linenums="1"
/*
 * Collect translated greetings into a single output file
 */
process collectGreetings {

    tag "${batch_name}"

    input:
    path input_files
    val batch_name
    val language

    output:
    path "${batch_name}-${language}-greetings.txt", emit: outfile
    path "${batch_name}-${language}-report.txt", emit: report

    script:
    count_greetings = input_files.size()
    """
    cat ${input_files} > '${batch_name}-${language}-greetings.txt'
    echo 'Number of greetings processed: ${count_greetings}' > '${batch_name}-${language}-report.txt'
    """
}
```

```groovy title="modules/cowpy.nf" linenums="1"
/*
 * Generate ASCII art with cowpy (https://github.com/jeffbuttars/cowpy)
 */
process cowpy {

    tag "${character}"

    container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
    conda 'conda-forge::cowpy==1.1.5'

    input:
    path input_file
    val character

    output:
    path "cowpy-${input_file}"

    script:
    """
    cat ${input_file} | cowpy -c "${character}" > cowpy-${input_file}
    """
}
```

A module can hold a script block in any language: the Python code of `translateGreeting` moves over exactly as it was.
The `container` and `conda` directives of `cowpy` move with the process, so the module carries its own software requirements.

Once that is done, delete the three process definitions from the workflow file.

#### 1.4.2. Add the include declarations

Add an include declaration for each of the new modules, under the one for `sayHello`:

=== "After"

    ```groovy title="main.nf" linenums="3" hl_lines="3-5"
    // Include modules
    include { sayHello } from './modules/sayHello.nf'
    include { translateGreeting } from './modules/translateGreeting.nf'
    include { collectGreetings } from './modules/collectGreetings.nf'
    include { cowpy } from './modules/cowpy.nf'
    ```

=== "Before"

    ```groovy title="main.nf" linenums="3"
    // Include modules
    include { sayHello } from './modules/sayHello.nf'
    ```

The workflow file now contains only the include declarations, the `params` block, the `workflow` block and the `output` block.

!!! tip "Skip the copy-pasting"

    If you're short on time, this command brings your directory to the end state of sections 1.4.1 and 1.4.2, by copying the four module files and the updated workflow file from the solutions:

    ```bash
    cp solutions/5-hello-modules/modules/*.nf modules/ && cp solutions/5-hello-modules/main-2.nf main.nf
    ```

#### 1.4.3. Run the workflow

Run the workflow with the `-resume` flag:

```bash
nextflow run main.nf -resume
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [special_ptolemy] revision: ceb815abde

    [4c/5eb228] sayHello (Celine)                  | 3 of 3, cached: 3 ✔
    [c1/c55d1e] translateGreeting (es-Raphael.txt) | 3 of 3, cached: 3 ✔
    [06/430547] collectGreetings (batch)           | 1 of 1, cached: 1 ✔
    [50/bbafd2] cowpy (turkey)                     | 1 of 1, cached: 1 ✔

    Outputs:

      /workspaces/training/hello-nextflow/results

      greetings:
        - batch/intermediates/Mateo.txt
        - batch/intermediates/Celine.txt
        - batch/intermediates/Raphael.txt

      translated:
        - batch/intermediates/es-Celine.txt
        - batch/intermediates/es-Mateo.txt
        - batch/intermediates/es-Raphael.txt

      collected: batch/intermediates/batch-es-greetings.txt

      batch_report: batch/batch-es-report.txt

      cowpy_art: batch/cowpy-batch-es-greetings.txt
    ```

Every task is retrieved from the cache once again, including `cowpy`, which still runs in its container.

### Takeaway

You know that the entry workflow file of a Nextflow project is conventionally called `main.nf`, and how to move processes into local modules without breaking the resumability of the workflow.

If you write another pipeline that calls on one of those processes, you just need to type one short `include` statement to use the relevant module.
This is better than copy-pasting the code, because if later you decide to improve the module, all your pipelines will inherit the improvements.

When you add a new step to a pipeline organized this way, you can write it as a module from the start: create the module file, include it in `main.nf`, and call it in the workflow.

### What's next?

Move the parameter values out of the workflow script.

---

## 2. Organize the configuration

With the code in modules, the next step is to organize the configuration: move the default parameter values into `nextflow.config`, put the test values in a `test` profile, and add a profile for running with Conda.

### 2.1. Move the default values to the configuration

The workflow script still sets default values for its parameters in the `params` block:

```groovy title="main.nf" linenums="9"
/*
 * Pipeline parameters
 */
params {
    input: Path = 'data/names.csv'
    language: String = 'es'
    batch: String = 'batch'
    character: String = 'turkey'
}
```

That was convenient while we were developing the workflow, but these values are not all the same kind of thing:

- `language = 'es'`, `batch = 'batch'` and `character = 'turkey'` are true defaults: sensible values for anyone who runs the pipeline without specifying them.
- `input = 'data/names.csv'` is a test value: it points to the example data provided with the course, which nobody would want to process by default in a real analysis.

True defaults belong in the `params` section of `nextflow.config`.
The test value doesn't belong among them; you'll give it its own place in section 2.2.
The workflow script keeps the parameter declarations and their types, so it is still easy to see at a glance which parameters the pipeline accepts.

#### 2.1.1. Remove the values from the workflow script

In the `params` block of the workflow script, make the following code change:

=== "After"

    ```groovy title="main.nf" linenums="9" hl_lines="5-8"
    /*
     * Pipeline parameters
     */
    params {
        input: Path
        language: String
        batch: String
        character: String
    }
    ```

=== "Before"

    ```groovy title="main.nf" linenums="9"
    /*
     * Pipeline parameters
     */
    params {
        input: Path = 'data/names.csv'
        language: String = 'es'
        batch: String = 'batch'
        character: String = 'turkey'
    }
    ```

#### 2.1.2. Add the default values to `nextflow.config`

Open `nextflow.config` and make the following change:

=== "After"

    ```groovy title="nextflow.config" linenums="1" hl_lines="3-10"
    docker.enabled = true

    /*
     * Pipeline parameters
     */
    params {
        language = 'es'
        batch = 'batch'
        character = 'turkey'
    }
    ```

=== "Before"

    ```groovy title="nextflow.config" linenums="1"
    docker.enabled = true
    ```

The `params` section of the configuration file uses `=` to assign values, without types: the types are declared in the workflow script.
There is no default value for `input`.

#### 2.1.3. Run the workflow without an input file

Try running the workflow the same way as before:

```bash
nextflow run main.nf -resume
```

??? failure "Command output"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [cheeky_euclid] revision: f1abb39e4d

    Parameter `input` is required but was not specified on the command line, params file, or config
    ```

Since `input` no longer has a default value, Nextflow stops before running anything and tells you that the parameter is missing.
`language`, `batch` and `character` did not cause an error, because their default values now come from `nextflow.config`.

The input has no default value on purpose.
No input file makes sense for every analysis, so whoever runs the pipeline has to say which file to process.

#### 2.1.4. Provide the input file on the command line

Provide the input file with `--input`:

```bash
nextflow run main.nf --input data/names.csv -resume
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [lethal_kare] revision: f1abb39e4d

    [f2/e0daf2] sayHello (Raphael)                 | 3 of 3, cached: 3 ✔
    [c1/c55d1e] translateGreeting (es-Raphael.txt) | 3 of 3, cached: 3 ✔
    [06/430547] collectGreetings (batch)           | 1 of 1, cached: 1 ✔
    [50/bbafd2] cowpy (turkey)                     | 1 of 1, cached: 1 ✔

    Outputs:

      /workspaces/training/hello-nextflow/results

      greetings:
        - batch/intermediates/Mateo.txt
        - batch/intermediates/Celine.txt
        - batch/intermediates/Raphael.txt

      translated:
        - batch/intermediates/es-Mateo.txt
        - batch/intermediates/es-Celine.txt
        - batch/intermediates/es-Raphael.txt

      collected: batch/intermediates/batch-es-greetings.txt

      batch_report: batch/batch-es-report.txt

      cowpy_art: batch/cowpy-batch-es-greetings.txt
    ```

The parameter values are the same as before, so every task is retrieved from the cache, and the outputs are still published under `results/batch/`.

### 2.2. Create a test profile

Typing the path to the input file every time gets tedious when you run the pipeline repeatedly to test it.
You may also want test runs to be easy to tell apart from real ones.
A `test` profile solves both: it bundles a test input file with any other values you want to use for testing, and you activate it with a single option.

#### 2.2.1. Add the profile to `nextflow.config`

In `nextflow.config`, make the following change:

=== "After"

    ```groovy title="nextflow.config" linenums="1" hl_lines="12-21"
    docker.enabled = true

    /*
     * Pipeline parameters
     */
    params {
        language = 'es'
        batch = 'batch'
        character = 'turkey'
    }

    /*
     * Profiles
     */
    profiles {
        test {
            params.input = 'data/names.csv'
            params.batch = 'test'
            params.character = 'tux'
        }
    }
    ```

=== "Before"

    ```groovy title="nextflow.config" linenums="1"
    docker.enabled = true

    /*
     * Pipeline parameters
     */
    params {
        language = 'es'
        batch = 'batch'
        character = 'turkey'
    }
    ```

The `test` profile sets `input` to the example data file.
It also overrides two default values: it names the batch `test`, so test outputs are published to `results/test/`, and it has the ASCII art drawn by the `tux` character.
When you activate the profile, its values take precedence over the ones in the `params` section.

#### 2.2.2. Run the workflow with the test profile

Run the workflow with `-profile test`:

```bash
nextflow run main.nf -profile test -resume
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [loving_brenner] revision: f1abb39e4d

    executor >  local (2)
    [9c/c67843] sayHello (Mateo)                   | 3 of 3, cached: 3 ✔
    [c1/c55d1e] translateGreeting (es-Raphael.txt) | 3 of 3, cached: 3 ✔
    [99/dc4374] collectGreetings (test)            | 1 of 1 ✔
    [c5/2440e4] cowpy (tux)                        | 1 of 1 ✔

    Outputs:

      /workspaces/training/hello-nextflow/results

      greetings:
        - test/intermediates/Celine.txt
        - test/intermediates/Mateo.txt
        - test/intermediates/Raphael.txt

      translated:
        - test/intermediates/es-Mateo.txt
        - test/intermediates/es-Celine.txt
        - test/intermediates/es-Raphael.txt

      collected: test/intermediates/test-es-greetings.txt

      batch_report: test/test-es-report.txt

      cowpy_art: test/cowpy-test-es-greetings.txt
    ```

The `sayHello` and `translateGreeting` tasks are retrieved from the cache, because the profile doesn't change their inputs.
`collectGreetings` and `cowpy` run again, because the batch name and the character changed.

The outputs are published under `results/test/`, separate from the outputs of your earlier runs in `results/batch/`, and the ASCII art is drawn by `tux`:

```console title="results/test/cowpy-test-es-greetings.txt"
 ________________
/ Hola, Celine!  \
| Hola, Raphael! |
\ Hola, Mateo!   /
 ----------------
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

The run shows how the profile combines with the defaults: it supplies the test input, overrides `batch` and `character`, and leaves `language` at its default value, `es`.

From now on, run the workflow with `-profile test`, or provide `--input` on the command line.

### 2.3. Add a `conda` profile

In Part 4, you gave the `cowpy` process both a container and a Conda package, and enabled Docker in `nextflow.config`.
The finished Hello pipeline in Nextflow Run is set up the same way: it uses Docker by default, and provides a `conda` profile to use Conda instead.

#### 2.3.1. Add the profile to `nextflow.config`

In the `profiles` section of `nextflow.config`, make the following change:

=== "After"

    ```groovy title="nextflow.config" linenums="15" hl_lines="7-10"
    profiles {
        test {
            params.input = 'data/names.csv'
            params.batch = 'test'
            params.character = 'tux'
        }
        conda {
            docker.enabled = false
            conda.enabled = true
        }
    }
    ```

=== "Before"

    ```groovy title="nextflow.config" linenums="15"
    profiles {
        test {
            params.input = 'data/names.csv'
            params.batch = 'test'
            params.character = 'tux'
        }
    }
    ```

The `conda` profile turns Docker off and turns Conda on.
Docker stays enabled by default, since the settings in a profile only apply when you activate it.

#### 2.3.2. Run the workflow with the `conda` profile

Activate both the `test` and `conda` profiles, separating their names with a comma:

```bash
nextflow run main.nf -profile test,conda -resume
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [reverent_church] revision: f1abb39e4d

    executor >  local (1)
    [9c/c67843] sayHello (Mateo)                   | 3 of 3, cached: 3 ✔
    [c1/c55d1e] translateGreeting (es-Raphael.txt) | 3 of 3, cached: 3 ✔
    [99/dc4374] collectGreetings (test)            | 1 of 1, cached: 1 ✔
    [1d/6735b3] cowpy (tux)                        | 1 of 1 ✔
    Creating env using conda: conda-forge::cowpy==1.1.5 [cache /workspaces/training/hello-nextflow/work/conda/env-898314d566668b6587ad714ae06b8520]

    Outputs:

      /workspaces/training/hello-nextflow/results

      greetings:
        - test/intermediates/Raphael.txt
        - test/intermediates/Celine.txt
        - test/intermediates/Mateo.txt

      translated:
        - test/intermediates/es-Celine.txt
        - test/intermediates/es-Mateo.txt
        - test/intermediates/es-Raphael.txt

      collected: test/intermediates/test-es-greetings.txt

      batch_report: test/test-es-report.txt

      cowpy_art: test/cowpy-test-es-greetings.txt
    ```

Nextflow creates a Conda environment from the package declared in the `conda` directive of `cowpy`, and runs the task in that environment.
The `cowpy` task runs again because its software environment changed.
The other processes don't declare any software, so their tasks are retrieved from the cache.

When you run the workflow without the `conda` profile, it uses Docker again.

### Takeaway

You know how to separate the parameter declarations in the workflow script from their values in the configuration.
True default values go in the `params` section of `nextflow.config`, and the test input, along with any values that should differ for test runs, goes in a `test` profile.
You also know how to add a profile that switches the pipeline from Docker to Conda.

Your project now has the standard layout of a Nextflow pipeline, the same as the finished Hello pipeline in [Nextflow Run](../nextflow_run/index.md): `main.nf` holds the entry workflow, `modules/` holds one module file per process, and `nextflow.config` holds the default parameter values and the `test` and `conda` profiles.
Anyone who knows how a Nextflow pipeline is organized can find their way around it.

### What's next?

Take a break!

You've completed the Hello Nextflow course.
Head to the [course summary](next_steps.md) to review what you learned, share your feedback, and choose where to go next.

---

## Bonus exercises

These optional exercises let you practice what you learned in this part.
Each one can be solved with what this part covered; open the hint if you get stuck.
The solutions are also in `solutions/5-hello-modules/bonus/`.

### Exercise 1: Add a profile for French greetings

Starting from your `main.nf` and `nextflow.config` at the end of this part, add a profile named `french` that sets the language to French.
Then run the workflow with both the `test` and `french` profiles, and check that the collected greetings are in French.

??? tip "Hint"

    Write it like the `test` profile, setting `params.language` instead of `params.input`.
    To activate two profiles, separate their names with a comma.

??? solution

    In `nextflow.config`, add the `french` profile next to the `test` profile:

    ```groovy title="nextflow.config" linenums="15" hl_lines="7-9"
    profiles {
        test {
            params.input = 'data/names.csv'
            params.batch = 'test'
            params.character = 'tux'
        }
        french {
            params.language = 'fr'
        }
        conda {
            docker.enabled = false
            conda.enabled = true
        }
    }
    ```

    Then run the workflow with both profiles:

    ```bash
    nextflow run main.nf -profile test,french -resume
    ```

    ```console title="Command output"
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [mighty_mercator] revision: f1abb39e4d

    executor >  local (5)
    [f2/e0daf2] sayHello (Raphael)               | 3 of 3, cached: 3 ✔
    [17/8eac23] translateGreeting (fr-Mateo.txt) | 3 of 3 ✔
    [fa/04dc6b] collectGreetings (test)          | 1 of 1 ✔
    [c8/397e0b] cowpy (tux)                      | 1 of 1 ✔

    Outputs:

      /workspaces/training/hello-nextflow/results

      greetings:
        - test/intermediates/Celine.txt
        - test/intermediates/Raphael.txt
        - test/intermediates/Mateo.txt

      translated:
        - test/intermediates/fr-Raphael.txt
        - test/intermediates/fr-Celine.txt
        - test/intermediates/fr-Mateo.txt

      collected: test/intermediates/test-fr-greetings.txt

      batch_report: test/test-fr-report.txt

      cowpy_art: test/cowpy-test-fr-greetings.txt
    ```

    The `sayHello` tasks are cached, and the translation runs again in French, followed by the steps that depend on it:

    ```console title="results/test/intermediates/test-fr-greetings.txt"
    Bonjour, Celine!
    Bonjour, Raphael!
    Bonjour, Mateo!
    ```

    Since the collected file, the report and the ASCII art are named after the language, the French results are published next to the Spanish ones instead of overwriting them:

    ```console title="Directory contents"
    results/test
    ├── cowpy-test-es-greetings.txt
    ├── cowpy-test-fr-greetings.txt
    ├── intermediates
    │   ├── Celine.txt
    │   ├── es-Celine.txt
    │   ├── es-Mateo.txt
    │   ├── es-Raphael.txt
    │   ├── fr-Celine.txt
    │   ├── fr-Mateo.txt
    │   ├── fr-Raphael.txt
    │   ├── Mateo.txt
    │   ├── Raphael.txt
    │   ├── test-es-greetings.txt
    │   └── test-fr-greetings.txt
    ├── test-es-report.txt
    └── test-fr-report.txt
    ```

### Exercise 2: Give each module its own folder

Larger projects often store each module in its own folder, as `modules/<name>/main.nf`, so that files related to the module can sit next to it.
Starting from your `main.nf` at the end of this part, move each module file into a folder of its own, with this layout, and update the workflow so it runs as before:

```console
modules
├── collectGreetings
│   └── main.nf
├── cowpy
│   └── main.nf
├── sayHello
│   └── main.nf
└── translateGreeting
    └── main.nf
```

??? tip "Hint"

    Nextflow finds each module from the path in its `include` declaration.

??? solution

    After moving the files, update the `include` declarations in `main.nf`:

    ```groovy title="main.nf" linenums="3" hl_lines="2-5"
    // Include modules
    include { sayHello } from './modules/sayHello/main.nf'
    include { translateGreeting } from './modules/translateGreeting/main.nf'
    include { collectGreetings } from './modules/collectGreetings/main.nf'
    include { cowpy } from './modules/cowpy/main.nf'
    ```

    Then run the workflow:

    ```bash
    nextflow run main.nf -profile test -resume
    ```

    ```console title="Command output"
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [evil_hypatia] revision: b5835a4414

    [9c/c67843] sayHello (Mateo)                 | 3 of 3, cached: 3 ✔
    [01/58d58b] translateGreeting (es-Mateo.txt) | 3 of 3, cached: 3 ✔
    [99/dc4374] collectGreetings (test)          | 1 of 1, cached: 1 ✔
    [c5/2440e4] cowpy (tux)                      | 1 of 1, cached: 1 ✔

    Outputs:

      /workspaces/training/hello-nextflow/results

      greetings:
        - test/intermediates/Raphael.txt
        - test/intermediates/Mateo.txt
        - test/intermediates/Celine.txt

      translated:
        - test/intermediates/es-Celine.txt
        - test/intermediates/es-Raphael.txt
        - test/intermediates/es-Mateo.txt

      collected: test/intermediates/test-es-greetings.txt

      batch_report: test/test-es-report.txt

      cowpy_art: test/cowpy-test-es-greetings.txt
    ```

    Every task is retrieved from the cache: moving a module to another file doesn't change the work it does.

---

## Quiz

<quiz>
What is the conventional name of the file that holds the entry workflow of a Nextflow project?
- [ ] `nextflow.nf`
- [ ] `workflow.nf`
- [x] `main.nf`
- [ ] `pipeline.nf`

Learn more: [1.1. Rename the workflow file to `main.nf`](#11-rename-the-workflow-file-to-mainnf)
</quiz>

<quiz>
What is a module in Nextflow?
- [ ] A workflow definition
- [ ] A configuration file
- [ ] A channel operator
- [x] A standalone file that can contain process definitions

Learn more: [1.3. Create a module for `sayHello()`](#13-create-a-module-for-sayhello)
</quiz>

<quiz>
What convention is typically used for storing module files?
- [x] In a `modules/` directory
- [ ] In a `bin/` directory
- [ ] In the same directory as the workflow
- [ ] In a `lib/` directory

Learn more: [1.2. Create a directory to store modules](#12-create-a-directory-to-store-modules)
</quiz>

<quiz>
What is the correct syntax to use a module?

- [ ] `#!groovy require { sayHello } from './modules/sayHello.nf'`
- [ ] `#!groovy import { sayHello } from './modules/sayHello.nf'`
- [ ] `#!groovy load { sayHello } from './modules/sayHello.nf'`
- [x] `#!groovy include { sayHello } from './modules/sayHello.nf'`

Learn more: [1.3.3. Add an include declaration before the `params` block](#133-add-an-include-declaration-before-the-params-block)
</quiz>

<quiz>
What do you need to change in a process whose script is written in Python when you move it into a module?
- [ ] Declare the module as a Python module
- [x] Nothing: the process moves over exactly as it was
- [ ] Add the Python interpreter to the `include` statement
- [ ] Rewrite the script in Bash

Learn more: [1.4.1. Move the process code into module files](#141-move-the-process-code-into-module-files)
</quiz>

<quiz>
What happens to the `-resume` functionality when using modules?
- [ ] It only works for local modules
- [ ] It no longer works
- [ ] It requires additional configuration
- [x] It works the same as before

Learn more: [1.3.4. Run the workflow](#134-run-the-workflow)
</quiz>

<quiz>
What are the benefits of using modules? (Select all that apply)
- [ ] Faster execution speed
- [x] Better organization of workflow code
- [x] Code reusability across workflows
- [x] Easier maintenance
</quiz>

<quiz>
Where should the path to the example input data used for testing go?
- [ ] In the `params` block of the workflow script, as a default value
- [x] In a `test` profile in `nextflow.config`
- [ ] In the `params` section of `nextflow.config`
- [ ] In the module that uses it

Learn more: [2.2. Create a test profile](#22-create-a-test-profile)
</quiz>

<quiz>
What does the `conda` profile in `nextflow.config` do when you activate it?
- [ ] It converts the containers into Conda packages
- [ ] It runs every process twice, once with Docker and once with Conda
- [x] It disables Docker and enables Conda, so processes use their `conda` packages instead of their containers
- [ ] It installs Conda in your environment

Learn more: [2.3. Add a `conda` profile](#23-add-a-conda-profile)
</quiz>
