# Part 4: Hello Modules

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/43Ot-f0iOME?si=y8lAedhEHWaTV4zd&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } See [the whole playlist](https://youtube.com/playlist?list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&si=eF7cLR62goy-lc6n) on the Nextflow YouTube channel.

:green_book: The video transcript is available [here](./transcripts/04_hello_modules.md).
///

This section covers how to organize your workflow code to make development and maintenance of your pipeline more efficient and sustainable.
Specifically, we are going to demonstrate how to use [**modules**](https://nextflow.io/docs/latest/module.html).

In Nextflow, a **module** is a standalone code file, often encapsulating a single process definition.
To use a module in a workflow, you just add a single-line `include` statement to your workflow code file; then you can integrate the process into the workflow the same way you normally would.
That makes it possible to reuse process definitions in multiple workflows without producing multiple copies of the code.

When we started developing our workflow, we wrote everything in one single code file.
Now we're going to put each piece of code in its proper place: the processes go into individual modules, and the parameter values go into the configuration file.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/modules.svg"
</figure>

This will make our code more shareable, flexible and maintainable.

??? info "How to begin from this section"

    This section of the course assumes you have completed Parts 1-3 of the [Hello Nextflow](./index.md) course, but if you are comfortable with the basics covered in those sections, you can start from here without doing anything special.

---

## 0. Warmup: Run `hello-modules.nf`

We're going to use the workflow script `hello-modules.nf` as a starting point.
It is equivalent to the script produced by working through Part 3 of this training course.

Just to make sure everything is working, run the script once before making any changes:

```bash
nextflow run hello-modules.nf
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `hello-modules.nf` [cranky_waddington] revision: b4393cf728

    executor >  local (7)
    [d5/968820] sayHello (3)          | 3 of 3 ✔
    [27/1b8c0b] translateGreeting (2) | 3 of 3 ✔
    [51/cd354d] collectGreetings      | 1 of 1 ✔

    Outputs:

      /workspaces/training/hello-nextflow/results

      greetings:
        - batch/intermediates/Mateo.txt
        - batch/intermediates/Celine.txt
        - batch/intermediates/Andrea.txt

      translated:
        - batch/intermediates/es-Mateo.txt
        - batch/intermediates/es-Andrea.txt
        - batch/intermediates/es-Celine.txt

      collected: batch/intermediates/batch-greetings.txt

      batch_report: batch/batch-report.txt
    ```

As previously, you will find the output files in the locations specified in the `output` block, under `results/batch/`, since we didn't set a value for `--batch`.

??? abstract "Directory contents"

    ```console
    results/batch
    ├── batch-report.txt
    └── intermediates
        ├── Andrea.txt
        ├── batch-greetings.txt
        ├── Celine.txt
        ├── es-Andrea.txt
        ├── es-Celine.txt
        ├── es-Mateo.txt
        └── Mateo.txt
    ```

If that worked for you, you're ready to learn how to modularize your workflow code.

---

## 1. Create a directory to store modules

It is best practice to store your modules in a specific directory.
You can call that directory anything you want, but the convention is to call it `modules/`.

```bash
mkdir modules
```

---

## 2. Create a module for `sayHello()`

In its simplest form, turning an existing process into a module is little more than a copy-paste operation.
We're going to create a file stub for the module, copy the relevant code over then delete it from the main workflow file.

Then all we'll need to do is add an `include` statement so that Nextflow will know to pull in the relevant code at runtime.

### 2.1. Create a file stub for the new module

Let's create an empty file for the module called `sayHello.nf`.

```bash
touch modules/sayHello.nf
```

This gives us a place to put the process code.

### 2.2. Move the `sayHello` process code to the module file

Copy the whole process definition over from the workflow file to the module file.

```groovy title="modules/sayHello.nf" linenums="1"
/*
 * Use echo to print a greeting to a file
 */
process sayHello {

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

### 2.3. Add an include declaration before the `params` block

The syntax for including a process from a module is fairly straightforward:

```groovy title="Syntax: include declaration"
include { <PROCESS_NAME> } from '<path_to_module>'
```

Let's insert that above the `params` block and fill it out appropriately.

=== "After"

    ```groovy title="hello-modules.nf" linenums="62" hl_lines="1 2"
    // Include modules
    include { sayHello } from './modules/sayHello.nf'

    /*
     * Pipeline parameters
     */
    params {
        input: Path = 'data/names.csv'
        language: String = 'es'
        batch: String = 'batch'
    }
    ```

=== "Before"

    ```groovy title="hello-modules.nf" linenums="62"
    /*
     * Pipeline parameters
     */
    params {
        input: Path = 'data/names.csv'
        language: String = 'es'
        batch: String = 'batch'
    }
    ```

You see we've filled in the process name, `sayHello`, and the path to the file containing the module code, `./modules/sayHello.nf`.

### 2.4. Run the workflow

We're running the workflow with essentially the same code and inputs as before, so let's run with the `-resume` flag and see what happens.

```bash
nextflow run hello-modules.nf -resume
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `hello-modules.nf` [irreverent_nightingale] revision: c301d4f9fa

    [d5/968820] sayHello (3)          | 3 of 3, cached: 3 ✔
    [27/1b8c0b] translateGreeting (2) | 3 of 3, cached: 3 ✔
    [51/cd354d] collectGreetings      | 1 of 1, cached: 1 ✔

    Outputs:

      /workspaces/training/hello-nextflow/results

      greetings:
        - batch/intermediates/Andrea.txt
        - batch/intermediates/Mateo.txt
        - batch/intermediates/Celine.txt

      translated:
        - batch/intermediates/es-Celine.txt
        - batch/intermediates/es-Andrea.txt
        - batch/intermediates/es-Mateo.txt

      collected: batch/intermediates/batch-greetings.txt

      batch_report: batch/batch-report.txt
    ```

This should run very quickly because everything is cached.
Feel free to check the published outputs.

Nextflow recognized that it's still all the same work to be done, even if the code is split up into multiple files.

### Takeaway

You know how to extract a process into a local module and you know doing this doesn't break the resumability of the workflow.

### What's next?

Practice making more modules.
Once you've done one, you can do a million more...
But let's just do two more for now.

---

## 3. Modularize the `translateGreeting()` process

### 3.1. Create a file stub for the new module

Create an empty file for the module called `translateGreeting.nf`.

```bash
touch modules/translateGreeting.nf
```

### 3.2. Move the `translateGreeting` process code to the module file

Copy the whole process definition over from the workflow file to the module file.

```groovy title="modules/translateGreeting.nf" linenums="1"
/*
 * Use Python to translate the greeting into another language
 */
process translateGreeting {

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

A module can hold a script block in any language: the Python code moves over exactly as it was.

Once that is done, delete the process definition from the workflow file.

### 3.3. Add an include declaration before the `params` block

Insert the include declaration above the `params` block and fill it out appropriately.

=== "After"

    ```groovy title="hello-modules.nf" linenums="24" hl_lines="3"
    // Include modules
    include { sayHello } from './modules/sayHello.nf'
    include { translateGreeting } from './modules/translateGreeting.nf'

    /*
     * Pipeline parameters
     */
    params {
        input: Path = 'data/names.csv'
        language: String = 'es'
        batch: String = 'batch'
    }
    ```

=== "Before"

    ```groovy title="hello-modules.nf" linenums="24"
    // Include modules
    include { sayHello } from './modules/sayHello.nf'

    /*
     * Pipeline parameters
     */
    params {
        input: Path = 'data/names.csv'
        language: String = 'es'
        batch: String = 'batch'
    }
    ```

This should start to look very familiar.

### 3.4. Run the workflow again

Run this with the `-resume` flag.

```bash
nextflow run hello-modules.nf -resume
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `hello-modules.nf` [cranky_shannon] revision: e8669b1308

    [d5/968820] sayHello (3)          | 3 of 3, cached: 3 ✔
    [27/1b8c0b] translateGreeting (2) | 3 of 3, cached: 3 ✔
    [51/cd354d] collectGreetings      | 1 of 1, cached: 1 ✔

    Outputs:

      /workspaces/training/hello-nextflow/results

      greetings:
        - batch/intermediates/Celine.txt
        - batch/intermediates/Andrea.txt
        - batch/intermediates/Mateo.txt

      translated:
        - batch/intermediates/es-Mateo.txt
        - batch/intermediates/es-Andrea.txt
        - batch/intermediates/es-Celine.txt

      collected: batch/intermediates/batch-greetings.txt

      batch_report: batch/batch-report.txt
    ```

This should still produce the same output as previously.

Two done, one more to go!

---

## 4. Modularize the `collectGreetings()` process

### 4.1. Create a file stub for the new module

Create an empty file for the module called `collectGreetings.nf`.

```bash
touch modules/collectGreetings.nf
```

### 4.2. Move the `collectGreetings` process code to the module file

Copy the whole process definition over from the workflow file to the module file.

```groovy title="modules/collectGreetings.nf" linenums="1"
/*
 * Collect translated greetings into a single output file
 */
process collectGreetings {

    input:
    path input_files
    val batch_name

    output:
    path "${batch_name}-greetings.txt", emit: outfile
    path "${batch_name}-report.txt", emit: report

    script:
    count_greetings = input_files.size()
    """
    cat ${input_files} > '${batch_name}-greetings.txt'
    echo 'Number of greetings processed: ${count_greetings}' > '${batch_name}-report.txt'
    """
}
```

Once that is done, delete the process definition from the workflow file.

### 4.3. Add an include declaration before the `params` block

Insert the include declaration above the `params` block and fill it out appropriately.

=== "After"

    ```groovy title="hello-modules.nf" linenums="3" hl_lines="4"
    // Include modules
    include { sayHello } from './modules/sayHello.nf'
    include { translateGreeting } from './modules/translateGreeting.nf'
    include { collectGreetings } from './modules/collectGreetings.nf'

    /*
     * Pipeline parameters
     */
    params {
        input: Path = 'data/names.csv'
        language: String = 'es'
        batch: String = 'batch'
    }
    ```

=== "Before"

    ```groovy title="hello-modules.nf" linenums="3"
    // Include modules
    include { sayHello } from './modules/sayHello.nf'
    include { translateGreeting } from './modules/translateGreeting.nf'

    /*
     * Pipeline parameters
     */
    params {
        input: Path = 'data/names.csv'
        language: String = 'es'
        batch: String = 'batch'
    }
    ```

Last one!

### 4.4. Run the workflow

Run this with the `-resume` flag.

```bash
nextflow run hello-modules.nf -resume
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `hello-modules.nf` [irreverent_shaw] revision: bd0f694fe5

    [d5/968820] sayHello (3)          | 3 of 3, cached: 3 ✔
    [27/1b8c0b] translateGreeting (2) | 3 of 3, cached: 3 ✔
    [51/cd354d] collectGreetings      | 1 of 1, cached: 1 ✔

    Outputs:

      /workspaces/training/hello-nextflow/results

      greetings:
        - batch/intermediates/Mateo.txt
        - batch/intermediates/Celine.txt
        - batch/intermediates/Andrea.txt

      translated:
        - batch/intermediates/es-Andrea.txt
        - batch/intermediates/es-Celine.txt
        - batch/intermediates/es-Mateo.txt

      collected: batch/intermediates/batch-greetings.txt

      batch_report: batch/batch-report.txt
    ```

This should still produce the same output as previously.

### Takeaway

You know how to modularize multiple processes in a workflow.

Congratulations, you've done all this work and absolutely nothing has changed to how the pipeline works!

Jokes aside, now your code is more modular, and if you decide to write another pipeline that calls on one of those processes, you just need to type one short `include` statement to use the relevant module.
This is better than copy-pasting the code, because if later you decide to improve the module, all your pipelines will inherit the improvements.

### What's next?

Move the parameter values out of the workflow script.

---

## 5. Move parameter values to the configuration

The workflow script still sets default values for its parameters in the `params` block:

```groovy title="hello-modules.nf" linenums="8"
/*
 * Pipeline parameters
 */
params {
    input: Path = 'data/names.csv'
    language: String = 'es'
    batch: String = 'batch'
}
```

That was convenient while we were developing the workflow, but these values are not all the same kind of thing:

- `language = 'es'` and `batch = 'batch'` are true defaults: sensible values for anyone who runs the pipeline without specifying them.
- `input = 'data/names.csv'` is a test value: it points to the example data provided with the course, which nobody would want to process by default in a real analysis.

As you saw in Nextflow Run, true defaults belong in the `params` section of `nextflow.config`, and test values belong in a `test` profile.
The workflow script keeps the parameter declarations and their types, so it is still easy to see at a glance which parameters the pipeline accepts.

### 5.1. Remove the values from the workflow script

In the `params` block of the workflow script, make the following code change:

=== "After"

    ```groovy title="hello-modules.nf" linenums="8" hl_lines="5-7"
    /*
     * Pipeline parameters
     */
    params {
        input: Path
        language: String
        batch: String
    }
    ```

=== "Before"

    ```groovy title="hello-modules.nf" linenums="8"
    /*
     * Pipeline parameters
     */
    params {
        input: Path = 'data/names.csv'
        language: String = 'es'
        batch: String = 'batch'
    }
    ```

### 5.2. Add the default value and the test profile to `nextflow.config`

Open `nextflow.config` and make the following change:

=== "After"

    ```groovy title="nextflow.config" linenums="1" hl_lines="3-18"
    docker.enabled = false

    /*
     * Pipeline parameters
     */
    params {
        language = 'es'
        batch = 'batch'
    }

    /*
     * Profiles
     */
    profiles {
        test {
            params.input = 'data/names.csv'
        }
    }
    ```

=== "Before"

    ```groovy title="nextflow.config" linenums="1"
    docker.enabled = false
    ```

The `params` section of the configuration file uses `=` to assign values, without types: the types are declared in the workflow script.

### 5.3. Run the workflow without the test profile

Try running the workflow the same way as before:

```bash
nextflow run hello-modules.nf -resume
```

??? failure "Command output"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `hello-modules.nf` [modest_austin] revision: 6024dbeca6

    Parameter `input` is required but was not specified on the command line, params file, or config
    ```

Since `input` no longer has a default value, Nextflow stops before running anything and tells you that the parameter is missing.
`language` and `batch` did not cause an error, because their default values now come from `nextflow.config`.

### 5.4. Run the workflow with the test profile

Add `-profile test` to provide the test value for `input`:

```bash
nextflow run hello-modules.nf -profile test -resume
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `hello-modules.nf` [jolly_einstein] revision: 6024dbeca6

    [d5/968820] sayHello (3)          | 3 of 3, cached: 3 ✔
    [27/1b8c0b] translateGreeting (2) | 3 of 3, cached: 3 ✔
    [51/cd354d] collectGreetings      | 1 of 1, cached: 1 ✔

    Outputs:

      /workspaces/training/hello-nextflow/results

      greetings:
        - batch/intermediates/Andrea.txt
        - batch/intermediates/Celine.txt
        - batch/intermediates/Mateo.txt

      translated:
        - batch/intermediates/es-Andrea.txt
        - batch/intermediates/es-Celine.txt
        - batch/intermediates/es-Mateo.txt

      collected: batch/intermediates/batch-greetings.txt

      batch_report: batch/batch-report.txt
    ```

The parameter values are the same as before, so every task is retrieved from the cache.

From now on, run the workflow with `-profile test`, or provide `--input` on the command line.

### Takeaway

You know how to separate the parameter declarations in the workflow script from their values in the configuration, keeping true default values in the `params` section and test values in a `test` profile.

### What's next?

Take a short break if you feel like it.

When you're ready, move on to [**Part 5: Hello Containers**](./05_hello_containers.md) to learn how to use containers to manage software dependencies more conveniently and reproducibly.

---

## Bonus exercises

These optional exercises let you practice what you learned in this part.
Each one can be solved with what this part covered; open the hint if you get stuck.
The solutions are also in `solutions/4-hello-modules/bonus/`.

### Exercise 1: Add a profile for French greetings

Starting from your `hello-modules.nf` and `nextflow.config` at the end of this part, add a profile named `french` that sets the language to French.
Then run the workflow with both the `test` and `french` profiles, and check that the collected greetings are in French.

??? tip "Hint"

    Write it like the `test` profile, setting `params.language` instead of `params.input`.
    To activate two profiles, separate their names with a comma.

??? solution

    In `nextflow.config`, add the `french` profile next to the `test` profile:

    ```groovy title="nextflow.config" linenums="14" hl_lines="5-7"
    profiles {
        test {
            params.input = 'data/names.csv'
        }
        french {
            params.language = 'fr'
        }
    }
    ```

    Then run the workflow with both profiles:

    ```bash
    nextflow run hello-modules.nf -profile test,french -resume
    ```

    ```console title="Command output"
    N E X T F L O W   ~  version 26.04.4

    Launching `hello-modules.nf` [agitated_kalman] revision: 6024dbeca6

    executor >  local (4)
    [d5/968820] sayHello (3)          | 3 of 3, cached: 3 ✔
    [ed/d7f3e3] translateGreeting (3) | 3 of 3 ✔
    [80/1d28d0] collectGreetings      | 1 of 1 ✔

    Outputs:

      /workspaces/training/hello-nextflow/results

      greetings:
        - batch/intermediates/Andrea.txt
        - batch/intermediates/Celine.txt
        - batch/intermediates/Mateo.txt

      translated:
        - batch/intermediates/fr-Celine.txt
        - batch/intermediates/fr-Andrea.txt
        - batch/intermediates/fr-Mateo.txt

      collected: batch/intermediates/batch-greetings.txt

      batch_report: batch/batch-report.txt
    ```

    The `sayHello` tasks are cached, and the translation runs again in French:

    ```console title="results/batch/intermediates/batch-greetings.txt"
    Bonjour, Celine!
    Bonjour, Andrea!
    Bonjour, Mateo!
    ```

### Exercise 2: Give each module its own folder

Larger projects often store each module in its own folder, as `modules/<name>/main.nf`, so that files related to the module can sit next to it.
Starting from your `hello-modules.nf` at the end of this part, move each module file into a folder of its own, with this layout, and update the workflow so it runs as before:

```console
modules
├── collectGreetings
│   └── main.nf
├── sayHello
│   └── main.nf
└── translateGreeting
    └── main.nf
```

??? tip "Hint"

    Only the path in each `include` declaration changes.

??? solution

    After moving the files, update the `include` declarations in `hello-modules.nf`:

    ```groovy title="hello-modules.nf" linenums="3" hl_lines="2-4"
    // Include modules
    include { sayHello } from './modules/sayHello/main.nf'
    include { translateGreeting } from './modules/translateGreeting/main.nf'
    include { collectGreetings } from './modules/collectGreetings/main.nf'
    ```

    Then run the workflow:

    ```bash
    nextflow run hello-modules.nf -profile test -resume
    ```

    ```console title="Command output"
    N E X T F L O W   ~  version 26.04.4

    Launching `hello-modules.nf` [romantic_baekeland] revision: 1a860c4bd9

    [d5/968820] sayHello (3)          | 3 of 3, cached: 3 ✔
    [27/1b8c0b] translateGreeting (2) | 3 of 3, cached: 3 ✔
    [51/cd354d] collectGreetings      | 1 of 1, cached: 1 ✔

    Outputs:

      /workspaces/training/hello-nextflow/results

      greetings:
        - batch/intermediates/Celine.txt
        - batch/intermediates/Mateo.txt
        - batch/intermediates/Andrea.txt

      translated:
        - batch/intermediates/es-Celine.txt
        - batch/intermediates/es-Mateo.txt
        - batch/intermediates/es-Andrea.txt

      collected: batch/intermediates/batch-greetings.txt

      batch_report: batch/batch-report.txt
    ```

    Every task is retrieved from the cache: moving a module to another file doesn't change the work it does.

---

## Quiz

<quiz>
What is a module in Nextflow?
- [ ] A configuration file
- [x] A standalone file that can contain process definitions
- [ ] A workflow definition
- [ ] A channel operator

Learn more: [2. Create a module for `sayHello()`](#2-create-a-module-for-sayhello)
</quiz>

<quiz>
What convention is typically used for storing module files?
- [ ] In the same directory as the workflow
- [ ] In a `bin/` directory
- [x] In a `modules/` directory
- [ ] In a `lib/` directory

Learn more: [1. Create a directory to store modules](#1-create-a-directory-to-store-modules)
</quiz>

<quiz>
What is the correct syntax to use a module?

- [ ] `#!groovy import { sayHello } from './modules/sayHello.nf'`
- [ ] `#!groovy require { sayHello } from './modules/sayHello.nf'`
- [x] `#!groovy include { sayHello } from './modules/sayHello.nf'`
- [ ] `#!groovy load { sayHello } from './modules/sayHello.nf'`

Learn more: [2.3. Add an include declaration](#23-add-an-include-declaration-before-the-params-block)
</quiz>

<quiz>
What do you need to change in a process whose script is written in Python when you move it into a module?
- [ ] Rewrite the script in Bash
- [ ] Add the Python interpreter to the `include` statement
- [x] Nothing: the process moves over exactly as it was
- [ ] Declare the module as a Python module

Learn more: [3.2. Move the `translateGreeting` process code to the module file](#32-move-the-translategreeting-process-code-to-the-module-file)
</quiz>

<quiz>
What happens to the `-resume` functionality when using modules?
- [ ] It no longer works
- [ ] It requires additional configuration
- [x] It works the same as before
- [ ] It only works for local modules
</quiz>

<quiz>
What are the benefits of using modules? (Select all that apply)
- [x] Code reusability across workflows
- [x] Easier maintenance
- [x] Better organization of workflow code
- [ ] Faster execution speed
</quiz>

<quiz>
Where should the path to the example input data used for testing go?
- [ ] In the `params` block of the workflow script, as a default value
- [ ] In the `params` section of `nextflow.config`
- [x] In a `test` profile in `nextflow.config`
- [ ] In the module that uses it

Learn more: [5. Move parameter values to the configuration](#5-move-parameter-values-to-the-configuration)
</quiz>
