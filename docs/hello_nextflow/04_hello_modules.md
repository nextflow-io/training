# Part 4: Hello Modules

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/Xxp_menS0E8?si=0AWnXB7xqHAzJdJV&amp;list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } See [the whole playlist](https://www.youtube.com/playlist?list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik) on the Nextflow YouTube channel.

:green_book: The video transcript is available [here](./transcripts/04_hello_modules.md).
///

This section covers how to organize your workflow code to make development and maintenance of your pipeline more efficient and sustainable.
Specifically, we are going to demonstrate how to use **modules**.

In Nextflow, a **module** is a single process definition that is encapsulated by itself in a standalone code file.
To use a module in a workflow, you just add a single-line import statement to your workflow code file; then you can integrate the process into the workflow the same way you normally would.
That makes it possible to reuse process definitions in multiple workflows without producing multiple copies of the code.

When we started developing our workflow, we wrote everything in one single code file.
Now we're going to move the processes out into individual modules.

<figure class="excalidraw">
--8<-- "docs/hello_nextflow/img/modules.svg"
</figure>

This will make our code more shareable, flexible and maintainable.

??? info "How to begin from this section"

    This section of the course assumes you have completed Parts 1-3 of the [Hello Nextflow](./index.md) course, but if you are comfortable with the basics covered in those sections, you can start from here without doing anything special.

---

## 0. Warmup: Run `hello-modules.nf`

We're going to use the workflow script `hello-modules.nf` as a starting point.
It is equivalent to the script produced by working through Part 3 of this training course, except we've changed the output destinations:

```groovy title="hello-modules.nf" linenums="37" hl_lines="3 7 11 15"
output {
    first_output {
        path 'hello_modules'
        mode 'copy'
    }
    uppercased {
        path 'hello_modules'
        mode 'copy'
    }
    collected {
        path 'hello_modules'
        mode 'copy'
    }
    batch_report {
        path 'hello_modules'
        mode 'copy'
    }
}
```

Just to make sure everything is working, run the script once before making any changes:

```bash
nextflow run hello-modules.nf
```

??? success "Command output"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-modules.nf` [hopeful_avogadro] DSL2 - revision: b09af1237d

    executor >  local (7)
    [0f/8795c9] sayHello (3)       [100%] 3 of 3 ✔
    [6a/eb2510] convertToUpper (3) [100%] 3 of 3 ✔
    [af/479117] collectGreetings   [100%] 1 of 1 ✔
    ```

As previously, you will find the output files in the directory specified in the `output` block (here, `results/hello_modules/`).

??? abstract "Directory contents"

    ```console
    results/hello_modules/
    ├── Bonjour-output.txt
    ├── COLLECTED-batch-output.txt
    ├── Hello-output.txt
    ├── Holà-output.txt
    ├── batch-report.txt
    ├── UPPER-Bonjour-output.txt
    ├── UPPER-Hello-output.txt
    └── UPPER-Holà-output.txt
    ```

If that worked for you, you're ready to learn how to modularize your workflow code.

---

## 1. Create a directory to store modules

It is best practice to store your modules in a specific directory.
You can call that directory anything you want, but the convention is to call it `modules/`.

```bash
mkdir modules
```

!!! tip

    Here we are showing you how to use **local modules**, meaning modules stored locally in the same repository as the rest of the workflow code, in contrast to remote modules, which are stored in other (remote) repositories.
    For more information about **remote modules**, see the [documentation](https://www.nextflow.io/docs/latest/module.html).

---

## 2. Create a module for `sayHello()`

In its simplest form, turning an existing process into a module is little more than a copy-paste operation.
We're going to create a file stub for the module, copy the relevant code over then delete it from the main workflow file.

Then all we'll need to do is add an import statement so that Nextflow will know to pull in the relevant code at runtime.

### 2.1. Create a file stub for the new module

Let's create an empty file for the module called `sayHello.nf`.

```bash
touch modules/sayHello.nf
```

This gives us a place to put the process code.

### 2.2. Move the `sayHello` process code to the module file

Copy the whole process definition over from the workflow file to the module file, making sure to copy over the `#!/usr/bin/env nextflow` shebang too.

```groovy title="modules/sayHello.nf" linenums="1"
#!/usr/bin/env nextflow

/*
 * Use echo to print 'Hello World!' to a file
 */
process sayHello {

    input:
    val greeting

    output:
    path "${greeting}-output.txt"

    script:
    """
    echo '${greeting}' > '${greeting}-output.txt'
    """
}
```

Once that is done, delete the process definition from the workflow file, but make sure to leave the shebang in place.

### 2.3. Add an import declaration before the workflow block

The syntax for importing a local module is fairly straightforward:

```groovy title="Syntax: Import declaration"
include { <MODULE_NAME> } from '<path_to_module>'
```

Let's insert that above the `params` block and fill it out appropriately.

=== "After"

    ```groovy title="hello-modules.nf" linenums="44" hl_lines="1 2"
    // Include modules
    include { sayHello } from './modules/sayHello.nf'

    /*
    * Pipeline parameters
    */
    params {
        greeting: Path = 'data/greetings.csv'
        batch: String = 'batch'
    }
    ```

=== "Before"

    ```groovy title="hello-modules.nf" linenums="44"
    /*
    * Pipeline parameters
    */
    params {
        greeting: Path = 'data/greetings.csv'
        batch: String = 'batch'
    }
    ```

You see we've filled in the module name, `sayHello`, and the path to the file containing the module code, `./modules/sayHello.nf`.

### 2.4. Run the workflow

We're running the workflow with essentially the same code and inputs as before, so let's run with the `-resume` flag and see what happens.

```bash
nextflow run hello-modules.nf -resume
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-modules.nf` [romantic_poisson] DSL2 - revision: 96edfa9ad3

    [f6/cc0107] sayHello (1)       | 3 of 3, cached: 3 ✔
    [3c/4058ba] convertToUpper (2) | 3 of 3, cached: 3 ✔
    [1a/bc5901] collectGreetings   | 1 of 1, cached: 1 ✔
    ```

This should run quickly very quickly because everything is cached.
Feel free to check the published outputs.

Nextflow recognized that it's still all the same work to be done, even if the code is split up into multiple files.

### Takeaway

You know how to extract a process into a local module and you know doing this doesn't break the resumability of the workflow.

### What's next?

Practice making more modules.
Once you've done one, you can do a million more...
But let's just do two more for now.

---

## 3. Modularize the `convertToUpper()` process

### 3.1. Create a file stub for the new module

Create an empty file for the module called `convertToUpper.nf`.

```bash
touch modules/convertToUpper.nf
```

### 3.2. Move the `convertToUpper` process code to the module file

Copy the whole process definition over from the workflow file to the module file, making sure to copy over the `#!/usr/bin/env nextflow` shebang too.

```groovy title="modules/convertToUpper.nf" linenums="1"
#!/usr/bin/env nextflow

/*
 * Use a text replacement tool to convert the greeting to uppercase
 */
process convertToUpper {

    input:
    path input_file

    output:
    path "UPPER-${input_file}"

    script:
    """
    cat '${input_file}' | tr '[a-z]' '[A-Z]' > 'UPPER-${input_file}'
    """
}
```

Once that is done, delete the process definition from the workflow file, but make sure to leave the shebang in place.

### 3.3. Add an import declaration before the `params` block

Insert the import declaration above the `params` block and fill it out appropriately.

=== "After"

    ```groovy title="hello-modules.nf" linenums="23" hl_lines="3"
    // Include modules
    include { sayHello } from './modules/sayHello.nf'
    include { convertToUpper } from './modules/convertToUpper.nf'

    /*
    * Pipeline parameters
    */
    params {
        greeting: Path = 'data/greetings.csv'
        batch: String = 'batch'
    }
    ```

=== "Before"

    ```groovy title="hello-modules.nf" linenums="23"
    // Include modules
    include { sayHello } from './modules/sayHello.nf'

    /*
    * Pipeline parameters
    */
    params {
        greeting: Path = 'data/greetings.csv'
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
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-modules.nf` [nauseous_heisenberg] DSL2 - revision: a04a9f2da0

    [c9/763d42] sayHello (3)       | 3 of 3, cached: 3 ✔
    [60/bc6831] convertToUpper (3) | 3 of 3, cached: 3 ✔
    [1a/bc5901] collectGreetings   | 1 of 1, cached: 1 ✔
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

Copy the whole process definition over from the workflow file to the module file, making sure to copy over the `#!/usr/bin/env nextflow` shebang too.

```groovy title="modules/collectGreetings.nf" linenums="1"
#!/usr/bin/env nextflow

/*
 * Collect uppercase greetings into a single output file
 */
process collectGreetings {

    input:
    path input_files
    val batch_name

    output:
    path "COLLECTED-${batch_name}-output.txt", emit: outfile
    path "${batch_name}-report.txt", emit: report

    script:
    count_greetings = input_files.size()
    """
    cat ${input_files} > 'COLLECTED-${batch_name}-output.txt'
    echo 'There were ${count_greetings} greetings in this batch.' > '${batch_name}-report.txt'
    """
}
```

Once that is done, delete the process definition from the workflow file, but make sure to leave the shebang in place.

### 4.3. Add an import declaration before the `params` block

Insert the import declaration above the `params` block and fill it out appropriately.

=== "After"

    ```groovy title="hello-modules.nf" linenums="3" hl_lines="4"
    // Include modules
    include { sayHello } from './modules/sayHello.nf'
    include { convertToUpper } from './modules/convertToUpper.nf'
    include { collectGreetings } from './modules/collectGreetings.nf'

    /*
    * Pipeline parameters
    */
    params {
        greeting: Path = 'data/greetings.csv'
        batch: String = 'batch'
    }
    ```

=== "Before"

    ```groovy title="hello-modules.nf" linenums="3"
    // Include modules
    include { sayHello } from './modules/sayHello.nf'
    include { convertToUpper } from './modules/convertToUpper.nf'

    /*
    * Pipeline parameters
    */
    params {
        greeting: Path = 'data/greetings.csv'
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
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-modules.nf` [friendly_coulomb] DSL2 - revision: 7aa2b9bc0f

    [f6/cc0107] sayHello (1)       | 3 of 3, cached: 3 ✔
    [3c/4058ba] convertToUpper (2) | 3 of 3, cached: 3 ✔
    [1a/bc5901] collectGreetings   | 1 of 1, cached: 1 ✔
    ```

This should still produce the same output as previously.

### Takeaway

You know how to modularize multiple processes in a workflow.

Congratulations, you've done all this work and absolutely nothing has changed to how the pipeline works!

Jokes aside, now your code is more modular, and if you decide to write another pipeline that calls on one of those processes, you just need to type one short import statement to use the relevant module.
This is better than copy-pasting the code, because if later you decide to improve the module, all your pipelines will inherit the improvements.

### What's next?

Take a short break if you feel like it.

When you're ready, move on to [**Part 5: Hello Containers**](./05_hello_containers.md) to learn how to use containers to manage software dependencies more conveniently and reproducibly.

---

## Quiz

<quiz>
What is a module in Nextflow?
- [ ] A configuration file
- [x] A standalone file containing a single process definition
- [ ] A workflow definition
- [ ] A channel operator
</quiz>

<quiz>
What is the recommended naming convention for module files?
- [ ] module_processName.nf
- [ ] processName_module.nf
- [x] processName.nf
- [ ] mod_processName.nf
</quiz>

<quiz>
Where should module files be stored?
- [ ] In the same directory as the workflow
- [ ] In a bin/ directory
- [x] In a modules/ directory
- [ ] In a lib/ directory
</quiz>

<quiz>
What is the correct syntax to import a module?
- [ ] import { processName } from './modules/processName.nf'
- [ ] require { processName } from './modules/processName.nf'
- [x] include { processName } from './modules/processName.nf'
- [ ] load { processName } from './modules/processName.nf'
</quiz>

<quiz>
What happens to the -resume functionality when using modules?
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
