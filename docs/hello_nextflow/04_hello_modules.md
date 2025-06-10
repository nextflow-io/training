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

When we started developing our workflow, we put everything in one single code file.

Putting processes into individual modules makes it possible to reuse process definitions in multiple workflows without producing multiple copies of the code.
This makes the code more shareable, flexible and maintainable.

!!!note

    It is also possible to encapsulate a section of a workflow as a 'subworkflow' that can be imported into a larger pipeline, but that is outside the scope of this course.

---

## 0. Warmup: Run `hello-modules.nf`

We're going to use the workflow script `hello-modules.nf` as a starting point.
It is equivalent to the script produced by working through Part 3 of this training course.

Just to make sure everything is working, run the script once before making any changes:

```bash
nextflow run hello-modules.nf
```

```console title="Output"
 N E X T F L O W   ~  version 25.04.3

Launching `hello-modules.nf` [festering_nobel] DSL2 - revision: eeca64cdb1

executor >  local (7)
[25/648bdd] sayHello (2)       | 3 of 3 ✔
[60/bc6831] convertToUpper (1) | 3 of 3 ✔
[1a/bc5901] collectGreetings   | 1 of 1 ✔
There were 3 greetings in this batch
```

As previously, you will find the output files in the `results` directory (specified by the `publishDir` directive).

```console title="Directory contents"
results
├── Bonjour-output.txt
├── COLLECTED-output.txt
├── COLLECTED-test-batch-output.txt
├── COLLECTED-trio-output.txt
├── Hello-output.txt
├── Holà-output.txt
├── UPPER-Bonjour-output.txt
├── UPPER-Hello-output.txt
└── UPPER-Holà-output.txt
```

!!! note

    There may also be a file named `output.txt` left over if you worked through Part 2 in the same environment.

If that worked for you, you're ready to learn how to modularize your workflow code.

---

## 1. Create a directory to store modules

It is best practice to store your modules in a specific directory.
You can call that directory anything you want, but the convention is to call it `modules/`.

```bash
mkdir modules
```

!!! note

    Here we are showing how to use local modules, meaning modules stored locally in the same repository as the rest of the workflow code, in contrast to remote modules, which are stored in other (remote) repositories. For more information about remote modules, see the [documentation](https://www.nextflow.io/docs/latest/module.html).

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

    publishDir 'results', mode: 'copy'

    input:
        val greeting

    output:
        path "${greeting}-output.txt"

    script:
    """
    echo '$greeting' > '$greeting-output.txt'
    """
}
```

Once that is done, delete the process definition from the workflow file, but make sure to leave the shebang in place.

### 2.3. Add an import declaration before the workflow block

The syntax for importing a local module is fairly straightforward:

```groovy title="Syntax: Import declaration"
include { <MODULE_NAME> } from '<path_to_module>'
```

Let's insert that above the workflow block and fill it out appropriately.

=== "After"

    ```groovy title="hello-modules.nf" linenums="50" hl_lines="1 2"
    // Include modules
    include { sayHello } from './modules/sayHello.nf'

    workflow {
    ```

=== "Before"

    ```groovy title="hello-modules.nf" linenums="50"
    workflow {
    ```

### 2.4. Run the workflow to verify that it does the same thing as before

We're running the workflow with essentially the same code and inputs as before, so let's run with the `-resume` flag and see what happens.

```bash
nextflow run hello-modules.nf -resume
```

This runs quickly very quickly because everything is cached.

```console title="Output" linenums="1"
 N E X T F L O W   ~  version 25.04.3

Launching `hello-modules.nf` [romantic_poisson] DSL2 - revision: 96edfa9ad3

[f6/cc0107] sayHello (1)       | 3 of 3, cached: 3 ✔
[3c/4058ba] convertToUpper (2) | 3 of 3, cached: 3 ✔
[1a/bc5901] collectGreetings   | 1 of 1, cached: 1 ✔
There were 3 greetings in this batch
```

Nextflow recognized that it's still all the same work to be done, even if the code is split up into multiple files.

### Takeaway

You know how to extract a process into a local module and you know doing this doesn't break the resumability of the workflow.

### What's next?

Practice making more modules.
Once you've done one, you can do a million modules...
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

    publishDir 'results', mode: 'copy'

    input:
        path input_file

    output:
        path "UPPER-${input_file}"

    script:
    """
    cat '$input_file' | tr '[a-z]' '[A-Z]' > 'UPPER-${input_file}'
    """
}
```

Once that is done, delete the process definition from the workflow file, but make sure to leave the shebang in place.

### 3.3. Add an import declaration before the workflow block

Insert the import declaration above the workflow block and fill it out appropriately.

=== "After"

    ```groovy title="hello-modules.nf" linenums="31" hl_lines="3"
    // Include modules
    include { sayHello } from './modules/sayHello.nf'
    include { convertToUpper } from './modules/convertToUpper.nf'

    workflow {
    ```

=== "Before"

    ```groovy title="hello-modules.nf" linenums="31"
    // Include modules
    include { sayHello } from './modules/sayHello.nf'

    workflow {
    ```

### 3.4. Run the workflow to verify that it does the same thing as before

Run this with the `-resume` flag.

```bash
nextflow run hello-modules.nf -resume
```

This should still produce the same output as previously.

```console title="Output" linenums="1"
 N E X T F L O W   ~  version 25.04.3

Launching `hello-modules.nf` [nauseous_heisenberg] DSL2 - revision: a04a9f2da0

[c9/763d42] sayHello (3)       | 3 of 3, cached: 3 ✔
[60/bc6831] convertToUpper (3) | 3 of 3, cached: 3 ✔
[1a/bc5901] collectGreetings   | 1 of 1, cached: 1 ✔
There were 3 greetings in this batch
```

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

    publishDir 'results', mode: 'copy'

    input:
        path input_files
        val batch_name

    output:
        path "COLLECTED-${batch_name}-output.txt" , emit: outfile
        val count_greetings , emit: count

    script:
        count_greetings = input_files.size()
    """
    cat ${input_files} > 'COLLECTED-${batch_name}-output.txt'
    """
}
```

Once that is done, delete the process definition from the workflow file, but make sure to leave the shebang in place.

### 4.3. Add an import declaration before the workflow block

Insert the import declaration above the workflow block and fill it out appropriately.

=== "After"

    ```groovy title="hello-modules.nf" linenums="9" hl_lines="4"
    // Include modules
    include { sayHello } from './modules/sayHello.nf'
    include { convertToUpper } from './modules/convertToUpper.nf'
    include { collectGreetings } from './modules/collectGreetings.nf'

    workflow {
    ```

=== "Before"

    ```groovy title="hello-modules.nf" linenums="9"
    // Include modules
    include { sayHello } from './modules/sayHello.nf'
    include { convertToUpper } from './modules/convertToUpper.nf'

    workflow {
    ```

### 4.4. Run the workflow to verify that it does the same thing as before

Run this with the `-resume` flag.

```bash
nextflow run hello-modules.nf -resume
```

This should still produce the same output as previously.

```console title="Output" linenums="1"
 N E X T F L O W   ~  version 25.04.3

Launching `hello-modules.nf` [friendly_coulomb] DSL2 - revision: 7aa2b9bc0f

[f6/cc0107] sayHello (1)       | 3 of 3, cached: 3 ✔
[3c/4058ba] convertToUpper (2) | 3 of 3, cached: 3 ✔
[1a/bc5901] collectGreetings   | 1 of 1, cached: 1 ✔
There were 3 greetings in this batch
```

### Takeaway

You know how to modularize multiple processes in a workflow.

Congratulations, you've done all this work and absolutely nothing has changed to how the pipeline works!

Jokes aside, now your code is more modular, and if you decide to write another pipeline that calls on one of those processes, you just need to type one short import statement to use the relevant module.
This is better than just copy-pasting the code, because if later you decide to improve the module, all your pipelines will inherit the improvements.

### What's next?

Take a short break if you feel like it.
When you're ready, move on to Part 5 to learn how to use containers to manage software dependencies more conveniently and reproducibly.
