# Part 4: Hello Modules

This section covers how to organize your workflow code to make development and maintenance of your pipeline more efficient and sustainable.
Specifically, we are going to demonstrate how to use **modules**.

In Nextflow, a **module** is a single process definition that is encapsulated by itself in a standalone code file.
To use a module in a workflow, you just add a single-line import statement to your workflow code file; then you can integrate the process into the workflow the same way you normally would.

When we started developing our workflow, we put everything in one single code file.

Putting processes into individual modules makes it possible to reuse process definitions in multiple workflows without producing multiple copies of the code.
This makes the code more shareable, flexible and maintainable.

!!!note

    It is also possible to encapsulate a section of a workflow as a 'subworkflow' that can be imported into a larger pipeline, but that is outside the scope of this training.

---

## 0. Warmup: Run `hello-modules.nf`

We're going to use the workflow script `hello-modules.nf` as a starting point.
It is equivalent to the script produced by working through Part 3 of this training course.

Just to make sure everything is working, run the script once before making any changes:

```bash
nextflow run hello-modules.nf
```

This should produce the following output:

```console title="Output"
 N E X T F L O W   ~  version 24.10.0

Launching `hello-modules.nf` [tender_becquerel] DSL2 - revision: f7cat8e223

executor >  local (7)
[bd/4bb541] sayHello (1)         [100%] 3 of 3 ✔
[85/b627e8] convertToUpper (3)   [100%] 3 of 3 ✔
[7d/f7961c] collectGreetings     [100%] 1 of 1 ✔
```

---

## 1. Create a directory to store modules

It is best practice to store your modules in a specific directory. You can call that directory anything you want, but the convention is to call it `modules/`.

```bash
mkdir modules
```

!!! note

    Here we are showing how to use local modules, meaning modules stored locally in the same repository as the rest of the workflow code, in contrast to remote modules, which are stored in other (remote) repositories. For more information about remote modules, see the [documentation](https://www.nextflow.io/docs/latest/module.html).

---

## 2. Modularize the `sayHello()` process

Turning an existing process into a module is little more than a copy-paste operation. We're going to create a file stub for the module, copy the relevant code over then delete it from the main workflow file.

Then all we need to do is add an import statement so that Nextflow will know to pull in the relevant code at runtime.

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

[TODO]
```

Once that is done, delete the process definition from the workflow file, but make sure to leave the shebang in place.

### 2.3. Add an import declaration before the workflow block

The syntax for importing a local module is fairly straightforward:

```groovy title="Import declaration syntax"
include { <MODULE_NAME> } from '<path_to_module>'
```

Let's insert that above the workflow block and fill it out appropriately.

_Before:_

```groovy title="hello-modules.nf" linenums="73"
workflow {
```

_After:_

```groovy title="hello-modules.nf" linenums="73"
// Include modules
include { sayHello } from './modules/sayHello.nf'

workflow {
```

### 2.4. Run the workflow to verify that it does the same thing as before

We're running the workflow with essentially the same code and inputs as before, so let's run with the `-resume` flag and see what happens.

```bash
nextflow run hello-modules.nf -resume
```

Sure enough, Nextflow recognizes that it's still all the same work to be done, even if the code is split up into multiple files.

```console title="Output"
[TODO]
```

So modularizing the code in the course of development does not break resumability!

### Takeaway

You know how to extract a process into a local module.

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

This gives us a place to put the process code.

### 3.2. Move the `convertToUpper` process code to the module file

Copy the whole process definition over from the workflow file to the module file, making sure to copy over the `#!/usr/bin/env nextflow` shebang too.

```groovy title="modules/convertToUpper.nf" linenums="1"
#!/usr/bin/env nextflow

[TODO]
```

Once that is done, delete the process definition from the workflow file, but make sure to leave the shebang in place.

### 3.3. Add an import declaration before the workflow block

Insert the import declaration above the workflow block and fill it out appropriately.

_Before:_

```groovy title="hello-modules.nf" linenums="73"
// Include modules
include { sayHello } from './modules/sayHello.nf'

workflow {
```

_After:_

```groovy title="hello-modules.nf" linenums="73"
// Include modules
include { sayHello } from './modules/sayHello.nf'
include { convertToUpper } from './modules/convertToUpper.nf'

workflow {
```

### 3.4. Run the workflow to verify that it does the same thing as before

Run this with the `-resume` flag.

```bash
nextflow run hello-modules.nf -resume
```

This should still produce the same output as previously.

---

## 4. Modularize the `collectGreetings()` process

### 4.1. Create a file stub for the new module

Create an empty file for the module called `collectGreetings.nf`.

```bash
touch modules/collectGreetings.nf
```

This gives us a place to put the process code.

### 4.2. Move the `collectGreetings` process code to the module file

Copy the whole process definition over from the workflow file to the module file, making sure to copy over the `#!/usr/bin/env nextflow` shebang too.

```groovy title="modules/collectGreetings.nf" linenums="1"
#!/usr/bin/env nextflow

[TODO]
```

Once that is done, delete the process definition from the workflow file, but make sure to leave the shebang in place.

### 4.3. Add an import declaration before the workflow block

Insert the import declaration above the workflow block and fill it out appropriately.

_Before:_

```groovy title="hello-modules.nf" linenums="73"
// Include modules
include { sayHello } from './modules/sayHello.nf'
include { convertToUpper } from './modules/convertToUpper.nf'

workflow {
```

_After:_

```groovy title="hello-modules.nf" linenums="73"
// Include modules
include { sayHello } from './modules/sayHello.nf'
include { convertToUpper } from './modules/convertToUpper.nf'
include { collectGreetings } from './modules/collectGreetings.nf'

workflow {
```

### 4.4. Run the workflow to verify that it does the same thing as before

Run this with the `-resume` flag.

```bash
nextflow run hello-modules.nf -resume
```

This should still produce the same output as previously.

Congratulations, you've done all this work and absolutely nothing has changed to how the pipeline works!

Jokes aside, now your code is more modular, and if you decide to write another pipeline that calls on one of those processes, you just need to type one short import statement to use the relevant module.
This is better than just copy-pasting the code, because if later you decide to improve the module, all your pipelines will inherit the improvements.

### Takeaway

You know how to modularize multiple processes in a workflow.

### What's next?

Take a short break if you feel like it.
When you're ready, move on to Part 5 to learn how to manage inputs and parameters with more flexibility and convenience.
