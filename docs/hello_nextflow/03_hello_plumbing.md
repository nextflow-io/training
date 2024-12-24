# Part 3: Hello Plumbing

[TODO]

---

## 0. Warmup: Run in terminal

[TODO]
Most real-world workflows involve more than one step. Here we introduce a second process that converts the text to uppercase (all-caps), using the classic UNIX one-liner:

```bash
tr '[a-z]' '[A-Z]'
```

We're going to run the command by itself in the terminal first to verify that it works as expected without any of the workflow code getting in the way of clarity, just like we did at the start with `echo 'Hello World'`.

### 0.1. Run the command in the terminal by itself

```bash
echo 'Hello World' | tr '[a-z]' '[A-Z]'
```

The output is simply the uppercase version of the text string:

```console title="Output"
HELLO WORLD
```

!!! note

    This is a very naive text replacement one-liner that does not account for accented letters, so for example 'Holà' will become 'HOLà'. This is expected.

### 0.2. Make the command take a file as input and write the output to a file

As previously, we want to output results to a dedicated file, which we name by prepending the original filename with `UPPER-`.

```bash
cat output.txt | tr '[a-z]' '[A-Z]' > UPPER-output.txt
```

Now the `HELLO WORLD` output is in the new output file, `UPPER-output.txt`.

---

## 1. Add a second step to the workflow

[TODO] INTRO

We're going to write a process that wraps the command we just ran in the terminal, then we'll add it to the workflow, setting it up to take the output of the `sayHello()` process as input.

### 1.1. Wrap the command in a new Nextflow process definition

We can model our new process on the first one, since we want to use all the same components.

```groovy title="hello-plumbing.nf" linenums="26"
/*
 * Use a text replace utility to convert the greeting to uppercase
 */
process convertToUpper {

    publishDir 'results', mode: 'copy'

    input:
        path input_file

    output:
        path "UPPER-${input_file}"

    script:
    """
    cat '$input_file' | tr '[a-z]' '[A-Z]' > UPPER-${input_file}
    """
}
```

Similarly to what we did for the output of the first process, we compose the second output filename based on the input filename.

!!! tip

    You MUST use double quotes around the output filename expression (NOT single quotes), otherwise it will fail.

### 1.2. Add a call to the new process in the workflow body

Don't forget we need to tell Nextflow to actually call the process we just created! To do that, we add it to the `workflow` body.

```groovy title="hello-plumbing.nf" linenums="44"
workflow {

    // create a channel for inputs from a CSV file
    greeting_ch = Channel.fromPath('data/greetings.csv')
                         .splitCsv()
                         .flatten()

    // emit a greeting
    sayHello(greeting_ch)

    // convert the greeting to uppercase
    convertToUpper()
}
```

Looking good! Now we just need to wire up the `convertToUpper` process call to run on the output of `sayHello`.

### 1.3. Pass the output of the first process to the second process

The output of the `sayHello` process is automatically provided as a channel called `sayHello.out`, so all we need to do is pass that as the input to the `convertToUpper` process.

```groovy title="hello-plumbing.nf" linenums="52"
// convert the greeting to uppercase
convertToUpper(sayHello.out)
```

For a simple case like this, that's all we need to do to connect two processes!

### 1.4. Run the same workflow command as before

Let's make sure this works:

```bash
nextflow run hello-plumbing.nf
```

Oh, how exciting! There is now an extra line in the log output, which corresponds to the new process we just added:

```console title="Output"
 N E X T F L O W   ~  version 24.10.0

 ┃ Launching `hello-plumbing.nf` [magical_brenner] DSL2 - revision: 0e18f34798

executor >  local (2)
[57/3836c0] sayHello (1)       [100%] 1 of 1 ✔
[ee/bb3cc8] convertToUpper (1) [100%] 1 of 1 ✔
```

You'll notice that this time the workflow produced two new work subdirectories; one per process call.
Check out the work directory of the call to the second process, where you should find two different output files listed. If you look carefully, you'll notice one of them (the output of the first process) has a little arrow icon on the right; that signifies it's a symbolic link.
It points to the location where that file lives in the work directory of the first process.
By default, Nextflow uses symbolic links to stage input files whenever possible, to avoid making duplicate copies.

!!! note

    All we did was connect the output of `sayHello` to the input of `convertToUpper` and the two processes could be run serially.
    Nextflow did the hard work of handling input and output files and passing them between the two commands for us.
    This is the power of channels in Nextflow, doing the busywork of connecting our pipeline steps together.

    What's more, Nextflow will automatically determine which call needs to be executed first based on how they're connected, so the order in which they're written in the workflow body does not matter.
    However, we do recommend you be kind to your collaborators and to your future self, and try to write them in a logical order!

### Takeaway

You know how to add a second step that takes the output of the first step as input.

### What's next?

[TODO]

---

## 2. [TODO] ADD A COLLECT STEP?

## 3. [TODO] MAKE THE COLLECT OPTIONAL?
