# Part 2: Hello Channels

[TODO] CHANNELS ARE X. LEARN HOW TO USE THEM TO PROVIDE INPUTS TO PROCESSES

---

## 0. Warmup: Run Hello Channels

[TODO] VERIFY THAT IT RUNS

## 1. Add variable inputs using a channel

So far, we've been emitting a greeting hardcoded into the process command.
Now we're going to add some flexibility by using an input variable, so that we can easily change the greeting.

This requires us to make a series of inter-related changes:

1. Tell the process about expected variable inputs using the `input:` block
2. Edit the process to use the input
3. Create a **channel** to pass input to the process (more on that in a minute)
4. Add the channel as input to the process call

### 1.1. Add an input definition to the process block

First we need to adapt the process definition to accept an input.

_Before:_

```groovy title="hello-world.nf" linenums="6"
process sayHello {

    publishDir 'results', mode: 'copy'

    output:
        path 'output.txt'
```

_After:_

```groovy title="hello-world.nf" linenums="6"
process sayHello {

    publishDir 'results', mode: 'copy'

    input:
        val greeting

    output:
        path 'output.txt'
```

### 1.2. Edit the process command to use the input variable

Now we swap the original hardcoded value for the input variable.

_Before:_

```groovy title="hello-world.nf" linenums="16"
"""
echo 'Hello World!' > output.txt
"""
```

_After:_

```groovy title="hello-world.nf" linenums="16"
"""
echo '$greeting' > output.txt
"""
```

### 1.3. Create an input channel

Now that our process expects an input, we need to set up that input in the workflow body.
This is where channels come in: Nextflow uses channels to feed inputs to processes and ferry data between processes that are connected together.

There are multiple ways to do this, but for now, we're just going to use the simplest possible channel, containing a single value.

We're going to create the channel using the `of()` channel factory, which sets up a simple value channel, and give it a hardcoded string to use as greeting by declaring `greeting_ch = Channel.of('Hello world!')`.

_Before:_

```groovy title="hello-world.nf" linenums="21"
workflow {

    // emit a greeting
    sayHello()
}
```

_After:_

```groovy title="hello-world.nf" linenums="21"
workflow {

    // create a channel for inputs
    greeting_ch = Channel.of('Hello world!')

    // emit a greeting
    sayHello()
}
```

### 1.4. Add the channel as input to the process call

Now we need to actually plug our newly created channel into the `sayHello()` process call.

_Before:_

```groovy title="hello-world.nf" linenums="26"
// emit a greeting
sayHello()
```

_After:_

```groovy title="hello-world.nf" linenums="26"
// emit a greeting
sayHello(greeting_ch)
```

### 1.5. Run the workflow command again

Let's run it!

```bash
nextflow run hello-world.nf
```

If you made all four edits correctly, you should get another successful execution:

```console title="Output"
 N E X T F L O W   ~  version 24.10.0

 ┃ Launching `hello-world.nf` [prickly_avogadro] DSL2 - revision: b58b6ab94b

executor >  local (1)
[1f/50efd5] sayHello (1) [100%] 1 of 1 ✔
```

Feel free to check the results directory to satisfy yourself that the outcome is still the same as previously; so far we're just progressively tweaking the internal plumbing to increase the flexibility of our workflow while achieving the same end result.

### Takeaway

You know how to use a simple channel to provide an input to a process.

### What's next?

Learn [TODO]

---

## 2. Modify the workflow to run on a batch of input values

Workflows typically run on batches of inputs that are meant to be processed in bulk, so we want to upgrade the workflow to accept multiple input values.

### 2.1. Load multiple greetings into the input channel

Conveniently, the `of()` channel factory we've been using is quite happy to accept more than one value, so we don't need to modify that at all; we just have to load more values into the channel.

_Before:_

```groovy title="hello-world.nf" linenums="46"
// create a channel for inputs
greeting_ch = Channel.of('Hello')
```

_After:_

```groovy title="hello-world.nf" linenums="46"
// create a channel for inputs
greeting_ch = Channel.of('Hello','Bonjour','Holà')
```

The documentation tells us this should work. Can it really be so simple?

### 2.2. Run the command and look at the log output

Let's try it.

```bash
nextflow run hello-world.nf
```

Well, it certainly seems to run just fine.

```console title="Output"
 N E X T F L O W   ~  version 24.10.0

 ┃ Launching `hello-world.nf` [lonely_pare] DSL2 - revision: b9f1d96905

executor >  local (3)
[3d/1fe62c] sayHello (2)       [100%] 3 of 3 ✔
```

However... This seems to indicate that '3 of 3' calls were made for the process, which is encouraging, but this only give us one subdirectory path. What's going on?

By default, the ANSI logging system writes the logging from multiple calls to the same process on the same line. Fortunately, we can disable that behavior.

### 2.3. Run the command again with the `-ansi-log false` option

To expand the logging to display one line per process call, just add `-ansi-log false` to the command.

```bash
nextflow run hello-world.nf -ansi-log false
```

This time we see all six work subdirectories listed in the output:

```console title="Output"
N E X T F L O W  ~  version 24.02.0-edge
Launching `hello-world.nf` [big_woese] DSL2 - revision: 53f20aeb70
[62/d81e63] Submitted process > sayHello (1)
[19/507af3] Submitted process > sayHello (2)
[8a/3126e6] Submitted process > sayHello (3)
```

That's much better; at least for a simple workflow.
For a complex workflow, or a large number of inputs, having the full list output to the terminal might get a bit overwhelming.

That being said, we have another problem. If you look in the `results` directory, there is only one file: `output.txt`!

```console title="Directory contents"
results
└── output.txt
```

What's up with that? Shouldn't we be expecting a separate file per input greeting, so three files in all?
Did all three greetings go into a single file?
You can check the contents of `output.txt`; you will find only one of the three.

You may recall that we hardcoded the output file name for the first process.
This was fine as long as there was only a single call made per process, but when we start processing multiple input values and publishing the outputs into the same directory of results, it becomes a problem.
For a given process, every call produces an output with the same file name, so Nextflow just overwrites the previous output file every time a new one is produced.

### 2.4. Ensure the output file names will be unique

Since we're going to be publishing all the outputs to the same results directory, we need to ensure they will have unique names.
Specifically, we need to modify the first process to generate a file name dynamically so that the final file names will be unique.

So how do we make the file names unique? A common way to do that is to use some unique piece of metadata as part of the file name.
Here, for convenience, we'll just use the greeting itself.

_Before:_

```groovy title="hello-world.nf" linenums="11"
process sayHello {

    publishDir 'results', mode: 'copy'

    input:
        val greeting

    output:
        path 'output.txt'

    script:
    """
    echo '$greeting' > output.txt
    """
}
```

_After:_

```groovy title="hello-world.nf" linenums="11"
process sayHello {

    publishDir 'results', mode: 'copy'

    input:
        val greeting

    output:
        path "${greeting}-output.txt"

    script:
    """
    echo '$greeting' > $greeting-output.txt
    """
}
```

!!! tip

    You MUST use double quotes around the output filename expression (NOT single quotes), otherwise it will fail.

This should produce a unique output file name for every call of each process.

### 2.5. Run the workflow and look at the results directory

Let's run it and check that it works.

```bash
nextflow run hello-world.nf
```

Reverting back to the summary view, the output looks like this again:

```console title="Output"
 N E X T F L O W   ~  version 24.10.0

 ┃ Launching `hello-world.nf` [jovial_mccarthy] DSL2 - revision: 53f20aeb70

executor >  local (3)
[03/f007f2] sayHello (1)       [100%] 3 of 3 ✔
```

But more importantly, now we have three new files in addition to the one we already had in the `results` directory:

```console title="Directory contents"
results
├── Bonjour-output.txt
├── Hello-output.txt
├── Holà-output.txt
└── output.txt
```

Success! Now we can add as many greetings as we like without worrying about output files being overwritten.

!!! note

    In practice, naming files based on the input data itself is almost always impractical.
    The better way to generate dynamic filenames is to use a samplesheet contain relevant metadata (such as unique sample IDs) and create a data structure called a 'map', which we pass to processes, and from which we can grab an appropriate identifier to generate the filenames.
    We'll show you how to do that later in this training course.

### Takeaway

You know how to feed a batch of multiple input elements through a channel.

### What's next?

Learn how to make the workflow take a file as its source of input values.

---

## 3. Use CLI parameters to supply input values

We want to be able to specify the input from the command line, since that is the piece that will almost always be different in subsequent runs of the workflow.
Good news: Nextflow has a built-in workflow parameter system called `params`, which makes it easy to declare and use CLI parameters.

### 3.1. Edit the input channel declaration to use a parameter

Here we replace the hardcoded input strings with `params.greeting` in the channel creation line.

_Before:_

```groovy title="hello-world.nf" linenums="23"
// create a channel for inputs
greeting_ch = Channel.of('Hello','Bonjour','Holà')
```

_After:_

```groovy title="hello-world.nf" linenums="23"
// create a channel for inputs
greeting_ch = Channel.of(params.greeting)
```

This automatically creates a parameter called `greeting` that you can use to provide a value in the command line.

### 3.2. Run the workflow again with the `--greeting` parameter

To provide a value for this parameter, simply add `--greeting <value>` to your command line. Let's start with using a single value.

```bash
nextflow run hello-world.nf --greeting 'Bonjour le monde!'
```

Running this should feel extremely familiar by now.

```console title="Output"
 N E X T F L O W   ~  version 24.10.0

 ┃ Launching `hello-world.nf` [cheesy_engelbart] DSL2 - revision: b58b6ab94b

executor >  local (1)
[1c/9b6dc9] sayHello (1) [100%] 1 of 1 ✔
```

Be sure to open up the output file to check that you now have the new version of the greeting. Voilà!

!!! tip

    It's helpful to distinguish Nextflow-level parameters from pipeline-level parameters.
    For parameters that apply to a pipeline, we use a double hyphen (`--`), whereas we use a single hyphen (`-`) for parameters that modify a specific Nextflow setting, _e.g._ the `-resume` feature we used earlier.

### 3.3. Set a default value for a command line parameter

In many cases, it makes sense to supply a default value for a given parameter so that you don't have to specify it for every run.

Let's initialize the `greeting` parameter with a default value by adding the parameter declaration before the workflow definition (with a comment block as a free bonus).

```groovy title="hello-world.nf" linenums="3"
/*
 * Pipeline parameters
 */
params.greeting = 'Holà mundo!'
```

!!! tip

    You can put the parameter declaration inside the workflow block if you prefer. Whatever you choose, try to group similar things in the same place so you don't end up with declarations all over the place.

### 3.4. Run the workflow again without specifying the parameter

Now that you have a default value set, you can run the workflow again without having to specify a value in the command line.

```bash
nextflow run hello-world.nf
```

The console output should look the same.

```console title="Output"
 N E X T F L O W   ~  version 24.10.0

 ┃ Launching `hello-world.nf` [wise_waddington] DSL2 - revision: 988fc779cf

executor >  local (1)
[c0/8b8332] sayHello (1) [100%] 1 of 1 ✔
```

Check the output in the results directory, and... Tadaa! It works!
Nextflow used the default value to name the output.

!!! note

    If you provide the parameter on the command line, the CLI value will override the default value. Feel free to test this out.

    ```bash
    nextflow run hello-world.nf --greeting 'Konnichiwa!'
    ```

    In Nextflow, there are multiple places where you can specify values for parameters.
    If the same parameter is set to different values in multiple places, Nexflow will determine what value to use based on the order of precedence that is described [here](https://www.nextflow.io/docs/latest/config.html).

### Takeaway

You know how to use CLI parameters to feed inputs to the workflow.

### What's next?

Learn how to make the workflow take a file as its source of input values.

---

## 4. Supply a batch of multiple values via the `params` system

We sneakily reverted to running on just one value there.
What if we want to run on a batch again, like we did earlier?

Common sense suggests we should be able to simply pass in an array of values instead of a single value. Right?

### 4.1. Switch the `params.greeting` value to an array of values

[TODO]

_Before:_

```groovy title="hello-world.nf" linenums="23"
/*
 * Pipeline parameters
 */
params.greeting = 'Holà mundo'
```

_After:_

```groovy title="hello-world.nf" linenums="23"
/*
 * Pipeline parameters
 */
params.greeting = ['Holà mundo','Konnichiwa,','Dobrý den']
```

### 4.2. Run the workflow

[TODO] OH NO IT DOES NOT WORK. SHOW ERROR. JUST RUNS ONCE, TRIES TO USE THE WHOLE ARRAY AS A SINGLE INPUT. NEED TO TRANSFORM HOW CONTENTS ARE ORGANIZED/PACKAGED IN THE CHANNEL

### 4.3. Use the `flatten()` operator

[TODO] INTRODUCE CONCEPT OF OPERATORS. "You can think of them as ways of transforming the contents of a channel in a variety of ways."

[TODO] LOOK AT DOCS, FIND FLATTEN. ADD TO CHANNEL CONSTRUCTION LIKE THIS

_Before:_

```groovy title="hello-world.nf" linenums="46"
// create a channel for inputs
greeting_ch = Channel.of(params.greeting)
```

_After:_

```groovy title="hello-world.nf" linenums="46"
// create a channel for inputs
greeting_ch = Channel.of(params.greeting)
                     .flatten()
```

### 4.4. Add `view()` to inspect channel contents [TODO]

We can inspect how each operator changes how the contents of a channel are organized using the `.view()` operator:

_Before:_

```groovy title="hello-world.nf" linenums="46"
// create a channel for inputs
greeting_ch = Channel.of(params.greeting)
                     .flatten()
```

_After:_

```groovy title="hello-world.nf" linenums="46"
// create a channel for inputs
greeting_ch = Channel.of(params.greeting)
                     .view{ "Before flatten: $it" }
                     .flatten()
                     .view{ "After flatten: $it" }
```

### 4.5. Run the workflow [TODO]

[TODO] YAY IT WORKS AND ALSO THE VIEW STATEMENTS SHOW US WHAT'S HAPPENING

As you can see from the view() statements, the `flatten()` operator has transformed the channel from containing arrays to containing individual elements. This can be useful when you want to process each item separately in your workflow.

!!! tip

    You can delete or comment out the `view()` statements before moving on.

    ```groovy title="hello-world.nf" linenums="46"
    // create a channel for inputs
    greeting_ch = Channel.of(params.greeting)
                         .flatten()
    ```

### Takeaway

You know how to use the flatten() operator to handle a batch of values passed in through the CLI parameter system, and use the view() directive to inspect channel contents before and after applying the operators.

### What's next?

Learn how to make the workflow take a file as its source of input values.

---

## 5. Modify the workflow to take a file as its source of input values

It's often the case that, when we want to run on a batch of multiple input elements, the input values are contained in a file.
As an example, we have provided you with a CSV file called `greetings.csv` in the `data/` directory, containing several greetings separated by commas.

```csv title="greetings.csv"
Hello,Bonjour,Holà
```

So we need to modify our workflow to read in the values from a file like that.

### 5.1. Switch the `params.greeting` to the CSV file

_Before:_

```groovy title="hello-world.nf" linenums="23"
/*
 * Pipeline parameters
 */
params.greeting = ['Holà mundo','Konnichiwa,','Dobrý den']
```

_After:_

```groovy title="hello-world.nf" linenums="23"
/*
 * Pipeline parameters
 */
params.greeting = 'data/greetings.csv'
```

### 5.1. Update the channel declaration to use the input file

Since we now want to use a file instead of a simple value as the input, we can't use the `of()` channel factory from before.
We need to switch to using a new channel factory, `fromPath()`, which has some built-in functionality for handling file paths.

_Before:_

```groovy title="hello-world.nf" linenums="46"
// create a channel for inputs
greeting_ch = Channel.of(params.greeting)
                     .flatten()
```

_After:_

```groovy title="hello-world.nf" linenums="46"
// create a channel for inputs from a CSV file
greeting_ch = Channel.fromPath(params.greeting)
```

### 5.2. Run it [TODO]

[TODO] OH NO IT BREAKS. SHOW ERROR. JUST RUNS ONCE, TRIES TO USE THE PATH ITSELF AS GREETING. NOT WHAT WE WANT. WE WANT TO READ IN THE CONTENTS OF THE FILE. SOUNDS LIKE WE NEED ANOTHER OPERATOR!

### 5.3. Add `splitCsv()` operator [TODO]

[TODO] LOOK AT OPERATOR DOCS, FIND SPLITCSV: "an 'operator' to transform that CSV file into channel contents".

To apply the operator, add it to the channel construction instruction like previously; and we're also going to include view statements while we're at it.

_Before:_

```groovy title="hello-world.nf" linenums="46"
// create a channel for inputs from a CSV file
greeting_ch = Channel.fromPath(params.greeting)
```

_After:_

```groovy title="hello-world.nf" linenums="46"
// create a channel for inputs from a CSV file
greeting_ch = Channel.fromPath(params.greeting)
                     .view{ "Before splitCsv: $it" }
                     .splitCsv()
                     .view{ "After splitCsv: $it" }
```

### 5.4. Run it [TODO]

[TODO] OH COME ON IT BREAKS AGAIN. SHOW ERROR. CAN PROBABLY ALREADY GUESS WHAT THE PROBLEM IS BUT HEY LET'S CHECK THOSE VIEW STATEMENTS. OH SEE, BRACKETS. CONFIRMS IT PASSED ALL THE ITEMS TOGETHER AS ONE ARRAY ELEMENT (INDICATED BY BRACKETS). BRACKETS IN THE OUTPUT FILE BREAK THE ECHO COMMAND. EVEN IF IT DIDN'T, THIS IS STILL NOT WHAT WE WANT. WE WANT TO BREAK UP THE PACKAGE FOR THE GREETINGS TO BE USED AS SEPARATE INPUT ITEMS.

### 5.5. Add `flatten()` operator [TODO]

[TODO] REMEMBER FLATTEN? WE LOVE FLATTEN

To apply the operator, add it to the channel construction instruction. Include another view() call.

_Before:_

```groovy title="hello-world.nf" linenums="46"
// create a channel for inputs from a CSV file
greeting_ch = Channel.fromPath(params.greeting)
                     .view{ "Before splitCsv: $it" }
                     .splitCsv()
                     .view{ "After splitCsv: $it" }
```

_After:_

```groovy title="hello-world.nf" linenums="46"
// create a channel for inputs from a CSV file
greeting_ch = Channel.fromPath(params.greeting)
                     .view{ "Before splitCsv: $it" }
                     .splitCsv()
                     .view{ "After splitCsv: $it" }
                     .flatten()
                     .view{ "After flatten: $it" }
```

### 5.6. Run it [TODO]

[TODO] THIS TIME IT WORKS, YAY

Looking at the outputs, we see each greeting was correctly extracted and processed through the workflow. We've achieved the same result as previously, but now we have a lot more flexibility to add more elements to the channel of greetings we want to process without modifying any code.

[TODO] NOTE THAT IF YOU ADD MORE LINES TO THE CSV EVERYTHING GETS PARSED AS INDIVIDUAL ITEMS. THAT'S WHAT FLATTEN IS DOING. LEARN MORE ABOUT PLUMBING LATER.

!!! note

    Be sure to remove the `.view()` operations before you continue.

    ```groovy title="hello-world.nf" linenums="46"
    // create a channel for inputs from a CSV file
    greeting_ch = Channel.fromPath('data/greetings.csv')
                         .splitCsv()
                         .flatten()
    ```

### Takeaway

You know how to use the splitCsv() and flatten() operators to handle a batch of values passed in through a file.

More generally, [TODO]

### What's next?

Take a break!

Don't worry if the channel factories and operators feel like a lot to grapple with the first time you encounter them.
You'll get more opportunities to practice using these components in various settings as you work through this training course.

When you're ready, move on to [TODO] (plumbing multiple processes)
