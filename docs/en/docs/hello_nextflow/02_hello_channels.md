# Part 2: Hello Channels

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/yDR66fzAMOg?si=y8lAedhEHWaTV4zd&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } See [the whole playlist](https://youtube.com/playlist?list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&si=eF7cLR62goy-lc6n) on the Nextflow YouTube channel.

:green_book: The video transcript is available [here](./transcripts/02_hello_channels.md).
///

In Part 1 of this course (Hello World), we showed you how to provide a variable input to a process by providing the input in the process call directly: `sayHello(params.name)`.
That was a deliberately simplified approach.
In practice, that approach has major limitations; namely that it only works for very simple cases where we only want to run the process once, on a single value.
In most realistic workflow use cases, we want to process multiple values (experimental data for multiple samples, for example), so we need a more sophisticated way to handle inputs.

That is what Nextflow [**channels**](https://nextflow.io/docs/latest/channel.html) are for.
Channels are queues designed to handle inputs efficiently and shuttle them from one step to another in multi-step workflows, while providing built-in parallelism and many additional benefits.

In this part of the course, you will learn how to use a channel to handle multiple inputs from a variety of different sources.
You will also learn to use [**operators**](https://nextflow.io/docs/latest/reference/operator.html) to transform channel contents as needed.

??? info "How to begin from this section"

    This section of the course assumes you have completed Part 1 of the [Hello Nextflow](./index.md) course, but if you are comfortable with the basics covered in that section, you can start from here without doing anything special.

---

## 0. Warmup: Run `hello-channels.nf`

We're going to use the workflow script `hello-channels.nf` as a starting point.
It is equivalent to the script produced by working through Part 1 of this training course.

Just to make sure everything is working, run the script once before making any changes:

```bash
nextflow run hello-channels.nf --name 'Channels'
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `hello-channels.nf` [fabulous_lamarck] revision: 07609cbd83

    executor >  local (1)
    [c7/c5367c] sayHello | 1 of 1 ✔

    Outputs:

      /workspaces/training/hello-nextflow/results

      greetings: output.txt
    ```

As previously, you will find the output file named `output.txt` in the `results` directory.

??? abstract "Directory contents"

    ```console title="results"
    results
    └── output.txt -> /workspaces/training/hello-nextflow/work/c7/c5367c.../output.txt
    ```

??? abstract "File contents"

    ```console title="results/output.txt"
    Hello, Channels!
    ```

If that worked for you, you're ready to learn about channels.

---

## 1. Provide variable inputs via a channel explicitly

We are going to create a **channel** to pass the variable input to the `sayHello()` process instead of relying on the implicit handling, which has certain limitations.

### 1.1. Create an input channel

There are a variety of [**channel factories**](https://nextflow.io/docs/latest/reference/channel.html) that we can use to set up a channel.
To keep things simple for now, we are going to use the most basic channel factory, called [`channel.of`](https://nextflow.io/docs/latest/reference/channel.html#of), which will create a channel containing a single value.
Functionally this will be similar to how we had it set up before, but instead of having Nextflow create a channel implicitly, we are doing this explicitly now.

This is the line of code we're going to use:

```console title="Syntax"
name_ch = channel.of('Mateo')
```

This creates a channel called `name_ch` using the `channel.of()` channel factory, which sets up a simple queue channel, and loads the string `'Mateo'` to use as the name value.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-pipeline-channel.svg"
</figure>

!!! note

    We are temporarily switching back to hardcoded strings instead of using a CLI parameter for the sake of readability. We'll go back to using CLI parameters once we've covered what's happening at the level of the channel.

In the workflow block, add the channel factory code:

=== "After"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="4 5"
    workflow {

        main:
        // create a channel for inputs
        name_ch = channel.of('Mateo')
        // emit a greeting
        sayHello(params.name)

        publish:
        greetings = sayHello.out
    }
    ```

=== "Before"

    ```groovy title="hello-channels.nf" linenums="27"
    workflow {

        main:
        // emit a greeting
        sayHello(params.name)

        publish:
        greetings = sayHello.out
    }
    ```

This is not yet functional since we haven't yet switched the input to the process call.

### 1.2. Add the channel as input to the process call

Now we need to actually plug our newly created channel into the `sayHello()` process call, replacing the CLI parameter which we were providing directly before.

In the workflow block, make the following code change:

=== "After"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="7"
    workflow {

        main:
        // create a channel for inputs
        name_ch = channel.of('Mateo')
        // emit a greeting
        sayHello(name_ch)

        publish:
        greetings = sayHello.out
    }
    ```

=== "Before"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="7"
    workflow {

        main:
        // create a channel for inputs
        name_ch = channel.of('Mateo')
        // emit a greeting
        sayHello(params.name)

        publish:
        greetings = sayHello.out
    }
    ```

This tells Nextflow to run the `sayHello` process on the contents of the `name_ch` channel.

Now our workflow is properly functional; it is the explicit equivalent of writing `sayHello('Mateo')`.

### 1.3. Run the workflow

Let's run it!

```bash
nextflow run hello-channels.nf
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `hello-channels.nf` [voluminous_leakey] revision: 857551c2e8

    executor >  local (1)
    [39/e7facf] sayHello (1) | 1 of 1 ✔

    Outputs:

      /workspaces/training/hello-nextflow/results

      greetings:
        - output.txt
    ```

If you made both edits correctly, you should get a successful execution.
You can check the results directory to see the greeting for the name we loaded into the channel.

??? abstract "File contents"

    ```console title="results/output.txt"
    Hello, Mateo!
    ```

So we've increased the flexibility of our workflow while keeping the same kind of result.
This may seem like we're writing more code for no tangible benefit, but the value will become clear as soon as we start handling more inputs.

As a preview of that, let's look at one more thing before we move on: one small but convenient benefit of using an explicit channel to manage data input.

### 1.4. Use `view()` to inspect the channel contents

Nextflow channels are built in a way that allows us to operate on their contents using operators, which we'll cover in detail later in this chapter.

For now, we're just going to show you how to use a super simple operator called [`view()`](https://www.nextflow.io/docs/latest/reference/operator.html#view) to inspect the contents of a channel.
You can think of `view()` as a debugging tool, like a `print()` statement in Python, or its equivalent in other languages.

Add this tiny line to the workflow block:

=== "After"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="6"
    workflow {

        main:
        // create a channel for inputs
        name_ch = channel.of('Mateo')
                            .view()
        // emit a greeting
        sayHello(name_ch)

        publish:
        greetings = sayHello.out
    }
    ```

=== "Before"

    ```groovy title="hello-channels.nf" linenums="27"
    workflow {

        main:
        // create a channel for inputs
        name_ch = channel.of('Mateo')
        // emit a greeting
        sayHello(name_ch)

        publish:
        greetings = sayHello.out
    }
    ```

The exact amount of whitespace doesn't matter; we're just lining up the `.view()` statement under the channel construction to make the code easier to read.

Now run the workflow again:

```bash
nextflow run hello-channels.nf
```

??? success "Command output"

    ```console hl_lines="7"
    N E X T F L O W   ~  version 26.04.4

    Launching `hello-channels.nf` [stupefied_heisenberg] revision: a7a1ebd855

    executor >  local (1)
    [41/a766f3] sayHello (1) | 1 of 1 ✔
    Mateo

    Outputs:

      /workspaces/training/hello-nextflow/results

      greetings:
        - output.txt
    ```

As you can see, this outputs the channel contents to the console.
Here we only have one element, but when we start loading multiple values into the channel in the next section, you'll see that this is set to output one element per line.

### Takeaway

You know how to use a basic channel factory to provide an input to a process.

### What's next?

Learn how to use channels to make the workflow iterate over multiple input values.

---

## 2. Modify the workflow to run on multiple input values

Workflows typically run on batches of inputs that are meant to be processed in bulk, so we want to upgrade the workflow to accept multiple input values.

### 2.1. Load multiple names into the input channel

Conveniently, the `channel.of()` channel factory we've been using is quite happy to accept more than one value, so we don't need to modify that at all.
We can just load multiple values into the channel.

Let's make them `'Mateo'`, `'Andrea'` and `'Celine'`.

#### 2.1.1. Add more names

In the workflow block, make the following code change:

=== "After"

    ```groovy title="hello-channels.nf" linenums="30" hl_lines="2"
    // create a channel for inputs
    name_ch = channel.of('Mateo','Andrea','Celine')
                        .view()
    ```

=== "Before"

    ```groovy title="hello-channels.nf" linenums="30" hl_lines="2"
    // create a channel for inputs
    name_ch = channel.of('Mateo')
                        .view()
    ```

The documentation tells us this should work. Can it really be so simple?

#### 2.1.2. Run the command and look at the log output

Let's try it.

```bash
nextflow run hello-channels.nf
```

??? success "Command output"

    ```console hl_lines="6-9"
    N E X T F L O W   ~  version 26.04.4

    Launching `hello-channels.nf` [determined_hoover] revision: cdd5bd6ff4

    executor >  local (3)
    [a7/df0f53] sayHello (3) | 3 of 3 ✔
    Mateo
    Andrea
    Celine

    Outputs:

      /workspaces/training/hello-nextflow/results

      greetings:
        - output.txt
        - output.txt
        - output.txt
    ```

It certainly seems to have run just fine.
The execution monitor shows that `3 of 3` calls were made for the `sayHello` process, and we see the three names enumerated by the `view()` statement, one per line as promised.

However, there is still only one output in the results directory:

??? abstract "Directory contents"

    ```console title="results"
    results
    └── output.txt -> /workspaces/training/hello-nextflow/work/.../output.txt
    ```

??? abstract "File contents"

    ```console title="results/output.txt"
    Hello, Celine!
    ```

You should see the greeting for one of the three names in there, though the one you got might be different from what is shown here.
Can you think of why that might be?

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-pipeline-channel-multi.svg"
</figure>

_In the diagram, the channel is represented in green, and the order of elements is represented like marbles in a pipe: the first one loaded is on the right, then the second one in the middle, then the third is on the left._

Looking back at the execution monitor, it gave us only one subdirectory path (`a7/df0f53`).
Let's have a look in there.

??? abstract "Directory contents"

    ```console hl_lines="9"
    work/a7/df0f53dbc94195e686e4222b3baa53/
    ├── .command.begin
    ├── .command.err
    ├── .command.log
    ├── .command.out
    ├── .command.run
    ├── .command.sh
    ├── .exitcode
    └── output.txt
    ```

??? abstract "File contents"

    ```console title="work/a7/df0f53dbc94195e686e4222b3baa53/output.txt"
    Hello, Celine!
    ```

That's one of the three greetings, but where are the other two? What's going on?

At this point, we need to tell you that by default, the ANSI logging system writes the logging from multiple calls to the same process on the same line.
So the status from all three calls to the `sayHello()` process are landing in the same spot.

Fortunately, we can disable that behavior to see the full list of process calls.

#### 2.1.3. Run the command again with the `-ansi-log false` option

To expand the logging to display one line per process call, add `-ansi-log false` to the command.

```bash
nextflow run hello-channels.nf -ansi-log false
```

??? success "Command output"

    ```console
    N E X T F L O W  ~  version 26.04.4
    Launching `hello-channels.nf` [amazing_mendel] - revision: cdd5bd6ff4
    Mateo
    Andrea
    Celine
    [c3/557357] Submitted process > sayHello (1)
    [24/a64e64] Submitted process > sayHello (2)
    [3c/1ef83b] Submitted process > sayHello (3)

    Outputs:

      /workspaces/training/hello-nextflow/results

      greetings:
        - output.txt
        - output.txt
        - output.txt
    ```

This time we see all three process runs and their associated work subdirectories listed in the output.

That's much better, at least for a simple workflow.
For a complex workflow, or a large number of inputs, having the full list output to the terminal would get a bit overwhelming.
That's why `-ansi-log false` is not the default behavior.

!!! tip

    The way the status is reported is a bit different between the two logging modes.
    In the condensed mode, Nextflow reports whether calls were completed successfully or not.
    In this expanded mode, it only reports that they were submitted.

Anyway, now that we have the subdirectories of each process call, we can look for the their logs and outputs.

??? abstract "Directory contents"

    ```console
    work/c3/5573571b6ea4e9c12bd95edf56fc94/
    ├── .command.begin
    ├── .command.err
    ├── .command.log
    ├── .command.out
    ├── .command.run
    ├── .command.sh
    ├── .exitcode
    └── output.txt
    ```

    ```console
    work/24/a64e64785030d9c2c7731ffcfa2469/
    ├── .command.begin
    ├── .command.err
    ├── .command.log
    ├── .command.out
    ├── .command.run
    ├── .command.sh
    ├── .exitcode
    └── output.txt
    ```

    ```console
    work/3c/1ef83bfb5e198624a40f13a767eef6/
    ├── .command.begin
    ├── .command.err
    ├── .command.log
    ├── .command.out
    ├── .command.run
    ├── .command.sh
    ├── .exitcode
    └── output.txt
    ```

??? abstract "File contents"

    ```txt title="work/c3/5573571b6ea4e9c12bd95edf56fc94/output.txt"
    Hello, Mateo!
    ```

    ```txt title="work/24/a64e64785030d9c2c7731ffcfa2469/output.txt"
    Hello, Andrea!
    ```

    ```txt title="work/3c/1ef83bfb5e198624a40f13a767eef6/output.txt"
    Hello, Celine!
    ```

This shows that all three processes ran successfully (yay).

That being said, we still have the problem that there is only one output file in the results directory.

You may recall that we hardcoded the output file name for the `sayHello` process, so all three calls produced a file called `output.txt`.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-channels-task-dirs.svg"
</figure>

As long as the output files stay in the work subdirectories, isolated from the other processes, that is okay.
But when they are published to the same results directory, whichever got published there first gets overwritten by the next one, and so on.

### 2.2. Ensure the output file names will be unique

We can continue publishing all the outputs to the same results directory, but we need to ensure they will have unique names.
Specifically, we need to modify the first process to generate a file name dynamically so that the final file names will be unique.

So how do we make the file names unique?
A common way to do that is to use some unique piece of metadata from the inputs (received from the input channel) as part of the output file name.
Here, for convenience, we'll just use the name itself since it's just a short string, and prepend it to the base output filename.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-pipeline-channel-multi-unique.svg"
</figure>

#### 2.2.1. Construct a dynamic output file name

In the process block, make the following code changes:

=== "After"

    ```groovy title="hello-channels.nf" linenums="6" hl_lines="7 11"
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

=== "Before"

    ```groovy title="hello-channels.nf" linenums="6" hl_lines="7 11"
    process sayHello {

        input:
        val name

        output:
        path 'output.txt'

        script:
        """
        echo 'Hello, ${name}!' > output.txt
        """
    }
    ```

Make sure to replace `output.txt` in both the output definition and in the `script:` command block.

!!! tip

    In the output definition, you MUST use double quotes around the output filename expression (NOT single quotes), otherwise it will fail.

This should produce a unique output file name every time the process is called, so that it can be distinguished from the outputs from other calls to the same process in the output directory.

#### 2.2.2. Run the workflow

Let's run it. Note that we're back to running with the default ANSI log settings.

```bash
nextflow run hello-channels.nf
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `hello-channels.nf` [grave_pare] revision: 80650ff4ed

    executor >  local (3)
    [12/0c7a46] sayHello (2) | 3 of 3 ✔
    Mateo
    Andrea
    Celine

    Outputs:

      /workspaces/training/hello-nextflow/results

      greetings:
        - Mateo.txt
        - Celine.txt
        - Andrea.txt
    ```

Reverting back to the summary view, the output is summarized on one line again.
Have a look at the `results` directory to see if all the greetings are there.

??? abstract "Directory contents"

    ```console
    results
    ├── Andrea.txt -> /workspaces/training/hello-nextflow/work/.../Andrea.txt
    ├── Celine.txt -> /workspaces/training/hello-nextflow/work/.../Celine.txt
    ├── Mateo.txt -> /workspaces/training/hello-nextflow/work/.../Mateo.txt
    └── output.txt -> /workspaces/training/hello-nextflow/work/.../output.txt
    ```

??? abstract "File contents"

    ```console title="Andrea.txt"
    Hello, Andrea!
    ```

    ```console title="Celine.txt"
    Hello, Celine!
    ```

    ```console title="Mateo.txt"
    Hello, Mateo!
    ```

Success! Now we can add as many names as we like without worrying about output files being overwritten.

!!! tip

    In practice, naming files based on the input data itself is almost always impractical.
    The better way to generate dynamic filenames is to pass metadata to a process along with the input files.
    The metadata is typically provided via a 'sample sheet' or equivalents.
    You'll learn how to do that later in your Nextflow training (see [Metadata side quest](../side_quests/metadata/index.md)).

### Takeaway

You know how to feed multiple input elements through a channel.

### What's next?

Learn to use an operator to transform the contents of a channel.

---

## 3. Provide multiple inputs via an array

We just showed you how to handle multiple input elements that were hardcoded directly in the channel factory.
What if we wanted to provide those multiple inputs in a different way?

For example, imagine we set up an input variable containing an array of elements like this:

`names_array = ['Mateo','Andrea','Celine']`

Can we load that into our output channel and expect it to work?

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-pipeline-multi-inputs-array.svg"
</figure>

Let's find out.

### 3.1. Provide an array of values as input to the channel

Common sense suggests we should be able to simply pass in an array of values instead of a single value.
Let's try it; we'll need to set up the input variable and load it into the channel factory.

#### 3.1.1. Set up the input variable

Let's take the `names_array` variable we just imagined and make it a reality by adding it to the workflow block:

=== "After"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="4 5"
    workflow {

        main:
        // declare an array of input names
        names_array = ['Mateo','Andrea','Celine']
        // create a channel for inputs
        name_ch = channel.of('Mateo','Andrea','Celine')
                            .view()
        // emit a greeting
        sayHello(name_ch)

        publish:
        greetings = sayHello.out
    }
    ```

=== "Before"

    ```groovy title="hello-channels.nf" linenums="27"
    workflow {

        main:
        // create a channel for inputs
        name_ch = channel.of('Mateo','Andrea','Celine')
                            .view()
        // emit a greeting
        sayHello(name_ch)

        publish:
        greetings = sayHello.out
    }
    ```

This is not yet functional, we've just added a declaration for the array.

#### 3.1.2. Set array of names as the input to the channel factory

Now we're going to replace the values `'Mateo','Andrea','Celine'` currently hardcoded in the channel factory with the `names_array` we just created.

In the workflow block, make the following change:

=== "After"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="7"
    workflow {

        main:
        // declare an array of input names
        names_array = ['Mateo','Andrea','Celine']
        // create a channel for inputs
        name_ch = channel.of(names_array)
                            .view()
        // emit a greeting
        sayHello(name_ch)

        publish:
        greetings = sayHello.out
    }
    ```

=== "Before"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="7"
    workflow {

        main:
        // declare an array of input names
        names_array = ['Mateo','Andrea','Celine']
        // create a channel for inputs
        name_ch = channel.of('Mateo','Andrea','Celine')
                            .view()
        // emit a greeting
        sayHello(name_ch)

        publish:
        greetings = sayHello.out
    }
    ```

This should be functional now.

#### 3.1.3. Run the workflow

Let's try running it:

```bash
nextflow run hello-channels.nf
```

??? failure "Command output"

    ```console hl_lines="7 11 16"
    N E X T F L O W   ~  version 26.04.4

    Launching `hello-channels.nf` [happy_dijkstra] revision: 5539b2ee39

    executor >  local (1)
    [d2/7d5016] sayHello (1) | 0 of 1 ✘
    [Mateo, Andrea, Celine]
    ERROR ~ Error executing process > 'sayHello (1)'

    Caused by:
      Missing output file(s) `[Mateo, Andrea, Celine].txt` expected by process `sayHello (1)`


    Command executed:

      echo 'Hello, [Mateo, Andrea, Celine]!' > '[Mateo, Andrea, Celine].txt'

    Command exit status:
      0

    Command output:
      (empty)

    Work dir:
      /workspaces/training/hello-nextflow/work/d2/7d5016c135c10e18df126e658ac26b

    Tip: you can replicate the issue by changing to the process work dir and entering the command `bash .command.run`

    -- Check '.nextflow.log' file for details
    ```

Oh no! There's an error!

Look at the output of `view()` and the error messages.

It looks like Nextflow tried to run a single process call, using `[Mateo, Andrea, Celine]` as a single string value, instead of using the three strings in the array as separate values.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-channels-array-fail.svg"
</figure>

So it's the 'packaging' that is causing the problem.
How do we get Nextflow to unpack the array and load the individual strings into the channel?

### 3.2. Use an operator to transform channel contents

This is where [**operators**](https://nextflow.io/docs/latest/reference/operator.html) come into play.
You've already used the `.view()` operator, which just looks at what's in there.
Now we're going to look at operators that allow us to act on the contents of a channel.

If you skim through the [list of operators](https://nextflow.io/docs/latest/reference/operator.html) in the Nextflow documentation, you'll find [`flatten()`](https://nextflow.io/docs/latest/reference/operator.html#flatten), which does exactly what we need: unpack the contents of an array and emit them as individual items.

#### 3.2.1. Add the `flatten()` operator

To apply the `flatten()` operator to our input channel, we append it to the channel factory declaration.

In the workflow block, make the following code change:

=== "After"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="9"
    workflow {

        main:
        // declare an array of input names
        names_array = ['Mateo','Andrea','Celine']
        // create a channel for inputs
        name_ch = channel.of(names_array)
                            .view()
                            .flatten()
        // emit a greeting
        sayHello(name_ch)

        publish:
        greetings = sayHello.out
    }
    ```

=== "Before"

    ```groovy title="hello-channels.nf" linenums="27"
    workflow {

        main:
        // declare an array of input names
        names_array = ['Mateo','Andrea','Celine']
        // create a channel for inputs
        name_ch = channel.of(names_array)
                            .view()
        // emit a greeting
        sayHello(name_ch)

        publish:
        greetings = sayHello.out
    }
    ```

Here we added the operator on the next line for readability, but you can add operators on the same line as the channel factory if you prefer, like this:
`name_ch = channel.of(names_array).view().flatten()`

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-channels-array-success.svg"
</figure>

#### 3.2.2. Refine the `view()` statement(s)

We could run this right away to test if it works, but while we're at it, we're going to refine how we inspect the channel contents.

We want to be able to contrast what the contents look like before and after the `flatten()` operator is applied, so we're going to add a second one, AND we're going to add a bit of code to get them labeled more clearly in the output.

In the workflow block, make the following code change:

=== "After"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="8-10"
    workflow {

        main:
        // declare an array of input names
        names_array = ['Mateo','Andrea','Celine']
        // create a channel for inputs
        name_ch = channel.of(names_array)
                            .view { names -> "Before flatten: $names" }
                            .flatten()
                            .view { name -> "After flatten: $name" }
        // emit a greeting
        sayHello(name_ch)

        publish:
        greetings = sayHello.out
    }
    ```

=== "Before"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="8"
    workflow {

        main:
        // declare an array of input names
        names_array = ['Mateo','Andrea','Celine']
        // create a channel for inputs
        name_ch = channel.of(names_array)
                            .view()
                            .flatten()
        // emit a greeting
        sayHello(name_ch)

        publish:
        greetings = sayHello.out
    }
    ```

You see we've added a second `.view` statement, and for each of them, we've replaced the empty parentheses (`()`) with curly braces containing some code, such as `#!groovy { name -> "After flatten: $name" }`.

These are called _closures_. The code they contain will be executed for each item in the channel.
We define a temporary variable for the inner value, here called `name`, which is only used within the scope of that closure.
The name is up to you, but it's good practice to name it after what we expect the channel to contain at that point: before `flatten()`, the channel contains the whole list of names, so the first closure calls it `names`.

In this example, `$name` represents each individual item loaded in the channel.
This will result in neatly labeled console output.

!!! info

    In some pipelines you may see a special variable called `$it` used inside operator closures.
    This is an _implicit_ variable that allows a short-hand access to the inner variable,
    without needing to define it with a `->`.

    We prefer to be explicit to aid code clarity, as such the `$it` syntax is discouraged and will slowly be phased out of the Nextflow language.

#### 3.2.3. Run the workflow

Finally, you can try running the workflow again!

```bash
nextflow run hello-channels.nf
```

??? success "Command output"

    ```console hl_lines="7-10"
    N E X T F L O W   ~  version 26.04.4

    Launching `hello-channels.nf` [stupefied_wing] revision: 4d592f7b48

    executor >  local (3)
    [25/82982d] sayHello (2) | 3 of 3 ✔
    Before flatten: [Mateo, Andrea, Celine]
    After flatten: Mateo
    After flatten: Andrea
    After flatten: Celine

    Outputs:

      /workspaces/training/hello-nextflow/results

      greetings:
        - Mateo.txt
        - Celine.txt
        - Andrea.txt
    ```

This time it works AND gives us the additional insight into what the contents of the channel look like before and after we run the `flatten()` operator.

- A single `Before flatten:` statement because at that point the channel contains one item, the original array.
- Three separate `After flatten:` statements, one for each name, which are now individual items in the channel.

Importantly, this means each item can now be processed separately by the workflow.

!!! tip

    It is technically possible to achieve the same results by using a different channel factory, [`channel.fromList`](https://nextflow.io/docs/latest/reference/channel.html#fromlist), which includes an implicit mapping step in its operation.
    Here we chose not to use that in order to demonstrate the use of an operator on a simple use case.

### Takeaway

You know how to use an operator like `flatten()` to transform the contents of a channel, and how to use the `view()` operator to inspect channel contents before and after applying an operator.

### What's next?

Learn how to make the workflow take a file as its source of input values.

---

## 4. Read input values from a CSV file

Realistically, we're rarely if ever going to start from an array of values.
Most likely, we'll have one or more files containing the data that needs to be processed, in some kind of structured format.

We've prepared a CSV file called `names.csv` under `data/`, mimicking the kind of columnar data you might want to process in a real data analysis.
Each row describes one person: their name, a two-letter language code and a cowpy character.

```csv title="data/names.csv" linenums="1"
Mateo,es,squirrel
Andrea,it,turtle
Celine,fr,moose
```

For now we only need the names, in the first column.
The other two columns aren't used in the main course; they're there for bonus exercises.

Our next task is to adapt our workflow to read in the values from this file.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-pipeline-multi-inputs-csv.svg"
</figure>

Let's see how we can make that happen.

### 4.1. Modify the script to expect a CSV file as the source of names

To get started, we're going to need to make two key changes to the script:

- Replace the `name` parameter with an `input` parameter that points to the CSV file
- Switch the channel factory to one designed to handle a file

#### 4.1.1. Replace the `name` parameter with an `input` parameter

In Part 1, we set up a `name` parameter to provide a single name on the command line.
The names now come from a file, so instead of reusing `name` for something different, we replace it with a new parameter called `input`, of type `Path`, which points to the CSV file.
That way, each parameter's name keeps describing what it holds.

Make the following edit to the parameter declaration:

=== "After"

    ```groovy title="hello-channels.nf" linenums="20" hl_lines="5"
    /*
     * Pipeline parameters
     */
    params {
        input: Path = 'data/names.csv'
    }
    ```

=== "Before"

    ```groovy title="hello-channels.nf" linenums="20" hl_lines="5"
    /*
     * Pipeline parameters
     */
    params {
        name: String = 'Friend'
    }
    ```

This assumes the file is in the `data/` directory relative to where you launch the workflow.
You'll learn how to deal with other data locations later in your Nextflow journey.

#### 4.1.2. Switch to a channel factory designed to handle a file

Since we now want to use a file instead of simple strings as the input, we can't use the `channel.of()` channel factory from before.
We need to switch to using a new channel factory, [`channel.fromPath()`](https://nextflow.io/docs/latest/reference/channel.html#frompath), which has some built-in functionality for handling file paths.

In the workflow block, make the following code change:

=== "After"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="4-8"
    workflow {

        main:
        // create a channel for inputs from a CSV file
        name_ch = channel.fromPath(params.input)
                            .view { names -> "Before flatten: $names" }
                            // .flatten()
                            // .view { name -> "After flatten: $name" }
        // emit a greeting
        sayHello(name_ch)

        publish:
        greetings = sayHello.out
    }
    ```

=== "Before"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="4-10"
    workflow {

        main:
        // declare an array of input names
        names_array = ['Mateo','Andrea','Celine']
        // create a channel for inputs
        name_ch = channel.of(names_array)
                            .view { names -> "Before flatten: $names" }
                            .flatten()
                            .view { name -> "After flatten: $name" }
        // emit a greeting
        sayHello(name_ch)

        publish:
        greetings = sayHello.out
    }
    ```

You'll notice we switched the channel input to the new `params.input`, and deleted the `names_array` declaration since we'll no longer need it.
We've also commented out the `flatten()` and the second `view()` statement.

#### 4.1.3. Run the workflow

Let's try running the workflow with the new channel factory and the input file.

```bash
nextflow run hello-channels.nf
```

??? failure "Command output"

    ```console hl_lines="5 6 10 15"
    N E X T F L O W   ~  version 26.04.4

    Launching `hello-channels.nf` [exotic_plateau] revision: 8a8ed28d3c

    [-        ] sayHello | 0 of 1
    Before flatten: /workspaces/training/hello-nextflow/data/names.csv
    ERROR ~ Error executing process > 'sayHello (1)'

    Caused by:
      File `/workspaces/training/hello-nextflow/data/names.csv.txt` is outside the scope of the process work directory: /workspaces/training/hello-nextflow/work/df/d3e28a682956a3e686acf8d4f5ee9d


    Command executed:

      echo 'Hello, /workspaces/training/hello-nextflow/data/names.csv!' > '/workspaces/training/hello-nextflow/data/names.csv.txt'

    Command exit status:
      -

    Command output:
      (empty)

    Work dir:
      /workspaces/training/hello-nextflow/work/df/d3e28a682956a3e686acf8d4f5ee9d

    Tip: view the complete command output by changing to the process work dir and entering the command `cat .command.out`

    -- Check '.nextflow.log' file for details
    ```

Oh no, it doesn't work. Have a look at the start of the console output and error message.
The `Command executed:` bit is especially helpful here.

This may look a little bit familiar.
It looks like Nextflow tried to run a single process call using the file path itself as a string value.
So it has resolved the file path correctly, but it didn't actually parse its contents, which is what we wanted.

How do we get Nextflow to open the file and load its contents into the channel?

Sounds like we need another [operator](https://nextflow.io/docs/latest/reference/operator.html)!

### 4.2. Use the `splitCsv()` operator to parse the file

Looking through the list of operators again, we find [`splitCsv()`](https://nextflow.io/docs/latest/reference/operator.html#splitcsv), which is designed to parse and split CSV-formatted text.

#### 4.2.1. Apply `splitCsv()` to the channel

To apply the operator, we append it to the channel factory line like previously.

In the workflow block, make the following code change to replace `flatten()` with `splitCsv()` (uncommented):

=== "After"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="6-8"
    workflow {

        main:
        // create a channel for inputs from a CSV file
        name_ch = channel.fromPath(params.input)
                            .view { csv -> "Before splitCsv: $csv" }
                            .splitCsv()
                            .view { row -> "After splitCsv: $row" }
        // emit a greeting
        sayHello(name_ch)

        publish:
        greetings = sayHello.out
    }
    ```

=== "Before"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="6-8"
    workflow {

        main:
        // create a channel for inputs from a CSV file
        name_ch = channel.fromPath(params.input)
                            .view { names -> "Before flatten: $names" }
                            // .flatten()
                            // .view { name -> "After flatten: $name" }
        // emit a greeting
        sayHello(name_ch)

        publish:
        greetings = sayHello.out
    }
    ```

As you can see, we also updated the before/after `view()` statements.
The closure variables are renamed to describe what the channel contains at each point: `csv` for the CSV file before `splitCsv()`, and `row` for each parsed row after it.

#### 4.2.2. Run the workflow again

Let's try running the workflow with the added CSV-parsing logic.

```bash
nextflow run hello-channels.nf
```

??? failure "Command output"

    ```console hl_lines="7-10 14 19"
    N E X T F L O W   ~  version 26.04.4

    Launching `hello-channels.nf` [elegant_brown] revision: fe77324844

    executor >  local (3)
    [46/9ce81d] sayHello (3) | 0 of 3 ✘
    Before splitCsv: /workspaces/training/hello-nextflow/data/names.csv
    After splitCsv: [Mateo, es, squirrel]
    After splitCsv: [Andrea, it, turtle]
    After splitCsv: [Celine, fr, moose]
    ERROR ~ Error executing process > 'sayHello (3)'

    Caused by:
      Missing output file(s) `[Celine, fr, moose].txt` expected by process `sayHello (3)`


    Command executed:

      echo 'Hello, [Celine, fr, moose]!' > '[Celine, fr, moose].txt'

    Command exit status:
      0

    Command output:
      (empty)

    Work dir:
      /workspaces/training/hello-nextflow/work/46/9ce81d06e91d7dbfae4eb960c11a91

    Tip: view the complete command output by changing to the process work dir and entering the command `cat .command.out`

    -- Check '.nextflow.log' file for details
    ```

Interestingly, this fails too, but with a different error.
This time Nextflow has parsed the contents of the file (yay!) but it has loaded each row as an array, and each array is an element in the channel.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-channels-split-fail.svg"
</figure>

We need to tell it to only take the first column in each row, which holds the name.
So how do we unpack this?

We've previously used `flatten()` to unpack the contents of a channel, but that wouldn't work here because flatten unpacks _everything_ (feel free to try it if you want to see for yourself).

Instead, we'll use another operator called `map()` that is really useful and pops up a lot in Nextflow pipelines.

### 4.3. Use the `map()` operator to extract the names

The [`map()`](https://nextflow.io/docs/latest/reference/operator.html#map) operator is a very handy little tool that allows us to do all kinds of mappings to the contents of a channel.

In this case, we're going to use it to extract that one element that we want from each row in our data file.
This is what the syntax looks like:

```groovy title="Syntax"
.map { row -> row[0] }
```

This means 'for each row in the channel, take the 0th (first) item it contains'.

So let's apply that to our CSV parsing.

#### 4.3.1. Apply `map()` to the channel

In the workflow block, make the following code change:

=== "After"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="9 10"
    workflow {

        main:
        // create a channel for inputs from a CSV file
        name_ch = channel.fromPath(params.input)
                            .view { csv -> "Before splitCsv: $csv" }
                            .splitCsv()
                            .view { row -> "After splitCsv: $row" }
                            .map { row -> row[0] }
                            .view { name -> "After map: $name" }
        // emit a greeting
        sayHello(name_ch)

        publish:
        greetings = sayHello.out
    }
    ```

=== "Before"

    ```groovy title="hello-channels.nf" linenums="27"
    workflow {

        main:
        // create a channel for inputs from a CSV file
        name_ch = channel.fromPath(params.input)
                            .view { csv -> "Before splitCsv: $csv" }
                            .splitCsv()
                            .view { row -> "After splitCsv: $row" }
        // emit a greeting
        sayHello(name_ch)

        publish:
        greetings = sayHello.out
    }
    ```

You see we added another `view()` call to confirm that the operator does what we expect.

#### 4.3.2. Run the workflow

Let's run this one more time:

```bash
nextflow run hello-channels.nf
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `hello-channels.nf` [nauseous_jones] revision: 289cdedebd

    executor >  local (3)
    [7a/74aa8d] sayHello (2) | 3 of 3 ✔
    Before splitCsv: /workspaces/training/hello-nextflow/data/names.csv
    After splitCsv: [Mateo, es, squirrel]
    After splitCsv: [Andrea, it, turtle]
    After splitCsv: [Celine, fr, moose]
    After map: Mateo
    After map: Andrea
    After map: Celine

    Outputs:

      /workspaces/training/hello-nextflow/results

      greetings:
        - Celine.txt
        - Andrea.txt
        - Mateo.txt
    ```

This time it should run without error.

Looking at the output of the `view()` statements, you see the following:

- A single `Before splitCsv:` statement: at that point the channel contains one item, the original file path.
- Three separate `After splitCsv:` statements: one for each row, each contained within an array that corresponds to that line in the file.
- Three separate `After map:` statements: one for each name, which are now individual elements in the channel.

_Note that the lines may appear in a different order in your output._

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-channels-split-and-map.svg"
</figure>

You can also look at the output files to verify that each name was correctly extracted and processed through the workflow.

We've achieved the same result as previously, but now we have a lot more flexibility to add more names to process by modifying an input file, without modifying any code.
You'll learn learn more sophisticated approaches for handling complex inputs in a later training.

### Takeaway

You know how to use the `.fromPath()` channel constructor and the operators `splitCsv()` and `map()` to read in a file of input values and handle them appropriately.

More generally, you have a basic understanding of how Nextflow uses **channels** to manage inputs to processes and **operators** to transform their contents.
You've also seen how channels handle parallel execution implicitly.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-channels-parallel.svg"
</figure>

### What's next?

Take a big break, you worked hard in this one!

When you're ready, move on to [**Part 3: Hello Workflow**](./03_hello_workflow.md) to learn how to add more steps and connect them together into a proper workflow.

---

## Bonus exercises

These optional exercises let you practice what you learned in this part.
Each one can be solved with what this part covered; open the hint if you get stuck.
The solutions are also in `solutions/2-hello-channels/bonus/`.
To run the solution to Exercise 1 directly, add `--input solutions/2-hello-channels/bonus/names-header.csv` to the command.

### Exercise 1: Read a CSV file with a header row

Many CSV files start with a header row that names the columns.
Create a copy of `data/names.csv` called `data/names-header.csv`, and add the header row `name,language,character` at the top.
Then adapt `hello-channels.nf` so that it reads this file and picks the name by its column name instead of its position.
The workflow should still greet the same three names.

??? tip "Hint"

    `splitCsv()` has a `header` option.
    With `header: true`, each row becomes a set of named values, so you can write `row.name` instead of `row[0]`.

??? solution

    ```csv title="data/names-header.csv"
    name,language,character
    Mateo,es,squirrel
    Andrea,it,turtle
    Celine,fr,moose
    ```

    Point the `input` parameter to the new file:

    ```groovy title="hello-channels.nf" linenums="23" hl_lines="2"
    params {
        input: Path = 'data/names-header.csv'
    }
    ```

    Then parse the header and pick the column by name:

    ```groovy title="hello-channels.nf" linenums="31" hl_lines="3 5"
        name_ch = channel.fromPath(params.input)
                            .view { csv -> "Before splitCsv: $csv" }
                            .splitCsv(header: true)
                            .view { row -> "After splitCsv: $row" }
                            .map { row -> row.name }
                            .view { name -> "After map: $name" }
    ```

    Run the workflow:

    ```bash
    nextflow run hello-channels.nf
    ```

    ```console title="Output (excerpt)"
    After splitCsv: [name:Mateo, language:es, character:squirrel]
    After splitCsv: [name:Andrea, language:it, character:turtle]
    After splitCsv: [name:Celine, language:fr, character:moose]
    After map: Mateo
    After map: Andrea
    After map: Celine
    ```

    Each row is now a set of named values, and the header row itself is not passed on as a row.

### Exercise 2: Greet only the French names

Starting from `hello-channels.nf` as you left it at the end of section 4, without the changes from Exercise 1, change the workflow so that it only greets the names whose language, in the second column of `data/names.csv`, is `fr`.
Only one greeting file should be produced.

??? tip "Hint"

    Use the [`filter()`](https://nextflow.io/docs/latest/reference/operator.html#filter) operator.
    Like `map()`, it takes a closure, and keeps only the elements for which the closure returns `true`.

??? solution

    Add a `filter()` step between `splitCsv()` and `map()`, while each row still holds all its columns:

    ```groovy title="hello-channels.nf" linenums="31" hl_lines="5"
        name_ch = channel.fromPath(params.input)
                            .view { csv -> "Before splitCsv: $csv" }
                            .splitCsv()
                            .view { row -> "After splitCsv: $row" }
                            .filter { row -> row[1] == 'fr' }
                            .map { row -> row[0] }
                            .view { name -> "After map: $name" }
    ```

    Run the workflow:

    ```bash
    nextflow run hello-channels.nf
    ```

    ```console title="Output (excerpt)"
    [42/6d5454] sayHello (1) | 1 of 1 ✔
    Before splitCsv: /workspaces/training/hello-nextflow/data/names.csv
    After splitCsv: [Mateo, es, squirrel]
    After splitCsv: [Andrea, it, turtle]
    After splitCsv: [Celine, fr, moose]
    After map: Celine
    ```

    All three rows are read, but only Celine's row passes the filter, so `sayHello` runs once and writes `Hello, Celine!`.

---

## Quiz

<quiz>
What is a channel in Nextflow?
- [ ] A file path specification
- [ ] A process definition
- [x] A queue-like structure for passing data between processes
- [ ] A configuration setting

Learn more: [1.1. Create an input channel](#11-create-an-input-channel)
</quiz>

<quiz>
What will this code output?

```groovy
channel.of('Mateo', 'Andrea', 'Celine')
    .view()
```

- [ ] `['Mateo', 'Andrea', 'Celine']` (a single list)
- [x] Each element on a separate line: `Mateo`, `Andrea`, `Celine`
- [ ] Nothing (channels don't print by default)
- [ ] An error (invalid syntax)

Learn more: [1.1. Create an input channel](#11-create-an-input-channel)
</quiz>

<quiz>
When a channel contains multiple values, how does Nextflow handle process execution?
- [ ] The process runs once with all values
- [x] The process runs once for each value in the channel
- [ ] The process runs only with the first value
- [ ] The process runs only with the last value

Learn more: [2. Modify the workflow to run on multiple input values](#2-modify-the-workflow-to-run-on-multiple-input-values)
</quiz>

<quiz>
What does the `flatten()` operator do?
- [ ] Combines multiple channels into one
- [ ] Sorts channel elements
- [x] Unpacks arrays into individual elements
- [ ] Removes duplicate elements

Learn more: [3.2.1. Add the `flatten()` operator](#321-add-the-flatten-operator)
</quiz>

<quiz>
What is the purpose of the `view()` operator?
- [ ] To filter channel contents
- [ ] To transform channel elements
- [x] To inspect and debug channel contents
- [ ] To save channel contents to a file

Learn more: [1.4. Use `view()` to inspect the channel contents](#14-use-view-to-inspect-the-channel-contents)
</quiz>

<quiz>
What does `splitCsv()` do?
- [ ] Creates a CSV file from channel contents
- [ ] Splits a string by commas
- [x] Parses a CSV file into arrays representing each row
- [ ] Merges multiple CSV files

Learn more: [4.2. Use the `splitCsv()` operator to parse the file](#42-use-the-splitcsv-operator-to-parse-the-file)
</quiz>

<quiz>
What is the purpose of the `map()` operator?
- [ ] To filter elements from a channel
- [ ] To combine multiple channels
- [x] To transform each element in a channel
- [ ] To count elements in a channel

Learn more: [4.3. Use the `map()` operator to extract the names](#43-use-the-map-operator-to-extract-the-names)
</quiz>

<quiz>
Why is it important to use dynamic output filenames when processing multiple inputs?
- [ ] To improve performance
- [ ] To reduce disk space
- [x] To prevent output files from overwriting each other
- [ ] To enable resume functionality

Learn more: [2.2. Ensure the output file names will be unique](#22-ensure-the-output-file-names-will-be-unique)
</quiz>
