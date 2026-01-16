# Part 2: Hello Channels

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/lJ41WMMm44M?si=xCItHLiOQWqoqBB9&amp;list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } See [the whole playlist](https://www.youtube.com/playlist?list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik) on the Nextflow YouTube channel.

:green_book: The video transcript is available [here](./transcripts/02_hello_channels.md).
///

In Part 1 of this course (Hello World), we showed you how to provide a variable input to a process by providing the input in the process call directly: `sayHello(params.input)`.
That was a deliberately simplified approach.
In practice, that approach has major limitations; namely that it only works for very simple cases where we only want to run the process once, on a single value.
In most realistic workflow use cases, we want to process multiple values (experimental data for multiple samples, for example), so we need a more sophisticated way to handle inputs.

That is what Nextflow **channels** are for.
Channels are queues designed to handle inputs efficiently and shuttle them from one step to another in multi-step workflows, while providing built-in parallelism and many additional benefits.

<!-- TODO: simple diagram for channels -->

In this part of the course, you will learn how to use a channel to handle multiple inputs from a variety of different sources.
You will also learn to use **operators** to transform channel contents as needed.

_For training on using channels to connect steps in a multi-step workflow, see Part 3 of this course._

---

## 0. Warmup: Run `hello-channels.nf`

We're going to use the workflow script `hello-channels.nf` as a starting point.
It is equivalent to the script produced by working through Part 1 of this training course, except we've changed the output destination:

```groovy title="hello-channels.nf" linenums="37" hl_lines="3"
output {
    first_output {
        path 'hello_channels'
        mode 'copy'
    }
}
```

Just to make sure everything is working, run the script once before making any changes:

```bash
nextflow run hello-channels.nf --input 'Hello Channels!'
```

??? success "Command output"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [wise_jennings] DSL2 - revision: b24f4902d6

    executor >  local (1)
    [6f/824bc1] process > sayHello [100%] 1 of 1 ✔
    ```

As previously, you will find the output file named `output.txt` in the `results/hello_channels` directory (as specified in the `output` block of the workflow script, shown above).

??? abstract "Directory contents"

    ```console title="results/hello_channels" hl_lines="2-3"
    results
    ├── hello_channels
    │   └── output.txt
    ├── hello_world
    │   └── output.txt
    └── output.txt -> /workspaces/training/hello-nextflow/work/8c/79499c11beea6e9d43605141f2817f/output.txt
    ```

??? abstract "File contents"

    ```console title="results/hello_channels/output.txt"
    Hello Channels!
    ```

If that worked for you, you're ready to learn about channels.

---

## 1. Provide variable inputs via a channel explicitly

We are going to create a **channel** to pass the variable input to the `sayHello()` process instead of relying on the implicit handling, which has certain limitations.

### 1.1. Create an input channel

There are a variety of **channel factories** that we can use to set up a channel.
To keep things simple for now, we are going to use the most basic channel factory, called `channel.of`, which will create a channel containing a single value.
Functionally this will be similar to how we had it set up before, but instead of having Nextflow create a channel implicitly, we are doing this explicitly now.

This is the line of code we're going to use:

```console title="Syntax"
greeting_ch = channel.of('Hello Channels!')
```

This creates a channel called `greeting_ch` using the `channel.of()` channel factory, which sets up a simple queue channel, and loads the string `'Hello Channels!'` to use as the greeting value.

!!! note

    We are temporarily switching back to hardcoded strings instead of using a CLI parameter for the sake of readability. We'll go back to using CLI parameters once we've covered what's happening at the level of the channel.

In the workflow block, add the channel factory code:

=== "After"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="4 5"
    workflow {

        main:
        // create a channel for inputs
        greeting_ch = channel.of('Hello Channels!')
        // emit a greeting
        sayHello(params.input)

        publish:
        first_output = sayHello.out
    }
    ```

=== "Before"

    ```groovy title="hello-channels.nf" linenums="27"
    workflow {

        main:
        // emit a greeting
        sayHello(params.input)

        publish:
        first_output = sayHello.out
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
        greeting_ch = channel.of('Hello Channels!')
        // emit a greeting
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

=== "Before"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="7"
    workflow {

        main:
        // create a channel for inputs
        greeting_ch = channel.of('Hello Channels!')
        // emit a greeting
        sayHello(params.input)

        publish:
        first_output = sayHello.out
    }
    ```

This tells Nextflow to run the `sayHello` process on the contents of the `greeting_ch` channel.

Now our workflow is properly functional; it is the explicit equivalent of writing `sayHello('Hello Channels!')`.

### 1.3. Run the workflow

Let's run it!

```bash
nextflow run hello-channels.nf
```

??? success "Command output"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [fabulous_crick] DSL2 - revision: 23e20f76e8

    executor >  local (1)
    [c0/4f1872] process > sayHello (1) [100%] 1 of 1 ✔
    ```

If you made both edits correctly, you should get a successful execution.
You can check the results directory to satisfy yourself that the outcome is still the same as previously.

??? abstract "File contents"

    ```console title="results/hello_channels/output.txt"
    Hello Channels!
    ```

So we've increased the flexibility of our workflow while achieving the same end result.
This may seem like we're writing more code for no tangible benefit, but the value will become clear as soon as we start handling more inputs.

As a preview of that, let's look at one more thing before we move on: one small but convenient benefit of using an explicit channel to manage data input.

### 1.4. Use `view()` to inspect the channel contents

Nextflow channels are built in a way that allows us to operate on their contents using operators, which we'll cover in detail later in this chapter.

For now, we're just going to show you how to use a super simple operator called [`view()`](https://www.nextflow.io/docs/latest/reference/operator.html#view) to inspect the contents of a channel.
You can think of `view()` as a debugging tool, like a `print()` statement in Python, or its equivalent in other languages.

Add this tiny line to the workflow block:

=== "After"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="7"
    workflow {

        main:
        // create a channel for inputs
        greeting_ch = channel.of('Hello Channels!')
                            .view()
        // emit a greeting
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

=== "Before"

    ```groovy title="hello-channels.nf" linenums="27"
    workflow {

        main:
        // create a channel for inputs
        greeting_ch = channel.of('Hello Channels!')
        // emit a greeting
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

The exact amount of spaces doesn't matter as long as it's a multiple of 4; we're just aiming to align the start of the `.view()` statement to the `.of()` part of the channel construction.

Now run the workflow again:

```bash
nextflow run hello-channels.nf
```

??? success "Command output"

    ```console hl_lines="7"
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [scruffy_shaw] DSL2 - revision: 2ede41e14a

    executor >  local (1)
    [ef/f7e40a] sayHello (1) [100%] 1 of 1 ✔
    Hello Channels!
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

### 2.1. Load multiple greetings into the input channel

Conveniently, the `channel.of()` channel factory we've been using is quite happy to accept more than one value, so we don't need to modify that at all.

We just have to load more values into the channel.

#### 2.1.1. Add more greetings

Before the workflow block, make the following code change:

=== "After"

    ```groovy title="hello-channels.nf" linenums="30" hl_lines="2"
    // create a channel for inputs
    greeting_ch = channel.of('Hello','Bonjour','Holà')
                        .view()
    ```

=== "Before"

    ```groovy title="hello-channels.nf" linenums="30" hl_lines="2"
    // create a channel for inputs
    greeting_ch = channel.of('Hello Channels')
                        .view()
    ```

The documentation tells us this should work. Can it really be so simple?

#### 2.1.2. Run the command and look at the log output

Let's try it.

```bash
nextflow run hello-channels.nf
```

??? success "Command output"

    ```console hl_lines="6"
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [amazing_crick] DSL2 - revision: 59a9a5888a

    executor >  local (3)
    [f4/c9962c] process > sayHello (1) [100%] 3 of 3 ✔
    Hello
    Bonjour
    Holà
    ```

It certainly seems to have run just fine.
The execution monitor shows that `3 of 3` calls were made for the `sayHello` process, and we see the three greetings enumerated by the `view()` statement, one per line as promised.

However, there is still only one output in the results directory:

??? abstract "Directory contents"

    ```console title="results/hello_channels" hl_lines="3"
    results
    ├── hello_channels
    │   └── output.txt
    ├── hello_world
    │   └── output.txt
    └── output.txt -> /workspaces/training/hello-nextflow/work/8c/79499c11beea6e9d43605141f2817f/output.txt
    ```

??? abstract "File contents"

    ```console title="results/hello_channels/output.txt"
    Holà
    ```

You should see one of the three greetings in there, but the one you got might be different from what is shown here.
Can you think of why that might be?

Looking back at the execution monitor, it gave us only one subdirectory path (`f4/c9962c`).
Let's have a look in there.

??? abstract "Directory contents"

    ```console hl_lines="9"
    work/f4/c9962ce91ef87480babcb86b2b9042/
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

    ```console title="work/f4/c9962ce91ef87480babcb86b2b9042/output.txt"
    Hello
    ```

That's not even the same greeting we got in the results directory! What's going on?

At this point, we need to tell you that by default, the ANSI logging system writes the logging from multiple calls to the same process on the same line.
So the status from all three calls to the sayHello() process are landing in the same spot.

Fortunately, we can disable that behavior to see the full list of process calls.

#### 2.1.3. Run the command again with the `-ansi-log false` option

To expand the logging to display one line per process call, add `-ansi-log false` to the command.

```bash
nextflow run hello-channels.nf -ansi-log false
```

??? success "Command output"

    ```console
     N E X T F L O W  ~  version 25.10.2
    Launching `hello-channels.nf` [desperate_monod] DSL2 - revision: 59a9a5888a
    Hello
    Bonjour
    Holà
    [23/871c7e] Submitted process > sayHello (2)
    [7f/21e2c2] Submitted process > sayHello (1)
    [f4/ea10a6] Submitted process > sayHello (3)
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
    work/23/871c7ec3642a898ecd5e6090d21300/
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
    work/7f/21e2c2f3cc8833ef3858b236e5575c/
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
    work/f4/ea10a680d5687596d3eaa3fcf69272/
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

    ```txt title="work/23/871c7ec3642a898ecd5e6090d21300/output.txt"
    Bonjour
    ```

    ```txt title="work/7f/21e2c2f3cc8833ef3858b236e5575c/output.txt"
    Hello
    ```

    ```txt title="work/f4/ea10a680d5687596d3eaa3fcf69272/output.txt"
    Holà
    ```

This shows that all three processes ran successfully (yay).

That being said, we still have the problem that there is only one output file in the results directory.

You may recall that we hardcoded the output file name for the `sayHello` process, so all three calls produced a file called `output.txt`.

As long as the output files stay in the work subdirectories, isolated from the other processes, that is okay.
But when they are published to the same results directory, whichever got copied there first gets overwritten by the next one, and so on.

### 2.2. Ensure the output file names will be unique

We can continue publishing all the outputs to the same results directory, but we need to ensure they will have unique names.
Specifically, we need to modify the first process to generate a file name dynamically so that the final file names will be unique.

So how do we make the file names unique?
A common way to do that is to use some unique piece of metadata from the inputs (received from the input channel) as part of the output file name.
Here, for convenience, we'll just use the greeting itself since it's just a short string, and prepend it to the base output filename.

#### 2.2.1. Construct a dynamic output file name

In the process block, make the following code changes:

=== "After"

    ```groovy title="hello-channels.nf" linenums="6" hl_lines="7 11"
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

=== "Before"

    ```groovy title="hello-channels.nf" linenums="6" hl_lines="7 11"
    process sayHello {

        input:
        val greeting

        output:
        path 'output.txt'

        script:
        """
        echo '${greeting}' > output.txt
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
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [sharp_minsky] DSL2 - revision: 16a291febe

    executor >  local (3)
    [e8/33ee64] sayHello (2) [100%] 3 of 3 ✔
    Hello
    Bonjour
    Holà
    ```

Reverting back to the summary view, the output is summarized on one line again.
Have a look at the `results` directory to see if all the output greetings are there.

??? abstract "Directory contents"

    ```console
    results/hello_channels/
    ├── Bonjour-output.txt
    ├── Hello-output.txt
    ├── Holà-output.txt
    └── output.txt
    ```

Yes! And they each have the expected contents.

??? abstract "File contents"

    ```console title="Bonjour-output.txt"
    Bonjour
    ```

    ```console title="Hello-output.txt"
    Hello
    ```

    ```console title="Holà-output.txt"
    Holà
    ```

Success! Now we can add as many greetings as we like without worrying about output files being overwritten.

!!! tip

    In practice, naming files based on the input data itself is almost always impractical.
    The better way to generate dynamic filenames is to pass metadata to a process along with the input files.
    The metadata is typically provided via a 'sample sheet' or equivalents.
    You'll learn how to do that later in your Nextflow training (see [Metadata side quest](../side_quests/metadata.md)).

### Takeaway

You know how to feed multiple input elements through a channel.

### What's next?

Learn to use an operator to transform the contents of a channel.

---

## 3. Use an operator to transform the contents of a channel

We just showed you how to handle multiple input elements that were hardcoded directly in the channel factory.
What if we wanted to provide those multiple inputs in a different way?

For example, imagine we set up an input variable containing an array of elements like this:

`greetings_array = ['Hello','Bonjour','Holà']`

Can we load that into our output channel and expect it to work? Let's find out.

### 3.1. Provide an array of values as input to the channel

Common sense suggests we should be able to simply pass in an array of values instead of a single value.
Let's try it; we'll need to set up the input variable and load it into the channel factory.

#### 3.1.1. Set up the input variable

Let's take the `greetings_array` variable we just imagined and make it a reality by adding it to the workflow block:

=== "After"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="4 5"
    workflow {

        main:
        // declare an array of input greetings
        greetings_array = ['Hello','Bonjour','Holà']
        // create a channel for inputs
        greeting_ch = channel.of('Hello','Bonjour','Holà')
                            .view()
        // emit a greeting
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

=== "Before"

    ```groovy title="hello-channels.nf" linenums="27"
    workflow {

        main:
        // create a channel for inputs
        greeting_ch = channel.of('Hello','Bonjour','Holà')
                            .view()
        // emit a greeting
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

This is not yet functional, we've just added a declaration for the array.

#### 3.1.2. Set array of greetings as the input to the channel factory

Now we're going to replace the values `'Hello','Bonjour','Holà'` currently hardcoded in the channel factory with the `greetings_array` we just created.

In the workflow block, make the following change:

=== "After"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="7"
    workflow {

        main:
        // declare an array of input greetings
        greetings_array = ['Hello','Bonjour','Holà']
        // create a channel for inputs
        greeting_ch = channel.of(greetings_array)
                            .view()
        // emit a greeting
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

=== "Before"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="7"
    workflow {

        main:
        // declare an array of input greetings
        greetings_array = ['Hello','Bonjour','Holà']
        // create a channel for inputs
        greeting_ch = channel.of('Hello','Bonjour','Holà')
                            .view()
        // emit a greeting
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
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
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [friendly_koch] DSL2 - revision: 97256837a7

    executor >  local (1)
    [a8/1f6ead] sayHello (1) | 0 of 1
    [Hello, Bonjour, Holà]
    ERROR ~ Error executing process > 'sayHello (1)'

    Caused by:
      Missing output file(s) `[Hello, Bonjour, Holà]-output.txt` expected by process `sayHello (1)`


    Command executed:

      echo '[Hello, Bonjour, Holà]' > '[Hello, Bonjour, Holà]-output.txt'

    Command exit status:
      0

    Command output:
      (empty)

    Work dir:
      /workspaces/training/hello-nextflow/work/a8/1f6ead5f3fa30a3c508e2e7cf83ffb

    Tip: you can replicate the issue by changing to the process work dir and entering the command `bash .command.run`

    -- Check '.nextflow.log' file for details
    ```

Oh no! There's an error!

Look at the output of `view()` and the error messages.

It looks like Nextflow tried to run a single process call, using `[Hello, Bonjour, Holà]` as a string value, instead of using the three strings in the array as separate values.

So it's the 'packaging' that is causing the problem.
How do we get Nextflow to unpack the array and load the individual strings into the channel?

### 3.2. Use an operator to transform channel contents

This is where **[operators](https://www.nextflow.io/docs/latest/reference/operator.html)** come into play.
You've already used the `.view()` operator, which just looks at what's in there.
Now we're going to look at operators that allow us to act on the contents of a channel.

If you skim through the [list of operators](https://www.nextflow.io/docs/latest/reference/operator.html) in the Nextflow documentation, you'll find [`flatten()`](https://www.nextflow.io/docs/latest/reference/operator.html#flatten), which does exactly what we need: unpack the contents of an array and emit them as individual items.

#### 3.2.1. Add the `flatten()` operator

To apply the `flatten()` operator to our input channel, we append it to the channel factory declaration.

In the workflow block, make the following code change:

=== "After"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="9"
    workflow {

        main:
        // declare an array of input greetings
        greetings_array = ['Hello','Bonjour','Holà']
        // create a channel for inputs
        greeting_ch = channel.of(greetings_array)
                            .view()
                            .flatten()
        // emit a greeting
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

=== "Before"

    ```groovy title="hello-channels.nf" linenums="27"
    workflow {

        main:
        // declare an array of input greetings
        greetings_array = ['Hello','Bonjour','Holà']
        // create a channel for inputs
        greeting_ch = channel.of(greetings_array)
                            .view()
        // emit a greeting
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

Here we added the operator on the next line for readability, but you can add operators on the same line as the channel factory if you prefer, like this:
`greeting_ch = channel.of(greetings_array).view().flatten()`

#### 3.2.2. Refine the `view()` statement(s)

We could run this right away to test if it works, but while we're at it, we're going to refine how we inspect the channel contents.

We want to be able to contrast what the contents look like before and after the `flatten()` operator is applied, so we're going to add a second one, AND we're going to add a bit of code to get them labeled more clearly in the output.

In the workflow block, make the following code change:

=== "After"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="8-10"
    workflow {

        main:
        // declare an array of input greetings
        greetings_array = ['Hello','Bonjour','Holà']
        // create a channel for inputs
        greeting_ch = channel.of(greetings_array)
                             .view { greeting -> "Before flatten: $greeting" }
                             .flatten()
                             .view { greeting -> "After flatten: $greeting" }
        // emit a greeting
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

=== "Before"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="8-9"
    workflow {

        main:
        // declare an array of input greetings
        greetings_array = ['Hello','Bonjour','Holà']
        // create a channel for inputs
        greeting_ch = channel.of(greetings_array)
                            .view()
                            .flatten()
        // emit a greeting
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

You see we've added a second `.view` statement, and for each of them, we've replaced the empty parentheses (`()`) with curly braces containing some code, such as `{ greeting -> "Before flatten: $greeting" }`.

These are called _closures_. The code they contain will be executed for each item in the channel.
We define a temporary variable for the inner value, here called `greeting` (but it could be any arbitrary name), which is only used within the scope of that closure.

In this example, `$greeting` represents each individual item loaded in the channel.
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
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [sleepy_gutenberg] DSL2 - revision: 1db4f760ee

    executor >  local (3)
    [b1/6a1e15] sayHello (2) [100%] 3 of 3 ✔
    Before flatten: [Hello, Bonjour, Holà]
    After flatten: Hello
    After flatten: Bonjour
    After flatten: Holà
    ```

This time it works AND gives us the additional insight into what the contents of the channel look like before and after we run the `flatten()` operator.

- - You see that we get a single `Before flatten:` statement because at that point the channel contains one item, the original array.
    Then we get three separate `After flatten:` statements, one for each greeting, which are now individual items in the channel.

Importantly, this means each item can now be processed separately by the workflow.

!!! tip

    It is technically possible to achieve the same results by using a different channel factory, [`channel.fromList`](https://nextflow.io/docs/latest/reference/channel.html#fromlist), which includes an implicit mapping step in its operation.
    Here we chose not to use that in order to demonstrate the use of an operator on a simple use case.

### Takeaway

You know how to use an operator like `flatten()` to transform the contents of a channel, and how to use the `view()` operator to inspect channel contents before and after applying an operator.

### What's next?

Learn how to make the workflow take a file as its source of input values.

---

## 4. Use an operator to parse input values from a CSV file

Realistically, we're rarely if ever going to start from an array of values.
Most likely, we'll have one or more files containing the data that needs to be processed, in some kind of structured format.

We've prepared a a CSV file called `greetings.csv` that contains several input greetings, mimicking the kind of columnar data you might want to process in a real data analysis, stored under `data/`.

```csv title="data/greetings.csv" linenums="1"
Hello,English,123
Bonjour,French,456
Holà,Spanish,789
```

Note that the numbers are not meaningful, they are just there for illustrative purposes.

Our next task then is to adapt our workflow to read in the values from this file.

### 4.1. Modify the script to expect a CSV file as the source of greetings

To get started, we're going to need to make two key changes to the script:

- Switch the input parameter to point to the CSV file
- Switch the channel factory to one designed to handle a file

#### 4.1.1. Switch the input parameter to point to the CSV file

Remember the `params.input` parameter we set up in Part 1?
We're going to update it to point to the CSV file containing our greetings.

Make the following edit to the parameter declaration:

=== "After"

    ```groovy title="hello-channels.nf" linenums="20" hl_lines="5"
    /*
    * Pipeline parameters
    */
    params {
        input: Path = 'data/greetings.csv'
    }
    ```

=== "Before"

    ```groovy title="hello-channels.nf" linenums="20" hl_lines="5"
    /*
     * Pipeline parameters
     */
    input: String = 'Holà mundo!'
    ```

This assumes the file is co-located with the workflow code.
You'll learn how to deal with other data locations later in your Nextflow journey.

#### 4.1.2. Switch to a channel factory designed to handle a file

Since we now want to use a file instead of simple strings as the input, we can't use the `channel.of()` channel factory from before.
We need to switch to using a new channel factory, [`channel.fromPath()`](https://www.nextflow.io/docs/latest/reference/channel.html#channel-path), which has some built-in functionality for handling file paths.

In the workflow block, make the following code change:

=== "After"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="4-8"
    workflow {

        main:
        // create a channel for inputs from a CSV file
        greeting_ch = channel.fromPath(params.input)
                            .view { greeting -> "Before flatten: $greeting" }
                            //.flatten()
                            //.view { greeting -> "After flatten: $greeting" }
        // emit a greeting
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

=== "Before"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="4-8"
    workflow {

        main:
        // declare an array of input greetings
        greetings_array = ['Hello','Bonjour','Holà']
        // create a channel for inputs
        greeting_ch = channel.of(greetings_array)
                            .view { greeting -> "Before flatten: $greeting" }
                            .flatten()
                            .view { greeting -> "After flatten: $greeting" }
        // emit a greeting
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

You'll notice we switched the channel input back to `param.input`, and deleted the `greetings_array` declaration since we'll no longer need it.
We've also commented out the `flatten()` and the second `view()` statement.

#### 4.1.3. Run the workflow

Let's try running the workflow with the new channel factory and the input file.

```bash
nextflow run hello-channels.nf
```

??? failure "Command output"

    ```console hl_lines="5 6 9 14"
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [peaceful_poisson] DSL2 - revision: a286c08ad5

    [-        ] sayHello [  0%] 0 of 1
    Before flatten: /workspaces/training/hello-nextflow/data/greetings.csv
    ERROR ~ Error executing process > 'sayHello (1)'

    Caused by:
      File `/workspaces/training/hello-nextflow/data/greetings.csv-output.txt` is outside the scope of the process work directory: /workspaces/training/hello-nextflow/work/30/e610cb4ea5ae8693f456ac3329c92f


    Command executed:

      echo '/workspaces/training/hello-nextflow/data/greetings.csv' > '/workspaces/training/hello-nextflow/data/greetings.csv-output.txt'

    Command exit status:
      -

    Command output:
      (empty)

    Work dir:
      /workspaces/training/hello-nextflow/work/30/e610cb4ea5ae8693f456ac3329c92f

    Tip: when you have fixed the problem you can continue the execution adding the option `-resume` to the run command line

    -- Check '.nextflow.log' file for details
    ```

Oh no, it doesn't work. Have a look at the start of the console output and error message.
The `Command executed:` bit is especially helpful here.

This may look a little bit familiar.
It looks like Nextflow tried to run a single process call using the file path itself as a string value.
So it has resolved the file path correctly, but it didn't actually parse its contents, which is what we wanted.

How do we get Nextflow to open the file and load its contents into the channel?

Sounds like we need another [operator](https://www.nextflow.io/docs/latest/reference/operator.html)!

### 4.2. Use the `splitCsv()` operator to parse the file

Looking through the list of operators again, we find [`splitCsv()`](https://www.nextflow.io/docs/latest/reference/operator.html#splitCsv), which is designed to parse and split CSV-formatted text.

#### 4.2.1. Apply `splitCsv()` to the channel

To apply the operator, we append it to the channel factory line like previously.

In the workflow block, make the following code change to replace `flatten()` with `splitcsv()` (uncommented):

=== "After"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="6-8"
    workflow {

        main:
        // create a channel for inputs from a CSV file
        greeting_ch = channel.fromPath(params.input)
                            .view { csv -> "Before splitCsv: $csv" }
                            .splitCsv()
                            .view { csv -> "After splitCsv: $csv" }
        // emit a greeting
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

=== "Before"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="6-8"
    workflow {

        main:
        // create a channel for inputs from a CSV file
        greeting_ch = channel.fromPath(params.input)
                            .view { greeting -> "Before flatten: $greeting" }
                            //.flatten()
                            //.view { greeting -> "After flatten: $greeting" }
        // emit a greeting
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

As you can see, we also updated the before/after `view()` statements.
Technically we could have used the same variable name (`greeting`) but we updated it to something more appropriate (`csv`) to make the code more readable by others.

#### 4.2.2. Run the workflow again

Let's try running the workflow with the added CSV-parsing logic.

```bash
nextflow run hello-channels.nf
```

??? failure "Command output"

    ```console hl_lines="7-11 14 19"
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [insane_fermat] DSL2 - revision: 8e62fcbeb1

    executor >  local (3)
    [24/76da2f] sayHello (2) [  0%] 0 of 3 ✘
    Before splitCsv: /workspaces/training/hello-nextflow/data/greetings.csv
    After splitCsv: [Hello, English, 123]
    After splitCsv: [Bonjour, French, 456]
    After splitCsv: [Holà, Spanish, 789]
    ERROR ~ Error executing process > 'sayHello (2)'

    Caused by:
      Missing output file(s) `[Bonjour, French, 456]-output.txt` expected by process `sayHello (2)`


    Command executed:

      echo '[Bonjour, French, 456]' > '[Bonjour, French, 456]-output.txt'

    Command exit status:
      0

    Command output:
      (empty)

    Work dir:
      /workspaces/training/hello-nextflow/work/24/76da2fcc4876b61632749f99e26a50

    Tip: you can try to figure out what's wrong by changing to the process work dir and showing the script file named `.command.sh`

    -- Check '.nextflow.log' file for details
    ```

Interestingly, this fails too, but with a different error.
This time Nextflow has parsed the contents of the file (yay!) but it has loaded each row as an array, and each array is an element in the channel.

We need to tell it to only take the first column in each row.
So how do we unpack this?

We've previously used `flatten()` to unpack the contents of a channel, but that wouldn't work here because flatten unpacks _everything_ (feel free to try it if you want to see for yourself).

Instead, we'll use another operator called `map()` that is really useful and pops up a lot in Nextflow pipelines.

### 4.3. Use the `map()` operator to extract the greetings

The [`map()`](https://www.nextflow.io/docs/latest/reference/operator.html#map) operator is a very handy little tool that allows us to do all kinds of mappings to the contents of a channel.

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
        greeting_ch = channel.fromPath(params.input)
                            .view { csv -> "Before splitCsv: $csv" }
                            .splitCsv()
                            .view { csv -> "After splitCsv: $csv" }
                            .map { item -> item[0] }
                            .view { csv -> "After map: $csv" }
        // emit a greeting
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

=== "Before"

    ```groovy title="hello-channels.nf" linenums="27"
    workflow {

        main:
        // create a channel for inputs from a CSV file
        greeting_ch = channel.fromPath(params.input)
                            .view { csv -> "Before splitCsv: $csv" }
                            .splitCsv()
                            .view { csv -> "After splitCsv: $csv" }
        // emit a greeting
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
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
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [focused_volhard] DSL2 - revision: de435e45be

    executor >  local (3)
    [54/6eebe3] sayHello (3) [100%] 3 of 3 ✔
    Before splitCsv: /workspaces/training/hello-nextflow/data/greetings.csv
    After splitCsv: [Hello, English, 123]
    After splitCsv: [Bonjour, French, 456]
    After splitCsv: [Holà, Spanish, 789]
    After map: Hello
    After map: Bonjour
    After map: Holà
    ```

This time it should run without error.

Looking at the output of the `view()` statements, you see the following:

- A single `Before splitCsv:` statement: at that point the channel contains one item, the original file path.
- Three separate `After splitCsv:` statements: one for each greeting, but each is contained within an array that corresponds to that line in the file.
- Three separate `After map:` statements: one for each greeting, which are now individual elements in the channel.

Note that the lines may appear in a different order in your output.

You can also look at the output files to verify that each greeting was correctly extracted and processed through the workflow.

We've achieved the same result as previously, but now we have a lot more flexibility to add more elements to the channel of greetings we want to process by modifying an input file, without modifying any code.
You'll learn learn more sophisticated approaches for handling complex inputs in a later training.

### Takeaway

You know how to use the `.fromPath()` channel constructor and the operators `splitCsv()` and `map()` to read in a file of input values and handle them appropriately.

More generally, you have a basic understanding of how Nextflow uses **channels** to manage inputs to processes and **operators** to transform their contents.

### What's next?

Take a big break, you worked hard in this one!

When you're ready, move on to [**Part 3: Hello Workflow**](./03_hello_workflow.md) to learn how to add more steps and connect them together into a proper workflow.
