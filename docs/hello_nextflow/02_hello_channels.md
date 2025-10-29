# Part 2: Hello Channels

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/lJ41WMMm44M?si=xCItHLiOQWqoqBB9&amp;list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } See [the whole playlist](https://www.youtube.com/playlist?list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik) on the Nextflow YouTube channel.

:green_book: The video transcript is available [here](./transcripts/02_hello_channels.md).
///

In Part 1 of this course (Hello World), we showed you how to provide a variable input to a process by providing the input in the process call directly: `sayHello(params.greeting)`.
That was a deliberately simplified approach.
In practice, that approach has major limitations; namely that it only works for very simple cases where we only want to run the process once, on a single value.
In most realistic workflow use cases, we want to process multiple values (experimental data for multiple samples, for example), so we need a more sophisticated way to handle inputs.

That is what Nextflow **channels** are for.
Channels are queues designed to handle inputs efficiently and shuttle them from one step to another in multi-step workflows, while providing built-in parallelism and many additional benefits.

In this part of the course, you will learn how to use a channel to handle multiple inputs from a variety of different sources.
You will also learn to use **operators** to transform channel contents as needed.

_For training on using channels to connect steps in a multi-step workflow, see Part 3 of this course._

---

## 0. Warmup: Run `hello-channels.nf`

We're going to use the workflow script `hello-channels.nf` as a starting point.
It is equivalent to the script produced by working through Part 1 of this training course.

Just to make sure everything is working, run the script once before making any changes:

```bash
nextflow run hello-channels.nf --greeting 'Hello Channels!'
```

```console title="Output"
 N E X T F L O W   ~  version 25.04.3

Launching `hello-channels.nf` [insane_lichterman] DSL2 - revision: c33d41f479

executor >  local (1)
[86/9efa08] sayHello | 1 of 1 ✔
```

As previously, you will find the output file named `output.txt` in the `results` directory (specified by the `publishDir` directive).

```console title="results/output.txt" linenums="1"
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

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="3 4"
    workflow {

        // create a channel for inputs
        greeting_ch = channel.of('Hello Channels!')

        // emit a greeting
        sayHello(params.greeting)
    }
    ```

=== "Before"

    ```groovy title="hello-channels.nf" linenums="27"
    workflow {

        // emit a greeting
        sayHello(params.greeting)
    }
    ```

This is not yet functional since we haven't yet switched the input to the process call.

### 1.2. Add the channel as input to the process call

Now we need to actually plug our newly created channel into the `sayHello()` process call, replacing the CLI parameter which we were providing directly before.

In the workflow block, make the following code change:

=== "After"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="7"
    workflow {

        // create a channel for inputs
        greeting_ch = channel.of('Hello Channels!')

        // emit a greeting
        sayHello(greeting_ch)
    }
    ```

=== "Before"

    ```groovy title="hello-channels.nf" linenums="27"
    workflow {

        // create a channel for inputs
        greeting_ch = channel.of('Hello Channels!')

        // emit a greeting
        sayHello(params.greeting)
    }
    ```

This tells Nextflow to run the `sayHello` process on the contents of the `greeting_ch` channel.

Now our workflow is properly functional; it is the explicit equivalent of writing `sayHello('Hello Channels!')`.

### 1.3. Run the workflow command again

Let's run it!

```bash
nextflow run hello-channels.nf
```

If you made both edits correctly, you should get another successful execution:

```console title="Output" linenums="1"
 N E X T F L O W   ~  version 25.04.3

Launching `hello-channels.nf` [nice_heisenberg] DSL2 - revision: 41b4aeb7e9

executor >  local (1)
[3b/f2b109] sayHello (1) | 1 of 1 ✔
```

You can check the results directory to satisfy yourself that the outcome is still the same as previously.

```console title="results/output.txt" linenums="1"
Hello Channels!
```

So far we're just progressively tweaking the code to increase the flexibility of our workflow while achieving the same end result.

!!! note

    This may seem like we're writing more code for no tangible benefit, but the value will become clear as soon as we start handling more inputs.

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

    ```groovy title="hello-channels.nf" linenums="29" hl_lines="2"
    // create a channel for inputs
    greeting_ch = channel.of('Hello','Bonjour','Holà')
    ```

=== "Before"

    ```groovy title="hello-channels.nf" linenums="29"
    // create a channel for inputs
    greeting_ch = channel.of('Hello Channels')
    ```

The documentation tells us this should work. Can it really be so simple?

#### 2.1.2. Run the command and look at the log output

Let's try it.

```bash
nextflow run hello-channels.nf
```

It certainly seems to run just fine:

```console title="Output" linenums="1"
 N E X T F L O W   ~  version 25.04.3

Launching `hello-channels.nf` [suspicious_lamport] DSL2 - revision: 778deadaea

executor >  local (3)
[cd/77a81f] sayHello (3) | 3 of 3 ✔
```

However... This seems to indicate that '3 of 3' calls were made for the process, which is encouraging, but this only shows us a single run of the process, with one subdirectory path (`cd/77a81f`).
What's going on?

By default, the ANSI logging system writes the logging from multiple calls to the same process on the same line.
Fortunately, we can disable that behavior to see the full list of process calls.

#### 2.1.3. Run the command again with the `-ansi-log false` option

To expand the logging to display one line per process call, add `-ansi-log false` to the command.

```bash
nextflow run hello-channels.nf -ansi-log false
```

This time we see all three process runs and their associated work subdirectories listed in the output:

```console title="Output" linenums="1"
N E X T F L O W  ~  version 25.04.3
Launching `hello-channels.nf` [pensive_poitras] DSL2 - revision: 778deadaea
[76/f61695] Submitted process > sayHello (1)
[6e/d12e35] Submitted process > sayHello (3)
[c1/097679] Submitted process > sayHello (2)
```

That's much better; at least for a simple workflow.
For a complex workflow, or a large number of inputs, having the full list output to the terminal might get a bit overwhelming, so you might not choose to use `-ansi-log false` in those cases.

!!! note

    The way the status is reported is a bit different between the two logging modes.
    In the condensed mode, Nextflow reports whether calls were completed successfully or not.
    In this expanded mode, it only reports that they were submitted.

That being said, we have another problem. If you look in the `results` directory, there is only one file: `output.txt`!

```console title="Directory contents"
results
└── output.txt
```

What's up with that? Shouldn't we be expecting a separate file per input greeting, so three files in all?
Did all three greetings go into a single file?

You can check the contents of `output.txt`; you will find only one of the three, containing one of the three greetings we provided.

```console title="output.txt" linenums="1"
Bonjour
```

You may recall that we hardcoded the output file name for the `sayHello` process, so all three calls produced a file called `output.txt`.
You can check the work subdirectories for each of the three processes; each of them contains a file called `output.txt` as expected.

As long as the output files stay there, isolated from the other processes, that is okay.
But when the `publishDir` directive copies each of them to the same `results` directory, whichever got copied there first gets overwritten by the next one, and so on.

### 2.2. Ensure the output file names will be unique

We can continue publishing all the outputs to the same results directory, but we need to ensure they will have unique names.
Specifically, we need to modify the first process to generate a file name dynamically so that the final file names will be unique.

So how do we make the file names unique?
A common way to do that is to use some unique piece of metadata from the inputs (received from the input channel) as part of the output file name.
Here, for convenience, we'll just use the greeting itself since it's just a short string, and prepend it to the base output filename.

#### 2.2.1. Construct a dynamic output file name

In the process block, make the following code changes:

=== "After"

    ```groovy title="hello-channels.nf" linenums="6" hl_lines="9 13"
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

=== "Before"

    ```groovy title="hello-channels.nf" linenums="6"
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

Make sure to replace `output.txt` in both the output definition and in the `script:` command block.

!!! tip

    In the output definition, you MUST use double quotes around the output filename expression (NOT single quotes), otherwise it will fail.

This should produce a unique output file name every time the process is called, so that it can be distinguished from the outputs from other iterations of the same process in the output directory.

#### 2.2.2. Run the workflow

Let's run it:

```bash
nextflow run hello-channels.nf
```

Reverting back to the summary view, the output looks like this again:

```console title="Output" linenums="1"
 N E X T F L O W   ~  version 25.04.3

Launching `hello-channels.nf` [astonishing_bell] DSL2 - revision: f57ff44a69

executor >  local (3)
[2d/90a2e2] sayHello (1) | 3 of 3 ✔
```

Importantly, now we have three new files in addition to the one we already had in the `results` directory:

```console title="Directory contents"
results
├── Bonjour-output.txt
├── Hello-output.txt
├── Holà-output.txt
└── output.txt
```

They each have the expected contents:

```console title="Bonjour-output.txt" linenums="1"
Bonjour
```

```console title="Hello-output.txt" linenums="1"
Hello
```

```console title="Holà-output.txt" linenums="1"
Holà
```

Success! Now we can add as many greetings as we like without worrying about output files being overwritten.

!!! note

    In practice, naming files based on the input data itself is almost always impractical.
    The better way to generate dynamic filenames is to pass metadata to a process along with the input files.
    The metadata is typically provided via a 'sample sheet' or equivalents.
    You'll learn how to do that later in your Nextflow training.

### Takeaway

You know how to feed multiple input elements through a channel.

### What's next?

Learn to use an operator to transform the contents of a channel.

---

## 3. Use an operator to transform the contents of a channel

In Nextflow, [operators](https://www.nextflow.io/docs/latest/reference/operator.html) allow us to transform the contents of a channel.

We just showed you how to handle multiple input elements that were hardcoded directly in the channel factory.
What if we wanted to provide those multiple inputs in a different form?

For example, imagine we set up an input variable containing an array of elements like this:

`greetings_array = ['Hello','Bonjour','Holà']`

Can we load that into our output channel and expect it to work? Let's find out.

### 3.1. Provide an array of values as input to the channel

Common sense suggests we should be able to simply pass in an array of values instead of a single value. Right?

#### 3.1.1. Set up the input variable

Let's take the `greetings_array` variable we just imagined and make it a reality by adding it to the workflow block:

=== "After"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="3 4"
    workflow {

        // declare an array of input greetings
        greetings_array = ['Hello','Bonjour','Holà']

        // create a channel for inputs
        greeting_ch = channel.of('Hello','Bonjour','Holà')
    ```

=== "Before"

    ```groovy title="hello-channels.nf" linenums="27"
    workflow {

        // create a channel for inputs
        greeting_ch = channel.of('Hello','Bonjour','Holà')
    ```

#### 3.1.2. Set array of greetings as the input to the channel factory

We're going to replace the values `'Hello','Bonjour','Holà'` currently hardcoded in the channel factory with the `greetings_array` we just created.

In the workflow block, make the following change:

=== "After"

    ```groovy title="hello-channels.nf" linenums="32" hl_lines="2"
        // create a channel for inputs
        greeting_ch = channel.of(greetings_array)
    ```

=== "Before"

    ```groovy title="hello-channels.nf" linenums="32"
        // create a channel for inputs
        greeting_ch = channel.of('Hello','Bonjour','Holà')
    ```

#### 3.1.3. Run the workflow

Let's try running this:

```bash
nextflow run hello-channels.nf
```

Oh no! Nextflow throws an error that starts like this:

```console title="Output" linenums="1"
 N E X T F L O W   ~  version 25.04.3

Launching `hello-channels.nf` [friendly_koch] DSL2 - revision: 97256837a7

executor >  local (1)
[22/57e015] sayHello (1) | 0 of 1
ERROR ~ Error executing process > 'sayHello (1)'

Caused by:
  Missing output file(s) `[Hello, Bonjour, Holà]-output.txt` expected by process `sayHello (1)`
```

It looks like Nextflow tried to run a single process call, using `[Hello, Bonjour, Holà]` as a string value, instead of using the three strings in the array as separate values.

How do we get Nextflow to unpack the array and load the individual strings into the channel?

### 3.2. Use an operator to transform channel contents

This is where **operators** come in.

If you skim through the [list of operators](https://www.nextflow.io/docs/latest/reference/operator.html) in the Nextflow documentation, you'll find [`flatten()`](https://www.nextflow.io/docs/latest/reference/operator.html#flatten), which does exactly what we need: unpack the contents of an array and emits them as individual items.

!!! note

    It is technically possible to achieve the same results by using a different channel factory, [`channel.fromList`](https://nextflow.io/docs/latest/reference/channel.html#fromlist), which includes an implicit mapping step in its operation.
    Here we chose not to use that in order to demonstrate the use of an operator on a fairly simple use case.

#### 3.2.1. Add the `flatten()` operator

To apply the `flatten()` operator to our input channel, we append it to the channel factory declaration.

In the workflow block, make the following code change:

=== "After"

    ```groovy title="hello-channels.nf" linenums="31" hl_lines="3"
        // create a channel for inputs
        greeting_ch = channel.of(greetings_array)
                             .flatten()
    ```

=== "Before"

    ```groovy title="hello-channels.nf" linenums="31"
        // create a channel for inputs
        greeting_ch = channel.of(greetings_array)

    ```

Here we added the operator on the next line for readability, but you can add operators on the same line as the channel factory if you prefer, like this: `greeting_ch = channel.of(greetings_array).flatten()`

#### 3.2.2. Add `view()` to inspect channel contents

We could run this right away to test if it works, but while we're at it, we're also going to add a couple of [`view()`](https://www.nextflow.io/docs/latest/reference/operator.html#view) operators, which allow us to inspect the contents of a channel.
You can think of `view()` as a debugging tool, like a `print()` statement in Python, or its equivalent in other languages.

In the workflow block, make the following code change:

=== "After"

    ```groovy title="hello-channels.nf" linenums="31" hl_lines="3-5"
        // create a channel for inputs
        greeting_ch = channel.of(greetings_array)
                             .view { greeting -> "Before flatten: $greeting" }
                             .flatten()
                             .view { greeting -> "After flatten: $greeting" }
    ```

=== "Before"

    ```groovy title="hello-channels.nf" linenums="31"
        // create a channel for inputs
        greeting_ch = channel.of(greetings_array)
                             .flatten()
    ```

We are using an operator _closure_ here - the curly brackets.
This code executes for each item in the channel.
We define a temporary variable for the inner value, here called `greeting` (it could be anything).
This variable is only used within the scope of that closure.

In this example, `$greeting` represents each individual item loaded in a channel.

!!! note "Note on `$it`"

    In some pipelines you may see a special variable called `$it` used inside operator closures.
    This is an _implicit_ variable that allows a short-hand access to the inner variable,
    without needing to define it with a `->`.

    We prefer to be explicit to aid code clarity, as such the `$it` syntax is discouraged and will slowly be phased out of the Nextflow language.

#### 3.2.3. Run the workflow

Finally, you can try running the workflow again!

```bash
nextflow run hello-channels.nf
```

This time it works AND gives us the additional insight into what the contents of the channel look like before and after we run the `flatten()` operator:

```console title="Output" linenums="1"
 N E X T F L O W   ~  version 25.04.3

Launching `hello-channels.nf` [tiny_elion] DSL2 - revision: 1d834f23d2

executor >  local (3)
[8e/bb08f3] sayHello (2) | 3 of 3 ✔
Before flatten: [Hello, Bonjour, Holà]
After flatten: Hello
After flatten: Bonjour
After flatten: Holà
```

You see that we get a single `Before flatten:` statement because at that point the channel contains one item, the original array.
Then we get three separate `After flatten:` statements, one for each greeting, which are now individual items in the channel.

Importantly, this means each item can now be processed separately by the workflow.

!!! tip

    You should delete or comment out the `view()` statements before moving on.

    ```groovy title="hello-channels.nf" linenums="31"
    // create a channel for inputs
    greeting_ch = channel.of(greetings_array)
                         .flatten()
    ```

    We left them in the `hello-channels-3.nf` solution file for reference purposes.

### Takeaway

You know how to use an operator like `flatten()` to transform the contents of a channel, and how to use the `view()` operator to inspect channel contents before and after applying an operator.

### What's next?

Learn how to make the workflow take a file as its source of input values.

---

## 4. Use an operator to parse input values from a CSV file

It's often the case that, when we want to run on multiple inputs, the input values are contained in a file.
As an example, we prepared a CSV file called `greetings.csv` containing several greetings, one on each line (like a column of data).

```csv title="greetings.csv" linenums="1"
Hello
Bonjour
Holà
```

So now we need to modify our workflow to read in the values from a file like that.

### 4.1. Modify the script to expect a CSV file as the source of greetings

To get started, we're going to need to make two key changes to the script:

- Switch the input parameter to point to the CSV file
- Switch to a channel factory designed to handle a file

#### 4.1.1. Switch the input parameter to point to the CSV file

Remember the `params.greeting` parameter we set up in Part 1?
We're going to update it to point to the CSV file containing our greetings.

Before the workflow block, make the following code change:

=== "After"

    ```groovy title="hello-channels.nf" linenums="25" hl_lines="4"
    /*
     * Pipeline parameters
     */
    params.greeting = 'greetings.csv'
    ```

=== "Before"

    ```groovy title="hello-channels.nf" linenums="25"
    /*
     * Pipeline parameters
     */
    params.greeting = ['Hello','Bonjour','Holà']
    ```

#### 4.1.2. Switch to a channel factory designed to handle a file

Since we now want to use a file instead of simple strings as the input, we can't use the `channel.of()` channel factory from before.
We need to switch to using a new channel factory, [`channel.fromPath()`](https://www.nextflow.io/docs/latest/reference/channel.html#channel-path), which has some built-in functionality for handling file paths.

In the workflow block, make the following code change:

=== "After"

    ```groovy title="hello-channels.nf" linenums="31" hl_lines="1 2"
        // create a channel for inputs from a CSV file
        greeting_ch = channel.fromPath(params.greeting)
    ```

=== "Before"

    ```groovy title="hello-channels.nf" linenums="31"
        // create a channel for inputs
        greeting_ch = channel.of(greetings_array)
                             .flatten()
    ```

#### 4.1.3. Run the workflow

Let's try running the workflow with the new channel factory and the input file.

```bash
nextflow run hello-channels.nf
```

Oh no, it doesn't work. Here's the start of the console output and error message:

```console title="Output" linenums="1"
 N E X T F L O W   ~  version 25.04.3

Launching `hello-channels.nf` [adoring_bhabha] DSL2 - revision: 8ce25edc39

[-        ] sayHello | 0 of 1
ERROR ~ Error executing process > 'sayHello (1)'

Caused by:
  File `/workspaces/training/hello-nextflow/data/greetings.csv-output.txt` is outside the scope of the process work directory: /workspaces/training/hello-nextflow/work/e3/c459b3c8f4029094cc778c89a4393d


Command executed:

  echo '/workspaces/training/hello-nextflow/data/greetings.csv' > '/workspaces/training/hello-nextflow/data/greetings.
```

The `Command executed:` bit (lines 13-15) is especially helpful here.

This may look a little bit familiar.
It looks like Nextflow tried to run a single process call using the file path itself as a string value.
So it has resolved the file path correctly, but it didn't actually parse its contents, which is what we wanted.

How do we get Nextflow to open the file and load its contents into the channel?

Sounds like we need another [operator](https://www.nextflow.io/docs/latest/reference/operator.html)!

### 4.2. Use the `splitCsv()` operator to parse the file

Looking through the list of operators again, we find [`splitCsv()`](https://www.nextflow.io/docs/latest/reference/operator.html#splitCsv), which is designed to parse and split CSV-formatted text.

#### 4.2.1. Apply `splitCsv()` to the channel

To apply the operator, we append it to the channel factory line like previously.

In the workflow block, make the following code change:

=== "After"

    ```groovy title="hello-channels.nf" linenums="31" hl_lines="3-5"
    // create a channel for inputs from a CSV file
    greeting_ch = channel.fromPath(params.greeting)
                         .view { csv -> "Before splitCsv: $csv" }
                         .splitCsv()
                         .view { csv -> "After splitCsv: $csv" }
    ```

=== "Before"

    ```groovy title="hello-channels.nf" linenums="31"
    // create a channel for inputs from a CSV file
    greeting_ch = channel.fromPath(params.greeting)

    ```

As you can see, we also include before/after view statements while we're at it.

#### 4.2.2. Run the workflow again

Let's try running the workflow with the added CSV-parsing logic.

```bash
nextflow run hello-channels.nf
```

Interestingly, this fails too, but with a different error. The console output and error starts like this:

```console title="Output" linenums="1"
 N E X T F L O W   ~  version 25.04.3

Launching `hello-channels.nf` [stoic_ride] DSL2 - revision: a0e5de507e

executor >  local (3)
[42/8fea64] sayHello (1) | 0 of 3
Before splitCsv: /workspaces/training/hello-nextflow/greetings.csv
After splitCsv: [Hello]
After splitCsv: [Bonjour]
After splitCsv: [Holà]
ERROR ~ Error executing process > 'sayHello (2)'

Caused by:
  Missing output file(s) `[Bonjour]-output.txt` expected by process `sayHello (2)`


Command executed:

  echo '[Bonjour]' > '[Bonjour]-output.txt'
```

This time Nextflow has parsed the contents of the file (yay!) but it's added brackets around the greetings.

Long story short, `splitCsv()` reads each line into an array, and each comma-separated value in the line becomes an element in the array.
So here it gives us three arrays containing one element each.

!!! note

    Even if this behavior feels inconvenient right now, it's going to be extremely useful later when we deal with input files with multiple columns of data.

We could solve this by using `flatten()`, which you already know.
However, there's another operator called `map()` that's more appropriate to use here and is really useful to know; it pops up a lot in Nextflow pipelines.

### 4.3. Use the `map()` operator to extract the greetings

The `map()` operator is a very handy little tool that allows us to do all kinds of mappings to the contents of a channel.

In this case, we're going to use it to extract that one element that we want from each line of our file.
This is what the syntax looks like:

```groovy title="Syntax"
.map { item -> item[0] }
```

This means 'for each element in the channel, take the first of any items it contains'.

So let's apply that to our CSV parsing.

#### 4.3.1. Apply `map()` to the channel

In the workflow block, make the following code change:

=== "After"

    ```groovy title="hello-channels.nf" linenums="31" hl_lines="6-8"
    // create a channel for inputs from a CSV file
    greeting_ch = channel.fromPath(params.greeting)
                         .view { csv -> "Before splitCsv: $csv" }
                         .splitCsv()
                         .view { csv -> "After splitCsv: $csv" }
                         .map { item -> item[0] }
                         .view { csv -> "After map: $csv" }
    ```

=== "Before"

    ```groovy title="hello-channels.nf" linenums="31"
    // create a channel for inputs from a CSV file
    greeting_ch = channel.fromPath(params.greeting)
                         .view { csv -> "Before splitCsv: $csv" }
                         .splitCsv()
                         .view { csv -> "After splitCsv: $csv" }

    ```

Once again we include another `view()` call to confirm that the operator does what we expect.

#### 4.3.2. Run the workflow one more time

Let's run it one more time:

```bash
nextflow run hello-channels.nf
```

This time it should run without error.

```console title="Output" linenums="1"
 N E X T F L O W   ~  version 25.04.3

Launching `hello-channels.nf` [tiny_heisenberg] DSL2 - revision: 845b471427

executor >  local (3)
[1a/1d19ab] sayHello (2) | 3 of 3 ✔
Before splitCsv: /workspaces/training/hello-nextflow/greetings.csv
After splitCsv: [Hello]
After splitCsv: [Bonjour]
After splitCsv: [Holà]
After map: Hello
After map: Bonjour
After map: Holà
```

Looking at the output of the `view()` statements, we see the following:

- A single `Before splitCsv:` statement: at that point the channel contains one item, the original file path.
- Three separate `After splitCsv:` statements: one for each greeting, but each is contained within an array that corresponds to that line in the file.
- Three separate `After map:` statements: one for each greeting, which are now individual elements in the channel.

You can also look at the output files to verify that each greeting was correctly extracted and processed through the workflow.

We've achieved the same result as previously, but now we have a lot more flexibility to add more elements to the channel of greetings we want to process by modifying an input file, without modifying any code.

!!! note

    Here we had all greetings on one line in the CSV file.
    You can try adding more columns to the CSV file and see what happens; for example, try the following:

    ```csv title="greetings.csv"
    Hello,English
    Bonjour,French
    Holà,Spanish
    ```

    You can also try replacing `.map { item -> item[0] }` with `.flatten()` and see what happens depending on how many lines and columns you have in the input file.

    You'll learn learn more advanced approaches for handling complex inputs in a later training.

### Takeaway

You know how to use the operators `splitCsv()` and `map()` to read in a file of input values and handle them appropriately.

More generally, you have a basic understanding of how Nextflow uses **channels** to manage inputs to processes and **operators** to transform their contents.

### What's next?

Take a big break, you worked hard in this one!
When you're ready, move on to Part 3 to learn how to add more steps and connect them together into a proper workflow.
