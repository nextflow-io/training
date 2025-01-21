# Part 2: Hello Channels

In Part 1 of this course (Hello World), we showed you how to provide a variable input to a process by providing the input in the process call directly: `sayHello(params.greet)`.
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
nextflow run hello-channels.nf --greet 'Hello Channels!'
```

```console title="Output"
 N E X T F L O W   ~  version 24.10.0

Launching `hello-channels.nf` [silly_moriondo] DSL2 - revision: a0dfbc86fe

executor >  local (1)
[03/b321a3] sayHello [100%] 1 of 1 ✔
```

---

## 1. Provide variable inputs via a channel explicitly

We are going to create a **channel** to pass the variable input to the `sayHello()` process instead of relying on the implicit handling, which has certain limitations.

### 1.1. Create an input channel

There are a variety of **channel factories** that we can use to set up a channel.
To keep things simple for now, we are going to use the simplest possible channel factory, which will create a simple channel containing a single value.
Functionally this be exactly equivalent to how we had it set up before.

This is the line of code to do it:

`greeting_ch = Channel.of('Hello world!')`

This creates a channel called `greeting_ch` using the `Channel.of()` channel factory, which sets up a simple value channel, and gives it the string `'Hello world!'` to use as the greeting value.

!!! note

    We are temporarily switching back to hardcoded strings instead of using a CLI parameter for the sake of simplicity. We'll go back to using CLI parameters once we've covered what's happening at the level of the channel.

In the workflow block, add the channel factory code:

_Before:_

```groovy title="hello-channels.nf" linenums="21"
workflow {

    // emit a greeting
    sayHello(params.greet)
}
```

_After:_

```groovy title="hello-channels.nf" linenums="21"
workflow {

    // create a channel for inputs
    greeting_ch = Channel.of('Hello world!')

    // emit a greeting
    sayHello(params.greet)
}
```

This is not yet functional since we haven't yet switched the input to the process call.

### 1.2. Add the channel as input to the process call

Now we need to actually plug our newly created channel into the `sayHello()` process call, replacing the CLI parameter which we were providing directly before.

In the workflow block, make the following code change:

_Before:_

```groovy title="hello-channels.nf" linenums="26"
workflow {

    // create a channel for inputs
    greeting_ch = Channel.of('Hello world!')

    // emit a greeting
    sayHello(params.greet)
}
```

_After:_

```groovy title="hello-channels.nf" linenums="26"
workflow {

    // create a channel for inputs
    greeting_ch = Channel.of('Hello world!')

    // emit a greeting
    sayHello(greeting_ch)
}
```

This tells Nextflow to run the `sayHello` process on the contents of the `greeting_ch` channel.

Now it's fully functional; it's the explicit equivalent of writing `sayHello('Hello world!')`.

### 1.3. Run the workflow command again

Let's run it!

```bash
nextflow run hello-channels.nf
```

If you made both edits correctly, you should get another successful execution:

```console title="Output"
 N E X T F L O W   ~  version 24.10.0

 ┃ Launching `hello-channels.nf` [prickly_avogadro] DSL2 - revision: b58b6ab94b

executor >  local (1)
[1f/50efd5] sayHello (1) [100%] 1 of 1 ✔
```

Feel free to check the results directory to satisfy yourself that the outcome is still the same as previously.
So far we're just progressively tweaking the code to increase the flexibility of our workflow while achieving the same end result.

!!! note

    This may seem like we're writing more code for no tangible benefit, but the value will become clear as soon as we start handling more complex inputs.

### Takeaway

You know how to use a simple channel to provide an input to a process.

### What's next?

Learn how to use channels to make the workflow iterate over multiple input values.

---

## 2. Modify the workflow to run on multiple input values

Workflows typically run on batches of inputs that are meant to be processed in bulk, so we want to upgrade the workflow to accept multiple input values.

### 2.1. Load multiple greetings into the input channel

Conveniently, the `Channel.of()` channel factory we've been using is quite happy to accept more than one value, so we don't need to modify that at all.
We just have to load more values into the channel.

#### 2.1.1. Add more greetings

In the workflow block, make the following code change:

_Before:_

```groovy title="hello-channels.nf" linenums="46"
// create a channel for inputs
greeting_ch = Channel.of('Hello')
```

_After:_

```groovy title="hello-channels.nf" linenums="46"
// create a channel for inputs
greeting_ch = Channel.of('Hello','Bonjour','Holà')
```

The documentation tells us this should work. Can it really be so simple?

#### 2.1.2. Run the command and look at the log output

Let's try it.

```bash
nextflow run hello-channels.nf
```

Well, it certainly seems to run just fine:

```console title="Output"
 N E X T F L O W   ~  version 24.10.0

 ┃ Launching `hello-channels.nf` [lonely_pare] DSL2 - revision: b9f1d96905

executor >  local (3)
[3d/1fe62c] sayHello (2)       [100%] 3 of 3 ✔
```

However... This seems to indicate that '3 of 3' calls were made for the process, which is encouraging, but this only shows us a single run of the process, with one subdirectory path (`3d/1fe62c...`). What's going on?

By default, the ANSI logging system writes the logging from multiple calls to the same process on the same line. Fortunately, we can disable that behavior.

#### 2.1.3. Run the command again with the `-ansi-log false` option

To expand the logging to display one line per process call, add `-ansi-log false` to the command.

```bash
nextflow run hello-channels.nf -ansi-log false
```

This time we see all three process runs and their associated work subdirectories listed in the output:

```console title="Output"
N E X T F L O W  ~  version 24.02.0-edge
Launching `hello-channels.nf` [big_woese] DSL2 - revision: 53f20aeb70
[62/d81e63] Submitted process > sayHello (1)
[19/507af3] Submitted process > sayHello (2)
[8a/3126e6] Submitted process > sayHello (3)
```

That's much better; at least for a simple workflow.
For a complex workflow, or a large number of inputs, having the full list output to the terminal might get a bit overwhelming, so you might not choose to use `-ansi-log false` in those cases.

That being said, we have another problem. If you look in the `results` directory, there is only one file: `output.txt`!

```console title="Directory contents"
results
└── output.txt
```

What's up with that? Shouldn't we be expecting a separate file per input greeting, so three files in all?
Did all three greetings go into a single file?
You can check the contents of `output.txt`; you will find only one of the three.

You may recall that we hardcoded the output file name for the `sayHello` process, so all three calls produced a file called `output.txt`.
You can check the work subdirectories for each of the three processes; each of them contains a file called `output.txt` as expected.

As long as the output files stay there, isolated from the other processes, that is okay.
But when the `publishDir` directive copies each of them to the same `results` directory, whichever got copied there first gets overwritten by the next one, and so on.

### 2.2. Ensure the output file names will be unique

We can continue publishing all the outputs to the same results directory, but we need to ensure they will have unique names.
Specifically, we need to modify the first process to generate a file name dynamically so that the final file names will be unique.

So how do we make the file names unique? A common way to do that is to use some unique piece of metadata from the inputs (received from the input channel) as part of the output file name.
Here, for convenience, we'll just use the greeting itself since it's just a short string, and prepend it to the base output filename.

#### 2.2.1. Construct a dynamic output file name

In the process block, make the following code changes:

_Before:_

```groovy title="hello-channels.nf" linenums="11"
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

```groovy title="hello-channels.nf" linenums="11"
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

Make sure to replace `output.txt` in both the output definition and in the `script:` command block.

!!! tip

    In the output definition, you MUST use double quotes around the output filename expression (NOT single quotes), otherwise it will fail.

This should produce a unique output file name every time the process is called, so that it can be distinguished from the outputs from other iterations of the same process in the output directory.

#### 2.2.2. Run the workflow and look at the results directory

Let's run it and check that it works.

```bash
nextflow run hello-channels.nf
```

Reverting back to the summary view, the output looks like this again:

```console title="Output"
 N E X T F L O W   ~  version 24.10.0

 ┃ Launching `hello-channels.nf` [jovial_mccarthy] DSL2 - revision: 53f20aeb70

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
    The better way to generate dynamic filenames is to pass metatdata to a process along with the input files. We can derive that metadata from a sample sheet as we're reading the input files.
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

`greetings = ['Hello','Bonjour','Holà']`

Can we load that into our output channel and expect it to work? Let's find out.

### 3.1. Provide an array of values as input to the channel

Common sense suggests we should be able to simply pass in an array of values instead of a single value. Right?

#### 3.1.1. Set up the input variable

Since we already have the `params.greet` declared, let's hijack that by changing its value to an array as follows:

_Before:_

```groovy title="hello-channels.nf" linenums="23"
/*
 * Pipeline parameters
 */
params.greet = 'Holà mundo'
```

_After:_

```groovy title="hello-channels.nf" linenums="23"
/*
 * Pipeline parameters
 */
params.greet = ['Hello','Bonjour','Holà']
```

#### 3.1.2. Set `params.greet` as the input to the channel factory

We replace the hardcoded values `'Hello','Bonjour','Holà'` with the `params.greet` that we just updated to contain the array `['Hello','Bonjour','Holà']`.

Modify the following code:

_Before:_

```groovy title="hello-channels.nf" linenums="23"
// create a channel for inputs
greeting_ch = Channel.of('Hello','Bonjour','Holà')
```

_After:_

```groovy title="hello-channels.nf" linenums="23"
// create a channel for inputs
greeting_ch = Channel.of(params.greet)
```

#### 3.1.3. Run the workflow

Let's try running this:

```bash
nextflow run hello-channels.nf
```

Oh no! Nextflow throws an error that starts like this:

```console title="Output"
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

    It is technically possible to achieve the same results by using a different channel factory, [`Channel.fromList`](https://nextflow.io/docs/latest/reference/channel.html#fromlist), which includes an implicit mapping step in its operation.
    Here we chose not to use that in order to demonstrate the use of an operator on a fairly simple use case.

#### 3.2.1. Add the `flatten()` operator

To apply the `flatten()` operator to our input channel, we simply append it to the channel factory declaration.

In the workflow block, make the following code change:

_Before:_

```groovy title="hello-channels.nf" linenums="29"
// create a channel for inputs
greeting_ch = Channel.of(params.greeting)
```

_After:_

```groovy title="hello-channels.nf" linenums="29"
// create a channel for inputs
greeting_ch = Channel.of(params.greeting)
                     .flatten()
```

Here we added the operator on the next line for readability, but you can add operators on the same line as the channel factory if you prefer.

#### 3.2.2. Add `view()` to inspect channel contents

We could run this right away to test if it works, but while we're at it, we're also going to add a couple of [`view()`](https://www.nextflow.io/docs/latest/reference/operator.html#view) directives, which allow us to inspect the contents of a channel.
You can think of `view()` as a debugging tool, like a `print()` statement in Python, if you're familiar with that.

In the workflow block, make the following code change:

_Before:_

```groovy title="hello-channels.nf" linenums="29"
// create a channel for inputs
greeting_ch = Channel.of(params.greeting)
                     .flatten()
```

_After:_

```groovy title="hello-channels.nf" linenums="29"
// create a channel for inputs
greeting_ch = Channel.of(params.greeting)
                     .view { "Before flatten: $it" }
                     .flatten()
                     .view { "After flatten: $it" }
```

Here `$it` is an implicit variable that represents each individual item loaded in a channel.

#### 3.2.3. Run the workflow

Finally, you can try running the workflow again!

```bash
nextflow run hello-channels.nf
```

This time it works AND gives us the additional insight into what the contents of the channel look like before and after we run the `flatten()` operator:

TODO UPDATE

```console title="Output"
Launching `hello-channels.nf` [irreverent_shaw] DSL2 - revision: b3a71cc376

executor >  local (3)
[20/8f32e4] sayHello (1) [100%] 3 of 3 ✔
Before flatten: [Holà mundo, Konnichiwa, Dobrý den]
After flatten: Holà mundo
After flatten: Konnichiwa
After flatten: Dobrý den
```

You see that we get a single `Before flatten:` statement because at that point the channel contains one item, the original array.
Then we get three separate `After flatten:` statements, one for each greeting, which are now individual items in the channel.

Importantly, this means each item can now be processed separately by the workflow.

!!! tip

    You should delete or comment out the `view()` statements before moving on.

    ```groovy title="hello-channels.nf" linenums="29"
    // create a channel for inputs
    greeting_ch = Channel.of(params.greeting)
                         .flatten()
    ```

### Takeaway

You know how to use an operator like `flatten()` to transform the contents of a channel, and how to use the `view()` directive to inspect channel contents before and after applying an operator.

### What's next?

Learn how to make the workflow take a file as its source of input values.

---

## 4. Use an operator to read in a file as the source of input values

It's often the case that, when we want to run on multiple inputs, the input values are contained in a file.
As an example, we prepared a CSV file called `greetings.csv` in the `data/` directory, containing several greetings separated by commas.

```csv title="greetings.csv"
Hello,Bonjour,Holà
```

So we need to modify our workflow to read in the values from a file like that.

### 4.1. Switch the input parameter to point to the CSV file

In the workflow block, make the following code change:

_Before:_

```groovy title="hello-channels.nf" linenums="23"
/*
 * Pipeline parameters
 */
params.greet = ['Hello','Bonjour','Holà']
```

_After:_

```groovy title="hello-channels.nf" linenums="23"
/*
 * Pipeline parameters
 */
params.greeting = 'data/greetings.csv'
```

### 4.2. Switch to a channel factory designed to handle a file

Since we now want to use a file instead of a simple value as the input, we can't use the `Channel.of()` channel factory from before.
We need to switch to using a new channel factory, `Channel.fromPath()`, which has some built-in functionality for handling file paths.

In the workflow block, make the following code change:

_Before:_

```groovy title="hello-channels.nf" linenums="29"
// create a channel for inputs
greeting_ch = Channel.of(params.greeting)
                     .flatten()
```

_After:_

```groovy title="hello-channels.nf" linenums="29"
// create a channel for inputs from a CSV file
greeting_ch = Channel.fromPath(params.greeting)
```

### 4.3. Run the workflow

Let's try running the workflow with the new channel factory and the input file.

```bash
nextflow run hello-channels.nf
```

Oh no, another error. This one starts like this:

```console title="Output"
ERROR ~ Error executing process > 'sayHello (1)'

Caused by:
  File `/workspace/gitpod/hello-nextflow/data/greetings.csv-output.txt` is outside the scope of the process work directory: /workspace/gitpod/hello-nextflow/work/cf/8f1c7f433eed3d79e2fa4d1ca81e10

Command executed:

  echo '/workspace/gitpod/hello-nextflow/data/greetings.csv' > '/workspace/gitpod/hello-nextflow/data/greetings.csv-output.txt'
```

The `Command executed:` bit is especially helpful here (you may need to scroll down a bit to find it).

Once again it looks like Nextflow tried to run a single process call, but using the file path itself as a string value.
So it has resolved the file path correctly, but it didn't open the file, which is what we wanted.

How do we get Nextflow to open the file and load its contents into the channel?

Sounds like we need another [operator](https://www.nextflow.io/docs/latest/reference/operator.html).

### 4.4. Add `splitCsv()` operator

Looking through the list of operators again, we find [`splitCsv()`](https://www.nextflow.io/docs/latest/reference/operator.html#splitCsv), which is designed to parse and split CSV-formatted text.

To apply the operator, add it to the channel construction instruction like previously; and we're also going to include view statements while we're at it.

In the workflow block, make the following code change:

_Before:_

```groovy title="hello-channels.nf" linenums="46"
// create a channel for inputs from a CSV file
greeting_ch = Channel.fromPath(params.greeting)
```

_After:_

```groovy title="hello-channels.nf" linenums="46"
// create a channel for inputs from a CSV file
greeting_ch = Channel.fromPath(params.greeting)
                     .view { "Before splitCsv: $it" }
                     .splitCsv()
                     .view { "After splitCsv: $it" }
```

### 4.5. Run the workflow again

Let's try running the workflow with the added CSV-parsing logic.

```bash
nextflow run hello-channels.nf
```

Sadly, this fails too. The console output and error starts like this:

```console title="Output"
Before splitCsv: /workspace/gitpod/hello-nextflow/data/greetings.csv
After splitCsv: [Hello, Bonjour, Holà]
ERROR ~ Error executing process > 'sayHello (1)'

Caused by:
  Missing output file(s) `[Hello, Bonjour, Holà]-output.txt` expected by process `sayHello (1)`
```

Okay, this looks a bit familiar. It looks like what happened earlier the first time we tried to run on an array of values.
So Nextflow has successfully loaded the file contents into the channel, but as a single item.
Once again we're going to need to split it up.

Can we use the same solution as that time?

### 4.6. Add `flatten()` operator

Remember `flatten()`? We love `flatten()`.

You know the drill now: we're going to add the operator to the channel construction instruction, and include another `view()` call.

In the workflow block, make the following code change:

_Before:_

```groovy title="hello-channels.nf" linenums="29"
// create a channel for inputs from a CSV file
greeting_ch = Channel.fromPath(params.greeting)
                     .view { "Before splitCsv: $it" }
                     .splitCsv()
                     .view { "After splitCsv: $it" }
```

_After:_

```groovy title="hello-channels.nf" linenums="29"
// create a channel for inputs from a CSV file
greeting_ch = Channel.fromPath(params.greeting)
                     .view { "Before splitCsv: $it" }
                     .splitCsv()
                     .view { "After splitCsv: $it" }
                     .flatten()
                     .view { "After flatten: $it" }
```

### 4.7. Run the workflow one more time

Let's run it one more time:

```bash
nextflow run hello-channels.nf
```

This time it should run without error!

```console title="Output"
executor >  local (3)
[99/657682] sayHello (3) [100%] 3 of 3 ✔
Before splitCsv: /workspace/gitpod/hello-nextflow/data/greetings.csv
After splitCsv: [Hello, Bonjour, Holà]
After flatten: Hello
After flatten: Bonjour
After flatten: Holà
```

Looking at the output of the `view()` statements, we see the following:

-   A single `Before splitCsv:` statement: at that point the channel contains one item, the original file path.
-   A single `After splitCsv:` statement: at that point the channel still contains only one item, an array containing the three values.
-   Three separate `After flatten:` statements: one for each greeting, which are now individual items in the channel.

We can also look at the output files, which show that each greeting was correctly extracted and processed through the workflow.

We've achieved the same result as previously, but now we have a lot more flexibility to add more elements to the channel of greetings we want to process by modifying an input file, without modifying any code.

!!! note

    Here we had all greetings on one line in the CSV file.
    You can try adding more lines to the CSV file and see what happens, with and without the `flatten()` operator.
    You'll learn how to handle more complex forms of inputs in a later training.

### Takeaway

You know how to use the operators `splitCsv()` and `flatten()` to read in a file of input values and handle them appropriately.

More generally, you have a basic understanding of how Nextflow uses **channels** to manage inputs to processes and **operators** to transform their contents.

### What's next?

Take a break!
When you're ready, move on to Part 3 to learn how to add more steps and connect them together into a proper workflow.
