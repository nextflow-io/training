# Part 3: Hello Workflow

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/zJP7cUYPEbA?si=Irl9nAQniDyICp2b&amp;list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } See [the whole playlist](https://www.youtube.com/playlist?list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik) on the Nextflow YouTube channel.

:green_book: The video transcript is available [here](./transcripts/03_hello_workflow.md).
///

Most real-world workflows involve more than one step.
In this training module, you'll learn how to connect processes together in a multi-step workflow.

This will teach you the Nextflow way of achieving the following:

1. Making data flow from one process to the next
2. Collecting outputs from multiple process calls into a single process call
3. Passing more than one input to a process
4. Handling multiple outputs coming out of a process

To demonstrate, we will continue building on the domain-agnostic Hello World example from Parts 1 and 2.
This time, we're going to make the following changes to our workflow to better reflect how people build actual workflows:

1. Add a second step that converts the greeting to uppercase.
2. Add a third step that collects all the transformed greetings and writes them into a single file.
3. Add a parameter to name the final output file and pass that as a secondary input to the collection step.
4. Make the collection step also output a simple statistic about what was processed.

---

## 0. Warmup: Run `hello-workflow.nf`

We're going to use the workflow script `hello-workflow.nf` as a starting point.
It is equivalent to the script produced by working through Part 2 of this training course.

Just to make sure everything is working, run the script once before making any changes:

```bash
nextflow run hello-workflow.nf
```

```console title="Output"
 N E X T F L O W   ~  version 25.04.3

Launching `hello-workflow.nf` [stupefied_sammet] DSL2 - revision: b9e466930b

executor >  local (3)
[2a/324ce6] sayHello (3) | 3 of 3 ✔
```

As previously, you will find the output files in the `results` directory (specified by the `publishDir` directive).

```console title="Directory contents"
results
├── Bonjour-output.txt
├── Hello-output.txt
└── Holà-output.txt
```

!!! note

    There may also be a file named `output.txt` left over if you worked through Part 2 in the same environment.

If that worked for you, you're ready to learn how to assemble a multi-step workflow.

---

## 1. Add a second step to the workflow

We're going to add a step to convert the greeting to uppercase.
To that end, we need to do three things:

- Define the command we're going to use to do the uppercase conversion.
- Write a new process that wraps the uppercasing command.
- Call the new process in the workflow block and set it up to take the output of the `sayHello()` process as input.

### 1.1. Define the uppercasing command and test it in the terminal

To do the conversion of the greetings to uppercase, we're going to use a classic UNIX tool called `tr` for 'text replacement', with the following syntax:

```bash title="Syntax"
tr '[a-z]' '[A-Z]'
```

This is a very naive text replacement one-liner that does not account for accented letters, so for example 'Holà' will become 'HOLà', but it will do a good enough job for demonstrating the Nextflow concepts and that's what matters.

To test it out, we can run the `echo 'Hello World'` command and pipe its output to the `tr` command:

```bash
echo 'Hello World' | tr '[a-z]' '[A-Z]' > UPPER-output.txt
```

The output is a text file called `UPPER-output.txt` that contains the uppercase version of the `Hello World` string:

```console title="UPPER-output.txt"
HELLO WORLD
```

That's basically what we're going to try to do with our workflow.

### 1.2. Write the uppercasing step as a Nextflow process

We can model our new process on the first one, since we want to use all the same components.

Add the following process definition to the workflow script:

```groovy title="hello-workflow.nf" linenums="22"
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

Here, we compose the second output filename based on the input filename, similarly to what we did originally for the output of the first process.

!!! note

    Nextflow will determine the order of operations based on the chaining of inputs and outputs, so the order of the process definitions in the workflow script does not matter.
    However, we do recommend you be kind to your collaborators and to your future self, and try to write them in a logical order for the sake of readability.

### 1.3. Add a call to the new process in the workflow block

Now we need to tell Nextflow to actually call the process that we just defined.

In the workflow block, make the following code change:

=== "After"

    ```groovy title="hello-workflow.nf" linenums="53" hl_lines="4 5"
        // emit a greeting
        sayHello(greeting_ch)

        // convert the greeting to uppercase
        convertToUpper()
    }
    ```

=== "Before"

    ```groovy title="hello-workflow.nf" linenums="53"
        // emit a greeting
        sayHello(greeting_ch)
    }
    ```

This is not yet functional because we have not specified what should be input to the `convertToUpper()` process.

### 1.4. Pass the output of the first process to the second process

Now we need to make the output of the `sayHello()` process flow into the `convertToUpper()` process.

Conveniently, Nextflow automatically packages the output of a process into a channel called `<process>.out`.
So the output of the `sayHello` process is a channel called `sayHello.out`, which we can plug straight into the call to `convertToUpper()`.

In the workflow block, make the following code change:

=== "After"

    ```groovy title="hello-workflow.nf" linenums="56" hl_lines="2"
        // convert the greeting to uppercase
        convertToUpper(sayHello.out)
    }
    ```

=== "Before"

    ```groovy title="hello-workflow.nf" linenums="56"
        // convert the greeting to uppercase
        convertToUpper()
    }
    ```

For a simple case like this (one output to one input), that's all we need to do to connect two processes!

### 1.5. Run the workflow again with `-resume`

Let's run this using the `-resume` flag, since we've already run the first step of the workflow successfully.

```bash
nextflow run hello-workflow.nf -resume
```

You should see the following output:

```console title="Output" linenums="1"
 N E X T F L O W   ~  version 25.04.3

Launching `hello-workflow.nf` [disturbed_darwin] DSL2 - revision: 4e252c048f

executor >  local (3)
[79/33b2f0] sayHello (2)       | 3 of 3, cached: 3 ✔
[b3/d52708] convertToUpper (3) | 3 of 3 ✔
```

There is now an extra line in the console output (line 7), which corresponds to the new process we just added.

Let's have a look inside the work directory of one of the calls to the second process.

```console title="Directory contents"
work/b3/d52708edba8b864024589285cb3445/
├── Bonjour-output.txt -> /workspaces/training/hello-nextflow/work/79/33b2f0af8438486258d200045bd9e8/Bonjour-output.txt
└── UPPER-Bonjour-output.txt
```

We find two output files: the output of the first process AND the output of the second.

The output of the first process is in there because Nextflow staged it there in order to have everything needed for execution within the same subdirectory.
However, it is actually a symbolic link pointing to the the original file in the subdirectory of the first process call.
By default, when running on a single machine as we're doing here, Nextflow uses symbolic links rather than copies to stage input and intermediate files.

You'll also find the final outputs in the `results` directory since we used the `publishDir` directive in the second process too.

```console title="Directory contents"
results
├── Bonjour-output.txt
├── Hello-output.txt
├── Holà-output.txt
├── UPPER-Bonjour-output.txt
├── UPPER-Hello-output.txt
└── UPPER-Holà-output.txt
```

Think about how all we did was connect the output of `sayHello` to the input of `convertToUpper` and the two processes could be run in series.
Nextflow did the hard work of handling individual input and output files and passing them between the two commands for us.

This is one of the reasons Nextflow channels are so powerful: they take care of the busywork involved in connecting workflow steps together.

### Takeaway

You know how to add a second step that takes the output of the first step as input.

### What's next?

Learn how to collect outputs from batched process calls and feed them into a single process.

---

## 2. Add a third step to collect all the greetings

When we use a process to apply a transformation to each of the elements in a channel, like we're doing here to the multiple greetings, we sometimes want to collect elements from the output channel of that process, and feed them into another process that performs some kind of analysis or summation.

In the next step we're simply going to write all the elements of a channel to a single file, using the UNIX `cat` command.

### 2.1. Define the collection command and test it in the terminal

The collection step we want to add to our workflow will use the `cat` command to concatenate multiple uppercased greetings into a single file.

Let's run the command by itself in the terminal to verify that it works as expected, just like we've done previously.

Run the following in your terminal:

```bash
echo 'Hello' | tr '[a-z]' '[A-Z]' > UPPER-Hello-output.txt
echo 'Bonjour' | tr '[a-z]' '[A-Z]' > UPPER-Bonjour-output.txt
echo 'Holà' | tr '[a-z]' '[A-Z]' > UPPER-Holà-output.txt
cat UPPER-Hello-output.txt UPPER-Bonjour-output.txt UPPER-Holà-output.txt > COLLECTED-output.txt
```

The output is a text file called `COLLECTED-output.txt` that contains the uppercase versions of the original greetings.

```console title="COLLECTED-output.txt"
HELLO
BONJOUR
HOLà
```

That is the result we want to achieve with our workflow.

### 2.2. Create a new process to do the collection step

Let's create a new process and call it `collectGreetings()`.
We can start writing it based on the previous one.

#### 2.2.1. Write the 'obvious' parts of the process

Add the following process definition to the workflow script:

```groovy title="hello-workflow.nf" linenums="41"
/*
 * Collect uppercase greetings into a single output file
 */
process collectGreetings {

    publishDir 'results', mode: 'copy'

    input:
        ???

    output:
        path "COLLECTED-output.txt"

    script:
    """
    ??? > 'COLLECTED-output.txt'
    """
}
```

This is what we can write with confidence based on what you've learned so far.
But this is not functional!
It leaves out the input definition(s) and the first half of the script command because we need to figure out how to write that.

#### 2.2.2. Define inputs to `collectGreetings()`

We need to collect the greetings from all the calls to the `convertToUpper()` process.
What do we know we can get from the previous step in the workflow?

The channel output by `convertToUpper()` will contain the paths to the individual files containing the uppercased greetings.
That amounts to one input slot; let's call it `input_files` for simplicity.

In the process block, make the following code change:

=== "After"

    ```groovy title="hello-workflow.nf" linenums="48" hl_lines="2"
            input:
                path input_files
    ```

=== "Before"

    ```groovy title="hello-workflow.nf" linenums="48"
            input:
                ???
    ```

Notice we use the `path` prefix even though we expect this to contain multiple files.
Nextflow doesn't mind, so it doesn't matter.

#### 2.2.3. Compose the concatenation command

This is where things could get a little tricky, because we need to be able to handle an arbitrary number of input files.
Specifically, we can't write the command up front, so we need to tell Nextflow how to compose it at runtime based on what inputs flow into the process.

In other words, if we have an input channel containing the element `[file1.txt, file2.txt, file3.txt]`, we need Nextflow to turn that into `cat file1.txt file2.txt file3.txt`.

Fortunately, Nextflow is quite happy to do that for us if we simply write `cat ${input_files}` in the script command.

In the process block, make the following code change:

=== "After"

    ```groovy title="hello-workflow.nf" linenums="54" hl_lines="3"
        script:
        """
        cat ${input_files} > 'COLLECTED-output.txt'
        """
    ```

=== "Before"

    ```groovy title="hello-workflow.nf" linenums="54"
        script:
        """
        ??? > 'COLLECTED-output.txt'
        """
    ```

In theory this should handle any arbitrary number of input files.

!!! tip

    Some command-line tools require providing an argument (like `-input`) for each input file.
    In that case, we would have to do a little bit of extra work to compose the command.
    You can see an example of this in the [Nextflow for Genomics](../../nf4_science/genomics/) training course.

<!--[ADD LINK to note above] -->

### 2.3. Add the collection step to the workflow

Now we should just need to call the collection process on the output of the uppercasing step.

#### 2.3.1. Connect the process calls

In the workflow block, make the following code change:

=== "After"

    ```groovy title="hello-workflow.nf" linenums="75" hl_lines="4 5"
        // convert the greeting to uppercase
        convertToUpper(sayHello.out)

        // collect all the greetings into one file
        collectGreetings(convertToUpper.out)
    }
    ```

=== "Before"

    ```groovy title="hello-workflow.nf" linenums="75"
        // convert the greeting to uppercase
        convertToUpper(sayHello.out)
    }
    ```

This connects the output of `convertToUpper()` to the input of `collectGreetings()`.

#### 2.3.2. Run the workflow with `-resume`

Let's try it.

```bash
nextflow run hello-workflow.nf -resume
```

It runs successfully, including the third step:

```console title="Output" linenums="1"
 N E X T F L O W   ~  version 25.04.3

Launching `hello-workflow.nf` [mad_gilbert] DSL2 - revision: 6acfd5e28d

executor >  local (3)
[79/33b2f0] sayHello (2)         | 3 of 3, cached: 3 ✔
[99/79394f] convertToUpper (3)   | 3 of 3, cached: 3 ✔
[47/50fe4a] collectGreetings (1) | 3 of 3 ✔
```

However, look at the number of calls for `collectGreetings()` on line 8.
We were only expecting one, but there are three.

And have a look at the contents of the final output file too:

```console title="results/COLLECTED-output.txt"
Holà
```

Oh no. The collection step was run individually on each greeting, which is NOT what we wanted.

We need to do something to tell Nextflow explicitly that we want that third step to run on all the elements in the channel output by `convertToUpper()`.

### 2.4. Use an operator to collect the greetings into a single input

Yes, once again the answer to our problem is an operator.

Specifically, we are going to use the aptly-named [`collect()`](https://www.nextflow.io/docs/latest/reference/operator.html#collect) operator.

#### 2.4.1. Add the `collect()` operator

This time it's going to look a bit different because we're not adding the operator in the context of a channel factory, but to an output channel.

We take the `convertToUpper.out` and append the `collect()` operator, which gives us `convertToUpper.out.collect()`.
We can plug that directly into the `collectGreetings()` process call.

In the workflow block, make the following code change:

=== "After"

    ```groovy title="hello-workflow.nf" linenums="78" hl_lines="2"
        // collect all the greetings into one file
        collectGreetings(convertToUpper.out.collect())
    }
    ```

=== "Before"

    ```groovy title="hello-workflow.nf" linenums="78"
        // collect all the greetings into one file
        collectGreetings(convertToUpper.out)
    }
    ```

#### 2.4.2. Add some `view()` statements

Let's also include a couple of `view()` statements to visualize the before and after states of the channel contents.

=== "After"

    ```groovy title="hello-workflow.nf" linenums="78" hl_lines="4 6"
        // collect all the greetings into one file
        collectGreetings(convertToUpper.out.collect())

        // optional view statements
        convertToUpper.out.view { greeting -> "Before collect: $greeting" }
        convertToUpper.out.collect().view { greeting -> "After collect: $greeting" }
    }
    ```

=== "Before"

    ```groovy title="hello-workflow.nf" linenums="78"
        // collect all the greetings into one file
        collectGreetings(convertToUpper.out.collect())
    }
    ```

The `view()` statements can go anywhere you want; we put them after the call for readability.

#### 2.4.3. Run the workflow again with `-resume`

Let's try it:

```bash
nextflow run hello-workflow.nf -resume
```

It runs successfully, although the log output may look a little messier than this (we cleaned it up for readability).

```console title="Output" linenums="1"
 N E X T F L O W   ~  version 25.04.3

Launching `hello-workflow.nf` [soggy_franklin] DSL2 - revision: bc8e1b2726

[d6/cdf466] sayHello (1)       | 3 of 3, cached: 3 ✔
[99/79394f] convertToUpper (2) | 3 of 3, cached: 3 ✔
[1e/83586c] collectGreetings   | 1 of 1 ✔
Before collect: /workspaces/training/hello-nextflow/work/b3/d52708edba8b864024589285cb3445/UPPER-Bonjour-output.txt
Before collect: /workspaces/training/hello-nextflow/work/99/79394f549e3040dfc2440f69ede1fc/UPPER-Hello-output.txt
Before collect: /workspaces/training/hello-nextflow/work/aa/56bfe7cf00239dc5badc1d04b60ac4/UPPER-Holà-output.txt
After collect: [/workspaces/training/hello-nextflow/work/b3/d52708edba8b864024589285cb3445/UPPER-Bonjour-output.txt, /workspaces/training/hello-nextflow/work/99/79394f549e3040dfc2440f69ede1fc/UPPER-Hello-output.txt, /workspaces/training/hello-nextflow/work/aa/56bfe7cf00239dc5badc1d04b60ac4/UPPER-Holà-output.txt]
```

This time the third step was only called once!

Looking at the output of the `view()` statements, we see the following:

- Three `Before collect:` statements, one for each greeting: at that point the file paths are individual items in the channel.
- A single `After collect:` statement: the three file paths are now packaged into a single element.

Have a look at the contents of the final output file too:

```console title="results/COLLECTED-output.txt"
BONJOUR
HELLO
HOLà
```

This time we have all three greetings in the final output file. Success! Remove the optional view calls to make the next outputs less verbose.

!!! note

    If you run this several times without `-resume`, you will see that the order of the greetings changes from one run to the next.
    This shows you that the order in which elements flow through process calls is not guaranteed to be consistent.

### Takeaway

You know how to collect outputs from a batch of process calls and feed them into a joint analysis or summation step.

### What's next?

Learn how to pass more than one input to a process.

---

## 3. Pass more than one input to a process in order to name the final output file uniquely

We want to be able to name the final output file something specific in order to process subsequent batches of greetings without overwriting the final results.

To that end, we're going to make the following refinements to the workflow:

- Modify the collector process to accept a user-defined name for the output file
- Add a command-line parameter to the workflow and pass it to the collector process

### 3.1. Modify the collector process to accept a user-defined name for the output file

We're going to need to declare the additional input and integrate it into the output file name.

#### 3.1.1. Declare the additional input in the process definition

Good news: we can declare as many input variables as we want.
Let's call this one `batch_name`.

In the process block, make the following code change:

=== "After"

    ```groovy title="hello-workflow.nf" linenums="48" hl_lines="3"
        input:
            path input_files
            val batch_name
    ```

=== "Before"

    ```groovy title="hello-workflow.nf" linenums="48"
        input:
            path input_files
    ```

You can set up your processes to expect as many inputs as you want.
Later on, you will learn how to manage required vs. optional inputs.

#### 3.1.2. Use the `batch_name` variable in the output file name

In the process block, make the following code change:

=== "After"

    ```groovy title="hello-workflow.nf" linenums="52" hl_lines="2 6"
        output:
            path "COLLECTED-${batch_name}-output.txt"

        script:
        """
        cat ${input_files} > 'COLLECTED-${batch_name}-output.txt'
        """
    ```

=== "Before"

    ```groovy title="hello-workflow.nf" linenums="52"
        output:
            path "COLLECTED-output.txt"

        script:
        """
        cat ${input_files} > 'COLLECTED-output.txt'
        """
    ```

This sets up the process to use the `batch_name` value to generate a specific filename for the final output of the workflow.

### 3.2. Add a `batch` command-line parameter

Now we need a way to supply the value for `batch_name` and feed it to the process call.

#### 3.2.1. Use `params` to set up the parameter

You already know how to use the `params` system to declare CLI parameters.
Let's use that to declare a `batch` parameter (with a default value because we are lazy).

In the pipeline parameters section, make the following code changes:

=== "After"

    ```groovy title="hello-workflow.nf" linenums="61" hl_lines="5"
    /*
     * Pipeline parameters
     */
    params.greeting = 'greetings.csv'
    params.batch = 'test-batch'
    ```

=== "Before"

    ```groovy title="hello-workflow.nf" linenums="61"
    /*
     * Pipeline parameters
     */
    params.greeting = 'greetings.csv'
    ```

Remember you can override that default value by specifying a value with `--batch` on the command line.

#### 3.2.2. Pass the `batch` parameter to the process

To provide the value of the parameter to the process, we need to add it in the process call.

In the workflow block, make the following code change:

=== "After"

    ```groovy title="hello-workflow.nf" linenums="80" hl_lines="2"
        // collect all the greetings into one file
        collectGreetings(convertToUpper.out.collect(), params.batch)
    ```

=== "Before"

    ```groovy title="hello-workflow.nf" linenums="80"
        // collect all the greetings into one file
        collectGreetings(convertToUpper.out.collect())
    ```

!!! warning

    You MUST provide the inputs to a process in the EXACT SAME ORDER as they are listed in the input definition block of the process.

### 3.3. Run the workflow

Let's try running this with a batch name on the command line.

```bash
nextflow run hello-workflow.nf -resume --batch trio
```

It runs successfully:

```console title="Output" linenums="1"
 N E X T F L O W   ~  version 25.04.3

Launching `hello-workflow.nf` [confident_rutherford] DSL2 - revision: bc58af409c

executor >  local (1)
[79/33b2f0] sayHello (2)       | 3 of 3, cached: 3 ✔
[99/79394f] convertToUpper (2) | 3 of 3, cached: 3 ✔
[b5/f19efe] collectGreetings   | 1 of 1 ✔
```

And produces the desired output:

```console title="bash"
cat results/COLLECTED-trio-output.txt
```

```console title="Output"
HELLO
BONJOUR
HOLà
```

Now, subsequent runs on other batches of inputs won't clobber previous results (as long as we specify the parameter appropriately).

### Takeaway

You know how to pass more than one input to a process.

### What's next?

Learn how to emit multiple outputs and handle them conveniently.

---

## 4. Add an output to the collector step

When a process produces only one output, it's easy to access it (in the workflow block) using the `<process>.out` syntax.
When there are two or more outputs, the default way to select a specific output is to use the corresponding (zero-based) index; for example, you would use `<process>.out[0]` to get the first output.
This is not terribly convenient; it's too easy to grab the wrong index.

Let's have a look at how we can select and use a specific output of a process when there are more than one.

For demonstration purposes, let's say we want to count and report the number of greetings that are being collected for a given batch of inputs.

To that end, we're going to make the following refinements to the workflow:

- Modify the process to count and output the number of greetings
- Once the process has run, select the count and report it using `view` (in the workflow block)

### 4.1. Modify the process to count and output the number of greetings

This will require two key changes to the process definition: we need a way to count the greetings, then we need to add that count to the `output` block of the process.

#### 4.1.1. Count the number of greetings collected

Conveniently, Nextflow lets us add arbitrary code in the `script:` block of the process definition, which comes in really handy for doing things like this.

That means we can use the built-in `size()` function to get the number of files in the `input_files` array.

In the `collectGreetings` process block, make the following code change:

=== "After"

    ```groovy title="hello-workflow.nf" linenums="55" hl_lines="2"
        script:
            count_greetings = input_files.size()
        """
        cat ${input_files} > 'COLLECTED-${batch_name}-output.txt'
        """
    ```

=== "Before"

    ```groovy title="hello-workflow.nf" linenums="55"
        script:
        """
        cat ${input_files} > 'COLLECTED-${batch_name}-output.txt'
        """
    ```

The `count_greetings` variable will be computed at runtime.

#### 4.1.2. Emit the count as a named output

In principle all we need to do is to add the `count_greetings` variable to the `output:` block.

However, while we're at it, we're also going to add some `emit:` tags to our output declarations. These will enable us to select the outputs by name instead of having to use positional indices.

In the process block, make the following code change:

=== "After"

    ```groovy title="hello-workflow.nf" linenums="52" hl_lines="2 3"
        output:
            path "COLLECTED-${batch_name}-output.txt" , emit: outfile
            val count_greetings , emit: count
    ```

=== "Before"

    ```groovy title="hello-workflow.nf" linenums="52"
        output:
            path "COLLECTED-${batch_name}-output.txt"
    ```

The `emit:` tags are optional, and we could have added a tag to only one of the outputs.
But as the saying goes, why not both?

### 4.2. Report the output at the end of the workflow

Now that we have two outputs coming out of the `collectGreetings` process, the `collectGreetings.out` output contains two channels:

- `collectGreetings.out.outfile` contains the final output file
- `collectGreetings.out.count` contains the count of greetings

We could send either or both of these to another process for further work. However, in the interest of wrapping this up, we're just going to use `view()` to demonstrate that we can access and report the count of greetings.

In the workflow block, make the following code change:

=== "After"

    ```groovy title="hello-workflow.nf" linenums="82" hl_lines="4 5"
        // collect all the greetings into one file
        collectGreetings(convertToUpper.out.collect(), params.batch)

        // emit a message about the size of the batch
        collectGreetings.out.count.view { num_greetings -> "There were $num_greetings greetings in this batch" }
    ```

=== "Before"

    ```groovy title="hello-workflow.nf" linenums="82"
        // collect all the greetings into one file
        collectGreetings(convertToUpper.out.collect(), params.batch)
    ```

!!! note

    There are a few other ways we could achieve a similar result, including some more elegant ones like the `count()` operator, but this allows us to show how to handle multiple outputs, which is what we care about.

### 4.3. Run the workflow

Let's try running this with the current batch of greetings.

```bash
nextflow run hello-workflow.nf -resume --batch trio
```

This runs successfully:

```console title="Output" linenums="1"
 N E X T F L O W   ~  version 25.04.3

Launching `hello-workflow.nf` [evil_sinoussi] DSL2 - revision: eeca64cdb1

[d6/cdf466] sayHello (1)       | 3 of 3, cached: 3 ✔
[99/79394f] convertToUpper (2) | 3 of 3, cached: 3 ✔
[9e/1dfda7] collectGreetings   | 1 of 1, cached: 1 ✔
There were 3 greetings in this batch
```

The last line (line 8) shows that we correctly retrieved the count of greetings processed.
Feel free to add more greetings to the CSV and see what happens.

### Takeaway

You know how to make a process emit a named output and how to access it from the workflow block.

More generally, you understand the key principles involved in connecting processes together in common ways.

### What's next?

Take an extra long break, you've earned it.
When you're ready, move on to Part 4 to learn how to modularize your code for better maintainability and code efficiency.
