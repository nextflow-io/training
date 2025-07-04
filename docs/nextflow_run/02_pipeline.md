# Part 2: Run pipelines

In Part 1 of this course (Run Basic Operations), we started with an example workflow that had only minimal features in order to keep the code complexity low.
However, most real-world pipelines use more sophisticated features in order to enable efficient processing of large amounts of data at scale, and apply multiple processing steps chained together by sometimes complex logic.

In this part of the training, we demonstrate key features of real-world pipelines through a set of example workflows that build on the original Hello World pipeline.

## 1. Processing input data from a file

The `hello.nf` workflow we ran in Part 1 used a command-line parameter (`--greeting`) to provide a single value at a time.
That was a deliberately simplified approach.

In a real-world pipeline, we typically want to process multiple data points (or data series) contained in one or more input files.
And wherever possible, we want to run the processing of independent data in parallel, to shorten the time spent waiting for analysis.

To enable this efficiently, Nextflow uses a system of queues called **channels**.

To demonstrate this, we've prepared a a CSV file called `greetings.csv` that contains several input greetings, mimicking the kind of columnar data you might want to process in a real data analysis.

```csv title="greetings.csv" linenums="1"
Hello,English,123
Bonjour,French,456
Holà,Spanish,789
```

_The numbers are not meaningful, they are just there for illustrative purposes._

And we've written an improved version of the original workflow, now called `channel.nf`, that will read in the CSV file, extract the greetings and write each of them to a separate file.

<!-- TODO: add diagram of operations -->

Let's run the workflow first, and we'll take a look at the relevant Nextflow code afterward.

### 1.1. Run the workflow

Run the following command in your terminal.

```bash
nextflow run channel.nf --input greetings.csv
```

This should run without error.

```console title="Output" linenums="1"
 N E X T F L O W   ~  version 25.04.3

Launching `channel.nf` [mighty_sammet] DSL2 - revision: 29fb5352b3

executor >  local (3)
[8e/0eb066] sayHello (2) [100%] 3 of 3 ✔
```

Excitingly, this seems to indicate that '3 of 3' calls were made for the process, which is encouraging, since there were three rows of data in the CSV we provided as input.
This suggests the sayHello() process was called three times, once on each input row.

### 1.2. Find the outputs in the `results` directory

Let's look at the 'results' directory to see if our workflow is still writing a copy of our outputs there.

````console title="results/" linenums="1"
results
├── Bonjour-output.txt
├── Hello-output.txt
└── Holà-output.txt

Yes! We see three output files with different names, conveniently enough.
(Spoiler: we changed the workflow to name the files differently.)

You can open each of them to satisfy yourself that they contain the appropriate greeting string.

```console title="results/Hello-output.txt"
Hello
````

```console title="results/Bonjour-output.txt"
Bonjour
```

```console title="results/Holà-output.txt"
Holà
```

This confirms each greeting in the input file has been processed appropriately.

### 1.3. Find the original outputs and logs

You may have noticed that the console output above referred to only one task directory.
Does that mean all three calls to `sayHello()` were executed within that one task directory?

Let's have a look inside that `8e/0eb066` task directory:

```console title="8e/0eb066"
work/8e/0eb066071cdb4123906b7b4ea8b047/
└── Bonjour-output.txt
```

No! We only find the output corresponding to one of the greetings (as well as the accessory files if we enable display of hidden files).

So what's going on here?

By default, the ANSI logging system writes the status information for all calls to the same process on the same line.
As a result, it only showed us one of the three task directory paths (`8e/0eb066`) in the console output.
There are two others that are not listed there.

We can modify the logging behavior to see the full list of process calls by adding the `-ansi-log false` to the command as follows:

```bash
nextflow run channel.nf --input greetings.csv -ansi-log false
```

This time we see all three process runs and their associated work subdirectories listed in the output:

```console title="Output" linenums="1"
N E X T F L O W  ~  version 25.04.3
Launching `channel.nf` [pedantic_hamilton] DSL2 - revision: 6bbc42e49f
[ab/1a8ece] Submitted process > sayHello (1)
[0d/2cae24] Submitted process > sayHello (2)
[b5/0df1d6] Submitted process > sayHello (3)
```

This confirms that the `sayHello()` process gets called three times, and a separate task directory is created for each one.

!!! note

    For a complex workflow, or a large number of inputs, having the full list output to the terminal might get a bit overwhelming, so you might prefer not to use `-ansi-log false` in those cases.

    Note also that the way the status is reported is a bit different between the two logging modes.
    In the condensed mode, Nextflow reports whether calls were completed successfully or not.
    In this expanded mode, it only reports that they were submitted.

If we look inside each of the task directories listed there, we can confirm that each one corresponds to one of the greetings, .

```console title="ab/1a8ece"
work/ab/1a8ece307e53f03fce689dde904b64/
└── Hello-output.txt
```

```console title="0d/2cae24"
work/0d/2cae2481a53593bc607077c80c9466/
└── Bonjour-output.txt
```

```console title="b5/0df1d6"
work/b5/0df1d642353269909c2ce23fc2a8fa/
└── Holà-output.txt
```

This confirms that each process call is executed in isolation from all the others.
That has many advantages, including avoiding collisions if the process produces any intermediate files with non-unique names.

### 1.4. Examine the code

So this version of the workflow is capable of reading in a CSV file of inputs, processing the inputs separately, and naming the outputs uniquely.

Let's take a look at what makes that possible in the workflow code.
Once again, we're not aiming to memorize code syntax, but to identify signature components of the workflow that provide important functionality.

```groovy title="channel.nf" linenums="1"
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

workflow {

    // create a channel for inputs from a CSV file
    greeting_ch = Channel.fromPath(params.input)
                        .splitCsv()
                        .map { line -> line[0] }

    // emit a greeting
    sayHello(greeting_ch)
}
```

#### 1.4.1. Load the inputs from the CSV

This is the most interesting part: how did we switch from taking a single value from the command-line, to taking a CSV file, parsing it and processing the individual greetings it contains?

In Nextflow, we do that with a **channel**: a construct designed to handle inputs efficiently and shuttle them from one step to another in multi-step workflows, while providing built-in parallelism and many additional benefits.

Let's break it down.

```groovy title="channel.nf" linenums="22"
workflow {

    // create a channel for inputs from a CSV file
    greeting_ch = Channel.fromPath(params.greeting)
                        .splitCsv()
                        .map { line -> line[0] }
```

This is where the magic happens, starting at line 25.
Here's what that line means in plain English:

Channel - create a **channel**, i.e. a queue that will hold the data,
.fromPath - from a filepath
(params.input) - provided with `--input` on the command line

In other words, that line tells Nextflow: take the filepath given with `--input` and get ready to treat its contents as input data.

Then the next two lines apply **operators** that do the actual parsing of the file and loading of the data into the appropriate data structure:

.splitCsv() - parse the CSV file into an array representing rows and columns
.map { line -> line[0] } - for each row (line), take only the element in the first column

So in practice, starting from the following CSV file:

```csv title="greetings.csv" linenums="1"
Hello,English,123
Bonjour,French,456
Holà,Spanish,789
```

We have transformed that into an array that looks like this:

```txt title="Array contents"
[[Hello,English,123],[Bonjour,French,456],[Holà,Spanish,789]]
```

And then we've taken the first element from each of the three rows and loaded them into a Nextflow channel that now contains: `Hello`, `Bonjour`, and `Holà`.

In other words, the result of this very short snippet of code is a channel called `greeting_ch` loaded with the three individual greetings from the CSV file, ready for processing.

#### 1.4.2. Call the process on each greeting

Next, in the last line of the workflow block, we provide the loaded `greeting_ch` channel as input to the `sayHello()` process.

```groovy title="channel.nf" linenums="28"
    sayHello(greeting_ch)
}
```

This tells Nextflow to run the process _individually_ on each element in the channel, i.e. on each greeting.

And because Nextflow is smart like that, it will run these process calls in parallel if possible, depending on the available computing infrastructure.

That is how you can achieve efficient and scalable processing of a lot of data (many samples, or data points, whatever is your unit of research) with comparatively very little code.

#### 1.4.3. Ensure the outputs are uniquely named

Finally, it's worth taking a quick look at how we get the output files to be named uniquely.

```groovy title="channel.nf" linenums="13"
    output:
        path "${greeting}-output.txt"

    script:
    """
    echo '$greeting' > '$greeting-output.txt'
    """
```

You see that, compared to the version of this process in `hello.nf`, the output declaration and the relevant bit of the command have changed to include the greeting value in the output file name.
This is one way to ensure that the output file names won't collide when they get published to the common `results` directory.

And that's the only change we've had to make inside the process declaration.

### Takeaway

You understand at a basic level how channels and operators enable us to process multiple inputs efficiently.

### What's next?

Discover how multi-step workflows are constructed and how they operate.

---

## 2. Multi-step workflows

Most real-world workflows involve more than one step.
Let's build on what we just learned about channels, and look at how Nextflow uses channels and operators to connect processes together in a multi-step workflow.

To that end, we provide you with an example workflow that chains together three separate steps and demonstrates the following:

1. Making data flow from one process to the next
2. Collecting outputs from multiple process calls into a single process call

Specifically, we made an expanded version of the workflow called `pipeline.nf` that takes each input greeting, converts it to uppercase, then collects all the uppercased greetings into a single output file.

<!-- TODO: add diagram of operations -->

As previously, we'll run the workflow first then look at the code to see what's changed.

### 2.1. Run the workflow

Run the following command in your terminal:

```bash
nextflow run pipeline.nf --input greetings.csv
```

Once again this should run successfully.

```console title="Output" linenums="1"
 N E X T F L O W   ~  version 25.04.3

Launching `pipeline.nf` [soggy_franklin] DSL2 - revision: bc8e1b2726

[d6/cdf466] sayHello (1)       | 3 of 3 ✔
[99/79394f] convertToUpper (2) | 3 of 3 ✔
[1e/83586c] collectGreetings   | 1 of 1 ✔
```

You see that as promised, multiple steps were run as part of the workflow; the first two (`sayHello` and `convertToUpper`) were presumably run on each individual greeting, and the third (`collectGreetings`) will have been run only once, on the outputs of all three of the `convertToUpper` calls.

### 2.2. Find the outputs

Let's verify that that is in fact what happened by taking a look in the `results` directory.

```console title="Directory contents"
results
├── Bonjour-output.txt
├── COLLECTED-output.txt
├── Hello-output.txt
├── Holà-output.txt
├── UPPER-Bonjour-output.txt
├── UPPER-Hello-output.txt
└── UPPER-Holà-output.txt
```

Look at the file names and check their contents to confirm that they are what you expect; for example:

```console title="bash"
cat results/COLLECTED-output.txt
```

```console title="Output"
HELLO
BONJOUR
HOLà
```

That is the expected final result of our multi-step pipeline.

### 2.3. Examine the code

Let's look at the code and see what we can tie back to what we just observed.

```groovy title="channel.nf" linenums="1"
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

/*
 * Collect uppercase greetings into a single output file
 */
process collectGreetings {

    publishDir 'results', mode: 'copy'

    input:
        path input_files

    output:
        path "COLLECTED-output.txt"

    script:
    """
    cat ${input_files} > 'COLLECTED-output.txt'
    """
}

workflow {

    // create a channel for inputs from a CSV file
    greeting_ch = Channel.fromPath(params.greeting)
                        .splitCsv()
                        .map { line -> line[0] }

    // emit a greeting
    sayHello(greeting_ch)

    // convert the greeting to uppercase
    convertToUpper(sayHello.out)

    // collect all the greetings into one file
    collectGreetings(convertToUpper.out.collect())
}

```

The most obvious difference compared to the previous version of the workflow is that now there are multiple process definitions, and correspondingly, several process calls in the workflow block.

#### 2.3.1. Multiple process definitions

In addition to the original `sayHello` process, we now also have `convertToUpper` and `collectGreetings`, which match the names of the processes we saw in the console output.

All three are structured in the same way and follow roughly the same logic.
We won't go into that in detail, but it shows how a process can be given additional parameters and emit multiple outputs.

#### 2.3.2. Processes are connected via channels

The really interesting thing to look at here is how the process calls are chained together in the workflow block.

```groovy title="channel.nf" linenums="69"
workflow {

    // create a channel for inputs from a CSV file
    greeting_ch = Channel.fromPath(params.input)
                        .splitCsv()
                        .map { line -> line[0] }

    // emit a greeting
    sayHello(greeting_ch)

    // convert the greeting to uppercase
    convertToUpper(sayHello.out)

    // collect all the greetings into one file
    collectGreetings(convertToUpper.out.collect())
}
```

You can see that the first process call, `sayHello(greeting_ch)`, is unchanged.

Then the next process call, to `convertToUpper`, _refers_ to the output of `sayHello` as `sayHello.out`:

```groovy title="channel.nf" linenums="79"
    // convert the greeting to uppercase
    convertToUpper(sayHello.out)
```

This tells Nextflow to provide `sayHello.out`, which represents a channel output by `sayHello()`, as an input to `convertToUpper`.

That is, at its simplest, how we shuttle data from one step to the next in Nextflow.

Then the next call is doing the same thing, with a twist:

```groovy title="channel.nf" linenums="82"
    // collect all the greetings into one file
    collectGreetings(convertToUpper.out.collect())
```

This one is a bit more complicated and deserves its own discussion.

#### 2.3.3. Operators provide additional wiring options

What we're seeing in `convertToUpper.out.collect()` is the use of another operator (like `splitCsv` and `map` in the previous section), called `collect()`.
This operator is used to collect the outputs from multiple calls to the same process (as when we run `sayHello` on multiple greetings independently) and package them into a single channel element.

This allows us to take all the separate uppercased greetings produced by the second step of the workflow and feed them all together to a single call in the third step of the pipeline.
If we didn't apply `collect()` to the output of `convertToUpper()` before feeding it to `collectGreetings()`, Nextflow would simply run `collectGreetings()` independently on each greeting, which would not achieve our goal.

<!-- TODO: add diagram of operations -->

There are many other operators available to apply transformations to the contents of channels between process calls.

This gives pipeline developers a lot of flexibility for customizing the flow logic of their pipeline.
The downside is that it can sometimes make it harder to decipher what the pipeline is doing.

### 2.4. Use the graph preview

One very helpful tool for understanding what a pipeline does, if it's not adequately documented, is the graph preview functionality available in VSCode. You can see this in the training environment by clicking on the small `DAG preview` link displayed just above the workflow block in any Nextflow script.

<!-- TODO: add screenshot? -->

This does not show operators, but it does give a useful representation of how process calls are connected and what are their inputs.

### Takeaway

You understand at a basic level how multi-step workflows are constructed using channels and operators and how they operate.

### What's next?

Learn how Nextflow pipelines can be modularized to promote code reuse and maintainability.

---

## 3. Modular code components

So far, all the workflows we've looked at have consisted of one single workflow file containing all the relevant code.

However, real-world pipelines typically benefit from being _modularized_, meaning that the code is split into different files.
This can make their development and maintenance more efficient and sustainable.

Here we are going to demonstrate the most common form of code modularity in Nextflow, which is the use of **modules**.

In Nextflow, a **module** is a single process definition that is encapsulated by itself in a standalone code file.
To use a module in a workflow, you just add a single-line import statement to your workflow code file; then you can integrate the process into the workflow the same way you normally would.

Putting processes into individual modules makes it possible to reuse process definitions in multiple workflows without producing multiple copies of the code.
This makes the code more shareable, flexible and maintainable.

We have of course once again prepared a suitable workflow for demonstration purposes, called `modular.nf`, along with a set of modules located in the `modules/` directory.

<!-- TODO: show directory contents -->

### 3.1. Examine the code

This time we're going to look at the code first, so let's open each of the files listed above.

We see that the processes and workflow logic are exactly the same as in the previous version of the workflow.
However, the process code is in the modules instead of being in the main workflow file, and there are now import statements in the workflow file telling Nextflow to pull them in at runtime.

```groovy title="hello-modules.nf" linenums="9" hl_lines="4"
// Include modules
include { sayHello } from './modules/sayHello.nf'
include { convertToUpper } from './modules/convertToUpper.nf'
include { collectGreetings } from './modules/collectGreetings.nf'

workflow {
```

You can look inside one of the modules to satisfy yourself that the process definition is unchanged; it's literally just been copy-pasted into a standalone file.

For example, this is the module containing the `sayHello` process:

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

So let's see what it looks like to run this new version.

### 3.2. Run the workflow

Run this command in your terminal, with the `-resume` flag:

```bash
nextflow run modular.nf --input greetings.csv -resume
```

Once again this should run successfully.

```console title="Output" linenums="1"
 N E X T F L O W   ~  version 25.04.3

Launching `modular.nf` [soggy_franklin] DSL2 - revision: bc8e1b2726

[j6/cdfa66] sayHello (1)       | 3 of 3, cached: ✔
[95/79484f] convertToUpper (2) | 3 of 3, cached: ✔
[5e/4358gc] collectGreetings   | 1 of 1, cached: ✔
```

You'll notice that the process executions all cached successfully, meaning that Nextflow recognized that it has already done the requested work, even though the code has been split up and the main workflow file has been renamed.

None of that matters to Nextflow; what matters is the job script that is generated once all the code has been pulled together and evaluated.

!!!note

    It is also possible to encapsulate a section of a workflow as a 'subworkflow' that can be imported into a larger pipeline, but that is outside the scope of this course.

    You can learn more about developing composable workflows in te Side Quest on [Workflows of Workflows](https://training.nextflow.io/latest/side_quests/workflows_of_workflows/).

### Takeaway

You know how processes can be stored in standalone modules to promote code reuse and improve maintainability.

### What's next?

Learn to use containers for managing software dependencies.

---

## 4. Using containerized software

So far the workflows we've been using as examples just needed to run very basic text procession operations using UNIX tools available in our environment.

However, real-world pipelines typically require specialized tools and packages that are not included by default in most environments.
Usually, you'd need to install these tools, manage their dependencies, and resolve any conflicts.

That is all very tedious and annoying.
A much better way to address this problem is to use **containers**.

A **container** is a lightweight, standalone, executable unit of software created from a container **image** that includes everything needed to run an application including code, system libraries and settings.

!!! note

    We teach this using the technology [Docker](https://www.docker.com/get-started/), but Nextflow supports [several other container technologies](https://www.nextflow.io/docs/latest/container.html#) as well.

### 4.1. Use a container directly

First, let's try interacting with a container directly.
This will help solidify your understanding of what containers are before we start using them in Nextflow.

#### 4.1.1. Pull the container image

To use a container, you usually download or "pull" a container image from a container registry, and then run the container image to create a container instance.

The general syntax is as follows:

```bash title="Syntax"
docker pull '<container>'
```

The `docker pull` part is the instruction to the container system to pull a container image from a repository.

The `'<container>'` part is the URI address of the container image.

As an example, let's pull a container image that contains [cowpy](https://github.com/jeffbuttars/cowpy), a python implementation of a tool called `cowsay` that generates ASCII art to display arbitrary text inputs in a fun way.

There are various repositories where you can find published containers.
We used the [Seqera Containers](https://seqera.io/containers/) service to generate this Docker container image from the `cowpy` Conda package: `'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'`.

Run the complete pull command:

```bash
docker pull 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
```

This gives you the following console output as the system downloads the image:

```console title="Output"
Unable to find image 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273' locally
131d6a1b707a8e65: Pulling from library/cowpy
dafa2b0c44d2: Pull complete
dec6b097362e: Pull complete
f88da01cff0b: Pull complete
4f4fb700ef54: Pull complete
92dc97a3ef36: Pull complete
403f74b0f85e: Pull complete
10b8c00c10a5: Pull complete
17dc7ea432cc: Pull complete
bb36d6c3110d: Pull complete
0ea1a16bbe82: Pull complete
030a47592a0a: Pull complete
622dd7f15040: Pull complete
895fb5d0f4df: Pull complete
Digest: sha256:fa50498b32534d83e0a89bb21fec0c47cc03933ac95c6b6587df82aaa9d68db3
Status: Downloaded newer image for community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273
community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273
```

Once the download is complete, you have a local copy of the container image.

You can also run a container interactively, which gives you a shell prompt inside the container and allows you to play with the command.

#### 4.1.2. Spin up the container

The general syntax is as follows:

```bash title="Syntax"
docker run --rm '<container>' [tool command]
```

The `docker run --rm '<container>'` part is the instruction to the container system to spin up a container instance from a container image and execute a command in it.
The `--rm` flag tells the system to shut down the container instance after the command has completed.

To run the container interactively, we add `-it` to the `docker run` command, and in this case we also add `-v` to mount a volume from the host system into the container using the following syntax:

```bash title="Syntax"
-v <outside_path>:<inside_path>
```

We need to do that because when you run a container, it is isolated from the host system by default.
This means that the container can't access any files on the host system unless you explicitly allow it to do so.

In our case `<outside_path>` will be the current working directory, so we can just use a dot (`.`), and `<inside_path>` is just a name we make up; let's call it `/data`.

Fully assembled, the container execution command looks like this:

```bash
docker run --rm -it -v .:/data 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
```

Run that command, and you should see your prompt change to something like `(base) root@b645838b3314:/tmp#`, which indicates that you are now inside the container.

You can verify this by running `ls` to list directory contents:

```bash
ls /
```

```console title="Output"
bin  boot  dev  data  etc  home  lib  lib64  media  mnt  opt  proc  root  run  sbin  srv  sys  tmp  usr  var
```

You can see that the filesystem inside the container is different from the filesystem on your host system.

Conveniently, the command we ran mounted the current working directory as a volume that is accessible under `/data` inside the container.

You can check that it works by listing the contents of `/data`:

```bash
ls /data
```

Depending on what part of this training you've done before, the output below may look slightly different, but you should now see the contents of the `data` directory.

<!-- TODO: show output: ls output needs to be updated
```console title="Output"

```
-->

This effectively established a tunnel through the container wall that you can use to access that part of your filesystem.

#### 4.1.3. Run the `cowpy` tool

Now that you are inside the container, you can run the `cowpy` command directly and give it some parameters.
For example, the tool documentation says we can set the character ('cowacter') with `-c`.

```bash
cowpy "Hello Containers" -c tux
```

Now the output shows the Linux penguin, Tux, because we specified the `-c tux` parameter.

```console title="Output"
 __________________
< Hello Containers >
 ------------------
   \
    \
        .--.
       |o_o |
       |:_/ |
      //   \ \
     (|     | )
    /'\_   _/`\
    \___)=(___/
```

Because you're inside the container, you can run the cowpy command as many times as you like, varying the input parameters, without having to bother with Docker commands.

!!! Tip

    Use the '-c' flag to pick a different character, including:
    `beavis`, `cheese`, `daemon`, `dragonandcow`, `ghostbusters`, `kitty`, `moose`, `milk`, `stegosaurus`, `turkey`, `turtle`, `tux`

This is neat. What would be even neater is if we could feed our `greetings.csv` as input into this.

#### 4.1.4. Run the `cowpy` tool

Good news: we can, since we have access to our files via the `/data` volume mount!

We can use `cat /data/greetings.csv | ` to pipe the contents of the CSV file into the `cowpy` command.

```bash
cat /data/greetings.csv | cowpy -c turkey
```

This produces the desired ASCII art of a turkey rattling off our example greetings:

```console title="Output"
 _________
/ Hello   \
| Bonjour |
\ Holà    /
 ---------
  \                                  ,+*^^*+___+++_
   \                           ,*^^^^              )
    \                       _+*                     ^**+_
     \                    +^       _ _++*+_+++_,         )
              _+^^*+_    (     ,+*^ ^          \+_        )
             {       )  (    ,(    ,_+--+--,      ^)      ^\
            { (\@)    } f   ,(  ,+-^ __*_*_  ^^\_   ^\       )
           {:;-/    (_+*-+^^^^^+*+*<_ _++_)_    )    )      /
          ( /  (    (        ,___    ^*+_+* )   <    <      \
           U _/     )    *--<  ) ^\-----++__)   )    )       )
            (      )  _(^)^^))  )  )\^^^^^))^*+/    /       /
          (      /  (_))_^)) )  )  ))^^^^^))^^^)__/     +^^
         (     ,/    (^))^))  )  ) ))^^^^^^^))^^)       _)
          *+__+*       (_))^)  ) ) ))^^^^^^))^^^^^)____*^
          \             \_)^)_)) ))^^^^^^^^^^))^^^^)
           (_             ^\__^^^^^^^^^^^^))^^^^^^^)
             ^\___            ^\__^^^^^^))^^^^^^^^)\\
                  ^^^^^\uuu/^^\uuu/^^^^\^\^\^\^\^\^\^\
                     ___) >____) >___   ^\_\_\_\_\_\_\)
                    ^^^//\\_^^//\\_^       ^(\_\_\_\)
                      ^^^ ^^ ^^^ ^
```

Feel free to play around with this command.
When you're done, exit the container using the `exit` command:

```bash
exit
```

You will find yourself back in your normal shell.

### 4.2. Use a container in a workflow

TODO: clone the content from hello_containers.md

### Takeaway

You understand what role containers play in managing software tool versions and ensuring reproducibility.

More generally, you have a basic understanding of the most common and most important components of real-world Nexflow pipelines.

### What's next?

Take another break!
TODO: finalize the transition text
