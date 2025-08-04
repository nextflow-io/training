# Part 5: Hello Containers

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/5PyOWjKnNmg?si=QinuAnFwFj-Z8CrO&amp;list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } See [the whole playlist](https://www.youtube.com/playlist?list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik) on the Nextflow YouTube channel.

:green_book: The video transcript is available [here](./transcripts/05_hello_containers.md).
///

In Parts 1-4 of this training course, you learned how to use the basic building blocks of Nextflow to assemble a simple workflow capable of processing some text, parallelizing execution if there were multiple inputs, and collecting the results for further processing.

However, you were limited to basic UNIX tools available in your environment.
Real-world tasks often require various tools and packages not included by default.
Typically, you'd need to install these tools, manage their dependencies, and resolve any conflicts.

That is all very tedious and annoying, so we're going to show you how to use **containers** to solve this problem much more conveniently.

A **container** is a lightweight, standalone, executable unit of software created from a container **image** that includes everything needed to run an application including code, system libraries and settings.

!!! note

    We'll be teaching this using the technology [Docker](https://www.docker.com/get-started/), but Nextflow supports [several other container technologies](https://www.nextflow.io/docs/latest/container.html#) as well.

---

## 0. Warmup: Run `hello-containers.nf`

We're going to use the workflow script `hello-containers.nf` as a starting point for the second section.
It is equivalent to the script produced by working through Part 4 of this training course.

Just to make sure everything is working, run the script once before making any changes:

```bash
nextflow run hello-containers.nf
```

This should produce the following output:

```console title="Output"
 N E X T F L O W   ~  version 25.04.3

Launching `hello-containers.nf` [tender_becquerel] DSL2 - revision: f7cat8e223

executor >  local (7)
[bd/4bb541] sayHello (1)         [100%] 3 of 3 ✔
[85/b627e8] convertToUpper (3)   [100%] 3 of 3 ✔
[7d/f7961c] collectGreetings     [100%] 1 of 1 ✔
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

If that worked for you, you're ready to learn how to use containers.

---

## 1. Use a container 'manually'

What we want to do is add a step to our workflow that will use a container for execution.

However, we are first going to go over some basic concepts and operations to solidify your understanding of what containers are before we start using them in Nextflow.

### 1.1. Pull the container image

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

### 1.2. Use the container to run `cowpy` as a one-off command

One very common way that people use containers is to run them directly, _i.e._ non-interactively.
This is great for running one-off commands.

The general syntax is as follows:

```bash title="Syntax"
docker run --rm '<container>' [tool command]
```

The `docker run --rm '<container>'` part is the instruction to the container system to spin up a container instance from a container image and execute a command in it.
The `--rm` flag tells the system to shut down the container instance after the command has completed.

The `[tool command]` syntax depends on the tool you are using and how the container is set up.
Let's just start with `cowpy`.

Fully assembled, the container execution command looks like this:

```bash
docker run --rm 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273' cowpy
```

Run it to produce the following output:

```console title="Output"
 ______________________________________________________
< Cowacter, eyes:default, tongue:False, thoughts:False >
 ------------------------------------------------------
     \   ^__^
      \  (oo)\_______
         (__)\       )\/\
           ||----w |
           ||     ||
```

The system spun up the container, ran the `cowpy` command with its parameters, sent the output to the console and finally, shut down the container instance.

### 1.3. Use the container to run `cowpy` interactively

You can also run a container interactively, which gives you a shell prompt inside the container and allows you to play with the command.

#### 1.3.1. Spin up the container

To run interactively, we just add `-it` to the `docker run` command.
Optionally, we can specify the shell we want to use inside the container by appending _e.g._ `/bin/bash` to the command.

```bash
docker run --rm -it 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273' /bin/bash
```

Notice that your prompt changes to something like `(base) root@b645838b3314:/tmp#`, which indicates that you are now inside the container.

You can verify this by running `ls` to list directory contents:

```bash
ls /
```

```console title="Output"
bin  boot  dev  etc  home  lib  lib64  media  mnt  opt  proc  root  run  sbin  srv  sys  tmp  usr  var
```

You can see that the filesystem inside the container is different from the filesystem on your host system.

!!! note

    When you run a container, it is isolated from the host system by default.
    This means that the container can't access any files on the host system unless you explicitly allow it to do so.

    You will learn how to do that in a minute.

#### 1.3.2. Run the desired tool command(s)

Now that you are inside the container, you can run the `cowpy` command directly and give it some parameters.
For example, the tool documentation says we can change the character ('cowacter') with `-c`.

```bash
cowpy "Hello Containers" -c tux
```

Now the output shows the Linux penguin, Tux, instead of the default cow, because we specified the `-c tux` parameter.

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
But since we don't have access to the filesystem, we can't.

Let's fix that.

#### 1.3.3. Exit the container

To exit the container, you can type `exit` at the prompt or use the ++ctrl+d++ keyboard shortcut.

```bash
exit
```

Your prompt should now be back to what it was before you started the container.

#### 1.3.4. Mount data into the container

When you run a container, it is isolated from the host system by default.
This means that the container can't access any files on the host system unless you explicitly allow it to do so.

One way to do this is to **mount** a **volume** from the host system into the container using the following syntax:

```bash title="Syntax"
-v <outside_path>:<inside_path>
```

In our case `<outside_path>` will be the current working directory, so we can just use a dot (`.`), and `<inside_path>` is just a name we make up; let's call it `/data`.

To mount a volume, we replace the paths and add the volume mounting argument to the docker run command as follows:

```bash
docker run --rm -it -v .:/data 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273' /bin/bash
```

This mounts the current working directory as a volume that will be accessible under `/data` inside the container.

You can check that it works by listing the contents of `/data`:

```bash
ls /data
```

Depending on what part of this training you've done before, the output below my look slightly different.

```console title="Output"
greetings.csv      hello-config.nf      hello-modules.nf   hello-world.nf  nextflow.config  solutions         work
hello-channels.nf  hello-containers.nf  hello-workflow.nf  modules         results          test-params.json
```

<!-- ls output may need to be updated -->

You can now see the contents of the `data` directory from inside the container, including the `greetings.csv` file.

This effectively established a tunnel through the container wall that you can use to access that part of your filesystem.

#### 1.3.5. Use the mounted data

Now that we have mounted the `data` directory into the container, we can use the `cowpy` command to display the contents of the `greetings.csv` file.

To do this, we'll use `cat /data/greetings.csv | ` to pipe the contents of the CSV file into the `cowpy` command.

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
When you're done, exit the container as previously:

```bash
exit
```

You will find yourself back in your normal shell.

### Takeaway

You know how to pull a container and run it either as a one-off or interactively. You also know how to make your data accessible from within your container, which lets you try any tool you're interested in on real data without having to install any software on your system.

### What's next?

Learn how to use containers for the execution of Nextflow processes.

---

## 2. Use containers in Nextflow

Nextflow has built-in support for running processes inside containers to let you run tools you don't have installed in your compute environment.
This means that you can use any container image you like to run your processes, and Nextflow will take care of pulling the image, mounting the data, and running the process inside it.

To demonstrate this, we are going to add a `cowpy` step to the pipeline we've been developing, after the `collectGreetings` step.

### 2.1. Write a `cowpy` module

#### 2.1.1. Create a file stub for the new module

Create an empty file for the module called `cowpy.nf`.

```bash
touch modules/cowpy.nf
```

This gives us a place to put the process code.

#### 2.1.2. Copy the `cowpy` process code in the module file

We can model our `cowpy` process on the other processes we've written previously.

```groovy title="modules/cowpy.nf" linenums="1"
#!/usr/bin/env nextflow

// Generate ASCII art with cowpy
process cowpy {

    publishDir 'results', mode: 'copy'

    input:
        path input_file
        val character

    output:
        path "cowpy-${input_file}"

    script:
    """
    cat $input_file | cowpy -c "$character" > cowpy-${input_file}
    """

}
```

The output will be a new text file containing the ASCII art generated by the `cowpy` tool.

### 2.2. Add cowpy to the workflow

Now we need to import the module and call the process.

#### 2.2.1. Import the `cowpy` process into `hello-containers.nf`

Insert the import declaration above the workflow block and fill it out appropriately.

=== "After"

    ```groovy title="hello-containers.nf" linenums="9" hl_lines="5"
    // Include modules
    include { sayHello } from './modules/sayHello.nf'
    include { convertToUpper } from './modules/convertToUpper.nf'
    include { collectGreetings } from './modules/collectGreetings.nf'
    include { cowpy } from './modules/cowpy.nf'

    workflow {
    ```

=== "Before"

    ```groovy title="hello-containers.nf" linenums="9"
    // Include modules
    include { sayHello } from './modules/sayHello.nf'
    include { convertToUpper } from './modules/convertToUpper.nf'
    include { collectGreetings } from './modules/collectGreetings.nf'

    workflow {
    ```

#### 2.2.2. Add a call to the `cowpy` process in the workflow

Let's connect the `cowpy()` process to the output of the `collectGreetings()` process, which as you may recall produces two outputs:

- `collectGreetings.out.outfile` contains the output file
- `collectGreetings.out.count` contains the count of greetings per batch

In the workflow block, make the following code change:

=== "After"

    ```groovy title="hello-containers.nf" linenums="28" hl_lines="7 8"
        // collect all the greetings into one file
        collectGreetings(convertToUpper.out.collect(), params.batch)

        // emit a message about the size of the batch
        collectGreetings.out.count.view{ num_greetings -> "There were $num_greetings greetings in this batch" }

        // generate ASCII art of the greetings with cowpy
        cowpy(collectGreetings.out.outfile, params.character)
    ```

=== "Before"

    ```groovy title="hello-containers.nf" linenums="28"
        // collect all the greetings into one file
        collectGreetings(convertToUpper.out.collect(), params.batch)

        // emit a message about the size of the batch
        collectGreetings.out.count.view{ num_greetings -> "There were $num_greetings greetings in this batch" }
    ```

Notice that we include a new CLI parameter, `params.character`, in order to specify which character we want to have say the greetings.

#### 2.2.3. Set a default value for `params.character`

We like to be lazy and skip typing parameters in our command lines.

=== "After"

    ```groovy title="hello-containers.nf" linenums="3" hl_lines="6"
    /*
     * Pipeline parameters
     */
    params.greeting = 'greetings.csv'
    params.batch = 'test-batch'
    params.character = 'turkey'
    ```

=== "Before"

    ```groovy title="hello-containers.nf" linenums="3"
    /*
     * Pipeline parameters
     */
    params.greeting = 'greetings.csv'
    params.batch = 'test-batch'
    ```

That should be all we need to make this work.

#### 2.2.4. Run the workflow to verify that it works

Run this with the `-resume` flag.

```bash
nextflow run hello-containers.nf -resume
```

Oh no, there's an error!

```console title="Output"
 N E X T F L O W   ~  version 25.04.3

Launching `hello-containers.nf` [special_lovelace] DSL2 - revision: 028a841db1

executor >  local (1)
[f6/cc0107] sayHello (1)       | 3 of 3, cached: 3 ✔
[2c/67a06b] convertToUpper (3) | 3 of 3, cached: 3 ✔
[1a/bc5901] collectGreetings   | 1 of 1, cached: 1 ✔
[b2/488871] cowpy             | 0 of 1
There were 3 greetings in this batch
ERROR ~ Error executing process > 'cowpy'

Caused by:
  Process `cowpy` terminated with an error exit status (127)

Command executed:

  cat COLLECTED-test-batch-output.txt | cowpy -c "turkey" > cowpy-COLLECTED-test-batch-output.txt

Command exit status:
  127

Command output:
  (empty)

Command error:
  .command.sh: line 2: cowpy: command not found

(trimmed output)
```

This error code, `error exit status (127)` means the executable we asked for was not found.

Of course, since we're calling the `cowpy` tool but we haven't actually specified a container yet.

### 2.3. Use a container to run it

We need to specify a container and tell Nextflow to use it for the `cowpy()` process.

#### 2.3.1. Specify a container for the `cowpy` process to use

Edit the `cowpy.nf` module to add the `container` directive to the process definition as follows:

=== "After"

    ```groovy title="modules/cowpy.nf" linenums="4" hl_lines="4"
    process cowpy {

        publishDir 'containers/results', mode: 'copy'
        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
    ```

=== "Before"

    ```groovy title="modules/cowpy.nf" linenums="4"
    process cowpy {

        publishDir 'containers/results', mode: 'copy'
    ```

This tells Nextflow that if the use of Docker is enabled, it should use the container image specified here to execute the process.

#### 2.3.2. Enable use of Docker via the `nextflow.config` file

Here we are going to slightly anticipate the topic of the next and last part of this course (Part 6), which covers configuration.

One of the main ways Nextflow offers for configuring workflow execution is to use a `nextflow.config` file. When such a file is present in the current directory, Nextflow will automatically load it in and apply any configuration it contains.

We provided a `nextflow.config` file with a single line of code that disables Docker: `docker.enabled = false`.

Now, let's switch that to `true` to enable Docker:

=== "After"

    ```console title="nextflow.config" linenums="1" hl_lines="1"
    docker.enabled = true
    ```

=== "Before"

    ```console title="nextflow.config" linenums="1"
    docker.enabled = false
    ```

!!! note

    It is possible to enable Docker execution from the command-line, on a per-run basis, using the `-with-docker <container>` parameter.
    However, that only allows us to specify one container for the entire workflow, whereas the approach we just showed you allows us to specify a different container per process.
    This is better for modularity, code maintenance and reproducibility.

#### 2.3.3. Run the workflow with Docker enabled

Run the workflow with the `-resume` flag:

```bash
nextflow run hello-containers.nf -resume
```

This time it does indeed work.

```console title="Output" linenums="1"
 N E X T F L O W   ~  version 25.04.3

Launching `hello-containers.nf` [elegant_brattain] DSL2 - revision: 028a841db1

executor >  local (1)
[95/fa0bac] sayHello (3)       | 3 of 3, cached: 3 ✔
[92/32533f] convertToUpper (3) | 3 of 3, cached: 3 ✔
[aa/e697a2] collectGreetings   | 1 of 1, cached: 1 ✔
[7f/caf718] cowpy              | 1 of 1 ✔
There were 3 greetings in this batch
```

You can find the cowpy'ed output in the `results` directory.

```console title="results/cowpy-COLLECTED-test-batch-output.txt"
 _________
/ HOLà    \
| HELLO   |
\ BONJOUR /
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

You see that the character is saying all the greetings, just as it did when we ran the `cowpy` command on the `greetings.csv` file from inside the container.

<!-- considering a side quest where we show how to use a conditional to skip the collect step if we want to emit the cowpy'ed greetings individually, and how to use metadata management to assign a specific character to each greeting, maybe do some cross products etc -->

#### 2.3.4. Inspect how Nextflow launched the containerized task

Let's take a look at the work subdirectory for one of the `cowpy` process calls to get a bit more insight on how Nextflow works with containers under the hood.

Check the output from your `nextflow run` command to find the path to the work subdirectory for the `cowpy` process.
Looking at what we got for the run shown above, the console log line for the `cowpy` process starts with `[7f/caf718]`.
That corresponds to the following truncated directory path: `work/7f/caf718`.
In it, you will find the `.command.run` file that contains all the commands Nextflow ran on your behalf in the course of executing the pipeline.

Open the `.command.run` file and search for `nxf_launch`; you should see something like this:

```bash
nxf_launch() {
    docker run -i --cpu-shares 1024 -e "NXF_TASK_WORKDIR" -v /workspaces/training/hello-nextflow/work:/workspaces/training/hello-nextflow/work -w "$NXF_TASK_WORKDIR" --name $NXF_BOXID community.wave.seqera.io/library/pip_cowpy:131d6a1b707a8e65 /bin/bash -ue /workspaces/training/hello-nextflow/work/7f/caf7189fca6c56ba627b75749edcb3/.command.sh
}
```

As you can see, Nextflow is using the `docker run` command to launch the process call.
It also mounts the corresponding work subdirectory into the container, sets the working directory inside the container accordingly, and runs our templated bash script in the `.command.sh` file.

All the hard work we had to do manually in the previous section is done for us by Nextflow!

### Takeaway

You know how to use containers in Nextflow to run processes.

### What's next?

Take a break!
When you're ready, move on to Part 6 to learn how to configure the execution of your pipeline to fit your infrastructure as well as manage configuration of inputs and parameters.
It's the very last part and then you're done!
