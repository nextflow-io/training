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
As you might imagine, that is a going to be very helpful for making your pipelines more reproducible.

Note that we'll be teaching this using [Docker](https://www.docker.com/get-started/), but keep in mind Nextflow supports [several other container technologies](https://www.nextflow.io/docs/latest/container.html#) as well.

??? info "How to begin from this section"

    This section of the course assumes you have completed Parts 1-4 of the [Hello Nextflow](./index.md) course and have a complete working pipeline.

    If you're starting the course from this point, you'll need to copy the `modules` directory over from the solutions:

    ```bash
    cp -r solutions/4-hello-modules/modules .
    ```

---

## 0. Warmup: Run `hello-containers.nf`

We're going to use the workflow script `hello-containers.nf` as a starting point.
It is equivalent to the script produced by working through Part 4 of this training course, except we've changed the output destinations:

```groovy title="hello-containers.nf" linenums="37" hl_lines="3 7 11 15"
output {
    first_output {
        path 'hello_containers'
        mode 'copy'
    }
    uppercased {
        path 'hello_containers'
        mode 'copy'
    }
    collected {
        path 'hello_containers'
        mode 'copy'
    }
    batch_report {
        path 'hello_containers'
        mode 'copy'
    }
}
```

Just to make sure everything is working, run the script once before making any changes:

```bash
nextflow run hello-containers.nf
```

??? success "Command output"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-containers.nf` [nice_escher] DSL2 - revision: d5dfdc9872

    executor > local (7)
    [5a/ec1fa1] sayHello (2) [100%] 3 of 3 ✔
    [30/32b5b8] convertToUpper (3) [100%] 3 of 3 ✔
    [d3/be01bc] collectGreetings [100%] 1 of 1 ✔

    ```

As previously, you will find the output files in the directory specified in the `output` block (`results/hello_containers/`).

??? abstract "Directory contents"

    ```console
    results/hello_containers/
    ├── Bonjour-output.txt
    ├── COLLECTED-batch-output.txt
    ├── Hello-output.txt
    ├── Holà-output.txt
    ├── batch-report.txt
    ├── UPPER-Bonjour-output.txt
    ├── UPPER-Hello-output.txt
    └── UPPER-Holà-output.txt
    ```

If that worked for you, you're ready to learn how to use containers.

---

## 1. Use a container 'manually'

What we want to do is add a step to our workflow that will use a container for execution.

However, we are first going to go over some basic concepts and operations to solidify your understanding of what containers are before we start using them in Nextflow.

### 1.1. Pull the container image

To use a container, you usually download or _pull_ a container image from a container registry, and then run the container image to create a container instance.

The general syntax is as follows:

```bash title="Syntax"
docker pull '<container>'
```

The `docker pull` part is the instruction to the container system to pull a container image from a repository.

The `'<container>'` part is the URI address of the container image.

As an example, let's pull a container image that contains [cowpy](https://github.com/jeffbuttars/cowpy), a python implementation of a tool called `cowsay` that generates ASCII art to display arbitrary text inputs in a fun way.

```txt title="Example"
 ________________________
< Are we having fun yet? >
 ------------------------
    \                                  ___-------___
     \                             _-~~             ~~-_
      \                         _-~                    /~-_
             /^\__/^\         /~  \                   /    \
           /|  O|| O|        /      \_______________/        \
          | |___||__|      /       /                \          \
          |          \    /      /                    \          \
          |   (_______) /______/                        \_________ \
          |         / /         \                      /            \
           \         \^\\         \                  /               \     /
             \         ||           \______________/      _-_       //\__//
               \       ||------_-~~-_ ------------- \ --/~   ~\    || __/
                 ~-----||====/~     |==================|       |/~~~~~
                  (_(__/  ./     /                    \_\      \.
                         (_(___/                         \_____)_)
```

There are various repositories where you can find published containers.
We used the [Seqera Containers](https://seqera.io/containers/) service to generate this Docker container image from the `cowpy` Conda package: `'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'`.

Run the complete pull command:

```bash
docker pull 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
```

??? success "Command output"

    ```console
    1.1.5--3db457ae1977a273: Pulling from library/cowpy
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
    c23bdb422167: Pull complete
    e1686ff32a11: Pull complete
    Digest: sha256:1ebc0043e8cafa61203bf42d29fd05bd14e7b4298e5e8cf986504c15f5aa4160
    Status: Downloaded newer image for community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273
    community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273
    ```

If you've never downloaded the image before, this may take a minute to complete.
Once it's done, you have a local copy of the container image.

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

Fully assembled, the container execution command looks like this; go ahead and run it.

```bash
docker run --rm 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273' cowpy
```

??? success "Command output"

    ```console
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

You can verify this by running `ls /` to list directory contents from the root of the filesystem:

```bash
ls /
```

??? abstract "Command output"

    ```console
    bin  boot  dev  etc  home  lib  lib64  media  mnt  opt  proc  root  run  sbin  srv  sys  tmp  usr  var
    ```

We use `ls` here instead of `tree` because the `tree` utility is not available in this container.
You can see that the filesystem inside the container is different from the filesystem on your host system.

One limitation of what we just did is that the container is completely isolated from the host system by default.
This means that the container can't access any files on the host system unless you explicitly allow it to do so.

We'll show you how to do that in a minute.

#### 1.3.2. Run the desired tool command(s)

Now that you are inside the container, you can run the `cowpy` command directly and give it some parameters.
For example, the tool documentation says we can change the character ('cowacter') with `-c`.

```bash
cowpy "Hello Containers" -c tux
```

??? success "Command output"

    ```console
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

Now the output shows the Linux penguin, Tux, instead of the default cow, because we specified the `-c tux` parameter.

Because you're inside the container, you can run the `cowpy` command as many times as you like, varying the input parameters, without having to bother with Docker commands.

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

As noted earlier, the container is isolated from the host system by default.

To allow the container to access the host filesystem, you can **mount** a **volume** from the host system into the container using the following syntax:

```bash title="Syntax"
-v <outside_path>:<inside_path>
```

In our case `<outside_path>` will be the current working directory, so we can just use a dot (`.`), and `<inside_path>` is just an alias we make up; let's call it `/my_project` (the inside path must be absolute).

To mount a volume, we replace the paths and add the volume mounting argument to the docker run command as follows:

```bash
docker run --rm -it -v .:/my_project 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273' /bin/bash
```

This mounts the current working directory as a volume that will be accessible under `/my_project` inside the container.

You can check that it works by listing the contents of `/my_project`:

```bash
ls /my_project
```

??? success "Command output"

    ```console
    data               hello-config.nf      hello-modules.nf   hello-world.nf  nextflow.config  solutions         work
    hello-channels.nf  hello-containers.nf  hello-workflow.nf  modules         results          test-params.json
    ```

You can now see the contents of the working directory from inside the container, including the `greetings.csv` file under `data/`.

This effectively established a tunnel through the container wall that you can use to access that part of your filesystem.

#### 1.3.5. Use the mounted data

Now that we have mounted the working directory into the container, we can use the `cowpy` command to display the contents of the `greetings.csv` file.

To do this, we'll use `cat /my_project/data/greetings.csv | ` to pipe the contents of the CSV file into the `cowpy` command.

```bash
cat /my_project/data/greetings.csv | cowpy -c turkey
```

??? success "Command output"

    ```console title="data/greetings.csv"
     ____________________
    / Hello,English,123  \
    | Bonjour,French,456 |
    \ Holà,Spanish,789   /
    --------------------
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

This produces the desired ASCII art of a turkey rattling off our example greetings!
Except here the turkey is repeating the full rows instead of just the greetings.
We already know our Nextflow workflow will do a better job!

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

<figure class="excalidraw">
--8<-- "docs/nextflow_run/img/hello-pipeline-cowpy.svg"
</figure>

Moo if you're ready to dive in!

### 2.1. Write a `cowpy` module

First, let's create the `cowpy` process module.

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

    input:
    path input_file
    val character

    output:
    path "cowpy-${input_file}"

    script:
    """
    cat ${input_file} | cowpy -c "${character}" > cowpy-${input_file}
    """

}
```

The process expects an `input_file` containing the greetings as well as a `character` value.

The output will be a new text file containing the ASCII art generated by the `cowpy` tool.

### 2.2. Add cowpy to the workflow

Now we need to import the module and call the process.

#### 2.2.1. Import the `cowpy` process into `hello-containers.nf`

Insert the import declaration above the workflow block and fill it out appropriately.

=== "After"

    ```groovy title="hello-containers.nf" linenums="3" hl_lines="5"
    // Include modules
    include { sayHello } from './modules/sayHello.nf'
    include { convertToUpper } from './modules/convertToUpper.nf'
    include { collectGreetings } from './modules/collectGreetings.nf'
    include { cowpy } from './modules/cowpy.nf'
    ```

=== "Before"

    ```groovy title="hello-containers.nf" linenums="3"
    // Include modules
    include { sayHello } from './modules/sayHello.nf'
    include { convertToUpper } from './modules/convertToUpper.nf'
    include { collectGreetings } from './modules/collectGreetings.nf'
    ```

Now the `cowpy` module is available to use in the workflow.

#### 2.2.2. Add a call to the `cowpy` process in the workflow

Let's connect the `cowpy()` process to the output of the `collectGreetings()` process, which as you may recall produces two outputs:

- `collectGreetings.out.outfile` contains the output file <--_what we want_
- `collectGreetings.out.report` contains the report file with the count of greetings per batch

In the workflow block, make the following code change:

=== "After"

    ```groovy title="hello-containers.nf" linenums="19" hl_lines="12-13"
        main:
        // create a channel for inputs from a CSV file
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        // emit a greeting
        sayHello(greeting_ch)
        // convert the greeting to uppercase
        convertToUpper(sayHello.out)
        // collect all the greetings into one file
        collectGreetings(convertToUpper.out.collect(), params.batch)
        // generate ASCII art of the greetings with cowpy
        cowpy(collectGreetings.out.outfile, params.character)
    ```

=== "Before"

    ```groovy title="hello-containers.nf" linenums="19"
        main:
        // create a channel for inputs from a CSV file
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        // emit a greeting
        sayHello(greeting_ch)
        // convert the greeting to uppercase
        convertToUpper(sayHello.out)
        // collect all the greetings into one file
        collectGreetings(convertToUpper.out.collect(), params.batch)
    ```

Notice that we declared a new CLI parameter, `params.character`, in order to specify which character we want to have say the greetings.

#### 2.2.3. Add the `character` parameter to the `params` block

This is technically optional but it's the recommended practice and it's an opportunity to set a default value for the character while we're at it.

=== "After"

    ```groovy title="hello-containers.nf" linenums="9" hl_lines="7"
    /*
    * Pipeline parameters
    */
    params {
        input: Path = 'data/greetings.csv'
        batch: String = 'batch'
        character: String = 'turkey'
    }
    ```

=== "Before"

    ```groovy title="hello-containers.nf" linenums="9"
    /*
    * Pipeline parameters
    */
    params {
        input: Path = 'data/greetings.csv'
        batch: String = 'batch'
    }
    ```

Now we can be lazy and skip typing the character parameter in our command lines.

#### 2.2.4. Update the workflow outputs

We need to update the workflow outputs to publish the output of the `cowpy` process.

##### 2.2.4.1. Update the `publish:` section

In the `workflow block`, make the following code change:

=== "After"

    ```groovy title="hello-containers.nf" linenums="34" hl_lines="6"
        publish:
        first_output = sayHello.out
        uppercased = convertToUpper.out
        collected = collectGreetings.out.outfile
        batch_report = collectGreetings.out.report
        cowpy_art = cowpy.out
    ```

=== "Before"

    ```groovy title="hello-containers.nf" linenums="34"
        publish:
        first_output = sayHello.out
        uppercased = convertToUpper.out
        collected = collectGreetings.out.outfile
        batch_report = collectGreetings.out.report
    ```

The `cowpy` process only produces one output so we can refer to it the usual way by appending `.out`.

But for now, let's finish updating the workflow-level outputs.

##### 2.2.4.2. Update the `output` block

We need to add the final `cowpy_art` output to the `output` block. While we're at it, let's also edit the publishing destinations since now our pipeline is complete and we know what outputs we really care about.

In the `output` block, make the following code changes:

=== "After"

    ```groovy title="hello-containers.nf" linenums="42" hl_lines="3 7 11 15 18-21"
    output {
        first_output {
            path 'hello_containers/intermediates'
            mode 'copy'
        }
        uppercased {
            path 'hello_containers/intermediates'
            mode 'copy'
        }
        collected {
            path 'hello_containers/intermediates'
            mode 'copy'
        }
        batch_report {
            path 'hello_containers'
            mode 'copy'
        }
        cowpy_art {
            path 'hello_containers'
            mode 'copy'
        }
    }
    ```

=== "Before"

    ```groovy title="hello-containers.nf" linenums="42" hl_lines="3 7 11 15"
    output {
        first_output {
            path 'hello_containers'
            mode 'copy'
        }
        uppercased {
            path 'hello_containers'
            mode 'copy'
        }
        collected {
            path 'hello_containers'
            mode 'copy'
        }
        batch_report {
            path 'hello_containers'
            mode 'copy'
        }
    }
    ```

Now the published outputs will be a bit more organized.

#### 2.2.5. Run the workflow

Just to recap, this is what we are aiming for:

<figure class="excalidraw">
--8<-- "docs/hello_nextflow/img/hello_pipeline_complete.svg"
</figure>

Do you think it's going to work?

Let's delete the previous published outputs to have a clean slate, and run the workflow with the `-resume` flag.

```bash
rm -r hello_containers/
nextflow run hello-containers.nf -resume
```

??? failure "Command output (edited for clarity)"

    ```console hl_lines="10 13 20-21 26-27"
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-containers.nf` [lonely_woese] DSL2 - revision: abf1dccf7f

    executor >  local (1)
    [c9/f5c686] sayHello (3)       [100%] 3 of 3, cached: 3 ✔
    [ef/3135a8] convertToUpper (3) [100%] 3 of 3, cached: 3 ✔
    [7f/f435e3] collectGreetings   [100%] 1 of 1, cached: 1 ✔
    [9b/02e776] cowpy              [  0%] 0 of 1 ✘
    ERROR ~ Error executing process > 'cowpy'

    Caused by:
      Process `cowpy` terminated with an error exit status (127)


    Command executed:

      cat COLLECTED-batch-output.txt | cowpy -c "turkey" > cowpy-COLLECTED-batch-output.txt

    Command exit status:
      127

    Command output:
      (empty)

    Command error:
      .command.sh: line 2: cowpy: command not found

    Work dir:
      /workspaces/training/hello-nextflow/work/9b/02e7761db848f82db3c3e59ff3a9b6

    Tip: when you have fixed the problem you can continue the execution adding the option `-resume` to the run command line

    -- Check '.nextflow.log' file for details
    ERROR ~ Cannot access first() element from an empty List

    -- Check '.nextflow.log' file for details
    ```

Oh no, there's an error!
The error code given by `error exit status (127)` means the executable we asked for was not found.

That makes sense, since we're calling the `cowpy` tool but we haven't actually specified a container yet (oops).

### 2.3. Use a container to run the `cowpy` process

We need to specify a container and tell Nextflow to use it for the `cowpy()` process.

#### 2.3.1. Specify a container for `cowpy`

We can use the same image we were using directly in the first section of this tutorial.

Edit the `cowpy.nf` module to add the `container` directive to the process definition as follows:

=== "After"

    ```groovy title="modules/cowpy.nf" linenums="4" hl_lines="3"
    process cowpy {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

        input:
        path input_file
        val character

        output:
        path "cowpy-${input_file}"

        script:
        """
        cat ${input_file} | cowpy -c "${character}" > cowpy-${input_file}
        """
    }
    ```

=== "Before"

    ```groovy title="modules/cowpy.nf" linenums="4"
    process cowpy {

        input:
        path input_file
        val character

        output:
        path "cowpy-${input_file}"

        script:
        """
        cat ${input_file} | cowpy -c "${character}" > cowpy-${input_file}
        """
    }
    ```

This tells Nextflow that _if the use of Docker is enabled_, it should use the container image specified here to execute the process.

#### 2.3.2. Enable use of Docker via the `nextflow.config` file

Notice we said _'if the use of Docker is enabled'_. By default, it is not, so we need to tell Nextflow it's allowed to use Docker.
To that end, we are going to slightly anticipate the topic of the next and last part of this course (Part 6), which covers configuration.

One of the main ways Nextflow offers for configuring workflow execution is to use a `nextflow.config` file.
When such a file is present in the current directory, Nextflow will automatically load it in and apply any configuration it contains.

We provided a `nextflow.config` file with a single line of code that explicitly disables Docker: `docker.enabled = false`.

Now, let's switch that to `true` to enable Docker:

=== "After"

    ```console title="nextflow.config" linenums="1" hl_lines="1"
    docker.enabled = true
    ```

=== "Before"

    ```console title="nextflow.config" linenums="1" hl_lines="1"
    docker.enabled = false
    ```

!!! tip

    It is possible to enable Docker execution from the command-line, on a per-run basis, using the `-with-docker <container>` parameter.
    However, that only allows us to specify one container for the entire workflow, whereas the approach we just showed you allows us to specify a different container per process.
    This is better for modularity, code maintenance and reproducibility.

#### 2.3.3. Run the workflow with Docker enabled

Run the workflow with the `-resume` flag:

```bash
nextflow run hello-containers.nf -resume
```

??? success "Command output"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-containers.nf` [drunk_perlman] DSL2 - revision: abf1dccf7f

    executor >  local (1)
    [c9/f5c686] sayHello (3)       [100%] 3 of 3, cached: 3 ✔
    [ef/3135a8] convertToUpper (3) [100%] 3 of 3, cached: 3 ✔
    [7f/f435e3] collectGreetings   [100%] 1 of 1, cached: 1 ✔
    [98/656c6c] cowpy              [100%] 1 of 1 ✔
    ```

This time it does indeed work!
As usual you can find the workflow outputs in the corresponding results directory, though this time they are a bit more neatly organized, with only the report and the final output at the top level, and all intermediate files shoved out of the way into a subdirectory.

??? abstract "Directory contents"

    ```console
    results/hello_containers/
    ├── cowpy-COLLECTED-batch-output.txt
    ├── intermediates
    │   ├── Bonjour-output.txt
    │   ├── COLLECTED-batch-output.txt
    │   ├── Hello-output.txt
    │   ├── Holà-output.txt
    │   ├── UPPER-Bonjour-output.txt
    │   ├── UPPER-Hello-output.txt
    │   └── UPPER-Holà-output.txt
    └── batch-report.txt
    ```

The final ASCII art output is in the `results/hello_containers/` directory, under the name `cowpy-COLLECTED-batch-output.txt`.

??? abstract "File contents"

    ```console title="results/hello_containers/cowpy-COLLECTED-batch-output.txt"
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

And there it is, our beautiful turkey saying the greetings as desired.

#### 2.3.4. Inspect how Nextflow launched the containerized task

As a final coda to this section, let's take a look at the work subdirectory for one of the `cowpy` process calls to get a bit more insight on how Nextflow works with containers under the hood.

Check the output from your `nextflow run` command to find the path to the work subdirectory for the `cowpy` process.
Looking at what we got for the run shown above, the console log line for the `cowpy` process starts with `[98/656c6c]`.
That corresponds to the following truncated directory path: `work/98/656c6c`.

In that directory, you will find the `.command.run` file that contains all the commands Nextflow ran on your behalf in the course of executing the pipeline.

??? abstract "File contents"

    ```console title="work/98/656c6c90cce1667c094d880f4b6dcc/.command.run"
    #!/bin/bash
    ### ---
    ### name: 'cowpy'
    ### container: 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
    ### outputs:
    ### - 'cowpy-COLLECTED-batch-output.txt'
    ### ...
    set -e
    set -u
    NXF_DEBUG=${NXF_DEBUG:=0}; [[ $NXF_DEBUG > 1 ]] && set -x
    NXF_ENTRY=${1:-nxf_main}


    nxf_sleep() {
      sleep $1 2>/dev/null || sleep 1;
    }

    nxf_date() {
        local ts=$(date +%s%3N);
        if [[ ${#ts} == 10 ]]; then echo ${ts}000
        elif [[ $ts == *%3N ]]; then echo ${ts/\%3N/000}
        elif [[ $ts == *3N ]]; then echo ${ts/3N/000}
        elif [[ ${#ts} == 13 ]]; then echo $ts
        else echo "Unexpected timestamp value: $ts"; exit 1
        fi
    }

    nxf_env() {
        echo '============= task environment ============='
        env | sort | sed "s/\(.*\)AWS\(.*\)=\(.\{6\}\).*/\1AWS\2=\3xxxxxxxxxxxxx/"
        echo '============= task output =================='
    }

    nxf_kill() {
        declare -a children
        while read P PP;do
            children[$PP]+=" $P"
        done < <(ps -e -o pid= -o ppid=)

        kill_all() {
            [[ $1 != $$ ]] && kill $1 2>/dev/null || true
            for i in ${children[$1]:=}; do kill_all $i; done
        }

        kill_all $1
    }

    nxf_mktemp() {
        local base=${1:-/tmp}
        mkdir -p "$base"
        if [[ $(uname) = Darwin ]]; then mktemp -d $base/nxf.XXXXXXXXXX
        else TMPDIR="$base" mktemp -d -t nxf.XXXXXXXXXX
        fi
    }

    nxf_fs_copy() {
      local source=$1
      local target=$2
      local basedir=$(dirname $1)
      mkdir -p $target/$basedir
      cp -fRL $source $target/$basedir
    }

    nxf_fs_move() {
      local source=$1
      local target=$2
      local basedir=$(dirname $1)
      mkdir -p $target/$basedir
      mv -f $source $target/$basedir
    }

    nxf_fs_rsync() {
      rsync -rRl $1 $2
    }

    nxf_fs_rclone() {
      rclone copyto $1 $2/$1
    }

    nxf_fs_fcp() {
      fcp $1 $2/$1
    }

    on_exit() {
        local last_err=$?
        local exit_status=${nxf_main_ret:=0}
        [[ ${exit_status} -eq 0 && ${nxf_unstage_ret:=0} -ne 0 ]] && exit_status=${nxf_unstage_ret:=0}
        [[ ${exit_status} -eq 0 && ${last_err} -ne 0 ]] && exit_status=${last_err}
        printf -- $exit_status > /workspaces/training/hello-nextflow/work/98/656c6c90cce1667c094d880f4b6dcc/.exitcode
        set +u
        docker rm $NXF_BOXID &>/dev/null || true
        exit $exit_status
    }

    on_term() {
        set +e
        docker stop $NXF_BOXID
    }

    nxf_launch() {
        docker run -i --cpu-shares 1024 -e "NXF_TASK_WORKDIR" -v /workspaces/training/hello-nextflow/work:/workspaces/training/hello-nextflow/work -w "$NXF_TASK_WORKDIR" --name $NXF_BOXID community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273 /bin/bash -ue /workspaces/training/hello-nextflow/work/98/656c6c90cce1667c094d880f4b6dcc/.command.sh
    }

    nxf_stage() {
        true
        # stage input files
        rm -f COLLECTED-batch-output.txt
        ln -s /workspaces/training/hello-nextflow/work/7f/f435e3f2cf95979b5f3d7647ae6696/COLLECTED-batch-output.txt COLLECTED-batch-output.txt
    }

    nxf_unstage_outputs() {
        true
    }

    nxf_unstage_controls() {
        true
    }

    nxf_unstage() {
        if [[ ${nxf_main_ret:=0} == 0 ]]; then
            (set -e -o pipefail; (nxf_unstage_outputs | tee -a .command.out) 3>&1 1>&2 2>&3 | tee -a .command.err)
            nxf_unstage_ret=$?
        fi
        nxf_unstage_controls
    }

    nxf_main() {
        trap on_exit EXIT
        trap on_term TERM INT USR2
        trap '' USR1

        [[ "${NXF_CHDIR:-}" ]] && cd "$NXF_CHDIR"
        export NXF_BOXID="nxf-$(dd bs=18 count=1 if=/dev/urandom 2>/dev/null | base64 | tr +/ 0A | tr -d '\r\n')"
        NXF_SCRATCH=''
        [[ $NXF_DEBUG > 0 ]] && nxf_env
        touch /workspaces/training/hello-nextflow/work/98/656c6c90cce1667c094d880f4b6dcc/.command.begin
        set +u
        set -u
        [[ $NXF_SCRATCH ]] && cd $NXF_SCRATCH
        export NXF_TASK_WORKDIR="$PWD"
        nxf_stage

        set +e
        (set -o pipefail; (nxf_launch | tee .command.out) 3>&1 1>&2 2>&3 | tee .command.err) &
        pid=$!
        wait $pid || nxf_main_ret=$?
        nxf_unstage
    }

    $NXF_ENTRY

    ```

If you search for `nxf_launch` in this file, you should see something like this:

```console
nxf_launch() {
    docker run -i --cpu-shares 1024 -e "NXF_TASK_WORKDIR" -v /workspaces/training/hello-nextflow/work:/workspaces/training/hello-nextflow/work -w "$NXF_TASK_WORKDIR" --name $NXF_BOXID community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273 /bin/bash -ue /workspaces/training/hello-nextflow/work/98/656c6c90cce1667c094d880f4b6dcc/.command.sh
}
```

As you can see, Nextflow is using the `docker run` command to launch the process call.
It also mounts the corresponding work subdirectory into the container, sets the working directory inside the container accordingly, and runs our templated bash script in the `.command.sh` file.

All the hard work we had to do manually in the first section? Nextflow does it for us behind the scenes!

```txt
 _______________________
< Hurray for robots...! >
 -----------------------
                                   ,-----.
                                   |     |
                                ,--|     |-.
                         __,----|  |     | |
                       ,;::     |  `_____' |
                       `._______|    i^i   |
                                `----| |---'| .
                           ,-------._| |== ||//
                           |       |_|P`.  /'/
                           `-------' 'Y Y/'/'
                                     .==\ /_\
   ^__^                             /   /'|  `i
   (oo)\_______                   /'   /  |   |
   (__)\       )\/\             /'    /   |   `i
       ||----w |           ___,;`----'.___L_,-'`\__
       ||     ||          i_____;----\.____i""\____\
```

### Takeaway

You know how to use containers in Nextflow to run processes.

### What's next?

Take a break!

When you're ready, move on to [**Part 6: Hello Config**](./06_hello_config.md) to learn how to configure the execution of your pipeline to fit your infrastructure as well as manage configuration of inputs and parameters.

It's the very last part, and then you'll be done with this course!

---

## Quiz

<quiz>
What is a container?
- [ ] A type of virtual machine
- [ ] A file compression format
- [x] A lightweight, standalone executable unit that includes everything needed to run an application
- [ ] A network protocol
</quiz>

<quiz>
What is the difference between a container image and a container instance?
- [ ] They are the same thing
- [x] An image is a template; an instance is a running container created from that image
- [ ] An instance is a template; an image is a running container
- [ ] Images are for Docker; instances are for Singularity
</quiz>

<quiz>
What does the `-v` flag do in a `docker run` command?
- [ ] Enables verbose output
- [ ] Validates the container
- [x] Mounts a volume from the host system into the container
- [ ] Specifies the version of the container

Learn more: [1.3.4. Mount data into the container](#134-mount-data-into-the-container)
</quiz>

<quiz>
Why do you need to mount volumes when using containers?
- [ ] To improve container performance
- [ ] To save disk space
- [x] Because containers are isolated from the host filesystem by default
- [ ] To enable networking

Learn more: [1.3.4. Mount data into the container](#134-mount-data-into-the-container)
</quiz>

<quiz>
How do you specify a container for a Nextflow process?
- [ ] docker 'container-uri'
- [ ] image 'container-uri'
- [x] container 'container-uri'
- [ ] use 'container-uri'

Learn more: [2.3.1. Specify a container for cowpy](#231-specify-a-container-for-cowpy)
</quiz>

<quiz>
How do you enable Docker in Nextflow?
- [ ] Add `--docker` to the command line
- [ ] Set `process.docker = true` in the config
- [x] Set `docker.enabled = true` in the config
- [ ] Docker is enabled by default

Learn more: [2.3.2. Enable use of Docker via the `nextflow.config` file](#232-enable-use-of-docker-via-the-nextflowconfig-file)
</quiz>

<quiz>
What does Nextflow automatically handle when running a process in a container? (Select all that apply)
- [x] Pulling the container image if needed
- [x] Mounting the work directory
- [x] Running the process script inside the container
- [x] Cleaning up the container instance after execution

Learn more: [2.3.4. Inspect how Nextflow launched the containerized task](#234-inspect-how-nextflow-launched-the-containerized-task)
</quiz>
