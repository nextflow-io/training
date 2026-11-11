# Part 4: Hello Containers

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/Xqr--bKEN9U?si=y8lAedhEHWaTV4zd&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } See [the whole playlist](https://youtube.com/playlist?list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&si=eF7cLR62goy-lc6n) on the Nextflow YouTube channel.

:green_book: The video transcript is available [here](./transcripts/05_hello_containers.md).
///

In Parts 1-3 of this training course, you learned how to use the basic building blocks of Nextflow to assemble a simple workflow capable of processing some text, parallelizing execution if there were multiple inputs, and collecting the results for further processing.

However, you were limited to basic UNIX tools available in your environment.
Real-world tasks often require various tools and packages not included by default.
Typically, you'd need to install these tools, manage their dependencies, and resolve any conflicts.

That is all very tedious and annoying, so we're going to show you how to use **containers** to solve this problem much more conveniently.

A **container** is a lightweight, standalone, executable unit of software created from a container **image** that includes everything needed to run an application including code, system libraries and settings.
As you might imagine, that is going to be very helpful for making your pipelines more reproducible.

Note that we'll be teaching this using [Docker](https://www.docker.com/get-started/), but keep in mind Nextflow supports [several other container technologies](https://nextflow.io/docs/latest/container.html) as well.

Here is where this part fits in the course:

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-course-map-part4.svg"
</figure>

??? info "How to begin from this section"

    This section of the course assumes you have completed Parts 1-3 of the [Hello Nextflow](./index.md) course, but if you are comfortable with the basics covered in those sections, you can start from here without doing anything special.

---

## 0. Warmup: Run `hello-containers.nf`

We're going to use the workflow script `hello-containers.nf` as a starting point.
It is equivalent to the script produced by working through Part 3 of this training course.

Just to make sure everything is working, run the script once before making any changes:

```bash
nextflow run hello-containers.nf
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `hello-containers.nf` [evil_poitras] revision: bb1df0f734

    executor >  local (7)
    [2d/7e3ffe] sayHello (2)          | 3 of 3 ✔
    [f4/e5528d] translateGreeting (1) | 3 of 3 ✔
    [93/a1e241] collectGreetings      | 1 of 1 ✔

    Outputs:

      /workspaces/training/hello-nextflow/results

      greetings:
        - batch/intermediates/Mateo.txt
        - batch/intermediates/Celine.txt
        - batch/intermediates/Raphael.txt

      translated:
        - batch/intermediates/es-Raphael.txt
        - batch/intermediates/es-Celine.txt
        - batch/intermediates/es-Mateo.txt

      collected: batch/intermediates/batch-es-greetings.txt

      batch_report: batch/batch-es-report.txt
    ```

As previously, you will find the output files in the locations specified in the `output` block, under `results/batch/`, since we didn't set a value for `--batch`.

??? abstract "Directory contents"

    ```console
    results/batch
    ├── batch-es-report.txt
    └── intermediates
        ├── batch-es-greetings.txt
        ├── Celine.txt
        ├── es-Celine.txt
        ├── es-Mateo.txt
        ├── es-Raphael.txt
        ├── Mateo.txt
        └── Raphael.txt
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
echo "Hello Containers" | cowpy -c tux
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

This is neat. What would be even neater is if we could feed our `names.csv` file as input into this.
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

To mount a volume, we replace the paths and add the volume mounting argument to the `docker run` command as follows:

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
    data               hello-containers.nf  hello-workflow.nf  nextflow.config  solutions
    hello-channels.nf  hello-modules.nf     hello-world.nf     results          work
    ```

You can now see the contents of the working directory from inside the container, including the `names.csv` file under `data/`.

This effectively established a tunnel through the container wall that you can use to access that part of your filesystem.

#### 1.3.5. Use the mounted data

Now that we have mounted the working directory into the container, we can use the `cowpy` command to display the contents of the `names.csv` file.

To do this, we'll use `cat /my_project/data/names.csv | ` to pipe the contents of the CSV file into the `cowpy` command.

```bash
cat /my_project/data/names.csv | cowpy -c turkey
```

??? success "Command output"

    ```console
     ___________________
    / Mateo,es,squirrel \
    | Raphael,it,turtle |
    \ Celine,fr,moose   /
     -------------------
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

This produces ASCII art of a turkey reciting the contents of our names file.
It repeats the full rows rather than greetings, but our Nextflow workflow will take care of that.

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
--8<-- "docs/en/docs/hello_nextflow/img/hello-pipeline-cowpy.svg"
</figure>

### 2.1. Add a `cowpy` process to the workflow script

We can model our `cowpy` process on the other processes we've written previously.

Add the following process definition to the workflow script, right under the `collectGreetings` process:

```groovy title="hello-containers.nf" linenums="80"
/*
 * Generate ASCII art with cowpy (https://github.com/jeffbuttars/cowpy)
 */
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

Now we need to call the process in the workflow, give it a character to use, and publish its output.

#### 2.2.1. Add a call to the `cowpy` process in the workflow

Let's connect the `cowpy()` process to the output of the `collectGreetings()` process, which as you may recall produces two outputs:

- `collectGreetings.out.outfile` contains the output file <--_what we want_
- `collectGreetings.out.report` contains the report file with the count of greetings per batch

In the workflow block, make the following code change:

=== "After"

    ```groovy title="hello-containers.nf" linenums="109" hl_lines="12-13"
        main:
        // create a channel for inputs from a CSV file
        name_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { row -> row[0] }
        // emit a greeting
        sayHello(name_ch)
        // translate the greeting
        translateGreeting(sayHello.out, params.language)
        // collect all the greetings into one file
        collectGreetings(translateGreeting.out.collect(), params.batch, params.language)
        // generate ASCII art of the greetings with cowpy
        cowpy(collectGreetings.out.outfile, params.character)
    ```

=== "Before"

    ```groovy title="hello-containers.nf" linenums="109"
        main:
        // create a channel for inputs from a CSV file
        name_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { row -> row[0] }
        // emit a greeting
        sayHello(name_ch)
        // translate the greeting
        translateGreeting(sayHello.out, params.language)
        // collect all the greetings into one file
        collectGreetings(translateGreeting.out.collect(), params.batch, params.language)
    ```

Notice that we used a new CLI parameter, `params.character`, in order to specify which character we want to have say the greetings.

#### 2.2.2. Declare the `character` parameter with a default value

Declare the new parameter in the `params` block of the workflow script, with a default value like the other parameters:

=== "After"

    ```groovy title="hello-containers.nf" linenums="98" hl_lines="8"
    /*
     * Pipeline parameters
     */
    params {
        input: Path = 'data/names.csv'
        language: String = 'es'
        batch: String = 'batch'
        character: String = 'turkey'
    }
    ```

=== "Before"

    ```groovy title="hello-containers.nf" linenums="98"
    /*
     * Pipeline parameters
     */
    params {
        input: Path = 'data/names.csv'
        language: String = 'es'
        batch: String = 'batch'
    }
    ```

Now we can skip typing the character parameter in our command lines.

#### 2.2.3. Update the workflow outputs

We need to update the workflow outputs to publish the output of the `cowpy` process.

##### 2.2.3.1. Update the `publish:` section

In the `workflow` block, make the following code change:

=== "After"

    ```groovy title="hello-containers.nf" linenums="124" hl_lines="6"
        publish:
        greetings = sayHello.out
        translated = translateGreeting.out
        collected = collectGreetings.out.outfile
        batch_report = collectGreetings.out.report
        cowpy_art = cowpy.out
    ```

=== "Before"

    ```groovy title="hello-containers.nf" linenums="124"
        publish:
        greetings = sayHello.out
        translated = translateGreeting.out
        collected = collectGreetings.out.outfile
        batch_report = collectGreetings.out.report
    ```

The `cowpy` process only produces one output so we can refer to it the usual way by appending `.out`.

##### 2.2.3.2. Update the `output` block

We need to add the final `cowpy_art` output to the `output` block.
It is a final result, so we publish it at the top level of the batch directory, next to the report, and as a copy like the other outputs.

In the `output` block, make the following code change:

=== "After"

    ```groovy title="hello-containers.nf" linenums="132" hl_lines="18-21"
    output {
        greetings {
            path "${params.batch}/intermediates"
            mode 'copy'
        }
        translated {
            path "${params.batch}/intermediates"
            mode 'copy'
        }
        collected {
            path "${params.batch}/intermediates"
            mode 'copy'
        }
        batch_report {
            path "${params.batch}"
            mode 'copy'
        }
        cowpy_art {
            path "${params.batch}"
            mode 'copy'
        }
    }
    ```

=== "Before"

    ```groovy title="hello-containers.nf" linenums="132"
    output {
        greetings {
            path "${params.batch}/intermediates"
            mode 'copy'
        }
        translated {
            path "${params.batch}/intermediates"
            mode 'copy'
        }
        collected {
            path "${params.batch}/intermediates"
            mode 'copy'
        }
        batch_report {
            path "${params.batch}"
            mode 'copy'
        }
    }
    ```

#### 2.2.4. Run the workflow

Do you think it's going to work?

Let's delete the previous published outputs to have a clean slate, and run the workflow with the `-resume` flag.

```bash
rm -r results/batch/
nextflow run hello-containers.nf -resume
```

??? failure "Command output"

    ```console hl_lines="10 13 20-21 26-27"
    N E X T F L O W   ~  version 26.04.4

    Launching `hello-containers.nf` [distraught_salas] revision: ef370ae430

    executor >  local (1)
    [ce/83663a] sayHello (1)          | 3 of 3, cached: 3 ✔
    [c2/d82633] translateGreeting (1) | 3 of 3, cached: 3 ✔
    [93/a1e241] collectGreetings      | 1 of 1, cached: 1 ✔
    [72/4be2f7] cowpy                 | 0 of 1 ✘
    ERROR ~ Error executing process > 'cowpy'

    Caused by:
      Process `cowpy` terminated with an error exit status (127)


    Command executed:

      cat batch-es-greetings.txt | cowpy -c "turkey" > cowpy-batch-es-greetings.txt

    Command exit status:
      127

    Command output:
      (empty)

    Command error:
      .command.sh: line 2: cowpy: command not found

    Work dir:
      /workspaces/training/hello-nextflow/work/72/4be2f70a8676312e8dde1b1f9522b4

    Tip: when you have fixed the problem you can continue the execution adding the option `-resume` to the run command line

     -- Check '.nextflow.log' file for details
    ```

Oh no, there's an error!
The error code given by `error exit status (127)` means the executable we asked for was not found.

That makes sense, since we're calling the `cowpy` tool but we haven't actually specified a container yet (oops).

### 2.3. Use a container to run the `cowpy` process

We need to specify a container and tell Nextflow to use it for the `cowpy()` process.

#### 2.3.1. Specify a container for `cowpy`

We can use the same image we were using directly in the first section of this tutorial.

Edit the `cowpy` process definition in `hello-containers.nf` to add the `container` directive as follows:

=== "After"

    ```groovy title="hello-containers.nf" linenums="83" hl_lines="3"
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

    ```groovy title="hello-containers.nf" linenums="83"
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
We do that in the `nextflow.config` file, which Nextflow loads automatically when it is present in the current directory; [Part 2 of Nextflow Run](../nextflow_run/02_configure_pipeline.md) covers it in detail.

The `nextflow.config` file we provided contains a single line, which explicitly disables Docker: `docker.enabled = false`.

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
    N E X T F L O W   ~  version 26.04.4

    Launching `hello-containers.nf` [drunk_spence] revision: bfe2e017c2

    executor >  local (1)
    [ce/83663a] sayHello (1)          | 3 of 3, cached: 3 ✔
    [f4/e5528d] translateGreeting (2) | 3 of 3, cached: 3 ✔
    [93/a1e241] collectGreetings      | 1 of 1, cached: 1 ✔
    [fd/d71f5c] cowpy                 | 1 of 1 ✔

    Outputs:

      /workspaces/training/hello-nextflow/results

      greetings:
        - batch/intermediates/Celine.txt
        - batch/intermediates/Mateo.txt
        - batch/intermediates/Raphael.txt

      translated:
        - batch/intermediates/es-Celine.txt
        - batch/intermediates/es-Raphael.txt
        - batch/intermediates/es-Mateo.txt

      collected: batch/intermediates/batch-es-greetings.txt

      batch_report: batch/batch-es-report.txt

      cowpy_art: batch/cowpy-batch-es-greetings.txt
    ```

This time it does indeed work!
As usual you can find the workflow outputs in the results directory, with the report and the final output at the top level of the batch directory, and the intermediate files in a subdirectory.

??? abstract "Directory contents"

    ```console
    results/batch
    ├── batch-es-report.txt
    ├── cowpy-batch-es-greetings.txt
    └── intermediates
        ├── batch-es-greetings.txt
        ├── Celine.txt
        ├── es-Celine.txt
        ├── es-Mateo.txt
        ├── es-Raphael.txt
        ├── Mateo.txt
        └── Raphael.txt
    ```

The final ASCII art output is in the `results/batch/` directory, under the name `cowpy-batch-es-greetings.txt`.

??? abstract "File contents"

    ```console title="results/batch/cowpy-batch-es-greetings.txt"
     ________________
    / Hola, Raphael! \
    | Hola, Celine!  |
    \ Hola, Mateo!   /
     ----------------
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

And there it is, our beautiful turkey saying the translated greetings as desired.

#### 2.3.4. Inspect how Nextflow launched the containerized task

As a final coda to this section, let's take a look at the work subdirectory for one of the `cowpy` process calls to get a bit more insight on how Nextflow works with containers under the hood.

Check the output from your `nextflow run` command to find the path to the work subdirectory for the `cowpy` process.
Looking at what we got for the run shown above, the console log line for the `cowpy` process starts with `[fd/d71f5c]`.
That corresponds to the following truncated directory path: `work/fd/d71f5c`.

In that directory, you will find the `.command.run` file that contains all the commands Nextflow ran on your behalf in the course of executing the pipeline.

??? abstract "File contents"

    ```console title="work/fd/d71f5c90fe222eedd0bb39f8ea8f15/.command.run"
    #!/bin/bash
    ### ---
    ### name: 'cowpy'
    ### container: 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
    ### outputs:
    ### - 'cowpy-batch-es-greetings.txt'
    ### ...
    set -e
    set -u
    NXF_DEBUG=${NXF_DEBUG:=0}; [[ $NXF_DEBUG > 1 ]] && set -x
    NXF_ENTRY=${1:-nxf_main}


    nxf_sleep() {
      sleep $1 2>/dev/null || sleep 1;
    }

    nxf_date() {
        local s n
        read s n < <(date '+%s %N')
        if [[ $n =~ ^[0-9]+$ ]]; then
            n=$(printf '%09d' "$((10#$n))")
            echo "$((s * 1000 + 10#${n:0:3}))"
        else
            echo "$((s * 1000))"
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
      local source=$1
      local target=$2
      mkdir -p $target
      rsync -rRl $source $target
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
        printf -- $exit_status > /workspaces/training/hello-nextflow/work/fd/d71f5c90fe222eedd0bb39f8ea8f15/.exitcode
        set +u
        docker rm $NXF_BOXID &>/dev/null || true
        exit $exit_status
    }

    on_term() {
        set +e
        docker stop $NXF_BOXID
    }

    nxf_launch() {
        docker run -i --cpu-shares 1024 -e "NXF_TASK_WORKDIR" -v /workspaces/training/hello-nextflow/work:/workspaces/training/hello-nextflow/work -w "$NXF_TASK_WORKDIR" --name $NXF_BOXID community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273 /bin/bash -ue /workspaces/training/hello-nextflow/work/fd/d71f5c90fe222eedd0bb39f8ea8f15/.command.sh
    }

    nxf_stage() {
        true
        # stage input files
        rm -f batch-es-greetings.txt
        ln -s /workspaces/training/hello-nextflow/work/93/a1e241bce7ab57b6d7fe663c99ef8b/batch-es-greetings.txt batch-es-greetings.txt
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
        touch /workspaces/training/hello-nextflow/work/fd/d71f5c90fe222eedd0bb39f8ea8f15/.command.begin
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
    docker run -i --cpu-shares 1024 -e "NXF_TASK_WORKDIR" -v /workspaces/training/hello-nextflow/work:/workspaces/training/hello-nextflow/work -w "$NXF_TASK_WORKDIR" --name $NXF_BOXID community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273 /bin/bash -ue /workspaces/training/hello-nextflow/work/fd/d71f5c90fe222eedd0bb39f8ea8f15/.command.sh
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

### 2.4. Add a Conda package as an alternative

Not everyone can use containers.
A pipeline can also supply its software through Conda instead, when the person running it enables Conda, as [Part 2 of Nextflow Run](../nextflow_run/02_configure_pipeline.md) shows.
For that to work, the process needs to declare a Conda package as well as a container.

The `cowpy` tool is available as the Conda package `conda-forge::cowpy==1.1.5`.
Add it to the `cowpy` process with the `conda` directive:

=== "After"

    ```groovy title="hello-containers.nf" linenums="83" hl_lines="4"
    process cowpy {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
        conda 'conda-forge::cowpy==1.1.5'

        input:
    ```

=== "Before"

    ```groovy title="hello-containers.nf" linenums="83"
    process cowpy {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

        input:
    ```

This doesn't replace the container: Nextflow uses whichever one matches the software packaging technology enabled in the configuration.
With `docker.enabled = true`, as in our `nextflow.config`, it keeps using the container.
In Part 5, you will add a profile that switches the pipeline to Conda.

!!! tip

    The [Seqera Containers](https://seqera.io/containers/) search shows the Conda package URI for a tool, even if you're not planning to build a container from it.

### Takeaway

You know how to use containers in Nextflow to run processes, and how to declare a Conda package as an alternative.

The pipeline is now complete, with all four steps:

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_pipeline_complete.svg"
</figure>

### What's next?

Take a short break if you feel like it.

When you're ready, move on to [**Part 5: Hello Modules**](./05_hello_modules.md) to organize your code into modules and your parameter values into the configuration, giving your pipeline the same layout as the finished Hello pipeline in Nextflow Run.

---

## Bonus exercises

These optional exercises let you practice what you learned in this part.
Each one can be solved with what this part covered; open the hint if you get stuck.
The solutions are also in `solutions/4-hello-containers/bonus/`.

### Exercise 1: Run the translation step in a container

The `translateGreeting` process still runs with the Python installed in your environment.
Starting from `hello-containers.nf` as you left it at the end of this part, make `translateGreeting` run in a container that provides Python.
Run the workflow, then check the `.command.run` file of one of the `translateGreeting` tasks to confirm which container it used.

??? tip "Hint"

    Add a `container` directive to the process, the same way you did for `cowpy`.
    Any image that provides Python 3 works, for example `python:3.12-slim`.

??? solution

    Add the `container` directive to the `translateGreeting` process:

    ```groovy title="hello-containers.nf" linenums="23" hl_lines="3"
    process translateGreeting {

        container 'python:3.12-slim'

        input:
    ```

    Then run the workflow:

    ```bash
    nextflow run hello-containers.nf
    ```

    The header of the `.command.run` file of a `translateGreeting` task names the container:

    ```console title="work/e3/8fd4702290fce24d0540d99a944c79/.command.run"
    #!/bin/bash
    ### ---
    ### name: 'translateGreeting (3)'
    ### container: 'python:3.12-slim'
    ### outputs:
    ### - 'es-Mateo.txt'
    ```

    Now both `translateGreeting` and `cowpy` run in containers, each with its own image.

### Exercise 2: Use a different character

The ASCII art uses the turkey, the default value of the `character` parameter.
Pick another character that `cowpy` offers, and run the workflow so that it draws that character instead.

??? tip "Hint"

    Section 1.3.2 lists characters you can pick with the `-c` flag.
    Set the `character` parameter on the command line, or change its default value in the `params` block of `hello-containers.nf`.

??? solution

    The list of characters in section 1.3.2 includes `stegosaurus`.
    Set it on the command line:

    ```bash
    nextflow run hello-containers.nf --character stegosaurus
    ```

    Or change the default value in the `params` block of `hello-containers.nf`:

    ```groovy title="hello-containers.nf" linenums="101" hl_lines="8"
    /*
     * Pipeline parameters
     */
    params {
        input: Path = 'data/names.csv'
        language: String = 'es'
        batch: String = 'batch'
        character: String = 'stegosaurus'
    }
    ```

    The ASCII art now shows a stegosaurus:

    ```console title="results/batch/cowpy-batch-es-greetings.txt"
     ________________
    / Hola, Celine!  \
    | Hola, Mateo!   |
    \ Hola, Raphael! /
     ----------------
    \                             .       .
     \                           / `.   .' "
      \                  .---.  <    > <    >  .---.
       \                 |    \  \ - ~ ~ - /  /    |
             _____          ..-~             ~-..-~
            |     |   \~~~\.'                    `./~~~/
           ---------   \__/                        \__/
          .'  O    \     /               /       \  "
         (_____,    `._.'               |         }  \/~~~/
          `----.          /       }     |        /    \__/
                `-.      |       /      |       /      `. ,~~|
                    ~-.__|      /_ - ~ ^|      /- _      `..-'
                         |     /        |     /     ~-.     `-. _  _  _
                         |_____|        |_____|         ~ - . _ _ _ _ _>
    ```

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
- [ ] `docker 'container-uri'`
- [ ] `image 'container-uri'`
- [x] `container 'container-uri'`
- [ ] `use 'container-uri'`

Learn more: [2.3.1. Specify a container for cowpy](#231-specify-a-container-for-cowpy)
</quiz>

<quiz>
Which `nextflow.config` setting enables Docker for your workflow?
- [ ] `#!groovy process.docker = true`
- [x] `#!groovy docker.enabled = true`
- [ ] `#!groovy container.engine = 'docker'`
- [ ] `#!groovy docker.activate = true`

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

<quiz>
What happens when a process declares both a `container` and a `conda` directive?
- [ ] Nextflow runs the process twice, once with each
- [ ] Nextflow fails with an error
- [x] Nextflow uses whichever one matches the software packaging technology enabled in the configuration
- [ ] The `conda` directive always takes precedence

Learn more: [2.4. Add a Conda package as an alternative](#24-add-a-conda-package-as-an-alternative)
</quiz>
