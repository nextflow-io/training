# Part 5: Hello Containers

In Parts 1-4 of this training course, you learned how to use the basic building blocks of Nextflow to assemble a simple pipeline capable of processing some text, parallelizing execution if there were multiple inputs and collecting the results for further processing.

However, you were limited to basic UNIX tools available in your environment.
Real-world tasks often require various tools and packages not included by default.
Typically, you'd need to install these tools, manage their dependencies, and resolve any conflicts.

That is all very tedious and annoying, so we're going to show you how to use **containers** to solve this problem much more conveniently.

A **container** is a lightweight, standalone, executable unit of software created from a container **image** that includes everything needed to run an application including code, system libraries and settings.

!!! Note

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
 N E X T F L O W   ~  version 24.10.0

Launching `hello-containers.nf` [tender_becquerel] DSL2 - revision: f7cat8e223

executor >  local (7)
[bd/4bb541] sayHello (1)         [100%] 3 of 3 ✔
[85/b627e8] convertToUpper (3)   [100%] 3 of 3 ✔
[7d/f7961c] collectGreetings     [100%] 1 of 1 ✔
```

Our goal will be to add a step to this workflow that will use a container for execution.
However, we are first going to go over some basic concepts and operations to solidify your understanding of what containers are before we start using them in Nextflow.

---

## 1. Use a container 'manually'

### 1.1. Pull the container image

To use a container, you usually download or "pull" a container image from a container registry, and then run the container image to create a container instance.

The general syntax is as follows:

```bash
docker pull '<container>'
```

The `docker pull` part is the instruction to the container system to pull a container image from a repository.

The `'<container>'` part is the URI address of the container image.

As an example, let's pull a container image that contains the [`cowsay` tool](https://pypi.org/project/cowsay/), which generates ASCII art to display arbitrary text inputs in a fun way.

There are various repositories where you can find published containers.
We looked in the [Seqera Containers](https://seqera.io/containers/) repository and found this `cowsay` container: `'community.wave.seqera.io/library/pip_cowsay:131d6a1b707a8e65'`.

The pull command becomes:

```bash
docker pull 'community.wave.seqera.io/library/pip_cowsay:131d6a1b707a8e65'
```

Running this gives you the following console output as the system downloads the image:

TODO OUTPUT

Once the download is complete, you have a local copy of the container image.

### 1.2. Use the container to run `cowsay` as a one-off command

One very common way that people use containers is to run them directly, _i.e._ non-interactively.
This is great for running one-off commands.

The general syntax is as follows:

```bash
docker run --rm '<container>' [tool command]
```

The `docker run --rm '<container>'` part is the instruction to the container system to spin up a container instance from a container image and execute a command in it.
The `--rm` flag tells the system to shut down the container instance after the command has completed.

The `[tool command]` syntax depends on the tool you are using and how the container is set up.
Here we will use `cowsay -t "Hello World"`.

Fully assembled, the container execution command looks like this:

```bash
docker run --rm 'community.wave.seqera.io/library/pip_cowsay:131d6a1b707a8e65' cowsay -t "Hello World"
```

Run it to produce the following output:

```console title="Output"
 _____________
< Hello World >
 -------------
        \   ^__^
         \  (oo)\_______
            (__)\       )\/\
                ||----w |
                ||     ||
```

The system spun up the container, ran the `cowsay` command with the parameters we specified, sent the output to the console and finally, shut down the container instance.

### 1.3. Use the container to run `cowsay` interactively

You can also run a container interactively, which gives you a shell prompt inside the container.

#### 1.3.1. Spin up the container

To run interactively, we just add `-it` to the `docker pull` command.
Optionally, we can specify the shell we want to use inside the container by appending _e.g._ `/bin/bash` to the command.

```bash
docker run --rm -it 'community.wave.seqera.io/library/pip_cowsay:131d6a1b707a8e65' /bin/bash
```

Notice that the prompt has changed to `(base) root@b645838b3314:/tmp#`, which indicates that you are now inside the container.

You can verify this by running `ls` to list directory contents:

```bash
ls
```

You can see that the filesystem inside the container is different from the filesystem on your host system:

```console title="Output"
(base) root@b645838b3314:/tmp# ls /
bin    dev    etc    home   lib    media  mnt    opt    proc   root   run    sbin   srv    sys    tmp    usr    var
```

#### 1.3.2. Run the desired tool command(s)

Now that you are inside the container, you can run the `cowsay` command directly.

```bash
cowsay -t "Hello World" -c tux
```

Now the output shows the Linux penguin, Tux, instead of the default cow, because we specified the `-c` parameter.

```console title="Output"
  ___________
| Hello World |
  ===========
                \
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

Because you're inside the container, you can run the cowsay command as many times as you like, varying the input parameters, without having to bother with docker commands.

!!! Tip

    Use the '-c' flag to pick a different character from this list:
    `beavis`, `cheese`, `cow`, `daemon`, `dragon`, `fox`, `ghostbusters`, `kitty`, `meow`, `miki`, `milk`, `octopus`, `pig`, `stegosaurus`, `stimpy`, `trex`, `turkey`, `turtle`, `tux`

#### 1.3.3. Exit the container

To exit the container, you can type `exit` at the prompt or use the ++ctrl+d++ keyboard shortcut.

```bash
exit
```

Your prompt should now be back to what it was before you started the container.

#### 1.3.4. Mounting data into containers

When you run a container, it is isolated from the host system by default.
This means that the container can't access any files on the host system unless you explicitly tell it to.
One way to do this is to **mount** a **volume** from the host system into the container.

To mount a volume, we add `-v <outside_path>:<inside_path>` to the `docker run command as follows:

```bash
docker run --rm -it -v data:/data 'community.wave.seqera.io/library/pip_cowsay:131d6a1b707a8e65' /bin/bash
```

This mounts the `data` directory (in the current working directory) as a volume that will be found under `/data` inside the container.

You can check that it works by listing the contents of `/data`:

```bash
ls /data
```

You should now be able to see the contents of the `data` directory from inside the container:

```console title="Output"
greetings.csv
```

This effectively established a tunnel through the container wall that you can use to access that part of your filesystem.

#### 1.3.5. Use the mounted data

Now that we have mounted the `data` directory into the container, we can use the `cowsay` command to display the contents of the `greetings.csv` file.

To do this, we'll use the syntax `-t "$(cat data/greetings.csv)"` to output the contents of the file into the `cowsay` command.

```bash
cowsay -t "$(cat /data/greetings.csv)" -c pig
```

This produces the desired ASCII art of the pig rattling off our example greetings:

```console title="Output"
  __________________
| Hello,Bonjour,Holà |
  ==================
                  \
                   \
                    \
                     \
                               ,.
                              (_|,.
                              ,' /, )_______   _
                          __j o``-'        `.'-)'
                          (")                 \'
                          `-j                |
                            `-._(           /
                               |_\  |--^.  /
                              /_]'|_| /_)_/
                                  /_]'  /_]'
```

Feel free to play around with this command.
When you're done, exit the container as previously:

```bash
exit
```

You will find yourself back in your normal shell.

### Takeaway

You know how to pull a container and run it either as a one-off or interactively. You also know how to make your data accessible from within your container, which lets you try any tool you're interested in without having to install any software on your system.

### What's next?

Learn how to use containers for the execution of Nextflow processes.

---

## 2. Use containers in Nextflow

Nextflow has built-in support for running processes inside containers to let you run tools you don't have installed in your compute environment.
This means that you can use any container image you like to run your processes, and Nextflow will take care of pulling the image, mounting the data, and running the process inside it.

To demonstrate this, we are going to add a `cowsay` step to the pipeline we've been developing, after the `collectGreetings` step.

### 2.1. Write a `cowSay` module

#### 2.1.1. Create a file stub for the new module

Create an empty file for the module called `cowSay.nf`.

```bash
touch modules/cowsay.nf
```

This gives us a place to put the process code.

#### 2.1.2. Copy the `cowSay` process code in the module file

We can model our `cowSay` process off of the processes we've written previously.

```groovy title="modules/cowSay.nf" linenums="1"
#!/usr/bin/env nextflow

// Generate ASCII art with cowsay
process cowSay {

    publishDir 'results', mode: 'copy'

    input:
        path input_file
        val character

    output:
        path "cowsay-${input_file}"

    script:
    """
    cowsay -c "$character" -t "\$(cat $input_file)" > cowsay-${input_file}
    """

}
```

The output will be a new text file containing the ASCII art generated by the `cowsay` tool.

### 2.2. Import the `cowSay` process into `hello-containers.nf`

Insert the import declaration above the workflow block and fill it out appropriately.

_Before:_

```groovy title="hello-containers.nf" linenums="73"
// Include modules
include { sayHello } from './modules/sayHello.nf'
include { convertToUpper } from './modules/convertToUpper.nf'
include { collectGreetings } from './modules/collectGreetings.nf'

workflow {
```

_After:_

```groovy title="hello-containers.nf" linenums="73"
// Include modules
include { sayHello } from './modules/sayHello.nf'
include { convertToUpper } from './modules/convertToUpper.nf'
include { collectGreetings } from './modules/collectGreetings.nf'
include { cowSay } from './modules/cowSay.nf'

workflow {
```

### 2.3 Add a call to the `cowSay` process in the workflow

Let's connect the `cowSay()` process to the output of the `collectGreetings()` process, which as you may recall produces two outputs:

-   `collectGreetings.out.outfile` contains the output file
-   `collectGreetings.out.count` contains the count of greetings per batch

In the workflow block, make the following code change:

_Before:_

```groovy title="hello-containers.nf" linenums="82"
    // collect all the greetings into one file
    collectGreetings(convertToUpper.out.collect(), params.batch)

    // emit a message about the size of the batch
    collectGreetings.out.count.view{ "There were $it greetings in this batch" }
```

_After:_

```groovy title="hello-containers.nf" linenums="82"
    // collect all the greetings into one file
    collectGreetings(convertToUpper.out.collect(), params.batch)

    // emit a message about the size of the batch
    collectGreetings.out.count.view{ "There were $it greetings in this batch" }

    // generate ASCII art of the greetings with cowSay
    cowSay(collectGreetings.out.outfile, params.character)
```

Notice that we include a new CLI parameter, `params.character`, in order to specify which character we want to have say the greetings.

### 2.4. Run the workflow to verify that it works

Run this with the `-resume` flag.

```bash
nextflow run hello-containers.nf -resume
```

Oh no, there's an error!

```console title="Output"
TODO
```

Of course, we're calling the `cowsay` tool but we haven't actually specified a container.

### 2.5. Specify a container for the process to use

Edit the `cowSay.nf` module to add the `container` directive as follows:

_Before:_

```groovy title="modules/cowSay.nf"
process cowSay {

    publishDir 'containers/results', mode: 'copy'
```

_After:_

```groovy title="modules/cowSay.nf"
process cowSay {

    publishDir 'containers/results', mode: 'copy'
    container 'community.wave.seqera.io/library/pip_cowsay:131d6a1b707a8e65'
```

This tells Nextflow that if Docker is available, it should use the container image specified here to execute the process.

### 2.6. Run the workflow again to verify that it works this time

Run this with the `-resume` flag and with `-with-docker`, which enables Docker execution on a per-command basis.
We will cover persistent configuration in the next section of this course.

```bash
nextflow run hello-containers.nf -resume -with-docker
```

This time it does indeed work!

```console title="Output"
TODO
```

You should find the cowsay'ed output in the `results` directory.

TODO UPDATE FILENAME AND CONTENT

```console title="results/cowsay-output-Bonjour.txt"
  _______
| Bonjour |
  =======
       \
        \
          ^__^
          (oo)\_______
          (__)\       )\/\
              ||----w |
              ||     ||
```

You see that the character is saying all the greetings.

TODO do we want to note that future course modules (either advanced or side quest) will show how to use a conditional to skip the collect step if we want to emit the cowsay'ed greetings individually, and how to use metadata management to assign a specific character to each greeting?

### 2.7. Inspect how Nextflow launched the containerized task

Let's take a look at the work subdirectory for one of the `cowSay` process calls to get a bit more insight on how Nextflow works with containers under the hood.

Check the output from your `nextflow run` command to find the call ID for the `cowsay` process.
Then navigate to the work subdirectory.
In it, you will find the `.command.run` file that contains all the commands Nextflow ran on your behalf in the course of executing the pipeline.

Open the `.command.run` file and search for `nxf_launch`; you should see something like this:

```bash
nxf_launch() {
    docker run -i --cpu-shares 1024 -e "NXF_TASK_WORKDIR" -v /workspace/gitpod/nf-training/hello-nextflow/work:/workspace/gitpod/nf-training/hello-nextflow/work -w "$NXF_TASK_WORKDIR" --name $NXF_BOXID community.wave.seqera.io/library/pip_cowsay:131d6a1b707a8e65 /bin/bash -ue /workspace/gitpod/nf-training/hello-nextflow/work/8c/738ac55b80e7b6170aa84a68412454/.command.sh
}
```

As you can see, Nextflow is using the `docker run` command to launch the process call.
It also mounts the corresponding work subdirectory into the container, sets the working directory inside the container accordingly, and runs our templated bash script in the `.command.sh` file.
All the hard work we learned about in the previous sections is done for us by Nextflow!

### Takeaway

You know how to use containers in Nextflow to run processes.

### What's next?

Take a break!
When you're ready, move on to Part 6 to learn how to configure the execution of your pipeline to fit your infrastructure as well as manage configuration of inputs and parameters.
