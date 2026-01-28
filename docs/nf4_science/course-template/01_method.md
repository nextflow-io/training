# Part 1: Method overview and manual testing

[High level summary of what the analysis aims to achieve.]

There are multiple valid methods for performing this type of analysis.
For this course, we are following the method described [here](url) by [authors] at [affiliation](url).

Our goal is therefore to develop a workflow that implements the following processing steps: [brief enumeration].

<figure class="excalidraw">
--8<-- "docs/nf4_science/[course-name]/img/[filename]"
</figure>

- **[TOOL]:** [Summary of analysis step]

However, before we dive into writing any workflow code, we are going to try out the commands manually on some test data.
The tools we need are not installed in the GitHub Codespaces environment, so we'll use them via containers.

For more detailed explanations about where to find containers and how they work, see the [Hello Containers](../../hello_nextflow/05_hello_containers.md) chapter in the beginners' course.

[TODO: consider making a note about Seqera containers]

!!! note

    Make sure you're in the `nf4-science/[course-name]` directory.
    The last part of the path shown when you type `pwd` should be `[course-name]`.

---

## 1. [Name of analysis step]

We're going to pull a container image that has [relevant tools] installed, spin it up interactively and run the [analysis] commands on one of the example data files.

### 1.1. Pull the `[tool]` container

First, we use `docker pull` to download the container image.

```bash
docker pull community.wave.seqera.io/library/[container]
```

This gives you the following console output as the system downloads the image:

```console title="Output"
[console output]
```

The image has been downloaded to your computer, but it's not active yet.
We need to actually spin it up (=make it run interactively).

### 1.2. Spin up the container

Let's spin up the container using `docker run -it` to run it interactively.

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/[container]
```

Your prompt will change to something like `(base) root@b645838b3314:/tmp#`, which indicates that you are now _inside_ the container.

The `-v ./data:/data` part of the command will enable us to access the contents of the `data/` directory from inside the container.

```bash
ls /data/[thing]
```

```console title="Output"
[console output]
```

You are now ready to try running [analysis or tool].

### 1.3. Run the [analysis] command

Now, we can run `[tool]` to [do analysis].

[optional explanations about command options]

```bash
[command]
```

This should run very quickly:

```console title="Output"
[console output]
```

[adapt depending on tool behavior] You can find the output files in the same directory as the original data:

```bash
ls /data/[path]
```

```console title="Output"
console output
```

[some explanation of the output; provide contents preview if applicable]

### 1.4. Save the output files

[adapt depending on tool behavior; delete if not necessary]

Output files that were created _inside_ the container will be inaccessible to future work, so let's move these to a new directory.

```bash
mkdir /data/[path]
mv [things]* /data/[path]
```

Now you will be able to use those files as inputs for subsequent analysis steps.

### 1.5. Exit the container

```bash
exit
```

[Concluding statement]

[For the first step: Now we're going to repeat this process for the other steps we plan to implement in our pipeline.]

---

## 2. [Name of second analysis step]

[repeat as needed]

[At final step: call out the final outputs]

---

### Takeaway

You have reviewed the analysis method and tested all the individual commands interactively in the relevant containers.

### What's next?

Learn how to wrap those same commands into a multi-step workflow that uses containers to execute the work.
