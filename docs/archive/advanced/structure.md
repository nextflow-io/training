---
title: Workflow Structure
description: Advanced Nextflow Training Workshop
---

# Workflow Structure

Nextflow includes a specific directory structure for workflows which can provide some features that can facilitate or enhance your code. In this section we will explore them.

First, let's move into the right directory:

```bash
cd structure
```

There are three directories in a Nextflow workflow repository that have a special purpose:

## `./bin`

The `bin` directory (if it exists) is always added to the `$PATH` for all tasks. If the tasks are performed on a remote machine, the directory is copied across to the new machine before the task begins. This Nextflow feature is designed to make it easy to include accessory scripts directly in the workflow without having to commit those scripts into the container. This feature also ensures that the scripts used inside the workflow move on the same revision schedule as the workflow itself.

It is important to know that Nextflow will take care of updating `$PATH` and ensuring the files are available wherever the task is running, but will not change the permissions of any files in that directory. If a file is called by a task as an executable, the workflow developer must ensure that the file has the correct permissions to be executed.

For example, let's say we have a small R script that produces a csv and a tsv:

```R linenums="1"
#!/usr/bin/env Rscript
library(tidyverse)

plot <- ggplot(mpg, aes(displ, hwy, colour = class)) + geom_point()
mtcars |> write_tsv("cars.tsv")
ggsave("cars.png", plot = plot)
```

We'd like to use this script in a simple workflow:

```groovy linenums="1"
process PlotCars {
    container 'rocker/tidyverse:latest'

    output:
    path("*.png"), emit: plot
    path("*.tsv"), emit: table

    script:
    """
    cars.R
    """
}

workflow {
    PlotCars()

    PlotCars.out.table.view { myfile -> "Found a tsv: $myfile" }
    PlotCars.out.plot.view { myfile -> "Found a png: $myfile" }
}
```

To do this, we can create the bin directory, write our R script into the directory. Finally, and crucially, we make the script executable. This is the code we used to create the `cars.R` script, no need to run it:

```bash
mkdir -p bin
cat << EOF > bin/cars.R
#!/usr/bin/env Rscript
library(tidyverse)

plot <- ggplot(mpg, aes(displ, hwy, colour = class)) + geom_point()
mtcars |> write_tsv("cars.tsv")
ggsave("cars.png", plot = plot)
EOF
chmod +x bin/cars.R
```

!!! warning

    Always ensure that your scripts are executable. The scripts will not be available to your Nextflow processes without this step.

Let's run the script and see what Nextflow is doing for us behind the scenes:

```bash linenums="1"
nextflow run .
```

and then inspect the `.command.run` file that Nextflow has generated

```bash
code work/*/*/.command.run
```

You'll notice a `nxf_container_env` bash function that appends our bin directory to `$PATH`:

```bash linenums="14"
nxf_container_env() {
cat << EOF
export PATH="\$PATH:/workspaces/training/nf-training-advanced/structure/bin"
EOF
}
```

When working on the cloud, Nextflow will also ensure that the bin directory is copied onto the virtual machine running your task in addition to the modification of `$PATH`.

!!! warning

    Always use a portable shebang line in your bin directory scripts.

    In the R script example shown above, I may have the `Rscript` program installed at (for example) `/opt/homebrew/bin/Rscript`. If I hard-code this path into my `cars.R`, everything will work when I'm testing locally outside of the docker container, but will fail when running with docker/singularity or in the cloud as the `Rscript` program may be installed in a different location in those contexts.

    It is __strongly__ recommended to use `#!/usr/bin/env` when setting the shebang for scripts in the `bin` directory to ensure maximum portability.

## `./templates`

If a process script block is becoming too long, it can be moved to a template file. The template file can then be imported into the process script block using the `template` method. This is useful for keeping the process block tidy and readable. Nextflow's use of `$` to indicate variables also allows for directly testing the template file by running it as a script.

The structure directory already contains an example template - a very simple python script. We can add a new process that uses this template:

```groovy linenums="14"
process SayHiTemplate {
    debug true
    input: val(name)
    script: template 'adder.py'
}
```

## `./lib`

In the previous chapter, we saw the addition of small helper Groovy functions to the `main.nf` file. It may at times be helpful to bundle functionality into a new Groovy class. Any classes defined in the `lib` directory are available for use in the workflow - both `main.nf` and any imported modules.

### Making a Metadata Class

!!! note
Using custom Groovy is considered a very advanced use case and you should not need it for the majority of workflows. The language server will complain about this but you can safely ignore it.

!!! exercise

    Create a new class in `./lib/Metadata.groovy` that extends the `HashMap` class and adds a `hi` method.

Let's consider an example where we create our own custom class to handle metadata. We can create a new class in `./lib/Metadata.groovy`. We'll extend the built-in `HashMap` class, and add a simple method to return a value:

```groovy linenums="1"
class Metadata extends HashMap {
    def hi() {
        return "Hello, workshop participants!"
    }
}
```

We can then use this class in our workflow:

```groovy linenums="1"
workflow {
    channel.of("Montreal")
        .map { new Metadata() }
        .view()
}
```

We can use the new `hi` method in the workflow:

```groovy linenums="1" hl_lines="3"
workflow {
    channel.of("Montreal")
        .map { new Metadata() }
        .view { metadata -> metadata.hi() }
}
```

At the moment, the `Metadata` class is not making use of the "Montreal" being passed into the closure. Let's change that by adding a constructor to the class:

```groovy linenums="1" hl_lines="3 7"
class Metadata extends HashMap {
    Metadata(String location) {
        this.location = location
    }

    def hi() {
        return this.location ? "Hello, from ${this.location}!" : "Hello, workshop participants!"
    }
}
```

Which we can use like so:

```groovy linenums="1" hl_lines="3"
workflow {
    channel.of("Montreal")
        .map { place -> new Metadata(place) }
        .view { metadata -> metadata.hi() }
}
```

We can also use this method when passing the object to a process:

```groovy linenums="1"
process UseMeta {
    input: val(meta)
    output: path("out.txt")
    script: "echo '${meta.hi()}' | tee out.txt"
}

workflow {
    place_ch = channel.of("Montreal")
        .map { place -> new Metadata(place) }

    UseMeta(place_ch)
        .view()
}
```

Why might this be helpful? You can add extra classes to the metadata which can be computed from the existing metadata. For example, we might want to add a method to get the adapter prefix into our Metadata class:

```groovy linenums="10"
def getAdapterStart() {
    this.adapter?.substring(0, 3)
}
```

Which we might use like so:

```groovy linenums="1" hl_lines="4 10"
process UseMeta {
    input: val(meta)
    output: path("out.txt")
    script: "echo '${meta.adapter} prefix is ${meta.getAdapterStart()}' | tee out.txt"
}

workflow {
    place_ch = channel.of("Montreal")
        .map { place -> new Metadata(place) }
        .map { metadata -> metadata + [adapter:"AACGTAGCTTGAC"] }

    UseMeta(place_ch)
        .view()
}
```

You might even want to reach out to external services such as a LIMS or the E-utilities API. Here we add a dummy "getSampleName()" method that reaches out to a public API:

```groovy linenums="13"
def getSampleName() {
    def get = new URL('https://postman-echo.com/get?sampleName=Fido').openConnection()
    def getRC = get.getResponseCode();
    if (getRC.equals(200)) {
        JsonSlurper jsonSlurper = new JsonSlurper()
        def json = jsonSlurper.parseText(get.getInputStream().getText())
        return json.args.sampleName
    }
}
```

This relies on jsonSlurper which isn't included by default. Import this by adding the following to the top of the file:

```groovy linenums="1"
import groovy.json.JsonSlurper
```

Which we can use like so:

```groovy linenums="1"
process UseMeta {
    input: val(meta)
    output: path("out.txt")
    script:
    "echo '${meta.adapter} prefix is ${meta.getAdapterStart()} with sampleName ${meta.getSampleName()}' | tee out.txt"
}
```

!!! note "Nextflow caching"

    When we start passing custom classes through the workflow, it's important to understand a little about the Nextflow caching mechanism. When a task is run, a unique hash is calculated based on the task name, the input files/values, and the input parameters. Our class extends from `HashMap`, which means that the hash will be calculated based on the contents of the `HashMap`. If we add a new method to the class, or amend a class method, this does not change the value of the objects in the hash, which means that the hash will not change.

!!! exercise

    Can you show changing a method in our `Metadata` class does not change the hash?

    ??? solution

        We might increase the length of the adapter prefix to 5 characters:

        ```groovy linenums="12"
            def getAdapterStart() {
                this.adapter?.substring(0, 5)
            }
        ```

        Changing this method and resuming the workflow will not change the hash, and the existing method will be used.

We are not limited to using or extending the built-in Groovy classes. Let's start by creating a `Dog` class in `./lib/Dog.groovy`:

```groovy linenums="1"
class Dog {
    String name
    Boolean isHungry = true
}
```

We can create a new dog at the beginning of the workflow:

```groovy linenums="1" hl_lines="3"
workflow {
    def dog = new Dog(name: "fido")
    log.info "Found a new dog: ${dog.name}"
}
```

We can pass objects of our class through channels. Here we take a channel of dog names and create a channel of dogs:

```groovy linenums="1" hl_lines="2-4"
workflow {
    channel.of("Argente", "Absolon", "Chowne")
        .map { name -> new Dog(name: name) }
        .view { dog -> "Found a new dog: ${dog.name}" }
}
```

If we try to use this new class in a resumed process, no caches will be used.

!!! exercise

    Show that the `Dog` class is not cached when resuming a workflow.

### Making a ValueObject

Nextflow has provided a decorator to help serialize your custom classes. By adding `@ValueObject` to the class definition, Nextflow will automatically serialize the class and cache it. This is useful if you want to pass a custom class through a channel, or if you want to use the class in a resumed workflow.

Let's add the decorator to our `Dog` class:

```groovy linenums="1"
import nextflow.io.ValueObject

@ValueObject
class Dog {
    String name
    Boolean isHungry = true
}
```

Lastly, we will need to register the class with Kryo, the Java serialization framework. Again, Nextflow provides a helper method to do this. We can add the following to the `main.nf` file:

```groovy linenums="1" hl_lines="2"
workflow {
    nextflow.util.KryoHelper.register(Dog)
```

!!! exercise

    Show that the `Dog` class can now be used in processes and cached correctly.
