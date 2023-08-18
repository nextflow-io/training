---
title: Workflow Structure
description: Advanced Nextflow Training Workshop
---

# Workflow Structure

There are three directories in a Nextflow workflow repository that have a special purpose:

## `./bin`
The `bin` directory (if it exists) is always added to the `$PATH` for all tasks. If the tasks are performed on a remote machine, the directory is copied across to the new machine before the task begins. It is important to know that Nextflow will take care of updating `$PATH` and ensuring the files are available wherever the task is running, but will not change the permissions of any files in that directory. If a file is called by a task as an executable, the workflow developer must ensure that the file has the correct permissions to be executed.

<!-- TODO: Add example of adding a script to the bin directory and using it in the workflow -->

## `./templates`
If a process script block is becoming too long, it can be moved to a template file. The template file can then be imported into the process script block using the `template` method. This is useful for keeping the process block tidy and readable. Nextflow's use of `$` to indicate variables also allows for directly testing the template file by running it as a script.

The chapter_05_structure directory already contains an example template - a very simple python script. We can add a new process that uses this template:

```groovy
process SayHiTemplate {
    debug true
    input: val(name)
    script: template 'adder.py'    
}
```

## `./lib`
In the previous chapter, we saw the addition of small helper Groovy functions to the `main.nf` file. It may at times be helpful to bundle functionality into a new Groovy class. Any classes defined in the `lib` directory are available for use in the workflow - both `main.nf` and any imported modules. 

Classes defined in `lib` can be used for a variety of purposes. For example, the [nf-core/rnaseq](https://github.com/nf-core/rnaseq/tree/master/lib) workflow uses five custom classes:
1) `NfcoreSchema.groovy` for parsing the schema.json file and validating the workflow parameters.
2) `NfcoreTemplate.groovy` for email templating and nf-core utility functions.
3) `Utils.groovy` for provision of a single `checkCondaChannels` method.
4) `WorkflowMain.groovy` for workflow setup and to call the `NfcoreTemplate` class.
5) `WorkflowRnaseq.groovy` for the workflow-specific functions.

The classes listed above all provide utility executed at the beginning of a workflow, and are generally used to "set up" the workflow. However, classes defined in `lib` can also be used to provide functionality to the workflow itself.

### Making a Metadata Class

Let's consider an example where we create our own custom class to handle metadata. We can create a new class in `./lib/Metadata.groovy`. We'll extend the built-in `HashMap` class, and add a simple method to return a value:

```groovy
class Metadata extends HashMap {
    def hi() {
        return "Hello, workshop participants!"
    }
}
```

We can then use this class in our workflow:

```groovy
workflow {
    Channel.of("Montreal") 
    | map { new Metadata() } 
    | view
}
```

We can use the new `hi` method in the workflow:

```groovy 
workflow {
    Channel.of("Montreal") 
    | map { new Metadata() } 
    | view { it.hi() }
}
```

At the moment, the `Metadata` class is not making use of the "Montreal" being passed into the closure. Let's change that by adding a constructor to the class:

```groovy
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

```groovy
workflow {
    Channel.of("Montreal") 
    | map { place -> new Metadata(place) } 
    | view { it.hi() }
}
```

We can also use this method when passing the object to a process:

```groovy
process UseMeta {
    input: val(meta)
    output: path("out.txt")
    script: "echo '${meta.hi()}' | tee out.txt"
}

workflow {   
    Channel.of("Montreal")
    | map { place -> new Metadata(place) } 
    | UseMeta
    | view
}
```

Why might this be helpful? You can add extra classes to the metadata which can be computed from the existing metadata. For exmaple, we might want want to grab the adapter prefix:

```groovy
def getAdapterStart() {
    this.adapter?.substring(0, 3)
}
```

Which we might use like so:

```groovy
process UseMeta {
    input: val(meta)
    output: path("out.txt")
    script: "echo '${meta.adapter} prefix is ${meta.getAdapterStart()}' | tee out.txt"
}

workflow {
    Channel.of("Montreal") 
    | map { place -> new Metadata(place) } 
    | map { it + [adapter:"AACGTAGCTTGAC"] } 
    | UseMeta
    | view
}
```

You might even want to reach out to external services such as a LIMS or e-utilis API. Here we add a dummy "getSampleName()" method that reaches out to a public API:

```groovy
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

Which we can use like so:

```groovy
process UseMeta {
    input: val(meta)
    output: path("out.txt")
    script: 
    "echo '${meta.adapter} prefix is ${meta.getAdapterStart()} with sampleName ${meta.getSampleName()}' | tee out.txt"
}
```

!!! note
    **Nextflow Caching**

    When we start passing custom classes through the workflow, it's important to understand a little about the Nextflow caching mechanism. When a task is run, a unique hash is calculated based on the task name, the input files/values, and the input parameters. Our class extends from `HashMap`, which means that the hash will be calculated based on the contents of the `HashMap`. If we add a new method to the class, or ammend a class method, this does not change the value of the objects in the hash, which means that the hash will not change.


!!! exercise

    Can you show changing a method in our `Metadata` class does not change the hash?

    !!! solution

        **Show answer**

        We might increase the length of the adapter prefix to 5 characters:

        ```{groovy}
            def getAdapterStart() {
                this.adapter?.substring(0, 5)
            }
        ```

        Changing this method and resuming the workflow will not change the hash, and the existing method will be used.

We are not limited to using or extending the built-in Groovy classes. Let's start by creating a `Dog` class in `./lib/Dog.groovy`:

```groovy
class Dog {
    String name
    Boolean isHungry = true
}
```

We can create a new dog at the beginning of the workflow:

```groovy
workflow {
    dog = new Dog(name: "fido")
    log.info "Found a new dog: $dog"
}
```

We can pass objects of our class through channels. Here we take a channel of dog names and create a channel of dogs:

```groovy
workflow {
    Channel.of("Argente", "Absolon", "Chowne")
    | map { new Dog(name: it) }
    | view
}
```

If we try to use this new class in a resumed process, no caches will be used.

!!! exercise

    Show that the `Dog` class is not cached when resuming a workflow.

### Making a ValueObject

Nextflow has provided a decorator to help serialize your custom classes. By adding `@ValueObject` to the class definition, Nextflow will automatically serialize the class and cache it. This is useful if you want to pass a custom class through a channel, or if you want to use the class in a resumed workflow.

Let's add the decorator to our `Dog` class:

```groovy
import nextflow.io.ValueObject

@ValueObject
class Dog {
    String name
    Boolean isHungry = true
}
```

Lastly, we will need to register the class with Kyro, the Java serialization framework. Again, Nextflow provides a helper method to do this. We can add the following to the `main.nf` file:

```groovy
import nextflow.util.KryoHelper

KryoHelper.register(Dog)
```

!!! exercise

    Show that the `Dog` class can now be used in processes and cached correctly.

