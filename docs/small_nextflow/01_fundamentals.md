# Part 1: Fundamentals

In this first part, we'll learn the building blocks of Nextflow by creating channels, defining our first process, and working with parameters and metadata.
We're going to build a workflow that produces a gallery of classified cat images, starting from the very beginning.

---

## Introduction

Nextflow offers you a way to iterate over a collection of files, so let's grab some files to iterate over.
We're going to write a workflow which produces a gallery of good and bad cats.
First things we're going to need are some cats.

We're starting with a (nearly) empty directory.
There is a `.stuff` directory that contains some bits and pieces to help us along during the workshop, but you can imagine that we're essentially starting from scratch.

### Fetch some cat images

The first thing we're going to need is some data.
I've created a small script that pulls from the Cat-As-A-Service API to give us some random cats.

Let's see what the script can do:

```bash
.stuff/cat_me.sh --help
```

To start, let's grab 4 cats.
By default, the script will save the images to `./data/pics`:

```bash
.stuff/cat_me.sh --count 4 --prefix data/pics
```

This will generate some example data.
It will look something like this:

```console title="Directory structure"
data
└── pics
    ├── 5n4MTAC6ld0bVeCe.jpg
    ├── 5n4MTAC6ld0bVeCe.txt
    ├── IOSNx33kkgkiPfaP.jpg
    ├── IOSNx33kkgkiPfaP.txt
    ├── IRuDgdPJZFA39dyf.jpg
    ├── IRuDgdPJZFA39dyf.txt
    ├── uq5KqqiF0qpgTQVA.jpg
    └── uq5KqqiF0qpgTQVA.txt
```

### Create a channel of images

Now let's iterate over those images in Nextflow.
To start, we'll just create a channel of those images.
We're not going to do anything with them, but to make sure that everything is working, we connect the channel to the `view` operator which takes the things in the channel (files in our case) and prints a String representation of those things to the command line.

Open `main.nf` and add the following:

```groovy title="main.nf" linenums="1"
#!/usr/bin/env nextflow

workflow {
    channel.fromPath("data/pics/*.{png,gif,jpg}")
    | view
}
```

Now run the workflow:

```bash
nextflow run main.nf
```

You should see output showing the paths to your cat images.

!!! tip "Glob patterns"

    The `{png,gif,jpg}` syntax is called brace expansion.
    It's a glob pattern that matches any file ending in `.png`, `.gif`, or `.jpg`.
    This is equivalent to writing three separate patterns: `*.png`, `*.gif`, and `*.jpg`.
    Using brace expansion keeps our code concise when we need to match multiple file extensions.

### Takeaway

You now know how to create a channel from files using glob patterns and view its contents.

### What's next?

Let's actually do something with those images by creating our first process.

---

## Channels and processes

Let's try actually doing something with those images.
We'll start with something simple - resizing images.
We'll need to download some software to do so.
Eventually, we'll talk about containers and reproducibility, but just to start, let's download the software to our machine:

```bash
sudo apt update
sudo apt-get install -y imagemagick
```

### Understanding process structure

We'll create a new "process" for our resize operation.
You can think of these processes as templates, or classes.
We'll connect the process to a channel and fire off a new "task" for each thing in the channel.

In our Resize process definition (indeed in most process definitions), there are three blocks - `input`, `output`, and `script`.
We connect these processes by channels and these blocks describe what we expect to get, what we expect to emit, and the work we want to do in-between.

### Create the Resize process

Update your `main.nf` to add a Resize process:

```groovy title="main.nf" linenums="1"
#!/usr/bin/env nextflow

workflow {
    images = channel.fromPath("data/pics/*.{png,gif,jpg}")

    images
    | Resize
    | view
}

process Resize {
    input: path(img)
    output: path("resized-*")
    script: "convert ${img} -resize 400x resized-${img.baseName}.png"
}
```

### Understanding the input block

The input block describes what we expect to take from the input channel.
The "things" in the channel can have a type.
The most common "types" are:

- `val` (values like Strings, integers, Maps, complex objects),
- `path` (paths like directories or files), or
- `tuple` (a collection of values and paths).

We'll see more of those later, but for the moment, the channel we created in the `workflow` block is a channel of files, so in our process definition, we'll say "I'm going to supply a channel of paths to this process, and as our process takes things from that channel to spawn a new task, we'll call the thing in the channel `img`."

### Understanding the output block

The output block describes what we want to emit into the process' output channel.
Again, we can describe the "type" of thing emitted - `val`, `path`, `tuple` and others.
For now, we'll promise to produce a file (or directory) that matches the glob pattern `resized-*`.

### Understanding the script block

The script block describes what work we want to do on each of the things in the channel - how we're going to transform each of the things we pull from the input channel into the files or values we promised to emit into the output channel.

By default, the script block will be rendered into a bash script, but you can use any interpreted language that makes sense to you - python, ruby, R, zsh, closure, whatever.
In this introductory workshop, we'll stick with the default bash.

We run the `magick` command from imagemagick which performs many types of manipulation.
In our case, we'll use the `-resize` argument to resize the image to a width of 400 pixels.
We also supply an output filename.
You'll notice that we use the `${img}` variable twice in our script block.
This `${img}` is the variable we defined in the input block.
For each iteration of our process (each task), the variable will be the path to our individual image.

For example, if the "thing" in the channel is the image `kitten.jpg`, then when Nextflow creates a new Resize task for this file, it will "render" our script block into bash, replacing the `${img}` variables with the path to produce this valid bash:

```bash
convert kitten.jpg -resize 400x resized-kitten.jpg
```

### Run the workflow

Now let's run our workflow!
We'll iterate over all of the images in `data/pics` (relative to our current location) and produce a channel of resized pictures that we then pipe into the `view` operator to print the channel contents to stdout.

```bash
nextflow run main.nf
```

### Takeaway

You now understand the three-part structure of a Nextflow process: input, output, and script blocks work together to transform data flowing through channels.

### What's next?

Let's explore how Nextflow executes each task in isolation.

---

## Investigate task execution

Every task in Nextflow is executed in its own unique work directory.
This directory isolation is a fundamental feature that ensures tasks cannot interfere with each other, even when running in parallel.

### Understanding work directories

The work directory path is calculated by constructing a hash of all task inputs.
This means that if you run the same task with the same inputs, Nextflow will recognize it and can reuse the cached results (we'll explore this with `-resume` later).

Let's explore where Nextflow actually ran our tasks.
Look at the output from your last workflow run - you'll see something like:

```console
executor >  local (4)
[a0/e7b2d4] Resize (1) | 4 of 4 ✔
```

That `a0/e7b2d4` is the hash prefix for the task directory.
Let's explore what's inside:

```bash
tree work
```

You'll see a directory structure like:

```console
work
└── a0
    └── e7b2d4a1f3c8e9b0a7f6d5c4b3a2e1f0
        ├── resized-5n4MTAC6ld0bVeCe.png
        └── ...
```

### Exploring task files

Each work directory contains several hidden files that Nextflow uses to track task execution.
Let's see them all:

```bash
tree -a work
```

Now you'll see additional files:

```console
work
└── a0
    └── e7b2d4a1f3c8e9b0a7f6d5c4b3a2e1f0
        ├── .command.begin
        ├── .command.err
        ├── .command.log
        ├── .command.out
        ├── .command.run
        ├── .command.sh
        ├── .exitcode
        └── resized-5n4MTAC6ld0bVeCe.png
```

These files serve different purposes:

- `.command.sh` - The actual script that was executed (with all variables resolved)
- `.command.run` - The wrapper script that Nextflow uses to execute the task
- `.command.out` - Standard output from the task
- `.command.err` - Standard error from the task
- `.command.log` - Combined stdout and stderr
- `.exitcode` - The exit code from the task (0 = success)
- `.command.begin` - Setup instructions run before the task

The most useful file for debugging is `.command.sh`.
Let's look at one:

```bash
cat work/a0/*/command.sh
```

You'll see the actual bash script that was executed with all Nextflow variables resolved:

```bash
convert 5n4MTAC6ld0bVeCe.jpg -resize 400x resized-5n4MTAC6ld0bVeCe.png
```

### Task isolation and idempotence

This isolation serves two critical purposes:

1. **Independence**: Tasks running in parallel cannot accidentally overwrite each other's files or interfere with each other's execution
2. **Idempotence**: Running the same task with the same inputs will produce the same outputs in the same location

**Idempotence** means that executing a task multiple times with identical inputs produces identical results.
This is crucial for reproducibility and for Nextflow's caching system.
Because the work directory is determined by hashing the inputs, identical inputs always map to the same directory, allowing Nextflow to detect when work can be reused.

### Takeaway

Each task runs in its own isolated directory to prevent interference between parallel executions.

### What's next?

Let's learn some convenient methods for working with file paths.

---

## Harmonization

One of the nice features of the `magick` utility is that it will also do file format conversion for us.
It will infer the format from the extension of the final argument.
For example, if we execute:

```bash
convert kitten.jpg -resize 400x resized-kitten.png
```

The `magick` utility will both resize the image and convert the jpg to png format.
Let's say we want to ensure that downstream in our workflow, we'd like to ensure all images are in the png format.
How might we modify our `script` block to replace the extension or pull out the file basename so that we can append the `.png` extension?

### The bash way (harder)

If you're a bash wizard, you might know that if you have a variable `$myFile` with the path to our file, you can replace the extension with this arcane incantation:

```bash
file=kitten.jpg
convert "$file" -resize 400x "${file%.*}.png"
```

Or perhaps you use the `basename` utility:

```bash
convert "$file" -resize 400x "$(basename "$file" .${file##*.}).png"
```

I love bash, but it's easy to forget this syntax or mistype it.

### The Nextflow way (easier)

Fortunately for us, when inside the script block the `img` variable is not a bash variable - it's a Nextflow variable, and Nextflow provides some convenience methods for operating on those path objects.
The full list is available in the [Nextflow stdlib documentation](https://www.nextflow.io/docs/latest/reference/stdlib-types.html#stdlib-types-path), but one handy method is `baseName`.

We can simply call `${img.baseName}` to return the file base name.
For example:

```groovy title="main.nf" hl_lines="14" linenums="1"
#!/usr/bin/env nextflow

workflow {
    images = channel.fromPath("data/pics/*.{png,gif,jpg}")

    images
    | Resize
    | view
}

process Resize {
    input: path(img)
    output: path("resized-*")
    script: "convert ${img} -resize 400x resized-${img.baseName}.png"
}
```

Update your workflow with this change and run it again:

```bash
nextflow run main.nf
```

### Takeaway

Nextflow path objects provide convenient methods like `baseName` that make file manipulation easier than bash string operations.

### What's next?

Let's make our workflow more flexible by adding parameters.

---

## Parameters

What if we want to make our workflow a little more flexible?
Let's pull out the width and expose it as a parameter to the user.

### Using parameters directly in the script

We could reference a parameter directly in the script block:

```groovy title="main.nf" hl_lines="14" linenums="1"
#!/usr/bin/env nextflow

workflow {
    images = channel.fromPath("data/pics/*.{png,gif,jpg}")

    images
    | Resize
    | view
}

process Resize {
    input: path(img)
    output: path("resized-*")
    script: "convert $img -resize ${params.width}x resized-${img.baseName}.png"
}
```

Now we can run with:

```bash
nextflow run main.nf --width 300
```

### Making inputs explicit (best practice)

This works, but it's considered best practice (and we'll see why in a bit) to make the inputs to a process explicit.
We can do this by adding a second input channel:

```groovy title="main.nf" hl_lines="8 16-17" linenums="1"
#!/usr/bin/env nextflow

workflow {
    images = channel.fromPath("data/pics/*.{png,gif,jpg}")

    Resize(images, params.width)
    | view
}

process Resize {
    input:
      path(img)
      val(width)
    output: path("resized-*")
    script: "convert $img -resize ${width}x resized-${img.baseName}.png"
}
```

The params object still works in the same way:

```bash
nextflow run main.nf --width 500
```

### Takeaway

Exposing parameters as explicit process inputs makes workflows more flexible and clearer about their dependencies.

### What's next?

Let's learn how to attach metadata to our files using tuples and maps.

---

## Extracting an ID

Great, but I'd like a way of retaining the original IDs.

### Understanding the map operator

The `map` operator is one of the most powerful tools in Nextflow.
It takes a collection of items in a channel and transforms them into a new collection of items.
The transformation is defined by a closure - a small piece of code that is evaluated "later" - during workflow execution.
Each item in the new channel is the result of applying the closure to the corresponding item in the original channel.

A closure is written as `{ input -> <expression> }` where to the left of the "stabby operator" `->`, you define the variable used to refer to the closure input, and then an expression or series of expressions. The last expression will be the return value of the closure. For map, the items in the resulting output channel are the collection of values returned by each invocation of the closure.

Let's use `map` to extract the ID from each filename:

```groovy title="main.nf" hl_lines="5-6" linenums="1"
#!/usr/bin/env nextflow

workflow {
    images = channel.fromPath("data/pics/*.{png,gif,jpg}")
    | map { img -> [img.baseName, img] }
    | view
}
```

Run this to see the structure:

```bash
nextflow run main.nf
```

### Using a metadata map

Better - we have the id extracted as a String.
What if we want to add other metadata later?
Let's turn it into a Map:

```groovy title="main.nf" hl_lines="5" linenums="1"
#!/usr/bin/env nextflow

workflow {
    images = channel.fromPath("data/pics/*.{png,gif,jpg}")
    | map { img -> [[id: img.baseName], img] }
    | view
}
```

### Update the process to handle tuples

Now we've changed the "shape" of the items in the channel, so we'll update the downstream process:

```groovy title="main.nf" hl_lines="8 13" linenums="1"
#!/usr/bin/env nextflow

workflow {
    images = channel.fromPath("data/pics/*.{png,gif,jpg}")
    | map { img -> [[id: img.baseName], img] }

    Resize(images, params.width)
    | view
}

process Resize {
    input:
        tuple val(meta), path(img)
        val(width)
    output: 
        tuple val(meta), path("resized-*")
    script: "convert $img -resize ${width}x resized-${img.baseName}.png"
}
```

Run the workflow and view the output:

```bash
nextflow run main.nf
```

### Takeaway

Using tuples with metadata maps allows you to carry important information alongside your files as they flow through the workflow.

### What's next?

Now that we understand the fundamentals, let's build a more complex workflow with classification and grouping operations.
