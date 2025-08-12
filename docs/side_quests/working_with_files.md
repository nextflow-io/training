# Working with Files

Bioinformatics workflows often involve processing large numbers of files. Nextflow provides powerful tools to handle files efficiently, helping you organize and process your data with minimal code.

In this side quest, we'll explore how Nextflow handles files, from basic file operations to more advanced techniques for working with file collections. You'll learn how to extract metadata from filenames - a common requirement in bioinformatics pipelines.

By the end of this side quest, you'll be able to:

- Create Path objects from file path strings using Nextflow's `file()` method
- Get file attributes such as name, extension, and path
- Use channels to automate file handling
- Extract sample metadata from filenames
- Group related files using pattern matching

---

## 0. Warmup

### 0.1. Prerequisites

Before taking on this side quest you should:

- Complete the [Hello Nextflow](../hello_nextflow/README.md) tutorial
- Understand basic Nextflow concepts (processes, channels, operators)
- Basic understanding of using sample-specific data (meta data): [Working with sample-specific data](./metadata.md)

### 0.2. Starting Point

Let's move into the project directory:

```bash
cd side-quests/working_with_files
```

You can set VSCode to focus on this directory:

```bash
code .
```

You'll find a simple workflow file (`file_operations.nf`) and a data directory containing some example files.

```console title="Directory contents"
> tree
.
├── data
│   ├── patientA_rep1_normal_R1_001.fastq.gz
│   ├── patientA_rep1_normal_R2_001.fastq.gz
│   ├── patientA_rep1_tumor_R1_001.fastq.gz
│   ├── patientA_rep1_tumor_R2_001.fastq.gz
│   ├── patientA_rep2_normal_R1_001.fastq.gz
│   ├── patientA_rep2_normal_R2_001.fastq.gz
│   ├── patientA_rep2_tumor_R1_001.fastq.gz
│   ├── patientA_rep2_tumor_R2_001.fastq.gz
│   ├── patientB_rep1_normal_R1_001.fastq.gz
│   ├── patientB_rep1_normal_R2_001.fastq.gz
│   ├── patientB_rep1_tumor_R1_001.fastq.gz
│   ├── patientB_rep1_tumor_R2_001.fastq.gz
│   ├── patientC_rep1_normal_R1_001.fastq.gz
│   ├── patientC_rep1_normal_R2_001.fastq.gz
│   ├── patientC_rep1_tumor_R1_001.fastq.gz
│   └── patientC_rep1_tumor_R2_001.fastq.gz
└── file_operations.nf
```

This directory contains paired-end sequencing data for three patients (A, B, C), with the typical `_R1_` and `_R2_` naming convention for forward and reverse reads. Each patient has normal and tumor tissue types, and patient A has two replicates.

### 0.3. Running the Workflow

Take a look at the workflow file:

```groovy title="file_operations.nf" linenums="1"
workflow {
    // Create a Path object from a string path
    myFile = 'data/patientA_rep1_normal_R1_001.fastq.gz'
    println "${myFile} is of class ${myFile.class}"
}
```

We have a mini-workflow that refers to a single file path in it's workflow, then prints it to the console.

Run the workflow:

```bash
nextflow run file_operations.nf
```

```console title="Starting Output"
 N E X T F L O W   ~  version 24.10.4

Launching `file_operations.nf` [romantic_chandrasekhar] DSL2 - revision: 5a4a89bc3a

data/patientA_rep1_normal_R1_001.fastq.gz is of class class java.lang.String
```

!!! note

    We printed the string path exactly as we wrote it. This is just text output - Nextflow hasn't done anything special with it yet.

---

## 1. Basic File Operations

### 1.1. Creating Path Objects

Let's start by understanding how to create [Path objects](https://www.nextflow.io/docs/latest/reference/stdlib-types.html#path) in Nextflow. In our workflow, we have a string path `data/patientA_rep1_normal_R1_001.fastq.gz`. This is just a plain string - Nextflow doesn't automatically recognize it as representing a file. To work with files properly in Nextflow, we need to convert string paths into proper Path objects using the `file()` method, which provides access to file properties and operations.

Edit the `file_operations.nf` file to include the following:

=== "After"

    ```groovy title="file_operations.nf" linenums="2" hl_lines="2"
    // Create a Path object from a string path
    myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')
    println "${myFile} is of class ${myFile.class}"
    ```

=== "Before"

    ```groovy title="file_operations.nf" linenums="2" hl_lines="2"
    // Create a Path object from a string path
    myFile = 'data/patientA_rep1_normal_R1_001.fastq.gz'
    println "${myFile} is of class ${myFile.class}"
    ```

!!! note

    **What is `.class`?**

    In Groovy (the language Nextflow uses), `.class` tells us what type of object we're working with. It's like asking "what kind of thing is this?" - whether it's a string, a number, a file, or something else. This will help us see the difference between a plain string and a Path object.

Run the workflow:

```bash
nextflow run file_operations.nf
```

```console title="Path object output"
 N E X T F L O W   ~  version 24.10.4

Launching `file_operations.nf` [kickass_coulomb] DSL2 - revision: 5af44b1b59

/Users/jonathan.manning/projects/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz is of class class sun.nio.fs.UnixPath
```

!!! note

    Now we see the full absolute path instead of the relative path we wrote. Nextflow has converted our string into a Path object and resolved it to the actual file location on the system.

    You should only see a single difference, the file path will now be absolute, like `/workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz`.

### 1.2. File Attributes

Why is this helpful? Well now Nextflow understands that `myFile` is a Path object and not just a string, we can access the various attributes of the Path object.

Let's update our workflow to print out the file attributes:

=== "After"

    ```groovy title="file_operations.nf" linenums="2" hl_lines="5-9"
    // Create a Path object from a string path
    myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

    // Print file attributes
    println "File object class: ${myFile.class.name}"
    println "File name: ${myFile.name}"
    println "Simple name: ${myFile.simpleName}"
    println "Extension: ${myFile.extension}"
    println "Parent directory: ${myFile.parent}"
    ```

=== "Before"

    ```groovy title="file_operations.nf" linenums="2" hl_lines="3"
    // Create a file object from a string path
    myFile = file('data/sampleA_rep1_normal_R1_001.fastq.gz')
    println "${myFile} is of class ${myFile.class}"
    ```

Run the workflow:

```bash
nextflow run file_operations.nf
```

You'll see various file attributes printed to the console:

```console title="File Attributes Output"
 N E X T F L O W   ~  version 24.10.4

Launching `file_operations.nf` [ecstatic_ampere] DSL2 - revision: f3fa3dcb48

File object class: sun.nio.fs.UnixPath
File name: patientA_rep1_normal_R1_001.fastq.gz
Simple name: patientA_rep1_normal_R1_001
Extension: gz
Parent directory: /workspaces/training/side-quests/working_with_files/data
```

Notice that the Path object class is `sun.nio.fs.UnixPath` - this is Nextflow's way of representing local files. As we'll see later, remote files will have different class names (like `nextflow.file.http.XPath` for HTTP files), but they all work exactly the same way and can be used identically in your workflows.

!!! note

    **The key difference:**

    - **Path string**: Just text that Nextflow treats as characters
    - **Path object**: A smart file reference that Nextflow can work with

    Think of it like this: a path string is like writing an address on paper, while a Path object is like having a GPS device that knows how to navigate and can tell you details about the journey.

### 1.3. Why proper file handling matters

The difference between strings and Path objects becomes critical when you start building actual workflows with processes. Let's work through this step by step to see what happens.

First, let's create a process that takes a `val` input and see what happens when we pass it a string. We've provided this example in `count_lines.nf`:

```groovy title="count_lines.nf" linenums="2" hl_lines="2 8"
// This will cause problems!
myFile = 'data/patientA_rep1_normal_R1_001.fastq.gz'

process COUNT_LINES {
    debug true

    input:
    val fastq_file

    script:
    """
    echo "Processing file: $fastq_file"
    gzip -dc $fastq_file | wc -l
    """
}

workflow {
    COUNT_LINES(myFile)
}
```

!!! note

    **About `debug = true`:**

    The `debug = true` directive in the process definition causes Nextflow to print the output from your script (like the line count "40") directly in the execution log. Without this, you would only see the process execution status but not the actual output from your script.

    For more information on debugging Nextflow processes, see the [Debugging Nextflow Workflows](../debugging.md) side quest.

Run this workflow to see the error:

```bash
nextflow run count_lines.nf
```

You'll get an error like this:

```console title="Val input with string error"
 N E X T F L O W   ~  version 24.10.2

Launching `file_operations.nf` [wise_poisson] DSL2 - revision: 9b7419747b

executor >  local (1)
[88/f8a197] COUNT_LINES [  0%] 0 of 1
ERROR ~ Error executing process > 'COUNT_LINES'

Caused by:
  Process `COUNT_LINES` terminated with an error exit status (1)

executor >  local (1)
[88/f8a197] COUNT_LINES [100%] 1 of 1, failed: 1 ✘
ERROR ~ Error executing process > 'COUNT_LINES'

Caused by:
  Process `COUNT_LINES` terminated with an error exit status (1)

Command executed:

  echo "Processing file: data/patientA_rep1_normal_R1_001.fastq.gz"
  wc -l data/patientA_rep1_normal_R1_001.fastq.gz

Command exit status:
  1

Command output:
  Processing file: data/patientA_rep1_normal_R1_001.fastq.gz

Command error:
  Processing file: data/patientA_rep1_normal_R1_001.fastq.gz
  wc: data/patientA_rep1_normal_R1_001.fastq.gz: open: No such file or directory

Work dir:
  /workspaces/training/side-quests/working_with_files/work/88/f8a19731c43c834c819ffd146c5713

Tip: you can replicate the issue by changing to the process work dir and entering the command `bash .command.run`
```

**What this error means:**

The process ran successfully, but the bash script failed because the file `data/patientA_rep1_normal_R1_001.fastq.gz` doesn't exist in the process working directory. When you use `val` input, Nextflow passes the string value through to your script, but it doesn't stage the actual file. The process tries to use the string as a file path, but the file isn't there.

Now let's fix this by changing the process to use a `path` input:

=== "After"

    ```groovy title="file_operations.nf" linenums="2" hl_lines="2 8"
    // Now using path input
    myFile = 'data/patientA_rep1_normal_R1_001.fastq.gz'

    process COUNT_LINES {
        debug true

        input:
        path fastq_file

        script:
        """
        echo "Processing file: $fastq_file"
        gzip -dc $fastq_file | wc -l
        """
    }

    workflow {
        COUNT_LINES(myFile)
    }
    ```

=== "Before"

    ```groovy title="file_operations.nf" linenums="2" hl_lines="2 8"
    // This will cause problems!
    myFile = 'data/patientA_rep1_normal_R1_001.fastq.gz'

    process COUNT_LINES {
        debug true

        input:
        val fastq_file

        script:
        """
        echo "Processing file: $fastq_file"
        gzip -dc $fastq_file | wc -l
        """
    }

    workflow {
        COUNT_LINES(myFile)
    }
    ```

Run this updated version and you'll get a different error:

```console title="Path input with string error"
 N E X T F L O W   ~  version 24.10.2

Launching `file_operations.nf` [mighty_poitras] DSL2 - revision: e996edfc53

[-        ] COUNT_LINES -
ERROR ~ Error executing process > 'COUNT_LINES'

Caused by:
  Not a valid path value: 'data/patientA_rep1_normal_R1_001.fastq.gz'

Tip: when you have fixed the problem you can continue the execution adding the option `-resume` to the run command line
```

**What this error means:**

This is much better! Nextflow immediately detected the problem and failed before even starting the process. When you use `path` input, Nextflow validates that you're passing actual file references, not just strings. It's telling you that `'data/patientA_rep1_normal_R1_001.fastq.gz'` is not a valid path value because it's a string, not a Path object.

Now let's fix this properly by using the `file()` method to create a Path object:

=== "After"

    ```groovy title="file_operations.nf" linenums="2" hl_lines="2 8"
    // This works correctly!
    myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

    process COUNT_LINES {
        debug true

        input:
        path fastq_file

        script:
        """
        echo "Processing file: $fastq_file"
        gzip -dc $fastq_file | wc -l
        """
    }

    workflow {
        COUNT_LINES(myFile)
    }
    ```

=== "Before"

    ```groovy title="file_operations.nf" linenums="2" hl_lines="2 8"
    // This will cause problems!
    myFile = 'data/patientA_rep1_normal_R1_001.fastq.gz'

    process COUNT_LINES {
        debug true

        input:
        path fastq_file

        script:
        """
        echo "Processing file: $fastq_file"
        gzip -dc $fastq_file | wc -l
        """
    }

    workflow {
        COUNT_LINES(myFile)
    }
    ```

Now when you run this, it should work correctly! The file will be staged into the process working directory and the `wc -l` command will succeed.

```bash
nextflow run count_lines.nf
```

You should see output like this:

```console title="Successful execution"
 N E X T F L O W   ~  version 24.10.2

Launching `count_lines.nf` [astonishing_tesla] DSL2 - revision: ee38f96485

executor >  local (1)
[8a/3f23d5] COUNT_LINES [100%] 1 of 1 ✔
Processing file: patientA_rep1_normal_R1_001.fastq.gz
      40
```

The process successfully:

- Staged the file into the working directory
- Decompressed the .gz file
- Counted the lines (40 lines in this case)
- Completed without errors

### Takeaway

- **Path strings vs Path objects**: Strings are just text, Path objects are smart file references
- The `file()` method converts a string path into a Path object that Nextflow can work with
- You can access file properties like `name`, `simpleName`, `extension`, and `parent` [using file attributes](https://www.nextflow.io/docs/latest/working-with-files.html#getting-file-attributes)
- Using Path objects instead of strings allows Nextflow to properly manage files in your workflow

**Process Input Outcomes**:
- **`val` input + string**: Process runs but fails at runtime when trying to access non-existent file
- **`path` input + string**: Process fails immediately with validation error (fail fast)
- **`path` input + Path object**: Process succeeds with proper file staging and execution

---

## 2. Using Remote Files

One of the key features of Nextflow is the ability to transparently switch between local files (on the same machine) to remote files accessible over the internet. To do this, all you need to do as a user is switch the file path from a normal file path (e.g. `/path/to/data`) to a file path with a remote protocol at the start. For example, replacing `/path/to/data` with `s3://path/to/data` will switch to using the S3 protocol. Many different protocols are supported:

- HTTP(S)/FTP (http://, https://, ftp://)
- Amazon S3 (s3://)
- Azure Blob Storage (az://)
- Google Cloud Storage (gs://)

To use these, simply replace the string and Nextflow will handle authentication and staging the files to the right place, downloading or uploading and all other file operations you would expect. We call this string the Uniform Resource Identifier (URI).

The key strength of this is we can switch between environments without changing any pipeline logic. For example, you can develop with a small, local test set before moving to a remote set of data by changing the URI.

### 2.1. Using a file from the internet

!!! warning

    Accessing remote data requires an internet connection!

In your workflow, you can replace the string path with an HTTPS one to download this file from the internet. We are going to swap the relative path of the FASTQ files with the remote one stored on the internet. This is the same data as we have been previously using.

=== "After"

    ```groovy title="file_operations.nf" linenums="2" hl_lines="2"
    // Using a remote file from the internet
    myFile = file('https://github.com/nextflow-io/training/blob/bb187e3bfdf4eec2c53b3b08d2b60fdd7003b763/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz')
    println "File object class: ${myFile.class.name}"
    println "${myFile}"
    ```

=== "Before"

    ```groovy title="file_operations.nf" linenums="2" hl_lines="2"
    // Create a file object from a string path
    myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')
    println "${myFile}"
    ```

!!! note

    HTTPS remote data does not accept globs because HTTPS cannot list multiple files, however other storage protocols such as blob storage (`s3://`, `az://`, `gs://`) can.

Run the workflow and it will automatically pull the data from the internet:

```bash
nextflow run file_operations.nf
```

```console title="Remote file output"
 N E X T F L O W   ~  version 24.10.5

Launching `file_operations.nf` [insane_swartz] DSL2 - revision: fff18abe6d

File object class: nextflow.file.http.XPath
https://github.com/nextflow-io/training/blob/bb187e3bfdf4eec2c53b3b08d2b60fdd7003b763/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz
```

In this example very little has changed! This shows how easy it is to switch between local and remote data using Nextflow.

To see the difference, check the working directory located at the hash value of the process. The file will be downloaded to a staging directory located within the work directory, then symlinked into the relevant process directory. If the file was used again, Nextflow would only download it once.

In this way, you can replace local with remote data without changing any pipeline logic.

### Takeaway

- Remote data is accessed using a URI
- Nextflow will automatically download and stage the data to the right place
- Do not write logic to download or upload remote files!

!!! note

    **Note on Object Types**: Notice that local files produce `sun.nio.fs.UnixPath` objects while remote files produce `nextflow.file.http.XPath` objects. Despite these different class names, both work exactly the same way and can be used identically in your workflows. This is a key feature of Nextflow - you can seamlessly switch between local and remote data sources without changing your code logic.

---

## 3. Using Channels for File Handling

While the `file()` method is useful for simple file operations, channels provide a more powerful way to handle multiple files in workflows. We could load our file into a channel using [`Channel.of()`](https://www.nextflow.io/docs/latest/reference/channel.html#of), but we have a much more convenient tool called [`Channel.fromPath()`](https://www.nextflow.io/docs/latest/reference/channel.html#frompath) which finds every file in a directory that matches a given pattern.

### 3.1. Reading Files with Channel.fromPath

Update your `file_operations.nf` file:

=== "After"

    ```groovy title="file_operations.nf" linenums="2" hl_lines="2"
    // Reading files with Channel.fromPath
    ch_fastq = Channel.fromPath('data/patientA_rep1_normal_R1_001.fastq.gz')
    ch_fastq.view { "Found file: $it" }

    // // Print file attributes
    // Comment these out for now, we'll come back to them!
    // println "File name: ${myFile.name}"
    // println "Simple name: ${myFile.simpleName}"
    // println "Extension: ${myFile.extension}"
    // println "Parent directory: ${myFile.parent}"
    ```

=== "Before"

    ```groovy title="file_operations.nf" linenums="2" hl_lines="5-8"
    // Create a Path object from a string path
    myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

    // Print file attributes
    println "File name: ${myFile.name}"
    println "Simple name: ${myFile.simpleName}"
    println "Extension: ${myFile.extension}"
    println "Parent directory: ${myFile.parent}"
    ```

Run the workflow:

```bash
nextflow run file_operations.nf
```

You'll see each file path being emitted as a separate element in the channel:

```console title="Channel.fromPath Output"
 N E X T F L O W   ~  version 24.10.4

Launching `file_operations.nf` [grave_meucci] DSL2 - revision: b09964a583

Found file: /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz

Note how Nextflow has grabbed the file we specified and turned it into a `file` object. `Channel.fromPath()` is a convenient way of creating a new channel populated by a list of files.

### 3.2. Viewing Channel Contents

In our first version, we use `.view()` to print the file name. Let's update our workflow to print out the file attributes:

=== "After"

    ```groovy title="file_operations.nf" linenums="2" hl_lines="3-8"
    // Reading files with Channel.fromPath
    ch_fastq = Channel.fromPath('data/patientA_rep1_normal_R1_001.fastq.gz')
    ch_fastq.view { myFile ->
        println "File object class: ${myFile.class.name}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    }
    ```

=== "Before"

    ```groovy title="file_operations.nf" linenums="2" hl_lines="3"
    // Reading files with Channel.fromPath
    ch_fastq = Channel.fromPath('data/patientA_rep1_normal_R1_001.fastq.gz')
    ch_fastq.view { myFile -> "Found file: $myFile" }

    // // Print file attributes
    // Comment these out for now, we'll come back to them!
    // println "File name: ${myFile.name}"
    // println "Simple name: ${myFile.simpleName}"
    // println "Extension: ${myFile.extension}"
    // println "Parent directory: ${myFile.parent}"
    ```

Run the workflow:

```bash
nextflow run file_operations.nf
```

```console title="Channel.fromPath Output"
 N E X T F L O W   ~  version 24.10.4

Launching `file_operations.nf` [furious_swanson] DSL2 - revision: c35c34950d

File object class: sun.nio.fs.UnixPath
File name: patientA_rep1_normal_R1_001.fastq.gz
Simple name: patientA_rep1_normal_R1_001
Extension: gz
Parent directory: /workspaces/training/side-quests/working_with_files/data
```

### 3.3. Using a glob to match multiple files

`Channel.fromPath()` can take a glob pattern as an argument, which will match all files in the directory that match the pattern. Let's grab both of the pair of FASTQs associated with this patient.

A glob pattern is a pattern that matches one or more characters in a string. The `*` wildcard is the most common glob pattern, which will match any character in it's place. To do this, we replace the full path with a `*` wildcard, which will match any character in it's place. In this case, we will replace the read number from `R1` to `R*`.

=== "After"

    ```groovy title="file_operations.nf" linenums="3"
    ch_fastq = Channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
    ```

=== "Before"

    ```groovy title="file_operations.nf" linenums="3"
    ch_fastq = Channel.fromPath('data/patientA_rep1_normal_R1_001.fastq.gz')
    ```

Run the workflow:

```bash
nextflow run file_operations.nf
```

```console title="Channel.fromPath Glob Output"
 N E X T F L O W   ~  version 24.10.4

Launching `file_operations.nf` [boring_sammet] DSL2 - revision: d2aa789c9a

File object class: sun.nio.fs.UnixPath
File name: patientA_rep1_normal_R1_001.fastq.gz
Simple name: patientA_rep1_normal_R1_001
Extension: gz
Parent directory: /workspaces/training/side-quests/working_with_files/data
File object class: sun.nio.fs.UnixPath
File name: patientA_rep1_normal_R2_001.fastq.gz
Simple name: patientA_rep1_normal_R2_001
Extension: gz
Parent directory: /workspaces/training/side-quests/working_with_files/data
```

Using this method, we could grab as many or as few files as we want just by changing the glob pattern. If we made it more generous, we could grab all the files in the `data` directory, but we'll come back to that later.

### Takeaway

- `Channel.fromPath()` creates a channel with files matching a pattern
- Each file is emitted as a separate element in the channel
- We can use a glob pattern to match multiple files

---

## 4. Extracting Patient Metadata from Filenames

One of the most common tasks in bioinformatics workflows is extracting metadata from filenames. This is especially important when working with sequencing data, where filenames often contain information about the sample, condition, replicate, and read number.

This isn't ideal - metadata should never be embedded in filenames, but it's a common reality. We want to extract that metadata in a standardised manner so we can use it later.

Let's explore how to extract metadata from our FASTQ filenames using Nextflow's powerful data transformation capabilities.

### 4.1. Basic Metadata Extraction

First, let's modify our workflow to extract metadata from the filenames.

First we will grab the simpleName of the file, which includes the metadata, and return with the file. Then, we will separate out the metadata by underscores using tokenize. Finally, we will use string handling to remove additional text like "rep" which aren't required right now.

=== "After"

    ```groovy title="file_operations.nf" linenums="4"
    ch_fastq.map { myFile ->
        [ myFile.simpleName, myFile ]
    }
    .view()
    ```

=== "Before"

    ```groovy title="file_operations.nf" linenums="4"
    ch_fastq.view {
        println "File name: ${it.name}"
        println "Simple name: ${it.simpleName}"
        println "Extension: ${it.extension}"
        println "Parent directory: ${it.parent}"
    }
    ```

```bash
nextflow run file_operations.nf
```

```console title="Sample Metadata Output"
 N E X T F L O W   ~  version 24.10.4

Launching `file_operations.nf` [furious_liskov] DSL2 - revision: dde7b5315e

[patientA_rep1_normal_R1_001, /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz]
[patientA_rep1_normal_R2_001, /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R2_001.fastq.gz]
```

Note how we have separated the patient `simpleName`, which includes the metadata, from the `file` object. This is useful if we want to use the patient metadata in a later process.

### 4.2. Extracting Metadata from Filenames

Our metadata is embedded in the filename, but it's not in a standard format. We need to split up the filename into it's components which are separated by underscores.

Nextflow includes a method called `tokenize()` which is perfect for this task.

=== "After"

    ```groovy title="file_operations.nf" linenums="4" hl_lines="2"
    ch_fastq.map { myFile ->
        [ myFile.simpleName.tokenize('_'), myFile ]
    }
    ```

=== "Before"

    ```groovy title="file_operations.nf" linenums="4" hl_lines="2"
    ch_fastq.map { myFile ->
        [ myFile.simpleName, myFile ]
    }
    ```

Once we run this, we should see the patient metadata as a list of strings, and the Path object as the second element in the tuple.

```bash
nextflow run file_operations.nf
```

```console title="Sample Tokenize Output"
 N E X T F L O W   ~  version 24.10.4

Launching `file_operations.nf` [gigantic_gauss] DSL2 - revision: a39baabb57

[[patientA, rep1, normal, R1, 001], /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz]
[[patientA, rep1, normal, R2, 001], /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R2_001.fastq.gz]
```

Success! We've broken down our patient information from a single string into a list of strings. We can now handle each part of the patient information separately.

### 4.3. Using a map to organise the data

Our meta data is just a flat list at the moment. It's easy to use but hard to read. What is the item at index 3? Can you tell without checking?

A [map](https://www.baeldung.com/groovy-maps) is Groovy's version of a key-value store. Every item has a key and a value and we can refer to each key to get the value. This will make our code much easier to read, i.e. we go from this:

```groovy
data = [patientA, 1, normal, R1]

println data[3]
```

to this:

```groovy
data = [sample: sampleA, replicate: 1, type: normal, readNum: 1]

println data.readNum
```

Let's convert our flat list into a map now.

=== "After"

    ```groovy title="file_operations.nf" linenums="4" hl_lines="3-11"
    ch_fastq.map { myFile ->
        def (sample, replicate, type, readNum) = myFile.simpleName.tokenize('_')
        [
          [
            id: sample,
            replicate: replicate,
            type: type,
            readNum: readNum,
          ],
          myFile
        ]
    }
    ```

=== "Before"

    ```groovy title="file_operations.nf" linenums="4" hl_lines="3"
    ch_fastq.map { myFile ->
        def (sample, replicate, type, readNum) = myFile.simpleName.tokenize('_')
        [ sample, replicate, type, readNum.readNum, myFile ]
    }
    ```

```bash
nextflow run file_operations.nf
```

```console title="Map Output"
 N E X T F L O W   ~  version 24.10.5

Launching `file_operations.nf` [infallible_swartz] DSL2 - revision: 7f4e68c0cb

[[id:patientA, replicate:rep1, type:normal, readNum:R2], /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R2_001.fastq.gz]
[[id:patientA, replicate:rep1, type:normal, readNum:R1], /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz]
```

We have converted our flat list into a map, and now we can refer to each bit of sample data by name instead of by index. This makes our code easier to read and more maintainable.

### 4.4. Simplifying the metadata

The metadata includes a lot of redundant information such as "rep". Let's remove this so replicate is just a number.

We can simplify the metadata by using the `.replace()` method on the replicate string to remove the "rep" prefix, and similarly remove the "R" prefix from the readNum value.

=== "After"

    ```groovy title="file_operations.nf" linenums="4" hl_lines="6 8"
    ch_fastq.map { myFile ->
        def (sample, replicate, type, readNum) = myFile.simpleName.tokenize('_')
        [
          [
            id: sample,
            replicate: replicate.replace('rep', ''),
            type: type,
            readNum: readNum.replace('R', ''),
          ],
          myFile
        ]
    }
    ```

=== "Before"

    ```groovy title="file_operations.nf" linenums="4"
    ch_fastq.map { myFile ->
        def (sample, replicate, type, readNum) = myFile.simpleName.tokenize('_')
        [
          [
            id: sample,
            replicate: replicate,
            type: type,
            readNum: readNum,
          ],
          myFile
        ]
    }
    .view()
    ```

```bash
nextflow run file_operations.nf
```

```console title="Simplified Metadata Output"
 N E X T F L O W   ~  version 24.10.4

Launching `file_operations.nf` [reverent_volta] DSL2 - revision: b3aac71fea

[[sample:sampleA, replicate:1, type:normal, readNum:2], /workspaces/training/side-quests/working_with_files/data/sampleA_rep1_normal_R2_001.fastq.gz]
[[sample:sampleA, replicate:1, type:normal, readNum:1], /workspaces/training/side-quests/working_with_files/data/sampleA_rep1_normal_R1_001.fastq.gz]
```

Nice! we have removed the "rep" from the replicate string.

### Takeaway

- We can handle filenames in Nextflow with the power of a full programming language
- We can treat the filenames as strings to extract relevant information
- Use of methods like `tokenize()` and `replace()` allows us to manipulate strings in the filename

Next up, we will look at how to handle paired-end reads.

---

## 5. Simplifying with Channel.fromFilePairs

Nextflow provides a specialized channel factory method for working with paired files: `Channel.fromFilePairs()`. This method automatically groups files that share a common prefix. This is particularly useful for paired-end sequencing data, where you have two files (e.g., R1 and R2) for each sample.

### 5.1. Basic Usage of fromFilePairs

Complete your `file_operations.nf` file with the following (we will comment out the map operation for now):

=== "After"

    ```groovy title="file_operations.nf" linenums="3" hl_lines="1"
    ch_fastq = Channel.fromFilePairs('data/sampleA_rep1_normal_R{1,2}_001.fastq.gz')
                      .view()
    ```

=== "Before"

    ```groovy title="file_operations.nf" linenums="3" hl_lines="1"
    ch_fastq = Channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
    ch_fastq.map { myFile ->
        def (sample, replicate, type, readNum) = myFile.simpleName.tokenize('_')
        [
            [
                id: sample,
                replicate: replicate.replace('rep', ''),
                type: type,
                readNum: readNum,
            ],
            myFile
        ]
    }
    .view()
    ```

Run the workflow:

```bash
nextflow run file_operations.nf
```

The output will show the paired files grouped together:

```console title="Channel.fromFilePairs Output"
 N E X T F L O W   ~  version 24.10.4

Launching `file_operations.nf` [chaotic_cuvier] DSL2 - revision: 472265a440

[patientA_rep1_normal_R, [/workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz, /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R2_001.fastq.gz]]
```

Note the difference in data structure. Rather than being a list of results, we have a single result in the format `id, [ fastq1, fastq2 ]`. Nextflow has done the hard work of extracting the patient name by examining the shared prefix and using it as a patient id.

### 5.2. Extract metadata from file pairs

We still need the metadata. Our `map` operation from before won't work because it doesn't match the data structure, but we can modify it to work. We already have access to the patient name in the `id` variable, so we can use that to extract the metadata without grabbing the `simpleName` from the Path object like before.

=== "After"

    ```groovy title="file_operations.nf" linenums="3" hl_lines="2-13"
    ch_fastq = Channel.fromFilePairs('data/patientA_rep1_normal_R{1,2}_001.fastq.gz')
    ch_fastq.map { id, fastqs ->
        def (sample, replicate, type, readNum) = id.tokenize('_')
        [
            [
                id: sample,
                replicate: replicate.replace('rep', ''),
                type: type
            ],
            fastqs
        ]
    }
    .view()
    ```

=== "Before"

    ```groovy title="file_operations.nf" linenums="3" hl_lines="3-11"
    ch_fastq = Channel.fromFilePairs('data/sampleA_rep1_normal_R{1,2}_001.fastq.gz')
    .view()
    ```

```bash
nextflow run file_operations.nf
```

```console title="File Pairs Output parsed"

 N E X T F L O W   ~  version 24.10.5

Launching `file_operations.nf` [prickly_stonebraker] DSL2 - revision: f62ab10a3f

[[id:patientA, replicate:1, type:normal], [/workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz, /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R2_001.fastq.gz]]
```

Notice, that this time we don't have a `readNum`. Because read 1 and read 2 are kept together, we do not need to track this in the meta data.

Well done! We have grabbed the metadata from the filenames and used them as values in the tuple.

### Takeaway

- [`Channel.fromFilePairs()` automatically finds and pairs related files](https://www.nextflow.io/docs/latest/reference/channel.html#fromfilepairs)
- The closure defines how to extract patient IDs from filenames
- This simplifies handling paired-end reads in your pipeline

---

## 6. Using File Operations in Processes

Now let's put all this together in a simple process to demonstrate how to use file operations in a Nextflow process.

### 6.1. Create the process

We'll keep it simple and make a process called `ANALYZE_READS` that takes in a tuple of metadata and a pair of fastq files and analyses them. We could imagine this is an alignment, or variant calling or any other step.

Add the following to the top of your `file_operations.nf` file:

=== "After"

    ```groovy title="file_operations.nf - process example" linenums="1" hl_lines="1-23"
    process ANALYZE_READS {
        tag "${meta.id}"

        publishDir "results/${meta.id}", mode: 'copy'

        input:
        tuple val(meta), path(fastqs)

        output:
        tuple val(meta.id), path("${meta.id}_stats.txt")

        script:
        """
        echo "Sample metadata: ${meta.id}" > ${meta.id}_stats.txt
        echo "Replicate: ${meta.replicate}" >> ${meta.id}_stats.txt
        echo "Type: ${meta.type}" >> ${meta.id}_stats.txt
        echo "Read 1: ${fastqs[0]}" >> ${meta.id}_stats.txt
        echo "Read 2: ${fastqs[1]}" >> ${meta.id}_stats.txt
        echo "Read 1 size: \$(gunzip -dc ${fastqs[0]} | wc -l | awk '{print \$1/4}') reads" >> ${meta.id}_stats.txt
        echo "Read 2 size: \$(gunzip -dc ${fastqs[1]} | wc -l | awk '{print \$1/4}') reads" >> ${meta.id}_stats.txt
        """
    }

    workflow {
    ```

=== "Before"

    ```groovy title="file_operations.nf" linenums="1"
    workflow {
    ```

!!! note

    We are calling our map '`meta`'. For a more in-depth introduction to meta maps, see [Patient Specific Data in Workflows](./metadata.md).

### 6.2. Implement the process in the workflow

Then implement the process in the workflow:

=== "After"

    ```groovy title="file_operations.nf" linenums="27" hl_lines="2 13"
        ch_fastq = Channel.fromFilePairs('data/patientA_rep1_normal_R{1,2}_001.fastq.gz')
        ch_samples = ch_fastq.map { id, fastqs ->
            def (sample, replicate, type, readNum) = id.tokenize('_')
            [
                [
                    id: sample,
                    replicate: replicate.replace('rep', ''),
                    type: type
                ],
                fastqs
            ]
        }
        ANALYZE_READS(ch_samples)
    }
    ```

=== "Before"

    ```groovy title="file_operations.nf" linenums="27"
        ch_fastq = Channel.fromFilePairs('data/patientA_rep1_normal_R{1,2}_001.fastq.gz')
        ch_fastq.map { id, fastqs ->
            def (sample, replicate, type, readNum) = id.tokenize('_')
            [
                [
                    id: sample,
                    replicate: replicate.replace('rep', ''),
                    type: type
                ],
                fastqs
            ]
        }
        .view()
    }
    ```

```bash
nextflow run file_operations.nf
```

```console title="ANALYZE_READS Output"
 N E X T F L O W   ~  version 24.10.5

Launching `./file_operations.nf` [shrivelled_cori] DSL2 - revision: b546a31769

executor >  local (1)
[b5/110360] process > ANALYZE_READS (patientA) [100%] 1 of 1 ✔
```

We should see the following files in the `results/patientA` directory:

```console title="Results Directory"
> tree results/patientA
results/patientA
└── patientA_stats.txt
```

The process took our inputs and created a new file with the patient metadata. Based on what you learned in hello-nextflow, what occurred in the working directory?

### 6.3. Include many more patients

Remember Channel.fromPath() accepts a _glob_ as input, which means it can accept any number of files that match the pattern. Therefore if we want to include all the patients we can just modify the input string to include more patients.

=== "After"

    ```groovy title="file_operations.nf" linenums="27"
    ch_fastq = Channel.fromFilePairs('data/*_R{1,2}_001.fastq.gz')
    ```

=== "Before"

    ```groovy title="file_operations.nf" linenums="27"
        ch_fastq = Channel.fromFilePairs('data/patientA_rep1_normal_R{1,2}_001.fastq.gz')
    ```

Run the pipeline now and see all the results:

```bash
nextflow run file_operations.nf
```

```console title="ANALYZE_READS Multiple Samples"
 N E X T F L O W   ~  version 24.10.5

Launching `./file_operations.nf` [big_stonebraker] DSL2 - revision: f7f9b8a76c

executor >  local (8)
[d5/441891] process > ANALYZE_READS (patientC) [100%] 8 of 8 ✔
```

Check the results directory now:

```console title="Results Directory"
> tree results
results
├── patientA
│   └── patientA_stats.txt
├── patientB
│   └── patientB_stats.txt
└── patientC
    └── patientC_stats.txt
```

See how we have analyzed all the patients in one go!

Wait, we have a problem. We have 2 replicates for patientA, but only 1 output file! We are overwriting the output file each time.

### 6.4. Make the published files unique

Since we have access to the patient metadata, we can use it to make the output files unique.

=== "After"

    ```groovy title="file_operations.nf" linenums="4"
    publishDir "results/${meta.type}/${meta.id}/${meta.replicate}", mode: 'copy'
    ```

=== "Before"

    ```groovy title="file_operations.nf" linenums="4"
    publishDir "results/${id}", mode: 'copy'
    ```

We have grabbed the metadata from the patients and used it to construct an output directory for each patient.

Run the pipeline now and see all the results. Remove the results directory first to give yourself a clean workspace:

```bash
rm -r results
nextflow run file_operations.nf
```

```console title="Results Directory"
 N E X T F L O W   ~  version 24.10.5

Launching `./file_operations.nf` [insane_swartz] DSL2 - revision: fff18abe6d

executor >  local (8)
[e3/449081] process > ANALYZE_READS (patientC) [100%] 8 of 8 ✔
```

Check the results directory now:

```console title="Results Directory"
> tree results
results/
├── normal
│   ├── patientA
│   │   ├── 1
│   │   │   └── patientA_stats.txt
│   │   └── 2
│   │       └── patientA_stats.txt
│   ├── patientB
│   │   └── 1
│   │       └── patientB_stats.txt
│   └── patientC
│       └── 1
│           └── patientC_stats.txt
└── tumor
    ├── patientA
    │   ├── 1
    │   │   └── patientA_stats.txt
    │   └── 2
    │       └── patientA_stats.txt
    ├── patientB
    │   └── 1
    │       └── patientB_stats.txt
    └── patientC
        └── 1
            └── patientC_stats.txt
```

See how using patient metadata as values gives us powerful flexibility in our pipeline. By propagating metadata alongside our data in tuples, we can:

1. Create organized output directories based on patient attributes
2. Make decisions in processes based on patient properties
3. Split, join, and recombine data based on metadata values

This pattern of keeping metadata explicit and attached to the data (rather than encoded in filenames) is a core best practice in Nextflow that enables building robust, maintainable bioinformatics workflows.

### Takeaway

- Processes can take files as input and produce new files as output
- The `publishDir` directive organizes outputs based on metadata values
- Metadata in tuples enables structured organization of results
- This approach creates maintainable workflows with clear data provenance



## Summary

In this side quest, you've learned how to work with files in Nextflow, from basic operations to more advanced techniques for handling collections of files. Here's a summary of what we covered:

1. **Basic File Operations**: We created Path objects with `file()` and accessed file attributes like name, path, and extension.

2. **Using Remote Files**: We learned how to transparently switch between local and remote files using URIs, demonstrating Nextflow's ability to handle files from various sources.

3. **Using Channels for File Handling**: We created channels from file patterns with `Channel.fromPath()` and viewed their file attributes with `.view()`.

4. **Extracting Patient Metadata from Filenames**: We used `tokenize()` to extract patient metadata from filenames and use it in our workflow.

5. **Simplifying with Channel.fromFilePairs**: We used `Channel.fromFilePairs()` to automatically pair files with a common prefix and use it in our workflow.

6. **Using File Operations in Processes**: We used `publishDir` to organize outputs based on patient metadata and `tuple` to pass metadata alongside data in our workflow.

These techniques will help you build more efficient and maintainable workflows, especially when working with large numbers of files with complex naming conventions.

### Key Concepts

- **Path Object Creation**

  ```groovy
  // Create a Path object from a string path
  myFile = file('path/to/file.txt')
  ```

- **File Attributes**

  ```groovy
  // Get file attributes
  println myFile.name       // file.txt
  println myFile.baseName   // file
  println myFile.extension  // txt
  println myFile.parent     // path/to
  ```

- **Channel Creation from Files**

  ```groovy
  // Create a channel from a file pattern
  ch_fastq = Channel.fromPath('data/*.fastq.gz')

  // Create a channel from paired files
  ch_pairs = Channel.fromFilePairs('data/*_R{1,2}_001.fastq.gz')
  ```

- **Extracting Metadata**

  ```groovy
  // Extract metadata with tokenize
  def name = file.name.tokenize('_')
  def patientId = name[0]
  def replicate = name[1].replace('rep', '')
  def type = name[2]
  def readNum = name[3]
  ```

- **Using remote files**

  ```groovy
  // Use a local file
  myFile = file('path/to/file.txt')

  // Use a file on FTP
  myFile = file('ftp://path/to/file.txt')

  // Use a file on HTTPS
  myFile = file('https://path/to/file.txt')

  // Use a file on S3
  myFile = file('s3://path/to/file.txt')

  // Use a file on Azure Blob Storage
  myFile = file('az://path/to/file.txt')

  // Use a file on Google Cloud Storage
  myFile = file('gs://path/to/file.txt')
  ```

## Resources

- [Nextflow Documentation: Working with Files](https://www.nextflow.io/docs/latest/working-with-files.html)
- [Channel.fromPath](https://www.nextflow.io/docs/latest/channel.html#frompath)
- [Channel.fromFilePairs](https://www.nextflow.io/docs/latest/channel.html#fromfilepairs)
