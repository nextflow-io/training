# File input processing

Scientific analysis workflows often involve processing large numbers of files.
Nextflow provides powerful tools to handle files efficiently, helping you organize and process your data with minimal code.

### Learning goals

In this side quest, we'll explore how Nextflow handles files, from basic file operations to more advanced techniques for working with file collections.
You'll learn how to extract metadata from filenames, which is a common requirement in scientific analysis pipelines.

By the end of this side quest, you'll be able to:

- Create Path objects from file path strings using Nextflow's `file()` method
- Access file attributes such as name, extension, and parent directory
- Handle both local and remote files transparently using URIs
- Use channels to automate file handling with `channel.fromPath()` and `channel.fromFilePairs()`
- Extract and structure metadata from filenames using string manipulation
- Group related files using pattern matching and glob expressions
- Integrate file operations into Nextflow processes with proper input handling
- Organize process outputs using metadata-driven directory structures

These skills will help you build workflows that can handle different kinds of file inputs with great flexibility.

### Prerequisites

Before taking on this side quest, you should:

- Have completed the [Hello Nextflow](../../hello_nextflow/) tutorial or equivalent beginner's course.
- Be comfortable using basic Nextflow concepts and mechanisms (processes, channels, operators)

<!-- I removed the suggestion to do the metamaps SQ first because that works more naturally after -->

---

## 0. Get started

#### Open the training codespace

If you haven't yet done so, make sure to open the training environment as described in the [Environment Setup](../envsetup/index.md).

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Move into the project directory

Let's move into the directory where the files for this tutorial are located.

```bash
cd side-quests/working_with_files
```

You can set VSCode to focus on this directory:

```bash
code .
```

#### Review the materials

You'll find a simple workflow file called `main.nf`, a `modules` directory containing two module files, and a `data` directory containing some example data files.

??? abstract "Directory contents"

    ```console
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
    ├── main.nf
    └── modules
        ├── analyze_reads.nf
        └── count_lines.nf
    ```

This directory contains paired-end sequencing data from three patients (A, B, C).

For each patient, we have samples that are of type `tumor` (typically originating from tumor biopsies) or `normal` (taken from healthy tissue or blood).
If you're not familiar with cancer analysis, just know that this corresponds to an experimental model that uses paired tumor/normal samples to perform contrastive analyses.

For patient A specifically, we have two sets of technical replicates (repeats).

The sequencing data files are named with a typical `_R1_` and `_R2_` convention for what are known as 'forward reads' and 'reverse reads'.

_Don't worry if you're not familiar with this experimental design, it's not critical for understanding this tutorial._

#### Review the assignment

Your challenge is to write a Nextflow workflow that will parse the samplesheet, extract basic metadata from the file naming structure, and use that metadata to organize the analysis and outputs appropriately.

<!-- TODO: give a bit more details, similar to how it's done in the Metadata side quest -->

#### Readiness checklist

Think you're ready to dive in?

- [ ] I understand the goal of this course and its prerequisites
- [ ] My codespace is up and running
- [ ] I've set my working directory appropriately
- [ ] I understand the assignment

If you can check all the boxes, you're good to go.

---

## 1. Basic file operations

### 1.1. Identify the type of an object with `.class`

Take a look at the workflow file `main.nf`:

```groovy title="main.nf" linenums="1"
#!/usr/bin/env nextflow

workflow {

    // Create a Path object from a string path
    myFile = 'data/patientA_rep1_normal_R1_001.fastq.gz'

    println "${myFile} is of class ${myFile.class}"
}
```

This is a mini-workflow (without any processes) that refers to a single file path in its workflow, then prints it to the console, along with its class.

??? info "What is `.class`?"

    In Nextflow, `.class` tells us what type of object we're working with. It's like asking "what kind of thing is this?" to find out whether it's a string, a number, a file, or something else.
    This will help us illustrate the difference between a plain string and a Path object in the next sections.

Let's run the workflow:

```bash
nextflow run main.nf
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [romantic_chandrasekhar] DSL2 - revision: 5a4a89bc3a

    data/patientA_rep1_normal_R1_001.fastq.gz is of class java.lang.String
    ```

As you can see, Nextflow printed the string path exactly as we wrote it.

This is just text output; Nextflow hasn't done anything special with it yet.
We've also confirmed that as far as Nextflow is concerned, this is only a string (of class `java.lang.String`).
That makes sense, since we haven't yet told Nextflow that it corresponds to a file.

### 1.2. Create a Path object with file()

We can tell Nextflow how to handle files by creating [Path objects](https://www.nextflow.io/docs/latest/reference/stdlib-types.html#path) from path strings.

In our workflow, we can convert the string path `data/patientA_rep1_normal_R1_001.fastq.gz` to a Path object using the `file()` method, which provides access to file properties and operations.

Edit the `main.nf` to wrap the string with `file()` as follows:

=== "After"

    ```groovy title="main.nf" linenums="5" hl_lines="2"
        // Create a Path object from a string path
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        println "${myFile} is of class ${myFile.class}"
    ```

=== "Before"

    ```groovy title="main.nf" linenums="5" hl_lines="2"
        // Create a Path object from a string path
        myFile = 'data/patientA_rep1_normal_R1_001.fastq.gz'

        println "${myFile} is of class ${myFile.class}"
    ```

Now run the workflow again:

```bash
nextflow run main.nf
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [kickass_coulomb] DSL2 - revision: 5af44b1b59

    /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz is of class class sun.nio.fs.UnixPath
    ```

This time, you see the full absolute path instead of the relative path we provided as input.

Nextflow has converted our string into a Path object and resolved it to the actual file location on the system.
The file path will now be absolute, as in `/workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz`.

Notice also that the Path object class is `sun.nio.fs.UnixPath`: this is Nextflow's way of representing local files.
As we'll see later, remote files will have different class names (such as `nextflow.file.http.XPath` for HTTP files), but they all work exactly the same way and can be used identically in your workflows.

!!! tip

    **The key difference:**

    - **Path string**: Just text that Nextflow treats as characters
    - **Path object**: A smart file reference that Nextflow can work with

    Think of it like this: a path string is like writing an address on paper, while a Path object is like having the address loaded in a GPS device that knows how to navigate to there and can tell you details about the journey.

### 1.3. Access file attributes

Why is this helpful? Well, now that Nextflow understands that `myFile` is a Path object and not just a string, we can access the various attributes of the Path object.

Let's update our workflow to print out the built-in file attributes:

=== "After"

    ```groovy title="main.nf" linenums="5" hl_lines="4-9"
        // Create a Path object from a string path
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Print file attributes
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    ```

=== "Before"

    ```groovy title="main.nf" linenums="5" hl_lines="4"
        // Create a Path object from a string path
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        println "${myFile} is of class ${myFile.class}"
    ```

Run the workflow:

```bash
nextflow run main.nf
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [ecstatic_ampere] DSL2 - revision: f3fa3dcb48

    File object class: sun.nio.fs.UnixPath
    File name: patientA_rep1_normal_R1_001.fastq.gz
    Simple name: patientA_rep1_normal_R1_001
    Extension: gz
    Parent directory: /workspaces/training/side-quests/working_with_files/data
    ```

You see the various file attributes printed to the console above.

### 1.4. Feed the file to a process

The difference between strings and Path objects becomes critical when you start building actual workflows with processes.
So far we've verified that Nextflow is now treating out input file as a file, but let's see if we can actually run something on that file in a process.

#### 1.4.1. Import the process and examine the code

We provide you with a pre-written process module called `COUNT_LINES` that takes a file input and counts how many lines are in it.

To use the process in the workflow, you just need to add an include statement before the workflow block:

=== "After"

    ```groovy title="main.nf" linenums="1" hl_lines="3"
    #!/usr/bin/env nextflow

    include { COUNT_LINES } from './modules/count_lines.nf'

    workflow {
    ```

=== "Before"

    ```groovy title="main.nf" linenums="1"
    #!/usr/bin/env nextflow

    workflow {
    ```

You can open the module file to examine its code:

```groovy title="modules/langid.nf" linenums="1" hl_lines="9 12"
#!/usr/bin/env nextflow

process COUNT_LINES {
    debug true

    input:
    val input_file

    script:
    """
    set -o pipefail
    echo "Processing file: $input_file"
    gzip -dc $input_file | wc -l
    """
}
```

As you can see, it's a fairly straightforward little script that unzips the file and counts how many lines it contains.

??? info "What does `debug true` do?"

    The `debug true` directive in the process definition causes Nextflow to print the output from your script (like the line count "40") directly in the execution log.
    Without this, you would only see the process execution status but not the actual output from your script.

    For more information on debugging Nextflow processes, see the [Debugging Nextflow Workflows](debugging.md) side quest.

#### 1.4.2. Add a call to `COUNT_LINES`

Now that the process is available to the workflow, we can add a call to the `COUNT_LINES` process to run it on the input file.

Make the following edits to the workflow:

=== "After"

    ```groovy title="main.nf" linenums="7" hl_lines="11-12"
        // Create a Path object from a string path
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Print file attributes
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"

        // Count the lines in the file
        COUNT_LINES(myFile)
    ```

=== "Before"

    ```groovy title="main.nf" linenums="7" hl_lines="4-9"
        // Create a Path object from a string path
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Print file attributes
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    ```

And now run the workflow:

```bash
nextflow run main.nf
```

??? success "Command output"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [cheeky_hypatia] DSL2 - revision: 281d13c414

    File object class: class sun.nio.fs.UnixPath
    File name: patientA_rep1_normal_R1_001.fastq.gz
    Simple name: patientA_rep1_normal_R1_001
    Extension: gz
    Parent directory: /workspaces/training/side-quests/working_with_files/data
    executor >  local (1)
    [e9/341c05] COUNT_LINES [100%] 1 of 1 ✔
    Processing file: /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz
    40
    ```

This shows we are able to operate on the file appropriately inside a process.

Specifically, Nextflow carried out the following operations successfully:

- Staged the file into the working directory
- Decompressed the .gz file
- Counted the lines (40 lines in this case)
- Completed without error

The key to this smooth operation is that we're explicitly telling Nextflow that our input is a file and should be treated as such.

### 1.5. Troubleshoot basic file input errors

This often trips up newcomers to Nextflow, so let's take a few minutes to look at what happens when you do it wrong.

There are two main places where you can get the file handling wrong: at the level of the workflow, and at the level of the process.

#### 1.5.1. Workflow-level error

Let's see what happens if we revert to treating the file as a string when we specify the input in the workflow block.

Make the following edits to the workflow, making sure to comment out the path-specific print statements:

=== "After"

    ```groovy title="main.nf" linenums="7" hl_lines="2 6-11"
        // Create a Path object from a string path
        myFile = 'data/patientA_rep1_normal_R1_001.fastq.gz'

        // Print file attributes
        println "File object class: ${myFile.class}"
        /*
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
        */

        // Count the lines in the file
        COUNT_LINES(myFile)
    ```

=== "Before"

    ```groovy title="main.nf" linenums="7" hl_lines="4-9"
        // Create a Path object from a string path
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Print file attributes
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"

        // Count the lines in the file
        COUNT_LINES(myFile)
    ```

And now run the workflow:

```bash
nextflow run main.nf
```

??? failure "Command output"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [friendly_goodall] DSL2 - revision: ae50609b20

    [-        ] COUNT_LINES -
    ERROR ~ Error executing process > 'COUNT_LINES'

    Caused by:
      Not a valid path value: 'data/patientA_rep1_normal_R1_001.fastq.gz'



    Tip: view the complete command output by changing to the process work dir and entering the command `cat .command.out`

    -- Check '.nextflow.log' file for details
    ```

This is the important bit:

```console
Not a valid path value: 'data/patientA_rep1_normal_R1_001.fastq.gz'
```

When you specify a `path` input, Nextflow validates that you're passing actual file references, not just strings.
This error is telling you that `'data/patientA_rep1_normal_R1_001.fastq.gz'` is not a valid path value because it's a string, not a Path object.

Nextflow immediately detected the problem and stopped before even starting the process.

#### 1.5.2. Process-level error

The other place we might forget to specify we want Nextflow to treat the input as a file is in the process definition.

Make the following edit to the module:

=== "After"

    ```groovy title="modules/count_lines.nf" linenums="3" hl_lines="7"
    process COUNT_LINES {
        debug true

        input:
        val input_file
    ```

=== "Before"

    ```groovy title="modules/count_lines.nf" linenums="3" hl_lines="7"
    process COUNT_LINES {
        debug true

        input:
        path input_file
    ```

And now run the workflow again:

```bash
nextflow run main.nf
```

??? failure "Command output"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [soggy_golick] DSL2 - revision: ae50609b20

    executor >  local (1)
    [b3/b3023c] COUNT_LINES [  0%] 0 of 1 ✘
    ERROR ~ Error executing process > 'COUNT_LINES'

    Caused by:
      Process `COUNT_LINES` terminated with an error exit status (1)


    Command executed:

      set -o pipefail
      echo "Processing file: data/patientA_rep1_normal_R1_001.fastq.gz"
      gzip -dc data/patientA_rep1_normal_R1_001.fastq.gz | wc -l

    Command exit status:
      1

    Command output:
      Processing file: data/patientA_rep1_normal_R1_001.fastq.gz
      0

    Command error:
      Processing file: data/patientA_rep1_normal_R1_001.fastq.gz
      gzip: data/patientA_rep1_normal_R1_001.fastq.gz: No such file or directory
      0

    Work dir:
      /workspaces/training/side-quests/working_with_files/work/b3/b3023cb2ccb986851301d8e369e79f

    Tip: when you have fixed the problem you can continue the execution adding the option `-resume` to the run command line

    -- Check '.nextflow.log' file for details
    ```

This shows a lot of details about the error because the process is set to output debugging information, as noted above.

These are the most relevant sections:

```console
Command executed:

  set -o pipefail
  echo "Processing file: data/patientA_rep1_normal_R1_001.fastq.gz"
  gzip -dc data/patientA_rep1_normal_R1_001.fastq.gz | wc -l
```

```console
Command error:
  Processing file: data/patientA_rep1_normal_R1_001.fastq.gz
  gzip: data/patientA_rep1_normal_R1_001.fastq.gz: No such file or directory
  0
```

This says the system couldn't find the file; however if you look up the path, there is a file by that name in that location.

When we ran this, Nextflow passed the string value through to the script, but it didn't _stage_ the actual file in the working directory.
So the process tried to use the relative string, `data/patientA_rep1_normal_R1_001.fastq.gz`, but that file doesn't exist within the process working directory.

Taken together, these two examples show you how important it is to tell Nextflow if an input should be handled as a file.

!!! note

    Make sure to go back and fix both intentional errors before continuing to the next section.

### Takeaway

- Path strings vs Path objects: Strings are just text, Path objects are smart file references
- The `file()` method converts a string path into a Path object that Nextflow can work with
- You can access file properties like `name`, `simpleName`, `extension`, and `parent` [using file attributes](https://www.nextflow.io/docs/latest/working-with-files.html#getting-file-attributes)
- Using Path objects instead of strings allows Nextflow to properly manage files in your workflow
- Process Input Outcomes: Proper file handling requires Path objects, not strings, to ensure files are correctly staged and accessible for use by processes.

---

## 2. Using remote files

One of the key features of Nextflow is the ability to switch seamlessly between local files (on the same machine) to remote files accessible over the internet.

If you're doing it right, you should never need to change the logic of your workflow to accommodate files coming from different locations.
All you need to do to use a remote file is to specify the appropriate prefix in the file path when you're supplying it to the workflow.

For example, `/path/to/data` has no prefix, indicating that it's a 'normal' local file path, whereas `s3://path/to/data` includes the `s3://` prefix, indicating that it's located in Amazon's S3 object storage.

Many different protocols are supported:

- HTTP(S)/FTP (http://, https://, ftp://)
- Amazon S3 (s3://)
- Azure Blob Storage (az://)
- Google Cloud Storage (gs://)

To use any of these, simply specify the relevant prefix in the string, which is then technically called a Uniform Resource Identifier (URI) instead of a file path.
Nextflow will handle authentication and staging the files to the right place, downloading or uploading and all other file operations you would expect.

The key strength of this system is that it enables us to switch between environments without changing any pipeline logic.
For example, you can develop with a small, local test set before switching to a full-scale test set located in remote storage simply by changing the URI.

### 2.1. Use a file from the internet

Let's test this out by switching the local path we're providing to our workflow with an HTTPS path pointing to a copy of the same data that is stored in Github.

!!! warning

    This will only work if you have an active internet connection.

Open `main.nf` again and change the input path as follows:

=== "After"

    ```groovy title="main.nf" linenums="2" hl_lines="2"
        // Using a remote file from the internet
        myFile = file('https://raw.github.com/nextflow-io/training/master/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz')

        // Print file attributes
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    ```

=== "Before"

    ```groovy title="main.nf" linenums="2" hl_lines="2"
        // Create a Path object from a string path
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Print file attributes
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    ```

Let's run the workflow:

```bash
nextflow run main.nf
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [insane_swartz] DSL2 - revision: fff18abe6d

    File object class: class nextflow.file.http.XPath
    File name: patientA_rep1_normal_R1_001.fastq.gz
    Simple name: patientA_rep1_normal_R1_001
    Extension: gz
    Parent directory: /nextflow-io/training/master/side-quests/working_with_files/data
    executor >  local (1)
    [8a/2ab7ca] COUNT_LINES [100%] 1 of 1 ✔
    Processing file: patientA_rep1_normal_R1_001.fastq.gz
    40
    ```

It works! You can see that very little has changed.

The one difference in the console output is that the path object class is now `nextflow.file.http.XPath`, whereas for the local path the class was `sun.nio.fs.UnixPath`.
You don't need to remember these classes; we just mention this to demonstrate that Nextflow identifies and handles the different locations appropriately.

Behind the scenes, Nextflow downloaded the file to a staging directory located within the work directory.
That staged file can then be treated as a local file and symlinked into the relevant process directory.

You can verify that that happened here by looking at the contents of the working directory located at the hash value of the process.

<!-- TODO (future) List work directory contents to show where the staging happens -->

Note that for larger files, the downloading step will take some extra time compared to running on local files.
However, Nextflow checks whether it already has a staged copy to avoid unnecessary downloads.
So if you run again on the same file and haven't deleted the staged file, Nextflow will use the staged copy.

This shows how easy it is to switch between local and remote data using Nextflow, which is a key feature of Nextflow.

!!! note

    The one important exception to this principle is that you can't use glob patterns or directory paths with HTTPS because HTTPS cannot list multiple files, so you must specify exact file URLs.
    However, other storage protocols such as blob storage (`s3://`, `az://`, `gs://`) can use both globs and directory paths.

    Here's how you could use glob patterns with cloud storage:

    ```groovy title="Cloud storage examples (not runnable in this environment)"
    // S3 with glob patterns - would match multiple files
    ch_s3_files = channel.fromPath('s3://my-bucket/data/*.fastq.gz')

    // Azure Blob Storage with glob patterns
    ch_azure_files = channel.fromPath('az://container/data/patient*_R{1,2}.fastq.gz')

    // Google Cloud Storage with glob patterns
    ch_gcs_files = channel.fromPath('gs://bucket/data/sample_*.fastq.gz')
    ```

    We'll show you how to work with globs in practice in the next section.

### 2.2. Switch back to the local file

We're going to go back to using our local example files for the rest of this side quest, so let's switch the workflow input back to the original file:

=== "After"

    ```groovy title="main.nf" linenums="2" hl_lines="2"
        // Create a Path object from a string path
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Print file attributes
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    ```

=== "Before"

    ```groovy title="main.nf" linenums="2" hl_lines="5-8"
        // Create a Path object from a string path
        myFile = file('https://raw.github.com/nextflow-io/training/master/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz')

        // Print file attributes
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    ```

### Takeaway

- Remote data is accessed using a URI (HTTP, FTP, S3, Azure, Google Cloud)
- Nextflow will automatically download and stage the data to the right place, as long as these paths are being fed to processes
- Do not write logic to download or upload remote files!
- Local and remote files produce different object types but work identically
- **Important**: HTTP/HTTPS only work with single files (no glob patterns)
- Cloud storage (S3, Azure, GCS) supports both single files and glob patterns
- You can seamlessly switch between local and remote data sources without changing code logic (as long as the protocol supports your required operations)

---

## 3. Using the `fromPath()` channel factory

So far we've been working with a single file at a time, but in Nextflow, we're typically going to want to create an input channel with multiple input files to process.

A naive way to do that would be to combine the `file()` method with [`channel.of()`](https://www.nextflow.io/docs/latest/reference/channel.html#of) like this:

```groovy title="Syntax example"
ch_files = channel.of([file('data/patientA_rep1_normal_R1_001.fastq.gz')],
                      [file('data/patientA_rep1_normal_R1_001.fastq.gz')])
```

That works, but it's clunky.

<!-- TODO (future) Discuss when it's good to use just file() vs channel.fromPath() -->

This is where [`channel.fromPath()`](https://www.nextflow.io/docs/latest/reference/channel.html#frompath) comes in: a convenient channel factory that bundles all the functionality we need to generate a channel from one or more static file strings as well as glob patterns.

### 3.1. Add the channel factory

Let's update our workflow to use `channel.fromPath`.

=== "After"

    ```groovy title="main.nf" linenums="7" hl_lines="1-3 6 12 15"
        // Load files with channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R1_001.fastq.gz')
        ch_files.view { myFile -> "Found file: $myFile" }

        // Print file attributes
        /* Comment these out for now, we'll come back to them!
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
        */

        // Count the lines in the file
        // COUNT_LINES(myFile)
    ```

=== "Before"

    ```groovy title="main.nf" linenums="7" hl_lines="1-2"
        // Create a Path object from a string path
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Print file attributes
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"

        // Count the lines in the file
        COUNT_LINES(myFile)
    ```

We've also commented out the code that prints out the attributes for now, and added a `.view` statement to print out just the filename instead.

Run the workflow:

```bash
nextflow run main.nf
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [grave_meucci] DSL2 - revision: b09964a583

    executor >  local (1)
    [5c/342c73] COUNT_LINES (1) [100%] 1 of 1 ✔
    Found file: /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz
    Processing file: patientA_rep1_normal_R1_001.fastq.gz
    40
    ```

As you can see, the file path is being loading as a `Path` type object in the channel.
This is similar to what `file()` would have done, except now we have a channel that we can load more files into if we want.

Using `channel.fromPath()` is a convenient way of creating a new channel populated by a list of files.

### 3.2. View attributes of files in channel

In our first pass at using the channel factory, we simplified the code and just printed out the file name.

Let's go back to printing out the full file attributes:

=== "After"

    ```groovy title="main.nf" linenums="7" hl_lines="3-9"
        // Loading files with channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R1_001.fastq.gz')
        ch_files.view { myFile ->
            println "File object class: ${myFile.class}"
            println "File name: ${myFile.name}"
            println "Simple name: ${myFile.simpleName}"
            println "Extension: ${myFile.extension}"
            println "Parent directory: ${myFile.parent}"
        }
    ```

=== "Before"

    ```groovy title="main.nf" linenums="7" hl_lines="3 6 12"
        // Load files with channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R1_001.fastq.gz')
        ch_files.view { myFile -> "Found file: $myFile" }

        // Print file attributes
        /* Comment these out for now, we'll come back to them!
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
        */
    ```

Since `myFile` is a proper Path object, we have access to all the same class attributes as before.

Run the workflow:

```bash
nextflow run main.nf
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [furious_swanson] DSL2 - revision: c35c34950d

    executor >  local (1)
    [9d/6701a6] COUNT_LINES (1) [100%] 1 of 1 ✔
    File object class: sun.nio.fs.UnixPath
    File name: patientA_rep1_normal_R1_001.fastq.gz
    Simple name: patientA_rep1_normal_R1_001
    Extension: gz
    Parent directory: /workspaces/training/side-quests/working_with_files/data
    Processing file: patientA_rep1_normal_R1_001.fastq.gz
    40
    ```

And there you are, same results as before but now we have the file in a channel, so we can add more.

### 3.3. Using a glob to match multiple files

There are several ways we could load more files into the channel.
Here we're going to show you how to use glob patterns, which are a convenient way to match and retrieve file and directory names based on wildcard characters.
The process of matching these patterns is called "globbing" or "filename expansion".

!!! note

    As noted previously, Nextflow supports globbing to manage input and output files in the majority of cases, except with HTTPS filepaths because HTTPS cannot list multiple files.

Let's say we want to retrieve both files in a pair of files associated with a given patient, `patientA`:

```console
patientA_rep1_normal_R1_001.fastq.gz
patientA_rep1_normal_R2_001.fastq.gz
```

Since the only difference between the filenames is the replicate number, _i.e._ the number after `R`, we can use the wildcard character `*` to stand in for the number as follows:

```console
patientA_rep1_normal_R*_001.fastq.gz
```

That is the glob pattern we need.

Now all we need to do is update the file path in the channel factory to use that glob pattern as follows:

=== "After"

    ```groovy title="main.nf" linenums="7"
      // Load files with channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
    ```

=== "Before"

    ```groovy title="main.nf" linenums="7"
      // Load files with channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R1_001.fastq.gz')
    ```

Nextflow will automatically recognize that this is a glob pattern and will handle it appropriately.

Run the workflow to test that out:

```bash
nextflow run main.nf
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [boring_sammet] DSL2 - revision: d2aa789c9a

    executor >  local (2)
    [3c/a65de5] COUNT_LINES (2) [100%] 2 of 2 ✔
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
    Processing file: patientA_rep1_normal_R1_001.fastq.gz
    40

    Processing file: patientA_rep1_normal_R2_001.fastq.gz
    40
    ```

As you can see, we now have two Path objects in our channel, which shows that Nextflow has done the filename expansion correctly, and has loaded and processed both files as expected.

Using this method, we can retrieve as many or as few files as we want just by changing the glob pattern. If we made it more generous, for example by replacing all the variable parts of the filenames by `*` (_e.g._ `data/patient*_rep*_*_R*_001.fastq.gz`) we could grab all the example files in the `data` directory.

### Takeaway

- `channel.fromPath()` creates a channel with files matching a pattern
- Each file is emitted as a separate element in the channel
- We can use a glob pattern to match multiple files
- Files are automatically converted to Path objects with full attributes
- The `.view()` method allows inspection of channel contents

---

## 4. Extracting basic metadata from filenames

In most scientific domains, it's very common to have metadata encoded in the names of the files that contain the data.
For example, in bioinformatics, files containing sequencing data are often named in a way that encodes information about the sample, condition, replicate, and read number.

If the filenames are constructed according to a consistent convention, you can extract that metadata in a standardized manner and use it in the course of your analysis.
That is a big 'if', of course, and you should be very cautious whenever you rely on filename structure; but the reality is that this approach is very widely used, so let's have a look at how it's done in Nextflow.

In the case of our example data, we know that the filenames include consistently structured metadata.
For example, the filename `patientA_rep1_normal_R2_001` encodes the following:

- patient ID: `patientA`
- replicate ID: `rep1`
- sample type: `normal` (as opposed to `tumor`)
- read set: `R1` (as opposed to `R2`)

We're going to modify our workflow to retrieve this information in three steps:

1. Retrieve the `simpleName` of the file, which includes the metadata
2. Separate the metadata using a method called `tokenize()`
3. Use a map to organize the metadata

!!! warning

    You should never encode sensitive information into filenames, such as patient names or other identifying characteristics, as that can compromise patient privacy or other relevant security restrictions.

### 4.1. Retrieve the `simpleName`

The `simpleName` is a file attribute that corresponds to the filename stripped of its path and extension.

Make the following edits to the workflow:

=== "After"

    ```groovy title="main.nf" linenums="7" hl_lines="3-6"
        // Load files with channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
        ch_files.map { myFile ->
            [ myFile.simpleName, myFile ]
        }
        .view()
    ```

=== "Before"

    ```groovy title="main.nf" linenums="7" hl_lines="3-9"
        // Load files with channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
        ch_files.view {
            println "File object class: ${myFile.class}"
            println "File name: ${it.name}"
            println "Simple name: ${it.simpleName}"
            println "Extension: ${it.extension}"
            println "Parent directory: ${it.parent}"
        }
    ```

This retrieves the `simpleName` and associates it with the full file object using a `map()` operation.

Run the workflow to test that it works:

```bash
nextflow run main.nf
```

??? success "Command output"

    ```console hl_lines="7-8"
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [suspicious_mahavira] DSL2 - revision: ae8edc4e48

    executor >  local (2)
    [e9/55774b] COUNT_LINES (2) [100%] 2 of 2 ✔
    [patientA_rep1_normal_R2_001, /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R2_001.fastq.gz]
    [patientA_rep1_normal_R1_001, /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz]
    Processing file: patientA_rep1_normal_R1_001.fastq.gz
    40

    Processing file: patientA_rep1_normal_R2_001.fastq.gz
    40
    ```

Each element in the channel is now a tuple containing the `simpleName` and the original file object.

### 4.2. Extract the metadata from the `simplename`

At this point, the metadata we want is embedded in the `simplename`, but we can't access individual items directly.
So we need to split the `simplename` into its components.
Fortunately, those components are simply separated by underscores in the original filename, so we can apply a common Nextflow method called `tokenize()` that is perfect for this task.

Make the following edits to the workflow:

=== "After"

    ```groovy title="main.nf" linenums="7" hl_lines="4"
        // Load files with channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
        ch_files.map { myFile ->
            [ myFile.simpleName.tokenize('_'), myFile ]
        }
    ```

=== "Before"

    ```groovy title="main.nf" linenums="7" hl_lines="4"
        // Load files with channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
        ch_files.map { myFile ->
            [ myFile.simpleName, myFile ]
        }
    ```

The `tokenize()` method will split the `simpleName` string wherever it finds underscores, and will return a list containing the substrings.

Run the workflow:

```bash
nextflow run main.nf
```

??? success "Command output"

    ```console hl_lines="7-8"
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [gigantic_gauss] DSL2 - revision: a39baabb57

    executor >  local (2)
    [e7/da2f4b] COUNT_LINES (2) [100%] 2 of 2 ✔
    [[patientA, rep1, normal, R2, 001], /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R2_001.fastq.gz]
    [[patientA, rep1, normal, R1, 001], /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz]
    Processing file: patientA_rep1_normal_R2_001.fastq.gz
    40

    Processing file: patientA_rep1_normal_R1_001.fastq.gz
    40
    ```

Now the tuple for each element in our channel contains the list of metadata (_e.g._ `[patientA, rep1, normal, R1, 001]`) and the original file object.

That's great!
We've broken down our patient information from a single string into a list of strings.
We can now handle each part of the patient information separately.

### 4.3. Use a map to organize the metadata

Our metadata is just a flat list at the moment.
It's easy enough to use but difficult to read.

```console
[patientA, rep1, normal, R1, 001]
```

What is the item at index 3? Can you tell without referring back to the original explanation of the metadata structure?

This is a great opportunity to use a key-value store, where every item has a set of keys and their associated values, so you can easily refer to each key to get the corresponding value.

In our example, that means going from this organization:

```groovy
data = [patientA, 1, normal, R1]

println data[3]
```

To this one:

```groovy
data = [id: patientA, replicate: 1, type: normal, readNum: 1]

println data.readNum
```

In Nextflow, that's called a [map](https://nextflow.io/docs/latest/script.html#maps).

Let's convert our flat list into a map now.
Make the following edits to the workflow:

=== "After"

    ```groovy title="main.nf" linenums="7" hl_lines="4-13"
        // Load files with channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
        ch_files.map { myFile ->
            def (patient, replicate, type, readNum) = myFile.simpleName.tokenize('_')
            [
              [
                id: patient,
                replicate: replicate.replace('rep', ''),
                type: type,
                readNum: readNum.replace('R', ''),
              ],
              myFile
            ]
        }
    ```

=== "Before"

    ```groovy title="main.nf" linenums="7" hl_lines="4"
        // Load files with channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
        ch_files.map { myFile ->
            [ myFile.simpleName.tokenize('_'), myFile ]
        }
    ```

<!-- TODO (future) Explain the map a little more? -->

While we're at it, we also simplified a couple of the metadata strings using a string replacement method called `replace()` to remove some characters that are unnecessary (_e.g._ `readNum.replace('rep', '')` to keep only the number from the replicate IDs).

Let's run the workflow again:

```bash
nextflow run main.nf
```

??? success "Command output"

    ```console hl_lines="7-8"
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [infallible_swartz] DSL2 - revision: 7f4e68c0cb

    executor >  local (2)
    [1b/e7fb27] COUNT_LINES (1) [100%] 2 of 2 ✔
    [[id:patientA, replicate:1, type:normal, readNum:2], /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R2_001.fastq.gz]
    [[id:patientA, replicate:1, type:normal, readNum:1], /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz]
    Processing file: patientA_rep1_normal_R2_001.fastq.gz
    40

    Processing file: patientA_rep1_normal_R1_001.fastq.gz
    40
    ```

Now the metadata is neatly labeled (_e.g._ `[id:patientA, replicate:1, type:normal, readNum:2]`) so it's a lot easier to tell what is what.

It'll also be a lot easier to actually make use of elements of metadata in the workflow, and will make our code easier to read and more maintainable.

### Takeaway

- We can handle filenames in Nextflow with the power of a full programming language
- We can treat the filenames as strings to extract relevant information
- Use of methods like `tokenize()` and `replace()` allows us to manipulate strings in the filename
- The `.map()` operation transforms channel elements while preserving structure
- Structured metadata (maps) makes code more readable and maintainable than positional lists

Next up, we will look at how to handle paired data files.

---

## 5. Handling paired data files

Many experimental designs produce paired data files that benefit from being handled in an explicitly paired way.
For example, in bioinformatics, sequencing data is often generated in the form of paired reads, meaning sequence strings that originate from the same fragment of DNA (often called 'forward' and 'reverse' because they are read from opposite ends).

That is the case of our example data, where R1 and R2 refer to the two sets of reads.

```console
data/patientA_rep1_normal_R1_001.fastq.gz
data/patientA_rep1_normal_R2_001.fastq.gz
```

Nextflow provides a specialized channel factory for working with paired files like this called `channel.fromFilePairs()`, which automatically groups files based on a shared naming pattern. That allows you to associate the paired files more tightly with less effort.

We're going to modify our workflow to take advantage of this.
It's going to take two steps:

1. Switch the channel factory to `channel.fromFilePairs()`
2. Extract and map the metadata

### 5.1. Switch the channel factory to `channel.fromFilePairs()`

To use `channel.fromFilePairs`, we need to specify the pattern that Nextflow should use to identify the two members in a pair.

Going back to our example data, we can formalize the naming pattern as follows:

```console
data/patientA_rep1_normal_R{1,2}_001.fastq.gz
```

This is similar to the glob pattern we used earlier, except this specifically enumerates the substrings (either `1` or `2` coming right after the R) that identify the two members of the pair.

Let's update the workflow `main.nf` accordingly:

=== "After"

    ```groovy title="main.nf" linenums="7" hl_lines="1-2"
        // Load files with channel.fromFilePairs
        ch_files = channel.fromFilePairs('data/patientA_rep1_normal_R{1,2}_001.fastq.gz')
        /* Comment out the mapping for now, we'll come back to it!
        ch_files.map { myFile ->
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
        */
        .view()
    ```

=== "Before"

    ```groovy title="main.nf" linenums="7" hl_lines="1-2"
        // Load files with channel.fromFilePairs
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
        ch_files.map { myFile ->
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

We've switched the channel factory and adapted the file matching pattern, and while we were at it, we commented out the map operation.
We'll add that back in later, with a few modifications.

Run the workflow to test it:

```bash
nextflow run main.nf
```

??? failure "Command output"

    ```console hl_lines="7-8"
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [angry_koch] DSL2 - revision: 44fdf66105

    [-        ] COUNT_LINES -
    [-        ] COUNT_LINES -
    [patientA_rep1_normal_R, [/workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz, /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R2_001.fastq.gz]]
    ERROR ~ Error executing process > 'COUNT_LINES (1)'

    Caused by:
      Not a valid path value: 'patientA_rep1_normal_R'



    Tip: when you have fixed the problem you can continue the execution adding the option `-resume` to the run command line

    -- Check '.nextflow.log' file for details
    ```

Uh-oh, this time the run failed!

The relevant bit of the error message is here:

```console
Not a valid path value: 'patientA_rep1_normal_R'
```

That is because we've changed the channel factory.
Until now, the original input channel only contained the file paths.
All the metadata manipulation we've been doing didn't actually affect the channel contents.

Now that we're using the `.fromFilePairs` channel factory, the contents of the resulting channel are different.
We see only one channel element, composed of a tuple containing two items: the part of the `simpleName` shared by the two files, which serves as an identifier, and a tuple containing the two file objects, in the format `id, [ file1, file2 ]`.

That's great, because Nextflow has done the hard work of extracting the patient name by examining the shared prefix and using it as a patient identifier.

However, it does break our current workflow.
If we wanted to still run `COUNT_LINES` the same way without changing the process, we would have to apply a mapping operation to extract the file paths.
But we're not going to do that, because our ultimate goal is to use a different process, `ANALYZE_READS`, that handles file pairs appropriately.

So let's simply comment out (or delete) the call to `COUNT_LINES` and move on.

=== "After"

    ```groovy title="main.nf" linenums="26" hl_lines="2"
        // Count the lines in the file
        // COUNT_LINES(ch_files)
    ```

=== "Before"

    ```groovy title="main.nf" linenums="26" hl_lines="2"
        // Count the lines in the file
        COUNT_LINES(ch_files)
    ```

You can also comment out or delete the `COUNT_LINES` include statement, but that will have no functional effect.

Now let's run the workflow again:

```bash
nextflow run main.nf
```

??? success "Command output"

    ```console hl_lines="7-8"
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [fabulous_davinci] DSL2 - revision: 22b53268dc

    [patientA_rep1_normal_R, [/workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz, /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R2_001.fastq.gz]]
    ```

Yay, this time the workflow succeeds!

However, we still need to get the rest of the metadata out of the `id` field.

### 5.2. Extract and organize metadata from file pairs

Our `map` operation from before won't work because it doesn't match the data structure, but we can modify it to work.

We already have access to the actual patient identifier in the string that `fromFilePairs()` used as an identifier, so we can use that to extract the metadata without getting the `simpleName` from the Path object like we did before.

Uncomment the map operation in the workflow and make the following edits:

=== "After"

    ```groovy title="main.nf" linenums="7" hl_lines="3-4 9 11 13"
        // Load files with channel.fromFilePairs
        ch_files = channel.fromFilePairs('data/patientA_rep1_normal_R{1,2}_001.fastq.gz')
        ch_files.map { id, files ->
            def (sample, replicate, type) = id.tokenize('_')
            [
                [
                    id: sample,
                    replicate: replicate.replace('rep', ''),
                    type: type
                ],
                files
            ]
        }
        .view()
    ```

=== "Before"

    ```groovy title="main.nf" linenums="7" hl_lines="3-5 11 13"
        // Load files with channel.fromFilePairs
        ch_files = channel.fromFilePairs('data/patientA_rep1_normal_R{1,2}_001.fastq.gz')
        /* Comment out the mapping for now, we'll come back to it!
        ch_files.map { myFile ->
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
        */
        .view()
    ```

This time the map starts from `id, files` instead of just `myFile`, and `tokenize()` is applied to `id` instead of to `myFile.simpleName`.

Notice also that we've dropped `readNum` from the `tokenize()` line; any substrings that we don't specifically name (starting from the left) will be silently dropped.
We can do this because the paired files are now tightly associated, so we no longer need `readNum` in the metadata map.

Let's run the workflow:

```bash
nextflow run main.nf
```

??? success "Command output"

    ```console

    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [prickly_stonebraker] DSL2 - revision: f62ab10a3f

    [[id:patientA, replicate:1, type:normal], [/workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz, /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R2_001.fastq.gz]]
    ```

And there it is: we have the metadata map (`[id:patientA, replicate:1, type:normal]`) in the first position of the output tuple, followed by the tuple of paired files, as intended.

Of course, this will only pick up and process that specific pair of files.
If you want to experiment with processing multiple pairs, you can try adding wildcards into the input pattern and see what happens.
Foe example, try using `data/patientA_rep1_*_R{1,2}_001.fastq.gz`

### Takeaway

- [`channel.fromFilePairs()` automatically finds and pairs related files](https://www.nextflow.io/docs/latest/reference/channel.html#fromfilepairs)
- This simplifies handling paired-end reads in your pipeline
- Paired files can be grouped as `[id, [file1, file2]]` tuples
- Metadata extraction can be done from the paired file ID rather than individual files

---

## 6. Using file operations in processes

Now let's put all this together in a simple process to reinforce how to use file operations inside a Nextflow process.

We provide you with a pre-written process module called `ANALYZE_READS` that takes a tuple of metadata and a pair of input files and analyzes them.
We could imagine this is doing sequence alignment, or variant calling or any other step that makes sense for this data type.

Let's get started.

### 6.1. Import the process and examine the code

To use this process in the workflow, we just need to add a module include statement before the workflow block.

Make the following edit to the workflow:

=== "After"

    ```groovy title="main.nf" linenums="1" hl_lines="3"
    #!/usr/bin/env nextflow

    include { ANALYZE_READS } from './modules/analyze_reads.nf'

    workflow {
    ```

=== "Before"

    ```groovy title="main.nf" linenums="1"
    #!/usr/bin/env nextflow

    workflow {
    ```

You can open the module file to examine its code:

```groovy title="modules/analyze_reads.nf - process example" linenums="1"
#!/usr/bin/env nextflow

process ANALYZE_READS {
    tag "${meta.id}"

    publishDir "results/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(files)

    output:
    tuple val(meta.id), path("${meta.id}_stats.txt")

    script:
    """
    echo "Sample metadata: ${meta.id}" > ${meta.id}_stats.txt
    echo "Replicate: ${meta.replicate}" >> ${meta.id}_stats.txt
    echo "Type: ${meta.type}" >> ${meta.id}_stats.txt
    echo "Read 1: ${files[0]}" >> ${meta.id}_stats.txt
    echo "Read 2: ${files[1]}" >> ${meta.id}_stats.txt
    echo "Read 1 size: \$(gunzip -dc ${files[0]} | wc -l | awk '{print \$1/4}') reads" >> ${meta.id}_stats.txt
    echo "Read 2 size: \$(gunzip -dc ${files[1]} | wc -l | awk '{print \$1/4}') reads" >> ${meta.id}_stats.txt
    """
}
```

!!! note

    We are calling our metadata map `meta` by convention.
    For a deeper dive into meta maps, see the [Metadata and meta maps](./metadata.md) side quest.

### 6.2. Call the process in the workflow

Now that the process is available to the workflow, we can add a call to the `ANALYZE_READS` process to run it.

To run it on our example data, we'll need to do two things:

1. Give a name to the remapped channel
2. Add a call to the process

#### 6.2.1. Name the remapped input channel

We previously applied the mapping manipulations directly to the input channel.
In order to feed the remapped contents to the `ANALYZE_READS` process (and do so in a way that is clear and easy to read) we want to create a new channel named `ch_samples`.

We can do that using the [`set`](https://www.nextflow.io/docs/latest/reference/operator.html#set) operator.

In the main workflow, replace the `.view()` operator with `.set { ch_samples }`, and add a line testing that we can refer to the channel by name.

=== "After"

    ```groovy title="main.nf" linenums="7" hl_lines="21 23-24"
        // Load files with channel.fromFilePairs
        ch_files = channel.fromFilePairs('data/patientA_rep1_normal_R{1,2}_001.fastq.gz')
        ch_files.map { id,  files ->
           def (sample, replicate, type, readNum) = id.tokenize('_')
           [
               [
                   id: sample,
                   replicate: replicate.replace('rep', ''),
                   type: type
               ],
               files
           ]
        }
            .set { ch_samples }

        // Temporary: peek into ch_samples
        ch_samples.view()
    ```

=== "Before"

    ```groovy title="main.nf" linenums="7" hl_lines="21"
        // Load files with channel.fromFilePairs
        ch_files = channel.fromFilePairs('data/patientA_rep1_normal_R{1,2}_001.fastq.gz')
        ch_files.map { id,  files ->
           def (sample, replicate, type, readNum) = id.tokenize('_')
           [
               [
                   id: sample,
                   replicate: replicate.replace('rep', ''),
                   type: type
               ],
               files
           ]
        }
        .view()
    }
    ```

Let's run this:

```bash
nextflow run main.nf
```

??? success "Command output"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [goofy_kirch] DSL2 - revision: 3313283e42

    [[id:patientA, replicate:1, type:normal], [/workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz, /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R2_001.fastq.gz]]
    ```

This confirms we can now refer to the channel by name.

#### 6.2.2. Call the process on the data

Now let's actually call the `ANALYZE_READS` process on the `ch_samples` channel.

In the main workflow, make the following code changes:

=== "After"

    ```groovy title="main.nf" linenums="23"
        // Temporary: peek into ch_samples
        ch_samples.view()
    ```

=== "Before"

    ```groovy title="main.nf" linenums="23"
        // Run the analysis
        ANALYZE_READS(ch_samples)
    ```

Let's run this:

```bash
nextflow run main.nf
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `./main.nf` [shrivelled_cori] DSL2 - revision: b546a31769

    executor >  local (1)
    [b5/110360] process > ANALYZE_READS (patientA) [100%] 1 of 1 ✔
    ```

This process is set up to publish its outputs to a `results` directory, so have a look in there.

??? abstract "Directory and file contents"

    ```console
    results
    └── patientA
        └── patientA_stats.txt
    ```

    ```txt title="patientA_stats.txt"
    Sample metadata: patientA
    Replicate: 1
    Type: normal
    Read 1: patientA_rep1_normal_R1_001.fastq.gz
    Read 2: patientA_rep1_normal_R2_001.fastq.gz
    Read 1 size: 10 reads
    Read 2 size: 10 reads
    ```

The process took our inputs and created a new file containing the patient metadata, as designed.
Splendid!

### 6.3. Include many more patients

Of course, this is just processing a single pair of files for a single patient, which is not exactly the kind of high throughput you're hoping to get with Nextflow.
You'll probably want to process a lot more data at a time.

Remember `channel.fromPath()` accepts a _glob_ as input, which means it can accept any number of files that match the pattern.
Therefore if we want to include all the patients, we can simply modify the input string to include more patients, as noted in passing earlier.

Let's pretend we want to be as greedy as possible.
Make the following edits to the workflow:

=== "After"

    ```groovy title="main.nf" linenums="7" hl_lines="2"
        // Load files with channel.fromFilePairs
        ch_files = channel.fromFilePairs('data/*_R{1,2}_001.fastq.gz')
    ```

=== "Before"

    ```groovy title="main.nf" linenums="7" hl_lines="2"
        // Load files with channel.fromFilePairs
        ch_files = channel.fromFilePairs('data/patientA_rep1_normal_R{1,2}_001.fastq.gz')
    ```

Run the pipeline again:

```bash
nextflow run main.nf
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `./main.nf` [big_stonebraker] DSL2 - revision: f7f9b8a76c

    executor >  local (8)
    [d5/441891] process > ANALYZE_READS (patientC) [100%] 8 of 8 ✔
    ```

The results directory should now contain results for all the available data.

??? abstract "Directory contents"

    ```console
    results
    ├── patientA
    │   └── patientA_stats.txt
    ├── patientB
    │   └── patientB_stats.txt
    └── patientC
        └── patientC_stats.txt
    ```

Success! We have analyzed all the patients in one go! Right?

Maybe not.
If you look more closely, we have a problem: we have two replicates for patientA, but only one output file!
We are overwriting the output file each time.

### 6.4. Make the published files unique

Since we have access to the patient metadata, we can use it to make the published files unique by including differentiating metadata, either in the directory structure or in the filenames themselves.

Make the following change to the workflow:

=== "After"

    ```groovy title="modules/analyze_reads.nf" linenums="6"
        publishDir "results/${meta.type}/${meta.id}/${meta.replicate}", mode: 'copy'
    ```

=== "Before"

    ```groovy title="modules/analyze_reads.nf" linenums="6"
        publishDir "results/${id}", mode: 'copy'
    ```

Here we show the option of using additional directory levels to account for sample types and replicates, but you could experiment with doing it at the filename level as well.

Now run the pipeline one more time, but be sure to remove the results directory first to give yourself a clean workspace:

```bash
rm -r results
nextflow run main.nf
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `./main.nf` [insane_swartz] DSL2 - revision: fff18abe6d

    executor >  local (8)
    [e3/449081] process > ANALYZE_READS (patientC) [100%] 8 of 8 ✔
    ```

Check the results directory now:

??? abstract "Directory contents"

    ```console
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

And there it is, all our metadata, neatly organized. That's success!

There's a lot more you can do once you have your metadata loaded into a map like this:

1. Create organized output directories based on patient attributes
2. Make decisions in processes based on patient properties
3. Split, join, and recombine data based on metadata values

This pattern of keeping metadata explicit and attached to the data (rather than encoded in filenames) is a core best practice in Nextflow that enables building robust, maintainable analysis workflows.
You can learn more about this in the [Metadata and meta maps](./metadata.md) side quest.

### Takeaway

- The `publishDir` directive can organize outputs based on metadata values
- Metadata in tuples enables structured organization of results
- This approach creates maintainable workflows with clear data provenance
- Processes can take tuples of metadata and files as input
- The `tag` directive provides process identification in execution logs
- Workflow structure separates channel creation from process execution

---

## Summary

In this side quest, you've learned how to work with files in Nextflow, from basic operations to more advanced techniques for handling collections of files.

Applying these techniques in your own work will enable you to build more efficient and maintainable workflows, especially when working with large numbers of files with complex naming conventions.

### Key patterns

1.  **Basic File Operations:** We created Path objects with `file()` and accessed file attributes like name, extension, and parent directory, learning the difference between strings and Path objects.

    - Create a Path object with `file()`

    ```groovy
    myFile = file('path/to/file.txt')
    ```

    - Get file attributes

    ```groovy
    println myFile.name       // file.txt
    println myFile.baseName   // file
    println myFile.extension  // txt
    println myFile.parent     // path/to
    ```

2.  **Using Remote Files**: We learned how to transparently switch between local and remote files using URIs, demonstrating Nextflow's ability to handle files from various sources without changing workflow logic.

    - Local file

    ```groovy
    myFile = file('path/to/file.txt')
    ```

    - FTP

    ```groovy
    myFile = file('ftp://path/to/file.txt')
    ```

    - HTTPS

    ```groovy
    myFile = file('https://path/to/file.txt')
    ```

    - Amazon S3

    ```groovy
    myFile = file('s3://path/to/file.txt')
    ```

    - Azure Blob Storage

    ```groovy
    myFile = file('az://path/to/file.txt')
    ```

    - Google Cloud Storage

    ```groovy
    myFile = file('gs://path/to/file.txt')
    ```

3.  **Loading files using the `fromPath()` channel factory:** We created channels from file patterns with `channel.fromPath()` and viewed their file attributes, including object types.

    - Create a channel from a file pattern

    ```groovy
     ch_files = channel.fromPath('data/*.fastq.gz')
    ```

    - Get file attributes

    ```groovy
     ch_files.view { myFile ->
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    }
    ```

4.  **Extracting Patient Metadata from Filenames:** We used `tokenize()` and `replace()` to extract and structure metadata from filenames, converting them to organized maps.

    ```groovy
    def name = file.name.tokenize('_')
    def patientId = name[0]
    def replicate = name[1].replace('rep', '')
    def type = name[2]
    def readNum = name[3].replace('R', '')
    ```

5.  **Simplifying with channel.fromFilePairs:** We used `channel.fromFilePairs()` to automatically pair related files and extract metadata from paired file IDs.

    ```groovy
    ch_pairs = channel.fromFilePairs('data/*_R{1,2}_001.fastq.gz')
    ```

6.  **Using File Operations in Processes:** We integrated file operations into Nextflow processes with proper input handling, using `publishDir` to organize outputs based on metadata.

    - Associate a meta map with the process inputs

    ```groovy
    ch_files = channel.fromFilePairs('data/patientA_rep1_normal_R{1,2}_001.fastq.gz')
    ch_files.map { id,  files ->
        def (sample, replicate, type, readNum) = id.tokenize('_')
        [
            [
                id: sample,
                replicate: replicate.replace('rep', ''),
                type: type
            ],
             files
        ]
    }
        .set { ch_samples }

    ANALYZE_READS(ch_samples)
    ```

    - Organize outputs based on metadata

    ```groovy
    publishDir "results/${meta.type}/${meta.id}/${meta.replicate}", mode: 'copy'
    ```

### Additional resources

- [Nextflow Documentation: Working with Files](https://www.nextflow.io/docs/latest/working-with-files.html)
- [channel.fromPath](https://www.nextflow.io/docs/latest/reference/channel.html#frompath)
- [channel.fromFilePairs](https://www.nextflow.io/docs/latest/reference/channel.html#fromfilepairs)

---

## What's next?

Return to the [menu of Side Quests](./index.md) or click the button in the bottom right of the page to move on to the next topic in the list.
