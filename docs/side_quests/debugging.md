# Debugging Workflows

Debugging is a critical skill that can save you hours of frustration and help you become a more effective Nextflow developer. Throughout your career, especially when you're starting out, you'll encounter bugs while building and maintaining your workflows. Learning systematic debugging approaches will help you identify and resolve issues quickly.

### Learning goals

In this side quest, we'll explore **systematic debugging techniques** for Nextflow workflows:

- **Syntax error debugging**: Using IDE features and Nextflow error messages effectively
- **Channel debugging**: Diagnosing data flow issues and channel structure problems
- **Process debugging**: Investigating execution failures and resource issues
- **Built-in debugging tools**: Leveraging Nextflow's preview mode, stub running, and work directories
- **Systematic approaches**: A four-phase methodology for efficient debugging

By the end, you'll have a robust debugging methodology that transforms frustrating error messages into clear roadmaps for solutions.

### Prerequisites

Before taking on this side quest, you should:

- Have completed the [Hello Nextflow](../hello_nextflow/README.md) tutorial or equivalent beginner's course.
- Be comfortable using basic Nextflow concepts and mechanisms (processes, channels, operators)

**Optional:** We recommend completing the [IDE Features for Nextflow Development](./ide_features.md) side quest first.
That covers comprehensive coverage of IDE features that support debugging (syntax highlighting, error detection, etc.), which we'll use heavily here.

---

## 0. Get started

#### Open the training codespace

If you haven't yet done so, make sure to open the training environment as described in the [Environment Setup](../envsetup/index.md).

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Move into the project directory

Let's move into the directory where the files for this tutorial are located.

```bash
cd side-quests/debugging
```

You can set VSCode to focus on this directory:

```bash
code .
```

#### Review the materials

You'll find a set of example workflows with various types of bugs that we'll use for practice:

```console title="Directory contents"
.
├── bad_bash_var.nf
├── bad_channel_shape.nf
├── bad_channel_shape_viewed_debug.nf
├── bad_channel_shape_viewed.nf
├── bad_number_inputs.nf
├── badpractice_syntax.nf
├── bad_resources.nf
├── bad_syntax.nf
├── buggy_workflow.nf
├── data
│   ├── sample_001.fastq.gz
│   ├── sample_002.fastq.gz
│   ├── sample_003.fastq.gz
│   ├── sample_004.fastq.gz
│   ├── sample_005.fastq.gz
│   └── sample_data.csv
├── exhausted.nf
├── invalid_process.nf
├── missing_output.nf
├── missing_software.nf
├── missing_software_with_stub.nf
├── nextflow.config
└── no_such_var.nf
```

These files represent common debugging scenarios you'll encounter in real-world development.

#### Review the assignment

Your challenge is to run each workflow, identify the error(s) and fix it.

<!-- TODO: give a bit more details, similar to how it's done in the Metadata side quest -->

#### Readiness checklist

Think you're ready to dive in?

- [ ] I understand the goal of this course and its prerequisites
- [ ] My codespace is up and running
- [ ] I've set my working directory appropriately
- [ ] I understand the assignment

If you can check all the boxes, you're good to go.

---

## 1. Syntax Errors

Syntax errors are the most common type of error you'll encounter when writing Nextflow code. They occur when the code does not conform to the expected syntax rules of the Nextflow DSL. These errors prevent your workflow from running at all, so it's important to learn how to identify and fix them quickly.

### 1.1. Missing braces

One of the most common syntax errors, and sometimes one of the more complex ones to debug is **missing or mismatched brackets**.

Let's start with a practical example.

#### Run the pipeline

```bash
nextflow run bad_syntax.nf
```

You'll see an error message like this:

```console title="Syntax error output"
 N E X T F L O W   ~  version 25.04.3

Launching `bad_syntax.nf` [stupefied_bhabha] DSL2 - revision: ca6327fad2

ERROR ~ Script compilation error
- file : /workspaces/training/side-quests/debugging/bad_syntax.nf
- cause: Unexpected input: '{' @ line 3, column 23.
   process PROCESS_FILES {
                         ^

1 error

NOTE: If this is the beginning of a process or workflow, there may be a syntax error in the body, such as a missing or extra comma, for which a more specific error message could not be produced.

 -- Check '.nextflow.log' file for details
```

**Key elements of syntax error messages:**

- **File location**: Shows exactly which file contains the error (`- file : /workspaces/training/side-quests/debugging/bad_syntax.nf`)
- **Error description**: Explains what the parser found that it didn't expect (`- cause: Unexpected input: '{'`)
- **Line and column**: Points to where the parser encountered the problem (`@ line 3, column 23.`)
- **Context**: Shows the problematic line with a caret (^) pointing to location of an unclosed brace (`process PROCESS_FILES {`)
- **Additional notes**: Provides hints about common causes

#### Check the code

Now, let's examine `bad_syntax.nf` to understand what's causing the error:

```groovy title="bad_syntax.nf" hl_lines="14" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    input:
    val sample_name

    output:
    path "${sample_name}_output.txt"

    script:
    """
    echo "Processing ${sample_name}" > ${sample_name}_output.txt
    """
// Missing closing brace for the process

workflow {

    // Create input channel
    input_ch = channel.of('sample1', 'sample2', 'sample3')

    // Call the process with the input channel
    PROCESS_FILES(input_ch)
}
```

For the purpose of this example we've left a comment for you to show where the error is. The Nextflow VSCode extension should also be giving you some hints about what might be wrong, putting the mismatched brace in red and highlighting the premature end of the file:

![Bad syntax](img/bad_syntax.png)

**Debugging strategy for bracket errors:**

1. Use VS Code's bracket matching (place cursor next to a bracket)
2. Check the Problems panel for bracket-related messages
3. Ensure each opening `{` has a corresponding closing `}`

#### Fix the code

Replace the comment with the missing closing brace:

=== "After"

    ```groovy title="bad_syntax.nf" hl_lines="14" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        """
        echo "Processing ${sample_name}" > ${sample_name}_output.txt
        """
    }  // Add the missing closing brace

    workflow {

        // Create input channel
        input_ch = channel.of('sample1', 'sample2', 'sample3')

        // Call the process with the input channel
        PROCESS_FILES(input_ch)
    }
    ```

=== "Before"

    ```groovy title="bad_syntax.nf" hl_lines="14" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        """
        echo "Processing ${sample_name}" > ${sample_name}_output.txt
        """
    // Missing closing brace for the process

    workflow {

        // Create input channel
        input_ch = channel.of('sample1', 'sample2', 'sample3')

        // Call the process with the input channel
        PROCESS_FILES(input_ch)
    }
    ```

#### Run the pipeline

Now run the workflow again to confirm it works:

```bash
nextflow run bad_syntax.nf
```

### 1.2. Using incorrect process keywords or directives

Another common syntax error is an **invalid process definition**. This can happen if you forget to define required blocks or use incorrect directives in the process definition.

#### Run the pipeline

```bash
nextflow run invalid_process.nf
```

You'll see an error like:

```console title="Invalid process keyword error"
 N E X T F L O W   ~  version 25.04.3

Launching `invalid_process.nf` [nasty_jepsen] DSL2 - revision: da9758d614

ERROR ~ Script compilation error
- file : /workspaces/training/side-quests/debugging/invalid_process.nf
- cause: Invalid process definition -- Unknown keyword `inputs` @ line 5, column 5.
       val sample_name
       ^

1 error


 -- Check '.nextflow.log' file for details
```

#### Check the code

Let's examine `invalid_process.nf` to see what's wrong:

```groovy title="invalid_process.nf" hl_lines="4" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    inputs:  // ERROR: Should be 'input' not 'inputs'
    val sample_name

    output:
    path "${sample_name}_output.txt"

    script:
    """
    echo "Processing ${sample_name}" > ${sample_name}_output.txt
    """
}

workflow {

    // Create input channel
    input_ch = channel.of('sample1', 'sample2', 'sample3')

    // Call the process with the input channel
    PROCESS_FILES(input_ch)
}
```

The error message was quite straightforward - we're using `inputs` instead of the correct `input` directive. You'll also see that the Nextflow VSCode exension is unhappy:

![Invalid process message](img/invalid_process_message.png)

#### Fix the code

Replace the incorrect keyword with the correct one by referencing [the documentation](https://www.nextflow.io/docs/latest/process.html#):

=== "After"

    ```groovy title="invalid_process.nf" hl_lines="4" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:  // Fixed: Changed 'inputs' to 'input'
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        """
        echo "Processing ${sample_name}" > ${sample_name}_output.txt
        """
    }

    workflow {

        // Create input channel
        input_ch = channel.of('sample1', 'sample2', 'sample3')

        // Call the process with the input channel
        PROCESS_FILES(input_ch)
    }
    ```

=== "Before"

    ```groovy title="invalid_process.nf" hl_lines="4" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        inputs:  // ERROR: Should be 'input' not 'inputs'
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        """
        echo "Processing ${sample_name}" > ${sample_name}_output.txt
        """
    }

    workflow {

        // Create input channel
        input_ch = channel.of('sample1', 'sample2', 'sample3')

        // Call the process with the input channel
        PROCESS_FILES(input_ch)
    }
    ```

#### Run the pipeline

Now run the workflow again to confirm it works:

```bash
nextflow run invalid_process.nf
```

### 1.3. Using bad variable names

The variable names you use in your script blocks must be valid, derived either from inputs or from groovy code inserted before the script. But when you're wrangling complexity at the start of pipeline development, it's easy to make mistakes in variable naming, and Nextflow will let you know quickly.

#### Run the pipeline

```bash
nextflow run no_such_var.nf
```

You should get a failure that looks like this:

```console title="No such variable error"
ERROR ~ Error executing process > 'PROCESS_FILES (3)'

Caused by:
  No such variable: undefined_var -- Check script 'no_such_var.nf' at line: 15


Source block:
  def output_prefix = "${sample_name}_processed"
  def timestamp = new Date().format("yyyy-MM-dd")
  """
  echo "Processing ${sample_name} on ${timestamp}" > ${output_prefix}.txt
  echo "Using undefined variable: ${undefined_var}" >> ${output_prefix}.txt  // ERROR: undefined_var not defined
  """

Tip: when you have fixed the problem you can continue the execution adding the option `-resume` to the run command line

 -- Check '.nextflow.log' file for details
```

#### Check the code

Let's examine `no_such_var.nf`:

```groovy title="no_such_var.nf" hl_lines="17" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    input:
    val sample_name

    output:
    path "${sample_name}_output.txt"

    script:
    // Define variables in Groovy code before the script
    def output_prefix = "${sample_name}_processed"
    def timestamp = new Date().format("yyyy-MM-dd")

    """
    echo "Processing ${sample_name} on ${timestamp}" > ${output_prefix}.txt
    echo "Using undefined variable: ${undefined_var}" >> ${output_prefix}.txt  // ERROR: undefined_var not defined
    """
}

workflow {
    input_ch = channel.of('sample1', 'sample2', 'sample3')
    PROCESS_FILES(input_ch)
}
```

The error message indicates that the variable is not recognized in the script template, and there you go- you should be able to see `${undefined_var}` used in the script block, but not defined elsewhere.

#### Fix the code

If you get a 'No such variable' error, you can fix it by either defining the variable (by correcting input variable names or editing groovy code before the script), or by removing it from the script block if it's not needed:

=== "After"

    ```groovy title="no_such_var.nf" hl_lines="15-17" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        // Define variables in Groovy code before the script
        def output_prefix = "${sample_name}_processed"
        def timestamp = new Date().format("yyyy-MM-dd")

        """
        echo "Processing ${sample_name} on ${timestamp}" > ${output_prefix}.txt
        """  // Removed the line with undefined_var
    }

    workflow {
        input_ch = channel.of('sample1', 'sample2', 'sample3')
        PROCESS_FILES(input_ch)
    }
    ```

=== "Before"

    ```groovy title="no_such_var.nf" hl_lines="17" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        // Define variables in Groovy code before the script
        def output_prefix = "${sample_name}_processed"
        def timestamp = new Date().format("yyyy-MM-dd")

        """
        echo "Processing ${sample_name} on ${timestamp}" > ${output_prefix}.txt
        echo "Using undefined variable: ${undefined_var}" >> ${output_prefix}.txt  // ERROR: undefined_var not defined
        """
    }

    workflow {
        input_ch = channel.of('sample1', 'sample2', 'sample3')
        PROCESS_FILES(input_ch)
    }
    ```

#### Run the pipeline

Now run the workflow again to confirm it works:

```bash
nextflow run no_such_var.nf
```

### 1.4. Bad use of Bash variables

Starting out in Nextflow, it can be difficult to understand the difference between Nextflow (Groovy) and Bash variables. This can generate another form of the bad variable error that appears when trying to use variables in the Bash content of the script block.

#### Run the pipeline

```bash
nextflow run bad_bash_var.nf
```

This throws the following error:

```console
ERROR ~ Error executing process > 'PROCESS_FILES (1)'

Caused by:
  No such variable: prefix -- Check script 'bad_bash_var.nf' at line: 11

```

#### Check the code

Let's examine `bad_bash_var.nf` to see what's causing the issue:

```groovy title="bad_bash_var.nf" hl_lines="13" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    input:
    val sample_name

    output:
    path "${sample_name}_output.txt"

    script:
    """
    prefix="${sample_name}_output"
    echo "Processing ${sample_name}" > ${prefix}.txt  # ERROR: ${prefix} is Groovy syntax, not Bash
    """
}
```

In this example, we're defining the `prefix` variable in Bash, but in a Nexflow process the `$` syntax we used to refer to it (`${prefix}`) is interpretes as a Groovy variable, not Bash. The variable doesn't exist in the Groovy context, so we get a 'no such variable' error.

#### Fix the code

If you want to use a Bash variable, you must escape the dollar sign like this:

=== "After"

    ```groovy title="bad_bash_var.nf" hl_lines="13" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        """
        prefix="${sample_name}_output"
        echo "Processing ${sample_name}" > \${prefix}.txt  # Fixed: Escaped the dollar sign
        """
    }

    workflow {
        input_ch = channel.of('sample1', 'sample2', 'sample3')
        PROCESS_FILES(input_ch)
    }
    ```

=== "Before"

    ```groovy title="bad_bash_var.nf" hl_lines="13" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        """
        prefix="${sample_name}_output"
        echo "Processing ${sample_name}" > ${prefix}.txt  # ERROR: ${prefix} is Groovy syntax, not Bash
        """
    }
    ```

This tells Nextflow to interpret this as a Bash variable.

#### Run the pipeline

Now run the workflow again to confirm it works:

```bash
nextflow run bad_bash_var.nf
```

!!! tip "Groovy vs Bash Variables"

    For simple variable manipulations like string concatenation or prefix/suffix operations, it's usually more readable to use Groovy variables in the script section rather than Bash variables in the script block:

    ```groovy linenums="1"
    script:
    def output_prefix = "${sample_name}_processed"
    def output_file = "${output_prefix}.txt"
    """
    echo "Processing ${sample_name}" > ${output_file}
    """
    ```

    This approach avoids the need to escape dollar signs and makes the code easier to read and maintain.

### 1.5. Non-Fatal Syntax Warnings

The Nextflow VSCode extension sometimes highlights issues that won't prevent execution but represent bad practices. For example, it's currently possible to define channels outside of the `workflow {}` block, but it's not good practice and the extension will highlight this as a potential issue. These warnings help you write better code even though they're not fatal errors.

#### Run the pipeline

```bash
nextflow run badpractice_syntax.nf
```

When you run this workflow, it will execute successfully:

```console title="Successful execution despite bad practice"
N E X T F L O W   ~  version 25.04.3

Launching `badpractice_syntax.nf` [peaceful_euler] DSL2 - revision: 7b2c9a1d45

executor >  local (3)
[a1/b2c3d4] process > PROCESS_FILES (1) [100%] 3 of 3 ✔
```

#### Check the code

Let's examine `badpractice_syntax.nf` to see what the VSCode extension is warning about:

```groovy title="badpractice_syntax.nf" hl_lines="3" linenums="1"
#!/usr/bin/env nextflow

input_ch = channel.of('sample1', 'sample2', 'sample3')  // WARNING: Channel defined outside workflow

process PROCESS_FILES {
    input:
    val sample_name

    output:
    path "${sample_name}_processed.txt"

    script:
    // Define variables in Groovy code before the script
    def output_prefix = "${sample_name}_processed"
    def timestamp = new Date().format("yyyy-MM-dd")

    """
    echo "Processing ${sample_name} on ${timestamp}" > ${output_prefix}.txt
    """
}

workflow {
    PROCESS_FILES(input_ch)
}
```

The VSCode extension will highlight the `input_ch` variable as being defined outside the workflow block, which is not recommended:

![Non-lethal syntax error](img/nonlethal.png)

This won't prevent execution but could lead to confusion or unexpected behavior in larger workflows.

#### Fix the code

Follow the VSCode extension's recommendation by moving the channel definition inside the workflow block:

=== "After"

    ```groovy title="badpractice_syntax.nf" hl_lines="21" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_processed.txt"

        script:
        // Define variables in Groovy code before the script
        def output_prefix = "${sample_name}_processed"
        def timestamp = new Date().format("yyyy-MM-dd")

        """
        echo "Processing ${sample_name} on ${timestamp}" > ${output_prefix}.txt
        """
    }

    workflow {
        input_ch = channel.of('sample1', 'sample2', 'sample3')  // Moved inside workflow block
        PROCESS_FILES(input_ch)
    }
    ```

=== "Before"

    ```groovy title="badpractice_syntax.nf" hl_lines="3" linenums="1"
    #!/usr/bin/env nextflow

    input_ch = channel.of('sample1', 'sample2', 'sample3')  // WARNING: Channel defined outside workflow

    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_processed.txt"

        script:
        // Define variables in Groovy code before the script
        def output_prefix = "${sample_name}_processed"
        def timestamp = new Date().format("yyyy-MM-dd")

        """
        echo "Processing ${sample_name} on ${timestamp}" > ${output_prefix}.txt
        """
    }

    workflow {
        PROCESS_FILES(input_ch)
    }
    ```

#### Run the pipeline

Run the workflow again to confirm it still works and the VSCode warning is resolved:

```bash
nextflow run badpractice_syntax.nf
```

Tighter restrictions on such things will likely become enforced in future Nextflow versions, so it's good practice to keep your input channels defined within the workflow block, and in general to follow any other recommendations the extension makes.

### Takeaway

You can systematically identify and fix syntax errors using Nextflow error messages and IDE visual indicators. Common syntax errors include missing braces, incorrect process keywords, undefined variables, and improper use of Bash vs. Nextflow variables. The VSCode extension helps catch many of these before runtime. With these syntax debugging skills in your toolkit, you'll be able to quickly resolve the most common Nextflow syntax errors and move on to tackling more complex runtime issues.

### What's next?

Learn to debug more complex channel structure errors that occur even when syntax is correct.

---

## 2. Channel Structure Errors

Channel structure errors are more subtle than syntax errors because the code is syntactically correct, but the data shapes don't match what processes expect. Nextflow will try to run the pipeline, but might find that the number of inputs doesn't match what it expects and fail. These errors typically only appear at runtime and require an understanding of the data flowing through your workflow.

!!! tip "Debugging Channels with `.view()`"

    Throughout this section, remember that you can use the `.view()` operator to inspect channel content at any point in your workflow. This is one of the most powerful debugging tools for understanding channel structure issues. We'll explore this technique in detail in section 2.4, but feel free to use it as you work through the examples.

    ```groovy
    my_channel.view()  // Shows what's flowing through the channel
    ```

### 2.1. Wrong Number of Input Channels

This error occurs when you pass a different number of channels than a process expects.

#### Run the pipeline

```bash
nextflow run bad_number_inputs.nf
```

```console title="Wrong number of channels error"
 N E X T F L O W   ~  version 25.04.3

Launching `bad_number_inputs.nf` [high_mendel] DSL2 - revision: 955705c51b

Process `PROCESS_FILES` declares 1 input channel but 2 were specified

 -- Check script 'bad_number_inputs.nf' at line: 23 or see '.nextflow.log' file for more details
```

#### Check the code

The error message clearly states that the process expects 1 input channel, but 2 were provided. Let's examine `bad_number_inputs.nf`:

```groovy title="bad_number_inputs.nf" hl_lines="5 23" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    input:
        val sample_name  // Process expects only 1 input

    output:
        path "${sample_name}_output.txt"

    script:
    """
    echo "Processing ${sample_name}" > ${sample_name}_output.txt
    """
}

workflow {

    // Create two separate channels
    samples_ch = channel.of('sample1', 'sample2', 'sample3')
    files_ch = channel.of('file1.txt', 'file2.txt', 'file3.txt')

    // ERROR: Passing 2 channels but process expects only 1
    PROCESS_FILES(samples_ch, files_ch)
}
```

You should see the mismatched `PROCESS_FILES` call, supplying multiple input channels when the process only defines one. The VSCode extension will also under line process call in red, and supply a diagnostic message when you mouse over:

![Incorrect number of args message](img/incorrect_num_args.png)

#### Fix the code

For this specific example, the process expects a single channel and doesn't require the second channel, so we can fix it by passing only the `samples_ch` channel:

=== "After"

    ```groovy title="bad_number_inputs.nf" hl_lines="23" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
            val sample_name  // Process expects only 1 input

        output:
            path "${sample_name}_output.txt"

        script:
        """
        echo "Processing ${sample_name}" > ${sample_name}_output.txt
        """
    }

    workflow {

        // Create two separate channels
        samples_ch = channel.of('sample1', 'sample2', 'sample3')
        files_ch = channel.of('file1.txt', 'file2.txt', 'file3.txt')

        // Fixed: Pass only the channel the process expects
        PROCESS_FILES(samples_ch)
    }
    ```

=== "Before"

    ```groovy title="bad_number_inputs.nf" hl_lines="5 23" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
            val sample_name  // Process expects only 1 input

        output:
            path "${sample_name}_output.txt"

        script:
        """
        echo "Processing ${sample_name}" > ${sample_name}_output.txt
        """
    }

    workflow {

        // Create two separate channels
        samples_ch = channel.of('sample1', 'sample2', 'sample3')
        files_ch = channel.of('file1.txt', 'file2.txt', 'file3.txt')

        // ERROR: Passing 2 channels but process expects only 1
        PROCESS_FILES(samples_ch, files_ch)
    }
    ```

#### Run the pipeline

```bash
nextflow run bad_number_inputs.nf
```

More commonly than this example, you might add additional inputs to a process and forget to update the workflow call accordingly, which can lead to this type of error. Fortunately, this is one of the easier-to-understand and fix errors, as the error message is quite clear about the mismatch.

### 2.2. Channel Exhaustion (Process Runs Fewer Times Than Expected)

Some channel structure errors are much more subtle and produce no errors at all. Probably the most common of these reflects a challenge that new Nextflow users face in understanding that queue channels can be exhausted and run out of items, meaning the workflow finishes prematurely.

#### Run the pipeline

```bash
nextflow run exhausted.nf
```

When you run this workflow, it will execute without error, processing a single sample:

```console title="Exhausted channel output"
 N E X T F L O W   ~  version 25.04.3

Launching `exhausted.nf` [extravagant_gauss] DSL2 - revision: 08cff7ba2a

executor >  local (1)
[bd/f61fff] PROCESS_FILES (1) [100%] 1 of 1 ✔
```

#### Check the code

Let's examine `exhausted.nf` to see if that's right:

```groovy title="exhausted.nf" hl_lines="23 24" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    input:
    val reference
    val sample_name

    output:
    path "${output_prefix}.txt"

    script:
    // Define variables in Groovy code before the script
    output_prefix = "${reference}_${sample_name}"
    def timestamp = new Date().format("yyyy-MM-dd")

    """
    echo "Processing ${sample_name} on ${timestamp}" > ${output_prefix}.txt
    """
}

workflow {

    reference_ch = channel.of('baseline_reference')
    input_ch = channel.of('sample1', 'sample2', 'sample3')

    PROCESS_FILES(reference_ch, input_ch)
}
```

The process only runs once instead of three times because the `reference_ch` channel is a queue channel that gets exhausted after the first process execution. When one channel is exhausted, the entire process stops, even if other channels still have items.

This is a common pattern where you have a single reference file that needs to be reused across multiple samples. The solution is to convert the reference channel to a value channel that can be reused indefinitely.

#### Fix the code

There are a couple of ways to address this depending on how many files are affected.

**Option 1**: You have a single reference file that you are re-using a lot. You can simply create a value channel type, which can be used over and over again. There are three ways to do this:

**1a** Use `channel.value()`:

```groovy title="exhausted.nf (fixed - Option 1a)" hl_lines="2" linenums="21"
workflow {
    reference_ch = channel.value('baseline_reference')  // Value channel can be reused
    input_ch = channel.of('sample1', 'sample2', 'sample3')

    PROCESS_FILES(reference_ch, input_ch)
}
```

**1b** Use the `first()` [operator](https://www.nextflow.io/docs/latest/reference/operator.html#first):

```groovy title="exhausted.nf (fixed - Option 1b)" hl_lines="2" linenums="21"
workflow {
    reference_ch = channel.of('baseline_reference').first()  // Convert to value channel
    input_ch = channel.of('sample1', 'sample2', 'sample3')

    PROCESS_FILES(reference_ch, input_ch)
}
```

**1c.** Use the `collect()` [operator](https://www.nextflow.io/docs/latest/reference/operator.html#collect):

```groovy title="exhausted.nf (fixed - Option 1c)" hl_lines="2" linenums="21"
workflow {
    reference_ch = channel.of('baseline_reference').collect()  // Convert to value channel
    input_ch = channel.of('sample1', 'sample2', 'sample3')

    PROCESS_FILES(reference_ch, input_ch)
}
```

**Option 2**: In more complex scenarios, perhaps where you have multiple reference files for all samples in the sample channel, you can use the `combine` operator to create a new channel that combines the two channels into tuples:

```groovy title="exhausted.nf (fixed - Option 2)" hl_lines="4" linenums="21"
workflow {
    reference_ch = channel.of('baseline_reference','other_reference')
    input_ch = channel.of('sample1', 'sample2', 'sample3')
    combined_ch = reference_ch.combine(input_ch)  // Creates cartesian product

    PROCESS_FILES(combined_ch)
}
```

The `.combine()` operator generates a cartesian product of the two channels, so each item in `reference_ch` will be paired with each item in `input_ch`. This allows the process to run for each sample while still using the reference.

This requires the process input to be adjusted. In our example, the start of the process definition would need to be adjusted as follows:

```groovy title="exhausted.nf (fixed - Option 2)" hl_lines="5" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    input:
        tuple val(reference), val(sample_name)
```

This approach may not be suitable in all situations.

#### Run the pipeline

Try one of the fixes above and run the workflow again:

```bash
nextflow run exhausted.nf
```

You should now see all three samples being processed instead of just one.

### 2.3. Wrong Channel Content Structure

When workflows reach a certain level of complexity, it can be a little difficult to keep track of the internal structures of each channel, and people commonly generate mismatches between what the process expects and what the channel actually contains. This is more subtle than the issue we discussed earlier, where the number of channels was incorrect. In this case, you can have the correct number of input channels, but the internal structure of one or more of those channels doesn't match what the process expects.

#### Run the pipeline

```bash
nextflow run bad_channel_shape.nf
```

You will see an error like this:

```console title="Channel structure error"
Launching `bad_channel_shape.nf` [hopeful_pare] DSL2 - revision: ffd66071a1

executor >  local (3)
executor >  local (3)
[3f/c2dcb3] PROCESS_FILES (3) [  0%] 0 of 3 ✘
ERROR ~ Error executing process > 'PROCESS_FILES (1)'

Caused by:
  Missing output file(s) `[sample1, file1.txt]_output.txt` expected by process `PROCESS_FILES (1)`


Command executed:

  echo "Processing [sample1, file1.txt]" > [sample1, file1.txt]_output.txt

Command exit status:
  0

Command output:
  (empty)

Work dir:
  /workspaces/training/side-quests/debugging/work/d6/1fb69d1d93300bbc9d42f1875b981e

Tip: when you have fixed the problem you can continue the execution adding the option `-resume` to the run command line

 -- Check '.nextflow.log' file for details
```

#### Check the code

The square brackets in the error message provide the clue here - the process is treating the tuple as a single value, which is not what we want. Let's examine `bad_channel_shape.nf`:

```groovy title="bad_channel_shape.nf" hl_lines="5 20-22" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    input:
        val sample_name  // Expects single value, gets tuple

    output:
        path "${sample_name}_output.txt"

    script:
    """
    echo "Processing ${sample_name}" > ${sample_name}_output.txt
    """
}

workflow {

    // Channel emits tuples, but process expects single values
    input_ch = channel.of(
      ['sample1', 'file1.txt'],
      ['sample2', 'file2.txt'],
      ['sample3', 'file3.txt']
    )
    PROCESS_FILES(input_ch)
}
```

You can see that we're generating a channel composed of tuples: `['sample1', 'file1.txt']`, but the process expects a single value, `val sample_name`. The command executed shows that the process is trying to create a file named `[sample3, file3.txt]_output.txt`, which is not the intended output.

#### Fix the code

To fix this, if the process requires both inputs we could adjust the process to accept a tuple:

=== "Option 1: Accept tuple in process"

    === "After"

        ```groovy title="bad_channel_shape.nf" hl_lines="5"  linenums="1"
        #!/usr/bin/env nextflow

        process PROCESS_FILES {
            input:
                tuple val(sample_name), path(file_name)  // Fixed: Accept tuple

            output:
                path "${sample_name}_output.txt"

            script:
            """
            echo "Processing ${sample_name}" > ${sample_name}_output.txt
            """
        }

        workflow {

            // Channel emits tuples, but process expects single values
            input_ch = channel.of(
              ['sample1', 'file1.txt'],
              ['sample2', 'file2.txt'],
              ['sample3', 'file3.txt']
            )
            PROCESS_FILES(input_ch)
        }
        ```

    === "Before"

        ```groovy title="bad_channel_shape.nf" hl_lines="5" linenums="1"
        #!/usr/bin/env nextflow

        process PROCESS_FILES {
            input:
                val sample_name  // Expects single value, gets tuple

            output:
                path "${sample_name}_output.txt"

            script:
            """
            echo "Processing ${sample_name}" > ${sample_name}_output.txt
            """
        }

        workflow {

            // Channel emits tuples, but process expects single values
            input_ch = channel.of(
              ['sample1', 'file1.txt'],
              ['sample2', 'file2.txt'],
              ['sample3', 'file3.txt']
            )
            PROCESS_FILES(input_ch)
        }
        ```

=== "Option 2: Extract first element"

    === "After"

        ```groovy title="bad_channel_shape.nf" hl_lines="9" linenums="16"
        workflow {

            // Channel emits tuples, but process expects single values
            input_ch = channel.of(
              ['sample1', 'file1.txt'],
              ['sample2', 'file2.txt'],
              ['sample3', 'file3.txt']
            )
            PROCESS_FILES(input_ch.map { it[0] })  // Fixed: Extract first element
        }
        ```

    === "Before"

        ```groovy title="bad_channel_shape.nf" hl_lines="9" linenums="16"
        workflow {

            // Channel emits tuples, but process expects single values
            input_ch = channel.of(
              ['sample1', 'file1.txt'],
              ['sample2', 'file2.txt'],
              ['sample3', 'file3.txt']
            )
            PROCESS_FILES(input_ch)
        }
        ```

#### Run the pipeline

Pick one of the solutions and re-run the workflow:

```bash
nextflow run bad_channel_shape.nf
```

### 2.4. Channel Debugging Techniques

#### Using `.view()` for Channel Inspection

The most powerful debugging tool for channels is the `.view()` operator. With `.view()`, you can understand the shape of your channels at all stages to help with debugging.

#### Run the pipeline

Run `bad_channel_shape_viewed.nf` to see this in action:

```bash
nextflow run bad_channel_shape_viewed.nf
```

You'll see output like this:

```console title="Channel debugging output"
 N E X T F L O W   ~  version 25.04.3

Launching `bad_channel_shape_viewed.nf` [maniac_poisson] DSL2 - revision: b4f24dc9da

executor >  local (3)
[c0/db76b3] PROCESS_FILES (3) [100%] 3 of 3 ✔
Channel content: [sample1, file1.txt]
Channel content: [sample2, file2.txt]
Channel content: [sample3, file3.txt]
After mapping: sample1
After mapping: sample2
After mapping: sample3
```

#### Check the code

Let's examine `bad_channel_shape_viewed.nf` to see how `.view()` is used:

```groovy title="bad_channel_shape_viewed.nf" linenums="16" hl_lines="9 11"
workflow {

    // Channel emits tuples, but process expects single values
    input_ch = channel.of(
      ['sample1', 'file1.txt'],
      ['sample2', 'file2.txt'],
      ['sample3', 'file3.txt']
    )
    .view { "Channel content: $it" }  // Debug: Show original channel content
    .map { tuple -> tuple[0] }        // Transform: Extract first element
    .view { "After mapping: $it" }    // Debug: Show transformed channel content

    PROCESS_FILES(input_ch)
}
```

#### Fix the code

To save you from using `.view()` operations excessively in future to understand channel content, it's advisable to add some comments to help:

```groovy title="bad_channel_shape_viewed.nf (with comments)" linenums="16" hl_lines="8 9"
workflow {

    // Channel emits tuples, but process expects single values
    input_ch = channel.of(
            ['sample1', 'file1.txt'],
            ['sample2', 'file2.txt'],
            ['sample3', 'file3.txt'],
        ) // [sample_name, file_name]
        .map { tuple -> tuple[0] } // sample_name

    PROCESS_FILES(input_ch)
}
```

This will become more important as your workflows grow in complexity and channel structure becomes more opaque.

#### Run the pipeline

```bash
nextflow run bad_channel_shape_viewed.nf
```

### Takeaway

Many channel structure errors can be created with valid Nextflow syntax. You can debug channel structure errors by understanding data flow, using `.view()` operators for inspection, and recognizing error message patterns like square brackets indicating unexpected tuple structures.

### What's next?

Learn about errors created by process definitions.

---

## 3. Process Structure Errors

Most of the errors you encounter related to processes will related to mistakes you have made in forming the command, or to issues related to the underlying software. That said, similarly to the channel issues above, you can make mistakes in the process definition that don't quality as syntax errors, but which will cause errors at run time.

### 3.1. Missing Output Files

One common error when writing processes is to do something that generates a mismatch between what the process expects and what is generated.

#### Run the pipeline

```bash
nextflow run missing_output.nf
```

You'll see an error like this:

```console title="Missing output files error"
 N E X T F L O W   ~  version 25.04.3

Launching `missing_output.nf` [zen_stone] DSL2 - revision: 37ff61f926

executor >  local (3)
executor >  local (3)
[fd/2642e9] process > PROCESS_FILES (2) [ 66%] 2 of 3, failed: 2
ERROR ~ Error executing process > 'PROCESS_FILES (3)'

Caused by:
  Missing output file(s) `sample3.txt` expected by process `PROCESS_FILES (3)`


Command executed:

  echo "Processing sample3" > sample3_output.txt

Command exit status:
  0

Command output:
  (empty)

Work dir:
  /workspaces/training/side-quests/debugging/work/02/9604d49fb8200a74d737c72a6c98ed

Tip: when you have fixed the problem you can continue the execution adding the option `-resume` to the run command line

 -- Check '.nextflow.log' file for details
```

#### Check the code

The error message indicates that the process expected to produce an output file named `sample3.txt`, but the script actually creates `sample3_output.txt`. Let's examine the process definition in `missing_output.nf`:

```groovy title="missing_output.nf" linenums="3" hl_lines="6 10"
process PROCESS_FILES {
    input:
    val sample_name

    output:
    path "${sample_name}.txt"  // Expects: sample3.txt

    script:
    """
    echo "Processing ${sample_name}" > ${sample_name}_output.txt  // Creates: sample3_output.txt
    """
}
```

You should see that there is a mismatch between the output file name in the `output:` block, and the one used in the script. This mismatch causes the process to fail. If you encounter this sort of error, go back and check that the outputs match between your process definition and your output block.

If the problem still isn't clear, check the work directory itself to identify the actual output files created:

```bash
❯ ls -h work/02/9604d49fb8200a74d737c72a6c98ed
sample3_output.txt
```

For this example this would highlight to us that a `_output` suffix is being incorporated into the output file name, contrary to our `output:` definition.

#### Fix the code

Fix the mismatch by making the output filename consistent:

=== "After"

    ```groovy title="missing_output.nf" hl_lines="6 10" linenums="3"
    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_output.txt"  // Fixed: Match the script output

        script:
        """
        echo "Processing ${sample_name}" > ${sample_name}_output.txt
        """
    }
    ```

=== "Before"

    ```groovy title="missing_output.nf" hl_lines="6 10" linenums="3"
    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}.txt"  // Expects: sample3.txt

        script:
        """
        echo "Processing ${sample_name}" > ${sample_name}_output.txt  // Creates: sample3_output.txt
        """
    }
    ```

#### Run the pipeline

```bash
nextflow run missing_output.nf
```

### 3.2. Missing software

Another class of errors occurs due to mistakes in software provisioning. `missing_software.nf` is a syntactically valid workflow, but it depends on some external software to provide the `cowpy` command it uses.

#### Run the pipeline

```bash
nextflow run missing_software.nf
```

You will see an error like this:

```console title="Missing software error" hl_lines="12 18"
ERROR ~ Error executing process > 'PROCESS_FILES (3)'

Caused by:
  Process `PROCESS_FILES (3)` terminated with an error exit status (127)


Command executed:

  cowpy sample3 > sample3_output.txt

Command exit status:
  127

Command output:
  (empty)

Command error:
  .command.sh: line 2: cowpy: command not found

Work dir:
  /workspaces/training/side-quests/debugging/work/82/42a5bfb60c9c6ee63ebdbc2d51aa6e

Tip: you can try to figure out what's wrong by changing to the process work directory and showing the script file named `.command.sh`

 -- Check '.nextflow.log' file for details
```

The process doesn't have access to the command we're specifying. Sometimes this is because a script is present in the workflow `bin` directory, but has not been made executable. Other times it is because the software is not installed in the container or environment where the workflow is running.

#### Check the code

Look out for that `127` exit code - it tells you exactly the problem. Let's examine `missing_software.nf`:

```groovy title="missing_software.nf" linenums="3" hl_lines="3"
process PROCESS_FILES {

    container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

    input:
    val sample_name

    output:
    path "${sample_name}_output.txt"

    script:
    """
    cowpy ${sample_name} > ${sample_name}_output.txt
    """
}
```

#### Fix the code

We've been a little disingenuous here, and there's actually nothing wrong with the code. We just need to specify the necessary configuration to run the process in such a way that it has access to the command in question. In this case the process has a container definition, so all we need to do is run the workflow with Docker enabled.

#### Run the pipeline

We've set up a Docker profile for you in `nextflow.config`, so you can run the workflow with:

```bash
nextflow run missing_software.nf -profile docker
```

This should run successfully now.

!!! note

    To learn more about how nextflow uses containers, go back to [Hello Nextflow](../hello_nextflow/05_hello_containers.md)

### 3.3. Bad resource configuration

In production usage, you'll be configuring resources on your processes. For example `memory` defines the maximum amount of memory available to your process, and if the process exceeds that, your scheduler will typically kill the process and return an exit code of `137`. We can't demonstrate that here because we're using the `local` executor, but we can show something similar with `time`.

#### Run the pipeline

`bad_resources.nf` has process configuration with an unrealistic bound on time of 1 millisecond:

```bash
nextflow run bad_resources.nf -profile docker
```

This gives us an error:

```console title="Resource time limit error"
ERROR ~ Error executing process > 'PROCESS_FILES (1)'

Caused by:
Process exceeded running time limit (1ms)
```

#### Check the code

Let's examine `bad_resources.nf`:

```groovy title="bad_resources.nf" linenums="3" hl_lines="3"
process PROCESS_FILES {

    time '1 ms'  // ERROR: Unrealistic time limit

    input:
    val sample_name

    output:
    path "${sample_name}_output.txt"

    script:
    """
    sleep 1  // Takes 1 second, but time limit is 1ms
    cowpy ${sample_name} > ${sample_name}_output.txt
    """
}
```

We know the process will take longer than a second (we've added a sleep in there to make sure), but the process is set to time out after 1 millisecond. Someone has been a little unrealistic with their configuration!

#### Fix the code

Increase the time limit to a realistic value:

=== "After"

    ```groovy title="bad_resources.nf" hl_lines="3" linenums="3"
    process PROCESS_FILES {

        time '100 s'  // Fixed: Realistic time limit

        input:
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        """
        sleep 1
        cowpy ${sample_name} > ${sample_name}_output.txt
        """
    }
    ```

=== "Before"

    ```groovy title="bad_resources.nf" hl_lines="3" linenums="3"
    process PROCESS_FILES {

        time '1 ms'  // ERROR: Unrealistic time limit

        input:
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        """
        sleep 1  // Takes 1 second, but time limit is 1ms
        cowpy ${sample_name} > ${sample_name}_output.txt
        """
    }
    ```

#### Run the pipeline

```bash
nextflow run bad_resources.nf -profile docker
```

If you make sure to read your error messages failures like this should not puzzle you for too long. But make sure you understand the resource requirements of the commands you are running so that you can configure your resource directives appropriately.

### 3.4. Process Debugging Techniques

When processes fail or behave unexpectedly, you need systematic techniques to investigate what went wrong. The work directory contains all the information you need to debug process execution.

#### Using Work Directory Inspection

The most powerful debugging tool for processes is examining the work directory. When a process fails, Nextflow creates a work directory for that specific process execution containing all the files needed to understand what happened.

#### Run the pipeline

Let's use the `missing_output.nf` example from earlier to demonstrate work directory inspection (re-generate an output naming mismatch if you need to):

```bash
nextflow run missing_output.nf
```

You'll see an error like this:

```console title="Missing output error"
  Missing output file(s) `sample2.txt` expected by process `PROCESS_FILES (2)`
```

#### Check the work directory

When you get this error, the work directory contains all the debugging information. Find the work directory path from the error message and examine its contents:

```bash
# Find the work directory from the error message
ls work/02/9604d49fb8200a74d737c72a6c98ed/
```

You can then examine the key files:

##### Check the Command Script

The `.command.sh` file shows exactly what command was executed:

```bash
# View the executed command
cat work/02/9604d49fb8200a74d737c72a6c98ed/.command.sh
```

This reveals:

- **Variable substitution**: Whether Nextflow variables were properly expanded
- **File paths**: Whether input files were correctly located
- **Command structure**: Whether the script syntax is correct

Common issues to look for:

- **Missing quotes**: Variables containing spaces need proper quoting
- **Wrong file paths**: Input files that don't exist or are in wrong locations
- **Incorrect variable names**: Typos in variable references
- **Missing environment setup**: Commands that depend on specific environments

##### Check Error Output

The `.command.err` file contains the actual error messages:

```bash
# View error output
cat work/02/9604d49fb8200a74d737c72a6c98ed/.command.err
```

This file will show:

- **Exit codes**: 127 (command not found), 137 (killed), etc.
- **Permission errors**: File access issues
- **Software errors**: Application-specific error messages
- **Resource errors**: Memory/time limit exceeded

##### Check Standard Output

The `.command.out` file shows what your command produced:

```bash
# View standard output
cat work/02/9604d49fb8200a74d737c72a6c98ed/.command.out
```

This helps verify:

- **Expected output**: Whether the command produced the right results
- **Partial execution**: Whether the command started but failed partway through
- **Debug information**: Any diagnostic output from your script

##### Check the Exit Code

The `.exitcode` file contains the exit code for the process:

```bash
# View exit code
cat work/*/*/.exitcode
```

Common exit codes and their meanings:

- **Exit code 127**: Command not found - check software installation
- **Exit code 137**: Process killed - check memory/time limits

##### Check File Existence

When processes fail due to missing output files, check what files were actually created:

```bash
# List all files in the work directory
ls -la work/02/9604d49fb8200a74d737c72a6c98ed/
```

This helps identify:

- **File naming mismatches**: Output files with different names than expected
- **Permission issues**: Files that couldn't be created
- **Path problems**: Files created in wrong directories

In our example earlier, this confirmed to us that while our expected `sample3.txt` wasn't present, `sample3_output.txt` was:

```bash
❯ ls -h work/02/9604d49fb8200a74d737c72a6c98ed
sample3_output.txt
```

### Takeaway

Process debugging requires examining work directories to understand what went wrong. Key files include `.command.sh` (the executed script), `.command.err` (error messages), and `.command.out` (standard output). Exit codes like 127 (command not found) and 137 (process killed) provide immediate diagnostic clues about the type of failure.

### What's next?

Learn about Nextflow's built-in debugging tools and systematic approaches to troubleshooting.

---

## 4. Built-in Debugging Tools and Advanced Techniques

Nextflow provides several powerful built-in tools for debugging and analyzing workflow execution. These tools help you understand what went wrong, where it went wrong, and how to fix it efficiently.

### 4.1. Real-time Process Output

Sometimes you need to see what's happening inside running processes. You can enable real-time process output, which shows you exactly what each task is doing as it executes.

#### Run the pipeline

`bad_channel_shape_viewed.nf` from our earlier examples printed channel content using `.view()`, but we can also use the `debug` directive to echo variables from within the process itself, which we demonstrate in `bad_channel_shape_viewed_debug.nf`. Run the workflow:

```bash
nextflow run bad_channel_shape_viewed_debug.nf
```

You will see output like this:

```console title="Real-time process output"
 N E X T F L O W   ~  version 25.04.3

Launching `bad_channel_shape_viewed_debug.nf` [agitated_crick] DSL2 - revision: ea3676d9ec

executor >  local (3)
[c6/2dac51] process > PROCESS_FILES (3) [100%] 3 of 3 ✔
Channel content: [sample1, file1.txt]
Channel content: [sample2, file2.txt]
Channel content: [sample3, file3.txt]
After mapping: sample1
After mapping: sample2
After mapping: sample3
Sample name inside process is sample2

Sample name inside process is sample1

Sample name inside process is sample3
```

#### Check the code

Let's examine `bad_channel_shape_viewed_debug.nf` to see how the `debug` directive works:

```groovy title="bad_channel_shape_viewed_debug.nf" linenums="3" hl_lines="2"
process PROCESS_FILES {
    debug true  // Enable real-time output

    input:
    val sample_name

    output:
    path "${sample_name}_output.txt"

    script:
    """
    echo "Sample name inside process is ${sample_name}"
    echo "Processing ${sample_name}" > ${sample_name}_output.txt
    """
}
```

The `debug` directive can be a quick and convenient way to understand the environment of a process.

### 4.2. Preview Mode

Sometimes you want to catch problems before any processes run. Nextflow provides a flag for this kind of proactive debugging: `-preview`.

#### Run the pipeline

The preview mode lets you test workflow logic without executing commands. This can be quite useful for quickly checking the structure of your workflow and ensuring that processes are connected correctly without running any actual commands.

For example, for our first syntax error from earlier:

```bash
nextflow run bad_syntax.nf -preview
```

!!! note

    If you fixed the file, reintroduce the syntax error by changing `input` to `inputs` before you run the command

You'll see output like this:

```console title="Preview mode output"
 N E X T F L O W   ~  version 25.04.3

Launching `bad_syntax.nf` [sick_fermi] DSL2 - revision: ca6327fad2

ERROR ~ Script compilation error
- file : /workspaces/training/side-quests/debugging/bad_syntax.nf
- cause: Unexpected input: '{' @ line 3, column 23.
   process PROCESS_FILES {
                         ^

1 error

NOTE: If this is the beginning of a process or workflow, there may be a syntax error in the body, such as a missing or extra comma, for which a more specific error message could not be produced.

 -- Check '.nextflow.log' file for details
```

Preview mode is particularly useful for catching syntax errors early without running any processes. It validates the workflow structure and process connections before execution.

### 4.3. Stub Running for Logic Testing

Sometimes errors are difficult to debug because commands take too long, require special software, or fail for complex reasons. Stub running lets you test workflow logic without executing the actual commands.

#### Run the pipeline

When you're developing a Nextflow process, you can use the `stub` directive to define 'dummy' commands that generate outputs of the correct form without running the real command. This approach is particularly valuable when you want to verify that your workflow logic is correct before dealing with the complexities of the actual software.

For example, remember our `missing_software.nf` from earlier? The one where we had missing software that prevented the workflow running until we added `-profile docker`? `missing_software_with_stub.nf` is a very similar workflow. If we run it in the same way, we will generate the same error:

```bash
nextflow run missing_software_with_stub.nf
```

```console title="Missing software error with stub" hl_lines="12 18"
ERROR ~ Error executing process > 'PROCESS_FILES (3)'

Caused by:
  Process `PROCESS_FILES (3)` terminated with an error exit status (127)


Command executed:

  cowpy sample3 > sample3_output.txt

Command exit status:
  127

Command output:
  (empty)

Command error:
  .command.sh: line 2: cowpy: command not found

Work dir:
  /workspaces/training/side-quests/debugging/work/82/42a5bfb60c9c6ee63ebdbc2d51aa6e

Tip: you can try to figure out what's wrong by changing to the process work directory and showing the script file named `.command.sh`

 -- Check '.nextflow.log' file for details
```

However, this workflow will not produce errors if we run it with `-stub-run`, even without the `docker` profile:

```bash
nextflow run missing_software_with_stub.nf -stub-run
```

#### Check the code

Let's examine `missing_software_with_stub.nf`:

```groovy title="missing_software.nf (with stub)" hl_lines="16-19" linenums="3"
process PROCESS_FILES {

    container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

    input:
    val sample_name

    output:
    path "${sample_name}_output.txt"

    script:
    """
    cowpy ${sample_name} > ${sample_name}_output.txt
    """

    stub:
    """
    touch ${sample_name}_output.txt
    """
}
```

Relative to `missing_software.nf`, this process has a `stub:` directive specifying a command to be used instead of the one specified in `script:`, in the event that that Nextflow is run in stub mode.

The `touch` command we're using here doesn't depend on any software or appropriate inputs, and will run in all situations, allowing us to debug workflow logic without worrying about the process internals.

**Stub running helps debug:**

- Channel structure and data flow
- Process connections and dependencies
- Parameter propagation
- Workflow logic without software dependencies

### 4.4. Systematic Debugging Approach

Now that you've learned individual debugging techniques - from trace files and work directories to preview mode, stub running, and resource monitoring - let's tie them together into a systematic methodology. Having a structured approach prevents you from getting overwhelmed by complex errors and ensures you don't miss important clues.

This methodology combines all the tools we've covered into an efficient workflow:

**Four-Phase Debugging Method:**

**Phase 1: Syntax Error Resolution (5 minutes)**

1. Check for red underlines in VSCode or your IDE
2. Run `nextflow run workflow.nf -preview` to identify syntax issues
3. Fix all syntax errors (missing braces, trailing commas, etc.)
4. Ensure the workflow parses successfully before proceeding

**Phase 2: Quick Assessment (5 minutes)**

1. Read runtime error messages carefully
2. Check if it's a runtime, logic, or resource error
3. Use preview mode to test basic workflow logic

**Phase 3: Detailed Investigation (15-30 minutes)**

1. Find the work directory of the failed task
2. Examine log files
3. Add `.view()` operators to inspect channels
4. Use `-stub-run` to test workflow logic without execution

**Phase 4: Fix and Validate (15 minutes)**

1. Make minimal targeted fixes
2. Test with resume: `nextflow run workflow.nf -resume`
3. Verify complete workflow execution

!!! tip "Using Resume for Efficient Debugging"

    Once you've identified a problem, you need an efficient way to test your fixes without wasting time re-running successful parts of your workflow. Nextflow's `-resume` functionality is invaluable for debugging.

    You will have encountered `-resume` if you've worked through [Hello Nextflow](../hello_nextflow/), and it's important that you make good use of it when debugging to save yourself waiting while the processes before your problem process run.

    **Resume debugging strategy:**

    1. Run workflow until failure
    2. Examine work directory for failed task
    3. Fix the specific issue
    4. Resume to test only the fix
    5. Repeat until workflow completes

#### Debugging Configuration Profile

To make this systematic approach even more efficient, you can create a dedicated debugging configuration that automatically enables all the tools you need:

```groovy title="nextflow.config (debug profile)" linenums="1"
profiles {
    debug {
        process {
            debug = true
            cleanup = false

            // Conservative resources for debugging
            maxForks = 1
            memory = '2.GB'
            cpus = 1
        }
    }
}
```

Then you can run the pipeline with this profile enabled:

```bash
nextflow run workflow.nf -profile debug
```

### 4.5. Practical Debugging Exercise

Now it's time to put the systematic debugging approach into practice. The workflow `buggy_workflow.nf` contains several common errors that represent the types of issues you'll encounter in real-world development.

!!! exercise

    Use the systematic debugging approach to identify and fix all errors in `buggy_workflow.nf`. This workflow attempts to process sample data from a CSV file but contains multiple intentional bugs representing common debugging scenarios.

    Start by running the workflow to see the first error:

    ```bash
    nextflow run buggy_workflow.nf
    ```

    Apply the four-phase debugging method you've learned:

    **Phase 1: Syntax Error Resolution**
    - Check for red underlines in VSCode or your IDE
    - Run `nextflow run workflow.nf -preview` to identify syntax issues
    - Fix all syntax errors (missing braces, trailing commas, etc.)
    - Ensure the workflow parses successfully before proceeding

    **Phase 2: Quick Assessment**
    - Read runtime error messages carefully
    - Identify whether errors are runtime, logic, or resource-related
    - Use `-preview` mode to test basic workflow logic

    **Phase 3: Detailed Investigation**
    - Examine work directories for failed tasks
    - Add `.view()` operators to inspect channels
    - Check log files in work directories
    - Use `-stub-run` to test workflow logic without execution

    **Phase 4: Fix and Validate**
    - Make targeted fixes
    - Use `-resume` to test fixes efficiently
    - Verify complete workflow execution

    **Debugging Tools at Your Disposal:**
    ```bash
    # Preview mode for syntax checking
    nextflow run buggy_workflow.nf -preview

    # Debug profile for detailed output
    nextflow run buggy_workflow.nf -profile debug

    # Stub running for logic testing
    nextflow run buggy_workflow.nf -stub-run

    # Resume after fixes
    nextflow run buggy_workflow.nf -resume
    ```

    ??? solution
        The `buggy_workflow.nf` contains 9 or 10 distinct errors (depending how you count) covering all major debugging categories. Here's a systematic breakdown of each error and how to fix it

        Let's start with those syntax errors:

        **Error 1: Syntax Error - Trailing Comma**
        ```groovy linenums="21"
        output:
            path "${sample_id}_result.txt",  // ERROR: Trailing comma
        ```
        **Fix:** Remove the trailing comma
        ```groovy linenums="21"
        output:
            path "${sample_id}_result.txt"
        ```

        **Error 2: Syntax Error - Missing Closing Brace**
        ```groovy linenums="24"
        script:
        """
        echo "Processing: ${sample}"
        cat ${input_file} > ${sample}_result.txt
        """
        // ERROR: Missing closing brace for processFiles process
        ```
        **Fix:** Add the missing closing brace
        ```groovy linenums="29"
        """
        echo "Processing: ${sample_id}"
        cat ${input_file} > ${sample_id}_result.txt
        """
        }  // Add missing closing brace
        ```

        **Error 3: Variable Name Error**
        ```groovy linenums="26"
        echo "Processing: ${sample}"     // ERROR: should be sample_id
        cat ${input_file} > ${sample}_result.txt  // ERROR: should be sample_id
        ```
        **Fix:** Use the correct input variable name
        ```groovy linenums="26"
        echo "Processing: ${sample_id}"
        cat ${input_file} > ${sample_id}_result.txt
        ```

        **Error 4: Undefined Variable Error**
        ```groovy linenums="87"
        heavy_ch = heavyProcess(sample_ids)  // ERROR: sample_ids undefined
        ```
        **Fix:** Use the correct channel and extract sample IDs
        ```groovy linenums="87"
        heavy_ch = heavyProcess(input_ch)
        ```

        At this point the workflow will run, but we'll still be getting errors (e.g. `Path value cannot be null` in `processFiles`), caused by bad channel structure.

        **Error 5: Channel Structure Error - Wrong Map Output**
        ```groovy linenums="83"
        .map { row -> row.sample_id }  // ERROR: processFiles expects tuple
        ```
        **Fix:** Return the tuple structure that processFiles expects
        ```groovy linenums="83"
        .map { row -> [row.sample_id, file(row.fastq_path)] }
        ```

        But this will break our for for running `heavyProcess()` above, so we'll need to use a map to pass just the sample IDs to that process:

        **Error 6: Bad channel structure for heavyProcess**
        ```groovy linenums="87"
        heavy_ch = heavyProcess(input_ch)  // ERROR: input_ch now has 2 elements per emission- heavyProcess only needs 1 (the first)
        ```
        **Fix:** Use the correct channel and extract sample IDs
        ```groovy linenums="87"
        heavy_ch = heavyProcess(input_ch.map{it[0]})
        ```

        Now we get a but further but receive an error about `No such variable: i`, because we didn't escape a Bash variable.

        **Error 7: Bash Variable Escaping Error**
        ```groovy linenums="48"
        echo "Heavy computation $i for ${sample_id}"  // ERROR: $i not escaped
        ```
        **Fix:** Escape the bash variable
        ```groovy linenums="48"
        echo "Heavy computation \${i} for ${sample_id}"
        ```

        Now we get `Process exceeded running time limit (1ms)`, so we fix the run time limit for the relevant process:

        **Error 8: Resource Configuration Error**
        ```groovy linenums="36"
        time '1 ms'  // ERROR: Unrealistic time limit
        ```
        **Fix:** Increase to a realistic time limit
        ```groovy linenums="36"
        time '100 s'
        ```

        Next we have a `Missing output file(s)` errror to resolve:

        **Error 9: Output File Name Mismatch**
        ```groovy linenums="49"
        done > ${sample_id}.txt  // ERROR: Wrong filename, should match output declaration
        ```
        **Fix:** Match the output declaration
        ```groovy linenums="49"
        done > ${sample_id}_heavy.txt
        ```

        The first two processes ran, but not the third.

        **Error 10: Output File Name Mismatch**
        ```groovy linenums="88"
        file_ch = channel.fromPath("*.txt") // Error: attempting to take input from the pwd rather than a process
        handleFiles(file_ch)
        ```
        **Fix:** Take the output from the previous process
        ```groovy linenums="88"
        handleFiles(heavyProcess.out)
        ```

        With that, the whole workflow should run.

        **Complete Corrected Workflow:**
        ```groovy linenums="1"
        #!/usr/bin/env nextflow

        /*
        * Buggy workflow for debugging exercises
        * This workflow contains several intentional bugs for learning purposes
        */

        // Parameters with missing validation
        params.input = 'data/sample_data.csv'
        params.output = 'results'

        /*
        * Process with input/output mismatch
        */
        process processFiles {
            publishDir "${params.output}/processed", mode: 'copy'

            input:
                tuple val(sample_id), path(input_file)

            output:
                path "${sample_id}_result.txt"

            script:
            """
            echo "Processing: ${sample_id}"
            cat ${input_file} > ${sample_id}_result.txt
            """
        }

        /*
        * Process with resource issues
        */
        process heavyProcess {
            publishDir "${params.output}/heavy", mode: 'copy'

            time '100 s'

            input:
                val sample_id

            output:
                path "${sample_id}_heavy.txt"

            script:
            """
            # Simulate heavy computation
            for i in {1..1000000}; do
                echo "Heavy computation \$i for ${sample_id}"
            done > ${sample_id}_heavy.txt
            """
        }

        /*
        * Process with file handling issues
        */
        process handleFiles {
            publishDir "${params.output}/files", mode: 'copy'

            input:
                path input_file

            output:
                path "processed_${input_file}"

            script:
            """
            if [ -f "${input_file}" ]; then
                cp ${input_file} processed_${input_file}
            fi
            """
        }

        /*
        * Main workflow with channel issues
        */
        workflow {

            // Channel with incorrect usage
            input_ch = channel
                .fromPath(params.input)
                .splitCsv(header: true)
                .map { row -> [row.sample_id, file(row.fastq_path)] }

            processed_ch = processFiles(input_ch)

            heavy_ch = heavyProcess(input_ch.map{it[0]})

            handleFiles(heavyProcess.out)
        }
        ```

**Error Categories Covered:**

- **Syntax errors**: Missing braces, trailing commas, undefined variables
- **Channel structure errors**: Wrong data shapes, undefined channels
- **Process errors**: Output file mismatches, variable escaping
- **Resource errors**: Unrealistic time limits

**Key Debugging Lessons:**

1. **Read error messages carefully** - they often point directly to the problem
2. **Use systematic approaches** - fix one error at a time and test with `-resume`
3. **Understand data flow** - channel structure errors are often the most subtle
4. **Check work directories** - when processes fail, the logs tell you exactly what went wrong

---

## Summary

In this side quest, you've learned a set of systematic techniques for debugging Nextflow workflows.
Applying these techniques in your own work will enable you to spend less time fighting your computer, solve problems faster and protect yourself from future issues.

### Key patterns

<!-- TODO: Can we add snippets of code below to illustrate? -->

**1. How to identify and fix syntax errors**:

- Interpreting Nextflow error messages and locating problems
- Common syntax errors: missing braces, incorrect keywords, undefined variables
- Distinguishing between Nextflow (Groovy) and Bash variables
- Using VS Code extension features for early error detection

**2. How to debug channel structure issues**:

- Understanding channel cardinality and exhaustion issues
- Debugging channel content structure mismatches
- Using `.view()` operators for channel inspection
- Recognizing error patterns like square brackets in output

**3. How to troubleshoot process execution problems**:

- Diagnosing missing output file errors
- Understanding exit codes (127 for missing software, 137 for memory issues)
- Investigating work directories and command files
- Configuring resources appropriately

**4. How to use Nextflow's built-in debugging tools**:

- Leveraging preview mode and real-time debugging
- Implementing stub running for logic testing
- Applying resume for efficient debugging cycles
- Following a four-phase systematic debugging methodology

!!! tip "Quick Debugging Reference"

    **Syntax errors?** → Check VSCode warnings, run `nextflow run workflow.nf -preview`

    **Channel issues?** → Use `.view()` to inspect content: `my_channel.view()`

    **Process failures?** → Check work directory files:

    - `.command.sh` - the executed script
    - `.command.err` - error messages
    - `.exitcode` - exit status (127 = command not found, 137 = killed)

    **Mysterious behavior?** → Run with `-stub-run` to test workflow logic

    **Made fixes?** → Use `-resume` to save time testing: `nextflow run workflow.nf -resume`

---

### Additional resources

<!-- TODO: Are there any specific articles to call out / link to? -->

Check out the [Nextflow documentation](https://www.nextflow.io/docs/latest/) for more advanced debugging features and best practices. You might want to:

- Add more comprehensive error handling to your workflows
- Write tests for edge cases and error conditions using nf-test
- Set up monitoring and logging for production workflows
- Learn about other debugging tools like profiling and performance analysis
- Explore more advanced debugging techniques for complex workflows

**Remember:** Effective debugging is a skill that improves with practice. The systematic methodology and comprehensive toolkit you've acquired here will serve you well throughout your Nextflow development journey.

---

## What's next?

Return to the [menu of Side Quests](./index.md) or click the button in the bottom right of the page to move on to the next topic in the list.
