# Part 3: Hello Workflow

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/_aO56V3iXGI?si=y8lAedhEHWaTV4zd&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } See [the whole playlist](https://youtube.com/playlist?list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&si=eF7cLR62goy-lc6n) on the Nextflow YouTube channel.

:green_book: The video transcript is available [here](./transcripts/03_hello_workflow.md).
///

Most real-world workflows involve more than one step.
In this training module, you'll learn how to connect processes together in a multi-step workflow.

This will teach you the Nextflow way of achieving the following:

1. Making data flow from one process to the next
2. Running a process written in another language
3. Passing more than one input to a process
4. Collecting outputs from multiple process calls into a single process call
5. Handling multiple outputs coming out of a process

To demonstrate, we will continue building on the Hello World example from Parts 1 and 2, and make the following changes to our workflow:

1. Add a second step, written in Python, that translates each greeting into another language.
2. Add a third step that collects all the translated greetings and writes them into a single file named after the batch and the language.
3. Make the collection step also report a simple statistic about what was processed.

Here is where this part fits in the course:

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-course-map-part3.svg"
</figure>

??? info "How to begin from this section"

    This section of the course assumes you have completed Parts 1-2 of the [Hello Nextflow](./index.md) course, but if you are comfortable with the basics covered in those sections, you can start from here without doing anything special.

---

## 0. Warmup: Run `hello-workflow.nf`

We're going to use the workflow script `hello-workflow.nf` as a starting point.
It is equivalent to the script produced by working through Part 2 of this training course, except we've removed the `view()` statements.

This diagram summarizes the current operation of the workflow.
It should look familiar, except now we're explicitly showing that the outputs of the process are packaged in a channel, just like the inputs were.
We're going to put that output channel to good use in a minute.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-workflow-channels.svg"
</figure>

Just to make sure everything is working, run the script once before making any changes:

```bash
nextflow run hello-workflow.nf
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `hello-workflow.nf` [mighty_becquerel] revision: 42ced9ea2d

    executor >  local (3)
    [1f/03f176] sayHello (Celine) | 3 of 3 ✔

    Outputs:

      /workspaces/training/hello-nextflow/results

      greetings:
        - hello/Raphael.txt
        - hello/Mateo.txt
        - hello/Celine.txt
    ```

As previously, you will find the output files in the `results/hello` directory, alongside the outputs of the previous parts.

??? abstract "Directory contents"

    ```console
    results
    ├── hello
    │   ├── Celine.txt
    │   ├── Mateo.txt
    │   ├── output.txt
    │   └── Raphael.txt
    └── output.txt -> /workspaces/training/hello-nextflow/work/.../output.txt
    ```

If that worked for you, you're ready to learn how to assemble a multi-step workflow.

---

## 1. Add a translation step to the workflow

We're going to add a step that translates each greeting into another language, for example from `Hello, Mateo!` to `Hola, Mateo!`.
This step is a small Python script, which shows that a process can run code in other languages than the shell.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-multistep.svg"
</figure>

We'll build it in three steps:

1. Test the translation code in the terminal.
2. Wrap it in a process that translates into Spanish, and connect it to the output of `sayHello`.
3. Make the language an input of the process, set by a parameter.

### 1.1. Test the translation code in the terminal

The translation is a simple text replacement: read the greeting, replace `Hello` with the word for hello in the target language, and write the result to a new file.

To test it, create a greeting file the same way `sayHello` does, then run the replacement with `python3`:

```bash
echo 'Hello, Mateo!' > Mateo.txt
python3 -c "greeting = open('Mateo.txt').read(); open('es-Mateo.txt', 'w').write(greeting.replace('Hello', 'Hola'))"
```

The output is a text file called `es-Mateo.txt` that contains the translated greeting.

??? abstract "File contents"

    ```console title="es-Mateo.txt"
    Hola, Mateo!
    ```

That's what we're going to do in our workflow, for every greeting.

### 1.2. Write the translation step as a Nextflow process

We start with a version that always translates into Spanish, and connect it to the first process.

#### 1.2.1. Define the process

Add the following process definition to the workflow script, right under the `sayHello` process:

```groovy title="hello-workflow.nf" linenums="22"
/*
 * Use Python to translate the greeting into Spanish
 */
process translateGreeting {

    tag "${input_file}"

    input:
    path input_file

    output:
    path "es-${input_file}"

    script:
    """
    #!/usr/bin/env python3

    with open('${input_file}') as greeting_file:
        greeting = greeting_file.read()

    with open('es-${input_file}', 'w') as translated_file:
        translated_file.write(greeting.replace('Hello', 'Hola'))
    """
}
```

The structure is the same as for `sayHello`, with two differences.
The `tag` directive, introduced in Part 2, labels each task with the name of the file it translates.

First, the input is a file, so it uses the `path` qualifier, and we compose the output filename from the input filename, similarly to what we did for the output of the first process.

Second, the script is written in Python.
The first line of the script, `#!/usr/bin/env python3`, is called a shebang: it tells the system which program should run the script.
Without it, Nextflow runs the script with the shell, as it did for `sayHello`.
With it, you can write the script of a process in any language that is available in the environment, such as Python, R or Perl.

The `${input_file}` placeholders work the same way as in a shell script: Nextflow replaces them with the actual file name before running the script.
Python never sees the placeholders, only the values.

??? tip "Longer scripts"

    Writing the code inline is convenient for a few lines.
    For a longer script, you can save it as a file in a directory named `bin/` next to your workflow script, make it executable, and call it by name in the `script` block.
    Nextflow adds the `bin/` directory to the `PATH` of every task automatically.

#### 1.2.2. Pass the output of the first process to the new process

Now we need to call the new process in the workflow block, and make the output of the `sayHello()` process flow into it.

Conveniently, Nextflow automatically packages the output of a process into a channel, as shown in the diagram in the warmup section.
We can refer to the output channel of a process as `<process>.out`.

So the output of the `sayHello` process is a channel called `sayHello.out`, which we can plug straight into the call to `translateGreeting()`.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-multistep-connector.svg"
</figure>

In the workflow block, make the following code change:

=== "After"

    ```groovy title="hello-workflow.nf" linenums="54" hl_lines="10-11"
    workflow {

        main:
        // create a channel for inputs from a CSV file
        name_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { row -> row[0] }
        // emit a greeting
        sayHello(name_ch)
        // translate the greeting
        translateGreeting(sayHello.out)

        publish:
        greetings = sayHello.out
    }
    ```

=== "Before"

    ```groovy title="hello-workflow.nf" linenums="54"
    workflow {

        main:
        // create a channel for inputs from a CSV file
        name_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { row -> row[0] }
        // emit a greeting
        sayHello(name_ch)

        publish:
        greetings = sayHello.out
    }
    ```

For a simple case like this (one output to one input), that's all we need to do to connect two processes!

#### 1.2.3. Set up the workflow output publishing

Finally, update the workflow outputs to publish the results from the second process too.

In the `publish:` section of the `workflow` block, make the following code change:

=== "After"

    ```groovy title="hello-workflow.nf" linenums="66" hl_lines="3"
        publish:
        greetings = sayHello.out
        translated = translateGreeting.out
    }
    ```

=== "Before"

    ```groovy title="hello-workflow.nf" linenums="66"
        publish:
        greetings = sayHello.out
    }
    ```

Then, in the `output` block, make the following code change:

=== "After"

    ```groovy title="hello-workflow.nf" linenums="71" hl_lines="6-9"
    output {
        greetings {
            path 'hello'
            mode 'copy'
        }
        translated {
            path 'hello'
            mode 'copy'
        }
    }
    ```

=== "Before"

    ```groovy title="hello-workflow.nf" linenums="71"
    output {
        greetings {
            path 'hello'
            mode 'copy'
        }
    }
    ```

Each output has its own entry, so you can set a different path for each one.
For now they all go to the same `hello` subdirectory; in section 2.5, you'll organize them by batch.

#### 1.2.4. Run the workflow with `-resume`

Run the workflow with the `-resume` flag, which reuses the results of steps that already ran successfully ([Part 1 of Nextflow Run](../nextflow_run/01_run_nextflow.md) covers it in detail):

```bash
nextflow run hello-workflow.nf -resume
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `hello-workflow.nf` [desperate_baekeland] revision: 18c770fedd

    executor >  local (3)
    [1f/03f176] sayHello (Celine)               | 3 of 3, cached: 3 ✔
    [f7/fb36e4] translateGreeting (Raphael.txt) | 3 of 3 ✔

    Outputs:

      /workspaces/training/hello-nextflow/results

      greetings:
        - hello/Raphael.txt
        - hello/Celine.txt
        - hello/Mateo.txt

      translated:
        - hello/es-Celine.txt
        - hello/es-Mateo.txt
        - hello/es-Raphael.txt
    ```

There is now an extra line in the console output that corresponds to the new process we just added, and the `sayHello` tasks were retrieved from the cache.

You'll find the translated greetings in the `results/hello` directory, as set in the `output` block.

??? abstract "File contents"

    ```console title="results/hello/es-Mateo.txt"
    Hola, Mateo!
    ```

It's also worth taking a look inside the work directory of one of the calls to the second process.

??? abstract "Directory contents"

    ```console
    work/2a/ab0718b375f74cf4c244ffcb0399ca/
    ├── Mateo.txt -> /workspaces/training/hello-nextflow/work/cf/a496aeed43fb0da1f62ca9ebedbb18/Mateo.txt
    └── es-Mateo.txt
    ```

Notice there are two output files: the output of the first process as well as the output of the second.

The output of the first process is in there because Nextflow **staged** it there in order to have everything needed for execution within the same subdirectory.
It is a symbolic link pointing to the original file in the subdirectory of the first process call.
By default, when running on a single machine as we're doing here, Nextflow uses symbolic links rather than copies to stage input and intermediate files.

The `.command.sh` file in the same directory shows the Python script that Nextflow ran, with the file names filled in:

```python title="work/2a/ab0718b375f74cf4c244ffcb0399ca/.command.sh"
#!/usr/bin/env python3

with open('Mateo.txt') as greeting_file:
    greeting = greeting_file.read()

with open('es-Mateo.txt', 'w') as translated_file:
    translated_file.write(greeting.replace('Hello', 'Hola'))
```

All we did was connect the output of `sayHello` to the input of `translateGreeting`, and Nextflow handled the individual files and passed them between the two steps for us.
This is one of the reasons Nextflow channels are so powerful: they take care of the busywork involved in connecting workflow steps together.

### 1.3. Make the language an input of the process

Our translation step always translates into Spanish.
To make it translate into other languages, we'll give the process a second input, `language`, set by a pipeline parameter.

#### 1.3.1. Add the `language` input and a dictionary of translations

Good news: we can declare as many inputs as we want in the process definition.

In the `translateGreeting` process definition, make the following code changes:

=== "After"

    ```groovy title="hello-workflow.nf" linenums="22" hl_lines="2 6 10 13 18 20-28 30-31 36-37"
    /*
     * Use Python to translate the greeting into another language
     */
    process translateGreeting {

        tag "${language}-${input_file}"

        input:
        path input_file
        val language

        output:
        path "${language}-${input_file}"

        script:
        """
        #!/usr/bin/env python3
        import sys

        translations = {
            'es': 'Hola',        # Spanish
            'fr': 'Bonjour',     # French
            'it': 'Ciao',        # Italian
            'nl': 'Hallo',       # Dutch
            'ar': 'Marhaban',    # Arabic
            'ja': 'Konnichiwa',  # Japanese
            'en': 'Hello',       # English
        }

        if '${language}' not in translations:
            sys.exit('Unsupported language code: ${language}')

        with open('${input_file}') as greeting_file:
            greeting = greeting_file.read()

        with open('${language}-${input_file}', 'w') as translated_file:
            translated_file.write(greeting.replace('Hello', translations['${language}']))
        """
    }
    ```

=== "Before"

    ```groovy title="hello-workflow.nf" linenums="22"
    /*
     * Use Python to translate the greeting into Spanish
     */
    process translateGreeting {

        tag "${input_file}"

        input:
        path input_file

        output:
        path "es-${input_file}"

        script:
        """
        #!/usr/bin/env python3

        with open('${input_file}') as greeting_file:
            greeting = greeting_file.read()

        with open('es-${input_file}', 'w') as translated_file:
            translated_file.write(greeting.replace('Hello', 'Hola'))
        """
    }
    ```

The new input, `language`, is a plain value, so it uses the `val` qualifier, like the input of `sayHello`.

The script now looks up the word for hello in a Python dictionary called `translations`, using the language as the key.
The keys are two-letter language codes, such as `es` for Spanish and `ja` for Japanese.
Nextflow fills in `${language}` before Python runs, so for Spanish, Python sees `translations['es']`.
If the code is not in the dictionary, the script exits with an error message (`sys.exit`) instead of a Python traceback.
The output file name also starts with the language code, so the outputs of different languages don't overwrite each other.
The tag now includes the language code too, matching the output file name; the language alone would give all the tasks of a run the same tag.

#### 1.3.2. Add a `language` parameter

You already know how to use the `params` system to declare CLI parameters.
Declare a `language` parameter, with Spanish (`es`) as the default value.

In the pipeline parameters section, make the following code change:

=== "After"

    ```groovy title="hello-workflow.nf" linenums="62" hl_lines="6"
    /*
     * Pipeline parameters
     */
    params {
        input: Path = 'data/names.csv'
        language: String = 'es'
    }
    ```

=== "Before"

    ```groovy title="hello-workflow.nf" linenums="62"
    /*
     * Pipeline parameters
     */
    params {
        input: Path = 'data/names.csv'
    }
    ```

#### 1.3.3. Pass the language to the process

To provide the value of the parameter to the process, add it to the process call.

In the workflow block, make the following code change:

=== "After"

    ```groovy title="hello-workflow.nf" linenums="79" hl_lines="2"
        // translate the greeting
        translateGreeting(sayHello.out, params.language)
    ```

=== "Before"

    ```groovy title="hello-workflow.nf" linenums="79"
        // translate the greeting
        translateGreeting(sayHello.out)
    ```

To provide multiple inputs to a process, you simply list them in the call parentheses, separated by commas.

!!! warning

    You MUST provide the inputs to the process in the EXACT SAME ORDER as they are listed in the input definition block of the process.

Notice that the two inputs are different kinds of things.
`sayHello.out` is a channel that carries three files, one per greeting, while `params.language` is a single value.
Nextflow calls the process once per file in the channel, and reuses the single value for every call.
So each of the three greetings gets translated into the same language.

??? note "What if the language comes from a channel?"

    The value is reused only because `params.language` is a plain value.
    If you put it in a channel instead, for example with `channel.of('es')`, it is consumed by the first call: the process runs only once, on the first greeting.
    Nextflow pairs up the elements of channel inputs one by one, and stops when one of the channels runs out.

#### 1.3.4. Run the workflow

Run the workflow again with `-resume`:

```bash
nextflow run hello-workflow.nf -resume
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `hello-workflow.nf` [pedantic_euler] revision: 861ea927ab

    executor >  local (3)
    [cf/a496ae] sayHello (Mateo)                 | 3 of 3, cached: 3 ✔
    [b8/7fc215] translateGreeting (es-Mateo.txt) | 3 of 3 ✔

    Outputs:

      /workspaces/training/hello-nextflow/results

      greetings:
        - hello/Celine.txt
        - hello/Mateo.txt
        - hello/Raphael.txt

      translated:
        - hello/es-Celine.txt
        - hello/es-Raphael.txt
        - hello/es-Mateo.txt
    ```

The results are the same as before, since the default language is Spanish.
However, the `translateGreeting` tasks ran again instead of being retrieved from the cache, because we changed the code of the process.

Now try another language:

```bash
nextflow run hello-workflow.nf -resume --language fr
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `hello-workflow.nf` [distracted_kalman] revision: 861ea927ab

    executor >  local (3)
    [cf/a496ae] sayHello (Mateo)                  | 3 of 3, cached: 3 ✔
    [54/7a4a37] translateGreeting (fr-Mateo.txt)  | 3 of 3 ✔

    Outputs:

      /workspaces/training/hello-nextflow/results

      greetings:
        - hello/Raphael.txt
        - hello/Mateo.txt
        - hello/Celine.txt

      translated:
        - hello/fr-Celine.txt
        - hello/fr-Raphael.txt
        - hello/fr-Mateo.txt
    ```

The `sayHello` tasks were cached, and the `translateGreeting` tasks ran again with the new input.

??? abstract "File contents"

    ```console title="results/hello/fr-Mateo.txt"
    Bonjour, Mateo!
    ```

Finally, try a language code that is not in the dictionary, such as `de` for German:

```bash
nextflow run hello-workflow.nf -resume --language de
```

??? failure "Command output"

    ```console hl_lines="44-45"
    N E X T F L O W   ~  version 26.04.4

    Launching `hello-workflow.nf` [astonishing_jones] revision: 861ea927ab

    executor >  local (3)
    [1f/03f176] sayHello (Celine)                  | 3 of 3, cached: 3 ✔
    [94/cb45cc] translateGreeting (de-Raphael.txt) | 0 of 3 ✘
    ERROR ~ Error executing process > 'translateGreeting (de-Mateo.txt)'

    Caused by:
      Process `translateGreeting (de-Mateo.txt)` terminated with an error exit status (1)


    Command executed:

      #!/usr/bin/env python3
      import sys

      translations = {
          'es': 'Hola',        # Spanish
          'fr': 'Bonjour',     # French
          'it': 'Ciao',        # Italian
          'nl': 'Hallo',       # Dutch
          'ar': 'Marhaban',    # Arabic
          'ja': 'Konnichiwa',  # Japanese
          'en': 'Hello',       # English
      }

      if 'de' not in translations:
          sys.exit('Unsupported language code: de')

      with open('Mateo.txt') as greeting_file:
          greeting = greeting_file.read()

      with open('de-Mateo.txt', 'w') as translated_file:
          translated_file.write(greeting.replace('Hello', translations['de']))

    Command exit status:
      1

    Command output:
      (empty)

    Command error:
      Unsupported language code: de

    Work dir:
      /workspaces/training/hello-nextflow/work/9e/94d39bec1a83bf0013f2fc54436a4b

    Tip: you can try to figure out what's wrong by changing to the process work dir and showing the script file named `.command.sh`

    -- Check '.nextflow.log' file for details
    ```

The run fails, and the error report shows the script as Nextflow ran it, with `de` filled in.
The message from `sys.exit` appears under `Command error:`, since it was written to standard error.
Nextflow also saves it in the `.command.err` file of the task directory listed under `Work dir:`.

### Takeaway

You know how to chain processes together by providing the output of one step as input to the next step, how to write the script of a process in another language, and how to pass more than one input to a process.

### What's next?

Learn how to collect outputs from batched process calls and feed them into a single process.

---

## 2. Add a third step to collect all the greetings

When we use a process to apply a transformation to each of the elements in a channel, like we're doing here to the multiple greetings, we sometimes want to collect elements from the output channel of that process, and feed them into another process that performs some kind of analysis or summation.

To demonstrate, we'll add a new step to our pipeline that collects all the translated greetings produced by the `translateGreeting` process and writes them out to a single file.
We'll also name that file after a batch name provided as a parameter and after the language, so that the results of different batches of inputs, and of different languages, don't overwrite each other.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-collect.svg"
</figure>

Not to spoil the surprise, but this is going to involve a very useful operator.

### 2.1. Define the collection command and test it in the terminal

The collection step we want to add to our workflow will use the `cat` command to concatenate multiple translated greetings into a single file.

Run the following in your terminal to verify that it works as expected:

```bash
echo 'Hola, Mateo!' > es-Mateo.txt
echo 'Hola, Raphael!' > es-Raphael.txt
echo 'Hola, Celine!' > es-Celine.txt
cat es-Mateo.txt es-Raphael.txt es-Celine.txt > greetings.txt
```

The output is a text file called `greetings.txt` that contains all three greetings.

??? abstract "File contents"

    ```console title="greetings.txt"
    Hola, Mateo!
    Hola, Raphael!
    Hola, Celine!
    ```

That is the result we want to achieve with our workflow.

### 2.2. Create a new process to do the collection step

Add the following process definition to the workflow script, right under the `translateGreeting` process:

```groovy title="hello-workflow.nf" linenums="62"
/*
 * Collect translated greetings into a single output file
 */
process collectGreetings {

    tag "${batch_name}"

    input:
    path input_files
    val batch_name
    val language

    output:
    path "${batch_name}-${language}-greetings.txt"

    script:
    """
    cat ${input_files} > '${batch_name}-${language}-greetings.txt'
    """
}
```

This process has three inputs:

- `input_files` will hold the files to collect.
  It uses the `path` qualifier even though we expect it to contain multiple files.
- `batch_name` is the name of the batch, which we use in the output file name, so that each batch gets its own output file.
- `language` is the language of the greetings, which we also use in the output file name.
  A batch names a group of samples, which you may translate into several languages; including the language in the file name keeps the results of different languages apart.

`batch_name` and `language` are plain values, so they use the `val` qualifier, like the `language` input of `translateGreeting` in section 1.3.
The process runs only once per run, but tagging it with the batch name makes runs of different batches easy to tell apart.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-collect-batch.svg"
</figure>

The script needs to handle an arbitrary number of input files, so we can't write the command up front.
In other words, if `input_files` contains `[file1.txt, file2.txt, file3.txt]`, we need Nextflow to turn that into `cat file1.txt file2.txt file3.txt`.

Fortunately, Nextflow is quite happy to do that for us if we simply write `cat ${input_files}` in the script command.

!!! tip

    Some command-line tools require providing an argument (like `-input`) for each input file.
    In that case, we would have to do a little bit of extra work to compose the command.
    You can see an example of this in the [Nextflow for Genomics](../nf4_science/genomics/index.md) training course.

### 2.3. Add the collection step to the workflow

Now we need a `batch` parameter, and a call to the collection process on the output of the translation step.
That output is also a channel, called `translateGreeting.out`.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-collect-connector.svg"
</figure>

#### 2.3.1. Add a `batch` parameter

Declare a `batch` parameter, with a default value, in the pipeline parameters section:

=== "After"

    ```groovy title="hello-workflow.nf" linenums="83" hl_lines="7"
    /*
     * Pipeline parameters
     */
    params {
        input: Path = 'data/names.csv'
        language: String = 'es'
        batch: String = 'batch'
    }
    ```

=== "Before"

    ```groovy title="hello-workflow.nf" linenums="83"
    /*
     * Pipeline parameters
     */
    params {
        input: Path = 'data/names.csv'
        language: String = 'es'
    }
    ```

#### 2.3.2. Connect the process calls

In the workflow block, make the following code change:

=== "After"

    ```groovy title="hello-workflow.nf" linenums="101" hl_lines="3-4"
        // translate the greeting
        translateGreeting(sayHello.out, params.language)
        // collect all the greetings into one file
        collectGreetings(translateGreeting.out, params.batch, params.language)
    ```

=== "Before"

    ```groovy title="hello-workflow.nf" linenums="101"
        // translate the greeting
        translateGreeting(sayHello.out, params.language)
    ```

This connects the output of `translateGreeting()` to the first input of `collectGreetings()`, and passes the batch name and the language as the second and third inputs.
The language comes from the same `params.language` parameter as for `translateGreeting()`, so the collected file is named after the language the greetings were translated into.

#### 2.3.3. Update the workflow outputs

In the `publish:` section of the `workflow` block, make the following code change:

=== "After"

    ```groovy title="hello-workflow.nf" linenums="106" hl_lines="4"
        publish:
        greetings = sayHello.out
        translated = translateGreeting.out
        collected = collectGreetings.out
    }
    ```

=== "Before"

    ```groovy title="hello-workflow.nf" linenums="106"
        publish:
        greetings = sayHello.out
        translated = translateGreeting.out
    }
    ```

Then, in the `output` block, make the following code change:

=== "After"

    ```groovy title="hello-workflow.nf" linenums="112" hl_lines="10-13"
    output {
        greetings {
            path 'hello'
            mode 'copy'
        }
        translated {
            path 'hello'
            mode 'copy'
        }
        collected {
            path 'hello'
            mode 'copy'
        }
    }
    ```

=== "Before"

    ```groovy title="hello-workflow.nf" linenums="112"
    output {
        greetings {
            path 'hello'
            mode 'copy'
        }
        translated {
            path 'hello'
            mode 'copy'
        }
    }
    ```

#### 2.3.4. Run the workflow with `-resume`

Let's try it.

```bash
nextflow run hello-workflow.nf -resume
```

??? success "Command output"

    ```console hl_lines="8"
    N E X T F L O W   ~  version 26.04.4

    Launching `hello-workflow.nf` [stoic_banach] revision: 512c90ea42

    executor >  local (3)
    [1f/03f176] sayHello (Celine)                  | 3 of 3, cached: 3 ✔
    [52/df8eff] translateGreeting (es-Raphael.txt) | 3 of 3, cached: 3 ✔
    [8a/eb1442] collectGreetings (batch)           | 3 of 3 ✔

    Outputs:

      /workspaces/training/hello-nextflow/results

      greetings:
        - hello/Raphael.txt
        - hello/Celine.txt
        - hello/Mateo.txt

      translated:
        - hello/es-Celine.txt
        - hello/es-Mateo.txt
        - hello/es-Raphael.txt

      collected:
        - hello/batch-es-greetings.txt
        - hello/batch-es-greetings.txt
        - hello/batch-es-greetings.txt
    ```

It runs successfully, including the third step.

However, look at the number of calls for `collectGreetings()` on the last line.
We were only expecting one, but there are three.

Now have a look at the contents of the final output file.

??? abstract "File contents"

    ```console title="results/hello/batch-es-greetings.txt"
    Hola, Celine!
    ```

Oh no. The collection step was run individually on each greeting, which is NOT what we wanted.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-collect-no-operator.svg"
</figure>

We need to do something to tell Nextflow explicitly that we want that third step to run on all the elements in the channel output by `translateGreeting()`.

### 2.4. Use an operator to collect the greetings into a single input

Yes, once again the answer to our problem is an operator.

Specifically, we are going to use the aptly-named [`collect()`](https://nextflow.io/docs/latest/reference/operator.html#collect) operator.

#### 2.4.1. Add the `collect()` operator

This time it's going to look a bit different because we're not adding the operator in the context of a channel factory; we're adding it to an output channel.

We take `translateGreeting.out` and append the `collect()` operator, which gives us `translateGreeting.out.collect()`.
We can plug that directly into the `collectGreetings()` process call.

In the workflow block, make the following code change:

=== "After"

    ```groovy title="hello-workflow.nf" linenums="103" hl_lines="2"
        // collect all the greetings into one file
        collectGreetings(translateGreeting.out.collect(), params.batch, params.language)
    ```

=== "Before"

    ```groovy title="hello-workflow.nf" linenums="103"
        // collect all the greetings into one file
        collectGreetings(translateGreeting.out, params.batch, params.language)
    ```

#### 2.4.2. Add some `view()` statements

Let's also include a couple of `view()` statements to visualize the before and after states of the channel contents.

=== "After"

    ```groovy title="hello-workflow.nf" linenums="103" hl_lines="4-6"
        // collect all the greetings into one file
        collectGreetings(translateGreeting.out.collect(), params.batch, params.language)

        // optional view statements
        translateGreeting.out.view { filepath -> "Before collect: $filepath" }
        translateGreeting.out.collect().view { filepaths -> "After collect: $filepaths" }
    ```

=== "Before"

    ```groovy title="hello-workflow.nf" linenums="103"
        // collect all the greetings into one file
        collectGreetings(translateGreeting.out.collect(), params.batch, params.language)
    ```

The `view()` statements can go anywhere you want; we put them right after the call for readability.
As in Part 2, the closure variables are named after what the channel contains: one file path per element before `collect()`, and a list of all the file paths after it.

#### 2.4.3. Run the workflow again with `-resume`

Let's try it:

```bash
nextflow run hello-workflow.nf -resume
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `hello-workflow.nf` [lonely_nightingale] revision: 98ed72b5e7

    executor >  local (1)
    [1f/03f176] sayHello (Celine)                  | 3 of 3, cached: 3 ✔
    [52/df8eff] translateGreeting (es-Raphael.txt) | 3 of 3, cached: 3 ✔
    [19/085ee1] collectGreetings (batch)           | 1 of 1 ✔
    Before collect: /workspaces/training/hello-nextflow/work/b8/7fc215ebc8eabfe57b3e3e07de1ceb/es-Mateo.txt
    Before collect: /workspaces/training/hello-nextflow/work/52/df8effab27c7ac0c7e9a4c9ebb7072/es-Raphael.txt
    Before collect: /workspaces/training/hello-nextflow/work/00/c5c01f2a09b6304128c99dd11fd228/es-Celine.txt
    After collect: [/workspaces/training/hello-nextflow/work/b8/7fc215ebc8eabfe57b3e3e07de1ceb/es-Mateo.txt, /workspaces/training/hello-nextflow/work/52/df8effab27c7ac0c7e9a4c9ebb7072/es-Raphael.txt, /workspaces/training/hello-nextflow/work/00/c5c01f2a09b6304128c99dd11fd228/es-Celine.txt]

    Outputs:

      /workspaces/training/hello-nextflow/results

      greetings:
        - hello/Raphael.txt
        - hello/Celine.txt
        - hello/Mateo.txt

      translated:
        - hello/es-Mateo.txt
        - hello/es-Raphael.txt
        - hello/es-Celine.txt

      collected: hello/batch-es-greetings.txt
    ```

It runs successfully, although the log output may look a little messier than this (we cleaned it up for readability).

This time the third step was only called once!
Looking at the output of the `view()` statements, we see the following:

- Three `Before collect:` statements, one for each greeting: at that point the file paths are individual items in the channel.
- A single `After collect:` statement: the three file paths are now packaged into a single element.

We can summarize that with the following diagram:

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-collect-WITH-operator.svg"
</figure>

Finally, you can have a look at the contents of the output file to satisfy yourself that everything worked correctly.

??? abstract "File contents"

    ```console title="results/hello/batch-es-greetings.txt"
    Hola, Celine!
    Hola, Raphael!
    Hola, Mateo!
    ```

This time we have all three greetings in the final output file. Success!

!!! note

    If you run this several times without `-resume`, you will see that the order of the greetings changes from one run to the next.
    This shows you that the order in which elements flow through process calls is not guaranteed to be consistent.

#### 2.4.4. Remove `view()` statements for readability

Before you move on to the next section, we recommend you delete the `view()` statements to avoid cluttering the console output.

=== "After"

    ```groovy title="hello-workflow.nf" linenums="103"
        // collect all the greetings into one file
        collectGreetings(translateGreeting.out.collect(), params.batch, params.language)
    ```

=== "Before"

    ```groovy title="hello-workflow.nf" linenums="103" hl_lines="4-6"
        // collect all the greetings into one file
        collectGreetings(translateGreeting.out.collect(), params.batch, params.language)

        // optional view statements
        translateGreeting.out.view { filepath -> "Before collect: $filepath" }
        translateGreeting.out.collect().view { filepaths -> "After collect: $filepaths" }
    ```

This is basically the reverse operation from point 2.4.2.

### 2.5. Use the batch name to organize the outputs

The `batch` parameter is also a good way to keep the published results of different batches apart.
Each output's `path` can include parameter values, so you can publish each batch to its own subdirectory, and group the intermediate files in a subdirectory of their own.

#### 2.5.1. Update the output paths

In the `output` block, make the following code change:

=== "After"

    ```groovy title="hello-workflow.nf" linenums="112" hl_lines="3 7 11"
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
    }
    ```

=== "Before"

    ```groovy title="hello-workflow.nf" linenums="112"
    output {
        greetings {
            path 'hello'
            mode 'copy'
        }
        translated {
            path 'hello'
            mode 'copy'
        }
        collected {
            path 'hello'
            mode 'copy'
        }
    }
    ```

The double quotes are required for Nextflow to insert the parameter value into the path, as in a process script.

This is how a pipeline author decides how the results are organized.
The person running the pipeline can change the location of the output directory, as [Part 2 of Nextflow Run](../nextflow_run/02_configure_pipeline.md) shows, but not the organization within it, except through parameters like `batch` that the author builds into the paths.

#### 2.5.2. Run the workflow with a batch name

Run the workflow with a batch name on the command line:

```bash
nextflow run hello-workflow.nf -resume --batch trio
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `hello-workflow.nf` [scruffy_franklin] revision: f3e169ca36

    executor >  local (1)
    [1f/03f176] sayHello (Celine)                 | 3 of 3, cached: 3 ✔
    [00/c5c01f] translateGreeting (es-Celine.txt) | 3 of 3, cached: 3 ✔
    [7d/4c5007] collectGreetings (trio)           | 1 of 1 ✔

    Outputs:

      /workspaces/training/hello-nextflow/results

      greetings:
        - trio/intermediates/Raphael.txt
        - trio/intermediates/Celine.txt
        - trio/intermediates/Mateo.txt

      translated:
        - trio/intermediates/es-Mateo.txt
        - trio/intermediates/es-Raphael.txt
        - trio/intermediates/es-Celine.txt

      collected: trio/intermediates/trio-es-greetings.txt
    ```

The collection step ran again, since its batch name input changed, and all the outputs of this batch are now under `results/trio/`.
As long as you specify a different batch name, subsequent runs on other batches of inputs won't clobber previous results.

### Takeaway

You know how to collect outputs from a batch of process calls and feed them into a joint analysis or summation step, and how to use a parameter to organize the published outputs.

To recap, this is what you've built so far:

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-collect.svg"
</figure>

### What's next?

Learn how to emit multiple outputs and handle them conveniently.

---

## 3. Add an output to the collector step

So far we've been using processes that only produced one output each.
We were able to access their respective outputs very conveniently using the `<process>.out` syntax, which we used both in the context of passing an output to the next process (e.g. `translateGreeting(sayHello.out, params.language)`) and in the context of the `publish:` section (e.g. `greetings = sayHello.out`).

What happens when a process produces more than one?
How do we handle the multiple outputs?
Can we select and use a specific output?

All excellent questions, and the short answer is yes we can!

Multiple outputs will be packaged into separate channels.
We can either choose to give those output channels names, which makes it easy to refer to them individually later on, or we can refer to them by index.

For demonstration purposes, let's say we want to count the number of greetings that are being collected for a given batch of inputs and report it in a file.

### 3.1. Modify the process to count and output the number of greetings

This will require two key changes to the process definition: we need a way to count the greetings and write a report file, then we need to add that report file to the `output` block of the process.

#### 3.1.1. Count the number of greetings collected

Conveniently, Nextflow lets us add arbitrary code in the `script:` block of the process definition, which comes in really handy for doing things like this.

That means we can use Nextflow's built-in `size()` function to get the number of files in the `input_files` list, and write the result to file with an `echo` command.

In the `collectGreetings` process block, make the following code changes:

=== "After"

    ```groovy title="hello-workflow.nf" linenums="77" hl_lines="2 5"
        script:
        count_greetings = input_files.size()
        """
        cat ${input_files} > '${batch_name}-${language}-greetings.txt'
        echo 'Number of greetings processed: ${count_greetings}' > '${batch_name}-${language}-report.txt'
        """
    ```

=== "Before"

    ```groovy title="hello-workflow.nf" linenums="77"
        script:
        """
        cat ${input_files} > '${batch_name}-${language}-greetings.txt'
        """
    ```

The `count_greetings` variable will be computed at runtime.

??? note "What if there is only one greeting?"

    If you run the workflow on a CSV file that contains a single name, the report claims a much larger number of greetings.
    With a single item in the channel, `input_files` holds one file instead of a list of files, and `size()` returns the size of that file in bytes.

    To make Nextflow always provide a list, declare how many files the input accepts with the `arity` option:

    ```groovy
    input:
    path input_files, arity: '1..*'
    ```

    A value of `'1..*'` means one or more files.

#### 3.1.2. Emit the report file and name outputs

In principle all we need to do is to add the report file to the `output:` block.

However, while we're at it, we're also going to add some `emit:` labels to our output declarations. These will enable us to select the outputs by name instead of having to use positional indices.

In the process block, make the following code change:

=== "After"

    ```groovy title="hello-workflow.nf" linenums="74" hl_lines="2 3"
        output:
        path "${batch_name}-${language}-greetings.txt", emit: outfile
        path "${batch_name}-${language}-report.txt", emit: report
    ```

=== "Before"

    ```groovy title="hello-workflow.nf" linenums="74"
        output:
        path "${batch_name}-${language}-greetings.txt"
    ```

The `emit:` labels are optional, and we could have added a label to only one of the outputs.
But as the saying goes, why not both?

!!! tip

    If you don't name the outputs of a process using `emit:`, you can still access them individually by using their respective (zero-based) index.
    For example, you would use `<process>.out[0]` to get the first output, `<process>.out[1]` to get the second output, and so on.

    We prefer naming outputs because otherwise, it's too easy to grab the wrong index by error, especially when the process produces a lot of outputs.

### 3.2. Update the workflow outputs

Now that we have two outputs coming out of the `collectGreetings` process, the `collectGreetings.out` output contains two channels:

- `collectGreetings.out.outfile` contains the final output file
- `collectGreetings.out.report` contains the report file

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-collect-report.svg"
</figure>

We need to update the workflow outputs accordingly.

#### 3.2.1. Update the `publish:` section

In the `workflow` block, make the following code change:

=== "After"

    ```groovy title="hello-workflow.nf" linenums="109" hl_lines="4 5"
        publish:
        greetings = sayHello.out
        translated = translateGreeting.out
        collected = collectGreetings.out.outfile
        batch_report = collectGreetings.out.report
    ```

=== "Before"

    ```groovy title="hello-workflow.nf" linenums="109" hl_lines="4"
        publish:
        greetings = sayHello.out
        translated = translateGreeting.out
        collected = collectGreetings.out
    ```

As you can see, referring to specific process outputs is now trivial.
When we go to add one more step to our pipeline in Part 4 (Containers), we'll be able to easily refer to `collectGreetings.out.outfile` and hand it to the new process (spoiler: the new process is called `cowpy`).

But for now, let's finish updating the workflow-level outputs.

#### 3.2.2. Update the `output` block

In the `output` block, make the following code change:

=== "After"

    ```groovy title="hello-workflow.nf" linenums="116" hl_lines="14-17"
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

=== "Before"

    ```groovy title="hello-workflow.nf" linenums="116"
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
    }
    ```

We don't need to update the `collected` output definition since that name hasn't changed.
We just need to add the new output, which we publish at the top level of the batch directory, since it is a final result rather than an intermediate file.

### 3.3. Run the workflow

Run the workflow with the same batch name as before:

```bash
nextflow run hello-workflow.nf -resume --batch trio
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `hello-workflow.nf` [loquacious_ride] revision: 375ebf79eb

    executor >  local (1)
    [1f/03f176] sayHello (Celine)                 | 3 of 3, cached: 3 ✔
    [00/c5c01f] translateGreeting (es-Celine.txt) | 3 of 3, cached: 3 ✔
    [c1/d3a120] collectGreetings (trio)           | 1 of 1 ✔

    Outputs:

      /workspaces/training/hello-nextflow/results

      greetings:
        - trio/intermediates/Celine.txt
        - trio/intermediates/Raphael.txt
        - trio/intermediates/Mateo.txt

      translated:
        - trio/intermediates/es-Raphael.txt
        - trio/intermediates/es-Mateo.txt
        - trio/intermediates/es-Celine.txt

      collected: trio/intermediates/trio-es-greetings.txt

      batch_report: trio/trio-es-report.txt
    ```

If you look in the `results/trio/` directory, you'll find the new report file, `trio-es-report.txt`.
Open it to verify that the workflow correctly reported the count of greetings that were processed.

??? abstract "File contents"

    ```txt title="trio-es-report.txt"
    Number of greetings processed: 3
    ```

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-collect-4-way.svg"
</figure>

Feel free to try the workflow on more names, for example with the `data/names-extended.csv` file, which adds two names to the original three.
Use `-resume`, and see which tasks Nextflow runs again.

### Takeaway

You know how to make a process emit multiple named outputs and how to handle them appropriately at the workflow level.

More generally, you understand the key principles involved in connecting processes together in common ways.

### What's next?

Take an extra long break, you've earned it.

When you're ready, move on to [**Part 4: Hello Containers**](./04_hello_containers.md) to learn how to use containers to manage software dependencies more conveniently and reproducibly.

---

## Bonus exercises

These optional exercises let you practice what you learned in this part.
Each one can be solved with what this part covered; open the hint if you get stuck.
The solutions are also in `solutions/3-hello-workflow/bonus/`.

### Exercise 1: Put the translations in their own folder

Starting from `hello-workflow.nf` as it stands at the end of this part, change where the translated greetings are published, so that they go to a `translations` subdirectory of the batch directory instead of `intermediates`.
Run the workflow with `--batch trio -resume` and check `results/trio/`.

??? tip "Hint"

    Where each output is published is set in the `output` block, not in the processes.

??? solution

    Change the `path` of the `translated` output:

    === "After"

        ```groovy title="hello-workflow.nf" linenums="121" hl_lines="2"
            translated {
                path "${params.batch}/translations"
                mode 'copy'
            }
        ```

    === "Before"

        ```groovy title="hello-workflow.nf" linenums="121" hl_lines="2"
            translated {
                path "${params.batch}/intermediates"
                mode 'copy'
            }
        ```

    Run the workflow:

    ```bash
    nextflow run hello-workflow.nf --batch trio -resume
    ```

    ```console title="Output (excerpt)"
    [1f/03f176] sayHello (Celine)                  | 3 of 3, cached: 3 ✔
    [52/df8eff] translateGreeting (es-Raphael.txt) | 3 of 3, cached: 3 ✔
    [c1/d3a120] collectGreetings (trio)            | 1 of 1, cached: 1 ✔

    Outputs:

      /workspaces/training/hello-nextflow/results

      greetings:
        - trio/intermediates/Mateo.txt
        - trio/intermediates/Celine.txt
        - trio/intermediates/Raphael.txt

      translated:
        - trio/translations/es-Mateo.txt
        - trio/translations/es-Celine.txt
        - trio/translations/es-Raphael.txt

      collected: trio/intermediates/trio-es-greetings.txt

      batch_report: trio/trio-es-report.txt
    ```

    Every task was cached, yet the translated files were published to the new location: Nextflow publishes the outputs of cached tasks too.
    The copies published to `results/trio/intermediates/` by earlier runs are still there.

### Exercise 2: Add a language, and predict what reruns

Starting from `hello-workflow.nf` as it stands at the end of this part, add Portuguese to the languages that `translateGreeting` can translate into, using the code `pt` and the greeting `Ola`.
Run the workflow with `--language pt -resume` before and after your change.

Then, before running each of the following commands, predict which steps will run again and which will be cached, and check your predictions:

1. `nextflow run hello-workflow.nf --language pt -resume` (right after adding Portuguese)
2. `nextflow run hello-workflow.nf --language es -resume`
3. `nextflow run hello-workflow.nf --language es --batch spanish -resume`

??? tip "Hint"

    Before your change, the run fails the same way as the `--language de` example in section 1.3.4.
    To predict what reruns, ask for each step whether its inputs or its code changed since a previous run.

??? solution

    Before the change, the run fails with `Unsupported language code: pt`.

    Add an entry to the dictionary:

    ```groovy title="hello-workflow.nf" linenums="48" hl_lines="2"
            'en': 'Hello',       # English
            'pt': 'Ola',         # Portuguese
        }
    ```

    1. **`--language pt`**: `sayHello` is cached; `translateGreeting` and `collectGreetings` run.
       The collected file contains the Portuguese greetings:

        ```console title="Output (excerpt)"
        [03/4bf37b] sayHello (Raphael)                | 3 of 3, cached: 3 ✔
        [63/700a01] translateGreeting (pt-Mateo.txt)  | 3 of 3 ✔
        [fd/bcdd3a] collectGreetings (batch)          | 1 of 1 ✔
        ```

        ```console title="results/batch/intermediates/batch-pt-greetings.txt"
        Ola, Raphael!
        Ola, Mateo!
        Ola, Celine!
        ```

    2. **`--language es`**: `translateGreeting` runs again, even though Spanish translations were computed in earlier runs.
       Adding Portuguese changed the code of `translateGreeting`, and a change to a process's code invalidates all of its cached results.
       `collectGreetings` runs again because its inputs changed: the files are new, and the language is different.
       Since the language is part of the collected file names, the Spanish results are published next to the Portuguese ones instead of replacing them.

        ```console title="Output (excerpt)"
        [1f/03f176] sayHello (Celine)                  | 3 of 3, cached: 3 ✔
        [00/fb475e] translateGreeting (es-Raphael.txt) | 3 of 3 ✔
        [3b/b95314] collectGreetings (batch)           | 1 of 1 ✔
        ```

    3. **`--batch spanish`**: only `collectGreetings` runs, because the batch name is one of its inputs; the other steps are cached.

        ```console title="Output (excerpt)"
        [03/4bf37b] sayHello (Raphael)               | 3 of 3, cached: 3 ✔
        [b4/a7125b] translateGreeting (es-Mateo.txt) | 3 of 3, cached: 3 ✔
        [29/cb8ddc] collectGreetings (spanish)       | 1 of 1 ✔
        ```

### Exercise 3: Collect both versions of the greetings

Starting from `hello-workflow.nf` as it stands at the end of this part, add a step that writes all the original greetings and all the translated greetings into a single file called `<batch>-<language>-bilingual.txt`, published at the top level of the batch directory.
Run it with `--batch trio -resume`.

??? tip "Hint"

    The new process needs two collections of files as inputs, one from `sayHello` and one from `translateGreeting`, plus the batch name and the language for the output file name, as in `collectGreetings`.
    Its script only needs a `cat` command.

??? solution

    Define the new process after `collectGreetings`:

    ```groovy title="hello-workflow.nf" linenums="86"
    /*
     * Combine the original and translated greetings into a single file
     */
    process collectBilingual {

        tag "${batch_name}"

        input:
        path greetings
        path translations
        val batch_name
        val language

        output:
        path "${batch_name}-${language}-bilingual.txt"

        script:
        """
        cat ${greetings} ${translations} > '${batch_name}-${language}-bilingual.txt'
        """
    }
    ```

    Call it in the workflow, collecting the outputs of both upstream processes, and publish its output:

    ```groovy title="hello-workflow.nf" linenums="128" hl_lines="3-4 11"
        // collect all the greetings into one file
        collectGreetings(translateGreeting.out.collect(), params.batch, params.language)
        // combine the original and translated greetings
        collectBilingual(sayHello.out.collect(), translateGreeting.out.collect(), params.batch, params.language)

        publish:
        greetings = sayHello.out
        translated = translateGreeting.out
        collected = collectGreetings.out.outfile
        batch_report = collectGreetings.out.report
        bilingual = collectBilingual.out
    }
    ```

    Add the matching entry to the `output` block:

    ```groovy title="hello-workflow.nf" linenums="154" hl_lines="5-8"
        batch_report {
            path "${params.batch}"
            mode 'copy'
        }
        bilingual {
            path "${params.batch}"
            mode 'copy'
        }
    }
    ```

    Run the workflow:

    ```bash
    nextflow run hello-workflow.nf --batch trio -resume
    ```

    ```console title="Output (excerpt)"
    [1f/03f176] sayHello (Celine)                  | 3 of 3, cached: 3 ✔
    [52/df8eff] translateGreeting (es-Raphael.txt) | 3 of 3, cached: 3 ✔
    [c1/d3a120] collectGreetings (trio)            | 1 of 1, cached: 1 ✔
    [7f/557284] collectBilingual (trio)            | 1 of 1 ✔
    ```

    Only the new step runs; everything upstream is cached.

    ```console title="results/trio/trio-es-bilingual.txt"
    Hello, Raphael!
    Hello, Mateo!
    Hello, Celine!
    Hola, Raphael!
    Hola, Celine!
    Hola, Mateo!
    ```

---

## Quiz

<quiz>
How do you access the output of a process in the workflow block?
- [x] `processName.out`
- [ ] `output.processName`
- [ ] `process.output`
- [ ] `get(processName)`

Learn more: [1.2.2. Pass the output of the first process to the new process](#122-pass-the-output-of-the-first-process-to-the-new-process)
</quiz>

<quiz>
Why does the task directory of a `translateGreeting` call also contain the greeting file produced by `sayHello`?
- [ ] The `output` block copied it there to publish it
- [ ] The `translateGreeting` process ran `sayHello` again to recreate it
- [ ] Both processes ran in the same task directory
- [x] Nextflow staged it there as an input, as a symbolic link to the original file

Learn more: [1.2.4. Run the workflow with `-resume`](#124-run-the-workflow-with-resume)
</quiz>

<quiz>
What does the line `#!/usr/bin/env python3` at the start of a process script do?
- [ ] It imports the Python standard library
- [x] It tells the system to run the script with Python instead of the shell
- [ ] It is a comment that Nextflow ignores
- [ ] It installs Python in the task environment

Learn more: [1.2.1. Define the process](#121-define-the-process)
</quiz>

<quiz>
A process is called with `translateGreeting(sayHello.out, params.language)`, where `sayHello.out` carries three files and `params.language` is a single value. How many times does the process run?
- [x] Three times, reusing the same language each time
- [ ] It fails because the inputs have different sizes
- [ ] Once per language
- [ ] Once, on all three files

Learn more: [1.3.3. Pass the language to the process](#133-pass-the-language-to-the-process)
</quiz>

<quiz>
When providing multiple inputs to a process, what must be true?
- [ ] Only two inputs can be provided at a time
- [ ] All inputs must be of the same type
- [ ] Inputs must be provided in alphabetical order
- [x] The order of inputs must match the order defined in the input block

Learn more: [1.3.3. Pass the language to the process](#133-pass-the-language-to-the-process)
</quiz>

<quiz>
Which operator should replace `???` to gather all outputs into a single list for the downstream process?

```groovy hl_lines="4"
workflow {
    names_ch = channel.of('Mateo', 'Raphael', 'Celine')
    sayHello(names_ch)
    gatherAll(sayHello.out.???)
}
```

- [ ] `flatten()`
- [ ] `join()`
- [x] `collect()`
- [ ] `mix()`

Learn more: [2.4. Use an operator to collect the greetings into a single input](#24-use-an-operator-to-collect-the-greetings-into-a-single-input)
</quiz>

<quiz>
When should you use the `collect()` operator?
- [ ] When you need to filter channel contents
- [x] When a downstream process needs all items from an upstream process
- [ ] When you want to split data across multiple processes
- [ ] When you want to process items in parallel

Learn more: [2.4. Use an operator to collect the greetings into a single input](#24-use-an-operator-to-collect-the-greetings-into-a-single-input)
</quiz>

<quiz>
How do you access a named output from a process?
- [ ] `output.processName.outputName`
- [ ] `processName.get(outputName)`
- [x] `processName.out.outputName`
- [ ] `processName.outputName`

Learn more: [3.1.2. Emit the report file and name outputs](#312-emit-the-report-file-and-name-outputs)
</quiz>

<quiz>
What is the correct syntax for naming an output in a process?
- [ ] `name: outputName`
- [x] `emit: outputName`
- [ ] `label: outputName`
- [ ] `output: outputName`

Learn more: [3.1.2. Emit the report file and name outputs](#312-emit-the-report-file-and-name-outputs)
</quiz>
