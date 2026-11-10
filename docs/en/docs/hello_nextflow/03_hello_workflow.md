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
2. Add a third step that collects all the translated greetings and writes them into a single file named after the batch.
3. Make the collection step also report a simple statistic about what was processed.

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

    Launching `hello-workflow.nf` [curious_albattani] revision: cef606f370

    executor >  local (3)
    [b0/9e0333] sayHello (1) | 3 of 3 ✔

    Outputs:

      /workspaces/training/hello-nextflow/results

      greetings:
        - hello/Raphael.txt
        - hello/Celine.txt
        - hello/Mateo.txt
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

```groovy title="hello-workflow.nf" linenums="20"
/*
 * Use Python to translate the greeting into Spanish
 */
process translateGreeting {

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

    ```groovy title="hello-workflow.nf" linenums="50" hl_lines="10-11"
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

    ```groovy title="hello-workflow.nf" linenums="50"
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

    ```groovy title="hello-workflow.nf" linenums="62" hl_lines="3"
        publish:
        greetings = sayHello.out
        translated = translateGreeting.out
    }
    ```

=== "Before"

    ```groovy title="hello-workflow.nf" linenums="62"
        publish:
        greetings = sayHello.out
    }
    ```

Then, in the `output` block, make the following code change:

=== "After"

    ```groovy title="hello-workflow.nf" linenums="67" hl_lines="6-9"
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

    ```groovy title="hello-workflow.nf" linenums="67"
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

Run the workflow with the `-resume` flag, which you used in Nextflow Run, since the first step of the workflow already ran successfully:

```bash
nextflow run hello-workflow.nf -resume
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `hello-workflow.nf` [suspicious_agnesi] revision: fa84af21f4

    executor >  local (3)
    [0d/1a60fc] sayHello (2)          | 3 of 3, cached: 3 ✔
    [c7/6a78e1] translateGreeting (3) | 3 of 3 ✔

    Outputs:

      /workspaces/training/hello-nextflow/results

      greetings:
        - hello/Celine.txt
        - hello/Mateo.txt
        - hello/Raphael.txt

      translated:
        - hello/es-Mateo.txt
        - hello/es-Raphael.txt
        - hello/es-Celine.txt
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
    work/45/3f7825238437f4f4467cd6a84ef195/
    ├── Mateo.txt -> /workspaces/training/hello-nextflow/work/b0/9e0333111ec99af27dede6f689dc5b/Mateo.txt
    └── es-Mateo.txt
    ```

Notice there are two output files: the output of the first process as well as the output of the second.

The output of the first process is in there because Nextflow **staged** it there in order to have everything needed for execution within the same subdirectory.
It is a symbolic link pointing to the original file in the subdirectory of the first process call.
By default, when running on a single machine as we're doing here, Nextflow uses symbolic links rather than copies to stage input and intermediate files.

The `.command.sh` file in the same directory shows the Python script that Nextflow ran, with the file names filled in:

```python title="work/45/3f7825238437f4f4467cd6a84ef195/.command.sh"
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

    ```groovy title="hello-workflow.nf" linenums="20" hl_lines="2 8 11 16 18-26 28-29 34-35"
    /*
     * Use Python to translate the greeting into another language
     */
    process translateGreeting {

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

    ```groovy title="hello-workflow.nf" linenums="20"
    /*
     * Use Python to translate the greeting into Spanish
     */
    process translateGreeting {

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

#### 1.3.2. Add a `language` parameter

You already know how to use the `params` system to declare CLI parameters.
Declare a `language` parameter, with Spanish (`es`) as the default value.

In the pipeline parameters section, make the following code change:

=== "After"

    ```groovy title="hello-workflow.nf" linenums="58" hl_lines="6"
    /*
     * Pipeline parameters
     */
    params {
        input: Path = 'data/names.csv'
        language: String = 'es'
    }
    ```

=== "Before"

    ```groovy title="hello-workflow.nf" linenums="58"
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

    ```groovy title="hello-workflow.nf" linenums="75" hl_lines="2"
        // translate the greeting
        translateGreeting(sayHello.out, params.language)
    ```

=== "Before"

    ```groovy title="hello-workflow.nf" linenums="75"
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

    Launching `hello-workflow.nf` [irreverent_dubinsky] revision: ebdc8255be

    executor >  local (3)
    [55/f9a2a6] sayHello (3)          | 3 of 3, cached: 3 ✔
    [44/9c4185] translateGreeting (1) | 3 of 3 ✔

    Outputs:

      /workspaces/training/hello-nextflow/results

      greetings:
        - hello/Celine.txt
        - hello/Mateo.txt
        - hello/Raphael.txt

      translated:
        - hello/es-Mateo.txt
        - hello/es-Raphael.txt
        - hello/es-Celine.txt
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

    Launching `hello-workflow.nf` [gloomy_maxwell] revision: ebdc8255be

    executor >  local (3)
    [b0/9e0333] sayHello (1)          | 3 of 3, cached: 3 ✔
    [e1/d0a3dd] translateGreeting (2) | 3 of 3 ✔

    Outputs:

      /workspaces/training/hello-nextflow/results

      greetings:
        - hello/Raphael.txt
        - hello/Mateo.txt
        - hello/Celine.txt

      translated:
        - hello/fr-Raphael.txt
        - hello/fr-Celine.txt
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

    Launching `hello-workflow.nf` [gigantic_rosalind] revision: 89eca833a5

    executor >  local (3)
    [cc/023c93] sayHello (1)          | 3 of 3, cached: 3 ✔
    [09/342f5a] translateGreeting (3) | 0 of 3 ✘
    ERROR ~ Error executing process > 'translateGreeting (3)'

    Caused by:
      Process `translateGreeting (3)` terminated with an error exit status (1)


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
      /workspaces/training/hello-nextflow/work/09/342f5abdee49ba5388f9cb6013731e

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
We'll also name that file after a batch name provided as a parameter, so that the results of different batches of inputs don't overwrite each other.

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

```groovy title="hello-workflow.nf" linenums="58"
/*
 * Collect translated greetings into a single output file
 */
process collectGreetings {

    input:
    path input_files
    val batch_name

    output:
    path "${batch_name}-greetings.txt"

    script:
    """
    cat ${input_files} > '${batch_name}-greetings.txt'
    """
}
```

This process has two inputs, like `translateGreeting`:

- `input_files` will hold the files to collect.
  It uses the `path` qualifier even though we expect it to contain multiple files.
- `batch_name` is a plain value that we use in the output file name, so that each batch gets its own output file.

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

    ```groovy title="hello-workflow.nf" linenums="76" hl_lines="7"
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

    ```groovy title="hello-workflow.nf" linenums="76"
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

    ```groovy title="hello-workflow.nf" linenums="94" hl_lines="3-4"
        // translate the greeting
        translateGreeting(sayHello.out, params.language)
        // collect all the greetings into one file
        collectGreetings(translateGreeting.out, params.batch)
    ```

=== "Before"

    ```groovy title="hello-workflow.nf" linenums="94"
        // translate the greeting
        translateGreeting(sayHello.out, params.language)
    ```

This connects the output of `translateGreeting()` to the first input of `collectGreetings()`, and passes the batch name as the second input, just like we did with the language.

#### 2.3.3. Update the workflow outputs

In the `publish:` section of the `workflow` block, make the following code change:

=== "After"

    ```groovy title="hello-workflow.nf" linenums="99" hl_lines="4"
        publish:
        greetings = sayHello.out
        translated = translateGreeting.out
        collected = collectGreetings.out
    }
    ```

=== "Before"

    ```groovy title="hello-workflow.nf" linenums="99"
        publish:
        greetings = sayHello.out
        translated = translateGreeting.out
    }
    ```

Then, in the `output` block, make the following code change:

=== "After"

    ```groovy title="hello-workflow.nf" linenums="105" hl_lines="10-13"
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

    ```groovy title="hello-workflow.nf" linenums="105"
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

    Launching `hello-workflow.nf` [backstabbing_lumiere] revision: d2c08a3285

    executor >  local (3)
    [55/f9a2a6] sayHello (3)          | 3 of 3, cached: 3 ✔
    [44/9c4185] translateGreeting (2) | 3 of 3, cached: 3 ✔
    [d9/315433] collectGreetings (3)  | 3 of 3 ✔

    Outputs:

      /workspaces/training/hello-nextflow/results

      greetings:
        - hello/Raphael.txt
        - hello/Celine.txt
        - hello/Mateo.txt

      translated:
        - hello/es-Mateo.txt
        - hello/es-Celine.txt
        - hello/es-Raphael.txt

      collected:
        - hello/batch-greetings.txt
        - hello/batch-greetings.txt
        - hello/batch-greetings.txt
    ```

It runs successfully, including the third step.

However, look at the number of calls for `collectGreetings()` on the last line.
We were only expecting one, but there are three.

Now have a look at the contents of the final output file.

??? abstract "File contents"

    ```console title="results/hello/batch-greetings.txt"
    Hola, Raphael!
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

    ```groovy title="hello-workflow.nf" linenums="96" hl_lines="2"
        // collect all the greetings into one file
        collectGreetings(translateGreeting.out.collect(), params.batch)
    ```

=== "Before"

    ```groovy title="hello-workflow.nf" linenums="96"
        // collect all the greetings into one file
        collectGreetings(translateGreeting.out, params.batch)
    ```

#### 2.4.2. Add some `view()` statements

Let's also include a couple of `view()` statements to visualize the before and after states of the channel contents.

=== "After"

    ```groovy title="hello-workflow.nf" linenums="96" hl_lines="4-6"
        // collect all the greetings into one file
        collectGreetings(translateGreeting.out.collect(), params.batch)

        // optional view statements
        translateGreeting.out.view { filepath -> "Before collect: $filepath" }
        translateGreeting.out.collect().view { filepaths -> "After collect: $filepaths" }
    ```

=== "Before"

    ```groovy title="hello-workflow.nf" linenums="96"
        // collect all the greetings into one file
        collectGreetings(translateGreeting.out.collect(), params.batch)
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

    Launching `hello-workflow.nf` [serene_panini] revision: 116be6008d

    executor >  local (1)
    [b0/9e0333] sayHello (1)          | 3 of 3, cached: 3 ✔
    [ae/3c0ca0] translateGreeting (2) | 3 of 3, cached: 3 ✔
    [db/6e3931] collectGreetings      | 1 of 1 ✔
    Before collect: /workspaces/training/hello-nextflow/work/44/9c4185f2ac241b19cb56a67526a966/es-Celine.txt
    Before collect: /workspaces/training/hello-nextflow/work/8c/174fc0d15b4f859d0914e3ae0e5fbd/es-Raphael.txt
    Before collect: /workspaces/training/hello-nextflow/work/ae/3c0ca00c751ea9330d7f103d4eae25/es-Mateo.txt
    After collect: [/workspaces/training/hello-nextflow/work/44/9c4185f2ac241b19cb56a67526a966/es-Celine.txt, /workspaces/training/hello-nextflow/work/8c/174fc0d15b4f859d0914e3ae0e5fbd/es-Raphael.txt, /workspaces/training/hello-nextflow/work/ae/3c0ca00c751ea9330d7f103d4eae25/es-Mateo.txt]

    Outputs:

      /workspaces/training/hello-nextflow/results

      greetings:
        - hello/Raphael.txt
        - hello/Mateo.txt
        - hello/Celine.txt

      translated:
        - hello/es-Celine.txt
        - hello/es-Raphael.txt
        - hello/es-Mateo.txt

      collected: hello/batch-greetings.txt
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

    ```console title="results/hello/batch-greetings.txt"
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

    ```groovy title="hello-workflow.nf" linenums="96"
        // collect all the greetings into one file
        collectGreetings(translateGreeting.out.collect(), params.batch)
    ```

=== "Before"

    ```groovy title="hello-workflow.nf" linenums="96" hl_lines="4-6"
        // collect all the greetings into one file
        collectGreetings(translateGreeting.out.collect(), params.batch)

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

    ```groovy title="hello-workflow.nf" linenums="105" hl_lines="3 7 11"
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

    ```groovy title="hello-workflow.nf" linenums="105"
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
As you saw in Nextflow Run, the person running the pipeline can change the location of the output directory, but not the organization within it, except through parameters like `batch` that the author builds into the paths.

#### 2.5.2. Run the workflow with a batch name

Run the workflow with a batch name on the command line:

```bash
nextflow run hello-workflow.nf -resume --batch trio
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `hello-workflow.nf` [nostalgic_noyce] revision: 45b5328138

    executor >  local (1)
    [55/f9a2a6] sayHello (3)          | 3 of 3, cached: 3 ✔
    [ae/3c0ca0] translateGreeting (1) | 3 of 3, cached: 3 ✔
    [93/8804a2] collectGreetings      | 1 of 1 ✔

    Outputs:

      /workspaces/training/hello-nextflow/results

      greetings:
        - trio/intermediates/Mateo.txt
        - trio/intermediates/Celine.txt
        - trio/intermediates/Raphael.txt

      translated:
        - trio/intermediates/es-Celine.txt
        - trio/intermediates/es-Raphael.txt
        - trio/intermediates/es-Mateo.txt

      collected: trio/intermediates/trio-greetings.txt
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

    ```groovy title="hello-workflow.nf" linenums="70" hl_lines="2 5"
        script:
        count_greetings = input_files.size()
        """
        cat ${input_files} > '${batch_name}-greetings.txt'
        echo 'Number of greetings processed: ${count_greetings}' > '${batch_name}-report.txt'
        """
    ```

=== "Before"

    ```groovy title="hello-workflow.nf" linenums="70"
        script:
        """
        cat ${input_files} > '${batch_name}-greetings.txt'
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

However, while we're at it, we're also going to add some `emit:` tags to our output declarations. These will enable us to select the outputs by name instead of having to use positional indices.

In the process block, make the following code change:

=== "After"

    ```groovy title="hello-workflow.nf" linenums="67" hl_lines="2 3"
        output:
        path "${batch_name}-greetings.txt", emit: outfile
        path "${batch_name}-report.txt", emit: report
    ```

=== "Before"

    ```groovy title="hello-workflow.nf" linenums="67"
        output:
        path "${batch_name}-greetings.txt"
    ```

The `emit:` tags are optional, and we could have added a tag to only one of the outputs.
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

    ```groovy title="hello-workflow.nf" linenums="102" hl_lines="4 5"
        publish:
        greetings = sayHello.out
        translated = translateGreeting.out
        collected = collectGreetings.out.outfile
        batch_report = collectGreetings.out.report
    ```

=== "Before"

    ```groovy title="hello-workflow.nf" linenums="102" hl_lines="4"
        publish:
        greetings = sayHello.out
        translated = translateGreeting.out
        collected = collectGreetings.out
    ```

As you can see, referring to specific process outputs is now trivial.
When we go to add one more step to our pipeline in Part 5 (Containers), we'll be able to easily refer to `collectGreetings.out.outfile` and hand it to the new process (spoiler: the new process is called `cowpy`).

But for now, let's finish updating the workflow-level outputs.

#### 3.2.2. Update the `output` block

In the `output` block, make the following code change:

=== "After"

    ```groovy title="hello-workflow.nf" linenums="109" hl_lines="14-17"
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

    ```groovy title="hello-workflow.nf" linenums="109"
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

    Launching `hello-workflow.nf` [chaotic_miescher] revision: b4393cf728

    executor >  local (1)
    [55/f9a2a6] sayHello (3)          | 3 of 3, cached: 3 ✔
    [44/9c4185] translateGreeting (3) | 3 of 3, cached: 3 ✔
    [06/567dd8] collectGreetings      | 1 of 1 ✔

    Outputs:

      /workspaces/training/hello-nextflow/results

      greetings:
        - trio/intermediates/Raphael.txt
        - trio/intermediates/Mateo.txt
        - trio/intermediates/Celine.txt

      translated:
        - trio/intermediates/es-Mateo.txt
        - trio/intermediates/es-Celine.txt
        - trio/intermediates/es-Raphael.txt

      collected: trio/intermediates/trio-greetings.txt

      batch_report: trio/trio-report.txt
    ```

If you look in the `results/trio/` directory, you'll find the new report file, `trio-report.txt`.
Open it to verify that the workflow correctly reported the count of greetings that were processed.

??? abstract "File contents"

    ```txt title="trio-report.txt"
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

When you're ready, move on to [**Part 4: Hello Modules**](./04_hello_modules.md) to learn how to modularize your code for better maintainability and code efficiency.

---

## Bonus exercises

These optional exercises let you practice what you learned in this part.
Each one can be solved with what this part covered; open the hint if you get stuck.
The solutions are also in `solutions/3-hello-workflow/bonus/`.

### Exercise 1: Sort the collected greetings

As the note in section 2.4.3 points out, the order of the greetings in the collected file can change from one run to the next.
Starting from `hello-workflow.nf` as it stands at the end of this part, modify `collectGreetings` so that the greetings in the collected file are always in the same order.
Run it with `--batch sorted` and check `results/sorted/intermediates/sorted-greetings.txt`.

??? tip "Hint"

    The script block can hold a shell pipeline: the `sort` command sorts its input lines alphabetically.

??? solution

    Pipe the concatenated files through `sort` before writing the collected file:

    ```groovy title="hello-workflow.nf" linenums="71" hl_lines="4"
        script:
        count_greetings = input_files.size()
        """
        cat ${input_files} | sort > '${batch_name}-greetings.txt'
        echo 'Number of greetings processed: ${count_greetings}' > '${batch_name}-report.txt'
        """
    ```

    ```bash
    nextflow run hello-workflow.nf --batch sorted
    ```

    The greetings are now in alphabetical order, whatever order the tasks finished in:

    ```console title="results/sorted/intermediates/sorted-greetings.txt"
    Hola, Celine!
    Hola, Mateo!
    Hola, Raphael!
    ```

### Exercise 2: Add a language

Starting from `hello-workflow.nf` as it stands at the end of this part, add Portuguese to the languages that `translateGreeting` can translate into.
Run the workflow with `--language pt --batch portuguese` and check the collected file.

??? tip "Hint"

    Only the `translations` dictionary in the Python script needs to change; use a romanized greeting such as `Ola`.

??? solution

    Add an entry to the dictionary:

    ```groovy title="hello-workflow.nf" linenums="44" hl_lines="2"
            'en': 'Hello',       # English
            'pt': 'Ola',         # Portuguese
        }
    ```

    ```bash
    nextflow run hello-workflow.nf --language pt --batch portuguese
    ```

    ```console title="results/portuguese/intermediates/portuguese-greetings.txt"
    Ola, Celine!
    Ola, Raphael!
    Ola, Mateo!
    ```

### Challenge: Translate each name into its own language

The second column of `data/names.csv` gives a language for each name.
Starting from `hello-workflow.nf` as it stands at the end of this part, translate each greeting into the language of its own row, instead of one language for the whole batch.
Run it with `--batch mixed` and check that the collected file contains `Hola, Mateo!`, `Ciao, Raphael!` and `Bonjour, Celine!`.

??? tip "Hint"

    Each item now needs to carry two values, the name and its language, through `sayHello` to `translateGreeting`.
    A channel element can be a list of values, and a process can accept and emit it with a `tuple`:

    ```groovy
    .map { row -> [row[0], row[1]] }       // channel elements: [name, language]

    input:
    tuple val(name), val(language)          // in sayHello

    output:
    tuple path("${name}.txt"), val(language)
    ```

    `translateGreeting` then takes `tuple path(input_file), val(language)` as its only input.

??? solution

    Keep both columns in the channel:

    ```groovy title="hello-workflow.nf" linenums="92" hl_lines="3"
        name_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { row -> [row[0], row[1]] }
    ```

    Make `sayHello` accept the pair and pass the language along with its output:

    ```groovy title="hello-workflow.nf" linenums="6" hl_lines="4 7"
    process sayHello {

        input:
        tuple val(name), val(language)

        output:
        tuple path("${name}.txt"), val(language)
    ```

    Make `translateGreeting` take that pair as a single input:

    ```groovy title="hello-workflow.nf" linenums="25" hl_lines="2"
        input:
        tuple path(input_file), val(language)
    ```

    Then call it with the output of `sayHello` alone, and remove the `language` parameter, which is no longer used:

    ```groovy title="hello-workflow.nf" linenums="95" hl_lines="2"
        // translate the greeting
        translateGreeting(sayHello.out)
    ```

    ```bash
    nextflow run hello-workflow.nf --batch mixed
    ```

    ```console title="results/mixed/intermediates/mixed-greetings.txt"
    Hola, Mateo!
    Ciao, Raphael!
    Bonjour, Celine!
    ```

    The `greetings` entries in the console output now show each file next to its language, since `sayHello` emits pairs; the files themselves are published as before.

---

## Quiz

<quiz>
How do you access the output of a process in the workflow block?
- [ ] `process.output`
- [ ] `output.processName`
- [x] `processName.out`
- [ ] `get(processName)`

Learn more: [1.2.2. Pass the output of the first process to the new process](#122-pass-the-output-of-the-first-process-to-the-new-process)
</quiz>

<quiz>
What determines the order of process execution in Nextflow?
- [ ] The order processes are written in the workflow block
- [ ] Alphabetical order by process name
- [x] Data dependencies between processes
- [ ] Random order for parallel execution

Learn more: [1.2.2. Pass the output of the first process to the new process](#122-pass-the-output-of-the-first-process-to-the-new-process)
</quiz>

<quiz>
What does the line `#!/usr/bin/env python3` at the start of a process script do?
- [ ] It is a comment that Nextflow ignores
- [x] It tells the system to run the script with Python instead of the shell
- [ ] It installs Python in the task environment
- [ ] It imports the Python standard library

Learn more: [1.2.1. Define the process](#121-define-the-process)
</quiz>

<quiz>
A process is called with `translateGreeting(sayHello.out, params.language)`, where `sayHello.out` carries three files and `params.language` is a single value. How many times does the process run?
- [ ] Once, on all three files
- [x] Three times, reusing the same language each time
- [ ] Once per language
- [ ] It fails because the inputs have different sizes

Learn more: [1.3.3. Pass the language to the process](#133-pass-the-language-to-the-process)
</quiz>

<quiz>
When providing multiple inputs to a process, what must be true?
- [ ] All inputs must be of the same type
- [ ] Inputs must be provided in alphabetical order
- [x] The order of inputs must match the order defined in the input block
- [ ] Only two inputs can be provided at a time

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
- [x] `collect()`
- [ ] `mix()`
- [ ] `join()`

Learn more: [2.4. Use an operator to collect the greetings into a single input](#24-use-an-operator-to-collect-the-greetings-into-a-single-input)
</quiz>

<quiz>
When should you use the `collect()` operator?
- [ ] When you want to process items in parallel
- [ ] When you need to filter channel contents
- [x] When a downstream process needs all items from an upstream process
- [ ] When you want to split data across multiple processes

Learn more: [2.4. Use an operator to collect the greetings into a single input](#24-use-an-operator-to-collect-the-greetings-into-a-single-input)
</quiz>

<quiz>
How do you access a named output from a process?
- [ ] `processName.outputName`
- [ ] `processName.get(outputName)`
- [x] `processName.out.outputName`
- [ ] `output.processName.outputName`

Learn more: [3.1.2. Emit the report file and name outputs](#312-emit-the-report-file-and-name-outputs)
</quiz>

<quiz>
What is the correct syntax for naming an output in a process?
- [ ] `name: outputName`
- [ ] `output: outputName`
- [x] `emit: outputName`
- [ ] `label: outputName`

Learn more: [3.1.2. Emit the report file and name outputs](#312-emit-the-report-file-and-name-outputs)
</quiz>
