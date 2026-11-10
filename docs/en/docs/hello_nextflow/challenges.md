# Challenges

The bonus exercises at the end of each part can be solved with what that part covers.
The challenges on this page go further on purpose: each one asks you to extend the Hello pipeline with a Nextflow feature that the course does not cover.

Each challenge starts from the end state of one part of the course.
Copy the solution file for that state to `hello-challenge.nf` in the `hello-nextflow` directory, as shown in each challenge, and make your changes in that file.
Your own course scripts stay as they are.

The hint for each challenge introduces the new feature and links to the Nextflow documentation.
Each solution ends with a pointer to the side quest that covers the feature in depth.

The solution files are in `solutions/challenges/`, one directory per challenge.
To run a solution directly, launch it from the `hello-nextflow` directory:

```bash
nextflow run solutions/challenges/challenge-2/hello-challenge.nf
```

---

## Challenge 1: Read a CSV file with a header row

Start from `hello-channels.nf` as it stands at the end of Part 2:

```bash
cp solutions/2-hello-channels/hello-channels-4.nf hello-challenge.nf
```

Many CSV files start with a header row that names the columns.
Create a copy of `data/names.csv` called `data/names-header.csv`, and add the header row `name,language,character` at the top.
Then adapt `hello-challenge.nf` so that it reads this file and picks the name by its column name instead of its position.
The workflow should still greet the same three names.

??? tip "Hint"

    The [`splitCsv()`](https://nextflow.io/docs/latest/reference/operator.html#splitcsv) operator has a `header` option.
    With `header: true`, Nextflow reads the column names from the first row, and returns each following row as a set of named values.
    You can then pick a value by its column name, as in `row.name`, instead of by its position, as in `row[0]`.

??? solution

    ```csv title="data/names-header.csv"
    name,language,character
    Mateo,es,squirrel
    Raphael,it,turtle
    Celine,fr,moose
    ```

    Point the `input` parameter to the new file:

    ```groovy title="hello-challenge.nf" linenums="23" hl_lines="2"
    params {
        input: Path = 'data/names-header.csv'
    }
    ```

    Then parse the header and pick the column by name:

    ```groovy title="hello-challenge.nf" linenums="31" hl_lines="3 5"
        name_ch = channel.fromPath(params.input)
                            .view { csv -> "Before splitCsv: $csv" }
                            .splitCsv(header: true)
                            .view { row -> "After splitCsv: $row" }
                            .map { row -> row.name }
                            .view { name -> "After map: $name" }
    ```

    Run the workflow:

    ```bash
    nextflow run hello-challenge.nf
    ```

    ```console title="Output (excerpt)"
    After splitCsv: [name:Mateo, language:es, character:squirrel]
    After splitCsv: [name:Raphael, language:it, character:turtle]
    After splitCsv: [name:Celine, language:fr, character:moose]
    After map: Mateo
    After map: Raphael
    After map: Celine
    ```

    Each row is now a set of named values, and the header row itself is not passed on as a row.

    The solution directory also contains `names-header.csv`.
    To run the solution directly, point the `input` parameter to that copy:

    ```bash
    nextflow run solutions/challenges/challenge-1/hello-challenge.nf --input solutions/challenges/challenge-1/names-header.csv
    ```

    To learn more about reading CSV files with a header row, see the [Splitting and Grouping](../side_quests/splitting_and_grouping/index.md#11-read-in-sample-data-with-splitcsv-and-create-meta-maps) side quest.

---

## Challenge 2: Greet only the French names

Start from `hello-channels.nf` as it stands at the end of Part 2:

```bash
cp solutions/2-hello-channels/hello-channels-4.nf hello-challenge.nf
```

Change the workflow so that it only greets the names whose language, in the second column of `data/names.csv`, is `fr`.
Only one greeting file should be produced.

??? tip "Hint"

    Use the [`filter()`](https://nextflow.io/docs/latest/reference/operator.html#filter) operator.
    Like `map()`, it takes a closure, and keeps only the elements for which the closure returns `true`.

??? solution

    Add a `filter()` step between `splitCsv()` and `map()`, while each row still holds all its columns:

    ```groovy title="hello-challenge.nf" linenums="31" hl_lines="5"
        name_ch = channel.fromPath(params.input)
                            .view { csv -> "Before splitCsv: $csv" }
                            .splitCsv()
                            .view { row -> "After splitCsv: $row" }
                            .filter { row -> row[1] == 'fr' }
                            .map { row -> row[0] }
                            .view { name -> "After map: $name" }
    ```

    Run the workflow:

    ```bash
    nextflow run hello-challenge.nf
    ```

    ```console title="Output (excerpt)"
    [3f/9c1f02] sayHello (1) | 1 of 1 ✔
    Before splitCsv: /workspaces/training/hello-nextflow/data/names.csv
    After splitCsv: [Mateo, es, squirrel]
    After splitCsv: [Raphael, it, turtle]
    After splitCsv: [Celine, fr, moose]
    After map: Celine
    ```

    All three rows are read, but only Celine's row passes the filter, so `sayHello` runs once and writes `Hello, Celine!`.

    To learn more about filtering channel contents, see the [Splitting and Grouping](../side_quests/splitting_and_grouping/index.md#21-filter-data-with-filter) side quest.

---

## Challenge 3: Translate each name into its own language

Start from `hello-workflow.nf` as it stands at the end of Part 3:

```bash
cp solutions/3-hello-workflow/hello-workflow-4.nf hello-challenge.nf
```

The second column of `data/names.csv` gives a language for each name.
Translate each greeting into the language of its own row, instead of one language for the whole batch.
Run it with `--batch mixed` and check that the collected file contains `Hola, Mateo!`, `Ciao, Raphael!` and `Bonjour, Celine!`.

??? tip "Hint"

    Each item now needs to carry two values, the name and its language, through `sayHello` to `translateGreeting`.
    A channel element can be a list of values, and a process can accept and emit such an element as a single [tuple](https://nextflow.io/docs/latest/process.html#input-tuples-tuple):

    ```groovy
    .map { row -> [row[0], row[1]] }       // channel elements: [name, language]

    input:
    tuple val(name), val(language)          // in sayHello

    output:
    tuple path("${name}.txt"), val(language)
    ```

    `translateGreeting` then takes `tuple path(input_file), val(language)` as its only input.

??? solution

    Make `sayHello` accept the pair and pass the language along with its output:

    ```groovy title="hello-challenge.nf" linenums="6" hl_lines="4 7"
    process sayHello {

        input:
        tuple val(name), val(language)

        output:
        tuple path("${name}.txt"), val(language)
    ```

    Make `translateGreeting` take that pair as a single input:

    ```groovy title="hello-challenge.nf" linenums="25" hl_lines="2"
        input:
        tuple path(input_file), val(language)
    ```

    In the workflow, keep both columns in the channel, and call `translateGreeting` with the output of `sayHello` alone:

    ```groovy title="hello-challenge.nf" linenums="89" hl_lines="1-2 4 6-8"
        // create a channel of [name, language] pairs from a CSV file
        person_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { row -> [row[0], row[1]] }
        // emit a greeting
        sayHello(person_ch)
        // translate each greeting into its own language
        translateGreeting(sayHello.out)
    ```

    The `language` parameter is no longer used, so remove it from the `params` block.

    `sayHello` now emits pairs, so publish only the file from each pair:

    ```groovy title="hello-challenge.nf" linenums="100" hl_lines="2"
        publish:
        greetings = sayHello.out.map { greeting_file, _language -> greeting_file }
    ```

    In the `map` closure, the underscore at the start of `_language` tells Nextflow that the variable is deliberately unused, so `nextflow lint` doesn't warn about it.

    Run the workflow:

    ```bash
    nextflow run hello-challenge.nf --batch mixed
    ```

    ??? success "Command output"

        ```console
        N E X T F L O W   ~  version 26.04.4

        Launching `hello-challenge.nf` [shrivelled_leibniz] revision: 6e10bef93d

        executor >  local (7)
        [a8/428799] sayHello (1)          | 3 of 3 ✔
        [95/2e24c9] translateGreeting (2) | 3 of 3 ✔
        [10/f0e0da] collectGreetings      | 1 of 1 ✔

        Outputs:

          /workspaces/training/hello-nextflow/results

          greetings:
            - mixed/intermediates/Celine.txt
            - mixed/intermediates/Raphael.txt
            - mixed/intermediates/Mateo.txt

          translated:
            - mixed/intermediates/fr-Celine.txt
            - mixed/intermediates/es-Mateo.txt
            - mixed/intermediates/it-Raphael.txt

          collected: mixed/intermediates/mixed-greetings.txt

          batch_report: mixed/mixed-report.txt
        ```

    Each translated file is named after its own language, and the collected file contains one greeting in each language:

    ```console title="results/mixed/intermediates/mixed-greetings.txt"
    Bonjour, Celine!
    Hola, Mateo!
    Ciao, Raphael!
    ```

    To learn more about passing several values together through a process, see the [Metadata and Meta Maps](../side_quests/metadata/index.md#14-group-everything-as-a-single-input-to-the-process) side quest.

---

## Challenge 4: Give each name its own character

Start from `hello-containers.nf` as it stands at the end of Part 5:

```bash
cp solutions/5-hello-containers/hello-containers-2.nf hello-challenge.nf
```

The copy uses the modules in `modules/` and the settings in `nextflow.config` from Parts 4 and 5.
This challenge changes the modules, which `hello-containers.nf` also uses, so make a copy of the `modules` directory first if you want to keep a working version of `hello-containers.nf`.

The third column of `data/names.csv` assigns a character to each name: `squirrel` for Mateo, `turtle` for Raphael and `moose` for Celine.
Make the workflow produce one ASCII art file per name, each drawn by the character from its row, instead of a single drawing of all the greetings.
Keep the collected file and the report as they are.

??? tip "Hint"

    Carry the character along with each greeting, in a [tuple](https://nextflow.io/docs/latest/process.html#input-tuples-tuple).
    Map each CSV row to `[row[0], row[2]]`, declare process inputs like `tuple val(name), val(character)` and outputs like `tuple path("${name}.txt"), val(character)`, and call `cowpy` on the output of `translateGreeting` rather than on the collected file.
    Where only the file is needed, remove the character with `map`.
    Challenge 3 uses the same approach.

??? solution

    In the workflow, build `[name, character]` pairs, pass them through the processes, and call `cowpy` on each translated greeting:

    ```groovy title="hello-challenge.nf" linenums="18" hl_lines="4-5 7 9 13-15 18-19"
    workflow {

        main:
        // create a channel of [name, character] pairs from a CSV file
        person_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { row -> [row[0], row[2]] }
        // emit a greeting
        sayHello(person_ch)
        // translate the greeting
        translateGreeting(sayHello.out, params.language)
        // collect all the greetings into one file
        collectGreetings(translateGreeting.out.map { greeting_file, _character -> greeting_file }.collect(), params.batch)
        // generate ASCII art of each greeting with its own character
        cowpy(translateGreeting.out)

        publish:
        greetings = sayHello.out.map { greeting_file, _character -> greeting_file }
        translated = translateGreeting.out.map { greeting_file, _character -> greeting_file }
        collected = collectGreetings.out.outfile
        batch_report = collectGreetings.out.report
        cowpy_art = cowpy.out
    }
    ```

    In the `map` closures, the underscore at the start of `_character` tells Nextflow that the variable is deliberately unused, so `nextflow lint` doesn't warn about it.

    The `character` parameter is no longer used, so you can remove it from the `params` block and from `nextflow.config`.

    In the modules, the inputs and outputs carry the character alongside the file:

    ```groovy title="modules/sayHello.nf" linenums="6" hl_lines="2 5"
        input:
        tuple val(name), val(character)

        output:
        tuple path("${name}.txt"), val(character)
    ```

    ```groovy title="modules/translateGreeting.nf" linenums="6" hl_lines="2 6"
        input:
        tuple path(input_file), val(character)
        val language

        output:
        tuple path("${language}-${input_file}"), val(character)
    ```

    ```groovy title="modules/cowpy.nf" linenums="7" hl_lines="2"
        input:
        tuple path(input_file), val(character)
    ```

    Run the workflow:

    ```bash
    nextflow run hello-challenge.nf -profile test
    ```

    ??? success "Command output"

        ```console
        N E X T F L O W   ~  version 26.04.4

        Launching `hello-challenge.nf` [stupefied_kowalevski] revision: 4c3e11a405

        executor >  local (10)
        [ab/d91fd4] sayHello (1)          | 3 of 3 ✔
        [2f/fe7a16] translateGreeting (3) | 3 of 3 ✔
        [fa/242a27] collectGreetings      | 1 of 1 ✔
        [57/9f025a] cowpy (2)             | 3 of 3 ✔

        Outputs:

          /workspaces/training/hello-nextflow/results

          greetings:
            - batch/intermediates/Raphael.txt
            - batch/intermediates/Celine.txt
            - batch/intermediates/Mateo.txt

          translated:
            - batch/intermediates/es-Celine.txt
            - batch/intermediates/es-Raphael.txt
            - batch/intermediates/es-Mateo.txt

          collected: batch/intermediates/batch-greetings.txt

          batch_report: batch/batch-report.txt

          cowpy_art:
            - batch/cowpy-es-Mateo.txt
            - batch/cowpy-es-Celine.txt
            - batch/cowpy-es-Raphael.txt
        ```

    `cowpy` now runs three times, once per name, and each file is drawn by that name's character:

    ```console title="results/batch/cowpy-es-Mateo.txt"
     ______________
    < Hola, Mateo! >
     --------------
      \
         \
                      _ _
           | \__/|  .~    ~.
           /oo `./      .'
          {o__,   \    {
            / .  . )    \
            `-` '-' \    }
           .(   _(   )_.'
          '---.~_ _ _|
    ```

    ```console title="results/batch/cowpy-es-Raphael.txt"
     ________________
    < Hola, Raphael! >
     ----------------
        \                                  ___-------___
         \                             _-~~             ~~-_
          \                         _-~                    /~-_
                 /^\__/^\         /~  \                   /    \
               /|  O|| O|        /      \_______________/        \
              | |___||__|      /       /                \          \
              |          \    /      /                    \          \
              |   (_______) /______/                        \_________ \
              |         / /         \                      /            \
               \         \^\\         \                  /               \     /
                 \         ||           \______________/      _-_       //\__//
                   \       ||------_-~~-_ ------------- \ --/~   ~\    || __/
                     ~-----||====/~     |==================|       |/~~~~~
                      (_(__/  ./     /                    \_\      \.
                             (_(___/                         \_____)_)
    ```

    ```console title="results/batch/cowpy-es-Celine.txt"
     _______________
    < Hola, Celine! >
     ---------------
      \
       \   \_\_    _/_/
        \      \__/
               (oo)\_______
               (__)\       )\/\
                   ||----w |
                   ||     ||
    ```

    The solution directory contains its own copy of the modified modules, so you can run the solution directly:

    ```bash
    nextflow run solutions/challenges/challenge-4/hello-challenge.nf -profile test
    ```

    To learn more about carrying metadata such as the character alongside each file, see the [Metadata and Meta Maps](../side_quests/metadata/index.md#15-use-a-meta-map-file-interface) side quest.

---

## What's next?

You've completed the Hello Nextflow course.
See [Next steps](next_steps.md) for where to go from here.
