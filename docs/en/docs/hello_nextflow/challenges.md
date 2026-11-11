# Challenges

The bonus exercises at the end of each part can be solved with what that part covers.
The challenges on this page go further on purpose: each one asks you to extend the Hello pipeline with a Nextflow feature that the course does not cover.

Each challenge starts from the end state of one part of the course.
Copy the solution file for that state to `hello-challenge.nf` in the `hello-nextflow` directory, as shown in each challenge, and make your changes in that file.
Your own course scripts stay as they are.

The hint for each challenge introduces the new feature and links to the Nextflow documentation.
Each solution ends with a pointer to the side quest that covers the feature in depth.

Here is what each challenge covers:

| Challenge                                                                                                              | Starts from the end of | New feature                                       | Side quest                                                               |
| ---------------------------------------------------------------------------------------------------------------------- | ---------------------- | ------------------------------------------------- | ------------------------------------------------------------------------ |
| [1. Read a CSV file with a header row](#challenge-1-read-a-csv-file-with-a-header-row)                                 | Part 2                 | `splitCsv()` with a header row                    | [Splitting and Grouping](../side_quests/splitting_and_grouping/index.md) |
| [2. Greet only the French names](#challenge-2-greet-only-the-french-names)                                             | Part 2                 | The `filter()` operator                           | [Splitting and Grouping](../side_quests/splitting_and_grouping/index.md) |
| [3. Greet each name in its own language](#challenge-3-greet-each-name-in-its-own-language)                             | Part 3                 | Tuples in channels and process inputs and outputs | [Metadata and Meta Maps](../side_quests/metadata/index.md)               |
| [4. Give each name its own character](#challenge-4-give-each-name-its-own-character)                                   | Part 5                 | Tuples across modules                             | [Metadata and Meta Maps](../side_quests/metadata/index.md)               |
| [5. Package the greeting steps as a reusable workflow](#challenge-5-package-the-greeting-steps-as-a-reusable-workflow) | Part 5                 | Named workflows with `take` and `emit`            | [Workflows of Workflows](../side_quests/workflows_of_workflows/index.md) |
| [6. Translate into two languages in one run](#challenge-6-translate-into-two-languages-in-one-run)                     | Part 5                 | Including a module under an alias                 | [Workflows of Workflows](../side_quests/workflows_of_workflows/index.md) |

The solution files are in `solutions/challenges/`, one directory per challenge.
To run a solution directly, launch it from the `hello-nextflow` directory.
Challenges 4 to 6 also need the `nextflow.config` from the end of Part 5 in that directory, which you already have if you completed the course; otherwise, copy it with `cp solutions/5-hello-modules/nextflow.config .` first.

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

## Challenge 3: Greet each name in its own language

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
    `collectGreetings` still takes a language for its file names; give it a label such as `'native'`, since each name is greeted in its own language.

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

    In the workflow, keep both columns in the channel, call `translateGreeting` with the output of `sayHello` alone, and give `collectGreetings` a fixed language label:

    ```groovy title="hello-challenge.nf" linenums="90" hl_lines="1-2 4 6-8 10"
        // create a channel of [name, language] pairs from a CSV file
        person_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { row -> [row[0], row[1]] }
        // emit a greeting
        sayHello(person_ch)
        // translate each greeting into its own language
        translateGreeting(sayHello.out)
        // collect all the greetings into one file
        collectGreetings(translateGreeting.out.collect(), params.batch, 'native')
    ```

    The greetings are now in several languages, so `params.language` no longer describes them.
    `collectGreetings` still needs a language for its file names, so it gets the label `'native'` instead.
    The `language` parameter is no longer used, so remove it from the `params` block.

    `sayHello` now emits pairs, so publish only the file from each pair:

    ```groovy title="hello-challenge.nf" linenums="101" hl_lines="2"
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

        Launching `hello-challenge.nf` [spontaneous_avogadro] revision: a822eaf8ec

        executor >  local (7)
        [04/fd2c8b] sayHello (1)          | 3 of 3 ✔
        [62/668f35] translateGreeting (1) | 3 of 3 ✔
        [cf/c8206c] collectGreetings      | 1 of 1 ✔

        Outputs:

          /workspaces/training/hello-nextflow/results

          greetings:
            - mixed/intermediates/Mateo.txt
            - mixed/intermediates/Celine.txt
            - mixed/intermediates/Raphael.txt

          translated:
            - mixed/intermediates/fr-Celine.txt
            - mixed/intermediates/it-Raphael.txt
            - mixed/intermediates/es-Mateo.txt

          collected: mixed/intermediates/mixed-native-greetings.txt

          batch_report: mixed/mixed-native-report.txt
        ```

    Each translated file is named after its own language, and the collected file contains one greeting in each language:

    ```console title="results/mixed/intermediates/mixed-native-greetings.txt"
    Bonjour, Celine!
    Ciao, Raphael!
    Hola, Mateo!
    ```

    To learn more about passing several values together through a process, see the [Metadata and Meta Maps](../side_quests/metadata/index.md#14-group-everything-as-a-single-input-to-the-process) side quest.

---

## Challenge 4: Give each name its own character

Start from `main.nf` and its modules as they stand at the end of Part 5:

```bash
cp solutions/5-hello-modules/main.nf hello-challenge.nf
cp -r modules challenge-modules
```

At the end of Part 5, the pipeline project consists of `main.nf`, the modules in `modules/`, and the parameter values and profiles in `nextflow.config`.
This challenge changes the modules, which `main.nf` and Challenges 5 and 6 also use, so the second command copies them to `challenge-modules/`.
Point the `include` declarations in `hello-challenge.nf` at `./challenge-modules/...` and edit only those copies, so that `main.nf` and the other challenges keep working.
As in Part 5 Exercise 2, Nextflow finds each module from the path in its `include` declaration.
`hello-challenge.nf` uses the same configuration as `main.nf`.

The third column of `data/names.csv` assigns a character to each name: `squirrel` for Mateo, `turtle` for Raphael and `moose` for Celine.
Make the workflow produce one ASCII art file per name, each drawn by the character from its row, instead of a single drawing of all the greetings.
Keep the collected file and the report as they are.

??? tip "Hint"

    Carry the character along with each greeting, in a [tuple](https://nextflow.io/docs/latest/process.html#input-tuples-tuple).
    Map each CSV row to `[row[0], row[2]]`, declare process inputs like `tuple val(name), val(character)` and outputs like `tuple path("${name}.txt"), val(character)`, and call `cowpy` on the output of `translateGreeting` rather than on the collected file.
    Where only the file is needed, remove the character with `map`.
    Challenge 3 uses the same approach.

??? solution

    In `hello-challenge.nf`, include the modules from `challenge-modules/`:

    ```groovy title="hello-challenge.nf" linenums="3" hl_lines="2-5"
    // Include modules
    include { sayHello } from './challenge-modules/sayHello.nf'
    include { translateGreeting } from './challenge-modules/translateGreeting.nf'
    include { collectGreetings } from './challenge-modules/collectGreetings.nf'
    include { cowpy } from './challenge-modules/cowpy.nf'
    ```

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
        collectGreetings(translateGreeting.out.map { greeting_file, _character -> greeting_file }.collect(), params.batch, params.language)
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

    The `character` parameter is no longer used, so you can remove it from the `params` block of `hello-challenge.nf`.
    Keep its values in `nextflow.config`, since `main.nf` still uses them.
    The `tux` value from the `test` profile has no effect on `hello-challenge.nf`: each drawing uses the character from its own row.

    In the modules in `challenge-modules/`, the inputs and outputs carry the character alongside the file:

    ```groovy title="challenge-modules/sayHello.nf" linenums="6" hl_lines="2 5"
        input:
        tuple val(name), val(character)

        output:
        tuple path("${name}.txt"), val(character)
    ```

    ```groovy title="challenge-modules/translateGreeting.nf" linenums="6" hl_lines="2 6"
        input:
        tuple path(input_file), val(character)
        val language

        output:
        tuple path("${language}-${input_file}"), val(character)
    ```

    ```groovy title="challenge-modules/cowpy.nf" linenums="9" hl_lines="2"
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

        Launching `hello-challenge.nf` [disturbed_stonebraker] revision: 987e9b3a29

        executor >  local (10)
        [cd/04bc7d] sayHello (2)          | 3 of 3 ✔
        [38/b93a35] translateGreeting (3) | 3 of 3 ✔
        [96/aa2a12] collectGreetings      | 1 of 1 ✔
        [ec/b3710e] cowpy (3)             | 3 of 3 ✔

        Outputs:

          /workspaces/training/hello-nextflow/results

          greetings:
            - test/intermediates/Mateo.txt
            - test/intermediates/Celine.txt
            - test/intermediates/Raphael.txt

          translated:
            - test/intermediates/es-Celine.txt
            - test/intermediates/es-Mateo.txt
            - test/intermediates/es-Raphael.txt

          collected: test/intermediates/test-es-greetings.txt

          batch_report: test/test-es-report.txt

          cowpy_art:
            - test/cowpy-es-Mateo.txt
            - test/cowpy-es-Celine.txt
            - test/cowpy-es-Raphael.txt
        ```

    `cowpy` now runs three times, once per name, and each file is drawn by that name's character:

    ```console title="results/test/cowpy-es-Mateo.txt"
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

    ```console title="results/test/cowpy-es-Raphael.txt"
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

    ```console title="results/test/cowpy-es-Celine.txt"
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

    The solution directory contains its own `challenge-modules/` with the modified modules, so you can run the solution directly:

    ```bash
    nextflow run solutions/challenges/challenge-4/hello-challenge.nf -profile test
    ```

    The [Metadata and Meta Maps](../side_quests/metadata/index.md) side quest covers this pattern in depth.
    Its [section 1.4](../side_quests/metadata/index.md#14-group-everything-as-a-single-input-to-the-process) groups several values into one tuple input, as you did here.
    [Section 1.5](../side_quests/metadata/index.md#15-use-a-meta-map-file-interface) then solves a problem much like this one, with a `cowpy` character per row, using a meta map: a single map that holds all of a sample's metadata next to its file, so that adding a field doesn't change the input and output declarations of every process.

---

## Challenge 5: Package the greeting steps as a reusable workflow

Start from `main.nf` as it stands at the end of Part 5:

```bash
cp solutions/5-hello-modules/main.nf hello-challenge.nf
```

The copy uses the same modules and configuration as `main.nf`.
This challenge doesn't change them.

Move the three greeting steps, `sayHello`, `translateGreeting` and `collectGreetings`, into a named workflow called `GREET`, in a new file `workflows/greet.nf`.
`GREET` takes the channel of names, the language and the batch name as inputs, and outputs the greetings, the translated greetings, the collected file and the report.
The entry workflow in `hello-challenge.nf` then calls `GREET`, passes the collected file to `cowpy`, and publishes the same outputs as before.

??? tip "Hint"

    A named workflow declares its inputs in a `take:` section and its outputs in an `emit:` section, around the `main:` section that holds its logic (see [Takes and emits](https://nextflow.io/docs/latest/workflow.html#takes-and-emits)):

    ```groovy
    workflow GREET {

        take:
        names_ch

        main:
        sayHello(names_ch)

        emit:
        greetings = sayHello.out
    }
    ```

    You [include](https://nextflow.io/docs/latest/reference/syntax.html#include) a named workflow the same way as a process, and its named outputs are available as `GREET.out.greetings` and so on.
    The paths in `include` declarations are relative to the file that contains them, so `workflows/greet.nf` includes the modules from `../modules/`.

??? solution

    Create `workflows/greet.nf` with the three greeting steps:

    ```groovy title="workflows/greet.nf" linenums="1" hl_lines="2-4 6 9-11 17 19 22-25"
    // Include modules
    include { sayHello } from '../modules/sayHello.nf'
    include { translateGreeting } from '../modules/translateGreeting.nf'
    include { collectGreetings } from '../modules/collectGreetings.nf'

    workflow GREET {

        take:
        names_ch
        language
        batch_name

        main:
        // emit a greeting
        sayHello(names_ch)
        // translate the greeting
        translateGreeting(sayHello.out, language)
        // collect all the greetings into one file
        collectGreetings(translateGreeting.out.collect(), batch_name, language)

        emit:
        greetings = sayHello.out
        translated = translateGreeting.out
        collected = collectGreetings.out.outfile
        report = collectGreetings.out.report
    }
    ```

    The `include` paths start with `../` because the file is in `workflows/`, one level below the `modules` directory.
    Inside `GREET`, the language and the batch name come from the `take:` inputs instead of `params`, and the language goes to both `translateGreeting` and `collectGreetings`.

    In `hello-challenge.nf`, include `GREET` instead of the three modules, call it, and publish its outputs:

    ```groovy title="hello-challenge.nf" linenums="1" hl_lines="4 25 27 30-33"
    #!/usr/bin/env nextflow

    // Include the greeting workflow and the cowpy module
    include { GREET } from './workflows/greet.nf'
    include { cowpy } from './modules/cowpy.nf'

    /*
     * Pipeline parameters
     */
    params {
        input: Path
        language: String
        batch: String
        character: String
    }

    workflow {

        main:
        // create a channel for inputs from a CSV file
        name_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { row -> row[0] }
        // greet, translate and collect
        GREET(name_ch, params.language, params.batch)
        // generate ASCII art of the greetings with cowpy
        cowpy(GREET.out.collected, params.character)

        publish:
        greetings = GREET.out.greetings
        translated = GREET.out.translated
        collected = GREET.out.collected
        batch_report = GREET.out.report
        cowpy_art = cowpy.out
    }
    ```

    The `output` block stays as it is.

    Only the entry workflow can publish outputs, so `GREET` emits everything that needs to be published, including the intermediate greetings and translations, and the entry workflow publishes them from `GREET.out`.

    Run the workflow:

    ```bash
    nextflow run hello-challenge.nf -profile test
    ```

    ??? success "Command output"

        ```console
        N E X T F L O W   ~  version 26.04.4

        Launching `hello-challenge.nf` [romantic_lavoisier] revision: 60ca931fc8

        executor >  local (8)
        [af/bfca39] GREET:sayHello (1)          | 3 of 3 ✔
        [ef/595f11] GREET:translateGreeting (3) | 3 of 3 ✔
        [74/8d690f] GREET:collectGreetings      | 1 of 1 ✔
        [d0/0236db] cowpy                       | 1 of 1 ✔

        Outputs:

          /workspaces/training/hello-nextflow/results

          greetings:
            - test/intermediates/Raphael.txt
            - test/intermediates/Celine.txt
            - test/intermediates/Mateo.txt

          translated:
            - test/intermediates/es-Raphael.txt
            - test/intermediates/es-Mateo.txt
            - test/intermediates/es-Celine.txt

          collected: test/intermediates/test-es-greetings.txt

          batch_report: test/test-es-report.txt

          cowpy_art: test/cowpy-test-es-greetings.txt
        ```

    The processes called inside `GREET` are listed as `GREET:sayHello` and so on, and the published outputs are the same as with `main.nf`.

    The solution directory contains its own copy of the modules and the configuration, so you can run the solution directly:

    ```bash
    nextflow run solutions/challenges/challenge-5/hello-challenge.nf -profile test
    ```

    To learn more about composing workflows, see the [Workflows of Workflows](../side_quests/workflows_of_workflows/index.md#12-make-the-workflow-composable) side quest.

---

## Challenge 6: Translate into two languages in one run

Start from `main.nf` as it stands at the end of Part 5:

```bash
cp solutions/5-hello-modules/main.nf hello-challenge.nf
```

The copy uses the same modules and configuration as `main.nf`.

`translateGreeting` takes the language as an input, so translating into a second language looks like a matter of calling it a second time.
However, a process can only be called once in a workflow.
A second call to `translateGreeting` fails with this error:

```console
Process 'translateGreeting' has been already used -- If you need to reuse the same component, include it with a different name or include it in a different workflow context
```

Add a `second_language` parameter with the default value `fr`, and make the workflow translate the greetings into both languages.
Collect the second-language greetings into their own file and report, and publish them alongside the existing outputs.
Keep `cowpy` on the greetings in the first language.

??? tip "Hint"

    Include the same process a second time under a different name with `include { ... as ... }` (see [Calling processes and workflows](https://nextflow.io/docs/latest/workflow.html#calling-processes-and-workflows)):

    ```groovy
    include { translateGreeting as translateSecond } from './modules/translateGreeting.nf'
    ```

    `translateSecond` is then a separate process that runs the same code, and you can call it with its own inputs.
    `collectGreetings` needs the same treatment.
    Since it includes the language in its file names, the second collection doesn't clash with the first.

??? solution

    Add the default value of the new parameter to the `params` section of `nextflow.config`:

    ```groovy title="nextflow.config" linenums="6" hl_lines="3"
    params {
        language = 'es'
        second_language = 'fr'
        batch = 'batch'
        character = 'turkey'
    }
    ```

    `main.nf` doesn't declare `second_language`, and still runs as before with this configuration.

    In `hello-challenge.nf`, include `translateGreeting` and `collectGreetings` a second time under aliases:

    ```groovy title="hello-challenge.nf" linenums="3" hl_lines="5-6"
    // Include modules
    include { sayHello } from './modules/sayHello.nf'
    include { translateGreeting } from './modules/translateGreeting.nf'
    include { collectGreetings } from './modules/collectGreetings.nf'
    include { translateGreeting as translateSecond } from './modules/translateGreeting.nf'
    include { collectGreetings as collectSecond } from './modules/collectGreetings.nf'
    include { cowpy } from './modules/cowpy.nf'
    ```

    Declare the new parameter:

    ```groovy title="hello-challenge.nf" linenums="14" hl_lines="4"
    params {
        input: Path
        language: String
        second_language: String
        batch: String
        character: String
    }
    ```

    In the workflow, call the aliased processes on the same greetings, and publish their outputs:

    ```groovy title="hello-challenge.nf" linenums="22" hl_lines="14-17 27-29"
    workflow {

        main:
        // create a channel for inputs from a CSV file
        name_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { row -> row[0] }
        // emit a greeting
        sayHello(name_ch)
        // translate the greeting
        translateGreeting(sayHello.out, params.language)
        // collect all the greetings into one file
        collectGreetings(translateGreeting.out.collect(), params.batch, params.language)
        // translate the greeting into the second language
        translateSecond(sayHello.out, params.second_language)
        // collect the second-language greetings into one file
        collectSecond(translateSecond.out.collect(), params.batch, params.second_language)
        // generate ASCII art of the greetings with cowpy
        cowpy(collectGreetings.out.outfile, params.character)

        publish:
        greetings = sayHello.out
        translated = translateGreeting.out
        collected = collectGreetings.out.outfile
        batch_report = collectGreetings.out.report
        cowpy_art = cowpy.out
        translated_second = translateSecond.out
        collected_second = collectSecond.out.outfile
        batch_report_second = collectSecond.out.report
    }
    ```

    Add the three new outputs to the `output` block, in the same directories as their first-language counterparts:

    ```groovy title="hello-challenge.nf" linenums="70" hl_lines="5-16"
        cowpy_art {
            path "${params.batch}"
            mode 'copy'
        }
        translated_second {
            path "${params.batch}/intermediates"
            mode 'copy'
        }
        collected_second {
            path "${params.batch}/intermediates"
            mode 'copy'
        }
        batch_report_second {
            path "${params.batch}"
            mode 'copy'
        }
    }
    ```

    Run the workflow:

    ```bash
    nextflow run hello-challenge.nf -profile test
    ```

    ??? success "Command output"

        ```console
        N E X T F L O W   ~  version 26.04.4

        Launching `hello-challenge.nf` [jovial_marconi] revision: 72805d564b

        executor >  local (12)
        [cc/272efd] sayHello (1)          | 3 of 3 ✔
        [c2/4b2b04] translateGreeting (1) | 3 of 3 ✔
        [da/4d2f05] collectGreetings      | 1 of 1 ✔
        [1a/fd2576] translateSecond (3)   | 3 of 3 ✔
        [33/7c130b] collectSecond         | 1 of 1 ✔
        [2e/f89690] cowpy                 | 1 of 1 ✔

        Outputs:

          /workspaces/training/hello-nextflow/results

          greetings:
            - test/intermediates/Mateo.txt
            - test/intermediates/Celine.txt
            - test/intermediates/Raphael.txt

          translated:
            - test/intermediates/es-Raphael.txt
            - test/intermediates/es-Celine.txt
            - test/intermediates/es-Mateo.txt

          collected: test/intermediates/test-es-greetings.txt

          batch_report: test/test-es-report.txt

          cowpy_art: test/cowpy-test-es-greetings.txt

          translated_second:
            - test/intermediates/fr-Celine.txt
            - test/intermediates/fr-Mateo.txt
            - test/intermediates/fr-Raphael.txt

          collected_second: test/intermediates/test-fr-greetings.txt

          batch_report_second: test/test-fr-report.txt
        ```

    `translateSecond` and `collectSecond` are listed as separate processes, and the second collected file contains the French greetings:

    ```console title="results/test/intermediates/test-fr-greetings.txt"
    Bonjour, Celine!
    Bonjour, Mateo!
    Bonjour, Raphael!
    ```

    Aliases work for workflows too.
    With the `GREET` workflow from Challenge 5, you could include it twice, as in `include { GREET as GREET_SECOND } from './workflows/greet.nf'`, and call each one with its own language.

    The solution directory contains its own copy of the modules and the configuration, so you can run the solution directly:

    ```bash
    nextflow run solutions/challenges/challenge-6/hello-challenge.nf -profile test
    ```

    To learn more about including workflows, with or without an alias, see the [Workflows of Workflows](../side_quests/workflows_of_workflows/index.md#key-patterns) side quest.

---

## What's next?

Return to the [course summary](next_steps.md#next-steps-to-build-your-skills) to choose a course for further study.
