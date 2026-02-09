# Teil 2: Hello für nf-core umschreiben

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

In diesem zweiten Teil des Hello nf-core Trainingskurses zeigen wir dir, wie du eine nf-core-kompatible Version der Pipeline erstellst, die im [Hello Nextflow](../hello_nextflow/index.md) Einsteigerkurs entwickelt wurde.

Du hast im ersten Abschnitt des Trainings bemerkt, dass nf-core-Pipelines einer ziemlich aufwendigen Struktur mit vielen zusätzlichen Dateien folgen.
All das von Grund auf neu zu erstellen wäre sehr mühsam, daher hat die nf-core-Community Werkzeuge entwickelt, um dies stattdessen aus einer Vorlage zu erstellen und den Prozess zu vereinfachen.

Wir zeigen dir, wie du diese Werkzeuge verwendest, um ein Pipeline-Gerüst zu erstellen und dann bestehenden 'regulären' Pipeline-Code auf das nf-core-Gerüst anzupassen.

Falls du mit der Hello-Pipeline nicht vertraut bist oder eine Auffrischung gebrauchen könntest, siehe [diese Infoseite](../info/hello_pipeline.md).

---

## 1. Ein neues Pipeline-Projekt erstellen

Zuerst erstellen wir das Gerüst für die neue Pipeline.

!!! note "Hinweis"

    Stelle sicher, dass du dich im Verzeichnis `hello-nf-core` in deinem Terminal befindest.

### 1.1. Das vorlagenbasierte Pipeline-Erstellungswerkzeug ausführen

Beginnen wir mit der Erstellung einer neuen Pipeline mit dem Befehl `nf-core pipelines create`.
Dies erstellt ein neues Pipeline-Gerüst unter Verwendung der nf-core-Basisvorlage, angepasst mit einem Pipeline-Namen, einer Beschreibung und einem\*einer Autor\*in.

```bash
nf-core pipelines create
```

Die Ausführung dieses Befehls öffnet eine Text User Interface (TUI) für die Pipeline-Erstellung:

<div style="text-align: center;">
    <iframe width="560" height="315" src="https://www.youtube.com/embed/VwjXNXONHlY?si=d0HkFSISnKn76TeI" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" allowfullscreen="" data-ruffle-polyfilled=""></iframe>
</div>

Diese TUI fordert dich auf, grundlegende Informationen über deine Pipeline anzugeben und bietet dir eine Auswahl an Funktionen, die du in dein Pipeline-Gerüst einbeziehen oder ausschließen kannst.

- Klicke auf dem Willkommensbildschirm auf **Let's go!**.
- Klicke auf dem Bildschirm `Choose pipeline type` auf **Custom**.
- Gib deine Pipeline-Details wie folgt ein (ersetze `< DEIN NAME >` durch deinen eigenen Namen) und klicke dann auf **Next**.

```
[ ] GitHub organisation: core
[ ] Workflow name: hello
[ ] A short description of your pipeline: A basic nf-core style version of Hello Nextflow
[ ] Name of the main author(s): < DEIN NAME >
```

- Setze auf dem Bildschirm Template features `Toggle all features` auf **off** und **aktiviere** dann selektiv die folgenden. Überprüfe deine Auswahl und klicke auf **Continue**.

```
[ ] Add testing profiles
[ ] Use nf-core components
[ ] Use nf-schema
[ ] Add configuration files
[ ] Add documentation
```

- Klicke auf dem Bildschirm `Final details` auf **Finish**. Warte, bis die Pipeline erstellt wurde, und klicke dann auf **Continue**.
- Klicke auf dem Bildschirm Create GitHub repository auf **Finish without creating a repo**. Dies zeigt Anweisungen zum späteren Erstellen eines GitHub-Repositorys an. Ignoriere diese und klicke auf **Close**.

Sobald die TUI geschlossen wird, solltest du die folgende Konsolenausgabe sehen.

??? success "Befehlsausgabe"

    ```console
                                              ,--./,-.
              ___     __   __   __   ___     /,-._.--~\
        |\ | |__  __ /  ` /  \ |__) |__         }  {
        | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                              `._,._,'

        nf-core/tools version 3.4.1 - https://nf-co.re


    INFO     Launching interactive nf-core pipeline creation tool.
    ```

Es gibt keine explizite Bestätigung in der Konsolenausgabe, dass die Pipeline-Erstellung funktioniert hat, aber du solltest ein neues Verzeichnis namens `core-hello` sehen.

Zeige den Inhalt des neuen Verzeichnisses an, um zu sehen, wie viel Arbeit du dir durch die Verwendung der Vorlage erspart hast.

```bash
tree core-hello
```

??? abstract "Verzeichnisinhalt"

    ```console
    core-hello/
    ├── assets
    │   ├── samplesheet.csv
    │   └── schema_input.json
    ├── conf
    │   ├── base.config
    │   ├── modules.config
    │   ├── test.config
    │   └── test_full.config
    ├── docs
    │   ├── output.md
    │   ├── README.md
    │   └── usage.md
    ├── main.nf
    ├── modules.json
    ├── nextflow.config
    ├── nextflow_schema.json
    ├── README.md
    ├── subworkflows
    │   ├── local
    │   │   └── utils_nfcore_hello_pipeline
    │   │       └── main.nf
    │   └── nf-core
    │       ├── utils_nextflow_pipeline
    │       │   ├── main.nf
    │       │   ├── meta.yml
    │       │   └── tests
    │       │       ├── main.function.nf.test
    │       │       ├── main.function.nf.test.snap
    │       │       ├── main.workflow.nf.test
    │       │       └── nextflow.config
    │       ├── utils_nfcore_pipeline
    │       │   ├── main.nf
    │       │   ├── meta.yml
    │       │   └── tests
    │       │       ├── main.function.nf.test
    │       │       ├── main.function.nf.test.snap
    │       │       ├── main.workflow.nf.test
    │       │       ├── main.workflow.nf.test.snap
    │       │       └── nextflow.config
    │       └── utils_nfschema_plugin
    │           ├── main.nf
    │           ├── meta.yml
    │           └── tests
    │               ├── main.nf.test
    │               ├── nextflow.config
    │               └── nextflow_schema.json
    └── workflows
        └── hello.nf

    14 directories, 34 files
    ```

Das sind eine Menge Dateien!

Hoffentlich erkennst du viele davon als dieselben, die wir bei der Untersuchung der `nf-core/demo`-Pipeline-Struktur kennengelernt haben.
Aber keine Sorge, wenn du dich noch etwas verloren fühlst; wir werden im Verlauf dieses Trainings gemeinsam die wichtigen Teile durchgehen.

!!! note "Hinweis"

    Ein wichtiger Unterschied im Vergleich zur `nf-core/demo`-Pipeline, die wir im ersten Teil dieses Trainings untersucht haben, ist, dass es kein `modules`-Verzeichnis gibt.
    Das liegt daran, dass wir uns nicht dafür entschieden haben, eines der Standard-nf-core-Module einzubeziehen.

### 1.2. Testen, dass das Gerüst funktionsfähig ist

Glaub es oder nicht, auch wenn du noch keine Module hinzugefügt hast, um echte Arbeit zu leisten, kann das Pipeline-Gerüst tatsächlich mit dem Testprofil ausgeführt werden, genauso wie wir die `nf-core/demo`-Pipeline ausgeführt haben.

```bash
nextflow run ./core-hello -profile docker,test --outdir core-hello-results
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.04.3

    Launching `./core-hello/main.nf` [scruffy_marconi] DSL2 - revision: b9e9b3b8de

    Downloading plugin nf-schema@2.5.1
    Input/output options
      input                     : https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv
      outdir                    : core-hello-results

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      trace_report_suffix       : 2025-11-21_04-47-18

    Core Nextflow options
      runName                   : scruffy_marconi
      containerEngine           : docker
      launchDir                 : /workspaces/training/hello-nf-core
      workDir                   : /workspaces/training/hello-nf-core/work
      projectDir                : /workspaces/training/hello-nf-core/core-hello
      userName                  : root
      profile                   : docker,test
      configFiles               : /workspaces/training/hello-nf-core/core-hello/nextflow.config

    !! Only displaying parameters that differ from the pipeline defaults !!
    ------------------------------------------------------
    -[core/hello] Pipeline completed successfully-
    ```

Dies zeigt dir, dass die gesamte grundlegende Verkabelung vorhanden ist.
Wo sind also die Ausgaben? Gibt es überhaupt welche?

Tatsächlich wurde ein neues Ergebnisverzeichnis namens `core-hello-results` erstellt, das die Standard-Ausführungsberichte enthält:

```bash
tree core-hello-results
```

??? abstract "Verzeichnisinhalt"

    ```console
    core-hello-results
    └── pipeline_info
        ├── execution_report_2025-11-21_04-47-18.html
        ├── execution_timeline_2025-11-21_04-47-18.html
        ├── execution_trace_2025-11-21_04-47-18.txt
        ├── hello_software_versions.yml
        ├── params_2025-11-21_04-47-18.json
        └── pipeline_dag_2025-11-21_04-47-18.html

    1 directory, 6 files
    ```

Du kannst einen Blick in die Berichte werfen, um zu sehen, was ausgeführt wurde, und die Antwort lautet: überhaupt nichts!

![leerer Ausführungs-Timeline-Bericht](./img/execution_timeline_empty.png)

Schauen wir uns an, was tatsächlich im Code steht.

### 1.3. Den Platzhalter-Workflow untersuchen

Wenn du in die Datei `main.nf` schaust, siehst du, dass sie einen Workflow namens `HELLO` aus `workflows/hello` importiert.

Dies entspricht dem `workflows/demo.nf`-Workflow, dem wir in Teil 1 begegnet sind, und dient als Platzhalter-Workflow für unseren gewünschten Workflow, mit bereits vorhandener nf-core-Funktionalität.

```groovy title="core-hello/workflows/hello.nf" linenums="1" hl_lines="15 17 19 35"
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow HELLO {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    main:

    ch_versions = channel.empty()

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'hello_software_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    emit:
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
```

Im Vergleich zu einem einfachen Nextflow-Workflow wie dem in [Hello Nextflow](../hello_nextflow/index.md) entwickelten, wirst du einige Dinge bemerken, die hier neu sind (hervorgehobene Zeilen oben):

- Der Workflow-Block hat einen Namen
- Workflow-Eingaben werden mit dem Schlüsselwort `take:` deklariert und die Kanal-Konstruktion wird in den übergeordneten Workflow verschoben
- Workflow-Inhalt wird innerhalb eines `main:`-Blocks platziert
- Ausgaben werden mit dem Schlüsselwort `emit:` deklariert

Dies sind optional Funktionen von Nextflow, die den Workflow **komponierbar** machen, was bedeutet, dass er von innerhalb eines anderen Workflows aufgerufen werden kann.

!!! note "Komponierbare Workflows im Detail"

    Die [Workflows of Workflows](../side_quests/workflows_of_workflows.md) Side Quest untersucht die Workflow-Komposition viel ausführlicher, einschließlich wie man mehrere Workflows zusammensetzt und komplexe Datenflüsse zwischen ihnen verwaltet. Wir führen die Komponierbarkeit hier ein, weil sie eine grundlegende Anforderung der nf-core-Template-Architektur ist, die verschachtelte Workflows verwendet, um Pipeline-Initialisierung, den Haupt-Analyse-Workflow und Abschlussaufgaben in separate, wiederverwendbare Komponenten zu organisieren.

Wir müssen die relevante Logik aus unserem gewünschten Workflow in diese Struktur einfügen.
Der erste Schritt dafür ist, unseren ursprünglichen Workflow komponierbar zu machen.

### Fazit

Du weißt jetzt, wie du ein Pipeline-Gerüst mit nf-core-Werkzeugen erstellst.

### Wie geht es weiter?

Lerne, wie du einen einfachen Workflow komponierbar machst als Vorbereitung darauf, ihn nf-core-kompatibel zu machen.

---

## 2. Den ursprünglichen Hello Nextflow Workflow komponierbar machen

Jetzt ist es Zeit, mit der Integration unseres Workflows in das nf-core-Gerüst zu beginnen.
Zur Erinnerung: Wir arbeiten mit dem Workflow aus unserem [Hello Nextflow](../hello_nextflow/index.md) Trainingskurs.

!!! tip "Tipp"

    Falls du mit dieser Pipeline nicht vertraut bist oder eine Auffrischung gebrauchen könntest, siehe [Die Hello-Pipeline](../info/hello_pipeline.md).

Wir stellen dir eine saubere, voll funktionsfähige Kopie des abgeschlossenen Hello Nextflow Workflows im Verzeichnis `original-hello` zusammen mit seinen Modulen und der Standard-CSV-Datei zur Verfügung, die als Eingabe erwartet wird.

```bash
tree original-hello/
```

??? abstract "Verzeichnisinhalt"

    ```console
    original-hello/
    ├── hello.nf
    ├── modules
    │   ├── collectGreetings.nf
    │   ├── convertToUpper.nf
    │   ├── cowpy.nf
    │   └── sayHello.nf
    └── nextflow.config

    1 directory, 6 files
    ```

Führe ihn gerne aus, um dich zu überzeugen, dass er funktioniert:

```bash
nextflow run original-hello/hello.nf
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.04.3

    Launching `original-hello/hello.nf` [goofy_babbage] DSL2 - revision: e9e72441e9

    executor >  local (8)
    [a4/081cec] sayHello (1)       | 3 of 3 ✔
    [e7/7e9058] convertToUpper (3) | 3 of 3 ✔
    [0c/17263b] collectGreetings   | 1 of 1 ✔
    [94/542280] cowpy              | 1 of 1 ✔
    ```

Öffnen wir die Workflow-Datei `hello.nf`, um den Code zu untersuchen, der unten vollständig gezeigt wird (ohne die Prozesse, die in Modulen sind):

```groovy title="original-hello/hello.nf" linenums="1"
#!/usr/bin/env nextflow

/*
* Pipeline parameters
*/
params.greeting = 'greetings.csv'
params.batch = 'test-batch'
params.character = 'turkey'

// Include modules
include { sayHello } from './modules/sayHello.nf'
include { convertToUpper } from './modules/convertToUpper.nf'
include { collectGreetings } from './modules/collectGreetings.nf'
include { cowpy } from './modules/cowpy.nf'

workflow {

  // create a channel for inputs from a CSV file
  greeting_ch = channel.fromPath(params.greeting)
                      .splitCsv()
                      .map { line -> line[0] }

  // emit a greeting
  sayHello(greeting_ch)

  // convert the greeting to uppercase
  convertToUpper(sayHello.out)

  // collect all the greetings into one file
  collectGreetings(convertToUpper.out.collect(), params.batch)

  // generate ASCII art of the greetings with cowpy
  cowpy(collectGreetings.out.outfile, params.character)
}
```

Wie du sehen kannst, wurde dieser Workflow als einfacher unbenannter Workflow geschrieben, der eigenständig ausgeführt werden kann.
Um ihn von innerhalb eines übergeordneten Workflows ausführbar zu machen, wie es die nf-core-Vorlage erfordert, müssen wir ihn **komponierbar** machen.

Gehen wir die notwendigen Änderungen einzeln durch.

### 2.1. Den Workflow benennen

Geben wir dem Workflow zunächst einen Namen, damit wir von einem übergeordneten Workflow darauf verweisen können.

=== "Danach"

    ```groovy title="original-hello/hello.nf" linenums="16"
    workflow HELLO {
    ```

=== "Vorher"

    ```groovy title="original-hello/hello.nf" linenums="16"
    workflow {
    ```

Für Workflow-Namen gelten die gleichen Konventionen wie für Modulnamen.

### 2.2. Kanal-Konstruktion durch `take` ersetzen

Ersetze nun die Kanal-Konstruktion durch eine einfache `take`-Anweisung, die erwartete Eingaben deklariert.

=== "Danach"

    ```groovy title="original-hello/hello.nf" linenums="18"
        take:
        // Kanal von Begrüßungen
        greeting_ch
    ```

=== "Vorher"

    ```groovy title="original-hello/hello.nf" linenums="18"
        // Einen Kanal für Eingaben aus einer CSV-Datei erstellen
        greeting_ch = channel.fromPath(params.greeting)
                            .splitCsv()
                            .map { line -> line[0] }
    ```

Dies überlässt die Details, wie die Eingaben bereitgestellt werden, dem übergeordneten Workflow.

Während wir dabei sind, können wir auch die Zeile `params.greeting = 'greetings.csv'` auskommentieren

=== "Danach"

    ```groovy title="original-hello/hello.nf" linenums="3" hl_lines="4"
        /*
        * Pipeline parameters
        */
        //params.greeting = 'greetings.csv'
        params.batch = 'test-batch'
        params.character = 'turkey'
    ```

=== "Vorher"

    ```groovy title="original-hello/hello.nf" linenums="3" hl_lines="4"
        /*
        * Pipeline parameters
        */
        params.greeting = 'greetings.csv'
        params.batch = 'test-batch'
        params.character = 'turkey'
    ```

!!! note "Hinweis"

    Wenn du die Nextflow Language Server Extension installiert hast, wird der Syntax-Checker deinen Code mit roten Wellenlinien markieren.
    Das liegt daran, dass du, wenn du eine `take:`-Anweisung einfügst, auch ein `main:` haben musst.

    Das fügen wir im nächsten Schritt hinzu.

### 2.3. Workflow-Operationen mit `main`-Anweisung einleiten

Füge als Nächstes eine `main`-Anweisung vor den restlichen Operationen hinzu, die im Hauptteil des Workflows aufgerufen werden.

=== "Danach"

    ```groovy title="original-hello/hello.nf" linenums="22" hl_lines="1"
        main:

        // Eine Begrüßung ausgeben
        sayHello(greeting_ch)

        // Die Begrüßung in Großbuchstaben umwandeln
        convertToUpper(sayHello.out)

        // Alle Begrüßungen in einer Datei sammeln
        collectGreetings(convertToUpper.out.collect(), params.batch)

        // ASCII-Art der Begrüßungen mit cowpy generieren
        cowpy(collectGreetings.out.outfile, params.character)
    ```

=== "Vorher"

    ```groovy title="original-hello/hello.nf" linenums="21"
        // Eine Begrüßung ausgeben
        sayHello(greeting_ch)

        // Die Begrüßung in Großbuchstaben umwandeln
        convertToUpper(sayHello.out)

        // Alle Begrüßungen in einer Datei sammeln
        collectGreetings(convertToUpper.out.collect(), params.batch)

        // ASCII-Art der Begrüßungen mit cowpy generieren
        cowpy(collectGreetings.out.outfile, params.character)
    ```

Dies sagt im Grunde 'das ist, was dieser Workflow _tut_'.

### 2.4. `emit`-Anweisung hinzufügen

Füge schließlich eine `emit`-Anweisung hinzu, die deklariert, was die finalen Ausgaben des Workflows sind.

```groovy title="original-hello/hello.nf" linenums="35"
    emit:
    cowpy_hellos = cowpy.out
```

Dies ist eine komplett neue Ergänzung zum Code im Vergleich zum ursprünglichen Workflow.

### 2.5. Zusammenfassung der abgeschlossenen Änderungen

Wenn du alle Änderungen wie beschrieben vorgenommen hast, sollte dein Workflow jetzt so aussehen:

```groovy title="original-hello/hello.nf" linenums="1" hl_lines="16 18-20 22 36-37"
#!/usr/bin/env nextflow

/*
* Pipeline parameters
*/
// params.greeting = 'greetings.csv'
params.batch = 'test-batch'
params.character = 'turkey'

// Include modules
include { sayHello } from './modules/sayHello.nf'
include { convertToUpper } from './modules/convertToUpper.nf'
include { collectGreetings } from './modules/collectGreetings.nf'
include { cowpy } from './modules/cowpy.nf'

workflow HELLO {

    take:
    // Kanal von Begrüßungen
    greeting_ch

    main:

    // Eine Begrüßung ausgeben
    sayHello(greeting_ch)

    // Die Begrüßung in Großbuchstaben umwandeln
    convertToUpper(sayHello.out)

    // Alle Begrüßungen in einer Datei sammeln
    collectGreetings(convertToUpper.out.collect(), params.batch)

    // ASCII-Art der Begrüßungen mit cowpy generieren
    cowpy(collectGreetings.out.outfile, params.character)

    emit:
    cowpy_hellos = cowpy.out
}
```

Dies beschreibt alles, was Nextflow braucht, AUSSER was in den Eingabekanal eingespeist werden soll.
Das wird im übergeordneten Workflow definiert, auch **Einstiegspunkt**-Workflow genannt.

### 2.6. Einen Dummy-Einstiegspunkt-Workflow erstellen

Bevor wir unseren komponierbaren Workflow in das komplexe nf-core-Gerüst integrieren, sollten wir überprüfen, ob er korrekt funktioniert.
Wir können einen einfachen Dummy-Einstiegspunkt-Workflow erstellen, um den komponierbaren Workflow isoliert zu testen.

Erstelle eine leere Datei namens `main.nf` im selben `original-hello`-Verzeichnis.

```bash
touch original-hello/main.nf
```

Kopiere den folgenden Code in die Datei `main.nf`.

```groovy title="original-hello/main.nf" linenums="1"
#!/usr/bin/env nextflow

// Den Workflow-Code aus der Datei hello.nf importieren
include { HELLO } from './hello.nf'

// Eingabeparameter deklarieren
params.greeting = 'greetings.csv'

workflow {
  // Einen Kanal für Eingaben aus einer CSV-Datei erstellen
  greeting_ch = channel.fromPath(params.greeting)
                      .splitCsv()
                      .map { line -> line[0] }

  // Den importierten Workflow auf dem Kanal von Begrüßungen aufrufen
  HELLO(greeting_ch)

  // Die vom Workflow ausgegebenen Ausgaben anzeigen
  HELLO.out.view { output -> "Output: $output" }
}
```

Hier gibt es zwei wichtige Beobachtungen:

- Die Syntax zum Aufrufen des importierten Workflows ist im Wesentlichen dieselbe wie die Syntax zum Aufrufen von Modulen.
- Alles, was mit dem Einlesen der Eingaben in den Workflow zu tun hat (Eingabeparameter und Kanal-Konstruktion), wird jetzt in diesem übergeordneten Workflow deklariert.

!!! note "Hinweis"

    Die Benennung der Einstiegspunkt-Workflow-Datei `main.nf` ist eine Konvention, keine Anforderung.

    Wenn du dieser Konvention folgst, kannst du die Angabe des Workflow-Dateinamens in deinem `nextflow run`-Befehl weglassen.
    Nextflow sucht automatisch nach einer Datei namens `main.nf` im Ausführungsverzeichnis.

    Du kannst die Einstiegspunkt-Workflow-Datei jedoch auch anders benennen, wenn du möchtest.
    Stelle in diesem Fall sicher, dass du den Workflow-Dateinamen in deinem `nextflow run`-Befehl angibst.

### 2.7. Testen, dass der Workflow läuft

Wir haben endlich alle Teile, die wir brauchen, um zu überprüfen, dass der komponierbare Workflow funktioniert.
Führen wir ihn aus!

```bash
nextflow run ./original-hello
```

Hier siehst du den Vorteil der Verwendung der `main.nf`-Namenskonvention.
Wenn wir den Einstiegspunkt-Workflow `something_else.nf` genannt hätten, hätten wir `nextflow run original-hello/something_else.nf` machen müssen.

Wenn du alle Änderungen korrekt vorgenommen hast, sollte dies bis zum Abschluss laufen.

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.04.3

    Launching `original-hello/main.nf` [friendly_wright] DSL2 - revision: 1ecd2d9c0a

    executor >  local (8)
    [24/c6c0d8] HELLO:sayHello (3)       | 3 of 3 ✔
    [dc/721042] HELLO:convertToUpper (3) | 3 of 3 ✔
    [48/5ab2df] HELLO:collectGreetings   | 1 of 1 ✔
    [e3/693b7e] HELLO:cowpy              | 1 of 1 ✔
    Output: /workspaces/training/hello-nf-core/work/e3/693b7e48dc119d0c54543e0634c2e7/cowpy-COLLECTED-test-batch-output.txt
    ```

Das bedeutet, wir haben unseren HELLO-Workflow erfolgreich auf komponierbar aktualisiert.

### Fazit

Du weißt, wie du einen Workflow komponierbar machst, indem du ihm einen Namen gibst und `take`-, `main`- und `emit`-Anweisungen hinzufügst, und wie du ihn von einem Einstiegspunkt-Workflow aufrufst.

### Wie geht es weiter?

Lerne, wie du einen einfachen komponierbaren Workflow auf das nf-core-Gerüst aufpfropfst.

---

## 3. Die aktualisierte Workflow-Logik in den Platzhalter-Workflow einfügen

Jetzt, da wir überprüft haben, dass unser komponierbarer Workflow korrekt funktioniert, kehren wir zum nf-core-Pipeline-Gerüst zurück, das wir in Abschnitt 1 erstellt haben.
Wir möchten den komponierbaren Workflow, den wir gerade entwickelt haben, in die nf-core-Template-Struktur integrieren, sodass das Endergebnis ungefähr so aussehen sollte.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nf-core/img/core-hello.svg"
</figure>

Wie machen wir das also? Schauen wir uns den aktuellen Inhalt des `HELLO`-Workflows in `core-hello/workflows/hello.nf` (dem nf-core-Gerüst) an.

```groovy title="core-hello/workflows/hello.nf" linenums="1"
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow HELLO {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    main:

    ch_versions = channel.empty()

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'hello_software_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    emit:
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
```

Insgesamt macht dieser Code sehr wenig außer etwas Verwaltungsarbeit, die mit der Erfassung der Version aller Software-Tools zu tun hat, die in der Pipeline ausgeführt werden.

Wir müssen den relevanten Code aus der komponierbaren Version des ursprünglichen Workflows hinzufügen, den wir in Abschnitt 2 entwickelt haben.

Wir werden dies in den folgenden Phasen angehen:

1. Die Module kopieren und Modul-Importe einrichten
2. Die `take`-Deklaration unverändert lassen
3. Die Workflow-Logik zum `main`-Block hinzufügen
4. Den `emit`-Block aktualisieren

!!! note "Hinweis"

    Wir werden die Versionserfassung für diesen ersten Durchgang ignorieren und uns später in diesem Training ansehen, wie man das verdrahtet.

### 3.1. Die Module kopieren und Modul-Importe einrichten

Die vier Prozesse aus unserem Hello Nextflow Workflow sind als Module in `original-hello/modules/` gespeichert.
Wir müssen diese Module in die nf-core-Projektstruktur (unter `core-hello/modules/local/`) kopieren und Import-Anweisungen zur nf-core-Workflow-Datei hinzufügen.

Kopieren wir zuerst die Modul-Dateien von `original-hello/` nach `core-hello/`:

```bash
mkdir -p core-hello/modules/local/
cp original-hello/modules/* core-hello/modules/local/.
```

Du solltest jetzt das Verzeichnis der Module unter `core-hello/` aufgelistet sehen.

```bash
tree core-hello/modules
```

??? abstract "Verzeichnisinhalt"

    ```console
    core-hello/modules
    └── local
        ├── collectGreetings.nf
        ├── convertToUpper.nf
        ├── cowpy.nf
        └── sayHello.nf

    1 directory, 4 files
    ```

Richten wir nun die Modul-Import-Anweisungen ein.

Dies waren die Import-Anweisungen im `original-hello/hello.nf`-Workflow:

```groovy title="original-hello/hello.nf" linenums="9"
// Include modules
include { sayHello } from './modules/sayHello.nf'
include { convertToUpper } from './modules/convertToUpper.nf'
include { collectGreetings } from './modules/collectGreetings.nf'
include { cowpy } from './modules/cowpy.nf'
```

Öffne die Datei `core-hello/workflows/hello.nf` und übertrage diese Import-Anweisungen wie unten gezeigt hinein.

=== "Danach"

    ```groovy title="core-hello/workflows/hello.nf" linenums="1" hl_lines="8-11"
    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    include { paramsSummaryMap       } from 'plugin/nf-schema'
    include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
    include { sayHello               } from '../modules/local/sayHello.nf'
    include { convertToUpper         } from '../modules/local/convertToUpper.nf'
    include { collectGreetings       } from '../modules/local/collectGreetings.nf'
    include { cowpy                  } from '../modules/local/cowpy.nf'
    ```

=== "Vorher"

    ```groovy title="core-hello/workflows/hello.nf" linenums="1"
    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    include { paramsSummaryMap       } from 'plugin/nf-schema'
    include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
    ```

Zwei weitere interessante Beobachtungen hier:

- Wir haben die Formatierung der Import-Anweisungen angepasst, um der nf-core-Stilkonvention zu folgen.
- Wir haben die relativen Pfade zu den Modulen aktualisiert, um widerzuspiegeln, dass sie jetzt auf einer anderen Verschachtelungsebene gespeichert sind.

### 3.2. Die `take`-Deklaration unverändert lassen

Das nf-core-Projekt hat viel vorgefertigte Funktionalität rund um das Konzept des Samplesheets, das typischerweise eine CSV-Datei mit Spaltendaten ist.
Da das im Wesentlichen das ist, was unsere `greetings.csv`-Datei ist, behalten wir die aktuelle `take`-Deklaration unverändert bei und aktualisieren einfach den Namen des Eingabekanals im nächsten Schritt.

```groovy title="core-hello/workflows/hello.nf" linenums="21"
    take:
    ch_samplesheet // channel: samplesheet read in from --input
```

Die Eingabeverarbeitung wird upstream von diesem Workflow durchgeführt (nicht in dieser Code-Datei).

### 3.3. Die Workflow-Logik zum `main`-Block hinzufügen

Jetzt, da unsere Module für den Workflow verfügbar sind, können wir die Workflow-Logik in den `main`-Block einfügen.

Zur Erinnerung, dies ist der relevante Code im ursprünglichen Workflow, der sich nicht viel geändert hat, als wir ihn komponierbar gemacht haben (wir haben nur die Zeile `main:` hinzugefügt):

```groovy title="original-hello/hello.nf" linenums="22"
    main:

    // Eine Begrüßung ausgeben
    sayHello(greeting_ch)

    // Die Begrüßung in Großbuchstaben umwandeln
    convertToUpper(sayHello.out)

    // Alle Begrüßungen in einer Datei sammeln
    collectGreetings(convertToUpper.out.collect(), params.batch)

    // ASCII-Art der Begrüßungen mit cowpy generieren
    cowpy(collectGreetings.out.outfile, params.character)
```

Wir müssen den Code, der nach `main:` kommt, in die neue Version des Workflows kopieren.

Es gibt bereits etwas Code darin, der mit der Erfassung der Versionen der Tools zu tun hat, die vom Workflow ausgeführt werden. Das lassen wir vorerst in Ruhe (wir werden uns später mit den Tool-Versionen befassen).
Wir behalten die `ch_versions = channel.empty()`-Initialisierung oben, fügen dann unsere Workflow-Logik ein und behalten den Versionszusammenstellungs-Code am Ende.
Diese Reihenfolge macht Sinn, weil in einer echten Pipeline die Prozesse Versionsinformationen ausgeben würden, die dem `ch_versions`-Kanal hinzugefügt würden, während der Workflow läuft.

=== "Danach"

    ```groovy title="core-hello/workflows/hello.nf" linenums="19" hl_lines="10-20"
    workflow HELLO {

        take:
        ch_samplesheet // channel: samplesheet read in from --input

        main:

        ch_versions = Channel.empty()

        // Eine Begrüßung ausgeben
        sayHello(greeting_ch)

        // Die Begrüßung in Großbuchstaben umwandeln
        convertToUpper(sayHello.out)

        // Alle Begrüßungen in einer Datei sammeln
        collectGreetings(convertToUpper.out.collect(), params.batch)

        // ASCII-Art der Begrüßungen mit cowpy generieren
        cowpy(collectGreetings.out.outfile, params.character)

        //
        // Collate and save software versions
        //
        softwareVersionsToYAML(ch_versions)
            .collectFile(
                storeDir: "${params.outdir}/pipeline_info",
                name:  'hello_software_'  + 'versions.yml',
                sort: true,
                newLine: true
            ).set { ch_collated_versions }


        emit:
        versions       = ch_versions                 // channel: [ path(versions.yml) ]

    }
    ```

=== "Vorher"

    ```groovy title="core-hello/workflows/hello.nf" linenums="19"
    workflow HELLO {

        take:
        ch_samplesheet // channel: samplesheet read in from --input
        main:

        ch_versions = Channel.empty()

        //
        // Collate and save software versions
        //
        softwareVersionsToYAML(ch_versions)
            .collectFile(
                storeDir: "${params.outdir}/pipeline_info",
                name:  'hello_software_'  + 'versions.yml',
                sort: true,
                newLine: true
            ).set { ch_collated_versions }


        emit:
        versions       = ch_versions                 // channel: [ path(versions.yml) ]

    }
    ```

Du wirst bemerken, dass wir auch eine Leerzeile vor `main:` hinzugefügt haben, um den Code lesbarer zu machen.

Das sieht großartig aus, aber wir müssen noch den Namen des Kanals aktualisieren, den wir an den `sayHello()`-Prozess übergeben, von `greeting_ch` zu `ch_samplesheet`, wie unten gezeigt, um zu dem zu passen, was unter dem Schlüsselwort `take:` geschrieben steht.

=== "Danach"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26"
        // Eine Begrüßung ausgeben (aktualisiert, um die nf-core-Konvention für Samplesheets zu verwenden)
        sayHello(ch_samplesheet)
    ```

=== "Vorher"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26"
        // Eine Begrüßung ausgeben
        sayHello(greeting_ch)
    ```

Jetzt ist die Workflow-Logik korrekt verdrahtet.

### 3.4. Den `emit`-Block aktualisieren

Schließlich müssen wir den `emit`-Block aktualisieren, um die Deklaration der finalen Ausgaben des Workflows einzuschließen.

=== "Danach"

    ```groovy title="core-hello/workflows/hello.nf" linenums="55" hl_lines="2"
        emit:
        cowpy_hellos   = cowpy.out
        versions       = ch_versions                 // channel: [ path(versions.yml) ]
    ```

=== "Vorher"

    ```groovy title="core-hello/workflows/hello.nf" linenums="55"
        emit:
        versions       = ch_versions                 // channel: [ path(versions.yml) ]
    ```

Dies schließt die Änderungen ab, die wir am HELLO-Workflow selbst vornehmen müssen.
An diesem Punkt haben wir die gesamte Code-Struktur erreicht, die wir implementieren wollten.

### Fazit

Du weißt, wie du die Kernstücke eines komponierbaren Workflows in einen nf-core-Platzhalter-Workflow einfügst.

### Wie geht es weiter?

Lerne, wie du die Eingabeverarbeitung im nf-core-Pipeline-Gerüst anpasst.

---

## 4. Die Eingabeverarbeitung anpassen

Jetzt, da wir unsere Workflow-Logik erfolgreich in das nf-core-Gerüst integriert haben, müssen wir noch ein kritisches Element ansprechen: sicherstellen, dass unsere Eingabedaten korrekt verarbeitet werden.
Die nf-core-Vorlage kommt mit ausgeklügelter Eingabeverarbeitung, die für komplexe Genomik-Datensätze entwickelt wurde, also müssen wir sie anpassen, um mit unserer einfacheren `greetings.csv`-Datei zu funktionieren.

### 4.1. Identifizieren, wo Eingaben verarbeitet werden

Der erste Schritt ist herauszufinden, wo die Eingabeverarbeitung durchgeführt wird.

Du erinnerst dich vielleicht, dass wir, als wir den Hello Nextflow Workflow umgeschrieben haben, um ihn komponierbar zu machen, die Eingabeparameter-Deklaration eine Ebene nach oben verschoben haben, in den `main.nf`-Einstiegspunkt-Workflow.
Schauen wir uns also den obersten `main.nf`-Einstiegspunkt-Workflow an, der als Teil des Pipeline-Gerüsts erstellt wurde:

```groovy title="core-hello/main.nf" linenums="1"
#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    core/hello
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/core/hello
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { HELLO  } from './workflows/hello'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_hello_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_hello_pipeline'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow CORE_HELLO {

    take:
    samplesheet // channel: samplesheet read in from --input

    main:

    //
    // WORKFLOW: Run pipeline
    //
    HELLO (
        samplesheet
    )
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input
    )

    //
    // WORKFLOW: Run main workflow
    //
    CORE_HELLO (
        PIPELINE_INITIALISATION.out.samplesheet
    )
    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.outdir,
        params.monochrome_logs,
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
```

Das nf-core-Projekt macht starken Gebrauch von verschachtelten Subworkflows, daher kann dieser Teil beim ersten Ansatz etwas verwirrend sein.

Was hier wichtig ist, ist, dass zwei Workflows definiert sind:

- `CORE_HELLO` ist ein dünner Wrapper zum Ausführen des HELLO-Workflows, den wir gerade in `core-hello/workflows/hello.nf` angepasst haben.
- Ein unbenannter Workflow, der `CORE_HELLO` sowie zwei andere Subworkflows aufruft, `PIPELINE_INITIALISATION` und `PIPELINE_COMPLETION`.

Hier ist ein Diagramm, wie sie zueinander in Beziehung stehen:

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nf-core/img/hello-nested-workflows.svg"
</figure>

Wichtig ist, dass wir auf dieser Ebene keinen Code finden können, der einen Eingabekanal konstruiert, nur Verweise auf ein Samplesheet, das über den Parameter `--input` bereitgestellt wird.

Ein bisschen Herumstöbern zeigt, dass die Eingabeverarbeitung vom Subworkflow `PIPELINE_INITIALISATION` durchgeführt wird, was passend ist, der aus `core-hello/subworkflows/local/utils_nfcore_hello_pipeline/main.nf` importiert wird.

Wenn wir diese Datei öffnen und nach unten scrollen, kommen wir zu diesem Code-Block:

```groovy title="core-hello/subworkflows/local/utils_nfcore_hello_pipeline/main.nf" linenums="76"
    //
    // Create channel from input file provided through params.input
    //

    channel
        .fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
        .map {
            meta, fastq_1, fastq_2 ->
                if (!fastq_2) {
                    return [ meta.id, meta + [ single_end:true ], [ fastq_1 ] ]
                } else {
                    return [ meta.id, meta + [ single_end:false ], [ fastq_1, fastq_2 ] ]
                }
        }
        .groupTuple()
        .map { samplesheet ->
            validateInputSamplesheet(samplesheet)
        }
        .map {
            meta, fastqs ->
                return [ meta, fastqs.flatten() ]
        }
        .set { ch_samplesheet }

    emit:
    samplesheet = ch_samplesheet
    versions    = ch_versions
```

Dies ist die Kanal-Factory, die das Samplesheet parst und es in einer Form weitergibt, die bereit ist, vom HELLO-Workflow konsumiert zu werden.

!!! note "Hinweis"

    Die obige Syntax ist etwas anders als das, was wir zuvor verwendet haben, aber im Grunde ist dies:

    ```groovy
    channel.<...>.set { ch_samplesheet }
    ```

    äquivalent zu diesem:

    ```groovy
    ch_samplesheet = channel.<...>
    ```

Dieser Code beinhaltet einige Parse- und Validierungsschritte, die sehr spezifisch für das Beispiel-Samplesheet sind, das mit der nf-core-Pipeline-Vorlage enthalten ist, das zum Zeitpunkt des Schreibens sehr domänenspezifisch und für unser einfaches Pipeline-Projekt nicht geeignet ist.

### 4.2. Den vorlagenbasierten Eingabekanal-Code ersetzen

Die gute Nachricht ist, dass die Bedürfnisse unserer Pipeline viel einfacher sind, sodass wir all das durch den Kanal-Konstruktionscode ersetzen können, den wir im ursprünglichen Hello Nextflow Workflow entwickelt haben.

Zur Erinnerung, so sah die Kanal-Konstruktion aus (wie im Lösungsverzeichnis zu sehen):

```groovy title="solutions/composable-hello/main.nf" linenums="10" hl_lines="4"
    // Einen Kanal für Eingaben aus einer CSV-Datei erstellen
    greeting_ch = channel.fromPath(params.greeting)
        .splitCsv()
        .map { line -> line[0] }
```

Wir müssen das also nur in den Initialisierungs-Workflow einfügen, mit kleinen Änderungen: Wir aktualisieren den Kanalnamen von `greeting_ch` zu `ch_samplesheet` und den Parameternamen von `params.greeting` zu `params.input` (siehe hervorgehobene Zeile).

=== "Danach"

    ```groovy title="core-hello/subworkflows/local/utils_nfcore_hello_pipeline/main.nf" linenums="76" hl_lines="5-7"
        //
        // Create channel from input file provided through params.input
        //

        ch_samplesheet = channel.fromPath(params.input)
            .splitCsv()
            .map { line -> line[0] }

        emit:
        samplesheet = ch_samplesheet
        versions    = ch_versions
    ```

=== "Vorher"

    ```groovy title="core-hello/subworkflows/local/utils_nfcore_hello_pipeline/main.nf" linenums="76" hl_lines="5-23"
        //
        // Create channel from input file provided through params.input
        //

        channel
            .fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
            .map {
                meta, fastq_1, fastq_2 ->
                    if (!fastq_2) {
                        return [ meta.id, meta + [ single_end:true ], [ fastq_1 ] ]
                    } else {
                        return [ meta.id, meta + [ single_end:false ], [ fastq_1, fastq_2 ] ]
                    }
            }
            .groupTuple()
            .map { samplesheet ->
                validateInputSamplesheet(samplesheet)
            }
            .map {
                meta, fastqs ->
                    return [ meta, fastqs.flatten() ]
            }
            .set { ch_samplesheet }

        emit:
        samplesheet = ch_samplesheet
        versions    = ch_versions
    ```

Das vervollständigt die Änderungen, die wir vornehmen müssen, damit die Eingabeverarbeitung funktioniert.

In seiner aktuellen Form wird uns dies nicht erlauben, die eingebauten Fähigkeiten von nf-core für Schema-Validierung zu nutzen, aber das können wir später hinzufügen.
Für jetzt konzentrieren wir uns darauf, es so einfach wie möglich zu halten, um zu etwas zu kommen, das wir erfolgreich auf Testdaten ausführen können.

### 4.3. Das Testprofil aktualisieren

Apropos Testdaten und Parameter, aktualisieren wir das Testprofil für diese Pipeline, um das `greetings.csv`-Mini-Samplesheet anstelle des Beispiel-Samplesheets zu verwenden, das in der Vorlage bereitgestellt wird.

Unter `core-hello/conf` finden wir zwei vorlagenbasierte Testprofile: `test.config` und `test_full.config`, die dazu gedacht sind, eine kleine Datenprobe und eine vollständige zu testen.
Angesichts des Zwecks unserer Pipeline gibt es nicht wirklich einen Grund, ein vollständiges Testprofil einzurichten, also kannst du `test_full.config` gerne ignorieren oder löschen.
Wir werden uns darauf konzentrieren, `test.config` so einzurichten, dass es auf unserer `greetings.csv`-Datei mit ein paar Standardparametern läuft.

#### 4.3.1. Die `greetings.csv`-Datei kopieren

Zuerst müssen wir die `greetings.csv`-Datei an einen geeigneten Ort in unserem Pipeline-Projekt kopieren.
Typischerweise werden kleine Testdateien im `assets`-Verzeichnis gespeichert, also kopieren wir die Datei aus unserem Arbeitsverzeichnis.

```bash
cp greetings.csv core-hello/assets/.
```

Jetzt ist die `greetings.csv`-Datei bereit, als Testeingabe verwendet zu werden.

#### 4.3.2. Die `test.config`-Datei aktualisieren

Jetzt können wir die `test.config`-Datei wie folgt aktualisieren:

=== "Danach"

    ```groovy title="core-hello/conf/test.config" linenums="21" hl_lines="6-10"
    params {
        config_profile_name        = 'Test profile'
        config_profile_description = 'Minimal test dataset to check pipeline function'

        // Input data
        input  = "${projectDir}/assets/greetings.csv"

        // Other parameters
        batch     = 'test'
        character = 'tux'
    }
    ```

=== "Vorher"

    ```groovy title="core-hello/conf/test.config" linenums="21" hl_lines="6-8"
    params {
        config_profile_name        = 'Test profile'
        config_profile_description = 'Minimal test dataset to check pipeline function'

        // Input data
        // TODO nf-core: Specify the paths to your test data on nf-core/test-datasets
        // TODO nf-core: Give any required params for the test so that command line flags are not needed
        input  = params.pipelines_testdata_base_path + 'viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv'
    }
    ```

Wichtige Punkte:

- **Verwendung von `${projectDir}`**: Dies ist eine implizite Nextflow-Variable, die auf das Verzeichnis zeigt, in dem sich das Haupt-Workflow-Skript befindet (das Pipeline-Root). Die Verwendung stellt sicher, dass der Pfad funktioniert, unabhängig davon, von wo aus die Pipeline ausgeführt wird.
- **Absolute Pfade**: Durch die Verwendung von `${projectDir}` erstellen wir einen absoluten Pfad, was wichtig für Testdaten ist, die mit der Pipeline ausgeliefert werden.
- **Testdaten-Speicherort**: nf-core-Pipelines speichern Testdaten typischerweise im `assets/`-Verzeichnis innerhalb des Pipeline-Repositorys für kleine Testdateien oder verweisen auf externe Testdatensätze für größere Dateien.

Und während wir dabei sind, verschärfen wir die Standard-Ressourcenlimits, um sicherzustellen, dass dies auf sehr einfachen Maschinen läuft (wie den minimalen VMs in Github Codespaces):

=== "Danach"

    ```groovy title="core-hello/config/test.config" linenums="13" hl_lines="3-4"
    process {
        resourceLimits = [
            cpus: 2,
            memory: '4.GB',
            time: '1.h'
        ]
    }
    ```

=== "Vorher"

    ```groovy title="core-hello/config/test.config" linenums="13" hl_lines="3-4"
    process {
        resourceLimits = [
            cpus: 4,
            memory: '15.GB',
            time: '1.h'
        ]
    }
    ```

Dies vervollständigt die Code-Änderungen, die wir vornehmen müssen.

### 4.4. Die Pipeline mit dem Testprofil ausführen

Das war viel, aber wir können endlich versuchen, die Pipeline auszuführen!
Beachte, dass wir `--validate_params false` zur Befehlszeile hinzufügen müssen, weil wir die Validierung noch nicht eingerichtet haben (das kommt später).

```bash
nextflow run core-hello --outdir core-hello-results -profile test,docker --validate_params false
```

Wenn du alle Änderungen korrekt vorgenommen hast, sollte dies bis zum Abschluss laufen.

??? success "Befehlsausgabe"

    ```console
     N E X T F L O W   ~  version 25.04.3

    Launching `core-hello/main.nf` [condescending_allen] DSL2 - revision: b9e9b3b8de

    Input/output options
      input                     : /workspaces/training/hello-nf-core/core-hello/assets/greetings.csv
      outdir                    : core-hello-results

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      validate_params           : false
      trace_report_suffix       : 2025-11-21_07-29-37

    Core Nextflow options
      runName                   : condescending_allen
      containerEngine           : docker
      launchDir                 : /workspaces/training/hello-nf-core
      workDir                   : /workspaces/training/hello-nf-core/work
      projectDir                : /workspaces/training/hello-nf-core/core-hello
      userName                  : root
      profile                   : test,docker
      configFiles               : /workspaces/training/hello-nf-core/core-hello/nextflow.config

    !! Only displaying parameters that differ from the pipeline defaults !!
    ------------------------------------------------------
    executor >  local (1)
    [ed/727b7e] CORE_HELLO:HELLO:sayHello (3)       [100%] 3 of 3 ✔
    [45/bb6096] CORE_HELLO:HELLO:convertToUpper (3) [100%] 3 of 3 ✔
    [81/7e2e34] CORE_HELLO:HELLO:collectGreetings   [100%] 1 of 1 ✔
    [96/9442a1] CORE_HELLO:HELLO:cowpy              [100%] 1 of 1 ✔
    -[core/hello] Pipeline completed successfully-
    ```

Wie du sehen kannst, hat dies die typische nf-core-Zusammenfassung am Anfang dank des Initialisierungs-Subworkflows erzeugt, und die Zeilen für jedes Modul zeigen jetzt die vollständigen PIPELINE:WORKFLOW:Modul-Namen.

### 4.5. Die Pipeline-Ausgaben finden

Die Frage ist jetzt: Wo sind die Ausgaben der Pipeline?
Und die Antwort ist ziemlich interessant: Es gibt jetzt zwei verschiedene Orte, an denen man nach den Ergebnissen suchen kann.

Wie du dich vielleicht von früher erinnerst, hat unsere erste Ausführung des neu erstellten Workflows ein Verzeichnis namens `core-hello-results/` erzeugt, das verschiedene Ausführungsberichte und Metadaten enthielt.

```bash
tree core-hello-results
```

??? abstract "Verzeichnisinhalt"

    ```console
    core-hello-results
    └── pipeline_info
        ├── execution_report_2025-11-21_04-47-18.html
        ├── execution_report_2025-11-21_07-29-37.html
        ├── execution_timeline_2025-11-21_04-47-18.html
        ├── execution_timeline_2025-11-21_07-29-37.html
        ├── execution_trace_2025-11-21_04-47-18.txt
        ├── execution_trace_2025-11-21_07-29-37.txt
        ├── hello_software_versions.yml
        ├── params_2025-11-21_04-47-13.json
        ├── params_2025-11-21_07-29-41.json
        └── pipeline_dag_2025-11-21_04-47-18.html
        └── pipeline_dag_2025-11-21_07-29-37.html

    1 directory, 12 files
    ```

Du siehst, wir haben einen weiteren Satz Ausführungsberichte zusätzlich zu denen, die wir von der ersten Ausführung bekommen haben, als der Workflow noch nur ein Platzhalter war.
Diesmal siehst du alle Aufgaben, die wie erwartet ausgeführt wurden.

![Ausführungs-Timeline-Bericht für die Hello-Pipeline](./img/execution_timeline_hello.png)

!!! note "Hinweis"

    Wieder wurden die Aufgaben nicht parallel ausgeführt, weil wir auf einer minimalistischen Maschine in Github Codespaces laufen.
    Um diese parallel laufen zu sehen, versuche die CPU-Zuweisung deines Codespace und die Ressourcenlimits in der Testkonfiguration zu erhöhen.

Das ist großartig, aber unsere tatsächlichen Pipeline-Ergebnisse sind nicht dort!

Folgendes ist passiert: Wir haben nichts an den Modulen selbst geändert, sodass die Ausgaben, die von Modul-Level-`publishDir`-Direktiven verarbeitet werden, immer noch in ein `results`-Verzeichnis gehen, wie in der ursprünglichen Pipeline angegeben.

```bash
tree results
```

??? abstract "Verzeichnisinhalt"

    ```console
    results
    ├── Bonjour-output.txt
    ├── COLLECTED-test-batch-output.txt
    ├── COLLECTED-test-output.txt
    ├── cowpy-COLLECTED-test-batch-output.txt
    ├── cowpy-COLLECTED-test-output.txt
    ├── Hello-output.txt
    ├── Holà-output.txt
    ├── UPPER-Bonjour-output.txt
    ├── UPPER-Hello-output.txt
    └── UPPER-Holà-output.txt

    0 directories, 10 files
    ```

Ah, da sind sie, vermischt mit den Ausgaben früherer Ausführungen der ursprünglichen Hello-Pipeline.

Wenn wir möchten, dass sie ordentlich organisiert sind wie die Ausgaben der Demo-Pipeline, müssen wir ändern, wie wir die zu veröffentlichenden Ausgaben einrichten.
Wir zeigen dir später in diesem Trainingskurs, wie das geht.

<!-- TODO: Update this once we've updated Hello Nextflow to use workflow-level outputs -->

Und da ist es! Es mag wie viel Arbeit erscheinen, um das gleiche Ergebnis wie die ursprüngliche Pipeline zu erreichen, aber du bekommst all diese schönen Berichte automatisch generiert, und du hast jetzt eine solide Grundlage, um zusätzliche Funktionen von nf-core zu nutzen, einschließlich Eingabevalidierung und einiger praktischer Metadaten-Verarbeitungsfähigkeiten, die wir in einem späteren Abschnitt behandeln werden.

---

### Fazit

Du weißt, wie du eine reguläre Nextflow-Pipeline in eine nf-core-Style-Pipeline umwandelst, indem du die nf-core-Vorlage verwendest.
Als Teil davon hast du gelernt, wie man einen Workflow komponierbar macht und wie man die Elemente der nf-core-Vorlage identifiziert, die am häufigsten angepasst werden müssen, wenn man eine benutzerdefinierte nf-core-Style-Pipeline entwickelt.

### Wie geht es weiter?

Mach eine Pause, das war harte Arbeit! Wenn du bereit bist, gehe weiter zu [Teil 3: Ein nf-core-Modul verwenden](./03_use_module.md), um zu lernen, wie du von der Community gepflegte Module aus dem nf-core/modules-Repository nutzen kannst.
