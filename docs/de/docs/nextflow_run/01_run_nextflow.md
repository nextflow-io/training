# Teil 1: Nextflow ausführen

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

In diesem Teil stellen wir die grundlegenden Konzepte zum Ausführen von Nextflow-Pipelines vor.
Wir beginnen mit einem einfachen Hello World-Workflow und arbeiten uns zu einer vollständigen mehrstufigen Pipeline vor, die mehrere Eingaben parallel mit Containern verarbeitet.

---

## 1. Hello World

Der Workflow `1-hello.nf` nimmt eine Begrüßung als Kommandozeilenargument entgegen und schreibt sie in eine Datei.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_world.svg"
</figure>

### 1.1. Den Workflow starten

Führe den folgenden Befehl in deinem Terminal aus.

```bash
nextflow run 1-hello.nf --input 'Hello World!'
```

??? success "Befehlsausgabe"

    ```console hl_lines="6"
    N E X T F L O W   ~  version 26.04.4

    Launching `1-hello.nf` [infallible_volhard] revision: 82d40dbf94

    executor >  local (1)
    [6d/740edd] sayHello | 1 of 1 ✔

    Outputs:

      /workspaces/training/nextflow-run/results

      first_output: 1-hello/Hello World!-output.txt
    ```

Die wichtigste Zeile in der Ausgabe ist die Prozessstatuszeile:

```console
[6d/740edd] sayHello | 1 of 1 ✔
```

Sie zeigt uns, dass der `sayHello`-Prozess einmal erfolgreich ausgeführt wurde.
Das Präfix `[6d/740edd]` ist ein abgekürzter Pfad zum Arbeitsverzeichnis der Aufgabe — dazu gleich mehr.
Der darauffolgende `Outputs:`-Block listet alle Dateien auf, die die Pipeline veröffentlicht hat, beschriftet entsprechend dem `output`-Block, der in [1.4](#14-optional-code-walkthrough) weiter unten behandelt wird.

### 1.2. Die Ausgabe finden

Dieser Workflow ist so konfiguriert, dass er seine Ausgabe in ein `results`-Verzeichnis veröffentlicht.
Nach der Ausführung solltest du die Ausgabe dort finden:

```console title="results/"
results
└── 1-hello
    └── Hello World!-output.txt
```

Öffne die Datei und überprüfe, ob sie `Hello World!` enthält.

### 1.3. Das `work/`-Verzeichnis erkunden

Im Hintergrund erstellt Nextflow für jeden Prozessaufruf ein eindeutiges Aufgabenverzeichnis innerhalb eines Verzeichnisses namens `work/`.
Der in der Konsolenausgabe angezeigte Hash (`[6d/740edd]`) ist der Pfad zu diesem Verzeichnis.

```bash
ls work/6d/740edd*
```

Darin findest du die Ausgabedatei sowie mehrere versteckte Log-Dateien:

- **`.command.sh`**: der genaue Befehl, den Nextflow ausgeführt hat
- **`.command.out`** / **`.command.err`**: stdout und stderr des Prozesses
- **`.command.log`**: kombinierte Log-Ausgabe
- **`.exitcode`**: der Exit-Code des Prozesses

Die Datei `.command.sh` ist beim Debuggen besonders nützlich — sie zeigt genau, was ausgeführt wurde.

### 1.4. Optional: Code-Erklärung

Das Verstehen des Codes ist nicht unbedingt notwendig, wenn du nur Pipelines ausführen möchtest. Wenn du neugierig bist, lohnt sich aber ein Blick.

??? optional "Klicke hier, um den Code zu dieser Übung zu erkunden"

    Öffnen wir `1-hello.nf` und schauen uns die wichtigsten Bestandteile an.

    ```groovy title="1-hello.nf" linenums="1"
    #!/usr/bin/env nextflow

    include { sayHello } from './modules/sayHello.nf'

    /*
     * Pipeline-Parameter
     */
    params {
        input: String
    }

    workflow {

        main:
        // Eine Begrüßung ausgeben
        sayHello(params.input)

        publish:
        first_output = sayHello.out
    }

    output {
        first_output {
            path '1-hello'
            mode 'copy'
        }
    }
    ```

    Wir sehen Folgendes:

    - eine `include`-Anweisung, die auf ein `process`-Modul verweist
    - einen `params`-Block, der Pipeline-Parameter definiert
    - einen `workflow`-Block, der die auszuführende Arbeit beschreibt
    - einen `output`-Block, der beschreibt, was mit den Ausgaben geschehen soll

    Schauen wir uns jeden dieser Teile der Reihe nach an.

    ### Das `process`-Modul

    Die `include`-Anweisung weist Nextflow an, etwas namens `sayHello` aus einer separaten Code-Datei zu laden.

    ```groovy title="1-hello.nf" linenums="3"
    include { sayHello } from './modules/sayHello.nf'
    ```

    In dieser Datei finden wir die Definition eines Prozesses namens `sayHello`:

    ```groovy title="modules/sayHello.nf" linenums="4"
    process sayHello {

        input:
        val greeting

        output:
        path "${greeting}-output.txt"

        script:
        """
        echo '${greeting}' > '${greeting}-output.txt'
        """
    }
    ```

    Ein **process** definiert einen einzelnen Schritt in der Pipeline.
    Er deklariert seine Eingaben, Ausgaben und das auszuführende Skript.
    Der `val`-Qualifier bedeutet, dass die Eingabe ein einfacher Wert ist (String, Zahl usw.).
    Der `path`-Qualifier bedeutet, dass die Ausgabe ein Dateipfad ist.

    Es ist möglich, die Prozessdefinition in der Haupt-Workflow-Datei zu schreiben. Wenn man sie jedoch in separaten Moduldateien hält, werden sie wiederverwendbar: Dasselbe Modul kann von mehreren Workflow-Skripten importiert werden.

    ### Der `params`-Block

    Der `params`-Block deklariert die Kommandozeilenparameter, die der Workflow akzeptiert:

    ```groovy title="1-hello.nf" linenums="8"
    params {
        input: String
    }
    ```

    Jeder hier deklarierte Parameter wird auf der Kommandozeile mit zwei Bindestrichen verfügbar (`--input`).
    Unterstützte Typen sind `String`, `Integer`, `Float`, `Boolean` und `Path`.

    !!! tip "Tipp"

        Workflow-Parameter verwenden immer zwei Bindestriche (`--input`), um sie von Nextflows eigenen CLI-Flags zu unterscheiden, die einen Bindestrich verwenden (z. B. `-resume`).

    ### Der `workflow`-Block

    Der **workflow**-Block definiert die Datenflusslogik: welche Prozesse in welcher Reihenfolge ausgeführt werden.

    ```groovy title="1-hello.nf" linenums="12"
    workflow {

        main:
        // Eine Begrüßung ausgeben
        sayHello(params.input)

        publish:
        first_output = sayHello.out
    }
    ```

    Hier wird nur ein Prozess aufgerufen, daher ist es sehr einfach. Realistischere Beispiele folgen später.

    Der `main:`-Abschnitt ruft den `sayHello`-Prozess mit dem `--input`-Wert auf.
    Der `publish:`-Abschnitt listet auf, welche Ausgaben in das Ergebnisverzeichnis kopiert werden sollen.

    ### Der `output`-Block

    Der `output`-Block am Ende der Datei gibt den Zielpfad und den Kopiermodus an.

    ```groovy title="1-hello.nf" linenums="22"
    output {
        first_output {
            path '1-hello'
            mode 'copy'
        }
    }
    ```

    Jeder benannte Eintrag entspricht einem `publish:`-Label im Workflow und ordnet ihn einem Unterverzeichnis unter `results/` zu.

### Fazit

Du weißt, wie du eine Nextflow-Pipeline ausführst und ihre Ausgaben findest, und du weißt, dass die Arbeit in Aufgabenverzeichnissen unter `work/` ausgeführt wird.

### Wie geht es weiter?

Erfahre, wie Nextflow mehrere Eingaben effizient verarbeitet.

---

## 2. Mehrere Eingaben verarbeiten

Pipelines in der Praxis verarbeiten typischerweise viele Datensätze, nicht nur einen.
Der Workflow `2-inputs.nf` liest aus einer CSV-Datei und führt `sayHello` einmal pro Zeile parallel aus.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-pipeline-multi-inputs-csv.svg"
</figure>

Führen wir zuerst den Workflow aus und schauen uns dann an, welchen Mechanismus Nextflow verwendet, um diese mehreren Eingaben zu verarbeiten.

### 2.1. Den Workflow ausführen

Führe den folgenden Befehl in deinem Terminal aus.

```bash
nextflow run 2-inputs.nf --input data/greetings.csv
```

??? success "Befehlsausgabe"

    ```console hl_lines="6"
    N E X T F L O W   ~  version 26.04.4

    Launching `2-inputs.nf` [nauseous_babbage] revision: b90778224d

    executor >  local (3)
    [66/de7844] sayHello (3) | 3 of 3 ✔

    Outputs:

      /workspaces/training/nextflow-run/results

      first_output:
        - 2-inputs/Bonjour-output.txt
        - 2-inputs/Hello-output.txt
        - 2-inputs/Hola-output.txt
    ```

Das `3 of 3` zeigt uns, dass der `sayHello`-Prozess dreimal aufgerufen wurde, einmal pro Zeile in der CSV-Datei.

Im `results`-Verzeichnis solltest du jetzt drei Ausgabedateien sehen, eine pro Begrüßung:

```console title="results/2-inputs/"
2-inputs
├── Bonjour-output.txt
├── Hello-output.txt
└── Hola-output.txt
```

Öffne eine der Ausgabedateien und überprüfe, ob sie eine Begrüßung enthält.

Die komprimierte Ausgabe oben zeigt eine einzelne Zusammenfassungszeile für `sayHello`, aber Nextflow hat tatsächlich drei separate Aufgabenausführungen gestartet — eine pro Zeile in der CSV-Datei — und sie parallel ausgeführt, sobald die Ressourcen deines Rechners es erlaubten.

Genau wie die einzelne Aufgabe, die du in [1.3](#13-explore-the-work-directory) erkundet hast, bekommt jede dieser drei Ausführungen ihr eigenes Aufgabenverzeichnis unter `work/`, vollständig von den anderen isoliert:

```console title="work/"
work
├── 2d/276c63.../
│   ├── .command.sh
│   └── Hola-output.txt
├── ab/007682.../
│   ├── .command.sh
│   └── Bonjour-output.txt
└── d9/2476082.../
    ├── .command.sh
    └── Hello-output.txt
```

Jede `.command.sh` enthält immer nur den Befehl für genau diese eine Begrüßung:

```console title="work/2d/276c63.../.command.sh"
#!/bin/bash -ue
echo 'Hola' > 'Hola-output.txt'
```

Diese Isolation macht die parallele Ausführung sicher: Drei gleichzeitig laufende Aufgaben teilen sich niemals ein Arbeitsverzeichnis. Daher kann nichts, was eine Aufgabe schreibt, mit dem kollidieren oder überschrieben werden, was eine andere Aufgabe schreibt — selbst wenn sie zufällig Dateien mit demselben Namen erzeugen.
Das ist auch der Grund, warum `-resume` (als nächstes behandelt) einzelne Aufgaben unabhängig voneinander cachen und wiederverwenden kann: Eingaben, Ausgaben und Logs jeder Aufgabe befinden sich vollständig in ihrem eigenen Verzeichnis, ohne dass etwas zwischen Aufgaben geteilt wird, das aus dem Takt geraten könnte.

### 2.2. Den Workflow erneut mit `-ansi-log false` ausführen

Standardmäßig fasst Nextflow die Ausgabe in einer einzelnen Zusammenfassungszeile pro Prozess zusammen.
Um jeden Prozessaufruf einzeln aufgelistet zu sehen, füge `-ansi-log false` hinzu:

```bash
nextflow run 2-inputs.nf --input data/greetings.csv -ansi-log false
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W  ~  version 26.04.4
    Launching `2-inputs.nf` [extravagant_bardeen] - revision: b90778224d
    [43/0bac1c] Submitted process > sayHello (1)
    [2d/99f604] Submitted process > sayHello (2)
    [6d/7578d7] Submitted process > sayHello (3)

    Outputs:

      /workspaces/training/nextflow-run/results

      first_output:
        - 2-inputs/Hello-output.txt
        - 2-inputs/Bonjour-output.txt
        - 2-inputs/Hola-output.txt
    ```

Dies zeigt alle drei Prozessaufrufe und das eindeutige Arbeitsunterverzeichnis, das für jeden erstellt wurde.

### 2.3. `-resume` verwenden, um abgeschlossene Arbeit zu überspringen

Wechsle nun zur erweiterten Eingabedatei, die zwei weitere Begrüßungen hinzufügt, und füge `-resume` zur Befehlszeile hinzu:

```bash
nextflow run 2-inputs.nf --input data/greetings-extended.csv -resume
```

??? success "Befehlsausgabe"

    ```console hl_lines="6"
    N E X T F L O W   ~  version 26.04.4

    Launching `2-inputs.nf` [adoring_mayer] revision: b90778224d

    executor >  local (2)
    [84/2f3067] sayHello (5) | 5 of 5, cached: 3 ✔

    Outputs:

      /workspaces/training/nextflow-run/results

      first_output:
        - 2-inputs/Bonjour-output.txt
        - 2-inputs/Ciao-output.txt
        - 2-inputs/Hello-output.txt
        - 2-inputs/Hola-output.txt
        - 2-inputs/Ola-output.txt
    ```

Nextflow hat nur die zwei neuen Eingaben verarbeitet.
Die drei Begrüßungen aus dem vorherigen Durchlauf wurden automatisch gecacht und wiederverwendet.

Das funktioniert auch, um die Ausführung von Prozessen zu überspringen, die in einer mehrstufigen Pipeline bereits erfolgreich abgeschlossen wurden.
Zum Beispiel wenn ein Pipeline-Durchlauf durch einen Systemfehler unterbrochen wurde oder wenn du einer Pipeline in der Entwicklung neue Schritte hinzugefügt hast.

Die `-resume`-Funktion ist besonders wertvoll bei langen Pipelines, wo die Wiederherstellung nach einem Fehler kritische Zeit und Ressourcen sparen kann.

### 2.4. Optional: Code-Erklärung

Das Verstehen des Codes ist nicht unbedingt notwendig, wenn du nur Pipelines ausführen möchtest. Wenn du neugierig bist, lohnt sich aber ein Blick.

??? optional "Klicke hier, um den Code zu dieser Übung zu erkunden"

    Die wichtigste Änderung in `2-inputs.nf` befindet sich im `main:`-Abschnitt des Workflows:

    ```groovy title="2-inputs.nf" linenums="14" hl_lines="3 4 5"
        main:
        // Einen Kanal für Eingaben aus einer CSV-Datei erstellen
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        // Eine Begrüßung ausgeben
        sayHello(greeting_ch)
    ```

    Was du hier siehst, nennt sich **channel**: ein Warteschlangen-Konstrukt, das Eingabedaten so verarbeitet, dass Operationen leicht parallelisiert werden können.

    - `channel.fromPath(params.input)` erstellt einen Kanal aus dem mit `--input` angegebenen Dateipfad
    - `.splitCsv()` parst die CSV-Datei in Zeilen
    - `#!groovy .map { line -> line[0] }` extrahiert die erste Spalte aus jeder Zeile

    Das Ergebnis ist ein Kanal, der `Hello`, `Bonjour` und `Hola` enthält.
    Wenn er an `sayHello(greeting_ch)` übergeben wird, ruft Nextflow den Prozess automatisch einmal pro Element auf und führt sie parallel aus, wenn Ressourcen verfügbar sind.

### Fazit

Du weißt, wie du mehrere Eingaben aus einer CSV-Datei parallel verarbeitest und wie du `-resume` verwendest, um abgeschlossene Arbeit nicht zu wiederholen.

### Wie geht es weiter?

Lerne, wie eine vollständige mehrstufige Pipeline Prozesse mithilfe von Kanälen miteinander verknüpft und wie du Container verwendest, um Analyse-Tools und ihre Abhängigkeiten zu verwalten.

---

## 3. Eine mehrstufige Pipeline ausführen

Bisher hast du einen einzelnen Prozess ausgeführt und ihn dann mehrmals parallel über eine Reihe von Eingaben laufen lassen.
Echte Pipelines gehen meist weiter: Sie verknüpfen mehrere Prozesse miteinander, wobei die Ausgabe eines Prozesses zur Eingabe des nächsten wird, und sind oft auf mehr als ein Softwareprogramm angewiesen.
Der Workflow `main.nf` kombiniert beides zu einer vollständigen Pipeline.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_pipeline_complete.svg"
</figure>

Jede Eingabe-Begrüßung durchläuft alle vier Schritte: `sayHello` schreibt sie in eine Datei, `convertToUpper` wandelt den Text in Großbuchstaben um, `collectGreetings` führt alle Ergebnisse in einer Datei zusammen, und `cowpy` erzeugt ASCII-Art aus der zusammengeführten Ausgabe mithilfe eines containerisierten Tools.
Nextflow verbindet diese Schritte mit Kanälen: Die Ausgabe eines Prozesses wird zur Eingabe des nächsten, sodass die gesamte Kette automatisch läuft, sobald Daten verfügbar sind — ohne dass du jeden Schritt manuell koordinieren musst.

Beachte, dass dieser Workflow Module verwendet: Jeder Prozess ist in einer eigenen Datei unter `modules/` definiert, und `main.nf` importiert sie mit `include`-Anweisungen, anstatt sie direkt zu definieren.
Das macht jeden Prozess in mehreren Workflows wiederverwendbar, ohne Code zu duplizieren. Mehr dazu findest du im Code-Erklärungsabschnitt weiter unten.

### 3.1. Den Workflow ausführen

Führe den folgenden Befehl in deinem Terminal aus.

```bash
nextflow run main.nf --input data/greetings.csv
```

Der `character`-Parameter ist in `nextflow.config` standardmäßig auf `turkey` gesetzt, sodass die ASCII-Art einen Truthahn verwendet, sofern du ihn nicht überschreibst (versuche `--character tux` hinzuzufügen).

??? success "Befehlsausgabe"

    ```console hl_lines="6 7 8 9"
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [nostalgic_brahmagupta] revision: ce74f81996

    executor >  local (8)
    [56/8499f6] sayHello (3)       | 3 of 3 ✔
    [cc/0ee42a] convertToUpper (3) | 3 of 3 ✔
    [eb/0f2e24] collectGreetings   | 1 of 1 ✔
    [b5/34e07f] cowpy              | 1 of 1 ✔

    Outputs:

      /workspaces/training/nextflow-run/results

      first_output:
        - batch/intermediates/Bonjour-output.txt
        - batch/intermediates/Hello-output.txt
        - batch/intermediates/Hola-output.txt

      uppercased:
        - batch/intermediates/UPPER-Bonjour-output.txt
        - batch/intermediates/UPPER-Hello-output.txt
        - batch/intermediates/UPPER-Hola-output.txt

      collected: batch/intermediates/COLLECTED-batch-output.txt

      batch_report: batch/batch-report.txt

      cowpy_art: batch/cowpy-COLLECTED-batch-output.txt
    ```

Vier Prozesse wurden ausgeführt, aber nicht gleich oft.
`sayHello` und `convertToUpper` liefen jeweils einmal pro Eingabe (3 of 3): Jede Begrüßung muss einzeln geschrieben und in Großbuchstaben umgewandelt werden.
`collectGreetings` und `cowpy` liefen jeweils nur einmal (1 of 1): Das Zusammenführen der Begrüßungen und das Erzeugen der ASCII-Art ergibt erst Sinn, wenn alle Einzelergebnisse vorliegen.
Diese Fächerform — mehrere parallele Aufgaben, die in eine kleinere Anzahl nachgelagerter Aufgaben münden — ist in echten Pipelines weit verbreitet.

Nextflow wartet nicht darauf, dass ein gesamter Schritt abgeschlossen ist, bevor der nächste beginnt.
Sobald eine `sayHello`-Ausgabe bereit ist, kann die entsprechende `convertToUpper`-Aufgabe starten. Aufgaben verschiedener Prozesse laufen also gleichzeitig, nicht in strikten Batches.
`collectGreetings` und `cowpy` müssen warten, da jeder von ihnen darauf angewiesen ist, dass alle vorgelagerten Ergebnisse verfügbar sind.

Das `results`-Verzeichnis spiegelt dieses Zusammenführen wider, zusammen mit dem, was der Pipeline-Autor zu veröffentlichen und wo zu speichern gewählt hat: Erinnere dich an den `output`-Block aus der Code-Erklärung in 1.4, der diese Struktur definiert.

```console title="results/"
results
└── batch
    ├── batch-report.txt
    ├── cowpy-COLLECTED-batch-output.txt
    └── intermediates
        ├── Bonjour-output.txt
        ├── COLLECTED-batch-output.txt
        ├── Hello-output.txt
        ├── Hola-output.txt
        ├── UPPER-Bonjour-output.txt
        ├── UPPER-Hello-output.txt
        └── UPPER-Hola-output.txt
```

Das oberste Verzeichnis ist nach dem `batch`-Parameter benannt, der standardmäßig `batch` lautet. In späteren Übungen wirst du sehen, wie er sich ändert.

Schau dir `cowpy-COLLECTED-batch-output.txt` für die ASCII-Art-Datei an.

??? abstract "Dateiinhalt"

    ```console title="results/batch/cowpy-COLLECTED-batch-output.txt"
     _________
    / HELLO   \
    | BONJOUR |
    \ HOLA    /
     ---------
      \                                  ,+*^^*+___+++_
       \                           ,*^^^^              )
        \                       _+*                     ^**+_
         \                    +^       _ _++*+_+++_,         )
                  _+^^*+_    (     ,+*^ ^          \+_        )
                 {       )  (    ,(    ,_+--+--,      ^)      ^\
                { (\@)    } f   ,(  ,+-^ __*_*_  ^^\_   ^\       )
               {:;-/    (_+*-+^^^^^+*+*<_ _++_)_    )    )      /
              ( /  (    (        ,___    ^*+_+* )   <    <      \
               U _/     )    *--<  ) ^\-----++__)   )    )       )
                (      )  _(^)^^))  )  )\^^^^^))^*+/    /       /
              (      /  (_))_^)) )  )  ))^^^^^))^^^)__/     +^^
             (     ,/    (^))^))  )  ) ))^^^^^^^))^^)       _)
              *+__+*       (_))^)  ) ) ))^^^^^^))^^^^^)____*^
              \             \_)^)_)) ))^^^^^^^^^^))^^^^)
               (_             ^\__^^^^^^^^^^^^))^^^^^^^)
                 ^\___            ^\__^^^^^^))^^^^^^^^)\\
                      ^^^^^\uuu/^^\uuu/^^^^\^\^\^\^\^\^\^\
                         ___) >____) >___   ^\_\_\_\_\_\_\)
                        ^^^//\\_^^//\\_^       ^(\_\_\_\)
                          ^^^ ^^ ^^^ ^
    ```

Genau wie in [2.1](#21-run-the-workflow) bekommt jede dieser 8 Aufgabenausführungen über alle vier Prozesse hinweg ihr eigenes Verzeichnis unter `work/`, vollständig von den anderen isoliert.
`collectGreetings` veranschaulicht gut, warum das wichtig ist: Es hängt von den Ausgaben aller drei `convertToUpper`-Aufgaben ab, die in drei verschiedenen Aufgabenverzeichnissen liegen. Daher legt Nextflow Symlinks zu diesen Dateien im eigenen Verzeichnis von `collectGreetings` an, anstatt direkt aus den Verzeichnissen der vorgelagerten Aufgaben zu lesen:

```console title="work/eb/0f2e24.../"
COLLECTED-batch-output.txt
UPPER-Bonjour-output.txt -> ../../69/e6c057.../UPPER-Bonjour-output.txt
UPPER-Hello-output.txt   -> ../../cc/ba7d19.../UPPER-Hello-output.txt
UPPER-Hola-output.txt    -> ../../cc/0ee42a.../UPPER-Hola-output.txt
batch-report.txt
.command.sh
```

Jede Aufgabe sieht immer nur die spezifischen Dateien, die sie benötigt — egal woher sie kommen — und niemals den internen Inhalt des Verzeichnisses einer anderen Aufgabe.
Über eine gesamte Pipeline hinweg ist es genau diese Isolation, die du bereits bei einem einzelnen Prozess in [2.1](#21-run-the-workflow) gesehen hast, die es Nextflow ermöglicht, jede Aufgabe aus jedem Prozess gleichzeitig und sicher auszuführen.

!!! note "Hinweis"

    Der `cowpy`-Schritt läuft innerhalb eines Docker-Containers, anstatt auf lokal installierter Software zu basieren.
    Ein Container verpackt eine Anwendung zusammen mit allem, was sie zum Ausführen benötigt. Du musst also keine Abhängigkeiten selbst installieren und verwalten, und die Pipeline verhält sich auf jedem Rechner, der den Container ausführen kann, gleich.
    Nextflow unterstützt auch Conda als Alternative zu Containern. Wie du zwischen beiden wechselst, erfährst du in [Teil 2](./02_configure_pipeline.md).

### 3.2. Optional: Code-Erklärung

Das Verstehen des Codes ist nicht unbedingt notwendig, wenn du nur Pipelines ausführen möchtest. Wenn du neugierig bist, lohnt sich aber ein Blick.

??? optional "Klicke hier, um den Code zu dieser Übung zu erkunden"

    ### Wie Daten von einem Schritt zum nächsten fließen

    Jeder Prozess gibt seinen Ausgabekanal an den nächsten weiter:

    ```groovy title="main.nf" linenums="19" hl_lines="7 8 9"
        main:
        // Einen Kanal für Eingaben aus einer CSV-Datei erstellen
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        sayHello(greeting_ch)
        convertToUpper(sayHello.out)
        collectGreetings(convertToUpper.out.collect(), params.batch)
        cowpy(collectGreetings.out.outfile, params.character)
    ```

    Das Muster `processName.out` verweist auf den Ausgabekanal eines Prozesses.

    Der `.collect()`-Operator sammelt alle einzelnen Ausgaben von `convertToUpper` in einem einzelnen Kanalelement, bevor sie an `collectGreetings` übergeben werden.

    ### Prozessmodule verwenden

    `main.nf` definiert keinen Prozesscode direkt.
    Stattdessen importiert es jeden Prozess aus seiner eigenen Datei unter `modules/`:

    ```groovy title="main.nf" linenums="3"
    include { sayHello } from './modules/sayHello.nf'
    include { convertToUpper } from './modules/convertToUpper.nf'
    include { collectGreetings } from './modules/collectGreetings.nf'
    include { cowpy } from './modules/cowpy.nf'
    ```

    Jede Moduldatei enthält eine einzelne Prozessdefinition, die genauso aufgebaut ist wie das `sayHello`-Modul in [1.4](#14-optional-code-walkthrough).
    Prozesse in separaten Dateien zu halten, macht sie in mehreren Workflows wiederverwendbar, ohne Code zu duplizieren.

    <figure class="excalidraw">
    --8<-- "docs/en/docs/hello_nextflow/img/modules.svg"
    </figure>

    ### Containerisierte Software verwenden

    Der `cowpy`-Prozess läuft innerhalb eines Docker-Containers, der in seiner Moduldatei angegeben ist:

    ```groovy title="modules/cowpy.nf" linenums="2" hl_lines="3"
    process cowpy {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

        input:
        path input_file
        val character

        output:
        path "cowpy-${input_file}"

        script:
        """
        cat ${input_file} | cowpy -c "${character}" > cowpy-${input_file}
        """
    }
    ```

    Nextflow lädt das Image automatisch herunter, führt das Skript innerhalb des Containers aus und räumt danach auf.
    Docker ist für dieses Projekt in `nextflow.config` aktiviert:

    ```groovy title="nextflow.config"
    docker.enabled = true
    ```

    Diese einzelne Zeile aktiviert Docker für jeden Prozess in der Pipeline, der einen Container angegeben hat.

### Fazit

Du hast eine vollständige mehrstufige Pipeline ausgeführt, die mehrere Eingaben parallel mit einem containerisierten Tool verarbeitet.

### Wie geht es weiter?

Weiter zu [Teil 2](./02_configure_pipeline.md), wo du lernst, wie du das Pipeline-Verhalten mit `nextflow.config` konfigurierst.

---

## Zusammenfassung

In diesem Teil hast du gelernt:

- Einen Nextflow-Workflow ausführen und seine Ausgaben finden
- Das `work/`-Verzeichnis und seine Log-Dateien erkunden
- Mehrere Eingaben aus einer CSV-Datei parallel verarbeiten
- `-resume` verwenden, um abgeschlossene Arbeit beim Hinzufügen neuer Eingaben zu überspringen
- Eine mehrstufige Pipeline ausführen, die ein containerisiertes Tool verwendet
