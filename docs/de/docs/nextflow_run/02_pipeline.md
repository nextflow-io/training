# Teil 2: Echte Pipelines ausführen

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

In Teil 1 dieses Kurses (Grundlegende Operationen ausführen) haben wir mit einem Beispiel-Workflow begonnen, der nur minimale Funktionen hatte, um die Code-Komplexität gering zu halten.
Zum Beispiel verwendete `1-hello.nf` einen Kommandozeilenparameter (`--input`), um jeweils einen einzelnen Wert bereitzustellen.

Die meisten realen Pipelines verwenden jedoch ausgefeiltere Funktionen, um eine effiziente Verarbeitung großer Datenmengen im großen Maßstab zu ermöglichen und mehrere Verarbeitungsschritte anzuwenden, die durch manchmal komplexe Logik miteinander verkettet sind.

In diesem Teil des Trainings demonstrieren wir wichtige Funktionen realer Pipelines, indem wir erweiterte Versionen der ursprünglichen Hello World Pipeline ausprobieren.

## 1. Eingabedaten aus einer Datei verarbeiten

In einer realen Pipeline möchten wir typischerweise mehrere Datenpunkte (oder Datenreihen) verarbeiten, die in einer oder mehreren Eingabedateien enthalten sind.
Und wo immer möglich, möchten wir die Verarbeitung unabhängiger Daten parallel ausführen, um die Wartezeit für die Analyse zu verkürzen.

Um zu demonstrieren, wie Nextflow dies macht, haben wir eine CSV-Datei namens `greetings.csv` vorbereitet, die mehrere Eingabebegrüßungen enthält und die Art von spaltenförmigen Daten nachahmt, die du in einer echten Datenanalyse verarbeiten möchtest.
Beachte, dass die Zahlen keine Bedeutung haben, sie sind nur zu Illustrationszwecken da.

```csv title="data/greetings.csv" linenums="1"
Hello,English,123
Bonjour,French,456
Holà,Spanish,789
```

Wir haben auch eine verbesserte Version des ursprünglichen Workflows geschrieben, jetzt `2a-inputs.nf` genannt, die die CSV-Datei einliest, die Begrüßungen extrahiert und jede von ihnen in eine separate Datei schreibt.

<figure class="excalidraw">
--8<-- "docs/en/docs/nextflow_run/img/hello-pipeline-multi-inputs.svg"
</figure>

Lass uns zuerst den Workflow ausführen, und danach schauen wir uns den relevanten Nextflow-Code an.

### 1.1. Den Workflow ausführen

Führe den folgenden Befehl in deinem Terminal aus.

```bash
nextflow run 2a-inputs.nf --input data/greetings.csv
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `2a-inputs.nf` [mighty_sammet] DSL2 - revision: 29fb5352b3

    executor >  local (3)
    [8e/0eb066] sayHello (2) [100%] 3 of 3 ✔
    ```

Erfreulicherweise scheint dies anzuzeigen, dass '3 of 3' Aufrufe für den Prozess gemacht wurden, was ermutigend ist, da es drei Datenzeilen in der CSV gab, die wir als Eingabe bereitgestellt haben.
Dies deutet darauf hin, dass der `sayHello()`-Prozess dreimal aufgerufen wurde, einmal für jede Eingabezeile.

### 1.2. Die veröffentlichten Ausgaben im `results`-Verzeichnis finden

Schauen wir uns das 'results'-Verzeichnis an, um zu sehen, ob unser Workflow immer noch eine Kopie unserer Ausgaben dort schreibt.

??? abstract "Verzeichnisinhalt"

    ```console linenums="1" hl_lines="4-7"
    results
    ├── 1-hello
    |   └── output.txt
    └── 2a-inputs
        ├── Bonjour-output.txt
        ├── Hello-output.txt
        └── Holà-output.txt
    ```

Ja! Wir sehen ein neues Verzeichnis namens `2a-inputs` mit drei Ausgabedateien mit unterschiedlichen Namen, praktischerweise.

Du kannst jede von ihnen öffnen, um dich zu überzeugen, dass sie die entsprechende Begrüßungszeichenkette enthalten.

??? abstract "Dateiinhalte"

    ```console title="results/2a-inputs/Hello-output.txt"
    Hello
    ```

    ```console title="results/2a-inputs/Bonjour-output.txt"
    Bonjour
    ```

    ```console title="results/2a-inputs/Holà-output.txt"
    Holà
    ```

Dies bestätigt, dass jede Begrüßung in der Eingabedatei angemessen verarbeitet wurde.

### 1.3. Die ursprünglichen Ausgaben und Logs finden

Dir ist vielleicht aufgefallen, dass die Konsolenausgabe oben nur auf ein Aufgabenverzeichnis verwies.
Bedeutet das, dass alle drei Aufrufe von `sayHello()` innerhalb dieses einen Aufgabenverzeichnisses ausgeführt wurden?

#### 1.3.1. Das im Terminal angegebene Aufgabenverzeichnis untersuchen

Schauen wir uns das `8e/0eb066`-Aufgabenverzeichnis an.

??? abstract "Verzeichnisinhalt"

    ```console title="8e/0eb066"
    work/8e/0eb066071cdb4123906b7b4ea8b047/
    └── Bonjour-output.txt
    ```

Wir finden nur die Ausgabe, die einer der Begrüßungen entspricht (sowie die Hilfsdateien, wenn wir die Anzeige versteckter Dateien aktivieren).

Was geht hier also vor?

Standardmäßig schreibt das ANSI-Logging-System die Statusinformationen für alle Aufrufe desselben Prozesses in dieselbe Zeile.
Daher zeigte es uns nur einen der drei Aufgabenverzeichnispfade (`8e/0eb066`) in der Konsolenausgabe.
Es gibt zwei weitere, die dort nicht aufgeführt sind.

#### 1.3.2. Das Terminal dazu bringen, mehr Details anzuzeigen

Wir können das Logging-Verhalten ändern, um die vollständige Liste der Prozessaufrufe zu sehen, indem wir `-ansi-log false` zum Befehl wie folgt hinzufügen:

```bash
nextflow run 2a-inputs.nf --input data/greetings.csv -ansi-log false
```

??? success "Befehlsausgabe"

    ```console linenums="1"
    N E X T F L O W  ~  version 25.10.2
    Launching `2a-inputs.nf` [pedantic_hamilton] DSL2 - revision: 6bbc42e49f
    [ab/1a8ece] Submitted process > sayHello (1)
    [0d/2cae24] Submitted process > sayHello (2)
    [b5/0df1d6] Submitted process > sayHello (3)
    ```

Diesmal sehen wir alle drei Prozessausführungen und ihre zugehörigen Work-Unterverzeichnisse in der Ausgabe aufgelistet.
Das Deaktivieren des ANSI-Loggings verhinderte auch, dass Nextflow Farben in der Terminalausgabe verwendete.

Beachte, dass die Art und Weise, wie der Status gemeldet wird, zwischen den beiden Logging-Modi etwas unterschiedlich ist.
Im kondensierten Modus meldet Nextflow, ob Aufrufe erfolgreich abgeschlossen wurden oder nicht.
In diesem erweiterten Modus meldet es nur, dass sie eingereicht wurden.

Dies bestätigt, dass der `sayHello()`-Prozess dreimal aufgerufen wird und für jeden ein separates Aufgabenverzeichnis erstellt wird.

Wenn wir in jedes der dort aufgelisteten Aufgabenverzeichnisse schauen, können wir überprüfen, dass jedes einer der Begrüßungen entspricht.

??? abstract "Verzeichnisinhalt"

    ```console title="ab/1a8ece"
    work/ab/1a8ece307e53f03fce689dde904b64/
    └── Hello-output.txt
    ```

    ```console title="0d/2cae24"
    work/0d/2cae2481a53593bc607077c80c9466/
    └── Bonjour-output.txt
    ```

    ```console title="b5/0df1d6"
    work/b5/0df1d642353269909c2ce23fc2a8fa/
    └── Holà-output.txt
    ```

Dies bestätigt, dass jeder Prozessaufruf isoliert von allen anderen ausgeführt wird.
Das hat viele Vorteile, einschließlich der Vermeidung von Kollisionen, wenn der Prozess Zwischendateien mit nicht eindeutigen Namen erzeugt.

!!! tip

    Für einen komplexen Workflow oder eine große Anzahl von Eingaben kann die vollständige Liste, die im Terminal ausgegeben wird, etwas überwältigend werden, daher verwenden Leute normalerweise `-ansi-log false` nicht im routinemäßigen Gebrauch.

### 1.4. Den Workflow-Code untersuchen

Diese Version des Workflows ist also in der Lage, eine CSV-Datei mit Eingaben einzulesen, die Eingaben separat zu verarbeiten und die Ausgaben eindeutig zu benennen.

Schauen wir uns an, was das im Workflow-Code möglich macht.

??? full-code "Vollständige Code-Datei"

    ```groovy title="2a-inputs.nf" linenums="1" hl_lines="31-33 35"
    #!/usr/bin/env nextflow

    /*
    * Use echo to print 'Hello World!' to a file
    */
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

    /*
    * Pipeline parameters
    */
    params {
        input: Path
    }

    workflow {

        main:
        // Einen Kanal für Eingaben aus einer CSV-Datei erstellen
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        // Eine Begrüßung ausgeben
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }

    output {
        first_output {
            path '2a-inputs'
            mode 'copy'
        }
    }
    ```

Auch hier musst du die Code-Syntax nicht auswendig lernen, aber es ist gut zu lernen, wichtige Komponenten des Workflows zu erkennen, die wichtige Funktionalität bieten.

#### 1.4.1. Die Eingabedaten aus der CSV laden

Dies ist der interessanteste Teil: Wie sind wir vom Nehmen eines einzelnen Werts von der Kommandozeile zum Nehmen einer CSV-Datei, deren Parsen und Verarbeiten der einzelnen darin enthaltenen Begrüßungen gewechselt?

In Nextflow machen wir das mit einem [**channel**](https://nextflow.io/docs/latest/channel.html): einem Warteschlangen-Konstrukt, das entwickelt wurde, um Eingaben effizient zu handhaben und sie von einem Schritt zum anderen in mehrstufigen Workflows zu transportieren, während es eingebaute Parallelität und viele zusätzliche Vorteile bietet.

Lass es uns aufschlüsseln.

```groovy title="2a-inputs.nf" linenums="29" hl_lines="3-5"
    main:
    // Einen Kanal für Eingaben aus einer CSV-Datei erstellen
    greeting_ch = channel.fromPath(params.input)
                        .splitCsv()
                        .map { line -> line[0] }
    // Eine Begrüßung ausgeben
    sayHello(greeting_ch)
```

Dieser Code erstellt einen Kanal namens `greeting_ch`, der die CSV-Datei liest, sie parst und die erste Spalte aus jeder Zeile extrahiert.
Das Ergebnis ist ein Kanal, der `Hello`, `Bonjour` und `Holà` enthält.

??? tip "Wie funktioniert das?"

    Hier ist, was diese Zeile auf Deutsch bedeutet:

    - `channel.fromPath` ist eine **Channel Factory**, die einen Kanal aus Dateipfad(en) erstellt
    - `(params.input)` gibt an, dass der Dateipfad durch `--input` auf der Kommandozeile bereitgestellt wird

    Mit anderen Worten, diese Zeile sagt Nextflow: Nimm den mit `--input` angegebenen Dateipfad und bereite dich darauf vor, seinen Inhalt als Eingabedaten zu behandeln.

    Dann wenden die nächsten beiden Zeilen **Operatoren** an, die das eigentliche Parsen der Datei und Laden der Daten in die entsprechende Datenstruktur durchführen:

    - `.splitCsv()` sagt Nextflow, die CSV-Datei in ein Array zu parsen, das Zeilen und Spalten repräsentiert
    - `.map { line -> line[0] }` sagt Nextflow, nur das Element in der ersten Spalte aus jeder Zeile zu nehmen

    In der Praxis haben wir also, ausgehend von der folgenden CSV-Datei:

    ```csv title="greetings.csv" linenums="1"
    Hello,English,123
    Bonjour,French,456
    Holà,Spanish,789
    ```

    Das in ein Array umgewandelt, das so aussieht:

    ```txt title="Array-Inhalt"
    [[Hello,English,123],[Bonjour,French,456],[Holà,Spanish,789]]
    ```

    Und dann haben wir das erste Element aus jeder der drei Zeilen genommen und sie in einen Nextflow-Kanal geladen, der jetzt enthält: `Hello`, `Bonjour` und `Holà`.

    Wenn du Kanäle und Operatoren im Detail verstehen möchtest, einschließlich wie man sie selbst schreibt, siehe [Hello Nextflow Teil 2: Hello Channels](../hello_nextflow/02_hello_channels.md#4-read-input-values-from-a-csv-file).

#### 1.4.2. Den Prozess für jede Begrüßung aufrufen

Als Nächstes stellen wir in der letzten Zeile des `main:`-Blocks des Workflows den geladenen `greeting_ch`-Kanal als Eingabe für den `sayHello()`-Prozess bereit.

```groovy title="2a-inputs.nf" linenums="29" hl_lines="7"
    main:
    // Einen Kanal für Eingaben aus einer CSV-Datei erstellen
    greeting_ch = channel.fromPath(params.input)
                        .splitCsv()
                        .map { line -> line[0] }
    // Eine Begrüßung ausgeben
    sayHello(greeting_ch)
```

Dies sagt Nextflow, den Prozess individuell auf jedem Element im Kanal auszuführen, _d.h._ auf jeder Begrüßung.
Und weil Nextflow so clever ist, wird es diese Prozessaufrufe parallel ausführen, wenn möglich, abhängig von der verfügbaren Recheninfrastruktur.

So kannst du eine effiziente und skalierbare Verarbeitung vieler Daten (viele Proben oder Datenpunkte, was auch immer deine Forschungseinheit ist) mit vergleichsweise sehr wenig Code erreichen.

#### 1.4.3. Wie die Ausgaben benannt werden

Schließlich lohnt es sich, einen kurzen Blick auf den Prozesscode zu werfen, um zu sehen, wie wir die Ausgabedateien eindeutig benennen lassen.

```groovy title="2a-inputs.nf" linenums="6" hl_lines="7 11"
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

Du siehst, dass im Vergleich zur Version dieses Prozesses in `1-hello.nf` die Ausgabedeklaration und der relevante Teil des Befehls geändert wurden, um den Begrüßungswert in den Ausgabedateinamen aufzunehmen.

Dies ist eine Möglichkeit sicherzustellen, dass die Ausgabedateinamen nicht kollidieren, wenn sie in das gemeinsame Ergebnisverzeichnis veröffentlicht werden.

Und das ist die einzige Änderung, die wir innerhalb der Prozessdeklaration vornehmen mussten!

### Fazit

Du verstehst auf grundlegender Ebene, wie Kanäle und Operatoren es uns ermöglichen, mehrere Eingaben effizient zu verarbeiten.

### Wie geht es weiter?

Entdecke, wie mehrstufige Workflows konstruiert werden und wie sie funktionieren.

---

## 2. Mehrstufige Workflows ausführen

Die meisten realen Workflows umfassen mehr als einen Schritt.
Lass uns auf dem aufbauen, was wir gerade über Kanäle gelernt haben, und schauen, wie Nextflow Kanäle und Operatoren verwendet, um Prozesse in einem mehrstufigen Workflow miteinander zu verbinden.

Zu diesem Zweck stellen wir dir einen Beispiel-Workflow zur Verfügung, der drei separate Schritte miteinander verkettet und Folgendes demonstriert:

1. Datenfluss von einem Prozess zum nächsten
2. Sammeln von Ausgaben aus mehreren Prozessaufrufen in einen einzelnen Prozessaufruf

Konkret haben wir eine erweiterte Version des Workflows namens `2b-multistep.nf` erstellt, die jede Eingabebegrüßung nimmt, sie in Großbuchstaben umwandelt und dann alle großgeschriebenen Begrüßungen in einer einzigen Ausgabedatei sammelt.

<figure class="excalidraw">
--8<-- "docs/en/docs/nextflow_run/img/hello-pipeline-multi-steps.svg"
</figure>

Wie zuvor führen wir zuerst den Workflow aus und schauen uns dann den Code an, um zu sehen, was neu ist.

### 2.1. Den Workflow ausführen

Führe den folgenden Befehl in deinem Terminal aus:

```bash
nextflow run 2b-multistep.nf --input data/greetings.csv
```

??? success "Befehlsausgabe"

    ```console linenums="1"
    N E X T F L O W   ~  version 25.10.2

    Launching `2b-multistep.nf` [soggy_franklin] DSL2 - revision: bc8e1b2726

    [d6/cdf466] sayHello (1)       | 3 of 3 ✔
    [99/79394f] convertToUpper (2) | 3 of 3 ✔
    [1e/83586c] collectGreetings   | 1 of 1 ✔
    ```

Du siehst, dass wie versprochen mehrere Schritte als Teil des Workflows ausgeführt wurden; die ersten beiden (`sayHello` und `convertToUpper`) wurden vermutlich auf jeder einzelnen Begrüßung ausgeführt, und der dritte (`collectGreetings`) wird nur einmal ausgeführt worden sein, auf den Ausgaben aller drei `convertToUpper`-Aufrufe.

### 2.2. Die Ausgaben finden

Lass uns überprüfen, dass das tatsächlich passiert ist, indem wir einen Blick in das `results`-Verzeichnis werfen.

??? abstract "Verzeichnisinhalt"

    ```console linenums="1" hl_lines="8-16"
    results
    ├── 1-hello
    |   └── output.txt
    ├── 2a-inputs
    |   ├── Bonjour-output.txt
    |   ├── Hello-output.txt
    |   └── Holà-output.txt
    └── 2b-multistep
        ├── COLLECTED-batch-output.txt
        ├── batch-report.txt
        └── intermediates
            ├── Bonjour-output.txt
            ├── Hello-output.txt
            ├── Holà-output.txt
            ├── UPPER-Bonjour-output.txt
            ├── UPPER-Hello-output.txt
            └── UPPER-Holà-output.txt

    ```

Wie du sehen kannst, haben wir ein neues Verzeichnis namens `2b-multistep`, und es enthält einige Dateien mehr als zuvor.
Einige der Dateien wurden in ein Unterverzeichnis namens `intermediates` gruppiert, während zwei Dateien auf der obersten Ebene liegen.

Diese beiden sind die Endergebnisse des mehrstufigen Workflows.
Nimm dir eine Minute Zeit, um die Dateinamen anzuschauen und ihre Inhalte zu überprüfen, um zu bestätigen, dass sie das sind, was du erwartest.

??? abstract "Dateiinhalte"

    ```txt title="results/2b-multistep/COLLECTED-batch-output.txt"
    HELLO
    BONJOUR
    HOLà
    ```

    ```txt title="results/2b-multistep/batch-report.txt"
    There were 3 greetings in this batch.
    ```

Die erste enthält unsere drei Begrüßungen, großgeschrieben und wie versprochen in einer einzigen Datei gesammelt.
Die zweite ist eine Berichtsdatei, die einige Informationen über die Ausführung zusammenfasst.

### 2.3. Den Code untersuchen

Schauen wir uns den Code an und identifizieren die wichtigsten Muster für mehrstufige Workflows.

??? full-code "Vollständige Code-Datei"

    ```groovy title="2b-multistep.nf" linenums="1" hl_lines="63 75-78 82-84"
    #!/usr/bin/env nextflow

    /*
    * Use echo to print 'Hello World!' to a file
    */
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

    /*
    * Use a text replacement tool to convert the greeting to uppercase
    */
    process convertToUpper {

        input:
        path input_file

        output:
        path "UPPER-${input_file}"

        script:
        """
        cat '${input_file}' | tr '[a-z]' '[A-Z]' > 'UPPER-${input_file}'
        """
    }

    /*
    * Collect uppercase greetings into a single output file
    */
    process collectGreetings {

        input:
        path input_files
        val batch_name

        output:
        path "COLLECTED-${batch_name}-output.txt", emit: outfile
        path "${batch_name}-report.txt", emit: report

        script:
        count_greetings = input_files.size()
        """
        cat ${input_files} > 'COLLECTED-${batch_name}-output.txt'
        echo 'There were ${count_greetings} greetings in this batch.' > '${batch_name}-report.txt'
        """
    }

    /*
    * Pipeline parameters
    */
    params {
        input: Path
        batch: String = 'batch'
    }

    workflow {

        main:
        // Einen Kanal für Eingaben aus einer CSV-Datei erstellen
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        // Eine Begrüßung ausgeben
        sayHello(greeting_ch)
        // Die Begrüßung in Großbuchstaben umwandeln
        convertToUpper(sayHello.out)
        // Alle Begrüßungen in eine Datei sammeln
        collectGreetings(convertToUpper.out.collect(), params.batch)

        publish:
        first_output = sayHello.out
        uppercased = convertToUpper.out
        collected = collectGreetings.out.outfile
        batch_report = collectGreetings.out.report
    }

    output {
        first_output {
            path '2b-multistep/intermediates'
            mode 'copy'
        }
        uppercased {
            path '2b-multistep/intermediates'
            mode 'copy'
        }
        collected {
            path '2b-multistep'
            mode 'copy'
        }
        batch_report {
            path '2b-multistep'
            mode 'copy'
        }
    }
    ```

Da passiert eine Menge, aber der offensichtlichste Unterschied im Vergleich zur vorherigen Version des Workflows ist, dass es jetzt mehrere Prozessdefinitionen gibt und entsprechend mehrere Prozessaufrufe im Workflow-Block.

Schauen wir genauer hin und sehen, ob wir die interessantesten Teile identifizieren können.

#### 2.3.1. Workflow-Struktur visualisieren

Wenn du VSCode mit der Nextflow-Erweiterung verwendest, kannst du ein hilfreiches Diagramm erhalten, wie die Prozesse verbunden sind, indem du auf den kleinen `DAG preview`-Link klickst, der direkt über dem Workflow-Block in jedem Nextflow-Skript angezeigt wird.

<figure class="excalidraw">
--8<-- "docs/en/docs/nextflow_run/img/DAG-multistep.svg"
</figure>

Dies gibt dir einen schönen Überblick darüber, wie die Prozesse verbunden sind und was sie produzieren.

Du siehst, dass wir zusätzlich zum ursprünglichen `sayHello`-Prozess jetzt auch `convertToUpper` und `collectGreetings` haben, die mit den Namen der Prozesse übereinstimmen, die wir in der Konsolenausgabe gesehen haben.
Die beiden neuen Prozessdefinitionen sind auf die gleiche Weise strukturiert wie der `sayHello`-Prozess, außer dass `collectGreetings` einen zusätzlichen Eingabeparameter namens `batch` nimmt und zwei Ausgaben produziert.

Wir werden nicht im Detail auf den Code für jeden eingehen, aber wenn du neugierig bist, kannst du die Details in [Teil 2 von Hello Nextflow](../hello_nextflow/03_hello_workflow.md) nachschlagen.

Für jetzt lass uns untersuchen, wie die Prozesse miteinander verbunden sind.

#### 2.3.2. Wie die Prozesse verbunden sind

Das wirklich Interessante hier ist, wie die Prozessaufrufe im `main:`-Block des Workflows miteinander verkettet sind.

```groovy title="2b-multistep.nf" linenums="68" hl_lines="9 11"
    main:
    // Einen Kanal für Eingaben aus einer CSV-Datei erstellen
    greeting_ch = channel.fromPath(params.input)
                        .splitCsv()
                        .map { line -> line[0] }
    // Eine Begrüßung ausgeben
    sayHello(greeting_ch)
    // Die Begrüßung in Großbuchstaben umwandeln
    convertToUpper(sayHello.out)
    // Alle Begrüßungen in eine Datei sammeln
    collectGreetings(convertToUpper.out.collect(), params.batch)
```

Du kannst sehen, dass der erste Prozessaufruf, `sayHello(greeting_ch)`, unverändert ist.
Dann verweist der nächste Prozessaufruf, zu `convertToUpper`, auf die Ausgabe von `sayHello` als `sayHello.out`.

Das Muster ist einfach: `processName.out` verweist auf den Ausgabekanal eines Prozesses, der direkt an den nächsten Prozess übergeben werden kann.
So transportieren wir Daten von einem Schritt zum nächsten in Nextflow.

#### 2.3.3. Ein Prozess kann mehrere Eingaben nehmen

Der dritte Prozessaufruf, zu `collectGreetings`, ist etwas anders.

```groovy title="2b-multistep.nf" linenums="77"
    // Alle Begrüßungen in eine Datei sammeln
    collectGreetings(convertToUpper.out.collect(), params.batch)
```

Du siehst, dass diesem Aufruf zwei Eingaben gegeben werden, `convertToUpper.out.collect()` und `params.batch`.
Wenn wir das `.collect()`-Bit vorerst ignorieren, können wir dies als `collectGreetings(input1, input2)` verallgemeinern.

Das passt zu den zwei Eingabedeklarationen im Prozessmodul:

```groovy title="2b-multistep.nf" linenums="40"
process collectGreetings {

    input:
    path input_files
    val batch_name
```

Wenn Nextflow dies parst, wird es die erste Eingabe im Aufruf `path input_files` zuweisen und die zweite `val batch_name`.

Jetzt weißt du also, dass ein Prozess mehrere Eingaben nehmen kann und wie der Aufruf im Workflow-Block aussieht.

Schauen wir uns jetzt diese erste Eingabe genauer an, `convertToUpper.out.collect()`.

#### 2.3.4. Was `collect()` im `collectGreetings`-Aufruf macht

Um die Ausgabe von `sayHello` an `convertToUpper` zu übergeben, haben wir einfach auf den Ausgabekanal von `sayHello` als `sayHello.out` verwiesen. Aber für den nächsten Schritt sehen wir einen Verweis auf `convertToUpper.out.collect()`.

Was ist dieses `collect()`-Bit und was macht es?

Es ist natürlich ein Operator. Genau wie die `splitCsv`- und `map`-Operatoren, denen wir früher begegnet sind.
Diesmal heißt der Operator `collect` und wird auf den von `convertToUpper` erzeugten Ausgabekanal angewendet.

Der `collect`-Operator wird verwendet, um die Ausgaben aus mehreren Aufrufen desselben Prozesses zu sammeln und sie in ein einzelnes Kanalelement zu verpacken.

Im Kontext dieses Workflows nimmt er die drei großgeschriebenen Begrüßungen im `convertToUpper.out`-Kanal (die drei separate Kanalelemente sind und normalerweise in separaten Aufrufen vom nächsten Prozess behandelt würden) und verpackt sie in ein einzelnes Element.
So bekommen wir alle Begrüßungen zurück in dieselbe Datei.

<figure class="excalidraw">
--8<-- "docs/en/docs/nextflow_run/img/with-collect-operator.svg"
</figure>

Im Gegensatz dazu würde Nextflow, wenn wir `collect()` nicht auf die Ausgabe von `convertToUpper()` anwenden würden, bevor wir sie an `collectGreetings()` übergeben, einfach `collectGreetings()` unabhängig auf jeder Begrüßung ausführen, was unser Ziel nicht erreichen würde.

<figure class="excalidraw">
--8<-- "docs/en/docs/nextflow_run/img/without-collect-operator.svg"
</figure>

Es gibt viele andere [Operatoren](https://nextflow.io/docs/latest/reference/operator.html), die verfügbar sind, um Transformationen auf den Inhalt von Kanälen zwischen Prozessaufrufen anzuwenden.

Dies gibt Pipeline-Entwickler\*innen viel Flexibilität für die Anpassung der Flusslogik ihrer Pipeline.
Der Nachteil ist, dass es manchmal schwieriger machen kann zu entschlüsseln, was die Pipeline macht.

#### 2.3.5. Ein Eingabeparameter kann einen Standardwert haben

Dir ist vielleicht aufgefallen, dass `collectGreetings` eine zweite Eingabe nimmt, `params.batch`:

```groovy title="2b-multistep.nf" linenums="77"
    // Alle Begrüßungen in eine Datei sammeln
    collectGreetings(convertToUpper.out.collect(), params.batch)
```

Dies übergibt einen CLI-Parameter namens `--batch` an den Workflow.
Als wir den Workflow jedoch früher gestartet haben, haben wir keinen `--batch`-Parameter angegeben.

Was geht hier vor?
Schau dir den `params`-Block an:

```groovy title="2b-multistep.nf" linenums="61" hl_lines="3"
params {
    input: Path
    batch: String = 'batch'
}
```

Es gibt einen Standardwert, der im Workflow konfiguriert ist, also müssen wir ihn nicht angeben.
Aber wenn wir einen auf der Kommandozeile angeben, wird der von uns angegebene Wert anstelle des Standards verwendet.

Probiere es aus:

```bash
nextflow run 2b-multistep.nf --input data/greetings.csv --batch test
```

??? success "Befehlsausgabe"

    ```console linenums="1"
    N E X T F L O W   ~  version 25.10.2

    Launching `2b-multistep.nf` [soggy_franklin] DSL2 - revision: bc8e1b2726

    [a5/cdff26] sayHello (1)       | 3 of 3 ✔
    [c5/78794f] convertToUpper (2) | 3 of 3 ✔
    [d3/b4d86c] collectGreetings   | 1 of 1 ✔
    ```

Du solltest neue finale Ausgaben sehen, die mit deinem benutzerdefinierten Batch-Namen benannt sind.

??? abstract "Verzeichnisinhalt"

    ```console linenums="1" hl_lines="10 12"
    results
    ├── 1-hello
    |   └── output.txt
    ├── 2a-inputs
    |   ├── Bonjour-output.txt
    |   ├── Hello-output.txt
    |   └── Holà-output.txt
    └── 2b-multistep
        ├── COLLECTED-batch-output.txt
        ├── COLLECTED-test-output.txt
        ├── batch-report.txt
        ├── test-report.txt
        └── intermediates
            ├── Bonjour-output.txt
            ├── Hello-output.txt
            ├── Holà-output.txt
            ├── UPPER-Bonjour-output.txt
            ├── UPPER-Hello-output.txt
            └── UPPER-Holà-output.txt
    ```

Dies ist ein Aspekt der Eingabekonfiguration, den wir in Teil 3 ausführlicher behandeln werden, aber für jetzt ist das Wichtige zu wissen, dass Eingabeparameter Standardwerte erhalten können.

#### 2.3.6. Ein Prozess kann mehrere Ausgaben produzieren

In der `collectGreetings`-Prozessdefinition sehen wir die folgenden Ausgabedeklarationen:

```groovy title="2b-multistep.nf" linenums="46"
    output:
    path "COLLECTED-${batch_name}-output.txt", emit: outfile
    path "${batch_name}-report.txt", emit: report
```

Auf die dann mit dem mit `emit:` angegebenen Namen im `publish:`-Block verwiesen wird:

```groovy title="2b-multistep.nf" linenums="80" hl_lines="4 5"
    publish:
    first_output = sayHello.out
    uppercased = convertToUpper.out
    collected = collectGreetings.out.outfile
    batch_report = collectGreetings.out.report
```

Dies macht es einfach, dann spezifische Ausgaben individuell an andere Prozesse im Workflow zu übergeben, in Kombination mit verschiedenen Operatoren.

#### 2.3.7. Veröffentlichte Ausgaben können organisiert werden

Im `output`-Block haben wir benutzerdefinierte Pfade verwendet, um Zwischenergebnisse zu gruppieren, um es einfacher zu machen, nur die finalen Ausgaben des Workflows herauszusuchen.

```groovy title="2b-multistep.nf" linenums="87" hl_lines="3 7 11 15"
output {
    first_output {
        path '2b-multistep/intermediates'
        mode 'copy'
    }
    uppercased {
        path '2b-multistep/intermediates'
        mode 'copy'
    }
    collected {
        path '2b-multistep'
        mode 'copy'
    }
    batch_report {
        path '2b-multistep'
        mode 'copy'
    }
}
```

Es gibt ausgefeiltere Möglichkeiten, veröffentlichte Ausgaben zu organisieren; wir werden einige davon im Teil über Konfiguration berühren.

!!! tip "Möchtest du mehr über das Erstellen von Workflows erfahren?"

    Für eine detaillierte Behandlung des Aufbaus mehrstufiger Workflows siehe [Hello Nextflow Teil 3: Hello Workflow](../hello_nextflow/03_hello_workflow.md).

### Fazit

Du verstehst auf grundlegender Ebene, wie mehrstufige Workflows mit Kanälen und Operatoren konstruiert werden und wie sie funktionieren.
Du hast auch gesehen, dass Prozesse mehrere Eingaben nehmen und mehrere Ausgaben produzieren können und dass diese auf strukturierte Weise veröffentlicht werden können.

### Wie geht es weiter?

Lerne, wie Nextflow-Pipelines modularisiert werden können, um Code-Wiederverwendung und Wartbarkeit zu fördern.

---

## 3. Modularisierte Pipelines ausführen

Bisher bestanden alle Workflows, die wir uns angeschaut haben, aus einer einzigen Workflow-Datei, die den gesamten relevanten Code enthielt.

Reale Pipelines profitieren jedoch typischerweise davon, _modularisiert_ zu werden, was bedeutet, dass der Code in verschiedene Dateien aufgeteilt wird.
Dies kann ihre Entwicklung und Wartung effizienter und nachhaltiger machen.

Hier werden wir die häufigste Form der Code-Modularität in Nextflow demonstrieren, nämlich die Verwendung von **Modulen**.

In Nextflow ist ein [**module**](https://nextflow.io/docs/latest/module.html) eine einzelne Prozessdefinition, die für sich allein in einer eigenständigen Code-Datei gekapselt ist.
Um ein Modul in einem Workflow zu verwenden, fügst du einfach eine einzeilige Import-Anweisung zu deiner Workflow-Code-Datei hinzu; dann kannst du den Prozess auf die gleiche Weise in den Workflow integrieren, wie du es normalerweise tun würdest.
Das macht es möglich, Prozessdefinitionen in mehreren Workflows wiederzuverwenden, ohne mehrere Kopien des Codes zu erstellen.

Bisher haben wir Workflows ausgeführt, die alle ihre Prozesse in einer monolithischen Code-Datei enthalten hatten.
Jetzt werden wir sehen, wie es aussieht, wenn die Prozesse in einzelnen Modulen gespeichert sind.

Wir haben natürlich wieder einen geeigneten Workflow für Demonstrationszwecke vorbereitet, genannt `2c-modules.nf`, zusammen mit einem Satz von Modulen, die sich im `modules/`-Verzeichnis befinden.

<figure class="excalidraw">
--8<-- "docs/en/docs/nextflow_run/img/modules.svg"
</figure>

??? abstract "Verzeichnisinhalt"

    ```console
    modules/
    ├── collectGreetings.nf
    ├── convertToUpper.nf
    ├── cowpy.nf
    └── sayHello.nf
    ```

Du siehst, es gibt vier Nextflow-Dateien, jede nach einem der Prozesse benannt.
Du kannst die `cowpy.nf`-Datei vorerst ignorieren; darauf kommen wir später.

### 3.1. Den Code untersuchen

Diesmal schauen wir uns zuerst den Code an.
Beginne damit, die `2c-modules.nf`-Workflow-Datei zu öffnen.

??? full-code "Vollständige Code-Datei"

    ```groovy title="2c-modules.nf" linenums="1"
    #!/usr/bin/env nextflow

    // Module einbinden
    include { sayHello } from './modules/sayHello.nf'
    include { convertToUpper } from './modules/convertToUpper.nf'
    include { collectGreetings } from './modules/collectGreetings.nf'

    /*
    * Pipeline parameters
    */
    params {
        input: Path
        batch: String = 'batch'
    }

    workflow {

        main:
        // Einen Kanal für Eingaben aus einer CSV-Datei erstellen
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        // Eine Begrüßung ausgeben
        sayHello(greeting_ch)
        // Die Begrüßung in Großbuchstaben umwandeln
        convertToUpper(sayHello.out)
        // Alle Begrüßungen in eine Datei sammeln
        collectGreetings(convertToUpper.out.collect(), params.batch)

        publish:
        first_output = sayHello.out
        uppercased = convertToUpper.out
        collected = collectGreetings.out.outfile
        batch_report = collectGreetings.out.report
    }

    output {
        first_output {
            path '2c-modules/intermediates'
            mode 'copy'
        }
        uppercased {
            path '2c-modules/intermediates'
            mode 'copy'
        }
        collected {
            path '2c-modules'
            mode 'copy'
        }
        batch_report {
            path '2c-modules'
            mode 'copy'
        }
    }
    ```

Du siehst, dass die Workflow-Logik genau die gleiche ist wie in der vorherigen Version des Workflows.
Der Prozesscode ist jedoch aus der Workflow-Datei verschwunden, und stattdessen gibt es `include`-Anweisungen, die auf separate Dateien unter `modules` verweisen.

```groovy title="hello-modules.nf" linenums="3"
// Module einbinden
include { sayHello } from './modules/sayHello.nf'
include { convertToUpper } from './modules/convertToUpper.nf'
include { collectGreetings } from './modules/collectGreetings.nf'
```

Öffne eine dieser Dateien und du wirst den Code für den entsprechenden Prozess finden.

??? full-code "Vollständige Code-Datei"

    ```groovy title="modules/sayHello.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
    * Use echo to print 'Hello World!' to a file
    */
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

Wie du sehen kannst, hat sich der Prozesscode nicht geändert; er wurde nur in eine individuelle Moduldatei kopiert, anstatt in der Haupt-Workflow-Datei zu sein.
Das Gleiche gilt für die anderen beiden Prozesse.

Schauen wir uns also an, wie es aussieht, diese neue Version auszuführen.

### 3.2. Den Workflow ausführen

Führe diesen Befehl in deinem Terminal aus, mit dem `-resume`-Flag:

```bash
nextflow run 2c-modules.nf --input data/greetings.csv -resume
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `2c-modules.nf` [soggy_franklin] DSL2 - revision: bc8e1b2726

    [d6/cdf466] sayHello (1)       | 3 of 3, cached: 3 ✔
    [99/79394f] convertToUpper (2) | 3 of 3, cached: 3 ✔
    [1e/83586c] collectGreetings   | 1 of 1, cached: 1 ✔
    ```

Du wirst bemerken, dass alle Prozessausführungen erfolgreich gecacht wurden, was bedeutet, dass Nextflow erkannt hat, dass es die angeforderte Arbeit bereits erledigt hat, obwohl der Code aufgeteilt und die Haupt-Workflow-Datei umbenannt wurde.

Nichts davon ist für Nextflow wichtig; was zählt, ist das Job-Skript, das generiert wird, sobald der gesamte Code zusammengezogen und ausgewertet wurde.

!!! tip

    Es ist auch möglich, einen Abschnitt eines Workflows als 'Subworkflow' zu kapseln, der in eine größere Pipeline importiert werden kann, aber das liegt außerhalb des Umfangs dieses Kurses.

    Du kannst mehr über die Entwicklung zusammensetzbarer Workflows in der Side Quest zu [Workflows of Workflows](https://training.nextflow.io/latest/side_quests/workflows_of_workflows/) erfahren.

### Fazit

Du weißt, wie Prozesse in eigenständigen Modulen gespeichert werden können, um Code-Wiederverwendung zu fördern und die Wartbarkeit zu verbessern.

### Wie geht es weiter?

Lerne, Container für die Verwaltung von Software-Abhängigkeiten zu verwenden.

---

## 4. Containerisierte Software verwenden

Bisher mussten die Workflows, die wir als Beispiele verwendet haben, nur sehr grundlegende Textverarbeitungsoperationen mit UNIX-Tools ausführen, die in unserer Umgebung verfügbar sind.

Reale Pipelines erfordern jedoch typischerweise spezialisierte Tools und Pakete, die nicht standardmäßig in den meisten Umgebungen enthalten sind.
Normalerweise müsstest du diese Tools installieren, ihre Abhängigkeiten verwalten und eventuelle Konflikte lösen.

Das ist alles sehr mühsam und nervig.
Ein viel besserer Weg, dieses Problem anzugehen, ist die Verwendung von **Containern**.

Ein **Container** ist eine leichtgewichtige, eigenständige, ausführbare Softwareeinheit, die aus einem Container-**Image** erstellt wird und alles enthält, was zum Ausführen einer Anwendung benötigt wird, einschließlich Code, Systembibliotheken und Einstellungen.

!!! Tip

    Wir lehren dies mit der Technologie [Docker](https://www.docker.com/get-started/), aber Nextflow unterstützt auch mehrere andere Container-Technologien.
    Du kannst mehr über die Nextflow-Unterstützung für Container [hier](https://nextflow.io/docs/latest/container.html) erfahren.

### 4.1. Einen Container direkt verwenden

Lass uns zuerst versuchen, direkt mit einem Container zu interagieren.
Dies wird helfen, dein Verständnis davon zu festigen, was Container sind, bevor wir anfangen, sie in Nextflow zu verwenden.

#### 4.1.1. Das Container-Image pullen

Um einen Container zu verwenden, lädst du normalerweise ein Container-Image aus einer Container-Registry herunter oder "pullst" es, und führst dann das Container-Image aus, um eine Container-Instanz zu erstellen.

Die allgemeine Syntax ist wie folgt:

```bash title="Syntax"
docker pull '<container>'
```

- `docker pull` ist die Anweisung an das Container-System, ein Container-Image aus einem Repository zu pullen.
- `'<container>'` ist die URI-Adresse des Container-Images.

Als Beispiel pullen wir ein Container-Image, das [cowpy](https://github.com/jeffbuttars/cowpy) enthält, eine Python-Implementierung eines Tools namens `cowsay`, das ASCII-Art generiert, um beliebige Texteingaben auf lustige Weise anzuzeigen.

Es gibt verschiedene Repositories, in denen du veröffentlichte Container finden kannst.
Wir haben den [Seqera Containers](https://seqera.io/containers/)-Service verwendet, um dieses Docker-Container-Image aus dem `cowpy`-Conda-Paket zu generieren: `'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'`.

Führe den vollständigen Pull-Befehl aus:

```bash
docker pull 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
```

??? success "Befehlsausgabe"

    ```console
    Unable to find image 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273' locally
    131d6a1b707a8e65: Pulling from library/cowpy
    dafa2b0c44d2: Pull complete
    dec6b097362e: Pull complete
    f88da01cff0b: Pull complete
    4f4fb700ef54: Pull complete
    92dc97a3ef36: Pull complete
    403f74b0f85e: Pull complete
    10b8c00c10a5: Pull complete
    17dc7ea432cc: Pull complete
    bb36d6c3110d: Pull complete
    0ea1a16bbe82: Pull complete
    030a47592a0a: Pull complete
    622dd7f15040: Pull complete
    895fb5d0f4df: Pull complete
    Digest: sha256:fa50498b32534d83e0a89bb21fec0c47cc03933ac95c6b6587df82aaa9d68db3
    Status: Downloaded newer image for community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273
    community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273
    ```

Dies sagt dem System, das angegebene Image herunterzuladen.
Sobald der Download abgeschlossen ist, hast du eine lokale Kopie des Container-Images.

#### 4.1.2. Den Container starten

Container können als einmaliger Befehl ausgeführt werden, aber du kannst sie auch interaktiv verwenden, was dir eine Shell-Eingabeaufforderung innerhalb des Containers gibt und es dir ermöglicht, mit dem Befehl zu spielen.

Die allgemeine Syntax ist wie folgt:

```bash title="Syntax"
docker run --rm '<container>' [tool command]
```

- `docker run --rm '<container>'` ist die Anweisung an das Container-System, eine Container-Instanz aus einem Container-Image zu starten und einen Befehl darin auszuführen.
- `--rm` sagt dem System, die Container-Instanz herunterzufahren, nachdem der Befehl abgeschlossen wurde.

Vollständig zusammengesetzt sieht der Container-Ausführungsbefehl so aus:

```bash
docker run --rm -it 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
```

Führe diesen Befehl aus, und du solltest sehen, dass sich deine Eingabeaufforderung zu etwas wie `(base) root@b645838b3314:/tmp#` ändert, was anzeigt, dass du jetzt innerhalb des Containers bist.

Du kannst dies überprüfen, indem du `ls` ausführst, um Verzeichnisinhalte aufzulisten:

```bash
ls /
```

??? success "Befehlsausgabe"

    ```console
    bin  boot  dev  etc  home  lib  lib64  media  mnt  opt  proc  root  run  sbin  srv  sys  tmp  usr  var
    ```

Du siehst, dass das Dateisystem innerhalb des Containers sich vom Dateisystem auf deinem Host-System unterscheidet.

!!! Tip

    Wenn du einen Container ausführst, ist er standardmäßig vom Host-System isoliert.
    Das bedeutet, dass der Container nicht auf Dateien auf dem Host-System zugreifen kann, es sei denn, du erlaubst es ausdrücklich, indem du angibst, dass du ein Volume als Teil des `docker run`-Befehls mounten möchtest, mit der folgenden Syntax:

    ```bash title="Syntax"
    -v <outside_path>:<inside_path>
    ```

    Dies etabliert effektiv einen Tunnel durch die Container-Wand, den du verwenden kannst, um auf diesen Teil deines Dateisystems zuzugreifen.

    Dies wird ausführlicher in [Teil 5 von Hello Nextflow](../hello_nextflow/05_hello_containers.md) behandelt.

#### 4.1.3. Das `cowpy`-Tool ausführen

Von innerhalb des Containers kannst du den `cowpy`-Befehl direkt ausführen.

```bash
cowpy "Hello Containers"
```

??? success "Befehlsausgabe"

    ```console
    ______________________________________________________
    < Hello Containers >
    ------------------------------------------------------
        \   ^__^
          \  (oo)\_______
            (__)\       )\/\
              ||----w |
              ||     ||
    ```

Dies erzeugt ASCII-Art des Standard-Kuh-Charakters (oder 'Cowacter') mit einer Sprechblase, die den von uns angegebenen Text enthält.

Jetzt, da du die grundlegende Verwendung getestet hast, kannst du versuchen, ihm einige Parameter zu geben.
Zum Beispiel sagt die Tool-Dokumentation, dass wir den Charakter mit `-c` setzen können.

```bash
cowpy "Hello Containers" -c tux
```

??? success "Befehlsausgabe"

    ```console
    __________________
    < Hello Containers >
    ------------------
      \
        \
            .--.
          |o_o |
          |:_/ |
          //   \ \
        (|     | )
        /'\_   _/`\
        \___)=(___/
    ```

Diesmal zeigt die ASCII-Art-Ausgabe den Linux-Pinguin Tux, weil wir den `-c tux`-Parameter angegeben haben.

Da du innerhalb des Containers bist, kannst du den cowpy-Befehl so oft ausführen, wie du möchtest, und die Eingabeparameter variieren, ohne dir Sorgen machen zu müssen, irgendwelche Bibliotheken auf deinem System selbst zu installieren.

??? tip "Andere verfügbare Charaktere"

    Verwende das '-c'-Flag, um einen anderen Charakter auszuwählen, einschließlich:

    `beavis`, `cheese`, `daemon`, `dragonandcow`, `ghostbusters`, `kitty`, `moose`, `milk`, `stegosaurus`, `turkey`, `turtle`, `tux`

Spiele gerne damit herum.
Wenn du fertig bist, verlasse den Container mit dem `exit`-Befehl:

```bash
exit
```

Du wirst dich wieder in deiner normalen Shell befinden.

### 4.2. Einen Container in einem Workflow verwenden

Wenn wir eine Pipeline ausführen, möchten wir in der Lage sein, Nextflow mitzuteilen, welchen Container es bei jedem Schritt verwenden soll, und wichtig ist, dass wir möchten, dass es all die Arbeit erledigt, die wir gerade gemacht haben: den Container pullen, ihn starten, den Befehl ausführen und den Container herunterfahren, wenn er fertig ist.

Gute Nachrichten: Das ist genau das, was Nextflow für uns tun wird.
Wir müssen nur einen Container für jeden Prozess angeben.

Um zu demonstrieren, wie dies funktioniert, haben wir eine weitere Version unseres Workflows erstellt, die `cowpy` auf der Datei der gesammelten Begrüßungen ausführt, die im dritten Schritt produziert wurde.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-pipeline-cowpy.svg"
</figure>

Dies sollte eine Datei ausgeben, die die ASCII-Art mit den drei Begrüßungen in der Sprechblase enthält.

#### 4.2.1. Den Code untersuchen

Der Workflow ist dem vorherigen sehr ähnlich, plus dem zusätzlichen Schritt, um `cowpy` auszuführen.

??? full-code "Vollständige Code-Datei"

    ```groovy title="2d-container.nf" linenums="1" hl_lines="7 15 32 39 59-62"
    #!/usr/bin/env nextflow

    // Module einbinden
    include { sayHello } from './modules/sayHello.nf'
    include { convertToUpper } from './modules/convertToUpper.nf'
    include { collectGreetings } from './modules/collectGreetings.nf'
    include { cowpy } from './modules/cowpy.nf'

    /*
    * Pipeline parameters
    */
    params {
        input: Path
        batch: String = 'batch'
        character: String
    }

    workflow {

        main:
        // Einen Kanal für Eingaben aus einer CSV-Datei erstellen
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        // Eine Begrüßung ausgeben
        sayHello(greeting_ch)
        // Die Begrüßung in Großbuchstaben umwandeln
        convertToUpper(sayHello.out)
        // Alle Begrüßungen in eine Datei sammeln
        collectGreetings(convertToUpper.out.collect(), params.batch)
        // ASCII-Art der Begrüßungen mit cowpy generieren
        cowpy(collectGreetings.out.outfile, params.character)

        publish:
        first_output = sayHello.out
        uppercased = convertToUpper.out
        collected = collectGreetings.out.outfile
        batch_report = collectGreetings.out.report
        cowpy_art = cowpy.out
    }

    output {
        first_output {
            path '2d-container/intermediates'
            mode 'copy'
        }
        uppercased {
            path '2d-container/intermediates'
            mode 'copy'
        }
        collected {
            path '2d-container/intermediates'
            mode 'copy'
        }
        batch_report {
            path '2d-container'
            mode 'copy'
        }
        cowpy_art {
            path '2d-container'
            mode 'copy'
        }
    }
    ```

Du siehst, dass dieser Workflow einen `cowpy`-Prozess aus einer Moduldatei importiert und ihn auf die Ausgabe des `collectGreetings()`-Aufrufs plus einen Eingabeparameter namens `params.character` aufruft.

```groovy title="2d-container.nf" linenums="31"
// ASCII-Art der Begrüßungen mit cowpy generieren
cowpy(collectGreetings.out.outfile, params.character)
```

Der `cowpy`-Prozess, der den cowpy-Befehl umschließt, um ASCII-Art zu generieren, ist im `cowpy.nf`-Modul definiert.

??? full-code "Vollständige Code-Datei"

    ```groovy title="modules/cowpy.nf" linenums="1"
    #!/usr/bin/env nextflow

    // ASCII-Art mit cowpy generieren (https://github.com/jeffbuttars/cowpy)
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

Der `cowpy`-Prozess benötigt zwei Eingaben: den Pfad zu einer Eingabedatei, die den Text enthält, der in die Sprechblase soll (`input_file`), und einen Wert für die Charaktervariable.

Wichtig ist, dass er auch die Zeile `container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'` enthält, die auf die Container-URI verweist, die wir früher verwendet haben.

#### 4.2.2. Überprüfen, dass Docker in der Konfiguration aktiviert ist

Wir werden Teil 3 dieses Trainingskurses leicht vorwegnehmen, indem wir die `nextflow.config`-Konfigurationsdatei einführen, die eine der Hauptmöglichkeiten ist, die Nextflow für die Konfiguration der Workflow-Ausführung bietet.
Wenn eine Datei namens `nextflow.config` im aktuellen Verzeichnis vorhanden ist, wird Nextflow sie automatisch laden und jede darin enthaltene Konfiguration anwenden.

Zu diesem Zweck haben wir eine `nextflow.config`-Datei mit einer einzigen Codezeile beigefügt, die Docker aktiviert.

```groovy title="nextflow.config" linenums="1"
docker.enabled = true
```

Diese Konfiguration sagt Nextflow, Docker für jeden Prozess zu verwenden, der einen kompatiblen Container angibt.

!!! tip

    Es ist technisch möglich, die Docker-Ausführung von der Kommandozeile aus zu aktivieren, auf Basis einer einzelnen Ausführung, mit dem `-with-docker <container>`-Parameter.
    Das erlaubt uns jedoch nur, einen Container für den gesamten Workflow anzugeben, während der Ansatz, den wir dir gerade gezeigt haben, es uns ermöglicht, einen anderen Container pro Prozess anzugeben.
    Letzteres ist viel besser für Modularität, Code-Wartung und Reproduzierbarkeit.

#### 4.2.3. Den Workflow ausführen

Nur zur Wiederholung, das ist, was wir gleich ausführen werden:

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_pipeline_complete.svg"
</figure>

Denkst du, es wird funktionieren?

Lass uns den Workflow mit dem `-resume`-Flag ausführen und angeben, dass wir möchten, dass der Charakter der Truthahn ist.

```bash
nextflow run 2d-container.nf --input data/greetings.csv --character turkey -resume
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `2d-container.nf` [elegant_brattain] DSL2 - revision: 028a841db1

    executor >  local (1)
    [95/fa0bac] sayHello (3)       | 3 of 3, cached: 3 ✔
    [92/32533f] convertToUpper (3) | 3 of 3, cached: 3 ✔
    [aa/e697a2] collectGreetings   | 1 of 1, cached: 1 ✔
    [7f/caf718] cowpy              | 1 of 1 ✔
    ```

Die ersten drei Schritte wurden gecacht, da wir sie bereits zuvor ausgeführt haben, aber der `cowpy`-Prozess ist neu, also wird dieser tatsächlich ausgeführt.

Du kannst die Ausgabe des `cowpy`-Schritts im `results`-Verzeichnis finden.

??? abstract "Dateiinhalte"

    ```console title="results/2d-container/cowpy-COLLECTED-batch-output.txt"
    _________
    / HOLà    \
    | HELLO   |
    \ BONJOUR /
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

Du siehst, dass der Charakter alle Begrüßungen sagt, da er auf der Datei der gesammelten großgeschriebenen Begrüßungen lief.

Wichtiger noch, wir konnten dies als Teil unserer Pipeline ausführen, ohne eine ordnungsgemäße Installation von cowpy und all seinen Abhängigkeiten durchführen zu müssen.
Und wir können die Pipeline jetzt mit Mitarbeiter\*innen teilen und sie auf ihrer Infrastruktur ausführen lassen, ohne dass sie auch etwas installieren müssen, abgesehen von Docker oder einer seiner Alternativen (wie Singularity/Apptainer), wie oben erwähnt.

#### 4.2.4. Untersuchen, wie Nextflow die containerisierte Aufgabe gestartet hat

Als abschließende Coda zu diesem Abschnitt werfen wir einen Blick in das Work-Unterverzeichnis für einen der `cowpy`-Prozessaufrufe, um etwas mehr Einblick zu bekommen, wie Nextflow unter der Haube mit Containern arbeitet.

Überprüfe die Ausgabe deines `nextflow run`-Befehls, um den Pfad zum Work-Unterverzeichnis für den `cowpy`-Prozess zu finden.
Wenn wir uns ansehen, was wir für die oben gezeigte Ausführung bekommen haben, beginnt die Konsolen-Log-Zeile für den `cowpy`-Prozess mit `[7f/caf718]`.
Das entspricht dem folgenden abgekürzten Verzeichnispfad: `work/7f/caf718`.

In diesem Verzeichnis findest du die `.command.run`-Datei, die alle Befehle enthält, die Nextflow in deinem Namen im Verlauf der Ausführung der Pipeline ausgeführt hat.

??? abstract "Dateiinhalte"

    ```console title="work/7f/caf71890cce1667c094d880f4b6dcc/.command.run"
    #!/bin/bash
    ### ---
    ### name: 'cowpy'
    ### container: 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
    ### outputs:
    ### - 'cowpy-COLLECTED-batch-output.txt'
    ### ...
    set -e
    set -u
    NXF_DEBUG=${NXF_DEBUG:=0}; [[ $NXF_DEBUG > 1 ]] && set -x
    NXF_ENTRY=${1:-nxf_main}


    nxf_sleep() {
      sleep $1 2>/dev/null || sleep 1;
    }

    nxf_date() {
        local ts=$(date +%s%3N);
        if [[ ${#ts} == 10 ]]; then echo ${ts}000
        elif [[ $ts == *%3N ]]; then echo ${ts/\%3N/000}
        elif [[ $ts == *3N ]]; then echo ${ts/3N/000}
        elif [[ ${#ts} == 13 ]]; then echo $ts
        else echo "Unexpected timestamp value: $ts"; exit 1
        fi
    }

    nxf_env() {
        echo '============= task environment ============='
        env | sort | sed "s/\(.*\)AWS\(.*\)=\(.\{6\}\).*/\1AWS\2=\3xxxxxxxxxxxxx/"
        echo '============= task output =================='
    }

    nxf_kill() {
        declare -a children
        while read P PP;do
            children[$PP]+=" $P"
        done < <(ps -e -o pid= -o ppid=)

        kill_all() {
            [[ $1 != $$ ]] && kill $1 2>/dev/null || true
            for i in ${children[$1]:=}; do kill_all $i; done
        }

        kill_all $1
    }

    nxf_mktemp() {
        local base=${1:-/tmp}
        mkdir -p "$base"
        if [[ $(uname) = Darwin ]]; then mktemp -d $base/nxf.XXXXXXXXXX
        else TMPDIR="$base" mktemp -d -t nxf.XXXXXXXXXX
        fi
    }

    nxf_fs_copy() {
      local source=$1
      local target=$2
      local basedir=$(dirname $1)
      mkdir -p $target/$basedir
      cp -fRL $source $target/$basedir
    }

    nxf_fs_move() {
      local source=$1
      local target=$2
      local basedir=$(dirname $1)
      mkdir -p $target/$basedir
      mv -f $source $target/$basedir
    }

    nxf_fs_rsync() {
      rsync -rRl $1 $2
    }

    nxf_fs_rclone() {
      rclone copyto $1 $2/$1
    }

    nxf_fs_fcp() {
      fcp $1 $2/$1
    }

    on_exit() {
        local last_err=$?
        local exit_status=${nxf_main_ret:=0}
        [[ ${exit_status} -eq 0 && ${nxf_unstage_ret:=0} -ne 0 ]] && exit_status=${nxf_unstage_ret:=0}
        [[ ${exit_status} -eq 0 && ${last_err} -ne 0 ]] && exit_status=${last_err}
        printf -- $exit_status > /workspaces/training/nextflow-run/work/7f/caf71890cce1667c094d880f4b6dcc/.exitcode
        set +u
        docker rm $NXF_BOXID &>/dev/null || true
        exit $exit_status
    }

    on_term() {
        set +e
        docker stop $NXF_BOXID
    }

    nxf_launch() {
        docker run -i --cpu-shares 1024 -e "NXF_TASK_WORKDIR" -v /workspaces/training/nextflow-run/work:/workspaces/training/nextflow-run/work -w "$NXF_TASK_WORKDIR" --name $NXF_BOXID community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273 /bin/bash -ue /workspaces/training/nextflow-run/work/7f/caf71890cce1667c094d880f4b6dcc/.command.sh
    }

    nxf_stage() {
        true
        # stage input files
        rm -f COLLECTED-batch-output.txt
        ln -s /workspaces/training/nextflow-run/work/7f/f435e3f2cf95979b5f3d7647ae6696/COLLECTED-batch-output.txt COLLECTED-batch-output.txt
    }

    nxf_unstage_outputs() {
        true
    }

    nxf_unstage_controls() {
        true
    }

    nxf_unstage() {
        if [[ ${nxf_main_ret:=0} == 0 ]]; then
            (set -e -o pipefail; (nxf_unstage_outputs | tee -a .command.out) 3>&1 1>&2 2>&3 | tee -a .command.err)
            nxf_unstage_ret=$?
        fi
        nxf_unstage_controls
    }

    nxf_main() {
        trap on_exit EXIT
        trap on_term TERM INT USR2
        trap '' USR1

        [[ "${NXF_CHDIR:-}" ]] && cd "$NXF_CHDIR"
        export NXF_BOXID="nxf-$(dd bs=18 count=1 if=/dev/urandom 2>/dev/null | base64 | tr +/ 0A | tr -d '\r\n')"
        NXF_SCRATCH=''
        [[ $NXF_DEBUG > 0 ]] && nxf_env
        touch /workspaces/training/nextflow-run/work/7f/caf71890cce1667c094d880f4b6dcc/.command.begin
        set +u
        set -u
        [[ $NXF_SCRATCH ]] && cd $NXF_SCRATCH
        export NXF_TASK_WORKDIR="$PWD"
        nxf_stage

        set +e
        (set -o pipefail; (nxf_launch | tee .command.out) 3>&1 1>&2 2>&3 | tee .command.err) &
        pid=$!
        wait $pid || nxf_main_ret=$?
        nxf_unstage
    }

    $NXF_ENTRY
    ```

Wenn du in dieser Datei nach `nxf_launch` suchst, solltest du so etwas sehen:

```console
nxf_launch() {
    docker run -i --cpu-shares 1024 -e "NXF_TASK_WORKDIR" -v /workspaces/training/nextflow-run/work:/workspaces/training/nextflow-run/work -w "$NXF_TASK_WORKDIR" --name $NXF_BOXID community.wave.seqera.io/library/pip_cowpy:131d6a1b707a8e65 /bin/bash -ue /workspaces/training/nextflow-run/work/7f/caf7189fca6c56ba627b75749edcb3/.command.sh
}
```

Dieser Launch-Befehl zeigt, dass Nextflow einen sehr ähnlichen `docker run`-Befehl verwendet, um den Prozessaufruf zu starten, wie wir es getan haben, als wir ihn manuell ausgeführt haben.
Er mountet auch das entsprechende Work-Unterverzeichnis in den Container, setzt das Arbeitsverzeichnis innerhalb des Containers entsprechend und führt unser vorlagenbasiertes Bash-Skript in der `.command.sh`-Datei aus.

Dies bestätigt, dass all die harte Arbeit, die wir im vorherigen Abschnitt manuell machen mussten, jetzt von Nextflow für uns erledigt wird!

### Fazit

Du verstehst, welche Rolle Container bei der Verwaltung von Software-Tool-Versionen und der Sicherstellung der Reproduzierbarkeit spielen.

Allgemeiner hast du ein grundlegendes Verständnis davon, was die Kernkomponenten realer Nextflow-Pipelines sind und wie sie organisiert sind.
Du kennst die Grundlagen, wie Nextflow mehrere Eingaben effizient verarbeiten kann, Workflows ausführt, die aus mehreren miteinander verbundenen Schritten bestehen, modulare Code-Komponenten nutzt und Container für größere Reproduzierbarkeit und Portabilität verwendet.

### Wie geht es weiter?

Mach eine weitere Pause! Das war ein großer Haufen Informationen darüber, wie Nextflow-Pipelines funktionieren.

Im letzten Abschnitt dieses Trainings werden wir tiefer in das Thema Konfiguration eintauchen.
Du wirst lernen, wie du die Ausführung deiner Pipeline an deine Infrastruktur anpasst sowie die Konfiguration von Eingaben und Parametern verwaltest.

---

## Quiz

<quiz>
Warum erstellt Nextflow ein separates Aufgabenverzeichnis für jeden Prozessaufruf?
- [ ] Um die Ausführungsgeschwindigkeit zu verbessern
- [ ] Um die Speichernutzung zu reduzieren
- [x] Um Ausführungen zu isolieren und Kollisionen zwischen Ausgaben zu vermeiden
- [ ] Um parallele Dateikompression zu ermöglichen

Mehr erfahren: [1.3. Die ursprünglichen Ausgaben und Logs finden](#13-die-ursprunglichen-ausgaben-und-logs-finden)
</quiz>

<quiz>
Was macht die `-ansi-log false`-Option beim Ausführen eines Workflows?
- [ ] Deaktiviert alle Konsolenausgaben
- [x] Entfernt Farbe aus der Ausgabe
- [x] Zeigt alle Aufgabenverzeichnispfade an, anstatt sie in einer Zeile zu kondensieren
- [ ] Aktiviert den ausführlichen Debug-Modus

Mehr erfahren: [1.3.2. Das Terminal dazu bringen, mehr Details anzuzeigen](#132-das-terminal-dazu-bringen-mehr-details-anzuzeigen)

Du kannst auch eine der folgenden Umgebungsvariablen verwenden, wenn du diesen Stil bevorzugst:

```bash
export NXF_ANSI_LOG=0
# oder
export NO_COLOR=1
```

</quiz>

<quiz>
Im Code `#!groovy channel.fromPath(params.input).splitCsv().map { line -> line[0] }`, was macht `#!groovy .map { line -> line[0] }`?
- [ ] Filtert leere Zeilen heraus
- [ ] Sortiert die Zeilen alphabetisch
- [x] Extrahiert die erste Spalte aus jeder CSV-Zeile
- [ ] Zählt die Anzahl der Zeilen

Mehr erfahren: [1.4.1. Die Eingabedaten aus der CSV laden](#141-die-eingabedaten-aus-der-csv-laden)
</quiz>

<quiz>
Warum ist es wichtig, den Eingabewert in Ausgabedateinamen einzuschließen (z.B. `#!groovy "${greeting}-output.txt"`)?
- [ ] Um die Verarbeitungsgeschwindigkeit zu verbessern
- [ ] Um die Resume-Funktionalität zu ermöglichen
- [x] Um zu verhindern, dass Ausgabedateien sich gegenseitig überschreiben, wenn mehrere Eingaben verarbeitet werden
- [ ] Um Dateien einfacher zu komprimieren

Mehr erfahren: [1.4.3. Wie die Ausgaben benannt werden](#143-wie-die-ausgaben-benannt-werden)
</quiz>

<quiz>
Was ist der Zweck der `include`-Anweisung in einem modularisierten Workflow?
- [ ] Prozesscode in die Workflow-Datei zu kopieren
- [x] Eine Prozessdefinition aus einer externen Moduldatei zu importieren
- [ ] Konfigurationseinstellungen einzuschließen
- [ ] Dokumentationskommentare hinzuzufügen

Mehr erfahren: [3. Modularisierte Pipelines ausführen](#3-modularisierte-pipelines-ausfuhren)
</quiz>

<quiz>
Wenn du einen Workflow modularisierst und ihn mit `-resume` ausführst, was passiert?
- [ ] Caching ist für modulare Prozesse deaktiviert
- [ ] Alle Aufgaben müssen erneut ausgeführt werden
- [x] Caching funktioniert normal basierend auf den generierten Job-Skripten
- [ ] Nur die Haupt-Workflow-Datei wird gecacht

Mehr erfahren: [3.2. Den Workflow ausführen](#32-den-workflow-ausfuhren)
</quiz>

<quiz>
Was gibt die `container`-Direktive in einer Prozessdefinition an?
- [ ] Das Arbeitsverzeichnis für den Prozess
- [ ] Die maximale Speicherzuweisung
- [x] Die Container-Image-URI, die für die Ausführung des Prozesses verwendet werden soll
- [ ] Das Ausgabedateiformat

Mehr erfahren: [4.2. Einen Container in einem Workflow verwenden](#42-einen-container-in-einem-workflow-verwenden)
</quiz>

<quiz>
Was enthält die `nxf_launch`-Funktion in der `.command.run`-Datei?
- [ ] Die Nextflow-Versionsinformationen
- [ ] Die Workflow-Parameter
- [x] Den `docker run`-Befehl mit Volume-Mounts und Container-Einstellungen
- [ ] Die Prozess-Eingabedeklarationen

Mehr erfahren: [4.2.4. Untersuchen, wie Nextflow die containerisierte Aufgabe gestartet hat](#424-untersuchen-wie-nextflow-die-containerisierte-aufgabe-gestartet-hat)
</quiz>

<quiz>
Was handhabt Nextflow automatisch beim Ausführen eines containerisierten Prozesses? (Wähle alle zutreffenden aus)
- [x] Das Container-Image bei Bedarf pullen
- [x] Das Work-Verzeichnis in den Container mounten
- [x] Das Prozess-Skript innerhalb des Containers ausführen
- [x] Die Container-Instanz nach der Ausführung aufräumen

Mehr erfahren: [4. Containerisierte Software verwenden](#4-containerisierte-software-verwenden)
</quiz>
