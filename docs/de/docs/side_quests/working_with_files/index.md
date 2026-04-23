# Verarbeitung von Datei-Eingaben

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Wissenschaftliche Analyse-Workflows beinhalten oft die Verarbeitung großer Mengen von Dateien.
Nextflow bietet leistungsstarke Werkzeuge, um Dateien effizient zu verwalten und dir zu helfen, deine Daten mit minimalem Code zu organisieren und zu verarbeiten.

### Lernziele

In dieser Side Quest erkunden wir, wie Nextflow mit Dateien umgeht – von grundlegenden Dateioperationen bis hin zu fortgeschritteneren Techniken für die Arbeit mit Dateisammlungen.
Du lernst, wie du Metadaten aus Dateinamen extrahierst, was eine häufige Anforderung in wissenschaftlichen Analyse-Pipelines ist.

Am Ende dieser Side Quest kannst du:

- Path-Objekte aus Dateipfad-Strings mit Nextflows `file()`-Methode erstellen
- Dateiattribute wie Name, Erweiterung und übergeordnetes Verzeichnis abrufen
- Lokale und entfernte Dateien transparent über URIs verarbeiten
- Kanäle zur Automatisierung der Dateiverarbeitung mit `channel.fromPath()` und `channel.fromFilePairs()` verwenden
- Metadaten aus Dateinamen mithilfe von String-Manipulation extrahieren und strukturieren
- Zusammengehörige Dateien mit Pattern-Matching und Glob-Ausdrücken gruppieren
- Dateioperationen mit korrekter Eingabeverarbeitung in Nextflow-Prozesse integrieren
- Prozessausgaben mithilfe metadatengesteuerter Verzeichnisstrukturen organisieren

Diese Fähigkeiten helfen dir, Workflows zu erstellen, die verschiedene Arten von Datei-Eingaben mit großer Flexibilität verarbeiten können.

### Voraussetzungen

Bevor du diese Side Quest in Angriff nimmst, solltest du:

- Das Tutorial [Hello Nextflow](../../hello_nextflow/) oder einen gleichwertigen Einsteigerkurs abgeschlossen haben.
- Mit grundlegenden Nextflow-Konzepten und -Mechanismen (Prozesse, Kanäle, Operatoren) vertraut sein.

<!-- I removed the suggestion to do the metamaps SQ first because that works more naturally after -->

---

## 0. Erste Schritte

#### Trainings-Codespace öffnen

Falls noch nicht geschehen, öffne die Trainingsumgebung wie in der [Umgebung einrichten](../envsetup/index.md) beschrieben.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### In das Projektverzeichnis wechseln

Wechseln wir in das Verzeichnis, in dem sich die Dateien für dieses Tutorial befinden.

```bash
cd side-quests/working_with_files
```

Du kannst VSCode so einstellen, dass es sich auf dieses Verzeichnis fokussiert:

```bash
code .
```

#### Die Materialien ansehen

Du findest eine einfache Workflow-Datei namens `main.nf`, ein `modules`-Verzeichnis mit zwei Moduldateien und ein `data`-Verzeichnis mit einigen Beispieldatendateien.

??? abstract "Verzeichnisinhalt"

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

Dieses Verzeichnis enthält Paired-End-Sequenzierungsdaten von drei Patient\*innen (A, B, C).

Für jede\*n Patient\*in haben wir Proben vom Typ `tumor` (typischerweise aus Tumorbiopsien stammend) oder `normal` (aus gesundem Gewebe oder Blut entnommen).
Falls du nicht mit der Krebsanalyse vertraut bist: Dieses Modell verwendet gepaarte Tumor/Normal-Proben für kontrastive Analysen.

Für Patient A haben wir speziell zwei Sätze technischer Replikate (Wiederholungen).

Die Sequenzierungsdatendateien sind nach der üblichen `_R1_`- und `_R2_`-Konvention für sogenannte 'Forward Reads' und 'Reverse Reads' benannt.

_Keine Sorge, wenn du mit diesem experimentellen Design nicht vertraut bist – es ist nicht entscheidend für das Verständnis dieses Tutorials._

#### Die Aufgabe verstehen

Deine Aufgabe ist es, einen Nextflow-Workflow zu schreiben, der:

1. **Eingabedateien lädt** mithilfe von Nextflows Dateiverarbeitungsmethoden
2. **Metadaten extrahiert** (Patienten-ID, Replikat, Probentyp) aus der Dateinamenstruktur
3. **Gepaarte Dateien (R1/R2) gruppiert** mithilfe von `channel.fromFilePairs()`
4. **Die Dateien verarbeitet** mit einem bereitgestellten Analysemodul
5. **Ausgaben organisiert** in einer Verzeichnisstruktur basierend auf den extrahierten Metadaten

#### Bereitschafts-Checkliste

Bereit zum Eintauchen?

- [ ] Ich verstehe das Ziel dieses Kurses und seine Voraussetzungen
- [ ] Mein Codespace läuft
- [ ] Ich habe mein Arbeitsverzeichnis entsprechend eingestellt
- [ ] Ich verstehe die Aufgabe

Wenn du alle Punkte abhaken kannst, kann es losgehen.

---

## 1. Grundlegende Dateioperationen

### 1.1. Den Typ eines Objekts mit `.class` bestimmen

Schau dir die Workflow-Datei `main.nf` an:

```groovy title="main.nf" linenums="1"
#!/usr/bin/env nextflow

workflow {

    // Ein Path-Objekt aus einem String-Pfad erstellen
    myFile = 'data/patientA_rep1_normal_R1_001.fastq.gz'

    println "${myFile} is of class ${myFile.class}"
}
```

Dies ist ein Mini-Workflow (ohne Prozesse), der auf einen einzelnen Dateipfad verweist, ihn dann zusammen mit seiner Klasse auf der Konsole ausgibt.

??? info "Was ist `.class`?"

    In Nextflow sagt uns `.class`, mit welchem Objekttyp wir arbeiten. Es ist wie die Frage „Was für ein Ding ist das?" – um herauszufinden, ob es ein String, eine Zahl, eine Datei oder etwas anderes ist.
    Das hilft uns, den Unterschied zwischen einem einfachen String und einem Path-Objekt in den nächsten Abschnitten zu veranschaulichen.

Führen wir den Workflow aus:

```bash
nextflow run main.nf
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [romantic_chandrasekhar] DSL2 - revision: 5a4a89bc3a

    data/patientA_rep1_normal_R1_001.fastq.gz is of class java.lang.String
    ```

Wie du siehst, hat Nextflow den String-Pfad genau so ausgegeben, wie wir ihn geschrieben haben.

Das ist nur Textausgabe; Nextflow hat noch nichts Besonderes damit gemacht.
Wir haben auch bestätigt, dass es sich aus Nextflows Sicht nur um einen String handelt (der Klasse `java.lang.String`).
Das macht Sinn, da wir Nextflow noch nicht mitgeteilt haben, dass es einer Datei entspricht.

### 1.2. Ein Path-Objekt mit file() erstellen

Wir können Nextflow mitteilen, wie es mit Dateien umgehen soll, indem wir [Path-Objekte](https://www.nextflow.io/docs/latest/reference/stdlib-types.html#path) aus Pfad-Strings erstellen.

In unserem Workflow können wir den String-Pfad `data/patientA_rep1_normal_R1_001.fastq.gz` mit der `file()`-Methode in ein Path-Objekt umwandeln, das Zugriff auf Dateieigenschaften und -operationen bietet.

Bearbeite `main.nf` und umschließe den String mit `file()` wie folgt:

=== "Danach"

    ```groovy title="main.nf" linenums="5" hl_lines="2"
        // Ein Path-Objekt aus einem String-Pfad erstellen
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        println "${myFile} is of class ${myFile.class}"
    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="5" hl_lines="2"
        // Ein Path-Objekt aus einem String-Pfad erstellen
        myFile = 'data/patientA_rep1_normal_R1_001.fastq.gz'

        println "${myFile} is of class ${myFile.class}"
    ```

Führe den Workflow erneut aus:

```bash
nextflow run main.nf
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [kickass_coulomb] DSL2 - revision: 5af44b1b59

    /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz is of class class sun.nio.fs.UnixPath
    ```

Diesmal siehst du den vollständigen absoluten Pfad anstelle des relativen Pfads, den wir als Eingabe angegeben haben.

Nextflow hat unseren String in ein Path-Objekt umgewandelt und es zum tatsächlichen Dateispeicherort auf dem System aufgelöst.
Der Dateipfad ist jetzt absolut, wie z. B. `/workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz`.

Beachte auch, dass die Path-Objekt-Klasse `sun.nio.fs.UnixPath` ist: So stellt Nextflow lokale Dateien dar.
Wie wir später sehen werden, haben entfernte Dateien andere Klassennamen (z. B. `nextflow.file.http.XPath` für HTTP-Dateien), aber sie funktionieren alle auf genau die gleiche Weise und können in deinen Workflows identisch verwendet werden.

!!! tip "Tipp"

    **Der wesentliche Unterschied:**

    - **Path-String**: Nur Text, den Nextflow als Zeichen behandelt
    - **Path-Objekt**: Eine intelligente Dateireferenz, mit der Nextflow arbeiten kann

    Stell es dir so vor: Ein Path-String ist wie eine Adresse auf Papier, während ein Path-Objekt wie eine in ein Navigationsgerät geladene Adresse ist, das weiß, wie man dorthin navigiert, und dir Details über die Route mitteilen kann.

### 1.3. Dateiattribute abrufen

Warum ist das hilfreich? Nun, da Nextflow versteht, dass `myFile` ein Path-Objekt und kein einfacher String ist, können wir auf die verschiedenen Attribute des Path-Objekts zugreifen.

Aktualisieren wir unseren Workflow, um die eingebauten Dateiattribute auszugeben:

=== "Danach"

    ```groovy title="main.nf" linenums="5" hl_lines="4-9"
        // Ein Path-Objekt aus einem String-Pfad erstellen
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Dateiattribute ausgeben
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="5" hl_lines="4"
        // Ein Path-Objekt aus einem String-Pfad erstellen
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        println "${myFile} is of class ${myFile.class}"
    ```

Workflow ausführen:

```bash
nextflow run main.nf
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [ecstatic_ampere] DSL2 - revision: f3fa3dcb48

    File object class: sun.nio.fs.UnixPath
    File name: patientA_rep1_normal_R1_001.fastq.gz
    Simple name: patientA_rep1_normal_R1_001
    Extension: gz
    Parent directory: /workspaces/training/side-quests/working_with_files/data
    ```

Du siehst die verschiedenen Dateiattribute, die oben auf der Konsole ausgegeben werden.

### 1.4. Die Datei an einen Prozess übergeben

Der Unterschied zwischen Strings und Path-Objekten wird entscheidend, wenn du anfängst, echte Workflows mit Prozessen zu erstellen.
Bisher haben wir bestätigt, dass Nextflow unsere Eingabedatei jetzt als Datei behandelt, aber lass uns sehen, ob wir tatsächlich etwas mit dieser Datei in einem Prozess ausführen können.

#### 1.4.1. Den Prozess importieren und den Code untersuchen

Wir stellen dir ein vorgefertigtes Prozessmodul namens `COUNT_LINES` zur Verfügung, das eine Dateieingabe entgegennimmt und zählt, wie viele Zeilen darin enthalten sind.

Um den Prozess im Workflow zu verwenden, musst du nur eine include-Anweisung vor dem Workflow-Block hinzufügen:

=== "Danach"

    ```groovy title="main.nf" linenums="1" hl_lines="3"
    #!/usr/bin/env nextflow

    include { COUNT_LINES } from './modules/count_lines.nf'

    workflow {
    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="1"
    #!/usr/bin/env nextflow

    workflow {
    ```

Du kannst die Moduldatei öffnen, um ihren Code zu untersuchen:

```groovy title="modules/count_lines.nf" linenums="1"
#!/usr/bin/env nextflow

process COUNT_LINES {
    debug true

    input:
    path input_file

    script:
    """
    set -o pipefail
    echo "Processing file: $input_file"
    gzip -dc $input_file | wc -l
    """
}
```

Wie du siehst, ist es ein recht einfaches Skript, das die Datei entpackt und zählt, wie viele Zeilen sie enthält.

??? info "Was macht `debug true`?"

    Die `debug true`-Direktive in der Prozessdefinition bewirkt, dass Nextflow die Ausgabe deines Skripts (wie die Zeilenanzahl „40") direkt im Ausführungsprotokoll ausgibt.
    Ohne diese Direktive würdest du nur den Prozessausführungsstatus sehen, aber nicht die tatsächliche Ausgabe deines Skripts.

    Weitere Informationen zum Debuggen von Nextflow-Prozessen findest du in der Side Quest [Debugging Nextflow Workflows](debugging.md).

#### 1.4.2. Einen Aufruf von `COUNT_LINES` hinzufügen

Da der Prozess jetzt im Workflow verfügbar ist, können wir einen Aufruf des `COUNT_LINES`-Prozesses hinzufügen, um ihn auf der Eingabedatei auszuführen.

Nimm folgende Änderungen am Workflow vor:

=== "Danach"

    ```groovy title="main.nf" linenums="7" hl_lines="11-12"
        // Ein Path-Objekt aus einem String-Pfad erstellen
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Dateiattribute ausgeben
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"

        // Zeilen in der Datei zählen
        COUNT_LINES(myFile)
    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="7" hl_lines="4-9"
        // Ein Path-Objekt aus einem String-Pfad erstellen
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Dateiattribute ausgeben
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    ```

Und jetzt den Workflow ausführen:

```bash
nextflow run main.nf
```

??? success "Befehlsausgabe"

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

Das zeigt, dass wir die Datei innerhalb eines Prozesses korrekt verarbeiten können.

Konkret hat Nextflow folgende Operationen erfolgreich durchgeführt:

- Die Datei in das work directory gestaged
- Die .gz-Datei dekomprimiert
- Die Zeilen gezählt (in diesem Fall 40 Zeilen)
- Ohne Fehler abgeschlossen

Der Schlüssel zu diesem reibungslosen Ablauf ist, dass wir Nextflow explizit mitteilen, dass unsere Eingabe eine Datei ist und als solche behandelt werden soll.

### 1.5. Grundlegende Datei-Eingabefehler beheben

Das ist eine häufige Stolperfalle für Nextflow-Einsteiger\*innen, also nehmen wir uns ein paar Minuten Zeit, um zu sehen, was passiert, wenn man es falsch macht.

Es gibt zwei Hauptstellen, an denen die Dateiverarbeitung schiefgehen kann: auf der Ebene des Workflows und auf der Ebene des Prozesses.

#### 1.5.1. Fehler auf Workflow-Ebene

Schauen wir uns an, was passiert, wenn wir die Datei wieder als String behandeln, wenn wir die Eingabe im Workflow-Block angeben.

Nimm folgende Änderungen am Workflow vor und stelle sicher, dass du die pfadspezifischen Print-Anweisungen auskommentierst:

=== "Danach"

    ```groovy title="main.nf" linenums="7" hl_lines="2 6-11"
        // Ein Path-Objekt aus einem String-Pfad erstellen
        myFile = 'data/patientA_rep1_normal_R1_001.fastq.gz'

        // Dateiattribute ausgeben
        println "File object class: ${myFile.class}"
        /*
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
        */

        // Zeilen in der Datei zählen
        COUNT_LINES(myFile)
    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="7" hl_lines="4-9"
        // Ein Path-Objekt aus einem String-Pfad erstellen
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Dateiattribute ausgeben
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"

        // Zeilen in der Datei zählen
        COUNT_LINES(myFile)
    ```

Und jetzt den Workflow ausführen:

```bash
nextflow run main.nf
```

??? failure "Befehlsausgabe"

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

Das ist der wichtige Teil:

```console
Not a valid path value: 'data/patientA_rep1_normal_R1_001.fastq.gz'
```

Wenn du eine `path`-Eingabe angibst, überprüft Nextflow, ob du tatsächliche Dateireferenzen übergibst und keine einfachen Strings.
Dieser Fehler teilt dir mit, dass `'data/patientA_rep1_normal_R1_001.fastq.gz'` kein gültiger Pfadwert ist, weil es ein String und kein Path-Objekt ist.

Nextflow hat das Problem sofort erkannt und gestoppt, bevor der Prozess überhaupt gestartet wurde.

#### 1.5.2. Fehler auf Prozess-Ebene

Die andere Stelle, an der wir vergessen könnten anzugeben, dass Nextflow die Eingabe als Datei behandeln soll, ist in der Prozessdefinition.

!!! warning "Warnung"

    Damit dieser Test korrekt funktioniert, behalte den Workflow in seinem fehlerhaften Zustand (mit einem einfachen String statt `file()`).
    In Kombination mit `val` im Prozess erzeugt das den unten gezeigten Fehler.

Nimm folgende Änderung am Modul vor:

=== "Danach"

    ```groovy title="modules/count_lines.nf" linenums="3" hl_lines="5"
    process COUNT_LINES {
        debug true

        input:
        val input_file
    ```

=== "Vorher"

    ```groovy title="modules/count_lines.nf" linenums="3" hl_lines="5"
    process COUNT_LINES {
        debug true

        input:
        path input_file
    ```

Und jetzt den Workflow erneut ausführen:

```bash
nextflow run main.nf
```

??? failure "Befehlsausgabe"

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

Das zeigt viele Details über den Fehler, da der Prozess so eingestellt ist, dass er Debug-Informationen ausgibt, wie oben erwähnt.

Die relevantesten Abschnitte sind:

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

Das bedeutet, dass das System die Datei nicht finden konnte; wenn du jedoch den Pfad nachschaust, gibt es dort eine Datei mit diesem Namen.

Als wir das ausgeführt haben, hat Nextflow den String-Wert an das Skript weitergegeben, aber die eigentliche Datei nicht in das work directory _gestaged_.
Der Prozess hat also versucht, den relativen String `data/patientA_rep1_normal_R1_001.fastq.gz` zu verwenden, aber diese Datei existiert nicht im work directory des Prozesses.

Zusammengenommen zeigen diese beiden Beispiele, wie wichtig es ist, Nextflow mitzuteilen, wenn eine Eingabe als Datei behandelt werden soll.

!!! note "Hinweis"

    Stelle sicher, dass du beide absichtlichen Fehler behebst, bevor du mit dem nächsten Abschnitt fortfährst.

### Fazit

- Path-Strings vs. Path-Objekte: Strings sind nur Text, Path-Objekte sind intelligente Dateireferenzen
- Die `file()`-Methode wandelt einen String-Pfad in ein Path-Objekt um, mit dem Nextflow arbeiten kann
- Du kannst auf Dateieigenschaften wie `name`, `simpleName`, `extension` und `parent` [mithilfe von Dateiattributen](https://www.nextflow.io/docs/latest/working-with-files.html#getting-file-attributes) zugreifen
- Die Verwendung von Path-Objekten statt Strings ermöglicht es Nextflow, Dateien in deinem Workflow korrekt zu verwalten
- Prozess-Eingabeergebnisse: Eine korrekte Dateiverarbeitung erfordert Path-Objekte, keine Strings, damit Dateien korrekt gestaged und für Prozesse zugänglich sind.

---

## 2. Entfernte Dateien verwenden

Eine der wichtigsten Funktionen von Nextflow ist die Möglichkeit, nahtlos zwischen lokalen Dateien (auf demselben Rechner) und entfernten Dateien, die über das Internet zugänglich sind, zu wechseln.

Wenn du es richtig machst, solltest du die Logik deines Workflows nie ändern müssen, um Dateien aus verschiedenen Speicherorten zu verarbeiten.
Alles, was du tun musst, um eine entfernte Datei zu verwenden, ist das entsprechende Präfix im Dateipfad anzugeben, wenn du ihn dem Workflow übergibst.

Zum Beispiel hat `/path/to/data` kein Präfix, was darauf hinweist, dass es ein „normaler" lokaler Dateipfad ist, während `s3://path/to/data` das Präfix `s3://` enthält, was darauf hinweist, dass es sich in Amazons S3-Objektspeicher befindet.

Viele verschiedene Protokolle werden unterstützt:

- HTTP(S)/FTP (http://, https://, ftp://)
- Amazon S3 (s3://)
- Azure Blob Storage (az://)
- Google Cloud Storage (gs://)

Um eines davon zu verwenden, gib einfach das entsprechende Präfix im String an, der dann technisch als Uniform Resource Identifier (URI) statt als Dateipfad bezeichnet wird.
Nextflow übernimmt die Authentifizierung und das Staging der Dateien an den richtigen Ort, einschließlich aller anderen Dateioperationen, die du erwarten würdest.

Die wichtigste Stärke dieses Systems ist, dass wir zwischen Umgebungen wechseln können, ohne die Pipeline-Logik zu ändern.
Du kannst zum Beispiel mit einem kleinen, lokalen Testdatensatz entwickeln und dann einfach durch Ändern der URI auf einen vollständigen Testdatensatz in einem entfernten Speicher umsteigen.

### 2.1. Eine Datei aus dem Internet verwenden

Testen wir das, indem wir den lokalen Pfad, den wir unserem Workflow übergeben, durch einen HTTPS-Pfad ersetzen, der auf eine Kopie derselben Daten verweist, die auf Github gespeichert ist.

!!! warning "Warnung"

    Das funktioniert nur, wenn du eine aktive Internetverbindung hast.

Öffne `main.nf` erneut und ändere den Eingabepfad wie folgt:

=== "Danach"

    ```groovy title="main.nf" linenums="2" hl_lines="2"
        // Eine entfernte Datei aus dem Internet verwenden
        myFile = file('https://raw.github.com/nextflow-io/training/master/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz')

        // Dateiattribute ausgeben
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="2" hl_lines="2"
        // Ein Path-Objekt aus einem String-Pfad erstellen
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Dateiattribute ausgeben
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    ```

Führen wir den Workflow aus:

```bash
nextflow run main.nf
```

??? success "Befehlsausgabe"

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

Es funktioniert! Du siehst, dass sich sehr wenig geändert hat.

Der einzige Unterschied in der Konsolenausgabe ist, dass die Path-Objekt-Klasse jetzt `nextflow.file.http.XPath` ist, während sie für den lokalen Pfad `sun.nio.fs.UnixPath` war.
Du musst dir diese Klassen nicht merken; wir erwähnen das nur, um zu zeigen, dass Nextflow die verschiedenen Speicherorte entsprechend identifiziert und verarbeitet.

Im Hintergrund hat Nextflow die Datei in ein Staging-Verzeichnis innerhalb des work directory heruntergeladen.
Diese gestagede Datei kann dann als lokale Datei behandelt und per Symlink in das entsprechende Prozessverzeichnis eingebunden werden.

Du kannst überprüfen, ob das passiert ist, indem du dir den Inhalt des work directory ansiehst, das sich am Hash-Wert des Prozesses befindet.

??? abstract "Work-Verzeichnis Inhalt"

    Wenn der Prozess-Hash `8a/2ab7ca` war, könntest du das work directory erkunden:

    ```console
    $ ls -la work/8a/2ab7ca*/
    total 16
    drwxr-xr-x  6 user  staff   192 Jan 28 10:00 .
    drwxr-xr-x  3 user  staff    96 Jan 28 10:00 ..
    -rw-r--r--  1 user  staff     0 Jan 28 10:00 .command.begin
    -rw-r--r--  1 user  staff   127 Jan 28 10:00 .command.sh
    lrwxr-xr-x  1 user  staff    89 Jan 28 10:00 patientA_rep1_normal_R1_001.fastq.gz -> /path/to/work/stage/.../patientA_rep1_normal_R1_001.fastq.gz
    ```

    Der Symlink verweist auf eine gestagede Kopie der entfernten Datei, die Nextflow automatisch heruntergeladen hat.

Beachte, dass der Download-Schritt bei größeren Dateien im Vergleich zur Ausführung mit lokalen Dateien etwas mehr Zeit in Anspruch nimmt.
Nextflow prüft jedoch, ob bereits eine gestagede Kopie vorhanden ist, um unnötige Downloads zu vermeiden.
Wenn du also erneut auf derselben Datei ausführst und die gestagede Datei nicht gelöscht hast, verwendet Nextflow die gestagede Kopie.

Das zeigt, wie einfach es ist, mit Nextflow zwischen lokalen und entfernten Daten zu wechseln – eine wichtige Funktion von Nextflow.

!!! note "Hinweis"

    Die eine wichtige Ausnahme von diesem Prinzip ist, dass du keine Glob-Muster oder Verzeichnispfade mit HTTPS verwenden kannst, da HTTPS keine mehreren Dateien auflisten kann. Du musst also genaue Datei-URLs angeben.
    Andere Speicherprotokolle wie Blob-Speicher (`s3://`, `az://`, `gs://`) können jedoch sowohl Globs als auch Verzeichnispfade verwenden.

    So könntest du Glob-Muster mit Cloud-Speicher verwenden:

    ```groovy title="Cloud storage examples (not runnable in this environment)"
    // S3 mit Glob-Mustern - würde mehrere Dateien matchen
    ch_s3_files = channel.fromPath('s3://my-bucket/data/*.fastq.gz')

    // Azure Blob Storage mit Glob-Mustern
    ch_azure_files = channel.fromPath('az://container/data/patient*_R{1,2}.fastq.gz')

    // Google Cloud Storage mit Glob-Mustern
    ch_gcs_files = channel.fromPath('gs://bucket/data/sample_*.fastq.gz')
    ```

    Im nächsten Abschnitt zeigen wir dir, wie du in der Praxis mit Globs arbeitest.

### 2.2. Zurück zur lokalen Datei wechseln

Für den Rest dieser Side Quest werden wir wieder unsere lokalen Beispieldateien verwenden, also wechseln wir die Workflow-Eingabe zurück zur ursprünglichen Datei:

=== "Danach"

    ```groovy title="main.nf" linenums="2" hl_lines="2"
        // Ein Path-Objekt aus einem String-Pfad erstellen
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Dateiattribute ausgeben
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="2" hl_lines="2"
        // Ein Path-Objekt aus einem String-Pfad erstellen
        myFile = file('https://raw.github.com/nextflow-io/training/master/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz')

        // Dateiattribute ausgeben
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    ```

### Fazit

- Auf entfernte Daten wird über eine URI zugegriffen (HTTP, FTP, S3, Azure, Google Cloud)
- Nextflow lädt die Daten automatisch herunter und stagt sie an den richtigen Ort, solange diese Pfade an Prozesse übergeben werden
- Schreibe keine Logik zum Herunterladen oder Hochladen entfernter Dateien!
- Lokale und entfernte Dateien erzeugen unterschiedliche Objekttypen, funktionieren aber identisch
- **Wichtig**: HTTP/HTTPS funktioniert nur mit einzelnen Dateien (keine Glob-Muster)
- Cloud-Speicher (S3, Azure, GCS) unterstützt sowohl einzelne Dateien als auch Glob-Muster
- Du kannst nahtlos zwischen lokalen und entfernten Datenquellen wechseln, ohne die Code-Logik zu ändern (solange das Protokoll die erforderlichen Operationen unterstützt)

---

## 3. Die `fromPath()`-Kanal-Factory verwenden

Bisher haben wir immer mit einer einzelnen Datei gearbeitet, aber in Nextflow werden wir typischerweise einen Eingabekanal mit mehreren Eingabedateien zur Verarbeitung erstellen wollen.

Eine naive Methode wäre, die `file()`-Methode mit [`channel.of()`](https://www.nextflow.io/docs/latest/reference/channel.html#of) zu kombinieren:

```groovy title="Syntax example"
ch_files = channel.of([file('data/patientA_rep1_normal_R1_001.fastq.gz')],
                      [file('data/patientA_rep1_normal_R1_001.fastq.gz')])
```

Das funktioniert, ist aber umständlich.

!!! tip "Tipp: Wann `file()` vs. `channel.fromPath()` verwenden"

    - Verwende `file()`, wenn du ein einzelnes Path-Objekt für die direkte Manipulation benötigst (prüfen, ob eine Datei existiert, ihre Attribute lesen oder an einen einzelnen Prozessaufruf übergeben)
    - Verwende `channel.fromPath()`, wenn du einen Kanal benötigst, der mehrere Dateien aufnehmen kann, insbesondere mit Glob-Mustern, oder wenn Dateien durch mehrere Prozesse fließen sollen

Hier kommt [`channel.fromPath()`](https://www.nextflow.io/docs/latest/reference/channel.html#frompath) ins Spiel: eine praktische Kanal-Factory, die alle Funktionen bündelt, die wir benötigen, um einen Kanal aus einem oder mehreren statischen Datei-Strings sowie Glob-Mustern zu erstellen.

### 3.1. Die Kanal-Factory hinzufügen

Aktualisieren wir unseren Workflow, um `channel.fromPath` zu verwenden.

=== "Danach"

    ```groovy title="main.nf" linenums="7" hl_lines="1-3"
        // Dateien mit channel.fromPath laden
        ch_files = channel.fromPath('data/patientA_rep1_normal_R1_001.fastq.gz')
        ch_files.view { myFile -> "Found file: $myFile" }

        // Dateiattribute ausgeben
        /* Diese vorerst auskommentieren, wir kommen später darauf zurück!
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
        */

        // Zeilen in der Datei zählen
        // COUNT_LINES(myFile)
    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="7" hl_lines="1-2"
        // Ein Path-Objekt aus einem String-Pfad erstellen
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Dateiattribute ausgeben
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"

        // Zeilen in der Datei zählen
        COUNT_LINES(myFile)
    ```

Wir haben auch den Code, der die Attribute ausgibt, vorerst auskommentiert und eine `.view`-Anweisung hinzugefügt, die stattdessen nur den Dateinamen ausgibt.

Workflow ausführen:

```bash
nextflow run main.nf
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [grave_meucci] DSL2 - revision: b09964a583

    Found file: /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz
    ```

Wie du siehst, wird der Dateipfad als `Path`-Objekt im Kanal geladen.
Das ist ähnlich wie bei `file()`, aber jetzt haben wir einen Kanal, in den wir bei Bedarf weitere Dateien laden können.

`channel.fromPath()` ist eine praktische Möglichkeit, einen neuen Kanal zu erstellen, der mit einer Liste von Dateien befüllt ist.

### 3.2. Attribute von Dateien im Kanal anzeigen

Bei unserem ersten Einsatz der Kanal-Factory haben wir den Code vereinfacht und nur den Dateinamen ausgegeben.

Gehen wir zurück zur Ausgabe der vollständigen Dateiattribute:

=== "Danach"

    ```groovy title="main.nf" linenums="7" hl_lines="3-9 12"
        // Dateien mit channel.fromPath laden
        ch_files = channel.fromPath('data/patientA_rep1_normal_R1_001.fastq.gz')
        ch_files.view { myFile ->
            println "File object class: ${myFile.class}"
            println "File name: ${myFile.name}"
            println "Simple name: ${myFile.simpleName}"
            println "Extension: ${myFile.extension}"
            println "Parent directory: ${myFile.parent}"
        }

        // Zeilen in der Datei zählen
        COUNT_LINES(ch_files)
    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="7" hl_lines="3"
        // Dateien mit channel.fromPath laden
        ch_files = channel.fromPath('data/patientA_rep1_normal_R1_001.fastq.gz')
        ch_files.view { myFile -> "Found file: $myFile" }

        // Zeilen in der Datei zählen
        // COUNT_LINES(ch_files)
    ```

Wir aktivieren auch den `COUNT_LINES`-Prozessaufruf wieder, um zu überprüfen, dass die Dateiverarbeitung mit unserem kanalbasierten Ansatz weiterhin korrekt funktioniert.

Workflow ausführen:

```bash
nextflow run main.nf
```

??? success "Befehlsausgabe"

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

Und da haben wir es – dieselben Ergebnisse wie zuvor, aber jetzt haben wir die Datei in einem Kanal, sodass wir weitere hinzufügen können.

### 3.3. Einen Glob verwenden, um mehrere Dateien zu matchen

Es gibt mehrere Möglichkeiten, weitere Dateien in den Kanal zu laden.
Hier zeigen wir dir, wie du Glob-Muster verwendest, eine praktische Methode zum Matchen und Abrufen von Datei- und Verzeichnisnamen basierend auf Platzhalterzeichen.
Der Prozess des Matchens dieser Muster wird als „Globbing" oder „Dateinamen-Expansion" bezeichnet.

!!! note "Hinweis"

    Wie bereits erwähnt, unterstützt Nextflow Globbing zur Verwaltung von Eingabe- und Ausgabedateien in den meisten Fällen, außer bei HTTPS-Dateipfaden, da HTTPS keine mehreren Dateien auflisten kann.

Angenommen, wir möchten beide Dateien in einem Dateipaar abrufen, das einem bestimmten Patienten `patientA` zugeordnet ist:

```console
patientA_rep1_normal_R1_001.fastq.gz
patientA_rep1_normal_R2_001.fastq.gz
```

Da der einzige Unterschied zwischen den Dateinamen die Replikatnummer ist, d. h. die Zahl nach `R`, können wir das Platzhalterzeichen `*` wie folgt als Ersatz für die Zahl verwenden:

```console
patientA_rep1_normal_R*_001.fastq.gz
```

Das ist das Glob-Muster, das wir benötigen.

Jetzt müssen wir nur noch den Dateipfad in der Kanal-Factory aktualisieren, um dieses Glob-Muster zu verwenden:

=== "Danach"

    ```groovy title="main.nf" linenums="7"
      // Dateien mit channel.fromPath laden
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="7"
      // Dateien mit channel.fromPath laden
        ch_files = channel.fromPath('data/patientA_rep1_normal_R1_001.fastq.gz')
    ```

Nextflow erkennt automatisch, dass es sich um ein Glob-Muster handelt, und verarbeitet es entsprechend.

Führe den Workflow aus, um das zu testen:

```bash
nextflow run main.nf
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [boring_sammet] DSL2 - revision: d2aa789c9a

    executor >  local (2)
    [3c/a65de5] COUNT_LINES (2) [100%] 2 of 2 ✔
    File object class: class sun.nio.fs.UnixPath
    File name: patientA_rep1_normal_R1_001.fastq.gz
    Simple name: patientA_rep1_normal_R1_001
    Extension: gz
    Parent directory: /workspaces/training/side-quests/working_with_files/data
    File object class: class sun.nio.fs.UnixPath
    File name: patientA_rep1_normal_R2_001.fastq.gz
    Simple name: patientA_rep1_normal_R2_001
    Extension: gz
    Parent directory: /workspaces/training/side-quests/working_with_files/data
    Processing file: patientA_rep1_normal_R1_001.fastq.gz
    40

    Processing file: patientA_rep1_normal_R2_001.fastq.gz
    40
    ```

Wie du siehst, haben wir jetzt zwei Path-Objekte in unserem Kanal, was zeigt, dass Nextflow die Dateinamen-Expansion korrekt durchgeführt und beide Dateien wie erwartet geladen und verarbeitet hat.

Mit dieser Methode können wir beliebig viele Dateien abrufen, indem wir einfach das Glob-Muster ändern. Wenn wir es großzügiger gestalten, zum Beispiel indem wir alle variablen Teile der Dateinamen durch `*` ersetzen (_z. B._ `data/patient*_rep*_*_R*_001.fastq.gz`), könnten wir alle Beispieldateien im `data`-Verzeichnis abrufen.

### Fazit

- `channel.fromPath()` erstellt einen Kanal mit Dateien, die einem Muster entsprechen
- Jede Datei wird als separates Element im Kanal ausgegeben
- Wir können ein Glob-Muster verwenden, um mehrere Dateien zu matchen
- Dateien werden automatisch in Path-Objekte mit vollständigen Attributen umgewandelt
- Die `.view()`-Methode ermöglicht die Inspektion von Kanalinhalten

---

## 4. Grundlegende Metadaten aus Dateinamen extrahieren

In den meisten wissenschaftlichen Bereichen ist es sehr üblich, Metadaten in den Namen der Datendateien zu kodieren.
In der Bioinformatik werden Dateien mit Sequenzierungsdaten beispielsweise oft so benannt, dass Informationen über die Probe, den Zustand, das Replikat und die Read-Nummer kodiert sind.

Wenn die Dateinamen nach einer konsistenten Konvention aufgebaut sind, kannst du diese Metadaten auf standardisierte Weise extrahieren und im Verlauf deiner Analyse verwenden.
Das ist natürlich ein großes „Wenn", und du solltest sehr vorsichtig sein, wenn du dich auf die Dateinamenstruktur verlässt; aber die Realität ist, dass dieser Ansatz sehr weit verbreitet ist, also schauen wir uns an, wie es in Nextflow gemacht wird.

Im Fall unserer Beispieldaten wissen wir, dass die Dateinamen konsistent strukturierte Metadaten enthalten.
Zum Beispiel kodiert der Dateiname `patientA_rep1_normal_R2_001` folgendes:

- Patienten-ID: `patientA`
- Replikat-ID: `rep1`
- Probentyp: `normal` (im Gegensatz zu `tumor`)
- Read-Set: `R1` (im Gegensatz zu `R2`)

Wir werden unseren Workflow in drei Schritten modifizieren, um diese Informationen abzurufen:

1. Den `simpleName` der Datei abrufen, der die Metadaten enthält
2. Die Metadaten mit einer Methode namens `tokenize()` trennen
3. Eine map verwenden, um die Metadaten zu organisieren

!!! warning "Warnung"

    Du solltest niemals sensible Informationen in Dateinamen kodieren, wie z. B. Patientennamen oder andere identifizierende Merkmale, da dies die Privatsphäre von Patient\*innen oder andere relevante Sicherheitsbeschränkungen gefährden kann.

### 4.1. Den `simpleName` abrufen

Der `simpleName` ist ein Dateiattribut, das dem Dateinamen ohne Pfad und Erweiterung entspricht.

Nimm folgende Änderungen am Workflow vor:

=== "Danach"

    ```groovy title="main.nf" linenums="7" hl_lines="3-6"
        // Dateien mit channel.fromPath laden
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
        ch_files.map { myFile ->
            [ myFile.simpleName, myFile ]
        }
        .view()
    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="7" hl_lines="3-9"
        // Dateien mit channel.fromPath laden
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
        ch_files.view { myFile ->
            println "File object class: ${myFile.class}"
            println "File name: ${myFile.name}"
            println "Simple name: ${myFile.simpleName}"
            println "Extension: ${myFile.extension}"
            println "Parent directory: ${myFile.parent}"
        }
    ```

Das ruft den `simpleName` ab und verknüpft ihn mit dem vollständigen Dateiobjekt mithilfe einer `map()`-Operation.

Führe den Workflow aus, um zu testen, ob es funktioniert:

```bash
nextflow run main.nf
```

??? success "Befehlsausgabe"

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

Jedes Element im Kanal ist jetzt ein Tupel, das den `simpleName` und das ursprüngliche Dateiobjekt enthält.

### 4.2. Die Metadaten aus dem `simpleName` extrahieren

An diesem Punkt sind die gewünschten Metadaten im `simpleName` eingebettet, aber wir können nicht direkt auf einzelne Elemente zugreifen.
Wir müssen also den `simpleName` in seine Bestandteile aufteilen.
Glücklicherweise sind diese Bestandteile im ursprünglichen Dateinamen einfach durch Unterstriche getrennt, sodass wir eine gängige Nextflow-Methode namens `tokenize()` anwenden können, die für diese Aufgabe perfekt geeignet ist.

Nimm folgende Änderungen am Workflow vor:

=== "Danach"

    ```groovy title="main.nf" linenums="7" hl_lines="4"
        // Dateien mit channel.fromPath laden
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
        ch_files.map { myFile ->
            [ myFile.simpleName.tokenize('_'), myFile ]
        }
    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="7" hl_lines="4"
        // Dateien mit channel.fromPath laden
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
        ch_files.map { myFile ->
            [ myFile.simpleName, myFile ]
        }
    ```

Die `tokenize()`-Methode teilt den `simpleName`-String überall dort auf, wo sie Unterstriche findet, und gibt eine Liste mit den Teilstrings zurück.

Workflow ausführen:

```bash
nextflow run main.nf
```

??? success "Befehlsausgabe"

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

Jetzt enthält das Tupel für jedes Element in unserem Kanal die Liste der Metadaten (_z. B._ `[patientA, rep1, normal, R1, 001]`) und das ursprüngliche Dateiobjekt.

Großartig!
Wir haben unsere Patienteninformationen von einem einzelnen String in eine Liste von Strings aufgeteilt.
Wir können jetzt jeden Teil der Patienteninformationen separat verarbeiten.

### 4.3. Eine map verwenden, um die Metadaten zu organisieren

Unsere Metadaten sind im Moment nur eine flache Liste.
Sie ist einfach genug zu verwenden, aber schwer zu lesen.

```console
[patientA, rep1, normal, R1, 001]
```

Was ist das Element an Index 3? Kannst du es sagen, ohne auf die ursprüngliche Erklärung der Metadatenstruktur zurückzugreifen?

Das ist eine gute Gelegenheit, einen Key-Value-Store zu verwenden, bei dem jedes Element einen Satz von Schlüsseln und ihren zugehörigen Werten hat, sodass du leicht auf jeden Schlüssel verweisen kannst, um den entsprechenden Wert zu erhalten.

In unserem Beispiel bedeutet das, von dieser Organisation:

```groovy
data = [patientA, 1, normal, R1]

println data[3]
```

Zu dieser hier:

```groovy
data = [id: patientA, replicate: 1, type: normal, readNum: 1]

println data.readNum
```

In Nextflow nennt man das eine [map](https://nextflow.io/docs/latest/script.html#maps).

Konvertieren wir jetzt unsere flache Liste in eine map.
Nimm folgende Änderungen am Workflow vor:

=== "Danach"

    ```groovy title="main.nf" linenums="7" hl_lines="4-13"
        // Dateien mit channel.fromPath laden
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

=== "Vorher"

    ```groovy title="main.nf" linenums="7" hl_lines="4"
        // Dateien mit channel.fromPath laden
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
        ch_files.map { myFile ->
            [ myFile.simpleName.tokenize('_'), myFile ]
        }
    ```

Die wichtigsten Änderungen hier sind:

- **Destrukturierungs-Zuweisung**: `def (patient, replicate, type, readNum) = ...` extrahiert die tokenisierten Werte in einer Zeile in benannte Variablen
- **Map-Literal-Syntax**: `[id: patient, replicate: ...]` erstellt eine map, bei der jeder Schlüssel (wie `id`) mit einem Wert (wie `patient`) verknüpft ist
- **Verschachtelte Struktur**: Die äußere Liste `[..., myFile]` verbindet die Metadaten-map mit dem ursprünglichen Dateiobjekt

Wir haben auch einige der Metadaten-Strings mit einer String-Ersetzungsmethode namens `replace()` vereinfacht, um einige unnötige Zeichen zu entfernen (_z. B._ `replicate.replace('rep', '')`, um nur die Zahl aus den Replikat-IDs zu behalten).

Führen wir den Workflow erneut aus:

```bash
nextflow run main.nf
```

??? success "Befehlsausgabe"

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

Jetzt sind die Metadaten übersichtlich beschriftet (_z. B._ `[id:patientA, replicate:1, type:normal, readNum:2]`), sodass es viel einfacher ist zu erkennen, was was ist.

Es wird auch viel einfacher sein, Metadaten-Elemente im Workflow tatsächlich zu nutzen, und unser Code wird leichter zu lesen und zu warten sein.

### Fazit

- Wir können Dateinamen in Nextflow mit der Leistung einer vollständigen Programmiersprache verarbeiten
- Wir können die Dateinamen als Strings behandeln, um relevante Informationen zu extrahieren
- Die Verwendung von Methoden wie `tokenize()` und `replace()` ermöglicht es uns, Strings im Dateinamen zu manipulieren
- Die `.map()`-Operation transformiert Kanalelemente und bewahrt dabei die Struktur
- Strukturierte Metadaten (maps) machen den Code lesbarer und wartbarer als Positionslisten

Als nächstes schauen wir uns an, wie man gepaarte Datendateien verarbeitet.

---

## 5. Gepaarte Datendateien verarbeiten

Viele experimentelle Designs erzeugen gepaarte Datendateien, die davon profitieren, explizit gepaart verarbeitet zu werden.
In der Bioinformatik werden Sequenzierungsdaten beispielsweise oft in Form von gepaarten Reads generiert, d. h. Sequenzstrings, die aus demselben DNA-Fragment stammen (oft als 'Forward' und 'Reverse' bezeichnet, weil sie von entgegengesetzten Enden gelesen werden).

Das ist der Fall bei unseren Beispieldaten, wo R1 und R2 auf die beiden Read-Sets verweisen.

```console
data/patientA_rep1_normal_R1_001.fastq.gz
data/patientA_rep1_normal_R2_001.fastq.gz
```

Nextflow bietet eine spezialisierte Kanal-Factory für die Arbeit mit solchen gepaarten Dateien namens `channel.fromFilePairs()`, die Dateien automatisch basierend auf einem gemeinsamen Benennungsmuster gruppiert. Das ermöglicht es dir, die gepaarten Dateien mit weniger Aufwand enger miteinander zu verknüpfen.

Wir werden unseren Workflow modifizieren, um davon zu profitieren.
Das erfordert zwei Schritte:

1. Die Kanal-Factory auf `channel.fromFilePairs()` umstellen
2. Die Metadaten extrahieren und mappen

### 5.1. Die Kanal-Factory auf `channel.fromFilePairs()` umstellen

Um `channel.fromFilePairs` zu verwenden, müssen wir das Muster angeben, das Nextflow verwenden soll, um die beiden Mitglieder eines Paares zu identifizieren.

Zurück zu unseren Beispieldaten können wir das Benennungsmuster wie folgt formalisieren:

```console
data/patientA_rep1_normal_R{1,2}_001.fastq.gz
```

Das ist ähnlich wie das Glob-Muster, das wir zuvor verwendet haben, außer dass es die Teilstrings (entweder `1` oder `2` direkt nach dem R) explizit aufzählt, die die beiden Mitglieder des Paares identifizieren.

Aktualisieren wir den Workflow `main.nf` entsprechend:

=== "Danach"

    ```groovy title="main.nf" linenums="7" hl_lines="1-2"
        // Dateien mit channel.fromFilePairs laden
        ch_files = channel.fromFilePairs('data/patientA_rep1_normal_R{1,2}_001.fastq.gz')
        /* Die Mapping-Operation vorerst auskommentieren, wir kommen später darauf zurück!
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

=== "Vorher"

    ```groovy title="main.nf" linenums="7" hl_lines="1-2"
        // Dateien mit channel.fromFilePairs laden
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

Wir haben die Kanal-Factory umgestellt und das Datei-Matching-Muster angepasst, und dabei auch die map-Operation auskommentiert.
Wir werden sie später mit einigen Änderungen wieder hinzufügen.

Führe den Workflow aus, um ihn zu testen:

```bash
nextflow run main.nf
```

??? failure "Befehlsausgabe"

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

Oh nein, diesmal ist die Ausführung fehlgeschlagen!

Der relevante Teil der Fehlermeldung ist hier:

```console
Not a valid path value: 'patientA_rep1_normal_R'
```

Das liegt daran, dass wir die Kanal-Factory geändert haben.
Bisher enthielt der ursprüngliche Eingabekanal nur die Dateipfade.
Alle Metadaten-Manipulationen, die wir durchgeführt haben, haben den Kanalinhalt nicht tatsächlich verändert.

Jetzt, da wir die `.fromFilePairs`-Kanal-Factory verwenden, ist der Inhalt des resultierenden Kanals anders.
Wir sehen nur ein Kanalelement, das aus einem Tupel mit zwei Elementen besteht: dem Teil des `simpleName`, der von den beiden Dateien geteilt wird und als Bezeichner dient, und einem Tupel mit den beiden Dateiobjekten im Format `id, [ file1, file2 ]`.

Das ist großartig, denn Nextflow hat die schwere Arbeit erledigt, den Patientennamen zu extrahieren, indem es das gemeinsame Präfix untersucht und es als Patientenbezeichner verwendet hat.

Allerdings bricht das unseren aktuellen Workflow.
Wenn wir `COUNT_LINES` weiterhin auf die gleiche Weise ausführen wollten, ohne den Prozess zu ändern, müssten wir eine Mapping-Operation anwenden, um die Dateipfade zu extrahieren.
Das werden wir aber nicht tun, da unser eigentliches Ziel ist, einen anderen Prozess, `ANALYZE_READS`, zu verwenden, der Dateipaare entsprechend verarbeitet.

Also kommentieren wir einfach den Aufruf von `COUNT_LINES` aus (oder löschen ihn) und machen weiter.

=== "Danach"

    ```groovy title="main.nf" linenums="26" hl_lines="2"
        // Zeilen in der Datei zählen
        // COUNT_LINES(ch_files)
    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="26" hl_lines="2"
        // Zeilen in der Datei zählen
        COUNT_LINES(ch_files)
    ```

Du kannst auch die `COUNT_LINES`-Include-Anweisung auskommentieren oder löschen, aber das hat keine funktionalen Auswirkungen.

Führen wir den Workflow erneut aus:

```bash
nextflow run main.nf
```

??? success "Befehlsausgabe"

    ```console hl_lines="5"
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [fabulous_davinci] DSL2 - revision: 22b53268dc

    [patientA_rep1_normal_R, [/workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz, /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R2_001.fastq.gz]]
    ```

Diesmal ist der Workflow erfolgreich!

Wir müssen jedoch noch die restlichen Metadaten aus dem `id`-Feld extrahieren.

### 5.2. Metadaten aus Dateipaaren extrahieren und organisieren

Unsere `map`-Operation von vorhin wird nicht funktionieren, weil sie nicht zur Datenstruktur passt, aber wir können sie anpassen.

Wir haben bereits Zugriff auf den eigentlichen Patientenbezeichner in dem String, den `fromFilePairs()` als Bezeichner verwendet hat, sodass wir diesen verwenden können, um die Metadaten zu extrahieren, ohne den `simpleName` aus dem Path-Objekt zu holen, wie wir es zuvor getan haben.

Kommentiere die map-Operation im Workflow aus und nimm folgende Änderungen vor:

=== "Danach"

    ```groovy title="main.nf" linenums="7" hl_lines="3-4 9 11 13"
        // Dateien mit channel.fromFilePairs laden
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

=== "Vorher"

    ```groovy title="main.nf" linenums="7" hl_lines="3-5 11 13"
        // Dateien mit channel.fromFilePairs laden
        ch_files = channel.fromFilePairs('data/patientA_rep1_normal_R{1,2}_001.fastq.gz')
        /* Die Mapping-Operation vorerst auskommentieren, wir kommen später darauf zurück!
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

Diesmal beginnt die map mit `id, files` statt nur `myFile`, und `tokenize()` wird auf `id` statt auf `myFile.simpleName` angewendet.

Beachte auch, dass wir `readNum` aus der `tokenize()`-Zeile entfernt haben; alle Teilstrings, die wir nicht explizit benennen (von links beginnend), werden stillschweigend verworfen.
Das können wir tun, weil die gepaarten Dateien jetzt eng miteinander verknüpft sind, sodass wir `readNum` nicht mehr in der Metadaten-map benötigen.

Führen wir den Workflow aus:

```bash
nextflow run main.nf
```

??? success "Befehlsausgabe"

    ```console

    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [prickly_stonebraker] DSL2 - revision: f62ab10a3f

    [[id:patientA, replicate:1, type:normal], [/workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz, /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R2_001.fastq.gz]]
    ```

Und da ist es: Wir haben die Metadaten-map (`[id:patientA, replicate:1, type:normal]`) an der ersten Position des Ausgabe-Tupels, gefolgt vom Tupel der gepaarten Dateien, wie beabsichtigt.

Natürlich wird das nur dieses spezifische Dateipaar aufnehmen und verarbeiten.
Wenn du mit der Verarbeitung mehrerer Paare experimentieren möchtest, kannst du versuchen, Platzhalter in das Eingabemuster einzufügen und sehen, was passiert.
Versuche zum Beispiel `data/patientA_rep1_*_R{1,2}_001.fastq.gz` zu verwenden.

### Fazit

- [`channel.fromFilePairs()` findet und paart verwandte Dateien automatisch](https://www.nextflow.io/docs/latest/reference/channel.html#fromfilepairs)
- Das vereinfacht die Verarbeitung von Paired-End-Reads in deiner Pipeline
- Gepaarte Dateien können als `[id, [file1, file2]]`-Tupel gruppiert werden
- Die Metadatenextraktion kann aus der gepaarten Datei-ID statt aus einzelnen Dateien erfolgen

---

## 6. Dateioperationen in Prozessen verwenden

Jetzt bringen wir das alles in einem einfachen Prozess zusammen, um zu festigen, wie man Dateioperationen innerhalb eines Nextflow-Prozesses verwendet.

Wir stellen dir ein vorgefertigtes Prozessmodul namens `ANALYZE_READS` zur Verfügung, das ein Tupel aus Metadaten und einem Paar Eingabedateien entgegennimmt und sie analysiert.
Wir könnten uns vorstellen, dass dies ein Sequenz-Alignment, Variant Calling oder einen anderen Schritt durchführt, der für diesen Datentyp sinnvoll ist.

Fangen wir an.

### 6.1. Den Prozess importieren und den Code untersuchen

Um diesen Prozess im Workflow zu verwenden, müssen wir nur eine Modul-Include-Anweisung vor dem Workflow-Block hinzufügen.

Nimm folgende Änderung am Workflow vor:

=== "Danach"

    ```groovy title="main.nf" linenums="1" hl_lines="3"
    #!/usr/bin/env nextflow

    include { ANALYZE_READS } from './modules/analyze_reads.nf'

    workflow {
    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="1"
    #!/usr/bin/env nextflow

    workflow {
    ```

Du kannst die Moduldatei öffnen, um ihren Code zu untersuchen:

```groovy title="modules/analyze_reads.nf - process example" linenums="1"
#!/usr/bin/env nextflow

process ANALYZE_READS {
    tag { meta.id }

    publishDir { "results/${meta.id}" }, mode: 'copy'

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

!!! note "Hinweis"

    Die `tag`- und `publishDir`-Direktiven verwenden Closure-Syntax (`{ ... }`) statt String-Interpolation (`"${...}"`).
    Das liegt daran, dass diese Direktiven auf Eingabevariablen (`meta`) verweisen, die erst zur Laufzeit verfügbar sind.
    Die Closure-Syntax verzögert die Auswertung, bis der Prozess tatsächlich ausgeführt wird.

!!! note "Hinweis"

    Wir nennen unsere Metadaten-map konventionsgemäß `meta`.
    Für einen tieferen Einblick in Meta-Maps, siehe die Side Quest [Metadata and meta maps](../metadata/).

### 6.2. Den Prozess im Workflow aufrufen

Da der Prozess jetzt im Workflow verfügbar ist, können wir einen Aufruf des `ANALYZE_READS`-Prozesses hinzufügen, um ihn auszuführen.

Um ihn auf unseren Beispieldaten auszuführen, müssen wir zwei Dinge tun:

1. Dem remappten Kanal einen Namen geben
2. Einen Aufruf des Prozesses hinzufügen

#### 6.2.1. Den remappten Eingabekanal benennen

Wir haben die Mapping-Manipulationen bisher direkt auf den Eingabekanal angewendet.
Um den remappten Inhalt an den `ANALYZE_READS`-Prozess zu übergeben (und das auf eine klare und leicht lesbare Weise), wollen wir einen neuen Kanal namens `ch_samples` erstellen.

Das können wir mit dem [`set`](https://www.nextflow.io/docs/latest/reference/operator.html#set)-Operator tun.

Ersetze im Haupt-Workflow den `.view()`-Operator durch `.set { ch_samples }` und füge eine Zeile hinzu, die testet, ob wir den Kanal beim Namen ansprechen können.

=== "Danach"

    ```groovy title="main.nf" linenums="7" hl_lines="14 16-17"
        // Dateien mit channel.fromFilePairs laden
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

        // Temporär: Blick in ch_samples
        ch_samples.view()
    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="7" hl_lines="14"
        // Dateien mit channel.fromFilePairs laden
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

Führen wir das aus:

```bash
nextflow run main.nf
```

??? success "Befehlsausgabe"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [goofy_kirch] DSL2 - revision: 3313283e42

    [[id:patientA, replicate:1, type:normal], [/workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz, /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R2_001.fastq.gz]]
    ```

Das bestätigt, dass wir den Kanal jetzt beim Namen ansprechen können.

#### 6.2.2. Den Prozess auf den Daten aufrufen

Jetzt rufen wir den `ANALYZE_READS`-Prozess tatsächlich auf dem `ch_samples`-Kanal auf.

Nimm folgende Codeänderungen im Haupt-Workflow vor:

=== "Danach"

    ```groovy title="main.nf" linenums="23"
        // Die Analyse ausführen
        ANALYZE_READS(ch_samples)
    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="23"
        // Temporär: Blick in ch_samples
        ch_samples.view()
    ```

Führen wir das aus:

```bash
nextflow run main.nf
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `./main.nf` [shrivelled_cori] DSL2 - revision: b546a31769

    executor >  local (1)
    [b5/110360] process > ANALYZE_READS (patientA) [100%] 1 of 1 ✔
    ```

Dieser Prozess ist so eingerichtet, dass er seine Ausgaben in ein `results`-Verzeichnis veröffentlicht, also schau dort nach.

??? abstract "Verzeichnis- und Dateiinhalte"

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

Der Prozess hat unsere Eingaben entgegengenommen und eine neue Datei mit den Patientenmetadaten erstellt, wie vorgesehen.
Ausgezeichnet!

### 6.3. Viele weitere Patient\*innen einbeziehen

Natürlich verarbeiten wir hier nur ein einzelnes Dateipaar für eine\*n einzelne\*n Patient\*in, was nicht gerade der Hochdurchsatz ist, den du dir von Nextflow erhoffst.
Du wirst wahrscheinlich viel mehr Daten auf einmal verarbeiten wollen.

Denke daran, dass `channel.fromPath()` einen _Glob_ als Eingabe akzeptiert, was bedeutet, dass es eine beliebige Anzahl von Dateien akzeptieren kann, die dem Muster entsprechen.
Wenn wir also alle Patient\*innen einbeziehen wollen, können wir einfach den Eingabe-String ändern, um mehr Patient\*innen einzuschließen, wie bereits kurz erwähnt.

Stellen wir uns vor, wir wollen so viele Daten wie möglich verarbeiten.
Nimm folgende Änderungen am Workflow vor:

=== "Danach"

    ```groovy title="main.nf" linenums="7" hl_lines="2"
        // Dateien mit channel.fromFilePairs laden
        ch_files = channel.fromFilePairs('data/*_R{1,2}_001.fastq.gz')
    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="7" hl_lines="2"
        // Dateien mit channel.fromFilePairs laden
        ch_files = channel.fromFilePairs('data/patientA_rep1_normal_R{1,2}_001.fastq.gz')
    ```

Führe die Pipeline erneut aus:

```bash
nextflow run main.nf
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `./main.nf` [big_stonebraker] DSL2 - revision: f7f9b8a76c

    executor >  local (8)
    [d5/441891] process > ANALYZE_READS (patientC) [100%] 8 of 8 ✔
    ```

Das `results`-Verzeichnis sollte jetzt Ergebnisse für alle verfügbaren Daten enthalten.

??? abstract "Verzeichnisinhalt"

    ```console
    results
    ├── patientA
    │   └── patientA_stats.txt
    ├── patientB
    │   └── patientB_stats.txt
    └── patientC
        └── patientC_stats.txt
    ```

Erfolg! Wir haben alle Patient\*innen auf einmal analysiert! Oder?

Vielleicht nicht.
Wenn du genauer hinschaust, haben wir ein Problem: Wir haben zwei Replikate für patientA, aber nur eine Ausgabedatei!
Wir überschreiben die Ausgabedatei jedes Mal.

### 6.4. Die veröffentlichten Dateien eindeutig machen

Da wir Zugriff auf die Patientenmetadaten haben, können wir sie verwenden, um die veröffentlichten Dateien eindeutig zu machen, indem wir differenzierende Metadaten entweder in die Verzeichnisstruktur oder in die Dateinamen selbst einbeziehen.

Nimm folgende Änderung am Workflow vor:

=== "Danach"

    ```groovy title="modules/analyze_reads.nf" linenums="6"
        publishDir { "results/${meta.type}/${meta.id}/${meta.replicate}" }, mode: 'copy'
    ```

=== "Vorher"

    ```groovy title="modules/analyze_reads.nf" linenums="6"
        publishDir { "results/${meta.id}" }, mode: 'copy'
    ```

Hier zeigen wir die Option, zusätzliche Verzeichnisebenen zu verwenden, um Probentypen und Replikate zu berücksichtigen, aber du könntest auch auf Dateinamenebene experimentieren.

Führe die Pipeline jetzt noch einmal aus, aber stelle sicher, dass du zuerst das `results`-Verzeichnis entfernst, um einen sauberen Arbeitsbereich zu haben:

```bash
rm -r results
nextflow run main.nf
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `./main.nf` [insane_swartz] DSL2 - revision: fff18abe6d

    executor >  local (8)
    [e3/449081] process > ANALYZE_READS (patientC) [100%] 8 of 8 ✔
    ```

Überprüfe jetzt das `results`-Verzeichnis:

??? abstract "Verzeichnisinhalt"

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

Und da ist es – alle unsere Metadaten, ordentlich organisiert. Das ist Erfolg!

Es gibt noch viel mehr, was du tun kannst, sobald du deine Metadaten in eine map geladen hast:

1. Organisierte Ausgabeverzeichnisse basierend auf Patientenattributen erstellen
2. Entscheidungen in Prozessen basierend auf Patienteneigenschaften treffen
3. Daten basierend auf Metadatenwerten aufteilen, zusammenführen und neu kombinieren

Dieses Muster, Metadaten explizit und an die Daten gebunden zu halten (anstatt sie in Dateinamen zu kodieren), ist eine grundlegende Best Practice in Nextflow, die den Aufbau robuster, wartbarer Analyse-Workflows ermöglicht.
Mehr dazu erfährst du in der Side Quest [Metadata and meta maps](../metadata/).

### Fazit

- Die `publishDir`-Direktive kann Ausgaben basierend auf Metadatenwerten organisieren
- Metadaten in Tupeln ermöglichen eine strukturierte Organisation der Ergebnisse
- Dieser Ansatz erstellt wartbare Workflows mit klarer Datenprovenienz
- Prozesse können Tupel aus Metadaten und Dateien als Eingabe entgegennehmen
- Die `tag`-Direktive liefert Prozessidentifikation in Ausführungsprotokollen
- Die Workflow-Struktur trennt die Kanalerstellung von der Prozessausführung

---

## Zusammenfassung

In dieser Side Quest hast du gelernt, wie man mit Dateien in Nextflow arbeitet – von grundlegenden Operationen bis hin zu fortgeschritteneren Techniken für die Verarbeitung von Dateisammlungen.

Die Anwendung dieser Techniken in deiner eigenen Arbeit ermöglicht es dir, effizientere und wartbarere Workflows zu erstellen, insbesondere wenn du mit großen Mengen von Dateien mit komplexen Benennungskonventionen arbeitest.

### Wichtige Muster

1.  **Grundlegende Dateioperationen:** Wir haben Path-Objekte mit `file()` erstellt und auf Dateiattribute wie Name, Erweiterung und übergeordnetes Verzeichnis zugegriffen und dabei den Unterschied zwischen Strings und Path-Objekten kennengelernt.

    - Ein Path-Objekt mit `file()` erstellen

    ```groovy
    myFile = file('path/to/file.txt')
    ```

    - Dateiattribute abrufen

    ```groovy
    println myFile.name       // file.txt
    println myFile.baseName   // file
    println myFile.extension  // txt
    println myFile.parent     // path/to
    ```

2.  **Entfernte Dateien verwenden**: Wir haben gelernt, wie man transparent zwischen lokalen und entfernten Dateien über URIs wechselt, und dabei Nextflows Fähigkeit demonstriert, Dateien aus verschiedenen Quellen zu verarbeiten, ohne die Workflow-Logik zu ändern.

    - Lokale Datei

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

3.  **Dateien mit der `fromPath()`-Kanal-Factory laden:** Wir haben Kanäle aus Dateimustern mit `channel.fromPath()` erstellt und ihre Dateiattribute einschließlich Objekttypen angezeigt.

    - Einen Kanal aus einem Dateimuster erstellen

    ```groovy
     ch_files = channel.fromPath('data/*.fastq.gz')
    ```

    - Dateiattribute abrufen

    ```groovy
     ch_files.view { myFile ->
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    }
    ```

4.  **Patientenmetadaten aus Dateinamen extrahieren:** Wir haben `tokenize()` und `replace()` verwendet, um Metadaten aus Dateinamen zu extrahieren und zu strukturieren und sie in organisierte maps umzuwandeln.

    ```groovy
    def name = file.name.tokenize('_')
    def patientId = name[0]
    def replicate = name[1].replace('rep', '')
    def type = name[2]
    def readNum = name[3].replace('R', '')
    ```

5.  **Vereinfachung mit channel.fromFilePairs:** Wir haben `channel.fromFilePairs()` verwendet, um verwandte Dateien automatisch zu paaren und Metadaten aus gepaarten Datei-IDs zu extrahieren.

    ```groovy
    ch_pairs = channel.fromFilePairs('data/*_R{1,2}_001.fastq.gz')
    ```

6.  **Dateioperationen in Prozessen verwenden:** Wir haben Dateioperationen mit korrekter Eingabeverarbeitung in Nextflow-Prozesse integriert und `publishDir` verwendet, um Ausgaben basierend auf Metadaten zu organisieren.

    - Eine Meta-Map mit den Prozesseingaben verknüpfen

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

    - Ausgaben basierend auf Metadaten organisieren

    ```groovy
    publishDir { "results/${meta.type}/${meta.id}/${meta.replicate}" }, mode: 'copy'
    ```

### Weitere Ressourcen

- [Nextflow-Dokumentation: Working with Files](https://www.nextflow.io/docs/latest/working-with-files.html)
- [channel.fromPath](https://www.nextflow.io/docs/latest/reference/channel.html#frompath)
- [channel.fromFilePairs](https://www.nextflow.io/docs/latest/reference/channel.html#fromfilepairs)

---

## Wie geht es weiter?

Kehre zum [Menü der Side Quests](../) zurück oder klicke auf die Schaltfläche unten rechts auf der Seite, um zum nächsten Thema in der Liste zu gelangen.
