# Datei-Eingabeverarbeitung

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Wissenschaftliche Analyse-Workflows beinhalten oft die Verarbeitung großer Mengen von Dateien.
Nextflow bietet leistungsstarke Tools, um Dateien effizient zu handhaben und hilft dir, deine Daten mit minimalem Code zu organisieren und zu verarbeiten.

### Lernziele

In dieser Side Quest werden wir erkunden, wie Nextflow mit Dateien umgeht, von grundlegenden Dateioperationen bis hin zu fortgeschritteneren Techniken für die Arbeit mit Dateisammlungen.
Du lernst, wie du Metadaten aus Dateinamen extrahierst, was eine häufige Anforderung in wissenschaftlichen Analyse-Pipelines ist.

Am Ende dieser Side Quest wirst du in der Lage sein:

- Path-Objekte aus Dateipfad-Strings mit Nextflows `file()`-Methode zu erstellen
- Auf Dateiattribute wie Name, Erweiterung und übergeordnetes Verzeichnis zuzugreifen
- Sowohl lokale als auch entfernte Dateien transparent mit URIs zu handhaben
- Kanäle zu verwenden, um die Dateiverarbeitung mit `channel.fromPath()` und `channel.fromFilePairs()` zu automatisieren
- Metadaten aus Dateinamen mit String-Manipulation zu extrahieren und zu strukturieren
- Zusammengehörige Dateien mit Pattern Matching und Glob-Ausdrücken zu gruppieren
- Dateioperationen mit korrekter Eingabeverarbeitung in Nextflow-Prozesse zu integrieren
- Prozess-Ausgaben mit metadatengesteuerten Verzeichnisstrukturen zu organisieren

Diese Fähigkeiten helfen dir, Workflows zu erstellen, die verschiedene Arten von Datei-Eingaben mit großer Flexibilität handhaben können.

### Voraussetzungen

Bevor du diese Side Quest beginnst, solltest du:

- Das [Hello Nextflow](../../hello_nextflow/)-Tutorial oder einen gleichwertigen Einsteigerkurs abgeschlossen haben
- Mit grundlegenden Nextflow-Konzepten und -Mechanismen vertraut sein (Prozesse, Kanäle, Operatoren)

---

## 0. Erste Schritte

#### Öffne den Training-Codespace

Falls du das noch nicht getan hast, stelle sicher, dass du die Trainingsumgebung wie in der [Umgebungseinrichtung](../envsetup/index.md) beschrieben öffnest.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Wechsle in das Projektverzeichnis

Lass uns in das Verzeichnis wechseln, in dem sich die Dateien für dieses Tutorial befinden.

```bash
cd side-quests/working_with_files
```

Du kannst VSCode so einstellen, dass es sich auf dieses Verzeichnis konzentriert:

```bash
code .
```

#### Überprüfe die Materialien

Du findest eine einfache Workflow-Datei namens `main.nf`, ein `modules`-Verzeichnis mit zwei Moduldateien und ein `data`-Verzeichnis mit einigen Beispieldateien.

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

Für jede\*n Patient\*in haben wir Proben vom Typ `tumor` (typischerweise aus Tumorbiopsien) oder `normal` (aus gesundem Gewebe oder Blut).
Falls du nicht mit Krebsanalysen vertraut bist: Dies entspricht einem experimentellen Modell, das gepaarte Tumor/Normal-Proben für kontrastive Analysen verwendet.

Speziell für Patient\*in A haben wir zwei Sätze technischer Replikate (Wiederholungen).

Die Sequenzierungsdateien sind mit einer typischen `_R1_`- und `_R2_`-Konvention benannt für sogenannte 'Forward Reads' und 'Reverse Reads'.

_Mach dir keine Sorgen, wenn du mit diesem experimentellen Design nicht vertraut bist, es ist nicht entscheidend für das Verständnis dieses Tutorials._

#### Überprüfe die Aufgabe

Deine Herausforderung ist es, einen Nextflow-Workflow zu schreiben, der:

1. **Lädt** Eingabedateien mit Nextflows Dateiverarbeitungsmethoden
2. **Extrahiert** Metadaten (Patient\*innen-ID, Replikat, Probentyp) aus der Dateinamenstruktur
3. **Gruppiert** gepaarte Dateien (R1/R2) mit `channel.fromFilePairs()`
4. **Verarbeitet** die Dateien mit einem bereitgestellten Analysemodul
5. **Organisiert** Ausgaben in einer Verzeichnisstruktur basierend auf den extrahierten Metadaten

#### Bereitschafts-Checkliste

Denkst du, du bist bereit loszulegen?

- [ ] Ich verstehe das Ziel dieses Kurses und seine Voraussetzungen
- [ ] Mein Codespace läuft
- [ ] Ich habe mein Arbeitsverzeichnis entsprechend eingestellt
- [ ] Ich verstehe die Aufgabe

Wenn du alle Kästchen abhaken kannst, kann es losgehen.

---

## 1. Grundlegende Dateioperationen

### 1.1. Identifiziere den Typ eines Objekts mit `.class`

Schau dir die Workflow-Datei `main.nf` an:

```groovy title="main.nf" linenums="1"
#!/usr/bin/env nextflow

workflow {

    // Ein Path-Objekt aus einem String-Pfad erstellen
    myFile = 'data/patientA_rep1_normal_R1_001.fastq.gz'

    println "${myFile} is of class ${myFile.class}"
}
```

Dies ist ein Mini-Workflow (ohne Prozesse), der auf einen einzelnen Dateipfad in seinem Workflow verweist und ihn dann auf der Konsole ausgibt, zusammen mit seiner Klasse.

??? info "Was ist `.class`?"

    In Nextflow sagt uns `.class`, mit welchem Objekttyp wir arbeiten. Es ist wie zu fragen "was für eine Art von Ding ist das?", um herauszufinden, ob es ein String, eine Zahl, eine Datei oder etwas anderes ist.
    Dies hilft uns, den Unterschied zwischen einem einfachen String und einem Path-Objekt in den nächsten Abschnitten zu veranschaulichen.

Lass uns den Workflow ausführen:

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

Dies ist nur Textausgabe; Nextflow hat noch nichts Besonderes damit gemacht.
Wir haben auch bestätigt, dass dies für Nextflow nur ein String ist (der Klasse `java.lang.String`).
Das macht Sinn, da wir Nextflow noch nicht gesagt haben, dass es einer Datei entspricht.

### 1.2. Erstelle ein Path-Objekt mit file()

Wir können Nextflow sagen, wie es mit Dateien umgehen soll, indem wir [Path-Objekte](https://www.nextflow.io/docs/latest/reference/stdlib-types.html#path) aus Pfad-Strings erstellen.

In unserem Workflow können wir den String-Pfad `data/patientA_rep1_normal_R1_001.fastq.gz` mit der `file()`-Methode in ein Path-Objekt umwandeln, das Zugriff auf Dateieigenschaften und -operationen bietet.

Bearbeite die `main.nf`, um den String wie folgt mit `file()` zu umschließen:

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

Führe den Workflow jetzt erneut aus:

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

Nextflow hat unseren String in ein Path-Objekt umgewandelt und ihn zum tatsächlichen Dateispeicherort im System aufgelöst.
Der Dateipfad ist jetzt absolut, wie in `/workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz`.

Beachte auch, dass die Path-Objekt-Klasse `sun.nio.fs.UnixPath` ist: Dies ist Nextflows Art, lokale Dateien darzustellen.
Wie wir später sehen werden, haben entfernte Dateien andere Klassennamen (wie `nextflow.file.http.XPath` für HTTP-Dateien), aber sie funktionieren alle genau gleich und können identisch in deinen Workflows verwendet werden.

!!! tip "Tipp"

    **Der Hauptunterschied:**

    - **Pfad-String**: Nur Text, den Nextflow als Zeichen behandelt
    - **Path-Objekt**: Eine intelligente Dateireferenz, mit der Nextflow arbeiten kann

    Stell es dir so vor: Ein Pfad-String ist wie eine auf Papier geschriebene Adresse, während ein Path-Objekt wie eine in ein GPS-Gerät geladene Adresse ist, das weiß, wie man dorthin navigiert und dir Details über die Reise mitteilen kann.

### 1.3. Greife auf Dateiattribute zu

Warum ist das hilfreich? Nun, jetzt, da Nextflow versteht, dass `myFile` ein Path-Objekt und nicht nur ein String ist, können wir auf die verschiedenen Attribute des Path-Objekts zugreifen.

Lass uns unseren Workflow aktualisieren, um die eingebauten Dateiattribute auszugeben:

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

Führe den Workflow aus:

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

Du siehst die verschiedenen Dateiattribute oben auf der Konsole ausgegeben.

### 1.4. Übergib die Datei an einen Prozess

Der Unterschied zwischen Strings und Path-Objekten wird entscheidend, wenn du beginnst, tatsächliche Workflows mit Prozessen zu erstellen.
Bisher haben wir überprüft, dass Nextflow unsere Eingabedatei jetzt als Datei behandelt, aber lass uns sehen, ob wir tatsächlich etwas mit dieser Datei in einem Prozess ausführen können.

#### 1.4.1. Importiere den Prozess und untersuche den Code

Wir stellen dir ein vorgefertigtes Prozessmodul namens `COUNT_LINES` zur Verfügung, das eine Dateieingabe nimmt und zählt, wie viele Zeilen sie enthält.

Um den Prozess im Workflow zu verwenden, musst du nur eine Include-Anweisung vor dem Workflow-Block hinzufügen:

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

Wie du siehst, ist es ein ziemlich einfaches kleines Skript, das die Datei entpackt und zählt, wie viele Zeilen sie enthält.

??? info "Was macht `debug true`?"

    Die `debug true`-Direktive in der Prozessdefinition veranlasst Nextflow, die Ausgabe deines Skripts (wie die Zeilenanzahl "40") direkt im Ausführungslog auszugeben.
    Ohne dies würdest du nur den Prozessausführungsstatus sehen, aber nicht die tatsächliche Ausgabe deines Skripts.

    Für weitere Informationen zum Debuggen von Nextflow-Workflows siehe die Side Quest [Debugging Nextflow Workflows](debugging.md).

#### 1.4.2. Füge einen Aufruf von `COUNT_LINES` hinzu

Jetzt, da der Prozess für den Workflow verfügbar ist, können wir einen Aufruf des `COUNT_LINES`-Prozesses hinzufügen, um ihn auf der Eingabedatei auszuführen.

Nimm die folgenden Änderungen am Workflow vor:

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

        // Die Zeilen in der Datei zählen
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

Und jetzt führe den Workflow aus:

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

Dies zeigt, dass wir in der Lage sind, die Datei innerhalb eines Prozesses angemessen zu verarbeiten.

Konkret hat Nextflow die folgenden Operationen erfolgreich durchgeführt:

- Die Datei in das Arbeitsverzeichnis bereitgestellt
- Die .gz-Datei dekomprimiert
- Die Zeilen gezählt (40 Zeilen in diesem Fall)
- Ohne Fehler abgeschlossen

Der Schlüssel zu diesem reibungslosen Ablauf ist, dass wir Nextflow explizit mitteilen, dass unsere Eingabe eine Datei ist und als solche behandelt werden sollte.

### 1.5. Behebe grundlegende Datei-Eingabefehler

Dies bringt Nextflow-Neulinge oft durcheinander, also lass uns ein paar Minuten nehmen, um zu sehen, was passiert, wenn du es falsch machst.

Es gibt zwei Hauptstellen, an denen du die Dateiverarbeitung falsch machen kannst: auf der Ebene des Workflows und auf der Ebene des Prozesses.

#### 1.5.1. Fehler auf Workflow-Ebene

Lass uns sehen, was passiert, wenn wir zur Behandlung der Datei als String zurückkehren, wenn wir die Eingabe im Workflow-Block angeben.

Nimm die folgenden Änderungen am Workflow vor und stelle sicher, dass du die pfadspezifischen Print-Anweisungen auskommentierst:

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

        // Die Zeilen in der Datei zählen
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

        // Die Zeilen in der Datei zählen
        COUNT_LINES(myFile)
    ```

Und jetzt führe den Workflow aus:

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

Dies ist der wichtige Teil:

```console
Not a valid path value: 'data/patientA_rep1_normal_R1_001.fastq.gz'
```

Wenn du eine `path`-Eingabe angibst, validiert Nextflow, dass du tatsächliche Dateireferenzen übergibst, nicht nur Strings.
Dieser Fehler sagt dir, dass `'data/patientA_rep1_normal_R1_001.fastq.gz'` kein gültiger Pfadwert ist, weil es ein String ist, kein Path-Objekt.

Nextflow hat das Problem sofort erkannt und gestoppt, bevor der Prozess überhaupt gestartet wurde.

#### 1.5.2. Fehler auf Prozess-Ebene

Die andere Stelle, an der wir vergessen könnten anzugeben, dass Nextflow die Eingabe als Datei behandeln soll, ist in der Prozessdefinition.

!!! warning "Behalte den Workflow-Fehler aus 1.5.1"

    Damit dieser Test korrekt funktioniert, behalte den Workflow in seinem fehlerhaften Zustand (mit einem einfachen String anstelle von `file()`).
    In Kombination mit `val` im Prozess erzeugt dies den unten gezeigten Fehler.

Nimm die folgende Änderung am Modul vor:

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

Und jetzt führe den Workflow erneut aus:

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

Dies zeigt viele Details über den Fehler, weil der Prozess so eingestellt ist, dass er Debugging-Informationen ausgibt, wie oben erwähnt.

Dies sind die relevantesten Abschnitte:

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

Dies sagt, dass das System die Datei nicht finden konnte; wenn du jedoch den Pfad nachschlägst, gibt es eine Datei mit diesem Namen an diesem Ort.

Als wir dies ausführten, übergab Nextflow den String-Wert an das Skript, aber es _bereitstellte_ die tatsächliche Datei nicht im Arbeitsverzeichnis.
Der Prozess versuchte also, den relativen String `data/patientA_rep1_normal_R1_001.fastq.gz` zu verwenden, aber diese Datei existiert nicht innerhalb des Prozess-Arbeitsverzeichnisses.

Zusammengenommen zeigen diese beiden Beispiele, wie wichtig es ist, Nextflow mitzuteilen, ob eine Eingabe als Datei behandelt werden soll.

!!! note "Hinweis"

    Stelle sicher, dass du beide absichtlichen Fehler behebst, bevor du mit dem nächsten Abschnitt fortfährst.

### Fazit

- Pfad-Strings vs. Path-Objekte: Strings sind nur Text, Path-Objekte sind intelligente Dateireferenzen
- Die `file()`-Methode wandelt einen String-Pfad in ein Path-Objekt um, mit dem Nextflow arbeiten kann
- Du kannst auf Dateieigenschaften wie `name`, `simpleName`, `extension` und `parent` [mit Dateiattributen](https://www.nextflow.io/docs/latest/working-with-files.html#getting-file-attributes) zugreifen
- Die Verwendung von Path-Objekten anstelle von Strings ermöglicht es Nextflow, Dateien in deinem Workflow ordnungsgemäß zu verwalten
- Prozess-Eingabe-Ergebnisse: Eine ordnungsgemäße Dateiverarbeitung erfordert Path-Objekte, keine Strings, um sicherzustellen, dass Dateien korrekt bereitgestellt und für die Verwendung durch Prozesse zugänglich sind

---

## 2. Verwendung entfernter Dateien

Eine der Hauptfunktionen von Nextflow ist die Fähigkeit, nahtlos zwischen lokalen Dateien (auf demselben Rechner) und entfernten Dateien, die über das Internet zugänglich sind, zu wechseln.

Wenn du es richtig machst, solltest du niemals die Logik deines Workflows ändern müssen, um Dateien von verschiedenen Orten zu verarbeiten.
Alles, was du tun musst, um eine entfernte Datei zu verwenden, ist das entsprechende Präfix im Dateipfad anzugeben, wenn du ihn dem Workflow bereitstellst.

Zum Beispiel hat `/path/to/data` kein Präfix, was anzeigt, dass es ein 'normaler' lokaler Dateipfad ist, während `s3://path/to/data` das `s3://`-Präfix enthält, was anzeigt, dass es sich in Amazons S3-Objektspeicher befindet.

Viele verschiedene Protokolle werden unterstützt:

- HTTP(S)/FTP (http://, https://, ftp://)
- Amazon S3 (s3://)
- Azure Blob Storage (az://)
- Google Cloud Storage (gs://)

Um eines davon zu verwenden, gib einfach das relevante Präfix im String an, der dann technisch als Uniform Resource Identifier (URI) anstelle eines Dateipfads bezeichnet wird.
Nextflow übernimmt die Authentifizierung und Bereitstellung der Dateien am richtigen Ort, das Herunterladen oder Hochladen und alle anderen Dateioperationen, die du erwarten würdest.

Die Hauptstärke dieses Systems ist, dass es uns ermöglicht, zwischen Umgebungen zu wechseln, ohne die Pipeline-Logik zu ändern.
Du kannst beispielsweise mit einem kleinen, lokalen Testset entwickeln, bevor du zu einem vollständigen Testset wechselst, das sich in entferntem Speicher befindet, indem du einfach die URI änderst.

### 2.1. Verwende eine Datei aus dem Internet

Lass uns das testen, indem wir den lokalen Pfad, den wir unserem Workflow bereitstellen, durch einen HTTPS-Pfad ersetzen, der auf eine Kopie derselben Daten zeigt, die in Github gespeichert sind.

!!! warning "Warnung"

    Dies funktioniert nur, wenn du eine aktive Internetverbindung hast.

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

Lass uns den Workflow ausführen:

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

Es funktioniert! Du kannst sehen, dass sich sehr wenig geändert hat.

Der eine Unterschied in der Konsolenausgabe ist, dass die Pfadobjektklasse jetzt `nextflow.file.http.XPath` ist, während für den lokalen Pfad die Klasse `sun.nio.fs.UnixPath` war.
Du musst dir diese Klassen nicht merken; wir erwähnen dies nur, um zu demonstrieren, dass Nextflow die verschiedenen Speicherorte identifiziert und angemessen behandelt.

Hinter den Kulissen hat Nextflow die Datei in ein Staging-Verzeichnis heruntergeladen, das sich im Arbeitsverzeichnis befindet.
Diese bereitgestellte Datei kann dann als lokale Datei behandelt und in das relevante Prozessverzeichnis verlinkt werden.

Du kannst überprüfen, dass dies geschehen ist, indem du dir den Inhalt des Arbeitsverzeichnisses am Hash-Wert des Prozesses ansiehst.

??? abstract "Work-Verzeichnis Inhalt"

    Wenn der Prozess-Hash `8a/2ab7ca` war, könntest du das Arbeitsverzeichnis erkunden:

    ```console
    $ ls -la work/8a/2ab7ca*/
    total 16
    drwxr-xr-x  6 user  staff   192 Jan 28 10:00 .
    drwxr-xr-x  3 user  staff    96 Jan 28 10:00 ..
    -rw-r--r--  1 user  staff     0 Jan 28 10:00 .command.begin
    -rw-r--r--  1 user  staff   127 Jan 28 10:00 .command.sh
    lrwxr-xr-x  1 user  staff    89 Jan 28 10:00 patientA_rep1_normal_R1_001.fastq.gz -> /path/to/work/stage/.../patientA_rep1_normal_R1_001.fastq.gz
    ```

    Der Symlink zeigt auf eine bereitgestellte Kopie der entfernten Datei, die Nextflow automatisch heruntergeladen hat.

Beachte, dass bei größeren Dateien der Download-Schritt im Vergleich zur Ausführung auf lokalen Dateien zusätzliche Zeit in Anspruch nimmt.
Nextflow prüft jedoch, ob es bereits eine bereitgestellte Kopie hat, um unnötige Downloads zu vermeiden.
Wenn du also erneut auf derselben Datei ausführst und die bereitgestellte Datei nicht gelöscht hast, verwendet Nextflow die bereitgestellte Kopie.

Dies zeigt, wie einfach es ist, mit Nextflow zwischen lokalen und entfernten Daten zu wechseln, was eine Hauptfunktion von Nextflow ist.

!!! note "Hinweis"

    Die eine wichtige Ausnahme von diesem Prinzip ist, dass du keine Glob-Muster oder Verzeichnispfade mit HTTPS verwenden kannst, da HTTPS nicht mehrere Dateien auflisten kann, sodass du exakte Datei-URLs angeben musst.
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

    Wir zeigen dir im nächsten Abschnitt, wie du in der Praxis mit Globs arbeitest.

### 2.2. Wechsle zurück zur lokalen Datei

Wir werden für den Rest dieser Side Quest wieder unsere lokalen Beispieldateien verwenden, also lass uns die Workflow-Eingabe zurück zur ursprünglichen Datei wechseln:

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
- Nextflow lädt die Daten automatisch herunter und stellt sie am richtigen Ort bereit, solange diese Pfade an Prozesse übergeben werden
- Schreibe keine Logik zum Herunterladen oder Hochladen entfernter Dateien!
- Lokale und entfernte Dateien erzeugen unterschiedliche Objekttypen, funktionieren aber identisch
- **Wichtig**: HTTP/HTTPS funktioniert nur mit einzelnen Dateien (keine Glob-Muster)
- Cloud-Speicher (S3, Azure, GCS) unterstützt sowohl einzelne Dateien als auch Glob-Muster
- Du kannst nahtlos zwischen lokalen und entfernten Datenquellen wechseln, ohne die Code-Logik zu ändern (solange das Protokoll deine erforderlichen Operationen unterstützt)

---

## 3. Verwendung der `fromPath()`-Kanal-Factory

Bisher haben wir mit jeweils einer Datei gearbeitet, aber in Nextflow möchten wir typischerweise einen Eingabekanal mit mehreren zu verarbeitenden Eingabedateien erstellen.

Eine naive Methode wäre, die `file()`-Methode mit [`channel.of()`](https://www.nextflow.io/docs/latest/reference/channel.html#of) wie folgt zu kombinieren:

```groovy title="Syntax example"
ch_files = channel.of([file('data/patientA_rep1_normal_R1_001.fastq.gz')],
                      [file('data/patientA_rep1_normal_R1_001.fastq.gz')])
```

Das funktioniert, ist aber umständlich.

!!! tip "Wann `file()` vs. `channel.fromPath()` verwenden"

    - Verwende `file()`, wenn du ein einzelnes Path-Objekt für direkte Manipulation benötigst (Überprüfung, ob eine Datei existiert, Lesen ihrer Attribute oder Übergabe an einen einzelnen Prozessaufruf)
    - Verwende `channel.fromPath()`, wenn du einen Kanal benötigst, der mehrere Dateien enthalten kann, insbesondere mit Glob-Mustern, oder wenn Dateien durch mehrere Prozesse fließen

Hier kommt [`channel.fromPath()`](https://www.nextflow.io/docs/latest/reference/channel.html#frompath) ins Spiel: eine praktische Kanal-Factory, die alle Funktionen bündelt, die wir benötigen, um einen Kanal aus einem oder mehreren statischen Datei-Strings sowie Glob-Mustern zu generieren.

### 3.1. Füge die Kanal-Factory hinzu

Lass uns unseren Workflow aktualisieren, um `channel.fromPath` zu verwenden.

=== "Danach"

    ```groovy title="main.nf" linenums="7" hl_lines="1-3"
        // Dateien mit channel.fromPath laden
        ch_files = channel.fromPath('data/patientA_rep1_normal_R1_001.fastq.gz')
        ch_files.view { myFile -> "Found file: $myFile" }

        // Dateiattribute ausgeben
        /* Kommentiere diese vorerst aus, wir kommen darauf zurück!
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
        */

        // Die Zeilen in der Datei zählen
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

        // Die Zeilen in der Datei zählen
        COUNT_LINES(myFile)
    ```

Wir haben auch den Code auskommentiert, der die Attribute ausgibt, und stattdessen eine `.view`-Anweisung hinzugefügt, um nur den Dateinamen auszugeben.

Führe den Workflow aus:

```bash
nextflow run main.nf
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [grave_meucci] DSL2 - revision: b09964a583

    Found file: /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz
    ```

Wie du siehst, wird der Dateipfad als `Path`-Typ-Objekt im Kanal geladen.
Dies ist ähnlich wie bei `file()`, außer dass wir jetzt einen Kanal haben, in den wir bei Bedarf mehr Dateien laden können.

Die Verwendung von `channel.fromPath()` ist eine bequeme Möglichkeit, einen neuen Kanal zu erstellen, der mit einer Liste von Dateien gefüllt ist.

### 3.2. Zeige Attribute von Dateien im Kanal an

In unserem ersten Versuch mit der Kanal-Factory haben wir den Code vereinfacht und nur den Dateinamen ausgegeben.

Lass uns zurückgehen und die vollständigen Dateiattribute ausgeben:

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

        // Die Zeilen in der Datei zählen
        COUNT_LINES(ch_files)
    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="7" hl_lines="3"
        // Dateien mit channel.fromPath laden
        ch_files = channel.fromPath('data/patientA_rep1_normal_R1_001.fastq.gz')
        ch_files.view { myFile -> "Found file: $myFile" }

        // Die Zeilen in der Datei zählen
        // COUNT_LINES(ch_files)
    ```

Wir aktivieren auch den `COUNT_LINES`-Prozessaufruf wieder, um zu überprüfen, dass die Dateiverarbeitung mit unserem kanalbasierten Ansatz weiterhin korrekt funktioniert.

Führe den Workflow aus:

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

Und da hast du es, dieselben Ergebnisse wie zuvor, aber jetzt haben wir die Datei in einem Kanal, sodass wir mehr hinzufügen können.

### 3.3. Verwendung eines Globs zum Matchen mehrerer Dateien

Es gibt mehrere Möglichkeiten, mehr Dateien in den Kanal zu laden.
Hier zeigen wir dir, wie du Glob-Muster verwendest, die eine bequeme Möglichkeit sind, Datei- und Verzeichnisnamen basierend auf Wildcard-Zeichen zu matchen und abzurufen.
Der Prozess des Matchens dieser Muster wird "Globbing" oder "Filename Expansion" genannt.

!!! note "Hinweis"

    Wie bereits erwähnt, unterstützt Nextflow Globbing zur Verwaltung von Eingabe- und Ausgabedateien in den meisten Fällen, außer bei HTTPS-Dateipfaden, da HTTPS nicht mehrere Dateien auflisten kann.

Nehmen wir an, wir möchten beide Dateien in einem Paar von Dateien abrufen, die mit einer bestimmten Patient\*in, `patientA`, verbunden sind:

```console
patientA_rep1_normal_R1_001.fastq.gz
patientA_rep1_normal_R2_001.fastq.gz
```

Da der einzige Unterschied zwischen den Dateinamen die Replikatnummer ist, _d.h._ die Zahl nach `R`, können wir das Wildcard-Zeichen `*` wie folgt als Platzhalter für die Zahl verwenden:

```console
patientA_rep1_normal_R*_001.fastq.gz
```

Das ist das Glob-Muster, das wir benötigen.

Jetzt müssen wir nur noch den Dateipfad in der Kanal-Factory aktualisieren, um dieses Glob-Muster wie folgt zu verwenden:

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

Nextflow erkennt automatisch, dass dies ein Glob-Muster ist, und behandelt es entsprechend.

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

Mit dieser Methode können wir so viele oder so wenige Dateien abrufen, wie wir möchten, indem wir einfach das Glob-Muster ändern. Wenn wir es großzügiger gestalten würden, zum Beispiel indem wir alle variablen Teile der Dateinamen durch `*` ersetzen (_z.B._ `data/patient*_rep*_*_R*_001.fastq.gz`), könnten wir alle Beispieldateien im `data`-Verzeichnis erfassen.

### Fazit

- `channel.fromPath()` erstellt einen Kanal mit Dateien, die einem Muster entsprechen
- Jede Datei wird als separates Element im Kanal ausgegeben
- Wir können ein Glob-Muster verwenden, um mehrere Dateien zu matchen
- Dateien werden automatisch in Path-Objekte mit vollständigen Attributen umgewandelt
- Die `.view()`-Methode ermöglicht die Inspektion von Kanalinhalten

---

## 4. Extraktion grundlegender Metadaten aus Dateinamen

In den meisten wissenschaftlichen Bereichen ist es sehr üblich, Metadaten in den Namen der Dateien zu kodieren, die die Daten enthalten.
In der Bioinformatik werden beispielsweise Dateien mit Sequenzierungsdaten oft so benannt, dass sie Informationen über die Probe, Bedingung, Replikat und Read-Nummer kodieren.

Wenn die Dateinamen nach einer konsistenten Konvention konstruiert sind, kannst du diese Metadaten standardisiert extrahieren und im Verlauf deiner Analyse verwenden.
Das ist natürlich ein großes 'wenn', und du solltest sehr vorsichtig sein, wenn du dich auf die Dateinamenstruktur verlässt; aber die Realität ist, dass dieser Ansatz sehr weit verbreitet ist, also lass uns einen Blick darauf werfen, wie das in Nextflow gemacht wird.

Im Fall unserer Beispieldaten wissen wir, dass die Dateinamen konsistent strukturierte Metadaten enthalten.
Zum Beispiel kodiert der Dateiname `patientA_rep1_normal_R2_001` Folgendes:

- Patient\*innen-ID: `patientA`
- Replikat-ID: `rep1`
- Probentyp: `normal` (im Gegensatz zu `tumor`)
- Read-Set: `R1` (im Gegensatz zu `R2`)

Wir werden unseren Workflow in drei Schritten modifizieren, um diese Informationen abzurufen:

1. Den `simpleName` der Datei abrufen, der die Metadaten enthält
2. Die Metadaten mit einer Methode namens `tokenize()` trennen
3. Eine Map verwenden, um die Metadaten zu organisieren

!!! warning "Warnung"

    Du solltest niemals sensible Informationen in Dateinamen kodieren, wie Patient\*innennamen oder andere identifizierende Merkmale, da dies die Privatsphäre von Patient\*innen oder andere relevante Sicherheitsbeschränkungen gefährden kann.

### 4.1. Rufe den `simpleName` ab

Der `simpleName` ist ein Dateiattribut, das dem Dateinamen ohne Pfad und Erweiterung entspricht.

Nimm die folgenden Änderungen am Workflow vor:

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

Dies ruft den `simpleName` ab und verknüpft ihn mit dem vollständigen Dateiobjekt unter Verwendung einer `map()`-Operation.

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

### 4.2. Extrahiere die Metadaten aus dem `simplename`

An diesem Punkt sind die Metadaten, die wir wollen, im `simplename` eingebettet, aber wir können nicht direkt auf einzelne Elemente zugreifen.
Wir müssen also den `simplename` in seine Komponenten aufteilen.
Glücklicherweise sind diese Komponenten im ursprünglichen Dateinamen einfach durch Unterstriche getrennt, sodass wir eine gängige Nextflow-Methode namens `tokenize()` anwenden können, die perfekt für diese Aufgabe ist.

Nimm die folgenden Änderungen am Workflow vor:

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

Die `tokenize()`-Methode teilt den `simpleName`-String überall dort, wo sie Unterstriche findet, und gibt eine Liste zurück, die die Teilstrings enthält.

Führe den Workflow aus:

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

Jetzt enthält das Tupel für jedes Element in unserem Kanal die Liste der Metadaten (_z.B._ `[patientA, rep1, normal, R1, 001]`) und das ursprüngliche Dateiobjekt.

Das ist großartig!
Wir haben unsere Patient\*inneninformationen von einem einzelnen String in eine Liste von Strings aufgeteilt.
Wir können jetzt jeden Teil der Patient\*inneninformationen separat handhaben.

### 4.3. Verwende eine Map zur Organisation der Metadaten

Unsere Metadaten sind im Moment nur eine flache Liste.
Sie ist einfach genug zu verwenden, aber schwer zu lesen.

```console
[patientA, rep1, normal, R1, 001]
```

Was ist das Element an Index 3? Kannst du es sagen, ohne auf die ursprüngliche Erklärung der Metadatenstruktur zurückzugreifen?

Dies ist eine großartige Gelegenheit, einen Key-Value-Store zu verwenden, bei dem jedes Element eine Reihe von Schlüsseln und ihre zugehörigen Werte hat, sodass du leicht auf jeden Schlüssel verweisen kannst, um den entsprechenden Wert zu erhalten.

In unserem Beispiel bedeutet das, von dieser Organisation:

```groovy
data = [patientA, 1, normal, R1]

println data[3]
```

Zu dieser:

```groovy
data = [id: patientA, replicate: 1, type: normal, readNum: 1]

println data.readNum
```

In Nextflow wird das eine [Map](https://nextflow.io/docs/latest/script.html#maps) genannt.

Lass uns jetzt unsere flache Liste in eine Map umwandeln.
Nimm die folgenden Änderungen am Workflow vor:

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

- **Destrukturierende Zuweisung**: `def (patient, replicate, type, readNum) = ...` extrahiert die tokenisierten Werte in benannte Variablen in einer Zeile
- **Map-Literal-Syntax**: `[id: patient, replicate: ...]` erstellt eine Map, bei der jeder Schlüssel (wie `id`) mit einem Wert (wie `patient`) verknüpft ist
- **Verschachtelte Struktur**: Die äußere Liste `[..., myFile]` paart die Metadaten-Map mit dem ursprünglichen Dateiobjekt

Wir haben auch ein paar der Metadaten-Strings mit einer String-Ersetzungsmethode namens `replace()` vereinfacht, um einige unnötige Zeichen zu entfernen (_z.B._ `replicate.replace('rep', '')`, um nur die Zahl aus den Replikat-IDs zu behalten).

Lass uns den Workflow erneut ausführen:

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

Jetzt sind die Metadaten ordentlich beschriftet (_z.B._ `[id:patientA, replicate:1, type:normal, readNum:2]`), sodass es viel einfacher ist zu erkennen, was was ist.

Es wird auch viel einfacher sein, Elemente der Metadaten im Workflow tatsächlich zu verwenden, und wird unseren Code lesbarer und wartbarer machen.

### Fazit

- Wir können Dateinamen in Nextflow mit der Kraft einer vollständigen Programmiersprache handhaben
- Wir können die Dateinamen als Strings behandeln, um relevante Informationen zu extrahieren
- Die Verwendung von Methoden wie `tokenize()` und `replace()` ermöglicht es uns, Strings im Dateinamen zu manipulieren
- Die `.map()`-Operation transformiert Kanalelemente unter Beibehaltung der Struktur
- Strukturierte Metadaten (Maps) machen Code lesbarer und wartbarer als positionelle Listen

Als Nächstes werden wir uns ansehen, wie man gepaarte Datendateien handhabt.

---

## 5. Handhabung gepaarter Datendateien

Viele experimentelle Designs erzeugen gepaarte Datendateien, die davon profitieren, explizit gepaart behandelt zu werden.
In der Bioinformatik werden beispielsweise Sequenzierungsdaten oft in Form von gepaarten Reads erzeugt, was Sequenzstrings bedeutet, die vom selben DNA-Fragment stammen (oft 'Forward' und 'Reverse' genannt, weil sie von entgegengesetzten Enden gelesen werden).

Das ist bei unseren Beispieldaten der Fall, wo R1 und R2 sich auf die beiden Read-Sets beziehen.

```console
data/patientA_rep1_normal_R1_001.fastq.gz
data/patientA_rep1_normal_R2_001.fastq.gz
```

Nextflow bietet eine spezialisierte Kanal-Factory für die Arbeit mit gepaarten Dateien wie dieser namens `channel.fromFilePairs()`, die automatisch Dateien basierend auf einem gemeinsamen Benennungsmuster gruppiert. Das ermöglicht es dir, die gepaarten Dateien mit weniger Aufwand enger zu verknüpfen.

Wir werden unseren Workflow modifizieren, um dies zu nutzen.
Es wird zwei Schritte dauern:

1. Die Kanal-Factory auf `channel.fromFilePairs()` umstellen
2. Die Metadaten extrahieren und mappen

### 5.1. Stelle die Kanal-Factory auf `channel.fromFilePairs()` um

Um `channel.fromFilePairs` zu verwenden, müssen wir das Muster angeben, das Nextflow verwenden soll, um die beiden Mitglieder eines Paares zu identifizieren.

Zurück zu unseren Beispieldaten können wir das Benennungsmuster wie folgt formalisieren:

```console
data/patientA_rep1_normal_R{1,2}_001.fastq.gz
```

Dies ist ähnlich dem Glob-Muster, das wir zuvor verwendet haben, außer dass dies speziell die Teilstrings aufzählt (entweder `1` oder `2`, die direkt nach dem R kommen), die die beiden Mitglieder des Paares identifizieren.

Lass uns den Workflow `main.nf` entsprechend aktualisieren:

=== "Danach"

    ```groovy title="main.nf" linenums="7" hl_lines="1-2"
        // Dateien mit channel.fromFilePairs laden
        ch_files = channel.fromFilePairs('data/patientA_rep1_normal_R{1,2}_001.fastq.gz')
        /* Kommentiere das Mapping vorerst aus, wir kommen darauf zurück!
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

Wir haben die Kanal-Factory umgestellt und das Datei-Matching-Muster angepasst, und dabei haben wir die Map-Operation auskommentiert.
Wir werden das später mit ein paar Modifikationen wieder hinzufügen.

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

Oh-oh, diesmal ist die Ausführung fehlgeschlagen!

Der relevante Teil der Fehlermeldung ist hier:

```console
Not a valid path value: 'patientA_rep1_normal_R'
```

Das liegt daran, dass wir die Kanal-Factory geändert haben.
Bis jetzt enthielt der ursprüngliche Eingabekanal nur die Dateipfade.
Alle Metadaten-Manipulationen, die wir durchgeführt haben, haben den Kanalinhalt nicht tatsächlich beeinflusst.

Jetzt, da wir die `.fromFilePairs`-Kanal-Factory verwenden, ist der Inhalt des resultierenden Kanals anders.
Wir sehen nur ein Kanalelement, das aus einem Tupel besteht, das zwei Elemente enthält: den Teil des `simpleName`, der von den beiden Dateien geteilt wird und als Identifikator dient, und ein Tupel, das die beiden Dateiobjekte enthält, im Format `id, [ file1, file2 ]`.

Das ist großartig, weil Nextflow die harte Arbeit geleistet hat, den Patient\*innennamen zu extrahieren, indem es das gemeinsame Präfix untersucht und es als Patient\*innen-Identifikator verwendet hat.

Es bricht jedoch unseren aktuellen Workflow. Wenn wir `COUNT_LINES` immer noch auf die gleiche Weise ausführen wollten, ohne den Prozess zu ändern, müssten wir eine Mapping-Operation anwenden, um die Dateipfade zu extrahieren.
Aber das werden wir nicht tun, weil unser ultimatives Ziel ist, einen anderen Prozess, `ANALYZE_READS`, zu verwenden, der Dateipaare angemessen handhabt.

Lass uns also einfach den Aufruf von `COUNT_LINES` auskommentieren (oder löschen) und weitermachen.

=== "Danach"

    ```groovy title="main.nf" linenums="26" hl_lines="2"
        // Die Zeilen in der Datei zählen
        // COUNT_LINES(ch_files)
    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="26" hl_lines="2"
        // Die Zeilen in der Datei zählen
        COUNT_LINES(ch_files)
    ```

Du kannst auch die `COUNT_LINES`-Include-Anweisung auskommentieren oder löschen, aber das wird keine funktionale Auswirkung haben.

Lass uns jetzt den Workflow erneut ausführen:

```bash
nextflow run main.nf
```

??? success "Befehlsausgabe"

    ```console hl_lines="5"
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [fabulous_davinci] DSL2 - revision: 22b53268dc

    [patientA_rep1_normal_R, [/workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz, /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R2_001.fastq.gz]]
    ```

Juhu, diesmal gelingt der Workflow!

Wir müssen jedoch noch den Rest der Metadaten aus dem `id`-Feld herausholen.

### 5.2. Extrahiere und organisiere Metadaten aus Dateipaaren

Unsere `map`-Operation von vorher funktioniert nicht, weil sie nicht zur Datenstruktur passt, aber wir können sie so modifizieren, dass sie funktioniert.

Wir haben bereits Zugriff auf den tatsächlichen Patient\*innen-Identifikator in dem String, den `fromFilePairs()` als Identifikator verwendet hat, sodass wir diesen verwenden können, um die Metadaten zu extrahieren, ohne den `simpleName` vom Path-Objekt zu holen, wie wir es zuvor getan haben.

Kommentiere die Map-Operation im Workflow aus und nimm die folgenden Änderungen vor:

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
        /* Kommentiere das Mapping vorerst aus, wir kommen darauf zurück!
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

Diesmal beginnt die Map mit `id, files` anstelle von nur `myFile`, und `tokenize()` wird auf `id` anstelle von `myFile.simpleName` angewendet.

Beachte auch, dass wir `readNum` aus der `tokenize()`-Zeile entfernt haben; alle Teilstrings, die wir nicht speziell benennen (von links beginnend), werden stillschweigend verworfen.
Wir können dies tun, weil die gepaarten Dateien jetzt eng verknüpft sind, sodass wir `readNum` nicht mehr in der Metadaten-Map benötigen.

Lass uns den Workflow ausführen:

```bash
nextflow run main.nf
```

??? success "Befehlsausgabe"

    ```console

    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [prickly_stonebraker] DSL2 - revision: f62ab10a3f

    [[id:patientA, replicate:1, type:normal], [/workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz, /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R2_001.fastq.gz]]
    ```

Und da ist es: Wir haben die Metadaten-Map (`[id:patientA, replicate:1, type:normal]`) an der ersten Position des Ausgabe-Tupels, gefolgt vom Tupel der gepaarten Dateien, wie beabsichtigt.

Natürlich wird dies nur dieses spezifische Dateipaar erfassen und verarbeiten.
Wenn du mit der Verarbeitung mehrerer Paare experimentieren möchtest, kannst du versuchen, Wildcards in das Eingabemuster einzufügen und sehen, was passiert.
Versuche zum Beispiel `data/patientA_rep1_*_R{1,2}_001.fastq.gz` zu verwenden.

### Fazit

- [`channel.fromFilePairs()` findet und paart automatisch zusammengehörige Dateien](https://www.nextflow.io/docs/latest/reference/channel.html#fromfilepairs)
- Dies vereinfacht die Handhabung von Paired-End-Reads in deiner Pipeline
- Gepaarte Dateien können als `[id, [file1, file2]]`-Tupel gruppiert werden
- Die Metadatenextraktion kann aus der gepaarten Datei-ID anstelle von einzelnen Dateien erfolgen

---

## 6. Verwendung von Dateioperationen in Prozessen

Lass uns jetzt all dies in einem einfachen Prozess zusammenfügen, um zu verstärken, wie man Dateioperationen innerhalb eines Nextflow-Prozesses verwendet.

Wir stellen dir ein vorgefertigtes Prozessmodul namens `ANALYZE_READS` zur Verfügung, das ein Tupel aus Metadaten und ein Paar Eingabedateien nimmt und sie analysiert.
Wir könnten uns vorstellen, dass dies Sequenz-Alignment, Variant Calling oder jeden anderen Schritt durchführt, der für diesen Datentyp sinnvoll ist.

Lass uns anfangen.

### 6.1. Importiere den Prozess und untersuche den Code

Um diesen Prozess im Workflow zu verwenden, müssen wir nur eine Modul-Include-Anweisung vor dem Workflow-Block hinzufügen.

Nimm die folgende Änderung am Workflow vor:

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

    Die `tag`- und `publishDir`-Direktiven verwenden Closure-Syntax (`{ ... }`) anstelle von String-Interpolation (`"${...}"`).
    Dies liegt daran, dass diese Direktiven auf Eingabevariablen (`meta`) verweisen, die erst zur Laufzeit verfügbar sind.
    Die Closure-Syntax verschiebt die Auswertung, bis der Prozess tatsächlich läuft.

!!! note "Hinweis"

    Wir nennen unsere Metadaten-Map `meta` nach Konvention.
    Für einen tieferen Einblick in Meta-Maps siehe die Side Quest [Metadata and meta maps](./metadata.md).

### 6.2. Rufe den Prozess im Workflow auf

Jetzt, da der Prozess für den Workflow verfügbar ist, können wir einen Aufruf des `ANALYZE_READS`-Prozesses hinzufügen, um ihn auszuführen.

Um ihn auf unseren Beispieldaten auszuführen, müssen wir zwei Dinge tun:

1. Dem remappten Kanal einen Namen geben
2. Einen Aufruf des Prozesses hinzufügen

#### 6.2.1. Benenne den remappten Eingabekanal

Wir haben die Mapping-Manipulationen zuvor direkt auf den Eingabekanal angewendet.
Um die remappten Inhalte an den `ANALYZE_READS`-Prozess zu übergeben (und dies auf eine klare und leicht lesbare Weise zu tun), möchten wir einen neuen Kanal namens `ch_samples` erstellen.

Wir können das mit dem [`set`](https://www.nextflow.io/docs/latest/reference/operator.html#set)-Operator tun.

Ersetze im Haupt-Workflow den `.view()`-Operator durch `.set { ch_samples }` und füge eine Zeile hinzu, die testet, dass wir auf den Kanal mit Namen verweisen können.

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

Lass uns das ausführen:

```bash
nextflow run main.nf
```

??? success "Befehlsausgabe"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [goofy_kirch] DSL2 - revision: 3313283e42

    [[id:patientA, replicate:1, type:normal], [/workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz, /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R2_001.fastq.gz]]
    ```

Dies bestätigt, dass wir jetzt auf den Kanal mit Namen verweisen können.

#### 6.2.2. Rufe den Prozess auf den Daten auf

Lass uns jetzt tatsächlich den `ANALYZE_READS`-Prozess auf dem `ch_samples`-Kanal aufrufen.

Nimm im Haupt-Workflow die folgenden Code-Änderungen vor:

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

Lass uns das ausführen:

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

Der Prozess hat unsere Eingaben genommen und eine neue Datei erstellt, die die Patient\*innen-Metadaten enthält, wie vorgesehen.
Großartig!

### 6.3. Beziehe viele weitere Patient\*innen ein

Natürlich verarbeitet dies nur ein einzelnes Dateipaar für eine\*n einzelne\*n Patient\*in, was nicht gerade der hohe Durchsatz ist, den du mit Nextflow erhoffst.
Du möchtest wahrscheinlich viel mehr Daten auf einmal verarbeiten.

Denk daran, dass `channel.fromPath()` einen _Glob_ als Eingabe akzeptiert, was bedeutet, dass es eine beliebige Anzahl von Dateien akzeptieren kann, die dem Muster entsprechen.
Wenn wir also alle Patient\*innen einbeziehen möchten, können wir einfach den Eingabe-String modifizieren, um mehr Patient\*innen einzubeziehen, wie zuvor beiläufig erwähnt.

Nehmen wir an, wir möchten so gierig wie möglich sein.
Nimm die folgenden Änderungen am Workflow vor:

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

Das Ergebnisverzeichnis sollte jetzt Ergebnisse für alle verfügbaren Daten enthalten.

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

Erfolg! Wir haben alle Patient\*innen in einem Durchgang analysiert! Richtig?

Vielleicht nicht.
Wenn du genauer hinschaust, haben wir ein Problem: Wir haben zwei Replikate für patientA, aber nur eine Ausgabedatei!
Wir überschreiben die Ausgabedatei jedes Mal.

### 6.4. Mache die veröffentlichten Dateien eindeutig

Da wir Zugriff auf die Patient\*innen-Metadaten haben, können wir sie verwenden, um die veröffentlichten Dateien eindeutig zu machen, indem wir differenzierende Metadaten entweder in der Verzeichnisstruktur oder in den Dateinamen selbst einbeziehen.

Nimm die folgende Änderung am Workflow vor:

=== "Danach"

    ```groovy title="modules/analyze_reads.nf" linenums="6"
        publishDir { "results/${meta.type}/${meta.id}/${meta.replicate}" }, mode: 'copy'
    ```

=== "Vorher"

    ```groovy title="modules/analyze_reads.nf" linenums="6"
        publishDir { "results/${meta.id}" }, mode: 'copy'
    ```

Hier zeigen wir die Option, zusätzliche Verzeichnisebenen zu verwenden, um Probentypen und Replikate zu berücksichtigen, aber du könntest auch experimentieren, dies auf Dateinamenebene zu tun.

Führe die Pipeline jetzt ein letztes Mal aus, aber stelle sicher, dass du zuerst das Ergebnisverzeichnis entfernst, um dir einen sauberen Arbeitsbereich zu geben:

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

Überprüfe jetzt das Ergebnisverzeichnis:

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

Und da ist es, alle unsere Metadaten, ordentlich organisiert. Das ist Erfolg!

Es gibt noch viel mehr, was du tun kannst, sobald du deine Metadaten in eine Map wie diese geladen hast:

1. Organisierte Ausgabeverzeichnisse basierend auf Patient\*innen-Attributen erstellen
2. Entscheidungen in Prozessen basierend auf Patient\*innen-Eigenschaften treffen
3. Daten basierend auf Metadatenwerten aufteilen, verbinden und neu kombinieren

Dieses Muster, Metadaten explizit zu halten und an die Daten anzuhängen (anstatt sie in Dateinamen zu kodieren), ist eine zentrale Best Practice in Nextflow, die den Aufbau robuster, wartbarer Analyse-Workflows ermöglicht.
Du kannst mehr darüber in der Side Quest [Metadata and meta maps](./metadata.md) erfahren.

### Fazit

- Die `publishDir`-Direktive kann Ausgaben basierend auf Metadatenwerten organisieren
- Metadaten in Tupeln ermöglichen eine strukturierte Organisation von Ergebnissen
- Dieser Ansatz schafft wartbare Workflows mit klarer Datenherkunft
- Prozesse können Tupel aus Metadaten und Dateien als Eingabe nehmen
- Die `tag`-Direktive bietet Prozessidentifikation in Ausführungslogs
- Die Workflow-Struktur trennt Kanalerstellung von Prozessausführung

---

## Zusammenfassung

In dieser Side Quest hast du gelernt, wie man mit Dateien in Nextflow arbeitet, von grundlegenden Operationen bis hin zu fortgeschritteneren Techniken für die Handhabung von Dateisammlungen.

Die Anwendung dieser Techniken in deiner eigenen Arbeit wird es dir ermöglichen, effizientere und wartbarere Workflows zu erstellen, insbesondere bei der Arbeit mit großen Mengen von Dateien mit komplexen Benennungskonventionen.

### Wichtige Muster

1.  **Grundlegende Dateioperationen:** Wir haben Path-Objekte mit `file()` erstellt und auf Dateiattribute wie Name, Erweiterung und übergeordnetes Verzeichnis zugegriffen, wobei wir den Unterschied zwischen Strings und Path-Objekten gelernt haben.

    - Ein Path-Objekt erstellen mit `file()`

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

2.  **Verwendung entfernter Dateien**: Wir haben gelernt, wie man transparent zwischen lokalen und entfernten Dateien mit URIs wechselt, was Nextflows Fähigkeit demonstriert, Dateien aus verschiedenen Quellen zu handhaben, ohne die Workflow-Logik zu ändern.

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

3.  **Laden von Dateien mit der `fromPath()`-Kanal-Factory:** Wir haben Kanäle aus Dateimustern mit `channel.fromPath()` erstellt und ihre Dateiattribute einschließlich Objekttypen angezeigt.

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

4.  **Extraktion von Patient\*innen-Metadaten aus Dateinamen:** Wir haben `tokenize()` und `replace()` verwendet, um Metadaten aus Dateinamen zu extrahieren und zu strukturieren und sie in organisierte Maps umzuwandeln.

    ```groovy
    def name = file.name.tokenize('_')
    def patientId = name[0]
    def replicate = name[1].replace('rep', '')
    def type = name[2]
    def readNum = name[3].replace('R', '')
    ```

5.  **Vereinfachung mit channel.fromFilePairs:** Wir haben `channel.fromFilePairs()` verwendet, um automatisch zusammengehörige Dateien zu paaren und Metadaten aus gepaarten Datei-IDs zu extrahieren.

    ```groovy
    ch_pairs = channel.fromFilePairs('data/*_R{1,2}_001.fastq.gz')
    ```

6.  **Verwendung von Dateioperationen in Prozessen:** Wir haben Dateioperationen mit korrekter Eingabeverarbeitung in Nextflow-Prozesse integriert und `publishDir` verwendet, um Ausgaben basierend auf Metadaten zu organisieren.

    - Eine Meta-Map mit den Prozess-Eingaben verknüpfen

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

### Zusätzliche Ressourcen

- [Nextflow-Dokumentation: Working with Files](https://www.nextflow.io/docs/latest/working-with-files.html)
- [channel.fromPath](https://www.nextflow.io/docs/latest/reference/channel.html#frompath)
- [channel.fromFilePairs](https://www.nextflow.io/docs/latest/reference/channel.html#fromfilepairs)

---

## Wie geht es weiter?

Kehre zum [Menü der Side Quests](./index.md) zurück oder klicke auf den Button unten rechts auf der Seite, um zum nächsten Thema in der Liste zu gelangen.
