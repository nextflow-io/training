# Verarbeitung von Datei-Eingaben

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Wissenschaftliche Analyse-Workflows beinhalten oft die Verarbeitung großer Mengen von Dateien.
Nextflow bietet leistungsstarke Werkzeuge, um Dateien effizient zu verarbeiten und hilft dir, deine Daten mit minimalem Code zu organisieren und zu verarbeiten.

### Lernziele

In dieser Side Quest werden wir erkunden, wie Nextflow mit Dateien umgeht, von grundlegenden Dateioperationen bis hin zu fortgeschritteneren Techniken für die Arbeit mit Dateisammlungen.
Du wirst lernen, wie man Metadaten aus Dateinamen extrahiert, was eine häufige Anforderung in wissenschaftlichen Analyse-Pipelines ist.

Am Ende dieser Side Quest kannst du:

- Path-Objekte aus Dateipfad-Strings mit der `file()`-Methode von Nextflow zu erstellen
- Auf Dateiattribute wie Name, Extension und übergeordnetes Verzeichnis zuzugreifen
- Sowohl lokale als auch entfernte Dateien transparent über URIs zu verarbeiten
- Channels zur Automatisierung der Dateiverarbeitung mit `channel.fromPath()` und `channel.fromFilePairs()` zu verwenden
- Metadaten aus Dateinamen mittels String-Manipulation zu extrahieren und zu strukturieren
- Zusammengehörige Dateien mittels Pattern Matching und Glob-Ausdrücken zu gruppieren
- Dateioperationen in Nextflow-Prozesse mit korrekter Input-Verarbeitung zu integrieren
- Prozess-Ausgaben mittels metadatengesteuerter Verzeichnisstrukturen zu organisieren

Diese Fähigkeiten helfen dir, Workflows zu erstellen, die verschiedene Arten von Datei-Eingaben mit großer Flexibilität verarbeiten können.

### Voraussetzungen

Bevor du diese Side Quest angehst, solltest du:

- Das [Hello Nextflow](../../hello_nextflow/)-Tutorial oder einen gleichwertigen Einsteigerkurs abgeschlossen haben.
- Mit grundlegenden Nextflow-Konzepten und -Mechanismen vertraut sein (Prozesse, Channels, Operatoren)

<!-- I removed the suggestion to do the metamaps SQ first because that works more naturally after -->

---

## 0. Erste Schritte

#### Öffne die Trainingsumgebung im Codespace

Falls du das noch nicht getan hast, stelle sicher, dass du die Trainingsumgebung wie in der [Umgebungseinrichtung](../envsetup/index.md) beschrieben öffnest.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Wechsle in das Projektverzeichnis

Lass uns in das Verzeichnis wechseln, in dem sich die Dateien für dieses Tutorial befinden.

```bash
cd side-quests/working_with_files
```

Du kannst VSCode so einstellen, dass es sich auf dieses Verzeichnis fokussiert:

```bash
code .
```

#### Überprüfe die Materialien

Du findest eine einfache Workflow-Datei namens `main.nf`, ein `modules`-Verzeichnis mit zwei Moduldateien und ein `data`-Verzeichnis mit einigen Beispieldateien.

??? abstract "Verzeichnisinhalte"

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

Dieses Verzeichnis enthält Paired-End-Sequenzierungsdaten von drei Patienten (A, B, C).

Für jeden Patienten haben wir Proben vom Typ `tumor` (typischerweise aus Tumorbiopsien) oder `normal` (aus gesundem Gewebe oder Blut).
Falls du mit Krebsanalysen nicht vertraut bist: Dies entspricht einem experimentellen Modell, das gepaarte Tumor/Normal-Proben für kontrastive Analysen verwendet.

Speziell für Patient A haben wir zwei Sets technischer Replikate (Wiederholungen).

Die Sequenzierungsdateien sind mit der typischen `_R1_`- und `_R2_`-Konvention für sogenannte 'Forward Reads' und 'Reverse Reads' benannt.

_Mach dir keine Sorgen, wenn du mit diesem experimentellen Design nicht vertraut bist, es ist nicht kritisch für das Verständnis dieses Tutorials._

#### Überprüfe die Aufgabenstellung

Deine Herausforderung besteht darin, einen Nextflow-Workflow zu schreiben, der:

1. Eingabedateien mit Nextflows Dateiverarbeitungsmethoden **lädt**
2. Metadaten (Patienten-ID, Replikat, Probentyp) aus der Dateinamenstruktur **extrahiert**
3. Gepaarte Dateien (R1/R2) mit `channel.fromFilePairs()` **gruppiert**
4. Die Dateien mit einem bereitgestellten Analysemodul **verarbeitet**
5. Ausgaben in einer Verzeichnisstruktur basierend auf den extrahierten Metadaten **organisiert**

#### Checkliste zur Bereitschaft

Denkst du, du bist bereit einzutauchen?

- [ ] Ich verstehe das Ziel dieses Kurses und seine Voraussetzungen
- [ ] Mein Codespace läuft
- [ ] Ich habe mein Arbeitsverzeichnis korrekt eingestellt
- [ ] Ich verstehe die Aufgabenstellung

Wenn du alle Kästchen abhaken kannst, bist du startklar.

---

## 1. Grundlegende Dateioperationen

### 1.1. Identifiziere den Typ eines Objekts mit `.class`

Wirf einen Blick auf die Workflow-Datei `main.nf`:

```groovy title="main.nf" linenums="1"
#!/usr/bin/env nextflow

workflow {

    // Create a Path object from a string path
    myFile = 'data/patientA_rep1_normal_R1_001.fastq.gz'

    println "${myFile} is of class ${myFile.class}"
}
```

Dies ist ein Mini-Workflow (ohne Prozesse), der auf einen einzelnen Dateipfad in seinem workflow verweist und ihn dann zusammen mit seiner Klasse auf der Konsole ausgibt.

??? info "Was ist `.class`?"

    In Nextflow sagt uns `.class`, mit welchem Objekttyp wir arbeiten. Es ist wie die Frage "Was für eine Art von Ding ist das?" zu stellen, um herauszufinden, ob es ein String, eine Zahl, eine Datei oder etwas anderes ist.
    Dies hilft uns, den Unterschied zwischen einem einfachen String und einem Path-Objekt in den nächsten Abschnitten zu illustrieren.

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

Dies ist nur eine Textausgabe; Nextflow hat noch nichts Besonderes damit gemacht.
Wir haben auch bestätigt, dass dies aus Sicht von Nextflow nur ein String (der Klasse `java.lang.String`) ist.
Das macht Sinn, da wir Nextflow noch nicht mitgeteilt haben, dass es sich um eine Datei handelt.

### 1.2. Erstelle ein Path-Objekt mit file()

Wir können Nextflow mitteilen, wie Dateien zu verarbeiten sind, indem wir [Path-Objekte](https://www.nextflow.io/docs/latest/reference/stdlib-types.html#path) aus Pfad-Strings erstellen.

In unserem Workflow können wir den String-Pfad `data/patientA_rep1_normal_R1_001.fastq.gz` mit der `file()`-Methode in ein Path-Objekt umwandeln, das Zugriff auf Dateieigenschaften und -operationen bietet.

Bearbeite `main.nf`, um den String wie folgt mit `file()` zu umschließen:

=== "Nachher"

    ```groovy title="main.nf" linenums="5" hl_lines="2"
        // Create a Path object from a string path
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        println "${myFile} is of class ${myFile.class}"
    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="5" hl_lines="2"
        // Create a Path object from a string path
        myFile = 'data/patientA_rep1_normal_R1_001.fastq.gz'

        println "${myFile} is of class ${myFile.class}"
    ```

Führe jetzt den Workflow erneut aus:

```bash
nextflow run main.nf
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [kickass_coulomb] DSL2 - revision: 5af44b1b59

    /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz is of class class sun.nio.fs.UnixPath
    ```

Dieses Mal siehst du den vollständigen absoluten Pfad statt des relativen Pfads, den wir als Eingabe bereitgestellt haben.

Nextflow hat unseren String in ein Path-Objekt umgewandelt und ihn zur tatsächlichen Dateiposition im System aufgelöst.
Der Dateipfad wird nun absolut sein, wie in `/workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz`.

Beachte auch, dass die Path-Objekt-Klasse `sun.nio.fs.UnixPath` ist: Dies ist Nextflows Art, lokale Dateien zu repräsentieren.
Wie wir später sehen werden, haben entfernte Dateien andere Klassennamen (wie `nextflow.file.http.XPath` für HTTP-Dateien), aber sie funktionieren alle genau gleich und können identisch in deinen Workflows verwendet werden.

!!! tip

    **Der Hauptunterschied:**

    - **Pfad-String**: Nur Text, den Nextflow als Zeichen behandelt
    - **Path-Objekt**: Eine intelligente Dateireferenz, mit der Nextflow arbeiten kann

    Stell es dir so vor: Ein Pfad-String ist wie eine auf Papier geschriebene Adresse, während ein Path-Objekt wie eine im GPS-Gerät geladene Adresse ist, die weiß, wie man dorthin navigiert und dir Details über die Route mitteilen kann.

### 1.3. Zugriff auf Dateiattribute

Warum ist das hilfreich? Nun, da Nextflow jetzt versteht, dass `myFile` ein Path-Objekt und nicht nur ein String ist, können wir auf die verschiedenen Attribute des Path-Objekts zugreifen.

Lass uns unseren Workflow aktualisieren, um die eingebauten Dateiattribute auszugeben:

=== "Nachher"

    ```groovy title="main.nf" linenums="5" hl_lines="4-9"
        // Create a Path object from a string path
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Print file attributes
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="5" hl_lines="4"
        // Create a Path object from a string path
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

Du siehst die verschiedenen Dateiattribute, die oben auf der Konsole ausgegeben werden.

### 1.4. Die Datei an einen Prozess übergeben

Der Unterschied zwischen Strings und Path-Objekten wird kritisch, wenn du anfängst, tatsächliche Workflows mit Prozessen zu erstellen.
Bisher haben wir überprüft, dass Nextflow unsere Eingabedatei nun als Datei behandelt, aber lass uns sehen, ob wir tatsächlich etwas mit dieser Datei in einem Prozess ausführen können.

#### 1.4.1. Importiere den Prozess und untersuche den Code

Wir stellen dir ein vorgefertigtes Prozessmodul namens `COUNT_LINES` zur Verfügung, das eine Dateieingabe nimmt und zählt, wie viele Zeilen sie enthält.

Um den Prozess im Workflow zu verwenden, musst du nur eine include-Anweisung vor dem workflow-Block hinzufügen:

=== "Nachher"

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

Wie du siehst, ist es ein ziemlich einfaches kleines Script, das die Datei entpackt und zählt, wie viele Zeilen sie enthält.

??? info "Was macht `debug true`?"

    Die `debug true`-Direktive in der Prozessdefinition veranlasst Nextflow, die Ausgabe deines Scripts (wie die Zeilenanzahl "40") direkt im Ausführungslog auszugeben.
    Ohne dies würdest du nur den Prozessausführungsstatus sehen, aber nicht die tatsächliche Ausgabe deines Scripts.

    Für weitere Informationen zum Debuggen von Nextflow-Workflows siehe die Side Quest [Debugging Nextflow Workflows](debugging.md).

#### 1.4.2. Füge einen Aufruf von `COUNT_LINES` hinzu

Nun, da der Prozess für den Workflow verfügbar ist, können wir einen Aufruf des `COUNT_LINES`-Prozesses hinzufügen, um ihn auf der Eingabedatei auszuführen.

Nimm folgende Änderungen am Workflow vor:

=== "Nachher"

    ```groovy title="main.nf" linenums="7" hl_lines="11-12"
        // Create a Path object from a string path
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Print file attributes
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"

        // Count the lines in the file
        COUNT_LINES(myFile)
    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="7" hl_lines="4-9"
        // Create a Path object from a string path
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Print file attributes
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

Konkret hat Nextflow folgende Operationen erfolgreich durchgeführt:

- Die Datei ins Arbeitsverzeichnis gestaged
- Die .gz-Datei dekomprimiert
- Die Zeilen gezählt (40 Zeilen in diesem Fall)
- Ohne Fehler abgeschlossen

Der Schlüssel zu dieser reibungslosen Operation ist, dass wir Nextflow explizit mitteilen, dass unsere Eingabe eine Datei ist und als solche behandelt werden soll.

### 1.5. Fehlerbehebung bei grundlegenden Dateieingabe-Fehlern

Dies führt oft zu Verwirrung bei Nextflow-Neulingen, also nehmen wir uns ein paar Minuten Zeit, um zu schauen, was passiert, wenn du es falsch machst.

Es gibt zwei Hauptstellen, an denen du die Dateiverarbeitung falsch machen kannst: auf Workflow-Ebene und auf Prozess-Ebene.

#### 1.5.1. Fehler auf Workflow-Ebene

Lass uns sehen, was passiert, wenn wir zur Behandlung der Datei als String zurückkehren, wenn wir die Eingabe im workflow-Block angeben.

Nimm folgende Änderungen am Workflow vor und stelle sicher, dass du die pfadspezifischen print-Anweisungen auskommentierst:

=== "Nachher"

    ```groovy title="main.nf" linenums="7" hl_lines="2 6-11"
        // Create a Path object from a string path
        myFile = 'data/patientA_rep1_normal_R1_001.fastq.gz'

        // Print file attributes
        println "File object class: ${myFile.class}"
        /*
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
        */

        // Count the lines in the file
        COUNT_LINES(myFile)
    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="7" hl_lines="4-9"
        // Create a Path object from a string path
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Print file attributes
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"

        // Count the lines in the file
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

Das ist der wichtige Teil:

```console
Not a valid path value: 'data/patientA_rep1_normal_R1_001.fastq.gz'
```

Wenn du eine `path`-Eingabe spezifizierst, validiert Nextflow, dass du tatsächliche Dateireferenzen übergibst, nicht nur Strings.
Dieser Fehler teilt dir mit, dass `'data/patientA_rep1_normal_R1_001.fastq.gz'` kein gültiger Pfadwert ist, weil es ein String ist, kein Path-Objekt.

Nextflow hat das Problem sofort erkannt und gestoppt, bevor der Prozess überhaupt gestartet wurde.

#### 1.5.2. Fehler auf Prozess-Ebene

Die andere Stelle, an der wir vergessen könnten zu spezifizieren, dass Nextflow die Eingabe als Datei behandeln soll, ist in der Prozessdefinition.

!!! warning "Behalte den Workflow-Fehler aus 1.5.1"

    Damit dieser Test korrekt funktioniert, belasse den Workflow in seinem fehlerhaften Zustand (mit einem einfachen String statt `file()`).
    In Kombination mit `val` im Prozess erzeugt dies den unten gezeigten Fehler.

Nimm folgende Änderung am Modul vor:

=== "Nachher"

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

Dies zeigt viele Details über den Fehler, da der Prozess so eingestellt ist, dass er Debugging-Informationen ausgibt, wie oben erwähnt.

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

Dies besagt, dass das System die Datei nicht finden konnte; wenn du jedoch den Pfad nachschlägst, gibt es dort eine Datei mit diesem Namen.

Als wir dies ausführten, hat Nextflow den String-Wert an das Script weitergegeben, aber es hat die tatsächliche Datei nicht ins Arbeitsverzeichnis _gestaged_.
Also versuchte der Prozess, den relativen String `data/patientA_rep1_normal_R1_001.fastq.gz` zu verwenden, aber diese Datei existiert nicht im Prozess-Arbeitsverzeichnis.

Zusammengenommen zeigen diese beiden Beispiele, wie wichtig es ist, Nextflow mitzuteilen, ob eine Eingabe als Datei verarbeitet werden soll.

!!! note

    Stelle sicher, dass du beide absichtlichen Fehler korrigierst, bevor du mit dem nächsten Abschnitt fortfährst.

### Zusammenfassung

- Pfad-Strings vs. Path-Objekte: Strings sind nur Text, Path-Objekte sind intelligente Dateireferenzen
- Die `file()`-Methode wandelt einen String-Pfad in ein Path-Objekt um, mit dem Nextflow arbeiten kann
- Du kannst auf Dateieigenschaften wie `name`, `simpleName`, `extension` und `parent` [über Dateiattribute](https://www.nextflow.io/docs/latest/working-with-files.html#getting-file-attributes) zugreifen
- Die Verwendung von Path-Objekten anstelle von Strings ermöglicht es Nextflow, Dateien in deinem Workflow ordnungsgemäß zu verwalten
- Prozess-Input-Ergebnisse: Eine ordnungsgemäße Dateiverarbeitung erfordert Path-Objekte, keine Strings, um sicherzustellen, dass Dateien korrekt gestaged und für die Verwendung durch Prozesse zugänglich sind.

---

## 2. Verwendung entfernter Dateien

Eines der Hauptmerkmale von Nextflow ist die Fähigkeit, nahtlos zwischen lokalen Dateien (auf demselben Rechner) und entfernten, über das Internet zugänglichen Dateien zu wechseln.

Wenn du es richtig machst, solltest du nie die Logik deines Workflows ändern müssen, um Dateien aus verschiedenen Speicherorten zu verarbeiten.
Alles, was du tun musst, um eine entfernte Datei zu verwenden, ist das entsprechende Präfix im Dateipfad anzugeben, wenn du ihn dem Workflow bereitstellst.

Zum Beispiel hat `/path/to/data` kein Präfix, was darauf hinweist, dass es sich um einen 'normalen' lokalen Dateipfad handelt, während `s3://path/to/data` das `s3://`-Präfix enthält, was darauf hinweist, dass es sich in Amazons S3-Objektspeicher befindet.

Viele verschiedene Protokolle werden unterstützt:

- HTTP(S)/FTP (http://, https://, ftp://)
- Amazon S3 (s3://)
- Azure Blob Storage (az://)
- Google Cloud Storage (gs://)

Um eines davon zu verwenden, gib einfach das entsprechende Präfix im String an, der dann technisch als Uniform Resource Identifier (URI) anstelle eines Dateipfads bezeichnet wird.
Nextflow übernimmt die Authentifizierung und das Staging der Dateien an den richtigen Ort, das Herunterladen oder Hochladen und alle anderen Dateioperationen, die du erwarten würdest.

Die Hauptstärke dieses Systems besteht darin, dass es uns ermöglicht, zwischen Umgebungen zu wechseln, ohne die Pipeline-Logik zu ändern.
Du kannst zum Beispiel mit einem kleinen lokalen Testdatensatz entwickeln, bevor du zu einem vollständigen Testdatensatz wechselst, der sich in entferntem Speicher befindet, indem du einfach den URI änderst.

### 2.1. Verwende eine Datei aus dem Internet

Lass uns dies testen, indem wir den lokalen Pfad, den wir unserem Workflow bereitstellen, durch einen HTTPS-Pfad ersetzen, der auf eine Kopie derselben Daten zeigt, die in Github gespeichert ist.

!!! warning

    Dies funktioniert nur, wenn du eine aktive Internetverbindung hast.

Öffne `main.nf` erneut und ändere den Eingabepfad wie folgt:

=== "Nachher"

    ```groovy title="main.nf" linenums="2" hl_lines="2"
        // Using a remote file from the internet
        myFile = file('https://raw.github.com/nextflow-io/training/master/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz')

        // Print file attributes
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="2" hl_lines="2"
        // Create a Path object from a string path
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Print file attributes
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

Der eine Unterschied in der Konsolenausgabe ist, dass die Pfad-Objekt-Klasse jetzt `nextflow.file.http.XPath` ist, während sie für den lokalen Pfad `sun.nio.fs.UnixPath` war.
Du musst dir diese Klassen nicht merken; wir erwähnen dies nur, um zu demonstrieren, dass Nextflow die verschiedenen Speicherorte angemessen identifiziert und verarbeitet.

Hinter den Kulissen hat Nextflow die Datei in ein Staging-Verzeichnis innerhalb des Arbeitsverzeichnisses heruntergeladen.
Diese gestaged Datei kann dann als lokale Datei behandelt und in das entsprechende Prozessverzeichnis per Symlink eingebunden werden.

Du kannst überprüfen, dass dies geschehen ist, indem du den Inhalt des Arbeitsverzeichnisses an der Hash-Position des Prozesses ansiehst.

??? abstract "Arbeitsverzeichnisinhalte"

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

    Der Symlink zeigt auf eine gestaged Kopie der entfernten Datei, die Nextflow automatisch heruntergeladen hat.

Beachte, dass bei größeren Dateien der Download-Schritt im Vergleich zur Ausführung auf lokalen Dateien zusätzliche Zeit in Anspruch nimmt.
Allerdings prüft Nextflow, ob es bereits eine gestaged Kopie hat, um unnötige Downloads zu vermeiden.
Wenn du also erneut auf derselben Datei ausführst und die gestaged Datei nicht gelöscht hast, wird Nextflow die gestaged Kopie verwenden.

Dies zeigt, wie einfach es ist, mit Nextflow zwischen lokalen und entfernten Daten zu wechseln, was ein Hauptmerkmal von Nextflow ist.

!!! note

    Die eine wichtige Ausnahme von diesem Prinzip ist, dass du keine Glob-Muster oder Verzeichnispfade mit HTTPS verwenden kannst, da HTTPS keine mehreren Dateien auflisten kann, sodass du exakte Datei-URLs angeben musst.
    Andere Speicherprotokolle wie Blob Storage (`s3://`, `az://`, `gs://`) können jedoch sowohl Globs als auch Verzeichnispfade verwenden.

    So könntest du Glob-Muster mit Cloud-Speicher verwenden:

    ```groovy title="Cloud-Speicher-Beispiele (in dieser Umgebung nicht ausführbar)"
    // S3 mit Glob-Mustern - würde mehrere Dateien matchen
    ch_s3_files = channel.fromPath('s3://my-bucket/data/*.fastq.gz')

    // Azure Blob Storage mit Glob-Mustern
    ch_azure_files = channel.fromPath('az://container/data/patient*_R{1,2}.fastq.gz')

    // Google Cloud Storage mit Glob-Mustern
    ch_gcs_files = channel.fromPath('gs://bucket/data/sample_*.fastq.gz')
    ```

    Wir zeigen dir im nächsten Abschnitt, wie du in der Praxis mit Globs arbeitest.

### 2.2. Wechsle zurück zur lokalen Datei

Wir werden für den Rest dieser Side Quest wieder unsere lokalen Beispieldateien verwenden, also lass uns den Workflow-Input zurück auf die ursprüngliche Datei umstellen:

=== "Nachher"

    ```groovy title="main.nf" linenums="2" hl_lines="2"
        // Create a Path object from a string path
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Print file attributes
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="2" hl_lines="2"
        // Create a Path object from a string path
        myFile = file('https://raw.github.com/nextflow-io/training/master/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz')

        // Print file attributes
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    ```

### Zusammenfassung

- Auf entfernte Daten wird über einen URI zugegriffen (HTTP, FTP, S3, Azure, Google Cloud)
- Nextflow lädt die Daten automatisch herunter und stellt sie am richtigen Ort bereit, solange diese Pfade an Prozesse übergeben werden
- Schreibe keine Logik zum Herunterladen oder Hochladen entfernter Dateien!
- Lokale und entfernte Dateien erzeugen verschiedene Objekttypen, funktionieren aber identisch
- **Wichtig**: HTTP/HTTPS funktionieren nur mit einzelnen Dateien (keine Glob-Muster)
- Cloud-Speicher (S3, Azure, GCS) unterstützt sowohl einzelne Dateien als auch Glob-Muster
- Du kannst nahtlos zwischen lokalen und entfernten Datenquellen wechseln, ohne Code-Logik zu ändern (solange das Protokoll deine erforderlichen Operationen unterstützt)

---

## 3. Verwendung der `fromPath()` Channel Factory

Bisher haben wir mit jeweils einer einzelnen Datei gearbeitet, aber in Nextflow möchten wir normalerweise einen Input-Channel mit mehreren zu verarbeitenden Eingabedateien erstellen.

Eine naive Art, dies zu tun, wäre die `file()`-Methode mit [`channel.of()`](https://www.nextflow.io/docs/latest/reference/channel.html#of) zu kombinieren, etwa so:

```groovy title="Syntaxbeispiel"
ch_files = channel.of([file('data/patientA_rep1_normal_R1_001.fastq.gz')],
                      [file('data/patientA_rep1_normal_R1_001.fastq.gz')])
```

Das funktioniert, ist aber umständlich.

!!! tip "Wann `file()` vs. `channel.fromPath()` verwenden"

    - Verwende `file()`, wenn du ein einzelnes Path-Objekt für direkte Manipulation benötigst (Prüfen, ob eine Datei existiert, Lesen ihrer Attribute oder Übergabe an einen einzelnen Prozessaufruf)
    - Verwende `channel.fromPath()`, wenn du einen Channel benötigst, der mehrere Dateien enthalten kann, insbesondere mit Glob-Mustern, oder wenn Dateien durch mehrere Prozesse fließen werden

Hier kommt [`channel.fromPath()`](https://www.nextflow.io/docs/latest/reference/channel.html#frompath) ins Spiel: eine bequeme Channel Factory, die alle Funktionalität bündelt, die wir benötigen, um einen Channel aus einem oder mehreren statischen Datei-Strings sowie Glob-Mustern zu generieren.

### 3.1. Füge die Channel Factory hinzu

Lass uns unseren Workflow aktualisieren, um `channel.fromPath` zu verwenden.

=== "Nachher"

    ```groovy title="main.nf" linenums="7" hl_lines="1-3"
        // Load files with channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R1_001.fastq.gz')
        ch_files.view { myFile -> "Found file: $myFile" }

        // Print file attributes
        /* Comment these out for now, we'll come back to them!
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
        */

        // Count the lines in the file
        // COUNT_LINES(myFile)
    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="7" hl_lines="1-2"
        // Create a Path object from a string path
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Print file attributes
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"

        // Count the lines in the file
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

Wie du sehen kannst, wird der Dateipfad als `Path`-Objekt im Channel geladen.
Dies ist ähnlich wie bei `file()`, nur dass wir jetzt einen Channel haben, in den wir bei Bedarf mehr Dateien laden können.

Die Verwendung von `channel.fromPath()` ist eine bequeme Möglichkeit, einen neuen Channel zu erstellen, der mit einer Liste von Dateien gefüllt ist.

### 3.2. Zeige Attribute von Dateien im Channel an

In unserem ersten Versuch mit der Channel Factory haben wir den Code vereinfacht und nur den Dateinamen ausgegeben.

Lass uns zurück zur Ausgabe der vollständigen Dateiattribute gehen:

=== "Nachher"

    ```groovy title="main.nf" linenums="7" hl_lines="3-9 12"
        // Load files with channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R1_001.fastq.gz')
        ch_files.view { myFile ->
            println "File object class: ${myFile.class}"
            println "File name: ${myFile.name}"
            println "Simple name: ${myFile.simpleName}"
            println "Extension: ${myFile.extension}"
            println "Parent directory: ${myFile.parent}"
        }

        // Count the lines in the file
        COUNT_LINES(ch_files)
    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="7" hl_lines="3"
        // Load files with channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R1_001.fastq.gz')
        ch_files.view { myFile -> "Found file: $myFile" }

        // Count the lines in the file
        // COUNT_LINES(ch_files)
    ```

Wir aktivieren auch den `COUNT_LINES`-Prozessaufruf erneut, um zu überprüfen, dass die Dateiverarbeitung mit unserem Channel-basierten Ansatz noch korrekt funktioniert.

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

Und da hast du es, dieselben Ergebnisse wie zuvor, aber jetzt haben wir die Datei in einem Channel, sodass wir mehr hinzufügen können.

### 3.3. Verwendung eines Globs zum Matchen mehrerer Dateien

Es gibt mehrere Möglichkeiten, mehr Dateien in den Channel zu laden.
Hier zeigen wir dir, wie man Glob-Muster verwendet, eine bequeme Möglichkeit, Datei- und Verzeichnisnamen basierend auf Wildcard-Zeichen zu matchen und abzurufen.
Der Prozess des Matchens dieser Muster wird "Globbing" oder "Filename Expansion" genannt.

!!! note

    Wie bereits erwähnt, unterstützt Nextflow Globbing zur Verwaltung von Input- und Output-Dateien in den meisten Fällen, außer bei HTTPS-Dateipfaden, da HTTPS keine mehreren Dateien auflisten kann.

Angenommen, wir möchten beide Dateien in einem Paar von Dateien abrufen, die mit einem bestimmten Patienten, `patientA`, assoziiert sind:

```console
patientA_rep1_normal_R1_001.fastq.gz
patientA_rep1_normal_R2_001.fastq.gz
```

Da der einzige Unterschied zwischen den Dateinamen die Replikat-Nummer ist, _d.h._ die Zahl nach `R`, können wir das Wildcard-Zeichen `*` anstelle der Zahl wie folgt verwenden:

```console
patientA_rep1_normal_R*_001.fastq.gz
```

Das ist das Glob-Muster, das wir benötigen.

Jetzt müssen wir nur noch den Dateipfad in der Channel Factory aktualisieren, um dieses Glob-Muster wie folgt zu verwenden:

=== "Nachher"

    ```groovy title="main.nf" linenums="7"
      // Load files with channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="7"
      // Load files with channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R1_001.fastq.gz')
    ```

Nextflow erkennt automatisch, dass dies ein Glob-Muster ist und wird es entsprechend verarbeiten.

Führe den Workflow aus, um dies zu testen:

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

Wie du sehen kannst, haben wir jetzt zwei Path-Objekte in unserem Channel, was zeigt, dass Nextflow die Dateinamen-Expansion korrekt durchgeführt und beide Dateien wie erwartet geladen und verarbeitet hat.

Mit dieser Methode können wir beliebig viele oder wenige Dateien abrufen, indem wir einfach das Glob-Muster ändern. Wenn wir es großzügiger gestalten würden, zum Beispiel indem wir alle variablen Teile der Dateinamen durch `*` ersetzen (_z.B._ `data/patient*_rep*_*_R*_001.fastq.gz`), könnten wir alle Beispieldateien im `data`-Verzeichnis erfassen.

### Zusammenfassung

- `channel.fromPath()` erstellt einen Channel mit Dateien, die einem Muster entsprechen
- Jede Datei wird als separates Element im Channel ausgegeben
- Wir können ein Glob-Muster verwenden, um mehrere Dateien zu matchen
- Dateien werden automatisch in Path-Objekte mit vollständigen Attributen umgewandelt
- Die `.view()`-Methode ermöglicht die Inspektion der Channel-Inhalte

---

## 4. Extraktion grundlegender Metadaten aus Dateinamen

In den meisten wissenschaftlichen Bereichen ist es sehr üblich, Metadaten in den Namen der Dateien zu kodieren, die die Daten enthalten.
Zum Beispiel werden in der Bioinformatik Dateien mit Sequenzierungsdaten oft so benannt, dass Informationen über die Probe, Bedingung, Replikat und Read-Nummer kodiert sind.

Wenn die Dateinamen nach einer konsistenten Konvention konstruiert sind, kannst du diese Metadaten auf standardisierte Weise extrahieren und im Verlauf deiner Analyse verwenden.
Das ist natürlich ein großes 'Wenn', und du solltest sehr vorsichtig sein, wenn du dich auf die Dateinamenstruktur verlässt; aber die Realität ist, dass dieser Ansatz sehr weit verbreitet ist, also schauen wir uns an, wie es in Nextflow gemacht wird.

Im Fall unserer Beispieldaten wissen wir, dass die Dateinamen konsistent strukturierte Metadaten enthalten.
Zum Beispiel kodiert der Dateiname `patientA_rep1_normal_R2_001` Folgendes:

- Patienten-ID: `patientA`
- Replikat-ID: `rep1`
- Probentyp: `normal` (im Gegensatz zu `tumor`)
- Read-Set: `R1` (im Gegensatz zu `R2`)

Wir werden unseren Workflow in drei Schritten modifizieren, um diese Informationen abzurufen:

1. Den `simpleName` der Datei abrufen, der die Metadaten enthält
2. Die Metadaten mit einer Methode namens `tokenize()` trennen
3. Eine Map verwenden, um die Metadaten zu organisieren

!!! warning

    Du solltest niemals sensible Informationen in Dateinamen kodieren, wie Patientennamen oder andere identifizierende Merkmale, da dies die Privatsphäre von Patienten oder andere relevante Sicherheitsbeschränkungen gefährden kann.

### 4.1. Abrufen des `simpleName`

Der `simpleName` ist ein Dateiattribut, das dem Dateinamen ohne Pfad und Extension entspricht.

Nimm folgende Änderungen am Workflow vor:

=== "Nachher"

    ```groovy title="main.nf" linenums="7" hl_lines="3-6"
        // Load files with channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
        ch_files.map { myFile ->
            [ myFile.simpleName, myFile ]
        }
        .view()
    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="7" hl_lines="3-9"
        // Load files with channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
        ch_files.view { myFile ->
            println "File object class: ${myFile.class}"
            println "File name: ${myFile.name}"
            println "Simple name: ${myFile.simpleName}"
            println "Extension: ${myFile.extension}"
            println "Parent directory: ${myFile.parent}"
        }
    ```

Dies ruft den `simpleName` ab und verknüpft ihn mit dem vollständigen Dateiobjekt über eine `map()`-Operation.

Führe den Workflow aus, um zu testen, dass es funktioniert:

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

Jedes Element im Channel ist jetzt ein Tupel, das den `simpleName` und das ursprüngliche Dateiobjekt enthält.

### 4.2. Extrahiere die Metadaten aus dem `simplename`

Zu diesem Zeitpunkt sind die Metadaten, die wir wollen, im `simplename` eingebettet, aber wir können nicht direkt auf einzelne Elemente zugreifen.
Also müssen wir den `simplename` in seine Komponenten aufteilen.
Glücklicherweise sind diese Komponenten im ursprünglichen Dateinamen einfach durch Unterstriche getrennt, sodass wir eine gängige Nextflow-Methode namens `tokenize()` anwenden können, die perfekt für diese Aufgabe ist.

Nimm folgende Änderungen am Workflow vor:

=== "Nachher"

    ```groovy title="main.nf" linenums="7" hl_lines="4"
        // Load files with channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
        ch_files.map { myFile ->
            [ myFile.simpleName.tokenize('_'), myFile ]
        }
    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="7" hl_lines="4"
        // Load files with channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
        ch_files.map { myFile ->
            [ myFile.simpleName, myFile ]
        }
    ```

Die `tokenize()`-Methode teilt den `simpleName`-String überall dort, wo sie Unterstriche findet, und gibt eine Liste mit den Teilstrings zurück.

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

Jetzt enthält das Tupel für jedes Element in unserem Channel die Liste der Metadaten (_z.B._ `[patientA, rep1, normal, R1, 001]`) und das ursprüngliche Dateiobjekt.

Das ist großartig!
Wir haben unsere Patienteninformationen von einem einzelnen String in eine Liste von Strings aufgeteilt.
Wir können nun jeden Teil der Patienteninformationen separat verarbeiten.

### 4.3. Verwende eine Map zur Organisation der Metadaten

Unsere Metadaten sind im Moment nur eine flache Liste.
Sie ist einfach genug zu verwenden, aber schwer zu lesen.

```console
[patientA, rep1, normal, R1, 001]
```

Was ist das Element an Index 3? Kannst du es sagen, ohne auf die ursprüngliche Erklärung der Metadatenstruktur zurückzugreifen?

Dies ist eine großartige Gelegen
