# Verarbeitung von Dateieingaben

Wissenschaftliche Analyse-Workflows beinhalten oft die Verarbeitung einer großen Anzahl von Dateien.
Nextflow bietet leistungsstarke Werkzeuge für die effiziente Handhabung von Dateien und hilft dir dabei, deine Daten mit minimalem Code zu organisieren und zu verarbeiten.

### Lernziele

In dieser Side Quest werden wir untersuchen, wie Nextflow mit Dateien umgeht, von grundlegenden Dateioperationen bis hin zu fortgeschritteneren Techniken für die Arbeit mit Dateisammlungen.
Du wirst lernen, wie du Metadaten aus Dateinamen extrahieren kannst, was eine häufige Anforderung in wissenschaftlichen Analyse-Pipelines ist.

Am Ende dieser Side Quest wirst du in der Lage sein:

- Path-Objekte aus Dateipfad-Strings mit Nextflows `file()`-Methode zu erstellen
- Auf Dateiattribute wie Name, Erweiterung und übergeordnetes Verzeichnis zuzugreifen
- Sowohl lokale als auch entfernte Dateien transparent über URIs zu verwalten
- Channels zur automatischen Dateiverarbeitung mit `channel.fromPath()` und `channel.fromFilePairs()` zu verwenden
- Metadaten aus Dateinamen mit String-Manipulation zu extrahieren und zu strukturieren
- Zusammengehörige Dateien mit Pattern Matching und Glob-Ausdrücken zu gruppieren
- Dateioperationen in Nextflow-Prozesse mit korrekter Eingabebehandlung zu integrieren
- Prozessausgaben mit metadatengesteuerten Verzeichnisstrukturen zu organisieren

Diese Fähigkeiten werden dir helfen, Workflows zu erstellen, die mit verschiedenen Arten von Dateieingaben mit großer Flexibilität umgehen können.

### Voraussetzungen

Bevor du diese Side Quest in Angriff nimmst, solltest du:

- Das Tutorial [Hello Nextflow](../../hello_nextflow/) oder einen gleichwertigen Anfängerkurs abgeschlossen haben.
- Mit der Verwendung grundlegender Nextflow-Konzepte und -Mechanismen (Prozesse, Channels, Operatoren) vertraut sein

<!-- I removed the suggestion to do the metamaps SQ first because that works more naturally after -->

---

## 0. Erste Schritte

#### Öffne den Trainings-Codespace

Wenn du es noch nicht getan hast, stelle sicher, dass du die Trainingsumgebung öffnest, wie in der [Umgebungseinrichtung](../envsetup/index.md) beschrieben.

[![In GitHub Codespaces öffnen](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Wechsle in das Projektverzeichnis

Lass uns in das Verzeichnis wechseln, in dem sich die Dateien für dieses Tutorial befinden.

```bash
cd side-quests/working_with_files
```

Du kannst VSCode auf dieses Verzeichnis fokussieren:

```bash
code .
```

#### Überprüfe die Materialien

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

Dieses Verzeichnis enthält Paired-End-Sequenzierungsdaten von drei Patienten (A, B, C).

Für jeden Patienten haben wir Proben, die vom Typ `tumor` (typischerweise aus Tumorbiopsien stammend) oder `normal` (aus gesundem Gewebe oder Blut entnommen) sind.
Wenn du mit Krebsanalyse nicht vertraut bist, solltest du wissen, dass dies einem experimentellen Modell entspricht, das gepaarte Tumor-/Normalproben verwendet, um kontrastive Analysen durchzuführen.

Für Patient A haben wir speziell zwei Sätze von technischen Replikaten (Wiederholungen).

Die Sequenzierungsdaten-Dateien werden mit einer typischen `_R1_`- und `_R2_`-Konvention für sogenannte "Forward-Reads" und "Reverse-Reads" benannt.

_Keine Sorge, wenn du mit diesem experimentellen Design nicht vertraut bist, es ist nicht entscheidend für das Verständnis dieses Tutorials._

#### Überprüfe die Aufgabe

Deine Herausforderung besteht darin, einen Nextflow-Workflow zu schreiben, der:

1. **Lädt** Eingabedateien mit Nextflow-Methoden zur Dateibehandlung
2. **Extrahiert** Metadaten (Patienten-ID, Replikat, Probentyp) aus der Dateinamenstruktur
3. **Gruppiert** zusammengehörige Dateien (R1/R2) mit `channel.fromFilePairs()`
4. **Verarbeitet** die Dateien mit einem bereitgestellten Analysemodul
5. **Organisiert** Ausgaben in einer Verzeichnisstruktur basierend auf den extrahierten Metadaten

#### Bereitschaftsliste

Denkst du, du bist bereit einzusteigen?

- [ ] Ich verstehe das Ziel dieses Kurses und seine Voraussetzungen
- [ ] Mein Codespace läuft
- [ ] Ich habe mein Arbeitsverzeichnis entsprechend eingestellt
- [ ] Ich verstehe die Aufgabe

Wenn du alle Kästchen ankreuzen kannst, bist du startklar.

---

## 1. Grundlegende Dateioperationen

### 1.1. Den Typ eines Objekts mit `.class` identifizieren

Werfen wir einen Blick auf die Workflow-Datei `main.nf`:

```groovy title="main.nf" linenums="1"
#!/usr/bin/env nextflow

workflow {

    // Create a Path object from a string path
    myFile = 'data/patientA_rep1_normal_R1_001.fastq.gz'

    println "${myFile} is of class ${myFile.class}"
}
```

Dies ist ein Mini-Workflow (ohne jegliche Prozesse), der auf einen einzelnen Dateipfad in seinem Workflow verweist und ihn dann zusammen mit seiner Klasse auf der Konsole ausgibt.

??? info "Was ist `.class`?"

    In Nextflow sagt uns `.class`, mit welcher Art von Objekt wir arbeiten. Es ist, als würde man fragen: "Was für eine Art von Ding ist das?", um herauszufinden, ob es sich um einen String, eine Zahl, eine Datei oder etwas anderes handelt.
    Dies wird uns helfen, den Unterschied zwischen einem einfachen String und einem Path-Objekt in den nächsten Abschnitten zu veranschaulichen.

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

Wie du sehen kannst, hat Nextflow den String-Pfad genau so ausgegeben, wie wir ihn geschrieben haben.

Dies ist nur eine Textausgabe; Nextflow hat noch nichts Besonderes damit gemacht.
Wir haben auch bestätigt, dass es sich aus Nextflows Sicht nur um einen String (der Klasse `java.lang.String`) handelt.
Das ergibt Sinn, da wir Nextflow noch nicht mitgeteilt haben, dass es einer Datei entspricht.

### 1.2. Ein Path-Objekt mit file() erstellen

Wir können Nextflow mitteilen, wie mit Dateien umzugehen ist, indem wir [Path-Objekte](https://www.nextflow.io/docs/latest/reference/stdlib-types.html#path) aus Pfad-Strings erstellen.

In unserem Workflow können wir den String-Pfad `data/patientA_rep1_normal_R1_001.fastq.gz` in ein Path-Objekt umwandeln, indem wir die `file()`-Methode verwenden, die Zugriff auf Dateieigenschaften und -operationen bietet.

Bearbeite die `main.nf`, um den String mit `file()` zu umschließen:

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

Dieses Mal siehst du den vollständigen absoluten Pfad anstelle des relativen Pfades, den wir als Eingabe angegeben haben.

Nextflow hat unseren String in ein Path-Objekt umgewandelt und ihn zum tatsächlichen Dateispeicherort im System aufgelöst.
Der Dateipfad wird jetzt absolut sein, wie in `/workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz`.

Beachte auch, dass die Path-Objektklasse `sun.nio.fs.UnixPath` ist: Dies ist Nextflows Art, lokale Dateien darzustellen.
Wie wir später sehen werden, haben entfernte Dateien andere Klassennamen (wie `nextflow.file.http.XPath` für HTTP-Dateien), aber sie funktionieren alle genau gleich und können in deinen Workflows identisch verwendet werden.

!!! tip

    **Der wichtigste Unterschied:**

    - **Pfad-String**: Nur Text, den Nextflow als Zeichen behandelt
    - **Path-Objekt**: Eine intelligente Dateireferenz, mit der Nextflow arbeiten kann

    Stell dir das so vor: Ein Pfad-String ist wie eine auf Papier geschriebene Adresse, während ein Path-Objekt wie eine in ein GPS-Gerät geladene Adresse ist, das weiß, wie man dorthin navigiert und dir Details über die Reise mitteilen kann.

### 1.3. Auf Dateiattribute zugreifen

Warum ist das hilfreich? Nun, da Nextflow jetzt versteht, dass `myFile` ein Path-Objekt und nicht nur ein String ist, können wir auf die verschiedenen Attribute des Path-Objekts zugreifen.

Aktualisieren wir unseren Workflow, um die eingebauten Dateiattribute auszugeben:

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

Du siehst die verschiedenen Dateiattribute, die auf der Konsole ausgegeben werden.

### 1.4. Die Datei an einen Prozess übergeben

Der Unterschied zwischen Strings und Path-Objekten wird entscheidend, wenn du beginnst, tatsächliche Workflows mit Prozessen zu erstellen.
Bisher haben wir überprüft, dass Nextflow unsere Eingabedatei jetzt als Datei behandelt, aber lassen wir sehen, ob wir tatsächlich etwas mit dieser Datei in einem Prozess ausführen können.

#### 1.4.1. Den Prozess importieren und den Code untersuchen

Wir stellen dir ein bereits geschriebenes Prozessmodul namens `COUNT_LINES` zur Verfügung, das eine Dateieingabe nimmt und zählt, wie viele Zeilen sie enthält.

Um den Prozess im Workflow zu verwenden, musst du nur eine Include-Anweisung vor dem Workflow-Block hinzufügen:

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

Wie du sehen kannst, handelt es sich um ein relativ unkompliziertes kleines Skript, das die Datei entpackt und zählt, wie viele Zeilen sie enthält.

??? info "Was bewirkt `debug true`?"

    Die `debug true`-Direktive in der Prozessdefinition bewirkt, dass Nextflow die Ausgabe deines Skripts (wie die Zeilenanzahl "40") direkt im Ausführungsprotokoll ausgibt.
    Ohne dies würdest du nur den Prozessausführungsstatus sehen, aber nicht die tatsächliche Ausgabe deines Skripts.

    Weitere Informationen zum Debuggen von Nextflow-Prozessen findest du in der Side Quest [Debugging Nextflow Workflows](debugging.md).

#### 1.4.2. Einen Aufruf zu `COUNT_LINES` hinzufügen

Jetzt, da der Prozess für den Workflow verfügbar ist, können wir einen Aufruf zum `COUNT_LINES`-Prozess hinzufügen, um ihn auf der Eingabedatei auszuführen.

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

Und führe jetzt den Workflow aus:

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

Dies zeigt, dass wir in der Lage sind, die Datei im Prozess angemessen zu verarbeiten.

Konkret hat Nextflow die folgenden Operationen erfolgreich durchgeführt:

- Die Datei ins Arbeitsverzeichnis gestellt
- Die .gz-Datei dekomprimiert
- Die Zeilen gezählt (40 Zeilen in diesem Fall)
- Ohne Fehler abgeschlossen

Der Schlüssel zu diesem reibungslosen Betrieb ist, dass wir Nextflow explizit mitteilen, dass unsere Eingabe eine Datei ist und als solche behandelt werden sollte.

### 1.5. Grundlegende Dateieingabefehler beheben

Dies verwirrt Neueinsteiger in Nextflow oft, also nehmen wir uns ein paar Minuten Zeit, um zu sehen, was passiert, wenn du es falsch machst.

Es gibt zwei Hauptstellen, an denen du die Dateiverarbeitung falsch machen kannst: auf der Ebene des Workflows und auf der Ebene des Prozesses.

#### 1.5.1. Workflow-Ebene-Fehler

Sehen wir, was passiert, wenn wir zur Behandlung der Datei als String zurückkehren, wenn wir die Eingabe im Workflow-Block spezifizieren.

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

Und führe jetzt den Workflow aus:

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

Der wichtige Teil ist:

```console
Not a valid path value: 'data/patientA_rep1_normal_R1_001.fastq.gz'
```

Wenn du eine `path`-Eingabe spezifizierst, überprüft Nextflow, dass du tatsächliche Dateireferenzen übergibst, nicht nur Strings.
Dieser Fehler teilt dir mit, dass `'data/patientA_rep1_normal_R1_001.fastq.gz'` kein gültiger Pfadwert ist, weil es ein String ist, kein Path-Objekt.

Nextflow hat das Problem sofort erkannt und angehalten, bevor es überhaupt den Prozess gestartet hat.

#### 1.5.2. Prozess-Ebene-Fehler

Die andere Stelle, an der wir vergessen könnten anzugeben, dass Nextflow die Eingabe als Datei behandeln soll, ist in der Prozessdefinition.

!!! warning "Behalte den Workflow-Fehler aus 1.5.1 bei"

    Damit dieser Test korrekt funktioniert, behalte den Workflow in seinem fehlerhaften Zustand bei (verwendet einen einfachen String anstelle von `file()`).
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

Und führe den Workflow jetzt erneut aus:

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

Dies besagt, dass das System die Datei nicht finden konnte; wenn du jedoch den Pfad überprüfst, gibt es eine Datei mit diesem Namen an diesem Speicherort.

Als wir dies ausführten, hat Nextflow den String-Wert an das Skript weitergegeben, aber es hat die eigentliche Datei nicht im Arbeitsverzeichnis _bereitgestellt_.
Also versuchte der Prozess, den relativen String `data/patientA_rep1_normal_R1_001.fastq.gz` zu verwenden, aber diese Datei existiert nicht innerhalb des Prozessarbeitsverzeichnisses.

Zusammen zeigen diese beiden Beispiele, wie wichtig es ist, Nextflow mitzuteilen, ob eine Eingabe als Datei behandelt werden soll.

!!! note

    Stelle sicher, dass du beide beabsichtigten Fehler behebst, bevor du zum nächsten Abschnitt übergehst.

### Zusammenfassung

- Pfad-Strings vs. Path-Objekte: Strings sind nur Text, Path-Objekte sind intelligente Dateireferenzen
- Die `file()`-Methode konvertiert einen String-Pfad in ein Path-Objekt, mit dem Nextflow arbeiten kann
- Du kannst auf Dateieigenschaften wie `name`, `simpleName`, `extension` und `parent` [über Dateiattribute](https://www.nextflow.io/docs/latest/working-with-files.html#getting-file-attributes) zugreifen
- Die Verwendung von Path-Objekten anstelle von Strings ermöglicht es Nextflow, Dateien in deinem Workflow richtig zu verwalten
- Ergebnisse der Prozesseingabe: Die richtige Dateibehandlung erfordert Path-Objekte, keine Strings, um sicherzustellen, dass Dateien korrekt bereitgestellt und für die Verwendung durch Prozesse zugänglich sind.

---

## 2. Verwendung von entfernten Dateien

Eine der Hauptfunktionen von Nextflow ist die Fähigkeit, nahtlos zwischen lokalen Dateien (auf demselben Rechner) und entfernten Dateien, die über das Internet zugänglich sind, zu wechseln.

Wenn du es richtig machst, solltest du niemals die Logik deines Workflows ändern müssen, um Dateien aus verschiedenen Quellen zu berücksichtigen.
Alles, was du tun musst, um eine entfernte Datei zu verwenden, ist, das entsprechende Präfix im Dateipfad anzugeben, wenn du ihn dem Workflow zur Verfügung stellst.

Zum Beispiel hat `/path/to/data` kein Präfix, was darauf hinweist, dass es sich um einen "normalen" lokalen Dateipfad handelt, während `s3://path/to/data` das `s3://`-Präfix enthält, was darauf hinweist, dass es sich in Amazons S3-Objektspeicher befindet.

Viele verschiedene Protokolle werden unterstützt:

- HTTP(S)/FTP (http://, https://, ftp://)
- Amazon S3 (s3://)
- Azure Blob Storage (az://)
- Google Cloud Storage (gs://)

Um eines davon zu verwenden, gib einfach das relevante Präfix im String an, der dann technisch gesehen ein Uniform Resource Identifier (URI) statt ein Dateipfad genannt wird.
Nextflow wird sich um Authentifizierung und Bereitstellung der Dateien am richtigen Ort kümmern, herunterladen oder hochladen und alle anderen Dateioperationen, die du erwarten würdest.

Die Hauptstärke dieses Systems ist, dass es uns ermöglicht, zwischen Umgebungen zu wechseln, ohne jegliche Pipeline-Logik zu ändern.
Du kannst zum Beispiel mit einem kleinen, lokalen Testsatz entwickeln, bevor du zu einem vollständigen Testsatz wechselst, der sich in einem entfernten Speicher befindet, einfach durch Ändern des URI.

### 2.1. Eine Datei aus dem Internet verwenden

Testen wir dies, indem wir den lokalen Pfad, den wir unserem Workflow zur Verfügung stellen, durch einen HTTPS-Pfad ersetzen, der auf eine Kopie derselben Daten zeigt, die in Github gespeichert sind.

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

Der einzige Unterschied in der Konsolenausgabe ist, dass die Pfadobjektklasse jetzt `nextflow.file.http.XPath` ist, während für den lokalen Pfad die Klasse `sun.nio.fs.UnixPath` war.
Du musst dir diese Klassen nicht merken; wir erwähnen dies nur, um zu demonstrieren, dass Nextflow die verschiedenen Orte identifiziert und angemessen behandelt.

Im Hintergrund hat Nextflow die Datei in ein Staging-Verzeichnis heruntergeladen, das sich im Arbeitsverzeichnis befindet.
Diese bereitgestellte Datei kann dann als lokale Datei behandelt und per Symlink in das relevante Prozessverzeichnis verknüpft werden.

Du kannst überprüfen, dass das hier passiert ist, indem du den Inhalt des Arbeitsverzeichnisses am Hash-Wert des Prozesses betrachtest.

??? abstract "Arbeitsverzeichnisinhalt"

    Wenn der Prozesshash `8a/2ab7ca` war, könntest du das Arbeitsverzeichnis untersuchen:

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

Beachte, dass bei größeren Dateien der Download-Schritt etwas mehr Zeit in Anspruch nehmen wird im Vergleich zur Ausführung mit lokalen Dateien.
Nextflow prüft jedoch, ob es bereits eine bereitgestellte Kopie hat, um unnötige Downloads zu vermeiden.
Wenn du also erneut dieselbe Datei ausführst und die bereitgestellte Datei nicht gelöscht hast, wird Nextflow die bereitgestellte Kopie verwenden.

Dies zeigt, wie einfach es ist, mit Nextflow zwischen lokalen und entfernten Daten zu wechseln, was eine Schlüsselfunktion von Nextflow ist.

!!! note

    Die eine wichtige Ausnahme zu diesem Prinzip ist, dass du keine Glob-Muster oder Verzeichnispfade mit HTTPS verwenden kannst, weil HTTPS keine mehreren Dateien auflisten kann, sodass du exakte Datei-URLs angeben musst.
    Andere Speicherprotokolle wie Blob-Speicher (`s3://`, `az://`, `gs://`) können jedoch sowohl Globs als auch Verzeichnispfade verwenden.

    Hier ist, wie du Glob-Muster mit Cloud-Speicher verwenden könntest:

    ```groovy title="Cloud-Speicher-Beispiele (nicht in dieser Umgebung ausführbar)"
    // S3 mit Glob-Mustern - würde mehrere Dateien abgleichen
    ch_s3_files = channel.fromPath('s3://my-bucket/data/*.fastq.gz')

    // Azure Blob Storage mit Glob-Mustern
    ch_azure_files = channel.fromPath('az://container/data/patient*_R{1,2}.fastq.gz')

    // Google Cloud Storage mit Glob-Mustern
    ch_gcs_files = channel.fromPath('gs://bucket/data/sample_*.fastq.gz')
    ```

    Wir werden dir im nächsten Abschnitt zeigen, wie du in der Praxis mit Globs arbeiten kannst.

### 2.2. Zurück zur lokalen Datei wechseln

Wir werden für den Rest dieser Side Quest zu unseren lokalen Beispieldateien zurückkehren, also wechseln wir die Workflow-Eingabe zurück zur ursprünglichen Datei:

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
- Schreibe keine Logik zum Herunterladen oder Hochladen von entfernten Dateien!
- Lokale und entfernte Dateien erzeugen unterschiedliche Objekttypen, funktionieren aber identisch
- **Wichtig**: HTTP/HTTPS funktioniert nur mit einzelnen Dateien (keine Glob-Muster)
- Cloud-Speicher (S3, Azure, GCS) unterstützt sowohl einzelne Dateien als auch Glob-Muster
- Du kannst nahtlos zwischen lokalen und entfernten Datenquellen wechseln, ohne den Code-Logik zu ändern (solange das Protokoll deine erforderlichen Operationen unterstützt)

---

## 3. Verwendung der `fromPath()` Channel Factory

Bisher haben wir mit jeweils einer Datei gearbeitet, aber in Nextflow werden wir typischerweise einen Input-Channel mit mehreren Eingabedateien zur Verarbeitung erstellen wollen.

Eine naive Art, dies zu tun, wäre, die `file()`-Methode mit [`channel.of()`](https://www.nextflow.io/docs/latest/reference/channel.html#of) zu kombinieren, wie hier:

```groovy title="Syntax-Beispiel"
ch_files = channel.of([file('data/patientA_rep1_normal_R1_001.fastq.gz')],
                      [file('data/patientA_rep1_normal_R1_001.fastq.gz')])
```

Das funktioniert, ist aber umständlich.

!!! tip "Wann `file()` vs `channel.fromPath()` verwenden"

    - Verwende `file()`, wenn du ein einzelnes Path-Objekt für direkte Manipulation benötigst (Überprüfen, ob eine Datei existiert, Lesen ihrer Attribute oder Übergeben an eine einzelne Prozessausführung)
    - Verwende `channel.fromPath()`, wenn du einen Channel benötigst, der mehrere Dateien aufnehmen kann, insbesondere mit Glob-Mustern, oder wenn Dateien durch mehrere Prozesse fließen werden

Hier kommt [`channel.fromPath()`](https://www.nextflow.io/docs/latest/reference/channel.html#frompath) ins Spiel: eine praktische Channel Factory, die alle Funktionen bündelt, die wir brauchen, um einen Channel aus einem oder mehreren statischen Datei-Strings sowie Glob-Mustern zu erzeugen.

### 3.1. Die Channel Factory hinzufügen

Aktualisieren wir unseren Workflow, um `channel.fromPath` zu verwenden.

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

Wir haben auch den Code, der die Attribute ausgibt, vorerst auskommentiert und eine `.view`-Anweisung hinzugefügt, um nur den Dateinamen auszugeben.

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

Wie du sehen kannst, wird der Dateipfad als Objekt vom Typ `Path` in den Channel geladen.
Dies ist ähnlich zu dem, was `file()` getan hätte, außer dass wir jetzt einen Channel haben, in den wir bei Bedarf mehr Dateien laden können.

Die Verwendung von `channel.fromPath()` ist eine bequeme Methode, um einen neuen Channel zu erstellen, der mit einer Liste von Dateien gefüllt ist.

### 3.2. Attribute von Dateien im Channel anzeigen

In unserem ersten Durchlauf mit der Channel Factory haben wir den Code vereinfacht und nur den Dateinamen ausgegeben.

Kehren wir zur Ausgabe der vollständigen Dateiattribute zurück:

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

Wir aktivieren auch wieder den Aufruf des `COUNT_LINES`-Prozesses, um zu überprüfen, dass die Dateiverarbeitung mit unserem Channel-basierten Ansatz immer noch korrekt funktioniert.

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

Und da haben wir es, die gleichen Ergebnisse wie zuvor, aber jetzt haben wir die Datei in einem Channel, sodass wir mehr hinzufügen können.

### 3.3. Einen Glob verwenden, um mehrere Dateien abzugleichen

Es gibt mehrere Möglichkeiten, wie wir mehr Dateien in den Channel laden könnten.
Hier werden wir dir zeigen, wie du Glob-Muster verwenden kannst, die eine praktische Methode sind, um Datei- und Verzeichnisnamen basierend auf Platzhalterzeichen abzugleichen und abzurufen.
Der Prozess des Abgleichens dieser Muster wird "Globbing" oder "Dateinamenserweiterung" genannt.

!!! note

    Wie bereits erwähnt, unterstützt Nextflow Globbing zur Verwaltung von Ein- und Ausgabedateien in der Mehrzahl der Fälle, außer bei HTTPS-Dateipfaden, da HTTPS keine mehreren Dateien auflisten kann.

Angenommen, wir möchten beide Dateien in einem Paar von Dateien abrufen, die einem bestimmten Patienten, `patientA`, zugeordnet sind:

```console
patientA_rep1_normal_R1_001.fastq.gz
patientA_rep1_normal_R2_001.fastq.gz
```

Da der einzige Unterschied zwischen den Dateinamen die Replikatnummer ist, _d.h._ die Zahl nach `R`, können wir das Platzhalterzeichen `*` verwenden, um für die Zahl einzuspringen:

```console
patientA_rep1_normal_R*_001.fastq.gz
```

Das ist das Glob-Muster, das wir brauchen.

Jetzt müssen wir nur noch den Dateipfad in der Channel Factory aktualisieren, um dieses Glob-Muster zu verwenden:

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

Nextflow wird automatisch erkennen, dass dies ein Glob-Muster ist, und es entsprechend behandeln.

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

Wie du sehen kannst, haben wir jetzt zwei Path-Objekte in unserem Channel, was zeigt, dass Nextflow die Dateinamenserweiterung korrekt durchgeführt hat und beide Dateien wie erwartet geladen und verarbeitet hat.

Mit dieser Methode können wir so viele oder so wenige Dateien abrufen, wie wir möchten, indem wir einfach das Glob-Muster ändern. Wenn wir es großzügiger gestalten würden, zum Beispiel indem wir alle variablen Teile der Dateinamen durch `*` ersetzen (_z.B._ `data/patient*_rep*_*_R*_001.fastq.gz`), könnten wir alle Beispieldateien im `data`-Verzeichnis erfassen.

### Zusammenfassung

- `channel.fromPath()` erstellt einen Channel mit Dateien, die einem Muster entsprechen
- Jede Datei wird als separates Element im Channel ausgegeben
- Wir können ein Glob-Muster verwenden, um mehrere Dateien abzugleichen
- Dateien werden automatisch in Path-Objekte mit vollständigen Attributen umgewandelt
- Die `.view()`-Methode ermöglicht die Inspektion von Channel-Inhalten

---

## 4. Extraktion grundlegender Metadaten aus Dateinamen

In den meisten wissenschaftlichen Bereichen ist es sehr üblich, dass Metadaten in den Namen der Dateien kodiert sind, die die Daten enthalten.
In der Bioinformatik werden beispielsweise Dateien mit Sequenzierungsdaten oft so benannt, dass Informationen über die Probe, den Zustand, das Replikat und die Read-Nummer kodiert sind.

Wenn die Dateinamen nach einer konsistenten Konvention erstellt werden, kannst du diese Metadaten standardisiert extrahieren und im Laufe deiner Analyse verwenden.
Das ist natürlich ein großes "Wenn", und du solltest sehr vorsichtig sein, wenn du dich auf die Dateinamenstruktur verlässt; aber die Realität ist, dass dieser Ansatz sehr weit verbreitet ist, also schauen wir uns an, wie er in Nextflow gemacht wird.

Im Fall unserer Beispieldaten wissen wir, dass die Dateinamen konsistent strukturierte Metadaten enthalten.
Zum Beispiel kodiert der Dateiname `patientA_rep1_normal_R2_001` Folgendes:

- Patienten-ID: `patientA`
- Replikat-ID: `rep1`
- Probentyp: `normal` (im Gegensatz zu `tumor`)
- Read-Set: `R1` (im Gegensatz zu `R2`)

Wir werden unseren Workflow modifizieren, um diese Informationen in drei Schritten abzurufen:

1. Den `simpleName` der Datei abrufen, der die Metadaten enthält
2. Die Metadaten mit einer Methode namens `tokenize()` trennen
3. Eine Map verwenden, um die Metadaten zu organisieren

!!! warning

    Du solltest niemals sensible Informationen in Dateinamen kodieren, wie Patientennamen oder andere identifizierende Merkmale, da dies die Patientenprivatsphäre oder andere relevante Sicherheitsbeschränkungen gefährden kann.

### 4.1. Den `simpleName` abrufen

Der `simpleName` ist ein Dateiattribut, das dem Dateinamen entspricht, der von seinem Pfad und seiner Erweiterung befreit wurde.

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

Dies ruft den `simpleName` ab und verknüpft ihn mit dem vollständigen Dateiobjekt mittels einer `map()`-Operation.

Führe den Workflow aus, um zu testen, dass er funktioniert:

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

### 4.2. Die Metadaten aus dem `simplename` extrahieren

An diesem Punkt sind die Metadaten, die wir wollen, im `simplename` eingebettet, aber wir können nicht direkt auf einzelne Elemente zugreifen.
Wir müssen also den `simplename` in seine Komponenten aufteilen.
Glücklicherweise sind diese Komponenten im ursprünglichen Dateinamen einfach durch Unterstriche getrennt, sodass wir eine häufige Nextflow-Methode namens `tokenize()` anwenden können, die perfekt für diese Aufgabe ist.

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

Die `tokenize()`-Methode wird den `simpleName`-String überall dort teilen, wo sie Unterstriche findet, und eine Liste mit den Teilstrings zurückgeben.

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
Wir haben unsere Patienteninformationen von einem einzelnen String in eine Liste von Strings zerlegt.
Wir können jetzt jeden Teil der Patienteninformationen separat behandeln.

### 4.3. Eine Map zur Organisation der Metadaten verwenden

Unsere Metadaten sind im Moment nur eine flache Liste.
Sie ist einfach zu verwenden, aber schwer zu lesen.

```console
[patientA, rep1, normal, R1, 001]
```

Was ist das Element bei Index 3? Kannst du es sagen, ohne auf die ursprüngliche Erklärung der Metadatenstruktur zurückzugreifen?

Das ist eine gute Gelegenheit, eine Map (Schlüssel-Wert-Speicher) zu verwenden, bei der jedes Element einen Satz von Schlüsseln und ihren zugehörigen Werten hat, sodass du einfach auf jeden Schlüssel verweisen kannst, um den entsprechenden Wert zu erhalten.

In unserem Beispiel bedeutet das, von dieser Organisation:

```groovy
data = [patientA, 1, normal, R1]

println data[3]
```

zu dieser zu wechseln:

```groovy
data = [id: patientA, replicate: 1, type: normal, readNum: 1]

println data.readNum
```

In Nextflow wird das als [Map](https://nextflow.io/docs/latest/script.html#maps) bezeichnet.

Lass uns unsere flache Liste jetzt in eine Map umwandeln.
Nimm folgende Änderungen am Workflow vor:

=== "Nachher"

    ```groovy title="main.nf" linenums="7" hl_lines="4-13"
        // Load files with channel.fromPath
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
        // Load files with channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
        ch_files.map { myFile ->
            [ myFile.simpleName.tokenize('_'), myFile ]
        }
    ```

Die wesentlichen Änderungen sind:

- **Destrukturierende Zuweisung**: `def (patient, replicate, type, readNum) = ...` extrahiert die tokenisierten Werte in benannte Variablen in einer Zeile
- **Map-Literal-Syntax**: `[id: patient, replicate: ...]` erstellt eine Map, in der jeder Schlüssel (wie `id`) mit einem Wert (wie `patient`) verknüpft ist
- **Verschachtelte Struktur**: Die äußere Liste `[..., myFile]` paart die Metadaten-Map mit dem ursprünglichen Dateiobjekt

Wir haben außerdem einige der Metadaten-Strings mit einer String-Ersetzungsmethode namens `replace()` vereinfacht, um einige unnötige Zeichen zu entfernen (_z.B._ `replicate.replace('rep', '')`, um nur die Nummer der Replikat-IDs zu behalten).

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

Jetzt sind die Metadaten übersichtlich beschriftet (_z.B._ `[id:patientA, replicate:1, type:normal, readNum:2]`), sodass es viel einfacher ist zu erkennen, was was ist.

Es wird auch viel einfacher sein, Metadatenelemente tatsächlich im Workflow zu nutzen, und macht unseren Code lesbarer und wartbarer.

### Zusammenfassung

- Wir können Dateinamen in Nextflow mit der vollen Leistung einer Programmiersprache verarbeiten
- Wir können die Dateinamen als Strings behandeln, um relevante Informationen zu extrahieren
- Die Verwendung von Methoden wie `tokenize()` und `replace()` ermöglicht es uns, Strings im Dateinamen zu manipulieren
- Die `.map()`-Operation transformiert Channel-Elemente unter Beibehaltung der Struktur
- Strukturierte Metadaten (Maps) machen Code lesbarer und wartbarer als Positionslisten

Als Nächstes schauen wir uns an, wie man mit gepaarten Datendateien umgeht.

---

## 5. Verarbeitung gepaarter Datendateien

Viele experimentelle Designs erzeugen gepaarte Datendateien, die von einer explizit gepaarten Handhabung profitieren.
Zum Beispiel werden in der Bioinformatik Sequenzierungsdaten oft in Form von gepaarten Reads erzeugt, d.h. Sequenz-Strings, die vom selben DNA-Fragment stammen (oft als 'Forward' und 'Reverse' bezeichnet, weil sie von entgegengesetzten Enden gelesen werden).

Das ist bei unseren Beispieldaten der Fall, wo R1 und R2 auf die beiden Read-Sets verweisen.

```console
data/patientA_rep1_normal_R1_001.fastq.gz
data/patientA_rep1_normal_R2_001.fastq.gz
```

Nextflow bietet eine spezialisierte Channel Factory für die Arbeit mit solchen gepaarten Dateien namens `channel.fromFilePairs()`, die Dateien automatisch basierend auf einem gemeinsamen Benennungsmuster gruppiert. Das ermöglicht es dir, die gepaarten Dateien mit weniger Aufwand enger zu verknüpfen.

Wir werden unseren Workflow modifizieren, um dies zu nutzen.
Das geschieht in zwei Schritten:

1. Die Channel Factory auf `channel.fromFilePairs()` umstellen
2. Die Metadaten extrahieren und mappen

### 5.1. Die Channel Factory auf `channel.fromFilePairs()` umstellen

Um `channel.fromFilePairs` zu verwenden, müssen wir das Muster angeben, das Nextflow verwenden soll, um die beiden Mitglieder eines Paares zu identifizieren.

Zurück zu unseren Beispieldaten können wir das Benennungsmuster wie folgt formalisieren:

```console
data/patientA_rep1_normal_R{1,2}_001.fastq.gz
```

Dies ähnelt dem Glob-Muster, das wir früher verwendet haben, außer dass es speziell die Teilstrings (entweder `1` oder `2` direkt nach dem R) auflistet, die die beiden Mitglieder des Paares identifizieren.

Aktualisieren wir den Workflow `main.nf` entsprechend:

=== "Nachher"

    ```groovy title="main.nf" linenums="7" hl_lines="1-2"
        // Load files with channel.fromFilePairs
        ch_files = channel.fromFilePairs('data/patientA_rep1_normal_R{1,2}_001.fastq.gz')
        /* Comment out the mapping for now, we'll come back to it!
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
        // Load files with channel.fromFilePairs
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

Wir haben die Channel Factory umgestellt und das Datei-Matching-Muster angepasst, und dabei die Map-Operation auskommentiert.
Wir werden sie später mit einigen Modifikationen wieder hinzufügen.

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

Oh je, dieses Mal ist die Ausführung fehlgeschlagen!

Der relevante Teil der Fehlermeldung ist hier:

```console
Not a valid path value: 'patientA_rep1_normal_R'
```

Das liegt daran, dass wir die Channel Factory geändert haben.
Bisher enthielt der ursprüngliche Input-Channel nur die Dateipfade.
Die gesamte Metadaten-Manipulation, die wir durchgeführt haben, hat die Channel-Inhalte nicht wirklich beeinflusst.

Jetzt, da wir die `.fromFilePairs`-Channel-Factory verwenden, sind die Inhalte des resultierenden Channels anders.
Wir sehen nur ein Channel-Element, bestehend aus einem Tupel mit zwei Einträgen: dem Teil des `simpleName`, der von beiden Dateien geteilt wird und als Identifikator dient, und einem Tupel mit den beiden Dateiobjekten im Format `id, [ file1, file2 ]`.

Das ist großartig, denn Nextflow hat die harte Arbeit übernommen, den Patientennamen zu extrahieren, indem es das gemeinsame Präfix untersucht und es als Patienten-Identifikator verwendet.

Allerdings bricht es unseren aktuellen Workflow.
Wenn wir `COUNT_LINES` noch auf die gleiche Weise ausführen wollten, ohne den Prozess zu ändern, müssten wir eine Mapping-Operation anwenden, um die Dateipfade zu extrahieren.
Aber das werden wir nicht tun, denn unser ultimatives Ziel ist es, einen anderen Prozess zu verwenden, `ANALYZE_READS`, der Dateipaare angemessen verarbeitet.

Lass uns also einfach den Aufruf von `COUNT_LINES` auskommentieren (oder löschen) und weitermachen.

=== "Nachher"

    ```groovy title="main.nf" linenums="26" hl_lines="2"
        // Count the lines in the file
        // COUNT_LINES(ch_files)
    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="26" hl_lines="2"
        // Count the lines in the file
        COUNT_LINES(ch_files)
    ```

Du kannst auch die `COUNT_LINES`-Include-Anweisung auskommentieren oder löschen, aber das wird keine funktionalen Auswirkungen haben.

Lass uns den Workflow jetzt erneut ausführen:

```bash
nextflow run main.nf
```

??? success "Befehlsausgabe"

    ```console hl_lines="5"
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [fabulous_davinci] DSL2 - revision: 22b53268dc

    [patientA_rep1_normal_R, [/workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz, /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R2_001.fastq.gz]]
    ```

Juhu, dieses Mal läuft der Workflow erfolgreich!

Allerdings müssen wir noch die restlichen Metadaten aus dem `id`-Feld extrahieren.

### 5.2. Metadaten aus Dateipaaren extrahieren und organisieren

Unsere vorherige `map`-Operation funktioniert nicht, weil sie nicht zur Datenstruktur passt, aber wir können sie anpassen.

Wir haben bereits Zugriff auf den tatsächlichen Patienten-Identifikator im String, den `fromFilePairs()` als Identifikator verwendet hat, sodass wir ihn zur Extraktion der Metadaten nutzen können, ohne den `simpleName` vom Path-Objekt zu holen, wie wir es zuvor getan haben.

Entferne die Kommentare bei der Map-Operation im Workflow und nimm folgende Änderungen vor:

=== "Nachher"

    ```groovy title="main.nf" linenums="7" hl_lines="3-4 9 11 13"
        // Load files with channel.fromFilePairs
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
        // Load files with channel.fromFilePairs
        ch_files = channel.fromFilePairs('data/patientA_rep1_normal_R{1,2}_001.fastq.gz')
        /* Comment out the mapping for now, we'll come back to it!
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

Dieses Mal beginnt die Map mit `id, files` statt nur `myFile`, und `tokenize()` wird auf `id` statt auf `myFile.simpleName` angewendet.

Beachte auch, dass wir `readNum` aus der `tokenize()`-Zeile entfernt haben; alle Teilstrings, die wir nicht explizit benennen (von links ausgehend), werden stillschweigend verworfen.
Wir können dies tun, weil die gepaarten Dateien jetzt eng verknüpft sind, sodass wir `readNum` nicht mehr in der Metadaten-Map benötigen.

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

Und da haben wir es: Wir haben die Metadaten-Map (`[id:patientA, replicate:1, type:normal]`) an erster Position des Ausgabe-Tupels, gefolgt vom Tupel der gepaarten Dateien, wie beabsichtigt.

Natürlich wird dies nur dieses spezifische Dateipaar erfassen und verarbeiten.
Wenn du mit der Verarbeitung mehrerer Paare experimentieren möchtest, kannst du versuchen, Wildcards in das Eingabemuster einzufügen und sehen, was passiert.
Versuche zum Beispiel `data/patientA_rep1_*_R{1,2}_001.fastq.gz` zu verwenden.

### Zusammenfassung

- [`channel.fromFilePairs()` findet und paart automatisch zusammengehörige Dateien](https://www.nextflow.io/docs/latest/reference/channel.html#fromfilepairs)
- Dies vereinfacht die Handhabung von Paired-End-Reads in deiner Pipeline
- Gepaarte Dateien können als `[id, [file1, file2]]`-Tupel gruppiert werden
- Die Extraktion von Metadaten kann aus der ID des Dateipaares statt aus einzelnen Dateien erfolgen

---

## 6. Verwendung von Dateioperationen in Prozessen

Lass uns das jetzt in einem einfachen Prozess zusammenführen, um zu vertiefen, wie Dateioperationen innerhalb eines Nextflow-Prozesses verwendet werden.

Wir stellen dir ein vorgefertigtes Prozessmodul namens `ANALYZE_READS` zur Verfügung, das ein Tupel aus Metadaten und einem Paar Eingabedateien nimmt und diese analysiert.
Wir könnten uns vorstellen, dass dies Sequenz-Alignment, Variant Calling oder jeden anderen Schritt durchführt, der für diesen Datentyp sinnvoll ist.

Los geht's.

### 6.1. Importiere den Prozess und untersuche den Code

Um diesen Prozess im Workflow zu verwenden, müssen wir nur eine Modul-Include-Anweisung vor dem Workflow-Block hinzufügen.

Nimm folgende Änderung am Workflow vor:

=== "Nachher"

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

!!! note

    Die `tag`- und `publishDir`-Direktiven verwenden Closure-Syntax (`{ ... }`) statt String-Interpolation (`"${...}"`).
    Das liegt daran, dass diese Direktiven auf Eingabevariablen (`meta`) verweisen, die erst zur Laufzeit verfügbar sind.
    Die Closure-Syntax schiebt die Auswertung auf, bis der Prozess tatsächlich läuft.

!!! note

    Wir nennen unsere Metadaten-Map konventionsgemäß `meta`.
    Für einen tieferen Einblick in Meta-Maps siehe die Side Quest [Metadata and meta maps](./metadata.md).

### 6.2. Rufe den Prozess im Workflow auf

Jetzt, da der Prozess für den Workflow verfügbar ist, können wir einen Aufruf des `ANALYZE_READS`-Prozesses hinzufügen, um ihn auszuführen.

Um ihn auf unseren Beispieldaten laufen zu lassen, müssen wir zwei Dinge tun:

1. Dem remappten Channel einen Namen geben
2. Einen Aufruf des Prozesses hinzufügen

#### 6.2.1. Benenne den remappten Input-Channel

Wir haben die Mapping-Manipulationen bisher direkt auf den Input-Channel angewendet.
Um die remappten Inhalte an den `ANALYZE_READS`-Prozess zu übergeben (und dies auf eine klare und leicht lesbare Weise zu tun), wollen wir einen neuen Channel namens `ch_samples` erstellen.

Das können wir mit dem [`set`](https://www.nextflow.io/docs/latest/reference/operator.html#set)-Operator tun.

Ersetze im Haupt-Workflow den `.view()`-Operator durch `.set { ch_samples }` und füge eine Zeile hinzu, die testet, dass wir den Channel per Namen referenzieren können.

=== "Nachher"

    ```groovy title="main.nf" linenums="7" hl_lines="14 16-17"
        // Load files with channel.fromFilePairs
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

        // Temporary: peek into ch_samples
        ch_samples.view()
    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="7" hl_lines="14"
        // Load files with channel.fromFilePairs
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

Dies bestätigt, dass wir den Channel nun per Namen referenzieren können.

#### 6.2.2. Rufe den Prozess auf den Daten auf

Lass uns nun tatsächlich den `ANALYZE_READS`-Prozess auf dem `ch_samples`-Channel aufrufen.

Nimm im Haupt-Workflow folgende Code-Änderungen vor:

=== "Nachher"

    ```groovy title="main.nf" linenums="23"
        // Run the analysis
        ANALYZE_READS(ch_samples)
    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="23"
        // Temporary: peek into ch_samples
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

Dieser Prozess ist so eingerichtet, dass er seine Ausgaben in ein `results`-Verzeichnis publiziert, also schau dort mal rein.

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

Der Prozess hat unsere Eingaben genommen und eine neue Datei mit den Patienten-Metadaten erstellt, wie vorgesehen.
Ausgezeichnet!

### 6.3. Viele weitere Patienten einbeziehen

Natürlich verarbeitet dies nur ein einzelnes Dateipaar für einen einzelnen Patienten, was nicht gerade der Durchsatz ist, den du mit Nextflow erhoffst.
Du wirst wahrscheinlich viel mehr Daten gleichzeitig verarbeiten wollen.

Erinnere dich, dass `channel.fromPath()` ein _Glob_ als Eingabe akzeptiert, was bedeutet, dass es beliebig viele Dateien akzeptieren kann, die zum Muster passen.
Wenn wir also alle Patienten einbeziehen möchten, können wir einfach den Eingabe-String modifizieren, um mehr Patienten einzubeziehen, wie bereits früher nebenbei erwähnt.

Tun wir so, als wollten wir so großzügig wie möglich sein.
Nimm folgende Änderungen am Workflow vor:

=== "Nachher"

    ```groovy title="main.nf" linenums="7" hl_lines="2"
        // Load files with channel.fromFilePairs
        ch_files = channel.fromFilePairs('data/*_R{1,2}_001.fastq.gz')
    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="7" hl_lines="2"
        // Load files with channel.fromFilePairs
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

??? abstract "Verzeichnisinhalte"

    ```console
    results
    ├── patientA
    │   └── patientA_stats.txt
    ├── patientB
    │   └── patientB_stats.txt
    └── patientC
        └── patientC_stats.txt
    ```

Erfolg! Wir haben alle Patienten in einem Durchgang analysiert! Richtig?

Vielleicht nicht.
Wenn du genauer hinschaust, haben wir ein Problem: Wir haben zwei Replikate für PatientA, aber nur eine Ausgabedatei!
Wir überschreiben die Ausgabedatei jedes Mal.

### 6.4. Die publizierten Dateien eindeutig machen

Da wir Zugriff auf die Patienten-Metadaten haben, können wir sie verwenden, um die publizierten Dateien eindeutig zu machen, indem wir unterscheidende Metadaten entweder in der Verzeichnisstruktur oder in den Dateinamen selbst einbeziehen.

Nimm folgende Änderung am Workflow vor:

=== "Nachher"

    ```groovy title="modules/analyze_reads.nf" linenums="6"
        publishDir { "results/${meta.type}/${meta.id}/${meta.replicate}" }, mode: 'copy'
    ```

=== "Vorher"

    ```groovy title="modules/analyze_reads.nf" linenums="6"
        publishDir { "results/${meta.id}" }, mode: 'copy'
    ```

Hier zeigen wir die Option, zusätzliche Verzeichnisebenen zu verwenden, um Probentypen und Replikate zu berücksichtigen, aber du könntest auch experimentieren, dies auf Dateinamen-Ebene zu tun.

Führe die Pipeline nun ein letztes Mal aus, aber stelle sicher, dass du zuerst das Ergebnisverzeichnis entfernst, um einen sauberen Arbeitsbereich zu haben:

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

??? abstract "Verzeichnisinhalte"

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

Und da haben wir es, all unsere Metadaten, ordentlich organisiert. Das ist ein Erfolg!

Es gibt noch viel mehr, was du tun kannst, sobald deine Metadaten in eine Map wie diese geladen sind:

1. Organisiere Ausgabeverzeichnisse basierend auf Patientenattributen
2. Treffe Entscheidungen in Prozessen basierend auf Patienteneigenschaften
3. Teile, verbinde und rekombiniere Daten basierend auf Metadatenwerten

Dieses Muster, Metadaten explizit zu halten und an die Daten anzuhängen (anstatt sie in Dateinamen zu kodieren), ist eine zentrale Best Practice in Nextflow, die den Aufbau robuster, wartbarer Analyse-Workflows ermöglicht.
Mehr darüber erfährst du in der Side Quest [Metadata and meta maps](./metadata.md).

### Zusammenfassung

- Die `publishDir`-Direktive kann Ausgaben basierend auf Metadatenwerten organisieren
- Metadaten in Tupeln ermöglichen strukturierte Organisation von Ergebnissen
- Dieser Ansatz schafft wartbare Workflows mit klarer Datenherkunft
- Prozesse können Tupel aus Metadaten und Dateien als Eingabe nehmen
- Die `tag`-Direktive liefert Prozess-Identifikation in Ausführungslogs
- Die Workflow-Struktur trennt Channel-Erstellung von Prozessausführung

---

## Zusammenfassung

In dieser Side Quest hast du gelernt, wie man mit Dateien in Nextflow arbeitet, von grundlegenden Operationen bis zu fortgeschritteneren Techniken für die Handhabung von Dateisammlungen.

Die Anwendung dieser Techniken in deiner eigenen Arbeit wird es dir ermöglichen, effizientere und wartbarere Workflows zu erstellen, insbesondere bei der Arbeit mit großen Mengen von Dateien mit komplexen Benennungskonventionen.

### Wichtige Muster

1.  **Grundlegende Dateioperationen:** Wir haben Path-Objekte mit `file()` erstellt und auf Dateiattribute wie Name, Extension und übergeordnetes Verzeichnis zugegriffen, wobei wir den Unterschied zwischen Strings und Path-Objekten gelernt haben.

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

2.  **Verwendung entfernter Dateien**: Wir haben gelernt, wie man transparent zwischen lokalen und entfernten Dateien mit URIs wechselt, was Nextflows Fähigkeit demonstriert, Dateien aus verschiedenen Quellen zu verarbeiten, ohne die Workflow-Logik zu ändern.

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

3.  **Dateien mit der `fromPath()`-Channel-Factory laden:** Wir haben Channels aus Dateimustern mit `channel.fromPath()` erstellt und ihre Dateiattribute einschließlich Objekttypen angesehen.

    - Einen Channel aus einem Dateimuster erstellen

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

4.  **Extraktion von Patienten-Metadaten aus Dateinamen:** Wir haben `tokenize()` und `replace()` verwendet, um Metadaten aus Dateinamen zu extrahieren und zu strukturieren und sie in organisierte Maps umzuwandeln.

    ```groovy
    def name = file.name.tokenize('_')
    def patientId = name[0]
    def replicate = name[1].replace('rep', '')
    def type = name[2]
    def readNum = name[3].replace('R', '')
    ```

5.  **Vereinfachung mit channel.fromFilePairs:** Wir haben `channel.fromFilePairs()` verwendet, um zusammengehörige Dateien automatisch zu paaren und Metadaten aus gepaarten Datei-IDs zu extrahieren.

    ```groovy
    ch_pairs = channel.fromFilePairs('data/*_R{1,2}_001.fastq.gz')
    ```

6.  **Verwendung von Dateioperationen in Prozessen:** Wir haben Dateioperationen in Nextflow-Prozesse mit korrekter Input-Verarbeitung integriert und `publishDir` verwendet, um Ausgaben basierend auf Metadaten zu organisieren.

    - Eine Meta-Map mit den Prozess-Inputs verknüpfen

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

Kehre zum [Menü der Side Quests](./index.md) zurück oder klicke auf die Schaltfläche unten rechts auf der Seite, um zum nächsten Thema in der Liste zu gelangen.
