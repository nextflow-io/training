# Wesentliche Nextflow-Skriptmuster

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Nextflow ist eine Programmiersprache, die auf der Java Virtual Machine läuft. Obwohl Nextflow auf [Groovy](http://groovy-lang.org/) basiert und viel von dessen Syntax teilt, ist Nextflow mehr als nur „Groovy mit Erweiterungen" – es ist eine eigenständige Sprache mit einer vollständig spezifizierten [Syntax](https://nextflow.io/docs/latest/reference/syntax.html) und [Standardbibliothek](https://nextflow.io/docs/latest/reference/stdlib.html).

Du kannst viel Nextflow schreiben, ohne über die grundlegende Syntax für Variablen, Maps und Listen hinauszugehen. Die meisten Nextflow-Tutorials konzentrieren sich auf die Workflow-Orchestrierung (Channels, Prozesse und Datenfluss), und damit kommst du überraschend weit.

Wenn du jedoch Daten manipulieren, komplexe Dateinamen parsen, bedingte Logik implementieren oder robuste Produktions-Workflows erstellen musst, hilft es, zwei verschiedene Aspekte deines Codes zu betrachten: **Dataflow** (Channels, Operatoren, Prozesse und Workflows) und **Scripting** (der Code innerhalb von Closures, Funktionen und Prozess-Scripts). Obwohl diese Unterscheidung etwas willkürlich ist – es ist alles Nextflow-Code – bietet sie ein nützliches mentales Modell, um zu verstehen, wann du deinen Pipeline orchestrierst und wann du Daten manipulierst. Die Beherrschung beider Aspekte verbessert deine Fähigkeit, klare, wartbare Workflows zu schreiben, dramatisch.

### Lernziele

Dieses Side Quest führt dich von grundlegenden Konzepten zu produktionsreifen Mustern.
Wir werden einen einfachen CSV-lesenden Workflow in einen ausgefeilten Bioinformatik-Pipeline verwandeln und ihn Schritt für Schritt durch realistische Herausforderungen weiterentwickeln:

- **Grenzen verstehen:** Unterscheide zwischen Dataflow-Operationen und Scripting und verstehe, wie sie zusammenarbeiten
- **Datenmanipulation:** Extrahiere, transformiere und wähle Maps und Collections mit leistungsstarken Operatoren aus
- **String-Verarbeitung:** Parse komplexe Dateibenennungsschemata mit Regex-Mustern und meistere Variable Interpolation
- **Wiederverwendbare Funktionen:** Extrahiere komplexe Logik in benannte Funktionen für sauberere, wartbarere Workflows
- **Dynamische Logik:** Erstelle Prozesse, die sich an verschiedene Eingabetypen anpassen, und verwende Closures für dynamische Ressourcenzuweisung

### Voraussetzungen

- Grundlegende Kenntnisse von Nextflow-Channels und Prozessen
- Verständnis der DSL2-Syntax
- Grundlegende Programmierkenntnisse (Variablen, Bedingungen, Schleifen)
- Übung mit der Ausführung von Nextflow-Workflows

### Erste Schritte

Für dieses Beispiel werden wir eine kleine Sammlung von Nextflow-Skripten verwenden. Bitte lade sie herunter, indem du diesen Befehl in deinem Terminal ausführst:

```bash
nextflow clone nextflow-io/training -r 25.10.2 -d side-quests/essential_scripting_patterns
cd side-quests/essential_scripting_patterns
```

## 1. Datenfluss vs. Scripting

Wie eingangs erwähnt, ist es beim Schreiben von Nextflow nützlich, ein mentales Modell von zwei verschiedenen Aspekten des Codes zu haben: **Datenfluss-Orchestrierung** (Channels, Operatoren, Prozesse) und **Scripting** (der Code in Closures, Funktionen und Prozess-Scripts).

- **Datenfluss:** Orchestrierung der Pipeline auf hoher Ebene

  - Channel-Erstellung und -Transformation
  - Ausführungsreihenfolge der Prozesse
  - Hauptpipeline-Logik (Verzweigungen, Filter, Zusammenfügungen)

- **Scripting:** Code, der bestimmte Daten manipuliert
  - Closures innerhalb von Operatoren (`map`, `filter`, `collect`)
  - Prozess-Script-Blöcke
  - Funktionen, die von anderen Scripten aufgerufen werden

Fangen wir an, indem wir ein einfaches Workflow-Skript erstellen, das eine CSV-Datei mit Metadaten zu Proben liest und die Metadaten verarbeitet. Nimm dir einen Moment, um zu verstehen, wie der Code arbeitet:

```groovy title="main.nf" linenums="1"
#!/usr/bin/env nextflow

workflow {
    // DATENFLUSS: Channel-Erstellung und -Transformation
    ch_samples = channel.fromPath("./data/samples.csv")
        .splitCsv(header: true)
        .map { row -> // SCRIPTING beginnt hier innerhalb der Closure
            def sample_meta = [
                id: row.sample_id.toLowerCase(),
                organism: row.organism,
                tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                depth: row.sequencing_depth.toInteger(),
                quality: row.quality_score.toDouble()
            ]
            def fastq_path = file(row.file_path)
            return tuple(sample_meta, fastq_path)
        } // SCRIPTING endet hier

    // DATENFLUSS: Kanalvisualisierung
    ch_samples.view()
}
```

Lass uns diesen Code durch die Linse des Datenfluss- vs. Scripting-Modells analysieren:

### 1.1. Closures und Bedingungen

Die obige Transformation `.map { row -> ... }` verwendet eine **Groovy-Closure**, die an den `.map()`-Operator übergeben wird. Eine Closure ist eine wiederverwendbare Codeblock, der Variablen deklarieren, Berechnungen durchführen und Werte zurückgeben kann. Das ist ein Kernkonzept beim Scripting in Nextflow.

Der Syntax ist:

```groovy title="closure-syntax"
{ Parameter ->
    // Codeblock
    return Wert // optional, gibt standardmäßig den letzten Ausdruck zurück
}
```

Die Parameter sind optional. Wenn du keine angibt, erhält die Closure standardmäßig einen Parameter namens `it`. So könnte das obige Beispiel auch geschrieben werden als:

```groovy title="implicit-it"
.map {
    def sample_meta = [
        id: it.sample_id.toLowerCase(),
        // ...andere Felder...
    ]
    // ...Rest der Closure...
}
```

Allerdings ist die explizite Benennung des Parameters (`row` in unserem Beispiel) oft klarer, besonders bei verschachtelter Logik.

Eine mächtige Eigenschaft von Closures ist ihre Fähigkeit, bedingte Logik zu verwenden. Innerhalb des Scripting-Teils können wir Standard-Programmierstrukturen wie `if/else`, den ternären Operator (`?:`), Schleifen und mehr verwenden.

Erweitern wir unser Beispiel, um eine bedingte Logik mit dem ternären Operator hinzuzufügen:

```groovy title="main.nf" linenums="6" hl_lines="6-10"
.map { row ->
    def sample_meta = [
        id: row.sample_id.toLowerCase(),
        organism: row.organism,
        tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
        depth: row.sequencing_depth.toInteger(),
        quality: row.quality_score.toDouble()
    ]
    def fastq_path = file(row.file_path)
    def priority = sample_meta.quality > 40 ? 'high' : 'normal'
    return tuple(sample_meta + [priority: priority], fastq_path)
}
```

Hier haben wir eine bedingte Logik hinzugefügt, die eine Prioritätskennzeichnung zur Probenmetadaten hinzufügt, basierend auf dem Qualitätswert. Der ternäre Operator (`?:`) ist eine kurze Form für `if/else` - wenn die Bedingung vor dem `?` wahr ist, wird der Wert nach dem `?` zurückgegeben, sonst der Wert nach dem `:`.

Wir verwenden auch die Syntax `sample_meta + [priority: priority]`, um ein neues Feld zum `sample_meta`-Map hinzuzufügen (Groovy erlaubt das Addieren von Maps).

Führe diesen aktualisierten Code aus:

```bash
nextflow run main.nf
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [cheeky_shirley] DSL2 - revision: 66c7a1a983
    [[id:sample_001, organism:human, tissue:liver, depth:30000000, quality:38.5, priority:normal], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_001_S1_L001_R1_001.fastq]
    [[id:sample_002, organism:mouse, tissue:brain, depth:25000000, quality:35.2, priority:normal], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_002_S2_L001_R1_001.fastq]
    [[id:sample_003, organism:human, tissue:kidney, depth:45000000, quality:42.1, priority:high], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_003_S3_L001_R1_001.fastq]
    ```

Du kannst sehen, dass die Probe mit der hohen Qualitätsbewertung als `priority:high` gekennzeichnet wurde, während die anderen als `priority:normal` gekennzeichnet wurden.

### 1.2. Collections vs. Channels

Eine häufige Verwirrung bei Nextflow entsteht durch den Unterschied zwischen Channel-Operatoren und Collection-Methoden. Beide können ähnlich benannt sein (z.B. `collect`), aber sie funktionieren unterschiedlich:

- **Channel-Operatoren:** Transformieren die Struktur eines Channels im Datenfluss (`collect`, `map`, `filter` etc. nach einem `.`)
- **Collection-Methoden:** Manipulieren die Daten innerhalb einer Collection (Liste, Map) im Scripting

Betrachte diese zwei Codeblöcke:

```groovy title="Channel Operator (Datenfluss)"
ch_samples = channel.of(1, 2, 3, 4)
    .collect()
    .view()  // [[1, 2, 3, 4]] - gibt eine Liste von allen Werten zurück
```

```groovy title="Collection Method (Scripting)"
ch_samples = channel.of(1, 2, 3, 4)
    .map { values ->
        values.collect { it * 2 }
    }
    .view()  // 2, 4, 6, 8 - verarbeitet jedes Element
```

Das ist ein häufiges Missverständnis! Die Channel-Version von `.collect()` (erste Beispiel) sammelt alle Werte in eine einzige Liste, wodurch der Channel seine asynchrone, streaming Natur verliert. Die Listen-Version von `.collect()` (zweites Beispiel) wendet eine Transformation auf jedes Element der Liste an.

!!! warning "Vorsicht bei gleichnamigen Methoden"

    Wenn du jemals einen Fehler wie "Method X is not applicable for class nextflow.Channel" siehst, könnte das daran liegen, dass du versuchst, eine Collection-Methode auf einen Channel anzuwenden (oder umgekehrt).

Für dieses Side Quest haben wir vorbereitete Module in `./modules/`. Lass uns einen einfachen Prozess integrieren, der `fastp` verwendet, um Qualitätsberichte für unsere Proben zu erstellen. Ändere dein Skript wie folgt:

```groovy title="main.nf" linenums="1"
include { FASTP } from './modules/fastp.nf'

workflow {
    ch_samples = channel.fromPath("./data/samples.csv")
        .splitCsv(header: true)
        .map { row ->
            def sample_meta = [
                id: row.sample_id.toLowerCase(),
                organism: row.organism,
                tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                depth: row.sequencing_depth.toInteger(),
                quality: row.quality_score.toDouble()
            ]
            def fastq_path = file(row.file_path)
            def priority = sample_meta.quality > 40 ? 'high' : 'normal'
            return tuple(sample_meta + [priority: priority], fastq_path)
        }

    ch_fastp = FASTP(ch_samples)
    ch_fastp.view()
}
```

Führe es aus:

```bash
nextflow run main.nf
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W  ~  version 25.10.2
    Launching `main.nf` [nice_kirch] DSL2 - revision: 66c7a1a983
    executor >  local (3)
    [0c/b8cc75] process > FASTP (1) [100%] 3 of 3 ✔
    [id:sample_001, organism:human, tissue:liver, depth:30000000, quality:38.5, priority:normal, fastp_json:/workspaces/training/side-quests/essential_scripting_patterns/results/reports/sample_001.fastp.json, fastp_html:/workspaces/training/side-quests/essential_scripting_patterns/results/reports/sample_001.fastp.html, fastp_log:/workspaces/training/side-quests/essential_scripting_patterns/results/logs/sample_001.fastp.log]
    [id:sample_002, organism:mouse, tissue:brain, depth:25000000, quality:35.2, priority:normal, fastp_json:/workspaces/training/side-quests/essential_scripting_patterns/results/reports/sample_002.fastp.json, fastp_html:/workspaces/training/side-quests/essential_scripting_patterns/results/reports/sample_002.fastp.html, fastp_log:/workspaces/training/side-quests/essential_scripting_patterns/results/logs/sample_002.fastp.log]
    [id:sample_003, organism:human, tissue:kidney, depth:45000000, quality:42.1, priority:high, fastp_json:/workspaces/training/side-quests/essential_scripting_patterns/results/reports/sample_003.fastp.json, fastp_html:/workspaces/training/side-quests/essential_scripting_patterns/results/reports/sample_003.fastp.html, fastp_log:/workspaces/training/side-quests/essential_scripting_patterns/results/logs/sample_003.fastp.log]
    ```

Der Prozess hat erfolgreich für jede Eingabeprobe einen FASTP-Bericht erstellt und die Metadaten mit den Ausgabedatei-Pfaden ergänzt.

Lass uns diesen Workflow weiter ausbauen, indem wir einen Prozess hinzufügen, der einfache Berichte für jede Probe generiert. Wir werden einen neuen Prozess aus der vorbereiteten Modul importieren:

```groovy title="main.nf" linenums="1" hl_lines="2"
include { FASTP } from './modules/fastp.nf'
include { GENERATE_REPORT } from './modules/generate_report.nf'

workflow {
    ch_samples = channel.fromPath("./data/samples.csv")
        .splitCsv(header: true)
        .map { row ->
            def sample_meta = [
                id: row.sample_id.toLowerCase(),
                organism: row.organism,
                tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                depth: row.sequencing_depth.toInteger(),
                quality: row.quality_score.toDouble()
            ]
            def fastq_path = file(row.file_path)
            def priority = sample_meta.quality > 40 ? 'high' : 'normal'
            return tuple(sample_meta + [priority: priority], fastq_path)
        }

    ch_fastp = FASTP(ch_samples)
    GENERATE_REPORT(ch_samples)
}
```

Jetzt haben wir zwei Prozesse: `FASTP` für die Qualitätskontrolle und `GENERATE_REPORT` für Berichtsgenerierung. Beachte, wie wir beide Prozesse mit dem gleichen Eingangs-Channel verbinden - dies ist ein Vorteil des Dataflow-Modells von Nextflow, das die Parallelisierung der Prozesse ohne explizite Parallelitätskontrolle oder Verzweigungslogik ermöglicht.

### Fazit

In diesem Abschnitt hast du gelernt:

- **Datenfluss vs. Scripting:** Den Unterschied zwischen Dataflow (Channel-Orchestrierung) und Scripting (Datenmanipulation) erkennen
- **Closures:** Codeblöcke für Channel-Operatoren schreiben und Parameter explizit benennen
- **Bedingungslogik:** Den ternären Operator (`?:`) für kurze if/else-Ausdrücke verwenden
- **Collections vs. Channels:** Den Unterschied zwischen Channel-Operatoren und Collection-Methoden verstehen

Diese Unterscheidung bildet die Grundlage für die Erstellung sauberer, wartbarer Workflows. Als nächstes werden wir uns mit fortgeschrittenen String-Verarbeitungstechniken befassen, um komplexe Dateinamen zu parsen.

---

## 2. Fortgeschrittene String-Verarbeitung

Bioinformatik-Workflows müssen oft mit komplexen Dateinamenskonventionen und -formaten arbeiten. In diesem Abschnitt werden wir reguläre Ausdrücke (Regex) und andere Methoden zur Stringverarbeitung verwenden, um Metadaten aus FASTQ-Dateinamen zu extrahieren.

### 2.1. Regex für das Parsen von Dateinamen

Ein häufiges Muster in der Bioinformatik ist die Extraktion von Probeninformationen aus standardisierten Dateinamen. Nehmen wir an, unsere FASTQ-Dateien folgen einer Illumina-ähnlichen Benennungskonvention:

```
SAMPLE_001_S1_L001_R1_001.fastq
```

wobei:

- `SAMPLE_001` ist die Proben-ID
- `S1` ist die Probennummer
- `L001` ist die Lane-Nummer
- `R1` gibt an, ob es sich um Read 1 oder Read 2 handelt (für PE-Sequenzierung)
- `001` ist eine fortlaufende Nummer für den Chunk

Diese Informationen aus dem Dateinamen zu extrahieren, ermöglicht eine bessere Organisation und Nachverfolgung von Proben in einer Pipeline. Erweitern wir unseren `.map()`-Operator, um Informationen aus dem Dateinamen zu extrahieren:

=== "Nach"

    ```groovy title="main.nf" linenums="11" hl_lines="6-12"
            def fastq_path = file(row.file_path)

            def m = (fastq_path.name =~ /^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$/)
            def file_meta = m ? [
                sample_num: m[0][2].toInteger(),
                lane: m[0][3],
                read: m[0][4],
                chunk: m[0][5]
            ] : [:]

            def priority = sample_meta.quality > 40 ? 'high' : 'normal'
            return tuple(sample_meta + file_meta + [priority: priority], fastq_path)
    ```

=== "Vor"

    ```groovy title="main.nf" linenums="11"
            def fastq_path = file(row.file_path)
            def priority = sample_meta.quality > 40 ? 'high' : 'normal'
            return tuple(sample_meta + [priority: priority], fastq_path)
    ```

Lass uns verstehen, was dieser Code tut:

1. `(fastq_path.name =~ /^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$/)` - Dies erstellt ein Matcher-Objekt mit einem Regex-Muster

   - `=~` ist der Regex-"Findall"-Operator in Groovy
   - Das Muster ist in Slashes eingeschlossen (`/pattern/`)
   - `()` erfassen Gruppen, auf die wir später zugreifen können

2. `def file_meta = m ? [...] : [:]` - Dies verwendet einen ternären Operator, um ein Map zu erstellen, wenn der Regex übereinstimmt, oder ein leeres Map, wenn nicht

3. `m[0][2].toInteger()` - Dies greift auf den Inhalt der zweiten erfassten Gruppe im ersten Match zu (und konvertiert ihn in einen Integer)

4. Schließlich fügen wir die extrahierten Metadaten zu unserem bestehenden `sample_meta`-Map hinzu

Führe diese aktualisierte Version aus:

```bash
nextflow run main.nf
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W  ~  version 25.10.2
    Launching `main.nf` [reverent_caravaggio] DSL2 - revision: 3a8c8a48b2
    executor >  local (3)
    [aa/e09a1e] process > FASTP (1)           [100%] 3 of 3 ✔
    [1b/3fb1a8] process > GENERATE_REPORT (1) [100%] 3 of 3 ✔
    ```

Schau dir einen der erzeugten Berichte an:

```bash
cat results/reports/sample_001_report.txt
```

??? success "Befehlsausgabe"

    ```console
    Processing /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_001_S1_L001_R1_001.fastq
    Sample: sample_001
    ```

Die Berichte sind etwas karg, aber wir können sie verbessern. Bearbeite `modules/generate_report.nf` um mehr Informationen einzubeziehen:

=== "Nach"

    ```groovy title="modules/generate_report.nf" linenums="1" hl_lines="12-16"
    process GENERATE_REPORT {
        publishDir "results/reports"

        input:
        tuple val(meta), path(reads)

        output:
        path "${meta.id}_report.txt"

        script:
        """
        echo "Processing ${reads}" > ${meta.id}_report.txt
        echo "Sample: ${meta.id}" >> ${meta.id}_report.txt
        echo "Organism: ${meta.organism}" >> ${meta.id}_report.txt
        echo "Quality Score: ${meta.quality}" >> ${meta.id}_report.txt
        echo "Lane: ${meta.lane}" >> ${meta.id}_report.txt
        echo "Read: ${meta.read}" >> ${meta.id}_report.txt
        """
    }
    ```

=== "Vor"

    ```groovy title="modules/generate_report.nf" linenums="1"
    process GENERATE_REPORT {
        publishDir "results/reports"

        input:
        tuple val(meta), path(reads)

        output:
        path "${meta.id}_report.txt"

        script:
        """
        echo "Processing ${reads}" > ${meta.id}_report.txt
        echo "Sample: ${meta.id}" >> ${meta.id}_report.txt
        """
    }
    ```

Führe es erneut aus und schaue dir einen Bericht an:

```bash
nextflow run main.nf
cat results/reports/sample_001_report.txt
```

??? success "Befehlsausgabe"

    ```console
    Processing /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_001_S1_L001_R1_001.fastq
    Sample: sample_001
    Organism: human
    Quality Score: 38.5
    Lane: 001
    Read: R1
    ```

Super! Jetzt enthält unser Bericht Informationen, die sowohl aus der CSV als auch aus dem Dateinamen extrahiert wurden.

### 2.2. Variable-Interpolation in Skripten

Eine wichtige Unterscheidung in Nextflow ist, wann verschiedene Arten von Variablen in Strings und Skripten interpoliert werden. Variablen in Nextflow können auf drei Arten verwendet werden:

1. **Nextflow-Variablen**: `${var}` wird zur Kompilierungszeit (während des Workflows) verarbeitet
2. **Shell-Umgebungsvariablen**: `\${VAR}` muss in Prozess-Scripts mit einem Backslash escaped werden
3. **Shell-Befehlsersetzung**: `\$(command)` muss mit einem Backslash escaped werden

Lass uns einen Workflow erstellen, der einige fortgeschrittenere Stringverarbeitungstechniken demonstriert:

```groovy title="advanced_string_processing.nf" linenums="1"
#!/usr/bin/env nextflow

include { FASTP } from './modules/fastp.nf'
include { GENERATE_REPORT } from './modules/generate_report.nf'

workflow {
    ch_samples = channel.fromPath("./data/samples.csv")
        .splitCsv(header: true)
        .map { row ->
            def sample_meta = [
                id: row.sample_id.toLowerCase(),
                organism: row.organism,
                tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                depth: row.sequencing_depth.toInteger(),
                quality: row.quality_score.toDouble()
            ]
            def fastq_path = file(row.file_path)

            def m = (fastq_path.name =~ /^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$/)
            def file_meta = m ? [
                sample_num: m[0][2].toInteger(),
                lane: m[0][3],
                read: m[0][4],
                chunk: m[0][5]
            ] : [:]

            def priority = sample_meta.quality > 40 ? 'high' : 'normal'
            return tuple(sample_meta + file_meta + [priority: priority], fastq_path)
        }

    ch_fastp = FASTP(ch_samples)
    GENERATE_REPORT(ch_samples)
}
```

Führe diesen Workflow aus und schau dir die erzeugten Berichte in `results/reports/` an. Sie sollten grundlegende Informationen über jede Probe enthalten.

<!-- TODO: den Befehl zum Ausführen hinzufügen -->

??? success "Befehlsausgabe"

    ```console
    <!-- TODO: Ausgabe -->
    ```

Aber was, wenn wir Informationen darüber hinzufügen wollen, wann und wo die Verarbeitung stattgefunden hat? Lass uns den Prozess ändern, um **Shell**-Variablen und ein wenig Befehlsersetzung zu verwenden, um den aktuellen Benutzer, Hostnamen und Datum in den Bericht aufzunehmen:

=== "Nach"

    ```groovy title="modules/generate_report.nf" linenums="10" hl_lines="5-7"
        script:
        """
        echo "Processing ${reads}" > ${meta.id}_report.txt
        echo "Sample: ${meta.id}" >> ${meta.id}_report.txt
        echo "Processed by: ${USER}" >> ${meta.id}_report.txt
        echo "Hostname: $(hostname)" >> ${meta.id}_report.txt
        echo "Date: $(date)" >> ${meta.id}_report.txt
        """
    ```

=== "Vor"

    ```groovy title="modules/generate_report.nf" linenums="10"
        script:
        """
        echo "Processing ${reads}" > ${meta.id}_report.txt
        echo "Sample: ${meta.id}" >> ${meta.id}_report.txt
        """
    ```

Wenn du dies ausführst, wirst du einen Fehler bemerken - Nextflow versucht, `${USER}` als Nextflow-Variable zu interpretieren, die nicht existiert.

??? failure "Befehlsausgabe"

    ```console
    Error modules/generate_report.nf:15:27: `USER` is not defined
    │  15 |     echo "Processed by: ${USER}" >> ${meta.id}_report.txt
    ╰     |                           ^^^^

    ERROR ~ Script compilation failed
    ```

Wir müssen es escapen, damit Bash es stattdessen verarbeiten kann.

Korrigiere dies, indem du die Shell-Variablen und Befehlsersetzungen mit einem Backslash (`\`) escapest:

=== "Nach"

    ```groovy title="modules/generate_report.nf" linenums="10" hl_lines="5-7"
        script:
        """
        echo "Processing ${reads}" > ${meta.id}_report.txt
        echo "Sample: ${meta.id}" >> ${meta.id}_report.txt
        echo "Processed by: \${USER}" >> ${meta.id}_report.txt
        echo "Hostname: \$(hostname)" >> ${meta.id}_report.txt
        echo "Date: \$(date)" >> ${meta.id}_report.txt
        """
    ```

=== "Vor"

    ```groovy title="modules/generate_report.nf" linenums="10"
        script:
        """
        echo "Processing ${reads}" > ${meta.id}_report.txt
        echo "Sample: ${meta.id}" >> ${meta.id}_report.txt
        echo "Processed by: ${USER}" >> ${meta.id}_report.txt
        echo "Hostname: $(hostname)" >> ${meta.id}_report.txt
        echo "Date: $(date)" >> ${meta.id}_report.txt
        """
    ```

Jetzt funktioniert es! Der Backslash (`\`) sagt Nextflow "interpretiere das nicht, gib es an Bash weiter."

### Fazit

In diesem Abschnitt hast du **String-Verarbeitungstechniken** gelernt:

- **Reguläre Ausdrücke für das Parsen von Dateien**: Verwendung des `=~`-Operators und Regex-Muster (`~/pattern/`), um Metadaten aus komplexen Dateibenennungskonventionen zu extrahieren
- **Dynamische Skripterstellung**: Verwendung von bedingter Logik (if/else, ternäre Operatoren), um verschiedene Skript-Strings basierend auf Eingabeeigenschaften zu generieren
- **Variable-Interpolation**: Verständnis, wann Nextflow Strings interpretiert vs. wann die Shell es tut
  - `${var}` - Nextflow-Variablen (zur Kompilierungszeit des Workflows von Nextflow interpoliert)
  - `\${var}` - Shell-Umgebungsvariablen (escaped, zur Laufzeit an Bash weitergegeben)
  - `\$(cmd)` - Shell-Befehlsersetzung (escaped, zur Laufzeit von Bash ausgeführt)

Diese String-Verarbeitungs- und Erzeugungsmuster sind unerlässlich für den Umgang mit den verschiedenen Dateiformaten und Namenskonventionen, denen du in realen Bioinformatik-Workflows begegnen wirst.

---

## 3. Wiederverwendbare Funktionen erstellen

Komplexe Workflow-Logik inline in Channel-Operatoren oder Prozessdefinitionen reduziert die Lesbarkeit und Wartbarkeit. **Funktionen** ermöglichen es dir, diese Logik in benannte, wiederverwendbare Komponenten zu extrahieren.

Unsere map-Operation ist lang und komplex geworden. Lass uns sie in eine wiederverwendbare Funktion mit dem Schlüsselwort `def` extrahieren.

Um zu veranschaulichen, wie das mit unserem bestehenden Workflow aussieht, nimm die folgende Änderung vor, indem du `def` verwendest, um eine wiederverwendbare Funktion namens `separateMetadata` zu definieren:

=== "Nach"

    ```groovy title="main.nf" linenums="1" hl_lines="4-24 29"
    include { FASTP } from './modules/fastp.nf'
    include { GENERATE_REPORT } from './modules/generate_report.nf'

    def separateMetadata(row) {
        def sample_meta = [
            id: row.sample_id.toLowerCase(),
            organism: row.organism,
            tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
            depth: row.sequencing_depth.toInteger(),
            quality: row.quality_score.toDouble()
        ]
        def fastq_path = file(row.file_path)

        def m = (fastq_path.name =~ /^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$/)
        def file_meta = m ? [
            sample_num: m[0][2].toInteger(),
            lane: m[0][3],
            read: m[0][4],
            chunk: m[0][5]
        ] : [:]

        def priority = sample_meta.quality > 40 ? 'high' : 'normal'
        return tuple(sample_meta + file_meta + [priority: priority], fastq_path)
    }

    workflow {
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map{ row -> separateMetadata(row) }

        ch_fastp = FASTP(ch_samples)
        GENERATE_REPORT(ch_samples)
    }
    ```

=== "Vor"

    ```groovy title="main.nf" linenums="1" hl_lines="7-27"
    include { FASTP } from './modules/fastp.nf'
    include { GENERATE_REPORT } from './modules/generate_report.nf'

    workflow {
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                def fastq_path = file(row.file_path)

                def m = (fastq_path.name =~ /^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$/)
                def file_meta = m ? [
                    sample_num: m[0][2].toInteger(),
                    lane: m[0][3],
                    read: m[0][4],
                    chunk: m[0][5]
                ] : [:]

                def priority = sample_meta.quality > 40 ? 'high' : 'normal'
                return tuple(sample_meta + file_meta + [priority: priority], fastq_path)
            }

        ch_fastp = FASTP(ch_samples)
        GENERATE_REPORT(ch_samples)
    }
    ```

Durch das Extrahieren dieser Logik in eine Funktion haben wir die eigentliche Workflow-Logik auf etwas viel Klareres reduziert:

```groovy title="minimal workflow"
    ch_samples = channel.fromPath("./data/samples.csv")
        .splitCsv(header: true)
        .map{ row -> separateMetadata(row) }

    ch_fastp = FASTP(ch_samples)
    GENERATE_REPORT(ch_samples)
```

Dies macht die Workflow-Logik viel einfacher zu lesen und auf einen Blick zu verstehen. Die Funktion `separateMetadata` kapselt die gesamte komplexe Logik zum Parsen und Anreichern von Metadaten und macht sie wiederverwendbar und testbar.

Führe den Workflow aus, um sicherzustellen, dass er immer noch funktioniert:

```bash
nextflow run main.nf
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [admiring_panini] DSL2 - revision: 8cc832e32f

    executor >  local (6)
    [8c/2e3f91] process > FASTP (3)           [100%] 3 of 3 ✔
    [7a/1b4c92] process > GENERATE_REPORT (3) [100%] 3 of 3 ✔
    ```

Die Ausgabe sollte zeigen, dass beide Prozesse erfolgreich abgeschlossen wurden. Der Workflow ist jetzt viel sauberer und leichter zu warten, wobei die gesamte komplexe Metadatenverarbeitungslogik in der Funktion `separateMetadata` gekapselt ist.

### Fazit

In diesem Abschnitt hast du **Funktionserstellung** gelernt:

- **Definieren von Funktionen mit `def`**: Das Schlüsselwort für die Erstellung benannter Funktionen (wie `def` in Python oder `function` in JavaScript)
- **Funktionsbereich**: Funktionen, die auf Skriptebene definiert sind, sind in deinem gesamten Nextflow-Workflow zugänglich
- **Rückgabewerte**: Funktionen geben automatisch den letzten Ausdruck zurück, oder verwenden explizit `return`
- **Saubererer Code**: Das Extrahieren komplexer Logik in Funktionen ist eine grundlegende Software-Engineering-Praxis in jeder Sprache

Als nächstes werden wir untersuchen, wie Closures in Prozessdirektiven für dynamische Ressourcenzuweisung verwendet werden können.

---

## 4. Dynamische Ressourcendirektiven mit Closures

Bisher haben wir Scripting im `script`-Block von Prozessen verwendet. Aber **Closures** (in Abschnitt 1.1 eingeführt) sind auch unglaublich nützlich in Prozessdirektiven, besonders für dynamische Ressourcenzuweisung. Fügen wir unserem FASTP-Prozess Ressourcendirektiven hinzu, die sich an die Probeneigenschaften anpassen.

### 4.1. Probenspezifische Ressourcenzuweisung

Derzeit verwendet unser FASTP-Prozess Standardressourcen. Machen wir ihn intelligenter, indem wir mehr CPUs für Proben mit hoher Tiefe zuweisen. Bearbeite `modules/fastp.nf`, um eine dynamische `cpus`-Direktive und eine statische `memory`-Direktive einzufügen:

=== "Nach"

    ```groovy title="modules/fastp.nf" linenums="1" hl_lines="4-5"
    process FASTP {
        container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

        cpus { meta.depth > 40000000 ? 2 : 1 }
        memory 2.GB

        input:
        tuple val(meta), path(reads)
    ```

=== "Vor"

    ```groovy title="modules/fastp.nf" linenums="1"
    process FASTP {
        container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

        input:
        tuple val(meta), path(reads)
    ```

Die Closure `{ meta.depth > 40000000 ? 2 : 1 }` verwendet den **ternären Operator** (in Abschnitt 1.1 behandelt) und wird für jede Aufgabe ausgewertet, was eine probespezifische Ressourcenzuweisung ermöglicht. Proben mit hoher Tiefe (>40M Reads) erhalten 2 CPUs, während andere 1 CPU erhalten.

!!! note "Zugriff auf Eingabevariablen in Direktiven"

    Die Closure kann auf beliebige Eingabevariablen (wie `meta` hier) zugreifen, da Nextflow diese Closures im Kontext jeder Task-Ausführung auswertet.

Führe den Workflow erneut aus mit der Option `-ansi-log false`, um die Task-Hashes besser sehen zu können.

```bash
nextflow run main.nf -ansi-log false
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W  ~  version 25.10.2
    Launching `main.nf` [fervent_albattani] DSL2 - revision: fa8f249759
    [bd/ff3d41] Submitted process > FASTP (2)
    [a4/a3aab2] Submitted process > FASTP (1)
    [48/6db0c9] Submitted process > FASTP (3)
    [ec/83439d] Submitted process > GENERATE_REPORT (3)
    [bd/15d7cc] Submitted process > GENERATE_REPORT (2)
    [42/699357] Submitted process > GENERATE_REPORT (1)
    ```

Du kannst den genauen `docker`-Befehl überprüfen, der ausgeführt wurde, um die CPU-Zuweisung für eine bestimmte Aufgabe zu sehen:

```console title="Check docker command"
cat work/48/6db0c9e9d8aa65e4bb4936cd3bd59e/.command.run | grep "docker run"
```

Du solltest so etwas wie folgendes sehen:

```bash title="docker command"
    docker run -i --cpu-shares 4096 --memory 2048m -e "NXF_TASK_WORKDIR" -v /workspaces/training/side-quests/essential_scripting_patterns:/workspaces/training/side-quests/essential_scripting_patterns -w "$NXF_TASK_WORKDIR" --name $NXF_BOXID community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690 /bin/bash -ue /workspaces/training/side-quests/essential_scripting_patterns/work/48/6db0c9e9d8aa65e4bb4936cd3bd59e/.command.sh
```

In diesem Beispiel haben wir ein Beispiel gewählt, das 2 CPUs angefordert hat (`--cpu-shares 2048`), weil es sich um eine Probe mit hoher Tiefe handelte, aber du solltest je nach Probentiefe unterschiedliche CPU-Zuweisungen sehen. Probiere dies auch für die anderen Aufgaben aus.

### 4.2. Wiederholungsstrategien

Ein weiteres leistungsstarkes Muster ist die Verwendung von `task.attempt` für Wiederholungsstrategien. Um zu zeigen, warum dies nützlich ist, werden wir zunächst die Speicherzuweisung für FASTP auf weniger reduzieren, als es benötigt. Ändere die `memory`-Direktive in `modules/fastp.nf` auf `1.GB`:

=== "Nach"

    ```groovy title="modules/fastp.nf" linenums="1" hl_lines="5"
    process FASTP {
        container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

        cpus { meta.depth > 40000000 ? 4 : 2 }
        memory 1.GB

        input:
        tuple val(meta), path(reads)
    ```

=== "Vor"

    ```groovy title="modules/fastp.nf" linenums="1" hl_lines="5"
    process FASTP {
        container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

        cpus { meta.depth > 40000000 ? 4 : 2 }
        memory 2.GB

        input:
        tuple val(meta), path(reads)
    ```

... und führe den Workflow erneut aus:

```bash
nextflow run main.nf
```

??? failure "Befehlsausgabe"

    ```console hl_lines="2 11"
    Command exit status:
      137

    Command output:
      (empty)

    Command error:
      Detecting adapter sequence for read1...
      No adapter detected for read1

      .command.sh: line 7:   101 Killed                  fastp --in1 SAMPLE_002_S2_L001_R1_001.fastq --out1 sample_002_trimmed.fastq.gz --json sample_002.fastp.json --html sample_002.fastp.html --thread 2
    ```

Dies zeigt an, dass der Prozess aufgrund der Überschreitung der Speicherlimits beendet wurde.

Das ist ein sehr häufiges Szenario in realen Workflows - manchmal weißt du einfach nicht, wie viel Speicher eine Aufgabe benötigen wird, bis du sie ausführst.

Um unseren Workflow robuster zu machen, können wir eine Wiederholungsstrategie implementieren, die die Speicherzuweisung bei jedem Versuch erhöht, wieder mit einer Groovy-Closure. Ändere die `memory`-Direktive, um den Basisspeicher mit `task.attempt` zu multiplizieren, und füge die Direktiven `errorStrategy 'retry'` und `maxRetries 2` hinzu:

=== "Nach"

    ```groovy title="modules/fastp.nf" linenums="1" hl_lines="5-7"
    process FASTP {
        container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

        cpus { meta.depth > 40000000 ? 4 : 2 }
        memory { 1.GB * task.attempt }
        errorStrategy 'retry'
        maxRetries 2

        input:
        tuple val(meta), path(reads)
    ```

=== "Vor"

    ```groovy title="modules/fastp.nf" linenums="1" hl_lines="5"
    process FASTP {
        container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

        cpus { meta.depth > 40000000 ? 4 : 2 }
        memory 2.GB

        input:
        tuple val(meta), path(reads)
    ```

Wenn der Prozess nun aufgrund unzureichenden Speichers fehlschlägt, wird Nextflow mit mehr Speicher wiederholen:

- Erster Versuch: 1 GB (task.attempt = 1)
- Zweiter Versuch: 2 GB (task.attempt = 2)

... und so weiter, bis zum `maxRetries`-Limit.

### Fazit

Dynamische Direktiven mit Closures ermöglichen dir:

- Ressourcen basierend auf Eingabeeigenschaften zuzuweisen
- Automatische Wiederholungsstrategien mit zunehmenden Ressourcen zu implementieren
- Mehrere Faktoren zu kombinieren (Metadaten, Versuchsnummer, Prioritäten)
- Bedingte Logik für komplexe Ressourcenberechnungen zu verwenden

Dies macht deine Workflows sowohl effizienter (keine Überzuweisung) als auch robuster (automatische Wiederholung mit mehr Ressourcen).

---

## 5. Bedingte Logik und Prozesssteuerung

Zuvor haben wir `.map()` mit Scripting verwendet, um Kanaldaten zu transformieren. Jetzt werden wir bedingte Logik verwenden, um zu steuern, welche Prozesse basierend auf Daten ausgeführt werden – wesentlich für flexible Workflows, die sich an verschiedene Probentypen anpassen können.

Nextflow's [Dataflow-Operatoren](https://www.nextflow.io/docs/latest/operator.html) nehmen Closures, die zur Laufzeit ausgewertet werden, und ermöglichen bedingte Logik, um Workflow-Entscheidungen basierend auf Kanalinhalt zu treffen.

### 5.1. Routing mit `.branch()`

Nehmen wir zum Beispiel an, dass unsere Sequenzierungsproben nur mit FASTP getrimmt werden müssen, wenn es sich um menschliche Proben mit einer Abdeckung über einem bestimmten Schwellenwert handelt. Mausproben oder Proben mit niedriger Abdeckung sollten stattdessen mit Trimgalore ausgeführt werden (dies ist ein konstruiertes Beispiel, aber es veranschaulicht den Punkt).

Wir haben einen einfachen Trimgalore-Prozess in `modules/trimgalore.nf` bereitgestellt. Schau ihn dir an, wenn du möchtest, aber die Details sind für diese Übung nicht wichtig. Der Hauptpunkt ist, dass wir Proben basierend auf ihren Metadaten routen wollen.

Füge den neuen Import aus `modules/trimgalore.nf` hinzu:

=== "Nach"

    ```groovy title="main.nf" linenums="1" hl_lines="2"
    include { FASTP } from './modules/fastp.nf'
    include { TRIMGALORE } from './modules/trimgalore.nf'
    ```

=== "Vor"

    ```groovy title="main.nf" linenums="1"
    include { FASTP } from './modules/fastp.nf'
    ```

... und ändere dann deinen `main.nf`-Workflow, um Proben basierend auf ihren Metadaten zu verzweigen und sie durch den entsprechenden Trimming-Prozess zu leiten, wie folgt:

=== "Nach"

    ```groovy title="main.nf" linenums="28" hl_lines="5-12"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row -> separateMetadata(row) }

        trim_branches = ch_samples
            .branch { meta, reads ->
                fastp: meta.organism == 'human' && meta.depth >= 30000000
                trimgalore: true
            }

        ch_fastp = FASTP(trim_branches.fastp)
        ch_trimgalore = TRIMGALORE(trim_branches.trimgalore)
        GENERATE_REPORT(ch_samples)
    ```

=== "Vor"

    ```groovy title="main.nf" linenums="28" hl_lines="5"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row -> separateMetadata(row) }

        ch_fastp = FASTP(ch_samples)
        GENERATE_REPORT(ch_samples)
    ```

Führe diesen modifizierten Workflow aus:

```bash
nextflow run main.nf
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [adoring_galileo] DSL2 - revision: c9e83aaef1

    executor >  local (6)
    [1d/0747ac] process > FASTP (2)           [100%] 2 of 2 ✔
    [cc/c44caf] process > TRIMGALORE (1)      [100%] 1 of 1 ✔
    [34/bd5a9f] process > GENERATE_REPORT (1) [100%] 3 of 3 ✔
    ```

Hier haben wir kleine, aber mächtige bedingte Ausdrücke innerhalb des `.branch{}`-Operators verwendet, um Proben basierend auf ihren Metadaten zu routen. Menschliche Proben mit hoher Abdeckung durchlaufen `FASTP`, während alle anderen Proben durch `TRIMGALORE` gehen.

### 5.2. Verwendung von `.filter()` mit Wahrheitswerten

Ein weiteres leistungsstarkes Muster zur Steuerung der Workflow-Ausführung ist der `.filter()`-Operator, der eine Closure verwendet, um zu bestimmen, welche Elemente in der Pipeline weitergehen sollen. In der Filter-Closure schreibst du **boolesche Ausdrücke**, die entscheiden, welche Elemente durchgelassen werden.

Nextflow (wie viele dynamische Sprachen) hat ein Konzept von **"Wahrheitswerten"** (truthiness), das bestimmt, welche Werte in booleschen Kontexten als `true` oder `false` ausgewertet werden:

- **Truthy**: Nicht-null-Werte, nicht-leere Strings, Zahlen ungleich Null, nicht-leere Collections
- **Falsy**: `null`, leere Strings `""`, Null `0`, leere Collections `[]` oder `[:]`, `false`

Das bedeutet, `meta.id` allein (ohne explizites `!= null`) prüft, ob die ID existiert und nicht leer ist. Lass uns das verwenden, um Proben herauszufiltern, die unsere Qualitätsanforderungen nicht erfüllen.

Füge Folgendes vor der Branch-Operation hinzu:

=== "Nach"

    ```groovy title="main.nf" linenums="28" hl_lines="5-11"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row -> separateMetadata(row) }

        // Ungültige oder minderwertige Proben herausfiltern
        ch_valid_samples = ch_samples
            .filter { meta, reads ->
                meta.id && meta.organism && meta.depth >= 25000000
            }

        trim_branches = ch_valid_samples
            .branch { meta, reads ->
                fastp: meta.organism == 'human' && meta.depth >= 30000000
                trimgalore: true
            }
    ```

=== "Vor"

    ```groovy title="main.nf" linenums="28" hl_lines="5"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row -> separateMetadata(row) }

        trim_branches = ch_samples
            .branch { meta, reads ->
                fastp: meta.organism == 'human' && meta.depth >= 30000000
                trimgalore: true
            }
    ```

Führe den Workflow erneut aus:

```bash
nextflow run main.nf
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W  ~  version 25.10.2
    Launching `main.nf` [lonely_williams] DSL2 - revision: d0b3f121ec
    [94/b48eac] Submitted process > FASTP (2)
    [2c/d2b28f] Submitted process > GENERATE_REPORT (2)
    [65/2e3be4] Submitted process > GENERATE_REPORT (1)
    [94/b48eac] NOTE: Process `FASTP (2)` terminated with an error exit status (137) -- Execution is retried (1)
    [3e/0d8664] Submitted process > TRIMGALORE (1)
    [6a/9137b0] Submitted process > FASTP (1)
    [6a/9137b0] NOTE: Process `FASTP (1)` terminated with an error exit status (137) -- Execution is retried (1)
    [83/577ac0] Submitted process > GENERATE_REPORT (3)
    [a2/5117de] Re-submitted process > FASTP (1)
    [1f/a1a4ca] Re-submitted process > FASTP (2)
    ```

Da wir einen Filter gewählt haben, der einige Proben ausschließt, wurden weniger Aufgaben ausgeführt.

Der Filterausdruck `meta.id && meta.organism && meta.depth >= 25000000` kombiniert Wahrheitswerte mit expliziten Vergleichen:

- `meta.id && meta.organism` prüft, ob beide Felder existieren und nicht leer sind (mit Wahrheitswerten)
- `meta.depth >= 25000000` stellt ausreichende Sequenzierungstiefe mit einem expliziten Vergleich sicher

!!! note "Wahrheitswerte in der Praxis"

    Der Ausdruck `meta.id && meta.organism` ist prägnanter als zu schreiben:
    ```groovy
    meta.id != null && meta.id != '' && meta.organism != null && meta.organism != ''
    ```

    Das macht die Filterlogik viel sauberer und leichter zu lesen.

### Fazit

In diesem Abschnitt hast du gelernt, bedingte Logik zur Steuerung der Workflow-Ausführung zu verwenden, indem du die Closure-Schnittstellen von Nextflow-Operatoren wie `.branch{}` und `.filter{}` nutzt und Wahrheitswerte verwendest, um prägnante bedingte Ausdrücke zu schreiben.

Unsere Pipeline leitet nun Proben intelligent durch entsprechende Prozesse, aber Produktions-Workflows müssen ungültige Daten elegant behandeln. Lass uns unseren Workflow robust gegen fehlende oder null Werte machen.

---

## 6. Safe Navigation und Elvis Operatoren

Unsere `separateMetadata`-Funktion geht derzeit davon aus, dass alle CSV-Felder vorhanden und gültig sind. Aber was passiert bei unvollständigen Daten? Finden wir es heraus.

### 6.1. Das Problem: Zugriff auf Eigenschaften, die nicht existieren

Nehmen wir an, wir möchten Unterstützung für optionale Sequenzierungslauf-Informationen hinzufügen. In einigen Laboren könnten Proben ein zusätzliches Feld für die Sequenzierungslauf-ID oder Batchnummer haben, aber unsere aktuelle CSV hat diese Spalte nicht. Lass uns trotzdem versuchen, darauf zuzugreifen.

Ändere die `separateMetadata`-Funktion, um ein run_id-Feld einzuschließen:

=== "Nach"

    ```groovy title="main.nf" linenums="5" hl_lines="9"
    def separateMetadata(row) {
        def sample_meta = [
            id: row.sample_id.toLowerCase(),
            organism: row.organism,
            tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
            depth: row.sequencing_depth.toInteger(),
            quality: row.quality_score.toDouble()
        ]
        def run_id = row.run_id.toUpperCase()
    ```

=== "Vor"

    ```groovy title="main.nf" linenums="5"
    def separateMetadata(row) {
        def sample_meta = [
            id: row.sample_id.toLowerCase(),
            organism: row.organism,
            tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
            depth: row.sequencing_depth.toInteger(),
            quality: row.quality_score.toDouble()
        ]
    ```

Führe nun den Workflow aus:

```bash
nextflow run main.nf
```

??? failure "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [trusting_torvalds] DSL2 - revision: b56fbfbce2

    ERROR ~ Cannot invoke method toUpperCase() on null object

    -- Check script 'main.nf' at line: 13 or see '.nextflow.log' file for more details
    ```

Dies stürzt mit einer NullPointerException ab.

Das Problem ist, dass `row.run_id` `null` zurückgibt, weil die Spalte `run_id` nicht in unserer CSV existiert. Wenn wir versuchen, `.toUpperCase()` auf `null` aufzurufen, stürzt es ab. Hier kommt der Safe-Navigation-Operator ins Spiel.

### 6.2. Safe Navigation Operator (`?.`)

Der Safe-Navigation-Operator (`?.`) gibt `null` zurück, anstatt eine Exception zu werfen, wenn er auf einen `null`-Wert aufgerufen wird. Wenn das Objekt vor `?.` `null` ist, wird der gesamte Ausdruck als `null` ausgewertet, ohne die Methode auszuführen.

Aktualisiere die Funktion, um Safe Navigation zu verwenden:

=== "Nach"

    ```groovy title="main.nf" linenums="4" hl_lines="9"
    def separateMetadata(row) {
        def sample_meta = [
            id: row.sample_id.toLowerCase(),
            organism: row.organism,
            tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
            depth: row.sequencing_depth.toInteger(),
            quality: row.quality_score.toDouble()
        ]
        def run_id = row.run_id?.toUpperCase()
    ```

=== "Vor"

    ```groovy title="main.nf" linenums="4" hl_lines="9"
    def separateMetadata(row) {
        def sample_meta = [
            id: row.sample_id.toLowerCase(),
            organism: row.organism,
            tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
            depth: row.sequencing_depth.toInteger(),
            quality: row.quality_score.toDouble()
        ]
        def run_id = row.run_id.toUpperCase()
    ```

Führe es noch einmal aus:

```bash
nextflow run main.nf
```

??? success "Befehlsausgabe"

    ```console
    <!-- TODO: Ausgabe -->
    ```

Kein Absturz mehr! Der Workflow behandelt nun das fehlende Feld elegant. Wenn `row.run_id` `null` ist, verhindert der `?.`-Operator den Aufruf von `.toUpperCase()`, und `run_id` wird `null`, anstatt eine Exception zu verursachen.

### 6.3. Elvis Operator (`?:`) für Standardwerte

Der Elvis-Operator (`?:`) stellt Standardwerte bereit, wenn die linke Seite "falsy" ist (wie zuvor erklärt). Er ist nach Elvis Presley benannt, weil `?:` wie seine berühmten Haare und Augen aussieht, wenn man es seitwärts betrachtet!

Jetzt, da wir Safe Navigation verwenden, wird `run_id` für Proben ohne dieses Feld `null` sein. Lass uns den Elvis-Operator verwenden, um einen Standardwert bereitzustellen und ihn zu unserem `sample_meta`-Map hinzuzufügen:

=== "Nach"

    ```groovy title="main.nf" linenums="5" hl_lines="9-10"
    def separateMetadata(row) {
        def sample_meta = [
            id: row.sample_id.toLowerCase(),
            organism: row.organism,
            tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
            depth: row.sequencing_depth.toInteger(),
            quality: row.quality_score.toDouble()
        ]
        def run_id = row.run_id?.toUpperCase() ?: 'UNSPECIFIED'
        sample_meta.run = run_id
    ```

=== "Vor"

    ```groovy title="main.nf" linenums="5" hl_lines="9"
    def separateMetadata(row) {
        def sample_meta = [
            id: row.sample_id.toLowerCase(),
            organism: row.organism,
            tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
            depth: row.sequencing_depth.toInteger(),
            quality: row.quality_score.toDouble()
        ]
        def run_id = row.run_id?.toUpperCase()
    ```

Füge auch einen `view()`-Operator im Workflow hinzu, um die Ergebnisse zu sehen:

=== "Nach"

    ```groovy title="main.nf" linenums="30" hl_lines="4"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map{ row -> separateMetadata(row) }
            .view()
    ```

=== "Vor"

    ```groovy title="main.nf" linenums="30"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map{ row -> separateMetadata(row) }
    ```

und führe den Workflow aus:

```bash
nextflow run main.nf
```

??? success "Befehlsausgabe"

    ```console
    [[id:sample_001, organism:human, tissue:liver, depth:30000000, quality:38.5, run:UNSPECIFIED, sample_num:1, lane:001, read:R1, chunk:001, priority:normal], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_001_S1_L001_R1_001.fastq]
    [[id:sample_002, organism:mouse, tissue:brain, depth:25000000, quality:35.2, run:UNSPECIFIED, sample_num:2, lane:001, read:R1, chunk:001, priority:normal], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_002_S2_L001_R1_001.fastq]
    [[id:sample_003, organism:human, tissue:kidney, depth:45000000, quality:42.1, run:UNSPECIFIED, sample_num:3, lane:001, read:R1, chunk:001, priority:high], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_003_S3_L001_R1_001.fastq]
    ```

Perfekt! Jetzt haben alle Proben ein `run`-Feld entweder mit ihrer tatsächlichen Lauf-ID (in Großbuchstaben) oder dem Standardwert 'UNSPECIFIED'. Die Kombination von `?.` und `?:` bietet sowohl Sicherheit (keine Abstürze) als auch sinnvolle Standardwerte.

Entferne den `.view()`-Operator jetzt, da wir bestätigt haben, dass es funktioniert.

!!! tip "Kombination von Safe Navigation und Elvis"

    Das Muster `value?.method() ?: 'default'` ist in Produktions-Workflows üblich:

    - `value?.method()` - Ruft die Methode sicher auf, gibt `null` zurück, wenn `value` `null` ist
    - `?: 'default'` - Bietet einen Fallback, wenn das Ergebnis `null` ist

    Dieses Muster behandelt fehlende/unvollständige Daten elegant.

Verwende diese Operatoren konsequent in Funktionen, Operator-Closures (`.map{}`, `.filter{}`), Prozess-Scripts und Konfigurationsdateien. Sie verhindern Abstürze bei der Verarbeitung von realen Daten.

### Fazit

- **Safe Navigation (`?.`)**: Verhindert Abstürze bei null-Werten - gibt null zurück anstatt eine Ausnahme zu werfen
- **Elvis-Operator (`?:`)**: Bietet Standardwerte - `value ?: 'default'`
- **Kombination**: `value?.method() ?: 'default'` ist das übliche Muster

Diese Operatoren machen Workflows widerstandsfähig gegenüber unvollständigen Daten - wesentlich für die Arbeit in der realen Welt.

---

## 7. Validierung mit `error()` und `log.warn`

Manchmal musst du den Workflow sofort stoppen, wenn Eingabeparameter ungültig sind. In Nextflow kannst du eingebaute Funktionen wie `error()` und `log.warn` sowie Standard-Programmierstrukturen wie `if`-Anweisungen und boolesche Logik verwenden, um Validierungslogik zu implementieren. Fügen wir unserem Workflow Validierung hinzu.

Erstelle eine Validierungsfunktion vor deinem Workflow-Block, rufe sie aus dem Workflow auf und ändere die Channel-Erstellung, um einen Parameter für den CSV-Dateipfad zu verwenden. Wenn der Parameter fehlt oder die Datei nicht existiert, rufe `error()` auf, um die Ausführung mit einer klaren Nachricht zu stoppen.

=== "Nach"

    ```groovy title="main.nf" linenums="1" hl_lines="5-20 23-24"
    include { FASTP } from './modules/fastp.nf'
    include { TRIMGALORE } from './modules/trimgalore.nf'
    include { GENERATE_REPORT } from './modules/generate_report.nf'

    def validateInputs() {
        // Prüfe, ob der Eingabeparameter bereitgestellt wird
        if (!params.input) {
            error("CSV-Dateipfad nicht angegeben. Bitte gib --input <file.csv> an")
        }

        // Prüfe, ob die CSV-Datei existiert
        if (!file(params.input).exists()) {
            error("CSV-Datei nicht gefunden: ${params.input}")
        }
    }
    ...
    workflow {
        validateInputs()
        ch_samples = channel.fromPath(params.input)
    ```

=== "Vor"

    ```groovy title="main.nf" linenums="1"
    include { FASTP } from './modules/fastp.nf'
    include { TRIMGALORE } from './modules/trimgalore.nf'
    include { GENERATE_REPORT } from './modules/generate_report.nf'

    ...
    workflow {
        ch_samples = channel.fromPath("./data/samples.csv")
    ```

Versuche nun, ohne die CSV-Datei auszuführen:

```bash
nextflow run main.nf
```

??? failure "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [confident_coulomb] DSL2 - revision: 07059399ed

    WARN: Access to undefined parameter `input` -- Initialise it to a default value eg. `params.input = some_value`
    CSV-Dateipfad nicht angegeben. Bitte gib --input <file.csv> an
    ```

Der Workflow stoppt sofort mit einer klaren Fehlermeldung, anstatt später mysteriös zu scheitern

Führe ihn jetzt mit einer nicht existierenden Datei aus:

```bash
nextflow run main.nf --input ./data/nonexistent.csv
```

??? failure "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [cranky_gates] DSL2 - revision: 26839ae3eb

    CSV-Datei nicht gefunden: ./data/nonexistent.csv
    ```

Führe ihn schließlich mit der richtigen Datei aus:

```bash
nextflow run main.nf --input ./data/samples.csv
```

??? success "Befehlsausgabe"

    ```console
    <!-- TODO: Ausgabe -->
    ```

Diesmal läuft es erfolgreich.

Du kannst auch Validierung innerhalb der `separateMetadata`-Funktion hinzufügen. Lass uns den nicht fatalen `log.warn` verwenden, um Warnungen für Proben mit geringer Sequenzierungstiefe auszugeben, aber den Workflow trotzdem weiterlaufen lassen:

=== "Nach"

    ```groovy title="main.nf" linenums="1" hl_lines="3-6"
        def priority = sample_meta.quality > 40 ? 'high' : 'normal'

        // Validiere, ob Daten Sinn machen
        if (sample_meta.depth < 30000000) {
            log.warn "Geringe Sequenzierungstiefe für ${sample_meta.id}: ${sample_meta.depth}"
        }

        return tuple(sample_meta + file_meta + [priority: priority], fastq_path)
    }
    ```

=== "Vor"

    ```groovy title="main.nf" linenums="1"
        def priority = sample_meta.quality > 40 ? 'high' : 'normal'

        return tuple(sample_meta + file_meta + [priority: priority], fastq_path)
    }
    ```

Führe den Workflow erneut mit der ursprünglichen CSV aus:

```bash
nextflow run main.nf --input ./data/samples.csv
```

??? warning "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [awesome_goldwasser] DSL2 - revision: a31662a7c1

    executor >  local (5)
    [ce/df5eeb] process > FASTP (2)           [100%] 2 of 2 ✔
    [-        ] process > TRIMGALORE          -
    [d1/7d2b4b] process > GENERATE_REPORT (3) [100%] 3 of 3 ✔
    WARN: Geringe Sequenzierungstiefe für sample_002: 25000000
    ```

Wir sehen eine Warnung über geringe Sequenzierungstiefe für eine der Proben.

### Fazit

- **`error()`**: Stoppt den Workflow sofort mit einer klaren Nachricht
- **`log.warn`**: Gibt Warnungen aus, ohne den Workflow zu stoppen
- **Frühzeitige Validierung**: Überprüft Eingaben vor der Verarbeitung, um schnell mit hilfreichen Fehlermeldungen zu scheitern
- **Validierungsfunktionen**: Erstelle wiederverwendbare Validierungslogik, die am Workflow-Start aufgerufen werden kann

Ordnungsgemäße Validierung macht Workflows robuster und benutzerfreundlicher, indem Probleme frühzeitig mit klaren Fehlermeldungen erkannt werden.

---

## 8. Workflow-Ereignishandler

Bis jetzt haben wir Code in unseren Workflow-Skripten und Prozessdefinitionen geschrieben. Aber es gibt noch ein weiteres wichtiges Feature, das du kennen solltest: Workflow-Ereignishandler.

Ereignishandler sind Closures, die zu bestimmten Zeitpunkten im Lebenszyklus deines Workflows ausgeführt werden. Sie sind perfekt, um Logging, Benachrichtigungen oder Aufräumoperationen hinzuzufügen. Diese Handler sollten in deinem Workflow-Skript neben deiner Workflow-Definition definiert werden.

### 8.1. Der `onComplete`-Handler

Der am häufigsten verwendete Ereignishandler ist `onComplete`, der ausgeführt wird, wenn dein Workflow endet (ob er erfolgreich war oder nicht). Fügen wir einen hinzu, um unsere Pipeline-Ergebnisse zusammenzufassen.

Füge den Ereignishandler zu deiner `main.nf`-Datei hinzu, innerhalb deiner Workflow-Definition:

=== "Nach"

    ```groovy title="main.nf" linenums="66" hl_lines="5-16"
        ch_fastp = FASTP(trim_branches.fastp)
        ch_trimgalore = TRIMGALORE(trim_branches.trimgalore)
        GENERATE_REPORT(ch_samples)

        workflow.onComplete = {
            println ""
            println "Pipeline-Ausführungszusammenfassung:"
            println "=========================="
            println "Abgeschlossen um: ${workflow.complete}"
            println "Dauer           : ${workflow.duration}"
            println "Erfolg          : ${workflow.success}"
            println "workDir         : ${workflow.workDir}"
            println "Exit-Status     : ${workflow.exitStatus}"
            println ""
        }
    }
    ```

=== "Vor"

    ```groovy title="main.nf" linenums="66" hl_lines="4"
        ch_fastp = FASTP(trim_branches.fastp)
        ch_trimgalore = TRIMGALORE(trim_branches.trimgalore)
        GENERATE_REPORT(ch_samples)
    }
    ```

Diese Closure wird ausgeführt, wenn der Workflow abgeschlossen ist. Innerhalb hast du Zugriff auf das `workflow`-Objekt, das nützliche Eigenschaften über die Ausführung bereitstellt.

Führe deinen Workflow aus, und du wirst sehen, dass diese Zusammenfassung am Ende erscheint!

```bash
nextflow run main.nf --input ./data/samples.csv -ansi-log false
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W  ~  version 25.10.2
    Launching `main.nf` [marvelous_boltzmann] DSL2 - revision: a31662a7c1
    WARN: Geringe Sequenzierungstiefe für sample_002: 25000000
    [9b/d48e40] Submitted process > FASTP (2)
    [6a/73867a] Submitted process > GENERATE_REPORT (2)
    [79/ad0ac5] Submitted process > GENERATE_REPORT (1)
    [f3/bda6cb] Submitted process > FASTP (1)
    [34/d5b52f] Submitted process > GENERATE_REPORT (3)

    Pipeline-Ausführungszusammenfassung:
    ==========================
    Abgeschlossen um: 2025-10-10T12:14:24.885384+01:00
    Dauer           : 2.9s
    Erfolg          : true
    workDir         : /workspaces/training/side-quests/essential_scripting_patterns/work
    Exit-Status     : 0
    ```

Machen wir ihn nützlicher, indem wir bedingte Logik hinzufügen:

=== "Nach"

    ```groovy title="main.nf" linenums="66" hl_lines="5-22"
        ch_fastp = FASTP(trim_branches.fastp)
        ch_trimgalore = TRIMGALORE(trim_branches.trimgalore)
        GENERATE_REPORT(ch_samples)

        workflow.onComplete = {
            println ""
            println "Pipeline-Ausführungszusammenfassung:"
            println "=========================="
            println "Abgeschlossen um: ${workflow.complete}"
            println "Dauer           : ${workflow.duration}"
            println "Erfolg          : ${workflow.success}"
            println "workDir         : ${workflow.workDir}"
            println "Exit-Status     : ${workflow.exitStatus}"
            println ""

            if (workflow.success) {
                println "✅ Pipeline erfolgreich abgeschlossen!"
            } else {
                println "❌ Pipeline fehlgeschlagen!"
                println "Fehler: ${workflow.errorMessage}"
            }
        }
    }
    ```

=== "Vor"

    ```groovy title="main.nf" linenums="66" hl_lines="5-16"
        ch_fastp = FASTP(trim_branches.fastp)
        ch_trimgalore = TRIMGALORE(trim_branches.trimgalore)
        GENERATE_REPORT(ch_samples)

        workflow.onComplete = {
            println ""
            println "Pipeline-Ausführungszusammenfassung:"
            println "=========================="
            println "Abgeschlossen um: ${workflow.complete}"
            println "Dauer           : ${workflow.duration}"
            println "Erfolg          : ${workflow.success}"
            println "workDir         : ${workflow.workDir}"
            println "Exit-Status     : ${workflow.exitStatus}"
            println ""
        }
    }
    ```

Jetzt erhalten wir eine noch informativere Zusammenfassung, einschließlich einer Erfolgs-/Fehlermeldung und des Ausgabeverzeichnisses, falls angegeben:

<!-- TODO: Befehl zum Ausführen hinzufügen -->

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W  ~  version 25.10.2
    Launching `main.nf` [boring_linnaeus] DSL2 - revision: a31662a7c1
    WARN: Geringe Sequenzierungstiefe für sample_002: 25000000
    [e5/242efc] Submitted process > FASTP (2)
    [3b/74047c] Submitted process > GENERATE_REPORT (3)
    [8a/7a57e6] Submitted process > GENERATE_REPORT (1)
    [a8/b1a31f] Submitted process > GENERATE_REPORT (2)
    [40/648429] Submitted process > FASTP (1)

    Pipeline-Ausführungszusammenfassung:
    ==========================
    Abgeschlossen um: 2025-10-10T12:16:00.522569+01:00
    Dauer           : 3.6s
    Erfolg          : true
    workDir         : /workspaces/training/side-quests/essential_scripting_patterns/work
    Exit-Status     : 0

    ✅ Pipeline erfolgreich abgeschlossen!
    ```

Du kannst die Zusammenfassung auch in eine Datei schreiben, indem du Dateioperationen verwendest:

```groovy title="main.nf - Schreiben der Zusammenfassung in eine Datei"
workflow {
    // ... dein Workflow-Code ...

    workflow.onComplete = {
        def summary = """
        Pipeline-Ausführungszusammenfassung
        ===========================
        Abgeschlossen: ${workflow.complete}
        Dauer        : ${workflow.duration}
        Erfolg       : ${workflow.success}
        Befehl       : ${workflow.commandLine}
        """

        println summary

        // In eine Log-Datei schreiben
        def log_file = file("${workflow.launchDir}/pipeline_summary.txt")
        log_file.text = summary
    }
}
```

### 8.2. Der `onError`-Handler

Neben `onComplete` gibt es einen weiteren Ereignishandler, den du verwenden kannst: `onError`, der nur ausgeführt wird, wenn der Workflow fehlschlägt:

```groovy title="main.nf - onError-Handler"
workflow {
    // ... dein Workflow-Code ...

    workflow.onError = {
        println "="* 50
        println "Pipeline-Ausführung fehlgeschlagen!"
        println "Fehlermeldung: ${workflow.errorMessage}"
        println "="* 50

        // Detailliertes Fehlerprotokoll schreiben
        def error_file = file("${workflow.launchDir}/error.log")
        error_file.text = """
        Workflow-Fehlerbericht
        =====================
        Zeit: ${new Date()}
        Fehler: ${workflow.errorMessage}
        Fehlerbericht: ${workflow.errorReport ?: 'Kein detaillierter Bericht verfügbar'}
        """

        println "Fehlerdetails geschrieben nach: ${error_file}"
    }
}
```

Du kannst mehrere Handler zusammen in deinem Workflow-Skript verwenden:

```groovy title="main.nf - Kombinierte Handler"
workflow {
    // ... dein Workflow-Code ...

    workflow.onError = {
        println "Workflow fehlgeschlagen: ${workflow.errorMessage}"
    }

    workflow.onComplete = {
        def duration_mins = workflow.duration.toMinutes().round(2)
        def status = workflow.success ? "ERFOLGREICH ✅" : "FEHLGESCHLAGEN ❌"

        println """
        Pipeline beendet: ${status}
        Dauer: ${duration_mins} Minuten
        """
    }
}
```

### Fazit

In diesem Abschnitt hast du Folgendes gelernt:

- **Ereignishandler-Closures**: Closures in deinem Workflow-Skript, die zu verschiedenen Lebenszykluspunkten ausgeführt werden
- **`onComplete`-Handler**: Für Ausführungszusammenfassungen und Ergebnisberichte
- **`onError`-Handler**: Für Fehlerbehandlung und Protokollierung von Fehlern
- **Workflow-Objekteigenschaften**: Zugriff auf `workflow.success`, `workflow.duration`, `workflow.errorMessage`, usw.

Ereignishandler zeigen, wie du die volle Macht der Nextflow-Sprache innerhalb deiner Workflow-Skripte nutzen kannst, um ausgeklügelte Logging- und Benachrichtigungsfunktionen hinzuzufügen.

---

## Zusammenfassung

Herzlichen Glückwunsch, du hast es geschafft!

In diesem Side Quest hast du eine umfassende Probenverarbeitungs-Pipeline erstellt, die sich von der grundlegenden Metadatenverarbeitung zu einem ausgeklügelten, produktionsreifen Workflow entwickelt hat.
Jeder Abschnitt baute auf dem vorherigen auf und zeigte, wie Programmierstrukturen einfache Workflows in leistungsfähige Datenverarbeitungssysteme verwandeln, mit folgenden Vorteilen:

- **Klarerer Code**: Das Verständnis von Datenfluss vs. Scripting hilft dir, besser organisierte Workflows zu schreiben
- **Robuste Handhabung**: Safe-Navigation und Elvis-Operatoren machen Workflows widerstandsfähig gegen fehlende Daten
- **Flexible Verarbeitung**: Bedingte Logik lässt deine Workflows verschiedene Probentypen angemessen verarbeiten
- **Adaptive Ressourcen**: Dynamische Direktiven optimieren die Ressourcennutzung basierend auf Eingabeeigenschaften

Diese Progression spiegelt die reale Evolution von Bioinformatik-Pipelines wider, von Forschungsprototypen, die einige Proben verarbeiten, bis hin zu Produktionssystemen, die Tausende von Proben über Labore und Einrichtungen hinweg verarbeiten.
Jede Herausforderung, die du gelöst hast, und jedes Muster, das du gelernt hast, spiegelt tatsächliche Probleme wider, mit denen Entwickler bei der Skalierung von Nextflow-Workflows konfrontiert sind.

Die Anwendung dieser Muster in deiner eigenen Arbeit wird es dir ermöglichen, robuste, produktionsreife Workflows zu erstellen.

### Schlüsselmuster

1.  **Datenfluss vs. Scripting:** Du hast gelernt, zwischen Datenflussoperationen (Channel-Orchestrierung) und Scripting (Code, der Daten manipuliert) zu unterscheiden, einschließlich der entscheidenden Unterschiede zwischen Operationen auf verschiedenen Typen wie `collect` auf Channel vs. Liste.

    - Datenfluss: Channel-Orchestrierung

    ```groovy
    channel.fromPath('*.fastq').splitCsv(header: true)
    ```

    - Scripting: Datenverarbeitung auf Collections

    ```groovy
    sample_data.collect { it.toUpperCase() }
    ```

2.  **Fortgeschrittene String-Verarbeitung**: Du hast reguläre Ausdrücke zum Parsen von Dateinamen, dynamische Skripterstellung in Prozessen und Variableninterpolation (Nextflow vs. Bash vs. Shell) gemeistert.

    - Mustererkennung

    ```groovy
    filename =~ ~/^(\w+)_(\w+)_(\d+)\.fastq$/
    ```

    - Funktion mit bedingter Rückgabe

    ```groovy
    def parseSample(filename) {
        def matcher = filename =~ pattern
        return matcher ? [valid: true, data: matcher[0]] : [valid: false]
    }
    ```

    - Dateisammlung zu Befehlsargumenten (im Prozess-Script-Block)

    ```groovy
    script:
    def file_args = input_files.collect { file -> "--input ${file}" }.join(' ')
    """
    analysis_tool ${file_args} --output results.txt
    """
    ```

3.  **Wiederverwendbare Funktionen erstellen**: Du hast gelernt, komplexe Logik in benannte Funktionen zu extrahieren, die von Channel-Operatoren aufgerufen werden können, wodurch Workflows lesbarer und wartbarer werden.

    - Eine benannte Funktion definieren

    ```groovy
    def separateMetadata(row) {
        def sample_meta = [ /* Code aus Platzgründen gekürzt */ ]
        def fastq_path = file(row.file_path)
        def m = (fastq_path.name =~ /^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$/)
        def file_meta = m ? [ /* Code aus Platzgründen gekürzt */ ] : [:]
        def priority = sample_meta.quality > 40 ? 'high' : 'normal'

        return tuple(sample_meta + file_meta + [priority: priority], fastq_path)
    }
    ```

    - Die benannte Funktion in einem Workflow aufrufen

    ```groovy
    workflow {
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map{ row -> separateMetadata(row) }

        ch_fastp = FASTP(ch_samples)
    }
    ```

4.  **Dynamische Ressourcendirektiven mit Closures**: Du hast die Verwendung von Closures in Prozessdirektiven für adaptive Ressourcenzuweisung basierend auf Eingabeeigenschaften untersucht.

    - Benannte Closures und Komposition

    ```groovy
    def enrichData = normalizeId >> addQualityCategory >> addFlags
    def processor = generalFunction.curry(fixedParam)
    ```

    - Closures mit Bereichszugriff

    ```groovy
    def collectStats = { data -> stats.count++; return data }
    ```

5.  **Bedingte Logik und Prozesssteuerung**: Du hast intelligentes Routing mit den Operatoren `.branch()` und `.filter()` hinzugefügt und dabei Wahrheitswerte für prägnante bedingte Ausdrücke genutzt.

    - Verwende `.branch()`, um Daten durch verschiedene Workflow-Zweige zu leiten

    ```groovy
    trim_branches = ch_samples
    .branch { meta, reads ->
        fastp: meta.organism == 'human' && meta.depth >= 30000000
        trimgalore: true
    }

    ch_fastp = FASTP(trim_branches.fastp)
    ch_trimgalore = TRIMGALORE(trim_branches.trimgalore)
    ```

    - Boolesche Auswertung mit Groovy Truth

    ```groovy
    if (sample.files) println "Hat Dateien"
    ```

    - Verwende `filter()`, um Daten mit 'Wahrheitswerten' zu filtern

    ```groovy
    ch_valid_samples = ch_samples
        .filter { meta, reads ->
            meta.id && meta.organism && meta.depth >= 25000000
        }
    ```

6.  **Safe Navigation und Elvis-Operatoren**: Du hast die Pipeline robust gegen fehlende Daten gemacht, indem du `?.` für nullsicheren Eigenschaftszugriff und `?:` für die Bereitstellung von Standardwerten verwendet hast.

    ```groovy
    def id = data?.sample?.id ?: 'unbekannt'
    ```

7.  **Validierung mit error() und log.warn**: Du hast gelernt, Eingaben frühzeitig zu validieren und schnell mit klaren Fehlermeldungen zu scheitern.

    ```groovy
    try {
        def errors = validateSample(sample)
        if (errors) throw new RuntimeException("Ungültig: ${errors.join(', ')}")
    } catch (Exception e) {
        println "Fehler: ${e.message}"
    }
    ```

8.  **Workflow-Ereignishandler**: Du hast gelernt, Workflow-Ereignishandler (`onComplete` und `onError`) für Logging, Benachrichtigungen und Lebenszyklus-Management zu verwenden.

    - Verwendung von `onComplete` für Logging und Benachrichtigung

    ```groovy
    workflow.onComplete = {
        println "Erfolg     : ${workflow.success}"
        println "Exit-Status: ${workflow.exitStatus}"

        if (workflow.success) {
            println "✅ Pipeline erfolgreich abgeschlossen!"
        } else {
            println "❌ Pipeline fehlgeschlagen!"
            println "Fehler: ${workflow.errorMessage}"
        }
    }
    ```

    - Verwendung von `onError`, um speziell im Fehlerfall Maßnahmen zu ergreifen

    ```groovy
    workflow.onError = {
        // Detailliertes Fehlerprotokoll schreiben
        def error_file = file("${workflow.launchDir}/error.log")
        error_file.text = """
        Zeit: ${new Date()}
        Fehler: ${workflow.errorMessage}
        Fehlerbericht: ${workflow.errorReport ?: 'Kein detaillierter Bericht verfügbar'}
        """

        println "Fehlerdetails geschrieben nach: ${error_file}"
    }
    ```

### Zusätzliche Ressourcen

- [Nextflow-Sprachreferenz](https://nextflow.io/docs/latest/reference/syntax.html)
- [Nextflow-Operatoren](https://www.nextflow.io/docs/latest/operator.html)
- [Nextflow-Skriptsyntax](https://www.nextflow.io/docs/latest/script.html)
- [Nextflow-Standardbibliothek](https://nextflow.io/docs/latest/reference/stdlib.html)

Schau dir diese Ressourcen an, wenn du fortgeschrittenere Funktionen erkunden musst.

Du wirst davon profitieren, deine Fähigkeiten zu üben und zu erweitern, um:

- Sauberere Workflows mit richtiger Trennung zwischen Datenfluss und Scripting zu schreiben
- Variable-Interpolation zu beherrschen, um häufige Fallstricke mit Nextflow-, Bash- und Shell-Variablen zu vermeiden
- Dynamische Ressourcendirektiven für effiziente, adaptive Workflows zu verwenden
- Dateisammlungen in korrekt formatierte Befehlszeilenargumente zu transformieren
- Verschiedene Dateibenennungskonventionen und Eingabeformate elegant zu behandeln, indem du Regex und String-Verarbeitung verwendest
- Wiederverwendbaren, wartbaren Code mit fortgeschrittenen Closure-Mustern und funktionaler Programmierung zu erstellen
- Komplexe Datensätze mit Collection-Operationen zu verarbeiten und zu organisieren
- Validierung, Fehlerbehandlung und Logging hinzuzufügen, um deine Workflows produktionsreif zu machen
- Workflow-Lebenszyklusmanagement mit Ereignishandlern zu implementieren

---

## Was kommt als nächstes?

Kehre zum [Menü der Side Quests](./index.md) zurück oder klicke auf die Schaltfläche unten rechts auf der Seite, um zum nächsten Thema in der Liste zu gelangen.
