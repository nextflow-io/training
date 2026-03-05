# Grundlegende Nextflow-Scripting-Muster

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Nextflow ist eine Programmiersprache, die auf der Java Virtual Machine läuft. Obwohl Nextflow auf [Groovy](http://groovy-lang.org/) aufbaut und viel von dessen Syntax teilt, ist Nextflow mehr als nur „Groovy mit Erweiterungen" – es ist eine eigenständige Sprache mit einer vollständig spezifizierten [Syntax](https://nextflow.io/docs/latest/reference/syntax.html) und [Standardbibliothek](https://nextflow.io/docs/latest/reference/stdlib.html).

Du kannst viel Nextflow schreiben, ohne über die grundlegende Syntax für Variablen, Maps und Listen hihinauszugehen. Die meisten Nextflow-Tutorials konzentrieren sich auf Workflow-Orchestrierung (Kanäle, Prozesse und Datenfluss), und damit kommst du überraschend weit.

Wenn du jedoch Daten manipulieren, komplexe Dateinamen parsen, bedingte Logik implementieren oder robuste Produktions-Workflows erstellen musst, hilft es, zwei verschiedene Aspekte deines Codes zu betrachten: **Datenfluss** (Kanäle, Operatoren, Prozesse und Workflows) und **Scripting** (der Code innerhalb von Closures, Funktionen und Prozess-Skripten). Obwohl diese Unterscheidung etwas willkürlich ist – es ist alles Nextflow-Code – bietet sie ein nützliches mentales Modell, um zu verstehen, wann du deine Pipeline orchestrierst und wann du Daten manipulierst. Die Beherrschung beider Aspekte verbessert deine Fähigkeit, klare, wartbare Workflows zu schreiben, dramatisch.

### Lernziele

Dieser Side Quest nimmt dich mit auf eine praktische Reise von grundlegenden Konzepten zu produktionsreifen Mustern.
Wir werden einen einfachen CSV-lesenden Workflow in eine anspruchsvolle Bioinformatik-Pipeline verwandeln und ihn Schritt für Schritt durch realistische Herausforderungen weiterentwickeln:

- **Grenzen verstehen:** Unterscheide zwischen Datenfluss-Operationen und Scripting und verstehe, wie sie zusammenarbeiten
- **Datenmanipulation:** Extrahiere, transformiere und unterteile Maps und Collections mit leistungsstarken Operatoren
- **String-Verarbeitung:** Parse komplexe Dateibenennungsschemata mit Regex-Mustern und beherrsche Variableninterpolation
- **Wiederverwendbare Funktionen:** Extrahiere komplexe Logik in benannte Funktionen für sauberere, wartbarere Workflows
- **Dynamische Logik:** Erstelle Prozesse, die sich an verschiedene Eingabetypen anpassen, und verwende Closures für dynamische Ressourcenzuweisung
- **Bedingtes Routing:** Route Proben intelligent durch verschiedene Prozesse basierend auf ihren Metadaten-Eigenschaften
- **Sichere Operationen:** Behandle fehlende Daten elegant mit null-sicheren Operatoren und validiere Eingaben mit klaren Fehlermeldungen
- **Konfigurationsbasierte Handler:** Verwende Workflow-Event-Handler für Logging, Benachrichtigungen und Lifecycle-Management

### Voraussetzungen

Bevor du diesen Side Quest beginnst, solltest du:

- Das [Hello Nextflow](../hello_nextflow/README.md)-Tutorial oder einen gleichwertigen Einsteigerkurs abgeschlossen haben
- Dich mit grundlegenden Nextflow-Konzepten und -Mechanismen wohlfühlen (Prozesse, Kanäle, Operatoren, Arbeiten mit Dateien, Metadaten)
- Grundlegende Vertrautheit mit gängigen Programmierkonzepten haben (Variablen, Maps, Listen)

Dieses Tutorial erklärt Programmierkonzepte, wenn wir auf sie stoßen, du brauchst also keine umfangreiche Programmiererfahrung.
Wir beginnen mit grundlegenden Konzepten und bauen zu fortgeschrittenen Mustern auf.

---

## 0. Erste Schritte

#### Öffne den Training-Codespace

Falls du es noch nicht getan hast, stelle sicher, dass du die Trainingsumgebung wie in [Umgebung einrichten](../envsetup/index.md) beschrieben öffnest.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Wechsle in das Projektverzeichnis

Lass uns in das Verzeichnis wechseln, in dem sich die Dateien für dieses Tutorial befinden.

```bash
cd side-quests/essential_scripting_patterns
```

#### Überprüfe die Materialien

Du findest eine Haupt-Workflow-Datei und ein `data`-Verzeichnis mit Beispieldaten.

```console title="Directory contents"
.
├── collect.nf
├── data
│   ├── samples.csv
│   └── sequences
│       ├── SAMPLE_001_S1_L001_R1_001.fastq
│       ├── SAMPLE_002_S2_L001_R1_001.fastq
│       └── SAMPLE_003_S3_L001_R1_001.fastq
├── main.nf
├── modules
│   ├── fastp.nf
│   ├── generate_report.nf
│   └── trimgalore.nf
└── nextflow.config
```

Unsere Beispiel-CSV enthält Informationen über biologische Proben, die je nach ihren Eigenschaften unterschiedlich verarbeitet werden müssen:

```console title="samples.csv"
sample_id,organism,tissue_type,sequencing_depth,file_path,quality_score
SAMPLE_001,human,liver,30000000,data/sequences/SAMPLE_001_S1_L001_R1_001.fastq,38.5
SAMPLE_002,mouse,brain,25000000,data/sequences/SAMPLE_002_S2_L001_R1_001.fastq,35.2
SAMPLE_003,human,kidney,45000000,data/sequences/SAMPLE_003_S3_L001_R1_001.fastq,42.1
```

Wir werden diesen realistischen Datensatz verwenden, um praktische Programmiertechniken zu erkunden, die dir in echten Bioinformatik-Workflows begegnen werden.

<!-- TODO: Can we make this more domain-agnostic? -->

<!-- TODO: add an assignment statement? #### Review the assignment -->

#### Bereitschafts-Checkliste

Denkst du, du bist bereit loszulegen?

- [ ] Ich verstehe das Ziel dieses Kurses und seine Voraussetzungen
- [ ] Mein Codespace läuft
- [ ] Ich habe mein Arbeitsverzeichnis entsprechend gesetzt
<!-- - [ ] I understand the assignment -->

Wenn du alle Kästchen abhaken kannst, kann es losgehen.

---

## 1. Datenfluss vs. Scripting: Die Grenzen verstehen

### 1.1. Identifizieren, was was ist

Beim Schreiben von Nextflow-Workflows ist es wichtig, zwischen **Datenfluss** (wie Daten durch Kanäle und Prozesse fließen) und **Scripting** (der Code, der Daten manipuliert und Entscheidungen trifft) zu unterscheiden. Lass uns einen Workflow erstellen, der zeigt, wie sie zusammenarbeiten.

#### 1.1.1. Grundlegender Nextflow-Workflow

Beginne mit einem einfachen Workflow, der nur die CSV-Datei liest (wir haben das bereits für dich in `main.nf` gemacht):

```groovy title="main.nf" linenums="1"
workflow {
    ch_samples = channel.fromPath("./data/samples.csv")
        .splitCsv(header: true)
        .view()
}
```

Der `workflow`-Block definiert unsere Pipeline-Struktur, während `channel.fromPath()` einen Kanal aus einem Dateipfad erstellt. Der `.splitCsv()`-Operator verarbeitet die CSV-Datei und konvertiert jede Zeile in eine Map-Datenstruktur.

Führe diesen Workflow aus, um die rohen CSV-Daten zu sehen:

```bash
nextflow run main.nf
```

??? success "Befehlsausgabe"

    ```console
    Launching `main.nf` [marvelous_tuckerman] DSL2 - revision: 6113e05c17

    [sample_id:SAMPLE_001, organism:human, tissue_type:liver, sequencing_depth:30000000, file_path:data/sequences/SAMPLE_001_S1_L001_R1_001.fastq, quality_score:38.5]
    [sample_id:SAMPLE_002, organism:mouse, tissue_type:brain, sequencing_depth:25000000, file_path:data/sequences/SAMPLE_002_S2_L001_R1_001.fastq, quality_score:35.2]
    [sample_id:SAMPLE_003, organism:human, tissue_type:kidney, sequencing_depth:45000000, file_path:data/sequences/SAMPLE_003_S3_L001_R1_001.fastq, quality_score:42.1]
    ```

#### 1.1.2. Hinzufügen des Map-Operators

Jetzt werden wir Scripting hinzufügen, um die Daten zu transformieren, indem wir den `.map()`-Operator verwenden, den du wahrscheinlich schon kennst. Dieser Operator nimmt eine 'Closure', in der wir Code schreiben können, um jedes Element zu transformieren.

!!! note "Hinweis"

    Eine **Closure** ist ein Codeblock, der herumgereicht und später ausgeführt werden kann. Stell sie dir wie eine Funktion vor, die du inline definierst. Closures werden mit geschweiften Klammern `{ }` geschrieben und können Parameter annehmen. Sie sind grundlegend dafür, wie Nextflow-Operatoren funktionieren, und wenn du schon eine Weile Nextflow schreibst, hast du sie vielleicht bereits verwendet, ohne es zu merken!

So sieht diese Map-Operation aus:

=== "Danach"

    ```groovy title="main.nf" linenums="2" hl_lines="3-6"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                return row
            }
            .view()
    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="2" hl_lines="3"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .view()
    ```

Das ist unsere erste **Closure** – eine anonyme Funktion, die du als Argument übergeben kannst (ähnlich wie Lambdas in Python oder Arrow-Funktionen in JavaScript). Closures sind essenziell für die Arbeit mit Nextflow-Operatoren.

Die Closure `{ row -> return row }` nimmt einen Parameter `row` (könnte jeder Name sein: `item`, `sample`, etc.).

Wenn der `.map()`-Operator jedes Kanal-Element verarbeitet, übergibt er dieses Element an deine Closure. Hier enthält `row` jeweils eine CSV-Zeile.

Wende diese Änderung an und führe den Workflow aus:

```bash
nextflow run main.nf
```

Du wirst die gleiche Ausgabe wie zuvor sehen, weil wir einfach die Eingabe unverändert zurückgeben. Das bestätigt, dass der Map-Operator korrekt funktioniert. Jetzt lass uns beginnen, die Daten zu transformieren.

#### 1.1.3. Erstellen einer Map-Datenstruktur

Jetzt werden wir **Scripting**-Logik innerhalb unserer Closure schreiben, um jede Datenzeile zu transformieren. Hier verarbeiten wir einzelne Datenelemente, anstatt den Datenfluss zu orchestrieren.

=== "Danach"

    ```groovy title="main.nf" linenums="2" hl_lines="4-12"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                // Scripting für Datentransformation
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                return sample_meta
            }
            .view()
    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="2" hl_lines="4"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                return row
            }
            .view()
    ```

Die `sample_meta`-Map ist eine Schlüssel-Wert-Datenstruktur (wie Dictionaries in Python, Objekte in JavaScript oder Hashes in Ruby), die zusammengehörige Informationen speichert: Proben-ID, Organismus, Gewebetyp, Sequenzierungstiefe und Qualitätswert.

Wir verwenden String-Manipulationsmethoden wie `.toLowerCase()` und `.replaceAll()`, um unsere Daten zu bereinigen, und Typkonvertierungsmethoden wie `.toInteger()` und `.toDouble()`, um String-Daten aus der CSV in die entsprechenden numerischen Typen zu konvertieren.

Wende diese Änderung an und führe den Workflow aus:

```bash
nextflow run main.nf
```

??? success "Befehlsausgabe"

    ```console
    [id:sample_001, organism:human, tissue:liver, depth:30000000, quality:38.5]
    [id:sample_002, organism:mouse, tissue:brain, depth:25000000, quality:35.2]
    [id:sample_003, organism:human, tissue:kidney, depth:45000000, quality:42.1]
    ```

#### 1.1.4. Hinzufügen bedingter Logik

Jetzt fügen wir mehr Scripting hinzu – diesmal mit einem ternären Operator, um Entscheidungen basierend auf Datenwerten zu treffen.

Nimm folgende Änderung vor:

=== "Danach"

    ```groovy title="main.nf" linenums="2" hl_lines="11-12"
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
                def priority = sample_meta.quality > 40 ? 'high' : 'normal'
                return sample_meta + [priority: priority]
            }
            .view()
    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="2" hl_lines="11"
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
                return sample_meta
            }
            .view()
    ```

Der ternäre Operator ist eine Kurzform für eine if/else-Anweisung, die dem Muster `bedingung ? wert_wenn_wahr : wert_wenn_falsch` folgt. Diese Zeile bedeutet: „Wenn die Qualität größer als 40 ist, verwende 'high', andernfalls verwende 'normal'". Sein Cousin, der **Elvis-Operator** (`?:`), liefert Standardwerte, wenn etwas null oder leer ist – wir werden dieses Muster später in diesem Tutorial erkunden.

Der Map-Additions-Operator `+` erstellt eine **neue Map**, anstatt die bestehende zu modifizieren. Diese Zeile erstellt eine neue Map, die alle Schlüssel-Wert-Paare aus `sample_meta` plus den neuen `priority`-Schlüssel enthält.

!!! Note "Hinweis"

    Modifiziere niemals Maps, die in Closures übergeben werden – erstelle immer neue mit `+` (zum Beispiel). In Nextflow fließen dieselben Daten oft gleichzeitig durch mehrere Operationen. Das Modifizieren einer Map an Ort und Stelle kann unvorhersehbare Nebeneffekte verursachen, wenn andere Operationen auf dasselbe Objekt verweisen. Das Erstellen neuer Maps stellt sicher, dass jede Operation ihre eigene saubere Kopie hat.

Führe den modifizierten Workflow aus:

```bash
nextflow run main.nf
```

??? success "Befehlsausgabe"

    ```console
    [id:sample_001, organism:human, tissue:liver, depth:30000000, quality:38.5, priority:normal]
    [id:sample_002, organism:mouse, tissue:brain, depth:25000000, quality:35.2, priority:normal]
    [id:sample_003, organism:human, tissue:kidney, depth:45000000, quality:42.1, priority:high]
    ```

Wir haben erfolgreich bedingte Logik hinzugefügt, um unsere Metadaten mit einer Prioritätsstufe basierend auf Qualitätswerten anzureichern.

#### 1.1.5. Unterteilen von Maps mit `.subMap()`

Während der `+`-Operator Schlüssel zu einer Map hinzufügt, musst du manchmal das Gegenteil tun – nur bestimmte Schlüssel extrahieren. Die `.subMap()`-Methode ist dafür perfekt.

Lass uns eine Zeile hinzufügen, um eine vereinfachte Version unserer Metadaten zu erstellen, die nur Identifikationsfelder enthält:

=== "Danach"

    ```groovy title="main.nf" linenums="2" hl_lines="12-15"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                // Scripting für Datentransformation
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                def id_only = sample_meta.subMap(['id', 'organism', 'tissue'])
                println "ID fields only: ${id_only}"

                def priority = sample_meta.quality > 40 ? 'high' : 'normal'
                return sample_meta + [priority: priority]
            }
            .view()
    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="2" hl_lines="12"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                // Scripting für Datentransformation
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                def priority = sample_meta.quality > 40 ? 'high' : 'normal'
                return sample_meta + [priority: priority]
            }
            .view()
    ```

Führe den modifizierten Workflow aus:

```bash
nextflow run main.nf
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [peaceful_cori] DSL2 - revision: 4cc4a8340f

    ID fields only: [id:sample_001, organism:human, tissue:liver]
    ID fields only: [id:sample_002, organism:mouse, tissue:brain]
    ID fields only: [id:sample_003, organism:human, tissue:kidney]
    [id:sample_001, organism:human, tissue:liver, depth:30000000, quality:38.5, priority:normal]
    [id:sample_002, organism:mouse, tissue:brain, depth:25000000, quality:35.2, priority:normal]
    [id:sample_003, organism:human, tissue:kidney, depth:45000000, quality:42.1, priority:high]
    ```

Dies zeigt sowohl die vollständigen Metadaten, die von der `view()`-Operation angezeigt werden, als auch die extrahierte Teilmenge, die wir mit `println` ausgegeben haben.

Die `.subMap()`-Methode nimmt eine Liste von Schlüsseln und gibt eine neue Map zurück, die nur diese Schlüssel enthält. Wenn ein Schlüssel in der ursprünglichen Map nicht existiert, wird er einfach nicht im Ergebnis enthalten.

Dies ist besonders nützlich, wenn du verschiedene Metadaten-Versionen für verschiedene Prozesse erstellen musst – einige benötigen möglicherweise vollständige Metadaten, während andere nur minimale Identifikationsfelder benötigen.

Entferne jetzt diese println-Anweisungen, um deinen Workflow in seinen vorherigen Zustand zurückzuversetzen, da wir sie nicht mehr benötigen.

!!! tip "Zusammenfassung der Map-Operationen"

    - **Schlüssel hinzufügen**: `map1 + [new_key: value]` - Erstellt neue Map mit zusätzlichen Schlüsseln
    - **Schlüssel extrahieren**: `map1.subMap(['key1', 'key2'])` - Erstellt neue Map mit nur angegebenen Schlüsseln
    - **Beide Operationen erstellen neue Maps** - Ursprüngliche Maps bleiben unverändert

#### 1.1.6. Kombinieren von Maps und Zurückgeben von Ergebnissen

Bisher haben wir nur das zurückgegeben, was die Nextflow-Community die 'Meta-Map' nennt, und wir haben die Dateien ignoriert, auf die sich diese Metadaten beziehen. Aber wenn du Nextflow-Workflows schreibst, möchtest du wahrscheinlich etwas mit diesen Dateien machen.

Lass uns eine Kanalstruktur ausgeben, die aus einem Tupel mit 2 Elementen besteht: der angereicherten Metadaten-Map und dem entsprechenden Dateipfad. Dies ist ein gängiges Muster in Nextflow zum Übergeben von Daten an Prozesse.

=== "Danach"

    ```groovy title="main.nf" linenums="2" hl_lines="12"
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
                def priority = sample_meta.quality > 40 ? 'high' : 'normal'
                return tuple( sample_meta + [priority: priority], file(row.file_path) )
            }
            .view()
    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="2" hl_lines="12"
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
                def priority = sample_meta.quality > 40 ? 'high' : 'normal'
                return sample_meta + [priority: priority]
            }
            .view()
    ```

Wende diese Änderung an und führe den Workflow aus:

```bash
nextflow run main.nf
```

??? success "Befehlsausgabe"

    ```console
    [[id:sample_001, organism:human, tissue:liver, depth:30000000, quality:38.5, priority:normal], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_001_S1_L001_R1_001.fastq]
    [[id:sample_002, organism:mouse, tissue:brain, depth:25000000, quality:35.2, priority:normal], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_002_S2_L001_R1_001.fastq]
    [[id:sample_003, organism:human, tissue:kidney, depth:45000000, quality:42.1, priority:high], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_003_S3_L001_R1_001.fastq]
    ```

Diese `[meta, file]`-Tupel-Struktur ist ein gängiges Muster in Nextflow zum Übergeben von Metadaten und zugehörigen Dateien an Prozesse.

!!! note "Hinweis"

    **Maps und Metadaten**: Maps sind grundlegend für die Arbeit mit Metadaten in Nextflow. Für eine detailliertere Erklärung der Arbeit mit Metadaten-Maps siehe den Side Quest [Arbeiten mit Metadaten](./metadata.md).

Unser Workflow demonstriert das Kernmuster: **Datenfluss-Operationen** (`workflow`, `channel.fromPath()`, `.splitCsv()`, `.map()`, `.view()`) orchestrieren, wie Daten durch die Pipeline fließen, während **Scripting** (Maps `[key: value]`, String-Methoden, Typkonvertierungen, ternäre Operatoren) innerhalb der `.map()`-Closure die Transformation einzelner Datenelemente übernimmt.

### 1.2. Verschiedene Typen verstehen: Channel vs. List

Bisher gut, wir können zwischen Datenfluss-Operationen und Scripting unterscheiden. Aber was ist, wenn derselbe Methodenname in beiden Kontexten existiert?

Ein perfektes Beispiel ist die `collect`-Methode, die sowohl für Kanal-Typen als auch für List-Typen in der Nextflow-Standardbibliothek existiert. Die `collect()`-Methode auf einer List transformiert jedes Element, während der `collect()`-Operator auf einem Kanal alle Kanal-Emissionen in einen Einzelelement-Kanal sammelt.

Lass uns das mit einigen Beispieldaten demonstrieren, indem wir uns zunächst daran erinnern, was der Kanal-`collect()`-Operator macht. Schau dir `collect.nf` an:

```groovy title="collect.nf" linenums="1"
def sample_ids = ['sample_001', 'sample_002', 'sample_003']

// channel.collect() - gruppiert mehrere Kanal-Emissionen in eine
ch_input = channel.fromList(sample_ids)
ch_input.view { sample -> "Individual channel item: ${sample}" }
ch_collected = ch_input.collect()
ch_collected.view { list -> "channel.collect() result: ${list} (${list.size()} items grouped into 1)" }
```

Schritte:

- Definiere eine Liste von Proben-IDs
- Erstelle einen Kanal mit `fromList()`, der jede Proben-ID separat ausgibt
- Gib jedes Element mit `view()` aus, während es durchfließt
- Sammle alle Elemente in eine einzelne Liste mit dem `collect()`-Operator des Kanals
- Gib das gesammelte Ergebnis (einzelnes Element mit allen Proben-IDs) mit einem zweiten `view()` aus

Wir haben die Struktur des Kanals geändert, aber wir haben die Daten selbst nicht geändert.

Führe den Workflow aus, um dies zu bestätigen:

```bash
nextflow run collect.nf
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `collect.nf` [loving_mendel] DSL2 - revision: e8d054a46e

    Individual channel item: sample_001
    Individual channel item: sample_002
    Individual channel item: sample_003
    channel.collect() result: [sample_001, sample_002, sample_003] (3 items grouped into 1)
    ```

`view()` gibt eine Ausgabe für jede Kanal-Emission zurück, also wissen wir, dass diese einzelne Ausgabe alle 3 ursprünglichen Elemente enthält, die in eine Liste gruppiert sind.

Jetzt lass uns die `collect`-Methode auf einer List in Aktion sehen. Modifiziere `collect.nf`, um die `collect`-Methode der List auf die ursprüngliche Liste der Proben-IDs anzuwenden:

=== "Danach"

    ```groovy title="main.nf" linenums="1" hl_lines="9-13"
    def sample_ids = ['sample_001', 'sample_002', 'sample_003']

    // channel.collect() - gruppiert mehrere Kanal-Emissionen in eine
    ch_input = channel.fromList(sample_ids)
    ch_input.view { sample -> "Individual channel item: ${sample}" }
    ch_collected = ch_input.collect()
    ch_collected.view { list -> "channel.collect() result: ${list} (${list.size()} items grouped into 1)" }

    // List.collect() - transformiert jedes Element, behält Struktur bei
    def formatted_ids = sample_ids.collect { id ->
        id.toUpperCase().replace('SAMPLE_', 'SPECIMEN_')
    }
    println "List.collect() result: ${formatted_ids} (${sample_ids.size()} items transformed into ${formatted_ids.size()})"
    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="1"
    def sample_ids = ['sample_001', 'sample_002', 'sample_003']

    // channel.collect() - gruppiert mehrere Kanal-Emissionen in eine
    ch_input = channel.fromList(sample_ids)
    ch_input.view { sample -> "Individual channel item: ${sample}" }
    ch_collected = ch_input.collect()
    ch_collected.view { list -> "channel.collect() result: ${list} (${list.size()} items grouped into 1)" }
    ```

In diesem neuen Snippet:

- Definieren wir eine neue Variable `formatted_ids`, die die `collect`-Methode der List verwendet, um jede Proben-ID in der ursprünglichen Liste zu transformieren
- Geben das Ergebnis mit `println` aus

Führe den modifizierten Workflow aus:

```bash
nextflow run collect.nf
```

??? success "Befehlsausgabe"

    ```console hl_lines="5"
    N E X T F L O W   ~  version 25.10.2

    Launching `collect.nf` [cheeky_stonebraker] DSL2 - revision: 2d5039fb47

    List.collect() result: [SPECIMEN_001, SPECIMEN_002, SPECIMEN_003] (3 items transformed into 3)
    Individual channel item: sample_001
    Individual channel item: sample_002
    Individual channel item: sample_003
    channel.collect() result: [sample_001, sample_002, sample_003] (3 items grouped into 1)
    ```

Diesmal haben wir die Struktur der Daten NICHT geändert, wir haben immer noch 3 Elemente in der Liste, aber wir HABEN jedes Element mit der `collect`-Methode der List transformiert, um eine neue Liste mit modifizierten Werten zu erzeugen. Dies ist ähnlich wie die Verwendung des `map`-Operators auf einem Kanal, aber es operiert auf einer List-Datenstruktur anstatt auf einem Kanal.

`collect` ist ein Extremfall, den wir hier verwenden, um einen Punkt zu verdeutlichen. Die wichtigste Lektion ist, dass du beim Schreiben von Workflows immer zwischen **Datenstrukturen** (Listen, Maps, etc.) und **Kanälen** (Datenfluss-Konstrukten) unterscheiden solltest. Operationen können Namen teilen, verhalten sich aber völlig unterschiedlich, je nachdem, auf welchem Typ sie aufgerufen werden.

### 1.3. Der Spread-Operator (`*.`) - Kurzform für Eigenschaftsextraktion

Verwandt mit der `collect`-Methode der List ist der Spread-Operator (`*.`), der eine prägnante Möglichkeit bietet, Eigenschaften aus Collections zu extrahieren. Er ist im Wesentlichen syntaktischer Zucker für ein gängiges `collect`-Muster.

Lass uns eine Demonstration zu unserer `collect.nf`-Datei hinzufügen:

=== "Danach"

    ```groovy title="collect.nf" linenums="1" hl_lines="15-18"
    def sample_ids = ['sample_001', 'sample_002', 'sample_003']

    // channel.collect() - gruppiert mehrere Kanal-Emissionen in eine
    ch_input = channel.fromList(sample_ids)
    ch_input.view { sample -> "Individual channel item: ${sample}" }
    ch_collected = ch_input.collect()
    ch_collected.view { list -> "channel.collect() result: ${list} (${list.size()} items grouped into 1)" }

    // List.collect() - transformiert jedes Element, behält Struktur bei
    def formatted_ids = sample_ids.collect { id ->
        id.toUpperCase().replace('SAMPLE_', 'SPECIMEN_')
    }
    println "List.collect() result: ${formatted_ids} (${sample_ids.size()} items transformed into ${formatted_ids.size()})"

    // Spread-Operator - prägnanter Eigenschaftszugriff
    def sample_data = [[id: 's1', quality: 38.5], [id: 's2', quality: 42.1], [id: 's3', quality: 35.2]]
    def all_ids = sample_data*.id
    println "Spread operator result: ${all_ids}"
    ```

=== "Vorher"

    ```groovy title="collect.nf" linenums="1"
    def sample_ids = ['sample_001', 'sample_002', 'sample_003']

    // channel.collect() - gruppiert mehrere Kanal-Emissionen in eine
    ch_input = channel.fromList(sample_ids)
    ch_input.view { sample -> "Individual channel item: ${sample}" }
    ch_collected = ch_input.collect()
    ch_collected.view { list -> "channel.collect() result: ${list} (${list.size()} items grouped into 1)" }

    // List.collect() - transformiert jedes Element, behält Struktur bei
    def formatted_ids = sample_ids.collect { id ->
        id.toUpperCase().replace('SAMPLE_', 'SPECIMEN_')
    }
    println "List.collect() result: ${formatted_ids} (${sample_ids.size()} items transformed into ${formatted_ids.size()})"
    ```

Führe den aktualisierten Workflow aus:

```bash title="Test spread operator"
nextflow run collect.nf
```

??? success "Befehlsausgabe"

    ```console hl_lines="6"
    N E X T F L O W   ~  version 25.10.2

    Launching `collect.nf` [cranky_galileo] DSL2 - revision: 5f3c8b2a91

    List.collect() result: [SPECIMEN_001, SPECIMEN_002, SPECIMEN_003] (3 items transformed into 3)
    Spread operator result: [s1, s2, s3]
    Individual channel item: sample_001
    Individual channel item: sample_002
    Individual channel item: sample_003
    channel.collect() result: [sample_001, sample_002, sample_003] (3 items grouped into 1)
    ```

Der Spread-Operator `*.` ist eine Kurzform für ein gängiges collect-Muster:

```groovy
// Diese sind äquivalent:
def ids = samples*.id
def ids = samples.collect { it.id }

// Funktioniert auch mit Methodenaufrufen:
def names = files*.getName()
def names = files.collect { it.getName() }
```

Der Spread-Operator ist besonders nützlich, wenn du eine einzelne Eigenschaft aus einer Liste von Objekten extrahieren musst – er ist lesbarer als das Ausschreiben der vollständigen `collect`-Closure.

!!! tip "Wann Spread vs. Collect verwenden"

    - **Verwende Spread (`*.`)** für einfachen Eigenschaftszugriff: `samples*.id`, `files*.name`
    - **Verwende collect** für Transformationen oder komplexe Logik: `samples.collect { it.id.toUpperCase() }`, `samples.collect { [it.id, it.quality > 40] }`

### Fazit

In diesem Abschnitt hast du gelernt:

- **Datenfluss vs. Scripting**: Kanal-Operatoren orchestrieren, wie Daten durch deine Pipeline fließen, während Scripting einzelne Datenelemente transformiert
- **Typen verstehen**: Derselbe Methodenname (wie `collect`) kann sich unterschiedlich verhalten, je nachdem, auf welchem Typ er aufgerufen wird (Channel vs. List)
- **Kontext ist wichtig**: Sei dir immer bewusst, ob du mit Kanälen (Datenfluss) oder Datenstrukturen (Scripting) arbeitest

Das Verstehen dieser Grenzen ist essenziell für Debugging, Dokumentation und das Schreiben wartbarer Workflows.

Als Nächstes werden wir tiefer in String-Verarbeitungsfähigkeiten eintauchen, die für den Umgang mit realen Daten unerlässlich sind.

---

## 2. String-Verarbeitung und dynamische Skript-Generierung

Die Beherrschung der String-Verarbeitung trennt fragile Workflows von robusten Pipelines. Dieser Abschnitt behandelt das Parsen komplexer Dateinamen, dynamische Skript-Generierung und Variableninterpolation.

### 2.1. Musterabgleich und reguläre Ausdrücke

Bioinformatik-Dateien haben oft komplexe Benennungskonventionen, die Metadaten kodieren. Lass uns diese automatisch mit Musterabgleich unter Verwendung regulärer Ausdrücke extrahieren.

Wir kehren zu unserem `main.nf`-Workflow zurück und fügen etwas Musterabgleich-Logik hinzu, um zusätzliche Probeninformationen aus Dateinamen zu extrahieren. Die FASTQ-Dateien in unserem Datensatz folgen Illumina-Benennungskonventionen mit Namen wie `SAMPLE_001_S1_L001_R1_001.fastq.gz`. Diese mögen kryptisch aussehen, aber sie kodieren tatsächlich nützliche Metadaten wie Proben-ID, Lane-Nummer und Read-Richtung. Wir werden Regex-Fähigkeiten verwenden, um diese Namen zu parsen.

Nimm folgende Änderung an deinem bestehenden `main.nf`-Workflow vor:

=== "Danach"

    ```groovy title="main.nf" linenums="4" hl_lines="10-21"
            .map { row ->
                // Scripting für Datentransformation
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
    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="4" hl_lines="10-11"
            .map { row ->
                // Scripting für Datentransformation
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                def priority = sample_meta.quality > 40 ? 'high' : 'normal'
                return tuple(sample_meta + [priority: priority], file(row.file_path))
            }
    ```

Dies demonstriert wichtige **String-Verarbeitungskonzepte**:

1. **Reguläre Ausdrücke** mit `~/pattern/`-Syntax - dies erstellt ein Regex-Muster, ohne Backslashes escapen zu müssen
2. **Musterabgleich** mit dem `=~`-Operator - dies versucht, einen String mit einem Regex-Muster abzugleichen
3. **Matcher-Objekte**, die Gruppen mit `[0][1]`, `[0][2]`, etc. erfassen - `[0]` bezieht sich auf den gesamten Match, `[1]`, `[2]`, etc. beziehen sich auf erfasste Gruppen in Klammern

Lass uns das Regex-Muster `^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$` aufschlüsseln:

| Muster              | Passt auf                                  | Erfasst                                |
| ------------------- | ------------------------------------------ | -------------------------------------- |
| `^(.+)`             | Probenname vom Anfang                      | Gruppe 1: Probenname                   |
| `_S(\d+)`           | Probennummer `_S1`, `_S2`, etc.            | Gruppe 2: Probennummer                 |
| `_L(\d{3})`         | Lane-Nummer `_L001`                        | Gruppe 3: Lane (3 Ziffern)             |
| `_(R[12])`          | Read-Richtung `_R1` oder `_R2`             | Gruppe 4: Read-Richtung                |
| `_(\d{3})`          | Chunk-Nummer `_001`                        | Gruppe 5: Chunk (3 Ziffern)            |
| `\.fastq(?:\.gz)?$` | Dateierweiterung `.fastq` oder `.fastq.gz` | Nicht erfasst (?: ist nicht-erfassend) |

Dies parst Illumina-Benennungskonventionen, um Metadaten automatisch zu extrahieren.

Führe den modifizierten Workflow aus:

```bash title="Test pattern matching"
nextflow run main.nf
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [clever_pauling] DSL2 - revision: 605d2058b4

    [[id:sample_001, organism:human, tissue:liver, depth:30000000, quality:38.5, sample_num:1, lane:001, read:R1, chunk:001, priority:normal], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_001_S1_L001_R1_001.fastq]
    [[id:sample_002, organism:mouse, tissue:brain, depth:25000000, quality:35.2, sample_num:2, lane:001, read:R1, chunk:001, priority:normal], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_002_S2_L001_R1_001.fastq]
    [[id:sample_003, organism:human, tissue:kidney, depth:45000000, quality:42.1, sample_num:3, lane:001, read:R1, chunk:001, priority:high], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_003_S3_L001_R1_001.fastq]
    ```

Dies zeigt die aus den Dateinamen angereicherten Metadaten.

### 2.2. Dynamische Skript-Generierung in Prozessen

Prozess-Skript-Blöcke sind im Wesentlichen mehrzeilige Strings, die an die Shell übergeben werden. Du kannst **bedingte Logik** (if/else, ternäre Operatoren) verwenden, um dynamisch verschiedene Skript-Strings basierend auf Eingabeeigenschaften zu generieren. Dies ist essenziell für den Umgang mit verschiedenen Eingabetypen – wie Single-End- vs. Paired-End-Sequenzierungs-Reads – ohne Prozessdefinitionen zu duplizieren.

Lass uns einen Prozess zu unserem Workflow hinzufügen, der dieses Muster demonstriert. Öffne `modules/fastp.nf` und schau es dir an:

```groovy title="modules/fastp.nf" linenums="1"
process FASTP {
    container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_trimmed*.fastq.gz"), emit: reads

    script:
    """
    fastp \\
        --in1 ${reads[0]} \\
        --in2 ${reads[1]} \\
        --out1 ${meta.id}_trimmed_R1.fastq.gz \\
        --out2 ${meta.id}_trimmed_R2.fastq.gz \\
        --json ${meta.id}.fastp.json \\
        --html ${meta.id}.fastp.html \\
        --thread $task.cpus
    """
}
```

Der Prozess nimmt FASTQ-Dateien als Eingabe und führt das `fastp`-Tool aus, um Adapter zu trimmen und Reads niedriger Qualität zu filtern. Leider hat die Person, die diesen Prozess geschrieben hat, die Single-End-Reads, die wir in unserem Beispieldatensatz haben, nicht berücksichtigt. Lass uns ihn zu unserem Workflow hinzufügen und sehen, was passiert:

Füge zunächst das Modul in der allerersten Zeile deines `main.nf`-Workflows ein:

```groovy title="main.nf" linenums="1"
include { FASTP } from './modules/fastp.nf'
```

Modifiziere dann den `workflow`-Block, um den `ch_samples`-Kanal mit dem `FASTP`-Prozess zu verbinden:

=== "Danach"

    ```groovy title="main.nf" linenums="25" hl_lines="27"
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
    }
    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="25" hl_lines="26"
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
                return [sample_meta + file_meta + [priority: priority], file(row.file_path)]
            }
            .view()
    }
    ```

Führe diesen modifizierten Workflow aus:

```bash
nextflow run main.nf
```

??? failure "Befehlsausgabe"

    ```console
    ERROR ~ Error executing process > 'FASTP (3)'

    Caused by:
      Process `FASTP (3)` terminated with an error exit status (255)


    Command executed:

      fastp \
          --in1 SAMPLE_003_S3_L001_R1_001.fastq \
          --in2 null \
          --out1 sample_003_trimmed_R1.fastq.gz \
          --out2 sample_003_trimmed_R2.fastq.gz \
          --json sample_003.fastp.json \
          --html sample_003.fastp.html \
          --thread 2

    Command exit status:
      255

    Command output:
      (empty)
    ```

Du kannst sehen, dass der Prozess versucht, `fastp` mit einem `null`-Wert für die zweite Eingabedatei auszuführen, was zum Fehlschlagen führt. Das liegt daran, dass unser Datensatz Single-End-Reads enthält, aber der Prozess fest codiert ist, um Paired-End-Reads zu erwarten (zwei Eingabedateien gleichzeitig).

Behebe dies, indem du bedingte Logik zum `script:`-Block des `FASTP`-Prozesses hinzufügst. Eine if/else-Anweisung prüft die Anzahl der Read-Dateien und passt den Befehl entsprechend an.

=== "Danach"

    ```groovy title="main.nf" linenums="10" hl_lines="3-27"
        script:
        // Einfache Single-End- vs. Paired-End-Erkennung
        def is_single = reads instanceof List ? reads.size() == 1 : true

        if (is_single) {
            def input_file = reads instanceof List ? reads[0] : reads
            """
            fastp \\
                --in1 ${input_file} \\
                --out1 ${meta.id}_trimmed.fastq.gz \\
                --json ${meta.id}.fastp.json \\
                --html ${meta.id}.fastp.html \\
                --thread $task.cpus
            """
        } else {
            """
            fastp \\
                --in1 ${reads[0]} \\
                --in2 ${reads[1]} \\
                --out1 ${meta.id}_trimmed_R1.fastq.gz \\
                --out2 ${meta.id}_trimmed_R2.fastq.gz \\
                --json ${meta.id}.fastp.json \\
                --html ${meta.id}.fastp.html \\
                --thread $task.cpus
            """
        }
    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="10" hl_lines="2-11"
            script:
            """
            fastp \\
                --in1 ${reads[0]} \\
                --in2 ${reads[1]} \\
                --out1 ${meta.id}_trimmed_R1.fastq.gz \\
                --out2 ${meta.id}_trimmed_R2.fastq.gz \\
                --json ${meta.id}.fastp.json \\
                --html ${meta.id}.fastp.html \\
                --thread $task.cpus
            """
        }
    ```

Jetzt kann der Workflow sowohl Single-End- als auch Paired-End-Reads elegant handhaben. Die bedingte Logik prüft die Anzahl der Eingabedateien und konstruiert den entsprechenden Befehl für `fastp`. Lass uns sehen, ob es funktioniert:

```bash
nextflow run main.nf
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [adoring_rosalind] DSL2 - revision: 04b1cd93e9

    executor >  local (3)
    [31/a8ad4d] process > FASTP (3) [100%] 3 of 3 ✔
    ```

Sieht gut aus! Wenn wir die tatsächlich ausgeführten Befehle überprüfen (passe für deinen Task-Hash an):

```console title="Check commands executed"
cat work/31/a8ad4d95749e685a6d842d3007957f/.command.sh
```

Wir können sehen, dass Nextflow korrekt den richtigen Befehl für Single-End-Reads ausgewählt hat:

```bash title=".command.sh"
#!/bin/bash -ue
fastp \
    --in1 SAMPLE_003_S3_L001_R1_001.fastq \
    --out1 sample_003_trimmed.fastq.gz \
    --json sample_003.fastp.json \
    --html sample_003.fastp.html \
    --thread 2
```

Eine weitere häufige Verwendung dynamischer Skript-Logik ist im [Nextflow for Science Genomics-Modul](../../nf4science/genomics/02_joint_calling) zu sehen. In diesem Modul kann der aufgerufene GATK-Prozess mehrere Eingabedateien annehmen, aber jede muss mit `-V` präfixiert werden, um eine korrekte Befehlszeile zu bilden. Der Prozess verwendet Scripting, um eine Collection von Eingabedateien (`all_gvcfs`) in die korrekten Befehlsargumente zu transformieren:

```groovy title="command line manipulation for GATK" linenums="1" hl_lines="2 5"
    script:
    def gvcfs_line = all_gvcfs.collect { gvcf -> "-V ${gvcf}" }.join(' ')
    """
    gatk GenomicsDBImport \
        ${gvcfs_line} \
        -L ${interval_list} \
        --genomicsdb-workspace-path ${cohort_name}_gdb
    """
```

Diese Muster der Verwendung von Scripting in Prozess-Skript-Blöcken sind extrem leistungsfähig und können in vielen Szenarien angewendet werden – vom Umgang mit variablen Eingabetypen bis zum Erstellen komplexer Befehlszeilenargumente aus Datei-Collections, wodurch deine Prozesse wirklich an die vielfältigen Anforderungen realer Daten anpassbar werden.

### 2.3. Variableninterpolation: Nextflow- und Shell-Variablen

Prozess-Skripte mischen Nextflow-Variablen, Shell-Variablen und Befehlssubstitutionen, jede mit unterschiedlicher Interpolationssyntax. Die Verwendung der falschen Syntax verursacht Fehler. Lass uns diese mit einem Prozess erkunden, der einen Verarbeitungsbericht erstellt.

Schau dir die Moduldatei `modules/generate_report.nf` an:

```groovy title="modules/generate_report.nf" linenums="1"
process GENERATE_REPORT {

    publishDir 'results/reports', mode: 'copy'

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

Dieser Prozess schreibt einen einfachen Bericht mit der Proben-ID und dem Dateinamen. Jetzt lass uns ihn ausführen, um zu sehen, was passiert, wenn wir verschiedene Arten von Variablen mischen müssen.

Füge den Prozess in deine `main.nf` ein und füge ihn zum Workflow hinzu:

=== "Danach"

    ```groovy title="main.nf" linenums="1" hl_lines="2 30"
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

=== "Vorher"

    ```groovy title="main.nf" linenums="1" hl_lines="1 10-29"
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
    }
    ```

Führe jetzt den Workflow aus und überprüfe die generierten Berichte in `results/reports/`. Sie sollten grundlegende Informationen über jede Probe enthalten.

<!-- TODO: add the run command -->

??? success "Befehlsausgabe"

    ```console
    <!-- TODO: output -->
    ```

Aber was, wenn wir Informationen darüber hinzufügen möchten, wann und wo die Verarbeitung stattgefunden hat? Lass uns den Prozess modifizieren, um **Shell**-Variablen und etwas Befehlssubstitution zu verwenden, um den aktuellen Benutzer, Hostnamen und das Datum in den Bericht aufzunehmen:

=== "Danach"

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

=== "Vorher"

    ```groovy title="modules/generate_report.nf" linenums="10"
        script:
        """
        echo "Processing ${reads}" > ${meta.id}_report.txt
        echo "Sample: ${meta.id}" >> ${meta.id}_report.txt
        """
    ```

Wenn du dies ausführst, wirst du einen Fehler bemerken – Nextflow versucht, `${USER}` als Nextflow-Variable zu interpretieren, die nicht existiert.

??? failure "Befehlsausgabe"

    ```console
    Error modules/generate_report.nf:15:27: `USER` is not defined
    │  15 |     echo "Processed by: ${USER}" >> ${meta.id}_report.txt
    ╰     |                           ^^^^

    ERROR ~ Script compilation failed
    ```

Wir müssen es escapen, damit Bash es stattdessen handhaben kann.

Behebe dies, indem du die Shell-Variablen und Befehlssubstitutionen mit einem Backslash (`\`) escapest:

=== "Danach"

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

=== "Vorher"

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

Jetzt funktioniert es! Der Backslash (`\`) sagt Nextflow „interpretiere das nicht, gib es an Bash weiter."

### Fazit

In diesem Abschnitt hast du **String-Verarbeitungs**-Techniken gelernt:

- **Reguläre Ausdrücke für Datei-Parsing**: Verwendung des `=~`-Operators und Regex-Mustern (`~/pattern/`), um Metadaten aus komplexen Dateibenennungskonventionen zu extrahieren
- **Dynamische Skript-Generierung**: Verwendung bedingter Logik (if/else, ternäre Operatoren), um verschiedene Skript-Strings basierend auf Eingabeeigenschaften zu generieren
- **Variableninterpolation**: Verstehen, wann Nextflow Strings interpretiert vs. wann die Shell es tut
  - `${var}` - Nextflow-Variablen (von Nextflow zur Workflow-Kompilierzeit interpoliert)
  - `\${var}` - Shell-Umgebungsvariablen (escaped, zur Laufzeit an Bash übergeben)
  - `\$(cmd)` - Shell-Befehlssubstitution (escaped, von Bash zur Laufzeit ausgeführt)

Diese String-Verarbeitungs- und Generierungsmuster sind essenziell für den Umgang mit den vielfältigen Dateiformaten und Benennungskonventionen, denen du in realen Bioinformatik-Workflows begegnen wirst.

---

## 3. Wiederverwendbare Funktionen erstellen

Komplexe Workflow-Logik inline in Kanal-Operatoren oder Prozessdefinitionen reduziert Lesbarkeit und Wartbarkeit. **Funktionen** ermöglichen es dir, diese Logik in benannte, wiederverwendbare Komponenten zu extrahieren.

Unsere Map-Operation ist lang und komplex geworden. Lass uns sie in eine wiederverwendbare Funktion mit dem `def`-Schlüsselwort extrahieren.

Um zu veranschaulichen, wie das mit unserem bestehenden Workflow aussieht, nimm die folgende Änderung vor und verwende `def`, um eine wiederverwendbare Funktion namens `separateMetadata` zu definieren:

=== "Danach"

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

=== "Vorher"

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

Durch das Extrahieren dieser Logik in eine Funktion haben wir die tatsächliche Workflow-Logik auf etwas viel Saubereres reduziert:

```groovy title="minimal workflow"
    ch_samples = channel.fromPath("./data/samples.csv")
        .splitCsv(header: true)
        .map{ row -> separateMetadata(row) }

    ch_fastp = FASTP(ch_samples)
    GENERATE_REPORT(ch_samples)
```

Dies macht die Workflow-Logik viel einfacher zu lesen und auf einen Blick zu verstehen. Die Funktion `separateMetadata` kapselt die gesamte komplexe Logik zum Parsen und Anreichern von Metadaten und macht sie wiederverwendbar und testbar.

Führe den Workflow aus, um sicherzustellen, dass er noch funktioniert:

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

Die Ausgabe sollte zeigen, dass beide Prozesse erfolgreich abgeschlossen wurden. Der Workflow ist jetzt viel sauberer und einfacher zu warten, wobei die gesamte komplexe Metadaten-Verarbeitungslogik in der `separateMetadata`-Funktion gekapselt ist.

### Fazit

In diesem Abschnitt hast du **Funktionserstellung** gelernt:

- **Funktionen mit `def` definieren**: Das Schlüsselwort zum Erstellen benannter Funktionen (wie `def` in Python oder `function` in JavaScript)
- **Funktionsbereich**: Auf Skriptebene definierte Funktionen sind in deinem gesamten Nextflow-Workflow zugänglich
- **Rückgabewerte**: Funktionen geben automatisch den letzten Ausdruck zurück oder verwenden explizites `return`
- **Saubererer Code**: Das Extrahieren komplexer Logik in Funktionen ist eine grundlegende Software-Engineering-Praxis in jeder Sprache

Als Nächstes werden wir erkunden, wie man Closures in Prozess-Direktiven für dynamische Ressourcenzuweisung verwendet.

---

## 4. Dynamische Ressourcen-Direktiven mit Closures

Bisher haben wir Scripting im `script`-Block von Prozessen verwendet. Aber **Closures** (eingeführt in Abschnitt 1.1) sind auch unglaublich nützlich in Prozess-Direktiven, besonders für dynamische Ressourcenzuweisung. Lass uns Ressourcen-Direktiven zu unserem FASTP-Prozess hinzufügen, die sich basierend auf den Probeneigenschaften anpassen.

### 4.1. Probenspezifische Ressourcenzuweisung

Derzeit verwendet unser FASTP-Prozess Standard-Ressourcen. Lass uns ihn intelligenter machen, indem wir mehr CPUs für Proben mit hoher Tiefe zuweisen. Bearbeite `modules/fastp.nf`, um eine dynamische `cpus`-Direktive und eine statische `memory`-Direktive einzufügen:

=== "Danach"

    ```groovy title="modules/fastp.nf" linenums="1" hl_lines="4-5"
    process FASTP {
        container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

        cpus { meta.depth > 40000000 ? 2 : 1 }
        memory 2.GB

        input:
        tuple val(meta), path(reads)
    ```

=== "Vorher"

    ```groovy title="modules/fastp.nf" linenums="1"
    process FASTP {
        container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

        input:
        tuple val(meta), path(reads)
    ```

Die Closure `{ meta.depth > 40000000 ? 2 : 1 }` verwendet den **ternären Operator** (behandelt in Abschnitt 1.1) und wird für jede Aufgabe ausgewertet, was probenspezifische Ressourcenzuweisung ermöglicht. Proben mit hoher Tiefe (>40M Reads) erhalten 2 CPUs, während andere 1 CPU erhalten.

!!! note "Zugriff auf Eingabevariablen in Direktiven"

    Die Closure kann auf alle Eingabevariablen (wie `meta` hier) zugreifen, weil Nextflow diese Closures im Kontext jeder Aufgabenausführung auswertet.

Führe den Workflow erneut mit der Option `-ansi-log false` aus, um die Task-Hashes leichter sehen zu können.

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

Du solltest etwas wie das sehen:

```bash title="docker command"
    docker run -i --cpu-shares 4096 --memory 2048m -e "NXF_TASK_WORKDIR" -v /workspaces/training/side-quests/essential_scripting_patterns:/workspaces/training/side-quests/essential_scripting_patterns -w "$NXF_TASK_WORKDIR" --name $NXF_BOXID community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690 /bin/bash -ue /workspaces/training/side-quests/essential_scripting_patterns/work/48/6db0c9e9d8aa65e4bb4936cd3bd59e/.command.sh
```

In diesem Beispiel haben wir ein Beispiel gewählt, das 2 CPUs angefordert hat (`--cpu-shares 2048`), weil es eine Probe mit hoher Tiefe war, aber du solltest unterschiedliche CPU-Zuweisungen je nach Probentiefe sehen. Probiere dies auch für die anderen Aufgaben aus.

### 4.2. Wiederholungsstrategien

Ein weiteres leistungsstarkes Muster ist die Verwendung von `task.attempt` für Wiederholungsstrategien. Um zu zeigen, warum dies nützlich ist, werden wir zunächst die Speicherzuweisung für FASTP auf weniger reduzieren, als es benötigt. Ändere die `memory`-Direktive in `modules/fastp.nf` auf `1.GB`:

=== "Danach"

    ```groovy title="modules/fastp.nf" linenums="1" hl_lines="5"
    process FASTP {
        container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

        cpus { meta.depth > 40000000 ? 4 : 2 }
        memory 1.GB

        input:
        tuple val(meta), path(reads)
    ```

=== "Vorher"

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

Dies zeigt an, dass der Prozess wegen Überschreitung der Speichergrenzen beendet wurde.

Dies ist ein sehr häufiges Szenario in realen Workflows – manchmal weißt du einfach nicht, wie viel Speicher eine Aufgabe benötigt, bis du sie ausführst.

Um unseren Workflow robuster zu machen, können wir eine Wiederholungsstrategie implementieren, die die Speicherzuweisung bei jedem Versuch erhöht, wieder einmal mit einer Groovy-Closure. Modifiziere die `memory`-Direktive, um den Basisspeicher mit `task.attempt` zu multiplizieren, und füge die Direktiven `errorStrategy 'retry'` und `maxRetries 2` hinzu:

=== "Danach"

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

=== "Vorher"

    ```groovy title="modules/fastp.nf" linenums="1" hl_lines="5"
    process FASTP {
        container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

        cpus { meta.depth > 40000000 ? 4 : 2 }
        memory 2.GB

        input:
        tuple val(meta), path(reads)
    ```

Jetzt wird Nextflow, wenn der Prozess aufgrund unzureichenden Speichers fehlschlägt, mit mehr Speicher wiederholen:

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

## 5. Bedingte Logik und Prozesskontrolle

Zuvor haben wir `.map()` mit Scripting verwendet, um Kanaldaten zu transformieren. Jetzt werden wir bedingte Logik verwenden, um zu steuern, welche Prozesse basierend auf Daten ausgeführt werden – essenziell für flexible Workflows, die sich an verschiedene Probentypen anpassen.

Nextflows [Datenfluss-Operatoren](https://www.nextflow.io/docs/latest/reference/operator.html) nehmen Closures, die zur Laufzeit ausgewertet werden, und ermöglichen bedingte Logik, um Workflow-Entscheidungen basierend auf Kanalinhalten zu treffen.

### 5.1. Routing mit `.branch()`

Nehmen wir zum Beispiel an, dass unsere Sequenzierungsproben nur mit FASTP getrimmt werden müssen, wenn es sich um menschliche Proben mit einer Abdeckung über einem bestimmten Schwellenwert handelt. Maus-Proben oder Proben mit niedriger Abdeckung sollten stattdessen mit Trimgalore ausgeführt werden (dies ist ein konstruiertes Beispiel, aber es veranschaulicht den Punkt).

Wir haben einen einfachen Trimgalore-Prozess in `modules/trimgalore.nf` bereitgestellt, schau ihn dir an, wenn du möchtest, aber die Details sind für diese Übung nicht wichtig. Der Schlüsselpunkt ist, dass wir Proben basierend auf ihren Metadaten routen möchten.

Füge das neue Modul aus `modules/trimgalore.nf` ein:

=== "Danach"

    ```groovy title="main.nf" linenums="1" hl_lines="2"
    include { FASTP } from './modules/fastp.nf'
    include { TRIMGALORE } from './modules/trimgalore.nf'
    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="1"
    include { FASTP } from './modules/fastp.nf'
    ```

... und modifiziere dann deinen `main.nf`-Workflow, um Proben basierend auf ihren Metadaten zu verzweigen und sie durch den entsprechenden Trimming-Prozess zu routen, so:

=== "Danach"

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

=== "Vorher"

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

Hier haben wir kleine, aber mächtige bedingte Ausdrücke innerhalb des `.branch{}`-Operators verwendet, um Proben basierend auf ihren Metadaten zu routen. Menschliche Proben mit hoher Abdeckung gehen durch `FASTP`, während alle anderen Proben durch `TRIMGALORE` gehen.

### 5.2. Verwendung von `.filter()` mit Truthiness

Ein weiteres leistungsstarkes Muster zur Steuerung der Workflow-Ausführung ist der `.filter()`-Operator, der eine Closure verwendet, um zu bestimmen, welche Elemente in der Pipeline fortfahren sollen. Innerhalb der Filter-Closure schreibst du **boolesche Ausdrücke**, die entscheiden, welche Elemente durchkommen.

Nextflow (wie viele dynamische Sprachen) hat ein Konzept von **"Truthiness"**, das bestimmt, welche Werte in booleschen Kontexten zu `true` oder `false` ausgewertet werden:

- **Truthy**: Nicht-null-Werte, nicht-leere Strings, Nicht-Null-Zahlen, nicht-leere Collections
- **Falsy**: `null`, leere Strings `""`, Null `0`, leere Collections `[]` oder `[:]`, `false`

Das bedeutet, dass `meta.id` allein (ohne explizites `!= null`) prüft, ob die ID existiert und nicht leer ist. Lass uns dies verwenden, um Proben herauszufiltern, die unsere Qualitätsanforderungen nicht erfüllen.

Füge Folgendes vor der Branch-Operation hinzu:

=== "Danach"

    ```groovy title="main.nf" linenums="28" hl_lines="5-11"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row -> separateMetadata(row) }

        // Ungültige oder qualitativ minderwertige Proben herausfiltern
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

=== "Vorher"

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

Der Filterausdruck `meta.id && meta.organism && meta.depth >= 25000000` kombiniert Truthiness mit expliziten Vergleichen:

- `meta.id && meta.organism` prüft, dass beide Felder existieren und nicht leer sind (mit Truthiness)
- `meta.depth >= 25000000` stellt ausreichende Sequenzierungstiefe mit einem expliziten Vergleich sicher

!!! note "Truthiness in der Praxis"

    Der Ausdruck `meta.id && meta.organism` ist prägnanter als zu schreiben:
    ```groovy
    meta.id != null && meta.id != '' && meta.organism != null && meta.organism != ''
    ```

    Dies macht die Filterlogik viel sauberer und einfacher zu lesen.

### Fazit

In diesem Abschnitt hast du gelernt, bedingte Logik zu verwenden, um die Workflow-Ausführung mit den Closure-Schnittstellen von Nextflow-Operatoren wie `.branch{}` und `.filter{}` zu steuern und Truthiness zu nutzen, um prägnante bedingte Ausdrücke zu schreiben.

Unsere Pipeline routet jetzt intelligent Proben durch geeignete Prozesse, aber Produktions-Workflows müssen ungültige Daten elegant handhaben. Lass uns unseren Workflow robust gegen fehlende oder null-Werte machen.

---

## 6. Sichere Navigation und Elvis-Operatoren

Unsere `separateMetadata`-Funktion geht derzeit davon aus, dass alle CSV-Felder vorhanden und gültig sind. Aber was passiert mit unvollständigen Daten? Lass es uns herausfinden.

### 6.1. Das Problem: Zugriff auf Eigenschaften, die nicht existieren

Nehmen wir an, wir möchten Unterstützung für optionale Sequenzierungs-Run-Informationen hinzufügen. In einigen Laboren haben Proben möglicherweise ein zusätzliches Feld für die Sequenzierungs-Run-ID oder Batch-Nummer, aber unsere aktuelle CSV hat diese Spalte nicht. Lass uns trotzdem versuchen, darauf zuzugreifen.

Modifiziere die `separateMetadata`-Funktion, um ein run_id-Feld einzuschließen:

=== "Danach"

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

=== "Vorher"

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

Führe jetzt den Workflow aus:

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

Das Problem ist, dass `row.run_id` `null` zurückgibt, weil die `run_id`-Spalte in unserer CSV nicht existiert. Wenn wir versuchen, `.toUpperCase()` auf `null` aufzurufen, stürzt es ab. Hier rettet der sichere Navigationsoperator den Tag.

### 6.2. Sicherer Navigationsoperator (`?.`)

Der sichere Navigationsoperator (`?.`) gibt `null` zurück, anstatt eine Exception zu werfen, wenn er auf einem `null`-Wert aufgerufen wird. Wenn das Objekt vor `?.` `null` ist, wird der gesamte Ausdruck zu `null` ausgewertet, ohne die Methode auszuführen.

Aktualisiere die Funktion, um sichere Navigation zu verwenden:

=== "Danach"

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

=== "Vorher"

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

Führe erneut aus:

```bash
nextflow run main.nf
```

??? success "Befehlsausgabe"

    ```console
    <!-- TODO: output -->
    ```

Kein Absturz! Der Workflow behandelt jetzt das fehlende Feld elegant. Wenn `row.run_id` `null` ist, verhindert der `?.`-Operator den `.toUpperCase()`-Aufruf, und `run_id` wird `null`, anstatt eine Exception zu verursachen.

### 6.3. Elvis-Operator (`?:`) für Standardwerte

Der Elvis-Operator (`?:`) liefert Standardwerte, wenn die linke Seite "falsy" ist (wie zuvor erklärt). Er ist nach Elvis Presley benannt, weil `?:` wie seine berühmten Haare und Augen aussieht, wenn man ihn seitlich betrachtet!

Jetzt, da wir sichere Navigation verwenden, wird `run_id` für Proben ohne dieses Feld `null` sein. Lass uns den Elvis-Operator verwenden, um einen Standardwert bereitzustellen und ihn zu unserer `sample_meta`-Map hinzuzufügen:

=== "Danach"

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

=== "Vorher"

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

=== "Danach"

    ```groovy title="main.nf" linenums="30" hl_lines="4"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map{ row -> separateMetadata(row) }
            .view()
    ```

=== "Vorher"

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

Perfekt! Jetzt haben alle Proben ein `run`-Feld mit entweder ihrer tatsächlichen Run-ID (in Großbuchstaben) oder dem Standardwert 'UNSPECIFIED'. Die Kombination von `?.` und `?:` bietet sowohl Sicherheit (keine Abstürze) als auch sinnvolle Standardwerte.

Entferne jetzt den `.view()`-Operator, nachdem wir bestätigt haben, dass es funktioniert.

!!! tip "Kombination von sicherer Navigation und Elvis"

    Das Muster `value?.method() ?: 'default'` ist in Produktions-Workflows üblich:

    - `value?.method()` - Ruft Methode sicher auf, gibt `null` zurück, wenn `value` `null` ist
    - `?: 'default'` - Liefert Fallback, wenn Ergebnis `null` ist

    Dieses Muster behandelt fehlende/unvollständige Daten elegant.

Verwende diese Operatoren konsistent in Funktionen, Operator-Closures (`.map{}`, `.filter{}`), Prozess-Skripten und Config-Dateien. Sie verhindern Abstürze beim Umgang mit realen Daten.

### Fazit

- **Sichere Navigation (`?.`)**: Verhindert Abstürze bei null-Werten - gibt null zurück statt Exception zu werfen
- **Elvis-Operator (`?:`)**: Liefert Standardwerte - `value ?: 'default'`
- **Kombination**: `value?.method() ?: 'default'` ist das gängige Muster

Diese Operatoren machen Workflows resilient gegen unvollständige Daten – essenziell für reale Arbeit.

---

## 7. Validierung mit `error()` und `log.warn`

Manchmal musst du den Workflow sofort stoppen, wenn Eingabeparameter ungültig sind. In Nextflow kannst du eingebaute Funktionen wie `error()` und `log.warn` sowie Standard-Programmierkonzepte wie `if`-Anweisungen und boolesche Logik verwenden, um Validierungslogik zu implementieren. Lass uns Validierung zu unserem Workflow hinzufügen.

Erstelle eine Validierungsfunktion vor deinem Workflow-Block, rufe sie vom Workflow auf und ändere die Kanal-Erstellung, um einen Parameter für den CSV-Dateipfad zu verwenden. Wenn der Parameter fehlt oder die Datei nicht existiert, rufe `error()` auf, um die Ausführung mit einer klaren Nachricht zu stoppen.

=== "Danach"

    ```groovy title="main.nf" linenums="1" hl_lines="5-20 23-24"
    include { FASTP } from './modules/fastp.nf'
    include { TRIMGALORE } from './modules/trimgalore.nf'
    include { GENERATE_REPORT } from './modules/generate_report.nf'

    def validateInputs() {
        // Prüfe, ob Eingabeparameter bereitgestellt wurde
        if (!params.input) {
            error("Input CSV file path not provided. Please specify --input <file.csv>")
        }

        // Prüfe, ob CSV-Datei existiert
        if (!file(params.input).exists()) {
            error("Input CSV file not found: ${params.input}")
        }
    }
    ...
    workflow {
        validateInputs()
        ch_samples = channel.fromPath(params.input)
    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="1"
    include { FASTP } from './modules/fastp.nf'
    include { TRIMGALORE } from './modules/trimgalore.nf'
    include { GENERATE_REPORT } from './modules/generate_report.nf'

    ...
    workflow {
        ch_samples = channel.fromPath("./data/samples.csv")
    ```

Versuche jetzt, ohne die CSV-Datei auszuführen:

```bash
nextflow run main.nf
```

??? failure "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [confident_coulomb] DSL2 - revision: 07059399ed

    WARN: Access to undefined parameter `input` -- Initialise it to a default value eg. `params.input = some_value`
    Input CSV file path not provided. Please specify --input <file.csv>
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

    Input CSV file not found: ./data/nonexistent.csv
    ```

Führe ihn schließlich mit der korrekten Datei aus:

```bash
nextflow run main.nf --input ./data/samples.csv
```

??? success "Befehlsausgabe"

    ```console
    <!-- TODO: output -->
    ```

Diesmal läuft er erfolgreich.

Du kannst auch Validierung innerhalb der `separateMetadata`-Funktion hinzufügen. Lass uns das nicht-fatale `log.warn` verwenden, um Warnungen für Proben mit niedriger Sequenzierungstiefe auszugeben, aber den Workflow trotzdem fortsetzen lassen:

=== "Danach"

    ```groovy title="main.nf" linenums="1" hl_lines="3-6"
        def priority = sample_meta.quality > 40 ? 'high' : 'normal'

        // Validiere, dass Daten Sinn ergeben
        if (sample_meta.depth < 30000000) {
            log.warn "Low sequencing depth for ${sample_meta.id}: ${sample_meta.depth}"
        }

        return tuple(sample_meta + file_meta + [priority: priority], fastq_path)
    }
    ```

=== "Vorher"

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
    WARN: Low sequencing depth for sample_002: 25000000
    ```

Wir sehen eine Warnung über niedrige Sequenzierungstiefe für eine der Proben.

### Fazit

- **`error()`**: Stoppt Workflow sofort mit klarer Nachricht
- **`log.warn`**: Gibt Warnungen aus, ohne Workflow zu stoppen
- **Frühe Validierung**: Prüfe Eingaben vor der Verarbeitung, um schnell mit hilfreichen Fehlern zu scheitern
- **Validierungsfunktionen**: Erstelle wiederverwendbare Validierungslogik, die beim Workflow-Start aufgerufen werden kann

Ordnungsgemäße Validierung macht Workflows robuster und benutzerfreundlicher, indem Probleme früh mit klaren Fehlermeldungen erkannt werden.

---

## 8. Workflow-Event-Handler

Bis jetzt haben wir Code in unseren Workflow-Skripten und Prozessdefinitionen geschrieben. Aber es gibt noch ein wichtiges Feature, das du kennen solltest: Workflow-Event-Handler.

Event-Handler sind Closures, die an bestimmten Punkten im Lebenszyklus deines Workflows ausgeführt werden. Sie sind perfekt zum Hinzufügen von Logging, Benachrichtigungen oder Aufräumoperationen. Diese Handler sollten in deinem Workflow-Skript neben deiner Workflow-Definition definiert werden.

### 8.1. Der `onComplete`-Handler

Der am häufigsten verwendete Event-Handler ist `onComplete`, der ausgeführt wird, wenn dein Workflow endet (egal ob erfolgreich oder fehlgeschlagen). Lass uns einen hinzufügen, um unsere Pipeline-Ergebnisse zusammenzufassen.

Füge den Event-Handler zu deiner `main.nf`-Datei hinzu, innerhalb deiner Workflow-Definition:

=== "Danach"

    ```groovy title="main.nf" linenums="66" hl_lines="5-16"
        ch_fastp = FASTP(trim_branches.fastp)
        ch_trimgalore = TRIMGALORE(trim_branches.trimgalore)
        GENERATE_REPORT(ch_samples)

        workflow.onComplete = {
            println ""
            println "Pipeline execution summary:"
            println "=========================="
            println "Completed at: ${workflow.complete}"
            println "Duration    : ${workflow.duration}"
            println "Success     : ${workflow.success}"
            println "workDir     : ${workflow.workDir}"
            println "exit status : ${workflow.exitStatus}"
            println ""
        }
    }
    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="66" hl_lines="4"
        ch_fastp = FASTP(trim_branches.fastp)
        ch_trimgalore = TRIMGALORE(trim_branches.trimgalore)
        GENERATE_REPORT(ch_samples)
    }
    ```

Diese Closure wird ausgeführt, wenn der Workflow abgeschlossen ist. Darin hast du Zugriff auf das `workflow`-Objekt, das nützliche Eigenschaften über die Ausführung bereitstellt.

Führe deinen Workflow aus und du wirst diese Zusammenfassung am Ende sehen!

```bash
nextflow run main.nf --input ./data/samples.csv -ansi-log false
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W  ~  version 25.10.2
    Launching `main.nf` [marvelous_boltzmann] DSL2 - revision: a31662a7c1
    WARN: Low sequencing depth for sample_002: 25000000
    [9b/d48e40] Submitted process > FASTP (2)
    [6a/73867a] Submitted process > GENERATE_REPORT (2)
    [79/ad0ac5] Submitted process > GENERATE_REPORT (1)
    [f3/bda6cb] Submitted process > FASTP (1)
    [34/d5b52f] Submitted process > GENERATE_REPORT (3)

    Pipeline execution summary:
    ==========================
    Completed at: 2025-10-10T12:14:24.885384+01:00
    Duration    : 2.9s
    Success     : true
    workDir     : /workspaces/training/side-quests/essential_scripting_patterns/work
    exit status : 0
    ```

Lass uns es nützlicher machen, indem wir bedingte Logik hinzufügen:

=== "Danach"

    ```groovy title="main.nf" linenums="66" hl_lines="5-22"
        ch_fastp = FASTP(trim_branches.fastp)
        ch_trimgalore = TRIMGALORE(trim_branches.trimgalore)
        GENERATE_REPORT(ch_samples)

        workflow.onComplete = {
            println ""
            println "Pipeline execution summary:"
            println "=========================="
            println "Completed at: ${workflow.complete}"
            println "Duration    : ${workflow.duration}"
            println "Success     : ${workflow.success}"
            println "workDir     : ${workflow.workDir}"
            println "exit status : ${workflow.exitStatus}"
            println ""

            if (workflow.success) {
                println "✅ Pipeline completed successfully!"
            } else {
                println "❌ Pipeline failed!"
                println "Error: ${workflow.errorMessage}"
            }
        }
    }
    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="66" hl_lines="5-16"
        ch_fastp = FASTP(trim_branches.fastp)
        ch_trimgalore = TRIMGALORE(trim_branches.trimgalore)
        GENERATE_REPORT(ch_samples)

        workflow.onComplete = {
            println ""
            println "Pipeline execution summary:"
            println "=========================="
            println "Completed at: ${workflow.complete}"
            println "Duration    : ${workflow.duration}"
            println "Success     : ${workflow.success}"
            println "workDir     : ${workflow.workDir}"
            println "exit status : ${workflow.exitStatus}"
            println ""
        }
    }
    ```

Jetzt erhalten wir eine noch informativere Zusammenfassung, einschließlich einer Erfolgs-/Fehlermeldung und des Ausgabeverzeichnisses, falls angegeben:

<!-- TODO: add run command -->

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W  ~  version 25.10.2
    Launching `main.nf` [boring_linnaeus] DSL2 - revision: a31662a7c1
    WARN: Low sequencing depth for sample_002: 25000000
    [e5/242efc] Submitted process > FASTP (2)
    [3b/74047c] Submitted process > GENERATE_REPORT (3)
    [8a/7a57e6] Submitted process > GENERATE_REPORT (1)
    [a8/b1a31f] Submitted process > GENERATE_REPORT (2)
    [40/648429] Submitted process > FASTP (1)

    Pipeline execution summary:
    ==========================
    Completed at: 2025-10-10T12:16:00.522569+01:00
    Duration    : 3.6s
    Success     : true
    workDir     : /workspaces/training/side-quests/essential_scripting_patterns/work
    exit status : 0

    ✅ Pipeline completed successfully!
    ```

Du kannst die Zusammenfassung auch mit Dateioperationen in eine Datei schreiben:

```groovy title="main.nf - Writing summary to file"
workflow {
    // ... dein Workflow-Code ...

    workflow.onComplete = {
        def summary = """
        Pipeline Execution Summary
        ===========================
        Completed: ${workflow.complete}
        Duration : ${workflow.duration}
        Success  : ${workflow.success}
        Command  : ${workflow.commandLine}
        """

        println summary

        // In Log-Datei schreiben
        def log_file = file("${workflow.launchDir}/pipeline_summary.txt")
        log_file.text = summary
    }
}
```

### 8.2. Der `onError`-Handler

Neben `onComplete` gibt es noch einen weiteren Event-Handler, den du verwenden kannst: `onError`, der nur ausgeführt wird, wenn der Workflow fehlschlägt:

```groovy title="main.nf - onError handler"
workflow {
    // ... dein Workflow-Code ...

    workflow.onError = {
        println "="* 50
        println "Pipeline execution failed!"
        println "Error message: ${workflow.errorMessage}"
        println "="* 50

        // Detailliertes Fehlerprotokoll schreiben
        def error_file = file("${workflow.launchDir}/error.log")
        error_file.text = """
        Workflow Error Report
        =====================
        Time: ${new Date()}
        Error: ${workflow.errorMessage}
        Error report: ${workflow.errorReport ?: 'No detailed report available'}
        """

        println "Error details written to: ${error_file}"
    }
}
```

Du kannst mehrere Handler zusammen in deinem Workflow-Skript verwenden:

```groovy title="main.nf - Combined handlers"
workflow {
    // ... dein Workflow-Code ...

    workflow.onError = {
        println "Workflow failed: ${workflow.errorMessage}"
    }

    workflow.onComplete = {
        def duration_mins = workflow.duration.toMinutes().round(2)
        def status = workflow.success ? "SUCCESS ✅" : "FAILED ❌"

        println """
        Pipeline finished: ${status}
        Duration: ${duration_mins} minutes
        """
    }
}
```

### Fazit

In diesem Abschnitt hast du gelernt:

- **Event-Handler-Closures**: Closures in deinem Workflow-Skript, die an verschiedenen Lebenszyklus-Punkten ausgeführt werden
- **`onComplete`-Handler**: Für Ausführungszusammenfassungen und Ergebnis-Reporting
- **`onError`-Handler**: Für Fehlerbehandlung und Protokollierung von Fehlern
- **Workflow-Objekt-Eigenschaften**: Zugriff auf `workflow.success`, `workflow.duration`, `workflow.errorMessage`, etc.

Event-Handler zeigen, wie du die volle Kraft der Nextflow-Sprache innerhalb deiner Workflow-Skripte nutzen kannst, um anspruchsvolle Logging- und Benachrichtigungsfähigkeiten hinzuzufügen.

---

## Zusammenfassung

Glückwunsch, du hast es geschafft!

In diesem Side Quest hast du eine umfassende Probenverarbeitungs-Pipeline erstellt, die sich von grundlegender Metadaten-Handhabung zu einem anspruchsvollen, produktionsreifen Workflow entwickelt hat.
Jeder Abschnitt baute auf dem vorherigen auf und demonstrierte, wie Programmierkonzepte einfache Workflows in leistungsstarke Datenverarbeitungssysteme verwandeln, mit folgenden Vorteilen:

- **Klarerer Code**: Das Verstehen von Datenfluss vs. Scripting hilft dir, organisiertere Workflows zu schreiben
- **Robuste Handhabung**: Sichere Navigation und Elvis-Operatoren machen Workflows resilient gegen fehlende Daten
- **Flexible Verarbeitung**: Bedingte Logik ermöglicht es deinen Workflows, verschiedene Probentypen angemessen zu verarbeiten
- **Adaptive Ressourcen**: Dynamische Direktiven optimieren die Ressourcennutzung basierend auf Eingabeeigenschaften

Diese Progression spiegelt die reale Evolution von Bioinformatik-Pipelines wider, von Forschungsprototypen, die einige Proben handhaben, zu Produktionssystemen, die Tausende von Proben über Labore und Institutionen hinweg verarbeiten.
Jede Herausforderung, die du gelöst hast, und jedes Muster, das du gelernt hast, spiegelt tatsächliche Probleme wider, denen Entwickler, denen Entwickler\*innen beim Skalieren von Nextflow-Workflows begegnen.
Die Anwendung dieser Muster in deiner eigenen Arbeit wird es dir ermöglichen, robuste, produktionsreife Workflows zu erstellen.

### Wichtige Muster

1.  **Datenfluss vs. Scripting:** Du hast gelernt, zwischen Datenfluss-Operationen (Kanal-Orchestrierung) und Scripting (Code, der Daten manipuliert) zu unterscheiden, einschließlich der entscheidenden Unterschiede zwischen Operationen auf verschiedenen Typen wie `collect` auf Channel vs. List.

    - Datenfluss: Kanal-Orchestrierung

    ```groovy
    channel.fromPath('*.fastq').splitCsv(header: true)
    ```

    - Scripting: Datenverarbeitung auf Collections

    ```groovy
    sample_data.collect { it.toUpperCase() }
    ```

2.  **Fortgeschrittene String-Verarbeitung**: Du hast reguläre Ausdrücke zum Parsen von Dateinamen, dynamische Skript-Generierung in Prozessen und Variableninterpolation (Nextflow vs. Bash vs. Shell) gemeistert.

    - Musterabgleich

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

    - Datei-Collection zu Befehlsargumenten (im Prozess-Skript-Block)

    ```groovy
    script:
    def file_args = input_files.collect { file -> "--input ${file}" }.join(' ')
    """
    analysis_tool ${file_args} --output results.txt
    """
    ```

3.  **Wiederverwendbare Funktionen erstellen**: Du hast gelernt, komplexe Logik in benannte Funktionen zu extrahieren, die von Kanal-Operatoren aufgerufen werden können, wodurch Workflows lesbarer und wartbarer werden.

    - Eine benannte Funktion definieren

    ```groovy
    def separateMetadata(row) {
        def sample_meta = [ /* code hidden for brevity */ ]
        def fastq_path = file(row.file_path)
        def m = (fastq_path.name =~ /^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$/)
        def file_meta = m ? [ /* code hidden for brevity */ ] : [:]
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

4.  **Dynamische Ressourcen-Direktiven mit Closures**: Du hast die Verwendung von Closures in Prozess-Direktiven für adaptive Ressourcenzuweisung basierend auf Eingabeeigenschaften erkundet.

    - Benannte Closures und Komposition

    ```groovy
    def enrichData = normalizeId >> addQualityCategory >> addFlags
    def processor = generalFunction.curry(fixedParam)
    ```

    - Closures mit Bereichszugriff

    ```groovy
    def collectStats = { data -> stats.count++; return data }
    ```

5.  **Bedingte Logik und Prozesskontrolle**: Du hast intelligentes Routing mit `.branch()`- und `.filter()`-Operatoren hinzugefügt und Truthiness für prägnante bedingte Ausdrücke genutzt.

    - Verwende `.branch()`, um Daten durch verschiedene Workflow-Zweige zu routen

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
    if (sample.files) println "Has files"
    ```

    - Verwende `filter()`, um Daten mit 'Truthiness' zu unterteilen

    ```groovy
    ch_valid_samples = ch_samples
        .filter { meta, reads ->
            meta.id && meta.organism && meta.depth >= 25000000
        }
    ```

6.  **Sichere Navigation und Elvis-Operatoren**: Du hast die Pipeline robust gegen fehlende Daten gemacht, indem du `?.` für null-sicheren Eigenschaftszugriff und `?:` zum Bereitstellen von Standardwerten verwendet hast.

    ```groovy
    def id = data?.sample?.id ?: 'unknown'
    ```

7.  **Validierung mit error() und log.warn**: Du hast gelernt, Eingaben früh zu validieren und schnell mit klaren Fehlermeldungen zu scheitern.

    ```groovy
    try {
        def errors = validateSample(sample)
        if (errors) throw new RuntimeException("Invalid: ${errors.join(', ')}")
    } catch (Exception e) {
        println "Error: ${e.message}"
    }
    ```

8.  **Konfigurations-Event-Handler**: Du hast gelernt, Workflow-Event-Handler (`onComplete` und `onError`) für Logging, Benachrichtigungen und Lifecycle-Management zu verwenden.

    - Verwendung von `onComplete` zum Protokollieren und Benachrichtigen

    ```groovy
    workflow.onComplete = {
        println "Success     : ${workflow.success}"
        println "exit status : ${workflow.exitStatus}"

        if (workflow.success) {
            println "✅ Pipeline completed successfully!"
        } else {
            println "❌ Pipeline failed!"
            println "Error: ${workflow.errorMessage}"
        }
    }
    ```

    - Verwendung von `onError`, um speziell im Fehlerfall Maßnahmen zu ergreifen

    ```groovy
    workflow.onError = {
        // Detailliertes Fehlerprotokoll schreiben
        def error_file = file("${workflow.launchDir}/error.log")
        error_file.text = """
        Time: ${new Date()}
        Error: ${workflow.errorMessage}
        Error report: ${workflow.errorReport ?: 'No detailed report available'}
        """

        println "Error details written to: ${error_file}"
    }
    ```

### Zusätzliche Ressourcen

- [Nextflow Language Reference](https://nextflow.io/docs/latest/reference/syntax.html)
- [Nextflow Operators](https://www.nextflow.io/docs/latest/operator.html)
- [Nextflow Script Syntax](https://www.nextflow.io/docs/latest/script.html)
- [Nextflow Standard Library](https://nextflow.io/docs/latest/reference/stdlib.html)

Schau dir diese Ressourcen unbedingt an, wenn du fortgeschrittenere Features erkunden musst.

Du wirst davon profitieren, zu üben und deine Fähigkeiten zu erweitern, um:

- Sauberere Workflows mit ordnungsgemäßer Trennung zwischen Datenfluss und Scripting zu schreiben
- Variableninterpolation zu beherrschen, um häufige Fallstricke mit Nextflow-, Bash- und Shell-Variablen zu vermeiden
- Dynamische Ressourcen-Direktiven für effiziente, adaptive Workflows zu verwenden
- Datei-Collections in korrekt formatierte Befehlszeilenargumente zu transformieren
- Verschiedene Dateibenennungskonventionen und Eingabeformate elegant mit Regex und String-Verarbeitung zu handhaben
- Wiederverwendbaren, wartbaren Code mit fortgeschrittenen Closure-Mustern und funktionaler Programmierung zu erstellen
- Komplexe Datensätze mit Collection-Operationen zu verarbeiten und zu organisieren
- Validierung, Fehlerbehandlung und Logging hinzuzufügen, um deine Workflows produktionsreif zu machen
- Workflow-Lifecycle-Management mit Event-Handlern zu implementieren

---

## Wie geht es weiter?

Kehre zum [Menü der Side Quests](./index.md) zurück oder klicke auf den Button unten rechts auf der Seite, um zum nächsten Thema in der Liste zu gelangen.
