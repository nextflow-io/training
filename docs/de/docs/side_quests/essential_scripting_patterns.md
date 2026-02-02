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
- **Bedingtes Routing:** Route Proben intelligent durch verschiedene Prozesse basierend auf ihren Metadaten-Eigenschaften
- **Sichere Operationen:** Gehe mit fehlenden Daten elegant um mit Null-Safe-Operatoren und validiere Eingaben mit klaren Fehlermeldungen
- **Konfigurationsbasierte Handler:** Verwende Workflow-Event-Handler für Logging, Benachrichtigungen und Lifecycle-Management

### Voraussetzungen

Bevor du dieses Side Quest angehst, solltest du:

- Das [Hello Nextflow](../hello_nextflow/README.md) Tutorial oder einen gleichwertigen Einsteigerkurs abgeschlossen haben.
- Dich wohl fühlen mit grundlegenden Nextflow-Konzepten und -Mechanismen (Prozesse, Channels, Operatoren, Arbeiten mit Dateien, Metadaten)
- Grundlegende Vertrautheit mit gängigen Programmierkonstrukten haben (Variablen, Maps, Listen)

Dieses Tutorial wird Programmierkonzepte erklären, sobald wir auf sie stoßen, du brauchst also keine umfangreiche Programmiererfahrung.
Wir beginnen mit grundlegenden Konzepten und bauen zu fortgeschrittenen Mustern auf.

---

## 0. Los geht's

#### Öffne den Training-Codespace

Falls noch nicht geschehen, stelle sicher, dass du die Trainingsumgebung wie in [Environment Setup](../envsetup/index.md) beschrieben öffnest.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Wechsle in das Projektverzeichnis

Lass uns in das Verzeichnis wechseln, in dem sich die Dateien für dieses Tutorial befinden.

```bash
cd side-quests/essential_scripting_patterns
```

#### Überprüfe die Materialien

Du findest eine Haupt-Workflow-Datei und ein `data`-Verzeichnis mit Beispieldateien.

```console title="Verzeichnisinhalt"
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

Wir werden diesen realistischen Datensatz verwenden, um praktische Programmiertechniken zu erkunden, denen du in echten Bioinformatik-Workflows begegnen wirst.

<!-- TODO: Can we make this more domain-agnostic? -->

<!-- TODO: add an assignment statement? #### Überprüfe die Aufgabe -->

#### Bereitschafts-Checkliste

Denkst du, du bist bereit einzutauchen?

- [ ] Ich verstehe das Ziel dieses Kurses und seine Voraussetzungen
- [ ] Mein Codespace ist betriebsbereit
- [ ] Ich habe mein Arbeitsverzeichnis entsprechend gesetzt
<!-- - [ ] Ich verstehe die Aufgabe -->

Wenn du alle Kästchen ankreuzen kannst, kann's losgehen.

---

## 1. Dataflow vs Scripting: Die Grenzen verstehen

### 1.1. Identifizieren, was was ist

Beim Schreiben von Nextflow-Workflows ist es wichtig, zwischen **Dataflow** (wie Daten durch Channels und Prozesse fließen) und **Scripting** (der Code, der Daten manipuliert und Entscheidungen trifft) zu unterscheiden. Lass uns einen Workflow erstellen, der zeigt, wie sie zusammenarbeiten.

#### 1.1.1. Grundlegender Nextflow-Workflow

Beginne mit einem einfachen Workflow, der nur die CSV-Datei liest (wir haben das bereits für dich in `main.nf` erledigt):

```groovy title="main.nf" linenums="1"
workflow {
    ch_samples = channel.fromPath("./data/samples.csv")
        .splitCsv(header: true)
        .view()
}
```

Der `workflow`-Block definiert unsere Pipeline-Struktur, während `channel.fromPath()` einen Channel aus einem Dateipfad erstellt. Der `.splitCsv()`-Operator verarbeitet die CSV-Datei und konvertiert jede Zeile in eine Map-Datenstruktur.

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

#### 1.1.2. Den Map-Operator hinzufügen

Jetzt werden wir Scripting hinzufügen, um die Daten zu transformieren, mit dem `.map()`-Operator, den du wahrscheinlich bereits kennst. Dieser Operator nimmt eine 'Closure', in der wir Code schreiben können, um jedes Element zu transformieren.

!!! note

    Eine **Closure** ist ein Codeblock, der herumgereicht und später ausgeführt werden kann. Denke daran wie an eine Funktion, die du inline definierst. Closures werden mit geschweiften Klammern `{ }` geschrieben und können Parameter annehmen. Sie sind grundlegend dafür, wie Nextflow-Operatoren funktionieren, und wenn du schon eine Weile Nextflow schreibst, hast du sie vielleicht bereits verwendet, ohne es zu merken!

So sieht diese Map-Operation aus:

=== "Nachher"

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

Das ist unsere erste **Closure** – eine anonyme Funktion, die du als Argument übergeben kannst (ähnlich wie Lambdas in Python oder Arrow Functions in JavaScript). Closures sind essentiell für die Arbeit mit Nextflow-Operatoren.

Die Closure `{ row -> return row }` nimmt einen Parameter `row` (könnte auch anders heißen: `item`, `sample`, etc.).

Wenn der `.map()`-Operator jedes Channel-Element verarbeitet, übergibt er dieses Element an deine Closure. Hier enthält `row` jeweils eine CSV-Zeile.

Wende diese Änderung an und führe den Workflow aus:

```bash
nextflow run main.nf
```

Du siehst die gleiche Ausgabe wie zuvor, weil wir einfach die Eingabe unverändert zurückgeben. Das bestätigt, dass der Map-Operator korrekt funktioniert. Jetzt fangen wir an, die Daten zu transformieren.

#### 1.1.3. Eine Map-Datenstruktur erstellen

Jetzt schreiben wir **Scripting**-Logik innerhalb unserer Closure, um jede Datenzeile zu transformieren. Hier verarbeiten wir einzelne Datenelemente, anstatt den Datenfluss zu orchestrieren.

=== "Nachher"

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

Die `sample_meta`-Map ist eine Schlüssel-Wert-Datenstruktur (wie Dictionaries in Python, Objekte in JavaScript oder Hashes in Ruby), die verwandte Informationen speichert: Proben-ID, Organismus, Gewebetyp, Sequenzierungstiefe und Qualitätswert.

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

#### 1.1.4. Bedingte Logik hinzufügen

Jetzt fügen wir mehr Scripting hinzu – diesmal mit einem ternären Operator, um Entscheidungen basierend auf Datenwerten zu treffen.

Nimm die folgende Änderung vor:

=== "Nachher"

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

Der ternäre Operator ist eine Kurzform für ein if/else-Statement, das dem Muster `Bedingung ? Wert_wenn_wahr : Wert_wenn_falsch` folgt. Diese Zeile bedeutet: „Wenn die Qualität größer als 40 ist, verwende 'high', ansonsten verwende 'normal'". Sein Cousin, der **Elvis-Operator** (`?:`), liefert Standardwerte, wenn etwas null oder leer ist – wir werden dieses Muster später in diesem Tutorial erkunden.

Der Map-Additions-Operator `+` erstellt eine **neue Map**, anstatt die bestehende zu modifizieren. Diese Zeile erstellt eine neue Map, die alle Schlüssel-Wert-Paare von `sample_meta` plus den neuen `priority`-Schlüssel enthält.

!!! Note

    Modifiziere niemals Maps, die in Closures übergeben werden – erstelle immer neue mit `+` (zum Beispiel). In Nextflow fließen die gleichen Daten oft gleichzeitig durch mehrere Operationen. Das In-Place-Modifizieren einer Map kann unvorhersehbare Nebeneffekte verursachen, wenn andere Operationen auf dasselbe Objekt verweisen. Das Erstellen neuer Maps stellt sicher, dass jede Operation ihre eigene saubere Kopie hat.

Führe den geänderten Workflow aus:

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

#### 1.1.5. Maps mit `.subMap()` aufteilen

Während der `+`-Operator Schlüssel zu einer Map hinzufügt, musst du manchmal das Gegenteil tun – nur bestimmte Schlüssel extrahieren. Die `.subMap()`-Methode ist dafür perfekt.

Fügen wir eine Zeile hinzu, um eine vereinfachte Version unserer Metadaten zu erstellen, die nur Identifikationsfelder enthält:

=== "Nachher"

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
                println "Nur ID-Felder: ${id_only}"

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

Führe den geänderten Workflow aus:

```bash
nextflow run main.nf
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [peaceful_cori] DSL2 - revision: 4cc4a8340f

    Nur ID-Felder: [id:sample_001, organism:human, tissue:liver]
    Nur ID-Felder: [id:sample_002, organism:mouse, tissue:brain]
    Nur ID-Felder: [id:sample_003, organism:human, tissue:kidney]
    [id:sample_001, organism:human, tissue:liver, depth:30000000, quality:38.5, priority:normal]
    [id:sample_002, organism:mouse, tissue:brain, depth:25000000, quality:35.2, priority:normal]
    [id:sample_003, organism:human, tissue:kidney, depth:45000000, quality:42.1, priority:high]
    ```

Dies zeigt sowohl die vollständigen Metadaten, die durch die `view()`-Operation angezeigt werden, als auch die extrahierte Teilmenge, die wir mit `println` ausgegeben haben.

Die `.subMap()`-Methode nimmt eine Liste von Schlüsseln und gibt eine neue Map zurück, die nur diese Schlüssel enthält. Wenn ein Schlüssel in der ursprünglichen Map nicht existiert, wird er einfach nicht im Ergebnis enthalten.

Dies ist besonders nützlich, wenn du verschiedene Metadaten-Versionen für verschiedene Prozesse erstellen musst – manche benötigen möglicherweise vollständige Metadaten, während andere nur minimale Identifikationsfelder brauchen.

Entferne jetzt diese println-Anweisungen, um deinen Workflow in seinen vorherigen Zustand zurückzuversetzen, da wir sie nicht mehr brauchen.

!!! tip "Zusammenfassung Map-Operationen"

    - **Schlüssel hinzufügen**: `map1 + [new_key: value]` - Erstellt neue Map mit zusätzlichen Schlüsseln
    - **Schlüssel extrahieren**: `map1.subMap(['key1', 'key2'])` - Erstellt neue Map mit nur angegebenen Schlüsseln
    - **Beide Operationen erstellen neue Maps** - Ursprüngliche Maps bleiben unverändert

#### 1.1.6. Maps kombinieren und Ergebnisse zurückgeben

Bisher haben wir nur das zurückgegeben, was die Nextflow-Community die 'Meta-Map' nennt, und wir haben die Dateien ignoriert, auf die sich diese Metadaten beziehen. Aber wenn du Nextflow-Workflows schreibst, möchtest du wahrscheinlich etwas mit diesen Dateien machen.

Lass uns eine Channel-Struktur ausgeben, die aus einem Tupel von 2 Elementen besteht: der angereicherten Metadaten-Map und dem entsprechenden Dateipfad. Dies ist ein gängiges Muster in Nextflow zum Übergeben von Daten an Prozesse.

=== "Nachher"

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

!!! note

    **Maps und Metadaten**: Maps sind grundlegend für die Arbeit mit Metadaten in Nextflow. Für eine detailliertere Erklärung zum Arbeiten mit Metadaten-Maps siehe das [Working with metadata](./metadata.md) Side Quest.

Unser Workflow demonstriert das Kernmuster: **Dataflow-Operationen** (`workflow`, `channel.fromPath()`, `.splitCsv()`, `.map()`, `.view()`) orchestrieren, wie Daten durch die Pipeline fließen, während **Scripting** (Maps `[key: value]`, String-Methoden, Typkonvertierungen, ternäre Operatoren) innerhalb der `.map()`-Closure die Transformation einzelner Datenelemente übernimmt.

### 1.2. Verschiedene Typen verstehen: Channel vs List

Bis jetzt gut, wir können zwischen Dataflow-Operationen und Scripting unterscheiden. Aber was ist, wenn derselbe Methodenname in beiden Kontexten existiert?

Ein perfektes Beispiel ist die `collect`-Methode, die sowohl für Channel-Typen als auch für List-Typen in der Nextflow-Standardbibliothek existiert. Die `collect()`-Methode auf einer List transformiert jedes Element, während der `collect()`-Operator auf einem Channel alle Channel-Emissionen in einem Single-Item-Channel sammelt.

Lass uns dies mit einigen Beispieldaten demonstrieren, beginnend mit einer Auffrischung, was der Channel-`collect()`-Operator macht. Schau dir `collect.nf` an:

```groovy title="collect.nf" linenums="1"
def sample_ids = ['sample_001', 'sample_002', 'sample_003']

// channel.collect() - gruppiert mehrere Channel-Emissionen zu einer
ch_input = channel.fromList(sample_ids)
ch_input.view { sample -> "Einzelnes Channel-Element: ${sample}" }
ch_collected = ch_input.collect()
ch_collected.view { list -> "channel.collect() Ergebnis: ${list} (${list.size()} Elemente zu 1 gruppiert)" }
```

Schritte:

- Definiere eine List von Proben-IDs
- Erstelle einen Channel mit `fromList()`, der jede Proben-ID separat emittiert
- Gib jedes Element mit `view()` aus, während es durchfließt
- Sammle alle Elemente in eine einzelne Liste mit dem Channel-`collect()`-Operator
- Gib das gesammelte Ergebnis (einzelnes Element, das alle Proben-IDs enthält) mit einem zweiten `view()` aus

Wir haben die Struktur des Channels geändert, aber wir haben die Daten selbst nicht geändert.

Führe den Workflow aus, um dies zu bestätigen:

```bash
nextflow run collect.nf
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `collect.nf` [loving_mendel] DSL2 - revision: e8d054a46e

    Einzelnes Channel-Element: sample_001
    Einzelnes Channel-Element: sample_002
    Einzelnes Channel-Element: sample_003
    channel.collect() Ergebnis: [sample_001, sample_002, sample_003] (3 Elemente zu 1 gruppiert)
    ```

`view()` gibt eine Ausgabe für jede Channel-Emission zurück, also wissen wir, dass diese einzelne Ausgabe alle 3 ursprünglichen Elemente gruppiert in eine Liste enthält.

Jetzt lass uns die `collect`-Methode auf einer List in Aktion sehen. Modifiziere `collect.nf`, um die List-`collect`-Methode auf die ursprüngliche Liste der Proben-IDs anzuwenden:

=== "Nachher"

    ```groovy title="main.nf" linenums="1" hl_lines="9-13"
    def sample_ids = ['sample_001', 'sample_002', 'sample_003']

    // channel.collect() - gruppiert mehrere Channel-Emissionen zu einer
    ch_input = channel.fromList(sample_ids)
    ch_input.view { sample -> "Einzelnes Channel-Element: ${sample}" }
    ch_collected = ch_input.collect()
    ch_collected.view { list -> "channel.collect() Ergebnis: ${list} (${list.size()} Elemente zu 1 gruppiert)" }

    // List.collect() - transformiert jedes Element, erhält Struktur
    def formatted_ids = sample_ids.collect { id ->
        id.toUpperCase().replace('SAMPLE_', 'SPECIMEN_')
    }
    println "List.collect() Ergebnis: ${formatted_ids} (${sample_ids.size()} Elemente in ${formatted_ids.size()} transformiert)"
    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="1"
    def sample_ids = ['sample_001', 'sample_002', 'sample_003']

    // channel.collect() - gruppiert mehrere Channel-Emissionen zu einer
    ch_input = channel.fromList(sample_ids)
    ch_input.view { sample -> "Einzelnes Channel-Element: ${sample}" }
    ch_collected = ch_input.collect()
    ch_collected.view { list -> "channel.collect() Ergebnis: ${list} (${list.size()} Elemente zu 1 gruppiert)" }
    ```

In diesem neuen Snippet:

- Definieren wir eine neue Variable `formatted_ids`, die die List-`collect`-Methode verwendet, um jede Proben-ID in der ursprünglichen Liste zu transformieren
- Geben das Ergebnis mit `println` aus

Führe den geänderten Workflow aus:

```bash
nextflow run collect.nf
```

??? success "Befehlsausgabe"

    ```console hl_lines="5"
    N E X T F L O W   ~  version 25.10.2

    Launching `collect.nf` [cheeky_stonebraker] DSL2 - revision: 2d5039fb47

    List.collect() Ergebnis: [SPECIMEN_001, SPECIMEN_002, SPECIMEN_003] (3 Elemente in 3 transformiert)
    Einzelnes Channel-Element: sample_001
    Einzelnes Channel-Element: sample_002
    Einzelnes Channel-Element: sample_003
    channel.collect() Ergebnis: [sample_001, sample_002, sample_003] (3 Elemente zu 1 gruppiert)
    ```

Diesmal haben wir NICHT die Struktur der Daten geändert, wir haben immer noch 3 Elemente in der Liste, aber wir HABEN jedes Element mit der List-`collect`-Methode transformiert, um eine neue Liste mit modifizierten Werten zu erzeugen. Dies ist ähnlich wie die Verwendung des `map`-Operators auf einem Channel, aber es arbeitet auf einer List-Datenstruktur anstatt einem Channel.

`collect` ist ein Extremfall, den wir hier verwenden, um einen Punkt zu machen. Die Kernlektion ist, dass du beim Schreiben von Workflows immer zwischen **Datenstrukturen** (Lists, Maps, etc.) und **Channels** (Dataflow-Konstrukten) unterscheiden solltest. Operationen können Namen teilen, verhalten sich aber völlig unterschiedlich, je nachdem, auf welchem Typ sie aufgerufen werden.

### 1.3. Der Spread-Operator (`*.`) - Kurzform für Eigenschaftsextraktion

Verwandt mit der List-`collect`-Methode ist der Spread-Operator (`*.`), der eine prägnante Möglichkeit bietet, Eigenschaften aus Collections zu extrahieren. Es ist im Wesentlichen syntaktischer Zucker für ein gängiges `collect`-Muster.

Fügen wir eine Demonstration zu unserer `collect.nf`-Datei hinzu:

=== "Nachher"

    ```groovy title="collect.nf" linenums="1" hl_lines="15-18"
    def sample_ids = ['sample_001', 'sample_002', 'sample_003']

    // channel.collect() - gruppiert mehrere Channel-Emissionen zu einer
    ch_input = channel.fromList(sample_ids)
    ch_input.view { sample -> "Einzelnes Channel-Element: ${sample}" }
    ch_collected = ch_input.collect()
    ch_collected.view { list -> "channel.collect() Ergebnis: ${list} (${list.size()} Elemente zu 1 gruppiert)" }

    // List.collect() - transformiert jedes Element, erhält Struktur
    def formatted_ids = sample_ids.collect { id ->
        id.toUpperCase().replace('SAMPLE_', 'SPECIMEN_')
    }
    println "List.collect() Ergebnis: ${formatted_ids} (${sample_ids.size()} Elemente in ${formatted_ids.size()} transformiert)"

    // Spread-Operator - prägnanter Eigenschaftszugriff
    def sample_data = [[id: 's1', quality: 38.5], [id: 's2', quality: 42.1], [id: 's3', quality: 35.2]]
    def all_ids = sample_data*.id
    println "Spread-Operator Ergebnis: ${all_ids}"
    ```

=== "Vorher"

    ```groovy title="collect.nf" linenums="1"
    def sample_ids = ['sample_001', 'sample_002', 'sample_003']

    // channel.collect() - gruppiert mehrere Channel-Emissionen zu einer
    ch_input = channel.fromList(sample_ids)
    ch_input.view { sample -> "Einzelnes Channel-Element: ${sample}" }
    ch_collected = ch_input.collect()
    ch_collected.view { list -> "channel.collect() Ergebnis: ${list} (${list.size()} Elemente zu 1 gruppiert)" }

    // List.collect() - transformiert jedes Element, erhält Struktur
    def formatted_ids = sample_ids.collect { id ->
        id.toUpperCase().replace('SAMPLE_', 'SPECIMEN_')
    }
    println "List.collect() Ergebnis: ${formatted_ids} (${sample_ids.size()} Elemente in ${formatted_ids.size()} transformiert)"
    ```

Führe den aktualisierten Workflow aus:

```bash title="Spread-Operator testen"
nextflow run collect.nf
```

??? success "Befehlsausgabe"

    ```console hl_lines="6"
    N E X T F L O W   ~  version 25.10.2

    Launching `collect.nf` [cranky_galileo] DSL2 - revision: 5f3c8b2a91

    List.collect() Ergebnis: [SPECIMEN_001, SPECIMEN_002, SPECIMEN_003] (3 Elemente in 3 transformiert)
    Spread-Operator Ergebnis: [s1, s2, s3]
    Einzelnes Channel-Element: sample_001
    Einzelnes Channel-Element: sample_002
    Einzelnes Channel-Element: sample_003
    channel.collect() Ergebnis: [sample_001, sample_002, sample_003] (3 Elemente zu 1 gruppiert)
    ```

Der Spread-Operator `*.` ist eine Kurzform für ein gängiges Collect-Muster:

```groovy
// Diese sind äquivalent:
def ids = samples*.id
def ids = samples.collect { it.id }

// Funktioniert auch mit Methodenaufrufen:
def names = files*.getName()
def names = files.collect { it.getName() }
```

Der Spread-Operator ist besonders nützlich, wenn du eine einzelne Eigenschaft aus einer Liste von Objekten extrahieren musst – er ist lesbarer als das Ausschreiben der vollständigen `collect`-Closure.

!!! tip "Wann Spread vs Collect verwenden"

    - **Verwende Spread (`*.`)** für einfachen Eigenschaftszugriff: `samples*.id`, `files*.name`
    - **Verwende collect** für Transformationen oder komplexe Logik: `samples.collect { it.id.toUpperCase() }`, `samples.collect { [it.id, it.quality > 40] }`

### Fazit

In diesem Abschnitt hast du gelernt:

- **Dataflow vs Scripting**: Channel-Operatoren orchestrieren, wie Daten durch deine Pipeline fließen, während Scripting einzelne Datenelemente transformiert
- **Typen verstehen**: Derselbe Methodenname (wie `collect`) kann sich je nach Typ, auf dem er aufgerufen wird, unterschiedlich verhalten (Channel vs List)
- **Kontext ist wichtig**: Sei dir immer bewusst, ob du mit Channels (Dataflow) oder Datenstrukturen (Scripting) arbeitest

Das Verstehen dieser Grenzen ist essentiell für Debugging, Dokumentation und das Schreiben wartbarer Workflows.

Als Nächstes tauchen wir tiefer in String-Verarbeitungsfähigkeiten ein, die essentiell für den Umgang mit realen Daten sind.

---

## 2. String-Verarbeitung und dynamische Script-Generierung

Die Beherrschung der String-Verarbeitung unterscheidet brüchige Workflows von robusten Pipelines. Dieser Abschnitt behandelt das Parsen komplexer Dateinamen, dynamische Script-Generierung und Variable Interpolation.

### 2.1. Pattern Matching und reguläre Ausdrücke

Bioinformatik-Dateien haben oft komplexe Namenskonventionen, die Metadaten kodieren. Lass uns dies automatisch mit Pattern Matching und regulären Ausdrücken extrahieren.

Wir kehren zu unserem `main.nf`-Workflow zurück und fügen etwas Pattern-Matching-Logik hinzu, um zusätzliche Probeninformationen aus Dateinamen zu extrahieren. Die FASTQ-Dateien in unserem Datensatz folgen Illumina-Namenskonventionen mit Namen wie `SAMPLE_001_S1_L001_R1_001.fastq.gz`. Diese mögen kryptisch aussehen, aber sie kodieren tatsächlich nützliche Metadaten wie Proben-ID, Lane-Nummer und Read-Richtung. Wir werden Regex-Fähigkeiten verwenden, um diese Namen zu parsen.

Nimm die folgende Änderung an deinem bestehenden `main.nf`-Workflow vor:

=== "Nachher"

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

1. **Reguläre Ausdrucks-Literale** mit `~/pattern/`-Syntax - dies erstellt ein Regex-Muster ohne die Notwendigkeit, Backslashes zu escapen
2. **Pattern Matching** mit dem `=~`-Operator - dies versucht, einen String gegen ein Regex-Muster zu matchen
3. **Matcher-Objekte**, die Gruppen mit `[0][1]`, `[0][2]`, etc. erfassen - `[0]` bezieht sich auf den gesamten Match, `[1]`, `[2]`, etc. beziehen sich auf erfasste Gruppen in Klammern

Lass uns das Regex-Muster `^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$` aufschlüsseln:

| Muster              | Passt                                      | Erfasst                                |
| ------------------- | ------------------------------------------ | -------------------------------------- |
| `^(.+)`             | Probenname vom Anfang                      | Gruppe 1: Probenname                   |
| `_S(\d+)`           | Probennummer `_S1`, `_S2`, etc.            | Gruppe 2: Probennummer                 |
| `_L(\d{3})`         | Lane-Nummer `_L001`                        | Gruppe 3: Lane (3 Ziffern)             |
| `_(R[12])`          | Read-Richtung `_R1` oder `_R2`             | Gruppe 4: Read-Richtung                |
| `_(\d{3})`          | Chunk-Nummer `_001`                        | Gruppe 5: Chunk (3 Ziffern)            |
| `\.fastq(?:\.gz)?$` | Dateierweiterung `.fastq` oder `.fastq.gz` | Nicht erfasst (?: ist nicht erfassend) |

Dies parst Illumina-Namenskonventionen, um Metadaten automatisch zu extrahieren.

Führe den geänderten Workflow aus:

```bash title="Pattern Matching testen"
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

Dies zeigt die Metadaten, die aus den Dateinamen angereichert wurden.

### 2.2. Dynamische Script-Generierung in Prozessen

Prozess-Script-Blöcke sind im Wesentlichen mehrzeilige Strings, die an die Shell übergeben werden. Du kannst **bedingte Logik** (if/else, ternäre Operatoren) verwenden, um dynamisch verschiedene Script-Strings basierend auf Eingabe-Eigenschaften zu generieren. Dies ist essentiell für den Umgang mit verschiedenen Eingabetypen – wie Single-End vs Paired-End-Sequenzierungs-Reads – ohne Prozessdefinitionen zu duplizieren.

Fügen wir einen Prozess zu unserem Workflow hinzu, der dieses Muster demonstriert. Öffne `modules/fastp.nf` und schau dir das an:

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

Der Prozess nimmt FASTQ-Dateien als Eingabe und führt das `fastp`-Tool aus, um Adapter zu trimmen und Reads niedriger Qualität zu filtern. Leider hat die Person, die diesen Prozess geschrieben hat, die Single-End-Reads nicht berücksichtigt, die wir in unserem Beispieldatensatz haben. Fügen wir ihn zu unserem Workflow hinzu und schauen, was passiert:

Füge zuerst das Modul in der allerersten Zeile deines `main.nf`-Workflows ein:

```groovy title="main.nf" linenums="1"
include { FASTP } from './modules/fastp.nf'
```

Dann modifiziere den `workflow`-Block, um den `ch_samples`-Channel mit dem `FASTP`-Prozess zu verbinden:

=== "Nachher"

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

Führe diesen geänderten Workflow aus:

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

Du siehst, dass der Prozess versucht, `fastp` mit einem `null`-Wert für die zweite Eingabedatei auszuführen, was dazu führt, dass er fehlschlägt. Dies liegt daran, dass unser Datensatz Single-End-Reads enthält, aber der Prozess fest codiert ist, um Paired-End-Reads zu erwarten (zwei Eingabedateien gleichzeitig).

Behebe dies, indem du bedingte Logik zum `FASTP`-Prozess-`script:`-Block hinzufügst. Ein if/else-Statement prüft die Anzahl der Read-Dateien und passt den Befehl entsprechend an.

=== "Nachher"

    ```groovy title="main.nf" linenums="10" hl_lines="3-27"
        script:
        // Einfache Single-End vs Paired-End-Erkennung
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

Jetzt kann der Workflow sowohl Single-End- als auch Paired-End-Reads elegant handhaben. Die bedingte Logik prüft die Anzahl der Eingabedateien und konstruiert den entsprechenden Befehl für `fastp`. Schauen wir, ob es funktioniert:

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

Sieht gut aus! Wenn wir die tatsächlichen Befehle überprüfen, die ausgeführt wurden (passe für deinen Task-Hash an):

```console title="Ausgeführte Befehle prüfen"
cat work/31/a8ad4d95749e685a6d842d3007957f/.command.sh
```

Können wir sehen, dass Nextflow korrekt den richtigen Befehl für Single-End-Reads ausgewählt hat:

```bash title=".command.sh"
#!/bin/bash -ue
fastp \
    --in1 SAMPLE_003_S3_L001_R1_001.fastq \
    --out1 sample_003_trimmed.fastq.gz \
    --json sample_003.fastp.json \
    --html sample_003.fastp.html \
    --thread 2
```

Eine weitere häufige Verwendung dynamischer Script-Logik ist im [Nextflow for Science Genomics-Modul](../../nf4science/genomics/02_joint_calling) zu sehen. In diesem Modul kann der aufgerufene GATK-Prozess mehrere Eingabedateien annehmen, aber jede muss mit `-V` präfixiert werden, um eine korrekte Befehlszeile zu bilden. Der Prozess verwendet Scripting, um eine Collection von Eingabedateien (`all_gvcfs`) in die korrekten Befehlsargumente zu transformieren:

```groovy title="Befehlszeilenmanipulation für GATK" linenums="1" hl_lines="2 5"
    script:
    def gvcfs_line = all_gvcfs.collect { gvcf -> "-V ${gvcf}" }.join(' ')
    """
    gatk GenomicsDBImport \
        ${gvcfs_line} \
        -L ${interval_list} \
        --genomicsdb-workspace-path ${cohort_name}_gdb
    """
```

Diese Muster der Verwendung von Scripting in Prozess-Script-Blöcken sind extrem leistungsfähig und können in vielen Szenarien angewendet werden – vom Umgang mit variablen Eingabetypen bis zum Erstellen komplexer Befehlszeilenargumente aus Datei-Collections, wodurch deine Prozesse wirklich an die vielfältigen Anforderungen realer Daten anpassbar werden.

### 2.3. Variable Interpolation: Nextflow- und Shell-Variablen

Prozess-Scripts mischen Nextflow-Variablen, Shell-Variablen und Command-Substitutions, jeweils mit unterschiedlicher Interpolations-Syntax. Die Verwendung der falschen Syntax verursacht Fehler. Lass uns diese mit einem Prozess erkunden, der einen Verarbeitungsbericht erstellt.

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
    echo "Verarbeite ${reads}" > ${meta.id}_report.txt
    echo "Probe: ${meta.id}" >> ${meta.id}_report.txt
    """
}
```

Dieser Prozess schreibt einen einfachen Bericht mit der Proben-ID und dem Dateinamen. Lass uns ihn ausführen, um zu sehen, was passiert, wenn wir verschiedene Variablentypen mischen müssen.

Füge den Prozess in deine `main.nf` ein und füge ihn zum Workflow hinzu:

=== "Nachher"

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

                def
