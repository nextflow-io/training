# Aufteilen und Gruppieren

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Nextflow bietet leistungsstarke Werkzeuge für die flexible Arbeit mit Daten. Eine wichtige Funktion ist das Aufteilen von Daten in verschiedene Streams und das anschließende Zusammenführen verwandter Elemente. Das ist besonders wertvoll in Bioinformatik-Workflows, wo du verschiedene Probentypen getrennt verarbeiten musst, bevor du die Ergebnisse für die Analyse zusammenführst.

Stell es dir wie das Sortieren von Post vor: Du trennst Briefe nach Zielort, verarbeitest jeden Stapel unterschiedlich und kombinierst dann die Elemente, die an dieselbe Person gehen, wieder. Nextflow verwendet spezielle Operatoren, um dies mit wissenschaftlichen Daten zu erreichen. Dieser Ansatz ist in der verteilten Datenverarbeitung und in Bioinformatik-Workflows auch als **scatter/gather**-Muster bekannt.

Das Kanalsystem von Nextflow ist das Herzstück dieser Flexibilität. Kanäle verbinden verschiedene Teile deines Workflows und ermöglichen den Datenfluss durch deine Analyse. Du kannst mehrere Kanäle aus einer einzigen Datenquelle erstellen, jeden Kanal unterschiedlich verarbeiten und Kanäle bei Bedarf wieder zusammenführen. Dieser Ansatz ermöglicht es dir, Workflows zu entwerfen, die die verzweigten und zusammenlaufenden Pfade komplexer Bioinformatik-Analysen auf natürliche Weise widerspiegeln.

### Lernziele

In dieser Side Quest lernst du, Daten mithilfe der Kanal-Operatoren von Nextflow aufzuteilen und zu gruppieren.
Wir beginnen mit einer CSV-Datei, die Probeninformationen und zugehörige Datendateien enthält, und manipulieren und reorganisieren diese Daten dann.

Am Ende dieser Side Quest kannst du Datenströme effektiv trennen und kombinieren, indem du folgende Techniken verwendest:

- Daten aus Dateien mit `splitCsv` lesen
- Daten mit `filter` und `map` filtern und transformieren
- Verwandte Daten mit `join` und `groupTuple` kombinieren
- Datenkombinationen mit `combine` für die Parallelverarbeitung erstellen
- Die Datenstruktur mit `subMap` und Deduplizierungsstrategien optimieren
- Wiederverwendbare Funktionen mit benannten Closures erstellen, um Kanalstrukturen zu manipulieren

Diese Fähigkeiten helfen dir, Workflows zu erstellen, die mehrere Eingabedateien und verschiedene Datentypen effizient verarbeiten können, während du eine saubere, wartbare Codestruktur beibehältst.

### Voraussetzungen

Bevor du diese Side Quest in Angriff nimmst, solltest du:

- Das Tutorial [Hello Nextflow](../hello_nextflow/README.md) oder einen gleichwertigen Einsteigerkurs abgeschlossen haben.
- Mit grundlegenden Nextflow-Konzepten und -Mechanismen vertraut sein (Prozesse, Kanäle, Operatoren, Arbeiten mit Dateien, Metadaten)

**Optional:** Wir empfehlen, zuerst die Side Quest [Metadata in workflows](../metadata/) abzuschließen.
Diese behandelt die Grundlagen des Lesens von CSV-Dateien mit `splitCsv` und das Erstellen von Meta-Maps, die wir hier intensiv verwenden werden.

---

## 0. Erste Schritte

#### Trainings-Codespace öffnen

Falls du das noch nicht getan hast, öffne die Trainingsumgebung wie in der [Umgebung einrichten](../envsetup/index.md) beschrieben.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### In das Projektverzeichnis wechseln

Wechseln wir in das Verzeichnis, in dem sich die Dateien für dieses Tutorial befinden.

```bash
cd side-quests/splitting_and_grouping
```

Du kannst VSCode so einstellen, dass es sich auf dieses Verzeichnis konzentriert:

```bash
code .
```

#### Die Materialien ansehen

Du findest eine Haupt-Workflow-Datei und ein `data`-Verzeichnis mit einem Samplesheet namens `samplesheet.csv`.

```console title="Directory contents"
.
├── data
│   └── samplesheet.csv
└── main.nf
```

Das Samplesheet enthält Informationen über Proben von verschiedenen Patient\*innen, einschließlich der Patienten-ID, der Wiederholungsnummer der Probe, des Typs (normal oder Tumor) und Pfade zu hypothetischen Datendateien (die eigentlich nicht existieren, aber wir tun so, als ob sie es täten).

```console title="samplesheet.csv"
id,repeat,type,bam
patientA,1,normal,patientA_rep1_normal.bam
patientA,1,tumor,patientA_rep1_tumor.bam
patientA,2,normal,patientA_rep2_normal.bam
patientA,2,tumor,patientA_rep2_tumor.bam
patientB,1,normal,patientB_rep1_normal.bam
patientB,1,tumor,patientB_rep1_tumor.bam
patientC,1,normal,patientC_rep1_normal.bam
patientC,1,tumor,patientC_rep1_tumor.bam
```

Dieses Samplesheet listet acht Proben von drei Patient\*innen (A, B, C) auf.

Für jede\*n Patient\*in haben wir Proben vom Typ `tumor` (typischerweise aus Tumorbiopsien stammend) oder `normal` (aus gesundem Gewebe oder Blut entnommen).
Falls du nicht mit der Krebsanalyse vertraut bist: Dies entspricht einem experimentellen Modell, das gepaarte Tumor/Normal-Proben für kontrastive Analysen verwendet.

Für Patient A haben wir speziell zwei Sätze technischer Replikate (Wiederholungen).

!!! note "Hinweis"

    Mach dir keine Sorgen, wenn du mit diesem experimentellen Design nicht vertraut bist – es ist nicht entscheidend für das Verständnis dieses Tutorials.

#### Die Aufgabe ansehen

Deine Aufgabe ist es, einen Nextflow-Workflow zu schreiben, der:

1. **Probendaten** aus einer CSV-Datei liest und mit Meta-Maps strukturiert
2. **Proben** basierend auf dem Typ (normal vs. Tumor) in verschiedene Kanäle trennt
3. **Übereinstimmende** Tumor/Normal-Paare nach Patienten-ID und Replikatnummer zusammenführt
4. **Proben** über genomische Intervalle für die Parallelverarbeitung verteilt
5. **Verwandte** Proben für die nachgelagerte Analyse wieder zusammenführt

Dies ist ein häufiges Bioinformatik-Muster, bei dem du Daten für die unabhängige Verarbeitung aufteilen und dann verwandte Elemente für die vergleichende Analyse wieder zusammenführen musst.

#### Bereitschafts-Checkliste

Bereit zum Loslegen?

- [ ] Ich verstehe das Ziel dieses Kurses und seine Voraussetzungen
- [ ] Mein Codespace läuft
- [ ] Ich habe mein Arbeitsverzeichnis entsprechend eingestellt
- [ ] Ich verstehe die Aufgabe

Wenn du alle Punkte abhaken kannst, kann es losgehen.

---

## 1. Probendaten einlesen

### 1.1. Probendaten mit `splitCsv` einlesen und Meta-Maps erstellen

Beginnen wir damit, die Probendaten mit `splitCsv` einzulesen und sie in das Meta-Map-Muster zu organisieren. In der `main.nf` siehst du, dass wir den Workflow bereits begonnen haben.

```groovy title="main.nf" linenums="1" hl_lines="2"
workflow {
    ch_samplesheet = channel.fromPath("./data/samplesheet.csv")
}
```

!!! note "Hinweis"

    In diesem Tutorial verwenden wir das Präfix `ch_` für alle Kanalvariablen, um klar anzuzeigen, dass es sich um Nextflow-Kanäle handelt.

Wenn du die Side Quest [Metadata in workflows](../metadata/) abgeschlossen hast, wirst du dieses Muster wiedererkennen. Wir verwenden `splitCsv`, um die CSV zu lesen, und strukturieren die Daten sofort mit einer Meta-Map, um Metadaten von Dateipfaden zu trennen.

!!! info "Info"

    In diesem Training begegnen uns zwei verschiedene Konzepte namens `map`:

    - **Datenstruktur**: Die Groovy-Map (entspricht Dictionaries/Hashes in anderen Sprachen), die Schlüssel-Wert-Paare speichert
    - **Kanal-Operator**: Der `.map()`-Operator, der Elemente in einem Kanal transformiert

    Wir werden im Kontext klarstellen, welches wir meinen, aber diese Unterscheidung ist wichtig, wenn du mit Nextflow arbeitest.

Wende diese Änderungen auf `main.nf` an:

=== "Danach"

    ```groovy title="main.nf" linenums="2" hl_lines="2-6"
        ch_samples = channel.fromPath("./data/samplesheet.csv")
            .splitCsv(header: true)
            .map{ row ->
              [[id:row.id, repeat:row.repeat, type:row.type], row.bam]
            }
            .view()
    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="2" hl_lines="1"
        ch_samplesheet = channel.fromPath("./data/samplesheet.csv")
    ```

Dies kombiniert die `splitCsv`-Operation (Lesen der CSV mit Kopfzeilen) und die `map`-Operation (Strukturierung der Daten als `[meta, file]`-Tupel) in einem Schritt. Wende diese Änderung an und führe die Pipeline aus:

```bash
nextflow run main.nf
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [deadly_mercator] DSL2 - revision: bd6b0224e9

    [[id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam]
    [[id:patientA, repeat:1, type:tumor], patientA_rep1_tumor.bam]
    [[id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam]
    [[id:patientA, repeat:2, type:tumor], patientA_rep2_tumor.bam]
    [[id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam]
    [[id:patientB, repeat:1, type:tumor], patientB_rep1_tumor.bam]
    [[id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam]
    [[id:patientC, repeat:1, type:tumor], patientC_rep1_tumor.bam]
    ```

Wir haben jetzt einen Kanal, bei dem jedes Element ein `[meta, file]`-Tupel ist – Metadaten getrennt von Dateipfaden. Diese Struktur ermöglicht es uns, unsere Arbeitslast basierend auf Metadatenfeldern aufzuteilen und zu gruppieren.

---

## 2. Daten filtern und transformieren

### 2.1. Daten mit `filter` filtern

Wir können den [`filter`-Operator](https://www.nextflow.io/docs/latest/operator.html#filter) verwenden, um die Daten basierend auf einer Bedingung zu filtern. Angenommen, wir möchten nur normale Proben verarbeiten. Wir können dies tun, indem wir die Daten basierend auf dem `type`-Feld filtern. Fügen wir dies vor dem `view`-Operator ein.

=== "Danach"

    ```groovy title="main.nf" linenums="2" hl_lines="6"
        ch_samples = channel.fromPath("./data/samplesheet.csv")
            .splitCsv(header: true)
            .map{ row ->
              [[id:row.id, repeat:row.repeat, type:row.type], row.bam]
            }
            .filter { meta, file -> meta.type == 'normal' }
            .view()
    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="2"
        ch_samples = channel.fromPath("./data/samplesheet.csv")
            .splitCsv(header: true)
            .map{ row ->
              [[id:row.id, repeat:row.repeat, type:row.type], row.bam]
            }
            .view()
    ```

Führe den Workflow erneut aus, um das gefilterte Ergebnis zu sehen:

```bash
nextflow run main.nf
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [admiring_brown] DSL2 - revision: 194d61704d

    [[id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam]
    [[id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam]
    [[id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam]
    [[id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam]
    ```

Wir haben die Daten erfolgreich gefiltert, sodass nur normale Proben enthalten sind. Lass uns kurz zusammenfassen, wie das funktioniert.

Der `filter`-Operator nimmt eine Closure, die auf jedes Element im Kanal angewendet wird. Wenn die Closure `true` zurückgibt, wird das Element eingeschlossen; wenn sie `false` zurückgibt, wird das Element ausgeschlossen.

In unserem Fall möchten wir nur Proben behalten, bei denen `meta.type == 'normal'`. Die Closure verwendet das Tupel `meta,file`, um auf jede Probe zu verweisen, greift mit `meta.type` auf den Probentyp zu und prüft, ob er gleich `'normal'` ist.

Dies wird mit der einzelnen Closure erreicht, die wir oben eingeführt haben:

```groovy title="main.nf" linenums="7"
    .filter { meta, file -> meta.type == 'normal' }
```

### 2.2. Separate gefilterte Kanäle erstellen

Derzeit wenden wir den Filter auf den direkt aus der CSV erstellten Kanal an, aber wir möchten diesen auf mehr als eine Weise filtern. Schreiben wir die Logik also um, um einen separaten gefilterten Kanal für normale Proben zu erstellen:

=== "Danach"

    ```groovy title="main.nf" linenums="2" hl_lines="6 8"
        ch_samples = channel.fromPath("./data/samplesheet.csv")
            .splitCsv(header: true)
            .map{ row ->
                [[id:row.id, repeat:row.repeat, type:row.type], row.bam]
            }
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
        ch_normal_samples
            .view()
    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="2"
        ch_samples = channel.fromPath("./data/samplesheet.csv")
            .splitCsv(header: true)
            .map{ row ->
              [[id:row.id, repeat:row.repeat, type:row.type], row.bam]
            }
            .filter { meta, file -> meta.type == 'normal' }
            .view()
    ```

Führe die Pipeline aus, um die Ergebnisse zu sehen:

```bash
nextflow run main.nf
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [trusting_poisson] DSL2 - revision: 639186ee74

    [[id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam]
    [[id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam]
    [[id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam]
    [[id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam]
    ```

Wir haben die Daten erfolgreich gefiltert und einen separaten Kanal für normale Proben erstellt.

Erstellen wir auch einen gefilterten Kanal für die Tumorproben:

=== "Danach"

    ```groovy title="main.nf" linenums="7" hl_lines="3-8"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
        ch_normal_samples
            .view{'Normal sample: ' + it}
        ch_tumor_samples
            .view{'Tumor sample: ' + it}
    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="7" hl_lines="3 4"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
        ch_normal_samples
            .view()
    ```

```bash
nextflow run main.nf
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [maniac_boltzmann] DSL2 - revision: 3636b6576b

    Tumor sample: [[id:patientA, repeat:1, type:tumor], patientA_rep1_tumor.bam]
    Tumor sample: [[id:patientA, repeat:2, type:tumor], patientA_rep2_tumor.bam]
    Normal sample: [[id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam]
    Normal sample: [[id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam]
    Normal sample: [[id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam]
    Normal sample: [[id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam]
    Tumor sample: [[id:patientB, repeat:1, type:tumor], patientB_rep1_tumor.bam]
    Tumor sample: [[id:patientC, repeat:1, type:tumor], patientC_rep1_tumor.bam]
    ```

Wir haben die normalen und Tumorproben in zwei verschiedene Kanäle getrennt und eine Closure verwendet, die an `view()` übergeben wird, um sie in der Ausgabe unterschiedlich zu kennzeichnen: `ch_tumor_samples.view{'Tumor sample: ' + it}`.

### Fazit

In diesem Abschnitt hast du gelernt:

- **Daten filtern**: Wie man Daten mit `filter` filtert
- **Daten aufteilen**: Wie man Daten basierend auf einer Bedingung in verschiedene Kanäle aufteilt
- **Daten anzeigen**: Wie man `view` verwendet, um Daten auszugeben und die Ausgabe verschiedener Kanäle zu kennzeichnen

Wir haben jetzt die normalen und Tumorproben in zwei verschiedene Kanäle getrennt. Als Nächstes werden wir die normalen und Tumorproben über das `id`-Feld zusammenführen.

---

## 3. Kanäle über Bezeichner zusammenführen

Im vorherigen Abschnitt haben wir die normalen und Tumorproben in zwei verschiedene Kanäle getrennt. Diese könnten unabhängig voneinander mit spezifischen Prozessen oder Workflows basierend auf ihrem Typ verarbeitet werden. Aber was passiert, wenn wir die normalen und Tumorproben desselben Patienten vergleichen möchten? An diesem Punkt müssen wir sie wieder zusammenführen und dabei sicherstellen, dass die Proben basierend auf ihrem `id`-Feld abgeglichen werden.

Nextflow enthält viele Methoden zum Kombinieren von Kanälen, aber in diesem Fall ist der am besten geeignete Operator [`join`](https://www.nextflow.io/docs/latest/operator.html#join). Wenn du mit SQL vertraut bist, verhält er sich wie die `JOIN`-Operation, bei der wir den Schlüssel angeben, nach dem zusammengeführt werden soll, und die Art des Joins.

### 3.1. `map` und `join` verwenden, um nach Patienten-ID zu kombinieren

Wenn wir die [`join`](https://www.nextflow.io/docs/latest/operator.html#join)-Dokumentation prüfen, sehen wir, dass standardmäßig zwei Kanäle basierend auf dem ersten Element in jedem Tupel zusammengeführt werden.

#### 3.1.1. Die Datenstruktur prüfen

Falls die Konsolenausgabe nicht mehr verfügbar ist, führen wir die Pipeline aus, um unsere Datenstruktur zu prüfen und zu sehen, wie wir sie ändern müssen, um nach dem `id`-Feld zusammenzuführen.

```bash
nextflow run main.nf
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [maniac_boltzmann] DSL2 - revision: 3636b6576b

    Tumor sample: [[id:patientA, repeat:1, type:tumor], patientA_rep1_tumor.bam]
    Tumor sample: [[id:patientA, repeat:2, type:tumor], patientA_rep2_tumor.bam]
    Normal sample: [[id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam]
    Normal sample: [[id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam]
    Normal sample: [[id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam]
    Normal sample: [[id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam]
    Tumor sample: [[id:patientB, repeat:1, type:tumor], patientB_rep1_tumor.bam]
    Tumor sample: [[id:patientC, repeat:1, type:tumor], patientC_rep1_tumor.bam]
    ```

Wir können sehen, dass das `id`-Feld das erste Element in jeder Meta-Map ist. Damit `join` funktioniert, sollten wir das `id`-Feld in jedem Tupel isolieren. Danach können wir einfach den `join`-Operator verwenden, um die beiden Kanäle zu kombinieren.

#### 3.1.2. Das `id`-Feld isolieren

Um das `id`-Feld zu isolieren, können wir den [`map`-Operator](https://www.nextflow.io/docs/latest/operator.html#map) verwenden, um ein neues Tupel mit dem `id`-Feld als erstem Element zu erstellen.

=== "Danach"

    ```groovy title="main.nf" linenums="7" hl_lines="3 6"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
            .map { meta, file -> [meta.id, meta, file] }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
            .map { meta, file -> [meta.id, meta, file] }
        ch_normal_samples
            .view{'Normal sample: ' + it}
        ch_tumor_samples
            .view{'Tumor sample: ' + it}
    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="7"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
        ch_normal_samples
            .view{'Normal sample: ' + it}
        ch_tumor_samples
            .view{'Tumor sample: ' + it}
    ```

```bash
nextflow run main.nf
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [mad_lagrange] DSL2 - revision: 9940b3f23d

    Tumor sample: [patientA, [id:patientA, repeat:1, type:tumor], patientA_rep1_tumor.bam]
    Tumor sample: [patientA, [id:patientA, repeat:2, type:tumor], patientA_rep2_tumor.bam]
    Normal sample: [patientA, [id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam]
    Normal sample: [patientA, [id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam]
    Tumor sample: [patientB, [id:patientB, repeat:1, type:tumor], patientB_rep1_tumor.bam]
    Tumor sample: [patientC, [id:patientC, repeat:1, type:tumor], patientC_rep1_tumor.bam]
    Normal sample: [patientB, [id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam]
    Normal sample: [patientC, [id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam]
    ```

Es mag subtil sein, aber du solltest erkennen können, dass das erste Element in jedem Tupel das `id`-Feld ist.

#### 3.1.3. Die beiden Kanäle kombinieren

Jetzt können wir den `join`-Operator verwenden, um die beiden Kanäle basierend auf dem `id`-Feld zu kombinieren.

Wir verwenden erneut `view`, um die zusammengeführten Ausgaben auszugeben.

=== "Danach"

    ```groovy title="main.nf" linenums="7" hl_lines="7-9"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
            .map { meta, file -> [meta.id, meta, file] }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
            .map { meta, file -> [meta.id, meta, file] }
        ch_joined_samples = ch_normal_samples
            .join(ch_tumor_samples)
        ch_joined_samples.view()
    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="7" hl_lines="7-10"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
            .map { meta, file -> [meta.id, meta, file] }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
            .map { meta, file -> [meta.id, meta, file] }
        ch_normal_samples
            .view{'Normal sample: ' + it}
        ch_tumor_samples
            .view{'Tumor sample: ' + it}
    ```

```bash
nextflow run main.nf
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [soggy_wiles] DSL2 - revision: 3bc1979889

    [patientA, [id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam, [id:patientA, repeat:1, type:tumor], patientA_rep1_tumor.bam]
    [patientA, [id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam, [id:patientA, repeat:2, type:tumor], patientA_rep2_tumor.bam]
    [patientB, [id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam, [id:patientB, repeat:1, type:tumor], patientB_rep1_tumor.bam]
    [patientC, [id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam, [id:patientC, repeat:1, type:tumor], patientC_rep1_tumor.bam]
    ```

Es ist etwas schwer zu erkennen, weil es so breit ist, aber du solltest sehen können, dass die Proben nach dem `id`-Feld zusammengeführt wurden. Jedes Tupel hat jetzt das Format:

- `id`: Die Proben-ID
- `normal_meta_map`: Die Metadaten der normalen Probe, einschließlich Typ, Replikat und Pfad zur BAM-Datei
- `normal_sample_file`: Die normale Probendatei
- `tumor_meta_map`: Die Metadaten der Tumorprobe, einschließlich Typ, Replikat und Pfad zur BAM-Datei
- `tumor_sample`: Die Tumorprobe, einschließlich Typ, Replikat und Pfad zur BAM-Datei

!!! warning "Warnung"

    Der `join`-Operator verwirft alle nicht übereinstimmenden Tupel. In diesem Beispiel haben wir sichergestellt, dass alle Proben für Tumor und Normal übereinstimmen, aber wenn das nicht der Fall ist, musst du den Parameter `remainder: true` verwenden, um die nicht übereinstimmenden Tupel zu behalten. Weitere Details findest du in der [Dokumentation](https://www.nextflow.io/docs/latest/operator.html#join).

Du weißt jetzt, wie du `map` verwendest, um ein Feld in einem Tupel zu isolieren, und wie du `join` verwendest, um Tupel basierend auf dem ersten Feld zu kombinieren.
Mit diesem Wissen können wir Kanäle erfolgreich basierend auf einem gemeinsamen Feld kombinieren.

Als Nächstes betrachten wir die Situation, in der du nach mehreren Feldern zusammenführen möchtest.

### 3.2. Nach mehreren Feldern zusammenführen

Wir haben 2 Replikate für sampleA, aber nur 1 für sampleB und sampleC. In diesem Fall konnten wir sie effektiv zusammenführen, indem wir das `id`-Feld verwendet haben, aber was würde passieren, wenn sie nicht synchron wären? Wir könnten die normalen und Tumorproben aus verschiedenen Replikaten durcheinanderbringen!

Um dies zu vermeiden, können wir nach mehreren Feldern zusammenführen. Es gibt tatsächlich mehrere Möglichkeiten, dies zu erreichen, aber wir konzentrieren uns auf die Erstellung eines neuen Zusammenführungsschlüssels, der sowohl die Proben-`id` als auch die `replicate`-Nummer enthält.

Beginnen wir mit der Erstellung eines neuen Zusammenführungsschlüssels. Wir können dies auf die gleiche Weise wie zuvor tun, indem wir den [`map`-Operator](https://www.nextflow.io/docs/latest/operator.html#map) verwenden, um ein neues Tupel mit den Feldern `id` und `repeat` als erstem Element zu erstellen.

=== "Danach"

    ```groovy title="main.nf" linenums="7" hl_lines="3 6"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
            .map { meta, file -> [[meta.id, meta.repeat], meta, file] }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
            .map { meta, file -> [[meta.id, meta.repeat], meta, file] }
    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="7" hl_lines="3 6"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
            .map { meta, file -> [meta.id, meta, file] }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
            .map { meta, file -> [meta.id, meta, file] }
    ```

Jetzt sollten wir sehen, dass die Zusammenführung stattfindet, aber sowohl die Felder `id` als auch `repeat` verwendet werden. Führe den Workflow aus:

```bash
nextflow run main.nf
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [prickly_wing] DSL2 - revision: 3bebf22dee

    [[patientA, 1], [id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam, [id:patientA, repeat:1, type:tumor], patientA_rep1_tumor.bam]
    [[patientA, 2], [id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam, [id:patientA, repeat:2, type:tumor], patientA_rep2_tumor.bam]
    [[patientB, 1], [id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam, [id:patientB, repeat:1, type:tumor], patientB_rep1_tumor.bam]
    [[patientC, 1], [id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam, [id:patientC, repeat:1, type:tumor], patientC_rep1_tumor.bam]
    ```

Beachte, dass wir ein Tupel aus zwei Elementen (Felder `id` und `repeat`) als erstes Element jedes zusammengeführten Ergebnisses haben. Dies zeigt, wie komplexe Elemente als Zusammenführungsschlüssel verwendet werden können, was eine recht präzise Übereinstimmung zwischen Proben aus denselben Bedingungen ermöglicht.

Wenn du weitere Möglichkeiten erkunden möchtest, nach verschiedenen Schlüsseln zusammenzuführen, schau dir die [Dokumentation des join-Operators](https://www.nextflow.io/docs/latest/operator.html#join) für weitere Optionen und Beispiele an.

### 3.3. `subMap` verwenden, um einen neuen Zusammenführungsschlüssel zu erstellen

Der vorherige Ansatz verliert die Feldnamen aus unserem Zusammenführungsschlüssel – die Felder `id` und `repeat` werden zu einer einfachen Werteliste. Um die Feldnamen für den späteren Zugriff beizubehalten, können wir die [`subMap`-Methode](<https://docs.groovy-lang.org/latest/html/groovy-jdk/java/util/Map.html#subMap(java.util.Collection)>) verwenden.

Die `subMap`-Methode extrahiert nur die angegebenen Schlüssel-Wert-Paare aus einer Map. Hier extrahieren wir nur die Felder `id` und `repeat`, um unseren Zusammenführungsschlüssel zu erstellen.

=== "Danach"

    ```groovy title="main.nf" linenums="7" hl_lines="3 6"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
            .map { meta, file -> [meta.subMap(['id', 'repeat']), meta, file] }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
            .map { meta, file -> [meta.subMap(['id', 'repeat']), meta, file] }
    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="7" hl_lines="3 6"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
            .map { meta, file -> [[meta.id, meta.repeat], meta, file] }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
            .map { meta, file -> [[meta.id, meta.repeat], meta, file] }
    ```

```bash
nextflow run main.nf
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [reverent_wing] DSL2 - revision: 847016c3b7

    [[id:patientA, repeat:1], [id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam, [id:patientA, repeat:1, type:tumor], patientA_rep1_tumor.bam]
    [[id:patientA, repeat:2], [id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam, [id:patientA, repeat:2, type:tumor], patientA_rep2_tumor.bam]
    [[id:patientB, repeat:1], [id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam, [id:patientB, repeat:1, type:tumor], patientB_rep1_tumor.bam]
    [[id:patientC, repeat:1], [id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam, [id:patientC, repeat:1, type:tumor], patientC_rep1_tumor.bam]
    ```

Jetzt haben wir einen neuen Zusammenführungsschlüssel, der nicht nur die Felder `id` und `repeat` enthält, sondern auch die Feldnamen beibehält, sodass wir später über den Namen darauf zugreifen können, z. B. `meta.id` und `meta.repeat`.

### 3.4. Eine benannte Closure in map verwenden

Um Duplikate zu vermeiden und Fehler zu reduzieren, können wir eine benannte Closure verwenden. Eine benannte Closure ermöglicht es uns, eine wiederverwendbare Funktion zu erstellen, die wir an mehreren Stellen aufrufen können.

Dazu definieren wir zunächst die Closure als neue Variable:

=== "Danach"

    ```groovy title="main.nf" linenums="2" hl_lines="7"
        ch_samples = channel.fromPath("./data/samplesheet.csv")
            .splitCsv(header: true)
            .map{ row ->
                [[id:row.id, repeat:row.repeat, type:row.type], row.bam]
            }

        getSampleIdAndReplicate = { meta, bam -> [ meta.subMap(['id', 'repeat']), meta, file(bam) ] }

        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="2"
        ch_samples = channel.fromPath("./data/samplesheet.csv")
            .splitCsv(header: true)
            .map{ row ->
                [[id:row.id, repeat:row.repeat, type:row.type], row.bam]
            }
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
    ```

Wir haben die Map-Transformation als benannte Variable definiert, die wir wiederverwenden können.

Beachte, dass wir den Dateipfad auch mit `file()` in ein Path-Objekt umwandeln, damit jeder Prozess, der diesen Kanal empfängt, die Datei korrekt verarbeiten kann (weitere Informationen findest du unter [Working with files](../working_with_files/)).

Implementieren wir die Closure in unserem Workflow:

=== "Danach"

    ```groovy title="main.nf" linenums="10" hl_lines="3 6"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
             .map ( getSampleIdAndReplicate )
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
             .map ( getSampleIdAndReplicate )

    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="10" hl_lines="3 6"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
            .map { meta, file -> [meta.subMap(['id', 'repeat']), meta, file] }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
            .map { meta, file -> [meta.subMap(['id', 'repeat']), meta, file] }
    ```

!!! note "Hinweis"

    Der `map`-Operator hat von `{ }` auf `( )` gewechselt, um die Closure als Argument zu übergeben. Das liegt daran, dass der `map`-Operator eine Closure als Argument erwartet und `{ }` zur Definition einer anonymen Closure verwendet wird. Beim Aufrufen einer benannten Closure verwende die `( )`-Syntax.

Führe den Workflow noch einmal aus, um zu prüfen, ob alles noch funktioniert:

```bash
nextflow run main.nf
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [angry_meninsky] DSL2 - revision: 2edc226b1d

    [[id:patientA, repeat:1], [id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam, [id:patientA, repeat:1, type:tumor], patientA_rep1_tumor.bam]
    [[id:patientA, repeat:2], [id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam, [id:patientA, repeat:2, type:tumor], patientA_rep2_tumor.bam]
    [[id:patientB, repeat:1], [id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam, [id:patientB, repeat:1, type:tumor], patientB_rep1_tumor.bam]
    [[id:patientC, repeat:1], [id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam, [id:patientC, repeat:1, type:tumor], patientC_rep1_tumor.bam]
    ```

Die Verwendung einer benannten Closure ermöglicht es uns, dieselbe Transformation an mehreren Stellen wiederzuverwenden, was das Fehlerrisiko verringert und den Code lesbarer und wartbarer macht.

### 3.5. Datenduplikate reduzieren

Wir haben viele duplizierte Daten in unserem Workflow. Jedes Element in den zusammengeführten Proben wiederholt die Felder `id` und `repeat`. Da diese Informationen bereits im Gruppierungsschlüssel verfügbar sind, können wir diese Redundanz vermeiden. Zur Erinnerung: Unsere aktuelle Datenstruktur sieht so aus:

```groovy
[
  [
    "id": "sampleC",
    "repeat": "1",
  ],
  [
    "id": "sampleC",
    "repeat": "1",
    "type": "normal",
  ],
  "sampleC_rep1_normal.bam"
  [
    "id": "sampleC",
    "repeat": "1",
    "type": "tumor",
  ],
  "sampleC_rep1_tumor.bam"
]
```

Da die Felder `id` und `repeat` im Gruppierungsschlüssel verfügbar sind, entfernen wir sie aus dem Rest jedes Kanalelements, um Duplikate zu vermeiden. Wir können dies tun, indem wir die `subMap`-Methode verwenden, um eine neue Map mit nur dem `type`-Feld zu erstellen. Dieser Ansatz ermöglicht es uns, alle notwendigen Informationen beizubehalten und gleichzeitig Redundanz in unserer Datenstruktur zu eliminieren.

=== "Danach"

    ```groovy title="main.nf" linenums="8" hl_lines="1"
        getSampleIdAndReplicate = { meta, bam -> [ meta.subMap(['id', 'repeat']), meta.subMap(['type']), file(bam) ] }
    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="8" hl_lines="1"
        getSampleIdAndReplicate = { meta, bam -> [ meta.subMap(['id', 'repeat']), meta, file(bam) ] }
    ```

Jetzt gibt die Closure ein Tupel zurück, bei dem das erste Element die Felder `id` und `repeat` enthält und das zweite Element nur das `type`-Feld enthält. Dies eliminiert Redundanz, indem die `id`- und `repeat`-Informationen einmal im Gruppierungsschlüssel gespeichert werden, während alle notwendigen Informationen erhalten bleiben.

Führe den Workflow aus, um zu sehen, wie das aussieht:

```bash
nextflow run main.nf
```

??? success "Befehlsausgabe"

    ```console
    [[id:patientA, repeat:1], [type:normal], /workspaces/training/side-quests/splitting_and_grouping/patientA_rep1_normal.bam, [type:tumor], /workspaces/training/side-quests/splitting_and_grouping/patientA_rep1_tumor.bam]
    [[id:patientA, repeat:2], [type:normal], /workspaces/training/side-quests/splitting_and_grouping/patientA_rep2_normal.bam, [type:tumor], /workspaces/training/side-quests/splitting_and_grouping/patientA_rep2_tumor.bam]
    [[id:patientB, repeat:1], [type:normal], /workspaces/training/side-quests/splitting_and_grouping/patientB_rep1_normal.bam, [type:tumor], /workspaces/training/side-quests/splitting_and_grouping/patientB_rep1_tumor.bam]
    [[id:patientC, repeat:1], [type:normal], /workspaces/training/side-quests/splitting_and_grouping/patientC_rep1_normal.bam, [type:tumor], /workspaces/training/side-quests/splitting_and_grouping/patientC_rep1_tumor.bam]
    ```

Wir können sehen, dass wir die Felder `id` und `repeat` nur einmal im Gruppierungsschlüssel angeben und das `type`-Feld in den Probendaten haben. Wir haben keine Informationen verloren, aber es ist uns gelungen, den Inhalt unseres Kanals kompakter zu gestalten.

### 3.6. Redundante Informationen entfernen

Wir haben oben duplizierte Informationen entfernt, aber es gibt noch andere redundante Informationen in unseren Kanälen.

Am Anfang haben wir die normalen und Tumorproben mit `filter` getrennt und sie dann basierend auf den Schlüsseln `id` und `repeat` zusammengeführt. Der `join`-Operator behält die Reihenfolge bei, in der Tupel zusammengeführt werden. In unserem Fall, mit normalen Proben auf der linken Seite und Tumorproben auf der rechten, behält der resultierende Kanal diese Struktur bei: `id, <normale Elemente>, <Tumor-Elemente>`.

Da wir die Position jedes Elements in unserem Kanal kennen, können wir die Struktur weiter vereinfachen, indem wir die Metadaten `[type:normal]` und `[type:tumor]` weglassen.

=== "Danach"

    ```groovy title="main.nf" linenums="8" hl_lines="1"
        getSampleIdAndReplicate = { meta, file -> [ meta.subMap(['id', 'repeat']), file ] }
    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="8" hl_lines="1"
        getSampleIdAndReplicate = { meta, file -> [ meta.subMap(['id', 'repeat']), meta.subMap(['type']), file ] }
    ```

Führe den Workflow erneut aus, um das Ergebnis zu sehen:

```bash
nextflow run main.nf
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [confident_leavitt] DSL2 - revision: a2303895bd

    [[id:patientA, repeat:1], patientA_rep1_normal.bam, patientA_rep1_tumor.bam]
    [[id:patientA, repeat:2], patientA_rep2_normal.bam, patientA_rep2_tumor.bam]
    [[id:patientB, repeat:1], patientB_rep1_normal.bam, patientB_rep1_tumor.bam]
    [[id:patientC, repeat:1], patientC_rep1_normal.bam, patientC_rep1_tumor.bam]
    ```

### Fazit

In diesem Abschnitt hast du gelernt:

- **Tupel manipulieren**: Wie man `map` verwendet, um ein Feld in einem Tupel zu isolieren
- **Tupel zusammenführen**: Wie man `join` verwendet, um Tupel basierend auf dem ersten Feld zu kombinieren
- **Zusammenführungsschlüssel erstellen**: Wie man `subMap` verwendet, um einen neuen Zusammenführungsschlüssel zu erstellen
- **Benannte Closures**: Wie man eine benannte Closure in map verwendet
- **Zusammenführung nach mehreren Feldern**: Wie man nach mehreren Feldern für eine präzisere Übereinstimmung zusammenführt
- **Datenstruktur optimieren**: Wie man die Kanalstruktur durch Entfernen redundanter Informationen vereinfacht

Du hast jetzt einen Workflow, der ein Samplesheet aufteilen, die normalen und Tumorproben filtern, sie nach Proben-ID und Replikatnummer zusammenführen und die Ergebnisse ausgeben kann.

Dies ist ein häufiges Muster in Bioinformatik-Workflows, bei dem du Proben oder andere Datentypen nach unabhängiger Verarbeitung zusammenführen musst – eine nützliche Fähigkeit. Als Nächstes schauen wir uns an, wie eine Probe mehrfach wiederholt werden kann.

## 4. Proben über Intervalle verteilen

Ein wichtiges Muster in Bioinformatik-Workflows ist die Verteilung der Analyse über genomische Regionen. Zum Beispiel kann die Variantenidentifikation parallelisiert werden, indem das Genom in Intervalle (wie Chromosomen oder kleinere Regionen) aufgeteilt wird. Diese Parallelisierungsstrategie verbessert die Pipeline-Effizienz erheblich, indem die Rechenlast auf mehrere Kerne oder Knoten verteilt wird, was die Gesamtausführungszeit reduziert.

Im folgenden Abschnitt zeigen wir, wie unsere Probendaten über mehrere genomische Intervalle verteilt werden. Wir paaren jede Probe mit jedem Intervall und ermöglichen so die parallele Verarbeitung verschiedener genomischer Regionen. Dadurch wird unsere Datenmenge um die Anzahl der Intervalle multipliziert, wodurch mehrere unabhängige Analyseeinheiten entstehen, die später wieder zusammengeführt werden können.

### 4.1. Proben über Intervalle mit `combine` verteilen

Beginnen wir mit der Erstellung eines Kanals mit Intervallen. Der Einfachheit halber verwenden wir nur 3 Intervalle, die wir manuell definieren. In einem echten Workflow könntest du diese aus einer Dateieingabe lesen oder sogar einen Kanal mit vielen Intervalldateien erstellen.

=== "Danach"

    ```groovy title="main.nf" linenums="17" hl_lines="2"
            .join(ch_tumor_samples)
        ch_intervals = channel.of('chr1', 'chr2', 'chr3')
    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="17" hl_lines="2"
            .join(ch_tumor_samples)
        ch_joined_samples.view()
    ```

Denk daran: Wir möchten jede Probe für jedes Intervall wiederholen. Dies wird manchmal als kartesisches Produkt der Proben und Intervalle bezeichnet. Wir können dies mit dem [`combine`-Operator](https://www.nextflow.io/docs/latest/operator.html#combine) erreichen. Dieser nimmt jedes Element aus Kanal 1 und wiederholt es für jedes Element in Kanal 2. Fügen wir einen combine-Operator zu unserem Workflow hinzu:

=== "Danach"

    ```groovy title="main.nf" linenums="18" hl_lines="3-5"
        ch_intervals = channel.of('chr1', 'chr2', 'chr3')

        ch_combined_samples = ch_joined_samples
            .combine(ch_intervals)
            .view()
    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="18"
        ch_intervals = channel.of('chr1', 'chr2', 'chr3')
    ```

Führen wir es aus und schauen, was passiert:

```bash
nextflow run main.nf
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [mighty_tesla] DSL2 - revision: ae013ab70b

    [[id:patientA, repeat:1], patientA_rep1_normal.bam, patientA_rep1_tumor.bam, chr1]
    [[id:patientA, repeat:1], patientA_rep1_normal.bam, patientA_rep1_tumor.bam, chr2]
    [[id:patientA, repeat:1], patientA_rep1_normal.bam, patientA_rep1_tumor.bam, chr3]
    [[id:patientA, repeat:2], patientA_rep2_normal.bam, patientA_rep2_tumor.bam, chr1]
    [[id:patientA, repeat:2], patientA_rep2_normal.bam, patientA_rep2_tumor.bam, chr2]
    [[id:patientA, repeat:2], patientA_rep2_normal.bam, patientA_rep2_tumor.bam, chr3]
    [[id:patientB, repeat:1], patientB_rep1_normal.bam, patientB_rep1_tumor.bam, chr1]
    [[id:patientB, repeat:1], patientB_rep1_normal.bam, patientB_rep1_tumor.bam, chr2]
    [[id:patientB, repeat:1], patientB_rep1_normal.bam, patientB_rep1_tumor.bam, chr3]
    [[id:patientC, repeat:1], patientC_rep1_normal.bam, patientC_rep1_tumor.bam, chr1]
    [[id:patientC, repeat:1], patientC_rep1_normal.bam, patientC_rep1_tumor.bam, chr2]
    [[id:patientC, repeat:1], patientC_rep1_normal.bam, patientC_rep1_tumor.bam, chr3]
    ```

Erfolg! Wir haben jede Probe für jedes einzelne Intervall in unserer 3-Intervall-Liste wiederholt. Wir haben die Anzahl der Elemente in unserem Kanal effektiv verdreifacht.

Es ist etwas schwer zu lesen, daher werden wir es im nächsten Abschnitt aufräumen.

### 4.2. Den Kanal organisieren

Wir können den `map`-Operator verwenden, um unsere Probendaten zu bereinigen und umzustrukturieren, damit sie leichter zu verstehen sind. Verschieben wir den Intervall-String in die Zusammenführungs-Map als erstes Element.

=== "Danach"

    ```groovy title="main.nf" linenums="20" hl_lines="3-9"
        ch_combined_samples = ch_joined_samples
            .combine(ch_intervals)
            .map { grouping_key, normal, tumor, interval ->
                [
                    grouping_key + [interval: interval],
                    normal,
                    tumor
                ]
            }
            .view()
    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="20"
        ch_combined_samples = ch_joined_samples
            .combine(ch_intervals)
            .view()
    ```

Schauen wir uns Schritt für Schritt an, was diese map-Operation macht.

Zuerst verwenden wir benannte Parameter, um den Code lesbarer zu machen. Durch die Verwendung der Namen `grouping_key`, `normal`, `tumor` und `interval` können wir auf die Elemente im Tupel nach Namen statt nach Index verweisen:

```groovy
        .map { grouping_key, normal, tumor, interval ->
```

Als Nächstes kombinieren wir den `grouping_key` mit dem `interval`-Feld. Der `grouping_key` ist eine Map mit den Feldern `id` und `repeat`. Wir erstellen eine neue Map mit dem `interval` und fügen sie mit Groovys Map-Addition (`+`) zusammen:

```groovy
                grouping_key + [interval: interval],
```

Schließlich geben wir dies als Tupel mit drei Elementen zurück: die kombinierte Metadaten-Map, die normale Probendatei und die Tumorprobendatei:

```groovy
            [
                grouping_key + [interval: interval],
                normal,
                tumor
            ]
```

Führen wir es erneut aus und prüfen den Kanalinhalt:

```bash
nextflow run main.nf
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [sad_hawking] DSL2 - revision: 1f6f6250cd

    [[id:patientA, repeat:1, interval:chr1], patientA_rep1_normal.bam, patientA_rep1_tumor.bam]
    [[id:patientA, repeat:1, interval:chr2], patientA_rep1_normal.bam, patientA_rep1_tumor.bam]
    [[id:patientA, repeat:1, interval:chr3], patientA_rep1_normal.bam, patientA_rep1_tumor.bam]
    [[id:patientA, repeat:2, interval:chr1], patientA_rep2_normal.bam, patientA_rep2_tumor.bam]
    [[id:patientA, repeat:2, interval:chr2], patientA_rep2_normal.bam, patientA_rep2_tumor.bam]
    [[id:patientA, repeat:2, interval:chr3], patientA_rep2_normal.bam, patientA_rep2_tumor.bam]
    [[id:patientB, repeat:1, interval:chr1], patientB_rep1_normal.bam, patientB_rep1_tumor.bam]
    [[id:patientB, repeat:1, interval:chr2], patientB_rep1_normal.bam, patientB_rep1_tumor.bam]
    [[id:patientB, repeat:1, interval:chr3], patientB_rep1_normal.bam, patientB_rep1_tumor.bam]
    [[id:patientC, repeat:1, interval:chr1], patientC_rep1_normal.bam, patientC_rep1_tumor.bam]
    [[id:patientC, repeat:1, interval:chr2], patientC_rep1_normal.bam, patientC_rep1_tumor.bam]
    [[id:patientC, repeat:1, interval:chr3], patientC_rep1_normal.bam, patientC_rep1_tumor.bam]
    ```

Die Verwendung von `map`, um Daten in die richtige Struktur zu bringen, kann knifflig sein, ist aber entscheidend für eine effektive Datenmanipulation.

Wir haben jetzt jede Probe über alle genomischen Intervalle wiederholt und mehrere unabhängige Analyseeinheiten erstellt, die parallel verarbeitet werden können. Aber was, wenn wir verwandte Proben wieder zusammenführen möchten? Im nächsten Abschnitt lernen wir, wie man Proben gruppiert, die gemeinsame Attribute teilen.

### Fazit

In diesem Abschnitt hast du gelernt:

- **Proben über Intervalle verteilen**: Wie man `combine` verwendet, um Proben über Intervalle zu wiederholen
- **Kartesische Produkte erstellen**: Wie man alle Kombinationen von Proben und Intervallen generiert
- **Kanalstruktur organisieren**: Wie man `map` verwendet, um Daten für bessere Lesbarkeit umzustrukturieren
- **Vorbereitung für die Parallelverarbeitung**: Wie man Daten für die verteilte Analyse einrichtet

## 5. Proben mit `groupTuple` aggregieren

In den vorherigen Abschnitten haben wir gelernt, wie man Daten aus einer Eingabedatei aufteilt und nach bestimmten Feldern filtert (in unserem Fall normale und Tumorproben). Aber das deckt nur eine Art der Zusammenführung ab. Was, wenn wir Proben nach einem bestimmten Attribut gruppieren möchten? Anstatt übereinstimmende Normal-Tumor-Paare zusammenzuführen, möchten wir vielleicht alle Proben von "sampleA" zusammen verarbeiten, unabhängig von ihrem Typ. Dieses Muster ist in Bioinformatik-Workflows üblich, wo du verwandte Proben aus Effizienzgründen getrennt verarbeiten möchtest, bevor du die Ergebnisse am Ende vergleichst oder kombinierst.

Nextflow enthält eingebaute Methoden dafür, die wichtigste, die wir uns ansehen werden, ist `groupTuple`.

Beginnen wir damit, alle unsere Proben zu gruppieren, die dieselben Felder `id` und `interval` haben. Dies wäre typisch für eine Analyse, bei der wir technische Replikate gruppieren, aber bedeutsam unterschiedliche Proben getrennt halten möchten.

Dazu sollten wir unsere Gruppierungsvariablen trennen, damit wir sie isoliert verwenden können.

Der erste Schritt ähnelt dem, was wir im vorherigen Abschnitt getan haben. Wir müssen unsere Gruppierungsvariable als erstes Element des Tupels isolieren. Denk daran: Unser erstes Element ist derzeit eine Map mit den Feldern `id`, `repeat` und `interval`:

```groovy title="main.nf" linenums="1"
{
  "id": "sampleA",
  "repeat": "1",
  "interval": "chr1"
}
```

Wir können die `subMap`-Methode von vorhin wiederverwenden, um unsere Felder `id` und `interval` aus der Map zu isolieren. Wie zuvor verwenden wir den `map`-Operator, um die `subMap`-Methode auf das erste Element des Tupels für jede Probe anzuwenden.

=== "Danach"

    ```groovy title="main.nf" linenums="20" hl_lines="11-19"
        ch_combined_samples = ch_joined_samples
            .combine(ch_intervals)
            .map { grouping_key, normal, tumor, interval ->
                [
                    grouping_key + [interval: interval],
                    normal,
                    tumor
                ]
            }

        ch_grouped_samples = ch_combined_samples
            .map { grouping_key, normal, tumor ->
                [
                    grouping_key.subMap('id', 'interval'),
                    normal,
                    tumor
                ]
              }
              .view()
    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="20" hl_lines="10"
        ch_combined_samples = ch_joined_samples
            .combine(ch_intervals)
            .map { grouping_key, normal, tumor, interval ->
                [
                    grouping_key + [interval: interval],
                    normal,
                    tumor
                ]
            }
            .view()
    ```

Führen wir es erneut aus und prüfen den Kanalinhalt:

```bash
nextflow run main.nf
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [hopeful_brenner] DSL2 - revision: 7f4f7fea76

    [[id:patientA, interval:chr1], patientA_rep1_normal.bam, patientA_rep1_tumor.bam]
    [[id:patientA, interval:chr2], patientA_rep1_normal.bam, patientA_rep1_tumor.bam]
    [[id:patientA, interval:chr3], patientA_rep1_normal.bam, patientA_rep1_tumor.bam]
    [[id:patientA, interval:chr1], patientA_rep2_normal.bam, patientA_rep2_tumor.bam]
    [[id:patientA, interval:chr2], patientA_rep2_normal.bam, patientA_rep2_tumor.bam]
    [[id:patientA, interval:chr3], patientA_rep2_normal.bam, patientA_rep2_tumor.bam]
    [[id:patientB, interval:chr1], patientB_rep1_normal.bam, patientB_rep1_tumor.bam]
    [[id:patientB, interval:chr2], patientB_rep1_normal.bam, patientB_rep1_tumor.bam]
    [[id:patientB, interval:chr3], patientB_rep1_normal.bam, patientB_rep1_tumor.bam]
    [[id:patientC, interval:chr1], patientC_rep1_normal.bam, patientC_rep1_tumor.bam]
    [[id:patientC, interval:chr2], patientC_rep1_normal.bam, patientC_rep1_tumor.bam]
    [[id:patientC, interval:chr3], patientC_rep1_normal.bam, patientC_rep1_tumor.bam]
    ```

Wir können sehen, dass wir die Felder `id` und `interval` erfolgreich isoliert haben, die Proben aber noch nicht gruppiert sind.

!!! note "Hinweis"

    Wir verwerfen hier das `replicate`-Feld. Das liegt daran, dass wir es für die weitere nachgelagerte Verarbeitung nicht benötigen. Versuche nach Abschluss dieses Tutorials, es einzuschließen, ohne die spätere Gruppierung zu beeinflussen!

Gruppieren wir jetzt die Proben nach diesem neuen Gruppierungselement mit dem [`groupTuple`-Operator](https://www.nextflow.io/docs/latest/operator.html#grouptuple).

=== "Danach"

    ```groovy title="main.nf" linenums="30" hl_lines="9"
        ch_grouped_samples = ch_combined_samples
            .map { grouping_key, normal, tumor ->
                [
                    grouping_key.subMap('id', 'interval'),
                    normal,
                    tumor
                ]
              }
              .groupTuple()
              .view()
    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="30"
        ch_grouped_samples = ch_combined_samples
            .map { grouping_key, normal, tumor ->
                [
                    grouping_key.subMap('id', 'interval'),
                    normal,
                    tumor
                ]
              }
              .view()
    ```

Das ist alles! Wir haben nur eine einzige Zeile Code hinzugefügt. Schauen wir, was passiert, wenn wir es ausführen:

```bash
nextflow run main.nf
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [friendly_jang] DSL2 - revision: a1bee1c55d

    [[id:patientA, interval:chr1], [patientA_rep1_normal.bam, patientA_rep2_normal.bam], [patientA_rep1_tumor.bam, patientA_rep2_tumor.bam]]
    [[id:patientA, interval:chr2], [patientA_rep1_normal.bam, patientA_rep2_normal.bam], [patientA_rep1_tumor.bam, patientA_rep2_tumor.bam]]
    [[id:patientA, interval:chr3], [patientA_rep1_normal.bam, patientA_rep2_normal.bam], [patientA_rep1_tumor.bam, patientA_rep2_tumor.bam]]
    [[id:patientB, interval:chr1], [patientB_rep1_normal.bam], [patientB_rep1_tumor.bam]]
    [[id:patientB, interval:chr2], [patientB_rep1_normal.bam], [patientB_rep1_tumor.bam]]
    [[id:patientB, interval:chr3], [patientB_rep1_normal.bam], [patientB_rep1_tumor.bam]]
    [[id:patientC, interval:chr1], [patientC_rep1_normal.bam], [patientC_rep1_tumor.bam]]
    [[id:patientC, interval:chr2], [patientC_rep1_normal.bam], [patientC_rep1_tumor.bam]]
    [[id:patientC, interval:chr3], [patientC_rep1_normal.bam], [patientC_rep1_tumor.bam]]
    ```

Beachte, dass sich unsere Datenstruktur geändert hat und die Dateien innerhalb jedes Kanalelements jetzt in Tupeln wie `[patientA_rep1_normal.bam, patientA_rep2_normal.bam]` enthalten sind. Das liegt daran, dass Nextflow beim Verwenden von `groupTuple` die einzelnen Dateien für jede Probe einer Gruppe kombiniert. Das ist wichtig zu beachten, wenn du die Daten nachgelagert verarbeiten möchtest.

!!! note "Hinweis"

    [`transpose`](https://www.nextflow.io/docs/latest/reference/operator.html#transpose) ist das Gegenteil von groupTuple. Es entpackt die Elemente in einem Kanal und flacht sie ab. Versuche, `transpose` hinzuzufügen und die oben durchgeführte Gruppierung rückgängig zu machen!

### Fazit

In diesem Abschnitt hast du gelernt:

- **Verwandte Proben gruppieren**: Wie man `groupTuple` verwendet, um Proben nach gemeinsamen Attributen zu aggregieren
- **Gruppierungsschlüssel isolieren**: Wie man `subMap` verwendet, um bestimmte Felder für die Gruppierung zu extrahieren
- **Gruppierte Datenstrukturen verarbeiten**: Wie man mit der verschachtelten Struktur arbeitet, die durch `groupTuple` erstellt wird
- **Technische Replikate verarbeiten**: Wie man Proben gruppiert, die dieselben experimentellen Bedingungen teilen

---

## Zusammenfassung

In dieser Side Quest hast du gelernt, wie man Daten mithilfe von Kanälen aufteilt und gruppiert.

Indem du die Daten beim Durchfließen der Pipeline modifizierst, kannst du eine skalierbare Pipeline erstellen, ohne Schleifen oder While-Anweisungen zu verwenden. Das bietet mehrere Vorteile gegenüber traditionelleren Ansätzen:

- Wir können ohne zusätzlichen Code auf beliebig viele oder wenige Eingaben skalieren
- Wir konzentrieren uns auf den Datenfluss durch die Pipeline, statt auf Iteration
- Wir können so komplex oder einfach sein wie nötig
- Die Pipeline wird deklarativer und konzentriert sich darauf, was passieren soll, statt wie es passieren soll
- Nextflow optimiert die Ausführung für uns, indem es unabhängige Operationen parallel ausführt

Die Beherrschung dieser Kanaloperationen ermöglicht es dir, flexible, skalierbare Pipelines zu erstellen, die komplexe Datenbeziehungen ohne Schleifen oder iterative Programmierung verwalten, sodass Nextflow die Ausführung optimieren und unabhängige Operationen automatisch parallelisieren kann.

### Wichtige Muster

1.  **Strukturierte Eingabedaten erstellen:** Ausgehend von einer CSV-Datei mit Meta-Maps (aufbauend auf Mustern aus [Metadata in workflows](../metadata/))

    ```groovy
    ch_samples = channel.fromPath("./data/samplesheet.csv")
        .splitCsv(header: true)
        .map{ row ->
          [[id:row.id, repeat:row.repeat, type:row.type], row.bam]
        }
    ```

2.  **Daten in separate Kanäle aufteilen:** Wir haben `filter` verwendet, um Daten basierend auf dem `type`-Feld in unabhängige Streams aufzuteilen

    ```groovy
    channel.filter { it.type == 'tumor' }
    ```

3.  **Übereinstimmende Proben zusammenführen:** Wir haben `join` verwendet, um verwandte Proben basierend auf den Feldern `id` und `repeat` wieder zusammenzuführen

    - Zwei Kanäle nach Schlüssel (erstes Element des Tupels) zusammenführen

    ```groovy
    tumor_ch.join(normal_ch)
    ```

    - Zusammenführungsschlüssel extrahieren und nach diesem Wert zusammenführen

    ```groovy
    tumor_ch.map { meta, file -> [meta.id, meta, file] }
        .join(
          normal_ch.map { meta, file -> [meta.id, meta, file] }
        )
    ```

    - Nach mehreren Feldern mit subMap zusammenführen

    ```groovy
    tumor_ch.map { meta, file -> [meta.subMap(['id', 'repeat']), meta, file] }
        .join(
          normal_ch.map { meta, file -> [meta.subMap(['id', 'repeat']), meta, file] }
        )
    ```

4.  **Über Intervalle verteilen:** Wir haben `combine` verwendet, um kartesische Produkte von Proben mit genomischen Intervallen für die Parallelverarbeitung zu erstellen.

    ```groovy
    samples_ch.combine(intervals_ch)
    ```

5.  **Nach Gruppierungsschlüsseln aggregieren:** Wir haben `groupTuple` verwendet, um nach dem ersten Element in jedem Tupel zu gruppieren und dabei Proben mit denselben Feldern `id` und `interval` zu sammeln und technische Replikate zusammenzuführen.

    ```groovy
    channel.groupTuple()
    ```

6.  **Datenstruktur optimieren:** Wir haben `subMap` verwendet, um bestimmte Felder zu extrahieren, und eine benannte Closure erstellt, um Transformationen wiederverwendbar zu machen.

    - Bestimmte Felder aus einer Map extrahieren

    ```groovy
    meta.subMap(['id', 'repeat'])
    ```

    - Benannte Closure für wiederverwendbare Transformationen verwenden

    ```groovy
    getSampleIdAndReplicate = { meta, file -> [meta.subMap(['id', 'repeat']), file] }
    channel.map(getSampleIdAndReplicate)
    ```

### Weitere Ressourcen

- [filter](https://www.nextflow.io/docs/latest/operator.html#filter)
- [map](https://www.nextflow.io/docs/latest/operator.html#map)
- [join](https://www.nextflow.io/docs/latest/operator.html#join)
- [groupTuple](https://www.nextflow.io/docs/latest/operator.html#grouptuple)
- [combine](https://www.nextflow.io/docs/latest/operator.html#combine)

---

## Wie geht es weiter?

Kehre zum [Menü der Side Quests](../) zurück oder klicke auf die Schaltfläche unten rechts auf der Seite, um zum nächsten Thema in der Liste zu wechseln.
