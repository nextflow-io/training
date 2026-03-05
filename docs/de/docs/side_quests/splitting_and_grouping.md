# Splitting und Grouping

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Nextflow bietet leistungsstarke Werkzeuge, um flexibel mit Daten zu arbeiten. Eine wichtige Fähigkeit ist das Aufteilen von Daten in verschiedene Streams und das anschließende Gruppieren zusammengehöriger Elemente. Dies ist besonders wertvoll in Bioinformatik-Workflows, wo du verschiedene Arten von Proben separat verarbeiten musst, bevor du die Ergebnisse für die Analyse kombinierst.

Stell es dir wie das Sortieren von Post vor: Du trennst Briefe nach Zielort, verarbeitest jeden Stapel unterschiedlich und kombinierst dann Elemente, die an dieselbe Person gehen. Nextflow verwendet spezielle Operatoren, um dies mit wissenschaftlichen Daten zu erreichen. Dieser Ansatz ist auch allgemein als **scatter/gather**-Muster in verteiltem Computing und Bioinformatik-Workflows bekannt.

Nextflows Channel-System steht im Mittelpunkt dieser Flexibilität. Channels verbinden verschiedene Teile deines Workflows und ermöglichen es den Daten, durch deine Analyse zu fließen. Du kannst mehrere Channels aus einer einzigen Datenquelle erstellen, jeden Channel unterschiedlich verarbeiten und dann Channels bei Bedarf wieder zusammenführen. Dieser Ansatz ermöglicht es dir, Workflows zu entwerfen, die auf natürliche Weise die verzweigenden und konvergierenden Pfade komplexer Bioinformatik-Analysen widerspiegeln.

### Lernziele

In dieser Side Quest lernst du, Daten mit Nextflows Channel-Operatoren aufzuteilen und zu gruppieren.
Wir beginnen mit einer CSV-Datei, die Probeninformationen und zugehörige Datendateien enthält, und manipulieren und reorganisieren diese Daten dann.

Am Ende dieser Side Quest kannst du Datenströme effektiv trennen und kombinieren, indem du die folgenden Techniken verwendest:

- Daten aus Dateien mit `splitCsv` einlesen
- Daten mit `filter` und `map` filtern und transformieren
- Zusammengehörige Daten mit `join` und `groupTuple` kombinieren
- Datenkombinationen mit `combine` für parallele Verarbeitung erstellen
- Datenstruktur mit `subMap` und Deduplizierungsstrategien optimieren
- Wiederverwendbare Funktionen mit benannten Closures erstellen, um dir bei der Manipulation von Channel-Strukturen zu helfen

Diese Fähigkeiten helfen dir, Workflows zu erstellen, die mehrere Eingabedateien und verschiedene Datentypen effizient verarbeiten können, während sie eine saubere, wartbare Code-Struktur beibehalten.

### Voraussetzungen

Bevor du diese Side Quest angehst, solltest du:

- Das [Hello Nextflow](../hello_nextflow/README.md)-Tutorial oder einen gleichwertigen Einsteigerkurs abgeschlossen haben.
- Dich mit grundlegenden Nextflow-Konzepten und -Mechanismen wohl fühlen (Prozesse, Channels, Operatoren, Arbeiten mit Dateien, Meta-Daten)

**Optional:** Wir empfehlen, zuerst die Side Quest [Metadaten in Workflows](./metadata.md) abzuschließen.
Diese behandelt die Grundlagen des Einlesens von CSV-Dateien mit `splitCsv` und der Erstellung von Meta-Maps, die wir hier intensiv nutzen werden.

---

## 0. Einstieg

#### Öffne den Training-Codespace

Falls du dies noch nicht getan hast, stelle sicher, dass du die Trainingsumgebung wie in [Environment Setup](../envsetup/index.md) beschrieben öffnest.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Wechsle in das Projektverzeichnis

Lass uns in das Verzeichnis wechseln, in dem sich die Dateien für dieses Tutorial befinden.

```bash
cd side-quests/splitting_and_grouping
```

Du kannst VSCode so einstellen, dass es sich auf dieses Verzeichnis konzentriert:

```bash
code .
```

#### Überprüfe die Materialien

Du findest eine Haupt-Workflow-Datei und ein `data`-Verzeichnis, das ein Samplesheet namens `samplesheet.csv` enthält.

```console title="Verzeichnisinhalt"
.
├── data
│   └── samplesheet.csv
└── main.nf
```

Das Samplesheet enthält Informationen über Proben von verschiedenen Patienten, einschließlich der Patienten-ID, der Probenwiederholungsnummer, des Typs (normal oder tumor) und Pfaden zu hypothetischen Datendateien (die nicht tatsächlich existieren, aber wir tun so, als ob).

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

Dieses Samplesheet listet acht Proben von drei Patienten (A, B, C) auf.

Für jeden Patienten haben wir Proben vom Typ `tumor` (typischerweise aus Tumorbiopsien stammend) oder `normal` (aus gesundem Gewebe oder Blut entnommen).
Falls du mit Krebs-Analysen nicht vertraut bist, solltest du nur wissen, dass dies einem experimentellen Modell entspricht, das gepaarte Tumor/Normal-Proben verwendet, um kontrastive Analysen durchzuführen.

Speziell für Patient A haben wir zwei Sätze technischer Replikate (Wiederholungen).

!!! note

    Mach dir keine Sorgen, wenn du mit diesem experimentellen Design nicht vertraut bist, es ist nicht entscheidend für das Verständnis dieses Tutorials.

#### Überprüfe die Aufgabe

Deine Herausforderung besteht darin, einen Nextflow-Workflow zu schreiben, der:

1. Probendaten aus einer CSV-Datei **einliest** und mit Meta-Maps strukturiert
2. Proben basierend auf dem Typ (normal vs. tumor) in verschiedene Channels **aufteilt**
3. Übereinstimmende Tumor/Normal-Paare nach Patienten-ID und Replikatnummer **verbindet**
4. Proben über genomische Intervalle für parallele Verarbeitung **verteilt**
5. Zusammengehörige Proben für die nachgelagerte Analyse wieder **gruppiert**

Dies stellt ein häufiges Bioinformatik-Muster dar, bei dem du Daten für unabhängige Verarbeitung aufteilen und dann zusammengehörige Elemente für vergleichende Analysen wieder kombinieren musst.

#### Bereitschafts-Checkliste

Glaubst du, du bist bereit loszulegen?

- [ ] Ich verstehe das Ziel dieses Kurses und seine Voraussetzungen
- [ ] Mein Codespace läuft
- [ ] Ich habe mein Arbeitsverzeichnis entsprechend eingestellt
- [ ] Ich verstehe die Aufgabe

Wenn du alle Kästchen ankreuzen kannst, kann es losgehen.

---

## 1. Probendaten einlesen

### 1.1. Probendaten mit `splitCsv` einlesen und Meta-Maps erstellen

Lass uns damit beginnen, die Probendaten mit `splitCsv` einzulesen und sie im Meta-Map-Muster zu organisieren. In der `main.nf` siehst du, dass wir bereits mit dem Workflow begonnen haben.

```groovy title="main.nf" linenums="1" hl_lines="2"
workflow {
    ch_samplesheet = channel.fromPath("./data/samplesheet.csv")
}
```

!!! note

    In diesem Tutorial werden wir das Präfix `ch_` für alle Channel-Variablen verwenden, um klar anzuzeigen, dass sie Nextflow-Channels sind.

Wenn du die Side Quest [Metadaten in Workflows](./metadata.md) abgeschlossen hast, wirst du dieses Muster wiedererkennen. Wir verwenden `splitCsv`, um die CSV zu lesen und strukturieren die Daten sofort mit einer Meta-Map, um Metadaten von Dateipfaden zu trennen.

!!! info

    In diesem Training werden wir zwei verschiedenen Konzepten begegnen, die `map` genannt werden:

    - **Datenstruktur**: Die Groovy-Map (entspricht Dictionaries/Hashes in anderen Sprachen), die Schlüssel-Wert-Paare speichert
    - **Channel-Operator**: Der `.map()`-Operator, der Elemente in einem Channel transformiert

    Wir werden im Kontext klären, welches wir meinen, aber diese Unterscheidung ist wichtig für die Arbeit mit Nextflow.

Wende diese Änderungen auf `main.nf` an:

=== "Nachher"

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

Dies kombiniert die `splitCsv`-Operation (Lesen der CSV mit Headern) und die `map`-Operation (Strukturierung der Daten als `[meta, file]`-Tupel) in einem Schritt. Wende diese Änderung an und führe die Pipeline aus:

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

Wir haben jetzt einen Channel, bei dem jedes Element ein `[meta, file]`-Tupel ist - Metadaten getrennt von Dateipfaden. Diese Struktur ermöglicht es uns, unsere Arbeitslast basierend auf Metadatenfeldern aufzuteilen und zu gruppieren.

---

## 2. Daten filtern und transformieren

### 2.1. Daten mit `filter` filtern

Wir können den [`filter`-Operator](https://www.nextflow.io/docs/latest/operator.html#filter) verwenden, um die Daten basierend auf einer Bedingung zu filtern. Nehmen wir an, wir wollen nur normale Proben verarbeiten. Wir können dies tun, indem wir die Daten basierend auf dem `type`-Feld filtern. Lass uns dies vor dem `view`-Operator einfügen.

=== "Nachher"

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

Wir haben die Daten erfolgreich gefiltert, um nur normale Proben einzuschließen. Lass uns zusammenfassen, wie dies funktioniert.

Der `filter`-Operator nimmt eine Closure, die auf jedes Element im Channel angewendet wird. Wenn die Closure `true` zurückgibt, wird das Element eingeschlossen; wenn sie `false` zurückgibt, wird das Element ausgeschlossen.

In unserem Fall wollen wir nur Proben behalten, bei denen `meta.type == 'normal'`. Die Closure verwendet das Tupel `meta,file`, um auf jede Probe zu verweisen, greift auf den Probentyp mit `meta.type` zu und prüft, ob er `'normal'` entspricht.

Dies wird mit der einzelnen Closure erreicht, die wir oben eingeführt haben:

```groovy title="main.nf" linenums="7"
    .filter { meta, file -> meta.type == 'normal' }
```

### 2.2. Separate gefilterte Channels erstellen

Derzeit wenden wir den Filter auf den direkt aus der CSV erstellten Channel an, aber wir wollen dies auf mehr als eine Weise filtern, also lass uns die Logik umschreiben, um einen separaten gefilterten Channel für normale Proben zu erstellen:

=== "Nachher"

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

Wir haben die Daten erfolgreich gefiltert und einen separaten Channel für normale Proben erstellt.

Lass uns auch einen gefilterten Channel für die Tumorproben erstellen:

=== "Nachher"

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

Wir haben die normalen und Tumorproben in zwei verschiedene Channels aufgeteilt und eine Closure, die an `view()` übergeben wird, verwendet, um sie in der Ausgabe unterschiedlich zu kennzeichnen: `ch_tumor_samples.view{'Tumor sample: ' + it}`.

### Zusammenfassung

In diesem Abschnitt hast du gelernt:

- **Daten filtern**: Wie man Daten mit `filter` filtert
- **Daten aufteilen**: Wie man Daten basierend auf einer Bedingung in verschiedene Channels aufteilt
- **Daten anzeigen**: Wie man `view` verwendet, um die Daten auszugeben und die Ausgabe verschiedener Channels zu kennzeichnen

Wir haben jetzt die normalen und Tumorproben in zwei verschiedene Channels aufgeteilt. Als nächstes werden wir die normalen und Tumorproben basierend auf dem `id`-Feld verbinden.

---

## 3. Channels nach Identifikatoren verbinden

Im vorherigen Abschnitt haben wir die normalen und Tumorproben in zwei verschiedene Channels aufgeteilt. Diese könnten unabhängig voneinander mit spezifischen Prozessen oder Workflows basierend auf ihrem Typ verarbeitet werden. Aber was passiert, wenn wir die normalen und Tumorproben desselben Patienten vergleichen wollen? An diesem Punkt müssen wir sie wieder zusammenführen und sicherstellen, dass wir die Proben basierend auf ihrem `id`-Feld abgleichen.

Nextflow enthält viele Methoden zum Kombinieren von Channels, aber in diesem Fall ist der am besten geeignete Operator [`join`](https://www.nextflow.io/docs/latest/operator.html#join). Wenn du mit SQL vertraut bist, verhält er sich wie die `JOIN`-Operation, bei der wir den Schlüssel zum Verbinden und den Typ des Joins angeben.

### 3.1. `map` und `join` verwenden, um basierend auf Patienten-ID zu kombinieren

Wenn wir die [`join`-Dokumentation](https://www.nextflow.io/docs/latest/operator.html#join) überprüfen, sehen wir, dass er standardmäßig zwei Channels basierend auf dem ersten Element in jedem Tupel verbindet.

#### 3.1.1. Datenstruktur überprüfen

Wenn du die Konsolenausgabe nicht mehr verfügbar hast, lass uns die Pipeline ausführen, um unsere Datenstruktur zu überprüfen und zu sehen, wie wir sie ändern müssen, um auf dem `id`-Feld zu verbinden.

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

Wir können sehen, dass das `id`-Feld das erste Element in jeder Meta-Map ist. Damit `join` funktioniert, sollten wir das `id`-Feld in jedem Tupel isolieren. Danach können wir einfach den `join`-Operator verwenden, um die beiden Channels zu kombinieren.

#### 3.1.2. Das `id`-Feld isolieren

Um das `id`-Feld zu isolieren, können wir den [`map`-Operator](https://www.nextflow.io/docs/latest/operator.html#map) verwenden, um ein neues Tupel mit dem `id`-Feld als erstes Element zu erstellen.

=== "Nachher"

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

Es mag subtil sein, aber du solltest sehen können, dass das erste Element in jedem Tupel das `id`-Feld ist.

#### 3.1.3. Die beiden Channels kombinieren

Jetzt können wir den `join`-Operator verwenden, um die beiden Channels basierend auf dem `id`-Feld zu kombinieren.

Noch einmal werden wir `view` verwenden, um die verbundenen Ausgaben auszugeben.

=== "Nachher"

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

Es ist etwas schwer zu erkennen, weil es so breit ist, aber du solltest sehen können, dass die Proben nach dem `id`-Feld verbunden wurden. Jedes Tupel hat jetzt das Format:

- `id`: Die Proben-ID
- `normal_meta_map`: Die Metadaten der normalen Probe einschließlich Typ, Replikat und Pfad zur BAM-Datei
- `normal_sample_file`: Die normale Probendatei
- `tumor_meta_map`: Die Metadaten der Tumorprobe einschließlich Typ, Replikat und Pfad zur BAM-Datei
- `tumor_sample`: Die Tumorprobe einschließlich Typ, Replikat und Pfad zur BAM-Datei

!!! warning

    Der `join`-Operator verwirft alle nicht übereinstimmenden Tupel. In diesem Beispiel haben wir sichergestellt, dass alle Proben für Tumor und Normal übereinstimmen, aber wenn dies nicht wahr ist, musst du den Parameter `remainder: true` verwenden, um die nicht übereinstimmenden Tupel zu behalten. Überprüfe die [Dokumentation](https://www.nextflow.io/docs/latest/operator.html#join) für weitere Details.

Jetzt weißt du, wie man `map` verwendet, um ein Feld in einem Tupel zu isolieren, und wie man `join` verwendet, um Tupel basierend auf dem ersten Feld zu kombinieren.
Mit diesem Wissen können wir erfolgreich Channels basierend auf einem gemeinsamen Feld kombinieren.

Als nächstes betrachten wir die Situation, in der du auf mehreren Feldern verbinden möchtest.

### 3.2. Auf mehreren Feldern verbinden

Wir haben 2 Replikate für sampleA, aber nur 1 für sampleB und sampleC. In diesem Fall konnten wir sie effektiv verbinden, indem wir das `id`-Feld verwendeten, aber was würde passieren, wenn sie nicht synchron wären? Wir könnten die normalen und Tumorproben verschiedener Replikate durcheinanderbringen!

Um dies zu vermeiden, können wir auf mehreren Feldern verbinden. Es gibt tatsächlich mehrere Möglichkeiten, dies zu erreichen, aber wir werden uns darauf konzentrieren, einen neuen Verbindungsschlüssel zu erstellen, der sowohl die Proben-`id` als auch die `replicate`-Nummer enthält.

Beginnen wir damit, einen neuen Verbindungsschlüssel zu erstellen. Wir können dies auf die gleiche Weise wie zuvor tun, indem wir den [`map`-Operator](https://www.nextflow.io/docs/latest/operator.html#map) verwenden, um ein neues Tupel mit den `id`- und `repeat`-Feldern als erstes Element zu erstellen.

=== "Nachher"

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

Jetzt sollten wir sehen, dass die Verbindung erfolgt, aber sowohl die `id`- als auch die `repeat`-Felder verwendet werden. Führe den Workflow aus:

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

Beachte, wie wir ein Tupel mit zwei Elementen (`id`- und `repeat`-Felder) als erstes Element jedes verbundenen Ergebnisses haben. Dies zeigt, wie komplexe Elemente als Verbindungsschlüssel verwendet werden können, was ziemlich komplizierte Übereinstimmungen zwischen Proben aus denselben Bedingungen ermöglicht.

Wenn du mehr Möglichkeiten erkunden möchtest, auf verschiedenen Schlüsseln zu verbinden, schau dir die [join-Operator-Dokumentation](https://www.nextflow.io/docs/latest/operator.html#join) für zusätzliche Optionen und Beispiele an.

### 3.3. `subMap` verwenden, um einen neuen Verbindungsschlüssel zu erstellen

Der vorherige Ansatz verliert die Feldnamen aus unserem Verbindungsschlüssel - die `id`- und `repeat`-Felder werden nur zu einer Liste von Werten. Um die Feldnamen für späteren Zugriff zu behalten, können wir die [`subMap`-Methode](<https://docs.groovy-lang.org/latest/html/groovy-jdk/java/util/Map.html#subMap(java.util.Collection)>) verwenden.

Die `subMap`-Methode extrahiert nur die angegebenen Schlüssel-Wert-Paare aus einer Map. Hier extrahieren wir nur die `id`- und `repeat`-Felder, um unseren Verbindungsschlüssel zu erstellen.

=== "Nachher"

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

Jetzt haben wir einen neuen Verbindungsschlüssel, der nicht nur die `id`- und `repeat`-Felder enthält, sondern auch die Feldnamen beibehält, sodass wir später über den Namen auf sie zugreifen können, z.B. `meta.id` und `meta.repeat`.

### 3.4. Eine benannte Closure in map verwenden

Um Duplikation zu vermeiden und Fehler zu reduzieren, können wir eine benannte Closure verwenden. Eine benannte Closure ermöglicht es uns, eine wiederverwendbare Funktion zu erstellen, die wir an mehreren Stellen aufrufen können.

Dazu definieren wir zuerst die Closure als neue Variable:

=== "Nachher"

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

Beachte, dass wir den Dateipfad auch mit `file()` in ein Path-Objekt konvertieren, sodass jeder Prozess, der diesen Channel empfängt, die Datei korrekt verarbeiten kann (für weitere Informationen siehe [Arbeiten mit Dateien](./working_with_files.md)).

Lass uns die Closure in unserem Workflow implementieren:

=== "Nachher"

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

!!! note

    Der `map`-Operator hat von der Verwendung von `{ }` zur Verwendung von `( )` gewechselt, um die Closure als Argument zu übergeben. Dies liegt daran, dass der `map`-Operator eine Closure als Argument erwartet und `{ }` verwendet wird, um eine anonyme Closure zu definieren. Beim Aufrufen einer benannten Closure verwende die `( )`-Syntax.

Führe den Workflow noch einmal aus, um zu überprüfen, ob alles noch funktioniert:

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

Die Verwendung einer benannten Closure ermöglicht es uns, dieselbe Transformation an mehreren Stellen wiederzuverwenden, das Fehlerrisiko zu reduzieren und den Code lesbarer und wartbarer zu machen.

### 3.5. Datenduplikation reduzieren

Wir haben viele duplizierte Daten in unserem Workflow. Jedes Element in den verbundenen Proben wiederholt die `id`- und `repeat`-Felder. Da diese Informationen bereits im Gruppierungsschlüssel verfügbar sind, können wir diese Redundanz vermeiden. Zur Erinnerung: Unsere aktuelle Datenstruktur sieht so aus:

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

Da die `id`- und `repeat`-Felder im Gruppierungsschlüssel verfügbar sind, lass uns sie aus dem Rest jedes Channel-Elements entfernen, um Duplikation zu vermeiden. Wir können dies tun, indem wir die `subMap`-Methode verwenden, um eine neue Map nur mit dem `type`-Feld zu erstellen. Dieser Ansatz ermöglicht es uns, alle notwendigen Informationen beizubehalten und gleichzeitig Redundanz in unserer Datenstruktur zu eliminieren.

=== "Nachher"

    ```groovy title="main.nf" linenums="8" hl_lines="1"
        getSampleIdAndReplicate = { meta, bam -> [ meta.subMap(['id', 'repeat']), meta.subMap(['type']), file(bam) ] }
    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="8" hl_lines="1"
        getSampleIdAndReplicate = { meta, bam -> [ meta.subMap(['id', 'repeat']), meta, file(bam) ] }
    ```

Jetzt gibt die Closure ein Tupel zurück, bei dem das erste Element die `id`- und `repeat`-Felder enthält und das zweite Element nur das `type`-Feld enthält. Dies eliminiert Redundanz, indem die `id`- und `repeat`-Informationen einmal im Gruppierungsschlüssel gespeichert werden, während alle notwendigen Informationen erhalten bleiben.

Führe den Workflow aus, um zu sehen, wie dies aussieht:

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

Wir können sehen, dass wir die `id`- und `repeat`-Felder nur einmal im Gruppierungsschlüssel angeben und das `type`-Feld in den Probendaten haben. Wir haben keine Informationen verloren, aber wir haben es geschafft, unsere Channel-Inhalte prägnanter zu machen.

### 3.6. Redundante Informationen entfernen

Wir haben oben duplizierte Informationen entfernt, aber wir haben immer noch einige andere redundante Informationen in unseren Channels.

Am Anfang haben wir die normalen und Tumorproben mit `filter` getrennt und sie dann basierend auf `id`- und `repeat`-Schlüsseln verbunden. Der `join`-Operator bewahrt die Reihenfolge, in der Tupel zusammengeführt werden. In unserem Fall, mit normalen Proben auf der linken Seite und Tumorproben auf der rechten, behält der resultierende Channel diese Struktur bei: `id, <normale Elemente>, <Tumor-Elemente>`.

Da wir die Position jedes Elements in unserem Channel kennen, können wir die Struktur weiter vereinfachen, indem wir die Metadaten `[type:normal]` und `[type:tumor]` weglassen.

=== "Nachher"

    ```groovy title="main.nf" linenums="8" hl_lines="1"
        getSampleIdAndReplicate = { meta, file -> [ meta.subMap(['id', 'repeat']), file ] }
    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="8" hl_lines="1"
        getSampleIdAndReplicate = { meta, file -> [ meta.subMap(['id', 'repeat']), meta.subMap(['type']), file ] }
    ```

Führe erneut aus, um das Ergebnis zu sehen:

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

### Zusammenfassung

In diesem Abschnitt hast du gelernt:

- **Tupel manipulieren**: Wie man `map` verwendet, um ein Feld in einem Tupel zu isolieren
- **Tupel verbinden**: Wie man `join` verwendet, um Tupel basierend auf dem ersten Feld zu kombinieren
- **Verbindungsschlüssel erstellen**: Wie man `subMap` verwendet, um einen neuen Verbindungsschlüssel zu erstellen
- **Benannte Closures**: Wie man eine benannte Closure in map verwendet
- **Verbinden auf mehreren Feldern**: Wie man auf mehreren Feldern verbindet für präzisere Übereinstimmung
- **Datenstruktur-Optimierung**: Wie man die Channel-Struktur durch Entfernen redundanter Informationen optimiert

Du hast jetzt einen Workflow, der ein Samplesheet aufteilen, die normalen und Tumorproben filtern, sie nach Proben-ID und Replikatnummer verbinden und dann die Ergebnisse ausgeben kann.

Dies ist ein häufiges Muster in Bioinformatik-Workflows, bei dem du Proben oder andere Arten von Daten nach unabhängiger Verarbeitung abgleichen musst, also ist es eine nützliche Fähigkeit. Als nächstes werden wir uns ansehen, wie man eine Probe mehrmals wiederholt.

## 4. Proben über Intervalle verteilen

Ein wichtiges Muster in Bioinformatik-Workflows ist die Verteilung der Analyse über genomische Regionen. Zum Beispiel kann Variant Calling parallelisiert werden, indem das Genom in Intervalle (wie Chromosomen oder kleinere Regionen) aufgeteilt wird. Diese Parallelisierungsstrategie verbessert die Pipeline-Effizienz erheblich, indem die Rechenlast auf mehrere Cores oder Nodes verteilt wird, was die Gesamtausführungszeit reduziert.

Im folgenden Abschnitt demonstrieren wir, wie du deine Probendaten über mehrere genomische Intervalle verteilst. Wir paaren jede Probe mit jedem Intervall, was die parallele Verarbeitung verschiedener genomischer Regionen ermöglicht. Dies multipliziert unsere Datensatzgröße mit der Anzahl der Intervalle und erstellt mehrere unabhängige Analyseeinheiten, die später wieder zusammengeführt werden können.

### 4.1. Proben mit `combine` über Intervalle verteilen

Beginnen wir damit, einen Channel von Intervallen zu erstellen. Um das Leben einfach zu halten, verwenden wir einfach 3 Intervalle, die wir manuell definieren. In einem echten Workflow könntest du diese aus einer Dateieingabe einlesen oder sogar einen Channel mit vielen Intervalldateien erstellen.

=== "Nachher"

    ```groovy title="main.nf" linenums="17" hl_lines="2"
            .join(ch_tumor_samples)
        ch_intervals = channel.of('chr1', 'chr2', 'chr3')
    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="17" hl_lines="2"
            .join(ch_tumor_samples)
        ch_joined_samples.view()
    ```

Denk daran, wir wollen jede Probe für jedes Intervall wiederholen. Dies wird manchmal als kartesisches Produkt der Proben und Intervalle bezeichnet. Wir können dies erreichen, indem wir den [`combine`-Operator](https://www.nextflow.io/docs/latest/operator.html#combine) verwenden. Dies nimmt jedes Element aus Channel 1 und wiederholt es für jedes Element in Channel 2. Lass uns einen combine-Operator zu unserem Workflow hinzufügen:

=== "Nachher"

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

Lass uns es jetzt ausführen und sehen, was passiert:

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

Erfolg! Wir haben jede Probe für jedes einzelne Intervall in unserer 3-Intervall-Liste wiederholt. Wir haben die Anzahl der Elemente in unserem Channel effektiv verdreifacht.

Es ist jedoch etwas schwer zu lesen, also werden wir es im nächsten Abschnitt aufräumen.

### 4.2. Den Channel organisieren

Wir können den `map`-Operator verwenden, um unsere Probendaten aufzuräumen und zu refaktorieren, damit sie leichter verständlich sind. Lass uns die Intervall-Zeichenkette zur Verbindungs-Map am ersten Element verschieben.

=== "Nachher"

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

Lass uns Schritt für Schritt aufschlüsseln, was diese Map-Operation macht.

Zuerst verwenden wir benannte Parameter, um den Code lesbarer zu machen. Durch die Verwendung der Namen `grouping_key`, `normal`, `tumor` und `interval` können wir auf die Elemente im Tupel nach Namen statt nach Index verweisen:

```groovy
        .map { grouping_key, normal, tumor, interval ->
```

Als nächstes kombinieren wir den `grouping_key` mit dem `interval`-Feld. Der `grouping_key` ist eine Map, die `id`- und `repeat`-Felder enthält. Wir erstellen eine neue Map mit dem `interval` und führen sie mit Groovys Map-Addition (`+`) zusammen:

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

Lass uns es erneut ausführen und die Channel-Inhalte überprüfen:

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

Die Verwendung von `map`, um deine Daten in die korrekte Struktur zu bringen, kann schwierig sein, aber es ist entscheidend für effektive Datenmanipulation.

Wir haben jetzt jede Probe über alle genomischen Intervalle hinweg wiederholt, was mehrere unabhängige Analyseeinheiten erstellt, die parallel verarbeitet werden können. Aber was, wenn wir zusammengehörige Proben wieder zusammenbringen wollen? Im nächsten Abschnitt lernen wir, wie man Proben gruppiert, die gemeinsame Attribute teilen.

### Zusammenfassung

In diesem Abschnitt hast du gelernt:

- **Proben über Intervalle verteilen**: Wie man `combine` verwendet, um Proben über Intervalle zu wiederholen
- **Kartesische Produkte erstellen**: Wie man alle Kombinationen von Proben und Intervallen generiert
- **Channel-Struktur organisieren**: Wie man `map` verwendet, um Daten für bessere Lesbarkeit umzustrukturieren
- **Vorbereitung für parallele Verarbeitung**: Wie man Daten für verteilte Analyse vorbereitet

## 5. Proben mit `groupTuple` aggregieren

In den vorherigen Abschnitten haben wir gelernt, wie man Daten aus einer Eingabedatei aufteilt und nach bestimmten Feldern filtert (in unserem Fall normale und Tumorproben). Aber das deckt nur eine Art der Verbindung ab. Was, wenn wir Proben nach einem bestimmten Attribut gruppieren wollen? Zum Beispiel, anstatt übereinstimmende Normal-Tumor-Paare zu verbinden, möchten wir vielleicht alle Proben von "sampleA" zusammen verarbeiten, unabhängig von ihrem Typ. Dieses Muster ist häufig in Bioinformatik-Workflows, in denen man zusammengehörige Proben aus Effizienzgründen getrennt verarbeiten möchte, bevor man die Ergebnisse am Ende vergleicht oder kombiniert.

Nextflow enthält eingebaute Methoden dafür, die wichtigste, die wir uns ansehen werden, ist `groupTuple`.

Beginnen wir damit, alle unsere Proben zu gruppieren, die die gleichen `id`- und `interval`-Felder haben. Dies wäre typisch für eine Analyse, bei der wir technische Replikate gruppieren, aber bedeutsam verschiedene Proben getrennt halten wollen.

Dazu sollten wir unsere Gruppierungsvariablen separieren, damit wir sie isoliert verwenden können.

Der erste Schritt ist ähnlich zu dem, was wir im vorherigen Abschnitt getan haben. Wir müssen unsere Gruppierungsvariable als erstes Element des Tupels isolieren. Denk daran, unser erstes Element ist derzeit eine Map der Felder `id`, `repeat` und `interval`:

```groovy title="main.nf" linenums="1"
{
  "id": "sampleA",
  "repeat": "1",
  "interval": "chr1"
}
```

Wir können die `subMap`-Methode von vorher wiederverwenden, um unsere `id`- und `interval`-Felder aus der Map zu isolieren. Wie zuvor verwenden wir den `map`-Operator, um die `subMap`-Methode auf das erste Element des Tupels für jede Probe anzuwenden.

=== "Nachher"

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

Lass uns es erneut ausführen und die Channel-Inhalte überprüfen:

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

Wir können sehen, dass wir die `id`- und `interval`-Felder erfolgreich isoliert haben, aber die Proben noch nicht gruppiert.

!!! note

    Wir verwerfen hier das `replicate`-Feld. Dies liegt daran, dass wir es für die weitere nachgelagerte Verarbeitung nicht benötigen. Nachdem du dieses Tutorial abgeschlossen hast, versuche, es einzuschließen, ohne die spätere Gruppierung zu beeinflussen!

Lass uns nun die Proben mit dem [`groupTuple`-Operator](https://www.nextflow.io/docs/latest/operator.html#grouptuple) nach diesem neuen Gruppierungselement gruppieren.

=== "Nachher"

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

Das ist alles! Wir haben nur eine einzige Zeile Code hinzugefügt. Lass uns sehen, was passiert, wenn wir es ausführen:

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

Beachte, dass sich unsere Datenstruktur geändert hat und die Dateien innerhalb jedes Channel-Elements jetzt in Tupeln wie `[patientA_rep1_normal.bam, patientA_rep2_normal.bam]` enthalten sind. Dies liegt daran, dass `groupTuple` die einzelnen Dateien für jede Probe einer Gruppe kombiniert, wenn wir es verwenden. Dies ist wichtig zu beachten, wenn man die Daten nachgelagert verarbeiten will.

!!! note

    [`transpose`](https://www.nextflow.io/docs/latest/reference/operator.html#transpose) ist das Gegenteil von groupTuple. Es entpackt die Elemente in einem Channel und flacht sie ab. Versuche, `transpose` hinzuzufügen und die Gruppierung rückgängig zu machen, die wir oben durchgeführt haben!

### Zusammenfassung

In diesem Abschnitt hast du gelernt:

- **Zusammengehörige Proben gruppieren**: Wie man `groupTuple` verwendet, um Proben nach gemeinsamen Attributen zu aggregieren
- **Gruppierungsschlüssel isolieren**: Wie man `subMap` verwendet, um bestimmte Felder für die Gruppierung zu extrahieren
- **Gruppierte Datenstrukturen handhaben**: Wie man mit der verschachtelten Struktur arbeitet, die von `groupTuple` erstellt wird
- **Technische Replikate handhaben**: Wie man Proben gruppiert, die die gleichen experimentellen Bedingungen teilen

---

## Zusammenfassung

In dieser Side Quest hast du gelernt, wie man Daten mithilfe von Channels aufteilt und gruppiert.

Indem du die Daten modifizierst, während sie durch die Pipeline fließen, kannst du eine skalierbare Pipeline ohne Schleifen oder While-Anweisungen erstellen, die mehrere Vorteile gegenüber traditionelleren Ansätzen bietet:

- Wir können auf beliebig viele oder wenige Eingaben skalieren, ohne zusätzlichen Code
- Wir konzentrieren uns auf den Datenfluss durch die Pipeline statt auf Iteration
- Wir können so komplex oder einfach wie nötig sein
- Die Pipeline wird deklarativer und konzentriert sich darauf, was passieren soll, statt darauf, wie es passieren soll
- Nextflow optimiert die Ausführung für uns, indem unabhängige Operationen parallel ausgeführt werden

Die Beherrschung dieser Channel-Operationen ermöglicht es dir, flexible, skalierbare Pipelines zu erstellen, die komplexe Datenbeziehungen handhaben, ohne auf Schleifen oder iterative Programmierung zurückzugreifen, wobei Nextflow die Ausführung optimiert und unabhängige Operationen automatisch parallelisiert.

### Wichtige Muster

1.  **Strukturierte Eingabedaten erstellen:** Ausgehend von einer CSV-Datei mit Meta-Maps (aufbauend auf Mustern aus [Metadaten in Workflows](./metadata.md))

    ```groovy
    ch_samples = channel.fromPath("./data/samplesheet.csv")
        .splitCsv(header: true)
        .map{ row ->
          [[id:row.id, repeat:row.repeat, type:row.type], row.bam]
        }
    ```

2.  **Daten in separate Channels aufteilen:** Wir haben `filter` verwendet, um Daten basierend auf dem `type`-Feld in unabhängige Streams aufzuteilen

    ```groovy
    channel.filter { it.type == 'tumor' }
    ```

3.  **Übereinstimmende Proben verbinden:** Wir haben `join` verwendet, um zusammengehörige Proben basierend auf `id`- und `repeat`-Feldern wieder zu kombinieren

    - Zwei Channels nach Schlüssel verbinden (erstes Element des Tupels)

    ```groovy
    tumor_ch.join(normal_ch)
    ```

    - Verbindungsschlüssel extrahieren und nach diesem Wert verbinden

    ```groovy
    tumor_ch.map { meta, file -> [meta.id, meta, file] }
        .join(
          normal_ch.map { meta, file -> [meta.id, meta, file] }
        )
    ```

    - Auf mehreren Feldern mit subMap verbinden

    ```groovy
    tumor_ch.map { meta, file -> [meta.subMap(['id', 'repeat']), meta, file] }
        .join(
          normal_ch.map { meta, file -> [meta.subMap(['id', 'repeat']), meta, file] }
        )
    ```

4.  **Über Intervalle verteilen:** Wir haben `combine` verwendet, um kartesische Produkte von Proben mit genomischen Intervallen für parallele Verarbeitung zu erstellen.

    ```groovy
    samples_ch.combine(intervals_ch)
    ```

5.  **Nach Gruppierungsschlüsseln aggregieren:** Wir haben `groupTuple` verwendet, um nach dem ersten Element in jedem Tupel zu gruppieren, wodurch Proben gesammelt wurden, die `id`- und `interval`-Felder teilen und technische Replikate zusammengeführt werden.

    ```groovy
    channel.groupTuple()
    ```

6.  **Datenstruktur optimieren:** Wir haben `subMap` verwendet, um bestimmte Felder zu extrahieren und eine benannte Closure erstellt, um Transformationen wiederverwendbar zu machen.

    - Bestimmte Felder aus einer Map extrahieren

    ```groovy
    meta.subMap(['id', 'repeat'])
    ```

    - Benannte Closure für wiederverwendbare Transformationen verwenden

    ```groovy
    getSampleIdAndReplicate = { meta, file -> [meta.subMap(['id', 'repeat']), file] }
    channel.map(getSampleIdAndReplicate)
    ```

### Zusätzliche Ressourcen

- [filter](https://www.nextflow.io/docs/latest/operator.html#filter)
- [map](https://www.nextflow.io/docs/latest/operator.html#map)
- [join](https://www.nextflow.io/docs/latest/operator.html#join)
- [groupTuple](https://www.nextflow.io/docs/latest/operator.html#grouptuple)
- [combine](https://www.nextflow.io/docs/latest/operator.html#combine)

---

## Wie geht es weiter?

Kehre zum [Menü der Side Quests](./index.md) zurück oder klicke auf die Schaltfläche unten rechts auf der Seite, um zum nächsten Thema in der Liste zu gelangen.
