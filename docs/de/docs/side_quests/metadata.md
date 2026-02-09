# Metadaten und Meta-Maps

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Bei jeder wissenschaftlichen Analyse arbeiten wir selten nur mit den reinen Datendateien.
Jede Datei kommt mit ihren eigenen zusätzlichen Informationen: was sie ist, woher sie stammt und was sie besonders macht.
Diese zusätzlichen Informationen nennen wir Metadaten.

Metadaten sind Daten, die andere Daten beschreiben.
Metadaten verfolgen wichtige Details über Dateien und experimentelle Bedingungen und helfen dabei, Analysen auf die einzigartigen Eigenschaften jedes Datensatzes zuzuschneiden.

Stell dir das wie einen Bibliothekskatalog vor: Während Bücher den eigentlichen Inhalt enthalten (Rohdaten), liefern die Katalogkarten wesentliche Informationen über jedes Buch – wann es veröffentlicht wurde, wer es geschrieben hat, wo man es findet (Metadaten).
In Nextflow-Pipelines können Metadaten verwendet werden, um:

- Dateispezifische Informationen durch den gesamten Workflow hindurch zu verfolgen
- Prozesse basierend auf Dateimerkmalen zu konfigurieren
- Zusammengehörige Dateien für gemeinsame Analysen zu gruppieren

### Lernziele

In dieser Side Quest werden wir erkunden, wie man Metadaten in Workflows handhabt.
Ausgehend von einem einfachen Datenblatt (in der Bioinformatik oft als Samplesheet bezeichnet), das grundlegende Dateiinformationen enthält, wirst du lernen, wie man:

- Dateimetadaten aus CSV-Dateien liest und verarbeitet
- Metadaten-Maps erstellt und manipuliert
- Neue Metadatenfelder während der Workflow-Ausführung hinzufügt
- Metadaten verwendet, um das Prozessverhalten anzupassen

Diese Fähigkeiten helfen dir, robustere und flexiblere Pipelines zu erstellen, die komplexe Dateibeziehungen und Verarbeitungsanforderungen handhaben können.

### Voraussetzungen

Bevor du diese Side Quest angehst, solltest du:

- Das [Hello Nextflow](../hello_nextflow/README.md)-Tutorial oder einen gleichwertigen Einsteigerkurs abgeschlossen haben.
- Mit grundlegenden Nextflow-Konzepten und -Mechanismen vertraut sein (Prozesse, Channels, Operatoren)

---

## 0. Erste Schritte

#### Öffne den Training-Codespace

Falls du es noch nicht getan hast, stelle sicher, dass du die Trainingsumgebung wie in der [Umgebungseinrichtung](../envsetup/index.md) beschrieben öffnest.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Wechsle in das Projektverzeichnis

Lass uns in das Verzeichnis wechseln, in dem sich die Dateien für dieses Tutorial befinden.

```bash
cd side-quests/metadata
```

Du kannst VSCode so einstellen, dass es sich auf dieses Verzeichnis fokussiert:

```bash
code .
```

#### Überprüfe die Materialien

Du findest eine Haupt-Workflow-Datei und ein `data`-Verzeichnis, das ein Datenblatt und eine Handvoll Datendateien enthält.

??? abstract "Verzeichnisinhalt"

    ```console
    .
    ├── data
    │   ├── bonjour.txt
    │   ├── ciao.txt
    │   ├── guten_tag.txt
    │   ├── hallo.txt
    │   ├── hello.txt
    │   ├── hola.txt
    │   ├── salut.txt
    │   └── datasheet.csv
    ├── main.nf
    └── nextflow.config
    ```

Der Workflow in der `main.nf`-Datei ist ein Stub, den du schrittweise zu einem vollständig funktionierenden Workflow erweitern wirst.

Das Datenblatt listet die Pfade zu den Datendateien und einige zugehörige Metadaten auf, organisiert in 3 Spalten:

- `id`: selbsterklärend, eine ID, die der Datei zugewiesen wurde
- `character`: ein Charaktername, den wir später verwenden werden, um verschiedene Kreaturen zu zeichnen
- `data`: Pfade zu `.txt`-Dateien, die Begrüßungen in verschiedenen Sprachen enthalten

```console title="datasheet.csv"
id,character,recording
sampleA,squirrel,/workspaces/training/side-quests/metadata/data/bonjour.txt
sampleB,tux,/workspaces/training/side-quests/metadata/data/guten_tag.txt
sampleC,sheep,/workspaces/training/side-quests/metadata/data/hallo.txt
sampleD,turkey,/workspaces/training/side-quests/metadata/data/hello.txt
sampleE,stegosaurus,/workspaces/training/side-quests/metadata/data/hola.txt
sampleF,moose,/workspaces/training/side-quests/metadata/data/salut.txt
sampleG,turtle,/workspaces/training/side-quests/metadata/data/ciao.txt
```

Jede Datendatei enthält Begrüßungstext in einer von fünf Sprachen (fr: Französisch, de: Deutsch, es: Spanisch, it: Italienisch, en: Englisch).

Wir stellen dir auch ein containerisiertes Sprachanalyse-Tool namens `langid` zur Verfügung.

#### Überprüfe die Aufgabe

Deine Herausforderung besteht darin, einen Nextflow-Workflow zu schreiben, der:

1. Die Sprache in jeder Datei automatisch **identifiziert**
2. Dateien nach Sprachfamilie **gruppiert** (germanische vs. romanische Sprachen)
3. Die Verarbeitung für jede Datei basierend auf ihrer Sprache und Metadaten **anpasst**
4. Ausgaben nach Sprachgruppe **organisiert**

Dies stellt ein typisches Workflow-Muster dar, bei dem dateispezifische Metadaten Verarbeitungsentscheidungen steuern; genau die Art von Problem, die Metadaten-Maps elegant lösen.

#### Bereitschafts-Checkliste

Denkst du, du bist bereit einzutauchen?

- [ ] Ich verstehe das Ziel dieses Kurses und seine Voraussetzungen
- [ ] Mein Codespace ist aktiv
- [ ] Ich habe mein Arbeitsverzeichnis entsprechend eingestellt
- [ ] Ich verstehe die Aufgabe

Wenn du alle Kästchen abhaken kannst, kannst du loslegen.

---

## 1. Metadaten aus einem Datenblatt laden

Öffne die `main.nf`-Workflow-Datei, um den Workflow-Stub zu untersuchen, den wir dir als Ausgangspunkt geben.

```groovy title="main.nf" linenums="1"
#!/usr/bin/env nextflow

workflow  {

    ch_datasheet = channel.fromPath("./data/datasheet.csv")

}
```

Du siehst, dass wir eine grundlegende Channel-Factory eingerichtet haben, um das Beispiel-Datenblatt als Datei zu laden, aber das liest noch nicht den Inhalt der Datei ein.
Lass uns damit beginnen, das hinzuzufügen.

### 1.1. Inhalt mit `splitCsv` einlesen

Wir müssen einen Operator wählen, der den Dateiinhalt mit minimalem Aufwand unsererseits passend parst.
Da unser Datenblatt im CSV-Format vorliegt, ist dies ein Job für den [`splitCsv`](https://www.nextflow.io/docs/latest/reference/operator.html#splitcsv)-Operator, der jede Zeile in der Datei als Element im Channel lädt.

Nimm die folgenden Änderungen vor, um eine `splitCsv()`-Operation zum Channel-Konstruktionscode hinzuzufügen, plus eine `view()`-Operation, um zu überprüfen, dass der Inhalt der Datei korrekt in den Channel geladen wird.

=== "Nachher"

    ```groovy title="main.nf" linenums="3" hl_lines="4-5"
    workflow  {

        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .view()

    }
    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="3"
    workflow  {

        ch_datasheet = channel.fromPath("./data/datasheet.csv")

    }
    ```

Beachte, dass wir die Option `header: true` verwenden, um Nextflow mitzuteilen, dass die erste Zeile der CSV-Datei die Header-Zeile ist.

Lass uns schauen, was dabei herauskommt, ja?
Führe den Workflow aus:

```bash
nextflow run main.nf
```

??? success "Befehlsausgabe"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [exotic_albattani] DSL2 - revision: c0d03cec83

    [id:sampleA, character:squirrel, recording:/workspaces/training/side-quests/metadata/data/bonjour.txt]
    [id:sampleB, character:tux, recording:/workspaces/training/side-quests/metadata/data/guten_tag.txt]
    [id:sampleC, character:sheep, recording:/workspaces/training/side-quests/metadata/data/hallo.txt]
    [id:sampleD, character:turkey, recording:/workspaces/training/side-quests/metadata/data/hello.txt]
    [id:sampleE, character:stegosaurus, recording:/workspaces/training/side-quests/metadata/data/hola.txt]
    [id:sampleF, character:moose, recording:/workspaces/training/side-quests/metadata/data/salut.txt]
    [id:sampleG, character:turtle, recording:/workspaces/training/side-quests/metadata/data/ciao.txt]
    ```

Wir können sehen, dass der Operator eine Map von Schlüssel-Wert-Paaren für jede Zeile in der CSV-Datei erstellt hat, mit den Spaltenüberschriften als Schlüssel für die entsprechenden Werte.

Jeder Map-Eintrag entspricht einer Spalte in unserem Datenblatt:

- `id`
- `character`
- `recording`

Das ist großartig! Es macht es einfach, auf spezifische Felder aus jeder Datei zuzugreifen.
Zum Beispiel könnten wir auf die Datei-ID mit `id` oder auf den txt-Dateipfad mit `recording` zugreifen.

??? info "(Optional) Mehr über Maps"

    In Groovy, der Programmiersprache, auf der Nextflow aufbaut, ist eine Map eine Schlüssel-Wert-Datenstruktur ähnlich wie Dictionaries in Python, Objekte in JavaScript oder Hashes in Ruby.

    Hier ist ein ausführbares Skript, das zeigt, wie du eine Map definieren und auf ihren Inhalt zugreifen kannst:

    ```groovy title="examples/map_demo.nf"
    #!/usr/bin/env nextflow

    // Erstelle eine einfache Map
    def my_map = [id:'sampleA', character:'squirrel']

    // Gib die gesamte Map aus
    println "map: ${my_map}"

    // Greife auf einzelne Werte mit Punkt-Notation zu
    println "id: ${my_map.id}"
    println "character: ${my_map.character}"
    ```

    Obwohl es keinen richtigen `workflow`-Block hat, kann Nextflow dies ausführen, als wäre es ein Workflow:

    ```bash
    nextflow run examples/map_demo.nf
    ```

    Und hier ist, was du in der Ausgabe erwarten kannst:

    ```console title="Ausgabe"
     N E X T F L O W   ~  version 25.10.2

    Launching `map_demo.nf` [cheesy_plateau] DSL2 - revision: fae5b8496e

    map: [id:sampleA, character:squirrel]
    id: sampleA
    character: squirrel
    ```

### 1.2. Spezifische Felder mit `map` auswählen

Nehmen wir an, wir möchten auf die `character`-Spalte aus dem Datenblatt zugreifen und sie ausgeben.
Wir können den Nextflow-`map`-Operator verwenden, um über jedes Element in unserem Channel zu iterieren und speziell den `character`-Eintrag aus dem Map-Objekt auszuwählen.

Nimm die folgenden Änderungen am Workflow vor:

=== "Nachher"

    ```groovy title="main.nf" linenums="3" hl_lines="5-7"
    workflow  {

        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map{ row ->
                row.character
            }
            .view()

    }
    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="3"
    workflow  {

        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .view()

    }
    ```

Führe jetzt den Workflow erneut aus:

```bash
nextflow run main.nf
```

??? success "Befehlsausgabe"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [exotic_albattani] DSL2 - revision: c0d03cec83

    squirrel
    tux
    sheep
    turkey
    stegosaurus
    moose
    turtle
    ```

Erfolg! Wir haben die Map-Struktur aus unserem Datenblatt genutzt, um auf die Werte aus einzelnen Spalten für jede Zeile zuzugreifen.

Jetzt, da wir das Datenblatt erfolgreich eingelesen haben und Zugriff auf die Daten in jeder Zeile haben, können wir beginnen, unsere Pipeline-Logik zu implementieren.

### 1.3. Die Metadaten in einer 'Meta-Map' organisieren

Im aktuellen Zustand des Workflows stehen die Eingabedateien (unter dem `recording`-Schlüssel) und die zugehörigen Metadaten (`id`, `character`) alle auf derselben Ebene, als wären sie alle in einem großen Beutel.
Die praktische Konsequenz ist, dass jeder Prozess, der diesen Channel konsumiert, mit dieser Struktur im Hinterkopf konfiguriert werden müsste:

```groovy
    input:
    tuple val(id), val(character), file(recording)
```

Das ist in Ordnung, solange sich die Anzahl der Spalten im Datenblatt nicht ändert.
Wenn du jedoch auch nur eine Spalte zum Datenblatt hinzufügst, wird die Form des Channels nicht mehr mit dem übereinstimmen, was der Prozess erwartet, und der Workflow wird Fehler produzieren.
Es macht den Prozess auch schwer mit anderen zu teilen, die möglicherweise leicht unterschiedliche Eingabedaten haben, und du musst möglicherweise Variablen in den Prozess hart codieren, die vom Script-Block nicht benötigt werden.

Um dieses Problem zu vermeiden, müssen wir einen Weg finden, die Channel-Struktur konsistent zu halten, unabhängig davon, wie viele Spalten das Datenblatt enthält.

Wir können das tun, indem wir alle Metadaten in einem Element innerhalb des Tupels sammeln, das wir die Metadaten-Map oder einfacher 'Meta-Map' nennen werden.

Nimm die folgenden Änderungen an der `map`-Operation vor:

=== "Nachher"

    ```groovy title="main.nf" linenums="5" hl_lines="4"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                    [ [id: row.id, character: row.character], row.recording ]
            }
            .view()
    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="5" hl_lines="4"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map{ row ->
                row.character
            }
            .view()
    ```

Wir haben unsere Channel-Elemente in ein Tupel umstrukturiert, das aus zwei Elementen besteht: der Meta-Map und dem entsprechenden Dateiobjekt.

Lass uns den Workflow ausführen:

```bash
nextflow run main.nf
```

??? success "Befehlsausgabe"

    ```console title="Meta-Map ansehen"
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [lethal_booth] DSL2 - revision: 0d8f844c07

    [[id:sampleA, character:squirrel], /workspaces/training/side-quests/metadata/data/bonjour.txt]
    [[id:sampleB, character:tux], /workspaces/training/side-quests/metadata/data/guten_tag.txt]
    [[id:sampleC, character:sheep], /workspaces/training/side-quests/metadata/data/hallo.txt]
    [[id:sampleD, character:turkey], /workspaces/training/side-quests/metadata/data/hello.txt]
    [[id:sampleE, character:stegosaurus], /workspaces/training/side-quests/metadata/data/hola.txt]
    [[id:sampleF, character:moose], /workspaces/training/side-quests/metadata/data/salut.txt]
    [[id:sampleG, character:turtle], /workspaces/training/side-quests/metadata/data/ciao.txt]
    ```

Jetzt enthält jedes Element im Channel zuerst die Metadaten-Map und zweitens das entsprechende Dateiobjekt:

```console title="Beispiel-Ausgabestruktur"
[
  [id:sampleA, character:squirrel],
  /workspaces/training/side-quests/metadata/data/bonjour.txt
]
```

Als Ergebnis macht das Hinzufügen weiterer Spalten im Datenblatt mehr Metadaten in der `meta`-Map verfügbar, ändert aber nicht die Channel-Form.
Dies ermöglicht es uns, Prozesse zu schreiben, die den Channel konsumieren, ohne die Metadaten-Elemente in die Eingabespezifikation hart codieren zu müssen:

```groovy title="Syntax-Beispiel"
    input:
    tuple val(meta), file(recording)
```

Dies ist eine weit verbreitete Konvention zur Organisation von Metadaten in Nextflow-Workflows.

### Fazit

In diesem Abschnitt hast du gelernt:

- **Warum Metadaten wichtig sind:** Das Behalten von Metadaten bei deinen Daten bewahrt wichtige Dateiinformationen durch den gesamten Workflow hindurch.
- **Wie man Datenblätter einliest:** Verwendung von `splitCsv` zum Lesen von CSV-Dateien mit Header-Informationen und Umwandlung von Zeilen in strukturierte Daten
- **Wie man eine Meta-Map erstellt:** Trennung von Metadaten von Dateidaten unter Verwendung der Tupel-Struktur `[ [id:wert, ...], datei ]`

---

## 2. Metadaten manipulieren

Jetzt, da wir unsere Metadaten geladen haben, lass uns etwas damit machen!

Wir werden ein Tool namens [`langid`](https://github.com/saffsd/langid.py) verwenden, um die Sprache zu identifizieren, die in der Aufnahmedatei jeder Kreatur enthalten ist.
Das Tool ist auf eine Reihe von Sprachen vortrainiert, und bei einem Textausschnitt gibt es eine Sprachvorhersage und einen zugehörigen Wahrscheinlichkeitswert aus, beide nach `stdout`.

### 2.1. Importiere den Prozess und untersuche den Code

Wir stellen dir ein vorgefertigtes Prozessmodul namens `IDENTIFY_LANGUAGE` zur Verfügung, das das `langid`-Tool umschließt, sodass du nur eine Include-Anweisung vor dem Workflow-Block hinzufügen musst.

Nimm die folgende Änderung am Workflow vor:

=== "Nachher"

    ```groovy title="main.nf" linenums="1" hl_lines="3"
    #!/usr/bin/env nextflow

    include { IDENTIFY_LANGUAGE } from './modules/langid.nf'

    workflow {
    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="1"
    #!/usr/bin/env nextflow

    workflow {
    ```

Du kannst die Moduldatei öffnen, um ihren Code zu untersuchen:

```groovy title="modules/langid.nf" linenums="1" hl_lines="9 12"
#!/usr/bin/env nextflow

// Verwende langid zur Vorhersage der Sprache jeder Eingabedatei
process IDENTIFY_LANGUAGE {

    container 'community.wave.seqera.io/library/pip_langid:b2269f456a5629ff'

    input:
    tuple val(meta), path(file)

    output:
    tuple val(meta), path(file), stdout

    script:
    """
    langid < ${file} -l en,de,fr,es,it | sed -E "s/.*\\('([a-z]+)'.*/\\1/" | tr -d '\\n'
    """
}
```

Wie du sehen kannst, verwendet die Eingabedefinition dieselbe `tuple val(meta), path(file)`-Struktur, die wir gerade auf unseren Eingabe-Channel angewendet haben.

Die Ausgabedefinition ist als Tupel mit einer ähnlichen Struktur wie die Eingabe strukturiert, außer dass sie auch `stdout` als drittes Element enthält.
Dieses `tuple val(meta), path(file), <output>`-Muster hält die Metadaten mit sowohl den Eingabedaten als auch den Ausgaben verbunden, während sie durch die Pipeline fließen.

Beachte, dass wir hier Nextflows [`stdout`](https://www.nextflow.io/docs/latest/process.html#outputs)-Ausgabequalifizierer verwenden, weil das Tool seine Ausgabe direkt auf die Konsole ausgibt, anstatt eine Datei zu schreiben; und wir verwenden `sed` in der Befehlszeile, um den Wahrscheinlichkeitswert zu entfernen, die Zeichenkette durch Entfernen von Zeilenumbruchzeichen zu bereinigen und nur die Sprachvorhersage zurückzugeben.

### 2.2. Füge einen Aufruf zu `IDENTIFY_LANGUAGE` hinzu

Jetzt, da der Prozess dem Workflow zur Verfügung steht, können wir einen Aufruf zum `IDENTIFY_LANGUAGE`-Prozess hinzufügen, um ihn auf dem Daten-Channel auszuführen.

Nimm die folgenden Änderungen am Workflow vor:

=== "Nachher"

    ```groovy title="main.nf" linenums="7" hl_lines="7-9"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [[id: row.id, character: row.character], row.recording]
            }

        // Führe langid aus, um die Sprache jeder Begrüßung zu identifizieren
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out.view()
    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="7" hl_lines="6"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [[id: row.id, character: row.character], row.recording]
            }
            .view()
    ```

Beachte, dass wir die ursprüngliche `.view()`-Operation in der Channel-Konstruktion entfernt haben.

Wir können jetzt den Workflow ausführen:

```bash
nextflow run main.nf
```

??? success "Befehlsausgabe"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [voluminous_mcnulty] DSL2 - revision: f9bcfebabb

    executor >  local (7)
    [4e/f722fe] IDENTIFY_LANGUAGE (7) [100%] 7 of 7 ✔
    [[id:sampleA, character:squirrel], /workspaces/training/side-quests/metadata/work/eb/f7148ebdd898fbe1136bec6a714acb/bonjour.txt, fr]
    [[id:sampleB, character:tux], /workspaces/training/side-quests/metadata/work/16/71d72410952c22cd0086d9bca03680/guten_tag.txt, de]
    [[id:sampleD, character:turkey], /workspaces/training/side-quests/metadata/work/c4/b7562adddc1cc0b7d414ec45d436eb/hello.txt, en]
    [[id:sampleC, character:sheep], /workspaces/training/side-quests/metadata/work/ea/04f5d979429e4455e14b9242fb3b45/hallo.txt, de]
    [[id:sampleF, character:moose], /workspaces/training/side-quests/metadata/work/5a/6c2b84bf8fadb98e28e216426be079/salut.txt, fr]
    [[id:sampleE, character:stegosaurus], /workspaces/training/side-quests/metadata/work/af/ee7c69bcab891c40d0529305f6b9e7/hola.txt, es]
    [[id:sampleG, character:turtle], /workspaces/training/side-quests/metadata/work/4e/f722fe47271ba7ebcd69afa42964ca/ciao.txt, it]
    ```

Ausgezeichnet! Wir haben jetzt eine Vorhersage darüber, welche Sprache jeder Charakter spricht.

Und wie bereits erwähnt, haben wir auch die Eingabedatei und die Meta-Map in die Ausgabe aufgenommen, was bedeutet, dass beide mit den neuen Informationen verbunden bleiben, die wir gerade produziert haben.
Dies wird sich im nächsten Schritt als nützlich erweisen.

!!! note

    Allgemeiner macht es dieses Muster, die Meta-Map mit Ergebnissen verbunden zu halten, einfacher, verwandte Ergebnisse zu verknüpfen, die dieselben Identifikatoren teilen.

    Wie du bereits gelernt hast, kannst du dich nicht auf die Reihenfolge der Elemente in Channels verlassen, um Ergebnisse über sie hinweg zuzuordnen.
    Stattdessen musst du Schlüssel verwenden, um Daten korrekt zuzuordnen, und Meta-Maps bieten eine ideale Struktur für diesen Zweck.

    Wir untersuchen diesen Anwendungsfall detailliert in der Side Quest [Aufteilen & Gruppieren](./splitting_and_grouping.md).

### 2.3. Erweitere Metadaten mit Prozessausgaben

Da die Ergebnisse, die wir gerade produziert haben, selbst eine Form von Metadaten über den Inhalt der Dateien sind, wäre es nützlich, sie zu unserer Meta-Map hinzuzufügen.

Allerdings wollen wir die bestehende Meta-Map nicht direkt modifizieren.
Aus technischer Sicht ist es _möglich_, das zu tun, aber es ist unsicher.

Stattdessen erstellen wir eine neue Meta-Map, die den Inhalt der bestehenden Meta-Map plus ein neues `lang: lang_id`-Schlüssel-Wert-Paar enthält, das die neuen Informationen hält, unter Verwendung des `+`-Operators (eine Groovy-Funktion).
Und wir werden dies mit einer [`map`](https://www.nextflow.io/docs/latest/operator.html#map)-Operation kombinieren, um die alte Map durch die neue zu ersetzen.

Hier sind die Änderungen, die du am Workflow vornehmen musst:

=== "Nachher"

    ```groovy title="main.nf" linenums="13" hl_lines="3-7"
        // Führe langid aus, um die Sprache jeder Begrüßung zu identifizieren
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out
            .map { meta, file, lang_id ->
                [meta + [lang: lang_id], file]
            }
            .view()
    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="13" hl_lines="3"
        // Führe langid aus, um die Sprache jeder Begrüßung zu identifizieren
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out.view()
    ```

Wenn du mit dem `+`-Operator noch nicht vertraut bist oder wenn dies verwirrend erscheint, nimm dir ein paar Minuten Zeit, um die detaillierte Erklärung unten durchzugehen.

??? info "Erstellung der neuen Meta-Map unter Verwendung des `+`-Operators"

    **Zuerst musst du wissen, dass wir den Inhalt von zwei Maps mit dem Groovy-Operator `+` zusammenführen können.**

    Nehmen wir an, wir haben die folgenden Maps:

    ```groovy
    map1 = [id: 'sampleA', character: 'squirrel']
    map2 = [lang: 'fr']
    ```

    Wir können sie so zusammenführen:

    ```groovy
    new_map = map1 + map2
    ```

    Der Inhalt von `new_map` wird sein:

    ```groovy
    [id: 'sampleA', character: 'squirrel', lang: 'fr']
    ```

    Großartig!

    **Aber was ist, wenn du ein Feld hinzufügen musst, das noch nicht Teil einer Map ist?**

    Nehmen wir an, du beginnst wieder von `map1`, aber die Sprachvorhersage ist nicht in ihrer eigenen Map (es gibt kein `map2`).
    Stattdessen wird sie in einer Variablen namens `lang_id` gehalten, und du weißt, dass du ihren Wert (`'fr'`) mit dem Schlüssel `lang` speichern möchtest.

    Du kannst tatsächlich folgendes tun:

    ```groovy
    new_map = [map1 + [lang: lang_id]]
    ```

    Hier erstellt `[lang: new_info]` spontan eine neue unbenannte Map, und `map1 + ` führt `map1` mit der neuen unbenannten Map zusammen, wodurch der gleiche `new_map`-Inhalt wie zuvor entsteht.

    Praktisch, oder?

    **Lass uns das jetzt in den Kontext einer Nextflow `channel.map()`-Operation übertragen.**

    Der Code wird zu:

    ```groovy
    .map { map1, lang_id ->
        [map1 + [lang: lang_id]]
    }
    ```

    Dies macht folgendes:

    - `map1, lang_id ->` nimmt die zwei Elemente im Tupel
    - `[map1 + [lang: lang_id]]` erstellt die neue Map wie oben detailliert

    Die Ausgabe ist eine einzelne unbenannte Map mit demselben Inhalt wie `new_map` in unserem obigen Beispiel.
    Wir haben also effektiv transformiert:

    ```groovy
    [id: 'sampleA', character: 'squirrel'], 'it'
    ```

    in:

    ```groovy
    [id: 'sampleA', character: 'squirrel', lang: 'fr']
    ```

    Hoffentlich kannst du sehen, dass, wenn wir `map1` in `meta` ändern, das im Grunde alles ist, was wir brauchen, um die Sprachvorhersage zu unserer Meta-Map in unserem Workflow hinzuzufügen.

    Außer einer Sache!

    Im Fall unseres Workflows **müssen wir auch die Anwesenheit des `file`-Objekts im Tupel berücksichtigen**, das aus `meta, file, lang_id` besteht.

    Also würde der Code hier zu:

    ```groovy
    .map { meta, file, lang_id ->
        [meta + [lang: lang_id], file]
    }
    ```

    Wenn du Schwierigkeiten hast zu verstehen, warum die `file` sich in der `map`-Operation zu bewegen scheint, stelle dir vor, dass statt `[meta + [lang: lang_id], file]` diese Zeile `[new_map, file]` lautet.
    Dies sollte deutlicher machen, dass wir einfach die `file` an ihrer ursprünglichen Position an zweiter Stelle im Tupel lassen. Wir haben nur den `new_info`-Wert genommen und ihn in die Map gefaltet, die an erster Stelle steht.

    **Und das bringt uns zurück zur `tuple val(meta), path(file)`-Channel-Struktur!**

Sobald du sicher bist, dass du verstehst, was dieser Code macht, führe den Workflow aus, um zu sehen, ob es funktioniert hat:

```bash
nextflow run main.nf -resume
```

??? success "Befehlsausgabe"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [cheeky_fermat] DSL2 - revision: d096281ee4

    [4e/f722fe] IDENTIFY_LANGUAGE (7) [100%] 7 of 7, cached: 7 ✔
    [[id:sampleA, character:squirrel, lang:fr], /workspaces/training/side-quests/metadata/work/eb/f7148ebdd898fbe1136bec6a714acb/bonjour.txt]
    [[id:sampleB, character:tux, lang:de], /workspaces/training/side-quests/metadata/work/16/71d72410952c22cd0086d9bca03680/guten_tag.txt]
    [[id:sampleC, character:sheep, lang:de], /workspaces/training/side-quests/metadata/work/ea/04f5d979429e4455e14b9242fb3b45/hallo.txt]
    [[id:sampleD, character:turkey, lang:en], /workspaces/training/side-quests/metadata/work/c4/b7562adddc1cc0b7d414ec45d436eb/hello.txt]
    [[id:sampleF, character:moose, lang:fr], /workspaces/training/side-quests/metadata/work/5a/6c2b84bf8fadb98e28e216426be079/salut.txt]
    [[id:sampleE, character:stegosaurus, lang:es], /workspaces/training/side-quests/metadata/work/af/ee7c69bcab891c40d0529305f6b9e7/hola.txt]
    [[id:sampleG, character:turtle, lang:it], /workspaces/training/side-quests/metadata/work/4e/f722fe47271ba7ebcd69afa42964ca/ciao.txt]
    ```

Ja, das stimmt!
Wir haben die Ausgabe des Prozesses von `meta, file, lang_id` sauber neu organisiert, sodass `lang_id` jetzt einer der Schlüssel in der Meta-Map ist, und die Tupel des Channels wieder dem `meta, file`-Modell entsprechen.

<!-- TODO (future) Should we also show how to remove a key using subMap?! Or note where to find that. -->

### 2.4. Weise eine Sprachgruppe mit bedingten Anweisungen zu

Jetzt, da wir unsere Sprachvorhersagen haben, lass uns die Informationen verwenden, um einige neue Gruppierungen zuzuweisen.

In unseren Beispieldaten können die von unseren Charakteren verwendeten Sprachen in germanische Sprachen (Englisch, Deutsch) und romanische Sprachen (Französisch, Spanisch, Italienisch) gruppiert werden.
Es könnte nützlich sein, diese Klassifizierung später in der Pipeline verfügbar zu haben, also lass uns diese Information zur Meta-Map hinzufügen.

Und gute Nachrichten, dies ist noch ein Fall, der sich perfekt zur Verwendung des `map`-Operators eignet!

Konkret werden wir eine Variable namens `lang_group` definieren, einfache bedingte Logik verwenden, um zu bestimmen, welchen Wert wir der `lang_group` für jedes Datenstück zuweisen sollen.

Die allgemeine Syntax wird so aussehen:

```groovy
.map { meta, file ->

    // bedingte Logik zur Definition von lang_group kommt hier hin

    [meta + [lang_group: lang_group], file]
}
```

Du kannst sehen, dass dies der spontanen Map-Zusammenführungsoperation sehr ähnlich ist, die wir im vorherigen Schritt verwendet haben.
Wir müssen nur die bedingten Anweisungen schreiben.

Hier ist die bedingte Logik, die wir anwenden wollen:

- Definiere eine Variable namens `lang_group` mit dem Standardwert `'unknown'`.
- Wenn `lang` entweder Deutsch (`'de'`) oder Englisch (`'en'`) ist, ändere `lang_group` zu `germanic`.
- Ansonsten, wenn `lang` in einer Liste enthalten ist, die Französisch (`'fr'`), Spanisch (`'es'`) und Italienisch (`'it'`) enthält, ändere `lang_group` zu `romance`.

Versuche es selbst zu schreiben, wenn du bereits weißt, wie man bedingte Anweisungen in Nextflow schreibt.

!!! tip

    Du kannst auf den Wert von `lang` innerhalb der Map-Operation mit `meta.lang` zugreifen.

Du solltest am Ende die folgenden Änderungen am Workflow vornehmen:

=== "Nachher"

    ```groovy title="main.nf" linenums="13" hl_lines="7-19"
        // Führe langid aus, um die Sprache jeder Begrüßung zu identifizieren
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out
            .map { meta, file, lang_id ->
                [meta + [lang: lang_id], file]
            }
            .map { meta, file ->

                def lang_group = "unknown"
                if (meta.lang.equals("de") || meta.lang.equals('en')) {
                    lang_group = "germanic"
                }
                else if (meta.lang in ["fr", "es", "it"]) {
                    lang_group = "romance"
                }

                [meta + [lang_group: lang_group], file]
            }
            .view()
    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="13" hl_lines="7"
        // Führe langid aus, um die Sprache jeder Begrüßung zu identifizieren
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out
            .map { meta, file, lang_id ->
                [meta + [lang: lang_id], file]
            }
            .view()
    ```

Hier sind die wichtigsten Punkte:

- Wir verwenden `def lang_group = "unknown"`, um die `lang_group`-Variable mit Standardwert auf `unknown` gesetzt zu erstellen.
- Wir verwenden eine `if {} else if {}`-Struktur für die bedingte Logik, mit alternativen `.equals()`-Tests für die zwei germanischen Sprachen und einem Test auf Existenz in einer Liste für die drei romanischen Sprachen.
- Wir verwenden die `meta + [lang_group:lang_group]`-Zusammenführungsoperation wie zuvor, um die aktualisierte Meta-Map zu generieren.

<!-- TODO (future) Add note/links to relevant docs in additional resources section -->

Sobald das alles Sinn macht, führe den Workflow erneut aus, um das Ergebnis zu sehen:

```bash
nextflow run main.nf -resume
```

??? success "Befehlsausgabe"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [wise_almeida] DSL2 - revision: 46778c3cd0

    [da/652cc6] IDENTIFY_LANGUAGE (7) [100%] 7 of 7, cached: 7 ✔
    [[id:sampleA, character:squirrel, lang:fr, lang_group:romance], /workspaces/training/side-quests/metadata/data/bonjour.txt]
    [[id:sampleB, character:tux, lang:de, lang_group:germanic], /workspaces/training/side-quests/metadata/data/guten_tag.txt]
    [[id:sampleC, character:sheep, lang:de, lang_group:germanic], /workspaces/training/side-quests/metadata/data/hallo.txt]
    [[id:sampleD, character:turkey, lang:en, lang_group:germanic], /workspaces/training/side-quests/metadata/data/hello.txt]
    [[id:sampleE, character:stegosaurus, lang:es, lang_group:romance], /workspaces/training/side-quests/metadata/data/hola.txt]
    [[id:sampleF, character:moose, lang:fr, lang_group:romance], /workspaces/training/side-quests/metadata/data/salut.txt]
    [[id:sampleG, character:turtle, lang:it, lang_group:romance], /workspaces/training/side-quests/metadata/data/ciao.txt]
    ```

Wie du sehen kannst, behalten die Channel-Elemente ihre `[meta, file]`-Struktur bei, aber die Meta-Map enthält jetzt diese neue Klassifizierung.

### Fazit

In diesem Abschnitt hast du gelernt, wie man:

- **Eingabe-Metadaten auf Ausgabe-Channels anwendet**: Das Kopieren von Metadaten auf diese Weise ermöglicht es uns, Ergebnisse später basierend auf Metadateninhalten zu verknüpfen.
- **Eigene Schlüssel erstellt**: Du hast zwei neue Schlüssel in deiner Meta-Map erstellt, die du mit `meta + [new_key:value]` in die bestehende Meta-Map zusammengeführt hast. Einen basierend auf einem berechneten Wert aus einem Prozess und einen basierend auf einer Bedingung, die du im `map`-Operator gesetzt hast.

Diese ermöglichen es dir, neue und bestehende Metadaten mit Dateien zu verknüpfen, während du durch deine Pipeline fortschreitest.
Selbst wenn du Metadaten nicht als Teil eines Prozesses verwendest, macht es das Beibehalten der Meta-Map mit den Daten wie hier einfach, alle relevanten Informationen zusammenzuhalten.

---

## 3. Meta-Map-Informationen in einem Prozess verwenden

Jetzt, da du weißt, wie man die Meta-Map erstellt und aktualisiert, können wir zum wirklich spaßigen Teil kommen: die Metadaten tatsächlich in einem Prozess verwenden.

Genauer gesagt, werden wir unserem Workflow einen zweiten Schritt hinzufügen, um jedes Tier als ASCII-Art zu zeichnen und es den aufgezeichneten Text in einer Sprechblase sagen zu lassen.
Wir werden dies mit einem Tool namens [`cowpy`](https://github.com/jeffbuttars/cowpy) tun.

??? info "Was macht `cowpy`?"

    `cowpy` ist ein Befehlszeilen-Tool, das ASCII-Art generiert, um beliebige Texteingaben auf lustige Weise anzuzeigen.
    Es ist eine Python-Implementierung des klassischen [cowsay](https://en.wikipedia.org/wiki/Cowsay)-Tools von Tony Monroe.

    ```console
    cowpy "Hello Nextflow"
    ```

    ```console
    ______________________________________________________
    < Hello Nextflow >
    ------------------------------------------------------
        \   ^__^
          \  (oo)\_______
            (__)\       )\/\
              ||----w |
              ||     ||
    ```

    Optional kannst du einen Charakter (oder 'Cowacter') auswählen, der anstelle der Standard-Kuh verwendet werden soll.

    ```console
    cowpy "Hello Nextflow" -c tux
    ```

    ```console
    __________________
    < Hello Nextflow >
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

Wenn du den Hello Nextflow-Kurs durchgearbeitet hast, hast du dieses Tool bereits in Aktion gesehen.
Wenn nicht, keine Sorge; wir werden alles abdecken, was du wissen musst, während wir weitermachen.

### 3.1. Importiere den Prozess und untersuche den Code

Wir stellen dir ein vorgefertigtes Prozessmodul namens `COWPY` zur Verfügung, das das `cowpy`-Tool umschließt, sodass du nur eine Include-Anweisung vor dem Workflow-Block hinzufügen musst.

Nimm die folgende Änderung am Workflow vor:

=== "Nachher"

    ```groovy title="main.nf" linenums="1" hl_lines="4"
    #!/usr/bin/env nextflow

    include { IDENTIFY_LANGUAGE } from './modules/langid.nf'
    include { COWPY } from './modules/cowpy.nf'

    workflow {
    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="1"
    #!/usr/bin/env nextflow

    include { IDENTIFY_LANGUAGE } from './modules/langid.nf'

    workflow {
    ```

Du kannst die Moduldatei öffnen, um ihren Code zu untersuchen:

```groovy title="modules/cowpy.nf" linenums="1"
#!/usr/bin/env nextflow

// Generiere ASCII-Art mit cowpy
process COWPY {

    publishDir "results/", mode: 'copy'

    container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

    input:
    path input_file
    val character

    output:
    path "cowpy-${input_file}"

    script:
    """
    cat ${input_file} | cowpy -c ${character} > cowpy-${input_file}
    """
}
```

Wie du sehen kannst, ist dieser Prozess derzeit so konzipiert, dass er eine Eingabedatei (die den anzuzeigenden Text enthält) und einen Wert entgegennimmt, der den Charakter angibt, der in ASCII-Art gezeichnet werden soll, normalerweise auf Workflow-Ebene durch einen Befehlszeilenparameter bereitgestellt.

### 3.2. Übergebe ein Meta-Map-Feld als Eingabe

Als wir das `cowpy`-Tool im Hello Nextflow-Kurs verwendet haben, haben wir einen Befehlszeilenparameter verwendet, um zu bestimmen, welcher Charakter zum Zeichnen des finalen Bildes verwendet werden soll.
Das machte Sinn, weil wir nur ein Bild pro Ausführung der Pipeline generierten.

In diesem Tutorial möchten wir jedoch ein geeignetes Bild für jedes Subjekt generieren, das wir verarbeiten, sodass die Verwendung eines Befehlszeilenparameters zu einschränkend wäre.

Gute Nachrichten: Wir haben eine `character`-Spalte in unserem Datenblatt und daher in unserer Meta-Map.
Lass uns diese verwenden, um den Charakter festzulegen, den der Prozess für jeden Eintrag verwenden soll.

Zu diesem Zweck müssen wir drei Dinge tun:

1. Dem Ausgabe-Channel, der aus dem vorherigen Prozess kommt, einen Namen geben, damit wir bequemer damit arbeiten können.
2. Bestimmen, wie man auf die interessierenden Informationen zugreift
3. Einen Aufruf zum zweiten Prozess hinzufügen und die Informationen entsprechend einspeisen.

Lass uns anfangen.

#### 3.2.1. Benenne den vorherigen Ausgabe-Channel

Wir haben die vorherigen Manipulationen direkt auf dem Ausgabe-Channel des ersten Prozesses angewendet, `IDENTIFY_LANGUAGE.out`.
Um den Inhalt dieses Channels an den nächsten Prozess weiterzugeben (und dies auf eine klare und leicht lesbare Weise zu tun), möchten wir ihm einen eigenen Namen geben, `ch_languages`.

Wir können das mit dem [`set`](https://www.nextflow.io/docs/latest/reference/operator.html#set)-Operator tun.

Ersetze im Haupt-Workflow den `.view()`-Operator durch `.set { ch_languages }` und füge eine Zeile hinzu, die testet, dass wir uns auf den Channel namentlich beziehen können.

=== "Nachher"

    ```groovy title="main.nf" linenums="14" hl_lines="19 21 22"
        // Führe langid aus, um die Sprache jeder Begrüßung zu identifizieren
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out
            .map { meta, file, lang_id ->
                [meta + [lang: lang_id], file]
            }
            .map { meta, file ->

                def lang_group = "unknown"
                if (meta.lang.equals("de") || meta.lang.equals('en')) {
                    lang_group = "germanic"
                }
                else if (meta.lang in ["fr", "es", "it"]) {
                    lang_group = "romance"
                }

                [meta + [lang_group: lang_group], file]
            }
            .set { ch_languages }

        // Temporär: Blick in ch_languages
        ch_languages.view()
    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="14" hl_lines="19"
        // Führe langid aus, um die Sprache jeder Begrüßung zu identifizieren
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out
            .map { meta, file, lang_id ->
                [meta + [lang: lang_id], file]
            }
            .map { meta, file ->

                def lang_group = "unknown"
                if (meta.lang.equals("de") || meta.lang.equals('en')) {
                    lang_group = "germanic"
                }
                else if (meta.lang in ["fr", "es", "it"]) {
                    lang_group = "romance"
                }

                [meta + [lang_group: lang_group], file]
            }
            .view()
    ```

Lass uns das ausführen:

```bash
nextflow run main.nf
```

??? success "Befehlsausgabe"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `./main.nf` [friendly_austin] DSL2 - revision: 3dbe460fd6

    [36/cca6a7] IDENTIFY_LANGUAGE (7) | 7 of 7 ✔
    [[id:sampleB, character:tux, lang:de, lang_group:germanic], /workspaces/training/side-quests/metadata/work/e2/6db2402d83cf72081bcd2d11784714/guten_tag.txt]
    [[id:sampleA, character:squirrel, lang:fr, lang_group:romance], /workspaces/training/side-quests/metadata/work/6c/114c818317d169457d6e7336d5d55b/bonjour.txt]
    [[id:sampleC, character:sheep, lang:de, lang_group:germanic], /workspaces/training/side-quests/metadata/work/55/68c69c5efb527f3604ddb3daab8057/hallo.txt]
    [[id:sampleD, character:turkey, lang:en, lang_group:germanic], /workspaces/training/side-quests/metadata/work/2a/4752055ccb5d1370b0ef9da41d3993/hello.txt]
    [[id:sampleE, character:stegosaurus, lang:es, lang_group:romance], /workspaces/training/side-quests/metadata/work/f4/fcd3186dc666d5d239ffa6c37d125d/hola.txt]
    [[id:sampleF, character:moose, lang:fr, lang_group:romance], /workspaces/training/side-quests/metadata/work/c3/3b2627f733f278a7088332a5806108/salut.txt]
    [[id:sampleG, character:turtle, lang:it, lang_group:romance], /workspaces/training/side-quests/metadata/work/36/cca6a7dbfa26ac24f9329787a32e9d/ciao.txt]
    ```

Dies bestätigt, dass wir uns jetzt namentlich auf den Channel beziehen können.

#### 3.2.2. Greife auf die Datei- und Charakter-Metadaten zu

Wir wissen aus der Betrachtung des Modulcodes, dass der `COWPY`-Prozess erwartet, eine Textdatei und einen `character`-Wert zu erhalten.
Um den Aufruf zum `COWPY`-Prozess zu schreiben, müssen wir nur wissen, wie man das entsprechende Dateiobjekt und die Metadaten aus jedem Element im Channel extrahiert.

Wie oft ist der einfachste Weg, das zu tun, eine `map`-Operation zu verwenden.

Unser Channel enthält Tupel, die als `[meta, file]` strukturiert sind, sodass wir direkt auf das `file`-Objekt zugreifen können, und wir können auf den `character`-Wert, der innerhalb der Meta-Map gespeichert ist, zugreifen, indem wir ihn als `meta.character` referenzieren.

Nimm im Haupt-Workflow die folgenden Code-Änderungen vor:

=== "Nachher"

    ```groovy title="main.nf" linenums="34"
        // Temporär: Zugriff auf Datei und Charakter
        ch_languages.map { meta, file -> file }.view { file -> "File: " + file }
        ch_languages.map { meta, file -> meta.character }.view { character -> "Character: " + character }
    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="34"
        // Temporär: Blick in ch_languages
        ch_languages.view()
    ```

Beachte, dass wir Closures (wie `{ file -> "File: " + file }`) verwenden, um die Ausgabe der `.view`-Operationen lesbarer zu machen.

Lass uns das ausführen:

```bash
nextflow run main.nf -resume
```

??? success "Befehlsausgabe"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `./main.nf` [cheesy_cantor] DSL2 - revision: 15af9c1ec7

    [43/05df08] IDENTIFY_LANGUAGE (7) [100%] 7 of 7, cached: 7 ✔
    Character: squirrel
    File: /workspaces/training/side-quests/metadata/work/8d/4b9498bbccb7a74f04e41877cdc3e5/bonjour.txt
    File: /workspaces/training/side-quests/metadata/work/d3/604274985406e40d79021dea658e60/guten_tag.txt
    Character: tux
    Character: turkey
    File: /workspaces/training/side-quests/metadata/work/d4/fafcc9415b61d2b0fea872e6a05e8a/hello.txt
    File: /workspaces/training/side-quests/metadata/work/02/468ac9efb27f636715e8144b37e9a7/hallo.txt
    Character: sheep
    Character: moose
    Character: stegosaurus
    File: /workspaces/training/side-quests/metadata/work/d4/61a7e1188b4f2742bc72004e226eca/salut.txt
    File: /workspaces/training/side-quests/metadata/work/ae/68364be238c11149c588bf6fc858b1/hola.txt
    File: /workspaces/training/side-quests/metadata/work/43/05df081af5d879ab52e5828fa0357e/ciao.txt
    Character: turtle
    ```

_Die Dateipfade und Charakterwerte können in deiner Ausgabe in einer anderen Reihenfolge herauskommen._

Dies bestätigt, dass wir auf die Datei und den Charakter für jedes Element im Channel zugreifen können.

#### 3.2.3. Rufe den `COWPY`-Prozess auf

Jetzt lass uns alles zusammenfügen und tatsächlich den `COWPY`-Prozess auf dem `ch_languages`-Channel aufrufen.

Nimm im Haupt-Workflow die folgenden Code-Änderungen vor:

=== "Nachher"

    ```groovy title="main.nf" linenums="34"
        // Führe cowpy aus, um ASCII-Art zu generieren
        COWPY(
            ch_languages.map { meta, file -> file },
            ch_languages.map { meta, file -> meta.character }
        )
    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="34"
        // Temporär: Zugriff auf Datei und Charakter
        ch_languages.map { meta, file -> [file, meta.character] }
            .view()
    ```

Du siehst, wir kopieren einfach die zwei Map-Operationen (minus die `.view()`-Anweisungen) als Eingaben für den Prozessaufruf.
Stelle nur sicher, dass du das Komma zwischen ihnen nicht vergisst!

Es ist etwas umständlich, aber wir werden sehen, wie man das im nächsten Abschnitt verbessern kann.

Lass uns das ausführen:

```bash
nextflow run main.nf -resume
```

??? success "Befehlsausgabe"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [suspicious_crick] DSL2 - revision: 25541014c5

    executor >  local (7)
    [43/05df08] IDENTIFY_LANGUAGE (7) [100%] 7 of 7, cached: 7 ✔
    [e7/317c18] COWPY (6)             [100%] 7 of 7 ✔
    ```

Wenn du ins Ergebnisverzeichnis schaust, solltest du die einzelnen Dateien sehen, die die ASCII-Art jeder Begrüßung enthalten, gesprochen vom entsprechenden Charakter.

??? abstract "Verzeichnis- und Beispieldateiinhalte"

    ```console
    results/
    ├── cowpy-bonjour.txt
    ├── cowpy-ciao.txt
    ├── cowpy-guten_tag.txt
    ├── cowpy-hallo.txt
    ├── cowpy-hello.txt
    ├── cowpy-hola.txt
    └── cowpy-salut.txt
    ```

    ```text title="results/cowpy-bonjour.txt"
     _________________
    / Bonjour         \
    \ Salut, à demain /
    -----------------
      \
        \
                      _ _
          | \__/|  .~    ~.
          /oo `./      .'
          {o__,   \    {
            / .  . )    \
            `-` '-' \    }
          .(   _(   )_.'
          '---.~_ _ _|
    ```

Dies zeigt, dass wir die Informationen in der Meta-Map verwenden konnten, um den Befehl im zweiten Schritt der Pipeline zu parametrisieren.

Allerdings war, wie oben bemerkt, ein Teil des beteiligten Codes etwas umständlich, da wir Metadaten entpacken mussten, während wir noch im Kontext des Workflow-Bodys waren.
Dieser Ansatz funktioniert gut, wenn man eine kleine Anzahl von Feldern aus der Meta-Map verwendet, würde aber schlecht skalieren, wenn wir viel mehr verwenden wollten.

Es gibt einen weiteren Operator namens `multiMap()`, der es uns ermöglicht, dies etwas zu rationalisieren, aber selbst dann ist es nicht ideal.

??? info "(Optional) Alternative Version mit `multiMap()`"

    Falls du dich fragst, wir konnten nicht einfach eine einzelne `map()`-Operation schreiben, die sowohl die `file` als auch den `character` ausgibt, weil das sie als Tupel zurückgeben würde.
    Wir mussten zwei separate `map()`-Operationen schreiben, um die `file`- und `character`-Elemente dem Prozess separat zu übergeben.

    Technisch gibt es eine andere Möglichkeit, dies durch eine einzelne Mapping-Operation zu tun, unter Verwendung des `multiMap()`-Operators, der mehrere Channels ausgeben kann.
    Zum Beispiel könntest du den Aufruf von `COWPY` oben durch folgenden Code ersetzen:

    === "Nachher"

        ```groovy title="main.nf" linenums="34"
            // Führe cowpy aus, um ASCII-Art zu generieren
            COWPY(
                ch_languages.multiMap { meta, file ->
                    file: file
                    character: meta.character
                }
            )
        ```

    === "Vorher"

        ```groovy title="main.nf" linenums="34"
            // Führe cowpy aus, um ASCII-Art zu generieren
            COWPY(
                ch_languages.map { meta, file -> file },
                ch_languages.map { meta, file -> meta.character }
            )
        ```

    Dies produziert exakt das gleiche Ergebnis.

In beiden Fällen ist es umständlich, dass wir etwas Entpacken auf Workflow-Ebene durchführen müssen.

Es wäre besser, wenn wir die gesamte Meta-Map in den Prozess einspeisen und dort auswählen könnten, was wir brauchen.

### 3.3. Die gesamte Meta-Map übergeben und verwenden

Der Sinn der Meta-Map ist es schließlich, alle Metadaten zusammen als ein Bündel zu übergeben.
Der einzige Grund, warum wir das oben nicht tun konnten, ist, dass der Prozess nicht für die Annahme einer Meta-Map eingerichtet ist.
Da wir aber den Prozesscode kontrollieren, können wir das ändern.

Lass uns den `COWPY`-Prozess modifizieren, um die `[meta, file]`-Tupelstruktur zu akzeptieren, die wir im ersten Prozess verwendet haben, damit wir den Workflow optimieren können.

Zu diesem Zweck müssen wir drei Dinge tun:

1. Die Eingabedefinitionen des `COWPY`-Prozessmoduls ändern
2. Den Prozessbefehl aktualisieren, um die Meta-Map zu verwenden
3. Den Prozessaufruf im Workflow-Body aktualisieren

Bereit? Los geht's!

#### 3.3.1. Die `COWPY`-Moduleingabe modifizieren

Nimm die folgenden Änderungen an der `cowpy.nf`-Moduldatei vor:

=== "Nachher"

    ```groovy title="cowpy.nf" linenums="10" hl_lines="2"
    input:
    tuple val(meta), path(input_file)
    ```

=== "Vorher"

    ```groovy title="cowpy.nf" linenums="10" hl_lines="2-3"
    input:
    path(input_file)
    val character
    ```

Dies ermöglicht es uns, die `[meta, file]`-Tupelstruktur zu verwenden, die wir zuvor im Tutorial behandelt haben.

Beachte, dass wir die Prozess-Ausgabedefinition nicht aktualisiert haben, um die Meta-Map auszugeben, um das Tutorial gestrafft zu halten, aber du kannst das gerne selbst als Übung nach dem Vorbild des `IDENTIFY_LANGUAGE`-Prozesses tun.

#### 3.3.2. Den Befehl aktualisieren, um das Meta-Map-Feld zu verwenden

Die gesamte Meta-Map ist jetzt innerhalb des Prozesses verfügbar, sodass wir direkt aus dem Befehlsblock heraus auf die darin enthaltenen Informationen verweisen können.

Nimm die folgenden Änderungen an der `cowpy.nf`-Moduldatei vor:

=== "Nachher"

    ```groovy title="cowpy.nf" linenums="16" hl_lines="3"
    script:
    """
    cat ${input_file} | cowpy -c ${meta.character} > cowpy-${input_file}
    """
    ```

=== "Vorher"

    ```groovy title="cowpy.nf" linenums="16" hl_lines="3"
    script:
    """
    cat ${input_file} | cowpy -c ${character} > cowpy-${input_file}
    """
    ```

Wir haben die Referenz auf den `character`-Wert, der zuvor als eigenständige Eingabe übergeben wurde, durch den in der Meta-Map gehaltenen Wert ersetzt, auf den wir mit `meta.character` verweisen.

Jetzt aktualisieren wir den Prozessaufruf entsprechend.

#### 3.3.3. Den Prozessaufruf aktualisieren und ausführen

Der Prozess erwartet jetzt, dass seine Eingabe die `[meta, file]`-Tupelstruktur verwendet, die der vorherige Prozess ausgibt, sodass wir einfach den gesamten `ch_languages`-Channel an den `COWPY`-Prozess übergeben können.

Nimm die folgenden Änderungen am Haupt-Workflow vor:

=== "Nachher"

    ```groovy title="main.nf" linenums="34" hl_lines="2"
    // Führe cowpy aus, um ASCII-Art zu generieren
    COWPY(ch_languages)
    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="34" hl_lines="3-4"
    // Führe cowpy aus, um ASCII-Art zu generieren
    COWPY(
        ch_languages.map { meta, file -> file },
        ch_languages.map { meta, file -> meta.character }
    )
    ```

Das vereinfacht den Aufruf erheblich!

Lass uns die Ergebnisse der vorherigen Ausführung löschen und es ausführen:

```bash
rm -r results
nextflow run main.nf
```

??? success "Befehlsausgabe"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [wise_sammet] DSL2 - revision: 99797b1e92

    executor >  local (14)
    [5d/dffd4e] process > IDENTIFY_LANGUAGE (7) [100%] 7 of 7 ✔
    [25/9243df] process > COWPY (7)             [100%] 7 of 7 ✔
    ```

Wenn du ins Ergebnisverzeichnis schaust, solltest du die gleichen Ausgaben wie zuvor sehen, d.h. einzelne Dateien, die die ASCII-Art jeder Begrüßung enthalten, gesprochen vom entsprechenden Charakter.

??? abstract "Verzeichnisinhalte"

    ```console
    ./results/
    ├── cowpy-bonjour.txt
    ├── cowpy-ciao.txt
    ├── cowpy-guten_tag.txt
    ├── cowpy-hallo.txt
    ├── cowpy-hello.txt
    ├── cowpy-hola.txt
    └── cowpy-salut.txt
    ```

Dies produziert also die gleichen Ergebnisse wie zuvor mit einfacherem Code.

Natürlich setzt dies voraus, dass du in der Lage bist, den Prozesscode zu modifizieren.
In manchen Fällen musst du dich möglicherweise auf bestehende Prozesse verlassen, die du nicht ändern darfst, was deine Optionen einschränkt.
Die gute Nachricht, wenn du planst Module aus dem [nf-core](https://nf-co.re/)-Projekt zu verwenden, ist, dass nf-core-Module alle standardmäßig die `[meta, file]`-Tupelstruktur verwenden.

### 3.4. Fehlerbehebung bei fehlenden erforderlichen Eingaben

Der `character`-Wert ist für den erfolgreichen Lauf des `COWPY`-Prozesses erforderlich.
Wenn wir keinen Standardwert in einer Konfigurationsdatei festlegen, MÜSSEN wir einen Wert im Datenblatt bereitstellen.

**Was passiert, wenn wir das nicht tun?**
Es hängt davon ab, was das Eingabedatenblatt enthält und welche Version des Workflows wir ausführen.

#### 3.4.1. Die Character-Spalte existiert, ist aber leer

Nehmen wir an, wir löschen den Character-Wert für einen der Einträge in unserem Datenblatt, um einen Datenerfassungsfehler zu simulieren:

```csv title="datasheet.csv" linenums="1" hl_lines="2"
id,character,recording
sampleA,,/workspaces/training/side-quests/metadata/data/bonjour.txt
sampleB,tux,/workspaces/training/side-quests/metadata/data/guten_tag.txt
sampleC,sheep,/workspaces/training/side-quests/metadata/data/hallo.txt
sampleD,turkey,/workspaces/training/side-quests/metadata/data/hello.txt
sampleE,stegosaurus,/workspaces/training/side-quests/metadata/data/hola.txt
sampleF,moose,/workspaces/training/side-quests/metadata/data/salut.txt
sampleG,turtle,/workspaces/training/side-quests/metadata/data/ciao.txt
```

Bei beiden Versionen des Workflows, die wir oben verwendet haben, wird der `character`-Schlüssel für alle Einträge erstellt, wenn das Datenblatt eingelesen wird, aber für `sampleA` wird der Wert eine leere Zeichenkette sein.

Dies wird einen Fehler verursachen.

??? failure "Befehlsausgabe"

    ```console hl_lines="8 11 16 28"
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [marvelous_hirsch] DSL2 - revision: 0dfeee3cc1

    executor >  local (9)
    [c1/c5dd4f] process > IDENTIFY_LANGUAGE (7) [ 85%] 6 of 7
    [d3/b7c415] process > COWPY (2)             [  0%] 0 of 6
    ERROR ~ Error executing process > 'COWPY (1)'

    Caused by:
      Process `COWPY (1)` terminated with an error exit status (2)


    Command executed:

      cat bonjour.txt | cowpy -c  > cowpy-bonjour.txt

    Command exit status:
      2

    Command output:
      (empty)

    Command error:
      usage: cowpy [-h] [-l] [-L] [-t] [-u] [-e EYES] [-c COWACTER] [-E] [-r] [-x]
                  [-C]
                  [msg ...]
      cowpy: error: argument -c/--cowacter: expected one argument

    Work dir:
      /workspaces/training/side-quests/metadata/work/ca/9d49796612a54dec5ed466063c809b

    Container:
      community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273

    Tip: you can try to figure out what's wrong by changing to the process work dir and showing the script file named `.command.sh`

    -- Check '.nextflow.log' file for details
    ```

Wenn Nextflow die `cowpy`-Befehlszeile für dieses Sample ausführt, wird `${meta.character}` mit einer leeren Zeichenkette in der `cowpy`-Befehlszeile gefüllt, sodass das `cowpy`-Tool einen Fehler wirft, der besagt, dass kein Wert für das `-c`-Argument bereitgestellt wurde.

#### 3.4.2. Die Character-Spalte existiert nicht im Datenblatt

Nehmen wir nun an, wir löschen die `character`-Spalte vollständig aus unserem Datenblatt:

```csv title="datasheet.csv" linenums="1"
id,recording
sampleA,/workspaces/training/side-quests/metadata/data/bonjour.txt
sampleB,/workspaces/training/side-quests/metadata/data/guten_tag.txt
sampleC,/workspaces/training/side-quests/metadata/data/hallo.txt
sampleD,/workspaces/training/side-quests/metadata/data/hello.txt
sampleE,/workspaces/training/side-quests/metadata/data/hola.txt
sampleF,/workspaces/training/side-quests/metadata/data/salut.txt
sampleG,/workspaces/training/side-quests/metadata/data/ciao.txt
```

In diesem Fall wird der `character`-Schlüssel beim Einlesen des Datenblatts überhaupt nicht erstellt.

##### 3.4.2.1. Wert wird auf Workflow-Ebene abgerufen

Wenn wir die Version des Codes verwenden, die wir in Abschnitt 3.2 geschrieben haben, wird Nextflow versuchen, VOR dem Aufruf des `COWPY`-Prozesses auf den `character`-Schlüssel in der Meta-Map zuzugreifen.

Es wird keine Elemente finden, die der Anweisung entsprechen, sodass `COWPY` überhaupt nicht ausgeführt wird.

??? success "Befehlsausgabe"

    ```console hl_lines="7"
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [desperate_montalcini] DSL2 - revision: 0dfeee3cc1

    executor >  local (7)
    [1a/df2544] process > IDENTIFY_LANGUAGE (7) [100%] 7 of 7 ✔
    [-        ] process > COWPY                 -
    ```

Aus Sicht von Nextflow lief dieser Workflow erfolgreich!
Allerdings werden keine der gewünschten Ausgaben produziert.

##### 3.4.2.2. Wert wird auf Prozessebene abgerufen

Wenn wir die Version aus Abschnitt 3.3 verwenden, wird Nextflow die gesamte Meta-Map an den `COWPY`-Prozess übergeben und versuchen, den Befehl auszuführen.

Dies wird einen Fehler verursachen, aber einen anderen als im ersten Fall.

??? failure "Befehlsausgabe"

    ```console hl_lines="8 11 16"
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [jovial_bohr] DSL2 - revision: eaaf375827

    executor >  local (9)
    [0d/ada9db] process > IDENTIFY_LANGUAGE (5) [ 85%] 6 of 7
    [06/28065f] process > COWPY (2)             [  0%] 0 of 6
    ERROR ~ Error executing process > 'COWPY (2)'

    Caused by:
      Process `COWPY (2)` terminated with an error exit status (1)


    Command executed:

      cat guten_tag.txt | cowpy -c null > cowpy-guten_tag.txt

    Command exit status:
      1

    Command output:
      (empty)

    Command error:
      Traceback (most recent call last):
        File "/opt/conda/bin/cowpy", line 10, in <module>
          sys.exit(main())
                  ~~~~^^
        File "/opt/conda/lib/python3.13/site-packages/cowpy/cow.py", line 1215, in main
          print(cow(eyes=args.eyes,
                ~~~^^^^^^^^^^^^^^^^
                tongue=args.tongue,
                ^^^^^^^^^^^^^^^^^^^
                thoughts=args.thoughts
                ^^^^^^^^^^^^^^^^^^^^^^
                    ).milk(msg)
                    ^
      TypeError: 'str' object is not callable

    Work dir:
      /workspaces/training/side-quests/metadata/work/06/28065f7d9fd7d22bba084aa941b6d6

    Container:
      community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273

    Tip: you can replicate the issue by changing to the process work dir and entering the command `bash .command.run`

    -- Check '.nextflow.log' file for details
    ```

Dies geschieht, weil `meta.character` nicht existiert, sodass unser Versuch, darauf zuzugreifen, `null` zurückgibt. Als Ergebnis setzt Nextflow buchstäblich `null` in die Befehlszeile ein, was natürlich vom `cowpy`-Tool nicht erkannt wird.

#### 3.4.3. Lösungen

Abgesehen von der Bereitstellung eines Standardwerts als Teil der Workflow-Konfiguration gibt es zwei Dinge, die wir tun können, um dies robuster zu handhaben:

1. Implementiere Eingabevalidierung in deinem Workflow, um sicherzustellen, dass das Datenblatt alle erforderlichen Informationen enthält. Eine [Einführung in die Eingabevalidierung](../hello_nf-core/05_input_validation.md) findest du im Hello nf-core Trainingskurs. <!-- TODO (future) pending a proper Validation side quest -->

2. Wenn du sicherstellen möchtest, dass jeder, der dein Prozessmodul verwendet, erforderliche Eingaben sofort identifizieren kann, kannst du die erforderliche Metadaten-Eigenschaft auch als explizite Eingabe definieren.

Hier ist ein Beispiel, wie das funktionieren würde.

Zunächst aktualisiere auf Prozessebene die Eingabedefinition wie folgt:

=== "Nachher"

    ```groovy title="cowpy.nf" linenums="12" hl_lines="2"
        input:
        tuple val(meta), val(character), path(input_file)
    ```

=== "Vorher"

    ```groovy title="cowpy.nf" linenums="12" hl_lines="2"
        input:
        tuple val(meta), path(input_file)
    ```

Verwende dann auf Workflow-Ebene eine Mapping-Operation, um die `character`-Eigenschaft aus den Metadaten zu extrahieren und sie zu einer expliziten Komponente des Eingabe-Tupels zu machen:

=== "Nachher"

    ```groovy title="main.nf" linenums="37" hl_lines="1"
        COWPY(ch_languages.map{meta, file -> [meta, meta.character, file]})
    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="37" hl_lines="1"
        COWPY(ch_languages)
    ```

Dieser Ansatz hat den Vorteil, dass er explizit zeigt, dass `character` erforderlich ist, und macht den Prozess einfacher in anderen Kontexten wiederverwendbar.

Dies unterstreicht ein wichtiges Designprinzip:

**Verwende die Meta-Map für optionale, beschreibende Informationen, aber extrahiere erforderliche Werte als explizite Eingaben.**

Die Meta-Map eignet sich hervorragend, um Channel-Strukturen sauber zu halten und willkürliche Channel-Strukturen zu vermeiden, aber für obligatorische Elemente, die direkt in einem Prozess referenziert werden, erzeugt das Extrahieren als explizite Eingaben robusteren und wartbareren Code.

### Fazit

In diesem Abschnitt hast du gelernt, wie du Metadaten nutzen kannst, um die Ausführung eines Prozesses anzupassen, indem du darauf entweder auf Workflow-Ebene oder auf Prozessebene zugreifst.

---

## Ergänzungsübung

Wenn du das Verwenden von Meta-Map-Informationen aus dem Inneren eines Prozesses üben möchtest, versuche andere Informationen aus der Meta-Map wie `lang` und `lang_group` zu verwenden, um anzupassen, wie die Ausgaben benannt und/oder organisiert werden.

Versuche zum Beispiel, den Code so zu modifizieren, dass dieses Ergebnis produziert wird:

```console title="Inhalt des Ergebnisverzeichnisses"
results/
├── germanic
│   ├── de-guten_tag.txt
│   ├── de-hallo.txt
│   └── en-hello.txt
└── romance
    ├── es-hola.txt
    ├── fr-bonjour.txt
    ├── fr-salut.txt
    └── it-ciao.txt
```

<!-- TODO (future) Provide worked out solution -->
<!-- the renaming should use the meta inside the process -->
<!-- the output org should use the meta in the workflow outputs -->

---

## Zusammenfassung

In dieser Side Quest hast du erkundet, wie du effektiv mit Metadaten in Nextflow-Workflows arbeiten kannst.

Dieses Muster, Metadaten explizit und mit den Daten verbunden zu halten, ist eine zentrale Best Practice in Nextflow und bietet mehrere Vorteile gegenüber dem Hardcodieren von Dateiinformationen:

- Dateimetadaten bleiben im gesamten Workflow mit den Dateien verknüpft
- Prozessverhalten kann pro Datei angepasst werden
- Ausgabeorganisation kann die Dateimetadaten widerspiegeln
- Dateiinformationen können während der Pipeline-Ausführung erweitert werden

Die Anwendung dieses Musters in deiner eigenen Arbeit wird es dir ermöglichen, robuste, wartbare Bioinformatik-Workflows zu erstellen.

### Wichtige Muster

1.  **Metadaten lesen und strukturieren:** CSV-Dateien lesen und organisierte Metadaten-Maps erstellen, die mit deinen Datendateien verknüpft bleiben.

    ```groovy
    channel.fromPath('datasheet.csv')
      .splitCsv(header: true)
      .map { row ->
          [ [id:row.id, character:row.character], row.recording ]
      }
    ```

2.  **Metadaten während des Workflows erweitern:** Neue Informationen zu deinen Metadaten hinzufügen, während deine Pipeline fortschreitet, indem du Prozessausgaben hinzufügst und Werte durch bedingte Logik ableitest.

    - Neue Schlüssel basierend auf Prozessausgaben hinzufügen

    ```groovy
    .map { meta, file, lang ->
      [ meta + [lang:lang], file ]
    }
    ```

    - Neue Schlüssel unter Verwendung einer bedingten Klausel hinzufügen

    ```groovy
    .map{ meta, file ->
        if ( meta.lang.equals("de") || meta.lang.equals('en') ){
            lang_group = "germanic"
        } else if ( meta.lang in ["fr", "es", "it"] ) {
            lang_group = "romance"
        } else {
            lang_group = "unknown"
        }
    }
    ```

3.  **Prozessverhalten anpassen:** Metadaten innerhalb des Prozesses verwenden.

    ```groovy
    cat $input_file | cowpy -c ${meta.character} > cowpy-${input_file}
    ```

### Zusätzliche Ressourcen

- [map](https://www.nextflow.io/docs/latest/operator.html#map)
- [stdout](https://www.nextflow.io/docs/latest/process.html#outputs)

---

## Was kommt als nächstes?

Kehre zum [Menü der Side Quests](./index.md) zurück oder klicke auf die Schaltfläche unten rechts auf der Seite, um zum nächsten Thema in der Liste zu gelangen.
