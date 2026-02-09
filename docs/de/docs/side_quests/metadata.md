# Metadaten und Meta-Maps

Bei jeder wissenschaftlichen Analyse arbeiten wir selten nur mit den reinen Datendateien.
Jede Datei kommt mit ihren eigenen zusätzlichen Informationen: was sie ist, woher sie stammt und was sie besonders macht.
Diese zusätzlichen Informationen nennen wir Metadaten.

Metadaten sind Daten, die andere Daten beschreiben.
Metadaten verfolgen wichtige Details über Dateien und experimentelle Bedingungen und helfen dabei, Analysen an die einzigartigen Eigenschaften jedes Datensatzes anzupassen.

Stell dir das wie einen Bibliothekskatalog vor: Während Bücher den eigentlichen Inhalt enthalten (Rohdaten), liefern die Katalogkarten wesentliche Informationen über jedes Buch – wann es veröffentlicht wurde, wer es geschrieben hat, wo man es findet (Metadaten).
In Nextflow-Pipelines können Metadaten verwendet werden, um:

- Dateispezifische Informationen durch den gesamten Workflow zu verfolgen
- Prozesse basierend auf Dateieigenschaften zu konfigurieren
- Zusammengehörige Dateien für gemeinsame Analysen zu gruppieren

### Lernziele

In diesem Side Quest werden wir erkunden, wie man mit Metadaten in Workflows umgeht.
Ausgehend von einer einfachen Datentabelle (in der Bioinformatik oft als Samplesheet bezeichnet), die grundlegende Dateiinformationen enthält, lernst du:

- Dateimetadaten aus CSV-Dateien zu lesen und zu parsen
- Metadaten-Maps zu erstellen und zu manipulieren
- Neue Metadatenfelder während der Workflow-Ausführung hinzuzufügen
- Metadaten zu verwenden, um das Prozessverhalten anzupassen

Diese Fähigkeiten helfen dir, robustere und flexiblere Pipelines zu erstellen, die komplexe Dateibeziehungen und Verarbeitungsanforderungen handhaben können.

### Voraussetzungen

Bevor du diesen Side Quest beginnst, solltest du:

- Das Tutorial [Hello Nextflow](../hello_nextflow/README.md) oder einen gleichwertigen Einsteigerkurs abgeschlossen haben.
- Dich sicher im Umgang mit grundlegenden Nextflow-Konzepten und -Mechanismen fühlen (Prozesse, Kanäle, Operatoren)

---

## 0. Erste Schritte

#### Öffne den Training-Codespace

Falls du das noch nicht getan hast, stelle sicher, dass du die Trainingsumgebung wie in [Umgebung einrichten](../envsetup/index.md) beschrieben öffnest.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Wechsle in das Projektverzeichnis

Lass uns in das Verzeichnis wechseln, in dem sich die Dateien für dieses Tutorial befinden.

```bash
cd side-quests/metadata
```

Du kannst VSCode so einstellen, dass es sich auf dieses Verzeichnis konzentriert:

```bash
code .
```

#### Überprüfe die Materialien

Du findest eine Haupt-Workflow-Datei und ein `data`-Verzeichnis, das eine Datentabelle und eine Handvoll Datendateien enthält.

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

Der Workflow in der `main.nf`-Datei ist ein Gerüst, das du schrittweise zu einem voll funktionsfähigen Workflow erweitern wirst.

Die Datentabelle listet die Pfade zu den Datendateien und einige zugehörige Metadaten auf, organisiert in 3 Spalten:

- `id`: selbsterklärend, eine ID, die der Datei gegeben wurde
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

Jede Datendatei enthält einen Begrüßungstext in einer von fünf Sprachen (fr: Französisch, de: Deutsch, es: Spanisch, it: Italienisch, en: Englisch).

Wir stellen dir auch ein containerisiertes Sprachanalyse-Tool namens `langid` zur Verfügung.

#### Überprüfe die Aufgabe

Deine Herausforderung besteht darin, einen Nextflow-Workflow zu schreiben, der:

1. Die Sprache in jeder Datei automatisch **identifiziert**
2. Dateien nach Sprachfamilie **gruppiert** (germanische vs. romanische Sprachen)
3. Die Verarbeitung für jede Datei basierend auf ihrer Sprache und ihren Metadaten **anpasst**
4. Ausgaben nach Sprachgruppe **organisiert**

Dies repräsentiert ein typisches Workflow-Muster, bei dem dateispezifische Metadaten Verarbeitungsentscheidungen steuern; genau die Art von Problem, die Metadaten-Maps elegant lösen.

#### Bereitschafts-Checkliste

Denkst du, du bist bereit loszulegen?

- [ ] Ich verstehe das Ziel dieses Kurses und seine Voraussetzungen
- [ ] Mein Codespace läuft
- [ ] Ich habe mein Arbeitsverzeichnis entsprechend eingestellt
- [ ] Ich verstehe die Aufgabe

Wenn du alle Kästchen abhaken kannst, kann es losgehen.

---

## 1. Metadaten aus einer Datentabelle laden

Öffne die `main.nf`-Workflow-Datei, um das Workflow-Gerüst zu untersuchen, das wir dir als Ausgangspunkt geben.

```groovy title="main.nf" linenums="1"
#!/usr/bin/env nextflow

workflow  {

    ch_datasheet = channel.fromPath("./data/datasheet.csv")

}
```

Du siehst, dass wir eine grundlegende Channel-Factory eingerichtet haben, um die Beispiel-Datentabelle als Datei zu laden, aber das liest noch nicht den Inhalt der Datei ein.
Lass uns damit beginnen, das hinzuzufügen.

### 1.1. Inhalt mit `splitCsv` einlesen

Wir müssen einen Operator wählen, der den Dateiinhalt mit minimalem Aufwand unsererseits angemessen parst.
Da unsere Datentabelle im CSV-Format vorliegt, ist dies ein Job für den [`splitCsv`](https://www.nextflow.io/docs/latest/reference/operator.html#splitcsv)-Operator, der jede Zeile in der Datei als Element im Kanal lädt.

Nimm die folgenden Änderungen vor, um eine `splitCsv()`-Operation zum Channel-Konstruktionscode hinzuzufügen, plus eine `view()`-Operation, um zu überprüfen, dass der Inhalt der Datei korrekt in den Kanal geladen wird.

=== "Danach"

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

Beachte, dass wir die Option `header: true` verwenden, um Nextflow mitzuteilen, dass die erste Zeile der CSV-Datei als Kopfzeile gelesen werden soll.

Lass uns schauen, was dabei herauskommt, oder?
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

Wir sehen, dass der Operator eine Map aus Schlüssel-Wert-Paaren für jede Zeile in der CSV-Datei konstruiert hat, mit den Spaltenüberschriften als Schlüssel für die entsprechenden Werte.

Jeder Map-Eintrag entspricht einer Spalte in unserer Datentabelle:

- `id`
- `character`
- `recording`

Das ist großartig! Es macht es einfach, auf bestimmte Felder aus jeder Datei zuzugreifen.
Zum Beispiel könnten wir auf die Datei-ID mit `id` oder den txt-Dateipfad mit `recording` zugreifen.

??? info "(Optional) Mehr über Maps"

    In Groovy, der Programmiersprache, auf der Nextflow aufbaut, ist eine Map eine Schlüssel-Wert-Datenstruktur ähnlich wie Dictionaries in Python, Objekte in JavaScript oder Hashes in Ruby.

    Hier ist ein ausführbares Skript, das zeigt, wie du eine Map definieren und auf ihren Inhalt zugreifen kannst:

    ```groovy title="examples/map_demo.nf"
    #!/usr/bin/env nextflow

    // Create a simple map
    def my_map = [id:'sampleA', character:'squirrel']

    // Print the whole map
    println "map: ${my_map}"

    // Access individual values using dot notation
    println "id: ${my_map.id}"
    println "character: ${my_map.character}"
    ```

    Obwohl es keinen richtigen `workflow`-Block hat, kann Nextflow dies ausführen, als wäre es ein Workflow:

    ```bash
    nextflow run examples/map_demo.nf
    ```

    Und hier ist, was du in der Ausgabe erwarten kannst:

    ```console title="Output"
     N E X T F L O W   ~  version 25.10.2

    Launching `map_demo.nf` [cheesy_plateau] DSL2 - revision: fae5b8496e

    map: [id:sampleA, character:squirrel]
    id: sampleA
    character: squirrel
    ```

### 1.2. Bestimmte Felder mit `map` auswählen

Nehmen wir an, wir möchten auf die `character`-Spalte aus der Datentabelle zugreifen und sie ausgeben.
Wir können den Nextflow-`map`-Operator verwenden, um über jedes Element in unserem Kanal zu iterieren und speziell den `character`-Eintrag aus dem Map-Objekt auszuwählen.

Nimm die folgenden Änderungen am Workflow vor:

=== "Danach"

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

Führe den Workflow jetzt erneut aus:

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

Erfolg! Wir haben die Map-Struktur, die aus unserer Datentabelle abgeleitet wurde, genutzt, um auf die Werte aus einzelnen Spalten für jede Zeile zuzugreifen.

Jetzt, da wir die Datentabelle erfolgreich eingelesen haben und Zugriff auf die Daten in jeder Zeile haben, können wir beginnen, unsere Pipeline-Logik zu implementieren.

### 1.3. Die Metadaten in einer 'Meta-Map' organisieren

Im aktuellen Zustand des Workflows befinden sich die Eingabedateien (unter dem `recording`-Schlüssel) und die zugehörigen Metadaten (`id`, `character`) alle auf der gleichen Ebene, als wären sie alle in einem großen Beutel.
Die praktische Konsequenz ist, dass jeder Prozess, der diesen Kanal konsumiert, mit dieser Struktur im Hinterkopf konfiguriert werden müsste:

```groovy
    input:
    tuple val(id), val(character), file(recording)
```

Das ist in Ordnung, solange sich die Anzahl der Spalten in der Datentabelle nicht ändert.
Wenn du jedoch auch nur eine Spalte zur Datentabelle hinzufügst, wird die Form des Kanals nicht mehr mit dem übereinstimmen, was der Prozess erwartet, und der Workflow wird Fehler produzieren.
Es macht den Prozess auch schwer mit anderen zu teilen, die möglicherweise leicht unterschiedliche Eingabedaten haben, und du musst möglicherweise Variablen in den Prozess hart codieren, die vom Script-Block nicht benötigt werden.

Um dieses Problem zu vermeiden, müssen wir einen Weg finden, die Kanalstruktur konsistent zu halten, unabhängig davon, wie viele Spalten die Datentabelle enthält.

Wir können das tun, indem wir alle Metadaten in einem Element innerhalb des Tupels sammeln, das wir die Metadaten-Map oder einfacher 'Meta-Map' nennen werden.

Nimm die folgenden Änderungen an der `map`-Operation vor:

=== "Danach"

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

Wir haben unsere Kanalelemente in ein Tupel umstrukturiert, das aus zwei Elementen besteht: der Meta-Map und dem entsprechenden Dateiobjekt.

Lass uns den Workflow ausführen:

```bash
nextflow run main.nf
```

??? success "Befehlsausgabe"

    ```console title="View meta map"
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

Jetzt enthält jedes Element im Kanal zuerst die Metadaten-Map und zweitens das entsprechende Dateiobjekt:

```console title="Example output structure"
[
  [id:sampleA, character:squirrel],
  /workspaces/training/side-quests/metadata/data/bonjour.txt
]
```

Dadurch wird das Hinzufügen weiterer Spalten in der Datentabelle mehr Metadaten in der `meta`-Map verfügbar machen, aber die Kanalform nicht ändern.
Dies ermöglicht es uns, Prozesse zu schreiben, die den Kanal konsumieren, ohne die Metadatenelemente in die Eingabespezifikation hart codieren zu müssen:

```groovy title="Syntax example"
    input:
    tuple val(meta), file(recording)
```

Dies ist eine weit verbreitete Konvention zur Organisation von Metadaten in Nextflow-Workflows.

### Fazit

In diesem Abschnitt hast du gelernt:

- **Warum Metadaten wichtig sind:** Das Beibehalten von Metadaten mit deinen Daten bewahrt wichtige Dateiinformationen während des gesamten Workflows.
- **Wie man Datentabellen einliest:** Verwendung von `splitCsv` zum Lesen von CSV-Dateien mit Kopfzeileninformationen und Transformation von Zeilen in strukturierte Daten
- **Wie man eine Meta-Map erstellt:** Trennung von Metadaten von Dateidaten unter Verwendung der Tupelstruktur `[ [id:value, ...], file ]`

---

## 2. Metadaten manipulieren

Jetzt, da wir unsere Metadaten geladen haben, lass uns etwas damit machen!

Wir werden ein Tool namens [`langid`](https://github.com/saffsd/langid.py) verwenden, um die Sprache zu identifizieren, die in der Aufnahmedatei jeder Kreatur enthalten ist.
Das Tool ist auf eine Reihe von Sprachen vortrainiert, und bei einem Textausschnitt gibt es eine Sprachvorhersage und einen zugehörigen Wahrscheinlichkeitswert aus, beide nach `stdout`.

### 2.1. Den Prozess importieren und den Code untersuchen

Wir stellen dir ein vorgeschriebenes Prozessmodul namens `IDENTIFY_LANGUAGE` zur Verfügung, das das `langid`-Tool umschließt, sodass du nur eine Include-Anweisung vor dem Workflow-Block hinzufügen musst.

Nimm die folgende Änderung am Workflow vor:

=== "Danach"

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

// Use langid to predict the language of each input file
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

Wie du sehen kannst, verwendet die Eingabedefinition die gleiche `tuple val(meta), path(file)`-Struktur, die wir gerade auf unseren Eingabekanal angewendet haben.

Die Ausgabedefinition ist als Tupel mit einer ähnlichen Struktur wie die Eingabe strukturiert, außer dass sie auch `stdout` als drittes Element enthält.
Dieses `tuple val(meta), path(file), <output>`-Muster hält die Metadaten sowohl mit den Eingabedaten als auch mit den Ausgaben verbunden, während sie durch die Pipeline fließen.

Beachte, dass wir hier den [`stdout`](https://www.nextflow.io/docs/latest/process.html#outputs)-Ausgabequalifizierer von Nextflow verwenden, weil das Tool seine Ausgabe direkt auf die Konsole ausgibt, anstatt eine Datei zu schreiben; und wir verwenden `sed` in der Befehlszeile, um den Wahrscheinlichkeitswert zu entfernen, den String durch Entfernen von Zeilenumbruchzeichen zu bereinigen und nur die Sprachvorhersage zurückzugeben.

### 2.2. Einen Aufruf zu `IDENTIFY_LANGUAGE` hinzufügen

Jetzt, da der Prozess für den Workflow verfügbar ist, können wir einen Aufruf zum `IDENTIFY_LANGUAGE`-Prozess hinzufügen, um ihn auf dem Datenkanal auszuführen.

Nimm die folgenden Änderungen am Workflow vor:

=== "Danach"

    ```groovy title="main.nf" linenums="7" hl_lines="7-9"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [[id: row.id, character: row.character], row.recording]
            }

        // Run langid to identify the language of each greeting
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

Beachte, dass wir die ursprüngliche `.view()`-Operation in der Kanalkonstruktion entfernt haben.

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

Ausgezeichnet! Wir haben jetzt eine Vorhersage, welche Sprache jeder Charakter spricht.

Und wie bereits erwähnt, haben wir auch die Eingabedatei und die Meta-Map in die Ausgabe aufgenommen, was bedeutet, dass beide mit den neuen Informationen, die wir gerade produziert haben, verbunden bleiben.
Dies wird sich im nächsten Schritt als nützlich erweisen.

!!! note

    Allgemeiner macht dieses Muster, die Meta-Map mit Ergebnissen verbunden zu halten, es einfacher, verwandte Ergebnisse zu verknüpfen, die dieselben Identifikatoren teilen.

    Wie du bereits gelernt hast, kannst du dich nicht auf die Reihenfolge der Elemente in Kanälen verlassen, um Ergebnisse über sie hinweg abzugleichen.
    Stattdessen musst du Schlüssel verwenden, um Daten korrekt zu verknüpfen, und Meta-Maps bieten eine ideale Struktur für diesen Zweck.

    Wir untersuchen diesen Anwendungsfall ausführlich im Side Quest [Splitting & Grouping](./splitting_and_grouping.md).

### 2.3. Metadaten mit Prozessausgaben erweitern

Da die Ergebnisse, die wir gerade produziert haben, selbst eine Form von Metadaten über den Inhalt der Dateien sind, wäre es nützlich, sie zu unserer Meta-Map hinzuzufügen.

Wir wollen jedoch die bestehende Meta-Map nicht direkt ändern.
Aus technischer Sicht ist es _möglich_, das zu tun, aber es ist unsicher.

Stattdessen erstellen wir eine neue Meta-Map, die den Inhalt der bestehenden Meta-Map plus ein neues `lang: lang_id`-Schlüssel-Wert-Paar enthält, das die neuen Informationen hält, unter Verwendung des `+`-Operators (ein Groovy-Feature).
Und wir kombinieren dies mit einer [`map`](https://www.nextflow.io/docs/latest/operator.html#map)-Operation, um die alte Map durch die neue zu ersetzen.

Hier sind die Änderungen, die du am Workflow vornehmen musst:

=== "Danach"

    ```groovy title="main.nf" linenums="13" hl_lines="3-7"
        // Run langid to identify the language of each greeting
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out
            .map { meta, file, lang_id ->
                [meta + [lang: lang_id], file]
            }
            .view()
    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="13" hl_lines="3"
        // Run langid to identify the language of each greeting
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out.view()
    ```

Wenn du noch nicht mit dem `+`-Operator vertraut bist oder wenn dies verwirrend erscheint, nimm dir ein paar Minuten Zeit, um die detaillierte Erklärung unten durchzugehen.

??? info "Erstellung der neuen Meta-Map mit dem `+`-Operator"

    **Zuerst musst du wissen, dass wir den Inhalt zweier Maps mit dem Groovy-Operator `+` zusammenführen können.**

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

    **Aber was, wenn du ein Feld hinzufügen musst, das noch nicht Teil einer Map ist?**

    Nehmen wir an, du beginnst wieder von `map1`, aber die Sprachvorhersage ist nicht in ihrer eigenen Map (es gibt keine `map2`).
    Stattdessen wird sie in einer Variablen namens `lang_id` gehalten, und du weißt, dass du ihren Wert (`'fr'`) mit dem Schlüssel `lang` speichern möchtest.

    Du kannst tatsächlich Folgendes tun:

    ```groovy
    new_map = [map1 + [lang: lang_id]]
    ```

    Hier erstellt `[lang: new_info]` spontan eine neue unbenannte Map, und `map1 + ` führt `map1` mit der neuen unbenannten Map zusammen, wodurch der gleiche `new_map`-Inhalt wie zuvor entsteht.

    Ordentlich, oder?

    **Jetzt lass uns das in den Kontext einer Nextflow-`channel.map()`-Operation übertragen.**

    Der Code wird zu:

    ```groovy
    .map { map1, lang_id ->
        [map1 + [lang: lang_id]]
    }
    ```

    Dies macht Folgendes:

    - `map1, lang_id ->` nimmt die zwei Elemente im Tupel
    - `[map1 + [lang: lang_id]]` erstellt die neue Map wie oben beschrieben

    Die Ausgabe ist eine einzelne unbenannte Map mit dem gleichen Inhalt wie `new_map` in unserem obigen Beispiel.
    Wir haben also effektiv transformiert:

    ```groovy
    [id: 'sampleA', character: 'squirrel'], 'it'
    ```

    in:

    ```groovy
    [id: 'sampleA', character: 'squirrel', lang: 'fr']
    ```

    Hoffentlich kannst du sehen, dass wenn wir `map1` in `meta` ändern, das im Grunde alles ist, was wir brauchen, um die Sprachvorhersage zu unserer Meta-Map in unserem Workflow hinzuzufügen.

    Außer einer Sache!

    Im Fall unseres Workflows **müssen wir auch die Anwesenheit des `file`-Objekts im Tupel berücksichtigen**, das aus `meta, file, lang_id` besteht.

    Der Code hier würde also zu:

    ```groovy
    .map { meta, file, lang_id ->
        [meta + [lang: lang_id], file]
    }
    ```

    Wenn du Schwierigkeiten hast zu verstehen, warum sich die `file` in der `map`-Operation zu bewegen scheint, stell dir vor, dass statt `[meta + [lang: lang_id], file]` diese Zeile `[new_map, file]` lautet.
    Dies sollte deutlicher machen, dass wir einfach die `file` an ihrer ursprünglichen Position an zweiter Stelle im Tupel belassen. Wir haben nur den `new_info`-Wert genommen und ihn in die Map gefaltet, die an erster Stelle steht.

    **Und das bringt uns zurück zur `tuple val(meta), path(file)`-Kanalstruktur!**

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
Wir haben die Ausgabe des Prozesses ordentlich von `meta, file, lang_id` reorganisiert, sodass `lang_id` jetzt einer der Schlüssel in der Meta-Map ist, und die Tupel des Kanals passen wieder zum `meta, file`-Modell.

<!-- TODO (future) Should we also show how to remove a key using subMap?! Or note where to find that. -->

### 2.4. Eine Sprachgruppe mit Bedingungen zuweisen

Jetzt, da wir unsere Sprachvorhersagen haben, lass uns die Informationen verwenden, um einige neue Gruppierungen zuzuweisen.

In unseren Beispieldaten können die von unseren Charakteren verwendeten Sprachen in germanische Sprachen (Englisch, Deutsch) und romanische Sprachen (Französisch, Spanisch, Italienisch) gruppiert werden.
Es könnte nützlich sein, diese Klassifizierung später in der Pipeline leicht verfügbar zu haben, also lass uns diese Information in der Meta-Map hinzufügen.

Und gute Nachrichten, dies ist ein weiterer Fall, der sich perfekt für die Verwendung des `map`-Operators eignet!

Konkret werden wir eine Variable namens `lang_group` definieren, etwas einfache bedingte Logik verwenden, um zu bestimmen, welchen Wert wir der `lang_group` für jedes Datenstück zuweisen sollen.

Die allgemeine Syntax wird so aussehen:

```groovy
.map { meta, file ->

    // conditional logic defining lang_group goes here

    [meta + [lang_group: lang_group], file]
}
```

Du siehst, dies ist sehr ähnlich zur spontanen Map-Zusammenführungsoperation, die wir im vorherigen Schritt verwendet haben.
Wir müssen nur die bedingten Anweisungen schreiben.

Hier ist die bedingte Logik, die wir anwenden möchten:

- Definiere eine Variable namens `lang_group` mit dem Standardwert `'unknown'`.
- Wenn `lang` entweder Deutsch (`'de'`) oder Englisch (`'en'`) ist, ändere `lang_group` zu `germanic`.
- Sonst, wenn `lang` in einer Liste enthalten ist, die Französisch (`'fr'`), Spanisch (`'es'`) und Italienisch (`'it'`) enthält, ändere `lang_group` zu `romance`.

Versuche selbst, es zu schreiben, wenn du bereits weißt, wie man bedingte Anweisungen in Nextflow schreibt.

!!! tip

    Du kannst auf den Wert von `lang` innerhalb der Map-Operation mit `meta.lang` zugreifen.

Du solltest am Ende die folgenden Änderungen am Workflow vornehmen:

=== "Danach"

    ```groovy title="main.nf" linenums="13" hl_lines="7-19"
        // Run langid to identify the language of each greeting
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
        // Run langid to identify the language of each greeting
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out
            .map { meta, file, lang_id ->
                [meta + [lang: lang_id], file]
            }
            .view()
    ```

Hier sind die wichtigsten Punkte:

- Wir verwenden `def lang_group = "unknown"`, um die `lang_group`-Variable mit dem Standardwert `unknown` zu erstellen.
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

Wie du sehen kannst, behalten die Kanalelemente ihre `[meta, file]`-Struktur bei, aber die Meta-Map enthält jetzt diese neue Klassifizierung.

### Fazit

In diesem Abschnitt hast du gelernt, wie man:

- **Eingabemetadaten auf Ausgabekanäle anwendet**: Das Kopieren von Metadaten auf diese Weise ermöglicht es uns, Ergebnisse später basierend auf Metadateninhalten zu verknüpfen.
- **Benutzerdefinierte Schlüssel erstellt**: Du hast zwei neue Schlüssel in deiner Meta-Map erstellt und sie mit `meta + [new_key:value]` in die bestehende Meta-Map eingefügt. Einen basierend auf einem berechneten Wert aus einem Prozess und einen basierend auf einer Bedingung, die du im `map`-Operator festgelegt hast.

Diese ermöglichen es dir, neue und bestehende Metadaten mit Dateien zu verknüpfen, während du durch deine Pipeline fortschreitest.
Selbst wenn du Metadaten nicht als Teil eines Prozesses verwendest, macht es das Beibehalten der Meta-Map mit den Daten auf diese Weise einfach, alle relevanten Informationen zusammenzuhalten.

---

## 3. Meta-Map-Informationen in einem Prozess verwenden

Jetzt, da du weißt, wie man die Meta-Map erstellt und aktualisiert, können wir zum wirklich spaßigen Teil kommen: die Metadaten tatsächlich in einem Prozess zu verwenden.

Genauer gesagt werden wir einen zweiten Schritt zu unserem Workflow hinzufügen, um jedes Tier als ASCII-Art zu zeichnen und es den aufgezeichneten Text in einer Sprechblase sagen zu lassen.
Wir werden dies mit einem Tool namens [`cowpy`](https://github.com/jeffbuttars/cowpy) tun.

??? info "Was macht `cowpy`?"

    `cowpy` ist ein Kommandozeilen-Tool, das ASCII-Art generiert, um beliebige Texteingaben auf unterhaltsame Weise anzuzeigen.
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

### 3.1. Den Prozess importieren und den Code untersuchen

Wir stellen dir ein vorgeschriebenes Prozessmodul namens `COWPY` zur Verfügung, das das `cowpy`-Tool umschließt, sodass du nur eine Include-Anweisung vor dem Workflow-Block hinzufügen musst.

Nimm die folgende Änderung am Workflow vor:

=== "Danach"

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

// Generate ASCII art with cowpy
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

Wie du sehen kannst, ist dieser Prozess derzeit so konzipiert, dass er eine Eingabedatei (die den anzuzeigenden Text enthält) und einen Wert, der den Charakter angibt, der in ASCII-Art gezeichnet werden soll, entgegennimmt, normalerweise auf Workflow-Ebene durch einen Befehlszeilenparameter bereitgestellt.

### 3.2. Ein Meta-Map-Feld als Eingabe übergeben

Als wir das `cowpy`-Tool im Hello Nextflow-Kurs verwendet haben, haben wir einen Befehlszeilenparameter verwendet, um zu bestimmen, welcher Charakter zum Zeichnen des endgültigen Bildes verwendet werden soll.
Das machte Sinn, weil wir nur ein Bild pro Durchlauf der Pipeline generierten.

In diesem Tutorial möchten wir jedoch ein geeignetes Bild für jedes Subjekt generieren, das wir verarbeiten, sodass die Verwendung eines Befehlszeilenparameters zu einschränkend wäre.

Gute Nachrichten: Wir haben eine `character`-Spalte in unserer Datentabelle und daher in unserer Meta-Map.
Lass uns das verwenden, um den Charakter festzulegen, den der Prozess für jeden Eintrag verwenden soll.

Dazu müssen wir drei Dinge tun:

1. Dem Ausgabekanal, der aus dem vorherigen Prozess kommt, einen Namen geben, damit wir bequemer damit arbeiten können.
2. Bestimmen, wie auf die interessierenden Informationen zugegriffen werden kann
3. Einen Aufruf zum zweiten Prozess hinzufügen und die Informationen entsprechend einspeisen.

Lass uns anfangen.

#### 3.2.1. Den vorherigen Ausgabekanal benennen

Wir haben die vorherigen Manipulationen direkt auf dem Ausgabekanal des ersten Prozesses angewendet, `IDENTIFY_LANGUAGE.out`.
Um den Inhalt dieses Kanals an den nächsten Prozess zu übergeben (und dies auf eine Weise zu tun, die klar und leicht zu lesen ist), möchten wir ihm einen eigenen Namen geben, `ch_languages`.

Wir können das mit dem [`set`](https://www.nextflow.io/docs/latest/reference/operator.html#set)-Operator tun.

Ersetze im Haupt-Workflow den `.view()`-Operator durch `.set { ch_languages }` und füge eine Zeile hinzu, die testet, dass wir auf den Kanal mit Namen verweisen können.

=== "Danach"

    ```groovy title="main.nf" linenums="14" hl_lines="19 21 22"
        // Run langid to identify the language of each greeting
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

        // Temporary: peek into ch_languages
        ch_languages.view()
    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="14" hl_lines="19"
        // Run langid to identify the language of each greeting
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

Dies bestätigt, dass wir jetzt auf den Kanal mit Namen verweisen können.

#### 3.2.2. Auf die Datei und Charakter-Metadaten zugreifen

Wir wissen aus dem Betrachten des Modulcodes, dass der `COWPY`-Prozess erwartet, eine Textdatei und einen `character`-Wert zu erhalten.
Um den Aufruf zum `COWPY`-Prozess zu schreiben, müssen wir nur wissen, wie wir das entsprechende Dateiobjekt und die Metadaten aus jedem Element im Kanal extrahieren können.

Wie oft ist der einfachste Weg, das zu tun, eine `map`-Operation zu verwenden.

Unser Kanal enthält Tupel, die als `[meta, file]` strukturiert sind, sodass wir direkt auf das `file`-Objekt zugreifen können, und wir können auf den `character`-Wert zugreifen, der in der Meta-Map gespeichert ist, indem wir ihn als `meta.character` bezeichnen.

Nimm im Haupt-Workflow die folgenden Codeänderungen vor:

=== "Danach"

    ```groovy title="main.nf" linenums="34"
        // Temporary: access the file and character
        ch_languages.map { meta, file -> file }.view { file -> "File: " + file }
        ch_languages.map { meta, file -> meta.character }.view { character -> "Character: " + character }
    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="34"
        // Temporary: peek into ch_languages
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

_Die Dateipfade und Charakterwerte können in deiner Ausgabe in einer anderen Reihenfolge erscheinen._

Dies bestätigt, dass wir auf die Datei und den Charakter für jedes Element im Kanal zugreifen können.

#### 3.2.3. Den `COWPY`-Prozess aufrufen

Jetzt lass uns alles zusammenfügen und tatsächlich den `COWPY`-Prozess auf dem `ch_languages`-Kanal aufrufen.

Nimm im Haupt-Workflow die folgenden Codeänderungen vor:

=== "Danach"

    ```groovy title="main.nf" linenums="34"
        // Run cowpy to generate ASCII art
        COWPY(
            ch_languages.map { meta, file -> file },
            ch_languages.map { meta, file -> meta.character }
        )
    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="34"
        // Temporary: access the file and character
        ch_languages.map { meta, file -> [file, meta.character] }
            .view()
    ```

Du siehst, wir kopieren einfach die zwei Map-Operationen (minus die `.view()`-Anweisungen) als Eingaben für den Prozessaufruf.
Stelle nur sicher, dass du das Komma zwischen ihnen nicht vergisst!

Es ist ein bisschen umständlich, aber wir werden im nächsten Abschnitt sehen, wie man das besser machen kann.

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

Wenn du im Ergebnisverzeichnis nachschaust, solltest du die einzelnen Dateien sehen, die die ASCII-Art jeder Begrüßung enthalten, gesprochen vom entsprechenden Charakter.

??? abstract "Verzeichnis- und Beispieldateiinhalt"

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

Wie jedoch oben erwähnt, war ein Teil des beteiligten Codes etwas umständlich, da wir Metadaten auspacken mussten, während wir uns noch im Kontext des Workflow-Bodys befanden.
Dieser Ansatz funktioniert gut für die Verwendung einer kleinen Anzahl von Feldern aus der Meta-Map, würde aber schlecht skalieren, wenn wir viel mehr verwenden wollten.

Es gibt einen anderen Operator namens `multiMap()`, der es uns ermöglicht, dies ein wenig zu optimieren, aber selbst dann ist es nicht ideal.

??? info "(Optional) Alternative Version mit `multiMap()`"

    Falls du dich fragst, wir konnten nicht einfach eine einzelne `map()`-Operation schreiben, die sowohl die `file` als auch den `character` ausgibt, weil das sie als Tupel zurückgeben würde.
    Wir mussten zwei separate `map()`-Operationen schreiben, um die `file`- und `character`-Elemente separat an den Prozess zu übergeben.

    Technisch gibt es einen anderen Weg, dies durch eine einzelne Mapping-Operation zu tun, unter Verwendung des `multiMap()`-Operators, der in der Lage ist, mehrere Kanäle auszugeben.
    Zum Beispiel könntest du den Aufruf zu `COWPY` oben durch den folgenden Code ersetzen:

    === "Danach"

        ```groovy title="main.nf" linenums="34"
            // Run cowpy to generate ASCII art
            COWPY(
                ch_languages.multiMap { meta, file ->
                    file: file
                    character: meta.character
                }
            )
        ```

    === "Vorher"

        ```groovy title="main.nf" linenums="34"
            // Run cowpy to generate ASCII art
            COWPY(
                ch_languages.map { meta, file -> file },
                ch_languages.map { meta, file -> meta.character }
            )
        ```

    Dies produziert genau das gleiche Ergebnis.

In beiden Fällen ist es umständlich, dass wir etwas Auspacken auf Workflow-Ebene machen müssen.

Es wäre besser, wenn wir die gesamte Meta-Map in den Prozess einspeisen und dort auswählen könnten, was wir brauchen.

### 3.3. Die gesamte Meta-Map übergeben und verwenden

Der Sinn der Meta-Map ist schließlich, alle Metadaten zusammen als Bündel zu übergeben.
Der einzige Grund, warum wir das oben nicht tun konnten, ist, dass der Prozess nicht so eingerichtet ist, dass er eine Meta-Map akzeptiert.
Aber da wir den Prozesscode kontrollieren, können wir das ändern.

Lass uns den `COWPY`-Prozess so modifizieren, dass er die `[meta, file]`-Tupelstruktur akzeptiert, die wir im ersten Prozess verwendet haben, damit wir den Workflow optimieren können.

Dazu müssen wir drei Dinge tun:

1. Die Eingabedefinitionen des `COWPY`-Prozessmoduls modifizieren
2. Den Prozessbefehl aktualisieren, um die Meta-Map zu verwenden
3. Den Prozessaufruf im Workflow-Body aktualisieren

Bereit? Los geht's!

#### 3.3.1. Die `COWPY`-Moduleingabe modifizieren

Nimm die folgenden Änderungen an der `cowpy.nf`-Moduldatei vor:

=== "Danach"

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

Dies ermöglicht es uns, die `[meta, file]`-Tupelstruktur zu verwenden, die wir früher im Tutorial behandelt haben.

Beachte, dass wir die Prozessausgabedefinition nicht aktualisiert haben, um die Meta-Map auszugeben, um das Tutorial übersichtlich zu halten, aber fühle dich frei, das selbst als Übung zu tun, indem du dem Modell des `IDENTIFY_LANGUAGE`-Prozesses folgst.

#### 3.3.2. Den Befehl aktualisieren, um das Meta-Map-Feld zu verwenden

Die gesamte Meta-Map ist jetzt innerhalb des Prozesses verfügbar, sodass wir direkt aus dem Befehlsblock auf die darin enthaltenen Informationen verweisen können.

Nimm die folgenden Änderungen an der `cowpy.nf`-Moduldatei vor:

=== "Danach"

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

Wir haben die Referenz auf den `character`-Wert, der zuvor als eigenständige Eingabe übergeben wurde, durch den Wert ersetzt, der in der Meta-Map gehalten wird, auf den wir mit `meta.character` verweisen.

Jetzt lass uns den Prozessaufruf entsprechend aktualisieren.

#### 3.3.3. Den Prozessaufruf aktualisieren und ausführen

Der Prozess erwartet jetzt, dass seine Eingabe die `[meta, file]`-Tupelstruktur verwendet, was der vorherige Prozess ausgibt, sodass wir einfach den gesamten `ch_languages`-Kanal an den `COWPY`-Prozess übergeben können.

Nimm die folgenden Änderungen am Haupt-Workflow vor:

=== "Danach"

    ```groovy title="main.nf" linenums="34" hl_lines="2"
    // Run cowpy to generate ASCII art
    COWPY(ch_languages)
    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="34" hl_lines="3-4"
    // Run cowpy to generate ASCII art
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

Wenn du im Ergebnisverzeichnis nachschaust, solltest du die gleichen Ausgaben wie zuvor sehen, _d.h._ einzelne Dateien, die die ASCII-Art jeder Begrüßung enthalten, gesprochen vom entsprechenden Charakter.

??? abstract "Verzeichnisinhalt"

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
In einigen Fällen musst du dich möglicherweise auf bestehende Prozesse verlassen, die du nicht frei modifizieren kannst, was deine Optionen einschränkt.
Die gute Nachricht, wenn du planst, Module aus dem [nf-core](https://nf-co.re/)-Projekt zu verwenden, ist, dass nf-core-Module alle so eingerichtet sind, dass sie die `[meta, file]`-Tupelstruktur als Standard verwenden.

### 3.4. Fehlerbehebung bei fehlenden erforderlichen Eingaben

Der `character`-Wert ist erforderlich, damit der `COWPY`-Prozess erfolgreich ausgeführt werden kann.
Wenn wir keinen Standardwert dafür in einer Konfigurationsdatei festlegen, MÜSSEN wir einen Wert dafür in der Datentabelle angeben.

**Was passiert, wenn wir das nicht tun?**
Es hängt davon ab, was die Eingabedatentabelle enthält und welche Version des Workflows wir ausführen.

#### 3.4.1. Die Charakterspalte existiert, ist aber leer

Nehmen wir an, wir löschen den Charakterwert für einen der Einträge in unserer Datentabelle, um einen Datenerfassungsfehler zu simulieren:

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

Für beide Versionen des Workflows, die wir oben verwendet haben, wird der `character`-Schlüssel für alle Einträge erstellt, wenn die Datentabelle eingelesen wird, aber für `sampleA` wird der Wert ein leerer String sein.

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

Wenn Nextflow die `cowpy`-Befehlszeile für diese Probe ausführt, wird `${meta.character}` mit einem leeren String in der `cowpy`-Befehlszeile gefüllt, sodass das `cowpy`-Tool einen Fehler wirft, der besagt, dass kein Wert für das `-c`-Argument angegeben wurde.

#### 3.4.2. Die Charakterspalte existiert nicht in der Datentabelle

Nehmen wir jetzt an, wir löschen die `character`-Spalte vollständig aus unserer Datentabelle:

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

In diesem Fall wird der `character`-Schlüssel überhaupt nicht erstellt, wenn die Datentabelle eingelesen wird.

##### 3.4.2.1. Wert auf Workflow-Ebene zugegriffen

Wenn wir die Version des Codes verwenden, die wir in Abschnitt 3.2 geschrieben haben, wird Nextflow versuchen, auf den `character`-Schlüssel in der Meta-Map zuzugreifen, BEVOR der `COWPY`-Prozess aufgerufen wird.

Es wird keine Elemente finden, die der Anweisung entsprechen, sodass es `COWPY` überhaupt nicht ausführen wird.

??? success "Befehlsausgabe"

    ```console hl_lines="7"
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [desperate_montalcini] DSL2 - revision: 0dfeee3cc1

    executor >  local (7)
    [1a/df2544] process > IDENTIFY_LANGUAGE (7) [100%] 7 of 7 ✔
    [-        ] process > COWPY                 -
    ```

Soweit es Nextflow betrifft, wurde dieser Workflow erfolgreich ausgeführt!
Allerdings werden keine der gewünschten Ausgaben produziert.

##### 3.4.2.2. Wert auf Prozessebene zugegriffen

Wenn wir die Version in Abschnitt 3.3 verwenden, wird Nextflow die gesamte Meta-Map an den `COWPY`-Prozess übergeben und versuchen, den Befehl auszuführen.

Dies wird einen Fehler verursachen, aber einen anderen im Vergleich zum ersten Fall.

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

Dies geschieht, weil `meta.character` nicht existiert, sodass unser Versuch, darauf zuzugreifen, `null` zurückgibt. Als Ergebnis fügt Nextflow buchstäblich `null` in die Befehlszeile ein, was natürlich vom `cowpy`-Tool nicht erkannt wird.

#### 3.4.3. Lösungen

Abgesehen von der Bereitstellung eines Standardwerts als Teil der Workflow-Konfiguration gibt es zwei Dinge, die wir tun können, um dies robuster zu handhaben:

1. Implementiere Eingabevalidierung in deinem Workflow, um sicherzustellen, dass die Datentabelle alle erforderlichen Informationen enthält. Du findest eine [Einführung in die Eingabevalidierung](../hello_nf-core/05_input_validation.md) im Hello nf-core-Trainingskurs. <!-- TODO (future) pending a proper Validation side quest -->

2. Wenn du sicherstellen möchtest, dass jeder, der dein Prozessmodul verwendet, erforderliche Eingaben sofort identifizieren kann, kannst du auch die erforderliche Metadateneigenschaft zu einer expliziten Eingabe machen.

Hier ist ein Beispiel, wie das funktionieren würde.

Zuerst aktualisiere auf Prozessebene die Eingabedefinition wie folgt:

=== "Danach"

    ```groovy title="cowpy.nf" linenums="12" hl_lines="2"
        input:
        tuple val(meta), val(character), path(input_file)
    ```

=== "Vorher"

    ```groovy title="cowpy.nf" linenums="12" hl_lines="2"
        input:
        tuple val(meta), path(input_file)
    ```

Dann verwende auf Workflow-Ebene eine Mapping-Operation, um die `character`-Eigenschaft aus den Metadaten zu extrahieren und sie zu einer expliziten Komponente des Eingabetupels zu machen:

=== "Danach"

    ```groovy title="main.nf" linenums="37" hl_lines="1"
        COWPY(ch_languages.map{meta, file -> [meta, meta.character, file]})
    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="37" hl_lines="1"
        COWPY(ch_languages)
    ```

Dieser Ansatz hat den Vorteil, dass er explizit zeigt, dass `character` erforderlich ist, und macht den Prozess einfacher in anderen Kontexten wiederverwendbar.

Dies hebt ein wichtiges Designprinzip hervor:

**Verwende die Meta-Map für optionale, beschreibende Informationen, aber extrahiere erforderliche Werte als explizite Eingaben.**

Die Meta-Map ist ausgezeichnet, um Kanalstrukturen sauber zu halten und willkürliche Kanalstrukturen zu verhindern, aber für obligatorische Elemente, die direkt in einem Prozess referenziert werden, schafft das Extrahieren als explizite Eingaben robusteren und wartbareren Code.

### Fazit

In diesem Abschnitt hast du gelernt, wie man Metadaten verwendet, um die Ausführung eines Prozesses anzupassen, indem man entweder auf Workflow-Ebene oder auf Prozessebene darauf zugreift.

---

## Ergänzende Übung

Wenn du das Verwenden von Meta-Map-Informationen innerhalb eines Prozesses üben möchtest, versuche, andere Informationen aus der Meta-Map wie `lang` und `lang_group` zu verwenden, um anzupassen, wie die Ausgaben benannt und/oder organisiert werden.

Versuche zum Beispiel, den Code so zu modifizieren, dass dieses Ergebnis produziert wird:

```console title="Results directory contents"
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

In diesem Side Quest hast du erkundet, wie man effektiv mit Metadaten in Nextflow-Workflows arbeitet.

Dieses Muster, Metadaten explizit und mit den Daten verbunden zu halten, ist eine zentrale Best Practice in Nextflow und bietet mehrere Vorteile gegenüber dem Hardcodieren von Dateiinformationen:

- Dateimetadaten bleiben während des gesamten Workflows mit Dateien verbunden
- Das Prozessverhalten kann pro Datei angepasst werden
- Die Ausgabeorganisation kann Dateimetadaten widerspiegeln
- Dateiinformationen können während der Pipeline-Ausführung erweitert werden

Die Anwendung dieses Musters in deiner eigenen Arbeit wird es dir ermöglichen, robuste, wartbare Bioinformatik-Workflows zu erstellen.

### Wichtige Muster

1.  **Metadaten lesen und strukturieren:** CSV-Dateien lesen und organisierte Metadaten-Maps erstellen, die mit deinen Datendateien verbunden bleiben.

    ```groovy
    channel.fromPath('datasheet.csv')
      .splitCsv(header: true)
      .map { row ->
          [ [id:row.id, character:row.character], row.recording ]
      }
    ```

2.  **Metadaten während des Workflows erweitern:** Neue Informationen zu deinen Metadaten hinzufügen, während deine Pipeline fortschreitet, indem Prozessausgaben hinzugefügt und Werte durch bedingte Logik abgeleitet werden.

    - Neue Schlüssel basierend auf Prozessausgabe hinzufügen

    ```groovy
    .map { meta, file, lang ->
      [ meta + [lang:lang], file ]
    }
    ```

    - Neue Schlüssel mit einer bedingten Klausel hinzufügen

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

## Wie geht es weiter?

Kehre zum [Menü der Side Quests](./index.md) zurück oder klicke auf den Button unten rechts auf der Seite, um zum nächsten Thema in der Liste zu gelangen.
