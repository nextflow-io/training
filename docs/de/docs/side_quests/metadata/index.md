# Metadaten und Meta-Maps

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

In jeder wissenschaftlichen Analyse arbeiten wir selten nur mit den rohen Datendateien.
Jede Datei bringt ihre eigenen zusätzlichen Informationen mit: was sie ist, woher sie stammt und was sie besonders macht.
Diese zusätzlichen Informationen nennen wir Metadaten.

Metadaten sind Daten, die andere Daten beschreiben.
Metadaten erfassen wichtige Details über Dateien und experimentelle Bedingungen und helfen dabei, Analysen an die einzigartigen Eigenschaften jedes Datensatzes anzupassen.

Stell dir das wie einen Bibliothekskatalog vor: Während Bücher den eigentlichen Inhalt enthalten (Rohdaten), liefern die Katalogkarten wesentliche Informationen über jedes Buch – wann es veröffentlicht wurde, wer es geschrieben hat, wo es zu finden ist (Metadaten).
In Nextflow-Pipelines können Metadaten verwendet werden, um:

- Dateispezifische Informationen durch den gesamten Workflow zu verfolgen
- Prozesse basierend auf Dateieigenschaften zu konfigurieren
- Verwandte Dateien für gemeinsame Analysen zu gruppieren

### Lernziele

In dieser Side Quest erkunden wir, wie man Metadaten in Workflows verarbeitet.
Ausgehend von einem einfachen Datenblatt (in der Bioinformatik oft als Samplesheet bezeichnet), das grundlegende Dateiinformationen enthält, lernst du:

- Datei-Metadaten aus CSV-Dateien zu lesen und zu verarbeiten
- Meta-Maps zu erstellen und zu bearbeiten
- Neue Metadatenfelder während der Workflow-Ausführung hinzuzufügen
- Metadaten zur Anpassung des Prozessverhaltens zu nutzen

Diese Fähigkeiten helfen dir, robustere und flexiblere Pipelines zu entwickeln, die komplexe Dateibeziehungen und Verarbeitungsanforderungen bewältigen können.

### Voraussetzungen

Bevor du diese Side Quest angehst, solltest du:

- Das Tutorial [Hello Nextflow](../hello_nextflow/README.md) oder einen gleichwertigen Einsteigerkurs abgeschlossen haben.
- Mit grundlegenden Nextflow-Konzepten und -Mechanismen (Prozesse, Kanäle, Operatoren) vertraut sein.

---

## 0. Erste Schritte

#### Öffne den Training-Codespace

Falls du das noch nicht getan hast, öffne die Trainingsumgebung wie in der [Umgebung einrichten](../envsetup/index.md)-Anleitung beschrieben.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Wechsle in das Projektverzeichnis

Wechseln wir in das Verzeichnis, in dem sich die Dateien für dieses Tutorial befinden.

```bash
cd side-quests/metadata
```

Du kannst VSCode so einstellen, dass es sich auf dieses Verzeichnis konzentriert:

```bash
code .
```

#### Schau dir die Materialien an

Du findest eine Haupt-Workflow-Datei und ein `data`-Verzeichnis mit einem Datenblatt und einigen Datendateien.

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

Der Workflow in der `main.nf`-Datei ist ein Grundgerüst, das du schrittweise zu einem vollständig funktionierenden Workflow ausbauen wirst.

Das Datenblatt listet die Pfade zu den Datendateien und einige zugehörige Metadaten auf, gegliedert in 3 Spalten:

- `id`: selbsterklärend, eine ID für die Datei
- `character`: ein Charaktername, den wir später verwenden, um verschiedene Figuren zu zeichnen
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

Außerdem stellen wir dir ein containerisiertes Sprachanalyse-Tool namens `langid` zur Verfügung.

#### Schau dir die Aufgabe an

Deine Aufgabe ist es, einen Nextflow-Workflow zu schreiben, der:

1. Die Sprache in jeder Datei **automatisch erkennt**
2. Dateien nach Sprachfamilie **gruppiert** (Germanisch vs. Romanisch)
3. Die Verarbeitung jeder Datei basierend auf ihrer Sprache und ihren Metadaten **anpasst**
4. Ausgaben nach Sprachgruppe **organisiert**

Dies ist ein typisches Workflow-Muster, bei dem dateispezifische Metadaten Verarbeitungsentscheidungen steuern – genau die Art von Problem, die Meta-Maps elegant lösen.

#### Bereitschafts-Checkliste

Bereit zum Eintauchen?

- [ ] Ich verstehe das Ziel dieses Kurses und seine Voraussetzungen
- [ ] Mein Codespace läuft
- [ ] Ich habe mein Arbeitsverzeichnis entsprechend eingestellt
- [ ] Ich verstehe die Aufgabe

Wenn du alle Punkte abhaken kannst, kann es losgehen.

---

## 1. Metadaten aus einem Datenblatt laden

Öffne die `main.nf`-Workflow-Datei, um das Workflow-Grundgerüst zu untersuchen, das wir dir als Ausgangspunkt geben.

```groovy title="main.nf" linenums="1"
#!/usr/bin/env nextflow

workflow  {

    ch_datasheet = channel.fromPath("./data/datasheet.csv")

}
```

Wir haben eine einfache Channel-Factory eingerichtet, um das Beispiel-Datenblatt als Datei zu laden, aber das liest den Inhalt der Datei noch nicht ein.
Fangen wir damit an.

### 1.1. Inhalt mit `splitCsv` einlesen

Wir müssen einen Operator wählen, der den Dateiinhalt mit minimalem Aufwand geeignet verarbeitet.
Da unser Datenblatt im CSV-Format vorliegt, ist das eine Aufgabe für den [`splitCsv`](https://www.nextflow.io/docs/latest/reference/operator.html#splitcsv)-Operator, der jede Zeile der Datei als Element im Kanal lädt.

Nimm die folgenden Änderungen vor, um eine `splitCsv()`-Operation zum Channel-Konstruktionscode hinzuzufügen, sowie eine `view()`-Operation, um zu prüfen, ob der Dateiinhalt korrekt in den Kanal geladen wird.

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

Beachte, dass wir die Option `header: true` verwenden, um Nextflow anzuweisen, die erste Zeile der CSV-Datei als Kopfzeile zu lesen.

Schauen wir mal, was dabei herauskommt!
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

Wir sehen, dass der Operator für jede Zeile in der CSV-Datei eine Map aus Schlüssel-Wert-Paaren erstellt hat, wobei die Spaltenüberschriften als Schlüssel für die entsprechenden Werte dienen.

Jeder Map-Eintrag entspricht einer Spalte in unserem Datenblatt:

- `id`
- `character`
- `recording`

Großartig! Das macht es einfach, auf bestimmte Felder jeder Datei zuzugreifen.
Zum Beispiel könnten wir mit `id` auf die Datei-ID oder mit `recording` auf den Pfad der txt-Datei zugreifen.

??? info "(Optional) Mehr über Maps"

    In Groovy, der Programmiersprache, auf der Nextflow aufbaut, ist eine Map eine Schlüssel-Wert-Datenstruktur, ähnlich wie Dictionaries in Python, Objekte in JavaScript oder Hashes in Ruby.

    Hier ist ein ausführbares Skript, das zeigt, wie du eine Map definieren und auf ihren Inhalt zugreifen kannst:

    ```groovy title="examples/map_demo.nf"
    #!/usr/bin/env nextflow

    // Eine einfache Map erstellen
    def my_map = [id:'sampleA', character:'squirrel']

    // Die gesamte Map ausgeben
    println "map: ${my_map}"

    // Auf einzelne Werte mit Punktnotation zugreifen
    println "id: ${my_map.id}"
    println "character: ${my_map.character}"
    ```

    Obwohl es keinen richtigen `workflow`-Block hat, kann Nextflow dies wie einen Workflow ausführen:

    ```bash
    nextflow run examples/map_demo.nf
    ```

    Und das kannst du in der Ausgabe erwarten:

    ```console title="Output"
     N E X T F L O W   ~  version 25.10.2

    Launching `map_demo.nf` [cheesy_plateau] DSL2 - revision: fae5b8496e

    map: [id:sampleA, character:squirrel]
    id: sampleA
    character: squirrel
    ```

### 1.2. Bestimmte Felder mit `map` auswählen

Angenommen, wir möchten die Spalte `character` aus dem Datenblatt abrufen und ausgeben.
Wir können den Nextflow-Operator `map` verwenden, um über jedes Element in unserem Kanal zu iterieren und gezielt den `character`-Eintrag aus dem Map-Objekt auszuwählen.

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

Führe den Workflow erneut aus:

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

Erfolg! Wir haben die Map-Struktur aus unserem Datenblatt genutzt, um auf die Werte einzelner Spalten für jede Zeile zuzugreifen.

Nachdem wir das Datenblatt erfolgreich eingelesen haben und auf die Daten in jeder Zeile zugreifen können, können wir mit der Implementierung unserer Pipeline-Logik beginnen.

### 1.3. Die Metadaten in einer 'Meta-Map' organisieren

Im aktuellen Zustand des Workflows stehen die Eingabedateien (unter dem Schlüssel `recording`) und die zugehörigen Metadaten (`id`, `character`) gleichberechtigt nebeneinander, als wären sie alle in einer großen Tasche.
Die praktische Konsequenz ist, dass jeder Prozess, der diesen Kanal verarbeitet, mit dieser Struktur im Hinterkopf konfiguriert werden müsste:

```groovy
    input:
    tuple val(id), val(character), file(recording)
```

Das ist in Ordnung, solange sich die Anzahl der Spalten im Datenblatt nicht ändert.
Wenn du jedoch auch nur eine Spalte zum Datenblatt hinzufügst, stimmt die Form des Kanals nicht mehr mit dem überein, was der Prozess erwartet, und der Workflow wird Fehler produzieren.
Außerdem ist der Prozess schwer mit anderen zu teilen, die möglicherweise leicht unterschiedliche Eingabedaten haben, und du könntest am Ende Variablen fest in den Prozess einprogrammieren müssen, die vom Skript-Block gar nicht benötigt werden.

Um dieses Problem zu vermeiden, müssen wir einen Weg finden, die Kanalstruktur konsistent zu halten, unabhängig davon, wie viele Spalten das Datenblatt enthält.

Das können wir tun, indem wir alle Metadaten in einem Element innerhalb des Tupels zusammenfassen, das wir die Metadaten-Map oder kurz 'Meta-Map' nennen.

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

Wir haben unsere Kanalelemente in ein Tupel aus zwei Elementen umstrukturiert: die Meta-Map und das entsprechende Dateiobjekt.

Führen wir den Workflow aus:

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

Jetzt enthält jedes Element im Kanal zuerst die Meta-Map und dann das entsprechende Dateiobjekt:

```console title="Example output structure"
[
  [id:sampleA, character:squirrel],
  /workspaces/training/side-quests/metadata/data/bonjour.txt
]
```

Das Hinzufügen weiterer Spalten im Datenblatt macht mehr Metadaten in der `meta`-Map verfügbar, ändert aber nicht die Kanalform.
Das ermöglicht es uns, Prozesse zu schreiben, die den Kanal verarbeiten, ohne die Metadaten-Elemente fest in die Eingabespezifikation einprogrammieren zu müssen:

```groovy title="Syntax example"
    input:
    tuple val(meta), file(recording)
```

Dies ist eine weit verbreitete Konvention zur Organisation von Metadaten in Nextflow-Workflows.

### Fazit

In diesem Abschnitt hast du gelernt:

- **Warum Metadaten wichtig sind:** Metadaten bei den Daten zu behalten bewahrt wichtige Dateiinformationen durch den gesamten Workflow.
- **Wie man Datenblätter einliest:** Mit `splitCsv` CSV-Dateien mit Kopfzeileninformationen lesen und Zeilen in strukturierte Daten umwandeln.
- **Wie man eine Meta-Map erstellt:** Metadaten von Dateidaten mit der Tupel-Struktur `[ [id:value, ...], file ]` trennen.

---

## 2. Metadaten bearbeiten

Jetzt, wo wir unsere Metadaten geladen haben, machen wir etwas damit!

Wir werden ein Tool namens [`langid`](https://github.com/saffsd/langid.py) verwenden, um die Sprache in der Aufnahmedatei jeder Figur zu identifizieren.
Das Tool ist auf einer Reihe von Sprachen vortrainiert und gibt bei einem Textausschnitt eine Sprachvorhersage und einen zugehörigen Wahrscheinlichkeitswert aus, beides nach `stdout`.

### 2.1. Den Prozess importieren und den Code untersuchen

Wir stellen dir ein vorgefertigtes Prozessmodul namens `IDENTIFY_LANGUAGE` zur Verfügung, das das `langid`-Tool kapselt. Du musst lediglich eine include-Anweisung vor dem Workflow-Block hinzufügen.

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

Du kannst die Moduldatei öffnen, um den Code zu untersuchen:

```groovy title="modules/langid.nf" linenums="1" hl_lines="9 12"
#!/usr/bin/env nextflow

// langid verwenden, um die Sprache jeder Eingabedatei vorherzusagen
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

Wie du siehst, verwendet die Eingabedefinition dieselbe `tuple val(meta), path(file)`-Struktur, die wir gerade auf unseren Eingabekanal angewendet haben.

Die Ausgabedefinition ist als Tupel mit einer ähnlichen Struktur wie die Eingabe aufgebaut, enthält aber zusätzlich `stdout` als drittes Element.
Dieses `tuple val(meta), path(file), <output>`-Muster hält die Metadaten sowohl mit den Eingabedaten als auch mit den Ausgaben verknüpft, während sie durch die Pipeline fließen.

Beachte, dass wir hier den [`stdout`](https://www.nextflow.io/docs/latest/process.html#outputs)-Ausgabe-Qualifier von Nextflow verwenden, weil das Tool seine Ausgabe direkt in die Konsole schreibt statt in eine Datei. Wir verwenden `sed` in der Befehlszeile, um den Wahrscheinlichkeitswert zu entfernen, den String durch Entfernen von Zeilenumbruchzeichen zu bereinigen und nur die Sprachvorhersage zurückzugeben.

### 2.2. Einen Aufruf von `IDENTIFY_LANGUAGE` hinzufügen

Jetzt, wo der Prozess im Workflow verfügbar ist, können wir einen Aufruf des `IDENTIFY_LANGUAGE`-Prozesses hinzufügen, um ihn auf dem Datenkanal auszuführen.

Nimm die folgenden Änderungen am Workflow vor:

=== "Danach"

    ```groovy title="main.nf" linenums="7" hl_lines="7-9"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [[id: row.id, character: row.character], row.recording]
            }

        // langid ausführen, um die Sprache jeder Begrüßung zu identifizieren
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

Wir können den Workflow jetzt ausführen:

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

Ausgezeichnet! Wir haben jetzt eine Vorhersage, welche Sprache jede Figur spricht.

Und wie bereits erwähnt, haben wir auch die Eingabedatei und die Meta-Map in der Ausgabe eingeschlossen, was bedeutet, dass beide mit den neu produzierten Informationen verknüpft bleiben.
Das wird im nächsten Schritt nützlich sein.

!!! note "Hinweis"

    Allgemeiner gesagt macht dieses Muster, die Meta-Map mit den Ergebnissen verknüpft zu halten, es einfacher, verwandte Ergebnisse zu verknüpfen, die dieselben Bezeichner teilen.

    Wie du bereits gelernt hast, kann man sich nicht auf die Reihenfolge der Elemente in Kanälen verlassen, um Ergebnisse kanalübergreifend zuzuordnen.
    Stattdessen musst du Schlüssel verwenden, um Daten korrekt zuzuordnen, und Meta-Maps bieten dafür eine ideale Struktur.

    Wir erkunden diesen Anwendungsfall ausführlich in der Side Quest [Splitting & Grouping](../splitting_and_grouping/).

### 2.3. Metadaten mit Prozessausgaben erweitern

Da die Ergebnisse, die wir gerade produziert haben, selbst eine Form von Metadaten über den Inhalt der Dateien sind, wäre es nützlich, sie zu unserer Meta-Map hinzuzufügen.

Wir wollen die bestehende Meta-Map jedoch nicht direkt verändern.
Technisch gesehen ist das zwar _möglich_, aber unsicher.

Stattdessen erstellen wir eine neue Meta-Map, die den Inhalt der bestehenden Meta-Map plus ein neues `lang: lang_id`-Schlüssel-Wert-Paar mit den neuen Informationen enthält, unter Verwendung des `+`-Operators (eine Groovy-Funktion).
Und wir kombinieren das mit einer [`map`](https://www.nextflow.io/docs/latest/operator.html#map)-Operation, um die alte Map durch die neue zu ersetzen.

Hier sind die Änderungen, die du am Workflow vornehmen musst:

=== "Danach"

    ```groovy title="main.nf" linenums="13" hl_lines="3-7"
        // langid ausführen, um die Sprache jeder Begrüßung zu identifizieren
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out
            .map { meta, file, lang_id ->
                [meta + [lang: lang_id], file]
            }
            .view()
    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="13" hl_lines="3"
        // langid ausführen, um die Sprache jeder Begrüßung zu identifizieren
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out.view()
    ```

Falls du mit dem `+`-Operator noch nicht vertraut bist oder das verwirrend erscheint, nimm dir ein paar Minuten, um die ausführliche Erklärung unten durchzulesen.

??? info "Erstellung der neuen Meta-Map mit dem `+`-Operator"

    **Zunächst musst du wissen, dass wir den Inhalt zweier Maps mit dem Groovy-Operator `+` zusammenführen können.**

    Angenommen, wir haben die folgenden Maps:

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

    Toll!

    **Aber was, wenn du ein Feld hinzufügen musst, das noch nicht Teil einer Map ist?**

    Angenommen, du beginnst wieder mit `map1`, aber die Sprachvorhersage ist nicht in ihrer eigenen Map (es gibt kein `map2`).
    Stattdessen ist sie in einer Variable namens `lang_id` gespeichert, und du weißt, dass du ihren Wert (`'fr'`) mit dem Schlüssel `lang` speichern möchtest.

    Du kannst tatsächlich Folgendes tun:

    ```groovy
    new_map = [map1 + [lang: lang_id]]
    ```

    Hier erstellt `[lang: new_info]` eine neue namenlose Map auf der Stelle, und `map1 + ` führt `map1` mit der neuen namenlosen Map zusammen, was denselben `new_map`-Inhalt wie zuvor erzeugt.

    Praktisch, oder?

    **Jetzt übertragen wir das in den Kontext einer Nextflow `channel.map()`-Operation.**

    Der Code wird zu:

    ```groovy
    .map { map1, lang_id ->
        [map1 + [lang: lang_id]]
    }
    ```

    Das macht Folgendes:

    - `map1, lang_id ->` nimmt die zwei Elemente im Tupel
    - `[map1 + [lang: lang_id]]` erstellt die neue Map wie oben beschrieben

    Die Ausgabe ist eine einzelne namenlose Map mit demselben Inhalt wie `new_map` in unserem Beispiel oben.
    Wir haben also effektiv transformiert:

    ```groovy
    [id: 'sampleA', character: 'squirrel'], 'fr'
    ```

    in:

    ```groovy
    [id: 'sampleA', character: 'squirrel', lang: 'fr']
    ```

    Hoffentlich siehst du, dass wenn wir `map1` durch `meta` ersetzen, das im Wesentlichen alles ist, was wir brauchen, um die Sprachvorhersage zu unserer Meta-Map in unserem Workflow hinzuzufügen.

    Bis auf eine Sache!

    In unserem Workflow **müssen wir auch das `file`-Objekt im Tupel berücksichtigen**, das aus `meta, file, lang_id` besteht.

    Der Code wird also zu:

    ```groovy
    .map { meta, file, lang_id ->
        [meta + [lang: lang_id], file]
    }
    ```

    Falls es schwer nachzuvollziehen ist, warum sich `file` in der `map`-Operation scheinbar bewegt, stell dir vor, dass statt `[meta + [lang: lang_id], file]` diese Zeile `[new_map, file]` lautet.
    Das sollte klarer machen, dass wir `file` einfach an seiner ursprünglichen zweiten Position im Tupel belassen. Wir haben nur den `new_info`-Wert in die Map an erster Position eingefaltet.

    **Und das bringt uns zurück zur `tuple val(meta), path(file)`-Kanalstruktur!**

Sobald du sicher bist, dass du verstehst, was dieser Code tut, führe den Workflow aus, um zu sehen, ob es funktioniert hat:

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
Wir haben die Ausgabe des Prozesses sauber von `meta, file, lang_id` umstrukturiert, sodass `lang_id` jetzt einer der Schlüssel in der Meta-Map ist und die Tupel des Kanals wieder dem `meta, file`-Modell entsprechen.

<!-- TODO (future) Should we also show how to remove a key using subMap?! Or note where to find that. -->

### 2.4. Eine Sprachgruppe mit Bedingungen zuweisen

Jetzt, wo wir unsere Sprachvorhersagen haben, nutzen wir die Informationen, um neue Gruppierungen zuzuweisen.

In unseren Beispieldaten können die von unseren Figuren verwendeten Sprachen in germanische Sprachen (Englisch, Deutsch) und romanische Sprachen (Französisch, Spanisch, Italienisch) eingeteilt werden.
Es könnte nützlich sein, diese Klassifizierung später in der Pipeline direkt verfügbar zu haben, also fügen wir diese Information zur Meta-Map hinzu.

Und gute Neuigkeiten: Das ist wieder ein Fall, der sich perfekt für den `map`-Operator eignet!

Konkret werden wir eine Variable namens `lang_group` definieren und einfache bedingte Logik verwenden, um zu bestimmen, welchen Wert `lang_group` für jedes Datenelement erhalten soll.

Die allgemeine Syntax wird so aussehen:

```groovy
.map { meta, file ->

    // Bedingte Logik zur Definition von lang_group kommt hier hin

    [meta + [lang_group: lang_group], file]
}
```

Du siehst, das ist sehr ähnlich wie die On-the-fly-Map-Zusammenführungsoperation, die wir im vorherigen Schritt verwendet haben.
Wir müssen nur die bedingten Anweisungen schreiben.

Hier ist die bedingte Logik, die wir anwenden möchten:

- Definiere eine Variable namens `lang_group` mit dem Standardwert `'unknown'`.
- Wenn `lang` entweder Deutsch (`'de'`) oder Englisch (`'en'`) ist, ändere `lang_group` zu `germanic`.
- Sonst wenn `lang` in einer Liste mit Französisch (`'fr'`), Spanisch (`'es'`) und Italienisch (`'it'`) enthalten ist, ändere `lang_group` zu `romance`.

Versuche es selbst zu schreiben, wenn du bereits weißt, wie man bedingte Anweisungen in Nextflow schreibt.

!!! tip "Tipp"

    Du kannst auf den Wert von `lang` innerhalb der map-Operation mit `meta.lang` zugreifen.

Du solltest am Ende die folgenden Änderungen am Workflow vornehmen:

=== "Danach"

    ```groovy title="main.nf" linenums="13" hl_lines="7-19"
        // langid ausführen, um die Sprache jeder Begrüßung zu identifizieren
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
        // langid ausführen, um die Sprache jeder Begrüßung zu identifizieren
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out
            .map { meta, file, lang_id ->
                [meta + [lang: lang_id], file]
            }
            .view()
    ```

Hier sind die wichtigsten Punkte:

- Wir verwenden `def lang_group = "unknown"`, um die Variable `lang_group` mit dem Standardwert `unknown` zu erstellen.
- Wir verwenden eine `if {} else if {}`-Struktur für die bedingte Logik, mit alternativen `.equals()`-Tests für die zwei germanischen Sprachen und einem Test auf Vorhandensein in einer Liste für die drei romanischen Sprachen.
- Wir verwenden die `meta + [lang_group:lang_group]`-Zusammenführungsoperation wie zuvor, um die aktualisierte Meta-Map zu erzeugen.

<!-- TODO (future) Add note/links to relevant docs in additional resources section -->

Sobald das alles Sinn ergibt, führe den Workflow erneut aus, um das Ergebnis zu sehen:

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

Wie du siehst, behalten die Kanalelemente ihre `[meta, file]`-Struktur, aber die Meta-Map enthält jetzt diese neue Klassifizierung.

### Fazit

In diesem Abschnitt hast du gelernt:

- **Eingabe-Metadaten auf Ausgabekanäle anwenden**: Das Kopieren von Metadaten auf diese Weise ermöglicht es uns, Ergebnisse später anhand von Metadateninhalten zu verknüpfen.
- **Benutzerdefinierte Schlüssel erstellen**: Du hast zwei neue Schlüssel in deiner Meta-Map erstellt und sie mit `meta + [new_key:value]` in die bestehende Meta-Map eingeführt. Einen basierend auf einem berechneten Wert aus einem Prozess und einen basierend auf einer Bedingung, die du im `map`-Operator gesetzt hast.

Das ermöglicht es dir, neue und bestehende Metadaten mit Dateien zu verknüpfen, während du durch deine Pipeline fortschreitest.
Auch wenn du Metadaten nicht als Teil eines Prozesses verwendest, macht es das Beibehalten der Meta-Map bei den Daten einfach, alle relevanten Informationen zusammenzuhalten.

---

## 3. Meta-Map-Informationen in einem Prozess verwenden

Jetzt, wo du weißt, wie man die Meta-Map erstellt und aktualisiert, kommen wir zum wirklich spaßigen Teil: die Metadaten tatsächlich in einem Prozess zu verwenden.

Konkret werden wir einen zweiten Schritt zu unserem Workflow hinzufügen, um jedes Tier als ASCII-Kunst zu zeichnen und es den aufgezeichneten Text in einer Sprechblase sagen zu lassen.
Dazu verwenden wir ein Tool namens [`cowpy`](https://github.com/jeffbuttars/cowpy).

??? info "Was macht `cowpy`?"

    `cowpy` ist ein Befehlszeilen-Tool, das ASCII-Kunst generiert, um beliebige Texteingaben auf unterhaltsame Weise darzustellen.
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

    Optional kannst du eine Figur (oder 'cowacter') auswählen, die anstelle der Standard-Kuh verwendet wird.

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
Falls nicht, keine Sorge; wir werden alles erklären, was du wissen musst.

### 3.1. Den Prozess importieren und den Code untersuchen

Wir stellen dir ein vorgefertigtes Prozessmodul namens `COWPY` zur Verfügung, das das `cowpy`-Tool kapselt. Du musst lediglich eine include-Anweisung vor dem Workflow-Block hinzufügen.

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

Du kannst die Moduldatei öffnen, um den Code zu untersuchen:

```groovy title="modules/cowpy.nf" linenums="1"
#!/usr/bin/env nextflow

// ASCII-Kunst mit cowpy generieren
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

Wie du siehst, ist dieser Prozess derzeit so konzipiert, dass er eine Eingabedatei (mit dem anzuzeigenden Text) und einen Wert entgegennimmt, der die Figur angibt, die als ASCII-Kunst gezeichnet werden soll – normalerweise auf Workflow-Ebene durch einen Befehlszeilenparameter bereitgestellt.

### 3.2. Ein Meta-Map-Feld als Eingabe übergeben

Als wir das `cowpy`-Tool im Hello Nextflow-Kurs verwendet haben, haben wir einen Befehlszeilenparameter verwendet, um zu bestimmen, welche Figur für das endgültige Bild verwendet werden soll.
Das war sinnvoll, weil wir pro Pipeline-Ausführung nur ein Bild generiert haben.

In diesem Tutorial möchten wir jedoch für jedes Subjekt, das wir verarbeiten, ein passendes Bild generieren, sodass ein Befehlszeilenparameter zu einschränkend wäre.

Gute Neuigkeiten: Wir haben eine `character`-Spalte in unserem Datenblatt und damit in unserer Meta-Map.
Nutzen wir das, um die Figur festzulegen, die der Prozess für jeden Eintrag verwenden soll.

Dazu müssen wir drei Dinge tun:

1. Dem Ausgabekanal des vorherigen Prozesses einen Namen geben, damit wir bequemer damit arbeiten können.
2. Bestimmen, wie wir auf die gewünschten Informationen zugreifen.
3. Einen Aufruf des zweiten Prozesses hinzufügen und die Informationen entsprechend einspeisen.

Fangen wir an.

#### 3.2.1. Den vorherigen Ausgabekanal benennen

Wir haben die vorherigen Manipulationen direkt auf dem Ausgabekanal des ersten Prozesses, `IDENTIFY_LANGUAGE.out`, angewendet.
Um den Inhalt dieses Kanals an den nächsten Prozess weiterzugeben (und das auf eine klare und leicht lesbare Weise), möchten wir ihm einen eigenen Namen geben: `ch_languages`.

Das können wir mit dem [`set`](https://www.nextflow.io/docs/latest/reference/operator.html#set)-Operator tun.

Ersetze im Haupt-Workflow den `.view()`-Operator durch `.set { ch_languages }` und füge eine Zeile hinzu, die testet, ob wir den Kanal beim Namen nennen können.

=== "Danach"

    ```groovy title="main.nf" linenums="14" hl_lines="19 21 22"
        // langid ausführen, um die Sprache jeder Begrüßung zu identifizieren
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

        // Temporär: einen Blick in ch_languages werfen
        ch_languages.view()
    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="14" hl_lines="19"
        // langid ausführen, um die Sprache jeder Begrüßung zu identifizieren
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

Führen wir das aus:

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

Das bestätigt, dass wir den Kanal jetzt beim Namen nennen können.

#### 3.2.2. Auf die Datei- und Charakter-Metadaten zugreifen

Aus dem Modulcode wissen wir, dass der `COWPY`-Prozess eine Textdatei und einen `character`-Wert erwartet.
Um den Aufruf des `COWPY`-Prozesses zu schreiben, müssen wir nur wissen, wie wir das entsprechende Dateiobjekt und die Metadaten aus jedem Element im Kanal extrahieren.

Wie so oft ist der einfachste Weg dafür eine `map`-Operation.

Unser Kanal enthält Tupel mit der Struktur `[meta, file]`, sodass wir direkt auf das `file`-Objekt zugreifen können und auf den in der Meta-Map gespeicherten `character`-Wert mit `meta.character` zugreifen können.

Nimm im Haupt-Workflow die folgenden Codeänderungen vor:

=== "Danach"

    ```groovy title="main.nf" linenums="34"
        // Temporär: auf Datei und Charakter zugreifen
        ch_languages.map { meta, file -> file }.view { file -> "File: " + file }
        ch_languages.map { meta, file -> meta.character }.view { character -> "Character: " + character }
    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="34"
        // Temporär: einen Blick in ch_languages werfen
        ch_languages.view()
    ```

Beachte, dass wir closures (wie `{ file -> "File: " + file }`) verwenden, um die Ausgabe der `.view`-Operationen lesbarer zu machen.

Führen wir das aus:

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

Das bestätigt, dass wir auf die Datei und den Charakter für jedes Element im Kanal zugreifen können.

#### 3.2.3. Den `COWPY`-Prozess aufrufen

Jetzt fügen wir alles zusammen und rufen den `COWPY`-Prozess tatsächlich auf dem `ch_languages`-Kanal auf.

Nimm im Haupt-Workflow die folgenden Codeänderungen vor:

=== "Danach"

    ```groovy title="main.nf" linenums="34"
        // cowpy ausführen, um ASCII-Kunst zu generieren
        COWPY(
            ch_languages.map { meta, file -> file },
            ch_languages.map { meta, file -> meta.character }
        )
    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="34"
        // Temporär: auf Datei und Charakter zugreifen
        ch_languages.map { meta, file -> [file, meta.character] }
            .view()
    ```

Du siehst, dass wir einfach die zwei map-Operationen (ohne die `.view()`-Anweisungen) als Eingaben für den Prozessaufruf kopieren.
Vergiss nur nicht das Komma zwischen ihnen!

Das ist etwas umständlich, aber wir werden im nächsten Abschnitt sehen, wie wir das verbessern können.

Führen wir das aus:

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

Wenn du im Ergebnisverzeichnis nachschaust, solltest du die einzelnen Dateien mit der ASCII-Kunst jeder Begrüßung sehen, gesprochen von der entsprechenden Figur.

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

Das zeigt, dass wir die Informationen in der Meta-Map verwenden konnten, um den Befehl im zweiten Schritt der Pipeline zu parametrisieren.

Wie oben angemerkt, war ein Teil des beteiligten Codes etwas umständlich, da wir Metadaten noch im Kontext des Workflow-Körpers entpacken mussten.
Dieser Ansatz funktioniert gut für die Verwendung einer kleinen Anzahl von Feldern aus der Meta-Map, würde aber schlecht skalieren, wenn wir viel mehr verwenden wollten.

Es gibt einen anderen Operator namens `multiMap()`, der uns das etwas vereinfachen lässt, aber auch das ist nicht ideal.

??? info "(Optional) Alternative Version mit `multiMap()`"

    Falls du dich fragst: Wir konnten keine einzelne `map()`-Operation schreiben, die sowohl `file` als auch `character` ausgibt, weil das sie als Tupel zurückgeben würde.
    Wir mussten zwei separate `map()`-Operationen schreiben, um die `file`- und `character`-Elemente separat an den Prozess zu übergeben.

    Technisch gibt es eine andere Möglichkeit, das durch eine einzelne Mapping-Operation zu tun, mit dem `multiMap()`-Operator, der mehrere Kanäle ausgeben kann.
    Du könntest zum Beispiel den Aufruf von `COWPY` oben durch folgenden Code ersetzen:

    === "Danach"

        ```groovy title="main.nf" linenums="34"
            // cowpy ausführen, um ASCII-Kunst zu generieren
            COWPY(
                ch_languages.multiMap { meta, file ->
                    file: file
                    character: meta.character
                }
            )
        ```

    === "Vorher"

        ```groovy title="main.nf" linenums="34"
            // cowpy ausführen, um ASCII-Kunst zu generieren
            COWPY(
                ch_languages.map { meta, file -> file },
                ch_languages.map { meta, file -> meta.character }
            )
        ```

    Das erzeugt genau dasselbe Ergebnis.

In beiden Fällen ist es umständlich, dass wir auf Workflow-Ebene etwas entpacken müssen.

Es wäre besser, wenn wir die gesamte Meta-Map an den Prozess übergeben und dort auswählen könnten, was wir brauchen.

### 3.3. Die gesamte Meta-Map übergeben und verwenden

Der Sinn der Meta-Map ist es schließlich, alle Metadaten zusammen als Paket weiterzugeben.
Der einzige Grund, warum wir das oben nicht tun konnten, ist, dass der Prozess nicht so eingerichtet ist, eine Meta-Map zu akzeptieren.
Da wir aber den Prozesscode kontrollieren, können wir das ändern.

Lass uns den `COWPY`-Prozess so modifizieren, dass er die `[meta, file]`-Tupel-Struktur akzeptiert, die wir im ersten Prozess verwendet haben, damit wir den Workflow vereinfachen können.

Dazu müssen wir drei Dinge tun:

1. Die Eingabedefinitionen des `COWPY`-Prozessmoduls ändern
2. Den Prozessbefehl aktualisieren, um die Meta-Map zu verwenden
3. Den Prozessaufruf im Workflow-Körper aktualisieren

Bereit? Los geht's!

#### 3.3.1. Die `COWPY`-Moduleingabe ändern

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

Das ermöglicht uns, die `[meta, file]`-Tupel-Struktur zu verwenden, die wir früher im Tutorial behandelt haben.

Beachte, dass wir die Prozessausgabedefinition nicht aktualisiert haben, um die Meta-Map auszugeben, um das Tutorial übersichtlich zu halten. Du kannst das aber gerne selbst als Übung nach dem Vorbild des `IDENTIFY_LANGUAGE`-Prozesses tun.

#### 3.3.2. Den Befehl aktualisieren, um das Meta-Map-Feld zu verwenden

Die gesamte Meta-Map ist jetzt innerhalb des Prozesses verfügbar, sodass wir direkt aus dem Befehlsblock auf die darin enthaltenen Informationen zugreifen können.

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

Wir haben den Verweis auf den `character`-Wert, der zuvor als eigenständige Eingabe übergeben wurde, durch den in der Meta-Map gespeicherten Wert ersetzt, auf den wir mit `meta.character` zugreifen.

Jetzt aktualisieren wir den Prozessaufruf entsprechend.

#### 3.3.3. Den Prozessaufruf aktualisieren und ausführen

Der Prozess erwartet jetzt, dass seine Eingabe die `[meta, file]`-Tupel-Struktur verwendet, was der vorherige Prozess ausgibt. Wir können also einfach den gesamten `ch_languages`-Kanal an den `COWPY`-Prozess übergeben.

Nimm die folgenden Änderungen am Haupt-Workflow vor:

=== "Danach"

    ```groovy title="main.nf" linenums="34" hl_lines="2"
    // cowpy ausführen, um ASCII-Kunst zu generieren
    COWPY(ch_languages)
    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="34" hl_lines="3-4"
    // cowpy ausführen, um ASCII-Kunst zu generieren
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

Wenn du im Ergebnisverzeichnis nachschaust, solltest du dieselben Ausgaben wie zuvor sehen, d.h. einzelne Dateien mit der ASCII-Kunst jeder Begrüßung, gesprochen von der entsprechenden Figur.

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

Das erzeugt also dieselben Ergebnisse wie zuvor, aber mit einfacherem Code.

Natürlich setzt das voraus, dass du den Prozesscode ändern kannst.
In manchen Fällen musst du dich auf bestehende Prozesse verlassen, die du nicht ändern darfst, was deine Möglichkeiten einschränkt.
Die gute Nachricht, wenn du planst, Module aus dem [nf-core](https://nf-co.re/)-Projekt zu verwenden, ist, dass nf-core-Module alle standardmäßig die `[meta, file]`-Tupel-Struktur verwenden.

### 3.4. Fehlende erforderliche Eingaben beheben

Der `character`-Wert ist erforderlich, damit der `COWPY`-Prozess erfolgreich ausgeführt werden kann.
Wenn wir keinen Standardwert in einer Konfigurationsdatei festlegen, MÜSSEN wir einen Wert dafür im Datenblatt angeben.

**Was passiert, wenn wir das nicht tun?**
Das hängt davon ab, was das Eingabe-Datenblatt enthält und welche Version des Workflows wir ausführen.

#### 3.4.1. Die Charakter-Spalte existiert, ist aber leer

Angenommen, wir löschen den Charakterwert für einen der Einträge in unserem Datenblatt, um einen Datenerfassungsfehler zu simulieren:

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

Bei beiden Versionen des Workflows, die wir oben verwendet haben, wird der `character`-Schlüssel für alle Einträge erstellt, wenn das Datenblatt eingelesen wird, aber für `sampleA` ist der Wert ein leerer String.

Das wird einen Fehler verursachen.

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

Wenn Nextflow die `cowpy`-Befehlszeile für diese Probe ausführt, wird `${meta.character}` mit einem leeren String in der `cowpy`-Befehlszeile gefüllt, sodass das `cowpy`-Tool einen Fehler ausgibt, der besagt, dass kein Wert für das `-c`-Argument angegeben wurde.

#### 3.4.2. Die Charakter-Spalte existiert nicht im Datenblatt

Angenommen, wir löschen die `character`-Spalte vollständig aus unserem Datenblatt:

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

Wenn wir die Version des Codes aus Abschnitt 3.2 verwenden, versucht Nextflow, auf den `character`-Schlüssel in der Meta-Map zuzugreifen, BEVOR der `COWPY`-Prozess aufgerufen wird.

Es findet keine Elemente, die der Anweisung entsprechen, sodass `COWPY` überhaupt nicht ausgeführt wird.

??? success "Befehlsausgabe"

    ```console hl_lines="7"
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [desperate_montalcini] DSL2 - revision: 0dfeee3cc1

    executor >  local (7)
    [1a/df2544] process > IDENTIFY_LANGUAGE (7) [100%] 7 of 7 ✔
    [-        ] process > COWPY                 -
    ```

Aus Nextflows Sicht wurde dieser Workflow erfolgreich ausgeführt!
Allerdings werden keine der gewünschten Ausgaben produziert.

##### 3.4.2.2. Wert wird auf Prozessebene abgerufen

Wenn wir die Version aus Abschnitt 3.3 verwenden, übergibt Nextflow die gesamte Meta-Map an den `COWPY`-Prozess und versucht, den Befehl auszuführen.

Das verursacht einen Fehler, aber einen anderen als im ersten Fall.

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

Das passiert, weil `meta.character` nicht existiert, sodass unser Versuch, darauf zuzugreifen, `null` zurückgibt. Als Ergebnis setzt Nextflow buchstäblich `null` in die Befehlszeile ein, was vom `cowpy`-Tool natürlich nicht erkannt wird.

#### 3.4.3. Lösungen

Abgesehen davon, einen Standardwert als Teil der Workflow-Konfiguration bereitzustellen, gibt es zwei Dinge, die wir tun können, um das robuster zu handhaben:

1. Implementiere eine Eingabevalidierung in deinem Workflow, um sicherzustellen, dass das Datenblatt alle erforderlichen Informationen enthält. Eine [Einführung in die Eingabevalidierung](../hello_nf-core/05_input_validation.md) findest du im Hello nf-core-Trainingskurs. <!-- TODO (future) pending a proper Validation side quest -->

2. Wenn du sicherstellen möchtest, dass jede Person, die dein Prozessmodul verwendet, erforderliche Eingaben sofort erkennen kann, kannst du die erforderliche Metadateneigenschaft auch als explizite Eingabe definieren.

Hier ist ein Beispiel, wie das funktionieren würde.

Aktualisiere zunächst auf Prozessebene die Eingabedefinition wie folgt:

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

Verwende dann auf Workflow-Ebene eine Mapping-Operation, um die `character`-Eigenschaft aus den Metadaten zu extrahieren und sie zu einer expliziten Komponente des Eingabe-Tupels zu machen:

=== "Danach"

    ```groovy title="main.nf" linenums="37" hl_lines="1"
        COWPY(ch_languages.map{meta, file -> [meta, meta.character, file]})
    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="37" hl_lines="1"
        COWPY(ch_languages)
    ```

Dieser Ansatz hat den Vorteil, dass er explizit zeigt, dass `character` erforderlich ist, und macht den Prozess einfacher in anderen Kontexten wiederzuverwenden.

Das verdeutlicht ein wichtiges Designprinzip:

**Verwende die Meta-Map für optionale, beschreibende Informationen, aber extrahiere erforderliche Werte als explizite Eingaben.**

Die Meta-Map eignet sich hervorragend, um Kanalstrukturen sauber zu halten und beliebige Kanalstrukturen zu vermeiden. Für obligatorische Elemente, die direkt in einem Prozess referenziert werden, schafft das Extrahieren als explizite Eingaben jedoch robusteren und wartbareren Code.

### Fazit

In diesem Abschnitt hast du gelernt, wie du Metadaten nutzt, um die Ausführung eines Prozesses anzupassen – entweder auf Workflow-Ebene oder auf Prozessebene.

---

## Ergänzende Übung

Wenn du das Verwenden von Meta-Map-Informationen innerhalb eines Prozesses üben möchtest, versuche andere Informationen aus der Meta-Map wie `lang` und `lang_group` zu verwenden, um anzupassen, wie die Ausgaben benannt und/oder organisiert werden.

Versuche zum Beispiel, den Code so zu ändern, dass dieses Ergebnis erzeugt wird:

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

In dieser Side Quest hast du erkundet, wie man effektiv mit Metadaten in Nextflow-Workflows arbeitet.

Dieses Muster, Metadaten explizit und mit den Daten verknüpft zu halten, ist eine zentrale Best Practice in Nextflow und bietet mehrere Vorteile gegenüber dem Fest-Einprogrammieren von Dateiinformationen:

- Datei-Metadaten bleiben während des gesamten Workflows mit den Dateien verknüpft
- Das Prozessverhalten kann pro Datei angepasst werden
- Die Ausgabeorganisation kann Datei-Metadaten widerspiegeln
- Dateiinformationen können während der Pipeline-Ausführung erweitert werden

Die Anwendung dieses Musters in deiner eigenen Arbeit ermöglicht es dir, robuste, wartbare bioinformatische Workflows zu entwickeln.

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

### Weitere Ressourcen

- [map](https://www.nextflow.io/docs/latest/operator.html#map)
- [stdout](https://www.nextflow.io/docs/latest/process.html#outputs)

---

## Wie geht es weiter?

Kehre zum [Menü der Side Quests](../) zurück oder klicke auf die Schaltfläche unten rechts auf der Seite, um zum nächsten Thema in der Liste zu wechseln.
