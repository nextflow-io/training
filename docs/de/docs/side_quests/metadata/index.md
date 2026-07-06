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
- Zu verstehen, warum die „Meta-Map + Datendatei"-Schnittstelle eine weit verbreitete Konvention ist
- Neue Metadatenfelder während der Workflow-Ausführung hinzuzufügen
- Metadaten zur Anpassung des Prozessverhaltens und zur Organisation von Ausgaben zu nutzen

Diese Fähigkeiten helfen dir, robustere und flexiblere Pipelines zu entwickeln, die komplexe Dateibeziehungen und Verarbeitungsanforderungen bewältigen können.

### Voraussetzungen

Bevor du diese Side Quest angehst, solltest du:

- Das Tutorial [Hello Nextflow](../../hello_nextflow/index.md) oder einen gleichwertigen Einsteigerkurs abgeschlossen haben.
- Mit grundlegenden Nextflow-Konzepten und -Mechanismen (Prozesse, Kanäle, Operatoren) vertraut sein.

---

## 0. Erste Schritte

#### Öffne den Training-Codespace

Falls du das noch nicht getan hast, öffne die Trainingsumgebung wie in der [Umgebung einrichten](../../envsetup/index.md)-Anleitung beschrieben.

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

Der Editor öffnet sich mit dem Projektverzeichnis im Fokus.

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

Wir werden ein Tool namens [`COWPY`](https://github.com/jeffbuttars/cowpy) verwenden, um ASCII-Kunst jeder Figur zu generieren, die ihre aufgezeichnete Begrüßung spricht.

??? info "Was macht `COWPY`?"

    `COWPY` ist ein Befehlszeilen-Tool, das ASCII-Kunst generiert, um beliebige Texteingaben auf unterhaltsame Weise darzustellen.
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

Außerdem verwenden wir ein Sprachanalyse-Tool namens `langid`, um zu erkennen, welche Sprache jede Figur spricht, und die Ausgaben der Pipeline entsprechend zu organisieren.

#### Schau dir die Aufgabe an

Deine Aufgabe ist es, einen Nextflow-Workflow zu schreiben, der:

1. **ASCII-Kunst** jeder Figur generiert
2. Ausgaben nach Sprachfamilie **organisiert** (Germanisch vs. Romanisch)

Dies ist ein typisches Workflow-Muster, bei dem dateispezifische Metadaten Verarbeitungsentscheidungen steuern – genau die Art von Problem, die Meta-Maps elegant lösen.

#### Bereitschafts-Checkliste

Bereit zum Eintauchen?

- [ ] Ich verstehe das Ziel dieses Kurses und seine Voraussetzungen
- [ ] Mein Codespace läuft
- [ ] Ich habe mein Arbeitsverzeichnis entsprechend eingestellt
- [ ] Ich verstehe die Aufgabe

Wenn du alle Punkte abhaken kannst, kann es losgehen.

---

## 1. Grundlegende Möglichkeiten zum Laden und Verwenden von Metadaten

Öffne die `main.nf`-Workflow-Datei, um das Workflow-Grundgerüst zu untersuchen, das wir dir als Ausgangspunkt geben.

```groovy title="main.nf" linenums="1"
#!/usr/bin/env nextflow

workflow  {
    main:
    ch_datasheet = channel.fromPath("./data/datasheet.csv")
        .splitCsv(header: true)
        .view()

    publish:
    cowpy_art = channel.empty()
}

output {
    cowpy_art {
    }
}
```

Der [`splitCsv`](https://www.nextflow.io/docs/latest/reference/operator.html#splitcsv)-Operator liest jede Zeile der Datei als Element im Kanal.
Das ist derselbe Ansatz, den wir in Hello Nextflow, unserem Einsteigerkurs, zum Laden von CSV-Daten verwenden.
Schau dir [diesen Abschnitt](../../hello_nextflow/02_hello_channels.md#4-read-input-values-from-a-csv-file) an, falls du eine Auffrischung brauchst.

Mit `header: true` wird die erste Zeile als Spaltenüberschriften behandelt, sodass jedes Element eine Map aus Schlüssel-Wert-Paaren wird, die nach Spaltennamen geordnet sind.

Da wir noch keine Prozesse auf den Daten ausführen, sind die `publish`- und `output`-Blöcke nur Platzhalter.

### 1.1. Den Workflow ausführen

Führe den Workflow aus, um zu sehen, wie der Kanalinhalt strukturiert ist, sobald alles geladen ist:

```bash
nextflow run main.nf
```

??? success "Befehlsausgabe"

    ```console
     N E X T F L O W   ~  version 25.10.4

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

Das macht es einfach, auf bestimmte Felder jeder Zeile zuzugreifen.
Zum Beispiel könnten wir mit `id` auf die Datei-ID oder mit `recording` auf den Pfad der txt-Datei zugreifen.

??? info "(Optional) Mehr über Groovy-Maps"

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
     N E X T F L O W   ~  version 25.10.4

    Launching `map_demo.nf` [cheesy_plateau] DSL2 - revision: fae5b8496e

    map: [id:sampleA, character:squirrel]
    id: sampleA
    character: squirrel
    ```

### 1.2. Ein bestimmtes Feld mit `map` auswählen

Wir verwenden den `map`-Operator, um über jedes Element in einem Kanal zu iterieren und gezielt das `character`-Feld auszuwählen, auf das wir per Punktnotation zugreifen können.

#### 1.2.1. Die map-Operation hinzufügen

Um auf die `character`-Spalte zuzugreifen, füge die `map`-Operation vor der `.view()`-Operation wie folgt hinzu:

=== "Danach"

    ```groovy title="main.nf" linenums="5" hl_lines="3-5"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                row.character
            }
            .view()
    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="5"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .view()
    ```

Diese Art, auf ein bestimmtes Feld zuzugreifen, wird in [diesem Abschnitt](../../hello_nextflow/02_hello_channels.md#43-use-the-map-operator-to-extract-the-greetings) von Hello Nextflow ausführlicher erklärt, falls du eine Auffrischung brauchst.

#### 1.2.2. Den Workflow ausführen

Führe den Workflow aus, um zu überprüfen, ob du die extrahierten Charakternamen siehst.

```bash
nextflow run main.nf
```

??? success "Befehlsausgabe"

    ```console
     N E X T F L O W   ~  version 25.10.4

    Launching `main.nf` [exotic_albattani] DSL2 - revision: c0d03cec83

    squirrel
    tux
    sheep
    turkey
    stegosaurus
    moose
    turtle
    ```

Das zeigt, dass wir auf die Werte der `character`-Spalte für jede Zeile zugreifen können.

Jetzt machen wir etwas mit diesen Daten: Wir verwenden die Felder `character` und `recording` zusammen, um mit `COWPY` ASCII-Kunst zu generieren.

### 1.3. Sub-Kanäle mit `multiMap` ausgeben

Wir stellen dir ein vorgefertigtes `COWPY`-Prozessmodul zur Verfügung. Zuerst musst du die Eingabeanforderungen des Prozesses untersuchen.

Du kannst die Datei öffnen, um den Prozess zu sehen:

```groovy title="modules/cowpy.nf" linenums="1" hl_lines="7 8"
// ASCII-Kunst mit cowpy generieren
process COWPY {

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

Wie du siehst, erwartet der Prozess zwei separate Eingaben: eine Aufnahmedatei und einen Charakternamen.
Wir haben Werte für beides, aber sie sind derzeit in jedem Element des Kanals gebündelt.

Eine Möglichkeit, mehrere Felder in separate Kanäle zu extrahieren, ist der [`multiMap`](https://www.nextflow.io/docs/latest/reference/operator.html#multimap)-Operator, der einen Kanal in mehrere benannte Sub-Kanäle in einer einzigen Operation aufteilt.

#### 1.3.1. Die multiMap-Operation hinzufügen

Ersetze die `map`-Operation durch `multiMap`:

=== "Danach"

    ```groovy title="main.nf" linenums="5" hl_lines="3-6"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .multiMap { row ->
                file: row.recording
                character: row.character
            }
    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="5" hl_lines="3-6"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                row.character
            }
            .view()
    ```

Der `multiMap`-Block definiert zwei benannte Sub-Kanäle (`file` und `character`) aus jeder Zeile, auf die wir als `ch_datasheet.file` und `ch_datasheet.character` zugreifen können.

#### 1.3.2. COWPY auf den Sub-Kanälen aufrufen

Füge nun das `COWPY`-Prozessmodul ein und übergib jeden Sub-Kanal als separates Argument:

=== "Danach"

    ```groovy title="main.nf" linenums="1" hl_lines="3 14"
    #!/usr/bin/env nextflow

    include { COWPY } from './modules/cowpy.nf'

    workflow {
        main:
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .multiMap { row ->
                file: row.recording
                character: row.character
            }

        COWPY(ch_datasheet.file, ch_datasheet.character)

        publish:
        cowpy_art = channel.empty()
    }

    output {
        cowpy_art {
        }
    }
    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="1"
    #!/usr/bin/env nextflow

    workflow {
        main:
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .multiMap { row ->
                file: row.recording
                character: row.character
            }

        publish:
        cowpy_art = channel.empty()
    }

    output {
        cowpy_art {
        }
    }
    ```

So können wir die zwei Felder separat übergeben, wie es `COWPY` erfordert.

#### 1.3.3. Die Ausgabeveröffentlichung einrichten

Füge abschließend die Ausgabe von `COWPY` zum `publish:`-Block hinzu:

=== "Danach"

    ```groovy title="main.nf" linenums="14" hl_lines="4"
        COWPY(ch_datasheet.file, ch_datasheet.character)

        publish:
        cowpy_art = COWPY.out
    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="14" hl_lines="4"
        COWPY(ch_datasheet.file, ch_datasheet.character)

        publish:
        cowpy_art = channel.empty()
    ```

So können wir die vom Workflow produzierten Ausgaben leicht einsehen.

#### 1.3.4. Den Workflow ausführen

Führe den Workflow aus, um zu überprüfen, dass `COWPY` mit den bereitgestellten Eingaben läuft:

```bash
nextflow run main.nf
```

??? success "Befehlsausgabe"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [clever_dijkstra] DSL2 - revision: a1b2c3d4e5

    executor >  local (7)
    [3a/f1c290] COWPY (7) [100%] 7 of 7 ✔
    ```

Wie du siehst, hat `COWPY` jede Datei mit der richtigen Figur verarbeitet.

??? abstract "Results-Verzeichnis Inhalt"

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

??? example "Inhalt von results/cowpy-guten_tag.txt"

    ```console
    $ cat results/cowpy-guten_tag.txt
     _____________________________
    / Guten Tag, wie geht es dir? \
    \ Auf Wiedersehen, bis morgen /
     -----------------------------
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

Dieser Ansatz funktioniert, hat aber eine Einschränkung: Wir mussten den Kanal in zwei separate Sub-Kanäle aufteilen.
Wenn wir mehr Felder an den Prozess übergeben wollten, müssten wir sie in weitere Sub-Kanäle aufteilen.
Das kann schnell unübersichtlich werden.

Gute Neuigkeiten: Es gibt einen einfacheren Weg.

### 1.4. Alles als einzelne Eingabe an den Prozess übergeben

Anstatt die Felder in separate Kanäle aufzuteilen, können wir den Prozess so aktualisieren, dass er alle Eingaben als einzelnes Tupel empfängt. Das vereinfacht den Prozessaufruf.

#### 1.4.1. Den COWPY-Prozess aktualisieren

Aktualisiere `COWPY`, um ein Tupel zu akzeptieren, das den drei Elementen jeder Zeile entspricht:

=== "Danach"

    ```groovy title="modules/cowpy.nf" linenums="1" hl_lines="7 10 14"
    // ASCII-Kunst mit cowpy generieren
    process COWPY {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

        input:
        tuple val(id), val(character), path(recording)

        output:
        path "cowpy-${recording}"

        script:
        """
        cat ${recording} | cowpy -c ${character} > cowpy-${recording}
        """
    }
    ```

=== "Vorher"

    ```groovy title="modules/cowpy.nf" linenums="1" hl_lines="7-8 11 15"
    // ASCII-Kunst mit cowpy generieren
    process COWPY {

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

Jetzt erwartet der Prozess nur noch eine Eingabe, die alles enthält, was wir ihm geben möchten.

#### 1.4.2. `map()` verwenden, um das Eingabe-Tupel zu erstellen

Wir müssen noch eine Mapping-Operation verwenden, um die Elemente aufzuzählen, die wir im Tupel an den Prozess übergeben möchten:

=== "Danach"

    ```groovy title="main.nf" linenums="5" hl_lines="3-5"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [row.id, row.character, row.recording]
            }
    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="5" hl_lines="3-6"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .multiMap { row ->
                file: row.recording
                character: row.character
            }
    ```

Du fragst dich vielleicht, warum wir nicht einfach die gesamte Groovy-Map aus `splitCsv` direkt übergeben können.
Der Grund: Wir müssen Nextflow explizit mitteilen, dass die Aufnahmedatei als Pfad behandelt werden muss (d.h. sie muss korrekt bereitgestellt werden).
Das geschieht auf der Ebene der Eingabeschnittstelle von `COWPY`, wo das `recording`-Element explizit als `path` deklariert ist.

#### 1.4.3. Den Prozessaufruf aktualisieren

Ersetzen wir abschließend die zwei separaten Eingaben im Prozessaufruf durch das einzelne Tupel, das wir gerade erstellt haben:

=== "Danach"

    ```groovy title="main.nf" linenums="5" hl_lines="7"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [row.id, row.character, row.recording]
            }

        COWPY(ch_datasheet)
    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="5" hl_lines="7"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [row.id, row.character, row.recording]
            }

        COWPY(ch_datasheet.file, ch_datasheet.character)
    ```

Das vereinfacht den Prozessaufruf etwas.

#### 1.4.4. Den Workflow ausführen

Führe den Workflow aus, um zu überprüfen, dass `COWPY` die Daten weiterhin korrekt verarbeitet:

```bash
nextflow run main.nf
```

??? success "Befehlsausgabe"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [pedantic_lovelace] DSL2 - revision: b2c3d4e5f6

    executor >  local (7)
    [5e/2a1b34] COWPY (7) [100%] 7 of 7 ✔
    ```

Die Ausgabe sind dieselben sieben `cowpy-*.txt`-Dateien wie zuvor, jetzt mit einem einfacheren Aufruf von `COWPY`.

??? abstract "Results-Verzeichnis Inhalt"

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

??? example "Inhalt von results/cowpy-guten_tag.txt"

    ```console
    $ cat results/cowpy-guten_tag.txt
     _____________________________
    / Guten Tag, wie geht es dir? \
    \ Auf Wiedersehen, bis morgen /
     -----------------------------
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

Das ist eine leichte Verbesserung gegenüber dem `multiMap`-Ansatz.
Aber wir mussten die ursprüngliche Groovy-Map immer noch entpacken, um das Eingabe-Tupel zu erstellen, und es gibt eine enge Kopplung zwischen dem Prozess und dem Datenblatt: Die `COWPY`-Eingabedefinition referenziert jetzt direkt die Spaltennamen `id`, `character` und `recording`.

```groovy
input:
tuple val(id), val(character), path(recording)
```

Wenn eine andere Person ein anders strukturiertes Datenblatt verwendet – mit zusätzlichen Spalten oder Spalten in einer anderen Reihenfolge – funktioniert dieser Prozess ohne Änderungen nicht.
Das macht den Prozess fragil, weil seine Eingabestruktur an die genaue Zusammensetzung des Datenblatts gebunden ist.

Um das zu lösen, brauchen wir eine Möglichkeit, alle Metadaten als Paket zu übergeben, ohne ihre genaue Struktur in die Prozessschnittstelle fest einzuprogrammieren.

### 1.5. Eine Meta-Map + Datei-Schnittstelle verwenden

Die Lösung besteht darin, zwei verschiedene Aspekte im Kanal zu trennen: die **Metadaten über eine Probe** und die **Datendatei** selbst.
Indem wir alle Metadaten in einer einzigen Map bündeln – der „Meta-Map" – erhalten wir ein konsistentes Tupel aus zwei Elementen, unabhängig davon, wie viele Metadatenspalten das Datenblatt enthält:

```groovy title="Syntax example"
input:
tuple val(meta), path(file)
```

Das Hinzufügen oder Entfernen von Spalten im Datenblatt ändert den Inhalt von `meta`, aber die Tupelform `[meta, file]` bleibt konstant.
Prozesse, die diese Struktur akzeptieren, müssen nicht wissen oder sich darum kümmern, wie viele Metadatenfelder vorhanden sind.

#### 1.5.1. Den Tupelinhalt in eine Meta-Map umstrukturieren

Lass uns die `map`-Operation umstrukturieren, um ein `[meta, file]`-Tupel zu erzeugen:

=== "Danach"

    ```groovy title="main.nf" linenums="5" hl_lines="4 6 8 11"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [[id: row.id, character: row.character], row.recording]
            }
            .view()

        // COWPY(ch_datasheet)  // Wird im nächsten Schritt aktualisiert

        publish:
        cowpy_art = channel.empty() // COWPY.out
    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="5" hl_lines="4 7"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [row.id, row.character, row.recording]
            }

        COWPY(ch_datasheet)

        publish:
        cowpy_art = COWPY.out
    ```

Du wirst bemerken, dass wir auch eine `view()`-Anweisung hinzugefügt, den `COWPY`-Aufruf auskommentiert und `COWPY.out` durch `channel.empty()` ersetzt haben, weil die Prozess-Eingabedefinition noch nicht zur neuen Struktur passt.

#### 1.5.2. Den Workflow ausführen, um den umstrukturierten Inhalt zu prüfen

Führe den Workflow aus, um die neue Kanalform zu sehen:

```bash
nextflow run main.nf
```

??? success "Befehlsausgabe"

    ```console title="View meta map"
     N E X T F L O W   ~  version 25.10.4

    Launching `main.nf` [lethal_booth] DSL2 - revision: 0d8f844c07

    [[id:sampleA, character:squirrel], /workspaces/training/side-quests/metadata/data/bonjour.txt]
    [[id:sampleB, character:tux], /workspaces/training/side-quests/metadata/data/guten_tag.txt]
    [[id:sampleC, character:sheep], /workspaces/training/side-quests/metadata/data/hallo.txt]
    [[id:sampleD, character:turkey], /workspaces/training/side-quests/metadata/data/hello.txt]
    [[id:sampleE, character:stegosaurus], /workspaces/training/side-quests/metadata/data/hola.txt]
    [[id:sampleF, character:moose], /workspaces/training/side-quests/metadata/data/salut.txt]
    [[id:sampleG, character:turtle], /workspaces/training/side-quests/metadata/data/ciao.txt]
    ```

Jedes Element im Kanal ist jetzt ein Tupel aus zwei Elementen: zuerst die Meta-Map, dann die Datei.

```console title="Example element structure"
[
  [id:sampleA, character:squirrel],
  /workspaces/training/side-quests/metadata/data/bonjour.txt
]
```

Wenn wir später eine `language`-Spalte zum Datenblatt hinzufügen, wird sie als `meta.language` verfügbar, ohne dass Änderungen an der Prozess-Eingabedefinition erforderlich sind.

#### 1.5.3. Den `COWPY`-Prozess aktualisieren, um die Meta-Map zu verwenden

Aktualisiere `COWPY`, um die `[meta, file]`-Tupelstruktur zu akzeptieren:

=== "Danach"

    ```groovy title="modules/cowpy.nf" linenums="1" hl_lines="7 10 14"
    // ASCII-Kunst mit cowpy generieren
    process COWPY {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

        input:
        tuple val(meta), path(input_file)

        output:
        path "cowpy-${input_file}"

        script:
        """
        cat ${input_file} | cowpy -c ${meta.character} > cowpy-${input_file}
        """
    }
    ```

=== "Vorher"

    ```groovy title="modules/cowpy.nf" linenums="1" hl_lines="7 10 14"
    // ASCII-Kunst mit cowpy generieren
    process COWPY {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

        input:
        tuple val(id), val(character), path(recording)

        output:
        path "cowpy-${recording}"

        script:
        """
        cat ${recording} | cowpy -c ${character} > cowpy-${recording}
        """
    }
    ```

Im Script-Block greift `meta.character` auf das `character`-Feld der Meta-Map zu.
Auf jedes Feld in der Meta-Map kann auf dieselbe Weise zugegriffen werden.

#### 1.5.4. Den Prozessaufruf aktualisieren

Stelle den `COWPY`-Aufruf wieder her und verbinde seine Ausgabe für die Veröffentlichung:

=== "Danach"

    ```groovy title="main.nf" linenums="5" hl_lines="7 10"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [[id: row.id, character: row.character], row.recording]
            }

        COWPY(ch_datasheet)

        publish:
        cowpy_art = COWPY.out
    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="5" hl_lines="6 8 11"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [[id: row.id, character: row.character], row.recording]
            }
            .view()

        // COWPY(ch_datasheet)  // Wird im nächsten Schritt aktualisiert

        publish:
        cowpy_art = channel.empty() // COWPY.out
    ```

Wir haben auch die Ausgabeveröffentlichung wiederhergestellt.

#### 1.5.5. Den Workflow ausführen

Führe den Workflow aus, um zu überprüfen, dass alles funktioniert:

```bash
nextflow run main.nf
```

??? success "Befehlsausgabe"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [wise_sammet] DSL2 - revision: 99797b1e92

    executor >  local (7)
    [5d/dffd4e] COWPY (7) [100%] 7 of 7 ✔
    ```

Das Ergebnisverzeichnis enthält jetzt die ASCII-Kunst-Dateien.

??? abstract "Verzeichnisinhalt"

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

??? example "Inhalt von results/cowpy-guten_tag.txt"

    ```console
    $ cat results/cowpy-guten_tag.txt
     _____________________________
    / Guten Tag, wie geht es dir? \
    \ Auf Wiedersehen, bis morgen /
     -----------------------------
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

Der Prozess empfängt jetzt alle Metadaten als Paket über `meta`, verwendet was er braucht (`meta.character`) und ignoriert den Rest.

Das ist die Standardschnittstelle, die alle [nf-core](https://nf-co.re/)-Module verwenden.
Das `tuple val(meta), path(file)`-Muster erscheint durchgängig in der nf-core-Modulbibliothek, weshalb Workflows, die diese Konvention übernehmen, nf-core-Module mit minimalem Aufwand einbinden können.

### Fazit

In diesem Abschnitt hast du gelernt:

- **Wie man Datenblätter einliest:** Mit `splitCsv` CSV-Dateien mit Kopfzeileninformationen lesen
- **Warum die Meta-Map-Konvention existiert:** Das Trennen von Metadaten und Datendateien in `[meta, file]`-Tupel hält die Kanalstruktur stabil, wenn sich das Datenblatt weiterentwickelt
- **Wie man Meta-Map-Felder innerhalb eines Prozesses verwendet:** Auf jedes Feld in der Meta-Map kann per Punktnotation im Script-Block zugegriffen werden

---

## 2. Weitere Metadaten-Manipulationen

Jetzt, wo die Meta-Map-Schnittstelle eingerichtet ist, können wir sie anreichern, während Daten durch die Pipeline fließen.

Wir werden ein Tool namens [`langid`](https://github.com/saffsd/langid.py) verwenden, um die Sprache in jeder Aufnahmedatei zu identifizieren.
Bei einem Textausschnitt gibt es eine Sprachvorhersage und einen Wahrscheinlichkeitswert nach `stdout` aus.

### 2.1. Einen Spracherkennungsschritt hinzufügen

Wir stellen dir ein vorgefertigtes Prozessmodul namens `IDENTIFY_LANGUAGE` zur Verfügung, das das `langid`-Tool kapselt.

Öffne die Moduldatei, um den Code zu untersuchen:

```groovy title="modules/langid.nf" linenums="1" hl_lines="7 10"
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

Die Eingabedefinition verwendet dieselbe `tuple val(meta), path(file)`-Struktur, die wir in Abschnitt 1 aufgebaut haben, sodass `ch_datasheet` direkt in diesen Prozess einfließen kann.

Die Ausgabe fügt `stdout` als drittes Element hinzu: Das erfasst die Sprachvorhersage, die `langid` auf der Konsole ausgibt.
Der `sed`-Befehl entfernt den Wahrscheinlichkeitswert und den abschließenden Zeilenumbruch und lässt nur den zweistelligen Sprachcode übrig.

#### 2.1.1. Einen Aufruf von `IDENTIFY_LANGUAGE` hinzufügen

Füge das `IDENTIFY_LANGUAGE`-Prozessmodul ein und rufe es auf dem Datenblatt-Kanal auf:

=== "Danach"

    ```groovy title="main.nf" linenums="1" hl_lines="4 14-16"
    #!/usr/bin/env nextflow

    include { COWPY } from './modules/cowpy.nf'
    include { IDENTIFY_LANGUAGE } from './modules/langid.nf'

    workflow {
        main:
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [[id: row.id, character: row.character], row.recording]
            }

        // langid ausführen, um die Sprache jeder Begrüßung zu identifizieren
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out.view()

        COWPY(ch_datasheet)
    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="1"
    #!/usr/bin/env nextflow

    include { COWPY } from './modules/cowpy.nf'

    workflow {
        main:
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [[id: row.id, character: row.character], row.recording]
            }

        COWPY(ch_datasheet)
    ```

Die Hauptausgabe dieses Prozesses ist nur ein String, daher gibt es keine Ausgabedateien zu veröffentlichen.
Stattdessen verwenden wir `IDENTIFY_LANGUAGE.out.view()`, um die Ergebnisse der Operation anzuzeigen.

#### 2.1.2. Den Workflow ausführen

Führe den Workflow aus, um die Sprachidentifikation zu erzeugen. Verwende `-resume`, um die `COWPY`-Aufgaben nicht erneut auszuführen:

```bash
nextflow run main.nf -resume
```

??? success "Befehlsausgabe"

    ```console
     N E X T F L O W   ~  version 25.10.4

    Launching `main.nf` [voluminous_mcnulty] DSL2 - revision: f9bcfebabb

    executor >  local (14)
    [5d/dffd4e] COWPY (7)             [100%] 7 of 7, cached: 7 ✔
    [4e/f722fe] IDENTIFY_LANGUAGE (7) [100%] 7 of 7 ✔
    [[id:sampleA, character:squirrel], /workspaces/training/side-quests/metadata/work/eb/f7148ebdd898fbe1136bec6a714acb/bonjour.txt, fr]
    [[id:sampleB, character:tux], /workspaces/training/side-quests/metadata/work/16/71d72410952c22cd0086d9bca03680/guten_tag.txt, de]
    [[id:sampleD, character:turkey], /workspaces/training/side-quests/metadata/work/c4/b7562adddc1cc0b7d414ec45d436eb/hello.txt, en]
    [[id:sampleC, character:sheep], /workspaces/training/side-quests/metadata/work/ea/04f5d979429e4455e14b9242fb3b45/hallo.txt, de]
    [[id:sampleF, character:moose], /workspaces/training/side-quests/metadata/work/5a/6c2b84bf8fadb98e28e216426be079/salut.txt, fr]
    [[id:sampleE, character:stegosaurus], /workspaces/training/side-quests/metadata/work/af/ee7c69bcab891c40d0529305f6b9e7/hola.txt, es]
    [[id:sampleG, character:turtle], /workspaces/training/side-quests/metadata/work/4e/f722fe47271ba7ebcd69afa42964ca/ciao.txt, it]
    ```

Wir haben jetzt eine Sprachvorhersage für jede Datei im Datensatz.

Das Ausgabe-Tupel besteht aus `[meta, file, lang_id]`, d.h. die Meta-Map und die Datei werden zusammen mit dem neuen Ergebnis weitergegeben.

!!! note "Hinweis"

    Dieses Muster, die Meta-Map mit den Ergebnissen verknüpft zu halten, macht es einfacher, Ergebnisse später kanalübergreifend zu verknüpfen.
    Man kann sich nicht auf die Reihenfolge der Elemente in Kanälen verlassen, um Daten korrekt zuzuordnen.
    Stattdessen müssen Schlüssel verwendet werden.
    Meta-Maps bieten dafür eine ideale Struktur.

    Dieser Anwendungsfall wird ausführlich in der Side Quest [Splitting & Grouping](../splitting_and_grouping/index.md) erkundet.

### 2.2. Metadaten mit Prozessausgaben anreichern

Die Sprachvorhersage ist selbst eine Form von Metadaten über den Inhalt der Datei.
Anstatt sie als separates Element zu behalten, falten wir sie zurück in die Meta-Map.

#### 2.2.1. Eine neue und erweiterte Meta-Map erstellen

Wir können eine neue Meta-Map erstellen, die die ursprüngliche ersetzt, indem wir den Groovy-`+`-Operator verwenden:

=== "Danach"

    ```groovy title="main.nf" linenums="14" hl_lines="3-7"
        // langid ausführen, um die Sprache jeder Begrüßung zu identifizieren
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out
            .map { meta, file, lang_id ->
                [meta + [lang: lang_id], file]
            }
            .view()
    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="14" hl_lines="3"
        // langid ausführen, um die Sprache jeder Begrüßung zu identifizieren
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out.view()
    ```

Das Herzstück dieser Operation ist `#!groovy meta + [lang: lang_id]`.

Dieser Code erstellt eine temporäre Map mit einem einzigen Schlüssel-Wert-Paar, das den Sprachcode enthält (`[lang: lang_id]`), und verwendet dann den Groovy-`+`-Operator, um sie mit der ursprünglichen `meta`-Map zu kombinieren, die die vorhandenen Metadaten enthält. Das ergibt eine neue, erweiterte Meta-Map.

Eine ausführlichere Erklärung findest du im Kasten unten.

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
    new_map = map1 + [lang: lang_id]
    ```

    Hier erstellt `[lang: lang_id]` eine neue namenlose Map auf der Stelle, und `map1 + ` führt `map1` mit der neuen namenlosen Map zusammen, was denselben `new_map`-Inhalt wie zuvor erzeugt.

    Praktisch, oder?

    **Jetzt übertragen wir das in den Kontext einer Nextflow `channel.map()`-Operation.**

    Der Code wird zu:

    ```groovy
    .map { map1, lang_id ->
        map1 + [lang: lang_id]
    }
    ```

    Das macht Folgendes:

    - `#!groovy map1, lang_id ->` nimmt die zwei Elemente im Tupel
    - `#!groovy map1 + [lang: lang_id]` erstellt die neue Map wie oben beschrieben

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

    Falls es schwer nachzuvollziehen ist, warum sich `file` in der `map`-Operation scheinbar bewegt, stell dir vor, dass statt `#!groovy [meta + [lang: lang_id], file]` diese Zeile `[new_map, file]` lautet.
    Das sollte klarer machen, dass wir `file` einfach an seiner ursprünglichen zweiten Position im Tupel belassen. Wir haben nur den `new_info`-Wert in die Map an erster Position eingefaltet.

    **Und das bringt uns zurück zur `tuple val(meta), path(file)`-Kanalstruktur!**

#### 2.2.2. Den Workflow ausführen

Sobald du sicher bist, dass du verstehst, was dieser Code tut, führe den Workflow aus, um zu sehen, ob es funktioniert hat:

```bash
nextflow run main.nf -resume
```

??? success "Befehlsausgabe"

    ```console
     N E X T F L O W   ~  version 25.10.4

    Launching `main.nf` [cheeky_fermat] DSL2 - revision: d096281ee4

    [5d/dffd4e] COWPY (7)             [100%] 7 of 7, cached: 7 ✔
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

!!! tip "Schlüssel aus einer Meta-Map entfernen"

    Du kannst einen Schlüssel aus einer Meta-Map mit der Groovy-Methode [`subMap`](https://docs.groovy-lang.org/latest/html/groovy-jdk/java/util/Map.html#subMap(java.util.Collection)) entfernen, die eine neue Map zurückgibt, die nur die angegebenen Schlüssel enthält:

    ```groovy
    meta.subMap(['id', 'character'])  // gibt eine Map nur mit 'id' und 'character' zurück
    ```

    Das ist nützlich, wenn ein nachgelagerter Prozess oder ein Modul nicht alle Felder benötigt, die sich in der Meta-Map angesammelt haben.

### 2.3. Eine Sprachgruppe mit Bedingungen zuweisen

Mit der Sprachvorhersage in der Meta-Map können wir weitere Metadaten daraus ableiten.
Die Sprachen in unserem Datensatz lassen sich in zwei Familien einteilen: Germanisch (Englisch, Deutsch) und Romanisch (Französisch, Spanisch, Italienisch).
Das Hinzufügen eines `lang_group`-Felds macht diese Klassifizierung für nachgelagerte Schritte verfügbar.

#### 2.3.1. Eine `map`-Operation mit der bedingten Logik hinzufügen

Wir verwenden eine zweite `map`-Operation mit bedingter Logik, um die Sprachfamilie zuzuweisen:

```groovy
.map { meta, file ->

    // Bedingte Logik zur Definition von lang_group kommt hier hin

    [meta + [lang_group: lang_group], file]
}
```

Hier ist die anzuwendende Logik:

- Beginne mit `lang_group = 'unknown'` als Standardwert.
- Wenn `meta.lang` `'de'` oder `'en'` ist, setze `lang_group` auf `'germanic'`.
- Sonst wenn `meta.lang` in `['fr', 'es', 'it']` enthalten ist, setze `lang_group` auf `'romance'`.

!!! tip "Tipp"

    Du kannst auf den Wert von `lang` innerhalb der map-Operation mit `meta.lang` zugreifen.

Nimm die folgenden Änderungen am Workflow vor:

=== "Danach"

    ```groovy title="main.nf" linenums="14" hl_lines="7-19 21"
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

        ch_languages.view()
    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="14" hl_lines="7"
        // langid ausführen, um die Sprache jeder Begrüßung zu identifizieren
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out
            .map { meta, file, lang_id ->
                [meta + [lang: lang_id], file]
            }
            .view()
    ```

Wichtige Punkte:

- `def lang_group = "unknown"` initialisiert die Variable mit einem sicheren Standardwert.
- Die `if / else if`-Struktur behandelt die zwei Sprachfamilien; alles andere bleibt `'unknown'`.
- `#!groovy .set { ch_languages }` gibt dem resultierenden Kanal einen Namen für den nächsten Schritt.

<!-- TODO (future) Add note/links to relevant docs in additional resources section -->

#### 2.3.2. Den Workflow ausführen

Führe den Workflow aus, um zu überprüfen, dass es funktioniert:

```bash
nextflow run main.nf -resume
```

??? success "Befehlsausgabe"

    ```console
     N E X T F L O W   ~  version 25.10.4

    Launching `main.nf` [wise_almeida] DSL2 - revision: 46778c3cd0

    [5d/dffd4e] COWPY (7)             [100%] 7 of 7, cached: 7 ✔
    [da/652cc6] IDENTIFY_LANGUAGE (7) [100%] 7 of 7, cached: 7 ✔
    [[id:sampleA, character:squirrel, lang:fr, lang_group:romance], /workspaces/training/side-quests/metadata/data/bonjour.txt]
    [[id:sampleB, character:tux, lang:de, lang_group:germanic], /workspaces/training/side-quests/metadata/data/guten_tag.txt]
    [[id:sampleC, character:sheep, lang:de, lang_group:germanic], /workspaces/training/side-quests/metadata/data/hallo.txt]
    [[id:sampleD, character:turkey, lang:en, lang_group:germanic], /workspaces/training/side-quests/metadata/data/hello.txt]
    [[id:sampleE, character:stegosaurus, lang:es, lang_group:romance], /workspaces/training/side-quests/metadata/data/hola.txt]
    [[id:sampleF, character:moose, lang:fr, lang_group:romance], /workspaces/training/side-quests/metadata/data/salut.txt]
    [[id:sampleG, character:turtle, lang:it, lang_group:romance], /workspaces/training/side-quests/metadata/data/ciao.txt]
    ```

Die Meta-Map enthält jetzt vier Felder: `id`, `character`, `lang` und `lang_group`.
Die Kanalstruktur ist weiterhin `[meta, file]`.

### 2.4. Metadaten zur Benennung und Organisation von Ausgaben verwenden

Mit `lang` und `lang_group` in der Meta-Map können wir sie verwenden, um einen Sprachcode zu den Ausgabedateinamen hinzuzufügen und die Dateien in Unterverzeichnisse nach Sprachfamilie zu organisieren.

Dazu sind drei Änderungen erforderlich: den `COWPY`-Prozess aktualisieren, um seine Ausgabe umzubenennen und `meta` in der Ausgabe einzuschließen, den `COWPY`-Aufruf auf `ch_languages` umstellen und den Output-Block aktualisieren, um den Unterverzeichnispfad anzugeben.

#### 2.4.1. Den `COWPY`-Prozess aktualisieren

Benenne die Ausgabedatei mit dem Sprachcode aus der Meta-Map um und füge `meta` zur Ausgabe hinzu, damit der Output-Block auf `lang_group` für die Unterverzeichnis-Zuordnung zugreifen kann:

=== "Danach"

    ```groovy title="modules/cowpy.nf" linenums="9" hl_lines="2 6"
        output:
        tuple val(meta), path("${meta.lang}-${input_file}")

        script:
        """
        cat ${input_file} | cowpy -c ${meta.character} > ${meta.lang}-${input_file}
        """
    ```

=== "Vorher"

    ```groovy title="modules/cowpy.nf" linenums="9" hl_lines="2 6"
        output:
        path "cowpy-${input_file}"

        script:
        """
        cat ${input_file} | cowpy -c ${meta.character} > cowpy-${input_file}
        """
    ```

Das zeigt, wie wir andere Metadatenfelder nutzen können, um das Verhalten eines Prozesses anzupassen, ohne die Eingabedefinition ändern zu müssen.

#### 2.4.2. Den `COWPY`-Aufruf auf `ch_languages` umstellen

Ersetze `COWPY(ch_datasheet)` durch `COWPY(ch_languages)`:

=== "Danach"

    ```groovy title="main.nf" linenums="32" hl_lines="3"
        .set { ch_languages }

        COWPY(ch_languages)

        publish:
        cowpy_art = COWPY.out
    }
    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="32" hl_lines="3 5"
        .set { ch_languages }

        ch_languages.view()

        COWPY(ch_datasheet)

        publish:
        cowpy_art = COWPY.out
    }
    ```

Wir entfernen auch die `ch_languages.view()`-Zeile, da wir den Kanalinhalt nicht mehr prüfen müssen.

#### 2.4.3. Den Output-Block aktualisieren

Füge eine `path`-Closure zum `output {}`-Block hinzu, um jede Datei in ihr Sprachgruppen-Unterverzeichnis zu leiten:

=== "Danach"

    ```groovy title="main.nf" linenums="40" hl_lines="3"
    output {
        cowpy_art {
            path { meta, file -> meta.lang_group }
        }
    }
    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="40" hl_lines="2 3"
    output {
        cowpy_art {
        }
    }
    ```

Das zeigt, wie wir Metadaten verwenden können, um Ausgaben mit großer Flexibilität zu organisieren.

#### 2.4.4. Die vollständige Pipeline ausführen

Lösche die vorherigen Ergebnisse und führe die vollständige Pipeline aus:

```bash
rm -r results
nextflow run main.nf
```

??? success "Befehlsausgabe"

    ```console
     N E X T F L O W   ~  version 25.10.4

    Launching `main.nf` [suspicious_crick] DSL2 - revision: 25541014c5

    executor >  local (14)
    [5d/dffd4e] IDENTIFY_LANGUAGE (7) [100%] 7 of 7 ✔
    [e7/317c18] COWPY (7)             [100%] 7 of 7 ✔
    ```

Das Ergebnisverzeichnis ist jetzt nach Sprachfamilie organisiert, wobei jede Datei nach der erkannten Sprache benannt ist:

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

Die `path`-Closure im `output {}`-Block empfängt jedes `[meta, file]`-Tupel und gibt `meta.lang_group` als Unterverzeichnisnamen zurück.
Der Dateiname selbst kommt von der Prozessausgabe (`#!groovy "${meta.lang}-${input_file}"`).
Beide Metadaten (Sprachcode und Sprachgruppe) stammen aus der in diesem Abschnitt aufgebauten angereicherten Meta-Map.

### Fazit

In diesem Abschnitt hast du gelernt:

- **Wie man die Meta-Map mit Prozessausgaben anreichert:** Das Hinzufügen neuer Schlüssel mit `#!groovy meta + [key: value]` erhält die `[meta, file]`-Kanalstruktur, während die Metadaten erweitert werden.
- **Wie man Metadaten aus Metadaten ableitet:** Bedingte Logik innerhalb einer `map`-Operation kann neue Felder aus bestehenden berechnen.
- **Wie man Metadaten zur Ausgabeorganisation verwendet:** Die `path`-Closure im `output {}`-Block kann aus der Meta-Map lesen, um Dateien in Unterverzeichnisse zu leiten.

---

## 3. Überlegungen zur Robustheit

Wenn Metadatenwerte das Prozessverhalten steuern, können fehlende oder unvollständige Daten Probleme verursachen, die schwer zu diagnostizieren sind.
Hier erfährst du, was zu erwarten ist und wie du damit umgehst.

### 3.1. Was passiert, wenn ein erforderliches Metadatenfeld fehlt

Der `character`-Wert ist erforderlich, damit der `COWPY`-Prozess ein gültiges Ergebnis erzeugt.
Das Fehlverhalten hängt davon ab, ob die Spalte im Datenblatt vorhanden, aber leer ist, oder ob sie ganz fehlt.

#### 3.1.1. Die Spalte existiert, aber ein Wert ist leer

Angenommen, ein Eintrag im Datenblatt hat ein leeres `character`-Feld:

```csv title="datasheet.csv" linenums="1" hl_lines="2"
id,character,recording
sampleA,,/workspaces/training/side-quests/metadata/data/bonjour.txt
sampleB,tux,/workspaces/training/side-quests/metadata/data/guten_tag.txt
...
```

Der `character`-Schlüssel wird für alle Einträge beim Einlesen des Datenblatts erstellt, aber `meta.character` für `sampleA` ist ein leerer String.
Wenn Nextflow `#!groovy ${meta.character}` in den Befehl einsetzt, erhält das `COWPY`-Tool ein leeres Argument für `-c` und schlägt fehl:

??? failure "Befehlsausgabe"

    ```console hl_lines="8 11 16 28"
     N E X T F L O W   ~  version 25.10.4

    Launching `main.nf` [marvelous_hirsch] DSL2 - revision: 0dfeee3cc1

    executor >  local (9)
    [c1/c5dd4f] process > IDENTIFY_LANGUAGE (7) [ 85%] 6 of 7
    [d3/b7c415] process > COWPY (2)             [  0%] 0 of 6
    ERROR ~ Error executing process > 'COWPY (1)'

    Caused by:
      Process `COWPY (1)` terminated with an error exit status (2)


    Command executed:

      cat bonjour.txt | cowpy -c  > fr-bonjour.txt

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

Die Fehlermeldung (`expected one argument`) weist auf das leere `-c`-Flag hin.
Ein Blick in die `.command.sh`-Datei im work-Verzeichnis bestätigt, dass der Befehl mit einem leeren Wert ausgeführt wurde.

#### 3.1.2. Die Spalte existiert nicht im Datenblatt

Wenn die `character`-Spalte vollständig fehlt:

```csv title="datasheet.csv" linenums="1"
id,recording
sampleA,/workspaces/training/side-quests/metadata/data/bonjour.txt
sampleB,/workspaces/training/side-quests/metadata/data/guten_tag.txt
...
```

Der `character`-Schlüssel wird in der Meta-Map nie erstellt.
Wenn der Prozess-Script `#!groovy ${meta.character}` auswertet, gibt der fehlende Schlüssel `null` zurück, und Nextflow setzt buchstäblich den String `null` in den Befehl ein:

??? failure "Befehlsausgabe"

    ```console hl_lines="8 11 16"
     N E X T F L O W   ~  version 25.10.4

    Launching `main.nf` [jovial_bohr] DSL2 - revision: eaaf375827

    executor >  local (9)
    [0d/ada9db] process > IDENTIFY_LANGUAGE (5) [ 85%] 6 of 7
    [06/28065f] process > COWPY (2)             [  0%] 0 of 6
    ERROR ~ Error executing process > 'COWPY (2)'

    Caused by:
      Process `COWPY (2)` terminated with an error exit status (1)


    Command executed:

      cat guten_tag.txt | cowpy -c null > de-guten_tag.txt

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

Das `cowpy -c null` im ausgeführten Befehl ist der diagnostische Hinweis.

### 3.2. Strategien für den Umgang mit fehlenden Metadaten

Es gibt zwei sich ergänzende Ansätze, um Workflows robuster gegen fehlende Metadaten zu machen.

**1. Eingabevalidierung**

Die zuverlässigste Lösung ist, das Datenblatt zu validieren, bevor die Verarbeitung beginnt, damit Probleme frühzeitig mit einer klaren Fehlermeldung erkannt werden, anstatt als kryptischer Prozessfehler mitten in der Ausführung aufzutauchen.
Das [Hello nf-core](../../hello_nf-core/05_input_validation.md)-Training zeigt, wie man Eingabevalidierung mit dem nf-schema-Plugin hinzufügt. <!-- TODO (future) pending a proper Validation side quest -->

**2. Explizite Prozesseingaben für erforderliche Werte**

Wenn die Prozessschnittstelle selbst kommunizieren soll, dass ein bestimmter Wert erforderlich ist, kannst du ihn als explizite Eingabe aus der Meta-Map extrahieren:

=== "Prozessdefinition"

    ```groovy title="modules/cowpy.nf" linenums="6"
    input:
    tuple val(meta), val(character), path(input_file)
    ```

=== "Workflow-Aufruf"

    ```groovy title="main.nf"
    COWPY(ch_languages.map { meta, file -> [meta, meta.character, file] })
    ```

Dieser Ansatz macht `character` zu einem sichtbaren, erforderlichen Teil des Prozessvertrags.
Wer das Modul liest, sieht sofort, dass ein Charakterwert angegeben werden muss.
Fehlt das Feld, schlägt der Workflow klar auf Kanalebene fehl, bevor der Prozess überhaupt ausgeführt wird.

Das verdeutlicht ein nützliches Designprinzip:

**Verwende die Meta-Map für optionale oder beschreibende Informationen; extrahiere erforderliche Werte als explizite Eingaben.**

Die Meta-Map hält Kanalstrukturen sauber und stabil, aber für Werte, die von einem Prozess wirklich benötigt werden, verbessert das Herauslösen als benannte Eingaben die Klarheit und macht das Modul in anderen Kontexten einfacher korrekt zu verwenden.

### Fazit

In diesem Abschnitt hast du gesehen:

- **Wie sich fehlende Metadaten äußern:** Ein leeres Feld erzeugt ein leeres Argument; ein fehlendes Feld erzeugt `null`, das buchstäblich in den Befehl eingesetzt wird.
- **Zwei sich ergänzende Strategien:** Eingabevalidierung, um Probleme frühzeitig zu erkennen, und explizite Prozesseingaben, um Anforderungen klar zu kommunizieren.

---

## Zusammenfassung

In dieser Side Quest hast du erkundet, wie man effektiv mit Metadaten in Nextflow-Workflows arbeitet.

Das „Meta-Map + Datendatei"-Tupel-Muster ist eine zentrale Konvention in Nextflow und bietet mehrere Vorteile gegenüber der Übergabe von Metadaten als einzelne Werte:

- Die Kanalstruktur bleibt stabil, wenn sich das Datenblatt weiterentwickelt
- Das Prozessverhalten kann pro Probe angepasst werden, ohne Feldnamen fest einzuprogrammieren
- Metadaten sind während der gesamten Pipeline für Benennung, Gruppierung und Organisation von Ausgaben verfügbar
- Module, die für diese Schnittstelle geschrieben wurden, sind austauschbar – einschließlich nf-core-Module

### Wichtige Muster

1.  **Metadaten lesen und strukturieren:** Eine CSV-Datei einlesen und eine Meta-Map erstellen.

    ```groovy
    channel.fromPath('datasheet.csv')
        .splitCsv(header: true)
        .map { row ->
            [ [id: row.id, character: row.character], row.recording ]
        }
    ```

2.  **Metadaten während des Workflows erweitern:** Neue Schlüssel aus Prozessausgaben oder abgeleiteter Logik hinzufügen.

    ```groovy
    // Aus einer Prozessausgabe
    .map { meta, file, lang ->
        [ meta + [lang: lang], file ]
    }

    // Aus bedingter Logik
    .map { meta, file ->
        def lang_group = "unknown"
        if (meta.lang in ["de", "en"]) { lang_group = "germanic" }
        else if (meta.lang in ["fr", "es", "it"]) { lang_group = "romance" }
        [ meta + [lang_group: lang_group], file ]
    }
    ```

3.  **Metadaten innerhalb eines Prozesses verwenden:** Auf jedes Feld per Punktnotation im Script-Block zugreifen.

    ```groovy
    cat ${input_file} | cowpy -c ${meta.character} > ${meta.lang}-${input_file}
    ```

4.  **Ausgaben nach Metadatenwert organisieren:** Die `path`-Closure im `output {}`-Block verwenden.

    ```groovy
    output {
        cowpy_art {
            path { meta, file -> meta.lang_group }
        }
    }
    ```

### Weitere Ressourcen

- [map-Operator](https://www.nextflow.io/docs/latest/operator.html#map)
- [multiMap-Operator](https://www.nextflow.io/docs/latest/reference/operator.html#multimap)
- [stdout-Ausgabe-Qualifier](https://www.nextflow.io/docs/latest/process.html#outputs)

---

## Wie geht es weiter?

Kehre zum [Menü der Side Quests](../index.md) zurück oder klicke auf die Schaltfläche unten rechts auf der Seite, um zum nächsten Thema in der Liste zu wechseln.
