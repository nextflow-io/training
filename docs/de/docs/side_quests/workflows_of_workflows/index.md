# Workflows of Workflows

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Wenn du eine Pipeline entwickelst, erstellst du oft ähnliche Abfolgen von Prozessen für verschiedene Datentypen oder Analyseschritte. Du könntest diese Prozessabfolgen kopieren und einfügen, was zu doppeltem Code führt, der schwer zu pflegen ist – oder du erstellst einen riesigen Workflow, der schwer zu verstehen und zu ändern ist.

Eine der mächtigsten Funktionen von Nextflow ist die Möglichkeit, komplexe Pipelines aus kleineren, wiederverwendbaren Workflow-Modulen zusammenzusetzen. Dieser modulare Ansatz macht Pipelines einfacher zu entwickeln, zu testen und zu pflegen.

### Lernziele

In diesem Side Quest erkunden wir, wie man Workflow-Module entwickelt, die separat getestet und verwendet werden können, diese Module zu einer größeren Pipeline zusammensetzt und den Datenfluss zwischen Modulen verwaltet.

Am Ende dieses Side Quests kannst du:

- Komplexe Pipelines in logische, wiederverwendbare Einheiten aufteilen
- Jedes Workflow-Modul unabhängig testen
- Workflows kombinieren, um neue Pipelines zu erstellen
- Gemeinsame Workflow-Module über verschiedene Pipelines hinweg teilen
- Deinen Code wartbarer und leichter verständlich machen

Diese Fähigkeiten helfen dir, komplexe Pipelines zu erstellen und dabei eine saubere, wartbare Code-Struktur beizubehalten.

### Voraussetzungen

Bevor du diesen Side Quest angehst, solltest du:

- Das Tutorial [Hello Nextflow](../../hello_nextflow/index.md) oder einen gleichwertigen Einsteigerkurs abgeschlossen haben.
- Mit grundlegenden Nextflow-Konzepten und -Mechanismen vertraut sein (Prozesse, Kanäle, Operatoren, Module)

---

## 0. Erste Schritte

#### Öffne den Training-Codespace

Falls du das noch nicht getan hast, öffne die Trainingsumgebung wie in der [Umgebung einrichten](../../envsetup/index.md) beschrieben.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Wechsle in das Projektverzeichnis

Wechsle in das Verzeichnis, in dem sich die Dateien für dieses Tutorial befinden.

```bash
cd side-quests/workflows_of_workflows
```

Du kannst VSCode so einstellen, dass es sich auf dieses Verzeichnis konzentriert:

```bash
code .
```

Der Editor öffnet sich mit dem Projektverzeichnis im Fokus.

#### Schau dir die Materialien an

Du findest ein `modules`-Verzeichnis mit Prozessdefinitionen, ein `workflows`-Verzeichnis mit zwei vorgefertigten Workflow-Skripten und eine `main.nf`-Datei, die du schrittweise aktualisieren wirst:

```console title="Directory contents"
├── main.nf
├── workflows/
│   ├── greeting.nf              # Standalone greeting workflow (to be made composable)
│   └── transform.nf             # Standalone transform workflow (to be made composable)
└── modules/
    ├── say_hello.nf             # Creates a greeting (from Hello Nextflow)
    ├── say_hello_upper.nf       # Converts to uppercase (from Hello Nextflow)
    ├── timestamp_greeting.nf    # Adds timestamps to greetings
    ├── validate_name.nf         # Validates input names
    └── reverse_text.nf          # Reverses text content
```

Das `modules/`-Verzeichnis enthält die einzelnen Prozessdefinitionen, und das `workflows/`-Verzeichnis enthält die zwei vorgefertigten Workflow-Skripte, mit denen du in diesem Side Quest arbeiten wirst.

#### Schau dir die Aufgabe an

Deine Aufgabe ist es, diese Module zu zwei separaten Workflows zusammenzusetzen, die wir dann zu einem Haupt-Workflow kombinieren:

- Ein `GREETING_WORKFLOW`, der Namen validiert, Begrüßungen erstellt und Zeitstempel hinzufügt
- Ein `TRANSFORM_WORKFLOW`, der Text in Großbuchstaben umwandelt und umkehrt

<!-- TODO: give a bit more details, similar to how it's done in the Metadata side quest -->

#### Bereitschafts-Checkliste

Bereit, loszulegen?

- [ ] Ich verstehe das Ziel dieses Kurses und seine Voraussetzungen
- [ ] Mein Codespace läuft
- [ ] Ich habe mein Arbeitsverzeichnis entsprechend eingestellt
- [ ] Ich verstehe die Aufgabe

Wenn du alle Punkte abhaken kannst, kann es losgehen.

---

## 1. Den Greeting Workflow zur Pipeline hinzufügen

Der Greeting Workflow validiert Namen und generiert Begrüßungen mit Zeitstempel.

### 1.1. Den Greeting Workflow ansehen und ausführen

Öffne `workflows/greeting.nf` und schau dir den Code an:

```groovy title="workflows/greeting.nf" linenums="1"
include { VALIDATE_NAME } from '../modules/validate_name'
include { SAY_HELLO } from '../modules/say_hello'
include { TIMESTAMP_GREETING } from '../modules/timestamp_greeting'

workflow {
    main:
    names_ch = channel.of('Alice', 'Bob', 'Charlie')

    // Prozesse verketten: validieren -> Begrüßung erstellen -> Zeitstempel hinzufügen
    validated_ch = VALIDATE_NAME(names_ch)
    greetings_ch = SAY_HELLO(validated_ch)
    timestamped_ch = TIMESTAMP_GREETING(greetings_ch)

    publish:
    greetings = greetings_ch
    timestamped = timestamped_ch
}

output {
    greetings {
    }
    timestamped {
    }
}
```

Dies ist ein vollständiger, in sich geschlossener Workflow mit der gleichen Struktur, die du im Tutorial 'Hello Nextflow' gesehen hast.
Er kodiert die Eingabenamen fest, verkettet drei Prozesse und veröffentlicht zwei Ausgaben.

Führe ihn aus, um zu überprüfen, ob alles funktioniert:

```bash
nextflow run workflows/greeting.nf
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W  ~  version 24.10.0
    Launching `workflows/greeting.nf` [peaceful_montalcini] DSL2 - revision: 90f61b7093
    executor >  local (9)
    [51/4f980f] process > VALIDATE_NAME (validating Bob)                    [100%] 3 of 3 ✔
    [2b/dd8dc2] process > SAY_HELLO (greeting Bob)                          [100%] 3 of 3 ✔
    [8e/882565] process > TIMESTAMP_GREETING (adding timestamp to greeting) [100%] 3 of 3 ✔
    ```

Damit er mit anderen Workflows kombinierbar wird, müssen einige Dinge geändert werden.

### 1.2. Den Workflow kombinierbar machen

Um einen Workflow kombinierbar zu machen, müssen vier Dinge geändert werden:
Der Workflow bekommt einen Namen, Eingaben werden in einen `take:`-Block verschoben, Ausgaben in einen `emit:`-Block,
und die eigenständigen `publish:`/`output {}`-Blöcke werden entfernt (sie gehören in den Entry Workflow).

Lass uns diese Änderungen Schritt für Schritt durchgehen.

#### 1.2.1. Den Workflow benennen

Gib dem Workflow einen Namen, damit er aus einem übergeordneten Workflow importiert werden kann.

=== "Danach"

    ```groovy title="workflows/greeting.nf" linenums="5" hl_lines="1"
    workflow GREETING_WORKFLOW {
    ```

=== "Vorher"

    ```groovy title="workflows/greeting.nf" linenums="5" hl_lines="1"
    workflow {
    ```

Mit einem Namen kann der Workflow in andere Skripte importiert werden.

#### 1.2.2. Eingaben mit `take:` deklarieren

Ersetze die fest kodierte Kanal-Deklaration durch einen `take:`-Block, der die erwarteten Eingaben des Workflows deklariert.
Der `take:`-Block kommt vor `main:`, und die Zeile `names_ch = channel.of(...)` wird entfernt.

=== "Danach"

    ```groovy title="workflows/greeting.nf" linenums="5" hl_lines="2 3 5"
    workflow GREETING_WORKFLOW {
        take:
        names_ch // Eingabekanal mit Namen

        main:
        // Prozesse verketten: validieren -> Begrüßung erstellen -> Zeitstempel hinzufügen
        validated_ch = VALIDATE_NAME(names_ch)
        greetings_ch = SAY_HELLO(validated_ch)
        timestamped_ch = TIMESTAMP_GREETING(greetings_ch)
    ```

=== "Vorher"

    ```groovy title="workflows/greeting.nf" linenums="5"
    workflow GREETING_WORKFLOW {
        main:
        names_ch = channel.of('Alice', 'Bob', 'Charlie')

        // Prozesse verketten: validieren -> Begrüßung erstellen -> Zeitstempel hinzufügen
        validated_ch = VALIDATE_NAME(names_ch)
        greetings_ch = SAY_HELLO(validated_ch)
        timestamped_ch = TIMESTAMP_GREETING(greetings_ch)
    ```

Der `take:`-Block deklariert den Kanal nur mit seinem Namen – was genau hineinfließt, wird vom übergeordneten Workflow festgelegt.

#### 1.2.3. Ausgaben mit `emit:` deklarieren

Ersetze den `publish:`-Abschnitt und entferne den `output {}`-Block, und ersetze sie durch einen `emit:`-Block, der die Ausgaben benennt.

=== "Danach"

    ```groovy title="workflows/greeting.nf" linenums="14" hl_lines="2 3 4"

        emit:
        greetings = greetings_ch // Ursprüngliche Begrüßungen
        timestamped = timestamped_ch // Begrüßungen mit Zeitstempel
    }
    ```

=== "Vorher"

    ```groovy title="workflows/greeting.nf" linenums="14"

        publish:
        greetings = greetings_ch
        timestamped = timestamped_ch
    }

    output {
        greetings {
        }
        timestamped {
        }
    }
    ```

Der `emit:`-Block stellt benannte Ausgaben bereit, auf die übergeordnete Workflows über `GREETING_WORKFLOW.out.greetings` und `GREETING_WORKFLOW.out.timestamped` zugreifen können.

#### 1.2.4. Das Ergebnis überprüfen und testen

Nach allen drei Änderungen sollte die vollständige Datei so aussehen:

```groovy title="workflows/greeting.nf" linenums="1" hl_lines="5 6 7 9 15 16 17"
include { VALIDATE_NAME } from '../modules/validate_name'
include { SAY_HELLO } from '../modules/say_hello'
include { TIMESTAMP_GREETING } from '../modules/timestamp_greeting'

workflow GREETING_WORKFLOW {
    take:
    names_ch // Eingabekanal mit Namen

    main:
    // Prozesse verketten: validieren -> Begrüßung erstellen -> Zeitstempel hinzufügen
    validated_ch = VALIDATE_NAME(names_ch)
    greetings_ch = SAY_HELLO(validated_ch)
    timestamped_ch = TIMESTAMP_GREETING(greetings_ch)

    emit:
    greetings = greetings_ch // Ursprüngliche Begrüßungen
    timestamped = timestamped_ch // Begrüßungen mit Zeitstempel
}
```

Versuche jetzt, ihn direkt auszuführen:

```bash
nextflow run workflows/greeting.nf
```

??? failure "Befehlsausgabe"

    ```console
    N E X T F L O W  ~  version 24.10.0
    Launching `workflows/greeting.nf` [high_brahmagupta] DSL2 - revision: 8f5857af25
    No entry workflow specified
    ```

Das führt ein wichtiges Konzept ein: den **Entry Workflow**.
Nextflow verwendet einen unbenannten `workflow {}`-Block als Einstiegspunkt, wenn du ein Skript direkt ausführst.
`GREETING_WORKFLOW` ist benannt, daher weiß Nextflow nicht, wie es ihn eigenständig ausführen soll.

Das ist beabsichtigt – kombinierbare Workflows sind dafür gedacht, von einem Entry Workflow aufgerufen zu werden, nicht direkt ausgeführt zu werden.
Die Lösung ist ein Entry Workflow in `main.nf`, der `GREETING_WORKFLOW` importiert und aufruft.

### 1.3. Den Haupt-Workflow aktualisieren und testen

Jetzt aktualisieren wir den Haupt-Workflow, um den Greeting Workflow aufzurufen.

#### 1.3.1. Den Greeting Workflow einbinden und aufrufen

Füge die `include`-Anweisung hinzu, aktualisiere den Workflow-Rumpf, um `GREETING_WORKFLOW` aufzurufen, und ersetze den `channel.empty()`-Platzhalter in `publish:`:

=== "Danach"

    ```groovy title="main.nf" linenums="1" hl_lines="1 7 8 11"
    include { GREETING_WORKFLOW } from './workflows/greeting'

    workflow {
        main:
        names = channel.of('Alice', 'Bob', 'Charlie')

        // Den Greeting Workflow ausführen
        GREETING_WORKFLOW(names)

        publish:
        greetings = GREETING_WORKFLOW.out.greetings
    }
    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="1"
    workflow {
        main:
        names = channel.of('Alice', 'Bob', 'Charlie')

        publish:
        greetings = channel.empty()
    }
    ```

Der Entry Workflow bleibt unbenannt, damit Nextflow ihn als Pipeline-Einstiegspunkt verwendet.

#### 1.3.2. Den output-Block aktualisieren

Füge eine `path`-Direktive hinzu, um veröffentlichte Begrüßungen in ein `greetings/`-Unterverzeichnis zu leiten:

=== "Danach"

    ```groovy title="main.nf" linenums="14" hl_lines="3"
    output {
        greetings {
            path 'greetings'
        }
    }
    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="14" hl_lines="2 3"
    output {
        greetings {
        }
    }
    ```

#### 1.3.3. Den Workflow ausführen

Führe den Workflow aus, um zu testen, ob er funktioniert:

```bash
nextflow run main.nf
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W  ~  version 24.10.0
    Launching `main.nf` [goofy_mayer] DSL2 - revision: 543f8742fe
    executor >  local (9)
    [05/3cc752] process > GREETING_WORKFLOW:VALIDATE_NAME (validating Char... [100%] 3 of 3 ✔
    [b1/b56ecf] process > GREETING_WORKFLOW:SAY_HELLO (greeting Charlie)      [100%] 3 of 3 ✔
    [ea/342168] process > GREETING_WORKFLOW:TIMESTAMP_GREETING (adding tim... [100%] 3 of 3 ✔
    ```

??? abstract "Verzeichnisinhalt"

    ```console
    results/
    └── greetings
        ├── Alice-output.txt
        ├── Bob-output.txt
        └── Charlie-output.txt
    ```

??? abstract "Dateiinhalt"

    ```console title="results/greetings/Alice-output.txt"
    Hello, Alice!
    ```

Die Begrüßungsdateien werden in `results/greetings/` veröffentlicht.
Der Haupt-Workflow ruft `GREETING_WORKFLOW` auf und leitet seine Ausgabe direkt an den `publish:`-Abschnitt weiter.

### Fazit

In diesem Abschnitt hast du mehrere wichtige Konzepte kennengelernt:

- **Benannte Workflows**: Einen benannten Workflow (`GREETING_WORKFLOW`) erstellen, der importiert und wiederverwendet werden kann
- **Workflow-Schnittstellen**: Klare Eingaben mit `take:` und Ausgaben mit `emit:` definieren, um einen kombinierbaren Workflow zu erstellen
- **Entry Points**: Verstehen, dass Nextflow einen unbenannten Entry Workflow benötigt, um ein Skript auszuführen
- **Workflow-Komposition**: Einen benannten Workflow innerhalb eines anderen Workflows importieren und verwenden
- **Workflow-Namespaces**: Auf Workflow-Ausgaben über den `.out`-Namespace zugreifen (`GREETING_WORKFLOW.out.greetings`)

Du hast jetzt einen funktionierenden Greeting Workflow, der:

- Einen Kanal mit Namen als Eingabe entgegennimmt
- Jeden Namen validiert
- Eine Begrüßung für jeden gültigen Namen erstellt
- Den Begrüßungen Zeitstempel hinzufügt
- Sowohl ursprüngliche als auch Begrüßungen mit Zeitstempel als Ausgaben bereitstellt

Dieser modulare Ansatz ermöglicht es dir, den Greeting Workflow unabhängig zu testen oder ihn als Komponente in größeren Pipelines zu verwenden.

---

## 2. Den Transform Workflow zur Pipeline hinzufügen

Der Transform Workflow wendet Texttransformationen auf die Begrüßungen mit Zeitstempel an.

### 2.1. Den Workflow ansehen und ausführen

Öffne `workflows/transform.nf` und schau dir den Code an:

```groovy title="workflows/transform.nf" linenums="1"
include { SAY_HELLO_UPPER } from '../modules/say_hello_upper'
include { REVERSE_TEXT } from '../modules/reverse_text'

workflow {
    main:
    input_ch = channel.fromPath('results/timestamped_*.txt')

    // Transformationen der Reihe nach anwenden
    upper_ch = SAY_HELLO_UPPER(input_ch)
    reversed_ch = REVERSE_TEXT(upper_ch)

    publish:
    upper = upper_ch
    reversed = reversed_ch
}

output {
    upper {
    }
    reversed {
    }
}
```

Dieser eigenständige Workflow liest Begrüßungsdateien mit Zeitstempel aus dem `results/`-Verzeichnis, das von `greeting.nf` erzeugt wurde, wandelt sie in Großbuchstaben um und kehrt dann den Text um.

Führe ihn aus, um zu überprüfen, ob er mit den Begrüßungsergebnissen aus Abschnitt 1.1 funktioniert:

```bash
nextflow run workflows/transform.nf
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W  ~  version 24.10.0
    Launching `workflows/transform.nf` [blissful_curie] DSL2 - revision: 4e7b1c9f02
    executor >  local (6)
    [3e/a14c29] process > SAY_HELLO_UPPER (converting t... [100%] 3 of 3 ✔
    [c8/51b9e3] process > REVERSE_TEXT (reversing UPPER... [100%] 3 of 3 ✔
    ```

Um ihn mit `GREETING_WORKFLOW` kombinierbar zu machen, gelten die gleichen drei Änderungen aus Abschnitt 1.2.

### 2.2. Ihn kombinierbar machen

Wende die gleichen drei Änderungen wie in Abschnitt 1.2 an: Benenne den Workflow, ersetze die fest kodierte Eingabe durch `take:`, und ersetze `publish:`/`output {}` durch `emit:`.

Die fertige Datei sollte so aussehen:

```groovy title="workflows/transform.nf" linenums="1" hl_lines="4 5 6 8 13 14 15"
include { SAY_HELLO_UPPER } from '../modules/say_hello_upper'
include { REVERSE_TEXT } from '../modules/reverse_text'

workflow TRANSFORM_WORKFLOW {
    take:
    input_ch // Eingabekanal mit Nachrichten

    main:
    // Transformationen der Reihe nach anwenden
    upper_ch = SAY_HELLO_UPPER(input_ch)
    reversed_ch = REVERSE_TEXT(upper_ch)

    emit:
    upper = upper_ch // Begrüßungen in Großbuchstaben
    reversed = reversed_ch // Umgekehrte Begrüßungen in Großbuchstaben
}
```

Der Transform Workflow ist jetzt kombinierbar und bereit, in den Haupt-Workflow importiert zu werden.

### 2.3. Den Haupt-Workflow aktualisieren und testen

Jetzt aktualisieren wir den Haupt-Workflow, um den Transform Workflow aufzurufen.

#### 2.3.1. Den Transform Workflow einbinden und aufrufen

Füge die include-Anweisung, einen Aufruf von `TRANSFORM_WORKFLOW` mit den Begrüßungen mit Zeitstempel und die zwei neuen `publish:`-Einträge hinzu:

=== "Danach"

    ```groovy title="main.nf" linenums="1" hl_lines="2 11 12 16 17"
    include { GREETING_WORKFLOW } from './workflows/greeting'
    include { TRANSFORM_WORKFLOW } from './workflows/transform'

    workflow {
        main:
        names = channel.of('Alice', 'Bob', 'Charlie')

        // Den Greeting Workflow ausführen
        GREETING_WORKFLOW(names)

        // Den Transform Workflow ausführen
        TRANSFORM_WORKFLOW(GREETING_WORKFLOW.out.timestamped)

        publish:
        greetings = GREETING_WORKFLOW.out.greetings
        upper = TRANSFORM_WORKFLOW.out.upper
        reversed = TRANSFORM_WORKFLOW.out.reversed
    }
    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="1"
    include { GREETING_WORKFLOW } from './workflows/greeting'

    workflow {
        main:
        names = channel.of('Alice', 'Bob', 'Charlie')

        // Den Greeting Workflow ausführen
        GREETING_WORKFLOW(names)

        publish:
        greetings = GREETING_WORKFLOW.out.greetings
    }
    ```

Damit wird der Transform Workflow auf die Begrüßungen mit Zeitstempel angewendet.

#### 2.3.2. Den output-Block aktualisieren

Füge `upper`- und `reversed`-Einträge zum `output {}`-Block hinzu, jeweils mit einer `path`-Direktive für das jeweilige Unterverzeichnis:

=== "Danach"

    ```groovy title="main.nf" linenums="20" hl_lines="5 6 7 8 9 10"
    output {
        greetings {
            path 'greetings'
        }
        upper {
            path 'upper'
        }
        reversed {
            path 'reversed'
        }
    }
    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="20" hl_lines="2 3 4 5"
    output {
        greetings {
            path 'greetings'
        }
    }
    ```

Damit werden die finalen Ausgaben in die entsprechenden Verzeichnisse veröffentlicht.

#### 2.3.3. Die vollständige Pipeline ausführen

Führe die Pipeline aus, um zu testen, ob alles funktioniert:

```bash
nextflow run main.nf
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W  ~  version 24.10.0
    Launching `main.nf` [sick_kimura] DSL2 - revision: 8dc45fc6a8
    executor >  local (15)
    [83/1b51f4] process > GREETING_WORKFLOW:VALIDATE_NAME (validating Alice)  [100%] 3 of 3 ✔
    [68/556150] process > GREETING_WORKFLOW:SAY_HELLO (greeting Alice)        [100%] 3 of 3 ✔
    [de/511abd] process > GREETING_WORKFLOW:TIMESTAMP_GREETING (adding tim... [100%] 3 of 3 ✔
    [cd/e6a7e0] process > TRANSFORM_WORKFLOW:SAY_HELLO_UPPER (converting t... [100%] 3 of 3 ✔
    [f0/74ba4a] process > TRANSFORM_WORKFLOW:REVERSE_TEXT (reversing UPPER... [100%] 3 of 3 ✔
    ```

??? abstract "Verzeichnisinhalt"

    ```console
    results/
    ├── greetings
    │   ├── Alice-output.txt
    │   ├── Bob-output.txt
    │   └── Charlie-output.txt
    ├── reversed
    │   ├── REVERSED-UPPER-timestamped_Alice-output.txt
    │   ├── REVERSED-UPPER-timestamped_Bob-output.txt
    │   └── REVERSED-UPPER-timestamped_Charlie-output.txt
    └── upper
        ├── UPPER-timestamped_Alice-output.txt
        ├── UPPER-timestamped_Bob-output.txt
        └── UPPER-timestamped_Charlie-output.txt
    ```

??? abstract "Dateiinhalt"

    ```console title="results/reversed/REVERSED-UPPER-timestamped_Alice-output.txt"
    !ECILA ,OLLEH ]04:50:71 60-30-5202[
    ```

Die Pipeline funktioniert von Anfang bis Ende: Die Begrüßung wurde in Großbuchstaben umgewandelt und umgekehrt.

### Fazit

Du solltest jetzt eine vollständige Pipeline haben, die:

- Namen durch den Greeting Workflow verarbeitet
- Die Begrüßungen mit Zeitstempel in den Transform Workflow einspeist
- Sowohl Großbuchstaben- als auch umgekehrte Versionen der Begrüßungen erzeugt

---

## Zusammenfassung

In diesem Side Quest haben wir das mächtige Konzept der Workflow-Komposition in Nextflow erkundet, das es uns ermöglicht, komplexe Pipelines aus kleineren, wiederverwendbaren Komponenten aufzubauen.

Dieser modulare Ansatz bietet gegenüber monolithischen Pipelines mehrere Vorteile:

- Jeder Workflow kann unabhängig entwickelt, getestet und debuggt werden
- Workflows können über verschiedene Pipelines hinweg wiederverwendet werden
- Die Gesamtstruktur der Pipeline wird lesbarer und wartbarer
- Änderungen an einem Workflow beeinflussen nicht zwangsläufig andere, solange die Schnittstellen konsistent bleiben
- Entry Points können so konfiguriert werden, dass sie verschiedene Teile deiner Pipeline nach Bedarf ausführen

_Es ist jedoch wichtig zu beachten, dass das Aufrufen von Workflows zwar ein bisschen wie das Aufrufen von Prozessen ist, aber nicht dasselbe. Du kannst einen Workflow beispielsweise nicht N-mal ausführen, indem du ihn mit einem Kanal der Größe N aufrufst – du müsstest einen Kanal der Größe N an den Workflow übergeben und intern iterieren._

Wenn du diese Techniken in deiner eigenen Arbeit anwendest, kannst du ausgefeiltere Nextflow-Pipelines erstellen, die komplexe Datenverarbeitungsaufgaben bewältigen und dabei wartbar und skalierbar bleiben.

### Wichtige Muster

1.  **Workflow-Struktur**: Wir haben klare Eingaben und Ausgaben für jeden Workflow mit der `take:`- und `emit:`-Syntax definiert, gut definierte Schnittstellen zwischen Komponenten erstellt und die Workflow-Logik im `main:`-Block eingebettet.

    ```groovy
    workflow EXAMPLE_WORKFLOW {
        take:
            // Eingabekanäle werden hier deklariert
            input_ch

        main:
            // Workflow-Logik kommt hier hin
            // Hier werden Prozesse aufgerufen und Kanäle verarbeitet
            result_ch = SOME_PROCESS(input_ch)

        emit:
            // Ausgabekanäle werden hier deklariert
            output_ch = result_ch
    }
    ```

2.  **Workflow-Importe:** Wir haben zwei unabhängige Workflow-Module erstellt und sie mit include-Anweisungen in eine Haupt-Pipeline importiert.

    - Einen einzelnen Workflow importieren

    ```groovy
    include { WORKFLOW_NAME } from './path/to/workflow'
    ```

    - Mehrere Workflows importieren

    ```groovy
    include { WORKFLOW_A; WORKFLOW_B } from './path/to/workflows'
    ```

    - Mit Alias importieren, um Namenskonflikte zu vermeiden

    ```groovy
    include { WORKFLOW_A as WORKFLOW_A_ALIAS } from './path/to/workflow'
    ```

3.  **Entry Points**: Nextflow benötigt einen unbenannten Entry Workflow, um zu wissen, wo die Ausführung beginnen soll. Dieser Entry Workflow ruft deine benannten Workflows auf.

    - Unbenannter Workflow (Entry Point)

    ```groovy
    workflow {
        // Dies ist der Einstiegspunkt, wenn das Skript ausgeführt wird
        NAMED_WORKFLOW(input_ch)
    }
    ```

    - Benannter Workflow (wird vom Entry Workflow aufgerufen)

    ```groovy
    workflow NAMED_WORKFLOW {
        // Muss vom Entry Workflow aufgerufen werden
    }
    ```

4.  **Datenflussverwaltung:** Wir haben gelernt, wie man auf Workflow-Ausgaben über die Namespace-Notation zugreift (`WORKFLOW_NAME.out.channel_name`) und sie an andere Workflows weitergibt.

    ```nextflow
    WORKFLOW_A(input_ch)
    WORKFLOW_B(WORKFLOW_A.out.some_channel)
    ```

### Weitere Ressourcen

- [Nextflow Workflow-Dokumentation](https://www.nextflow.io/docs/latest/workflow.html)
- [Referenz für Kanal-Operatoren](https://www.nextflow.io/docs/latest/operator.html)
- [Dokumentation zur Error Strategy](https://www.nextflow.io/docs/latest/process.html#errorstrategy)

---

## Wie geht es weiter?

Kehre zum [Menü der Side Quests](../index.md) zurück oder klicke auf die Schaltfläche unten rechts auf der Seite, um zum nächsten Thema in der Liste zu wechseln.
