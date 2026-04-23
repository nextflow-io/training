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

- Das Tutorial [Hello Nextflow](../hello_nextflow/README.md) oder einen gleichwertigen Einsteigerkurs abgeschlossen haben.
- Mit grundlegenden Nextflow-Konzepten und -Mechanismen vertraut sein (Prozesse, Kanäle, Operatoren, Module)

---

## 0. Erste Schritte

#### Öffne den Training-Codespace

Falls du das noch nicht getan hast, öffne die Trainingsumgebung wie in der [Umgebung einrichten](../envsetup/index.md) beschrieben.

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

#### Schau dir die Materialien an

Du findest ein `modules`-Verzeichnis mit mehreren Prozessdefinitionen, die auf dem aufbauen, was du in 'Hello Nextflow' gelernt hast:

```console title="Directory contents"
modules/
├── say_hello.nf             # Creates a greeting (from Hello Nextflow)
├── say_hello_upper.nf       # Converts to uppercase (from Hello Nextflow)
├── timestamp_greeting.nf    # Adds timestamps to greetings
├── validate_name.nf         # Validates input names
└── reverse_text.nf          # Reverses text content
```

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

## 1. Den Greeting Workflow erstellen

Beginnen wir mit der Erstellung eines Workflows, der Namen validiert und Begrüßungen mit Zeitstempel generiert.

### 1.1. Die Workflow-Struktur erstellen

```bash title="Create workflow directory and file"
mkdir -p workflows
touch workflows/greeting.nf
```

### 1.2. Den ersten (Sub-)Workflow-Code hinzufügen

Füge diesen Code zu `workflows/greeting.nf` hinzu:

```groovy title="workflows/greeting.nf" linenums="1"
include { VALIDATE_NAME } from '../modules/validate_name'
include { SAY_HELLO } from '../modules/say_hello'
include { TIMESTAMP_GREETING } from '../modules/timestamp_greeting'

workflow {

    names_ch = channel.of('Alice', 'Bob', 'Charlie')

    // Prozesse verketten: validieren -> Begrüßung erstellen -> Zeitstempel hinzufügen
    validated_ch = VALIDATE_NAME(names_ch)
    greetings_ch = SAY_HELLO(validated_ch)
    timestamped_ch = TIMESTAMP_GREETING(greetings_ch)
}
```

Dies ist ein vollständiger Workflow mit einer ähnlichen Struktur wie die, die du im Tutorial 'Hello Nextflow' gesehen hast, und den wir unabhängig testen können. Lass uns das jetzt ausprobieren:

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

Das funktioniert wie erwartet, aber damit es kombinierbar wird, müssen wir einige Dinge ändern.

### 1.3. Den Workflow kombinierbar machen

Kombinierbare Workflows unterscheiden sich in einigen Punkten von denen, die du im Tutorial 'Hello Nextflow' gesehen hast:

- Der workflow-Block muss einen Namen haben
- Eingaben werden mit dem Schlüsselwort `take:` deklariert
- Der Workflow-Inhalt wird im `main:`-Block platziert
- Ausgaben werden mit dem Schlüsselwort `emit:` deklariert

Lass uns den Greeting Workflow aktualisieren, um dieser Struktur zu entsprechen. Ändere den Code wie folgt:

<!-- TODO: switch to before/after tabs -->

```groovy title="workflows/greeting.nf" linenums="1" hl_lines="6 7 9 15 16 17"
include { VALIDATE_NAME } from '../modules/validate_name'
include { SAY_HELLO } from '../modules/say_hello'
include { TIMESTAMP_GREETING } from '../modules/timestamp_greeting'

workflow GREETING_WORKFLOW {
    take:
        names_ch        // Eingabekanal mit Namen

    main:
        // Prozesse verketten: validieren -> Begrüßung erstellen -> Zeitstempel hinzufügen
        validated_ch = VALIDATE_NAME(names_ch)
        greetings_ch = SAY_HELLO(validated_ch)
        timestamped_ch = TIMESTAMP_GREETING(greetings_ch)

    emit:
        greetings = greetings_ch      // Ursprüngliche Begrüßungen
        timestamped = timestamped_ch  // Begrüßungen mit Zeitstempel
}
```

Du siehst, dass der Workflow jetzt einen Namen hat und einen `take:`- und `emit:`-Block besitzt – das sind die Verbindungen, die wir verwenden werden, um einen übergeordneten Workflow zusammenzusetzen.
Der Workflow-Inhalt wird auch im `main:`-Block platziert. Beachte außerdem, dass wir die Deklaration des Eingabekanals `names_ch` entfernt haben, da er jetzt als Argument an den Workflow übergeben wird.

Lass uns den Workflow erneut testen, um zu sehen, ob er wie erwartet funktioniert:

```bash
nextflow run workflows/greeting.nf
```

??? failure "Befehlsausgabe"

    ```console
    N E X T F L O W  ~  version 24.10.0
    Launching `workflows/greeting.nf` [high_brahmagupta] DSL2 - revision: 8f5857af25
    No entry workflow specified
    ```

Das weist dich auf ein weiteres neues Konzept hin: einen 'Entry Workflow'. Der Entry Workflow ist der Workflow, der aufgerufen wird, wenn du ein Nextflow-Skript ausführst. Standardmäßig verwendet Nextflow einen unbenannten Workflow als Entry Workflow, wenn einer vorhanden ist – und genau das hast du bisher gemacht, mit workflow-Blöcken, die so beginnen:

```groovy title="hello.nf" linenums="1"
workflow {
```

Unser Greeting Workflow hat jedoch keinen unbenannten Workflow, sondern einen benannten:

```groovy title="workflows/greeting.nf" linenums="1"
workflow GREETING_WORKFLOW {
```

Deshalb hat Nextflow einen Fehler ausgegeben und nicht das getan, was wir wollten.

Wir haben die `take:`/`emit:`-Syntax nicht hinzugefügt, um den Workflow direkt aufzurufen – wir haben es getan, damit wir ihn mit anderen Workflows kombinieren können. Die Lösung besteht darin, ein Hauptskript mit einem unbenannten Entry Workflow zu erstellen, der unseren benannten Workflow importiert und aufruft.

### 1.4. Den Haupt-Workflow erstellen und testen

Jetzt erstellen wir einen Haupt-Workflow, der den `greeting`-Workflow importiert und verwendet.

Erstelle `main.nf`:

```groovy title="main.nf" linenums="1"
include { GREETING_WORKFLOW } from './workflows/greeting'

workflow {
    names = channel.of('Alice', 'Bob', 'Charlie')
    GREETING_WORKFLOW(names)

    GREETING_WORKFLOW.out.greetings.view { "Original: $it" }
    GREETING_WORKFLOW.out.timestamped.view { "Timestamped: $it" }
}

```

Beachte, dass unser workflow-Einstiegspunkt in dieser Datei unbenannt ist – das liegt daran, dass wir ihn als Entry Workflow verwenden werden.

Führe das aus und schau dir die Ausgabe an:

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
    Original: /workspaces/training/side_quests/workflows_of_workflows/work/bb/c8aff3df0ebc15a4d7d35f736db44c/Alice-output.txt
    Original: /workspaces/training/side_quests/workflows_of_workflows/work/fb/fa877776e8a5d90b537b1bcd3b6f5b/Bob-output.txt
    Original: /workspaces/training/side_quests/workflows_of_workflows/work/b1/b56ecf938fda8bcbec211847c8f0be/Charlie-output.txt
    Timestamped: /workspaces/training/side_quests/workflows_of_workflows/work/06/877bc909f140bbf8223343450cea36/timestamped_Alice-output.txt
    Timestamped: /workspaces/training/side_quests/workflows_of_workflows/work/aa/bd31b71cdb745b7c155ca7f8837b8a/timestamped_Bob-output.txt
    Timestamped: /workspaces/training/side_quests/workflows_of_workflows/work/ea/342168d4ba04cc899a89c56cbfd9b0/timestamped_Charlie-output.txt
    ```

Es funktioniert! Wir haben den benannten Greeting Workflow in einen Haupt-Workflow mit einem unbenannten Entry-`workflow`-Block eingebettet. Der Haupt-Workflow verwendet den `GREETING_WORKFLOW` fast (aber nicht ganz) wie einen Prozess und übergibt den `names`-Kanal als Argument.

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

## 2. Den Transform Workflow hinzufügen

Jetzt erstellen wir einen Workflow, der Texttransformationen auf die Begrüßungen anwendet.

### 2.1. Die Workflow-Datei erstellen

```bash
touch workflows/transform.nf
```

### 2.2. Den Workflow-Code hinzufügen

Füge diesen Code zu `workflows/transform.nf` hinzu:

```groovy title="workflows/transform.nf" linenums="1"
include { SAY_HELLO_UPPER } from '../modules/say_hello_upper'
include { REVERSE_TEXT } from '../modules/reverse_text'

workflow TRANSFORM_WORKFLOW {
    take:
        input_ch         // Eingabekanal mit Nachrichten

    main:
        // Transformationen der Reihe nach anwenden
        upper_ch = SAY_HELLO_UPPER(input_ch)
        reversed_ch = REVERSE_TEXT(upper_ch)

    emit:
        upper = upper_ch        // Begrüßungen in Großbuchstaben
        reversed = reversed_ch  // Umgekehrte Begrüßungen in Großbuchstaben
}
```

Wir wiederholen die Erklärung der kombinierbaren Syntax hier nicht, aber beachte, dass der benannte Workflow wieder mit einem `take:`- und `emit:`-Block deklariert wird und der Workflow-Inhalt im `main:`-Block platziert ist.

### 2.3. Den Haupt-Workflow aktualisieren

Aktualisiere `main.nf`, um beide Workflows zu verwenden:

```groovy title="main.nf" linenums="1"
include { GREETING_WORKFLOW } from './workflows/greeting'
include { TRANSFORM_WORKFLOW } from './workflows/transform'

workflow {
    names = channel.of('Alice', 'Bob', 'Charlie')

    // Den Greeting Workflow ausführen
    GREETING_WORKFLOW(names)

    // Den Transform Workflow ausführen
    TRANSFORM_WORKFLOW(GREETING_WORKFLOW.out.timestamped)

    // Ergebnisse anzeigen
    TRANSFORM_WORKFLOW.out.upper.view { "Uppercase: $it" }
    TRANSFORM_WORKFLOW.out.reversed.view { "Reversed: $it" }
}
```

Die vollständige Pipeline ausführen:

```bash
nextflow run main.nf
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W  ~  version 24.10.0
    Launching `main.nf` [sick_kimura] DSL2 - revision: 8dc45fc6a8
    executor >  local (13)
    executor >  local (15)
    [83/1b51f4] process > GREETING_WORKFLOW:VALIDATE_NAME (validating Alice)  [100%] 3 of 3 ✔
    [68/556150] process > GREETING_WORKFLOW:SAY_HELLO (greeting Alice)        [100%] 3 of 3 ✔
    [de/511abd] process > GREETING_WORKFLOW:TIMESTAMP_GREETING (adding tim... [100%] 3 of 3 ✔
    [cd/e6a7e0] process > TRANSFORM_WORKFLOW:SAY_HELLO_UPPER (converting t... [100%] 3 of 3 ✔
    [f0/74ba4a] process > TRANSFORM_WORKFLOW:REVERSE_TEXT (reversing UPPER... [100%] 3 of 3 ✔
    Uppercase: /workspaces/training/side_quests/workflows_of_workflows/work/a0/d4f5df4d6344604498fa47a6084a11/UPPER-timestamped_Bob-output.txt
    Uppercase: /workspaces/training/side_quests/workflows_of_workflows/work/69/b5e37f6c79c2fd38adb75d0eca8f87/UPPER-timestamped_Charlie-output.txt
    Uppercase: /workspaces/training/side_quests/workflows_of_workflows/work/cd/e6a7e0b17e7d5a2f71bb8123cd53a7/UPPER-timestamped_Alice-output.txt
    Reversed: /workspaces/training/side_quests/workflows_of_workflows/work/7a/7a222f7957b35d1d121338566a24ac/REVERSED-UPPER-timestamped_Bob-output.txt
    Reversed: /workspaces/training/side_quests/workflows_of_workflows/work/46/8d19af62e33a5a6417c773496e0f90/REVERSED-UPPER-timestamped_Charlie-output.txt
    Reversed: /workspaces/training/side_quests/workflows_of_workflows/work/f0/74ba4a10d9ef5c82f829d1c154d0f6/REVERSED-UPPER-timestamped_Alice-output.txt
    ```

Wenn du dir eine dieser umgekehrten Dateien anschaust, siehst du, dass es die Großbuchstabenversion der Begrüßung umgekehrt ist:

```bash
cat /workspaces/training/side_quests/workflows_of_workflows/work/f0/74ba4a10d9ef5c82f829d1c154d0f6/REVERSED-UPPER-timestamped_Alice-output.txt
```

```console title="Reversed file content"
!ECILA ,OLLEH ]04:50:71 60-30-5202[
```

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

Wenn du diese Techniken in deiner eigenen Arbeit anwendest, kannst du ausgefeiltere Nextflow-Pipelines erstellen, die komplexe Bioinformatik-Aufgaben bewältigen und dabei wartbar und skalierbar bleiben.

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

Kehre zum [Menü der Side Quests](../) zurück oder klicke auf die Schaltfläche unten rechts auf der Seite, um zum nächsten Thema in der Liste zu wechseln.
