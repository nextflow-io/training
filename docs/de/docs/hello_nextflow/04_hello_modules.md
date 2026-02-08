# Teil 4: Hello Modules

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/43Ot-f0iOME?si=y8lAedhEHWaTV4zd&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1&amp;cc_lang_pref=de" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } Siehe [die gesamte Playlist](https://youtube.com/playlist?list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&si=eF7cLR62goy-lc6n) auf dem Nextflow YouTube-Kanal.

:green_book: Das Videotranskript ist [hier](./transcripts/04_hello_modules.md) verfügbar.
///

Dieser Abschnitt behandelt, wie du deinen Workflow-Code organisierst, um die Entwicklung und Wartung deiner Pipeline effizienter und nachhaltiger zu gestalten.
Konkret werden wir demonstrieren, wie man [**Module**](https://nextflow.io/docs/latest/module.html) verwendet.

In Nextflow ist ein **Modul** eine eigenständige Codedatei, die oft eine einzelne process-Definition enthält.
Um ein Modul in einem Workflow zu verwenden, fügst du einfach eine einzeilige `include`-Anweisung zu deiner Workflow-Codedatei hinzu; dann kannst du den process genauso in den Workflow integrieren, wie du es normalerweise tun würdest.
Das ermöglicht es, process-Definitionen in mehreren Workflows wiederzuverwenden, ohne mehrere Kopien des Codes zu produzieren.

Als wir anfingen, unseren Workflow zu entwickeln, haben wir alles in einer einzigen Codedatei geschrieben.
Jetzt werden wir die processes in einzelne Module auslagern.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/modules.svg"
</figure>

Dies wird unseren Code besser teilbar, flexibler und wartbarer machen.

??? info "Wie du von diesem Abschnitt aus beginnen kannst"

    Dieser Abschnitt des Kurses setzt voraus, dass du die Teile 1-3 des [Hello Nextflow](./index.md)-Kurses abgeschlossen hast, aber wenn du mit den dort behandelten Grundlagen vertraut bist, kannst du von hier aus starten, ohne etwas Besonderes tun zu müssen.

---

## 0. Aufwärmübung: `hello-modules.nf` ausführen

Wir werden das Workflow-Skript `hello-modules.nf` als Ausgangspunkt verwenden.
Es entspricht dem Skript, das durch das Durcharbeiten von Teil 3 dieses Trainingskurses erstellt wurde, außer dass wir die Ausgabeziele geändert haben:

```groovy title="hello-modules.nf" linenums="37" hl_lines="3 7 11 15"
output {
    first_output {
        path 'hello_modules'
        mode 'copy'
    }
    uppercased {
        path 'hello_modules'
        mode 'copy'
    }
    collected {
        path 'hello_modules'
        mode 'copy'
    }
    batch_report {
        path 'hello_modules'
        mode 'copy'
    }
}
```

Um sicherzustellen, dass alles funktioniert, führe das Skript einmal aus, bevor du Änderungen vornimmst:

```bash
nextflow run hello-modules.nf
```

??? success "Befehlsausgabe"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-modules.nf` [hopeful_avogadro] DSL2 - revision: b09af1237d

    executor >  local (7)
    [0f/8795c9] sayHello (3)       [100%] 3 of 3 ✔
    [6a/eb2510] convertToUpper (3) [100%] 3 of 3 ✔
    [af/479117] collectGreetings   [100%] 1 of 1 ✔
    ```

Wie zuvor findest du die Ausgabedateien in dem im `output`-Block angegebenen Verzeichnis (hier `results/hello_modules/`).

??? abstract "Verzeichnisinhalte"

    ```console
    results/hello_modules/
    ├── Bonjour-output.txt
    ├── COLLECTED-batch-output.txt
    ├── Hello-output.txt
    ├── Holà-output.txt
    ├── batch-report.txt
    ├── UPPER-Bonjour-output.txt
    ├── UPPER-Hello-output.txt
    └── UPPER-Holà-output.txt
    ```

Wenn das bei dir funktioniert hat, bist du bereit zu lernen, wie du deinen Workflow-Code modularisierst.

---

## 1. Ein Verzeichnis zum Speichern von Modulen erstellen

Es ist eine bewährte Praxis, deine Module in einem bestimmten Verzeichnis zu speichern.
Du kannst dieses Verzeichnis beliebig nennen, aber die Konvention ist, es `modules/` zu nennen.

```bash
mkdir modules
```

---

## 2. Ein Modul für `sayHello()` erstellen

In seiner einfachsten Form ist das Umwandeln eines bestehenden process in ein Modul kaum mehr als eine Kopier-und-Einfüge-Operation.
Wir werden einen Dateistub für das Modul erstellen, den relevanten Code hinüberkopieren und ihn dann aus der Haupt-Workflow-Datei löschen.

Dann müssen wir nur noch eine `include`-Anweisung hinzufügen, damit Nextflow weiß, dass es den relevanten Code zur Laufzeit einziehen soll.

### 2.1. Einen Dateistub für das neue Modul erstellen

Lass uns eine leere Datei für das Modul namens `sayHello.nf` erstellen.

```bash
touch modules/sayHello.nf
```

Dies gibt uns einen Platz, um den process-Code zu platzieren.

### 2.2. Den `sayHello`-process-Code in die Moduldatei verschieben

Kopiere die gesamte process-Definition von der Workflow-Datei in die Moduldatei.

```groovy title="modules/sayHello.nf" linenums="1"
/*
 * Verwende echo, um 'Hello World!' in eine Datei zu schreiben
 */
process sayHello {

    input:
    val greeting

    output:
    path "${greeting}-output.txt"

    script:
    """
    echo '${greeting}' > '${greeting}-output.txt'
    """
}
```

Sobald das erledigt ist, lösche die process-Definition aus der Workflow-Datei.

### 2.3. Eine include-Deklaration vor dem workflow-Block hinzufügen

Die Syntax für das Einbinden eines process aus einem Modul ist ziemlich unkompliziert:

```groovy title="Syntax: include declaration"
include { <PROCESS_NAME> } from '<path_to_module>'
```

Lass uns das oberhalb des `params`-Blocks einfügen und es entsprechend ausfüllen.

=== "Nachher"

    ```groovy title="hello-modules.nf" linenums="44" hl_lines="1 2"
    // Include modules
    include { sayHello } from './modules/sayHello.nf'

    /*
    * Pipeline-Parameter
    */
    params {
        greeting: Path = 'data/greetings.csv'
        batch: String = 'batch'
    }
    ```

=== "Vorher"

    ```groovy title="hello-modules.nf" linenums="44"
    /*
    * Pipeline-Parameter
    */
    params {
        greeting: Path = 'data/greetings.csv'
        batch: String = 'batch'
    }
    ```

Du siehst, wir haben den Prozessnamen `sayHello` und den Pfad zu der Datei, die den Modulcode enthält, `./modules/sayHello.nf`, ausgefüllt.

### 2.4. Den Workflow ausführen

Wir führen den Workflow mit im Wesentlichen demselben Code und denselben Eingaben wie zuvor aus, also lass ihn uns mit dem `-resume`-Flag ausführen und sehen, was passiert.

```bash
nextflow run hello-modules.nf -resume
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-modules.nf` [romantic_poisson] DSL2 - revision: 96edfa9ad3

    [f6/cc0107] sayHello (1)       | 3 of 3, cached: 3 ✔
    [3c/4058ba] convertToUpper (2) | 3 of 3, cached: 3 ✔
    [1a/bc5901] collectGreetings   | 1 of 1, cached: 1 ✔
    ```

Dies sollte sehr schnell laufen, weil alles im Cache ist.
Du kannst gerne die veröffentlichten Ausgaben überprüfen.

Nextflow hat erkannt, dass es immer noch dieselbe Arbeit zu erledigen gibt, auch wenn der Code auf mehrere Dateien aufgeteilt ist.

### Fazit

Du weißt, wie du einen process in ein lokales Modul extrahierst und weißt, dass dies die Wiederaufnahmefähigkeit des Workflows nicht beeinträchtigt.

### Wie geht es weiter?

Übe das Erstellen weiterer Module.
Sobald du eines gemacht hast, kannst du eine Million weitere machen...
Aber lass uns vorerst nur noch zwei weitere machen.

---

## 3. Den `convertToUpper()`-process modularisieren

### 3.1. Einen Dateistub für das neue Modul erstellen

Erstelle eine leere Datei für das Modul namens `convertToUpper.nf`.

```bash
touch modules/convertToUpper.nf
```

### 3.2. Den `convertToUpper`-process-Code in die Moduldatei verschieben

Kopiere die gesamte process-Definition von der Workflow-Datei in die Moduldatei.

```groovy title="modules/convertToUpper.nf" linenums="1"
/*
 * Verwende ein Textersetzungstool, um die Begrüßung in Großbuchstaben umzuwandeln
 */
process convertToUpper {

    input:
    path input_file

    output:
    path "UPPER-${input_file}"

    script:
    """
    cat '${input_file}' | tr '[a-z]' '[A-Z]' > 'UPPER-${input_file}'
    """
}
```

Sobald das erledigt ist, lösche die process-Definition aus der Workflow-Datei.

### 3.3. Eine include-Deklaration vor dem `params`-Block hinzufügen

Füge die include-Deklaration oberhalb des `params`-Blocks ein und fülle sie entsprechend aus.

=== "Nachher"

    ```groovy title="hello-modules.nf" linenums="23" hl_lines="3"
    // Include modules
    include { sayHello } from './modules/sayHello.nf'
    include { convertToUpper } from './modules/convertToUpper.nf'

    /*
    * Pipeline-Parameter
    */
    params {
        greeting: Path = 'data/greetings.csv'
        batch: String = 'batch'
    }
    ```

=== "Vorher"

    ```groovy title="hello-modules.nf" linenums="23"
    // Include modules
    include { sayHello } from './modules/sayHello.nf'

    /*
    * Pipeline-Parameter
    */
    params {
        greeting: Path = 'data/greetings.csv'
        batch: String = 'batch'
    }
    ```

Dies sollte jetzt sehr vertraut aussehen.

### 3.4. Den Workflow erneut ausführen

Führe dies mit dem `-resume`-Flag aus.

```bash
nextflow run hello-modules.nf -resume
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-modules.nf` [nauseous_heisenberg] DSL2 - revision: a04a9f2da0

    [c9/763d42] sayHello (3)       | 3 of 3, cached: 3 ✔
    [60/bc6831] convertToUpper (3) | 3 of 3, cached: 3 ✔
    [1a/bc5901] collectGreetings   | 1 of 1, cached: 1 ✔
    ```

Dies sollte immer noch dieselbe Ausgabe wie zuvor produzieren.

Zwei erledigt, noch eines übrig!

---

## 4. Den `collectGreetings()`-process modularisieren

### 4.1. Einen Dateistub für das neue Modul erstellen

Erstelle eine leere Datei für das Modul namens `collectGreetings.nf`.

```bash
touch modules/collectGreetings.nf
```

### 4.2. Den `collectGreetings`-process-Code in die Moduldatei verschieben

Kopiere die gesamte process-Definition von der Workflow-Datei in die Moduldatei.

```groovy title="modules/collectGreetings.nf" linenums="1"
/*
 * Großbuchstaben-Begrüßungen in einer einzelnen Ausgabedatei sammeln
 */
process collectGreetings {

    input:
    path input_files
    val batch_name

    output:
    path "COLLECTED-${batch_name}-output.txt", emit: outfile
    path "${batch_name}-report.txt", emit: report

    script:
    count_greetings = input_files.size()
    """
    cat ${input_files} > 'COLLECTED-${batch_name}-output.txt'
    echo 'There were ${count_greetings} greetings in this batch.' > '${batch_name}-report.txt'
    """
}
```

Sobald das erledigt ist, lösche die process-Definition aus der Workflow-Datei.

### 4.3. Eine include-Deklaration vor dem `params`-Block hinzufügen

Füge die include-Deklaration oberhalb des `params`-Blocks ein und fülle sie entsprechend aus.

=== "Nachher"

    ```groovy title="hello-modules.nf" linenums="3" hl_lines="4"
    // Include modules
    include { sayHello } from './modules/sayHello.nf'
    include { convertToUpper } from './modules/convertToUpper.nf'
    include { collectGreetings } from './modules/collectGreetings.nf'

    /*
    * Pipeline-Parameter
    */
    params {
        greeting: Path = 'data/greetings.csv'
        batch: String = 'batch'
    }
    ```

=== "Vorher"

    ```groovy title="hello-modules.nf" linenums="3"
    // Include modules
    include { sayHello } from './modules/sayHello.nf'
    include { convertToUpper } from './modules/convertToUpper.nf'

    /*
    * Pipeline-Parameter
    */
    params {
        greeting: Path = 'data/greetings.csv'
        batch: String = 'batch'
    }
    ```

Das letzte!

### 4.4. Den Workflow ausführen

Führe dies mit dem `-resume`-Flag aus.

```bash
nextflow run hello-modules.nf -resume
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-modules.nf` [friendly_coulomb] DSL2 - revision: 7aa2b9bc0f

    [f6/cc0107] sayHello (1)       | 3 of 3, cached: 3 ✔
    [3c/4058ba] convertToUpper (2) | 3 of 3, cached: 3 ✔
    [1a/bc5901] collectGreetings   | 1 of 1, cached: 1 ✔
    ```

Dies sollte immer noch dieselbe Ausgabe wie zuvor produzieren.

### Fazit

Du weißt, wie du mehrere processes in einem Workflow modularisierst.

Herzlichen Glückwunsch, du hast all diese Arbeit geleistet und absolut nichts hat sich daran geändert, wie die Pipeline funktioniert!

Scherz beiseite, jetzt ist dein Code modularer, und wenn du dich entscheidest, eine andere Pipeline zu schreiben, die einen dieser processes aufruft, musst du nur eine kurze `include`-Anweisung eingeben, um das relevante Modul zu verwenden.
Das ist besser als den Code zu kopieren und einzufügen, denn wenn du dich später entscheidest, das Modul zu verbessern, werden alle deine Pipelines die Verbesserungen erben.

### Wie geht es weiter?

Mach eine kurze Pause, wenn du möchtest.

Wenn du bereit bist, gehe zu [**Teil 5: Hello Container**](./05_hello_containers.md), um zu lernen, wie du Container verwendest, um Softwareabhängigkeiten bequemer und reproduzierbarer zu verwalten.

---

## Quiz

<quiz>
Was ist ein Modul in Nextflow?
- [ ] Eine Konfigurationsdatei
- [x] Eine eigenständige Datei, die process-Definitionen enthalten kann
- [ ] Eine Workflow-Definition
- [ ] Ein channel-Operator

Mehr erfahren: [2. Ein Modul für `sayHello()` erstellen](#2-create-a-module-for-sayhello)
</quiz>

<quiz>
Welche Konvention wird typischerweise für das Speichern von Moduldateien verwendet?
- [ ] Im selben Verzeichnis wie der Workflow
- [ ] In einem `bin/`-Verzeichnis
- [x] In einem `modules/`-Verzeichnis
- [ ] In einem `lib/`-Verzeichnis

Mehr erfahren: [1. Ein Verzeichnis zum Speichern von Modulen erstellen](#1-create-a-directory-to-store-modules)
</quiz>

<quiz>
Was ist die korrekte Syntax, um ein Modul zu verwenden?

- [ ] `#!groovy import { SAYHELLO } from './modules/sayhello.nf'`
- [ ] `#!groovy require { SAYHELLO } from './modules/sayhello.nf'`
- [x] `#!groovy include { SAYHELLO } from './modules/sayhello.nf'`
- [ ] `#!groovy load { SAYHELLO } from './modules/sayhello.nf'`

Mehr erfahren: [2.3. Eine include-Deklaration hinzufügen](#23-add-an-include-declaration-before-the-workflow-block)
</quiz>

<quiz>
Was passiert mit der `-resume`-Funktionalität bei der Verwendung von Modulen?
- [ ] Sie funktioniert nicht mehr
- [ ] Sie erfordert zusätzliche Konfiguration
- [x] Sie funktioniert genauso wie zuvor
- [ ] Sie funktioniert nur für lokale Module
</quiz>

<quiz>
Was sind die Vorteile der Verwendung von Modulen? (Wähle alle zutreffenden aus)
- [x] Wiederverwendbarkeit von Code über Workflows hinweg
- [x] Einfachere Wartung
- [x] Bessere Organisation von Workflow-Code
- [ ] Schnellere Ausführungsgeschwindigkeit
</quiz>
