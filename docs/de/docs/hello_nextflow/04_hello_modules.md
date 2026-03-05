# Teil 4: Hello Modules

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/43Ot-f0iOME?si=y8lAedhEHWaTV4zd&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } Siehe [die gesamte Playlist](https://youtube.com/playlist?list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&si=eF7cLR62goy-lc6n) auf dem Nextflow YouTube-Kanal.

:green_book: Das Video-Transkript ist [hier](./transcripts/04_hello_modules.md) verfügbar.
///

Dieser Abschnitt behandelt, wie du deinen Workflow-Code organisierst, um die Entwicklung und Wartung deiner Pipeline effizienter und nachhaltiger zu gestalten.
Konkret zeigen wir dir, wie du [**Module**](https://nextflow.io/docs/latest/module.html) verwendest.

In Nextflow ist ein **Modul** eine eigenständige Code-Datei, die oft eine einzelne Prozess-Definition enthält.
Um ein Modul in einem Workflow zu verwenden, fügst du einfach eine einzeilige `include`-Anweisung zu deiner Workflow-Code-Datei hinzu; dann kannst du den Prozess auf die gleiche Weise in den Workflow integrieren wie normalerweise.
Das macht es möglich, Prozess-Definitionen in mehreren Workflows wiederzuverwenden, ohne mehrere Kopien des Codes zu erstellen.

Als wir mit der Entwicklung unseres Workflows begonnen haben, haben wir alles in einer einzigen Code-Datei geschrieben.
Jetzt werden wir die Prozesse in einzelne Module auslagern.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/modules.svg"
</figure>

Das macht unseren Code besser teilbar, flexibler und wartbarer.

??? info "Wie du von diesem Abschnitt aus beginnst"

    Dieser Abschnitt des Kurses setzt voraus, dass du die Teile 1-3 des [Hello Nextflow](./index.md)-Kurses abgeschlossen hast. Wenn du jedoch mit den dort behandelten Grundlagen vertraut bist, kannst du auch hier beginnen, ohne etwas Besonderes tun zu müssen.

---

## 0. Aufwärmen: Führe `hello-modules.nf` aus

Wir werden das Workflow-Skript `hello-modules.nf` als Ausgangspunkt verwenden.
Es entspricht dem Skript, das durch die Bearbeitung von Teil 3 dieses Trainingskurses erstellt wurde, außer dass wir die Ausgabeziele geändert haben:

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

Wie zuvor findest du die Ausgabedateien im Verzeichnis, das im `output`-Block angegeben ist (hier `results/hello_modules/`).

??? abstract "Verzeichnisinhalt"

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

## 1. Erstelle ein Verzeichnis zum Speichern von Modulen

Es ist Best Practice, deine Module in einem bestimmten Verzeichnis zu speichern.
Du kannst dieses Verzeichnis beliebig benennen, aber die Konvention ist, es `modules/` zu nennen.

```bash
mkdir modules
```

---

## 2. Erstelle ein Modul für `sayHello()`

In seiner einfachsten Form ist die Umwandlung eines bestehenden Prozesses in ein Modul kaum mehr als eine Kopier-Einfügen-Operation.
Wir werden eine Datei-Vorlage für das Modul erstellen, den relevanten Code hinüberkopieren und ihn dann aus der Haupt-Workflow-Datei löschen.

Dann müssen wir nur noch eine `include`-Anweisung hinzufügen, damit Nextflow weiß, dass es den relevanten Code zur Laufzeit einbinden soll.

### 2.1. Erstelle eine Datei-Vorlage für das neue Modul

Lass uns eine leere Datei für das Modul namens `sayHello.nf` erstellen.

```bash
touch modules/sayHello.nf
```

Das gibt uns einen Ort, an dem wir den Prozess-Code ablegen können.

### 2.2. Verschiebe den `sayHello`-Prozess-Code in die Modul-Datei

Kopiere die gesamte Prozess-Definition von der Workflow-Datei in die Modul-Datei.

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

Sobald das erledigt ist, lösche die Prozess-Definition aus der Workflow-Datei.

### 2.3. Füge eine include-Deklaration vor dem Workflow-Block hinzu

Die Syntax zum Einbinden eines Prozesses aus einem Modul ist ziemlich einfach:

```groovy title="Syntax: include declaration"
include { <PROCESS_NAME> } from '<path_to_module>'
```

Lass uns das oberhalb des `params`-Blocks einfügen und entsprechend ausfüllen.

=== "Danach"

    ```groovy title="hello-modules.nf" linenums="44" hl_lines="1 2"
    // Module einbinden
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

Du siehst, wir haben den Prozessnamen `sayHello` und den Pfad zur Datei mit dem Modul-Code `./modules/sayHello.nf` ausgefüllt.

### 2.4. Führe den Workflow aus

Wir führen den Workflow mit im Wesentlichen dem gleichen Code und den gleichen Eingaben wie zuvor aus, also lass uns mit dem `-resume`-Flag ausführen und sehen, was passiert.

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

Das sollte sehr schnell laufen, weil alles gecacht ist.
Du kannst gerne die veröffentlichten Ausgaben überprüfen.

Nextflow hat erkannt, dass es immer noch die gleiche Arbeit zu erledigen gibt, auch wenn der Code auf mehrere Dateien aufgeteilt ist.

### Fazit

Du weißt, wie du einen Prozess in ein lokales Modul extrahierst, und du weißt, dass dies die Wiederaufnahmefähigkeit des Workflows nicht beeinträchtigt.

### Wie geht es weiter?

Übe das Erstellen weiterer Module.
Wenn du eines gemacht hast, kannst du eine Million mehr machen...
Aber lass uns vorerst nur noch zwei weitere machen.

---

## 3. Modularisiere den `convertToUpper()`-Prozess

### 3.1. Erstelle eine Datei-Vorlage für das neue Modul

Erstelle eine leere Datei für das Modul namens `convertToUpper.nf`.

```bash
touch modules/convertToUpper.nf
```

### 3.2. Verschiebe den `convertToUpper`-Prozess-Code in die Modul-Datei

Kopiere die gesamte Prozess-Definition von der Workflow-Datei in die Modul-Datei.

```groovy title="modules/convertToUpper.nf" linenums="1"
/*
 * Verwende ein Text-Ersetzungstool, um die Begrüßung in Großbuchstaben umzuwandeln
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

Sobald das erledigt ist, lösche die Prozess-Definition aus der Workflow-Datei.

### 3.3. Füge eine include-Deklaration vor dem `params`-Block hinzu

Füge die include-Deklaration oberhalb des `params`-Blocks ein und fülle sie entsprechend aus.

=== "Danach"

    ```groovy title="hello-modules.nf" linenums="23" hl_lines="3"
    // Module einbinden
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
    // Module einbinden
    include { sayHello } from './modules/sayHello.nf'

    /*
    * Pipeline-Parameter
    */
    params {
        greeting: Path = 'data/greetings.csv'
        batch: String = 'batch'
    }
    ```

Das sollte dir langsam sehr vertraut vorkommen.

### 3.4. Führe den Workflow erneut aus

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

Das sollte immer noch die gleiche Ausgabe wie zuvor erzeugen.

Zwei geschafft, noch einer!

---

## 4. Modularisiere den `collectGreetings()`-Prozess

### 4.1. Erstelle eine Datei-Vorlage für das neue Modul

Erstelle eine leere Datei für das Modul namens `collectGreetings.nf`.

```bash
touch modules/collectGreetings.nf
```

### 4.2. Verschiebe den `collectGreetings`-Prozess-Code in die Modul-Datei

Kopiere die gesamte Prozess-Definition von der Workflow-Datei in die Modul-Datei.

```groovy title="modules/collectGreetings.nf" linenums="1"
/*
 * Sammle Großbuchstaben-Begrüßungen in einer einzigen Ausgabedatei
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

Sobald das erledigt ist, lösche die Prozess-Definition aus der Workflow-Datei.

### 4.3. Füge eine include-Deklaration vor dem `params`-Block hinzu

Füge die include-Deklaration oberhalb des `params`-Blocks ein und fülle sie entsprechend aus.

=== "Danach"

    ```groovy title="hello-modules.nf" linenums="3" hl_lines="4"
    // Module einbinden
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
    // Module einbinden
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

Der letzte!

### 4.4. Führe den Workflow aus

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

Das sollte immer noch die gleiche Ausgabe wie zuvor erzeugen.

### Fazit

Du weißt, wie du mehrere Prozesse in einem Workflow modularisierst.

Herzlichen Glückwunsch, du hast all diese Arbeit geleistet und absolut nichts hat sich daran geändert, wie die Pipeline funktioniert!

Spaß beiseite, jetzt ist dein Code modularer, und wenn du dich entscheidest, eine weitere Pipeline zu schreiben, die einen dieser Prozesse aufruft, musst du nur eine kurze `include`-Anweisung eingeben, um das entsprechende Modul zu verwenden.
Das ist besser als Code zu kopieren und einzufügen, denn wenn du später beschließt, das Modul zu verbessern, werden alle deine Pipelines die Verbesserungen erben.

### Wie geht es weiter?

Mach eine kurze Pause, wenn du möchtest.

Wenn du bereit bist, gehe weiter zu [**Teil 5: Hello Containers**](./05_hello_containers.md), um zu lernen, wie du Container verwendest, um Software-Abhängigkeiten bequemer und reproduzierbarer zu verwalten.

---

## Quiz

<quiz>
Was ist ein Modul in Nextflow?
- [ ] Eine Konfigurationsdatei
- [x] Eine eigenständige Datei, die Prozess-Definitionen enthalten kann
- [ ] Eine Workflow-Definition
- [ ] Ein Kanal-Operator

Mehr erfahren: [2. Erstelle ein Modul für `sayHello()`](#2-create-a-module-for-sayhello)
</quiz>

<quiz>
Welche Konvention wird typischerweise zum Speichern von Modul-Dateien verwendet?
- [ ] Im gleichen Verzeichnis wie der Workflow
- [ ] In einem `bin/`-Verzeichnis
- [x] In einem `modules/`-Verzeichnis
- [ ] In einem `lib/`-Verzeichnis

Mehr erfahren: [1. Erstelle ein Verzeichnis zum Speichern von Modulen](#1-create-a-directory-to-store-modules)
</quiz>

<quiz>
Was ist die korrekte Syntax, um ein Modul zu verwenden?

- [ ] `#!groovy import { SAYHELLO } from './modules/sayhello.nf'`
- [ ] `#!groovy require { SAYHELLO } from './modules/sayhello.nf'`
- [x] `#!groovy include { SAYHELLO } from './modules/sayhello.nf'`
- [ ] `#!groovy load { SAYHELLO } from './modules/sayhello.nf'`

Mehr erfahren: [2.3. Füge eine include-Deklaration hinzu](#23-add-an-include-declaration-before-the-workflow-block)
</quiz>

<quiz>
Was passiert mit der `-resume`-Funktionalität bei Verwendung von Modulen?
- [ ] Sie funktioniert nicht mehr
- [ ] Sie erfordert zusätzliche Konfiguration
- [x] Sie funktioniert genauso wie zuvor
- [ ] Sie funktioniert nur für lokale Module
</quiz>

<quiz>
Was sind die Vorteile der Verwendung von Modulen? (Wähle alle zutreffenden aus)
- [x] Code-Wiederverwendbarkeit über Workflows hinweg
- [x] Einfachere Wartung
- [x] Bessere Organisation des Workflow-Codes
- [ ] Schnellere Ausführungsgeschwindigkeit
</quiz>
