# Teil 1: Hello World

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/tOukLxWCHiA?si=y8lAedhEHWaTV4zd&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } Sieh dir [die gesamte Playlist](https://youtube.com/playlist?list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&si=eF7cLR62goy-lc6n) auf dem Nextflow YouTube-Kanal an.

:green_book: Das Video-Transkript findest du [hier](./transcripts/01_hello_world.md).
///

Im ersten Teil des Hello Nextflow Trainings steigen wir mit einem sehr einfachen, domänenunabhängigen Hello World Beispiel ein, das wir schrittweise erweitern, um die Verwendung grundlegender Nextflow-Logik und -Komponenten zu demonstrieren.

??? info "Was ist ein Hello World Beispiel?"

    Ein "Hello World!" ist ein minimalistisches Beispiel, das die grundlegende Syntax und Struktur einer Programmiersprache oder eines Software-Frameworks demonstrieren soll.
    Das Beispiel besteht typischerweise darin, die Phrase "Hello, World!" auf einem Ausgabegerät wie der Konsole oder dem Terminal auszugeben oder in eine Datei zu schreiben.

---

## 0. Aufwärmen: Ein Hello World Beispiel direkt ausführen

Lass uns das mit einem einfachen Befehl demonstrieren, den wir direkt im Terminal ausführen, um zu zeigen, was er macht, bevor wir ihn in Nextflow einbinden.

!!! tip "Tipp"

    Denk daran, dass du dich jetzt im Verzeichnis `hello-nextflow/` befinden solltest, wie auf der Seite [Erste Schritte](00_orientation.md) beschrieben.

### 0.1. Das Terminal zum Sprechen bringen

Führe den folgenden Befehl in deinem Terminal aus.

```bash
echo 'Hello World!'
```

??? success "Befehlsausgabe"

    ```console
    Hello World!
    ```

Dies gibt den Text 'Hello World' direkt im Terminal aus.

### 0.2. Die Ausgabe in eine Datei schreiben

Pipelines lesen hauptsächlich Daten aus Dateien und schreiben Ergebnisse in andere Dateien. Lass uns den Befehl daher so ändern, dass die Textausgabe in eine Datei geschrieben wird, um das Beispiel etwas relevanter zu machen.

```bash
echo 'Hello World!' > output.txt
```

??? success "Befehlsausgabe"

    ```console

    ```

Dies gibt nichts im Terminal aus.

### 0.3. Die Ausgabe finden

Der Text 'Hello World' sollte sich jetzt in der von uns angegebenen Ausgabedatei namens `output.txt` befinden.
Du kannst sie im Datei-Explorer öffnen oder über die Kommandozeile mit dem `cat`-Befehl, zum Beispiel.

??? abstract "Dateiinhalt"

    ```console title="output.txt" linenums="1"
    Hello World!
    ```

Das werden wir jetzt versuchen, mit unserem allerersten Nextflow-Workflow nachzubilden.

### Fazit

Du weißt jetzt, wie du einen einfachen Befehl im Terminal ausführst, der Text ausgibt, und optional, wie du die Ausgabe in eine Datei schreiben lässt.

### Wie geht es weiter?

Finde heraus, wie das als Nextflow-Workflow aussehen würde.

---

## 1. Das Skript untersuchen und ausführen

Wir stellen dir ein voll funktionsfähiges, wenn auch minimalistisches Workflow-Skript namens `hello-world.nf` zur Verfügung, das dasselbe tut wie zuvor (schreibt 'Hello World!'), aber mit Nextflow.

Zum Einstieg öffnen wir das Workflow-Skript, damit du ein Gefühl für seine Struktur bekommst.
Dann führen wir es aus und suchen nach den Ausgaben.

### 1.1. Den Code untersuchen

Du findest das Skript `hello-world.nf` in deinem aktuellen Verzeichnis, das `hello-nextflow` sein sollte. Öffne es im Editor-Fenster.

??? full-code "Vollständige Code-Datei"

    ```groovy title="hello-world.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
    * echo verwenden, um 'Hello World!' in eine Datei zu schreiben
    */
    process sayHello {

        output:
        path 'output.txt'

        script:
        """
        echo 'Hello World!' > output.txt
        """
    }

    workflow {

        main:
        // Eine Begrüßung ausgeben
        sayHello()
    }
    ```

Ein Nextflow-Workflow-Skript enthält typischerweise eine oder mehrere [**process**](https://nextflow.io/docs/latest/process.html)-Definitionen und den [**workflow**](https://nextflow.io/docs/latest/workflow.html) selbst, plus einige optionale Blöcke (hier nicht vorhanden), die wir später einführen werden.

Jeder **process** beschreibt, welche Operation(en) der entsprechende Schritt in der Pipeline ausführen soll, während der **workflow** die Datenfluss-Logik beschreibt, die die verschiedenen Schritte verbindet.

Wir schauen uns zuerst den **process**-Block genauer an, dann den **workflow**-Block.

#### 1.1.1. Die `process`-Definition

Der erste Codeblock beschreibt einen **process**.

Die Prozess-Definition beginnt mit dem Schlüsselwort `process`, gefolgt vom Prozessnamen und schließlich dem Prozess-Body, der durch geschweifte Klammern begrenzt wird.
Der Prozess-Body muss einen Script-Block enthalten, der den auszuführenden Befehl angibt, was alles sein kann, was du in einem Kommandozeilen-Terminal ausführen könntest.

```groovy title="hello-world.nf" linenums="3"
/*
* echo verwenden, um 'Hello World!' in eine Datei zu schreiben
*/
process sayHello {

    output:
    path 'output.txt'

    script:
    """
    echo 'Hello World!' > output.txt
    """
}
```

Hier haben wir einen **process** namens `sayHello`, der seine **output** in eine Datei namens `output.txt` schreibt.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_world.svg"
</figure>

Dies ist eine sehr minimale Prozess-Definition, die nur eine `output`-Definition und das auszuführende `script` enthält.

Die `output`-Definition enthält den Qualifier `path`, der Nextflow mitteilt, dass dies als Pfad behandelt werden soll (umfasst sowohl Verzeichnispfade als auch Dateien).
Ein anderer häufiger Qualifier ist `val`.

Wichtig ist, dass die Output-Definition nicht _bestimmt_, welche Ausgabe erstellt wird.
Sie _deklariert_ lediglich, was die erwartete Ausgabe ist, damit Nextflow danach suchen kann, sobald die Ausführung abgeschlossen ist.
Dies ist notwendig, um zu überprüfen, dass der Befehl erfolgreich ausgeführt wurde, und um die Ausgabe bei Bedarf an nachfolgende Prozesse weiterzugeben. Ausgaben, die nicht mit dem übereinstimmen, was im Output-Block deklariert ist, werden nicht an nachfolgende Prozesse weitergegeben.

!!! warning "Warnung"

    Dieses Beispiel ist anfällig, weil wir den Ausgabedateinamen an zwei separaten Stellen fest codiert haben (im Script- und im Output-Block).
    Wenn wir das eine ändern, aber nicht das andere, wird das Skript fehlschlagen.
    Später lernst du Möglichkeiten kennen, Variablen zu verwenden, um dieses Problem zu vermeiden.

In einer realen Pipeline enthält ein Prozess normalerweise zusätzliche Blöcke wie Direktiven und Eingaben, die wir gleich einführen werden.

#### 1.1.2. Die `workflow`-Definition

Der zweite Codeblock beschreibt den **workflow** selbst.
Die Workflow-Definition beginnt mit dem Schlüsselwort `workflow`, gefolgt von einem optionalen Namen, dann dem Workflow-Body, der durch geschweifte Klammern begrenzt wird.

Hier haben wir einen **workflow**, der aus einem `main:`-Block besteht (der sagt 'dies ist der Hauptteil des Workflows'), der einen Aufruf des `sayHello`-Prozesses enthält.

```groovy title="hello-world.nf" linenums="17"
workflow {

    main:
    // Eine Begrüßung ausgeben
    sayHello()
}
```

Dies ist eine sehr minimale **workflow**-Definition.
In einer realen Pipeline enthält der Workflow typischerweise mehrere Aufrufe von **processes**, die durch **channels** verbunden sind, und die Prozesse erwarten eine oder mehrere variable **input(s)**.

Du lernst später in diesem Trainingsmodul, wie du variable Eingaben hinzufügst; und du lernst in Teil 3 dieses Kurses, wie du weitere Prozesse hinzufügst und sie durch Kanäle verbindest.

!!! tip "Tipp"

    Technisch gesehen ist die Zeile `main:` für einfache Workflows wie diesen nicht erforderlich, daher kannst du auf Workflows stoßen, die sie nicht haben.
    Aber wir werden sie brauchen, um Workflow-Level-Ausgaben nutzen zu können, also können wir sie gleich von Anfang an einbinden.

### 1.2. Den Workflow ausführen

Code anzuschauen ist nicht annähernd so unterhaltsam wie ihn auszuführen, also lass uns das in der Praxis ausprobieren.

#### 1.2.1. Den Workflow starten und die Ausführung überwachen

Führe im Terminal den folgenden Befehl aus:

```bash
nextflow run hello-world.nf
```

??? success "Befehlsausgabe"

    ```console hl_lines="7"
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-world.nf` [goofy_torvalds] DSL2 - revision: c33d41f479

    executor >  local (1)
    [65/7be2fa] sayHello | 1 of 1 ✔
    ```

Wenn deine Konsolenausgabe ungefähr so aussieht, dann herzlichen Glückwunsch, du hast gerade deinen ersten Nextflow-Workflow ausgeführt!

Die wichtigste Ausgabe hier ist die letzte Zeile, die in der obigen Ausgabe hervorgehoben ist:

```console
[65/7be2fa] sayHello | 1 of 1 ✔
```

Dies sagt uns, dass der `sayHello`-Prozess erfolgreich einmal ausgeführt wurde (`1 of 1 ✔`).

Wichtig ist, dass diese Zeile dir auch sagt, wo du die Ausgabe des `sayHello`-Prozessaufrufs findest.
Schauen wir uns das jetzt an.

#### 1.2.2. Die Ausgabe und Logs im `work`-Verzeichnis finden

Wenn du Nextflow zum ersten Mal in einem bestimmten Verzeichnis ausführst, erstellt es ein Verzeichnis namens `work`, in das es alle Dateien (und alle Symlinks) schreibt, die im Laufe der Ausführung generiert werden.

Innerhalb des `work`-Verzeichnisses organisiert Nextflow Ausgaben und Logs pro Prozessaufruf.
Für jeden Prozessaufruf erstellt Nextflow ein verschachteltes Unterverzeichnis, das mit einem Hash benannt wird, um es eindeutig zu machen, wo es alle notwendigen Eingaben bereitstellt (standardmäßig mit Symlinks), Hilfsdateien schreibt und Logs und alle Ausgaben des Prozesses ausgibt.

Der Pfad zu diesem Unterverzeichnis wird in gekürzter Form in eckigen Klammern in der Konsolenausgabe angezeigt.
Wenn wir uns ansehen, was wir für den oben gezeigten Lauf erhalten haben, beginnt die Konsolen-Log-Zeile für den sayHello-Prozess mit `[65/7be2fa]`. Das entspricht dem folgenden Verzeichnispfad: `work/65/7be2fa7be2fad5e71e5f49998f795677fd68`

Schauen wir uns an, was dort drin ist.

??? abstract "Verzeichnisinhalt"

    ```console
    work
    └── 65
        └── 7be2fad5e71e5f49998f795677fd68
            ├── .command.begin
            ├── .command.err
            ├── .command.log
            ├── .command.out
            ├── .command.run
            ├── .command.sh
            ├── .exitcode
            └── output.txt
    ```

??? question "Siehst du nicht dasselbe?"

    Die genauen Unterverzeichnisnamen werden auf deinem System anders sein.

    Wenn du den Inhalt des Task-Unterverzeichnisses im VSCode-Datei-Explorer durchsuchst, siehst du alle Dateien sofort.
    Die Log-Dateien sind jedoch im Terminal auf unsichtbar gesetzt. Wenn du also `ls` oder `tree` verwenden möchtest, um sie anzuzeigen, musst du die entsprechende Option zum Anzeigen unsichtbarer Dateien setzen.

    ```bash
    tree -a work
    ```

Das Erste, was du dir ansehen möchtest, ist die tatsächliche Ausgabe des Workflows, d.h. die `output.txt`-Datei, die vom `sayHello`-Prozess erzeugt wurde.
Öffne sie und du wirst die `Hello World!`-Begrüßung finden, die der Sinn unseres minimalistischen Workflows war.

??? abstract "Dateiinhalt"

    ```console title="output.txt"
    Hello World!
    ```

Es hat funktioniert!

Zugegeben, das mag nach viel Wrapper-Code für ein so kleines Ergebnis erscheinen, aber der Wert all dieses Wrapper-Codes wird deutlicher, sobald wir anfangen, Eingabedateien zu lesen und mehrere Schritte miteinander zu verbinden.

Trotzdem schauen wir uns auch die anderen Dateien in diesem Verzeichnis an. Das sind Hilfs- und Log-Dateien, die von Nextflow als Teil der Task-Ausführung erzeugt wurden.

- **`.command.begin`**: Metadaten zum Beginn der Ausführung des Prozessaufrufs
- **`.command.err`**: Fehlermeldungen (`stderr`), die vom Prozessaufruf ausgegeben wurden
- **`.command.log`**: Vollständige Log-Ausgabe des Prozessaufrufs
- **`.command.out`**: Reguläre Ausgabe (`stdout`) des Prozessaufrufs
- **`.command.run`**: Vollständiges Skript, das von Nextflow ausgeführt wurde, um den Prozessaufruf auszuführen
- **`.command.sh`**: Der Befehl, der tatsächlich vom Prozessaufruf ausgeführt wurde
- **`.exitcode`**: Der Exit-Code, der aus dem Befehl resultierte

Die Datei `.command.sh` ist besonders nützlich, weil sie dir den Hauptbefehl zeigt, den Nextflow ausgeführt hat, ohne die gesamte Buchhaltung und Task-/Umgebungseinrichtung.

??? abstract "Dateiinhalt"

    ```console title=".command.sh"
    #!/bin/bash -ue
    echo 'Hello World!' > output.txt
    ```

Dies entspricht dem, was wir zuvor manuell ausgeführt haben.

In diesem Fall ist es sehr einfach, weil der Prozessbefehl fest codiert war, aber später im Kurs wirst du Prozessbefehle sehen, die eine Interpolation von Variablen beinhalten.
Das macht es besonders wertvoll, genau sehen zu können, wie Nextflow den Code interpretiert hat und welcher Befehl erzeugt wurde, wenn du einen fehlgeschlagenen Lauf debuggst.

### 1.3. Den Workflow erneut ausführen

Versuche, den Workflow ein paar Mal erneut auszuführen, und schaue dir dann die Task-Verzeichnisse unter `work/` an.

??? abstract "Verzeichnisinhalt"

    ```console
    work
    ├── 0f
    │   └── 52b7e07b0e274a80843fca48ed21b8
    │       ├── .command.begin
    │       ├── .command.err
    │       ├── .command.log
    │       ├── .command.out
    │       ├── .command.run
    │       ├── .command.sh
    │       ├── .exitcode
    │       └── output.txt
    ├── 65
        └── 7be2fad5e71e5f49998f795677fd68
    │   │   ├── .command.begin
    │   │   ├── .command.err
    │   │   ├── .command.log
    │   │   ├── .command.out
    │   │   ├── .command.run
    │   │   ├── .command.sh
    │   │   ├── .exitcode
    │   │   └── output.txt
    │   └── e029f2e75305874a9ab263d21ebc2c
    │       ├── .command.begin
    │       ├── .command.err
    │       ├── .command.log
    │       ├── .command.out
    │       ├── .command.run
    │       ├── .command.sh
    │       ├── .exitcode
    │       └── output.txt
    ├── 6c
    │   └── d4fd787e0b01b3c82e85696c297500
    │       ├── .command.begin
    │       ├── .command.err
    │       ├── .command.log
    │       ├── .command.out
    │       ├── .command.run
    │       ├── .command.sh
    │       ├── .exitcode
    │       └── output.txt
    └── e8
        └── ab99fad46ade52905ec973ff39bb80
            ├── .command.begin
            ├── .command.err
            ├── .command.log
            ├── .command.out
            ├── .command.run
            ├── .command.sh
            ├── .exitcode
            └── output.txt
    ```

Du siehst, dass für jeden Lauf ein neues Unterverzeichnis mit einem vollständigen Satz von Ausgabe- und Log-Dateien erstellt wurde.
Dies zeigt dir, dass das mehrfache Ausführen desselben Workflows die Ergebnisse früherer Läufe nicht überschreibt.

### Fazit

Du weißt, wie du ein einfaches Nextflow-Skript entschlüsselst, es ausführst und die Ausgabe sowie relevante Log-Dateien im Work-Verzeichnis findest.

### Wie geht es weiter?

Lerne, wie du die Workflow-Ausgaben an einem bequemeren Ort veröffentlichst.

---

## 2. Ausgaben veröffentlichen

Wie du gerade gelernt hast, ist die von unserer Pipeline erzeugte Ausgabe in einem Arbeitsverzeichnis mehrere Ebenen tief vergraben.
Dies geschieht absichtlich; Nextflow hat die Kontrolle über dieses Verzeichnis und wir sollen nicht damit interagieren.
Das macht es jedoch unpraktisch, Ausgaben abzurufen, die uns wichtig sind.

Glücklicherweise bietet Nextflow eine Möglichkeit, Ausgaben in einem bestimmten Verzeichnis zu veröffentlichen, indem [Workflow-Output-Definitionen](https://nextflow.io/docs/latest/workflow.html#workflow-outputs) verwendet werden.

### 2.1. Grundlegende Verwendung

Dies erfordert zwei neue Code-Teile:

1. Einen `publish:`-Block innerhalb des `workflow`-Bodys, der Prozess-Ausgaben deklariert.
2. Einen `output`-Block im Skript, der Ausgabeoptionen wie Modus und Speicherort angibt.

#### 2.1.1. Die Ausgabe des `sayHello`-Prozesses deklarieren

Wir müssen einen `publish:`-Block zum Workflow-Body hinzufügen (dieselbe Art von Code-Element wie der `main:`-Block) und die Ausgabe des `sayHello()`-Prozesses auflisten.

Füge in der Workflow-Skript-Datei `hello-world.nf` die folgenden Codezeilen hinzu:

=== "Danach"

    ```groovy title="hello-world.nf" linenums="17" hl_lines="7-8"
    workflow {

        main:
        // Eine Begrüßung ausgeben
        sayHello()

        publish:
        first_output = sayHello.out
    }
    ```

=== "Vorher"

    ```groovy title="hello-world.nf" linenums="17"
    workflow {

        main:
        // Eine Begrüßung ausgeben
        sayHello()
    }
    ```

Du siehst, dass wir auf die Ausgabe des Prozesses einfach mit `sayHello().out` verweisen und ihr einen beliebigen Namen zuweisen können, `first_output`.

#### 2.1.2. Einen `output:`-Block zum Skript hinzufügen

Jetzt müssen wir nur noch den `output:`-Block hinzufügen, in dem der Ausgabeverzeichnispfad angegeben wird. Beachte, dass dieser neue Block **außerhalb** und **unterhalb** des `workflow`-Blocks im Skript steht.

Füge in der Workflow-Skript-Datei `hello-world.nf` die folgenden Codezeilen hinzu:

=== "Danach"

    ```groovy title="hello-world.nf" linenums="17" hl_lines="11-15"
    workflow {

        main:
        // Eine Begrüßung ausgeben
        sayHello()

        publish:
        first_output = sayHello.out
    }

    output {
        first_output {
            path '.'
        }
    }
    ```

=== "Vorher"

    ```groovy title="hello-world.nf" linenums="17"
    workflow {

        main:
        // Eine Begrüßung ausgeben
        sayHello()

        publish:
        first_output = sayHello.out
    }
    ```

Wir können dies verwenden, um spezifische Pfade für alle im `workflow`-Block deklarierten Prozess-Ausgaben zuzuweisen.
Später lernst du Möglichkeiten kennen, ausgefeilte Ausgabeverzeichnisstrukturen zu generieren, aber vorerst codieren wir der Einfachheit halber nur einen minimalen Pfad fest.

#### 2.1.3. Den Workflow ausführen

Führe jetzt das modifizierte Workflow-Skript aus:

```bash
nextflow run hello-world.nf
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-world.nf` [jovial_mayer] DSL2 - revision: 35bd3425e5

    executor >  local (1)
    [9f/48ef97] sayHello | 1 of 1 ✔
    ```

Die Terminalausgabe sollte vertraut aussehen. Äußerlich hat sich nichts geändert.

Überprüfe jedoch deinen Datei-Explorer: Diesmal hat Nextflow ein neues Verzeichnis namens `results/` erstellt.

??? abstract "Verzeichnisinhalt"

    ```console hl_lines="10-11 22"
    .
    ├── greetings.csv
    ├── hello-channels.nf
    ├── hello-config.nf
    ├── hello-containers.nf
    ├── hello-modules.nf
    ├── hello-workflow.nf
    ├── hello-world.nf
    ├── nextflow.config
    ├── results
    │   └── output.txt -> /workspaces/training/hello-nextflow/work/9f/48ef97f110b0dbd83635d7cbe288d2/output.txt
    ├── solutions
    │   ├── 1-hello-world
    │   ├── 2-hello-channels
    │   ├── 3-hello-workflow
    │   ├── 4-hello-modules
    │   ├── 5-hello-containers
    │   └── 6-hello-config
    ├── test-params.json
    └── work
        ├── 65
        └── 9f
    ```

Im `results`-Verzeichnis finden wir einen symbolischen Link zur `output.txt`, die im Work-Verzeichnis durch den gerade ausgeführten Befehl erzeugt wurde.

Dies ermöglicht es uns, Ausgabedateien einfach abzurufen, ohne durch das Work-Unterverzeichnis graben zu müssen.

### 2.2. Einen benutzerdefinierten Speicherort festlegen

Einen Standardspeicherort zu haben ist großartig, aber du möchtest vielleicht anpassen, wo die Ergebnisse gespeichert werden und wie sie organisiert sind.

Zum Beispiel möchtest du deine Ausgaben möglicherweise in Unterverzeichnissen organisieren.
Der einfachste Weg, dies zu tun, besteht darin, für jede Ausgabe einen spezifischen Ausgabepfad zuzuweisen.

#### 2.2.1. Den Ausgabepfad ändern

Auch hier ist das Ändern des Veröffentlichungsverhaltens für eine bestimmte Ausgabe wirklich unkompliziert.
Um einen benutzerdefinierten Speicherort festzulegen, bearbeite einfach den `path` entsprechend:

=== "Danach"

    ```groovy title="hello-world.nf" linenums="27" hl_lines="3"
    output {
        first_output {
            path 'hello_world'
        }
    }
    ```

=== "Vorher"

    ```groovy title="hello-world.nf" linenums="27" hl_lines="3"
    output {
        first_output {
            path '.'
        }
    }
    ```

Da dies auf der Ebene der einzelnen Ausgabe festgelegt wird, kannst du verschiedene Speicherorte und Unterverzeichnisse nach deinen Bedürfnissen angeben.

#### 2.2.2. Den Workflow erneut ausführen

Probieren wir es aus.

```bash
nextflow run hello-world.nf
```

??? success "Befehlsausgabe"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-world.nf` [tiny_shaw] DSL2 - revision: 757723adc1

    executor >  local (1)
    [8c/79499c] process > sayHello [100%] 1 of 1 ✔
    ```

Diesmal wird das Ergebnis unter dem angegebenen Unterverzeichnis geschrieben.

??? abstract "Verzeichnisinhalt"

    ```console hl_lines="2-3"
    results/
    ├── hello_world
    │   └── output.txt -> /workspaces/training/hello-nextflow/work/8c/79499c2e506b79e2e01acb808d9d12/output.txt
    └── output.txt -> /workspaces/training/hello-nextflow/work/65/f56f2cd75df1352e106fcdd084b97b/output.txt
    ```

Du siehst, dass das Ergebnis der vorherigen Ausführung noch da ist.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_world_output.svg"
</figure>

Du kannst so viele Verschachtelungsebenen verwenden, wie du möchtest.
Es ist auch möglich, den Prozessnamen oder andere Variablen zu verwenden, um die Verzeichnisse zu benennen, die zur Organisation der Ergebnisse verwendet werden, und es ist möglich, den Standardnamen des obersten Ausgabeverzeichnisses zu ändern (der durch das CLI-Flag `-o` oder die Konfigurationsvariable `outputDir` gesteuert wird).
Wir werden diese Optionen später im Training behandeln.

### 2.3. Den Veröffentlichungsmodus auf Kopieren setzen

Standardmäßig werden die Ausgaben als symbolische Links aus dem `work`-Verzeichnis veröffentlicht.
Das bedeutet, dass es nur eine einzige Datei im Dateisystem gibt.

Das ist großartig, wenn du mit sehr großen Dateien arbeitest, von denen du nicht mehrere Kopien speichern möchtest.
Wenn du jedoch das Work-Verzeichnis irgendwann löschst (wir werden Aufräumoperationen in Kürze behandeln), verlierst du den Zugriff auf die Datei.
Du musst also einen Plan haben, um Kopien aller wichtigen Dateien an einem sicheren Ort zu speichern.

Eine einfache Option ist, den Veröffentlichungsmodus für die Ausgaben, die dir wichtig sind, auf Kopieren umzustellen.

#### 2.3.1. Die Modus-Direktive hinzufügen

Dieser Teil ist wirklich unkompliziert.
Füge einfach `mode 'copy'` zur entsprechenden Workflow-Level-Output-Definition hinzu:

=== "Danach"

    ```groovy title="hello-world.nf" linenums="27" hl_lines="4"
    output {
        first_output {
            path 'hello_world'
            mode 'copy'
        }
    }
    ```

=== "Vorher"

    ```groovy title="hello-world.nf" linenums="27"
    output {
        first_output {
            path 'hello_world'
        }
    }
    ```

Dies setzt den Veröffentlichungsmodus für diese spezifische Ausgabe.

#### 2.3.2. Den Workflow erneut ausführen

Probieren wir es aus.

```bash
nextflow run hello-world.nf
```

??? success "Befehlsausgabe"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-world.nf` [tiny_shaw] DSL2 - revision: 757723adc1

    executor >  local (1)
    [df/521638] process > sayHello [100%] 1 of 1 ✔
    ```

Diesmal ist die Datei, wenn du dir die Ergebnisse ansiehst, eine echte Kopie anstatt nur ein Symlink.

??? abstract "Verzeichnisinhalt"

    ```console hl_lines="3"
    results/
    ├── hello_world
    │   └── output.txt
    └── output.txt -> /workspaces/training/hello-nextflow/work/65/f56f2cd75df1352e106fcdd084b97b/output.txt
    ```

Da auch dies auf der Ebene der einzelnen Ausgabe festgelegt wird, ermöglicht es dir, den Veröffentlichungsmodus granular festzulegen.
Dies wird besonders praktisch sein, wenn wir später zu mehrstufigen Pipelines übergehen, wo du beispielsweise nur finale Ausgaben kopieren und Zwischenausgaben als Symlinks belassen möchtest.

Wie bereits erwähnt, gibt es andere, ausgefeiltere Optionen zur Steuerung, wie Ausgaben veröffentlicht werden.
Wir zeigen dir zu gegebener Zeit auf deiner Nextflow-Reise, wie du sie verwendest.

### 2.4. Hinweis zu prozess-level `publishDir`-Direktiven

Bis vor kurzem war die etablierte Methode zum Veröffentlichen von Ausgaben, dies auf der Ebene jedes einzelnen Prozesses mit einer `publishDir`-Direktive zu tun.

Um das zu erreichen, was wir gerade für die Ausgaben des `sayHello`-Prozesses getan haben, hätten wir stattdessen die folgende Zeile zur Prozess-Definition hinzugefügt:

```groovy title="hello-world.nf" linenums="6" hl_lines="3"
process sayHello {

    publishDir 'results/hello_world', mode: 'copy'

    output:
    path 'output.txt'

    script:
    """
    echo 'Hello World!' > output.txt
    """
}
```

Du wirst dieses Code-Muster immer noch überall in älteren Nextflow-Pipelines und Prozess-Modulen finden, daher ist es wichtig, sich dessen bewusst zu sein.
Wir empfehlen jedoch nicht, es in neuen Arbeiten zu verwenden, da es in zukünftigen Versionen der Nextflow-Sprache irgendwann nicht mehr erlaubt sein wird.

### Fazit

Du weißt, wie du Workflow-Ausgaben an einem bequemeren Ort veröffentlichst.

### Wie geht es weiter?

Lerne, eine variable Eingabe über einen Kommandozeilenparameter bereitzustellen und Standardwerte effektiv zu nutzen.

---

## 3. Eine variable Eingabe verwenden, die über die Kommandozeile übergeben wird

In seinem aktuellen Zustand verwendet unser Workflow eine Begrüßung, die fest in den Prozessbefehl codiert ist.
Wir möchten etwas Flexibilität hinzufügen, indem wir eine Eingabevariable verwenden, damit wir die Begrüßung zur Laufzeit leichter ändern können.

Dies erfordert drei Änderungen an unserem Skript:

1. Den Prozess ändern, um eine variable Eingabe zu erwarten
2. Einen Kommandozeilenparameter einrichten, um Benutzereingaben zu erfassen
3. Die Eingabe im Workflow-Body an den Prozess übergeben

Lass uns diese Änderungen nacheinander vornehmen.

### 3.1. Den `sayHello`-Prozess ändern, um eine variable Eingabe zu erwarten

Wir müssen die Prozess-Definition bearbeiten, um (1) eine Eingabevariable zu akzeptieren und (2) diese Variable in der Kommandozeile zu verwenden.

#### 3.1.1. Einen Input-Block zur Prozess-Definition hinzufügen

Zuerst passen wir die Prozess-Definition an, um eine Eingabe namens `greeting` zu akzeptieren.

Nimm im Prozess-Block die folgende Code-Änderung vor:

=== "Danach"

    ```groovy title="hello-world.nf" linenums="6" hl_lines="3-4"
    process sayHello {

        input:
        val greeting

        output:
        path 'output.txt'
    ```

=== "Vorher"

    ```groovy title="hello-world.nf" linenums="6"
    process sayHello {

        output:
        path 'output.txt'
    ```

Die Variable `greeting` ist mit `val` präfixiert, um Nextflow mitzuteilen, dass es sich um einen Wert handelt (nicht um einen Pfad).

#### 3.1.2. Den Prozessbefehl bearbeiten, um die Eingabevariable zu verwenden

Jetzt tauschen wir den ursprünglichen fest codierten Wert gegen den Wert der Eingabevariable aus, die wir erwarten zu erhalten.

Nimm im Prozess-Block die folgende Code-Änderung vor:

=== "Danach"

    ```groovy title="hello-world.nf" linenums="14" hl_lines="3"
    script:
    """
    echo '${greeting}' > output.txt
    """
    ```

=== "Vorher"

    ```groovy title="hello-world.nf" linenums="14" hl_lines="3"
    script:
    """
    echo 'Hello World!' > output.txt
    """
    ```

Das `$`-Symbol und die geschweiften Klammern (`{ }`) sagen Nextflow, dass dies ein Variablenname ist, der durch den tatsächlichen Eingabewert ersetzt werden muss (=interpoliert).

!!! tip "Tipp"

    Die geschweiften Klammern (`{ }`) waren in früheren Versionen von Nextflow technisch optional, daher siehst du möglicherweise ältere Workflows, in denen dies als `echo '$greeting' > output.txt` geschrieben ist.

Jetzt, da der `sayHello()`-Prozess bereit ist, eine variable Eingabe zu akzeptieren, brauchen wir eine Möglichkeit, dem Prozessaufruf auf Workflow-Ebene einen Eingabewert bereitzustellen.

### 3.2. Einen Kommandozeilenparameter einrichten, um Benutzereingaben zu erfassen

Wir könnten einfach eine Eingabe direkt fest codieren, indem wir den Prozessaufruf `sayHello('Hello World!')` machen.
Wenn wir jedoch echte Arbeit mit unserem Workflow erledigen, möchten wir in der Lage sein, seine Eingaben von der Kommandozeile aus zu steuern, damit wir so etwas tun können:

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_world_input.svg"
</figure>

Glücklicherweise hat Nextflow ein eingebautes Workflow-Parametersystem namens [`params`](https://nextflow.io/docs/latest/config.html#params), das es einfach macht, CLI-Parameter zu deklarieren und zu verwenden.

Die allgemeine Syntax besteht darin, `params.<parameter_name>` zu deklarieren, um Nextflow mitzuteilen, dass ein `--<parameter_name>`-Parameter auf der Kommandozeile erwartet wird.

Hier möchten wir einen Parameter namens `--input` erstellen, also müssen wir `params.input` irgendwo im Workflow deklarieren.
Prinzipiell können wir es überall schreiben; aber da wir es dem `sayHello()`-Prozessaufruf geben möchten, können wir es dort direkt einfügen, indem wir `sayHello(params.input)` schreiben.

Nimm im Workflow-Block die folgende Code-Änderung vor:

=== "Danach"

    ```groovy title="hello-world.nf" linenums="23" hl_lines="2"
    // Eine Begrüßung ausgeben
    sayHello(params.input)
    ```

=== "Vorher"

    ```groovy title="hello-world.nf" linenums="23" hl_lines="2"
    // Eine Begrüßung ausgeben
    sayHello()
    ```

Dies sagt Nextflow, den `sayHello`-Prozess mit dem Wert auszuführen, der über den `--input`-Parameter bereitgestellt wird.

Tatsächlich haben wir die Schritte (2) und (3), die zu Beginn des Abschnitts skizziert wurden, in einem Rutsch erledigt.

### 3.3. Den Workflow-Befehl ausführen

Lass es uns ausführen!

```bash
nextflow run hello-world.nf --input 'Bonjour le monde!'
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-world.nf` [elated_lavoisier] DSL2 - revision: 7c031b42ea

    executor >  local (1)
    [4b/654319] sayHello | 1 of 1 ✔
    ```

Wenn du all diese Änderungen korrekt vorgenommen hast, solltest du eine weitere erfolgreiche Ausführung erhalten.

Stelle sicher, dass du die Ausgabedatei öffnest, um zu überprüfen, dass du jetzt die neue Version der Begrüßung hast.

??? abstract "Dateiinhalt"

    ```console title="results/hello_world/output.txt"
    Bonjour le monde!
    ```

Et voilà!

Beachte, wie die neue Ausführung die im `results`-Verzeichnis veröffentlichte Ausgabedatei überschrieben hat.
Die Ergebnisse der vorherigen Läufe sind jedoch immer noch in den Task-Verzeichnissen unter `work` erhalten.

!!! tip "Tipp"

    Du kannst Nextflow-Level-Parameter leicht von Pipeline-Level-Parametern unterscheiden.

    - Parameter, die für eine Pipeline gelten, nehmen immer einen doppelten Bindestrich (`--`).
    - Parameter, die eine Nextflow-Einstellung ändern, _z.B._ die `-resume`-Funktion, die wir zuvor verwendet haben, nehmen einen einfachen Bindestrich (`-`).

### 3.4. Standardwerte für Kommandozeilenparameter verwenden

Ok, das war praktisch, aber in vielen Fällen macht es Sinn, einen Standardwert für einen bestimmten Parameter anzugeben, damit du ihn nicht bei jedem Lauf angeben musst.

#### 3.4.1. Einen Standardwert für den CLI-Parameter festlegen

Geben wir dem `input`-Parameter einen Standardwert, indem wir ihn vor der Workflow-Definition deklarieren.

```groovy title="hello-world.nf" linenums="20"
/*
 * Pipeline-Parameter
 */
params {
    input: String = 'Holà mundo!'
}
```

Wie du siehst, können wir den Typ der Eingabe angeben, den der Workflow erwartet (Nextflow 25.10.2 und später).
Die Syntax ist `name: Type = default_value`.
Unterstützte Typen sind `String`, `Integer`, `Float`, `Boolean` und `Path`.

!!! info "Info"

    In älteren Workflows siehst du möglicherweise, dass der gesamte `params`-Block nur als `input = 'Holà mundo!'` geschrieben ist.

Wenn du deiner Pipeline weitere Parameter hinzufügst, solltest du sie alle zu diesem Block hinzufügen, unabhängig davon, ob du ihnen einen Standardwert geben musst oder nicht.
Dies macht es einfach, alle konfigurierbaren Parameter auf einen Blick zu finden.

#### 3.4.2. Den Workflow erneut ausführen, ohne den Parameter anzugeben

Jetzt, da du einen Standardwert festgelegt hast, kannst du den Workflow erneut ausführen, ohne einen Wert in der Kommandozeile angeben zu müssen.

```bash
nextflow run hello-world.nf
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-world.nf` [determined_edison] DSL2 - revision: 3539118582

    executor >  local (1)
    [72/394147] sayHello | 1 of 1 ✔
    ```

Die Ausgabe wird am selben Ort wie zuvor sein, aber der Inhalt sollte mit dem neuen Text aktualisiert sein.

??? abstract "Dateiinhalt"

    ```console title="results/hello_world/output.txt"
    Holà mundo!
    ```

Nextflow hat den Standardwert des Begrüßungsparameters verwendet, um die Ausgabe zu erstellen.

#### 3.4.3. Den Standardwert überschreiben

Wenn du den Parameter auf der Kommandozeile angibst, überschreibt der CLI-Wert den Standardwert.

Probiere es aus:

```bash
nextflow run hello-world.nf --input 'Konnichiwa!'
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-world.nf` [elegant_faraday] DSL2 - revision: 3539118582

    executor >  local (1)
    [6f/a12a91] sayHello | 1 of 1 ✔
    ```

Noch einmal solltest du die entsprechende aktualisierte Ausgabe in deinem Ergebnisverzeichnis finden.

??? abstract "Dateiinhalt"

    ```console title="results/hello_world/output.txt"
    Konnichiwa!
    ```

!!! note "Hinweis"

    In Nextflow gibt es mehrere Stellen, an denen du Werte für Parameter angeben kannst.
    Wenn derselbe Parameter an mehreren Stellen auf unterschiedliche Werte gesetzt wird, bestimmt Nextflow, welcher Wert verwendet werden soll, basierend auf der Rangfolge, die [hier](https://www.nextflow.io/docs/latest/config.html) beschrieben ist.

    Wir werden dies in Teil 6 (Konfiguration) ausführlicher behandeln.

### Fazit

Du weißt, wie du eine einfache variable Eingabe verwendest, die zur Laufzeit über einen Kommandozeilenparameter bereitgestellt wird, sowie wie du Standardwerte einrichtest, verwendest und überschreibst.

### Wie geht es weiter?

Lerne, wie du Ausführungen bequemer verwaltest.

---

## 4. Workflow-Ausführungen verwalten

Zu wissen, wie man Workflows startet und Ausgaben abruft, ist großartig, aber du wirst schnell feststellen, dass es einige andere Aspekte der Workflow-Verwaltung gibt, die dir das Leben erleichtern werden, besonders wenn du deine eigenen Workflows entwickelst.

Hier zeigen wir dir, wie du die [`-resume`](https://nextflow.io/docs/latest/cache-and-resume.html)-Funktion verwendest, wenn du denselben Workflow erneut starten musst, wie du das Log vergangener Ausführungen mit [`nextflow log`](https://nextflow.io/docs/latest/reference/cli.html#log) inspizierst und wie du ältere Work-Verzeichnisse mit [`nextflow clean`](https://nextflow.io/docs/latest/reference/cli.html#clean) löschst.

### 4.1. Einen Workflow mit `-resume` erneut starten

Manchmal möchtest du eine Pipeline erneut ausführen, die du bereits zuvor gestartet hast, ohne Schritte zu wiederholen, die bereits erfolgreich abgeschlossen wurden.

Nextflow hat eine Option namens [`-resume`](https://nextflow.io/docs/latest/cache-and-resume.html), die dir dies ermöglicht.
Konkret werden in diesem Modus alle Prozesse übersprungen, die bereits mit genau demselben Code, denselben Einstellungen und denselben Eingaben ausgeführt wurden.
Das bedeutet, dass Nextflow nur Prozesse ausführt, die du seit dem letzten Lauf hinzugefügt oder geändert hast, oder denen du neue Einstellungen oder Eingaben bereitstellst.

Es gibt zwei Hauptvorteile dabei:

- Wenn du mitten in der Entwicklung deiner Pipeline bist, kannst du schneller iterieren, da du nur den/die Prozess(e) ausführen musst, an dem/denen du aktiv arbeitest, um deine Änderungen zu testen.
- Wenn du eine Pipeline in Produktion ausführst und etwas schief geht, kannst du in vielen Fällen das Problem beheben und die Pipeline erneut starten, und sie wird ab dem Punkt des Fehlers fortgesetzt, was dir viel Zeit und Rechenleistung sparen kann.

Um es zu verwenden, füge einfach `-resume` zu deinem Befehl hinzu und führe ihn aus:

```bash
nextflow run hello-world.nf -resume
```

??? success "Befehlsausgabe"

    ```console hl_lines="5"
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-world.nf` [golden_cantor] DSL2 - revision: 35bd3425e5

    [62/49a1f8] sayHello | 1 of 1, cached: 1 ✔
    ```

Die Konsolenausgabe sollte vertraut aussehen, aber es gibt eine Sache, die etwas anders ist als zuvor.

Suche nach dem `cached:`-Teil, der in der Prozess-Statuszeile (Zeile 5) hinzugefügt wurde, was bedeutet, dass Nextflow erkannt hat, dass es diese Arbeit bereits erledigt hat und einfach das Ergebnis des vorherigen erfolgreichen Laufs wiederverwendet.

Du kannst auch sehen, dass der Work-Unterverzeichnis-Hash derselbe ist wie beim vorherigen Lauf.
Nextflow zeigt dir buchstäblich auf die vorherige Ausführung und sagt "Das habe ich dort drüben schon gemacht."

!!! tip "Tipp"

    Wenn du eine Pipeline mit `resume` erneut ausführst, überschreibt Nextflow keine Dateien, die außerhalb des Work-Verzeichnisses von Ausführungen veröffentlicht wurden, die zuvor erfolgreich ausgeführt wurden.

### 4.2. Das Log vergangener Ausführungen inspizieren

Ob du eine neue Pipeline entwickelst oder Pipelines in Produktion ausführst, irgendwann wirst du wahrscheinlich Informationen über vergangene Läufe nachschlagen müssen.
So machst du das.

Wann immer du einen Nextflow-Workflow startest, wird eine Zeile in eine Log-Datei namens `history` geschrieben, die sich in einem versteckten Verzeichnis namens `.nextflow` im aktuellen Arbeitsverzeichnis befindet.

??? abstract "Dateiinhalt"

    ```txt title=".nextflow/history" linenums="1"
    2025-07-04 19:27:09	1.8s	wise_watson	OK	3539118582ccde68dde471cc2c66295c	a02c9c46-c3c7-4085-9139-d1b9b5b194c8	nextflow run 1-hello.nf --input 'Hello World'
    2025-07-04 19:27:20	2.9s	spontaneous_blackwell	OK	3539118582ccde68dde471cc2c66295c	59a5db23-d83c-4c02-a54e-37ddb73a337e	nextflow run 1-hello.nf --input Bonjour
    2025-07-04 19:27:31	1.8s	gigantic_yonath	OK	3539118582ccde68dde471cc2c66295c	5acaa83a-6ad6-4509-bebc-cb25d5d7ddd0	nextflow run 1-hello.nf --input 'Dobry den'
    2025-07-04 19:27:45	2.4s	backstabbing_swartz	OK	3539118582ccde68dde471cc2c66295c	5f4b3269-5b53-404a-956c-cac915fbb74e	nextflow run 1-hello.nf --input Konnichiwa
    2025-07-04 19:27:57	2.1s	goofy_wilson	OK	3539118582ccde68dde471cc2c66295c	5f4b3269-5b53-404a-956c-cac915fbb74e	nextflow run 1-hello.nf --input Konnichiwa -resume
    ```

Diese Datei gibt dir den Zeitstempel, den Laufnamen, den Status, die Revisions-ID, die Sitzungs-ID und die vollständige Kommandozeile für jeden Nextflow-Lauf, der aus dem aktuellen Arbeitsverzeichnis gestartet wurde.

Eine bequemere Möglichkeit, auf diese Informationen zuzugreifen, ist die Verwendung des Befehls `nextflow log`.

```bash
nextflow log
```

??? success "Befehlsausgabe"

    ```console linenums="1"
    TIMESTAMP               DURATION        RUN NAME                STATUS  REVISION ID     SESSION ID                              COMMAND
    2025-07-04 19:27:09     1.8s            wise_watson             OK       3539118582     a02c9c46-c3c7-4085-9139-d1b9b5b194c8    nextflow run 1-hello.nf --input 'Hello World'
    2025-07-04 19:27:20     2.9s            spontaneous_blackwell   OK       3539118582     59a5db23-d83c-4c02-a54e-37ddb73a337e    nextflow run 1-hello.nf --input Bonjour
    2025-07-04 19:27:31     1.8s            gigantic_yonath         OK       3539118582     5acaa83a-6ad6-4509-bebc-cb25d5d7ddd0    nextflow run 1-hello.nf --input 'Dobry den'
    2025-07-04 19:27:45     2.4s            backstabbing_swartz     OK       3539118582     5f4b3269-5b53-404a-956c-cac915fbb74e    nextflow run 1-hello.nf --input Konnichiwa
    2025-07-04 19:27:57     2.1s            goofy_wilson            OK       3539118582     5f4b3269-5b53-404a-956c-cac915fbb74e    nextflow run 1-hello.nf --input Konnichiwa -resume
    ```

Dies gibt den Inhalt der Log-Datei im Terminal aus, ergänzt um eine Kopfzeile.

Du wirst bemerken, dass sich die Sitzungs-ID jedes Mal ändert, wenn du einen neuen `nextflow run`-Befehl ausführst, AUSSER wenn du die `-resume`-Option verwendest.
In diesem Fall bleibt die Sitzungs-ID gleich.

Nextflow verwendet die Sitzungs-ID, um Lauf-Caching-Informationen unter dem `cache`-Verzeichnis zu gruppieren, das sich ebenfalls unter `.nextflow` befindet.

### 4.3. Ältere Work-Verzeichnisse löschen

Während des Entwicklungsprozesses wirst du typischerweise deinen Entwurf der Pipeline viele Male ausführen, was zu einer Ansammlung vieler Dateien über viele Unterverzeichnisse hinweg führen kann.

Glücklicherweise enthält Nextflow einen hilfreichen `clean`-Unterbefehl, der automatisch die Work-Unterverzeichnisse für vergangene Läufe löschen kann, die dich nicht mehr interessieren.

#### 4.3.1. Löschkriterien bestimmen

Es gibt mehrere [Optionen](https://www.nextflow.io/docs/latest/reference/cli.html#clean), um zu bestimmen, was gelöscht werden soll.

Hier zeigen wir dir ein Beispiel, das alle Unterverzeichnisse von Läufen vor einem bestimmten Lauf löscht, der anhand seines Laufnamens angegeben wird.

Suche den letzten erfolgreichen Lauf, bei dem du `-resume` nicht verwendet hast; in unserem Fall war der Laufname `golden_cantor`.

Der Laufname ist die maschinell generierte zweiteilige Zeichenfolge, die in eckigen Klammern in der `Launching (...)`-Konsolenausgabezeile angezeigt wird.
Du kannst auch das Nextflow-Log verwenden, um einen Lauf anhand seines Zeitstempels und/oder seiner Kommandozeile nachzuschlagen.

#### 4.3.2. Einen Testlauf durchführen

Zuerst verwenden wir das Testlauf-Flag `-n`, um zu überprüfen, was bei dem Befehl gelöscht wird:

```bash
nextflow clean -before golden_cantor -n
```

??? success "Befehlsausgabe"

    ```console
    Would remove /workspaces/training/hello-nextflow/work/a3/7be2fad5e71e5f49998f795677fd68
    ```

Deine Ausgabe wird andere Task-Verzeichnisnamen haben und möglicherweise eine andere Anzahl von Zeilen, aber sie sollte ähnlich wie das Beispiel aussehen.

Wenn du keine Zeilen ausgegeben siehst, hast du entweder keinen gültigen Laufnamen angegeben oder es gibt keine vergangenen Läufe zum Löschen. Stelle sicher, dass du `golden_cantor` im Beispielbefehl durch den entsprechenden neuesten Laufnamen in deinem Log ersetzt.

#### 4.3.3. Mit der Löschung fortfahren

Wenn die Ausgabe wie erwartet aussieht und du mit der Löschung fortfahren möchtest, führe den Befehl mit dem `-f`-Flag anstelle von `-n` erneut aus:

```bash
nextflow clean -before golden_cantor -f
```

??? success "Befehlsausgabe"

    ```console
    Removed /workspaces/training/hello-nextflow/work/a3/7be2fad5e71e5f49998f795677fd68
    ```

Die Ausgabe sollte ähnlich wie zuvor sein, aber jetzt mit 'Removed' anstelle von 'Would remove'.
Beachte, dass dies die zweistelligen Unterverzeichnisse (wie `a3/` oben) nicht entfernt, aber deren Inhalt leert.

!!! Warning "Warnung"

    Das Löschen von Work-Unterverzeichnissen aus vergangenen Läufen entfernt sie aus dem Cache von Nextflow und löscht alle Ausgaben, die in diesen Verzeichnissen gespeichert wurden.
    Das bedeutet, dass es die Fähigkeit von Nextflow unterbricht, die Ausführung fortzusetzen, ohne die entsprechenden Prozesse erneut auszuführen.

    Du bist dafür verantwortlich, alle Ausgaben zu speichern, die dir wichtig sind oder auf die du dich verlassen möchtest! Das ist der Hauptgrund, warum wir es vorziehen, den `copy`-Modus anstelle des `symlink`-Modus für die `publish`-Direktive zu verwenden.

### Fazit

Du weißt, wie du Ausgaben in einem bestimmten Verzeichnis veröffentlichst, eine Pipeline erneut startest, ohne Schritte zu wiederholen, die bereits auf identische Weise ausgeführt wurden, und den Befehl `nextflow clean` verwendest, um alte Work-Verzeichnisse aufzuräumen.

Allgemeiner weißt du, wie du einen einfachen Nextflow-Workflow interpretierst, seine Ausführung verwaltest und Ausgaben abrufst.

### Wie geht es weiter?

Mach eine kleine Pause, du hast sie dir verdient!

Wenn du bereit bist, gehe weiter zu [**Teil 2: Hello Channels**](./02_hello_channels.md), um zu lernen, wie du Kanäle verwendest, um Eingaben in deinen Workflow einzuspeisen, was es dir ermöglicht, Nextflows eingebauten Datenfluss-Parallelismus und andere leistungsstarke Funktionen zu nutzen.

---

## Quiz

<quiz>
Was sind die minimal erforderlichen Komponenten eines Nextflow-Prozesses?
- [ ] Nur Input- und Output-Blöcke
- [x] Output- und Script-Blöcke
- [ ] Input-, Output- und Script-Blöcke
- [ ] Nur ein Script-Block

Mehr erfahren: [1.1.1. Die process-Definition](#111-the-process-definition)
</quiz>

<quiz>
Was ist der Zweck des Output-Blocks in einem Prozess?
- [ ] Ergebnisse auf der Konsole ausgeben
- [ ] Dateien im Work-Verzeichnis speichern
- [x] Erwartete Ausgaben des Prozesses deklarieren
- [ ] Umgebungsvariablen definieren

Mehr erfahren: [1.1.1. Die process-Definition](#111-the-process-definition)
</quiz>

<quiz>
Welcher Befehl wird verwendet, um einen Nextflow-Workflow auszuführen?
- [ ] `nextflow start`
- [ ] `nextflow execute`
- [x] `nextflow run`
- [ ] `nextflow launch`
</quiz>

<quiz>
Wenn du dir das Work-Verzeichnis einer Aufgabe ansiehst, welche Datei enthält den tatsächlich ausgeführten Befehl?

```
work/a3/7be2fa.../
├── .command.begin
├── .command.err
├── .command.log
├── .command.out
├── .command.run
├── .command.sh
├── .exitcode
└── output.txt
```

- [ ] `.command.run`
- [x] `.command.sh`
- [ ] `.command.log`
- [ ] `.command.out`

Mehr erfahren: [1.2.2. Die Ausgabe und Logs im `work`-Verzeichnis finden](#122-find-the-output-and-logs-in-the-work-directory)
</quiz>

<quiz>
Was bewirkt das `-resume`-Flag?
- [ ] Startet den Workflow von vorne
- [ ] Pausiert den Workflow
- [x] Überspringt Prozesse, die bereits erfolgreich abgeschlossen wurden
- [ ] Erstellt ein Backup des Workflows

Mehr erfahren: [4.1. Einen Workflow mit `-resume` erneut starten](#41-re-launch-a-workflow-with--resume)
</quiz>

<quiz>
Was ist der Standardmodus für die Veröffentlichung von Workflow-Ausgaben?
- [ ] Dateien in das Ausgabeverzeichnis kopieren
- [x] Symbolische Links im Ausgabeverzeichnis erstellen
- [ ] Dateien in das Ausgabeverzeichnis verschieben
- [ ] Dateien im Ausgabeverzeichnis komprimieren

Mehr erfahren: [2.3. Den Veröffentlichungsmodus auf Kopieren setzen](#23-set-the-publish-mode-to-copy)
</quiz>

<quiz>
Wie übergibst du einen Parameterwert von der Kommandozeile an einen Nextflow-Workflow?
- [ ] `-parameter value`
- [ ] `--parameter:value`
- [x] `--parameter value`
- [ ] `-p parameter=value`

Mehr erfahren: [3.2. Einen Kommandozeilenparameter einrichten, um Benutzereingaben zu erfassen](#32-set-up-a-command-line-parameter-to-capture-user-input)
</quiz>

<quiz>
Wie referenzierst du eine Variable innerhalb eines Nextflow-Script-Blocks?
- [ ] Verwende `%variable%`-Syntax
- [x] Verwende `#!groovy ${variable}`-Syntax
- [ ] Verwende `{{variable}}`-Syntax
- [ ] Verwende `[variable]`-Syntax
</quiz>
