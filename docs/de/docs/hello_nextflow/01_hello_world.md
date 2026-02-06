# Teil 1: Hello World

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<!--
<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/8X2hHI-9vms?si=F0t9LFYLjAWoyRXj&amp;list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } Siehe [die gesamte Playlist](https://www.youtube.com/playlist?list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik) auf dem Nextflow YouTube-Kanal.

:green_book: Das Videotranskript ist [hier](./transcripts/01_hello_world.md) verfügbar.
///
-->

In diesem ersten Teil des Hello Nextflow Trainingskurses steigen wir mit einem sehr einfachen, domänenunabhängigen Hello World Beispiel ein. Dieses bauen wir schrittweise aus, um die Verwendung grundlegender Nextflow-Logik und -Komponenten zu demonstrieren.

??? info "Was ist ein Hello World Beispiel?"

    Ein "Hello World!" ist ein minimalistisches Beispiel, das die grundlegende Syntax und Struktur einer Programmiersprache oder eines Software-Frameworks demonstrieren soll.
    Das Beispiel besteht typischerweise darin, den Text "Hello, World!" auf einem Ausgabegerät wie der Konsole oder dem Terminal auszugeben oder in eine Datei zu schreiben.

---

## 0. Aufwärmen: Ein Hello World Beispiel direkt ausführen

Lass uns dies mit einem einfachen Befehl demonstrieren, den wir direkt im Terminal ausführen, um zu zeigen, was er macht, bevor wir ihn in Nextflow einbetten.

!!! tip "Tipp"

    Denke daran, dass du dich jetzt im Verzeichnis `hello-nextflow/` befinden solltest, wie auf der Seite [Erste Schritte](00_orientation.md) beschrieben.

### 0.1. Das Terminal zur Begrüßung auffordern

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

Bei der Ausführung von Pipelines geht es hauptsächlich darum, Daten aus Dateien zu lesen und Ergebnisse in andere Dateien zu schreiben. Lass uns den Befehl so modifizieren, dass die Textausgabe in eine Datei geschrieben wird, um das Beispiel etwas relevanter zu machen.

```bash
echo 'Hello World!' > output.txt
```

??? success "Befehlsausgabe"

    ```console

    ```

Dies gibt nichts im Terminal aus.

### 0.3. Die Ausgabe finden

Der Text 'Hello World' sollte sich jetzt in der von uns angegebenen Ausgabedatei namens `output.txt` befinden.
Du kannst sie im Datei-Explorer öffnen oder zum Beispiel mit dem Dienstprogramm `cat` von der Befehlszeile aus.

??? abstract "Dateiinhalt"

    ```console title="output.txt" linenums="1"
    Hello World!
    ```

Dies werden wir mit unserem allerersten Nextflow-Workflow zu replizieren versuchen.

### Fazit

Du weißt jetzt, wie man einen einfachen Befehl im Terminal ausführt, der etwas Text ausgibt, und optional, wie man ihn dazu bringt, die Ausgabe in eine Datei zu schreiben.

### Wie geht es weiter?

Finde heraus, wie das als Nextflow-Workflow aussehen würde.

---

## 1. Das Skript untersuchen und ausführen

Wir stellen dir ein voll funktionsfähiges, wenn auch minimalistisches Workflow-Skript namens `hello-world.nf` zur Verfügung, das dasselbe macht wie zuvor (es gibt 'Hello World!' aus), aber mit Nextflow.

Um dir den Einstieg zu erleichtern, öffnen wir das Workflow-Skript, damit du ein Gefühl für seine Struktur bekommst.
Dann führen wir es aus und suchen nach seinen Ausgaben.

### 1.1. Den Code untersuchen

Du findest das Skript `hello-world.nf` in deinem aktuellen Verzeichnis, das `hello-nextflow` sein sollte. Öffne es im Editor-Bereich.

??? full-code "Vollständige Codedatei"

    ```groovy title="hello-world.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
    * Verwende echo, um 'Hello World!' in eine Datei zu schreiben
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

Jeder **process** beschreibt, welche Operation(en) der entsprechende Schritt in der Pipeline ausführen soll, während der **workflow** die Datenflusslogik beschreibt, die die verschiedenen Schritte verbindet.

Wir werden uns zuerst den **process**-Block genauer ansehen, dann werden wir uns den **workflow**-Block ansehen.

#### 1.1.1. Die `process`-Definition

Der erste Codeblock beschreibt einen **process**.

Die Prozessdefinition beginnt mit dem Schlüsselwort `process`, gefolgt vom Prozessnamen und schließlich dem Prozesskörper, der von geschweiften Klammern begrenzt wird.
Der Prozesskörper muss einen script-Block enthalten, der den auszuführenden Befehl angibt, der alles sein kann, was du in einem Befehlszeilen-Terminal ausführen könntest.

```groovy title="hello-world.nf" linenums="3"
/*
* Verwende echo, um 'Hello World!' in eine Datei zu schreiben
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

Dies ist eine sehr minimale Prozessdefinition, die nur eine `output`-Definition und das auszuführende `script` enthält.

Die `output`-Definition enthält den `path`-Qualifier, der Nextflow mitteilt, dass dies als Pfad behandelt werden soll (umfasst sowohl Verzeichnispfade als auch Dateien).
Ein weiterer häufiger Qualifier ist `val`.

Wichtig ist, dass die Output-Definition nicht _bestimmt_, welche Ausgabe erstellt wird.
Sie _deklariert_ lediglich, was die erwartete Ausgabe ist, damit Nextflow nach Abschluss der Ausführung danach suchen kann.
Dies ist notwendig, um zu überprüfen, ob der Befehl erfolgreich ausgeführt wurde, und um die Ausgabe bei Bedarf an nachgelagerte Prozesse weiterzugeben. Ausgaben, die nicht mit dem übereinstimmen, was im output-Block deklariert ist, werden nicht an nachgelagerte Prozesse weitergegeben.

!!! warning "Warnung"

    Dieses Beispiel ist fragil, weil wir den Ausgabedateinamen an zwei verschiedenen Stellen (im script- und im output-Block) fest eincodiert haben.
    Wenn wir einen ändern, aber nicht den anderen, wird das Skript nicht mehr funktionieren.
    Später wirst du lernen, wie du Variablen verwenden kannst, um dieses Problem zu vermeiden.

In einer realen Pipeline enthält ein Prozess normalerweise zusätzliche Blöcke wie Direktiven und Eingaben, die wir gleich einführen werden.

#### 1.1.2. Die `workflow`-Definition

Der zweite Codeblock beschreibt den **workflow** selbst.
Die Workflow-Definition beginnt mit dem Schlüsselwort `workflow`, gefolgt von einem optionalen Namen, dann dem Workflow-Körper, der von geschweiften Klammern begrenzt wird.

Hier haben wir einen **workflow**, der aus einem `main:`-Block (der sagt 'dies ist der Hauptteil des Workflows') besteht, der einen Aufruf des `sayHello`-Prozesses enthält.

```groovy title="hello-world.nf" linenums="17"
workflow {

    main:
    // Eine Begrüßung ausgeben
    sayHello()
}
```

Dies ist eine sehr minimale **workflow**-Definition.
In einer realen Pipeline enthält der Workflow typischerweise mehrere Aufrufe von **Prozessen**, die durch **Channels** verbunden sind, und die Prozesse erwarten eine oder mehrere variable **Eingabe(n)**.

Variable Eingaben behandeln wir später in diesem Modul. In Teil 3 lernst du, wie du mehrere Prozesse hinzufügst und durch Channels verbindest.

!!! tip "Tipp"

    Technisch gesehen ist die `main:`-Zeile für einfache Workflows wie diesen nicht erforderlich, daher könntest du auf Workflows stoßen, die sie nicht haben.
    Aber wir werden sie brauchen, um Workflow-Level-Outputs zu nutzen, also können wir sie gleich von Anfang an einbeziehen.

### 1.2. Den Workflow ausführen

Code anzuschauen macht nicht annähernd so viel Spaß wie ihn auszuführen, also lass uns das in der Praxis ausprobieren.

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

Dies sagt uns, dass der `sayHello`-Prozess einmal erfolgreich ausgeführt wurde (`1 of 1 ✔`).

Wichtig ist, dass diese Zeile dir auch sagt, wo du die Ausgabe des `sayHello`-Prozessaufrufs findest.
Lass uns das jetzt anschauen.

#### 1.2.2. Die Ausgabe und Logs im `work`-Verzeichnis finden

Wenn du Nextflow zum ersten Mal in einem bestimmten Verzeichnis ausführst, erstellt es ein Verzeichnis namens `work`, in das es alle Dateien (und alle Symlinks) schreibt, die im Laufe der Ausführung generiert werden.

Innerhalb des `work`-Verzeichnisses organisiert Nextflow Ausgaben und Logs pro Prozessaufruf.
Für jeden Prozessaufruf erstellt Nextflow ein verschachteltes Unterverzeichnis, das mit einem Hash benannt ist, um es eindeutig zu machen, in dem es alle notwendigen Eingaben bereitstellt (standardmäßig mit Symlinks), Hilfsdateien schreibt und Logs sowie alle Ausgaben des Prozesses ausgibt.

Der Pfad zu diesem Unterverzeichnis wird in verkürzter Form in eckigen Klammern in der Konsolenausgabe angezeigt.
Wenn wir uns ansehen, was wir für den oben gezeigten Lauf bekommen haben, beginnt die Konsolenlogzeile für den sayHello-Prozess mit `[65/7be2fa]`. Das entspricht dem folgenden Verzeichnispfad: `work/65/7be2fa7be2fad5e71e5f49998f795677fd68`

Lass uns einen Blick darauf werfen, was dort drin ist.

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
    Die Logdateien sind jedoch so eingestellt, dass sie im Terminal unsichtbar sind. Wenn du also `ls` oder `tree` verwenden möchtest, um sie anzuzeigen, musst du die entsprechende Option zum Anzeigen unsichtbarer Dateien setzen.

    ```bash
    tree -a work
    ```

Das Erste, was du dir ansehen möchtest, ist die tatsächliche Ausgabe des Workflows, d.h. die Datei `output.txt`, die vom `sayHello`-Prozess erzeugt wurde.
Öffne sie und du findest die Begrüßung `Hello World!`, was der Sinn unseres minimalistischen Workflows war.

??? abstract "Dateiinhalt"

    ```console title="output.txt"
    Hello World!
    ```

Es hat funktioniert!

Zugegeben, das mag wie viel Wrapper-Code für ein so winziges Ergebnis erscheinen, aber der Wert all dieses Wrapper-Codes wird offensichtlicher, sobald wir beginnen, Eingabedateien zu lesen und mehrere Schritte miteinander zu verbinden.

Davon abgesehen, lass uns auch die anderen Dateien in diesem Verzeichnis ansehen. Das sind Hilfs- und Logdateien, die Nextflow im Rahmen der Task-Ausführung erstellt hat.

- **`.command.begin`**: Metadaten zum Beginn der Ausführung des Prozessaufrufs
- **`.command.err`**: Fehlermeldungen (`stderr`), die vom Prozessaufruf ausgegeben wurden
- **`.command.log`**: Vollständige Logausgabe, die vom Prozessaufruf ausgegeben wurde
- **`.command.out`**: Reguläre Ausgabe (`stdout`) des Prozessaufrufs
- **`.command.run`**: Vollständiges Skript, das von Nextflow ausgeführt wurde, um den Prozessaufruf auszuführen
- **`.command.sh`**: Der Befehl, der tatsächlich vom Prozessaufruf ausgeführt wurde
- **`.exitcode`**: Der Exit-Code, der sich aus dem Befehl ergab

Die Datei `.command.sh` ist besonders nützlich, weil sie dir den Hauptbefehl zeigt, den Nextflow ausgeführt hat, ohne all die Buchhaltung und Task-/Umgebungseinrichtung.

??? abstract "Dateiinhalt"

    ```console title=".command.sh"
    #!/bin/bash -ue
    echo 'Hello World!' > output.txt
    ```

Dies entspricht dem, was wir zuvor manuell ausgeführt haben.

In diesem Fall ist es sehr einfach, weil der Prozessbefehl fest eincodiert war, aber später im Kurs wirst du Prozessbefehle sehen, die eine gewisse Interpolation von Variablen beinhalten.
Das macht es besonders wertvoll, genau sehen zu können, wie Nextflow den Code interpretiert hat und welcher Befehl erzeugt wurde, wenn du einen fehlgeschlagenen Lauf beheben möchtest.

### 1.3. Den Workflow erneut ausführen

Versuche, den Workflow ein paar Mal erneut auszuführen, und sieh dir dann die Task-Verzeichnisse unter `work/` an.

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

Du siehst, dass für jeden Lauf ein neues Unterverzeichnis mit einem vollständigen Satz von Ausgabe- und Logdateien erstellt wurde.
Dies zeigt dir, dass das mehrmalige Ausführen desselben Workflows die Ergebnisse vorheriger Läufe nicht überschreibt.

### Fazit

Du weißt, wie man ein einfaches Nextflow-Skript entschlüsselt, es ausführt und die Ausgabe sowie die relevanten Logdateien im work-Verzeichnis findet.

### Wie geht es weiter?

Lerne, wie du die Workflow-Ausgaben an einen bequemeren Ort veröffentlichen kannst.

---

## 2. Ausgaben veröffentlichen

Wie du gerade gelernt hast, ist die von unserer Pipeline erzeugte Ausgabe mehrere Ebenen tief in einem Arbeitsverzeichnis vergraben.
Dies geschieht absichtlich; Nextflow kontrolliert dieses Verzeichnis und wir sollen nicht damit interagieren.
Das macht es jedoch unpraktisch, Ausgaben abzurufen, die uns wichtig sind.

Glücklicherweise bietet Nextflow eine Möglichkeit, Ausgaben mit [Workflow-Output-Definitionen](https://www.nextflow.io/docs/latest/workflow.html#workflow-outputs) in ein bestimmtes Verzeichnis zu veröffentlichen.

### 2.1. Grundlegende Verwendung

Dies wird zwei neue Code-Elemente beinhalten:

1. Einen `publish:`-Block innerhalb des `workflow`-Körpers, der Prozessausgaben deklariert.
2. Einen `output`-Block zum Skript, der Ausgabeoptionen wie Modus und Speicherort angibt.

#### 2.1.1. Die Ausgabe des `sayHello`-Prozesses deklarieren

Wir müssen einen `publish:`-Block zum Workflow-Körper hinzufügen (dasselbe Code-Element wie der `main:`-Block) und die Ausgabe des `sayHello()`-Prozesses auflisten.

Füge in der Workflow-Skriptdatei `hello-world.nf` die folgenden Codezeilen hinzu:

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

Du siehst, dass wir uns auf die Ausgabe des Prozesses einfach mit `sayHello().out` beziehen und ihr einen beliebigen Namen `first_output` zuweisen können.

#### 2.1.2. Einen `output:`-Block zum Skript hinzufügen

Jetzt müssen wir nur noch den `output:`-Block hinzufügen, in dem der Pfad des Ausgabeverzeichnisses angegeben wird. Beachte, dass dieser neue Block **außerhalb** und **unterhalb** des `workflow`-Blocks innerhalb des Skripts sitzt.

Füge in der Workflow-Skriptdatei `hello-world.nf` die folgenden Codezeilen hinzu:

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

Wir können dies verwenden, um allen im `workflow`-Block deklarierten Prozessausgaben spezifische Pfade zuzuweisen.
Später wirst du Möglichkeiten kennenlernen, ausgeklügelte Ausgabeverzeichnisstrukturen zu generieren, aber jetzt codieren wir der Einfachheit halber einen minimalen Pfad fest ein.

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

Im Verzeichnis `results` finden wir einen symbolischen Link zur `output.txt`, die im work-Verzeichnis durch den gerade ausgeführten Befehl erzeugt wurde.

Dies ermöglicht es uns, Ausgabedateien einfach abzurufen, ohne das work-Unterverzeichnis durchsuchen zu müssen.

### 2.2. Einen benutzerdefinierten Speicherort festlegen

Ein Standardspeicherort ist großartig, aber du möchtest vielleicht anpassen, wo die Ergebnisse gespeichert werden und wie sie organisiert sind.

Du möchtest zum Beispiel deine Ausgaben in Unterverzeichnisse organisieren.
Der einfachste Weg, das zu tun, ist, pro Ausgabe einen spezifischen Ausgabepfad zuzuweisen.

#### 2.2.1. Den Ausgabepfad modifizieren

Auch hier ist das Modifizieren des Veröffentlichungsverhaltens für eine bestimmte Ausgabe wirklich einfach.
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

Da dies auf der Ebene der einzelnen Ausgabe festgelegt wird, kannst du verschiedene Speicherorte und Unterverzeichnisse entsprechend deinen Bedürfnissen angeben.

#### 2.2.2. Den Workflow erneut ausführen

Lass es uns ausprobieren.

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

Diesmal wird das Ergebnis in das angegebene Unterverzeichnis geschrieben.

??? abstract "Verzeichnisinhalt"

    ```console hl_lines="2-3"
    results/
    ├── hello_world
    │   └── output.txt -> /workspaces/training/hello-nextflow/work/8c/79499c2e506b79e2e01acb808d9d12/output.txt
    └── output.txt -> /workspaces/training/hello-nextflow/work/65/f56f2cd75df1352e106fcdd084b97b/output.txt
    ```

Du siehst, dass das Ergebnis der vorherigen Ausführung immer noch da ist.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_world_output.svg"
</figure>

Du kannst so viele Verschachtelungsebenen verwenden, wie du möchtest.
Es ist auch möglich, den Prozessnamen oder andere Variablen zu verwenden, um die Verzeichnisse zu benennen, die zur Organisation der Ergebnisse verwendet werden, und es ist möglich, den Standardnamen des Ausgabeverzeichnisses der obersten Ebene zu ändern (der durch das `-o`-CLI-Flag oder die Konfigurationsvariable `outputDir` gesteuert wird).
Wir werden diese Optionen später im Training behandeln.

### 2.3. Den Veröffentlichungsmodus auf Kopieren setzen

Standardmäßig werden die Ausgaben als symbolische Links aus dem `work`-Verzeichnis veröffentlicht.
Das bedeutet, dass es nur eine einzige Datei im Dateisystem gibt.

Das ist großartig, wenn du mit sehr großen Dateien arbeitest, für die du nicht mehrere Kopien speichern möchtest.
Wenn du jedoch irgendwann das work-Verzeichnis löschst (wir werden Bereinigungsoperationen in Kürze behandeln), verlierst du den Zugriff auf die Datei.
Du musst also einen Plan haben, um Kopien aller wichtigen Dateien an einem sicheren Ort zu speichern.

Eine einfache Option ist, den Veröffentlichungsmodus für die Ausgaben, die dir wichtig sind, auf Kopieren umzustellen.

#### 2.3.1. Die mode-Direktive hinzufügen

Dieser Teil ist wirklich einfach.
Füge einfach `mode 'copy'` zur relevanten Workflow-Output-Definition hinzu:

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

Lass es uns ausprobieren.

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

Diesmal, wenn du dir die Ergebnisse ansiehst, ist die Datei eine echte Kopie anstatt nur ein Symlink.

??? abstract "Verzeichnisinhalt"

    ```console hl_lines="3"
    results/
    ├── hello_world
    │   └── output.txt
    └── output.txt -> /workspaces/training/hello-nextflow/work/65/f56f2cd75df1352e106fcdd084b97b/output.txt
    ```

Da dies auch auf der Ebene der einzelnen Ausgabe festgelegt wird, kannst du den Veröffentlichungsmodus sehr granular festlegen.
Dies wird besonders praktisch sein, wenn wir zu mehrstufigen Pipelines übergehen, wo du vielleicht nur endgültige Ausgaben kopieren und Zwischenausgaben als Symlinks belassen möchtest.

Wie bereits erwähnt, gibt es andere, ausgeklügeltere Optionen zur Steuerung, wie Ausgaben veröffentlicht werden.
Wir werden dir zeigen, wie du sie beim Arbeiten mit Nextflow verwendest.

### 2.4. Hinweis zu `publishDir`-Direktiven auf Prozessebene

Bis vor kurzem war die etablierte Methode, Ausgaben zu veröffentlichen, dies auf der Ebene jedes einzelnen Prozesses mit einer `publishDir`-Direktive zu tun.

Um das zu erreichen, was wir gerade für die Ausgaben des `sayHello`-Prozesses getan haben, hätten wir stattdessen die folgende Zeile zur Prozessdefinition hinzugefügt:

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

Du wirst dieses Codemuster immer noch überall in älteren Nextflow-Pipelines und Prozessmodulen finden, daher ist es wichtig, es zu kennen.
Wir empfehlen jedoch nicht, es in neuer Arbeit zu verwenden, da es in zukünftigen Versionen der Nextflow-Sprache eventuell nicht mehr erlaubt sein wird.

### Fazit

Du weißt, wie man Workflow-Ausgaben an einen bequemeren Ort veröffentlicht.

### Wie geht es weiter?

Lerne, wie du eine variable Eingabe über einen Befehlszeilenparameter bereitstellst und Standardwerte effektiv nutzt.

---

## 3. Eine variable Eingabe verwenden, die über die Befehlszeile übergeben wird

In seinem aktuellen Zustand verwendet unser Workflow eine in den Prozessbefehl fest eincodierte Begrüßung.
Wir möchten etwas Flexibilität hinzufügen, indem wir eine Eingabevariable verwenden, damit wir die Begrüßung zur Laufzeit leichter ändern können.

Dies erfordert drei Änderungen an unserem Skript:

1. Den Prozess so ändern, dass er eine variable Eingabe erwartet
2. Einen Befehlszeilenparameter einrichten, um Benutzereingaben zu erfassen
3. Die Eingabe an den Prozess im Workflow-Körper übergeben

Lass uns diese Änderungen nacheinander vornehmen.

### 3.1. Den `sayHello`-Prozess ändern, um eine variable Eingabe zu erwarten

Wir müssen die Prozessdefinition bearbeiten, um (1) eine Eingabevariable zu akzeptieren und (2) diese Variable in der Befehlszeile zu verwenden.

#### 3.1.1. Einen input-Block zur Prozessdefinition hinzufügen

Zuerst passen wir die Prozessdefinition an, um eine Eingabe namens `greeting` zu akzeptieren.

Nimm im Prozessblock die folgende Codeänderung vor:

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

Der Variable `greeting` wird `val` vorangestellt, um Nextflow mitzuteilen, dass es sich um einen Wert handelt (nicht um einen Pfad).

#### 3.1.2. Den Prozessbefehl bearbeiten, um die Eingabevariable zu verwenden

Jetzt tauschen wir den ursprünglichen fest eincodierten Wert gegen den Wert der Eingabevariable aus, die wir erwarten zu erhalten.

Nimm im Prozessblock die folgende Codeänderung vor:

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

    Die geschweiften Klammern (`{ }`) waren in früheren Versionen von Nextflow technisch optional, daher siehst du vielleicht ältere Workflows, in denen dies als `echo '$greeting' > output.txt` geschrieben ist.

Jetzt, da der `sayHello()`-Prozess bereit ist, eine variable Eingabe zu akzeptieren, brauchen wir eine Möglichkeit, dem Prozessaufruf auf Workflow-Ebene einen Eingabewert bereitzustellen.

### 3.2. Einen Befehlszeilenparameter einrichten, um Benutzereingaben zu erfassen

Wir könnten einfach eine Eingabe direkt fest eincodieren, indem wir den Prozessaufruf `sayHello('Hello World!')` machen.
Wenn wir jedoch echte Arbeit mit unserem Workflow erledigen, werden wir seine Eingaben von der Befehlszeile aus steuern wollen, sodass wir so etwas machen können:

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_world_input.svg"
</figure>

Glücklicherweise hat Nextflow ein eingebautes Workflow-Parameter-System namens [`params`](https://nextflow.io/docs/latest/config.html#params), das es einfach macht, CLI-Parameter zu deklarieren und zu verwenden.

Die allgemeine Syntax ist, `params.<parameter_name>` zu deklarieren, um Nextflow mitzuteilen, dass es einen `--<parameter_name>`-Parameter auf der Befehlszeile erwartet.

Hier wollen wir einen Parameter namens `--input` erstellen, also müssen wir irgendwo im Workflow `params.input` deklarieren.
Im Prinzip können wir es überall hinschreiben; aber da wir es dem `sayHello()`-Prozessaufruf geben wollen, können wir es dort direkt einfügen, indem wir `sayHello(params.input)` schreiben.

Nimm im Workflow-Block die folgende Codeänderung vor:

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

Dies sagt Nextflow, den `sayHello`-Prozess mit dem über den `--input`-Parameter bereitgestellten Wert auszuführen.

Tatsächlich haben wir die oben skizzierten Schritte (2) und (3) in einem Zug erledigt.

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

Wenn du alle diese Änderungen richtig gemacht hast, solltest du eine weitere erfolgreiche Ausführung erhalten.

Öffne unbedingt die Ausgabedatei, um zu überprüfen, dass du jetzt die neue Version der Begrüßung hast.

??? abstract "Dateiinhalt"

    ```console title="results/hello_world/output.txt"
    Bonjour le monde!
    ```

Voilà!

Beachte, dass die neue Ausführung die im `results`-Verzeichnis veröffentlichte Ausgabedatei überschrieben hat.
Die Ergebnisse der vorherigen Läufe sind jedoch immer noch in den Task-Verzeichnissen unter `work` erhalten.

!!! tip "Tipp"

    Du kannst Nextflow-Level-Parameter leicht von Pipeline-Level-Parametern unterscheiden.

    - Parameter, die für eine Pipeline gelten, nehmen immer einen doppelten Bindestrich (`--`).
    - Parameter, die eine Nextflow-Einstellung ändern, _z.B._ die `-resume`-Funktion, die wir zuvor verwendet haben, nehmen einen einzelnen Bindestrich (`-`).

### 3.4. Standardwerte für Befehlszeilenparameter verwenden

Ok, das war praktisch, aber in vielen Fällen ist es sinnvoll, einen Standardwert für einen bestimmten Parameter anzugeben, damit du ihn nicht bei jedem Lauf angeben musst.

#### 3.4.1. Einen Standardwert für den CLI-Parameter festlegen

Lass uns dem `input`-Parameter einen Standardwert geben, indem wir ihn vor der Workflow-Definition deklarieren.

```groovy title="hello-world.nf" linenums="20"
/*
 * Pipeline-Parameter
 */
params {
    input: String = 'Holà mundo!'
}
```

Wie du siehst, können wir den Typ der Eingabe angeben, die der Workflow erwartet (Nextflow 25.10.2 und später).
Die Syntax ist `name: Type = default_value`.
Unterstützte Typen sind `String`, `Integer`, `Float`, `Boolean` und `Path`.

!!! info "Info"

    In älteren Workflows siehst du vielleicht, dass der gesamte `params`-Block nur als `input = 'Holà mundo!'` geschrieben ist.

Wenn du deiner Pipeline weitere Parameter hinzufügst, solltest du sie alle zu diesem Block hinzufügen, unabhängig davon, ob du ihnen einen Standardwert geben musst.
So kannst du alle konfigurierbaren Parameter auf einen Blick leicht finden.

#### 3.4.2. Den Workflow erneut ausführen, ohne den Parameter anzugeben

Jetzt, da du einen Standardwert festgelegt hast, kannst du den Workflow erneut ausführen, ohne einen Wert in der Befehlszeile angeben zu müssen.

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

Die Ausgabe wird an derselben Stelle wie zuvor sein, aber der Inhalt sollte mit dem neuen Text aktualisiert sein.

??? abstract "Dateiinhalt"

    ```console title="results/hello_world/output.txt"
    Holà mundo!
    ```

Nextflow hat den Standardwert des greeting-Parameters verwendet, um die Ausgabe zu erstellen.

#### 3.4.3. Den Standardwert überschreiben

Wenn du den Parameter auf der Befehlszeile angibst, überschreibt der CLI-Wert den Standardwert.

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

Auch hier solltest du die entsprechend aktualisierte Ausgabe in deinem results-Verzeichnis finden.

??? abstract "Dateiinhalt"

    ```console title="results/hello_world/output.txt"
    Konnichiwa!
    ```

!!! note "Hinweis"

    In Nextflow gibt es mehrere Stellen, an denen du Werte für Parameter angeben kannst.
    Wenn derselbe Parameter an mehreren Stellen auf verschiedene Werte gesetzt wird, bestimmt Nextflow anhand der [hier](https://www.nextflow.io/docs/latest/config.html) beschriebenen Rangfolge, welcher Wert verwendet wird.

    Wir werden dies in Teil 6 (Konfiguration) ausführlicher behandeln.

### Fazit

Du weißt, wie man eine einfache variable Eingabe verwendet, die zur Laufzeit über einen Befehlszeilenparameter bereitgestellt wird, sowie wie man Standardwerte einrichtet, verwendet und überschreibt.

### Wie geht es weiter?

Lerne, wie du Ausführungen bequemer verwalten kannst.

---

## 4. Workflow-Ausführungen verwalten

Zu wissen, wie man Workflows startet und Ausgaben abruft, ist großartig, aber du wirst schnell feststellen, dass es einige andere Aspekte der Workflow-Verwaltung gibt, die dein Leben einfacher machen werden, besonders wenn du deine eigenen Workflows entwickelst.

Hier zeigen wir dir, wie du die [`-resume`](https://nextflow.io/docs/latest/cache-and-resume.html)-Funktion verwendest, wenn du denselben Workflow erneut starten musst, wie du das Log vergangener Ausführungen mit [`nextflow log`](https://nextflow.io/docs/latest/reference/cli.html#log) inspizierst und wie du ältere work-Verzeichnisse mit [`nextflow clean`](https://nextflow.io/docs/latest/reference/cli.html#clean) löschst.

<!-- Any other cool options we should include? Added log -->

### 4.1. Einen Workflow mit `-resume` erneut starten

Manchmal möchtest du eine Pipeline, die du bereits zuvor gestartet hast, erneut ausführen, ohne Schritte zu wiederholen, die bereits erfolgreich abgeschlossen wurden.

Nextflow hat eine Option namens [`-resume`](https://nextflow.io/docs/latest/cache-and-resume.html), die dir dies ermöglicht.
Konkret werden in diesem Modus alle Prozesse, die bereits mit genau demselben Code, denselben Einstellungen und Eingaben ausgeführt wurden, übersprungen.
Das bedeutet, dass Nextflow nur Prozesse ausführt, die du seit dem letzten Lauf hinzugefügt oder geändert hast, oder denen du neue Einstellungen oder Eingaben bereitstellst.

Es gibt zwei Hauptvorteile:

- Wenn du mitten in der Entwicklung deiner Pipeline bist, kannst du schneller iterieren, da du nur den/die Prozess(e) ausführen musst, an dem/denen du aktiv arbeitest, um deine Änderungen zu testen.
- Wenn du eine Pipeline in der Produktion ausführst und etwas schief geht, kannst du in vielen Fällen das Problem beheben und die Pipeline neu starten, und sie wird von der Fehlerstelle aus weiterlaufen, was dir viel Zeit und Rechenleistung sparen kann.

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

Die Konsolenausgabe sollte vertraut aussehen, aber es gibt eine Sache, die im Vergleich zu vorher ein wenig anders ist.

Achte auf das `cached:`-Bit, das in der Prozessstatuszeile (Zeile 5) hinzugefügt wurde, was bedeutet, dass Nextflow erkannt hat, dass es diese Arbeit bereits erledigt hat und einfach das Ergebnis vom vorherigen erfolgreichen Lauf wiederverwendet hat.

Du kannst auch sehen, dass der work-Unterverzeichnis-Hash derselbe ist wie beim vorherigen Lauf.
Nextflow zeigt buchstäblich auf die vorherige Ausführung und sagt "Das habe ich schon dort drüben gemacht."

!!! tip "Tipp"

    Wenn du eine Pipeline mit `resume` erneut ausführst, überschreibt Nextflow keine Dateien, die außerhalb des work-Verzeichnisses von zuvor erfolgreich ausgeführten Ausführungen veröffentlicht wurden.

### 4.2. Das Log vergangener Ausführungen inspizieren

Ob du eine neue Pipeline entwickelst oder Pipelines in der Produktion ausführst, irgendwann musst du wahrscheinlich Informationen über vergangene Läufe nachschlagen.
So geht's.

Immer wenn du einen Nextflow-Workflow startest, wird eine Zeile in eine Logdatei namens `history` geschrieben, unter einem versteckten Verzeichnis namens `.nextflow` im aktuellen Arbeitsverzeichnis.

??? abstract "Dateiinhalt"

    ```txt title=".nextflow/history" linenums="1"
    2025-07-04 19:27:09	1.8s	wise_watson	OK	3539118582ccde68dde471cc2c66295c	a02c9c46-c3c7-4085-9139-d1b9b5b194c8	nextflow run 1-hello.nf --input 'Hello World'
    2025-07-04 19:27:20	2.9s	spontaneous_blackwell	OK	3539118582ccde68dde471cc2c66295c	59a5db23-d83c-4c02-a54e-37ddb73a337e	nextflow run 1-hello.nf --input Bonjour
    2025-07-04 19:27:31	1.8s	gigantic_yonath	OK	3539118582ccde68dde471cc2c66295c	5acaa83a-6ad6-4509-bebc-cb25d5d7ddd0	nextflow run 1-hello.nf --input 'Dobry den'
    2025-07-04 19:27:45	2.4s	backstabbing_swartz	OK	3539118582ccde68dde471cc2c66295c	5f4b3269-5b53-404a-956c-cac915fbb74e	nextflow run 1-hello.nf --input Konnichiwa
    2025-07-04 19:27:57	2.1s	goofy_wilson	OK	3539118582ccde68dde471cc2c66295c	5f4b3269-5b53-404a-956c-cac915fbb74e	nextflow run 1-hello.nf --input Konnichiwa -resume
    ```

Diese Datei gibt dir den Zeitstempel, Laufnamen, Status, Revisions-ID, Session-ID und die vollständige Befehlszeile für jeden Nextflow-Lauf, der aus dem aktuellen Arbeitsverzeichnis gestartet wurde.

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

Dies gibt den Inhalt der Logdatei im Terminal aus, ergänzt um eine Kopfzeile.

Du wirst bemerken, dass sich die Session-ID jedes Mal ändert, wenn du einen neuen `nextflow run`-Befehl ausführst, AUSSER wenn du die `-resume`-Option verwendest.
In diesem Fall bleibt die Session-ID gleich.

Nextflow verwendet die Session-ID, um Lauf-Caching-Informationen unter dem Verzeichnis `cache` zu gruppieren, das sich ebenfalls unter `.nextflow` befindet.

### 4.3. Ältere work-Verzeichnisse löschen

Während des Entwicklungsprozesses wirst du deine Entwurfs-Pipeline typischerweise sehr oft ausführen, was zur Ansammlung vieler Dateien in vielen Unterverzeichnissen führen kann.

Glücklicherweise enthält Nextflow einen hilfreichen `clean`-Unterbefehl, der automatisch die work-Unterverzeichnisse für vergangene Läufe löschen kann, die dich nicht mehr interessieren.

#### 4.3.1. Löschkriterien bestimmen

Es gibt mehrere [Optionen](https://www.nextflow.io/docs/latest/reference/cli.html#clean), um zu bestimmen, was gelöscht werden soll.

Hier zeigen wir dir ein Beispiel, das alle Unterverzeichnisse von Läufen vor einem bestimmten Lauf löscht, der anhand seines Laufnamens angegeben wird.

Suche den aktuellsten erfolgreichen Lauf, bei dem du `-resume` nicht verwendet hast; in unserem Fall war der Laufname `golden_cantor`.

Der Laufname ist die maschinengenerierte zweiteilige Zeichenkette, die in eckigen Klammern in der Konsolenausgabezeile `Launching (...)` angezeigt wird.
Du kannst auch das Nextflow-Log verwenden, um einen Lauf basierend auf seinem Zeitstempel und/oder seiner Befehlszeile nachzuschlagen.

#### 4.3.2. Einen Probelauf durchführen

Zuerst verwenden wir das Probelauf-Flag `-n`, um zu überprüfen, was bei dem Befehl gelöscht wird:

```bash
nextflow clean -before golden_cantor -n
```

??? success "Befehlsausgabe"

    ```console
    Would remove /workspaces/training/hello-nextflow/work/a3/7be2fad5e71e5f49998f795677fd68
    ```

Deine Ausgabe wird andere Task-Verzeichnisnamen haben und möglicherweise eine andere Anzahl von Zeilen, aber sie sollte ähnlich wie das Beispiel aussehen.

Wenn du keine Zeilen in der Ausgabe siehst, hast du entweder keinen gültigen Laufnamen angegeben oder es gibt keine vergangenen Läufe zu löschen. Stelle sicher, dass du `golden_cantor` im Beispielbefehl durch den entsprechenden letzten Laufnamen in deinem Log ersetzt.

#### 4.3.3. Mit dem Löschen fortfahren

Wenn die Ausgabe wie erwartet aussieht und du mit dem Löschen fortfahren möchtest, führe den Befehl mit dem `-f`-Flag anstelle von `-n` erneut aus:

```bash
nextflow clean -before golden_cantor -f
```

??? success "Befehlsausgabe"

    ```console
    Removed /workspaces/training/hello-nextflow/work/a3/7be2fad5e71e5f49998f795677fd68
    ```

Die Ausgabe sollte ähnlich wie zuvor sein, aber jetzt mit 'Removed' anstelle von 'Would remove'.
Beachte, dass dies nicht die zweistelligen Unterverzeichnisse (wie `a3/` oben) entfernt, aber deren Inhalte leert.

!!! warning "Warnung"

    Das Löschen von work-Unterverzeichnissen aus vergangenen Läufen entfernt sie aus dem Nextflow-Cache und löscht alle Ausgaben, die in diesen Verzeichnissen gespeichert waren.
    Das bedeutet, dass Nextflows Fähigkeit, die Ausführung fortzusetzen, ohne die entsprechenden Prozesse erneut auszuführen, beeinträchtigt wird.

    Du bist verantwortlich für das Speichern aller Ausgaben, die dir wichtig sind oder auf die du dich verlassen willst! Das ist der Hauptgrund, warum wir den `copy`-Modus anstelle des `symlink`-Modus für die `publish`-Direktive bevorzugen.

### Fazit

Du weißt, wie man Ausgaben in ein bestimmtes Verzeichnis veröffentlicht, eine Pipeline erneut startet, ohne Schritte zu wiederholen, die bereits auf identische Weise ausgeführt wurden, und den `nextflow clean`-Befehl verwendet, um alte work-Verzeichnisse zu bereinigen.

Allgemeiner gesagt, weißt du, wie man einen einfachen Nextflow-Workflow interpretiert, seine Ausführung verwaltet und Ausgaben abruft.

### Wie geht es weiter?

Mach eine kleine Pause, du hast sie dir verdient!

Wenn du bereit bist, gehe zu [**Teil 2: Hello Channels**](./02_hello_channels.md), um zu lernen, wie du Channels verwendest, um Eingaben in deinen Workflow einzuspeisen, was es dir ermöglicht, Nextflows eingebaute Datenflussparallelisierung und andere leistungsstarke Funktionen zu nutzen.

---

## Quiz

<quiz>
Was sind die mindestens erforderlichen Komponenten eines Nextflow-Prozesses?
- [ ] Nur input- und output-Blöcke
- [x] Output- und script-Blöcke
- [ ] Input-, output- und script-Blöcke
- [ ] Nur ein script-Block

Mehr erfahren: [1.1.1. Die process-Definition](#111-die-process-definition)
</quiz>

<quiz>
Was ist der Zweck des output-Blocks in einem Prozess?
- [ ] Ergebnisse auf der Konsole ausgeben
- [ ] Dateien im work-Verzeichnis speichern
- [x] Erwartete Ausgaben des Prozesses deklarieren
- [ ] Umgebungsvariablen definieren

Mehr erfahren: [1.1.1. Die process-Definition](#111-die-process-definition)
</quiz>

<quiz>
Welcher Befehl wird verwendet, um einen Nextflow-Workflow auszuführen?
- [ ] `nextflow start`
- [ ] `nextflow execute`
- [x] `nextflow run`
- [ ] `nextflow launch`
</quiz>

<quiz>
Welche Datei im work-Verzeichnis eines Tasks enthält den tatsächlich ausgeführten Befehl?

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

Mehr erfahren: [1.2.2. Die Ausgabe und Logs im `work`-Verzeichnis finden](#122-die-ausgabe-und-logs-im-work-verzeichnis-finden)
</quiz>

<quiz>
Was macht das `-resume`-Flag?
- [ ] Startet den Workflow von Anfang an neu
- [ ] Pausiert den Workflow
- [x] Überspringt Prozesse, die bereits erfolgreich abgeschlossen wurden
- [ ] Erstellt ein Backup des Workflows

Mehr erfahren: [4.1. Einen Workflow mit `-resume` erneut starten](#41-einen-workflow-mit--resume-erneut-starten)
</quiz>

<quiz>
Was ist der Standardmodus für die Veröffentlichung von Workflow-Ausgaben?
- [ ] Dateien in das Ausgabeverzeichnis kopieren
- [x] Symbolische Links im Ausgabeverzeichnis erstellen
- [ ] Dateien in das Ausgabeverzeichnis verschieben
- [ ] Dateien im Ausgabeverzeichnis komprimieren

Mehr erfahren: [2.3. Den Veröffentlichungsmodus auf Kopieren setzen](#23-den-veröffentlichungsmodus-auf-kopieren-setzen)
</quiz>

<quiz>
Wie übergibst du einen Parameterwert an einen Nextflow-Workflow von der Befehlszeile?
- [ ] `-parameter value`
- [ ] `--parameter:value`
- [x] `--parameter value`
- [ ] `-p parameter=value`

Mehr erfahren: [3.2. Einen Befehlszeilenparameter einrichten, um Benutzereingaben zu erfassen](#32-einen-befehlszeilenparameter-einrichten-um-benutzereingaben-zu-erfassen)
</quiz>

<quiz>
Wie referenzierst du eine Variable innerhalb eines Nextflow-script-Blocks?
- [ ] Mit `%variable%`-Syntax
- [x] Mit `#!groovy ${variable}`-Syntax
- [ ] Mit `{{variable}}`-Syntax
- [ ] Mit `[variable]`-Syntax
</quiz>
