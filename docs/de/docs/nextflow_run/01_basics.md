# Teil 1: Grundlegende Operationen

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Im ersten Teil des Nextflow Run-Trainingskurses steigen wir mit einem sehr grundlegenden, fachunabhängigen Hello World-Beispiel in das Thema ein, das wir verwenden werden, um wesentliche Operationen zu demonstrieren und auf die entsprechenden Nextflow-Code-Komponenten hinzuweisen.

??? info "Was ist ein Hello World-Beispiel?"

    Ein "Hello World!" ist ein minimalistisches Beispiel, das die grundlegende Syntax und Struktur einer Programmiersprache oder eines Software-Frameworks demonstrieren soll.
    Das Beispiel besteht typischerweise darin, die Phrase "Hello, World!" auf dem Ausgabegerät wie der Konsole oder dem Terminal auszugeben oder in eine Datei zu schreiben.

---

## 1. Ein Hello World direkt ausführen

Demonstrieren wir dieses Konzept mit einem einfachen Befehl, den wir direkt im Terminal ausführen, um zu zeigen, was er tut, bevor wir ihn in Nextflow einpacken.

!!! tip "Tipp"

    Denke daran, dass du dich jetzt im Verzeichnis `nextflow-run/` befinden solltest, wie auf der Seite [Erste Schritte](00_orientation.md) beschrieben.

### 1.1. Das Terminal "Hallo" sagen lassen

Führe folgenden Befehl in deinem Terminal aus.

```bash
echo 'Hello World!'
```

??? success "Befehlsausgabe"

    ```console
    Hello World!
    ```

Dies gibt den Text 'Hello World' direkt im Terminal aus.

### 1.2. Die Ausgabe in eine Datei schreiben

Beim Ausführen von Pipelines geht es hauptsächlich darum, Daten aus Dateien zu lesen und Ergebnisse in andere Dateien zu schreiben, also modifizieren wir den Befehl so, dass die Textausgabe in eine Datei geschrieben wird, um das Beispiel etwas relevanter zu machen.

```bash
echo 'Hello World!' > output.txt
```

??? success "Befehlsausgabe"

    ```console

    ```

Dies gibt nichts im Terminal aus.

### 1.3. Die Ausgabe finden

Der Text 'Hello World' sollte jetzt in der Ausgabedatei sein, die wir angegeben haben, namens `output.txt`.
Du kannst sie im Datei-Explorer öffnen oder von der Befehlszeile aus mit dem `cat`-Tool zum Beispiel.

??? abstract "Dateiinhalt"

    ```console title="output.txt" linenums="1"
    Hello World!
    ```

Das ist es, was wir mit unserem allerersten Nextflow-Workflow replizieren werden.

### Erkenntnisse

Du weißt jetzt, wie man einen einfachen Befehl im Terminal ausführt, der etwas Text ausgibt, und optional, wie man ihn die Ausgabe in eine Datei schreiben lässt.

### Was kommt als Nächstes?

Finde heraus, was es braucht, um einen Nextflow-Workflow auszuführen, der dasselbe Ergebnis erzielt.

---

## 2. Den Workflow ausführen

Wir stellen dir ein Workflow-Script namens `1-hello.nf` zur Verfügung, das einen Eingabe-Gruß über ein Befehlszeilenargument namens `--input` nimmt und eine Textdatei produziert, die diesen Gruß enthält.

Wir werden uns den Code noch nicht ansehen; zuerst schauen wir uns an, wie es aussieht, ihn auszuführen.

### 2.1. Den Workflow starten und die Ausführung überwachen

Führe im Terminal folgenden Befehl aus:

```bash
nextflow run 1-hello.nf --input 'Hello World!'
```

??? success "Befehlsausgabe"

    ```console hl_lines="6"
    N E X T F L O W   ~  version 25.10.2

    Launching `1-hello.nf` [goofy_torvalds] DSL2 - revision: c33d41f479

    executor >  local (1)
    [a3/7be2fa] sayHello | 1 of 1 ✔
    ```

Wenn deine Konsolenausgabe ungefähr so aussieht, dann herzlichen Glückwunsch, du hast gerade deinen ersten Nextflow-Workflow ausgeführt!

Die wichtigste Ausgabe hier ist die letzte Zeile, die in der obigen Ausgabe hervorgehoben ist:

```console
[a3/7be2fa] sayHello | 1 of 1 ✔
```

Das sagt uns, dass der `sayHello`-Prozess einmal erfolgreich ausgeführt wurde (`1 of 1 ✔`).

Das ist großartig, aber du fragst dich vielleicht: wo ist die Ausgabe?

### 2.2. Die Ausgabedatei im `results`-Verzeichnis finden

Dieser Workflow ist so konfiguriert, dass er seine Ausgabe in ein results-Verzeichnis veröffentlicht.
Wenn du dein aktuelles Verzeichnis anschaust, wirst du sehen, dass Nextflow beim Ausführen des Workflows ein neues Verzeichnis namens `results` erstellt hat, sowie ein Unterverzeichnis namens `1-hello` darunter, das eine Datei namens `output.txt` enthält.

```console title="results/"
results
└── 1-hello
    └── output.txt
```

Öffne die Datei; der Inhalt sollte mit dem String übereinstimmen, den du in der Befehlszeile angegeben hast.

```console title="results/1-hello/output.txt" linenums="1"
Hello World!
```

Das ist großartig, unser Workflow hat getan, was er sollte!

Beachte jedoch, dass das 'veröffentlichte' Ergebnis eine Kopie (oder in einigen Fällen ein symbolischer Link) der tatsächlichen Ausgabe ist, die Nextflow produziert hat, als es den Workflow ausführte.

Also werden wir jetzt unter die Haube schauen, um zu sehen, wo Nextflow die Arbeit tatsächlich ausgeführt hat.

!!! warning "Warnung"

    Nicht alle Workflows werden so eingerichtet sein, dass sie Ausgaben in ein results-Verzeichnis veröffentlichen, und/oder die Verzeichnisnamen und -struktur können unterschiedlich sein.
    Etwas später in diesem Abschnitt zeigen wir dir, wie du herausfindest, wo dieses Verhalten festgelegt ist.

### 2.3. Die ursprüngliche Ausgabe und Protokolle im `work/`-Verzeichnis finden

Wenn du einen Workflow ausführst, erstellt Nextflow ein eigenes 'Aufgabenverzeichnis' für jeden einzelnen Aufruf jedes Prozesses im Workflow (=jeden Schritt in der Pipeline).
Für jeden wird es die notwendigen Eingaben bereitstellen, die relevanten Anweisungen ausführen und Ausgaben und Protokolldateien in dieses eine Verzeichnis schreiben, das automatisch mit einem Hash benannt wird, um es einzigartig zu machen.

Alle diese Aufgabenverzeichnisse leben unter einem Verzeichnis namens `work` in deinem aktuellen Verzeichnis (wo du den Befehl ausführst).

Das mag verwirrend klingen, also schauen wir uns an, wie das in der Praxis aussieht.

Wenn wir auf die Konsolenausgabe für den zuvor ausgeführten Workflow zurückgehen, hatten wir diese Zeile:

```console
[a3/7be2fa] sayHello | 1 of 1 ✔
```

Siehst du, wie die Zeile mit `[a3/7be2fa]` beginnt?
Das ist eine abgekürzte Form des Aufgabenverzeichnispfads für diesen einen Prozessaufruf und sagt dir, wo du die Ausgabe des `sayHello`-Prozessaufrufs innerhalb des `work/`-Verzeichnispfads findest.

Du kannst den vollständigen Pfad finden, indem du folgenden Befehl eingibst (ersetze `a3/7be2fa` durch das, was du in deinem eigenen Terminal siehst) und die Tab-Taste drückst, um den Pfad automatisch zu vervollständigen, oder ein Sternchen hinzufügst:

```bash
ls work/a3/7be2fa*
```

Dies sollte den vollständigen Verzeichnispfad ergeben: `work/a3/7be2fa7be2fad5e71e5f49998f795677fd68`

Schauen wir uns an, was dort drin ist.

??? abstract "Verzeichnisinhalt"

    ```console
    work
    └── a3
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

    Wenn du den Inhalt des Aufgabenunterverzeichnisses im VSCode-Datei-Explorer durchsuchst, siehst du alle Dateien sofort.
    Die Protokolldateien sind jedoch so eingestellt, dass sie im Terminal unsichtbar sind, also wenn du `ls` oder `tree` verwenden möchtest, um sie anzuzeigen, musst du die entsprechende Option zum Anzeigen unsichtbarer Dateien setzen.

    ```bash
    tree -a work
    ```

Du solltest sofort die `output.txt`-Datei erkennen, die tatsächlich die ursprüngliche Ausgabe des `sayHello`-Prozesses ist, die in das `results`-Verzeichnis veröffentlicht wurde.
Wenn du sie öffnest, wirst du den `Hello World!`-Gruß wieder finden.

```console title="work/a3/7be2fa7be2fad5e71e5f49998f795677fd68/output.txt"
Hello World!
```

Was ist also mit all den anderen Dateien?

Das sind die Hilfs- und Protokolldateien, die Nextflow als Teil der Aufgabenausführung geschrieben hat:

- **`.command.begin`**: Sentinel-Datei, die erstellt wird, sobald die Aufgabe gestartet wird.
- **`.command.err`**: Fehlermeldungen (`stderr`), die vom Prozessaufruf ausgegeben werden
- **`.command.log`**: Vollständige Protokollausgabe, die vom Prozessaufruf ausgegeben wird
- **`.command.out`**: Reguläre Ausgabe (`stdout`) vom Prozessaufruf
- **`.command.run`**: Vollständiges Script, das von Nextflow ausgeführt wird, um den Prozessaufruf auszuführen
- **`.command.sh`**: Der Befehl, der tatsächlich vom Prozessaufruf ausgeführt wurde
- **`.exitcode`**: Der Exit-Code, der aus dem Befehl resultiert

Die `.command.sh`-Datei ist besonders nützlich, weil sie dir den Hauptbefehl zeigt, den Nextflow ausgeführt hat, ohne die gesamte Buchführung und Aufgaben-/Umgebungseinrichtung.

```console title="work/a3/7be2fa7be2fad5e71e5f49998f795677fd68/command.sh"
#!/bin/bash -ue
echo 'Hello World!' > output.txt

```

Das bestätigt also, dass der Workflow denselben Befehl zusammengestellt hat, den wir zuvor direkt auf der Befehlszeile ausgeführt haben.

Wenn etwas schief geht und du Fehler beheben musst, kann es nützlich sein, das `command.sh`-Script anzusehen, um genau zu überprüfen, welchen Befehl Nextflow basierend auf den Workflow-Anweisungen, Variableninterpolation usw. zusammengestellt hat.

### 2.4. Den Workflow mit verschiedenen Grüßen erneut ausführen

Versuche, den Workflow einige Male mit verschiedenen Werten für das `--input`-Argument erneut auszuführen, und schaue dann die Aufgabenverzeichnisse an.

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
    ├── 67
    │   ├── 134e6317f90726c6c17ad53234a32b
    │   │   ├── ...
    │   │   └── output.txt
    │   └── e029f2e75305874a9ab263d21ebc2c
    │       ├── ...
    │       └── output.txt
    ...
    ```

Du siehst, dass für jeden Lauf ein neues Unterverzeichnis mit einem vollständigen Satz von Ausgabe- und Protokolldateien erstellt wurde.

Im Gegensatz dazu, wenn du das `results`-Verzeichnis anschaust, gibt es immer noch nur einen Satz von Ergebnissen, und der Inhalt der Ausgabedatei entspricht dem, was du zuletzt ausgeführt hast.

??? abstract "Verzeichnisinhalt"

    ```console title="results/"
    results
    └── 1-hello
        └── output.txt
    ```

Das zeigt dir, dass die veröffentlichten Ergebnisse durch nachfolgende Ausführungen überschrieben werden, während die Aufgabenverzeichnisse unter `work/` erhalten bleiben.

### Erkenntnisse

Du weißt, wie man ein einfaches Nextflow-Script ausführt, seine Ausführung überwacht und seine Ausgaben findet.

### Was kommt als Nächstes?

Lerne, wie man ein grundlegendes Nextflow-Script liest und identifiziert, wie seine Komponenten mit seiner Funktionalität zusammenhängen.

---

## 3. Das Hello World Workflow-Starter-Script untersuchen

Was wir dort getan haben, war im Grunde, das Workflow-Script wie eine Black Box zu behandeln.
Jetzt, da wir gesehen haben, was es tut, öffnen wir die Box und schauen hinein.

Unser Ziel hier ist nicht, die Syntax von Nextflow-Code auswendig zu lernen, sondern eine grundlegende Intuition zu entwickeln, was die Hauptkomponenten sind und wie sie organisiert sind.

### 3.1. Die Gesamtstruktur des Codes untersuchen

Du findest das `1-hello.nf`-Script in deinem aktuellen Verzeichnis, das `nextflow-run` sein sollte. Öffne es im Editor-Bereich.

??? full-code "Vollständige Code-Datei"

    ```groovy title="1-hello.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
    * Verwende echo, um 'Hello World!' in eine Datei zu schreiben
    */
    process sayHello {

        input:
        val greeting

        output:
        path 'output.txt'

        script:
        """
        echo '${greeting}' > output.txt
        """
    }

    /*
    * Pipeline-Parameter
    */
    params {
        input: String
    }

    workflow {

        main:
        // Eine Begrüßung ausgeben
        sayHello(params.input)

        publish:
        first_output = sayHello.out
    }

    output {
        first_output {
            path '1-hello'
            mode 'copy'
        }
    }
    ```

Ein Nextflow-Workflow-Script enthält typischerweise eine oder mehrere **Prozess**-Definitionen, den **Workflow** selbst und einige optionale Blöcke wie **params** und **output**.

Jeder **Prozess** beschreibt, welche Operation(en) der entsprechende Schritt in der Pipeline ausführen soll, während der **Workflow** die Datenflusslogik beschreibt, die die verschiedenen Schritte verbindet.

Schauen wir uns zuerst den **Prozess**-Block genauer an, dann schauen wir uns den **Workflow**-Block an.

### 3.2. Die `process`-Definition

Der erste Code-Block beschreibt einen **Prozess**.
Die Prozessdefinition beginnt mit dem Schlüsselwort `process`, gefolgt vom Prozessnamen und schließlich dem Prozess-Body, der durch geschweifte Klammern begrenzt wird.
Der Prozess-Body muss einen Script-Block enthalten, der den auszuführenden Befehl angibt, der alles sein kann, was du in einem Befehlszeilen-Terminal ausführen könntest.

```groovy title="1-hello.nf" linenums="3"
/*
* Verwende echo, um eine Begrüßung in eine Datei zu schreiben
*/
process sayHello {

    input:
    val greeting

    output:
    path 'output.txt'

    script:
    """
    echo '${greeting}' > output.txt
    """
}
```

Hier haben wir einen **Prozess** namens `sayHello`, der eine **Eingabe**-Variable namens `greeting` nimmt und seine **Ausgabe** in eine Datei namens `output.txt` schreibt.

<figure class="excalidraw">
--8<-- "docs/en/docs/nextflow_run/img/sayhello_with_input.svg"
</figure>

Dies ist eine sehr minimale Prozessdefinition, die nur eine `input`-Definition, eine `output`-Definition und das auszuführende `script` enthält.

Die `input`-Definition enthält den `val`-Qualifier, der Nextflow mitteilt, einen Wert irgendeiner Art zu erwarten (kann ein String, eine Zahl, was auch immer sein).

Die `output`-Definition enthält den `path`-Qualifier, der Nextflow mitteilt, dass dies als Pfad behandelt werden soll (umfasst sowohl Verzeichnispfade als auch Dateien).

### 3.3. Die `workflow`-Definition

Der zweite Code-Block beschreibt den **Workflow** selbst.
Die Workflow-Definition beginnt mit dem Schlüsselwort `workflow`, gefolgt von einem optionalen Namen, dann dem Workflow-Body, der durch geschweifte Klammern begrenzt wird.

Hier haben wir einen **Workflow**, der aus einem `main:`-Block und einem `publish:`-Block besteht.
Der `main:`-Block ist der Haupt-Body des Workflows und der `publish:`-Block listet die Ausgaben auf, die in das `results`-Verzeichnis veröffentlicht werden sollen.

```groovy title="1-hello.nf" linenums="27"
workflow {

    main:
    // Eine Begrüßung ausgeben
    sayHello(params.input)

    publish:
    first_output = sayHello.out
}
```

In diesem Fall enthält der `main:`-Block einen Aufruf des `sayHello`-Prozesses und gibt ihm eine Eingabe namens `params.input` als Gruß.

Wie wir gleich ausführlicher besprechen werden, enthält `params.input` den Wert, den wir dem `--input`-Parameter in unserer Befehlszeile gegeben haben.

Der `publish:`-Block listet die Ausgabe des `sayHello()`-Prozessaufrufs auf, die er als `sayHello.out` bezeichnet und den Namen `first_output` gibt (das kann alles sein, was der Workflow-Autor möchte).

Dies ist eine sehr minimale **Workflow**-Definition.
In einer realen Pipeline enthält der Workflow typischerweise mehrere Aufrufe von **Prozessen**, die durch **Channels** verbunden sind, und es können Standardwerte für die variablen Eingaben eingerichtet sein.

Wir werden das in Teil 2 des Kurses behandeln.
Für jetzt schauen wir uns genauer an, wie unser Workflow Eingaben und Ausgaben handhabt.

### 3.4. Das `params`-System der Befehlszeilenparameter

Das `params.input`, das wir dem `sayHello()`-Prozessaufruf geben, ist ein raffiniertes Stück Nextflow-Code und lohnt sich, eine extra Minute darüber zu verbringen.

Wie oben erwähnt, übergeben wir so den Wert des `--input`-Befehlszeilenparameters an den `sayHello()`-Prozessaufruf.
Tatsächlich reicht es aus, einfach `params.someParameterName` zu deklarieren, um dem Workflow einen Parameter namens `--someParameterName` von der Befehlszeile zu geben.

Hier haben wir diese Parameterdeklaration formalisiert, indem wir einen `params`-Block eingerichtet haben, der den Typ der Eingabe angibt, die der Workflow erwartet (Nextflow 25.10.2 und später).

```groovy title="1-hello.nf" linenums="20"
/*
 * Pipeline-Parameter
 */
params {
    input: String
}
```

Unterstützte Typen umfassen `String`, `Integer`, `Float`, `Boolean` und `Path`.

!!! tip "Tipp"

    Workflow-Parameter, die mit dem `params`-System deklariert werden, nehmen immer zwei Bindestriche auf der Befehlszeile (`--`).
    Das unterscheidet sie von Nextflow-Level-Parametern, die nur einen Bindestrich nehmen (`-`).

### 3.5. Die `publish`-Direktive

Am anderen Ende des Workflows haben wir bereits einen Blick auf den `publish:`-Block geworfen.
Das ist eine Hälfte des Ausgabehandhabungssystems; die andere Hälfte ist der `output`-Block, der sich unten befindet.

```groovy title="1-hello.nf" linenums="37"
output {
    first_output {
        path '1-hello'
        mode 'copy'
    }
}
```

Das gibt an, dass die `first_output`-Ausgabe, die im `publish:`-Block aufgelistet ist, in ein Unterverzeichnis namens `1-hello` unter dem Standard-`results`-Ausgabeverzeichnis kopiert werden soll.

Die Zeile `mode 'copy'` überschreibt das Standardverhalten des Systems, das darin besteht, einen symbolischen Link (oder Symlink) auf die ursprüngliche Datei im `work/`-Verzeichnis zu erstellen, anstatt einer richtigen Kopie.

Es gibt mehr Optionen als hier angezeigt, um das Veröffentlichungsverhalten zu steuern; wir werden später einige behandeln.
Du wirst auch sehen, dass wenn ein Workflow mehrere Ausgaben generiert, jede auf diese Weise im `output`-Block aufgelistet wird.

??? info "Ältere Syntax zum Veröffentlichen von Ausgaben mit `publishDir`"

    Bis vor kurzem war die etablierte Methode zum Veröffentlichen von Ausgaben, es auf der Ebene jedes einzelnen Prozesses mit einer `publishDir`-Direktive zu tun.

    Du wirst dieses Code-Muster überall in älteren Nextflow-Pipelines und Prozessmodulen finden, daher ist es wichtig, davon zu wissen.

    Anstatt einen `publish:`-Block im Workflow und einen `output`-Block auf oberster Ebene zu haben, würdest du eine `publishDir`-Zeile in der `sayHello`-Prozessdefinition sehen:

    ```groovy title="Syntax-Beispiel" linenums="1" hl_lines="3"
    process sayHello {

        publishDir 'results/1-hello', mode: 'copy'

        output:
        path 'output.txt'

        script:
        """
        echo 'Hello World!' > output.txt
        """
    }
    ```

    Wir empfehlen jedoch nicht, dies in neuer Arbeit zu verwenden, da es in zukünftigen Versionen der Nextflow-Sprache irgendwann nicht mehr erlaubt sein wird.

### Erkenntnisse

Du weißt jetzt, wie ein einfacher Nextflow-Workflow strukturiert ist und wie die grundlegenden Komponenten mit seiner Funktionalität zusammenhängen.

### Was kommt als Nächstes?

Lerne, deine Workflow-Ausführungen bequem zu verwalten.

---

## 4. Workflow-Ausführungen verwalten

Zu wissen, wie man Workflows startet und Ausgaben abruft, ist großartig, aber du wirst schnell feststellen, dass es einige andere Aspekte des Workflow-Managements gibt, die dir das Leben erleichtern werden.

Hier zeigen wir dir, wie du die `resume`-Funktion nutzt, wenn du denselben Workflow neu starten musst, wie du die Ausführungsprotokolle mit `nextflow log` überprüfst und wie du ältere Work-Verzeichnisse mit `nextflow clean` löschst.

### 4.1. Einen Workflow mit `-resume` neu starten

Manchmal möchtest du eine Pipeline erneut ausführen, die du zuvor bereits gestartet hast, ohne Arbeit zu wiederholen, die bereits erfolgreich abgeschlossen wurde.

Nextflow hat eine Option namens `-resume`, die dir das ermöglicht.
Konkret werden in diesem Modus alle Prozesse, die bereits mit exakt demselben Code, denselben Einstellungen und Eingaben ausgeführt wurden, übersprungen.
Das bedeutet, Nextflow wird nur Prozesse ausführen, die du seit dem letzten Lauf hinzugefügt oder geändert hast, oder denen du neue Einstellungen oder Eingaben gibst.

Es gibt zwei Hauptvorteile dabei:

- Wenn du gerade dabei bist, eine Pipeline zu entwickeln, kannst du schneller iterieren, da du nur die Prozesse ausführen musst, an denen du aktiv arbeitest, um deine Änderungen zu testen.
- Wenn du eine Pipeline in der Produktion ausführst und etwas schief geht, kannst du in vielen Fällen das Problem beheben und die Pipeline neu starten, und sie wird vom Punkt des Fehlers aus fortgesetzt, was dir viel Zeit und Rechenleistung sparen kann.

Um es zu verwenden, füge einfach `-resume` zu deinem Befehl hinzu und führe ihn aus:

```bash
nextflow run 1-hello.nf --input 'Hello World!' -resume
```

??? success "Befehlsausgabe"

    ```console linenums="1"
    N E X T F L O W   ~  version 25.10.2

    Launching `1-hello.nf` [tiny_noyce] DSL2 - revision: c33d41f479

    [a3/7be2fa] sayHello | 1 of 1, cached: 1 ✔
    ```

Die Konsolenausgabe sollte vertraut aussehen, aber es gibt eine Sache, die ein bisschen anders ist als vorher.

Suche nach dem `cached:`-Teil, der in der Prozessstatuszeile (Zeile 5) hinzugefügt wurde, was bedeutet, dass Nextflow erkannt hat, dass es diese Arbeit bereits erledigt hat, und einfach das Ergebnis des vorherigen erfolgreichen Laufs wiederverwendet hat.

Du kannst auch sehen, dass der Work-Unterverzeichnis-Hash derselbe ist wie im vorherigen Lauf.
Nextflow zeigt dir buchstäblich auf die vorherige Ausführung und sagt "Das habe ich schon dort drüben gemacht."

!!! tip "Tipp"

    Wenn du eine Pipeline mit `resume` erneut ausführst, überschreibt Nextflow keine Dateien, die außerhalb des Work-Verzeichnisses von Ausführungen veröffentlicht wurden, die zuvor erfolgreich ausgeführt wurden.

### 4.2. Das Protokoll vergangener Ausführungen überprüfen

Immer wenn du einen Nextflow-Workflow startest, wird eine Zeile in eine Protokolldatei namens `history` geschrieben, unter einem versteckten Verzeichnis namens `.nextflow` im aktuellen Arbeitsverzeichnis.

??? abstract "Dateiinhalt"

    ```txt title=".nextflow/history" linenums="1"
    2025-07-04 19:27:09	1.8s	wise_watson	OK	3539118582ccde68dde471cc2c66295c	a02c9c46-c3c7-4085-9139-d1b9b5b194c8	nextflow run 1-hello.nf --input 'Hello World'
    2025-07-04 19:27:20	2.9s	spontaneous_blackwell	OK	3539118582ccde68dde471cc2c66295c	59a5db23-d83c-4c02-a54e-37ddb73a337e	nextflow run 1-hello.nf --input Bonjour
    ...
    ```

Diese Datei gibt dir den Zeitstempel, Laufnamen, Status, Revisions-ID, Sitzungs-ID und die vollständige Befehlszeile für jeden Nextflow-Lauf, der aus dem aktuellen Arbeitsverzeichnis gestartet wurde.

Eine bequemere Möglichkeit, auf diese Informationen zuzugreifen, ist der Befehl `nextflow log`.

```bash
nextflow log
```

??? success "Befehlsausgabe"

    ```console linenums="1"
    TIMESTAMP               DURATION        RUN NAME                STATUS  REVISION ID     SESSION ID                              COMMAND
    2025-07-04 19:27:09     1.8s            wise_watson             OK       3539118582     a02c9c46-c3c7-4085-9139-d1b9b5b194c8    nextflow run 1-hello.nf --input 'Hello World'
    ...
    ```

Das gibt den Inhalt der Protokolldatei im Terminal aus, ergänzt durch eine Kopfzeile.

Du wirst bemerken, dass sich die Sitzungs-ID ändert, wenn du einen neuen `nextflow run`-Befehl ausführst, AUSSER wenn du die `-resume`-Option verwendest.
In diesem Fall bleibt die Sitzungs-ID gleich.

Nextflow verwendet die Sitzungs-ID, um Lauf-Caching-Informationen unter dem `cache`-Verzeichnis zu gruppieren, das sich ebenfalls unter `.nextflow` befindet.

### 4.3. Ältere Work-Verzeichnisse löschen

Wenn du viele Pipelines ausführst, können sich sehr viele Dateien über viele Unterverzeichnisse ansammeln.
Da die Unterverzeichnisse zufällig benannt sind, ist es schwierig, anhand ihrer Namen zu erkennen, welche älteren gegenüber neueren Läufen entsprechen.

Glücklicherweise enthält Nextflow einen hilfreichen `clean`-Unterbefehl, der automatisch die Work-Unterverzeichnisse für vergangene Läufe löschen kann, die du nicht mehr benötigst.

#### 4.3.1. Löschkriterien festlegen

Es gibt mehrere [Optionen](https://www.nextflow.io/docs/latest/reference/cli.html#clean), um festzulegen, was gelöscht werden soll.

Hier zeigen wir dir ein Beispiel, das alle Unterverzeichnisse von Läufen vor einem bestimmten Lauf löscht, der über seinen Laufnamen angegeben wird.

Schau dir den letzten erfolgreichen Lauf an, bei dem du nicht `-resume` verwendet hast; in unserem Fall war der Laufname `backstabbing_swartz`.

Der Laufname ist der maschinengenerierte zweiteilige String, der in eckigen Klammern in der `Launching (...)`-Konsolenausgabezeile angezeigt wird.
Du kannst auch das Nextflow-Protokoll verwenden, um einen Lauf basierend auf seinem Zeitstempel und/oder seiner Befehlszeile nachzuschlagen.

#### 4.3.2. Einen Probelauf durchführen

Zuerst verwenden wir das Probelauf-Flag `-n`, um zu überprüfen, was bei dem Befehl gelöscht wird:

```bash
nextflow clean -before backstabbing_swartz -n
```

??? success "Befehlsausgabe"

    ```console
    Would remove /workspaces/training/hello-nextflow/work/eb/1a5de36637b475afd88fca7f79e024
    Would remove /workspaces/training/hello-nextflow/work/6b/19b0e002ea13486d3a0344c336c1d0
    Would remove /workspaces/training/hello-nextflow/work/45/9a6dd7ab771f93003d040956282883
    ```

Deine Ausgabe wird andere Aufgabenverzeichnisnamen haben und möglicherweise eine andere Anzahl von Zeilen, aber sie sollte ähnlich wie das Beispiel aussehen.

Wenn du keine Zeilen in der Ausgabe siehst, hast du entweder keinen gültigen Laufnamen angegeben oder es gibt keine vergangenen Läufe zum Löschen. Stelle sicher, dass du `backstabbing_swartz` im Beispielbefehl durch den entsprechenden letzten Laufnamen in deinem Protokoll änderst.

#### 4.3.3. Mit dem Löschen fortfahren

Wenn die Ausgabe wie erwartet aussieht und du mit dem Löschen fortfahren möchtest, führe den Befehl mit dem `-f`-Flag anstelle von `-n` erneut aus:

```bash
nextflow clean -before backstabbing_swartz -f
```

??? success "Befehlsausgabe"

    ```console
    Removed /workspaces/training/hello-nextflow/work/eb/1a5de36637b475afd88fca7f79e024
    Removed /workspaces/training/hello-nextflow/work/6b/19b0e002ea13486d3a0344c336c1d0
    Removed /workspaces/training/hello-nextflow/work/45/9a6dd7ab771f93003d040956282883
    ```

Die Ausgabe sollte ähnlich wie vorher sein, aber jetzt sagt sie 'Removed' anstelle von 'Would remove'.
Beachte, dass dies nicht die zweistelligen Unterverzeichnisse (wie `eb/` oben) entfernt, aber ihren Inhalt leert.

!!! warning "Warnung"

    Das Löschen von Work-Unterverzeichnissen von vergangenen Läufen entfernt sie aus Nextflows Cache und löscht alle Ausgaben, die in diesen Verzeichnissen gespeichert waren.
    Das bedeutet, es bricht Nextflows Fähigkeit, die Ausführung fortzusetzen, ohne die entsprechenden Prozesse erneut auszuführen.

    Du bist verantwortlich dafür, alle Ausgaben zu sichern, die dir wichtig sind! Das ist der Hauptgrund, warum wir den `copy`-Modus anstelle des `symlink`-Modus für die `publish`-Direktive bevorzugen.

### Erkenntnisse

Du weißt, wie du eine Pipeline neu startest, ohne Schritte zu wiederholen, die bereits identisch ausgeführt wurden, das Ausführungsprotokoll überprüfst und den `nextflow clean`-Befehl verwendest, um alte Work-Verzeichnisse zu bereinigen.

### Was kommt als Nächstes?

Mach eine kleine Pause! Du hast gerade die Bausteine der Nextflow-Syntax und grundlegende Nutzungsanweisungen aufgenommen.

Im nächsten Abschnitt dieses Trainings werden wir uns vier sukzessiv realistischere Versionen der Hello World-Pipeline ansehen, die demonstrieren werden, wie Nextflow dir ermöglicht, mehrere Eingaben effizient zu verarbeiten, Workflows auszuführen, die aus mehreren miteinander verbundenen Schritten bestehen, modulare Code-Komponenten zu nutzen und Container für größere Reproduzierbarkeit und Portabilität zu verwenden.

---

## Quiz

<quiz>
Was repräsentiert `[a3/7be2fa]` in der Konsolenausgabezeile `[a3/7be2fa] SAYHELLO | 1 of 1 ✔`?
- [ ] Die Prozess-Versionsnummer
- [ ] Ein eindeutiger Lauf-Identifikator
- [x] Der abgekürzte Pfad zum Work-Verzeichnis der Aufgabe
- [ ] Die Prüfsumme der Ausgabedatei

Mehr erfahren: [2.3. Die ursprüngliche Ausgabe und Protokolle im `work/`-Verzeichnis finden](#23-die-ursprungliche-ausgabe-und-protokolle-im-work-verzeichnis-finden)
</quiz>

<quiz>
Was ist der Zweck der `.command.sh`-Datei in einem Aufgabenverzeichnis?
- [ ] Sie speichert die Konfigurationseinstellungen der Aufgabe
- [x] Sie zeigt den tatsächlichen Befehl, der vom Prozess ausgeführt wurde
- [ ] Sie enthält Fehlermeldungen von fehlgeschlagenen Aufgaben
- [ ] Sie listet Eingabedateien auf, die für die Aufgabe bereitgestellt wurden

Mehr erfahren: [2.3. Die ursprüngliche Ausgabe und Protokolle im `work/`-Verzeichnis finden](#23-die-ursprungliche-ausgabe-und-protokolle-im-work-verzeichnis-finden)
</quiz>

<quiz>
Was passiert mit veröffentlichten Ergebnissen, wenn du einen Workflow ohne `-resume` erneut ausführst?
- [ ] Sie werden in separaten Verzeichnissen mit Zeitstempel aufbewahrt
- [x] Sie werden von der neuen Ausführung überschrieben
- [ ] Nextflow verhindert das Überschreiben und schlägt fehl
- [ ] Sie werden automatisch gesichert

Mehr erfahren: [2.4. Den Workflow mit verschiedenen Grüßen erneut ausführen](#24-den-workflow-mit-verschiedenen-grussen-erneut-ausfuhren)
</quiz>

<quiz>
Was zeigt diese Konsolenausgabe an?

```console
[skipped  ] process > sayHello (1) [100%] 1 of 1, cached: 1 ✔
```

- [ ] Die Aufgabe ist fehlgeschlagen und wurde übersprungen
- [ ] Die Aufgabe wartet in einer Warteschlange
- [x] Nextflow hat Ergebnisse von einer vorherigen identischen Ausführung wiederverwendet
- [ ] Die Aufgabe wurde manuell abgebrochen

Mehr erfahren: [4.1. Einen Workflow mit `-resume` neu starten](#41-einen-workflow-mit--resume-neu-starten)
</quiz>

<quiz>
Wo speichert Nextflow den Ausführungsverlauf, den der Befehl `nextflow log` anzeigt?
- [ ] Im results-Verzeichnis
- [ ] Im work-Verzeichnis
- [x] In der `.nextflow/history`-Datei
- [ ] In `nextflow.config`

Mehr erfahren: [4.2. Das Protokoll vergangener Ausführungen überprüfen](#42-das-protokoll-vergangener-ausfuhrungen-uberprufen)
</quiz>

<quiz>
Was ist der Zweck des `params`-Blocks in einer Workflow-Datei?
- [ ] Um Prozess-Ressourcenanforderungen zu definieren
- [ ] Um den Executor zu konfigurieren
- [x] Um Workflow-Eingabeparameter zu deklarieren und zu typisieren
- [ ] Um Ausgabe-Veröffentlichungsoptionen anzugeben

Mehr erfahren: [3.4. Das params-System der Befehlszeilenparameter](#34-das-params-system-der-befehlszeilenparameter)
</quiz>

<quiz>
Was macht `mode 'copy'` im `output`-Block des Workflows?
- [ ] Erstellt ein Backup des work-Verzeichnisses
- [x] Erstellt eine vollständige Kopie der Dateien anstelle von symbolischen Links
- [ ] Kopiert das Workflow-Script zu results
- [ ] Aktiviert inkrementelles Datei-Kopieren

Mehr erfahren: [3.5. Die publish-Direktive](#35-die-publish-direktive)
</quiz>

<quiz>
Was ist das empfohlene Flag für den `nextflow clean`-Befehl, bevor tatsächlich Dateien gelöscht werden?
- [x] `-n` (Probelauf), um eine Vorschau dessen zu sehen, was gelöscht würde
- [ ] `-v` (verbose), um detaillierte Ausgabe zu sehen
- [ ] `-a` (all), um alle Verzeichnisse auszuwählen
- [ ] `-q` (quiet), um Warnungen zu unterdrücken

Mehr erfahren: [4.3. Ältere Work-Verzeichnisse löschen](#43-altere-work-verzeichnisse-loschen)
</quiz>
