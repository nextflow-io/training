# Teil 1: Grundlegende Operationen ausführen

Im ersten Teil des Nextflow Run Trainingskurses steigen wir mit einem sehr einfachen, domänenunabhängigen Hello World Beispiel ein, das wir nutzen werden, um grundlegende Operationen zu demonstrieren und die entsprechenden Nextflow-Code-Komponenten aufzuzeigen.

??? info "Was ist ein Hello World Beispiel?"

    Ein "Hello World!" ist ein minimalistisches Beispiel, das die grundlegende Syntax und Struktur einer Programmiersprache oder eines Software-Frameworks demonstrieren soll.
    Das Beispiel besteht typischerweise darin, die Phrase "Hello, World!" auf einem Ausgabegerät wie der Konsole oder dem Terminal auszugeben oder in eine Datei zu schreiben.

---

## 1. Ein Hello World direkt ausführen

Lass uns dieses Konzept mit einem einfachen Befehl demonstrieren, den wir direkt im Terminal ausführen, um zu zeigen, was er macht, bevor wir ihn in Nextflow einbinden.

!!! tip

    Denk daran, dass du dich jetzt im Verzeichnis `nextflow-run/` befinden solltest, wie auf der Seite [Erste Schritte](00_orientation.md) beschrieben.

### 1.1. Das Terminal "Hallo" sagen lassen

Führe den folgenden Befehl in deinem Terminal aus.

```bash
echo 'Hello World!'
```

??? success "Befehlsausgabe"

    ```console
    Hello World!
    ```

Dies gibt den Text 'Hello World' direkt im Terminal aus.

### 1.2. Die Ausgabe in eine Datei schreiben

Das Ausführen von Pipelines beinhaltet hauptsächlich das Lesen von Daten aus Dateien und das Schreiben von Ergebnissen in andere Dateien. Lass uns daher den Befehl so ändern, dass die Textausgabe in eine Datei geschrieben wird, um das Beispiel etwas relevanter zu machen.

```bash
echo 'Hello World!' > output.txt
```

??? success "Befehlsausgabe"

    ```console

    ```

Dies gibt nichts im Terminal aus.

### 1.3. Die Ausgabe finden

Der Text 'Hello World' sollte sich jetzt in der von uns angegebenen Ausgabedatei namens `output.txt` befinden.
Du kannst sie im Datei-Explorer öffnen oder von der Kommandozeile aus mit dem `cat`-Befehl, zum Beispiel.

??? abstract "Dateiinhalt"

    ```console title="output.txt" linenums="1"
    Hello World!
    ```

Das ist es, was wir mit unserem allerersten Nextflow-Workflow zu replizieren versuchen werden.

### Fazit

Du weißt jetzt, wie du einen einfachen Befehl im Terminal ausführst, der Text ausgibt, und optional, wie du die Ausgabe in eine Datei schreiben lässt.

### Wie geht es weiter?

Finde heraus, was nötig ist, um einen Nextflow-Workflow auszuführen, der das gleiche Ergebnis erzielt.

---

## 2. Den Workflow ausführen

Wir stellen dir ein Workflow-Skript namens `1-hello.nf` zur Verfügung, das eine Eingabebegrüßung über ein Kommandozeilenargument namens `--input` entgegennimmt und eine Textdatei mit dieser Begrüßung erzeugt.

Wir werden uns den Code noch nicht ansehen; schauen wir uns zunächst an, wie es aussieht, ihn auszuführen.

### 2.1. Den Workflow starten und die Ausführung überwachen

Führe im Terminal den folgenden Befehl aus.

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

Dies sagt uns, dass der **process** `sayHello` erfolgreich einmal ausgeführt wurde (`1 of 1 ✔`).

Das ist großartig, aber du fragst dich vielleicht: Wo ist die Ausgabe?

### 2.2. Die Ausgabedatei im Verzeichnis `results` finden

Dieser Workflow ist so konfiguriert, dass er seine Ausgabe in ein Ergebnisverzeichnis veröffentlicht.
Wenn du dir dein aktuelles Verzeichnis ansiehst, wirst du sehen, dass Nextflow beim Ausführen des Workflows ein neues Verzeichnis namens `results` erstellt hat, sowie ein Unterverzeichnis namens `1-hello` darunter, das eine Datei namens `output.txt` enthält.

```console title="results/"
results
└── 1-hello
    └── output.txt
```

Öffne die Datei; der Inhalt sollte mit der Zeichenkette übereinstimmen, die du in der Kommandozeile angegeben hast.

```console title="results/1-hello/output.txt" linenums="1"
Hello World!
```

Das ist großartig, unser Workflow hat getan, was er sollte!

### 2.3. Die Ergebnisse in einem anderen Verzeichnis speichern

Standardmäßig speichert Nextflow Pipeline-Ausgaben in einem Verzeichnis namens `results` in deinem aktuellen Pfad.
Um zu ändern, wohin deine Dateien veröffentlicht werden, verwende das CLI-Flag `-output-dir` (oder kurz `-o`)

!!! danger

    Beachte, dass `--input` zwei Bindestriche hat und `-output-dir` einen!
    Das liegt daran, dass `--input` ein Pipeline-_Parameter_ ist und `-output-dir` ein Nextflow-CLI-Flag auf Kernebene.
    Mehr dazu später.

```bash
nextflow run 1-hello.nf --input 'Hello World!' -output-dir hello_results
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `1-hello.nf` [hungry_celsius] DSL2 - revision: f048d6ea78

    executor >  local (1)
    [a3/1e1535] sayHello [100%] 1 of 1 ✔
    ```

Du solltest sehen, dass deine Ausgaben jetzt in einem Verzeichnis namens `hello_results` statt `results` veröffentlicht werden:

```console title="hello_results/"
hello_results
└── 1-hello
    └── output.txt
```

Die Dateien in diesem Verzeichnis sind genau die gleichen wie zuvor, nur das oberste Verzeichnis ist anders.
Beachte jedoch in beiden Fällen, dass das 'veröffentlichte' Ergebnis eine Kopie (oder in manchen Fällen ein symbolischer Link) der tatsächlichen Ausgabe ist, die Nextflow bei der Ausführung des Workflows erzeugt hat.

Jetzt werden wir einen Blick unter die Haube werfen, um zu sehen, wo Nextflow die Arbeit tatsächlich ausgeführt hat.

!!! Warning

    Nicht alle Workflows sind so eingerichtet, dass sie Ausgaben in ein Ergebnisverzeichnis veröffentlichen, und/oder die Verzeichnisnamen und -struktur können unterschiedlich sein.
    Etwas weiter in diesem Abschnitt zeigen wir dir, wie du herausfindest, wo dieses Verhalten festgelegt ist.

### 2.4. Die ursprüngliche Ausgabe und Logs im Verzeichnis `work/` finden

Wenn du einen Workflow ausführst, erstellt Nextflow für jeden einzelnen Aufruf jedes **process** im Workflow (=jeden Schritt in der Pipeline) ein eigenes 'Aufgabenverzeichnis'.
Für jedes wird es die notwendigen Eingaben bereitstellen, die relevanten Anweisungen ausführen und Ausgaben sowie Logdateien in diesem einen Verzeichnis schreiben, das automatisch mit einem Hash benannt wird, um es eindeutig zu machen.

Alle diese Aufgabenverzeichnisse befinden sich in einem Verzeichnis namens `work` innerhalb deines aktuellen Verzeichnisses (wo du den Befehl ausführst).

Das mag verwirrend klingen, also schauen wir uns an, wie das in der Praxis aussieht.

Zurück zur Konsolenausgabe für den Workflow, den wir zuvor ausgeführt haben, hatten wir diese Zeile:

```console
[a3/1e1535] sayHello [100%] 1 of 1 ✔
```

Siehst du, wie die Zeile mit `[a3/1e1535]` beginnt?
Das ist eine verkürzte Form des Aufgabenverzeichnispfads für diesen einen Prozessaufruf und sagt dir, wo du die Ausgabe des `sayHello`-Prozessaufrufs innerhalb des `work/`-Verzeichnispfads findest.

Du kannst den vollständigen Pfad finden, indem du den folgenden Befehl eingibst (ersetze `a3/1e1535` durch das, was du in deinem eigenen Terminal siehst) und die Tab-Taste drückst, um den Pfad automatisch zu vervollständigen, oder ein Sternchen hinzufügst:

```bash
ls work/a3/1e1535*
```

Dies sollte den vollständigen Verzeichnispfad ergeben: `work/a3/1e153543b0a7f9d2c4735ddb4ab231`

Schauen wir uns an, was sich darin befindet.

??? abstract "Verzeichnisinhalt"

    ```console
    work
    ├── a3
    │   └── 1e153543b0a7f9d2c4735ddb4ab231
    │       ├── .command.begin
    │       ├── .command.err
    │       ├── .command.log
    │       ├── .command.out
    │       ├── .command.run
    │       ├── .command.sh
    │       ├── .exitcode
    │       └── output.txt
    └── a4
        └── aa3694b8808bdcc1135ef4a1187a4d
            ├── .command.begin
            ├── .command.err
            ├── .command.log
            ├── .command.out
            ├── .command.run
            ├── .command.sh
            ├── .exitcode
            └── output.txt
    ```

??? question "Siehst du nicht das Gleiche?"

    Die genauen Unterverzeichnisnamen werden auf deinem System unterschiedlich sein.

    Wenn du den Inhalt des Aufgaben-Unterverzeichnisses im VSCode-Datei-Explorer durchsuchst, siehst du alle Dateien sofort.
    Die Logdateien sind jedoch im Terminal als unsichtbar eingestellt. Wenn du also `ls` oder `tree` verwenden möchtest, um sie anzuzeigen, musst du die entsprechende Option zum Anzeigen unsichtbarer Dateien setzen.

    ```bash
    tree -a work
    ```

Es gibt zwei Verzeichnissätze in `work/`, von den zwei verschiedenen Pipeline-Ausführungen, die wir gemacht haben.
Jede Aufgabenausführung bekommt ihr eigenes, isoliertes Verzeichnis zum Arbeiten.
In diesem Fall hat die Pipeline beide Male das Gleiche getan, daher sind die Inhalte jedes Aufgabenverzeichnisses identisch.

Du solltest die Datei `output.txt` sofort erkennen, die tatsächlich die ursprüngliche Ausgabe des `sayHello`-Prozesses ist, die in das `results`-Verzeichnis veröffentlicht wurde.
Wenn du sie öffnest, findest du wieder die `Hello World!`-Begrüßung.

```console title="work/a3/1e153543b0a7f9d2c4735ddb4ab231/output.txt"
Hello World!
```

Was ist also mit all den anderen Dateien?

Dies sind die Hilfs- und Logdateien, die Nextflow als Teil der Aufgabenausführung geschrieben hat:

- **`.command.begin`**: Sentinel-Datei, die erstellt wird, sobald die Aufgabe gestartet wird.
- **`.command.err`**: Fehlermeldungen (`stderr`), die vom Prozessaufruf ausgegeben wurden
- **`.command.log`**: Vollständige Logausgabe, die vom Prozessaufruf ausgegeben wurde
- **`.command.out`**: Reguläre Ausgabe (`stdout`) vom Prozessaufruf
- **`.command.run`**: Vollständiges Skript, das von Nextflow ausgeführt wurde, um den Prozessaufruf auszuführen
- **`.command.sh`**: Der Befehl, der tatsächlich vom Prozessaufruf ausgeführt wurde
- **`.exitcode`**: Der Exit-Code, der aus dem Befehl resultierte

Die Datei `.command.sh` ist besonders nützlich, weil sie dir den Hauptbefehl zeigt, den Nextflow ausgeführt hat, ohne die gesamte Buchhaltung und Aufgaben-/Umgebungseinrichtung.

```console title="work/a3/1e153543b0a7f9d2c4735ddb4ab231/.command.sh"
#!/bin/bash -ue
echo 'Hello World!' > output.txt

```

Dies bestätigt also, dass der Workflow denselben Befehl zusammengestellt hat, den wir zuvor direkt in der Kommandozeile ausgeführt haben.

Wenn etwas schiefgeht und du herausfinden musst, was passiert ist, kann es nützlich sein, das `command.sh`-Skript anzusehen, um genau zu prüfen, welchen Befehl Nextflow basierend auf den Workflow-Anweisungen, Variableninterpolation usw. zusammengestellt hat.

### 2.5. Den Workflow mit verschiedenen Begrüßungen erneut ausführen

Versuche, den Workflow mehrmals mit verschiedenen Werten für das `--input`-Argument erneut auszuführen, und schaue dir dann die Aufgabenverzeichnisse an.

??? abstract "Verzeichnisinhalt"

    ```console
    work/
    ├── 09
    │   └── 5ea8665939daf6f04724286c9b3c8a
    │       ├── .command.begin
    │       ├── .command.err
    │       ├── .command.log
    │       ├── .command.out
    │       ├── .command.run
    │       ├── .command.sh
    │       ├── .exitcode
    │       └── output.txt
    ├── 92
    │   └── ceb95e05d87621c92a399da9bd2067
    │       ├── .command.begin
    │       ├── .command.err
    │       ├── .command.log
    │       ├── .command.out
    │       ├── .command.run
    │       ├── .command.sh
    │       ├── .exitcode
    │       └── output.txt
    ├── 93
    │   └── 6708dbc20c7efdc6769cbe477061ec
    │       ├── .command.begin
    │       ├── .command.err
    │       ├── .command.log
    │       ├── .command.out
    │       ├── .command.run
    │       ├── .command.sh
    │       ├── .exitcode
    │       └── output.txt
    ├── a3
    │   └── 1e153543b0a7f9d2c4735ddb4ab231
    │       ├── .command.begin
    │       ├── .command.err
    │       ├── .command.log
    │       ├── .command.out
    │       ├── .command.run
    │       ├── .command.sh
    │       ├── .exitcode
    │       └── output.txt
    └── a4
        └── aa3694b8808bdcc1135ef4a1187a4d
            ├── .command.begin
            ├── .command.err
            ├── .command.log
            ├── .command.out
            ├── .command.run
            ├── .command.sh
            ├── .exitcode
            └── output.txt
    ```

Du siehst, dass für jede Ausführung ein neues Unterverzeichnis mit einem vollständigen Satz von Ausgabe- und Logdateien erstellt wurde.

Im Gegensatz dazu gibt es im `results`-Verzeichnis immer noch nur einen Ergebnissatz, und der Inhalt der Ausgabedatei entspricht dem, was du zuletzt ausgeführt hast.

??? abstract "Verzeichnisinhalt"

    ```console title="results/"
    results
    └── 1-hello
        └── output.txt
    ```

Dies zeigt dir, dass die veröffentlichten Ergebnisse durch nachfolgende Ausführungen überschrieben werden, während die Aufgabenverzeichnisse unter `work/` erhalten bleiben.

### Fazit

Du weißt, wie du ein einfaches Nextflow-Skript ausführst, seine Ausführung überwachst und seine Ausgaben findest.

### Wie geht es weiter?

Lerne, wie du ein einfaches Nextflow-Skript liest und erkennst, wie seine Komponenten mit seiner Funktionalität zusammenhängen.

---

## 3. Das Hello World Workflow-Starterskript untersuchen

Was wir dort gemacht haben, war im Grunde, das Workflow-Skript wie eine Black Box zu behandeln.
Jetzt, da wir gesehen haben, was es tut, lass uns die Box öffnen und hineinschauen.

Unser Ziel hier ist nicht, die Syntax von Nextflow-Code auswendig zu lernen, sondern ein grundlegendes Verständnis dafür zu entwickeln, was die Hauptkomponenten sind und wie sie organisiert sind.

### 3.1. Die allgemeine Code-Struktur untersuchen

Du findest das Skript `1-hello.nf` in deinem aktuellen Verzeichnis, das `nextflow-run` sein sollte. Öffne es im Editor-Fenster.

??? full-code "Vollständige Code-Datei"

    ```groovy title="1-hello.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
    * echo verwenden, um 'Hello World!' in eine Datei zu schreiben
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
        // emit a greeting
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

Ein Nextflow-Workflow-Skript enthält typischerweise eine oder mehrere **process**-Definitionen, den **workflow** selbst und einige optionale Blöcke wie **params** und **output**.

Jeder **process** beschreibt, welche Operation(en) der entsprechende Schritt in der Pipeline ausführen soll, während der **workflow** die Datenflusslogik beschreibt, die die verschiedenen Schritte verbindet.

Schauen wir uns zunächst den **process**-Block genauer an, dann den **workflow**-Block.

### 3.2. Die `process`-Definition

Der erste Codeblock beschreibt einen [**process**](https://nextflow.io/docs/latest/process.html).
Die Prozessdefinition beginnt mit dem Schlüsselwort `process`, gefolgt vom Prozessnamen und schließlich dem Prozesskörper, der durch geschweifte Klammern begrenzt wird.
Der Prozesskörper muss einen script-Block enthalten, der den auszuführenden Befehl spezifiziert, was alles sein kann, was du in einem Kommandozeilen-Terminal ausführen könntest.

```groovy title="1-hello.nf" linenums="3"
/*
* echo verwenden, um eine Begrüßung in eine Datei zu schreiben
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

Hier haben wir einen **process** namens `sayHello`, der eine **input**-Variable namens `greeting` entgegennimmt und seine **output** in eine Datei namens `output.txt` schreibt.

<figure class="excalidraw">
--8<-- "docs/en/docs/nextflow_run/img/sayhello_with_input.svg"
</figure>

Dies ist eine sehr minimale Prozessdefinition, die nur eine `input`-Definition, eine `output`-Definition und das auszuführende `script` enthält.

Die `input`-Definition enthält den Qualifier `val`, der Nextflow mitteilt, dass ein Wert irgendeiner Art erwartet wird (kann eine Zeichenkette, eine Zahl, was auch immer sein).

Die `output`-Definition enthält den Qualifier `path`, der Nextflow mitteilt, dass dies als Pfad behandelt werden soll (umfasst sowohl Verzeichnispfade als auch Dateien).

### 3.3. Die `workflow`-Definition

Der zweite Codeblock beschreibt den [**workflow**](https://nextflow.io/docs/latest/workflow.html) selbst.
Die Workflow-Definition beginnt mit dem Schlüsselwort `workflow`, gefolgt von einem optionalen Namen, dann dem Workflow-Körper, der durch geschweifte Klammern begrenzt wird.

Hier haben wir einen **workflow**, der aus einem `main:`-Block und einem `publish:`-Block besteht.
Der `main:`-Block ist der Hauptkörper des Workflows und der `publish:`-Block listet die Ausgaben auf, die im `results`-Verzeichnis veröffentlicht werden sollen.

```groovy title="1-hello.nf" linenums="27"
workflow {

    main:
    // emit a greeting
    sayHello(params.input)

    publish:
    first_output = sayHello.out
}
```

In diesem Fall enthält der `main:`-Block einen Aufruf des `sayHello`-Prozesses und gibt ihm eine Eingabe namens `params.input`, die als Begrüßung verwendet werden soll.

Wie wir gleich ausführlicher besprechen werden, enthält `params.input` den Wert, den wir dem `--input`-Parameter in unserer Kommandozeile gegeben haben.

Der `publish:`-Block listet die Ausgabe des `sayHello()`-Prozessaufrufs auf, auf die er als `sayHello.out` verweist und der er den Namen `first_output` gibt (dies kann alles sein, was der Workflow-Autor möchte).

Dies ist eine sehr minimale **workflow**-Definition.
In einer realen Pipeline enthält der Workflow typischerweise mehrere Aufrufe von **processes**, die durch **channels** verbunden sind, und es können Standardwerte für die variablen Eingaben eingerichtet sein.

Darauf gehen wir in Teil 2 des Kurses ein.
Schauen wir uns zunächst genauer an, wie unser Workflow Eingaben und Ausgaben verarbeitet.

### 3.4. Das `params`-System für Kommandozeilenparameter

Das `params.input`, das wir dem `sayHello()`-Prozessaufruf zur Verfügung stellen, ist ein nettes Stück Nextflow-Code und ist es wert, eine zusätzliche Minute darauf zu verwenden.

Wie oben erwähnt, ist das die Art und Weise, wie wir den Wert des `--input`-Kommandozeilenparameters an den `sayHello()`-Prozessaufruf übergeben.
Tatsächlich reicht es aus, einfach `params.someParameterName` zu deklarieren, um dem Workflow einen Parameter namens `--someParameterName` von der Kommandozeile aus zu geben.

Hier haben wir diese Parameterdeklaration formalisiert, indem wir einen `params`-Block eingerichtet haben, der den Typ der Eingabe spezifiziert, die der Workflow erwartet (Nextflow 25.10.2 und später).

```groovy title="1-hello.nf" linenums="20"
/*
 * Pipeline-Parameter
 */
params {
    input: String
}
```

Unterstützte Typen sind `String`, `Integer`, `Float`, `Boolean` und `Path`.
Um mehr zu erfahren, siehe [Workflow parameters](https://nextflow.io/docs/latest/config.html#workflow-parameters) in der Nextflow-Referenzdokumentation.

!!! tip

    Denk daran, dass _Workflow_-Parameter, die mit dem `params`-System deklariert werden, in der Kommandozeile immer zwei Bindestriche (`--`) haben.
    Dies unterscheidet sie von _Nextflow-Level_-CLI-Flags, die nur einen Bindestrich (`-`) haben.

### 3.5. Die `publish`-Direktive

Am anderen Ende des Workflows haben wir bereits einen Blick auf den `publish:`-Block geworfen.
Das ist die eine Hälfte des Ausgabeverwaltungssystems; die andere Hälfte ist der `output`-Block, der sich weiter unten befindet.

```groovy title="1-hello.nf" linenums="37"
output {
    first_output {
        path '1-hello'
        mode 'copy'
    }
}
```

Dies spezifiziert, dass die `first_output`-Ausgabe, die im `publish:`-Block aufgeführt ist, in ein Unterverzeichnis namens `1-hello` unter dem Standard-`results`-Ausgabeverzeichnis kopiert werden soll.

Die Zeile `mode 'copy'` überschreibt das Standardverhalten des Systems, das darin besteht, einen symbolischen Link (oder Symlink) zur Originaldatei im `work/`-Verzeichnis zu erstellen, anstatt einer echten Kopie.

Es gibt mehr Optionen als hier angezeigt, um das Veröffentlichungsverhalten zu steuern; einige davon werden wir später behandeln.
Du wirst auch sehen, dass wenn ein Workflow mehrere Ausgaben erzeugt, jede auf diese Weise im `output`-Block aufgeführt wird.

Um mehr zu erfahren, siehe [Publishing outputs](https://nextflow.io/docs/latest/workflow.html#publishing-outputs) in der Nextflow-Referenzdokumentation.

??? info "Ältere Syntax zum Veröffentlichen von Ausgaben mit `publishDir`"

    Bis vor kurzem war die etablierte Methode zum Veröffentlichen von Ausgaben, dies auf der Ebene jedes einzelnen Prozesses mit einer `publishDir`-Direktive zu tun.

    Du wirst dieses Code-Muster immer noch überall in älteren Nextflow-Pipelines und Prozessmodulen finden, daher ist es wichtig, sich dessen bewusst zu sein.

    Anstatt einen `publish:`-Block im Workflow und einen `output`-Block auf oberster Ebene zu haben, würdest du eine `publishDir`-Zeile in der `sayHello`-Prozessdefinition sehen:

    ```groovy title="Syntax example" linenums="1" hl_lines="3"
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

    Wir empfehlen jedoch nicht, dies in neuen Arbeiten zu verwenden, da es in zukünftigen Versionen der Nextflow-Sprache irgendwann nicht mehr erlaubt sein wird.

### Fazit

Du weißt jetzt, wie ein einfacher Nextflow-Workflow strukturiert ist und wie die grundlegenden Komponenten mit seiner Funktionalität zusammenhängen.

### Wie geht es weiter?

Lerne, deine Workflow-Ausführungen bequem zu verwalten.

---

## 4. Workflow-Ausführungen verwalten

Zu wissen, wie man Workflows startet und Ausgaben abruft, ist großartig, aber du wirst schnell feststellen, dass es einige andere Aspekte der Workflow-Verwaltung gibt, die dir das Leben leichter machen werden.

Hier zeigen wir dir, wie du die `resume`-Funktion nutzen kannst, wenn du denselben Workflow erneut starten musst, wie du die Ausführungslogs mit `nextflow log` inspizierst und wie du ältere Work-Verzeichnisse mit `nextflow clean` löschst.

### 4.1. Einen Workflow mit `-resume` erneut starten

Manchmal möchtest du eine Pipeline erneut ausführen, die du bereits zuvor gestartet hast, ohne Arbeit zu wiederholen, die bereits erfolgreich abgeschlossen wurde.

Nextflow hat eine Option namens `-resume`, die dir dies ermöglicht.
Konkret werden in diesem Modus alle Prozesse übersprungen, die bereits mit genau demselben Code, denselben Einstellungen und Eingaben ausgeführt wurden.
Das bedeutet, dass Nextflow nur Prozesse ausführt, die du seit dem letzten Lauf hinzugefügt oder geändert hast, oder denen du neue Einstellungen oder Eingaben zur Verfügung stellst.

Es gibt zwei Hauptvorteile dabei:

- Wenn du gerade eine Pipeline entwickelst, kannst du schneller iterieren, da du nur den/die Prozess(e) ausführen musst, an dem/denen du aktiv arbeitest, um deine Änderungen zu testen.
- Wenn du eine Pipeline in Produktion ausführst und etwas schiefgeht, kannst du in vielen Fällen das Problem beheben und die Pipeline erneut starten, und sie wird ab dem Punkt des Fehlers fortgesetzt, was dir viel Zeit und Rechenleistung sparen kann.

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

Die Konsolenausgabe sollte vertraut aussehen, aber es gibt eine Sache, die etwas anders ist als zuvor.

Suche nach dem `cached:`-Teil, der in der Prozessstatuszeile (Zeile 5) hinzugefügt wurde, was bedeutet, dass Nextflow erkannt hat, dass es diese Arbeit bereits erledigt hat und einfach das Ergebnis aus der vorherigen erfolgreichen Ausführung wiederverwendet hat.

Du kannst auch sehen, dass der Hash des Work-Unterverzeichnisses derselbe ist wie in der vorherigen Ausführung.
Nextflow zeigt dir buchstäblich auf die vorherige Ausführung und sagt "Das habe ich dort drüben schon gemacht."

!!! tip

    Wenn du eine Pipeline mit `resume` erneut ausführst, überschreibt Nextflow keine Dateien, die außerhalb des Work-Verzeichnisses von Ausführungen veröffentlicht wurden, die zuvor erfolgreich ausgeführt wurden.

    Um mehr zu erfahren, siehe [Cache and resume](https://nextflow.io/docs/latest/cache-and-resume.html) in der Nextflow-Referenzdokumentation.

### 4.2. Das Log vergangener Ausführungen inspizieren

Wann immer du einen Nextflow-Workflow startest, wird eine Zeile in eine Logdatei namens `history` geschrieben, die sich in einem versteckten Verzeichnis namens `.nextflow` im aktuellen Arbeitsverzeichnis befindet.

??? abstract "Dateiinhalt"

    ```txt title=".nextflow/history" linenums="1"
    2025-07-04 19:27:09	1.8s	wise_watson	OK	3539118582ccde68dde471cc2c66295c	a02c9c46-c3c7-4085-9139-d1b9b5b194c8	nextflow run 1-hello.nf --input 'Hello World'
    2025-07-04 19:27:20	2.9s	spontaneous_blackwell	OK	3539118582ccde68dde471cc2c66295c	59a5db23-d83c-4c02-a54e-37ddb73a337e	nextflow run 1-hello.nf --input Bonjour
    2025-07-04 19:27:31	1.8s	gigantic_yonath	OK	3539118582ccde68dde471cc2c66295c	5acaa83a-6ad6-4509-bebc-cb25d5d7ddd0	nextflow run 1-hello.nf --input 'Dobry den'
    2025-07-04 19:27:45	2.4s	backstabbing_swartz	OK	3539118582ccde68dde471cc2c66295c	5f4b3269-5b53-404a-956c-cac915fbb74e	nextflow run 1-hello.nf --input Konnichiwa
    2025-07-04 19:27:57	2.1s	goofy_wilson	OK	3539118582ccde68dde471cc2c66295c	5f4b3269-5b53-404a-956c-cac915fbb74e	nextflow run 1-hello.nf --input Konnichiwa -resume
    ```

Diese Datei gibt dir den Zeitstempel, den Run-Namen, den Status, die Revisions-ID, die Session-ID und die vollständige Kommandozeile für jeden Nextflow-Lauf, der aus dem aktuellen Arbeitsverzeichnis gestartet wurde.

Eine bequemere Möglichkeit, auf diese Informationen zuzugreifen, ist die Verwendung des Befehls [`nextflow log`](https://nextflow.io/docs/latest/reference/cli.html#log).

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

Nextflow verwendet die Session-ID, um Run-Caching-Informationen unter dem `cache`-Verzeichnis zu gruppieren, das sich ebenfalls unter `.nextflow` befindet.

### 4.3. Ältere Work-Verzeichnisse löschen

Wenn du viele Pipelines ausführst, kannst du am Ende sehr viele Dateien über viele Unterverzeichnisse hinweg ansammeln.
Da die Unterverzeichnisse zufällig benannt sind, ist es schwierig, anhand ihrer Namen zu erkennen, welche älter vs. neuere Läufe sind.

Glücklicherweise enthält Nextflow einen hilfreichen Befehl namens [`nextflow clean`](https://www.nextflow.io/docs/latest/reference/cli.html#clean), der automatisch die Work-Unterverzeichnisse für vergangene Läufe löschen kann, die dich nicht mehr interessieren.

#### 4.3.1. Löschkriterien festlegen

Es gibt mehrere Optionen, um zu bestimmen, was gelöscht werden soll, die du in der oben verlinkten Dokumentation erkunden kannst.
Hier zeigen wir dir ein Beispiel, das alle Unterverzeichnisse von Läufen vor einem bestimmten Lauf löscht, der über seinen Run-Namen angegeben wird.

Suche den letzten erfolgreichen Lauf, bei dem du `-resume` nicht verwendet hast; in unserem Fall war der Run-Name `backstabbing_swartz`.

Der Run-Name ist die maschinengenerierte zweiteilige Zeichenkette, die in eckigen Klammern in der `Launching (...)`-Konsolenausgabezeile angezeigt wird.
Du kannst auch das Nextflow-Log verwenden, um einen Lauf basierend auf seinem Zeitstempel und/oder seiner Kommandozeile nachzuschlagen.

#### 4.3.2. Einen Probelauf durchführen

Zuerst verwenden wir das Dry-Run-Flag `-n`, um zu überprüfen, was bei dem Befehl gelöscht würde:

```bash
nextflow clean -before backstabbing_swartz -n
```

??? success "Befehlsausgabe"

    ```console
    Would remove /workspaces/training/hello-nextflow/work/eb/1a5de36637b475afd88fca7f79e024
    Would remove /workspaces/training/hello-nextflow/work/6b/19b0e002ea13486d3a0344c336c1d0
    Would remove /workspaces/training/hello-nextflow/work/45/9a6dd7ab771f93003d040956282883
    ```

Deine Ausgabe wird andere Aufgabenverzeichnisnamen haben und möglicherweise eine andere Anzahl von Zeilen, sollte aber ähnlich wie das Beispiel aussehen.

Wenn du keine Zeilen ausgegeben siehst, hast du entweder keinen gültigen Run-Namen angegeben oder es gibt keine vergangenen Läufe zum Löschen. Stelle sicher, dass du `backstabbing_swartz` im Beispielbefehl durch den entsprechenden neuesten Run-Namen in deinem Log ersetzt.

#### 4.3.3. Mit der Löschung fortfahren

Wenn die Ausgabe wie erwartet aussieht und du mit der Löschung fortfahren möchtest, führe den Befehl mit dem `-f`-Flag anstelle von `-n` erneut aus:

```bash
nextflow clean -before backstabbing_swartz -f
```

??? success "Befehlsausgabe"

    ```console
    Removed /workspaces/training/hello-nextflow/work/eb/1a5de36637b475afd88fca7f79e024
    Removed /workspaces/training/hello-nextflow/work/6b/19b0e002ea13486d3a0344c336c1d0
    Removed /workspaces/training/hello-nextflow/work/45/9a6dd7ab771f93003d040956282883
    ```

Die Ausgabe sollte ähnlich wie zuvor sein, aber jetzt mit 'Removed' statt 'Would remove'.
Beachte, dass dies die zweistelligen Unterverzeichnisse (wie `eb/` oben) nicht entfernt, aber deren Inhalt leert.

!!! Warning

    Das Löschen von Work-Unterverzeichnissen aus vergangenen Läufen entfernt sie aus dem Cache von Nextflow und löscht alle Ausgaben, die in diesen Verzeichnissen gespeichert waren.
    Das bedeutet, es bricht die Fähigkeit von Nextflow, die Ausführung fortzusetzen, ohne die entsprechenden Prozesse erneut auszuführen.

    Du bist dafür verantwortlich, alle Ausgaben zu speichern, die dir wichtig sind! Das ist der Hauptgrund, warum wir es vorziehen, den `copy`-Modus anstelle des `symlink`-Modus für die `publish`-Direktive zu verwenden.

### Fazit

Du weißt, wie du eine Pipeline erneut startest, ohne Schritte zu wiederholen, die bereits auf identische Weise ausgeführt wurden, das Ausführungslog inspizierst und den `nextflow clean`-Befehl verwendest, um alte Work-Verzeichnisse aufzuräumen.

### Wie geht es weiter?

Mach eine kleine Pause! Du hast gerade die Bausteine der Nextflow-Syntax und grundlegende Nutzungsanweisungen aufgenommen.

Im nächsten Abschnitt dieses Trainings werden wir uns vier sukzessive realistischere Versionen der Hello World Pipeline ansehen, die demonstrieren werden, wie Nextflow es dir ermöglicht, mehrere Eingaben effizient zu verarbeiten, Workflows auszuführen, die aus mehreren miteinander verbundenen Schritten bestehen, modulare Code-Komponenten zu nutzen und Container für größere Reproduzierbarkeit und Portabilität zu verwenden.

---

## Quiz

<quiz>
In der Konsolenausgabezeile `[a3/7be2fa] SAYHELLO | 1 of 1 ✔`, was repräsentiert `[a3/7be2fa]`?
- [ ] Die Prozessversionsnummer
- [ ] Eine eindeutige Run-Kennung
- [x] Den verkürzten Pfad zum Work-Verzeichnis der Aufgabe
- [ ] Die Prüfsumme der Ausgabedatei

Mehr erfahren: [2.3. Die ursprüngliche Ausgabe und Logs im Verzeichnis `work/` finden](#23-find-the-original-output-and-logs-in-the-work-directory)
</quiz>

<quiz>
Was ist der Zweck der Datei `.command.sh` in einem Aufgabenverzeichnis?
- [ ] Sie speichert die Konfigurationseinstellungen der Aufgabe
- [x] Sie zeigt den tatsächlichen Befehl, der vom Prozess ausgeführt wurde
- [ ] Sie enthält Fehlermeldungen von fehlgeschlagenen Aufgaben
- [ ] Sie listet Eingabedateien auf, die für die Aufgabe bereitgestellt wurden

Mehr erfahren: [2.3. Die ursprüngliche Ausgabe und Logs im Verzeichnis `work/` finden](#23-find-the-original-output-and-logs-in-the-work-directory)
</quiz>

<quiz>
Was passiert mit veröffentlichten Ergebnissen, wenn du einen Workflow ohne `-resume` erneut ausführst?
- [ ] Sie werden in separaten zeitgestempelten Verzeichnissen aufbewahrt
- [x] Sie werden durch die neue Ausführung überschrieben
- [ ] Nextflow verhindert das Überschreiben und schlägt fehl
- [ ] Sie werden automatisch gesichert

Mehr erfahren: [2.4. Den Workflow mit verschiedenen Begrüßungen erneut ausführen](#24-re-run-the-workflow-with-different-greetings)
</quiz>

<quiz>
Was zeigt diese Konsolenausgabe an?

```console
[skipped  ] process > sayHello (1) [100%] 1 of 1, cached: 1 ✔
```

- [ ] Die Aufgabe ist fehlgeschlagen und wurde übersprungen
- [ ] Die Aufgabe wartet in einer Warteschlange
- [x] Nextflow hat Ergebnisse aus einer vorherigen identischen Ausführung wiederverwendet
- [ ] Die Aufgabe wurde manuell abgebrochen

Mehr erfahren: [4.1. Einen Workflow mit `-resume` erneut starten](#41-re-launch-a-workflow-with--resume)
</quiz>

<quiz>
Wo speichert Nextflow den Ausführungsverlauf, den der Befehl `nextflow log` anzeigt?
- [ ] Im Ergebnisverzeichnis
- [ ] Im Work-Verzeichnis
- [x] In der Datei `.nextflow/history`
- [ ] In `nextflow.config`

Mehr erfahren: [4.2. Das Log vergangener Ausführungen inspizieren](#42-inspect-the-log-of-past-executions)
</quiz>

<quiz>
Was ist der Zweck des `params`-Blocks in einer Workflow-Datei?
- [ ] Prozess-Ressourcenanforderungen zu definieren
- [ ] Den Executor zu konfigurieren
- [x] Workflow-Eingabeparameter zu deklarieren und zu typisieren
- [ ] Ausgabeveröffentlichungsoptionen zu spezifizieren

Mehr erfahren: [3.4. Das `params`-System für Kommandozeilenparameter](#34-the-params-system-of-command-line-parameters)
</quiz>

<quiz>
Was bewirkt `mode 'copy'` im `output`-Block des Workflows?
- [ ] Erstellt ein Backup des Work-Verzeichnisses
- [x] Erstellt vollständige Kopien von Dateien anstelle von symbolischen Links
- [ ] Kopiert das Workflow-Skript in die Ergebnisse
- [ ] Aktiviert inkrementelles Dateikopieren

Mehr erfahren: [3.5. Die `publish`-Direktive](#35-the-publish-directive)
</quiz>

<quiz>
Was ist das empfohlene Flag, das mit dem Befehl `nextflow clean` verwendet werden sollte, bevor tatsächlich Dateien gelöscht werden?
- [x] `-n` (Dry-Run), um eine Vorschau dessen zu sehen, was gelöscht würde
- [ ] `-v` (Verbose), um detaillierte Ausgaben zu sehen
- [ ] `-a` (All), um alle Verzeichnisse auszuwählen
- [ ] `-q` (Quiet), um Warnungen zu unterdrücken

Mehr erfahren: [4.3. Ältere Work-Verzeichnisse löschen](#43-delete-older-work-directories)
</quiz>
