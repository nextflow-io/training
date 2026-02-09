# Teil 1: Grundlegende Operationen ausführen

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Im ersten Teil des Nextflow for Bioimaging Trainings verwenden wir ein sehr einfaches, domänenunabhängiges Hello World Beispiel, um grundlegende Operationen zu demonstrieren und die entsprechenden Nextflow-Code-Komponenten zu zeigen.

## 1. Den Workflow ausführen

Wir stellen dir ein Workflow-Skript namens `hello-world.nf` zur Verfügung, das eine Eingabe über ein Kommandozeilenargument namens `--greeting` entgegennimmt und eine Textdatei mit dieser Begrüßung erzeugt.
Wir schauen uns den Code noch nicht an; zuerst sehen wir uns an, wie die Ausführung aussieht.

### 1.1. Den Workflow starten und die Ausführung überwachen

Führe im Terminal folgenden Befehl aus:

```bash
nextflow run hello-world.nf --greeting 'Hello World!'
```

Deine Konsolenausgabe sollte ungefähr so aussehen:

```console title="Output" linenums="1"
 N E X T F L O W   ~  version 25.04.3

Launching `hello-world.nf` [goofy_torvalds] DSL2 - revision: c33d41f479

executor >  local (1)
[a3/7be2fa] sayHello | 1 of 1 ✔
```

Glückwunsch, du hast gerade deinen ersten Nextflow-Workflow ausgeführt!

Die wichtigste Ausgabe ist hier die letzte Zeile (Zeile 6):

```console title="Output" linenums="6"
[a3/7be2fa] sayHello | 1 of 1 ✔
```

Das zeigt uns, dass der `sayHello`-Prozess erfolgreich einmal ausgeführt wurde (`1 of 1 ✔`).

Das ist großartig, aber du fragst dich vielleicht: Wo ist die Ausgabe?

### 1.2. Die Ausgabedatei im `results`-Verzeichnis finden

Dieser Workflow ist so konfiguriert, dass er seine Ausgabe in einem Verzeichnis namens `results` veröffentlicht.
Wenn du dir dein aktuelles Verzeichnis ansiehst, wirst du sehen, dass Nextflow beim Ausführen des Workflows ein neues Verzeichnis namens `results` erstellt hat, das eine Datei namens `output.txt` enthält.

```console title="results/" linenums="1"
results
└── output.txt
```

Öffne die Datei; der Inhalt sollte mit der Begrüßung übereinstimmen, die du in der Kommandozeile angegeben hast.

<details>
  <summary>Dateiinhalt</summary>

```console title="results/output.txt" linenums="1"
Hello World!
```

</details>

Super, unser Workflow hat getan, was er sollte!

Beachte jedoch, dass das 'veröffentlichte' Ergebnis eine Kopie (oder in manchen Fällen ein Symlink) der tatsächlichen Ausgabe ist, die Nextflow bei der Ausführung des Workflows erzeugt hat.

Jetzt schauen wir unter die Haube, um zu sehen, wo Nextflow die Arbeit tatsächlich ausgeführt hat.

!!! warning "Warnung"

    Nicht alle Workflows sind so eingerichtet, dass sie Ausgaben in einem results-Verzeichnis veröffentlichen, und/oder der Verzeichnisname kann anders sein.
    Etwas weiter unten in diesem Abschnitt zeigen wir dir, wie du herausfindest, wo dieses Verhalten festgelegt ist.

### 1.3. Die ursprüngliche Ausgabe und Logs im `work/`-Verzeichnis finden

Wenn du einen Workflow ausführst, erstellt Nextflow für jeden einzelnen Aufruf jedes Prozesses im Workflow (=jeden Schritt in der Pipeline) ein eigenes 'Aufgabenverzeichnis'.
Für jedes dieser Verzeichnisse werden die notwendigen Eingaben bereitgestellt, die relevanten Anweisungen ausgeführt und Ausgaben sowie Logdateien in diesem einen Verzeichnis geschrieben, das automatisch mit einem Hash benannt wird, um es eindeutig zu machen.

Alle diese Aufgabenverzeichnisse befinden sich in einem Verzeichnis namens `work` innerhalb deines aktuellen Verzeichnisses (wo du den Befehl ausführst).

Das klingt vielleicht verwirrend, also schauen wir uns an, wie das in der Praxis aussieht.

Zurück zur Konsolenausgabe des Workflows, den wir vorhin ausgeführt haben, hatten wir diese Zeile:

```console title="Auszug der Befehlsausgabe" linenums="6"
[a3/7be2fa] sayHello | 1 of 1 ✔
```

Siehst du, wie die Zeile mit `[a3/7be2fa]` beginnt?
Das ist eine gekürzte Form des Aufgabenverzeichnispfads für diesen einen Prozessaufruf und zeigt dir, wo du die Ausgabe des `sayHello`-Prozessaufrufs innerhalb des `work/`-Verzeichnispfads findest.

Du kannst den vollständigen Pfad finden, indem du folgenden Befehl eingibst (ersetze `a3/7be2fa` mit dem, was du in deinem eigenen Terminal siehst) und die Tab-Taste drückst, um den Pfad zu vervollständigen, oder ein Sternchen hinzufügst:

```bash
tree work/a3/7be2fa*
```

Das sollte den vollständigen Verzeichnispfad ergeben: `work/a3/7be2fa7be2fad5e71e5f49998f795677fd68`

Schauen wir uns an, was da drin ist.

!!! Tip "Tipp"

    Wenn du den Inhalt des Aufgaben-Unterverzeichnisses im VSCode-Datei-Explorer durchsuchst, siehst du alle Dateien sofort.
    Die Logdateien sind jedoch im Terminal als unsichtbar eingestellt. Wenn du also `ls` oder `tree` verwenden möchtest, um sie anzuzeigen, musst du die entsprechende Option zum Anzeigen unsichtbarer Dateien setzen.

    ```bash
    tree -a work
    ```

Die genauen Unterverzeichnisnamen werden auf deinem System anders sein.

<details>
  <summary>Verzeichnisinhalt</summary>

```console title="work/"
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

</details>

Du solltest sofort die `output.txt`-Datei erkennen, die tatsächlich die ursprüngliche Ausgabe des `sayHello`-Prozesses ist, die in das `results`-Verzeichnis veröffentlicht wurde.
Wenn du sie öffnest, findest du wieder die `Hello World!`-Begrüßung.

<details>
  <summary>Dateiinhalt von output.txt</summary>

```console title="work/a3/7be2fa7be2fad5e71e5f49998f795677fd68/output.txt" linenums="1"
Hello World!
```

</details>

Was ist also mit all den anderen Dateien?

Das sind die Hilfs- und Logdateien, die Nextflow als Teil der Aufgabenausführung geschrieben hat:

- **`.command.begin`**: Sentinel-Datei, die erstellt wird, sobald die Aufgabe gestartet wird.
- **`.command.err`**: Fehlermeldungen (`stderr`), die vom Prozessaufruf ausgegeben wurden
- **`.command.log`**: Vollständige Logausgabe des Prozessaufrufs
- **`.command.out`**: Reguläre Ausgabe (`stdout`) des Prozessaufrufs
- **`.command.run`**: Vollständiges Skript, das von Nextflow ausgeführt wurde, um den Prozessaufruf auszuführen
- **`.command.sh`**: Der Befehl, der tatsächlich vom Prozessaufruf ausgeführt wurde
- **`.exitcode`**: Der Exit-Code, der aus dem Befehl resultierte

Die `.command.sh`-Datei ist besonders nützlich, weil sie dir den Hauptbefehl zeigt, den Nextflow ausgeführt hat, ohne die gesamte Buchhaltung und Aufgaben-/Umgebungseinrichtung.

<details>
  <summary>Dateiinhalt</summary>

```console title="work/a3/7be2fa7be2fad5e71e5f49998f795677fd68/.command.sh" linenums="1"
#!/bin/bash -ue
echo 'Hello World!' > output.txt

```

</details>

!!! Tip "Tipp"

    Wenn etwas schiefgeht und du herausfinden musst, was passiert ist, kann es nützlich sein, das `command.sh`-Skript anzuschauen, um genau zu prüfen, welchen Befehl Nextflow basierend auf den Workflow-Anweisungen, Variableninterpolation usw. zusammengestellt hat.

### 1.4. Optionale Übung: Mit verschiedenen Begrüßungen erneut ausführen

Versuche, den Workflow mehrmals mit verschiedenen Werten für das `--greeting`-Argument auszuführen, und schaue dir dann sowohl den Inhalt des `results/`-Verzeichnisses als auch die Aufgabenverzeichnisse an.

Beobachte, wie die Ausgaben und Logs der isolierten Aufgabenverzeichnisse erhalten bleiben, während der Inhalt des `results`-Verzeichnisses durch die Ausgabe nachfolgender Ausführungen überschrieben wird.

### Fazit

Du weißt, wie du ein einfaches Nextflow-Skript ausführst, seine Ausführung überwachst und seine Ausgaben findest.

### Wie geht es weiter?

Lerne, wie du ein einfaches Nextflow-Skript liest und erkennst, wie seine Komponenten mit seiner Funktionalität zusammenhängen.

---

## 2. Das Hello World Workflow-Starterskript untersuchen

Was wir dort gemacht haben, war im Grunde, das Workflow-Skript wie eine Black Box zu behandeln.
Jetzt, wo wir gesehen haben, was es tut, öffnen wir die Box und schauen hinein.

_Das Ziel hier ist nicht, die Syntax von Nextflow-Code auswendig zu lernen, sondern ein grundlegendes Verständnis dafür zu entwickeln, was die Hauptkomponenten sind und wie sie organisiert sind._

### 2.1. Die allgemeine Code-Struktur untersuchen

Öffnen wir das `hello-world.nf`-Skript im Editor-Bereich.

<details>
  <summary>Code</summary>

```groovy title="hello-world.nf" linenums="1"
#!/usr/bin/env nextflow

/*
 * Verwende echo, um eine Begrüßung in eine Datei zu schreiben
 */
process sayHello {

    publishDir 'results', mode: 'copy'

    input:
    val greeting

    output:
    path 'output.txt'

    script:
    """
    echo '$greeting' > output.txt
    """
}

workflow {

    // Eine Begrüßung ausgeben
    sayHello(params.greeting)
}
```

</details>

Ein Nextflow-Skript besteht aus zwei Haupttypen von Kernkomponenten: einem oder mehreren **Prozessen** und dem **Workflow** selbst.
Jeder **Prozess** beschreibt, welche Operation(en) der entsprechende Schritt in der Pipeline ausführen soll, während der **Workflow** die Datenflusslogik beschreibt, die die verschiedenen Schritte verbindet.

Schauen wir uns zuerst den **Prozess**-Block genauer an, dann den **Workflow**-Block.

### 2.2. Die `process`-Definition

Der erste Codeblock beschreibt einen **Prozess**.
Die Prozessdefinition beginnt mit dem Schlüsselwort `process`, gefolgt vom Prozessnamen und schließlich dem Prozesskörper, der durch geschweifte Klammern begrenzt wird.
Der Prozesskörper muss einen script-Block enthalten, der den auszuführenden Befehl angibt, was alles sein kann, was du in einem Kommandozeilen-Terminal ausführen könntest.

Hier haben wir einen **Prozess** namens `sayHello`, der eine **Eingabe**-Variable namens `greeting` entgegennimmt und seine **Ausgabe** in eine Datei namens `output.txt` schreibt.

<details>
  <summary>Code</summary>

```groovy title="hello-world.nf" linenums="3"
/*
 * Verwende echo, um eine Begrüßung in eine Datei zu schreiben
 */
process sayHello {

    publishDir 'results', mode: 'copy'

    input:
    val greeting

    output:
    path 'output.txt'

    script:
    """
    echo '$greeting' > output.txt
    """
}
```

</details>

Das ist eine sehr minimale Prozessdefinition, die nur eine `input`-Definition, eine `output`-Definition und das auszuführende `script` enthält.

Die `input`-Definition enthält den `val`-Qualifier, der Nextflow mitteilt, dass ein Wert irgendeiner Art erwartet wird (kann ein String, eine Zahl oder was auch immer sein).

Die `output`-Definition enthält den `path`-Qualifier, der Nextflow mitteilt, dass dies als Pfad behandelt werden soll (umfasst sowohl Verzeichnispfade als auch Dateien).

!!! Tip "Tipp"

    Die output-Definition _bestimmt_ nicht, welche Ausgabe erstellt wird.
    Sie _deklariert_ lediglich, wo die erwartete(n) Ausgabedatei(en) zu finden ist/sind, damit Nextflow danach suchen kann, sobald die Ausführung abgeschlossen ist.

    Dies ist notwendig, um zu überprüfen, dass der Befehl erfolgreich ausgeführt wurde, und um die Ausgabe bei Bedarf an nachgelagerte Prozesse weiterzugeben.
    Erzeugte Ausgaben, die nicht mit dem übereinstimmen, was im output-Block deklariert ist, werden nicht an nachgelagerte Prozesse weitergegeben.

In einer realen Pipeline enthält ein Prozess normalerweise zusätzliche Informationen wie Prozessdirektiven, die wir gleich einführen werden.

### 2.3. Die `workflow`-Definition

Der zweite Codeblock beschreibt den **Workflow** selbst.
Die Workflow-Definition beginnt mit dem Schlüsselwort `workflow`, gefolgt von einem optionalen Namen, dann dem Workflow-Körper, der durch geschweifte Klammern begrenzt wird.

Hier haben wir einen **Workflow**, der aus einem Aufruf des `sayHello`-Prozesses besteht, der eine Eingabe `params.greeting` entgegennimmt, die den Wert enthält, den wir dem `--greeting`-Parameter gegeben haben.

```groovy title="hello-world.nf" linenums="22"
workflow {

    // Eine Begrüßung ausgeben
    sayHello(params.greeting)
}
```

Das ist eine sehr minimale **Workflow**-Definition.
In einer realen Pipeline enthält der Workflow typischerweise mehrere Aufrufe von **Prozessen**, die durch **Kanäle** verbunden sind, und es können Standardwerte für die variablen Eingaben eingerichtet sein.

Das werden wir in Aktion sehen, wenn wir nf-core/molkart in Teil 2 des Kurses ausführen.

### 2.4. Das `params`-System für Kommandozeilenparameter

Das `params.greeting`, das wir dem `sayHello()`-Prozessaufruf übergeben, ist ein cleveres Stück Nextflow-Code und ist es wert, eine zusätzliche Minute darauf zu verwenden.

Wie oben erwähnt, übergeben wir damit den Wert des `--greeting`-Kommandozeilenparameters an den `sayHello()`-Prozessaufruf.
Tatsächlich ermöglicht uns die einfache Deklaration von `params.someParameterName`, dem Workflow einen Parameter namens `--someParameterName` von der Kommandozeile aus zu geben.

!!! Tip "Tipp"

    Diese Workflow-Parameter, die mit dem `params`-System deklariert werden, nehmen immer zwei Bindestriche (`--`).
    Das unterscheidet sie von Nextflow-Parametern, die nur einen Bindestrich (`-`) nehmen.

### Fazit

Du weißt jetzt, wie ein einfacher Nextflow-Workflow strukturiert ist und wie die grundlegenden Komponenten mit seiner Funktionalität zusammenhängen.

### Wie geht es weiter?

Lerne, deine Workflow-Ausführungen bequem zu verwalten.

---

## 3. Workflow-Ausführungen verwalten

Zu wissen, wie man Workflows startet und Ausgaben abruft, ist großartig, aber du wirst schnell feststellen, dass es ein paar andere Aspekte der Workflow-Verwaltung gibt, die dir das Leben leichter machen.

Hier zeigen wir dir, wie du die `resume`-Funktion nutzen kannst, wenn du denselben Workflow erneut starten musst, wie du die Ausführungslogs mit `nextflow log` inspizierst und wie du ältere work-Verzeichnisse mit `nextflow clean` löschst.

### 3.1. Einen Workflow mit `-resume` erneut starten

Manchmal möchtest du eine Pipeline erneut ausführen, die du bereits zuvor gestartet hast, ohne Arbeit zu wiederholen, die bereits erfolgreich abgeschlossen wurde.

Nextflow hat eine Option namens `-resume`, die dir das ermöglicht.
Konkret werden in diesem Modus alle Prozesse übersprungen, die bereits mit genau demselben Code, denselben Einstellungen und Eingaben ausgeführt wurden.
Das bedeutet, dass Nextflow nur Prozesse ausführt, die du seit dem letzten Lauf hinzugefügt oder geändert hast, oder denen du neue Einstellungen oder Eingaben gibst.

Das hat zwei wesentliche Vorteile:

- Wenn du gerade eine Pipeline entwickelst, kannst du schneller iterieren, da du nur den/die Prozess(e) ausführen musst, an dem/denen du aktiv arbeitest, um deine Änderungen zu testen.
- Wenn du eine Pipeline in Produktion ausführst und etwas schiefgeht, kannst du in vielen Fällen das Problem beheben und die Pipeline erneut starten, und sie wird ab dem Fehlerpunkt fortgesetzt, was dir viel Zeit und Rechenleistung sparen kann.

Um es zu verwenden, füge einfach `-resume` zu deinem Befehl hinzu und führe ihn aus:

```bash
nextflow run hello-world.nf --greeting 'Hello World!' -resume
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.04.3

    Launching `hello-world.nf` [tiny_noyce] DSL2 - revision: c33d41f479

    [a3/7be2fa] process > sayHello [100%] 1 of 1, cached: 1 ✔
    ```

Achte auf das `cached:`-Bit, das in der Prozessstatuszeile (Zeile 5) hinzugefügt wurde, was bedeutet, dass Nextflow erkannt hat, dass es diese Arbeit bereits erledigt hat und einfach das Ergebnis des vorherigen erfolgreichen Laufs wiederverwendet.

Du kannst auch sehen, dass der work-Unterverzeichnis-Hash derselbe ist wie beim vorherigen Lauf.
Nextflow zeigt dir buchstäblich auf die vorherige Ausführung und sagt "Das habe ich dort drüben schon gemacht."

!!! Tip "Tipp"

    Wenn du eine Pipeline mit `resume` erneut ausführst, überschreibt Nextflow keine Dateien, die von einem Prozessaufruf, der zuvor erfolgreich ausgeführt wurde, in ein `publishDir`-Verzeichnis geschrieben wurden.

### 3.2. Das Log vergangener Ausführungen inspizieren

Wann immer du einen Nextflow-Workflow startest, wird eine Zeile in eine Logdatei namens `history` geschrieben, die sich in einem versteckten Verzeichnis namens `.nextflow` im aktuellen Arbeitsverzeichnis befindet.

Eine bequemere Möglichkeit, auf diese Informationen zuzugreifen, ist die Verwendung des `nextflow log`-Befehls.

```bash
nextflow log
```

Dies gibt den Inhalt der Logdatei im Terminal aus und zeigt dir den Zeitstempel, den Laufnamen, den Status und die vollständige Kommandozeile für jeden Nextflow-Lauf, der vom aktuellen Arbeitsverzeichnis aus gestartet wurde.

### 3.3. Ältere work-Verzeichnisse löschen

Während des Entwicklungsprozesses führst du deine Entwurfs-Pipelines typischerweise sehr oft aus, was zu einer Ansammlung sehr vieler Dateien über viele Unterverzeichnisse hinweg führen kann.
Da die Unterverzeichnisse zufällig benannt werden, ist es schwierig, anhand ihrer Namen zu erkennen, welche älter und welche aktueller sind.

Nextflow enthält einen praktischen `clean`-Unterbefehl, der automatisch die work-Unterverzeichnisse für vergangene Läufe löschen kann, die dich nicht mehr interessieren, mit mehreren [Optionen](https://www.nextflow.io/docs/latest/reference/cli.html#clean), um zu steuern, was gelöscht wird.

Du kannst das Nextflow-Log verwenden, um einen Lauf anhand seines Zeitstempels und/oder seiner Kommandozeile nachzuschlagen, und dann `nextflow clean -before <run_name> -f` verwenden, um work-Verzeichnisse von früheren Läufen zu löschen.

!!! Warning "Warnung"

    Das Löschen von work-Unterverzeichnissen vergangener Läufe entfernt sie aus dem Cache von Nextflow und löscht alle Ausgaben, die in diesen Verzeichnissen gespeichert waren.
    Das bedeutet, es bricht die Fähigkeit von Nextflow, die Ausführung fortzusetzen, ohne die entsprechenden Prozesse erneut auszuführen.

    Du bist dafür verantwortlich, alle Ausgaben zu speichern, die dir wichtig sind oder auf die du dich verlassen möchtest! Wenn du die `publishDir`-Direktive für diesen Zweck verwendest, stelle sicher, dass du den `copy`-Modus verwendest, nicht den `symlink`-Modus.

### Fazit

Du weißt, wie du eine Pipeline erneut startest, ohne Schritte zu wiederholen, die bereits auf identische Weise ausgeführt wurden, das Ausführungslog inspizierst und den `nextflow clean`-Befehl verwendest, um alte work-Verzeichnisse aufzuräumen.

### Wie geht es weiter?

Jetzt, wo du grundlegende Nextflow-Operationen verstehst, bist du bereit, eine echte Bioimaging-Pipeline mit nf-core/molkart auszuführen.
