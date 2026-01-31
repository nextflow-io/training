# Teil 1: Grundlegende Operationen ausführen

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Im ersten Teil des Nextflow für Bioimaging Trainingskurses verwenden wir ein sehr einfaches, domänenunabhängiges Hello World Beispiel, um wesentliche Operationen zu demonstrieren und auf die entsprechenden Nextflow Code-Komponenten hinzuweisen.

## 1. Den Workflow ausführen

Wir stellen dir ein Workflow-Skript namens `hello-world.nf` zur Verfügung, das eine Eingabe über ein Befehlszeilenargument namens `--greeting` entgegennimmt und eine Textdatei mit dieser Begrüßung erzeugt.
Wir werden uns den Code noch nicht ansehen; zuerst schauen wir uns an, wie es aussieht, wenn man ihn ausführt.

### 1.1. Den Workflow starten und die Ausführung überwachen

Führe im Terminal folgenden Befehl aus:

```bash
nextflow run hello-world.nf --greeting 'Hello World!'
```

Deine Konsolenausgabe sollte in etwa so aussehen:

```console title="Ausgabe" linenums="1"
 N E X T F L O W   ~  version 25.04.3

Launching `hello-world.nf` [goofy_torvalds] DSL2 - revision: c33d41f479

executor >  local (1)
[a3/7be2fa] sayHello | 1 of 1 ✔
```

Glückwunsch, du hast gerade deinen ersten Nextflow Workflow ausgeführt!

Die wichtigste Ausgabe hier ist die letzte Zeile (Zeile 6):

```console title="Ausgabe" linenums="6"
[a3/7be2fa] sayHello | 1 of 1 ✔
```

Dies zeigt uns, dass der `sayHello` Prozess einmal erfolgreich ausgeführt wurde (`1 of 1 ✔`).

Das ist großartig, aber du fragst dich vielleicht: wo ist die Ausgabe?

### 1.2. Die Ausgabedatei im `results` Verzeichnis finden

Dieser Workflow ist so konfiguriert, dass er seine Ausgabe in ein Verzeichnis namens `results` veröffentlicht.
Wenn du dir dein aktuelles Verzeichnis ansiehst, wirst du sehen, dass Nextflow beim Ausführen des Workflows ein neues Verzeichnis namens `results` erstellt hat, das eine Datei namens `output.txt` enthält.

```console title="results/" linenums="1"
results
└── output.txt
```

Öffne die Datei; der Inhalt sollte mit der Begrüßung übereinstimmen, die du auf der Befehlszeile angegeben hast.

<details>
  <summary>Dateiinhalt</summary>

```console title="results/output.txt" linenums="1"
Hello World!
```

</details>

Das ist großartig, unser Workflow hat getan, was er sollte!

Beachte jedoch, dass das 'veröffentlichte' Ergebnis eine Kopie (oder in manchen Fällen ein Symlink) der tatsächlichen Ausgabe ist, die Nextflow bei der Ausführung des Workflows erzeugt hat.

Jetzt werden wir also unter die Haube schauen, um zu sehen, wo Nextflow tatsächlich die Arbeit ausgeführt hat.

!!! warning "Warnung"

    Nicht alle Workflows sind so eingerichtet, dass sie Ausgaben in ein results Verzeichnis veröffentlichen, und/oder der Verzeichnisname kann unterschiedlich sein.
    Etwas weiter unten in diesem Abschnitt zeigen wir dir, wie du herausfindest, wo dieses Verhalten festgelegt ist.

### 1.3. Die ursprüngliche Ausgabe und Logs im `work/` Verzeichnis finden

Wenn du einen Workflow ausführst, erstellt Nextflow für jeden einzelnen Aufruf jedes Prozesses im Workflow (=jeden Schritt in der Pipeline) ein eigenes 'Aufgabenverzeichnis'.
Für jedes wird es die notwendigen Eingaben bereitstellen, die relevanten Anweisung(en) ausführen und Ausgaben sowie Logdateien innerhalb dieses einen Verzeichnisses schreiben, das automatisch mit einem Hash benannt wird, um es eindeutig zu machen.

Alle diese Aufgabenverzeichnisse befinden sich unter einem Verzeichnis namens `work` in deinem aktuellen Verzeichnis (wo du den Befehl ausführst).

Das klingt vielleicht verwirrend, also schauen wir uns an, wie das in der Praxis aussieht.

Zurück zur Konsolenausgabe für den Workflow, den wir vorhin ausgeführt haben, hatten wir diese Zeile:

```console title="Auszug der Befehlsausgabe" linenums="6"
[a3/7be2fa] sayHello | 1 of 1 ✔
```

Siehst du, wie die Zeile mit `[a3/7be2fa]` beginnt?
Das ist eine gekürzte Form des Aufgabenverzeichnispfads für diesen einen Prozessaufruf und zeigt dir, wo du die Ausgabe des `sayHello` Prozessaufrufs innerhalb des `work/` Verzeichnispfads findest.

Du kannst den vollständigen Pfad finden, indem du folgenden Befehl eingibst (ersetze `a3/7be2fa` mit dem, was du in deinem eigenen Terminal siehst) und die Tab-Taste drückst, um den Pfad zu vervollständigen oder ein Sternchen hinzuzufügen:

```bash
tree work/a3/7be2fa*
```

Dies sollte den vollständigen Verzeichnispfad ergeben: `work/a3/7be2fa7be2fad5e71e5f49998f795677fd68`

Schauen wir uns an, was sich darin befindet.

!!! Tip "Tipp"

    Wenn du den Inhalt des Aufgabenunterverzeichnisses im VSCode Datei-Explorer durchsuchst, siehst du alle Dateien sofort.
    Allerdings sind die Logdateien so eingestellt, dass sie im Terminal unsichtbar sind. Wenn du also `ls` oder `tree` verwenden möchtest, um sie anzuzeigen, musst du die entsprechende Option zum Anzeigen unsichtbarer Dateien setzen.

    ```bash
    tree -a work
    ```

Die genauen Unterverzeichnisnamen werden auf deinem System unterschiedlich sein.

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

Du solltest sofort die `output.txt` Datei erkennen, die tatsächlich die ursprüngliche Ausgabe des `sayHello` Prozesses ist, die in das `results` Verzeichnis veröffentlicht wurde.
Wenn du sie öffnest, wirst du die `Hello World!` Begrüßung wieder finden.

<details>
  <summary>Dateiinhalt von output.txt</summary>

```console title="work/a3/7be2fa7be2fad5e71e5f49998f795677fd68/output.txt" linenums="1"
Hello World!
```

</details>

Was ist also mit all diesen anderen Dateien?

Dies sind die Hilfs- und Logdateien, die Nextflow als Teil der Aufgabenausführung geschrieben hat:

- **`.command.begin`**: Sentinel-Datei, die erstellt wird, sobald die Aufgabe gestartet wird.
- **`.command.err`**: Fehlermeldungen (`stderr`), die vom Prozessaufruf ausgegeben wurden
- **`.command.log`**: Vollständige Logausgabe, die vom Prozessaufruf ausgegeben wurde
- **`.command.out`**: Reguläre Ausgabe (`stdout`) des Prozessaufrufs
- **`.command.run`**: Vollständiges Skript, das von Nextflow ausgeführt wurde, um den Prozessaufruf auszuführen
- **`.command.sh`**: Der Befehl, der tatsächlich vom Prozessaufruf ausgeführt wurde
- **`.exitcode`**: Der Exit-Code, der aus dem Befehl resultierte

Die `.command.sh` Datei ist besonders nützlich, weil sie dir den Hauptbefehl zeigt, den Nextflow ausgeführt hat, ohne die gesamte Buchführung und Aufgaben-/Umgebungseinrichtung.

<details>
  <summary>Dateiinhalt</summary>

```console title="work/a3/7be2fa7be2fad5e71e5f49998f795677fd68/.command.sh" linenums="1"
#!/bin/bash -ue
echo 'Hello World!' > output.txt

```

</details>

!!! Tip "Tipp"

    Wenn etwas schiefgeht und du Fehlersuche betreiben musst, kann es nützlich sein, das `command.sh` Skript anzuschauen, um genau zu prüfen, welchen Befehl Nextflow basierend auf den Workflow-Anweisungen, Variableninterpolation usw. zusammengestellt hat.

### 1.4. Optionale Übung: erneut mit verschiedenen Begrüßungen ausführen

Versuche, den Workflow mehrmals mit verschiedenen Werten für das `--greeting` Argument erneut auszuführen, und schaue dir dann sowohl den Inhalt des `results/` Verzeichnisses als auch die Aufgabenverzeichnisse an.

Beobachte, wie die Ausgaben und Logs isolierter Aufgabenverzeichnisse erhalten bleiben, während der Inhalt des `results` Verzeichnisses von der Ausgabe nachfolgender Ausführungen überschrieben wird.

### Zusammenfassung

Du weißt, wie man ein einfaches Nextflow Skript ausführt, seine Ausführung überwacht und seine Ausgaben findet.

### Was kommt als Nächstes?

Lerne, wie man ein grundlegendes Nextflow Skript liest und erkennt, wie seine Komponenten mit seiner Funktionalität zusammenhängen.

---

## 2. Das Hello World Workflow Starter-Skript untersuchen

Was wir dort gemacht haben, war im Grunde, das Workflow-Skript wie eine Black Box zu behandeln.
Jetzt, da wir gesehen haben, was es tut, lass uns die Box öffnen und hineinschauen.

_Das Ziel hier ist nicht, die Syntax von Nextflow Code auswendig zu lernen, sondern eine grundlegende Intuition dafür zu entwickeln, was die Hauptkomponenten sind und wie sie organisiert sind._

### 2.1. Die gesamte Code-Struktur untersuchen

Lass uns das `hello-world.nf` Skript im Editor-Fenster öffnen.

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

    // gib eine Begrüßung aus
    sayHello(params.greeting)
}
```

</details>

Ein Nextflow Skript umfasst zwei Haupttypen von Kernkomponenten: einen oder mehrere **Prozesse** und den **Workflow** selbst.
Jeder **Prozess** beschreibt, welche Operation(en) der entsprechende Schritt in der Pipeline durchführen soll, während der **Workflow** die Datenflusslogik beschreibt, die die verschiedenen Schritte verbindet.

Schauen wir uns zuerst den **Prozess**-Block genauer an, dann betrachten wir den **Workflow**-Block.

### 2.2. Die `process` Definition

Der erste Codeblock beschreibt einen **Prozess**.
Die Prozessdefinition beginnt mit dem Schlüsselwort `process`, gefolgt vom Prozessnamen und schließlich dem Prozesskörper, der durch geschweifte Klammern begrenzt wird.
Der Prozesskörper muss einen script-Block enthalten, der den auszuführenden Befehl angibt, was alles sein kann, was du in einem Befehlszeilenterminal ausführen könntest.

Hier haben wir einen **Prozess** namens `sayHello`, der eine **Eingabe**-Variable namens `greeting` nimmt und seine **Ausgabe** in eine Datei namens `output.txt` schreibt.

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

Dies ist eine sehr minimale Prozessdefinition, die nur eine `input` Definition, eine `output` Definition und das `script` zur Ausführung enthält.

Die `input` Definition enthält den `val` Qualifier, der Nextflow mitteilt, dass ein Wert irgendeiner Art erwartet wird (kann ein String, eine Zahl, was auch immer sein).

Die `output` Definition enthält den `path` Qualifier, der Nextflow mitteilt, dass dies als Pfad behandelt werden soll (umfasst sowohl Verzeichnispfade als auch Dateien).

!!! Tip "Tipp"

    Die Ausgabedefinition _bestimmt_ nicht, welche Ausgabe erstellt wird.
    Sie _deklariert_ einfach, wo die erwarteten Ausgabedatei(en) zu finden sind, damit Nextflow danach suchen kann, sobald die Ausführung abgeschlossen ist.

    Dies ist notwendig, um zu überprüfen, dass der Befehl erfolgreich ausgeführt wurde, und um die Ausgabe bei Bedarf an nachgelagerte Prozesse weiterzugeben.
    Erzeugte Ausgabe, die nicht mit dem übereinstimmt, was im output-Block deklariert ist, wird nicht an nachgelagerte Prozesse weitergegeben.

In einer realen Pipeline enthält ein Prozess normalerweise zusätzliche Informationen wie Prozessdirektiven, die wir in Kürze einführen werden.

### 2.3. Die `workflow` Definition

Der zweite Codeblock beschreibt den **Workflow** selbst.
Die Workflow-Definition beginnt mit dem Schlüsselwort `workflow`, gefolgt von einem optionalen Namen, dann dem Workflow-Körper, der durch geschweifte Klammern begrenzt wird.

Hier haben wir einen **Workflow**, der aus einem Aufruf des `sayHello` Prozesses besteht, der eine Eingabe `params.greeting` nimmt, die den Wert enthält, den wir dem `--greeting` Parameter gegeben haben.

```groovy title="hello-world.nf" linenums="22"
workflow {

    // gib eine Begrüßung aus
    sayHello(params.greeting)
}
```

Dies ist eine sehr minimale **Workflow**-Definition.
In einer realen Pipeline enthält der Workflow typischerweise mehrere Aufrufe von **Prozessen**, die durch **Channels** verbunden sind, und es können Standardwerte für die variablen Eingaben eingerichtet sein.

Wir werden dies in Aktion sehen, wenn wir nf-core/molkart in Teil 2 des Kurses ausführen.

### 2.4. Das `params` System von Befehlszeilenparametern

Das `params.greeting`, das wir dem `sayHello()` Prozessaufruf übergeben, ist ein cleveres Stück Nextflow Code und ist es wert, eine zusätzliche Minute darauf zu verwenden.

Wie oben erwähnt, geben wir so den Wert des `--greeting` Befehlszeilenparameters an den `sayHello()` Prozessaufruf weiter.
Tatsächlich ermöglicht die einfache Deklaration von `params.someParameterName` uns, dem Workflow einen Parameter namens `--someParameterName` von der Befehlszeile aus zu geben.

!!! Tip "Tipp"

    Diese Workflow-Parameter, die mit dem `params` System deklariert werden, verwenden immer zwei Bindestriche (`--`).
    Dies unterscheidet sie von Nextflow-Parametern, die nur einen Bindestrich (`-`) verwenden.

### Zusammenfassung

Du weißt jetzt, wie ein einfacher Nextflow Workflow strukturiert ist und wie die grundlegenden Komponenten mit seiner Funktionalität zusammenhängen.

### Was kommt als Nächstes?

Lerne, deine Workflow-Ausführungen bequem zu verwalten.

---

## 3. Workflow-Ausführungen verwalten

Zu wissen, wie man Workflows startet und Ausgaben abruft, ist großartig, aber du wirst schnell feststellen, dass es ein paar andere Aspekte der Workflow-Verwaltung gibt, die dir das Leben erleichtern werden.

Hier zeigen wir dir, wie du die `resume` Funktion nutzen kannst, wenn du denselben Workflow erneut starten musst, wie du die Ausführungslogs mit `nextflow log` einsehen kannst und wie du ältere work Verzeichnisse mit `nextflow clean` löschen kannst.

### 3.1. Einen Workflow mit `-resume` erneut starten

Manchmal wirst du eine Pipeline erneut ausführen wollen, die du bereits zuvor gestartet hast, ohne bereits erfolgreich abgeschlossene Arbeiten erneut durchzuführen.

Nextflow hat eine Option namens `-resume`, die dir dies ermöglicht.
Konkret werden in diesem Modus alle Prozesse übersprungen, die bereits mit exakt demselben Code, denselben Einstellungen und Eingaben ausgeführt wurden.
Das bedeutet, dass Nextflow nur Prozesse ausführt, die du seit dem letzten Lauf hinzugefügt oder geändert hast, oder denen du neue Einstellungen oder Eingaben zur Verfügung stellst.

Dies hat zwei Hauptvorteile:

- Wenn du gerade eine Pipeline entwickelst, kannst du schneller iterieren, da du nur den/die Prozess(e) ausführen musst, an dem/denen du aktiv arbeitest, um deine Änderungen zu testen.
- Wenn du eine Pipeline im Produktionsbetrieb ausführst und etwas schiefgeht, kannst du in vielen Fällen das Problem beheben und die Pipeline erneut starten, und sie wird vom Fehlerpunkt aus weiterlaufen, was dir viel Zeit und Rechenleistung sparen kann.

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

Achte auf den `cached:` Teil, der in der Prozessstatuszeile (Zeile 5) hinzugefügt wurde, was bedeutet, dass Nextflow erkannt hat, dass es diese Arbeit bereits erledigt hat und einfach das Ergebnis aus dem vorherigen erfolgreichen Lauf wiederverwendet.

Du kannst auch sehen, dass der work Unterverzeichnis-Hash derselbe ist wie im vorherigen Lauf.
Nextflow zeigt dir buchstäblich auf die vorherige Ausführung und sagt "Das habe ich bereits dort drüben gemacht."

!!! Tip "Tipp"

    Wenn du eine Pipeline mit `resume` erneut ausführst, überschreibt Nextflow keine Dateien, die von einem Prozessaufruf, der zuvor erfolgreich ausgeführt wurde, in ein `publishDir` Verzeichnis geschrieben wurden.

### 3.2. Das Log vergangener Ausführungen einsehen

Jedes Mal, wenn du einen Nextflow Workflow startest, wird eine Zeile in eine Logdatei namens `history` geschrieben, unter einem versteckten Verzeichnis namens `.nextflow` im aktuellen Arbeitsverzeichnis.

Eine bequemere Möglichkeit, auf diese Informationen zuzugreifen, ist die Verwendung des `nextflow log` Befehls.

```bash
nextflow log
```

Dies gibt den Inhalt der Logdatei im Terminal aus und zeigt dir den Zeitstempel, Laufnamen, Status und die vollständige Befehlszeile für jeden Nextflow Lauf, der aus dem aktuellen Arbeitsverzeichnis heraus gestartet wurde.

### 3.3. Ältere work Verzeichnisse löschen

Während des Entwicklungsprozesses führst du deine Entwurfs-Pipelines typischerweise sehr oft aus, was zu einer Ansammlung von sehr vielen Dateien über viele Unterverzeichnisse führen kann.
Da die Unterverzeichnisse zufällig benannt sind, ist es schwierig, anhand ihrer Namen zu erkennen, welche älter und welche aktueller sind.

Nextflow enthält einen praktischen `clean` Unterbefehl, der automatisch die work Unterverzeichnisse für vergangene Läufe löschen kann, die dich nicht mehr interessieren, mit mehreren [Optionen](https://www.nextflow.io/docs/latest/reference/cli.html#clean), um zu steuern, was gelöscht wird.

Du kannst das Nextflow Log verwenden, um einen Lauf anhand seines Zeitstempels und/oder seiner Befehlszeile nachzuschlagen, und dann `nextflow clean -before <run_name> -f` verwenden, um work Verzeichnisse von früheren Läufen zu löschen.

!!! Warning "Warnung"

    Das Löschen von work Unterverzeichnissen aus vergangenen Läufen entfernt sie aus dem Nextflow Cache und löscht alle Ausgaben, die in diesen Verzeichnissen gespeichert waren.
    Das bedeutet, dass es Nextflows Fähigkeit zerstört, die Ausführung fortzusetzen, ohne die entsprechenden Prozesse erneut auszuführen.

    Du bist dafür verantwortlich, alle Ausgaben zu speichern, die dir wichtig sind oder auf die du dich verlassen möchtest! Wenn du die `publishDir` Direktive für diesen Zweck verwendest, stelle sicher, dass du den `copy` Modus verwendest, nicht den `symlink` Modus.

### Zusammenfassung

Du weißt, wie man eine Pipeline erneut startet, ohne Schritte zu wiederholen, die bereits auf identische Weise ausgeführt wurden, das Ausführungslog einsieht und den `nextflow clean` Befehl verwendest, um alte work Verzeichnisse aufzuräumen.

### Was kommt als Nächstes?

Jetzt, da du grundlegende Nextflow Operationen verstehst, bist du bereit, eine echte Bioimaging Pipeline mit nf-core/molkart auszuführen.
