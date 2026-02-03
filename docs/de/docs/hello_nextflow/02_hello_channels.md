# Teil 2: Hello Channels

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<!--
<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/lJ41WMMm44M?si=xCItHLiOQWqoqBB9&amp;list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } Siehe [die gesamte Playlist](https://www.youtube.com/playlist?list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik) auf dem Nextflow YouTube-Kanal.

:green_book: Das Videotranskript ist [hier](./transcripts/02_hello_channels.md) verfügbar.
///
-->

In Teil 1 dieses Kurses (Hello World) haben wir dir gezeigt, wie du eine variable Eingabe an einen process übergibst, indem du sie direkt im process-Aufruf angibst: `sayHello(params.input)`.
Das war ein bewusst vereinfachter Ansatz.
In der Praxis hat dieser Ansatz wesentliche Einschränkungen; er funktioniert nämlich nur für sehr einfache Fälle, bei denen wir den process nur einmal mit einem einzelnen Wert ausführen möchten.
In den meisten realistischen Workflow-Anwendungsfällen wollen wir mehrere Werte verarbeiten (z.B. experimentelle Daten für mehrere Proben), daher brauchen wir einen ausgefeilteren Weg, um Eingaben zu handhaben.

Dafür sind Nextflow [**Kanäle**](https://nextflow.io/docs/latest/channel.html) da.
Kanäle sind Warteschlangen, die entwickelt wurden, um Eingaben effizient zu handhaben und sie von einem Schritt zum nächsten in mehrstufigen Workflows zu transportieren, während sie eingebaute Parallelität und viele zusätzliche Vorteile bieten.

In diesem Teil des Kurses lernst du, wie du einen Kanal verwendest, um mehrere Eingaben aus verschiedenen Quellen zu handhaben.
Du lernst auch, [**Operatoren**](https://nextflow.io/docs/latest/reference/operator.html) zu verwenden, um Kanal-Inhalte nach Bedarf zu transformieren.

??? info "Wie du von diesem Abschnitt aus beginnen kannst"

    Dieser Abschnitt des Kurses setzt voraus, dass du Teil 1 des [Hello Nextflow](./index.md)-Kurses abgeschlossen hast, aber wenn du mit den dort behandelten Grundlagen vertraut bist, kannst du von hier aus starten, ohne etwas Besonderes tun zu müssen.

---

## 0. Aufwärmen: `hello-channels.nf` ausführen

Wir werden das Workflow-Skript `hello-channels.nf` als Ausgangspunkt verwenden.
Es entspricht dem Skript, das durch das Durcharbeiten von Teil 1 dieses Trainingskurses erstellt wurde, außer dass wir das Ausgabeziel geändert haben:

```groovy title="hello-channels.nf" linenums="37" hl_lines="3"
output {
    first_output {
        path 'hello_channels'
        mode 'copy'
    }
}
```

Um sicherzustellen, dass alles funktioniert, führe das Skript einmal aus, bevor du Änderungen vornimmst:

```bash
nextflow run hello-channels.nf --input 'Hello Channels!'
```

??? success "Befehlsausgabe"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [wise_jennings] DSL2 - revision: b24f4902d6

    executor >  local (1)
    [6f/824bc1] process > sayHello [100%] 1 of 1 ✔
    ```

Wie zuvor findest du die Ausgabedatei `output.txt` im Verzeichnis `results/hello_channels` (wie im `output`-Block des Workflow-Skripts oben angegeben).

??? abstract "Verzeichnisinhalte"

    ```console title="results/hello_channels" hl_lines="2-3"
    results
    ├── hello_channels
    │   └── output.txt
    ├── hello_world
    │   └── output.txt
    └── output.txt -> /workspaces/training/hello-nextflow/work/8c/79499c11beea6e9d43605141f2817f/output.txt
    ```

??? abstract "Dateiinhalte"

    ```console title="results/hello_channels/output.txt"
    Hello Channels!
    ```

Wenn das bei dir funktioniert hat, bist du bereit, etwas über Kanäle zu lernen.

---

## 1. Variable Eingaben explizit über einen Kanal bereitstellen

Wir werden einen **Kanal** erstellen, um die variable Eingabe an den `sayHello()`-process zu übergeben, anstatt uns auf die implizite Handhabung zu verlassen, die gewisse Einschränkungen hat.

### 1.1. Einen Eingabe-Kanal erstellen

Es gibt verschiedene [**channel factories**](https://nextflow.io/docs/latest/reference/channel.html), die wir verwenden können, um einen Kanal einzurichten.
Um die Dinge vorerst einfach zu halten, werden wir die grundlegendste channel factory namens [`channel.of`](https://nextflow.io/docs/latest/reference/channel.html#of) verwenden, die einen Kanal mit einem einzelnen Wert erstellt.
Funktional wird dies ähnlich sein wie zuvor, aber anstatt Nextflow einen Kanal implizit erstellen zu lassen, tun wir dies jetzt explizit.

Dies ist die Codezeile, die wir verwenden werden:

```console title="Syntax"
greeting_ch = channel.of('Hello Channels!')
```

Dies erstellt einen Kanal namens `greeting_ch` unter Verwendung der `channel.of()` channel factory, die einen einfachen queue channel einrichtet, und lädt den String `'Hello Channels!'` als Begrüßungswert.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-pipeline-channel.svg"
</figure>

!!! note "Hinweis"

    Wir wechseln vorübergehend zurück zu fest codierten Strings, anstatt einen CLI-Parameter zu verwenden, um die Lesbarkeit zu verbessern. Wir werden zurück zu CLI-Parametern wechseln, sobald wir behandelt haben, was auf Ebene des Kanals passiert.

Füge im workflow-Block den channel factory-Code hinzu:

=== "Danach"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="4 5"
    workflow {

        main:
        // Einen Kanal für Eingaben erstellen
        greeting_ch = channel.of('Hello Channels!')
        // Eine Begrüßung ausgeben
        sayHello(params.input)

        publish:
        first_output = sayHello.out
    }
    ```

=== "Vorher"

    ```groovy title="hello-channels.nf" linenums="27"
    workflow {

        main:
        // Eine Begrüßung ausgeben
        sayHello(params.input)

        publish:
        first_output = sayHello.out
    }
    ```

Dies ist noch nicht funktionsfähig, da wir die Eingabe für den process-Aufruf noch nicht umgestellt haben.

### 1.2. Den Kanal als Eingabe zum process-Aufruf hinzufügen

Jetzt müssen wir unseren neu erstellten Kanal tatsächlich in den `sayHello()`-process-Aufruf einbinden und den CLI-Parameter ersetzen, den wir zuvor direkt bereitgestellt haben.

Nimm im workflow-Block die folgende Codeänderung vor:

=== "Danach"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="7"
    workflow {

        main:
        // Einen Kanal für Eingaben erstellen
        greeting_ch = channel.of('Hello Channels!')
        // Eine Begrüßung ausgeben
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

=== "Vorher"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="7"
    workflow {

        main:
        // Einen Kanal für Eingaben erstellen
        greeting_ch = channel.of('Hello Channels!')
        // Eine Begrüßung ausgeben
        sayHello(params.input)

        publish:
        first_output = sayHello.out
    }
    ```

Dies sagt Nextflow, den `sayHello`-process auf den Inhalt des `greeting_ch`-Kanals auszuführen.

Jetzt ist unser Workflow richtig funktionsfähig; es ist das explizite Äquivalent zum Schreiben von `sayHello('Hello Channels!')`.

### 1.3. Den Workflow ausführen

Lass ihn laufen!

```bash
nextflow run hello-channels.nf
```

??? success "Befehlsausgabe"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [fabulous_crick] DSL2 - revision: 23e20f76e8

    executor >  local (1)
    [c0/4f1872] process > sayHello (1) [100%] 1 of 1 ✔
    ```

Wenn du beide Änderungen korrekt vorgenommen hast, solltest du eine erfolgreiche Ausführung erhalten.
Du kannst das Ergebnisverzeichnis überprüfen, um dich zu vergewissern, dass das Ergebnis noch dasselbe ist wie zuvor.

??? abstract "Dateiinhalte"

    ```console title="results/hello_channels/output.txt"
    Hello Channels!
    ```

Wir haben also die Flexibilität unseres Workflows erhöht, während wir dasselbe Endergebnis erzielen.
Das mag so erscheinen, als würden wir mehr Code für keinen greifbaren Nutzen schreiben, aber der Wert wird deutlich, sobald wir anfangen, mehr Eingaben zu handhaben.

Als Vorschau darauf schauen wir uns noch eine Sache an, bevor wir weitermachen: einen kleinen, aber praktischen Vorteil der Verwendung eines expliziten Kanals zur Verwaltung von Dateneingaben.

### 1.4. `view()` verwenden, um Kanal-Inhalte zu inspizieren

Nextflow-Kanäle sind so aufgebaut, dass wir mit Operatoren auf ihre Inhalte einwirken können, was wir später in diesem Kapitel ausführlich behandeln werden.

Vorerst zeigen wir dir nur, wie du einen super einfachen Operator namens [`view()`](https://www.nextflow.io/docs/latest/reference/operator.html#view) verwendest, um die Inhalte eines Kanals zu inspizieren.
Du kannst `view()` als Debugging-Werkzeug betrachten, wie eine `print()`-Anweisung in Python oder deren Äquivalent in anderen Sprachen.

Füge diese kleine Zeile zum workflow-Block hinzu:

=== "Danach"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="7"
    workflow {

        main:
        // Einen Kanal für Eingaben erstellen
        greeting_ch = channel.of('Hello Channels!')
                             .view()
        // Eine Begrüßung ausgeben
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

=== "Vorher"

    ```groovy title="hello-channels.nf" linenums="27"
    workflow {

        main:
        // Einen Kanal für Eingaben erstellen
        greeting_ch = channel.of('Hello Channels!')
        // Eine Begrüßung ausgeben
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

Die genaue Anzahl der Leerzeichen spielt keine Rolle, solange es ein Vielfaches von 4 ist; wir zielen nur darauf ab, den Anfang der `.view()`-Anweisung mit dem `.of()`-Teil der Kanal-Konstruktion auszurichten.

Führe jetzt den Workflow erneut aus:

```bash
nextflow run hello-channels.nf
```

??? success "Befehlsausgabe"

    ```console hl_lines="7"
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [scruffy_shaw] DSL2 - revision: 2ede41e14a

    executor >  local (1)
    [ef/f7e40a] sayHello (1) [100%] 1 of 1 ✔
    Hello Channels!
    ```

Wie du sehen kannst, gibt dies die Kanal-Inhalte in der Konsole aus.
Hier haben wir nur ein Element, aber wenn wir im nächsten Abschnitt anfangen, mehrere Werte in den Kanal zu laden, wirst du sehen, dass dies so eingestellt ist, dass ein Element pro Zeile ausgegeben wird.

### Fazit

Du weißt, wie du eine grundlegende channel factory verwendest, um eine Eingabe an einen process zu übergeben.

### Wie geht es weiter?

Lerne, wie du Kanäle verwendest, damit der Workflow über mehrere Eingabewerte iteriert.

---

## 2. Den Workflow so modifizieren, dass er mit mehreren Eingabewerten läuft

Workflows werden typischerweise auf Stapel von Eingaben ausgeführt, die in großen Mengen verarbeitet werden sollen, daher wollen wir den Workflow aufrüsten, um mehrere Eingabewerte zu akzeptieren.

### 2.1. Mehrere Begrüßungen in den Eingabe-Kanal laden

Praktischerweise ist die `channel.of()` channel factory, die wir verwendet haben, durchaus bereit, mehr als einen Wert zu akzeptieren, daher müssen wir das überhaupt nicht ändern.
Wir können einfach mehrere Werte in den Kanal laden.

Machen wir sie zu `'Hello'`, `'Bonjour'` und `'Holà'`.

#### 2.1.1. Weitere Begrüßungen hinzufügen

Nimm vor dem workflow-Block die folgende Codeänderung vor:

=== "Danach"

    ```groovy title="hello-channels.nf" linenums="30" hl_lines="2"
    // Einen Kanal für Eingaben erstellen
    greeting_ch = channel.of('Hello','Bonjour','Holà')
                         .view()
    ```

=== "Vorher"

    ```groovy title="hello-channels.nf" linenums="30" hl_lines="2"
    // Einen Kanal für Eingaben erstellen
    greeting_ch = channel.of('Hello Channels')
                         .view()
    ```

Die Dokumentation sagt uns, dass das funktionieren sollte. Kann es wirklich so einfach sein?

#### 2.1.2. Den Befehl ausführen und die Log-Ausgabe betrachten

Lass uns es versuchen.

```bash
nextflow run hello-channels.nf
```

??? success "Befehlsausgabe"

    ```console hl_lines="6"
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [amazing_crick] DSL2 - revision: 59a9a5888a

    executor >  local (3)
    [f4/c9962c] process > sayHello (1) [100%] 3 of 3 ✔
    Hello
    Bonjour
    Holà
    ```

Es scheint tatsächlich gut gelaufen zu sein.
Der Ausführungsmonitor zeigt, dass `3 von 3` Aufrufe für den `sayHello`-process gemacht wurden, und wir sehen die drei Begrüßungen, die durch die `view()`-Anweisung aufgelistet werden, eine pro Zeile wie versprochen.

Allerdings gibt es im Ergebnisverzeichnis immer noch nur eine Ausgabe:

??? abstract "Verzeichnisinhalte"

    ```console title="results/hello_channels" hl_lines="3"
    results
    ├── hello_channels
    │   └── output.txt
    ├── hello_world
    │   └── output.txt
    └── output.txt -> /workspaces/training/hello-nextflow/work/8c/79499c11beea6e9d43605141f2817f/output.txt
    ```

??? abstract "Dateiinhalte"

    ```console title="results/hello_channels/output.txt"
    Holà
    ```

Du solltest eine der drei Begrüßungen darin sehen, aber die, die du bekommen hast, könnte anders sein als hier gezeigt.
Kannst du dir vorstellen, warum das so sein könnte?

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-pipeline-channel-multi.svg"
</figure>

_Im Diagramm wird der Kanal in Grün dargestellt, und die Reihenfolge der Elemente wird wie Murmeln in einer Röhre dargestellt: das zuerst geladene ist rechts, dann das zweite in der Mitte, dann das dritte links._

Wenn wir zurück auf den Ausführungsmonitor schauen, hat er uns nur einen Unterverzeichnispfad gegeben (`f4/c9962c`).
Lass uns dort hineinschauen.

??? abstract "Verzeichnisinhalte"

    ```console hl_lines="9"
    work/f4/c9962ce91ef87480babcb86b2b9042/
    ├── .command.begin
    ├── .command.err
    ├── .command.log
    ├── .command.out
    ├── .command.run
    ├── .command.sh
    ├── .exitcode
    └── output.txt
    ```

??? abstract "Dateiinhalte"

    ```console title="work/f4/c9962ce91ef87480babcb86b2b9042/output.txt"
    Hello
    ```

Das ist nicht einmal dieselbe Begrüßung, die wir im Ergebnisverzeichnis bekommen haben! Was geht hier vor?

An diesem Punkt müssen wir dir sagen, dass das ANSI-Logging-System standardmäßig das Logging von mehreren Aufrufen desselben process auf derselben Zeile schreibt.
Also landen die Statusmeldungen von allen drei Aufrufen des sayHello()-process an derselben Stelle.

Glücklicherweise können wir dieses Verhalten deaktivieren, um die vollständige Liste der process-Aufrufe zu sehen.

#### 2.1.3. Den Befehl erneut mit der Option `-ansi-log false` ausführen

Um das Logging zu erweitern und eine Zeile pro process-Aufruf anzuzeigen, füge `-ansi-log false` zum Befehl hinzu.

```bash
nextflow run hello-channels.nf -ansi-log false
```

??? success "Befehlsausgabe"

    ```console
     N E X T F L O W  ~  version 25.10.2
    Launching `hello-channels.nf` [desperate_monod] DSL2 - revision: 59a9a5888a
    Hello
    Bonjour
    Holà
    [23/871c7e] Submitted process > sayHello (2)
    [7f/21e2c2] Submitted process > sayHello (1)
    [f4/ea10a6] Submitted process > sayHello (3)
    ```

Diesmal sehen wir alle drei process-Ausführungen und ihre zugehörigen work-Unterverzeichnisse in der Ausgabe aufgelistet.

Das ist viel besser, zumindest für einen einfachen Workflow.
Für einen komplexen Workflow oder eine große Anzahl von Eingaben würde die vollständige Liste, die im Terminal ausgegeben wird, etwas überwältigend werden.
Deshalb ist `-ansi-log false` nicht das Standardverhalten.

!!! tip "Tipp"

    Die Art und Weise, wie der Status gemeldet wird, ist zwischen den beiden Logging-Modi etwas unterschiedlich.
    Im komprimierten Modus meldet Nextflow, ob Aufrufe erfolgreich abgeschlossen wurden oder nicht.
    In diesem erweiterten Modus wird nur gemeldet, dass sie eingereicht wurden.

Wie auch immer, jetzt da wir die Unterverzeichnisse jedes process-Aufrufs haben, können wir nach ihren Logs und Ausgaben suchen.

??? abstract "Verzeichnisinhalte"

    ```console
    work/23/871c7ec3642a898ecd5e6090d21300/
    ├── .command.begin
    ├── .command.err
    ├── .command.log
    ├── .command.out
    ├── .command.run
    ├── .command.sh
    ├── .exitcode
    └── output.txt
    ```

    ```console
    work/7f/21e2c2f3cc8833ef3858b236e5575c/
    ├── .command.begin
    ├── .command.err
    ├── .command.log
    ├── .command.out
    ├── .command.run
    ├── .command.sh
    ├── .exitcode
    └── output.txt
    ```

    ```console
    work/f4/ea10a680d5687596d3eaa3fcf69272/
    ├── .command.begin
    ├── .command.err
    ├── .command.log
    ├── .command.out
    ├── .command.run
    ├── .command.sh
    ├── .exitcode
    └── output.txt
    ```

??? abstract "Dateiinhalte"

    ```txt title="work/23/871c7ec3642a898ecd5e6090d21300/output.txt"
    Bonjour
    ```

    ```txt title="work/7f/21e2c2f3cc8833ef3858b236e5575c/output.txt"
    Hello
    ```

    ```txt title="work/f4/ea10a680d5687596d3eaa3fcf69272/output.txt"
    Holà
    ```

Dies zeigt, dass alle drei processes erfolgreich ausgeführt wurden (juhu).

Allerdings haben wir immer noch das Problem, dass es nur eine Ausgabedatei im Ergebnisverzeichnis gibt.

Du erinnerst dich vielleicht, dass wir den Ausgabedateinamen für den `sayHello`-process fest codiert haben, sodass alle drei Aufrufe eine Datei namens `output.txt` produziert haben.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-channels-task-dirs.svg"
</figure>

Solange die Ausgabedateien in den work-Unterverzeichnissen bleiben, isoliert von den anderen processes, ist das in Ordnung.
Aber wenn sie in dasselbe Ergebnisverzeichnis veröffentlicht werden, wird die, die zuerst dort kopiert wurde, von der nächsten überschrieben, und so weiter.

### 2.2. Sicherstellen, dass die Ausgabedateinamen eindeutig sind

Wir können weiterhin alle Ausgaben in dasselbe Ergebnisverzeichnis veröffentlichen, aber wir müssen sicherstellen, dass sie eindeutige Namen haben werden.
Genauer gesagt müssen wir den ersten process so modifizieren, dass er einen Dateinamen dynamisch generiert, sodass die endgültigen Dateinamen eindeutig sein werden.

Wie machen wir also die Dateinamen eindeutig?
Ein gängiger Weg, das zu tun, ist, ein eindeutiges Stück Metadaten aus den Eingaben (die vom Eingabe-Kanal empfangen wurden) als Teil des Ausgabedateinamens zu verwenden.
Hier werden wir der Einfachheit halber nur die Begrüßung selbst verwenden, da sie nur ein kurzer String ist, und sie dem Basis-Ausgabedateinamen voranstellen.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-pipeline-channel-multi-unique.svg"
</figure>

#### 2.2.1. Einen dynamischen Ausgabedateinamen konstruieren

Nimm im process-Block die folgenden Codeänderungen vor:

=== "Danach"

    ```groovy title="hello-channels.nf" linenums="6" hl_lines="7 11"
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

=== "Vorher"

    ```groovy title="hello-channels.nf" linenums="6" hl_lines="7 11"
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

Stelle sicher, dass du `output.txt` sowohl in der Ausgabedefinition als auch im `script:`-Befehlsblock ersetzt.

!!! tip "Tipp"

    In der Ausgabedefinition MUSST du doppelte Anführungszeichen um den Ausgabedateinamensausdruck verwenden (NICHT einfache Anführungszeichen), sonst schlägt es fehl.

Dies sollte jedes Mal einen eindeutigen Ausgabedateinamen produzieren, wenn der process aufgerufen wird, sodass er von den Ausgaben anderer Aufrufe desselben process im Ausgabeverzeichnis unterschieden werden kann.

#### 2.2.2. Den Workflow ausführen

Lass ihn laufen. Beachte, dass wir wieder mit den Standard-ANSI-Log-Einstellungen ausführen.

```bash
nextflow run hello-channels.nf
```

??? success "Befehlsausgabe"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [sharp_minsky] DSL2 - revision: 16a291febe

    executor >  local (3)
    [e8/33ee64] sayHello (2) [100%] 3 of 3 ✔
    Hello
    Bonjour
    Holà
    ```

Zurück in der Zusammenfassungsansicht wird die Ausgabe wieder auf einer Zeile zusammengefasst.
Schau in das `results`-Verzeichnis, um zu sehen, ob alle Ausgabebegrüßungen da sind.

??? abstract "Verzeichnisinhalte"

    ```console
    results/hello_channels/
    ├── Bonjour-output.txt
    ├── Hello-output.txt
    ├── Holà-output.txt
    └── output.txt
    ```

Ja! Und sie haben jeweils die erwarteten Inhalte.

??? abstract "Dateiinhalte"

    ```console title="Bonjour-output.txt"
    Bonjour
    ```

    ```console title="Hello-output.txt"
    Hello
    ```

    ```console title="Holà-output.txt"
    Holà
    ```

Erfolg! Jetzt können wir so viele Begrüßungen hinzufügen, wie wir möchten, ohne uns Sorgen zu machen, dass Ausgabedateien überschrieben werden.

!!! tip "Tipp"

    In der Praxis ist das Benennen von Dateien basierend auf den Eingabedaten selbst fast immer unpraktisch.
    Der bessere Weg, dynamische Dateinamen zu generieren, ist, Metadaten zusammen mit den Eingabedateien an einen process zu übergeben.
    Die Metadaten werden typischerweise über ein 'Sample Sheet' oder Äquivalente bereitgestellt.
    Du wirst später in deinem Nextflow-Training lernen, wie du das machst (siehe [Metadaten-Side Quest](../side_quests/metadata.md)).

### Fazit

Du weißt, wie du mehrere Eingabeelemente durch einen Kanal führst.

### Wie geht es weiter?

Lerne, einen Operator zu verwenden, um die Inhalte eines Kanals zu transformieren.

---

## 3. Mehrere Eingaben über ein Array bereitstellen

Wir haben dir gerade gezeigt, wie du mehrere Eingabeelemente handhabst, die direkt in der channel factory fest codiert waren.
Was ist, wenn wir diese mehreren Eingaben auf eine andere Weise bereitstellen wollten?

Stell dir zum Beispiel vor, wir richten eine Eingabevariable ein, die ein Array von Elementen enthält wie dieses:

`greetings_array = ['Hello','Bonjour','Holà']`

Können wir das in unseren Ausgabe-Kanal laden und erwarten, dass es funktioniert?

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-pipeline-multi-inputs-array.svg"
</figure>

Lass uns das herausfinden.

### 3.1. Ein Array von Werten als Eingabe für den Kanal bereitstellen

Der gesunde Menschenverstand legt nahe, dass wir einfach ein Array von Werten anstelle eines einzelnen Wertes übergeben können sollten.
Lass uns es versuchen; wir müssen die Eingabevariable einrichten und sie in die channel factory laden.

#### 3.1.1. Die Eingabevariable einrichten

Lass uns die `greetings_array`-Variable, die wir uns gerade vorgestellt haben, Realität werden lassen, indem wir sie zum workflow-Block hinzufügen:

=== "Danach"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="4 5"
    workflow {

        main:
        // Ein Array von Eingabe-Begrüßungen deklarieren
        greetings_array = ['Hello','Bonjour','Holà']
        // Einen Kanal für Eingaben erstellen
        greeting_ch = channel.of('Hello','Bonjour','Holà')
                             .view()
        // Eine Begrüßung ausgeben
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

=== "Vorher"

    ```groovy title="hello-channels.nf" linenums="27"
    workflow {

        main:
        // Einen Kanal für Eingaben erstellen
        greeting_ch = channel.of('Hello','Bonjour','Holà')
                             .view()
        // Eine Begrüßung ausgeben
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

Dies ist noch nicht funktionsfähig, wir haben nur eine Deklaration für das Array hinzugefügt.

#### 3.1.2. Das Array von Begrüßungen als Eingabe für die channel factory setzen

Jetzt werden wir die Werte `'Hello','Bonjour','Holà'`, die derzeit in der channel factory fest codiert sind, durch das `greetings_array` ersetzen, das wir gerade erstellt haben.

Nimm im workflow-Block die folgende Änderung vor:

=== "Danach"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="7"
    workflow {

        main:
        // Ein Array von Eingabe-Begrüßungen deklarieren
        greetings_array = ['Hello','Bonjour','Holà']
        // Einen Kanal für Eingaben erstellen
        greeting_ch = channel.of(greetings_array)
                             .view()
        // Eine Begrüßung ausgeben
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

=== "Vorher"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="7"
    workflow {

        main:
        // Ein Array von Eingabe-Begrüßungen deklarieren
        greetings_array = ['Hello','Bonjour','Holà']
        // Einen Kanal für Eingaben erstellen
        greeting_ch = channel.of('Hello','Bonjour','Holà')
                             .view()
        // Eine Begrüßung ausgeben
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

Dies sollte jetzt funktionsfähig sein.

#### 3.1.3. Den Workflow ausführen

Lass uns versuchen, ihn auszuführen:

```bash
nextflow run hello-channels.nf
```

??? failure "Befehlsausgabe"

    ```console hl_lines="7 11 16"
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [friendly_koch] DSL2 - revision: 97256837a7

    executor >  local (1)
    [a8/1f6ead] sayHello (1) | 0 of 1
    [Hello, Bonjour, Holà]
    ERROR ~ Error executing process > 'sayHello (1)'

    Caused by:
      Missing output file(s) `[Hello, Bonjour, Holà]-output.txt` expected by process `sayHello (1)`


    Command executed:

      echo '[Hello, Bonjour, Holà]' > '[Hello, Bonjour, Holà]-output.txt'

    Command exit status:
      0

    Command output:
      (empty)

    Work dir:
      /workspaces/training/hello-nextflow/work/a8/1f6ead5f3fa30a3c508e2e7cf83ffb

    Tip: you can replicate the issue by changing to the process work dir and entering the command `bash .command.run`

    -- Check '.nextflow.log' file for details
    ```

Oh nein! Es gibt einen Fehler!

Schau dir die Ausgabe von `view()` und die Fehlermeldungen an.

Es sieht so aus, als hätte Nextflow versucht, einen einzelnen process-Aufruf auszuführen, wobei `[Hello, Bonjour, Holà]` als Stringwert verwendet wurde, anstatt die drei Strings im Array als separate Werte zu verwenden.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-channels-array-fail.svg"
</figure>

Es ist also die 'Verpackung', die das Problem verursacht.
Wie bringen wir Nextflow dazu, das Array zu entpacken und die einzelnen Strings in den Kanal zu laden?

### 3.2. Einen Operator verwenden, um Kanal-Inhalte zu transformieren

Hier kommen [**Operatoren**](https://www.nextflow.io/docs/latest/reference/operator.html) ins Spiel.
Du hast bereits den `.view()`-Operator verwendet, der nur hineinschaut, was drin ist.
Jetzt werden wir uns Operatoren ansehen, die es uns ermöglichen, auf die Inhalte eines Kanals einzuwirken.

Wenn du durch die [Liste der Operatoren](https://www.nextflow.io/docs/latest/reference/operator.html) in der Nextflow-Dokumentation blätterst, wirst du [`flatten()`](https://www.nextflow.io/docs/latest/reference/operator.html#flatten) finden, das genau das tut, was wir brauchen: den Inhalt eines Arrays entpacken und sie als einzelne Elemente ausgeben.

#### 3.2.1. Den `flatten()`-Operator hinzufügen

Um den `flatten()`-Operator auf unseren Eingabe-Kanal anzuwenden, hängen wir ihn an die channel factory-Deklaration an.

Nimm im workflow-Block die folgende Codeänderung vor:

=== "Danach"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="9"
    workflow {

        main:
        // Ein Array von Eingabe-Begrüßungen deklarieren
        greetings_array = ['Hello','Bonjour','Holà']
        // Einen Kanal für Eingaben erstellen
        greeting_ch = channel.of(greetings_array)
                             .view()
                             .flatten()
        // Eine Begrüßung ausgeben
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

=== "Vorher"

    ```groovy title="hello-channels.nf" linenums="27"
    workflow {

        main:
        // Ein Array von Eingabe-Begrüßungen deklarieren
        greetings_array = ['Hello','Bonjour','Holà']
        // Einen Kanal für Eingaben erstellen
        greeting_ch = channel.of(greetings_array)
                             .view()
        // Eine Begrüßung ausgeben
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

Hier haben wir den Operator für die Lesbarkeit in der nächsten Zeile hinzugefügt, aber du kannst Operatoren auch auf derselben Zeile wie die channel factory hinzufügen, wenn du möchtest, so:
`greeting_ch = channel.of(greetings_array).view().flatten()`

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-channels-array-success.svg"
</figure>

#### 3.2.2. Die `view()`-Anweisung(en) verfeinern

Wir könnten dies jetzt sofort ausführen, um zu testen, ob es funktioniert, aber während wir dabei sind, werden wir verfeinern, wie wir die Kanal-Inhalte inspizieren.

Wir wollen den Inhalt vor und nach der Anwendung des `flatten()`-Operators vergleichen können, also werden wir ein zweites hinzufügen, UND wir werden etwas Code hinzufügen, um sie in der Ausgabe deutlicher zu beschriften.

Nimm im workflow-Block die folgende Codeänderung vor:

=== "Danach"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="8-10"
    workflow {

        main:
        // Ein Array von Eingabe-Begrüßungen deklarieren
        greetings_array = ['Hello','Bonjour','Holà']
        // Einen Kanal für Eingaben erstellen
        greeting_ch = channel.of(greetings_array)
                             .view { greeting -> "Vor flatten: $greeting" }
                             .flatten()
                             .view { greeting -> "Nach flatten: $greeting" }
        // Eine Begrüßung ausgeben
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

=== "Vorher"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="8-9"
    workflow {

        main:
        // Ein Array von Eingabe-Begrüßungen deklarieren
        greetings_array = ['Hello','Bonjour','Holà']
        // Einen Kanal für Eingaben erstellen
        greeting_ch = channel.of(greetings_array)
                             .view()
                             .flatten()
        // Eine Begrüßung ausgeben
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

Du siehst, wir haben eine zweite `.view`-Anweisung hinzugefügt, und für jede von ihnen haben wir die leeren Klammern (`()`) durch geschweifte Klammern mit etwas Code ersetzt, wie `{ greeting -> "Vor flatten: $greeting" }`.

Diese werden _closures_ genannt. Der Code, den sie enthalten, wird für jedes Element im Kanal ausgeführt.
Wir definieren eine temporäre Variable für den inneren Wert, hier `greeting` genannt (aber es könnte ein beliebiger Name sein), die nur innerhalb des Gültigkeitsbereichs dieser closure verwendet wird.

In diesem Beispiel repräsentiert `$greeting` jedes einzelne Element, das in den Kanal geladen wurde.
Dies führt zu sauber beschrifteter Konsolenausgabe.

!!! info

    In einigen Pipelines siehst du möglicherweise eine spezielle Variable namens `$it`, die innerhalb von Operator-closures verwendet wird.
    Dies ist eine _implizite_ Variable, die einen Kurzform-Zugriff auf die innere Variable ermöglicht,
    ohne sie mit einem `->` definieren zu müssen.

    Wir bevorzugen die explizite Form, um die Klarheit des Codes zu unterstützen. Die `$it`-Syntax wird nicht empfohlen und wird langsam aus der Nextflow-Sprache auslaufen.

#### 3.2.3. Den Workflow ausführen

Schließlich kannst du versuchen, den Workflow erneut auszuführen!

```bash
nextflow run hello-channels.nf
```

??? success "Befehlsausgabe"

    ```console hl_lines="7-10"
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [sleepy_gutenberg] DSL2 - revision: 1db4f760ee

    executor >  local (3)
    [b1/6a1e15] sayHello (2) [100%] 3 of 3 ✔
    Vor flatten: [Hello, Bonjour, Holà]
    Nach flatten: Hello
    Nach flatten: Bonjour
    Nach flatten: Holà
    ```

Diesmal funktioniert es UND gibt uns den zusätzlichen Einblick, wie die Inhalte des Kanals vor und nach der Ausführung des `flatten()`-Operators aussehen.

- Eine einzelne `Vor flatten:`-Anweisung, weil zu diesem Zeitpunkt der Kanal ein Element enthält, das ursprüngliche Array.
- Drei separate `Nach flatten:`-Anweisungen, eine für jede Begrüßung, die jetzt einzelne Elemente im Kanal sind.

Wichtig ist, dass dies bedeutet, dass jedes Element jetzt separat vom Workflow verarbeitet werden kann.

!!! tip "Tipp"

    Es ist technisch möglich, dieselben Ergebnisse zu erzielen, indem man eine andere channel factory verwendet, [`channel.fromList`](https://nextflow.io/docs/latest/reference/channel.html#fromlist), die einen impliziten Mapping-Schritt in ihrer Operation enthält.
    Hier haben wir uns entschieden, das nicht zu verwenden, um die Verwendung eines Operators an einem einfachen Anwendungsfall zu demonstrieren.

### Fazit

Du weißt, wie du einen Operator wie `flatten()` verwendest, um die Inhalte eines Kanals zu transformieren, und wie du den `view()`-Operator verwendest, um Kanal-Inhalte vor und nach der Anwendung eines Operators zu inspizieren.

### Wie geht es weiter?

Lerne, wie du den Workflow so einrichtest, dass er eine Datei als Quelle von Eingabewerten verwendet.

---

## 4. Eingabewerte aus einer CSV-Datei lesen

Realistischerweise werden wir selten, wenn überhaupt, von einem Array von Werten ausgehen.
Höchstwahrscheinlich haben wir eine oder mehrere Dateien, die die zu verarbeitenden Daten enthalten, in irgendeiner Art von strukturiertem Format.

Wir haben eine CSV-Datei namens `greetings.csv` vorbereitet, die mehrere Eingabebegrüßungen enthält und die Art von tabellarischen Daten nachahmt, die du in einer echten Datenanalyse verarbeiten möchtest, gespeichert unter `data/`.
(Die Zahlen sind nicht bedeutsam, sie sind nur zu Illustrationszwecken da.)

```csv title="data/greetings.csv" linenums="1"
Hello,English,123
Bonjour,French,456
Holà,Spanish,789
```

Unsere nächste Aufgabe ist es, unseren Workflow anzupassen, um die Werte aus dieser Datei einzulesen.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-pipeline-multi-inputs-csv.svg"
</figure>

Lass uns sehen, wie wir das bewerkstelligen können.

### 4.1. Das Skript so modifizieren, dass es eine CSV-Datei als Quelle von Begrüßungen erwartet

Um loszulegen, müssen wir zwei wichtige Änderungen am Skript vornehmen:

- Den Eingabeparameter so umstellen, dass er auf die CSV-Datei zeigt
- Die channel factory auf eine umstellen, die für das Handhaben einer Datei ausgelegt ist

#### 4.1.1. Den Eingabeparameter so umstellen, dass er auf die CSV-Datei zeigt

Erinnerst du dich an den `params.input`-Parameter, den wir in Teil 1 eingerichtet haben?
Wir werden ihn aktualisieren, um auf die CSV-Datei zu zeigen, die unsere Begrüßungen enthält.

Nimm die folgende Änderung an der Parameterdeklaration vor:

=== "Danach"

    ```groovy title="hello-channels.nf" linenums="20" hl_lines="5"
    /*
    * Pipeline-Parameter
    */
    params {
        input: Path = 'data/greetings.csv'
    }
    ```

=== "Vorher"

    ```groovy title="hello-channels.nf" linenums="20" hl_lines="5"
    /*
     * Pipeline-Parameter
     */
    input: String = 'Holà mundo!'
    ```

Dies setzt voraus, dass die Datei am selben Ort wie der Workflow-Code liegt.
Wie du mit anderen Datenspeicherorten umgehst, lernst du in späteren Kursen.

#### 4.1.2. Zu einer channel factory wechseln, die für das Handhaben einer Datei ausgelegt ist

Da wir jetzt eine Datei anstelle von einfachen Strings als Eingabe verwenden wollen, können wir die `channel.of()` channel factory von vorher nicht verwenden.
Wir müssen zu einer neuen channel factory wechseln, [`channel.fromPath()`](https://nextflow.io/docs/latest/reference/channel.html#frompath), die einige eingebaute Funktionalitäten für das Handhaben von Dateipfaden hat.

Nimm im workflow-Block die folgende Codeänderung vor:

=== "Danach"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="4-8"
    workflow {

        main:
        // Einen Kanal für Eingaben aus einer CSV-Datei erstellen
        greeting_ch = channel.fromPath(params.input)
                             .view { greeting -> "Vor flatten: $greeting" }
                             // .flatten()
                             // .view { greeting -> "Nach flatten: $greeting" }
        // Eine Begrüßung ausgeben
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

=== "Vorher"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="4-8"
    workflow {

        main:
        // Ein Array von Eingabe-Begrüßungen deklarieren
        greetings_array = ['Hello','Bonjour','Holà']
        // Einen Kanal für Eingaben erstellen
        greeting_ch = channel.of(greetings_array)
                             .view { greeting -> "Vor flatten: $greeting" }
                             .flatten()
                             .view { greeting -> "Nach flatten: $greeting" }
        // Eine Begrüßung ausgeben
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

Du wirst bemerken, dass wir die Kanal-Eingabe zurück zu `param.input` gewechselt und die `greetings_array`-Deklaration gelöscht haben, da wir sie nicht mehr brauchen werden.
Wir haben auch `flatten()` und die zweite `view()`-Anweisung auskommentiert.

#### 4.1.3. Den Workflow ausführen

Lass uns versuchen, den Workflow mit der neuen channel factory und der Eingabedatei auszuführen.

```bash
nextflow run hello-channels.nf
```

??? failure "Befehlsausgabe"

    ```console hl_lines="5 6 9 14"
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [peaceful_poisson] DSL2 - revision: a286c08ad5

    [-        ] sayHello [  0%] 0 of 1
    Vor flatten: /workspaces/training/hello-nextflow/data/greetings.csv
    ERROR ~ Error executing process > 'sayHello (1)'

    Caused by:
      File `/workspaces/training/hello-nextflow/data/greetings.csv-output.txt` is outside the scope of the process work directory: /workspaces/training/hello-nextflow/work/30/e610cb4ea5ae8693f456ac3329c92f


    Command executed:

      echo '/workspaces/training/hello-nextflow/data/greetings.csv' > '/workspaces/training/hello-nextflow/data/greetings.csv-output.txt'

    Command exit status:
      -

    Command output:
      (empty)

    Work dir:
      /workspaces/training/hello-nextflow/work/30/e610cb4ea5ae8693f456ac3329c92f

    Tip: when you have fixed the problem you can continue the execution adding the option `-resume` to the run command line

    -- Check '.nextflow.log' file for details
    ```

Oh nein, es funktioniert nicht. Schau dir den Anfang der Konsolenausgabe und die Fehlermeldung an.
Der Teil `Command executed:` ist hier besonders hilfreich.

Das sieht vielleicht ein bisschen vertraut aus.
Es sieht so aus, als hätte Nextflow versucht, einen einzelnen process-Aufruf auszuführen, wobei der Dateipfad selbst als Stringwert verwendet wurde.
Also hat es den Dateipfad korrekt aufgelöst, aber es hat nicht wirklich seinen Inhalt geparst, was wir wollten.

Wie bringen wir Nextflow dazu, die Datei zu öffnen und ihren Inhalt in den Kanal zu laden?

Klingt so, als bräuchten wir einen weiteren [Operator](https://www.nextflow.io/docs/latest/reference/operator.html)!

### 4.2. Den `splitCsv()`-Operator verwenden, um die Datei zu parsen

Wenn wir wieder durch die Liste der Operatoren schauen, finden wir [`splitCsv()`](https://www.nextflow.io/docs/latest/reference/operator.html#splitcsv), das dafür ausgelegt ist, CSV-formatierten Text zu parsen und zu teilen.

#### 4.2.1. `splitCsv()` auf den Kanal anwenden

Um den Operator anzuwenden, hängen wir ihn wie zuvor an die channel factory-Zeile an.

Nimm im workflow-Block die folgende Codeänderung vor, um `flatten()` durch `splitcsv()` zu ersetzen (nicht auskommentiert):

=== "Danach"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="6-8"
    workflow {

        main:
        // Einen Kanal für Eingaben aus einer CSV-Datei erstellen
        greeting_ch = channel.fromPath(params.input)
                             .view { csv -> "Vor splitCsv: $csv" }
                             .splitCsv()
                             .view { csv -> "Nach splitCsv: $csv" }
        // Eine Begrüßung ausgeben
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

=== "Vorher"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="6-8"
    workflow {

        main:
        // Einen Kanal für Eingaben aus einer CSV-Datei erstellen
        greeting_ch = channel.fromPath(params.input)
                             .view { greeting -> "Vor flatten: $greeting" }
                             // .flatten()
                             // .view { greeting -> "Nach flatten: $greeting" }
        // Eine Begrüßung ausgeben
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

Wie du sehen kannst, haben wir auch die vorher/nachher `view()`-Anweisungen aktualisiert.
Technisch hätten wir denselben Variablennamen (`greeting`) verwenden können, aber wir haben ihn auf etwas Passenderes (`csv`) aktualisiert, um den Code für andere lesbarer zu machen.

#### 4.2.2. Den Workflow erneut ausführen

Lass uns versuchen, den Workflow mit der hinzugefügten CSV-Parsing-Logik auszuführen.

```bash
nextflow run hello-channels.nf
```

??? failure "Befehlsausgabe"

    ```console hl_lines="7-11 14 19"
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [insane_fermat] DSL2 - revision: 8e62fcbeb1

    executor >  local (3)
    [24/76da2f] sayHello (2) [  0%] 0 of 3 ✘
    Vor splitCsv: /workspaces/training/hello-nextflow/data/greetings.csv
    Nach splitCsv: [Hello, English, 123]
    Nach splitCsv: [Bonjour, French, 456]
    Nach splitCsv: [Holà, Spanish, 789]
    ERROR ~ Error executing process > 'sayHello (2)'

    Caused by:
      Missing output file(s) `[Bonjour, French, 456]-output.txt` expected by process `sayHello (2)`


    Command executed:

      echo '[Bonjour, French, 456]' > '[Bonjour, French, 456]-output.txt'

    Command exit status:
      0

    Command output:
      (empty)

    Work dir:
      /workspaces/training/hello-nextflow/work/24/76da2fcc4876b61632749f99e26a50

    Tip: you can try to figure out what's wrong by changing to the process work dir and showing the script file named `.command.sh`

    -- Check '.nextflow.log' file for details
    ```

Interessanterweise schlägt dies auch fehl, aber mit einem anderen Fehler.
Diesmal hat Nextflow den Inhalt der Datei geparst (juhu!), aber es hat jede Zeile als Array geladen, und jedes Array ist ein Element im Kanal.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-channels-split-fail.svg"
</figure>

Wir müssen ihm sagen, dass es nur die erste Spalte in jeder Zeile nehmen soll.
Wie entpacken wir das also?

Wir haben zuvor `flatten()` verwendet, um die Inhalte eines Kanals zu entpacken, aber das würde hier nicht funktionieren, weil flatten _alles_ entpackt (versuche es gerne, wenn du es selbst sehen möchtest).

Stattdessen werden wir einen anderen Operator namens `map()` verwenden, der wirklich nützlich ist und in Nextflow-Pipelines häufig vorkommt.

### 4.3. Den `map()`-Operator verwenden, um die Begrüßungen zu extrahieren

Der [`map()`](https://www.nextflow.io/docs/latest/reference/operator.html#map)-Operator ist ein sehr praktisches kleines Werkzeug, das es uns ermöglicht, alle Arten von Mappings auf die Inhalte eines Kanals durchzuführen.

In diesem Fall werden wir ihn verwenden, um das eine Element zu extrahieren, das wir aus jeder Zeile in unserer Datendatei wollen.
So sieht die Syntax aus:

```groovy title="Syntax"
.map { row -> row[0] }
```

Das bedeutet 'für jede Zeile im Kanal, nimm das 0. (erste) Element, das sie enthält'.

Also wenden wir das auf unser CSV-Parsing an.

#### 4.3.1. `map()` auf den Kanal anwenden

Nimm im workflow-Block die folgende Codeänderung vor:

=== "Danach"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="9 10"
    workflow {

        main:
        // Einen Kanal für Eingaben aus einer CSV-Datei erstellen
        greeting_ch = channel.fromPath(params.input)
                             .view { csv -> "Vor splitCsv: $csv" }
                             .splitCsv()
                             .view { csv -> "Nach splitCsv: $csv" }
                             .map { item -> item[0] }
                             .view { csv -> "Nach map: $csv" }
        // Eine Begrüßung ausgeben
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

=== "Vorher"

    ```groovy title="hello-channels.nf" linenums="27"
    workflow {

        main:
        // Einen Kanal für Eingaben aus einer CSV-Datei erstellen
        greeting_ch = channel.fromPath(params.input)
                             .view { csv -> "Vor splitCsv: $csv" }
                             .splitCsv()
                             .view { csv -> "Nach splitCsv: $csv" }
        // Eine Begrüßung ausgeben
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

Du siehst, wir haben einen weiteren `view()`-Aufruf hinzugefügt, um zu bestätigen, dass der Operator das tut, was wir erwarten.

#### 4.3.2. Den Workflow ausführen

Lass uns das noch einmal ausführen:

```bash
nextflow run hello-channels.nf
```

??? success "Befehlsausgabe"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [focused_volhard] DSL2 - revision: de435e45be

    executor >  local (3)
    [54/6eebe3] sayHello (3) [100%] 3 of 3 ✔
    Vor splitCsv: /workspaces/training/hello-nextflow/data/greetings.csv
    Nach splitCsv: [Hello, English, 123]
    Nach splitCsv: [Bonjour, French, 456]
    Nach splitCsv: [Holà, Spanish, 789]
    Nach map: Hello
    Nach map: Bonjour
    Nach map: Holà
    ```

Diesmal sollte es ohne Fehler laufen.

Wenn du dir die Ausgabe der `view()`-Anweisungen anschaust, siehst
