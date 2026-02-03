# Teil 3: Hello Workflow

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<!--
<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/zJP7cUYPEbA?si=Irl9nAQniDyICp2b&amp;list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } Siehe [die gesamte Playlist](https://www.youtube.com/playlist?list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik) auf dem Nextflow YouTube-Kanal.

:green_book: Das Videotranskript ist [hier](./transcripts/03_hello_workflow.md) verfügbar.
///
-->

Die meisten realen Workflows umfassen mehr als einen Schritt.
In diesem Trainingsmodul lernst du, wie du Processes in einem mehrstufigen Workflow verbindest.

Dies wird dir den Nextflow-Weg beibringen, um Folgendes zu erreichen:

1. Daten von einem Process zum nächsten fließen lassen
2. Ausgaben von mehreren Process-Aufrufen in einen einzelnen Process-Aufruf sammeln
3. Zusätzliche Parameter an einen Process übergeben
4. Mehrere Ausgaben aus einem Process handhaben

Zur Demonstration werden wir weiterhin auf dem domänenunabhängigen Hello World-Beispiel aus den Teilen 1 und 2 aufbauen.
Diesmal werden wir die folgenden Änderungen an unserem Workflow vornehmen, um besser widerzuspiegeln, wie Leute tatsächliche Workflows erstellen:

1. Einen zweiten Schritt hinzufügen, der die Begrüßung in Großbuchstaben umwandelt.
2. Einen dritten Schritt hinzufügen, der alle transformierten Begrüßungen sammelt und sie in eine einzelne Datei schreibt.
3. Einen Parameter hinzufügen, um die endgültige Ausgabedatei zu benennen, und diesen als sekundäre Eingabe an den Sammelschritt übergeben.
4. Den Sammelschritt auch eine einfache Statistik über das Verarbeitete berichten lassen.

??? info "Wie du von diesem Abschnitt aus beginnen kannst"

    Dieser Abschnitt des Kurses setzt voraus, dass du die Teile 1-2 des [Hello Nextflow](./index.md)-Kurses abgeschlossen hast, aber wenn du mit den dort behandelten Grundlagen vertraut bist, kannst du von hier aus starten, ohne etwas Besonderes tun zu müssen.

---

## 0. Aufwärmen: `hello-workflow.nf` ausführen

Wir werden das Workflow-Skript `hello-workflow.nf` als Ausgangspunkt verwenden.
Es entspricht dem Skript, das durch das Durcharbeiten von Teil 2 dieses Trainingskurses erstellt wurde, außer dass wir die `view()`-Anweisungen entfernt und das Ausgabeziel geändert haben:

```groovy title="hello-workflow.nf" linenums="37" hl_lines="3"
output {
    first_output {
        path 'hello_workflow'
        mode 'copy'
    }
}
```

Dieses Diagramm fasst die aktuelle Funktionsweise des Workflows zusammen.
Es sollte vertraut aussehen, außer dass wir jetzt explizit zeigen, dass die Ausgaben des Process in einem Kanal verpackt werden, genau wie die Eingaben.
Wir werden diesen Ausgabekanal gleich gut nutzen.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-workflow-channels.svg"
</figure>

Um sicherzustellen, dass alles funktioniert, führe das Skript einmal aus, bevor du Änderungen vornimmst:

```bash
nextflow run hello-workflow.nf
```

??? success "Befehlsausgabe"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-workflow.nf` [admiring_lamarr] DSL2 - revision: 4d4053520d

    executor >  local (3)
    [b1/5826b5] process > sayHello (2) [100%] 3 of 3 ✔
    ```

Wie zuvor findest du die Ausgabedateien an dem im `output`-Block angegebenen Ort.
Für dieses Kapitel ist es unter `results/hello_workflow/`.

??? abstract "Verzeichnisinhalte"

    ```console
    results/hello_workflow
    ├── Bonjour-output.txt
    ├── Hello-output.txt
    └── Holà-output.txt
    ```

Wenn das bei dir funktioniert hat, bist du bereit zu lernen, wie man einen mehrstufigen Workflow zusammenstellt.

---

## 1. Einen zweiten Schritt zum Workflow hinzufügen

Wir werden einen Schritt hinzufügen, um jede Begrüßung in Großbuchstaben umzuwandeln.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-multistep.svg"
</figure>

Dazu müssen wir drei Dinge tun:

- Den Befehl definieren, den wir für die Großbuchstaben-Umwandlung verwenden werden.
- Einen neuen Process schreiben, der den Großbuchstaben-Befehl umhüllt.
- Den neuen Process im Workflow-Block aufrufen und so einrichten, dass er die Ausgabe des `sayHello()`-Process als Eingabe nimmt.

### 1.1. Den Großbuchstaben-Befehl definieren und im Terminal testen

Für die Umwandlung der Begrüßungen in Großbuchstaben werden wir ein klassisches UNIX-Tool namens `tr` für 'text replacement' verwenden, mit der folgenden Syntax:

```bash title="Syntax"
tr '[a-z]' '[A-Z]'
```

Dies ist ein sehr naiver Textersetzungs-Einzeiler, der keine akzentuierten Buchstaben berücksichtigt, also wird zum Beispiel 'Holà' zu 'HOLà', aber er wird für die Demonstration der Nextflow-Konzepte gut genug funktionieren, und das ist, was zählt.

Um es auszuprobieren, können wir den Befehl `echo 'Hello World'` ausführen und seine Ausgabe an den `tr`-Befehl pipen:

```bash
echo 'Hello World' | tr '[a-z]' '[A-Z]' > UPPER-output.txt
```

Die Ausgabe ist eine Textdatei namens `UPPER-output.txt`, die die Großbuchstaben-Version des `Hello World`-Strings enthält.

??? abstract "Dateiinhalte"

    ```console title="UPPER-output.txt"
    HELLO WORLD
    ```

Das ist im Grunde das, was wir mit unserem Workflow versuchen werden zu tun.

### 1.2. Den Großbuchstaben-Schritt als Nextflow Process schreiben

Wir können unseren neuen Process nach dem ersten modellieren, da wir alle dieselben Komponenten verwenden wollen.

Füge die folgende Process-Definition zum Workflow-Skript hinzu, direkt unter dem ersten:

```groovy title="hello-workflow.nf" linenums="20"
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
    cat '$input_file' | tr '[a-z]' '[A-Z]' > 'UPPER-${input_file}'
    """
}
```

In diesem Fall setzen wir den zweiten Ausgabedateinamen basierend auf dem Eingabedateinamen zusammen, ähnlich wie wir es ursprünglich für die Ausgabe des ersten Process getan haben.

### 1.3. Einen Aufruf des neuen Process im Workflow-Block hinzufügen

Jetzt müssen wir Nextflow sagen, dass es den Process, den wir gerade definiert haben, tatsächlich aufrufen soll.

Nimm im Workflow-Block die folgende Codeänderung vor:

=== "Nachher"

    ```groovy title="hello-workflow.nf" linenums="44" hl_lines="10-11"
    workflow {

        main:
        // Einen Kanal für Eingaben aus einer CSV-Datei erstellen
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        // Eine Begrüßung ausgeben
        sayHello(greeting_ch)
        // Die Begrüßung in Großbuchstaben umwandeln
        convertToUpper()

        publish:
        first_output = sayHello.out
    }
    ```

=== "Vorher"

    ```groovy title="hello-workflow.nf" linenums="44"
    workflow {

        main:
        // Einen Kanal für Eingaben aus einer CSV-Datei erstellen
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        // Eine Begrüßung ausgeben
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

Dies ist noch nicht funktionsfähig, weil wir nicht angegeben haben, was in den `convertToUpper()`-Process eingegeben werden soll.

### 1.4. Die Ausgabe des ersten Process an den zweiten Process übergeben

Jetzt müssen wir die Ausgabe des `sayHello()`-Process in den `convertToUpper()`-Process fließen lassen.

Praktischerweise verpackt Nextflow automatisch die Ausgabe eines Process in einen Kanal, wie im Diagramm im Aufwärmabschnitt gezeigt.
Wir können auf den Ausgabekanal eines Process als `<process>.out` verweisen.

Also ist die Ausgabe des `sayHello`-Process ein Kanal namens `sayHello.out`, den wir direkt in den Aufruf von `convertToUpper()` einstecken können.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-multistep-connector.svg"
</figure>

Nimm im Workflow-Block die folgende Codeänderung vor:

=== "Nachher"

    ```groovy title="hello-workflow.nf" linenums="53" hl_lines="2"
        // Die Begrüßung in Großbuchstaben umwandeln
        convertToUpper(sayHello.out)
    ```

=== "Vorher"

    ```groovy title="hello-workflow.nf" linenums="53" hl_lines="2"
        // Die Begrüßung in Großbuchstaben umwandeln
        convertToUpper()
    ```

Für einen einfachen Fall wie diesen (eine Ausgabe zu einer Eingabe) ist das alles, was wir tun müssen, um zwei Processes zu verbinden!

### 1.5. Die Workflow-Ausgabeveröffentlichung einrichten

Zum Schluss aktualisieren wir die Workflow-Ausgaben, um auch die Ergebnisse des zweiten Process zu veröffentlichen.

#### 1.5.1. Den `publish:`-Abschnitt des `workflow`-Blocks aktualisieren

Nimm im `workflow`-Block die folgende Codeänderung vor:

=== "Nachher"

    ```groovy title="hello-workflow.nf" linenums="56" hl_lines="3"
        publish:
        first_output = sayHello.out
        uppercased = convertToUpper.out
    }
    ```

=== "Vorher"

    ```groovy title="hello-workflow.nf" linenums="56"
        publish:
        first_output = sayHello.out
    }
    ```

Die Logik ist dieselbe wie zuvor.

#### 1.5.2. Den `output`-Block aktualisieren

Nimm im `output`-Block die folgende Codeänderung vor:

=== "Nachher"

    ```groovy title="hello-workflow.nf" linenums="61" hl_lines="6-9"
    output {
        first_output {
            path 'hello_workflow'
            mode 'copy'
        }
        uppercased {
            path 'hello_workflow'
            mode 'copy'
        }
    }
    ```

=== "Vorher"

    ```groovy title="hello-workflow.nf" linenums="61"
    output {
        first_output {
            path 'hello_workflow'
            mode 'copy'
        }
    }
    ```

Wieder ist die Logik dieselbe wie zuvor.

Dies zeigt dir, dass du die Ausgabeeinstellungen auf einer sehr granularen Ebene kontrollieren kannst, für jede einzelne Ausgabe.
Versuche gerne, die Pfade oder den Veröffentlichungsmodus für einen der Processes zu ändern, um zu sehen, was passiert.

Natürlich bedeutet das, dass wir hier einige Informationen wiederholen, was unpraktisch werden könnte, wenn wir den Speicherort für alle Ausgaben auf dieselbe Weise aktualisieren wollten.
Später im Kurs lernst du, wie du diese Einstellungen für mehrere Ausgaben auf strukturierte Weise konfigurierst.

### 1.6. Den Workflow mit `-resume` ausführen

Lass uns das mit dem `-resume`-Flag testen, da wir den ersten Schritt des Workflows bereits erfolgreich ausgeführt haben.

```bash
nextflow run hello-workflow.nf -resume
```

??? success "Befehlsausgabe"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-workflow.nf` [high_cantor] DSL2 - revision: d746983511

    executor >  local (3)
    [ab/816321] process > sayHello (3)       [100%] 3 of 3, cached: 3 ✔
    [e0/ecf81b] process > convertToUpper (3) [100%] 3 of 3 ✔
    ```

In der Konsolenausgabe gibt es jetzt eine zusätzliche Zeile, die dem neuen Process entspricht, den wir gerade hinzugefügt haben.

Du findest die Ausgaben im Verzeichnis `results/hello_workflow`, wie im `output`-Block festgelegt.

??? abstract "Verzeichnisinhalte"

    ```console
    results/hello_workflow/
    ├── Bonjour-output.txt
    ├── Hello-output.txt
    ├── Holà-output.txt
    ├── UPPER-Bonjour-output.txt
    ├── UPPER-Hello-output.txt
    └── UPPER-Holà-output.txt
    ```

Das ist praktisch! Aber es lohnt sich trotzdem, einen Blick in das Work-Verzeichnis eines der Aufrufe des zweiten Process zu werfen.

??? abstract "Verzeichnisinhalte"

    ```console
    work/e0/ecf81b4cacc648b9b994218d5b29d7/
    ├── Holà-output.txt -> /workspaces/training/hello-nextflow/work/ab/81632178cd37e9e815959278808819/Holà-output.txt
    └── UPPER-Holà-output.txt
    ```

Beachte, dass es zwei `*-output`-Dateien gibt: die Ausgabe des ersten Process sowie die Ausgabe des zweiten.

Die Ausgabe des ersten Process ist dort, weil Nextflow sie dort **gestageed** hat, um alles Notwendige für die Ausführung innerhalb desselben Unterverzeichnisses zu haben.

Es handelt sich jedoch tatsächlich um einen symbolischen Link, der auf die Originaldatei im Unterverzeichnis des ersten Process-Aufrufs zeigt.
Standardmäßig verwendet Nextflow, wenn es auf einer einzelnen Maschine wie hier ausgeführt wird, symbolische Links anstelle von Kopien, um Eingabe- und Zwischendateien zu stagen.

Bevor du weitermachst, denke darüber nach, wie wir nur die Ausgabe von `sayHello` mit der Eingabe von `convertToUpper` verbunden haben und die beiden Processes nacheinander ausgeführt werden konnten.
Nextflow hat die harte Arbeit erledigt, einzelne Eingabe- und Ausgabedateien zu handhaben und sie zwischen den beiden Befehlen für uns weiterzugeben.

Das ist einer der Gründe, warum Nextflow Kanäle so mächtig sind: Sie kümmern sich um die Routinearbeit, die beim Verbinden von Workflow-Schritten anfällt.

### Fazit

Du weißt, wie du Processes verketten kannst, indem du die Ausgabe eines Schritts als Eingabe für den nächsten Schritt bereitstellst.

### Wie geht es weiter?

Lerne, wie du Ausgaben von gestapelten Process-Aufrufen sammelst und sie in einen einzelnen Process einspeist.

---

## 2. Einen dritten Schritt hinzufügen, um alle Begrüßungen zu sammeln

Wenn wir einen Process verwenden, um eine Transformation auf jedes der Elemente in einem Kanal anzuwenden, wie wir es hier mit den mehreren Begrüßungen tun, wollen wir manchmal Elemente aus dem Ausgabekanal dieses Process sammeln und sie in einen anderen Process einspeisen, der eine Art Analyse oder Zusammenfassung durchführt.

Zur Demonstration werden wir einen neuen Schritt zu unserer Pipeline hinzufügen, der alle vom `convertToUpper`-Process produzierten Großbuchstaben-Begrüßungen sammelt und sie in eine einzelne Datei schreibt.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-collect.svg"
</figure>

Um die Überraschung nicht zu verderben, aber das wird einen sehr nützlichen Operator beinhalten.

### 2.1. Den Sammelbefehl definieren und im Terminal testen

Der Sammelschritt, den wir zu unserem Workflow hinzufügen wollen, wird den `cat`-Befehl verwenden, um mehrere Großbuchstaben-Begrüßungen in eine einzelne Datei zu verketten.

Lass uns den Befehl für sich allein im Terminal ausführen, um zu überprüfen, dass er wie erwartet funktioniert, genau wie wir es zuvor getan haben.

Führe Folgendes in deinem Terminal aus:

```bash
echo 'Hello' | tr '[a-z]' '[A-Z]' > UPPER-Hello-output.txt
echo 'Bonjour' | tr '[a-z]' '[A-Z]' > UPPER-Bonjour-output.txt
echo 'Holà' | tr '[a-z]' '[A-Z]' > UPPER-Holà-output.txt
cat UPPER-Hello-output.txt UPPER-Bonjour-output.txt UPPER-Holà-output.txt > COLLECTED-output.txt
```

Die Ausgabe ist eine Textdatei namens `COLLECTED-output.txt`, die die Großbuchstaben-Versionen der ursprünglichen Begrüßungen enthält.

??? abstract "Dateiinhalte"

    ```console title="COLLECTED-output.txt"
    HELLO
    BONJOUR
    HOLà
    ```

Das ist das Ergebnis, das wir mit unserem Workflow erreichen wollen.

### 2.2. Einen neuen Process für den Sammelschritt erstellen

Lass uns einen neuen Process erstellen und ihn `collectGreetings()` nennen.
Wir können beginnen, ihn basierend auf dem zu schreiben, was wir bisher gesehen haben.

#### 2.2.1. Die 'offensichtlichen' Teile des Process schreiben

Füge die folgende Process-Definition zum Workflow-Skript hinzu:

```groovy title="hello-workflow.nf" linenums="37"
/*
 * Großbuchstaben-Begrüßungen in einer einzelnen Ausgabedatei sammeln
 */
process collectGreetings {

    input:
    ???

    output:
    path "COLLECTED-output.txt"

    script:
    """
    cat ??? > 'COLLECTED-output.txt'
    """
}
```

Das ist, was wir mit Zuversicht basierend auf dem, was du bisher gelernt hast, schreiben können.
Aber das ist nicht funktionsfähig!
Es lässt die Eingabedefinition(en) und die erste Hälfte des Skriptbefehls aus, weil wir herausfinden müssen, wie wir das schreiben.

#### 2.2.2. Eingaben für `collectGreetings()` definieren

Wir müssen die Begrüßungen von allen Aufrufen des `convertToUpper()`-Process sammeln.
Was wissen wir, dass wir vom vorherigen Schritt im Workflow bekommen können?

Der von `convertToUpper()` ausgegebene Kanal wird die Pfade zu den einzelnen Dateien enthalten, die die Großbuchstaben-Begrüßungen enthalten.
Das entspricht einem Eingabeslot; nennen wir ihn der Einfachheit halber `input_files`.

Nimm im Process-Block die folgende Codeänderung vor:

=== "Nachher"

    ```groovy title="hello-workflow.nf" linenums="42" hl_lines="2"
          input:
          path input_files
    ```

=== "Vorher"

    ```groovy title="hello-workflow.nf" linenums="42" hl_lines="2"
          input:
          ???
    ```

Beachte, dass wir das `path`-Präfix verwenden, obwohl wir erwarten, dass dies mehrere Dateien enthält.

#### 2.2.3. Den Verkettungsbefehl zusammenstellen

Hier könnten die Dinge etwas knifflig werden, weil wir in der Lage sein müssen, eine beliebige Anzahl von Eingabedateien zu handhaben.
Konkret können wir den Befehl nicht im Voraus schreiben, also müssen wir Nextflow sagen, wie es ihn zur Laufzeit basierend auf den Eingaben, die in den Process fließen, zusammenstellen soll.

Mit anderen Worten, wenn wir einen Eingabekanal haben, der das Element `[file1.txt, file2.txt, file3.txt]` enthält, brauchen wir Nextflow, um das in `cat file1.txt file2.txt file3.txt` umzuwandeln.

Glücklicherweise macht Nextflow das gerne für uns, wenn wir einfach `cat ${input_files}` in den Skriptbefehl schreiben.

Nimm im Process-Block die folgende Codeänderung vor:

=== "Nachher"

    ```groovy title="hello-workflow.nf" linenums="54" hl_lines="3"
        script:
        """
        cat ${input_files} > 'COLLECTED-output.txt'
        """
    ```

=== "Vorher"

    ```groovy title="hello-workflow.nf" linenums="54"
        script:
        """
        cat ??? > 'COLLECTED-output.txt'
        """
    ```

Theoretisch sollte dies jede beliebige Anzahl von Eingabedateien handhaben können.

!!! tip "Tipp"

    Einige Kommandozeilenwerkzeuge erfordern die Angabe eines Arguments (wie `-input`) für jede Eingabedatei.
    In diesem Fall müssten wir ein bisschen zusätzliche Arbeit leisten, um den Befehl zusammenzustellen.
    Du kannst ein Beispiel dafür im [Nextflow for Genomics](../../nf4_science/genomics/)-Trainingskurs sehen.

### 2.3. Den Sammelschritt zum Workflow hinzufügen

Jetzt sollten wir nur noch den Sammel-Process auf der Ausgabe des Großbuchstaben-Schritts aufrufen müssen.
Das ist auch ein Kanal, genannt `convertToUpper.out`.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-collect-connector.svg"
</figure>

#### 2.3.1. Die Process-Aufrufe verbinden

Nimm im Workflow-Block die folgende Codeänderung vor:

=== "Nachher"

    ```groovy title="hello-workflow.nf" linenums="75" hl_lines="4 5"
        // Die Begrüßung in Großbuchstaben umwandeln
        convertToUpper(sayHello.out)

        // Alle Begrüßungen in einer Datei sammeln
        collectGreetings(convertToUpper.out)
    }
    ```

=== "Vorher"

    ```groovy title="hello-workflow.nf" linenums="75"
        // Die Begrüßung in Großbuchstaben umwandeln
        convertToUpper(sayHello.out)
    }
    ```

Dies verbindet die Ausgabe von `convertToUpper()` mit der Eingabe von `collectGreetings()`.

#### 2.3.2. Den Workflow mit `-resume` ausführen

Lass uns es versuchen.

```bash
nextflow run hello-workflow.nf -resume
```

??? success "Befehlsausgabe"

    ```console hl_lines="8"
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-workflow.nf` [mad_gilbert] DSL2 - revision: 6acfd5e28d

    executor >  local (3)
    [79/33b2f0] sayHello (2)         | 3 of 3, cached: 3 ✔
    [99/79394f] convertToUpper (3)   | 3 of 3, cached: 3 ✔
    [47/50fe4a] collectGreetings (1) | 3 of 3 ✔
    ```

Es läuft erfolgreich, einschließlich des dritten Schritts.

Schau dir jedoch die Anzahl der Aufrufe für `collectGreetings()` in der letzten Zeile an.
Wir haben nur einen erwartet, aber es sind drei.

Schau dir jetzt den Inhalt der endgültigen Ausgabedatei an.

??? abstract "Dateiinhalte"

    ```console title="results/COLLECTED-output.txt"
    Holà
    ```

Oh nein. Der Sammelschritt wurde einzeln auf jede Begrüßung ausgeführt, was NICHT das ist, was wir wollten.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-collect-no-operator.svg"
</figure>

Wir müssen etwas tun, um Nextflow explizit zu sagen, dass wir wollen, dass dieser dritte Schritt auf allen Elementen im vom `convertToUpper()` ausgegebenen Kanal läuft.

### 2.4. Einen Operator verwenden, um die Begrüßungen in eine einzelne Eingabe zu sammeln

Ja, wieder einmal ist die Antwort auf unser Problem ein Operator.

Konkret werden wir den treffend benannten [`collect()`](https://nextflow.io/docs/latest/reference/operator.html#collect)-Operator verwenden.

#### 2.4.1. Den `collect()`-Operator hinzufügen

Diesmal wird es etwas anders aussehen, weil wir den Operator nicht im Kontext einer Channel Factory hinzufügen; wir fügen ihn zu einem Ausgabekanal hinzu.

Wir nehmen das `convertToUpper.out` und hängen den `collect()`-Operator an, was uns `convertToUpper.out.collect()` gibt.
Wir können das direkt in den `collectGreetings()`-Process-Aufruf einstecken.

Nimm im Workflow-Block die folgende Codeänderung vor:

=== "Nachher"

    ```groovy title="hello-workflow.nf" linenums="73" hl_lines="2"
        // Alle Begrüßungen in einer Datei sammeln
        collectGreetings(convertToUpper.out.collect())
    }
    ```

=== "Vorher"

    ```groovy title="hello-workflow.nf" linenums="73" hl_lines="2"
        // Alle Begrüßungen in einer Datei sammeln
        collectGreetings(convertToUpper.out)
    }
    ```

#### 2.4.2. Einige `view()`-Anweisungen hinzufügen

Lass uns auch ein paar `view()`-Anweisungen einfügen, um die Zustände des Kanal-Inhalts vor und nach zu visualisieren.

=== "Nachher"

    ```groovy title="hello-workflow.nf" linenums="73" hl_lines="4-6"
        // Alle Begrüßungen in einer Datei sammeln
        collectGreetings(convertToUpper.out.collect())

        // Optionale view-Anweisungen
        convertToUpper.out.view { contents -> "Before collect: $contents" }
        convertToUpper.out.collect().view { contents -> "After collect: $contents" }
    }
    ```

=== "Vorher"

    ```groovy title="hello-workflow.nf" linenums="73"
        // Alle Begrüßungen in einer Datei sammeln
        collectGreetings(convertToUpper.out.collect())
    }
    ```

Die `view()`-Anweisungen können überall platziert werden, wo du möchtest; wir haben sie für die Lesbarkeit direkt nach dem Aufruf platziert.

#### 2.4.3. Den Workflow erneut mit `-resume` ausführen

Lass uns es versuchen:

```bash
nextflow run hello-workflow.nf -resume
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-workflow.nf` [soggy_franklin] DSL2 - revision: bc8e1b2726

    [d6/cdf466] sayHello (1)       | 3 of 3, cached: 3 ✔
    [99/79394f] convertToUpper (2) | 3 of 3, cached: 3 ✔
    [1e/83586c] collectGreetings   | 1 of 1 ✔
    Before collect: /workspaces/training/hello-nextflow/work/b3/d52708edba8b864024589285cb3445/UPPER-Bonjour-output.txt
    Before collect: /workspaces/training/hello-nextflow/work/99/79394f549e3040dfc2440f69ede1fc/UPPER-Hello-output.txt
    Before collect: /workspaces/training/hello-nextflow/work/aa/56bfe7cf00239dc5badc1d04b60ac4/UPPER-Holà-output.txt
    After collect: [/workspaces/training/hello-nextflow/work/b3/d52708edba8b864024589285cb3445/UPPER-Bonjour-output.txt, /workspaces/training/hello-nextflow/work/99/79394f549e3040dfc2440f69ede1fc/UPPER-Hello-output.txt, /workspaces/training/hello-nextflow/work/aa/56bfe7cf00239dc5badc1d04b60ac4/UPPER-Holà-output.txt]
    ```

Es läuft erfolgreich, obwohl die Log-Ausgabe möglicherweise etwas unordentlicher aussieht als dies (wir haben es für die Lesbarkeit aufgeräumt).

Diesmal wurde der dritte Schritt nur einmal aufgerufen!
Wenn du dir die Ausgabe der `view()`-Anweisungen anschaust, siehst du Folgendes:

- Drei `Before collect:`-Anweisungen, eine für jede Begrüßung: zu diesem Zeitpunkt sind die Dateipfade einzelne Elemente im Kanal.
- Eine einzelne `After collect:`-Anweisung: die drei Dateipfade sind jetzt in ein einzelnes Element verpackt.

Wir können das mit dem folgenden Diagramm zusammenfassen:

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-collect-WITH-operator.svg"
</figure>

Zum Schluss kannst du dir den Inhalt der Ausgabedatei anschauen, um dich selbst davon zu überzeugen, dass alles korrekt funktioniert hat.

??? abstract "Dateiinhalte"

    ```console title="results/COLLECTED-output.txt"
    BONJOUR
    HELLO
    HOLà
    ```

Diesmal haben wir alle drei Begrüßungen in der endgültigen Ausgabedatei. Erfolg!

!!! note "Hinweis"

    Wenn du dies mehrmals ohne `-resume` ausführst, wirst du sehen, dass sich die Reihenfolge der Begrüßungen von einem Lauf zum nächsten ändert.
    Dies zeigt dir, dass die Reihenfolge, in der Elemente durch Process-Aufrufe fließen, nicht garantiert konsistent ist.

#### 2.4.4. `view()`-Anweisungen für die Lesbarkeit entfernen

Bevor du zum nächsten Abschnitt übergehst, empfehlen wir, die `view()`-Anweisungen zu löschen, um die Konsolenausgabe nicht zu überladen.

=== "Nachher"

    ```groovy title="hello-workflow.nf" linenums="73"
        // Alle Begrüßungen in einer Datei sammeln
        collectGreetings(convertToUpper.out.collect())
    ```

=== "Vorher"

    ```groovy title="hello-workflow.nf" linenums="73" hl_lines="4-6"
        // Alle Begrüßungen in einer Datei sammeln
        collectGreetings(convertToUpper.out.collect())

        // Optionale view-Anweisungen
        convertToUpper.out.view { contents -> "Before collect: $contents" }
        convertToUpper.out.collect().view { contents -> "After collect: $contents" }
    ```

Dies ist im Grunde die umgekehrte Operation von Punkt 2.4.2.

### Fazit

Du weißt, wie du Ausgaben aus einem Stapel von Process-Aufrufen sammelst und sie in einen gemeinsamen Analyse- oder Zusammenfassungsschritt einspeist.

Zur Wiederholung, das ist, was du bisher erstellt hast:

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-collect.svg"
</figure>

### Wie geht es weiter?

Lerne, wie du zusätzliche Parameter an einen Process übergibst.

---

## 3. Zusätzliche Parameter an einen Process übergeben

Wir wollen in der Lage sein, der endgültigen Ausgabedatei einen bestimmten Namen zu geben, um nachfolgende Stapel von Begrüßungen verarbeiten zu können, ohne die endgültigen Ergebnisse zu überschreiben.

Zu diesem Zweck werden wir die folgenden Verfeinerungen am Workflow vornehmen:

- Den Sammel-Process so modifizieren, dass er einen benutzerdefinierten Namen für die Ausgabedatei akzeptiert (`batch_name`)
- Einen Kommandozeilenparameter zum Workflow hinzufügen (`--batch`) und ihn an den Sammel-Process übergeben

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-collect-batch.svg"
</figure>

### 3.1. Den Sammel-Process modifizieren

Wir müssen die zusätzliche Eingabe deklarieren und sie in den Ausgabedateinamen integrieren.

#### 3.1.1. Die zusätzliche Eingabe deklarieren

Gute Nachricht: Wir können so viele Eingabevariablen deklarieren, wie wir wollen, in der Process-Definition.
Nennen wir diese `batch_name`.

Nimm im Process-Block die folgende Codeänderung vor:

=== "Nachher"

    ```groovy title="hello-workflow.nf" linenums="42" hl_lines="3"
        input:
        path input_files
        val batch_name
    ```

=== "Vorher"

    ```groovy title="hello-workflow.nf" linenums="42"
        input:
        path input_files
    ```

Du kannst deine Processes so einrichten, dass sie so viele Eingaben erwarten, wie du möchtest.
Im Moment sind diese alle als erforderliche Eingaben eingerichtet; du _musst_ einen Wert angeben, damit der Workflow funktioniert.

Wie man erforderliche und optionale Eingaben handhabt, lernst du später auf deiner Nextflow-Reise.

#### 3.1.2. Die `batch_name`-Variable im Ausgabedateinamen verwenden

Wir können die Variable auf dieselbe Weise in den Ausgabedateinamen einfügen, wie wir zuvor dynamische Dateinamen zusammengestellt haben.

Nimm im Process-Block die folgende Codeänderung vor:

=== "Nachher"

    ```groovy title="hello-workflow.nf" linenums="46" hl_lines="2 6"
        output:
        path "COLLECTED-${batch_name}-output.txt"

        script:
        """
        cat ${input_files} > 'COLLECTED-${batch_name}-output.txt'
        """
    ```

=== "Vorher"

    ```groovy title="hello-workflow.nf" linenums="46" hl_lines="2 6"
        output:
        path "COLLECTED-output.txt"

        script:
        """
        cat ${input_files} > 'COLLECTED-output.txt'
        """
    ```

Dies richtet den Process ein, den `batch_name`-Wert zu verwenden, um einen spezifischen Dateinamen für die endgültige Ausgabe des Workflows zu generieren.

### 3.2. Einen `batch`-Kommandozeilenparameter hinzufügen

Jetzt brauchen wir einen Weg, den Wert für `batch_name` zu liefern und ihn an den Process-Aufruf zu übergeben.

#### 3.2.1. `params` verwenden, um den Parameter einzurichten

Du weißt bereits, wie du das `params`-System verwendest, um CLI-Parameter zu deklarieren.
Lass uns das verwenden, um einen `batch`-Parameter zu deklarieren (mit einem Standardwert, weil wir faul sind).

Nimm im Abschnitt für Pipeline-Parameter die folgenden Codeänderungen vor:

=== "Nachher"

    ```groovy title="hello-workflow.nf" linenums="55" hl_lines="6"
    /*
     * Pipeline-Parameter
     */
    params {
        input: Path = 'data/greetings.csv'
        batch: String = 'batch'
    }
    ```

=== "Vorher"

    ```groovy title="hello-workflow.nf" linenums="55"
    /*
     * Pipeline-Parameter
     */
    params {
        input: Path = 'data/greetings.csv'
    }
    ```

Genau wie wir es für `--input` demonstriert haben, kannst du diesen Standardwert überschreiben, indem du einen Wert mit `--batch` auf der Kommandozeile angibst.

#### 3.2.2. Den `batch`-Parameter an den Process übergeben

Um den Wert des Parameters an den Process zu übergeben, müssen wir ihn im Process-Aufruf hinzufügen.

Nimm im Workflow-Block die folgende Codeänderung vor:

=== "Nachher"

    ```groovy title="hello-workflow.nf" linenums="74" hl_lines="2"
        // Alle Begrüßungen in einer Datei sammeln
        collectGreetings(convertToUpper.out.collect(), params.batch)
    ```

=== "Vorher"

    ```groovy title="hello-workflow.nf" linenums="74" hl_lines="2"
        // Alle Begrüßungen in einer Datei sammeln
        collectGreetings(convertToUpper.out.collect())
    ```

Du siehst, dass um mehrere Eingaben an einen Process zu übergeben, du sie einfach in den Aufruf-Klammern auflistest, durch Kommas getrennt.

!!! warning "Warnung"

    Du MUSST die Eingaben an den Process in der EXAKT GLEICHEN REIHENFOLGE übergeben, wie sie im Eingabedefinitionsblock des Process aufgelistet sind.

### 3.3. Den Workflow ausführen

Lass uns versuchen, dies mit einem Stapelnamen auf der Kommandozeile auszuführen.

```bash
nextflow run hello-workflow.nf -resume --batch trio
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-workflow.nf` [confident_rutherford] DSL2 - revision: bc58af409c

    executor >  local (1)
    [79/33b2f0] sayHello (2)       | 3 of 3, cached: 3 ✔
    [99/79394f] convertToUpper (2) | 3 of 3, cached: 3 ✔
    [b5/f19efe] collectGreetings   | 1 of 1 ✔
    ```

Es läuft erfolgreich und produziert die gewünschte Ausgabe:

??? abstract "Dateiinhalte"

    ```console title="results/COLLECTED-trio-output.txt"
    HELLO
    BONJOUR
    HOLà
    ```

Solange wir den Parameter entsprechend angeben, werden nachfolgende Läufe auf anderen Eingabestapeln keine vorherigen Ergebnisse überschreiben.

### Fazit

Du weißt, wie du mehr als eine Eingabe an einen Process übergibst.

### Wie geht es weiter?

Lerne, wie du mehrere Ausgaben ausgibst und sie bequem handhabst.

---

## 4. Eine Ausgabe zum Sammelschritt hinzufügen

Bisher haben wir Processes verwendet, die jeweils nur eine Ausgabe produzierten.
Wir konnten auf ihre jeweiligen Ausgaben sehr bequem mit der `<process>.out`-Syntax zugreifen, die wir sowohl im Kontext der Übergabe einer Ausgabe an den nächsten Process (z.B. `convertToUpper(sayHello.out)`) als auch im Kontext des `publish:`-Abschnitts (z.B. `first_output = sayHello.out`) verwendet haben.

Was passiert, wenn ein Process mehr als eine Ausgabe produziert?
Wie handhaben wir die mehreren Ausgaben?
Können wir eine bestimmte Ausgabe auswählen und verwenden?

Alles ausgezeichnete Fragen, und die kurze Antwort ist ja, das können wir!

Mehrere Ausgaben werden in separate Kanäle verpackt.
Wir können entweder wählen, diesen Ausgabekanälen Namen zu geben, was es einfach macht, sie später einzeln zu referenzieren, oder wir können sie per Index referenzieren.

Zur Demonstration sagen wir, dass wir die Anzahl der Begrüßungen zählen wollen, die für einen gegebenen Eingabestapel gesammelt werden, und sie in einer Datei berichten wollen.

### 4.1. Den Process so modifizieren, dass er die Anzahl der Begrüßungen zählt und ausgibt

Dies erfordert zwei wichtige Änderungen an der Process-Definition: wir brauchen einen Weg, die Begrüßungen zu zählen und eine Berichtsdatei zu schreiben, dann müssen wir diese Berichtsdatei zum `output`-Block des Process hinzufügen.

#### 4.1.1. Die Anzahl der gesammelten Begrüßungen zählen

Praktischerweise erlaubt Nextflow uns, beliebigen Code im `script:`-Block der Process-Definition hinzuzufügen, was wirklich praktisch ist, um solche Dinge zu tun.

Das bedeutet, wir können Nextflows eingebaute `size()`-Funktion verwenden, um die Anzahl der Dateien im `input_files`-Array zu erhalten, und das Ergebnis mit einem `echo`-Befehl in eine Datei schreiben.

Nimm im `collectGreetings`-Process-Block die folgenden Codeänderungen vor:

=== "Nachher"

    ```groovy title="hello-workflow.nf" linenums="55" hl_lines="2 5"
        script:
        count_greetings = input_files.size()
        """
        cat ${input_files} > 'COLLECTED-${batch_name}-output.txt'
        echo 'There were ${count_greetings} greetings in this batch.' > '${batch_name}-report.txt'
        """
    ```

=== "Vorher"

    ```groovy title="hello-workflow.nf" linenums="55"
        script:
        """
        cat ${input_files} > 'COLLECTED-${batch_name}-output.txt'
        """
    ```

Die `count_greetings`-Variable wird zur Laufzeit berechnet.

#### 4.1.2. Die Berichtsdatei ausgeben und Ausgaben benennen

Im Prinzip müssen wir nur die Berichtsdatei zum `output:`-Block hinzufügen.

Aber während wir dabei sind, werden wir auch einige `emit:`-Tags zu unseren Ausgabedeklarationen hinzufügen. Diese werden es uns ermöglichen, die Ausgaben per Namen auszuwählen, anstatt Positionsindizes verwenden zu müssen.

Nimm im Process-Block die folgende Codeänderung vor:

=== "Nachher"

    ```groovy title="hello-workflow.nf" linenums="46" hl_lines="2 3"
        output:
        path "COLLECTED-${batch_name}-output.txt", emit: outfile
        path "${batch_name}-report.txt", emit: report
    ```

=== "Vorher"

    ```groovy title="hello-workflow.nf" linenums="46"
        output:
        path "COLLECTED-${batch_name}-output.txt"
    ```

Die `emit:`-Tags sind optional, und wir hätten ein Tag nur zu einer der Ausgaben hinzufügen können.
Aber wie das Sprichwort sagt, warum nicht beides?

!!! tip "Tipp"

    Wenn du die Ausgaben eines Process nicht mit `emit:` benennst, kannst du trotzdem einzeln auf sie zugreifen, indem du ihren jeweiligen (nullbasierten) Index verwendest.
    Zum Beispiel würdest du `<process>.out[0]` verwenden, um die erste Ausgabe zu bekommen, `<process>.out[1]` für die zweite Ausgabe, und so weiter.

    Wir bevorzugen es, Ausgaben zu benennen, weil es sonst zu einfach ist, versehentlich den falschen Index zu greifen, besonders wenn der Process viele Ausgaben produziert.

### 4.2. Die Workflow-Ausgaben aktualisieren

Jetzt, da wir zwei Ausgaben aus dem `collectGreetings`-Process haben, enthält die `collectGreetings.out`-Ausgabe zwei Kanäle:

- `collectGreetings.out.outfile` enthält die endgültige Ausgabedatei
- `collectGreetings.out.report` enthält die Berichtsdatei

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-collect-report.svg"
</figure>

Wir müssen die Workflow-Ausgaben entsprechend aktualisieren.

#### 4.2.1. Den `publish:`-Abschnitt aktualisieren

Nimm im `workflow`-Block die folgende Codeänderung vor:

=== "Nachher"

    ```groovy title="hello-workflow.nf" linenums="80" hl_lines="4 5"
        publish:
        first_output = sayHello.out
        uppercased = convertToUpper.out
        collected = collectGreetings.out.outfile
        batch_report = collectGreetings.out.report
    ```

=== "Vorher"

    ```groovy title="hello-workflow.nf" linenums="80" hl_lines="4"
        publish:
        first_output = sayHello.out
        uppercased = convertToUpper.out
        collected = collectGreetings.out
    ```

Wie du sehen kannst, ist das Referenzieren spezifischer Process-Ausgaben jetzt trivial.
Wenn wir in Teil 5 (Containers) einen weiteren Schritt zu unserer Pipeline hinzufügen, werden wir in der Lage sein, leicht auf `collectGreetings.out.outfile` zu verweisen und es dem neuen Process zu übergeben (Spoiler: der neue Process heißt `cowpy`).

Aber für jetzt aktualisieren wir die Workflow-Level-Ausgaben fertig.

#### 4.2.2. Den `output`-Block aktualisieren

Nimm im `output`-Block die folgende Codeänderung vor:

=== "Nachher"

    ```groovy title="hello-workflow.nf" linenums="86" hl_lines="14-17"
    output {
        first_output {
            path 'hello_workflow'
            mode 'copy'
        }
        uppercased {
            path 'hello_workflow'
            mode 'copy'
        }
        collected {
            path 'hello_workflow'
            mode 'copy'
        }
        batch_report {
            path 'hello_workflow'
            mode 'copy'
        }
    }
    ```

=== "Vorher"

    ```groovy title="hello-workflow.nf" linenums="80"
    output {
        first_output {
            path 'hello_workflow'
            mode 'copy'
        }
        uppercased {
            path 'hello_workflow'
            mode 'copy'
        }
        collected {
            path 'hello_workflow'
            mode 'copy'
        }
    }
    ```

Wir müssen die `collected`-Ausgabedefinition nicht aktualisieren, da sich dieser Name nicht geändert hat.
Wir müssen nur die neue Ausgabe hinzufügen.

### 4.3. Den Workflow ausführen

Lass uns versuchen, dies mit dem aktuellen Stapel von Begrüßungen auszuführen.

```bash
nextflow run hello-workflow.nf -resume --batch trio
```

??? success "Befehlsausgabe"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-workflow.nf` [ecstatic_wilson] DSL2 - revision: c80285f8c8

    executor >  local (1)
    [c5/4c6ca9] sayHello (3)       [100%] 3 of 3, cached: 3 ✔
    [0e/6cbc59] convertToUpper (3) [100%] 3 of 3, cached: 3 ✔
    [02/61ead2] collectGreetings   [100%] 1 of 1 ✔
    ```

Wenn du ins Verzeichnis `results/hello_workflow/` schaust, findest du die neue Berichtsdatei `trio-report.txt`.
Öffne sie, um zu überprüfen, dass der Workflow korrekt die Anzahl der verarbeiteten Begrüßungen berichtet hat.

??? abstract "Dateiinhalte"

    ```txt title="trio-report.txt"
    There were 3 greetings in this batch.
    ```

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-collect-4-way.svg"
</figure>

Du kannst gerne weitere Begrüßungen zur CSV hinzufügen und testen, was passiert.

### Fazit

Du weißt, wie du einen Process mehrere benannte Ausgaben ausgeben lässt und wie du sie auf Workflow-Ebene angemessen handhabst.

Allgemeiner verstehst du die wichtigsten Prinzipien, die beim Verbinden von Processes auf gängige Weise involviert sind.

### Wie geht es weiter?

Mach eine extra lange Pause, du hast sie verdient.

Wenn du bereit bist, gehe zu [**Teil 4: Hello Modules**](./04_hello_modules.md), um zu lernen, wie du deinen Code für bessere Wartbarkeit und Code-Effizienz modularisierst.

---

## Quiz

<quiz>
Wie greifst du auf die Ausgabe eines Process im Workflow-Block zu?
- [ ] `process.output`
- [ ] `output.processName`
- [x] `processName.out`
- [ ] `get(processName)`

Mehr erfahren: [1.4. Die Ausgabe des ersten Process an den zweiten Process übergeben](#14-die-ausgabe-des-ersten-process-an-den-zweiten-process-ubergeben)
</quiz>

<quiz>
Was bestimmt die Reihenfolge der Process-Ausführung in Nextflow?
- [ ] Die Reihenfolge, in der Processes im Workflow-Block geschrieben sind
- [ ] Alphabetische Reihenfolge nach Process-Namen
- [x] Datenabhängigkeiten zwischen Processes
- [ ] Zufällige Reihenfolge für parallele Ausführung

Mehr erfahren: [1.4. Die Ausgabe des ersten Process an den zweiten Process übergeben](#14-die-ausgabe-des-ersten-process-an-den-zweiten-process-ubergeben)
</quiz>

<quiz>
Welcher Operator sollte `???` ersetzen, um alle Ausgaben in eine einzelne Liste für den nachfolgenden Process zu sammeln?

```groovy hl_lines="4"
workflow {
    greetings_ch = Channel.of('Hello', 'Bonjour', 'Hola')
    SAYHELLO(greetings_ch)
    GATHER_ALL(SAYHELLO.out.???)
}
```

- [ ] `flatten()`
- [x] `collect()`
- [ ] `mix()`
- [ ] `join()`

Mehr erfahren: [2.4. Einen Operator verwenden, um die Begrüßungen in eine einzelne Eingabe zu sammeln](#24-einen-operator-verwenden-um-die-begrussungen-in-eine-einzelne-eingabe-zu-sammeln)
</quiz>

<quiz>
Wann solltest du den `collect()`-Operator verwenden?
- [ ] Wenn du Elemente parallel verarbeiten möchtest
- [ ] Wenn du Kanal-Inhalte filtern musst
- [x] Wenn ein nachfolgender Process alle Elemente von einem vorgelagerten Process benötigt
- [ ] Wenn du Daten auf mehrere Processes aufteilen möchtest

Mehr erfahren: [2.4. Einen Operator verwenden, um die Begrüßungen in eine einzelne Eingabe zu sammeln](#24-einen-operator-verwenden-um-die-begrussungen-in-eine-einzelne-eingabe-zu-sammeln)
</quiz>

<quiz>
Wie greifst du auf eine benannte Ausgabe von einem Process zu?
- [ ] `processName.outputName`
- [ ] `processName.get(outputName)`
- [x] `processName.out.outputName`
- [ ] `output.processName.outputName`

Mehr erfahren: [4.1.2. Die Berichtsdatei ausgeben und Ausgaben benennen](#412-die-berichtsdatei-ausgeben-und-ausgaben-benennen)
</quiz>

<quiz>
Was ist die korrekte Syntax, um eine Ausgabe in einem Process zu benennen?
- [ ] `name: outputName`
- [ ] `output: outputName`
- [x] `emit: outputName`
- [ ] `label: outputName`

Mehr erfahren: [4.1.2. Die Berichtsdatei ausgeben und Ausgaben benennen](#412-die-berichtsdatei-ausgeben-und-ausgaben-benennen)
</quiz>

<quiz>
Was muss bei der Übergabe mehrerer Eingaben an einen Process zutreffen?
- [ ] Alle Eingaben müssen vom gleichen Typ sein
- [ ] Eingaben müssen in alphabetischer Reihenfolge übergeben werden
- [x] Die Reihenfolge der Eingaben muss mit der im Eingabeblock definierten Reihenfolge übereinstimmen
- [ ] Nur zwei Eingaben können gleichzeitig übergeben werden

Mehr erfahren: [3. Zusätzliche Parameter an einen Process übergeben](#3-zusatzliche-parameter-an-einen-process-ubergeben)
</quiz>
