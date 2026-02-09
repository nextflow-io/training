# Teil 3: Hello Workflow

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/_aO56V3iXGI?si=y8lAedhEHWaTV4zd&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } Siehe [die gesamte Playlist](https://youtube.com/playlist?list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&si=eF7cLR62goy-lc6n) auf dem Nextflow YouTube-Kanal.

:green_book: Das Video-Transkript ist [hier](./transcripts/03_hello_workflow.md) verfügbar.
///

Die meisten realen Workflows umfassen mehr als einen Schritt.
In diesem Trainingsmodul lernst du, wie du Prozesse in einem mehrstufigen Workflow miteinander verbindest.

Du lernst die Nextflow-Methode, um Folgendes zu erreichen:

1. Daten von einem Prozess zum nächsten fließen lassen
2. Ausgaben von mehreren Prozessaufrufen in einem einzigen Prozessaufruf sammeln
3. Zusätzliche Parameter an einen Prozess übergeben
4. Mehrere Ausgaben verarbeiten, die aus einem Prozess kommen

Zur Demonstration bauen wir weiter auf dem domänenunabhängigen Hello World-Beispiel aus Teil 1 und 2 auf.
Diesmal nehmen wir folgende Änderungen an unserem Workflow vor, um besser zu zeigen, wie echte Workflows aufgebaut werden:

1. Einen zweiten Schritt hinzufügen, der die Begrüßung in Großbuchstaben umwandelt.
2. Einen dritten Schritt hinzufügen, der alle umgewandelten Begrüßungen sammelt und in eine einzige Datei schreibt.
3. Einen Parameter hinzufügen, um die finale Ausgabedatei zu benennen, und diesen als zweite Eingabe an den Sammelschritt übergeben.
4. Den Sammelschritt auch eine einfache Statistik über die verarbeiteten Daten ausgeben lassen.

??? info "Wie du von diesem Abschnitt aus beginnst"

    Dieser Abschnitt des Kurses setzt voraus, dass du Teil 1-2 des [Hello Nextflow](./index.md)-Kurses abgeschlossen hast. Wenn du jedoch mit den dort behandelten Grundlagen vertraut bist, kannst du von hier aus beginnen, ohne etwas Besonderes tun zu müssen.

---

## 0. Aufwärmen: `hello-workflow.nf` ausführen

Wir verwenden das Workflow-Skript `hello-workflow.nf` als Ausgangspunkt.
Es entspricht dem Skript, das durch die Bearbeitung von Teil 2 dieses Trainingskurses entsteht, außer dass wir die `view()`-Anweisungen entfernt und das Ausgabeziel geändert haben:

```groovy title="hello-workflow.nf" linenums="37" hl_lines="3"
output {
    first_output {
        path 'hello_workflow'
        mode 'copy'
    }
}
```

Dieses Diagramm fasst die aktuelle Funktionsweise des Workflows zusammen.
Es sollte vertraut aussehen, außer dass wir jetzt explizit zeigen, dass die Ausgaben des Prozesses in einem Kanal verpackt werden, genau wie die Eingaben.
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

??? abstract "Verzeichnisinhalt"

    ```console
    results/hello_workflow
    ├── Bonjour-output.txt
    ├── Hello-output.txt
    └── Holà-output.txt
    ```

Wenn das bei dir funktioniert hat, bist du bereit zu lernen, wie man einen mehrstufigen Workflow zusammenstellt.

---

## 1. Einen zweiten Schritt zum Workflow hinzufügen

Wir fügen einen Schritt hinzu, um jede Begrüßung in Großbuchstaben umzuwandeln.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-multistep.svg"
</figure>

Dazu müssen wir drei Dinge tun:

- Den Befehl definieren, den wir für die Großbuchstabenumwandlung verwenden werden.
- Einen neuen Prozess schreiben, der den Großbuchstabenbefehl umschließt.
- Den neuen Prozess im Workflow-Block aufrufen und so einrichten, dass er die Ausgabe des `sayHello()`-Prozesses als Eingabe nimmt.

### 1.1. Den Großbuchstabenbefehl definieren und im Terminal testen

Für die Umwandlung der Begrüßungen in Großbuchstaben verwenden wir ein klassisches UNIX-Tool namens `tr` für 'text replacement' mit folgender Syntax:

```bash title="Syntax"
tr '[a-z]' '[A-Z]'
```

Dies ist eine sehr naive Textersetzung, die keine Akzentbuchstaben berücksichtigt, sodass zum Beispiel 'Holà' zu 'HOLà' wird, aber es reicht aus, um die Nextflow-Konzepte zu demonstrieren, und darauf kommt es an.

Um es zu testen, können wir den Befehl `echo 'Hello World'` ausführen und seine Ausgabe an den `tr`-Befehl weiterleiten:

```bash
echo 'Hello World' | tr '[a-z]' '[A-Z]' > UPPER-output.txt
```

Die Ausgabe ist eine Textdatei namens `UPPER-output.txt`, die die Großbuchstabenversion des `Hello World`-Strings enthält.

??? abstract "Dateiinhalt"

    ```console title="UPPER-output.txt"
    HELLO WORLD
    ```

Das ist im Grunde das, was wir mit unserem Workflow erreichen wollen.

### 1.2. Den Großbuchstabenschritt als Nextflow-Prozess schreiben

Wir können unseren neuen Prozess am ersten orientieren, da wir alle gleichen Komponenten verwenden wollen.

Füge die folgende Prozessdefinition zum Workflow-Skript hinzu, direkt unter dem ersten:

```groovy title="hello-workflow.nf" linenums="20"
/*
 * Ein Textersetzungstool verwenden, um die Begrüßung in Großbuchstaben umzuwandeln
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

In diesem Fall setzen wir den zweiten Ausgabedateinamen basierend auf dem Eingabedateinamen zusammen, ähnlich wie wir es ursprünglich für die Ausgabe des ersten Prozesses getan haben.

### 1.3. Einen Aufruf des neuen Prozesses im Workflow-Block hinzufügen

Jetzt müssen wir Nextflow mitteilen, dass es den gerade definierten Prozess tatsächlich aufrufen soll.

Nimm im Workflow-Block folgende Codeänderung vor:

=== "Danach"

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

Dies ist noch nicht funktionsfähig, da wir nicht angegeben haben, was in den `convertToUpper()`-Prozess eingegeben werden soll.

### 1.4. Die Ausgabe des ersten Prozesses an den zweiten Prozess übergeben

Jetzt müssen wir die Ausgabe des `sayHello()`-Prozesses in den `convertToUpper()`-Prozess fließen lassen.

Praktischerweise verpackt Nextflow die Ausgabe eines Prozesses automatisch in einen Kanal, wie im Diagramm im Aufwärmabschnitt gezeigt.
Wir können auf den Ausgabekanal eines Prozesses als `<process>.out` verweisen.

Die Ausgabe des `sayHello`-Prozesses ist also ein Kanal namens `sayHello.out`, den wir direkt in den Aufruf von `convertToUpper()` einstecken können.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-multistep-connector.svg"
</figure>

Nimm im Workflow-Block folgende Codeänderung vor:

=== "Danach"

    ```groovy title="hello-workflow.nf" linenums="53" hl_lines="2"
        // Die Begrüßung in Großbuchstaben umwandeln
        convertToUpper(sayHello.out)
    ```

=== "Vorher"

    ```groovy title="hello-workflow.nf" linenums="53" hl_lines="2"
        // Die Begrüßung in Großbuchstaben umwandeln
        convertToUpper()
    ```

Für einen einfachen Fall wie diesen (eine Ausgabe zu einer Eingabe) ist das alles, was wir tun müssen, um zwei Prozesse zu verbinden!

### 1.5. Die Workflow-Ausgabeveröffentlichung einrichten

Lass uns abschließend die Workflow-Ausgaben aktualisieren, um auch die Ergebnisse des zweiten Prozesses zu veröffentlichen.

#### 1.5.1. Den `publish:`-Abschnitt des `workflow`-Blocks aktualisieren

Nimm im `workflow`-Block folgende Codeänderung vor:

=== "Danach"

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

Nimm im `output`-Block folgende Codeänderung vor:

=== "Danach"

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

Auch hier ist die Logik dieselbe wie zuvor.

Dies zeigt dir, dass du die Ausgabeeinstellungen auf sehr granularer Ebene steuern kannst, für jede einzelne Ausgabe.
Du kannst gerne versuchen, die Pfade oder den Veröffentlichungsmodus für einen der Prozesse zu ändern, um zu sehen, was passiert.

Natürlich bedeutet das, dass wir hier einige Informationen wiederholen, was unpraktisch werden könnte, wenn wir den Speicherort für alle Ausgaben auf die gleiche Weise aktualisieren wollten.
Später im Kurs lernst du, wie du diese Einstellungen für mehrere Ausgaben strukturiert konfigurieren kannst.

### 1.6. Den Workflow mit `-resume` ausführen

Lass uns dies mit dem `-resume`-Flag testen, da wir den ersten Schritt des Workflows bereits erfolgreich ausgeführt haben.

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

Es gibt jetzt eine zusätzliche Zeile in der Konsolenausgabe, die dem neuen Prozess entspricht, den wir gerade hinzugefügt haben.

Du findest die Ausgaben im Verzeichnis `results/hello_workflow`, wie im `output`-Block festgelegt.

??? abstract "Verzeichnisinhalt"

    ```console
    results/hello_workflow/
    ├── Bonjour-output.txt
    ├── Hello-output.txt
    ├── Holà-output.txt
    ├── UPPER-Bonjour-output.txt
    ├── UPPER-Hello-output.txt
    └── UPPER-Holà-output.txt
    ```

Das ist praktisch! Aber es lohnt sich trotzdem, einen Blick in das Arbeitsverzeichnis eines der Aufrufe des zweiten Prozesses zu werfen.

??? abstract "Verzeichnisinhalt"

    ```console
    work/e0/ecf81b4cacc648b9b994218d5b29d7/
    ├── Holà-output.txt -> /workspaces/training/hello-nextflow/work/ab/81632178cd37e9e815959278808819/Holà-output.txt
    └── UPPER-Holà-output.txt
    ```

Beachte, dass es zwei `*-output`-Dateien gibt: die Ausgabe des ersten Prozesses sowie die Ausgabe des zweiten.

Die Ausgabe des ersten Prozesses ist dort, weil Nextflow sie dort **bereitgestellt** hat, um alles Notwendige für die Ausführung im selben Unterverzeichnis zu haben.

Es handelt sich jedoch tatsächlich um einen symbolischen Link, der auf die Originaldatei im Unterverzeichnis des ersten Prozessaufrufs zeigt.
Standardmäßig verwendet Nextflow beim Ausführen auf einem einzelnen Rechner, wie wir es hier tun, symbolische Links anstelle von Kopien, um Eingabe- und Zwischendateien bereitzustellen.

Denke nun darüber nach, wie wir nur die Ausgabe von `sayHello` mit der Eingabe von `convertToUpper` verbunden haben und die beiden Prozesse in Serie ausgeführt werden konnten.
Nextflow hat die harte Arbeit übernommen, einzelne Eingabe- und Ausgabedateien zu verarbeiten und sie zwischen den beiden Befehlen für uns weiterzugeben.

Das ist einer der Gründe, warum Nextflow-Kanäle so leistungsstark sind: Sie übernehmen die Routinearbeit beim Verbinden von Workflow-Schritten.

### Fazit

Du weißt, wie du Prozesse miteinander verketten kannst, indem du die Ausgabe eines Schritts als Eingabe für den nächsten Schritt bereitstellst.

### Wie geht es weiter?

Lerne, wie du Ausgaben von Batch-Prozessaufrufen sammelst und in einen einzelnen Prozess einspeist.

---

## 2. Einen dritten Schritt zum Sammeln aller Begrüßungen hinzufügen

Wenn wir einen Prozess verwenden, um eine Transformation auf jedes der Elemente in einem Kanal anzuwenden, wie wir es hier mit den mehreren Begrüßungen tun, möchten wir manchmal Elemente aus dem Ausgabekanal dieses Prozesses sammeln und sie in einen anderen Prozess einspeisen, der eine Art Analyse oder Summierung durchführt.

Zur Demonstration fügen wir unserem Pipeline einen neuen Schritt hinzu, der alle von `convertToUpper` erzeugten Großbuchstaben-Begrüßungen sammelt und in eine einzige Datei schreibt.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-collect.svg"
</figure>

Um die Überraschung nicht zu verderben: Dies wird einen sehr nützlichen Operator beinhalten.

### 2.1. Den Sammelbefehl definieren und im Terminal testen

Der Sammelschritt, den wir unserem Workflow hinzufügen möchten, verwendet den `cat`-Befehl, um mehrere Großbuchstaben-Begrüßungen in eine einzige Datei zu verketten.

Lass uns den Befehl allein im Terminal ausführen, um zu überprüfen, dass er wie erwartet funktioniert, genau wie wir es zuvor getan haben.

Führe Folgendes in deinem Terminal aus:

```bash
echo 'Hello' | tr '[a-z]' '[A-Z]' > UPPER-Hello-output.txt
echo 'Bonjour' | tr '[a-z]' '[A-Z]' > UPPER-Bonjour-output.txt
echo 'Holà' | tr '[a-z]' '[A-Z]' > UPPER-Holà-output.txt
cat UPPER-Hello-output.txt UPPER-Bonjour-output.txt UPPER-Holà-output.txt > COLLECTED-output.txt
```

Die Ausgabe ist eine Textdatei namens `COLLECTED-output.txt`, die die Großbuchstabenversionen der ursprünglichen Begrüßungen enthält.

??? abstract "Dateiinhalt"

    ```console title="COLLECTED-output.txt"
    HELLO
    BONJOUR
    HOLà
    ```

Das ist das Ergebnis, das wir mit unserem Workflow erreichen wollen.

### 2.2. Einen neuen Prozess für den Sammelschritt erstellen

Lass uns einen neuen Prozess erstellen und ihn `collectGreetings()` nennen.
Wir können mit dem Schreiben beginnen, basierend auf dem, was wir bisher gesehen haben.

#### 2.2.1. Die 'offensichtlichen' Teile des Prozesses schreiben

Füge die folgende Prozessdefinition zum Workflow-Skript hinzu:

```groovy title="hello-workflow.nf" linenums="37"
/*
 * Großbuchstaben-Begrüßungen in einer einzigen Ausgabedatei sammeln
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

Das ist, was wir mit Zuversicht schreiben können, basierend auf dem, was du bisher gelernt hast.
Aber das ist nicht funktionsfähig!
Es lässt die Eingabedefinition(en) und die erste Hälfte des Skriptbefehls aus, weil wir herausfinden müssen, wie wir das schreiben.

#### 2.2.2. Eingaben für `collectGreetings()` definieren

Wir müssen die Begrüßungen von allen Aufrufen des `convertToUpper()`-Prozesses sammeln.
Was wissen wir, dass wir vom vorherigen Schritt im Workflow erhalten können?

Der von `convertToUpper()` ausgegebene Kanal enthält die Pfade zu den einzelnen Dateien mit den Großbuchstaben-Begrüßungen.
Das entspricht einem Eingabeslot; nennen wir ihn der Einfachheit halber `input_files`.

Nimm im Prozessblock folgende Codeänderung vor:

=== "Danach"

    ```groovy title="hello-workflow.nf" linenums="42" hl_lines="2"
          input:
          path input_files
    ```

=== "Vorher"

    ```groovy title="hello-workflow.nf" linenums="42" hl_lines="2"
          input:
          ???
    ```

Beachte, dass wir das Präfix `path` verwenden, obwohl wir erwarten, dass dies mehrere Dateien enthält.

#### 2.2.3. Den Verkettungsbefehl zusammensetzen

Hier könnte es etwas knifflig werden, weil wir eine beliebige Anzahl von Eingabedateien verarbeiten müssen.
Konkret können wir den Befehl nicht im Voraus schreiben, also müssen wir Nextflow mitteilen, wie es ihn zur Laufzeit basierend auf den eingehenden Eingaben zusammensetzen soll.

Mit anderen Worten: Wenn wir einen Eingabekanal mit dem Element `[file1.txt, file2.txt, file3.txt]` haben, muss Nextflow daraus `cat file1.txt file2.txt file3.txt` machen.

Glücklicherweise ist Nextflow gerne bereit, das für uns zu tun, wenn wir einfach `cat ${input_files}` im Skriptbefehl schreiben.

Nimm im Prozessblock folgende Codeänderung vor:

=== "Danach"

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

Theoretisch sollte dies jede beliebige Anzahl von Eingabedateien verarbeiten.

!!! tip

    Einige Kommandozeilen-Tools erfordern die Angabe eines Arguments (wie `-input`) für jede Eingabedatei.
    In diesem Fall müssten wir etwas zusätzliche Arbeit leisten, um den Befehl zusammenzusetzen.
    Ein Beispiel dafür findest du im [Nextflow for Genomics](../../nf4_science/genomics/)-Trainingskurs.

### 2.3. Den Sammelschritt zum Workflow hinzufügen

Jetzt sollten wir nur noch den Sammelprozess auf die Ausgabe des Großbuchstabenschritts aufrufen müssen.
Das ist auch ein Kanal, genannt `convertToUpper.out`.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-collect-connector.svg"
</figure>

#### 2.3.1. Die Prozessaufrufe verbinden

Nimm im Workflow-Block folgende Codeänderung vor:

=== "Danach"

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

Lass es uns versuchen.

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

Schaue dir jedoch die Anzahl der Aufrufe für `collectGreetings()` in der letzten Zeile an.
Wir haben nur einen erwartet, aber es sind drei.

Schaue dir nun den Inhalt der finalen Ausgabedatei an.

??? abstract "Dateiinhalt"

    ```console title="results/COLLECTED-output.txt"
    Holà
    ```

Oh nein. Der Sammelschritt wurde einzeln auf jede Begrüßung angewendet, was NICHT das war, was wir wollten.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-collect-no-operator.svg"
</figure>

Wir müssen etwas tun, um Nextflow explizit mitzuteilen, dass wir möchten, dass dieser dritte Schritt auf allen Elementen im Kanal ausgeführt wird, der von `convertToUpper()` ausgegeben wird.

### 2.4. Einen Operator verwenden, um die Begrüßungen in einer einzigen Eingabe zu sammeln

Ja, wieder einmal ist die Antwort auf unser Problem ein Operator.

Konkret werden wir den treffend benannten [`collect()`](https://nextflow.io/docs/latest/reference/operator.html#collect)-Operator verwenden.

#### 2.4.1. Den `collect()`-Operator hinzufügen

Diesmal wird es etwas anders aussehen, weil wir den Operator nicht im Kontext einer Kanalfabrik hinzufügen; wir fügen ihn einem Ausgabekanal hinzu.

Wir nehmen `convertToUpper.out` und hängen den `collect()`-Operator an, was uns `convertToUpper.out.collect()` gibt.
Das können wir direkt in den `collectGreetings()`-Prozessaufruf einstecken.

Nimm im Workflow-Block folgende Codeänderung vor:

=== "Danach"

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

Lass uns auch ein paar `view()`-Anweisungen einfügen, um die Zustände des Kanalinhalts vor und nach der Transformation zu visualisieren.

=== "Danach"

    ```groovy title="hello-workflow.nf" linenums="73" hl_lines="4-6"
        // Alle Begrüßungen in einer Datei sammeln
        collectGreetings(convertToUpper.out.collect())

        // Optionale view-Anweisungen
        convertToUpper.out.view { contents -> "Vor collect: $contents" }
        convertToUpper.out.collect().view { contents -> "Nach collect: $contents" }
    }
    ```

=== "Vorher"

    ```groovy title="hello-workflow.nf" linenums="73"
        // Alle Begrüßungen in einer Datei sammeln
        collectGreetings(convertToUpper.out.collect())
    }
    ```

Die `view()`-Anweisungen können überall platziert werden; wir haben sie direkt nach dem Aufruf platziert, um die Lesbarkeit zu verbessern.

#### 2.4.3. Den Workflow erneut mit `-resume` ausführen

Lass es uns versuchen:

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
    Vor collect: /workspaces/training/hello-nextflow/work/b3/d52708edba8b864024589285cb3445/UPPER-Bonjour-output.txt
    Vor collect: /workspaces/training/hello-nextflow/work/99/79394f549e3040dfc2440f69ede1fc/UPPER-Hello-output.txt
    Vor collect: /workspaces/training/hello-nextflow/work/aa/56bfe7cf00239dc5badc1d04b60ac4/UPPER-Holà-output.txt
    Nach collect: [/workspaces/training/hello-nextflow/work/b3/d52708edba8b864024589285cb3445/UPPER-Bonjour-output.txt, /workspaces/training/hello-nextflow/work/99/79394f549e3040dfc2440f69ede1fc/UPPER-Hello-output.txt, /workspaces/training/hello-nextflow/work/aa/56bfe7cf00239dc5badc1d04b60ac4/UPPER-Holà-output.txt]
    ```

Es läuft erfolgreich, obwohl die Log-Ausgabe etwas unübersichtlicher aussehen kann (wir haben sie zur besseren Lesbarkeit bereinigt).

Diesmal wurde der dritte Schritt nur einmal aufgerufen!
Wenn wir uns die Ausgabe der `view()`-Anweisungen ansehen, sehen wir Folgendes:

- Drei `Vor collect:`-Anweisungen, eine für jede Begrüßung: zu diesem Zeitpunkt sind die Dateipfade einzelne Elemente im Kanal.
- Eine einzige `Nach collect:`-Anweisung: die drei Dateipfade sind jetzt in einem einzigen Element verpackt.

Wir können das mit folgendem Diagramm zusammenfassen:

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-collect-WITH-operator.svg"
</figure>

Schließlich kannst du dir den Inhalt der Ausgabedatei ansehen, um dich zu überzeugen, dass alles korrekt funktioniert hat.

??? abstract "Dateiinhalt"

    ```console title="results/COLLECTED-output.txt"
    BONJOUR
    HELLO
    HOLà
    ```

Diesmal haben wir alle drei Begrüßungen in der finalen Ausgabedatei. Erfolg!

!!! note

    Wenn du dies mehrmals ohne `-resume` ausführst, wirst du sehen, dass sich die Reihenfolge der Begrüßungen von einem Lauf zum nächsten ändert.
    Dies zeigt dir, dass die Reihenfolge, in der Elemente durch Prozessaufrufe fließen, nicht garantiert konsistent ist.

#### 2.4.4. `view()`-Anweisungen zur besseren Lesbarkeit entfernen

Bevor du zum nächsten Abschnitt übergehst, empfehlen wir dir, die `view()`-Anweisungen zu löschen, um die Konsolenausgabe nicht zu überladen.

=== "Danach"

    ```groovy title="hello-workflow.nf" linenums="73"
        // Alle Begrüßungen in einer Datei sammeln
        collectGreetings(convertToUpper.out.collect())
    ```

=== "Vorher"

    ```groovy title="hello-workflow.nf" linenums="73" hl_lines="4-6"
        // Alle Begrüßungen in einer Datei sammeln
        collectGreetings(convertToUpper.out.collect())

        // Optionale view-Anweisungen
        convertToUpper.out.view { contents -> "Vor collect: $contents" }
        convertToUpper.out.collect().view { contents -> "Nach collect: $contents" }
    ```

Dies ist im Grunde die umgekehrte Operation von Punkt 2.4.2.

### Fazit

Du weißt, wie du Ausgaben von einem Batch von Prozessaufrufen sammelst und in einen gemeinsamen Analyse- oder Summierungsschritt einspeist.

Zur Zusammenfassung, das hast du bisher aufgebaut:

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-collect.svg"
</figure>

### Wie geht es weiter?

Lerne, wie du mehr als eine Eingabe an einen Prozess übergibst.

---

## 3. Zusätzliche Parameter an einen Prozess übergeben

Wir möchten die finale Ausgabedatei spezifisch benennen können, um nachfolgende Batches von Begrüßungen zu verarbeiten, ohne die finalen Ergebnisse zu überschreiben.

Zu diesem Zweck nehmen wir folgende Verfeinerungen am Workflow vor:

- Den Sammelprozess so ändern, dass er einen benutzerdefinierten Namen für die Ausgabedatei akzeptiert (`batch_name`)
- Einen Kommandozeilenparameter zum Workflow hinzufügen (`--batch`) und ihn an den Sammelprozess übergeben

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-collect-batch.svg"
</figure>

### 3.1. Den Sammelprozess ändern

Wir müssen die zusätzliche Eingabe deklarieren und sie in den Ausgabedateinamen integrieren.

#### 3.1.1. Die zusätzliche Eingabe deklarieren

Gute Nachrichten: Wir können so viele Eingabevariablen deklarieren, wie wir in der Prozessdefinition möchten.
Nennen wir diese `batch_name`.

Nimm im Prozessblock folgende Codeänderung vor:

=== "Danach"

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

Du kannst deine Prozesse so einrichten, dass sie so viele Eingaben erwarten, wie du möchtest.
Im Moment sind diese alle als erforderliche Eingaben eingerichtet; du _musst_ einen Wert angeben, damit der Workflow funktioniert.

Du wirst später in deiner Nextflow-Reise lernen, wie du erforderliche vs. optionale Eingaben verwaltest.

#### 3.1.2. Die `batch_name`-Variable im Ausgabedateinamen verwenden

Wir können die Variable auf die gleiche Weise in den Ausgabedateinamen einfügen, wie wir zuvor dynamische Dateinamen zusammengesetzt haben.

Nimm im Prozessblock folgende Codeänderung vor:

=== "Danach"

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

Dies richtet den Prozess so ein, dass er den `batch_name`-Wert verwendet, um einen spezifischen Dateinamen für die finale Ausgabe des Workflows zu generieren.

### 3.2. Einen `batch`-Kommandozeilenparameter hinzufügen

Jetzt brauchen wir eine Möglichkeit, den Wert für `batch_name` bereitzustellen und ihn an den Prozessaufruf weiterzugeben.

#### 3.2.1. `params` verwenden, um den Parameter einzurichten

Du weißt bereits, wie du das `params`-System verwendest, um CLI-Parameter zu deklarieren.
Lass uns das verwenden, um einen `batch`-Parameter zu deklarieren (mit einem Standardwert, weil wir faul sind).

Nimm im Abschnitt Pipeline-Parameter folgende Codeänderungen vor:

=== "Danach"

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

#### 3.2.2. Den `batch`-Parameter an den Prozess übergeben

Um den Wert des Parameters an den Prozess zu übergeben, müssen wir ihn im Prozessaufruf hinzufügen.

Nimm im Workflow-Block folgende Codeänderung vor:

=== "Danach"

    ```groovy title="hello-workflow.nf" linenums="74" hl_lines="2"
        // Alle Begrüßungen in einer Datei sammeln
        collectGreetings(convertToUpper.out.collect(), params.batch)
    ```

=== "Vorher"

    ```groovy title="hello-workflow.nf" linenums="74" hl_lines="2"
        // Alle Begrüßungen in einer Datei sammeln
        collectGreetings(convertToUpper.out.collect())
    ```

Du siehst, dass du, um mehrere Eingaben an einen Prozess zu übergeben, sie einfach in den Aufrufklammern auflistest, getrennt durch Kommas.

!!! warning

    Du MUSST die Eingaben an den Prozess in GENAU DERSELBEN REIHENFOLGE angeben, wie sie im Eingabedefinitionsblock des Prozesses aufgelistet sind.

### 3.3. Den Workflow ausführen

Lass uns versuchen, dies mit einem Batch-Namen auf der Kommandozeile auszuführen.

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

Es läuft erfolgreich und erzeugt die gewünschte Ausgabe:

??? abstract "Dateiinhalt"

    ```console title="results/COLLECTED-trio-output.txt"
    HELLO
    BONJOUR
    HOLà
    ```

Solange wir den Parameter entsprechend angeben, werden nachfolgende Läufe auf anderen Batches von Eingaben frühere Ergebnisse nicht überschreiben.

### Fazit

Du weißt, wie du mehr als eine Eingabe an einen Prozess übergibst.

### Wie geht es weiter?

Lerne, wie du mehrere Ausgaben ausgibst und sie bequem verarbeitest.

---

## 4. Eine Ausgabe zum Sammelschritt hinzufügen

Bisher haben wir Prozesse verwendet, die jeweils nur eine Ausgabe erzeugt haben.
Wir konnten auf ihre jeweiligen Ausgaben sehr bequem mit der `<process>.out`-Syntax zugreifen, die wir sowohl im Kontext der Übergabe einer Ausgabe an den nächsten Prozess (z.B. `convertToUpper(sayHello.out)`) als auch im Kontext des `publish:`-Abschnitts (z.B. `first_output = sayHello.out`) verwendet haben.

Was passiert, wenn ein Prozess mehr als eine erzeugt?
Wie gehen wir mit den mehreren Ausgaben um?
Können wir eine bestimmte Ausgabe auswählen und verwenden?

Alles ausgezeichnete Fragen, und die kurze Antwort ist: Ja, das können wir!

Mehrere Ausgaben werden in separate Kanäle verpackt.
Wir können entweder wählen, diesen Ausgabekanälen Namen zu geben, was es einfach macht, später individuell auf sie zu verweisen, oder wir können über ihren Index auf sie verweisen.

Zu Demonstrationszwecken nehmen wir an, wir möchten die Anzahl der Begrüßungen zählen, die für einen bestimmten Batch von Eingaben gesammelt werden, und dies in einer Datei berichten.

### 4.1. Den Prozess ändern, um die Anzahl der Begrüßungen zu zählen und auszugeben

Dies erfordert zwei wesentliche Änderungen an der Prozessdefinition: Wir brauchen eine Möglichkeit, die Begrüßungen zu zählen und eine Berichtsdatei zu schreiben, dann müssen wir diese Berichtsdatei zum `output`-Block des Prozesses hinzufügen.

#### 4.1.1. Die Anzahl der gesammelten Begrüßungen zählen

Praktischerweise erlaubt uns Nextflow, beliebigen Code im `script:`-Block der Prozessdefinition hinzuzufügen, was sehr praktisch für solche Dinge ist.

Das bedeutet, wir können Nextflows eingebaute `size()`-Funktion verwenden, um die Anzahl der Dateien im `input_files`-Array zu erhalten, und das Ergebnis mit einem `echo`-Befehl in eine Datei schreiben.

Nimm im `collectGreetings`-Prozessblock folgende Codeänderungen vor:

=== "Danach"

    ```groovy title="hello-workflow.nf" linenums="55" hl_lines="2 5"
        script:
        count_greetings = input_files.size()
        """
        cat ${input_files} > 'COLLECTED-${batch_name}-output.txt'
        echo 'Es gab ${count_greetings} Begrüßungen in diesem Batch.' > '${batch_name}-report.txt'
        """
    ```

=== "Vorher"

    ```groovy title="hello-workflow.nf" linenums="55"
        script:
        """
        cat ${input_files} > 'COLLECTED-${batch_name}-output.txt'
        """
    ```

Die Variable `count_greetings` wird zur Laufzeit berechnet.

#### 4.1.2. Die Berichtsdatei ausgeben und Ausgaben benennen

Im Prinzip müssen wir nur die Berichtsdatei zum `output:`-Block hinzufügen.

Während wir dabei sind, fügen wir jedoch auch einige `emit:`-Tags zu unseren Ausgabedeklarationen hinzu. Diese ermöglichen es uns, die Ausgaben nach Namen auszuwählen, anstatt positionelle Indizes verwenden zu müssen.

Nimm im Prozessblock folgende Codeänderung vor:

=== "Danach"

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
Aber wie das Sprichwort sagt: Warum nicht beides?

!!! tip

    Wenn du die Ausgaben eines Prozesses nicht mit `emit:` benennst, kannst du trotzdem individuell auf sie zugreifen, indem du ihren jeweiligen (nullbasierten) Index verwendest.
    Zum Beispiel würdest du `<process>.out[0]` verwenden, um die erste Ausgabe zu erhalten, `<process>.out[1]`, um die zweite Ausgabe zu erhalten, und so weiter.

    Wir bevorzugen es, Ausgaben zu benennen, weil es sonst zu einfach ist, versehentlich den falschen Index zu verwenden, besonders wenn der Prozess viele Ausgaben erzeugt.

### 4.2. Die Workflow-Ausgaben aktualisieren

Jetzt, da wir zwei Ausgaben aus dem `collectGreetings`-Prozess haben, enthält die `collectGreetings.out`-Ausgabe zwei Kanäle:

- `collectGreetings.out.outfile` enthält die finale Ausgabedatei
- `collectGreetings.out.report` enthält die Berichtsdatei

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-collect-report.svg"
</figure>

Wir müssen die Workflow-Ausgaben entsprechend aktualisieren.

#### 4.2.1. Den `publish:`-Abschnitt aktualisieren

Nimm im `workflow`-Block folgende Codeänderung vor:

=== "Danach"

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

Wie du sehen kannst, ist das Verweisen auf spezifische Prozessausgaben jetzt trivial.
Wenn wir in Teil 5 (Container) einen weiteren Schritt zu unserer Pipeline hinzufügen, können wir einfach auf `collectGreetings.out.outfile` verweisen und es an den neuen Prozess übergeben (Spoiler: der neue Prozess heißt `cowpy`).

Aber für jetzt lass uns die Aktualisierung der Workflow-Level-Ausgaben abschließen.

#### 4.2.2. Den `output`-Block aktualisieren

Nimm im `output`-Block folgende Codeänderung vor:

=== "Danach"

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

Lass uns versuchen, dies mit dem aktuellen Batch von Begrüßungen auszuführen.

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

Wenn du im Verzeichnis `results/hello_workflow/` nachsiehst, findest du die neue Berichtsdatei `trio-report.txt`.
Öffne sie, um zu überprüfen, dass der Workflow die Anzahl der verarbeiteten Begrüßungen korrekt gemeldet hat.

??? abstract "Dateiinhalt"

    ```txt title="trio-report.txt"
    Es gab 3 Begrüßungen in diesem Batch.
    ```

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-collect-4-way.svg"
</figure>

Du kannst gerne weitere Begrüßungen zur CSV hinzufügen und testen, was passiert.

### Fazit

Du weißt, wie du einen Prozess mehrere benannte Ausgaben ausgeben lässt und wie du sie auf Workflow-Ebene angemessen verarbeitest.

Allgemeiner verstehst du die Schlüsselprinzipien, die beim Verbinden von Prozessen auf gängige Weise beteiligt sind.

### Wie geht es weiter?

Mache eine extra lange Pause, du hast sie dir verdient.

Wenn du bereit bist, gehe weiter zu [**Teil 4: Hello Modules**](./04_hello_modules.md), um zu lernen, wie du deinen Code für bessere Wartbarkeit und Code-Effizienz modularisierst.

---

## Quiz

<quiz>
Wie greifst du auf die Ausgabe eines Prozesses im Workflow-Block zu?
- [ ] `process.output`
- [ ] `output.processName`
- [x] `processName.out`
- [ ] `get(processName)`

Mehr erfahren: [1.4. Die Ausgabe des ersten Prozesses an den zweiten Prozess übergeben](#14-die-ausgabe-des-ersten-prozesses-an-den-zweiten-prozess-ubergeben)
</quiz>

<quiz>
Was bestimmt die Reihenfolge der Prozessausführung in Nextflow?
- [ ] Die Reihenfolge, in der Prozesse im Workflow-Block geschrieben sind
- [ ] Alphabetische Reihenfolge nach Prozessname
- [x] Datenabhängigkeiten zwischen Prozessen
- [ ] Zufällige Reihenfolge für parallele Ausführung

Mehr erfahren: [1.4. Die Ausgabe des ersten Prozesses an den zweiten Prozess übergeben](#14-die-ausgabe-des-ersten-prozesses-an-den-zweiten-prozess-ubergeben)
</quiz>

<quiz>
Welcher Operator sollte `???` ersetzen, um alle Ausgaben in einer einzigen Liste für den nachgelagerten Prozess zu sammeln?

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

Mehr erfahren: [2.4. Einen Operator verwenden, um die Begrüßungen in einer einzigen Eingabe zu sammeln](#24-einen-operator-verwenden-um-die-begrussungen-in-einer-einzigen-eingabe-zu-sammeln)
</quiz>

<quiz>
Wann solltest du den `collect()`-Operator verwenden?
- [ ] Wenn du Elemente parallel verarbeiten möchtest
- [ ] Wenn du Kanalinhalte filtern musst
- [x] Wenn ein nachgelagerter Prozess alle Elemente von einem vorgelagerten Prozess benötigt
- [ ] Wenn du Daten über mehrere Prozesse verteilen möchtest

Mehr erfahren: [2.4. Einen Operator verwenden, um die Begrüßungen in einer einzigen Eingabe zu sammeln](#24-einen-operator-verwenden-um-die-begrussungen-in-einer-einzigen-eingabe-zu-sammeln)
</quiz>

<quiz>
Wie greifst du auf eine benannte Ausgabe eines Prozesses zu?
- [ ] `processName.outputName`
- [ ] `processName.get(outputName)`
- [x] `processName.out.outputName`
- [ ] `output.processName.outputName`

Mehr erfahren: [4.1.2. Die Berichtsdatei ausgeben und Ausgaben benennen](#412-die-berichtsdatei-ausgeben-und-ausgaben-benennen)
</quiz>

<quiz>
Was ist die korrekte Syntax zum Benennen einer Ausgabe in einem Prozess?
- [ ] `name: outputName`
- [ ] `output: outputName`
- [x] `emit: outputName`
- [ ] `label: outputName`

Mehr erfahren: [4.1.2. Die Berichtsdatei ausgeben und Ausgaben benennen](#412-die-berichtsdatei-ausgeben-und-ausgaben-benennen)
</quiz>

<quiz>
Was muss wahr sein, wenn du mehrere Eingaben an einen Prozess übergibst?
- [ ] Alle Eingaben müssen vom gleichen Typ sein
- [ ] Eingaben müssen in alphabetischer Reihenfolge angegeben werden
- [x] Die Reihenfolge der Eingaben muss mit der im Eingabeblock definierten Reihenfolge übereinstimmen
- [ ] Es können nur zwei Eingaben gleichzeitig angegeben werden

Mehr erfahren: [3. Mehr als eine Eingabe an einen Prozess übergeben](#3-zusatzliche-parameter-an-einen-prozess-ubergeben)
</quiz>
