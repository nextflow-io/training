heit, eine Map zu verwenden, d.h. eine Datenstruktur, bei der jedes Element einen Satz von Schlüsseln und zugehörigen Werten hat, sodass du leicht auf jeden Schlüssel verweisen kannst, um den entsprechenden Wert zu erhalten.

In unserem Beispiel bedeutet das, von dieser Organisation zu wechseln:

```groovy
data = [patientA, 1, normal, R1]

println data[3]
```

Zu dieser hier:

```groovy
data = [id: patientA, replicate: 1, type: normal, readNum: 1]

println data.readNum
```

In Nextflow wird das eine [Map](https://nextflow.io/docs/latest/script.html#maps) genannt.

Lass uns unsere flache Liste jetzt in eine Map umwandeln.
Nimm folgende Änderungen am Workflow vor:

=== "Nachher"

    ```groovy title="main.nf" linenums="7" hl_lines="4-13"
        // Load files with channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
        ch_files.map { myFile ->
            def (patient, replicate, type, readNum) = myFile.simpleName.tokenize('_')
            [
              [
                id: patient,
                replicate: replicate.replace('rep', ''),
                type: type,
                readNum: readNum.replace('R', ''),
              ],
              myFile
            ]
        }
    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="7" hl_lines="4"
        // Load files with channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
        ch_files.map { myFile ->
            [ myFile.simpleName.tokenize('_'), myFile ]
        }
    ```

Die wichtigsten Änderungen hier sind:

- **Destrukturierende Zuweisung**: `def (patient, replicate, type, readNum) = ...` extrahiert die tokenisierten Werte in einer Zeile in benannte Variablen
- **Map-Literal-Syntax**: `[id: patient, replicate: ...]` erstellt eine Map, bei der jeder Schlüssel (wie `id`) mit einem Wert (wie `patient`) assoziiert ist
- **Verschachtelte Struktur**: Die äußere Liste `[..., myFile]` paart die Metadaten-Map mit dem ursprünglichen Dateiobjekt

Wir haben auch ein paar der Metadaten-Strings mit einer String-Ersetzungsmethode namens `replace()` vereinfacht, um unnötige Zeichen zu entfernen (_z.B._ `replicate.replace('rep', '')`, um nur die Zahl aus den Replikat-IDs zu behalten).

Lass uns den Workflow erneut ausführen:

```bash
nextflow run main.nf
```

??? success "Befehlsausgabe"

    ```console hl_lines="7-8"
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [infallible_swartz] DSL2 - revision: 7f4e68c0cb

    executor >  local (2)
    [1b/e7fb27] COUNT_LINES (1) [100%] 2 of 2 ✔
    [[id:patientA, replicate:1, type:normal, readNum:2], /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R2_001.fastq.gz]
    [[id:patientA, replicate:1, type:normal, readNum:1], /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz]
    Processing file: patientA_rep1_normal_R2_001.fastq.gz
    40

    Processing file: patientA_rep1_normal_R1_001.fastq.gz
    40
    ```

Jetzt sind die Metadaten sauber beschriftet (_z.B._ `[id:patientA, replicate:1, type:normal, readNum:2]`), sodass es viel einfacher ist zu erkennen, was was ist.

Es wird auch viel einfacher sein, Elemente der Metadaten im Workflow tatsächlich zu nutzen, und macht unseren Code lesbarer und wartbarer.

### Zusammenfassung

- Wir können Dateinamen in Nextflow mit der vollen Kraft einer Programmiersprache verarbeiten
- Wir können die Dateinamen als Strings behandeln, um relevante Informationen zu extrahieren
- Die Verwendung von Methoden wie `tokenize()` und `replace()` erlaubt es uns, Strings im Dateinamen zu manipulieren
- Die `.map()`-Operation transformiert Channel-Elemente unter Beibehaltung der Struktur
- Strukturierte Metadaten (Maps) machen Code lesbarer und wartbarer als positionsbasierte Listen

Als Nächstes schauen wir uns an, wie man mit gepaarten Datendateien umgeht.

---

## 5. Verarbeitung gepaarter Datendateien

Viele experimentelle Designs erzeugen gepaarte Datendateien, die von einer explizit gepaarten Handhabung profitieren.
Zum Beispiel werden in der Bioinformatik Sequenzierungsdaten oft in Form von gepaarten Reads erzeugt, d.h. Sequenz-Strings, die vom selben DNA-Fragment stammen (oft als 'Forward' und 'Reverse' bezeichnet, weil sie von entgegengesetzten Enden gelesen werden).

Das ist bei unseren Beispieldaten der Fall, wo R1 und R2 auf die beiden Read-Sets verweisen.

```console
data/patientA_rep1_normal_R1_001.fastq.gz
data/patientA_rep1_normal_R2_001.fastq.gz
```

Nextflow bietet eine spezialisierte Channel Factory für die Arbeit mit solchen gepaarten Dateien namens `channel.fromFilePairs()`, die Dateien basierend auf einem gemeinsamen Benennungsmuster automatisch gruppiert. Das ermöglicht es dir, die gepaarten Dateien enger mit weniger Aufwand zu assoziieren.

Wir werden unseren Workflow modifizieren, um dies zu nutzen.
Es wird zwei Schritte dauern:

1. Die Channel Factory auf `channel.fromFilePairs()` umstellen
2. Die Metadaten extrahieren und abbilden

### 5.1. Umstellung der Channel Factory auf `channel.fromFilePairs()`

Um `channel.fromFilePairs` zu verwenden, müssen wir das Muster angeben, das Nextflow verwenden soll, um die beiden Mitglieder eines Paares zu identifizieren.

Zurück zu unseren Beispieldaten können wir das Benennungsmuster wie folgt formalisieren:

```console
data/patientA_rep1_normal_R{1,2}_001.fastq.gz
```

Dies ist ähnlich wie das Glob-Muster, das wir früher verwendet haben, außer dass hier die Substrings (entweder `1` oder `2`, die direkt nach dem R kommen) speziell aufgezählt werden, die die beiden Mitglieder des Paares identifizieren.

Lass uns den Workflow `main.nf` entsprechend aktualisieren:

=== "Nachher"

    ```groovy title="main.nf" linenums="7" hl_lines="1-2"
        // Load files with channel.fromFilePairs
        ch_files = channel.fromFilePairs('data/patientA_rep1_normal_R{1,2}_001.fastq.gz')
        /* Comment out the mapping for now, we'll come back to it!
        ch_files.map { myFile ->
            def (sample, replicate, type, readNum) = myFile.simpleName.tokenize('_')
            [
                [
                    id: sample,
                    replicate: replicate.replace('rep', ''),
                    type: type,
                    readNum: readNum,
                ],
                myFile
            ]
        }
        */
        .view()
    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="7" hl_lines="1-2"
        // Load files with channel.fromFilePairs
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
        ch_files.map { myFile ->
            def (sample, replicate, type, readNum) = myFile.simpleName.tokenize('_')
            [
                [
                    id: sample,
                    replicate: replicate.replace('rep', ''),
                    type: type,
                    readNum: readNum,
                ],
                myFile
            ]
        }
        .view()
    ```

Wir haben die Channel Factory umgestellt und das Datei-Matching-Muster angepasst, und während wir dabei waren, haben wir die Map-Operation auskommentiert.
Wir werden sie später mit ein paar Modifikationen wieder hinzufügen.

Führe den Workflow aus, um ihn zu testen:

```bash
nextflow run main.nf
```

??? failure "Befehlsausgabe"

    ```console hl_lines="7-8"
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [angry_koch] DSL2 - revision: 44fdf66105

    [-        ] COUNT_LINES -
    [-        ] COUNT_LINES -
    [patientA_rep1_normal_R, [/workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz, /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R2_001.fastq.gz]]
    ERROR ~ Error executing process > 'COUNT_LINES (1)'

    Caused by:
      Not a valid path value: 'patientA_rep1_normal_R'



    Tip: when you have fixed the problem you can continue the execution adding the option `-resume` to the run command line

    -- Check '.nextflow.log' file for details
    ```

Oh-oh, dieses Mal ist die Ausführung fehlgeschlagen!

Der relevante Teil der Fehlermeldung ist hier:

```console
Not a valid path value: 'patientA_rep1_normal_R'
```

Das liegt daran, dass wir die Channel Factory geändert haben.
Bis jetzt enthielt der ursprüngliche Input-Channel nur die Dateipfade.
Alle Metadaten-Manipulationen, die wir vorgenommen haben, haben die Channel-Inhalte eigentlich nicht beeinflusst.

Jetzt, da wir die `.fromFilePairs`-Channel-Factory verwenden, sind die Inhalte des resultierenden Channels anders.
Wir sehen nur ein Channel-Element, bestehend aus einem Tupel mit zwei Elementen: dem Teil des `simpleName`, der von beiden Dateien geteilt wird und als Identifikator dient, und einem Tupel mit den beiden Dateiobjekten, im Format `id, [ file1, file2 ]`.

Das ist großartig, weil Nextflow die harte Arbeit geleistet hat, den Patientennamen zu extrahieren, indem es das gemeinsame Präfix untersucht und es als Patienten-Identifikator verwendet hat.

Es bricht jedoch unseren aktuellen Workflow. Wenn wir `COUNT_LINES` immer noch auf die gleiche Weise ausführen wollten, ohne den Prozess zu ändern, müssten wir eine Mapping-Operation anwenden, um die Dateipfade zu extrahieren.
Aber wir werden das nicht tun, weil unser letztendliches Ziel darin besteht, einen anderen Prozess, `ANALYZE_READS`, zu verwenden, der Dateipaare angemessen verarbeitet.

Also kommentieren wir einfach den Aufruf von `COUNT_LINES` aus (oder löschen ihn) und fahren fort.

=== "Nachher"

    ```groovy title="main.nf" linenums="26" hl_lines="2"
        // Count the lines in the file
        // COUNT_LINES(ch_files)
    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="26" hl_lines="2"
        // Count the lines in the file
        COUNT_LINES(ch_files)
    ```

Du kannst auch die `COUNT_LINES`-Include-Anweisung auskommentieren oder löschen, aber das hat keine funktionale Auswirkung.

Lass uns den Workflow jetzt erneut ausführen:

```bash
nextflow run main.nf
```

??? success "Befehlsausgabe"

    ```console hl_lines="5"
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [fabulous_davinci] DSL2 - revision: 22b53268dc

    [patientA_rep1_normal_R, [/workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz, /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R2_001.fastq.gz]]
    ```

Juhu, dieses Mal gelingt der Workflow!

Wir müssen jedoch noch die restlichen Metadaten aus dem `id`-Feld extrahieren.

### 5.2. Extraktion und Organisation von Metadaten aus Dateipaaren

Unsere `map`-Operation von vorher funktioniert nicht, weil sie nicht zur Datenstruktur passt, aber wir können sie anpassen, damit sie funktioniert.

Wir haben bereits Zugriff auf den eigentlichen Patienten-Identifikator im String, den `fromFilePairs()` als Identifikator verwendet hat, sodass wir diesen verwenden können, um die Metadaten zu extrahieren, ohne den `simpleName` vom Path-Objekt zu holen, wie wir es vorher getan haben.

Kommentiere die Map-Operation im Workflow aus und nimm folgende Änderungen vor:

=== "Nachher"

    ```groovy title="main.nf" linenums="7" hl_lines="3-4 9 11 13"
        // Load files with channel.fromFilePairs
        ch_files = channel.fromFilePairs('data/patientA_rep1_normal_R{1,2}_001.fastq.gz')
        ch_files.map { id, files ->
            def (sample, replicate, type) = id.tokenize('_')
            [
                [
                    id: sample,
                    replicate: replicate.replace('rep', ''),
                    type: type
                ],
                files
            ]
        }
        .view()
    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="7" hl_lines="3-5 11 13"
        // Load files with channel.fromFilePairs
        ch_files = channel.fromFilePairs('data/patientA_rep1_normal_R{1,2}_001.fastq.gz')
        /* Comment out the mapping for now, we'll come back to it!
        ch_files.map { myFile ->
            def (sample, replicate, type, readNum) = myFile.simpleName.tokenize('_')
            [
                [
                    id: sample,
                    replicate: replicate.replace('rep', ''),
                    type: type,
                    readNum: readNum,
                ],
                myFile
            ]
        }
        */
        .view()
    ```

Dieses Mal beginnt die Map mit `id, files` anstelle von nur `myFile`, und `tokenize()` wird auf `id` anstelle von `myFile.simpleName` angewendet.

Beachte auch, dass wir `readNum` aus der `tokenize()`-Zeile entfernt haben; alle Substrings, die wir nicht speziell benennen (von links beginnend), werden stillschweigend verworfen.
Wir können das tun, weil die gepaarten Dateien jetzt eng assoziiert sind, sodass wir `readNum` nicht mehr in der Metadaten-Map benötigen.

Lass uns den Workflow ausführen:

```bash
nextflow run main.nf
```

??? success "Befehlsausgabe"

    ```console

    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [prickly_stonebraker] DSL2 - revision: f62ab10a3f

    [[id:patientA, replicate:1, type:normal], [/workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz, /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R2_001.fastq.gz]]
    ```

Und da ist es: Wir haben die Metadaten-Map (`[id:patientA, replicate:1, type:normal]`) an der ersten Position des Ausgabe-Tupels, gefolgt vom Tupel der gepaarten Dateien, wie beabsichtigt.

Natürlich wird dies nur dieses spezifische Dateipaar erfassen und verarbeiten.
Wenn du mit der Verarbeitung mehrerer Paare experimentieren möchtest, kannst du versuchen, Wildcards in das Eingabemuster einzufügen und sehen, was passiert.
Versuche zum Beispiel `data/patientA_rep1_*_R{1,2}_001.fastq.gz` zu verwenden.

### Zusammenfassung

- [`channel.fromFilePairs()` findet und paart automatisch zusammengehörige Dateien](https://www.nextflow.io/docs/latest/reference/channel.html#fromfilepairs)
- Dies vereinfacht die Handhabung von Paired-End-Reads in deiner Pipeline
- Gepaarte Dateien können als `[id, [file1, file2]]`-Tupel gruppiert werden
- Die Metadatenextraktion kann von der ID des gepaarten Dateipaares statt von einzelnen Dateien erfolgen

---

## 6. Verwendung von Dateioperationen in Prozessen

Lass uns jetzt all dies in einem einfachen Prozess zusammenführen, um zu verstärken, wie Dateioperationen innerhalb eines Nextflow-Prozesses verwendet werden.

Wir stellen dir ein vorgefertigtes Prozessmodul namens `ANALYZE_READS` zur Verfügung, das ein Tupel aus Metadaten und einem Paar von Eingabedateien nimmt und sie analysiert.
Wir könnten uns vorstellen, dass dies Sequenz-Alignment, Variant Calling oder jeden anderen Schritt durchführt, der für diesen Datentyp sinnvoll ist.

Lass uns anfangen.

### 6.1. Importiere den Prozess und untersuche den Code

Um diesen Prozess im Workflow zu verwenden, müssen wir nur eine Modul-Include-Anweisung vor dem workflow-Block hinzufügen.

Nimm folgende Änderung am Workflow vor:

=== "Nachher"

    ```groovy title="main.nf" linenums="1" hl_lines="3"
    #!/usr/bin/env nextflow

    include { ANALYZE_READS } from './modules/analyze_reads.nf'

    workflow {
    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="1"
    #!/usr/bin/env nextflow

    workflow {
    ```

Du kannst die Moduldatei öffnen, um ihren Code zu untersuchen:

```groovy title="modules/analyze_reads.nf - Prozessbeispiel" linenums="1"
#!/usr/bin/env nextflow

process ANALYZE_READS {
    tag { meta.id }

    publishDir { "results/${meta.id}" }, mode: 'copy'

    input:
    tuple val(meta), path(files)

    output:
    tuple val(meta.id), path("${meta.id}_stats.txt")

    script:
    """
    echo "Sample metadata: ${meta.id}" > ${meta.id}_stats.txt
    echo "Replicate: ${meta.replicate}" >> ${meta.id}_stats.txt
    echo "Type: ${meta.type}" >> ${meta.id}_stats.txt
    echo "Read 1: ${files[0]}" >> ${meta.id}_stats.txt
    echo "Read 2: ${files[1]}" >> ${meta.id}_stats.txt
    echo "Read 1 size: \$(gunzip -dc ${files[0]} | wc -l | awk '{print \$1/4}') reads" >> ${meta.id}_stats.txt
    echo "Read 2 size: \$(gunzip -dc ${files[1]} | wc -l | awk '{print \$1/4}') reads" >> ${meta.id}_stats.txt
    """
}
```

!!! note

    Die `tag`- und `publishDir`-Direktiven verwenden Closure-Syntax (`{ ... }`) anstelle von String-Interpolation (`"${...}"`).
    Das liegt daran, dass diese Direktiven auf Eingabevariablen (`meta`) verweisen, die erst zur Laufzeit verfügbar sind.
    Die Closure-Syntax verschiebt die Auswertung bis zur tatsächlichen Ausführung des Prozesses.

!!! note

    Wir nennen unsere Metadaten-Map per Konvention `meta`.
    Für einen tieferen Einblick in Meta-Maps siehe die Side Quest [Metadata and meta maps](./metadata.md).

### 6.2. Aufruf des Prozesses im Workflow

Jetzt, da der Prozess für den Workflow verfügbar ist, können wir einen Aufruf des `ANALYZE_READS`-Prozesses hinzufügen, um ihn auszuführen.

Um ihn auf unseren Beispieldaten auszuführen, müssen wir zwei Dinge tun:

1. Dem remapped Channel einen Namen geben
2. Einen Aufruf des Prozesses hinzufügen

#### 6.2.1. Benennung des remapped Input-Channels

Wir haben die Mapping-Manipulationen bisher direkt auf den Input-Channel angewendet.
Um die remapped Inhalte an den `ANALYZE_READS`-Prozess zu übergeben (und dies auf eine Weise zu tun, die klar und leicht zu lesen ist), wollen wir einen neuen Channel namens `ch_samples` erstellen.

Wir können das mit dem [`set`](https://www.nextflow.io/docs/latest/reference/operator.html#set)-Operator tun.

Ersetze im Haupt-Workflow den `.view()`-Operator durch `.set { ch_samples }` und füge eine Zeile hinzu, die testet, dass wir auf den Channel mit Namen verweisen können.

=== "Nachher"

    ```groovy title="main.nf" linenums="7" hl_lines="14 16-17"
        // Load files with channel.fromFilePairs
        ch_files = channel.fromFilePairs('data/patientA_rep1_normal_R{1,2}_001.fastq.gz')
        ch_files.map { id,  files ->
           def (sample, replicate, type, readNum) = id.tokenize('_')
           [
               [
                   id: sample,
                   replicate: replicate.replace('rep', ''),
                   type: type
               ],
               files
           ]
        }
            .set { ch_samples }

        // Temporary: peek into ch_samples
        ch_samples.view()
    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="7" hl_lines="14"
        // Load files with channel.fromFilePairs
        ch_files = channel.fromFilePairs('data/patientA_rep1_normal_R{1,2}_001.fastq.gz')
        ch_files.map { id,  files ->
           def (sample, replicate, type, readNum) = id.tokenize('_')
           [
               [
                   id: sample,
                   replicate: replicate.replace('rep', ''),
                   type: type
               ],
               files
           ]
        }
        .view()
    }
    ```

Lass uns dies ausführen:

```bash
nextflow run main.nf
```

??? success "Befehlsausgabe"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [goofy_kirch] DSL2 - revision: 3313283e42

    [[id:patientA, replicate:1, type:normal], [/workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz, /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R2_001.fastq.gz]]
    ```

Dies bestätigt, dass wir jetzt auf den Channel mit Namen verweisen können.

#### 6.2.2. Aufruf des Prozesses auf den Daten

Lass uns jetzt tatsächlich den `ANALYZE_READS`-Prozess auf dem `ch_samples`-Channel aufrufen.

Nimm im Haupt-Workflow folgende Codeänderungen vor:

=== "Nachher"

    ```groovy title="main.nf" linenums="23"
        // Run the analysis
        ANALYZE_READS(ch_samples)
    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="23"
        // Temporary: peek into ch_samples
        ch_samples.view()
    ```

Lass uns dies ausführen:

```bash
nextflow run main.nf
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `./main.nf` [shrivelled_cori] DSL2 - revision: b546a31769

    executor >  local (1)
    [b5/110360] process > ANALYZE_READS (patientA) [100%] 1 of 1 ✔
    ```

Dieser Prozess ist so eingerichtet, dass er seine Ausgaben in ein `results`-Verzeichnis publiziert, also schaue dort nach.

??? abstract "Verzeichnis- und Dateiinhalte"

    ```console
    results
    └── patientA
        └── patientA_stats.txt
    ```

    ```txt title="patientA_stats.txt"
    Sample metadata: patientA
    Replicate: 1
    Type: normal
    Read 1: patientA_rep1_normal_R1_001.fastq.gz
    Read 2: patientA_rep1_normal_R2_001.fastq.gz
    Read 1 size: 10 reads
    Read 2 size: 10 reads
    ```

Der Prozess hat unsere Eingaben genommen und eine neue Datei erstellt, die die Patienten-Metadaten enthält, wie vorgesehen.
Hervorragend!

### 6.3. Viele weitere Patienten einbeziehen

Natürlich verarbeitet dies nur ein einzelnes Dateipaar für einen einzelnen Patienten, was nicht gerade der hohe Durchsatz ist, den du mit Nextflow erhoffst.
Du wirst wahrscheinlich viel mehr Daten auf einmal verarbeiten wollen.

Denke daran, dass `channel.fromPath()` einen _Glob_ als Eingabe akzeptiert, was bedeutet, dass es beliebig viele Dateien akzeptieren kann, die dem Muster entsprechen.
Wenn wir also alle Patienten einbeziehen wollen, können wir einfach den Input-String so modifizieren, dass er mehr Patienten einschließt, wie früher beiläufig erwähnt.

Lass uns so tun, als wollten wir so großzügig wie möglich sein.
Nimm folgende Änderungen am Workflow vor:

=== "Nachher"

    ```groovy title="main.nf" linenums="7" hl_lines="2"
        // Load files with channel.fromFilePairs
        ch_files = channel.fromFilePairs('data/*_R{1,2}_001.fastq.gz')
    ```

=== "Vorher"

    ```groovy title="main.nf" linenums="7" hl_lines="2"
        // Load files with channel.fromFilePairs
        ch_files = channel.fromFilePairs('data/patientA_rep1_normal_R{1,2}_001.fastq.gz')
    ```

Führe die Pipeline erneut aus:

```bash
nextflow run main.nf
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `./main.nf` [big_stonebraker] DSL2 - revision: f7f9b8a76c

    executor >  local (8)
    [d5/441891] process > ANALYZE_READS (patientC) [100%] 8 of 8 ✔
    ```

Das results-Verzeichnis sollte jetzt Ergebnisse für alle verfügbaren Daten enthalten.

??? abstract "Verzeichnisinhalte"

    ```console
    results
    ├── patientA
    │   └── patientA_stats.txt
    ├── patientB
    │   └── patientB_stats.txt
    └── patientC
        └── patientC_stats.txt
    ```

Erfolg! Wir haben alle Patienten auf einmal analysiert! Richtig?

Vielleicht nicht.
Wenn du genauer hinschaust, haben wir ein Problem: Wir haben zwei Replikate für patientA, aber nur eine Ausgabedatei!
Wir überschreiben die Ausgabedatei jedes Mal.

### 6.4. Die publizierten Dateien eindeutig machen

Da wir Zugriff auf die Patienten-Metadaten haben, können wir sie verwenden, um die publizierten Dateien eindeutig zu machen, indem wir differenzierende Metadaten entweder in der Verzeichnisstruktur oder in den Dateinamen selbst einbeziehen.

Nimm folgende Änderung am Workflow vor:

=== "Nachher"

    ```groovy title="modules/analyze_reads.nf" linenums="6"
        publishDir { "results/${meta.type}/${meta.id}/${meta.replicate}" }, mode: 'copy'
    ```

=== "Vorher"

    ```groovy title="modules/analyze_reads.nf" linenums="6"
        publishDir { "results/${meta.id}" }, mode: 'copy'
    ```

Hier zeigen wir die Option, zusätzliche Verzeichnisebenen zu verwenden, um Probentypen und Replikate zu berücksichtigen, aber du könntest auch experimentieren, dies auf Dateinamenebene zu tun.

Führe nun die Pipeline ein letztes Mal aus, aber stelle sicher, dass du zuerst das results-Verzeichnis entfernst, um dir einen sauberen Arbeitsbereich zu verschaffen:

```bash
rm -r results
nextflow run main.nf
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `./main.nf` [insane_swartz] DSL2 - revision: fff18abe6d

    executor >  local (8)
    [e3/449081] process > ANALYZE_READS (patientC) [100%] 8 of 8 ✔
    ```

Überprüfe jetzt das results-Verzeichnis:

??? abstract "Verzeichnisinhalte"

    ```console
    results/
    ├── normal
    │   ├── patientA
    │   │   ├── 1
    │   │   │   └── patientA_stats.txt
    │   │   └── 2
    │   │       └── patientA_stats.txt
    │   ├── patientB
    │   │   └── 1
    │   │       └── patientB_stats.txt
    │   └── patientC
    │       └── 1
    │           └── patientC_stats.txt
    └── tumor
        ├── patientA
        │   ├── 1
        │   │   └── patientA_stats.txt
        │   └── 2
        │       └── patientA_stats.txt
        ├── patientB
        │   └── 1
        │       └── patientB_stats.txt
        └── patientC
            └── 1
                └── patientC_stats.txt
    ```

Und da ist es, all unsere Metadaten, sauber organisiert. Das ist Erfolg!

Es gibt noch viel mehr, was du tun kannst, sobald du deine Metadaten in eine Map wie diese geladen hast:

1. Organisierte Ausgabeverzeichnisse basierend auf Patientenattributen erstellen
2. Entscheidungen in Prozessen basierend auf Patienteneigenschaften treffen
3. Daten basierend auf Metadatenwerten aufteilen, verbinden und rekombinieren

Dieses Muster, Metadaten explizit und an die Daten gebunden zu halten (anstatt in Dateinamen kodiert), ist eine zentrale Best Practice in Nextflow, die den Aufbau robuster, wartbarer Analyse-Workflows ermöglicht.
Mehr darüber kannst du in der Side Quest [Metadata and meta maps](./metadata.md) erfahren.

### Zusammenfassung

- Die `publishDir`-Direktive kann Ausgaben basierend auf Metadatenwerten organisieren
- Metadaten in Tupeln ermöglichen eine strukturierte Organisation von Ergebnissen
- Dieser Ansatz schafft wartbare Workflows mit klarer Datenherkunft
- Prozesse können Tupel aus Metadaten und Dateien als Eingabe nehmen
- Die `tag`-Direktive bietet Prozessidentifikation in Ausführungslogs
- Die Workflow-Struktur trennt Channel-Erstellung von Prozessausführung

---

## Zusammenfassung

In dieser Side Quest hast du gelernt, wie man mit Dateien in Nextflow arbeitet, von grundlegenden Operationen bis hin zu fortgeschritteneren Techniken für die Handhabung von Dateisammlungen.

Die Anwendung dieser Techniken in deiner eigenen Arbeit wird es dir ermöglichen, effizientere und wartbarere Workflows zu erstellen, besonders beim Arbeiten mit großen Mengen von Dateien mit komplexen Benennungskonventionen.

### Wichtige Muster

1.  **Grundlegende Dateioperationen:** Wir erstellten Path-Objekte mit `file()` und griffen auf Dateiattribute wie Name, Extension und übergeordnetes Verzeichnis zu, wobei wir den Unterschied zwischen Strings und Path-Objekten kennenlernten.

    - Ein Path-Objekt erstellen mit `file()`

    ```groovy
    myFile = file('path/to/file.txt')
    ```

    - Dateiattribute abrufen

    ```groovy
    println myFile.name       // file.txt
    println myFile.baseName   // file
    println myFile.extension  // txt
    println myFile.parent     // path/to
    ```

2.  **Verwendung entfernter Dateien**: Wir lernten, wie man transparent zwischen lokalen und entfernten Dateien mittels URIs wechselt und demonstrierten Nextflows Fähigkeit, Dateien aus verschiedenen Quellen zu verarbeiten, ohne die Workflow-Logik zu ändern.

    - Lokale Datei

    ```groovy
    myFile = file('path/to/file.txt')
    ```

    - FTP

    ```groovy
    myFile = file('ftp://path/to/file.txt')
    ```

    - HTTPS

    ```groovy
    myFile = file('https://path/to/file.txt')
    ```

    - Amazon S3

    ```groovy
    myFile = file('s3://path/to/file.txt')
    ```

    - Azure Blob Storage

    ```groovy
    myFile = file('az://path/to/file.txt')
    ```

    - Google Cloud Storage

    ```groovy
    myFile = file('gs://path/to/file.txt')
    ```

3.  **Laden von Dateien mit der `fromPath()`-Channel-Factory:** Wir erstellten Channels aus Dateimustern mit `channel.fromPath()` und zeigten ihre Dateiattribute an, einschließlich Objekttypen.

    - Einen Channel aus einem Dateimuster erstellen

    ```groovy
     ch_files = channel.fromPath('data/*.fastq.gz')
    ```

    - Dateiattribute abrufen

    ```groovy
     ch_files.view { myFile ->
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    }
    ```

4.  **Extraktion von Patienten-Metadaten aus Dateinamen:** Wir verwendeten `tokenize()` und `replace()`, um Metadaten aus Dateinamen zu extrahieren und zu strukturieren und sie in organisierte Maps umzuwandeln.

    ```groovy
    def name = file.name.tokenize('_')
    def patientId = name[0]
    def replicate = name[1].replace('rep', '')
    def type = name[2]
    def readNum = name[3].replace('R', '')
    ```

5.  **Vereinfachung mit channel.fromFilePairs:** Wir verwendeten `channel.fromFilePairs()`, um automatisch zusammengehörige Dateien zu paaren und Metadaten aus IDs gepaarter Dateien zu extrahieren.

    ```groovy
    ch_pairs = channel.fromFilePairs('data/*_R{1,2}_001.fastq.gz')
    ```

6.  **Verwendung von Dateioperationen in Prozessen:** Wir integrierten Dateioperationen in Nextflow-Prozesse mit korrekter Input-Verarbeitung und verwendeten `publishDir`, um Ausgaben basierend auf Metadaten zu organisieren.

    - Eine Meta-Map mit den Prozess-Inputs assoziieren

    ```groovy
    ch_files = channel.fromFilePairs('data/patientA_rep1_normal_R{1,2}_001.fastq.gz')
    ch_files.map { id,  files ->
        def (sample, replicate, type, readNum) = id.tokenize('_')
        [
            [
                id: sample,
                replicate: replicate.replace('rep', ''),
                type: type
            ],
             files
        ]
    }
        .set { ch_samples }

    ANALYZE_READS(ch_samples)
    ```

    - Ausgaben basierend auf Metadaten organisieren

    ```groovy
    publishDir { "results/${meta.type}/${meta.id}/${meta.replicate}" }, mode: 'copy'
    ```

### Zusätzliche Ressourcen

- [Nextflow-Dokumentation: Arbeiten mit Dateien](https://www.nextflow.io/docs/latest/working-with-files.html)
- [channel.fromPath](https://www.nextflow.io/docs/latest/reference/channel.html#frompath)
- [channel.fromFilePairs](https://www.nextflow.io/docs/latest/reference/channel.html#fromfilepairs)

---

## Wie geht es weiter?

Kehre zum [Menü der Side Quests](./index.md) zurück oder klicke auf die Schaltfläche unten rechts auf der Seite, um mit dem nächsten Thema in der Liste fortzufahren.
