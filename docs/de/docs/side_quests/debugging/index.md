# Workflows debuggen

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Debugging ist eine wichtige Fähigkeit, die dir stundenlange Frustration ersparen und dich zu einer effektiveren Nextflow-Entwickler\*in machen kann. Im Laufe deiner Karriere, besonders am Anfang, wirst du beim Entwickeln und Pflegen von Workflows auf Bugs stoßen. Systematische Debugging-Ansätze helfen dir, Probleme schnell zu erkennen und zu lösen.

### Lernziele

In dieser Side Quest erkunden wir **systematische Debugging-Techniken** für Nextflow-Workflows:

- **Syntaxfehler debuggen**: IDE-Funktionen und Nextflow-Fehlermeldungen effektiv nutzen
- **Kanal-Debugging**: Datenflussprobleme und Kanalstrukturprobleme diagnostizieren
- **Prozess-Debugging**: Ausführungsfehler und Ressourcenprobleme untersuchen
- **Eingebaute Debugging-Tools**: Nextflows Preview-Modus, Stub-Ausführung und work directories nutzen
- **Systematische Ansätze**: Eine Vier-Phasen-Methodik für effizientes Debugging

Am Ende hast du eine robuste Debugging-Methodik, die frustrierende Fehlermeldungen in klare Lösungswege verwandelt.

### Voraussetzungen

Bevor du diese Side Quest angehst, solltest du:

- Das Tutorial [Hello Nextflow](../hello_nextflow/README.md) oder einen gleichwertigen Einsteigerkurs abgeschlossen haben.
- Mit grundlegenden Nextflow-Konzepten und -Mechanismen (Prozesse, Kanäle, Operatoren) vertraut sein.

**Optional:** Wir empfehlen, zuerst die Side Quest [IDE Features for Nextflow Development](../dev_environment/) abzuschließen.
Diese behandelt umfassend IDE-Funktionen, die das Debugging unterstützen (Syntaxhervorhebung, Fehlererkennung usw.), die wir hier intensiv nutzen werden.

---

## 0. Erste Schritte

#### Trainings-Codespace öffnen

Falls noch nicht geschehen, öffne die Trainingsumgebung wie in der [Umgebung einrichten](../envsetup/index.md) beschrieben.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### In das Projektverzeichnis wechseln

Wechsle in das Verzeichnis, in dem sich die Dateien für dieses Tutorial befinden.

```bash
cd side-quests/debugging
```

Du kannst VSCode so einstellen, dass es sich auf dieses Verzeichnis fokussiert:

```bash
code .
```

#### Die Materialien ansehen

Du findest eine Reihe von Beispiel-Workflows mit verschiedenen Arten von Bugs, die wir zum Üben verwenden:

??? abstract "Verzeichnisinhalt"

    ```console
    .
    ├── bad_bash_var.nf
    ├── bad_channel_shape.nf
    ├── bad_channel_shape_viewed_debug.nf
    ├── bad_channel_shape_viewed.nf
    ├── bad_number_inputs.nf
    ├── badpractice_syntax.nf
    ├── bad_resources.nf
    ├── bad_syntax.nf
    ├── buggy_workflow.nf
    ├── data
    │   ├── sample_001.fastq.gz
    │   ├── sample_002.fastq.gz
    │   ├── sample_003.fastq.gz
    │   ├── sample_004.fastq.gz
    │   ├── sample_005.fastq.gz
    │   └── sample_data.csv
    ├── exhausted.nf
    ├── invalid_process.nf
    ├── missing_output.nf
    ├── missing_software.nf
    ├── missing_software_with_stub.nf
    ├── nextflow.config
    └── no_such_var.nf
    ```

Diese Dateien repräsentieren häufige Debugging-Szenarien, denen du in der realen Entwicklung begegnen wirst.

#### Die Aufgabe verstehen

Deine Aufgabe ist es, jeden Workflow auszuführen, den/die Fehler zu identifizieren und sie zu beheben.

Für jeden fehlerhaften Workflow:

1. **Führe den Workflow aus** und beobachte den Fehler
2. **Analysiere die Fehlermeldung**: Was teilt dir Nextflow mit?
3. **Finde das Problem** im Code anhand der gegebenen Hinweise
4. **Behebe den Bug** und überprüfe, ob deine Lösung funktioniert
5. **Setze die Datei zurück**, bevor du zum nächsten Abschnitt gehst (verwende `git checkout <filename>`)

Die Übungen beginnen mit einfachen Syntaxfehlern und werden zu subtileren Laufzeitproblemen.
Die Lösungen werden direkt im Text besprochen, aber versuche, jede Aufgabe selbst zu lösen, bevor du weiterliest.

#### Bereitschafts-Checkliste

Bereit zum Eintauchen?

- [ ] Ich verstehe das Ziel dieses Kurses und seine Voraussetzungen
- [ ] Mein Codespace läuft
- [ ] Ich habe mein Arbeitsverzeichnis entsprechend eingestellt
- [ ] Ich verstehe die Aufgabe

Wenn du alle Punkte abhaken kannst, kann es losgehen.

---

## 1. Syntaxfehler

Syntaxfehler sind die häufigste Art von Fehlern, die du beim Schreiben von Nextflow-Code begegnen wirst. Sie treten auf, wenn der Code nicht den erwarteten Syntaxregeln des Nextflow-DSL entspricht. Diese Fehler verhindern, dass dein Workflow überhaupt ausgeführt wird, daher ist es wichtig zu lernen, wie man sie schnell erkennt und behebt.

### 1.1. Fehlende geschweifte Klammern

Einer der häufigsten Syntaxfehler, und manchmal einer der komplexeren beim Debuggen, sind **fehlende oder nicht übereinstimmende Klammern**.

Fangen wir mit einem praktischen Beispiel an.

#### Pipeline ausführen

```bash
nextflow run bad_syntax.nf
```

??? failure "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_syntax.nf` [stupefied_bhabha] DSL2 - revision: ca6327fad2

    Error bad_syntax.nf:24:1: Unexpected input: '<EOF>'

    ERROR ~ Script compilation failed

     -- Check '.nextflow.log' file for details
    ```

**Wichtige Elemente von Syntaxfehlermeldungen:**

- **Datei und Position**: Zeigt, welche Datei und welche Zeile/Spalte den Fehler enthält (`bad_syntax.nf:24:1`)
- **Fehlerbeschreibung**: Erklärt, was der Parser gefunden hat, das er nicht erwartet hat (`Unexpected input: '<EOF>'`)
- **EOF-Indikator**: Die `<EOF>`-Meldung (End Of File) zeigt an, dass der Parser das Dateiende erreicht hat, während er noch mehr Inhalt erwartet hat – ein klassisches Zeichen für nicht geschlossene geschweifte Klammern

#### Den Code prüfen

Schauen wir uns `bad_syntax.nf` an, um zu verstehen, was den Fehler verursacht:

```groovy title="bad_syntax.nf" hl_lines="14" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    input:
    val sample_name

    output:
    path "${sample_name}_output.txt"

    script:
    """
    echo "Processing ${sample_name}" > ${sample_name}_output.txt
    """
// Fehlende schließende geschweifte Klammer für den Prozess

workflow {

    // Eingabekanal erstellen
    input_ch = channel.of('sample1', 'sample2', 'sample3')

    // Prozess mit dem Eingabekanal aufrufen
    PROCESS_FILES(input_ch)
}
```

Für dieses Beispiel haben wir einen Kommentar hinterlassen, der zeigt, wo der Fehler liegt. Die Nextflow-VSCode-Erweiterung sollte dir ebenfalls Hinweise geben, indem sie die nicht übereinstimmende geschweifte Klammer rot markiert und das vorzeitige Dateiende hervorhebt:

![Bad syntax](img/bad_syntax.png)

**Debugging-Strategie für Klammerfehler:**

1. VS Codes Klammerabgleich verwenden (Cursor neben eine Klammer setzen)
2. Das Problems-Panel auf klammerbezogene Meldungen prüfen
3. Sicherstellen, dass jede öffnende `{` eine entsprechende schließende `}` hat

#### Den Code korrigieren

Ersetze den Kommentar durch die fehlende schließende geschweifte Klammer:

=== "Danach"

    ```groovy title="bad_syntax.nf" hl_lines="14" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        """
        echo "Processing ${sample_name}" > ${sample_name}_output.txt
        """
    }  // Fehlende schließende geschweifte Klammer hinzufügen

    workflow {

        // Eingabekanal erstellen
        input_ch = channel.of('sample1', 'sample2', 'sample3')

        // Prozess mit dem Eingabekanal aufrufen
        PROCESS_FILES(input_ch)
    }
    ```

=== "Vorher"

    ```groovy title="bad_syntax.nf" hl_lines="14" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        """
        echo "Processing ${sample_name}" > ${sample_name}_output.txt
        """
    // Fehlende schließende geschweifte Klammer für den Prozess

    workflow {

        // Eingabekanal erstellen
        input_ch = channel.of('sample1', 'sample2', 'sample3')

        // Prozess mit dem Eingabekanal aufrufen
        PROCESS_FILES(input_ch)
    }
    ```

#### Pipeline ausführen

Führe den Workflow erneut aus, um zu bestätigen, dass er funktioniert:

```bash
nextflow run bad_syntax.nf
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_syntax.nf` [insane_faggin] DSL2 - revision: 961938ee2b

    executor >  local (3)
    [48/cd7f54] PROCESS_FILES (1) | 3 of 3 ✔
    ```

### 1.2. Falsche Prozess-Keywords oder Direktiven verwenden

Ein weiterer häufiger Syntaxfehler ist eine **ungültige Prozessdefinition**. Das kann passieren, wenn du vergisst, erforderliche Blöcke zu definieren, oder falsche Direktiven in der Prozessdefinition verwendest.

#### Pipeline ausführen

```bash
nextflow run invalid_process.nf
```

??? failure "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `invalid_process.nf` [nasty_jepsen] DSL2 - revision: da9758d614

    Error invalid_process.nf:3:1: Invalid process definition -- check for missing or out-of-order section labels
    │   3 | process PROCESS_FILES {
    │     | ^^^^^^^^^^^^^^^^^^^^^^^
    │   4 |     inputs:
    │   5 |     val sample_name
    │   6 |
    ╰   7 |     output:

    ERROR ~ Script compilation failed

     -- Check '.nextflow.log' file for details
    ```

#### Den Code prüfen

Der Fehler weist auf eine „ungültige Prozessdefinition" hin und zeigt den Kontext rund um das Problem. In den Zeilen 3–7 sehen wir `inputs:` in Zeile 4, was das Problem ist. Schauen wir uns `invalid_process.nf` an:

```groovy title="invalid_process.nf" hl_lines="4" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    inputs:  // FEHLER: Sollte 'input' sein, nicht 'inputs'
    val sample_name

    output:
    path "${sample_name}_output.txt"

    script:
    """
    echo "Processing ${sample_name}" > ${sample_name}_output.txt
    """
}

workflow {

    // Eingabekanal erstellen
    input_ch = channel.of('sample1', 'sample2', 'sample3')

    // Prozess mit dem Eingabekanal aufrufen
    PROCESS_FILES(input_ch)
}
```

In Zeile 4 des Fehlerkontexts sehen wir das Problem: Wir verwenden `inputs` statt der korrekten `input`-Direktive. Die Nextflow-VSCode-Erweiterung wird das ebenfalls markieren:

![Invalid process message](img/invalid_process_message.png)

#### Den Code korrigieren

Ersetze das falsche Keyword durch das richtige, indem du [die Dokumentation](https://www.nextflow.io/docs/latest/process.html#) zu Rate ziehst:

=== "Danach"

    ```groovy title="invalid_process.nf" hl_lines="4" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:  // Korrigiert: 'inputs' zu 'input' geändert
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        """
        echo "Processing ${sample_name}" > ${sample_name}_output.txt
        """
    }

    workflow {

        // Eingabekanal erstellen
        input_ch = channel.of('sample1', 'sample2', 'sample3')

        // Prozess mit dem Eingabekanal aufrufen
        PROCESS_FILES(input_ch)
    }
    ```

=== "Vorher"

    ```groovy title="invalid_process.nf" hl_lines="4" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        inputs:  // FEHLER: Sollte 'input' sein, nicht 'inputs'
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        """
        echo "Processing ${sample_name}" > ${sample_name}_output.txt
        """
    }

    workflow {

        // Eingabekanal erstellen
        input_ch = channel.of('sample1', 'sample2', 'sample3')

        // Prozess mit dem Eingabekanal aufrufen
        PROCESS_FILES(input_ch)
    }
    ```

#### Pipeline ausführen

Führe den Workflow erneut aus, um zu bestätigen, dass er funktioniert:

```bash
nextflow run invalid_process.nf
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `invalid_process.nf` [silly_fermi] DSL2 - revision: 961938ee2b

    executor >  local (3)
    [b7/76cd9d] PROCESS_FILES (2) | 3 of 3 ✔
    ```

### 1.3. Ungültige Variablennamen verwenden

Die Variablennamen, die du in deinen script-Blöcken verwendest, müssen gültig sein – entweder aus Eingaben abgeleitet oder aus Groovy-Code, der vor dem Skript eingefügt wird. Wenn du aber zu Beginn der Pipeline-Entwicklung mit Komplexität jonglierst, passieren leicht Fehler bei der Benennung von Variablen, und Nextflow wird dich schnell darauf hinweisen.

#### Pipeline ausführen

```bash
nextflow run no_such_var.nf
```

??? failure "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `no_such_var.nf` [gloomy_meninsky] DSL2 - revision: 0c4d3bc28c

    Error no_such_var.nf:17:39: `undefined_var` is not defined
    │  17 |     echo "Using undefined variable: ${undefined_var}" >> ${output_pref
    ╰     |                                       ^^^^^^^^^^^^^

    ERROR ~ Script compilation failed

     -- Check '.nextflow.log' file for details
    ```

Der Fehler wird zur Kompilierzeit abgefangen und zeigt direkt auf die undefinierte Variable in Zeile 17, mit einem Caret-Zeichen, das genau auf das Problem hinweist.

#### Den Code prüfen

Schauen wir uns `no_such_var.nf` an:

```groovy title="no_such_var.nf" hl_lines="17" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    input:
    val sample_name

    output:
    path "${sample_name}_processed.txt"

    script:
    // Variablen in Groovy-Code vor dem Skript definieren
    def output_prefix = "${sample_name}_processed"
    def timestamp = new Date().format("yyyy-MM-dd")

    """
    echo "Processing ${sample_name} on ${timestamp}" > ${output_prefix}.txt
    echo "Using undefined variable: ${undefined_var}" >> ${output_prefix}.txt  // FEHLER: undefined_var ist nicht definiert
    """
}

workflow {
    input_ch = channel.of('sample1', 'sample2', 'sample3')
    PROCESS_FILES(input_ch)
}
```

Die Fehlermeldung zeigt an, dass die Variable im Skript-Template nicht erkannt wird – und tatsächlich siehst du `${undefined_var}` im script-Block, aber nirgendwo definiert.

#### Den Code korrigieren

Bei einem 'No such variable'-Fehler kannst du ihn beheben, indem du die Variable definierst (durch Korrektur der Eingabevariablennamen oder Bearbeitung des Groovy-Codes vor dem Skript), oder indem du sie aus dem script-Block entfernst, wenn sie nicht benötigt wird:

=== "Danach"

    ```groovy title="no_such_var.nf" hl_lines="15-17" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        // Variablen in Groovy-Code vor dem Skript definieren
        def output_prefix = "${sample_name}_processed"
        def timestamp = new Date().format("yyyy-MM-dd")

        """
        echo "Processing ${sample_name} on ${timestamp}" > ${output_prefix}.txt
        """  // Zeile mit undefined_var entfernt
    }

    workflow {
        input_ch = channel.of('sample1', 'sample2', 'sample3')
        PROCESS_FILES(input_ch)
    }
    ```

=== "Vorher"

    ```groovy title="no_such_var.nf" hl_lines="17" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        // Variablen in Groovy-Code vor dem Skript definieren
        def output_prefix = "${sample_name}_processed"
        def timestamp = new Date().format("yyyy-MM-dd")

        """
        echo "Processing ${sample_name} on ${timestamp}" > ${output_prefix}.txt
        echo "Using undefined variable: ${undefined_var}" >> ${output_prefix}.txt  // FEHLER: undefined_var ist nicht definiert
        """
    }

    workflow {
        input_ch = channel.of('sample1', 'sample2', 'sample3')
        PROCESS_FILES(input_ch)
    }
    ```

#### Pipeline ausführen

Führe den Workflow erneut aus, um zu bestätigen, dass er funktioniert:

```bash
nextflow run no_such_var.nf
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `no_such_var.nf` [suspicious_venter] DSL2 - revision: 6ba490f7c5

    executor >  local (3)
    [21/237300] PROCESS_FILES (2) | 3 of 3 ✔
    ```

### 1.4. Falsche Verwendung von Bash-Variablen

Am Anfang mit Nextflow kann es schwierig sein, den Unterschied zwischen Nextflow-(Groovy-)Variablen und Bash-Variablen zu verstehen. Das kann eine weitere Form des Variablenfehlers erzeugen, der auftritt, wenn man versucht, Variablen im Bash-Inhalt des script-Blocks zu verwenden.

#### Pipeline ausführen

```bash
nextflow run bad_bash_var.nf
```

??? failure "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_bash_var.nf` [infallible_mandelbrot] DSL2 - revision: 0853c11080

    Error bad_bash_var.nf:13:42: `prefix` is not defined
    │  13 |     echo "Processing ${sample_name}" > ${prefix}.txt
    ╰     |                                          ^^^^^^

    ERROR ~ Script compilation failed

     -- Check '.nextflow.log' file for details
    ```

#### Den Code prüfen

Der Fehler zeigt auf Zeile 13, wo `${prefix}` verwendet wird. Schauen wir uns `bad_bash_var.nf` an, um zu sehen, was das Problem verursacht:

```groovy title="bad_bash_var.nf" hl_lines="13" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    input:
    val sample_name

    output:
    path "${sample_name}_output.txt"

    script:
    """
    prefix="${sample_name}_output"
    echo "Processing ${sample_name}" > ${prefix}.txt  # FEHLER: ${prefix} ist Groovy-Syntax, nicht Bash
    """
}
```

In diesem Beispiel definieren wir die Variable `prefix` in Bash, aber in einem Nextflow-Prozess wird die `$`-Syntax, die wir verwenden, um darauf zu verweisen (`${prefix}`), als Groovy-Variable interpretiert, nicht als Bash. Die Variable existiert nicht im Groovy-Kontext, daher erhalten wir einen 'no such variable'-Fehler.

#### Den Code korrigieren

Wenn du eine Bash-Variable verwenden möchtest, musst du das Dollarzeichen so escapen:

=== "Danach"

    ```groovy title="bad_bash_var.nf" hl_lines="13" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        """
        prefix="${sample_name}_output"
        echo "Processing ${sample_name}" > \${prefix}.txt  # Korrigiert: Dollarzeichen escaped
        """
    }

    workflow {
        input_ch = channel.of('sample1', 'sample2', 'sample3')
        PROCESS_FILES(input_ch)
    }
    ```

=== "Vorher"

    ```groovy title="bad_bash_var.nf" hl_lines="13" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        """
        prefix="${sample_name}_output"
        echo "Processing ${sample_name}" > ${prefix}.txt  # FEHLER: ${prefix} ist Groovy-Syntax, nicht Bash
        """
    }
    ```

Das teilt Nextflow mit, dies als Bash-Variable zu interpretieren.

#### Pipeline ausführen

Führe den Workflow erneut aus, um zu bestätigen, dass er funktioniert:

```bash
nextflow run bad_bash_var.nf
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_bash_var.nf` [naughty_franklin] DSL2 - revision: 58c1c83709

    executor >  local (3)
    [4e/560285] PROCESS_FILES (2) | 3 of 3 ✔
    ```

!!! tip "Groovy- vs. Bash-Variablen"

    Für einfache Variablenoperationen wie String-Verkettung oder Präfix-/Suffix-Operationen ist es in der Regel lesbarer, Groovy-Variablen im script-Abschnitt zu verwenden, anstatt Bash-Variablen im script-Block:

    ```groovy linenums="1"
    script:
    def output_prefix = "${sample_name}_processed"
    def output_file = "${output_prefix}.txt"
    """
    echo "Processing ${sample_name}" > ${output_file}
    """
    ```

    Dieser Ansatz vermeidet das Escapen von Dollarzeichen und macht den Code leichter lesbar und wartbar.

### 1.5. Anweisungen außerhalb des Workflow-Blocks

Die Nextflow-VSCode-Erweiterung hebt Probleme mit der Codestruktur hervor, die Fehler verursachen. Ein häufiges Beispiel ist das Definieren von Kanälen außerhalb des `workflow {}`-Blocks – das wird jetzt als Syntaxfehler erzwungen.

#### Pipeline ausführen

```bash
nextflow run badpractice_syntax.nf
```

??? failure "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `badpractice_syntax.nf` [intergalactic_colden] DSL2 - revision: 5e4b291bde

    Error badpractice_syntax.nf:3:1: Statements cannot be mixed with script declarations -- move statements into a process or workflow
    │   3 | input_ch = channel.of('sample1', 'sample2', 'sample3')
    ╰     | ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    ERROR ~ Script compilation failed

     -- Check '.nextflow.log' file for details
    ```

Die Fehlermeldung zeigt das Problem klar: Anweisungen (wie Kanaldefinitionen) können nicht mit Skriptdeklarationen außerhalb eines workflow- oder process-Blocks gemischt werden.

#### Den Code prüfen

Schauen wir uns `badpractice_syntax.nf` an, um zu sehen, was den Fehler verursacht:

```groovy title="badpractice_syntax.nf" hl_lines="3" linenums="1"
#!/usr/bin/env nextflow

input_ch = channel.of('sample1', 'sample2', 'sample3')  // FEHLER: Kanal außerhalb des workflow-Blocks definiert

process PROCESS_FILES {
    input:
    val sample_name

    output:
    path "${sample_name}_processed.txt"

    script:
    // Variablen in Groovy-Code vor dem Skript definieren
    def output_prefix = "${sample_name}_processed"
    def timestamp = new Date().format("yyyy-MM-dd")

    """
    echo "Processing ${sample_name} on ${timestamp}" > ${output_prefix}.txt
    """
}

workflow {
    PROCESS_FILES(input_ch)
}
```

Die VSCode-Erweiterung hebt die Variable `input_ch` ebenfalls hervor, da sie außerhalb des workflow-Blocks definiert ist:

![Non-lethal syntax error](img/nonlethal.png)

#### Den Code korrigieren

Verschiebe die Kanaldefinition in den workflow-Block:

=== "Danach"

    ```groovy title="badpractice_syntax.nf" hl_lines="21" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_processed.txt"

        script:
        // Variablen in Groovy-Code vor dem Skript definieren
        def output_prefix = "${sample_name}_processed"
        def timestamp = new Date().format("yyyy-MM-dd")

        """
        echo "Processing ${sample_name} on ${timestamp}" > ${output_prefix}.txt
        """
    }

    workflow {
        input_ch = channel.of('sample1', 'sample2', 'sample3')  // In den workflow-Block verschoben
        PROCESS_FILES(input_ch)
    }
    ```

=== "Vorher"

    ```groovy title="badpractice_syntax.nf" hl_lines="3" linenums="1"
    #!/usr/bin/env nextflow

    input_ch = channel.of('sample1', 'sample2', 'sample3')  // FEHLER: Kanal außerhalb des workflow-Blocks definiert

    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_processed.txt"

        script:
        // Variablen in Groovy-Code vor dem Skript definieren
        def output_prefix = "${sample_name}_processed"
        def timestamp = new Date().format("yyyy-MM-dd")

        """
        echo "Processing ${sample_name} on ${timestamp}" > ${output_prefix}.txt
        """
    }

    workflow {
        PROCESS_FILES(input_ch)
    }
    ```

#### Pipeline ausführen

Führe den Workflow erneut aus, um zu bestätigen, dass die Korrektur funktioniert:

```bash
nextflow run badpractice_syntax.nf
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `badpractice_syntax.nf` [naughty_ochoa] DSL2 - revision: 5e4b291bde

    executor >  local (3)
    [6a/84a608] PROCESS_FILES (2) | 3 of 3 ✔
    ```

Definiere deine Eingabekanäle immer innerhalb des workflow-Blocks und befolge generell alle anderen Empfehlungen der Erweiterung.

### Fazit

Du kannst Syntaxfehler systematisch identifizieren und beheben, indem du Nextflow-Fehlermeldungen und visuelle IDE-Indikatoren nutzt. Häufige Syntaxfehler sind fehlende geschweifte Klammern, falsche Prozess-Keywords, undefinierte Variablen und die falsche Verwendung von Bash- vs. Nextflow-Variablen. Die VSCode-Erweiterung hilft dabei, viele dieser Fehler vor der Laufzeit zu erkennen. Mit diesen Syntaxdebugging-Fähigkeiten in deinem Werkzeugkasten kannst du die häufigsten Nextflow-Syntaxfehler schnell beheben und dich komplexeren Laufzeitproblemen widmen.

### Wie geht es weiter?

Lerne, komplexere Kanalstrukturfehler zu debuggen, die auch dann auftreten, wenn die Syntax korrekt ist.

---

## 2. Kanalstrukturfehler

Kanalstrukturfehler sind subtiler als Syntaxfehler, weil der Code syntaktisch korrekt ist, aber die Datenformen nicht mit dem übereinstimmen, was Prozesse erwarten. Nextflow versucht, die Pipeline auszuführen, stellt aber möglicherweise fest, dass die Anzahl der Eingaben nicht mit den Erwartungen übereinstimmt, und schlägt fehl. Diese Fehler treten typischerweise erst zur Laufzeit auf und erfordern ein Verständnis der Daten, die durch deinen Workflow fließen.

!!! tip "Kanäle mit `.view()` debuggen"

    Denke in diesem Abschnitt daran, dass du den `.view()`-Operator verwenden kannst, um den Kanalinhalt an jedem Punkt in deinem Workflow zu inspizieren. Das ist eines der mächtigsten Debugging-Tools zum Verstehen von Kanalstrukturproblemen. Wir erkunden diese Technik in Abschnitt 2.4 im Detail, aber du kannst sie gerne schon beim Durcharbeiten der Beispiele verwenden.

    ```groovy
    my_channel.view()  // Zeigt, was durch den Kanal fließt
    ```

### 2.1. Falsche Anzahl von Eingabekanälen

Dieser Fehler tritt auf, wenn du eine andere Anzahl von Kanälen übergibst, als ein Prozess erwartet.

#### Pipeline ausführen

```bash
nextflow run bad_number_inputs.nf
```

??? failure "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_number_inputs.nf` [happy_swartz] DSL2 - revision: d83e58dcd3

    Error bad_number_inputs.nf:23:5: Incorrect number of call arguments, expected 1 but received 2
    │  23 |     PROCESS_FILES(samples_ch, files_ch)
    ╰     |     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    ERROR ~ Script compilation failed

     -- Check '.nextflow.log' file for details
    ```

#### Den Code prüfen

Die Fehlermeldung besagt klar, dass der Aufruf 1 Argument erwartet hat, aber 2 erhalten hat, und zeigt auf Zeile 23. Schauen wir uns `bad_number_inputs.nf` an:

```groovy title="bad_number_inputs.nf" hl_lines="5 23" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    input:
        val sample_name  // Prozess erwartet nur 1 Eingabe

    output:
        path "${sample_name}_output.txt"

    script:
    """
    echo "Processing ${sample_name}" > ${sample_name}_output.txt
    """
}

workflow {

    // Zwei separate Kanäle erstellen
    samples_ch = channel.of('sample1', 'sample2', 'sample3')
    files_ch = channel.of('file1.txt', 'file2.txt', 'file3.txt')

    // FEHLER: 2 Kanäle übergeben, aber Prozess erwartet nur 1
    PROCESS_FILES(samples_ch, files_ch)
}
```

Du siehst den nicht übereinstimmenden `PROCESS_FILES`-Aufruf, der mehrere Eingabekanäle liefert, obwohl der Prozess nur einen definiert. Die VSCode-Erweiterung unterstreicht den Prozessaufruf ebenfalls rot und zeigt eine Diagnosemeldung beim Darüberfahren mit der Maus:

![Incorrect number of args message](img/incorrect_num_args.png)

#### Den Code korrigieren

Für dieses spezifische Beispiel erwartet der Prozess einen einzelnen Kanal und benötigt den zweiten Kanal nicht, daher können wir es beheben, indem wir nur den `samples_ch`-Kanal übergeben:

=== "Danach"

    ```groovy title="bad_number_inputs.nf" hl_lines="23" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
            val sample_name  // Prozess erwartet nur 1 Eingabe

        output:
            path "${sample_name}_output.txt"

        script:
        """
        echo "Processing ${sample_name}" > ${sample_name}_output.txt
        """
    }

    workflow {

        // Zwei separate Kanäle erstellen
        samples_ch = channel.of('sample1', 'sample2', 'sample3')
        files_ch = channel.of('file1.txt', 'file2.txt', 'file3.txt')

        // Korrigiert: Nur den Kanal übergeben, den der Prozess erwartet
        PROCESS_FILES(samples_ch)
    }
    ```

=== "Vorher"

    ```groovy title="bad_number_inputs.nf" hl_lines="5 23" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
            val sample_name  // Prozess erwartet nur 1 Eingabe

        output:
            path "${sample_name}_output.txt"

        script:
        """
        echo "Processing ${sample_name}" > ${sample_name}_output.txt
        """
    }

    workflow {

        // Zwei separate Kanäle erstellen
        samples_ch = channel.of('sample1', 'sample2', 'sample3')
        files_ch = channel.of('file1.txt', 'file2.txt', 'file3.txt')

        // FEHLER: 2 Kanäle übergeben, aber Prozess erwartet nur 1
        PROCESS_FILES(samples_ch, files_ch)
    }
    ```

#### Pipeline ausführen

```bash
nextflow run bad_number_inputs.nf
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_number_inputs.nf` [big_euler] DSL2 - revision: e302bd87be

    executor >  local (3)
    [48/497f7b] PROCESS_FILES (3) | 3 of 3 ✔
    ```

Häufiger als in diesem Beispiel wirst du zusätzliche Eingaben zu einem Prozess hinzufügen und vergessen, den workflow-Aufruf entsprechend zu aktualisieren, was zu dieser Art von Fehler führen kann. Glücklicherweise ist das einer der leichter verständlichen und behebbaren Fehler, da die Fehlermeldung die Diskrepanz klar beschreibt.

### 2.2. Kanalerschöpfung (Prozess läuft seltener als erwartet)

Einige Kanalstrukturfehler sind viel subtiler und erzeugen überhaupt keine Fehler. Wahrscheinlich der häufigste davon spiegelt eine Herausforderung wider, mit der neue Nextflow-Nutzer\*innen konfrontiert sind: zu verstehen, dass queue channels erschöpft werden und keine Elemente mehr haben können, was dazu führt, dass der Workflow vorzeitig endet.

#### Pipeline ausführen

```bash
nextflow run exhausted.nf
```

??? success "Befehlsausgabe"

```console title="Exhausted channel output"
 N E X T F L O W   ~  version 25.10.2

Launching `exhausted.nf` [extravagant_gauss] DSL2 - revision: 08cff7ba2a

executor >  local (1)
[bd/f61fff] PROCESS_FILES (1) [100%] 1 of 1 ✔
```

Dieser Workflow wird ohne Fehler abgeschlossen, verarbeitet aber nur eine einzige Probe!

#### Den Code prüfen

Schauen wir uns `exhausted.nf` an, um zu sehen, ob das korrekt ist:

```groovy title="exhausted.nf" hl_lines="23 24" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    input:
    val reference
    val sample_name

    output:
    path "${output_prefix}.txt"

    script:
    // Variablen in Groovy-Code vor dem Skript definieren
    output_prefix = "${reference}_${sample_name}"
    def timestamp = new Date().format("yyyy-MM-dd")

    """
    echo "Processing ${sample_name} on ${timestamp}" > ${output_prefix}.txt
    """
}

workflow {

    reference_ch = channel.of('baseline_reference')
    input_ch = channel.of('sample1', 'sample2', 'sample3')

    PROCESS_FILES(reference_ch, input_ch)
}
```

Der Prozess läuft nur einmal statt dreimal, weil der `reference_ch`-Kanal ein queue channel ist, der nach der ersten Prozessausführung erschöpft wird. Wenn ein Kanal erschöpft ist, stoppt der gesamte Prozess, auch wenn andere Kanäle noch Elemente haben.

Das ist ein häufiges Muster, bei dem du eine einzelne Referenzdatei hast, die für mehrere Proben wiederverwendet werden muss. Die Lösung besteht darin, den Referenzkanal in einen value channel umzuwandeln, der unbegrenzt wiederverwendet werden kann.

#### Den Code korrigieren

Es gibt je nach Anzahl der betroffenen Dateien verschiedene Möglichkeiten, das zu beheben.

**Option 1**: Du hast eine einzelne Referenzdatei, die du häufig wiederverwendest. Du kannst einfach einen value channel-Typ erstellen, der immer wieder verwendet werden kann. Es gibt drei Möglichkeiten:

**1a** `channel.value()` verwenden:

```groovy title="exhausted.nf (fixed - Option 1a)" hl_lines="2" linenums="21"
workflow {
    reference_ch = channel.value('baseline_reference')  // Value channel kann wiederverwendet werden
    input_ch = channel.of('sample1', 'sample2', 'sample3')

    PROCESS_FILES(reference_ch, input_ch)
}
```

**1b** Den `first()`-[Operator](https://www.nextflow.io/docs/latest/reference/operator.html#first) verwenden:

```groovy title="exhausted.nf (fixed - Option 1b)" hl_lines="2" linenums="21"
workflow {
    reference_ch = channel.of('baseline_reference').first()  // In value channel umwandeln
    input_ch = channel.of('sample1', 'sample2', 'sample3')

    PROCESS_FILES(reference_ch, input_ch)
}
```

**1c.** Den `collect()`-[Operator](https://www.nextflow.io/docs/latest/reference/operator.html#collect) verwenden:

```groovy title="exhausted.nf (fixed - Option 1c)" hl_lines="2" linenums="21"
workflow {
    reference_ch = channel.of('baseline_reference').collect()  // In value channel umwandeln
    input_ch = channel.of('sample1', 'sample2', 'sample3')

    PROCESS_FILES(reference_ch, input_ch)
}
```

**Option 2**: In komplexeren Szenarien, z. B. wenn du mehrere Referenzdateien für alle Proben im Probenkanal hast, kannst du den `combine`-Operator verwenden, um einen neuen Kanal zu erstellen, der die beiden Kanäle zu Tupeln kombiniert:

```groovy title="exhausted.nf (fixed - Option 2)" hl_lines="4" linenums="21"
workflow {
    reference_ch = channel.of('baseline_reference','other_reference')
    input_ch = channel.of('sample1', 'sample2', 'sample3')
    combined_ch = reference_ch.combine(input_ch)  // Erstellt kartesisches Produkt

    PROCESS_FILES(combined_ch)
}
```

Der `.combine()`-Operator erzeugt ein kartesisches Produkt der beiden Kanäle, sodass jedes Element in `reference_ch` mit jedem Element in `input_ch` gepaart wird. Das ermöglicht es dem Prozess, für jede Probe zu laufen und dabei die Referenz zu verwenden.

Dafür muss die Prozesseingabe angepasst werden. In unserem Beispiel müsste der Anfang der Prozessdefinition wie folgt angepasst werden:

```groovy title="exhausted.nf (fixed - Option 2)" hl_lines="5" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    input:
        tuple val(reference), val(sample_name)
```

Dieser Ansatz ist nicht in allen Situationen geeignet.

#### Pipeline ausführen

Probiere eine der obigen Korrekturen aus und führe den Workflow erneut aus:

```bash
nextflow run exhausted.nf
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `exhausted.nf` [maniac_leavitt] DSL2 - revision: f372a56a7d

    executor >  local (3)
    [80/0779e9] PROCESS_FILES (3) | 3 of 3 ✔
    ```

Du solltest jetzt sehen, dass alle drei Proben verarbeitet werden, nicht nur eine.

### 2.3. Falsche Kanalinhaltsstruktur

Wenn Workflows eine gewisse Komplexität erreichen, kann es schwierig sein, die internen Strukturen jedes Kanals im Blick zu behalten, und es kommt häufig zu Diskrepanzen zwischen dem, was der Prozess erwartet, und dem, was der Kanal tatsächlich enthält. Das ist subtiler als das zuvor besprochene Problem, bei dem die Anzahl der Kanäle falsch war. In diesem Fall kann die korrekte Anzahl von Eingabekanälen vorhanden sein, aber die interne Struktur eines oder mehrerer dieser Kanäle stimmt nicht mit dem überein, was der Prozess erwartet.

#### Pipeline ausführen

```bash
nextflow run bad_channel_shape.nf
```

??? failure "Befehlsausgabe"

    ```console
    Launching `bad_channel_shape.nf` [hopeful_pare] DSL2 - revision: ffd66071a1

    executor >  local (3)
    executor >  local (3)
    [3f/c2dcb3] PROCESS_FILES (3) [  0%] 0 of 3 ✘
    ERROR ~ Error executing process > 'PROCESS_FILES (1)'

    Caused by:
      Missing output file(s) `[sample1, file1.txt]_output.txt` expected by process `PROCESS_FILES (1)`


    Command executed:

      echo "Processing [sample1, file1.txt]" > [sample1, file1.txt]_output.txt

    Command exit status:
      0

    Command output:
      (empty)

    Work dir:
      /workspaces/training/side-quests/debugging/work/d6/1fb69d1d93300bbc9d42f1875b981e

    Tip: when you have fixed the problem you can continue the execution adding the option `-resume` to the run command line

    -- Check '.nextflow.log' file for details
    ```

#### Den Code prüfen

Die eckigen Klammern in der Fehlermeldung geben den Hinweis – der Prozess behandelt das Tupel als einzelnen Wert, was nicht gewollt ist. Schauen wir uns `bad_channel_shape.nf` an:

```groovy title="bad_channel_shape.nf" hl_lines="5 20-22" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    input:
        val sample_name  // Erwartet einzelnen Wert, erhält Tupel

    output:
        path "${sample_name}_output.txt"

    script:
    """
    echo "Processing ${sample_name}" > ${sample_name}_output.txt
    """
}

workflow {

    // Kanal gibt Tupel aus, aber Prozess erwartet einzelne Werte
    input_ch = channel.of(
      ['sample1', 'file1.txt'],
      ['sample2', 'file2.txt'],
      ['sample3', 'file3.txt']
    )
    PROCESS_FILES(input_ch)
}
```

Du siehst, dass wir einen Kanal aus Tupeln erzeugen: `['sample1', 'file1.txt']`, aber der Prozess erwartet einen einzelnen Wert, `val sample_name`. Der ausgeführte Befehl zeigt, dass der Prozess versucht, eine Datei namens `[sample3, file3.txt]_output.txt` zu erstellen, was nicht die beabsichtigte Ausgabe ist.

#### Den Code korrigieren

Um das zu beheben: Wenn der Prozess beide Eingaben benötigt, könnten wir den Prozess so anpassen, dass er ein Tupel akzeptiert:

=== "Option 1: Tupel im Prozess akzeptieren"

    === "Danach"

        ```groovy title="bad_channel_shape.nf" hl_lines="5"  linenums="1"
        #!/usr/bin/env nextflow

        process PROCESS_FILES {
            input:
                tuple val(sample_name), val(file_name)  // Korrigiert: Tupel akzeptieren

            output:
                path "${sample_name}_output.txt"

            script:
            """
            echo "Processing ${sample_name}" > ${sample_name}_output.txt
            """
        }

        workflow {

            // Kanal gibt Tupel aus, aber Prozess erwartet einzelne Werte
            input_ch = channel.of(
              ['sample1', 'file1.txt'],
              ['sample2', 'file2.txt'],
              ['sample3', 'file3.txt']
            )
            PROCESS_FILES(input_ch)
        }
        ```

    === "Vorher"

        ```groovy title="bad_channel_shape.nf" hl_lines="5" linenums="1"
        #!/usr/bin/env nextflow

        process PROCESS_FILES {
            input:
                val sample_name  // Erwartet einzelnen Wert, erhält Tupel

            output:
                path "${sample_name}_output.txt"

            script:
            """
            echo "Processing ${sample_name}" > ${sample_name}_output.txt
            """
        }

        workflow {

            // Kanal gibt Tupel aus, aber Prozess erwartet einzelne Werte
            input_ch = channel.of(
              ['sample1', 'file1.txt'],
              ['sample2', 'file2.txt'],
              ['sample3', 'file3.txt']
            )
            PROCESS_FILES(input_ch)
        }
        ```

=== "Option 2: Erstes Element extrahieren"

    === "Danach"

        ```groovy title="bad_channel_shape.nf" hl_lines="9" linenums="16"
        workflow {

            // Kanal gibt Tupel aus, aber Prozess erwartet einzelne Werte
            input_ch = channel.of(
              ['sample1', 'file1.txt'],
              ['sample2', 'file2.txt'],
              ['sample3', 'file3.txt']
            )
            PROCESS_FILES(input_ch.map { it[0] })  // Korrigiert: Erstes Element extrahieren
        }
        ```

    === "Vorher"

        ```groovy title="bad_channel_shape.nf" hl_lines="9" linenums="16"
        workflow {

            // Kanal gibt Tupel aus, aber Prozess erwartet einzelne Werte
            input_ch = channel.of(
              ['sample1', 'file1.txt'],
              ['sample2', 'file2.txt'],
              ['sample3', 'file3.txt']
            )
            PROCESS_FILES(input_ch)
        }
        ```

#### Pipeline ausführen

Wähle eine der Lösungen und führe den Workflow erneut aus:

```bash
nextflow run bad_channel_shape.nf
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_channel_shape.nf` [clever_thompson] DSL2 - revision: 8cbcae3746

    executor >  local (3)
    [bb/80a958] PROCESS_FILES (2) | 3 of 3 ✔
    ```

### 2.4. Kanal-Debugging-Techniken

#### `.view()` zur Kanalinspektion verwenden

Das mächtigste Debugging-Tool für Kanäle ist der `.view()`-Operator. Mit `.view()` kannst du die Form deiner Kanäle in allen Phasen verstehen und beim Debuggen helfen.

#### Pipeline ausführen

Führe `bad_channel_shape_viewed.nf` aus, um das in Aktion zu sehen:

```bash
nextflow run bad_channel_shape_viewed.nf
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_channel_shape_viewed.nf` [maniac_poisson] DSL2 - revision: b4f24dc9da

    executor >  local (3)
    [c0/db76b3] PROCESS_FILES (3) [100%] 3 of 3 ✔
    Channel content: [sample1, file1.txt]
    Channel content: [sample2, file2.txt]
    Channel content: [sample3, file3.txt]
    After mapping: sample1
    After mapping: sample2
    After mapping: sample3
    ```

#### Den Code prüfen

Schauen wir uns `bad_channel_shape_viewed.nf` an, um zu sehen, wie `.view()` verwendet wird:

```groovy title="bad_channel_shape_viewed.nf" linenums="16" hl_lines="9 11"
workflow {

    // Kanal gibt Tupel aus, aber Prozess erwartet einzelne Werte
    input_ch = channel.of(
      ['sample1', 'file1.txt'],
      ['sample2', 'file2.txt'],
      ['sample3', 'file3.txt']
    )
    .view { "Channel content: $it" }  // Debug: Ursprünglichen Kanalinhalt anzeigen
    .map { tuple -> tuple[0] }        // Transformation: Erstes Element extrahieren
    .view { "After mapping: $it" }    // Debug: Transformierten Kanalinhalt anzeigen

    PROCESS_FILES(input_ch)
}
```

#### Den Code verbessern

Um zu vermeiden, dass du in Zukunft übermäßig viele `.view()`-Operationen verwenden musst, um den Kanalinhalt zu verstehen, empfiehlt es sich, Kommentare hinzuzufügen:

```groovy title="bad_channel_shape_viewed.nf (with comments)" linenums="16" hl_lines="8 9"
workflow {

    // Kanal gibt Tupel aus, aber Prozess erwartet einzelne Werte
    input_ch = channel.of(
            ['sample1', 'file1.txt'],
            ['sample2', 'file2.txt'],
            ['sample3', 'file3.txt'],
        ) // [sample_name, file_name]
        .map { tuple -> tuple[0] } // sample_name

    PROCESS_FILES(input_ch)
}
```

Das wird wichtiger, wenn deine Workflows komplexer werden und die Kanalstruktur undurchsichtiger wird.

#### Pipeline ausführen

```bash
nextflow run bad_channel_shape_viewed.nf
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_channel_shape_viewed.nf` [marvelous_koch] DSL2 - revision: 03e79cdbad

    executor >  local (3)
    [ff/d67cec] PROCESS_FILES (2) | 3 of 3 ✔
    Channel content: [sample1, file1.txt]
    Channel content: [sample2, file2.txt]
    Channel content: [sample3, file3.txt]
    After mapping: sample1
    After mapping: sample2
    After mapping: sample3
    ```

### Fazit

Viele Kanalstrukturfehler können mit gültigem Nextflow-Syntax entstehen. Du kannst Kanalstrukturfehler debuggen, indem du den Datenfluss verstehst, `.view()`-Operatoren zur Inspektion verwendest und Fehlermeldungsmuster erkennst, wie eckige Klammern, die auf unerwartete Tupelstrukturen hinweisen.

### Wie geht es weiter?

Lerne über Fehler, die durch Prozessdefinitionen entstehen.

---

## 3. Prozessstrukturfehler

Die meisten Fehler, die du bei Prozessen begegnest, hängen mit Fehlern beim Formulieren des Befehls oder mit Problemen der zugrunde liegenden Software zusammen. Ähnlich wie bei den Kanalproblemen oben kannst du jedoch Fehler in der Prozessdefinition machen, die keine Syntaxfehler sind, aber zur Laufzeit Fehler verursachen.

### 3.1. Fehlende Ausgabedateien

Ein häufiger Fehler beim Schreiben von Prozessen ist eine Diskrepanz zwischen dem, was der Prozess erwartet, und dem, was tatsächlich erzeugt wird.

#### Pipeline ausführen

```bash
nextflow run missing_output.nf
```

??? failure "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `missing_output.nf` [zen_stone] DSL2 - revision: 37ff61f926

    executor >  local (3)
    executor >  local (3)
    [fd/2642e9] process > PROCESS_FILES (2) [ 66%] 2 of 3, failed: 2
    ERROR ~ Error executing process > 'PROCESS_FILES (3)'

    Caused by:
      Missing output file(s) `sample3.txt` expected by process `PROCESS_FILES (3)`


    Command executed:

      echo "Processing sample3" > sample3_output.txt

    Command exit status:
      0

    Command output:
      (empty)

    Work dir:
      /workspaces/training/side-quests/debugging/work/02/9604d49fb8200a74d737c72a6c98ed

    Tip: when you have fixed the problem you can continue the execution adding the option `-resume` to the run command line

    -- Check '.nextflow.log' file for details
    ```

#### Den Code prüfen

Die Fehlermeldung zeigt, dass der Prozess erwartet hat, eine Ausgabedatei namens `sample3.txt` zu erzeugen, aber das Skript tatsächlich `sample3_output.txt` erstellt. Schauen wir uns die Prozessdefinition in `missing_output.nf` an:

```groovy title="missing_output.nf" linenums="3" hl_lines="6 10"
process PROCESS_FILES {
    input:
    val sample_name

    output:
    path "${sample_name}.txt"  // Erwartet: sample3.txt

    script:
    """
    echo "Processing ${sample_name}" > ${sample_name}_output.txt  // Erstellt: sample3_output.txt
    """
}
```

Du siehst, dass es eine Diskrepanz zwischen dem Ausgabedateinamen im `output:`-Block und dem im Skript verwendeten gibt. Diese Diskrepanz führt zum Scheitern des Prozesses. Wenn du auf diese Art von Fehler stößt, überprüfe, ob die Ausgaben zwischen deiner Prozessdefinition und deinem output-Block übereinstimmen.

Wenn das Problem noch nicht klar ist, überprüfe das work directory selbst, um die tatsächlich erstellten Ausgabedateien zu identifizieren:

```bash
❯ ls -h work/02/9604d49fb8200a74d737c72a6c98ed
sample3_output.txt
```

In diesem Beispiel würde uns das zeigen, dass ein `_output`-Suffix in den Ausgabedateinamen eingebaut wird, entgegen unserer `output:`-Definition.

#### Den Code korrigieren

Behebe die Diskrepanz, indem du den Ausgabedateinamen konsistent machst:

=== "Danach"

    ```groovy title="missing_output.nf" hl_lines="6 10" linenums="3"
    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_output.txt"  // Korrigiert: Skriptausgabe angepasst

        script:
        """
        echo "Processing ${sample_name}" > ${sample_name}_output.txt
        """
    }
    ```

=== "Vorher"

    ```groovy title="missing_output.nf" hl_lines="6 10" linenums="3"
    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}.txt"  // Erwartet: sample3.txt

        script:
        """
        echo "Processing ${sample_name}" > ${sample_name}_output.txt  // Erstellt: sample3_output.txt
        """
    }
    ```

#### Pipeline ausführen

```bash
nextflow run missing_output.nf
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `missing_output.nf` [elated_hamilton] DSL2 - revision: 961938ee2b

    executor >  local (3)
    [16/1c437c] PROCESS_FILES (3) | 3 of 3 ✔
    ```

### 3.2. Fehlende Software

Eine weitere Klasse von Fehlern entsteht durch Fehler bei der Software-Bereitstellung. `missing_software.nf` ist ein syntaktisch gültiger Workflow, der aber auf externe Software angewiesen ist, um den `cowpy`-Befehl bereitzustellen.

#### Pipeline ausführen

```bash
nextflow run missing_software.nf
```

??? failure "Befehlsausgabe"

    ```console hl_lines="12 18"
    ERROR ~ Error executing process > 'PROCESS_FILES (3)'

    Caused by:
      Process `PROCESS_FILES (3)` terminated with an error exit status (127)


    Command executed:

      cowpy sample3 > sample3_output.txt

    Command exit status:
      127

    Command output:
      (empty)

    Command error:
      .command.sh: line 2: cowpy: command not found

    Work dir:
      /workspaces/training/side-quests/debugging/work/82/42a5bfb60c9c6ee63ebdbc2d51aa6e

    Tip: you can try to figure out what's wrong by changing to the process work directory and showing the script file named `.command.sh`

    -- Check '.nextflow.log' file for details
    ```

Der Prozess hat keinen Zugriff auf den angegebenen Befehl. Manchmal liegt das daran, dass ein Skript im `bin`-Verzeichnis des Workflows vorhanden ist, aber nicht ausführbar gemacht wurde. In anderen Fällen ist die Software nicht im Container oder in der Umgebung installiert, in der der Workflow läuft.

#### Den Code prüfen

Achte auf den Exit-Code `127` – er sagt dir genau, was das Problem ist. Schauen wir uns `missing_software.nf` an:

```groovy title="missing_software.nf" linenums="3" hl_lines="3"
process PROCESS_FILES {

    container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

    input:
    val sample_name

    output:
    path "${sample_name}_output.txt"

    script:
    """
    cowpy ${sample_name} > ${sample_name}_output.txt
    """
}
```

#### Den Code korrigieren

Wir waren hier etwas unehrlich – am Code selbst ist eigentlich nichts falsch. Wir müssen nur die notwendige Konfiguration angeben, um den Prozess so auszuführen, dass er Zugriff auf den betreffenden Befehl hat. In diesem Fall hat der Prozess eine Container-Definition, also müssen wir den Workflow nur mit aktiviertem Docker ausführen.

#### Pipeline ausführen

Wir haben ein Docker-Profil in `nextflow.config` für dich eingerichtet, sodass du den Workflow ausführen kannst mit:

```bash
nextflow run missing_software.nf -profile docker
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `missing_software.nf` [awesome_stonebraker] DSL2 - revision: 0296d12839

    executor >  local (3)
    [38/ab20d1] PROCESS_FILES (1) | 3 of 3 ✔
    ```

!!! note "Hinweis"

    Um mehr darüber zu erfahren, wie Nextflow Container verwendet, siehe [Hello Nextflow](../hello_nextflow/05_hello_containers.md)

### 3.3. Fehlerhafte Ressourcenkonfiguration

Im Produktionsbetrieb wirst du Ressourcen für deine Prozesse konfigurieren. Zum Beispiel definiert `memory` die maximale Speichermenge, die deinem Prozess zur Verfügung steht. Wenn der Prozess diese überschreitet, beendet dein Scheduler den Prozess typischerweise und gibt einen Exit-Code von `137` zurück. Das können wir hier nicht demonstrieren, weil wir den `local`-Executor verwenden, aber wir können etwas Ähnliches mit `time` zeigen.

#### Pipeline ausführen

`bad_resources.nf` hat eine Prozesskonfiguration mit einer unrealistischen Zeitbegrenzung von 1 Millisekunde:

```bash
nextflow run bad_resources.nf -profile docker
```

??? failure "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_resources.nf` [disturbed_elion] DSL2 - revision: 27d2066e86

    executor >  local (3)
    [c0/ded8e1] PROCESS_FILES (3) | 0 of 3 ✘
    ERROR ~ Error executing process > 'PROCESS_FILES (2)'

    Caused by:
      Process exceeded running time limit (1ms)

    Command executed:

      cowpy sample2 > sample2_output.txt

    Command exit status:
      -

    Command output:
      (empty)

    Work dir:
      /workspaces/training/side-quests/debugging/work/53/f0a4cc56d6b3dc2a6754ff326f1349

    Container:
      community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273

    Tip: you can replicate the issue by changing to the process work dir and entering the command `bash .command.run`

     -- Check '.nextflow.log' file for details
    ```

#### Den Code prüfen

Schauen wir uns `bad_resources.nf` an:

```groovy title="bad_resources.nf" linenums="3" hl_lines="3"
process PROCESS_FILES {

    time '1 ms'  // FEHLER: Unrealistische Zeitbegrenzung

    input:
    val sample_name

    output:
    path "${sample_name}_output.txt"

    script:
    """
    sleep 1  // Dauert 1 Sekunde, aber Zeitlimit ist 1ms
    cowpy ${sample_name} > ${sample_name}_output.txt
    """
}
```

Wir wissen, dass der Prozess länger als eine Sekunde dauern wird (wir haben ein sleep eingebaut, um das sicherzustellen), aber der Prozess ist so eingestellt, dass er nach 1 Millisekunde abbricht. Jemand war mit seiner Konfiguration etwas unrealistisch!

#### Den Code korrigieren

Erhöhe das Zeitlimit auf einen realistischen Wert:

=== "Danach"

    ```groovy title="bad_resources.nf" hl_lines="3" linenums="3"
    process PROCESS_FILES {

        time '100 s'  // Korrigiert: Realistisches Zeitlimit

        input:
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        """
        sleep 1
        cowpy ${sample_name} > ${sample_name}_output.txt
        """
    }
    ```

=== "Vorher"

    ```groovy title="bad_resources.nf" hl_lines="3" linenums="3"
    process PROCESS_FILES {

        time '1 ms'  // FEHLER: Unrealistische Zeitbegrenzung

        input:
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        """
        sleep 1  // Dauert 1 Sekunde, aber Zeitlimit ist 1ms
        cowpy ${sample_name} > ${sample_name}_output.txt
        """
    }
    ```

#### Pipeline ausführen

```bash
nextflow run bad_resources.nf -profile docker
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_resources.nf` [friendly_mcclintock] DSL2 - revision: 381567d2c1

    executor >  local (3)
    [c2/9b4c41] PROCESS_FILES (3) | 3 of 3 ✔
    ```

Wenn du deine Fehlermeldungen sorgfältig liest, sollten dich solche Fehler nicht lange beschäftigen. Stelle aber sicher, dass du die Ressourcenanforderungen der von dir ausgeführten Befehle verstehst, damit du deine Ressourcen-Direktiven entsprechend konfigurieren kannst.

### 3.4. Prozess-Debugging-Techniken

Wenn Prozesse fehlschlagen oder sich unerwartet verhalten, brauchst du systematische Techniken, um zu untersuchen, was schiefgelaufen ist. Das work directory enthält alle Informationen, die du zum Debuggen der Prozessausführung benötigst.

#### Work Directory-Inspektion verwenden

Das mächtigste Debugging-Tool für Prozesse ist die Untersuchung des work directory. Wenn ein Prozess fehlschlägt, erstellt Nextflow ein work directory für diese spezifische Prozessausführung, das alle Dateien enthält, die du brauchst, um zu verstehen, was passiert ist.

#### Pipeline ausführen

Verwenden wir das `missing_output.nf`-Beispiel von früher, um die work directory-Inspektion zu demonstrieren (erzeuge bei Bedarf erneut eine Ausgabebenennungsdiskrepanz):

```bash
nextflow run missing_output.nf
```

??? failure "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `missing_output.nf` [irreverent_payne] DSL2 - revision: 3d5117f7e2

    executor >  local (3)
    [5d/d544a4] PROCESS_FILES (2) | 0 of 3 ✘
    ERROR ~ Error executing process > 'PROCESS_FILES (1)'

    Caused by:
      Missing output file(s) `sample1.txt` expected by process `PROCESS_FILES (1)`

    Command executed:

      echo "Processing sample1" > sample1_output.txt

    Command exit status:
      0

    Command output:
      (empty)

    Work dir:
      /workspaces/training/side-quests/debugging/work/1e/2011154d0b0f001cd383d7364b5244

    Tip: you can replicate the issue by changing to the process work dir and entering the command `bash .command.run`

     -- Check '.nextflow.log' file for details
    ```

#### Das work directory prüfen

Wenn du diesen Fehler erhältst, enthält das work directory alle Debugging-Informationen. Finde den work directory-Pfad aus der Fehlermeldung und untersuche seinen Inhalt:

```bash
# Work directory aus der Fehlermeldung finden
ls work/02/9604d49fb8200a74d737c72a6c98ed/
```

Du kannst dann die wichtigsten Dateien untersuchen:

##### Das Befehlsskript prüfen

Die `.command.sh`-Datei zeigt genau, welcher Befehl ausgeführt wurde:

```bash
# Ausgeführten Befehl anzeigen
cat work/02/9604d49fb8200a74d737c72a6c98ed/.command.sh
```

Das zeigt:

- **Variablensubstitution**: Ob Nextflow-Variablen korrekt expandiert wurden
- **Dateipfade**: Ob Eingabedateien korrekt gefunden wurden
- **Befehlsstruktur**: Ob die Skriptsyntax korrekt ist

Häufige Probleme, auf die du achten solltest:

- **Fehlende Anführungszeichen**: Variablen mit Leerzeichen benötigen korrekte Anführungszeichen
- **Falsche Dateipfade**: Eingabedateien, die nicht existieren oder am falschen Ort sind
- **Falsche Variablennamen**: Tippfehler in Variablenreferenzen
- **Fehlende Umgebungseinrichtung**: Befehle, die von bestimmten Umgebungen abhängen

##### Fehlerausgabe prüfen

Die `.command.err`-Datei enthält die tatsächlichen Fehlermeldungen:

```bash
# Fehlerausgabe anzeigen
cat work/02/9604d49fb8200a74d737c72a6c98ed/.command.err
```

Diese Datei zeigt:

- **Exit-Codes**: 127 (Befehl nicht gefunden), 137 (beendet) usw.
- **Berechtigungsfehler**: Probleme beim Dateizugriff
- **Softwarefehler**: Anwendungsspezifische Fehlermeldungen
- **Ressourcenfehler**: Speicher-/Zeitlimit überschritten

##### Standardausgabe prüfen

Die `.command.out`-Datei zeigt, was dein Befehl produziert hat:

```bash
# Standardausgabe anzeigen
cat work/02/9604d49fb8200a74d737c72a6c98ed/.command.out
```

Das hilft zu überprüfen:

- **Erwartete Ausgabe**: Ob der Befehl die richtigen Ergebnisse produziert hat
- **Teilausführung**: Ob der Befehl gestartet, aber auf halbem Weg fehlgeschlagen ist
- **Debug-Informationen**: Jegliche Diagnosausgabe aus deinem Skript

##### Den Exit-Code prüfen

Die `.exitcode`-Datei enthält den Exit-Code des Prozesses:

```bash
# Exit-Code anzeigen
cat work/*/*/.exitcode
```

Häufige Exit-Codes und ihre Bedeutungen:

- **Exit-Code 127**: Befehl nicht gefunden – Software-Installation prüfen
- **Exit-Code 137**: Prozess beendet – Speicher-/Zeitlimits prüfen

##### Dateiexistenz prüfen

Wenn Prozesse aufgrund fehlender Ausgabedateien fehlschlagen, prüfe, welche Dateien tatsächlich erstellt wurden:

```bash
# Alle Dateien im work directory auflisten
ls -la work/02/9604d49fb8200a74d737c72a6c98ed/
```

Das hilft zu identifizieren:

- **Dateibenennungsdiskrepanzen**: Ausgabedateien mit anderen Namen als erwartet
- **Berechtigungsprobleme**: Dateien, die nicht erstellt werden konnten
- **Pfadprobleme**: Dateien, die in falschen Verzeichnissen erstellt wurden

In unserem früheren Beispiel bestätigte uns das, dass zwar unsere erwartete `sample3.txt` nicht vorhanden war, aber `sample3_output.txt` schon:

```bash
❯ ls -h work/02/9604d49fb8200a74d737c72a6c98ed
sample3_output.txt
```

### Fazit

Prozess-Debugging erfordert die Untersuchung von work directories, um zu verstehen, was schiefgelaufen ist. Wichtige Dateien sind `.command.sh` (das ausgeführte Skript), `.command.err` (Fehlermeldungen) und `.command.out` (Standardausgabe). Exit-Codes wie 127 (Befehl nicht gefunden) und 137 (Prozess beendet) liefern sofortige Diagnosehinweise auf die Art des Fehlers.

### Wie geht es weiter?

Lerne über Nextflows eingebaute Debugging-Tools und systematische Ansätze zur Fehlerbehebung.

---

## 4. Eingebaute Debugging-Tools und fortgeschrittene Techniken

Nextflow bietet mehrere leistungsstarke eingebaute Tools zum Debuggen und Analysieren der Workflow-Ausführung. Diese Tools helfen dir zu verstehen, was schiefgelaufen ist, wo es schiefgelaufen ist und wie du es effizient beheben kannst.

### 4.1. Echtzeit-Prozessausgabe

Manchmal musst du sehen, was in laufenden Prozessen passiert. Du kannst die Echtzeit-Prozessausgabe aktivieren, die dir genau zeigt, was jede Aufgabe während der Ausführung tut.

#### Pipeline ausführen

`bad_channel_shape_viewed.nf` aus unseren früheren Beispielen hat Kanalinhalt mit `.view()` ausgegeben, aber wir können auch die `debug`-Direktive verwenden, um Variablen aus dem Prozess selbst zu echoen, was wir in `bad_channel_shape_viewed_debug.nf` demonstrieren. Führe den Workflow aus:

```bash
nextflow run bad_channel_shape_viewed_debug.nf
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_channel_shape_viewed_debug.nf` [agitated_crick] DSL2 - revision: ea3676d9ec

    executor >  local (3)
    [c6/2dac51] process > PROCESS_FILES (3) [100%] 3 of 3 ✔
    Channel content: [sample1, file1.txt]
    Channel content: [sample2, file2.txt]
    Channel content: [sample3, file3.txt]
    After mapping: sample1
    After mapping: sample2
    After mapping: sample3
    Sample name inside process is sample2

    Sample name inside process is sample1

    Sample name inside process is sample3
    ```

#### Den Code prüfen

Schauen wir uns `bad_channel_shape_viewed_debug.nf` an, um zu sehen, wie die `debug`-Direktive funktioniert:

```groovy title="bad_channel_shape_viewed_debug.nf" linenums="3" hl_lines="2"
process PROCESS_FILES {
    debug true  // Echtzeit-Ausgabe aktivieren

    input:
    val sample_name

    output:
    path "${sample_name}_output.txt"

    script:
    """
    echo "Sample name inside process is ${sample_name}"
    echo "Processing ${sample_name}" > ${sample_name}_output.txt
    """
}
```

Die `debug`-Direktive kann eine schnelle und praktische Möglichkeit sein, die Umgebung eines Prozesses zu verstehen.

### 4.2. Preview-Modus

Manchmal möchtest du Probleme erkennen, bevor Prozesse ausgeführt werden. Nextflow bietet ein Flag für diese Art von proaktivem Debugging: `-preview`.

#### Pipeline ausführen

Der Preview-Modus ermöglicht es dir, die Workflow-Logik zu testen, ohne Befehle auszuführen. Das kann sehr nützlich sein, um schnell die Struktur deines Workflows zu überprüfen und sicherzustellen, dass Prozesse korrekt verbunden sind, ohne tatsächliche Befehle auszuführen.

!!! note "Hinweis"

    Wenn du `bad_syntax.nf` früher korrigiert hast, führe den Syntaxfehler wieder ein, indem du die schließende geschweifte Klammer nach dem script-Block entfernst, bevor du diesen Befehl ausführst.

Führe diesen Befehl aus:

```bash
nextflow run bad_syntax.nf -preview
```

??? failure "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_syntax.nf` [magical_mercator] DSL2 - revision: 550b9a8873

    Error bad_syntax.nf:24:1: Unexpected input: '<EOF>'

    ERROR ~ Script compilation failed

     -- Check '.nextflow.log' file for details
    ```

Der Preview-Modus ist besonders nützlich, um Syntaxfehler frühzeitig zu erkennen, ohne Prozesse auszuführen. Er validiert die Workflow-Struktur und Prozessverbindungen vor der Ausführung.

### 4.3. Stub-Ausführung zum Testen der Logik

Manchmal sind Fehler schwer zu debuggen, weil Befehle zu lange dauern, spezielle Software erfordern oder aus komplexen Gründen fehlschlagen. Die Stub-Ausführung ermöglicht es dir, die Workflow-Logik zu testen, ohne die eigentlichen Befehle auszuführen.

#### Pipeline ausführen

Wenn du einen Nextflow-Prozess entwickelst, kannst du die `stub`-Direktive verwenden, um „Dummy"-Befehle zu definieren, die Ausgaben der richtigen Form erzeugen, ohne den echten Befehl auszuführen. Dieser Ansatz ist besonders wertvoll, wenn du die Korrektheit deiner Workflow-Logik überprüfen möchtest, bevor du dich mit den Komplexitäten der eigentlichen Software befasst.

Erinnerst du dich an unser `missing_software.nf` von früher? Das, bei dem fehlende Software die Ausführung des Workflows verhinderte, bis wir `-profile docker` hinzufügten? `missing_software_with_stub.nf` ist ein sehr ähnlicher Workflow. Wenn wir ihn auf die gleiche Weise ausführen, erhalten wir denselben Fehler:

```bash
nextflow run missing_software_with_stub.nf
```

??? failure "Befehlsausgabe"

    ```console hl_lines="12 18"
    ERROR ~ Error executing process > 'PROCESS_FILES (3)'

    Caused by:
      Process `PROCESS_FILES (3)` terminated with an error exit status (127)


    Command executed:

      cowpy sample3 > sample3_output.txt

    Command exit status:
      127

    Command output:
      (empty)

    Command error:
      .command.sh: line 2: cowpy: command not found

    Work dir:
      /workspaces/training/side-quests/debugging/work/82/42a5bfb60c9c6ee63ebdbc2d51aa6e

    Tip: you can try to figure out what's wrong by changing to the process work directory and showing the script file named `.command.sh`

    -- Check '.nextflow.log' file for details
    ```

Dieser Workflow erzeugt jedoch keine Fehler, wenn wir ihn mit `-stub-run` ausführen, auch ohne das `docker`-Profil:

```bash
nextflow run missing_software_with_stub.nf -stub-run
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `missing_software_with_stub.nf` [astonishing_shockley] DSL2 - revision: f1f4f05d7d

    executor >  local (3)
    [b5/2517a3] PROCESS_FILES (3) | 3 of 3 ✔
    ```

#### Den Code prüfen

Schauen wir uns `missing_software_with_stub.nf` an:

```groovy title="missing_software.nf (with stub)" hl_lines="16-19" linenums="3"
process PROCESS_FILES {

    container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

    input:
    val sample_name

    output:
    path "${sample_name}_output.txt"

    script:
    """
    cowpy ${sample_name} > ${sample_name}_output.txt
    """

    stub:
    """
    touch ${sample_name}_output.txt
    """
}
```

Im Vergleich zu `missing_software.nf` hat dieser Prozess eine `stub:`-Direktive, die einen Befehl angibt, der anstelle des in `script:` angegebenen verwendet wird, wenn Nextflow im Stub-Modus ausgeführt wird.

Der `touch`-Befehl, den wir hier verwenden, hängt von keiner Software oder geeigneten Eingaben ab und läuft in allen Situationen, sodass wir die Workflow-Logik debuggen können, ohne uns um die Prozess-Interna zu kümmern.

**Die Stub-Ausführung hilft beim Debuggen von:**

- Kanalstruktur und Datenfluss
- Prozessverbindungen und -abhängigkeiten
- Parameter-Weitergabe
- Workflow-Logik ohne Software-Abhängigkeiten

### 4.4. Systematischer Debugging-Ansatz

Jetzt, wo du einzelne Debugging-Techniken gelernt hast – von Trace-Dateien und work directories bis hin zu Preview-Modus, Stub-Ausführung und Ressourcenüberwachung – lass uns diese zu einer systematischen Methodik zusammenfassen. Ein strukturierter Ansatz verhindert, dass du von komplexen Fehlern überwältigt wirst, und stellt sicher, dass du keine wichtigen Hinweise übersiehst.

Diese Methodik kombiniert alle behandelten Tools zu einem effizienten Workflow:

**Vier-Phasen-Debugging-Methode:**

**Phase 1: Syntaxfehler beheben (5 Minuten)**

1. Auf rote Unterstreichungen in VSCode oder deiner IDE prüfen
2. `nextflow run workflow.nf -preview` ausführen, um Syntaxprobleme zu identifizieren
3. Alle Syntaxfehler beheben (fehlende geschweifte Klammern, abschließende Kommas usw.)
4. Sicherstellen, dass der Workflow erfolgreich geparst wird, bevor du weitermachst

**Phase 2: Schnelle Bewertung (5 Minuten)**

1. Laufzeit-Fehlermeldungen sorgfältig lesen
2. Prüfen, ob es sich um einen Laufzeit-, Logik- oder Ressourcenfehler handelt
3. Preview-Modus verwenden, um grundlegende Workflow-Logik zu testen

**Phase 3: Detaillierte Untersuchung (15–30 Minuten)**

1. Das work directory der fehlgeschlagenen Aufgabe finden
2. Log-Dateien untersuchen
3. `.view()`-Operatoren hinzufügen, um Kanäle zu inspizieren
4. `-stub-run` verwenden, um Workflow-Logik ohne Ausführung zu testen

**Phase 4: Korrigieren und validieren (15 Minuten)**

1. Minimale, gezielte Korrekturen vornehmen
2. Mit resume testen: `nextflow run workflow.nf -resume`
3. Vollständige Workflow-Ausführung überprüfen

!!! tip "Resume für effizientes Debugging nutzen"

    Sobald du ein Problem identifiziert hast, brauchst du eine effiziente Möglichkeit, deine Korrekturen zu testen, ohne Zeit damit zu verschwenden, erfolgreiche Teile deines Workflows erneut auszuführen. Nextflows `-resume`-Funktionalität ist beim Debugging unverzichtbar.

    Du bist `-resume` begegnet, wenn du [Hello Nextflow](../hello_nextflow/) durchgearbeitet hast, und es ist wichtig, dass du es beim Debugging gut nutzt, um nicht warten zu müssen, während die Prozesse vor deinem Problemprozess erneut ausgeführt werden.

    **Resume-Debugging-Strategie:**

    1. Workflow bis zum Fehler ausführen
    2. Work directory der fehlgeschlagenen Aufgabe untersuchen
    3. Das spezifische Problem beheben
    4. Mit resume fortfahren, um nur die Korrektur zu testen
    5. Wiederholen, bis der Workflow abgeschlossen ist

#### Debugging-Konfigurationsprofil

Um diesen systematischen Ansatz noch effizienter zu machen, kannst du eine dedizierte Debugging-Konfiguration erstellen, die automatisch alle benötigten Tools aktiviert:

```groovy title="nextflow.config (debug profile)" linenums="1"
profiles {
    debug {
        process {
            debug = true
            cleanup = false

            // Konservative Ressourcen für das Debugging
            maxForks = 1
            memory = '2.GB'
            cpus = 1
        }
    }
}
```

Dann kannst du die Pipeline mit diesem aktivierten Profil ausführen:

```bash
nextflow run workflow.nf -profile debug
```

Dieses Profil aktiviert die Echtzeit-Ausgabe, bewahrt work directories und begrenzt die Parallelisierung für einfacheres Debugging.

### 4.5. Praktische Debugging-Übung

Jetzt ist es Zeit, den systematischen Debugging-Ansatz in die Praxis umzusetzen. Der Workflow `buggy_workflow.nf` enthält mehrere häufige Fehler, die die Arten von Problemen repräsentieren, denen du in der realen Entwicklung begegnen wirst.

!!! exercise "Übung"

    Verwende den systematischen Debugging-Ansatz, um alle Fehler in `buggy_workflow.nf` zu identifizieren und zu beheben. Dieser Workflow versucht, Probendaten aus einer CSV-Datei zu verarbeiten, enthält aber mehrere absichtliche Bugs, die häufige Debugging-Szenarien repräsentieren.

    Beginne damit, den Workflow auszuführen, um den ersten Fehler zu sehen:

    ```bash
    nextflow run buggy_workflow.nf
    ```

    ??? failure "Befehlsausgabe"

        ```console
        N E X T F L O W   ~  version 25.10.2

        Launching `buggy_workflow.nf` [wise_ramanujan] DSL2 - revision: d51a8e83fd

        ERROR ~ Range [11, 12) out of bounds for length 11

         -- Check '.nextflow.log' file for details
        ```

        Dieser kryptische Fehler weist auf ein Parsing-Problem um Zeile 11–12 im `params{}`-Block hin. Der v2-Parser erkennt strukturelle Probleme frühzeitig.

    Wende die Vier-Phasen-Debugging-Methode an, die du gelernt hast:

    **Phase 1: Syntaxfehler beheben**
    - Auf rote Unterstreichungen in VSCode oder deiner IDE prüfen
    - `nextflow run workflow.nf -preview` ausführen, um Syntaxprobleme zu identifizieren
    - Alle Syntaxfehler beheben (fehlende geschweifte Klammern, abschließende Kommas usw.)
    - Sicherstellen, dass der Workflow erfolgreich geparst wird, bevor du weitermachst

    **Phase 2: Schnelle Bewertung**
    - Laufzeit-Fehlermeldungen sorgfältig lesen
    - Identifizieren, ob Fehler laufzeitbezogen, logikbezogen oder ressourcenbezogen sind
    - `-preview`-Modus verwenden, um grundlegende Workflow-Logik zu testen

    **Phase 3: Detaillierte Untersuchung**
    - Work directories für fehlgeschlagene Aufgaben untersuchen
    - `.view()`-Operatoren hinzufügen, um Kanäle zu inspizieren
    - Log-Dateien in work directories prüfen
    - `-stub-run` verwenden, um Workflow-Logik ohne Ausführung zu testen

    **Phase 4: Korrigieren und validieren**
    - Gezielte Korrekturen vornehmen
    - `-resume` verwenden, um Korrekturen effizient zu testen
    - Vollständige Workflow-Ausführung überprüfen

    **Verfügbare Debugging-Tools:**
    ```bash
    # Preview-Modus für Syntaxprüfung
    nextflow run buggy_workflow.nf -preview

    # Debug-Profil für detaillierte Ausgabe
    nextflow run buggy_workflow.nf -profile debug

    # Stub-Ausführung zum Testen der Logik
    nextflow run buggy_workflow.nf -stub-run

    # Resume nach Korrekturen
    nextflow run buggy_workflow.nf -resume
    ```

    ??? solution "Lösung"
        Der `buggy_workflow.nf` enthält 9 oder 10 verschiedene Fehler (je nachdem, wie man zählt), die alle wichtigen Debugging-Kategorien abdecken. Hier ist eine systematische Aufschlüsselung jedes Fehlers und wie man ihn behebt.

        Beginnen wir mit den Syntaxfehlern:

        **Fehler 1: Syntaxfehler – Abschließendes Komma**
        ```groovy linenums="21"
        output:
            path "${sample_id}_result.txt",  // FEHLER: Abschließendes Komma
        ```
        **Korrektur:** Abschließendes Komma entfernen
        ```groovy linenums="21"
        output:
            path "${sample_id}_result.txt"
        ```

        **Fehler 2: Syntaxfehler – Fehlende schließende geschweifte Klammer**
        ```groovy linenums="24"
        script:
        """
        echo "Processing: ${sample}"
        cat ${input_file} > ${sample}_result.txt
        """
        // FEHLER: Fehlende schließende geschweifte Klammer für den processFiles-Prozess
        ```
        **Korrektur:** Fehlende schließende geschweifte Klammer hinzufügen
        ```groovy linenums="29"
        """
        echo "Processing: ${sample_id}"
        cat ${input_file} > ${sample_id}_result.txt
        """
        }  // Fehlende schließende geschweifte Klammer hinzufügen
        ```

        **Fehler 3: Variablennamenfehler**
        ```groovy linenums="26"
        echo "Processing: ${sample}"     // FEHLER: Sollte sample_id sein
        cat ${input_file} > ${sample}_result.txt  // FEHLER: Sollte sample_id sein
        ```
        **Korrektur:** Korrekten Eingabevariablennamen verwenden
        ```groovy linenums="26"
        echo "Processing: ${sample_id}"
        cat ${input_file} > ${sample_id}_result.txt
        ```

        **Fehler 4: Undefinierter Variablenfehler**
        ```groovy linenums="87"
        heavy_ch = heavyProcess(sample_ids)  // FEHLER: sample_ids ist undefiniert
        ```
        **Korrektur:** Korrekten Kanal verwenden und Proben-IDs extrahieren
        ```groovy linenums="87"
        heavy_ch = heavyProcess(input_ch)
        ```

        An diesem Punkt wird der Workflow ausgeführt, aber wir erhalten noch Fehler (z. B. `Path value cannot be null` in `processFiles`), verursacht durch eine fehlerhafte Kanalstruktur.

        **Fehler 5: Kanalstrukturfehler – Falsche map-Ausgabe**
        ```groovy linenums="83"
        .map { row -> row.sample_id }  // FEHLER: processFiles erwartet Tupel
        ```
        **Korrektur:** Die Tupelstruktur zurückgeben, die processFiles erwartet
        ```groovy linenums="83"
        .map { row -> [row.sample_id, file(row.fastq_path)] }
        ```

        Das bricht aber unsere Lösung für `heavyProcess()` oben, also müssen wir eine map verwenden, um nur die Proben-IDs an diesen Prozess zu übergeben:

        **Fehler 6: Fehlerhafte Kanalstruktur für heavyProcess**
        ```groovy linenums="87"
        heavy_ch = heavyProcess(input_ch)  // FEHLER: input_ch hat jetzt 2 Elemente pro Emission – heavyProcess braucht nur 1 (das erste)
        ```
        **Korrektur:** Korrekten Kanal verwenden und Proben-IDs extrahieren
        ```groovy linenums="87"
        heavy_ch = heavyProcess(input_ch.map{it[0]})
        ```

        Jetzt kommen wir weiter, erhalten aber einen Fehler über `No such variable: i`, weil wir eine Bash-Variable nicht escaped haben.

        **Fehler 7: Bash-Variablen-Escaping-Fehler**
        ```groovy linenums="48"
        echo "Heavy computation $i for ${sample_id}"  // FEHLER: $i nicht escaped
        ```
        **Korrektur:** Bash-Variable escapen
        ```groovy linenums="48"
        echo "Heavy computation \${i} for ${sample_id}"
        ```

        Jetzt erhalten wir `Process exceeded running time limit (1ms)`, also korrigieren wir das Zeitlimit für den betreffenden Prozess:

        **Fehler 8: Ressourcenkonfigurationsfehler**
        ```groovy linenums="36"
        time '1 ms'  // FEHLER: Unrealistisches Zeitlimit
        ```
        **Korrektur:** Auf ein realistisches Zeitlimit erhöhen
        ```groovy linenums="36"
        time '100 s'
        ```

        Als nächstes haben wir einen `Missing output file(s)`-Fehler zu beheben:

        **Fehler 9: Ausgabedateiname stimmt nicht überein**
        ```groovy linenums="49"
        done > ${sample_id}.txt  // FEHLER: Falscher Dateiname, sollte mit der output-Deklaration übereinstimmen
        ```
        **Korrektur:** Mit der output-Deklaration übereinstimmen
        ```groovy linenums="49"
        done > ${sample_id}_heavy.txt
        ```

        Die ersten beiden Prozesse liefen, aber nicht der dritte.

        **Fehler 10: Ausgabedateiname stimmt nicht überein**
        ```groovy linenums="88"
        file_ch = channel.fromPath("*.txt") // FEHLER: Versucht, Eingabe aus dem aktuellen Verzeichnis zu nehmen, statt aus einem Prozess
        handleFiles(file_ch)
        ```
        **Korrektur:** Die Ausgabe des vorherigen Prozesses verwenden
        ```groovy linenums="88"
        handleFiles(heavyProcess.out)
        ```

        Damit sollte der gesamte Workflow ausgeführt werden.

        **Vollständig korrigierter Workflow:**
        ```groovy linenums="1"
        #!/usr/bin/env nextflow

        /*
        * Fehlerhafter Workflow für Debugging-Übungen
        * Dieser Workflow enthält mehrere absichtliche Bugs zu Lernzwecken
        */

        params{
            // Parameter ohne Validierung
            input: Path = 'data/sample_data.csv'
            output: String = 'results'
        }

        /*
        * Prozess mit Eingabe-/Ausgabe-Diskrepanz
        */
        process processFiles {
            publishDir "${params.output}/processed", mode: 'copy'

            input:
                tuple val(sample_id), path(input_file)

            output:
                path "${sample_id}_result.txt"

            script:
            """
            echo "Processing: ${sample_id}"
            cat ${input_file} > ${sample_id}_result.txt
            """
        }

        /*
        * Prozess mit Ressourcenproblemen
        */
        process heavyProcess {
            publishDir "${params.output}/heavy", mode: 'copy'

            time '100 s'

            input:
                val sample_id

            output:
                path "${sample_id}_heavy.txt"

            script:
            """
            # Schwere Berechnung simulieren
            for i in {1..1000000}; do
                echo "Heavy computation \$i for ${sample_id}"
            done > ${sample_id}_heavy.txt
            """
        }

        /*
        * Prozess mit Dateiverarbeitungsproblemen
        */
        process handleFiles {
            publishDir "${params.output}/files", mode: 'copy'

            input:
                path input_file

            output:
                path "processed_${input_file}"

            script:
            """
            if [ -f "${input_file}" ]; then
                cp ${input_file} processed_${input_file}
            fi
            """
        }

        /*
        * Haupt-Workflow mit Kanalproblemen
        */
        workflow {

            // Kanal mit falscher Verwendung
            input_ch = channel
                .fromPath(params.input)
                .splitCsv(header: true)
                .map { row -> [row.sample_id, file(row.fastq_path)] }

            processed_ch = processFiles(input_ch)

            heavy_ch = heavyProcess(input_ch.map{it[0]})

            handleFiles(heavyProcess.out)
        }
        ```

**Abgedeckte Fehlerkategorien:**

- **Syntaxfehler**: Fehlende geschweifte Klammern, abschließende Kommas, undefinierte Variablen
- **Kanalstrukturfehler**: Falsche Datenformen, undefinierte Kanäle
- **Prozessfehler**: Ausgabedatei-Diskrepanzen, Variablen-Escaping
- **Ressourcenfehler**: Unrealistische Zeitlimits

**Wichtige Debugging-Lektionen:**

1. **Fehlermeldungen sorgfältig lesen** – sie zeigen oft direkt auf das Problem
2. **Systematische Ansätze verwenden** – einen Fehler nach dem anderen beheben und mit `-resume` testen
3. **Datenfluss verstehen** – Kanalstrukturfehler sind oft die subtilsten
4. **Work directories prüfen** – wenn Prozesse fehlschlagen, sagen dir die Logs genau, was schiefgelaufen ist

---

## Zusammenfassung

In dieser Side Quest hast du eine Reihe systematischer Techniken zum Debuggen von Nextflow-Workflows gelernt.
Wenn du diese Techniken in deiner eigenen Arbeit anwendest, wirst du weniger Zeit damit verbringen, mit deinem Computer zu kämpfen, Probleme schneller lösen und dich vor zukünftigen Problemen schützen.

### Wichtige Muster

**1. Syntaxfehler identifizieren und beheben**:

- Nextflow-Fehlermeldungen interpretieren und Probleme lokalisieren
- Häufige Syntaxfehler: fehlende geschweifte Klammern, falsche Keywords, undefinierte Variablen
- Unterschied zwischen Nextflow-(Groovy-)Variablen und Bash-Variablen
- VS Code-Erweiterungsfunktionen für frühzeitige Fehlererkennung nutzen

```groovy
// Fehlende geschweifte Klammer – auf rote Unterstreichungen in der IDE achten
process FOO {
    script:
    """
    echo "hello"
    """
// } <-- fehlt!

// Falsches Keyword
inputs:  // Sollte 'input:' sein

// Undefinierte Variable – mit Backslash für Bash-Variablen escapen
echo "${undefined_var}"      // Nextflow-Variable (Fehler, wenn nicht definiert)
echo "\${bash_var}"          // Bash-Variable (escaped)
```

**2. Kanalstrukturprobleme debuggen**:

- Kanalkardinaliät und Erschöpfungsprobleme verstehen
- Diskrepanzen in der Kanalinhaltsstruktur debuggen
- `.view()`-Operatoren zur Kanalinspektion verwenden
- Fehlermuster wie eckige Klammern in der Ausgabe erkennen

```groovy
// Kanalinhalt inspizieren
my_channel.view { "Content: $it" }

// Queue channel in value channel umwandeln (verhindert Erschöpfung)
reference_ch = channel.value('ref.fa')
// oder
reference_ch = channel.of('ref.fa').first()
```

**3. Prozessausführungsprobleme beheben**:

- Fehlende Ausgabedatei-Fehler diagnostizieren
- Exit-Codes verstehen (127 für fehlende Software, 137 für Speicherprobleme)
- Work directories und Befehlsdateien untersuchen
- Ressourcen angemessen konfigurieren

```bash
# Prüfen, was tatsächlich ausgeführt wurde
cat work/ab/cdef12/.command.sh

# Fehlerausgabe prüfen
cat work/ab/cdef12/.command.err

# Exit-Code 127 = Befehl nicht gefunden
# Exit-Code 137 = beendet (Speicher-/Zeitlimit)
```

**4. Nextflows eingebaute Debugging-Tools nutzen**:

- Preview-Modus und Echtzeit-Debugging nutzen
- Stub-Ausführung für Logiktests implementieren
- Resume für effiziente Debugging-Zyklen anwenden
- Eine systematische Vier-Phasen-Debugging-Methodik befolgen

!!! tip "Schnelle Debugging-Referenz"

    **Syntaxfehler?** → VSCode-Warnungen prüfen, `nextflow run workflow.nf -preview` ausführen

    **Kanalprobleme?** → `.view()` verwenden, um Inhalt zu inspizieren: `my_channel.view()`

    **Prozessfehler?** → Work directory-Dateien prüfen:

    - `.command.sh` – das ausgeführte Skript
    - `.command.err` – Fehlermeldungen
    - `.exitcode` – Exit-Status (127 = Befehl nicht gefunden, 137 = beendet)

    **Mysteriöses Verhalten?** → Mit `-stub-run` ausführen, um Workflow-Logik zu testen

    **Korrekturen vorgenommen?** → `-resume` verwenden, um Zeit beim Testen zu sparen: `nextflow run workflow.nf -resume`

---

### Weitere Ressourcen

- [Nextflow-Fehlerbehebungshandbuch](https://www.nextflow.io/docs/latest/troubleshooting.html): Offizielle Fehlerbehebungsdokumentation
- [Nextflow-Kanäle verstehen](https://www.nextflow.io/docs/latest/channel.html): Tiefgehende Einführung in Kanaltypen und -verhalten
- [Prozess-Direktiven-Referenz](https://www.nextflow.io/docs/latest/process.html#directives): Alle verfügbaren Prozesskonfigurationsoptionen
- [nf-test](https://www.nf-test.com/): Test-Framework für Nextflow-Pipelines
- [Nextflow Slack-Community](https://www.nextflow.io/slack-invite.html): Hilfe von der Community erhalten

Für Produktions-Workflows solltest du Folgendes in Betracht ziehen:

- [Seqera Platform](https://seqera.io/platform/) für Monitoring und Debugging im großen Maßstab einrichten
- [Wave-Container](https://seqera.io/wave/) für reproduzierbare Software-Umgebungen verwenden

**Denk daran:** Effektives Debugging ist eine Fähigkeit, die sich mit der Praxis verbessert. Die systematische Methodik und das umfassende Toolkit, das du hier erworben hast, werden dir auf deiner gesamten Nextflow-Entwicklungsreise gute Dienste leisten.

---

## Wie geht es weiter?

Kehre zum [Menü der Side Quests](../) zurück oder klicke auf die Schaltfläche unten rechts auf der Seite, um zum nächsten Thema in der Liste zu wechseln.
