# Debugging von Workflows

Debugging ist eine wichtige Fähigkeit, die dir Stunden der Frustration ersparen und dich zu einem effektiveren Nextflow-Entwickler machen kann. Im Laufe deiner Karriere, besonders wenn du gerade anfängst, wirst du beim Erstellen und Pflegen deiner Workflows auf Fehler stoßen. Systematische Debugging-Ansätze zu lernen hilft dir, Probleme schnell zu identifizieren und zu beheben.

### Lernziele

In dieser Side Quest erkunden wir **systematische Debugging-Techniken** für Nextflow-Workflows:

- **Debugging von Syntaxfehlern**: Effektive Nutzung von IDE-Features und Nextflow-Fehlermeldungen
- **Kanal-Debugging**: Diagnose von Datenfluss-Problemen und Kanal-Strukturproblemen
- **Prozess-Debugging**: Untersuchung von Ausführungsfehlern und Ressourcenproblemen
- **Eingebaute Debugging-Tools**: Nutzung von Nextflows Preview-Modus, Stub-Running und Work-Verzeichnissen
- **Systematische Ansätze**: Eine vierphasige Methodik für effizientes Debugging

Am Ende wirst du eine robuste Debugging-Methodik haben, die frustrierende Fehlermeldungen in klare Lösungswege verwandelt.

### Voraussetzungen

Bevor du diese Side Quest beginnst, solltest du:

- Das [Hello Nextflow](../hello_nextflow/README.md) Tutorial oder einen gleichwertigen Einsteigerkurs abgeschlossen haben.
- Dich mit grundlegenden Nextflow-Konzepten und -Mechanismen (Prozesse, Kanäle, Operatoren) wohlfühlen

**Optional:** Wir empfehlen, zuerst die Side Quest [IDE Features for Nextflow Development](./ide_features.md) abzuschließen.
Diese bietet umfassende Informationen zu IDE-Features, die das Debugging unterstützen (Syntax-Highlighting, Fehlererkennung usw.), die wir hier intensiv nutzen werden.

---

## 0. Erste Schritte

#### Öffne den Training-Codespace

Falls noch nicht geschehen, öffne die Trainingsumgebung wie in [Umgebung einrichten](../envsetup/index.md) beschrieben.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Wechsle in das Projektverzeichnis

Wechseln wir in das Verzeichnis, in dem sich die Dateien für dieses Tutorial befinden.

```bash
cd side-quests/debugging
```

Du kannst VSCode auf dieses Verzeichnis fokussieren:

```bash
code .
```

#### Überprüfe die Materialien

Du findest eine Reihe von Beispiel-Workflows mit verschiedenen Arten von Fehlern, die wir zum Üben verwenden werden:

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

#### Überprüfe die Aufgabe

Deine Herausforderung ist es, jeden Workflow auszuführen, die Fehler zu identifizieren und zu beheben.

Für jeden fehlerhaften Workflow:

1. **Führe den Workflow aus** und beobachte den Fehler
2. **Analysiere die Fehlermeldung**: Was sagt dir Nextflow?
3. **Lokalisiere das Problem** im Code anhand der bereitgestellten Hinweise
4. **Behebe den Fehler** und überprüfe, ob deine Lösung funktioniert
5. **Setze die Datei zurück**, bevor du zum nächsten Abschnitt übergehst (verwende `git checkout <dateiname>`)

Die Übungen gehen von einfachen Syntaxfehlern zu subtileren Laufzeitproblemen über.
Lösungen werden inline diskutiert, aber versuche, jede Übung selbst zu lösen, bevor du weiterliest.

#### Bereitschafts-Checkliste

Denkst du, du bist bereit loszulegen?

- [ ] Ich verstehe das Ziel dieses Kurses und seine Voraussetzungen
- [ ] Mein Codespace läuft
- [ ] Ich habe mein Arbeitsverzeichnis entsprechend gesetzt
- [ ] Ich verstehe die Aufgabe

Wenn du alle Kästchen abhaken kannst, kann es losgehen.

---

## 1. Syntaxfehler

Syntaxfehler sind die häufigste Art von Fehlern, denen du beim Schreiben von Nextflow-Code begegnen wirst. Sie treten auf, wenn der Code nicht den erwarteten Syntaxregeln der Nextflow-DSL entspricht. Diese Fehler verhindern, dass dein Workflow überhaupt läuft, daher ist es wichtig zu lernen, wie man sie schnell identifiziert und behebt.

### 1.1. Fehlende Klammern

Einer der häufigsten Syntaxfehler, und manchmal einer der komplexeren zu debuggen, sind **fehlende oder nicht übereinstimmende Klammern**.

Beginnen wir mit einem praktischen Beispiel.

#### Führe die Pipeline aus

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

**Schlüsselelemente von Syntaxfehler-Meldungen:**

- **Datei und Position**: Zeigt, welche Datei und Zeile/Spalte den Fehler enthält (`bad_syntax.nf:24:1`)
- **Fehlerbeschreibung**: Erklärt, was der Parser gefunden hat, das er nicht erwartet hat (`Unexpected input: '<EOF>'`)
- **EOF-Indikator**: Die `<EOF>` (End Of File) Meldung zeigt an, dass der Parser das Ende der Datei erreicht hat, während er noch mehr Inhalt erwartet - ein klassisches Zeichen für nicht geschlossene Klammern

#### Überprüfe den Code

Schauen wir uns nun `bad_syntax.nf` an, um zu verstehen, was den Fehler verursacht:

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
// Missing closing brace for the process

workflow {

    // Create input channel
    input_ch = channel.of('sample1', 'sample2', 'sample3')

    // Call the process with the input channel
    PROCESS_FILES(input_ch)
}
```

Für dieses Beispiel haben wir einen Kommentar hinterlassen, der zeigt, wo der Fehler ist. Die Nextflow VSCode-Extension sollte dir auch einige Hinweise geben, was falsch sein könnte, indem sie die nicht übereinstimmende Klammer rot markiert und das vorzeitige Ende der Datei hervorhebt:

![Bad syntax](img/bad_syntax.png)

**Debugging-Strategie für Klammerfehler:**

1. Verwende VS Codes Klammer-Matching (platziere den Cursor neben einer Klammer)
2. Überprüfe das Problems-Panel auf klammerbezogene Meldungen
3. Stelle sicher, dass jede öffnende `{` eine entsprechende schließende `}` hat

#### Behebe den Code

Ersetze den Kommentar durch die fehlende schließende Klammer:

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
    }  // Add the missing closing brace

    workflow {

        // Create input channel
        input_ch = channel.of('sample1', 'sample2', 'sample3')

        // Call the process with the input channel
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
    // Missing closing brace for the process

    workflow {

        // Create input channel
        input_ch = channel.of('sample1', 'sample2', 'sample3')

        // Call the process with the input channel
        PROCESS_FILES(input_ch)
    }
    ```

#### Führe die Pipeline aus

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

### 1.2. Verwendung falscher Prozess-Schlüsselwörter oder Direktiven

Ein weiterer häufiger Syntaxfehler ist eine **ungültige Prozessdefinition**. Dies kann passieren, wenn du vergisst, erforderliche Blöcke zu definieren oder falsche Direktiven in der Prozessdefinition verwendest.

#### Führe die Pipeline aus

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

#### Überprüfe den Code

Der Fehler zeigt eine "Invalid process definition" an und zeigt den Kontext um das Problem. Wenn wir uns die Zeilen 3-7 ansehen, können wir `inputs:` in Zeile 4 sehen, was das Problem ist. Schauen wir uns `invalid_process.nf` an:

```groovy title="invalid_process.nf" hl_lines="4" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    inputs:  // ERROR: Should be 'input' not 'inputs'
    val sample_name

    output:
    path "${sample_name}_output.txt"

    script:
    """
    echo "Processing ${sample_name}" > ${sample_name}_output.txt
    """
}

workflow {

    // Create input channel
    input_ch = channel.of('sample1', 'sample2', 'sample3')

    // Call the process with the input channel
    PROCESS_FILES(input_ch)
}
```

Wenn wir uns Zeile 4 im Fehlerkontext ansehen, können wir das Problem erkennen: Wir verwenden `inputs` statt der korrekten `input`-Direktive. Die Nextflow VSCode-Extension wird dies ebenfalls markieren:

![Invalid process message](img/invalid_process_message.png)

#### Behebe den Code

Ersetze das falsche Schlüsselwort durch das korrekte, indem du [die Dokumentation](https://www.nextflow.io/docs/latest/process.html#) konsultierst:

=== "Danach"

    ```groovy title="invalid_process.nf" hl_lines="4" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:  // Fixed: Changed 'inputs' to 'input'
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        """
        echo "Processing ${sample_name}" > ${sample_name}_output.txt
        """
    }

    workflow {

        // Create input channel
        input_ch = channel.of('sample1', 'sample2', 'sample3')

        // Call the process with the input channel
        PROCESS_FILES(input_ch)
    }
    ```

=== "Vorher"

    ```groovy title="invalid_process.nf" hl_lines="4" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        inputs:  // ERROR: Should be 'input' not 'inputs'
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        """
        echo "Processing ${sample_name}" > ${sample_name}_output.txt
        """
    }

    workflow {

        // Create input channel
        input_ch = channel.of('sample1', 'sample2', 'sample3')

        // Call the process with the input channel
        PROCESS_FILES(input_ch)
    }
    ```

#### Führe die Pipeline aus

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

### 1.3. Verwendung ungültiger Variablennamen

Die Variablennamen, die du in deinen Script-Blöcken verwendest, müssen gültig sein und entweder aus Eingaben oder aus Groovy-Code stammen, der vor dem Script eingefügt wurde. Aber wenn du zu Beginn der Pipeline-Entwicklung mit Komplexität jonglierst, ist es leicht, Fehler bei der Variablenbenennung zu machen, und Nextflow wird dich schnell darauf hinweisen.

#### Führe die Pipeline aus

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

Der Fehler wird zur Kompilierzeit erkannt und zeigt direkt auf die undefinierte Variable in Zeile 17, mit einem Caret, das genau anzeigt, wo das Problem ist.

#### Überprüfe den Code

Schauen wir uns `no_such_var.nf` an:

```groovy title="no_such_var.nf" hl_lines="17" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    input:
    val sample_name

    output:
    path "${sample_name}_processed.txt"

    script:
    // Define variables in Groovy code before the script
    def output_prefix = "${sample_name}_processed"
    def timestamp = new Date().format("yyyy-MM-dd")

    """
    echo "Processing ${sample_name} on ${timestamp}" > ${output_prefix}.txt
    echo "Using undefined variable: ${undefined_var}" >> ${output_prefix}.txt  // ERROR: undefined_var not defined
    """
}

workflow {
    input_ch = channel.of('sample1', 'sample2', 'sample3')
    PROCESS_FILES(input_ch)
}
```

Die Fehlermeldung zeigt an, dass die Variable im Script-Template nicht erkannt wird, und da ist es - du solltest `${undefined_var}` im Script-Block verwendet sehen, aber nirgendwo anders definiert.

#### Behebe den Code

Wenn du einen 'No such variable'-Fehler erhältst, kannst du ihn beheben, indem du entweder die Variable definierst (durch Korrektur von Eingabevariablennamen oder Bearbeitung von Groovy-Code vor dem Script) oder sie aus dem Script-Block entfernst, wenn sie nicht benötigt wird:

=== "Danach"

    ```groovy title="no_such_var.nf" hl_lines="15-17" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        // Define variables in Groovy code before the script
        def output_prefix = "${sample_name}_processed"
        def timestamp = new Date().format("yyyy-MM-dd")

        """
        echo "Processing ${sample_name} on ${timestamp}" > ${output_prefix}.txt
        """  // Removed the line with undefined_var
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
        // Define variables in Groovy code before the script
        def output_prefix = "${sample_name}_processed"
        def timestamp = new Date().format("yyyy-MM-dd")

        """
        echo "Processing ${sample_name} on ${timestamp}" > ${output_prefix}.txt
        echo "Using undefined variable: ${undefined_var}" >> ${output_prefix}.txt  // ERROR: undefined_var not defined
        """
    }

    workflow {
        input_ch = channel.of('sample1', 'sample2', 'sample3')
        PROCESS_FILES(input_ch)
    }
    ```

#### Führe die Pipeline aus

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

Zu Beginn mit Nextflow kann es schwierig sein, den Unterschied zwischen Nextflow- (Groovy-) und Bash-Variablen zu verstehen. Dies kann eine weitere Form des Fehlers mit ungültigen Variablen erzeugen, der auftritt, wenn man versucht, Variablen im Bash-Inhalt des Script-Blocks zu verwenden.

#### Führe die Pipeline aus

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

#### Überprüfe den Code

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

In diesem Beispiel definieren wir die `prefix`-Variable in Bash, aber in einem Nextflow-Prozess wird die `$`-Syntax, die wir verwendet haben, um darauf zu verweisen (`${prefix}`), als Groovy-Variable interpretiert, nicht als Bash. Die Variable existiert nicht im Groovy-Kontext, daher erhalten wir einen 'no such variable'-Fehler.

#### Behebe den Code

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
        echo "Processing ${sample_name}" > \${prefix}.txt  # Behoben: Dollarzeichen escaped
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

Dies sagt Nextflow, dies als Bash-Variable zu interpretieren.

#### Führe die Pipeline aus

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

!!! tip "Groovy vs Bash-Variablen"

    Für einfache Variablenmanipulationen wie String-Verkettung oder Präfix-/Suffix-Operationen ist es normalerweise lesbarer, Groovy-Variablen im Script-Abschnitt zu verwenden statt Bash-Variablen im Script-Block:

    ```groovy linenums="1"
    script:
    def output_prefix = "${sample_name}_processed"
    def output_file = "${output_prefix}.txt"
    """
    echo "Processing ${sample_name}" > ${output_file}
    """
    ```

    Dieser Ansatz vermeidet die Notwendigkeit, Dollarzeichen zu escapen, und macht den Code einfacher zu lesen und zu warten.

### 1.5. Anweisungen außerhalb des Workflow-Blocks

Die Nextflow VSCode-Extension hebt Probleme mit der Code-Struktur hervor, die Fehler verursachen werden. Ein häufiges Beispiel ist die Definition von Kanälen außerhalb des `workflow {}`-Blocks - dies wird jetzt als Syntaxfehler erzwungen.

#### Führe die Pipeline aus

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

Die Fehlermeldung zeigt das Problem klar an: Anweisungen (wie Kanaldefinitionen) können nicht mit Script-Deklarationen außerhalb eines Workflow- oder Prozess-Blocks gemischt werden.

#### Überprüfe den Code

Schauen wir uns `badpractice_syntax.nf` an, um zu sehen, was den Fehler verursacht:

```groovy title="badpractice_syntax.nf" hl_lines="3" linenums="1"
#!/usr/bin/env nextflow

input_ch = channel.of('sample1', 'sample2', 'sample3')  // ERROR: Channel defined outside workflow

process PROCESS_FILES {
    input:
    val sample_name

    output:
    path "${sample_name}_processed.txt"

    script:
    // Define variables in Groovy code before the script
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

Die VSCode-Extension wird auch die `input_ch`-Variable hervorheben, die außerhalb des Workflow-Blocks definiert ist:

![Non-lethal syntax error](img/nonlethal.png)

#### Behebe den Code

Verschiebe die Kanaldefinition in den Workflow-Block:

=== "Danach"

    ```groovy title="badpractice_syntax.nf" hl_lines="21" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_processed.txt"

        script:
        // Define variables in Groovy code before the script
        def output_prefix = "${sample_name}_processed"
        def timestamp = new Date().format("yyyy-MM-dd")

        """
        echo "Processing ${sample_name} on ${timestamp}" > ${output_prefix}.txt
        """
    }

    workflow {
        input_ch = channel.of('sample1', 'sample2', 'sample3')  // Moved inside workflow block
        PROCESS_FILES(input_ch)
    }
    ```

=== "Vorher"

    ```groovy title="badpractice_syntax.nf" hl_lines="3" linenums="1"
    #!/usr/bin/env nextflow

    input_ch = channel.of('sample1', 'sample2', 'sample3')  // ERROR: Channel defined outside workflow

    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_processed.txt"

        script:
        // Define variables in Groovy code before the script
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

#### Führe die Pipeline aus

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

Halte deine Eingabekanäle innerhalb des Workflow-Blocks definiert und folge im Allgemeinen allen anderen Empfehlungen, die die Extension macht.

### Fazit

Du kannst Syntaxfehler systematisch identifizieren und beheben, indem du Nextflow-Fehlermeldungen und visuelle IDE-Indikatoren verwendest. Häufige Syntaxfehler umfassen fehlende Klammern, falsche Prozess-Schlüsselwörter, undefinierte Variablen und unsachgemäße Verwendung von Bash- vs. Nextflow-Variablen. Die VSCode-Extension hilft, viele davon vor der Laufzeit zu erkennen. Mit diesen Syntax-Debugging-Fähigkeiten in deinem Werkzeugkasten wirst du in der Lage sein, die häufigsten Nextflow-Syntaxfehler schnell zu beheben und dich komplexeren Laufzeitproblemen zuzuwenden.

### Wie geht es weiter?

Lerne, komplexere Kanal-Strukturfehler zu debuggen, die auch bei korrekter Syntax auftreten.

---

## 2. Kanal-Strukturfehler

Kanal-Strukturfehler sind subtiler als Syntaxfehler, weil der Code syntaktisch korrekt ist, aber die Datenformen nicht mit dem übereinstimmen, was Prozesse erwarten. Nextflow wird versuchen, die Pipeline auszuführen, könnte aber feststellen, dass die Anzahl der Eingaben nicht mit den Erwartungen übereinstimmt und fehlschlagen. Diese Fehler treten typischerweise nur zur Laufzeit auf und erfordern ein Verständnis der Daten, die durch deinen Workflow fließen.

!!! tip "Debugging von Kanälen mit `.view()`"

    Denke während dieses Abschnitts daran, dass du den `.view()`-Operator verwenden kannst, um Kanalinhalte an jedem Punkt in deinem Workflow zu inspizieren. Dies ist eines der mächtigsten Debugging-Tools zum Verstehen von Kanal-Strukturproblemen. Wir werden diese Technik in Abschnitt 2.4 im Detail erkunden, aber fühle dich frei, sie zu verwenden, während du durch die Beispiele arbeitest.

    ```groovy
    my_channel.view()  // Shows what's flowing through the channel
    ```

### 2.1. Falsche Anzahl von Eingabekanälen

Dieser Fehler tritt auf, wenn du eine andere Anzahl von Kanälen übergibst, als ein Prozess erwartet.

#### Führe die Pipeline aus

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

#### Überprüfe den Code

Die Fehlermeldung besagt klar, dass der Aufruf 1 Argument erwartet, aber 2 erhalten hat, und zeigt auf Zeile 23. Schauen wir uns `bad_number_inputs.nf` an:

```groovy title="bad_number_inputs.nf" hl_lines="5 23" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    input:
        val sample_name  // Process expects only 1 input

    output:
        path "${sample_name}_output.txt"

    script:
    """
    echo "Processing ${sample_name}" > ${sample_name}_output.txt
    """
}

workflow {

    // Create two separate channels
    samples_ch = channel.of('sample1', 'sample2', 'sample3')
    files_ch = channel.of('file1.txt', 'file2.txt', 'file3.txt')

    // ERROR: Passing 2 channels but process expects only 1
    PROCESS_FILES(samples_ch, files_ch)
}
```

Du solltest den nicht übereinstimmenden `PROCESS_FILES`-Aufruf sehen, der mehrere Eingabekanäle liefert, wenn der Prozess nur einen definiert. Die VSCode-Extension wird auch den Prozessaufruf rot unterstreichen und eine Diagnosemeldung liefern, wenn du mit der Maus darüber fährst:

![Incorrect number of args message](img/incorrect_num_args.png)

#### Behebe den Code

Für dieses spezifische Beispiel erwartet der Prozess einen einzelnen Kanal und benötigt den zweiten Kanal nicht, also können wir es beheben, indem wir nur den `samples_ch`-Kanal übergeben:

=== "Danach"

    ```groovy title="bad_number_inputs.nf" hl_lines="23" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
            val sample_name  // Process expects only 1 input

        output:
            path "${sample_name}_output.txt"

        script:
        """
        echo "Processing ${sample_name}" > ${sample_name}_output.txt
        """
    }

    workflow {

        // Create two separate channels
        samples_ch = channel.of('sample1', 'sample2', 'sample3')
        files_ch = channel.of('file1.txt', 'file2.txt', 'file3.txt')

        // Fixed: Pass only the channel the process expects
        PROCESS_FILES(samples_ch)
    }
    ```

=== "Vorher"

    ```groovy title="bad_number_inputs.nf" hl_lines="5 23" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
            val sample_name  // Process expects only 1 input

        output:
            path "${sample_name}_output.txt"

        script:
        """
        echo "Processing ${sample_name}" > ${sample_name}_output.txt
        """
    }

    workflow {

        // Create two separate channels
        samples_ch = channel.of('sample1', 'sample2', 'sample3')
        files_ch = channel.of('file1.txt', 'file2.txt', 'file3.txt')

        // ERROR: Passing 2 channels but process expects only 1
        PROCESS_FILES(samples_ch, files_ch)
    }
    ```

#### Führe die Pipeline aus

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

Häufiger als in diesem Beispiel könntest du zusätzliche Eingaben zu einem Prozess hinzufügen und vergessen, den Workflow-Aufruf entsprechend zu aktualisieren, was zu dieser Art von Fehler führen kann. Glücklicherweise ist dies einer der leichter zu verstehenden und zu behebenden Fehler, da die Fehlermeldung ziemlich klar über die Nichtübereinstimmung ist.

### 2.2. Kanal-Erschöpfung (Prozess läuft seltener als erwartet)

Einige Kanal-Strukturfehler sind viel subtiler und erzeugen überhaupt keine Fehler. Wahrscheinlich der häufigste davon spiegelt eine Herausforderung wider, der neue Nextflow-Nutzer\*innen gegenüberstehen, wenn sie verstehen, dass Queue-Kanäle erschöpft werden und keine Elemente mehr haben können, was bedeutet, dass der Workflow vorzeitig endet.

#### Führe die Pipeline aus

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

Dieser Workflow wird ohne Fehler abgeschlossen, aber er verarbeitet nur eine einzelne Probe!

#### Überprüfe den Code

Schauen wir uns `exhausted.nf` an, um zu sehen, ob das richtig ist:

```groovy title="exhausted.nf" hl_lines="23 24" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    input:
    val reference
    val sample_name

    output:
    path "${output_prefix}.txt"

    script:
    // Define variables in Groovy code before the script
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

Der Prozess läuft nur einmal statt dreimal, weil der `reference_ch`-Kanal ein Queue-Kanal ist, der nach der ersten Prozessausführung erschöpft ist. Wenn ein Kanal erschöpft ist, stoppt der gesamte Prozess, auch wenn andere Kanäle noch Elemente haben.

Dies ist ein häufiges Muster, bei dem du eine einzelne Referenzdatei hast, die über mehrere Proben hinweg wiederverwendet werden muss. Die Lösung besteht darin, den Referenzkanal in einen Value-Kanal umzuwandeln, der unbegrenzt wiederverwendet werden kann.

#### Behebe den Code

Es gibt ein paar Möglichkeiten, dies zu beheben, abhängig davon, wie viele Dateien betroffen sind.

**Option 1**: Du hast eine einzelne Referenzdatei, die du oft wiederverwendest. Du kannst einfach einen Value-Kanal-Typ erstellen, der immer wieder verwendet werden kann. Es gibt drei Möglichkeiten, dies zu tun:

**1a** Verwende `channel.value()`:

```groovy title="exhausted.nf (fixed - Option 1a)" hl_lines="2" linenums="21"
workflow {
    reference_ch = channel.value('baseline_reference')  // Value channel can be reused
    input_ch = channel.of('sample1', 'sample2', 'sample3')

    PROCESS_FILES(reference_ch, input_ch)
}
```

**1b** Verwende den `first()`-[Operator](https://www.nextflow.io/docs/latest/reference/operator.html#first):

```groovy title="exhausted.nf (fixed - Option 1b)" hl_lines="2" linenums="21"
workflow {
    reference_ch = channel.of('baseline_reference').first()  // Convert to value channel
    input_ch = channel.of('sample1', 'sample2', 'sample3')

    PROCESS_FILES(reference_ch, input_ch)
}
```

**1c.** Verwende den `collect()`-[Operator](https://www.nextflow.io/docs/latest/reference/operator.html#collect):

```groovy title="exhausted.nf (fixed - Option 1c)" hl_lines="2" linenums="21"
workflow {
    reference_ch = channel.of('baseline_reference').collect()  // Convert to value channel
    input_ch = channel.of('sample1', 'sample2', 'sample3')

    PROCESS_FILES(reference_ch, input_ch)
}
```

**Option 2**: In komplexeren Szenarien, vielleicht wo du mehrere Referenzdateien für alle Proben im Probenkanal hast, kannst du den `combine`-Operator verwenden, um einen neuen Kanal zu erstellen, der die beiden Kanäle zu Tupeln kombiniert:

```groovy title="exhausted.nf (fixed - Option 2)" hl_lines="4" linenums="21"
workflow {
    reference_ch = channel.of('baseline_reference','other_reference')
    input_ch = channel.of('sample1', 'sample2', 'sample3')
    combined_ch = reference_ch.combine(input_ch)  // Creates cartesian product

    PROCESS_FILES(combined_ch)
}
```

Der `.combine()`-Operator erzeugt ein kartesisches Produkt der beiden Kanäle, sodass jedes Element in `reference_ch` mit jedem Element in `input_ch` gepaart wird. Dies ermöglicht es dem Prozess, für jede Probe zu laufen, während die Referenz weiterhin verwendet wird.

Dies erfordert eine Anpassung der Prozesseingabe. In unserem Beispiel müsste der Anfang der Prozessdefinition wie folgt angepasst werden:

```groovy title="exhausted.nf (fixed - Option 2)" hl_lines="5" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    input:
        tuple val(reference), val(sample_name)
```

Dieser Ansatz ist möglicherweise nicht in allen Situationen geeignet.

#### Führe die Pipeline aus

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

Du solltest jetzt sehen, dass alle drei Proben verarbeitet werden, statt nur einer.

### 2.3. Falsche Kanal-Inhaltsstruktur

Wenn Workflows ein gewisses Maß an Komplexität erreichen, kann es etwas schwierig sein, die internen Strukturen jedes Kanals im Auge zu behalten, und Menschen erzeugen häufig Nichtübereinstimmungen zwischen dem, was der Prozess erwartet, und dem, was der Kanal tatsächlich enthält. Dies ist subtiler als das Problem, das wir früher besprochen haben, wo die Anzahl der Kanäle falsch war. In diesem Fall kannst du die richtige Anzahl von Eingabekanälen haben, aber die interne Struktur eines oder mehrerer dieser Kanäle stimmt nicht mit dem überein, was der Prozess erwartet.

#### Führe die Pipeline aus

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

#### Überprüfe den Code

Die eckigen Klammern in der Fehlermeldung liefern hier den Hinweis - der Prozess behandelt das Tupel als einzelnen Wert, was nicht das ist, was wir wollen. Schauen wir uns `bad_channel_shape.nf` an:

```groovy title="bad_channel_shape.nf" hl_lines="5 20-22" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    input:
        val sample_name  // Expects single value, gets tuple

    output:
        path "${sample_name}_output.txt"

    script:
    """
    echo "Processing ${sample_name}" > ${sample_name}_output.txt
    """
}

workflow {

    // Channel emits tuples, but process expects single values
    input_ch = channel.of(
      ['sample1', 'file1.txt'],
      ['sample2', 'file2.txt'],
      ['sample3', 'file3.txt']
    )
    PROCESS_FILES(input_ch)
}
```

Du kannst sehen, dass wir einen Kanal erzeugen, der aus Tupeln besteht: `['sample1', 'file1.txt']`, aber der Prozess erwartet einen einzelnen Wert, `val sample_name`. Der ausgeführte Befehl zeigt, dass der Prozess versucht, eine Datei namens `[sample3, file3.txt]_output.txt` zu erstellen, was nicht die beabsichtigte Ausgabe ist.

#### Behebe den Code

Um dies zu beheben, könnten wir, wenn der Prozess beide Eingaben benötigt, den Prozess anpassen, um ein Tupel zu akzeptieren:

=== "Option 1: Tupel im Prozess akzeptieren"

    === "Danach"

        ```groovy title="bad_channel_shape.nf" hl_lines="5"  linenums="1"
        #!/usr/bin/env nextflow

        process PROCESS_FILES {
            input:
                tuple val(sample_name), val(file_name)  // Fixed: Accept tuple

            output:
                path "${sample_name}_output.txt"

            script:
            """
            echo "Processing ${sample_name}" > ${sample_name}_output.txt
            """
        }

        workflow {

            // Channel emits tuples, but process expects single values
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
                val sample_name  // Expects single value, gets tuple

            output:
                path "${sample_name}_output.txt"

            script:
            """
            echo "Processing ${sample_name}" > ${sample_name}_output.txt
            """
        }

        workflow {

            // Channel emits tuples, but process expects single values
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

            // Channel emits tuples, but process expects single values
            input_ch = channel.of(
              ['sample1', 'file1.txt'],
              ['sample2', 'file2.txt'],
              ['sample3', 'file3.txt']
            )
            PROCESS_FILES(input_ch.map { it[0] })  // Fixed: Extract first element
        }
        ```

    === "Vorher"

        ```groovy title="bad_channel_shape.nf" hl_lines="9" linenums="16"
        workflow {

            // Channel emits tuples, but process expects single values
            input_ch = channel.of(
              ['sample1', 'file1.txt'],
              ['sample2', 'file2.txt'],
              ['sample3', 'file3.txt']
            )
            PROCESS_FILES(input_ch)
        }
        ```

#### Führe die Pipeline aus

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

#### Verwendung von `.view()` zur Kanal-Inspektion

Das mächtigste Debugging-Tool für Kanäle ist der `.view()`-Operator. Mit `.view()` kannst du die Form deiner Kanäle in allen Phasen verstehen, um beim Debugging zu helfen.

#### Führe die Pipeline aus

Führe `bad_channel_shape_viewed.nf` aus, um dies in Aktion zu sehen:

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

#### Überprüfe den Code

Schauen wir uns `bad_channel_shape_viewed.nf` an, um zu sehen, wie `.view()` verwendet wird:

```groovy title="bad_channel_shape_viewed.nf" linenums="16" hl_lines="9 11"
workflow {

    // Channel emits tuples, but process expects single values
    input_ch = channel.of(
      ['sample1', 'file1.txt'],
      ['sample2', 'file2.txt'],
      ['sample3', 'file3.txt']
    )
    .view { "Channel content: $it" }  // Debug: Show original channel content
    .map { tuple -> tuple[0] }        // Transform: Extract first element
    .view { "After mapping: $it" }    // Debug: Show transformed channel content

    PROCESS_FILES(input_ch)
}
```

#### Behebe den Code

Um dich davor zu bewahren, in Zukunft übermäßig `.view()`-Operationen zu verwenden, um Kanalinhalte zu verstehen, ist es ratsam, einige Kommentare hinzuzufügen:

```groovy title="bad_channel_shape_viewed.nf (with comments)" linenums="16" hl_lines="8 9"
workflow {

    // Channel emits tuples, but process expects single values
    input_ch = channel.of(
            ['sample1', 'file1.txt'],
            ['sample2', 'file2.txt'],
            ['sample3', 'file3.txt'],
        ) // [sample_name, file_name]
        .map { tuple -> tuple[0] } // sample_name

    PROCESS_FILES(input_ch)
}
```

Dies wird wichtiger, wenn deine Workflows an Komplexität zunehmen und die Kanalstruktur undurchsichtiger wird.

#### Führe die Pipeline aus

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

Viele Kanal-Strukturfehler können mit gültiger Nextflow-Syntax erstellt werden. Du kannst Kanal-Strukturfehler debuggen, indem du den Datenfluss verstehst, `.view()`-Operatoren zur Inspektion verwendest und Fehlermeldungsmuster wie eckige Klammern erkennst, die auf unerwartete Tupelstrukturen hinweisen.

### Wie geht es weiter?

Lerne über Fehler, die durch Prozessdefinitionen verursacht werden.

---

## 3. Prozess-Strukturfehler

Die meisten Fehler, denen du im Zusammenhang mit Prozessen begegnest, beziehen sich auf Fehler, die du beim Formulieren des Befehls gemacht hast, oder auf Probleme mit der zugrunde liegenden Software. Allerdings kannst du, ähnlich wie bei den Kanalproblemen oben, Fehler in der Prozessdefinition machen, die nicht als Syntaxfehler gelten, aber zur Laufzeit Fehler verursachen.

### 3.1. Fehlende Ausgabedateien

Ein häufiger Fehler beim Schreiben von Prozessen ist, etwas zu tun, das eine Nichtübereinstimmung zwischen dem erzeugt, was der Prozess erwartet, und dem, was generiert wird.

#### Führe die Pipeline aus

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

#### Überprüfe den Code

Die Fehlermeldung zeigt an, dass der Prozess erwartet hat, eine Ausgabedatei namens `sample3.txt` zu erzeugen, aber das Script tatsächlich `sample3_output.txt` erstellt. Schauen wir uns die Prozessdefinition in `missing_output.nf` an:

```groovy title="missing_output.nf" linenums="3" hl_lines="6 10"
process PROCESS_FILES {
    input:
    val sample_name

    output:
    path "${sample_name}.txt"  // Expects: sample3.txt

    script:
    """
    echo "Processing ${sample_name}" > ${sample_name}_output.txt  // Creates: sample3_output.txt
    """
}
```

Du solltest sehen, dass es eine Nichtübereinstimmung zwischen dem Ausgabedateinamen im `output:`-Block und dem im Script verwendeten gibt. Diese Nichtübereinstimmung führt dazu, dass der Prozess fehlschlägt. Wenn du auf diese Art von Fehler stößt, gehe zurück und überprüfe, ob die Ausgaben zwischen deiner Prozessdefinition und deinem Ausgabeblock übereinstimmen.

Wenn das Problem immer noch nicht klar ist, überprüfe das Work-Verzeichnis selbst, um die tatsächlich erstellten Ausgabedateien zu identifizieren:

```bash
❯ ls -h work/02/9604d49fb8200a74d737c72a6c98ed
sample3_output.txt
```

Für dieses Beispiel würde dies uns verdeutlichen, dass ein `_output`-Suffix in den Ausgabedateinamen eingebaut wird, entgegen unserer `output:`-Definition.

#### Behebe den Code

Behebe die Nichtübereinstimmung, indem du den Ausgabedateinamen konsistent machst:

=== "Danach"

    ```groovy title="missing_output.nf" hl_lines="6 10" linenums="3"
    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_output.txt"  // Fixed: Match the script output

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
        path "${sample_name}.txt"  // Expects: sample3.txt

        script:
        """
        echo "Processing ${sample_name}" > ${sample_name}_output.txt  // Creates: sample3_output.txt
        """
    }
    ```

#### Führe die Pipeline aus

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

Eine weitere Fehlerklasse tritt aufgrund von Fehlern bei der Software-Bereitstellung auf. `missing_software.nf` ist ein syntaktisch gültiger Workflow, aber er hängt von externer Software ab, um den `cowpy`-Befehl bereitzustellen, den er verwendet.

#### Führe die Pipeline aus

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

Der Prozess hat keinen Zugriff auf den Befehl, den wir angeben. Manchmal liegt dies daran, dass ein Script im Workflow-`bin`-Verzeichnis vorhanden ist, aber nicht ausführbar gemacht wurde. Andere Male liegt es daran, dass die Software nicht im Container oder in der Umgebung installiert ist, in der der Workflow läuft.

#### Überprüfe den Code

Achte auf den `127`-Exit-Code - er sagt dir genau das Problem. Schauen wir uns `missing_software.nf` an:

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

#### Behebe den Code

Wir waren hier ein wenig unaufrichtig, und es ist eigentlich nichts mit dem Code falsch. Wir müssen nur die notwendige Konfiguration angeben, um den Prozess so auszuführen, dass er Zugriff auf den fraglichen Befehl hat. In diesem Fall hat der Prozess eine Container-Definition, also müssen wir nur den Workflow mit aktiviertem Docker ausführen.

#### Führe die Pipeline aus

Wir haben ein Docker-Profil für dich in `nextflow.config` eingerichtet, sodass du den Workflow ausführen kannst mit:

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

!!! note

    Um mehr darüber zu erfahren, wie Nextflow Container verwendet, siehe [Hello Nextflow](../hello_nextflow/05_hello_containers.md)

### 3.3. Schlechte Ressourcenkonfiguration

In der Produktionsnutzung wirst du Ressourcen für deine Prozesse konfigurieren. Zum Beispiel definiert `memory` die maximale Menge an Speicher, die deinem Prozess zur Verfügung steht, und wenn der Prozess diese überschreitet, wird dein Scheduler typischerweise den Prozess beenden und einen Exit-Code von `137` zurückgeben. Wir können das hier nicht demonstrieren, weil wir den `local`-Executor verwenden, aber wir können etwas Ähnliches mit `time` zeigen.

#### Führe die Pipeline aus

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

#### Überprüfe den Code

Schauen wir uns `bad_resources.nf` an:

```groovy title="bad_resources.nf" linenums="3" hl_lines="3"
process PROCESS_FILES {

    time '1 ms'  // ERROR: Unrealistic time limit

    input:
    val sample_name

    output:
    path "${sample_name}_output.txt"

    script:
    """
    sleep 1  // Takes 1 second, but time limit is 1ms
    cowpy ${sample_name} > ${sample_name}_output.txt
    """
}
```

Wir wissen, dass der Prozess länger als eine Sekunde dauern wird (wir haben einen Sleep eingefügt, um sicherzustellen), aber der Prozess ist so eingestellt, dass er nach 1 Millisekunde ein Timeout hat. Jemand war ein wenig unrealistisch mit seiner Konfiguration!

#### Behebe den Code

Erhöhe das Zeitlimit auf einen realistischen Wert:

=== "Danach"

    ```groovy title="bad_resources.nf" hl_lines="3" linenums="3"
    process PROCESS_FILES {

        time '100 s'  // Fixed: Realistic time limit

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

        time '1 ms'  // ERROR: Unrealistic time limit

        input:
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        """
        sleep 1  // Takes 1 second, but time limit is 1ms
        cowpy ${sample_name} > ${sample_name}_output.txt
        """
    }
    ```

#### Führe die Pipeline aus

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

Wenn du sicherstellst, dass du deine Fehlermeldungen liest, sollten Fehler wie dieser dich nicht zu lange verwirren. Aber stelle sicher, dass du die Ressourcenanforderungen der Befehle verstehst, die du ausführst, damit du deine Ressourcen-Direktiven angemessen konfigurieren kannst.

### 3.4. Prozess-Debugging-Techniken

Wenn Prozesse fehlschlagen oder sich unerwartet verhalten, benötigst du systematische Techniken, um zu untersuchen, was schiefgelaufen ist. Das Work-Verzeichnis enthält alle Informationen, die du zum Debuggen der Prozessausführung benötigst.

#### Verwendung der Work-Verzeichnis-Inspektion

Das mächtigste Debugging-Tool für Prozesse ist die Untersuchung des Work-Verzeichnisses. Wenn ein Prozess fehlschlägt, erstellt Nextflow ein Work-Verzeichnis für diese spezifische Prozessausführung, das alle Dateien enthält, die zum Verstehen dessen, was passiert ist, benötigt werden.

#### Führe die Pipeline aus

Verwenden wir das `missing_output.nf`-Beispiel von früher, um die Work-Verzeichnis-Inspektion zu demonstrieren (erzeuge bei Bedarf erneut eine Ausgabenamen-Nichtübereinstimmung):

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

#### Überprüfe das Work-Verzeichnis

Wenn du diesen Fehler erhältst, enthält das Work-Verzeichnis alle Debugging-Informationen. Finde den Work-Verzeichnis-Pfad aus der Fehlermeldung und untersuche seinen Inhalt:

```bash
# Finde das Work-Verzeichnis aus der Fehlermeldung
ls work/02/9604d49fb8200a74d737c72a6c98ed/
```

Du kannst dann die Schlüsseldateien untersuchen:

##### Überprüfe das Befehls-Script

Die `.command.sh`-Datei zeigt genau, welcher Befehl ausgeführt wurde:

```bash
# Zeige den ausgeführten Befehl an
cat work/02/9604d49fb8200a74d737c72a6c98ed/.command.sh
```

Dies zeigt:

- **Variablensubstitution**: Ob Nextflow-Variablen richtig erweitert wurden
- **Dateipfade**: Ob Eingabedateien korrekt lokalisiert wurden
- **Befehlsstruktur**: Ob die Script-Syntax korrekt ist

Häufige Probleme, auf die man achten sollte:

- **Fehlende Anführungszeichen**: Variablen mit Leerzeichen benötigen richtige Anführungszeichen
- **Falsche Dateipfade**: Eingabedateien, die nicht existieren oder an falschen Orten sind
- **Falsche Variablennamen**: Tippfehler in Variablenreferenzen
- **Fehlendes Umgebungs-Setup**: Befehle, die von spezifischen Umgebungen abhängen

##### Überprüfe die Fehlerausgabe

Die `.command.err`-Datei enthält die tatsächlichen Fehlermeldungen:

```bash
# Zeige Fehlerausgabe an
cat work/02/9604d49fb8200a74d737c72a6c98ed/.command.err
```

Diese Datei zeigt:

- **Exit-Codes**: 127 (Befehl nicht gefunden), 137 (beendet), usw.
- **Berechtigungsfehler**: Dateizugriffsprobleme
- **Software-Fehler**: Anwendungsspezifische Fehlermeldungen
- **Ressourcenfehler**: Speicher-/Zeitlimit überschritten

##### Überprüfe die Standardausgabe

Die `.command.out`-Datei zeigt, was dein Befehl produziert hat:

```bash
# Zeige Standardausgabe an
cat work/02/9604d49fb8200a74d737c72a6c98ed/.command.out
```

Dies hilft zu überprüfen:

- **Erwartete Ausgabe**: Ob der Befehl die richtigen Ergebnisse produziert hat
- **Teilweise Ausführung**: Ob der Befehl gestartet, aber auf halbem Weg fehlgeschlagen ist
- **Debug-Informationen**: Jegliche Diagnoseausgabe aus deinem Script

##### Überprüfe den Exit-Code

Die `.exitcode`-Datei enthält den Exit-Code für den Prozess:

```bash
# Zeige Exit-Code an
cat work/*/*/.exitcode
```

Häufige Exit-Codes und ihre Bedeutungen:

- **Exit-Code 127**: Befehl nicht gefunden - überprüfe Software-Installation
- **Exit-Code 137**: Prozess beendet - überprüfe Speicher-/Zeitlimits

##### Überprüfe Dateiexistenz

Wenn Prozesse aufgrund fehlender Ausgabedateien fehlschlagen, überprüfe, welche Dateien tatsächlich erstellt wurden:

```bash
# Liste alle Dateien im Work-Verzeichnis auf
ls -la work/02/9604d49fb8200a74d737c72a6c98ed/
```

Dies hilft zu identifizieren:

- **Dateinamen-Nichtübereinstimmungen**: Ausgabedateien mit anderen Namen als erwartet
- **Berechtigungsprobleme**: Dateien, die nicht erstellt werden konnten
- **Pfadprobleme**: Dateien, die in falschen Verzeichnissen erstellt wurden

In unserem früheren Beispiel bestätigte uns dies, dass während unsere erwartete `sample3.txt` nicht vorhanden war, `sample3_output.txt` vorhanden war:

```bash
❯ ls -h work/02/9604d49fb8200a74d737c72a6c98ed
sample3_output.txt
```

### Fazit

Prozess-Debugging erfordert die Untersuchung von Work-Verzeichnissen, um zu verstehen, was schiefgelaufen ist. Schlüsseldateien umfassen `.command.sh` (das ausgeführte Script), `.command.err` (Fehlermeldungen) und `.command.out` (Standardausgabe). Exit-Codes wie 127 (Befehl nicht gefunden) und 137 (Prozess beendet) liefern sofortige diagnostische Hinweise auf die Art des Fehlers.

### Wie geht es weiter?

Lerne über Nextflows eingebaute Debugging-Tools und systematische Ansätze zur Fehlerbehebung.

---

## 4. Eingebaute Debugging-Tools und fortgeschrittene Techniken

Nextflow bietet mehrere leistungsstarke eingebaute Tools zum Debuggen und Analysieren der Workflow-Ausführung. Diese Tools helfen dir zu verstehen, was schiefgelaufen ist, wo es schiefgelaufen ist und wie man es effizient behebt.

### 4.1. Echtzeit-Prozessausgabe

Manchmal musst du sehen, was in laufenden Prozessen passiert. Du kannst Echtzeit-Prozessausgabe aktivieren, die dir genau zeigt, was jede Aufgabe während der Ausführung tut.

#### Führe die Pipeline aus

`bad_channel_shape_viewed.nf` aus unseren früheren Beispielen hat Kanalinhalte mit `.view()` ausgegeben, aber wir können auch die `debug`-Direktive verwenden, um Variablen aus dem Prozess selbst auszugeben, was wir in `bad_channel_shape_viewed_debug.nf` demonstrieren. Führe den Workflow aus:

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

#### Überprüfe den Code

Schauen wir uns `bad_channel_shape_viewed_debug.nf` an, um zu sehen, wie die `debug`-Direktive funktioniert:

```groovy title="bad_channel_shape_viewed_debug.nf" linenums="3" hl_lines="2"
process PROCESS_FILES {
    debug true  // Enable real-time output

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

Die `debug`-Direktive kann eine schnelle und bequeme Möglichkeit sein, die Umgebung eines Prozesses zu verstehen.

### 4.2. Preview-Modus

Manchmal möchtest du Probleme erkennen, bevor Prozesse laufen. Nextflow bietet ein Flag für diese Art von proaktivem Debugging: `-preview`.

#### Führe die Pipeline aus

Der Preview-Modus ermöglicht es dir, Workflow-Logik zu testen, ohne Befehle auszuführen. Dies kann sehr nützlich sein, um schnell die Struktur deines Workflows zu überprüfen und sicherzustellen, dass Prozesse korrekt verbunden sind, ohne tatsächliche Befehle auszuführen.

!!! note

    Wenn du `bad_syntax.nf` früher behoben hast, führe den Syntaxfehler wieder ein, indem du die schließende Klammer nach dem Script-Block entfernst, bevor du diesen Befehl ausführst.

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

Der Preview-Modus ist besonders nützlich, um Syntaxfehler früh zu erkennen, ohne Prozesse auszuführen. Er validiert die Workflow-Struktur und Prozessverbindungen vor der Ausführung.

### 4.3. Stub-Running zum Testen der Logik

Manchmal sind Fehler schwer zu debuggen, weil Befehle zu lange dauern, spezielle Software benötigen oder aus komplexen Gründen fehlschlagen. Stub-Running ermöglicht es dir, Workflow-Logik zu testen, ohne die tatsächlichen Befehle auszuführen.

#### Führe die Pipeline aus

Wenn du einen Nextflow-Prozess entwickelst, kannst du die `stub`-Direktive verwenden, um 'Dummy'-Befehle zu definieren, die Ausgaben in der richtigen Form generieren, ohne den echten Befehl auszuführen. Dieser Ansatz ist besonders wertvoll, wenn du überprüfen möchtest, dass deine Workflow-Logik korrekt ist, bevor du dich mit den Komplexitäten der tatsächlichen Software befasst.

Erinnerst du dich zum Beispiel an unser `missing_software.nf` von früher? Das, bei dem wir fehlende Software hatten, die verhinderte, dass der Workflow lief, bis wir `-profile docker` hinzufügten? `missing_software_with_stub.nf` ist ein sehr ähnlicher Workflow. Wenn wir ihn auf die gleiche Weise ausführen, werden wir denselben Fehler erzeugen:

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

Dieser Workflow wird jedoch keine Fehler erzeugen, wenn wir ihn mit `-stub-run` ausführen, auch ohne das `docker`-Profil:

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

#### Überprüfe den Code

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

Im Vergleich zu `missing_software.nf` hat dieser Prozess eine `stub:`-Direktive, die einen Befehl spezifiziert, der anstelle des in `script:` angegebenen verwendet werden soll, falls Nextflow im Stub-Modus ausgeführt wird.

Der `touch`-Befehl, den wir hier verwenden, hängt nicht von Software oder geeigneten Eingaben ab und wird in allen Situationen laufen, sodass wir Workflow-Logik debuggen können, ohne uns um die Prozess-Interna zu sorgen.

**Stub-Running hilft beim Debuggen:**

- Kanalstruktur und Datenfluss
- Prozessverbindungen und Abhängigkeiten
- Parameter-Propagierung
- Workflow-Logik ohne Software-Abhängigkeiten

### 4.4. Systematischer Debugging-Ansatz

Jetzt, da du einzelne Debugging-Techniken gelernt hast - von Trace-Dateien und Work-Verzeichnissen bis zu Preview-Modus, Stub-Running und Ressourcen-Monitoring - lass uns sie zu einer systematischen Methodik zusammenfügen. Ein strukturierter Ansatz verhindert, dass du von komplexen Fehlern überwältigt wirst, und stellt sicher, dass du keine wichtigen Hinweise übersiehst.

Diese Methodik kombiniert alle Tools, die wir behandelt haben, zu einem effizienten Workflow:

**Vier-Phasen-Debugging-Methode:**

**Phase 1: Syntaxfehler-Behebung (5 Minuten)**

1. Überprüfe auf rote Unterstreichungen in VSCode oder deiner IDE
2. Führe `nextflow run workflow.nf -preview` aus, um Syntaxprobleme zu identifizieren
3. Behebe alle Syntaxfehler (fehlende Klammern, nachgestellte Kommas usw.)
4. Stelle sicher, dass der Workflow erfolgreich geparst wird, bevor du fortfährst

**Phase 2: Schnelle Bewertung (5 Minuten)**

1. Lies Laufzeit-Fehlermeldungen sorgfältig
2. Überprüfe, ob es sich um einen Laufzeit-, Logik- oder Ressourcenfehler handelt
3. Verwende den Preview-Modus, um grundlegende Workflow-Logik zu testen

**Phase 3: Detaillierte Untersuchung (15-30 Minuten)**

1. Finde das Work-Verzeichnis der fehlgeschlagenen Aufgabe
2. Untersuche Log-Dateien
3. Füge `.view()`-Operatoren hinzu, um Kanäle zu inspizieren
4. Verwende `-stub-run`, um Workflow-Logik ohne Ausführung zu testen

**Phase 4: Beheben und Validieren (15 Minuten)**

1. Mache minimale gezielte Korrekturen
2. Teste mit Resume: `nextflow run workflow.nf -resume`
3. Überprüfe vollständige Workflow-Ausführung

!!! tip "Verwendung von Resume für effizientes Debugging"

    Sobald du ein Problem identifiziert hast, benötigst du eine effiziente Möglichkeit, deine Korrekturen zu testen, ohne Zeit mit dem erneuten Ausführen erfolgreicher Teile deines Workflows zu verschwenden. Nextflows `-resume`-Funktionalität ist beim Debugging von unschätzbarem Wert.

    Du wirst `-resume` begegnet sein, wenn du [Hello Nextflow](../hello_nextflow/) durchgearbeitet hast, und es ist wichtig, dass du es beim Debugging gut nutzt, um dir das Warten zu ersparen, während die Prozesse vor deinem Problemprozess laufen.

    **Resume-Debugging-Strategie:**

    1. Führe Workflow bis zum Fehler aus
    2. Untersuche Work-Verzeichnis für fehlgeschlagene Aufgabe
    3. Behebe das spezifische Problem
    4. Resume, um nur die Korrektur zu testen
    5. Wiederhole, bis der Workflow abgeschlossen ist

#### Debugging-Konfigurationsprofil

Um diesen systematischen Ansatz noch effizienter zu machen, kannst du eine dedizierte Debugging-Konfiguration erstellen, die automatisch alle Tools aktiviert, die du benötigst:

```groovy title="nextflow.config (debug profile)" linenums="1"
profiles {
    debug {
        process {
            debug = true
            cleanup = false

            // Conservative resources for debugging
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

Dieses Profil aktiviert Echtzeit-Ausgabe, bewahrt Work-Verzeichnisse und begrenzt Parallelisierung für einfacheres Debugging.

### 4.5. Praktische Debugging-Übung

Jetzt ist es an der Zeit, den systematischen Debugging-Ansatz in die Praxis umzusetzen. Der Workflow `buggy_workflow.nf` enthält mehrere häufige Fehler, die die Arten von Problemen repräsentieren, denen du in der realen Entwicklung begegnen wirst.

!!! exercise

    Verwende den systematischen Debugging-Ansatz, um alle Fehler in `buggy_workflow.nf` zu identifizieren und zu beheben. Dieser Workflow versucht, Probendaten aus einer CSV-Datei zu verarbeiten, enthält aber mehrere absichtliche Fehler, die häufige Debugging-Szenarien repräsentieren.

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

        Dieser kryptische Fehler zeigt ein Parsing-Problem um Zeile 11-12 im `params{}`-Block an. Der v2-Parser erkennt strukturelle Probleme früh.

    Wende die Vier-Phasen-Debugging-Methode an, die du gelernt hast:

    **Phase 1: Syntaxfehler-Behebung**
    - Überprüfe auf rote Unterstreichungen in VSCode oder deiner IDE
    - Führe `nextflow run workflow.nf -preview` aus, um Syntaxprobleme zu identifizieren
    - Behebe alle Syntaxfehler (fehlende Klammern, nachgestellte Kommas usw.)
    - Stelle sicher, dass der Workflow erfolgreich geparst wird, bevor du fortfährst

    **Phase 2: Schnelle Bewertung**
    - Lies Laufzeit-Fehlermeldungen sorgfältig
    - Identifiziere, ob Fehler laufzeit-, logik- oder ressourcenbezogen sind
    - Verwende `-preview`-Modus, um grundlegende Workflow-Logik zu testen

    **Phase 3: Detaillierte Untersuchung**
    - Untersuche Work-Verzeichnisse für fehlgeschlagene Aufgaben
    - Füge `.view()`-Operatoren hinzu, um Kanäle zu inspizieren
    - Überprüfe Log-Dateien in Work-Verzeichnissen
    - Verwende `-stub-run`, um Workflow-Logik ohne Ausführung zu testen

    **Phase 4: Beheben und Validieren**
    - Mache gezielte Korrekturen
    - Verwende `-resume`, um Korrekturen effizient zu testen
    - Überprüfe vollständige Workflow-Ausführung

    **Debugging-Tools zur Verfügung:**
    ```bash
    # Preview-Modus für Syntaxprüfung
    nextflow run buggy_workflow.nf -preview

    # Debug-Profil für detaillierte Ausgabe
    nextflow run buggy_workflow.nf -profile debug

    # Stub-Running für Logik-Tests
    nextflow run buggy_workflow.nf -stub-run

    # Resume nach Korrekturen
    nextflow run buggy_workflow.nf -resume
    ```

    ??? solution
        `buggy_workflow.nf` enthält 9 oder 10 verschiedene Fehler (je nachdem, wie man zählt), die alle wichtigen Debugging-Kategorien abdecken. Hier ist eine systematische Aufschlüsselung jedes Fehlers und wie man ihn behebt

        Beginnen wir mit den Syntaxfehlern:

        **Fehler 1: Syntaxfehler - Nachgestelltes Komma**
        ```groovy linenums="21"
        output:
            path "${sample_id}_result.txt",  // ERROR: Trailing comma
        ```
        **Korrektur:** Entferne das nachgestellte Komma
        ```groovy linenums="21"
        output:
            path "${sample_id}_result.txt"
        ```

        **Fehler 2: Syntaxfehler - Fehlende schließende Klammer**
        ```groovy linenums="24"
        script:
        """
        echo "Processing: ${sample}"
        cat ${input_file} > ${sample}_result.txt
        """
        // ERROR: Missing closing brace for processFiles process
        ```
        **Korrektur:** Füge die fehlende schließende Klammer hinzu
        ```groovy linenums="29"
        """
        echo "Processing: ${sample_id}"
        cat ${input_file} > ${sample_id}_result.txt
        """
        }  // Add missing closing brace
        ```

        **Fehler 3: Variablennamenfehler**
        ```groovy linenums="26"
        echo "Processing: ${sample}"     // ERROR: should be sample_id
        cat ${input_file} > ${sample}_result.txt  // ERROR: should be sample_id
        ```
        **Korrektur:** Verwende den korrekten Eingabevariablennamen
        ```groovy linenums="26"
        echo "Processing: ${sample_id}"
        cat ${input_file} > ${sample_id}_result.txt
        ```

        **Fehler 4: Undefinierter Variablenfehler**
        ```groovy linenums="87"
        heavy_ch = heavyProcess(sample_ids)  // ERROR: sample_ids undefined
        ```
        **Korrektur:** Verwende den korrekten Kanal und extrahiere Proben-IDs
        ```groovy linenums="87"
        heavy_ch = heavyProcess(input_ch)
        ```

        An diesem Punkt wird der Workflow laufen, aber wir werden immer noch Fehler erhalten (z.B. `Path value cannot be null` in `processFiles`), verursacht durch schlechte Kanalstruktur.

        **Fehler 5: Kanal-Strukturfehler - Falsche Map-Ausgabe**
        ```groovy linenums="83"
        .map { row -> row.sample_id }  // ERROR: processFiles expects tuple
        ```
        **Korrektur:** Gib die Tupelstruktur zurück, die processFiles erwartet
        ```groovy linenums="83"
        .map { row -> [row.sample_id, file(row.fastq_path)] }
        ```

        Aber dies wird unsere Korrektur für das Ausführen von `heavyProcess()` oben brechen, also müssen wir eine Map verwenden, um nur die Proben-IDs an diesen Prozess zu übergeben:

        **Fehler 6: Schlechte Kanalstruktur für heavyProcess**
        ```groovy linenums="87"
        heavy_ch = heavyProcess(input_ch)  // ERROR: input_ch now has 2 elements per emission- heavyProcess only needs 1 (the first)
        ```
        **Korrektur:** Verwende den korrekten Kanal und extrahiere Proben-IDs
        ```groovy linenums="87"
        heavy_ch = heavyProcess(input_ch.map{it[0]})
        ```

        Jetzt kommen wir etwas weiter, erhalten aber einen Fehler über `No such variable: i`, weil wir eine Bash-Variable nicht escaped haben.

        **Fehler 7: Bash-Variablen-Escaping-Fehler**
        ```groovy linenums="48"
        echo "Heavy computation $i for ${sample_id}"  // ERROR: $i not escaped
        ```
        **Korrektur:** Escape die Bash-Variable
        ```groovy linenums="48"
        echo "Heavy computation \${i} for ${sample_id}"
        ```

        Jetzt erhalten wir `Process exceeded running time limit (1ms)`, also beheben wir das Laufzeitlimit für den relevanten Prozess:

        **Fehler 8: Ressourcen-Konfigurationsfehler**
        ```groovy linenums="36"
        time '1 ms'  // ERROR: Unrealistic time limit
        ```
        **Korrektur:** Erhöhe auf ein realistisches Zeitlimit
        ```groovy linenums="36"
        time '100 s'
        ```

        Als Nächstes haben wir einen `Missing output file(s)`-Fehler zu beheben:

        **Fehler 9: Ausgabedateiname-Nichtübereinstimmung**
        ```groovy linenums="49"
        done > ${sample_id}.txt  // ERROR: Wrong filename, should match output declaration
        ```
        **Korrektur:** Passe die Ausgabe-Deklaration an
        ```groovy linenums="49"
        done > ${sample_id}_heavy.txt
        ```

        Die ersten beiden Prozesse liefen, aber nicht der dritte.

        **Fehler 10: Ausgabedateiname-Nichtübereinstimmung**
        ```groovy linenums="88"
        file_ch = channel.fromPath("*.txt") // Error: attempting to take input from the pwd rather than a process
        handleFiles(file_ch)
        ```
        **Korrektur:** Nimm die Ausgabe vom vorherigen Prozess
        ```groovy linenums="88"
        handleFiles(heavyProcess.out)
        ```

        Damit sollte der gesamte Workflow laufen.

        **Vollständiger korrigierter Workflow:**
        ```groovy linenums="1"
        #!/usr/bin/env nextflow

        /*
        * Fehlerhafter Workflow für Debugging-Übungen
        * Dieser Workflow enthält mehrere absichtliche Fehler zu Lernzwecken
        */

        params{
            // Parameters with missing validation
            input: Path = 'data/sample_data.csv'
            output: String = 'results'
        }

        /*
        * Prozess mit Eingabe-/Ausgabe-Nichtübereinstimmung
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
        * Prozess mit Dateibehandlungsproblemen
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

            // Channel with incorrect usage
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

- **Syntaxfehler**: Fehlende Klammern, nachgestellte Kommas, undefinierte Variablen
- **Kanal-Strukturfehler**: Falsche Datenformen, undefinierte Kanäle
- **Prozessfehler**: Ausgabedatei-Nichtübereinstimmungen, Variablen-Escaping
- **Ressourcenfehler**: Unrealistische Zeitlimits

**Wichtige Debugging-Lektionen:**

1. **Lies Fehlermeldungen sorgfältig** - sie zeigen oft direkt auf das Problem
2. **Verwende systematische Ansätze** - behebe einen Fehler nach dem anderen und teste mit `-resume`
3. **Verstehe den Datenfluss** - Kanal-Strukturfehler sind oft die subtilsten
4. **Überprüfe Work-Verzeichnisse** - wenn Prozesse fehlschlagen, sagen dir die Logs genau, was schiefgelaufen ist

---

## Zusammenfassung

In dieser Side Quest hast du eine Reihe systematischer Techniken zum Debuggen von Nextflow-Workflows gelernt.
Die Anwendung dieser Techniken in deiner eigenen Arbeit wird es dir ermöglichen, weniger Zeit im Kampf mit deinem Computer zu verbringen, Probleme schneller zu lösen und dich vor zukünftigen Problemen zu schützen.

### Wichtige Muster

**1. Wie man Syntaxfehler identifiziert und behebt**:

- Interpretation von Nextflow-Fehlermeldungen und Lokalisierung von Problemen
- Häufige Syntaxfehler: fehlende Klammern, falsche Schlüsselwörter, undefinierte Variablen
- Unterscheidung zwischen Nextflow- (Groovy-) und Bash-Variablen
- Verwendung von VS Code-Extension-Features zur frühen Fehlererkennung

```groovy
// Missing brace - look for red underlines in IDE
process FOO {
    script:
    """
    echo "hello"
    """
// } <-- missing!

// Wrong keyword
inputs:  // Should be 'input:'

// Undefined variable - escape with backslash for Bash variables
echo "${undefined_var}"      // Nextflow variable (error if not defined)
echo "\${bash_var}"          // Bash variable (escaped)
```

**2. Wie man Kanal-Strukturprobleme debuggt**:

- Verstehen von Kanal-Kardinalität und Erschöpfungsproblemen
- Debugging von Kanal-Inhaltsstruktur-Nichtübereinstimmungen
- Verwendung von `.view()`-Operatoren zur Kanal-Inspektion
- Erkennung von Fehlermustern wie eckigen Klammern in der Ausgabe

```groovy
// Inspect channel content
my_channel.view { "Content: $it" }

// Convert queue to value channel (prevents exhaustion)
reference_ch = channel.value('ref.fa')
// or
reference_ch = channel.of('ref.fa').first()
```

**3. Wie man Prozessausführungsprobleme behebt**:

- Diagnose von Fehlern bei fehlenden Ausgabedateien
- Verstehen von Exit-Codes (127 für fehlende Software, 137 für Speicherprobleme)
- Untersuchung von Work-Verzeichnissen und Befehlsdateien
- Angemessene Konfiguration von Ressourcen

```bash
# Überprüfe, was tatsächlich ausgeführt wurde
cat work/ab/cdef12/.command.sh

# Überprüfe Fehlerausgabe
cat work/ab/cdef12/.command.err

# Exit-Code 127 = Befehl nicht gefunden
# Exit-Code 137 = beendet (Speicher-/Zeitlimit)
```

**4. Wie man Nextflows eingebaute Debugging-Tools verwendet**:

- Nutzung von Preview-Modus und Echtzeit-Debugging
- Implementierung von Stub-Running für Logik-Tests
- Anwendung von Resume für effiziente Debugging-Zyklen
- Befolgung einer vierphasigen systematischen Debugging-Methodik

!!! tip "Schnelle Debugging-Referenz"

    **Syntaxfehler?** → Überprüfe VSCode-Warnungen, führe `nextflow run workflow.nf -preview` aus

    **Kanalprobleme?** → Verwende `.view()` zur Inhaltsinspektion: `my_channel.view()`

    **Prozessfehler?** → Überprüfe Work-Verzeichnis-Dateien:

    - `.command.sh` - das ausgeführte Script
    - `.command.err` - Fehlermeldungen
    - `.exitcode` - Exit-Status (127 = Befehl nicht gefunden, 137 = beendet)

    **Mysteriöses Verhalten?** → Führe mit `-stub-run` aus, um Workflow-Logik zu testen

    **Korrekturen gemacht?** → Verwende `-resume`, um Zeit beim Testen zu sparen: `nextflow run workflow.nf -resume`

---

### Zusätzliche Ressourcen

- [Nextflow Troubleshooting-Guide](https://www.nextflow.io/docs/latest/troubleshooting.html): Offizielle Troubleshooting-Dokumentation
- [Understanding Nextflow channels](https://www.nextflow.io/docs/latest/channel.html): Tiefgehende Betrachtung von Kanaltypen und -verhalten
- [Process directives reference](https://www.nextflow.io/docs/latest/process.html#directives): Alle verfügbaren Prozesskonfigurationsoptionen
- [nf-test](https://www.nf-test.com/): Testing-Framework für Nextflow-Pipelines
- [Nextflow Slack-Community](https://www.nextflow.io/slack-invite.html): Hilfe von der Community erhalten

Für Produktions-Workflows erwäge:

- Einrichtung von [Seqera Platform](https://seqera.io/platform/) für Monitoring und Debugging im großen Maßstab
- Verwendung von [Wave-Containern](https://seqera.io/wave/) für reproduzierbare Software-Umgebungen

**Denke daran:** Effektives Debugging ist eine Fähigkeit, die sich mit der Praxis verbessert. Die systematische Methodik und das umfassende Toolkit, das du hier erworben hast, werden dir während deiner gesamten Nextflow-Entwicklungsreise gute Dienste leisten.

---

## Wie geht es weiter?

Kehre zum [Menü der Side Quests](./index.md) zurück oder klicke auf den Button unten rechts auf der Seite, um zum nächsten Thema in der Liste zu gelangen.
