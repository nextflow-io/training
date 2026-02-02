# Debugging von Workflows

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Debugging ist eine kritische Fähigkeit, die dir Stunden der Frustration ersparen und dich zu einer effektiveren Nextflow-Entwickler\*in machen kann. Im Laufe deiner Karriere, besonders wenn du gerade anfängst, wirst du beim Erstellen und Warten deiner Workflows auf Bugs stoßen. Systematische Debugging-Ansätze zu lernen wird dir helfen, Probleme schnell zu identifizieren und zu lösen.

### Lernziele

In dieser Side Quest werden wir **systematische Debugging-Techniken** für Nextflow-Workflows erkunden:

- **Debugging von Syntaxfehlern**: Effektive Nutzung von IDE-Features und Nextflow-Fehlermeldungen
- **Channel-Debugging**: Diagnose von Datenfluss-Problemen und Channel-Strukturproblemen
- **Process-Debugging**: Untersuchung von Ausführungsfehlern und Ressourcenproblemen
- **Integrierte Debugging-Tools**: Nutzung von Nextflows Preview-Modus, Stub-Running und Work-Verzeichnissen
- **Systematische Ansätze**: Eine Vier-Phasen-Methodik für effizientes Debugging

Am Ende wirst du eine robuste Debugging-Methodik haben, die frustrierende Fehlermeldungen in klare Wegweiser für Lösungen verwandelt.

### Voraussetzungen

Bevor du diese Side Quest angehst, solltest du:

- Das [Hello Nextflow](../hello_nextflow/README.md) Tutorial oder einen gleichwertigen Einsteigerkurs abgeschlossen haben
- Mit grundlegenden Nextflow-Konzepten und -Mechanismen (Processes, Channels, Operatoren) vertraut sein

**Optional:** Wir empfehlen, zuerst die Side Quest [IDE Features for Nextflow Development](./ide_features.md) zu absolvieren.
Diese bietet eine umfassende Abdeckung von IDE-Features, die Debugging unterstützen (Syntaxhervorhebung, Fehlererkennung usw.), die wir hier intensiv nutzen werden.

---

## 0. Erste Schritte

#### Öffne den Training-Codespace

Falls du es noch nicht getan hast, öffne die Trainingsumgebung wie im [Environment Setup](../envsetup/index.md) beschrieben.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Wechsle in das Projektverzeichnis

Lass uns in das Verzeichnis wechseln, in dem sich die Dateien für dieses Tutorial befinden.

```bash
cd side-quests/debugging
```

Du kannst VSCode so einstellen, dass es sich auf dieses Verzeichnis fokussiert:

```bash
code .
```

#### Überprüfe die Materialien

Du findest eine Reihe von Beispiel-Workflows mit verschiedenen Arten von Bugs, die wir zum Üben verwenden werden:

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

Diese Dateien repräsentieren gängige Debugging-Szenarien, denen du in der realen Entwicklung begegnen wirst.

#### Überprüfe die Aufgabenstellung

Deine Herausforderung ist es, jeden Workflow auszuführen, die Fehler zu identifizieren und sie zu beheben.

Für jeden fehlerhaften Workflow:

1. **Führe den Workflow aus** und beobachte den Fehler
2. **Analysiere die Fehlermeldung**: Was sagt dir Nextflow?
3. **Lokalisiere das Problem** im Code anhand der bereitgestellten Hinweise
4. **Behebe den Bug** und verifiziere, dass deine Lösung funktioniert
5. **Setze die Datei zurück**, bevor du zum nächsten Abschnitt übergehst (verwende `git checkout <filename>`)

Die Übungen steigern sich von einfachen Syntaxfehlern zu subtileren Laufzeitproblemen.
Lösungen werden inline diskutiert, aber versuche jede selbst zu lösen, bevor du weiterliest.

#### Bereitschaftscheckliste

Denkst du, du bist bereit einzutauchen?

- [ ] Ich verstehe das Ziel dieses Kurses und seine Voraussetzungen
- [ ] Mein Codespace ist aktiv
- [ ] Ich habe mein Arbeitsverzeichnis entsprechend eingestellt
- [ ] Ich verstehe die Aufgabenstellung

Wenn du alle Kästchen abhaken kannst, kannst du loslegen.

---

## 1. Syntaxfehler

Syntaxfehler sind die häufigsten Fehler, denen du beim Schreiben von Nextflow-Code begegnen wirst. Sie treten auf, wenn der Code nicht den erwarteten Syntaxregeln der Nextflow-DSL entspricht. Diese Fehler verhindern, dass dein Workflow überhaupt läuft, daher ist es wichtig zu lernen, wie man sie schnell identifiziert und behebt.

### 1.1. Fehlende Klammern

Einer der häufigsten Syntaxfehler, und manchmal einer der komplexeren beim Debugging, sind **fehlende oder nicht übereinstimmende Klammern**.

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

**Schlüsselelemente von Syntaxfehlermeldungen:**

- **Datei und Position**: Zeigt, welche Datei und welche Zeile/Spalte den Fehler enthält (`bad_syntax.nf:24:1`)
- **Fehlerbeschreibung**: Erklärt, was der Parser gefunden hat, das er nicht erwartet hat (`Unexpected input: '<EOF>'`)
- **EOF-Indikator**: Die `<EOF>` (End Of File) Nachricht zeigt an, dass der Parser das Ende der Datei erreicht hat, während er noch mehr Inhalt erwartet - ein klassisches Zeichen für nicht geschlossene Klammern

#### Überprüfe den Code

Jetzt schauen wir uns `bad_syntax.nf` an, um zu verstehen, was den Fehler verursacht:

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
// Fehlende schließende Klammer für den Process

workflow {

    // Eingabe-Channel erstellen
    input_ch = channel.of('sample1', 'sample2', 'sample3')

    // Den Process mit dem Eingabe-Channel aufrufen
    PROCESS_FILES(input_ch)
}
```

Für dieses Beispiel haben wir einen Kommentar hinterlassen, der dir zeigt, wo der Fehler ist. Die Nextflow VSCode Extension sollte dir auch einige Hinweise geben, was falsch sein könnte, indem sie die nicht übereinstimmende Klammer rot markiert und das vorzeitige Ende der Datei hervorhebt:

![Bad syntax](img/bad_syntax.png)

**Debugging-Strategie für Klammerfehler:**

1. Verwende VS Codes Klammerabgleich (platziere den Cursor neben einer Klammer)
2. Überprüfe das Problems-Panel auf klammerbezogene Meldungen
3. Stelle sicher, dass jede öffnende `{` eine entsprechende schließende `}` hat

#### Behebe den Code

Ersetze den Kommentar mit der fehlenden schließenden Klammer:

=== "Nach"

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
    }  // Fehlende schließende Klammer hinzufügen

    workflow {

        // Eingabe-Channel erstellen
        input_ch = channel.of('sample1', 'sample2', 'sample3')

        // Den Process mit dem Eingabe-Channel aufrufen
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
    // Fehlende schließende Klammer für den Process

    workflow {

        // Eingabe-Channel erstellen
        input_ch = channel.of('sample1', 'sample2', 'sample3')

        // Den Process mit dem Eingabe-Channel aufrufen
        PROCESS_FILES(input_ch)
    }
    ```

#### Führe die Pipeline aus

Führe den Workflow jetzt erneut aus, um zu bestätigen, dass er funktioniert:

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

### 1.2. Verwendung falscher Process-Schlüsselwörter oder Direktiven

Ein weiterer häufiger Syntaxfehler ist eine **ungültige Process-Definition**. Dies kann passieren, wenn du vergisst, erforderliche Blöcke zu definieren oder falsche Direktiven in der Process-Definition verwendest.

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

Der Fehler zeigt eine "Invalid process definition" an und zeigt den Kontext um das Problem herum. Wenn wir uns die Zeilen 3-7 ansehen, können wir `inputs:` in Zeile 4 sehen, was das Problem ist. Schauen wir uns `invalid_process.nf` an:

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

    // Eingabe-Channel erstellen
    input_ch = channel.of('sample1', 'sample2', 'sample3')

    // Den Process mit dem Eingabe-Channel aufrufen
    PROCESS_FILES(input_ch)
}
```

Wenn wir uns Zeile 4 im Fehlerkontext ansehen, können wir das Problem erkennen: Wir verwenden `inputs` anstelle der korrekten `input`-Direktive. Die Nextflow VSCode Extension wird dies ebenfalls kennzeichnen:

![Invalid process message](img/invalid_process_message.png)

#### Behebe den Code

Ersetze das falsche Schlüsselwort durch das korrekte, indem du [die Dokumentation](https://www.nextflow.io/docs/latest/process.html#) konsultierst:

=== "Nach"

    ```groovy title="invalid_process.nf" hl_lines="4" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:  // Behoben: 'inputs' zu 'input' geändert
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        """
        echo "Processing ${sample_name}" > ${sample_name}_output.txt
        """
    }

    workflow {

        // Eingabe-Channel erstellen
        input_ch = channel.of('sample1', 'sample2', 'sample3')

        // Den Process mit dem Eingabe-Channel aufrufen
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

        // Eingabe-Channel erstellen
        input_ch = channel.of('sample1', 'sample2', 'sample3')

        // Den Process mit dem Eingabe-Channel aufrufen
        PROCESS_FILES(input_ch)
    }
    ```

#### Führe die Pipeline aus

Führe den Workflow jetzt erneut aus, um zu bestätigen, dass er funktioniert:

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

Die Variablennamen, die du in deinen Script-Blöcken verwendest, müssen gültig sein und entweder aus Eingaben oder aus Groovy-Code vor dem Script abgeleitet werden. Aber wenn du zu Beginn der Pipeline-Entwicklung mit Komplexität jonglierst, ist es leicht, Fehler bei der Variablenbenennung zu machen, und Nextflow wird dich schnell darauf hinweisen.

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

Der Fehler wird zur Compile-Zeit abgefangen und zeigt direkt auf die undefinierte Variable in Zeile 17, mit einem Caret, das genau anzeigt, wo das Problem ist.

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
    // Variablen im Groovy-Code vor dem Script definieren
    def output_prefix = "${sample_name}_processed"
    def timestamp = new Date().format("yyyy-MM-dd")

    """
    echo "Processing ${sample_name} on ${timestamp}" > ${output_prefix}.txt
    echo "Using undefined variable: ${undefined_var}" >> ${output_prefix}.txt  // FEHLER: undefined_var nicht definiert
    """
}

workflow {
    input_ch = channel.of('sample1', 'sample2', 'sample3')
    PROCESS_FILES(input_ch)
}
```

Die Fehlermeldung zeigt an, dass die Variable im Script-Template nicht erkannt wird, und dort siehst du es - `${undefined_var}` wird im Script-Block verwendet, aber nicht anderweitig definiert.

#### Behebe den Code

Wenn du einen 'No such variable'-Fehler erhältst, kannst du ihn beheben, indem du entweder die Variable definierst (durch Korrektur von Eingabevariablennamen oder Bearbeitung des Groovy-Codes vor dem Script) oder sie aus dem Script-Block entfernst, wenn sie nicht benötigt wird:

=== "Nach"

    ```groovy title="no_such_var.nf" hl_lines="15-17" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        // Variablen im Groovy-Code vor dem Script definieren
        def output_prefix = "${sample_name}_processed"
        def timestamp = new Date().format("yyyy-MM-dd")

        """
        echo "Processing ${sample_name} on ${timestamp}" > ${output_prefix}.txt
        """  // Die Zeile mit undefined_var entfernt
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
        // Variablen im Groovy-Code vor dem Script definieren
        def output_prefix = "${sample_name}_processed"
        def timestamp = new Date().format("yyyy-MM-dd")

        """
        echo "Processing ${sample_name} on ${timestamp}" > ${output_prefix}.txt
        echo "Using undefined variable: ${undefined_var}" >> ${output_prefix}.txt  // FEHLER: undefined_var nicht definiert
        """
    }

    workflow {
        input_ch = channel.of('sample1', 'sample2', 'sample3')
        PROCESS_FILES(input_ch)
    }
    ```

#### Führe die Pipeline aus

Führe den Workflow jetzt erneut aus, um zu bestätigen, dass er funktioniert:

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

### 1.4. Schlechte Verwendung von Bash-Variablen

Am Anfang mit Nextflow kann es schwierig sein, den Unterschied zwischen Nextflow (Groovy) und Bash-Variablen zu verstehen. Dies kann eine weitere Form des Variablenfehlers erzeugen, der auftritt, wenn man versucht, Variablen im Bash-Inhalt des Script-Blocks zu verwenden.

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
    echo "Processing ${sample_name}" > ${prefix}.txt  # ERROR: ${prefix} is Groovy syntax, not Bash
    """
}
```

In diesem Beispiel definieren wir die `prefix`-Variable in Bash, aber in einem Nextflow-Process wird die `$`-Syntax, die wir verwendet haben, um darauf zu verweisen (`${prefix}`), als Groovy-Variable interpretiert, nicht als Bash. Die Variable existiert nicht im Groovy-Kontext, also erhalten wir einen 'no such variable'-Fehler.

#### Behebe den Code

Wenn du eine Bash-Variable verwenden möchtest, musst du das Dollarzeichen so escapen:

=== "Nach"

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
        echo "Processing ${sample_name}" > \${prefix}.txt  # Fixed: Escaped the dollar sign
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
        echo "Processing ${sample_name}" > ${prefix}.txt  # ERROR: ${prefix} is Groovy syntax, not Bash
        """
    }
    ```

Das sagt Nextflow, dies als Bash-Variable zu interpretieren.

#### Führe die Pipeline aus

Führe den Workflow jetzt erneut aus, um zu bestätigen, dass er funktioniert:

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

    Für einfache Variablenmanipulationen wie String-Verkettung oder Präfix-/Suffix-Operationen ist es normalerweise lesbarer, Groovy-Variablen im Script-Abschnitt zu verwenden anstatt Bash-Variablen im Script-Block:

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

Die Nextflow VSCode Extension hebt Probleme mit der Code-Struktur hervor, die Fehler verursachen werden. Ein häufiges Beispiel ist das Definieren von Channels außerhalb des `workflow {}`-Blocks - dies wird jetzt als Syntaxfehler durchgesetzt.

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

Die Fehlermeldung weist klar auf das Problem hin: Anweisungen (wie Channel-Definitionen) können nicht außerhalb eines Workflow- oder Process-Blocks mit Script-Deklarationen gemischt werden.

#### Überprüfe den Code

Schauen wir uns `badpractice_syntax.nf` an, um zu sehen, was den Fehler verursacht:

```groovy title="badpractice_syntax.nf" hl_lines="3" linenums="1"
#!/usr/bin/env nextflow

input_ch = channel.of('sample1', 'sample2', 'sample3')  // FEHLER: Channel außerhalb des Workflows definiert

process PROCESS_FILES {
    input:
    val sample_name

    output:
    path "${sample_name}_processed.txt"

    script:
    // Variablen im Groovy-Code vor dem Script definieren
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

Die VSCode Extension wird auch die `input_ch`-Variable hervorheben, da sie außerhalb des Workflow-Blocks definiert ist:

![Non-lethal syntax error](img/nonlethal.png)

#### Behebe den Code

Verschiebe die Channel-Definition in den Workflow-Block:

=== "Nach"

    ```groovy title="badpractice_syntax.nf" hl_lines="21" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_processed.txt"

        script:
        // Variablen im Groovy-Code vor dem Script definieren
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

    input_ch = channel.of('sample1', 'sample2', 'sample3')  // FEHLER: Channel außerhalb des Workflows definiert

    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_processed.txt"

        script:
        // Variablen im Groovy-Code vor dem Script definieren
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

Halte deine Input-Channels im Workflow-Block definiert und folge im Allgemeinen allen anderen Empfehlungen, die die Extension macht.

### Zusammenfassung

Du kannst Syntaxfehler systematisch identifizieren und beheben, indem du Nextflow-Fehlermeldungen und IDE-visuelle Indikatoren verwendest. Häufige Syntaxfehler umfassen fehlende Klammern, falsche Process-Schlüsselwörter, undefinierte Variablen und unsachgemäße Verwendung von Bash- vs. Nextflow-Variablen. Die VSCode Extension hilft, viele davon vor der Laufzeit zu erkennen. Mit diesen Syntax-Debugging-Fähigkeiten kannst du die häufigsten Nextflow-Syntaxfehler schnell beheben und dich komplexeren Laufzeitproblemen widmen.

### Was kommt als Nächstes?

Lerne, komplexere Channel-Strukturfehler zu debuggen, die auftreten, selbst wenn die Syntax korrekt ist.

---

## 2. Channel-Strukturfehler

Channel-Strukturfehler sind subtiler als Syntaxfehler, weil der Code syntaktisch korrekt ist, aber die Datenformen nicht zu dem passen, was Processes erwarten. Nextflow wird versuchen, die Pipeline auszuführen, könnte aber feststellen, dass die Anzahl der Eingaben nicht übereinstimmt und fehlschlagen. Diese Fehler treten typischerweise nur zur Laufzeit auf und erfordern ein Verständnis der Daten, die durch deinen Workflow fließen.

!!! tip "Debugging von Channels mit `.view()`"

    Denke während dieses Abschnitts daran, dass du den `.view()`-Operator verwenden kannst, um Channel-Inhalte an jedem Punkt in deinem Workflow zu inspizieren. Dies ist eines der mächtigsten Debugging-Tools zum Verstehen von Channel-Strukturproblemen. Wir werden diese Technik in Abschnitt 2.4 im Detail erkunden, aber fühle dich frei, sie zu verwenden, während du durch die Beispiele arbeitest.

    ```groovy
    my_channel.view()  // Shows what's flowing through the channel
    ```

### 2.1. Falsche Anzahl von Input-Channels

Dieser Fehler tritt auf, wenn du eine andere Anzahl von Channels übergibst, als ein Process erwartet.

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

Die Fehlermeldung gibt klar an, dass der Aufruf 1 Argument erwartet, aber 2 erhalten hat, und zeigt auf Zeile 23. Schauen wir uns `bad_number_inputs.nf` an:

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

    // Zwei separate Channels erstellen
    samples_ch = channel.of('sample1', 'sample2', 'sample3')
    files_ch = channel.of('file1.txt', 'file2.txt', 'file3.txt')

    // FEHLER: 2 Channels übergeben, aber Process erwartet nur 1
    PROCESS_FILES(samples_ch, files_ch)
}
```

Du solltest den nicht übereinstimmenden `PROCESS_FILES`-Aufruf sehen, der mehrere Input-Channels bereitstellt, wenn der Process nur einen definiert. Die VSCode Extension wird auch den Process-Aufruf rot unterstreichen und eine Diagnosemeldung liefern, wenn du mit der Maus darüberfährst:

![Incorrect number of args message](img/incorrect_num_args.png)

#### Behebe den Code

Für dieses spezifische Beispiel erwartet der Process einen einzelnen Channel und benötigt den zweiten Channel nicht, also können wir es beheben, indem wir nur den `samples_ch`-Channel übergeben:

=== "Nach"

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

        // Zwei separate Channels erstellen
        samples_ch = channel.of('sample1', 'sample2', 'sample3')
        files_ch = channel.of('file1.txt', 'file2.txt', 'file3.txt')

        // Behoben: Nur den Channel übergeben, den der Process erwartet
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

        // Zwei separate Channels erstellen
        samples_ch = channel.of('sample1', 'sample2', 'sample3')
        files_ch = channel.of('file1.txt', 'file2.txt', 'file3.txt')

        // FEHLER: 2 Channels übergeben, aber Process erwartet nur 1
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

Häufiger als in diesem Beispiel könntest du zusätzliche Eingaben zu einem Process hinzufügen und vergessen, den Workflow-Aufruf entsprechend zu aktualisieren, was zu dieser Art von Fehler führen kann. Glücklicherweise ist dies einer der leichter zu verstehenden und zu behebenden Fehler, da die Fehlermeldung recht klar über die Diskrepanz ist.

### 2.2. Channel-Erschöpfung (Process läuft weniger oft als erwartet)

Einige Channel-Strukturfehler sind viel subtiler und erzeugen überhaupt keine Fehler. Wahrscheinlich am häufigsten davon reflektiert eine Herausforderung, der neue Nextflow-Nutzer\*innen beim Verständnis gegenüberstehen, dass Queue-Channels erschöpft werden und keine Items mehr haben können, was bedeutet, dass der Workflow vorzeitig endet.

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
    // Variablen im Groovy-Code vor dem Script definieren
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

Der Process läuft nur einmal statt dreimal, weil der `reference_ch`-Channel ein Queue-Channel ist, der nach der ersten Process-Ausführung erschöpft ist. Wenn ein Channel erschöpft ist, stoppt der gesamte Process, selbst wenn andere Channels noch Items haben.

Dies ist ein häufiges Muster, bei dem du eine einzelne Referenzdatei hast, die über mehrere Proben wiederverwendet werden muss. Die Lösung besteht darin, den Referenz-Channel in einen Value-Channel umzuwandeln, der unbegrenzt wiederverwendet werden kann.

#### Behebe den Code

Es gibt ein paar Möglichkeiten, dies zu adressieren, abhängig davon, wie viele Dateien betroffen sind.

**Option 1**: Du hast eine einzelne Referenzdatei, die du häufig wiederverwendest. Du kannst einfach einen Value-Channel-Typ erstellen, der immer wieder verwendet werden kann. Es gibt drei Wege, dies zu tun:

**1a** Verwende `channel.value()`:

```groovy title="exhausted.nf (behoben - Option 1a)" hl_lines="2" linenums="21"
workflow {
    reference_ch = channel.value('baseline_reference')  // Value channel can be reused
    input_ch = channel.of('sample1', 'sample2', 'sample3')

    PROCESS_FILES(reference_ch, input_ch)
}
```

**1b** Verwende den `first()` [Operator](https://www.nextflow.io/docs/latest/reference/operator.html#first):

```groovy title="exhausted.nf (behoben - Option 1b)" hl_lines="2" linenums="21"
workflow {
    reference_ch = channel.of('baseline_reference').first()  // Convert to value channel
    input_ch = channel.of('sample1', 'sample2', 'sample3')

    PROCESS_FILES(reference_ch, input_ch)
}
```

**1c.** Verwende den `collect()` [Operator](https://www.nextflow.io/docs/latest/reference/operator.html#collect):

```groovy title="exhausted.nf (behoben - Option 1c)" hl_lines="2" linenums="21"
workflow {
    reference_ch = channel.of('baseline_reference').collect()  // Convert to value channel
    input_ch = channel.of('sample1', 'sample2', 'sample3')

    PROCESS_FILES(reference_ch, input_ch)
}
```

**Option 2**: In komplexeren Szenarien, vielleicht wo du mehrere Referenzdateien für alle Proben im Proben-Channel hast, kannst du den `combine`-Operator verwenden, um einen neuen Channel zu erstellen, der die beiden Channels zu Tupeln kombiniert:

```groovy title="exhausted.nf (behoben - Option 2)" hl_lines="4" linenums="21"
workflow {
    reference_ch = channel.of('baseline_reference','other_reference')
    input_ch = channel.of('sample1', 'sample2', 'sample3')
    combined_ch = reference_ch.combine(input_ch)  // Creates cartesian product

    PROCESS_FILES(combined_ch)
}
```

Der `.combine()`-Operator erzeugt ein kartesisches Produkt der beiden Channels, sodass jedes Item in `reference_ch` mit jedem Item in `input_ch` gepaart wird. Dies ermöglicht es dem Process, für jede Probe zu laufen und dabei trotzdem die Referenz zu verwenden.

Dies erfordert, dass die Process-Eingabe angepasst wird. In unserem Beispiel müsste der Beginn der Process-Definition wie folgt angepasst werden:

```groovy title="exhausted.nf (behoben - Option 2)" hl_lines="5" linenums="1"
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

Du solltest jetzt sehen, dass alle drei Proben verarbeitet werden, anstatt nur eine.

### 2.3. Falsche Channel-Inhaltsstruktur

Wenn Workflows ein gewisses Maß an Komplexität erreichen, kann es etwas schwierig sein, den Überblick über die internen Strukturen jedes Channels zu behalten, und Menschen erzeugen häufig Diskrepanzen zwischen dem, was der Process erwartet, und dem, was der Channel tatsächlich enthält. Dies ist subtiler als das Problem, das wir früher besprochen haben, wo die Anzahl der Channels falsch war. In diesem Fall kannst du die richtige Anzahl von Input-Channels haben, aber die interne Struktur eines oder mehrerer dieser Channels stimmt nicht mit dem überein, was der Process erwartet.

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

Die eckigen Klammern in der Fehlermeldung liefern hier den Hinweis - der Process behandelt das Tupel als einen einzelnen Wert, was nicht das ist, was wir wollen. Schauen wir uns `bad_channel_shape.nf` an:

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

Du kannst sehen, dass wir einen Channel erstellen, der aus Tupeln besteht: `['sample1', 'file1.txt']`, aber der Process einen einzelnen Wert erwartet, `val sample_name`. Der ausgeführte Befehl zeigt, dass der Process versucht, eine Datei mit dem Namen `[sample3, file3.txt]_output.txt` zu erstellen, was nicht die beabsichtigte Ausgabe ist.

#### Behebe den Code

Um dies zu beheben, könnten wir, wenn der Process beide Eingaben benötigt, den Process so anpassen, dass er ein Tupel akzeptiert:

=== "Option 1: Tupel im Process akzeptieren"

    === "Nach"

        ```groovy title="bad_channel_shape.nf" hl_lines="5"  linenums="1"
        #!/usr/bin/env nextflow

        process PROCESS_FILES {
            input:
                tuple val(sample_name), val(file_name)  // Behoben: Tuple akzeptieren

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

    === "Nach"

        ```groovy title="bad_channel_shape.nf" hl_lines="9" linenums="16"
        workflow {

            // Channel emits tuples, but process expects single values
            input_ch = channel.of(
              ['sample1', 'file1.txt'],
              ['sample2', 'file2.txt'],
              ['sample3', 'file3.txt']
            )
            PROCESS_FILES(input_ch.map { it[0] })  // Behoben: Erstes Element extrahieren
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

### 2.4. Channel-Debugging-Techniken

#### Verwendung von `.view()` zur Channel-Inspektion

Das mächtigste Debugging-Tool für Channels ist der `.view()`-Operator. Mit `.view()` kannst du die Form deiner Channels in allen Phasen verstehen, um beim Debugging zu helfen.

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

Um zu vermeiden, dass du in Zukunft übermäßig `.view()`-Operationen verwenden musst, um Channel-Inhalte zu verstehen, ist es ratsam, einige Kommentare hinzuzufügen:

```groovy title="bad_channel_shape_viewed.nf (mit Kommentaren)" linenums="16" hl_lines="8 9"
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

Dies wird wichtiger, wenn deine Workflows an Komplexität zunehmen und die Channel-Struktur undurchsichtiger wird.

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

### Zusammenfassung

Viele Channel-Strukturfehler können mit gültiger Nextflow-Syntax erstellt werden. Du kannst Channel-Strukturfehler debuggen, indem du den Datenfluss verstehst, `.view()`-Operatoren zur Inspektion verwendest und Fehlermuster wie eckige Klammern erkennst, die auf unerwartete Tupelstrukturen hinweisen.

### Was kommt als Nächstes?

Lerne über Fehler, die durch Process-Definitionen erzeugt werden.

---

## 3. Process-Strukturfehler

Die meisten Fehler, denen du im Zusammenhang mit Processes begegnest, werden sich auf Fehler beziehen, die du bei der Formung des Befehls gemacht hast, oder auf Probleme im Zusammenhang mit der zugrunde liegenden Software. Allerdings kannst du, ähnlich wie bei den Channel-Problemen oben, Fehler in der Process-Definition machen, die nicht als Syntaxfehler qualifizieren, aber die zur Laufzeit Fehler verursachen werden.

### 3.1. Fehlende Ausgabedateien

Ein häufiger Fehler beim Schreiben von Processes ist es, etwas zu tun, das eine Diskrepanz zwischen dem, was der Process erwartet, und dem, was generiert wird, erzeugt.

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

Die Fehlermeldung zeigt an, dass der Process erwartet hat, eine Ausgabedatei namens `sample3.txt` zu erzeugen, aber das Script tatsächlich `sample3_output.txt` erstellt. Schauen wir uns die Process-Definition in `missing_output.nf` an:

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

Du solltest sehen, dass es eine Diskrepanz zwischen dem Ausgabedateinamen im `output:`-Block und dem im Script verwendeten gibt. Diese Diskrepanz führt dazu, dass der Process fehlschlägt. Wenn du auf diese Art von Fehler stößt, gehe zurück und überprüfe, dass die Ausgaben zwischen deiner Process-Definition und deinem Output-Block übereinstimmen.

Wenn das Problem immer noch nicht klar ist, überprüfe das Work-Verzeichnis selbst, um die tatsächlich erstellten Ausgabedateien zu identifizieren:

```bash
❯ ls -h work/02/9604d49fb8200a74d737c72a6c98ed
sample3_output.txt
```

Für dieses Beispiel würde uns dies verdeutlichen, dass ein `_output`-Suffix in den Ausgabedateinamen einbezogen wird, entgegen unserer `output:`-Definition.

#### Behebe den Code

Behebe die Diskrepanz, indem du den Ausgabedateinamen konsistent machst:

=== "Nach"

    ```groovy title="missing_output.nf" hl_lines="6 10" linenums="3"
    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_output.txt"  // Behoben: Mit der Script-Ausgabe übereinstimmen

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
