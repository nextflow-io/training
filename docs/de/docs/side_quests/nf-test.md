# Testen mit nf-test

Die Möglichkeit, systematisch zu testen, dass jeder Teil deines Workflows das tut, was er soll, ist entscheidend für Reproduzierbarkeit und langfristige Wartbarkeit und kann während des Entwicklungsprozesses eine große Hilfe sein.

Lass uns kurz darüber sprechen, warum Testen so wichtig ist. Wenn du einen Workflow entwickelst, ist eines der ersten Dinge, die du tust, dir Testdaten zu besorgen, von denen du weißt, dass sie gültig sind und ein Ergebnis liefern sollten. Du fügst den ersten Prozess zur Pipeline hinzu und verbindest ihn mit deinen Eingaben, damit er funktioniert. Um zu prüfen, ob alles funktioniert, führst du ihn mit den Testdaten aus. Wenn das funktioniert, gehst du zum nächsten Prozess über und führst die Testdaten erneut aus. Du wiederholst diesen Prozess, bis du eine Pipeline hast, mit der du zufrieden bist.

Dann fügst du vielleicht einen einfachen wahr/falsch-Parameter wie `--skip_process` hinzu. Jetzt musst du die Pipeline zweimal ausführen, einmal mit jedem Parameter, um sicherzustellen, dass sie wie erwartet funktioniert. Aber warte, wie prüfen wir, ob `--skip_process` den Prozess tatsächlich überspringt? Wir müssen in den Ausgaben nachsehen oder die Log-Dateien überprüfen! Das ist mühsam und fehleranfällig.

Während du deine Pipeline entwickelst, wird sie schnell so komplex, dass manuelles Testen jeder Iteration langsam und fehleranfällig ist. Wenn du einen Fehler findest, wird es außerdem sehr schwierig sein, genau zu bestimmen, woher in deiner Pipeline der Fehler kommt. Hier kommt das Testen ins Spiel.

Testen ermöglicht es dir, systematisch zu überprüfen, dass jeder Teil deiner Pipeline wie erwartet funktioniert. Die Vorteile gut geschriebener Tests für Entwickler\*innen sind enorm:

- **Vertrauen**: Da die Tests die gesamte Pipeline abdecken, kannst du sicher sein, dass eine Änderung nichts anderes beeinträchtigt
- **Zuverlässigkeit**: Wenn mehrere Entwickler\*innen an der Pipeline arbeiten, wissen sie, dass die anderen die Pipeline und jede Komponente nicht kaputt gemacht haben
- **Transparenz**: Die Tests zeigen, wo eine Pipeline fehlschlägt, und erleichtern es, das Problem aufzuspüren. Sie fungieren auch als eine Form der Dokumentation und zeigen, wie ein Prozess oder Workflow ausgeführt wird
- **Geschwindigkeit**: Da die Tests automatisiert sind, können sie sehr schnell und wiederholt ausgeführt werden. Du kannst schnell iterieren mit weniger Angst, neue Fehler einzuführen

Es gibt viele verschiedene Arten von Tests, die wir schreiben können:

1. **Modul-Tests**: Für einzelne Prozesse
2. **Workflow-Tests**: Für einen einzelnen Workflow
3. **Pipeline-Tests**: Für die Pipeline als Ganzes
4. **Performance-Tests**: Für die Geschwindigkeit und Effizienz der Pipeline
5. **Stress-Tests**: Bewertung der Pipeline-Leistung unter extremen Bedingungen, um ihre Grenzen zu bestimmen

Das Testen einzelner Prozesse ist analog zu Unit-Tests in anderen Sprachen. Das Testen des Workflows oder der gesamten Pipeline ist analog zu sogenannten Integrationstests in anderen Sprachen, bei denen wir die Interaktionen der Komponenten testen.

[**nf-test**](https://www.nf-test.com/) ist ein Tool, mit dem du Modul-, Workflow- und Pipeline-Tests schreiben kannst. Kurz gesagt, ermöglicht es dir, systematisch zu überprüfen, dass jeder einzelne Teil der Pipeline wie erwartet funktioniert, _isoliert_.

### Lernziele

In dieser Side Quest lernst du, nf-test zu verwenden, um einen Workflow-Test für die Pipeline sowie Modul-Tests für die drei Prozesse zu schreiben, die sie aufruft.

Am Ende dieser Side Quest wirst du die folgenden Techniken effektiv einsetzen können:

- nf-test in deinem Projekt initialisieren
- Modul- und Workflow-Tests generieren
- Gängige Arten von Assertions hinzufügen
- Verstehen, wann Snapshots vs. Content-Assertions verwendet werden
- Tests für ein gesamtes Projekt ausführen

Diese Fähigkeiten helfen dir, eine umfassende Teststrategie in deinen Pipeline-Projekten zu implementieren und sicherzustellen, dass sie robuster und wartbarer sind.

### Voraussetzungen

Bevor du diese Side Quest angehst, solltest du:

- Das [Hello Nextflow](../hello_nextflow/README.md) Tutorial oder einen gleichwertigen Einsteigerkurs abgeschlossen haben
- Mit grundlegenden Nextflow-Konzepten und -Mechanismen vertraut sein (Prozesse, Kanäle, Operatoren, Arbeiten mit Dateien, Metadaten)

---

## 0. Erste Schritte

#### Öffne den Training-Codespace

Falls du es noch nicht getan hast, stelle sicher, dass du die Trainingsumgebung wie in [Umgebung einrichten](../envsetup/index.md) beschrieben öffnest.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Wechsle in das Projektverzeichnis

Lass uns in das Verzeichnis wechseln, in dem sich die Dateien für dieses Tutorial befinden.

```bash
cd side-quests/nf-test
```

Du kannst VSCode so einstellen, dass es sich auf dieses Verzeichnis konzentriert:

```bash
code .
```

#### Überprüfe die Materialien

Du findest eine Haupt-Workflow-Datei und eine CSV-Datei namens `greetings.csv`, die die Eingabe für die Pipeline enthält.

```console title="Directory contents"
.
├── greetings.csv
└── main.nf
```

Eine detaillierte Beschreibung der Dateien findest du im [Warmup von Hello Nextflow](../hello_nextflow/00_orientation.md).

Der Workflow, den wir testen werden, ist eine Teilmenge des Hello-Workflows, der in [Hello Workflow](../hello_nextflow/03_hello_workflow.md) erstellt wurde.

??? example "Was macht der Hello Nextflow Workflow?"

    Falls du das [Hello Nextflow](../hello_nextflow/index.md) Training nicht gemacht hast, hier ist ein kurzer Überblick darüber, was dieser einfache Workflow macht.

    Der Workflow nimmt eine CSV-Datei mit Begrüßungen, führt vier aufeinanderfolgende Transformationsschritte durch und gibt eine einzelne Textdatei aus, die ein ASCII-Bild eines lustigen Charakters enthält, der die Begrüßungen sagt.

    Die vier Schritte sind als Nextflow-Prozesse (`sayHello`, `convertToUpper`, `collectGreetings` und `cowpy`) implementiert, die in separaten Moduldateien gespeichert sind.

    1. **`sayHello`:** Schreibt jede Begrüßung in ihre eigene Ausgabedatei (z.B. "Hello-output.txt")
    2. **`convertToUpper`:** Konvertiert jede Begrüßung in Großbuchstaben (z.B. "HELLO")
    3. **`collectGreetings`:** Sammelt alle Großbuchstaben-Begrüßungen in einer einzigen Batch-Datei
    4. **`cowpy`:** Generiert ASCII-Art mit dem `cowpy`-Tool

    Die Ergebnisse werden in einem Verzeichnis namens `results/` veröffentlicht, und die endgültige Ausgabe der Pipeline (wenn sie mit Standardparametern ausgeführt wird) ist eine Klartextdatei mit ASCII-Art eines Charakters, der die großgeschriebenen Begrüßungen sagt.

    In dieser Side Quest verwenden wir eine Zwischenform des Hello-Workflows, die nur die ersten beiden Prozesse enthält.

Die Teilmenge, mit der wir arbeiten werden, besteht aus zwei Prozessen: `sayHello` und `convertToUpper`.
Den vollständigen Workflow-Code siehst du unten.

??? example "Workflow-Code"

    ```groovy title="main.nf"
    /*
    * Pipeline-Parameter
    */
    params.input_file = "greetings.csv"

    /*
    * Verwende echo, um 'Hello World!' auf die Standardausgabe zu drucken
    */
    process sayHello {

        publishDir 'results', mode: 'copy'

        input:
            val greeting

        output:
            path "${greeting}-output.txt"

        script:
        """
        echo '$greeting' > '$greeting-output.txt'
        """
    }

    /*
    * Verwende ein Text-Ersetzungs-Utility, um die Begrüßung in Großbuchstaben zu konvertieren
    */
    process convertToUpper {

        publishDir 'results', mode: 'copy'

        input:
            path input_file

        output:
            path "UPPER-${input_file}"

        script:
        """
        cat '$input_file' | tr '[a-z]' '[A-Z]' > UPPER-${input_file}
        """
    }

    workflow {

        // create a channel for inputs from a CSV file
        greeting_ch = channel.fromPath(params.input_file).splitCsv().flatten()

        // emit a greeting
        sayHello(greeting_ch)

        // convert the greeting to uppercase
        convertToUpper(sayHello.out)
    }
    ```

#### Führe den Workflow aus

Lass uns den Workflow ausführen, um sicherzustellen, dass er wie erwartet funktioniert.

```bash
nextflow run main.nf
```

```console title="Result of running the workflow"
 N E X T F L O W   ~  version 24.10.2

Launching `main.nf` [soggy_linnaeus] DSL2 - revision: bbf79d5c31

executor >  local (6)
[f7/c3be66] sayHello (3)       | 3 of 3 ✔
[cd/e15303] convertToUpper (3) | 3 of 3 ✔
```

GLÜCKWUNSCH! Du hast gerade einen Test ausgeführt!

"Warte, was? Ich habe nur den Workflow ausgeführt und er hat funktioniert! Wie ist das ein Test?"

Gute Frage!

Lass uns aufschlüsseln, was gerade passiert ist.

Du hast den Workflow mit den Standardparametern ausgeführt, bestätigt, dass er funktioniert hat, und bist mit den Ergebnissen zufrieden. Das ist die Essenz des Testens. Wenn du das Hello Nextflow Training durchgearbeitet hast, wirst du bemerkt haben, dass wir jeden Abschnitt immer damit begonnen haben, den Workflow auszuführen, den wir als Ausgangspunkt verwendeten, um zu bestätigen, dass alles korrekt eingerichtet ist.

Software-Tests machen im Wesentlichen diesen Prozess für uns.

#### Überprüfe die Aufgabe

Deine Herausforderung besteht darin, diesem Workflow standardisierte Tests mit nf-test hinzuzufügen, um es einfach zu machen, zu überprüfen, dass jeder Teil weiterhin wie erwartet funktioniert, falls weitere Änderungen vorgenommen werden.

#### Bereitschafts-Checkliste

Denkst du, du bist bereit einzutauchen?

- [ ] Ich verstehe das Ziel dieses Kurses und seine Voraussetzungen
- [ ] Mein Codespace läuft
- [ ] Ich habe mein Arbeitsverzeichnis entsprechend eingestellt
- [ ] Ich habe den Workflow erfolgreich ausgeführt
- [ ] Ich verstehe die Aufgabe

Wenn du alle Kästchen ankreuzen kannst, bist du startklar.

---

## 1. `nf-test` initialisieren

Das `nf-test`-Paket bietet einen Initialisierungsbefehl, der einige Dinge einrichtet, damit wir mit der Entwicklung von Tests für unser Projekt beginnen können.

```bash
nf-test init
```

Dies sollte die folgende Ausgabe erzeugen:

```bash
🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr

Project configured. Configuration is stored in nf-test.config
```

Es erstellt auch ein `tests`-Verzeichnis, das einen Konfigurationsdatei-Stub enthält.

### 1.1. Generiere einen nf-test

`nf-test` kommt mit einer Reihe von Tools zum Erstellen von nf-test-Dateien, die uns den Großteil der Arbeit ersparen. Diese fallen unter den Unterbefehl `generate`. Lass uns einen Test für die Pipeline generieren:

```bash
nf-test generate pipeline main.nf
```

```console title="Output"
> nf-test generate pipeline main.nf

Load source file '/workspaces/training/side-quests/nf-test/main.nf'
Wrote pipeline test file '/workspaces/training/side-quests/nf-test/tests/main.nf.test

SUCCESS: Generated 1 test files.
```

Dies erstellt eine `main.nf.test`-Datei im `tests`-Verzeichnis. Dies ist unsere Pipeline-Test-Datei. Wenn du `tree tests/` ausführst, solltest du so etwas sehen:

```console title="Test directory contents"
tests/
├── main.nf.test
└── nextflow.config
```

Die `main.nf.test`-Datei ist unsere Pipeline-Test-Datei. Lass uns sie öffnen und den Inhalt ansehen.

```groovy title="tests/main.nf.test"
nextflow_pipeline {

    name "Test Workflow main.nf"
    script "main.nf"

    test("Should run without failures") {

        when {
            params {
                // define parameters here. Example:
                // outdir = "tests/results"
            }
        }

        then {
            assert workflow.success
        }

    }

}
```

Wir nehmen uns einen Moment Zeit, um die Struktur der Testdatei zu verstehen.

Der `nextflow_pipeline`-Block ist der Einstiegspunkt für alle Pipeline-Tests. Er enthält Folgendes:

- `name`: Der Name des Tests
- `script`: Der Pfad zum Pipeline-Skript

Der `test`-Block ist der eigentliche Test. Er enthält Folgendes:

- `when`: Die Bedingungen, unter denen der Test ausgeführt werden soll. Dies umfasst die Parameter, die zum Ausführen der Pipeline verwendet werden
- `then`: Die Assertions, die gemacht werden sollen. Dies umfasst die erwarteten Ergebnisse der Pipeline

In einfachem Deutsch liest sich die Logik des Tests wie folgt:
"**Wenn** diese _Parameter_ dieser _Pipeline_ bereitgestellt werden, **dann** erwarten wir diese Ergebnisse zu sehen."

Dies ist kein funktionaler Test, wir werden im nächsten Abschnitt demonstrieren, wie man ihn in einen solchen verwandelt.

### Eine Anmerkung zu Testnamen

Im obigen Beispiel haben wir den Standardnamen "Should run without failures" verwendet, der für einen grundlegenden Test geeignet ist, der nur prüft, ob die Pipeline erfolgreich läuft. Wenn wir jedoch spezifischere Testfälle hinzufügen, sollten wir aussagekräftigere Namen verwenden, die angeben, was wir tatsächlich testen. Zum Beispiel:

- "Should convert input to uppercase" - beim Testen spezifischer Funktionalität
- "Should handle empty input gracefully" - beim Testen von Grenzfällen
- "Should respect max memory parameter" - beim Testen von Ressourcenbeschränkungen
- "Should create expected output files" - beim Testen der Dateigenerierung

Gute Testnamen sollten:

1. Mit "Should" beginnen, um klarzustellen, was das erwartete Verhalten ist
2. Die spezifische Funktionalität oder das Szenario beschreiben, das getestet wird
3. Klar genug sein, dass du bei einem Testfehler weißt, welche Funktionalität kaputt ist

Wenn wir später weitere Assertions und spezifische Testfälle hinzufügen, werden wir diese aussagekräftigeren Namen verwenden, um klarzustellen, was jeder Test überprüft.

### 1.2. Führe den Test aus

Lass uns den Test ausführen, um zu sehen, was passiert.

```bash
nf-test test tests/main.nf.test
```

```console title="nf-test pipeline fail"
> nf-test test tests/main.nf.test

🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Workflow main.nf

  Test [693ba951] 'Should run without failures' FAILED (4.652s)

  Assertion failed:

  assert workflow.success
         |        |
         workflow false

  Nextflow stdout:

  ERROR ~ No such file or directory: /workspaces/training/side-quests/nf-test/.nf-test/tests/693ba951a20fec36a5a9292ed1cc8a9f/greetings.csv

   -- Check '/workspaces/training/side-quests/nf-test/.nf-test/tests/693ba951a20fec36a5a9292ed1cc8a9f/meta/nextflow.log' file for details
  Nextflow stderr:

FAILURE: Executed 1 tests in 4.679s (1 failed)
```

Der Test schlägt fehl! Was ist passiert?

1. nf-test hat versucht, die Pipeline wie sie ist auszuführen, unter Verwendung der Einstellungen im `when`-Block:

```groovy title="tests/main.nf.test"
when {
    params {
        // define parameters here. Example:
        // outdir = "tests/results"
    }
}
```

2. nf-test hat den Status der Pipeline überprüft und ihn mit dem `when`-Block verglichen:

```groovy title="tests/main.nf.test"
then {
    assert workflow.success
}
```

Beachte, wie nf-test gemeldet hat, dass die Pipeline fehlgeschlagen ist, und die Fehlermeldung von Nextflow bereitgestellt hat:

```console title="Error"
ERROR ~ No such file or directory: /workspaces/training/side-quests/nf-test/.nf-test/tests/693ba951a20fec36a5a9292ed1cc8a9f/greetings.csv
```

Was war also das Problem? Denk daran, dass die Pipeline eine greetings.csv-Datei im Projektverzeichnis hat. Wenn nf-test die Pipeline ausführt, wird es nach dieser Datei suchen, kann sie aber nicht finden. Die Datei ist da, was passiert? Nun, wenn wir uns den Pfad ansehen, können wir sehen, dass der Test im Pfad `./nf-test/tests/longHashString/` stattfindet. Genau wie Nextflow erstellt nf-test ein neues Verzeichnis für jeden Test, um alles isoliert zu halten. Die Datendatei befindet sich nicht dort, also müssen wir den Pfad zur Datei im ursprünglichen Test korrigieren.

Lass uns zur Testdatei zurückgehen und den Pfad zur Datei im `when`-Block ändern.

Du fragst dich vielleicht, wie wir im Test auf das Stammverzeichnis der Pipeline zeigen werden. Da dies eine häufige Situation ist, hat nf-test eine Reihe globaler Variablen, die wir verwenden können, um uns das Leben zu erleichtern. Die vollständige Liste findest du [hier](https://www.nf-test.com/docs/testcases/global_variables/), aber in der Zwischenzeit verwenden wir die Variable `projectDir`, die das Stammverzeichnis des Pipeline-Projekts bedeutet.

_Vorher:_

```groovy title="tests/main.nf.test" linenums="1" hl_lines="3 4"
when {
    params {
        // define parameters here. Example:
        // outdir = "tests/results"
    }
}
```

_Nachher:_

```groovy title="tests/main.nf.test" linenums="1" hl_lines="3"
when {
    params {
        input_file = "${projectDir}/greetings.csv"
    }
}
```

Lass uns den Test erneut ausführen, um zu sehen, ob er funktioniert.

```bash title="nf-test pipeline pass"
nf-test test tests/main.nf.test
```

```console title="Pipeline passes"
> nf-test test tests/main.nf.test

🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Workflow main.nf

  Test [1d4aaf12] 'Should run without failures' PASSED (1.619s)


SUCCESS: Executed 1 tests in 1.626s
```

Erfolg! Die Pipeline läuft erfolgreich und der Test ist bestanden. Führe ihn so oft aus, wie du möchtest, und du wirst immer das gleiche Ergebnis erhalten!

Standardmäßig ist die Nextflow-Ausgabe verborgen, aber um dich selbst zu überzeugen, dass nf-test definitiv den Workflow ausführt, kannst du das `--verbose`-Flag verwenden:

```bash
nf-test test tests/main.nf.test --verbose
```

```console title="Pipeline runs all processes"
> nf-test test tests/main.nf.test

🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Workflow main.nf

  Test [693ba951] 'Should run without failures'
    > Nextflow 24.10.4 is available - Please consider updating your version to it
    > N E X T F L O W  ~  version 24.10.0
    > Launching `/workspaces/training/side-quests/nf-test/main.nf` [zen_ampere] DSL2 - revision: bbf79d5c31
    > [2b/61e453] Submitted process > sayHello (2)
    > [31/4e1606] Submitted process > sayHello (1)
    > [bb/5209ee] Submitted process > sayHello (3)
    > [83/83db6f] Submitted process > convertToUpper (2)
    > [9b/3428b1] Submitted process > convertToUpper (1)
    > [ca/0ba51b] Submitted process > convertToUpper (3)
    PASSED (5.206s)


SUCCESS: Executed 1 tests in 5.239s
```

### 1.3. Füge Assertions hinzu

Eine einfache Überprüfung besteht darin, sicherzustellen, dass unsere Pipeline alle erwarteten Prozesse ausführt und keine stillschweigend überspringt. Denk daran, dass unsere Pipeline 6 Prozesse ausführt, einen namens `sayHello` und einen namens `convertToUpper` für jede der 3 Begrüßungen.

Lass uns unserem Test eine Assertion hinzufügen, um zu überprüfen, dass die Pipeline die erwartete Anzahl von Prozessen ausführt. Wir werden auch unseren Testnamen aktualisieren, um besser widerzuspiegeln, was wir testen.

**Vorher:**

```groovy title="tests/main.nf.test" linenums="1" hl_lines="1"
    test("Should run without failures") {

        when {
            params {
                input_file = "${projectDir}/greetings.csv"
            }
        }

        then {
            assert workflow.success
        }

    }
```

**Nachher:**

```groovy title="tests/main.nf.test" linenums="1" hl_lines="1 11"
    test("Should run successfully with correct number of processes") {

        when {
            params {
                input_file = "${projectDir}/greetings.csv"
            }
        }

        then {
            assert workflow.success
            assert workflow.trace.tasks().size() == 6
        }

    }
```

Der Testname spiegelt jetzt besser wider, was wir tatsächlich überprüfen - nicht nur, dass die Pipeline ohne Fehler läuft, sondern dass sie die erwartete Anzahl von Prozessen ausführt.

Lass uns den Test erneut ausführen, um zu sehen, ob er funktioniert.

```bash title="nf-test pipeline pass"
nf-test test tests/main.nf.test
```

```console title="Pipeline passes with assertions"
🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Workflow main.nf

  Test [1d4aaf12] 'Should run successfully with correct number of processes' PASSED (1.567s)


SUCCESS: Executed 1 tests in 1.588s
```

Erfolg! Die Pipeline läuft erfolgreich und der Test ist bestanden. Jetzt haben wir begonnen, die Details der Pipeline zu testen, sowie den Gesamtstatus.

### 1.4. Teste die Ausgabe

Lass uns unserem Test eine Assertion hinzufügen, um zu überprüfen, dass die Ausgabedatei erstellt wurde. Wir fügen sie als separaten Test mit einem informativen Namen hinzu, um die Ergebnisse leichter interpretierbar zu machen.

**Vorher:**

```groovy title="tests/main.nf.test" linenums="1" hl_lines="14"
    test("Should run successfully with correct number of processes") {

        when {
            params {
                input_file = "${projectDir}/greetings.csv"
            }
        }

        then {
            assert workflow.success
            assert workflow.trace.tasks().size() == 6
        }

    }
```

**Nachher:**

```groovy title="tests/main.nf.test" linenums="1" hl_lines="14-33"
    test("Should run successfully with correct number of processes") {

        when {
            params {
                input_file = "${projectDir}/greetings.csv"
            }
        }

        then {
            assert workflow.success
            assert workflow.trace.tasks().size() == 6
        }

    }

    test("Should produce correct output files") {

        when {
            params {
                input_file = "${projectDir}/greetings.csv"
            }
        }

        then {
            assert file("$launchDir/results/Bonjour-output.txt").exists()
            assert file("$launchDir/results/Hello-output.txt").exists()
            assert file("$launchDir/results/Holà-output.txt").exists()
            assert file("$launchDir/results/UPPER-Bonjour-output.txt").exists()
            assert file("$launchDir/results/UPPER-Hello-output.txt").exists()
            assert file("$launchDir/results/UPPER-Holà-output.txt").exists()
        }

    }
```

Führe den Test erneut aus, um zu sehen, ob er funktioniert.

```bash title="nf-test pipeline pass"
nf-test test tests/main.nf.test
```

```console title="Pipeline passes with file assertions"
> nf-test test tests/main.nf.test

🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Workflow main.nf

  Test [f0e08a68] 'Should run successfully with correct number of processes' PASSED (8.144s)
  Test [d7e32a32] 'Should produce correct output files' PASSED (6.994s)


SUCCESS: Executed 2 tests in 15.165s
```

Erfolg! Die Tests sind bestanden, weil die Pipeline erfolgreich abgeschlossen wurde, die richtige Anzahl von Prozessen ausgeführt wurde und die Ausgabedateien erstellt wurden. Dies sollte dir auch zeigen, wie nützlich es ist, diese informativen Namen für deine Tests bereitzustellen.

Das ist nur die Oberfläche, wir können weiterhin Assertions schreiben, um die Details der Pipeline zu überprüfen, aber lass uns vorerst zum Testen der Interna der Pipeline übergehen.

### Fazit

Du weißt, wie man einen nf-test für eine Pipeline schreibt.

### Wie geht es weiter?

Lerne, wie man einen Nextflow-Prozess testet.

---

## 2. Teste einen Nextflow-Prozess

Wir müssen nicht für jeden Teil der Pipeline Tests schreiben, aber je mehr Tests wir haben, desto umfassender können wir über die Pipeline sein und desto sicherer können wir sein, dass sie wie erwartet funktioniert. In diesem Abschnitt werden wir beide Prozesse in der Pipeline als einzelne Einheiten testen.

### 2.1. Teste den `sayHello`-Prozess

Beginnen wir mit dem `sayHello`-Prozess.

Lass uns den `nf-test generate`-Befehl erneut verwenden, um Tests für den Prozess zu generieren.

```bash
nf-test generate process main.nf
```

```console title="Output"
> nf-test generate process main.nf

Load source file '/workspaces/training/side-quests/nf-test/main.nf'
Wrote process test file '/workspaces/training/side-quests/nf-test/tests/main.sayhello.nf.test
Wrote process test file '/workspaces/training/side-quests/nf-test/tests/main.converttoupper.nf.test

SUCCESS: Generated 2 test files.
```

Konzentrieren wir uns vorerst auf den `sayhello`-Prozess in der `main.sayhello.nf.test`-Datei.

Lass uns die Datei öffnen und den Inhalt ansehen.

```groovy title="tests/main.sayhello.nf.test"
nextflow_process {

    name "Test Process sayHello"
    script "main.nf"
    process "sayHello"

    test("Should run without failures") {

        when {
            params {
                // define parameters here. Example:
                // outdir = "tests/results"
            }
            process {
                """
                // define inputs of the process here. Example:
                // input[0] = file("test-file.txt")
                """
            }
        }

        then {
            assert process.success
            assert snapshot(process.out).match()
        }

    }

}
```

Wie zuvor beginnen wir mit den Testdetails, gefolgt von den `when`- und `then`-Blöcken. Wir haben jedoch auch einen zusätzlichen `process`-Block, der es uns ermöglicht, die Eingaben für den Prozess zu definieren.

Lass uns den Test ausführen, um zu sehen, ob er funktioniert.

```bash title="nf-test pipeline pass"
nf-test test tests/main.sayhello.nf.test
```

```console title="Process test fails"
> nf-test test tests/main.sayhello.nf.test

🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Process sayHello

  Test [1eaad118] 'Should run without failures' FAILED (4.876s)

  Assertion failed:

  assert process.success
         |       |
         |       false
         sayHello

  Nextflow stdout:

  Process `sayHello` declares 1 input but was called with 0 arguments
  Nextflow stderr:

FAILURE: Executed 1 tests in 4.884s (1 failed)
```

Der Test schlägt fehl, weil der `sayHello`-Prozess 1 Eingabe deklariert, aber mit 0 Argumenten aufgerufen wurde. Lass uns das beheben, indem wir dem Prozess eine Eingabe hinzufügen. Denk daran aus [Hello Workflow](../hello_nextflow/03_hello_workflow.md) (und dem Warmup-Abschnitt oben), dass unser `sayHello`-Prozess eine einzelne Werteingabe nimmt, die wir bereitstellen müssen. Wir sollten auch den Testnamen korrigieren, um besser widerzuspiegeln, was wir testen.

**Vorher:**

```groovy title="tests/main.sayhello.nf.test" linenums="1" hl_lines="1 10 11"
    test("Should run without failures") {

        when {
            params {
                // define parameters here. Example:
                // outdir = "tests/results"
            }
            process {
                """
                // define inputs of the process here. Example:
                // input[0] = file("test-file.txt")
                """
            }
        }

        then {
            assert process.success
            assert snapshot(process.out).match()
        }

    }
```

**Nachher:**

```groovy title="tests/main.sayhello.nf.test" linenums="1" hl_lines="1 10"
    test("Should run without failures and produce correct output") {

        when {
            params {
                // define parameters here. Example:
                // outdir = "tests/results"
            }
            process {
                """
                input[0] = "hello"
                """
            }
        }

        then {
            assert process.success
            assert snapshot(process.out).match()
        }

    }
```

Lass uns den Test erneut ausführen, um zu sehen, ob er funktioniert.

```console title="nf-test pipeline pass"
> nf-test test tests/main.sayhello.nf.test

🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Process sayHello

  Test [f91a1bcd] 'Should run without failures and produce correct output' PASSED (1.604s)
  Snapshots:
    1 created [Should run without failures and produce correct output]


Snapshot Summary:
  1 created

SUCCESS: Executed 1 tests in 1.611s
```

Erfolg! Der Test ist bestanden, weil der `sayHello`-Prozess erfolgreich ausgeführt wurde und die Ausgabe erstellt wurde.

### 2.2. Sieh dir den vom Test erstellten Snapshot an

Wenn wir uns die `tests/main.sayhello.nf.test`-Datei ansehen, können wir sehen, dass sie eine Methode `snapshot()` im Assertion-Block verwendet:

```groovy title="tests/main.sayhello.nf.test"
assert snapshot(process.out).match()
```

Dies teilt nf-test mit, einen Snapshot der Ausgabe des `sayHello`-Prozesses zu erstellen. Lass uns den Inhalt der Snapshot-Datei ansehen.

```console title="Snapshot file contents"
code tests/main.sayhello.nf.test.snap
```

Wir werden es hier nicht drucken, aber du solltest eine JSON-Datei sehen, die Details des Prozesses und der Prozessausgaben enthält. Insbesondere können wir eine Zeile sehen, die so aussieht:

```json title="Snapshot file contents"
"0": [
    "hello-output.txt:md5,b1946ac92492d2347c6235b4d2611184"
]
```

Dies repräsentiert die vom `sayHello`-Prozess erstellten Ausgaben, die wir explizit testen. Wenn wir den Test erneut ausführen, wird das Programm überprüfen, dass die neue Ausgabe mit der ursprünglich aufgezeichneten Ausgabe übereinstimmt. Dies ist eine schnelle, einfache Möglichkeit zu testen, dass sich Prozessausgaben nicht ändern, weshalb nf-test dies als Standard bereitstellt.

!!!warning

    Das bedeutet, wir müssen sicher sein, dass die Ausgabe, die wir im ursprünglichen Lauf aufzeichnen, korrekt ist!

Wenn sich im Laufe der zukünftigen Entwicklung etwas im Code ändert, das dazu führt, dass die Ausgabe anders ist, wird der Test fehlschlagen und wir müssen bestimmen, ob die Änderung erwartet ist oder nicht.

- Wenn sich herausstellt, dass etwas im Code kaputt gegangen ist, müssen wir es reparieren, mit der Erwartung, dass der reparierte Code den Test besteht
- Wenn es eine erwartete Änderung ist (z.B. das Tool wurde verbessert und die Ergebnisse sind besser), dann müssen wir den Snapshot aktualisieren, um die neue Ausgabe als Referenz zum Abgleich zu akzeptieren. nf-test hat einen Parameter `--update-snapshot` für diesen Zweck

Wir können den Test erneut ausführen und sehen, dass der Test bestehen sollte:

```console title="nf-test process pass with snapshot"
> nf-test test tests/main.sayhello.nf.test

🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Process sayHello

  Test [f91a1bcd] 'Should run without failures and produce correct output' PASSED (1.675s)


SUCCESS: Executed 1 tests in 1.685s
```

Erfolg! Der Test ist bestanden, weil der `sayHello`-Prozess erfolgreich ausgeführt wurde und die Ausgabe mit dem Snapshot übereinstimmte.

### 2.3. Alternative zu Snapshots: Direkte Content-Assertions

Während Snapshots großartig sind, um Änderungen in der Ausgabe zu erfassen, möchtest du manchmal spezifischen Inhalt überprüfen, ohne so streng zu sein, dass die gesamte Datei übereinstimmen muss. Zum Beispiel:

- Wenn Teile der Ausgabe sich ändern könnten (Zeitstempel, zufällige IDs usw.), aber bestimmte Schlüsselinhalte vorhanden sein müssen
- Wenn du nach spezifischen Mustern oder Werten in der Ausgabe suchen möchtest
- Wenn du den Test expliziter darüber machen möchtest, was Erfolg ausmacht

So könnten wir unseren Test modifizieren, um spezifischen Inhalt zu überprüfen:

**Vorher:**

```groovy title="tests/main.sayhello.nf.test" linenums="1" hl_lines="1 5 6 17"
    test("Should run without failures and produce correct output") {

        when {
            params {
                // define parameters here. Example:
                // outdir = "tests/results"
            }
            process {
                """
                input[0] = "hello"
                """
            }
        }

        then {
            assert process.success
            assert snapshot(process.out).match()
        }

    }
```

**Nachher:**

```groovy title="tests/main.sayhello.nf.test" linenums="1" hl_lines="1 5 16 17"
    test("Should run without failures and contain expected greeting") {

        when {
            params {
                // define parameters here
            }
            process {
                """
                input[0] = "hello"
                """
            }
        }

        then {
            assert process.success
            assert path(process.out[0][0]).readLines().contains('hello')
            assert !path(process.out[0][0]).readLines().contains('HELLO')
        }

    }
```

Beachte, dass nf-test die Prozessausgaben als Liste von Listen sieht, also holt `process.out[0][0]` den ersten Teil des ersten Kanalelements (oder 'Emission') von diesem Prozess.

Dieser Ansatz:

- Macht klar, was genau wir in der Ausgabe erwarten
- Ist widerstandsfähiger gegen irrelevante Änderungen in der Ausgabe
- Bietet bessere Fehlermeldungen, wenn Tests fehlschlagen
- Ermöglicht komplexere Validierungen (Regex-Muster, numerische Vergleiche usw.)

Lass uns den Test ausführen, um zu sehen, ob er funktioniert.

```bash title="nf-test pipeline pass"
nf-test test tests/main.sayhello.nf.test
```

```console title="Process test fails"
> nf-test test tests/main.sayhello.nf.test

🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Process sayHello

  Test [58df4e4b] 'Should run without failures and contain expected greeting' PASSED (7.196s)


SUCCESS: Executed 1 tests in 7.208s
```

### 2.4. Teste den `convertToUpper`-Prozess

Lass uns die `tests/main.converttoupper.nf.test`-Datei öffnen und den Inhalt ansehen:

```groovy title="tests/main.converttoupper.nf.test"
nextflow_process {

    name "Test Process convertToUpper"
    script "main.nf"
    process "convertToUpper"

    test("Should run without failures") {

        when {
            params {
                // define parameters here. Example:
                // outdir = "tests/results"
            }
            process {
                """
                // define inputs of the process here. Example:
                // input[0] = file("test-file.txt")
                """
            }
        }

        then {
            assert process.success
            assert snapshot(process.out).match()
        }

    }

}
```

Dies ist ein ähnlicher Test wie der `sayHello`-Prozess, aber er testet den `convertToUpper`-Prozess. Wir wissen, dass dieser fehlschlagen wird, weil genau wie bei `sayHello` der `convertToUpper`-Prozess eine einzelne Pfadeingabe nimmt, wir aber keine angegeben haben.

Wir müssen jetzt eine einzelne Eingabedatei für den convertToUpper-Prozess bereitstellen, die etwas Text enthält, den wir in Großbuchstaben konvertieren möchten. Es gibt viele Möglichkeiten, wie wir das tun könnten:

- Wir könnten eine dedizierte Datei zum Testen erstellen
- Wir könnten die vorhandene data/greetings.csv-Datei wiederverwenden
- Wir könnten sie spontan innerhalb des Tests erstellen

Lass uns vorerst die vorhandene data/greetings.csv-Datei wiederverwenden, indem wir das Beispiel verwenden, das wir beim Pipeline-Test verwendet haben. Wie zuvor können wir den Test benennen, um besser widerzuspiegeln, was wir testen, aber dieses Mal lassen wir ihn den Inhalt 'snapshotten', anstatt nach spezifischen Strings zu suchen (wie wir es im anderen Prozess getan haben).

**Vorher:**

```groovy title="tests/main.converttoupper.nf.test" linenums="1" hl_lines="1 10 11"
    test("Should run without failures") {

        when {
            params {
                // define parameters here. Example:
                // outdir = "tests/results"
            }
            process {
                """
                // define inputs of the process here. Example:
                // input[0] = file("test-file.txt")
                """
            }
        }

        then {
            assert process.success
            assert snapshot(process.out).match()
        }

    }
```

**Nachher:**

```groovy title="tests/main.converttoupper.nf.test" linenums="1" hl_lines="1 10"
    test("Should run without failures and produce correct output") {

        when {
            params {
                // define parameters here. Example:
                // outdir = "tests/results"
            }
            process {
                """
                input[0] = "${projectDir}/greetings.csv"
                """
            }
        }

        then {
            assert process.success
            assert snapshot(process.out).match()
        }

    }
```

Und führe den Test aus!

```bash title="nf-test pipeline pass"
nf-test test tests/main.converttoupper.nf.test
```

```console title="nf-test process convertToUpper pass"
> nf-test test tests/main.converttoupper.nf.test

🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Process convertToUpper

  Test [c59b6044] 'Should run without failures and produce correct output' PASSED (1.755s)
  Snapshots:
    1 created [Should run without failures and produce correct output]


Snapshot Summary:
  1 created

SUCCESS: Executed 1 tests in 1.764s
```

Beachte, wir haben eine Snapshot-Datei für den `convertToUpper`-Prozess unter `tests/main.converttoupper.nf.test.snap` erstellt. Wenn wir den Test erneut ausführen, sollten wir sehen, dass nf-test erneut besteht.

```bash title="nf-test process convertToUpper pass"
nf-test test tests/main.converttoupper.nf.test
```

```console title="nf-test process convertToUpper pass"
> nf-test test tests/main.converttoupper.nf.test

🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Process convertToUpper

  Test [c59b6044] 'Should run without failures and produce correct output' PASSED (1.798s)


SUCCESS: Executed 1 tests in 1.811s
```

### Fazit

Du weißt, wie man Tests für einen Nextflow-Prozess schreibt und ausführt.

### Wie geht es weiter?

Lerne, wie man Tests für alles auf einmal ausführt!

## 3. Führe Tests für das gesamte Repository aus

nf-test auf jeder Komponente auszuführen ist in Ordnung, aber mühsam und fehleranfällig. Können wir nicht einfach alles auf einmal testen?

Ja, das können wir!

Lass uns nf-test auf dem gesamten Repo ausführen.

### 3.1. Führe nf-test auf dem gesamten Repo aus

Wir können nf-test auf dem gesamten Repo ausführen, indem wir den Befehl `nf-test test` ausführen.

```bash
nf-test test .
```

Beachte, wir verwenden nur das `.`, um alles von unserem aktuellen Verzeichnis aus auszuführen. Dies wird jeden Test einschließen!

```console title="nf-test repo pass"
> nf-test test .

🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Process convertToUpper

  Test [3d26d9af] 'Should run without failures and produce correct output' PASSED (4.155s)

Test Workflow main.nf

  Test [f183df37] 'Should run successfully with correct number of processes' PASSED (3.33s)
  Test [d7e32a32] 'Should produce correct output files' PASSED (3.102s)

Test Process sayHello

  Test [58df4e4b] 'Should run without failures and contain expected greeting' PASSED (2.614s)


SUCCESS: Executed 4 tests in 13.481s
```

Sieh dir das an! Wir haben 4 Tests ausgeführt, 1 für jeden Prozess und 2 für die gesamte Pipeline mit einem einzigen Befehl. Stell dir vor, wie mächtig das bei einer großen Codebasis ist!

---

## Zusammenfassung

In dieser Side Quest hast du gelernt, die Funktionen von nf-test zu nutzen, um Tests für einzelne Prozesse sowie End-to-End-Tests für die gesamte Pipeline zu erstellen und auszuführen.
Du bist dir jetzt der beiden Hauptansätze zur Ausgabevalidierung bewusst, Snapshots und direkte Content-Assertions, und wann du welchen verwenden solltest.
Du weißt auch, wie man Tests entweder einzeln oder für ein gesamtes Projekt ausführt.

Die Anwendung dieser Techniken in deiner eigenen Arbeit wird es dir ermöglichen sicherzustellen, dass:

- Dein Code wie erwartet funktioniert
- Änderungen bestehende Funktionalität nicht kaputt machen
- Andere Entwickler\*innen mit Vertrauen beitragen können
- Probleme schnell identifiziert und behoben werden können
- Ausgabeinhalte den Erwartungen entsprechen

### Wichtige Muster

1. Pipeline-Tests:
   - Grundlegende Erfolgstests
   - Prozesszahl-Überprüfung
   - Ausgabedatei-Existenzprüfungen
2. Prozess-Tests
3. Zwei Ansätze zur Ausgabevalidierung:
   - Verwendung von Snapshots für vollständige Ausgabeüberprüfung
   - Verwendung direkter Content-Assertions für spezifische Inhaltsprüfungen
4. Ausführen aller Tests in einem Repository mit einem einzigen Befehl

### Zusätzliche Ressourcen

Sieh dir die [nf-test-Dokumentation](https://www.nf-test.com/) für erweiterte Testfunktionen und Best Practices an. Du möchtest vielleicht:

- Umfassendere Assertions zu deinen Tests hinzufügen
- Tests für Grenzfälle und Fehlerbedingungen schreiben
- Continuous Integration einrichten, um Tests automatisch auszuführen
- Mehr über andere Arten von Tests wie Workflow- und Modul-Tests erfahren
- Erweiterte Content-Validierungstechniken erkunden

**Denk daran:** Tests sind lebende Dokumentation darüber, wie sich dein Code verhalten sollte. Je mehr Tests du schreibst und je spezifischer deine Assertions sind, desto sicherer kannst du in die Zuverlässigkeit deiner Pipeline sein.

---

## Wie geht es weiter?

Kehre zum [Menü der Side Quests](./index.md) zurück oder klicke auf den Button unten rechts auf der Seite, um zum nächsten Thema in der Liste zu gelangen.
