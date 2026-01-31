# Teil 4: Tests hinzufügen

Im ersten Teil dieses Kurses hast du eine Varianten-Calling-Pipeline erstellt, die vollständig linear war und die Daten jeder Probe unabhängig von den anderen verarbeitet hat.

Im zweiten Teil haben wir dir gezeigt, wie du Channels und Channel-Operatoren verwendest, um gemeinsames Varianten-Calling mit GATK zu implementieren.

Im dritten Teil haben wir die Pipeline modularisiert.

In diesem Teil des Trainings zeigen wir dir, wie du [**nf-test**](https://www.nf-test.com/) verwendest, ein Test-Framework, das sich gut in Nextflow integriert und es einfach macht, sowohl Tests auf Modul-Ebene als auch auf Workflow-Ebene zu deiner Pipeline hinzuzufügen. Um diesem Teil des Trainings zu folgen, solltest du Teil 1, Teil 2 und Teil 3 abgeschlossen haben, sowie die [nf-test Side Quest](../../side_quests/nf-test.md), die die Grundlagen von nf-test behandelt und erklärt, warum Testen wichtig ist.

---

## 0. Aufwärmen

!!! note "Hinweis"

    Stelle sicher, dass du im richtigen Arbeitsverzeichnis bist:
    `cd /workspaces/training/nf4-science/genomics`

Wenn du die vorherigen Teile dieses Trainingskurses durchgearbeitet hast, solltest du eine funktionierende Version der Genomik-Pipeline mit der entsprechenden Modulverzeichnisstruktur haben.

??? abstract "Verzeichnisinhalte"

    ```console
    modules/
    ├── gatk
    │   ├── haplotypecaller
    │   │   └── main.nf
    │   └── jointgenotyping
    │       └── main.nf
    └── samtools
        └── index
            └── main.nf
    ```

Dieses Modulverzeichnis findest du im `solutions`-Verzeichnis, falls du es benötigst.

Wir beginnen mit demselben Workflow wie in Teil 3, den wir dir in der Datei `genomics-4.nf` bereitgestellt haben. Genau wie bei der [nf-test Side Quest](../../side_quests/nf-test.md) werden wir den drei Prozessen in dieser Pipeline verschiedene Arten von Tests hinzufügen sowie einen Test auf Workflow-Ebene.

### 0.1. Prüfen, ob der Workflow läuft

Bevor wir mit dem Hinzufügen von Tests beginnen, stelle sicher, dass der Workflow wie erwartet läuft.

```bash
nextflow run genomics-4.nf -resume
```

Das sollte inzwischen sehr vertraut aussehen, wenn du diesen Trainingskurs von Anfang an durchgearbeitet hast.

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `genomics-4.nf` [gloomy_poincare] DSL2 - revision: 43203316e0

    executor >  local (7)
    [18/89dfa4] SAMTOOLS_INDEX (1)       | 3 of 3 ✔
    [30/b2522b] GATK_HAPLOTYPECALLER (2) | 3 of 3 ✔
    [a8/d2c189] GATK_JOINTGENOTYPING     | 1 of 1 ✔
    ```

Wie zuvor wird es jetzt ein `work`-Verzeichnis und ein `results_genomics`-Verzeichnis in deinem Projektverzeichnis geben. Diese Ergebnisse werden wir später beim Testen tatsächlich verwenden. Aber ab jetzt werden wir das `nf-test`-Paket verwenden, um die Pipeline zu testen.

### 0.2. `nf-test` initialisieren

Wie bei der [nf-test Side Quest](../../side_quests/nf-test.md) müssen wir das `nf-test`-Paket initialisieren.

```bash
nf-test init
```

??? success "Befehlsausgabe"

    ```console
    🚀 nf-test 0.9.3
    https://www.nf-test.com
    (c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr

    Project configured. Configuration is stored in nf-test.config
    ```

??? abstract "nf-test.config-Inhalte"

    ```groovy title="nf-test.config"
    config {

        testsDir "tests"
        workDir ".nf-test"
        configFile "tests/nextflow.config"
        profile ""

    }
    ```

Es erstellt auch ein `tests`-Verzeichnis mit einer Konfigurationsdatei-Vorlage.

### Zusammenfassung

Jetzt sind wir bereit, Tests für unsere Genomik-Pipeline zu schreiben.

### Wie geht es weiter?

Schreibe grundlegende Tests, die prüfen, ob die Prozessaufrufe erfolgreich waren und die richtigen Ausgaben produziert haben.

---

## 1. Einen Prozess auf Erfolg und übereinstimmende Ausgaben testen

Wir beginnen mit dem Testen des `SAMTOOLS_INDEX`-Prozesses, der Indexdateien für BAM-Dateien erstellt, um effizienten Direktzugriff zu ermöglichen. Dies ist ein guter erster Testfall, weil:

1. Er eine einzelne, klar definierte Eingabe hat (eine BAM-Datei)
2. Er eine vorhersagbare Ausgabe produziert (eine BAI-Indexdatei)
3. Die Ausgabe bei identischen Eingaben identisch sein sollte

### 1.1. Eine Test-Datei-Vorlage generieren

Generiere zunächst eine Test-Datei-Vorlage:

```bash
nf-test generate process modules/samtools/index/main.nf
```

??? success "Befehlsausgabe"

    ```console
    Load source file '/workspaces/training/nf4-science/genomics/modules/samtools/index/main.nf'
    Wrote process test file '/workspaces/training/nf4-science/genomics/tests/modules/samtools/index/main.nf.test'

    SUCCESS: Generated 1 test files.
    ```

Dies erstellt eine Datei im selben Verzeichnis wie `main.nf`.
Du kannst im Datei-Explorer zu dem Verzeichnis navigieren und die Datei öffnen, die folgenden Code enthalten sollte:

```groovy title="tests/modules/samtools/index/main.nf.test" linenums="1"
nextflow_process {

    name "Test Process SAMTOOLS_INDEX"
    script "modules/samtools/index/main.nf"
    process "SAMTOOLS_INDEX"

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

Die anfänglichen Assertions sollten aus der [nf-test Side Quest](../../side_quests/nf-test.md) bekannt sein:

- `assert process.success` besagt, dass wir erwarten, dass der Prozess erfolgreich läuft und ohne Fehler abgeschlossen wird.
- `snapshot(process.out).match()` besagt, dass wir erwarten, dass das Ergebnis des Laufs mit dem Ergebnis eines vorherigen Laufs identisch ist (falls zutreffend).
  Wir besprechen dies später ausführlicher.

Mit diesem Ausgangspunkt müssen wir die richtigen Test-Eingaben für den samtools index-Prozess hinzufügen und gegebenenfalls Parameter.

### 1.2. Die Test-Datei verschieben und den Script-Pfad aktualisieren

Bevor wir mit dem Ausfüllen des Tests beginnen, müssen wir die Datei an ihren endgültigen Speicherort verschieben. Ein Teil des Grundes, warum wir ein Verzeichnis für jedes Modul hinzugefügt haben, ist, dass wir jetzt Tests in einem `tests`-Verzeichnis zusammen mit der `main.nf`-Datei jedes Moduls ablegen können. Erstelle dieses Verzeichnis und verschiebe die Test-Datei dorthin.

```bash
mkdir -p modules/samtools/index/tests
mv tests/modules/samtools/index/main.nf.test modules/samtools/index/tests/
```

Jetzt können wir den `script`-Abschnitt der Test-Datei zu einem relativen Pfad vereinfachen:

=== "Nachher"

    ```groovy title="modules/samtools/index/tests/main.nf.test" linenums="3" hl_lines="2"
    name "Test Process SAMTOOLS_INDEX"
    script "../main.nf"
    process "SAMTOOLS_INDEX"
    ```

=== "Vorher"

    ```groovy title="modules/samtools/index/tests/main.nf.test" linenums="3" hl_lines="2"
    name "Test Process SAMTOOLS_INDEX"
    script "modules/samtools/index/main.nf"
    process "SAMTOOLS_INDEX"
    ```

Dies teilt dem Test mit, wo die `main.nf`-Datei des Moduls zu finden ist, ohne den vollständigen Pfad angeben zu müssen.

### 1.3. Test-Eingaben für SAMTOOLS_INDEX bereitstellen

Die Vorlagendatei enthält einen Platzhalter, den wir durch eine tatsächliche Test-Eingabe ersetzen müssen, die für die Eingabe von `samtools index` geeignet ist. Die passende Eingabe ist eine BAM-Datei, die wir im Verzeichnis `data/bam` verfügbar haben.

=== "Nachher"

    ```groovy title="modules/samtools/index/tests/main.nf.test" linenums="14"
    process {
        """
        input[0] = file("${projectDir}/data/bam/reads_son.bam")
        """
    }
    ```

=== "Vorher"

    ```groovy title="modules/samtools/index/tests/main.nf.test" linenums="14"
    process {
        """
        // define inputs of the process here. Example:
        // input[0] = file("test-file.txt")
        """
    }
    ```

### 1.4. Den Test basierend auf der Funktionalität benennen

Wie wir zuvor gelernt haben, ist es eine gute Praxis, den Test so umzubenennen, dass er im Kontext des Tests sinnvoll ist.

=== "Nachher"

    ```groovy title="modules/samtools/index/tests/main.nf.test" linenums="7"
    test("Should index reads_son.bam correctly") {
    ```

    Dies nimmt eine beliebige Zeichenkette, wir könnten also alles eingeben, was wir wollen.
    Hier wählen wir, uns auf den Dateinamen und sein Format zu beziehen.

=== "Vorher"

    ```groovy title="modules/samtools/index/tests/main.nf.test" linenums="7"
    test("Should run without failures") {
    ```

### 1.5. Den Test ausführen und die Ausgabe untersuchen

Führe den Test aus:

```bash
nf-test test modules/samtools/index/tests/main.nf.test
```

??? success "Befehlsausgabe"

    ```console
    🚀 nf-test 0.9.3
    https://www.nf-test.com
    (c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


    Test Process SAMTOOLS_INDEX

      Test [625e39ee] 'Should index reads_son.bam correctly' PASSED (7.717s)
      Snapshots:
        1 created [Should index reads_son.bam correctly]


    Snapshot Summary:
      1 created

    SUCCESS: Executed 1 tests in 7.727s
    ```

Wie wir zuvor gelernt haben, hat dies die grundlegende Assertion über den Erfolg des Prozesses verifiziert und eine Snapshot-Datei basierend auf der Ausgabe des Prozesses erstellt. Wir können den Inhalt der Snapshot-Datei in der Datei `tests/modules/samtools/index/tests/main.nf.test.snap` sehen:

```json title="modules/samtools/index/tests/main.nf.test.snap" linenums="1"
{
  "Should index reads_son.bam correctly": {
    "content": [
      {
        "0": [
          [
            "reads_son.bam:md5,af5956d9388ba017944bef276b71d809",
            "reads_son.bam.bai:md5,a2ca7b84998218ee77eff14af8eb8ca2"
          ]
        ]
      }
    ],
    "meta": {
      "nf-test": "0.9.3",
      "nextflow": "25.10.2"
    },
    "timestamp": "2026-01-27T15:09:48.394063389"
  }
}
```

Wir können den Test auch erneut ausführen und sehen, dass er besteht, weil die Ausgabe mit dem Snapshot identisch ist:

```bash
nf-test test modules/samtools/index/tests/main.nf.test
```

??? success "Befehlsausgabe"

    ```console
    🚀 nf-test 0.9.3
    https://www.nf-test.com
    (c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


    Test Process SAMTOOLS_INDEX

      Test [625e39ee] 'Should index reads_son.bam correctly' PASSED (7.938s)


    SUCCESS: Executed 1 tests in 7.987s
    ```

### 1.6. Weitere Tests zu `SAMTOOLS_INDEX` hinzufügen

Manchmal ist es nützlich, eine Reihe verschiedener Eingabedateien zu testen, um sicherzustellen, dass wir auf verschiedene potenzielle Probleme testen. Füge Tests für die BAM-Dateien von Mutter und Vater im Trio aus unseren Testdaten hinzu.

```groovy
    test("Should index reads_mother.bam correctly") {

        when {
            params {
                // define parameters here
            }
            process {
                """
                input[0] = file("${projectDir}/data/bam/reads_mother.bam")
                """
            }
        }

        then {
            assert process.success
            assert snapshot(process.out).match()
        }

    }

    test("Should index reads_father.bam correctly") {

        when {
            params {
                // define parameters here
            }
            process {
                """
                input[0] = file("${projectDir}/data/bam/reads_father.bam")
                """
            }
        }

        then {
            assert process.success
            assert snapshot(process.out).match()
        }

    }
```

Dann kannst du den Test erneut ausführen:

```bash
nf-test test modules/samtools/index/tests/main.nf.test
```

??? success "Befehlsausgabe"

    ```console
    🚀 nf-test 0.9.3
    https://www.nf-test.com
    (c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


    Test Process SAMTOOLS_INDEX

      Test [625e39ee] 'Should index reads_son.bam correctly' PASSED (7.185s)
      Test [a8b28f36] 'Should index reads_mother.bam correctly' PASSED (6.576s)
      Test [c15852a1] 'Should index reads_father.bam correctly' PASSED (6.31s)
      Snapshots:
        2 created [Should index reads_father.bam correctly, Should index reads_mother.bam correctly]


    Snapshot Summary:
      2 created

    SUCCESS: Executed 3 tests in 20.117s
    ```

Beachte die Warnung, die sich auf die Wirkung des Parameters `--update-snapshot` bezieht.

!!! note "Hinweis"

    Hier verwenden wir Testdaten, die wir zuvor verwendet haben, um die wissenschaftlichen Ausgaben der Pipeline zu demonstrieren.
    Wenn wir geplant hätten, diese Tests in einer Produktionsumgebung zu betreiben, hätten wir kleinere Eingaben zu Testzwecken generiert.

    Im Allgemeinen ist es wichtig, Unit-Tests so leicht wie möglich zu halten, indem man die kleinsten notwendigen und ausreichenden Datenmengen zur Bewertung der Prozessfunktionalität verwendet, da sich sonst die Gesamtlaufzeit erheblich summieren kann.
    Eine Test-Suite, die zu lange dauert, um regelmäßig ausgeführt zu werden, ist eine Test-Suite, die wahrscheinlich im Interesse der Schnelligkeit übersprungen wird.

### Zusammenfassung

Du hast deinen ersten Modultest für einen Genomik-Prozess geschrieben und verifiziert, dass `SAMTOOLS_INDEX` korrekt Indexdateien für verschiedene BAM-Dateien erstellt. Die Test-Suite stellt sicher, dass:

1. Der Prozess erfolgreich läuft
2. Indexdateien erstellt werden
3. Die Ausgaben über Läufe hinweg konsistent sind
4. Der Prozess für alle Proben-BAM-Dateien funktioniert

### Wie geht es weiter?

Lerne, wie du Tests für andere Prozesse in unserem Genomik-Workflow schreibst, indem du die setup-Methode verwendest, um verkettete Prozesse zu handhaben. Wir werden auch prüfen, ob Ausgaben, insbesondere unsere VCF-Dateien, erwartete Variantenaufrufe enthalten.

---

## 2. Tests zu einem verketteten Prozess hinzufügen und auf Inhalte testen

Um `GATK_HAPLOTYPECALLER` zu testen, müssen wir dem Prozess die Ausgabe von `SAMTOOLS_INDEX` als Eingabe bereitstellen. Wir könnten das tun, indem wir `SAMTOOLS_INDEX` ausführen, seine Ausgaben abrufen und sie mit den Testdaten für den Workflow speichern. Das ist tatsächlich der empfohlene Ansatz für eine ausgereifte Pipeline, aber nf-test bietet einen alternativen Ansatz mit der `setup`-Methode.

Mit der setup-Methode können wir den `SAMTOOLS_INDEX`-Prozess als Teil des Test-Setups auslösen und dann seine Ausgabe als Eingabe für `GATK_HAPLOTYPECALLER` verwenden. Das hat einen Preis: Wir müssen den `SAMTOOLS_INDEX`-Prozess jedes Mal ausführen, wenn wir den Test für `GATK_HAPLOTYPECALLER` durchführen. Vielleicht entwickeln wir jedoch noch den Workflow und wollen keine Testdaten vorab generieren, die wir später möglicherweise ändern müssen. Der `SAMTOOLS_INDEX`-Prozess ist auch sehr schnell, sodass die Vorteile des Vorab-Generierens und Speicherns seiner Ausgaben vielleicht vernachlässigbar sind. So funktioniert die setup-Methode.

### 2.1. Die Test-Datei generieren und platzieren

Wie zuvor generieren wir zunächst die Datei-Vorlage:

```bash
nf-test generate process modules/gatk/haplotypecaller/main.nf
```

??? success "Befehlsausgabe"

    ```console
    Load source file '/workspaces/training/nf4-science/genomics/modules/gatk/haplotypecaller/main.nf'
    Wrote process test file '/workspaces/training/nf4-science/genomics/tests/modules/gatk/haplotypecaller/main.nf.test'

    SUCCESS: Generated 1 test files.
    ```

Dies erzeugt die folgende Test-Vorlage:

```groovy title="tests/modules/gatk/haplotypecaller/main.nf.test" linenums="1"
nextflow_process {

    name "Test Process GATK_HAPLOTYPECALLER"
    script "modules/gatk/haplotypecaller/main.nf"
    process "GATK_HAPLOTYPECALLER"

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

### 2.2. Die Test-Datei verschieben und den Script-Pfad aktualisieren

Wir erstellen ein Verzeichnis für die Test-Datei zusammen mit der `main.nf`-Datei des Moduls:

```bash
mkdir -p modules/gatk/haplotypecaller/tests
```

Und verschieben die Test-Vorlagendatei dorthin:

```bash
mv tests/modules/gatk/haplotypecaller/main.nf.test modules/gatk/haplotypecaller/tests/
```

Vergiss nicht, den Script-Pfad zu aktualisieren:

=== "Nachher"

    ```groovy title="modules/gatk/haplotypecaller/tests/main.nf.test" linenums="3" hl_lines="2"
        name "Test Process GATK_HAPLOTYPECALLER"
        script "../main.nf"
        process "GATK_HAPLOTYPECALLER"
    ```

=== "Vorher"

    ```groovy title="modules/gatk/haplotypecaller/tests/main.nf.test" linenums="3" hl_lines="2"
        name "Test Process GATK_HAPLOTYPECALLER"
        script "modules/gatk/haplotypecaller/main.nf"
        process "GATK_HAPLOTYPECALLER"
    ```

### 2.3. Eingaben mit der setup-Methode bereitstellen

Wir fügen vor dem `when`-Block einen `setup`-Block ein, in dem wir einen Lauf des `SAMTOOLS_INDEX`-Prozesses auf einer unserer ursprünglichen Eingabedateien auslösen können. Denke auch daran, wie zuvor den Testnamen in etwas Aussagekräftiges zu ändern.

=== "Nachher"

    ```groovy title="modules/gatk/haplotypecaller/tests/main.nf.test" linenums="7" hl_lines="1-12"
        test("Should call son's haplotype correctly") {

            setup {
                run("SAMTOOLS_INDEX") {
                    script "../../../samtools/index/main.nf"
                    process {
                        """
                        input[0] =  file("${projectDir}/data/bam/reads_son.bam")
                        """
                    }
                }
            }

            when {
    ```

=== "Vorher"

    ```groovy title="modules/gatk/haplotypecaller/tests/main.nf.test" linenums="7"  hl_lines="1"
    test("Should run without failures") {

        when {
    ```

Dann können wir auf die Ausgabe dieses Prozesses im `when`-Block verweisen, wo wir die Test-Eingaben spezifizieren:

```groovy title="modules/gatk/haplotypecaller/tests/main.nf.test" linenums="20"
        when {
            params {
                // define parameters here
            }
            process {
                """
                input[0] = SAMTOOLS_INDEX.out
                input[1] = file("${projectDir}/data/ref/ref.fasta")
                input[2] = file("${projectDir}/data/ref/ref.fasta.fai")
                input[3] = file("${projectDir}/data/ref/ref.dict")
                input[4] = file("${projectDir}/data/ref/intervals.bed")
                """
            }
        }
```

Führe diese Änderung durch und führe den Test erneut aus:

```bash
nf-test test modules/gatk/haplotypecaller/tests/main.nf.test
```

??? success "Befehlsausgabe"

    ```console
    🚀 nf-test 0.9.3
    https://www.nf-test.com
    (c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


    Test Process GATK_HAPLOTYPECALLER

      Test [c5156c2b] 'Should call son's haplotype correctly' PASSED (40.53s)
      Snapshots:
        1 created [Should call son's haplotype correctly]


    Snapshot Summary:
      1 created

    SUCCESS: Executed 1 tests in 40.555s
    ```

Es erzeugt auch wie zuvor eine Snapshot-Datei.

### 2.4. Erneut ausführen und Fehler beobachten

Interessanterweise wird der Test fehlschlagen, wenn du genau denselben Befehl noch einmal ausführst.

```bash
nf-test test modules/gatk/haplotypecaller/tests/main.nf.test
```

??? failure "Befehlsausgabe"

    ```console
    🚀 nf-test 0.9.3
    https://www.nf-test.com
    (c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


    Test Process GATK_HAPLOTYPECALLER

      Test [c5156c2b] 'Should call son's haplotype correctly' FAILED (40.123s)

      java.lang.RuntimeException: Different Snapshot:
      [                                                                                           [
          {                                                                                           {
              "0": [                                                                                      "0": [
                  "reads_son.bam.g.vcf:md5,069316cdd4328542ffc6ae247b1dac39"                         |                 "reads_son.bam.g.vcf:md5,005f1a13ee39f11b0fc9bea094850eac"
              ],                                                                                          ],
              "1": [                                                                                      "1": [
                  "reads_son.bam.g.vcf.idx:md5,dc36c18f2afdc546f41e68b2687e9334"                     |                 "reads_son.bam.g.vcf.idx:md5,dbad4b76a4b90c158ffc9c9740764242"
              ],                                                                                          ],
              "idx": [                                                                                    "idx": [
                  "reads_son.bam.g.vcf.idx:md5,dc36c18f2afdc546f41e68b2687e9334"                     |                 "reads_son.bam.g.vcf.idx:md5,dbad4b76a4b90c158ffc9c9740764242"
              ],                                                                                          ],
              "vcf": [                                                                                    "vcf": [
                  "reads_son.bam.g.vcf:md5,069316cdd4328542ffc6ae247b1dac39"                         |                 "reads_son.bam.g.vcf:md5,005f1a13ee39f11b0fc9bea094850eac"
              ]                                                                                           ]
          }                                                                                           }
      ]                                                                                           ]

      Nextflow stdout:

      Nextflow stderr:


        Obsolete snapshots can only be checked if all tests of a file are executed successful.


    FAILURE: Executed 1 tests in 40.156s (1 failed)
    ```

Die Fehlermeldung sagt dir, dass es Unterschiede zwischen den Snapshots für die beiden Läufe gab; insbesondere sind die md5sum-Werte für die VCF-Dateien unterschiedlich.

Warum? Um es kurz zu machen: Das HaplotypeCaller-Tool fügt einen Zeitstempel in den VCF-Header ein, der jedes Mal anders ist (per Definition).
Folglich können wir nicht einfach erwarten, dass die Dateien identische md5sums haben, selbst wenn sie in Bezug auf die Variantenaufrufe selbst identischen Inhalt haben.

Wie gehen wir damit um?

### 2.5. Eine Inhaltsassertion-Methode verwenden, um eine bestimmte Variante zu prüfen

Eine Möglichkeit, das Problem zu lösen, besteht darin, eine [andere Art von Assertion](https://nf-co.re/docs/contributing/tutorials/nf-test_assertions) zu verwenden.
In diesem Fall werden wir auf spezifischen Inhalt prüfen, anstatt Identität zu behaupten.
Genauer gesagt lassen wir das Tool die Zeilen der VCF-Datei lesen und prüfen auf das Vorhandensein bestimmter Zeilen.

In der Praxis ersetzen wir die zweite Assertion im `then`-Block wie folgt:

=== "Nachher"

    ```groovy title="modules/gatk/haplotypecaller/tests/main.nf.test" linenums="35" hl_lines="3 4"
            then {
                assert process.success
                assert path(process.out[0][0]).readLines().contains('#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	reads_son')
                assert path(process.out[0][0]).readLines().contains('20_10037292_10066351	3277	.	G	<NON_REF>	.	.	END=3282	GT:DP:GQ:MIN_DP:PL	0/0:25:72:24:0,72,719')
            }
    ```

=== "Vorher"

    ```groovy title="modules/gatk/haplotypecaller/tests/main.nf.test" linenums="35" hl_lines="3"
    then {
        assert process.success
        assert snapshot(process.out).match()
    }
    ```

Hier lesen wir den vollständigen Inhalt der VCF-Ausgabedatei ein und suchen nach einer Inhaltsübereinstimmung, was bei einer kleinen Testdatei in Ordnung ist, aber du würdest das nicht bei einer größeren Datei tun wollen.
Stattdessen könntest du dich entscheiden, bestimmte Zeilen einzulesen.

Dieser Ansatz erfordert eine sorgfältigere Auswahl dessen, was wir als "Signal" zum Testen verwenden möchten.
Auf der positiven Seite kann es verwendet werden, um mit großer Präzision zu testen, ob ein Analyse-Tool konsistent "schwierige" Merkmale (wie seltene Varianten) identifizieren kann, während es sich weiterentwickelt.

### 2.6. Erneut ausführen und Erfolg beobachten

Sobald wir den Test auf diese Weise geändert haben, können wir den Test mehrmals ausführen, und er wird konsistent bestehen.

```bash
nf-test test modules/gatk/haplotypecaller/tests/main.nf.test
```

??? success "Befehlsausgabe"

    ```console
    🚀 nf-test 0.9.3
    https://www.nf-test.com
    (c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


    Test Process GATK_HAPLOTYPECALLER

      Test [c5156c2b] 'Should call son's haplotype correctly' PASSED (40.53s)


    SUCCESS: Executed 1 tests in 40.555s
    ```

### 2.7. Weitere Tests hinzufügen

Füge ähnliche Tests für die Mutter- und Vaterproben hinzu:

```groovy title="modules/gatk/haplotypecaller/tests/main.nf.test" linenums="43"
    test("Should call mother's haplotype correctly") {

        setup {
            run("SAMTOOLS_INDEX") {
                script "../../../samtools/index/main.nf"
                process {
                    """
                    input[0] =  file("${projectDir}/data/bam/reads_mother.bam")
                    """
                }
            }
        }

        when {
            params {
                // define parameters here
            }
            process {
                """
                input[0] = SAMTOOLS_INDEX.out
                input[1] = file("${projectDir}/data/ref/ref.fasta")
                input[2] = file("${projectDir}/data/ref/ref.fasta.fai")
                input[3] = file("${projectDir}/data/ref/ref.dict")
                input[4] = file("${projectDir}/data/ref/intervals.bed")
                """
            }
        }

        then {
            assert process.success
            assert path(process.out[0][0]).readLines().contains('#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	reads_mother')
            assert path(process.out[0][0]).readLines().contains('20_10037292_10066351	3277	.	G	<NON_REF>	.	.	END=3278	GT:DP:GQ:MIN_DP:PL	0/0:38:99:37:0,102,1530')
        }
    }

    test("Should call father's haplotype correctly") {

        setup {
            run("SAMTOOLS_INDEX") {
                script "../../../samtools/index/main.nf"
                process {
                    """
                    input[0] =  file("${projectDir}/data/bam/reads_father.bam")
                    """
                }
            }
        }

        when {
            params {
                // define parameters here
            }
            process {
                """
                input[0] = SAMTOOLS_INDEX.out
                input[1] = file("${projectDir}/data/ref/ref.fasta")
                input[2] = file("${projectDir}/data/ref/ref.fasta.fai")
                input[3] = file("${projectDir}/data/ref/ref.dict")
                input[4] = file("${projectDir}/data/ref/intervals.bed")
                """
            }
        }

        then {
            assert process.success
            assert path(process.out[0][0]).readLines().contains('#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	reads_father')
            assert path(process.out[0][0]).readLines().contains('20_10037292_10066351	3277	.	G	<NON_REF>	.	.	END=3281	GT:DP:GQ:MIN_DP:PL	0/0:44:99:42:0,120,1800')
        }
    }
```

### 2.8. Den Testbefehl ausführen

```bash
nf-test test modules/gatk/haplotypecaller/tests/main.nf.test
```

??? success "Befehlsausgabe"

    ```console
    🚀 nf-test 0.9.3
    https://www.nf-test.com
    (c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


    Test Process GATK_HAPLOTYPECALLER

      Test [c5156c2b] 'Should call son's haplotype correctly' PASSED (40.53s)
      Test [10de94a8] 'Should call mother's haplotype correctly' PASSED (41.47s)
      Test [c0386fc7] 'Should call father's haplotype correctly' PASSED (45.556s)


    SUCCESS: Executed 3 tests in 127.586s
    ```

Das vervollständigt den grundlegenden Testplan für diesen zweiten Schritt in der Pipeline. Weiter zum dritten und letzten Test auf Modulebene!

### Zusammenfassung

Du hast gelernt, wie man:

1. Prozesse testet, die von Ausgaben anderer Prozesse abhängen
2. Spezifische genomische Varianten in VCF-Ausgabedateien verifiziert
3. Nicht-deterministische Ausgaben handhabt, indem man spezifischen Inhalt prüft
4. Variantenaufruf über mehrere Proben hinweg testet

### Wie geht es weiter?

Lerne, wie man Tests schreibt, die vorab generierte Testdaten für den Joint-Genotyping-Schritt verwenden.

---

## 3. Vorab generierte Testdaten verwenden

Für den Joint-Genotyping-Schritt verwenden wir einen anderen Ansatz - die Verwendung vorab generierter Testdaten. Dies ist oft vorzuziehen für:

1. Komplexe Prozesse mit mehreren Abhängigkeiten
2. Prozesse, die lange laufen
3. Prozesse, die Teil einer stabilen Produktions-Pipeline sind

### 3.1. Testdaten generieren

Untersuche die Ergebnisse, die wir zu Beginn dieses Abschnitts generiert haben:

```bash
tree results_genomics/
```

```console title="Inhalt des Ergebnisverzeichnisses"
results_genomics/
├── family_trio.joint.vcf
├── family_trio.joint.vcf.idx
├── gvcf
│   ├── reads_father.bam.g.vcf -> /workspaces/training/nf4-science/genomics/work/30/b2522b83c63baff8c3cf75704512a2/reads_father.bam.g.vcf
│   ├── reads_father.bam.g.vcf.idx -> /workspaces/training/nf4-science/genomics/work/30/b2522b83c63baff8c3cf75704512a2/reads_father.bam.g.vcf.idx
│   ├── reads_mother.bam.g.vcf -> /workspaces/training/nf4-science/genomics/work/f6/be2efa58e625d08cf8d0da1d0e9f09/reads_mother.bam.g.vcf
│   ├── reads_mother.bam.g.vcf.idx -> /workspaces/training/nf4-science/genomics/work/f6/be2efa58e625d08cf8d0da1d0e9f09/reads_mother.bam.g.vcf.idx
│   ├── reads_son.bam.g.vcf -> /workspaces/training/nf4-science/genomics/work/fe/2f22d56aa16ed45f8bc419312894f6/reads_son.bam.g.vcf
│   └── reads_son.bam.g.vcf.idx -> /workspaces/training/nf4-science/genomics/work/fe/2f22d56aa16ed45f8bc419312894f6/reads_son.bam.g.vcf.idx
└── indexed_bam
    ├── reads_father.bam -> /workspaces/training/nf4-science/genomics/work/42/a3bf19dbfaf1f3672b16a5d5e6a8be/reads_father.bam
    ├── reads_father.bam.bai -> /workspaces/training/nf4-science/genomics/work/cf/289c2d264f496d60a69e3e9ba6463e/reads_father.bam.bai
    ├── reads_mother.bam -> /workspaces/training/nf4-science/genomics/work/af/f31a6ade82cc0cf853c4f61c8bc473/reads_mother.bam
    ├── reads_mother.bam.bai -> /workspaces/training/nf4-science/genomics/work/18/89dfa40a3def17e45421e54431a126/reads_mother.bam.bai
    ├── reads_son.bam -> /workspaces/training/nf4-science/genomics/work/9f/9615dd553d6f13d8bec4f006ac395f/reads_son.bam
    └── reads_son.bam.bai -> /workspaces/training/nf4-science/genomics/work/4d/cb384a97db5687cc9daab002017c7c/reads_son.bam.bai

2 directories, 14 files
```

Der Joint-Genotyping-Schritt benötigt die von den Haplotype-Caller-Schritten erzeugten VCF-Dateien als Eingaben, zusammen mit den Indizes. Also lass uns die Ergebnisse, die wir haben, in das Test-Verzeichnis des `jointgenotyping`-Moduls kopieren.

```bash
mkdir -p modules/gatk/jointgenotyping/tests/inputs/
cp results_genomics/gvcf/*.g.vcf results_genomics/gvcf/*.g.vcf.idx modules/gatk/jointgenotyping/tests/inputs/
```

Jetzt können wir diese Dateien als Eingaben für den Test verwenden, den wir für den Joint-Genotyping-Schritt schreiben werden.

### 3.2. Die Test-Datei-Vorlage generieren

Wie zuvor generieren wir zunächst die Datei-Vorlage:

```bash
nf-test generate process modules/gatk/jointgenotyping/main.nf
```

??? success "Befehlsausgabe"

    ```console
    Load source file '/workspaces/training/nf4-science/genomics/modules/gatk/jointgenotyping/main.nf'
    Wrote process test file '/workspaces/training/nf4-science/genomics/tests/modules/gatk/jointgenotyping/main.nf.test'

    SUCCESS: Generated 1 test files.
    ```

Dies erzeugt die folgende Test-Vorlage:

```groovy title="tests/modules/gatk/jointgenotyping/main.nf.test" linenums="1"
nextflow_process {

    name "Test Process GATK_JOINTGENOTYPING"
    script "modules/gatk/jointgenotyping/main.nf"
    process "GATK_JOINTGENOTYPING"

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

### 3.3. Die Test-Datei verschieben und den Script-Pfad aktualisieren

Diesmal haben wir bereits ein Verzeichnis für Tests zusammen mit der `main.nf`-Datei des Moduls, sodass wir die Test-Vorlagendatei dorthin verschieben können:

```bash
mv tests/modules/gatk/jointgenotyping/main.nf.test modules/gatk/jointgenotyping/tests/
```

Und vergiss nicht, den Script-Pfad zu aktualisieren:

=== "Nachher"

    ```groovy title="modules/gatk/jointgenotyping/tests/main.nf.test" linenums="3" hl_lines="2"
    name "Test Process GATK_JOINTGENOTYPING"
    script "../main.nf"
    process "GATK_JOINTGENOTYPING"
    ```

=== "Vorher"

    ```groovy title="modules/gatk/jointgenotyping/tests/main.nf.test" linenums="3" hl_lines="2"
    name "Test Process GATK_JOINTGENOTYPING"
    script "modules/gatk/jointgenotyping/main.nf"
    process "GATK_JOINTGENOTYPING"
    ```

### 3.4. Eingaben bereitstellen

Fülle die Eingaben basierend auf den Prozess-Eingabedefinitionen aus und benenne den Test entsprechend um:

```groovy title="modules/gatk/jointgenotyping/tests/main.nf.test" linenums="7"
    test("Should call trio's joint genotype correctly") {

        when {
            params {
                // define parameters here
            }
            process {
                """
                input[0] = [
                    file("${projectDir}/modules/gatk/jointgenotyping/tests/inputs/reads_father.bam.g.vcf"),
                    file("${projectDir}/modules/gatk/jointgenotyping/tests/inputs/reads_mother.bam.g.vcf"),
                    file("${projectDir}/modules/gatk/jointgenotyping/tests/inputs/reads_son.bam.g.vcf")
                ]
                input[1] = [
                    file("${projectDir}/modules/gatk/jointgenotyping/tests/inputs/reads_father.bam.g.vcf.idx"),
                    file("${projectDir}/modules/gatk/jointgenotyping/tests/inputs/reads_mother.bam.g.vcf.idx"),
                    file("${projectDir}/modules/gatk/jointgenotyping/tests/inputs/reads_son.bam.g.vcf.idx")
                ]
                input[2] = file("${projectDir}/data/ref/intervals.bed")
                input[3] = "family_trio"
                input[4] = file("${projectDir}/data/ref/ref.fasta")
                input[5] = file("${projectDir}/data/ref/ref.fasta.fai")
                input[6] = file("${projectDir}/data/ref/ref.dict")
                """
            }
        }
```

### 3.5. Inhaltsassertions verwenden

Die Ausgabe des Joint-Genotyping-Schritts ist wieder eine VCF-Datei, also werden wir wieder eine Inhaltsassertion verwenden.

```groovy title="modules/gatk/jointgenotyping/tests/main.nf.test" linenums="25"
    then {
        assert process.success
        assert path(process.out[0][0]).readLines().contains('#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	reads_father	reads_mother	reads_son')
        assert path(process.out[0][0]).readLines().contains('20_10037292_10066351	3480	.	C	CT	1625.89	.	AC=5;AF=0.833;AN=6;BaseQRankSum=0.220;DP=85;ExcessHet=0.0000;FS=2.476;MLEAC=5;MLEAF=0.833;MQ=60.00;MQRankSum=0.00;QD=21.68;ReadPosRankSum=-1.147e+00;SOR=0.487	GT:AD:DP:GQ:PL	0/1:15,16:31:99:367,0,375	1/1:0,18:18:54:517,54,0	1/1:0,26:26:78:756,78,0')
    }
```

Durch die Überprüfung des Inhalts einer bestimmten Variante in der Ausgabedatei verifiziert dieser Test, dass:

1. Der Joint-Genotyping-Prozess erfolgreich läuft
2. Die Ausgabe-VCF alle drei Proben in der richtigen Reihenfolge enthält
3. Eine spezifische Variante korrekt aufgerufen wird mit:
   - Genauen Genotypen für jede Probe (0/1 für Vater, 1/1 für Mutter und Sohn)
   - Korrekten Read-Tiefen und Genotyp-Qualitäten
   - Populationsweiten Statistiken wie Allelfrequenz (AF=0.833)

Wir haben nicht die gesamte Datei als Snapshot gespeichert, aber indem wir eine spezifische Variante überprüfen, können wir sicher sein, dass der Joint-Genotyping-Prozess wie erwartet funktioniert.

### 3.6. Den Test ausführen

```bash
nf-test test modules/gatk/jointgenotyping/tests/main.nf.test
```

??? success "Befehlsausgabe"

    ```console
    🚀 nf-test 0.9.3
    https://www.nf-test.com
    (c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


    Test Process GATK_JOINTGENOTYPING

      Test [ac2067de] 'Should call trio's joint genotype correctly' PASSED (53.827s)


    SUCCESS: Executed 1 tests in 53.837s
    ```

Der Test besteht und verifiziert, dass unser Joint-Genotyping-Prozess korrekt:

1. Individuelle Proben-VCFs kombiniert
2. Gemeinsames Varianten-Calling durchführt
3. Eine Multi-Proben-VCF mit konsistenten Genotyp-Aufrufen über Läufe hinweg produziert

### Zusammenfassung

Du weißt, wie man:

- Zuvor generierte Ergebnisse als Eingaben für Tests verwendet
- Tests mit vorab generierten Testdaten schreibt

### Wie geht es weiter?

Füge einen Test auf Workflow-Ebene hinzu, um zu verifizieren, dass die gesamte Varianten-Calling-Pipeline end-to-end funktioniert.

---

## 4. Einen Test auf Workflow-Ebene hinzufügen

Jetzt testen wir die vollständige Varianten-Calling-Pipeline, von BAM-Dateien bis zu Joint-Genotypen. Dies verifiziert, dass:

1. Alle Prozesse korrekt zusammenarbeiten
2. Daten ordnungsgemäß zwischen den Schritten fließen
3. Abschließende Variantenaufrufe konsistent sind

### 4.1. Den Workflow-Test generieren

Generiere eine Test-Datei für die vollständige Pipeline:

```bash
nf-test generate pipeline genomics-4.nf
```

??? success "Befehlsausgabe"

    ```console
    Load source file '/workspaces/training/nf4-science/genomics/genomics-4.nf'
    Wrote pipeline test file '/workspaces/training/nf4-science/genomics/tests/genomics-4.nf.test'

    SUCCESS: Generated 1 test files.
    ```

Dies erstellt eine grundlegende Test-Vorlage:

```groovy title="tests/genomics-4.nf.test" linenums="1"
nextflow_pipeline {

    name "Test Workflow genomics-4.nf"
    script "genomics-4.nf"

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

Korrigiere einfach den Namen zu etwas Aussagekräftigem (du wirst bald sehen, warum das nützlich ist).

=== "Nachher"

    ```groovy title="tests/genomics-4.nf.test" linenums="6" hl_lines="1"
        test("Should run the pipeline without failures") {
    ```

=== "Vorher"

    ```groovy title="tests/genomics-4.nf.test" linenums="6" hl_lines="1"
        test("Should run without failures") {
    ```

!!! note "Hinweis"

    In diesem Fall kann die Test-Datei dort bleiben, wo `nf-test` sie erstellt hat.

### 4.2. Eingabeparameter spezifizieren

Wir müssen noch Eingaben spezifizieren, was auf Workflow-Ebene etwas anders gemacht wird als bei Tests auf Modulebene.
Es gibt mehrere Möglichkeiten, dies zu tun, einschließlich durch Angabe eines Profils.
Ein einfacherer Weg ist jedoch, einen `params {}`-Block in der `nextflow.config`-Datei einzurichten, die `nf-test init` ursprünglich im `tests`-Verzeichnis erstellt hat.

```groovy title="tests/nextflow.config" linenums="1"
/*
========================================================================================
    Nextflow-Konfigurationsdatei zum Ausführen von Tests
========================================================================================
*/

// Ausgabeverzeichnis für Workflow-Ausgaben
outputDir = 'results_genomics'

/*
 * Pipeline-Parameter
 */

params {
    // Primäre Eingabe (Datei mit Eingabedateien, eine pro Zeile)
    reads_bam = "${projectDir}/data/sample_bams.txt"

    // Zusatzdateien
    reference = "${projectDir}/data/ref/ref.fasta"
    reference_index = "${projectDir}/data/ref/ref.fasta.fai"
    reference_dict = "${projectDir}/data/ref/ref.dict"
    intervals = "${projectDir}/data/ref/intervals.bed"

    // Basisname für abschließende Ausgabedatei
    cohort_name = "family_trio"
}
```

Wenn wir den Test ausführen, übernimmt `nf-test` diese Konfigurationsdatei und zieht die Eingaben entsprechend heran.

### 4.3. Den Workflow-Test ausführen

```bash
nf-test test tests/genomics-4.nf.test
```

??? success "Befehlsausgabe"

    ```console
    🚀 nf-test 0.9.3
    https://www.nf-test.com
    (c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


    Test Workflow genomics-4.nf

      Test [1b4c6936] 'Should run the pipeline without failures' PASSED (171.019s)


    SUCCESS: Executed 1 tests in 171.056s
    ```

Der Test besteht und bestätigt, dass unsere vollständige Varianten-Calling-Pipeline:

1. Alle Proben erfolgreich verarbeitet
2. Alle Schritte korrekt verkettet

### 4.4. ALLE Tests ausführen

nf-test hat noch einen weiteren Trick auf Lager. Wir können alle Tests auf einmal ausführen! Ändere die `nf-test.config`-Datei so, dass nf-test in jedem Verzeichnis nach nf-test-Dateien sucht. Du kannst dies tun, indem du den `testsDir`-Parameter änderst:

=== "Nachher"

    ```groovy title="nf-test.config" linenums="1" hl_lines="3"
    config {

        testsDir "."
        workDir ".nf-test"
        configFile "tests/nextflow.config"
        profile ""

    }
    ```

=== "Vorher"

    ```groovy title="nf-test.config" linenums="1" hl_lines="3"
    config {

        testsDir "tests"
        workDir ".nf-test"
        configFile "tests/nextflow.config"
        profile ""

    }
    ```

Jetzt können wir einfach nf-test ausführen und es wird _jeden einzelnen Test_ in unserem Repository ausführen:

```bash
nf-test test
```

??? success "Befehlsausgabe"

    ```console
    🚀 nf-test 0.9.3
    https://www.nf-test.com
    (c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


    Test Process GATK_HAPLOTYPECALLER

      Test [c5156c2b] 'Should call son's haplotype correctly' PASSED (39.947s)
      Test [10de94a8] 'Should call mother's haplotype correctly' PASSED (43.17s)
      Test [c0386fc7] 'Should call father's haplotype correctly' PASSED (44.244s)

    Test Process GATK_JOINTGENOTYPING

      Test [ac2067de] 'Should call trio's joint genotype correctly' PASSED (61.129s)

    Test Process SAMTOOLS_INDEX

      Test [625e39ee] 'Should index reads_son.bam correctly' PASSED (8.671s)
      Test [a8b28f36] 'Should index reads_mother.bam correctly' PASSED (8.518s)
      Test [c15852a1] 'Should index reads_father.bam correctly' PASSED (5.378s)

    Test Workflow genomics-4.nf

      Test [1b4c6936] 'Should run the pipeline without failures' PASSED (169.714s)


    SUCCESS: Executed 8 tests in 380.801s
    ```

8 Tests mit 1 Befehl! Wir haben lange damit verbracht, viele Tests zu konfigurieren, aber als es darum ging, sie auszuführen, war es sehr schnell und einfach. Du kannst sehen, wie nützlich das bei der Wartung einer großen Pipeline ist, die Hunderte verschiedener Elemente enthalten könnte. Wir investieren Zeit in das Schreiben von Tests einmal, damit wir Zeit sparen können, wenn wir sie viele Male ausführen.

Außerdem können wir das automatisieren! Stell dir vor, Tests laufen jedes Mal, wenn du oder ein Kollege versucht, neuen Code hinzuzufügen. So stellen wir sicher, dass unsere Pipelines einen hohen Standard aufrechterhalten.

## Zusammenfassung

Du weißt jetzt, wie man mehrere Arten von Tests für deine Genomik-Pipeline mit nf-test schreibt und ausführt. Dieses Test-Framework hilft sicherzustellen, dass dein Varianten-Calling-Workflow konsistente, zuverlässige Ergebnisse über verschiedene Umgebungen hinweg produziert und wenn du Code-Änderungen vornimmst.

Du hast gelernt, kritische Komponenten zu testen wie:

- Den `SAMTOOLS_INDEX`-Prozess, der BAM-Dateien für Varianten-Calling vorbereitet
- Den `GATK_HAPLOTYPECALLER`-Prozess, der Varianten in einzelnen Proben identifiziert
- Den `GATK_JOINTGENOTYPING`-Prozess, der Variantenaufrufe über eine Kohorte hinweg kombiniert

Du hast auch verschiedene Teststrategien implementiert, die spezifisch für Genomik-Daten sind:

- Verifizieren, dass VCF-Dateien erwartete Variantenaufrufe enthalten, trotz nicht-deterministischer Elemente wie Zeitstempel
- Testen mit einem Familien-Trio-Datensatz, um die ordnungsgemäße Variantenidentifikation über verwandte Proben hinweg sicherzustellen
- Überprüfen auf spezifische genomische Koordinaten und Varianteninformationen in deinen Ausgabedateien

Diese Testfähigkeiten sind wesentlich für die Entwicklung robuster Bioinformatik-Pipelines, die zuverlässig Genomdaten verarbeiten und genaue Variantenaufrufe produzieren können. Während du weiterhin mit Nextflow für Genomik-Analysen arbeitest, wird dir diese Test-Grundlage helfen, qualitativ hochwertigen Code aufrechtzuerhalten, der vertrauenswürdige wissenschaftliche Ergebnisse produziert.
