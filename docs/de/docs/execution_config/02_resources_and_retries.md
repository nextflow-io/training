# Teil 2: Rechenressourcen und Fehler verwalten

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

In [Teil 1](./01_packaging_and_execution.md) hast du angepasst, wo und wie die Aufgaben einer Pipeline ausgeführt werden.
Hier passt du an, wie viele Rechenressourcen jede Aufgabe erhält und was passiert, wenn eine Aufgabe trotz deiner besten Schätzung fehlschlägt.

---

## 1. Rechenressourcen zuweisen

Standardmäßig weist Nextflow jedem Prozess über die `cpus`-Direktive eine einzelne CPU zu und setzt kein Speicherlimit, sofern du keines festlegst:

```groovy title="Built-in configuration"
process {
    cpus = 1
}
```

Aus [Nextflow Run](../nextflow_run/index.md) weißt du bereits, dass die Konfiguration dieser Pipeline `memory` für alle Prozesse auf 1 GB setzt.
Aber woher weißt du, welche Werte du für deine eigenen Pipelines verwenden solltest?

### 1.1. Einen Ressourcennutzungsbericht erstellen

In [Nextflow Run](../nextflow_run/02_configure_pipeline.md) hast du bereits einen Ausführungsbericht mit `-with-report` erstellt.
Genau dieser Bericht zeigt dir, wie viel CPU und Arbeitsspeicher deine Prozesse tatsächlich benötigen: Führe den Workflow mit einigen Standardzuweisungen aus, erfasse die tatsächliche Nutzung und passe sie dann entsprechend an.

```bash
nextflow run main.nf -with-report report-config-1.html
```

Der Bericht ist eine HTML-Datei, die du im Browser öffnen kannst.
Er schlüsselt Laufzeit und Ressourcennutzung pro Prozess auf, einschließlich des prozentualen Anteils der tatsächlich genutzten zugewiesenen Ressourcen.
So sieht es für `cowpy` mit den aktuellen Standardwerten aus (1 CPU, 1 GB Arbeitsspeicher):

| Metrik                      | Wert   |
| --------------------------- | ------ |
| CPU-Nutzung                 | 116%   |
| Maximaler Speicherverbrauch | 6,4 MB |
| Zugewiesener Speicher       | 1 GB   |

`cowpy` nutzt deutlich weniger als 1% seiner 1-GB-Zuweisung. Der `%cpu`-Wert über 100% bedeutet lediglich, dass es kurzzeitig mehr als eine CPU-Einheit an Rechenleistung innerhalb des Containers nutzt.

Die vollständige Liste der verfügbaren Funktionen findest du unter [Reports](https://nextflow.io/docs/latest/reports.html).

### 1.2. Ressourcenzuweisungen für einen bestimmten Prozess festlegen

Der obige Bericht zeigt, dass `cowpy` gut innerhalb seiner aktuellen Zuweisung liegt. Angenommen, du möchtest ihm trotzdem mehr Spielraum geben, zum Beispiel weil du in der Produktion größere Eingaben erwartest.
Du kannst die Standardwerte für einen einzelnen Prozess mit `withName` überschreiben.

=== "Danach"

    ```groovy title="nextflow.config" linenums="6" hl_lines="5-6"
    process {
        memory = 1.GB
        withName: 'cowpy' {
            conda = 'conda-forge::cowpy==1.1.5'
            memory = 2.GB
            cpus = 2
        }
    }
    ```

=== "Vorher"

    ```groovy title="nextflow.config" linenums="6"
    process {
        memory = 1.GB
        withName: 'cowpy' {
            conda = 'conda-forge::cowpy==1.1.5'
        }
    }
    ```

Mit dieser Einstellung fordert jeder Prozess 1 GB Arbeitsspeicher und eine einzelne CPU an – außer `cowpy`, das 2 GB und 2 CPUs anfordert (zusätzlich zur `conda`-Einstellung aus [Teil 1](./01_packaging_and_execution.md)).

!!! info "Info"

    Wenn dein Rechner wenige CPUs hat und du eine hohe Anzahl pro Prozess zuweist, können sich Aufgaben hintereinander einreihen, da Nextflow nicht mehr CPUs anfordert als verfügbar sind.

Führe es erneut mit einem anderen Berichtsdateinamen aus, damit du vorher und nachher vergleichen kannst.

```bash
nextflow run main.nf -with-report report-config-2.html
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [voluminous_venter] revision: c3c85dec78

    executor >  local (8)
    [a1/0e96d4] sayHello (1)       | 3 of 3 ✔
    [a3/7173a3] convertToUpper (2) | 3 of 3 ✔
    [4f/a8ae3d] collectGreetings   | 1 of 1 ✔
    [91/3724f8] cowpy              | 1 of 1 ✔

    Outputs:

      /workspaces/training/execution-config/results

      first_output:
        - full_pipeline/intermediates/Bonjour-output.txt
        - full_pipeline/intermediates/Hola-output.txt
        - full_pipeline/intermediates/Hello-output.txt

      uppercased:
        - full_pipeline/intermediates/UPPER-Bonjour-output.txt
        - full_pipeline/intermediates/UPPER-Hello-output.txt
        - full_pipeline/intermediates/UPPER-Hola-output.txt

      collected: full_pipeline/intermediates/COLLECTED-batch-output.txt

      batch_report: full_pipeline/batch-report.txt

      cowpy_art: full_pipeline/cowpy-COLLECTED-batch-output.txt
    ```

Vergleich der beiden Berichte für `cowpy`:

| Metrik                      | Vorher (1 CPU, 1 GB) | Danach (2 CPUs, 2 GB) |
| --------------------------- | -------------------- | --------------------- |
| Maximaler Speicherverbrauch | 6,4 MB               | 6,4 MB                |
| CPU-Nutzung                 | 116%                 | 118%                  |

Die doppelte Zuweisung hat die tatsächliche Nutzung überhaupt nicht verändert. Das zeigt, dass die ursprünglichen 1 GB / 1 CPU für diese einfache Arbeitslast bereits großzügig waren.
Bei einer echten Pipeline, die nicht-triviale Daten verarbeitet, würden sich die Zahlen zwischen den Prozessen deutlich unterscheiden – genau deshalb erstellst du ein Profil, bevor du entscheidest, was du zuweist, anstatt zu raten.

### 1.3. Ressourcenlimits festlegen

Je nach Recheninfrastruktur kann es harte Beschränkungen geben, was du anfordern kannst, zum Beispiel eine clusterweite Obergrenze.
Die `resourceLimits`-Direktive ermöglicht es dir, diese Limits festzulegen:

```groovy title="Syntax example"
process {
    resourceLimits = [
        memory: 750.GB,
        cpus: 200,
        time: 30.d
    ]
}
```

Nextflow übersetzt diese in das, was der Ziel-Executor erwartet.
Wenn ein Prozess mehr als das Limit anfordert, wird die Anfrage begrenzt statt abgelehnt.

!!! warning "Warnung"

    Das kannst du in der Trainingsumgebung nicht ausführen, da es HPC-Infrastruktur erfordert, um eine Wirkung zu haben.

??? info "Institutionelle Referenzkonfigurationen"

    Das nf-core-Projekt pflegt eine [Sammlung von Konfigurationsdateien](https://nf-co.re/configs/), die von Institutionen weltweit geteilt werden und eine breite Palette von HPC- und Cloud-Executors abdecken.
    Sie sind ein nützlicher Ausgangspunkt, unabhängig davon, ob deine eigene Institution darunter ist.

### Fazit

Du weißt, wie du einen Profiling-Bericht erstellst, um die Ressourcennutzung zu bewerten, Ressourcenzuweisungen für einen bestimmten Prozess überschreibst und Zuweisungen mit `resourceLimits` begrenzt.

### Wie geht es weiter?

Lerne, wie du eine Pipeline dazu bringst, sich automatisch zu erholen, wenn eine Aufgabe fehlschlägt – unabhängig davon, ob deine Ressourcenschätzung richtig war.

---

## 2. Aufgabenfehler mit Wiederholungen behandeln

Profiling zeigt dir, was ein Prozess meistens benötigt, aber echte Arbeitslasten variieren: Eine Zuweisung, die für die meisten Eingaben ausreicht, kann für eine ungewöhnlich große Eingabe zu knapp sein, und Schätzungen können schlicht falsch sein.
Anstatt einen einzelnen fehlgeschlagenen Task den gesamten Lauf zum Absturz zu bringen, kann Nextflow eine fehlgeschlagene Aufgabe automatisch wiederholen und ihr dabei optional bei jedem Versuch mehr Ressourcen geben.

### 2.1. Eine fehlgeschlagene Aufgabe automatisch wiederholen

Um das in Aktion zu sehen, setze die Speicherzuweisung von `cowpy` absichtlich unter das, was es tatsächlich benötigt: Aus [1.1](#11-generate-a-resource-utilization-report) weißt du, dass es bei etwa 6,4 MB seinen Höchstwert erreicht, also sollten 6 MB knapp zu wenig sein.

=== "Danach"

    ```groovy title="nextflow.config" linenums="6" hl_lines="5-7"
    process {
        memory = 1.GB
        withName: 'cowpy' {
            conda = 'conda-forge::cowpy==1.1.5'
            memory = 6.MB
            errorStrategy = 'retry'
            maxRetries = 2
        }
    }
    ```

=== "Vorher"

    ```groovy title="nextflow.config" linenums="6"
    process {
        memory = 1.GB
        withName: 'cowpy' {
            conda = 'conda-forge::cowpy==1.1.5'
            memory = 2.GB
            cpus = 2
        }
    }
    ```

`errorStrategy` teilt Nextflow mit, was zu tun ist, wenn eine Aufgabe fehlschlägt: `'retry'` sendet die Aufgabe erneut ein, anstatt die gesamte Pipeline zu stoppen.
`maxRetries` begrenzt, wie viele zusätzliche Versuche sie bekommt, bevor Nextflow aufgibt.

```bash
nextflow run main.nf
```

??? failure "Befehlsausgabe (gekürzt)"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [desperate_brazil] revision: c3c85dec78

    executor >  local (10)
    [67/fe1f49] sayHello (1)       | 3 of 3 ✔
    [8a/f13335] convertToUpper (1) | 3 of 3 ✔
    [39/1b24ed] collectGreetings   | 1 of 1 ✔
    [7a/d5eb6f] cowpy              | 0 of 1, retries: 2 ✘
    [9d/b79eb3] NOTE: Process `cowpy` terminated with an error exit status (137) -- Execution is retried (1)
    [6d/1d9d84] NOTE: Process `cowpy` terminated with an error exit status (137) -- Execution is retried (2)
    ERROR ~ Error executing process > 'cowpy'

    Caused by:
      Process `cowpy` terminated with an error exit status (137)

    Command executed:
      cat COLLECTED-batch-output.txt | cowpy -c "turkey" > cowpy-COLLECTED-batch-output.txt

    Command exit status:
      137

    Command output:
      (empty)

    Command error:
      /usr/local/bin/_activate_current_env.sh: line 35:    14 Killed                  micromamba activate "${ENV_NAME:-base}"

    Work dir:
      /workspaces/training/execution-config/work/7a/d5eb6feeac0eed18d95d3da7a7aeb4

    Tip: when you have fixed the problem you can continue the execution adding the option `-resume` to the run command line

    -- Check '.nextflow.log' file for details
    ```

Exit-Code 137 ist das Standardsignal für einen Out-of-Memory-Kill: Der Container hatte nicht genug Arbeitsspeicher, um `cowpy` überhaupt auszuführen.
Nextflow hat die Aufgabe zweimal wiederholt – insgesamt drei Versuche, entsprechend `maxRetries = 2`.
Da sich die Speicherzuweisung zwischen den Versuchen nie geändert hat, ist jeder Versuch an derselben Grenze gescheitert. Sobald die Wiederholungen erschöpft sind, meldet Nextflow den Fehler vollständig und stoppt die Pipeline mit einem Nicht-Null-Status.

Wiederholungen allein beheben nichts, wenn sich die zugrunde liegende Ursache zwischen den Versuchen nicht ändert.

### 2.2. Ressourcen bei jeder Wiederholung erhöhen

Innerhalb einer Prozess-Direktive enthält `task.attempt` die aktuelle Versuchsnummer, beginnend bei 1.
Du kannst sie in einem closure verwenden, um eine Ressourcenzuweisung bei jeder Wiederholung zu skalieren.

=== "Danach"

    ```groovy title="nextflow.config" linenums="6" hl_lines="5 7"
    process {
        memory = 1.GB
        withName: 'cowpy' {
            conda = 'conda-forge::cowpy==1.1.5'
            memory = { 6.MB * task.attempt }
            errorStrategy = 'retry'
            maxRetries = 3
        }
    }
    ```

=== "Vorher"

    ```groovy title="nextflow.config" linenums="6"
    process {
        memory = 1.GB
        withName: 'cowpy' {
            conda = 'conda-forge::cowpy==1.1.5'
            memory = 6.MB
            errorStrategy = 'retry'
            maxRetries = 2
        }
    }
    ```

Führe den Workflow erneut aus:

```bash
nextflow run main.nf
```

??? success "Befehlsausgabe (gekürzt)"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [grave_joliot] revision: c3c85dec78

    executor >  local (9)
    [0f/211b8a] sayHello (2)       | 3 of 3 ✔
    [26/301ea2] convertToUpper (3) | 3 of 3 ✔
    [22/5a895b] collectGreetings   | 1 of 1 ✔
    [e1/beee86] cowpy              | 1 of 1, retries: 1 ✔
    [b6/7aed6a] NOTE: Process `cowpy` terminated with an error exit status (137) -- Execution is retried (1)

    Outputs:

      /workspaces/training/execution-config/results

      first_output:
        - full_pipeline/intermediates/Hola-output.txt
        - full_pipeline/intermediates/Hello-output.txt
        - full_pipeline/intermediates/Bonjour-output.txt

      uppercased:
        - full_pipeline/intermediates/UPPER-Hola-output.txt
        - full_pipeline/intermediates/UPPER-Hello-output.txt
        - full_pipeline/intermediates/UPPER-Bonjour-output.txt

      collected: full_pipeline/intermediates/COLLECTED-batch-output.txt

      batch_report: full_pipeline/batch-report.txt

      cowpy_art: full_pipeline/cowpy-COLLECTED-batch-output.txt
    ```

Der erste Versuch schlägt bei 6 MB noch fehl, aber die Wiederholung läuft mit 12 MB (`6.MB * 2`) und ist erfolgreich. Die Pipeline wird mit allen veröffentlichten Ausgaben abgeschlossen.

!!! warning "Warnung"

    Die Konsolenausgabe enthält weiterhin eine `NOTE:`-Zeile, die den fehlgeschlagenen ersten Versuch meldet, obwohl die Pipeline insgesamt erfolgreich war: Nextflow protokolliert jede Wiederholung einzeln, aber ein wiederholter Fehler beeinflusst das Gesamtergebnis nicht.
    Prüfe die `Outputs:`-Zusammenfassung oder den Exit-Status des Befehls, um zu bestätigen, ob der Lauf tatsächlich erfolgreich war.

Weitere fortgeschrittene Wiederholungsmuster, einschließlich der Skalierung basierend auf dem aufgetretenen Fehler, findest du unter [Dynamic computing resources](https://nextflow.io/docs/latest/process.html#dynamic-task-resources) in der Nextflow-Dokumentation.

### Fazit

Du weißt, wie du eine Pipeline dazu bringst, fehlgeschlagene Aufgaben automatisch zu wiederholen, und wie du Ressourcenzuweisungen bei jeder Wiederholung mit `task.attempt` skalierst.

### Wie geht es weiter?

Weiter zu [Teil 3](./03_profiles.md), wo du lernst, wie du solche Konfigurationen in umschaltbare Profile bündelst.

---

## Zusammenfassung

In diesem Teil hast du gelernt, wie du:

- Einen Ressourcen-Profiling-Bericht erstellst und prozessspezifische Ressourcenzuweisungen festlegst
- Ressourcenanfragen mit `resourceLimits` begrenzt
- Eine fehlgeschlagene Aufgabe mit `errorStrategy` und `maxRetries` automatisch wiederholst
- Eine Ressourcenzuweisung bei jeder Wiederholung mit `task.attempt` erhöhst
