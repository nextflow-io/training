# Teil 1: An deine Rechenumgebung anpassen

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>


In [Nextflow Run](../nextflow_run/index.md) hast du die Eingaben, Parameter und Ausgaben einer Pipeline konfiguriert.
Dieser Kurs behandelt die andere Seite: wie du die Ausführung einer Pipeline an die jeweilige Rechenumgebung anpasst, ohne den Workflow-Code zu ändern.

!!! example "Szenario"

    Du hast deine Pipeline auf deinem Laptop mit Docker entwickelt und getestet.
    Jetzt musst du sie weitergeben: Eine Kollegin hat nur Conda eingerichtet, und der HPC-Cluster deiner Institution erwartet, dass Jobs über seinen eigenen Scheduler mit eigenen Ressourcenlimits eingereicht werden.
    Das alles sollte keine Änderungen an der Pipeline selbst erfordern.

Derselbe Pipeline-Code kann an all diesen Orten laufen, weil nichts davon fest im Workflow verankert ist.
Software-Paketierung, Ausführungsplattform und Ressourcenzuweisung werden alle über Konfiguration gesteuert, die über dem Code liegt – und genau das behandelt dieser Kurs: wie du dieselbe Pipeline durch Änderungen an der Konfiguration, nicht am Code, an eine neue Umgebung anpasst.

---

## 1. Eine Software-Paketierungstechnologie auswählen

In [Nextflow Run](../nextflow_run/index.md) hast du ein `conda`-Profil gesehen, das in `nextflow.config` bereits als Alternative zu Docker eingerichtet war.
Hier baust du denselben Wechsel selbst nach und siehst, was nötig ist, damit ein Prozess tatsächlich mit Conda nutzbar ist.

### 1.1. Docker deaktivieren und Conda aktivieren

Setze `docker.enabled` auf `false` und füge eine Direktive hinzu, die Conda aktiviert.

=== "Danach"

    ```groovy title="nextflow.config" linenums="1" hl_lines="1-2"
    docker.enabled = false
    conda.enabled = true
    ```

=== "Vorher"

    ```groovy title="nextflow.config" linenums="1"
    docker.enabled = true
    ```

Damit kann Nextflow Conda-Umgebungen für jeden Prozess erstellen und verwenden, für den ein Conda-Paket angegeben ist.
Der `cowpy`-Prozess hat noch keines, also fügen wir eines hinzu – vollständig über die Konfiguration.

### 1.2. Ein Conda-Paket über die Konfiguration hinzufügen

Eine `conda`-Direktive kann in der Prozessdefinition selbst gesetzt werden, genauso wie `container` bereits in `modules/cowpy.nf` gesetzt ist. Das muss aber nicht so sein: Mit `withName` kannst du sie stattdessen aus der Konfiguration heraus setzen, beschränkt auf den `cowpy`-Prozess.

=== "Danach"

    ```groovy title="nextflow.config" linenums="6" hl_lines="3-5"
    process {
        memory = 1.GB
        withName: 'cowpy' {
            conda = 'conda-forge::cowpy==1.1.5'
        }
    }
    ```

=== "Vorher"

    ```groovy title="nextflow.config" linenums="6"
    process {
        memory = 1.GB
    }
    ```

Das ersetzt nicht die `container`-Direktive, die bereits im Pipeline-Code steht, sondern fügt eine Alternative daneben hinzu – ohne diesen Code überhaupt anzufassen.

!!! tip "Tipp"

    Die [Seqera Containers](https://seqera.io/containers/)-Suche ist eine praktische Möglichkeit, den Conda-Paket-URI für ein bestimmtes Tool nachzuschlagen, auch wenn du nicht planst, daraus einen Container zu bauen.

### 1.3. Den Workflow ausführen, um Conda zu testen

```bash
nextflow run main.nf --batch conda
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [romantic_pike] revision: c3c85dec78

    executor >  local (8)
    [6d/d48030] sayHello (2)       | 3 of 3 ✔
    [f5/7a9d76] convertToUpper (1) | 3 of 3 ✔
    [1c/79b693] collectGreetings   | 1 of 1 ✔
    Creating env using conda: conda-forge::cowpy==1.1.5 [cache /workspaces/training/config-exec/work/conda/env-898314d566668b6587ad714ae06b8520]
    [bb/64b67c] cowpy              | 1 of 1 ✔

    Outputs:

      /workspaces/training/config-exec/results

      first_output:
        - full_pipeline/intermediates/Hola-output.txt
        - full_pipeline/intermediates/Bonjour-output.txt
        - full_pipeline/intermediates/Hello-output.txt

      uppercased:
        - full_pipeline/intermediates/UPPER-Hello-output.txt
        - full_pipeline/intermediates/UPPER-Bonjour-output.txt
        - full_pipeline/intermediates/UPPER-Hola-output.txt

      collected: full_pipeline/intermediates/COLLECTED-conda-output.txt

      batch_report: full_pipeline/conda-report.txt

      cowpy_art: full_pipeline/cowpy-COLLECTED-conda-output.txt
    ```

Das erzeugt dieselbe Ausgabe wie die Ausführung mit Docker, auch wenn die Mechanik im Hintergrund anders ist: Nextflow lädt das Conda-Paket herunter und baut daraus eine Umgebung, anstatt ein Container-Image zu pullen.

!!! info "Info"

    Das Erstellen einer neuen Conda-Umgebung kann beim ersten Mal etwas länger dauern als das Pullen eines Containers. Das hier verwendete Paket ist jedoch klein, daher sollte es schnell gehen.

Wechsle jetzt für den Rest dieses Kurses zurück zu Docker.

```groovy title="nextflow.config" linenums="1"
docker.enabled = true
```

??? tip "Docker und Conda kombinieren"

    Da diese Einstellungen pro Prozess gelten, kannst du sie kombinieren: Manche Prozesse verwenden Docker, andere Conda, je nachdem, was für das jeweilige Tool verfügbar ist.
    Wenn sowohl eine `container`-Direktive (im Pipeline-Code) als auch eine `conda`-Direktive (hier, aus der Konfiguration) für denselben Prozess gesetzt sind und beide Paketierungssysteme aktiviert sind, priorisiert Nextflow Container.

### Fazit

Du weißt jetzt, wie du konfigurierst, welche Software-Paketierungstechnologie ein Prozess verwenden soll, und wie du zwischen Docker und Conda wechselst.

### Wie geht es weiter?

Lerne, wie du die Ausführungsplattform änderst, die Nextflow zum Ausführen deiner Aufgaben verwendet.

---

## 2. Eine Ausführungsplattform auswählen

Jede Pipeline, die du bisher ausgeführt hast, hat den lokalen Executor verwendet: Jede Aufgabe läuft auf demselben Rechner wie Nextflow selbst.
Nextflow prüft die verfügbaren CPUs und den Arbeitsspeicher und hält Aufgaben zurück, bis genügend Ressourcen frei sind.

Der lokale Executor ist praktisch, skaliert aber nicht über einen einzelnen Rechner hinaus.
Nextflow unterstützt [viele weitere Ausführungs-Backends](https://nextflow.io/docs/latest/executor.html), darunter HPC-Scheduler (Slurm, LSF, SGE, PBS und andere) sowie Cloud-Plattformen (AWS Batch, Google Cloud Batch, Azure Batch, Kubernetes und mehr).

### 2.1. Ein anderes Backend auswählen

Der Executor wird durch eine Prozess-Direktive namens `executor` festgelegt.
Standardmäßig ist er `local`, daher gilt implizit folgendes:

```groovy title="Built-in configuration"
process {
    executor = 'local'
}
```

Um ein anderes Backend zu verwenden, setze die Direktive auf den gewünschten Executor.

```groovy title="nextflow.config"
process {
    executor = 'slurm'
}
```

!!! warning "Warnung"

    Die Trainingsumgebung ist nicht mit einem HPC-Cluster verbunden, daher kannst du das hier nicht ausprobieren.

### 2.2. Backend-spezifische Syntax wird abstrahiert

Die meisten HPC-Plattformen erfordern, dass Job-Einreichungen Ressourcenanforderungen angeben – wie CPUs, Arbeitsspeicher und einen Queue-Namen – in ihrer eigenen Syntax.
Dieselbe Anforderung für 8 CPUs und 4 GB RAM auf einer Queue namens `my-science-work` sieht je nach Scheduler völlig unterschiedlich aus.

??? abstract "Beispiele"

    ```bash title="Config for SLURM / submit using sbatch"
    #SBATCH -o /path/to/my/task/directory/my-task-1.log
    #SBATCH --no-requeue
    #SBATCH -c 8
    #SBATCH --mem 4096M
    #SBATCH -p my-science-work
    ```

    ```bash title="Config for PBS / submit using qsub"
    #PBS -o /path/to/my/task/directory/my-task-1.log
    #PBS -j oe
    #PBS -q my-science-work
    #PBS -l nodes=1:ppn=8
    #PBS -l mem=4gb
    ```

    ```bash title="Config for SGE / submit using qsub"
    #$ -o /path/to/my/task/directory/my-task-1.log
    #$ -j y
    #$ -terse
    #$ -notify
    #$ -q my-science-work
    #$ -l slots=8
    #$ -l h_rss=4096M,mem_free=4096M
    ```

Nextflow abstrahiert das alles: Du gibst standardisierte Eigenschaften wie `cpus`, `memory` und `queue` einmal an (siehe [Prozess-Direktiven](https://nextflow.io/docs/latest/reference/process.html#process-directives) für die vollständige Liste), und Nextflow übersetzt sie zur Laufzeit in die passenden Backend-spezifischen Skripte.

### 2.3. Was Nextflow tatsächlich ausführt

Diese Übersetzung ist nicht nur eine Konfigurationsdatei-Bequemlichkeit: Sie wird durch etwas Konkretes unterstützt, das du dir jetzt ansehen kannst, auch mit dem lokalen Executor.
In [Nextflow Run, Abschnitt 1.3](../nextflow_run/01_run_nextflow.md#13-explore-the-work-directory) hast du in ein Aufgabenverzeichnis unter `work/` geschaut und `.command.sh` gefunden – den genauen Befehl, den Nextflow ausgeführt hat.
Dasselbe Verzeichnis enthält auch eine Datei, die du noch nicht angeschaut hast: `.command.run`.

```bash
cat work/0a/0df4a1*/.command.run
```

??? success "Befehlsausgabe (Auszug)"

    ```console
    #!/bin/bash
    ### ---
    ### name: 'convertToUpper (3)'
    ### container: 'null'
    ### outputs:
    ### - 'UPPER-Bonjour-output.txt'
    ### ...
    set -e
    set -u
    ...
    nxf_launch() {
        /bin/bash -ue /workspaces/training/nextflow-run/work/0a/0df4a1028c2001758b1841cff92fc7/.command.sh
    }
    ...
    ```

`.command.run` ist das eigentliche Skript, das Nextflow zur Ausführung übergibt.
Es umhüllt `.command.sh` mit allem, was zum Ausführen benötigt wird: Umgebungseinrichtung, Ein-/Ausgabe-Staging und die Rückmeldung des Ergebnisses an Nextflow.
Mit dem `local`-Executor führt Nextflow dieses Skript einfach auf demselben Rechner aus.

Genau das ändert sich, wenn du einen anderen `executor` setzt.
Für einen HPC-Scheduler wie Slurm oder PBS generiert Nextflow dasselbe Wrapper-Skript, fügt den scheduler-spezifischen Header hinzu, den du in [2.2](#22-backend-specific-syntax-is-abstracted-away) gesehen hast (übersetzt aus deinen `cpus`-, `memory`- und `queue`-Einstellungen), und übergibt das Ergebnis an den eigenen Einreichungsbefehl des Schedulers, zum Beispiel `sbatch` für Slurm.
Danach fragt Nextflow den Scheduler nach dem Job-Status, anstatt einen lokalen Prozess direkt zu beobachten.
Cloud-Batch-Backends funktionieren etwas anders, da sie über API-Aufrufe statt über einen Einreichungsbefehl gesteuert werden, aber die grundlegende Idee ist dieselbe: Dasselbe Aufgaben-Skript wird ausgeführt, nur wie es gestartet und verfolgt wird, ändert sich.

### Fazit

Du weißt jetzt, wie du den Executor änderst, um verschiedene Recheninfrastrukturen anzusprechen, dass Nextflow die Backend-spezifische Einreichungssyntax abstrahiert, und was tatsächlich im Hintergrund passiert, wenn eine Aufgabe auf einem anderen Backend ausgeführt wird.

### Wie geht es weiter?

Weiter zu [Teil 2](./02_resources_and_retries.md), wo du lernst, wie du Rechenressourcen profilierst und zuweist und Aufgabenfehler mit Wiederholungsversuchen behandelst.

---

## Zusammenfassung

In diesem Teil hast du gelernt, wie du:

- Die Software-Paketierungstechnologie zwischen Docker und Conda wechselst
- Eine `conda`-Direktive zu einer Prozessdefinition hinzufügst
- Die Ausführungsplattform mit der `executor`-Direktive änderst
- Inspizierst, was Nextflow tatsächlich für eine Aufgabe generiert und ausführt, und wie sich das je nach Executor unterscheidet
