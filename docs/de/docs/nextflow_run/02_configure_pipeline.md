# Teil 2: Die Pipeline konfigurieren

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

In [Teil 1](./01_run_nextflow.md) hast du eine vollständige mehrstufige Pipeline ausgeführt, die mehrere Eingaben parallel mit Containern verarbeitet.
Jetzt schauen wir uns an, wie du das Verhalten der Pipeline mit `nextflow.config` konfigurierst: zuerst durch die Untersuchung der Konfigurationsdatei, die wir dir bereits mitgegeben haben, dann durch die Erkundung weiterer Möglichkeiten, Konfiguration bereitzustellen, und schließlich durch die Steuerung, wie und wo Ausgaben veröffentlicht werden.

---

## 1. Die Hauptkonfigurationsdatei untersuchen

Nextflow liest `nextflow.config` automatisch aus dem Arbeitsverzeichnis und wendet die Einstellungen auf jeden Lauf an.

Wir stellen dir eine Konfigurationsdatei bereit, die vier Bereiche abdeckt: Software-Paketierung, Prozesseinstellungen, Pipeline-Parameter und Ausführungsprofile.

??? full-code "nextflow.config"

    ```groovy title="nextflow.config" linenums="1"
    /*
     * Software-Paketierung
     */
    docker.enabled = true

    /*
     * Prozesseinstellungen
     */
    process {
        cpus = 1
        memory = 1.GB
    }

    /*
     * Pipeline-Parameter
     */
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }

    /*
     * Profile
     */
    profiles {
        test {
            params.input = 'data/greetings.csv'
            params.batch = 'test'
            params.character = 'tux'
        }
        conda {
            docker.enabled = false
            conda.enabled = true
        }
    }
    ```

Gehen wir jeden Bereich durch und nutzen dann Profile, indem wir die Pipeline mit einem davon ausführen.

!!! note "Hinweis"

    Diese Konfiguration deckt die lokale Ausführung auf einem einzelnen Rechner ab.
    Nextflow unterstützt auch HPC-Scheduler (SLURM, PBS, LSF) und Cloud-Executors (AWS Batch, Google Cloud Batch, Azure Batch), die alle über denselben `nextflow.config`-Mechanismus konfiguriert werden.
    Einen vollständigen Überblick über diese Optionen findest du in [Teil 1: An deine Rechenumgebung anpassen](../config_exec/01_packaging_and_execution.md) im Kurs [Configure Execution](../config_exec/index.md).

### 1.1. Software-Paketierung

Software-Paketierung beschreibt, wie Nextflow die tatsächlichen Tools bereitstellt, die deine Prozesse benötigen – ob das ein Container-Image, eine Conda-Umgebung oder etwas anderes ist.

```groovy title="nextflow.config" linenums="1"
/*
 * Software-Paketierung
 */
docker.enabled = true
```

Diese Zeile aktiviert Docker für jeden Prozess.
Jeder Prozess, der eine `container`-Direktive deklariert, läuft innerhalb des angegebenen Images.

### 1.2. Prozesseinstellungen

Denk daran, dass ein Prozess ein einzelner Schritt in deiner Pipeline ist, wie `sayHello` oder `cowpy`.
Nextflow lässt dich eine Reihe von Dingen konfigurieren, wie jeder Prozess tatsächlich ausgeführt wird: wie viel CPU und Arbeitsspeicher er bekommt, welchen Container oder welche Conda-Umgebung er verwendet, und mehr.

```groovy title="nextflow.config" linenums="6"
/*
 * Prozesseinstellungen
 */
process {
    cpus = 1
    memory = 1.GB
}
```

Das begrenzt jeden Prozess auf eine einzelne CPU und 1 GB Arbeitsspeicher.

Nextflow lässt dich auch unterschiedliche Werte für einzelne benannte Prozesse oder Gruppen von Prozessen festlegen. Wie das geht, lernst du in [Teil 2: Rechenressourcen und Fehler verwalten](../config_exec/02_resources_and_retries.md#12-set-resource-allocations-for-a-specific-process) des Kurses [Configure Execution](../config_exec/index.md).

### 1.3. Pipeline-Parameter

Parameter sind die Befehlszeileneingaben der Pipeline – dieselben `--input`-, `--batch`- und `--character`-Flags, die du bereits direkt auf der Befehlszeile gesetzt hast.
Standardwerte hier festzulegen bedeutet, dass du sie nicht jedes Mal eintippen musst. Wie du später in diesem Teil sehen wirst, gibt es noch ein paar andere Möglichkeiten, sie bereitzustellen.

```groovy title="nextflow.config" linenums="14"
/*
 * Pipeline-Parameter
 */
params {
    input = 'data/greetings.csv'
    batch = 'batch'
    character = 'turkey'
}
```

Diese Standardwerte greifen, wenn ein Parameter nicht auf der Befehlszeile angegeben wird. So funktioniert `nextflow run main.nf` ohne Flags trotzdem.

### 1.4. Profile

Profile ermöglichen es dir, eine Reihe von Einstellungen unter einem einzigen Namen zu bündeln, sodass du mit einem einzigen Flag zwischen ganzen Konfigurationen wechseln kannst, anstatt jedes Mal Werte manuell zu ändern.

```groovy title="nextflow.config" linenums="23"
/*
 * Profile
 */
profiles {
    test {
        params.input = 'data/greetings.csv'
        params.batch = 'test'
        params.character = 'tux'
    }
    conda {
        docker.enabled = false
        conda.enabled = true
    }
}
```

Das `test`-Profil überschreibt drei Parameter, um die Pipeline mit einem kleinen, klar definierten Eingabedatensatz auszuführen. Jede nf-core-Pipeline wird mit einem solchen Profil für die schnelle Validierung ausgeliefert – eine Konvention, die sich auch in deinen eigenen Pipelines lohnt.

Das `conda`-Profil wechselt die Software-Paketierung von Docker zu Conda.

Du aktivierst ein Profil, indem du `-profile <name>` auf der Befehlszeile übergibst.

Probieren wir das `test`-Profil aus.

```bash
nextflow run main.nf -profile test
```

??? success "Befehlsausgabe"

    ```console hl_lines="6 7 8 9"
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [reverent_heisenberg] revision: ce74f81996

    executor >  local (8)
    [3d/8a12c7] sayHello (3)       | 3 of 3 ✔
    [e7/e0934f] convertToUpper (2) | 3 of 3 ✔
    [1d/616569] collectGreetings   | 1 of 1 ✔
    [44/7d46cf] cowpy              | 1 of 1 ✔

    Outputs:

      /workspaces/training/nextflow-run/results

      first_output:
        - test/intermediates/Bonjour-output.txt
        - test/intermediates/Hello-output.txt
        - test/intermediates/Hola-output.txt

      uppercased:
        - test/intermediates/UPPER-Bonjour-output.txt
        - test/intermediates/UPPER-Hello-output.txt
        - test/intermediates/UPPER-Hola-output.txt

      collected: test/intermediates/COLLECTED-test-output.txt

      batch_report: test/test-report.txt

      cowpy_art: test/cowpy-COLLECTED-test-output.txt
    ```

Die Pipeline läuft mit `batch = 'test'` und `character = 'tux'`.
Schau dir `results/test/` an: Der Batch-Name ist jetzt Teil des Verzeichnispfads selbst, und die ASCII-Kunst zeigt den Tux-Pinguin statt eines Truthahns.

!!! note "Hinweis"

    Du kannst mehrere Profile gleichzeitig aktivieren und `nextflow config -profile <name>,<name>` verwenden, um das vollständig aufgelöste Ergebnis zu sehen, bevor du irgendetwas ausführst.
    Das Kombinieren von Profilen und wie Nextflow Konflikte zwischen ihnen auflöst, wird ausführlich in [Teil 3: Profile zum Wechseln von Konfigurationen verwenden](../config_exec/03_profiles.md) des Kurses [Configure Execution](../config_exec/index.md) behandelt.

### Fazit

Du weißt jetzt, was die häufigsten Elemente einer `nextflow.config`-Datei bewirken und wie du ein Profil aktivierst.

### Wie geht es weiter?

Lerne ein paar weitere Möglichkeiten kennen, Konfigurationswerte bereitzustellen, ohne die Haupt-`nextflow.config`-Datei zu ändern – nützlich für die Konfiguration einzelner Läufe und für die Weitergabe einer genauen Einstellungssammlung an andere.

---

## 2. Konfiguration über ergänzende Dateien bereitstellen

Standardwerte in `nextflow.config` festzulegen funktioniert gut für Werte, die sich selten ändern.
Nextflow bietet dir außerdem zwei gezieltere Mechanismen: eine laufspezifische Konfigurationsdatei zur Anpassung der Ausführung an eine bestimmte Umgebung und eine Parameterdatei zur Weitergabe einer genauen Eingabewertmenge an Mitarbeiter\*innen.

### 2.1. Eine laufspezifische Konfigurationsdatei verwenden

Angenommen, du verschiebst die Pipeline auf einen Rechner ohne Docker und möchtest jedem Prozess mehr Ressourcen geben.
Erstelle eine neue Konfigurationsdatei mit nur den Überschreibungen, die du brauchst:

```groovy title="custom.config" linenums="1"
process {
    cpus = 2
    memory = 2.GB
}

docker.enabled = false
conda.enabled = true
```

Übergib sie zusammen mit deiner Hauptpipeline mit `-c`:

```bash
nextflow run main.nf -c custom.config
```

??? success "Befehlsausgabe"

    ```console hl_lines="6 7 8 9"
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [exotic_cray] revision: ce74f81996

    executor >  local (8)
    [77/e02315] sayHello (1)       | 3 of 3 ✔
    [a6/ccf44b] convertToUpper (3) | 3 of 3 ✔
    [56/fd1296] collectGreetings   | 1 of 1 ✔
    [2c/205a94] cowpy              | 1 of 1 ✔

    Outputs:

      /workspaces/training/nextflow-run/results

      first_output:
        - batch/intermediates/Bonjour-output.txt
        - batch/intermediates/Hello-output.txt
        - batch/intermediates/Hola-output.txt

      uppercased:
        - batch/intermediates/UPPER-Bonjour-output.txt
        - batch/intermediates/UPPER-Hello-output.txt
        - batch/intermediates/UPPER-Hola-output.txt

      collected: batch/intermediates/COLLECTED-batch-output.txt

      batch_report: batch/batch-report.txt

      cowpy_art: batch/cowpy-COLLECTED-batch-output.txt
    ```

Nextflow fügt `custom.config` über die eigene `nextflow.config` der Pipeline zusammen, sodass jeder Prozess jetzt 2 CPUs und 2 GB Arbeitsspeicher statt der Standardwerte bekommt und über Conda statt Docker läuft.
`cowpy` ist der einzige Prozess mit einem deklarierten Conda-Paket neben seinem Container, daher ist er derjenige, für den Nextflow tatsächlich eine Umgebung erstellt:

```console
Creating env using conda: conda-forge::cowpy==1.1.5 [cache /path/to/work/conda/env-898314d566668b6587ad714ae06b8520]
```

Eine kleine Datei, die nur Ressourcenzuweisung und Paketierung überschreibt, ohne Pipeline-Parameter anzufassen, ist genau das Muster, das nf-core-Pipelines von institutionellen Konfigurationen erwarten.
Schau dir das Repository [nf-core/configs](https://github.com/nf-core/configs) für Beispiele aus der Praxis an.

Das gibt dir eine flexible Möglichkeit, eine Pipeline an eine neue Umgebung anzupassen, ohne deine normale Konfiguration zu verändern.

### 2.2. Eine Parameterdatei verwenden

Angenommen, du musst stattdessen eine genaue Menge von Laufparametern mit Mitarbeiter\*innen teilen oder sie für eine Veröffentlichung festhalten.

Nextflow erlaubt dir, [Parameterdateien](https://nextflow.io/docs/latest/config.html#parameter-file) im YAML- oder JSON-Format bereitzustellen. Das ist eine einfachere Möglichkeit, eine genaue, reproduzierbare Wertemenge zu verteilen.

Eine Parameterdatei namens `test-params.yaml` ist bereits in deinem Arbeitsverzeichnis vorhanden:

```yaml title="test-params.yaml" linenums="1"
input: "data/greetings.csv"
batch: "yaml"
character: "stegosaurus"
```

Die Syntax verwendet Doppelpunkte (`:`), anstelle der Gleichheitszeichen (`=`) aus `nextflow.config`, da diese Datei reines YAML und kein Groovy ist.

!!! info "Info"

    Eine JSON-Version, `test-params.json`, ist ebenfalls vorhanden. Probiere sie gerne selbst aus; die Syntax zum Übergeben ist identisch.

Übergib die Datei mit `-params-file`:

```bash
nextflow run main.nf -params-file test-params.yaml
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [sharp_faraday] revision: ce74f81996

    executor >  local (8)
    [1c/9ff63e] sayHello (1)       | 3 of 3 ✔
    [3b/bb5691] convertToUpper (2) | 3 of 3 ✔
    [cd/2c1f6e] collectGreetings   | 1 of 1 ✔
    [89/c333bc] cowpy              | 1 of 1 ✔

    Outputs:

      /workspaces/training/nextflow-run/results

      first_output:
        - yaml/intermediates/Bonjour-output.txt
        - yaml/intermediates/Hello-output.txt
        - yaml/intermediates/Hola-output.txt

      uppercased:
        - yaml/intermediates/UPPER-Bonjour-output.txt
        - yaml/intermediates/UPPER-Hello-output.txt
        - yaml/intermediates/UPPER-Hola-output.txt

      collected: yaml/intermediates/COLLECTED-yaml-output.txt

      batch_report: yaml/yaml-report.txt

      cowpy_art: yaml/cowpy-COLLECTED-yaml-output.txt
    ```

??? abstract "Dateiinhalt"

    ```console title="results/yaml/cowpy-COLLECTED-yaml-output.txt"
     _________
    / BONJOUR \
    | HOLA    |
    \ HELLO   /
     ---------
    \                             .       .
     \                           / `.   .' "
      \                  .---.  <    > <    >  .---.
       \                 |    \  \ - ~ ~ - /  /    |
             _____          ..-~             ~-..-~
            |     |   \~~~\.'                    `./~~~/
           ---------   \__/                        \__/
          .'  O    \     /               /       \  "
         (_____,    `._.'               |         }  \/~~~/
          `----.          /       }     |        /    \__/
                `-.      |       /      |       /      `. ,~~|
                    ~-.__|      /_ - ~ ^|      /- _      `..-'
                        |     /        |     /     ~-.     `-. _  _  _
                        |_____|        |_____|         ~ - . _ _ _ _ _>
    ```

Eine Parameterdatei ist besonders wertvoll, wenn eine Pipeline mehr als eine Handvoll Parameter hat: Sie ermöglicht es dir, alle auf einmal bereitzustellen, ohne eine lange Befehlszeile oder Änderungen am Workflow-Skript, und lässt sich leicht zusammen mit deinen Ergebnissen weitergeben.

### Fazit

Du kennst jetzt zwei weitere Möglichkeiten, Konfiguration bereitzustellen: eine laufspezifische Konfigurationsdatei zur Anpassung der Ausführung an eine neue Umgebung und eine Parameterdatei zur Weitergabe genauer, reproduzierbarer Eingabewerte.

### Wie geht es weiter?

Lerne, wie du steuerst, wie und wo die Ausgaben deiner Pipeline veröffentlicht werden.

---

## 3. Pipeline-Ausgaben verwalten

Pipeline-Autor\*innen entscheiden im Code, wie Ausgaben organisiert werden. Du musst diesen Code jedoch nicht anfassen, um zu steuern, wo sie landen oder wie sie dorthin gelangen.
Nextflow bietet dir dafür Möglichkeiten auf Konfigurationsebene: ein Basisausgabeverzeichnis festlegen und wählen, ob Dateien kopiert oder verlinkt werden.

### 3.1. Das Ausgabeverzeichnis anpassen

Standardmäßig veröffentlicht Nextflow Ausgaben unter `results/`.
Mit `-output-dir` (oder der Kurzform `-o`) kannst du einen anderen Ort angeben:

```bash
nextflow run main.nf -output-dir outputs
```

??? success "Befehlsausgabe"

    ```console hl_lines="6 7 8 9"
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [serene_kimura] revision: ce74f81996

    executor >  local (8)
    [31/df5c15] sayHello (3)       | 3 of 3 ✔
    [08/8bb2e5] convertToUpper (1) | 3 of 3 ✔
    [e5/5814da] collectGreetings   | 1 of 1 ✔
    [cf/8ab8c6] cowpy              | 1 of 1 ✔

    Outputs:

      /workspaces/training/nextflow-run/outputs

      first_output:
        - batch/intermediates/Bonjour-output.txt
        - batch/intermediates/Hello-output.txt
        - batch/intermediates/Hola-output.txt

      uppercased:
        - batch/intermediates/UPPER-Bonjour-output.txt
        - batch/intermediates/UPPER-Hello-output.txt
        - batch/intermediates/UPPER-Hola-output.txt

      collected: batch/intermediates/COLLECTED-batch-output.txt

      batch_report: batch/batch-report.txt

      cowpy_art: batch/cowpy-COLLECTED-batch-output.txt
    ```

??? abstract "Verzeichnisinhalt"

    ```console
    outputs/batch
    ├── batch-report.txt
    ├── cowpy-COLLECTED-batch-output.txt
    └── intermediates
        ├── Bonjour-output.txt
        ├── COLLECTED-batch-output.txt
        ├── Hello-output.txt
        ├── Hola-output.txt
        ├── UPPER-Bonjour-output.txt
        ├── UPPER-Hello-output.txt
        └── UPPER-Hola-output.txt
    ```

Die Ausgaben landen jetzt unter `outputs/batch/` statt im Standard `results/batch/`.
Der eigene Code der Pipeline bestimmt weiterhin die Struktur innerhalb dieses Basisverzeichnisses, wie die Unterverzeichnisse `batch/` und `intermediates/`. `-output-dir` steuert nur, wo diese Struktur beginnt.

`-output-dir` ist eigentlich nur eine Befehlszeilen-Abkürzung für die Konfigurationsoption `outputDir`, die überall dort verwendet werden kann, wo Konfiguration möglich ist: direkt in `nextflow.config`, innerhalb eines Profils oder in einer `-c`-Overlay-Datei wie der, die du früher in diesem Teil verwendet hast.
Dieses Beispiel zeigt dieselbe Einstellung direkt in `nextflow.config` statt auf der Befehlszeile:

```groovy title="nextflow.config"
outputDir = 'outputs'
```

Die vollständige Liste der Orte, an denen eine Konfigurationsoption wie diese stehen kann, findest du unter [Configuration file](https://nextflow.io/docs/latest/config.html) in der Nextflow-Referenz.

### 3.2. Wählen, wie Ausgaben veröffentlicht werden

Standardmäßig veröffentlicht Nextflow Ausgaben als Symlinks, die auf die Speicherorte der Ausgaben unter `work/` zeigen, keine echten Kopien:

```console
$ ls -l results/batch/intermediates/Hello-output.txt
lrwxr-xr-x  ...  Hello-output.txt -> /workspaces/training/nextflow-run/work/b7/b8c4d1.../Hello-output.txt
```

Pipeline-Autor\*innen können den Veröffentlichungsmodus für jeden einzelnen Prozess im Workflow-Code auf `'copy'` oder `'move'` setzen.
Typischerweise tun sie das für die finalen Ausgaben der Pipeline, während das Standard-`'symlink'`-Verhalten für Zwischendateien beibehalten wird, die gelöscht werden können, sobald die gesamte Pipeline durchgelaufen ist.

Das vermeidet doppelte Daten auf der Festplatte, bedeutet aber, dass du die Aufgabenverzeichnisse unter `work/` nicht löschen kannst, ohne den Link zu unterbrechen und die Möglichkeit zu verlieren, `-resume` zu verwenden.
Wenn du möchtest, dass alle Ausgabedateien richtig kopiert werden, setze [`workflow.output.mode`](https://nextflow.io/docs/latest/reference/config.html#workflow) in deiner Pipeline-Konfiguration auf `'copy'`. (Anders als `-output-dir` gibt es dafür kein Befehlszeilen-Flag; das geht nur über die Konfiguration.)

Probiere es aus, indem du es in `nextflow.config` setzt:

=== "Danach"

    ```groovy title="nextflow.config" hl_lines="7"
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }

    workflow.output.mode = 'copy'
    ```

=== "Vorher"

    ```groovy title="nextflow.config"
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }
    ```

Führe dann die Pipeline aus und ändere den Batch-Namen, damit du den Unterschied in den Ausgaben sehen kannst:

```bash
nextflow run main.nf --batch withmode
```

??? success "Befehlsausgabe"

    ```console hl_lines="6 7 8 9"
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [angry_noether] revision: ce74f81996

    executor >  local (8)
    [41/2b478d] sayHello (3)       | 3 of 3 ✔
    [bf/dd2840] convertToUpper (2) | 3 of 3 ✔
    [ea/364e97] collectGreetings   | 1 of 1 ✔
    [10/76fe7b] cowpy              | 1 of 1 ✔

    Outputs:

      /workspaces/training/nextflow-run/results

      first_output:
        - withmode/intermediates/Bonjour-output.txt
        - withmode/intermediates/Hello-output.txt
        - withmode/intermediates/Hola-output.txt

      uppercased:
        - withmode/intermediates/UPPER-Bonjour-output.txt
        - withmode/intermediates/UPPER-Hello-output.txt
        - withmode/intermediates/UPPER-Hola-output.txt

      collected: withmode/intermediates/COLLECTED-withmode-output.txt

      batch_report: withmode/withmode-report.txt

      cowpy_art: withmode/cowpy-COLLECTED-withmode-output.txt
    ```

Schau dir eine der Ausgabedateien wie zuvor an:

```console
$ ls -l results/withmode/intermediates/Hello-output.txt
-rw-r--r--  ...  Hello-output.txt
```

Jetzt ist es eine echte, unabhängige Datei, die auch dann verfügbar bleibt, wenn `work/` bereinigt wird.

!!! warning "Warnung"

    Die Einstellung `workflow.output.mode` füllt nur einen Standardwert für Ausgaben aus, für die im Pipeline-Code noch kein Modus gesetzt ist.
    Sie kann einen vom Autor fest kodierten Modus nicht überschreiben, egal was du einstellst.

### Fazit

Du weißt jetzt, wie du das Basisausgabeverzeichnis anpasst und zwischen kopierten und verlinkten Ausgaben wählst – alles ohne den Code der Pipeline anzufassen.

### Wie geht es weiter?

Weiter zu [Teil 3](./03_manage_executions.md), wo du lernst, wie du den Verlauf vergangener Läufe einsehen, Ausführungsberichte erstellen und alte Arbeitsverzeichnisse bereinigen kannst.

---

## Zusammenfassung

In diesem Teil hast du gelernt:

- Das Verhalten der Pipeline mit `nextflow.config` und Profilen konfigurieren
- Konfiguration über eine laufspezifische Konfigurationsdatei oder eine Parameterdatei bereitstellen
- Das Ausgabeverzeichnis anpassen und zwischen kopierten und verlinkten Ausgaben wählen
