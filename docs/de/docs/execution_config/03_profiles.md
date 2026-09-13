# Teil 3: Profile verwenden, um Konfigurationen zu wechseln

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

In [Teil 1](./01_packaging_and_execution.md) und [Teil 2](./02_resources_and_retries.md) hast du einige Konfigurationsoptionen gesammelt: Software-Paketierung, Ausführungsplattform und Ressourcenzuweisungen.
In der Praxis möchtest du oft zwischen ganzen Sätzen dieser Optionen wechseln, je nachdem, wo du arbeitest – zum Beispiel ein Laptop für die Entwicklung und ein HPC-Cluster für den Produktionsbetrieb.

Nextflow ermöglicht es dir, beliebig viele [Profile](https://nextflow.io/docs/latest/config.html#profiles) mit unterschiedlichen Konfigurationen einzurichten und zur Laufzeit mit einem einzigen Flag eines (oder mehrere) auszuwählen.

Du hast bereits eines verwendet: Das `test`-Profil aus [Nextflow Run](../nextflow_run/index.md) überschreibt die Eingabeparameter mit einem kleinen, klar definierten Satz.
Jetzt erstellst du deine eigenen Infrastrukturprofile und kombinierst sie damit.

---

## 1. Profile für verschiedene Umgebungen erstellen

### 1.1. Profile einrichten

Füge zwei Profile zu `nextflow.config` hinzu: eines für die Ausführung auf einem normalen Laptop mit Docker und eines für einen universitären HPC-Cluster mit einem Slurm-Scheduler und Conda.

=== "Danach"

    ```groovy title="nextflow.config" linenums="35" hl_lines="10-19"
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
        my_laptop {
            process.executor = 'local'
            docker.enabled = true
        }
        univ_hpc {
            process.executor = 'slurm'
            conda.enabled = true
            process.resourceLimits = [
                memory: 750.GB,
                cpus: 200,
                time: 30.d
            ]
        }
    }
    ```

=== "Vorher"

    ```groovy title="nextflow.config" linenums="35"
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

Das `univ_hpc`-Profil legt außerdem Ressourcenlimits fest, da dies auf gemeinsam genutzter HPC-Infrastruktur typischerweise erforderlich ist.

### 1.2. Den Workflow mit einem Profil ausführen

Wähle ein Profil zur Laufzeit mit `-profile` aus.

```bash
nextflow run main.nf -profile my_laptop
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [nice_kimura] revision: c3c85dec78

    executor >  local (8)
    [df/b86b2f] sayHello (2)       | 3 of 3 ✔
    [b2/e2cb9c] convertToUpper (3) | 3 of 3 ✔
    [af/2ecff7] collectGreetings   | 1 of 1 ✔
    [c6/713c96] cowpy              | 1 of 1 ✔

    Outputs:

      /workspaces/training/execution-config/results

      first_output:
        - full_pipeline/intermediates/Bonjour-output.txt
        - full_pipeline/intermediates/Hello-output.txt
        - full_pipeline/intermediates/Hola-output.txt

      uppercased:
        - full_pipeline/intermediates/UPPER-Bonjour-output.txt
        - full_pipeline/intermediates/UPPER-Hello-output.txt
        - full_pipeline/intermediates/UPPER-Hola-output.txt

      collected: full_pipeline/intermediates/COLLECTED-batch-output.txt

      batch_report: full_pipeline/batch-report.txt

      cowpy_art: full_pipeline/cowpy-COLLECTED-batch-output.txt
    ```

!!! warning "Warnung"

    Das `univ_hpc`-Profil lässt sich in der Trainingsumgebung nicht ausführen, da kein Slurm-Scheduler verfügbar ist.

Wenn du weitere Einstellungen findest, die immer zusammengehören, füge sie dem entsprechenden Profil hinzu.
Du kannst auch weitere Profile erstellen, um beliebige andere Kombinationen zu gruppieren.

### 1.3. Mit mehreren Profilen ausführen

Profile schließen sich nicht gegenseitig aus.
Du kannst mehrere gleichzeitig aktivieren mit `-profile <profil1>,<profil2>`.
Kombiniere `my_laptop` mit dem `test`-Profil, das du bereits aus Nextflow Run kennst.

```bash
nextflow run main.nf -profile my_laptop,test
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [romantic_volhard] revision: c3c85dec78

    executor >  local (8)
    [22/81de4f] sayHello (3)       | 3 of 3 ✔
    [93/e16b04] convertToUpper (1) | 3 of 3 ✔
    [ca/3e6a91] collectGreetings   | 1 of 1 ✔
    [16/0cc7e3] cowpy              | 1 of 1 ✔

    Outputs:

      /workspaces/training/execution-config/results

      first_output:
        - full_pipeline/intermediates/Hola-output.txt
        - full_pipeline/intermediates/Hello-output.txt
        - full_pipeline/intermediates/Bonjour-output.txt

      uppercased:
        - full_pipeline/intermediates/UPPER-Hello-output.txt
        - full_pipeline/intermediates/UPPER-Bonjour-output.txt
        - full_pipeline/intermediates/UPPER-Hola-output.txt

      collected: full_pipeline/intermediates/COLLECTED-test-output.txt

      batch_report: full_pipeline/test-report.txt

      cowpy_art: full_pipeline/cowpy-COLLECTED-test-output.txt
    ```

Die einzelnen Dateinamen übernehmen korrekt `batch = 'test'` aus dem `test`-Profil (`COLLECTED-test-output.txt` usw.).

Wenn du Profile kombinierst, die dieselbe Option setzen, löst Nextflow den Konflikt mit dem Wert, den es zuletzt liest – also dem, der später in der Datei steht.
Wenn die widersprüchlichen Einstellungen aus völlig unterschiedlichen Konfigurationsquellen stammen, gilt die standardmäßige [Rangfolge](https://www.nextflow.io/docs/latest/config.html).

### Fazit

Du weißt jetzt, wie du Profile definierst, die infrastrukturspezifische Konfigurationen bündeln, zur Laufzeit mit `-profile` eines auswählst, mehrere Profile in einem einzigen Lauf kombinierst und wie Nextflow Konflikte auflöst, wenn mehr als ein Profil dieselbe Option setzt.

### Wie geht es weiter?

Lerne, wie du die vollständig aufgelöste Konfiguration überprüfst, bevor du irgendetwas ausführst.

---

## 2. Die aufgelöste Konfiguration prüfen

Du hast `nextflow config -profile test` bereits in [Nextflow Run](../nextflow_run/02_configure_pipeline.md) verwendet, um zu prüfen, was ein einzelnes Profil auflöst.
Dieser Befehl wird besonders nützlich, sobald du mehrere Profile kombinierst: Wie du gerade gesehen hast, kann es schwierig sein, von Hand herauszufinden, welcher Wert tatsächlich gewinnt, wenn zwei Profile dieselbe Option setzen.
Der Befehl `nextflow config` löst das alles für dich auf, ohne die Pipeline auszuführen.

### 2.1. Die Standardkonfiguration auflösen

```bash
nextflow config
```

??? success "Befehlsausgabe"

    ```groovy
    params {
       input = 'data/greetings.csv'
       batch = 'batch'
       character = 'turkey'
    }

    docker {
       enabled = true
    }

    process {
       memory = '1 GB'
       withName:cowpy {
          conda = 'conda-forge::cowpy==1.1.5'
          memory = '2 GB'
          cpus = 2
       }
    }
    ```

Das ist genau das, was gelten würde, wenn du die Pipeline ohne zusätzliche Flags ausführst.

### 2.2. Die Konfiguration mit aktivierten Profilen auflösen

Füge dieselben Profile hinzu, die du auch für einen echten Lauf verwenden würdest.

```bash
nextflow config -profile my_laptop,test
```

??? success "Befehlsausgabe"

    ```groovy
    params {
       input = 'data/greetings.csv'
       batch = 'test'
       character = 'tux'
    }

    docker {
       enabled = true
    }

    process {
       memory = '1 GB'
       withName:cowpy {
          conda = 'conda-forge::cowpy==1.1.5'
          memory = '2 GB'
          cpus = 2
       }
       executor = 'local'
    }
    ```

Der Vergleich der beiden Ausgaben bestätigt, was sich geändert hat: `params.batch`, `params.character` und `process.executor` spiegeln alle die Profile `my_laptop,test` wider.
Das ist besonders wertvoll bei Pipelines mit vielen Konfigurationsebenen, bei denen es mühsam und fehleranfällig wäre, die aufgelösten Einstellungen von Hand zu ermitteln.

### Fazit

Du weißt jetzt, wie du `nextflow config` verwendest, um die vollständig aufgelöste Konfiguration für jede Kombination von Profilen zu prüfen, bevor du irgendetwas ausführst.

### Wie geht es weiter?

Du hast die wesentlichen Grundlagen der Konfiguration von Nextflow-Pipelines abgedeckt.
Unter [Kursübersicht](next_steps.md) erfährst du, wie es von hier aus weitergeht.

---

## Zusammenfassung

In diesem Teil hast du gelernt:

- Profile zu definieren, die infrastrukturspezifische Konfigurationen bündeln
- Mehrere Profile in einem einzigen Lauf zu kombinieren und zu verstehen, wie Konflikte zwischen ihnen aufgelöst werden
- `nextflow config` zu verwenden, um die vollständig aufgelöste Konfiguration zu prüfen
