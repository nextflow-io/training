# Teil 2: Pipeline-Ausführung konfigurieren

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

In [Teil 1](./01_run_demo.md) hast du die nf-core/demo Pipeline gefunden und mit ihrem Test-Profil ausgeführt.
Jetzt schauen wir uns an, wie du die Pipeline-Ausführung konfigurierst: Parameter setzen, Validierung verstehen und Ressourcenzuweisung sowie Tool-Argumente anpassen.

Wie in [Hello Config](../hello_nextflow/06_hello_config.md) erklärt, möchten wir die Daten und die Ausführungsweise unserer Pipeline ändern können, ohne den Pipeline-Code selbst anzupassen.
Dafür unterstützt Nextflow mehrere Möglichkeiten zur Steuerung der Pipeline-Konfiguration, was zunächst etwas überwältigend wirken kann.

Das nf-core-Projekt legt Konventionen für die Organisation von Konfigurationselementen fest und unterscheidet auf oberster Ebene zwei Arten von Konfiguration: **Pipeline-Parameter** und **Konfiguration** im engeren Sinne.

- **Pipeline-Parameter** (über das `params`-System gesetzt) umfassen typischerweise Eingabedateien, Flags für das Tool-Verhalten und Analyseparameter.
- **Konfiguration** im engeren Sinne bezieht sich auf die Logistik der Pipeline-Ausführung, also den Executor, die Zuweisung von Rechenressourcen usw.

<figure class="excalidraw">
    --8<-- "docs/en/docs/nfcore_use/img/params_vs_config.excalidraw.svg"
</figure>

Fangen wir mit den Pipeline-Parametern an und schauen uns danach die Konfiguration im engeren Sinne an.

---

## 1. Pipeline-Parameter

Für alle nf-core Pipelines kannst du eine vollständige Liste der Pipeline-Parameter direkt über die Kommandozeile abrufen, indem du das Flag `--help` verwendest – das selbst ein Pipeline-Parameter ist.

### 1.1. Parameterliste mit `--help` abrufen

Führe den Hilfe-Befehl für die Demo-Pipeline aus:

```bash
nextflow run nf-core/demo --help
```

??? success "Befehlsausgabe"

    ```console
     N E X T F L O W   ~  version 26.04.4

    Launching `https://github.com/nf-core/demo` [adoring_meucci] revision: 32893afef8 [master]


    ------------------------------------------------------
                                            ,--./,-.
            ___     __   __   __   ___     /,-._.--~'
      |\ | |__  __ /  ` /  \ |__) |__         }  {
      | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                            `._,._,'
      nf-core/demo 1.2.0
    ------------------------------------------------------
    Typical pipeline command:

      nextflow run nf-core/demo -profile <docker/singularity/.../institute> --input samplesheet.csv --outdir <OUTDIR>


    Input/output options
      --input                       [string] Path to a metadata file containing information about the samples in the experiment.
      --outdir                      [string] The output directory where the results will be saved. You have to use absolute paths to storage on Cloud infrastructure.
      --email                       [string] Email address for completion summary.
      --multiqc_title               [string] MultiQC report title. Printed as page header, used for filename if not otherwise specified.

    Reference genome options
      --genome                      [string] Name of iGenomes reference.
      --fasta                       [string] Path to FASTA genome file.

    Process skipping options
      --skip_trim                   [boolean] Skip trimming fastq files with seqtk

    Generic options
      --multiqc_methods_description [string]          Custom MultiQC yaml file containing HTML including a methods description.
      --help                        [boolean, string] Display the help message.
      --help_full                   [boolean]         Display the full detailed help message.
      --show_hidden                 [boolean]         Display hidden parameters in the help message (only works when --help or --help_full are provided).
     !! Hiding 19 param(s), use the `--showHidden` parameter to show them !!
    ------------------------------------------------------

    * The pipeline
        https://doi.org/10.5281/zenodo.12192442

    * The nf-core framework
        https://doi.org/10.1038/s41587-020-0439-x

    * Software dependencies
        https://github.com/nf-core/demo/blob/master/CITATIONS.md
    ```

Wie du siehst, gruppiert die Ausgabe die Parameter in Kategorien (Input/output options, Reference genome options usw.) mit Typ und Beschreibung für jeden Parameter.

Diese Kategorisierung wird durch eine Schema-Datei festgelegt, auf die wir weiter unten eingehen.
Bei einfachen Nextflow-Pipelines funktioniert `--help` nur, wenn der Entwickler bzw. die Entwicklerin es manuell implementiert hat.

!!! tip "Tipp"

    Verwende `--help --show_hidden`, um zusätzliche Parameter zu sehen, die standardmäßig ausgeblendet sind, wie z. B. `--publish_dir_mode` oder `--monochrome_logs`.

### 1.2. Parameterwerte setzen

Wie in [Hello Config](../hello_nextflow/06_hello_config.md) beschrieben, kannst du Parameterwerte auf der Kommandozeile mit `--param_name` setzen oder eine Reihe von Parametern in einer YAML-Datei sammeln und diese mit `-params-file` übergeben.
Beide Ansätze funktionieren bei nf-core Pipelines gleich.

Um zum Beispiel den Trimming-Schritt zu überspringen, setzen wir den booleschen Parameter `skip_trim` auf `true`.
In deinem Arbeitsverzeichnis liegt eine Params-Datei namens `my_params.yml` mit diesem Wert:

```yaml title="my_params.yml"
skip_trim: true
```

Übergib sie mit `-params-file`:

```bash
nextflow run nf-core/demo -profile test,docker --outdir demo-results-notrim -params-file my_params.yml
```

??? success "Befehlsausgabe"

    ```console
     N E X T F L O W   ~  version 26.04.4

    Launching `https://github.com/nf-core/demo` [focused_heisenberg] revision: 32893afef8 [master]


    ------------------------------------------------------
                                            ,--./,-.
            ___     __   __   __   ___     /,-._.--~'
      |\ | |__  __ /  ` /  \ |__) |__         }  {
      | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                            `._,._,'
      nf-core/demo 1.2.0
    ------------------------------------------------------

    Input/output options
      input                     : https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv
      outdir                    : demo-results-notrim

    Process skipping options
      skip_trim                 : true

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      trace_report_suffix       : 2026-07-03_22-08-47

    Core Nextflow options
      revision                  : master
      runName                   : focused_heisenberg
      containerEngine           : docker
      launchDir                 : /workspaces/training/nfcore-use
      workDir                   : /workspaces/training/nfcore-use/work
      projectDir                : /workspaces/.nextflow/assets/.repos/nf-core/demo/clones/32893afef8076a03a2767a020b3f0cab2e0b40b2
      userName                  : root
      profile                   : test,docker
      configFiles               : /workspaces/.nextflow/assets/.repos/nf-core/demo/clones/32893afef8076a03a2767a020b3f0cab2e0b40b2/nextflow.config

    !! Only displaying parameters that differ from the pipeline defaults !!
    ------------------------------------------------------

    * The pipeline
        https://doi.org/10.5281/zenodo.12192442

    * The nf-core framework
        https://doi.org/10.1038/s41587-020-0439-x

    * Software dependencies
        https://github.com/nf-core/demo/blob/master/CITATIONS.md

    executor >  local (5)
    [7a/f3599e] NFCORE_DEMO:DEMO:FASTQC (SAMPLE3_SE) | 3 of 3 ✔
    [b0/2f0bdc] NFCORE_DEMO:DEMO:COWPY               | 1 of 1 ✔
    [c3/3c2278] NFCORE_DEMO:DEMO:MULTIQC (demo)      | 1 of 1 ✔
    -[nf-core/demo] Pipeline completed successfully-
    ```

Der `SEQTK_TRIM`-Prozess erscheint nicht mehr in der Ausgabe.

!!! warning "Warnung: Wichtige Einschränkungen bei Parameter-Eingaben"

    **Boolesche Parameter auf der Kommandozeile setzen**

    Ab Nextflow Version 26.04 werden alle auf der Kommandozeile übergebenen Werte als Strings typisiert.
    Bei einem booleschen Parameter wie `skip_trim` wird das Übergeben als einfaches Flag (`--skip_trim`) oder als `--skip_trim true` als **String** `"true"` ausgewertet, was die Schema-Validierung fehlschlagen lässt:

    ```console
    * --skip_trim (true): Value is [string] but should be [boolean]
    ```

    Um einen booleschen Parameter auf einen echten `true`/`false`-Wert zu setzen, verwende wie oben gezeigt eine `-params-file` oder setze ihn in einer Konfigurationsdatei.
    String-, Integer- und Dateipfad-Parameter sind davon nicht betroffen und können weiterhin direkt auf der Kommandozeile gesetzt werden.
    Dieser Kurs verwendet dieses Muster durchgehend für boolesche Parameter.

    **Eigene Konfigurationsdateien verwenden**

    Obwohl es technisch möglich ist, Pipeline-Parameter in einer eigenen Konfigurationsdatei zu setzen, die mit `-c` übergeben wird, überschreibt dies möglicherweise nicht die bereits in der `nextflow.config` der Pipeline gesetzten Standardwerte – abhängig von Nextflows Konfigurationsprioritätsregeln.
    Die Verwendung von `--param_name` auf der Kommandozeile oder `-params-file` ist zuverlässiger, da diese immer Vorrang haben.

    Als Faustregel gilt: Wenn ein Parameter in der `--help`-Ausgabe erscheint, setze ihn über die Kommandozeile oder eine Params-Datei statt über eine Konfigurationsdatei.

### 1.3. Parameter-Validierung

Interessant zu wissen: Der `--help`-Befehl funktioniert bei allen nf-core Pipelines, weil das nf-core-Projekt von Entwickler\*innen verlangt, alle Pipeline-Parameter formal in einer JSON-Schema-Datei (`nextflow_schema.json`) zu definieren.
Dieses Schema erfasst Typ, Beschreibung, Standardwert und Gruppierung jedes Parameters.

Neben der Bereitstellung der `--help`-Ausgabe ermöglicht die Schema-Datei auch eine automatische Validierung beim Start.
Das bedeutet, dass Nextflow prüfen kann, ob jeder übergebene Parameter existiert und einen geeigneten Wert erhalten hat (richtiger Typ, innerhalb des erlaubten Wertebereichs usw.).

Wir gehen darauf in [Abschnitt zur Eingabe-Validierung](../nfcore_build/04_input_validation.md) genauer ein, aber du kannst es bereits in Aktion sehen, indem du der Demo-Pipeline ungültige Parameter-Eingaben gibst.

#### 1.3.1. Nicht erkannte Parameter

Versuche, einen nicht existierenden Parameter zu übergeben:

```bash
nextflow run nf-core/demo -profile test,docker --outdir demo-results --foobar "invalid"
```

Die Konsolenausgabe enthält eine Warnung:

```console
WARN: The following invalid input values have been detected:

* --foobar: invalid
```

Die Pipeline läuft weiter, aber die Warnung macht dich sofort darauf aufmerksam, dass `--foobar` kein bekannter Parameter ist.
Das soll deine Aufmerksamkeit auf nicht-kritische Tippfehler lenken, wie z. B. `--outDir` statt `--outdir`, was dir helfen kann, Zeit und Rechenressourcen zu sparen.

#### 1.3.2. Ungültige Parameterwerte

Die Validierung prüft auch Parameter**werte**.
Der Parameter `--skip_trim` ist ein boolesches Flag. Das Übergeben eines String-Werts lässt die Pipeline sofort fehlschlagen:

```bash
nextflow run nf-core/demo -profile test,docker --outdir demo-results --skip_trim yes
```

```console
ERROR ~ Validation of pipeline parameters failed!

 -- Check '.nextflow.log' file for details
The following invalid input values have been detected:

* --skip_trim (yes): Value is [string] but should be [boolean]

 -- Check script '/workspaces/.nextflow/assets/.repos/nf-core/demo/clones/32893afef8076a03a2767a020b3f0cab2e0b40b2/subworkflows/nf-core/utils_nfschema_plugin/main.nf' at line: 68 or see '.nextflow.log' file for more details
```

Die Pipeline stoppt, bevor irgendein Prozess ausgeführt wird, und bewahrt dich so vor einer fehlgeschlagenen oder falschen Ausführung.
Wie in [1.2](#12-set-parameter-values) erwähnt, sollten boolesche Parameter auf einen echten `true`/`false`-Wert in einer Params-Datei gesetzt werden, statt sie auf der Kommandozeile zu übergeben, da Kommandozeilenwerte als Strings typisiert werden.

### 1.4. Eingabe-Validierung

Dieselbe Validierungslogik kann auch verwendet werden, um die Gültigkeit von Eingabedateien zu prüfen.
Wenn eine Pipeline beispielsweise ein Samplesheet als Hauptdateneingabe erwartet (was bei vielen, wenn nicht den meisten nf-core Pipelines der Fall ist), kann der Entwickler bzw. die Entwicklerin ein Eingabe-Schema bereitstellen (getrennt vom Parameter-Schema), das beschreibt, wie die Eingabedatei strukturiert sein soll.

Zur Laufzeit kann Nextflow dann prüfen, ob die bereitgestellte Eingabedatei gültig ist.

Wir gehen auch darauf in [Abschnitt zur Eingabe-Validierung](../nfcore_build/04_input_validation.md) genauer ein, aber du kannst es bereits in Aktion sehen, indem du der Demo-Pipeline ein ungültiges Eingabe-Samplesheet gibst.

Die `nf-core/demo` Pipeline erwartet eine CSV-Datei mit den Spalten `sample`, `fastq_1` und `fastq_2`.
Dies ist in einer Schema-Datei (`assets/schema_input.json`) definiert, die die erwartete Struktur, Spaltentypen und Einschränkungen festlegt.

??? abstract "Schema-Datei für Eingaben"

    ```json title="assets/schema_input.json"
    {
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        "$id": "https://raw.githubusercontent.com/nf-core/demo/master/assets/schema_input.json",
        "title": "nf-core/demo pipeline - params.input schema",
        "description": "Schema for the file provided with params.input",
        "type": "array",
        "items": {
            "type": "object",
            "properties": {
                "sample": {
                    "type": "string",
                    "pattern": "^\\S+$",
                    "errorMessage": "Sample name must be provided and cannot contain spaces",
                    "meta": ["id"]
                },
                "fastq_1": {
                    "type": "string",
                    "format": "file-path",
                    "exists": true,
                    "pattern": "^([\\S\\s]*\\/)?[^\\s\\/]+\\.f(ast)?q\\.gz$",
                    "errorMessage": "FastQ file for reads 1 must be provided, cannot contain spaces and must have extension '.fq.gz' or '.fastq.gz'"
                },
                "fastq_2": {
                    "type": "string",
                    "format": "file-path",
                    "exists": true,
                    "pattern": "^([\\S\\s]*\\/)?[^\\s\\/]+\\.f(ast)?q\\.gz$",
                    "errorMessage": "FastQ file for reads 2 cannot contain spaces and must have extension '.fq.gz' or '.fastq.gz'"
                }
            },
            "required": ["sample", "fastq_1"]
        }
    }
    ```

Das Schema legt fest, dass `sample` und `fastq_1` erforderlich sind, während `fastq_2` optional ist (unterstützt sowohl Paired-End- als auch Single-End-Daten).
Dateipfade werden auf Existenz und Erweiterungsmuster geprüft.

Um dies zu demonstrieren, stellen wir ein fehlerhaftes Samplesheet namens `malformed_samplesheet.csv` in deinem Arbeitsverzeichnis bereit:

```csv title="malformed_samplesheet.csv"
sample,fastq_2
SAMPLE1,/not/a/real/file.fastq.gz
```

Diesem Samplesheet fehlt die erforderliche Spalte `fastq_1` und es enthält einen nicht existierenden Dateipfad in `fastq_2`.

Führe die Demo-Pipeline mit `malformed_samplesheet.csv` als Eingabe aus:

```bash
nextflow run nf-core/demo -profile test,docker --outdir demo-results --input malformed_samplesheet.csv
```

```console
ERROR ~ Validation of pipeline parameters failed!

 -- Check '.nextflow.log' file for details
The following invalid input values have been detected:

* --input (malformed_samplesheet.csv): Validation of file failed:
    -> Entry 1: Error for field 'fastq_2' (/not/a/real/file.fastq.gz): the file or directory '/not/a/real/file.fastq.gz' does not exist (FastQ file for reads 2 cannot contain spaces and must have extension '.fq.gz' or '.fastq.gz')
    -> Entry 1: Missing required field(s): fastq_1

 -- Check script '/workspaces/.nextflow/assets/.repos/nf-core/demo/clones/32893afef8076a03a2767a020b3f0cab2e0b40b2/subworkflows/nf-core/utils_nfschema_plugin/main.nf' at line: 68 or see '.nextflow.log' file for more details
```

Wie du siehst, schlägt die Pipeline sofort fehl und meldet **alle** Validierungsfehler auf einmal.
nf-schema stoppt nicht beim ersten Fehler – es sammelt alle Probleme und listet sie zusammen auf, damit du alles auf einmal beheben kannst, statt Probleme einzeln zu entdecken.

Jeder Fehler identifiziert den genauen Eintrag und das Feld, das das Problem verursacht hat. So kannst du dein Samplesheet korrigieren und die Pipeline dann mit der Gewissheit neu starten, dass sie nicht an einem späteren Punkt fehlschlägt, wenn Nextflow tatsächlich auf den Dateipfad zugreift.

Für Entwickler\*innen wird all das in [Teil 4 von Build with nf-core](../nfcore_build/04_input_validation.md) ausführlicher behandelt.

### Fazit

Du weißt jetzt, wie du mit `--help` eine vollständige Liste der Parameter einer Pipeline abrufst, sie über die Kommandozeile oder eine Params-Datei setzt und wie Nextflow sowohl Parameterwerte als auch Eingabedateien anhand der Pipeline-Schemas validiert.

### Wie geht es weiter?

Lerne mehr über die andere Art der Konfiguration: wie die Pipeline ausgeführt wird, einschließlich Ressourcenzuweisung und Tool-Argumente.

---

## 2. Konfiguration

Konfiguration im engeren Sinne steuert **wie** die Pipeline ausgeführt wird: Ressourcenzuweisung, Tool-spezifische Argumente, wo Aufgaben ausgeführt werden und welches Software-Packaging-System verwendet wird.

nf-core Pipelines enthalten Standardkonfiguration in `nextflow.config` und dem Verzeichnis `conf/`.
Bevor du etwas überschreibst, ist es hilfreich zu wissen, wo die Standardwerte zu finden sind.

### 2.1. Konfigurationsdateien erkunden

Du hast in [Teil 1](./01_run_demo.md) bereits gesehen, dass der Pipeline-Quellcode unter `$NXF_HOME/assets` liegt.
Verwende den `pipelines`-Symlink, den du in [Teil 1](./01_run_demo.md) erstellt hast, und liste die Konfigurationsdateien auf, um zu sehen, was verfügbar ist:

```bash
ls pipelines/nf-core/demo/conf/
```

```console
base.config
containers_conda_lock_files_amd64.config
containers_conda_lock_files_arm64.config
containers_docker_amd64.config
containers_docker_arm64.config
containers_singularity_https_amd64.config
containers_singularity_https_arm64.config
containers_singularity_oras_amd64.config
containers_singularity_oras_arm64.config
igenomes.config
igenomes_ignored.config
modules.config
test.config
test_full.config
```

<figure class="excalidraw">
--8<-- "docs/en/docs/nfcore_use/img/nfcore_config_files.excalidraw.svg"
</figure>

Die wichtigsten Konfigurationsdateien sind:

- **`conf/base.config`**: Definiert Ressourcen-Labels (`process_low`, `process_medium`, `process_high`), die Prozessen CPUs, Arbeitsspeicher und Zeit zuweisen. Wenn ein Prozess mehr Ressourcen als erwartet verwendet, kommen die Standardwerte von hier.
- **`conf/modules.config`**: Legt prozessspezifische Tool-Argumente (`ext.args`) und Einstellungen für die Ausgabe-Veröffentlichung (`publishDir`) fest. Öffne diese Datei, um zu sehen, welche Argumente jedes Tool standardmäßig erhält.
- **`conf/test.config`**: Das Test-Profil, das du in [Teil 1](./01_run_demo.md) verwendet hast. Es begrenzt Ressourcen über `resourceLimits` und setzt ein Test-Samplesheet. Wird mit `-profile test` aktiviert.
  Es gibt auch eine `conf/test_full.config` für die Ausführung mit einem vollständigen Test-Datensatz, was für Benchmarking nützlich ist.

Die zentrale `nextflow.config` lädt all das oben Genannte und setzt die entsprechenden Standardwerte für alles.

Wenn du Einstellungen aus diesen Dateien ändern möchtest, bearbeite keine dieser Dateien direkt.
Erstelle stattdessen deine eigene Konfigurationsdatei und übergib sie mit `-c`.
Die von dir angegebenen Werte überschreiben die in den anderen Dateien gesetzten Standardwerte.

Lass uns das in der Praxis ausprobieren.

### 2.2. Prozessressourcen und Tool-Argumente anpassen

nf-core Module unterstützen zwei gängige Arten der Konfigurationsüberschreibung: **Ressourcenzuweisung** (CPUs, Arbeitsspeicher, Zeit) und **Tool-Argumente** über `ext.args`.

Viele Kommandozeilen-Tools haben Argumente, die nicht häufig genug verwendet werden, um als Pipeline-Parameter verfügbar gemacht zu werden.
Die `ext.args`-Konvention ermöglicht es dir, diese Argumente über eine Konfigurationsdatei an das zugrunde liegende Tool zu übergeben.

Die in deinem Arbeitsverzeichnis bereitgestellte Datei `custom.config` demonstriert beide Überschreibungen:

```groovy title="custom.config" linenums="1"
process {
    withName: 'FASTQC' {
        cpus = 2
        memory = 4.GB
    }
    withName: 'SEQTK_TRIM' {
        ext.args = '-b 5'
    }
}
```

Der erste Block überschreibt die Ressourcenzuweisung für `FASTQC`.
Standardmäßig verwendet `FASTQC` das Label `process_medium` aus `base.config`, das 6 CPUs und 36 GB Arbeitsspeicher zuweist; hier begrenzen wir es auf 2 CPUs und 4 GB.

Der zweite Block übergibt ein zusätzliches Argument an `SEQTK_TRIM` über `ext.args`.
Das Flag `-b 5` weist `seqtk trimfq` an, zusätzlich zum Qualitäts-Trimming 5 Basen vom Anfang jedes Reads zu entfernen.

Führe die Pipeline mit dieser Konfiguration aus:

```bash
nextflow run nf-core/demo -profile test,docker --outdir demo-results-custom -c custom.config
```

??? success "Befehlsausgabe"

    ```console
    executor >  local (8)
    [95/b32876] NFCORE_DEMO:DEMO:FASTQC (SAMPLE1_PE)     | 3 of 3 ✔
    [17/428668] NFCORE_DEMO:DEMO:SEQTK_TRIM (SAMPLE1_PE) | 3 of 3 ✔
    [cf/85991a] NFCORE_DEMO:DEMO:COWPY                   | 1 of 1 ✔
    [3c/94a7a0] NFCORE_DEMO:DEMO:MULTIQC (demo)          | 1 of 1 ✔
    -[nf-core/demo] Pipeline completed successfully-
    ```

Das Flag `-c` fügt deine Konfiguration zur integrierten Konfiguration der Pipeline hinzu.

Um zu überprüfen, ob die `ext.args`-Überschreibung wirksam war, suche den Hash des `SEQTK_TRIM`-Arbeitsverzeichnisses aus der Ausführungsausgabe (z. B. `work/17/428668...`) und prüfe die darin enthaltene Datei `.command.sh`:

```bash
cat work/17/428668/.command.sh
```

??? success "Befehlsausgabe"

    ```console
    #!/usr/bin/env bash -e -u -o pipefail
    printf "%s\n" sample1_R1.fastq.gz sample1_R2.fastq.gz | while read f;
    do
        seqtk \
            trimfq \
            -b 5 \
            $f \
            | gzip --no-name > SAMPLE1_PE_$(basename $f)
    done
    ...
    ```

Du solltest `-b 5` im `seqtk trimfq`-Befehl sehen.

Wichtig zu wissen über `ext.args`: Wenn ein Modul bereits einen Standardwert gesetzt hat, **ersetzt** dein Wert diesen vollständig, anstatt ihn zu ergänzen.
Zum Beispiel hat `FASTQC` standardmäßig `ext.args = '--quiet'` in `conf/modules.config` gesetzt:

```groovy title="conf/modules.config" linenums="21" hl_lines="2"
    withName: FASTQC {
        ext.args   = '--quiet'
        publishDir = [
            path: { "${params.outdir}/fastqc/${meta.id}" },
            mode: params.publish_dir_mode,
            pattern: "*.{html,json}",
        ]
    }
```

Wenn du `ext.args = '--kmers 8'` für `FASTQC` setzt, wird das Flag `--quiet` nicht mehr angewendet.
Um beides beizubehalten, setze `ext.args = '--quiet --kmers 8'`.

Du solltest immer die Standardkonfiguration eines Moduls prüfen, bevor du `ext.args` überschreibst.

### Fazit

Du weißt jetzt, wo die Konfigurationsstandards von nf-core Pipelines zu finden sind und wie du Ressourcenzuweisungen und Tool-Argumente mit einer eigenen Konfigurationsdatei überschreibst.

### Wie geht es weiter?

Weiter zu [Teil 3](./03_run_production_pipeline.md), wo du das Gelernte auf eine echte Produktions-Pipeline anwendest.

---

## Zusammenfassung

In diesem Teil hast du gelernt, wie du:

- Hilfe abrufst, Parameter setzt und die Parameter- und Eingabe-Validierung verstehst
- Ressourcenzuweisung und Tool-Argumente über Konfigurationsdateien anpasst
