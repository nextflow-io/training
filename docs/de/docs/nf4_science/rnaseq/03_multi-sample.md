# Teil 3: Implementierung für mehrere Proben mit Paired-End-Daten

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Zuvor hast du eine Pipeline zur Variantenerkennung pro Probe erstellt, die die Daten jeder Probe unabhängig verarbeitet hat.
In diesem Teil des Kurses werden wir unseren einfachen Workflow auf die nächste Stufe heben, indem wir ihn in ein leistungsstarkes Batch-Automatisierungswerkzeug verwandeln, das eine beliebige Anzahl von Proben verarbeiten kann.
Und dabei werden wir ihn auch auf Paired-End-Daten umstellen, die in neueren Studien häufiger vorkommen.

??? info "Wie du von diesem Abschnitt aus beginnst"

    Dieser Abschnitt des Kurses setzt voraus, dass du [Teil 1: Methodenübersicht](./01_method.md), [Teil 2: Einzelproben-Implementierung](./02_single-sample.md) abgeschlossen hast und eine funktionierende `rnaseq.nf`-Pipeline mit ausgefüllten Moduldateien besitzt.

    Falls du Teil 2 nicht abgeschlossen hast oder für diesen Teil neu beginnen möchtest, kannst du die Lösung von Teil 2 als Ausgangspunkt verwenden.
    Führe diese Befehle im Verzeichnis `nf4-science/rnaseq/` aus:

    ```bash
    cp solutions/part2/rnaseq-2.nf rnaseq.nf
    cp solutions/part2/modules/fastqc.nf modules/
    cp solutions/part2/modules/trim_galore.nf modules/
    cp solutions/part2/modules/hisat2_align.nf modules/
    cp solutions/part2/nextflow.config .
    ```

    Dies gibt dir einen vollständigen Workflow zur Verarbeitung einzelner Proben.
    Du kannst testen, ob er erfolgreich läuft:

    ```bash
    nextflow run rnaseq.nf -profile test
    ```

## Aufgabe

In diesem Teil des Kurses werden wir den Workflow erweitern, um Folgendes zu tun:

1. Probeninformationen aus einem CSV-Samplesheet lesen
2. QC, Trimming und Alignment pro Probe für alle Proben parallel ausführen
3. Alle QC-Berichte in einem umfassenden MultiQC-Bericht aggregieren

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/rnaseq/img/rnaseq-wf-03.svg"
</figure>

Dies automatisiert die Schritte aus dem zweiten Abschnitt von [Teil 1: Methodenübersicht](./01_method.md#2-multi-sample-qc-aggregation), wo du diese Befehle manuell in ihren Containern ausgeführt hast.

## Lektionsplan

Wir haben dies in drei Phasen unterteilt:

1. **Den Workflow so anpassen, dass er mehrere Eingabeproben akzeptiert.**
   Dies umfasst den Wechsel von einem einzelnen Dateipfad zu einem CSV-Samplesheet, dessen Parsen mit `splitCsv()` und das Ausführen aller vorhandenen Prozesse auf mehreren Proben.
2. **Umfassende QC-Berichtsgenerierung hinzufügen.**
   Dies führt den Operator `collect()` ein, um Ausgaben über Proben hinweg zu aggregieren, und fügt einen MultiQC-Prozess hinzu, um einen kombinierten Bericht zu erstellen.
3. **Auf Paired-End-RNAseq-Daten umstellen.**
   Dies umfasst die Anpassung von Prozessen für Paired-End-Eingaben (unter Verwendung von Tupeln), das Erstellen von Paired-End-Modulen und das Einrichten eines separaten Testprofils.

Dies implementiert die in [Teil 1: Methodenübersicht](./01_method.md) beschriebene Methode (zweiter Abschnitt zum Multi-Sample-Anwendungsfall) und baut direkt auf dem in Teil 2 erstellten Workflow auf.

!!! tip "Tipp"

     Stelle sicher, dass du im richtigen Arbeitsverzeichnis bist:
     `cd /workspaces/training/nf4-science/rnaseq`

---

## 1. Den Workflow so anpassen, dass er mehrere Eingabeproben akzeptiert

Um auf mehreren Proben zu laufen, müssen wir ändern, wie wir die Eingabe verwalten: Anstatt einen einzelnen Dateipfad bereitzustellen, lesen wir Probeninformationen aus einer CSV-Datei.

Wir stellen eine CSV-Datei zur Verfügung, die Proben-IDs und FASTQ-Dateipfade im Verzeichnis `data/` enthält.

```csv title="data/single-end.csv" linenums="1"
sample_id,fastq_path
ENCSR000COQ1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ1_1.fastq.gz
ENCSR000COQ2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ2_1.fastq.gz
ENCSR000COR1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR1_1.fastq.gz
ENCSR000COR2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR2_1.fastq.gz
ENCSR000CPO1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO1_1.fastq.gz
ENCSR000CPO2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO2_1.fastq.gz
```

Diese CSV-Datei enthält eine Kopfzeile, die die Spalten benennt.

Beachte, dass dies immer noch Single-End-Read-Daten sind.

!!! warning "Warnung"

    Die Dateipfade in der CSV sind absolute Pfade, die zu deiner Umgebung passen müssen.
    Falls du dies nicht in der von uns bereitgestellten Trainingsumgebung ausführst, musst du die Pfade an dein System anpassen.

### 1.1. Die primäre Eingabe im Testprofil auf eine CSV-Datei mit Dateipfaden ändern

Zuerst müssen wir das Testprofil in `nextflow.config` aktualisieren, um den CSV-Dateipfad anstelle des einzelnen FASTQ-Pfads bereitzustellen.

=== "Danach"

    ```groovy title="nextflow.config" linenums="1" hl_lines="5"
    docker.enabled = true

    profiles {
        test {
            params.input = "${projectDir}/data/single-end.csv"
            params.hisat2_index_zip = "${projectDir}/data/genome_index.tar.gz"
        }
    }
    ```

=== "Vorher"

    ```groovy title="nextflow.config" linenums="1"
    docker.enabled = true

    profiles {
        test {
            params.input = "${projectDir}/data/reads/ENCSR000COQ1_1.fastq.gz"
            params.hisat2_index_zip = "${projectDir}/data/genome_index.tar.gz"
        }
    }
    ```

Als Nächstes müssen wir die Kanalerstellung aktualisieren, um aus dieser CSV zu lesen.

### 1.2. Die Channel-Factory aktualisieren, um CSV-Eingaben zu parsen

Wir müssen den Inhalt der Datei in den Kanal laden anstatt nur des Dateipfads selbst.

Wir können dies mit demselben Muster tun, das wir in [Teil 2 von Hello Nextflow](../../hello_nextflow/02_hello_channels.md#42-use-the-splitcsv-operator-to-parse-the-file) verwendet haben: Anwenden des Operators [`splitCsv()`](https://nextflow.io/docs/latest/reference/operator.html#splitcsv) zum Parsen der Datei, dann eine `map`-Operation, um den FASTQ-Dateipfad aus jeder Zeile zu extrahieren.

=== "Danach"

    ```groovy title="rnaseq.nf" linenums="19" hl_lines="4-7"
    workflow {

        main:
        // Eingabekanal aus dem Inhalt einer CSV-Datei erstellen
        read_ch = channel.fromPath(params.input)
            .splitCsv(header: true)
            .map { row -> file(row.fastq_path) }
    ```

=== "Vorher"

    ```groovy title="rnaseq.nf" linenums="19"
    workflow {

        main:
        // Eingabekanal aus einem Dateipfad erstellen
        read_ch = channel.fromPath(params.input)
    ```

Etwas Neues im Vergleich zu dem, was du im Hello Nextflow-Kurs kennengelernt hast, ist, dass diese CSV eine Kopfzeile hat, also fügen wir `#!groovy header: true` zum Aufruf von `splitCsv()` hinzu.
Das erlaubt uns, Spalten nach Namen in der `map`-Operation zu referenzieren: `#!groovy row.fastq_path` extrahiert den Dateipfad aus der Spalte `fastq_path` jeder Zeile.

Die Eingabeverarbeitung ist aktualisiert und der Workflow ist bereit zum Testen.

### 1.3. Den Workflow ausführen

Der Workflow liest jetzt Probeninformationen aus einer CSV-Datei und verarbeitet alle Proben parallel.

```bash
nextflow run rnaseq.nf -profile test
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 24.10.0

    Launching `rnaseq.nf` [golden_curry] DSL2 - revision: 2a5ba5be1e

    executor >  local (18)
    [07/3ff9c5] FASTQC (6)       [100%] 6 of 6 ✔
    [cc/16859f] TRIM_GALORE (6)  [100%] 6 of 6 ✔
    [68/4c27b5] HISAT2_ALIGN (6) [100%] 6 of 6 ✔
    ```

Diesmal wird jeder Schritt 6 Mal ausgeführt, einmal für jede Probe in der CSV-Datei.

Das war alles, was nötig war, um den Workflow auf mehreren Dateien laufen zu lassen.
Nextflow kümmert sich um die gesamte Parallelisierung für uns.

### Fazit

Du weißt jetzt, wie du von einer Einzeldatei-Eingabe zu einer CSV-basierten Multi-Sample-Eingabe wechselst, die Nextflow parallel verarbeitet.

### Wie geht es weiter?

Füge einen QC-Berichtsaggregationsschritt hinzu, der Metriken aus allen Proben kombiniert.

---

## 2. Vor-Verarbeitungs-QC-Metriken in einem einzelnen MultiQC-Bericht aggregieren

All dies erzeugt viele QC-Berichte, und wir möchten nicht durch einzelne Berichte graben müssen.
Dies ist der perfekte Punkt, um einen MultiQC-Berichtsaggregationsschritt einzufügen.

Erinnere dich an den `multiqc`-Befehl aus [Teil 1](01_method.md):

```bash
multiqc . -n <output_name>.html
```

Der Befehl durchsucht das aktuelle Verzeichnis nach erkannten QC-Ausgabedateien und aggregiert sie in einem einzelnen HTML-Bericht.
Die Container-URI war `community.wave.seqera.io/library/pip_multiqc:a3c26f6199d64b7c`.

Wir müssen einen zusätzlichen Parameter einrichten, die Eingaben vorbereiten, den Prozess schreiben, ihn einbinden und die Ausgabeverarbeitung aktualisieren.

### 2.1. Die Eingaben einrichten

Der MultiQC-Prozess benötigt einen Berichtsnamen-Parameter und die gesammelten QC-Ausgaben aus allen vorherigen Schritten gebündelt.

#### 2.1.1. Einen Parameter `report_id` hinzufügen

Füge einen Parameter hinzu, um den Ausgabebericht zu benennen.

=== "Danach"

    ```groovy title="rnaseq.nf" linenums="11" hl_lines="8-9"
    params {
        // Primäre Eingabe
        input: Path

        // Referenzgenom-Archiv
        hisat2_index_zip: Path

        // Bericht-ID
        report_id: String
    }
    ```

=== "Vorher"

    ```groovy title="rnaseq.nf" linenums="11"
    params {
        // Primäre Eingabe
        input: Path

        // Referenzgenom-Archiv
        hisat2_index_zip: Path
    }
    ```

Füge die Standard-Bericht-ID zum Testprofil hinzu:

=== "Danach"

    ```groovy title="nextflow.config" linenums="1" hl_lines="7"
    docker.enabled = true

    profiles {
        test {
            params.input = "${projectDir}/data/single-end.csv"
            params.hisat2_index_zip = "${projectDir}/data/genome_index.tar.gz"
            params.report_id = "all_single-end"
        }
    }
    ```

=== "Vorher"

    ```groovy title="nextflow.config" linenums="1"
    docker.enabled = true

    profiles {
        test {
            params.input = "${projectDir}/data/single-end.csv"
            params.hisat2_index_zip = "${projectDir}/data/genome_index.tar.gz"
        }
    }
    ```

Als Nächstes müssen wir die Eingaben für den MultiQC-Prozess vorbereiten.

#### 2.1.2. QC-Ausgaben aus vorherigen Schritten sammeln und kombinieren

Wir müssen dem Prozess `MULTIQC` alle QC-bezogenen Ausgaben aus den vorherigen Schritten gebündelt geben.

Dafür verwenden wir den Operator `.mix()`, der mehrere Kanäle in einen einzigen aggregiert.
Wir beginnen mit `channel.empty()` und mischen alle Ausgabekanäle ein, die wir kombinieren möchten.
Dies ist sauberer, als `.mix()` direkt an einen der Ausgabekanäle zu ketten, da es alle Eingaben symmetrisch behandelt.

In unserem Workflow sind die zu aggregierenden QC-bezogenen Ausgaben:

- `FASTQC.out.zip`
- `FASTQC.out.html`
- `TRIM_GALORE.out.trimming_reports`
- `TRIM_GALORE.out.fastqc_reports`
- `HISAT2_ALIGN.out.log`

Wir mischen sie in einen einzigen Kanal, dann verwenden wir `.collect()`, um die Berichte über alle Proben hinweg in eine einzige Liste zu aggregieren.

Füge diese Zeilen zum Workflow-Body nach dem Aufruf von `HISAT2_ALIGN` hinzu:

=== "Danach"

    ```groovy title="rnaseq.nf" linenums="38" hl_lines="4-13"
        // Alignment zu einem Referenzgenom
        HISAT2_ALIGN(TRIM_GALORE.out.trimmed_reads, file(params.hisat2_index_zip))

        // Umfassende QC-Berichtsgenerierung
        multiqc_files_ch = channel.empty().mix(
            FASTQC.out.zip,
            FASTQC.out.html,
            TRIM_GALORE.out.trimming_reports,
            TRIM_GALORE.out.fastqc_reports,
            HISAT2_ALIGN.out.log,
        )
        multiqc_files_list = multiqc_files_ch.collect()
    ```

=== "Vorher"

    ```groovy title="rnaseq.nf" linenums="38"
        // Alignment zu einem Referenzgenom
        HISAT2_ALIGN(TRIM_GALORE.out.trimmed_reads, file(params.hisat2_index_zip))
    ```

Die Verwendung von Zwischenvariablen macht jeden Schritt klar: `multiqc_files_ch` enthält alle einzelnen QC-Dateien, die in einen Kanal gemischt wurden, und `multiqc_files_list` ist das gesammelte Bündel, das bereit ist, an MultiQC übergeben zu werden.

### 2.2. Den QC-Aggregationsprozess schreiben und im Workflow aufrufen

Wie zuvor müssen wir die Prozessdefinition ausfüllen, das Modul importieren und den Prozessaufruf hinzufügen.

#### 2.2.1. Das Modul für den QC-Aggregationsprozess ausfüllen

Öffne `modules/multiqc.nf` und untersuche die Gliederung der Prozessdefinition.

Fülle die Prozessdefinition selbst aus, indem du die oben bereitgestellten Informationen verwendest, und überprüfe dann deine Arbeit anhand der Lösung im Tab "Danach" unten.

=== "Vorher"

    ```groovy title="modules/multiqc.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
     * QC-Berichte mit MultiQC aggregieren
     */
    process MULTIQC {

        container

        input:

        output:

        script:
        """

        """
    }
    ```

=== "Danach"

    ```groovy title="modules/multiqc.nf" linenums="1" hl_lines="8 11 12 15 16 20"
    #!/usr/bin/env nextflow

    /*
     * QC-Berichte mit MultiQC aggregieren
     */
    process MULTIQC {

        container "community.wave.seqera.io/library/pip_multiqc:a3c26f6199d64b7c"

        input:
        path '*'
        val output_name

        output:
        path "${output_name}.html", emit: report
        path "${output_name}_data", emit: data

        script:
        """
        multiqc . -n ${output_name}.html
        """
    }
    ```

Dieser Prozess verwendet `#!groovy path '*'` als Eingabequalifizierer für die QC-Dateien.
Der Platzhalter `'*'` teilt Nextflow mit, alle gesammelten Dateien in das Arbeitsverzeichnis zu stagen, ohne spezifische Namen zu erfordern.
Die Eingabe `val output_name` ist ein String, der den Berichtsdateinamen steuert.

Der Befehl `multiqc .` durchsucht das aktuelle Verzeichnis (wo alle gestagten QC-Dateien sind) und generiert den Bericht.

Sobald du dies abgeschlossen hast, ist der Prozess einsatzbereit.

#### 2.2.2. Das Modul einbinden

Füge die Import-Anweisung zu `rnaseq.nf` hinzu:

=== "Danach"

    ```groovy title="rnaseq.nf" linenums="3" hl_lines="5"
    // Modul-INCLUDE-Anweisungen
    include { FASTQC } from './modules/fastqc.nf'
    include { TRIM_GALORE } from './modules/trim_galore.nf'
    include { HISAT2_ALIGN } from './modules/hisat2_align.nf'
    include { MULTIQC } from './modules/multiqc.nf'
    ```

=== "Vorher"

    ```groovy title="rnaseq.nf" linenums="3"
    // Modul-INCLUDE-Anweisungen
    include { FASTQC } from './modules/fastqc.nf'
    include { TRIM_GALORE } from './modules/trim_galore.nf'
    include { HISAT2_ALIGN } from './modules/hisat2_align.nf'
    ```

Füge nun den Prozessaufruf zum Workflow hinzu.

#### 2.2.3. Den Prozessaufruf hinzufügen

Übergib die gesammelten QC-Dateien und die Bericht-ID an den Prozess `MULTIQC`:

=== "Danach"

    ```groovy title="rnaseq.nf" linenums="49" hl_lines="2"
        multiqc_files_list = multiqc_files_ch.collect()
        MULTIQC(multiqc_files_list, params.report_id)
    ```

=== "Vorher"

    ```groovy title="rnaseq.nf" linenums="49"
        multiqc_files_list = multiqc_files_ch.collect()
    ```

Der MultiQC-Prozess ist jetzt in den Workflow eingebunden.

### 2.3. Die Ausgabeverarbeitung aktualisieren

Wir müssen die MultiQC-Ausgaben zur Publish-Deklaration hinzufügen und konfigurieren, wohin sie gehen.

#### 2.3.1. Publish-Ziele für die MultiQC-Ausgaben hinzufügen

Füge die MultiQC-Ausgaben zum Abschnitt `publish:` hinzu:

=== "Danach"

    ```groovy title="rnaseq.nf" linenums="52" hl_lines="9-10"
        publish:
        fastqc_zip = FASTQC.out.zip
        fastqc_html = FASTQC.out.html
        trimmed_reads = TRIM_GALORE.out.trimmed_reads
        trimming_reports = TRIM_GALORE.out.trimming_reports
        trimming_fastqc = TRIM_GALORE.out.fastqc_reports
        bam = HISAT2_ALIGN.out.bam
        align_log = HISAT2_ALIGN.out.log
        multiqc_report = MULTIQC.out.report
        multiqc_data = MULTIQC.out.data
    }
    ```

=== "Vorher"

    ```groovy title="rnaseq.nf" linenums="52"
        publish:
        fastqc_zip = FASTQC.out.zip
        fastqc_html = FASTQC.out.html
        trimmed_reads = TRIM_GALORE.out.trimmed_reads
        trimming_reports = TRIM_GALORE.out.trimming_reports
        trimming_fastqc = TRIM_GALORE.out.fastqc_reports
        bam = HISAT2_ALIGN.out.bam
        align_log = HISAT2_ALIGN.out.log
    }
    ```

Als Nächstes müssen wir Nextflow mitteilen, wohin diese Ausgaben sollen.

#### 2.3.2. Die neuen Ausgabeziele konfigurieren

Füge Einträge für die MultiQC-Ziele im Block `output {}` hinzu und veröffentliche sie in einem Unterverzeichnis `multiqc/`:

=== "Danach"

    ```groovy title="rnaseq.nf" linenums="82" hl_lines="7-12"
        align_log {
            path 'align'
        }
        multiqc_report {
            path 'multiqc'
        }
        multiqc_data {
            path 'multiqc'
        }
    }
    ```

=== "Vorher"

    ```groovy title="rnaseq.nf" linenums="82"
        align_log {
            path 'align'
        }
    }
    ```

Die Ausgabekonfiguration ist vollständig.

### 2.4. Den Workflow ausführen

Wir verwenden `-resume`, damit die vorherigen Verarbeitungsschritte gecacht werden und nur der neue MultiQC-Schritt läuft.

```bash
nextflow run rnaseq.nf -profile test -resume
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 24.10.0

    Launching `rnaseq.nf` [modest_pare] DSL2 - revision: fc724d3b49

    executor >  local (1)
    [07/3ff9c5] FASTQC (6)       [100%] 6 of 6, cached: 6 ✔
    [2c/8d8e1e] TRIM_GALORE (5)  [100%] 6 of 6, cached: 6 ✔
    [a4/7f9c44] HISAT2_ALIGN (6) [100%] 6 of 6, cached: 6 ✔
    [56/e1f102] MULTIQC          [100%] 1 of 1 ✔
    ```

Ein einzelner Aufruf von MULTIQC wurde nach den gecachten Prozessaufrufen hinzugefügt.

Du findest die MultiQC-Ausgaben im Ergebnisverzeichnis.

```bash
tree -L 2 results/multiqc
```

```console title="Ausgabe"
results/multiqc
├── all_single-end_data
│   ├── cutadapt_filtered_reads_plot.txt
│   ├── cutadapt_trimmed_sequences_plot_3_Counts.txt
│   ├── cutadapt_trimmed_sequences_plot_3_Obs_Exp.txt
│   ├── fastqc_adapter_content_plot.txt
│   ├── fastqc_overrepresented_sequences_plot.txt
│   ├── fastqc_per_base_n_content_plot.txt
│   ├── fastqc_per_base_sequence_quality_plot.txt
│   ├── fastqc_per_sequence_gc_content_plot_Counts.txt
│   ├── fastqc_per_sequence_gc_content_plot_Percentages.txt
│   ├── fastqc_per_sequence_quality_scores_plot.txt
│   ├── fastqc_sequence_counts_plot.txt
│   ├── fastqc_sequence_duplication_levels_plot.txt
│   ├── fastqc_sequence_length_distribution_plot.txt
│   ├── fastqc-status-check-heatmap.txt
│   ├── fastqc_top_overrepresented_sequences_table.txt
│   ├── hisat2_se_plot.txt
│   ├── multiqc_citations.txt
│   ├── multiqc_cutadapt.txt
│   ├── multiqc_data.json
│   ├── multiqc_fastqc.txt
│   ├── multiqc_general_stats.txt
│   ├── multiqc_hisat2.txt
│   ├── multiqc.log
│   ├── multiqc_software_versions.txt
│   └── multiqc_sources.txt
└── all_single-end.html
```

Die letzte Datei `all_single-end.html` ist der vollständige aggregierte Bericht, praktisch in eine einfach zu durchsuchende HTML-Datei verpackt.

### Fazit

Du weißt jetzt, wie du Ausgaben aus mehreren Kanälen sammelst, sie mit `.mix()` und `.collect()` bündelst und an einen Aggregationsprozess übergibst.

### Wie geht es weiter?

Passe den Workflow an, um Paired-End-RNAseq-Daten zu verarbeiten.

---

## 3. Verarbeitung von Paired-End-RNAseq-Daten ermöglichen

Momentan kann unser Workflow nur Single-End-RNAseq-Daten verarbeiten.
Es wird zunehmend üblich, Paired-End-RNAseq-Daten zu sehen, also möchten wir in der Lage sein, diese zu verarbeiten.

Den Workflow völlig unabhängig vom Datentyp zu machen, würde etwas fortgeschrittenere Nextflow-Sprachfeatures erfordern, also werden wir das hier nicht tun, aber wir können eine Paired-End-Verarbeitungsversion erstellen, um zu demonstrieren, was angepasst werden muss.

### 3.1. Den Workflow kopieren und die Eingaben aktualisieren

Wir beginnen damit, die Single-End-Workflow-Datei zu kopieren und sie für Paired-End-Daten zu aktualisieren.

#### 3.1.1. Die Workflow-Datei kopieren

Erstelle eine Kopie der Workflow-Datei, um sie als Ausgangspunkt für die Paired-End-Version zu verwenden.

```bash
cp rnaseq.nf rnaseq_pe.nf
```

Aktualisiere nun die Parameter und die Eingabeverarbeitung in der neuen Datei.

#### 3.1.2. Ein Paired-End-Testprofil hinzufügen

Wir stellen eine zweite CSV-Datei zur Verfügung, die Proben-IDs und gepaarte FASTQ-Dateipfade im Verzeichnis `data/` enthält.

```csv title="data/paired-end.csv" linenums="1"
sample_id,fastq_1,fastq_2
ENCSR000COQ1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ1_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ1_2.fastq.gz
ENCSR000COQ2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ2_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ2_2.fastq.gz
ENCSR000COR1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR1_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR1_2.fastq.gz
ENCSR000COR2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR2_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR2_2.fastq.gz
ENCSR000CPO1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO1_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO1_2.fastq.gz
ENCSR000CPO2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO2_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO2_2.fastq.gz
```

Füge ein Profil `test_pe` zu `nextflow.config` hinzu, das auf diese Datei zeigt und eine Paired-End-Bericht-ID verwendet.

=== "Danach"

    ```groovy title="nextflow.config" linenums="1" hl_lines="9-13"
    docker.enabled = true

    profiles {
        test {
            params.input = "${projectDir}/data/single-end.csv"
            params.hisat2_index_zip = "${projectDir}/data/genome_index.tar.gz"
            params.report_id = "all_single-end"
        }
        test_pe {
            params.input = "${projectDir}/data/paired-end.csv"
            params.hisat2_index_zip = "${projectDir}/data/genome_index.tar.gz"
            params.report_id = "all_paired-end"
        }
    }
    ```

=== "Vorher"

    ```groovy title="nextflow.config" linenums="1"
    docker.enabled = true

    profiles {
        test {
            params.input = "${projectDir}/data/single-end.csv"
            params.hisat2_index_zip = "${projectDir}/data/genome_index.tar.gz"
            params.report_id = "all_single-end"
        }
    }
    ```

Das Testprofil für Paired-End-Daten ist bereit.

#### 3.1.3. Die Channel-Factory aktualisieren

Der Operator `.map()` muss beide FASTQ-Dateipfade greifen und sie als Liste zurückgeben.

=== "Danach"

    ```groovy title="rnaseq_pe.nf" linenums="25" hl_lines="4"
        // Eingabekanal aus dem Inhalt einer CSV-Datei erstellen
        read_ch = channel.fromPath(params.input)
            .splitCsv(header: true)
            .map { row -> [file(row.fastq_1), file(row.fastq_2)] }
    ```

=== "Vorher"

    ```groovy title="rnaseq_pe.nf" linenums="25"
        // Eingabekanal aus dem Inhalt einer CSV-Datei erstellen
        read_ch = channel.fromPath(params.input)
            .splitCsv(header: true)
            .map { row -> file(row.fastq_path) }
    ```

Die Eingabeverarbeitung ist für Paired-End-Daten konfiguriert.

### 3.2. Das FASTQC-Modul für Paired-End-Daten anpassen

Kopiere das Modul, um eine Paired-End-Version zu erstellen:

```bash
cp modules/fastqc.nf modules/fastqc_pe.nf
```

Die Eingabe des FASTQC-Prozesses muss sich nicht ändern — wenn Nextflow eine Liste von zwei Dateien erhält, stagt es beide und `reads` expandiert zu beiden Dateinamen.
Die einzige notwendige Änderung ist im Ausgabeblock: Da wir jetzt zwei FastQC-Berichte pro Probe erhalten, wechseln wir von `simpleName`-basierten Mustern zu Platzhaltern.

=== "Danach"

    ```groovy title="modules/fastqc_pe.nf" linenums="10" hl_lines="2 3"
        output:
        path "*_fastqc.zip", emit: zip
        path "*_fastqc.html", emit: html
    ```

=== "Vorher"

    ```groovy title="modules/fastqc_pe.nf" linenums="10"
        output:
        path "${reads.simpleName}_fastqc.zip", emit: zip
        path "${reads.simpleName}_fastqc.html", emit: html
    ```

Dies verallgemeinert den Prozess auf eine Weise, die ihn in die Lage versetzt, entweder Single-End- oder Paired-End-Daten zu verarbeiten.

Aktualisiere den Import in `rnaseq_pe.nf`, um die Paired-End-Version zu verwenden:

=== "Danach"

    ```groovy title="rnaseq_pe.nf" linenums="4" hl_lines="1"
    include { FASTQC } from './modules/fastqc_pe.nf'
    ```

=== "Vorher"

    ```groovy title="rnaseq_pe.nf" linenums="4"
    include { FASTQC } from './modules/fastqc.nf'
    ```

Das FASTQC-Modul und sein Import sind für Paired-End-Daten aktualisiert.

### 3.3. Das TRIM_GALORE-Modul für Paired-End-Daten anpassen

Kopiere das Modul, um eine Paired-End-Version zu erstellen:

```bash
cp modules/trim_galore.nf modules/trim_galore_pe.nf
```

Dieses Modul benötigt umfangreichere Änderungen:

- Die Eingabe ändert sich von einem einzelnen Pfad zu einem Tupel aus zwei Pfaden
- Der Befehl fügt das Flag `--paired` hinzu und nimmt beide Read-Dateien
- Die Ausgabe ändert sich, um die Paired-End-Benennungskonventionen von Trim Galore widerzuspiegeln und separate FastQC-Berichte für jede Read-Datei zu erzeugen

=== "Danach"

    ```groovy title="modules/trim_galore_pe.nf" linenums="8" hl_lines="2 5 7 8 12"
        input:
        tuple path(read1), path(read2)

        output:
        tuple path("*_val_1.fq.gz"), path("*_val_2.fq.gz"), emit: trimmed_reads
        path "*_trimming_report.txt", emit: trimming_reports
        path "*_val_1_fastqc.{zip,html}", emit: fastqc_reports_1
        path "*_val_2_fastqc.{zip,html}", emit: fastqc_reports_2

        script:
        """
        trim_galore --fastqc --paired ${read1} ${read2}
        """
    ```

=== "Vorher"

    ```groovy title="modules/trim_galore_pe.nf" linenums="8"
        input:
        path reads

        output:
        path "${reads.simpleName}_trimmed.fq.gz", emit: trimmed_reads
        path "${reads}_trimming_report.txt", emit: trimming_reports
        path "${reads.simpleName}_trimmed_fastqc.{zip,html}", emit: fastqc_reports

        script:
        """
        trim_galore --fastqc ${reads}
        """
    ```

Aktualisiere den Import in `rnaseq_pe.nf`:

=== "Danach"

    ```groovy title="rnaseq_pe.nf" linenums="5" hl_lines="1"
    include { TRIM_GALORE } from './modules/trim_galore_pe.nf'
    ```

=== "Vorher"

    ```groovy title="rnaseq_pe.nf" linenums="5"
    include { TRIM_GALORE } from './modules/trim_galore.nf'
    ```

Das TRIM_GALORE-Modul und sein Import sind für Paired-End-Daten aktualisiert.

### 3.4. Das HISAT2_ALIGN-Modul für Paired-End-Daten anpassen

Kopiere das Modul, um eine Paired-End-Version zu erstellen:

```bash
cp modules/hisat2_align.nf modules/hisat2_align_pe.nf
```

Dieses Modul benötigt ähnliche Änderungen:

- Die Eingabe ändert sich von einem einzelnen Pfad zu einem Tupel aus zwei Pfaden
- Der HISAT2-Befehl ändert sich von `-U` (unpaired) zu `-1` und `-2` (paired) Read-Argumenten
- Alle Verwendungen von `reads.simpleName` ändern sich zu `read1.simpleName`, da wir jetzt auf ein spezifisches Mitglied des Paares verweisen

=== "Danach"

    ```groovy title="modules/hisat2_align_pe.nf" linenums="8" hl_lines="2 6 7 12 13 14"
        input:
        tuple path(read1), path(read2)
        path index_zip

        output:
        path "${read1.simpleName}.bam", emit: bam
        path "${read1.simpleName}.hisat2.log", emit: log

        script:
        """
        tar -xzvf ${index_zip}
        hisat2 -x ${index_zip.simpleName} -1 ${read1} -2 ${read2} \
            --new-summary --summary-file ${read1.simpleName}.hisat2.log | \
            samtools view -bS -o ${read1.simpleName}.bam
        """
    ```

=== "Vorher"

    ```groovy title="modules/hisat2_align_pe.nf" linenums="8"
        input:
        path reads
        path index_zip

        output:
        path "${reads.simpleName}.bam", emit: bam
        path "${reads.simpleName}.hisat2.log", emit: log

        script:
        """
        tar -xzvf ${index_zip}
        hisat2 -x ${index_zip.simpleName} -U ${reads} \
            --new-summary --summary-file ${reads.simpleName}.hisat2.log | \
            samtools view -bS -o ${reads.simpleName}.bam
        """
    ```

Aktualisiere den Import in `rnaseq_pe.nf`:

=== "Danach"

    ```groovy title="rnaseq_pe.nf" linenums="6" hl_lines="1"
    include { HISAT2_ALIGN } from './modules/hisat2_align_pe.nf'
    ```

=== "Vorher"

    ```groovy title="rnaseq_pe.nf" linenums="6"
    include { HISAT2_ALIGN } from './modules/hisat2_align.nf'
    ```

Das HISAT2_ALIGN-Modul und sein Import sind für Paired-End-Daten aktualisiert.

### 3.5. Die MultiQC-Aggregation für Paired-End-Ausgaben aktualisieren

Der Paired-End-Prozess `TRIM_GALORE` erzeugt jetzt zwei separate FastQC-Berichtskanäle (`fastqc_reports_1` und `fastqc_reports_2`) anstelle von einem.
Aktualisiere den Block `.mix()` in `rnaseq_pe.nf`, um beide einzuschließen:

=== "Danach"

    ```groovy title="rnaseq_pe.nf" linenums="40" hl_lines="5 6"
        multiqc_files_ch = channel.empty().mix(
            FASTQC.out.zip,
            FASTQC.out.html,
            TRIM_GALORE.out.trimming_reports,
            TRIM_GALORE.out.fastqc_reports_1,
            TRIM_GALORE.out.fastqc_reports_2,
            HISAT2_ALIGN.out.log,
        )
    ```

=== "Vorher"

    ```groovy title="rnaseq_pe.nf" linenums="40"
        multiqc_files_ch = channel.empty().mix(
            FASTQC.out.zip,
            FASTQC.out.html,
            TRIM_GALORE.out.trimming_reports,
            TRIM_GALORE.out.fastqc_reports,
            HISAT2_ALIGN.out.log,
        )
    ```

Die MultiQC-Aggregation enthält jetzt beide Sätze von Paired-End-FastQC-Berichten.

### 3.6. Die Ausgabeverarbeitung für Paired-End-Ausgaben aktualisieren

Der Abschnitt `publish:` und der Block `output {}` müssen ebenfalls die zwei separaten FastQC-Berichtskanäle aus dem Paired-End-Prozess `TRIM_GALORE` widerspiegeln.

Aktualisiere den Abschnitt `publish:` in `rnaseq_pe.nf`:

=== "Danach"

    ```groovy title="rnaseq_pe.nf" linenums="52" hl_lines="6-7"
        publish:
        fastqc_zip = FASTQC.out.zip
        fastqc_html = FASTQC.out.html
        trimmed_reads = TRIM_GALORE.out.trimmed_reads
        trimming_reports = TRIM_GALORE.out.trimming_reports
        trimming_fastqc_1 = TRIM_GALORE.out.fastqc_reports_1
        trimming_fastqc_2 = TRIM_GALORE.out.fastqc_reports_2
        bam = HISAT2_ALIGN.out.bam
        align_log = HISAT2_ALIGN.out.log
        multiqc_report = MULTIQC.out.report
        multiqc_data = MULTIQC.out.data
    }
    ```

=== "Vorher"

    ```groovy title="rnaseq_pe.nf" linenums="52"
        publish:
        fastqc_zip = FASTQC.out.zip
        fastqc_html = FASTQC.out.html
        trimmed_reads = TRIM_GALORE.out.trimmed_reads
        trimming_reports = TRIM_GALORE.out.trimming_reports
        trimming_fastqc = TRIM_GALORE.out.fastqc_reports
        bam = HISAT2_ALIGN.out.bam
        align_log = HISAT2_ALIGN.out.log
        multiqc_report = MULTIQC.out.report
        multiqc_data = MULTIQC.out.data
    }
    ```

Aktualisiere die entsprechenden Einträge im Block `output {}`:

=== "Danach"

    ```groovy title="rnaseq_pe.nf" linenums="77" hl_lines="4-9"
        trimming_reports {
            path 'trimming'
        }
        trimming_fastqc_1 {
            path 'trimming'
        }
        trimming_fastqc_2 {
            path 'trimming'
        }
        bam {
            path 'align'
        }
    ```

=== "Vorher"

    ```groovy title="rnaseq_pe.nf" linenums="77"
        trimming_reports {
            path 'trimming'
        }
        trimming_fastqc {
            path 'trimming'
        }
        bam {
            path 'align'
        }
    ```

Der Paired-End-Workflow ist jetzt vollständig aktualisiert und bereit zum Ausführen.

### 3.7. Den Workflow ausführen

Wir verwenden nicht `-resume`, da dies nicht cachen würde, und es gibt doppelt so viele Daten zu verarbeiten wie zuvor, aber es sollte trotzdem in weniger als einer Minute abgeschlossen sein.

```bash
nextflow run rnaseq_pe.nf -profile test_pe
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 24.10.0

    Launching `rnaseq_pe.nf` [reverent_kare] DSL2 - revision: 9c376cc219

    executor >  local (19)
    [c5/cbde15] FASTQC (5)       [100%] 6 of 6 ✔
    [e4/fa2784] TRIM_GALORE (5)  [100%] 6 of 6 ✔
    [3a/e23049] HISAT2_ALIGN (5) [100%] 6 of 6 ✔
    [e6/a3ccd9] MULTIQC          [100%] 1 of 1 ✔
    ```

Jetzt haben wir zwei leicht divergierende Versionen unseres Workflows, eine für Single-End-Read-Daten und eine für Paired-End-Daten.
Der nächste logische Schritt wäre, den Workflow so zu gestalten, dass er beide Datentypen spontan akzeptiert, was außerhalb des Umfangs dieses Kurses liegt, aber wir könnten das in einem Folgekurs angehen.

---

### Fazit

Du weißt jetzt, wie du einen Einzelproben-Workflow anpasst, um die Verarbeitung mehrerer Proben zu parallelisieren, einen umfassenden QC-Bericht zu generieren und den Workflow für Paired-End-Read-Daten anzupassen.

### Wie geht es weiter?

Klopf dir selbst auf die Schulter! Du hast den Nextflow for RNAseq-Kurs abgeschlossen.

Gehe weiter zur abschließenden [Kurszusammenfassung](./next_steps.md), um zu überprüfen, was du gelernt hast, und herauszufinden, was als Nächstes kommt.
