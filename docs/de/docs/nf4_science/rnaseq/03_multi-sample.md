# Teil 3: Implementierung für mehrere Proben mit Paired-End-Daten

In diesem letzten Teil des Kurses werden wir unseren einfachen Workflow auf die nächste Stufe heben, indem wir ihn in ein leistungsstarkes Batch-Automatisierungstool verwandeln, das beliebig viele Proben verarbeiten kann.
Und während wir dabei sind, werden wir ihn auch so anpassen, dass er Paired-End-Daten erwartet, die in neueren Studien häufiger vorkommen.

Wir machen das in drei Schritten:

1. Den Workflow so anpassen, dass er mehrere Eingabeproben akzeptiert und die Ausführung parallelisiert
2. Umfassende QC-Berichtserstellung hinzufügen
3. Auf Paired-End-RNAseq-Daten umstellen

---

## 1. Den Workflow so anpassen, dass er mehrere Eingabeproben akzeptiert und die Ausführung parallelisiert

Wir müssen ändern, wie wir die Eingabe verwalten.

### 1.1. Die primäre Eingabe in eine CSV-Datei mit Dateipfaden ändern statt einer einzelnen Datei

Wir stellen eine CSV-Datei bereit, die Proben-IDs und FASTQ-Dateipfade im Verzeichnis `data/` enthält.
Diese CSV-Datei enthält eine Kopfzeile.
Beachte, dass die FASTQ-Dateipfade absolute Pfade sind.

```csv title="data/single-end.csv" linenums="1"
sample_id,fastq_path
ENCSR000COQ1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ1_1.fastq.gz
ENCSR000COQ2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ2_1.fastq.gz
ENCSR000COR1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR1_1.fastq.gz
ENCSR000COR2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR2_1.fastq.gz
ENCSR000CPO1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO1_1.fastq.gz
ENCSR000CPO2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO2_1.fastq.gz
```

Lass uns den primären Eingabeparameter in `input_csv` umbenennen und den Standardwert auf den Pfad zur Datei `single-end.csv` ändern.

```groovy title="rnaseq.nf" linenums="13"
params {
    // Primary input
    input_csv: Path = "data/single-end.csv"

    // Reference genome archive
    hisat2_index_zip: Path = "data/genome_index.tar.gz"
}
```

### 1.2. Die Eingabekanal-Factory aktualisieren, um eine CSV als Eingabe zu verarbeiten

Wir wollen den Inhalt der Datei in den Kanal laden, anstatt nur den Dateipfad selbst. Dazu verwenden wir den Operator `.splitCsv()`, um das CSV-Format zu parsen, und dann den Operator `.map()`, um die spezifische Information zu extrahieren, die wir wollen (den FASTQ-Dateipfad).

```groovy title="rnaseq.nf" linenums="16"
    // Create input channel from the contents of a CSV file
    read_ch = channel.fromPath(params.input_csv)
        .splitCsv(header:true)
        .map { row -> file(row.fastq_path) }
```

### 1.3. Den Workflow ausführen, um zu testen, dass er funktioniert

```bash
nextflow run rnaseq.nf
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

Diesmal sehen wir, dass jeder Schritt 6-mal ausgeführt wird, für jede der 6 Datendateien, die wir bereitgestellt haben.

Das war alles, was nötig war, um den Workflow auf mehreren Dateien laufen zu lassen!
Nextflow kümmert sich um die gesamte Parallelisierung für uns.

---

## 2. Pre-Processing-QC-Metriken in einem einzigen MultiQC-Bericht aggregieren

All das erzeugt viele QC-Berichte, und wir wollen nicht durch einzelne Berichte graben müssen.
Das ist der perfekte Punkt, um einen MultiQC-Berichtsaggregationsschritt einzufügen!

### 2.1. Ein Modul für den QC-Aggregationsprozess erstellen

Lass uns eine Moduldatei namens `modules/multiqc.nf` erstellen, um den Prozess `MULTIQC` unterzubringen:

```bash
touch modules/multiqc.nf
```

Öffne die Datei im Code-Editor und kopiere den folgenden Code hinein:

```groovy title="modules/multiqc.nf" linenums="1"
#!/usr/bin/env nextflow

process MULTIQC {

    container "community.wave.seqera.io/library/pip_multiqc:a3c26f6199d64b7c"
    publishDir "results/multiqc", mode: 'symlink'

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

### 2.2. Das Modul in die Workflow-Datei importieren

Füge die Anweisung `include { MULTIQC } from './modules/multiqc.nf'` zur Datei `rnaseq.nf` hinzu:

```groovy title="rnaseq.nf" linenums="3"
// Module INCLUDE statements
include { FASTQC } from './modules/fastqc.nf'
include { TRIM_GALORE } from './modules/trim_galore.nf'
include { HISAT2_ALIGN } from './modules/hisat2_align.nf'
include { MULTIQC } from './modules/multiqc.nf'
```

### 2.3. Einen Parameter `report_id` hinzufügen und ihm einen sinnvollen Standardwert geben

```groovy title="rnaseq.nf" linenums="9"
params {
    // Primary input
    input_csv: Path = "data/single-end.csv"

    // Reference genome archive
    hisat2_index_zip: Path = "data/genome_index.tar.gz"

    // Report ID
    report_id: String = "all_single-end"
}
```

### 2.4. Den Prozess mit den Ausgaben der vorherigen Schritte aufrufen

Wir müssen dem Prozess `MULTIQC` alle QC-bezogenen Ausgaben aus den vorherigen Schritten geben.

Dafür verwenden wir den Operator `.mix()`, der mehrere Kanäle zu einem einzigen aggregiert.

Wenn wir vier Prozesse namens A, B, C und D mit jeweils einem einfachen Kanal `.out` hätten, würde die Syntax so aussehen: `A.out.mix( B.out, C.out, D.out )`. Wie du siehst, wendest du ihn auf den ersten der Kanäle an, die du kombinieren möchtest (egal welchen), und fügst einfach alle anderen, durch Kommas getrennt, in den Klammern hinzu.

Im Fall unseres Workflows haben wir folgende Ausgaben zu aggregieren:

- `FASTQC.out.zip`
- `FASTQC.out.html`
- `TRIM_GALORE.out.trimming_reports`
- `TRIM_GALORE.out.fastqc_reports`
- `HISAT2_ALIGN.out.log`

Das Syntaxbeispiel wird also zu:

```groovy title="Applying .mix() in the MULTIQC call"
        FASTQC.out.zip.mix(
        FASTQC.out.html,
        TRIM_GALORE.out.trimming_reports,
        TRIM_GALORE.out.fastqc_reports,
        HISAT2_ALIGN.out.log
        )
```

Das sammelt QC-Berichte pro Probe.
Aber da wir sie über alle Proben hinweg aggregieren wollen, müssen wir den Operator `collect()` hinzufügen, um die Berichte für alle Proben in einen einzigen Aufruf von `MULTIQC` zu ziehen.
Und wir müssen ihm auch den Parameter `report_id` geben.

Das ergibt Folgendes:

```groovy title="The completed MULTIQC call" linenums="33"
    // Comprehensive QC report generation
    MULTIQC(
        FASTQC.out.zip.mix(
        FASTQC.out.html,
        TRIM_GALORE.out.trimming_reports,
        TRIM_GALORE.out.fastqc_reports,
        HISAT2_ALIGN.out.log
        ).collect(),
        params.report_id
    )
```

Im Kontext des vollständigen Workflow-Blocks sieht es so aus:

```groovy title="rnaseq.nf" linenums="18"
workflow {
    // Create input channel from the contents of a CSV file
    read_ch = channel.fromPath(params.input_csv)
        .splitCsv(header:true)
        .map { row -> file(row.fastq_path) }

    /// Initial quality control
    FASTQC(read_ch)

    // Adapter trimming and post-trimming QC
    TRIM_GALORE(read_ch)

    // Alignment to a reference genome
    HISAT2_ALIGN(TRIM_GALORE.out.trimmed_reads, file (params.hisat2_index_zip))

    // Comprehensive QC report generation
    MULTIQC(
        FASTQC.out.zip.mix(
        FASTQC.out.html,
        TRIM_GALORE.out.trimming_reports,
        TRIM_GALORE.out.fastqc_reports,
        HISAT2_ALIGN.out.log
        ).collect(),
        params.report_id
    )
}
```

### 2.5. Den Workflow ausführen, um zu testen, dass er funktioniert

```bash
nextflow run rnaseq.nf -resume
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

Diesmal sehen wir einen einzigen Aufruf von MULTIQC, der nach den gecachten Prozessaufrufen hinzugefügt wurde:

Du findest die Ausgaben unter `results/trimming`, wie im Prozess `TRIM_GALORE` durch die Direktive `publishDir` angegeben.

```bash
tree -L 2 results/multiqc
```

```console title="Output"
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

Die letzte Datei `all_single-end.html` ist der vollständige aggregierte Bericht, praktisch in einer einzigen, leicht zu durchsuchenden HTML-Datei verpackt.

---

## 3. Verarbeitung von Paired-End-RNAseq-Daten ermöglichen

Momentan kann unser Workflow nur Single-End-RNAseq-Daten verarbeiten.
Es ist zunehmend üblich, Paired-End-RNAseq-Daten zu sehen, also wollen wir das auch verarbeiten können.

Den Workflow vollständig unabhängig vom Datentyp zu machen, würde etwas fortgeschrittenere Nextflow-Sprachfunktionen erfordern, also werden wir das hier nicht tun, aber wir können eine Paired-End-Verarbeitungsversion erstellen, um zu zeigen, was angepasst werden muss.

### 3.1. Eine Kopie des Workflows namens `rnaseq_pe.nf` erstellen

```bash
cp rnaseq.nf rnaseq_pe.nf
```

### 3.2. Den Standard-`input_csv` so ändern, dass er auf die Paired-End-Daten zeigt

Wir stellen eine zweite CSV-Datei bereit, die Proben-IDs und gepaarte FASTQ-Dateipfade im Verzeichnis `data/` enthält

```csv title="data/paired-end.csv" linenums="1"
sample_id,fastq_1,fastq_2
ENCSR000COQ1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ1_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ1_2.fastq.gz
ENCSR000COQ2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ2_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ2_2.fastq.gz
ENCSR000COR1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR1_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR1_2.fastq.gz
ENCSR000COR2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR2_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR2_2.fastq.gz
ENCSR000CPO1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO1_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO1_2.fastq.gz
ENCSR000CPO2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO2_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO2_2.fastq.gz
```

Lass uns den Standard für `input_csv` auf den Pfad zur Datei `paired-end.csv` ändern.

```groovy title="rnaseq_pe.nf" linenums="15"
params {
    // Primary input
    input_csv: Path = "data/paired-end.csv"

    // Reference genome archive
    hisat2_index_zip: Path = "data/genome_index.tar.gz"

    // Report ID
    report_id: String = "all_single-end"
}
```

### 3.3. Die Kanal-Factory aktualisieren

Wir müssen dem Operator `.map()` sagen, dass er jetzt beide FASTQ-Dateipfade greifen soll.

Also wird `row -> file(row.fastq_path)` zu `row -> [file(row.fastq_1), file(row.fastq_2)]`

```groovy title="rnaseq_pe.nf" linenums="19"
    // Create input channel from the contents of a CSV file
    read_ch = channel.fromPath(params.input_csv)
        .splitCsv(header:true)
        .map { row -> [file(row.fastq_1), file(row.fastq_2)] }
```

### 3.4. Eine Paired-End-Version des FASTQC-Prozesses erstellen

Lass uns eine Kopie des Moduls erstellen, damit wir beide Versionen zur Hand haben.

```bash
cp modules/fastqc.nf modules/fastqc_pe.nf
```

Öffne die neue Moduldatei `fastqc_pe.nf` im Code-Editor und nimm folgende Codeänderungen vor:

- Ändere `fastqc $reads` zu `fastqc ${reads}` im Block `script` (Zeile 17), damit die Eingabe `reads` entpackt wird, da es jetzt ein Tupel aus zwei Pfaden statt eines einzelnen Pfads ist.
- Ersetze `${reads.simpleName}` durch einen Platzhalter (`*`), um zu vermeiden, dass die Ausgabedateien einzeln behandelt werden müssen.

```groovy title="modules/fastqc_pe.nf" linenums="8"
    input:
    path reads

    output:
    path "*_fastqc.zip", emit: zip
    path "*_fastqc.html", emit: html

    script:
    """
    fastqc ${reads}
    """
```

Technisch gesehen verallgemeinert dies den Prozess `FASTQC` so, dass er entweder Single-End- oder Paired-End-RNAseq-Daten verarbeiten kann.

Aktualisiere abschließend die Modul-Import-Anweisung, um die Paired-End-Version des Moduls zu verwenden.

```groovy title="rnaseq_pe.nf" linenums="4"
include { FASTQC } from './modules/fastqc_pe.nf'
```

### 3.5. Eine Paired-End-Version des TRIM_GALORE-Prozesses erstellen

Erstelle eine Kopie des Moduls, damit wir beide Versionen zur Hand haben.

```bash
cp modules/trim_galore.nf modules/trim_galore_pe.nf
```

Öffne die neue Moduldatei `trim_galore_pe.nf` im Code-Editor und nimm folgende Codeänderungen vor:

- Ändere die Eingabedeklaration von `path reads` zu `tuple path(read1), path(read2)`
- Aktualisiere den Befehl im Block `script`, ersetze `$reads` durch `--paired ${read1} ${read2}`
- Aktualisiere die Ausgabedeklarationen, um die hinzugefügten Dateien und unterschiedlichen Namenskonventionen widerzuspiegeln, verwende Platzhalter, um zu vermeiden, dass alles aufgelistet werden muss.

```groovy title="modules/trim_galore_pe.nf" linenums="8"
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

Aktualisiere abschließend die Modul-Import-Anweisung, um die Paired-End-Version des Moduls zu verwenden.

```groovy title="rnaseq_pe.nf" linenums="5"
include { TRIM_GALORE } from './modules/trim_galore_pe.nf'
```

### 3.6. Den Aufruf des MULTIQC-Prozesses aktualisieren, um zwei Berichte von TRIM_GALORE zu erwarten

Der Prozess `TRIM_GALORE` erzeugt jetzt einen zusätzlichen Ausgabekanal, also müssen wir diesen an MultiQC weiterleiten.

Ersetze `TRIM_GALORE.out.fastqc_reports,` durch `TRIM_GALORE.out.fastqc_reports_1,` plus `TRIM_GALORE.out.fastqc_reports_2,`:

```groovy title="rnaseq_pe.nf" linenums="33"
    // Comprehensive QC report generation
    MULTIQC(
        FASTQC.out.zip.mix(
        FASTQC.out.html,
        TRIM_GALORE.out.trimming_reports,
        TRIM_GALORE.out.fastqc_reports_1,
        TRIM_GALORE.out.fastqc_reports_2,
        HISAT2_ALIGN.out.log
        ).collect(),
        params.report_id
    )
```

Während wir bei MultiQC sind, lass uns auch den Standardwert des Parameters `report_id` von `"all_single-end"` zu `"all_paired-end"` aktualisieren.

```groovy title="rnaseq_pe.nf" linenums="9"
params {
    // Primary input
    input_csv: Path = "data/paired-end.csv"

    // Reference genome archive
    hisat2_index_zip: Path = "data/genome_index.tar.gz"

    // Report ID
    report_id: String = "all_paired-end"
}
```

### 3.7. Eine Paired-End-Version des HISAT2_ALIGN-Prozesses erstellen

Erstelle eine Kopie des Moduls, damit wir beide Versionen zur Hand haben.

```bash
cp modules/hisat2_align.nf modules/hisat2_align_pe.nf
```

Öffne die neue Moduldatei `hisat2_align_pe.nf` im Code-Editor und nimm folgende Codeänderungen vor:

- Ändere die Eingabedeklaration von `path reads` zu `tuple path(read1), path(read2)`
- Aktualisiere den Befehl im Block `script`, ersetze `-U $reads` durch `-1 ${read1} -2 ${read2}`
- Ersetze alle Vorkommen von `${reads.simpleName}` durch `${read1.simpleName}` im Befehl im Block `script` sowie in den Ausgabedeklarationen.

```groovy title="modules/hisat2_align_pe.nf" linenums="8"
    input:
    tuple path(read1), path(read2)
    path index_zip

    output:
    path "${read1.simpleName}.bam", emit: bam
    path "${read1.simpleName}.hisat2.log", emit: log

    script:
    """
    tar -xzvf $index_zip
    hisat2 -x ${index_zip.simpleName} -1 ${read1} -2 ${read2} \
        --new-summary --summary-file ${read1.simpleName}.hisat2.log | \
        samtools view -bS -o ${read1.simpleName}.bam
    """
```

Aktualisiere abschließend die Modul-Import-Anweisung, um die Paired-End-Version des Moduls zu verwenden.

```groovy title="rnaseq_pe.nf" linenums="5"
include { HISAT2_ALIGN } from './modules/hisat2_align_pe.nf'
```

### 3.8. Den Workflow ausführen, um zu testen, dass er funktioniert

Wir verwenden nicht `-resume`, da dies nicht cachen würde, und es gibt doppelt so viele Daten zu verarbeiten wie zuvor, aber es sollte trotzdem in unter einer Minute abgeschlossen sein.

```bash
nextflow run rnaseq_pe.nf
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

Und das war's! Jetzt haben wir zwei leicht unterschiedliche Versionen unseres Workflows, eine für Single-End-Read-Daten und eine für Paired-End-Daten.
Der nächste logische Schritt wäre, den Workflow so zu gestalten, dass er beide Datentypen spontan akzeptiert, was außerhalb des Umfangs dieses Kurses liegt, aber wir könnten das in einem Folgekurs angehen.

---

### Fazit

Du weißt jetzt, wie du einen Workflow für eine einzelne Probe anpasst, um die Verarbeitung mehrerer Proben zu parallelisieren, einen umfassenden QC-Bericht zu erstellen und den Workflow bei Bedarf für Paired-End-Read-Daten anzupassen.

### Wie geht es weiter?

Herzlichen Glückwunsch, du hast den Nextflow For RNAseq Mini-Kurs abgeschlossen! Feiere deinen Erfolg und gönn dir eine wohlverdiente Pause!

Als Nächstes bitten wir dich, eine sehr kurze Umfrage über deine Erfahrung mit diesem Trainingskurs auszufüllen. Danach leiten wir dich zu einer Seite mit Links zu weiteren Trainingsressourcen und hilfreichen Links weiter.
