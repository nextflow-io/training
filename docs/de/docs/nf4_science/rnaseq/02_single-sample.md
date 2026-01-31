# Teil 2: Implementierung für eine einzelne Probe

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

In diesem Teil des Kurses werden wir den einfachstmöglichen Workflow schreiben, der alle Befehle aus Teil 1 zusammenfasst, um ihre Ausführung zu automatisieren. Dabei werden wir zunächst nur eine Probe gleichzeitig verarbeiten.

Wir werden dies in drei Schritten tun:

1. Einen einstufigen Workflow schreiben, der den ersten QC-Schritt ausführt
2. Adapter-Trimming und Post-Trimming-QC hinzufügen
3. Alignment zum Referenzgenom hinzufügen

!!! warning "Voraussetzung"

    Du musst Teil 1 des Kurses durcharbeiten, bevor du mit dieser Lektion beginnst.
    Insbesondere durch das Durcharbeiten der Abschnitte 2.1-3 wird die Genom-Indexdatei (`data/genome_index.tar.gz`) erstellt, die für den Alignment-Schritt in dieser Lektion benötigt wird.

---

## 1. Einen einstufigen Workflow schreiben, der die erste QC ausführt

Beginnen wir damit, einen einfachen Workflow zu schreiben, der das FastQC-Tool auf eine FASTQ-Datei mit Single-End-RNAseq-Reads anwendet.

Wir stellen dir eine Workflow-Datei, `rnaseq.nf`, zur Verfügung, die die Hauptteile des Workflows skizziert.

```groovy title="rnaseq.nf" linenums="1"
#!/usr/bin/env nextflow

// Modul-INCLUDE-Anweisungen

/*
 * Pipeline parameters
 */

// Primäre Eingabe

workflow {

    // Eingabe-Channel erstellen

    // Prozesse aufrufen

}
```

Beachte, dass dieser Workflow-Code zwar korrekt, aber noch nicht funktionsfähig ist; sein Zweck ist es lediglich, als Gerüst zu dienen, das du zum Schreiben des eigentlichen Workflows verwenden wirst.

### 1.1. Ein Verzeichnis zum Speichern von Modulen erstellen

Wir werden eigenständige Module für jeden Prozess erstellen, um sie einfacher verwalten und wiederverwenden zu können. Erstellen wir also ein Verzeichnis zum Speichern.

```bash
mkdir modules
```

### 1.2. Ein Modul für den Prozess zur Erfassung von QC-Metriken erstellen

Erstellen wir eine Moduldatei namens `modules/fastqc.nf`, um den `FASTQC`-Prozess unterzubringen:

```bash
touch modules/fastqc.nf
```

Öffne die Datei im Code-Editor und kopiere den folgenden Code hinein:

```groovy title="modules/fastqc.nf" linenums="1"
#!/usr/bin/env nextflow

process FASTQC {

    container "community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18"
    publishDir "results/fastqc", mode: 'symlink'

    input:
    path reads

    output:
    path "${reads.simpleName}_fastqc.zip", emit: zip
    path "${reads.simpleName}_fastqc.html", emit: html

    script:
    """
    fastqc $reads
    """
}
```

Du solltest alle Teile aus dem wiedererkennen, was du in Teil 1 & Teil 2 dieser Trainingsreihe gelernt hast; die einzige nennenswerte Änderung ist, dass wir diesmal `mode: symlink` für die `publishDir`-Direktive verwenden und einen Parameter zur Definition von `publishDir` verwenden.

!!! note "Hinweis"

    Obwohl die Datendateien, die wir hier verwenden, sehr klein sind, können sie in der Genomik sehr groß werden. Zu Demonstrationszwecken in der Trainingsumgebung verwenden wir den 'symlink'-Veröffentlichungsmodus, um unnötige Dateikopien zu vermeiden. Du solltest dies in deinen finalen Workflows nicht tun, da du Ergebnisse verlierst, wenn du dein `work`-Verzeichnis aufräumst.

### 1.3. Das Modul in die Workflow-Datei importieren

Füge die Anweisung `include { FASTQC } from './modules/fastqc.nf'` zur Datei `rnaseq.nf` hinzu:

```groovy title="rnaseq.nf" linenums="3"
// Modul-INCLUDE-Anweisungen
include { FASTQC } from './modules/fastqc.nf'
```

### 1.4. Eine Eingabedeklaration hinzufügen

Deklariere einen Eingabeparameter mit einem Standardwert:

```groovy title="rnaseq.nf" linenums="10"
params {
    // Primäre Eingabe
    reads: Path = "data/reads/ENCSR000COQ1_1.fastq.gz"
}
```

### 1.5. Einen Eingabekanal im Workflow-Block erstellen

Verwende eine einfache `.fromPath()`-Channel-Factory, um den Eingabekanal zu erstellen:

```groovy title="rnaseq.nf" linenums="13"
workflow {

    // Eingabe-Channel aus einem Dateipfad erstellen
    read_ch = channel.fromPath(params.reads)

    // Prozesse aufrufen

}
```

### 1.6. Den `FASTQC`-Prozess auf dem Eingabekanal aufrufen

```groovy title="rnaseq.nf" linenums="13"
workflow {

    // Eingabe-Channel aus einem Dateipfad erstellen
    read_ch = channel.fromPath(params.reads)

    // Initiale Qualitätskontrolle
    FASTQC(read_ch)

}
```

### 1.7. Den Workflow ausführen, um zu testen, ob er funktioniert

Wir könnten den `--reads`-Parameter verwenden, um eine Eingabe von der Befehlszeile aus anzugeben, aber während der Entwicklung können wir faul sein und einfach den eingerichteten Teststandard verwenden.

```bash
nextflow run rnaseq.nf
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 24.10.0

    Launching `rnaseq.nf` [fabulous_snyder] DSL2 - revision: 3394c725ee

    executor >  local (1)
    [d6/d94c3a] FASTQC (1) [100%] 1 of 1 ✔
    ```

Dies sollte sehr schnell ausgeführt werden, wenn du Teil 1 durchgearbeitet hast und den Container bereits heruntergeladen hast.
Wenn du ihn übersprungen hast, wird Nextflow den Container für dich herunterladen; du musst nichts dafür tun, aber du musst möglicherweise bis zu einer Minute warten.

Du findest die Ausgaben unter `results/fastqc`, wie im `FASTQC`-Prozess durch die `publishDir`-Direktive angegeben.

```bash
ls results/fastqc
```

```console title="Output"
ENCSR000COQ1_1_fastqc.html  ENCSR000COQ1_1_fastqc.zip
```

---

## 2. Adapter-Trimming und Post-Trimming-Qualitätskontrolle hinzufügen

Wir werden den Trim_Galore-Wrapper verwenden, der Cutadapt für das Trimming selbst und FastQC für die Post-Trimming-Qualitätskontrolle bündelt.

### 2.1. Ein Modul für den Trimming- und QC-Prozess erstellen

Erstellen wir eine Moduldatei namens `modules/trim_galore.nf`, um den `TRIM_GALORE`-Prozess unterzubringen:

```bash
touch modules/trim_galore.nf
```

Öffne die Datei im Code-Editor und kopiere den folgenden Code hinein:

```groovy title="modules/trim_galore.nf" linenums="1"
#!/usr/bin/env nextflow

process TRIM_GALORE {

    container "community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18"
    publishDir "results/trimming", mode: 'symlink'

    input:
    path reads

    output:
    path "${reads.simpleName}_trimmed.fq.gz", emit: trimmed_reads
    path "${reads}_trimming_report.txt", emit: trimming_reports
    path "${reads.simpleName}_trimmed_fastqc.{zip,html}", emit: fastqc_reports

    script:
    """
    trim_galore --fastqc $reads
    """
}
```

### 2.2. Das Modul in die Workflow-Datei importieren

Füge die Anweisung `include { TRIM_GALORE } from './modules/trim_galore.nf'` zur Datei `rnaseq.nf` hinzu:

```groovy title="rnaseq.nf" linenums="3"
// Modul-INCLUDE-Anweisungen
include { FASTQC } from './modules/fastqc.nf'
include { TRIM_GALORE } from './modules/trim_galore.nf'
```

### 2.3. Den Prozess auf dem Eingabekanal aufrufen

```groovy title="rnaseq.nf" linenums="14"
workflow {

    // Eingabe-Channel aus einem Dateipfad erstellen
    read_ch = channel.fromPath(params.reads)

    // Initiale Qualitätskontrolle
    FASTQC(read_ch)

    // Adapter-Trimming und Post-Trimming-QC
    TRIM_GALORE(read_ch)
}
```

### 2.4. Den Workflow ausführen, um zu testen, ob er funktioniert

```bash
nextflow run rnaseq.nf
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 24.10.0

    Launching `rnaseq.nf` [fabulous_snyder] DSL2 - revision: 3394c725ee

    executor >  local (1)
    [d6/d94c3a] FASTQC (1) [100%] 1 of 1 ✔
    [c2/e4a9bb] TRIM_GALORE (1)  [100%] 1 of 1 ✔
    ```

Dies sollte ebenfalls sehr schnell ausgeführt werden, da wir mit einer so kleinen Eingabedatei arbeiten.

Du findest die Ausgaben unter `results/trimming`, wie im `TRIM_GALORE`-Prozess durch die `publishDir`-Direktive angegeben.

```bash
ls results/trimming
```

```console title="Output"
ENCSR000COQ1_1.fastq.gz_trimming_report.txt  ENCSR000COQ1_1_trimmed_fastqc.zip
ENCSR000COQ1_1_trimmed_fastqc.html           ENCSR000COQ1_1_trimmed.fq.gz
```

---

## 3. Die Reads zum Referenzgenom alignieren

Schließlich können wir den Genom-Alignment-Schritt mit Hisat2 ausführen, der auch FastQC-ähnliche Qualitätskontrollmetriken ausgibt.

### 3.1. Ein Modul für den HiSat2-Prozess erstellen

Erstellen wir eine Moduldatei namens `modules/hisat2_align.nf`, um den `HISAT2_ALIGN`-Prozess unterzubringen:

```bash
touch modules/hisat2_align.nf
```

Öffne die Datei im Code-Editor und kopiere den folgenden Code hinein:

```groovy title="modules/hisat2_align.nf" linenums="1"
#!/usr/bin/env nextflow

process HISAT2_ALIGN {

    container "community.wave.seqera.io/library/hisat2_samtools:5e49f68a37dc010e"
    publishDir "results/align", mode: 'symlink'

    input:
    path reads
    path index_zip

    output:
    path "${reads.simpleName}.bam", emit: bam
    path "${reads.simpleName}.hisat2.log", emit: log

    script:
    """
    tar -xzvf $index_zip
    hisat2 -x ${index_zip.simpleName} -U $reads \
        --new-summary --summary-file ${reads.simpleName}.hisat2.log | \
        samtools view -bS -o ${reads.simpleName}.bam
    """
}
```

### 3.2. Das Modul in die Workflow-Datei importieren

Füge die Anweisung `include { HISAT2_ALIGN } from './modules/hisat2_align.nf'` zur Datei `rnaseq.nf` hinzu:

```groovy title="rnaseq.nf" linenums="3"
// Modul-INCLUDE-Anweisungen
include { FASTQC } from './modules/fastqc.nf'
include { TRIM_GALORE } from './modules/trim_galore.nf'
include { HISAT2_ALIGN } from './modules/hisat2_align.nf'
```

### 3.3. Eine Parameterdeklaration hinzufügen, um den Genomindex bereitzustellen

Deklariere einen Eingabeparameter mit einem Standardwert:

```groovy title="rnaseq.nf" linenums="8"
params {
    // Primäre Eingabe
    reads: Path = "data/reads/ENCSR000COQ1_1.fastq.gz"

    // Referenzgenom-Archiv
    hisat2_index_zip: Path = "data/genome_index.tar.gz"
}
```

### 3.4. Den `HISAT2_ALIGN`-Prozess auf den getrimmten Reads aufrufen, die von `TRIM_GALORE` ausgegeben wurden

Die getrimmten Reads befinden sich im `TRIM_GALORE.out.trimmed_reads`-Channel, der vom vorherigen Schritt ausgegeben wurde.

Zusätzlich verwenden wir `file (params.hisat2_index_zip)`, um dem Hisat2-Tool das gezippte Genom-Index-Tarball bereitzustellen.

```groovy title="rnaseq.nf" linenums="16"
workflow {

    // Eingabe-Channel aus einem Dateipfad erstellen
    read_ch = channel.fromPath(params.reads)

    // Initiale Qualitätskontrolle
    FASTQC(read_ch)

    // Adapter-Trimming und Post-Trimming-QC
    TRIM_GALORE(read_ch)

    // Alignment zum Referenzgenom
    HISAT2_ALIGN(TRIM_GALORE.out.trimmed_reads, file (params.hisat2_index_zip))
}
```

### 3.5. Den Workflow ausführen, um zu testen, ob er funktioniert

```bash
nextflow run rnaseq.nf
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 24.10.0

    Launching `rnaseq.nf` [extravagant_khorana] DSL2 - revision: 701b41bd16

    executor >  local (3)
    [e4/d15ad4] FASTQC (1)       [100%] 1 of 1 ✔
    [c6/12b2be] TRIM_GALORE (1)  [100%] 1 of 1 ✔
    [c6/7a9f13] HISAT2_ALIGN (1) [100%] 1 of 1 ✔
    ```

Du findest die Ausgaben unter `results/align`, wie im `HISAT2_ALIGN`-Prozess durch die `publishDir`-Direktive angegeben.

```bash
ls results/align
```

```console title="Output"
ENCSR000COQ1_1_trimmed.bam  ENCSR000COQ1_1_trimmed.hisat2.log
```

Dies vervollständigt die grundlegende Verarbeitung, die wir auf jede Probe anwenden müssen.

_Wir werden die MultiQC-Report-Aggregation in Teil 2 hinzufügen, nachdem wir den Workflow so angepasst haben, dass er mehrere Proben gleichzeitig akzeptiert._

---

### Zusammenfassung

Du weißt jetzt, wie du alle Kernschritte zur individuellen Verarbeitung von Single-End-RNAseq-Proben zusammenfasst.

### Wie geht es weiter?

Lerne, wie du den Workflow anpasst, um mehrere Proben parallel zu verarbeiten, QC-Reports über alle Schritte für alle Proben zu aggregieren und die Ausführung des Workflows mit Paired-End-RNAseq-Daten zu ermöglichen.
