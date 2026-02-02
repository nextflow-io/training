# Teil 1: Methodenübersicht und manuelles Testen

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Es gibt mehrere gültige Methoden zur Verarbeitung und Analyse von Bulk-RNAseq-Daten.
Für diesen Kurs folgen wir der Methode, die [hier](https://www.bioinformatics.babraham.ac.uk/training/RNASeq_Course/Analysing%20RNA-Seq%20data%20Exercise.pdf) von Dr. Simon Andrews und Dr. Laura Biggins am [Babraham Institute](https://www.babraham.ac.uk/) beschrieben wird.

Unser Ziel ist es, einen Workflow zu entwickeln, der die folgenden Verarbeitungsschritte implementiert: initiale Qualitätskontrolle der Reads in einer Bulk-RNAseq-Probe durchführen, Adaptersequenzen aus den Reads trimmen, die Reads auf ein Referenzgenom alignieren und einen umfassenden Qualitätskontroll-Bericht (QC) erstellen.

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/rnaseq/img/preprocess.svg"
</figure>

- **FASTQC:** QC der Read-Daten vor dem Trimmen mit FastQC durchführen
- **TRIM_GALORE:** Adaptersequenzen trimmen und QC nach dem Trimmen mit Trim Galore durchführen (bündelt Cutadapt und FastQC)
- **HISAT2_ALIGN:** Reads mit Hisat2 auf das Referenzgenom alignieren
- **MULTIQC:** Einen umfassenden QC-Bericht mit MultiQC generieren

Bevor wir jedoch mit dem Schreiben von Workflow-Code beginnen, werden wir die Befehle manuell an einigen Testdaten ausprobieren.
Die benötigten Tools sind in der GitHub Codespaces-Umgebung nicht installiert, daher werden wir sie über Container verwenden (siehe [Hello Containers](../../hello_nextflow/05_hello_containers.md)).

!!! note "Hinweis"

     Stelle sicher, dass du dich im Verzeichnis `nf4-science/rnaseq` befindest. Der letzte Teil des Pfads, der angezeigt wird, wenn du `pwd` eingibst, sollte `rnaseq` sein.

---

## 1. Initiale QC und Adapter-Trimming

Wir werden ein Container-Image herunterladen, das sowohl `fastqc` als auch `trim_galore` installiert hat, es interaktiv starten und die Trimming- und QC-Befehle auf einer der Beispieldateien ausführen.

### 1.1. Container herunterladen

```bash
docker pull community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18
```

Dies gibt dir die folgende Konsolenausgabe, während das System das Image herunterlädt:

??? success "Befehlsausgabe"

    ```console
    0.6.10--1bf8ca4e1967cd18: Pulling from library/trim-galore
    dafa2b0c44d2: Pull complete
    dec6b097362e: Pull complete
    f88da01cff0b: Pull complete
    4f4fb700ef54: Pull complete
    92dc97a3ef36: Pull complete
    403f74b0f85e: Pull complete
    10b8c00c10a5: Pull complete
    17dc7ea432cc: Pull complete
    bb36d6c3110d: Pull complete
    0ea1a16bbe82: Pull complete
    030a47592a0a: Pull complete
    32ec762be2d0: Pull complete
    d2cb90387285: Pull complete
    Digest: sha256:4f00e7b2a09f3c8d8a9ce955120e177152fb1e56f63a2a6e186088b1250d9907
    Status: Downloaded newer image for community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18
    community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18
    ```

### 1.2. Container interaktiv starten

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18
```

<!--
??? success "Befehlsausgabe"

    ```console

    ```
-->

Dein Prompt ändert sich zu etwas wie `(base) root@b645838b3314:/tmp#`, was anzeigt, dass du dich jetzt innerhalb des Containers befindest.

Der Teil `-v ./data:/data` des Befehls ermöglicht es uns, auf den Inhalt des Verzeichnisses `data/` von innerhalb des Containers zuzugreifen.

```bash
ls /data/reads
```

??? success "Befehlsausgabe"

    ```console
    ENCSR000COQ1_1.fastq.gz  ENCSR000COQ2_2.fastq.gz  ENCSR000COR2_1.fastq.gz  ENCSR000CPO1_2.fastq.gz
    ENCSR000COQ1_2.fastq.gz  ENCSR000COR1_1.fastq.gz  ENCSR000COR2_2.fastq.gz  ENCSR000CPO2_1.fastq.gz
    ENCSR000COQ2_1.fastq.gz  ENCSR000COR1_2.fastq.gz  ENCSR000CPO1_1.fastq.gz  ENCSR000CPO2_2.fastq.gzO
    ```

### 1.3. Ersten `fastqc`-Befehl ausführen

Lass uns `fastqc` ausführen, um Qualitätskontroll-Metriken der Read-Daten zu sammeln.

```bash
fastqc /data/reads/ENCSR000COQ1_1.fastq.gz
```

??? success "Befehlsausgabe"

    ```console
    application/gzip
    Started analysis of ENCSR000COQ1_1.fastq.gz
    Approx 5% complete for ENCSR000COQ1_1.fastq.gz
    Approx 10% complete for ENCSR000COQ1_1.fastq.gz
    Approx 15% complete for ENCSR000COQ1_1.fastq.gz
    Approx 20% complete for ENCSR000COQ1_1.fastq.gz
    Approx 25% complete for ENCSR000COQ1_1.fastq.gz
    Approx 30% complete for ENCSR000COQ1_1.fastq.gz
    Approx 35% complete for ENCSR000COQ1_1.fastq.gz
    Approx 40% complete for ENCSR000COQ1_1.fastq.gz
    Approx 45% complete for ENCSR000COQ1_1.fastq.gz
    Approx 50% complete for ENCSR000COQ1_1.fastq.gz
    Approx 55% complete for ENCSR000COQ1_1.fastq.gz
    Approx 60% complete for ENCSR000COQ1_1.fastq.gz
    Approx 65% complete for ENCSR000COQ1_1.fastq.gz
    Approx 70% complete for ENCSR000COQ1_1.fastq.gz
    Approx 75% complete for ENCSR000COQ1_1.fastq.gz
    Approx 80% complete for ENCSR000COQ1_1.fastq.gz
    Approx 85% complete for ENCSR000COQ1_1.fastq.gz
    Approx 90% complete for ENCSR000COQ1_1.fastq.gz
    Approx 95% complete for ENCSR000COQ1_1.fastq.gz
    Analysis complete for ENCSR000COQ1_1.fastq.gz
    ```

Dies sollte sehr schnell ablaufen.
Du findest die Ausgabedateien im selben Verzeichnis wie die ursprünglichen Daten:

```bash
ls /data/reads/ENCSR000COQ1_1_fastqc*
```

<!-- switch to tree -->

```console title="Ausgabe"
/data/reads/ENCSR000COQ1_1_fastqc.html  /data/reads/ENCSR000COQ1_1_fastqc.zip
```

### 1.4. Adaptersequenzen mit `trim_galore` trimmen

Jetzt führen wir `trim_galore` aus, das Cutadapt und FastQC bündelt, um die Adaptersequenzen zu trimmen und Post-Trimming-QC-Metriken zu sammeln.

```bash
trim_galore --fastqc /data/reads/ENCSR000COQ1_1.fastq.gz
```

Das Flag `--fastqc` bewirkt, dass der Befehl automatisch einen QC-Sammelschritt ausführt, nachdem das Trimmen abgeschlossen ist.

_Die Ausgabe ist sehr ausführlich, daher ist das Folgende gekürzt._

??? success "Befehlsausgabe"

    ```console
    Multicore support not enabled. Proceeding with single-core trimming.
    Path to Cutadapt set as: 'cutadapt' (default)
    Cutadapt seems to be working fine (tested command 'cutadapt --version')
    Cutadapt version: 4.9
    single-core operation.
    igzip command line interface 2.31.0
    igzip detected. Using igzip for decompressing

    <...>

    Analysis complete for ENCSR000COQ1_1_trimmed.fq.gz
    ```

Du findest die Ausgabedateien im Arbeitsverzeichnis:

```bash
ls ENCSR000COQ1_1*
```

```console title="Ausgabe"
ENCSR000COQ1_1.fastq.gz_trimming_report.txt  ENCSR000COQ1_1_trimmed_fastqc.html
ENCSR000COQ1_1_trimmed.fq.gz                 ENCSR000COQ1_1_trimmed_fastqc.zip
```

### 1.5. Ausgabedateien in das Dateisystem außerhalb des Containers verschieben

Alles, was im Container verbleibt, wird für zukünftige Arbeiten nicht zugänglich sein, also verschieben wir diese in ein neues Verzeichnis.

```bash
mkdir /data/trimmed
mv ENCSR000COQ1_1* /data/trimmed
```

### 1.6. Container beenden

```bash
exit
```

---

## 2. Reads auf das Referenzgenom alignieren

Wir werden ein Container-Image herunterladen, das `hisat2` installiert hat, es interaktiv starten und den Alignment-Befehl ausführen, um die RNAseq-Daten auf ein Referenzgenom zu alignieren.

### 2.1. `hisat2`-Container herunterladen

```bash
docker pull community.wave.seqera.io/library/hisat2_samtools:5e49f68a37dc010e
```

??? success "Befehlsausgabe"

    ```console
    Unable to find image 'community.wave.seqera.io/library/hisat2_samtools:5e49f68a37dc010e' locally
    5e49f68a37dc010e: Pulling from library/hisat2_samtools
    dafa2b0c44d2: Already exists
    dec6b097362e: Already exists
    f88da01cff0b: Already exists
    4f4fb700ef54: Already exists
    92dc97a3ef36: Already exists
    403f74b0f85e: Already exists
    10b8c00c10a5: Already exists
    17dc7ea432cc: Already exists
    bb36d6c3110d: Already exists
    0ea1a16bbe82: Already exists
    030a47592a0a: Already exists
    e74ed5dd390b: Pull complete
    abfcf0185e51: Pull complete
    Digest: sha256:29d8e1a3172a2bdde7be813f7ebec22d331388194a7c0de872b4ccca4bed8f45
    Status: Downloaded newer image for community.wave.seqera.io/library/hisat2_samtools:5e49f68a37dc010e
    ```

### 2.2. `hisat2`-Container interaktiv starten

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/hisat2_samtools:5e49f68a37dc010e
```

Der Befehl ist derselbe wie zuvor, mit der entsprechenden Container-URI ausgetauscht.

### 2.3. Hisat2-Genom-Indexdateien erstellen

Hisat2 benötigt die Genomreferenz in einem sehr spezifischen Format und kann nicht einfach die Datei `genome.fa` im FASTA-Format verarbeiten, die wir bereitstellen. Daher nutzen wir diese Gelegenheit, um die entsprechenden Ressourcen zu erstellen.

```bash
hisat2-build /data/genome.fa genome_index
```

Die Ausgabe ist sehr ausführlich, daher ist das Folgende gekürzt:

<!-- TODO: switch to full output -->

??? success "Befehlsausgabe"

    ```console
    Settings:
      Output files: "genome_index.*.ht2"
    <...>
    Total time for call to driver() for forward index: 00:00:16
    ```

Dies erstellt mehrere Genom-Indexdateien, die du im Arbeitsverzeichnis finden kannst.

```bash
ls genome_index.*
```

```console title="Ausgabe"
genome_index.1.ht2  genome_index.3.ht2  genome_index.5.ht2  genome_index.7.ht2
genome_index.2.ht2  genome_index.4.ht2  genome_index.6.ht2  genome_index.8.ht2
```

Wir werden diese gleich verwenden, aber zuerst erstellen wir einen gzippten Tarball mit diesen Genom-Indexdateien; wir werden sie später benötigen und das Generieren dieser Dateien ist normalerweise nichts, was wir als Teil eines Workflows tun möchten.

```bash
tar -czvf /data/genome_index.tar.gz genome_index.*
```

Dies speichert einen Tarball `genome_index.tar.gz`, der die Genom-Indexdateien enthält, im Verzeichnis `data/` auf unserem Dateisystem, was in Teil 2 dieses Kurses nützlich sein wird.

### 2.4. `hisat2`-Befehl ausführen

Jetzt können wir den Alignment-Befehl ausführen, der den Alignment-Schritt mit `hisat2` durchführt und dann die Ausgabe an `samtools` weiterleitet, um die Ausgabe als BAM-Datei zu schreiben.

Die Read-Dateneingabe ist die Datei `/data/trimmed/ENCSR000COQ1_1_trimmed.fq.gz`, die wir im vorherigen Schritt mit `trim_galore` generiert haben.

```bash
hisat2 -x genome_index -U /data/trimmed/ENCSR000COQ1_1_trimmed.fq.gz \
    --new-summary --summary-file ENCSR000COQ1_1_trimmed.hisat2.log | \
    samtools view -bS -o ENCSR000COQ1_1_trimmed.bam
```

??? success "Befehlsausgabe"

    ```console
    HISAT2 summary stats:
            Total reads: 27816
                    Aligned 0 time: 1550 (5.57%)
                    Aligned 1 time: 25410 (91.35%)
                    Aligned >1 times: 856 (3.08%)
            Overall alignment rate: 94.43%
    ```

Dies läuft fast sofort ab, weil es sich um eine sehr kleine Testdatei handelt.
In der Realität könnte dies bei entsprechender Größe viel länger dauern.

Wieder kannst du die Ausgabedateien im Arbeitsverzeichnis finden:

```bash
ls ENCSR000COQ1_1*
```

```console title="Ausgabe"
ENCSR000COQ1_1_trimmed.bam  ENCSR000COQ1_1_trimmed.hisat2.log
```

### 2.5. Ausgabedateien in das Dateisystem außerhalb des Containers verschieben

```bash
mkdir /data/aligned
mv ENCSR000COQ1_1* /data/aligned
```

### 2.6. Container beenden

```bash
exit
```

---

## 3. Umfassenden QC-Bericht generieren

Wir werden ein Container-Image herunterladen, das `multiqc` installiert hat, es interaktiv starten und einen Berichtsgenerierungsbefehl auf den Vorher/Nachher-FastQC-Berichtsdateien ausführen.

### 3.1. `multiqc`-Container herunterladen

```bash
docker pull community.wave.seqera.io/library/pip_multiqc:a3c26f6199d64b7c
```

??? success "Befehlsausgabe"

    ```console
    ad8f247edb55897c: Pulling from library/pip_multiqc
    dafa2b0c44d2: Already exists
    dec6b097362e: Already exists
    f88da01cff0b: Already exists
    4f4fb700ef54: Already exists
    92dc97a3ef36: Already exists
    403f74b0f85e: Already exists
    10b8c00c10a5: Already exists
    17dc7ea432cc: Already exists
    bb36d6c3110d: Already exists
    0ea1a16bbe82: Already exists
    030a47592a0a: Already exists
    3f229294c69a: Pull complete
    5a5ad47fd84c: Pull complete
    Digest: sha256:0ebb1d9605395a7df49ad0eb366b21f46afd96a5090376b0d8941cf5294a895a
    Status: Downloaded newer image for community.wave.seqera.io/library/pip_multiqc:a3c26f6199d64b7c
    community.wave.seqera.io/library/pip_multiqc:a3c26f6199d64b7c
    ```

### 3.2. `multiqc`-Container interaktiv starten

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/pip_multiqc:a3c26f6199d64b7c
```

### 3.3. `multiqc`-Befehl ausführen

```bash
multiqc /data/reads /data/trimmed /data/aligned -n ENCSR000COQ1_1_QC
```

??? success "Befehlsausgabe"

    ```console

    /// MultiQC 🔍 v1.27.1

          file_search | Search path: /data/reads
          file_search | Search path: /data/trimmed
          file_search | Search path: /data/aligned
            searching | ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 100% 20/20
                hisat2 | Found 1 reports
              cutadapt | Found 1 reports
                fastqc | Found 1 reports
        write_results | Data        : ENCSR000COQ1_1_QC_data
        write_results | Report      : ENCSR000COQ1_1_QC.html
              multiqc | MultiQC complete
    ```

MultiQC kann Verzeichnisse nach kompatiblen QC-Berichten durchsuchen und aggregiert alles, was es findet.

Hier sehen wir, dass das Tool alle drei QC-Berichte gefunden hat, die wir generiert haben: die initiale QC, die wir mit `fastqc` durchgeführt haben, den Post-Trimming-Bericht von `cutadapt` (erstellt über `trim_galore`) und die Post-Alignment-QC, die von `hisat2` produziert wurde.

Die Ausgabedateien befinden sich wieder im Arbeitsverzeichnis:

```bash
ls ENCSR000COQ1_1_QC*
```

```console title="Ausgabe"
ENCSR000COQ1_1_QC.html

ENCSR000COQ1_1_QC_data:
cutadapt_filtered_reads_plot.txt                     fastqc_top_overrepresented_sequences_table.txt
cutadapt_trimmed_sequences_plot_3_Counts.txt         hisat2_se_plot.txt
cutadapt_trimmed_sequences_plot_3_Obs_Exp.txt        multiqc.log
fastqc-status-check-heatmap.txt                      multiqc_citations.txt
fastqc_adapter_content_plot.txt                      multiqc_cutadapt.txt
fastqc_per_base_n_content_plot.txt                   multiqc_data.json
fastqc_per_base_sequence_quality_plot.txt            multiqc_fastqc.txt
fastqc_per_sequence_gc_content_plot_Counts.txt       multiqc_general_stats.txt
fastqc_per_sequence_gc_content_plot_Percentages.txt  multiqc_hisat2.txt
fastqc_per_sequence_quality_scores_plot.txt          multiqc_software_versions.txt
fastqc_sequence_counts_plot.txt                      multiqc_sources.txt
fastqc_sequence_duplication_levels_plot.txt
```

### 3.4. Ausgabedateien in das Dateisystem außerhalb des Containers verschieben

```bash
mkdir /data/final_qc
mv ENCSR000COQ1_1_QC** /data/final_qc
```

### 3.5. Container beenden

```bash
exit
```

---

### Zusammenfassung

Du hast alle einzelnen Befehle interaktiv in den entsprechenden Containern getestet.

### Wie geht es weiter?

Lerne, wie du dieselben Befehle in einen mehrstufigen Workflow verpackst, der Container zur Ausführung der Arbeit verwendet.
