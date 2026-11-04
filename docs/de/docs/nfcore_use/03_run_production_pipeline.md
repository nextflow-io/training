# Teil 3: Eine Produktions-Pipeline ausführen

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

In [Teil 2](./02_configure_execution.md) hast du gelernt, wie du Parameter setzt und die Konfiguration für nf-core/demo anpasst.
Jetzt wenden wir das Gelernte auf eine echte Produktions-Pipeline an: nf-core/rnaseq.

---

## 1. nf-core/rnaseq herunterladen und ausführen

Bisher haben wir `nf-core/demo` verwendet, eine minimale Pipeline, die für das Training entwickelt wurde.
Jetzt laden wir eine echte Produktions-Pipeline herunter und führen sie mit ihrem Test-Profil aus.

Die `nf-core/rnaseq` Pipeline führt die wichtigsten Schritte der Bulk-RNA-Sequenzierungsanalyse durch: Qualitätskontrolle, Adapter-Trimming, Read-Alignment und Quantifizierung auf Genebene.
Sie ist wahrscheinlich die am häufigsten verwendete nf-core Pipeline überhaupt.

### 1.1. Die Pipeline herunterladen

Führe den folgenden Befehl aus, um sie herunterzuladen.

```bash
nextflow pull nf-core/rnaseq
```

??? success "Befehlsausgabe"

    ```console
    Checking nf-core/rnaseq ...
     downloaded from https://github.com/nf-core/rnaseq.git - revision: 1f03b53ef7 [master]
    ```

Die Pipeline ist jetzt lokal gecacht und bereit zur Ausführung.

### 1.2. Das Test-Profil ausführen

Führe sie mit dem Test-Profil und Docker aus:

```bash
nextflow run nf-core/rnaseq -profile test,docker --outdir rnaseq-results
```

??? failure "Befehlsausgabe"

    ```console
     N E X T F L O W   ~  version 26.04.4

    Launching `https://github.com/nf-core/rnaseq` [suspicious_dijkstra] revision: 1f03b53ef7 [master]

    ------------------------------------------------------
                                            ,--./,-.
            ___     __   __   __   ___     /,-._.--~'
      |\ | |__  __ /  ` /  \ |__) |__         }  {
      | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                            `._,._,'
      nf-core/rnaseq 3.26.0
    ------------------------------------------------------
    ...
    [-        ] NFCORE_RNASEQ:RNASEQ:FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS:FQ_LINT  -
    Plus 47 more processes waiting for tasks…

    Execution cancelled -- Finishing pending tasks before exit
    -[nf-core/rnaseq] Pipeline completed with errors-
    ERROR ~ Error executing process > 'NFCORE_RNASEQ:RNASEQ:FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS:FQ_LINT (WT_REP2)'

    Caused by:
      Process requirement exceeds available memory -- req: 12 GB; avail: 7.7 GB


    Command executed:

      fq lint \
          --disable-validator P001 \
          SRR6357072_1.fastq.gz SRR6357072_2.fastq.gz > WT_REP2.fq_lint.txt

    Command exit status:
      -

    Command output:
      (empty)

    Work dir:
      /workspaces/training/nfcore-use/work/xx/xxxxxxxxxxxxxxxxxxxxxx

    Container:
      quay.io/biocontainers/fq:0.12.0--h9ee0642_0

    Tip: you can replicate the issue by changing to the process work dir and entering the command `bash .command.run`

     -- Check '.nextflow.log' file for details
    ERROR ~ Pipeline failed. Please refer to troubleshooting docs for common issues: https://nf-co.re/docs/running/troubleshooting
    ```

Die entscheidende Zeile in dieser Fehlermeldung ist:

```console
Process requirement exceeds available memory -- req: 12 GB; avail: 7.7 GB
```

Die Standard-Codespaces-Maschine hat 8 GB RAM – das ist auch der typische Standardwert für Docker Desktop.
Die Pipeline fordert 12 GB für den `FQ_LINT`-Prozess an – mehr als die Maschine bereitstellen kann.

Diese 12 GB stammen aus dem `process_low`-Ressourcen-Label, das in `conf/base.config` definiert ist:

```groovy title="conf/base.config"
withLabel:process_low {
    cpus   = { 2     * task.attempt }
    memory = { 12.GB * task.attempt }
    time   = { 4.h   * task.attempt }
}
```

Eine Möglichkeit wäre, einen größeren Maschinentyp zu verwenden. Für Testzwecke möchten wir jedoch auf beliebiger verfügbarer Hardware ausführen können.
Der bessere Ansatz ist, die Standard-Ressourcen in einer eigenen Konfigurationsdatei zu überschreiben.

### 1.3. Mit einer eigenen Konfiguration erneut ausführen

Wir stellen dir eine eigene Konfigurationsdatei zur Verfügung, die die label-basierten Ressourcen-Standardwerte überschreibt.

??? full-code "laptop.config"

    ```groovy title="laptop.config"
    process {
        withLabel: 'process_low' {
            cpus   = 2
            memory = 6.GB
        }
        withLabel: 'process_medium' {
            cpus   = 4
            memory = 6.GB
        }
        withLabel: 'process_high' {
            cpus   = 6
            memory = 6.GB
        }
        withLabel: 'process_high_memory' {
            memory = 6.GB
        }
    }
    ```

In [Teil 2](./02_configure_execution.md) haben wir `withName:` verwendet, um einen einzelnen Prozess nach Name anzusprechen.
Hier verwenden wir `withLabel:`, um alle Prozesse mit einem bestimmten Label auf einmal anzusprechen.

Diese Datei ist bereits in deinem Arbeitsverzeichnis vorhanden.
Übergib sie mit `-c`, um die Überschreibungen anzuwenden:

```bash
nextflow run nf-core/rnaseq -profile test,docker -c laptop.config --outdir rnaseq-results
```

??? success "Befehlsausgabe (Pipeline startet)"

    ```console
     N E X T F L O W   ~  version 26.04.4

    Launching `https://github.com/nf-core/rnaseq` [romantic_faraday] revision: 1f03b53ef7 [master]

    ------------------------------------------------------
                                            ,--./,-.
            ___     __   __   __   ___     /,-._.--~'
      |\ | |__  __ /  ` /  \ |__) |__         }  {
      | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                            `._,._,'
      nf-core/rnaseq 3.26.0
    ------------------------------------------------------
    ...
    executor >  local (7)
    [xx/xxxxxx] NFCORE_RNASEQ:RNASEQ:FQ_LINT (RAP1_IAA_30M_REP1)   | 3 of 5
    [xx/xxxxxx] NFCORE_RNASEQ:RNASEQ:FASTQC (RAP1_IAA_30M_REP1)    | 2 of 5
    ...
    ```

Die Pipeline läuft jetzt, und du kannst beobachten, wie die Aufgaben eine nach der anderen abgeschlossen werden.
Auf diesem minimalen Testdatensatz wird sie in 15–20 Minuten fertig sein und dabei über 200 Aufgaben ausführen.

Echte RNA-seq-Experimente umfassen typischerweise Dutzende von Proben und laufen stunden- oder tagelang.
Nextflow unterstützt HPC-Scheduler (SLURM, PBS, LSF) und Cloud-Plattformen (AWS, Google Cloud, Azure), die die Laufzeit durch die Verteilung der Arbeit auf viele Knoten erheblich reduzieren können.
Die Einrichtung dieser Umgebungen ist jedoch mit erheblichem Aufwand verbunden.

Die Seqera Platform (entwickelt von den Schöpfern von Nextflow) bietet eine webbasierte Oberfläche zum Starten von Nextflow-Pipelines auf HPC- oder Cloud-Infrastruktur – entweder deiner eigenen oder einer für dich verwalteten –, mit Rechen- und Datenverwaltungsfunktionen, die das Ausführen von Pipelines in großem Maßstab vereinfachen.

!!! tip "Tipp"

    Wissenschaftler\*innen an Hochschulen können die Seqera Platform kostenlos über das [Seqera Academic Program](https://seqera.io/academic-program/) nutzen.

### Fazit

Du hast `nf-core/rnaseq` heruntergeladen, gesehen, wie nf-core-Ressourcen-Labels funktionieren, und gelernt, sie mit einer eigenen Konfigurationsdatei zu überschreiben.
Noch wichtiger: Du hast verstanden, warum die lokale Ausführung eher ein Ausgangspunkt als ein Ziel für Analysen in echtem Maßstab ist.

### Wie geht es weiter?

Du hast die Grundlagen der Ausführung von nf-core-Pipelines kennengelernt.
Unter [Nächste Schritte](next_steps.md) erfährst du, wie es weitergeht.

---

## Zusammenfassung

In diesem Teil hast du gelernt, wie du:

- Eine Produktions-Pipeline (nf-core/rnaseq) herunterlädst und ausführst und ihre Standard-Ressourcen-Labels überschreibst
