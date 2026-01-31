# Orientierung

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Die Trainingsumgebung enthält alle Software, Code und Daten, die zum Durcharbeiten dieses Trainingskurses notwendig sind, sodass du nichts selbst installieren musst.
Allerdings benötigst du einen (kostenlosen) Account zum Anmelden, und du solltest dir ein paar Minuten Zeit nehmen, um dich mit der Benutzeroberfläche vertraut zu machen.

Falls du dies noch nicht getan hast, absolviere bitte den Mini-Kurs [Environment Setup](../../envsetup/), bevor du fortfährst.

## Bereitgestellte Materialien

Während dieses Trainingskurses arbeiten wir im Verzeichnis `nf4-science/rnaseq/`, in das du wechseln musst, wenn du den Trainingsarbeitsbereich öffnest.
Dieses Verzeichnis enthält alle Code-Dateien, Testdaten und zusätzlichen Dateien, die du benötigen wirst.

Erkunde gerne den Inhalt dieses Verzeichnisses; am einfachsten geht das mit dem Datei-Explorer auf der linken Seite des Trainingsarbeitsbereichs in der VSCode-Oberfläche.
Alternativ kannst du den Befehl `tree` verwenden.
Während des Kurses verwenden wir die Ausgabe von `tree`, um die Verzeichnisstruktur und den Inhalt in lesbarer Form darzustellen, manchmal mit geringfügigen Änderungen zur besseren Lesbarkeit.

Hier erstellen wir ein Inhaltsverzeichnis bis zur zweiten Ebene:

```bash
tree . -L 3
```

??? success "Verzeichnisinhalt"

    ```console
    rnaseq
    ├── data
    │   ├── genome.fa
    │   ├── paired-end.csv
    │   ├── reads
    │   │   ├── ENCSR000COQ1_1.fastq.gz
    │   │   ├── ENCSR000COQ1_2.fastq.gz
    │   │   ├── ENCSR000COQ2_1.fastq.gz
    │   │   ├── ENCSR000COQ2_2.fastq.gz
    │   │   ├── ENCSR000COR1_1.fastq.gz
    │   │   ├── ENCSR000COR1_2.fastq.gz
    │   │   ├── ENCSR000COR2_1.fastq.gz
    │   │   ├── ENCSR000COR2_2.fastq.gz
    │   │   ├── ENCSR000CPO1_1.fastq.gz
    │   │   ├── ENCSR000CPO1_2.fastq.gz
    │   │   ├── ENCSR000CPO2_1.fastq.gz
    │   │   └── ENCSR000CPO2_2.fastq.gz
    │   └── single-end.csv
    ├── nextflow.config
    ├── rnaseq.nf
    └── solutions
        ├── modules
        │   ├── fastqc.nf
        │   ├── fastqc_pe.nf
        │   ├── hisat2_align.nf
        │   ├── hisat2_align_pe.nf
        │   ├── multiqc.nf
        │   ├── trim_galore.nf
        │   └── trim_galore_pe.nf
        ├── rnaseq-2.1.nf
        ├── rnaseq-2.2.nf
        ├── rnaseq-2.3.nf
        ├── rnaseq-3.1.nf
        ├── rnaseq-3.2.nf
        └── rnaseq_pe-3.3.nf
    ```

!!!note "Hinweis"

    Keine Sorge, wenn das viel erscheint; wir gehen die relevanten Teile bei jedem Schritt des Kurses durch.
    Dies soll dir nur einen Überblick verschaffen.

**Hier ist eine Zusammenfassung dessen, was du zum Einstieg wissen solltest:**

- **Die Datei `rnaseq.nf`** ist die Grundstruktur des Workflow-Scripts, das wir entwickeln werden.

- **Die Datei `nextflow.config`** ist eine Konfigurationsdatei, die minimale Umgebungseigenschaften festlegt. Du kannst sie vorerst ignorieren.

- **Das Verzeichnis `data`** enthält Eingabedaten und zugehörige Ressourcen:

  - _Ein Referenzgenom_ namens `genome.fa`, das aus einer kleinen Region des menschlichen Chromosoms 20 besteht (aus hg19/b37).
  - _RNAseq-Daten_, die auf eine kleine Region reduziert wurden, um die Dateigrößen klein zu halten, im Verzeichnis `reads/`.
  - _CSV-Dateien_, die die IDs und Pfade der Beispieldatendateien auflisten, zur Stapelverarbeitung.

- **Das Verzeichnis `solutions`** enthält die fertigen Workflow-Scripts und Module, die aus jedem Schritt des Kurses resultieren.
  Sie sind als Referenz gedacht, um deine Arbeit zu überprüfen und eventuelle Probleme zu beheben.
  Die Nummer im Dateinamen entspricht dem Schritt des relevanten Kursteils.

!!!tip "Tipp"

    Falls du aus irgendeinem Grund dieses Verzeichnis verlässt, kannst du jederzeit diesen Befehl ausführen, um dorthin zurückzukehren:

    ```bash
    cd /workspaces/training/nf4-science/rnaseq
    ```

Um nun mit dem Kurs zu beginnen, klicke auf den Pfeil in der unteren rechten Ecke dieser Seite.
