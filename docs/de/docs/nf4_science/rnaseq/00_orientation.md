# Orientierung

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Die Trainingsumgebung enthält alle Software, Code und Daten, die du für diesen Kurs benötigst. Du musst also nichts selbst installieren.
Du brauchst jedoch einen (kostenlosen) Account zum Einloggen und solltest dir ein paar Minuten Zeit nehmen, um dich mit der Oberfläche vertraut zu machen.

Falls du es noch nicht getan hast, absolviere bitte den Mini-Kurs [Umgebung einrichten](../../envsetup/), bevor du weitermachst.

## Bereitgestellte Materialien

In diesem Kurs arbeiten wir im Verzeichnis `nf4-science/rnaseq/`, in das du wechseln musst, wenn du den Trainings-Workspace öffnest.
Dieses Verzeichnis enthält alle Code-Dateien, Testdaten und zusätzlichen Dateien, die du benötigst.

Erkunde gerne den Inhalt dieses Verzeichnisses. Am einfachsten geht das über den Datei-Explorer auf der linken Seite des Trainings-Workspace in der VSCode-Oberfläche.
Alternativ kannst du den Befehl `tree` verwenden.
Im gesamten Kurs nutzen wir die Ausgabe von `tree`, um die Verzeichnisstruktur und Inhalte in lesbarer Form darzustellen, manchmal mit kleinen Anpassungen für bessere Übersichtlichkeit.

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

    Keine Sorge, wenn das viel erscheint. Wir gehen die relevanten Teile bei jedem Schritt des Kurses durch.
    Dies soll dir nur einen Überblick geben.

**Hier ist eine Zusammenfassung dessen, was du zum Start wissen solltest:**

- **Die Datei `rnaseq.nf`** ist der Entwurf des Workflow-Skripts, das wir entwickeln werden.

- **Die Datei `nextflow.config`** ist eine Konfigurationsdatei, die minimale Umgebungseigenschaften festlegt. Du kannst sie vorerst ignorieren.

- **Das Verzeichnis `data`** enthält Eingabedaten und zugehörige Ressourcen:

  - _Ein Referenzgenom_ namens `genome.fa`, das aus einem kleinen Bereich des menschlichen Chromosoms 20 besteht (von hg19/b37).
  - _RNAseq-Daten_, die auf einen kleinen Bereich reduziert wurden, um die Dateigrößen klein zu halten, im Verzeichnis `reads/`.
  - _CSV-Dateien_, die die IDs und Pfade der Beispieldateien auflisten, für die Stapelverarbeitung.

- **Das Verzeichnis `solutions`** enthält die fertigen Workflow-Skripte und Module, die aus jedem Schritt des Kurses resultieren.
  Sie dienen als Referenz, um deine Arbeit zu überprüfen und Probleme zu beheben.
  Die Nummer im Dateinamen entspricht dem Schritt des jeweiligen Kursteils.

!!!tip "Tipp"

    Falls du aus irgendeinem Grund dieses Verzeichnis verlässt, kannst du jederzeit mit diesem Befehl zurückkehren:

    ```bash
    cd /workspaces/training/nf4-science/rnaseq
    ```

Um mit dem Kurs zu beginnen, klicke auf den Pfeil unten rechts auf dieser Seite.
