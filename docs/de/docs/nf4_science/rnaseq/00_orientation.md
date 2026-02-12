# Erste Schritte

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

## Eine Trainingsumgebung starten

Um die vorgefertigte Umgebung zu nutzen, die wir auf GitHub Codespaces bereitstellen, klicke auf den Button "Open in GitHub Codespaces" unten. Für andere Optionen siehe [Umgebungsoptionen](../../envsetup/index.md).

Wir empfehlen, die Trainingsumgebung in einem neuen Browser-Tab oder -Fenster zu öffnen (verwende Rechtsklick, Strg-Klick oder Cmd-Klick, je nach deiner Ausstattung), damit du weiterlesen kannst, während die Umgebung lädt.
Du musst diese Anleitung parallel geöffnet halten, um den Kurs durchzuarbeiten.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

### Grundlagen der Umgebung

Diese Trainingsumgebung enthält alle Software, Code und Daten, die zum Durcharbeiten des Trainingskurses notwendig sind, sodass du nichts selbst installieren musst.

Der Codespace ist mit einer VSCode-Oberfläche eingerichtet, die einen Dateisystem-Explorer, einen Code-Editor und ein Terminal-Shell umfasst.
Alle Anweisungen während des Kurses (z.B. 'öffne die Datei', 'bearbeite den Code' oder 'führe diesen Befehl aus') beziehen sich auf diese drei Teile der VSCode-Oberfläche, sofern nicht anders angegeben.

Falls du diesen Kurs selbstständig durcharbeitest, mache dich bitte mit den [Grundlagen der Umgebung](../../envsetup/01_setup.md) vertraut, um weitere Details zu erfahren.

### Versionsanforderungen

Dieses Training ist für Nextflow 25.10.2 oder später **mit AKTIVIERTEM v2-Syntax-Parser** konzipiert.
Falls du eine lokale oder benutzerdefinierte Umgebung verwendest, stelle bitte sicher, dass du die korrekten Einstellungen verwendest, wie [hier](../../info/nxf_versions.md) dokumentiert.

## Bereit zum Arbeiten

Sobald dein Codespace läuft, musst du zwei Dinge tun, bevor du ins Training eintauchst: dein Arbeitsverzeichnis für diesen spezifischen Kurs festlegen und einen Blick auf die bereitgestellten Materialien werfen.

### Das Arbeitsverzeichnis festlegen

Standardmäßig öffnet sich der Codespace mit dem Arbeitsverzeichnis im Stammverzeichnis aller Trainingskurse, aber für diesen Kurs arbeiten wir im Verzeichnis `nf4-science/rnaseq/`.

Wechsle jetzt das Verzeichnis, indem du diesen Befehl im Terminal ausführst:

```bash
cd nf4-science/rnaseq/
```

Du kannst VSCode so einstellen, dass es sich auf dieses Verzeichnis fokussiert, sodass nur die relevanten Dateien im Datei-Explorer in der Seitenleiste angezeigt werden:

```bash
code .
```

!!! tip "Tipp"

    Falls du aus irgendeinem Grund dieses Verzeichnis verlässt (z.B. dein Codespace geht in den Ruhezustand), kannst du jederzeit den vollständigen Pfad verwenden, um dorthin zurückzukehren, vorausgesetzt du arbeitest in der GitHub Codespaces Trainingsumgebung:

    ```bash
    cd /workspaces/training/nf4-science/rnaseq
    ```

Schauen wir uns nun den Inhalt an.

### Die bereitgestellten Materialien erkunden

Du kannst den Inhalt dieses Verzeichnisses mit dem Datei-Explorer auf der linken Seite des Trainingsarbeitsbereichs erkunden.
Alternativ kannst du den Befehl `tree` verwenden.

Während des Kurses verwenden wir die Ausgabe von `tree`, um die Verzeichnisstruktur und den Inhalt in lesbarer Form darzustellen, manchmal mit geringfügigen Änderungen zur besseren Lesbarkeit.

Hier erstellen wir ein Inhaltsverzeichnis bis zur dritten Ebene:

```bash
tree . -L 3
```

??? abstract "Verzeichnisinhalt"

    ```console
    .
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

Klicke auf die farbige Box, um den Abschnitt zu erweitern und seinen Inhalt anzuzeigen.
Wir verwenden solche ausklappbaren Abschnitte, um erwartete Befehlsausgaben sowie Verzeichnis- und Dateiinhalte auf kompakte Weise darzustellen.

- **Die Datei `rnaseq.nf`** ist eine Grundstruktur für ein Workflow-Script, das du im Verlauf des Kurses aufbauen wirst.

- **Das Verzeichnis `modules`** enthält Grundstrukturen für Prozessmodule, die du während des Kurses ausfüllen wirst.

- **Die Datei `nextflow.config`** ist eine Konfigurationsdatei, die minimale Umgebungseigenschaften festlegt.
  Du kannst sie vorerst ignorieren.

- **Das Verzeichnis `data`** enthält Eingabedaten und zugehörige Ressourcen, die später im Kurs beschrieben werden.

- **Das Verzeichnis `solutions`** enthält die fertigen Workflow-Scripts und Module, die aus jedem Schritt des Kurses resultieren.
  Sie sind als Referenz gedacht, um deine Arbeit zu überprüfen und eventuelle Probleme zu beheben.
  Die Lösung aus Teil 2 kann als Ausgangspunkt für Teil 3 verwendet werden.

## Bereitschafts-Checkliste

Denkst du, du bist bereit loszulegen?

- [ ] Ich verstehe das Ziel dieses Kurses und seine Voraussetzungen
- [ ] Meine Umgebung ist eingerichtet und läuft
- [ ] Ich habe mein Arbeitsverzeichnis entsprechend festgelegt

Wenn du alle Kästchen abhaken kannst, bist du startklar.

**Um mit [Teil 1: Methodenübersicht](./01_method.md) fortzufahren, klicke auf den Pfeil in der unteren rechten Ecke dieser Seite.**
