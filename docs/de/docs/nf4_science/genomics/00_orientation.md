# Erste Schritte

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

## Eine Trainingsumgebung starten

Um die vorgefertigte Umgebung zu nutzen, die wir auf GitHub Codespaces bereitstellen, klicke unten auf den Button „Open in GitHub Codespaces". Für andere Optionen siehe [Umgebungsoptionen](../../envsetup/index.md).

Wir empfehlen, die Trainingsumgebung in einem neuen Browser-Tab oder -Fenster zu öffnen (verwende Rechtsklick, Strg+Klick oder Cmd+Klick, je nach deinem Gerät), damit du diese Anleitung parallel zur ladenden Umgebung lesen kannst.
Du musst diese Anleitung während des Kurses geöffnet lassen.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

### Grundlagen der Umgebung

Diese Trainingsumgebung enthält alle Software, den Code und die Daten, die für diesen Trainingskurs notwendig sind, sodass du nichts selbst installieren musst.

Der Codespace ist mit einer VSCode-Oberfläche eingerichtet, die einen Datei-Explorer, einen Code-Editor und ein Terminal-Shell enthält.
Alle Anweisungen während des Kurses (z. B. „öffne die Datei", „bearbeite den Code" oder „führe diesen Befehl aus") beziehen sich auf diese drei Teile der VSCode-Oberfläche, sofern nicht anders angegeben.

Falls du diesen Kurs selbstständig durcharbeitest, mach dich bitte mit den [Grundlagen der Umgebung](../../envsetup/01_setup.md) vertraut, um weitere Details zu erfahren.

### Versionsanforderungen

Dieses Training ist für Nextflow 25.10.2 oder neuer **mit AKTIVIERTEM v2-Syntax-Parser** konzipiert.
Falls du eine lokale oder benutzerdefinierte Umgebung verwendest, stelle bitte sicher, dass du die korrekten Einstellungen verwendest, wie [hier](../../info/nxf_versions.md) dokumentiert.

## Vorbereitung zur Arbeit

Sobald dein Codespace läuft, musst du zwei Dinge tun, bevor du ins Training einsteigst: Setze dein Arbeitsverzeichnis für diesen spezifischen Kurs und wirf einen Blick auf die bereitgestellten Materialien.

### Das Arbeitsverzeichnis festlegen

Standardmäßig öffnet sich der Codespace mit dem Arbeitsverzeichnis am Root aller Trainingskurse, aber für diesen Kurs arbeiten wir im Verzeichnis `nf4-science/genomics/`.

Wechsle jetzt das Verzeichnis, indem du diesen Befehl im Terminal ausführst:

```bash
cd nf4-science/genomics/
```

Du kannst VSCode so einstellen, dass es sich auf dieses Verzeichnis fokussiert, sodass nur die relevanten Dateien im Datei-Explorer in der Seitenleiste angezeigt werden:

```bash
code .
```

!!! tip "Tipp"

    Falls du aus irgendeinem Grund aus diesem Verzeichnis herausnavigierst (z. B. wenn dein Codespace in den Ruhezustand geht), kannst du jederzeit mit dem vollständigen Pfad zurückkehren, vorausgesetzt du arbeitest in der GitHub Codespaces Trainingsumgebung:

    ```bash
    cd /workspaces/training/nf4-science/genomics
    ```

Schauen wir uns nun den Inhalt an.

### Die bereitgestellten Materialien erkunden

Du kannst den Inhalt dieses Verzeichnisses mit dem Datei-Explorer auf der linken Seite des Trainings-Workspace erkunden.
Alternativ kannst du den Befehl `tree` verwenden.

Im Laufe des Kurses nutzen wir die Ausgabe von `tree`, um die Verzeichnisstruktur und den Inhalt in lesbarer Form darzustellen, manchmal mit kleinen Änderungen zur besseren Übersicht.

Hier erstellen wir ein Inhaltsverzeichnis bis zur zweiten Ebene:

```bash
tree . -L 2
```

??? abstract "Verzeichnisinhalt"

    ```console
    .
    ├── data
    │   ├── bam
    │   ├── ref
    │   ├── sample_bams.txt
    │   └── samplesheet.csv
    ├── genomics.nf
    ├── modules
    │   ├── gatk_haplotypecaller.nf
    │   └── samtools_index.nf
    ├── nextflow.config
    └── solutions
        ├── modules
        ├── nf-test.config
        ├── part2
        └── tests

    8 directories, 8 files
    ```

Klicke auf das farbige Feld, um den Abschnitt zu erweitern und seinen Inhalt anzuzeigen.
Wir verwenden klappbare Abschnitte wie diesen, um erwartete Befehlsausgaben sowie Verzeichnis- und Dateiinhalte kompakt darzustellen.

- **Die Datei `genomics.nf`** ist ein Workflow-Skript, das du im Verlauf des Kurses aufbaust.

- **Das Verzeichnis `modules`** enthält Skeleton-Moduldateien, die du während des Kurses ausfüllst.

- **Die Datei `nextflow.config`** ist eine Konfigurationsdatei, die minimale Umgebungseigenschaften festlegt.
  Du kannst sie vorerst ignorieren.

- **Das Verzeichnis `data`** enthält Eingabedaten und zugehörige Ressourcen, die später im Kurs beschrieben werden.

- **Das Verzeichnis `solutions`** enthält fertige Moduldateien und eine Lösung für Teil 2, die als Ausgangspunkt für Teil 3 dienen kann.
  Sie sind als Referenz gedacht, um deine Arbeit zu überprüfen und eventuelle Probleme zu beheben.

## Checkliste zur Bereitschaft

Glaubst du, du bist bereit einzusteigen?

- [ ] Ich verstehe das Ziel dieses Kurses und seine Voraussetzungen
- [ ] Meine Umgebung läuft
- [ ] Ich habe mein Arbeitsverzeichnis entsprechend gesetzt

Wenn du alle Kästchen abhaken kannst, kann es losgehen.

**Um mit [Teil 1: Methodenübersicht und manuelle Tests](./01_method.md) fortzufahren, klicke auf den Pfeil in der unteren rechten Ecke dieser Seite.**
