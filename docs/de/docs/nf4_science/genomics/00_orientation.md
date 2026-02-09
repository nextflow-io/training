# Erste Schritte

## Eine Trainingsumgebung starten

Um die vorgefertigte Umgebung zu nutzen, die wir auf GitHub Codespaces bereitstellen, klicke auf die Schaltfläche „In GitHub Codespaces öffnen" unten. Für andere Optionen siehe [Umgebungsoptionen](../../envsetup/index.md).

Wir empfehlen, die Trainingsumgebung in einem neuen Browser-Tab oder -Fenster zu öffnen (verwende Rechtsklick, Strg-Klick oder Cmd-Klick, je nach deinem Gerät), damit du weiterlesen kannst, während die Umgebung lädt.
Du musst diese Anleitung parallel geöffnet halten, um den Kurs durchzuarbeiten.

[![In GitHub Codespaces öffnen](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

### Grundlagen der Umgebung

Diese Trainingsumgebung enthält alle Software, den Code und die Daten, die für den Kurs notwendig sind. Du musst also nichts selbst installieren.

Der Codespace ist mit einer VSCode-Oberfläche eingerichtet, die einen Dateisystem-Explorer, einen Code-Editor und ein Terminal-Shell umfasst.
Alle Anweisungen im Kurs (z. B. „öffne die Datei", „bearbeite den Code" oder „führe diesen Befehl aus") beziehen sich auf diese drei Teile der VSCode-Oberfläche, sofern nicht anders angegeben.

Wenn du diesen Kurs selbstständig durcharbeitest, mache dich bitte mit den [Grundlagen der Umgebung](../../envsetup/01_setup.md) vertraut, um weitere Details zu erfahren.

### Versionsanforderungen

Dieses Training ist für Nextflow 25.10.2 oder neuer **mit AKTIVIERTEM v2-Syntax-Parser** konzipiert.
Wenn du eine lokale oder benutzerdefinierte Umgebung verwendest, stelle bitte sicher, dass du die korrekten Einstellungen verwendest, wie [hier](../../info/nxf_versions.md) dokumentiert.

## Bereit zum Arbeiten

Sobald dein Codespace läuft, musst du zwei Dinge tun, bevor du ins Training einsteigst: Setze dein Arbeitsverzeichnis für diesen spezifischen Kurs und wirf einen Blick auf die bereitgestellten Materialien.

### Das Arbeitsverzeichnis festlegen

Standardmäßig öffnet sich der Codespace mit dem Arbeitsverzeichnis im Stammverzeichnis aller Trainingskurse, aber für diesen Kurs arbeiten wir im Verzeichnis `nf4-science/genomics/`.

Wechsle jetzt das Verzeichnis, indem du diesen Befehl im Terminal ausführst:

```bash
cd nf4-science/genomics/
```

Du kannst VSCode so einstellen, dass es sich auf dieses Verzeichnis fokussiert, sodass nur die relevanten Dateien in der Dateisystem-Seitenleiste angezeigt werden:

```bash
code .
```

!!! tip "Tipp"

    Falls du aus irgendeinem Grund dieses Verzeichnis verlässt (z. B. wenn dein Codespace in den Ruhezustand geht), kannst du jederzeit den vollständigen Pfad verwenden, um dorthin zurückzukehren, vorausgesetzt, du arbeitest in der GitHub Codespaces Trainingsumgebung:

    ```bash
    cd /workspaces/training/nf4-science/genomics
    ```

Schauen wir uns nun den Inhalt an.

### Die bereitgestellten Materialien erkunden

Du kannst den Inhalt dieses Verzeichnisses über den Datei-Explorer auf der linken Seite des Trainings-Arbeitsbereichs erkunden.
Alternativ kannst du den Befehl `tree` verwenden.

Im gesamten Kurs verwenden wir die Ausgabe von `tree`, um die Verzeichnisstruktur und den Inhalt in lesbarer Form darzustellen, manchmal mit kleineren Anpassungen zur besseren Übersichtlichkeit.

Hier erzeugen wir ein Inhaltsverzeichnis bis zur zweiten Ebene:

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
Wir verwenden ausklappbare Abschnitte wie diesen, um erwartete Befehlsausgaben sowie Verzeichnis- und Dateiinhalte auf kompakte Weise darzustellen.

- **Die Datei `genomics.nf`** ist ein Workflow-Skript, das du im Laufe des Kurses aufbauen wirst.

- **Das Verzeichnis `modules`** enthält Modul-Gerüstdateien, die du während des Kurses ausfüllen wirst.

- **Die Datei `nextflow.config`** ist eine Konfigurationsdatei, die minimale Umgebungseigenschaften festlegt.
  Du kannst sie vorerst ignorieren.

- **Das Verzeichnis `data`** enthält Eingabedaten und zugehörige Ressourcen, die später im Kurs beschrieben werden.

- **Das Verzeichnis `solutions`** enthält fertige Moduldateien und eine Lösung für Teil 2, die als Ausgangspunkt für Teil 3 dienen kann.
  Sie sind als Referenz gedacht, um deine Arbeit zu überprüfen und Probleme zu beheben.

## Bereitschafts-Checkliste

Denkst du, du bist bereit loszulegen?

- [ ] Ich verstehe das Ziel dieses Kurses und seine Voraussetzungen
- [ ] Meine Umgebung läuft
- [ ] Ich habe mein Arbeitsverzeichnis entsprechend festgelegt

Wenn du alle Kästchen abhaken kannst, bist du startklar.

**Um mit [Teil 1: Methodenübersicht und manuelles Testen](./01_method.md) fortzufahren, klicke auf den Pfeil in der unteren rechten Ecke dieser Seite.**
