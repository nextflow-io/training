# Erste Schritte

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

## Eine Trainingsumgebung starten

Um die vorgefertigte Umgebung zu verwenden, die wir auf GitHub Codespaces bereitstellen, klicke auf die Schaltfläche "Open in GitHub Codespaces" unten. Für andere Optionen siehe [Umgebungsoptionen](../envsetup/index.md).

Wir empfehlen, die Trainingsumgebung in einem neuen Browser-Tab oder -Fenster zu öffnen (verwende Rechtsklick, Strg-Klick oder Cmd-Klick je nach Gerät), damit du weiterlesen kannst, während die Umgebung lädt.
Du musst diese Anleitung parallel geöffnet halten, um den Kurs durchzuarbeiten.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

### Grundlagen der Umgebung

Diese Trainingsumgebung enthält alle Software, den Code und die Daten, die zum Durcharbeiten des Kurses notwendig sind, sodass du nichts selbst installieren musst.

Der Codespace ist mit einer VSCode-Oberfläche eingerichtet, die einen Dateisystem-Explorer, einen Code-Editor und eine Terminal-Shell enthält.
Alle Anweisungen während des Kurses (z.B. 'öffne die Datei', 'bearbeite den Code' oder 'führe diesen Befehl aus') beziehen sich auf diese drei Teile der VSCode-Oberfläche, sofern nicht anders angegeben.

Wenn du diesen Kurs selbstständig durcharbeitest, mache dich bitte mit den [Grundlagen der Umgebung](../envsetup/01_setup.md) für weitere Details vertraut.

### Versionsanforderungen

Dieses Training ist für Nextflow 25.10.2 oder später **mit aktiviertem v2-Syntax-Parser** konzipiert.
Wenn du eine lokale oder benutzerdefinierte Umgebung verwendest, stelle bitte sicher, dass du die richtigen Einstellungen wie [hier](../info/nxf_versions.md) dokumentiert verwendest.

## Bereit zum Arbeiten

Sobald dein Codespace läuft, musst du zwei Dinge tun, bevor du in das Training eintauchst: dein Arbeitsverzeichnis für diesen speziellen Kurs festlegen und dir die bereitgestellten Materialien ansehen.

### Das Arbeitsverzeichnis festlegen

Standardmäßig öffnet sich der Codespace mit dem Arbeitsverzeichnis am Stamm aller Trainingskurse, aber für diesen Kurs werden wir im Verzeichnis `nextflow-run/` arbeiten.

Wechsle jetzt das Verzeichnis, indem du diesen Befehl im Terminal ausführst:

```bash
cd nextflow-run/
```

Du kannst VSCode so einstellen, dass es sich auf dieses Verzeichnis konzentriert, sodass nur die relevanten Dateien in der Datei-Explorer-Seitenleiste angezeigt werden:

```bash
code .
```

!!! tip "Tipp"

    Wenn du aus irgendeinem Grund dieses Verzeichnis verlässt (z.B. dein Codespace geht in den Schlafmodus), kannst du immer den vollständigen Pfad verwenden, um dorthin zurückzukehren, vorausgesetzt du arbeitest in der GitHub Codespaces-Trainingsumgebung:

    ```bash
    cd /workspaces/training/nextflow-run
    ```

Schauen wir uns jetzt den Inhalt an.

### Die bereitgestellten Materialien erkunden

Du kannst den Inhalt dieses Verzeichnisses erkunden, indem du den Datei-Explorer auf der linken Seite des Trainingsarbeitsbereichs verwendest.
Alternativ kannst du den Befehl `tree` verwenden.

Während des gesamten Kurses verwenden wir die Ausgabe von `tree`, um Verzeichnisstruktur und -inhalt in lesbarer Form darzustellen, manchmal mit kleineren Anpassungen zur Klarheit.

Hier generieren wir ein Inhaltsverzeichnis bis zur zweiten Ebene:

```bash
tree . -L 2
```

??? abstract "Verzeichnisinhalt"

    ```console
    .
    ├── 1-hello.nf
    ├── 2a-inputs.nf
    ├── 2b-multistep.nf
    ├── 2c-modules.nf
    ├── 2d-container.nf
    ├── 3-main.nf
    ├── data
    │   └── greetings.csv
    ├── modules
    │   ├── collectGreetings.nf
    │   ├── convertToUpper.nf
    │   ├── cowpy.nf
    │   └── sayHello.nf
    ├── nextflow.config
    ├── solutions
    │   ├── 3-main.nf
    │   ├── modules
    │   └── nextflow.config
    ├── test-params.json
    └── test-params.yaml
    ```

Klicke auf das farbige Feld, um den Abschnitt zu erweitern und seinen Inhalt anzuzeigen.
Wir verwenden solche einklappbaren Abschnitte, um erwartete Befehlsausgaben sowie Verzeichnis- und Dateiinhalte kompakt darzustellen.

- **Die `.nf`-Dateien** sind Workflow-Scripts, die basierend darauf nummeriert sind, in welchem Teil des Kurses sie verwendet werden.

- **Die Datei `nextflow.config`** ist eine Konfigurationsdatei, die minimale Umgebungseigenschaften festlegt.
  Du kannst sie vorerst ignorieren.

- **Die Datei `greetings.csv`** unter `data/` enthält Eingabedaten, die wir im größten Teil des Kurses verwenden werden. Sie wird in Teil 2 (Pipelines ausführen) beschrieben, wenn wir sie zum ersten Mal einführen.

- **Die `test-params.*`-Dateien** sind Konfigurationsdateien, die wir in Teil 3 (Konfiguration) verwenden werden. Du kannst sie vorerst ignorieren.

- **Das `solutions`-Verzeichnis** enthält den Endzustand des Workflows und seiner Hilfsdateien (Config und Module), die sich aus dem Abschluss des Kurses ergeben.
  Sie sind als Referenz gedacht, um deine Arbeit zu überprüfen und eventuelle Probleme zu beheben.

## Bereitschafts-Checkliste

Denkst du, du bist bereit einzutauchen?

- [ ] Ich verstehe das Ziel dieses Kurses und seine Voraussetzungen
- [ ] Meine Umgebung ist eingerichtet und läuft
- [ ] Ich habe mein Arbeitsverzeichnis entsprechend festgelegt

Wenn du alle Kästchen abhaken kannst, bist du startklar.

**Um zu [Teil 1: Grundlegende Operationen](./01_basics.md) fortzufahren, klicke auf den Pfeil in der unteren rechten Ecke dieser Seite.**
