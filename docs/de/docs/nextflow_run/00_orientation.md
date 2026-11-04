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

Dieser Kurs erfordert Nextflow 25.10.2 oder später, mit aktiviertem v2-Syntax-Parser (Standard ab 25.10+).
Wenn du eine lokale oder benutzerdefinierte Umgebung verwendest, stelle bitte sicher, dass du die richtigen Einstellungen wie [hier](../info/nxf_versions.md) dokumentiert verwendest.

## Bereit zum Arbeiten

Sobald dein Codespace läuft, musst du zwei Dinge tun, bevor du in das Training eintauchst: dein Arbeitsverzeichnis festlegen und dir die bereitgestellten Materialien ansehen.

### Das Arbeitsverzeichnis festlegen

Standardmäßig öffnet sich der Codespace am Stamm aller Trainingskurse.
Wechsle für diesen Kurs in das Verzeichnis `nextflow-run/`:

```bash
cd nextflow-run/
```

Dann stelle VSCode so ein, dass es sich auf dieses Verzeichnis konzentriert, sodass nur die relevanten Dateien in der Datei-Explorer-Seitenleiste angezeigt werden:

```bash
code .
```

!!! tip "Tipp"

    Wenn du aus irgendeinem Grund dieses Verzeichnis verlässt (z.B. dein Codespace geht in den Schlafmodus), kannst du immer den vollständigen Pfad verwenden, um dorthin zurückzukehren, vorausgesetzt du arbeitest in der GitHub Codespaces-Trainingsumgebung:

    ```bash
    cd /workspaces/training/nextflow-run
    ```

### Die bereitgestellten Materialien erkunden

Du kannst die Kursmaterialien mit dem Datei-Explorer auf der linken Seite oder mit dem Befehl `tree` erkunden.
Führe folgenden Befehl im Terminal aus, um die vollständige Struktur zu sehen:

```bash
tree . -L 2
```

??? abstract "Verzeichnisinhalt"

    ```console
    .
    ├── 1-hello.nf
    ├── 2-inputs.nf
    ├── data
    │   ├── greetings-extended.csv
    │   └── greetings.csv
    ├── main.nf
    ├── modules
    │   ├── collectGreetings.nf
    │   ├── convertToUpper.nf
    │   ├── cowpy.nf
    │   └── sayHello.nf
    └── nextflow.config
    ```

Die **`.nf`-Dateien** sind Workflow-Scripts mit zunehmender Komplexität, die in dieser Reihenfolge im Kurs verwendet werden.

Das **`data/`-Verzeichnis** enthält die CSV-Eingabedateien, die wir ab Abschnitt 2 verwenden werden.

Das **`modules/`-Verzeichnis** enthält die Prozessdefinitionen, die von `main.nf` verwendet werden.

Die **`nextflow.config`-Datei** ist eine Konfigurationsdatei, die minimale Umgebungseigenschaften festlegt. Du kannst sie vorerst ignorieren; wir gehen in Abschnitt 4 darauf ein.

## Bereitschafts-Checkliste

Denkst du, du bist bereit einzutauchen?

- [ ] Ich verstehe das Ziel dieses Kurses und seine Voraussetzungen
- [ ] Meine Umgebung ist eingerichtet und läuft
- [ ] Ich habe mein Arbeitsverzeichnis entsprechend festgelegt

Wenn du alle Kästchen abhaken kannst, bist du startklar.

**Um zu [Teil 1: Nextflow ausführen](./01_run_nextflow.md) fortzufahren, klicke auf den Pfeil in der unteren rechten Ecke dieser Seite.**
