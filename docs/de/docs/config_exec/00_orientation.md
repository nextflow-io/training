# Erste Schritte

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

## Trainingsumgebung starten

Um die von uns auf GitHub Codespaces bereitgestellte Umgebung zu nutzen, klicke auf die Schaltfläche „Open in GitHub Codespaces" unten. Weitere Optionen findest du unter [Umgebungsoptionen](../envsetup/index.md).

Wir empfehlen, die Trainingsumgebung in einem neuen Browser-Tab oder -Fenster zu öffnen (je nach Gerät per Rechtsklick, Strg+Klick oder Cmd+Klick), damit du weiterlesen kannst, während die Umgebung lädt.
Du solltest diese Anleitung parallel geöffnet halten, um den Kurs durchzuarbeiten.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

### Grundlagen der Umgebung

Diese Trainingsumgebung enthält die gesamte Software, den Code und die Daten, die du für den Kurs benötigst – du musst also nichts selbst installieren.

Der Codespace ist mit einer VSCode-Oberfläche eingerichtet, die einen Datei-Explorer, einen Code-Editor und ein Terminal enthält.
Alle Anweisungen im Kurs (z. B. „Datei öffnen", „Code bearbeiten" oder „diesen Befehl ausführen") beziehen sich auf diese drei Teile der VSCode-Oberfläche, sofern nicht anders angegeben.

Wenn du diesen Kurs eigenständig durcharbeitest, mach dich bitte mit den [Grundlagen der Umgebung](../envsetup/01_setup.md) für weitere Details vertraut.

### Versionsanforderungen

Dieser Kurs erfordert Nextflow 25.10.2 oder höher, mit aktiviertem v2-Syntax-Parser (Standard ab 25.10+).
Wenn du eine lokale oder benutzerdefinierte Umgebung verwendest, stelle sicher, dass du die richtigen Einstellungen verwendest, wie [hier](../info/nxf_versions.md) dokumentiert.

## Bereit zum Arbeiten

Sobald dein Codespace läuft, gibt es zwei Dinge zu tun, bevor du einsteigst: das Arbeitsverzeichnis festlegen und einen Blick auf die bereitgestellten Materialien werfen.

### Arbeitsverzeichnis festlegen

Standardmäßig öffnet sich der Codespace im Stammverzeichnis aller Trainingskurse.
Wechsle für diesen Kurs in das Verzeichnis `config-exec/`:

```bash
cd config-exec/
```

Dann weise VSCode an, sich auf dieses Verzeichnis zu fokussieren, damit nur die relevanten Dateien in der Datei-Explorer-Seitenleiste angezeigt werden:

```bash
code .
```

!!! tip "Tipp"

    Falls du dieses Verzeichnis aus irgendeinem Grund verlässt (z. B. wenn dein Codespace in den Ruhezustand wechselt), kannst du jederzeit den vollständigen Pfad verwenden, um dorthin zurückzukehren – vorausgesetzt, du arbeitest in der GitHub Codespaces-Trainingsumgebung:

    ```bash
    cd /workspaces/training/config-exec
    ```

### Bereitgestellte Materialien erkunden

Du kannst die Kursmaterialien über den Datei-Explorer auf der linken Seite oder mit dem Befehl `tree` erkunden.
Führe folgenden Befehl im Terminal aus, um die vollständige Struktur anzuzeigen:

```bash
tree . -L 2
```

??? abstract "Verzeichnisinhalt"

    ```console
    .
    ├── data
    │   └── greetings.csv
    ├── main.nf
    ├── modules
    │   ├── collectGreetings.nf
    │   ├── convertToUpper.nf
    │   ├── cowpy.nf
    │   └── sayHello.nf
    └── nextflow.config
    ```

Die Dateien **`main.nf`** und **`modules/`** sind dieselbe mehrstufige Pipeline aus [Nextflow Run](../nextflow_run/index.md), und die Datei **`nextflow.config`** ist dieselbe Konfiguration, die du dort bereits gesehen hast.
Du wirst beide im Laufe dieser Übungen erweitern.

Das Verzeichnis **`data/`** enthält die CSV-Eingabedatei, aus der die Pipeline liest.

## Bereitschafts-Checkliste

Bereit zum Einsteigen?

- [ ] Ich verstehe das Ziel dieses Kurses und seine Voraussetzungen
- [ ] Meine Umgebung ist eingerichtet und läuft
- [ ] Ich habe mein Arbeitsverzeichnis entsprechend festgelegt

Wenn du alle Punkte abhaken kannst, kann es losgehen.

**Um zu [Teil 1: An deine Rechenumgebung anpassen](./01_packaging_and_execution.md) zu gelangen, klicke auf den Pfeil in der unteren rechten Ecke dieser Seite.**
