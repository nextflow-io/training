# Erste Schritte

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/PIjOdFaYwWA?si=gZxlXgkVxuLEzOsC" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } Sieh dir [die gesamte Playlist](https://youtube.com/playlist?list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&si=eF7cLR62goy-lc6n) auf dem Nextflow YouTube-Kanal an.

:green_book: Das Video-Transkript ist [hier](./transcripts/00_orientation.md) verfügbar.
///

!!! tip

    Die YouTube-Videos haben Superkräfte!

    - :fontawesome-solid-closed-captioning: Hochwertige (manuell kuratierte) Untertitel. Schalte sie mit dem :material-subtitles: Symbol ein
    - :material-bookmark: Video-Kapitel in der Timeline, die den Seitenüberschriften entsprechen.

## Eine Trainingsumgebung starten

Um die vorgefertigte Umgebung zu nutzen, die wir auf GitHub Codespaces bereitstellen, klicke auf die Schaltfläche „Open in GitHub Codespaces" unten. Für andere Optionen siehe [Umgebungsoptionen](../envsetup/index.md).

Wir empfehlen, die Trainingsumgebung in einem neuen Browser-Tab oder -Fenster zu öffnen (verwende Rechtsklick, Strg-Klick oder Cmd-Klick, je nach deiner Ausstattung), damit du weiterlesen kannst, während die Umgebung lädt.
Du musst diese Anleitung parallel geöffnet halten, um den Kurs durchzuarbeiten.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

### Grundlagen der Umgebung

Diese Trainingsumgebung enthält alle Software, den Code und die Daten, die zum Durcharbeiten des Trainingskurses erforderlich sind, sodass du nichts selbst installieren musst.

Der Codespace ist mit einer VSCode-Oberfläche eingerichtet, die einen Dateisystem-Explorer, einen Code-Editor und ein Terminal-Shell umfasst.
Alle Anweisungen während des Kurses (z. B. „öffne die Datei", „bearbeite den Code" oder „führe diesen Befehl aus") beziehen sich auf diese drei Teile der VSCode-Oberfläche, sofern nicht anders angegeben.

Wenn du diesen Kurs selbstständig durcharbeitest, mache dich bitte mit den [Grundlagen der Umgebung](../envsetup/01_setup.md) vertraut, um weitere Details zu erhalten.

### Versionsanforderungen

Dieses Training ist für Nextflow 25.10.2 oder später **mit AKTIVIERTEM v2-Syntax-Parser** konzipiert.
Wenn du eine lokale oder benutzerdefinierte Umgebung verwendest, stelle bitte sicher, dass du die korrekten Einstellungen verwendest, wie [hier](../info/nxf_versions.md) dokumentiert.

## Bereit zum Arbeiten

Sobald dein Codespace läuft, musst du zwei Dinge tun, bevor du ins Training eintauchst: Setze dein Arbeitsverzeichnis für diesen spezifischen Kurs und wirf einen Blick auf die bereitgestellten Materialien.

### Das Arbeitsverzeichnis festlegen

Standardmäßig öffnet sich der Codespace mit dem Arbeitsverzeichnis im Stammverzeichnis aller Trainingskurse, aber für diesen Kurs arbeiten wir im Verzeichnis `hello-nextflow/`.

Wechsle jetzt das Verzeichnis, indem du diesen Befehl im Terminal ausführst:

```bash
cd hello-nextflow/
```

Du kannst VSCode so einstellen, dass es sich auf dieses Verzeichnis konzentriert, sodass nur die relevanten Dateien in der Dateisystem-Explorer-Seitenleiste angezeigt werden:

```bash
code .
```

!!! tip

    Wenn du aus irgendeinem Grund dieses Verzeichnis verlässt (z. B. wenn dein Codespace in den Ruhezustand geht), kannst du jederzeit den vollständigen Pfad verwenden, um dorthin zurückzukehren, vorausgesetzt, du arbeitest in der Github Codespaces-Trainingsumgebung:

    ```bash
    cd /workspaces/training/hello-nextflow
    ```

Schauen wir uns nun den Inhalt an.

### Die bereitgestellten Materialien erkunden

Du kannst den Inhalt dieses Verzeichnisses mit dem Datei-Explorer auf der linken Seite des Trainingsarbeitsbereichs erkunden.
Alternativ kannst du den Befehl `tree` verwenden.

Im gesamten Kurs verwenden wir die Ausgabe von `tree`, um die Verzeichnisstruktur und den Inhalt in lesbarer Form darzustellen, manchmal mit kleineren Änderungen zur besseren Übersichtlichkeit.

Hier erzeugen wir ein Inhaltsverzeichnis bis zur zweiten Ebene:

```bash
tree . -L 2
```

??? abstract "Verzeichnisinhalt"

    ```console
    .
    ├── data
    │   └── greetings.csv
    ├── hello-channels.nf
    ├── hello-config.nf
    ├── hello-containers.nf
    ├── hello-modules.nf
    ├── hello-workflow.nf
    ├── hello-world.nf
    ├── nextflow.config
    ├── solutions
    │   ├── 1-hello-world
    │   ├── 2-hello-channels
    │   ├── 3-hello-workflow
    │   ├── 4-hello-modules
    │   ├── 5-hello-containers
    │   └── 6-hello-config
    ├── test-params.json
    └── test-params.yaml
    ```

Klicke auf das farbige Feld, um den Abschnitt zu erweitern und seinen Inhalt anzuzeigen.
Wir verwenden solche ausklappbaren Abschnitte, um die erwartete Befehlsausgabe auf kompakte Weise einzubeziehen.

- **Die `.nf`-Dateien** sind Workflow-Skripte, die nach dem Teil des Kurses benannt sind, in dem sie verwendet werden.

- **Die Datei `nextflow.config`** ist eine Konfigurationsdatei, die minimale Umgebungseigenschaften festlegt.
  Du kannst sie vorerst ignorieren.

- **Die Datei `greetings.csv`** unter `data/` enthält Eingabedaten, die wir in den meisten Teilen des Kurses verwenden werden. Sie wird in Teil 2 (Kanäle) beschrieben, wenn wir sie zum ersten Mal einführen.

- **Die `test-params.*`-Dateien** sind Konfigurationsdateien, die wir in Teil 6 (Konfiguration) verwenden werden. Du kannst sie vorerst ignorieren.

- **Das Verzeichnis `solutions`** enthält die fertigen Workflow-Skripte, die aus jedem Schritt des Kurses resultieren.
  Sie sind als Referenz gedacht, um deine Arbeit zu überprüfen und eventuelle Probleme zu beheben.

## Bereitschafts-Checkliste

Denkst du, du bist bereit einzutauchen?

- [ ] Ich verstehe das Ziel dieses Kurses und seine Voraussetzungen
- [ ] Meine Umgebung ist eingerichtet und läuft
- [ ] Ich habe mein Arbeitsverzeichnis entsprechend festgelegt

Wenn du alle Kästchen abhaken kannst, bist du startklar.

**Um mit [Teil 1: Hello World](./01_hello_world.md) fortzufahren, klicke auf den Pfeil in der unteren rechten Ecke dieser Seite.**
