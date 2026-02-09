# Orientierung

Diese Orientierung setzt voraus, dass du die Trainingsumgebung bereits durch Klicken auf den Button „Open in GitHub Codespaces" geöffnet hast.
Falls nicht, tue dies bitte jetzt, idealerweise in einem zweiten Browserfenster oder -tab, damit du auf diese Anleitung zurückgreifen kannst.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://github.com/codespaces/new?hide_repo_select=true&ref=master&repo=290790519&skip_quickstart=true&machine=premiumLinux&devcontainer_path=.devcontainer%2Fdevcontainer.json)

!!!warning "Anforderung an die Maschinengröße"

    Stelle sicher, dass du beim Erstellen deines Codespace für diesen Trainingskurs eine **8-Kern-Maschine** auswählst. Die Bioimaging-Workflows benötigen zusätzliche Rechenressourcen.

## GitHub Codespaces

Die GitHub Codespaces-Umgebung enthält alle Software, den Code und die Daten, die für diesen Trainingskurs erforderlich sind, sodass du nichts selbst installieren musst.
Du benötigst jedoch ein (kostenloses) GitHub-Konto zum Anmelden. Falls du mit der Oberfläche nicht vertraut bist, solltest du dir ein paar Minuten Zeit nehmen, um dich damit vertraut zu machen, indem du den Mini-Kurs [GitHub Codespaces Orientation](../../envsetup/index.md) absolvierst.

## Docker-Images vorab herunterladen

Sobald du deinen Codespace geöffnet hast, lass uns alle Docker-Images herunterladen, die wir für diesen Trainingskurs benötigen.
Das spart später Zeit und sorgt für eine reibungslose Ausführung der Workflows.

Öffne einen neuen Terminal-Tab und führe folgenden Befehl aus:

```bash
nextflow run nf-core/molkart -profile docker,test -stub -resume --outdir results
```

Dieser Befehl lädt alle notwendigen Docker-Images im Hintergrund herunter.
Du kannst mit dem Rest der Orientierung fortfahren, während dies läuft.

!!!tip

    Das Flag `-stub` ermöglicht es der Pipeline, schnell zu laufen, ohne echte Daten zu verarbeiten – perfekt zum Herunterladen von Images. Du kannst den Fortschritt im Terminal-Tab verfolgen.

## Arbeitsverzeichnis

Während dieses Trainingskurses arbeiten wir im Verzeichnis `nf4-science/imaging/`.

Wechsle jetzt das Verzeichnis, indem du diesen Befehl im Terminal ausführst:

```bash
cd nf4-science/imaging/
```

!!!tip

    Falls du aus irgendeinem Grund dieses Verzeichnis verlässt, kannst du jederzeit mit dem vollständigen Pfad dorthin zurückkehren, vorausgesetzt, du arbeitest in der GitHub Codespaces-Trainingsumgebung:

    ```bash
    cd /workspaces/training/nf4-science/imaging
    ```

**Um nun mit dem Kurs zu beginnen, klicke auf den Pfeil unten rechts auf dieser Seite.**
