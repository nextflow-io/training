# Erste Schritte

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

## Trainingsumgebung starten

Um die von uns auf GitHub Codespaces bereitgestellte Umgebung zu nutzen, klicke auf den Button „Open in GitHub Codespaces" unten. Weitere Optionen findest du unter [Umgebungsoptionen](../envsetup/index.md).

Wir empfehlen, die Trainingsumgebung in einem neuen Browser-Tab oder -Fenster zu öffnen (je nach Gerät mit Rechtsklick, Strg+Klick oder Cmd+Klick), damit du weiterlesen kannst, während die Umgebung lädt.
Du solltest diese Anleitung parallel geöffnet halten, um den Kurs durchzuarbeiten.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

### Grundlagen der Umgebung

Diese Trainingsumgebung enthält die gesamte Software, den Code und die Daten, die du für den Kurs benötigst – du musst nichts selbst installieren.

Der Codespace ist mit einer VSCode-Oberfläche eingerichtet, die einen Datei-Explorer, einen Code-Editor und ein Terminal enthält.
Alle Anweisungen im Kurs (z. B. „Datei öffnen", „Code bearbeiten" oder „diesen Befehl ausführen") beziehen sich auf diese drei Teile der VSCode-Oberfläche, sofern nicht anders angegeben.

Wenn du diesen Kurs eigenständig durcharbeitest, mach dich bitte mit den [Grundlagen der Umgebung](../envsetup/01_setup.md) vertraut.

### Versionsanforderungen

Dieses Training funktioniert mit Nextflow 25.10.2 oder höher **mit dem v2-Syntax-Parser**, der ab Nextflow 26.04 standardmäßig verwendet wird.
In unserer Trainingsumgebung musst du nichts weiter tun: Sie läuft mit Nextflow 26.04.4 und dem v2-Parser. Wenn du eine lokale oder benutzerdefinierte Umgebung verwendest, lies die [Versionshinweise](../info/nxf_versions.md).

!!! warning "nf-core/demo erfordert Nextflow 25.10.4 oder höher"

    Die in Teil 1 verwendete `nf-core/demo`-Pipeline setzt ihre eigene Mindestversion von Nextflow voraus (`>=25.10.4`), die strenger ist als die allgemeine Trainingsanforderung von 25.10.2.
    Unsere Trainingsumgebung erfüllt diese Anforderung bereits. Wenn du eine lokale oder benutzerdefinierte Umgebung verwendest, stelle sicher, dass du Nextflow 25.10.4 oder höher nutzt.

Dieses Training erfordert außerdem **nf-core tools 4.0.2**.
Wenn du eine andere Version der nf-core-Tools verwendest, kann es schwierig sein, dem Kurs zu folgen.

Mit dem Befehl `nf-core --version` kannst du prüfen, welche Version in deiner Umgebung installiert ist.

!!! warning "Kompatibilität mit dem v2-Parser"

    Viele nf-core-Pipelines unterstützen den v2-Syntax-Parser noch nicht.
    Wenn du eine andere nf-core-Pipeline als die in diesem Kurs verwendeten ausführst und dabei Fehler auftreten, musst du möglicherweise zum v1-Parser wechseln: `export NXF_SYNTAX_PARSER=v1`.
    Weitere Details findest du in den [Versionshinweisen](../info/nxf_versions.md).

## Bereit zum Arbeiten

Sobald dein Codespace läuft, gibt es zwei Dinge, die du vor dem Einstieg ins Training erledigen musst: das Arbeitsverzeichnis für diesen Kurs festlegen und einen Blick auf die bereitgestellten Materialien werfen.

### Arbeitsverzeichnis festlegen

Standardmäßig öffnet sich der Codespace mit dem Arbeitsverzeichnis im Stammverzeichnis aller Trainingskurse. Für diesen Kurs arbeiten wir jedoch im Verzeichnis `nfcore-use/`.

Wechsle jetzt das Verzeichnis, indem du diesen Befehl im Terminal ausführst:

```bash
cd nfcore-use/
```

!!! tip "Tipp"

    Falls du dieses Verzeichnis aus irgendeinem Grund verlässt (z. B. weil dein Codespace in den Ruhezustand wechselt), kannst du jederzeit mit dem vollständigen Pfad dorthin zurückkehren – vorausgesetzt, du arbeitest in der GitHub Codespaces-Trainingsumgebung:

    ```bash
    cd /workspaces/training/nfcore-use
    ```

Erkunde als Nächstes den Inhalt dieses Verzeichnisses.

### Bereitgestellte Materialien erkunden

Du kannst den Inhalt dieses Verzeichnisses über den Datei-Explorer auf der linken Seite des Trainings-Workspaces erkunden.
Alternativ kannst du den Befehl `tree` verwenden.

```bash
tree .
```

??? abstract "Verzeichnisinhalt"

    ```console
    .
    ├── custom.config
    ├── laptop.config
    ├── malformed_samplesheet.csv
    └── my_params.yml
    ```

- **Die Datei `laptop.config`** ist eine Konfigurationsdatei, die wir in Abschnitt 4 verwenden, um die Ressourcennutzung beim lokalen Ausführen einer produktionsreifen Pipeline zu begrenzen.
  Du kannst sie bis dahin ignorieren.
- **Die Dateien `my_params.yml`, `malformed_samplesheet.csv` und `custom.config`** werden in Teil 2 verwendet, um das Festlegen von Parametern aus einer Datei, die Eingabevalidierung und prozessebene Konfigurationsüberschreibungen zu demonstrieren.
  Du kannst auch diese bis dahin ignorieren.

## Bereitschafts-Checkliste

Bereit zum Einstieg?

- [ ] Ich verstehe das Ziel dieses Kurses und seine Voraussetzungen
- [ ] Meine Umgebung läuft
- [ ] Ich verwende nf-core tools 4.0.2 (prüfen mit `nf-core --version`)
- [ ] Ich habe mein Arbeitsverzeichnis entsprechend gesetzt

Wenn du alle Punkte abhaken kannst, kann es losgehen.

**Um zu Teil 1 zu gelangen, klicke auf den Pfeil in der unteren rechten Ecke dieser Seite.**
