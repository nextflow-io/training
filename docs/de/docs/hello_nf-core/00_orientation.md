# Erste Schritte

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

## Eine Trainingsumgebung starten

Um die vorgefertigte Umgebung zu nutzen, die wir auf GitHub Codespaces bereitstellen, klicke auf die Schaltfläche "Open in GitHub Codespaces" unten. Für andere Optionen siehe [Umgebungsoptionen](../envsetup/index.md).

Wir empfehlen, die Trainingsumgebung in einem neuen Browser-Tab oder -Fenster zu öffnen (verwende Rechtsklick, Strg-Klick oder Cmd-Klick, je nach deinem System), damit du weiterlesen kannst, während die Umgebung lädt.
Du musst diese Anleitung parallel geöffnet halten, um den Kurs durchzuarbeiten.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

### Grundlagen der Umgebung

Diese Trainingsumgebung enthält alle Software, Code und Daten, die zum Durcharbeiten des Trainingskurses notwendig sind, sodass du nichts selbst installieren musst.

Der Codespace ist mit einer VSCode-Oberfläche eingerichtet, die einen Dateisystem-Explorer, einen Code-Editor und eine Terminal-Shell umfasst.
Alle Anweisungen, die während des Kurses gegeben werden (z.B. 'öffne die Datei', 'bearbeite den Code' oder 'führe diesen Befehl aus'), beziehen sich auf diese drei Teile der VSCode-Oberfläche, sofern nicht anders angegeben.

Wenn du diesen Kurs selbstständig durcharbeitest, mache dich bitte mit den [Grundlagen der Umgebung](../envsetup/01_setup.md) vertraut, um weitere Details zu erfahren.

### Versionsanforderungen

Dieses Training ist für **Nextflow 25.10.2** oder neuer **mit DEAKTIVIERTEM v2-Syntax-Parser** konzipiert.

#### Wenn du unsere Trainingsumgebung verwendest:

Du MUSST den folgenden Befehl ausführen, bevor du fortfährst:

```bash
export NXF_SYNTAX_PARSER=v1
```

#### Wenn du eine lokale oder benutzerdefinierte Umgebung verwendest:

Bitte stelle sicher, dass du die korrekten Einstellungen verwendest, wie [hier](../info/nxf_versions.md) dokumentiert.

Das Training erfordert zusätzlich **nf-core tools 3.4.1**.
Wenn du eine andere Version der nf-core-Tools verwendest, könntest du Schwierigkeiten haben, dem Training zu folgen.

Du kannst überprüfen, welche Version in deiner Umgebung installiert ist, indem du den Befehl `nf-core --version` verwendest.

## Bereit zum Arbeiten

Sobald dein Codespace läuft, musst du zwei Dinge tun, bevor du ins Training einsteigst: Setze dein Arbeitsverzeichnis für diesen spezifischen Kurs und wirf einen Blick auf die bereitgestellten Materialien.

### Das Arbeitsverzeichnis setzen

Standardmäßig öffnet der Codespace mit dem Arbeitsverzeichnis im Stammverzeichnis aller Trainingskurse, aber für diesen Kurs werden wir im Verzeichnis `hello-nf-core/` arbeiten.

Wechsle jetzt das Verzeichnis, indem du diesen Befehl im Terminal ausführst:

```bash
cd hello-nf-core/
```

!!! tip "Tipp"

    Solltest du aus irgendeinem Grund aus diesem Verzeichnis herauswechseln (z.B. wenn dein Codespace in den Ruhemodus geht), kannst du immer den vollständigen Pfad verwenden, um dorthin zurückzukehren, vorausgesetzt du arbeitest in der GitHub Codespaces Trainingsumgebung:

    ```bash
    cd /workspaces/training/hello-nf-core
    ```

Lass uns nun einen Blick auf den Inhalt dieses Verzeichnisses werfen.

### Die bereitgestellten Materialien erkunden

Du kannst den Inhalt dieses Verzeichnisses mit dem Datei-Explorer auf der linken Seite des Trainings-Arbeitsbereichs erkunden.
Alternativ kannst du den Befehl `tree` verwenden.

Im Verlauf des Kurses verwenden wir die Ausgabe von `tree`, um die Verzeichnisstruktur und -inhalte in einer lesbaren Form darzustellen, manchmal mit kleinen Änderungen zur besseren Übersicht.

Hier erzeugen wir ein Inhaltsverzeichnis bis zur zweiten Ebene:

```bash
tree . -L 2
```

??? abstract "Verzeichnisinhalte"

    ```console
    .
    ├── greetings.csv
    ├── original-hello
    │   ├── hello.nf
    │   ├── modules
    │   └── nextflow.config
    └── solutions
        ├── composable-hello
        ├── core-hello-part2
        ├── core-hello-part3
        ├── core-hello-part4
        ├── core-hello-part5
        └── core-hello-start
    ```

Klicke auf das farbige Feld, um den Abschnitt zu erweitern und seinen Inhalt anzuzeigen.
Wir verwenden aufklappbare Abschnitte wie diesen, um die erwartete Befehlsausgabe auf prägnante Weise einzubinden.

- **Die Datei `greetings.csv`** ist eine CSV-Datei mit einigen minimalen Spaltendaten, die wir zu Testzwecken verwenden.

- **Das Verzeichnis `original-hello`** enthält eine Kopie des Quellcodes, der durch das Durcharbeiten der vollständigen Hello Nextflow Trainingsreihe entsteht (mit aktiviertem Docker).

- **Das Verzeichnis `solutions`** enthält die fertigen Workflow-Skripte, die aus jedem Schritt des Kurses resultieren.
  Sie sind als Referenz gedacht, um deine Arbeit zu überprüfen und eventuell auftretende Probleme zu beheben.

## Bereitschafts-Checkliste

Denkst du, du bist bereit loszulegen?

- [ ] Ich verstehe das Ziel dieses Kurses und seine Voraussetzungen
- [ ] Meine Umgebung ist eingerichtet und läuft
- [ ] Ich habe sichergestellt, dass der Syntax-Parser auf **v1** gesetzt ist
- [ ] Ich habe mein Arbeitsverzeichnis entsprechend gesetzt

Wenn du alle Kästchen abhaken kannst, kann es losgehen.

**Um mit Teil 1 fortzufahren, klicke auf den Pfeil unten rechts auf dieser Seite.**
