# Orientierung

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Die GitHub Codespaces-Umgebung enthält alle Software, Code und Daten, die für dieses Training notwendig sind, sodass du nichts selbst installieren musst.
Du benötigst jedoch einen (kostenlosen) Account zum Einloggen und solltest dir ein paar Minuten Zeit nehmen, um dich mit der Benutzeroberfläche vertraut zu machen.

Falls du das noch nicht getan hast, folge bitte [diesem Link](../../envsetup/), bevor du fortfährst.

## Bereitgestellte Materialien

Während dieses Trainings werden wir im Verzeichnis `side-quests/` arbeiten.
Dieses Verzeichnis enthält alle Code-Dateien, Testdaten und zusätzlichen Dateien, die du benötigst.

Erkunde gerne den Inhalt dieses Verzeichnisses; am einfachsten geht das mit dem Datei-Explorer auf der linken Seite des GitHub Codespaces-Arbeitsbereichs.
Alternativ kannst du den Befehl `tree` verwenden.
Im Verlauf des Kurses nutzen wir die Ausgabe von `tree`, um Verzeichnisstrukturen und -inhalte in einer lesbaren Form darzustellen, manchmal mit kleinen Anpassungen zur besseren Übersichtlichkeit.

Hier erstellen wir ein Inhaltsverzeichnis bis zur zweiten Ebene:

```bash
tree . -L 2
```

Wenn du dies innerhalb von `side-quests` ausführst, solltest du folgende Ausgabe sehen:

```console title="Verzeichnisinhalt"
.
├── metadata
├── nf-core
├── nf-test
├── solutions
├── splitting_and_grouping
└── workflows_of_workflows
```

**Hier ist eine Zusammenfassung dessen, was du zum Einstieg wissen solltest:**

- **Jedes Verzeichnis entspricht einer einzelnen Side Quest.**
  Deren Inhalte werden auf der entsprechenden Seite der Side Quest detailliert beschrieben.

- **Das Verzeichnis `solutions`** enthält die vollständigen Workflow- und/oder Modul-Skripte, die aus den verschiedenen Schritten jeder Side Quest resultieren.
  Sie dienen als Referenz, um deine Arbeit zu überprüfen und eventuell auftretende Probleme zu beheben.

!!!tip

    Falls du aus irgendeinem Grund dieses Verzeichnis verlässt, kannst du jederzeit mit diesem Befehl zurückkehren:

    ```bash
    cd /workspaces/training/side-quests
    ```

Um nun mit dem Kurs zu beginnen, klicke auf den Pfeil in der unteren rechten Ecke dieser Seite.
