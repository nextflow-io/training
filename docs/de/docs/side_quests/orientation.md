# Orientierung

Die GitHub Codespaces-Umgebung enthält alle Software, Code und Daten, die für diesen Trainingskurs notwendig sind, sodass du nichts selbst installieren musst.
Du benötigst jedoch ein (kostenloses) Konto zum Einloggen und solltest dir ein paar Minuten Zeit nehmen, um dich mit der Oberfläche vertraut zu machen.

Falls du dies noch nicht getan hast, folge bitte [diesem Link](../../envsetup/), bevor du weitermachst.

## Bereitgestellte Materialien

Während dieses Trainingskurses arbeiten wir im Verzeichnis `side-quests/`.
Dieses Verzeichnis enthält alle Code-Dateien, Testdaten und zusätzlichen Dateien, die du benötigen wirst.

Erkunde gerne den Inhalt dieses Verzeichnisses. Am einfachsten geht das über den Datei-Explorer auf der linken Seite des GitHub Codespaces-Arbeitsbereichs.
Alternativ kannst du den Befehl `tree` verwenden.
Im Verlauf des Kurses nutzen wir die Ausgabe von `tree`, um die Verzeichnisstruktur und -inhalte in lesbarer Form darzustellen, manchmal mit kleinen Anpassungen zur besseren Übersichtlichkeit.

Hier erzeugen wir ein Inhaltsverzeichnis bis zur zweiten Ebene:

```bash
tree . -L 2
```

Wenn du dies innerhalb von `side-quests` ausführst, solltest du folgende Ausgabe sehen:

```console title="Directory contents"
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
  Deren Inhalte werden auf der Seite der jeweiligen Side Quest detailliert beschrieben.

- **Das Verzeichnis `solutions`** enthält die fertigen Workflow- und/oder Modul-Skripte, die sich aus den verschiedenen Schritten jeder Side Quest ergeben.
  Sie dienen als Referenz, um deine Arbeit zu überprüfen und eventuelle Probleme zu beheben.

!!!tip

    Falls du aus irgendeinem Grund dieses Verzeichnis verlässt, kannst du jederzeit mit diesem Befehl zurückkehren:

    ```bash
    cd /workspaces/training/side-quests
    ```

Um nun mit dem Kurs zu beginnen, klicke auf den Pfeil unten rechts auf dieser Seite.
