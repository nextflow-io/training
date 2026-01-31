# Orientierung

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Die Trainingsumgebung enthält alle Software, den Code und die Daten, die für diesen Trainingskurs notwendig sind, sodass du nichts selbst installieren musst.
Du benötigst jedoch einen (kostenlosen) Account, um dich anzumelden, und solltest dir ein paar Minuten Zeit nehmen, um dich mit der Oberfläche vertraut zu machen.

Falls du dies noch nicht getan hast, folge bitte [diesem Link](../../../envsetup/), bevor du weiter fortfährst.

## Bereitgestellte Materialien

Während dieses Trainingskurses arbeiten wir im Verzeichnis `nf4-science/genomics/`, in das du wechseln musst, wenn du den Trainings-Workspace öffnest.
Dieses Verzeichnis enthält alle Code-Dateien, Testdaten und Zusatzdateien, die du benötigst.

Du kannst gerne den Inhalt dieses Verzeichnisses erkunden; am einfachsten geht das mit dem Datei-Explorer auf der linken Seite des Trainings-Workspace in der VSCode-Oberfläche.
Alternativ kannst du den Befehl `tree` verwenden.
Im Laufe des Kurses nutzen wir die Ausgabe von `tree`, um die Verzeichnisstruktur und den Inhalt in lesbarer Form darzustellen, manchmal mit kleinen Änderungen zur besseren Übersicht.

Hier erstellen wir ein Inhaltsverzeichnis bis zur zweiten Ebene:

```bash
tree . -L 2
```

Wenn du dies innerhalb von `nf4-science/genomics` ausführst, solltest du folgende Ausgabe sehen:

```console title="Verzeichnisinhalt"

.
├── data
│   ├── bam
│   ├── ref
│   ├── sample_bams.txt
│   └── samplesheet.csv
├── genomics-1.nf
├── genomics-2.nf
├── genomics-3.nf
├── genomics-4.nf
├── nextflow.config
└── solutions
    ├── modules
    ├── nf-test.config
    └── tests

6 directories, 8 files

```

!!!note "Hinweis"

    Mach dir keine Sorgen, falls das viel erscheint; wir werden die relevanten Teile bei jedem Schritt des Kurses durchgehen.
    Dies soll dir nur einen Überblick geben.

**Hier ist eine Zusammenfassung dessen, was du für den Einstieg wissen solltest:**

- **Die `.nf`-Dateien** sind Workflow-Skripte, die nach dem entsprechenden Teil des Kurses benannt sind, in dem sie verwendet werden.

- **Die Datei `nextflow.config`** ist eine Konfigurationsdatei, die minimale Umgebungseigenschaften festlegt.
  Du kannst sie vorerst ignorieren.

- **Das Verzeichnis `data`** enthält Eingabedaten und zugehörige Ressourcen, die später im Kurs beschrieben werden.

- **Das Verzeichnis `solutions`** enthält Moduldateien und Testkonfigurationen, die aus den Teilen 3 und 4 des Kurses resultieren.
  Sie sind als Referenz gedacht, um deine Arbeit zu überprüfen und eventuelle Probleme zu beheben.

!!!tip "Tipp"

    Falls du aus irgendeinem Grund aus diesem Verzeichnis herausnavigierst, kannst du jederzeit mit diesem Befehl zurückkehren:

    ```bash
    cd /workspaces/training/nf4-science/genomics
    ```

Um den Kurs zu beginnen, klicke auf den Pfeil in der unteren rechten Ecke dieser Seite.
