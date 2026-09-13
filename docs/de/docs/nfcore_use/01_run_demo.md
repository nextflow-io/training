# Teil 1: Eine Demo-Pipeline ausführen

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

In diesem ersten Teil des Kurses „Use nf-core" zeigen wir dir, wie du eine nf-core-Pipeline findest und mit ihrem integrierten Test-Profil ausprobierst.

Wir verwenden eine Pipeline namens nf-core/demo, die vom nf-core-Projekt als Teil seiner Sammlung von Pipelines für Demonstrations- und Trainingszwecke gepflegt wird.

Stelle sicher, dass dein Arbeitsverzeichnis auf `nfcore-use/` gesetzt ist, wie auf der Seite [Erste Schritte](./00_orientation.md) beschrieben.

---

## 1. Die nf-core/demo-Pipeline finden und herunterladen

Beginnen wir damit, die nf-core/demo-Pipeline auf der Projektwebsite unter [nf-co.re](https://nf-co.re) zu finden. Diese Website bündelt alle Informationen: allgemeine Dokumentation und Hilfeartikel, Dokumentation für jede Pipeline, Blogbeiträge, Veranstaltungsankündigungen und vieles mehr.

### 1.1. Die Pipeline auf der Website finden

Öffne in deinem Webbrowser [https://nf-co.re/pipelines/](https://nf-co.re/pipelines/) und gib `demo` in die Suchleiste ein.

![Suchergebnisse](./img/search-results.png)

Klicke auf den Pipeline-Namen `demo`, um zur Dokumentationsseite der Pipeline zu gelangen.

Jede veröffentlichte Pipeline hat eine eigene Seite mit folgenden Dokumentationsabschnitten:

- **Introduction:** Eine Einführung und Übersicht der Pipeline
- **Usage:** Beschreibungen zur Ausführung der Pipeline
- **Parameters:** Gruppierte Pipeline-Parameter mit Beschreibungen
- **Output:** Beschreibungen und Beispiele der erwarteten Ausgabedateien
- **Results:** Beispielausgaben aus dem vollständigen Test-Datensatz
- **Releases & Statistics:** Versionsverlauf und Statistiken der Pipeline

Bevor du eine neue Pipeline verwendest, solltest du die Dokumentation sorgfältig lesen, um zu verstehen, was sie tut und wie sie konfiguriert werden muss.

Schau jetzt nach und versuche herauszufinden:

- Welche Tools die Pipeline ausführt (Tab: `Introduction`)
- Welche Eingaben und Parameter die Pipeline akzeptiert oder erfordert (Tab: `Parameters`)
- Welche Ausgaben die Pipeline erzeugt (Tab: `Output`)

#### 1.1.1. Pipeline-Übersicht

Der Tab `Introduction` bietet eine Übersicht der Pipeline, einschließlich einer visuellen Darstellung (auch „Subway Map" genannt) und einer Liste der Tools, die als Teil der Pipeline ausgeführt werden.

![Pipeline Subway Map](./img/nf-core-demo-subway-cropped.png)

1. Read QC ([FASTQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
2. Adapter- und Qualitäts-Trimming ([SEQTK_TRIM](https://github.com/lh3/seqtk))
3. QC-Darstellung für Rohdaten ([MULTIQC](http://multiqc.info/))
4. Eine humorvolle Textnachricht von einer Kuh generieren ([COWPY](https://github.com/jeffbuttars/cowpy))

#### 1.1.2. Beispiel-Befehlszeile

Die Dokumentation enthält außerdem eine Beispiel-Eingabedatei (weiter unten besprochen) und eine Beispiel-Befehlszeile.

```bash
nextflow run nf-core/demo \
  -profile <docker/singularity/.../institute> \
  --input samplesheet.csv \
  --outdir <OUTDIR>
```

Du wirst bemerken, dass der Beispielbefehl KEINE Workflow-Datei angibt, sondern nur den Verweis auf das Pipeline-Repository `nf-core/demo`.

Wenn Nextflow so aufgerufen wird, geht es davon aus, dass der Code auf eine bestimmte Weise organisiert ist.
Lass uns den Code herunterladen, damit wir diese Struktur untersuchen können.

### 1.2. Den Pipeline-Code herunterladen

Nachdem wir festgestellt haben, dass die Pipeline für unsere Zwecke geeignet erscheint, wollen wir sie ausprobieren.
Glücklicherweise macht es Nextflow einfach, Pipelines aus korrekt formatierten Repositories herunterzuladen, ohne etwas manuell herunterladen zu müssen.

#### 1.2.1. `nextflow pull` verwenden

Kehren wir zum Terminal zurück und führen Folgendes aus:

```bash
nextflow pull nf-core/demo
```

??? success "Befehlsausgabe"

    ```console
    Checking nf-core/demo ...
     downloaded from https://github.com/nf-core/demo.git - revision: 32893afef8 [master]
    ```

Nextflow führt einen `pull` des Pipeline-Codes durch, d. h. es lädt das vollständige Repository auf dein lokales Laufwerk herunter.

Das funktioniert mit jeder Nextflow-Pipeline, die entsprechend auf GitHub eingerichtet ist – nicht nur mit nf-core-Pipelines.
nf-core ist jedoch die größte Open-Source-Sammlung von Nextflow-Pipelines.

#### 1.2.2. `nextflow list` verwenden

Du kannst Nextflow dazu bringen, dir eine Liste der auf diese Weise heruntergeladenen Pipelines anzuzeigen:

```bash
nextflow list
```

??? success "Befehlsausgabe"

    ```console
    nf-core/demo
    ```

Du kannst weitere Pipelines herunterladen, um zu sehen, wie sie aufgelistet werden, wenn du mehr als eine hast.

#### 1.2.3. Herausfinden, wo die Pipeline gespeichert wurde

Du wirst bemerken, dass sich die Dateien nicht in deinem aktuellen Arbeitsverzeichnis befinden.
Standardmäßig speichert Nextflow heruntergeladene Pipelines unter `$NXF_HOME/assets`.

Um herauszufinden, wo eine bestimmte Pipeline gespeichert ist, frage Nextflow direkt:

```bash
nextflow info nf-core/demo
```

??? success "Befehlsausgabe"

    ```console
     project name: nf-core/demo
     repository  : https://github.com/nf-core/demo
     local path  : /workspaces/.nextflow/assets/.repos/nf-core/demo
     main script : main.nf
     description : An nf-core demo pipeline
     revisions   :
       TEMPLATE
       bumper
       dev
       fix-nxfversion
       manually-merge-3_0_2
     > master (default)
       nf-core-template-merge-2.13.2.dev0
       nf-core-template-merge-2.14.0
       nf-core-template-merge-2.14.1
       nf-core-template-merge-3.0.0
       nf-core-template-merge-3.0.1
       nf-core-template-merge-3.0.2
       nf-core-template-merge-3.1.0
       nf-core-template-merge-3.1.2
       nf-core-template-merge-3.2.0
       nf-core-template-merge-3.2.1
       nf-core-template-merge-3.3.1
       nf-core-template-merge-3.3.2
       nf-core-template-merge-4.0.0
       nf-core-template-merge-4.0.3
       nf-core-template-merge-4.1.0
       nf-core-template-merge-4.1.0-2
       patch
       1.0.0 [t]
       1.0.1 [t]
       1.0.2 [t]
       1.1.0 [t]
     > 1.2.0 [t]
    ```

!!! info "Info"

    Der vollständige Pfad kann auf deinem System abweichen, wenn du nicht unsere Trainingsumgebung verwendest.

Nextflow speichert den heruntergeladenen Quellcode bewusst „aus dem Weg", nach dem Prinzip, dass diese Pipelines eher wie Bibliotheken verwendet werden sollten als Code, mit dem du direkt interagierst.

Im Hintergrund speichert Nextflow jede heruntergeladene Pipeline als git-Repository unter `$NXF_HOME/assets/.repos/` und checkt den Code für jede Revision in ein Unterverzeichnis `clones/<commit>/` aus.
Da `.repos` ein verstecktes Verzeichnis ist, erscheint ein einfaches `tree -L 2 $NXF_HOME/assets/` leer.

#### 1.2.4. Einen Symlink erstellen, um einfach auf den Quellcode zuzugreifen

Wir werden den Code nicht im Detail durchgehen, aber lass uns kurz einen Blick darauf werfen, um ein Gefühl für die Gesamtorganisation zu bekommen.

Um das Durchsuchen des Pipeline-Quellcodes zu erleichtern, erstelle einen symbolischen Link, der auf die ausgecheckte Kopie der Pipeline zeigt:

```bash
mkdir -p pipelines/nf-core
ln -s "$(echo $NXF_HOME/assets/.repos/nf-core/demo/clones/*/)" pipelines/nf-core/demo
```

Damit erstellst du eine Verknüpfung, mit der du den Code mit `tree -L 2 pipelines/nf-core/demo` erkunden oder Dateien direkt öffnen kannst.

#### 1.2.5. Übersicht der Code-Organisation

Du kannst entweder `tree` verwenden oder den Datei-Explorer nutzen, um das Verzeichnis `nf-core/demo` zu finden und zu öffnen.

```bash
tree -L 1 pipelines/nf-core/demo
```

??? abstract "Verzeichnisinhalt"

    ```console
    pipelines/nf-core/demo
    ├── assets
    ├── CHANGELOG.md
    ├── CITATIONS.md
    ├── CODE_OF_CONDUCT.md
    ├── conf
    ├── docs
    ├── LICENSE
    ├── main.nf
    ├── modules
    ├── modules.json
    ├── nextflow.config
    ├── nextflow_schema.json
    ├── nf-test.config
    ├── README.md
    ├── ro-crate-metadata.json
    ├── subworkflows
    ├── tests
    ├── tower.yml
    └── workflows

    7 directories, 12 files
    ```

Wie du siehst, ist dort einiges los – aber das meiste davon musst du dir nicht merken.

Kurz gesagt: Auf der obersten Ebene findest du eine README-Datei mit zusammenfassenden Informationen sowie Hilfsdateien, die Projektinformationen wie Lizenz, Beitragsrichtlinien, Zitierhinweise und Verhaltenskodex enthalten.
Die detaillierte Pipeline-Dokumentation befindet sich im Verzeichnis `docs`.
All diese Inhalte werden verwendet, um die Webseiten auf der nf-core-Website automatisch zu generieren – sie sind also immer auf dem neuesten Stand des Codes.

Für den Rest können wir drei funktionale Gruppen von Code-Dateien unterscheiden:

1. Pipeline-Code-Komponenten (`main.nf`, `workflows`, `subworkflows`, `modules`)
2. Pipeline-Konfiguration
3. Pipeline-Parameter / Eingaben und Validierung

Wir werden die Pipeline-Code-Komponenten in diesem Teil des Kurses nicht durchgehen, aber wir werden Elemente der Konfiguration und Validierung ansprechen, die für dich als Endnutzer\*in von nf-core-Pipelines relevant sein dürften.

!!! tip "Tipp"

    Du kannst den Quellcode jeder nf-core-Pipeline auch auf GitHub durchsuchen, z. B. [github.com/nf-core/demo](https://github.com/nf-core/demo).
    Jede nf-core-Pipeline folgt demselben Verzeichnislayout. Wenn du die Struktur einmal kennst, findest du Konfigurationsdateien, Module und Workflows für jede Pipeline auf die gleiche Weise.

Jetzt aber zur Ausführung der Pipeline!

### Fazit

Du weißt jetzt, wie du eine Pipeline über die nf-core-Website findest und eine lokale Kopie des Quellcodes herunterlädst.

### Wie geht es weiter?

Lerne, wie du eine nf-core-Pipeline mit minimalem Aufwand ausprobierst.

---

## 2. Die Pipeline mit ihrem Test-Profil ausprobieren

Praktischerweise wird jede nf-core-Pipeline mit einem Test-Profil geliefert.
Dies ist ein minimaler Satz von Konfigurationseinstellungen, mit denen die Pipeline einen kleinen Test-Datensatz aus dem Repository [nf-core/test-datasets](https://github.com/nf-core/test-datasets) verwenden kann.
Es ist eine großartige Möglichkeit, eine Pipeline schnell in kleinem Maßstab auszuprobieren.

!!! tip "Tipp"

    Das Konfigurationsprofil-System von Nextflow ermöglicht es dir, einfach zwischen verschiedenen Container-Engines oder Ausführungsumgebungen zu wechseln.
    Weitere Details findest du unter [Hello Nextflow Teil 6: Konfiguration](../hello_nextflow/06_hello_config.md).

### 2.1. Das Test-Profil untersuchen

Es ist gute Praxis, vor der Ausführung zu prüfen, was das Test-Profil einer Pipeline festlegt.
Das `test`-Profil für `nf-core/demo` befindet sich in der Konfigurationsdatei `conf/test.config`.
Du findest es lokal im Pipeline-Quellcode, den `nextflow pull` heruntergeladen hat, über den in Abschnitt 1.2.4 erstellten `pipelines`-Symlink:

```bash
code pipelines/nf-core/demo/conf/test.config
```

Hier ist der Inhalt dieser Datei:

```groovy title="conf/test.config" linenums="1" hl_lines="8 26"
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Nextflow config file for running minimal tests
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Defines input files and everything required to run a fast and simple pipeline test.

    Use as follows:
        nextflow run nf-core/demo -profile test,<docker/singularity> --outdir <OUTDIR>

----------------------------------------------------------------------------------------
*/

process {
    resourceLimits = [
        cpus: 2,
        memory: '4.GB',
        time: '1.h',
    ]
}

params {
    config_profile_name        = 'Test profile'
    config_profile_description = 'Minimal test dataset to check pipeline function'

    // Eingabedaten
    input                      = 'https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv'
}
```

Du wirst sofort bemerken, dass der Kommentarblock am Anfang ein Verwendungsbeispiel enthält, das zeigt, wie die Pipeline mit diesem Test-Profil ausgeführt wird.

```groovy title="conf/test.config" linenums="7"
    Use as follows:
        nextflow run nf-core/demo -profile test,<docker/singularity> --outdir <OUTDIR>
```

Wir müssen nur angeben, was zwischen den spitzen Klammern im Beispielbefehl steht: `<docker/singularity>` und `<OUTDIR>`.

Zur Erinnerung: `<docker/singularity>` bezieht sich auf die Wahl des Container-Systems. Alle nf-core-Pipelines sind so konzipiert, dass sie mit Containern (Docker, Singularity usw.) verwendet werden können, um Reproduzierbarkeit zu gewährleisten und Probleme bei der Software-Installation zu vermeiden.
Wir müssen also angeben, ob wir Docker oder Singularity zum Testen der Pipeline verwenden möchten.

Der Teil `--outdir <OUTDIR>` bezieht sich auf das Verzeichnis, in das Nextflow die Ausgaben der Pipeline schreibt.
Wir müssen einen Namen dafür angeben, den wir uns einfach ausdenken können.
Falls es noch nicht existiert, erstellt Nextflow es zur Laufzeit für uns.

Im Abschnitt nach dem Kommentarblock zeigt das Test-Profil, was für das Testen vorkonfiguriert wurde: Besonders wichtig ist, dass der Parameter `input` bereits auf einen Test-Datensatz zeigt, sodass wir keine eigenen Daten bereitstellen müssen.
Wenn du dem Link zur vorkonfigurierten Eingabe folgst, siehst du, dass es sich um eine CSV-Datei handelt, die Probenbezeichner und Dateipfade für mehrere experimentelle Proben enthält.

```csv title="samplesheet_test_illumina_amplicon.csv"
sample,fastq_1,fastq_2
SAMPLE1_PE,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample1_R1.fastq.gz,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample1_R2.fastq.gz
SAMPLE2_PE,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample2_R1.fastq.gz,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample2_R2.fastq.gz
SAMPLE3_SE,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample1_R1.fastq.gz,
SAMPLE3_SE,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample2_R1.fastq.gz,
```

Dies wird als Samplesheet bezeichnet und ist die häufigste Form der Eingabe für nf-core-Pipelines.
Mach dir keine Sorgen, wenn du mit den Datenformaten und -typen nicht vertraut bist – das ist für das Folgende nicht wichtig.

Wir haben jetzt alles, was wir brauchen, um die Pipeline auszuprobieren.

### 2.2. Die Pipeline ausführen

Wie oben erwähnt, können wir den Beispiel-Testbefehl fast unverändert verwenden. Wir müssen nur angeben, welches Software-Paketierungssystem wir verwenden möchten, und wie das Ausgabeverzeichnis heißen soll.
Hier verwenden wir Docker als Container-System und `demo-results` als Verzeichnisnamen.

Damit können wir den Testbefehl ausführen:

```bash
nextflow run nf-core/demo -profile test,docker --outdir demo-results
```

??? success "Befehlsausgabe"

    ```console
     N E X T F L O W   ~  version 26.04.4

    Downloading plugin nf-schema@2.7.2
    Launching `https://github.com/nf-core/demo` [cranky_curry] revision: 32893afef8 [master]


    ------------------------------------------------------
                                            ,--./,-.
            ___     __   __   __   ___     /,-._.--~'
      |\ | |__  __ /  ` /  \ |__) |__         }  {
      | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                            `._,._,'
      nf-core/demo 1.2.0
    ------------------------------------------------------

    Input/output options
      input                     : https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv
      outdir                    : demo-results

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      trace_report_suffix       : 2026-07-03_21-31-35

    Core Nextflow options
      revision                  : master
      runName                   : cranky_curry
      containerEngine           : docker
      launchDir                 : /workspaces/training/nfcore-use
      workDir                   : /workspaces/training/nfcore-use/work
      projectDir                : /workspaces/.nextflow/assets/.repos/nf-core/demo/clones/32893afef8076a03a2767a020b3f0cab2e0b40b2
      userName                  : root
      profile                   : test,docker
      configFiles               : /workspaces/.nextflow/assets/.repos/nf-core/demo/clones/32893afef8076a03a2767a020b3f0cab2e0b40b2/nextflow.config

    !! Only displaying parameters that differ from the pipeline defaults !!
    ------------------------------------------------------

    * The pipeline
        https://doi.org/10.5281/zenodo.12192442

    * The nf-core framework
        https://doi.org/10.1038/s41587-020-0439-x

    * Software dependencies
        https://github.com/nf-core/demo/blob/master/CITATIONS.md

    executor >  local (8)
    [ca/5b0f3e] NFCORE_DEMO:DEMO:FASTQC (SAMPLE3_SE)     | 3 of 3 ✔
    [b7/cb6812] NFCORE_DEMO:DEMO:SEQTK_TRIM (SAMPLE3_SE) | 3 of 3 ✔
    [ff/6ebd98] NFCORE_DEMO:DEMO:COWPY                   | 1 of 1 ✔
    [09/bbd1b4] NFCORE_DEMO:DEMO:MULTIQC (demo)          | 1 of 1 ✔
    -[nf-core/demo] Pipeline completed successfully-
    ```

Wenn deine Ausgabe damit übereinstimmt – herzlichen Glückwunsch! Du hast gerade deine erste nf-core-Pipeline ausgeführt.

Du wirst bemerken, dass die Konsolenausgabe viel umfangreicher ist als bei einer einfachen Nextflow-Pipeline.
Es gibt einen Header mit einer Zusammenfassung der Pipeline-Version, Eingaben und Ausgaben sowie einigen Konfigurationselementen.

!!! info "Info"

    Deine Ausgabe zeigt andere Zeitstempel, Ausführungsnamen und Dateipfade, aber die Gesamtstruktur und die Prozessausführung sollten ähnlich sein.

Beachte die Zeile nahe dem Anfang der Ausgabe:

```console
Launching `https://github.com/nf-core/demo` [cranky_curry] revision: 32893afef8 [master]
```

Diese zeigt dir, welche Revision der Pipeline verwendet wurde.
Da wir keine Version angegeben haben, verwendete Nextflow den neuesten Commit auf `master`.
Für reproduzierbare Ausführungen solltest du eine bestimmte Version mit dem Flag `-r` festlegen:

```bash
nextflow run nf-core/demo -r 1.2.0 -profile test,docker --outdir demo-results
```

So wird sichergestellt, dass jedes Mal derselbe Pipeline-Code verwendet wird, unabhängig von neuen Commits oder Releases.
In diesem Training lassen wir `-r` der Einfachheit halber weg, aber in der Produktion solltest du es immer angeben.

Weiter zur Ausführungsausgabe: Schauen wir uns die Zeilen an, die uns zeigen, welche Prozesse ausgeführt wurden:

```console
executor >  local (8)
[ca/5b0f3e] NFCORE_DEMO:DEMO:FASTQC (SAMPLE3_SE)     | 3 of 3 ✔
[b7/cb6812] NFCORE_DEMO:DEMO:SEQTK_TRIM (SAMPLE3_SE) | 3 of 3 ✔
[ff/6ebd98] NFCORE_DEMO:DEMO:COWPY                   | 1 of 1 ✔
[09/bbd1b4] NFCORE_DEMO:DEMO:MULTIQC (demo)          | 1 of 1 ✔
-[nf-core/demo] Pipeline completed successfully-
```

Das zeigt uns, dass vier Prozesse ausgeführt wurden, entsprechend den vier Tools, die auf der Pipeline-Dokumentationsseite der nf-core-Website aufgeführt sind: `FASTQC`, `SEQTK_TRIM`, `MULTIQC` und `COWPY`.

Die vollständigen Prozessnamen, wie sie hier angezeigt werden – z. B. `NFCORE_DEMO:DEMO:MULTIQC` – sind länger als das, was du möglicherweise im einführenden Hello Nextflow-Material gesehen hast.
Sie enthalten die Namen ihrer übergeordneten Workflows und spiegeln die Modularität des Pipeline-Codes wider.
Wenn du lernen möchtest, nf-core-artige Pipelines selbst zu entwickeln, schau dir den Kurs [Build with nf-core](../nfcore_build/index.md) an.

### 2.3. Die Ausgaben der Pipeline untersuchen

Schauen wir uns abschließend das Verzeichnis `demo-results` an, das von der Pipeline erzeugt wurde.

```bash
tree -L 2 demo-results
```

??? abstract "Verzeichnisinhalt"

    ```console
    demo-results
    ├── cowpy
    │   └── cowpy.txt
    ├── fastqc
    │   ├── SAMPLE1_PE
    │   ├── SAMPLE2_PE
    │   └── SAMPLE3_SE
    ├── fq
    │   ├── SAMPLE1_PE
    │   ├── SAMPLE2_PE
    │   └── SAMPLE3_SE
    ├── multiqc
    │   ├── multiqc_data
    │   └── multiqc_report.html
    └── pipeline_info
        ├── execution_report_2026-07-03_21-31-35.html
        ├── execution_timeline_2026-07-03_21-31-35.html
        ├── execution_trace_2026-07-03_21-31-35.txt
        ├── nf_core_demo_software_mqc_versions.yml
        ├── params_2026-07-03_21-31-43.json
        └── pipeline_dag_2026-07-03_21-31-35.html

    12 directories, 8 files
    ```

Das mag viel erscheinen.
Um mehr über die Ausgaben der `nf-core/demo`-Pipeline zu erfahren, schau dir die [Dokumentationsseite](https://nf-co.re/demo/1.2.0/docs/output/) an.

An dieser Stelle ist es wichtig zu beobachten, dass die Ergebnisse nach Modul organisiert sind, und es gibt zusätzlich ein Verzeichnis namens `pipeline_info`, das verschiedene zeitgestempelte Berichte über die Pipeline-Ausführung enthält.

Die Datei `execution_timeline_*` zeigt dir zum Beispiel, welche Prozesse in welcher Reihenfolge ausgeführt wurden und wie lange sie gedauert haben:

![Ausführungs-Timeline-Bericht](./img/execution_timeline.png)

!!! info "Info"

    Hier wurden die Aufgaben nicht parallel ausgeführt, da wir auf einem minimalistischen Rechner in Github Codespaces arbeiten.
    Um sie parallel laufen zu sehen, versuche die CPU-Zuweisung deines Codespace und die Ressourcenlimits in der Testkonfiguration zu erhöhen.

Diese Berichte werden automatisch für alle nf-core-Pipelines generiert.

### Fazit

Du weißt jetzt, wie du eine nf-core-Pipeline mit ihrem integrierten Test-Profil ausführst und wo du ihre Ausgaben findest.

### Wie geht es weiter?

Weiter zu [Teil 2](./02_configure_execution.md), wo du lernst, wie du die Pipeline-Ausführung konfigurierst.

---

## Zusammenfassung

In diesem Teil hast du gelernt:

- Eine nf-core-Pipeline zu finden und herunterzuladen sowie ihre Code-Struktur zu untersuchen
- Eine Pipeline mit ihrem integrierten Test-Profil auszuführen
