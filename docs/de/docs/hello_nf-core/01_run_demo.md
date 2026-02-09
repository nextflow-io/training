# Teil 1: Eine Demo-Pipeline ausführen

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Im ersten Teil des Hello nf-core Trainingskurses zeigen wir dir, wie du eine nf-core Pipeline findest und ausprobierst, verstehst, wie der Code organisiert ist, und erkennst, wie er sich von reinem Nextflow-Code unterscheidet, wie in [Hello Nextflow](../hello_nextflow/index.md) gezeigt.

Wir werden eine Pipeline namens nf-core/demo verwenden, die vom nf-core-Projekt als Teil seines Pipeline-Inventars zur Demonstration von Codestruktur und Tool-Operationen gepflegt wird.

Stelle sicher, dass dein Arbeitsverzeichnis auf `hello-nf-core/` gesetzt ist, wie auf der Seite [Erste Schritte](./00_orientation.md) beschrieben.

---

## 1. Die nf-core/demo Pipeline finden und abrufen

Beginnen wir damit, die nf-core/demo Pipeline auf der Projektwebsite unter [nf-co.re](https://nf-co.re) zu finden, die alle Informationen zentralisiert, wie: allgemeine Dokumentation und Hilfeartikel, Dokumentation für jede der Pipelines, Blogbeiträge, Veranstaltungsankündigungen und so weiter.

### 1.1. Die Pipeline auf der Website finden

Gehe in deinem Webbrowser zu [https://nf-co.re/pipelines/](https://nf-co.re/pipelines/) und tippe `demo` in die Suchleiste.

![Suchergebnisse](./img/search-results.png)

Klicke auf den Pipeline-Namen `demo`, um zur Pipeline-Dokumentationsseite zu gelangen.

Jede veröffentlichte Pipeline hat eine eigene Seite mit folgenden Dokumentationsabschnitten:

- **Introduction:** Eine Einführung und Übersicht der Pipeline
- **Usage:** Beschreibungen zur Ausführung der Pipeline
- **Parameters:** Gruppierte Pipeline-Parameter mit Beschreibungen
- **Output:** Beschreibungen und Beispiele der erwarteten Ausgabedateien
- **Results:** Beispiel-Ausgabedateien, die aus dem vollständigen Testdatensatz generiert wurden
- **Releases & Statistics:** Pipeline-Versionshistorie und Statistiken

Wann immer du erwägst, eine neue Pipeline zu übernehmen, solltest du zuerst die Pipeline-Dokumentation sorgfältig lesen, um zu verstehen, was sie tut und wie sie konfiguriert werden sollte, bevor du versuchst, sie auszuführen.

Schau dir das jetzt an und versuche herauszufinden:

- Welche Tools die Pipeline ausführen wird (Prüfe den Tab: `Introduction`)
- Welche Eingaben und Parameter die Pipeline akzeptiert oder benötigt (Prüfe den Tab: `Parameters`)
- Was sind die Ausgaben, die von der Pipeline produziert werden (Prüfe den Tab: `Output`)

#### 1.1.1. Pipeline-Übersicht

Der Tab `Introduction` bietet eine Übersicht der Pipeline, einschließlich einer visuellen Darstellung (genannt U-Bahn-Plan) und einer Liste von Tools, die als Teil der Pipeline ausgeführt werden.

![Pipeline U-Bahn-Plan](./img/nf-core-demo-subway-cropped.png)

1. Read QC (FASTQC)
2. Adapter and quality trimming (SEQTK_TRIM)
3. Present QC for raw reads (MULTIQC)

#### 1.1.2. Beispiel-Befehlszeile

Die Dokumentation bietet auch eine Beispiel-Eingabedatei (weiter unten ausführlicher besprochen) und eine Beispiel-Befehlszeile.

```bash
nextflow run nf-core/demo \
  -profile <docker/singularity/.../institute> \
  --input samplesheet.csv \
  --outdir <OUTDIR>
```

Du wirst bemerken, dass der Beispielbefehl KEINE Workflow-Datei angibt, sondern nur die Referenz zum Pipeline-Repository `nf-core/demo`.

Wenn auf diese Weise aufgerufen, nimmt Nextflow an, dass der Code auf eine bestimmte Weise organisiert ist.
Lass uns den Code abrufen, damit wir diese Struktur untersuchen können.

### 1.2. Den Pipeline-Code abrufen

Nachdem wir festgestellt haben, dass die Pipeline für unsere Zwecke geeignet zu sein scheint, lass sie uns ausprobieren.
Glücklicherweise macht es Nextflow einfach, Pipelines aus korrekt formatierten Repositories abzurufen, ohne etwas manuell herunterladen zu müssen.

Kehren wir zum Terminal zurück und führen Folgendes aus:

```bash
nextflow pull nf-core/demo
```

??? success "Befehlsausgabe"

    ```console
    Checking nf-core/demo ...
    downloaded from https://github.com/nf-core/demo.git - revision: 04060b4644 [master]
    ```

Nextflow führt einen `pull` des Pipeline-Codes durch, was bedeutet, dass es das vollständige Repository auf deine lokale Festplatte herunterlädt.

Um es klarzustellen: Du kannst dies mit jeder Nextflow-Pipeline tun, die entsprechend in GitHub eingerichtet ist, nicht nur mit nf-core Pipelines.
Allerdings ist nf-core die größte Open-Source-Sammlung von Nextflow-Pipelines.

Du kannst Nextflow eine Liste der Pipelines geben lassen, die du auf diese Weise abgerufen hast:

```bash
nextflow list
```

??? success "Befehlsausgabe"

    ```console
    nf-core/demo
    ```

Du wirst bemerken, dass sich die Dateien nicht in deinem aktuellen Arbeitsverzeichnis befinden.
Standardmäßig speichert Nextflow sie in `$NXF_HOME/assets`.

```bash
tree -L 2 $NXF_HOME/assets/
```

```console title="Verzeichnisinhalt"
/workspaces/.nextflow/assets/
└── nf-core
    └── demo

2 directories, 0 files
```

!!! note

    Der vollständige Pfad kann auf deinem System abweichen, wenn du nicht unsere Trainingsumgebung verwendest.

Nextflow hält den heruntergeladenen Quellcode absichtlich 'aus dem Weg' nach dem Prinzip, dass diese Pipelines eher wie Bibliotheken verwendet werden sollten als Code, mit dem du direkt interagieren würdest.

Für die Zwecke dieses Trainings möchten wir jedoch in der Lage sein, herumzustöbern und zu sehen, was darin ist.
Um das zu erleichtern, erstellen wir einen symbolischen Link zu diesem Ort von unserem aktuellen Arbeitsverzeichnis aus.

```bash
ln -s $NXF_HOME/assets pipelines
```

Dies erstellt eine Verknüpfung, die es einfacher macht, den gerade heruntergeladenen Code zu erkunden.

```bash
tree -L 2 pipelines
```

```console title="Verzeichnisinhalt"
pipelines
└── nf-core
    └── demo

2 directories, 0 files
```

Jetzt können wir bei Bedarf leichter in den Quellcode schauen.

Aber zuerst, lass uns unsere erste nf-core Pipeline ausführen!

### Fazit

Du weißt jetzt, wie du eine Pipeline über die nf-core Website findest und eine lokale Kopie des Quellcodes abrufst.

### Wie geht es weiter?

Lerne, wie du eine nf-core Pipeline mit minimalem Aufwand ausprobierst.

---

## 2. Die Pipeline mit ihrem Testprofil ausprobieren

Praktischerweise kommt jede nf-core Pipeline mit einem Testprofil.
Dies ist ein minimaler Satz von Konfigurationseinstellungen für die Pipeline, um mit einem kleinen Testdatensatz zu laufen, der im [nf-core/test-datasets](https://github.com/nf-core/test-datasets) Repository gehostet wird.
Es ist eine großartige Möglichkeit, eine Pipeline schnell im kleinen Maßstab auszuprobieren.

!!! note

    Das Konfigurationsprofil-System von Nextflow ermöglicht es dir, einfach zwischen verschiedenen Container-Engines oder Ausführungsumgebungen zu wechseln.
    Für weitere Details siehe [Hello Nextflow Teil 6: Konfiguration](../hello_nextflow/06_hello_config.md).

### 2.1. Das Testprofil untersuchen

Es ist gute Praxis zu prüfen, was das Testprofil einer Pipeline spezifiziert, bevor man sie ausführt.
Das `test`-Profil für `nf-core/demo` befindet sich in der Konfigurationsdatei `conf/test.config` und ist unten gezeigt.

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
        cpus: 4,
        memory: '4.GB',
        time: '1.h'
    ]
}

params {
    config_profile_name        = 'Test profile'
    config_profile_description = 'Minimal test dataset to check pipeline function'

    // Input data
    input  = 'https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv'

}
```

Du wirst sofort bemerken, dass der Kommentarblock oben ein Verwendungsbeispiel enthält, das zeigt, wie die Pipeline mit diesem Testprofil ausgeführt wird.

```groovy title="conf/test.config" linenums="7"
Use as follows:
        nextflow run nf-core/demo -profile test,<docker/singularity> --outdir <OUTDIR>
```

Die einzigen Dinge, die wir angeben müssen, sind das, was im Beispielbefehl zwischen spitzen Klammern gezeigt wird: `<docker/singularity>` und `<OUTDIR>`.

Zur Erinnerung: `<docker/singularity>` bezieht sich auf die Wahl des Container-Systems. Alle nf-core Pipelines sind so konzipiert, dass sie mit Containern (Docker, Singularity usw.) verwendbar sind, um Reproduzierbarkeit zu gewährleisten und Softwareinstallationsprobleme zu eliminieren.
Wir müssen also angeben, ob wir Docker oder Singularity verwenden möchten, um die Pipeline zu testen.

Der Teil `--outdir <OUTDIR>` bezieht sich auf das Verzeichnis, in das Nextflow die Ausgaben der Pipeline schreiben wird.
Wir müssen einen Namen dafür angeben, den wir einfach erfinden können.
Wenn es noch nicht existiert, wird Nextflow es zur Laufzeit für uns erstellen.

Weiter zum Abschnitt nach dem Kommentarblock: Das Testprofil zeigt uns, was für Tests vorkonfiguriert wurde: vor allem ist der Parameter `input` bereits so eingestellt, dass er auf einen Testdatensatz zeigt, sodass wir keine eigenen Daten bereitstellen müssen.
Wenn du dem Link zur vorkonfigurierten Eingabe folgst, wirst du sehen, dass es sich um eine CSV-Datei handelt, die Probenidentifikatoren und Dateipfade für mehrere experimentelle Proben enthält.

```csv title="samplesheet_test_illumina_amplicon.csv"
sample,fastq_1,fastq_2
SAMPLE1_PE,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample1_R1.fastq.gz,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample1_R2.fastq.gz
SAMPLE2_PE,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample2_R1.fastq.gz,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample2_R2.fastq.gz
SAMPLE3_SE,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample1_R1.fastq.gz,
SAMPLE3_SE,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample2_R1.fastq.gz,
```

Dies wird als Samplesheet bezeichnet und ist die häufigste Form der Eingabe für nf-core Pipelines.

!!! note

    Mach dir keine Sorgen, wenn du mit den Datenformaten und -typen nicht vertraut bist, das ist für das Folgende nicht wichtig.

Dies bestätigt also, dass wir alles haben, was wir brauchen, um die Pipeline auszuprobieren.

### 2.2. Die Pipeline ausführen

Entscheiden wir uns, Docker für das Container-System zu verwenden und `demo-results` als Ausgabeverzeichnis, und wir sind bereit, den Testbefehl auszuführen:

```bash
nextflow run nf-core/demo -profile docker,test --outdir demo-results
```

??? success "Befehlsausgabe"

    ```console
     N E X T F L O W   ~  version 25.04.3

    Launching `https://github.com/nf-core/demo` [magical_pauling] DSL2 - revision: db7f526ce1 [master]


    ------------------------------------------------------
                                            ,--./,-.
            ___     __   __   __   ___     /,-._.--~'
      |\ | |__  __ /  ` /  \ |__) |__         }  {
      | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                            `._,._,'
      nf-core/demo 1.0.2
    ------------------------------------------------------
    Input/output options
      input                     : https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv
      outdir                    : demo-results

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      trace_report_suffix       : 2025-11-21_04-57-41

    Core Nextflow options
      revision                  : master
      runName                   : magical_pauling
      containerEngine           : docker
      launchDir                 : /workspaces/training/hello-nf-core
      workDir                   : /workspaces/training/hello-nf-core/work
      projectDir                : /workspaces/.nextflow/assets/nf-core/demo
      userName                  : root
      profile                   : docker,test
      configFiles               : /workspaces/.nextflow/assets/nf-core/demo/nextflow.config

    !! Only displaying parameters that differ from the pipeline defaults !!
    ------------------------------------------------------
    * The pipeline
        https://doi.org/10.5281/zenodo.12192442

    * The nf-core framework
        https://doi.org/10.1038/s41587-020-0439-x

    * Software dependencies
        https://github.com/nf-core/demo/blob/master/CITATIONS.md


    executor >  local (7)
    [ff/a6976b] NFCORE_DEMO:DEMO:FASTQC (SAMPLE3_SE)     | 3 of 3 ✔
    [39/731ab7] NFCORE_DEMO:DEMO:SEQTK_TRIM (SAMPLE3_SE) | 3 of 3 ✔
    [7c/78d96e] NFCORE_DEMO:DEMO:MULTIQC                 | 1 of 1 ✔
    -[nf-core/demo] Pipeline completed successfully-
    ```

Wenn deine Ausgabe damit übereinstimmt, herzlichen Glückwunsch! Du hast gerade deine erste nf-core Pipeline ausgeführt.

Du wirst bemerken, dass es viel mehr Konsolenausgabe gibt als beim Ausführen einer einfachen Nextflow-Pipeline.
Es gibt einen Header, der eine Zusammenfassung der Pipeline-Version, Eingaben und Ausgaben sowie einige Konfigurationselemente enthält.

!!! note

    Deine Ausgabe wird unterschiedliche Zeitstempel, Ausführungsnamen und Dateipfade zeigen, aber die Gesamtstruktur und Prozessausführung sollte ähnlich sein.

Weiter zur Ausführungsausgabe: Schauen wir uns die Zeilen an, die uns sagen, welche Prozesse ausgeführt wurden:

```console
    [ff/a6976b] NFCORE_DEMO:DEMO:FASTQC (SAMPLE3_SE)     | 3 of 3 ✔
    [39/731ab7] NFCORE_DEMO:DEMO:SEQTK_TRIM (SAMPLE3_SE) | 3 of 3 ✔
    [7c/78d96e] NFCORE_DEMO:DEMO:MULTIQC                 | 1 of 1 ✔
```

Dies sagt uns, dass drei Prozesse ausgeführt wurden, entsprechend den drei Tools, die auf der Pipeline-Dokumentationsseite der nf-core Website gezeigt werden: FASTQC, SEQTK_TRIM und MULTIQC.

Die vollständigen Prozessnamen, wie sie hier gezeigt werden, wie `NFCORE_DEMO:DEMO:MULTIQC`, sind länger als das, was du möglicherweise im einführenden Hello Nextflow Material gesehen hast.
Diese enthalten die Namen ihrer übergeordneten Workflows und spiegeln die Modularität des Pipeline-Codes wider.
Wir werden gleich etwas ausführlicher darauf eingehen.

### 2.3. Die Ausgaben der Pipeline untersuchen

Schauen wir uns schließlich das Verzeichnis `demo-results` an, das von der Pipeline produziert wurde.

```bash
tree -L 2 demo-results
```

??? abstract "Verzeichnisinhalt"

    ```console
    demo-results
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
    │   ├── multiqc_plots
    │   └── multiqc_report.html
    └── pipeline_info
        ├── execution_report_2025-11-21_04-57-41.html
        ├── execution_timeline_2025-11-21_04-57-41.html
        ├── execution_trace_2025-11-21_04-57-41.txt
        ├── nf_core_demo_software_mqc_versions.yml
        ├── params_2025-11-21_04-57-46.json
        └── pipeline_dag_2025-11-21_04-57-41.html
    ```

Das mag nach viel aussehen.
Um mehr über die Ausgaben der `nf-core/demo` Pipeline zu erfahren, schau dir ihre [Dokumentationsseite](https://nf-co.re/demo/1.0.2/docs/output/) an.

In diesem Stadium ist es wichtig zu beobachten, dass die Ergebnisse nach Modul organisiert sind, und es gibt zusätzlich ein Verzeichnis namens `pipeline_info`, das verschiedene mit Zeitstempel versehene Berichte über die Pipeline-Ausführung enthält.

Zum Beispiel zeigt dir die Datei `execution_timeline_*`, welche Prozesse ausgeführt wurden, in welcher Reihenfolge und wie lange sie zur Ausführung benötigt haben:

![Ausführungs-Timeline-Bericht](./img/execution_timeline.png)

!!! note

    Hier wurden die Aufgaben nicht parallel ausgeführt, weil wir auf einer minimalistischen Maschine in Github Codespaces laufen.
    Um diese parallel laufen zu sehen, versuche die CPU-Zuweisung deines Codespace und die Ressourcenlimits in der Testkonfiguration zu erhöhen.

Diese Berichte werden automatisch für alle nf-core Pipelines generiert.

### Fazit

Du weißt, wie du eine nf-core Pipeline mit ihrem eingebauten Testprofil ausführst und wo du ihre Ausgaben findest.

### Wie geht es weiter?

Lerne, wie der Pipeline-Code organisiert ist.

---

## 3. Die Pipeline-Codestruktur untersuchen

Nachdem wir die Pipeline erfolgreich als Benutzer\*innen ausgeführt haben, wechseln wir nun die Perspektive und schauen uns an, wie nf-core Pipelines intern strukturiert sind.

Das nf-core-Projekt setzt strenge Richtlinien durch, wie Pipelines strukturiert sind und wie der Code organisiert, konfiguriert und dokumentiert wird.
Zu verstehen, wie das alles organisiert ist, ist der erste Schritt zur Entwicklung deiner eigenen nf-core-kompatiblen Pipelines, was wir in Teil 2 dieses Kurses angehen werden.

Schauen wir uns an, wie der Pipeline-Code im `nf-core/demo` Repository organisiert ist, unter Verwendung des Symlinks `pipelines`, den wir zuvor erstellt haben.

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
    ```

Da ist eine Menge los, also gehen wir das Schritt für Schritt an.

Zunächst bemerken wir, dass du auf der obersten Ebene eine README-Datei mit zusammenfassenden Informationen sowie Zusatzdateien finden kannst, die Projektinformationen wie Lizenzierung, Beitragsrichtlinien, Zitation und Verhaltenskodex zusammenfassen.
Detaillierte Pipeline-Dokumentation befindet sich im Verzeichnis `docs`.
All dieser Inhalt wird verwendet, um die Webseiten auf der nf-core Website programmatisch zu generieren, sodass sie immer mit dem Code auf dem neuesten Stand sind.

Nun, für den Rest werden wir unsere Erkundung in drei Phasen unterteilen:

1. Pipeline-Code-Komponenten (`main.nf`, `workflows`, `subworkflows`, `modules`)
2. Pipeline-Konfiguration
3. Eingaben und Validierung

Beginnen wir mit den Pipeline-Code-Komponenten.
Wir werden uns auf die Dateihierarchie und strukturelle Organisation konzentrieren, anstatt in den Code innerhalb einzelner Dateien einzutauchen.

### 3.1. Pipeline-Code-Komponenten

Die Standard-nf-core Pipeline-Code-Organisation folgt einer modularen Struktur, die darauf ausgelegt ist, die Wiederverwendung von Code zu maximieren, wie in [Hello Modules](../hello_nextflow/04_hello_modules.md), Teil 4 des [Hello Nextflow](../hello_nextflow/index.md) Kurses eingeführt, obwohl dies in echter nf-core-Manier mit etwas zusätzlicher Komplexität implementiert ist.
Insbesondere machen nf-core Pipelines reichlich Gebrauch von Subworkflows, d.h. Workflow-Skripten, die von einem übergeordneten Workflow importiert werden.

Das mag etwas abstrakt klingen, also schauen wir uns an, wie dies in der Praxis in der `nf-core/demo` Pipeline verwendet wird.

!!! note

    Wir werden nicht den tatsächlichen Code durchgehen, _wie_ diese modularen Komponenten verbunden sind, weil es einige zusätzliche Komplexität im Zusammenhang mit der Verwendung von Subworkflows gibt, die verwirrend sein kann, und das Verständnis davon ist in diesem Stadium des Trainings nicht notwendig.
    Vorerst konzentrieren wir uns auf die Gesamtorganisation und Logik.

#### 3.1.1. Allgemeine Übersicht

So sehen die Beziehungen zwischen den relevanten Code-Komponenten für die `nf-core/demo` Pipeline aus:

<figure class="excalidraw">
    --8<-- "docs/en/docs/hello_nf-core/img/nf-core_demo_code_organization.svg"
</figure>

Es gibt ein sogenanntes _Entrypoint_-Skript namens `main.nf`, das als Wrapper für zwei Arten von verschachtelten Workflows fungiert: den Workflow, der die eigentliche Analyselogik enthält, der sich unter `workflows/` befindet und `demo.nf` heißt, und eine Reihe von Housekeeping-Workflows, die sich unter `subworkflows/` befinden.
Der `demo.nf` Workflow ruft **Module** auf, die sich unter `modules/` befinden; diese enthalten die **Prozesse**, die die eigentlichen Analyseschritte durchführen werden.

!!! note

    Subworkflows sind nicht auf Housekeeping-Funktionen beschränkt, und sie können Prozessmodule verwenden.

    Die hier gezeigte `nf-core/demo` Pipeline ist zufällig auf der einfacheren Seite des Spektrums, aber andere nf-core Pipelines (wie `nf-core/rnaseq`) nutzen Subworkflows, die an der eigentlichen Analyse beteiligt sind.

Schauen wir uns nun diese Komponenten der Reihe nach an.

#### 3.1.2. Das Entrypoint-Skript: `main.nf`

Das `main.nf` Skript ist der Entrypoint, von dem Nextflow startet, wenn wir `nextflow run nf-core/demo` ausführen.
Das bedeutet, wenn du `nextflow run nf-core/demo` ausführst, um die Pipeline zu starten, findet und führt Nextflow automatisch das `main.nf` Skript aus.
Dies funktioniert für jede Nextflow-Pipeline, die dieser konventionellen Benennung und Struktur folgt, nicht nur für nf-core Pipelines.

Die Verwendung eines Entrypoint-Skripts macht es einfach, standardisierte 'Housekeeping'-Subworkflows vor und nach dem eigentlichen Analyseskript auszuführen.
Wir werden diese durchgehen, nachdem wir den eigentlichen Analyse-Workflow und seine Module überprüft haben.

#### 3.1.3. Das Analyseskript: `workflows/demo.nf`

Der `workflows/demo.nf` Workflow ist der Ort, an dem die zentrale Logik der Pipeline gespeichert ist.
Er ist ähnlich strukturiert wie ein normaler Nextflow-Workflow, außer dass er so konzipiert ist, dass er von einem übergeordneten Workflow aufgerufen wird, was einige zusätzliche Funktionen erfordert.
Wir werden die relevanten Unterschiede im nächsten Teil dieses Kurses behandeln, wenn wir die Konvertierung der einfachen Hello-Pipeline aus Hello Nextflow in eine nf-core-kompatible Form angehen.

Der `demo.nf` Workflow ruft **Module** auf, die sich unter `modules/` befinden, die wir als Nächstes überprüfen werden.

!!! note

    Einige nf-core Analyse-Workflows zeigen zusätzliche Verschachtelungsebenen, indem sie Subworkflows auf niedrigerer Ebene aufrufen.
    Dies wird hauptsächlich verwendet, um zwei oder mehr Module, die häufig zusammen verwendet werden, in leicht wiederverwendbare Pipeline-Segmente zu verpacken.
    Du kannst einige Beispiele sehen, indem du verfügbare [nf-core Subworkflows](https://nf-co.re/subworkflows/) auf der nf-core Website durchsuchst.

    Wenn das Analyseskript Subworkflows verwendet, werden diese unter dem Verzeichnis `subworkflows/` gespeichert.

#### 3.1.4. Die Module

Die Module sind der Ort, an dem der Prozesscode lebt, wie in [Teil 4 des Hello Nextflow Trainingskurses](../hello_nextflow/04_hello_modules.md) beschrieben.

Im nf-core-Projekt werden Module unter Verwendung einer mehrstufigen verschachtelten Struktur organisiert, die sowohl ihre Herkunft als auch ihren Inhalt widerspiegelt.
Auf der obersten Ebene werden Module entweder als `nf-core` oder `local` (nicht Teil des nf-core-Projekts) unterschieden und dann weiter in ein Verzeichnis platziert, das nach dem/den Tool(s) benannt ist, das/die sie umschließen.
Wenn das Tool zu einem Toolkit gehört (d.h. einem Paket, das mehrere Tools enthält), gibt es eine Zwischenverzeichnisebene, die nach dem Toolkit benannt ist.

Du kannst dies in der Praxis auf die `nf-core/demo` Pipeline-Module angewendet sehen:

```bash
tree -L 3 pipelines/nf-core/demo/modules
```

??? abstract "Verzeichnisinhalt"

    ```console
    pipelines/nf-core/demo/modules
    └── nf-core
        ├── fastqc
        │   ├── environment.yml
        │   ├── main.nf
        │   ├── meta.yml
        │   └── tests
        ├── multiqc
        │   ├── environment.yml
        │   ├── main.nf
        │   ├── meta.yml
        │   └── tests
        └── seqtk
            └── trim

    7 directories, 6 files
    ```

Hier siehst du, dass die Module `fastqc` und `multiqc` auf der obersten Ebene innerhalb der `nf-core` Module sitzen, während das Modul `trim` unter dem Toolkit sitzt, zu dem es gehört, `seqtk`.
In diesem Fall gibt es keine `local` Module.

Die Modul-Codedatei, die den Prozess beschreibt, heißt immer `main.nf` und wird von Tests und `.yml`-Dateien begleitet, die wir vorerst ignorieren werden.

Zusammengenommen sind der Entrypoint-Workflow, der Analyse-Workflow und die Module ausreichend, um die 'interessanten' Teile der Pipeline auszuführen.
Wir wissen jedoch, dass es dort auch Housekeeping-Subworkflows gibt, also schauen wir uns diese jetzt an.

#### 3.1.5. Die Housekeeping-Subworkflows

Wie Module werden Subworkflows in `local` und `nf-core` Verzeichnisse unterschieden, und jeder Subworkflow hat seine eigene verschachtelte Verzeichnisstruktur mit seinem eigenen `main.nf` Skript, Tests und `.yml`-Datei.

```bash
tree -L 3 pipelines/nf-core/demo/subworkflows
```

??? abstract "Verzeichnisinhalt"

    ```console
    pipelines/nf-core/demo/subworkflows
    ├── local
    │   └── utils_nfcore_demo_pipeline
    │       └── main.nf
    └── nf-core
        ├── utils_nextflow_pipeline
        │   ├── main.nf
        │   ├── meta.yml
        │   └── tests
        ├── utils_nfcore_pipeline
        │   ├── main.nf
        │   ├── meta.yml
        │   └── tests
        └── utils_nfschema_plugin
            ├── main.nf
            ├── meta.yml
            └── tests

    9 directories, 7 files
    ```

Wie oben erwähnt, enthält die `nf-core/demo` Pipeline keine analysespezifischen Subworkflows, sodass alle Subworkflows, die wir hier sehen, sogenannte 'Housekeeping'- oder 'Utility'-Workflows sind, wie durch das Präfix `utils_` in ihren Namen gekennzeichnet.
Diese Subworkflows sind es, die unter anderem den schicken nf-core Header in der Konsolenausgabe produzieren.

!!! tip

    Abgesehen von ihrem Benennungsmuster ist ein weiterer Hinweis darauf, dass diese Subworkflows keine wirklich analysebezogene Funktion ausführen, dass sie überhaupt keine Prozesse aufrufen.

Dies vervollständigt die Zusammenfassung der Kern-Code-Komponenten, die die `nf-core/demo` Pipeline ausmachen.
Schauen wir uns nun die verbleibenden Elemente an, über die du ein wenig wissen solltest, bevor du in die Entwicklung eintauchst: Pipeline-Konfiguration und Eingabevalidierung.

### 3.2. Pipeline-Konfiguration

Du hast zuvor gelernt, dass Nextflow viele Optionen zur Konfiguration der Pipeline-Ausführung bietet, sei es in Bezug auf Eingaben und Parameter, Rechenressourcen und andere Aspekte der Orchestrierung.
Das nf-core-Projekt wendet hochstandardisierte Richtlinien für die Pipeline-Konfiguration an, die darauf abzielen, auf Nextflows flexiblen Anpassungsoptionen auf eine Weise aufzubauen, die größere Konsistenz und Wartbarkeit über Pipelines hinweg bietet.

Die zentrale Konfigurationsdatei `nextflow.config` wird verwendet, um Standardwerte für Parameter und andere Konfigurationsoptionen festzulegen.
Die Mehrheit dieser Konfigurationsoptionen wird standardmäßig angewendet, während andere (z.B. Software-Abhängigkeitsprofile) als optionale Profile enthalten sind.

Es gibt mehrere zusätzliche Konfigurationsdateien, die im Ordner `conf` gespeichert sind und die standardmäßig oder optional als Profile zur Konfiguration hinzugefügt werden können:

- `base.config`: Eine 'leere Tafel'-Konfigurationsdatei, geeignet für die allgemeine Verwendung in den meisten Hochleistungsrechenumgebungen. Dies definiert breite Bins der Ressourcennutzung, die beispielsweise bequem auf Module angewendet werden können.
- `modules.config`: Zusätzliche Modul-Direktiven und Argumente.
- `test.config`: Ein Profil zum Ausführen der Pipeline mit minimalen Testdaten, das wir verwendet haben, als wir die Demo-Pipeline ausgeführt haben.
- `test_full.config`: Ein Profil zum Ausführen der Pipeline mit einem vollständigen Testdatensatz.

Wir werden einige dieser Dateien später im Kurs berühren.

### 3.3. Eingaben und Validierung

Wie wir zuvor festgestellt haben, als wir das Testprofil der `nf-core/demo` Pipeline untersucht haben, ist sie so konzipiert, dass sie als Eingabe ein Samplesheet mit Dateipfaden und Probenidentifikatoren nimmt.
Die Dateipfade verweisen auf echte Daten, die sich im `nf-core/test-datasets` Repository befinden.

Ein Beispiel-Samplesheet wird auch unter dem Verzeichnis `assets` bereitgestellt, obwohl die Pfade in diesem nicht echt sind.

```csv title="assets/samplesheet.csv" linenums="1"
sample,fastq_1,fastq_2
SAMPLE_PAIRED_END,/path/to/fastq/files/AEG588A1_S1_L002_R1_001.fastq.gz,/path/to/fastq/files/AEG588A1_S1_L002_R2_001.fastq.gz
SAMPLE_SINGLE_END,/path/to/fastq/files/AEG588A4_S4_L003_R1_001.fastq.gz,

```

Dieses spezielle Samplesheet ist ziemlich einfach, aber einige Pipelines laufen auf Samplesheets, die komplexer sind, mit viel mehr Metadaten, die mit den primären Eingaben verbunden sind.

Leider, weil diese Dateien schwer mit dem Auge zu überprüfen sein können, ist eine unsachgemäße Formatierung von Eingabedaten eine sehr häufige Quelle für Pipeline-Fehler.
Ein verwandtes Problem ist, wenn Parameter falsch bereitgestellt werden.

Die Lösung für diese Probleme besteht darin, automatisierte Validierungsprüfungen für alle Eingabedateien durchzuführen, um sicherzustellen, dass sie die erwarteten Arten von Informationen enthalten, korrekt formatiert sind, und für Parameter, um sicherzustellen, dass sie vom erwarteten Typ sind.
Dies wird Eingabevalidierung genannt und sollte idealerweise durchgeführt werden, _bevor_ versucht wird, eine Pipeline auszuführen, anstatt darauf zu warten, dass die Pipeline fehlschlägt, um herauszufinden, dass es ein Problem mit den Eingaben gab.

Genau wie bei der Konfiguration ist das nf-core-Projekt sehr meinungsstark in Bezug auf Eingabevalidierung und empfiehlt die Verwendung des [nf-schema Plugins](https://nextflow-io.github.io/nf-schema/latest/), eines Nextflow-Plugins, das umfassende Validierungsfähigkeiten für Nextflow-Pipelines bietet.

Wir werden dieses Thema in Teil 5 dieses Kurses ausführlicher behandeln.
Sei dir vorerst bewusst, dass es zwei JSON-Dateien für diesen Zweck gibt: `nextflow_schema.json` und `assets/schema_input.json`.

Die `nextflow_schema.json` ist eine Datei, die verwendet wird, um Informationen über die Pipeline-Parameter einschließlich Typ, Beschreibung und Hilfetext in einem maschinenlesbaren Format zu speichern.
Dies wird für verschiedene Zwecke verwendet, einschließlich automatisierter Parametervalidierung, Hilfetextgenerierung und interaktiver Parameterformular-Darstellung in UI-Schnittstellen.

Die `schema_input.json` ist eine Datei, die verwendet wird, um die Struktur des Eingabe-Samplesheets zu definieren.
Jede Spalte kann einen Typ, ein Muster, eine Beschreibung und einen Hilfetext in einem maschinenlesbaren Format haben.
Das Schema wird für verschiedene Zwecke verwendet, einschließlich automatisierter Validierung und Bereitstellung hilfreicher Fehlermeldungen.

### Fazit

Du weißt, was die Hauptkomponenten einer nf-core Pipeline sind und wie der Code organisiert ist; wo sich die Hauptelemente der Konfiguration befinden; und du bist dir bewusst, wofür Eingabevalidierung da ist.

### Wie geht es weiter?

Mach eine Pause! Das war viel. Wenn du dich erfrischt und bereit fühlst, gehe zum nächsten Abschnitt über, um das Gelernte anzuwenden und eine nf-core-kompatible Pipeline zu schreiben.

!!! tip

    Wenn du lernen möchtest, wie man Workflows mit Subworkflows komponiert, bevor du zum nächsten Teil übergehst, schau dir die [Workflows of Workflows](../side_quests/workflows_of_workflows.md) Side Quest an.
