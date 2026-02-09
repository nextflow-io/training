# Teil 2: nf-core/molkart ausführen

In Teil 1 haben wir einen einfachen Hello-World-Workflow ausgeführt, um die Grundlagen der Nextflow-Ausführung zu verstehen.
Jetzt führen wir eine echte Bioimaging-Pipeline aus: **nf-core/molkart**.

Diese Pipeline verarbeitet Molecular-Cartography-Daten zur räumlichen Transkriptomik von Resolve Bioscience.
Die Nextflow-Muster, die du hier lernst, gelten jedoch für jede nf-core-Pipeline oder Produktions-Workflow.

## 1. nf-core-Pipelines verstehen

Bevor wir die Pipeline ausführen, lass uns verstehen, was nf-core ist und warum es für die Ausführung von Workflows wichtig ist.

### 1.1. Was ist nf-core?

[nf-core](https://nf-co.re/) ist eine Community-gesteuerte Sammlung hochwertiger Nextflow-Pipelines.
Alle nf-core-Pipelines folgen derselben Struktur und denselben Konventionen. Das bedeutet: Wenn du eine ausführen kannst, kannst du alle ausführen.

Hauptmerkmale von nf-core-Pipelines:

- **Standardisierte Struktur**: Alle Pipelines haben konsistente Parameternamen und Verwendungsmuster
- **Integrierte Testdaten**: Jede Pipeline enthält Testprofile zur schnellen Validierung
- **Umfassende Dokumentation**: Detaillierte Nutzungsanweisungen und Parameterbeschreibungen
- **Qualitätskontrolle**: Automatisierte QC-Berichte mit MultiQC
- **Container-Unterstützung**: Vorgefertigte Container für Reproduzierbarkeit

!!! tip "Möchtest du mehr über nf-core erfahren?"

    Für eine ausführliche Einführung in die nf-core-Pipeline-Entwicklung schau dir den Kurs [Hello nf-core](../../hello_nf-core/index.md) an.
    Er behandelt, wie man nf-core-Pipelines von Grund auf erstellt und anpasst.

### 1.2. Die molkart-Pipeline

![nf-core/molkart-Pipeline](img/molkart.png)

Die [nf-core/molkart](https://nf-co.re/molkart)-Pipeline verarbeitet Bilddaten zur räumlichen Transkriptomik in mehreren Schritten:

1. **Bildvorverarbeitung**: Auffüllen von Gittermustern und optional Kontrastverbesserung
2. **Zellsegmentierung**: Mehrere Algorithmusoptionen (Cellpose, Mesmer, ilastik, Stardist)
3. **Spot-Zuordnung**: Transkript-Spots segmentierten Zellen zuordnen
4. **Qualitätskontrolle**: Umfassende QC-Berichte generieren

Die wichtigsten Ausgaben sind:

- Zell-nach-Transkript-Zähltabellen
- Segmentierungsmasken
- MultiQC-Qualitätskontrollbericht

---

## 2. molkart mit Testdaten ausführen

Bevor wir beginnen, klonen wir das molkart-Repository lokal, damit wir den Code inspizieren können:

```bash
cd /workspaces/training/nf4-science/imaging
git clone --branch 1.2.0 --depth 1 https://github.com/nf-core/molkart
```

Dies erstellt ein `molkart/`-Verzeichnis mit dem vollständigen Pipeline-Quellcode.

!!! note "Warum klonen wir lokal?"

    Normalerweise würdest du nf-core-Pipelines direkt von GitHub mit `nextflow run nf-core/molkart -r 1.2.0` ausführen.
    Nextflow lädt die angeforderte Pipeline-Version automatisch für dich nach `$HOME/.nextflow/assets/nf-core/molkart` herunter und führt sie von dort aus.
    Für dieses Training klonen wir die Pipeline jedoch in ein anderes lokales Verzeichnis, damit wir den Code leichter inspizieren können.

### 2.1. Container-Anforderungen verstehen

Bevor wir die vollständige Pipeline ausführen, lass uns lernen, warum Container für nf-core-Pipelines unverzichtbar sind.

Versuchen wir, die Pipeline mit dem Testdatensatz und den Parametern aus der molkart-Testkonfiguration auszuführen:

```bash
nextflow run ./molkart \
  --input 'data/samplesheet.csv' \
  --mindagap_tilesize 90 \
  --mindagap_boxsize 7 \
  --mindagap_loopnum 100 \
  --clahe_pyramid_tile 368 \
  --segmentation_method "mesmer,cellpose,stardist" \
  --outdir results
```

Lass uns diese Parameter aufschlüsseln:

- `--input`: Pfad zum Samplesheet mit Sample-Metadaten
- `--mindagap_tilesize`, `--mindagap_boxsize`, `--mindagap_loopnum`: Parameter zum Auffüllen von Gittermustern
- `--clahe_pyramid_tile`: Kernel-Größe für Kontrastverbesserung
- `--segmentation_method`: Welche(r) Algorithmus/Algorithmen für die Zellsegmentierung verwendet werden soll(en)
- `--outdir`: Wo die Ergebnisse gespeichert werden sollen

!!! Warning "Dieser Befehl wird fehlschlagen - das ist beabsichtigt!"

    Wir führen dies absichtlich ohne Container aus, um zu zeigen, warum sie benötigt werden.

Nach einigen Momenten siehst du einen Fehler wie diesen:

??? failure "Befehlsausgabe"

    ```console
    ERROR ~ Error executing process > 'NFCORE_MOLKART:MOLKART:MINDAGAP_DUPLICATEFINDER (mem_only)'

    Caused by:
      Process `NFCORE_MOLKART:MOLKART:MINDAGAP_DUPLICATEFINDER (mem_only)` terminated with an error exit status (127)

    Command executed:

      duplicate_finder.py \
          spots.txt \
          90

    Command exit status:
      127

    Command error:
      .command.sh: line 3: duplicate_finder.py: command not found
    ```

**Was passiert hier?**

Der Fehler `command not found` (Exit-Status 127) bedeutet, dass Nextflow versucht hat, `duplicate_finder.py` auszuführen, es aber nicht auf deinem System finden konnte.
Das liegt daran, dass:

1. Die Pipeline erwartet, dass spezialisierte Bioinformatik-Software installiert ist
2. Diese Tools (wie `duplicate_finder.py`, `apply_clahe.dask.py` usw.) nicht Teil von Standard-Linux-Distributionen sind
3. Ohne Container versucht Nextflow, Befehle direkt auf deinem lokalen Rechner auszuführen

**Woher sollen diese Tools kommen?**

Lass uns eines der Prozessmodule inspizieren, um zu sehen, wie es seine Software-Anforderungen deklariert.

Öffne das CLAHE-Vorverarbeitungsmodul:

```bash
code molkart/modules/local/clahe/main.nf
```

Schau dir Zeile 5 an - du siehst:

```groovy
container 'ghcr.io/schapirolabor/molkart-local:v0.0.4'
```

Diese Zeile sagt Nextflow: "Um diesen Prozess auszuführen, verwende das Docker-Image `ghcr.io/schapirolabor/molkart-local:v0.0.4`, das alle erforderlichen Tools enthält."

Jeder Prozess deklariert, welches Container-Image seine benötigten Tools bereitstellt.
Nextflow verwendet diese Container jedoch nur, wenn du es ihm sagst!

**Die Lösung: Docker in der Konfiguration aktivieren**

### 2.2. Docker konfigurieren und die Pipeline starten

Um Docker zu aktivieren, müssen wir `docker.enabled` von `false` auf `true` in der Datei `nextflow.config` ändern.

Öffne die Konfigurationsdatei:

```bash
code nextflow.config
```

Ändere `docker.enabled = false` zu `docker.enabled = true`:

```groovy
docker.enabled = true
process {
    resourceLimits = [
        cpus: 2,
        memory: '7.GB',
    ]
}
```

Führe die Pipeline jetzt erneut mit demselben Befehl aus:

```bash
nextflow run ./molkart \
  --input 'data/samplesheet.csv' \
  --mindagap_tilesize 90 \
  --mindagap_boxsize 7 \
  --mindagap_loopnum 100 \
  --clahe_pyramid_tile 368 \
  --segmentation_method "cellpose,mesmer,stardist" \
  --outdir results
```

Dieses Mal wird Nextflow:

1. Die Einstellung `docker.enabled = true` aus der Konfiguration lesen
2. Die erforderlichen Docker-Images herunterladen (nur beim ersten Mal)
3. Jeden Prozess innerhalb seines angegebenen Containers ausführen
4. Erfolgreich ausgeführt werden, weil alle Tools innerhalb der Container verfügbar sind

!!! Tip "Warum Container wichtig sind"

    Die meisten nf-core-Pipelines **erfordern** Containerisierung (Docker, Singularity, Podman usw.), weil:

    - Sie spezialisierte Bioinformatik-Software verwenden, die in Standardumgebungen nicht verfügbar ist
    - Container Reproduzierbarkeit gewährleisten - exakt dieselben Software-Versionen laufen überall
    - Du nicht manuell Dutzende von Tools und ihre Abhängigkeiten installieren musst

    Weitere Details zu Containern in Nextflow findest du unter [Hello Containers](../../hello_nextflow/05_hello_containers.md) aus dem Hello-Nextflow-Training.

### 2.3. Ausführung überwachen

Während die Pipeline läuft, siehst du eine Ausgabe ähnlich dieser:

??? success "Befehlsausgabe"

    ```console
    Nextflow 25.04.8 is available - Please consider updating your version to it

    N E X T F L O W   ~  version 25.04.3

    Launching `https://github.com/nf-core/molkart` [soggy_kalam] DSL2 - revision: 5e54b29cb3 [dev]


    ------------------------------------------------------
                                            ,--./,-.
            ___     __   __   __   ___     /,-._.--~'
      |\ | |__  __ /  ` /  \ |__) |__         }  {
      | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                            `._,._,'
      nf-core/molkart 1.2.0dev
    ------------------------------------------------------
    Segmentation methods and options
      segmentation_method       : mesmer,cellpose,stardist

    Image preprocessing
      mindagap_boxsize          : 7
      mindagap_loopnum          : 100
      clahe_kernel              : 25
      mindagap_tilesize         : 90
      clahe_pyramid_tile        : 368

    Input/output options
      input                     : https://raw.githubusercontent.com/nf-core/test-datasets/molkart/test_data/samplesheets/samplesheet_membrane.csv
      outdir                    : results

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      trace_report_suffix       : 2025-10-18_22-22-21

    Core Nextflow options
      revision                  : dev
      runName                   : soggy_kalam
      containerEngine           : docker
      launchDir                 : /workspaces/training/nf4-science/imaging
      workDir                   : /workspaces/training/nf4-science/imaging/work
      projectDir                : /workspaces/.nextflow/assets/nf-core/molkart
      userName                  : root
      profile                   : docker,test
      configFiles               :

    !! Only displaying parameters that differ from the pipeline defaults !!
    ------------------------------------------------------
    * The pipeline
        https://doi.org/10.5281/zenodo.10650748

    * The nf-core framework
        https://doi.org/10.1038/s41587-020-0439-x

    * Software dependencies
        https://github.com/nf-core/molkart/blob/master/CITATIONS.md

    executor >  local (22)
    [c1/da5009] NFCORE_MOLKART:MOLKART:MINDAGAP_MINDAGAP (mem_only)        [100%] 2 of 2 ✔
    [73/8f5e8a] NFCORE_MOLKART:MOLKART:CLAHE (mem_only)                    [100%] 2 of 2 ✔
    [ec/8f84d5] NFCORE_MOLKART:MOLKART:CREATE_STACK (mem_only)             [100%] 1 of 1 ✔
    [a2/99349b] NFCORE_MOLKART:MOLKART:MINDAGAP_DUPLICATEFINDER (mem_only) [100%] 1 of 1 ✔
    [95/c9b4b1] NFCORE_MOLKART:MOLKART:DEEPCELL_MESMER (mem_only)          [100%] 1 of 1 ✔
    [d4/1ebd1e] NFCORE_MOLKART:MOLKART:STARDIST (mem_only)                 [100%] 1 of 1 ✔
    [3e/3c0736] NFCORE_MOLKART:MOLKART:CELLPOSE (mem_only)                 [100%] 1 of 1 ✔
    [a0/415c6a] NFCORE_MOLKART:MOLKART:MASKFILTER (mem_only)               [100%] 3 of 3 ✔
    [14/a830c9] NFCORE_MOLKART:MOLKART:SPOT2CELL (mem_only)                [100%] 3 of 3 ✔
    [b5/391836] NFCORE_MOLKART:MOLKART:CREATE_ANNDATA (mem_only)           [100%] 3 of 3 ✔
    [77/aed558] NFCORE_MOLKART:MOLKART:MOLKARTQC (mem_only)                [100%] 3 of 3 ✔
    [e6/b81475] NFCORE_MOLKART:MOLKART:MULTIQC                             [100%] 1 of 1 ✔
    -[nf-core/molkart] Pipeline completed successfully-
    Completed at: 19-Oct-2025 22:23:01
    Duration    : 2m 52s
    CPU hours   : 0.1
    Succeeded   : 22
    ```

Beachte, wie diese Ausgabe detaillierter ist als unser Hello-World-Beispiel, weil die Pipeline den nf-core-Konventionen folgt:

- Die Pipeline zeigt ihre Version und ihr Logo
- Konfigurationsparameter werden angezeigt
- Mehrere Prozesse laufen parallel (angezeigt durch mehrere Prozesszeilen)
- Prozessnamen enthalten den vollständigen Modulpfad (z.B. `NFCORE_MOLKART:MOLKART:MINDAGAP_MINDAGAP`)

### 2.4. Prozessausführung verstehen

Die Executor-Zeile `executor > local (22)` sagt dir:

- **executor**: Welche Rechenumgebung verwendet wird (`local` = dein Rechner)
- **(22)**: Gesamtzahl der gestarteten Aufgaben

Jede Prozesszeile zeigt:

- **Hash** (`[1a/2b3c4d]`): Work-Verzeichnis-Identifikator (wie zuvor)
- **Prozessname**: Vollständiger Modulpfad und Prozessname
- **Eingabe-Identifikator**: Sample-Name in Klammern
- **Fortschritt**: Prozentsatz abgeschlossen und Anzahl (z.B. `1 of 1 ✔`)

### Fazit

Du weißt, wie man eine nf-core-Pipeline mit Testdaten startet und ihre Ausführungsausgabe interpretiert.

### Wie geht es weiter?

Lerne, wo du die Ergebnisse findest und wie du sie interpretierst.

---

## 3. Ausgaben finden und untersuchen

Wenn die Pipeline erfolgreich abgeschlossen ist, siehst du eine Abschlussmeldung und eine Ausführungszusammenfassung.

### 3.1. Ergebnisverzeichnis finden

Standardmäßig schreiben nf-core-Pipelines Ausgaben in ein Verzeichnis, das durch den Parameter `outdir` angegeben wird, den wir auf `results/` gesetzt haben.

Liste den Inhalt auf:

```bash
tree results/
```

Du solltest mehrere Unterverzeichnisse sehen:

```console title="results/"
results/
├── anndata/
├── clahe/
├── mindagap/
├── molkartqc/
├── multiqc/
├── pipeline_info/
├── segmentation/
├── spot2cell/
└── stack/
```

Jedes Unterverzeichnis enthält Ausgaben aus einem bestimmten Stadium der Pipeline:

- **mindagap/**: Gitteraufgefüllte Bilder aus dem MindaGap-Vorverarbeitungsschritt
- **clahe/**: Kontrastverbesserte Bilder aus der CLAHE-Vorverarbeitung
- **stack/**: Mehrkanal-Bildstapel, die für die Segmentierung erstellt wurden
- **segmentation/**: Segmentierungsergebnisse verschiedener Algorithmen (cellpose/, mesmer/, stardist/, filtered_masks/)
- **spot2cell/**: Zell-nach-Transkript-Zähltabellen
- **anndata/**: AnnData-Objekte mit Zell-nach-Transkript-Matrizen und räumlichen Koordinaten
- **molkartqc/**: Qualitätskontrollmetriken für Spot-Zuordnung
- **multiqc/**: Umfassender Qualitätskontrollbericht
- **pipeline_info/**: Ausführungsberichte und Logs

### 3.2. MultiQC-Bericht untersuchen

Der MultiQC-Bericht ist eine umfassende HTML-Datei, die Qualitätsmetriken aus allen Pipeline-Schritten aggregiert.

Öffne den Bericht im Dateibrowser und klicke dann auf die Schaltfläche "Show Preview", um ihn direkt in VS Code gerendert zu sehen.

Der Bericht enthält:

- Allgemeine Statistiken für alle Samples
- Vorverarbeitungsmetriken
- Segmentierungsqualitätsmetriken
- Anzahl der erkannten Zellen und Spots

!!! Tip

    MultiQC-Berichte sind typischerweise in allen nf-core-Pipelines enthalten.
    Sie bieten immer einen Überblick über die Pipeline-Ausführung und Datenqualität.

### 3.3. Zell-nach-Transkript-Tabellen untersuchen

Die wichtigste wissenschaftliche Ausgabe ist die Zell-nach-Transkript-Zähltabelle.
Diese sagt dir, wie viele von jedem Transkript in jeder Zelle erkannt wurden.

Navigiere zum spot2cell-Verzeichnis:

```bash
ls results/spot2cell/
```

Du findest Dateien wie:

- `cellxgene_mem_only_cellpose.csv`: Zell-nach-Transkript-Tabelle mit Cellpose-Segmentierung
- `cellxgene_mem_only_mesmer.csv`: Zell-nach-Transkript-Tabelle mit Mesmer-Segmentierung
- `cellxgene_mem_only_stardist.csv`: Zell-nach-Transkript-Tabelle mit Stardist-Segmentierung

Wir haben in diesem Testdatensatz nur 1 Sample ausgeführt, aber in einem echten Experiment hätten wir diese Tabellen für jedes Sample.
Beachte, wie Nextflow mehrere Segmentierungsmethoden parallel verarbeiten kann, was den Vergleich von Ergebnissen erleichtert.

### 3.4. Ausführungsberichte ansehen

Nextflow generiert automatisch mehrere Ausführungsberichte.

Überprüfe das pipeline_info-Verzeichnis:

```bash
ls results/pipeline_info/
```

Wichtige Dateien:

- **execution_report.html**: Timeline und Ressourcennutzungsvisualisierung
- **execution_timeline.html**: Gantt-Diagramm der Prozessausführung
- **execution_trace.txt**: Detaillierte Aufgabenausführungsmetriken
- **pipeline_dag.html**: Gerichteter azyklischer Graph, der die Workflow-Struktur zeigt

Öffne den Ausführungsbericht, um die Ressourcennutzung zu sehen:

```bash
code results/pipeline_info/execution_report.html
```

Dies zeigt:

- Wie lange jeder Prozess gedauert hat
- CPU- und Speichernutzung
- Welche Aufgaben gecacht vs. ausgeführt wurden

!!! Tip

    Diese Berichte sind unglaublich nützlich zur Optimierung der Ressourcenzuweisung und zur Fehlerbehebung bei Leistungsproblemen.

### Fazit

Du weißt, wie man Pipeline-Ausgaben findet, Qualitätskontrollberichte untersucht und auf Ausführungsmetriken zugreift.

### Wie geht es weiter?

Lerne über das Work-Verzeichnis und wie Nextflow Zwischendateien verwaltet.

---

## 4. Das Work-Verzeichnis erkunden

Genau wie bei unserem Hello-World-Beispiel findet die gesamte eigentliche Arbeit im `work/`-Verzeichnis statt.

### 4.1. Work-Verzeichnisstruktur verstehen

Das Work-Verzeichnis enthält ein Unterverzeichnis für jede ausgeführte Aufgabe.
Für diese Pipeline mit 12 Aufgaben gibt es 12 Work-Unterverzeichnisse.

Liste das Work-Verzeichnis auf:

```bash
ls -d work/*/*/ | head -5
```

Dies zeigt die ersten 5 Aufgabenverzeichnisse.

### 4.2. Ein Aufgabenverzeichnis inspizieren

Wähle einen der Segmentierungsprozess-Hashes aus der Konsolenausgabe (z.B. `[3m/4n5o6p]`) und schau hinein:

```bash
ls -la work/3m/4n5o6p*/
```

Du siehst:

- **.command.\*-Dateien**: Nextflow-Ausführungsskripte und Logs (wie zuvor)
- **Bereitgestellte Eingabedateien**: Symlinks zu den tatsächlichen Eingabedateien
- **Ausgabedateien**: Segmentierungsmasken, Zwischenergebnisse usw.

Der Hauptunterschied zu Hello World:

- Echte Pipelines stellen große Eingabedateien bereit (Bilder, Referenzdaten)
- Ausgabedateien können recht groß sein (Segmentierungsmasken, verarbeitete Bilder)
- Mehrere Eingabe- und Ausgabedateien pro Aufgabe

!!! Tip

    Wenn ein Prozess fehlschlägt, kannst du zu seinem Work-Verzeichnis navigieren, `.command.err` auf Fehlermeldungen untersuchen und sogar `.command.sh` manuell erneut ausführen, um das Problem zu debuggen.

### 4.3. Work-Verzeichnis aufräumen

Das Work-Verzeichnis kann über mehrere Pipeline-Ausführungen hinweg recht groß werden.
Wie wir in Teil 1 gelernt haben, kannst du `nextflow clean` verwenden, um Work-Verzeichnisse von alten Ausführungen zu entfernen.

Für nf-core-Pipelines mit großen Zwischendateien ist es jedoch besonders wichtig, regelmäßig aufzuräumen.

### Fazit

Du verstehst, wie nf-core-Pipelines ihre Work-Verzeichnisse organisieren und wie man einzelne Aufgaben zum Debuggen inspiziert.

### Wie geht es weiter?

Lerne über den Nextflow-Cache und wie man fehlgeschlagene Pipeline-Ausführungen fortsetzt.

---

## 5. Eine Pipeline-Ausführung fortsetzen

Eine der mächtigsten Funktionen von Nextflow ist die Fähigkeit, eine Pipeline vom Punkt des Fehlers aus fortzusetzen.

### 5.1. Der Cache-Mechanismus

Wenn du eine Pipeline mit `-resume` ausführst, macht Nextflow Folgendes:

1. Überprüft den Cache für jede Aufgabe
2. Wenn Eingaben, Code und Parameter identisch sind, wird das gecachte Ergebnis wiederverwendet
3. Führt nur Aufgaben erneut aus, die sich geändert haben oder fehlgeschlagen sind

Dies ist unverzichtbar für lang laufende Pipelines, bei denen Fehler spät in der Ausführung auftreten können.

### 5.2. Resume mit molkart ausprobieren

Führe denselben Befehl erneut aus, aber füge `-resume` hinzu:

```bash
nextflow run ./molkart \
  --input 'data/samplesheet.csv' \
  --mindagap_tilesize 90 \
  --mindagap_boxsize 7 \
  --mindagap_loopnum 100 \
  --clahe_pyramid_tile 368 \
  --segmentation_method "cellpose" \
  --outdir results \
  -resume
```

Du solltest eine Ausgabe wie diese sehen: <!-- TODO: vollständige Ausgabe -->

```console
executor >  local (0)
[1a/2b3c4d] NFCORE_MOLKART:MOLKART:MINDAGAP_MINDAGAP (mem_only)        [100%] 2 of 2, cached: 2 ✔
[5e/6f7g8h] NFCORE_MOLKART:MOLKART:CLAHE (mem_only)                    [100%] 2 of 2, cached: 2 ✔
[7f/8g9h0i] NFCORE_MOLKART:MOLKART:CREATE_STACK (mem_only)             [100%] 1 of 1, cached: 1 ✔
[9h/0i1j2k] NFCORE_MOLKART:MOLKART:MINDAGAP_DUPLICATEFINDER (mem_only) [100%] 1 of 1, cached: 1 ✔
[2k/3l4m5n] NFCORE_MOLKART:MOLKART:CELLPOSE (mem_only)                 [100%] 1 of 1, cached: 1 ✔
...
```

Beachte `cached: 2` oder `cached: 1` für jeden Prozess - nichts wurde erneut ausgeführt!

### 5.3. Wann Resume nützlich ist

Resume ist besonders wertvoll, wenn:

- Eine Pipeline aufgrund von Ressourcenlimits fehlschlägt (Speicher voll, Zeitlimit überschritten)
- Du nachgelagerte Prozesse ändern musst, ohne vorgelagerte Schritte erneut auszuführen
- Deine Netzwerkverbindung während des Datendownloads abbricht
- Du zusätzliche Ausgaben hinzufügen möchtest, ohne die Berechnung zu wiederholen

!!! Warning

    Resume funktioniert nur, wenn du die Eingabedaten, den Pipeline-Code oder die Parameter nicht geändert hast.
    Wenn du eines davon änderst, wird Nextflow betroffene Aufgaben korrekt erneut ausführen.

### Fazit

Du weißt, wie man `-resume` verwendet, um Pipelines effizient erneut auszuführen, ohne erfolgreiche Aufgaben zu wiederholen.

### Wie geht es weiter?

Jetzt, da du nf-core/molkart mit Testdaten ausführen kannst, bist du bereit zu lernen, wie du sie für deine eigenen Datensätze konfigurierst.
