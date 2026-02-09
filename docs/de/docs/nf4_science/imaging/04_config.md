# Teil 4: Konfiguration

In den Teilen 1-3 haben wir gelernt, wie man Nextflow ausführt, eine nf-core-Pipeline verwendet und Eingaben mit Parameterdateien und Samplesheets verwaltet.
Jetzt werden wir erkunden, wie man Pipelines für verschiedene Rechenumgebungen mit **Konfigurationsdateien** und **Profilen** konfiguriert.

## Lernziele

Am Ende dieses Teils wirst du in der Lage sein:

- Zu verstehen, wie Nextflow Konfigurationen aus mehreren Quellen auflöst
- Integrierte nf-core-Profile für Container und Tests zu verwenden
- Eigene Profile für verschiedene Rechenumgebungen zu erstellen
- Ressourcenanforderungen mithilfe von Prozess-Labels anzupassen
- Ressourcenlimits in eingeschränkten Umgebungen zu verwalten
- Aufgelöste Konfigurationen mit `nextflow config` zu inspizieren

---

## 1. Nextflow-Konfiguration verstehen

### 1.1. Was ist eine Konfigurationsdatei?

Nextflow verwendet Konfigurationsdateien, um **Workflow-Logik** (was zu tun ist) von **Ausführungseinstellungen** (wie und wo es zu tun ist) zu trennen.

Konfigurationsdateien steuern:

- Container-Engines (Docker, Singularity, Conda)
- Rechenressourcen (CPUs, Arbeitsspeicher, Zeit)
- Ausführungsplattformen (lokal, HPC, Cloud)
- Pipeline-Parameter

### 1.2. Konfigurationspriorität

Nextflow lädt Konfigurationen aus mehreren Quellen, wobei spätere Quellen frühere überschreiben:

1. **Pipeline-Konfiguration**: `nextflow.config` im Pipeline-Repository
2. **Verzeichnis-Konfiguration**: `nextflow.config` in deinem aktuellen Arbeitsverzeichnis
3. **Benutzer-Konfiguration**: `~/.nextflow/config`
4. **Kommandozeile**: Direkt übergebene Parameter und Optionen

Dieser mehrschichtige Ansatz ermöglicht es dir, Standardwerte in der Pipeline zu behalten, sie mit benutzerspezifischen Einstellungen zu überschreiben und schnelle Anpassungen über die Kommandozeile vorzunehmen.

### 1.3. Unsere aktuelle Konfiguration

Schauen wir uns die Konfiguration an, die wir bisher verwendet haben:

```groovy title="nextflow.config"
docker.enabled = true
process {
    resourceLimits = [
        cpus: 2,
        memory: '7.GB',
    ]
}

```

Kommentiere die Zeile `docker.enabled = true` aus Teil 2 aus oder ändere sie zurück, und finde heraus, wie wir das gleiche Ergebnis mit einem Profil in molkart erreichen können.

---

## 2. Profile verwenden

### 2.1. Was sind Profile?

Profile sind benannte Konfigurationssets, die mit dem Flag `-profile` über den Befehl `nextflow run` aktiviert werden können.
Sie erleichtern den Wechsel zwischen verschiedenen Rechenszenarien, ohne Konfigurationsdateien bearbeiten zu müssen.

Alle nf-core-Pipelines enthalten eine Reihe von Standardprofilen, die wir nutzen können.

### 2.2. Integrierte Profile inspizieren

Schauen wir sie uns in der Datei `molkart/nextflow.config` an, die mit der Pipeline-Codebasis verknüpft ist:

```bash
code molkart/nextflow.config
```

Suche nach dem `profiles`-Block:

```groovy title="molkart/nextflow.config (excerpt)"
profiles {
    docker {
        docker.enabled          = true
        singularity.enabled     = false
        conda.enabled           = false
    }
    singularity {
        singularity.enabled     = true
        singularity.autoMounts  = true
        docker.enabled          = false
    }
    conda {
        conda.enabled           = true
        docker.enabled          = false
        conda.channels          = ['conda-forge', 'bioconda']
    }
}
```

Gängige Container-Profile:

- `docker`: Docker-Container verwenden (am häufigsten für lokale Entwicklung)
- `singularity`: Singularity/Apptainer verwenden (häufig auf HPC)
- `conda`: Conda-Umgebungen verwenden
- `apptainer`: Apptainer-Container verwenden

### 2.3. Erneute Ausführung mit Profilen statt nextflow.config

Nachdem wir die Docker-Konfiguration in unserer lokalen `nextflow.config`-Datei deaktiviert haben und Profile verstehen, führen wir die Pipeline nun mit dem Flag `-profile` erneut aus.

Zuvor in Teil 3 haben wir eine `params.yaml`-Datei mit unseren benutzerdefinierten Parametern erstellt.
Wir können diese nun mit dem integrierten Docker-Profil kombinieren:

```bash
nextflow run ./molkart \
  -profile docker \
  -params-file params.yaml \
  -resume
```

Schauen wir uns an, was jedes Flag bewirkt:

- `-profile docker`: Aktiviert das Docker-Profil aus molkarts `nextflow.config`, das `docker.enabled = true` setzt
- `-params-file params.yaml`: Lädt alle Pipeline-Parameter aus unserer YAML-Datei
- `-resume`: Verwendet zwischengespeicherte Ergebnisse aus vorherigen Ausführungen wieder

Da wir `-resume` verwenden, prüft Nextflow, ob sich seit der letzten Ausführung etwas geändert hat.
Wenn Parameter, Eingaben und Code gleich sind, werden alle Aufgaben aus dem Cache abgerufen und die Pipeline wird fast sofort abgeschlossen.

```console title="Output (excerpt)"
executor >  local (12)
...
[1a/2b3c4d] NFCORE_MOLKART:MOLKART:MINDAGAP_MINDAGAP (mem_only)   [100%] 2 of 2, cached: 2 ✔
[5e/6f7g8h] NFCORE_MOLKART:MOLKART:CLAHE (mem_only)               [100%] 2 of 2, cached: 2 ✔
...
-[nf-core/molkart] Pipeline completed successfully-
```

Beachte, dass alle Prozesse `cached: 2` oder `cached: 1` anzeigen - nichts wurde erneut ausgeführt!

### 2.4. Test-Profile

Test-Profile bieten schnelle Möglichkeiten, Standard-Eingabeparameter und Datendateien anzugeben, damit du überprüfen kannst, ob die Pipeline funktioniert.
nf-core-Pipelines enthalten immer mindestens zwei Test-Profile:

- `test`: Kleiner Datensatz mit schnellen Parametern für schnelle Tests
- `test_full`: Umfassenderer Test mit größeren Daten

Schauen wir uns das `test`-Profil in molkart genauer an, das mit der `includeConfig`-Direktive eingebunden wird:

```groovy title="molkart/nextflow.config (excerpt)"
profiles {
  ...
    test      { includeConfig 'conf/test.config'      }
}
```

Das bedeutet, dass Nextflow immer dann, wenn wir die Pipeline mit `-profile test` ausführen, die Konfiguration aus `conf/test.config` lädt.

```groovy title="molkart/conf/test.config (excerpt)"
params {
    config_profile_name        = 'Test profile'
    config_profile_description = 'Minimal test dataset to check pipeline function'

    input = 'https://raw.githubusercontent.com/nf-core/test-datasets/molkart/test_data/samplesheets/samplesheet_membrane.csv'
    mindagap_tilesize = 90
    mindagap_boxsize = 7
    mindagap_loopnum = 100
    clahe_pyramid_tile = 368
    segmentation_method = "mesmer,cellpose,stardist"
}

process {
    resourceLimits = [
        cpus: 4,
        memory: '15.GB',
        time: '1.h'
    ]
}
```

Beachte, dass dieses Profil die gleichen Parameter enthält, die wir zuvor in unserer `params.yaml`-Datei verwendet haben.

Du kannst mehrere Profile aktivieren, indem du sie durch Kommas trennst.
Nutzen wir das, um unsere Pipeline ohne unsere Parameterdatei zu testen:

```bash
nextflow run ./molkart -profile docker,test --outdir results -resume
```

Dies kombiniert:

- `docker`: Docker-Container aktivieren
- `test`: Test-Datensatz und -Parameter verwenden

Profile werden von links nach rechts angewendet, sodass spätere Profile frühere überschreiben, wenn sie die gleichen Werte setzen.

### Fazit

nf-core-Pipelines enthalten integrierte Profile für Container, Tests und spezielle Umgebungen.
Du kannst mehrere Profile kombinieren, um die benötigte Konfiguration zu erstellen.

### Wie geht es weiter?

Lerne, wie du eigene Profile für verschiedene Rechenumgebungen erstellst.

---

## 3. Eigene Profile erstellen

### 3.1. Profile für den Wechsel zwischen lokaler Entwicklung und Ausführung auf HPC erstellen

Erstellen wir eigene Profile für zwei Szenarien:

1. Lokale Entwicklung mit Docker
2. Universitäts-HPC mit Slurm-Scheduler und Singularity

Füge Folgendes zu deiner `nextflow.config` hinzu:

```groovy title="nextflow.config"
profiles {
    local_dev {
        docker.enabled = true
        process.executor = 'local'
    }

    hpc_cluster {
        singularity.enabled = true
        process.executor = 'slurm'
        process.queue = 'standard_queue'
        singularity.cacheDir = '/shared/containers'
    }
}
```

Jetzt kannst du einfach zwischen Umgebungen wechseln:

```bash
# Für lokale Entwicklung
nextflow run ./molkart -profile local_dev --input data/samplesheet.csv --outdir results

# Für HPC (wenn verfügbar)
nextflow run ./molkart -profile hpc_cluster --input data/samplesheet.csv --outdir results
```

!!! note "Hinweis"

    Wir können das HPC-Profil in dieser Trainingsumgebung nicht testen, da wir keinen Zugriff auf einen Slurm-Scheduler haben.
    Aber dies zeigt, wie du es für den realen Einsatz konfigurieren würdest.

### 3.2. `nextflow config` verwenden, um Konfigurationen zu inspizieren

Der Befehl `nextflow config` zeigt die vollständig aufgelöste Konfiguration an, ohne die Pipeline auszuführen.

Standardkonfiguration anzeigen:

```bash
nextflow config ./molkart
```

Konfiguration mit einem bestimmten Profil anzeigen:

```bash
nextflow config -profile local_dev ./molkart
```

Dies ist äußerst nützlich für:

- Debugging von Konfigurationsproblemen
- Verstehen, welche Werte tatsächlich verwendet werden
- Überprüfen, wie mehrere Profile interagieren

### Fazit

Eigene Profile ermöglichen es dir, mit einem einzigen Kommandozeilen-Flag zwischen verschiedenen Rechenumgebungen zu wechseln.
Verwende `nextflow config`, um die aufgelöste Konfiguration vor der Ausführung zu inspizieren.

### Wie geht es weiter?

Lerne, wie du Ressourcenanforderungen für einzelne Prozesse mithilfe des Prozess-Label-Systems von nf-core anpasst.

---

## 4. Ressourcenanforderungen anpassen

### 4.1. Prozess-Labels in nf-core-Pipelines verstehen

Der Einfachheit halber verwenden nf-core-Pipelines [**Prozess-Labels**](https://www.nextflow.io/docs/latest/reference/process.html#process-label), um die Ressourcenzuweisung über alle Pipelines hinweg zu standardisieren.
Jeder Prozess wird mit einem Label wie `process_low`, `process_medium` oder `process_high` versehen, um niedrige, mittlere bzw. hohe Rechenressourcenanforderungen zu beschreiben.
Diese Labels werden in einer der Konfigurationsdateien im Verzeichnis `conf/` der Pipeline in spezifische Ressourcenanforderungen umgewandelt.

```groovy title="molkart/conf/base.config (excerpt)"
process {
    cpus   = { 1      * task.attempt }
    memory = { 6.GB   * task.attempt }
    time   = { 4.h    * task.attempt }

    errorStrategy = { task.exitStatus in ((130..145) + 104) ? 'retry' : 'finish' }
    maxRetries    = 1
    maxErrors     = '-1'

    withLabel:process_single {
        cpus   = { 1                   }
        memory = { 6.GB * task.attempt }
        time   = { 4.h  * task.attempt }
    }
    withLabel:process_low {
        cpus   = { 2     * task.attempt }
        memory = { 12.GB * task.attempt }
        time   = { 4.h   * task.attempt }
    }
    withLabel:process_medium {
        cpus   = { 6     * task.attempt }
        memory = { 36.GB * task.attempt }
        time   = { 8.h   * task.attempt }
    }
    withLabel:process_high {
        cpus   = { 12    * task.attempt }
        memory = { 72.GB * task.attempt }
        time   = { 16.h  * task.attempt }
    }
}
```

Beachte den Multiplikator `task.attempt` - dieser ermöglicht es nachfolgenden Aufgabenwiederholungen, mehr Ressourcen anzufordern, wenn die Pipeline mit `process.maxRetries > 1` konfiguriert ist.

### 4.2. Ressourcen für bestimmte Prozesse überschreiben

Für eine feinkörnige Kontrolle kannst du einzelne Prozesse nach Namen ansprechen:

```groovy title="nextflow.config"
process {
    withName: 'NFCORE_MOLKART:MOLKART:CELLPOSE' {
        cpus   = 16
        memory = 32.GB
    }
}
```

Wenn wir versuchen, diese Pipeline mit der obigen Überschreibung auszuführen, wird der `CELLPOSE`-Prozess 16 CPUs und 32 GB Arbeitsspeicher anfordern statt der durch sein Label definierten Standardwerte.
Dies führt in unserer aktuellen Umgebung zu einem Fehler, da wir nicht so viel RAM zur Verfügung haben.
Im nächsten Abschnitt lernen wir, wie man solche Fehler verhindert.

!!! tip "Tipp"

    Um Prozessnamen zu finden, schaue dir die Pipeline-Ausführungsausgabe an oder prüfe `.nextflow.log`.
    Prozessnamen folgen dem Muster `WORKFLOW:SUBWORKFLOW:PROCESS`.

### Fazit

nf-core-Pipelines verwenden Prozess-Labels, um die Ressourcenzuweisung zu standardisieren.
Du kannst Ressourcen nach Label (betrifft mehrere Prozesse) oder nach Namen (betrifft einen bestimmten Prozess) überschreiben.

### Wie geht es weiter?

Lerne, wie du Ressourcenlimits in eingeschränkten Umgebungen wie GitHub Codespaces verwaltest.

---

## 5. Ressourcen in eingeschränkten Umgebungen verwalten

### 5.1. Das Problem der Ressourcenlimits

Wenn wir versuchen würden, molkart mit einem Prozess auszuführen, der 16 CPUs und 32 GB Arbeitsspeicher anfordert (wie in Abschnitt 4.2 gezeigt), würde dies in unserer aktuellen Umgebung fehlschlagen, da wir nicht so viele Ressourcen zur Verfügung haben.
In einer Cluster-Umgebung mit größeren Knoten würden solche Anforderungen an den Scheduler übermittelt.

In eingeschränkten Umgebungen wie GitHub Codespaces würde Nextflow ohne Limits die Ausführung von Prozessen verweigern, die die verfügbaren Ressourcen überschreiten.

### 5.2. Ressourcenlimits festlegen

Die Direktive `resourceLimits` begrenzt Ressourcenanforderungen auf angegebene Werte:

```groovy title="nextflow.config"
process {
    resourceLimits = [ cpus: 2, memory: 7.GB ]
}
```

Dies teilt Nextflow mit: "Wenn ein Prozess mehr als 2 CPUs oder 7 GB Arbeitsspeicher anfordert, begrenze es stattdessen auf diese Limits."

### 5.3. Ressourcenlimits zu eigenen Profilen hinzufügen

Aktualisiere deine eigenen Profile, um geeignete Limits einzuschließen:

```groovy title="nextflow.config"
profiles {
    local_dev {
        docker.enabled = true
        process.executor = 'local'
        process.resourceLimits = [
            cpus: 2,
            memory: 7.GB
        ]
    }

    hpc_cluster {
        singularity.enabled = true
        process.executor = 'slurm'
        process.queue = 'batch'
        process.resourceLimits = [
            cpus: 32,
            memory: 128.GB,
            time: 24.h
        ]
    }
}
```

!!! warning "Warnung"

    Zu niedrige Ressourcenlimits können dazu führen, dass Prozesse fehlschlagen oder langsam laufen.
    Die Pipeline muss möglicherweise weniger speicherintensive Algorithmen verwenden oder Daten in kleineren Blöcken verarbeiten.

### Fazit

Verwende `resourceLimits`, um Pipelines in ressourcenbeschränkten Umgebungen auszuführen, indem du Prozessressourcenanforderungen begrenzt.
Verschiedene Profile können unterschiedliche Limits haben, die für ihre Umgebung geeignet sind.

### Wie geht es weiter?

Du hast das Kerntraining Nextflow für Bioimaging abgeschlossen!

---

## Zusammenfassung

Du verstehst jetzt, wie man Nextflow-Pipelines für verschiedene Rechenumgebungen konfiguriert.

Wichtige Fähigkeiten, die du gelernt hast:

- **Konfigurationspriorität**: Wie Nextflow Einstellungen aus mehreren Quellen auflöst
- **nf-core-Profile**: Verwendung integrierter Profile für Container, Tests und Hilfsprogramme
- **Eigene Profile**: Erstellen eigener Profile für verschiedene Umgebungen
- **Prozess-Labels**: Verstehen und Überschreiben von Ressourcenanforderungen nach Label
- **Ressourcenlimits**: Verwaltung eingeschränkter Umgebungen mit `resourceLimits`
- **Konfigurationsinspektion**: Verwendung von `nextflow config` zum Debuggen und Überprüfen von Einstellungen

Diese Konfigurationsfähigkeiten sind auf jede Nextflow-Pipeline übertragbar und helfen dir, Workflows effizient auf lokalen Rechnern, HPC-Clustern und Cloud-Plattformen auszuführen.

### Wie geht es weiter?

Herzlichen Glückwunsch zum Abschluss des Kurses Nextflow für Bioimaging!

Nächste Schritte:

- Fülle die Kursumfrage aus, um Feedback zu geben
- Schau dir [Hello Nextflow](../hello_nextflow/index.md) an, um mehr über die Entwicklung von Workflows zu erfahren
- Erkunde [Hello nf-core](../hello_nf-core/index.md), um tiefer in die nf-core-Tools einzutauchen
- Durchsuche weitere Kurse in den [Trainingssammlungen](../training_collections/index.md)
