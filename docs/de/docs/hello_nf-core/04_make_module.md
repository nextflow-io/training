# Teil 4: Ein nf-core-Modul erstellen

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

In diesem vierten Teil des Hello nf-core Trainingskurses zeigen wir dir, wie du ein nf-core-Modul erstellst, indem du die wichtigsten Konventionen anwendest, die Module portabel und wartbar machen.

Das nf-core-Projekt stellt einen Befehl (`nf-core modules create`) bereit, der automatisch korrekt strukturierte Modul-Templates generiert, ähnlich wie wir es in Teil 2 für den Workflow verwendet haben.
Zu Lehrzwecken werden wir jedoch zunächst manuell vorgehen: Wir transformieren das lokale `cowpy`-Modul in deiner `core-hello`-Pipeline Schritt für Schritt in ein nf-core-Modul.
Danach zeigen wir dir, wie du die template-basierte Modul-Erstellung nutzen kannst, um in Zukunft effizienter zu arbeiten.

??? info "Wie du von diesem Abschnitt aus beginnst"

    Dieser Abschnitt setzt voraus, dass du [Teil 3: Ein nf-core-Modul verwenden](./03_use_module.md) abgeschlossen hast und das `CAT_CAT`-Modul in deine Pipeline integriert hast.

    Falls du Teil 3 nicht abgeschlossen hast oder für diesen Teil neu beginnen möchtest, kannst du die `core-hello-part3`-Lösung als Ausgangspunkt verwenden.
    Führe diese Befehle im `hello-nf-core/`-Verzeichnis aus:

    ```bash
    cp -r solutions/core-hello-part3 core-hello
    cd core-hello
    ```

    Dies gibt dir eine Pipeline mit dem bereits integrierten `CAT_CAT`-Modul.
    Du kannst testen, dass sie erfolgreich läuft, indem du folgenden Befehl ausführst:

    ```bash
    nextflow run . --outdir core-hello-results -profile test,docker --validate_params false
    ```

---

## 1. `cowpy` in ein nf-core-Modul transformieren

In diesem Abschnitt wenden wir nf-core-Konventionen auf das lokale `cowpy`-Modul in deiner `core-hello`-Pipeline an und transformieren es in ein Modul, das den nf-core-Community-Standards folgt.

Dies ist der aktuelle Code für das `cowpy`-Prozessmodul:

```groovy title="core-hello/modules/local/cowpy.nf" linenums="1"
#!/usr/bin/env nextflow

// ASCII-Art mit cowpy generieren (https://github.com/jeffbuttars/cowpy)
process cowpy {

    publishDir 'results', mode: 'copy'

    container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
    conda 'conda-forge::cowpy==1.1.5'

    input:
    path input_file
    val character

    output:
    path "cowpy-${input_file}"

    script:
    """
    cat $input_file | cowpy -c "$character" > cowpy-${input_file}
    """
}
```

Wir werden schrittweise die folgenden nf-core-Konventionen anwenden:

1. **Den Prozessnamen in `COWPY` großschreiben**, um der Konvention zu folgen.
2. **`COWPY` aktualisieren, um Metadaten-Tupel zu verwenden**, um Sample-Metadaten durch den Workflow zu propagieren.
3. **Tool-Argument-Konfiguration mit `ext.args` zentralisieren**, um die Vielseitigkeit des Moduls zu erhöhen und gleichzeitig die Schnittstelle minimal zu halten.
4. **Ausgabe-Benennung mit `ext.prefix` standardisieren**, um Konsistenz zu fördern.
5. **Die Publishing-Konfiguration zentralisieren**, um Konsistenz zu fördern.

Nach jedem Schritt führen wir die Pipeline aus, um zu testen, dass alles wie erwartet funktioniert.

!!! warning "Arbeitsverzeichnis"

    Stelle sicher, dass du dich im `core-hello`-Verzeichnis (deinem Pipeline-Root) befindest für alle Dateibearbeitungen und Befehlsausführungen in diesem Abschnitt.

    ```bash
    cd core-hello
    ```

### 1.1. Den Prozessnamen großschreiben

Dies ist rein eine stilistische Konvention (es gibt keine technische Rechtfertigung), aber da es die Norm für nf-core-Module ist, sollten wir uns daran halten.

Wir müssen drei Änderungen vornehmen:

1. Den Prozessnamen im Modul aktualisieren
2. Die Modul-Import-Anweisung im Workflow-Header aktualisieren
3. Den Prozessaufruf und die emit-Deklaration im Workflow-Body aktualisieren

Los geht's!

#### 1.1.1. Den Prozessnamen im Modul aktualisieren

Öffne die `cowpy.nf`-Moduldatei (unter `core-hello/modules/local/`) und ändere den Prozessnamen in Großbuchstaben:

=== "Danach"

    ```groovy title="core-hello/modules/local/cowpy.nf" linenums="3" hl_lines="2"
    // ASCII-Art mit cowpy generieren (https://github.com/jeffbuttars/cowpy)
    process COWPY {
    ```

=== "Vorher"

    ```groovy title="core-hello/modules/local/cowpy.nf" linenums="3" hl_lines="2"
    // ASCII-Art mit cowpy generieren (https://github.com/jeffbuttars/cowpy)
    process cowpy {
    ```

In diesem Fall ist die Großschreibung völlig unkompliziert.

Wenn der Prozessname aus mehreren Wörtern bestünde, zum Beispiel wenn wir einen Prozess namens MyCowpyTool ursprünglich in Camel Case hätten, wäre die nf-core-Konvention, Unterstriche zur Trennung zu verwenden, was MY_COWPY_TOOL ergäbe.

#### 1.1.2. Die Modul-Import-Anweisung aktualisieren

Prozessnamen sind case-sensitiv, daher müssen wir jetzt, da wir den Prozessnamen geändert haben, die Modul-Import-Anweisung entsprechend im Workflow-Header von `hello.nf` aktualisieren:

=== "Danach"

    ```groovy title="core-hello/workflows/hello.nf" linenums="1" hl_lines="10"
    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        MODULE / SUBWORKFLOWS / FUNKTIONEN IMPORTIEREN
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    include { paramsSummaryMap       } from 'plugin/nf-schema'
    include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
    include { sayHello               } from '../modules/local/sayHello.nf'
    include { convertToUpper         } from '../modules/local/convertToUpper.nf'
    include { COWPY                  } from '../modules/local/cowpy/main.nf'
    include { CAT_CAT                } from '../modules/nf-core/cat/cat/main'
    ```

=== "Vorher"

    ```groovy title="core-hello/workflows/hello.nf" linenums="1" hl_lines="10"
    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        MODULE / SUBWORKFLOWS / FUNKTIONEN IMPORTIEREN
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    include { paramsSummaryMap       } from 'plugin/nf-schema'
    include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
    include { sayHello               } from '../modules/local/sayHello.nf'
    include { convertToUpper         } from '../modules/local/convertToUpper.nf'
    include { cowpy                  } from '../modules/local/cowpy/main.nf'
    include { CAT_CAT                } from '../modules/nf-core/cat/cat/main'
    ```

Wir könnten einen Alias in der Import-Anweisung verwenden, um zu vermeiden, dass wir Aufrufe des Prozesses aktualisieren müssen, aber das würde den Sinn der Übernahme der Großschreibungskonvention etwas verfehlen.

#### 1.1.3. Den Prozessaufruf und die emit-Deklaration aktualisieren

Also aktualisieren wir jetzt die beiden Referenzen auf den Prozess im Workflow-Block von `hello.nf`:

=== "Danach"

    ```groovy title="core-hello/workflows/hello.nf" linenums="43" hl_lines="2 17"
    // ASCII-Art der Begrüßungen mit cowpy generieren
    COWPY(CAT_CAT.out.file_out)

    //
    // Software-Versionen zusammenstellen und speichern
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'hello_software_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    emit:
    cowpy_hellos   = COWPY.out.cowpy_output
    versions       = ch_versions
    ```

=== "Vorher"

    ```groovy title="core-hello/workflows/hello.nf" linenums="43" hl_lines="2 17"
    // ASCII-Art der Begrüßungen mit cowpy generieren
    cowpy(CAT_CAT.out.file_out)

    //
    // Software-Versionen zusammenstellen und speichern
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'hello_software_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    emit:
    cowpy_hellos   = cowpy.out.cowpy_output
    versions       = ch_versions
    ```

Stelle sicher, dass du **beide** Änderungen vornimmst, sonst erhältst du einen Fehler, wenn du dies ausführst.

#### 1.1.4. Die Pipeline ausführen, um sie zu testen

Lass uns den Workflow ausführen, um zu testen, dass nach diesen Änderungen alles korrekt funktioniert.

```bash
nextflow run . --outdir core-hello-results -profile test,docker --validate_params false
```

??? success "Befehlsausgabe"

    ```console
     N E X T F L O W   ~  version 25.04.3

    Launching `./main.nf` [elegant_plateau] DSL2 - revision: b9e9b3b8de

    Input/output options
      input                     : /workspaces/training/hello-nf-core/core-hello/assets/greetings.csv
      outdir                    : core-hello-results

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      validate_params           : false
      trace_report_suffix       : 2026-01-06_04-51-29

    Core Nextflow options
      runName                   : elegant_plateau
      containerEngine           : docker
      launchDir                 : /workspaces/training/hello-nf-core/core-hello
      workDir                   : /workspaces/training/hello-nf-core/core-hello/work
      projectDir                : /workspaces/training/hello-nf-core/core-hello
      userName                  : root
      profile                   : test,docker
      configFiles               : /workspaces/training/hello-nf-core/core-hello/nextflow.config

    !! Only displaying parameters that differ from the pipeline defaults !!
    ------------------------------------------------------
    executor >  local (8)
    [7b/66ceb5] CORE_HELLO:HELLO:sayHello (3)       | 3 of 3 ✔
    [8e/1bafb9] CORE_HELLO:HELLO:convertToUpper (3) | 3 of 3 ✔
    [bb/203575] CORE_HELLO:HELLO:CAT_CAT (test)     | 1 of 1 ✔
    [39/715489] CORE_HELLO:HELLO:COWPY              | 1 of 1 ✔
    -[core/hello] Pipeline completed successfully-
    ```

Gut, das funktioniert! Jetzt machen wir weiter mit substanzielleren Änderungen.

### 1.2. `COWPY` aktualisieren, um Metadaten-Tupel zu verwenden

In der aktuellen Version der `core-hello`-Pipeline extrahieren wir die Datei aus dem Output-Tupel von `CAT_CAT`, um sie an `COWPY` zu übergeben, wie in der oberen Hälfte des Diagramms unten gezeigt.

<figure class="excalidraw">
    --8<-- "docs/en/docs/hello_nf-core/img/cowpy-inputs.svg"
</figure>

Es wäre besser, wenn `COWPY` Metadaten-Tupel direkt akzeptieren würde, sodass Metadaten durch den Workflow fließen können, wie in der unteren Hälfte des Diagramms gezeigt.

Dazu müssen wir folgende Änderungen vornehmen:

1. Die Input- und Output-Definitionen aktualisieren
2. Den Prozessaufruf im Workflow aktualisieren
3. Den emit-Block im Workflow aktualisieren

Sobald wir das alles erledigt haben, führen wir die Pipeline aus, um zu testen, dass alles noch wie zuvor funktioniert.

#### 1.2.1. Die Input- und Output-Definitionen aktualisieren

Kehre zur `cowpy.nf`-Moduldatei zurück und modifiziere sie, um Metadaten-Tupel zu akzeptieren, wie unten gezeigt.

=== "Danach"

    ```groovy title="core-hello/modules/local/cowpy.nf" linenums="11" hl_lines="2 6"
        input:
            tuple val(meta), path(input_file)
            val character

        output:
            tuple val(meta), path("cowpy-${input_file}"), emit: cowpy_output
    ```

=== "Vorher"

    ```groovy title="core-hello/modules/local/cowpy.nf" linenums="11" hl_lines="2 6"
        input:
            path input_file
            val character

        output:
            path "cowpy-${input_file}"
    ```

Wie du sehen kannst, haben wir sowohl die **Haupteingabe** als auch die **Ausgabe** in ein Tupel geändert, das dem `tuple val(meta), path(input_file)`-Muster folgt, das in Teil 3 dieses Trainings eingeführt wurde.
Für die Ausgabe haben wir diese Gelegenheit auch genutzt, um `emit: cowpy_output` hinzuzufügen, um dem Ausgabekanal einen beschreibenden Namen zu geben.

Jetzt, da wir geändert haben, was der Prozess erwartet, müssen wir aktualisieren, was wir ihm im Prozessaufruf bereitstellen.

#### 1.2.2. Den Prozessaufruf im Workflow aktualisieren

Die gute Nachricht ist, dass diese Änderung den Prozessaufruf vereinfachen wird.
Jetzt, da die Ausgabe von `CAT_CAT` und die Eingabe von `COWPY` die gleiche 'Form' haben, d.h. beide aus einer `tuple val(meta), path(input_file)`-Struktur bestehen, können wir sie einfach direkt verbinden, anstatt die Datei explizit aus der Ausgabe des `CAT_CAT`-Prozesses extrahieren zu müssen.

Öffne die `hello.nf`-Workflow-Datei (unter `core-hello/workflows/`) und aktualisiere den Aufruf von `COWPY` wie unten gezeigt.

=== "Danach"

    ```groovy title="core-hello/workflows/hello.nf" linenums="43" hl_lines="2"
        // ASCII-Art der Begrüßungen mit cowpy generieren
        COWPY(CAT_CAT.out.file_out, params.character)
    ```

=== "Vorher"

    ```groovy title="core-hello/workflows/hello.nf" linenums="43" hl_lines="1-2 5"
        // Datei aus dem Tupel extrahieren, da cowpy noch keine Metadaten verwendet
        ch_for_cowpy = CAT_CAT.out.file_out.map{ meta, file -> file }

        // ASCII-Art der Begrüßungen mit cowpy generieren
        COWPY(ch_for_cowpy, params.character)
    ```

Wir rufen jetzt `COWPY` direkt auf `CAT_CAT.out.file_out` auf.

Dadurch müssen wir den `ch_for_cowpy`-Kanal nicht mehr konstruieren, sodass diese Zeile (und ihre Kommentarzeile) vollständig gelöscht werden kann.

#### 1.2.3. Den emit-Block im Workflow aktualisieren

Da `COWPY` jetzt eine benannte Ausgabe `cowpy_output` ausgibt, können wir den `emit:`-Block der `hello.nf`-Workflow-Datei aktualisieren, um diese zu verwenden.

=== "Danach"

    ```groovy title="core-hello/workflows/hello.nf" linenums="60" hl_lines="2"
        emit:
        cowpy_hellos   = COWPY.out.cowpy_output
        versions       = ch_versions
    ```

=== "Vorher"

    ```groovy title="core-hello/workflows/hello.nf" linenums="60" hl_lines="2"
        emit:
        cowpy_hellos   = COWPY.out
        versions       = ch_versions
    ```

Dies ist technisch nicht erforderlich, aber es ist gute Praxis, benannte Ausgaben zu verwenden, wann immer möglich.

#### 1.2.4. Die Pipeline ausführen, um sie zu testen

Lass uns den Workflow ausführen, um zu testen, dass nach diesen Änderungen alles korrekt funktioniert.

```bash
nextflow run . --outdir core-hello-results -profile test,docker --validate_params false
```

??? success "Befehlsausgabe"

    ```console
     N E X T F L O W   ~  version 25.04.3

    Launching `./main.nf` [modest_saha] DSL2 - revision: b9e9b3b8de

    Downloading plugin nf-schema@2.5.1
    Input/output options
      input                     : /workspaces/training/hello-nf-core/core-hello/assets/greetings.csv
      outdir                    : core-hello-results

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      validate_params           : false
      trace_report_suffix       : 2025-12-27_06-16-55

    Core Nextflow options
      runName                   : modest_saha
      containerEngine           : docker
      launchDir                 : /workspaces/training/hello-nf-core/core-hello
      workDir                   : /workspaces/training/hello-nf-core/core-hello/work
      projectDir                : /workspaces/training/hello-nf-core/core-hello
      userName                  : root
      profile                   : test,docker
      configFiles               : /workspaces/training/hello-nf-core/core-hello/nextflow.config

    !! Only displaying parameters that differ from the pipeline defaults !!
    ------------------------------------------------------
    executor >  local (8)
    [a8/447993] CORE_HELLO:HELLO:sayHello (3)       [100%] 3 of 3 ✔
    [00/1fc59c] CORE_HELLO:HELLO:convertToUpper (3) [100%] 3 of 3 ✔
    [57/ac800d] CORE_HELLO:HELLO:CAT_CAT (test)     [100%] 1 of 1 ✔
    [b7/092f2b] CORE_HELLO:HELLO:COWPY              [100%] 1 of 1 ✔
    -[core/hello] Pipeline completed successfully-
    ```

Die Pipeline sollte erfolgreich laufen, wobei Metadaten jetzt von `CAT_CAT` durch `COWPY` fließen.

Damit ist abgeschlossen, was wir tun mussten, um `COWPY` Metadaten-Tupel verarbeiten zu lassen.
Schauen wir uns jetzt an, was wir sonst noch tun können, um die nf-core-Modul-Muster zu nutzen.

### 1.3. Tool-Argument-Konfiguration mit `ext.args` zentralisieren

In seinem aktuellen Zustand erwartet der `COWPY`-Prozess einen Wert für den `character`-Parameter.
Daher müssen wir jedes Mal, wenn wir den Prozess aufrufen, einen Wert angeben, selbst wenn wir mit den vom Tool gesetzten Standardwerten zufrieden wären.
Für `COWPY` ist dies zugegebenermaßen kein großes Problem, aber für Tools mit vielen optionalen Parametern kann es ziemlich umständlich werden.

Das nf-core-Projekt empfiehlt die Verwendung einer Nextflow-Funktion namens [`ext.args`](https://www.nextflow.io/docs/latest/reference/process.html#ext), um Tool-Argumente bequemer über Konfigurationsdateien zu verwalten.

Anstatt Prozess-Inputs für jede Tool-Option zu deklarieren, schreibst du das Modul so, dass es `ext.args` beim Aufbau seiner Befehlszeile referenziert.
Dann ist es nur noch eine Frage der Einrichtung der `ext.args`-Variable, um die Argumente und Werte zu halten, die du in der `modules.config`-Datei verwenden möchtest, die Konfigurationsdetails für alle Module konsolidiert.
Nextflow wird diese Argumente mit ihren Werten zur Laufzeit in die Tool-Befehlszeile einfügen.

Lass uns diesen Ansatz auf das `COWPY`-Modul anwenden.
Wir müssen folgende Änderungen vornehmen:

1. Das `COWPY`-Modul aktualisieren
2. `ext.args` in der `modules.config`-Datei konfigurieren
3. Den `hello.nf`-Workflow aktualisieren

Sobald wir das alles erledigt haben, führen wir die Pipeline aus, um zu testen, dass alles noch wie zuvor funktioniert.

#### 1.3.1. Das `COWPY`-Modul aktualisieren

Los geht's.
Öffne die `cowpy.nf`-Moduldatei (unter `core-hello/modules/local/`) und modifiziere sie, um `ext.args` zu referenzieren, wie unten gezeigt.

=== "Danach"

    ```groovy title="modules/local/cowpy.nf" linenums="1" hl_lines="18 20"
    #!/usr/bin/env nextflow

    // ASCII-Art mit cowpy generieren (https://github.com/jeffbuttars/cowpy)
    process COWPY {

        publishDir 'results', mode: 'copy'

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
        conda 'conda-forge::cowpy==1.1.5'

        input:
            tuple val(meta), path(input_file)

        output:
            tuple val(meta), path("cowpy-${input_file}"), emit: cowpy_output

        script:
        def args = task.ext.args ?: ''
        """
        cat $input_file | cowpy $args > cowpy-${input_file}
        """
    }
    ```

=== "Vorher"

    ```groovy title="core-hello/modules/local/cowpy.nf" linenums="1" hl_lines="13 20"
    #!/usr/bin/env nextflow

    // ASCII-Art mit cowpy generieren (https://github.com/jeffbuttars/cowpy)
    process COWPY {

        publishDir 'results', mode: 'copy'

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
        conda 'conda-forge::cowpy==1.1.5'

        input:
            tuple val(meta), path(input_file)
            val character

        output:
            tuple val(meta), path("cowpy-${input_file}"), emit: cowpy_output

        script:
        """
        cat $input_file | cowpy -c "$character" > cowpy-${input_file}
        """
    }
    ```

Du kannst sehen, dass wir drei Änderungen vorgenommen haben.

1. **Im `input:`-Block haben wir die `val character`-Eingabe entfernt.**
   In Zukunft werden wir dieses Argument über die `ext.args`-Konfiguration bereitstellen, wie weiter unten beschrieben.

2. **Im `script:`-Block haben wir die Zeile `def args = task.ext.args ?: ''` hinzugefügt.**
   Diese Zeile verwendet den `?:`-Operator, um den Wert der `args`-Variable zu bestimmen: den Inhalt von `task.ext.args`, wenn er nicht leer ist, oder einen leeren String, wenn er leer ist.
   Beachte, dass wir zwar allgemein auf `ext.args` verweisen, dieser Code aber `task.ext.args` referenzieren muss, um die Modul-Level-`ext.args`-Konfiguration herauszuziehen.

3. **In der Befehlszeile haben wir `-c "$character"` durch `$args` ersetzt.**
   Hier wird Nextflow alle in `ext.args` in der `modules.config`-Datei gesetzten Tool-Argumente einfügen.

Dadurch ist die Modul-Schnittstelle jetzt einfacher: Sie erwartet nur die wesentlichen Metadaten- und Datei-Inputs.

!!! note

    Der `?:`-Operator wird oft als 'Elvis-Operator' bezeichnet, weil er wie ein seitliches Elvis-Presley-Gesicht aussieht, wobei das `?`-Zeichen die Welle in seinem Haar symbolisiert.

#### 1.3.2. `ext.args` in der `modules.config`-Datei konfigurieren

Jetzt, da wir die `character`-Deklaration aus dem Modul entfernt haben, müssen wir sie zu `ext.args` in der `modules.config`-Konfigurationsdatei hinzufügen.

Konkret werden wir diesen kleinen Code-Block zum `process {}`-Block hinzufügen:

```groovy title="Hinzuzufügender Code"
withName: 'COWPY' {
    ext.args = { "-c ${params.character}" }
}
```

Die `withName:`-Syntax weist diese Konfiguration nur dem `COWPY`-Prozess zu, und `ext.args = { "-c ${params.character}" }` komponiert einfach einen String, der den Wert des `character`-Parameters enthält.
Beachte die Verwendung von geschweiften Klammern, die Nextflow anweisen, den Wert des Parameters zur Laufzeit auszuwerten.

Macht das Sinn? Fügen wir es hinzu.

Öffne `conf/modules.config` und füge den Konfigurationscode innerhalb des `process {}`-Blocks hinzu, wie unten gezeigt.

=== "Danach"

    ```groovy title="core-hello/conf/modules.config" linenums="13" hl_lines="8-10"
    process {
        publishDir = [
            path: { "${params.outdir}/${task.process.tokenize(':')[-1].tokenize('_')[0].toLowerCase()}" },
            mode: params.publish_dir_mode,
            saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
        ]

        withName: 'COWPY' {
            ext.args = { "-c ${params.character}" }
        }
    }
    ```

=== "Vorher"

    ```groovy title="core-hello/conf/modules.config" linenums="13"
    process {
        publishDir = [
            path: { "${params.outdir}/${task.process.tokenize(':')[-1].tokenize('_')[0].toLowerCase()}" },
            mode: params.publish_dir_mode,
            saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
        ]
    }
    ```

Hoffentlich kannst du dir vorstellen, dass alle Module in einer Pipeline ihre `ext.args` in dieser Datei spezifiziert haben, mit folgenden Vorteilen:

- Die **Modul-Schnittstelle bleibt einfach** - Sie akzeptiert nur die wesentlichen Metadaten- und Datei-Inputs
- Die **Pipeline stellt immer noch `params.character` bereit** - Endnutzer\*innen können es weiterhin wie zuvor konfigurieren
- Das **Modul ist jetzt portabel** - Es kann in anderen Pipelines wiederverwendet werden, ohne einen spezifischen Parameternamen zu erwarten
- Die Konfiguration ist **zentralisiert** in `modules.config`, was die Workflow-Logik sauber hält

Indem wir die `modules.config`-Datei als Ort verwenden, an dem alle Pipelines die Konfiguration pro Modul zentralisieren, machen wir unsere Module über verschiedene Pipelines hinweg wiederverwendbarer.

#### 1.3.3. Den `hello.nf`-Workflow aktualisieren

Da das `COWPY`-Modul den `character`-Parameter nicht mehr als Eingabe benötigt, müssen wir den Workflow-Aufruf entsprechend aktualisieren.

Öffne die `hello.nf`-Workflow-Datei (unter `core-hello/workflows/`) und aktualisiere den Aufruf von `COWPY` wie unten gezeigt.

=== "Danach"

    ```groovy title="core-hello/workflows/hello.nf" linenums="39" hl_lines="2"
        // ASCII-Art der Begrüßungen mit cowpy generieren
        COWPY(CAT_CAT.out.file_out)
    ```

=== "Vorher"

    ```groovy title="core-hello/workflows/hello.nf" linenums="39" hl_lines="2"
        // ASCII-Art der Begrüßungen mit cowpy generieren
        COWPY(CAT_CAT.out.file_out, params.character)
    ```

Der Workflow-Code ist jetzt sauberer: Wir müssen `params.character` nicht direkt an den Prozess übergeben.
Die Modul-Schnittstelle wird minimal gehalten, was sie portabler macht, während die Pipeline die explizite Option weiterhin über die Konfiguration bereitstellt.

#### 1.3.4. Die Pipeline ausführen, um sie zu testen

Lass uns testen, dass der Workflow noch wie erwartet funktioniert, indem wir einen anderen Charakter angeben, um zu überprüfen, dass die `ext.args`-Konfiguration funktioniert.

Führe diesen Befehl mit `kosh` aus, einer der... enigmatischeren Optionen:

```bash
nextflow run . --outdir core-hello-results -profile test,docker --validate_params false --character kosh
```

??? success "Befehlsausgabe"

    ```console
     N E X T F L O W   ~  version 25.04.3

    Launching `./main.nf` [exotic_planck] DSL2 - revision: b9e9b3b8de

    Input/output options
      input                     : /workspaces/training/hello-nf-core/core-hello/assets/greetings.csv
      outdir                    : core-hello-results

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      validate_params           : false
      trace_report_suffix       : 2025-12-27_06-23-13

    Core Nextflow options
      runName                   : exotic_planck
      containerEngine           : docker
      launchDir                 : /workspaces/training/hello-nf-core/core-hello
      workDir                   : /workspaces/training/hello-nf-core/core-hello/work
      projectDir                : /workspaces/training/hello-nf-core/core-hello
      userName                  : root
      profile                   : test,docker
      configFiles               : /workspaces/training/hello-nf-core/core-hello/nextflow.config

    !! Only displaying parameters that differ from the pipeline defaults !!
    ------------------------------------------------------
    executor >  local (8)
    [13/9e3c0e] CORE_HELLO:HELLO:sayHello (3)       [100%] 3 of 3 ✔
    [e2/5b0ee5] CORE_HELLO:HELLO:convertToUpper (3) [100%] 3 of 3 ✔
    [b6/4fb569] CORE_HELLO:HELLO:CAT_CAT (test)     [100%] 1 of 1 ✔
    [38/eb29ea] CORE_HELLO:HELLO:COWPY              [100%] 1 of 1 ✔
    -[core/hello] Pipeline completed successfully-
    ```

Dies sollte wie zuvor erfolgreich laufen.

Lass uns überprüfen, dass die `ext.args`-Konfiguration funktioniert hat, indem wir die Ausgabe überprüfen.
Finde die Ausgabe im Dateibrowser oder verwende den Task-Hash (der `38/eb29ea`-Teil im obigen Beispiel), um die Ausgabedatei anzusehen:

```bash
cat work/38/eb29ea*/cowpy-test.txt
```

??? success "Befehlsausgabe"

    ```console
    _________
    / HELLO   \
    | HOLà    |
    \ BONJOUR /
    ---------
        \
        \
          \
      ___       _____     ___
    /   \     /    /|   /   \
    |     |   /    / |  |     |
    |     |  /____/  |  |     |
    |     |  |    |  |  |     |
    |     |  | {} | /   |     |
    |     |  |____|/    |     |
    |     |    |==|     |     |
    |      \___________/      |
    |                         |
    |                         |
    ```

Du solltest die ASCII-Art mit dem `kosh`-Charakter sehen, was bestätigt, dass die `ext.args`-Konfiguration funktioniert hat!

??? info "(Optional) Die Befehlsdatei inspizieren"

    Wenn du genau sehen möchtest, wie die Konfiguration angewendet wurde, kannst du die `.command.sh`-Datei inspizieren:

    ```bash
    cat work/38/eb29ea*/.command.sh
    ```

    Du wirst den `cowpy`-Befehl mit dem `-c kosh`-Argument sehen:

    ```console
    #!/usr/bin/env bash
    ...
    cat test.txt | cowpy -c kosh > cowpy-test.txt
    ```

    Dies zeigt, dass die `.command.sh`-Datei korrekt basierend auf der `ext.args`-Konfiguration generiert wurde.

Nimm dir einen Moment Zeit, um darüber nachzudenken, was wir hier erreicht haben.
Dieser Ansatz hält die Modul-Schnittstelle auf wesentliche Daten fokussiert (Dateien, Metadaten und alle obligatorischen pro-Sample-Parameter), während Optionen, die das Verhalten des Tools steuern, separat über die Konfiguration behandelt werden.

Dies mag für ein einfaches Tool wie `cowpy` unnötig erscheinen, kann aber einen großen Unterschied für Datenanalyse-Tools machen, die viele optionale Argumente haben.

Um die Vorteile dieses Ansatzes zusammenzufassen:

- **Saubere Schnittstelle**: Das Modul konzentriert sich auf wesentliche Dateneingaben (Metadaten und Dateien)
- **Flexibilität**: Nutzer\*innen können Tool-Argumente über die Konfiguration angeben, einschließlich sample-spezifischer Werte
- **Konsistenz**: Alle nf-core-Module folgen diesem Muster
- **Portabilität**: Module können ohne fest codierte Tool-Optionen wiederverwendet werden
- **Keine Workflow-Änderungen**: Das Hinzufügen oder Ändern von Tool-Optionen erfordert keine Aktualisierung des Workflow-Codes

!!! note

    Das `ext.args`-System hat leistungsstarke zusätzliche Fähigkeiten, die hier nicht behandelt werden, einschließlich des dynamischen Wechselns von Argumentwerten basierend auf Metadaten. Siehe die [nf-core-Modul-Spezifikationen](https://nf-co.re/docs/guidelines/components/modules) für weitere Details.

### 1.4. Ausgabe-Benennung mit `ext.prefix` standardisieren

Jetzt, da wir dem `COWPY`-Prozess Zugriff auf die Metamap gegeben haben, können wir beginnen, ein weiteres nützliches nf-core-Muster zu nutzen: die Benennung von Ausgabedateien basierend auf Metadaten.

Hier werden wir eine Nextflow-Funktion namens `ext.prefix` verwenden, die es uns ermöglicht, die Ausgabedatei-Benennung über Module hinweg mit `meta.id` (dem in der Metamap enthaltenen Identifier) zu standardisieren, während wir Module bei Bedarf individuell konfigurieren können.

Dies wird ähnlich sein wie das, was wir mit `ext.args` gemacht haben, mit einigen Unterschieden, die wir im Verlauf detailliert erläutern werden.

Lass uns diesen Ansatz auf das `COWPY`-Modul anwenden.
Wir müssen folgende Änderungen vornehmen:

1. Das `COWPY`-Modul aktualisieren
2. `ext.prefix` in der `modules.config`-Datei konfigurieren

(Keine Änderungen am Workflow erforderlich.)

Sobald wir das erledigt haben, führen wir die Pipeline aus, um zu testen, dass alles noch wie zuvor funktioniert.

#### 1.4.1. Das `COWPY`-Modul aktualisieren

Öffne die `cowpy.nf`-Moduldatei (unter `core-hello/modules/local/`) und modifiziere sie, um `ext.prefix` zu referenzieren, wie unten gezeigt.

=== "Danach"

    ```groovy title="core-hello/modules/local/cowpy.nf" linenums="1" hl_lines="2 6 8"
        output:
            tuple val(meta), path("${prefix}.txt"), emit: cowpy_output

        script:
        def args = task.ext.args ?: ''
        prefix = task.ext.prefix ?: "${meta.id}"
        """
        cat $input_file | cowpy $args > ${prefix}.txt
        """
    }
    ```

=== "Vorher"

    ```groovy title="core-hello/modules/local/cowpy.nf" linenums="1" hl_lines="2 7"
        output:
            tuple val(meta), path("cowpy-${input_file}"), emit: cowpy_output

        script:
        def args = task.ext.args ?: ''
        """
        cat $input_file | cowpy $args > cowpy-${input_file}
        """
    }
    ```

Du kannst sehen, dass wir drei Änderungen vorgenommen haben.

1. **Im `script:`-Block haben wir die Zeile `prefix = task.ext.prefix ?: "${meta.id}"` hinzugefügt.**
   Diese Zeile verwendet den `?:`-Operator, um den Wert der `prefix`-Variable zu bestimmen: den Inhalt von `task.ext.prefix`, wenn er nicht leer ist, oder den Identifier aus der Metamap (`meta.id`), wenn er leer ist.
   Beachte, dass wir zwar allgemein auf `ext.prefix` verweisen, dieser Code aber `task.ext.prefix` referenzieren muss, um die Modul-Level-`ext.prefix`-Konfiguration herauszuziehen.

2. **In der Befehlszeile haben wir `cowpy-${input_file}` durch `${prefix}.txt` ersetzt.**
   Hier wird Nextflow den Wert von `prefix` einfügen, der durch die obige Zeile bestimmt wurde.

3. **Im `output:`-Block haben wir `path("cowpy-${input_file}")` durch `path("${prefix}.txt")` ersetzt.**
   Dies wiederholt einfach, was der Dateipfad gemäß dem sein wird, was in der Befehlszeile geschrieben ist.

Dadurch wird der Ausgabedateiname jetzt mit einem sinnvollen Standard (dem Identifier aus der Metamap) kombiniert mit der entsprechenden Dateiformaterweiterung konstruiert.

#### 1.4.2. `ext.prefix` in der `modules.config`-Datei konfigurieren

In diesem Fall ist der sinnvolle Standard für unseren Geschmack nicht ausreichend ausdrucksstark; wir möchten ein benutzerdefiniertes Benennungsmuster verwenden, das den Tool-Namen enthält, `cowpy-<id>.txt`, wie wir es zuvor hatten.

Wir werden das tun, indem wir `ext.prefix` in `modules.config` konfigurieren, genau wie wir es für den `character`-Parameter mit `ext.args` getan haben, außer dass diesmal der `withName: 'COWPY' {}`-Block bereits existiert und wir nur die folgende Zeile hinzufügen müssen:

```groovy title="Hinzuzufügender Code"
ext.prefix = { "cowpy-${meta.id}" }
```

Dies wird den gewünschten String komponieren.
Beachte, dass wir wieder geschweifte Klammern verwenden, diesmal um Nextflow anzuweisen, den Wert von `meta.id` zur Laufzeit auszuwerten.

Fügen wir es hinzu.

Öffne `conf/modules.config` und füge den Konfigurationscode innerhalb des `process {}`-Blocks hinzu, wie unten gezeigt.

=== "Danach"

    ```groovy title="core-hello/conf/modules.config" linenums="21" hl_lines="3"
        withName: 'COWPY' {
            ext.args = { "-c ${params.character}" }
            ext.prefix = { "cowpy-${meta.id}" }
        }
    ```

=== "Vorher"

    ```groovy title="core-hello/conf/modules.config" linenums="21"
        withName: 'COWPY' {
            ext.args = { "-c ${params.character}" }
        }
    ```

Falls du dich fragst: Die `ext.prefix`-Closure hat Zugriff auf das richtige Stück Metadaten, weil die Konfiguration im Kontext der Prozessausführung ausgewertet wird, wo Metadaten verfügbar sind.

#### 1.4.3. Die Pipeline ausführen, um sie zu testen

Lass uns testen, dass der Workflow noch wie erwartet funktioniert.

```bash
nextflow run . --outdir core-hello-results -profile test,docker --validate_params false
```

??? success "Befehlsausgabe"

    ```console
     N E X T F L O W   ~  version 25.04.3

    Launching `./main.nf` [admiring_turing] DSL2 - revision: b9e9b3b8de

    Input/output options
      input                     : /workspaces/training/hello-nf-core/core-hello/assets/greetings.csv
      outdir                    : core-hello-results

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      validate_params           : false
      trace_report_suffix       : 2025-12-27_06-29-02

    Core Nextflow options
      runName                   : admiring_turing
      containerEngine           : docker
      launchDir                 : /workspaces/training/hello-nf-core/core-hello
      workDir                   : /workspaces/training/hello-nf-core/core-hello/work
      projectDir                : /workspaces/training/hello-nf-core/core-hello
      userName                  : root
      profile                   : test,docker
      configFiles               : /workspaces/training/hello-nf-core/core-hello/nextflow.config

    !! Only displaying parameters that differ from the pipeline defaults !!
    ------------------------------------------------------
    executor >  local (8)
    [b2/e08524] CORE_HELLO:HELLO:sayHello (3)       [100%] 3 of 3 ✔
    [13/88939f] CORE_HELLO:HELLO:convertToUpper (3) [100%] 3 of 3 ✔
    [23/4554e1] CORE_HELLO:HELLO:CAT_CAT (test)     [100%] 1 of 1 ✔
    [a3/c6cbe9] CORE_HELLO:HELLO:COWPY              [100%] 1 of 1 ✔
    -[core/hello] Pipeline completed successfully-
    ```

Schau dir die Ausgabe im Ergebnisverzeichnis an.
Du solltest die cowpy-Ausgabedatei mit der gleichen Benennung wie zuvor sehen: `cowpy-test.txt`, basierend auf dem Standard-Batch-Namen.

??? abstract "Verzeichnisinhalt"

    ```console hl_lines="3"
    results
    ├── Bonjour-output.txt
    ├── cowpy-test.txt
    ├── Hello-output.txt
    ├── Holà-output.txt
    ├── UPPER-Bonjour-output.txt
    ├── UPPER-Hello-output.txt
    └── UPPER-Holà-output.txt
    ```

Du kannst gerne die `ext.prefix`-Konfiguration in `conf/modules.config` ändern, um dich zu überzeugen, dass du das Benennungsmuster ändern kannst, ohne Änderungen am Modul- oder Workflow-Code vornehmen zu müssen.

Alternativ kannst du dies auch erneut mit einem anderen `--batch`-Parameter ausführen, der auf der Befehlszeile angegeben wird, um dich zu überzeugen, dass dieser Teil weiterhin spontan anpassbar ist.

Dies demonstriert, wie `ext.prefix` es dir ermöglicht, deine bevorzugte Benennungskonvention beizubehalten und gleichzeitig die Modul-Schnittstelle flexibel zu halten.

Um die Vorteile dieses Ansatzes zusammenzufassen:

- **Standardisierte Benennung**: Ausgabedateien werden typischerweise mit Sample-IDs aus Metadaten benannt
- **Konfigurierbar**: Nutzer\*innen können die Standardbenennung bei Bedarf überschreiben
- **Konsistent**: Alle nf-core-Module folgen diesem Muster
- **Vorhersagbar**: Es ist einfach zu wissen, wie Ausgabedateien heißen werden

Ziemlich gut, oder?
Nun, es gibt noch eine wichtige Änderung, die wir vornehmen müssen, um unser Modul an die nf-core-Richtlinien anzupassen.

### 1.5. Die Publishing-Konfiguration zentralisieren

Du hast vielleicht bemerkt, dass wir Ausgaben in zwei verschiedene Verzeichnisse veröffentlicht haben:

- **`results`** — Das ursprüngliche Ausgabeverzeichnis, das wir von Anfang an für unsere lokalen Module verwendet haben, individuell mit pro-Modul-`publishDir`-Direktiven gesetzt;
- **`core-hello-results`** — Das Ausgabeverzeichnis, das mit `--outdir` auf der Befehlszeile gesetzt wurde, das die nf-core-Logs und die von `CAT_CAT` veröffentlichten Ergebnisse erhalten hat.

Das ist unordentlich und suboptimal; es wäre besser, einen Ort für alles zu haben.
Natürlich könnten wir in jedes unserer lokalen Module gehen und die `publishDir`-Direktive manuell aktualisieren, um das `core-hello-results`-Verzeichnis zu verwenden, aber was ist beim nächsten Mal, wenn wir das Ausgabeverzeichnis ändern möchten?

Einzelne Module Publishing-Entscheidungen treffen zu lassen, ist eindeutig nicht der richtige Weg, besonders in einer Welt, in der dasselbe Modul in vielen verschiedenen Pipelines verwendet werden könnte, von Leuten, die unterschiedliche Bedürfnisse oder Präferenzen haben.
Wir möchten in der Lage sein zu kontrollieren, wo Ausgaben auf der Ebene der Workflow-Konfiguration veröffentlicht werden.

"Hey", könntest du sagen, "`CAT_CAT` sendet seine Ausgaben an `--outdir`. Vielleicht sollten wir seine `publishDir`-Direktive kopieren?"

Ja, das ist eine großartige Idee.

Außer dass es keine `publishDir`-Direktive hat. (Schau dir den Modul-Code an.)

Das liegt daran, dass nf-core-Pipelines die Kontrolle auf Workflow-Ebene zentralisieren, indem sie `publishDir` in `conf/modules.config` konfigurieren, anstatt in einzelnen Modulen.
Konkret deklariert das nf-core-Template eine Standard-`publishDir`-Direktive (mit einer vordefinierten Verzeichnisstruktur), die für alle Module gilt, es sei denn, eine überschreibende Direktive wird bereitgestellt.

Klingt das nicht großartig? Könnte es sein, dass wir, um diese Standard-Direktive zu nutzen, nur die aktuelle `publishDir`-Direktive aus unseren lokalen Modulen entfernen müssen?

Lass uns das an `COWPY` ausprobieren, um zu sehen, was passiert, dann schauen wir uns den Code für die Standardkonfiguration an, um zu verstehen, wie es funktioniert.

Schließlich werden wir demonstrieren, wie man das Standardverhalten bei Bedarf überschreibt.

#### 1.5.1. Die `publishDir`-Direktive aus `COWPY` entfernen

Lass uns das tun.
Öffne die `cowpy.nf`-Moduldatei (unter `core-hello/modules/local/`) und entferne die `publishDir`-Direktive wie unten gezeigt.

=== "Danach"

    ```groovy title="core-hello/modules/local/cowpy.nf (Auszug)" linenums="1"
    #!/usr/bin/env nextflow

    // ASCII-Art mit cowpy generieren (https://github.com/jeffbuttars/cowpy)
    process COWPY {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
        conda 'conda-forge::cowpy==1.1.5'
    ```

=== "Vorher"

    ```groovy title="core-hello/modules/local/cowpy.nf (Auszug)" linenums="1" hl_lines="6"
    #!/usr/bin/env nextflow

    // ASCII-Art mit cowpy generieren (https://github.com/jeffbuttars/cowpy)
    process COWPY {

        publishDir 'results', mode: 'copy'

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
        conda 'conda-forge::cowpy==1.1.5'

    ```

Das war's!

#### 1.5.2. Die Pipeline ausführen, um sie zu testen

Schauen wir uns an, was passiert, wenn wir die Pipeline jetzt ausführen.

```bash
nextflow run . --outdir core-hello-results -profile test,docker --validate_params false
```

??? success "Befehlsausgabe"

    ```console
     N E X T F L O W   ~  version 25.04.3

    Launching `./main.nf` [silly_caravaggio] DSL2 - revision: b9e9b3b8de

    Input/output options
      input                     : /workspaces/training/hello-nf-core/core-hello/assets/greetings.csv
      outdir                    : core-hello-results

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      validate_params           : false
      trace_report_suffix       : 2025-12-27_06-35-56

    Core Nextflow options
      runName                   : silly_caravaggio
      containerEngine           : docker
      launchDir                 : /workspaces/training/hello-nf-core/core-hello
      workDir                   : /workspaces/training/hello-nf-core/core-hello/work
      projectDir                : /workspaces/training/hello-nf-core/core-hello
      userName                  : root
      profile                   : test,docker
      configFiles               : /workspaces/training/hello-nf-core/core-hello/nextflow.config

    !! Only displaying parameters that differ from the pipeline defaults !!
    ------------------------------------------------------
    executor >  local (8)
    [db/39978e] CORE_HELLO:HELLO:sayHello (3)       [100%] 3 of 3 ✔
    [b5/bf6a8d] CORE_HELLO:HELLO:convertToUpper (3) [100%] 3 of 3 ✔
    [b7/c61842] CORE_HELLO:HELLO:CAT_CAT (test)     [100%] 1 of 1 ✔
    [46/5839d6] CORE_HELLO:HELLO:COWPY              [100%] 1 of 1 ✔
    -[core/hello] Pipeline completed successfully-
    ```

Schau dir dein aktuelles Arbeitsverzeichnis an.
Jetzt enthält `core-hello-results` auch die Ausgaben des `COWPY`-Moduls.

??? abstract "Verzeichnisinhalt"

    ```console hl_lines="4-5"
    core-hello-results/
    ├── cat
    │   └── test.txt
    ├── cowpy
    │   └── cowpy-test.txt
    └── pipeline_info
        ├── execution_report_2025-12-27_06-16-55.html
        ├── execution_report_2025-12-27_06-23-13.html
        ├── execution_report_2025-12-27_06-29-02.html
        ├── execution_report_2025-12-27_06-35-56.html
        ├── execution_timeline_2025-12-27_06-16-55.html
        ├── execution_timeline_2025-12-27_06-23-13.html
        ├── execution_timeline_2025-12-27_06-29-02.html
        ├── execution_timeline_2025-12-27_06-35-56.html
        ├── execution_trace_2025-12-27_06-16-55.txt
        ├── execution_trace_2025-12-27_06-23-13.txt
        ├── execution_trace_2025-12-27_06-29-02.txt
        ├── execution_trace_2025-12-27_06-35-56.txt
        ├── hello_software_versions.yml
        ├── params_2025-12-27_06-17-00.json
        ├── params_2025-12-27_06-23-17.json
        ├── params_2025-12-27_06-29-07.json
        ├── params_2025-12-27_06-36-01.json
        ├── pipeline_dag_2025-12-27_06-16-55.html
        ├── pipeline_dag_2025-12-27_06-23-13.html
        ├── pipeline_dag_2025-12-27_06-29-02.html
        └── pipeline_dag_2025-12-27_06-35-56.html
    ```

Du kannst sehen, dass Nextflow diese Hierarchie von Verzeichnissen basierend auf den Namen des Workflows und des Moduls erstellt hat.

Der verantwortliche Code befindet sich in der `conf/modules.config`-Datei.
Dies ist die Standard-`publishDir`-Konfiguration, die Teil des nf-core-Templates ist und für alle Prozesse gilt:

```groovy
process {
    publishDir = [
        path: { "${params.outdir}/${task.process.tokenize(':')[-1].tokenize('_')[0].toLowerCase()}" },
        mode: params.publish_dir_mode,
        saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    ]
}
```

Das mag kompliziert aussehen, also schauen wir uns jede der drei Komponenten an:

- **`path:`** Bestimmt das Ausgabeverzeichnis basierend auf dem Prozessnamen.
  Der vollständige Name eines Prozesses, der in `task.process` enthalten ist, umfasst die Hierarchie von Workflow- und Modul-Importen (wie `CORE_HELLO:HELLO:CAT_CAT`).
  Die `tokenize`-Operationen entfernen diese Hierarchie, um nur den Prozessnamen zu erhalten, nehmen dann den ersten Teil vor einem Unterstrich (falls zutreffend) und konvertieren ihn in Kleinbuchstaben.
  Dies bestimmt, dass die Ergebnisse von `CAT_CAT` nach `${params.outdir}/cat/` veröffentlicht werden.
- **`mode:`** Steuert, wie Dateien veröffentlicht werden (copy, symlink, etc.).
  Dies ist über den `params.publish_dir_mode`-Parameter konfigurierbar.
- **`saveAs:`** Filtert, welche Dateien veröffentlicht werden.
  Dieses Beispiel schließt `versions.yml`-Dateien aus, indem es `null` für sie zurückgibt, was verhindert, dass sie veröffentlicht werden.

Dies bietet eine konsistente Logik für die Organisation von Ausgaben.

Die Ausgabe sieht noch besser aus, wenn alle Module in einer Pipeline diese Konvention übernehmen, also kannst du gerne die `publishDir`-Direktiven aus den anderen Modulen in deiner Pipeline löschen.
Dieser Standard wird auch auf Module angewendet, die wir nicht explizit modifiziert haben, um den nf-core-Richtlinien zu folgen.

Allerdings könntest du entscheiden, dass du deine Eingaben anders organisieren möchtest, und die gute Nachricht ist, dass es einfach ist, dies zu tun.

#### 1.5.3. Den Standard überschreiben

Um die Standard-`publishDir`-Direktive zu überschreiben, kannst du einfach deine eigenen Direktiven zur `conf/modules.config`-Datei hinzufügen.

Zum Beispiel könntest du den Standard für einen einzelnen Prozess mit dem `withName:`-Selektor überschreiben, wie in diesem Beispiel, wo wir eine benutzerdefinierte `publishDir`-Direktive für den 'COWPY'-Prozess hinzufügen.

```groovy title="core-hello/conf/modules.config" linenums="13" hl_lines="8-10"
process {
    publishDir = [
        path: { "${params.outdir}/${task.process.tokenize(':')[-1].tokenize('_')[0].toLowerCase()}" },
        mode: params.publish_dir_mode,
        saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    ]

    withName: 'COWPY' {
        ext.args = { "-c ${params.character}" }
        publishDir = [
            path: 'my_custom_results'
        ]
    }
}
```

Wir werden diese Änderung nicht tatsächlich vornehmen, aber du kannst gerne damit spielen und sehen, welche Logik du implementieren kannst.

Der Punkt ist, dass dieses System dir das Beste aus beiden Welten gibt: Konsistenz standardmäßig und die Flexibilität, die Konfiguration bei Bedarf anzupassen.

Zusammenfassend erhältst du:

- **Single Source of Truth**: Die gesamte Publishing-Konfiguration befindet sich in `modules.config`
- **Nützlicher Standard**: Prozesse funktionieren out-of-the-box ohne pro-Modul-Konfiguration
- **Einfache Anpassung**: Publishing-Verhalten in der Konfiguration überschreiben, nicht im Modul-Code
- **Portable Module**: Module codieren keine Ausgabeorte fest

Dies vervollständigt die Reihe von nf-core-Modul-Features, die du unbedingt lernen solltest, aber es gibt andere, über die du in den [nf-core-Modul-Spezifikationen](https://nf-co.re/docs/guidelines/components/modules) lesen kannst.

### Fazit

Du weißt jetzt, wie du lokale Module anpasst, um nf-core-Konventionen zu folgen:

- Entwirf deine Module so, dass sie Metadaten-Tupel akzeptieren und propagieren;
- Verwende `ext.args`, um Modul-Schnittstellen minimal und portabel zu halten;
- Verwende `ext.prefix` für konfigurierbare, standardisierte Ausgabedatei-Benennung;
- Übernimm die standardmäßige zentralisierte `publishDir`-Direktive für eine konsistente Ergebnisverzeichnisstruktur.

### Wie geht es weiter?

Lerne, wie du die eingebauten template-basierten Tools von nf-core verwendest, um Module auf einfache Weise zu erstellen.

---

## 2. Ein Modul mit den nf-core-Tools erstellen

Jetzt, da du die nf-core-Modul-Muster durch manuelle Anwendung gelernt hast, schauen wir uns an, wie du Module in der Praxis erstellen würdest.

### 2.1. Ein Modul-Gerüst aus einem Template generieren

Ähnlich wie bei der Erstellung von Pipelines stellt das nf-core-Projekt Tools bereit, um korrekt strukturierte Module basierend auf einem Template zu generieren, mit all diesen Mustern von Anfang an eingebaut.

#### 2.1.1. Den Modul-Erstellungsbefehl ausführen

Der `nf-core modules create`-Befehl generiert ein Modul-Template, das bereits allen Konventionen folgt, die du gelernt hast.

Lass uns eine neue Version des `COWPY`-Moduls mit einem minimalen Template erstellen, indem wir diesen Befehl ausführen:

```bash
nf-core modules create --empty-template COWPY
```

Das `--empty-template`-Flag erstellt ein sauberes Starter-Template ohne zusätzlichen Code, was es einfacher macht, die wesentliche Struktur zu sehen.

Der Befehl läuft interaktiv und führt dich durch die Einrichtung.
Er sucht automatisch Tool-Informationen aus Paket-Repositories wie Bioconda und bio.tools, um Metadaten vorab auszufüllen.

Du wirst zu mehreren Konfigurationsoptionen aufgefordert:

- **Autor\*innen-Informationen**: Dein GitHub-Benutzername für die Zuordnung
- **Ressourcen-Label**: Ein vordefinierter Satz von Rechenanforderungen.
  Das nf-core-Projekt stellt Standard-Labels wie `process_single` für leichtgewichtige Tools und `process_high` für anspruchsvolle Tools bereit.
  Diese Labels helfen bei der Verwaltung der Ressourcenzuweisung über verschiedene Ausführungsumgebungen hinweg.
- **Metadaten-Anforderung**: Ob das Modul sample-spezifische Informationen über eine `meta`-Map benötigt (normalerweise ja für Datenverarbeitungsmodule).

Das Tool übernimmt die Komplexität der Suche nach Paketinformationen und der Einrichtung der Struktur, sodass du dich auf die Implementierung der spezifischen Logik des Tools konzentrieren kannst.

#### 2.1.2. Das Modul-Gerüst untersuchen

Das Tool erstellt eine vollständige Modulstruktur in `modules/local/` (oder `modules/nf-core/`, wenn du im nf-core/modules-Repository bist):

??? abstract "Verzeichnisinhalt"

    ```console
    modules/local/cowpy
    ├── environment.yml
    ├── main.nf
    ├── meta.yml
    └── tests
        └── main.nf.test
    ```

Jede Datei dient einem bestimmten Zweck:

- **`main.nf`**: Prozessdefinition mit allen eingebauten nf-core-Mustern
- **`meta.yml`**: Modul-Dokumentation, die Inputs, Outputs und das Tool beschreibt
- **`environment.yml`**: Conda-Umgebungsspezifikation für Abhängigkeiten
- **`tests/main.nf.test`**: nf-test-Testfälle zur Validierung, dass das Modul funktioniert

!!! tip "Mehr über Testing erfahren"

    Die generierte Testdatei verwendet nf-test, ein Testing-Framework für Nextflow-Pipelines und -Module. Um zu lernen, wie man diese Tests schreibt und ausführt, siehe die [nf-test Side Quest](../side_quests/nf-test.md).

Die generierte `main.nf` enthält alle Muster, die du gerade gelernt hast, plus einige zusätzliche Features:

```groovy title="modules/local/cowpy/main.nf" hl_lines="11 21 22"
process COWPY {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
        'biocontainers/YOUR-TOOL-HERE' }"

    input:
    tuple val(meta), path(input)        // Muster 1: Metadaten-Tupel ✓

    output:
    tuple val(meta), path("*"), emit: output
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''                // Muster 2: ext.args ✓
    def prefix = task.ext.prefix ?: "${meta.id}"  // Muster 3: ext.prefix ✓

    """
    // Füge hier deinen Tool-Befehl hinzu

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        COWPY: \$(cowpy --version)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo $args
    touch ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        COWPY: \$(cowpy --version)
    END_VERSIONS
    """
}
```

Beachte, wie alle Muster, die du oben manuell angewendet hast, bereits vorhanden sind!

Das Template enthält auch mehrere zusätzliche nf-core-Konventionen.
Einige davon funktionieren out-of-the-box, während andere Platzhalter sind, die wir ausfüllen müssen, wie unten beschrieben.

**Features, die sofort funktionieren:**

- **`tag "$meta.id"`**: Fügt Sample-ID zu Prozessnamen in Logs hinzu für einfacheres Tracking
- **`label 'process_single'`**: Ressourcen-Label zur Konfiguration von CPU/Speicher-Anforderungen
- **`when:`-Block**: Ermöglicht bedingte Ausführung über `task.ext.when`-Konfiguration

Diese Features sind bereits funktional und machen Module wartbarer.

**Platzhalter, die wir unten anpassen werden:**

- **`input:`- und `output:`-Blöcke**: Generische Deklarationen, die wir an unser Tool anpassen werden
- **`script:`-Block**: Enthält einen Kommentar, wo wir den `cowpy`-Befehl hinzufügen werden
- **`stub:`-Block**: Template, das wir aktualisieren werden, um die korrekten Ausgaben zu produzieren
- **Container und Umgebung**: Platzhalter, die wir mit Paketinformationen füllen werden

Die nächsten Abschnitte führen durch die Vervollständigung dieser Anpassungen.

### 2.2. Container und Conda-Umgebung einrichten

Die nf-core-Richtlinien erfordern, dass wir sowohl einen Container als auch eine Conda-Umgebung als Teil des Moduls spezifizieren.

#### 2.2.1. Container

Für den Container kannst du [Seqera Containers](https://seqera.io/containers/) verwenden, um automatisch einen Container aus jedem Conda-Paket zu erstellen, einschließlich conda-forge-Paketen.
In diesem Fall verwenden wir denselben vorgefertigten Container wie zuvor.

Der Standardcode bietet an, zwischen Docker und Singularity umzuschalten, aber wir werden diese Zeile vereinfachen und nur den Docker-Container angeben, den wir oben von Seqera Containers erhalten haben.

=== "Danach"

```groovy title="modules/local/cowpy/main.nf" linenums="3" hl_lines="6"
process COWPY {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273"
```

=== "Vorher"

```groovy title="modules/local/cowpy/main.nf" linenums="3" hl_lines="6"
process COWPY {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
        'biocontainers/YOUR-TOOL-HERE' }"
```

#### 2.2.2. Conda-Umgebung

Für die Conda-Umgebung spezifiziert der Modul-Code `conda "${moduleDir}/environment.yml"`, was bedeutet, dass sie in der `environment.yml`-Datei konfiguriert werden sollte.

Das Modul-Erstellungstool hat uns gewarnt, dass es das `cowpy`-Paket nicht in Bioconda finden konnte (dem primären Kanal für Bioinformatik-Tools).
Allerdings ist `cowpy` in conda-forge verfügbar, sodass du die `environment.yml` so vervollständigen kannst:

=== "Danach"

    ```yaml title="modules/local/cowpy/environment.yml"  linenums="1" hl_lines="1 3 5"
    name: COWPY
    channels:
      - conda-forge
    dependencies:
      - cowpy=1.1.5
    ```

=== "Vorher"

    ```yaml title="modules/local/cowpy/environment.yml" linenums="1"
    ---
    # yaml-language-server: $schema=https://raw.githubusercontent.com/nf-core/modules/master/modules/environment-schema.json
    channels:
      - conda-forge
      - bioconda
    dependencies:
      # TODO nf-core: Erforderliche Conda-Pakete auflisten.
      #               Software MUSS an Kanal (z.B. "bioconda"), Version (z.B. "1.10") gepinnt werden.
      #               Für Conda MUSS der Build (z.B. "h9402c20_2") AUSGESCHLOSSEN werden, um Installation auf verschiedenen Betriebssystemen zu unterstützen.
      - "YOUR-TOOL-HERE"
    ```

Für die Einreichung bei nf-core müssten wir den Standards genauer folgen, aber für unsere eigene Verwendung können wir den Code auf diese Weise vereinfachen.

!!! tip "Bioconda vs conda-forge-Pakete"

    - **Bioconda-Pakete**: Erhalten automatisch BioContainers, die gebrauchsfertige Container bereitstellen
    - **conda-forge-Pakete**: Können Seqera Containers verwenden, um Container on-demand aus dem Conda-Rezept zu erstellen

    Die meisten Bioinformatik-Tools sind in Bioconda, aber für conda-forge-Tools bietet Seqera Containers eine einfache Lösung für die Containerisierung.

### 2.3. Die `COWPY`-Logik einfügen

Jetzt aktualisieren wir die Code-Elemente, die spezifisch für das sind, was der `COWPY`-Prozess tut: die Inputs und Outputs und den script-Block.

#### 2.3.1. Inputs und Outputs

Das generierte Template enthält generische Input- und Output-Deklarationen, die du für dein spezifisches Tool anpassen musst.
Wenn wir auf unser manuelles `COWPY`-Modul aus Abschnitt 1 zurückblicken, können wir das als Leitfaden verwenden.

Aktualisiere die Input- und Output-Blöcke:

=== "Danach"

    ```groovy title="modules/local/cowpy/main.nf" linenums="8" hl_lines="2 5"
    input:
    tuple val(meta), path(input_file)

    output:
    tuple val(meta), path("${prefix}.txt"), emit: cowpy_output
    path "versions.yml"           , emit: versions
    ```

=== "Vorher"

    ```groovy title="modules/local/cowpy/main.nf" linenums="8" hl_lines="2 5"
    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path("*"), emit: output
    path "versions.yml"           , emit: versions
    ```

Dies spezifiziert:

- Den Input-Datei-Parameternamen (`input_file` statt generischem `input`)
- Den Ausgabedateinamen mit dem konfigurierbaren Präfix-Muster (`${prefix}.txt` statt Wildcard `*`)
- Einen beschreibenden emit-Namen (`cowpy_output` statt generischem `output`)

Wenn du den Nextflow-Language-Server zur Syntaxvalidierung verwendest, wird der `${prefix}`-Teil zu diesem Zeitpunkt als Fehler markiert, weil wir ihn noch nicht zum script-Block hinzugefügt haben.
Lass uns das jetzt tun.

#### 2.3.2. Der script-Block

Das Template bietet einen Kommentar-Platzhalter im script-Block, wo du den tatsächlichen Tool-Befehl hinzufügen solltest.

Basierend auf dem Modul, das wir früher manuell geschrieben haben, sollten wir folgende Änderungen vornehmen:

=== "Danach"

    ```groovy title="modules/local/cowpy/main.nf" linenums="15" hl_lines="3 6"
    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    cat $input_file | cowpy $args > ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        COWPY: \$(cowpy --version)
    END_VERSIONS
    """
    ```

=== "Vorher"

    ```groovy title="modules/local/cowpy/main.nf" linenums="15" hl_lines="6"
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    // Füge hier deinen Tool-Befehl hinzu

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        COWPY: \$(cowpy --version)
    END_VERSIONS
    """
    ```

Wichtige Änderungen:

- Ändere `def prefix` zu nur `prefix` (ohne `def`), um es im output-Block zugänglich zu machen
- Ersetze den Kommentar durch den tatsächlichen `cowpy`-Befehl, der sowohl `$args` als auch `${prefix}.txt` verwendet

Beachte, dass wir, wenn wir nicht bereits die Arbeit getan hätten, die `ext.args`- und `ext.prefix`-Konfiguration für den `COWPY`-Prozess zur `modules.config`-Datei hinzuzufügen, dies jetzt tun müssten.

#### 2.3.3. Den stub-Block implementieren

Im Nextflow-Kontext ermöglicht ein [stub](https://www.nextflow.io/docs/latest/process.html#stub)-Block die Definition eines leichtgewichtigen Dummy-Skripts, das für schnelles Prototyping und Testen der Logik einer Pipeline verwendet wird, ohne den tatsächlichen Befehl auszuführen.

Mach dir keine Sorgen, wenn das mysteriös erscheint; wir fügen dies der Vollständigkeit halber ein, aber du kannst den stub-Abschnitt auch einfach löschen, wenn du dich nicht damit befassen möchtest, da er völlig optional ist.

=== "Danach"

    ```groovy title="modules/local/cowpy/main.nf" linenums="27" hl_lines="3 6"
    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        COWPY: \$(cowpy --version)
    END_VERSIONS
    """
    ```

=== "Vorher"

    ```groovy title="modules/local/cowpy/main.nf" linenums="27" hl_lines="3 6"
    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        COWPY: \$(cowpy --version)
    END_VERSIONS
    """
    ```

Wichtige Änderungen:

- Ändere `def prefix` zu nur `prefix`, um dem script-Block zu entsprechen
- Entferne die `echo $args`-Zeile (die nur Template-Platzhalter-Code war)
- Der stub erstellt eine leere `${prefix}.txt`-Datei, die dem entspricht, was der script-Block produziert

Dies ermöglicht es dir, Workflow-Logik und Dateihandhabung zu testen, ohne auf die tatsächliche Ausführung des Tools zu warten.

Sobald du die Umgebungseinrichtung (Abschnitt 2.2), Inputs/Outputs (Abschnitt 2.3.1), script-Block (Abschnitt 2.3.2) und stub-Block (Abschnitt 2.3.3) abgeschlossen hast, ist das Modul bereit zum Testen!

### 2.4. Das neue `COWPY`-Modul einsetzen und die Pipeline ausführen

Alles, was wir tun müssen, um diese neue Version des `COWPY`-Moduls auszuprobieren, ist die Import-Anweisung in der `hello.nf`-Workflow-Datei zu ändern, um auf die neue Datei zu zeigen.

=== "Danach"

    ```groovy title="workflows/hello.nf" linenums="1" hl_lines="10"
    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        MODULE / SUBWORKFLOWS / FUNKTIONEN IMPORTIEREN
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    include { paramsSummaryMap       } from 'plugin/nf-schema'
    include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
    include { sayHello               } from '../modules/local/sayHello.nf'
    include { convertToUpper         } from '../modules/local/convertToUpper.nf'
    include { COWPY                  } from '../modules/local/cowpy/main.nf'
    include { CAT_CAT                } from '../modules/nf-core/cat/cat/main'
    ```

=== "Vorher"

    ```groovy title="modules/local/cowpy/main.nf" linenums="1" hl_lines="10"
    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        MODULE / SUBWORKFLOWS / FUNKTIONEN IMPORTIEREN
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    include { paramsSummaryMap       } from 'plugin/nf-schema'
    include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
    include { sayHello               } from '../modules/local/sayHello.nf'
    include { convertToUpper         } from '../modules/local/convertToUpper.nf'
    include { COWPY                  } from '../modules/local/cowpy.nf'
    include { CAT_CAT                } from '../modules/nf-core/cat/cat/main'
    ```

Lass uns die Pipeline ausführen, um sie zu testen.

```bash
nextflow run . --outdir core-hello-results -profile test,docker --validate_params false
```

??? success "Befehlsausgabe"

    ```console hl_lines="33"
      N E X T F L O W   ~  version 25.04.3

    Launching `./main.nf` [prickly_neumann] DSL2 - revision: b9e9b3b8de

    Input/output options
      input                     : /workspaces/training/hello-nf-core/core-hello/assets/greetings.csv
      outdir                    : core-hello-results

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      validate_params           : false
      trace_report_suffix       : 2025-12-27_08-23-51

    Core Nextflow options
      runName                   : prickly_neumann
      containerEngine           : docker
      launchDir                 : /workspaces/training/hello-nf-core/core-hello
      workDir                   : /workspaces/training/hello-nf-core/core-hello/work
      projectDir                : /workspaces/training/hello-nf-core/core-hello
      userName                  : root
      profile                   : test,docker
      configFiles               : /workspaces/training/hello-nf-core/core-hello/nextflow.config

    !! Only displaying parameters that differ from the pipeline defaults !!
    ------------------------------------------------------
    executor >  local (8)
    [e9/008ede] CORE_HELLO:HELLO:sayHello (3)       [100%] 3 of 3 ✔
    [f0/d70cfe] CORE_HELLO:HELLO:convertToUpper (3) [100%] 3 of 3 ✔
    [be/0ecc58] CORE_HELLO:HELLO:CAT_CAT (test)     [100%] 1 of 1 ✔
    [11/8e082f] CORE_HELLO:HELLO:COWPY (test)       [100%] 1 of 1 ✔
    -[core/hello] Pipeline completed successfully-
    ```

Dies produziert die gleichen Ergebnisse wie zuvor.

### Fazit

Du weißt jetzt, wie du die eingebauten nf-core-Tools verwendest, um Module effizient mit Templates zu erstellen, anstatt alles von Grund auf zu schreiben.

### Wie geht es weiter?

Lerne, welche Vorteile es hat, Module zu nf-core beizutragen, und welche Hauptschritte und Anforderungen damit verbunden sind.

---

## 3. Module zu nf-core beitragen

Das [nf-core/modules](https://github.com/nf-core/modules)-Repository begrüßt Beiträge von gut getesteten, standardisierten Modulen.

### 3.1. Warum beitragen?

Das Beitragen deiner Module zu nf-core:

- Macht deine Tools für die gesamte nf-core-Community über den Modul-Katalog unter [nf-co.re/modules](https://nf-co.re/modules) verfügbar
- Stellt laufende Community-Wartung und Verbesserungen sicher
- Bietet Qualitätssicherung durch Code-Review und automatisiertes Testing
- Gibt deiner Arbeit Sichtbarkeit und Anerkennung

### 3.2. Checkliste für Beitragende

Um ein Modul zu nf-core beizutragen, musst du folgende Schritte durchlaufen:

1. Prüfe, ob es bereits unter [nf-co.re/modules](https://nf-co.re/modules) existiert
2. Forke das [nf-core/modules](https://github.com/nf-core/modules)-Repository
3. Verwende `nf-core modules create`, um das Template zu generieren
4. Fülle die Modul-Logik und Tests aus
5. Teste mit `nf-core modules test tool/subtool`
6. Linte mit `nf-core modules lint tool/subtool`
7. Reiche einen Pull Request ein

Für detaillierte Anweisungen siehe das [nf-core Components Tutorial](https://nf-co.re/docs/tutorials/nf-core_components/components).

### 3.3. Ressourcen

- **Components Tutorial**: [Vollständiger Leitfaden zur Erstellung und zum Beitragen von Modulen](https://nf-co.re/docs/tutorials/nf-core_components/components)
- **Modul-Spezifikationen**: [Technische Anforderungen und Richtlinien](https://nf-co.re/docs/guidelines/components/modules)
- **Community-Support**: [nf-core Slack](https://nf-co.re/join) - Tritt dem `#modules`-Kanal bei

### Fazit

Du weißt jetzt, wie du nf-core-Module erstellst! Du hast die vier Schlüsselmuster gelernt, die Module portabel und wartbar machen:

- **Metadaten-Tupel** propagieren Metadaten durch den Workflow
- **`ext.args`** vereinfacht Modul-Schnittstellen, indem optionale Argumente über die Konfiguration behandelt werden
- **`ext.prefix`** standardisiert die Ausgabedatei-Benennung
- **Zentralisiertes Publishing** über `publishDir`, konfiguriert in `modules.config` statt fest codiert in Modulen

Durch die schrittweise Transformation von `COWPY` hast du ein tiefes Verständnis dieser Muster entwickelt, was dich befähigt, mit nf-core-Modulen zu arbeiten, sie zu debuggen und zu erstellen.
In der Praxis wirst du `nf-core modules create` verwenden, um korrekt strukturierte Module mit diesen Mustern von Anfang an zu generieren.

Schließlich hast du gelernt, wie du Module zur nf-core-Community beiträgst, wodurch Tools für Forscher\*innen weltweit verfügbar werden, während du von laufender Community-Wartung profitierst.

### Wie geht es weiter?

Wenn du bereit bist, fahre fort mit [Teil 5: Input-Validierung](./05_input_validation.md), um zu lernen, wie du schema-basierte Input-Validierung zu deiner Pipeline hinzufügst.
