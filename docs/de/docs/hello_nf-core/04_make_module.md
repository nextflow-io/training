# Teil 4: Ein nf-core-Modul erstellen

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

In diesem vierten Teil des Hello nf-core Trainingskurses zeigen wir dir, wie du ein nf-core-Modul erstellst, indem du die wichtigsten Konventionen anwendest, die Module portabel und wartbar machen.

Das nf-core-Projekt stellt einen Befehl (`nf-core modules create`) bereit, der automatisch korrekt strukturierte Modulvorlagen generiert, ähnlich wie bei dem, was wir für den Workflow in Teil 2 verwendet haben.
Zu Lehrzwecken werden wir jedoch damit beginnen, es manuell zu machen: Wir transformieren das lokale `cowpy`-Modul in deiner `core-hello`-Pipeline Schritt für Schritt in ein Modul im nf-core-Stil.
Danach zeigen wir dir, wie du die vorlagenbasierte Modulerstellung nutzt, um in Zukunft effizienter zu arbeiten.

??? info "Wie du von diesem Abschnitt aus beginnst"

    Dieser Abschnitt setzt voraus, dass du [Teil 3: Ein nf-core-Modul verwenden](./03_use_module.md) abgeschlossen hast und das `CAT_CAT`-Modul in deine Pipeline integriert hast.

    Wenn du Teil 3 nicht abgeschlossen hast oder für diesen Teil neu beginnen möchtest, kannst du die Lösung `core-hello-part3` als Ausgangspunkt verwenden.
    Führe diese Befehle im Verzeichnis `hello-nf-core/` aus:

    ```bash
    cp -r solutions/core-hello-part3 core-hello
    cd core-hello
    ```

    Dies gibt dir eine Pipeline mit dem bereits integrierten `CAT_CAT`-Modul.
    Du kannst testen, ob sie erfolgreich läuft, indem du den folgenden Befehl ausführst:

    ```bash
    nextflow run . --outdir core-hello-results -profile test,docker --validate_params false
    ```

---

## 1. `cowpy` in ein nf-core-Modul transformieren

In diesem Abschnitt wenden wir nf-core-Konventionen auf das lokale `cowpy`-Modul in deiner `core-hello`-Pipeline an und transformieren es in ein Modul, das den Community-Standards von nf-core folgt.

Dies ist der aktuelle Code für das `cowpy`-Prozessmodul:

```groovy title="core-hello/modules/local/cowpy.nf" linenums="1"
#!/usr/bin/env nextflow

// ASCII-Kunst mit cowpy generieren (https://github.com/jeffbuttars/cowpy)
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

Wir werden die folgenden nf-core-Konventionen schrittweise anwenden:

1. **Den Prozessnamen zu `COWPY` in Großbuchstaben ändern**, um der Konvention zu folgen.
2. **`COWPY` aktualisieren, um Metadaten-Tupel zu verwenden**, um Probenmetadaten durch den Workflow zu propagieren.
3. **Tool-Argument-Konfiguration mit `ext.args` zentralisieren**, um die Vielseitigkeit des Moduls zu erhöhen und gleichzeitig die Schnittstelle minimal zu halten.
4. **Ausgabe-Benennung mit `ext.prefix` standardisieren**, um Konsistenz zu fördern.
5. **Die Publishing-Konfiguration zentralisieren**, um Konsistenz zu fördern.

Nach jedem Schritt führen wir die Pipeline aus, um zu testen, dass alles wie erwartet funktioniert.

!!! warning "Arbeitsverzeichnis"

    Stelle sicher, dass du dich im Verzeichnis `core-hello` (dein Pipeline-Stammverzeichnis) befindest, für alle Dateibearbeitungen und Befehlsausführungen in diesem Abschnitt.

    ```bash
    cd core-hello
    ```

### 1.1. Den Prozessnamen in Großbuchstaben ändern

Dies ist rein eine stilistische Konvention (es gibt keine technische Rechtfertigung), aber da es die Norm für nf-core-Module ist, sollten wir uns daran halten.

Wir müssen drei Änderungssets vornehmen:

1. Den Prozessnamen im Modul aktualisieren
2. Die Modul-Import-Anweisung im Workflow-Header aktualisieren
3. Den Prozessaufruf und die Emit-Deklaration im Workflow-Body aktualisieren

Los geht's!

#### 1.1.1. Den Prozessnamen im Modul aktualisieren

Öffne die Moduldatei `cowpy.nf` (unter `core-hello/modules/local/`) und ändere den Prozessnamen in Großbuchstaben:

=== "Nachher"

    ```groovy title="core-hello/modules/local/cowpy.nf" linenums="3" hl_lines="2"
    // ASCII-Kunst mit cowpy generieren (https://github.com/jeffbuttars/cowpy)
    process COWPY {
    ```

=== "Vorher"

    ```groovy title="core-hello/modules/local/cowpy.nf" linenums="3" hl_lines="2"
    // ASCII-Kunst mit cowpy generieren (https://github.com/jeffbuttars/cowpy)
    process cowpy {
    ```

In diesem Fall ist die Großschreibung völlig unkompliziert.

Wenn der Prozessname aus mehreren Wörtern bestünde, zum Beispiel wenn wir einen Prozess namens MyCowpyTool ursprünglich in CamelCase hätten, wäre die nf-core-Konvention, Unterstriche zu verwenden, um sie zu trennen, was MY_COWPY_TOOL ergäbe.

#### 1.1.2. Die Modul-Import-Anweisung aktualisieren

Prozessnamen sind case-sensitiv, daher müssen wir, nachdem wir den Prozessnamen geändert haben, die Modul-Import-Anweisung entsprechend im Workflow-Header von `hello.nf` aktualisieren:

=== "Nachher"

    ```groovy title="core-hello/workflows/hello.nf" linenums="1" hl_lines="10"
    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
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
        IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    include { paramsSummaryMap       } from 'plugin/nf-schema'
    include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
    include { sayHello               } from '../modules/local/sayHello.nf'
    include { convertToUpper         } from '../modules/local/convertToUpper.nf'
    include { cowpy                  } from '../modules/local/cowpy/main.nf'
    include { CAT_CAT                } from '../modules/nf-core/cat/cat/main'
    ```

Wir könnten einen Alias in der Import-Anweisung verwenden, um zu vermeiden, dass wir Aufrufe des Prozesses aktualisieren müssen, aber das würde den Sinn der Übernahme der Großbuchstaben-Konvention etwas verfehlen.

#### 1.1.3. Den Prozessaufruf und die Emit-Deklaration aktualisieren

Also aktualisieren wir jetzt die beiden Referenzen auf den Prozess im Workflow-Block von `hello.nf`:

=== "Nachher"

    ```groovy title="core-hello/workflows/hello.nf" linenums="43" hl_lines="2 17"
    // ASCII-Kunst der Begrüßungen mit cowpy generieren
    COWPY(CAT_CAT.out.file_out)

    //
    // Software-Versionen sammeln und speichern
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
    // ASCII-Kunst der Begrüßungen mit cowpy generieren
    cowpy(CAT_CAT.out.file_out)

    //
    // Software-Versionen sammeln und speichern
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

Gut, das funktioniert! Lass uns nun zu substanzielleren Änderungen übergehen.

### 1.2. `COWPY` aktualisieren, um Metadaten-Tupel zu verwenden

In der aktuellen Version der `core-hello`-Pipeline extrahieren wir die Datei aus dem Ausgabe-Tupel von `CAT_CAT`, um sie an `COWPY` zu übergeben, wie in der oberen Hälfte des Diagramms unten gezeigt.

<figure class="excalidraw">
    --8<-- "docs/en/docs/hello_nf-core/img/cowpy-inputs.svg"
</figure>

Es wäre besser, wenn `COWPY` Metadaten-Tupel direkt akzeptieren würde, sodass Metadaten durch den Workflow fließen können, wie in der unteren Hälfte des Diagramms gezeigt.

Zu diesem Zweck müssen wir die folgenden Änderungen vornehmen:

1. Die Eingabe- und Ausgabedefinitionen aktualisieren
2. Den Prozessaufruf im Workflow aktualisieren
3. Den Emit-Block im Workflow aktualisieren

Nachdem wir all das getan haben, führen wir die Pipeline aus, um zu testen, dass alles noch wie zuvor funktioniert.

#### 1.2.1. Die Eingabe- und Ausgabedefinitionen aktualisieren

Kehre zur Moduldatei `cowpy.nf` zurück und ändere sie so, dass sie Metadaten-Tupel akzeptiert, wie unten gezeigt.

=== "Nachher"

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

Wie du sehen kannst, haben wir sowohl die **Haupteingabe** als auch die **Ausgabe** in ein Tupel geändert, das dem Muster `tuple val(meta), path(input_file)` folgt, das in Teil 3 dieses Trainings eingeführt wurde.
Für die Ausgabe haben wir auch diese Gelegenheit genutzt, um `emit: cowpy_output` hinzuzufügen, um dem Ausgabekanal einen beschreibenden Namen zu geben.

Jetzt, da wir geändert haben, was der Prozess erwartet, müssen wir aktualisieren, was wir ihm im Prozessaufruf bereitstellen.

#### 1.2.2. Den Prozessaufruf im Workflow aktualisieren

Die gute Nachricht ist, dass diese Änderung den Prozessaufruf vereinfachen wird.
Da die Ausgabe von `CAT_CAT` und die Eingabe von `COWPY` jetzt die gleiche 'Form' haben, d.h. beide aus einer `tuple val(meta), path(input_file)`-Struktur bestehen, können wir sie einfach direkt verbinden, anstatt die Datei explizit aus der Ausgabe des `CAT_CAT`-Prozesses extrahieren zu müssen.

Öffne die Workflow-Datei `hello.nf` (unter `core-hello/workflows/`) und aktualisiere den Aufruf von `COWPY` wie unten gezeigt.

=== "Nachher"

    ```groovy title="core-hello/workflows/hello.nf" linenums="43" hl_lines="2"
        // ASCII-Kunst der Begrüßungen mit cowpy generieren
        COWPY(CAT_CAT.out.file_out, params.character)
    ```

=== "Vorher"

    ```groovy title="core-hello/workflows/hello.nf" linenums="43" hl_lines="1-2 5"
        // Die Datei aus dem Tupel extrahieren, da cowpy noch keine Metadaten verwendet
        ch_for_cowpy = CAT_CAT.out.file_out.map{ meta, file -> file }

        // ASCII-Kunst der Begrüßungen mit cowpy generieren
        COWPY(ch_for_cowpy, params.character)
    ```

Wir rufen jetzt `COWPY` direkt mit `CAT_CAT.out.file_out` auf.

Dadurch müssen wir den `ch_for_cowpy`-Kanal nicht mehr konstruieren, sodass diese Zeile (und ihre Kommentarzeile) vollständig gelöscht werden kann.

#### 1.2.3. Den Emit-Block im Workflow aktualisieren

Da `COWPY` jetzt eine benannte Ausgabe, `cowpy_output`, ausgibt, können wir den `emit:`-Block des `hello.nf`-Workflows aktualisieren, um diese zu verwenden.

=== "Nachher"

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

Dies ist technisch nicht erforderlich, aber es ist gute Praxis, auf benannte Ausgaben zu verweisen, wann immer möglich.

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

Die Pipeline sollte erfolgreich laufen, wobei Metadaten nun von `CAT_CAT` durch `COWPY` fließen.

Damit ist abgeschlossen, was wir tun mussten, um `COWPY` Metadaten-Tupel verarbeiten zu lassen.
Schauen wir uns nun an, was wir sonst noch tun können, um die nf-core-Modulmuster zu nutzen.

### 1.3. Tool-Argument-Konfiguration mit `ext.args` zentralisieren

In seinem aktuellen Zustand erwartet der `COWPY`-Prozess, einen Wert für den `character`-Parameter zu erhalten.
Daher müssen wir jedes Mal, wenn wir den Prozess aufrufen, einen Wert angeben, auch wenn wir mit den vom Tool gesetzten Standardwerten zufrieden wären.
Für `COWPY` ist dies zugegebenermaßen kein großes Problem, aber für Tools mit vielen optionalen Parametern kann es ziemlich umständlich werden.

Das nf-core-Projekt empfiehlt die Verwendung einer Nextflow-Funktion namens [`ext.args`](https://www.nextflow.io/docs/latest/reference/process.html#ext), um Tool-Argumente bequemer über Konfigurationsdateien zu verwalten.

Anstatt Prozesseingaben für jede Tool-Option zu deklarieren, schreibst du das Modul so, dass es `ext.args` bei der Konstruktion seiner Befehlszeile referenziert.
Dann ist es nur noch eine Frage der Einrichtung der `ext.args`-Variable, um die Argumente und Werte zu halten, die du in der Datei `modules.config` verwenden möchtest, die Konfigurationsdetails für alle Module konsolidiert.
Nextflow fügt diese Argumente mit ihren Werten zur Laufzeit in die Tool-Befehlszeile ein.

Lass uns diesen Ansatz auf das `COWPY`-Modul anwenden.
Wir müssen die folgenden Änderungen vornehmen:

1. Das `COWPY`-Modul aktualisieren
2. `ext.args` in der Datei `modules.config` konfigurieren
3. Den `hello.nf`-Workflow aktualisieren

Nachdem wir all das getan haben, führen wir die Pipeline aus, um zu testen, dass alles noch wie zuvor funktioniert.

#### 1.3.1. Das `COWPY`-Modul aktualisieren

Lass es uns tun.
Öffne die Moduldatei `cowpy.nf` (unter `core-hello/modules/local/`) und ändere sie so, dass sie auf `ext.args` verweist, wie unten gezeigt.

=== "Nachher"

    ```groovy title="modules/local/cowpy.nf" linenums="1" hl_lines="18 20"
    #!/usr/bin/env nextflow

    // ASCII-Kunst mit cowpy generieren (https://github.com/jeffbuttars/cowpy)
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

    // ASCII-Kunst mit cowpy generieren (https://github.com/jeffbuttars/cowpy)
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

1. **Im `input:`-Block haben wir die Eingabe `val character` entfernt.**
   Von nun an werden wir dieses Argument über die `ext.args`-Konfiguration bereitstellen, wie weiter unten beschrieben.

2. **Im `script:`-Block haben wir die Zeile `def args = task.ext.args ?: ''` hinzugefügt.**
   Diese Zeile verwendet den `?:`-Operator, um den Wert der `args`-Variable zu bestimmen: den Inhalt von `task.ext.args`, wenn er nicht leer ist, oder eine leere Zeichenkette, wenn er es ist.
   Beachte, dass wir zwar allgemein auf `ext.args` verweisen, dieser Code jedoch `task.ext.args` referenzieren muss, um die `ext.args`-Konfiguration auf Modulebene herauszuziehen.

3. **In der Befehlszeile haben wir `-c "$character"` durch `$args` ersetzt.**
   Hier fügt Nextflow alle Tool-Argumente ein, die in `ext.args` in der Datei `modules.config` gesetzt sind.

Dadurch ist die Modulschnittstelle jetzt einfacher: Sie erwartet nur die wesentlichen Metadaten- und Dateieingaben.

!!! note

    Der `?:`-Operator wird oft 'Elvis-Operator' genannt, weil er wie ein seitwärts gedrehtes Elvis-Presley-Gesicht aussieht, wobei das `?`-Zeichen die Welle in seinen Haaren symbolisiert.

#### 1.3.2. `ext.args` in der Datei `modules.config` konfigurieren

Nachdem wir die `character`-Deklaration aus dem Modul entfernt haben, müssen wir sie zu `ext.args` in der Konfigurationsdatei `modules.config` hinzufügen.

Konkret werden wir dieses kleine Code-Stück zum `process {}`-Block hinzufügen:

```groovy title="Hinzuzufügender Code"
withName: 'COWPY' {
    ext.args = { "-c ${params.character}" }
}
```

Die `withName:`-Syntax weist diese Konfiguration nur dem `COWPY`-Prozess zu, und `ext.args = { "-c ${params.character}" }` komponiert einfach eine Zeichenkette, die den Wert des `character`-Parameters enthält.
Beachte die Verwendung der geschweiften Klammern, die Nextflow anweisen, den Wert des Parameters zur Laufzeit zu evaluieren.

Macht das Sinn? Lass es uns hinzufügen.

Öffne `conf/modules.config` und füge den Konfigurationscode innerhalb des `process {}`-Blocks hinzu, wie unten gezeigt.

=== "Nachher"

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

Hoffentlich kannst du dir vorstellen, dass alle Module in einer Pipeline ihre `ext.args` in dieser Datei spezifiziert haben, mit den folgenden Vorteilen:

- Die **Modulschnittstelle bleibt einfach** - Sie akzeptiert nur die wesentlichen Metadaten- und Dateieingaben
- Die **Pipeline stellt immer noch `params.character` bereit** - Endbenutzer können es weiterhin wie zuvor konfigurieren
- Das **Modul ist jetzt portabel** - Es kann in anderen Pipelines wiederverwendet werden, ohne einen bestimmten Parameternamen zu erwarten
- Die Konfiguration ist **zentralisiert** in `modules.config`, wodurch die Workflow-Logik sauber bleibt

Indem wir die Datei `modules.config` als den Ort verwenden, an dem alle Pipelines die Konfiguration pro Modul zentralisieren, machen wir unsere Module über verschiedene Pipelines hinweg wiederverwendbarer.

#### 1.3.3. Den `hello.nf`-Workflow aktualisieren

Da das `COWPY`-Modul den `character`-Parameter nicht mehr als Eingabe benötigt, müssen wir den Workflow-Aufruf entsprechend aktualisieren.

Öffne die Workflow-Datei `hello.nf` (unter `core-hello/workflows/`) und aktualisiere den Aufruf von `COWPY` wie unten gezeigt.

=== "Nachher"

    ```groovy title="core-hello/workflows/hello.nf" linenums="39" hl_lines="2"
        // ASCII-Kunst der Begrüßungen mit cowpy generieren
        COWPY(CAT_CAT.out.file_out)
    ```

=== "Vorher"

    ```groovy title="core-hello/workflows/hello.nf" linenums="39" hl_lines="2"
        // ASCII-Kunst der Begrüßungen mit cowpy generieren
        COWPY(CAT_CAT.out.file_out, params.character)
    ```

Der Workflow-Code ist jetzt sauberer: Wir müssen `params.character` nicht direkt an den Prozess übergeben.
Die Modulschnittstelle wird minimal gehalten, was sie portabler macht, während die Pipeline die explizite Option weiterhin über die Konfiguration bereitstellt.

#### 1.3.4. Die Pipeline ausführen, um sie zu testen

Lass uns testen, dass der Workflow immer noch wie erwartet funktioniert, indem wir einen anderen Charakter angeben, um zu überprüfen, dass die `ext.args`-Konfiguration funktioniert.

Führe diesen Befehl mit `kosh` aus, einer der... rätselhafteren Optionen:

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
Finde die Ausgabe im Dateibrowser oder verwende den Aufgaben-Hash (der Teil `38/eb29ea` im obigen Beispiel), um die Ausgabedatei anzusehen:

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

    Du wirst den `cowpy`-Befehl mit dem Argument `-c kosh` sehen:

    ```console
    #!/usr/bin/env bash
    ...
    cat test.txt | cowpy -c kosh > cowpy-test.txt
    ```

    Dies zeigt, dass die `.command.sh`-Datei basierend auf der `ext.args`-Konfiguration korrekt generiert wurde.

Nimm dir einen Moment Zeit, um darüber nachzudenken, was wir hier erreicht haben.
Dieser Ansatz hält die Modulschnittstelle auf wesentliche Daten (Dateien, Metadaten und alle obligatorischen probenbezogenen Parameter) fokussiert, während Optionen, die das Verhalten des Tools steuern, separat über die Konfiguration gehandhabt werden.

Dies mag für ein einfaches Tool wie `cowpy` unnötig erscheinen, aber es kann für Datenanalyse-Tools mit vielen optionalen Argumenten einen großen Unterschied machen.

Zusammenfassung der Vorteile dieses Ansatzes:

- **Saubere Schnittstelle**: Das Modul konzentriert sich auf wesentliche Dateneingaben (Metadaten und Dateien)
- **Flexibilität**: Tool-Argumente können über die Konfiguration angegeben werden, einschließlich probenspezifischer Werte
- **Konsistenz**: Alle nf-core-Module folgen diesem Muster
- **Portabilität**: Module können ohne fest codierte Tool-Optionen wiederverwendet werden
- **Keine Workflow-Änderungen**: Das Hinzufügen oder Ändern von Tool-Optionen erfordert keine Aktualisierung des Workflow-Codes

!!! note

    Das `ext.args`-System hat leistungsstarke zusätzliche Funktionen, die hier nicht behandelt werden, einschließlich des dynamischen Wechselns von Argumentwerten basierend auf Metadaten. Siehe die [nf-core-Modulspezifikationen](https://nf-co.re/docs/guidelines/components/modules) für weitere Details.

### 1.4. Ausgabe-Benennung mit `ext.prefix` standardisieren

Nachdem wir dem `COWPY`-Prozess nun Zugriff auf die Metamap gegeben haben, können wir ein weiteres nützliches nf-core-Muster nutzen: die Benennung von Ausgabedateien basierend auf Metadaten.

Hier werden wir eine Nextflow-Funktion namens `ext.prefix` verwenden, die es uns ermöglicht, die Benennung von Ausgabedateien über Module hinweg mit `meta.id` (der Kennung, die in der Metamap enthalten ist) zu standardisieren, während wir Module bei Bedarf weiterhin individuell konfigurieren können.

Dies wird ähnlich sein wie das, was wir mit `ext.args` gemacht haben, mit ein paar Unterschieden, die wir im Laufe der Zeit detailliert erläutern werden.

Lass uns diesen Ansatz auf das `COWPY`-Modul anwenden.
Wir müssen die folgenden Änderungen vornehmen:

1. Das `COWPY`-Modul aktualisieren
2. `ext.prefix` in der Datei `modules.config` konfigurieren

(Keine Änderungen am Workflow erforderlich.)

Nachdem wir das getan haben, führen wir die Pipeline aus, um zu testen, dass alles noch wie zuvor funktioniert.

#### 1.4.1. Das `COWPY`-Modul aktualisieren

Öffne die Moduldatei `cowpy.nf` (unter `core-hello/modules/local/`) und ändere sie so, dass sie auf `ext.prefix` verweist, wie unten gezeigt.

=== "Nachher"

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
   Diese Zeile verwendet den `?:`-Operator, um den Wert der `prefix`-Variable zu bestimmen: den Inhalt von `task.ext.prefix`, wenn er nicht leer ist, oder die Kennung aus der Metamap (`meta.id`), wenn er es ist.
   Beachte, dass wir zwar allgemein auf `ext.prefix` verweisen, dieser Code jedoch `task.ext.prefix` referenzieren muss, um die `ext.prefix`-Konfiguration auf Modulebene herauszuziehen.

2. **In der Befehlszeile haben wir `cowpy-${input_file}` durch `${prefix}.txt` ersetzt.**
   Hier fügt Nextflow den durch die obige Zeile bestimmten Wert von `prefix` ein.

3. **Im `output:`-Block haben wir `path("cowpy-${input_file}")` durch `path("${prefix}.txt")` ersetzt.**
   Dies wiederholt einfach, wie der Dateipfad gemäß dem in der Befehlszeile Geschriebenen sein wird.

Dadurch wird der Name der Ausgabedatei nun mit einem sinnvollen Standard (der Kennung aus der Metamap) kombiniert mit der entsprechenden Dateiformaterweiterung konstruiert.

#### 1.4.2. `ext.prefix` in der Datei `modules.config` konfigurieren

In diesem Fall ist der sinnvolle Standard nach unserem Geschmack nicht ausreichend ausdrucksstark; wir möchten ein benutzerdefiniertes Benennungsmuster verwenden, das den Tool-Namen enthält, `cowpy-<id>.txt`, wie wir es zuvor hatten.

Wir tun dies, indem wir `ext.prefix` in `modules.config` konfigurieren, genau wie wir es für den `character`-Parameter mit `ext.args` getan haben, außer dass diesmal der Block `withName: 'COWPY' {}` bereits existiert und wir nur die folgende Zeile hinzufügen müssen:

```groovy title="Hinzuzufügender Code"
ext.prefix = { "cowpy-${meta.id}" }
```

Dies wird die gewünschte Zeichenkette komponieren.
Beachte, dass wir erneut geschweifte Klammern verwenden, diesmal um Nextflow anzuweisen, den Wert von `meta.id` zur Laufzeit zu evaluieren.

Lass es uns hinzufügen.

Öffne `conf/modules.config` und füge den Konfigurationscode innerhalb des `process {}`-Blocks hinzu, wie unten gezeigt.

=== "Nachher"

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

Lass uns testen, dass der Workflow immer noch wie erwartet funktioniert.

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

Schaue dir die Ausgabe im Ergebnisverzeichnis an.
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

Ändere gerne die `ext.prefix`-Konfiguration in `conf/modules.config`, um dich zu überzeugen, dass du das Benennungsmuster ändern kannst, ohne Änderungen am Modul- oder Workflow-Code vornehmen zu müssen.

Alternativ kannst du dies auch erneut mit einem anderen auf der Befehlszeile angegebenen `--batch`-Parameter ausführen, um dich zu überzeugen, dass dieser Teil weiterhin spontan anpassbar ist.

Dies zeigt, wie `ext.prefix` es dir ermöglicht, deine bevorzugte Benennungskonvention beizubehalten und gleichzeitig die Modulschnittstelle flexibel zu halten.

Zusammenfassung der Vorteile dieses Ansatzes:

- **Standardisierte Benennung**: Ausgabedateien werden typischerweise mit Proben-IDs aus Metadaten benannt
- **Konfigurierbar**: Die Standardbenennung kann bei Bedarf überschrieben werden
- **Konsistent**: Alle nf-core-Module folgen diesem Muster
- **Vorhersagbar**: Es ist einfach zu wissen, wie Ausgabedateien heißen werden

Ziemlich gut, oder?
Nun, es gibt noch eine wichtige Änderung, die wir vornehmen müssen, um unser Modul an die nf-core-Richtlinien anzupassen.

### 1.5. Die Publishing-Konfiguration zentralisieren

Du hast vielleicht bemerkt, dass wir Ausgaben in zwei verschiedene Verzeichnisse veröffentlicht haben:

- **`results`** — Das ursprüngliche Ausgabeverzeichnis, das wir von Anfang an für unsere lokalen Module verwendet haben, individuell mit per-Modul `publishDir`-Direktiven gesetzt;
- **`core-hello-results`** — Das Ausgabeverzeichnis, das mit `--outdir` auf der Befehlszeile gesetzt wurde, das die nf-core-Logs und die von `CAT_CAT` veröffentlichten Ergebnisse erhalten hat.

Dies ist unordentlich und suboptimal; es wäre besser, einen Ort für alles zu haben.
Natürlich könnten wir in jedes unserer lokalen Module gehen und die `publishDir`-Direktive manuell aktualisieren, um das Verzeichnis `core-hello-results` zu verwenden, aber was ist beim nächsten Mal, wenn wir das Ausgabeverzeichnis ändern?

Einzelne Module Publishing-Entscheidungen treffen zu lassen, ist eindeutig nicht der richtige Weg, insbesondere in einer Welt, in der dasselbe Modul in vielen verschiedenen Pipelines von Personen verwendet werden könnte, die unterschiedliche Bedürfnisse oder Präferenzen haben.
Wir möchten in der Lage sein, auf Ebene der Workflow-Konfiguration zu kontrollieren, wo Ausgaben veröffentlicht werden.

"Hey", könntest du sagen, "`CAT_CAT` sendet seine Ausgaben an `--outdir`. Vielleicht sollten wir dessen `publishDir`-Direktive kopieren?"

Ja, das ist eine großartige Idee.

Außer dass es keine `publishDir`-Direktive hat. (Schau mal, sieh dir den Modulcode an.)

Das liegt daran, dass nf-core-Pipelines die Kontrolle auf Workflow-Ebene zentralisieren, indem sie `publishDir` in `conf/modules.config` konfigurieren, anstatt in einzelnen Modulen.
Konkret deklariert die nf-core-Vorlage eine Standard-`publishDir`-Direktive (mit einer vordefinierten Verzeichnisstruktur), die für alle Module gilt, es sei denn, eine überschreibende Direktive wird bereitgestellt.

Klingt das nicht großartig? Könnte es sein, dass wir, um diese Standarddirektive zu nutzen, nur die aktuelle `publishDir`-Direktive aus unseren lokalen Modulen entfernen müssen?

Lass uns das bei `COWPY` ausprobieren, um zu sehen, was passiert, dann schauen wir uns den Code für die Standardkonfiguration an, um zu verstehen, wie sie funktioniert.

Schließlich werden wir demonstrieren, wie man das Standardverhalten bei Bedarf überschreibt.

#### 1.5.1. Die `publishDir`-Direktive aus `COWPY` entfernen

Lass uns das tun.
Öffne die Moduldatei `cowpy.nf` (unter `core-hello/modules/local/`) und entferne die `publishDir`-Direktive, wie unten gezeigt.

=== "Nachher"

    ```groovy title="core-hello/modules/local/cowpy.nf (Auszug)" linenums="1"
    #!/usr/bin/env nextflow

    // ASCII-Kunst mit cowpy generieren (https://github.com/jeffbuttars/cowpy)
    process COWPY {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
        conda 'conda-forge::cowpy==1.1.5'
    ```

=== "Vorher"

    ```groovy title="core-hello/modules/local/cowpy.nf (Auszug)" linenums="1" hl_lines="6"
    #!/usr/bin/env nextflow

    // ASCII-Kunst mit cowpy generieren (https://github.com/jeffbuttars/cowpy)
    process COWPY {

        publishDir 'results', mode: 'copy'

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
        conda 'conda-forge::cowpy==1.1.5'

    ```

Das war's!

#### 1.5.2. Die Pipeline ausführen, um sie zu testen

Lass uns sehen, was passiert, wenn wir die Pipeline jetzt ausführen.

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

Schaue dir dein aktuelles Arbeitsverzeichnis an.
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

Der verantwortliche Code befindet sich in der Datei `conf/modules.config`.
Dies ist die Standard-`publishDir`-Konfiguration, die Teil der nf-core-Vorlage ist und für alle Prozesse gilt:

```groovy
process {
    publishDir = [
        path: { "${params.outdir}/${task.process.tokenize(':')[-1].tokenize('_')[0].toLowerCase()}" },
        mode: params.publish_dir_mode,
        saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    ]
}
```

Dies mag kompliziert aussehen, also schauen wir uns jede der drei Komponenten an:

- **`path:`** Bestimmt das Ausgabeverzeichnis basierend auf dem Prozessnamen.
  Der vollständige Name eines in `task.process` enthaltenen Prozesses umfasst die Hierarchie von Workflow- und Modulimporten (wie `CORE_HELLO:HELLO:CAT_CAT`).
  Die `tokenize`-Operationen entfernen diese Hierarchie, um nur den Prozessnamen zu erhalten, nehmen dann den ersten Teil vor einem Unterstrich (falls zutreffend) und konvertieren ihn in Kleinbuchstaben.
  Dies bestimmt, dass die Ergebnisse von `CAT_CAT` nach `${params.outdir}/cat/` veröffentlicht werden.
- **`mode:`** Steuert, wie Dateien veröffentlicht werden (copy, symlink usw.).
  Dies ist über den Parameter `params.publish_dir_mode` konfigurierbar.
- **`saveAs:`** Filtert, welche Dateien veröffentlicht werden sollen.
  Dieses Beispiel schließt `versions.yml`-Dateien aus, indem es `null` für sie zurückgibt und verhindert, dass sie veröffentlicht werden.

Dies bietet eine konsistente Logik zur Organisation von Ausgaben.

Die Ausgabe sieht noch besser aus, wenn alle Module in einer Pipeline diese Konvention übernehmen, also lösche gerne die `publishDir`-Direktiven aus den anderen Modulen in deiner Pipeline.
Dieser Standard wird auch auf Module angewendet, die wir nicht explizit geändert haben, um den nf-core-Richtlinien zu folgen.

Allerdings kannst du entscheiden, dass du deine Eingaben anders organisieren möchtest, und die gute Nachricht ist, dass dies einfach zu tun ist.

#### 1.5.3. Den Standard überschreiben

Um die Standard-`publishDir`-Direktive zu überschreiben, kannst du einfach deine eigenen Direktiven zur Datei `conf/modules.config` hinzufügen.

Du könntest beispielsweise den Standard für einen einzelnen Prozess mit dem `withName:`-Selektor überschreiben, wie in diesem Beispiel, wo wir eine benutzerdefinierte `publishDir`-Direktive für den 'COWPY'-Prozess hinzufügen.

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

Wir werden diese Änderung nicht tatsächlich vornehmen, aber spiele gerne damit herum und sieh, welche Logik du implementieren kannst.

Der Punkt ist, dass dieses System dir das Beste aus beiden Welten gibt: Konsistenz standardmäßig und die Flexibilität, die Konfiguration bei Bedarf anzupassen.

Zusammenfassend erhältst du:

- **Single Source of Truth**: Alle Publishing-Konfigurationen befinden sich in `modules.config`
- **Nützlicher Standard**: Prozesse funktionieren sofort ohne modulspezifische Konfiguration
- **Einfache Anpassung**: Überschreibe das Publishing-Verhalten in der Konfiguration, nicht im Modulcode
- **Portable Module**: Module enthalten keine fest codierten Ausgabeorte

Dies vervollständigt die Reihe der nf-core-Modul-Features, die du unbedingt lernen solltest, aber es gibt weitere, über die du in den [nf-core Modul-Spezifikationen](https://nf-co.re/docs/guidelines/components/modules) nachlesen kannst.

### Zusammenfassung

Du weißt jetzt, wie du lokale Module anpassen kannst, um nf-core-Konventionen zu folgen:

- Entwirf deine Module so, dass sie Metadaten-Tupel akzeptieren und weitergeben;
- Verwende `ext.args`, um Modul-Interfaces minimal und portabel zu halten;
- Verwende `ext.prefix` für konfigurierbare, standardisierte Ausgabedatei-Benennung;
- Übernimm die standardmäßig zentralisierte `publishDir`-Direktive für eine konsistente Ergebnis-Verzeichnisstruktur.

### Was kommt als Nächstes?

Erfahre, wie du die integrierten vorlagenbasierten Tools von nf-core verwenden kannst, um Module auf einfache Weise zu erstellen.

---

## 2. Ein Modul mit den nf-core-Tools erstellen

Nachdem du die nf-core-Modulmuster durch manuelle Anwendung gelernt hast, schauen wir uns an, wie du Module in der Praxis erstellen würdest.

### 2.1. Ein Modul-Gerüst aus einer Vorlage generieren

Ähnlich wie bei der Erstellung von Pipelines bietet das nf-core-Projekt Tools, um korrekt strukturierte Module basierend auf einer Vorlage zu generieren, mit allen diesen Mustern von Anfang an integriert.

#### 2.1.1. Den Modul-Erstellungsbefehl ausführen

Der Befehl `nf-core modules create` generiert eine Modulvorlage, die bereits alle Konventionen befolgt, die du gelernt hast.

Erstellen wir eine neue Version des `COWPY`-Moduls mit einer minimalen Vorlage, indem wir diesen Befehl ausführen:

```bash
nf-core modules create --empty-template COWPY
```

Das Flag `--empty-template` erstellt eine saubere Startvorlage ohne zusätzlichen Code, wodurch die wesentliche Struktur leichter zu erkennen ist.

Der Befehl wird interaktiv ausgeführt und führt dich durch die Einrichtung.
Er sucht automatisch nach Tool-Informationen aus Paket-Repositories wie Bioconda und bio.tools, um Metadaten vorab auszufüllen.

Du wirst nach mehreren Konfigurationsoptionen gefragt:

- **Autoreninformationen**: Dein GitHub-Benutzername für die Zuordnung
- **Ressourcen-Label**: Ein vordefinierter Satz von Rechenanforderungen.
  Das nf-core-Projekt stellt Standard-Labels wie `process_single` für leichtgewichtige Tools und `process_high` für anspruchsvolle bereit.
  Diese Labels helfen bei der Verwaltung der Ressourcenzuweisung in verschiedenen Ausführungsumgebungen.
- **Metadaten-Anforderung**: Ob das Modul probenspezifische Informationen über eine `meta`-Map benötigt (normalerweise ja für Datenverarbeitungsmodule).

Das Tool übernimmt die Komplexität der Paketinformationssuche und der Struktureinrichtung, sodass du dich auf die Implementierung der spezifischen Tool-Logik konzentrieren kannst.

#### 2.1.2. Das Modul-Gerüst untersuchen

Das Tool erstellt eine vollständige Modulstruktur in `modules/local/` (oder `modules/nf-core/`, wenn du dich im nf-core/modules-Repository befindest):

??? abstract "Verzeichnisinhalt"

    ```console
    modules/local/cowpy
    ├── environment.yml
    ├── main.nf
    ├── meta.yml
    └── tests
        └── main.nf.test
    ```

Jede Datei erfüllt einen bestimmten Zweck:

- **`main.nf`**: Prozessdefinition mit allen integrierten nf-core-Mustern
- **`meta.yml`**: Moduldokumentation, die Eingaben, Ausgaben und das Tool beschreibt
- **`environment.yml`**: Conda-Umgebungsspezifikation für Abhängigkeiten
- **`tests/main.nf.test`**: nf-test-Testfälle zur Validierung der Modulfunktionalität

!!! tip "Mehr über Tests erfahren"

    Die generierte Testdatei verwendet nf-test, ein Test-Framework für Nextflow-Pipelines und -Module. Um zu lernen, wie man diese Tests schreibt und ausführt, siehe die [nf-test-Nebenquest](../side_quests/nf-test.md).

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
    tuple val(meta), path(input)        // Pattern 1: Metadata tuples ✓

    output:
    tuple val(meta), path("*"), emit: output
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''                // Pattern 2: ext.args ✓
    def prefix = task.ext.prefix ?: "${meta.id}"  // Pattern 3: ext.prefix ✓

    """
    // Add your tool command here

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

Beachte, dass alle Muster, die du oben manuell angewendet hast, bereits vorhanden sind!

Die Vorlage enthält auch mehrere zusätzliche nf-core-Konventionen.
Einige davon funktionieren sofort, während andere Platzhalter sind, die wir ausfüllen müssen, wie unten beschrieben.

**Features, die sofort funktionieren:**

- **`tag "$meta.id"`**: Fügt die Proben-ID zu Prozessnamen in den Logs hinzu, um die Nachverfolgung zu erleichtern
- **`label 'process_single'`**: Ressourcen-Label zur Konfiguration von CPU/Speicher-Anforderungen
- **`when:`-Block**: Ermöglicht bedingte Ausführung über die `task.ext.when`-Konfiguration

Diese Features sind bereits funktionsfähig und machen Module wartbarer.

**Platzhalter, die wir unten anpassen:**

- **`input:`- und `output:`-Blöcke**: Generische Deklarationen, die wir an unser Tool anpassen
- **`script:`-Block**: Enthält einen Kommentar, wo wir den `cowpy`-Befehl einfügen
- **`stub:`-Block**: Vorlage, die wir aktualisieren, um die korrekten Ausgaben zu erzeugen
- **Container und Umgebung**: Platzhalter, die wir mit Paketinformationen ausfüllen

Die nächsten Abschnitte führen durch die Vervollständigung dieser Anpassungen.

### 2.2. Container und Conda-Umgebung einrichten

Die nf-core-Richtlinien erfordern, dass wir sowohl einen Container als auch eine Conda-Umgebung als Teil des Moduls angeben.

#### 2.2.1. Container

Für den Container kannst du [Seqera Containers](https://seqera.io/containers/) verwenden, um automatisch einen Container aus jedem Conda-Paket zu erstellen, einschließlich conda-forge-Paketen.
In diesem Fall verwenden wir denselben vorgefertigten Container wie zuvor.

Der Standardcode bietet einen Wechsel zwischen Docker und Singularity an, aber wir werden diese Zeile vereinfachen und einfach den Docker-Container angeben, den wir von Seqera Containers oben erhalten haben.

=== "Nachher"

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

Für die Conda-Umgebung gibt der Modulcode `conda "${moduleDir}/environment.yml"` an, was bedeutet, dass sie in der Datei `environment.yml` konfiguriert werden sollte.

Das Modul-Erstellungstool hat uns gewarnt, dass es das `cowpy`-Paket in Bioconda (dem primären Kanal für Bioinformatik-Tools) nicht finden konnte.
Allerdings ist `cowpy` in conda-forge verfügbar, sodass du die `environment.yml` wie folgt vervollständigen kannst:

=== "Nachher"

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
      # TODO nf-core: List required Conda package(s).
      #               Software MUST be pinned to channel (i.e. "bioconda"), version (i.e. "1.10").
      #               For Conda, the build (i.e. "h9402c20_2") must be EXCLUDED to support installation on different operating systems.
      - "YOUR-TOOL-HERE"
    ```

Für eine Einreichung bei nf-core müssten wir die Standards genauer einhalten, aber für unseren eigenen Gebrauch können wir den Code auf diese Weise vereinfachen.

!!! tip "Bioconda vs conda-forge Pakete"

    - **Bioconda-Pakete**: Erhalten automatisch BioContainers, die sofort einsatzbereite Container bereitstellen
    - **conda-forge-Pakete**: Können Seqera Containers verwenden, um Container on-demand aus dem Conda-Rezept zu erstellen

    Die meisten Bioinformatik-Tools befinden sich in Bioconda, aber für conda-forge-Tools bietet Seqera Containers eine einfache Lösung für die Containerisierung.

### 2.3. Die `COWPY`-Logik einbauen

Aktualisieren wir nun die Code-Elemente, die spezifisch für das sind, was der `COWPY`-Prozess tut: die Eingaben und Ausgaben sowie den Script-Block.

#### 2.3.1. Eingaben und Ausgaben

Die generierte Vorlage enthält generische Eingabe- und Ausgabedeklarationen, die du für dein spezifisches Tool anpassen musst.
Wenn wir auf unser manuelles `COWPY`-Modul aus Abschnitt 1 zurückblicken, können wir das als Leitfaden verwenden.

Aktualisiere die Eingabe- und Ausgabeblöcke:

=== "Nachher"

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

- Den Eingabedatei-Parameternamen (`input_file` statt dem generischen `input`)
- Den Ausgabedateinamen mit dem konfigurierbaren Präfix-Muster (`${prefix}.txt` statt Wildcard `*`)
- Einen beschreibenden Emit-Namen (`cowpy_output` statt dem generischen `output`)

Wenn du den Nextflow Language Server zur Syntaxvalidierung verwendest, wird der Teil `${prefix}` in dieser Phase als Fehler markiert, weil wir ihn noch nicht zum Script-Block hinzugefügt haben.
Kommen wir jetzt dazu.

#### 2.3.2. Der Script-Block

Die Vorlage bietet einen Kommentar-Platzhalter im Script-Block, wo du den eigentlichen Tool-Befehl einfügen solltest.

Basierend auf dem Modul, das wir zuvor manuell geschrieben haben, sollten wir die folgenden Änderungen vornehmen:

=== "Nachher"

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
    // Add your tool command here

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        COWPY: \$(cowpy --version)
    END_VERSIONS
    """
    ```

Wichtige Änderungen:

- Ändere `def prefix` zu nur `prefix` (ohne `def`), um es im Output-Block zugänglich zu machen
- Ersetze den Kommentar durch den eigentlichen `cowpy`-Befehl, der sowohl `$args` als auch `${prefix}.txt` verwendet

Beachte, dass wir, wenn wir nicht bereits die Arbeit des Hinzufügens der `ext.args`- und `ext.prefix`-Konfiguration für den `COWPY`-Prozess zur Datei `modules.config` erledigt hätten, dies jetzt tun müssten.

#### 2.3.3. Den Stub-Block implementieren

Im Nextflow-Kontext erlaubt ein [Stub](https://www.nextflow.io/docs/latest/process.html#stub)-Block, ein leichtgewichtiges Dummy-Skript zu definieren, das für schnelles Prototyping und Testen der Pipeline-Logik verwendet wird, ohne den eigentlichen Befehl auszuführen.

<!-- TODO (future) This is super glossed over but should really be explained or at least link out to an explanation about stubs (the reference doc isn't terribly helpful either). Right now this is likely to be mostly meaningless to anyone who doesn't already know about stubs. -->

Mach dir nicht zu viele Sorgen, wenn dies mysteriös erscheint; wir fügen dies der Vollständigkeit halber hinzu, aber du kannst den Stub-Abschnitt auch einfach löschen, wenn du dich nicht damit befassen möchtest, da er komplett optional ist.

=== "Nachher"

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

- Ändere `def prefix` zu nur `prefix`, um es dem Script-Block anzupassen
- Entferne die Zeile `echo $args` (die nur Vorlagen-Platzhaltercode war)
- Der Stub erstellt eine leere `${prefix}.txt`-Datei, die dem entspricht, was der Script-Block erzeugt

Dies ermöglicht es dir, die Workflow-Logik und Dateibehandlung zu testen, ohne auf die Ausführung des eigentlichen Tools warten zu müssen.

Sobald du die Umgebungseinrichtung (Abschnitt 2.2), Eingaben/Ausgaben (Abschnitt 2.3.1), den Script-Block (Abschnitt 2.3.2) und den Stub-Block (Abschnitt 2.3.3) abgeschlossen hast, ist das Modul bereit zum Testen!

### 2.4. Das neue `COWPY`-Modul einsetzen und die Pipeline ausführen

Alles, was wir tun müssen, um diese neue Version des `COWPY`-Moduls auszuprobieren, ist die Import-Anweisung in der Workflow-Datei `hello.nf` auf die neue Datei umzuleiten.

=== "Nachher"

    ```groovy title="workflows/hello.nf" linenums="1" hl_lines="10"
    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
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
        IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    include { paramsSummaryMap       } from 'plugin/nf-schema'
    include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
    include { sayHello               } from '../modules/local/sayHello.nf'
    include { convertToUpper         } from '../modules/local/convertToUpper.nf'
    include { COWPY                  } from '../modules/local/cowpy.nf'
    include { CAT_CAT                } from '../modules/nf-core/cat/cat/main'
    ```

Führen wir die Pipeline aus, um sie zu testen.

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

Dies erzeugt dieselben Ergebnisse wie zuvor.

### Zusammenfassung

Du weißt jetzt, wie du die integrierten nf-core-Tools verwenden kannst, um Module effizient mit Vorlagen zu erstellen, anstatt alles von Grund auf zu schreiben.

### Was kommt als Nächstes?

Erfahre, welche Vorteile es hat, Module zu nf-core beizutragen, und welche Hauptschritte und Anforderungen damit verbunden sind.

---

## 3. Module zu nf-core beitragen

Das [nf-core/modules](https://github.com/nf-core/modules)-Repository begrüßt Beiträge gut getesteter, standardisierter Module.

### 3.1. Warum beitragen?

Das Beitragen deiner Module zu nf-core:

- Macht deine Tools der gesamten nf-core-Community über den Modulkatalog unter [nf-co.re/modules](https://nf-co.re/modules) zugänglich
- Gewährleistet fortlaufende Community-Wartung und Verbesserungen
- Bietet Qualitätssicherung durch Code-Review und automatisierte Tests
- Gibt deiner Arbeit Sichtbarkeit und Anerkennung

### 3.2. Checkliste für Beitragende

Um ein Modul zu nf-core beizutragen, musst du die folgenden Schritte durchlaufen:

1. Prüfe, ob es bereits unter [nf-co.re/modules](https://nf-co.re/modules) existiert
2. Forke das [nf-core/modules](https://github.com/nf-core/modules)-Repository
3. Verwende `nf-core modules create`, um die Vorlage zu generieren
4. Fülle die Modullogik und Tests aus
5. Teste mit `nf-core modules test tool/subtool`
6. Linte mit `nf-core modules lint tool/subtool`
7. Reiche einen Pull Request ein

Für detaillierte Anweisungen siehe das [nf-core Components Tutorial](https://nf-co.re/docs/tutorials/nf-core_components/components).

### 3.3. Ressourcen

- **Components Tutorial**: [Vollständige Anleitung zum Erstellen und Beitragen von Modulen](https://nf-co.re/docs/tutorials/nf-core_components/components)
- **Modul-Spezifikationen**: [Technische Anforderungen und Richtlinien](https://nf-co.re/docs/guidelines/components/modules)
- **Community-Support**: [nf-core Slack](https://nf-co.re/join) - Tritt dem Kanal `#modules` bei

### Zusammenfassung

Du weißt jetzt, wie man nf-core-Module erstellt! Du hast die vier Schlüsselmuster gelernt, die Module portabel und wartbar machen:

- **Metadaten-Tupel** propagieren Metadaten durch den Workflow
- **`ext.args`** vereinfacht Modul-Interfaces durch die Handhabung optionaler Argumente über Konfiguration
- **`ext.prefix`** standardisiert die Benennung von Ausgabedateien
- **Zentralisiertes Publishing** über `publishDir`, konfiguriert in `modules.config` statt fest codiert in Modulen

Durch die schrittweise Transformation von `COWPY` hast du ein tiefes Verständnis dieser Muster entwickelt, das dich befähigt, mit nf-core-Modulen zu arbeiten, sie zu debuggen und zu erstellen.
In der Praxis wirst du `nf-core modules create` verwenden, um korrekt strukturierte Module mit diesen Mustern von Anfang an zu generieren.

Schließlich hast du gelernt, wie du Module zur nf-core-Community beitragen kannst, um Tools Forschenden weltweit zur Verfügung zu stellen und gleichzeitig von der fortlaufenden Community-Wartung zu profitieren.

### Was kommt als Nächstes?

Wenn du bereit bist, fahre mit [Teil 5: Eingabevalidierung](./05_input_validation.md) fort, um zu lernen, wie du Schema-basierte Eingabevalidierung zu deiner Pipeline hinzufügst.
