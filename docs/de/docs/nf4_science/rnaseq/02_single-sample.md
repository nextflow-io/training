# Teil 2: Implementierung für eine einzelne Probe

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

In diesem Teil des Kurses werden wir einen Workflow entwickeln, der alle Befehle aus Teil 1 zusammenfasst, um ihre Ausführung zu automatisieren. Dabei werden wir zunächst nur eine Probe gleichzeitig verarbeiten.

!!! warning "Voraussetzung"

    Du musst [Teil 1: Methodenübersicht](./01_method.md) durcharbeiten, bevor du mit dieser Lektion beginnst.
    Insbesondere durch das Durcharbeiten von Abschnitt 1.2.3 wird die Genom-Indexdatei (`data/genome_index.tar.gz`) erstellt, die für den Alignment-Schritt in dieser Lektion benötigt wird.

## Aufgabe

In diesem Teil des Kurses werden wir einen Workflow entwickeln, der Folgendes tut:

1. Qualitätskontrolle (FastQC) auf Eingabe-Reads ausführen
2. Adapter trimmen und Post-Trimming-QC ausführen (Trim Galore)
3. Getrimmte Reads zum Referenzgenom alignieren (HISAT2)

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/rnaseq/img/rnaseq-wf-02.svg"
</figure>

Dies automatisiert die Schritte aus dem ersten Abschnitt von [Teil 1: Methodenübersicht](./01_method.md#1-single-sample-processing), wo du diese Befehle manuell in ihren Containern ausgeführt hast.

Als Ausgangspunkt stellen wir dir eine Workflow-Datei, `rnaseq.nf`, zur Verfügung, die die Hauptteile des Workflows skizziert, sowie vier Moduldateien im Verzeichnis `modules/` (`fastqc.nf`, `trim_galore.nf`, `hisat2_align.nf` und `multiqc.nf`), die die Struktur jedes Prozesses skizzieren.

??? full-code "Gerüstdateien"

    ```groovy title="rnaseq.nf"
    #!/usr/bin/env nextflow

    // Modul-INCLUDE-Anweisungen

    /*
     * Pipeline parameters
     */

    // Primäre Eingabe

    workflow {

        main:
        // Eingabe-Channel erstellen

        // Prozesse aufrufen

        publish:
        // Zu veröffentlichende Ausgaben deklarieren
    }

    output {
        // Veröffentlichungsziele konfigurieren
    }
    ```

    ```groovy title="modules/fastqc.nf"
    #!/usr/bin/env nextflow

    /*
     * Run FastQC on input reads
     */
    process FASTQC {

        container

        input:

        output:

        script:
        """

        """
    }
    ```

    ```groovy title="modules/trim_galore.nf"
    #!/usr/bin/env nextflow

    /*
     * Trim adapters and run post-trimming QC
     */
    process TRIM_GALORE {

        container

        input:

        output:

        script:
        """

        """
    }
    ```

    ```groovy title="modules/hisat2_align.nf"
    #!/usr/bin/env nextflow

    /*
     * Align reads to a reference genome
     */
    process HISAT2_ALIGN {

        container

        input:

        output:

        script:
        """

        """
    }
    ```

    ```groovy title="modules/multiqc.nf"
    #!/usr/bin/env nextflow

    /*
     * Aggregate QC reports with MultiQC
     */
    process MULTIQC {

        container

        input:

        output:

        script:
        """

        """
    }
    ```

Diese Dateien sind nicht funktionsfähig; ihr Zweck ist es lediglich, als Gerüst zu dienen, das du mit den interessanten Teilen des Codes ausfüllen kannst.

## Lektionsplan

Um den Entwicklungsprozess lehrreicher zu gestalten, haben wir dies in drei Phasen unterteilt:

1. **Einen einstufigen Workflow schreiben, der den ersten QC-Schritt ausführt.**
   Dies umfasst das Einrichten eines CLI-Parameters, das Erstellen eines Eingabekanals, das Schreiben eines Prozessmoduls und das Konfigurieren der Ausgabeveröffentlichung.
2. **Adapter-Trimming und Post-Trimming-QC hinzufügen.**
   Dies führt das Verketten von Prozessen ein, indem die Ausgabe eines Prozesses mit der Eingabe eines anderen verbunden wird.
3. **Alignment zum Referenzgenom hinzufügen.**
   Dies behandelt die Handhabung zusätzlicher Referenzeingaben und die Arbeit mit komprimierten Archiven.

Jeder Schritt konzentriert sich auf einen bestimmten Aspekt der Workflow-Entwicklung.

!!! tip "Tipp"

     Stelle sicher, dass du im richtigen Arbeitsverzeichnis bist:
     `cd /workspaces/training/nf4-science/rnaseq`

---

## 1. Einen einstufigen Workflow schreiben, der die erste QC ausführt

Dieser erste Schritt konzentriert sich auf die Grundlagen: Laden einer FASTQ-Datei und Ausführen der Qualitätskontrolle darauf.

Erinnere dich an den `fastqc`-Befehl aus [Teil 1](01_method.md):

```bash
fastqc <reads>
```

Der Befehl nimmt eine FASTQ-Datei als Eingabe und erzeugt einen Qualitätskontrollbericht als `.zip`-Archiv und eine `.html`-Zusammenfassung.
Die Container-URI war `community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18`.

Wir werden diese Informationen nehmen und in Nextflow in drei Phasen einbinden:

1. Die Eingabe einrichten
2. Den QC-Prozess schreiben und im Workflow aufrufen
3. Die Ausgabebehandlung konfigurieren

### 1.1. Die Eingabe einrichten

Wir müssen einen Eingabeparameter deklarieren, ein Testprofil erstellen, um einen praktischen Standardwert bereitzustellen, und einen Eingabekanal erstellen.

#### 1.1.1. Eine Eingabeparameterdeklaration hinzufügen

Deklariere in `rnaseq.nf` unter dem Abschnitt `Pipeline parameters` einen Parameter namens `reads` mit dem Typ `Path`.

=== "Danach"

    ```groovy title="rnaseq.nf" linenums="5" hl_lines="4-7"
    /*
     * Pipeline parameters
     */
    params {
        // Primäre Eingabe
        input: Path
    }
    ```

=== "Vorher"

    ```groovy title="rnaseq.nf" linenums="5"
    /*
     * Pipeline parameters
     */

    // Primäre Eingabe
    ```

Das richtet den CLI-Parameter ein, aber wir wollen nicht jedes Mal den Dateipfad eingeben, wenn wir den Workflow während der Entwicklung ausführen.
Es gibt mehrere Optionen, um einen Standardwert bereitzustellen; hier verwenden wir ein Testprofil.

#### 1.1.2. Ein Testprofil mit einem Standardwert in `nextflow.config` erstellen

Ein Testprofil bietet praktische Standardwerte zum Ausprobieren eines Workflows, ohne Eingaben auf der Befehlszeile anzugeben.
Dies ist eine gängige Konvention im Nextflow-Ökosystem (siehe [Hello Config](../../hello_nextflow/06_hello_config.md) für weitere Details).

Füge einen `profiles`-Block zu `nextflow.config` mit einem `test`-Profil hinzu, das den `reads`-Parameter auf eine der Test-FASTQ-Dateien setzt.

=== "Danach"

    ```groovy title="nextflow.config" linenums="1" hl_lines="3-7"
    docker.enabled = true

    profiles {
        test {
            params.input = "${projectDir}/data/reads/ENCSR000COQ1_1.fastq.gz"
        }
    }
    ```

=== "Vorher"

    ```groovy title="nextflow.config" linenums="1"
    docker.enabled = true
    ```

Hier verwenden wir `#!groovy ${projectDir}`, eine eingebaute Nextflow-Variable, die auf das Verzeichnis zeigt, in dem sich das Workflow-Skript befindet.
Dies macht es einfach, auf Datendateien und andere Ressourcen zu verweisen, ohne absolute Pfade fest zu codieren.

Der Parameter hat jetzt einen praktischen Standardwert. Als Nächstes müssen wir einen Kanal daraus erstellen.

#### 1.1.3. Den Eingabekanal einrichten

Erstelle im Workflow-Block einen Eingabekanal aus dem Parameterwert mit der `.fromPath`-Channel-Factory (wie in [Hello Channels](../../hello_nextflow/02_hello_channels.md) verwendet).

=== "Danach"

    ```groovy title="rnaseq.nf" linenums="13" hl_lines="4-5"
    workflow {

        main:
        // Eingabe-Channel aus einem Dateipfad erstellen
        read_ch = channel.fromPath(params.input)

        // Prozesse aufrufen

        publish:
        // Zu veröffentlichende Ausgaben deklarieren
    }
    ```

=== "Vorher"

    ```groovy title="rnaseq.nf" linenums="13"
    workflow {

        main:
        // Eingabe-Channel erstellen

        // Prozesse aufrufen

        publish:
        // Zu veröffentlichende Ausgaben deklarieren
    }
    ```

Als Nächstes müssen wir den Prozess erstellen, um QC auf dieser Eingabe auszuführen.

### 1.2. Den QC-Prozess schreiben und im Workflow aufrufen

Wir müssen die Prozessdefinition in der Moduldatei ausfüllen, sie mit einer Include-Anweisung in den Workflow importieren und sie auf der Eingabe aufrufen.

#### 1.2.1. Das Modul für den QC-Prozess ausfüllen

Öffne `modules/fastqc.nf` und untersuche die Gliederung der Prozessdefinition.
Du solltest die Hauptstrukturelemente wiedererkennen; falls nicht, lies [Hello Nextflow](../../hello_nextflow/01_hello_world.md) zur Auffrischung.

Fülle die Prozessdefinition selbst aus, indem du die oben bereitgestellten Informationen verwendest, und überprüfe dann deine Arbeit anhand der Lösung im Tab "Danach" unten.

=== "Vorher"

    ```groovy title="modules/fastqc.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
     * Run FastQC on input reads
     */
    process FASTQC {

        container

        input:

        output:

        script:
        """

        """
    }
    ```

=== "Danach"

    ```groovy title="modules/fastqc.nf" linenums="1" hl_lines="8 11 14 15 19"
    #!/usr/bin/env nextflow

    /*
     * Run FastQC on input reads
     */
    process FASTQC {

        container "community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18"

        input:
        path reads

        output:
        path "${reads.simpleName}_fastqc.zip", emit: zip
        path "${reads.simpleName}_fastqc.html", emit: html

        script:
        """
        fastqc ${reads}
        """
    }
    ```

Der `simpleName`-Accessor entfernt alle Erweiterungen vom Dateinamen, sodass `ENCSR000COQ1_1.fastq.gz` zu `ENCSR000COQ1_1` wird.
Wir verwenden die `emit:`-Syntax, um jedem Ausgabekanal Namen zuzuweisen, was nützlich sein wird, um Ausgaben in den Publish-Block zu leiten.

Sobald du dies abgeschlossen hast, ist der Prozess fertig.
Um ihn im Workflow zu verwenden, musst du das Modul importieren und einen Prozessaufruf hinzufügen.

#### 1.2.2. Das Modul einbinden

Füge in `rnaseq.nf` eine `include`-Anweisung hinzu, um den Prozess für den Workflow verfügbar zu machen:

=== "Danach"

    ```groovy title="rnaseq.nf" linenums="3" hl_lines="2"
    // Modul-INCLUDE-Anweisungen
    include { FASTQC } from './modules/fastqc.nf'
    ```

=== "Vorher"

    ```groovy title="rnaseq.nf" linenums="3"
    // Modul-INCLUDE-Anweisungen
    ```

Der Prozess ist jetzt im Workflow-Scope verfügbar.

#### 1.2.3. Den QC-Prozess auf der Eingabe aufrufen

Füge einen Aufruf von `FASTQC` im Workflow-Block hinzu und übergebe den Eingabekanal als Argument.

=== "Danach"

    ```groovy title="rnaseq.nf" linenums="14" hl_lines="7-8"
    workflow {

        main:
        // Eingabe-Channel aus einem Dateipfad erstellen
        read_ch = channel.fromPath(params.input)

        // Initiale Qualitätskontrolle
        FASTQC(read_ch)

        publish:
        // Zu veröffentlichende Ausgaben deklarieren
    }
    ```

=== "Vorher"

    ```groovy title="rnaseq.nf" linenums="14"
    workflow {

        main:
        // Eingabe-Channel aus einem Dateipfad erstellen
        read_ch = channel.fromPath(params.input)

        // Prozesse aufrufen

        publish:
        // Zu veröffentlichende Ausgaben deklarieren
    }
    ```

Der Workflow lädt jetzt die Eingabe und führt den QC-Prozess darauf aus.
Als Nächstes müssen wir konfigurieren, wie die Ausgabe veröffentlicht wird.

### 1.3. Die Ausgabebehandlung konfigurieren

Wir müssen deklarieren, welche Prozessausgaben veröffentlicht werden sollen, und angeben, wohin sie gehen sollen.

#### 1.3.1. Ausgaben im `publish:`-Abschnitt deklarieren

Der `publish:`-Abschnitt innerhalb des Workflow-Blocks deklariert, welche Prozessausgaben veröffentlicht werden sollen.
Weise die Ausgaben von `FASTQC` benannten Zielen zu.

=== "Danach"

    ```groovy title="rnaseq.nf" linenums="23" hl_lines="2-3"
        publish:
        fastqc_zip = FASTQC.out.zip
        fastqc_html = FASTQC.out.html
    }
    ```

=== "Vorher"

    ```groovy title="rnaseq.nf" linenums="23"
        publish:
        // Zu veröffentlichende Ausgaben deklarieren
    }
    ```

Als Nächstes müssen wir Nextflow mitteilen, wohin die veröffentlichten Ausgaben gelegt werden sollen.

#### 1.3.2. Die Ausgabeziele im `output {}`-Block konfigurieren

Der `output {}`-Block befindet sich außerhalb des Workflows und gibt an, wohin jedes benannte Ziel veröffentlicht wird.
Konfiguriere beide Ziele so, dass sie in ein `fastqc/`-Unterverzeichnis veröffentlicht werden.

=== "Danach"

    ```groovy title="rnaseq.nf" linenums="28" hl_lines="2-7"
    output {
        fastqc_zip {
            path 'fastqc'
        }
        fastqc_html {
            path 'fastqc'
        }
    }
    ```

=== "Vorher"

    ```groovy title="rnaseq.nf" linenums="28"
    output {
        // Veröffentlichungsziele konfigurieren
    }
    ```

!!! note "Hinweis"

    Standardmäßig veröffentlicht Nextflow Ausgabedateien als symbolische Links, was unnötige Duplizierung vermeidet.
    Obwohl die Datendateien, die wir hier verwenden, sehr klein sind, können sie in der Genomik sehr groß werden.
    Symlinks funktionieren nicht mehr, wenn du dein `work`-Verzeichnis aufräumst. Für Produktions-Workflows möchtest du daher möglicherweise den Standard-Veröffentlichungsmodus auf `'copy'` überschreiben.

### 1.4. Den Workflow ausführen

An diesem Punkt haben wir einen einstufigen QC-Workflow, der voll funktionsfähig sein sollte.

Wir führen mit `-profile test` aus, um den im Testprofil eingerichteten Standardwert zu verwenden und zu vermeiden, den Pfad auf der Befehlszeile schreiben zu müssen.

```bash
nextflow run rnaseq.nf -profile test
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `rnaseq.nf` [mad_lorenz] DSL2 - revision: 5846a164d2

    executor >  local (1)
    [7b/8ee79e] FASTQC (1) | 1 of 1 ✔
    ```

Dies sollte sehr schnell ausgeführt werden, wenn du Teil 1 durchgearbeitet hast und den Container bereits heruntergeladen hast.
Wenn du ihn übersprungen hast, wird Nextflow den Container für dich herunterladen; du musst nichts dafür tun, aber du musst möglicherweise bis zu einer Minute warten.

Du kannst die Ausgaben im Ergebnisverzeichnis überprüfen.

```bash
ls results/fastqc
```

```console title="Output"
ENCSR000COQ1_1_fastqc.html  ENCSR000COQ1_1_fastqc.zip
```

Die QC-Berichte für die Probe sind jetzt im `fastqc/`-Unterverzeichnis veröffentlicht.

### Fazit

Du weißt jetzt, wie du ein Modul mit einem Prozess erstellst, es in einen Workflow importierst, es mit einem Eingabekanal aufrufst und die Ergebnisse mit dem Workflow-Level-Output-Block veröffentlichst.

### Wie geht es weiter?

Füge Adapter-Trimming mit Post-Trimming-QC als zweiten Schritt im Workflow hinzu.

---

## 2. Adapter-Trimming und Post-Trimming-Qualitätskontrolle hinzufügen

Jetzt, da wir die initiale QC eingerichtet haben, können wir den Adapter-Trimming-Schritt mit seiner eingebauten Post-Trimming-QC hinzufügen.

Erinnere dich an den `trim_galore`-Befehl aus [Teil 1](01_method.md):

```bash
trim_galore --fastqc <reads>
```

Der Befehl trimmt Adapter aus einer FASTQ-Datei und führt FastQC auf der getrimmten Ausgabe aus.
Er erzeugt getrimmte Reads, einen Trimming-Bericht und FastQC-Berichte für die getrimmten Reads.
Die Container-URI war `community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18`.

Wir müssen nur die Prozessdefinition schreiben, sie importieren, sie im Workflow aufrufen und die Ausgabebehandlung aktualisieren.

### 2.1. Den Trimming-Prozess schreiben und im Workflow aufrufen

Wie zuvor müssen wir die Prozessdefinition ausfüllen, das Modul importieren und den Prozessaufruf hinzufügen.

#### 2.1.1. Das Modul für den Trimming-Prozess ausfüllen

Öffne `modules/trim_galore.nf` und untersuche die Gliederung der Prozessdefinition.

Fülle die Prozessdefinition selbst aus, indem du die oben bereitgestellten Informationen verwendest, und überprüfe dann deine Arbeit anhand der Lösung im Tab "Danach" unten.

=== "Vorher"

    ```groovy title="modules/trim_galore.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
     * Trim adapters and run post-trimming QC
     */
    process TRIM_GALORE {

        container

        input:

        output:

        script:
        """

        """
    }
    ```

=== "Danach"

    ```groovy title="modules/trim_galore.nf" linenums="1" hl_lines="8 11 14 15 16 20"
    #!/usr/bin/env nextflow

    /*
     * Trim adapters and run post-trimming QC
     */
    process TRIM_GALORE {

        container "community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18"

        input:
        path reads

        output:
        path "${reads.simpleName}_trimmed.fq.gz", emit: trimmed_reads
        path "${reads}_trimming_report.txt", emit: trimming_reports
        path "${reads.simpleName}_trimmed_fastqc.{zip,html}", emit: fastqc_reports

        script:
        """
        trim_galore --fastqc ${reads}
        """
    }
    ```

Dieser Prozess hat drei benannte Ausgaben: die getrimmten Reads, die in den Alignment-Schritt einfließen, den Trimming-Bericht und die Post-Trimming-FastQC-Berichte.
Das `--fastqc`-Flag weist Trim Galore an, automatisch FastQC auf der getrimmten Ausgabe auszuführen.

#### 2.1.2. Das Modul einbinden

Aktualisiere `rnaseq.nf`, um das neue Modul zu importieren:

=== "Danach"

    ```groovy title="rnaseq.nf" linenums="3" hl_lines="3"
    // Modul-INCLUDE-Anweisungen
    include { FASTQC } from './modules/fastqc.nf'
    include { TRIM_GALORE } from './modules/trim_galore.nf'
    ```

=== "Vorher"

    ```groovy title="rnaseq.nf" linenums="3"
    // Modul-INCLUDE-Anweisungen
    include { FASTQC } from './modules/fastqc.nf'
    ```

Als Nächstes fügen wir den Prozessaufruf zum Workflow hinzu.

#### 2.1.3. Den Trimming-Prozess auf der Eingabe aufrufen

Füge den Prozessaufruf im Workflow-Block hinzu:

=== "Danach"

    ```groovy title="rnaseq.nf" linenums="15" hl_lines="10-11"
    workflow {

        main:
        // Eingabe-Channel aus einem Dateipfad erstellen
        read_ch = channel.fromPath(params.input)

        // Initiale Qualitätskontrolle
        FASTQC(read_ch)

        // Adapter-Trimming und Post-Trimming-QC
        TRIM_GALORE(read_ch)

        publish:
        fastqc_zip = FASTQC.out.zip
        fastqc_html = FASTQC.out.html
    }
    ```

=== "Vorher"

    ```groovy title="rnaseq.nf" linenums="15"
    workflow {

        main:
        // Eingabe-Channel aus einem Dateipfad erstellen
        read_ch = channel.fromPath(params.input)

        // Initiale Qualitätskontrolle
        FASTQC(read_ch)

        publish:
        fastqc_zip = FASTQC.out.zip
        fastqc_html = FASTQC.out.html
    }
    ```

Der Trimming-Prozess ist jetzt in den Workflow eingebunden.

### 2.2. Die Ausgabebehandlung aktualisieren

Wir müssen die Trimming-Ausgaben zur Publish-Deklaration hinzufügen und konfigurieren, wohin sie gehen.

#### 2.2.1. Publish-Ziele für die Trimming-Ausgaben hinzufügen

Füge die Trimming-Ausgaben zum `publish:`-Abschnitt hinzu:

=== "Danach"

    ```groovy title="rnaseq.nf" linenums="27" hl_lines="4-6"
        publish:
        fastqc_zip = FASTQC.out.zip
        fastqc_html = FASTQC.out.html
        trimmed_reads = TRIM_GALORE.out.trimmed_reads
        trimming_reports = TRIM_GALORE.out.trimming_reports
        trimming_fastqc = TRIM_GALORE.out.fastqc_reports
    }
    ```

=== "Vorher"

    ```groovy title="rnaseq.nf" linenums="27"
        publish:
        fastqc_zip = FASTQC.out.zip
        fastqc_html = FASTQC.out.html
    }
    ```

Als Nächstes müssen wir Nextflow mitteilen, wohin diese Ausgaben gelegt werden sollen.

#### 2.2.2. Die neuen Ausgabeziele konfigurieren

Füge Einträge für die Trimming-Ziele im `output {}`-Block hinzu und veröffentliche sie in ein `trimming/`-Unterverzeichnis:

=== "Danach"

    ```groovy title="rnaseq.nf" linenums="35" hl_lines="8-16"
    output {
        fastqc_zip {
            path 'fastqc'
        }
        fastqc_html {
            path 'fastqc'
        }
        trimmed_reads {
            path 'trimming'
        }
        trimming_reports {
            path 'trimming'
        }
        trimming_fastqc {
            path 'trimming'
        }
    }
    ```

=== "Vorher"

    ```groovy title="rnaseq.nf" linenums="35"
    output {
        fastqc_zip {
            path 'fastqc'
        }
        fastqc_html {
            path 'fastqc'
        }
    }
    ```

Die Ausgabekonfiguration ist abgeschlossen.

### 2.3. Den Workflow ausführen

Der Workflow umfasst jetzt sowohl initiale QC als auch Adapter-Trimming.

```bash
nextflow run rnaseq.nf -profile test
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `rnaseq.nf` [gloomy_becquerel] DSL2 - revision: bb11055736

    executor >  local (2)
    [f6/c8ef2e] FASTQC (1)      | 1 of 1 ✔
    [58/c58d8a] TRIM_GALORE (1) | 1 of 1 ✔
    ```

Dies sollte ebenfalls sehr schnell ausgeführt werden, da wir mit einer so kleinen Eingabedatei arbeiten.

Du findest die Trimming-Ausgaben im Ergebnisverzeichnis.

```bash
ls results/trimming
```

```console title="Output"
ENCSR000COQ1_1.fastq.gz_trimming_report.txt  ENCSR000COQ1_1_trimmed_fastqc.zip
ENCSR000COQ1_1_trimmed_fastqc.html           ENCSR000COQ1_1_trimmed.fq.gz
```

Die Trimming-Ausgaben und Post-Trimming-QC-Berichte befinden sich jetzt im `trimming/`-Unterverzeichnis.

### Fazit

Du weißt jetzt, wie du einen zweiten Verarbeitungsschritt hinzufügst, der unabhängig auf derselben Eingabe läuft und mehrere benannte Ausgaben erzeugt.

### Wie geht es weiter?

Füge den Alignment-Schritt hinzu, der sich an die getrimmten Reads-Ausgabe anschließt.

---

## 3. Alignment zum Referenzgenom hinzufügen

Schließlich können wir den Genom-Alignment-Schritt mit HISAT2 hinzufügen.

Erinnere dich an den Alignment-Befehl aus [Teil 1](01_method.md):

```bash
hisat2 -x <genome_index> -U <reads> \
    --new-summary --summary-file <reads>.hisat2.log | \
    samtools view -bS -o <reads>.bam
```

Der Befehl aligniert Reads zu einem Referenzgenom und konvertiert die Ausgabe in das BAM-Format.
Er benötigt ein vorgefertigtes Genom-Index-Archiv und erzeugt eine BAM-Datei und ein Alignment-Zusammenfassungslog.
Die Container-URI war `community.wave.seqera.io/library/hisat2_samtools:5e49f68a37dc010e`.

Dieser Prozess benötigt eine zusätzliche Eingabe (das Genom-Index-Archiv), also müssen wir das zuerst einrichten und dann den Prozess schreiben und verdrahten.

### 3.1. Die Eingaben einrichten

Wir müssen einen Parameter für das Genom-Index-Archiv deklarieren.

#### 3.1.1. Einen Parameter für den Genomindex hinzufügen

Füge eine Parameterdeklaration für das Genom-Index-Archiv in `rnaseq.nf` hinzu:

=== "Danach"

    ```groovy title="rnaseq.nf" linenums="11" hl_lines="5-6"
    params {
        // Primäre Eingabe
        input: Path

        // Referenzgenom-Archiv
        hisat2_index_zip: Path
    }
    ```

=== "Vorher"

    ```groovy title="rnaseq.nf" linenums="11"
    params {
        // Primäre Eingabe
        input: Path
    }
    ```

#### 3.1.2. Den Genomindex-Standard zum Testprofil hinzufügen

Genau wie wir es für `reads` in Abschnitt 1.1.2 getan haben, füge einen Standardwert für den Genomindex zum Testprofil in `nextflow.config` hinzu:

=== "Danach"

    ```groovy title="nextflow.config" linenums="1" hl_lines="6"
    docker.enabled = true

    profiles {
        test {
            params.input = "${projectDir}/data/reads/ENCSR000COQ1_1.fastq.gz"
            params.hisat2_index_zip = "${projectDir}/data/genome_index.tar.gz"
        }
    }
    ```

=== "Vorher"

    ```groovy title="nextflow.config" linenums="1"
    docker.enabled = true

    profiles {
        test {
            params.input = "${projectDir}/data/reads/ENCSR000COQ1_1.fastq.gz"
        }
    }
    ```

Der Parameter ist bereit; jetzt können wir den Alignment-Prozess erstellen.

### 3.2. Den Alignment-Prozess schreiben und im Workflow aufrufen

Wie zuvor müssen wir die Prozessdefinition ausfüllen, das Modul importieren und den Prozessaufruf hinzufügen.

#### 3.2.1. Das Modul für den Alignment-Prozess ausfüllen

Öffne `modules/hisat2_align.nf` und untersuche die Gliederung der Prozessdefinition.

Fülle die Prozessdefinition selbst aus, indem du die oben bereitgestellten Informationen verwendest, und überprüfe dann deine Arbeit anhand der Lösung im Tab "Danach" unten.

=== "Vorher"

    ```groovy title="modules/hisat2_align.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
     * Align reads to a reference genome
     */
    process HISAT2_ALIGN {

        container

        input:

        output:

        script:
        """

        """
    }
    ```

=== "Danach"

    ```groovy title="modules/hisat2_align.nf" linenums="1" hl_lines="8 11 12 15 16 20 21 22 23"
    #!/usr/bin/env nextflow

    /*
     * Align reads to a reference genome
     */
    process HISAT2_ALIGN {

        container "community.wave.seqera.io/library/hisat2_samtools:5e49f68a37dc010e"

        input:
        path reads
        path index_zip

        output:
        path "${reads.simpleName}.bam", emit: bam
        path "${reads.simpleName}.hisat2.log", emit: log

        script:
        """
        tar -xzvf ${index_zip}
        hisat2 -x ${index_zip.simpleName} -U ${reads} \
            --new-summary --summary-file ${reads.simpleName}.hisat2.log | \
            samtools view -bS -o ${reads.simpleName}.bam
        """
    }
    ```

Dieser Prozess nimmt zwei Eingaben: die Reads und das Genom-Index-Archiv.
Der Script-Block extrahiert zuerst den Index aus dem Archiv und führt dann das HISAT2-Alignment aus, das in `samtools view` geleitet wird, um die Ausgabe in das BAM-Format zu konvertieren.
Der `simpleName`-Accessor auf `index_zip` extrahiert den Basisnamen des Archivs (`genome_index`), um ihn als Index-Präfix zu verwenden.

#### 3.2.2. Das Modul einbinden

Aktualisiere `rnaseq.nf`, um das neue Modul zu importieren:

=== "Danach"

    ```groovy title="rnaseq.nf" linenums="3" hl_lines="4"
    // Modul-INCLUDE-Anweisungen
    include { FASTQC } from './modules/fastqc.nf'
    include { TRIM_GALORE } from './modules/trim_galore.nf'
    include { HISAT2_ALIGN } from './modules/hisat2_align.nf'
    ```

=== "Vorher"

    ```groovy title="rnaseq.nf" linenums="3"
    // Modul-INCLUDE-Anweisungen
    include { FASTQC } from './modules/fastqc.nf'
    include { TRIM_GALORE } from './modules/trim_galore.nf'
    ```

Als Nächstes fügen wir den Prozessaufruf zum Workflow hinzu.

#### 3.2.3. Den Alignment-Prozess aufrufen

Die getrimmten Reads befinden sich im `TRIM_GALORE.out.trimmed_reads`-Kanal, der vom vorherigen Schritt ausgegeben wurde.
Wir verwenden `#!groovy file(params.hisat2_index_zip)`, um das Genom-Index-Archiv bereitzustellen.

=== "Danach"

    ```groovy title="rnaseq.nf" linenums="19" hl_lines="14-15"
    workflow {

        main:
        // Eingabe-Channel aus einem Dateipfad erstellen
        read_ch = channel.fromPath(params.input)

        // Initiale Qualitätskontrolle
        FASTQC(read_ch)

        // Adapter-Trimming und Post-Trimming-QC
        TRIM_GALORE(read_ch)

        // Alignment zum Referenzgenom
        HISAT2_ALIGN(TRIM_GALORE.out.trimmed_reads, file(params.hisat2_index_zip))
    ```

=== "Vorher"

    ```groovy title="rnaseq.nf" linenums="19"
    workflow {

        main:
        // Eingabe-Channel aus einem Dateipfad erstellen
        read_ch = channel.fromPath(params.input)

        // Initiale Qualitätskontrolle
        FASTQC(read_ch)

        // Adapter-Trimming und Post-Trimming-QC
        TRIM_GALORE(read_ch)
    ```

Der Alignment-Prozess ist jetzt in den Workflow eingebunden.

### 3.3. Die Ausgabebehandlung aktualisieren

Wir müssen die Alignment-Ausgaben zur Publish-Deklaration hinzufügen und konfigurieren, wohin sie gehen.

#### 3.3.1. Publish-Ziele für die Alignment-Ausgaben hinzufügen

Füge die Alignment-Ausgaben zum `publish:`-Abschnitt hinzu:

=== "Danach"

    ```groovy title="rnaseq.nf" linenums="34" hl_lines="7-8"
        publish:
        fastqc_zip = FASTQC.out.zip
        fastqc_html = FASTQC.out.html
        trimmed_reads = TRIM_GALORE.out.trimmed_reads
        trimming_reports = TRIM_GALORE.out.trimming_reports
        trimming_fastqc = TRIM_GALORE.out.fastqc_reports
        bam = HISAT2_ALIGN.out.bam
        align_log = HISAT2_ALIGN.out.log
    }
    ```

=== "Vorher"

    ```groovy title="rnaseq.nf" linenums="34"
        publish:
        fastqc_zip = FASTQC.out.zip
        fastqc_html = FASTQC.out.html
        trimmed_reads = TRIM_GALORE.out.trimmed_reads
        trimming_reports = TRIM_GALORE.out.trimming_reports
        trimming_fastqc = TRIM_GALORE.out.fastqc_reports
    }
    ```

Als Nächstes müssen wir Nextflow mitteilen, wohin diese Ausgaben gelegt werden sollen.

#### 3.3.2. Die neuen Ausgabeziele konfigurieren

Füge Einträge für die Alignment-Ziele im `output {}`-Block hinzu und veröffentliche sie in ein `align/`-Unterverzeichnis:

=== "Danach"

    ```groovy title="rnaseq.nf" linenums="44" hl_lines="17-22"
    output {
        fastqc_zip {
            path 'fastqc'
        }
        fastqc_html {
            path 'fastqc'
        }
        trimmed_reads {
            path 'trimming'
        }
        trimming_reports {
            path 'trimming'
        }
        trimming_fastqc {
            path 'trimming'
        }
        bam {
            path 'align'
        }
        align_log {
            path 'align'
        }
    }
    ```

=== "Vorher"

    ```groovy title="rnaseq.nf" linenums="44"
    output {
        fastqc_zip {
            path 'fastqc'
        }
        fastqc_html {
            path 'fastqc'
        }
        trimmed_reads {
            path 'trimming'
        }
        trimming_reports {
            path 'trimming'
        }
        trimming_fastqc {
            path 'trimming'
        }
    }
    ```

Die Ausgabekonfiguration ist abgeschlossen.

### 3.4. Den Workflow ausführen

Der Workflow umfasst jetzt alle drei Verarbeitungsschritte: QC, Trimming und Alignment.

```bash
nextflow run rnaseq.nf -profile test
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `rnaseq.nf` [elated_stonebraker] DSL2 - revision: e8e57d0cdd

    executor >  local (3)
    [e8/fa29d6] FASTQC (1)       | 1 of 1 ✔
    [ca/ffdde2] TRIM_GALORE (1)  | 1 of 1 ✔
    [b6/1c6ca3] HISAT2_ALIGN (1) | 1 of 1 ✔
    ```

Du findest die Alignment-Ausgaben im Ergebnisverzeichnis.

```bash
ls results/align
```

```console title="Output"
ENCSR000COQ1_1_trimmed.bam  ENCSR000COQ1_1_trimmed.hisat2.log
```

Dies vervollständigt die grundlegende Verarbeitung, die wir auf jede Probe anwenden müssen.

_Wir werden die MultiQC-Report-Aggregation in Teil 3 hinzufügen, nachdem wir den Workflow so angepasst haben, dass er mehrere Proben gleichzeitig akzeptiert._

---

### Fazit

Du weißt jetzt, wie du alle Kernschritte zur individuellen Verarbeitung von Single-End-RNAseq-Proben zusammenfasst.

### Wie geht es weiter?

Mach eine Pause! Das war eine Menge.

Wenn du dich erfrischt fühlst, gehe weiter zu [Teil 3](./03_multi-sample.md), wo du lernst, wie du den Workflow anpasst, um mehrere Proben parallel zu verarbeiten, QC-Reports über alle Schritte für alle Proben zu aggregieren und die Ausführung des Workflows mit Paired-End-RNAseq-Daten zu ermöglichen.
