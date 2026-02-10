# Teil 2: Variantenerkennung pro Probe

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

In Teil 1 hast du die Samtools- und GATK-Befehle manuell in ihren jeweiligen Containern getestet.
Jetzt werden wir diese Befehle in einen Nextflow-Workflow einbinden.

## Aufgabe

In diesem Teil des Kurses entwickeln wir einen Workflow, der Folgendes tut:

1. Eine Indexdatei für jede BAM-Eingabedatei mit [Samtools](https://www.htslib.org/) generieren
2. GATK HaplotypeCaller auf jede BAM-Eingabedatei anwenden, um Variantenaufrufe pro Probe im VCF-Format (Variant Call Format) zu generieren

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/genomics/img/hello-gatk-1.svg"
</figure>

Dies repliziert die Schritte aus Teil 1, wo du diese Befehle manuell in ihren Containern ausgeführt hast.

Als Ausgangspunkt stellen wir dir eine Workflow-Datei `genomics.nf` zur Verfügung, die die Hauptteile des Workflows skizziert, sowie zwei Moduldateien, samtools_index.nf und gatk_haplotypecaller.nf, die die Struktur der Module umreißen.
Diese Dateien sind nicht funktionsfähig; ihr Zweck ist es lediglich, als Gerüst zu dienen, das du mit den interessanten Teilen des Codes ausfüllen kannst.

## Lernplan

Um den Entwicklungsprozess lehrreicher zu gestalten, haben wir dies in vier Schritte unterteilt:

1. **Einen einstufigen Workflow schreiben, der Samtools index auf eine BAM-Datei anwendet.**
   Dies umfasst das Erstellen eines Moduls, das Importieren und Aufrufen im Workflow.
2. **Einen zweiten Prozess hinzufügen, um GATK HaplotypeCaller auf die indizierte BAM-Datei anzuwenden.**
   Dies führt die Verkettung von Prozessausgaben zu Eingaben ein und behandelt Hilfsdateien.
3. **Den Workflow anpassen, um auf einem Batch von Proben zu laufen.**
   Dies behandelt parallele Ausführung und führt Tupel ein, um zugehörige Dateien zusammenzuhalten.
4. **Den Workflow so anpassen, dass er eine Textdatei mit einem Batch von Eingabedateien akzeptiert.**
   Dies demonstriert ein gängiges Muster zur Bereitstellung von Eingaben in großen Mengen.

Jeder Schritt konzentriert sich auf einen spezifischen Aspekt der Workflow-Entwicklung.

---

## 1. Einen einstufigen Workflow schreiben, der Samtools index auf eine BAM-Datei anwendet

Dieser erste Schritt konzentriert sich auf die Grundlagen: Laden einer BAM-Datei und Generieren eines Index dafür.

Erinnere dich an den `samtools index`-Befehl aus [Teil 1](01_method.md):

```bash
samtools index '<input_bam>'
```

Der Befehl nimmt eine BAM-Datei als Eingabe und erzeugt eine `.bai`-Indexdatei daneben.
Die Container-URI war `community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464`.

Wir werden diese Informationen nehmen und in Nextflow in drei Stufen einbinden:

1. Die Eingabe einrichten
2. Den Indizierungsprozess schreiben und im Workflow aufrufen
3. Die Ausgabeverarbeitung konfigurieren

### 1.1. Die Eingabe einrichten

Wir müssen einen Eingabeparameter deklarieren, ein Testprofil erstellen, um einen praktischen Standardwert bereitzustellen, und einen Eingabekanal erstellen.

#### 1.1.1. Eine Eingabeparameterdeklaration hinzufügen

In der Haupt-Workflow-Datei `genomics.nf` deklariere unter dem Abschnitt `Pipeline parameters` einen CLI-Parameter namens `reads_bam`.

=== "Danach"

    ```groovy title="genomics.nf" linenums="5" hl_lines="4-7"
    /*
     * Pipeline parameters
     */
    params {
        // Primary input
        reads_bam: Path
    }
    ```

=== "Vorher"

    ```groovy title="genomics.nf" linenums="5"
    /*
     * Pipeline parameters
     */

    // Primary input
    ```

Das richtet den CLI-Parameter ein, aber wir wollen nicht jedes Mal den Dateipfad eingeben, wenn wir den Workflow während der Entwicklung ausführen.
Es gibt mehrere Optionen, um einen Standardwert bereitzustellen; hier verwenden wir ein Testprofil.

#### 1.1.2. Ein Testprofil mit einem Standardwert in `nextflow.config` erstellen

Ein Testprofil bietet praktische Standardwerte zum Ausprobieren eines Workflows, ohne Eingaben auf der Kommandozeile anzugeben.
Dies ist eine gängige Konvention im Nextflow-Ökosystem (siehe [Hello Config](../../hello_nextflow/06_hello_config.md) für weitere Details).

Füge einen `profiles`-Block zu `nextflow.config` mit einem `test`-Profil hinzu, das den Parameter `reads_bam` auf eine der Test-BAM-Dateien setzt.

=== "Danach"

    ```groovy title="nextflow.config" linenums="1" hl_lines="3-7"
    docker.enabled = true

    profiles {
        test {
            params.reads_bam = "${projectDir}/data/bam/reads_mother.bam"
        }
    }
    ```

=== "Vorher"

    ```groovy title="nextflow.config" linenums="1"
    docker.enabled = true
    ```

Hier verwenden wir `${projectDir}`, eine eingebaute Nextflow-Variable, die auf das Verzeichnis zeigt, in dem sich das Workflow-Skript befindet.
Dies macht es einfach, auf Datendateien und andere Ressourcen zu verweisen, ohne absolute Pfade fest zu codieren.

#### 1.1.3. Den Eingabekanal einrichten

Erstelle im Workflow-Block einen Eingabekanal aus dem Parameterwert mit der `.fromPath`-Kanal-Factory (wie in [Hello Channels](../../hello_nextflow/02_hello_channels.md) verwendet).

=== "Danach"

    ```groovy title="genomics.nf" linenums="13" hl_lines="4-5"
    workflow {

        main:
        // Eingabekanal erstellen (einzelne Datei über CLI-Parameter)
        reads_ch = channel.fromPath(params.reads_bam)
    ```

=== "Vorher"

    ```groovy title="genomics.nf" linenums="13"
    workflow {

        main:
        // Create input channel
    ```

Jetzt müssen wir den Prozess erstellen, um die Indizierung auf diese Eingabe anzuwenden.

### 1.2. Den Indizierungsprozess schreiben und im Workflow aufrufen

Wir müssen die Prozessdefinition in der Moduldatei schreiben, sie mit einer Include-Anweisung in den Workflow importieren und sie auf die Eingabe anwenden.

#### 1.2.1. Das Modul für den Indizierungsprozess ausfüllen

Öffne `modules/samtools_index.nf` und untersuche die Gliederung der Prozessdefinition.
Du solltest die strukturellen Hauptelemente erkennen; falls nicht, lies [Hello Nextflow](../../hello_nextflow/01_hello_world.md) zur Auffrischung.

Fülle die Prozessdefinition selbst aus, indem du die oben bereitgestellten Informationen verwendest, und überprüfe dann deine Arbeit anhand der Lösung im Tab "Danach" unten.

=== "Vorher"

    ```groovy title="modules/samtools_index.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
     * Generate BAM index file
     */
    process SAMTOOLS_INDEX {

        container

        input:

        output:

        script:
        """

        """
    }
    ```

=== "Danach"

    ```groovy title="modules/samtools_index.nf" linenums="1" hl_lines="8 11 14 18"
    #!/usr/bin/env nextflow

    /*
     * BAM-Indexdatei generieren
     */
    process SAMTOOLS_INDEX {

        container 'community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464'

        input:
        path input_bam

        output:
        path "${input_bam}.bai"

        script:
        """
        samtools index '$input_bam'
        """
    }
    ```

Sobald du dies abgeschlossen hast, ist der Prozess fertig.
Um ihn im Workflow zu verwenden, musst du das Modul importieren und einen Prozessaufruf hinzufügen.

#### 1.2.2. Das Modul einbinden

Füge in `genomics.nf` eine `include`-Anweisung hinzu, um den Prozess für den Workflow verfügbar zu machen:

=== "Danach"

    ```groovy title="genomics.nf" linenums="3" hl_lines="2"
    // Module INCLUDE statements
    include { SAMTOOLS_INDEX } from './modules/samtools_index.nf'
    ```

=== "Vorher"

    ```groovy title="genomics.nf" linenums="3"
    // Module INCLUDE statements
    ```

Der Prozess ist jetzt im Workflow-Scope verfügbar.

#### 1.2.3. Den Indizierungsprozess auf die Eingabe anwenden

Fügen wir nun einen Aufruf von `SAMTOOLS_INDEX` im Workflow-Block hinzu und übergeben den Eingabekanal als Argument.

=== "Danach"

    ```groovy title="genomics.nf" linenums="14" hl_lines="7-8"
    workflow {

        main:
        // Eingabekanal erstellen (einzelne Datei über CLI-Parameter)
        reads_ch = channel.fromPath(params.reads_bam)

        // Indexdatei für Eingabe-BAM-Datei erstellen
        SAMTOOLS_INDEX(reads_ch)
    ```

=== "Vorher"

    ```groovy title="genomics.nf" linenums="14"
    workflow {

        main:
        // Eingabekanal erstellen (einzelne Datei über CLI-Parameter)
        reads_ch = channel.fromPath(params.reads_bam)

        // Call processes
    ```

Der Workflow lädt jetzt die Eingabe und führt den Indizierungsprozess darauf aus.
Als Nächstes müssen wir konfigurieren, wie die Ausgabe veröffentlicht wird.

### 1.3. Die Ausgabeverarbeitung konfigurieren

Wir müssen deklarieren, welche Prozessausgaben veröffentlicht werden sollen, und angeben, wohin sie gehen sollen.

#### 1.3.1. Eine Ausgabe im Abschnitt `publish:` deklarieren

Der Abschnitt `publish:` innerhalb des Workflow-Blocks deklariert, welche Prozessausgaben veröffentlicht werden sollen.
Weise die Ausgabe von `SAMTOOLS_INDEX` einem benannten Ziel namens `bam_index` zu.

=== "Danach"

    ```groovy title="genomics.nf" linenums="22" hl_lines="2"
        publish:
        bam_index = SAMTOOLS_INDEX.out
    }
    ```

=== "Vorher"

    ```groovy title="genomics.nf" linenums="22"
        publish:
        // Declare outputs to publish
    }
    ```

Jetzt müssen wir Nextflow mitteilen, wohin die veröffentlichte Ausgabe gelegt werden soll.

#### 1.3.2. Das Ausgabeziel im Block `output {}` konfigurieren

Der Block `output {}` befindet sich außerhalb des Workflows und gibt an, wo jedes benannte Ziel veröffentlicht wird.
Fügen wir ein Ziel für `bam_index` hinzu, das in ein Unterverzeichnis `bam/` veröffentlicht.

=== "Danach"

    ```groovy title="genomics.nf" linenums="26" hl_lines="2-4"
    output {
        bam_index {
            path 'bam'
        }
    }
    ```

=== "Vorher"

    ```groovy title="genomics.nf" linenums="26"
    output {
        // Configure publish targets
    }
    ```

!!! note "Hinweis"

    Standardmäßig veröffentlicht Nextflow Ausgabedateien als symbolische Links, was unnötige Duplizierung vermeidet.
    Obwohl die Datendateien, die wir hier verwenden, sehr klein sind, können sie in der Genomik sehr groß werden.
    Symlinks funktionieren nicht mehr, wenn du dein `work`-Verzeichnis aufräumst. Für Produktions-Workflows möchtest du daher möglicherweise den Standard-Veröffentlichungsmodus auf `'copy'` überschreiben.

### 1.4. Den Workflow ausführen

An diesem Punkt haben wir einen einstufigen Indizierungs-Workflow, der voll funktionsfähig sein sollte. Testen wir, ob er funktioniert!

Wir können ihn mit `-profile test` ausführen, um den im Testprofil eingerichteten Standardwert zu verwenden und zu vermeiden, den Pfad auf der Kommandozeile schreiben zu müssen.

```bash
nextflow run genomics.nf -profile test
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics.nf` [reverent_sinoussi] DSL2 - revision: 41d43ad7fe

    executor >  local (1)
    [2a/e69536] SAMTOOLS_INDEX (1) | 1 of 1 ✔
    ```

Du kannst überprüfen, dass die Indexdatei korrekt generiert wurde, indem du im Arbeitsverzeichnis oder im Ergebnisverzeichnis nachsiehst.

??? abstract "Work-Verzeichnis Inhalt"

    ```console
    work/2a/e695367b2f60df09cf826b07192dc3
    ├── reads_mother.bam -> /workspaces/training/nf4-science/genomics/data/bam/reads_mother.bam
    └── reads_mother.bam.bai
    ```

??? abstract "Results-Verzeichnis Inhalt"

    ```console
    results/
    └── bam/
        └── reads_mother.bam.bai -> ...
    ```

Da ist sie!

### Fazit

Du weißt, wie man ein Modul mit einem Prozess erstellt, es in einen Workflow importiert, es mit einem Eingabekanal aufruft und die Ergebnisse veröffentlicht.

### Wie geht es weiter?

Füge einen zweiten Schritt hinzu, der die Ausgabe des Indizierungsprozesses nimmt und sie verwendet, um Variantenerkennung durchzuführen.

---

## 2. Einen zweiten Prozess hinzufügen, um GATK HaplotypeCaller auf die indizierte BAM-Datei anzuwenden

Jetzt, da wir einen Index für unsere Eingabedatei haben, können wir mit der Einrichtung des Variantenerkennungsschritts fortfahren.

Erinnere dich an den `gatk HaplotypeCaller`-Befehl aus [Teil 1](01_method.md):

```bash
gatk HaplotypeCaller \
        -R /data/ref/ref.fasta \
        -I /data/bam/reads_mother.bam \
        -O reads_mother.vcf \
        -L /data/ref/intervals.bed
```

Der Befehl nimmt eine BAM-Datei (`-I`), ein Referenzgenom (`-R`) und eine Intervalldatei (`-L`) und erzeugt eine VCF-Datei (`-O`) zusammen mit ihrem Index.
Das Tool erwartet auch, dass der BAM-Index, der Referenzindex und das Referenzwörterbuch zusammen mit ihren jeweiligen Dateien liegen.
Die Container-URI war `community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867`.

Wir folgen denselben drei Stufen wie zuvor:

1. Die Eingaben einrichten
2. Den Variantenerkennungsprozess schreiben und im Workflow aufrufen
3. Die Ausgabeverarbeitung konfigurieren

### 2.1. Die Eingaben einrichten

Der Variantenerkennungsschritt erfordert mehrere zusätzliche Eingabedateien.
Wir müssen Parameter dafür deklarieren, Standardwerte zum Testprofil hinzufügen und Variablen erstellen, um sie zu laden.

#### 2.1.1. Parameterdeklarationen für Hilfseingaben hinzufügen

Da unser neuer Prozess eine Handvoll zusätzlicher Dateien erwartet, füge Parameterdeklarationen dafür in `genomics.nf` unter dem Abschnitt `Pipeline parameters` hinzu:

=== "Danach"

    ```groovy title="genomics.nf" linenums="9" hl_lines="5-9"
    params {
        // Primary input
        reads_bam: Path

        // Hilfsdateien
        reference: Path
        reference_index: Path
        reference_dict: Path
        intervals: Path
    }
    ```

=== "Vorher"

    ```groovy title="genomics.nf" linenums="9"
    params {
        // Primary input
        reads_bam: Path
    }
    ```

Wie zuvor stellen wir Standardwerte über das Testprofil bereit, anstatt sie inline anzugeben.

#### 2.1.2. Standardwerte für Hilfsdateien zum Testprofil hinzufügen

Genau wie wir es für `reads_bam` in Abschnitt 1.1.2 getan haben, füge Standardwerte für die Hilfsdateien zum Testprofil in `nextflow.config` hinzu:

=== "Danach"

    ```groovy title="nextflow.config" linenums="4" hl_lines="3-6"
    test {
        params.reads_bam = "${projectDir}/data/bam/reads_mother.bam"
        params.reference = "${projectDir}/data/ref/ref.fasta"
        params.reference_index = "${projectDir}/data/ref/ref.fasta.fai"
        params.reference_dict = "${projectDir}/data/ref/ref.dict"
        params.intervals = "${projectDir}/data/ref/intervals.bed"
    }
    ```

=== "Vorher"

    ```groovy title="nextflow.config" linenums="4"
    test {
        params.reads_bam = "${projectDir}/data/bam/reads_mother.bam"
    }
    ```

Jetzt müssen wir Variablen erstellen, die diese Dateipfade zur Verwendung im Workflow laden.

#### 2.1.3. Variablen für die Hilfsdateien erstellen

Füge Variablen für die Hilfsdateipfade innerhalb des Workflow-Blocks hinzu:

=== "Danach"

    ```groovy title="genomics.nf" linenums="21" hl_lines="7-11"
    workflow {

        main:
        // Eingabekanal erstellen (einzelne Datei über CLI-Parameter)
        reads_ch = channel.fromPath(params.reads_bam)

        // Dateipfade für die Hilfsdateien laden (Referenz und Intervalle)
        ref_file        = file(params.reference)
        ref_index_file  = file(params.reference_index)
        ref_dict_file   = file(params.reference_dict)
        intervals_file  = file(params.intervals)

        // Indexdatei für Eingabe-BAM-Datei erstellen
        SAMTOOLS_INDEX(reads_ch)
    ```

=== "Vorher"

    ```groovy title="genomics.nf" linenums="21"
    workflow {

        main:
        // Eingabekanal erstellen (einzelne Datei über CLI-Parameter)
        reads_ch = channel.fromPath(params.reads_bam)

        // Indexdatei für Eingabe-BAM-Datei erstellen
        SAMTOOLS_INDEX(reads_ch)
    ```

Die `file()`-Syntax teilt Nextflow explizit mit, diese Eingaben als Dateipfade zu behandeln.
Mehr darüber erfährst du in der Side Quest [Working with files](../../side_quests/working_with_files.md).

### 2.2. Den Variantenerkennungsprozess schreiben und im Workflow aufrufen

Wir müssen die Prozessdefinition in der Moduldatei schreiben, sie mit einer Include-Anweisung in den Workflow importieren und sie auf die Eingabe-Reads plus die Ausgabe des Indizierungsschritts und die Hilfsdateien anwenden.

#### 2.2.1. Das Modul für den Variantenerkennungsprozess ausfüllen

Öffne `modules/gatk_haplotypecaller.nf` und untersuche die Gliederung der Prozessdefinition.

Fülle die Prozessdefinition selbst aus, indem du die oben bereitgestellten Informationen verwendest, und überprüfe dann deine Arbeit anhand der Lösung im Tab "Danach" unten.

=== "Vorher"

    ```groovy title="modules/gatk_haplotypecaller.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
     * Call variants with GATK HaplotypeCaller
     */
    process GATK_HAPLOTYPECALLER {

        container

        input:

        output:

        script:
        """

        """
    }
    ```

=== "Danach"

    ```groovy title="modules/gatk_haplotypecaller.nf" linenums="1" hl_lines="8 11-16 19-20 24-28"
    #!/usr/bin/env nextflow

    /*
     * Varianten mit GATK HaplotypeCaller aufrufen
     */
    process GATK_HAPLOTYPECALLER {

        container "community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867"

        input:
        path input_bam
        path input_bam_index
        path ref_fasta
        path ref_index
        path ref_dict
        path interval_list

        output:
        path "${input_bam}.vcf"     , emit: vcf
        path "${input_bam}.vcf.idx" , emit: idx

        script:
        """
        gatk HaplotypeCaller \
            -R ${ref_fasta} \
            -I ${input_bam} \
            -O ${input_bam}.vcf \
            -L ${interval_list}
        """
    }
    ```

Du wirst bemerken, dass dieser Prozess mehr Eingaben hat, als der GATK-Befehl selbst benötigt.
GATK weiß aufgrund von Namenskonventionen, wo es nach der BAM-Indexdatei und den Hilfsdateien des Referenzgenoms suchen muss, aber Nextflow ist domänenunabhängig und kennt diese Konventionen nicht.
Wir müssen sie explizit auflisten, damit Nextflow sie zur Laufzeit im Arbeitsverzeichnis bereitstellt; andernfalls wird GATK einen Fehler über fehlende Dateien ausgeben.

Ebenso listen wir die Indexdatei der Ausgabe-VCF (`"${input_bam}.vcf.idx"`) explizit auf, damit Nextflow sie für nachfolgende Schritte verfolgt.
Wir verwenden die `emit:`-Syntax, um jedem Ausgabekanal einen Namen zuzuweisen, was nützlich wird, wenn wir die Ausgaben in den Publish-Block einbinden.

Sobald du dies abgeschlossen hast, ist der Prozess fertig.
Um ihn im Workflow zu verwenden, musst du das Modul importieren und einen Prozessaufruf hinzufügen.

#### 2.2.2. Das neue Modul importieren

Aktualisiere `genomics.nf`, um das neue Modul zu importieren:

=== "Danach"

    ```groovy title="genomics.nf" linenums="3" hl_lines="3"
    // Module INCLUDE statements
    include { SAMTOOLS_INDEX } from './modules/samtools_index.nf'
    include { GATK_HAPLOTYPECALLER } from './modules/gatk_haplotypecaller.nf'
    ```

=== "Vorher"

    ```groovy title="genomics.nf" linenums="3"
    // Module INCLUDE statements
    include { SAMTOOLS_INDEX } from './modules/samtools_index.nf'
    ```

Der Prozess ist jetzt im Workflow-Scope verfügbar.

#### 2.2.3. Den Prozessaufruf hinzufügen

Füge den Prozessaufruf im Workflow-Body unter `main:` hinzu:

=== "Danach"

    ```groovy title="genomics.nf" linenums="33" hl_lines="4-12"
        // Indexdatei für Eingabe-BAM-Datei erstellen
        SAMTOOLS_INDEX(reads_ch)

        // Varianten aus der indizierten BAM-Datei aufrufen
        GATK_HAPLOTYPECALLER(
            reads_ch,
            SAMTOOLS_INDEX.out,
            ref_file,
            ref_index_file,
            ref_dict_file,
            intervals_file
        )
    ```

=== "Vorher"

    ```groovy title="genomics.nf" linenums="33"
        // Indexdatei für Eingabe-BAM-Datei erstellen
        SAMTOOLS_INDEX(reads_ch)
    ```

Du solltest die `*.out`-Syntax aus der Hello Nextflow-Trainingsreihe erkennen; wir sagen Nextflow, dass es den von `SAMTOOLS_INDEX` ausgegebenen Kanal nehmen und in den `GATK_HAPLOTYPECALLER`-Prozessaufruf einstecken soll.

!!! note "Hinweis"

    Beachte, dass die Eingaben im Aufruf des Prozesses in genau derselben Reihenfolge bereitgestellt werden, wie sie im Eingabeblock des Prozesses aufgelistet sind.
    In Nextflow sind Eingaben positionsabhängig, was bedeutet, dass du _dieselbe Reihenfolge_ einhalten musst; und natürlich muss es dieselbe Anzahl von Elementen geben.

### 2.3. Die Ausgabeverarbeitung konfigurieren

Wir müssen die neuen Ausgaben zur Publish-Deklaration hinzufügen und konfigurieren, wohin sie gehen.

#### 2.3.1. Publish-Ziele für die Variantenerkennungsausgaben hinzufügen

Füge die VCF- und Index-Ausgaben zum Abschnitt `publish:` hinzu:

=== "Danach"

    ```groovy title="genomics.nf" linenums="45" hl_lines="3-4"
        publish:
        bam_index = SAMTOOLS_INDEX.out
        vcf = GATK_HAPLOTYPECALLER.out.vcf
        vcf_idx = GATK_HAPLOTYPECALLER.out.idx
    }
    ```

=== "Vorher"

    ```groovy title="genomics.nf" linenums="45"
        publish:
        bam_index = SAMTOOLS_INDEX.out
    }
    ```

Jetzt müssen wir Nextflow mitteilen, wohin die neuen Ausgaben gelegt werden sollen.

#### 2.3.2. Die neuen Ausgabeziele konfigurieren

Füge Einträge für die Ziele `vcf` und `vcf_idx` im Block `output {}` hinzu und veröffentliche beide in ein Unterverzeichnis `vcf/`:

=== "Danach"

    ```groovy title="genomics.nf" linenums="51" hl_lines="5-10"
    output {
        bam_index {
            path 'bam'
        }
        vcf {
            path 'vcf'
        }
        vcf_idx {
            path 'vcf'
        }
    }
    ```

=== "Vorher"

    ```groovy title="genomics.nf" linenums="49"
    output {
        bam_index {
            path 'bam'
        }
    }
    ```

Die VCF und ihr Index werden als separate Ziele veröffentlicht, die beide in das Unterverzeichnis `vcf/` gehen.

### 2.4. Den Workflow ausführen

Führe den erweiterten Workflow aus und füge diesmal `-resume` hinzu, damit wir den Indizierungsschritt nicht erneut ausführen müssen.

```bash
nextflow run genomics.nf -profile test -resume
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics.nf` [grave_volta] DSL2 - revision: 4790abc96a

    executor >  local (1)
    [2a/e69536] SAMTOOLS_INDEX (1)       | 1 of 1, cached: 1 ✔
    [53/e18e98] GATK_HAPLOTYPECALLER (1) | 1 of 1 ✔
    ```

Wenn wir uns jetzt die Konsolenausgabe ansehen, sehen wir die beiden aufgelisteten Prozesse.

Der erste Prozess wurde dank des Cachings übersprungen, wie erwartet, während der zweite Prozess ausgeführt wurde, da er brandneu ist.

Du findest die Ausgabedateien im Ergebnisverzeichnis (als symbolische Links zum Arbeitsverzeichnis).

??? abstract "Verzeichnisinhalt"

    ```console
    results/
    ├── bam/
    │   └── reads_mother.bam.bai -> ...
    └── vcf/
        ├── reads_mother.bam.vcf -> ...
        └── reads_mother.bam.vcf.idx -> ...
    ```

Wenn du die VCF-Datei öffnest, solltest du denselben Inhalt sehen wie in der Datei, die du durch direktes Ausführen des GATK-Befehls im Container generiert hast.

??? abstract "Dateiinhalt"

    ```console title="reads_mother.bam.vcf" linenums="26"
    #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	reads_mother
    20_10037292_10066351	3480	.	C	CT	503.03	.	AC=2;AF=1.00;AN=2;DP=23;ExcessHet=0.0000;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=27.95;SOR=1.179	GT:AD:DP:GQ:PL	1/1:0,18:18:54:517,54,0
    20_10037292_10066351	3520	.	AT	A	609.03	.	AC=2;AF=1.00;AN=2;DP=18;ExcessHet=0.0000;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=33.83;SOR=0.693	GT:AD:DP:GQ:PL	1/1:0,18:18:54:623,54,0
    20_10037292_10066351	3529	.	T	A	155.64	.	AC=1;AF=0.500;AN=2;BaseQRankSum=-0.544;DP=21;ExcessHet=0.0000;FS=1.871;MLEAC=1;MLEAF=0.500;MQ=60.00;MQRankSum=0.000;QD=7.78;ReadPosRankSum=-1.158;SOR=1.034	GT:AD:DP:GQ:PL	0/1:12,8:20:99:163,0,328
    ```

Dies ist die Ausgabe, die wir für jede Probe in unserer Studie generieren möchten.

### Fazit

Du weißt, wie man einen zweistufigen modularen Workflow erstellt, der echte Analysearbeit leistet und mit Eigenheiten von Genomik-Dateiformaten wie den Hilfsdateien umgehen kann.

### Wie geht es weiter?

Mache den Workflow so, dass er mehrere Proben in großen Mengen verarbeitet.

---

## 3. Den Workflow anpassen, um auf einem Batch von Proben zu laufen

Es ist schön und gut, einen Workflow zu haben, der die Verarbeitung einer einzelnen Probe automatisieren kann, aber was ist, wenn du 1000 Proben hast?
Musst du ein Bash-Skript schreiben, das durch alle deine Proben loopt?

Nein, zum Glück nicht! Mache einfach eine kleine Anpassung am Code und Nextflow wird das auch für dich erledigen.

### 3.1. Die Eingabe aktualisieren, um drei Proben aufzulisten

Um auf mehreren Proben zu laufen, aktualisiere das Testprofil, um ein Array von Dateipfaden anstelle eines einzelnen bereitzustellen.
Dies ist eine schnelle Möglichkeit, die Multi-Proben-Ausführung zu testen; im nächsten Schritt wechseln wir zu einem skalierbareren Ansatz mit einer Datei von Eingaben.

Kommentiere zuerst die Typannotation in der Parameterdeklaration aus, da Arrays keine typisierten Deklarationen verwenden können:

=== "Danach"

    ```groovy title="genomics.nf" linenums="10" hl_lines="1-2"
        // Primäre Eingabe (Array von drei Proben)
        reads_bam //: Path
    ```

=== "Vorher"

    ```groovy title="genomics.nf" linenums="10"
        // Primary input
        reads_bam: Path
    ```

Aktualisiere dann das Testprofil, um alle drei Proben aufzulisten:

=== "Danach"

    ```groovy title="nextflow.config" linenums="4" hl_lines="2-6"
    test {
        params.reads_bam = [
            "${projectDir}/data/bam/reads_mother.bam",
            "${projectDir}/data/bam/reads_father.bam",
            "${projectDir}/data/bam/reads_son.bam"
        ]
        params.reference = "${projectDir}/data/ref/ref.fasta"
        params.reference_index = "${projectDir}/data/ref/ref.fasta.fai"
        params.reference_dict = "${projectDir}/data/ref/ref.dict"
        params.intervals = "${projectDir}/data/ref/intervals.bed"
    }
    ```

=== "Vorher"

    ```groovy title="nextflow.config" linenums="4"
    test {
        params.reads_bam = "${projectDir}/data/bam/reads_mother.bam"
        params.reference = "${projectDir}/data/ref/ref.fasta"
        params.reference_index = "${projectDir}/data/ref/ref.fasta.fai"
        params.reference_dict = "${projectDir}/data/ref/ref.dict"
        params.intervals = "${projectDir}/data/ref/intervals.bed"
    }
    ```

Die Kanal-Factory im Workflow-Body (`.fromPath`) akzeptiert mehrere Dateipfade genauso gut wie einen einzelnen, daher sind keine weiteren Änderungen erforderlich.

### 3.2. Den Workflow ausführen

Versuche jetzt, den Workflow auszuführen, nachdem die Verkabelung eingerichtet ist, um auf allen drei Testproben zu laufen.

```bash
nextflow run genomics.nf -profile test -resume
```

Lustige Sache: Dies _könnte funktionieren_, ODER es _könnte fehlschlagen_. Hier ist zum Beispiel ein Lauf, der erfolgreich war:

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics.nf` [peaceful_yalow] DSL2 - revision: a256d113ad

    executor >  local (6)
    [4f/7071b0] SAMTOOLS_INDEX (3)       | 3 of 3, cached: 1 ✔
    [7a/89bc43] GATK_HAPLOTYPECALLER (2) | 3 of 3, cached: 1 ✔
    ```

Wenn dein Workflow-Lauf erfolgreich war, führe ihn erneut aus, bis du einen Fehler wie diesen erhältst:

??? failure "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics.nf` [loving_pasteur] DSL2 - revision: d2a8e63076

    executor >  local (4)
    [01/eea165] SAMTOOLS_INDEX (2)       | 3 of 3, cached: 1 ✔
    [a5/fa9fd0] GATK_HAPLOTYPECALLER (3) | 1 of 3, cached: 1
    ERROR ~ Error executing process > 'GATK_HAPLOTYPECALLER (2)'

    Caused by:
      Process `GATK_HAPLOTYPECALLER (2)` terminated with an error exit status (2)

    Command executed:

      gatk HaplotypeCaller         -R ref.fasta         -I reads_father.bam         -O reads_father.bam.vcf         -L intervals.bed

    Command exit status:
      2

    Command error:
      ...
      A USER ERROR has occurred: Traversal by intervals was requested but some input files are not indexed.
      ...
    ```

Wenn du dir die GATK-Befehlsfehlerausgabe ansiehst, wird es eine Zeile wie diese geben:

```console
A USER ERROR has occurred: Traversal by intervals was requested but some input files are not indexed.
```

Nun, das ist seltsam, wenn man bedenkt, dass wir die BAM-Dateien im ersten Schritt des Workflows explizit indiziert haben. Könnte etwas mit der Verkabelung nicht stimmen?

### 3.3. Das Problem beheben

Wir werden die Arbeitsverzeichnisse inspizieren und den Operator `view()` verwenden, um herauszufinden, was schief gelaufen ist.

#### 3.3.1. Die Arbeitsverzeichnisse für die relevanten Aufrufe überprüfen

Schaue in das Arbeitsverzeichnis für den fehlgeschlagenen `GATK_HAPLOTYPECALLER`-Prozessaufruf, der in der Konsolenausgabe aufgeführt ist.

??? abstract "Verzeichnisinhalt"

    ```console
    work/a5/fa9fd0994b6beede5fb9ea073596c2
    ├── intervals.bed -> /workspaces/training/nf4-science/genomics/data/ref/intervals.bed
    ├── reads_father.bam.bai -> /workspaces/training/nf4-science/genomics/work/01/eea16597bd6e810fb4cf89e60f8c2d/reads_father.bam.bai
    ├── reads_son.bam -> /workspaces/training/nf4-science/genomics/data/bam/reads_son.bam
    ├── reads_son.bam.vcf
    ├── reads_son.bam.vcf.idx
    ├── ref.dict -> /workspaces/training/nf4-science/genomics/data/ref/ref.dict
    ├── ref.fasta -> /workspaces/training/nf4-science/genomics/data/ref/ref.fasta
    └── ref.fasta.fai -> /workspaces/training/nf4-science/genomics/data/ref/ref.fasta.fai
    ```

Achte besonders auf die Namen der BAM-Datei und des BAM-Index, die in diesem Verzeichnis aufgeführt sind: `reads_son.bam` und `reads_father.bam.bai`.

Was zum Teufel? Nextflow hat eine Indexdatei im Arbeitsverzeichnis dieses Prozessaufrufs bereitgestellt, aber es ist die falsche. Wie konnte das passieren?

#### 3.3.2. Den [view()-Operator](https://www.nextflow.io/docs/latest/reference/operator.html#view) verwenden, um Kanalinhalte zu inspizieren

Füge diese beiden Zeilen im Workflow-Body vor dem `GATK_HAPLOTYPECALLER`-Prozessaufruf hinzu, um den Inhalt des Kanals anzuzeigen:

=== "Danach"

    ```groovy title="genomics.nf" hl_lines="3-5"
        SAMTOOLS_INDEX(reads_ch)

        // temporäre Diagnose
        reads_ch.view()
        SAMTOOLS_INDEX.out.view()

        // Varianten aus der indizierten BAM-Datei aufrufen
        GATK_HAPLOTYPECALLER(
    ```

=== "Vorher"

    ```groovy title="genomics.nf"
        SAMTOOLS_INDEX(reads_ch)

        // Varianten aus der indizierten BAM-Datei aufrufen
        GATK_HAPLOTYPECALLER(
    ```

Führe dann den Workflow-Befehl erneut aus.

```bash
nextflow run genomics.nf -profile test
```

Auch hier kann dies erfolgreich sein oder fehlschlagen. So sieht die Ausgabe der beiden `.view()`-Aufrufe für einen fehlgeschlagenen Lauf aus:

```console
/workspaces/training/nf4-science/genomics/data/bam/reads_mother.bam
/workspaces/training/nf4-science/genomics/data/bam/reads_father.bam
/workspaces/training/nf4-science/genomics/data/bam/reads_son.bam
/workspaces/training/nf4-science/genomics/work/9c/53492e3518447b75363e1cd951be4b/reads_father.bam.bai
/workspaces/training/nf4-science/genomics/work/cc/37894fffdf6cc84c3b0b47f9b536b7/reads_son.bam.bai
/workspaces/training/nf4-science/genomics/work/4d/dff681a3d137ba7d9866e3d9307bd0/reads_mother.bam.bai
```

Die ersten drei Zeilen entsprechen dem Eingabekanal und die zweiten dem Ausgabekanal.
Du kannst sehen, dass die BAM-Dateien und Indexdateien für die drei Proben nicht in derselben Reihenfolge aufgelistet sind!

!!! note "Hinweis"

    Wenn du einen Nextflow-Prozess auf einem Kanal aufrufst, der mehrere Elemente enthält, wird Nextflow versuchen, die Ausführung so weit wie möglich zu parallelisieren, und wird Ausgaben in der Reihenfolge sammeln, in der sie verfügbar werden.
    Die Konsequenz ist, dass die entsprechenden Ausgaben möglicherweise in einer anderen Reihenfolge gesammelt werden als die ursprünglichen Eingaben eingegeben wurden.

Wie derzeit geschrieben, geht unser Workflow-Skript davon aus, dass die Indexdateien aus dem Indizierungsschritt in derselben Mutter/Vater/Sohn-Reihenfolge herauskommen, wie die Eingaben gegeben wurden.
Aber das ist nicht garantiert, weshalb manchmal (aber nicht immer) die falschen Dateien im zweiten Schritt gepaart werden.

Um dies zu beheben, müssen wir sicherstellen, dass die BAM-Dateien und ihre Indexdateien zusammen durch die Kanäle reisen.

!!! tip "Tipp"

    Die `view()`-Anweisungen im Workflow-Code tun nichts, daher ist es kein Problem, sie drin zu lassen.
    Sie werden jedoch deine Konsolenausgabe überladen, daher empfehlen wir, sie zu entfernen, wenn du mit der Fehlerbehebung fertig bist.

### 3.4. Den Workflow aktualisieren, um die Indexdateien korrekt zu handhaben

Die Lösung besteht darin, jede BAM-Datei mit ihrem Index in ein Tupel zu bündeln und dann den nachgelagerten Prozess und die Workflow-Verkabelung entsprechend zu aktualisieren.

#### 3.4.1. Die Ausgabe des SAMTOOLS_INDEX-Moduls in ein Tupel ändern

Der einfachste Weg, um sicherzustellen, dass eine BAM-Datei und ihr Index eng verbunden bleiben, besteht darin, sie zusammen in ein Tupel zu packen, das aus der Indexaufgabe herauskommt.

!!! note "Hinweis"

    Ein **Tupel** ist eine endliche, geordnete Liste von Elementen, die häufig verwendet wird, um mehrere Werte von einer Funktion zurückzugeben. Tupel sind besonders nützlich, um mehrere Eingaben oder Ausgaben zwischen Prozessen zu übergeben und dabei ihre Zuordnung und Reihenfolge beizubehalten.

Aktualisiere die Ausgabe in `modules/samtools_index.nf`, um die BAM-Datei einzuschließen:

=== "Danach"

    ```groovy title="modules/samtools_index.nf" linenums="14" hl_lines="2"
        output:
        tuple path(input_bam), path("${input_bam}.bai")
    ```

=== "Vorher"

    ```groovy title="modules/samtools_index.nf" linenums="14"
        output:
        path "${input_bam}.bai"
    ```

Auf diese Weise wird jede Indexdatei eng mit ihrer ursprünglichen BAM-Datei gekoppelt, und die Gesamtausgabe des Indizierungsschritts wird ein einzelner Kanal sein, der Dateipaare enthält.

#### 3.4.2. Die Eingabe des GATK_HAPLOTYPECALLER-Moduls ändern, um ein Tupel zu akzeptieren

Da wir die 'Form' der Ausgabe des ersten Prozesses geändert haben, müssen wir die Eingabedefinition des zweiten Prozesses entsprechend aktualisieren.

Aktualisiere `modules/gatk_haplotypecaller.nf`:

=== "Danach"

    ```groovy title="modules/gatk_haplotypecaller.nf" linenums="11" hl_lines="2"
        input:
        tuple path(input_bam), path(input_bam_index)
    ```

=== "Vorher"

    ```groovy title="modules/gatk_haplotypecaller.nf" linenums="11"
        input:
        path input_bam
        path input_bam_index
    ```

Jetzt müssen wir den Workflow aktualisieren, um die neue Tupelstruktur im Prozessaufruf und den Publish-Zielen widerzuspiegeln.

#### 3.4.3. Den Aufruf von GATK_HAPLOTYPECALLER im Workflow aktualisieren

Wir müssen dem `GATK_HAPLOTYPECALLER`-Prozess nicht mehr den ursprünglichen `reads_ch` bereitstellen, da die BAM-Datei jetzt in die Kanalausgabe von `SAMTOOLS_INDEX` gebündelt ist.

Aktualisiere den Aufruf in `genomics.nf`:

=== "Danach"

    ```groovy title="genomics.nf" linenums="42" hl_lines="2"
        GATK_HAPLOTYPECALLER(
            SAMTOOLS_INDEX.out,
    ```

=== "Vorher"

    ```groovy title="genomics.nf" linenums="42"
        GATK_HAPLOTYPECALLER(
            reads_ch,
            SAMTOOLS_INDEX.out,
    ```

Schließlich müssen wir die Publish-Ziele aktualisieren, um die neue Ausgabestruktur widerzuspiegeln.

#### 3.4.4. Das Publish-Ziel für die indizierte BAM-Ausgabe aktualisieren

Da die SAMTOOLS_INDEX-Ausgabe jetzt ein Tupel ist, das sowohl die BAM-Datei als auch ihren Index enthält, benenne das Publish-Ziel von `bam_index` in `indexed_bam` um, um seinen Inhalt besser widerzuspiegeln:

=== "Danach"

    ```groovy title="genomics.nf" linenums="46" hl_lines="2 8"
        publish:
        indexed_bam = SAMTOOLS_INDEX.out
        vcf = GATK_HAPLOTYPECALLER.out.vcf
        vcf_idx = GATK_HAPLOTYPECALLER.out.idx
    }

    output {
        indexed_bam {
            path 'bam'
        }
        vcf {
            path 'vcf'
        }
        vcf_idx {
            path 'vcf'
        }
    }
    ```

=== "Vorher"

    ```groovy title="genomics.nf" linenums="46"
        publish:
        bam_index = SAMTOOLS_INDEX.out
        vcf = GATK_HAPLOTYPECALLER.out.vcf
        vcf_idx = GATK_HAPLOTYPECALLER.out.idx
    }

    output {
        bam_index {
            path 'bam'
        }
        vcf {
            path 'vcf'
        }
        vcf_idx {
            path 'vcf'
        }
    }
    ```

Mit diesen Änderungen ist garantiert, dass die BAM und ihr Index zusammen reisen, sodass die Paarung immer korrekt ist.

### 3.5. Den korrigierten Workflow ausführen

Führe den Workflow erneut aus, um sicherzustellen, dass dies in Zukunft zuverlässig funktioniert.

```bash
nextflow run genomics.nf -profile test
```

Diesmal (und jedes Mal) sollte alles korrekt laufen:

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics.nf` [special_goldstine] DSL2 - revision: 4cbbf6ea3e

    executor >  local (6)
    [d6/10c2c4] SAMTOOLS_INDEX (1)       | 3 of 3 ✔
    [88/1783aa] GATK_HAPLOTYPECALLER (2) | 3 of 3 ✔
    ```

Das Ergebnisverzeichnis enthält jetzt sowohl BAM- als auch BAI-Dateien für jede Probe (aus dem Tupel) zusammen mit den VCF-Ausgaben:

??? abstract "Results-Verzeichnis Inhalt"

    ```console
    results/
    ├── bam/
    │   ├── reads_father.bam -> ...
    │   ├── reads_father.bam.bai -> ...
    │   ├── reads_mother.bam -> ...
    │   ├── reads_mother.bam.bai -> ...
    │   ├── reads_son.bam -> ...
    │   └── reads_son.bam.bai -> ...
    └── vcf/
        ├── reads_father.bam.vcf -> ...
        ├── reads_father.bam.vcf.idx -> ...
        ├── reads_mother.bam.vcf -> ...
        ├── reads_mother.bam.vcf.idx -> ...
        ├── reads_son.bam.vcf -> ...
        └── reads_son.bam.vcf.idx -> ...
    ```

Durch das Bündeln zugehöriger Dateien in Tupel haben wir sichergestellt, dass die richtigen Dateien immer zusammen durch den Workflow reisen.
Der Workflow verarbeitet jetzt zuverlässig eine beliebige Anzahl von Proben, aber sie einzeln in der Config aufzulisten ist nicht sehr skalierbar.
Im nächsten Schritt wechseln wir zum Lesen von Eingaben aus einer Datei.

### Fazit

Du weißt, wie du deinen Workflow auf mehreren Proben (unabhängig) laufen lässt.

### Wie geht es weiter?

Mache es einfacher, Proben in großen Mengen zu handhaben.

---

## 4. Den Workflow so anpassen, dass er eine Textdatei mit einem Batch von Eingabedateien akzeptiert

Eine sehr gängige Methode, um mehrere Dateneingabedateien für einen Workflow bereitzustellen, ist dies mit einer Textdatei zu tun, die die Dateipfade enthält.
Es kann so einfach sein wie eine Textdatei, die einen Dateipfad pro Zeile auflistet und sonst nichts, oder die Datei kann zusätzliche Metadaten enthalten, in diesem Fall wird sie oft als Samplesheet bezeichnet.

Hier zeigen wir dir, wie du den einfachen Fall machst.

### 4.1. Die bereitgestellte Textdatei mit den Eingabedateipfaden untersuchen

Wir haben bereits eine Textdatei erstellt, die die Eingabedateipfade auflistet, genannt `sample_bams.txt`, die du im Verzeichnis `data/` findest.

```txt title="sample_bams.txt"
/workspaces/training/nf4-science/genomics/data/bam/reads_mother.bam
/workspaces/training/nf4-science/genomics/data/bam/reads_father.bam
/workspaces/training/nf4-science/genomics/data/bam/reads_son.bam
```

Wie du sehen kannst, haben wir einen Dateipfad pro Zeile aufgelistet, und es sind absolute Pfade.

!!! note "Hinweis"

    Die Dateien, die wir hier verwenden, befinden sich nur auf dem lokalen Dateisystem deines GitHub Codespaces, aber wir könnten auch auf Dateien im Cloud-Speicher zeigen.
    Wenn du nicht die bereitgestellte Codespaces-Umgebung verwendest, musst du möglicherweise die Dateipfade an dein lokales Setup anpassen.

### 4.2. Den Parameter und das Testprofil aktualisieren

Ändere den Parameter `reads_bam` so, dass er auf die Datei `sample_bams.txt` zeigt, anstatt einzelne Proben aufzulisten.

Stelle die Typannotation im params-Block wieder her (da es wieder ein einzelner Pfad ist):

=== "Danach"

    ```groovy title="genomics.nf" linenums="10" hl_lines="1-2"
        // Primäre Eingabe (Datei mit Eingabedateien, eine pro Zeile)
        reads_bam: Path
    ```

=== "Vorher"

    ```groovy title="genomics.nf" linenums="10"
        // Primäre Eingabe (Array von drei Proben)
        reads_bam
    ```

Aktualisiere dann das Testprofil, um auf die Textdatei zu zeigen:

=== "Danach"

    ```groovy title="nextflow.config" linenums="4" hl_lines="2"
    test {
        params.reads_bam = "${projectDir}/data/sample_bams.txt"
        params.reference = "${projectDir}/data/ref/ref.fasta"
        params.reference_index = "${projectDir}/data/ref/ref.fasta.fai"
        params.reference_dict = "${projectDir}/data/ref/ref.dict"
        params.intervals = "${projectDir}/data/ref/intervals.bed"
    }
    ```

=== "Vorher"

    ```groovy title="nextflow.config" linenums="4"
    test {
        params.reads_bam = [
            "${projectDir}/data/bam/reads_mother.bam",
            "${projectDir}/data/bam/reads_father.bam",
            "${projectDir}/data/bam/reads_son.bam"
        ]
        params.reference = "${projectDir}/data/ref/ref.fasta"
        params.reference_index = "${projectDir}/data/ref/ref.fasta.fai"
        params.reference_dict = "${projectDir}/data/ref/ref.dict"
        params.intervals = "${projectDir}/data/ref/intervals.bed"
    }
    ```

Die Liste der Dateien lebt nicht mehr im Code, was ein großer Schritt in die richtige Richtung ist.

### 4.3. Die Kanal-Factory aktualisieren, um Zeilen aus einer Datei zu lesen

Derzeit behandelt unsere Eingabekanal-Factory alle Dateien, die wir ihr geben, als die Dateneingaben, die wir dem Indizierungsprozess zuführen möchten.
Da wir ihr jetzt eine Datei geben, die Eingabedateipfade auflistet, müssen wir ihr Verhalten ändern, um die Datei zu parsen und die darin enthaltenen Dateipfade als Dateneingaben zu behandeln.

Wir können dies mit demselben Muster tun, das wir in [Teil 2 von Hello Nextflow](../../hello_nextflow/02_hello_channels.md#42-use-the-splitcsv-operator-to-parse-the-file) verwendet haben: Anwenden des Operators [`splitCsv()`](https://nextflow.io/docs/latest/reference/operator.html#splitcsv), um die Datei zu parsen, dann eine `map`-Operation, um das erste Feld jeder Zeile auszuwählen.

=== "Danach"

    ```groovy title="genomics.nf" linenums="24" hl_lines="1-4"
        // Eingabekanal aus einer CSV-Datei erstellen, die Eingabedateipfade auflistet
        reads_ch = Channel.fromPath(params.reads_bam)
                .splitCsv()
                .map { line -> file(line[0]) }
    ```

=== "Vorher"

    ```groovy title="genomics.nf" linenums="24"
        // Eingabekanal erstellen (einzelne Datei über CLI-Parameter)
        reads_ch = channel.fromPath(params.reads_bam)
    ```

Technisch könnten wir dies einfacher mit dem Operator [`.splitText()`](https://www.nextflow.io/docs/latest/reference/operator.html#operator-splittext) tun, da unsere Eingabedatei derzeit nur Dateipfade enthält.
Durch die Verwendung des vielseitigeren `splitCsv`-Operators (ergänzt durch `map`) können wir unseren Workflow jedoch zukunftssicher machen, falls wir uns entscheiden, Metadaten zur Datei mit Dateipfaden hinzuzufügen.

!!! tip "Tipp"

    Wenn du nicht sicher bist, dass du verstehst, was die Operatoren hier tun, ist dies eine weitere großartige Gelegenheit, den `.view()`-Operator zu verwenden, um zu sehen, wie die Kanalinhalte vor und nach der Anwendung aussehen.

### 4.4. Den Workflow ausführen

Führe den Workflow noch einmal aus. Dies sollte dasselbe Ergebnis wie zuvor liefern, oder?

```bash
nextflow run genomics.nf -profile test -resume
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics.nf` [sick_albattani] DSL2 - revision: 46d84642f6

    [18/23b4bb] SAMTOOLS_INDEX (1)       | 3 of 3, cached: 3 ✔
    [12/f727bb] GATK_HAPLOTYPECALLER (3) | 3 of 3, cached: 3 ✔
    ```

Ja! Tatsächlich erkennt Nextflow korrekt, dass die Prozessaufrufe genau gleich sind, und macht sich nicht einmal die Mühe, alles erneut auszuführen, da wir mit `-resume` liefen.

Und das war's! Unser einfacher Variantenerkennungs-Workflow hat alle grundlegenden Funktionen, die wir wollten.

### Fazit

Du weißt, wie man einen mehrstufigen modularen Workflow erstellt, um eine BAM-Datei zu indizieren und Variantenerkennung pro Probe mit GATK anzuwenden.

Allgemeiner hast du gelernt, wie man wesentliche Nextflow-Komponenten und -Logik verwendet, um eine einfache Genomik-Pipeline zu erstellen, die echte Arbeit leistet und die Eigenheiten von Genomik-Dateiformaten und Tool-Anforderungen berücksichtigt.

### Wie geht es weiter?

Feiere deinen Erfolg und mache eine extra lange Pause!

Im nächsten Teil dieses Kurses lernst du, wie du diesen einfachen Variantenerkennungs-Workflow pro Probe transformierst, um gemeinsame Variantenerkennung auf die Daten anzuwenden.
