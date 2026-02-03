# Teil 1: Variantenerkennung pro Probe

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Im ersten Teil dieses Kurses zeigen wir dir, wie du eine einfache Varianten-Pipeline erstellst, die GATK-Variantenerkennung auf einzelne Sequenzierungsproben anwendet.

### Methodenübersicht

Variantenerkennung ist eine genomische Analysemethode, die darauf abzielt, Variationen in einer Genomsequenz relativ zu einem Referenzgenom zu identifizieren.
Hier werden wir Tools und Methoden verwenden, die für die Erkennung kurzer Varianten entwickelt wurden, _d.h._ SNPs und Indels.

![GATK pipeline](img/gatk-pipeline.png)

Eine vollständige Varianten-Pipeline umfasst typischerweise viele Schritte, einschließlich Mapping zum Referenzgenom (manchmal als Genom-Alignment bezeichnet) sowie Variantenfilterung und -priorisierung.
Der Einfachheit halber konzentrieren wir uns in diesem Teil des Kurses nur auf den Variantenerkennungsteil.

### Datensatz

Wir stellen die folgenden Daten und zugehörige Ressourcen bereit:

- **Ein Referenzgenom**, das aus einer kleinen Region des menschlichen Chromosoms 20 (aus hg19/b37) und seinen Zusatzdateien (Index und Sequenzwörterbuch) besteht.
- **Drei Gesamtgenom-Sequenzierungsproben**, die einem Familientrio (Mutter, Vater und Sohn) entsprechen, die auf einen kleinen Datenausschnitt auf Chromosom 20 reduziert wurden, um die Dateigrößen klein zu halten.
  Dies sind Illumina-Short-Read-Sequenzierungsdaten, die bereits auf das Referenzgenom gemappt wurden und im [BAM](https://samtools.github.io/hts-specs/SAMv1.pdf)-Format (Binary Alignment Map, eine komprimierte Version von SAM, Sequence Alignment Map) bereitgestellt werden.
- **Eine Liste genomischer Intervalle**, d.h. Koordinaten auf dem Genom, wo unsere Proben Daten haben, die für die Variantenerkennung geeignet sind, bereitgestellt im BED-Format.

### Workflow

In diesem Teil des Kurses werden wir einen Workflow entwickeln, der Folgendes tut:

1. Generiere eine Indexdatei für jede BAM-Eingabedatei mit [Samtools](https://www.htslib.org/)
2. Führe GATK HaplotypeCaller auf jeder BAM-Eingabedatei aus, um Variantenaufrufe pro Probe im VCF-Format (Variant Call Format) zu generieren

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/genomics/img/hello-gatk-1.svg"
</figure>

!!! note "Hinweis"

    Indexdateien sind ein häufiges Merkmal bioinformatischer Dateiformate; sie enthalten Informationen über die Struktur der Hauptdatei, die es Tools wie GATK ermöglichen, auf eine Teilmenge der Daten zuzugreifen, ohne die gesamte Datei lesen zu müssen.
    Dies ist wichtig, weil diese Dateien sehr groß werden können.

---

## 0. Aufwärmen: Teste die Samtools- und GATK-Befehle interaktiv

Zuerst möchten wir die Befehle manuell ausprobieren, bevor wir versuchen, sie in einen Workflow zu verpacken.
Die Tools, die wir benötigen (Samtools und GATK), sind nicht in der GitHub Codespaces-Umgebung installiert, also werden wir sie über Container verwenden (siehe [Hello Containers](../../hello_nextflow/05_hello_containers.md)).

!!! note "Hinweis"

     Stelle sicher, dass du im Verzeichnis `nf4-science/genomics` bist, sodass der letzte Teil des Pfads, der angezeigt wird, wenn du `pwd` eingibst, `genomics` ist.

### 0.1. Indexiere eine BAM-Eingabedatei mit Samtools

Wir werden einen Samtools-Container herunterladen, ihn interaktiv starten und den Befehl `samtools index` auf eine der BAM-Dateien anwenden.

#### 0.1.1. Lade den Samtools-Container herunter

```bash
docker pull community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464
```

<!--
??? success "Befehlsausgabe"

    ```console

    ```
-->

#### 0.1.2. Starte den Samtools-Container interaktiv

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464
```

<!--
??? success "Befehlsausgabe"

    ```console

    ```
-->

#### 0.1.3. Führe den Indexierungsbefehl aus

Die [Samtools-Dokumentation](https://www.htslib.org/doc/samtools-index.html) gibt uns die Befehlszeile, um eine BAM-Datei zu indexieren.

Wir müssen nur die Eingabedatei angeben; das Tool generiert automatisch einen Namen für die Ausgabe, indem es `.bai` an den Eingabedateinamen anhängt.

```bash
samtools index /data/bam/reads_mother.bam
```

Dies sollte sofort abgeschlossen sein, und du solltest jetzt eine Datei namens `reads_mother.bam.bai` im selben Verzeichnis wie die ursprüngliche BAM-Eingabedatei sehen.

??? abstract "Verzeichnisinhalt"

    ```console
    data/bam/
    ├── reads_father.bam
    ├── reads_mother.bam
    ├── reads_mother.bam.bai
    └── reads_son.bam
    ```

#### 0.1.4. Verlasse den Samtools-Container

```bash
exit
```

### 0.2. Rufe Varianten mit GATK HaplotypeCaller auf

Wir werden einen GATK-Container herunterladen, ihn interaktiv starten und den Befehl `gatk HaplotypeCaller` auf die BAM-Datei anwenden, die wir gerade indexiert haben.

#### 0.2.1. Lade den GATK-Container herunter

```bash
docker pull community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867
```

<!--
??? success "Befehlsausgabe"

    ```console

    ```
-->

#### 0.2.2. Starte den GATK-Container interaktiv

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867
```

<!--
??? success "Befehlsausgabe"

    ```console

    ```
-->

#### 0.2.3. Führe den Variantenerkennungsbefehl aus

Die [GATK-Dokumentation](https://gatk.broadinstitute.org/hc/en-us/articles/21905025322523-HaplotypeCaller) gibt uns die Befehlszeile, um Variantenerkennung auf einer BAM-Datei durchzuführen.

Wir müssen die BAM-Eingabedatei (`-I`) sowie das Referenzgenom (`-R`), einen Namen für die Ausgabedatei (`-O`) und eine Liste genomischer Intervalle zur Analyse (`-L`) angeben.

Wir müssen jedoch nicht den Pfad zur Indexdatei angeben; das Tool sucht automatisch danach im selben Verzeichnis, basierend auf der etablierten Namens- und Platzierungskonvention.
Dasselbe gilt für die Zusatzdateien des Referenzgenoms (Index- und Sequenzwörterbuchdateien, `*.fai` und `*.dict`).

```bash
gatk HaplotypeCaller \
        -R /data/ref/ref.fasta \
        -I /data/bam/reads_mother.bam \
        -O reads_mother.vcf \
        -L /data/ref/intervals.bed
```

<!--
??? success "Befehlsausgabe"

    ```console

    ```
-->

Die Ausgabedatei `reads_mother.vcf` wird in deinem Arbeitsverzeichnis im Container erstellt, sodass du sie im VS Code-Datei-Explorer nicht sehen wirst, es sei denn, du änderst den Ausgabedateipfad.
Es ist jedoch eine kleine Testdatei, sodass du sie mit `cat` öffnen und den Inhalt anzeigen kannst.
Wenn du ganz nach oben zum Anfang der Datei scrollst, findest du einen Header, der aus vielen Zeilen mit Metadaten besteht, gefolgt von einer Liste von Variantenaufrufen, einer pro Zeile.

```console title="reads_mother.vcf" linenums="26"
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	reads_mother
20_10037292_10066351	3480	.	C	CT	503.03	.	AC=2;AF=1.00;AN=2;DP=23;ExcessHet=0.0000;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=27.95;SOR=1.179	GT:AD:DP:GQ:PL	1/1:0,18:18:54:517,54,0
20_10037292_10066351	3520	.	AT	A	609.03	.	AC=2;AF=1.00;AN=2;DP=18;ExcessHet=0.0000;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=33.83;SOR=0.693	GT:AD:DP:GQ:PL	1/1:0,18:18:54:623,54,0
20_10037292_10066351	3529	.	T	A	155.64	.	AC=1;AF=0.500;AN=2;BaseQRankSum=-0.544;DP=21;ExcessHet=0.0000;FS=1.871;MLEAC=1;MLEAF=0.500;MQ=60.00;MQRankSum=0.000;QD=7.78;ReadPosRankSum=-1.158;SOR=1.034	GT:AD:DP:GQ:PL	0/1:12,8:20:99:163,0,328
```

Jede Zeile beschreibt eine mögliche Variante, die in den Sequenzierungsdaten der Probe identifiziert wurde. Für Anleitungen zur Interpretation des VCF-Formats siehe [diesen hilfreichen Artikel](https://www.ebi.ac.uk/training/online/courses/human-genetic-variation-introduction/variant-identification-and-analysis/understanding-vcf-format/).

Die Ausgabe-VCF-Datei wird von einer Indexdatei namens `reads_mother.vcf.idx` begleitet, die automatisch von GATK erstellt wurde.
Sie hat dieselbe Funktion wie die BAM-Indexdatei, um Tools zu ermöglichen, Teilmengen von Daten zu suchen und abzurufen, ohne die gesamte Datei zu laden.

#### 0.2.4. Verlasse den GATK-Container

```bash
exit
```

### Zusammenfassung

Du weißt, wie du die Samtools-Indexierung und GATK-Variantenerkennungsbefehle in ihren jeweiligen Containern testest.

### Was kommt als Nächstes?

Lerne, wie du diese Befehle in einen zweistufigen Workflow verpackst, der Container zur Ausführung der Arbeit verwendet.

---

## 1. Schreibe einen einstufigen Workflow, der Samtools index auf einer BAM-Datei ausführt

Wir stellen dir eine Workflow-Datei, `genomics-1.nf`, zur Verfügung, die die Hauptteile des Workflows umreißt.
Sie ist nicht funktionsfähig; ihr Zweck ist es nur, als Skelett zu dienen, das du zum Schreiben des tatsächlichen Workflows verwenden wirst.

### 1.1. Definiere den Indexierungsprozess

Beginnen wir mit dem Schreiben eines Prozesses, den wir `SAMTOOLS_INDEX` nennen werden, der die Indexierungsoperation beschreibt.

```groovy title="genomics-1.nf" linenums="9"
/*
 * Generiere BAM-Indexdatei
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

Du solltest alle Teile aus dem wiedererkennen, was du in Teil 1 & Teil 2 dieser Trainingsreihe gelernt hast.

Dieser Prozess wird von uns verlangen, einen Dateipfad über die Eingabe `input_bam` zu übergeben, also richten wir das als Nächstes ein.

### 1.2. Füge eine Eingabeparameterdeklaration hinzu

Am Anfang der Datei, unter dem Abschnitt `Pipeline parameters`, deklarieren wir einen CLI-Parameter namens `reads_bam` und geben ihm einen Standardwert.
Auf diese Weise können wir faul sein und die Eingabe nicht angeben, wenn wir den Befehl zum Starten der Pipeline eingeben (für Entwicklungszwecke).

```groovy title="genomics-1.nf" linenums="3"
/*
 * Pipeline parameters
 */
params {
    // Primäre Eingabe
    reads_bam: Path = "${projectDir}/data/bam/reads_mother.bam"
}
```

Jetzt haben wir einen Prozess bereit, sowie einen Parameter, um ihm eine Eingabe zum Ausführen zu geben, also verbinden wir diese Dinge zusammen.

!!! note "Hinweis"

    `${projectDir}` ist eine eingebaute Nextflow-Variable, die auf das Verzeichnis zeigt, in dem sich das aktuelle Nextflow-Workflow-Skript (`genomics-1.nf`) befindet.

    Dies macht es einfach, auf Dateien, Datenverzeichnisse und andere Ressourcen zu verweisen, die im Workflow-Repository enthalten sind, ohne absolute Pfade fest zu codieren.

### 1.3. Füge einen workflow-Block hinzu, um SAMTOOLS_INDEX auszuführen

Im `workflow`-Block müssen wir einen **channel** einrichten, um die Eingabe an den `SAMTOOLS_INDEX`-Prozess zu übergeben; dann können wir den Prozess selbst aufrufen, um auf den Inhalten dieses Channels zu laufen.

```groovy title="genomics-1.nf" linenums="24"
workflow {

    main:
    // Erstelle Eingabekanal (einzelne Datei über CLI-Parameter)
    reads_ch = channel.fromPath(params.reads_bam)

    // Erstelle Indexdatei für Eingabe-BAM-Datei
    SAMTOOLS_INDEX(reads_ch)

    publish:
    bam_index = SAMTOOLS_INDEX.out
}
```

Der workflow-Block hat zwei Abschnitte:

- `main:` enthält die Channel-Operationen und Prozessaufrufe
- `publish:` deklariert, welche Ausgaben veröffentlicht werden sollen, und weist sie benannten Zielen zu

Du wirst bemerken, dass wir dieselbe `.fromPath`-Channel-Factory verwenden, wie wir sie in [Hello Channels](../../hello_nextflow/02_hello_channels.md) verwendet haben.
Tatsächlich machen wir etwas sehr Ähnliches.
Der Unterschied ist, dass wir Nextflow sagen, nur den Dateipfad selbst in den Channel als Eingabeelement zu laden, anstatt seinen Inhalt einzulesen.

### 1.4. Füge einen output-Block hinzu, um zu definieren, wo Ergebnisse veröffentlicht werden

Nach dem workflow-Block fügen wir einen `output`-Block hinzu, der angibt, wo die Workflow-Ausgaben veröffentlicht werden sollen.

```groovy title="genomics-1.nf" linenums="37"
output {
    bam_index {
        path '.'
    }
}
```

Jedes benannte Ziel aus dem `publish:`-Abschnitt (wie `bam_index`) erhält seinen eigenen Block, in dem du den Ausgabepfad relativ zum Basis-Ausgabeverzeichnis konfigurieren kannst.

!!! note "Hinweis"

    Obwohl die Datendateien, die wir hier verwenden, sehr klein sind, können sie in der Genomik sehr groß werden.
    Standardmäßig erstellt Nextflow symbolische Links zu den Ausgabedateien im Veröffentlichungsverzeichnis, was unnötige Dateikopien vermeidet.
    Du kannst dieses Verhalten mit der Option `mode` ändern (z.B. `mode 'copy'`), um stattdessen tatsächliche Kopien zu erstellen.
    Beachte, dass Symlinks kaputt gehen, wenn du dein `work`-Verzeichnis aufräumst, sodass du für Produktions-Workflows möglicherweise `mode 'copy'` verwenden möchtest.

### 1.5. Konfiguriere das Ausgabeverzeichnis

Das Basis-Ausgabeverzeichnis wird über die Konfigurationsoption `outputDir` festgelegt. Füge sie zu `nextflow.config` hinzu:

=== "Danach"

    ```groovy title="nextflow.config" hl_lines="2"
    docker.enabled = true
    outputDir = 'results_genomics'
    ```

=== "Vorher"

    ```groovy title="nextflow.config"
    docker.enabled = true
    ```

### 1.6. Führe den Workflow aus, um zu überprüfen, dass der Indexierungsschritt funktioniert

Lass uns den Workflow ausführen! Zur Erinnerung: Wir müssen keine Eingabe in der Befehlszeile angeben, weil wir einen Standardwert für die Eingabe festgelegt haben, als wir den Eingabeparameter deklariert haben.

```bash
nextflow run genomics-1.nf
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics-1.nf` [reverent_sinoussi] DSL2 - revision: 41d43ad7fe

    executor >  local (1)
    [2a/e69536] SAMTOOLS_INDEX (1) | 1 of 1 ✔
    ```

Du kannst überprüfen, dass die Indexdatei korrekt generiert wurde, indem du im work-Verzeichnis oder im results-Verzeichnis nachsiehst.

??? abstract "Work-Verzeichnisinhalt"

    ```console
    work/2a/e695367b2f60df09cf826b07192dc3
    ├── reads_mother.bam -> /workspaces/training/nf4-science/genomics/data/bam/reads_mother.bam
    └── reads_mother.bam.bai
    ```

??? abstract "Results-Verzeichnisinhalt"

    ```console
    results_genomics/
    └── reads_mother.bam.bai -> */87/908bba*/reads_mother.bam.bai
    ```

Da ist sie!

### Zusammenfassung

Du weißt, wie du ein Genomik-Tool in einen einstufigen Nextflow-Workflow verpackst und es mit einem Container ausführen lässt.

### Was kommt als Nächstes?

Füge einen zweiten Schritt hinzu, der die Ausgabe des ersten verwendet.

---

## 2. Füge einen zweiten Prozess hinzu, um GATK HaplotypeCaller auf der indexierten BAM-Datei auszuführen

Jetzt, da wir einen Index für unsere Eingabedatei haben, können wir mit dem Einrichten des Variantenerkennungsschritts fortfahren, der der interessante Teil des Workflows ist.

### 2.1. Definiere den Variantenerkennungsprozess

Lass uns einen Prozess schreiben, den wir `GATK_HAPLOTYPECALLER` nennen werden, der die Variantenerkennungsoperation beschreibt.

```groovy title="genomics-1.nf" linenums="44"
/*
 * Rufe Varianten mit GATK HaplotypeCaller auf
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

Du wirst bemerken, dass wir hier eine neue Syntax (`emit:`) eingeführt haben, um jeden unserer Ausgabekanäle eindeutig zu benennen, und die Gründe dafür werden bald klar werden.

Dieser Befehl benötigt einige weitere Eingaben, weil GATK mehr Informationen benötigt, um die Analyse durchzuführen, verglichen mit einem einfachen Indexierungsjob.
Aber du wirst feststellen, dass im Eingabeblock noch mehr Eingaben definiert sind, als im GATK-Befehl aufgelistet sind. Warum ist das so?

!!! note "Hinweis"

    GATK weiß, dass es nach der BAM-Indexdatei und den Zusatzdateien des Referenzgenoms suchen muss, weil es sich der Konventionen bewusst ist, die diese Dateien umgeben.
    Nextflow ist jedoch so konzipiert, dass es domänenunabhängig ist und nichts über bioinformatische Dateiformatanforderungen weiß.

Wir müssen Nextflow explizit mitteilen, dass es diese Dateien zur Laufzeit im Arbeitsverzeichnis bereitstellen muss; andernfalls wird es das nicht tun, und GATK wird (zu Recht) einen Fehler werfen, dass die Indexdateien fehlen.

Ebenso müssen wir die Indexdatei der Ausgabe-VCF (die Datei `"${input_bam}.vcf.idx"`) explizit auflisten, damit Nextflow weiß, dass es diese Datei im Auge behalten muss, falls sie in nachfolgenden Schritten benötigt wird.

### 2.2. Füge Definitionen für Zusatzeingaben hinzu

Da unser neuer Prozess erwartet, dass einige zusätzliche Dateien bereitgestellt werden, richten wir einige CLI-Parameter für sie unter dem Abschnitt `Pipeline parameters` ein, zusammen mit einigen Standardwerten (aus denselben Gründen wie zuvor).

```groovy title="genomics-1.nf" linenums="8"
    // Zusatzdateien
    reference: Path = "${projectDir}/data/ref/ref.fasta"
    reference_index: Path = "${projectDir}/data/ref/ref.fasta.fai"
    reference_dict: Path = "${projectDir}/data/ref/ref.dict"
    intervals: Path = "${projectDir}/data/ref/intervals.bed"
```

### 2.3. Erstelle Variablen, um die Zusatzdateipfade zu halten

Während Hauptdateneingaben dynamisch durch Channels geleitet werden, gibt es zwei Ansätze für die Handhabung von Zusatzdateien. Der empfohlene Ansatz ist das Erstellen expliziter Channels, was den Datenfluss klarer und konsistenter macht. Alternativ kann die Funktion file() zum Erstellen von Variablen für einfachere Fälle verwendet werden, insbesondere wenn du dieselbe Datei in mehreren Prozessen referenzieren musst - beachte jedoch, dass dies immer noch implizit Channels erstellt. <!-- TODO: Klarstellen: Ist das immer noch notwendig bei typisierten Eingaben? -->

Füge dies zum workflow-Block hinzu (nach der Erstellung von `reads_ch`, innerhalb des `main:`-Abschnitts):

```groovy title="genomics-1.nf" linenums="79"
    // Lade die Dateipfade für die Zusatzdateien (Referenz und Intervalle)
    ref_file        = file(params.reference)
    ref_index_file  = file(params.reference_index)
    ref_dict_file   = file(params.reference_dict)
    intervals_file  = file(params.intervals)
```

Dies macht die Zusatzdateipfade verfügbar, um sie als Eingabe für alle Prozesse bereitzustellen, die sie benötigen.

### 2.4. Füge einen Aufruf zum workflow-Block hinzu, um GATK_HAPLOTYPECALLER auszuführen

Jetzt, da wir unseren zweiten Prozess eingerichtet haben und alle Eingaben und Zusatzdateien bereit und verfügbar sind, können wir einen Aufruf zum `GATK_HAPLOTYPECALLER`-Prozess im workflow-Body hinzufügen.

```groovy title="genomics-1.nf" linenums="88"
    // Rufe Varianten aus der indexierten BAM-Datei auf
    GATK_HAPLOTYPECALLER(
        reads_ch,
        SAMTOOLS_INDEX.out,
        ref_file,
        ref_index_file,
        ref_dict_file,
        intervals_file
    )
```

Du solltest die `*.out`-Syntax aus Teil 1 dieser Trainingsreihe wiedererkennen; wir sagen Nextflow, dass es die von `SAMTOOLS_INDEX` ausgegebene Channel-Ausgabe nimmt und diese in den `GATK_HAPLOTYPECALLER`-Prozessaufruf einsteckt.

!!! note "Hinweis"

    Du wirst bemerken, dass die Eingaben im Aufruf des Prozesses in genau derselben Reihenfolge bereitgestellt werden, in der sie im Eingabeblock des Prozesses aufgelistet sind.
    In Nextflow sind Eingaben positionsabhängig, was bedeutet, dass du _derselben Reihenfolge folgen musst_; und natürlich muss es die gleiche Anzahl von Elementen geben.

### 2.5. Aktualisiere den publish-Abschnitt und den output-Block

Wir müssen den `publish:`-Abschnitt aktualisieren, um die VCF-Ausgaben einzuschließen, und entsprechende Ziele im `output`-Block hinzufügen.

```groovy title="genomics-1.nf" linenums="99"
    publish:
    bam_index = SAMTOOLS_INDEX.out
    vcf = GATK_HAPLOTYPECALLER.out.vcf
    vcf_idx = GATK_HAPLOTYPECALLER.out.idx
}

output {
    bam_index {
        path '.'
    }
    vcf {
        path '.'
    }
    vcf_idx {
        path '.'
    }
}
```

### 2.6. Führe den Workflow aus, um zu überprüfen, dass der Variantenerkennungsschritt funktioniert

Lass uns den erweiterten Workflow mit `-resume` ausführen, sodass wir den Indexierungsschritt nicht erneut ausführen müssen.

```bash
nextflow run genomics-1.nf -resume
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics-1.nf` [grave_volta] DSL2 - revision: 4790abc96a

    executor >  local (1)
    [2a/e69536] SAMTOOLS_INDEX (1)       | 1 of 1, cached: 1 ✔
    [53/e18e98] GATK_HAPLOTYPECALLER (1) | 1 of 1 ✔
    ```

Wenn wir jetzt die Konsolenausgabe betrachten, sehen wir die beiden aufgelisteten Prozesse.

Der erste Prozess wurde dank des Cachings übersprungen, wie erwartet, während der zweite Prozess ausgeführt wurde, da er brandneu ist.

Du findest die Ausgabedateien im results-Verzeichnis (als symbolische Links zum work-Verzeichnis).

??? abstract "Verzeichnisinhalt"

    ```console
    results_genomics/
    ├── reads_mother.bam.bai -> */87/908bba*/reads_mother.bam.bai
    ├── reads_mother.bam.vcf -> */cf/36f756*/reads_mother.bam.vcf
    └── reads_mother.bam.vcf.idx -> */cf/36f756*/reads_mother.bam.vcf.idx
    ```

Wenn du die VCF-Datei öffnest, solltest du denselben Inhalt sehen wie in der Datei, die du durch direktes Ausführen des GATK-Befehls im Container generiert hast.

```console title="reads_mother.bam.vcf" linenums="26"
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	reads_mother
20_10037292_10066351	3480	.	C	CT	503.03	.	AC=2;AF=1.00;AN=2;DP=23;ExcessHet=0.0000;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=27.95;SOR=1.179	GT:AD:DP:GQ:PL	1/1:0,18:18:54:517,54,0
20_10037292_10066351	3520	.	AT	A	609.03	.	AC=2;AF=1.00;AN=2;DP=18;ExcessHet=0.0000;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=33.83;SOR=0.693	GT:AD:DP:GQ:PL	1/1:0,18:18:54:623,54,0
20_10037292_10066351	3529	.	T	A	155.64	.	AC=1;AF=0.500;AN=2;BaseQRankSum=-0.544;DP=21;ExcessHet=0.0000;FS=1.871;MLEAC=1;MLEAF=0.500;MQ=60.00;MQRankSum=0.000;QD=7.78;ReadPosRankSum=-1.158;SOR=1.034	GT:AD:DP:GQ:PL	0/1:12,8:20:99:163,0,328
```

Dies ist die Ausgabe, die uns wichtig ist, für jede Probe in unserer Studie zu generieren.

### Zusammenfassung

Du weißt, wie man einen sehr einfachen zweistufigen Workflow erstellt, der echte Analysearbeit leistet und in der Lage ist, mit Eigenheiten genomischer Dateiformate wie den Zusatzdateien umzugehen.

### Was kommt als Nächstes?

Mache den Workflow so, dass er mehrere Proben in großen Mengen verarbeitet.

---

## 3. Passe den Workflow an, um auf mehreren Proben zu laufen

Es ist schön und gut, einen Workflow zu haben, der die Verarbeitung einer einzelnen Probe automatisieren kann, aber was ist, wenn du 1000 Proben hast?
Musst du ein Bash-Skript schreiben, das durch alle deine Proben läuft?

Nein, zum Glück! Mache einfach eine kleine Anpassung am Code und Nextflow wird das auch für dich erledigen.

### 3.1. Verwandle die Eingabeparameterdeklaration in ein Array, das die drei Proben auflistet

Lass uns diesen Standarddateipfad in der Eingabe-BAM-Dateideklaration in ein Array verwandeln, das Dateipfade für unsere drei Testproben auflistet, oben unter dem Abschnitt `Pipeline parameters`.

=== "Danach"

    ```groovy title="genomics-1.nf" linenums="7"
    // Primäre Eingabe (Array von drei Proben)
        reads_bam = [
            "${projectDir}/data/bam/reads_mother.bam",
            "${projectDir}/data/bam/reads_father.bam",
            "${projectDir}/data/bam/reads_son.bam"
        ]
    ```

=== "Vorher"

    ```groovy title="genomics-1.nf" linenums="7"
        // Primäre Eingabe
        reads_bam: Path = "${projectDir}/data/bam/reads_mother.bam"
    ```

!!! note "Hinweis"

    Wenn du typisierte Parameterdeklarationen verwendest (wie `reads_bam: Path`), kannst du keinen Array-Wert zuweisen.
    Für Arrays lasse die Typannotation weg.

Und das ist tatsächlich alles, was wir tun müssen, weil die Channel-Factory, die wir im workflow-Body verwenden (`.fromPath`), genauso gerne mehrere Dateipfade akzeptiert, um sie in den Eingabekanal zu laden, wie sie einen einzelnen geladen hat.

!!! note "Hinweis"

    Normalerweise würdest du die Liste der Proben nicht in deine Workflow-Datei fest codieren wollen, aber wir tun das hier, um die Dinge einfach zu halten.
    Wir werden später in dieser Trainingsreihe elegantere Möglichkeiten für die Handhabung von Eingaben präsentieren.

### 3.2. Führe den Workflow aus, um zu überprüfen, dass er auf allen drei Proben läuft

Lass uns versuchen, den Workflow jetzt auszuführen, da die Verkabelung so eingerichtet ist, dass sie auf allen drei Testproben läuft.

```bash
nextflow run genomics-1.nf -resume
```

Lustige Sache: Dies _könnte funktionieren_, ODER es _könnte fehlschlagen_. Hier ist zum Beispiel eine Ausführung, die erfolgreich war:

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics-1.nf` [peaceful_yalow] DSL2 - revision: a256d113ad

    executor >  local (6)
    [4f/7071b0] SAMTOOLS_INDEX (3)       | 3 of 3, cached: 1 ✔
    [7a/89bc43] GATK_HAPLOTYPECALLER (2) | 3 of 3, cached: 1 ✔
    ```

Wenn deine Workflow-Ausführung erfolgreich war, führe sie erneut aus, bis du einen Fehler wie diesen erhältst:

??? failure "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics-1.nf` [loving_pasteur] DSL2 - revision: d2a8e63076

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

Wenn du die GATK-Befehlsfehlerausgabe betrachtest, wird es eine Zeile wie diese geben:

```console
A USER ERROR has occurred: Traversal by intervals was requested but some input files are not indexed.
```

Nun, das ist seltsam, wenn man bedenkt, dass wir die BAM-Dateien explizit im ersten Schritt des Workflows indexiert haben. Könnte es ein Problem mit der Verkabelung geben?

#### 3.2.1. Überprüfe die work-Verzeichnisse für die relevanten Aufrufe

Lass uns einen Blick in das work-Verzeichnis für den fehlgeschlagenen `GATK_HAPLOTYPECALLER`-Prozessaufruf werfen, der in der Konsolenausgabe aufgelistet ist.

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

Achte besonders auf die Namen der BAM-Datei und des BAM-Index, die in diesem Verzeichnis aufgelistet sind: `reads_son.bam` und `reads_father.bam.bai`.

Was zum Teufel? Nextflow hat eine Indexdatei in diesem Prozessaufruf-work-Verzeichnis bereitgestellt, aber es ist die falsche. Wie konnte das passieren?

#### 3.2.2. Verwende den [view()-Operator](https://www.nextflow.io/docs/latest/reference/operator.html#view), um Channel-Inhalte zu inspizieren

Füge diese beiden Zeilen im workflow-Body vor dem `GATK_HAPLOTYPER`-Prozessaufruf hinzu:

```groovy title="genomics-1.nf" linenums="84"
    // temporäre Diagnose
    reads_ch.view()
    SAMTOOLS_INDEX.out.view()
```

Führe dann den Workflow-Befehl erneut aus.

```bash
nextflow run genomics-1.nf
```

Wieder einmal kann dies erfolgreich sein oder fehlschlagen. Hier ist eine erfolgreiche Ausführung:

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics-1.nf` [fervent_pasteur] DSL2 - revision: a256d113ad

    /workspaces/training/nf4-science/genomics/data/bam/reads_mother.bam
    /workspaces/training/nf4-science/genomics/data/bam/reads_father.bam
    /workspaces/training/nf4-science/genomics/data/bam/reads_son.bam
    executor >  local (6)
    [4f/7071b0] SAMTOOLS_INDEX (3)       | 3 of 3 ✔
    /workspaces/training/nf4-science/genomics/work/b4/45a376f0e724be1dc626a6807f73d8/reads_mother.bam.bai
    /workspaces/training/nf4-science/genomics/work/4f/7071b082b45dd85b1c9b6b3b32cb69/reads_father.bam.bai
    /workspaces/training/nf4-science/genomics/work/3c/331645a9e20e67edae10da5ba17c7b/reads_son.bam.bai
    [a2/dbd8d5] GATK_HAPLOTYPECALLER (3) | 3 of 3 ✔
    ```

Und hier ist eine fehlgeschlagene:

??? failure "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics-1.nf` [angry_hamilton] DSL2 - revision: a256d113ad

    /workspaces/training/nf4-science/genomics/data/bam/reads_mother.bam
    /workspaces/training/nf4-science/genomics/data/bam/reads_father.bam
    /workspaces/training/nf4-science/genomics/data/bam/reads_son.bam
    executor >  local (6)
    [4f/7071b0] SAMTOOLS_INDEX (3)       | 3 of 3 ✔
    /workspaces/training/nf4-science/genomics/work/4f/7071b082b45dd85b1c9b6b3b32cb69/reads_father.bam.bai
    /workspaces/training/nf4-science/genomics/work/3c/331645a9e20e67edae10da5ba17c7b/reads_son.bam.bai
    /workspaces/training/nf4-science/genomics/work/b4/45a376f0e724be1dc626a6807f73d8/reads_mother.bam.bai
    [a3/cf3a89] GATK_HAPLOTYPECALLER (3) | 1 of 3
    ERROR ~ Error executing process > 'GATK_HAPLOTYPECALLER (3)'
    ...
    ```

Du musst sie möglicherweise mehrmals ausführen, damit sie erneut fehlschlägt.
Dieser Fehler wird sich nicht konsistent reproduzieren, weil er von einer gewissen Variabilität in den Ausführungszeiten der einzelnen Prozessaufrufe abhängt.

So sieht die Ausgabe der beiden `.view()`-Aufrufe aus, die wir für eine fehlgeschlagene Ausführung hinzugefügt haben:

```console
/workspaces/training/nf4-science/genomics/data/bam/reads_mother.bam
/workspaces/training/nf4-science/genomics/data/bam/reads_father.bam
/workspaces/training/nf4-science/genomics/data/bam/reads_son.bam
/workspaces/training/nf4-science/genomics/work/9c/53492e3518447b75363e1cd951be4b/reads_father.bam.bai
/workspaces/training/nf4-science/genomics/work/cc/37894fffdf6cc84c3b0b47f9b536b7/reads_son.bam.bai
/workspaces/training/nf4-science/genomics/work/4d/dff681a3d137ba7d9866e3d9307bd0/reads_mother.bam.bai
```

Die ersten drei Zeilen entsprechen dem Eingabekanal und die zweite dem Ausgabekanal.
Du kannst sehen, dass die BAM-Dateien und Indexdateien für die drei Proben nicht in derselben Reihenfolge aufgelistet sind!

!!! note "Hinweis"

    Wenn du einen Nextflow-Prozess auf einem Channel aufrufst, der mehrere Elemente enthält, wird Nextflow versuchen, die Ausführung so weit wie möglich zu parallelisieren, und wird Ausgaben in der Reihenfolge sammeln, in der sie verfügbar werden.
    Die Folge ist, dass die entsprechenden Ausgaben in einer anderen Reihenfolge gesammelt werden können als die ursprünglichen Eingaben eingegeben wurden.

Wie derzeit geschrieben, geht unser Workflow-Skript davon aus, dass die Indexdateien aus dem Indexierungsschritt in derselben Mutter/Vater/Sohn-Reihenfolge herauskommen, wie die Eingaben gegeben wurden.
Aber das ist nicht garantiert der Fall, weshalb manchmal (aber nicht immer) die falschen Dateien im zweiten Schritt gepaart werden.

Um dies zu beheben, müssen wir sicherstellen, dass die BAM-Dateien und ihre Indexdateien zusammen durch die Channels reisen.

!!! tip "Tipp"

    Die `view()`-Anweisungen im Workflow-Code tun nichts, sodass es kein Problem ist, sie drin zu lassen.
    Sie werden jedoch deine Konsolenausgabe überladen, daher empfehlen wir, sie zu entfernen, wenn du mit der Fehlerbehebung fertig bist.

### 3.3. Ändere die Ausgabe des SAMTOOLS_INDEX-Prozesses in ein Tupel, das die Eingabedatei und ihren Index zusammenhält

Der einfachste Weg, um sicherzustellen, dass eine BAM-Datei und ihr Index eng verbunden bleiben, besteht darin, sie zusammen in ein Tupel zu verpacken, das aus der Indexaufgabe herauskommt.

!!! note "Hinweis"

    Ein **Tupel** ist eine endliche, geordnete Liste von Elementen, die häufig zum Zurückgeben mehrerer Werte aus einer Funktion verwendet wird. Tupel sind besonders nützlich für die Übergabe mehrerer Ein- oder Ausgaben zwischen Prozessen, während ihre Zuordnung und Reihenfolge erhalten bleibt.

Zuerst ändern wir die Ausgabe des `SAMTOOLS_INDEX`-Prozesses, um die BAM-Datei in ihre Ausgabedeklaration aufzunehmen.

=== "Danach"

    ```groovy title="genomics-1.nf" linenums="20"
        output:
        tuple path(input_bam), path("${input_bam}.bai")
    ```

=== "Vorher"

    ```groovy title="genomics-1.nf" linenums="20"
        output:
        path "${input_bam}.bai"
    ```

Auf diese Weise wird jede Indexdatei eng mit ihrer ursprünglichen BAM-Datei gekoppelt, und die Gesamtausgabe des Indexierungsschritts wird ein einzelner Channel sein, der Dateipaare enthält.

### 3.4. Ändere die Eingabe zum GATK_HAPLOTYPECALLER-Prozess zu einem Tupel

Da wir die 'Form' der Ausgabe des ersten Prozesses im Workflow geändert haben, müssen wir die Eingabedefinition des zweiten Prozesses aktualisieren, damit sie übereinstimmt.

Konkret, wo wir zuvor zwei separate Eingabepfade im Eingabeblock des `GATK_HAPLOTYPECALLER`-Prozesses deklariert haben, deklarieren wir jetzt eine einzelne Eingabe, die der Struktur des von `SAMTOOLS_INDEX` ausgegebenen Tupels entspricht.

=== "Danach"

    ```groovy title="genomics-1.nf" linenums="51"
        input:
        tuple path(input_bam), path(input_bam_index)
    ```

=== "Vorher"

    ```groovy title="genomics-1.nf" linenums="51"
        input:
        path input_bam
        path input_bam_index
    ```

Natürlich müssen wir, da wir jetzt die Form der Eingaben geändert haben, die `GATK_HAPLOTYPECALLER` erwartet, den Prozessaufruf entsprechend im workflow-Body aktualisieren.

### 3.5. Aktualisiere den Aufruf von GATK_HAPLOTYPECALLER im workflow-Block

Wir müssen nicht mehr das ursprüngliche `reads_ch` an den `GATK_HAPLOTYPECALLER`-Prozess übergeben, da die BAM-Datei jetzt in die Channel-Ausgabe von `SAMTOOLS_INDEX` gebündelt ist.

Als Ergebnis können wir diese Zeile einfach löschen.

=== "Danach"

    ```groovy title="genomics-1.nf" linenums="88"
    GATK_HAPLOTYPECALLER(
        SAMTOOLS_INDEX.out,
    ```

=== "Vorher"

    ```groovy title="genomics-1.nf" linenums="88"
    GATK_HAPLOTYPECALLER(
        reads_ch,
        SAMTOOLS_INDEX.out,
    ```

Das ist die gesamte Neu-Verkabelung, die notwendig ist, um das Index-Mismatch-Problem zu lösen.

### 3.6. Aktualisiere den publish-Abschnitt und den output-Block für das Tupel

Da `SAMTOOLS_INDEX.out` jetzt ein Tupel ist, das sowohl die BAM als auch ihren Index enthält, werden beide Dateien zusammen veröffentlicht.
Wir benennen das Ziel von `bam_index` in `indexed_bam` um, um widerzuspiegeln, dass es jetzt beide Dateien enthält.

=== "Danach"

    ```groovy title="genomics-1.nf" hl_lines="2"
        publish:
        indexed_bam = SAMTOOLS_INDEX.out
    ```

=== "Vorher"

    ```groovy title="genomics-1.nf"
        publish:
        bam_index = SAMTOOLS_INDEX.out
    ```

Wir müssen auch den output-Block aktualisieren, um den neuen Zielnamen zu verwenden:

=== "Danach"

    ```groovy title="genomics-1.nf" hl_lines="2"
    output {
        indexed_bam {
            path '.'
        }
    ```

=== "Vorher"

    ```groovy title="genomics-1.nf"
    output {
        bam_index {
            path '.'
        }
    ```

### 3.7. Führe den Workflow aus, um zu überprüfen, dass er jedes Mal korrekt auf allen drei Proben funktioniert

Natürlich liegt der Beweis im Pudding, also lass uns den Workflow einige Male erneut ausführen, um sicherzustellen, dass dies zukünftig zuverlässig funktioniert.

```bash
nextflow run genomics-1.nf
```

Dieses Mal (und jedes Mal) sollte alles korrekt laufen:

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics-1.nf` [special_goldstine] DSL2 - revision: 4cbbf6ea3e

    executor >  local (6)
    [d6/10c2c4] SAMTOOLS_INDEX (1)       | 3 of 3 ✔
    [88/1783aa] GATK_HAPLOTYPECALLER (2) | 3 of 3 ✔
    ```

Das results-Verzeichnis enthält jetzt sowohl BAM- als auch BAI-Dateien für jede Probe (aus dem Tupel), zusammen mit den VCF-Ausgaben:

??? abstract "Results-Verzeichnisinhalt"

    ```console
    results_genomics/
    ├── reads_father.bam -> */60/e2614c*/reads_father.bam
    ├── reads_father.bam.bai -> */60/e2614c*/reads_father.bam.bai
    ├── reads_father.bam.vcf -> */b8/91b3c8*/reads_father.bam.vcf
    ├── reads_father.bam.vcf.idx -> */b8/91b3c8*/reads_father.bam.vcf.idx
    ├── reads_mother.bam -> */3e/fededc*/reads_mother.bam
    ├── reads_mother.bam.bai -> */3e/fededc*/reads_mother.bam.bai
    ├── reads_mother.bam.vcf -> */32/5ca037*/reads_mother.bam.vcf
    ├── reads_mother.bam.vcf.idx -> */32/5ca037*/reads_mother.bam.vcf.idx
    ├── reads_son.bam -> */3c/36d1c2*/reads_son.bam
    ├── reads_son.bam.bai -> */3c/36d1c2*/reads_son.bam.bai
    ├── reads_son.bam.vcf -> */d7/a6b046*/reads_son.bam.vcf
    └── reads_son.bam.vcf.idx -> */d7/a6b046*/reads_son.bam.vcf.idx
    ```

Wenn du möchtest, kannst du `.view()` erneut verwenden, um einen Blick darauf zu werfen, wie die Inhalte des `SAMTOOLS_INDEX`-Ausgabekanals aussehen:

```groovy title="genomics-1.nf" linenums="92"
SAMTOOLS_INDEX.out.view()
```

Du wirst sehen, dass der Channel die drei erwarteten Tupel enthält (Dateipfade zur besseren Lesbarkeit gekürzt).

```console title="Ausgabe"
[*/60/e2614c*/reads_father.bam, */60/e2614c*/reads_father.bam.bai]
[*/3e/fededc*/reads_mother.bam, */3e/fededc*/reads_mother.bam.bai]
[*/3c/36d1c2*/reads_son.bam, */3c/36d1c2*/reads_son.bam.bai]
```

Das wird zukünftig viel sicherer sein.

### Zusammenfassung

Du weißt, wie du deinen Workflow auf mehreren Proben (unabhängig) ausführen lässt.

### Was kommt als Nächstes?

Mache es einfacher, Proben in großen Mengen zu handhaben.

---

## 4. Lasse den Workflow eine Textdatei akzeptieren, die mehrere Eingabedateien enthält

Eine sehr häufige Methode, um mehrere Dateneingabedateien an einen Workflow zu übergeben, besteht darin, dies mit einer Textdatei zu tun, die die Dateipfade enthält.
Es kann so einfach sein wie eine Textdatei, die einen Dateipfad pro Zeile auflistet und sonst nichts, oder die Datei kann zusätzliche Metadaten enthalten, in welchem Fall sie oft als Samplesheet bezeichnet wird.

Hier zeigen wir dir, wie du den einfachen Fall machst.

### 4.1. Untersuche die bereitgestellte Textdatei, die die Eingabedateipfade auflistet

Wir haben bereits eine Textdatei erstellt, die die Eingabedateipfade auflistet, genannt `sample_bams.txt`, die du im Verzeichnis `data/` finden kannst.

```txt title="sample_bams.txt"
/workspaces/training/nf4-science/genomics/data/bam/reads_mother.bam
/workspaces/training/nf4-science/genomics/data/bam/reads_father.bam
/workspaces/training/nf4-science/genomics/data/bam/reads_son.bam
```

Wie du sehen kannst, haben wir einen Dateipfad pro Zeile aufgelistet, und es sind absolute Pfade.

!!! note "Hinweis"

    Die Dateien, die wir hier verwenden, befinden sich nur auf dem lokalen Dateisystem deines GitHub Codespaces, aber wir könnten auch auf Dateien im Cloud-Speicher verweisen.

### 4.2. Aktualisiere den Parameterstandardwert

Lass uns den Standardwert für unseren `reads_bam`-Eingabeparameter ändern, um auf die `sample_bams.txt`-Datei zu zeigen.

=== "Danach"

    ```groovy title="genomics-1.nf" linenums="7"
        // Primäre Eingabe (Datei mit Eingabedateien, eine pro Zeile)
        reads_bam: Path = "${projectDir}/data/sample_bams.txt"
    ```

=== "Vorher"

    ```groovy title="genomics-1.nf" linenums="7"
    // Primäre Eingabe (Array von drei Proben)
        reads_bam = [
            "${projectDir}/data/bam/reads_mother.bam",
            "${projectDir}/data/bam/reads_father.bam",
            "${projectDir}/data/bam/reads_son.bam"
        ]
    ```

Auf diese Weise können wir weiterhin faul sein, aber die Liste der Dateien lebt nicht mehr im Workflow-Code selbst, was ein großer Schritt in die richtige Richtung ist.

### 4.3. Aktualisiere die Channel-Factory, um Zeilen aus einer Datei zu lesen

Derzeit behandelt unsere Eingabe-Channel-Factory alle Dateien, die wir ihr geben, als die Dateneingaben, die wir dem Indexierungsprozess zuführen möchten.
Da wir ihr jetzt eine Datei geben, die Eingabedateipfade auflistet, müssen wir ihr Verhalten ändern, um die Datei zu parsen und die darin enthaltenen Dateipfade als Dateneingaben zu behandeln.

Glücklicherweise können wir das sehr einfach tun, indem wir einfach den [`.splitText()`-Operator](https://www.nextflow.io/docs/latest/reference/operator.html#operator-splittext) zum Channel-Konstruktionsschritt hinzufügen.

=== "Danach"

    ```groovy title="genomics-1.nf" linenums="68"
        // Erstelle Eingabekanal aus einer Textdatei, die Eingabedateipfade auflistet
        reads_ch = channel.fromPath(params.reads_bam).splitText()
    ```

=== "Vorher"

    ```groovy title="genomics-1.nf" linenums="68"
        // Erstelle Eingabekanal (einzelne Datei über CLI-Parameter)
        reads_ch = channel.fromPath(params.reads_bam)
    ```

!!! tip "Tipp"

    Dies ist eine weitere großartige Gelegenheit, den `.view()`-Operator zu verwenden, um zu sehen, wie die Channel-Inhalte vor und nach der Anwendung eines Operators aussehen.

### 4.4. Führe den Workflow aus, um zu überprüfen, dass er korrekt funktioniert

Lass uns den Workflow noch einmal ausführen. Dies sollte dasselbe Ergebnis wie zuvor produzieren, richtig?

```bash
nextflow run genomics-1.nf -resume
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics-1.nf` [sick_albattani] DSL2 - revision: 46d84642f6

    [18/23b4bb] SAMTOOLS_INDEX (1)       | 3 of 3, cached: 3 ✔
    [12/f727bb] GATK_HAPLOTYPECALLER (3) | 3 of 3
