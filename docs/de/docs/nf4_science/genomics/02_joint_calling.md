# Teil 2: Joint Calling auf einer Kohorte

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Im ersten Teil dieses Kurses hast du eine Variant-Calling-Pipeline erstellt, die komplett linear war und die Daten jeder Probe unabhängig von den anderen verarbeitet hat.
In einem echten Genomik-Anwendungsfall musst du jedoch normalerweise die Variant Calls mehrerer Proben gemeinsam betrachten.

In diesem zweiten Teil zeigen wir dir, wie du Channels und Channel-Operatoren verwendest, um Joint Variant Calling mit GATK zu implementieren, aufbauend auf der Pipeline aus Teil 1.

### Übersicht über die Methode

Die GATK-Variant-Calling-Methode, die wir im ersten Teil dieses Kurses verwendet haben, hat einfach Variant Calls pro Probe generiert.
Das ist in Ordnung, wenn du dir nur die Varianten jeder Probe isoliert ansehen möchtest, aber das liefert nur begrenzte Informationen.
Es ist oft interessanter zu betrachten, wie sich Variant Calls über mehrere Proben hinweg unterscheiden, und dafür bietet GATK eine alternative Methode namens Joint Variant Calling, die wir hier demonstrieren.

Joint Variant Calling beinhaltet die Generierung einer speziellen Art von Variant-Ausgabe namens GVCF (für Genomic VCF) für jede Probe, dann die Kombination der GVCF-Daten aller Proben und schließlich die Durchführung einer statistischen Analyse namens 'Joint Genotyping'.

![Joint-Analyse](img/joint-calling.png)

Das Besondere an einer Proben-GVCF ist, dass sie Einträge enthält, die Sequenzdatenstatistiken über alle Positionen im Zielbereich des Genoms zusammenfassen, nicht nur die Positionen, an denen das Programm Hinweise auf Variation gefunden hat.
Das ist entscheidend für die Joint-Genotyping-Berechnung ([weiterführende Informationen](https://gatk.broadinstitute.org/hc/en-us/articles/360035890431-The-logic-of-joint-calling-for-germline-short-variants)).

Die GVCF wird von GATK HaplotypeCaller erzeugt, demselben Tool, das wir in Teil 1 verwendet haben, mit einem zusätzlichen Parameter (`-ERC GVCF`).
Die Kombination der GVCFs erfolgt mit GATK GenomicsDBImport, das die Calls pro Probe in einen Datenspeicher (analog zu einer Datenbank) kombiniert, dann wird die eigentliche 'Joint Genotyping'-Analyse mit GATK GenotypeGVCFs durchgeführt.

### Workflow

Um zusammenzufassen, in diesem Teil des Kurses werden wir einen Workflow entwickeln, der Folgendes tut:

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/genomics/img/hello-gatk-2.svg"
</figure>

1. Generiere eine Indexdatei für jede BAM-Eingabedatei mit Samtools
2. Führe GATK HaplotypeCaller auf jeder BAM-Eingabedatei aus, um eine GVCF mit genomischen Variant Calls pro Probe zu generieren
3. Sammle alle GVCFs und kombiniere sie in einen GenomicsDB-Datenspeicher
4. Führe Joint Genotyping auf dem kombinierten GVCF-Datenspeicher aus, um eine VCF auf Kohorten-Ebene zu erzeugen

Wir wenden dies auf denselben Datensatz wie in Teil 1 an.

---

## 0. Aufwärmen: Samtools und GATK direkt ausführen

Genau wie zuvor möchten wir die Befehle manuell ausprobieren, bevor wir versuchen, sie in einen Workflow zu verpacken.

!!! note

     Stelle sicher, dass du im richtigen Arbeitsverzeichnis bist:
     `cd /workspaces/training/nf4-science/genomics`

### 0.1. Eine BAM-Eingabedatei mit Samtools indizieren

Dieser erste Schritt ist derselbe wie in Teil 1, sollte sich also sehr vertraut anfühlen, aber dieses Mal müssen wir es für alle drei Proben tun.

!!! note

    Wir haben technisch gesehen bereits Indexdateien für die drei Proben durch unsere Pipeline generiert, also könnten wir diese aus dem Ergebnisverzeichnis holen. Es ist jedoch sauberer, dies einfach manuell zu wiederholen, und es dauert nur eine Minute.

#### 0.1.1. Den Samtools-Container interaktiv starten

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464
```

<!--
??? success "Befehlsausgabe"

    ```console

    ```
-->

#### 0.1.2. Den Indizierungsbefehl für die drei Proben ausführen

```bash
samtools index /data/bam/reads_mother.bam
samtools index /data/bam/reads_father.bam
samtools index /data/bam/reads_son.bam
```

Genau wie zuvor sollte dies die Indexdateien im selben Verzeichnis wie die entsprechenden BAM-Dateien erzeugen.

??? abstract "Verzeichnisinhalt"

    ```console
    data/bam/
    ├── reads_father.bam
    ├── reads_father.bam.bai
    ├── reads_mother.bam
    ├── reads_mother.bam.bai
    ├── reads_son.bam
    └── reads_son.bam.bai
    ```

Jetzt, da wir Indexdateien für alle drei Proben haben, können wir mit der Generierung der GVCFs für jede von ihnen fortfahren.

#### 0.1.3. Den Samtools-Container verlassen

```bash
exit
```

### 0.2. Varianten mit GATK HaplotypeCaller im GVCF-Modus aufrufen

Dieser zweite Schritt ist dem sehr ähnlich, was wir in Teil 1 getan haben: Hello Genomics, aber wir werden GATK jetzt im 'GVCF-Modus' ausführen.

#### 0.2.1. Den GATK-Container interaktiv starten

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867
```

<!--
??? success "Befehlsausgabe"

    ```console

    ```
-->

#### 0.2.2. Den Variant-Calling-Befehl mit der GVCF-Option ausführen

Um eine genomische VCF (GVCF) zu erzeugen, fügen wir die Option `-ERC GVCF` zum Basisbefehl hinzu, was den GVCF-Modus von HaplotypeCaller aktiviert.

Wir ändern auch die Dateierweiterung für die Ausgabedatei von `.vcf` zu `.g.vcf`.
Dies ist technisch gesehen keine Anforderung, aber es ist eine stark empfohlene Konvention.

```bash
gatk HaplotypeCaller \
        -R /data/ref/ref.fasta \
        -I /data/bam/reads_mother.bam \
        -O reads_mother.g.vcf \
        -L /data/ref/intervals.bed \
        -ERC GVCF
```

<!--
??? success "Befehlsausgabe"

    ```console

    ```
-->

Dies erzeugt die GVCF-Ausgabedatei `reads_mother.g.vcf` im aktuellen Arbeitsverzeichnis im Container.

Wenn du sie mit `cat` anzeigst, um den Inhalt zu sehen, wirst du feststellen, dass sie viel länger ist als die entsprechende VCF, die wir in Teil 1 generiert haben. Du kannst nicht einmal zum Anfang der Datei hochscrollen, und die meisten Zeilen sehen ganz anders aus als das, was wir in der VCF in Teil 1 gesehen haben.

```console title="Ausgabe" linenums="1674"
20_10037292_10066351    14714   .       T       <NON_REF>       .       .       END=14718       GT:DP:GQ:MIN_DP:PL       0/0:37:99:37:0,99,1192
20_10037292_10066351    14719   .       T       <NON_REF>       .       .       END=14719       GT:DP:GQ:MIN_DP:PL       0/0:36:82:36:0,82,1087
20_10037292_10066351    14720   .       T       <NON_REF>       .       .       END=14737       GT:DP:GQ:MIN_DP:PL       0/0:42:99:37:0,100,1160
```

Diese repräsentieren Nicht-Varianten-Regionen, in denen der Variant Caller keine Hinweise auf Variation gefunden hat, daher hat er einige Statistiken erfasst, die sein Vertrauen in das Fehlen von Variation beschreiben. Dies ermöglicht es, zwischen zwei sehr unterschiedlichen Fällen zu unterscheiden: (1) es gibt qualitativ hochwertige Daten, die zeigen, dass die Probe homozygot-Referenz ist, und (2) es sind nicht genügend gute Daten verfügbar, um eine Bestimmung vorzunehmen.

In einer GVCF gibt es typischerweise viele solcher Nicht-Varianten-Zeilen, mit einer kleineren Anzahl von Varianten-Einträgen, die dazwischen verstreut sind. Versuche, `head -176` auf der GVCF auszuführen, um nur die ersten 176 Zeilen der Datei zu laden, um einen tatsächlichen Variant Call zu finden.

```console title="Ausgabe" linenums="174"
20_10037292_10066351    3479    .       T       <NON_REF>       .       .       END=3479        GT:DP:GQ:MIN_DP:PL       0/0:34:36:34:0,36,906
20_10037292_10066351    3480    .       C       CT,<NON_REF>    503.03  .       DP=23;ExcessHet=0.0000;MLEAC=2,0;MLEAF=1.00,0.00;RAW_MQandDP=82800,23    GT:AD:DP:GQ:PL:SB       1/1:0,18,0:18:54:517,54,0,517,54,517:0,0,7,11
20_10037292_10066351    3481    .       T       <NON_REF>       .       .       END=3481        GT:DP:GQ:MIN_DP:PL       0/0:21:51:21:0,51,765
```

Die zweite Zeile zeigt den ersten Varianten-Eintrag in der Datei, der der ersten Variante in der VCF-Datei entspricht, die wir uns in Teil 1 angesehen haben.

Genau wie die ursprüngliche VCF wird auch die GVCF-Ausgabedatei von einer Indexdatei begleitet, genannt `reads_mother.g.vcf.idx`.

#### 0.2.3. Den Prozess auf den anderen beiden Proben wiederholen

Um den Joint-Genotyping-Schritt zu testen, benötigen wir GVCFs für alle drei Proben, also generieren wir diese jetzt manuell.

```bash
gatk HaplotypeCaller \
        -R /data/ref/ref.fasta \
        -I /data/bam/reads_father.bam \
        -O reads_father.g.vcf \
        -L /data/ref/intervals.bed \
        -ERC GVCF
```

<!--
??? success "Befehlsausgabe"

    ```console

    ```
-->

```bash
gatk HaplotypeCaller \
        -R /data/ref/ref.fasta \
        -I /data/bam/reads_son.bam \
        -O reads_son.g.vcf \
        -L /data/ref/intervals.bed \
        -ERC GVCF
```

<!--
??? success "Befehlsausgabe"

    ```console

    ```
-->

Sobald dies abgeschlossen ist, solltest du drei Dateien mit der Endung `.g.vcf` in deinem aktuellen Verzeichnis haben (eine pro Probe) und ihre jeweiligen Indexdateien mit der Endung `.g.vcf.idx`.

### 0.3. Joint Genotyping ausführen

Jetzt, da wir alle GVCFs haben, können wir endlich den Joint-Genotyping-Ansatz zur Generierung von Variant Calls für eine Kohorte von Proben ausprobieren.
Zur Erinnerung: Es ist eine zweistufige Methode, die darin besteht, die Daten aus allen GVCFs in einen Datenspeicher zu kombinieren und dann die eigentliche Joint-Genotyping-Analyse durchzuführen, um die finale VCF der gemeinsam aufgerufenen Varianten zu generieren.

#### 0.3.1. Alle Proben-GVCFs kombinieren

Dieser erste Schritt verwendet ein weiteres GATK-Tool namens GenomicsDBImport, um die Daten aus allen GVCFs in einen GenomicsDB-Datenspeicher zu kombinieren.

```bash
gatk GenomicsDBImport \
    -V reads_mother.g.vcf \
    -V reads_father.g.vcf \
    -V reads_son.g.vcf \
    -L /data/ref/intervals.bed \
    --genomicsdb-workspace-path family_trio_gdb
```

<!--
??? success "Befehlsausgabe"

    ```console

    ```
-->

Die Ausgabe dieses Schritts ist effektiv ein Verzeichnis, das eine Reihe weiterer verschachtelter Verzeichnisse enthält, die die kombinierten Variantendaten in Form mehrerer verschiedener Dateien enthalten.
Du kannst darin herumstöbern, aber du wirst schnell sehen, dass dieses Datenspeicherformat nicht dafür gedacht ist, direkt von Menschen gelesen zu werden.

!!! note

    GATK enthält Tools, die es ermöglichen, Variant-Call-Daten aus dem Datenspeicher bei Bedarf zu inspizieren und zu extrahieren.

#### 0.3.2. Die eigentliche Joint-Genotyping-Analyse ausführen

Dieser zweite Schritt verwendet noch ein weiteres GATK-Tool namens GenotypeGVCFs, um Variantenstatistiken und individuelle Genotypen im Lichte der über alle Proben in der Kohorte verfügbaren Daten neu zu berechnen.

```bash
gatk GenotypeGVCFs \
    -R /data/ref/ref.fasta \
    -V gendb://family_trio_gdb \
    -O family_trio.vcf
```

<!--
??? success "Befehlsausgabe"

    ```console

    ```
-->

Dies erzeugt die VCF-Ausgabedatei `family_trio.vcf` im aktuellen Arbeitsverzeichnis im Container.
Es ist eine weitere relativ kleine Datei, sodass du diese Datei mit `cat` anzeigen kannst, um ihren Inhalt zu sehen, und nach oben scrollen kannst, um die ersten Variantenzeilen zu finden.

```console title="family_trio.vcf" linenums="40"
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  reads_father    reads_mother    reads_son
20_10037292_10066351    3480    .       C       CT      1625.89 .       AC=5;AF=0.833;AN=6;BaseQRankSum=0.220;DP=85;ExcessHet=0.0000;FS=2.476;MLEAC=5;MLEAF=0.833;MQ=60.00;MQRankSum=0.00;QD=21.68;ReadPosRankSum=-1.147e+00;SOR=0.487    GT:AD:DP:GQ:PL  0/1:15,16:31:99:367,0,375       1/1:0,18:18:54:517,54,0 1/1:0,26:26:78:756,78,0
20_10037292_10066351    3520    .       AT      A       1678.89 .       AC=5;AF=0.833;AN=6;BaseQRankSum=1.03;DP=80;ExcessHet=0.0000;FS=2.290;MLEAC=5;MLEAF=0.833;MQ=60.00;MQRankSum=0.00;QD=22.39;ReadPosRankSum=0.701;SOR=0.730 GT:AD:DP:GQ:PL   0/1:18,13:31:99:296,0,424       1/1:0,18:18:54:623,54,0 1/1:0,26:26:78:774,78,0
20_10037292_10066351    3529    .       T       A       154.29  .       AC=1;AF=0.167;AN=6;BaseQRankSum=-5.440e-01;DP=104;ExcessHet=0.0000;FS=1.871;MLEAC=1;MLEAF=0.167;MQ=60.00;MQRankSum=0.00;QD=7.71;ReadPosRankSum=-1.158e+00;SOR=1.034       GT:AD:DP:GQ:PL  0/0:44,0:44:99:0,112,1347       0/1:12,8:20:99:163,0,328        0/0:39,0:39:99:0,105,1194
```

Dies sieht eher wie die ursprüngliche VCF aus, die wir in Teil 1 generiert haben, außer dass wir diesmal Informationen auf Genotyp-Ebene für alle drei Proben haben.
Die letzten drei Spalten in der Datei sind die Genotyp-Blöcke für die Proben, aufgelistet in alphabetischer Reihenfolge.

Wenn wir uns die aufgerufenen Genotypen für unser Test-Familien-Trio für die allererste Variante ansehen, sehen wir, dass der Vater heterozygot-variant (`0/1`) ist und Mutter und Sohn beide homozygot-variant (`1/1`) sind.

Das ist letztendlich die Information, die wir aus dem Datensatz extrahieren möchten! Also lass uns all dies in einen Nextflow-Workflow packen, damit wir dies im großen Maßstab tun können.

#### 0.3.3. Den GATK-Container verlassen

```bash
exit
```

### Erkenntnis

Du weißt, wie du die einzelnen Befehle, die am Joint Variant Calling beteiligt sind, im Terminal ausführst, um zu überprüfen, dass sie die gewünschten Informationen erzeugen.

### Was kommt als Nächstes?

Diese Befehle in eine tatsächliche Pipeline verpacken.

---

## 1. Den Variant-Calling-Schritt pro Probe modifizieren, um eine GVCF zu erzeugen

Die gute Nachricht ist, dass wir nicht ganz von vorne anfangen müssen, da wir bereits einen Workflow geschrieben haben, der einen Teil dieser Arbeit in Teil 1 erledigt.
Allerdings erzeugt diese Pipeline VCF-Dateien, wohingegen wir jetzt GVCF-Dateien wollen, um das Joint Genotyping durchzuführen.
Wir müssen also damit beginnen, den GVCF-Variant-Calling-Modus zu aktivieren und die Dateierweiterung der Ausgabe zu aktualisieren.

!!! note

    Der Einfachheit halber werden wir mit einer neuen Kopie des GATK-Workflows arbeiten, wie er am Ende von Teil 1 steht, aber unter einem anderen Namen: `genomics-2.nf`.

### 1.1. HaplotypeCaller mitteilen, eine GVCF auszugeben, und die Ausgabeerweiterung aktualisieren

Öffnen wir die Datei `genomics-2.nf` im Code-Editor.
Sie sollte sehr vertraut aussehen, aber du kannst sie gerne ausführen, wenn du dich davon überzeugen möchtest, dass sie wie erwartet läuft.

Wir werden damit beginnen, zwei Änderungen vorzunehmen:

- Füge den Parameter `-ERC GVCF` zum GATK-HaplotypeCaller-Befehl hinzu;
- Aktualisiere den Ausgabedateipfad, um die entsprechende `.g.vcf`-Erweiterung zu verwenden, gemäß der GATK-Konvention.

Stelle sicher, dass du einen Backslash (`\`) am Ende der vorherigen Zeile hinzufügst, wenn du `-ERC GVCF` hinzufügst.

=== "Danach"

    ```groovy title="genomics-2.nf" linenums="56" hl_lines="4 6"
        """
        gatk HaplotypeCaller \
            -R ${ref_fasta} \
            -I ${input_bam} \
            -O ${input_bam}.g.vcf \
            -L ${interval_list} \
            -ERC GVCF
        """
    ```

=== "Davor"

    ```groovy title="genomics-2.nf" linenums="56" hl_lines="4"
        """
        gatk HaplotypeCaller \
            -R ${ref_fasta} \
            -I ${input_bam} \
            -O ${input_bam}.vcf \
            -L ${interval_list}
        """
    ```

Und das ist alles, was nötig ist, um HaplotypeCaller dazu zu bringen, GVCFs statt VCFs zu generieren, oder?

### 1.2. Die Pipeline ausführen, um zu überprüfen, dass du GVCFs generieren kannst

Der Nextflow-Ausführungsbefehl ist derselbe wie zuvor, bis auf den Workflow-Dateinamen selbst.
Stelle sicher, dass du ihn entsprechend aktualisierst.

```bash
nextflow run genomics-2.nf
```

??? failure "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics-2.nf` [crazy_venter] DSL2 - revision: a2d6f6f09f

    executor >  local (6)
    [f1/8d8486] SAMTOOLS_INDEX (1)       | 3 of 3 ✔
    [72/3249ca] GATK_HAPLOTYPECALLER (3) | 0 of 3
    ERROR ~ Error executing process > 'GATK_HAPLOTYPECALLER (2)'

    Caused by:
      Missing output file(s) `reads_son.bam.vcf` expected by process `GATK_HAPLOTYPECALLER (2)`

    Command executed:

      gatk HaplotypeCaller         -R ref.fasta         -I reads_son.bam         -O reads_son.bam.g.vcf         -L intervals.bed         -ERC GVCF
    ```

Und die Ausgabe ist... ganz rot! Oh nein.

Der ausgeführte Befehl ist korrekt, also hatten wir recht, dass das ausreichte, um das Verhalten des GATK-Tools zu ändern.
Aber schau dir diese Zeile über die fehlende Ausgabedatei an. Fällt dir etwas auf?

Das stimmt, wir haben vergessen, Nextflow mitzuteilen, dass es einen neuen Dateinamen erwarten soll. Ups.

### 1.3. Die Ausgabedateierweiterung auch im Prozessausgabe-Block aktualisieren

Denn es reicht nicht aus, nur die Dateierweiterung im Tool-Befehl selbst zu ändern, du musst Nextflow auch mitteilen, dass sich der erwartete Ausgabedateiname geändert hat.

=== "Danach"

    ```groovy title="genomics-2.nf" linenums="50" hl_lines="2 3"
        output:
        path "${input_bam}.g.vcf"     , emit: vcf
        path "${input_bam}.g.vcf.idx" , emit: idx
    ```

=== "Davor"

    ```groovy title="genomics-2.nf" linenums="50" hl_lines="2 3"
        output:
        path "${input_bam}.vcf"     , emit: vcf
        path "${input_bam}.vcf.idx" , emit: idx
    ```

### 1.4. Die Publish-Ziele für die neuen GVCF-Ausgaben aktualisieren

Da wir jetzt GVCFs statt VCFs erzeugen, sollten wir den `publish:`-Abschnitt des Workflows aktualisieren, um aussagekräftigere Namen zu verwenden.
Wir werden die GVCF-Dateien auch der Übersichtlichkeit halber in ihr eigenes Unterverzeichnis organisieren.

=== "Danach"

    ```groovy title="genomics-2.nf" linenums="88" hl_lines="3 4"
        publish:
        indexed_bam = SAMTOOLS_INDEX.out
        gvcf = GATK_HAPLOTYPECALLER.out.vcf
        gvcf_idx = GATK_HAPLOTYPECALLER.out.idx
    ```

=== "Davor"

    ```groovy title="genomics-2.nf" linenums="88"
        publish:
        indexed_bam = SAMTOOLS_INDEX.out
        vcf = GATK_HAPLOTYPECALLER.out.vcf
        vcf_idx = GATK_HAPLOTYPECALLER.out.idx
    ```

### 1.5. Den Output-Block für die neue Verzeichnisstruktur aktualisieren

Wir müssen auch den `output`-Block aktualisieren, um die GVCF-Dateien in ein `gvcf`-Unterverzeichnis zu legen.

=== "Danach"

    ```groovy title="genomics-2.nf" linenums="94" hl_lines="3 5 6 8 9"
    output {
        indexed_bam {
            path 'indexed_bam'
        }
        gvcf {
            path 'gvcf'
        }
        gvcf_idx {
            path 'gvcf'
        }
    }
    ```

=== "Davor"

    ```groovy title="genomics-2.nf" linenums="94"
    output {
        indexed_bam {
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

### 1.6. Die Pipeline erneut ausführen

Lass uns sie diesmal mit `-resume` ausführen.

```bash
nextflow run genomics-2.nf -resume
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics-2.nf` [nostalgic_franklin] DSL2 - revision: f2c0a93c6a

    executor >  local (3)
    [cc/fbc705] SAMTOOLS_INDEX (3)       | 3 of 3, cached: 3 ✔
    [27/0d7eb9] GATK_HAPLOTYPECALLER (2) | 3 of 3 ✔
    ```

Diesmal funktioniert es.

Die Nextflow-Ausgabe selbst sieht nicht anders aus (verglichen mit einem erfolgreichen Lauf im normalen VCF-Modus), aber jetzt können wir die `.g.vcf`-Dateien und ihre jeweiligen Indexdateien für alle drei Proben finden, organisiert in Unterverzeichnissen.

??? abstract "Verzeichnisinhalt (Symlinks gekürzt)"

    ```console
    results_genomics/
    ├── gvcf/
    │   ├── reads_father.bam.g.vcf -> */27/0d7eb9*/reads_father.bam.g.vcf
    │   ├── reads_father.bam.g.vcf.idx -> */27/0d7eb9*/reads_father.bam.g.vcf.idx
    │   ├── reads_mother.bam.g.vcf -> */e4/4ed55e*/reads_mother.bam.g.vcf
    │   ├── reads_mother.bam.g.vcf.idx -> */e4/4ed55e*/reads_mother.bam.g.vcf.idx
    │   ├── reads_son.bam.g.vcf -> */08/e95962*/reads_son.bam.g.vcf
    │   └── reads_son.bam.g.vcf.idx -> */08/e95962*/reads_son.bam.g.vcf.idx
    └── indexed_bam/
        ├── reads_father.bam -> */9a/c7a873*/reads_father.bam
        ├── reads_father.bam.bai -> */9a/c7a873*/reads_father.bam.bai
        ├── reads_mother.bam -> */f1/8d8486*/reads_mother.bam
        ├── reads_mother.bam.bai -> */f1/8d8486*/reads_mother.bam.bai
        ├── reads_son.bam -> */cc/fbc705*/reads_son.bam
        └── reads_son.bam.bai -> */cc/fbc705*/reads_son.bam.bai
    ```

Wenn du eine der GVCF-Dateien öffnest und durchblätterst, kannst du überprüfen, dass GATK HaplotypeCaller GVCF-Dateien wie angefordert erzeugt hat.

### Erkenntnis

Okay, diese hier war minimal in Bezug auf das Nextflow-Lernen...
Aber es war eine schöne Gelegenheit, die Bedeutung des Prozessausgabe-Blocks zu wiederholen!

### Was kommt als Nächstes?

Lerne, die Inhalte eines Channels zu sammeln und sie als einzelne Eingabe an den nächsten Prozess weiterzugeben.

---

## 2. Die GVCF-Daten über alle Proben hinweg sammeln und kombinieren

Wir müssen jetzt die Daten aus allen Proben-GVCFs in eine Form kombinieren, die die Joint-Genotyping-Analyse unterstützt, die wir durchführen möchten.

### 2.1. Den Prozess definieren, der die GVCFs kombiniert

Zur Erinnerung an das, was wir früher im Aufwärmabschnitt getan haben: Das Kombinieren der GVCFs ist eine Aufgabe für das GATK-Tool GenomicsDBImport, das einen Datenspeicher im sogenannten GenomicsDB-Format erzeugt.

Lass uns einen neuen Prozess schreiben, um zu definieren, wie das funktionieren wird, basierend auf dem Befehl, den wir früher im Aufwärmabschnitt verwendet haben.

```groovy title="genomics-2.nf" linenums="66"
/*
 * Kombiniere GVCFs in einen GenomicsDB-Datenspeicher
 */
process GATK_GENOMICSDB {

    container "community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867"

    input:
    path all_gvcfs
    path all_idxs
    path interval_list
    val cohort_name

    output:
    path "${cohort_name}_gdb"

    script:
    """
    gatk GenomicsDBImport \
        -V ${all_gvcfs} \
        -L ${interval_list} \
        --genomicsdb-workspace-path ${cohort_name}_gdb
    """
}
```

Was denkst du, sieht vernünftig aus?

Lass es uns verdrahten und sehen, was passiert.

### 2.2. Einen `cohort_name`-Parameter mit einem Standardwert hinzufügen

Wir müssen einen beliebigen Namen für die Kohorte angeben.
Später in der Trainingsreihe wirst du lernen, wie du Proben-Metadaten für so etwas verwendest, aber vorerst deklarieren wir einfach einen CLI-Parameter mit `params` und geben ihm der Einfachheit halber einen Standardwert.

```groovy title="genomics-2.nf" linenums="16"
    // Basisname für finale Ausgabedatei
    cohort_name: String = "family_trio"
```

### 2.3. Die Ausgaben von GATK_HAPLOTYPECALLER über die Proben hinweg sammeln

Wenn wir einfach den Ausgabe-Channel vom `GATK_HAPLOTYPECALLER`-Prozess so wie er ist einstecken würden, würde Nextflow den Prozess für jede Proben-GVCF separat aufrufen.
Wir möchten jedoch alle drei GVCFs (und ihre Indexdateien) so bündeln, dass Nextflow sie alle zusammen an einen einzelnen Prozessaufruf übergibt.

Gute Nachrichten: Wir können das mit dem `collect()`-Channel-Operator tun. Fügen wir die folgenden Zeilen zum `workflow`-Body hinzu, direkt nach dem Aufruf von GATK_HAPLOTYPECALLER:

```groovy title="genomics-2.nf" linenums="118"
// Sammle Variant-Calling-Ausgaben über Proben hinweg
all_gvcfs_ch = GATK_HAPLOTYPECALLER.out.vcf.collect()
all_idxs_ch = GATK_HAPLOTYPECALLER.out.idx.collect()
```

Erscheint das etwas kompliziert? Lass uns das aufschlüsseln und in einfache Sprache übersetzen.

1. Wir nehmen den Ausgabe-Channel vom `GATK_HAPLOTYPECALLER`-Prozess, der mit der `.out`-Eigenschaft referenziert wird.
2. Jedes 'Element', das aus dem Channel kommt, ist ein Paar von Dateien: die GVCF und ihre Indexdatei, in dieser Reihenfolge, weil das die Reihenfolge ist, in der sie im Prozessausgabe-Block aufgelistet sind. Praktischerweise können wir, weil wir in der letzten Session die Ausgaben dieses Prozesses benannt haben (mit `emit:`), die GVCFs einerseits durch Hinzufügen von `.vcf` und die Indexdateien andererseits durch Hinzufügen von `.idx` nach der `.out`-Eigenschaft herauspicken. Hätten wir diese Ausgaben nicht benannt, hätten wir sie mit `.out[0]` bzw. `.out[1]` referenzieren müssen.
3. Wir hängen den `collect()`-Channel-Operator an, um alle GVCF-Dateien zusammen in ein einzelnes Element in einem neuen Channel namens `all_gvcfs_ch` zu bündeln, und machen dasselbe mit den Indexdateien, um den neuen Channel namens `all_idxs_ch` zu bilden.

!!! tip

    Wenn du dir schwer vorstellen kannst, was hier genau passiert, denke daran, dass du den `view()`-Operator verwenden kannst, um den Inhalt von Channels vor und nach der Anwendung von Channel-Operatoren zu inspizieren.

Die resultierenden `all_gvcfs_ch`- und `all_idxs_ch`-Channels sind das, was wir in den `GATK_GENOMICSDB`-Prozess einstecken werden, den wir gerade geschrieben haben.

!!! note

    Falls du dich gefragt hast: Wir sammeln die GVCFs und ihre Indexdateien separat, weil der GATK-GenomicsDBImport-Befehl nur die GVCF-Dateipfade sehen möchte. Glücklicherweise müssen wir uns nicht um die Reihenfolge der Dateien kümmern wie bei BAMs und ihrem Index in Teil 1, da Nextflow alle Dateien zusammen zur Ausführung bereitstellt.

### 2.4. Einen Aufruf zum Workflow-Block hinzufügen, um GATK_GENOMICSDB auszuführen

Wir haben einen Prozess und wir haben Eingabe-Channels. Wir müssen nur den Prozessaufruf hinzufügen.

```groovy title="genomics-2.nf" linenums="122"
    // Kombiniere GVCFs in einen GenomicsDB-Datenspeicher
    GATK_GENOMICSDB(
        all_gvcfs_ch,
        all_idxs_ch,
        intervals_file,
        params.cohort_name
    )
```

Ok, alles ist verdrahtet.

### 2.5. Den Workflow ausführen

Lass uns sehen, ob das funktioniert.

```bash
nextflow run genomics-2.nf -resume
```

??? failure "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics-2.nf` [disturbed_bell] DSL2 - revision: 57942246cc

    executor >  local (1)
    [f1/8d8486] SAMTOOLS_INDEX (1)       | 3 of 3, cached: 3 ✔
    [e4/4ed55e] GATK_HAPLOTYPECALLER (2) | 3 of 3, cached: 3 ✔
    [51/d350ea] GATK_GENOMICSDB          | 0 of 1
    ERROR ~ Error executing process > 'GATK_GENOMICSDB'

    Caused by:
      Process `GATK_GENOMICSDB` terminated with an error exit status (1)

    Command executed:

      gatk GenomicsDBImport         -V reads_son.bam.g.vcf reads_father.bam.g.vcf reads_mother.bam.g.vcf         -L intervals.bed         --genomicsdb-workspace-path family_trio_gdb
    ```

Es läuft ziemlich schnell, da wir mit `-resume` laufen, aber es schlägt fehl!

Ah. Auf der positiven Seite sehen wir, dass Nextflow den `GATK_GENOMICSDB`-Prozess aufgenommen hat und ihn speziell nur einmal aufgerufen hat.
Das deutet darauf hin, dass der `collect()`-Ansatz bis zu einem gewissen Grad funktioniert hat.
Aber, und das ist ein großes Aber, der Prozessaufruf ist fehlgeschlagen.

Wenn wir uns die Konsolenausgabe oben genauer ansehen, können wir sehen, dass der ausgeführte Befehl nicht korrekt ist.

Kannst du den Fehler erkennen?
Schau dir diesen Teil an: `-V reads_father.bam.g.vcf reads_son.bam.g.vcf reads_mother.bam.g.vcf`

Wir haben `gatk GenomicsDBImport` mehrere GVCF-Dateien für ein einzelnes `-V`-Argument gegeben, aber das Tool erwartet ein separates `-V`-Argument für jede GVCF-Datei.

Zur Erinnerung, das war der Befehl, den wir im Container ausgeführt haben:

```bash
gatk GenomicsDBImport \
    -V reads_mother.g.vcf \
    -V reads_father.g.vcf \
    -V reads_son.g.vcf \
    -L /data/ref/intervals.bed \
    --genomicsdb-workspace-path family_trio_gdb
```

Das bedeutet also, dass wir unser Bündel von GVCF-Dateien irgendwie in einen korrekt formatierten Befehlsstring umwandeln müssen.

### 2.6. Eine Befehlszeile mit einem separaten `-V`-Argument für jede Eingabe-GVCF konstruieren

Hier kommt es uns zugute, dass Nextflow auf Groovy basiert, denn es wird uns ermöglichen, einige ziemlich einfache String-Manipulationen zu verwenden, um den notwendigen Befehlsstring zu konstruieren.

Speziell unter Verwendung dieser Syntax: `all_gvcfs.collect { gvcf -> "-V ${gvcf}" }.join(' ')`

Noch einmal, lass uns das in seine Komponenten aufschlüsseln.

1. Zuerst nehmen wir den Inhalt des `all_gvcfs`-Eingabe-Channels und wenden `.collect()` darauf an (genau wie vorher).
2. Das ermöglicht es uns, jeden einzelnen GVCF-Dateipfad im Bündel an die **Closure**, `{ gvcf -> "-V ${gvcf}" }`, zu übergeben, wobei `gvcf` sich auf diesen GVCF-Dateipfad bezieht.
   Die Closure ist eine Mini-Funktion, die wir verwenden, um `-V ` dem Dateipfad voranzustellen, in der Form von `"-V ${gvcf}"`.
3. Dann verwenden wir `.join(' ')`, um alle drei Strings mit einem einzigen Leerzeichen als Trennzeichen zu verketten.

Mit einem konkreten Beispiel sieht es so aus:

1. Wir haben drei Dateien:

   `[A.ext, B.ext, C.ext]`

2. Die Closure modifiziert jede einzelne, um die Strings zu erstellen:

   `"-V A.ext", "-V B.ext", "-V C.ext"`

3. Die `.join(' ')`-Operation generiert den finalen String:

   `"-V A.ext -V B.ext -V C.ext"`

Sobald wir diesen String haben, können wir ihn einer lokalen Variable `gvcfs_line` zuweisen, die mit dem `def`-Schlüsselwort definiert ist:

`def gvcfs_line = all_gvcfs.collect { gvcf -> "-V ${gvcf}" }.join(' ')`

Ok, also haben wir unser String-Manipulations-Ding. Wo platzieren wir es?

Wir möchten, dass dies irgendwo innerhalb der Prozessdefinition steht, weil wir es machen möchten, _nachdem_ wir die GVCF-Dateipfade in den Prozess kanalisiert haben.
Das liegt daran, dass Nextflow sie als Dateipfade sehen muss, um die Dateien selbst korrekt zur Ausführung bereitzustellen.

Aber _wo_ im Prozess können wir das hinzufügen?

Fun Fact: Du kannst beliebigen Code nach `script:` und vor den `"""` hinzufügen!

Großartig, dann fügen wir unsere String-Manipulations-Zeile dort hinzu und aktualisieren den `gatk GenomicsDBImport`-Befehl, um den verketteten String zu verwenden, den sie erzeugt.

=== "Danach"

    ```groovy title="genomics-2.nf" linenums="87"  hl_lines="2 5"
        script:
        def gvcfs_line = all_gvcfs.collect { gvcf -> "-V ${gvcf}" }.join(' ')
        """
        gatk GenomicsDBImport \
            ${gvcfs_line} \
            -L ${interval_list} \
            --genomicsdb-workspace-path ${cohort_name}_gdb
        """
    ```

=== "Davor"

    ```groovy title="genomics-2.nf" linenums="87"  hl_lines="4"
        script:
        """
        gatk GenomicsDBImport \
            -V ${all_gvcfs} \
            -L ${interval_list} \
            --genomicsdb-workspace-path ${cohort_name}_gdb
        """
    ```

Das sollte alles sein, was nötig ist, um die Eingaben korrekt an `gatk GenomicsDBImport` zu übergeben.

!!! tip

    Wenn du den `gatk GenomicsDBImport`-Befehl aktualisierst, stelle sicher, dass du das `-V `-Präfix entfernst, wenn du die `${gvcfs_line}`-Variable einsetzt.

### 2.7. Den Workflow ausführen, um zu überprüfen, dass er die GenomicsDB-Ausgabe wie erwartet generiert

Alles klar, lass uns sehen, ob das das Problem behoben hat.

```bash
nextflow run genomics-2.nf -resume
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics-2.nf` [peaceful_gates] DSL2 - revision: ca0bf847ed

    executor >  local (1)
    [cc/fbc705] SAMTOOLS_INDEX (3)       | 3 of 3, cached: 3 ✔
    [27/0d7eb9] GATK_HAPLOTYPECALLER (2) | 3 of 3, cached: 3 ✔
    [76/d13861] GATK_GENOMICSDB          | 1 of 1 ✔
    ```

Aha! Es scheint jetzt zu funktionieren.

Die ersten beiden Schritte wurden erfolgreich übersprungen, und der dritte Schritt hat diesmal einwandfrei funktioniert.
Der GenomicsDB-Datenspeicher wird im Work-Verzeichnis erstellt, aber nicht in den Ergebnissen veröffentlicht, da es nur ein Zwischenformat ist, das wir für Joint Genotyping verwenden werden.

Übrigens mussten wir nichts Besonderes tun, um damit umzugehen, dass die Ausgabe ein Verzeichnis anstelle einer einzelnen Datei ist.

### Erkenntnis

Jetzt weißt du, wie du Ausgaben aus einem Channel sammelst und sie als einzelne Eingabe an einen anderen Prozess bündelst.
Du weißt auch, wie du eine Befehlszeile konstruierst, um Eingaben an ein bestimmtes Tool mit der entsprechenden Syntax zu übergeben.

### Was kommt als Nächstes?

Lerne, wie du einen zweiten Befehl zum selben Prozess hinzufügst.

---

## 3. Den Joint-Genotyping-Schritt als Teil desselben Prozesses ausführen

Jetzt, da wir die kombinierten genomischen Variant Calls haben, können wir das Joint-Genotyping-Tool ausführen, das die finale Ausgabe erzeugt, die uns tatsächlich interessiert: die VCF der Variant Calls auf Kohorten-Ebene.

Aus logistischen Gründen entscheiden wir uns, das Joint Genotyping innerhalb desselben Prozesses einzuschließen.

### 3.1. Den Prozess von GATK_GENOMICSDB in GATK_JOINTGENOTYPING umbenennen

Da der Prozess mehr als ein Tool ausführt, ändern wir seinen Namen, um sich auf die Gesamtoperation zu beziehen, anstatt auf einen einzelnen Tool-Namen.

=== "Danach"

    ```groovy title="genomics-2.nf"
    /*
     * Kombiniere GVCFs in einen GenomicsDB-Datenspeicher und führe Joint Genotyping aus, um Calls auf Kohorten-Ebene zu erzeugen
     */
    process GATK_JOINTGENOTYPING {
    ```

=== "Davor"

    ```groovy title="genomics-2.nf"
    /*
     * Kombiniere GVCFs in einen GenomicsDB-Datenspeicher
     */
    process GATK_GENOMICSDB {
    ```

Denke daran, deine Prozessnamen so beschreibend wie möglich zu halten, um die Lesbarkeit für deine Kollegen – und dein zukünftiges Ich – zu maximieren!

### 3.2. Den Joint-Genotyping-Befehl zum GATK_JOINTGENOTYPING-Prozess hinzufügen

Füge einfach den zweiten Befehl nach dem ersten innerhalb des Script-Abschnitts hinzu.

=== "Danach"

    ```groovy title="genomics-2.nf" linenums="89"  hl_lines="6-10"
        """
        gatk GenomicsDBImport \
            ${gvcfs_line} \
            -L ${interval_list} \
            --genomicsdb-workspace-path ${cohort_name}_gdb

        gatk GenotypeGVCFs \
            -R ${ref_fasta} \
            -V gendb://${cohort_name}_gdb \
            -L ${interval_list} \
            -O ${cohort_name}.joint.vcf
        """
    ```

=== "Davor"

    ```groovy title="genomics-2.nf" linenums="89"
        """
        gatk GenomicsDBImport \
            ${gvcfs_line} \
            -L ${interval_list} \
            --genomicsdb-workspace-path ${cohort_name}_gdb
        """
    ```

Die beiden Befehle werden seriell ausgeführt, auf die gleiche Weise, wie sie ausgeführt würden, wenn wir sie manuell im Terminal ausführen würden.

### 3.3. Die Referenzgenomdateien zu den GATK_JOINTGENOTYPING-Prozess-Eingabedefinitionen hinzufügen

Der zweite Befehl benötigt die Referenzgenomdateien, also müssen wir diese zu den Prozesseingaben hinzufügen.

=== "Danach"

    ```groovy title="genomics-2.nf" linenums="78"  hl_lines="6-8"
        input:
        path all_gvcfs
        path all_idxs
        path interval_list
        val cohort_name
        path ref_fasta
        path ref_index
        path ref_dict
    ```

=== "Davor"

    ```groovy title="genomics-2.nf" linenums="78"
        input:
        path all_gvcfs
        path all_idxs
        path interval_list
        val cohort_name
    ```

Es mag lästig erscheinen, diese auszutippen, aber denke daran, du tippst sie nur einmal, und dann kannst du den Workflow eine Million Mal ausführen. Lohnt sich das?

### 3.4. Die Prozessausgabedefinition aktualisieren, um die VCF mit Variant Calls auf Kohorten-Ebene auszugeben

Wir kümmern uns nicht wirklich darum, den GenomicsDB-Datenspeicher zu speichern, der nur ein Zwischenformat ist, das nur aus logistischen Gründen existiert, also können wir ihn einfach aus dem Output-Block entfernen, wenn wir möchten.

Die Ausgabe, an der wir tatsächlich interessiert sind, ist die VCF, die vom Joint-Genotyping-Befehl erzeugt wird.

=== "Danach"

    ```groovy title="genomics-2.nf" linenums="87" hl_lines="2 3"
        output:
        path "${cohort_name}.joint.vcf"     , emit: vcf
        path "${cohort_name}.joint.vcf.idx" , emit: idx
    ```

=== "Davor"

    ```groovy title="genomics-2.nf" linenums="87"
        output:
        path "${cohort_name}_gdb"
    ```

Wir sind fast fertig!

### 3.5. Den Prozessaufruf von GATK_GENOMICSDB in GATK_JOINTGENOTYPING aktualisieren

Vergessen wir nicht, den Prozessaufruf im Workflow-Body von GATK_GENOMICSDB in GATK_JOINTGENOTYPING umzubenennen. Und während wir dabei sind, sollten wir auch die Referenzgenomdateien als Eingaben hinzufügen, da wir sie dem Joint-Genotyping-Tool bereitstellen müssen.

=== "Danach"

    ```groovy title="genomics-2.nf" linenums="126"
    // Kombiniere GVCFs in einen GenomicsDB-Datenspeicher und wende Joint Genotyping an
    GATK_JOINTGENOTYPING(
        all_gvcfs_ch,
        all_idxs_ch,
        intervals_file,
        params.cohort_name,
        ref_file,
        ref_index_file,
        ref_dict_file
    )
    ```

=== "Davor"

    ```groovy title="genomics-2.nf" linenums="126"
    // Kombiniere GVCFs in einen GenomicsDB-Datenspeicher
    GATK_GENOMICSDB(
        all_gvcfs_ch,
        all_idxs_ch,
        intervals_file,
        params.cohort_name
    )
    ```

Jetzt ist der Prozess vollständig verdrahtet.

### 3.6. Die gemeinsame VCF zum Publish-Abschnitt hinzufügen

Wir müssen die gemeinsamen VCF-Ausgaben vom neuen Prozess veröffentlichen.
Füge diese Zeilen zum `publish:`-Abschnitt des Workflows hinzu:

```groovy title="genomics-2.nf" linenums="145"
    joint_vcf = GATK_JOINTGENOTYPING.out.vcf
    joint_vcf_idx = GATK_JOINTGENOTYPING.out.idx
```

### 3.7. Die gemeinsamen VCF-Ziele zum Output-Block hinzufügen

Füge schließlich Output-Ziele für die gemeinsamen VCF-Dateien hinzu.
Wir werden sie im Root des Ergebnisverzeichnisses platzieren, da dies die finale Ausgabe ist.

```groovy title="genomics-2.nf" linenums="157"
    joint_vcf {
        path '.'
    }
    joint_vcf_idx {
        path '.'
    }
```

Jetzt sollte alles vollständig verdrahtet sein.

### 3.8. Den Workflow ausführen

Schließlich können wir den modifizierten Workflow ausführen...

```bash
nextflow run genomics-2.nf -resume
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics-2.nf` [crazy_marconi] DSL2 - revision: 5da9afc841

    executor >  local (1)
    [9a/c7a873] SAMTOOLS_INDEX (2)       | 3 of 3, cached: 3 ✔
    [e4/4ed55e] GATK_HAPLOTYPECALLER (2) | 3 of 3, cached: 3 ✔
    [a6/7cc8ed] GATK_JOINTGENOTYPING     | 1 of 1 ✔
    ```

Und es funktioniert!

Du findest die finale Ausgabedatei `family_trio.joint.vcf` (und ihren Dateiindex) im Ergebnisverzeichnis.

??? abstract "Verzeichnisinhalt (Symlinks gekürzt)"

    ```console
    results_genomics/
    ├── family_trio.joint.vcf -> */a6/7cc8ed*/family_trio.joint.vcf
    ├── family_trio.joint.vcf.idx -> */a6/7cc8ed*/family_trio.joint.vcf.idx
    ├── gvcf/
    │   ├── reads_father.bam.g.vcf -> */27/0d7eb9*/reads_father.bam.g.vcf
    │   ├── reads_father.bam.g.vcf.idx -> */27/0d7eb9*/reads_father.bam.g.vcf.idx
    │   ├── reads_mother.bam.g.vcf -> */e4/4ed55e*/reads_mother.bam.g.vcf
    │   ├── reads_mother.bam.g.vcf.idx -> */e4/4ed55e*/reads_mother.bam.g.vcf.idx
    │   ├── reads_son.bam.g.vcf -> */08/e95962*/reads_son.bam.g.vcf
    │   └── reads_son.bam.g.vcf.idx -> */08/e95962*/reads_son.bam.g.vcf.idx
    └── indexed_bam/
        ├── reads_father.bam -> */9a/c7a873*/reads_father.bam
        ├── reads_father.bam.bai -> */9a/c7a873*/reads_father.bam.bai
        ├── reads_mother.bam -> */f1/8d8486*/reads_mother.bam
        ├── reads_mother.bam.bai -> */f1/8d8486*/reads_mother.bam.bai
        ├── reads_son.bam -> */cc/fbc705*/reads_son.bam
        └── reads_son.bam.bai -> */cc/fbc705*/reads_son.bam.bai
    ```

Wenn du der skeptische Typ bist, kannst du auf die gemeinsame VCF-Datei klicken, um sie zu öffnen, und überprüfen, dass der Workflow dieselben Variant Calls generiert hat, die du erhalten hast, indem du die Tools zu Beginn dieses Abschnitts manuell ausgeführt hast.

```console title="family_trio.joint.vcf" linenums="40"
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	reads_father	reads_mother	reads_son
20_10037292_10066351	3480	.	C	CT	1625.89	.	AC=5;AF=0.833;AN=6;BaseQRankSum=0.220;DP=85;ExcessHet=0.0000;FS=2.476;MLEAC=5;MLEAF=0.833;MQ=60.00;MQRankSum=0.00;QD=21.68;ReadPosRankSum=-1.147e+00;SOR=0.487	GT:AD:DP:GQ:PL	0/1:15,16:31:99:367,0,375	1/1:0,18:18:54:517,54,0	1/1:0,26:26:78:756,78,0
20_10037292_10066351	3520	.	AT	A	1678.89	.	AC=5;AF=0.833;AN=6;BaseQRankSum=1.03;DP=80;ExcessHet=0.0000;FS=2.290;MLEAC=5;MLEAF=0.833;MQ=60.00;MQRankSum=0.00;QD=22.39;ReadPosRankSum=0.701;SOR=0.730	GT:AD:DP:GQ:PL	0/1:18,13:31:99:296,0,424	1/1:0,18:18:54:623,54,0	1/1:0,26:26:78:774,78,0
20_10037292_10066351	3529	.	T	A	154.29	.	AC=1;AF=0.167;AN=6;BaseQRankSum=-5.440e-01;DP=104;ExcessHet=0.0000;FS=1.871;MLEAC=1;MLEAF=0.167;MQ=60.00;MQRankSum=0.00;QD=7.71;ReadPosRankSum=-1.158e+00;SOR=1.034	GT:AD:DP:GQ:PL	0/0:44,0:44:99:0,112,1347	0/1:12,8:20:99:163,0,328	0/0:39,0:39:99:0,105,1194
```

Du hast jetzt einen automatisierten, vollständig reproduzierbaren Joint-Variant-Calling-Workflow!

!!! note

    Bedenke, dass die Datendateien, die wir dir gegeben haben, nur einen winzigen Teil von Chromosom 20 abdecken.
    Die tatsächliche Größe eines Variant-Callsets würde in Millionen von Varianten gezählt werden.
    Deshalb verwenden wir nur winzige Teilmengen von Daten zu Trainingszwecken!

### Erkenntnis

Du weißt, wie du einige gängige
