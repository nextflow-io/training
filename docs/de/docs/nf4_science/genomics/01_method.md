# Teil 1: Methodenübersicht und manuelles Testen

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Variant Calling ist eine genomische Analysemethode, die darauf abzielt, Variationen in einer Genomsequenz im Vergleich zu einem Referenzgenom zu identifizieren.
Hier werden wir Tools und Methoden verwenden, die für das Calling von kurzen Keimbahnvarianten entwickelt wurden, _d. h._ SNPs und Indels, in Gesamtgenom-Sequenzierungsdaten.

![GATK-Pipeline](img/gatk-pipeline.png)

Eine vollständige Variant-Calling-Pipeline umfasst typischerweise viele Schritte, einschließlich Mapping auf die Referenz (manchmal als Genom-Alignment bezeichnet) sowie Filterung und Priorisierung von Varianten.
Der Einfachheit halber konzentrieren wir uns in diesem Kurs nur auf den Variant-Calling-Teil.

### Methoden

Wir zeigen dir zwei Möglichkeiten, Variant Calling auf Gesamtgenom-Sequenzierungsproben anzuwenden, um Keimbahn-SNPs und -Indels zu identifizieren.
Zunächst beginnen wir mit einem einfachen **probenweisen Ansatz**, der Varianten unabhängig für jede Probe aufruft.
Dann zeigen wir dir einen ausgefeilteren **Joint-Calling-Ansatz**, der mehrere Proben gemeinsam analysiert und genauere und informativere Ergebnisse liefert.

Bevor wir uns mit dem Schreiben von Workflow-Code für einen der beiden Ansätze befassen, werden wir die Befehle manuell an einigen Testdaten ausprobieren.

### Datensatz

Wir stellen die folgenden Daten und zugehörigen Ressourcen bereit:

- **Ein Referenzgenom**, das aus einer kleinen Region des menschlichen Chromosoms 20 (aus hg19/b37) und seinen Hilfsdateien (Index und Sequenzwörterbuch) besteht.
- **Drei Gesamtgenom-Sequenzierungsproben**, die einem Familien-Trio (Mutter, Vater und Sohn) entsprechen und auf einen kleinen Datenausschnitt auf Chromosom 20 reduziert wurden, um die Dateigrößen klein zu halten.
  Dies sind Illumina-Short-Read-Sequenzierungsdaten, die bereits auf das Referenzgenom gemappt wurden und im [BAM](https://samtools.github.io/hts-specs/SAMv1.pdf)-Format (Binary Alignment Map, eine komprimierte Version von SAM, Sequence Alignment Map) bereitgestellt werden.
- **Eine Liste genomischer Intervalle**, d. h. Koordinaten auf dem Genom, wo unsere Proben Daten haben, die für das Calling von Varianten geeignet sind, bereitgestellt im BED-Format.

### Software

Die beiden Haupttools sind [Samtools](https://www.htslib.org/), ein weit verbreitetes Toolkit zur Manipulation von Sequenz-Alignment-Dateien, und [GATK](https://gatk.broadinstitute.org/) (Genome Analysis Toolkit), eine Sammlung von Tools zur Variantenerkennung, die am Broad Institute entwickelt wurde.

Diese Tools sind nicht in der GitHub-Codespaces-Umgebung installiert, daher verwenden wir sie über Container (siehe [Hello Containers](../../hello_nextflow/05_hello_containers.md)).

!!! note "Hinweis"

    Stelle sicher, dass du dich im Verzeichnis `nf4-science/genomics` befindest, sodass der letzte Teil des Pfads, der angezeigt wird, wenn du `pwd` eingibst, `genomics` ist.

---

## 1. Probenweises Variant Calling

Probenweises Variant Calling verarbeitet jede Probe unabhängig: Der Variant Caller untersucht die Sequenzierungsdaten für jeweils eine Probe und identifiziert Positionen, an denen sich die Probe von der Referenz unterscheidet.

In diesem Abschnitt testen wir die beiden Befehle, aus denen der probenweise Variant-Calling-Ansatz besteht: Indexierung einer BAM-Datei mit Samtools und Calling von Varianten mit GATK HaplotypeCaller.
Dies sind die Befehle, die wir in Teil 2 dieses Kurses in einen Nextflow-Workflow einbinden werden.

1. Erzeuge eine Indexdatei für eine BAM-Eingabedatei mit [Samtools](https://www.htslib.org/)
2. Führe den GATK HaplotypeCaller auf der indexierten BAM-Datei aus, um probenweise Variantenaufrufe im VCF-Format (Variant Call Format) zu generieren

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/genomics/img/hello-gatk-1.svg"
</figure>

Wir beginnen damit, die beiden Befehle an nur einer Probe zu testen.

### 1.1. Indexiere eine BAM-Eingabedatei mit Samtools

Indexdateien sind ein häufiges Merkmal bioinformatischer Dateiformate; sie enthalten Informationen über die Struktur der Hauptdatei, die es Tools wie GATK ermöglichen, auf eine Teilmenge der Daten zuzugreifen, ohne die gesamte Datei lesen zu müssen.
Dies ist wichtig, da diese Dateien sehr groß werden können.

BAM-Dateien werden oft ohne Index bereitgestellt, daher ist der erste Schritt in vielen Analyse-Workflows, einen Index mit `samtools index` zu generieren.

Wir werden einen Samtools-Container herunterladen, ihn interaktiv starten und den Befehl `samtools index` auf einer der BAM-Dateien ausführen.

#### 1.1.1. Lade den Samtools-Container herunter

Führe den Befehl `docker pull` aus, um das Samtools-Container-Image herunterzuladen:

```bash
docker pull community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464
```

??? success "Befehlsausgabe"

    ```console
    1.20--b5dfbd93de237464: Pulling from library/samtools
    6360b3717211: Pull complete
    2ec3f7ad9b3c: Pull complete
    7716ca300600: Pull complete
    4f4fb700ef54: Pull complete
    8c61d418774c: Pull complete
    03dae77ff45c: Pull complete
    aab7f787139d: Pull complete
    4f4fb700ef54: Pull complete
    837d55536720: Pull complete
    897362c12ca7: Pull complete
    3893cbe24e91: Pull complete
    d1b61e94977b: Pull complete
    c72ff66fb90f: Pull complete
    0e0388f29b6d: Pull complete
    Digest: sha256:bbfc45b4f228975bde86cba95e303dd94ecf2fdacea5bfb2e2f34b0d7b141e41
    Status: Downloaded newer image for community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464
    community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464
    ```

Falls du dieses Image noch nicht heruntergeladen hast, kann es eine Minute dauern, bis der Vorgang abgeschlossen ist.
Sobald er fertig ist, hast du eine lokale Kopie des Container-Images.

#### 1.1.2. Starte den Samtools-Container interaktiv

Um den Container interaktiv auszuführen, verwende `docker run` mit den Flags `-it`.
Die Option `-v ./data:/data` bindet das lokale Verzeichnis `data` in den Container ein, sodass die Tools auf die Eingabedateien zugreifen können.

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464
```

Deine Eingabeaufforderung ändert sich zu etwas wie `(base) root@a1b2c3d4e5f6:/tmp#`, was anzeigt, dass du dich jetzt im Container befindest.
Die Datendateien sind unter `/data` zugänglich.

#### 1.1.3. Führe den Indexierungsbefehl aus

Die [Samtools-Dokumentation](https://www.htslib.org/doc/samtools-index.html) gibt uns die Befehlszeile, um eine BAM-Datei zu indexieren.

Wir müssen nur die Eingabedatei angeben; das Tool generiert automatisch einen Namen für die Ausgabe, indem es `.bai` an den Eingabedateinamen anhängt.

```bash
samtools index /data/bam/reads_mother.bam
```

??? abstract "Verzeichnisinhalt"

    ```console
    data/bam/
    ├── reads_father.bam
    ├── reads_mother.bam
    ├── reads_mother.bam.bai
    └── reads_son.bam
    ```

Du solltest jetzt eine Datei namens `reads_mother.bam.bai` im selben Verzeichnis wie die ursprüngliche BAM-Eingabedatei sehen.

#### 1.1.4. Verlasse den Samtools-Container

Um den Container zu verlassen, gib `exit` ein.

```bash
exit
```

Deine Eingabeaufforderung sollte jetzt wieder so sein wie vor dem Start des Containers.

### 1.2. Rufe Varianten mit GATK HaplotypeCaller auf

Wir werden einen GATK-Container herunterladen, ihn interaktiv starten und den Befehl `gatk HaplotypeCaller` auf der BAM-Datei ausführen, die wir gerade indexiert haben.

#### 1.2.1. Lade den GATK-Container herunter

Führe den Befehl `docker pull` aus, um das GATK-Container-Image herunterzuladen:

```bash
docker pull community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867
```

??? success "Befehlsausgabe"

    Einige Layer zeigen `Already exists`, weil sie mit dem Samtools-Container-Image geteilt werden, das wir zuvor heruntergeladen haben.

    ```console
    4.5.0.0--730ee8817e436867: Pulling from library/gatk4
    6360b3717211: Already exists
    2ec3f7ad9b3c: Already exists
    7716ca300600: Already exists
    4f4fb700ef54: Already exists
    8c61d418774c: Already exists
    03dae77ff45c: Already exists
    aab7f787139d: Already exists
    4f4fb700ef54: Already exists
    837d55536720: Already exists
    897362c12ca7: Already exists
    3893cbe24e91: Already exists
    d1b61e94977b: Already exists
    e5c558f54708: Pull complete
    087cce32d294: Pull complete
    Digest: sha256:e33413b9100f834fcc62fd5bc9edc1e881e820aafa606e09301eac2303d8724b
    Status: Downloaded newer image for community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867
    community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867
    ```

Dies sollte schneller gehen als der erste Download, da die beiden Container-Images die meisten ihrer Layer teilen.

#### 1.2.2. Starte den GATK-Container interaktiv

Starte den GATK-Container interaktiv mit eingebundenem Datenverzeichnis, genau wie wir es für Samtools getan haben.

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867
```

Deine Eingabeaufforderung ändert sich, um anzuzeigen, dass du dich jetzt im GATK-Container befindest.

#### 1.2.3. Führe den Variant-Calling-Befehl aus

Die [GATK-Dokumentation](https://gatk.broadinstitute.org/hc/en-us/articles/21905025322523-HaplotypeCaller) gibt uns die Befehlszeile, um Variant Calling auf einer BAM-Datei durchzuführen.

Wir müssen die BAM-Eingabedatei (`-I`) sowie das Referenzgenom (`-R`), einen Namen für die Ausgabedatei (`-O`) und eine Liste genomischer Intervalle zur Analyse (`-L`) angeben.

Wir müssen jedoch nicht den Pfad zur Indexdatei angeben; das Tool sucht automatisch danach im selben Verzeichnis, basierend auf der etablierten Namens- und Speicherortkonvention.
Dasselbe gilt für die Hilfsdateien des Referenzgenoms (Index- und Sequenzwörterbuchdateien, `*.fai` und `*.dict`).

```bash
gatk HaplotypeCaller \
        -R /data/ref/ref.fasta \
        -I /data/bam/reads_mother.bam \
        -O reads_mother.vcf \
        -L /data/ref/intervals.bed
```

??? success "Befehlsausgabe"

    Das Tool erzeugt ausführliche Logging-Ausgaben. Die hervorgehobenen Zeilen bestätigen den erfolgreichen Abschluss.

    ```console hl_lines="37 51 56 57"
    Using GATK jar /opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar
    Running:
        java -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Dsamjdk.compression_level=2 -jar /opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar HaplotypeCaller -R /data/ref/ref.fasta -I /data/bam/reads_mother.bam -O reads_mother.vcf -L /data/ref/intervals.bed
    00:27:50.687 INFO  NativeLibraryLoader - Loading libgkl_compression.so from jar:file:/opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar!/com/intel/gkl/native/libgkl_compression.so
    00:27:50.854 INFO  HaplotypeCaller - ------------------------------------------------------------
    00:27:50.858 INFO  HaplotypeCaller - The Genome Analysis Toolkit (GATK) v4.5.0.0
    00:27:50.858 INFO  HaplotypeCaller - For support and documentation go to https://software.broadinstitute.org/gatk/
    00:27:50.858 INFO  HaplotypeCaller - Executing as root@a1fe8ff42d07 on Linux v6.10.14-linuxkit amd64
    00:27:50.858 INFO  HaplotypeCaller - Java runtime: OpenJDK 64-Bit Server VM v17.0.11-internal+0-adhoc..src
    00:27:50.859 INFO  HaplotypeCaller - Start Date/Time: February 8, 2026 at 12:27:50 AM GMT
    00:27:50.859 INFO  HaplotypeCaller - ------------------------------------------------------------
    00:27:50.859 INFO  HaplotypeCaller - ------------------------------------------------------------
    00:27:50.861 INFO  HaplotypeCaller - HTSJDK Version: 4.1.0
    00:27:50.861 INFO  HaplotypeCaller - Picard Version: 3.1.1
    00:27:50.861 INFO  HaplotypeCaller - Built for Spark Version: 3.5.0
    00:27:50.862 INFO  HaplotypeCaller - HTSJDK Defaults.COMPRESSION_LEVEL : 2
    00:27:50.862 INFO  HaplotypeCaller - HTSJDK Defaults.USE_ASYNC_IO_READ_FOR_SAMTOOLS : false
    00:27:50.862 INFO  HaplotypeCaller - HTSJDK Defaults.USE_ASYNC_IO_WRITE_FOR_SAMTOOLS : true
    00:27:50.863 INFO  HaplotypeCaller - HTSJDK Defaults.USE_ASYNC_IO_WRITE_FOR_TRIBBLE : false
    00:27:50.864 INFO  HaplotypeCaller - Deflater: IntelDeflater
    00:27:50.864 INFO  HaplotypeCaller - Inflater: IntelInflater
    00:27:50.864 INFO  HaplotypeCaller - GCS max retries/reopens: 20
    00:27:50.864 INFO  HaplotypeCaller - Requester pays: disabled
    00:27:50.865 INFO  HaplotypeCaller - Initializing engine
    00:27:50.991 INFO  FeatureManager - Using codec BEDCodec to read file file:///data/ref/intervals.bed
    00:27:51.016 INFO  IntervalArgumentCollection - Processing 6369 bp from intervals
    00:27:51.029 INFO  HaplotypeCaller - Done initializing engine
    00:27:51.040 INFO  NativeLibraryLoader - Loading libgkl_utils.so from jar:file:/opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar!/com/intel/gkl/native/libgkl_utils.so
    00:27:51.042 INFO  NativeLibraryLoader - Loading libgkl_smithwaterman.so from jar:file:/opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar!/com/intel/gkl/native/libgkl_smithwaterman.so
    00:27:51.042 INFO  SmithWatermanAligner - Using AVX accelerated SmithWaterman implementation
    00:27:51.046 INFO  HaplotypeCallerEngine - Disabling physical phasing, which is supported only for reference-model confidence output
    00:27:51.063 INFO  NativeLibraryLoader - Loading libgkl_pairhmm_omp.so from jar:file:/opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar!/com/intel/gkl/native/libgkl_pairhmm_omp.so
    00:27:51.085 INFO  IntelPairHmm - Flush-to-zero (FTZ) is enabled when running PairHMM
    00:27:51.086 INFO  IntelPairHmm - Available threads: 10
    00:27:51.086 INFO  IntelPairHmm - Requested threads: 4
    00:27:51.086 INFO  PairHMM - Using the OpenMP multi-threaded AVX-accelerated native PairHMM implementation
    00:27:51.128 INFO  ProgressMeter - Starting traversal
    00:27:51.136 INFO  ProgressMeter -        Current Locus  Elapsed Minutes     Regions Processed   Regions/Minute
    00:27:51.882 WARN  InbreedingCoeff - InbreedingCoeff will not be calculated at position 20_10037292_10066351:3480 and possibly subsequent; at least 10 samples must have called genotypes
    00:27:52.969 INFO  HaplotypeCaller - 7 read(s) filtered by: MappingQualityReadFilter
    0 read(s) filtered by: MappingQualityAvailableReadFilter
    0 read(s) filtered by: MappedReadFilter
    0 read(s) filtered by: NotSecondaryAlignmentReadFilter
    0 read(s) filtered by: NotDuplicateReadFilter
    0 read(s) filtered by: PassesVendorQualityCheckReadFilter
    0 read(s) filtered by: NonZeroReferenceLengthAlignmentReadFilter
    0 read(s) filtered by: GoodCigarReadFilter
    0 read(s) filtered by: WellformedReadFilter
    7 total reads filtered out of 1867 reads processed
    00:27:52.971 INFO  ProgressMeter - 20_10037292_10066351:13499              0.0                    35           1145.7
    00:27:52.971 INFO  ProgressMeter - Traversal complete. Processed 35 total regions in 0.0 minutes.
    00:27:52.976 INFO  VectorLoglessPairHMM - Time spent in setup for JNI call : 0.003346916
    00:27:52.976 INFO  PairHMM - Total compute time in PairHMM computeLogLikelihoods() : 0.045731709
    00:27:52.977 INFO  SmithWatermanAligner - Total compute time in native Smith-Waterman : 0.02 sec
    00:27:52.981 INFO  HaplotypeCaller - Shutting down engine
    [February 8, 2026 at 12:27:52 AM GMT] org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCaller done. Elapsed time: 0.04 minutes.
    Runtime.totalMemory()=203423744
    ```

Die Ausgabedatei `reads_mother.vcf` wird in deinem Arbeitsverzeichnis im Container erstellt, sodass du sie im VS Code-Datei-Explorer nicht sehen wirst, es sei denn, du änderst den Ausgabedateipfad.
Es ist jedoch eine kleine Testdatei, sodass du sie mit `cat` öffnen und den Inhalt anzeigen kannst.
Wenn du ganz nach oben zum Anfang der Datei scrollst, findest du einen Header, der aus vielen Zeilen mit Metadaten besteht, gefolgt von einer Liste von Variantenaufrufen, einer pro Zeile.

??? abstract "Dateiinhalt"

    ```console title="reads_mother.vcf" linenums="26"
    #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	reads_mother
    20_10037292_10066351	3480	.	C	CT	503.03	.	AC=2;AF=1.00;AN=2;DP=23;ExcessHet=0.0000;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=27.95;SOR=1.179	GT:AD:DP:GQ:PL	1/1:0,18:18:54:517,54,0
    20_10037292_10066351	3520	.	AT	A	609.03	.	AC=2;AF=1.00;AN=2;DP=18;ExcessHet=0.0000;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=33.83;SOR=0.693	GT:AD:DP:GQ:PL	1/1:0,18:18:54:623,54,0
    20_10037292_10066351	3529	.	T	A	155.64	.	AC=1;AF=0.500;AN=2;BaseQRankSum=-0.544;DP=21;ExcessHet=0.0000;FS=1.871;MLEAC=1;MLEAF=0.500;MQ=60.00;MQRankSum=0.000;QD=7.78;ReadPosRankSum=-1.158;SOR=1.034	GT:AD:DP:GQ:PL	0/1:12,8:20:99:163,0,328
    ```

Jede Zeile beschreibt eine mögliche Variante, die in den Sequenzierungsdaten der Probe identifiziert wurde. Für eine Anleitung zur Interpretation des VCF-Formats siehe [diesen hilfreichen Artikel](https://www.ebi.ac.uk/training/online/courses/human-genetic-variation-introduction/variant-identification-and-analysis/understanding-vcf-format/).

Die Ausgabe-VCF-Datei wird von einer Indexdatei namens `reads_mother.vcf.idx` begleitet, die automatisch von GATK erstellt wurde.
Sie hat dieselbe Funktion wie die BAM-Indexdatei: Sie ermöglicht es Tools, Teilmengen von Daten zu suchen und abzurufen, ohne die gesamte Datei laden zu müssen.

#### 1.2.4. Verlasse den GATK-Container

Um den Container zu verlassen, gib `exit` ein.

```bash
exit
```

Deine Eingabeaufforderung sollte wieder normal sein.
Damit ist der Test des probenweisen Variant Callings abgeschlossen.

---

## 2. Joint Calling auf einer Kohorte

Der Variant-Calling-Ansatz, den wir gerade verwendet haben, generiert Variantenaufrufe pro Probe.
Das ist in Ordnung, um Varianten von jeder Probe isoliert zu betrachten, liefert aber begrenzte Informationen.
Es ist oft interessanter zu sehen, wie sich Variantenaufrufe über mehrere Proben hinweg unterscheiden.
GATK bietet für diesen Zweck eine alternative Methode namens Joint Variant Calling.

Joint Variant Calling beinhaltet die Erzeugung einer speziellen Art von Variantenausgabe namens GVCF (für Genomic VCF) für jede Probe, dann das Kombinieren der GVCF-Daten von allen Proben und das Ausführen einer statistischen Analyse namens „Joint Genotyping".

![Joint-Analyse](img/joint-calling.png)

Das Besondere an einem GVCF einer Probe ist, dass es Datensätze enthält, die Sequenzdatenstatistiken über alle Positionen im Zielbereich des Genoms zusammenfassen, nicht nur die Positionen, an denen das Programm Hinweise auf Variation gefunden hat.
Dies ist entscheidend für die Joint-Genotyping-Berechnung ([weiterführende Literatur](https://gatk.broadinstitute.org/hc/en-us/articles/360035890431-The-logic-of-joint-calling-for-germline-short-variants)).

Das GVCF wird von GATK HaplotypeCaller erzeugt, demselben Tool, das wir gerade getestet haben, mit einem zusätzlichen Parameter (`-ERC GVCF`).
Das Kombinieren der GVCFs erfolgt mit GATK GenomicsDBImport, das die probenweisen Aufrufe in einen Datenspeicher (analog zu einer Datenbank) kombiniert.
Die eigentliche „Joint-Genotyping"-Analyse wird dann mit GATK GenotypeGVCFs durchgeführt.

Hier testen wir die Befehle, die zum Generieren von GVCFs und zum Ausführen von Joint Genotyping benötigt werden.
Dies sind die Befehle, die wir in Teil 3 dieses Kurses in einen Nextflow-Workflow einbinden werden.

1. Erzeuge eine Indexdatei für jede BAM-Eingabedatei mit Samtools
2. Führe den GATK HaplotypeCaller auf jeder BAM-Eingabedatei aus, um ein GVCF mit probenweisen genomischen Variantenaufrufen zu generieren
3. Sammle alle GVCFs und kombiniere sie in einen GenomicsDB-Datenspeicher
4. Führe Joint Genotyping auf dem kombinierten GVCF-Datenspeicher aus, um ein VCF auf Kohortenebene zu erzeugen

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/genomics/img/hello-gatk-2.svg"
</figure>

Wir müssen jetzt alle diese Befehle testen, beginnend mit der Indexierung aller drei BAM-Dateien.

### 2.1. Indexiere BAM-Dateien für alle drei Proben

Im ersten Abschnitt oben haben wir nur eine BAM-Datei indexiert.
Jetzt müssen wir alle drei Proben indexieren, damit GATK HaplotypeCaller sie verarbeiten kann.

#### 2.1.1. Starte den Samtools-Container interaktiv

Wir haben das Samtools-Container-Image bereits heruntergeladen, sodass wir es direkt starten können:

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464
```

Deine Eingabeaufforderung ändert sich, um anzuzeigen, dass du dich im Container befindest, mit dem Datenverzeichnis wie zuvor eingebunden.

#### 2.1.2. Führe den Indexierungsbefehl auf allen drei Proben aus

Führe den Indexierungsbefehl auf jeder der drei BAM-Dateien aus:

```bash
samtools index /data/bam/reads_mother.bam
samtools index /data/bam/reads_father.bam
samtools index /data/bam/reads_son.bam
```

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

Dies sollte die Indexdateien im selben Verzeichnis wie die entsprechenden BAM-Dateien erzeugen.

#### 2.1.3. Verlasse den Samtools-Container

Um den Container zu verlassen, gib `exit` ein.

```bash
exit
```

Deine Eingabeaufforderung sollte wieder normal sein.

### 2.2. Generiere GVCFs für alle drei Proben

Um den Joint-Genotyping-Schritt auszuführen, benötigen wir GVCFs für alle drei Proben.

#### 2.2.1. Starte den GATK-Container interaktiv

Wir haben das GATK-Container-Image bereits früher heruntergeladen, sodass wir es direkt starten können:

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867
```

Deine Eingabeaufforderung ändert sich, um anzuzeigen, dass du dich im GATK-Container befindest.

#### 2.2.2. Führe den Variant-Calling-Befehl mit der GVCF-Option aus

Um ein genomisches VCF (GVCF) zu erzeugen, fügen wir die Option `-ERC GVCF` zum Basisbefehl hinzu, die den GVCF-Modus des HaplotypeCaller aktiviert.

Wir ändern auch die Dateierweiterung für die Ausgabedatei von `.vcf` zu `.g.vcf`.
Dies ist technisch keine Anforderung, aber eine dringend empfohlene Konvention.

```bash
gatk HaplotypeCaller \
        -R /data/ref/ref.fasta \
        -I /data/bam/reads_mother.bam \
        -O reads_mother.g.vcf \
        -L /data/ref/intervals.bed \
        -ERC GVCF
```

??? success "Befehlsausgabe"

    ```console hl_lines="39 53 58 59"
    Using GATK jar /opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar
    Running:
        java -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Dsamjdk.compression_level=2 -jar /opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar HaplotypeCaller -R /data/ref/ref.fasta -I /data/bam/reads_mother.bam -O reads_mother.g.vcf -L /data/ref/intervals.bed -ERC GVCF
    00:28:03.593 INFO  NativeLibraryLoader - Loading libgkl_compression.so from jar:file:/opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar!/com/intel/gkl/native/libgkl_compression.so
    00:28:03.765 INFO  HaplotypeCaller - ------------------------------------------------------------
    00:28:03.768 INFO  HaplotypeCaller - The Genome Analysis Toolkit (GATK) v4.5.0.0
    00:28:03.768 INFO  HaplotypeCaller - For support and documentation go to https://software.broadinstitute.org/gatk/
    00:28:03.768 INFO  HaplotypeCaller - Executing as root@8515e5a0598e on Linux v6.10.14-linuxkit amd64
    00:28:03.768 INFO  HaplotypeCaller - Java runtime: OpenJDK 64-Bit Server VM v17.0.11-internal+0-adhoc..src
    00:28:03.769 INFO  HaplotypeCaller - Start Date/Time: February 8, 2026 at 12:28:03 AM GMT
    00:28:03.769 INFO  HaplotypeCaller - ------------------------------------------------------------
    00:28:03.770 INFO  HaplotypeCaller - ------------------------------------------------------------
    00:28:03.772 INFO  HaplotypeCaller - HTSJDK Version: 4.1.0
    00:28:03.773 INFO  HaplotypeCaller - Picard Version: 3.1.1
    00:28:03.773 INFO  HaplotypeCaller - Built for Spark Version: 3.5.0
    00:28:03.773 INFO  HaplotypeCaller - HTSJDK Defaults.COMPRESSION_LEVEL : 2
    00:28:03.773 INFO  HaplotypeCaller - HTSJDK Defaults.USE_ASYNC_IO_READ_FOR_SAMTOOLS : false
    00:28:03.773 INFO  HaplotypeCaller - HTSJDK Defaults.USE_ASYNC_IO_WRITE_FOR_SAMTOOLS : true
    00:28:03.774 INFO  HaplotypeCaller - HTSJDK Defaults.USE_ASYNC_IO_WRITE_FOR_TRIBBLE : false
    00:28:03.774 INFO  HaplotypeCaller - Deflater: IntelDeflater
    00:28:03.774 INFO  HaplotypeCaller - Inflater: IntelInflater
    00:28:03.775 INFO  HaplotypeCaller - GCS max retries/reopens: 20
    00:28:03.775 INFO  HaplotypeCaller - Requester pays: disabled
    00:28:03.776 INFO  HaplotypeCaller - Initializing engine
    00:28:03.896 INFO  FeatureManager - Using codec BEDCodec to read file file:///data/ref/intervals.bed
    00:28:03.919 INFO  IntervalArgumentCollection - Processing 6369 bp from intervals
    00:28:03.934 INFO  HaplotypeCaller - Done initializing engine
    00:28:03.935 INFO  HaplotypeCallerEngine - Tool is in reference confidence mode and the annotation, the following changes will be made to any specified annotations: 'StrandBiasBySample' will be enabled. 'ChromosomeCounts', 'FisherStrand', 'StrandOddsRatio' and 'QualByDepth' annotations have been disabled
    00:28:03.943 INFO  NativeLibraryLoader - Loading libgkl_utils.so from jar:file:/opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar!/com/intel/gkl/native/libgkl_utils.so
    00:28:03.945 INFO  NativeLibraryLoader - Loading libgkl_smithwaterman.so from jar:file:/opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar!/com/intel/gkl/native/libgkl_smithwaterman.so
    00:28:03.946 INFO  SmithWatermanAligner - Using AVX accelerated SmithWaterman implementation
    00:28:03.955 INFO  HaplotypeCallerEngine - Standard Emitting and Calling confidence set to -0.0 for reference-model confidence output
    00:28:03.956 INFO  HaplotypeCallerEngine - All sites annotated with PLs forced to true for reference-model confidence output
    00:28:03.972 INFO  NativeLibraryLoader - Loading libgkl_pairhmm_omp.so from jar:file:/opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar!/com/intel/gkl/native/libgkl_pairhmm_omp.so
    00:28:03.993 INFO  IntelPairHmm - Flush-to-zero (FTZ) is enabled when running PairHMM
    00:28:03.994 INFO  IntelPairHmm - Available threads: 10
    00:28:03.994 INFO  IntelPairHmm - Requested threads: 4
    00:28:03.994 INFO  PairHMM - Using the OpenMP multi-threaded AVX-accelerated native PairHMM implementation
    00:28:04.044 INFO  ProgressMeter - Starting traversal
    00:28:04.070 INFO  ProgressMeter -        Current Locus  Elapsed Minutes     Regions Processed   Regions/Minute
    00:28:04.874 WARN  InbreedingCoeff - InbreedingCoeff will not be calculated at position 20_10037292_10066351:3480 and possibly subsequent; at least 10 samples must have called genotypes
    00:28:06.535 INFO  HaplotypeCaller - 7 read(s) filtered by: MappingQualityReadFilter
    0 read(s) filtered by: MappingQualityAvailableReadFilter
    0 read(s) filtered by: MappedReadFilter
    0 read(s) filtered by: NotSecondaryAlignmentReadFilter
    0 read(s) filtered by: NotDuplicateReadFilter
    0 read(s) filtered by: PassesVendorQualityCheckReadFilter
    0 read(s) filtered by: NonZeroReferenceLengthAlignmentReadFilter
    0 read(s) filtered by: GoodCigarReadFilter
    0 read(s) filtered by: WellformedReadFilter
    7 total reads filtered out of 1867 reads processed
    00:28:06.537 INFO  ProgressMeter - 20_10037292_10066351:13499              0.0                    35            851.6
    00:28:06.538 INFO  ProgressMeter - Traversal complete. Processed 35 total regions in 0.0 minutes.
    00:28:06.543 INFO  VectorLoglessPairHMM - Time spent in setup for JNI call : 0.003648749
    00:28:06.544 INFO  PairHMM - Total compute time in PairHMM computeLogLikelihoods() : 0.031498916
    00:28:06.544 INFO  SmithWatermanAligner - Total compute time in native Smith-Waterman : 0.02 sec
    00:28:06.547 INFO  HaplotypeCaller - Shutting down engine
    [February 8, 2026 at 12:28:06 AM GMT] org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCaller done. Elapsed time: 0.05 minutes.
    Runtime.totalMemory()=281018368
    ```

Dies erstellt die GVCF-Ausgabedatei `reads_mother.g.vcf` im aktuellen Arbeitsverzeichnis im Container.

Wenn du sie mit `cat` anzeigst, um den Inhalt zu sehen, wirst du feststellen, dass sie viel länger ist als das entsprechende VCF, das wir in Abschnitt 1 generiert haben. Du kannst nicht einmal zum Anfang der Datei hochscrollen, und die meisten Zeilen sehen ganz anders aus als das, was wir im VCF gesehen haben.

??? abstract "Dateiinhalt"

    ```console title="reads_mother.g.vcf" linenums="1674"
    20_10037292_10066351    14714   .       T       <NON_REF>       .       .       END=14718       GT:DP:GQ:MIN_DP:PL       0/0:37:99:37:0,99,1192
    20_10037292_10066351    14719   .       T       <NON_REF>       .       .       END=14719       GT:DP:GQ:MIN_DP:PL       0/0:36:82:36:0,82,1087
    20_10037292_10066351    14720   .       T       <NON_REF>       .       .       END=14737       GT:DP:GQ:MIN_DP:PL       0/0:42:99:37:0,100,1160
    ```

Diese repräsentieren Nicht-Varianten-Regionen, in denen der Variant Caller keine Hinweise auf Variation gefunden hat, sodass er einige Statistiken erfasst hat, die sein Vertrauensniveau in die Abwesenheit von Variation beschreiben.
Dies ermöglicht es, zwischen zwei sehr unterschiedlichen Fällen zu unterscheiden: (1) Es gibt qualitativ hochwertige Daten, die zeigen, dass die Probe homozygot-Referenz ist, und (2) Es sind nicht genügend gute Daten verfügbar, um eine Bestimmung in die eine oder andere Richtung vorzunehmen.

In einem GVCF gibt es typischerweise viele solcher Nicht-Varianten-Zeilen, mit einer kleineren Anzahl von Variantendatensätzen, die dazwischen verstreut sind.
Versuche, `head -176` auf dem GVCF auszuführen, um nur die ersten 176 Zeilen der Datei zu laden, um einen tatsächlichen Variantenaufruf zu finden.

??? abstract "Dateiinhalt"

    ```console title="reads_mother.g.vcf" linenums="174"
    20_10037292_10066351    3479    .       T       <NON_REF>       .       .       END=3479        GT:DP:GQ:MIN_DP:PL       0/0:34:36:34:0,36,906
    20_10037292_10066351    3480    .       C       CT,<NON_REF>    503.03  .       DP=23;ExcessHet=0.0000;MLEAC=2,0;MLEAF=1.00,0.00;RAW_MQandDP=82800,23    GT:AD:DP:GQ:PL:SB       1/1:0,18,0:18:54:517,54,0,517,54,517:0,0,7,11
    20_10037292_10066351    3481    .       T       <NON_REF>       .       .       END=3481        GT:DP:GQ:MIN_DP:PL       0/0:21:51:21:0,51,765
    ```

Die zweite Zeile zeigt den ersten Variantendatensatz in der Datei, der der ersten Variante in der VCF-Datei entspricht, die wir uns zuvor angesehen haben.

Genau wie das ursprüngliche VCF wird auch die Ausgabe-GVCF-Datei von einer Indexdatei namens `reads_mother.g.vcf.idx` begleitet.

#### 2.2.3. Wiederhole den Vorgang für die anderen beiden Proben

Generiere GVCFs für die verbleibenden zwei Proben, indem du die folgenden Befehle nacheinander ausführst.

```bash
gatk HaplotypeCaller \
        -R /data/ref/ref.fasta \
        -I /data/bam/reads_father.bam \
        -O reads_father.g.vcf \
        -L /data/ref/intervals.bed \
        -ERC GVCF
```

```bash
gatk HaplotypeCaller \
        -R /data/ref/ref.fasta \
        -I /data/bam/reads_son.bam \
        -O reads_son.g.vcf \
        -L /data/ref/intervals.bed \
        -ERC GVCF
```

Sobald dies abgeschlossen ist, solltest du drei Dateien mit der Endung `.g.vcf` in deinem aktuellen Verzeichnis haben (eine pro Probe) und ihre jeweiligen Indexdateien mit der Endung `.g.vcf.idx`.

Aber verlasse den Container noch nicht!
Wir werden denselben Container im nächsten Schritt verwenden.

### 2.3. Führe Joint Genotyping aus

Jetzt, da wir alle GVCFs haben, können wir den Joint-Genotyping-Ansatz zur Generierung von Variantenaufrufen für eine Kohorte von Proben ausprobieren.
Es ist eine zweistufige Methode, die darin besteht, die Daten aus allen GVCFs in einen Datenspeicher zu kombinieren und dann die eigentliche Joint-Genotyping-Analyse auszuführen, um das endgültige VCF mit gemeinsam aufgerufenen Varianten zu generieren.

#### 2.3.1. Kombiniere alle probenweisen GVCFs

Dieser erste Schritt verwendet ein weiteres GATK-Tool namens GenomicsDBImport, um die Daten aus allen GVCFs in einen GenomicsDB-Datenspeicher zu kombinieren.

```bash
gatk GenomicsDBImport \
    -V reads_mother.g.vcf \
    -V reads_father.g.vcf \
    -V reads_son.g.vcf \
    -L /data/ref/intervals.bed \
    --genomicsdb-workspace-path family_trio_gdb
```

??? success "Befehlsausgabe"

    ```console hl_lines="33 36 37 39 40"
    Using GATK jar /opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar
    Running:
        java -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Dsamjdk.compression_level=2 -jar /opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar GenomicsDBImport -V reads_mother.g.vcf -V reads_father.g.vcf -V reads_son.g.vcf -L /data/ref/intervals.bed --genomicsdb-workspace-path family_trio_gdb
    00:28:20.772 INFO  NativeLibraryLoader - Loading libgkl_compression.so from jar:file:/opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar!/com/intel/gkl/native/libgkl_compression.so
    00:28:20.914 INFO  GenomicsDBImport - ------------------------------------------------------------
    00:28:20.917 INFO  GenomicsDBImport - The Genome Analysis Toolkit (GATK) v4.5.0.0
    00:28:20.917 INFO  GenomicsDBImport - For support and documentation go to https://software.broadinstitute.org/gatk/
    00:28:20.917 INFO  GenomicsDBImport - Executing as root@8515e5a0598e on Linux v6.10.14-linuxkit amd64
    00:28:20.917 INFO  GenomicsDBImport - Java runtime: OpenJDK 64-Bit Server VM v17.0.11-internal+0-adhoc..src
    00:28:20.918 INFO  GenomicsDBImport - Start Date/Time: February 8, 2026 at 12:28:20 AM GMT
    00:28:20.918 INFO  GenomicsDBImport - ------------------------------------------------------------
    00:28:20.918 INFO  GenomicsDBImport - ------------------------------------------------------------
    00:28:20.920 INFO  GenomicsDBImport - HTSJDK Version: 4.1.0
    00:28:20.921 INFO  GenomicsDBImport - Picard Version: 3.1.1
    00:28:20.921 INFO  GenomicsDBImport - Built for Spark Version: 3.5.0
    00:28:20.922 INFO  GenomicsDBImport - HTSJDK Defaults.COMPRESSION_LEVEL : 2
    00:28:20.923 INFO  GenomicsDBImport - HTSJDK Defaults.USE_ASYNC_IO_READ_FOR_SAMTOOLS : false
    00:28:20.923 INFO  GenomicsDBImport - HTSJDK Defaults.USE_ASYNC_IO_WRITE_FOR_SAMTOOLS : true
    00:28:20.923 INFO  GenomicsDBImport - HTSJDK Defaults.USE_ASYNC_IO_WRITE_FOR_TRIBBLE : false
    00:28:20.923 INFO  GenomicsDBImport - Deflater: IntelDeflater
    00:28:20.924 INFO  GenomicsDBImport - Inflater: IntelInflater
    00:28:20.924 INFO  GenomicsDBImport - GCS max retries/reopens: 20
    00:28:20.924 INFO  GenomicsDBImport - Requester pays: disabled
    00:28:20.925 INFO  GenomicsDBImport - Initializing engine
    00:28:21.144 INFO  FeatureManager - Using codec BEDCodec to read file file:///data/ref/intervals.bed
    00:28:21.152 INFO  IntervalArgumentCollection - Processing 6369 bp from intervals
    00:28:21.157 INFO  GenomicsDBImport - Done initializing engine
    00:28:21.287 INFO  GenomicsDBLibLoader - GenomicsDB native library version : 1.5.1-84e800e
    00:28:21.290 INFO  GenomicsDBImport - Vid Map JSON file will be written to /tmp/family_trio_gdb/vidmap.json
    00:28:21.290 INFO  GenomicsDBImport - Callset Map JSON file will be written to /tmp/family_trio_gdb/callset.json
    00:28:21.291 INFO  GenomicsDBImport - Complete VCF Header will be written to /tmp/family_trio_gdb/vcfheader.vcf
    00:28:21.291 INFO  GenomicsDBImport - Importing to workspace - /tmp/family_trio_gdb
    00:28:21.453 INFO  GenomicsDBImport - Importing batch 1 with 3 samples
    00:28:21.757 INFO  GenomicsDBImport - Importing batch 1 with 3 samples
    00:28:21.859 INFO  GenomicsDBImport - Importing batch 1 with 3 samples
    00:28:21.979 INFO  GenomicsDBImport - Done importing batch 1/1
    00:28:21.988 INFO  GenomicsDBImport - Import completed!
    00:28:21.988 INFO  GenomicsDBImport - Shutting down engine
    [February 8, 2026 at 12:28:21 AM GMT] org.broadinstitute.hellbender.tools.genomicsdb.GenomicsDBImport done. Elapsed time: 0.02 minutes.
    Runtime.totalMemory()=305135616
    ```

Die Ausgabe dieses Schritts ist effektiv ein Verzeichnis, das eine Reihe weiterer verschachtelter Verzeichnisse enthält, die die kombinierten Variantendaten in Form mehrerer verschiedener Dateien enthalten.
Du kannst darin herumstöbern, aber du wirst schnell sehen, dass dieses Datenspeicherformat nicht dazu gedacht ist, direkt von Menschen gelesen zu werden.

!!! note "Hinweis"

    GATK enthält Tools, die es ermöglichen, Variantenaufruf-Daten aus dem Datenspeicher nach Bedarf zu inspizieren und zu extrahieren.

#### 2.3.2. Führe die eigentliche Joint-Genotyping-Analyse aus

Dieser zweite Schritt verwendet ein weiteres GATK-Tool namens GenotypeGVCFs, um Variantenstatistiken und individuelle Genotypen im Lichte der über alle Proben in der Kohorte verfügbaren Daten neu zu berechnen.

```bash
gatk GenotypeGVCFs \
    -R /data/ref/ref.fasta \
    -V gendb://family_trio_gdb \
    -O family_trio.vcf
```

??? success "Befehlsausgabe"

    ```console hl_lines="30 35 37 38"
    Using GATK jar /opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar
    Running:
        java -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Dsamjdk.compression_level=2 -jar /opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar GenotypeGVCFs -R /data/ref/ref.fasta -V gendb://family_trio_gdb -O family_trio.vcf
    00:28:24.625 INFO  NativeLibraryLoader - Loading libgkl_compression.so from jar:file:/opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar!/com/intel/gkl/native/libgkl_compression.so
    00:28:24.798 INFO  GenotypeGVCFs - ------------------------------------------------------------
    00:28:24.801 INFO  GenotypeGVCFs - The Genome Analysis Toolkit (GATK) v4.5.0.0
    00:28:24.801 INFO  GenotypeGVCFs - For support and documentation go to https://software.broadinstitute.org/gatk/
    00:28:24.801 INFO  GenotypeGVCFs - Executing as root@8515e5a0598e on Linux v6.10.14-linuxkit amd64
    00:28:24.801 INFO  GenotypeGVCFs - Java runtime: OpenJDK 64-Bit Server VM v17.0.11-internal+0-adhoc..src
    00:28:24.802 INFO  GenotypeGVCFs - Start Date/Time: February 8, 2026 at 12:28:24 AM GMT
    00:28:24.802 INFO  GenotypeGVCFs - ------------------------------------------------------------
    00:28:24.802 INFO  GenotypeGVCFs - ------------------------------------------------------------
    00:28:24.804 INFO  GenotypeGVCFs - HTSJDK Version: 4.1.0
    00:28:24.804 INFO  GenotypeGVCFs - Picard Version: 3.1.1
    00:28:24.804 INFO  GenotypeGVCFs - Built for Spark Version: 3.5.0
    00:28:24.805 INFO  GenotypeGVCFs - HTSJDK Defaults.COMPRESSION_LEVEL : 2
    00:28:24.805 INFO  GenotypeGVCFs - HTSJDK Defaults.USE_ASYNC_IO_READ_FOR_SAMTOOLS : false
    00:28:24.806 INFO  GenotypeGVCFs - HTSJDK Defaults.USE_ASYNC_IO_WRITE_FOR_SAMTOOLS : true
    00:28:24.806 INFO  GenotypeGVCFs - HTSJDK Defaults.USE_ASYNC_IO_WRITE_FOR_TRIBBLE : false
    00:28:24.806 INFO  GenotypeGVCFs - Deflater: IntelDeflater
    00:28:24.806 INFO  GenotypeGVCFs - Inflater: IntelInflater
    00:28:24.807 INFO  GenotypeGVCFs - GCS max retries/reopens: 20
    00:28:24.807 INFO  GenotypeGVCFs - Requester pays: disabled
    00:28:24.808 INFO  GenotypeGVCFs - Initializing engine
    00:28:25.023 INFO  GenomicsDBLibLoader - GenomicsDB native library version : 1.5.1-84e800e
    00:28:25.081 INFO  NativeGenomicsDB - pid=162 tid=163 No valid combination operation found for INFO field InbreedingCoeff  - the field will NOT be part of INFO fields in the generated VCF records
    00:28:25.082 INFO  NativeGenomicsDB - pid=162 tid=163 No valid combination operation found for INFO field MLEAC  - the field will NOT be part of INFO fields in the generated VCF records
    00:28:25.082 INFO  NativeGenomicsDB - pid=162 tid=163 No valid combination operation found for INFO field MLEAF  - the field will NOT be part of INFO fields in the generated VCF records
    00:28:25.109 INFO  GenotypeGVCFs - Done initializing engine
    00:28:25.184 INFO  ProgressMeter - Starting traversal
    00:28:25.187 INFO  ProgressMeter -        Current Locus  Elapsed Minutes    Variants Processed  Variants/Minute
    00:28:25.446 WARN  InbreedingCoeff - InbreedingCoeff will not be calculated at position 20_10037292_10066351:3480 and possibly subsequent; at least 10 samples must have called genotypes
    GENOMICSDB_TIMER,GenomicsDB iterator next() timer,Wall-clock time(s),0.15034835899999904,Cpu time(s),0.1355218420000006
    00:28:26.189 INFO  ProgressMeter - 20_10037292_10066351:13953              0.0                  3390         202994.0
    00:28:26.190 INFO  ProgressMeter - Traversal complete. Processed 3390 total variants in 0.0 minutes.
    00:28:26.194 INFO  GenotypeGVCFs - Shutting down engine
    [February 8, 2026 at 12:28:26 AM GMT] org.broadinstitute.hellbender.tools.walkers.GenotypeGVCFs done. Elapsed time: 0.03 minutes.
    Runtime.totalMemory()=296747008
    ```

Dies erstellt die VCF-Ausgabedatei `family_trio.vcf` im aktuellen Arbeitsverzeichnis im Container.
Es ist eine weitere relativ kleine Datei, sodass du diese Datei mit `cat` anzeigen kannst, um ihren Inhalt zu sehen, und nach oben scrollen kannst, um die ersten paar Variantenzeilen zu finden.

??? abstract "Dateiinhalt"

    ```console title="family_trio.vcf" linenums="40"
    #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  reads_father    reads_mother    reads_son
    20_10037292_10066351    3480    .       C       CT      1625.89 .       AC=5;AF=0.833;AN=6;BaseQRankSum=0.220;DP=85;ExcessHet=0.0000;FS=2.476;MLEAC=5;MLEAF=0.833;MQ=60.00;MQRankSum=0.00;QD=21.68;ReadPosRankSum=-1.147e+00;SOR=0.487    GT:AD:DP:GQ:PL  0/1:15,16:31:99:367,0,375       1/1:0,18:18:54:517,54,0 1/1:0,26:26:78:756,78,0
    20_10037292_10066351    3520    .       AT      A       1678.89 .       AC=5;AF=0.833;AN=6;BaseQRankSum=1.03;DP=80;ExcessHet=0.0000;FS=2.290;MLEAC=5;MLEAF=0.833;MQ=60.00;MQRankSum=0.00;QD=22.39;ReadPosRankSum=0.701;SOR=0.730 GT:AD:DP:GQ:PL   0/1:18,13:31:99:296,0,424       1/1:0,18:18:54:623,54,0 1/1:0,26:26:78:774,78,0
    20_10037292_10066351    3529    .       T       A       154.29  .       AC=1;AF=0.167;AN=6;BaseQRankSum=-5.440e-01;DP=104;ExcessHet=0.0000;FS=1.871;MLEAC=1;MLEAF=0.167;MQ=60.00;MQRankSum=0.00;QD=7.71;ReadPosRankSum=-1.158e+00;SOR=1.034       GT:AD:DP:GQ:PL  0/0:44,0:44:99:0,112,1347       0/1:12,8:20:99:163,0,328        0/0:39,0:39:99:0,105,1194
    ```

Dies sieht ähnlich aus wie das VCF, das wir zuvor generiert haben, außer dass wir diesmal Genotyp-Informationen für alle drei Proben haben.
Die letzten drei Spalten in der Datei sind die Genotyp-Blöcke für die Proben, aufgelistet in alphabetischer Reihenfolge.

Wenn wir uns die aufgerufenen Genotypen für unser Test-Familien-Trio für die allererste Variante ansehen, sehen wir, dass der Vater heterozygot-variant (`0/1`) ist und die Mutter und der Sohn beide homozygot-variant (`1/1`) sind.

Das ist letztendlich die Information, die wir aus dem Datensatz extrahieren möchten!

#### 2.3.3. Verlasse den GATK-Container

Um den Container zu verlassen, gib `exit` ein.

```bash
exit
```

Deine Eingabeaufforderung sollte wieder normal sein.
Damit ist das manuelle Testen der Variant-Calling-Befehle abgeschlossen.

---

### Fazit

Du weißt, wie du die Samtools-Indexierungs- und GATK-Variant-Calling-Befehle in ihren jeweiligen Containern testest, einschließlich der Generierung von GVCFs und der Ausführung von Joint Genotyping auf mehreren Proben.

### Wie geht es weiter?

Lerne, wie du dieselben Befehle in Workflows einbindest, die Container verwenden, um die Arbeit auszuführen.
