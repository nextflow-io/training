# Parte 1: Panoramica del metodo e test manuale

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduzione assistita da IA - [ulteriori informazioni e suggerimenti per miglioramenti](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Il variant calling è un metodo di analisi genomica che mira a identificare variazioni in una sequenza genomica rispetto a un genoma di riferimento.
Qui utilizzeremo strumenti e metodi progettati per identificare varianti germinali brevi, _cioè_ SNP e indel, in dati di sequenziamento dell'intero genoma.

![Pipeline GATK](img/gatk-pipeline.png)

Una pipeline completa di variant calling comprende tipicamente molti passaggi, incluso il mapping al riferimento (a volte chiamato allineamento genomico) e il filtraggio e la prioritizzazione delle varianti.
Per semplicità, in questo corso ci concentreremo solo sulla parte di variant calling.

### Metodi

Vi mostreremo due modi per applicare il variant calling a campioni di sequenziamento dell'intero genoma per identificare SNP e indel germinali.
Inizieremo con un semplice **approccio per campione** che identifica le varianti indipendentemente per ciascun campione.
Poi vi mostreremo un approccio più sofisticato di **joint calling** che analizza più campioni insieme, producendo risultati più accurati e informativi.

Prima di addentrarci nella scrittura di qualsiasi codice per il flusso di lavoro per entrambi gli approcci, testeremo i comandi manualmente su alcuni dati di test.

### Dataset

Forniamo i seguenti dati e risorse correlate:

- **Un genoma di riferimento** costituito da una piccola regione del cromosoma umano 20 (da hg19/b37) e i suoi file accessori (indice e dizionario della sequenza).
- **Tre campioni di sequenziamento dell'intero genoma** corrispondenti a un trio familiare (madre, padre e figlio), che sono stati ridotti a una piccola porzione di dati sul cromosoma 20 per mantenere le dimensioni dei file ridotte.
  Si tratta di dati di sequenziamento Illumina a lettura breve che sono già stati mappati sul genoma di riferimento, forniti in formato [BAM](https://samtools.github.io/hts-specs/SAMv1.pdf) (Binary Alignment Map, una versione compressa di SAM, Sequence Alignment Map).
- **Una lista di intervalli genomici**, cioè coordinate sul genoma dove i nostri campioni hanno dati adatti per identificare varianti, fornita in formato BED.

### Software

I due strumenti principali coinvolti sono [Samtools](https://www.htslib.org/), un toolkit ampiamente utilizzato per manipolare file di allineamento di sequenze, e [GATK](https://gatk.broadinstitute.org/) (Genome Analysis Toolkit), un insieme di strumenti per la scoperta di varianti sviluppato al Broad Institute.

Questi strumenti non sono installati nell'ambiente GitHub Codespaces, quindi li utilizzeremo tramite container (vedere [Hello Containers](../../hello_nextflow/05_hello_containers.md)).

!!! note

     Assicuratevi di trovarvi nella directory `nf4-science/genomics` in modo che l'ultima parte del percorso mostrato quando digitate `pwd` sia `genomics`.

---

## 1. Variant calling per campione

Il variant calling per campione elabora ogni campione indipendentemente: il variant caller esamina i dati di sequenziamento per un campione alla volta e identifica le posizioni in cui il campione differisce dal riferimento.

In questa sezione testiamo i due comandi che costituiscono l'approccio di variant calling per campione: indicizzazione di un file BAM con Samtools e identificazione delle varianti con GATK HaplotypeCaller.
Questi sono i comandi che avvolgeremo in un flusso di lavoro Nextflow nella Parte 2 di questo corso.

1. Generare un file indice per un file BAM di input utilizzando [Samtools](https://www.htslib.org/)
2. Eseguire GATK HaplotypeCaller sul file BAM indicizzato per generare varianti per campione in VCF (Variant Call Format)

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/genomics/img/hello-gatk-1.svg"
</figure>

Iniziamo testando i due comandi su un solo campione.

### 1.1. Indicizzare un file BAM di input con Samtools

I file indice sono una caratteristica comune dei formati di file bioinformatici; contengono informazioni sulla struttura del file principale che consentono a strumenti come GATK di accedere a un sottoinsieme dei dati senza dover leggere l'intero file.
Questo è importante a causa delle dimensioni che questi file possono raggiungere.

I file BAM sono spesso forniti senza un indice, quindi il primo passo in molti flussi di lavoro di analisi è generarne uno utilizzando `samtools index`.

Scaricheremo un container Samtools, lo avvieremo in modo interattivo ed eseguiremo il comando `samtools index` su uno dei file BAM.

#### 1.1.1. Scaricare il container Samtools

Eseguite il comando `docker pull` per scaricare l'immagine del container Samtools:

```bash
docker pull community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464
```

??? success "Output del comando"

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

Se non avete scaricato questa immagine in precedenza, potrebbe volerci un minuto per completare.
Una volta terminato, avrete una copia locale dell'immagine del container.

#### 1.1.2. Avviare il container Samtools in modo interattivo

Per eseguire il container in modo interattivo, utilizzate `docker run` con i flag `-it`.
L'opzione `-v ./data:/data` monta la directory locale `data` nel container in modo che gli strumenti possano accedere ai file di input.

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464
```

Il vostro prompt cambia in qualcosa come `(base) root@a1b2c3d4e5f6:/tmp#`, indicando che siete ora all'interno del container.
I file di dati sono accessibili in `/data`.

#### 1.1.3. Eseguire il comando di indicizzazione

La [documentazione di Samtools](https://www.htslib.org/doc/samtools-index.html) ci fornisce la riga di comando da eseguire per indicizzare un file BAM.

Dobbiamo solo fornire il file di input; lo strumento genererà automaticamente un nome per l'output aggiungendo `.bai` al nome del file di input.

```bash
samtools index /data/bam/reads_mother.bam
```

??? abstract "Directory contents"

    ```console
    data/bam/
    ├── reads_father.bam
    ├── reads_mother.bam
    ├── reads_mother.bam.bai
    └── reads_son.bam
    ```

Dovreste ora vedere un file chiamato `reads_mother.bam.bai` nella stessa directory del file BAM di input originale.

#### 1.1.4. Uscire dal container Samtools

Per uscire dal container, digitate `exit`.

```bash
exit
```

Il vostro prompt dovrebbe ora essere tornato a quello che era prima di avviare il container.

### 1.2. Identificare varianti con GATK HaplotypeCaller

Scaricheremo un container GATK, lo avvieremo in modo interattivo ed eseguiremo il comando `gatk HaplotypeCaller` sul file BAM che abbiamo appena indicizzato.

#### 1.2.1. Scaricare il container GATK

Eseguite il comando `docker pull` per scaricare l'immagine del container GATK:

```bash
docker pull community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867
```

??? success "Output del comando"

    Alcuni layer mostrano `Already exists` perché sono condivisi con l'immagine del container Samtools che abbiamo scaricato in precedenza.

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

Questo dovrebbe essere più veloce del primo download perché le due immagini del container condividono la maggior parte dei loro layer.

#### 1.2.2. Avviare il container GATK in modo interattivo

Avviate il container GATK in modo interattivo con la directory data montata, proprio come abbiamo fatto per Samtools.

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867
```

Il vostro prompt cambia per indicare che siete ora all'interno del container GATK.

#### 1.2.3. Eseguire il comando di variant calling

La [documentazione GATK](https://gatk.broadinstitute.org/hc/en-us/articles/21905025322523-HaplotypeCaller) ci fornisce la riga di comando da eseguire per effettuare il variant calling su un file BAM.

Dobbiamo fornire il file BAM di input (`-I`) così come il genoma di riferimento (`-R`), un nome per il file di output (`-O`) e una lista di intervalli genomici da analizzare (`-L`).

Tuttavia, non dobbiamo specificare il percorso del file indice; lo strumento lo cercherà automaticamente nella stessa directory, in base alla convenzione di denominazione e co-locazione stabilita.
Lo stesso vale per i file accessori del genoma di riferimento (file indice e dizionario della sequenza, `*.fai` e `*.dict`).

```bash
gatk HaplotypeCaller \
        -R /data/ref/ref.fasta \
        -I /data/bam/reads_mother.bam \
        -O reads_mother.vcf \
        -L /data/ref/intervals.bed
```

??? success "Output del comando"

    Lo strumento produce output di log dettagliato. Le righe evidenziate confermano il completamento con successo.

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

Il file di output `reads_mother.vcf` viene creato all'interno della vostra directory di lavoro nel container, quindi non lo vedrete nell'explorer di file di VS Code a meno che non cambiate il percorso del file di output.
Tuttavia, è un file di test piccolo, quindi potete usare `cat` per aprirlo e visualizzarne il contenuto.
Se scorrete fino all'inizio del file, troverete un'intestazione composta da molte righe di metadati, seguita da una lista di varianti identificate, una per riga.

??? abstract "File contents"

    ```console title="reads_mother.vcf" linenums="26"
    #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	reads_mother
    20_10037292_10066351	3480	.	C	CT	503.03	.	AC=2;AF=1.00;AN=2;DP=23;ExcessHet=0.0000;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=27.95;SOR=1.179	GT:AD:DP:GQ:PL	1/1:0,18:18:54:517,54,0
    20_10037292_10066351	3520	.	AT	A	609.03	.	AC=2;AF=1.00;AN=2;DP=18;ExcessHet=0.0000;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=33.83;SOR=0.693	GT:AD:DP:GQ:PL	1/1:0,18:18:54:623,54,0
    20_10037292_10066351	3529	.	T	A	155.64	.	AC=1;AF=0.500;AN=2;BaseQRankSum=-0.544;DP=21;ExcessHet=0.0000;FS=1.871;MLEAC=1;MLEAF=0.500;MQ=60.00;MQRankSum=0.000;QD=7.78;ReadPosRankSum=-1.158;SOR=1.034	GT:AD:DP:GQ:PL	0/1:12,8:20:99:163,0,328
    ```

Ogni riga descrive una possibile variante identificata nei dati di sequenziamento del campione. Per una guida sull'interpretazione del formato VCF, consultate [questo utile articolo](https://www.ebi.ac.uk/training/online/courses/human-genetic-variation-introduction/variant-identification-and-analysis/understanding-vcf-format/).

Il file VCF di output è accompagnato da un file indice chiamato `reads_mother.vcf.idx` che è stato creato automaticamente da GATK.
Ha la stessa funzione del file indice BAM, permettere agli strumenti di cercare e recuperare sottoinsiemi di dati senza caricare l'intero file.

#### 1.2.4. Uscire dal container GATK

Per uscire dal container, digitate `exit`.

```bash
exit
```

Il vostro prompt dovrebbe essere tornato alla normalità.
Questo conclude il test del variant calling per campione.

---

## 2. Joint calling su una coorte

L'approccio di variant calling che abbiamo appena utilizzato genera varianti per campione.
Questo va bene per esaminare le varianti da ciascun campione in isolamento, ma fornisce informazioni limitate.
È spesso più interessante osservare come le varianti identificate differiscono tra più campioni.
GATK offre un metodo alternativo chiamato joint variant calling per questo scopo.

Il joint variant calling comporta la generazione di un tipo speciale di output di varianti chiamato GVCF (per Genomic VCF) per ciascun campione, quindi la combinazione dei dati GVCF di tutti i campioni e l'esecuzione di un'analisi statistica di 'joint genotyping'.

![Analisi congiunta](img/joint-calling.png)

Ciò che è speciale nel GVCF di un campione è che contiene record che riassumono le statistiche dei dati di sequenza per tutte le posizioni nell'area mirata del genoma, non solo le posizioni in cui il programma ha trovato evidenze di variazione.
Questo è fondamentale per il calcolo del joint genotyping ([ulteriori informazioni](https://gatk.broadinstitute.org/hc/en-us/articles/360035890431-The-logic-of-joint-calling-for-germline-short-variants)).

Il GVCF è prodotto da GATK HaplotypeCaller, lo stesso strumento che abbiamo appena testato, con un parametro aggiuntivo (`-ERC GVCF`).
La combinazione dei GVCF viene effettuata con GATK GenomicsDBImport, che combina le varianti per campione in un data store (analogo a un database).
L'analisi di 'joint genotyping' vera e propria viene quindi eseguita con GATK GenotypeGVCFs.

Qui testiamo i comandi necessari per generare i GVCF ed eseguire il joint genotyping.
Questi sono i comandi che avvolgeremo in un flusso di lavoro Nextflow nella Parte 3 di questo corso.

1. Generare un file indice per ciascun file BAM di input utilizzando Samtools
2. Eseguire GATK HaplotypeCaller su ciascun file BAM di input per generare un GVCF di varianti genomiche per campione
3. Raccogliere tutti i GVCF e combinarli in un data store GenomicsDB
4. Eseguire il joint genotyping sul data store GVCF combinato per produrre un VCF a livello di coorte

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/genomics/img/hello-gatk-2.svg"
</figure>

Dobbiamo ora testare tutti questi comandi, iniziando con l'indicizzazione di tutti e tre i file BAM.

### 2.1. Indicizzare i file BAM per tutti e tre i campioni

Nella prima sezione sopra, abbiamo indicizzato solo un file BAM.
Ora dobbiamo indicizzare tutti e tre i campioni in modo che GATK HaplotypeCaller possa elaborarli.

#### 2.1.1. Avviare il container Samtools in modo interattivo

Abbiamo già scaricato l'immagine del container Samtools, quindi possiamo avviarla direttamente:

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464
```

Il vostro prompt cambia per indicare che siete all'interno del container, con la directory data montata come prima.

#### 2.1.2. Eseguire il comando di indicizzazione su tutti e tre i campioni

Eseguite il comando di indicizzazione su ciascuno dei tre file BAM:

```bash
samtools index /data/bam/reads_mother.bam
samtools index /data/bam/reads_father.bam
samtools index /data/bam/reads_son.bam
```

??? abstract "Directory contents"

    ```console
    data/bam/
    ├── reads_father.bam
    ├── reads_father.bam.bai
    ├── reads_mother.bam
    ├── reads_mother.bam.bai
    ├── reads_son.bam
    └── reads_son.bam.bai
    ```

Questo dovrebbe produrre i file indice nella stessa directory dei corrispondenti file BAM.

#### 2.1.3. Uscire dal container Samtools

Per uscire dal container, digitate `exit`.

```bash
exit
```

Il vostro prompt dovrebbe essere tornato alla normalità.

### 2.2. Generare i GVCF per tutti e tre i campioni

Per eseguire il passaggio di joint genotyping, abbiamo bisogno dei GVCF per tutti e tre i campioni.

#### 2.2.1. Avviare il container GATK in modo interattivo

Abbiamo già scaricato l'immagine del container GATK in precedenza, quindi possiamo avviarla direttamente:

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867
```

Il vostro prompt cambia per indicare che siete all'interno del container GATK.

#### 2.2.2. Eseguire il comando di variant calling con l'opzione GVCF

Per produrre un VCF genomico (GVCF), aggiungiamo l'opzione `-ERC GVCF` al comando base, che attiva la modalità GVCF di HaplotypeCaller.

Cambiamo anche l'estensione del file per il file di output da `.vcf` a `.g.vcf`.
Tecnicamente questo non è un requisito, ma è una convenzione fortemente raccomandata.

```bash
gatk HaplotypeCaller \
        -R /data/ref/ref.fasta \
        -I /data/bam/reads_mother.bam \
        -O reads_mother.g.vcf \
        -L /data/ref/intervals.bed \
        -ERC GVCF
```

??? success "Output del comando"

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

Questo crea il file di output GVCF `reads_mother.g.vcf` nella directory di lavoro corrente nel container.

Se usate `cat` per visualizzarne il contenuto, vedrete che è molto più lungo del VCF equivalente che abbiamo generato nella sezione 1. Non potete nemmeno scorrere fino all'inizio del file, e la maggior parte delle righe appare abbastanza diversa da ciò che abbiamo visto nel VCF.

??? abstract "File contents"

    ```console title="reads_mother.g.vcf" linenums="1674"
    20_10037292_10066351    14714   .       T       <NON_REF>       .       .       END=14718       GT:DP:GQ:MIN_DP:PL       0/0:37:99:37:0,99,1192
    20_10037292_10066351    14719   .       T       <NON_REF>       .       .       END=14719       GT:DP:GQ:MIN_DP:PL       0/0:36:82:36:0,82,1087
    20_10037292_10066351    14720   .       T       <NON_REF>       .       .       END=14737       GT:DP:GQ:MIN_DP:PL       0/0:42:99:37:0,100,1160
    ```

Queste rappresentano regioni non varianti in cui il variant caller non ha trovato evidenze di variazione, quindi ha catturato alcune statistiche che descrivono il suo livello di confidenza nell'assenza di variazione.
Questo rende possibile distinguere tra due casi molto diversi: (1) ci sono dati di buona qualità che mostrano che il campione è omozigote per il riferimento, e (2) non ci sono abbastanza dati di buona qualità disponibili per fare una determinazione in un modo o nell'altro.

In un GVCF, ci sono tipicamente molte di queste righe non varianti, con un numero minore di record di varianti sparse tra di esse.
Provate a eseguire `head -176` sul GVCF per caricare solo le prime 176 righe del file per trovare una variante effettiva.

??? abstract "File contents"

    ```console title="reads_mother.g.vcf" linenums="174"
    20_10037292_10066351    3479    .       T       <NON_REF>       .       .       END=3479        GT:DP:GQ:MIN_DP:PL       0/0:34:36:34:0,36,906
    20_10037292_10066351    3480    .       C       CT,<NON_REF>    503.03  .       DP=23;ExcessHet=0.0000;MLEAC=2,0;MLEAF=1.00,0.00;RAW_MQandDP=82800,23    GT:AD:DP:GQ:PL:SB       1/1:0,18,0:18:54:517,54,0,517,54,517:0,0,7,11
    20_10037292_10066351    3481    .       T       <NON_REF>       .       .       END=3481        GT:DP:GQ:MIN_DP:PL       0/0:21:51:21:0,51,765
    ```

La seconda riga mostra il primo record di variante nel file, che corrisponde alla prima variante nel file VCF che abbiamo esaminato in precedenza.

Proprio come il VCF originale, anche il file di output GVCF è accompagnato da un file indice, chiamato `reads_mother.g.vcf.idx`.

#### 2.2.3. Ripetere il processo sugli altri due campioni

Generate i GVCF per i due campioni rimanenti eseguendo i comandi sottostanti, uno dopo l'altro.

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

Una volta completato, dovreste avere tre file che terminano con `.g.vcf` nella vostra directory corrente (uno per campione) e i rispettivi file indice che terminano con `.g.vcf.idx`.

Ma non uscite dal container!
Utilizzeremo lo stesso container nel passaggio successivo.

### 2.3. Eseguire il joint genotyping

Ora che abbiamo tutti i GVCF, possiamo provare l'approccio di joint genotyping per generare varianti per una coorte di campioni.
È un metodo in due fasi che consiste nel combinare i dati di tutti i GVCF in un data store, quindi eseguire l'analisi di joint genotyping vera e propria per generare il VCF finale delle varianti identificate congiuntamente.

#### 2.3.1. Combinare tutti i GVCF per campione

Questo primo passaggio utilizza un altro strumento GATK, chiamato GenomicsDBImport, per combinare i dati di tutti i GVCF in un data store GenomicsDB.

```bash
gatk GenomicsDBImport \
    -V reads_mother.g.vcf \
    -V reads_father.g.vcf \
    -V reads_son.g.vcf \
    -L /data/ref/intervals.bed \
    --genomicsdb-workspace-path family_trio_gdb
```

??? success "Output del comando"

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

L'output di questo passaggio è effettivamente una directory contenente un insieme di ulteriori directory nidificate che contengono i dati di varianti combinati sotto forma di più file diversi.
Potete esplorarlo ma vedrete rapidamente che questo formato di data store non è pensato per essere letto direttamente dagli esseri umani.

!!! note

    GATK include strumenti che rendono possibile ispezionare ed estrarre dati di varianti dal data store secondo necessità.

#### 2.3.2. Eseguire l'analisi di joint genotyping vera e propria

Questo secondo passaggio utilizza un altro strumento GATK, chiamato GenotypeGVCFs, per ricalcolare le statistiche delle varianti e i genotipi individuali alla luce dei dati disponibili in tutti i campioni della coorte.

```bash
gatk GenotypeGVCFs \
    -R /data/ref/ref.fasta \
    -V gendb://family_trio_gdb \
    -O family_trio.vcf
```

??? success "Output del comando"

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

Questo crea il file di output VCF `family_trio.vcf` nella directory di lavoro corrente nel container.
È un altro file ragionevolmente piccolo quindi potete usare `cat` per visualizzarne il contenuto e scorrere verso l'alto per trovare le prime righe di varianti.

??? abstract "File contents"

    ```console title="family_trio.vcf" linenums="40"
    #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  reads_father    reads_mother    reads_son
    20_10037292_10066351    3480    .       C       CT      1625.89 .       AC=5;AF=0.833;AN=6;BaseQRankSum=0.220;DP=85;ExcessHet=0.0000;FS=2.476;MLEAC=5;MLEAF=0.833;MQ=60.00;MQRankSum=0.00;QD=21.68;ReadPosRankSum=-1.147e+00;SOR=0.487    GT:AD:DP:GQ:PL  0/1:15,16:31:99:367,0,375       1/1:0,18:18:54:517,54,0 1/1:0,26:26:78:756,78,0
    20_10037292_10066351    3520    .       AT      A       1678.89 .       AC=5;AF=0.833;AN=6;BaseQRankSum=1.03;DP=80;ExcessHet=0.0000;FS=2.290;MLEAC=5;MLEAF=0.833;MQ=60.00;MQRankSum=0.00;QD=22.39;ReadPosRankSum=0.701;SOR=0.730 GT:AD:DP:GQ:PL   0/1:18,13:31:99:296,0,424       1/1:0,18:18:54:623,54,0 1/1:0,26:26:78:774,78,0
    20_10037292_10066351    3529    .       T       A       154.29  .       AC=1;AF=0.167;AN=6;BaseQRankSum=-5.440e-01;DP=104;ExcessHet=0.0000;FS=1.871;MLEAC=1;MLEAF=0.167;MQ=60.00;MQRankSum=0.00;QD=7.71;ReadPosRankSum=-1.158e+00;SOR=1.034       GT:AD:DP:GQ:PL  0/0:44,0:44:99:0,112,1347       0/1:12,8:20:99:163,0,328        0/0:39,0:39:99:0,105,1194
    ```

Questo appare simile al VCF che abbiamo generato in precedenza, tranne che questa volta abbiamo informazioni a livello di genotipo per tutti e tre i campioni.
Le ultime tre colonne nel file sono i blocchi di genotipo per i campioni, elencati in ordine alfabetico.

Se osserviamo i genotipi identificati per il nostro trio familiare di test per la primissima variante, vediamo che il padre è eterozigote-variante (`0/1`), e la madre e il figlio sono entrambi omozigoti-variante (`1/1`).

Questa è in definitiva l'informazione che stiamo cercando di estrarre dal dataset!

#### 2.3.3. Uscire dal container GATK

Per uscire dal container, digitate `exit`.

```bash
exit
```

Il vostro prompt dovrebbe essere tornato alla normalità.
Questo conclude il test manuale dei comandi di variant calling.

---

### Takeaway

Sapete come testare i comandi di indicizzazione Samtools e variant calling GATK nei rispettivi container, incluso come generare i GVCF ed eseguire il joint genotyping su più campioni.

### Cosa c'è dopo?

Imparate come avvolgere quegli stessi comandi in flussi di lavoro che utilizzano container per eseguire il lavoro.
