# Parte 1: Panoramica del metodo e test manuali

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduzione assistita da IA - [scopri di più e suggerisci miglioramenti](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Esistono molteplici metodi validi per elaborare e analizzare dati RNAseq in bulk.
Per questo corso, seguiamo il metodo descritto [qui](https://www.bioinformatics.babraham.ac.uk/training/RNASeq_Course/Analysing%20RNA-Seq%20data%20Exercise.pdf) dai Dott. Simon Andrews e Laura Biggins presso il [Babraham Institute](https://www.babraham.ac.uk/).

Il nostro obiettivo è sviluppare un workflow che implementi le seguenti fasi di elaborazione: eseguire il controllo qualità iniziale sulle read in un campione RNAseq in bulk, rimuovere le sequenze degli adapter dalle read, allineare le read a un genoma di riferimento e produrre un report completo di controllo qualità (QC).

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/rnaseq/img/preprocess.svg"
</figure>

- **FASTQC:** Eseguire il QC sui dati delle read prima del trimming utilizzando FastQC
- **TRIM_GALORE:** Rimuovere le sequenze degli adapter ed eseguire il QC dopo il trimming utilizzando Trim Galore (che include Cutadapt e FastQC)
- **HISAT2_ALIGN:** Allineare le read al genoma di riferimento utilizzando Hisat2
- **MULTIQC:** Generare un report QC completo utilizzando MultiQC

Tuttavia, prima di iniziare a scrivere qualsiasi codice del workflow, proveremo i comandi manualmente su alcuni dati di test.
Gli strumenti necessari non sono installati nell'ambiente GitHub Codespaces, quindi li utilizzeremo tramite container (vedere [Hello Containers](../../hello_nextflow/05_hello_containers.md)).

!!! note "Nota"

     Verificare di trovarsi nella directory `nf4-science/rnaseq`. L'ultima parte del percorso mostrato quando si digita `pwd` dovrebbe essere `rnaseq`.

---

## 1. QC iniziale e rimozione degli adapter

Otterremo un'immagine container che ha sia `fastqc` che `trim_galore` installati, la avvieremo in modo interattivo ed eseguiremo i comandi di trimming e QC su uno dei file di dati di esempio.

### 1.1. Ottenere il container

```bash
docker pull community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18
```

Questo produce il seguente output della console mentre il sistema scarica l'immagine:

??? success "Output del comando"

    ```console
    0.6.10--1bf8ca4e1967cd18: Pulling from library/trim-galore
    dafa2b0c44d2: Pull complete
    dec6b097362e: Pull complete
    f88da01cff0b: Pull complete
    4f4fb700ef54: Pull complete
    92dc97a3ef36: Pull complete
    403f74b0f85e: Pull complete
    10b8c00c10a5: Pull complete
    17dc7ea432cc: Pull complete
    bb36d6c3110d: Pull complete
    0ea1a16bbe82: Pull complete
    030a47592a0a: Pull complete
    32ec762be2d0: Pull complete
    d2cb90387285: Pull complete
    Digest: sha256:4f00e7b2a09f3c8d8a9ce955120e177152fb1e56f63a2a6e186088b1250d9907
    Status: Downloaded newer image for community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18
    community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18
    ```

### 1.2. Avviare il container in modo interattivo

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18
```

<!--
??? success "Output del comando"

    ```console

    ```
-->

Il prompt cambierà in qualcosa come `(base) root@b645838b3314:/tmp#`, il che indica che ci si trova all'interno del container.

La parte `-v ./data:/data` del comando consentirà di accedere ai contenuti della directory `data/` dall'interno del container.

```bash
ls /data/reads
```

??? success "Output del comando"

    ```console
    ENCSR000COQ1_1.fastq.gz  ENCSR000COQ2_2.fastq.gz  ENCSR000COR2_1.fastq.gz  ENCSR000CPO1_2.fastq.gz
    ENCSR000COQ1_2.fastq.gz  ENCSR000COR1_1.fastq.gz  ENCSR000COR2_2.fastq.gz  ENCSR000CPO2_1.fastq.gz
    ENCSR000COQ2_1.fastq.gz  ENCSR000COR1_2.fastq.gz  ENCSR000CPO1_1.fastq.gz  ENCSR000CPO2_2.fastq.gzO
    ```

### 1.3. Eseguire il primo comando `fastqc`

Eseguiamo `fastqc` per raccogliere le metriche di controllo qualità sui dati delle read.

```bash
fastqc /data/reads/ENCSR000COQ1_1.fastq.gz
```

??? success "Output del comando"

    ```console
    application/gzip
    Started analysis of ENCSR000COQ1_1.fastq.gz
    Approx 5% complete for ENCSR000COQ1_1.fastq.gz
    Approx 10% complete for ENCSR000COQ1_1.fastq.gz
    Approx 15% complete for ENCSR000COQ1_1.fastq.gz
    Approx 20% complete for ENCSR000COQ1_1.fastq.gz
    Approx 25% complete for ENCSR000COQ1_1.fastq.gz
    Approx 30% complete for ENCSR000COQ1_1.fastq.gz
    Approx 35% complete for ENCSR000COQ1_1.fastq.gz
    Approx 40% complete for ENCSR000COQ1_1.fastq.gz
    Approx 45% complete for ENCSR000COQ1_1.fastq.gz
    Approx 50% complete for ENCSR000COQ1_1.fastq.gz
    Approx 55% complete for ENCSR000COQ1_1.fastq.gz
    Approx 60% complete for ENCSR000COQ1_1.fastq.gz
    Approx 65% complete for ENCSR000COQ1_1.fastq.gz
    Approx 70% complete for ENCSR000COQ1_1.fastq.gz
    Approx 75% complete for ENCSR000COQ1_1.fastq.gz
    Approx 80% complete for ENCSR000COQ1_1.fastq.gz
    Approx 85% complete for ENCSR000COQ1_1.fastq.gz
    Approx 90% complete for ENCSR000COQ1_1.fastq.gz
    Approx 95% complete for ENCSR000COQ1_1.fastq.gz
    Analysis complete for ENCSR000COQ1_1.fastq.gz
    ```

Questo dovrebbe essere eseguito molto rapidamente.
È possibile trovare i file di output nella stessa directory dei dati originali:

```bash
ls /data/reads/ENCSR000COQ1_1_fastqc*
```

<!-- switch to tree -->

```console title="Output"
/data/reads/ENCSR000COQ1_1_fastqc.html  /data/reads/ENCSR000COQ1_1_fastqc.zip
```

### 1.4. Rimuovere le sequenze degli adapter con `trim_galore`

Ora eseguiamo `trim_galore`, che include Cutadapt e FastQC, per rimuovere le sequenze degli adapter e raccogliere le metriche QC post-trimming.

```bash
trim_galore --fastqc /data/reads/ENCSR000COQ1_1.fastq.gz
```

Il flag `--fastqc` fa sì che il comando esegua automaticamente una fase di raccolta QC dopo il completamento del trimming.

_L'output è molto dettagliato, quindi quanto segue è abbreviato._

??? success "Output del comando"

    ```console
    Multicore support not enabled. Proceeding with single-core trimming.
    Path to Cutadapt set as: 'cutadapt' (default)
    Cutadapt seems to be working fine (tested command 'cutadapt --version')
    Cutadapt version: 4.9
    single-core operation.
    igzip command line interface 2.31.0
    igzip detected. Using igzip for decompressing

    <...>

    Analysis complete for ENCSR000COQ1_1_trimmed.fq.gz
    ```

È possibile trovare i file di output nella directory di lavoro:

```bash
ls ENCSR000COQ1_1*
```

```console title="Output"
ENCSR000COQ1_1.fastq.gz_trimming_report.txt  ENCSR000COQ1_1_trimmed_fastqc.html
ENCSR000COQ1_1_trimmed.fq.gz                 ENCSR000COQ1_1_trimmed_fastqc.zip
```

### 1.5. Spostare i file di output nel filesystem esterno al container

Tutto ciò che rimane all'interno del container sarà inaccessibile per il lavoro futuro, quindi spostiamo questi file in una nuova directory.

```bash
mkdir /data/trimmed
mv ENCSR000COQ1_1* /data/trimmed
```

### 1.6. Uscire dal container

```bash
exit
```

---

## 2. Allineare le read al genoma di riferimento

Otterremo un'immagine container che ha `hisat2` installato, la avvieremo in modo interattivo ed eseguiremo il comando di allineamento per allineare i dati RNAseq a un genoma di riferimento.

### 2.1. Ottenere il container `hisat2`

```bash
docker pull community.wave.seqera.io/library/hisat2_samtools:5e49f68a37dc010e
```

??? success "Output del comando"

    ```console
    Unable to find image 'community.wave.seqera.io/library/hisat2_samtools:5e49f68a37dc010e' locally
    5e49f68a37dc010e: Pulling from library/hisat2_samtools
    dafa2b0c44d2: Already exists
    dec6b097362e: Already exists
    f88da01cff0b: Already exists
    4f4fb700ef54: Already exists
    92dc97a3ef36: Already exists
    403f74b0f85e: Already exists
    10b8c00c10a5: Already exists
    17dc7ea432cc: Already exists
    bb36d6c3110d: Already exists
    0ea1a16bbe82: Already exists
    030a47592a0a: Already exists
    e74ed5dd390b: Pull complete
    abfcf0185e51: Pull complete
    Digest: sha256:29d8e1a3172a2bdde7be813f7ebec22d331388194a7c0de872b4ccca4bed8f45
    Status: Downloaded newer image for community.wave.seqera.io/library/hisat2_samtools:5e49f68a37dc010e
    ```

### 2.2. Avviare il container `hisat2` in modo interattivo

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/hisat2_samtools:5e49f68a37dc010e
```

Il comando è lo stesso di prima, con l'URI del container pertinente sostituito.

### 2.3. Creare i file di indice del genoma per Hisat2

Hisat2 richiede che il riferimento del genoma sia fornito in un formato molto specifico e non può semplicemente utilizzare il file FASTA `genome.fa` che forniamo, quindi coglieremo questa opportunità per creare le risorse pertinenti.

```bash
hisat2-build /data/genome.fa genome_index
```

L'output è molto dettagliato, quindi quanto segue è abbreviato:

<!-- TODO: switch to full output -->

??? success "Output del comando"

    ```console
    Settings:
      Output files: "genome_index.*.ht2"
    <...>
    Total time for call to driver() for forward index: 00:00:16
    ```

Questo crea molteplici file di indice del genoma, che è possibile trovare nella directory di lavoro.

```bash
ls genome_index.*
```

```console title="Output"
genome_index.1.ht2  genome_index.3.ht2  genome_index.5.ht2  genome_index.7.ht2
genome_index.2.ht2  genome_index.4.ht2  genome_index.6.ht2  genome_index.8.ht2
```

Li utilizzeremo tra un momento, ma prima generiamo un tarball compresso con questi file di indice del genoma; ne avremo bisogno più avanti e generarli non è tipicamente qualcosa che vogliamo fare come parte di un workflow.

```bash
tar -czvf /data/genome_index.tar.gz genome_index.*
```

Questo memorizza un tarball `genome_index.tar.gz` contenente i file di indice del genoma nella directory `data/` sul nostro filesystem, il che tornerà utile nella Parte 2 di questo corso.

### 2.4. Eseguire il comando `hisat2`

Ora possiamo eseguire il comando di allineamento, che esegue la fase di allineamento con `hisat2` e poi invia l'output a `samtools` per scrivere l'output come file BAM.

L'input dei dati delle read è il file `/data/trimmed/ENCSR000COQ1_1_trimmed.fq.gz` che abbiamo generato con `trim_galore` nella fase precedente.

```bash
hisat2 -x genome_index -U /data/trimmed/ENCSR000COQ1_1_trimmed.fq.gz \
    --new-summary --summary-file ENCSR000COQ1_1_trimmed.hisat2.log | \
    samtools view -bS -o ENCSR000COQ1_1_trimmed.bam
```

??? success "Output del comando"

    ```console
    HISAT2 summary stats:
            Total reads: 27816
                    Aligned 0 time: 1550 (5.57%)
                    Aligned 1 time: 25410 (91.35%)
                    Aligned >1 times: 856 (3.08%)
            Overall alignment rate: 94.43%
    ```

Questo viene eseguito quasi istantaneamente perché è un file di test molto piccolo.
Su scala reale potrebbe richiedere molto più tempo.

Ancora una volta è possibile trovare i file di output nella directory di lavoro:

```bash
ls ENCSR000COQ1_1*
```

```console title="Output"
ENCSR000COQ1_1_trimmed.bam  ENCSR000COQ1_1_trimmed.hisat2.log
```

### 2.5. Spostare i file di output nel filesystem esterno al container

```bash
mkdir /data/aligned
mv ENCSR000COQ1_1* /data/aligned
```

### 2.6. Uscire dal container

```bash
exit
```

---

## 3. Generare un report QC completo

Otterremo un'immagine container che ha `multiqc` installato, la avvieremo in modo interattivo ed eseguiremo un comando di generazione del report sui file di report FastQC prima/dopo.

### 3.1. Ottenere il container `multiqc`

```bash
docker pull community.wave.seqera.io/library/pip_multiqc:a3c26f6199d64b7c
```

??? success "Output del comando"

    ```console
    ad8f247edb55897c: Pulling from library/pip_multiqc
    dafa2b0c44d2: Already exists
    dec6b097362e: Already exists
    f88da01cff0b: Already exists
    4f4fb700ef54: Already exists
    92dc97a3ef36: Already exists
    403f74b0f85e: Already exists
    10b8c00c10a5: Already exists
    17dc7ea432cc: Already exists
    bb36d6c3110d: Already exists
    0ea1a16bbe82: Already exists
    030a47592a0a: Already exists
    3f229294c69a: Pull complete
    5a5ad47fd84c: Pull complete
    Digest: sha256:0ebb1d9605395a7df49ad0eb366b21f46afd96a5090376b0d8941cf5294a895a
    Status: Downloaded newer image for community.wave.seqera.io/library/pip_multiqc:a3c26f6199d64b7c
    community.wave.seqera.io/library/pip_multiqc:a3c26f6199d64b7c
    ```

### 3.2. Avviare il container `multiqc` in modo interattivo

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/pip_multiqc:a3c26f6199d64b7c
```

### 3.3. Eseguire il comando `multiqc`

```bash
multiqc /data/reads /data/trimmed /data/aligned -n ENCSR000COQ1_1_QC
```

??? success "Output del comando"

    ```console

    /// MultiQC 🔍 v1.27.1

          file_search | Search path: /data/reads
          file_search | Search path: /data/trimmed
          file_search | Search path: /data/aligned
            searching | ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 100% 20/20
                hisat2 | Found 1 reports
              cutadapt | Found 1 reports
                fastqc | Found 1 reports
        write_results | Data        : ENCSR000COQ1_1_QC_data
        write_results | Report      : ENCSR000COQ1_1_QC.html
              multiqc | MultiQC complete
    ```

MultiQC è in grado di cercare nelle directory i report QC compatibili e aggrega tutto ciò che trova.

Qui vediamo che lo strumento ha trovato tutti e tre i report QC che abbiamo generato: il QC iniziale eseguito con `fastqc`, il report post-trimming da `cutadapt` (prodotto tramite `trim_galore`) e il QC post-allineamento prodotto da `hisat2`.

I file di output sono ancora una volta nella directory di lavoro:

```bash
ls ENCSR000COQ1_1_QC*
```

```console title="Output"
ENCSR000COQ1_1_QC.html

ENCSR000COQ1_1_QC_data:
cutadapt_filtered_reads_plot.txt                     fastqc_top_overrepresented_sequences_table.txt
cutadapt_trimmed_sequences_plot_3_Counts.txt         hisat2_se_plot.txt
cutadapt_trimmed_sequences_plot_3_Obs_Exp.txt        multiqc.log
fastqc-status-check-heatmap.txt                      multiqc_citations.txt
fastqc_adapter_content_plot.txt                      multiqc_cutadapt.txt
fastqc_per_base_n_content_plot.txt                   multiqc_data.json
fastqc_per_base_sequence_quality_plot.txt            multiqc_fastqc.txt
fastqc_per_sequence_gc_content_plot_Counts.txt       multiqc_general_stats.txt
fastqc_per_sequence_gc_content_plot_Percentages.txt  multiqc_hisat2.txt
fastqc_per_sequence_quality_scores_plot.txt          multiqc_software_versions.txt
fastqc_sequence_counts_plot.txt                      multiqc_sources.txt
fastqc_sequence_duplication_levels_plot.txt
```

### 3.4. Spostare i file di output nel filesystem esterno al container

```bash
mkdir /data/final_qc
mv ENCSR000COQ1_1_QC** /data/final_qc
```

### 3.5. Uscire dal container

```bash
exit
```

---

### Takeaway

Sono stati testati tutti i comandi individuali in modo interattivo nei container pertinenti.

### Passi successivi

Imparare come inserire gli stessi comandi in un workflow multi-fase che utilizza container per eseguire il lavoro.
