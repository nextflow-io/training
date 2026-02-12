# Parte 3: Implementazione con campioni multipli paired-end

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduzione assistita da IA - [scopri di più e suggerisci miglioramenti](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

In precedenza, avete costruito una pipeline di variant calling per campione singolo che elaborava i dati di ciascun campione in modo indipendente.
In questa parte del corso, porteremo il nostro semplice workflow al livello successivo trasformandolo in un potente strumento di automazione batch per gestire un numero arbitrario di campioni.
E mentre ci siamo, lo aggiorneremo anche per accettare dati paired-end, che sono più comuni negli studi più recenti.

??? info "Come iniziare da questa sezione"

    Questa sezione del corso presuppone che abbiate completato la [Parte 1: Panoramica del metodo](./01_method.md), la [Parte 2: Implementazione con campione singolo](./02_single-sample.md) e che abbiate una pipeline `rnaseq.nf` funzionante con i file dei moduli compilati.

    Se non avete completato la Parte 2 o volete iniziare da zero per questa parte, potete usare la soluzione della Parte 2 come punto di partenza.
    Eseguite questi comandi dall'interno della directory `nf4-science/rnaseq/`:

    ```bash
    cp solutions/part2/rnaseq-2.nf rnaseq.nf
    cp solutions/part2/modules/fastqc.nf modules/
    cp solutions/part2/modules/trim_galore.nf modules/
    cp solutions/part2/modules/hisat2_align.nf modules/
    cp solutions/part2/nextflow.config .
    ```

    Questo vi fornisce un workflow completo di elaborazione per campione singolo.
    Potete verificare che funzioni correttamente:

    ```bash
    nextflow run rnaseq.nf -profile test
    ```

## Compito

In questa parte del corso, estenderemo il workflow per fare quanto segue:

1. Leggere le informazioni sui campioni da un samplesheet CSV
2. Eseguire QC per campione, trimming e allineamento su tutti i campioni in parallelo
3. Aggregare tutti i report QC in un report MultiQC completo

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/rnaseq/img/rnaseq-wf-03.svg"
</figure>

Questo automatizza i passaggi della seconda sezione della [Parte 1: Panoramica del metodo](./01_method.md#2-multi-sample-qc-aggregation), dove avete eseguito questi comandi manualmente nei loro container.

## Piano della lezione

Abbiamo suddiviso questo in tre fasi:

1. **Far accettare al workflow campioni di input multipli.**
   Questo copre il passaggio da un singolo percorso di file a un samplesheet CSV, l'analisi con `splitCsv()` e l'esecuzione di tutti i processi esistenti su campioni multipli.
2. **Aggiungere la generazione di report QC completi.**
   Questo introduce l'operatore `collect()` per aggregare gli output attraverso i campioni e aggiunge un processo MultiQC per produrre un report combinato.
3. **Passare a dati RNAseq paired-end.**
   Questo copre l'adattamento dei processi per input paired-end (usando tuple), la creazione di moduli paired-end e la configurazione di un profilo di test separato.

Questo implementa il metodo descritto nella [Parte 1: Panoramica del metodo](./01_method.md) (seconda sezione che copre il caso d'uso multi-campione) e si basa direttamente sul workflow prodotto dalla Parte 2.

!!! tip "Suggerimento"

     Assicuratevi di essere nella directory di lavoro corretta:
     `cd /workspaces/training/nf4-science/rnaseq`

---

## 1. Far accettare al workflow campioni di input multipli

Per eseguire su campioni multipli, dobbiamo cambiare il modo in cui gestiamo l'input: invece di fornire un singolo percorso di file, leggeremo le informazioni sui campioni da un file CSV.

Forniamo un file CSV contenente ID dei campioni e percorsi dei file FASTQ nella directory `data/`.

```csv title="data/single-end.csv" linenums="1"
sample_id,fastq_path
ENCSR000COQ1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ1_1.fastq.gz
ENCSR000COQ2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ2_1.fastq.gz
ENCSR000COR1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR1_1.fastq.gz
ENCSR000COR2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR2_1.fastq.gz
ENCSR000CPO1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO1_1.fastq.gz
ENCSR000CPO2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO2_1.fastq.gz
```

Questo file CSV include una riga di intestazione che nomina le colonne.

Si noti che questi sono ancora dati read single-end.

!!! warning "Avviso"

    I percorsi dei file nel CSV sono percorsi assoluti che devono corrispondere al vostro ambiente.
    Se non state eseguendo questo nell'ambiente di formazione che forniamo, dovrete aggiornare i percorsi per corrispondere al vostro sistema.

### 1.1. Modificare l'input primario in un CSV di percorsi di file nel profilo di test

Prima, dobbiamo aggiornare il profilo di test in `nextflow.config` per fornire il percorso del file CSV invece del singolo percorso FASTQ.

=== "Dopo"

    ```groovy title="nextflow.config" linenums="1" hl_lines="5"
    docker.enabled = true

    profiles {
        test {
            params.input = "${projectDir}/data/single-end.csv"
            params.hisat2_index_zip = "${projectDir}/data/genome_index.tar.gz"
        }
    }
    ```

=== "Prima"

    ```groovy title="nextflow.config" linenums="1"
    docker.enabled = true

    profiles {
        test {
            params.input = "${projectDir}/data/reads/ENCSR000COQ1_1.fastq.gz"
            params.hisat2_index_zip = "${projectDir}/data/genome_index.tar.gz"
        }
    }
    ```

Successivamente, dovremo aggiornare la creazione del canale per leggere da questo CSV.

### 1.2. Aggiornare la factory del canale per analizzare l'input CSV

Dobbiamo caricare il contenuto del file nel canale invece del solo percorso del file.

Possiamo farlo usando lo stesso pattern che abbiamo usato nella [Parte 2 di Hello Nextflow](../../hello_nextflow/02_hello_channels.md#42-use-the-splitcsv-operator-to-parse-the-file): applicando l'operatore [`splitCsv()`](https://nextflow.io/docs/latest/reference/operator.html#splitcsv) per analizzare il file, poi un'operazione `map` per estrarre il percorso del file FASTQ da ogni riga.

=== "Dopo"

    ```groovy title="rnaseq.nf" linenums="19" hl_lines="4-7"
    workflow {

        main:
        // Crea canale di input dal contenuto di un file CSV
        read_ch = channel.fromPath(params.input)
            .splitCsv(header: true)
            .map { row -> file(row.fastq_path) }
    ```

=== "Prima"

    ```groovy title="rnaseq.nf" linenums="19"
    workflow {

        main:
        // Crea canale di input da un percorso di file
        read_ch = channel.fromPath(params.input)
    ```

Una cosa nuova rispetto a quanto incontrato nel corso Hello Nextflow è che questo CSV ha una riga di intestazione, quindi aggiungiamo `#!groovy header: true` alla chiamata `splitCsv()`.
Questo ci permette di riferirci alle colonne per nome nell'operazione `map`: `#!groovy row.fastq_path` estrae il percorso del file dalla colonna `fastq_path` di ogni riga.

La gestione dell'input è aggiornata e il workflow è pronto per essere testato.

### 1.3. Eseguire il workflow

Il workflow ora legge le informazioni sui campioni da un file CSV ed elabora tutti i campioni in parallelo.

```bash
nextflow run rnaseq.nf -profile test
```

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 24.10.0

    Launching `rnaseq.nf` [golden_curry] DSL2 - revision: 2a5ba5be1e

    executor >  local (18)
    [07/3ff9c5] FASTQC (6)       [100%] 6 of 6 ✔
    [cc/16859f] TRIM_GALORE (6)  [100%] 6 of 6 ✔
    [68/4c27b5] HISAT2_ALIGN (6) [100%] 6 of 6 ✔
    ```

Questa volta ogni passaggio viene eseguito 6 volte, una volta per ciascuno dei 6 campioni nel file CSV.

Questo è tutto ciò che è servito per far eseguire il workflow su file multipli.
Nextflow gestisce tutto il parallelismo per noi.

### Takeaway

Sapete come passare da un input a file singolo a un input multi-campione basato su CSV che Nextflow elabora in parallelo.

### Cosa seguirà?

Aggiungere un passaggio di aggregazione dei report QC che combina le metriche di tutti i campioni.

---

## 2. Aggregare le metriche QC di pre-elaborazione in un singolo report MultiQC

Tutto questo produce molti report QC, e non vogliamo dover esaminare i singoli report.
Questo è il punto perfetto per inserire un passaggio di aggregazione dei report MultiQC.

Ricordate il comando `multiqc` dalla [Parte 1](01_method.md):

```bash
multiqc . -n <output_name>.html
```

Il comando scansiona la directory corrente per file di output QC riconosciuti e li aggrega in un singolo report HTML.
L'URI del container era `community.wave.seqera.io/library/pip_multiqc:a3c26f6199d64b7c`.

Dobbiamo configurare un parametro aggiuntivo, preparare gli input, scrivere il processo, collegarlo e aggiornare la gestione degli output.

### 2.1. Configurare gli input

Il processo MultiQC necessita di un parametro per il nome del report e degli output QC raccolti da tutti i passaggi precedenti raggruppati insieme.

#### 2.1.1. Aggiungere un parametro `report_id`

Aggiungere un parametro per nominare il report di output.

=== "Dopo"

    ```groovy title="rnaseq.nf" linenums="11" hl_lines="8-9"
    params {
        // Input primario
        input: Path

        // Reference genome archive
        hisat2_index_zip: Path

        // Report ID
        report_id: String
    }
    ```

=== "Prima"

    ```groovy title="rnaseq.nf" linenums="11"
    params {
        // Input primario
        input: Path

        // Reference genome archive
        hisat2_index_zip: Path
    }
    ```

Aggiungere il valore predefinito del report ID al profilo di test:

=== "Dopo"

    ```groovy title="nextflow.config" linenums="1" hl_lines="7"
    docker.enabled = true

    profiles {
        test {
            params.input = "${projectDir}/data/single-end.csv"
            params.hisat2_index_zip = "${projectDir}/data/genome_index.tar.gz"
            params.report_id = "all_single-end"
        }
    }
    ```

=== "Prima"

    ```groovy title="nextflow.config" linenums="1"
    docker.enabled = true

    profiles {
        test {
            params.input = "${projectDir}/data/single-end.csv"
            params.hisat2_index_zip = "${projectDir}/data/genome_index.tar.gz"
        }
    }
    ```

Successivamente, dovremo preparare gli input per il processo MultiQC.

#### 2.1.2. Raccogliere e combinare gli output QC dai passaggi precedenti

Dobbiamo dare al processo `MULTIQC` tutti gli output relativi al QC dai passaggi precedenti raggruppati insieme.

Per questo, usiamo l'operatore `.mix()`, che aggrega più canali in uno solo.
Partiamo da `channel.empty()` e mescoliamo tutti i canali di output che vogliamo combinare.
Questo è più pulito che concatenare `.mix()` direttamente su uno dei canali di output, perché tratta tutti gli input in modo simmetrico.

Nel nostro workflow, gli output relativi al QC da aggregare sono:

- `FASTQC.out.zip`
- `FASTQC.out.html`
- `TRIM_GALORE.out.trimming_reports`
- `TRIM_GALORE.out.fastqc_reports`
- `HISAT2_ALIGN.out.log`

Li mescoliamo in un singolo canale, poi usiamo `.collect()` per aggregare i report attraverso tutti i campioni in una singola lista.

Aggiungere queste righe al corpo del workflow dopo la chiamata `HISAT2_ALIGN`:

=== "Dopo"

    ```groovy title="rnaseq.nf" linenums="38" hl_lines="4-13"
        // Allineamento a un genoma di riferimento
        HISAT2_ALIGN(TRIM_GALORE.out.trimmed_reads, file(params.hisat2_index_zip))

        // Generazione di report QC completi
        multiqc_files_ch = channel.empty().mix(
            FASTQC.out.zip,
            FASTQC.out.html,
            TRIM_GALORE.out.trimming_reports,
            TRIM_GALORE.out.fastqc_reports,
            HISAT2_ALIGN.out.log,
        )
        multiqc_files_list = multiqc_files_ch.collect()
    ```

=== "Prima"

    ```groovy title="rnaseq.nf" linenums="38"
        // Allineamento a un genoma di riferimento
        HISAT2_ALIGN(TRIM_GALORE.out.trimmed_reads, file(params.hisat2_index_zip))
    ```

L'uso di variabili intermedie rende ogni passaggio chiaro: `multiqc_files_ch` contiene tutti i singoli file QC mescolati in un canale, e `multiqc_files_list` è il bundle raccolto pronto per essere passato a MultiQC.

### 2.2. Scrivere il processo di aggregazione QC e chiamarlo nel workflow

Come prima, dobbiamo compilare la definizione del processo, importare il modulo e aggiungere la chiamata al processo.

#### 2.2.1. Compilare il modulo per il processo di aggregazione QC

Aprire `modules/multiqc.nf` ed esaminare lo schema della definizione del processo.

Procedete e compilate la definizione del processo da soli usando le informazioni fornite sopra, poi verificate il vostro lavoro confrontandolo con la soluzione nella scheda "Dopo" qui sotto.

=== "Prima"

    ```groovy title="modules/multiqc.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
     * Aggrega report QC con MultiQC
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

=== "Dopo"

    ```groovy title="modules/multiqc.nf" linenums="1" hl_lines="8 11 12 15 16 20"
    #!/usr/bin/env nextflow

    /*
     * Aggrega report QC con MultiQC
     */
    process MULTIQC {

        container "community.wave.seqera.io/library/pip_multiqc:a3c26f6199d64b7c"

        input:
        path '*'
        val output_name

        output:
        path "${output_name}.html", emit: report
        path "${output_name}_data", emit: data

        script:
        """
        multiqc . -n ${output_name}.html
        """
    }
    ```

Questo processo usa `#!groovy path '*'` come qualificatore di input per i file QC.
Il carattere jolly `'*'` dice a Nextflow di mettere in stage tutti i file raccolti nella directory di lavoro senza richiedere nomi specifici.
L'input `val output_name` è una stringa che controlla il nome del file del report.

Il comando `multiqc .` scansiona la directory corrente (dove sono tutti i file QC in stage) e genera il report.

Una volta completato questo, il processo è pronto per essere usato.

#### 2.2.2. Includere il modulo

Aggiungere l'istruzione di importazione a `rnaseq.nf`:

=== "Dopo"

    ```groovy title="rnaseq.nf" linenums="3" hl_lines="5"
    // Istruzioni INCLUDE dei moduli
    include { FASTQC } from './modules/fastqc.nf'
    include { TRIM_GALORE } from './modules/trim_galore.nf'
    include { HISAT2_ALIGN } from './modules/hisat2_align.nf'
    include { MULTIQC } from './modules/multiqc.nf'
    ```

=== "Prima"

    ```groovy title="rnaseq.nf" linenums="3"
    // Istruzioni INCLUDE dei moduli
    include { FASTQC } from './modules/fastqc.nf'
    include { TRIM_GALORE } from './modules/trim_galore.nf'
    include { HISAT2_ALIGN } from './modules/hisat2_align.nf'
    ```

Ora aggiungere la chiamata al processo nel workflow.

#### 2.2.3. Aggiungere la chiamata al processo

Passare i file QC raccolti e il report ID al processo `MULTIQC`:

=== "Dopo"

    ```groovy title="rnaseq.nf" linenums="49" hl_lines="2"
        multiqc_files_list = multiqc_files_ch.collect()
        MULTIQC(multiqc_files_list, params.report_id)
    ```

=== "Prima"

    ```groovy title="rnaseq.nf" linenums="49"
        multiqc_files_list = multiqc_files_ch.collect()
    ```

Il processo MultiQC è ora collegato nel workflow.

### 2.3. Aggiornare la gestione degli output

Dobbiamo aggiungere gli output MultiQC alla dichiarazione di pubblicazione e configurare dove vanno.

#### 2.3.1. Aggiungere target di pubblicazione per gli output MultiQC

Aggiungere gli output MultiQC alla sezione `publish:`:

=== "Dopo"

    ```groovy title="rnaseq.nf" linenums="52" hl_lines="9-10"
        publish:
        fastqc_zip = FASTQC.out.zip
        fastqc_html = FASTQC.out.html
        trimmed_reads = TRIM_GALORE.out.trimmed_reads
        trimming_reports = TRIM_GALORE.out.trimming_reports
        trimming_fastqc = TRIM_GALORE.out.fastqc_reports
        bam = HISAT2_ALIGN.out.bam
        align_log = HISAT2_ALIGN.out.log
        multiqc_report = MULTIQC.out.report
        multiqc_data = MULTIQC.out.data
    }
    ```

=== "Prima"

    ```groovy title="rnaseq.nf" linenums="52"
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

Successivamente, dovremo dire a Nextflow dove mettere questi output.

#### 2.3.2. Configurare i nuovi target di output

Aggiungere voci per i target MultiQC nel blocco `output {}`, pubblicandoli in una sottodirectory `multiqc/`:

=== "Dopo"

    ```groovy title="rnaseq.nf" linenums="82" hl_lines="7-12"
        align_log {
            path 'align'
        }
        multiqc_report {
            path 'multiqc'
        }
        multiqc_data {
            path 'multiqc'
        }
    }
    ```

=== "Prima"

    ```groovy title="rnaseq.nf" linenums="82"
        align_log {
            path 'align'
        }
    }
    ```

La configurazione degli output è completa.

### 2.4. Eseguire il workflow

Usiamo `-resume` in modo che i passaggi di elaborazione precedenti siano messi in cache e venga eseguito solo il nuovo passaggio MultiQC.

```bash
nextflow run rnaseq.nf -profile test -resume
```

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 24.10.0

    Launching `rnaseq.nf` [modest_pare] DSL2 - revision: fc724d3b49

    executor >  local (1)
    [07/3ff9c5] FASTQC (6)       [100%] 6 of 6, cached: 6 ✔
    [2c/8d8e1e] TRIM_GALORE (5)  [100%] 6 of 6, cached: 6 ✔
    [a4/7f9c44] HISAT2_ALIGN (6) [100%] 6 of 6, cached: 6 ✔
    [56/e1f102] MULTIQC          [100%] 1 of 1 ✔
    ```

Una singola chiamata a MULTIQC è stata aggiunta dopo le chiamate ai processi in cache.

Potete trovare gli output MultiQC nella directory dei risultati.

```bash
tree -L 2 results/multiqc
```

```console title="Output"
results/multiqc
├── all_single-end_data
│   ├── cutadapt_filtered_reads_plot.txt
│   ├── cutadapt_trimmed_sequences_plot_3_Counts.txt
│   ├── cutadapt_trimmed_sequences_plot_3_Obs_Exp.txt
│   ├── fastqc_adapter_content_plot.txt
│   ├── fastqc_overrepresented_sequences_plot.txt
│   ├── fastqc_per_base_n_content_plot.txt
│   ├── fastqc_per_base_sequence_quality_plot.txt
│   ├── fastqc_per_sequence_gc_content_plot_Counts.txt
│   ├── fastqc_per_sequence_gc_content_plot_Percentages.txt
│   ├── fastqc_per_sequence_quality_scores_plot.txt
│   ├── fastqc_sequence_counts_plot.txt
│   ├── fastqc_sequence_duplication_levels_plot.txt
│   ├── fastqc_sequence_length_distribution_plot.txt
│   ├── fastqc-status-check-heatmap.txt
│   ├── fastqc_top_overrepresented_sequences_table.txt
│   ├── hisat2_se_plot.txt
│   ├── multiqc_citations.txt
│   ├── multiqc_cutadapt.txt
│   ├── multiqc_data.json
│   ├── multiqc_fastqc.txt
│   ├── multiqc_general_stats.txt
│   ├── multiqc_hisat2.txt
│   ├── multiqc.log
│   ├── multiqc_software_versions.txt
│   └── multiqc_sources.txt
└── all_single-end.html
```

Quell'ultimo file `all_single-end.html` è il report aggregato completo, comodamente confezionato in un unico file HTML facile da consultare.

### Takeaway

Sapete come raccogliere output da canali multipli, raggrupparli con `.mix()` e `.collect()`, e passarli a un processo di aggregazione.

### Cosa seguirà?

Adattare il workflow per gestire dati RNAseq paired-end.

---

## 3. Abilitare l'elaborazione di dati RNAseq paired-end

Attualmente il nostro workflow può gestire solo dati RNAseq single-end.
È sempre più comune vedere dati RNAseq paired-end, quindi vogliamo essere in grado di gestirli.

Rendere il workflow completamente agnostico del tipo di dati richiederebbe l'uso di funzionalità del linguaggio Nextflow leggermente più avanzate, quindi non lo faremo qui, ma possiamo creare una versione per l'elaborazione paired-end per dimostrare cosa deve essere adattato.

### 3.1. Copiare il workflow e aggiornare gli input

Iniziamo copiando il file del workflow single-end e aggiornandolo per i dati paired-end.

#### 3.1.1. Copiare il file del workflow

Creare una copia del file del workflow da usare come punto di partenza per la versione paired-end.

```bash
cp rnaseq.nf rnaseq_pe.nf
```

Ora aggiornare i parametri e la gestione degli input nel nuovo file.

#### 3.1.2. Aggiungere un profilo di test paired-end

Forniamo un secondo file CSV contenente ID dei campioni e percorsi di file FASTQ paired nella directory `data/`.

```csv title="data/paired-end.csv" linenums="1"
sample_id,fastq_1,fastq_2
ENCSR000COQ1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ1_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ1_2.fastq.gz
ENCSR000COQ2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ2_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ2_2.fastq.gz
ENCSR000COR1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR1_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR1_2.fastq.gz
ENCSR000COR2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR2_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR2_2.fastq.gz
ENCSR000CPO1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO1_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO1_2.fastq.gz
ENCSR000CPO2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO2_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO2_2.fastq.gz
```

Aggiungere un profilo `test_pe` a `nextflow.config` che punta a questo file e usa un report ID paired-end.

=== "Dopo"

    ```groovy title="nextflow.config" linenums="1" hl_lines="9-13"
    docker.enabled = true

    profiles {
        test {
            params.input = "${projectDir}/data/single-end.csv"
            params.hisat2_index_zip = "${projectDir}/data/genome_index.tar.gz"
            params.report_id = "all_single-end"
        }
        test_pe {
            params.input = "${projectDir}/data/paired-end.csv"
            params.hisat2_index_zip = "${projectDir}/data/genome_index.tar.gz"
            params.report_id = "all_paired-end"
        }
    }
    ```

=== "Prima"

    ```groovy title="nextflow.config" linenums="1"
    docker.enabled = true

    profiles {
        test {
            params.input = "${projectDir}/data/single-end.csv"
            params.hisat2_index_zip = "${projectDir}/data/genome_index.tar.gz"
            params.report_id = "all_single-end"
        }
    }
    ```

Il profilo di test per i dati paired-end è pronto.

#### 3.1.3. Aggiornare la factory del canale

L'operatore `.map()` deve acquisire ora entrambi i percorsi dei file FASTQ e restituirli come lista.

=== "Dopo"

    ```groovy title="rnaseq_pe.nf" linenums="25" hl_lines="4"
        // Crea canale di input dal contenuto di un file CSV
        read_ch = channel.fromPath(params.input)
            .splitCsv(header: true)
            .map { row -> [file(row.fastq_1), file(row.fastq_2)] }
    ```

=== "Prima"

    ```groovy title="rnaseq_pe.nf" linenums="25"
        // Crea canale di input dal contenuto di un file CSV
        read_ch = channel.fromPath(params.input)
            .splitCsv(header: true)
            .map { row -> file(row.fastq_path) }
    ```

La gestione degli input è configurata per i dati paired-end.

### 3.2. Adattare il modulo FASTQC per i dati paired-end

Copiare il modulo per creare una versione paired-end:

```bash
cp modules/fastqc.nf modules/fastqc_pe.nf
```

L'input del processo FASTQC non ha bisogno di cambiare — quando Nextflow riceve una lista di due file, mette in stage entrambi e `reads` si espande in entrambi i nomi dei file.
L'unico cambiamento necessario è nel blocco di output: poiché ora otteniamo due report FastQC per campione, passiamo da pattern basati su `simpleName` a caratteri jolly.

=== "Dopo"

    ```groovy title="modules/fastqc_pe.nf" linenums="10" hl_lines="2 3"
        output:
        path "*_fastqc.zip", emit: zip
        path "*_fastqc.html", emit: html
    ```

=== "Prima"

    ```groovy title="modules/fastqc_pe.nf" linenums="10"
        output:
        path "${reads.simpleName}_fastqc.zip", emit: zip
        path "${reads.simpleName}_fastqc.html", emit: html
    ```

Questo generalizza il processo in modo da renderlo in grado di gestire sia dati single-end che paired-end.

Aggiornare l'importazione in `rnaseq_pe.nf` per usare la versione paired-end:

=== "Dopo"

    ```groovy title="rnaseq_pe.nf" linenums="4" hl_lines="1"
    include { FASTQC } from './modules/fastqc_pe.nf'
    ```

=== "Prima"

    ```groovy title="rnaseq_pe.nf" linenums="4"
    include { FASTQC } from './modules/fastqc.nf'
    ```

Il modulo FASTQC e la sua importazione sono aggiornati per i dati paired-end.

### 3.3. Adattare il modulo TRIM_GALORE per i dati paired-end

Copiare il modulo per creare una versione paired-end:

```bash
cp modules/trim_galore.nf modules/trim_galore_pe.nf
```

Questo modulo necessita di cambiamenti più sostanziali:

- L'input cambia da un singolo percorso a una tupla di due percorsi
- Il comando aggiunge il flag `--paired` e prende entrambi i file read
- L'output cambia per riflettere le convenzioni di denominazione paired-end di Trim Galore, producendo report FastQC separati per ciascun file read

=== "Dopo"

    ```groovy title="modules/trim_galore_pe.nf" linenums="8" hl_lines="2 5 7 8 12"
        input:
        tuple path(read1), path(read2)

        output:
        tuple path("*_val_1.fq.gz"), path("*_val_2.fq.gz"), emit: trimmed_reads
        path "*_trimming_report.txt", emit: trimming_reports
        path "*_val_1_fastqc.{zip,html}", emit: fastqc_reports_1
        path "*_val_2_fastqc.{zip,html}", emit: fastqc_reports_2

        script:
        """
        trim_galore --fastqc --paired ${read1} ${read2}
        """
    ```

=== "Prima"

    ```groovy title="modules/trim_galore_pe.nf" linenums="8"
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
    ```

Aggiornare l'importazione in `rnaseq_pe.nf`:

=== "Dopo"

    ```groovy title="rnaseq_pe.nf" linenums="5" hl_lines="1"
    include { TRIM_GALORE } from './modules/trim_galore_pe.nf'
    ```

=== "Prima"

    ```groovy title="rnaseq_pe.nf" linenums="5"
    include { TRIM_GALORE } from './modules/trim_galore.nf'
    ```

Il modulo TRIM_GALORE e la sua importazione sono aggiornati per i dati paired-end.

### 3.4. Adattare il modulo HISAT2_ALIGN per i dati paired-end

Copiare il modulo per creare una versione paired-end:

```bash
cp modules/hisat2_align.nf modules/hisat2_align_pe.nf
```

Questo modulo necessita di cambiamenti simili:

- L'input cambia da un singolo percorso a una tupla di due percorsi
- Il comando HISAT2 cambia da `-U` (unpaired) agli argomenti read `-1` e `-2` (paired)
- Tutti gli usi di `reads.simpleName` cambiano in `read1.simpleName` poiché ora facciamo riferimento a un membro specifico della coppia

=== "Dopo"

    ```groovy title="modules/hisat2_align_pe.nf" linenums="8" hl_lines="2 6 7 12 13 14"
        input:
        tuple path(read1), path(read2)
        path index_zip

        output:
        path "${read1.simpleName}.bam", emit: bam
        path "${read1.simpleName}.hisat2.log", emit: log

        script:
        """
        tar -xzvf ${index_zip}
        hisat2 -x ${index_zip.simpleName} -1 ${read1} -2 ${read2} \
            --new-summary --summary-file ${read1.simpleName}.hisat2.log | \
            samtools view -bS -o ${read1.simpleName}.bam
        """
    ```

=== "Prima"

    ```groovy title="modules/hisat2_align_pe.nf" linenums="8"
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
    ```

Aggiornare l'importazione in `rnaseq_pe.nf`:

=== "Dopo"

    ```groovy title="rnaseq_pe.nf" linenums="6" hl_lines="1"
    include { HISAT2_ALIGN } from './modules/hisat2_align_pe.nf'
    ```

=== "Prima"

    ```groovy title="rnaseq_pe.nf" linenums="6"
    include { HISAT2_ALIGN } from './modules/hisat2_align.nf'
    ```

Il modulo HISAT2_ALIGN e la sua importazione sono aggiornati per i dati paired-end.

### 3.5. Aggiornare l'aggregazione MultiQC per gli output paired-end

Il processo `TRIM_GALORE` paired-end ora produce due canali di report FastQC separati (`fastqc_reports_1` e `fastqc_reports_2`) invece di uno.
Aggiornare il blocco `.mix()` in `rnaseq_pe.nf` per includere entrambi:

=== "Dopo"

    ```groovy title="rnaseq_pe.nf" linenums="40" hl_lines="5 6"
        multiqc_files_ch = channel.empty().mix(
            FASTQC.out.zip,
            FASTQC.out.html,
            TRIM_GALORE.out.trimming_reports,
            TRIM_GALORE.out.fastqc_reports_1,
            TRIM_GALORE.out.fastqc_reports_2,
            HISAT2_ALIGN.out.log,
        )
    ```

=== "Prima"

    ```groovy title="rnaseq_pe.nf" linenums="40"
        multiqc_files_ch = channel.empty().mix(
            FASTQC.out.zip,
            FASTQC.out.html,
            TRIM_GALORE.out.trimming_reports,
            TRIM_GALORE.out.fastqc_reports,
            HISAT2_ALIGN.out.log,
        )
    ```

L'aggregazione MultiQC ora include entrambi i set di report FastQC paired-end.

### 3.6. Aggiornare la gestione degli output per gli output paired-end

La sezione `publish:` e il blocco `output {}` devono anche riflettere i due canali di report FastQC separati dal processo `TRIM_GALORE` paired-end.

Aggiornare la sezione `publish:` in `rnaseq_pe.nf`:

=== "Dopo"

    ```groovy title="rnaseq_pe.nf" linenums="52" hl_lines="6-7"
        publish:
        fastqc_zip = FASTQC.out.zip
        fastqc_html = FASTQC.out.html
        trimmed_reads = TRIM_GALORE.out.trimmed_reads
        trimming_reports = TRIM_GALORE.out.trimming_reports
        trimming_fastqc_1 = TRIM_GALORE.out.fastqc_reports_1
        trimming_fastqc_2 = TRIM_GALORE.out.fastqc_reports_2
        bam = HISAT2_ALIGN.out.bam
        align_log = HISAT2_ALIGN.out.log
        multiqc_report = MULTIQC.out.report
        multiqc_data = MULTIQC.out.data
    }
    ```

=== "Prima"

    ```groovy title="rnaseq_pe.nf" linenums="52"
        publish:
        fastqc_zip = FASTQC.out.zip
        fastqc_html = FASTQC.out.html
        trimmed_reads = TRIM_GALORE.out.trimmed_reads
        trimming_reports = TRIM_GALORE.out.trimming_reports
        trimming_fastqc = TRIM_GALORE.out.fastqc_reports
        bam = HISAT2_ALIGN.out.bam
        align_log = HISAT2_ALIGN.out.log
        multiqc_report = MULTIQC.out.report
        multiqc_data = MULTIQC.out.data
    }
    ```

Aggiornare le voci corrispondenti nel blocco `output {}`:

=== "Dopo"

    ```groovy title="rnaseq_pe.nf" linenums="77" hl_lines="4-9"
        trimming_reports {
            path 'trimming'
        }
        trimming_fastqc_1 {
            path 'trimming'
        }
        trimming_fastqc_2 {
            path 'trimming'
        }
        bam {
            path 'align'
        }
    ```

=== "Prima"

    ```groovy title="rnaseq_pe.nf" linenums="77"
        trimming_reports {
            path 'trimming'
        }
        trimming_fastqc {
            path 'trimming'
        }
        bam {
            path 'align'
        }
    ```

Il workflow paired-end è ora completamente aggiornato e pronto per essere eseguito.

### 3.7. Eseguire il workflow

Non usiamo `-resume` poiché questo non verrebbe messo in cache, e c'è il doppio dei dati da elaborare rispetto a prima, ma dovrebbe comunque completarsi in meno di un minuto.

```bash
nextflow run rnaseq_pe.nf -profile test_pe
```

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 24.10.0

    Launching `rnaseq_pe.nf` [reverent_kare] DSL2 - revision: 9c376cc219

    executor >  local (19)
    [c5/cbde15] FASTQC (5)       [100%] 6 of 6 ✔
    [e4/fa2784] TRIM_GALORE (5)  [100%] 6 of 6 ✔
    [3a/e23049] HISAT2_ALIGN (5) [100%] 6 of 6 ✔
    [e6/a3ccd9] MULTIQC          [100%] 1 of 1 ✔
    ```

Ora abbiamo due versioni leggermente divergenti del nostro workflow, una per dati read single-end e una per dati paired-end.
Il passo logico successivo sarebbe rendere il workflow in grado di accettare entrambi i tipi di dati al volo, che è al di fuori dello scopo di questo corso, ma potremmo affrontarlo in un seguito.

---

### Takeaway

Sapete come adattare un workflow per campione singolo per parallelizzare l'elaborazione di campioni multipli, generare un report QC completo e adattare il workflow per utilizzare dati read paired-end.

### Cosa seguirà?

Datevi una grande pacca sulla spalla! Avete completato il corso Nextflow per RNAseq.

Passate al [riepilogo finale del corso](./next_steps.md) per rivedere ciò che avete imparato e scoprire cosa viene dopo.
