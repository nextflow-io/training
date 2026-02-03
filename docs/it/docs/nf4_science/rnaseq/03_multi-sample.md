# Parte 3: Implementazione con campioni multipli paired-end

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduzione assistita da IA - [scopri di più e suggerisci miglioramenti](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

In questa parte finale del corso, porteremo il nostro semplice workflow al livello successivo trasformandolo in un potente strumento di automazione batch per gestire un numero arbitrario di campioni.
E mentre ci siamo, lo modificheremo anche per accettare dati paired-end, che sono più comuni negli studi più recenti.

Lo faremo in tre fasi:

1. Far accettare al workflow campioni di input multipli e parallelizzare l'esecuzione
2. Aggiungere la generazione di report QC completi
3. Passare a dati RNAseq paired-end

---

## 1. Far accettare al workflow campioni di input multipli e parallelizzare l'esecuzione

Dovremo cambiare il modo in cui gestiamo l'input.

### 1.1. Modificare l'input primario in un CSV di percorsi di file invece di un singolo file

Forniamo un file CSV contenente ID dei campioni e percorsi dei file FASTQ nella directory `data/`.
Questo file CSV include una riga di intestazione.
Si noti che i percorsi dei file FASTQ sono percorsi assoluti.

```csv title="data/single-end.csv" linenums="1"
sample_id,fastq_path
ENCSR000COQ1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ1_1.fastq.gz
ENCSR000COQ2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ2_1.fastq.gz
ENCSR000COR1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR1_1.fastq.gz
ENCSR000COR2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR2_1.fastq.gz
ENCSR000CPO1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO1_1.fastq.gz
ENCSR000CPO2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO2_1.fastq.gz
```

Rinominiamo il parametro di input primario in `input_csv` e modifichiamo il valore predefinito nel percorso del file `single-end.csv`.

```groovy title="rnaseq.nf" linenums="13"
params {
    // Input primario
    input_csv: Path = "data/single-end.csv"

    // Reference genome archive
    hisat2_index_zip: Path = "data/genome_index.tar.gz"
}
```

### 1.2. Aggiornare la factory del canale di input per gestire un CSV come input

Vogliamo caricare il contenuto del file nel canale invece del solo percorso del file, quindi usiamo l'operatore `.splitCsv()` per analizzare il formato CSV, poi l'operatore `.map()` per acquisire l'informazione specifica che vogliamo (il percorso del file FASTQ).

```groovy title="rnaseq.nf" linenums="16"
    // Crea canale di input dal contenuto di un file CSV
    read_ch = channel.fromPath(params.input_csv)
        .splitCsv(header:true)
        .map { row -> file(row.fastq_path) }
```

### 1.3. Eseguire il workflow per verificare che funzioni

```bash
nextflow run rnaseq.nf
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

Questa volta vediamo che ogni passaggio viene eseguito 6 volte, su ciascuno dei 6 file di dati che abbiamo fornito.

Questo è tutto ciò che è servito per far eseguire il workflow su file multipli!
Nextflow gestisce tutto il parallelismo per noi.

---

## 2. Aggregare le metriche QC di pre-elaborazione in un singolo report MultiQC

Tutto questo produce molti report QC, e non vogliamo dover esaminare i singoli report.
Questo è il punto perfetto per inserire un passaggio di aggregazione dei report MultiQC!

### 2.1. Creare un modulo per il processo di aggregazione QC

Creiamo un file di modulo chiamato `modules/multiqc.nf` per ospitare il processo `MULTIQC`:

```bash
touch modules/multiqc.nf
```

Aprire il file nell'editor di codice e copiarvi il seguente codice:

```groovy title="modules/multiqc.nf" linenums="1"
#!/usr/bin/env nextflow

process MULTIQC {

    container "community.wave.seqera.io/library/pip_multiqc:a3c26f6199d64b7c"
    publishDir "results/multiqc", mode: 'symlink'

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

### 2.2. Importare il modulo nel file del workflow

Aggiungere l'istruzione `include { MULTIQC } from './modules/multiqc.nf'` al file `rnaseq.nf`:

```groovy title="rnaseq.nf" linenums="3"
// Module INCLUDE statements
include { FASTQC } from './modules/fastqc.nf'
include { TRIM_GALORE } from './modules/trim_galore.nf'
include { HISAT2_ALIGN } from './modules/hisat2_align.nf'
include { MULTIQC } from './modules/multiqc.nf'
```

### 2.3. Aggiungere un parametro `report_id` e dargli un valore predefinito sensato

```groovy title="rnaseq.nf" linenums="9"
params {
    // Input primario
    input_csv: Path = "data/single-end.csv"

    // Reference genome archive
    hisat2_index_zip: Path = "data/genome_index.tar.gz"

    // Report ID
    report_id: String = "all_single-end"
}
```

### 2.4. Chiamare il processo sugli output dei passaggi precedenti

Dobbiamo dare al processo `MULTIQC` tutti gli output relativi al QC dai passaggi precedenti.

Per questo, useremo l'operatore `.mix()`, che aggrega più canali in uno solo.

Se avessimo quattro processi chiamati A, B, C e D con un semplice canale `.out` ciascuno, la sintassi sarebbe così: `A.out.mix( B.out, C.out, D.out )`. Come si potete vedere, lo si applica al primo dei canali che si vogliono combinare (non importa quale) e si aggiungono tutti gli altri, separati da virgole, nelle parentesi che seguono.

Nel caso del nostro workflow, abbiamo i seguenti output da aggregare:

- `FASTQC.out.zip`
- `FASTQC.out.html`
- `TRIM_GALORE.out.trimming_reports`
- `TRIM_GALORE.out.fastqc_reports`
- `HISAT2_ALIGN.out.log`

Quindi l'esempio di sintassi diventa:

```groovy title="Applicare .mix() nella chiamata MULTIQC"
        FASTQC.out.zip.mix(
        FASTQC.out.html,
        TRIM_GALORE.out.trimming_reports,
        TRIM_GALORE.out.fastqc_reports,
        HISAT2_ALIGN.out.log
        )
```

Questo raccoglierà i report QC per campione.
Ma poiché vogliamo aggregarli attraverso tutti i campioni, dobbiamo aggiungere l'operatore `collect()` per raccogliere i report di tutti i campioni in una singola chiamata a `MULTIQC`.
E dobbiamo anche dargli il parametro `report_id`.

Questo ci dà quanto segue:

```groovy title="La chiamata MULTIQC completa" linenums="33"
    // Comprehensive QC report generation
    MULTIQC(
        FASTQC.out.zip.mix(
        FASTQC.out.html,
        TRIM_GALORE.out.trimming_reports,
        TRIM_GALORE.out.fastqc_reports,
        HISAT2_ALIGN.out.log
        ).collect(),
        params.report_id
    )
```

Nel contesto del blocco del workflow completo, finisce per apparire così:

```groovy title="rnaseq.nf" linenums="18"
workflow {
    // Crea canale di input dal contenuto di un file CSV
    read_ch = channel.fromPath(params.input_csv)
        .splitCsv(header:true)
        .map { row -> file(row.fastq_path) }

    /// Initial quality control
    FASTQC(read_ch)

    // Adapter trimming and post-trimming QC
    TRIM_GALORE(read_ch)

    // Alignment to a reference genome
    HISAT2_ALIGN(TRIM_GALORE.out.trimmed_reads, file (params.hisat2_index_zip))

    // Comprehensive QC report generation
    MULTIQC(
        FASTQC.out.zip.mix(
        FASTQC.out.html,
        TRIM_GALORE.out.trimming_reports,
        TRIM_GALORE.out.fastqc_reports,
        HISAT2_ALIGN.out.log
        ).collect(),
        params.report_id
    )
}
```

### 2.5. Eseguire il workflow per verificare che funzioni

```bash
nextflow run rnaseq.nf -resume
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

Questa volta vediamo una singola chiamata a MULTIQC aggiunta dopo le chiamate ai processi in cache:

È possibile trovare gli output in `results/trimming` come specificato nel processo `TRIM_GALORE` dalla direttiva `publishDir`.

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

---

## 3. Abilitare l'elaborazione di dati RNAseq paired-end

Attualmente il nostro workflow può gestire solo dati RNAseq single-end.
È sempre più comune vedere dati RNAseq paired-end, quindi vogliamo essere in grado di gestirli.

Rendere il workflow completamente agnostico del tipo di dati richiederebbe l'uso di funzionalità del linguaggio Nextflow leggermente più avanzate, quindi non lo faremo qui, ma possiamo creare una versione per l'elaborazione paired-end per dimostrare cosa deve essere adattato.

### 3.1. Creare una copia del workflow chiamata `rnaseq_pe.nf`

```bash
cp rnaseq.nf rnaseq_pe.nf
```

### 3.2. Modificare l'`input_csv` predefinito per puntare ai dati paired-end

Forniamo un secondo file CSV contenente ID dei campioni e percorsi di file FASTQ paired nella directory `data/`

```csv title="data/paired-end.csv" linenums="1"
sample_id,fastq_1,fastq_2
ENCSR000COQ1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ1_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ1_2.fastq.gz
ENCSR000COQ2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ2_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ2_2.fastq.gz
ENCSR000COR1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR1_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR1_2.fastq.gz
ENCSR000COR2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR2_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR2_2.fastq.gz
ENCSR000CPO1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO1_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO1_2.fastq.gz
ENCSR000CPO2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO2_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO2_2.fastq.gz
```

Modifichiamo il valore predefinito di `input_csv` nel percorso del file `paired-end.csv`.

```groovy title="rnaseq_pe.nf" linenums="15"
params {
    // Input primario
    input_csv: Path = "data/paired-end.csv"

    // Reference genome archive
    hisat2_index_zip: Path = "data/genome_index.tar.gz"

    // Report ID
    report_id: String = "all_single-end"
}
```

### 3.3. Aggiornare la factory del canale

Dobbiamo dire all'operatore `.map()` di acquisire ora entrambi i percorsi dei file FASTQ.

Quindi `row -> file(row.fastq_path)` diventa `row -> [file(row.fastq_1), file(row.fastq_2)]`

```groovy title="rnaseq_pe.nf" linenums="19"
    // Crea canale di input dal contenuto di un file CSV
    read_ch = channel.fromPath(params.input_csv)
        .splitCsv(header:true)
        .map { row -> [file(row.fastq_1), file(row.fastq_2)] }
```

### 3.4. Creare una versione paired-end del processo FASTQC

Creiamo una copia del modulo in modo da avere entrambe le versioni a disposizione.

```bash
cp modules/fastqc.nf modules/fastqc_pe.nf
```

Aprire il nuovo file di modulo `fastqc_pe.nf` nell'editor di codice e apportare le seguenti modifiche al codice:

- Modificare `fastqc $reads` in `fastqc ${reads}` nel blocco `script` (riga 17) in modo che l'input `reads` venga decompresso, poiché ora è una tupla di due percorsi invece di un singolo percorso.
- Sostituire `${reads.simpleName}` con un carattere jolly (`*`) per evitare di dover gestire i file di output individualmente.

```groovy title="modules/fastqc_pe.nf" linenums="8"
    input:
    path reads

    output:
    path "*_fastqc.zip", emit: zip
    path "*_fastqc.html", emit: html

    script:
    """
    fastqc ${reads}
    """
```

Tecnicamente questo generalizza il processo `FASTQC` in modo da renderlo in grado di gestire sia dati RNAseq single-end che paired-end.

Infine, aggiornare l'istruzione di importazione del modulo per usare la versione paired-end del modulo.

```groovy title="rnaseq_pe.nf" linenums="4"
include { FASTQC } from './modules/fastqc_pe.nf'
```

### 3.5. Creare una versione paired-end del processo TRIM_GALORE

Creare una copia del modulo in modo da avere entrambe le versioni a disposizione.

```bash
cp modules/trim_galore.nf modules/trim_galore_pe.nf
```

Aprire il nuovo file di modulo `trim_galore_pe.nf` nell'editor di codice e apportare le seguenti modifiche al codice:

- Modificare la dichiarazione di input da `path reads` a `tuple path(read1), path(read2)`
- Aggiornare il comando nel blocco `script`, sostituendo `$reads` con `--paired ${read1} ${read2}`
- Aggiornare le dichiarazioni di output per riflettere i file aggiunti e le diverse convenzioni di denominazione, usando caratteri jolly per evitare di dover elencare tutto.

```groovy title="modules/trim_galore_pe.nf" linenums="8"
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

Infine, aggiornare l'istruzione di importazione del modulo per usare la versione paired-end del modulo.

```groovy title="rnaseq_pe.nf" linenums="5"
include { TRIM_GALORE } from './modules/trim_galore_pe.nf'
```

### 3.6. Aggiornare la chiamata al processo MULTIQC per aspettarsi due report da TRIM_GALORE

Il processo `TRIM_GALORE` ora produce un canale di output aggiuntivo, quindi dobbiamo fornirlo a MultiQC.

Sostituire `TRIM_GALORE.out.fastqc_reports,` con `TRIM_GALORE.out.fastqc_reports_1,` più `TRIM_GALORE.out.fastqc_reports_2,`:

```groovy title="rnaseq_pe.nf" linenums="33"
    // Comprehensive QC report generation
    MULTIQC(
        FASTQC.out.zip.mix(
        FASTQC.out.html,
        TRIM_GALORE.out.trimming_reports,
        TRIM_GALORE.out.fastqc_reports_1,
        TRIM_GALORE.out.fastqc_reports_2,
        HISAT2_ALIGN.out.log
        ).collect(),
        params.report_id
    )
```

Mentre siamo su MultiQC, aggiorniamo anche il valore predefinito del parametro `report_id` da `"all_single-end"` a `"all_paired-end"`.

```groovy title="rnaseq_pe.nf" linenums="9"
params {
    // Input primario
    input_csv: Path = "data/paired-end.csv"

    // Reference genome archive
    hisat2_index_zip: Path = "data/genome_index.tar.gz"

    // Report ID
    report_id: String = "all_paired-end"
}
```

### 3.7. Creare una versione paired-end del processo HISAT2_ALIGN

Creare una copia del modulo in modo da avere entrambe le versioni a disposizione.

```bash
cp modules/hisat2_align.nf modules/hisat2_align_pe.nf
```

Aprire il nuovo file di modulo `hisat2_align_pe.nf` nell'editor di codice e apportare le seguenti modifiche al codice:

- Modificare la dichiarazione di input da `path reads` a `tuple path(read1), path(read2)`
- Aggiornare il comando nel blocco `script`, sostituendo `-U $reads` con `-1 ${read1} -2 ${read2}`
- Sostituire tutte le istanze di `${reads.simpleName}` con `${read1.simpleName}` nel comando nel blocco `script` così come nelle dichiarazioni di output.

```groovy title="modules/hisat2_align_pe.nf" linenums="8"
    input:
    tuple path(read1), path(read2)
    path index_zip

    output:
    path "${read1.simpleName}.bam", emit: bam
    path "${read1.simpleName}.hisat2.log", emit: log

    script:
    """
    tar -xzvf $index_zip
    hisat2 -x ${index_zip.simpleName} -1 ${read1} -2 ${read2} \
        --new-summary --summary-file ${read1.simpleName}.hisat2.log | \
        samtools view -bS -o ${read1.simpleName}.bam
    """
```

Infine, aggiornare l'istruzione di importazione del modulo per usare la versione paired-end del modulo.

```groovy title="rnaseq_pe.nf" linenums="5"
include { HISAT2_ALIGN } from './modules/hisat2_align_pe.nf'
```

### 3.8. Eseguire il workflow per verificare che funzioni

Non usiamo `-resume` poiché questo non verrebbe messo in cache, e c'è il doppio dei dati da elaborare rispetto a prima, ma dovrebbe comunque completarsi in meno di un minuto.

```bash
nextflow run rnaseq_pe.nf
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

Ed è tutto! Ora abbiamo due versioni leggermente divergenti del nostro workflow, una per dati single-end read e una per dati paired-end.
Il passo logico successivo sarebbe rendere il workflow in grado di accettare entrambi i tipi di dati al volo, che è al di fuori dello scopo di questo corso, ma potremmo affrontarlo in un seguito.

---

### Takeaway

Sapete come adattare un workflow per campione singolo per parallelizzare l'elaborazione di campioni multipli, generare un report QC completo e adattare il workflow per utilizzare dati read paired-end se necessario.

### Cosa seguirà?

Congratulazioni, ha completato il mini-corso Nextflow Per RNAseq! Celebrate il vostro successo e prendetevi una meritata pausa!

Successivamente, vi chiediamo di completare un breve sondaggio sulla vostra esperienza con questo corso di formazione, poi vi porteremo a una pagina con collegamenti a ulteriori risorse di formazione e collegamenti utili.
