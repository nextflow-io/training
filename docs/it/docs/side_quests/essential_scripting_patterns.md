---
title: Pattern di scripting essenziali
description: Imparare tecniche di programmazione avanzate in Nextflow.
weight: 1200
---

# Pattern di scripting essenziali

Questa lezione ti guiderà attraverso i modelli di programmazione essenziali che sono fondamentali per creare flussi di lavoro Nextflow efficaci. Tratteremo modelli per gestire gli input di dati, trasformare i valori, controllare la logica del flusso di lavoro, allocare risorse dinamicamente e altro ancora.

## Obiettivi di apprendimento

- Comprendere le differenze tra i paradigmi di flusso dati e scripting
- Applicare closure, operatori ternari e altre tecniche di Groovy
- Padroneggiare le tecniche per manipolare i metadati ed estrarre informazioni dai file
- Utilizzare espressioni regolari ed elaborazione di stringhe per analizzare i nomi dei file
- Creare funzioni riutilizzabili per logiche complesse
- Implementare l'allocazione dinamica delle risorse e strategie di ripetizione
- Aggiungere logica condizionale per controllare l'esecuzione del flusso di lavoro
- Scrivere codice robusto utilizzando operatori di navigazione sicura ed Elvis
- Validare gli input con messaggi di errore chiari
- Utilizzare gestori di eventi per gestire il completamento del flusso di lavoro

## Prerequisiti

- Comprensione di base dei flussi di lavoro Nextflow
- Familiarità con la sintassi DSL2
- Conoscenza di base dei canali e dei processi
- Ambiente di sviluppo configurato con Nextflow v23.04.0 o successivo
- Accesso a Docker (o Conda) per i container software

## Come iniziare

Questo tutorial presuppone che tu abbia una conoscenza di base della sintassi di Nextflow e delle operazioni sui canali.

Per iniziare, naviga nella cartella `side-quests/essential_scripting_patterns`:

```bash
cd side-quests/essential_scripting_patterns
```

Il repository contiene diversi file:

- `main.nf`: Il flusso di lavoro principale
- `modules/fastp.nf`: Modulo per l'elaborazione della qualità con FASTP
- `modules/generate_report.nf`: Modulo per generare report
- `modules/trimgalore.nf`: Modulo per un programma di trimming alternativo
- `data/samples.csv`: Un CSV con informazioni sui campioni
- `data/sequences/*.fastq`: File fastq di esempio

Esaminiamo i file:

```bash
cat main.nf
```

Dovresti vedere un flusso di lavoro minimo che utilizza inclusioni di moduli:

```groovy title="main.nf" linenums="1"
include { FASTP } from './modules/fastp.nf'
include { GENERATE_REPORT } from './modules/generate_report.nf'

workflow {
    ch_samples = channel.fromPath("./data/samples.csv")
        .splitCsv(header: true)
        .map { row ->
            tuple(
                [id: row.sample_id],
                file(row.file_path)
            )
        }

    ch_fastp = FASTP(ch_samples)
    GENERATE_REPORT(ch_samples)
}
```

Anche i moduli sono semplici:

```bash
cat modules/fastp.nf
```

```groovy title="modules/fastp.nf" linenums="1"
process FASTP {
    container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}_trimmed.fastq.gz"), emit: fastq
    tuple val(meta), path("${meta.id}.fastp.json"), emit: json

    script:
    """
    fastp --in1 ${reads} --out1 ${meta.id}_trimmed.fastq.gz --json ${meta.id}.fastp.json --html ${meta.id}.fastp.html
    """
}
```

```bash
cat modules/generate_report.nf
```

```groovy title="modules/generate_report.nf" linenums="1"
process GENERATE_REPORT {
    container 'community.wave.seqera.io/library/bash:5.0'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}_report.txt"), emit: report

    script:
    """
    echo "Sample ID: ${meta.id}" > ${meta.id}_report.txt
    echo "File processed: ${reads}" >> ${meta.id}_report.txt
    echo "Process command: fastp --in1 ${reads} --out1 ${meta.id}_trimmed.fastq.gz" >> ${meta.id}_report.txt
    """
}
```

Infine, diamo un'occhiata ai dati di esempio:

```bash
cat data/samples.csv
```

```csv title="data/samples.csv"
sample_id,organism,tissue_type,sequencing_depth,quality_score,file_path
SAMPLE_001,human,liver,30000000,38.5,./data/sequences/SAMPLE_001_S1_L001_R1_001.fastq
SAMPLE_002,mouse,brain,25000000,35.2,./data/sequences/SAMPLE_002_S2_L001_R1_001.fastq
SAMPLE_003,human,kidney,45000000,42.1,./data/sequences/SAMPLE_003_S3_L001_R1_001.fastq
```

Eseguiamo prima il flusso di lavoro base per assicurarci che funzioni:

```bash
nextflow run main.nf
```

Dovresti vedere il flusso di lavoro eseguirsi e completarsi con successo.

Ora iniziamo a migliorare il flusso di lavoro con pattern di scripting avanzati!

---

## 1. Comprendere il Flusso di Dati vs lo Scripting

Nextflow utilizza **due paradigmi di programmazione distinti**:

1. **Flusso di dati**: L'orchestrazione di canali attraverso operatori (`.map`, `.filter`, `.branch`)
2. **Scripting**: Codice Groovy eseguito all'interno di closure o blocchi di script nei processi

Entrambi sono cruciali, ma funzionano in modo diverso, quindi questa prima sezione chiarirà le differenze.

### 1.1. Closure e operatori ternari

Un concetto fondamentale in Nextflow è la **closure** - un blocco di codice che può essere passato come oggetto ed eseguito in seguito. Le closure sono essenziali per gli operatori dei canali (`map`, `filter`, ecc.).

Prendiamo il nostro attuale operatore `.map` e miglioriamolo per estrarre più metadati:

```groovy
.map { row ->
    tuple(
        [id: row.sample_id],
        file(row.file_path)
    )
}
```

Il blocco `{ row -> ... }` è una closure che riceve un argomento `row`.

Ora, miglioriamolo per estrarre più metadati dal CSV:

=== "After"

    ```groovy title="main.nf" linenums="7" hl_lines="2-8"
        .map { row ->
            def sample_meta = [
                id: row.sample_id.toLowerCase(),
                organism: row.organism,
                tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                depth: row.sequencing_depth.toInteger(),
                quality: row.quality_score.toDouble()
            ]
            tuple(sample_meta, file(row.file_path))
        }
    ```

=== "Before"

    ```groovy title="main.nf" linenums="7" hl_lines="2-5"
        .map { row ->
            tuple(
                [id: row.sample_id],
                file(row.file_path)
            )
        }
    ```

Ora, modifichiamo `main.nf` per utilizzare questi metadati aggiuntivi. Prima, cambiamo il modo in cui generiamo il report. Etichetteremo i campioni di alta qualità come "alta priorità" usando un **operatore ternario**.

=== "After"

    ```groovy title="main.nf" linenums="7" hl_lines="9-10"
        .map { row ->
            def sample_meta = [
                id: row.sample_id.toLowerCase(),
                organism: row.organism,
                tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                depth: row.sequencing_depth.toInteger(),
                quality: row.quality_score.toDouble()
            ]
            def priority = sample_meta.quality > 40 ? 'high' : 'normal'
            tuple(sample_meta + [priority: priority], file(row.file_path))
        }
    ```

=== "Before"

    ```groovy title="main.nf" linenums="7" hl_lines="9"
        .map { row ->
            def sample_meta = [
                id: row.sample_id.toLowerCase(),
                organism: row.organism,
                tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                depth: row.sequencing_depth.toInteger(),
                quality: row.quality_score.toDouble()
            ]
            tuple(sample_meta, file(row.file_path))
        }
    ```

L'operatore ternario `condizione ? valore_se_vera : valore_se_falsa` è un modo conciso per scrivere un'istruzione if/else. Se `sample_meta.quality > 40`, allora `priority` sarà `'high'`; altrimenti, sarà `'normal'`.

Modifichiamo anche il processo `GENERATE_REPORT` per includere i metadati aggiuntivi:

=== "After"

    ```groovy title="modules/generate_report.nf" linenums="11" hl_lines="2-7"
    script:
    """
    echo "Sample ID: ${meta.id}" > ${meta.id}_report.txt
    echo "Organism: ${meta.organism}" >> ${meta.id}_report.txt
    echo "Tissue: ${meta.tissue}" >> ${meta.id}_report.txt
    echo "Sequencing depth: ${meta.depth}" >> ${meta.id}_report.txt
    echo "Quality score: ${meta.quality}" >> ${meta.id}_report.txt
    echo "Priority: ${meta.priority}" >> ${meta.id}_report.txt
    echo "File processed: ${reads}" >> ${meta.id}_report.txt
    echo "Process command: fastp --in1 ${reads} --out1 ${meta.id}_trimmed.fastq.gz" >> ${meta.id}_report.txt
    """
    ```

=== "Before"

    ```groovy title="modules/generate_report.nf" linenums="11" hl_lines="2-4"
    script:
    """
    echo "Sample ID: ${meta.id}" > ${meta.id}_report.txt
    echo "File processed: ${reads}" >> ${meta.id}_report.txt
    echo "Process command: fastp --in1 ${reads} --out1 ${meta.id}_trimmed.fastq.gz" >> ${meta.id}_report.txt
    """
    ```

Ora eseguiamo il flusso di lavoro aggiornato:

```bash
nextflow run main.nf
```

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [jolly_fourier] DSL2 - revision: d3e76a7fce

    executor >  local (6)
    [a3/6b9e80] process > FASTP (sample_003)           [100%] 3 of 3 ✔
    [29/54f4b6] process > GENERATE_REPORT (sample_001) [100%] 3 of 3 ✔
    ```

Possiamo verificare il risultato esaminando uno dei file di report:

```console
cat work/29/54f4b6b0eb90fed9e3a673b8e47629/sample_001_report.txt
```

Ora dovresti vedere tutti i metadati inclusi:

```
Sample ID: sample_001
Organism: human
Tissue: liver
Sequencing depth: 30000000
Quality score: 38.5
Priority: normal
File processed: SAMPLE_001_S1_L001_R1_001.fastq
Process command: fastp --in1 SAMPLE_001_S1_L001_R1_001.fastq --out1 sample_001_trimmed.fastq.gz
```

### 1.2. La collezione 'meta' vs 'reads'

In Nextflow, canali, tuple e collezioni (come mappe e liste) sono strutture dati fondamentali. Esiste una differenza importante tra **collezioni nel flusso dati** (canali/operatori) e **collezioni nei blocchi di script**.

Per dimostrare questo, modifichiamo il nostro flusso di lavoro per estrarre metadati dai nomi dei file FASTQ utilizzando espressioni regolari. È comune che i nomi dei file FASTQ seguano una convenzione come: `SAMPLE_001_S1_L001_R1_001.fastq`.

In questo formato:

- `SAMPLE_001`: ID del campione
- `S1`: Numero del campione
- `L001`: Numero della lane
- `R1`: Numero della read
- `001`: Chunk o frammento

Estraiamo questi metadati dal nome del file FASTQ:

=== "After"

    ```groovy title="main.nf" linenums="7" hl_lines="12-20"
        .map { row ->
            def sample_meta = [
                id: row.sample_id.toLowerCase(),
                organism: row.organism,
                tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                depth: row.sequencing_depth.toInteger(),
                quality: row.quality_score.toDouble()
            ]
            def priority = sample_meta.quality > 40 ? 'high' : 'normal'
            def fastq_path = file(row.file_path)

            // Extract metadata from filename using regex
            def m = (fastq_path.name =~ /^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$/)
            def file_meta = m ? [
                sample_num: m[0][2].toInteger(),
                lane: m[0][3],
                read: m[0][4],
                chunk: m[0][5]
            ] : [:]

            tuple(sample_meta + file_meta + [priority: priority], fastq_path)
        }
    ```

=== "Before"

    ```groovy title="main.nf" linenums="7" hl_lines="10"
        .map { row ->
            def sample_meta = [
                id: row.sample_id.toLowerCase(),
                organism: row.organism,
                tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                depth: row.sequencing_depth.toInteger(),
                quality: row.quality_score.toDouble()
            ]
            def priority = sample_meta.quality > 40 ? 'high' : 'normal'
            tuple(sample_meta + [priority: priority], file(row.file_path))
        }
    ```

C'è molto da analizzare qui!

1. `fastq_path.name` ottiene il nome del file (senza il percorso)
2. L'operatore `=~` è per il pattern matching regex
3. Il pattern `/^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$/` cattura i componenti
4. `m ? ... : ...` è un operatore ternario che gestisce il caso in cui il nome del file non corrisponde al pattern
5. `m[0][2]` accede al secondo gruppo di cattura (gli indici iniziano da 1 per i gruppi)
6. La funzione `toInteger()` converte la stringa catturata in un numero intero
7. `sample_meta + file_meta + [priority: priority]` unisce le tre mappe in una

Modifichiamo anche il processo `GENERATE_REPORT` per includere questi nuovi metadati:

=== "After"

    ```groovy title="modules/generate_report.nf" linenums="11" hl_lines="7-10"
    script:
    """
    echo "Sample ID: ${meta.id}" > ${meta.id}_report.txt
    echo "Organism: ${meta.organism}" >> ${meta.id}_report.txt
    echo "Tissue: ${meta.tissue}" >> ${meta.id}_report.txt
    echo "Sequencing depth: ${meta.depth}" >> ${meta.id}_report.txt
    echo "Quality score: ${meta.quality}" >> ${meta.id}_report.txt
    echo "Sample number: ${meta.sample_num}" >> ${meta.id}_report.txt
    echo "Lane: ${meta.lane}" >> ${meta.id}_report.txt
    echo "Read: ${meta.read}" >> ${meta.id}_report.txt
    echo "Chunk: ${meta.chunk}" >> ${meta.id}_report.txt
    echo "Priority: ${meta.priority}" >> ${meta.id}_report.txt
    echo "File processed: ${reads}" >> ${meta.id}_report.txt
    echo "Process command: fastp --in1 ${reads} --out1 ${meta.id}_trimmed.fastq.gz" >> ${meta.id}_report.txt
    """
    ```

=== "Before"

    ```groovy title="modules/generate_report.nf" linenums="11" hl_lines="7"
    script:
    """
    echo "Sample ID: ${meta.id}" > ${meta.id}_report.txt
    echo "Organism: ${meta.organism}" >> ${meta.id}_report.txt
    echo "Tissue: ${meta.tissue}" >> ${meta.id}_report.txt
    echo "Sequencing depth: ${meta.depth}" >> ${meta.id}_report.txt
    echo "Quality score: ${meta.quality}" >> ${meta.id}_report.txt
    echo "Priority: ${meta.priority}" >> ${meta.id}_report.txt
    echo "File processed: ${reads}" >> ${meta.id}_report.txt
    echo "Process command: fastp --in1 ${reads} --out1 ${meta.id}_trimmed.fastq.gz" >> ${meta.id}_report.txt
    """
    ```

Eseguiamo di nuovo il flusso di lavoro:

```bash
nextflow run main.nf
```

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [jolly_hawking] DSL2 - revision: cd0a5b0d29

    executor >  local (6)
    [b3/1cb89f] process > FASTP (sample_001)           [100%] 3 of 3 ✔
    [e6/c2f254] process > GENERATE_REPORT (sample_001) [100%] 3 of 3 ✔
    ```

Verifichiamo uno dei file di report aggiornati:

```console
cat work/e6/c2f2542ec23ee6e7f0aa9c66c12e30/sample_001_report.txt
```

Ora dovresti vedere i metadati aggiuntivi estratti dal nome del file:

```
Sample ID: sample_001
Organism: human
Tissue: liver
Sequencing depth: 30000000
Quality score: 38.5
Sample number: 1
Lane: 001
Read: R1
Chunk: 001
Priority: normal
File processed: SAMPLE_001_S1_L001_R1_001.fastq
Process command: fastp --in1 SAMPLE_001_S1_L001_R1_001.fastq --out1 sample_001_trimmed.fastq.gz
```

### 1.3. Operazioni di collezione nelle closure vs manipolazione dei canali

Un punto spesso fonte di confusione in Nextflow riguarda le operazioni sulle collezioni all'interno delle closure. Metodi come `collect()` funzionano in modo diverso a seconda del contesto:

1. Nei **canali Nextflow**: `channel.collect()` è un operatore che raccoglie tutti gli elementi di un canale in una singola lista
2. Nelle **liste e mappe Groovy**: `list.collect {...}` applica una funzione a ciascun elemento e restituisce una nuova lista

Modifichiamo il processo FASTP per illustrare questa differenza:

```groovy title="modules/fastp.nf" linenums="11" hl_lines="3-5"
script:
"""
# Dimostrazione di operazioni su collezioni nello script
def options = ["--in1", reads, "--out1", "${meta.id}_trimmed.fastq.gz", "--json", "${meta.id}.fastp.json", "--html", "${meta.id}.fastp.html"]
fastp ${options.join(' ')}
"""
```

Questo esempio usa `options.join(' ')` per unire gli elementi della lista in una singola stringa con spazi tra di loro.

Tuttavia, questo genererà un errore perché stiamo cercando di eseguire codice Groovy all'interno di un blocco di script bash. Cambiamolo per spostare la logica di collezione fuori dal blocco di script:

=== "After"

    ```groovy title="modules/fastp.nf" linenums="11" hl_lines="1-3"
    def options = ["--in1", reads, "--out1", "${meta.id}_trimmed.fastq.gz", "--json", "${meta.id}.fastp.json", "--html", "${meta.id}.fastp.html"]
    def cmd = "fastp ${options.join(' ')}"

    script:
    """
    $cmd
    """
    ```

=== "Before"

    ```groovy title="modules/fastp.nf" linenums="11" hl_lines="2-3"
    script:
    """
    fastp --in1 ${reads} --out1 ${meta.id}_trimmed.fastq.gz --json ${meta.id}.fastp.json --html ${meta.id}.fastp.html
    """
    ```

Esegui di nuovo il flusso di lavoro e dovrebbe funzionare! Questo dimostra l'uso dello scripting Groovy per la manipolazione delle collezioni prima di passare il comando al blocco di script.

```bash
nextflow run main.nf
```

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [compassionate_shaw] DSL2 - revision: 3471dc57d9

    executor >  local (6)
    [32/8e94af] process > FASTP (sample_001)           [100%] 3 of 3 ✔
    [6e/e3e56d] process > GENERATE_REPORT (sample_001) [100%] 3 of 3 ✔
    ```

### Takeaway

In questa sezione, abbiamo imparato le differenze tra il **flusso di dati** (operazioni sui canali) e lo **scripting** (codice all'interno di closure e blocchi di script). Abbiamo usato:

- **Closure** per estrarre e trasformare i metadati dai dati di input
- **Operatori ternari** (`condizione ? valore_se_vera : valore_se_falsa`) per logica condizionale compatta
- **Espressioni regolari** con l'operatore `=~` per estrarre componenti dai nomi dei file
- **Manipolazione delle collezioni** come `join()` per creare stringhe di comandi

Queste tecniche sono fondamentali per scrivere flussi di lavoro Nextflow che siano puliti, efficaci e facili da mantenere.

---

## 2. Elaborazione delle stringhe per gestire nomi di file e metadati

L'elaborazione delle stringhe è un compito comune nei flussi di lavoro bioinformatici, specialmente quando si estraggono metadati dai nomi dei file o si generano script dinamicamente. Nextflow ha potenti capacità di manipolazione delle stringhe che sono cruciali per flussi di lavoro robusti.

### 2.1. Espressioni regolari per analizzare i nomi dei file

Abbiamo già usato espressioni regolari per estrarre metadati dai nomi dei file FASTQ. Comprendiamo meglio come funziona l'operatore di pattern matching `=~`.

Quando scrivi `x =~ /pattern/`:

1. Crea un oggetto `java.util.regex.Matcher`
2. Se valutato in un contesto booleano, verifica se c'è **qualche corrispondenza**
3. Se assegnato a una variabile, puoi accedere ai **gruppi di cattura**

Vediamo come potremmo complicare la nostra regex fastq per gestire più variazioni:

=== "After"

    ```groovy title="main.nf" linenums="14" hl_lines="2"
            // Extract metadata from filename using regex
            def m = (fastq_path.name =~ /^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$/)
            def file_meta = m ? [
                sample_num: m[0][2].toInteger(),
                lane: m[0][3],
                read: m[0][4],
                chunk: m[0][5]
            ] : [:]
    ```

=== "Before"

    ```groovy title="main.nf" linenums="14" hl_lines="2"
            // Extract metadata from filename using regex
            def m = (fastq_path.name =~ /^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq$/)
            def file_meta = m ? [
                sample_num: m[0][2].toInteger(),
                lane: m[0][3],
                read: m[0][4],
                chunk: m[0][5]
            ] : [:]
    ```

Ora, `/^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$/` corrisponderà sia a `.fastq` che a `.fastq.gz`.

!!! note "Gruppi di cattura nelle regex" - Le parentesi `(...)` definiscono "gruppi di cattura" - `m[0]` è la corrispondenza completa - `m[0][1]`, `m[0][2]`, ecc. sono i gruppi di cattura (a partire da 1)

Se avessi nomi di file con convenzioni diverse, potresti usare l'operatore OR (`|`) nella tua regex:

```groovy title="esempio di regex"
def m = (filename =~ /^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$/ |
                      /^(.+)_(\d{6})_([ACGT]+)_L(\d{3})_(R[12])\.fastq(?:\.gz)?$/)
```

### 2.2. Generazione dinamica di script

Un'altra potente applicazione dell'elaborazione delle stringhe è la generazione dinamica di script bash basati su metadati o input. Questo è particolarmente utile per la logica condizionale all'interno dei processi.

Modifichiamo il processo `GENERATE_REPORT` per generare diversi tipi di report basati sulla priorità:

=== "After"

    ```groovy title="modules/generate_report.nf" linenums="10" hl_lines="3-11"
    script:
    // Generate different report content based on priority
    def report_content = meta.priority == 'high' ?
    """
    echo "CAMPIONE AD ALTA PRIORITÀ" > ${meta.id}_report.txt
    echo "===============================================" >> ${meta.id}_report.txt
    """ :
    """
    echo "Campione standard" > ${meta.id}_report.txt
    echo "---------------------------------------------" >> ${meta.id}_report.txt
    """

    """
    ${report_content}
    echo "Sample ID: ${meta.id}" >> ${meta.id}_report.txt
    echo "Organism: ${meta.organism}" >> ${meta.id}_report.txt
    echo "Tissue: ${meta.tissue}" >> ${meta.id}_report.txt
    echo "Sequencing depth: ${meta.depth}" >> ${meta.id}_report.txt
    echo "Quality score: ${meta.quality}" >> ${meta.id}_report.txt
    echo "Sample number: ${meta.sample_num}" >> ${meta.id}_report.txt
    echo "Lane: ${meta.lane}" >> ${meta.id}_report.txt
    echo "Read: ${meta.read}" >> ${meta.id}_report.txt
    echo "Chunk: ${meta.chunk}" >> ${meta.id}_report.txt
    echo "Priority: ${meta.priority}" >> ${meta.id}_report.txt
    echo "File processed: ${reads}" >> ${meta.id}_report.txt
    echo "Process command: fastp --in1 ${reads} --out1 ${meta.id}_trimmed.fastq.gz" >> ${meta.id}_report.txt
    """
    ```

=== "Before"

    ```groovy title="modules/generate_report.nf" linenums="11" hl_lines="2-13"
    script:
    """
    echo "Sample ID: ${meta.id}" > ${meta.id}_report.txt
    echo "Organism: ${meta.organism}" >> ${meta.id}_report.txt
    echo "Tissue: ${meta.tissue}" >> ${meta.id}_report.txt
    echo "Sequencing depth: ${meta.depth}" >> ${meta.id}_report.txt
    echo "Quality score: ${meta.quality}" >> ${meta.id}_report.txt
    echo "Sample number: ${meta.sample_num}" >> ${meta.id}_report.txt
    echo "Lane: ${meta.lane}" >> ${meta.id}_report.txt
    echo "Read: ${meta.read}" >> ${meta.id}_report.txt
    echo "Chunk: ${meta.chunk}" >> ${meta.id}_report.txt
    echo "Priority: ${meta.priority}" >> ${meta.id}_report.txt
    echo "File processed: ${reads}" >> ${meta.id}_report.txt
    echo "Process command: fastp --in1 ${reads} --out1 ${meta.id}_trimmed.fastq.gz" >> ${meta.id}_report.txt
    """
    ```

Eseguiamo il flusso di lavoro per vedere il risultato:

```bash
nextflow run main.nf
```

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [condescending_northcutt] DSL2 - revision: 1a3c16a96f

    executor >  local (6)
    [95/e633a2] process > FASTP (sample_001)           [100%] 3 of 3 ✔
    [a8/e4c214] process > GENERATE_REPORT (sample_001) [100%] 3 of 3 ✔
    ```

Esaminiamo i file di report generati - cerchiamo specificamente quello ad alta priorità:

```console
find work -name "sample_003_report.txt" -exec cat {} \;
```

Dovresti vedere l'intestazione speciale:

```
CAMPIONE AD ALTA PRIORITÀ
===============================================
Sample ID: sample_003
Organism: human
...
```

### 2.3. Interpolazione delle variabili: quando Nextflow valuta vs quando bash valuta

Un punto sottile ma cruciale da comprendere è quando avviene l'interpolazione delle variabili:

1. `${var}` - Interpolato da Nextflow durante la compilazione dello script
2. `\${var}` - Escaped, passato letteralmente a bash come `${var}` (per le variabili d'ambiente bash)

Possiamo vedere questo in azione aggiornando il processo `GENERATE_REPORT` per utilizzare una variabile d'ambiente:

=== "After"

    ```groovy title="modules/generate_report.nf" linenums="11" hl_lines="15"
    script:
    // Generate different report content based on priority
    def report_content = meta.priority == 'high' ?
    """
    echo "CAMPIONE AD ALTA PRIORITÀ" > ${meta.id}_report.txt
    echo "===============================================" >> ${meta.id}_report.txt
    """ :
    """
    echo "Campione standard" > ${meta.id}_report.txt
    echo "---------------------------------------------" >> ${meta.id}_report.txt
    """

    """
    ${report_content}
    echo "Sample ID: ${meta.id}" >> ${meta.id}_report.txt
    echo "Processed by: \${USER}" >> ${meta.id}_report.txt
    echo "Organism: ${meta.organism}" >> ${meta.id}_report.txt
    echo "Tissue: ${meta.tissue}" >> ${meta.id}_report.txt
    echo "Sequencing depth: ${meta.depth}" >> ${meta.id}_report.txt
    echo "Quality score: ${meta.quality}" >> ${meta.id}_report.txt
    echo "Sample number: ${meta.sample_num}" >> ${meta.id}_report.txt
    echo "Lane: ${meta.lane}" >> ${meta.id}_report.txt
    echo "Read: ${meta.read}" >> ${meta.id}_report.txt
    echo "Chunk: ${meta.chunk}" >> ${meta.id}_report.txt
    echo "Priority: ${meta.priority}" >> ${meta.id}_report.txt
    echo "File processed: ${reads}" >> ${meta.id}_report.txt
    echo "Process command: fastp --in1 ${reads} --out1 ${meta.id}_trimmed.fastq.gz" >> ${meta.id}_report.txt
    """
    ```

=== "Before"

    ```groovy title="modules/generate_report.nf" linenums="11"
    script:
    // Generate different report content based on priority
    def report_content = meta.priority == 'high' ?
    """
    echo "CAMPIONE AD ALTA PRIORITÀ" > ${meta.id}_report.txt
    echo "===============================================" >> ${meta.id}_report.txt
    """ :
    """
    echo "Campione standard" > ${meta.id}_report.txt
    echo "---------------------------------------------" >> ${meta.id}_report.txt
    """

    """
    ${report_content}
    echo "Sample ID: ${meta.id}" >> ${meta.id}_report.txt
    echo "Organism: ${meta.organism}" >> ${meta.id}_report.txt
    echo "Tissue: ${meta.tissue}" >> ${meta.id}_report.txt
    echo "Sequencing depth: ${meta.depth}" >> ${meta.id}_report.txt
    echo "Quality score: ${meta.quality}" >> ${meta.id}_report.txt
    echo "Sample number: ${meta.sample_num}" >> ${meta.id}_report.txt
    echo "Lane: ${meta.lane}" >> ${meta.id}_report.txt
    echo "Read: ${meta.read}" >> ${meta.id}_report.txt
    echo "Chunk: ${meta.chunk}" >> ${meta.id}_report.txt
    echo "Priority: ${meta.priority}" >> ${meta.id}_report.txt
    echo "File processed: ${reads}" >> ${meta.id}_report.txt
    echo "Process command: fastp --in1 ${reads} --out1 ${meta.id}_trimmed.fastq.gz" >> ${meta.id}_report.txt
    """
    ```

Se eseguiamo il flusso di lavoro così com'è, fallirà:

```bash
nextflow run main.nf
```

??? failure "Output del comando"

    ```console
    Error modules/generate_report.nf:15:27: `USER` is not defined
    │  15 |     echo "Processed by: ${USER}" >> ${meta.id}_report.txt
    ╰     |                           ^^^^

    ERROR ~ Script compilation failed
    ```

Il problema è che Nextflow cerca di interpretare `${USER}` come una variabile Nextflow, ma non c'è nessuna variabile chiamata `USER`. Dobbiamo fare l'escape del `$` in modo che venga passato a bash, che ha una variabile d'ambiente `USER`:

```groovy title="modules/generate_report.nf" linenums="15" hl_lines="1"
    echo "Processed by: \${USER}" >> ${meta.id}_report.txt
```

Ora funziona! Il backslash (`\`) dice a Nextflow "non interpretare questo, passalo a Bash."

### Takeaway

In questa sezione, ha imparato tecniche di **elaborazione delle stringhe**:

- **Espressioni regolari per l'analisi dei file**: Utilizzare l'operatore `=~` e pattern regex (`~/pattern/`) per estrarre metadati da convenzioni di denominazione dei file complesse
- **Generazione dinamica di script**: Utilizzare logica condizionale (if/else, operatori ternari) per generare diverse stringhe di script basate sulle caratteristiche di input
- **Interpolazione di variabili**: Comprendere quando Nextflow interpreta le stringhe vs quando lo fa la shell
  - `${var}` - Variabili Nextflow (interpolate da Nextflow al momento della compilazione del workflow)
  - `\${var}` - Variabili d'ambiente della shell (escaped, passate a bash a runtime)
  - `\$(cmd)` - Sostituzione di comandi shell (escaped, eseguita da bash a runtime)

Questi pattern di elaborazione e generazione di stringhe sono essenziali per gestire i diversi formati di file e convenzioni di denominazione che incontrerai nei workflow bioinformatici del mondo reale.

---

## 3. Creazione di funzioni riutilizzabili

La logica di workflow complessa inserita negli operatori di canale o nelle definizioni di processo riduce la leggibilità e la manutenibilità. Le **funzioni** ti permettono di estrarre questa logica in componenti nominati e riutilizzabili.

La nostra operazione map è diventata lunga e complessa. Estraiamola in una funzione riutilizzabile usando la keyword `def`.

Per illustrare come appare con il nostro workflow esistente, fai la modifica qui sotto, usando `def` per definire una funzione riutilizzabile chiamata `separateMetadata`:

=== "After"

    ```groovy title="main.nf" linenums="1" hl_lines="4-24 29"
    include { FASTP } from './modules/fastp.nf'
    include { GENERATE_REPORT } from './modules/generate_report.nf'

    def separateMetadata(row) {
        def sample_meta = [
            id: row.sample_id.toLowerCase(),
            organism: row.organism,
            tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
            depth: row.sequencing_depth.toInteger(),
            quality: row.quality_score.toDouble()
        ]
        def fastq_path = file(row.file_path)

        def m = (fastq_path.name =~ /^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$/)
        def file_meta = m ? [
            sample_num: m[0][2].toInteger(),
            lane: m[0][3],
            read: m[0][4],
            chunk: m[0][5]
        ] : [:]

        def priority = sample_meta.quality > 40 ? 'high' : 'normal'
        return tuple(sample_meta + file_meta + [priority: priority], fastq_path)
    }

    workflow {
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map{ row -> separateMetadata(row) }

        ch_fastp = FASTP(ch_samples)
        GENERATE_REPORT(ch_samples)
    }
    ```

=== "Before"

    ```groovy title="main.nf" linenums="1" hl_lines="7-27"
    include { FASTP } from './modules/fastp.nf'
    include { GENERATE_REPORT } from './modules/generate_report.nf'

    workflow {
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                def fastq_path = file(row.file_path)

                def m = (fastq_path.name =~ /^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$/)
                def file_meta = m ? [
                    sample_num: m[0][2].toInteger(),
                    lane: m[0][3],
                    read: m[0][4],
                    chunk: m[0][5]
                ] : [:]

                def priority = sample_meta.quality > 40 ? 'high' : 'normal'
                return tuple(sample_meta + file_meta + [priority: priority], fastq_path)
            }

        ch_fastp = FASTP(ch_samples)
        GENERATE_REPORT(ch_samples)
    }
    ```

Estraendo questa logica in una funzione, abbiamo ridotto la logica effettiva del workflow a qualcosa di molto più pulito:

```groovy title="workflow minimo"
    ch_samples = channel.fromPath("./data/samples.csv")
        .splitCsv(header: true)
        .map{ row -> separateMetadata(row) }

    ch_fastp = FASTP(ch_samples)
    GENERATE_REPORT(ch_samples)
```

Questo rende la logica del workflow molto più facile da leggere e capire a colpo d'occhio. La funzione `separateMetadata` incapsula tutta la logica complessa per analizzare e arricchire i metadati, rendendola riutilizzabile e testabile.

Esegui il workflow per assicurarti che funzioni ancora:

```bash
nextflow run main.nf
```

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [admiring_panini] DSL2 - revision: 8cc832e32f

    executor >  local (6)
    [8c/2e3f91] process > FASTP (3)           [100%] 3 of 3 ✔
    [7a/1b4c92] process > GENERATE_REPORT (3) [100%] 3 of 3 ✔
    ```

L'output dovrebbe mostrare che entrambi i processi si completano con successo. Il workflow è ora molto più pulito e facile da mantenere, con tutta la logica complessa di elaborazione dei metadati incapsulata nella funzione `separateMetadata`.

### Takeaway

In questa sezione, hai imparato la **creazione di funzioni**:

- **Definizione di funzioni con `def`**: La keyword per creare funzioni con nome (come `def` in Python o `function` in JavaScript)
- **Scope delle funzioni**: Le funzioni definite a livello di script sono accessibili in tutto il tuo workflow Nextflow
- **Valori di ritorno**: Le funzioni restituiscono automaticamente l'ultima espressione, o usano un esplicito `return`
- **Codice più pulito**: Estrarre logica complessa in funzioni è una pratica fondamentale di ingegneria del software in qualsiasi linguaggio

Successivamente, esploreremo come usare le closure nelle direttive di processo per l'allocazione dinamica delle risorse.

---

## 4. Direttive di risorse dinamiche con Closure

Finora abbiamo usato lo scripting nel blocco `script` dei processi. Ma le **closure** (introdotte nella Sezione 1.1) sono anche incredibilmente utili nelle direttive dei processi, specialmente per l'allocazione dinamica delle risorse. Aggiungiamo direttive di risorse al nostro processo FASTP che si adattano in base alle caratteristiche del campione.

### 4.1. Allocazione di risorse specifiche per campione

Attualmente, il nostro processo FASTP utilizza risorse predefinite. Rendiamolo più intelligente allocando più CPU per i campioni ad alta profondità. Modifica `modules/fastp.nf` per includere una direttiva `cpus` dinamica e una direttiva `memory` statica:

=== "After"

    ```groovy title="modules/fastp.nf" linenums="1" hl_lines="4-5"
    process FASTP {
        container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

        cpus { meta.depth > 40000000 ? 2 : 1 }
        memory 2.GB

        input:
        tuple val(meta), path(reads)
    ```

=== "Before"

    ```groovy title="modules/fastp.nf" linenums="1"
    process FASTP {
        container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

        input:
        tuple val(meta), path(reads)
    ```

La closure `{ meta.depth > 40000000 ? 2 : 1 }` usa l'**operatore ternario** (trattato nella Sezione 1.1) e viene valutata per ogni task, permettendo l'allocazione di risorse per campione. I campioni ad alta profondità (>40M reads) ottengono 2 CPU, mentre gli altri ottengono 1 CPU.

!!! note "Accesso alle variabili di input nelle direttive"

    La closure può accedere a qualsiasi variabile di input (come `meta` qui) perché Nextflow valuta queste closure nel contesto di ogni esecuzione di task.

Esegui di nuovo il workflow con l'opzione `-ansi-log false` per rendere più facile vedere gli hash dei task.

```bash
nextflow run main.nf -ansi-log false
```

??? success "Output del comando"

    ```console
    N E X T F L O W  ~  version 25.10.2
    Launching `main.nf` [fervent_albattani] DSL2 - revision: fa8f249759
    [bd/ff3d41] Submitted process > FASTP (2)
    [a4/a3aab2] Submitted process > FASTP (1)
    [48/6db0c9] Submitted process > FASTP (3)
    [ec/83439d] Submitted process > GENERATE_REPORT (3)
    [bd/15d7cc] Submitted process > GENERATE_REPORT (2)
    [42/699357] Submitted process > GENERATE_REPORT (1)
    ```

Puoi controllare l'esatto comando `docker` che è stato eseguito per vedere l'allocazione CPU per qualsiasi task dato:

```console title="Controlla comando docker"
cat work/48/6db0c9e9d8aa65e4bb4936cd3bd59e/.command.run | grep "docker run"
```

Dovresti vedere qualcosa del genere:

```bash title="comando docker"
    docker run -i --cpu-shares 4096 --memory 2048m -e "NXF_TASK_WORKDIR" -v /workspaces/training/side-quests/essential_scripting_patterns:/workspaces/training/side-quests/essential_scripting_patterns -w "$NXF_TASK_WORKDIR" --name $NXF_BOXID community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690 /bin/bash -ue /workspaces/training/side-quests/essential_scripting_patterns/work/48/6db0c9e9d8aa65e4bb4936cd3bd59e/.command.sh
```

In questo esempio abbiamo scelto un esempio che ha richiesto 2 CPU (`--cpu-shares 2048`), perché era un campione ad alta profondità, ma dovresti vedere diverse allocazioni di CPU a seconda della profondità del campione. Prova questo anche per gli altri task.

### 4.2. Strategie di ritentativo

Un altro pattern potente è l'uso di `task.attempt` per strategie di ritentativo. Per mostrare perché questo è utile, inizieremo riducendo l'allocazione di memoria per FASTP a meno di quanto ne abbia bisogno. Cambia la direttiva `memory` in `modules/fastp.nf` a `1.GB`:

=== "After"

    ```groovy title="modules/fastp.nf" linenums="1" hl_lines="5"
    process FASTP {
        container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

        cpus { meta.depth > 40000000 ? 4 : 2 }
        memory 1.GB

        input:
        tuple val(meta), path(reads)
    ```

=== "Before"

    ```groovy title="modules/fastp.nf" linenums="1" hl_lines="5"
    process FASTP {
        container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

        cpus { meta.depth > 40000000 ? 4 : 2 }
        memory 2.GB

        input:
        tuple val(meta), path(reads)
    ```

... ed esegui di nuovo il workflow:

```bash
nextflow run main.nf
```

??? failure "Output del comando"

    ```console hl_lines="2 11"
    Command exit status:
      137

    Command output:
      (empty)

    Command error:
      Detecting adapter sequence for read1...
      No adapter detected for read1

      .command.sh: line 7:   101 Killed                  fastp --in1 SAMPLE_002_S2_L001_R1_001.fastq --out1 sample_002_trimmed.fastq.gz --json sample_002.fastp.json --html sample_002.fastp.html --thread 2
    ```

Questo indica che il processo è stato terminato per aver superato i limiti di memoria.

Questo è uno scenario molto comune nei workflow del mondo reale - a volte semplicemente non sai quanta memoria avrà bisogno un task finché non lo esegui.

Per rendere il nostro workflow più robusto, possiamo implementare una strategia di ritentativo che aumenta l'allocazione di memoria ad ogni tentativo, ancora una volta usando una closure Groovy. Modifica la direttiva `memory` per moltiplicare la memoria base per `task.attempt`, e aggiungi le direttive `errorStrategy 'retry'` e `maxRetries 2`:

=== "After"

    ```groovy title="modules/fastp.nf" linenums="1" hl_lines="5-7"
    process FASTP {
        container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

        cpus { meta.depth > 40000000 ? 4 : 2 }
        memory { 1.GB * task.attempt }
        errorStrategy 'retry'
        maxRetries 2

        input:
        tuple val(meta), path(reads)
    ```

=== "Before"

    ```groovy title="modules/fastp.nf" linenums="1" hl_lines="5"
    process FASTP {
        container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

        cpus { meta.depth > 40000000 ? 4 : 2 }
        memory 2.GB

        input:
        tuple val(meta), path(reads)
    ```

Ora se il processo fallisce a causa di memoria insufficiente, Nextflow ritenterà con più memoria:

- Primo tentativo: 1 GB (task.attempt = 1)
- Secondo tentativo: 2.GB (task.attempt = 2)

... e così via, fino al limite `maxRetries`.

### Takeaway

Le direttive dinamiche con closure ti permettono di:

- Allocare risorse in base alle caratteristiche di input
- Implementare strategie di ritentativo automatiche con risorse crescenti
- Combinare molteplici fattori (metadati, numero di tentativi, priorità)
- Utilizzare logica condizionale per calcoli di risorse complessi

Questo rende i tuoi workflow sia più efficienti (senza sovra-allocazione) che più robusti (ritentativo automatico con più risorse).

---

## 5. Logica condizionale e controllo del processo

In precedenza, abbiamo usato `.map()` con lo scripting per trasformare i dati del canale. Ora useremo la logica condizionale per controllare quali processi eseguire in base ai dati — essenziale per flussi di lavoro flessibili che si adattano a diversi tipi di campione.

Gli [operatori di flusso dati](https://www.nextflow.io/docs/latest/reference/operator.html) di Nextflow prendono closure valutate a runtime, abilitando la logica condizionale per guidare le decisioni del workflow basate sul contenuto del canale.

### 5.1. Routing con `.branch()`

Ad esempio, immaginiamo che i nostri campioni di sequenziamento debbano essere trimmati con FASTP solo se sono campioni umani con una copertura sopra una certa soglia. I campioni di topo o campioni a bassa copertura dovrebbero invece essere eseguiti con Trimgalore (questo è un esempio artificiale, ma illustra il punto).

Abbiamo fornito un semplice processo Trimgalore in `modules/trimgalore.nf`, dai un'occhiata se vuoi, ma i dettagli non sono importanti per questo esercizio. Il punto chiave è che vogliamo instradare i campioni in base ai loro metadati.

Includi il nuovo modulo da `modules/trimgalore.nf`:

=== "After"

    ```groovy title="main.nf" linenums="1" hl_lines="2"
    include { FASTP } from './modules/fastp.nf'
    include { TRIMGALORE } from './modules/trimgalore.nf'
    ```

=== "Before"

    ```groovy title="main.nf" linenums="1"
    include { FASTP } from './modules/fastp.nf'
    ```

... e poi modifica il tuo workflow `main.nf` per suddividere i campioni in base ai loro metadati e instradarli attraverso il processo di trimming appropriato, così:

=== "After"

    ```groovy title="main.nf" linenums="28" hl_lines="5-12"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row -> separateMetadata(row) }

        trim_branches = ch_samples
            .branch { meta, reads ->
                fastp: meta.organism == 'human' && meta.depth >= 30000000
                trimgalore: true
            }

        ch_fastp = FASTP(trim_branches.fastp)
        ch_trimgalore = TRIMGALORE(trim_branches.trimgalore)
        GENERATE_REPORT(ch_samples)
    ```

=== "Before"

    ```groovy title="main.nf" linenums="28" hl_lines="5"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row -> separateMetadata(row) }

        ch_fastp = FASTP(ch_samples)
        GENERATE_REPORT(ch_samples)
    ```

Esegui questo workflow modificato:

```bash
nextflow run main.nf
```

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [adoring_galileo] DSL2 - revision: c9e83aaef1

    executor >  local (6)
    [1d/0747ac] process > FASTP (2)           [100%] 2 of 2 ✔
    [cc/c44caf] process > TRIMGALORE (1)      [100%] 1 of 1 ✔
    [34/bd5a9f] process > GENERATE_REPORT (1) [100%] 3 of 3 ✔
    ```

Qui, abbiamo usato piccole ma potenti espressioni condizionali all'interno dell'operatore `.branch{}` per instradare i campioni in base ai loro metadati. I campioni umani con alta copertura passano attraverso `FASTP`, mentre tutti gli altri campioni passano attraverso `TRIMGALORE`.

### 5.2. Uso di `.filter()` con Truthiness

Un altro pattern potente per controllare l'esecuzione del workflow è l'operatore `.filter()`, che utilizza una closure per determinare quali elementi dovrebbero continuare nel pipeline. All'interno della closure di filtraggio, scriverai **espressioni booleane** che decidono quali elementi passano.

Nextflow (come molti linguaggi dinamici) ha un concetto di **"truthiness"** che determina quali valori vengono valutati come `true` o `false` in contesti booleani:

- **Truthy**: Valori non-null, stringhe non vuote, numeri non-zero, collezioni non vuote
- **Falsy**: `null`, stringhe vuote `""`, zero `0`, collezioni vuote `[]` o `[:]`, `false`

Questo significa che `meta.id` da solo (senza un esplicito `!= null`) controlla se l'ID esiste e non è vuoto. Usiamo questo per filtrare i campioni che non soddisfano i nostri requisiti di qualità.

Aggiungi quanto segue prima dell'operazione branch:

=== "After"

    ```groovy title="main.nf" linenums="28" hl_lines="5-11"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row -> separateMetadata(row) }

        // Filtra campioni invalidi o di bassa qualità
        ch_valid_samples = ch_samples
            .filter { meta, reads ->
                meta.id && meta.organism && meta.depth >= 25000000
            }

        trim_branches = ch_valid_samples
            .branch { meta, reads ->
                fastp: meta.organism == 'human' && meta.depth >= 30000000
                trimgalore: true
            }
    ```

=== "Before"

    ```groovy title="main.nf" linenums="28" hl_lines="5"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row -> separateMetadata(row) }

        trim_branches = ch_samples
            .branch { meta, reads ->
                fastp: meta.organism == 'human' && meta.depth >= 30000000
                trimgalore: true
            }
    ```

Esegui di nuovo il workflow:

```bash
nextflow run main.nf
```

??? success "Output del comando"

    ```console
    N E X T F L O W  ~  version 25.10.2
    Launching `main.nf` [lonely_williams] DSL2 - revision: d0b3f121ec
    [94/b48eac] Submitted process > FASTP (2)
    [2c/d2b28f] Submitted process > GENERATE_REPORT (2)
    [65/2e3be4] Submitted process > GENERATE_REPORT (1)
    [94/b48eac] NOTE: Process `FASTP (2)` terminated with an error exit status (137) -- Execution is retried (1)
    [3e/0d8664] Submitted process > TRIMGALORE (1)
    [6a/9137b0] Submitted process > FASTP (1)
    [6a/9137b0] NOTE: Process `FASTP (1)` terminated with an error exit status (137) -- Execution is retried (1)
    [83/577ac0] Submitted process > GENERATE_REPORT (3)
    [a2/5117de] Re-submitted process > FASTP (1)
    [1f/a1a4ca] Re-submitted process > FASTP (2)
    ```

Poiché abbiamo scelto un filtro che esclude alcuni campioni, sono stati eseguiti meno task.

L'espressione di filtro `meta.id && meta.organism && meta.depth >= 25000000` combina truthiness con confronti espliciti:

- `meta.id && meta.organism` verifica che entrambi i campi esistano e non siano vuoti (usando truthiness)
- `meta.depth >= 25000000` garantisce una profondità di sequenziamento sufficiente con un confronto esplicito

!!! note "Truthiness in pratica"

    L'espressione `meta.id && meta.organism` è più concisa che scrivere:
    ```groovy
    meta.id != null && meta.id != '' && meta.organism != null && meta.organism != ''
    ```

    Questo rende la logica di filtraggio molto più pulita e facile da leggere.

### Takeaway

In questa sezione, hai imparato a usare la logica condizionale per controllare l'esecuzione del workflow usando le interfacce closure degli operatori Nextflow come `.branch{}` e `.filter{}`, sfruttando truthiness per scrivere espressioni condizionali concise.

Il nostro pipeline ora instrada intelligentemente i campioni attraverso processi appropriati, ma i workflow di produzione devono gestire elegantemente i dati invalidi. Rendiamo il nostro workflow robusto contro valori mancanti o nulli.

---

## 6. Operatori di navigazione sicura ed Elvis

La nostra funzione `separateMetadata` attualmente presuppone che tutti i campi CSV siano presenti e validi. Ma cosa succede con dati incompleti? Scopriamolo.

### 6.1. Il Problema: Accedere a proprietà che non esistono

Diciamo che vogliamo aggiungere supporto per informazioni opzionali sulla sequenza di esecuzione. In alcuni laboratori, i campioni potrebbero avere un campo aggiuntivo per l'ID di sequenza o il numero di lotto, ma il nostro CSV attuale non ha questa colonna. Proviamo ad accedervi comunque.

Modifica la funzione `separateMetadata` per includere un campo run_id:

=== "After"

    ```groovy title="main.nf" linenums="5" hl_lines="9"
    def separateMetadata(row) {
        def sample_meta = [
            id: row.sample_id.toLowerCase(),
            organism: row.organism,
            tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
            depth: row.sequencing_depth.toInteger(),
            quality: row.quality_score.toDouble()
        ]
        def run_id = row.run_id.toUpperCase()
    ```

=== "Before"

    ```groovy title="main.nf" linenums="5"
    def separateMetadata(row) {
        def sample_meta = [
            id: row.sample_id.toLowerCase(),
            organism: row.organism,
            tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
            depth: row.sequencing_depth.toInteger(),
            quality: row.quality_score.toDouble()
        ]
    ```

Ora esegui il workflow:

```bash
nextflow run main.nf
```

??? failure "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [trusting_torvalds] DSL2 - revision: b56fbfbce2

    ERROR ~ Cannot invoke method toUpperCase() on null object

    -- Check script 'main.nf' at line: 13 or see '.nextflow.log' file for more details
    ```

Questo fallisce con una NullPointerException.

Il problema è che `row.run_id` restituisce `null` perché la colonna `run_id` non esiste nel nostro CSV. Quando proviamo a chiamare `.toUpperCase()` su `null`, si blocca. È qui che l'operatore di navigazione sicura salva la situazione.

### 6.2. Operatore di navigazione sicura (`?.`)

L'operatore di navigazione sicura (`?.`) restituisce `null` invece di lanciare un'eccezione quando viene chiamato su un valore `null`. Se l'oggetto prima di `?.` è `null`, l'intera espressione viene valutata come `null` senza eseguire il metodo.

Aggiorna la funzione per usare la navigazione sicura:

=== "After"

    ```groovy title="main.nf" linenums="4" hl_lines="9"
    def separateMetadata(row) {
        def sample_meta = [
            id: row.sample_id.toLowerCase(),
            organism: row.organism,
            tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
            depth: row.sequencing_depth.toInteger(),
            quality: row.quality_score.toDouble()
        ]
        def run_id = row.run_id?.toUpperCase()
    ```

=== "Before"

    ```groovy title="main.nf" linenums="4" hl_lines="9"
    def separateMetadata(row) {
        def sample_meta = [
            id: row.sample_id.toLowerCase(),
            organism: row.organism,
            tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
            depth: row.sequencing_depth.toInteger(),
            quality: row.quality_score.toDouble()
        ]
        def run_id = row.run_id.toUpperCase()
    ```

Esegui di nuovo:

```bash
nextflow run main.nf
```

??? success "Output del comando"

    ```console
    <!-- TODO: output -->
    ```

Nessun crash! Il workflow ora gestisce con eleganza il campo mancante. Quando `row.run_id` è `null`, l'operatore `?.` impedisce la chiamata `.toUpperCase()`, e `run_id` diventa `null` invece di causare un'eccezione.

### 6.3. Operatore Elvis (`?:`) per valori predefiniti

L'operatore Elvis (`?:`) fornisce valori predefiniti quando il lato sinistro è "falsy" (come spiegato in precedenza). È chiamato così per Elvis Presley perché `?:` assomiglia ai suoi famosi capelli e occhi quando visti di lato!

Ora che stiamo usando la navigazione sicura, `run_id` sarà `null` per i campioni senza quel campo. Usiamo l'operatore Elvis per fornire un valore predefinito e aggiungerlo alla nostra mappa `sample_meta`:

=== "After"

    ```groovy title="main.nf" linenums="5" hl_lines="9-10"
    def separateMetadata(row) {
        def sample_meta = [
            id: row.sample_id.toLowerCase(),
            organism: row.organism,
            tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
            depth: row.sequencing_depth.toInteger(),
            quality: row.quality_score.toDouble()
        ]
        def run_id = row.run_id?.toUpperCase() ?: 'NON_SPECIFICATO'
        sample_meta.run = run_id
    ```

=== "Before"

    ```groovy title="main.nf" linenums="5" hl_lines="9"
    def separateMetadata(row) {
        def sample_meta = [
            id: row.sample_id.toLowerCase(),
            organism: row.organism,
            tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
            depth: row.sequencing_depth.toInteger(),
            quality: row.quality_score.toDouble()
        ]
        def run_id = row.run_id?.toUpperCase()
    ```

Aggiungi anche un operatore `view()` nel workflow per vedere i risultati:

=== "After"

    ```groovy title="main.nf" linenums="30" hl_lines="4"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map{ row -> separateMetadata(row) }
            .view()
    ```

=== "Before"

    ```groovy title="main.nf" linenums="30"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map{ row -> separateMetadata(row) }
    ```

ed esegui il workflow:

```bash
nextflow run main.nf
```

??? success "Output del comando"

    ```console
    [[id:sample_001, organism:human, tissue:liver, depth:30000000, quality:38.5, run:NON_SPECIFICATO, sample_num:1, lane:001, read:R1, chunk:001, priority:normal], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_001_S1_L001_R1_001.fastq]
    [[id:sample_002, organism:mouse, tissue:brain, depth:25000000, quality:35.2, run:NON_SPECIFICATO, sample_num:2, lane:001, read:R1, chunk:001, priority:normal], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_002_S2_L001_R1_001.fastq]
    [[id:sample_003, organism:human, tissue:kidney, depth:45000000, quality:42.1, run:NON_SPECIFICATO, sample_num:3, lane:001, read:R1, chunk:001, priority:high], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_003_S3_L001_R1_001.fastq]
    ```

Perfetto! Ora tutti i campioni hanno un campo `run` con il loro ID di esecuzione effettivo (in maiuscolo) o il valore predefinito 'NON_SPECIFICATO'. La combinazione di `?.` e `?:` fornisce sia sicurezza (nessun crash) che valori predefiniti sensati.

Rimuovi l'operatore `.view()` ora che abbiamo confermato che funziona.

!!! tip "Combinare navigazione sicura ed Elvis"

    Il pattern `value?.method() ?: 'default'` è comune nei workflow di produzione:

    - `value?.method()` - Chiama il metodo in sicurezza, restituisce `null` se `value` è `null`
    - `?: 'default'` - Fornisce un fallback se il risultato è `null`

    Questo pattern gestisce con eleganza i dati mancanti/incompleti.

Usa questi operatori in modo coerente in funzioni, closure di operatori (`.map{}`, `.filter{}`), script di processo e file di configurazione. Prevengono crash quando si gestiscono dati del mondo reale.

### Takeaway

- **Navigazione sicura (`?.`)**: Previene crash su valori nulli - restituisce null invece di lanciare un'eccezione
- **Operatore Elvis (`?:`)**: Fornisce valori predefiniti - `value ?: 'default'`
- **Combinazione**: `value?.method() ?: 'default'` è il pattern comune

Questi operatori rendono i workflow resistenti a dati incompleti - essenziali per il lavoro nel mondo reale.

---

## 7. Validazione con `error()` e `log.warn`

A volte è necessario fermare immediatamente il workflow se i parametri di input non sono validi. In Nextflow, puoi usare funzioni integrate come `error()` e `log.warn`, così come costrutti di programmazione standard come istruzioni `if` e logica booleana, per implementare la logica di validazione. Aggiungiamo la validazione al nostro workflow.

Crea una funzione di validazione prima del tuo blocco workflow, chiamala dal workflow, e cambia la creazione del canale per usare un parametro per il percorso del file CSV. Se il parametro manca o il file non esiste, chiama `error()` per fermare l'esecuzione con un messaggio chiaro.

=== "After"

    ```groovy title="main.nf" linenums="1" hl_lines="5-20 23-24"
    include { FASTP } from './modules/fastp.nf'
    include { TRIMGALORE } from './modules/trimgalore.nf'
    include { GENERATE_REPORT } from './modules/generate_report.nf'

    def validateInputs() {
        // Controlla che il parametro di input sia fornito
        if (!params.input) {
            error("Percorso file CSV di input non fornito. Per favore specifica --input <file.csv>")
        }

        // Controlla che il file CSV esista
        if (!file(params.input).exists()) {
            error("File CSV di input non trovato: ${params.input}")
        }
    }
    ...
    workflow {
        validateInputs()
        ch_samples = channel.fromPath(params.input)
    ```

=== "Before"

    ```groovy title="main.nf" linenums="1"
    include { FASTP } from './modules/fastp.nf'
    include { TRIMGALORE } from './modules/trimgalore.nf'
    include { GENERATE_REPORT } from './modules/generate_report.nf'

    ...
    workflow {
        ch_samples = channel.fromPath("./data/samples.csv")
    ```

Ora prova a eseguire senza il file CSV:

```bash
nextflow run main.nf
```

??? failure "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [confident_coulomb] DSL2 - revision: 07059399ed

    WARN: Access to undefined parameter `input` -- Initialise it to a default value eg. `params.input = some_value`
    Percorso file CSV di input non fornito. Per favore specifica --input <file.csv>
    ```

Il workflow si ferma immediatamente con un messaggio di errore chiaro invece di fallire misteriosamente più tardi

Ora eseguilo con un file inesistente:

```bash
nextflow run main.nf --input ./data/nonexistent.csv
```

??? failure "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [cranky_gates] DSL2 - revision: 26839ae3eb

    File CSV di input non trovato: ./data/nonexistent.csv
    ```

Infine, eseguilo con il file corretto:

```bash
nextflow run main.nf --input ./data/samples.csv
```

??? success "Output del comando"

    ```console
    <!-- TODO: output -->
    ```

Questa volta viene eseguito con successo.

Puoi anche aggiungere la validazione all'interno della funzione `separateMetadata`. Usiamo il non-fatal `log.warn` per emettere avvisi per i campioni con bassa profondità di sequenziamento, ma permettiamo comunque al workflow di continuare:

=== "After"

    ```groovy title="main.nf" linenums="1" hl_lines="3-6"
        def priority = sample_meta.quality > 40 ? 'high' : 'normal'

        // Validazione che i dati abbiano senso
        if (sample_meta.depth < 30000000) {
            log.warn "Bassa profondità di sequenziamento per ${sample_meta.id}: ${sample_meta.depth}"
        }

        return tuple(sample_meta + file_meta + [priority: priority], fastq_path)
    }
    ```

=== "Before"

    ```groovy title="main.nf" linenums="1"
        def priority = sample_meta.quality > 40 ? 'high' : 'normal'

        return tuple(sample_meta + file_meta + [priority: priority], fastq_path)
    }
    ```

Esegui di nuovo il workflow con il CSV originale:

```bash
nextflow run main.nf --input ./data/samples.csv
```

??? warning "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [awesome_goldwasser] DSL2 - revision: a31662a7c1

    executor >  local (5)
    [ce/df5eeb] process > FASTP (2)           [100%] 2 of 2 ✔
    [-        ] process > TRIMGALORE          -
    [d1/7d2b4b] process > GENERATE_REPORT (3) [100%] 3 of 3 ✔
    WARN: Bassa profondità di sequenziamento per sample_002: 25000000
    ```

Vediamo un avviso riguardo alla bassa profondità di sequenziamento per uno dei campioni.

### Takeaway

- **`error()`**: Ferma il workflow immediatamente con un messaggio chiaro
- **`log.warn`**: Emette avvisi senza fermare il workflow
- **Validazione precoce**: Controlla gli input prima dell'elaborazione per fallire velocemente con errori utili
- **Funzioni di validazione**: Crea logica di validazione riutilizzabile che può essere chiamata all'inizio del workflow

Una validazione appropriata rende i workflow più robusti e facili da usare rilevando i problemi presto con messaggi di errore chiari.

---

## 8. Gestori di eventi del workflow

Finora abbiamo scritto codice nei nostri script di workflow e definizioni di processo. Ma c'è un'altra importante caratteristica che dovresti conoscere: i gestori di eventi del workflow.

I gestori di eventi sono closure che vengono eseguite in momenti specifici del ciclo di vita del tuo workflow. Sono perfetti per aggiungere logging, notifiche o operazioni di pulizia. Questi gestori dovrebbero essere definiti nel tuo script di workflow accanto alla definizione del workflow.

### 8.1. Il gestore `onComplete`

Il gestore di eventi più comunemente utilizzato è `onComplete`, che viene eseguito quando il tuo workflow termina (sia che abbia avuto successo o sia fallito). Aggiungiamone uno per riassumere i risultati del nostro pipeline.

Aggiungi il gestore di eventi al tuo file `main.nf`, all'interno della tua definizione di workflow:

=== "After"

    ```groovy title="main.nf" linenums="66" hl_lines="5-16"
        ch_fastp = FASTP(trim_branches.fastp)
        ch_trimgalore = TRIMGALORE(trim_branches.trimgalore)
        GENERATE_REPORT(ch_samples)

        workflow.onComplete = {
            println ""
            println "Riepilogo dell'esecuzione del pipeline:"
            println "=========================="
            println "Completato a: ${workflow.complete}"
            println "Durata      : ${workflow.duration}"
            println "Successo    : ${workflow.success}"
            println "workDir     : ${workflow.workDir}"
            println "stato uscita: ${workflow.exitStatus}"
            println ""
        }
    }
    ```

=== "Before"

    ```groovy title="main.nf" linenums="66" hl_lines="4"
        ch_fastp = FASTP(trim_branches.fastp)
        ch_trimgalore = TRIMGALORE(trim_branches.trimgalore)
        GENERATE_REPORT(ch_samples)
    }
    ```

Questa closure viene eseguita quando il workflow si completa. All'interno, hai accesso all'oggetto `workflow` che fornisce proprietà utili sull'esecuzione.

Esegui il tuo workflow e vedrai questo riepilogo apparire alla fine!

```bash
nextflow run main.nf --input ./data/samples.csv -ansi-log false
```

??? success "Output del comando"

    ```console
    N E X T F L O W  ~  version 25.10.2
    Launching `main.nf` [marvelous_boltzmann] DSL2 - revision: a31662a7c1
    WARN: Bassa profondità di sequenziamento per sample_002: 25000000
    [9b/d48e40] Submitted process > FASTP (2)
    [6a/73867a] Submitted process > GENERATE_REPORT (2)
    [79/ad0ac5] Submitted process > GENERATE_REPORT (1)
    [f3/bda6cb] Submitted process > FASTP (1)
    [34/d5b52f] Submitted process > GENERATE_REPORT (3)

    Riepilogo dell'esecuzione del pipeline:
    ==========================
    Completato a: 2025-10-10T12:14:24.885384+01:00
    Durata      : 2.9s
    Successo    : true
    workDir     : /workspaces/training/side-quests/essential_scripting_patterns/work
    stato uscita: 0
    ```

Rendiamolo più utile aggiungendo una logica condizionale:

=== "After"

    ```groovy title="main.nf" linenums="66" hl_lines="5-22"
        ch_fastp = FASTP(trim_branches.fastp)
        ch_trimgalore = TRIMGALORE(trim_branches.trimgalore)
        GENERATE_REPORT(ch_samples)

        workflow.onComplete = {
            println ""
            println "Riepilogo dell'esecuzione del pipeline:"
            println "=========================="
            println "Completato a: ${workflow.complete}"
            println "Durata      : ${workflow.duration}"
            println "Successo    : ${workflow.success}"
            println "workDir     : ${workflow.workDir}"
            println "stato uscita: ${workflow.exitStatus}"
            println ""

            if (workflow.success) {
                println "✅ Pipeline completato con successo!"
            } else {
                println "❌ Pipeline fallito!"
                println "Errore: ${workflow.errorMessage}"
            }
        }
    }
    ```

=== "Before"

    ```groovy title="main.nf" linenums="66" hl_lines="5-16"
        ch_fastp = FASTP(trim_branches.fastp)
        ch_trimgalore = TRIMGALORE(trim_branches.trimgalore)
        GENERATE_REPORT(ch_samples)

        workflow.onComplete = {
            println ""
            println "Riepilogo dell'esecuzione del pipeline:"
            println "=========================="
            println "Completato a: ${workflow.complete}"
            println "Durata      : ${workflow.duration}"
            println "Successo    : ${workflow.success}"
            println "workDir     : ${workflow.workDir}"
            println "stato uscita: ${workflow.exitStatus}"
            println ""
        }
    }
    ```

Ora otteniamo un riepilogo ancora più informativo, incluso un messaggio di successo/fallimento e la directory di output se specificata:

<!-- TODO: add run command -->

??? success "Output del comando"

    ```console
    N E X T F L O W  ~  version 25.10.2
    Launching `main.nf` [boring_linnaeus] DSL2 - revision: a31662a7c1
    WARN: Bassa profondità di sequenziamento per sample_002: 25000000
    [e5/242efc] Submitted process > FASTP (2)
    [3b/74047c] Submitted process > GENERATE_REPORT (3)
    [8a/7a57e6] Submitted process > GENERATE_REPORT (1)
    [a8/b1a31f] Submitted process > GENERATE_REPORT (2)
    [40/648429] Submitted process > FASTP (1)

    Riepilogo dell'esecuzione del pipeline:
    ==========================
    Completato a: 2025-10-10T12:16:00.522569+01:00
    Durata      : 3.6s
    Successo    : true
    workDir     : /workspaces/training/side-quests/essential_scripting_patterns/work
    stato uscita: 0

    ✅ Pipeline completato con successo!
    ```

Puoi anche scrivere il riepilogo su un file usando operazioni di file:

```groovy title="main.nf - Scrittura del riepilogo su file"
workflow {
    // ... il tuo codice di workflow ...

    workflow.onComplete = {
        def summary = """
        Riepilogo dell'esecuzione del pipeline
        ===========================
        Completato: ${workflow.complete}
        Durata    : ${workflow.duration}
        Successo  : ${workflow.success}
        Comando   : ${workflow.commandLine}
        """

        println summary

        // Scrivi su un file di log
        def log_file = file("${workflow.launchDir}/pipeline_summary.txt")
        log_file.text = summary
    }
}
```

### 8.2. Il gestore `onError`

Oltre a `onComplete`, c'è un altro gestore di eventi che puoi utilizzare: `onError`, che viene eseguito solo se il workflow fallisce:

```groovy title="main.nf - gestore onError"
workflow {
    // ... il tuo codice di workflow ...

    workflow.onError = {
        println "="* 50
        println "L'esecuzione del pipeline è fallita!"
        println "Messaggio di errore: ${workflow.errorMessage}"
        println "="* 50

        // Scrivi log di errore dettagliato
        def error_file = file("${workflow.launchDir}/error.log")
        error_file.text = """
        Rapporto di errore del workflow
        =====================
        Ora: ${new Date()}
        Errore: ${workflow.errorMessage}
        Rapporto di errore: ${workflow.errorReport ?: 'Nessun rapporto dettagliato disponibile'}
        """

        println "Dettagli dell'errore scritti in: ${error_file}"
    }
}
```

Puoi utilizzare più gestori insieme nel tuo script di workflow:

```groovy title="main.nf - Gestori combinati"
workflow {
    // ... il tuo codice di workflow ...

    workflow.onError = {
        println "Workflow fallito: ${workflow.errorMessage}"
    }

    workflow.onComplete = {
        def duration_mins = workflow.duration.toMinutes().round(2)
        def status = workflow.success ? "SUCCESSO ✅" : "FALLITO ❌"

        println """
        Pipeline terminato: ${status}
        Durata: ${duration_mins} minuti
        """
    }
}
```

### Takeaway

In questa sezione, hai imparato:

- **Closure di gestori di eventi**: Closure nel tuo script di workflow che vengono eseguite in diversi punti del ciclo di vita
- **Gestore `onComplete`**: Per riepiloghi di esecuzione e reportistica dei risultati
- **Gestore `onError`**: Per gestione degli errori e logging dei fallimenti
- **Proprietà dell'oggetto workflow**: Accedere a `workflow.success`, `workflow.duration`, `workflow.errorMessage`, ecc.

I gestori di eventi mostrano come puoi usare tutta la potenza del linguaggio Nextflow nei tuoi script di workflow per aggiungere sofisticate capacità di logging e notifica.

---

## Riepilogo

Congratulazioni, ce l'hai fatta!

Durante questa quest secondaria, hai costruito un pipeline completo di elaborazione di campioni che è evoluto dalla gestione di base dei metadati a un workflow sofisticato e pronto per la produzione.
Ogni sezione ha costruito sulla precedente, dimostrando come i costrutti di programmazione trasformano semplici workflow in potenti sistemi di elaborazione dati, con i seguenti benefici:

- **Codice più chiaro**: Comprendere il flusso di dati vs lo scripting ti aiuta a scrivere workflow più organizzati
- **Gestione robusta**: Gli operatori di navigazione sicura ed Elvis rendono i workflow resistenti ai dati mancanti
- **Elaborazione flessibile**: La logica condizionale permette ai tuoi workflow di elaborare diversi tipi di campione in modo appropriato
- **Risorse adattive**: Le direttive dinamiche ottimizzano l'utilizzo delle risorse in base alle caratteristiche di input

Questa progressione riflette l'evoluzione nel mondo reale dei pipeline bioinformatici, dai prototipi di ricerca che gestiscono pochi campioni ai sistemi di produzione che elaborano migliaia di campioni tra laboratori e istituzioni.
Ogni sfida che hai risolto e ogni pattern che hai imparato riflette problemi reali che gli sviluppatori affrontano quando scalano i workflow Nextflow.

Applicare questi pattern nel tuo lavoro ti permetterà di costruire workflow robusti e pronti per la produzione.

### Pattern chiave

1.  **Flusso di dati vs Scripting:** Hai imparato a distinguere tra operazioni di flusso dati (orchestrazione di canali) e scripting (codice che manipola i dati), incluse le differenze cruciali tra operazioni su diversi tipi come `collect` su Channel vs List.

    - Flusso di dati: orchestrazione di canali

    ```groovy
    channel.fromPath('*.fastq').splitCsv(header: true)
    ```

    - Scripting: elaborazione dati su collezioni

    ```groovy
    sample_data.collect { it.toUpperCase() }
    ```

2.  **Elaborazione avanzata delle stringhe**: Hai padroneggiato le espressioni regolari per analizzare i nomi dei file, la generazione dinamica di script nei processi e l'interpolazione delle variabili (Nextflow vs Bash vs Shell).

    - Pattern matching

    ```groovy
    filename =~ ~/^(\w+)_(\w+)_(\d+)\.fastq$/
    ```

    - Funzione con ritorno condizionale

    ```groovy
    def parseSample(filename) {
        def matcher = filename =~ pattern
        return matcher ? [valid: true, data: matcher[0]] : [valid: false]
    }
    ```

    - Collezione di file ad argomenti di comando (nel blocco script di processo)

    ```groovy
    script:
    def file_args = input_files.collect { file -> "--input ${file}" }.join(' ')
    """
    analysis_tool ${file_args} --output results.txt
    """
    ```

3.  **Creazione di funzioni riutilizzabili**: Hai imparato a estrarre logica complessa in funzioni con nome che possono essere chiamate dagli operatori di canale, rendendo i workflow più leggibili e manutenibili.

    - Definire una funzione con nome

    ```groovy
    def separateMetadata(row) {
        def sample_meta = [ /* codice nascosto per brevità */ ]
        def fastq_path = file(row.file_path)
        def m = (fastq_path.name =~ /^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$/)
        def file_meta = m ? [ /* codice nascosto per brevità */ ] : [:]
        def priority = sample_meta.quality > 40 ? 'high' : 'normal'

        return tuple(sample_meta + file_meta + [priority: priority], fastq_path)
    }
    ```

    - Chiamare la funzione con nome in un workflow

    ```groovy
    workflow {
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map{ row -> separateMetadata(row) }

        ch_fastp = FASTP(ch_samples)
    }
    ```

4.  **Direttive di risorse dinamiche con Closure**: Hai esplorato l'uso di closure nelle direttive di processo per un'allocazione adattiva delle risorse basata sulle caratteristiche di input.

    - Closure con nome e composizione

    ```groovy
    def enrichData = normalizeId >> addQualityCategory >> addFlags
    def processor = generalFunction.curry(fixedParam)
    ```

    - Closure con accesso all'ambito

    ```groovy
    def collectStats = { data -> stats.count++; return data }
    ```

5.  **Logica condizionale e controllo del processo**: Hai aggiunto routing intelligente usando gli operatori `.branch()` e `.filter()`, sfruttando truthiness per espressioni condizionali concise.

    - Usa `.branch()` per instradare i dati attraverso diverse diramazioni del workflow

    ```groovy
    trim_branches = ch_samples
    .branch { meta, reads ->
        fastp: meta.organism == 'human' && meta.depth >= 30000000
        trimgalore: true
    }

    ch_fastp = FASTP(trim_branches.fastp)
    ch_trimgalore = TRIMGALORE(trim_branches.trimgalore)
    ```

    - Valutazione booleana con Groovy Truth

    ```groovy
    if (sample.files) println "Ha dei file"
    ```

    - Usa `filter()` per sottoinsiemi di dati con 'truthiness'

    ```groovy
    ch_valid_samples = ch_samples
        .filter { meta, reads ->
            meta.id && meta.organism && meta.depth >= 25000000
        }
    ```

6.  **Operatori di navigazione sicura ed Elvis**: Hai reso il pipeline robusto contro i dati mancanti usando `?.` per un accesso sicuro alle proprietà e `?:` per fornire valori predefiniti.

    ```groovy
    def id = data?.sample?.id ?: 'sconosciuto'
    ```

7.  **Validazione con error() e log.warn**: Hai imparato a validare gli input presto e fallire velocemente con messaggi di errore chiari.

    ```groovy
    try {
        def errors = validateSample(sample)
        if (errors) throw new RuntimeException("Non valido: ${errors.join(', ')}")
    } catch (Exception e) {
        println "Errore: ${e.message}"
    }
    ```

8.  **Gestori di eventi di configurazione**: Hai imparato a usare i gestori di eventi del workflow (`onComplete` e `onError`) per logging, notifiche e gestione del ciclo di vita.

    - Usare `onComplete` per loggare e notificare

    ```groovy
    workflow.onComplete = {
        println "Successo    : ${workflow.success}"
        println "stato uscita: ${workflow.exitStatus}"

        if (workflow.success) {
            println "✅ Pipeline completato con successo!"
        } else {
            println "❌ Pipeline fallito!"
            println "Errore: ${workflow.errorMessage}"
        }
    }
    ```

    - Usare `onError` per prendere azioni specificamente in caso di fallimento

    ```groovy
    workflow.onError = {
        // Scrivi log di errore dettagliato
        def error_file = file("${workflow.launchDir}/error.log")
        error_file.text = """
        Ora: ${new Date()}
        Errore: ${workflow.errorMessage}
        Rapporto di errore: ${workflow.errorReport ?: 'Nessun rapporto dettagliato disponibile'}
        """

        println "Dettagli dell'errore scritti in: ${error_file}"
    }
    ```

### Risorse aggiuntive

- [Riferimento del linguaggio Nextflow](https://nextflow.io/docs/latest/reference/syntax.html)
- [Operatori Nextflow](https://www.nextflow.io/docs/latest/operator.html)
- [Sintassi degli script Nextflow](https://www.nextflow.io/docs/latest/script.html)
- [Libreria standard di Nextflow](https://nextflow.io/docs/latest/reference/stdlib.html)

Assicurati di controllare queste risorse quando hai bisogno di esplorare funzionalità più avanzate.

Trarrai beneficio dalla pratica e dall'espansione delle tue abilità per:

- Scrivere workflow più puliti con una corretta separazione tra flusso dati e scripting
- Padroneggiare l'interpolazione delle variabili per evitare insidie comuni con variabili Nextflow, Bash e shell
- Utilizzare direttive di risorse dinamiche per workflow efficienti e adattivi
- Trasformare collezioni di file in argomenti di linea di comando formattati correttamente
- Gestire con eleganza diverse convenzioni di denominazione file e formati di input usando regex ed elaborazione delle stringhe
- Costruire codice riutilizzabile e manutenibile usando pattern avanzati di closure e programmazione funzionale
- Elaborare e organizzare dataset complessi utilizzando operazioni di collezione
- Aggiungere validazione, gestione degli errori e logging per rendere i tuoi workflow pronti per la produzione
- Implementare la gestione del ciclo di vita del workflow con gestori di eventi

---

## Cosa c'è dopo?

Torna al [menu delle Side Quests](./index.md) o clicca il pulsante in basso a destra della pagina per passare all'argomento successivo della lista.
