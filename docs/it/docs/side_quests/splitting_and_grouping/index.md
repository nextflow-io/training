# Suddivisione e Raggruppamento

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduzione assistita da IA - [scopri di più e suggerisci miglioramenti](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Nextflow fornisce strumenti potenti per lavorare con i dati in modo flessibile. Una capacità fondamentale è la suddivisione dei dati in diversi flussi e il successivo raggruppamento degli elementi correlati. Questo è particolarmente utile nei flussi di lavoro bioinformatici dove è necessario elaborare diversi tipi di campioni separatamente prima di combinare i risultati per l'analisi.

Pensatelo come smistare la posta: separate le lettere per destinazione, elaborate ogni pila in modo diverso, poi ricombinate gli elementi diretti alla stessa persona. Nextflow utilizza operatori speciali per realizzare questo con dati scientifici. Questo approccio è anche comunemente noto come pattern **scatter/gather** nel calcolo distribuito e nei flussi di lavoro bioinformatici.

Il sistema di canali di Nextflow è al cuore di questa flessibilità. I canali collegano diverse parti del vostro flusso di lavoro, permettendo ai dati di scorrere attraverso l'analisi. Potete creare più canali da una singola sorgente di dati, elaborare ogni canale in modo diverso, e poi unire i canali quando necessario. Questo approccio vi permette di progettare flussi di lavoro che rispecchiano naturalmente i percorsi ramificati e convergenti delle analisi bioinformatiche complesse.

### Obiettivi di apprendimento

In questa missione secondaria, imparerete a suddividere e raggruppare i dati usando gli operatori dei canali di Nextflow.
Inizieremo con un file CSV contenente informazioni sui campioni e i file di dati associati, poi manipoleremo e riorganizzeremo questi dati.

Al termine di questa missione secondaria, sarete in grado di separare e combinare flussi di dati efficacemente, usando le seguenti tecniche:

- Leggere dati da file usando `splitCsv`
- Filtrare e trasformare i dati con `filter` e `map`
- Combinare dati correlati usando `join` e `groupTuple`
- Creare combinazioni di dati con `combine` per l'elaborazione parallela
- Ottimizzare la struttura dei dati usando `subMap` e strategie di deduplicazione
- Costruire funzioni riutilizzabili con closure con nome per manipolare le strutture dei canali

Queste competenze vi aiuteranno a costruire flussi di lavoro in grado di gestire più file di input e diversi tipi di dati in modo efficiente, mantenendo al contempo una struttura del codice pulita e manutenibile.

### Prerequisiti

Prima di affrontare questa missione secondaria, dovreste:

- Aver completato il tutorial [Hello Nextflow](../hello_nextflow/README.md) o un corso equivalente per principianti.
- Essere a proprio agio con i concetti e i meccanismi di base di Nextflow (processi, canali, operatori, lavoro con i file, metadati)

**Facoltativo:** Raccomandiamo di completare prima la missione secondaria [Metadata in workflows](../metadata/).
Questa copre i fondamentali della lettura di file CSV con `splitCsv` e della creazione di meta map, che useremo ampiamente qui.

---

## 0. Iniziamo

#### Aprite il codespace di formazione

Se non lo avete ancora fatto, assicuratevi di aprire l'ambiente di formazione come descritto in [Environment Setup](../envsetup/index.md).

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Spostatevi nella directory del progetto

Spostiamoci nella directory dove si trovano i file per questo tutorial.

```bash
cd side-quests/splitting_and_grouping
```

Potete impostare VSCode per concentrarsi su questa directory:

```bash
code .
```

#### Esaminate i materiali

Troverete un file principale del flusso di lavoro e una directory `data` contenente un samplesheet chiamato `samplesheet.csv`.

```console title="Directory contents"
.
├── data
│   └── samplesheet.csv
└── main.nf
```

Il samplesheet contiene informazioni sui campioni di diversi pazienti, incluso l'ID del paziente, il numero di replica del campione, il tipo (normale o tumorale) e i percorsi a file di dati ipotetici (che in realtà non esistono, ma faremo finta che esistano).

```console title="samplesheet.csv"
id,repeat,type,bam
patientA,1,normal,patientA_rep1_normal.bam
patientA,1,tumor,patientA_rep1_tumor.bam
patientA,2,normal,patientA_rep2_normal.bam
patientA,2,tumor,patientA_rep2_tumor.bam
patientB,1,normal,patientB_rep1_normal.bam
patientB,1,tumor,patientB_rep1_tumor.bam
patientC,1,normal,patientC_rep1_normal.bam
patientC,1,tumor,patientC_rep1_tumor.bam
```

Questo samplesheet elenca otto campioni di tre pazienti (A, B, C).

Per ogni paziente, abbiamo campioni di tipo `tumor` (tipicamente provenienti da biopsie tumorali) o `normal` (prelevati da tessuto sano o sangue).
Se non avete familiarità con l'analisi del cancro, sappiate semplicemente che questo corrisponde a un modello sperimentale che utilizza campioni tumorali/normali accoppiati per eseguire analisi contrastive.

Per il paziente A in particolare, abbiamo due serie di repliche tecniche.

!!! note "Nota"

    Non preoccupatevi se non avete familiarità con questo disegno sperimentale, non è fondamentale per comprendere questo tutorial.

#### Esaminate il compito

La vostra sfida è scrivere un flusso di lavoro Nextflow che:

1. **Legga** i dati dei campioni da un file CSV e li strutturi con meta map
2. **Separi** i campioni in diversi canali in base al tipo (normale vs tumorale)
3. **Unisca** le coppie tumorale/normale corrispondenti per ID paziente e numero di replica
4. **Distribuisca** i campioni su intervalli genomici per l'elaborazione parallela
5. **Raggruppi** i campioni correlati per l'analisi a valle

Questo rappresenta un pattern bioinformatico comune dove è necessario suddividere i dati per l'elaborazione indipendente, poi ricombinare gli elementi correlati per l'analisi comparativa.

#### Lista di controllo per la preparazione

Pensate di essere pronti a tuffarvi?

- [ ] Comprendo l'obiettivo di questo corso e i suoi prerequisiti
- [ ] Il mio codespace è attivo e funzionante
- [ ] Ho impostato correttamente la mia directory di lavoro
- [ ] Comprendo il compito

Se potete spuntare tutte le caselle, siete pronti a partire.

---

## 1. Leggere i dati dei campioni

### 1.1. Leggere i dati dei campioni con `splitCsv` e creare meta map

Iniziamo leggendo i dati dei campioni con `splitCsv` e organizzandoli nel pattern della meta map. Nel file `main.nf`, vedrete che abbiamo già avviato il flusso di lavoro.

```groovy title="main.nf" linenums="1" hl_lines="2"
workflow {
    ch_samplesheet = channel.fromPath("./data/samplesheet.csv")
}
```

!!! note "Nota"

    In tutto questo tutorial, useremo il prefisso `ch_` per tutte le variabili dei canali per indicare chiaramente che sono canali Nextflow.

Se avete completato la missione secondaria [Metadata in workflows](../metadata/), riconoscerete questo pattern. Useremo `splitCsv` per leggere il CSV e strutturare immediatamente i dati con una meta map per separare i metadati dai percorsi dei file.

!!! info "Info"

    In questa formazione incontreremo due concetti diversi chiamati `map`:

    - **Struttura dati**: La map Groovy (equivalente ai dizionari/hash in altri linguaggi) che memorizza coppie chiave-valore
    - **Operatore del canale**: L'operatore `.map()` che trasforma gli elementi in un canale

    Chiariremo quale intendiamo nel contesto, ma questa distinzione è importante da capire quando si lavora con Nextflow.

Applicate queste modifiche a `main.nf`:

=== "Dopo"

    ```groovy title="main.nf" linenums="2" hl_lines="2-6"
        ch_samples = channel.fromPath("./data/samplesheet.csv")
            .splitCsv(header: true)
            .map{ row ->
              [[id:row.id, repeat:row.repeat, type:row.type], row.bam]
            }
            .view()
    ```

=== "Prima"

    ```groovy title="main.nf" linenums="2" hl_lines="1"
        ch_samplesheet = channel.fromPath("./data/samplesheet.csv")
    ```

Questo combina l'operazione `splitCsv` (lettura del CSV con intestazioni) e l'operazione `map` (strutturazione dei dati come tuple `[meta, file]`) in un unico passaggio. Applicate quella modifica ed eseguite la pipeline:

```bash
nextflow run main.nf
```

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [deadly_mercator] DSL2 - revision: bd6b0224e9

    [[id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam]
    [[id:patientA, repeat:1, type:tumor], patientA_rep1_tumor.bam]
    [[id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam]
    [[id:patientA, repeat:2, type:tumor], patientA_rep2_tumor.bam]
    [[id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam]
    [[id:patientB, repeat:1, type:tumor], patientB_rep1_tumor.bam]
    [[id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam]
    [[id:patientC, repeat:1, type:tumor], patientC_rep1_tumor.bam]
    ```

Abbiamo ora un canale dove ogni elemento è una tupla `[meta, file]` - i metadati separati dai percorsi dei file. Questa struttura ci permette di suddividere e raggruppare il nostro carico di lavoro in base ai campi dei metadati.

---

## 2. Filtrare e trasformare i dati

### 2.1. Filtrare i dati con `filter`

Possiamo usare l'[operatore `filter`](https://www.nextflow.io/docs/latest/operator.html#filter) per filtrare i dati in base a una condizione. Supponiamo di voler elaborare solo i campioni normali. Possiamo farlo filtrando i dati in base al campo `type`. Inseriamo questo prima dell'operatore `view`.

=== "Dopo"

    ```groovy title="main.nf" linenums="2" hl_lines="6"
        ch_samples = channel.fromPath("./data/samplesheet.csv")
            .splitCsv(header: true)
            .map{ row ->
              [[id:row.id, repeat:row.repeat, type:row.type], row.bam]
            }
            .filter { meta, file -> meta.type == 'normal' }
            .view()
    ```

=== "Prima"

    ```groovy title="main.nf" linenums="2"
        ch_samples = channel.fromPath("./data/samplesheet.csv")
            .splitCsv(header: true)
            .map{ row ->
              [[id:row.id, repeat:row.repeat, type:row.type], row.bam]
            }
            .view()
    ```

Eseguite di nuovo il flusso di lavoro per vedere il risultato filtrato:

```bash
nextflow run main.nf
```

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [admiring_brown] DSL2 - revision: 194d61704d

    [[id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam]
    [[id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam]
    [[id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam]
    [[id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam]
    ```

Abbiamo filtrato con successo i dati per includere solo i campioni normali. Ricapitoliamo come funziona.

L'operatore `filter` accetta una closure che viene applicata a ogni elemento nel canale. Se la closure restituisce `true`, l'elemento viene incluso; se restituisce `false`, l'elemento viene escluso.

Nel nostro caso, vogliamo mantenere solo i campioni dove `meta.type == 'normal'`. La closure usa la tupla `meta,file` per riferirsi a ogni campione, accede al tipo di campione con `meta.type`, e verifica se è uguale a `'normal'`.

Questo viene realizzato con la singola closure che abbiamo introdotto sopra:

```groovy title="main.nf" linenums="7"
    .filter { meta, file -> meta.type == 'normal' }
```

### 2.2. Creare canali filtrati separati

Attualmente stiamo applicando il filtro al canale creato direttamente dal CSV, ma vogliamo filtrarlo in più modi, quindi riscriviamo la logica per creare un canale filtrato separato per i campioni normali:

=== "Dopo"

    ```groovy title="main.nf" linenums="2" hl_lines="6 8"
        ch_samples = channel.fromPath("./data/samplesheet.csv")
            .splitCsv(header: true)
            .map{ row ->
                [[id:row.id, repeat:row.repeat, type:row.type], row.bam]
            }
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
        ch_normal_samples
            .view()
    ```

=== "Prima"

    ```groovy title="main.nf" linenums="2"
        ch_samples = channel.fromPath("./data/samplesheet.csv")
            .splitCsv(header: true)
            .map{ row ->
              [[id:row.id, repeat:row.repeat, type:row.type], row.bam]
            }
            .filter { meta, file -> meta.type == 'normal' }
            .view()
    ```

Eseguite la pipeline per vedere i risultati:

```bash
nextflow run main.nf
```

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [trusting_poisson] DSL2 - revision: 639186ee74

    [[id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam]
    [[id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam]
    [[id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam]
    [[id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam]
    ```

Abbiamo filtrato con successo i dati e creato un canale separato per i campioni normali.

Creiamo anche un canale filtrato per i campioni tumorali:

=== "Dopo"

    ```groovy title="main.nf" linenums="7" hl_lines="3-8"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
        ch_normal_samples
            .view{'Normal sample: ' + it}
        ch_tumor_samples
            .view{'Tumor sample: ' + it}
    ```

=== "Prima"

    ```groovy title="main.nf" linenums="7" hl_lines="3 4"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
        ch_normal_samples
            .view()
    ```

```bash
nextflow run main.nf
```

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [maniac_boltzmann] DSL2 - revision: 3636b6576b

    Tumor sample: [[id:patientA, repeat:1, type:tumor], patientA_rep1_tumor.bam]
    Tumor sample: [[id:patientA, repeat:2, type:tumor], patientA_rep2_tumor.bam]
    Normal sample: [[id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam]
    Normal sample: [[id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam]
    Normal sample: [[id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam]
    Normal sample: [[id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam]
    Tumor sample: [[id:patientB, repeat:1, type:tumor], patientB_rep1_tumor.bam]
    Tumor sample: [[id:patientC, repeat:1, type:tumor], patientC_rep1_tumor.bam]
    ```

Abbiamo separato i campioni normali e tumorali in due canali diversi, e usato una closure fornita a `view()` per etichettarli diversamente nell'output: `ch_tumor_samples.view{'Tumor sample: ' + it}`.

### Takeaway

In questa sezione, avete imparato:

- **Filtrare i dati**: Come filtrare i dati con `filter`
- **Suddividere i dati**: Come suddividere i dati in diversi canali in base a una condizione
- **Visualizzare i dati**: Come usare `view` per stampare i dati ed etichettare l'output da diversi canali

Abbiamo ora separato i campioni normali e tumorali in due canali diversi. Successivamente, uniremo i campioni normali e tumorali sul campo `id`.

---

## 3. Unire i canali tramite identificatori

Nella sezione precedente, abbiamo separato i campioni normali e tumorali in due canali diversi. Questi potrebbero essere elaborati indipendentemente usando processi o flussi di lavoro specifici in base al loro tipo. Ma cosa succede quando vogliamo confrontare i campioni normali e tumorali dello stesso paziente? A questo punto, dobbiamo unirli di nuovo assicurandoci di abbinare i campioni in base al campo `id`.

Nextflow include molti metodi per combinare i canali, ma in questo caso l'operatore più appropriato è [`join`](https://www.nextflow.io/docs/latest/operator.html#join). Se avete familiarità con SQL, si comporta come l'operazione `JOIN`, dove specifichiamo la chiave su cui unire e il tipo di unione da eseguire.

### 3.1. Usare `map` e `join` per combinare in base all'ID paziente

Se consultiamo la documentazione di [`join`](https://www.nextflow.io/docs/latest/operator.html#join), possiamo vedere che per impostazione predefinita unisce due canali in base al primo elemento di ogni tupla.

#### 3.1.1. Verificare la struttura dei dati

Se non avete ancora l'output della console disponibile, eseguiamo la pipeline per verificare la struttura dei dati e vedere come dobbiamo modificarla per unire sul campo `id`.

```bash
nextflow run main.nf
```

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [maniac_boltzmann] DSL2 - revision: 3636b6576b

    Tumor sample: [[id:patientA, repeat:1, type:tumor], patientA_rep1_tumor.bam]
    Tumor sample: [[id:patientA, repeat:2, type:tumor], patientA_rep2_tumor.bam]
    Normal sample: [[id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam]
    Normal sample: [[id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam]
    Normal sample: [[id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam]
    Normal sample: [[id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam]
    Tumor sample: [[id:patientB, repeat:1, type:tumor], patientB_rep1_tumor.bam]
    Tumor sample: [[id:patientC, repeat:1, type:tumor], patientC_rep1_tumor.bam]
    ```

Possiamo vedere che il campo `id` è il primo elemento in ogni meta map. Affinché `join` funzioni, dovremmo isolare il campo `id` in ogni tupla. Dopodiché, possiamo semplicemente usare l'operatore `join` per combinare i due canali.

#### 3.1.2. Isolare il campo `id`

Per isolare il campo `id`, possiamo usare l'[operatore `map`](https://www.nextflow.io/docs/latest/operator.html#map) per creare una nuova tupla con il campo `id` come primo elemento.

=== "Dopo"

    ```groovy title="main.nf" linenums="7" hl_lines="3 6"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
            .map { meta, file -> [meta.id, meta, file] }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
            .map { meta, file -> [meta.id, meta, file] }
        ch_normal_samples
            .view{'Normal sample: ' + it}
        ch_tumor_samples
            .view{'Tumor sample: ' + it}
    ```

=== "Prima"

    ```groovy title="main.nf" linenums="7"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
        ch_normal_samples
            .view{'Normal sample: ' + it}
        ch_tumor_samples
            .view{'Tumor sample: ' + it}
    ```

```bash
nextflow run main.nf
```

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [mad_lagrange] DSL2 - revision: 9940b3f23d

    Tumor sample: [patientA, [id:patientA, repeat:1, type:tumor], patientA_rep1_tumor.bam]
    Tumor sample: [patientA, [id:patientA, repeat:2, type:tumor], patientA_rep2_tumor.bam]
    Normal sample: [patientA, [id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam]
    Normal sample: [patientA, [id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam]
    Tumor sample: [patientB, [id:patientB, repeat:1, type:tumor], patientB_rep1_tumor.bam]
    Tumor sample: [patientC, [id:patientC, repeat:1, type:tumor], patientC_rep1_tumor.bam]
    Normal sample: [patientB, [id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam]
    Normal sample: [patientC, [id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam]
    ```

Potrebbe essere sottile, ma dovreste riuscire a vedere che il primo elemento in ogni tupla è il campo `id`.

#### 3.1.3. Combinare i due canali

Ora possiamo usare l'operatore `join` per combinare i due canali in base al campo `id`.

Ancora una volta, useremo `view` per stampare gli output uniti.

=== "Dopo"

    ```groovy title="main.nf" linenums="7" hl_lines="7-9"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
            .map { meta, file -> [meta.id, meta, file] }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
            .map { meta, file -> [meta.id, meta, file] }
        ch_joined_samples = ch_normal_samples
            .join(ch_tumor_samples)
        ch_joined_samples.view()
    ```

=== "Prima"

    ```groovy title="main.nf" linenums="7" hl_lines="7-10"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
            .map { meta, file -> [meta.id, meta, file] }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
            .map { meta, file -> [meta.id, meta, file] }
        ch_normal_samples
            .view{'Normal sample: ' + it}
        ch_tumor_samples
            .view{'Tumor sample: ' + it}
    ```

```bash
nextflow run main.nf
```

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [soggy_wiles] DSL2 - revision: 3bc1979889

    [patientA, [id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam, [id:patientA, repeat:1, type:tumor], patientA_rep1_tumor.bam]
    [patientA, [id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam, [id:patientA, repeat:2, type:tumor], patientA_rep2_tumor.bam]
    [patientB, [id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam, [id:patientB, repeat:1, type:tumor], patientB_rep1_tumor.bam]
    [patientC, [id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam, [id:patientC, repeat:1, type:tumor], patientC_rep1_tumor.bam]
    ```

È un po' difficile da leggere perché è molto largo, ma dovreste riuscire a vedere che i campioni sono stati uniti in base al campo `id`. Ogni tupla ha ora il formato:

- `id`: L'ID del campione
- `normal_meta_map`: I metadati del campione normale inclusi tipo, replica e percorso al file bam
- `normal_sample_file`: Il file del campione normale
- `tumor_meta_map`: I metadati del campione tumorale inclusi tipo, replica e percorso al file bam
- `tumor_sample`: Il campione tumorale inclusi tipo, replica e percorso al file bam

!!! warning "Avviso"

    L'operatore `join` scarterà le tuple non abbinate. In questo esempio, ci siamo assicurati che tutti i campioni fossero abbinati per tumorale e normale, ma se questo non è vero dovete usare il parametro `remainder: true` per mantenere le tuple non abbinate. Consultate la [documentazione](https://www.nextflow.io/docs/latest/operator.html#join) per maggiori dettagli.

Ora sapete come usare `map` per isolare un campo in una tupla, e come usare `join` per combinare tuple in base al primo campo.
Con questa conoscenza, possiamo combinare con successo i canali in base a un campo condiviso.

Successivamente, considereremo la situazione in cui si vuole unire su più campi.

### 3.2. Unire su più campi

Abbiamo 2 repliche per il campione A, ma solo 1 per i campioni B e C. In questo caso siamo riusciti a unirli efficacemente usando il campo `id`, ma cosa succederebbe se fossero fuori sincronia? Potremmo confondere i campioni normali e tumorali di repliche diverse!

Per evitare questo, possiamo unire su più campi. In realtà ci sono più modi per ottenere questo risultato, ma ci concentreremo sulla creazione di una nuova chiave di unione che includa sia l'`id` del campione che il numero di `replicate`.

Iniziamo creando una nuova chiave di unione. Possiamo farlo nello stesso modo di prima usando l'[operatore `map`](https://www.nextflow.io/docs/latest/operator.html#map) per creare una nuova tupla con i campi `id` e `repeat` come primo elemento.

=== "Dopo"

    ```groovy title="main.nf" linenums="7" hl_lines="3 6"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
            .map { meta, file -> [[meta.id, meta.repeat], meta, file] }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
            .map { meta, file -> [[meta.id, meta.repeat], meta, file] }
    ```

=== "Prima"

    ```groovy title="main.nf" linenums="7" hl_lines="3 6"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
            .map { meta, file -> [meta.id, meta, file] }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
            .map { meta, file -> [meta.id, meta, file] }
    ```

Ora dovremmo vedere che l'unione avviene usando sia i campi `id` che `repeat`. Eseguite il flusso di lavoro:

```bash
nextflow run main.nf
```

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [prickly_wing] DSL2 - revision: 3bebf22dee

    [[patientA, 1], [id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam, [id:patientA, repeat:1, type:tumor], patientA_rep1_tumor.bam]
    [[patientA, 2], [id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam, [id:patientA, repeat:2, type:tumor], patientA_rep2_tumor.bam]
    [[patientB, 1], [id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam, [id:patientB, repeat:1, type:tumor], patientB_rep1_tumor.bam]
    [[patientC, 1], [id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam, [id:patientC, repeat:1, type:tumor], patientC_rep1_tumor.bam]
    ```

Notate come abbiamo una tupla di due elementi (campi `id` e `repeat`) come primo elemento di ogni risultato unito. Questo dimostra come elementi complessi possano essere usati come chiave di unione, consentendo un abbinamento piuttosto intricato tra campioni delle stesse condizioni.

Se volete esplorare altri modi per unire su chiavi diverse, consultate la [documentazione dell'operatore join](https://www.nextflow.io/docs/latest/operator.html#join) per opzioni ed esempi aggiuntivi.

### 3.3. Usare `subMap` per creare una nuova chiave di unione

L'approccio precedente perde i nomi dei campi dalla nostra chiave di unione - i campi `id` e `repeat` diventano solo un elenco di valori. Per mantenere i nomi dei campi per un accesso successivo, possiamo usare il [metodo `subMap`](<https://docs.groovy-lang.org/latest/html/groovy-jdk/java/util/Map.html#subMap(java.util.Collection)>).

Il metodo `subMap` estrae solo le coppie chiave-valore specificate da una map. Qui estrarremo solo i campi `id` e `repeat` per creare la nostra chiave di unione.

=== "Dopo"

    ```groovy title="main.nf" linenums="7" hl_lines="3 6"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
            .map { meta, file -> [meta.subMap(['id', 'repeat']), meta, file] }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
            .map { meta, file -> [meta.subMap(['id', 'repeat']), meta, file] }
    ```

=== "Prima"

    ```groovy title="main.nf" linenums="7" hl_lines="3 6"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
            .map { meta, file -> [[meta.id, meta.repeat], meta, file] }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
            .map { meta, file -> [[meta.id, meta.repeat], meta, file] }
    ```

```bash
nextflow run main.nf
```

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [reverent_wing] DSL2 - revision: 847016c3b7

    [[id:patientA, repeat:1], [id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam, [id:patientA, repeat:1, type:tumor], patientA_rep1_tumor.bam]
    [[id:patientA, repeat:2], [id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam, [id:patientA, repeat:2, type:tumor], patientA_rep2_tumor.bam]
    [[id:patientB, repeat:1], [id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam, [id:patientB, repeat:1, type:tumor], patientB_rep1_tumor.bam]
    [[id:patientC, repeat:1], [id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam, [id:patientC, repeat:1, type:tumor], patientC_rep1_tumor.bam]
    ```

Ora abbiamo una nuova chiave di unione che non solo include i campi `id` e `repeat` ma mantiene anche i nomi dei campi in modo da potervi accedere in seguito per nome, ad esempio `meta.id` e `meta.repeat`.

### 3.4. Usare una closure con nome in map

Per evitare duplicazioni e ridurre gli errori, possiamo usare una closure con nome. Una closure con nome ci permette di creare una funzione riutilizzabile che possiamo chiamare in più punti.

Per farlo, prima definiamo la closure come una nuova variabile:

=== "Dopo"

    ```groovy title="main.nf" linenums="2" hl_lines="7"
        ch_samples = channel.fromPath("./data/samplesheet.csv")
            .splitCsv(header: true)
            .map{ row ->
                [[id:row.id, repeat:row.repeat, type:row.type], row.bam]
            }

        getSampleIdAndReplicate = { meta, bam -> [ meta.subMap(['id', 'repeat']), meta, file(bam) ] }

        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
    ```

=== "Prima"

    ```groovy title="main.nf" linenums="2"
        ch_samples = channel.fromPath("./data/samplesheet.csv")
            .splitCsv(header: true)
            .map{ row ->
                [[id:row.id, repeat:row.repeat, type:row.type], row.bam]
            }
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
    ```

Abbiamo definito la trasformazione map come una variabile con nome che possiamo riutilizzare.

Notate che convertiamo anche il percorso del file in un oggetto Path usando `file()` in modo che qualsiasi processo che riceve questo canale possa gestire il file correttamente (per maggiori informazioni vedere [Working with files](../working_with_files/)).

Implementiamo la closure nel nostro flusso di lavoro:

=== "Dopo"

    ```groovy title="main.nf" linenums="10" hl_lines="3 6"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
             .map ( getSampleIdAndReplicate )
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
             .map ( getSampleIdAndReplicate )

    ```

=== "Prima"

    ```groovy title="main.nf" linenums="10" hl_lines="3 6"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
            .map { meta, file -> [meta.subMap(['id', 'repeat']), meta, file] }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
            .map { meta, file -> [meta.subMap(['id', 'repeat']), meta, file] }
    ```

!!! note "Nota"

    L'operatore `map` è passato dall'uso di `{ }` all'uso di `( )` per passare la closure come argomento. Questo perché l'operatore `map` si aspetta una closure come argomento e `{ }` viene usato per definire una closure anonima. Quando si chiama una closure con nome, usate la sintassi `( )`.

Eseguite di nuovo il flusso di lavoro per verificare che tutto funzioni ancora:

```bash
nextflow run main.nf
```

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [angry_meninsky] DSL2 - revision: 2edc226b1d

    [[id:patientA, repeat:1], [id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam, [id:patientA, repeat:1, type:tumor], patientA_rep1_tumor.bam]
    [[id:patientA, repeat:2], [id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam, [id:patientA, repeat:2, type:tumor], patientA_rep2_tumor.bam]
    [[id:patientB, repeat:1], [id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam, [id:patientB, repeat:1, type:tumor], patientB_rep1_tumor.bam]
    [[id:patientC, repeat:1], [id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam, [id:patientC, repeat:1, type:tumor], patientC_rep1_tumor.bam]
    ```

Usare una closure con nome ci permette di riutilizzare la stessa trasformazione in più punti, riducendo il rischio di errori e rendendo il codice più leggibile e manutenibile.

### 3.5. Ridurre la duplicazione dei dati

Abbiamo molti dati duplicati nel nostro flusso di lavoro. Ogni elemento nei campioni uniti ripete i campi `id` e `repeat`. Poiché queste informazioni sono già disponibili nella chiave di raggruppamento, possiamo evitare questa ridondanza. Come promemoria, la nostra struttura dati attuale è simile a questa:

```groovy
[
  [
    "id": "sampleC",
    "repeat": "1",
  ],
  [
    "id": "sampleC",
    "repeat": "1",
    "type": "normal",
  ],
  "sampleC_rep1_normal.bam"
  [
    "id": "sampleC",
    "repeat": "1",
    "type": "tumor",
  ],
  "sampleC_rep1_tumor.bam"
]
```

Poiché i campi `id` e `repeat` sono disponibili nella chiave di raggruppamento, rimuoviamoli dal resto di ogni elemento del canale per evitare duplicazioni. Possiamo farlo usando il metodo `subMap` per creare una nuova map con solo il campo `type`. Questo approccio ci permette di mantenere tutte le informazioni necessarie eliminando la ridondanza nella nostra struttura dati.

=== "Dopo"

    ```groovy title="main.nf" linenums="8" hl_lines="1"
        getSampleIdAndReplicate = { meta, bam -> [ meta.subMap(['id', 'repeat']), meta.subMap(['type']), file(bam) ] }
    ```

=== "Prima"

    ```groovy title="main.nf" linenums="8" hl_lines="1"
        getSampleIdAndReplicate = { meta, bam -> [ meta.subMap(['id', 'repeat']), meta, file(bam) ] }
    ```

Ora la closure restituisce una tupla dove il primo elemento contiene i campi `id` e `repeat`, e il secondo elemento contiene solo il campo `type`. Questo elimina la ridondanza memorizzando le informazioni `id` e `repeat` una sola volta nella chiave di raggruppamento, mantenendo al contempo tutte le informazioni necessarie.

Eseguite il flusso di lavoro per vedere come appare:

```bash
nextflow run main.nf
```

??? success "Output del comando"

    ```console
    [[id:patientA, repeat:1], [type:normal], /workspaces/training/side-quests/splitting_and_grouping/patientA_rep1_normal.bam, [type:tumor], /workspaces/training/side-quests/splitting_and_grouping/patientA_rep1_tumor.bam]
    [[id:patientA, repeat:2], [type:normal], /workspaces/training/side-quests/splitting_and_grouping/patientA_rep2_normal.bam, [type:tumor], /workspaces/training/side-quests/splitting_and_grouping/patientA_rep2_tumor.bam]
    [[id:patientB, repeat:1], [type:normal], /workspaces/training/side-quests/splitting_and_grouping/patientB_rep1_normal.bam, [type:tumor], /workspaces/training/side-quests/splitting_and_grouping/patientB_rep1_tumor.bam]
    [[id:patientC, repeat:1], [type:normal], /workspaces/training/side-quests/splitting_and_grouping/patientC_rep1_normal.bam, [type:tumor], /workspaces/training/side-quests/splitting_and_grouping/patientC_rep1_tumor.bam]
    ```

Possiamo vedere che indichiamo i campi `id` e `repeat` solo una volta nella chiave di raggruppamento e abbiamo il campo `type` nei dati del campione. Non abbiamo perso alcuna informazione ma siamo riusciti a rendere il contenuto del nostro canale più conciso.

### 3.6. Rimuovere le informazioni ridondanti

Abbiamo rimosso le informazioni duplicate sopra, ma abbiamo ancora altre informazioni ridondanti nei nostri canali.

All'inizio, abbiamo separato i campioni normali e tumorali usando `filter`, poi li abbiamo uniti in base alle chiavi `id` e `repeat`. L'operatore `join` preserva l'ordine in cui le tuple vengono unite, quindi nel nostro caso, con i campioni normali sul lato sinistro e i campioni tumorali sul lato destro, il canale risultante mantiene questa struttura: `id, <elementi normali>, <elementi tumorali>`.

Poiché conosciamo la posizione di ogni elemento nel nostro canale, possiamo semplificare ulteriormente la struttura eliminando i metadati `[type:normal]` e `[type:tumor]`.

=== "Dopo"

    ```groovy title="main.nf" linenums="8" hl_lines="1"
        getSampleIdAndReplicate = { meta, file -> [ meta.subMap(['id', 'repeat']), file ] }
    ```

=== "Prima"

    ```groovy title="main.nf" linenums="8" hl_lines="1"
        getSampleIdAndReplicate = { meta, file -> [ meta.subMap(['id', 'repeat']), meta.subMap(['type']), file ] }
    ```

Eseguite di nuovo per vedere il risultato:

```bash
nextflow run main.nf
```

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [confident_leavitt] DSL2 - revision: a2303895bd

    [[id:patientA, repeat:1], patientA_rep1_normal.bam, patientA_rep1_tumor.bam]
    [[id:patientA, repeat:2], patientA_rep2_normal.bam, patientA_rep2_tumor.bam]
    [[id:patientB, repeat:1], patientB_rep1_normal.bam, patientB_rep1_tumor.bam]
    [[id:patientC, repeat:1], patientC_rep1_normal.bam, patientC_rep1_tumor.bam]
    ```

### Takeaway

In questa sezione, avete imparato:

- **Manipolare le tuple**: Come usare `map` per isolare un campo in una tupla
- **Unire le tuple**: Come usare `join` per combinare tuple in base al primo campo
- **Creare chiavi di unione**: Come usare `subMap` per creare una nuova chiave di unione
- **Closure con nome**: Come usare una closure con nome in map
- **Unione su più campi**: Come unire su più campi per un abbinamento più preciso
- **Ottimizzazione della struttura dati**: Come semplificare la struttura del canale rimuovendo le informazioni ridondanti

Avete ora un flusso di lavoro che può suddividere un samplesheet, filtrare i campioni normali e tumorali, unirli per ID campione e numero di replica, poi stampare i risultati.

Questo è un pattern comune nei flussi di lavoro bioinformatici dove è necessario abbinare campioni o altri tipi di dati dopo l'elaborazione indipendente, quindi è una competenza utile. Successivamente, vedremo come ripetere un campione più volte.

## 4. Distribuire i campioni su intervalli

Un pattern chiave nei flussi di lavoro bioinformatici è la distribuzione dell'analisi su regioni genomiche. Ad esempio, la chiamata delle varianti può essere parallelizzata dividendo il genoma in intervalli (come cromosomi o regioni più piccole). Questa strategia di parallelizzazione migliora significativamente l'efficienza della pipeline distribuendo il carico computazionale su più core o nodi, riducendo il tempo di esecuzione complessivo.

Nella sezione seguente, dimostreremo come distribuire i dati dei nostri campioni su più intervalli genomici. Abbineremo ogni campione a ogni intervallo, consentendo l'elaborazione parallela di diverse regioni genomiche. Questo moltiplicherà la dimensione del nostro dataset per il numero di intervalli, creando più unità di analisi indipendenti che potranno essere riunite in seguito.

### 4.1. Distribuire i campioni sugli intervalli usando `combine`

Iniziamo creando un canale di intervalli. Per semplicità, useremo solo 3 intervalli che definiremo manualmente. In un flusso di lavoro reale, potreste leggerli da un file di input o persino creare un canale con molti file di intervalli.

=== "Dopo"

    ```groovy title="main.nf" linenums="17" hl_lines="2"
            .join(ch_tumor_samples)
        ch_intervals = channel.of('chr1', 'chr2', 'chr3')
    ```

=== "Prima"

    ```groovy title="main.nf" linenums="17" hl_lines="2"
            .join(ch_tumor_samples)
        ch_joined_samples.view()
    ```

Ricordate, vogliamo ripetere ogni campione per ogni intervallo. Questo è talvolta chiamato prodotto cartesiano dei campioni e degli intervalli. Possiamo ottenerlo usando l'[operatore `combine`](https://www.nextflow.io/docs/latest/operator.html#combine). Questo prenderà ogni elemento dal canale 1 e lo ripeterà per ogni elemento nel canale 2. Aggiungiamo un operatore combine al nostro flusso di lavoro:

=== "Dopo"

    ```groovy title="main.nf" linenums="18" hl_lines="3-5"
        ch_intervals = channel.of('chr1', 'chr2', 'chr3')

        ch_combined_samples = ch_joined_samples
            .combine(ch_intervals)
            .view()
    ```

=== "Prima"

    ```groovy title="main.nf" linenums="18"
        ch_intervals = channel.of('chr1', 'chr2', 'chr3')
    ```

Eseguiamolo e vediamo cosa succede:

```bash
nextflow run main.nf
```

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [mighty_tesla] DSL2 - revision: ae013ab70b

    [[id:patientA, repeat:1], patientA_rep1_normal.bam, patientA_rep1_tumor.bam, chr1]
    [[id:patientA, repeat:1], patientA_rep1_normal.bam, patientA_rep1_tumor.bam, chr2]
    [[id:patientA, repeat:1], patientA_rep1_normal.bam, patientA_rep1_tumor.bam, chr3]
    [[id:patientA, repeat:2], patientA_rep2_normal.bam, patientA_rep2_tumor.bam, chr1]
    [[id:patientA, repeat:2], patientA_rep2_normal.bam, patientA_rep2_tumor.bam, chr2]
    [[id:patientA, repeat:2], patientA_rep2_normal.bam, patientA_rep2_tumor.bam, chr3]
    [[id:patientB, repeat:1], patientB_rep1_normal.bam, patientB_rep1_tumor.bam, chr1]
    [[id:patientB, repeat:1], patientB_rep1_normal.bam, patientB_rep1_tumor.bam, chr2]
    [[id:patientB, repeat:1], patientB_rep1_normal.bam, patientB_rep1_tumor.bam, chr3]
    [[id:patientC, repeat:1], patientC_rep1_normal.bam, patientC_rep1_tumor.bam, chr1]
    [[id:patientC, repeat:1], patientC_rep1_normal.bam, patientC_rep1_tumor.bam, chr2]
    [[id:patientC, repeat:1], patientC_rep1_normal.bam, patientC_rep1_tumor.bam, chr3]
    ```

Ottimo! Abbiamo ripetuto ogni campione per ogni singolo intervallo nella nostra lista di 3 intervalli. Abbiamo effettivamente triplicato il numero di elementi nel nostro canale.

È un po' difficile da leggere però, quindi nella prossima sezione lo sistemeremo.

### 4.2. Organizzare il canale

Possiamo usare l'operatore `map` per riordinare e rifattorizzare i dati dei nostri campioni in modo che siano più facili da capire. Spostiamo la stringa degli intervalli nella map di unione come primo elemento.

=== "Dopo"

    ```groovy title="main.nf" linenums="20" hl_lines="3-9"
        ch_combined_samples = ch_joined_samples
            .combine(ch_intervals)
            .map { grouping_key, normal, tumor, interval ->
                [
                    grouping_key + [interval: interval],
                    normal,
                    tumor
                ]
            }
            .view()
    ```

=== "Prima"

    ```groovy title="main.nf" linenums="20"
        ch_combined_samples = ch_joined_samples
            .combine(ch_intervals)
            .view()
    ```

Analizziamo passo per passo cosa fa questa operazione map.

Prima, usiamo parametri con nome per rendere il codice più leggibile. Usando i nomi `grouping_key`, `normal`, `tumor` e `interval`, possiamo riferirci agli elementi nella tupla per nome invece che per indice:

```groovy
        .map { grouping_key, normal, tumor, interval ->
```

Successivamente, combiniamo la `grouping_key` con il campo `interval`. La `grouping_key` è una map contenente i campi `id` e `repeat`. Creiamo una nuova map con l'`interval` e le uniamo usando l'addizione di map di Groovy (`+`):

```groovy
                grouping_key + [interval: interval],
```

Infine, restituiamo questo come una tupla con tre elementi: la map dei metadati combinata, il file del campione normale e il file del campione tumorale:

```groovy
            [
                grouping_key + [interval: interval],
                normal,
                tumor
            ]
```

Eseguiamolo di nuovo e verifichiamo il contenuto del canale:

```bash
nextflow run main.nf
```

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [sad_hawking] DSL2 - revision: 1f6f6250cd

    [[id:patientA, repeat:1, interval:chr1], patientA_rep1_normal.bam, patientA_rep1_tumor.bam]
    [[id:patientA, repeat:1, interval:chr2], patientA_rep1_normal.bam, patientA_rep1_tumor.bam]
    [[id:patientA, repeat:1, interval:chr3], patientA_rep1_normal.bam, patientA_rep1_tumor.bam]
    [[id:patientA, repeat:2, interval:chr1], patientA_rep2_normal.bam, patientA_rep2_tumor.bam]
    [[id:patientA, repeat:2, interval:chr2], patientA_rep2_normal.bam, patientA_rep2_tumor.bam]
    [[id:patientA, repeat:2, interval:chr3], patientA_rep2_normal.bam, patientA_rep2_tumor.bam]
    [[id:patientB, repeat:1, interval:chr1], patientB_rep1_normal.bam, patientB_rep1_tumor.bam]
    [[id:patientB, repeat:1, interval:chr2], patientB_rep1_normal.bam, patientB_rep1_tumor.bam]
    [[id:patientB, repeat:1, interval:chr3], patientB_rep1_normal.bam, patientB_rep1_tumor.bam]
    [[id:patientC, repeat:1, interval:chr1], patientC_rep1_normal.bam, patientC_rep1_tumor.bam]
    [[id:patientC, repeat:1, interval:chr2], patientC_rep1_normal.bam, patientC_rep1_tumor.bam]
    [[id:patientC, repeat:1, interval:chr3], patientC_rep1_normal.bam, patientC_rep1_tumor.bam]
    ```

Usare `map` per trasformare i dati nella struttura corretta può essere complicato, ma è fondamentale per una manipolazione efficace dei dati.

Abbiamo ora ogni campione ripetuto su tutti gli intervalli genomici, creando più unità di analisi indipendenti che possono essere elaborate in parallelo. Ma cosa succede se vogliamo riunire i campioni correlati? Nella prossima sezione, impareremo come raggruppare i campioni che condividono attributi comuni.

### Takeaway

In questa sezione, avete imparato:

- **Distribuire i campioni sugli intervalli**: Come usare `combine` per ripetere i campioni sugli intervalli
- **Creare prodotti cartesiani**: Come generare tutte le combinazioni di campioni e intervalli
- **Organizzare la struttura del canale**: Come usare `map` per ristrutturare i dati per una migliore leggibilità
- **Preparazione all'elaborazione parallela**: Come impostare i dati per l'analisi distribuita

## 5. Aggregare i campioni usando `groupTuple`

Nelle sezioni precedenti, abbiamo imparato come suddividere i dati da un file di input e filtrare per campi specifici (nel nostro caso campioni normali e tumorali). Ma questo copre solo un singolo tipo di unione. Cosa succede se vogliamo raggruppare i campioni per un attributo specifico? Ad esempio, invece di unire coppie normale-tumorale corrispondenti, potremmo voler elaborare tutti i campioni di "sampleA" insieme indipendentemente dal loro tipo. Questo pattern è comune nei flussi di lavoro bioinformatici dove potreste voler elaborare campioni correlati separatamente per ragioni di efficienza prima di confrontare o combinare i risultati alla fine.

Nextflow include metodi integrati per farlo, il principale che esamineremo è `groupTuple`.

Iniziamo raggruppando tutti i nostri campioni che hanno gli stessi campi `id` e `interval`; questo sarebbe tipico di un'analisi dove volessimo raggruppare le repliche tecniche ma mantenere separati i campioni significativamente diversi.

Per farlo, dovremmo separare le nostre variabili di raggruppamento in modo da poterle usare in isolamento.

Il primo passo è simile a quello che abbiamo fatto nella sezione precedente. Dobbiamo isolare la nostra variabile di raggruppamento come primo elemento della tupla. Ricordate, il nostro primo elemento è attualmente una map dei campi `id`, `repeat` e `interval`:

```groovy title="main.nf" linenums="1"
{
  "id": "sampleA",
  "repeat": "1",
  "interval": "chr1"
}
```

Possiamo riutilizzare il metodo `subMap` di prima per isolare i campi `id` e `interval` dalla map. Come prima, useremo l'operatore `map` per applicare il metodo `subMap` al primo elemento della tupla per ogni campione.

=== "Dopo"

    ```groovy title="main.nf" linenums="20" hl_lines="11-19"
        ch_combined_samples = ch_joined_samples
            .combine(ch_intervals)
            .map { grouping_key, normal, tumor, interval ->
                [
                    grouping_key + [interval: interval],
                    normal,
                    tumor
                ]
            }

        ch_grouped_samples = ch_combined_samples
            .map { grouping_key, normal, tumor ->
                [
                    grouping_key.subMap('id', 'interval'),
                    normal,
                    tumor
                ]
              }
              .view()
    ```

=== "Prima"

    ```groovy title="main.nf" linenums="20" hl_lines="10"
        ch_combined_samples = ch_joined_samples
            .combine(ch_intervals)
            .map { grouping_key, normal, tumor, interval ->
                [
                    grouping_key + [interval: interval],
                    normal,
                    tumor
                ]
            }
            .view()
    ```

Eseguiamolo di nuovo e verifichiamo il contenuto del canale:

```bash
nextflow run main.nf
```

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [hopeful_brenner] DSL2 - revision: 7f4f7fea76

    [[id:patientA, interval:chr1], patientA_rep1_normal.bam, patientA_rep1_tumor.bam]
    [[id:patientA, interval:chr2], patientA_rep1_normal.bam, patientA_rep1_tumor.bam]
    [[id:patientA, interval:chr3], patientA_rep1_normal.bam, patientA_rep1_tumor.bam]
    [[id:patientA, interval:chr1], patientA_rep2_normal.bam, patientA_rep2_tumor.bam]
    [[id:patientA, interval:chr2], patientA_rep2_normal.bam, patientA_rep2_tumor.bam]
    [[id:patientA, interval:chr3], patientA_rep2_normal.bam, patientA_rep2_tumor.bam]
    [[id:patientB, interval:chr1], patientB_rep1_normal.bam, patientB_rep1_tumor.bam]
    [[id:patientB, interval:chr2], patientB_rep1_normal.bam, patientB_rep1_tumor.bam]
    [[id:patientB, interval:chr3], patientB_rep1_normal.bam, patientB_rep1_tumor.bam]
    [[id:patientC, interval:chr1], patientC_rep1_normal.bam, patientC_rep1_tumor.bam]
    [[id:patientC, interval:chr2], patientC_rep1_normal.bam, patientC_rep1_tumor.bam]
    [[id:patientC, interval:chr3], patientC_rep1_normal.bam, patientC_rep1_tumor.bam]
    ```

Possiamo vedere che abbiamo isolato con successo i campi `id` e `interval`, ma non abbiamo ancora raggruppato i campioni.

!!! note "Nota"

    Stiamo scartando il campo `replicate` qui. Questo perché non ne abbiamo bisogno per l'elaborazione a valle. Dopo aver completato questo tutorial, provate a includerlo senza influenzare il raggruppamento successivo!

Raggruppiamo ora i campioni per questo nuovo elemento di raggruppamento, usando l'[operatore `groupTuple`](https://www.nextflow.io/docs/latest/operator.html#grouptuple).

=== "Dopo"

    ```groovy title="main.nf" linenums="30" hl_lines="9"
        ch_grouped_samples = ch_combined_samples
            .map { grouping_key, normal, tumor ->
                [
                    grouping_key.subMap('id', 'interval'),
                    normal,
                    tumor
                ]
              }
              .groupTuple()
              .view()
    ```

=== "Prima"

    ```groovy title="main.nf" linenums="30"
        ch_grouped_samples = ch_combined_samples
            .map { grouping_key, normal, tumor ->
                [
                    grouping_key.subMap('id', 'interval'),
                    normal,
                    tumor
                ]
              }
              .view()
    ```

È tutto! Abbiamo aggiunto solo una riga di codice. Vediamo cosa succede quando lo eseguiamo:

```bash
nextflow run main.nf
```

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [friendly_jang] DSL2 - revision: a1bee1c55d

    [[id:patientA, interval:chr1], [patientA_rep1_normal.bam, patientA_rep2_normal.bam], [patientA_rep1_tumor.bam, patientA_rep2_tumor.bam]]
    [[id:patientA, interval:chr2], [patientA_rep1_normal.bam, patientA_rep2_normal.bam], [patientA_rep1_tumor.bam, patientA_rep2_tumor.bam]]
    [[id:patientA, interval:chr3], [patientA_rep1_normal.bam, patientA_rep2_normal.bam], [patientA_rep1_tumor.bam, patientA_rep2_tumor.bam]]
    [[id:patientB, interval:chr1], [patientB_rep1_normal.bam], [patientB_rep1_tumor.bam]]
    [[id:patientB, interval:chr2], [patientB_rep1_normal.bam], [patientB_rep1_tumor.bam]]
    [[id:patientB, interval:chr3], [patientB_rep1_normal.bam], [patientB_rep1_tumor.bam]]
    [[id:patientC, interval:chr1], [patientC_rep1_normal.bam], [patientC_rep1_tumor.bam]]
    [[id:patientC, interval:chr2], [patientC_rep1_normal.bam], [patientC_rep1_tumor.bam]]
    [[id:patientC, interval:chr3], [patientC_rep1_normal.bam], [patientC_rep1_tumor.bam]]
    ```

Notate che la struttura dei dati è cambiata e all'interno di ogni elemento del canale i file sono ora contenuti in tuple come `[patientA_rep1_normal.bam, patientA_rep2_normal.bam]`. Questo perché quando usiamo `groupTuple`, Nextflow combina i singoli file per ogni campione di un gruppo. È importante ricordarlo quando si cerca di gestire i dati a valle.

!!! note "Nota"

    [`transpose`](https://www.nextflow.io/docs/latest/reference/operator.html#transpose) è l'opposto di groupTuple. Spacchetta gli elementi in un canale e li appiattisce. Provate ad aggiungere `transpose` e annullare il raggruppamento che abbiamo eseguito sopra!

### Takeaway

In questa sezione, avete imparato:

- **Raggruppare campioni correlati**: Come usare `groupTuple` per aggregare campioni per attributi comuni
- **Isolare le chiavi di raggruppamento**: Come usare `subMap` per estrarre campi specifici per il raggruppamento
- **Gestire strutture dati raggruppate**: Come lavorare con la struttura annidata creata da `groupTuple`
- **Gestione delle repliche tecniche**: Come raggruppare campioni che condividono le stesse condizioni sperimentali

---

## Riepilogo

In questa missione secondaria, avete imparato come suddividere e raggruppare i dati usando i canali.

Modificando i dati mentre scorrono attraverso la pipeline, potete costruire una pipeline scalabile senza usare cicli o istruzioni while, offrendo diversi vantaggi rispetto agli approcci più tradizionali:

- Possiamo scalare a quanti o pochi input vogliamo senza codice aggiuntivo
- Ci concentriamo sulla gestione del flusso di dati attraverso la pipeline, invece dell'iterazione
- Possiamo essere complessi o semplici quanto richiesto
- La pipeline diventa più dichiarativa, concentrandosi su cosa dovrebbe accadere piuttosto che su come dovrebbe accadere
- Nextflow ottimizzerà l'esecuzione per noi eseguendo operazioni indipendenti in parallelo

Padroneggiare queste operazioni sui canali vi permetterà di costruire pipeline flessibili e scalabili che gestiscono relazioni dati complesse senza ricorrere a cicli o programmazione iterativa, consentendo a Nextflow di ottimizzare l'esecuzione e parallelizzare automaticamente le operazioni indipendenti.

### Pattern chiave

1.  **Creare dati di input strutturati:** Partendo da un file CSV con meta map (basandosi sui pattern di [Metadata in workflows](../metadata/))

    ```groovy
    ch_samples = channel.fromPath("./data/samplesheet.csv")
        .splitCsv(header: true)
        .map{ row ->
          [[id:row.id, repeat:row.repeat, type:row.type], row.bam]
        }
    ```

2.  **Suddividere i dati in canali separati:** Abbiamo usato `filter` per dividere i dati in flussi indipendenti in base al campo `type`

    ```groovy
    channel.filter { it.type == 'tumor' }
    ```

3.  **Unire campioni corrispondenti:** Abbiamo usato `join` per ricombinare campioni correlati in base ai campi `id` e `repeat`

    - Unire due canali per chiave (primo elemento della tupla)

    ```groovy
    tumor_ch.join(normal_ch)
    ```

    - Estrarre la chiave di unione e unire per questo valore

    ```groovy
    tumor_ch.map { meta, file -> [meta.id, meta, file] }
        .join(
          normal_ch.map { meta, file -> [meta.id, meta, file] }
        )
    ```

    - Unire su più campi usando subMap

    ```groovy
    tumor_ch.map { meta, file -> [meta.subMap(['id', 'repeat']), meta, file] }
        .join(
          normal_ch.map { meta, file -> [meta.subMap(['id', 'repeat']), meta, file] }
        )
    ```

4.  **Distribuire sugli intervalli:** Abbiamo usato `combine` per creare prodotti cartesiani di campioni con intervalli genomici per l'elaborazione parallela.

    ```groovy
    samples_ch.combine(intervals_ch)
    ```

5.  **Aggregare per chiavi di raggruppamento:** Abbiamo usato `groupTuple` per raggruppare per il primo elemento in ogni tupla, raccogliendo così i campioni che condividono i campi `id` e `interval` e unendo le repliche tecniche.

    ```groovy
    channel.groupTuple()
    ```

6.  **Ottimizzare la struttura dati:** Abbiamo usato `subMap` per estrarre campi specifici e creato una closure con nome per rendere le trasformazioni riutilizzabili.

    - Estrarre campi specifici da una map

    ```groovy
    meta.subMap(['id', 'repeat'])
    ```

    - Usare una closure con nome per trasformazioni riutilizzabili

    ```groovy
    getSampleIdAndReplicate = { meta, file -> [meta.subMap(['id', 'repeat']), file] }
    channel.map(getSampleIdAndReplicate)
    ```

### Risorse aggiuntive

- [filter](https://www.nextflow.io/docs/latest/operator.html#filter)
- [map](https://www.nextflow.io/docs/latest/operator.html#map)
- [join](https://www.nextflow.io/docs/latest/operator.html#join)
- [groupTuple](https://www.nextflow.io/docs/latest/operator.html#grouptuple)
- [combine](https://www.nextflow.io/docs/latest/operator.html#combine)

---

## Cosa c'è dopo?

Tornate al [menu delle missioni secondarie](../) o cliccate il pulsante in basso a destra della pagina per passare all'argomento successivo nell'elenco.
