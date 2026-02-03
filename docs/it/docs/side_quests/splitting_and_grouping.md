# Divisione e Raggruppamento

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduzione assistita da IA - [scopri di più e suggerisci miglioramenti](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Nextflow fornisce strumenti potenti per lavorare con i dati in modo flessibile. Una capacità chiave è la divisione dei dati in flussi diversi e il successivo raggruppamento degli elementi correlati. Questo è particolarmente prezioso nei workflow di bioinformatica dove è necessario processare diversi tipi di campioni separatamente prima di combinare i risultati per l'analisi.

Pensatelo come smistare la posta: separate le lettere per destinazione, processate ogni pila in modo diverso, poi ricombinate gli elementi destinati alla stessa persona. Nextflow utilizza operatori speciali per realizzare questo con dati scientifici. Questo approccio è anche comunemente noto come pattern **scatter/gather** nel computing distribuito e nei workflow di bioinformatica.

Il sistema di canali di Nextflow è al centro di questa flessibilità. I canali connettono diverse parti del workflow, permettendo ai dati di fluire attraverso l'analisi. È possibile creare canali multipli da una singola fonte di dati, processare ciascun canale in modo diverso, e poi unire i canali quando necessario. Questo approccio permette di progettare workflow che rispecchiano naturalmente i percorsi ramificati e convergenti delle analisi bioinformatiche complesse.

### Obiettivi di apprendimento

In questa missione secondaria, imparerete a dividere e raggruppare dati usando gli operatori di canale di Nextflow.
Inizieremo con un file CSV contenente informazioni sui campioni e file di dati associati, poi manipoleremo e riorganizzeremo questi dati.

Alla fine di questa missione secondaria, sarete in grado di separare e combinare flussi di dati in modo efficace, utilizzando le seguenti tecniche:

- Leggere dati da file usando `splitCsv`
- Filtrare e trasformare dati con `filter` e `map`
- Combinare dati correlati usando `join` e `groupTuple`
- Creare combinazioni di dati con `combine` per l'elaborazione parallela
- Ottimizzare la struttura dati usando `subMap` e strategie di deduplicazione
- Costruire funzioni riutilizzabili con closure nominate per aiutarvi a manipolare le strutture dei canali

Queste competenze vi aiuteranno a costruire workflow che possono gestire file di input multipli e diversi tipi di dati in modo efficiente, mantenendo una struttura di codice pulita e manutenibile.

### Prerequisiti

Prima di intraprendere questa missione secondaria, dovreste:

- Aver completato il tutorial [Hello Nextflow](../hello_nextflow/README.md) o un corso equivalente per principianti.
- Essere a proprio agio con i concetti e meccanismi base di Nextflow (processi, canali, operatori, lavoro con file, meta data)

**Opzionale:** Raccomandiamo di completare prima la missione secondaria [Metadata nei workflow](./metadata.md).
Questa copre i fondamenti della lettura di file CSV con `splitCsv` e della creazione di meta map, che useremo intensamente qui.

---

## 0. Iniziare

#### Aprire il codespace di formazione

Se non lo avete ancora fatto, assicuratevi di aprire l'ambiente di formazione come descritto nella [Configurazione dell'Ambiente](../envsetup/index.md).

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Spostarsi nella directory del progetto

Spostiamoci nella directory dove si trovano i file per questo tutorial.

```bash
cd side-quests/splitting_and_grouping
```

Potete impostare VSCode per concentrarsi su questa directory:

```bash
code .
```

#### Rivedere i materiali

Troverete un file workflow principale e una directory `data` contenente un samplesheet chiamato `samplesheet.csv`.

```console title="Contenuti della directory"
.
├── data
│   └── samplesheet.csv
└── main.nf
```

Il samplesheet contiene informazioni sui campioni di diversi pazienti, inclusi l'ID paziente, il numero di ripetizione del campione, il tipo (normale o tumorale), e i percorsi ai file di dati ipotetici (che in realtà non esistono, ma faremo finta che esistano).

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

Questo samplesheet elenca otto campioni da tre pazienti (A, B, C).

Per ciascun paziente, abbiamo campioni di tipo `tumor` (tipicamente originati da biopsie tumorali) o `normal` (prelevati da tessuto sano o sangue).
Se non avete familiarità con l'analisi del cancro, sappiate semplicemente che questo corrisponde a un modello sperimentale che utilizza campioni appaiati tumore/normale per eseguire analisi contrastive.

Per il paziente A specificatamente, abbiamo due set di repliche tecniche.

!!! note

    Non preoccupatevi se non avete familiarità con questo disegno sperimentale, non è critico per comprendere questo tutorial.

#### Rivedere l'incarico

La vostra sfida è scrivere un workflow Nextflow che:

1. **Legga** i dati dei campioni da un file CSV e li strutturi con meta map
2. **Separi** i campioni in canali diversi in base al tipo (normale vs tumore)
3. **Unisca** le coppie tumore/normale corrispondenti per ID paziente e numero di replica
4. **Distribuisca** i campioni attraverso intervalli genomici per l'elaborazione parallela
5. **Raggruppi** i campioni correlati per l'analisi downstream

Questo rappresenta un pattern bioinformatico comune dove è necessario dividere i dati per l'elaborazione indipendente, poi ricombinare gli elementi correlati per l'analisi comparativa.

#### Checklist di preparazione

Pensate di essere pronti per iniziare?

- [ ] Comprendo l'obiettivo di questo corso e i suoi prerequisiti
- [ ] Il mio codespace è attivo e funzionante
- [ ] Ho impostato la directory di lavoro appropriatamente
- [ ] Comprendo l'incarico

Se potete spuntare tutte le caselle, siete pronti per partire.

---

## 1. Leggere i dati dei campioni

### 1.1. Leggere i dati dei campioni con `splitCsv` e creare meta map

Iniziamo leggendo i dati dei campioni con `splitCsv` e organizzandoli nel pattern meta map. Nel `main.nf`, vedrete che abbiamo già iniziato il workflow.

```groovy title="main.nf" linenums="1" hl_lines="2"
workflow {
    ch_samplesheet = channel.fromPath("./data/samplesheet.csv")
}
```

!!! note

    In tutto questo tutorial, useremo il prefisso `ch_` per tutte le variabili canale per indicare chiaramente che sono canali Nextflow.

Se avete completato la missione secondaria [Metadata nei workflow](./metadata.md), riconoscerete questo pattern. Useremo `splitCsv` per leggere il CSV e strutturare immediatamente i dati con una meta map per separare i metadata dai percorsi dei file.

!!! info

    Incontreremo due concetti diversi chiamati `map` in questa formazione:

    - **Struttura dati**: La map Groovy (equivalente a dizionari/hash in altri linguaggi) che memorizza coppie chiave-valore
    - **Operatore di canale**: L'operatore `.map()` che trasforma elementi in un canale

    Chiariremo quale intendiamo nel contesto, ma questa distinzione è importante da comprendere quando si lavora con Nextflow.

Applicate questi cambiamenti al `main.nf`:

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

Questo combina l'operazione `splitCsv` (lettura del CSV con intestazioni) e l'operazione `map` (strutturazione dei dati come tuple `[meta, file]`) in un solo passaggio. Applicate quella modifica ed eseguite il pipeline:

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

Ora abbiamo un canale dove ogni elemento è una tupla `[meta, file]` - metadata separati dai percorsi dei file. Questa struttura ci permette di dividere e raggruppare il nostro carico di lavoro in base ai campi dei metadata.

---

## 2. Filtrare e trasformare i dati

### 2.1. Filtrare i dati con `filter`

Possiamo usare l'[operatore `filter`](https://www.nextflow.io/docs/latest/operator.html#filter) per filtrare i dati in base a una condizione. Diciamo che vogliamo processare solo campioni normali. Possiamo farlo filtrando i dati in base al campo `type`. Inseriamo questo prima dell'operatore `view`.

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

Eseguite di nuovo il workflow per vedere il risultato filtrato:

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

L'operatore `filter` prende una closure che viene applicata a ciascun elemento nel canale. Se la closure restituisce `true`, l'elemento viene incluso; se restituisce `false`, l'elemento viene escluso.

Nel nostro caso, vogliamo mantenere solo i campioni dove `meta.type == 'normal'`. La closure usa la tupla `meta,file` per riferirsi a ciascun campione, accede al tipo di campione con `meta.type`, e verifica se è uguale a `'normal'`.

Questo viene realizzato con la singola closure che abbiamo introdotto sopra:

```groovy title="main.nf" linenums="7"
    .filter { meta, file -> meta.type == 'normal' }
```

### 2.2. Creare canali filtrati separati

Attualmente stiamo applicando il filtro al canale creato direttamente dal CSV, ma vogliamo filtrare questo in più modi di uno, quindi riscriviamo la logica per creare un canale filtrato separato per i campioni normali:

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

Eseguite il pipeline per vedere i risultati:

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

### Sintesi

In questa sezione, avete imparato:

- **Filtrare i dati**: Come filtrare i dati con `filter`
- **Dividere i dati**: Come dividere i dati in canali diversi in base a una condizione
- **Visualizzare i dati**: Come usare `view` per stampare i dati ed etichettare l'output da canali diversi

Abbiamo ora separato i campioni normali e tumorali in due canali diversi. Successivamente, uniremo i campioni normali e tumorali sul campo `id`.

---

## 3. Unire i canali per identificatori

Nella sezione precedente, abbiamo separato i campioni normali e tumorali in due canali diversi. Questi potrebbero essere processati indipendentemente usando processi o workflow specifici in base al loro tipo. Ma cosa succede quando vogliamo confrontare i campioni normali e tumorali dello stesso paziente? A questo punto, dobbiamo unirli assicurandoci di abbinare i campioni in base al loro campo `id`.

Nextflow include molti metodi per combinare i canali, ma in questo caso l'operatore più appropriato è [`join`](https://www.nextflow.io/docs/latest/operator.html#join). Se avete familiarità con SQL, agisce come l'operazione `JOIN`, dove specifichiamo la chiave su cui unire e il tipo di join da eseguire.

### 3.1. Usare `map` e `join` per combinare in base all'ID paziente

Se controlliamo la documentazione di [`join`](https://www.nextflow.io/docs/latest/operator.html#join), possiamo vedere che per impostazione predefinita unisce due canali in base al primo elemento in ciascuna tupla.

#### 3.1.1. Controllare la struttura dati

Se non avete ancora l'output della console disponibile, eseguiamo il pipeline per controllare la nostra struttura dati e vedere come dobbiamo modificarla per unire sul campo `id`.

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

Possiamo vedere che il campo `id` è il primo elemento in ciascuna meta map. Perché `join` funzioni, dovremmo isolare il campo `id` in ciascuna tupla. Dopo di che, possiamo semplicemente usare l'operatore `join` per combinare i due canali.

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

Potrebbe essere sottile, ma dovreste essere in grado di vedere che il primo elemento in ciascuna tupla è il campo `id`.

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

È un po' difficile da leggere perché è molto largo, ma dovreste essere in grado di vedere che i campioni sono stati uniti per il campo `id`. Ciascuna tupla ora ha il formato:

- `id`: L'ID del campione
- `normal_meta_map`: I metadata del campione normale inclusi tipo, replica e percorso al file bam
- `normal_sample_file`: Il file del campione normale
- `tumor_meta_map`: I metadata del campione tumorale inclusi tipo, replica e percorso al file bam
- `tumor_sample`: Il campione tumorale inclusi tipo, replica e percorso al file bam

!!! warning

    L'operatore `join` scarterà qualsiasi tupla non abbinata. In questo esempio, ci siamo assicurati che tutti i campioni fossero abbinati per tumore e normale ma se questo non è vero dovete usare il parametro `remainder: true` per mantenere le tuple non abbinate. Controllate la [documentazione](https://www.nextflow.io/docs/latest/operator.html#join) per maggiori dettagli.

Quindi ora sapete come usare `map` per isolare un campo in una tupla, e come usare `join` per combinare tuple in base al primo campo.
Con questa conoscenza, possiamo combinare con successo i canali in base a un campo condiviso.

Successivamente, considereremo la situazione in cui volete unire su campi multipli.

### 3.2. Unire su campi multipli

Abbiamo 2 repliche per il sampleA, ma solo 1 per sampleB e sampleC. In questo caso siamo stati in grado di unirli efficacemente usando il campo `id`, ma cosa succederebbe se fossero fuori sincronizzazione? Potremmo mescolare i campioni normali e tumorali da diverse repliche!

Per evitare questo, possiamo unire su campi multipli. Ci sono in realtà diversi modi per ottenere questo ma ci concentreremo sulla creazione di una nuova chiave di unione che include sia il campo `id` che il numero di `replicate`.

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

Ora dovremmo vedere l'unione avvenire ma usando sia i campi `id` che `repeat`. Eseguite il workflow:

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

Notate come abbiamo una tupla di due elementi (campi `id` e `repeat`) come primo elemento di ciascun risultato unito. Questo dimostra come elementi complessi possono essere usati come chiave di unione, abilitando abbinamenti abbastanza intricati tra campioni delle stesse condizioni.

Se volete esplorare altri modi per unire su chiavi diverse, controllate la [documentazione dell'operatore join](https://www.nextflow.io/docs/latest/operator.html#join) per opzioni ed esempi aggiuntivi.

### 3.3. Usare `subMap` per creare una nuova chiave di unione

L'approccio precedente perde i nomi dei campi dalla nostra chiave di unione - i campi `id` e `repeat` diventano solo una lista di valori. Per mantenere i nomi dei campi per l'accesso successivo, possiamo usare il [metodo `subMap`](<https://docs.groovy-lang.org/latest/html/groovy-jdk/java/util/Map.html#subMap(java.util.Collection)>).

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

Ora abbiamo una nuova chiave di unione che non solo include i campi `id` e `repeat` ma mantiene anche i nomi dei campi così possiamo accedervi successivamente per nome, ad esempio `meta.id` e `meta.repeat`.

### 3.4. Usare una closure nominata in map

Per evitare duplicazioni e ridurre errori, possiamo usare una closure nominata. Una closure nominata ci permette di creare una funzione riutilizzabile che possiamo chiamare in più posti.

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

Abbiamo definito la trasformazione map come variabile nominata che possiamo riutilizzare.

Notate che convertiamo anche il percorso del file in un oggetto Path usando `file()` così che qualsiasi processo che riceve questo canale possa gestire il file correttamente (per maggiori informazioni vedere [Lavorare con i file](./working_with_files.md)).

Implementiamo la closure nel nostro workflow:

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

!!! note

    L'operatore `map` è passato dall'usare `{ }` all'usare `( )` per passare la closure come argomento. Questo perché l'operatore `map` si aspetta una closure come argomento e `{ }` viene usato per definire una closure anonima. Quando si chiama una closure nominata, usare la sintassi `( )`.

Eseguite di nuovo il workflow per verificare che tutto funzioni ancora:

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

Usare una closure nominata ci permette di riutilizzare la stessa trasformazione in più posti, riducendo il rischio di errori e rendendo il codice più leggibile e manutenibile.

### 3.5. Ridurre la duplicazione dei dati

Abbiamo molti dati duplicati nel nostro workflow. Ciascun elemento nei campioni uniti ripete i campi `id` e `repeat`. Dato che questa informazione è già disponibile nella chiave di raggruppamento, possiamo evitare questa ridondanza. Come promemoria, la nostra struttura dati attuale appare così:

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

Dato che i campi `id` e `repeat` sono disponibili nella chiave di raggruppamento, rimuoviamoli dal resto di ciascun elemento del canale per evitare la duplicazione. Possiamo farlo usando il metodo `subMap` per creare una nuova map con solo il campo `type`. Questo approccio ci permette di mantenere tutte le informazioni necessarie eliminando la ridondanza nella nostra struttura dati.

=== "Dopo"

    ```groovy title="main.nf" linenums="8" hl_lines="1"
        getSampleIdAndReplicate = { meta, bam -> [ meta.subMap(['id', 'repeat']), meta.subMap(['type']), file(bam) ] }
    ```

=== "Prima"

    ```groovy title="main.nf" linenums="8" hl_lines="1"
        getSampleIdAndReplicate = { meta, bam -> [ meta.subMap(['id', 'repeat']), meta, file(bam) ] }
    ```

Ora la closure restituisce una tupla dove il primo elemento contiene i campi `id` e `repeat`, e il secondo elemento contiene solo il campo `type`. Questo elimina la ridondanza memorizzando le informazioni `id` e `repeat` una volta nella chiave di raggruppamento, mantenendo tutte le informazioni necessarie.

Eseguite il workflow per vedere come appare:

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

Possiamo vedere che dichiariamo i campi `id` e `repeat` solo una volta nella chiave di raggruppamento e abbiamo il campo `type` nei dati del campione. Non abbiamo perso alcuna informazione ma siamo riusciti a rendere i contenuti del nostro canale più succinti.

### 3.6. Rimuovere informazioni ridondanti

Abbiamo rimosso le informazioni duplicate sopra, ma abbiamo ancora altre informazioni ridondanti nei nostri canali.

All'inizio, abbiamo separato i campioni normali e tumorali usando `filter`, poi li abbiamo uniti in base alle chiavi `id` e `repeat`. L'operatore `join` preserva l'ordine in cui le tuple vengono unite, quindi nel nostro caso, con i campioni normali sul lato sinistro e i campioni tumorali sulla destra, il canale risultante mantiene questa struttura: `id, <elementi normali>, <elementi tumorali>`.

Dato che conosciamo la posizione di ciascun elemento nel nostro canale, possiamo semplificare ulteriormente la struttura eliminando i metadata `[type:normal]` e `[type:tumor]`.

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

### Sintesi

In questa sezione, avete imparato:

- **Manipolare le Tuple**: Come usare `map` per isolare un campo in una tupla
- **Unire le Tuple**: Come usare `join` per combinare tuple in base al primo campo
- **Creare Chiavi di Unione**: Come usare `subMap` per creare una nuova chiave di unione
- **Closure Nominate**: Come usare una closure nominata in map
- **Unione su Campi Multipli**: Come unire su campi multipli per abbinamenti più precisi
- **Ottimizzazione della Struttura Dati**: Come semplificare la struttura del canale rimuovendo informazioni ridondanti

Ora avete un workflow che può dividere un samplesheet, filtrare i campioni normali e tumorali, unirli per ID campione e numero di replica, poi stampare i risultati.

Questo è un pattern comune nei workflow di bioinformatica dove è necessario abbinare campioni o altri tipi di dati dopo l'elaborazione indipendente, quindi è una competenza utile. Successivamente, vedremo come ripetere un campione più volte.

## 4. Distribuire campioni su intervalli

Un pattern chiave nei workflow di bioinformatica è la distribuzione dell'analisi attraverso regioni genomiche. Ad esempio, il variant calling può essere parallelizzato dividendo il genoma in intervalli (come cromosomi o regioni più piccole). Questa strategia di parallelizzazione migliora significativamente l'efficienza del pipeline distribuendo il carico computazionale attraverso core o nodi multipli, riducendo il tempo complessivo di esecuzione.

Nella sezione seguente, dimostreremo come distribuire i nostri dati dei campioni attraverso intervalli genomici multipli. Abbineremo ciascun campione con ogni intervallo, permettendo l'elaborazione parallela di diverse regioni genomiche. Questo moltiplicherà la dimensione del nostro dataset per il numero di intervalli, creando unità di analisi indipendenti multiple che possono essere riportate insieme successivamente.

### 4.1. Distribuire campioni su intervalli usando `combine`

Iniziamo creando un canale di intervalli. Per mantenere la vita semplice, useremo solo 3 intervalli che definiremo manualmente. In un workflow reale, potreste leggerli da un file di input o persino creare un canale con molti file di intervalli.

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

Ora ricordate, vogliamo ripetere ciascun campione per ogni intervallo. Questo è a volte chiamato prodotto Cartesiano dei campioni e intervalli. Possiamo ottenere questo usando l'[operatore `combine`](https://www.nextflow.io/docs/latest/operator.html#combine). Questo prenderà ogni elemento dal canale 1 e lo ripeterà per ogni elemento nel canale 2. Aggiungiamo un operatore combine al nostro workflow:

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

Ora eseguiamolo e vediamo cosa succede:

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

Successo! Abbiamo ripetuto ogni campione per ogni singolo intervallo nella nostra lista di 3 intervalli. Abbiamo effettivamente triplicato il numero di elementi nel nostro canale.

È un po' difficile da leggere però, quindi nella prossima sezione lo riordineremo.

### 4.2. Organizzare il canale

Possiamo usare l'operatore `map` per riordinare e ristrutturare i dati del nostro campione così che siano più facili da comprendere. Spostiamo la stringa degli intervalli alla map di unione al primo elemento.

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

Analizziamo cosa fa questa operazione map passo dopo passo.

Prima, usiamo parametri nominati per rendere il codice più leggibile. Usando i nomi `grouping_key`, `normal`, `tumor` e `interval`, possiamo riferirci agli elementi nella tupla per nome invece che per indice:

```groovy
        .map { grouping_key, normal, tumor, interval ->
```

Successivamente, combiniamo la `grouping_key` con il campo `interval`. La `grouping_key` è una map contenente i campi `id` e `repeat`. Creiamo una nuova map con l'`interval` e le uniamo usando l'addizione di map Groovy (`+`):

```groovy
                grouping_key + [interval: interval],
```

Infine, restituiamo questo come tupla con tre elementi: la map di metadata combinata, il file del campione normale, e il file del campione tumorale:

```groovy
            [
                grouping_key + [interval: interval],
                normal,
                tumor
            ]
```

Eseguiamolo di nuovo e controlliamo i contenuti del canale:

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

Usare `map` per coercere i dati nella struttura corretta può essere complicato, ma è cruciale per la manipolazione efficace dei dati.

Ora abbiamo ogni campione ripetuto attraverso tutti gli intervalli genomici, creando unità di analisi indipendenti multiple che possono essere processate in parallelo. Ma cosa succede se vogliamo riportare insieme campioni correlati? Nella prossima sezione, impareremo come raggruppare campioni che condividono attributi comuni.

### Sintesi

In questa sezione, avete imparato:

- **Distribuire campioni su intervalli**: Come usare `combine` per ripetere campioni su intervalli
- **Creare prodotti Cartesiani**: Come generare tutte le combinazioni di campioni e intervalli
- **Organizzare la struttura del canale**: Come usare `map` per ristrutturare i dati per una migliore leggibilità
- **Preparazione all'elaborazione parallela**: Come impostare i dati per l'analisi distribuita

## 5. Aggregare campioni usando `groupTuple`

Nelle sezioni precedenti, abbiamo imparato come dividere i dati da un file di input e filtrare per campi specifici (nel nostro caso campioni normali e tumorali). Ma questo copre solo un singolo tipo di unione. Cosa succede se vogliamo raggruppare campioni per un attributo specifico? Ad esempio, invece di unire coppie tumore-normale abbinate, potremmo voler processare tutti i campioni da "sampleA" insieme indipendentemente dal loro tipo. Questo pattern è comune nei workflow di bioinformatica dove potreste voler processare campioni correlati separatamente per ragioni di efficienza prima di confrontare o combinare i risultati alla fine.

Nextflow include metodi integrati per fare questo, il principale che vedremo è `groupTuple`.

Iniziamo raggruppando tutti i nostri campioni che hanno gli stessi campi `id` e `interval`, questo sarebbe tipico di un'analisi dove vogliamo raggruppare repliche tecniche ma mantenere campioni significativamente diversi separati.

Per farlo, dovremmo separare le nostre variabili di raggruppamento così possiamo usarle in isolamento.

Il primo passo è simile a quello che abbiamo fatto nella sezione precedente. Dobbiamo isolare la nostra variabile di raggruppamento come primo elemento della tupla. Ricordate, il nostro primo elemento è attualmente una map dei campi `id`, `repeat` e `interval`:

```groovy title="main.nf" linenums="1"
{
  "id": "sampleA",
  "repeat": "1",
  "interval": "chr1"
}
```

Possiamo riutilizzare il metodo `subMap` da prima per isolare i nostri campi `id` e `interval` dalla map. Come prima, useremo l'operatore `map` per applicare il metodo `subMap` al primo elemento della tupla per ciascun campione.

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

Eseguiamolo di nuovo e controlliamo i contenuti del canale:

```bash
nextflow run main.nf
```

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [hop
