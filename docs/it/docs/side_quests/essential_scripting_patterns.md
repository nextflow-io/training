# Pattern di Scripting Essenziali in Nextflow

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduzione assistita da IA - [scopri di più e suggerisci miglioramenti](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Nextflow è un linguaggio di programmazione che viene eseguito sulla Java Virtual Machine. Sebbene Nextflow sia costruito su [Groovy](http://groovy-lang.org/) e condivida gran parte della sua sintassi, Nextflow è più di un semplice "Groovy con estensioni" -- è un linguaggio autonomo con una [sintassi](https://nextflow.io/docs/latest/reference/syntax.html) e una [libreria standard](https://nextflow.io/docs/latest/reference/stdlib.html) completamente specificate.

Potete scrivere molto codice Nextflow senza avventurarvi oltre la sintassi di base per variabili, mappe e liste. La maggior parte dei tutorial su Nextflow si concentra sull'orchestrazione del flusso di lavoro (canali, processi e flusso di dati), e potete arrivare sorprendentemente lontano con solo questo.

Tuttavia, quando avete bisogno di manipolare dati, analizzare nomi di file complessi, implementare logica condizionale o costruire flussi di lavoro robusti per la produzione, è utile pensare a due aspetti distinti del vostro codice: **dataflow** (canali, operatori, processi e workflow) e **scripting** (il codice all'interno di closure, funzioni e script di processo). Sebbene questa distinzione sia in qualche modo arbitraria—è tutto codice Nextflow—fornisce un modello mentale utile per capire quando state orchestrando la vostra pipeline rispetto a quando state manipolando dati. Padroneggiare entrambi migliora notevolmente la vostra capacità di scrivere flussi di lavoro chiari e manutenibili.

### Obiettivi di apprendimento

Questa side quest vi accompagna in un viaggio pratico dai concetti di base ai pattern pronti per la produzione.
Trasformeremo un semplice flusso di lavoro che legge CSV in una sofisticata pipeline bioinformatica, evolvendola passo dopo passo attraverso sfide realistiche:

- **Comprendere i confini:** Distinguere tra operazioni di dataflow e scripting, e capire come lavorano insieme
- **Manipolazione dei dati:** Estrarre, trasformare e creare sottoinsiemi di mappe e collezioni usando operatori potenti
- **Elaborazione di stringhe:** Analizzare schemi complessi di denominazione dei file con pattern regex e padroneggiare l'interpolazione di variabili
- **Funzioni riutilizzabili:** Estrarre logica complessa in funzioni con nome per flussi di lavoro più puliti e manutenibili
- **Logica dinamica:** Costruire processi che si adattano a diversi tipi di input e usare closure per l'allocazione dinamica delle risorse
- **Routing condizionale:** Instradare intelligentemente i campioni attraverso diversi processi in base alle loro caratteristiche di metadati
- **Operazioni sicure:** Gestire con eleganza i dati mancanti con operatori null-safe e validare gli input con messaggi di errore chiari
- **Handler basati su configurazione:** Usare handler di eventi del flusso di lavoro per logging, notifiche e gestione del ciclo di vita

### Prerequisiti

Prima di intraprendere questa side quest, dovreste:

- Aver completato il tutorial [Hello Nextflow](../hello_nextflow/README.md) o un corso equivalente per principianti.
- Essere a vostro agio nell'uso di concetti e meccanismi di base di Nextflow (processi, canali, operatori, lavorare con file, metadati)
- Avere una familiarità di base con costrutti di programmazione comuni (variabili, mappe, liste)

Questo tutorial spiegherà i concetti di programmazione man mano che li incontriamo, quindi non avete bisogno di un'esperienza di programmazione estesa.
Inizieremo con concetti fondamentali e costruiremo fino a pattern avanzati.

---

## 0. Iniziamo

#### Aprite il codespace di formazione

Se non l'avete ancora fatto, assicuratevi di aprire l'ambiente di formazione come descritto in [Configurazione dell'Ambiente](../envsetup/index.md).

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Spostatevi nella directory del progetto

Spostiamoci nella directory dove si trovano i file per questo tutorial.

```bash
cd side-quests/essential_scripting_patterns
```

#### Esaminate i materiali

Troverete un file di flusso di lavoro principale e una directory `data` contenente file di dati di esempio.

```console title="Directory contents"
.
├── collect.nf
├── data
│   ├── samples.csv
│   └── sequences
│       ├── SAMPLE_001_S1_L001_R1_001.fastq
│       ├── SAMPLE_002_S2_L001_R1_001.fastq
│       └── SAMPLE_003_S3_L001_R1_001.fastq
├── main.nf
├── modules
│   ├── fastp.nf
│   ├── generate_report.nf
│   └── trimgalore.nf
└── nextflow.config
```

Il nostro CSV di esempio contiene informazioni su campioni biologici che necessitano di elaborazioni diverse in base alle loro caratteristiche:

```console title="samples.csv"
sample_id,organism,tissue_type,sequencing_depth,file_path,quality_score
SAMPLE_001,human,liver,30000000,data/sequences/SAMPLE_001_S1_L001_R1_001.fastq,38.5
SAMPLE_002,mouse,brain,25000000,data/sequences/SAMPLE_002_S2_L001_R1_001.fastq,35.2
SAMPLE_003,human,kidney,45000000,data/sequences/SAMPLE_003_S3_L001_R1_001.fastq,42.1
```

Useremo questo dataset realistico per esplorare tecniche di programmazione pratiche che incontrerete nei flussi di lavoro bioinformatici reali.

#### Checklist di preparazione

Pensate di essere pronti per iniziare?

- [ ] Comprendo l'obiettivo di questo corso e i suoi prerequisiti
- [ ] Il mio codespace è attivo e funzionante
- [ ] Ho impostato la mia directory di lavoro in modo appropriato

Se potete spuntare tutte le caselle, siete pronti per partire.

---

## 1. Dataflow vs Scripting: Comprendere i Confini

### 1.1. Identificare Cosa è Cosa

Quando scrivete flussi di lavoro Nextflow, è importante distinguere tra **dataflow** (come i dati si muovono attraverso canali e processi) e **scripting** (il codice che manipola i dati e prende decisioni). Costruiamo un flusso di lavoro che dimostri come lavorano insieme.

#### 1.1.1. Flusso di Lavoro Nextflow di Base

Iniziamo con un semplice flusso di lavoro che legge solo il file CSV (l'abbiamo già fatto per voi in `main.nf`):

```groovy title="main.nf" linenums="1"
workflow {
    ch_samples = channel.fromPath("./data/samples.csv")
        .splitCsv(header: true)
        .view()
}
```

Il blocco `workflow` definisce la struttura della nostra pipeline, mentre `channel.fromPath()` crea un canale da un percorso di file. L'operatore `.splitCsv()` elabora il file CSV e converte ogni riga in una struttura dati di tipo mappa.

Eseguite questo flusso di lavoro per vedere i dati CSV grezzi:

```bash
nextflow run main.nf
```

??? success "Output del comando"

    ```console
    Launching `main.nf` [marvelous_tuckerman] DSL2 - revision: 6113e05c17

    [sample_id:SAMPLE_001, organism:human, tissue_type:liver, sequencing_depth:30000000, file_path:data/sequences/SAMPLE_001_S1_L001_R1_001.fastq, quality_score:38.5]
    [sample_id:SAMPLE_002, organism:mouse, tissue_type:brain, sequencing_depth:25000000, file_path:data/sequences/SAMPLE_002_S2_L001_R1_001.fastq, quality_score:35.2]
    [sample_id:SAMPLE_003, organism:human, tissue_type:kidney, sequencing_depth:45000000, file_path:data/sequences/SAMPLE_003_S3_L001_R1_001.fastq, quality_score:42.1]
    ```

#### 1.1.2. Aggiungere l'Operatore Map

Ora aggiungeremo scripting per trasformare i dati, usando l'operatore `.map()` che probabilmente conoscete già. Questo operatore prende una 'closure' dove possiamo scrivere codice per trasformare ogni elemento.

!!! note "Nota"

    Una **closure** è un blocco di codice che può essere passato in giro ed eseguito successivamente. Pensatela come una funzione che definite inline. Le closure sono scritte con parentesi graffe `{ }` e possono prendere parametri. Sono fondamentali per il funzionamento degli operatori Nextflow e se avete scritto Nextflow per un po', potreste averle già usate senza rendervene conto!

Ecco come appare quell'operazione map:

=== "Dopo"

    ```groovy title="main.nf" linenums="2" hl_lines="3-6"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                return row
            }
            .view()
    ```

=== "Prima"

    ```groovy title="main.nf" linenums="2" hl_lines="3"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .view()
    ```

Questa è la nostra prima **closure** - una funzione anonima che potete passare come argomento (simile alle lambda in Python o alle arrow function in JavaScript). Le closure sono essenziali per lavorare con gli operatori Nextflow.

La closure `{ row -> return row }` prende un parametro `row` (potrebbe essere qualsiasi nome: `item`, `sample`, ecc.).

Quando l'operatore `.map()` elabora ogni elemento del canale, passa quell'elemento alla vostra closure. Qui, `row` contiene una riga CSV alla volta.

Applicate questa modifica ed eseguite il flusso di lavoro:

```bash
nextflow run main.nf
```

Vedrete lo stesso output di prima, perché stiamo semplicemente restituendo l'input invariato. Questo conferma che l'operatore map funziona correttamente. Ora iniziamo a trasformare i dati.

#### 1.1.3. Creare una Struttura Dati Map

Ora scriveremo logica di **scripting** all'interno della nostra closure per trasformare ogni riga di dati. Qui è dove elaboriamo singoli elementi di dati piuttosto che orchestrare il flusso di dati.

=== "Dopo"

    ```groovy title="main.nf" linenums="2" hl_lines="4-12"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                // Scripting per la trasformazione dei dati
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                return sample_meta
            }
            .view()
    ```

=== "Prima"

    ```groovy title="main.nf" linenums="2" hl_lines="4"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                return row
            }
            .view()
    ```

La mappa `sample_meta` è una struttura dati chiave-valore (come i dizionari in Python, gli oggetti in JavaScript o gli hash in Ruby) che memorizza informazioni correlate: ID del campione, organismo, tipo di tessuto, profondità di sequenziamento e punteggio di qualità.

Usiamo metodi di manipolazione delle stringhe come `.toLowerCase()` e `.replaceAll()` per pulire i nostri dati, e metodi di conversione di tipo come `.toInteger()` e `.toDouble()` per convertire i dati stringa dal CSV nei tipi numerici appropriati.

Applicate questa modifica ed eseguite il flusso di lavoro:

```bash
nextflow run main.nf
```

??? success "Output del comando"

    ```console
    [id:sample_001, organism:human, tissue:liver, depth:30000000, quality:38.5]
    [id:sample_002, organism:mouse, tissue:brain, depth:25000000, quality:35.2]
    [id:sample_003, organism:human, tissue:kidney, depth:45000000, quality:42.1]
    ```

#### 1.1.4. Aggiungere Logica Condizionale

Ora aggiungiamo più scripting - questa volta usando un operatore ternario per prendere decisioni basate sui valori dei dati.

Fate la seguente modifica:

=== "Dopo"

    ```groovy title="main.nf" linenums="2" hl_lines="11-12"
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
                def priority = sample_meta.quality > 40 ? 'high' : 'normal'
                return sample_meta + [priority: priority]
            }
            .view()
    ```

=== "Prima"

    ```groovy title="main.nf" linenums="2" hl_lines="11"
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
                return sample_meta
            }
            .view()
    ```

L'operatore ternario è una scorciatoia per un'istruzione if/else che segue il pattern `condizione ? valore_se_vero : valore_se_falso`. Questa riga significa: "Se la qualità è maggiore di 40, usa 'high', altrimenti usa 'normal'". Il suo cugino, l'**operatore Elvis** (`?:`), fornisce valori predefiniti quando qualcosa è null o vuoto - esploreremo quel pattern più avanti in questo tutorial.

L'operatore di addizione di mappe `+` crea una **nuova mappa** piuttosto che modificare quella esistente. Questa riga crea una nuova mappa che contiene tutte le coppie chiave-valore da `sample_meta` più la nuova chiave `priority`.

!!! Note "Nota"

    Non modificate mai le mappe passate nelle closure - create sempre nuove usando `+` (per esempio). In Nextflow, gli stessi dati spesso fluiscono attraverso più operazioni simultaneamente. Modificare una mappa in-place può causare effetti collaterali imprevedibili quando altre operazioni fanno riferimento a quello stesso oggetto. Creare nuove mappe assicura che ogni operazione abbia la propria copia pulita.

Eseguite il flusso di lavoro modificato:

```bash
nextflow run main.nf
```

??? success "Output del comando"

    ```console
    [id:sample_001, organism:human, tissue:liver, depth:30000000, quality:38.5, priority:normal]
    [id:sample_002, organism:mouse, tissue:brain, depth:25000000, quality:35.2, priority:normal]
    [id:sample_003, organism:human, tissue:kidney, depth:45000000, quality:42.1, priority:high]
    ```

Abbiamo aggiunto con successo logica condizionale per arricchire i nostri metadati con un livello di priorità basato sui punteggi di qualità.

#### 1.1.5. Creare Sottoinsiemi di Mappe con `.subMap()`

Mentre l'operatore `+` aggiunge chiavi a una mappa, a volte dovete fare l'opposto - estrarre solo chiavi specifiche. Il metodo `.subMap()` è perfetto per questo.

Aggiungiamo una riga per creare una versione semplificata dei nostri metadati che contiene solo campi di identificazione:

=== "Dopo"

    ```groovy title="main.nf" linenums="2" hl_lines="12-15"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                // Scripting per la trasformazione dei dati
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                def id_only = sample_meta.subMap(['id', 'organism', 'tissue'])
                println "ID fields only: ${id_only}"

                def priority = sample_meta.quality > 40 ? 'high' : 'normal'
                return sample_meta + [priority: priority]
            }
            .view()
    ```

=== "Prima"

    ```groovy title="main.nf" linenums="2" hl_lines="12"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                // Scripting per la trasformazione dei dati
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                def priority = sample_meta.quality > 40 ? 'high' : 'normal'
                return sample_meta + [priority: priority]
            }
            .view()
    ```

Eseguite il flusso di lavoro modificato:

```bash
nextflow run main.nf
```

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [peaceful_cori] DSL2 - revision: 4cc4a8340f

    ID fields only: [id:sample_001, organism:human, tissue:liver]
    ID fields only: [id:sample_002, organism:mouse, tissue:brain]
    ID fields only: [id:sample_003, organism:human, tissue:kidney]
    [id:sample_001, organism:human, tissue:liver, depth:30000000, quality:38.5, priority:normal]
    [id:sample_002, organism:mouse, tissue:brain, depth:25000000, quality:35.2, priority:normal]
    [id:sample_003, organism:human, tissue:kidney, depth:45000000, quality:42.1, priority:high]
    ```

Questo mostra sia i metadati completi visualizzati dall'operazione `view()` che il sottoinsieme estratto che abbiamo stampato con `println`.

Il metodo `.subMap()` prende una lista di chiavi e restituisce una nuova mappa contenente solo quelle chiavi. Se una chiave non esiste nella mappa originale, semplicemente non viene inclusa nel risultato.

Questo è particolarmente utile quando dovete creare diverse versioni di metadati per processi diversi - alcuni potrebbero aver bisogno di metadati completi mentre altri necessitano solo di campi di identificazione minimi.

Ora rimuovete quelle istruzioni println per ripristinare il vostro flusso di lavoro allo stato precedente, poiché non ne abbiamo bisogno andando avanti.

!!! tip "Riepilogo Operazioni su Mappe"

    - **Aggiungere chiavi**: `map1 + [new_key: value]` - Crea nuova mappa con chiavi aggiuntive
    - **Estrarre chiavi**: `map1.subMap(['key1', 'key2'])` - Crea nuova mappa con solo le chiavi specificate
    - **Entrambe le operazioni creano nuove mappe** - Le mappe originali rimangono invariate

#### 1.1.6. Combinare Mappe e Restituire Risultati

Finora, abbiamo restituito solo quella che la comunità Nextflow chiama 'meta map', e abbiamo ignorato i file a cui quei metadati si riferiscono. Ma se state scrivendo flussi di lavoro Nextflow, probabilmente volete fare qualcosa con quei file.

Produciamo una struttura di canale composta da una tupla di 2 elementi: la mappa di metadati arricchita e il percorso del file corrispondente. Questo è un pattern comune in Nextflow per passare dati ai processi.

=== "Dopo"

    ```groovy title="main.nf" linenums="2" hl_lines="12"
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
                def priority = sample_meta.quality > 40 ? 'high' : 'normal'
                return tuple( sample_meta + [priority: priority], file(row.file_path) )
            }
            .view()
    ```

=== "Prima"

    ```groovy title="main.nf" linenums="2" hl_lines="12"
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
                def priority = sample_meta.quality > 40 ? 'high' : 'normal'
                return sample_meta + [priority: priority]
            }
            .view()
    ```

Applicate questa modifica ed eseguite il flusso di lavoro:

```bash
nextflow run main.nf
```

??? success "Output del comando"

    ```console
    [[id:sample_001, organism:human, tissue:liver, depth:30000000, quality:38.5, priority:normal], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_001_S1_L001_R1_001.fastq]
    [[id:sample_002, organism:mouse, tissue:brain, depth:25000000, quality:35.2, priority:normal], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_002_S2_L001_R1_001.fastq]
    [[id:sample_003, organism:human, tissue:kidney, depth:45000000, quality:42.1, priority:high], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_003_S3_L001_R1_001.fastq]
    ```

Questa struttura di tupla `[meta, file]` è un pattern comune in Nextflow per passare sia metadati che file associati ai processi.

!!! note "Nota"

    **Mappe e Metadati**: Le mappe sono fondamentali per lavorare con i metadati in Nextflow. Per una spiegazione più dettagliata su come lavorare con mappe di metadati, consultate la side quest [Lavorare con i metadati](./metadata.md).

Il nostro flusso di lavoro dimostra il pattern fondamentale: le **operazioni di dataflow** (`workflow`, `channel.fromPath()`, `.splitCsv()`, `.map()`, `.view()`) orchestrano come i dati si muovono attraverso la pipeline, mentre lo **scripting** (mappe `[key: value]`, metodi di stringhe, conversioni di tipo, operatori ternari) all'interno della closure `.map()` gestisce la trasformazione dei singoli elementi di dati.

### 1.2. Comprendere Tipi Diversi: Channel vs List

Finora tutto bene, possiamo distinguere tra operazioni di dataflow e scripting. Ma che dire quando lo stesso nome di metodo esiste in entrambi i contesti?

Un esempio perfetto è il metodo `collect`, che esiste sia per i tipi di canale che per i tipi List nella libreria standard di Nextflow. Il metodo `collect()` su una List trasforma ogni elemento, mentre l'operatore `collect()` su un canale raccoglie tutte le emissioni del canale in un canale a singolo elemento.

Dimostriamolo con alcuni dati di esempio, iniziando col rinfrescarci su cosa fa l'operatore `collect()` del canale. Date un'occhiata a `collect.nf`:

```groovy title="collect.nf" linenums="1"
def sample_ids = ['sample_001', 'sample_002', 'sample_003']

// channel.collect() - raggruppa più emissioni di canale in una
ch_input = channel.fromList(sample_ids)
ch_input.view { sample -> "Individual channel item: ${sample}" }
ch_collected = ch_input.collect()
ch_collected.view { list -> "channel.collect() result: ${list} (${list.size()} items grouped into 1)" }
```

Passaggi:

- Definire una List di ID campione
- Creare un canale con `fromList()` che emette ogni ID campione separatamente
- Stampare ogni elemento con `view()` mentre fluisce attraverso
- Raccogliere tutti gli elementi in una singola lista con l'operatore `collect()` del canale
- Stampare il risultato raccolto (singolo elemento contenente tutti gli ID campione) con un secondo `view()`

Abbiamo cambiato la struttura del canale, ma non abbiamo cambiato i dati stessi.

Eseguite il flusso di lavoro per confermarlo:

```bash
nextflow run collect.nf
```

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `collect.nf` [loving_mendel] DSL2 - revision: e8d054a46e

    Individual channel item: sample_001
    Individual channel item: sample_002
    Individual channel item: sample_003
    channel.collect() result: [sample_001, sample_002, sample_003] (3 items grouped into 1)
    ```

`view()` restituisce un output per ogni emissione del canale, quindi sappiamo che questo singolo output contiene tutti e 3 gli elementi originali raggruppati in una lista.

Ora vediamo il metodo `collect` su una List in azione. Modificate `collect.nf` per applicare il metodo `collect` della List alla lista originale di ID campione:

=== "Dopo"

    ```groovy title="main.nf" linenums="1" hl_lines="9-13"
    def sample_ids = ['sample_001', 'sample_002', 'sample_003']

    // channel.collect() - raggruppa più emissioni di canale in una
    ch_input = channel.fromList(sample_ids)
    ch_input.view { sample -> "Individual channel item: ${sample}" }
    ch_collected = ch_input.collect()
    ch_collected.view { list -> "channel.collect() result: ${list} (${list.size()} items grouped into 1)" }

    // List.collect() - trasforma ogni elemento, preserva la struttura
    def formatted_ids = sample_ids.collect { id ->
        id.toUpperCase().replace('SAMPLE_', 'SPECIMEN_')
    }
    println "List.collect() result: ${formatted_ids} (${sample_ids.size()} items transformed into ${formatted_ids.size()})"
    ```

=== "Prima"

    ```groovy title="main.nf" linenums="1"
    def sample_ids = ['sample_001', 'sample_002', 'sample_003']

    // channel.collect() - raggruppa più emissioni di canale in una
    ch_input = channel.fromList(sample_ids)
    ch_input.view { sample -> "Individual channel item: ${sample}" }
    ch_collected = ch_input.collect()
    ch_collected.view { list -> "channel.collect() result: ${list} (${list.size()} items grouped into 1)" }
    ```

In questo nuovo snippet:

- Definiamo una nuova variabile `formatted_ids` che usa il metodo `collect` della List per trasformare ogni ID campione nella lista originale
- Stampiamo il risultato usando `println`

Eseguite il flusso di lavoro modificato:

```bash
nextflow run collect.nf
```

??? success "Output del comando"

    ```console hl_lines="5"
    N E X T F L O W   ~  version 25.10.2

    Launching `collect.nf` [cheeky_stonebraker] DSL2 - revision: 2d5039fb47

    List.collect() result: [SPECIMEN_001, SPECIMEN_002, SPECIMEN_003] (3 items transformed into 3)
    Individual channel item: sample_001
    Individual channel item: sample_002
    Individual channel item: sample_003
    channel.collect() result: [sample_001, sample_002, sample_003] (3 items grouped into 1)
    ```

Questa volta, NON abbiamo cambiato la struttura dei dati, abbiamo ancora 3 elementi nella lista, ma ABBIAMO trasformato ogni elemento usando il metodo `collect` della List per produrre una nuova lista con valori modificati. Questo è simile all'uso dell'operatore `map` su un canale, ma sta operando su una struttura dati List piuttosto che su un canale.

`collect` è un caso estremo che stiamo usando qui per sottolineare un punto. La lezione chiave è che quando scrivete flussi di lavoro, distinguete sempre tra **strutture dati** (List, Map, ecc.) e **canali** (costrutti di dataflow). Le operazioni possono condividere nomi ma comportarsi in modo completamente diverso a seconda del tipo su cui vengono chiamate.

### 1.3. L'Operatore Spread (`*.`) - Scorciatoia per l'Estrazione di Proprietà

Correlato al metodo `collect` della List è l'operatore spread (`*.`), che fornisce un modo conciso per estrarre proprietà dalle collezioni. È essenzialmente zucchero sintattico per un pattern `collect` comune.

Aggiungiamo una dimostrazione al nostro file `collect.nf`:

=== "Dopo"

    ```groovy title="collect.nf" linenums="1" hl_lines="15-18"
    def sample_ids = ['sample_001', 'sample_002', 'sample_003']

    // channel.collect() - raggruppa più emissioni di canale in una
    ch_input = channel.fromList(sample_ids)
    ch_input.view { sample -> "Individual channel item: ${sample}" }
    ch_collected = ch_input.collect()
    ch_collected.view { list -> "channel.collect() result: ${list} (${list.size()} items grouped into 1)" }

    // List.collect() - trasforma ogni elemento, preserva la struttura
    def formatted_ids = sample_ids.collect { id ->
        id.toUpperCase().replace('SAMPLE_', 'SPECIMEN_')
    }
    println "List.collect() result: ${formatted_ids} (${sample_ids.size()} items transformed into ${formatted_ids.size()})"

    // Operatore spread - accesso conciso alle proprietà
    def sample_data = [[id: 's1', quality: 38.5], [id: 's2', quality: 42.1], [id: 's3', quality: 35.2]]
    def all_ids = sample_data*.id
    println "Spread operator result: ${all_ids}"
    ```

=== "Prima"

    ```groovy title="collect.nf" linenums="1"
    def sample_ids = ['sample_001', 'sample_002', 'sample_003']

    // channel.collect() - raggruppa più emissioni di canale in una
    ch_input = channel.fromList(sample_ids)
    ch_input.view { sample -> "Individual channel item: ${sample}" }
    ch_collected = ch_input.collect()
    ch_collected.view { list -> "channel.collect() result: ${list} (${list.size()} items grouped into 1)" }

    // List.collect() - trasforma ogni elemento, preserva la struttura
    def formatted_ids = sample_ids.collect { id ->
        id.toUpperCase().replace('SAMPLE_', 'SPECIMEN_')
    }
    println "List.collect() result: ${formatted_ids} (${sample_ids.size()} items transformed into ${formatted_ids.size()})"
    ```

Eseguite il flusso di lavoro aggiornato:

```bash title="Test spread operator"
nextflow run collect.nf
```

??? success "Output del comando"

    ```console hl_lines="6"
    N E X T F L O W   ~  version 25.10.2

    Launching `collect.nf` [cranky_galileo] DSL2 - revision: 5f3c8b2a91

    List.collect() result: [SPECIMEN_001, SPECIMEN_002, SPECIMEN_003] (3 items transformed into 3)
    Spread operator result: [s1, s2, s3]
    Individual channel item: sample_001
    Individual channel item: sample_002
    Individual channel item: sample_003
    channel.collect() result: [sample_001, sample_002, sample_003] (3 items grouped into 1)
    ```

L'operatore spread `*.` è una scorciatoia per un pattern collect comune:

```groovy
// Questi sono equivalenti:
def ids = samples*.id
def ids = samples.collect { it.id }

// Funziona anche con chiamate di metodi:
def names = files*.getName()
def names = files.collect { it.getName() }
```

L'operatore spread è particolarmente utile quando dovete estrarre una singola proprietà da una lista di oggetti - è più leggibile che scrivere la closure `collect` completa.

!!! tip "Quando Usare Spread vs Collect"

    - **Usate spread (`*.`)** per accesso semplice alle proprietà: `samples*.id`, `files*.name`
    - **Usate collect** per trasformazioni o logica complessa: `samples.collect { it.id.toUpperCase() }`, `samples.collect { [it.id, it.quality > 40] }`

### Takeaway

In questa sezione, avete imparato:

- **Dataflow vs scripting**: Gli operatori di canale orchestrano come i dati fluiscono attraverso la vostra pipeline, mentre lo scripting trasforma i singoli elementi di dati
- **Comprendere i tipi**: Lo stesso nome di metodo (come `collect`) può comportarsi diversamente a seconda del tipo su cui viene chiamato (Channel vs List)
- **Il contesto conta**: Siate sempre consapevoli se state lavorando con canali (dataflow) o strutture dati (scripting)

Comprendere questi confini è essenziale per il debugging, la documentazione e la scrittura di flussi di lavoro manutenibili.

Successivamente ci addentreremo più a fondo nelle capacità di elaborazione delle stringhe, che sono essenziali per gestire dati del mondo reale.

---

## 2. Elaborazione di Stringhe e Generazione Dinamica di Script

Padroneggiare l'elaborazione delle stringhe separa i flussi di lavoro fragili dalle pipeline robuste. Questa sezione copre l'analisi di nomi di file complessi, la generazione dinamica di script e l'interpolazione di variabili.

### 2.1. Pattern Matching ed Espressioni Regolari

I file bioinformatici hanno spesso convenzioni di denominazione complesse che codificano metadati. Estraiamoli automaticamente usando il pattern matching con espressioni regolari.

Torneremo al nostro flusso di lavoro `main.nf` e aggiungeremo della logica di pattern matching per estrarre informazioni aggiuntive sui campioni dai nomi dei file. I file FASTQ nel nostro dataset seguono convenzioni di denominazione in stile Illumina con nomi come `SAMPLE_001_S1_L001_R1_001.fastq.gz`. Questi potrebbero sembrare criptici, ma in realtà codificano metadati utili come ID campione, numero di corsia e direzione di lettura. Useremo le capacità regex per analizzare questi nomi.

Fate la seguente modifica al vostro flusso di lavoro `main.nf` esistente:

=== "Dopo"

    ```groovy title="main.nf" linenums="4" hl_lines="10-21"
            .map { row ->
                // Scripting per la trasformazione dei dati
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
    ```

=== "Prima"

    ```groovy title="main.nf" linenums="4" hl_lines="10-11"
            .map { row ->
                // Scripting per la trasformazione dei dati
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                def priority = sample_meta.quality > 40 ? 'high' : 'normal'
                return tuple(sample_meta + [priority: priority], file(row.file_path))
            }
    ```

Questo dimostra **concetti chiave di elaborazione delle stringhe**:

1. **Letterali di espressioni regolari** usando la sintassi `~/pattern/` - questo crea un pattern regex senza dover fare l'escape dei backslash
2. **Pattern matching** con l'operatore `=~` - questo tenta di far corrispondere una stringa a un pattern regex
3. **Oggetti Matcher** che catturano gruppi con `[0][1]`, `[0][2]`, ecc. - `[0]` si riferisce all'intera corrispondenza, `[1]`, `[2]`, ecc. si riferiscono ai gruppi catturati tra parentesi

Analizziamo il pattern regex `^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$`:

| Pattern             | Corrisponde                            | Cattura                             |
| ------------------- | -------------------------------------- | ----------------------------------- |
| `^(.+)`             | Nome campione dall'inizio              | Gruppo 1: nome campione             |
| `_S(\d+)`           | Numero campione `_S1`, `_S2`, ecc.     | Gruppo 2: numero campione           |
| `_L(\d{3})`         | Numero corsia `_L001`                  | Gruppo 3: corsia (3 cifre)          |
| `_(R[12])`          | Direzione lettura `_R1` o `_R2`        | Gruppo 4: direzione lettura         |
| `_(\d{3})`          | Numero chunk `_001`                    | Gruppo 5: chunk (3 cifre)           |
| `\.fastq(?:\.gz)?$` | Estensione file `.fastq` o `.fastq.gz` | Non catturato (?: è non-catturante) |

Questo analizza le convenzioni di denominazione in stile Illumina per estrarre metadati automaticamente.

Eseguite il flusso di lavoro modificato:

```bash title="Test pattern matching"
nextflow run main.nf
```

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [clever_pauling] DSL2 - revision: 605d2058b4

    [[id:sample_001, organism:human, tissue:liver, depth:30000000, quality:38.5, sample_num:1, lane:001, read:R1, chunk:001, priority:normal], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_001_S1_L001_R1_001.fastq]
    [[id:sample_002, organism:mouse, tissue:brain, depth:25000000, quality:35.2, sample_num:2, lane:001, read:R1, chunk:001, priority:normal], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_002_S2_L001_R1_001.fastq]
    [[id:sample_003, organism:human, tissue:kidney, depth:45000000, quality:42.1, sample_num:3, lane:001, read:R1, chunk:001, priority:high], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_003_S3_L001_R1_001.fastq]
    ```

Questo mostra i metadati arricchiti dai nomi dei file.

### 2.2. Generazione Dinamica di Script nei Processi

I blocchi script dei processi sono essenzialmente stringhe multi-linea che vengono passate alla shell. Potete usare **logica condizionale** (if/else, operatori ternari) per generare dinamicamente stringhe di script diverse in base alle caratteristiche dell'input. Questo è essenziale per gestire tipi di input diversi—come letture single-end vs paired-end—senza duplicare le definizioni dei processi.

Aggiungiamo un processo al nostro flusso di lavoro che dimostri questo pattern. Aprite `modules/fastp.nf` e date un'occhiata:

```groovy title="modules/fastp.nf" linenums="1"
process FASTP {
    container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_trimmed*.fastq.gz"), emit: reads

    script:
    """
    fastp \\
        --in1 ${reads[0]} \\
        --in2 ${reads[1]} \\
        --out1 ${meta.id}_trimmed_R1.fastq.gz \\
        --out2 ${meta.id}_trimmed_R2.fastq.gz \\
        --json ${meta.id}.fastp.json \\
        --html ${meta.id}.fastp.html \\
        --thread $task.cpus
    """
}
```

Il processo prende file FASTQ come input ed esegue lo strumento `fastp` per tagliare gli adattatori e filtrare le letture di bassa qualità. Sfortunatamente, la persona che ha scritto questo processo non ha previsto le letture single-end che abbiamo nel nostro dataset di esempio. Aggiungiamolo al nostro flusso di lavoro e vediamo cosa succede:

Prima, includete il modulo alla primissima riga del vostro flusso di lavoro `main.nf`:

```groovy title="main.nf" linenums="1"
include { FASTP } from './modules/fastp.nf'
```

Poi modificate il blocco `workflow` per connettere il canale `ch_samples` al processo `FASTP`:

=== "Dopo"

    ```groovy title="main.nf" linenums="25" hl_lines="27"
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
    }
    ```

=== "Prima"

    ```groovy title="main.nf" linenums="25" hl_lines="26"
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
                return [sample_meta + file_meta + [priority: priority], file(row.file_path)]
            }
            .view()
    }
    ```

Eseguite questo flusso di lavoro modificato:

```bash
nextflow run main.nf
```

??? failure "Output del comando"

    ```console
    ERROR ~ Error executing process > 'FASTP (3)'

    Caused by:
      Process `FASTP (3)` terminated with an error exit status (255)


    Command executed:

      fastp \
          --in1 SAMPLE_003_S3_L001_R1_001.fastq \
          --in2 null \
          --out1 sample_003_trimmed_R1.fastq.gz \
          --out2 sample_003_trimmed_R2.fastq.gz \
          --json sample_003.fastp.json \
          --html sample_003.fastp.html \
          --thread 2

    Command exit status:
      255

    Command output:
      (empty)
    ```

Potete vedere che il processo sta cercando di eseguire `fastp` con un valore `null` per il secondo file di input, il che lo fa fallire. Questo perché il nostro dataset contiene letture single-end, ma il processo è codificato per aspettarsi letture paired-end (due file di input alla volta).

Risolviamo questo aggiungendo logica condizionale al blocco `script:` del processo `FASTP`. Un'istruzione if/else controlla il conteggio dei file di lettura e adatta il comando di conseguenza.

=== "Dopo"

    ```groovy title="main.nf" linenums="10" hl_lines="3-27"
        script:
        // Semplice rilevamento single-end vs paired-end
        def is_single = reads instanceof List ? reads.size() == 1 : true

        if (is_single) {
            def input_file = reads instanceof List ? reads[0] : reads
            """
            fastp \\
                --in1 ${input_file} \\
                --out1 ${meta.id}_trimmed.fastq.gz \\
                --json ${meta.id}.fastp.json \\
                --html ${meta.id}.fastp.html \\
                --thread $task.cpus
            """
        } else {
            """
            fastp \\
                --in1 ${reads[0]} \\
                --in2 ${reads[1]} \\
                --out1 ${meta.id}_trimmed_R1.fastq.gz \\
                --out2 ${meta.id}_trimmed_R2.fastq.gz \\
                --json ${meta.id}.fastp.json \\
                --html ${meta.id}.fastp.html \\
                --thread $task.cpus
            """
        }
    ```

=== "Prima"

    ```groovy title="main.nf" linenums="10" hl_lines="2-11"
            script:
            """
            fastp \\
                --in1 ${reads[0]} \\
                --in2 ${reads[1]} \\
                --out1 ${meta.id}_trimmed_R1.fastq.gz \\
                --out2 ${meta.id}_trimmed_R2.fastq.gz \\
                --json ${meta.id}.fastp.json \\
                --html ${meta.id}.fastp.html \\
                --thread $task.cpus
            """
        }
    ```

Ora il flusso di lavoro può gestire con eleganza sia letture single-end che paired-end. La logica condizionale controlla il numero di file di input e costruisce il comando appropriato per `fastp`. Vediamo se funziona:

```bash
nextflow run main.nf
```

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [adoring_rosalind] DSL2 - revision: 04b1cd93e9

    executor >  local (3)
    [31/a8ad4d] process > FASTP (3) [100%] 3 of 3 ✔
    ```

Sembra buono! Se controlliamo i comandi effettivi che sono stati eseguiti (personalizzate per il vostro hash di attività):

```console title="Check commands executed"
cat work/31/a8ad4d95749e685a6d842d3007957f/.command.sh
```

Possiamo vedere che Nextflow ha correttamente scelto il comando giusto per le letture single-end:

```bash title=".command.sh"
#!/bin/bash -ue
fastp \
    --in1 SAMPLE_003_S3_L001_R1_001.fastq \
    --out1 sample_003_trimmed.fastq.gz \
    --json sample_003.fastp.json \
    --html sample_003.fastp.html \
    --thread 2
```

Un altro uso comune della logica dinamica degli script può essere visto in [il modulo Genomics di Nextflow for Science](../../nf4science/genomics/02_joint_calling). In quel modulo, il processo GATK chiamato può prendere più file di input, ma ognuno deve essere prefissato con `-V` per formare una riga di comando corretta. Il processo usa lo scripting per trasformare una collezione di file di input (`all_gvcfs`) negli argomenti di comando corretti:

```groovy title="command line manipulation for GATK" linenums="1" hl_lines="2 5"
    script:
    def gvcfs_line = all_gvcfs.collect { gvcf -> "-V ${gvcf}" }.join(' ')
    """
    gatk GenomicsDBImport \
        ${gvcfs_line} \
        -L ${interval_list} \
        --genomicsdb-workspace-path ${cohort_name}_gdb
    """
```

Questi pattern di uso dello scripting nei blocchi script dei processi sono estremamente potenti e possono essere applicati in molti scenari - dalla gestione di tipi di input variabili alla costruzione di argomenti complessi da linea di comando da collezioni di file, rendendo i vostri processi veramente adattabili ai requisiti diversi dei dati del mondo reale.

### 2.3. Interpolazione di Variabili: Variabili Nextflow e Shell

Gli script dei processi mescolano variabili Nextflow, variabili shell e sostituzioni di comandi, ognuna con sintassi di interpolazione diversa. Usare la sintassi sbagliata causa errori. Esploriamo questi con un processo che crea un report di elaborazione.

Date un'occhiata al file del modulo `modules/generate_report.nf`:

```groovy title="modules/generate_report.nf" linenums="1"
process GENERATE_REPORT {

    publishDir 'results/reports', mode: 'copy'

    input:
    tuple val(meta), path(reads)

    output:
    path "${meta.id}_report.txt"

    script:
    """
    echo "Processing ${reads}" > ${meta.id}_report.txt
    echo "Sample: ${meta.id}" >> ${meta.id}_report.txt
    """
}
```

Questo processo scrive un semplice report con l'ID del campione e il nome del file. Ora eseguiamolo per vedere cosa succede quando dobbiamo mescolare diversi tipi di variabili.

Includete il processo nel vostro `main.nf` e aggiungetelo al flusso di lavoro:

=== "Dopo"

    ```groovy title="main.nf" linenums="1" hl_lines="2 30"
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

=== "Prima"

    ```groovy title="main.nf" linenums="1" hl_lines="1 10-29"
    include { FASTP } from './modules/fastp.nf'

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
    }
    ```

Ora eseguite il flusso di lavoro e controllate i report generati in `results/reports/`. Dovrebbero contenere informazioni di base su ogni campione.

??? success "Output del comando"

    ```console
    <!-- TODO: output -->
    ```

Ma cosa succede se vogliamo aggiungere informazioni su quando e dove è avvenuta l'elaborazione? Modifichiamo il processo per usare variabili **shell** e un po' di sostituzione di comandi per includere l'utente corrente, l'hostname e la data nel report:

=== "Dopo"

    ```groovy title="modules/generate_report.nf" linenums="10" hl_lines="5-7"
        script:
        """
        echo "Processing ${reads}" > ${meta.id}_report.txt
        echo "Sample: ${meta.id}" >> ${meta.id}_report.txt
        echo "Processed by: ${USER}" >> ${meta.id}_report.txt
        echo "Hostname: $(hostname)" >> ${meta.id}_report.txt
        echo "Date: $(date)" >> ${meta.id}_report.txt
        """
    ```

=== "Prima"

    ```groovy title="modules/generate_report.nf" linenums="10"
        script:
        """
        echo "Processing ${reads}" > ${meta.id}_report.txt
        echo "Sample: ${meta.id}" >> ${meta.id}_report.txt
        """
    ```

Se eseguite questo, noterete un errore - Nextflow cerca di interpretare `${USER}` come una variabile Nextflow che non esiste.

??? failure "Output del comando"

    ```console
    Error modules/generate_report.nf:15:27: `USER` is not defined
    │  15 |     echo "Processed by: ${USER}" >> ${meta.id}_report.txt
    ╰     |                           ^^^^

    ERROR ~ Script compilation failed
    ```

Dobbiamo fare l'escape in modo che Bash possa gestirlo invece.

Risolviamo questo facendo l'escape delle variabili shell e delle sostituzioni di comandi con un backslash (`\`):

=== "Dopo"

    ```groovy title="modules/generate_report.nf" linenums="10" hl_lines="5-7"
        script:
        """
        echo "Processing ${reads}" > ${meta.id}_report.txt
        echo "Sample: ${meta.id}" >> ${meta.id}_report.txt
        echo "Processed by: \${USER}" >> ${meta.id}_report.txt
        echo "Hostname: \$(hostname)" >> ${meta.id}_report.txt
        echo "Date: \$(date)" >> ${meta.id}_report.txt
        """
    ```

=== "Prima"

    ```groovy title="modules/generate_report.nf" linenums="10"
        script:
        """
        echo "Processing ${reads}" > ${meta.id}_report.txt
        echo "Sample: ${meta.id}" >> ${meta.id}_report.txt
        echo "Processed by: ${USER}" >> ${meta.id}_report.txt
        echo "Hostname: $(hostname)" >> ${meta.id}_report.txt
        echo "Date: $(date)" >> ${meta.id}_report.txt
        """
    ```

Ora funziona! Il backslash (`\`) dice a Nextflow "non interpretare questo, passalo a Bash."

### Takeaway

In questa sezione, avete imparato tecniche di **elaborazione delle stringhe**:

- **Espressioni regolari per l'analisi dei file**: Usare l'operatore `=~` e pattern regex (`~/pattern/`) per estrarre metadati da convenzioni di denominazione dei file complesse
- **Generazione dinamica di script**: Usare logica condizionale (if/else, operatori ternari) per generare stringhe di script diverse in base alle caratteristiche dell'input
- **Interpolazione di variabili**: Capire quando Nextflow interpreta le stringhe vs quando lo fa la shell
  - `${var}` - Variabili Nextflow (interpolate da Nextflow al momento della compilazione del flusso di lavoro)
  - `\${var}` - Variabili d'ambiente shell (con escape, passate a bash al runtime)
  - `\$(cmd)` - Sostituzione di comandi shell (con escape, eseguita da bash al runtime)

Questi pattern di elaborazione e generazione di stringhe sono essenziali per gestire i diversi formati di file e convenzioni di denominazione che incontrerete nei flussi di lavoro bioinformatici del mondo reale.

---

## 3. Creare Funzioni Riutilizzabili

La logica complessa del flusso di lavoro inline negli operatori di canale o nelle definizioni dei processi riduce la leggibilità e la manutenibilità. Le **funzioni** vi permettono di estrarre questa logica in componenti con nome e riutilizzabili.

La nostra operazione map è diventata lunga e complessa. Estraiamola in una funzione riutilizzabile usando la parola chiave `def`.

Per illustrare come appare con il nostro flusso di lavoro esistente, fate la modifica qui sotto, usando `def` per definire una funzione riutilizzabile chiamata `separateMetadata`:

=== "Dopo"

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

=== "Prima"

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

Estraendo questa logica in una funzione, abbiamo ridotto la logica effettiva del flusso di lavoro a qualcosa di molto più pulito:

```groovy title="minimal workflow"
    ch_samples = channel.fromPath("./data/samples.csv")
        .splitCsv(header: true)
        .map{ row -> separateMetadata(row) }

    ch_fastp = FASTP(ch_samples)
    GENERATE_REPORT(ch_samples)
```

Questo rende la logica del flusso di lavoro molto più facile da leggere e capire a colpo d'occhio. La funzione `separateMetadata` incapsula tutta la logica complessa per analizzare e arricchire i metadati, rendendola riutilizzabile e testabile.

Eseguite il flusso di lavoro per assicurarvi che funzioni ancora:

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

L'output dovrebbe mostrare entrambi i processi completati con successo. Il flusso di lavoro è ora molto più pulito e facile da mantenere, con tutta la logica complessa di elaborazione dei metadati incapsulata nella funzione `separateMetadata`.

### Takeaway

In questa sezione, avete imparato la **creazione di funzioni**:

- **Definire funzioni con `def`**: La parola chiave per creare funzioni con nome (come `def` in Python o `function` in JavaScript)
- **Scope delle funzioni**: Le funzioni definite a livello di script sono accessibili in tutto il vostro flusso di lavoro Nextflow
- **Valori di ritorno**: Le funzioni restituiscono automaticamente l'ultima espressione, o usano `return` esplicito
- **Codice più pulito**: Estrarre logica complessa in funzioni è una pratica fondamentale di ingegneria del software in qualsiasi linguaggio

Successivamente, esploreremo come usare le closure nelle direttive dei processi per l'allocazione dinamica delle risorse.

---

## 4. Direttive di Risorse Dinamiche con Closure

Finora abbiamo usato lo scripting nel blocco `script` dei processi. Ma le **closure** (introdotte nella Sezione 1.1) sono anche incredibilmente utili nelle direttive dei processi, specialmente per l'allocazione dinamica delle risorse. Aggiungiamo direttive di risorse al nostro processo FASTP che si adattano in base alle caratteristiche del campione.

### 4.1. Allocazione di risorse specifica per campione

Attualmente, il nostro processo FASTP usa risorse predefinite. Rendiamolo più intelligente allocando più CPU per campioni ad alta profondità. Modificate `modules/fastp.nf` per includere una direttiva `cpus` dinamica e una direttiva `memory` statica:

=== "Dopo"

    ```groovy title="modules/fastp.nf" linenums="1" hl_lines="4-5"
    process FASTP {
        container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

        cpus { meta.depth > 40000000 ? 2 : 1 }
        memory 2.GB

        input:
        tuple val(meta), path(reads)
    ```

=== "Prima"

    ```groovy title="modules/fastp.nf" linenums="1"
    process FASTP {
        container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

        input:
        tuple val(meta), path(reads)
    ```

La closure `{ meta.depth > 40000000 ? 2 : 1 }` usa l'**operatore ternario** (trattato nella Sezione 1.1) e viene valutata per ogni attività, permettendo l'allocazione di risorse per campione. I campioni ad alta profondità (>40M letture) ottengono 2 CPU, mentre gli altri ottengono 1 CPU.

!!! note "Accedere alle Variabili di Input nelle Direttive"

    La closure può accedere a qualsiasi variabile di input (come `meta` qui) perché Nextflow valuta queste closure nel contesto di ogni esecuzione di attività.

Eseguite nuovamente il flusso di lavoro con l'opzione `-ansi-log false` per rendere più facile vedere gli hash delle attività.

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

Potete controllare l'esatto comando `docker` che è stato eseguito per vedere l'allocazione di CPU per una data attività:

```console title="Check docker command"
cat work/48/6db0c9e9d8aa65e4bb4936cd3bd59e/.command.run | grep "docker run"
```

Dovreste vedere qualcosa come:

```bash title="docker command"
    docker run -i --cpu-shares 4096 --memory 2048m -e "NXF_TASK_WORKDIR" -v /workspaces/training/side-quests/essential_scripting_patterns:/workspaces/training/side-quests/essential_scripting_patterns -w "$NXF_TASK_WORKDIR" --name $NXF_BOXID community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690 /bin/bash -ue /workspaces/training/side-quests/essential_scripting_patterns/work/48/6db0c9e9d8aa65e4bb4936cd3bd59e/.command.sh
```

In questo esempio abbiamo scelto un esempio che ha richiesto 2 CPU (`--cpu-shares 2048`), perché era un campione ad alta profondità, ma dovreste vedere allocazioni di CPU diverse a seconda della profondità del campione. Provate questo anche per le altre attività.

### 4.2. Strategie di retry

Un altro pattern potente è usare `task.attempt` per strategie di retry. Per mostrare perché questo è utile, inizieremo riducendo l'allocazione di memoria a FASTP a meno di quanto necessario. Cambiate la direttiva `memory` in `modules/fastp.nf` a `1.GB`:

=== "Dopo"

    ```groovy title="modules/fastp.nf" linenums="1" hl_lines="5"
    process FASTP {
        container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

        cpus { meta.depth > 40000000 ? 4 : 2 }
        memory 1.GB

        input:
        tuple val(meta), path(reads)
    ```

=== "Prima"

    ```groovy title="modules/fastp.nf" linenums="1" hl_lines="5"
    process FASTP {
        container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

        cpus { meta.depth > 40000000 ? 4 : 2 }
        memory 2.GB

        input:
        tuple val(meta), path(reads)
    ```

... ed eseguite nuovamente il flusso di lavoro:

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

Questo è uno scenario molto comune nei flussi di lavoro del mondo reale - a volte semplicemente non sapete quanta memoria avrà bisogno un'attività finché non la eseguite.

Per rendere il nostro flusso di lavoro più robusto, possiamo implementare una strategia di retry che aumenta l'allocazione di memoria ad ogni tentativo, ancora una volta usando una closure Groovy. Modificate la direttiva `memory` per moltiplicare la memoria base per `task.attempt`, e aggiungete le direttive `errorStrategy 'retry'` e `maxRetries 2`:

=== "Dopo"

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

=== "Prima"

    ```groovy title="modules/fastp.nf" linenums="1" hl_lines="5"
    process FASTP {
        container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

        cpus { meta.depth > 40000000 ? 4 : 2 }
        memory 2.GB

        input:
        tuple val(meta), path(reads)
    ```

Ora se il processo fallisce a causa di memoria insufficiente, Nextflow riproverà con più memoria:

- Primo tentativo: 1 GB (task.attempt = 1)
- Secondo tentativo: 2.GB (task.attempt = 2)

... e così via, fino al limite `maxRetries`.

### Takeaway

Le direttive dinamiche con closure vi permettono di:

- Allocare risorse in base alle caratteristiche dell'input
- Implementare strategie di retry automatiche con risorse crescenti
- Combinare più fattori (metadati, numero di tentativi, priorità)
- Usare logica condizionale per calcoli complessi delle risorse

Questo rende i vostri flussi di lavoro sia più efficienti (non sovra-allocando) che più robusti (retry automatico con più risorse).

---

## 5. Logica Condizionale e Controllo dei Processi

In precedenza, abbiamo usato `.map()` con scripting per trasformare i dati del canale. Ora useremo la logica condizionale per controllare quali processi vengono eseguiti in base ai dati—essenziale per flussi di lavoro flessibili che si adattano a diversi tipi di campioni.

Gli [operatori di dataflow](https://www.nextflow.io/docs/latest/reference/operator.html) di Nextflow prendono closure valutate al runtime, abilitando la logica condizionale per guidare le decisioni del flusso di lavoro in base al contenuto del canale.

### 5.1. Routing con `.branch()`

Per esempio, fingiamo che i nostri campioni di sequenziamento debbano essere tagliati con FASTP solo se sono campioni umani con una copertura sopra una certa soglia. I campioni di topo o i campioni a bassa copertura dovrebbero essere eseguiti con Trimgalore invece (questo è un esempio forzato, ma illustra il punto).

Abbiamo fornito un semplice processo Trimgalore in `modules/trimgalore.nf`, date un'occhiata se volete, ma i dettagli non sono importanti per questo esercizio. Il punto chiave è che vogliamo instradare i campioni in base ai loro metadati.

Includete il nuovo modulo da `modules/trimgalore.nf`:

=== "Dopo"

    ```groovy title="main.nf" linenums="1" hl_lines="2"
    include { FASTP } from './modules/fastp.nf'
    include { TRIMGALORE } from './modules/trimgalore.nf'
    ```

=== "Prima"

    ```groovy title="main.nf" linenums="1"
    include { FASTP } from './modules/fastp.nf'
    ```

... e poi modificate il vostro flusso di lavoro `main.nf` per ramificare i campioni in base ai loro metadati e instradarli attraverso il processo di trimming appropriato, così:

=== "Dopo"

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

=== "Prima"

    ```groovy title="main.nf" linenums="28" hl_lines="5"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row -> separateMetadata(row) }

        ch_fastp = FASTP(ch_samples)
        GENERATE_REPORT(ch_samples)
    ```

Eseguite questo flusso di lavoro modificato:

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

Qui, abbiamo usato espressioni condizionali piccole ma potenti all'interno dell'operatore `.branch{}` per instradare i campioni in base ai loro metadati. I campioni umani con alta copertura passano attraverso `FASTP`, mentre tutti gli altri campioni passano attraverso `TRIMGALORE`.

### 5.2. Usare `.filter()` con Truthiness

Un altro pattern potente per controllare l'esecuzione del flusso di lavoro è l'operatore `.filter()`, che usa una closure per determinare quali elementi dovrebbero continuare lungo la pipeline. All'interno della closure filter, scriverete **espressioni booleane** che decidono quali elementi passano.

Nextflow (come molti linguaggi dinamici) ha un concetto di **"truthiness"** che determina quali valori vengono valutati come `true` o `false` in contesti booleani:

- **Truthy**: Valori non-null, stringhe non vuote, numeri non-zero, collezioni non vuote
- **Falsy**: `null`, stringhe vuote `""`, zero `0`, collezioni vuote `[]` o `[:]`, `false`

Questo significa che `meta.id` da solo (senza esplicito `!= null`) controlla se l'ID esiste e non è vuoto. Usiamo questo per filtrare i campioni che non soddisfano i nostri requisiti di qualità.

Aggiungete quanto segue prima dell'operazione branch:

=== "Dopo"

    ```groovy title="main.nf" linenums="28" hl_lines="5-11"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row -> separateMetadata(row) }

        // Filtra campioni non validi o di bassa qualità
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

=== "Prima"

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

Eseguite nuovamente il flusso di lavoro:

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

Poiché abbiamo scelto un filtro che esclude alcuni campioni, sono state eseguite meno attività.

L'espressione filter `meta.id && meta.organism && meta.depth >= 25000000` combina truthiness con confronti espliciti:

- `meta.id && meta.organism` controlla che entrambi i campi esistano e non siano vuoti (usando truthiness)
- `meta.depth >= 25000000` assicura una profondità di sequenziamento sufficiente con un confronto esplicito

!!! note "Truthiness in Pratica"

    L'espressione `meta.id && meta.organism` è più concisa che scrivere:
    ```groovy
    meta.id != null && meta.id != '' && meta.organism != null && meta.organism != ''
    ```

    Questo rende la logica di filtraggio molto più pulita e facile da leggere.

### Takeaway

In questa sezione, avete imparato a usare la logica condizionale per controllare l'esecuzione del flusso di lavoro usando le interfacce closure degli operatori Nextflow come `.branch{}` e `.filter{}`, sfruttando la truthiness per scrivere espressioni condizionali concise.

La nostra pipeline ora instrada intelligentemente i campioni attraverso processi appropriati, ma i flussi di lavoro di produzione devono gestire con eleganza i dati non validi. Rendiamo il nostro flusso di lavoro robusto contro valori mancanti o null.

---

## 6. Operatori di Navigazione Sicura ed Elvis

La nostra funzione `separateMetadata` attualmente presume che tutti i campi CSV siano presenti e validi. Ma cosa succede con dati incompleti? Scopriamolo.

### 6.1. Il Problema: Accedere a Proprietà che Non Esistono

Diciamo che vogliamo aggiungere supporto per informazioni opzionali sulla corsa di sequenziamento. In alcuni laboratori, i campioni potrebbero avere un campo aggiuntivo per l'ID della corsa di sequenziamento o il numero di batch, ma il nostro CSV attuale non ha questa colonna. Proviamo ad accedervi comunque.

Modificate la funzione `separateMetadata` per includere un campo run_id:

=== "Dopo"

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

=== "Prima"

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

Ora eseguite il flusso di lavoro:

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

Questo va in crash con una NullPointerException.

Il problema è che `row.run_id` restituisce `null` perché la colonna `run_id` non esiste nel nostro CSV. Quando proviamo a chiamare `.toUpperCase()` su `null`, va in crash. Qui è dove l'operatore di navigazione sicura salva la situazione.

### 6.2. Operatore di Navigazione Sicura (`?.`)

L'operatore di navigazione sicura (`?.`) restituisce `null` invece di lanciare un'eccezione quando chiamato su un valore `null`. Se l'oggetto prima di `?.` è `null`, l'intera espressione viene valutata come `null` senza eseguire il metodo.

Aggiornate la funzione per usare la navigazione sicura:

=== "Dopo"

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

=== "Prima"

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

Eseguite di nuovo:

```bash
nextflow run main.nf
```

??? success "Output del comando"

    ```console
    <!-- TODO: output -->
    ```

Nessun crash! Il flusso di lavoro ora gestisce con eleganza il campo mancante. Quando `row.run_id` è `null`, l'operatore `?.` previene la chiamata `.toUpperCase()`, e `run_id` diventa `null` invece di causare un'eccezione.

### 6.3. Operatore Elvis (`?:`) per Valori Predefiniti

L'operatore Elvis (`?:`) fornisce valori predefiniti quando il lato sinistro è "falsy" (come spiegato in precedenza). È chiamato così per Elvis Presley perché `?:` assomiglia ai suoi famosi capelli e occhi quando visto di lato!

Ora che stiamo usando la navigazione sicura, `run_id` sarà `null` per i campioni senza quel campo. Usiamo l'operatore Elvis per fornire un valore predefinito e aggiungiamolo alla nostra mappa `sample_meta`:

=== "Dopo"

    ```groovy title="main.nf" linenums="5" hl_lines="9-10"
    def separateMetadata(row) {
        def sample_meta = [
            id: row.sample_id.toLowerCase(),
            organism: row.organism,
            tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
            depth: row.sequencing_depth.toInteger(),
            quality: row.quality_score.toDouble()
        ]
        def run_id = row.run_id?.toUpperCase() ?: 'UNSPECIFIED'
        sample_meta.run = run_id
    ```

=== "Prima"

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

Aggiungete anche un operatore `view()` nel flusso di lavoro per vedere i risultati:

=== "Dopo"

    ```groovy title="main.nf" linenums="30" hl_lines="4"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map{ row -> separateMetadata(row) }
            .view()
    ```

=== "Prima"

    ```groovy title="main.nf" linenums="30"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map{ row -> separateMetadata(row) }
    ```

ed eseguite il flusso di lavoro:

```bash
nextflow run main.nf
```

??? success "Output del comando"

    ```console
    [[id:sample_001, organism:human, tissue:liver, depth:30000000, quality:38.5, run:UNSPECIFIED, sample_num:1, lane:001, read:R1, chunk:001, priority:normal], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_001_S1_L001_R1_001.fastq]
    [[id:sample_002, organism:mouse, tissue:brain, depth:25000000, quality:35.2, run:UNSPECIFIED, sample_num:2, lane:001, read:R1, chunk:001, priority:normal], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_002_S2_L001_R1_001.fastq]
    [[id:sample_003, organism:human, tissue:kidney, depth:45000000, quality:42.1, run:UNSPECIFIED, sample_num:3, lane:001, read:R1, chunk:001, priority:high], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_003_S3_L001_R1_001.fastq]
    ```

Perfetto! Ora tutti i campioni hanno un campo `run` con il loro ID di corsa effettivo (in maiuscolo) o il valore predefinito 'UNSPECIFIED'. La combinazione di `?.` e `?:` fornisce sia sicurezza (nessun crash) che valori predefiniti sensati.

Togliete l'operatore `.view()` ora che abbiamo confermato che funziona.

!!! tip "Combinare Navigazione Sicura ed Elvis"

    Il pattern `value?.method() ?: 'default'` è comune nei flussi di lavoro di produzione:

    - `value?.method()` - Chiama il metodo in modo sicuro, restituisce `null` se `value` è `null`
    - `?: 'default'` - Fornisce un fallback se il risultato è `null`

    Questo pattern gestisce con eleganza dati mancanti/incompleti.

Usate questi operatori in modo coerente in funzioni, closure di operatori (`.map{}`, `.filter{}`), script di processi e file di configurazione. Prevengono crash quando si gestiscono dati del mondo reale.

### Takeaway

- **Navigazione sicura (`?.`)**: Previene crash su valori null - restituisce null invece di lanciare un'eccezione
- **Operatore Elvis (`?:`)**: Fornisce valori predefiniti - `value ?: 'default'`
- **Combinazione**: `value?.method() ?: 'default'` è il pattern comune

Questi operatori rendono i flussi di lavoro resilienti a dati incompleti - essenziale per il lavoro nel mondo reale.

---

## 7. Validazione con `error()` e `log.warn`

A volte dovete fermare immediatamente il flusso di lavoro se i parametri di input non sono validi. In Nextflow, potete usare funzioni integrate come `error()` e `log.warn`, così come costrutti di programmazione standard come istruzioni `if` e logica booleana, per implementare logica di validazione. Aggiungiamo validazione al nostro flusso di lavoro.

Create una funzione di validazione prima del vostro blocco workflow, chiamatela dal workflow, e cambiate la creazione del canale per usare un parametro per il percorso del file CSV. Se il parametro manca o il file non esiste, chiamate `error()` per fermare l'esecuzione con un messaggio chiaro.

=== "Dopo"

    ```groovy title="main.nf" linenums="1" hl_lines="5-20 23-24"
    include { FASTP } from './modules/fastp.nf'
    include { TRIMGALORE } from './modules/trimgalore.nf'
    include { GENERATE_REPORT } from './modules/generate_report.nf'

    def validateInputs() {
        // Controlla che il parametro di input sia fornito
        if (!params.input) {
            error("Input CSV file path not provided. Please specify --input <file.csv>")
        }

        // Controlla che il file CSV esista
        if (!file(params.input).exists()) {
            error("Input CSV file not found: ${params.input}")
        }
    }
    ...
    workflow {
        validateInputs()
        ch_samples = channel.fromPath(params.input)
    ```

=== "Prima"

    ```groovy title="main.nf" linenums="1"
    include { FASTP } from './modules/fastp.nf'
    include { TRIMGALORE } from './modules/trimgalore.nf'
    include { GENERATE_REPORT } from './modules/generate_report.nf'

    ...
    workflow {
        ch_samples = channel.fromPath("./data/samples.csv")
    ```

Ora provate a eseguire senza il file CSV:

```bash
nextflow run main.nf
```

??? failure "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [confident_coulomb] DSL2 - revision: 07059399ed

    WARN: Access to undefined parameter `input` -- Initialise it to a default value eg. `params.input = some_value`
    Input CSV file path not provided. Please specify --input <file.csv>
    ```

Il flusso di lavoro si ferma immediatamente con un messaggio di errore chiaro invece di fallire misteriosamente più tardi

Ora eseguitelo con un file non esistente:

```bash
nextflow run main.nf --input ./data/nonexistent.csv
```

??? failure "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [cranky_gates] DSL2 - revision: 26839ae3eb

    Input CSV file not found: ./data/nonexistent.csv
    ```

Infine, eseguitelo con il file corretto:

```bash
nextflow run main.nf --input ./data/samples.csv
```

??? success "Output del comando"

    ```console
    <!-- TODO: output -->
    ```

Questa volta viene eseguito con successo.

Potete anche aggiungere validazione all'interno della funzione `separateMetadata`. Usiamo il non-fatale `log.warn` per emettere avvisi per campioni con bassa profondità di sequenziamento, ma permettendo comunque al flusso di lavoro di continuare:

=== "Dopo"

    ```groovy title="main.nf" linenums="1" hl_lines="3-6"
        def priority = sample_meta.quality > 40 ? 'high' : 'normal'

        // Valida che i dati abbiano senso
        if (sample_meta.depth < 30000000) {
            log.warn "Low sequencing depth for ${sample_meta.id}: ${sample_meta.depth}"
        }

        return tuple(sample_meta + file_meta + [priority: priority], fastq_path)
    }
    ```

=== "Prima"

    ```groovy title="main.nf" linenums="1"
        def priority = sample_meta.quality > 40 ? 'high' : 'normal'

        return tuple(sample_meta + file_meta + [priority: priority], fastq_path)
    }
    ```

Eseguite nuovamente il flusso di lavoro con il CSV originale:

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
    WARN: Low sequencing depth for sample_002: 25000000
    ```

Vediamo un avviso sulla bassa profondità di sequenziamento per uno dei campioni.

### Takeaway

- **`error()`**: Ferma immediatamente il flusso di lavoro con un messaggio chiaro
- **`log.warn`**: Emette avvisi senza fermare il flusso di lavoro
- **Validazione precoce**: Controllate gli input prima dell'elaborazione per fallire velocemente con errori utili
- **Funzioni di validazione**: Create logica di validazione riutilizzabile che può essere chiamata all'avvio del flusso di lavoro

Una validazione appropriata rende i flussi di lavoro più robusti e user-friendly catturando i problemi precocemente con messaggi di errore chiari.

---

## 8. Handler di Eventi del Flusso di Lavoro

Fino ad ora, abbiamo scritto codice nei nostri script di flusso di lavoro e nelle definizioni dei processi. Ma c'è un'altra caratteristica importante che dovreste conoscere: gli handler di eventi del flusso di lavoro.

Gli handler di eventi sono closure che vengono eseguite in punti specifici del ciclo di vita del vostro flusso di lavoro. Sono perfetti per aggiungere logging, notifiche o operazioni di pulizia. Questi handler dovrebbero essere definiti nel vostro script di flusso di lavoro insieme alla vostra definizione di workflow.

### 8.1. L'Handler `onComplete`

L'handler di eventi più comunemente usato è `onComplete`, che viene eseguito quando il vostro flusso di lavoro termina (sia che abbia avuto successo o fallito). Aggiungiamone uno per riassumere i risultati della nostra pipeline.

Aggiungete l'handler di eventi al vostro file `main.nf`, all'interno della vostra definizione di workflow:

=== "Dopo"

    ```groovy title="main.nf" linenums="66" hl_lines="5-16"
        ch_fastp = FASTP(trim_branches.fastp)
        ch_trimgalore = TRIMGALORE(trim_branches.trimgalore)
        GENERATE_REPORT(ch_samples)

        workflow.onComplete = {
            println ""
            println "Pipeline execution summary:"
            println "=========================="
            println "Completed at: ${workflow.complete}"
            println "Duration    : ${workflow.duration}"
            println "Success     : ${workflow.success}"
            println "workDir     : ${workflow.workDir}"
            println "exit status : ${workflow.exitStatus}"
            println ""
        }
    }
    ```

=== "Prima"

    ```groovy title="main.nf" linenums="66" hl_lines="4"
        ch_fastp = FASTP(trim_branches.fastp)
        ch_trimgalore = TRIMGALORE(trim_branches.trimgalore)
        GENERATE_REPORT(ch_samples)
    }
    ```

Questa closure viene eseguita quando il flusso di lavoro si completa. All'interno, avete accesso all'oggetto `workflow` che fornisce proprietà utili sull'esecuzione.

Eseguite il vostro flusso di lavoro e vedrete questo riepilogo apparire alla fine!

```bash
nextflow run main.nf --input ./data/samples.csv -ansi-log false
```

??? success "Output del comando"

    ```console
    N E X T F L O W  ~  version 25.10.2
    Launching `main.nf` [marvelous_boltzmann] DSL2 - revision: a31662a7c1
    WARN: Low sequencing depth for sample_002: 25000000
    [9b/d48e40] Submitted process > FASTP (2)
    [6a/73867a] Submitted process > GENERATE_REPORT (2)
    [79/ad0ac5] Submitted process > GENERATE_REPORT (1)
    [f3/bda6cb] Submitted process > FASTP (1)
    [34/d5b52f] Submitted process > GENERATE_REPORT (3)

    Pipeline execution summary:
    ==========================
    Completed at: 2025-10-10T12:14:24.885384+01:00
    Duration    : 2.9s
    Success     : true
    workDir     : /workspaces/training/side-quests/essential_scripting_patterns/work
    exit status : 0
    ```

Rendiamolo più utile aggiungendo logica condizionale:

=== "Dopo"

    ```groovy title="main.nf" linenums="66" hl_lines="5-22"
        ch_fastp = FASTP(trim_branches.fastp)
        ch_trimgalore = TRIMGALORE(trim_branches.trimgalore)
        GENERATE_REPORT(ch_samples)

        workflow.onComplete = {
            println ""
            println "Pipeline execution summary:"
            println "=========================="
            println "Completed at: ${workflow.complete}"
            println "Duration    : ${workflow.duration}"
            println "Success     : ${workflow.success}"
            println "workDir     : ${workflow.workDir}"
            println "exit status : ${workflow.exitStatus}"
            println ""

            if (workflow.success) {
                println "✅ Pipeline completed successfully!"
            } else {
                println "❌ Pipeline failed!"
                println "Error: ${workflow.errorMessage}"
            }
        }
    }
    ```

=== "Prima"

    ```groovy title="main.nf" linenums="66" hl_lines="5-16"
        ch_fastp = FASTP(trim_branches.fastp)
        ch_trimgalore = TRIMGALORE(trim_branches.trimgalore)
        GENERATE_REPORT(ch_samples)

        workflow.onComplete = {
            println ""
            println "Pipeline execution summary:"
            println "=========================="
            println "Completed at: ${workflow.complete}"
            println "Duration    : ${workflow.duration}"
            println "Success     : ${workflow.success}"
            println "workDir     : ${workflow.workDir}"
            println "exit status : ${workflow.exitStatus}"
            println ""
        }
    }
    ```

Ora otteniamo un riepilogo ancora più informativo, incluso un messaggio di successo/fallimento e la directory di output se specificata:

??? success "Output del comando"

    ```console
    N E X T F L O W  ~  version 25.10.2
    Launching `main.nf` [boring_linnaeus] DSL2 - revision: a31662a7c1
    WARN: Low sequencing depth for sample_002: 25000000
    [e5/242efc] Submitted process > FASTP (2)
    [3b/74047c] Submitted process > GENERATE_REPORT (3)
    [8a/7a57e6] Submitted process > GENERATE_REPORT (1)
    [a8/b1a31f] Submitted process > GENERATE_REPORT (2)
    [40/648429] Submitted process > FASTP (1)

    Pipeline execution summary:
    ==========================
    Completed at: 2025-10-10T12:16:00.522569+01:00
    Duration    : 3.6s
    Success     : true
    workDir     : /workspaces/training/side-quests/essential_scripting_patterns/work
    exit status : 0

    ✅ Pipeline completed successfully!
    ```

Potete anche scrivere il riepilogo in un file usando operazioni sui file:

```groovy title="main.nf - Writing summary to file"
workflow {
    // ... il vostro codice del flusso di lavoro ...

    workflow.onComplete = {
        def summary = """
        Pipeline Execution Summary
        ===========================
        Completed: ${workflow.complete}
        Duration : ${workflow.duration}
        Success  : ${workflow.success}
        Command  : ${workflow.commandLine}
        """

        println summary

        // Scrivi su un file di log
        def log_file = file("${workflow.launchDir}/pipeline_summary.txt")
        log_file.text = summary
    }
}
```

### 8.2. L'Handler `onError`

Oltre a `onComplete`, c'è un altro handler di eventi che potete usare: `onError`, che viene eseguito solo se il flusso di lavoro fallisce:

```groovy title="main.nf - onError handler"
workflow {
    // ... il vostro codice del flusso di lavoro ...

    workflow.onError = {
        println "="* 50
        println "Pipeline execution failed!"
        println "Error message: ${workflow.errorMessage}"
        println "="* 50

        // Scrivi log di errore dettagliato
        def error_file = file("${workflow.launchDir}/error.log")
        error_file.text = """
        Workflow Error Report
        =====================
        Time: ${new Date()}
        Error: ${workflow.errorMessage}
        Error report: ${workflow.errorReport ?: 'No detailed report available'}
        """

        println "Error details written to: ${error_file}"
    }
}
```

Potete usare più handler insieme nel vostro script di flusso di lavoro:

```groovy title="main.nf - Combined handlers"
workflow {
    // ... il vostro codice del flusso di lavoro ...

    workflow.onError = {
        println "Workflow failed: ${workflow.errorMessage}"
    }

    workflow.onComplete = {
        def duration_mins = workflow.duration.toMinutes().round(2)
        def status = workflow.success ? "SUCCESS ✅" : "FAILED ❌"

        println """
        Pipeline finished: ${status}
        Duration: ${duration_mins} minutes
        """
    }
}
```

### Takeaway

In questa sezione, avete imparato:

- **Closure handler di eventi**: Closure nel vostro script di flusso di lavoro che vengono eseguite in diversi punti del ciclo di vita
- **Handler `onComplete`**: Per riepiloghi di esecuzione e reporting dei risultati
- **Handler `onError`**: Per gestione degli errori e logging dei fallimenti
- **Proprietà dell'oggetto workflow**: Accesso a `workflow.success`, `workflow.duration`, `workflow.errorMessage`, ecc.

Gli handler di eventi mostrano come potete usare la piena potenza del linguaggio Nextflow all'interno dei vostri script di flusso di lavoro per aggiungere capacità sofisticate di logging e notifica.

---

## Riepilogo

Congratulazioni, ce l'avete fatta!

In questa side quest, avete costruito una pipeline completa di elaborazione campioni che si è evoluta dalla gestione di base dei metadati a un flusso di lavoro sofisticato e pronto per la produzione.
Ogni sezione si è basata sulla precedente, dimostrando come i costrutti di programmazione trasformano semplici flussi di lavoro in potenti sistemi di elaborazione dati, con i seguenti benefici:

- **Codice più chiaro**: Comprendere dataflow vs scripting vi aiuta a scrivere flussi di lavoro più organizzati
- **Gestione robusta**: Navigazione sicura e operatori Elvis rendono i flussi di lavoro resilienti a dati mancanti
- **Elaborazione flessibile**: La logica condizionale permette ai vostri flussi di lavoro di elaborare appropriatamente diversi tipi di campioni
- **Risorse adattive**: Le direttive dinamiche ottimizzano l'uso delle risorse in base alle caratteristiche dell'input

Questa progressione rispecchia l'evoluzione nel mondo reale delle pipeline bioinformatiche, da prototipi di ricerca che gestiscono pochi campioni a sistemi di produzione che elaborano migliaia di campioni attraverso laboratori e istituzioni.
Ogni sfida che avete risolto e pattern che avete imparato riflette problemi reali che gli sviluppatori affrontano quando scalano i flussi di lavoro Nextflow.

Applicare questi pattern nel vostro lavoro vi permetterà di costruire flussi di lavoro robusti e pronti per la produzione.

### Pattern chiave

1.  **Dataflow vs Scripting:** Avete imparato a distinguere tra operazioni di dataflow (orchestrazione di canali) e scripting (codice che manipola dati), incluse le differenze cruciali tra operazioni su tipi diversi come `collect` su Channel vs List.

    - Dataflow: orchestrazione di canali

    ```groovy
    channel.fromPath('*.fastq').splitCsv(header: true)
    ```

    - Scripting: elaborazione dati su collezioni

    ```groovy
    sample_data.collect { it.toUpperCase() }
    ```

2.  **Elaborazione Avanzata di Stringhe**: Avete padroneggiato le espressioni regolari per analizzare nomi di file, la generazione dinamica di script nei processi e l'interpolazione di variabili (Nextflow vs Bash vs Shell).

    - Pattern matching

    ```groovy
    filename =~ ~/^(\w+)_(\w+)_(\d+)\.fastq$/
    ```

    - Funzione con return condizionale

    ```groovy
    def parseSample(filename) {
        def matcher = filename =~ pattern
        return matcher ? [valid: true, data: matcher[0]] : [valid: false]
    }
    ```

    - Collezione di file ad argomenti di comando (nel blocco script del processo)

    ```groovy
    script:
    def file_args = input_files.collect { file -> "--input ${file}" }.join(' ')
    """
    analysis_tool ${file_args} --output results.txt
    """
    ```

3.  **Creare Funzioni Riutilizzabili**: Avete imparato a estrarre logica complessa in funzioni con nome che possono essere chiamate dagli operatori di canale, rendendo i flussi di lavoro più leggibili e manutenibili.

    - Definire una funzione con nome

    ```groovy
    def separateMetadata(row) {
        def sample_meta = [ /* code hidden for brevity */ ]
        def fastq_path = file(row.file_path)
        def m = (fastq_path.name =~ /^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$/)
        def file_meta = m ? [ /* code hidden for brevity */ ] : [:]
        def priority = sample_meta.quality > 40 ? 'high' : 'normal'

        return tuple(sample_meta + file_meta + [priority: priority], fastq_path)
    }
    ```

    - Chiamare la funzione con nome in un flusso di lavoro

    ```groovy
    workflow {
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map{ row -> separateMetadata(row) }

        ch_fastp = FASTP(ch_samples)
    }
    ```

4.  **Direttive di Risorse Dinamiche con Closure**: Avete esplorato l'uso di closure nelle direttive dei processi per l'allocazione adattiva delle risorse in base alle caratteristiche dell'input.

    - Closure con nome e composizione

    ```groovy
    def enrichData = normalizeId >> addQualityCategory >> addFlags
    def processor = generalFunction.curry(fixedParam)
    ```

    - Closure con accesso allo scope

    ```groovy
    def collectStats = { data -> stats.count++; return data }
    ```

5.  **Logica Condizionale e Controllo dei Processi**: Avete aggiunto routing intelligente usando gli operatori `.branch()` e `.filter()`, sfruttando la truthiness per espressioni condizionali concise.

    - Usare `.branch()` per instradare dati attraverso diversi rami del flusso di lavoro

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
    if (sample.files) println "Has files"
    ```

    - Usare `filter()` per creare sottoinsiemi di dati con 'truthiness'

    ```groovy
    ch_valid_samples = ch_samples
        .filter { meta, reads ->
            meta.id && meta.organism && meta.depth >= 25000000
        }
    ```

6.  **Operatori di Navigazione Sicura ed Elvis**: Avete reso la pipeline robusta contro dati mancanti usando `?.` per accesso null-safe alle proprietà e `?:` per fornire valori predefiniti.

    ```groovy
    def id = data?.sample?.id ?: 'unknown'
    ```

7.  **Validazione con error() e log.warn**: Avete imparato a validare gli input precocemente e fallire velocemente con messaggi di errore chiari.

    ```groovy
    try {
        def errors = validateSample(sample)
        if (errors) throw new RuntimeException("Invalid: ${errors.join(', ')}")
    } catch (Exception e) {
        println "Error: ${e.message}"
    }
    ```

8.  **Handler di Eventi di Configurazione**: Avete imparato a usare gli handler di eventi del flusso di lavoro (`onComplete` e `onError`) per logging, notifiche e gestione del ciclo di vita.

    - Usare `onComplete` per logging e notifiche

    ```groovy
    workflow.onComplete = {
        println "Success     : ${workflow.success}"
        println "exit status : ${workflow.exitStatus}"

        if (workflow.success) {
            println "✅ Pipeline completed successfully!"
        } else {
            println "❌ Pipeline failed!"
            println "Error: ${workflow.errorMessage}"
        }
    }
    ```

    - Usare `onError` per intraprendere azioni specificamente in caso di fallimento

    ```groovy
    workflow.onError = {
        // Scrivi log di errore dettagliato
        def error_file = file("${workflow.launchDir}/error.log")
        error_file.text = """
        Time: ${new Date()}
        Error: ${workflow.errorMessage}
        Error report: ${workflow.errorReport ?: 'No detailed report available'}
        """

        println "Error details written to: ${error_file}"
    }
    ```

### Risorse aggiuntive

- [Riferimento Linguaggio Nextflow](https://nextflow.io/docs/latest/reference/syntax.html)
- [Operatori Nextflow](https://www.nextflow.io/docs/latest/operator.html)
- [Sintassi Script Nextflow](https://www.nextflow.io/docs/latest/script.html)
- [Libreria Standard Nextflow](https://nextflow.io/docs/latest/reference/stdlib.html)

Assicuratevi di consultare queste risorse quando dovete esplorare funzionalità più avanzate.

Trarrete beneficio dalla pratica e dall'espansione delle vostre competenze per:

- Scrivere flussi di lavoro più puliti con una corretta separazione tra dataflow e scripting
- Padroneggiare l'interpolazione di variabili per evitare insidie comuni con variabili Nextflow, Bash e shell
- Usare direttive di risorse dinamiche per flussi di lavoro efficienti e adattivi
- Trasformare collezioni di file in argomenti da linea di comando formattati correttamente
- Gestire con eleganza diverse convenzioni di denominazione dei file e formati di input usando regex ed elaborazione di stringhe
- Costruire codice riutilizzabile e manutenibile usando pattern avanzati di closure e programmazione funzionale
- Elaborare e organizzare dataset complessi usando operazioni su collezioni
- Aggiungere validazione, gestione degli errori e logging per rendere i vostri flussi di lavoro pronti per la produzione
- Implementare gestione del ciclo di vita del flusso di lavoro con handler di eventi

---

## Cosa c'è dopo?

Tornate al [menu delle Side Quest](./index.md) o cliccate il pulsante in basso a destra della pagina per passare al prossimo argomento nella lista.
