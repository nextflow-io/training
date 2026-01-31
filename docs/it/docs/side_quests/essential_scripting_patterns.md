# Pattern di Scripting Essenziali in Nextflow

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduzione assistita da IA - [scopri di più e suggerisci miglioramenti](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Nextflow è un linguaggio di programmazione che viene eseguito sulla Java Virtual Machine. Sebbene Nextflow sia costruito su [Groovy](http://groovy-lang.org/) e ne condivida gran parte della sintassi, Nextflow è più di un semplice "Groovy con estensioni" -- è un linguaggio autonomo con una [sintassi](https://nextflow.io/docs/latest/reference/syntax.html) completamente specificata e una [libreria standard](https://nextflow.io/docs/latest/reference/stdlib.html).

È possibile scrivere molto codice Nextflow senza andare oltre la sintassi di base per variabili, mappe e liste. La maggior parte dei tutorial su Nextflow si concentra sull'orchestrazione dei workflow (canali, processi e flusso di dati), e si può andare sorprendentemente lontano con solo questo.

Tuttavia, quando è necessario manipolare dati, analizzare nomi di file complessi, implementare logica condizionale o costruire workflow di produzione robusti, è utile pensare a due aspetti distinti del codice: **dataflow** (canali, operatori, processi e workflow) e **scripting** (il codice all'interno di closure, funzioni e script di processo). Sebbene questa distinzione sia in qualche modo arbitraria—è tutto codice Nextflow—fornisce un modello mentale utile per comprendere quando si sta orchestrando la pipeline rispetto a quando si stanno manipolando i dati. Padroneggiare entrambi migliora notevolmente la capacità di scrivere workflow chiari e manutenibili.

### Obiettivi di apprendimento

Questa missione secondaria ti accompagna in un percorso pratico dai concetti di base ai pattern pronti per la produzione.
Trasformeremo un semplice workflow di lettura CSV in una sofisticata pipeline bioinformatica, evolvendola passo dopo passo attraverso sfide realistiche:

- **Comprendere i confini:** Distinguere tra operazioni di dataflow e scripting, e capire come lavorano insieme
- **Manipolazione dei dati:** Estrarre, trasformare e creare sottoinsiemi di mappe e collezioni utilizzando operatori potenti
- **Elaborazione delle stringhe:** Analizzare schemi di denominazione dei file complessi con pattern regex e padroneggiare l'interpolazione delle variabili
- **Funzioni riutilizzabili:** Estrarre logica complessa in funzioni nominate per workflow più puliti e manutenibili
- **Logica dinamica:** Costruire processi che si adattano a diversi tipi di input e utilizzare closure per l'allocazione dinamica delle risorse
- **Routing condizionale:** Instradare intelligentemente i campioni attraverso diversi processi in base alle loro caratteristiche di metadati
- **Operazioni sicure:** Gestire in modo elegante i dati mancanti con operatori null-safe e validare gli input con messaggi di errore chiari
- **Handler basati sulla configurazione:** Utilizzare handler di eventi del workflow per logging, notifiche e gestione del ciclo di vita

### Prerequisiti

Prima di intraprendere questa missione secondaria, dovreste:

- Aver completato il tutorial [Hello Nextflow](../hello_nextflow/README.md) o un corso equivalente per principianti.
- Essere a proprio agio nell'uso di concetti e meccanismi di base di Nextflow (processi, canali, operatori, lavoro con file, metadati)
- Avere familiarità di base con costrutti di programmazione comuni (variabili, mappe, liste)

Questo tutorial spiegherà i concetti di programmazione man mano che li incontriamo, quindi non è necessaria un'esperienza di programmazione approfondita.
Inizieremo con concetti fondamentali e costruiremo fino a pattern avanzati.

---

## 0. Iniziare

#### Aprire il codespace di formazione

Se non lo ha ancora fatto, assicuratevi di aprire l'ambiente di formazione come descritto in [Configurazione Ambiente](../envsetup/index.md).

[![Apri in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Spostarsi nella directory del progetto

Spostiamoci nella directory dove si trovano i file per questo tutorial.

```bash
cd side-quests/essential_scripting_patterns
```

#### Rivedere i materiali

Troverà un file workflow principale e una directory `data` contenente file di dati di esempio.

```console title="Contenuto della directory"
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

Il nostro CSV di esempio contiene informazioni sui campioni biologici che necessitano di elaborazioni diverse in base alle loro caratteristiche:

```console title="samples.csv"
sample_id,organism,tissue_type,sequencing_depth,file_path,quality_score
SAMPLE_001,human,liver,30000000,data/sequences/SAMPLE_001_S1_L001_R1_001.fastq,38.5
SAMPLE_002,mouse,brain,25000000,data/sequences/SAMPLE_002_S2_L001_R1_001.fastq,35.2
SAMPLE_003,human,kidney,45000000,data/sequences/SAMPLE_003_S3_L001_R1_001.fastq,42.1
```

Utilizzeremo questo dataset realistico per esplorare tecniche di programmazione pratiche che incontrerà in workflow bioinformatici reali.

<!-- TODO: Can we make this more domain-agnostic? -->

<!-- TODO: add an assignment statement? #### Review the assignment -->

#### Checklist di preparazione

Pensa di essere pronto per iniziare?

- [ ] Comprendo l'obiettivo di questo corso e i suoi prerequisiti
- [ ] Il mio codespace è attivo e funzionante
- [ ] Ho impostato la mia directory di lavoro in modo appropriato
<!-- - [ ] I understand the assignment -->

Se può spuntare tutte le caselle, è pronto per partire.

---

## 1. Dataflow vs Scripting: Comprendere i Confini

### 1.1. Identificare Cosa è Cosa

Quando si scrivono workflow Nextflow, è importante distinguere tra **dataflow** (come i dati si muovono attraverso canali e processi) e **scripting** (il codice che manipola i dati e prende decisioni). Costruiamo un workflow che dimostri come lavorano insieme.

#### 1.1.1. Workflow Nextflow di Base

Iniziamo con un semplice workflow che legge solo il file CSV (lo abbiamo già fatto per voi in `main.nf`):

```groovy title="main.nf" linenums="1"
workflow {
    ch_samples = channel.fromPath("./data/samples.csv")
        .splitCsv(header: true)
        .view()
}
```

Il blocco `workflow` definisce la struttura della nostra pipeline, mentre `channel.fromPath()` crea un canale da un percorso di file. L'operatore `.splitCsv()` elabora il file CSV e converte ogni riga in una struttura dati di tipo mappa.

Esegua questo workflow per vedere i dati CSV grezzi:

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

Ora aggiungeremo scripting per trasformare i dati, utilizzando l'operatore `.map()` che probabilmente conosce già. Questo operatore prende una 'closure' dove possiamo scrivere codice per trasformare ogni elemento.

!!! note

    Una **closure** è un blocco di codice che può essere passato in giro ed eseguito successivamente. Pensate ad essa come a una funzione definita inline. Le closure sono scritte con parentesi graffe `{ }` e possono prendere parametri. Sono fondamentali per come funzionano gli operatori Nextflow e se avete scritto Nextflow per un po', potrebbe averle già utilizzate senza rendersene conto!

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

Questa è la nostra prima **closure** - una funzione anonima che si può passare come argomento (simile alle lambda in Python o alle arrow function in JavaScript). Le closure sono essenziali per lavorare con gli operatori Nextflow.

La closure `{ row -> return row }` prende un parametro `row` (potrebbe essere qualsiasi nome: `item`, `sample`, ecc.).

Quando l'operatore `.map()` elabora ogni elemento del canale, passa quell'elemento alla closure. Qui, `row` contiene una riga CSV alla volta.

Applicate questa modifica ed eseguite il workflow:

```bash
nextflow run main.nf
```

Vedrà lo stesso output di prima, perché stiamo semplicemente restituendo l'input invariato. Questo conferma che l'operatore map funziona correttamente. Ora iniziamo a trasformare i dati.

#### 1.1.3. Creare una Struttura Dati Map

Ora scriveremo logica di **scripting** all'interno della nostra closure per trasformare ogni riga di dati. Questo è dove elaboriamo i singoli elementi di dati piuttosto che orchestrare il flusso di dati.

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

La mappa `sample_meta` è una struttura dati chiave-valore (come i dizionari in Python, gli oggetti in JavaScript o gli hash in Ruby) che memorizza informazioni correlate: ID campione, organismo, tipo di tessuto, profondità di sequenziamento e punteggio di qualità.

Utilizziamo metodi di manipolazione delle stringhe come `.toLowerCase()` e `.replaceAll()` per pulire i nostri dati, e metodi di conversione di tipo come `.toInteger()` e `.toDouble()` per convertire i dati stringa dal CSV nei tipi numerici appropriati.

Applicate questa modifica ed eseguite il workflow:

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

Ora aggiungiamo più scripting - questa volta utilizzando un operatore ternario per prendere decisioni in base ai valori dei dati.

Effettui la seguente modifica:

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

!!! Note

    Non modifichi mai le mappe passate nelle closure - crei sempre nuove mappe utilizzando `+` (per esempio). In Nextflow, gli stessi dati spesso fluiscono attraverso più operazioni simultaneamente. Modificare una mappa sul posto può causare effetti collaterali imprevedibili quando altre operazioni fanno riferimento a quello stesso oggetto. Creare nuove mappe assicura che ogni operazione abbia la propria copia pulita.

Esegua il workflow modificato:

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

Mentre l'operatore `+` aggiunge chiavi a una mappa, a volte è necessario fare il contrario - estrarre solo chiavi specifiche. Il metodo `.subMap()` è perfetto per questo.

Aggiungiamo una riga per creare una versione semplificata dei nostri metadati che contiene solo i campi di identificazione:

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
                println "Solo campi ID: ${id_only}"

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

Esegua il workflow modificato:

```bash
nextflow run main.nf
```

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [peaceful_cori] DSL2 - revision: 4cc4a8340f

    Solo campi ID: [id:sample_001, organism:human, tissue:liver]
    Solo campi ID: [id:sample_002, organism:mouse, tissue:brain]
    Solo campi ID: [id:sample_003, organism:human, tissue:kidney]
    [id:sample_001, organism:human, tissue:liver, depth:30000000, quality:38.5, priority:normal]
    [id:sample_002, organism:mouse, tissue:brain, depth:25000000, quality:35.2, priority:normal]
    [id:sample_003, organism:human, tissue:kidney, depth:45000000, quality:42.1, priority:high]
    ```

Questo mostra sia i metadati completi visualizzati dall'operazione `view()` che il sottoinsieme estratto che abbiamo stampato con `println`.

Il metodo `.subMap()` prende una lista di chiavi e restituisce una nuova mappa contenente solo quelle chiavi. Se una chiave non esiste nella mappa originale, semplicemente non viene inclusa nel risultato.

Questo è particolarmente utile quando è necessario creare diverse versioni di metadati per diversi processi - alcuni potrebbero aver bisogno di metadati completi mentre altri necessitano solo di campi di identificazione minimali.

Ora rimuova quelle istruzioni println per ripristinare il workflow allo stato precedente, poiché non ne abbiamo bisogno andando avanti.

!!! tip "Riepilogo Operazioni su Mappe"

    - **Aggiungere chiavi**: `map1 + [new_key: value]` - Crea nuova mappa con chiavi aggiuntive
    - **Estrarre chiavi**: `map1.subMap(['key1', 'key2'])` - Crea nuova mappa con solo le chiavi specificate
    - **Entrambe le operazioni creano nuove mappe** - Le mappe originali rimangono invariate

#### 1.1.6. Combinare Mappe e Restituire Risultati

Finora, abbiamo restituito solo quella che la comunità Nextflow chiama 'meta map', e abbiamo ignorato i file a cui quei metadati si riferiscono. Ma se state scrivendo workflow Nextflow, probabilmente vuole fare qualcosa con quei file.

Produciamo una struttura di canale comprendente una tupla di 2 elementi: la mappa di metadati arricchita e il percorso di file corrispondente. Questo è un pattern comune in Nextflow per passare dati ai processi.

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

Applicate questa modifica ed eseguite il workflow:

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

!!! note

    **Mappe e Metadati**: Le mappe sono fondamentali per lavorare con i metadati in Nextflow. Per una spiegazione più dettagliata del lavoro con mappe di metadati, veda la missione secondaria [Lavorare con i metadati](./metadata.md).

Il nostro workflow dimostra il pattern centrale: **operazioni di dataflow** (`workflow`, `channel.fromPath()`, `.splitCsv()`, `.map()`, `.view()`) orchestrano come i dati si muovono attraverso la pipeline, mentre **scripting** (mappe `[key: value]`, metodi di stringhe, conversioni di tipo, operatori ternari) all'interno della closure `.map()` gestisce la trasformazione dei singoli elementi di dati.

### 1.2. Comprendere Diversi Tipi: Channel vs List

Finora, tutto bene, possiamo distinguere tra operazioni di dataflow e scripting. Ma che dire quando lo stesso nome di metodo esiste in entrambi i contesti?

Un esempio perfetto è il metodo `collect`, che esiste sia per i tipi di canale che per i tipi List nella libreria standard Nextflow. Il metodo `collect()` su una List trasforma ogni elemento, mentre l'operatore `collect()` su un canale raccoglie tutte le emissioni del canale in un canale a singolo elemento.

Dimostriamo questo con alcuni dati di esempio, iniziando col rinfrescarci su cosa fa l'operatore `collect()` del canale. Dia un'occhiata a `collect.nf`:

```groovy title="collect.nf" linenums="1"
def sample_ids = ['sample_001', 'sample_002', 'sample_003']

// channel.collect() - raggruppa più emissioni di canale in una
ch_input = channel.fromList(sample_ids)
ch_input.view { sample -> "Elemento individuale del canale: ${sample}" }
ch_collected = ch_input.collect()
ch_collected.view { list -> "Risultato channel.collect(): ${list} (${list.size()} elementi raggruppati in 1)" }
```

Passaggi:

- Definire una List di ID campione
- Creare un canale con `fromList()` che emette ogni ID campione separatamente
- Stampare ogni elemento con `view()` mentre fluisce attraverso
- Raccogliere tutti gli elementi in una singola lista con l'operatore `collect()` del canale
- Stampare il risultato raccolto (singolo elemento contenente tutti gli ID campione) con un secondo `view()`

Abbiamo cambiato la struttura del canale, ma non abbiamo cambiato i dati stessi.

Esegua il workflow per confermarlo:

```bash
nextflow run collect.nf
```

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `collect.nf` [loving_mendel] DSL2 - revision: e8d054a46e

    Elemento individuale del canale: sample_001
    Elemento individuale del canale: sample_002
    Elemento individuale del canale: sample_003
    Risultato channel.collect(): [sample_001, sample_002, sample_003] (3 elementi raggruppati in 1)
    ```

`view()` restituisce un output per ogni emissione del canale, quindi sappiamo che questo singolo output contiene tutti e 3 gli elementi originali raggruppati in una lista.

Ora vediamo il metodo `collect` su una List in azione. Modifichi `collect.nf` per applicare il metodo `collect` della List alla lista originale di ID campione:

=== "Dopo"

    ```groovy title="main.nf" linenums="1" hl_lines="9-13"
    def sample_ids = ['sample_001', 'sample_002', 'sample_003']

    // channel.collect() - raggruppa più emissioni di canale in una
    ch_input = channel.fromList(sample_ids)
    ch_input.view { sample -> "Elemento individuale del canale: ${sample}" }
    ch_collected = ch_input.collect()
    ch_collected.view { list -> "Risultato channel.collect(): ${list} (${list.size()} elementi raggruppati in 1)" }

    // List.collect() - trasforma ogni elemento, preserva la struttura
    def formatted_ids = sample_ids.collect { id ->
        id.toUpperCase().replace('SAMPLE_', 'SPECIMEN_')
    }
    println "Risultato List.collect(): ${formatted_ids} (${sample_ids.size()} elementi trasformati in ${formatted_ids.size()})"
    ```

=== "Prima"

    ```groovy title="main.nf" linenums="1"
    def sample_ids = ['sample_001', 'sample_002', 'sample_003']

    // channel.collect() - raggruppa più emissioni di canale in una
    ch_input = channel.fromList(sample_ids)
    ch_input.view { sample -> "Elemento individuale del canale: ${sample}" }
    ch_collected = ch_input.collect()
    ch_collected.view { list -> "Risultato channel.collect(): ${list} (${list.size()} elementi raggruppati in 1)" }
    ```

In questo nuovo snippet:

- Definiamo una nuova variabile `formatted_ids` che utilizza il metodo `collect` della List per trasformare ogni ID campione nella lista originale
- Stampiamo il risultato usando `println`

Esegua il workflow modificato:

```bash
nextflow run collect.nf
```

??? success "Output del comando"

    ```console hl_lines="5"
    N E X T F L O W   ~  version 25.10.2

    Launching `collect.nf` [cheeky_stonebraker] DSL2 - revision: 2d5039fb47

    Risultato List.collect(): [SPECIMEN_001, SPECIMEN_002, SPECIMEN_003] (3 elementi trasformati in 3)
    Elemento individuale del canale: sample_001
    Elemento individuale del canale: sample_002
    Elemento individuale del canale: sample_003
    Risultato channel.collect(): [sample_001, sample_002, sample_003] (3 elementi raggruppati in 1)
    ```

Questa volta, NON abbiamo cambiato la struttura dei dati, abbiamo ancora 3 elementi nella lista, ma ABBIAMO trasformato ogni elemento usando il metodo `collect` della List per produrre una nuova lista con valori modificati. Questo è simile all'uso dell'operatore `map` su un canale, ma sta operando su una struttura dati List piuttosto che su un canale.

`collect` è un caso estremo che stiamo usando qui per fare un punto. La lezione chiave è che quando si scrivono workflow, si deve sempre distinguere tra **strutture dati** (List, Map, ecc.) e **canali** (costrutti di dataflow). Le operazioni possono condividere nomi ma comportarsi in modo completamente diverso a seconda del tipo su cui vengono chiamate.

### 1.3. L'Operatore Spread (`*.`) - Scorciatoia per l'Estrazione di Proprietà

Correlato al metodo `collect` della List è l'operatore spread (`*.`), che fornisce un modo conciso per estrarre proprietà dalle collezioni. È essenzialmente zucchero sintattico per un pattern `collect` comune.

Aggiungiamo una dimostrazione al nostro file `collect.nf`:

=== "Dopo"

    ```groovy title="collect.nf" linenums="1" hl_lines="15-18"
    def sample_ids = ['sample_001', 'sample_002', 'sample_003']

    // channel.collect() - raggruppa più emissioni di canale in una
    ch_input = channel.fromList(sample_ids)
    ch_input.view { sample -> "Elemento individuale del canale: ${sample}" }
    ch_collected = ch_input.collect()
    ch_collected.view { list -> "Risultato channel.collect(): ${list} (${list.size()} elementi raggruppati in 1)" }

    // List.collect() - trasforma ogni elemento, preserva la struttura
    def formatted_ids = sample_ids.collect { id ->
        id.toUpperCase().replace('SAMPLE_', 'SPECIMEN_')
    }
    println "Risultato List.collect(): ${formatted_ids} (${sample_ids.size()} elementi trasformati in ${formatted_ids.size()})"

    // Operatore spread - accesso conciso alle proprietà
    def sample_data = [[id: 's1', quality: 38.5], [id: 's2', quality: 42.1], [id: 's3', quality: 35.2]]
    def all_ids = sample_data*.id
    println "Risultato operatore spread: ${all_ids}"
    ```

=== "Prima"

    ```groovy title="collect.nf" linenums="1"
    def sample_ids = ['sample_001', 'sample_002', 'sample_003']

    // channel.collect() - raggruppa più emissioni di canale in una
    ch_input = channel.fromList(sample_ids)
    ch_input.view { sample -> "Elemento individuale del canale: ${sample}" }
    ch_collected = ch_input.collect()
    ch_collected.view { list -> "Risultato channel.collect(): ${list} (${list.size()} elementi raggruppati in 1)" }

    // List.collect() - trasforma ogni elemento, preserva la struttura
    def formatted_ids = sample_ids.collect { id ->
        id.toUpperCase().replace('SAMPLE_', 'SPECIMEN_')
    }
    println "Risultato List.collect(): ${formatted_ids} (${sample_ids.size()} elementi trasformati in ${formatted_ids.size()})"
    ```

Esegua il workflow aggiornato:

```bash title="Testare operatore spread"
nextflow run collect.nf
```

??? success "Output del comando"

    ```console hl_lines="6"
    N E X T F L O W   ~  version 25.10.2

    Launching `collect.nf` [cranky_galileo] DSL2 - revision: 5f3c8b2a91

    Risultato List.collect(): [SPECIMEN_001, SPECIMEN_002, SPECIMEN_003] (3 elementi trasformati in 3)
    Risultato operatore spread: [s1, s2, s3]
    Elemento individuale del canale: sample_001
    Elemento individuale del canale: sample_002
    Elemento individuale del canale: sample_003
    Risultato channel.collect(): [sample_001, sample_002, sample_003] (3 elementi raggruppati in 1)
    ```

L'operatore spread `*.` è una scorciatoia per un pattern collect comune:

```groovy
// Questi sono equivalenti:
def ids = samples*.id
def ids = samples.collect { it.id }

// Funziona anche con chiamate a metodi:
def names = files*.getName()
def names = files.collect { it.getName() }
```

L'operatore spread è particolarmente utile quando è necessario estrarre una singola proprietà da una lista di oggetti - è più leggibile che scrivere la closure `collect` completa.

!!! tip "Quando Usare Spread vs Collect"

    - **Usare spread (`*.`)** per accesso semplice alle proprietà: `samples*.id`, `files*.name`
    - **Usare collect** per trasformazioni o logica complessa: `samples.collect { it.id.toUpperCase() }`, `samples.collect { [it.id, it.quality > 40] }`

### Takeaway

In questa sezione, ha imparato:

- **Dataflow vs scripting**: Gli operatori di canale orchestrano come i dati fluiscono attraverso la pipeline, mentre lo scripting trasforma i singoli elementi di dati
- **Comprendere i tipi**: Lo stesso nome di metodo (come `collect`) può comportarsi in modo diverso a seconda del tipo su cui viene chiamato (Channel vs List)
- **Il contesto conta**: Sia sempre consapevole se state lavorando con canali (dataflow) o strutture dati (scripting)

Comprendere questi confini è essenziale per il debug, la documentazione e la scrittura di workflow manutenibili.

Successivamente approfondiremo le capacità di elaborazione delle stringhe, che sono essenziali per gestire dati del mondo reale.

---

## 2. Elaborazione di Stringhe e Generazione Dinamica di Script

Padroneggiare l'elaborazione delle stringhe separa i workflow fragili dalle pipeline robuste. Questa sezione copre l'analisi di nomi di file complessi, la generazione dinamica di script e l'interpolazione di variabili.

### 2.1. Pattern Matching ed Espressioni Regolari

I file bioinformatici hanno spesso convenzioni di denominazione complesse che codificano metadati. Estraiamoli automaticamente usando il pattern matching con espressioni regolari.

Torneremo al nostro workflow `main.nf` e aggiungeremo un po' di logica di pattern matching per estrarre informazioni aggiuntive sui campioni dai nomi dei file. I file FASTQ nel nostro dataset seguono convenzioni di denominazione in stile Illumina con nomi come `SAMPLE_001_S1_L001_R1_001.fastq.gz`. Questi potrebbero sembrare criptici, ma in realtà codificano metadati utili come ID campione, numero di lane e direzione di lettura. Useremo le capacità regex per analizzare questi nomi.

Effettui la seguente modifica al workflow `main.nf` esistente:

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

1. **Letterali di espressioni regolari** usando la sintassi `~/pattern/` - questo crea un pattern regex senza bisogno di escape dei backslash
2. **Pattern matching** con l'operatore `=~` - questo tenta di abbinare una stringa a un pattern regex
3. **Oggetti matcher** che catturano gruppi con `[0][1]`, `[0][2]`, ecc. - `[0]` si riferisce all'intera corrispondenza, `[1]`, `[2]`, ecc. si riferiscono ai gruppi catturati tra parentesi

Analizziamo il pattern regex `^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$`:

| Pattern             | Corrisponde                            | Cattura                             |
| ------------------- | -------------------------------------- | ----------------------------------- |
| `^(.+)`             | Nome campione dall'inizio              | Gruppo 1: nome campione             |
| `_S(\d+)`           | Numero campione `_S1`, `_S2`, ecc.     | Gruppo 2: numero campione           |
| `_L(\d{3})`         | Numero lane `_L001`                    | Gruppo 3: lane (3 cifre)            |
| `_(R[12])`          | Direzione lettura `_R1` o `_R2`        | Gruppo 4: direzione lettura         |
| `_(\d{3})`          | Numero chunk `_001`                    | Gruppo 5: chunk (3 cifre)           |
| `\.fastq(?:\.gz)?$` | Estensione file `.fastq` o `.fastq.gz` | Non catturato (?: è non-catturante) |

Questo analizza le convenzioni di denominazione in stile Illumina per estrarre automaticamente i metadati.

Esegua il workflow modificato:

```bash title="Testare pattern matching"
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

I blocchi script dei processi sono essenzialmente stringhe multi-linea che vengono passate alla shell. Si potete usare **logica condizionale** (if/else, operatori ternari) per generare dinamicamente diverse stringhe di script in base alle caratteristiche dell'input. Questo è essenziale per gestire diversi tipi di input—come letture single-end vs paired-end—senza duplicare le definizioni dei processi.

Aggiungiamo un processo al nostro workflow che dimostri questo pattern. Apra `modules/fastp.nf` e dia un'occhiata:

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

Il processo prende file FASTQ come input ed esegue lo strumento `fastp` per tagliare adapter e filtrare letture di bassa qualità. Sfortunatamente, la persona che ha scritto questo processo non ha previsto le letture single-end che abbiamo nel nostro dataset di esempio. Aggiungiamolo al nostro workflow e vediamo cosa succede:

Per prima cosa, includa il modulo alla prima riga del workflow `main.nf`:

```groovy title="main.nf" linenums="1"
include { FASTP } from './modules/fastp.nf'
```

Quindi modifichi il blocco `workflow` per collegare il canale `ch_samples` al processo `FASTP`:

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

Esegua questo workflow modificato:

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

Si potete vedere che il processo sta cercando di eseguire `fastp` con un valore `null` per il secondo file di input, il che lo fa fallire. Questo perché il nostro dataset contiene letture single-end, ma il processo è codificato per aspettarsi letture paired-end (due file di input alla volta).

Corregga questo aggiungendo logica condizionale al blocco `script:` del processo `FASTP`. Un'istruzione if/else controlla il conteggio dei file di lettura e regola il comando di conseguenza.

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

Ora il workflow può gestire con eleganza sia letture single-end che paired-end. La logica condizionale controlla il numero di file di input e costruisce il comando appropriato per `fastp`. Vediamo se funziona:

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

Sembra buono! Se controlliamo i comandi effettivi che sono stati eseguiti (personalizzi per il suo hash di attività):

```console title="Controllare comandi eseguiti"
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

Un altro uso comune della logica dinamica degli script può essere visto in [il modulo Genomics di Nextflow for Science](../../nf4science/genomics/02_joint_calling). In quel modulo, il processo GATK chiamato può prendere più file di input, ma ciascuno deve essere prefissato con `-V` per formare una riga di comando corretta. Il processo usa scripting per trasformare una collezione di file di input (`all_gvcfs`) negli argomenti di comando corretti:

```groovy title="manipolazione riga di comando per GATK" linenums="1" hl_lines="2 5"
    script:
    def gvcfs_line = all_gvcfs.collect { gvcf -> "-V ${gvcf}" }.join(' ')
    """
    gatk GenomicsDBImport \
        ${gvcfs_line} \
        -L ${interval_list} \
        --genomicsdb-workspace-path ${cohort_name}_gdb
    """
```

Questi pattern di utilizzo dello scripting nei blocchi script dei processi sono estremamente potenti e possono essere applicati in molti scenari - dalla gestione di tipi di input variabili alla costruzione di argomenti di riga di comando complessi da collezioni di file, rendendo i processi veramente adattabili ai requisiti diversificati dei dati del mondo reale.

### 2.3. Interpolazione di Variabili: Variabili Nextflow e Shell

Gli script dei processi mescolano variabili Nextflow, variabili shell e sostituzioni di comandi, ciascuna con sintassi di interpolazione diversa. Usare la sintassi sbagliata causa errori. Esploriamo questi con un processo che crea un report di elaborazione.

Dia un'occhiata al file modulo `modules/generate_report.nf`:

```groovy title="modules/generate_report.nf" linenums="1"
process GENERATE_REPORT {

    publishDir 'results/reports', mode: 'copy'

    input:
    tuple val(meta), path(reads)

    output:
    path "${meta.id}_report.txt"

    script:
    """
    echo "Elaborazione ${reads}" > ${meta.id}_report.txt
    echo "Campione: ${meta.id}" >> ${meta.id}_report.txt
    """
}
```

Questo processo scrive un semplice report con l'ID campione e il nome del file. Ora eseguiamolo per vedere cosa succede quando dobbiamo mescolare diversi tipi di variabili.

Includa il processo nel suo `main.nf` e lo aggiunga al workflow:

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

Ora eseguite il workflow e controllate i report generati in `results/reports/`. Dovrebbero contenere informazioni di base su ogni campione.

<!-- TODO: add the run command -->

??? success "Output del comando"

    ```console
    <!-- TODO: output -->
    ```

Ma cosa succederebbe se volessimo aggiungere informazioni su quando e dove è avvenuta l'elaborazione? Modifichiamo il processo per utilizzare variabili **shell** e un po' di sostituzione di comandi per includere l'utente corrente, l'hostname e la data nel report:

=== "Dopo"

    ```groovy title="modules/generate_report.nf" linenums="10" hl_lines="5-7"
        script:
        """
        echo "Elaborazione ${reads}" > ${meta.id}_report.txt
        echo "Campione: ${meta.id}" >> ${meta.id}_report.txt
        echo "Elaborato da: ${USER}" >> ${meta.id}_report.txt
        echo "Hostname: $(hostname)" >> ${meta.id}_report.txt
        echo "Data: $(date)" >> ${meta.id}_report.txt
        """
    ```

=== "Prima"

    ```groovy title="modules/generate_report.nf" linenums="10"
        script:
        """
        echo "Elaborazione ${reads}" > ${meta.id}_report.txt
        echo "Campione: ${meta.id}" >> ${meta.id}_report.txt
        """
    ```

Se esegue questo, noterà un errore - Nextflow cerca di interpretare `${USER}` come una variabile Nextflow che non esiste.

??? failure "Output del comando"

    ```console
    Error modules/generate_report.nf:15:27: `USER` is not defined
    │  15 |     echo "Elaborato da: ${USER}" >> ${meta.id}_report.txt
    ╰     |                           ^^^^

    ERROR ~ Script compilation failed
    ```

Dobbiamo fare l'escape così che Bash possa gestirlo invece.

Corregga questo facendo l'escape delle variabili shell e delle sostituzioni di comandi con un backslash (`\`):

=== "Dopo"

    ```groovy title="modules/generate_report.nf" linenums="10" hl_lines="5-7"
        script:
        """
        echo "Elaborazione ${reads}" > ${meta.id}_report.txt
        echo "Campione: ${meta.id}" >> ${meta.id}_report.txt
        echo "Elaborato da: \${USER}" >> ${meta.id}_report.txt
        echo "Hostname: \$(hostname)" >> ${meta.id}_report.txt
        echo "Data: \$(date)" >> ${meta.id}_report.txt
        """
    ```

=== "Prima"

    ```groovy title="modules/generate_report.nf" linenums="10"
        script:
        """
        echo "Elaborazione ${reads}" > ${meta.id}_report.txt
        echo "Campione: ${meta.id}" >> ${meta.id}_report.txt
        echo "Elaborato da: ${USER}" >> ${meta.id}_report.txt
        echo "Hostname: $(hostname)" >> ${meta.id}_report.txt
        echo "Data: $(date)" >> ${meta.id}_report.txt
        """
    ```

Ora funziona! Il backslash (`\`) dice a Nextflow "non interpretare questo, passalo a Bash."

### Takeaway

In questa sezione, ha imparato tecniche di **elaborazione delle stringhe**:

- **Espressioni regolari per l'analisi dei file**: Utilizzare l'operatore `=~` e pattern regex (`~/pattern/`) per estrarre metadati da convenzioni di denominazione dei file complesse
- **Generazione dinamica di script**: Utilizz
