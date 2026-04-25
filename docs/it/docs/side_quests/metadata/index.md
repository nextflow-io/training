# Metadati e meta map

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduzione assistita da IA - [scopri di più e suggerisci miglioramenti](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

In qualsiasi analisi scientifica, raramente lavoriamo solo con i file di dati grezzi.
Ogni file porta con sé informazioni aggiuntive: cosa contiene, da dove proviene e cosa lo rende speciale.
Queste informazioni extra sono ciò che chiamiamo metadati.

I metadati sono dati che descrivono altri dati.
I metadati tracciano dettagli importanti sui file e sulle condizioni sperimentali, e aiutano ad adattare le analisi alle caratteristiche uniche di ogni dataset.

Pensateci come a un catalogo di biblioteca: mentre i libri contengono il contenuto effettivo (dati grezzi), le schede catalografiche forniscono informazioni essenziali su ogni libro — quando è stato pubblicato, chi lo ha scritto, dove trovarlo (metadati).
Nei pipeline Nextflow, i metadati possono essere utilizzati per:

- Tracciare informazioni specifiche dei file durante il flusso di lavoro
- Configurare i processi in base alle caratteristiche dei file
- Raggruppare file correlati per analisi congiunte

### Obiettivi di apprendimento

In questa side quest, esploreremo come gestire i metadati nei flussi di lavoro.
Partendo da un semplice foglio dati (spesso chiamato samplesheet in bioinformatica) contenente informazioni di base sui file, imparerete come:

- Leggere e analizzare i metadati dei file da file CSV
- Creare e manipolare meta map
- Aggiungere nuovi campi di metadati durante l'esecuzione del flusso di lavoro
- Utilizzare i metadati per personalizzare il comportamento dei processi

Queste competenze vi aiuteranno a costruire pipeline più robusti e flessibili, in grado di gestire relazioni complesse tra file e requisiti di elaborazione.

### Prerequisiti

Prima di affrontare questa side quest, dovreste:

- Aver completato il tutorial [Hello Nextflow](../hello_nextflow/README.md) o un corso equivalente per principianti.
- Essere a proprio agio con i concetti e i meccanismi di base di Nextflow (processi, canali, operatori)

---

## 0. Iniziamo

#### Aprite il codespace di formazione

Se non l'avete ancora fatto, assicuratevi di aprire l'ambiente di formazione come descritto nella sezione [Configurazione dell'ambiente](../envsetup/index.md).

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Spostatevi nella directory del progetto

Spostiamoci nella directory dove si trovano i file per questo tutorial.

```bash
cd side-quests/metadata
```

Potete impostare VSCode per concentrarsi su questa directory:

```bash
code .
```

#### Esaminate i materiali

Troverete un file principale del flusso di lavoro e una directory `data` contenente un foglio dati e alcuni file di dati.

??? abstract "Contenuto della directory"

    ```console
    .
    ├── data
    │   ├── bonjour.txt
    │   ├── ciao.txt
    │   ├── guten_tag.txt
    │   ├── hallo.txt
    │   ├── hello.txt
    │   ├── hola.txt
    │   ├── salut.txt
    │   └── datasheet.csv
    ├── main.nf
    └── nextflow.config
    ```

Il flusso di lavoro nel file `main.nf` è uno stub che espanderete gradualmente in un flusso di lavoro completamente funzionante.

Il foglio dati elenca i percorsi dei file di dati e alcuni metadati associati, organizzati in 3 colonne:

- `id`: autoesplicativo, un ID assegnato al file
- `character`: il nome di un personaggio, che useremo in seguito per disegnare diverse creature
- `data`: percorsi ai file `.txt` che contengono saluti in diverse lingue

```console title="datasheet.csv"
id,character,recording
sampleA,squirrel,/workspaces/training/side-quests/metadata/data/bonjour.txt
sampleB,tux,/workspaces/training/side-quests/metadata/data/guten_tag.txt
sampleC,sheep,/workspaces/training/side-quests/metadata/data/hallo.txt
sampleD,turkey,/workspaces/training/side-quests/metadata/data/hello.txt
sampleE,stegosaurus,/workspaces/training/side-quests/metadata/data/hola.txt
sampleF,moose,/workspaces/training/side-quests/metadata/data/salut.txt
sampleG,turtle,/workspaces/training/side-quests/metadata/data/ciao.txt
```

Ogni file di dati contiene del testo di saluto in una delle cinque lingue (fr: francese, de: tedesco, es: spagnolo, it: italiano, en: inglese).

Vi forniremo anche uno strumento di analisi linguistica containerizzato chiamato `langid`.

#### Esaminate il compito

La vostra sfida è scrivere un flusso di lavoro Nextflow che:

1. **Identifichi** automaticamente la lingua in ogni file
2. **Raggruppi** i file per famiglia linguistica (lingue germaniche vs lingue romanze)
3. **Personalizzi** l'elaborazione per ogni file in base alla sua lingua e ai suoi metadati
4. **Organizzi** gli output per gruppo linguistico

Questo rappresenta un tipico schema di flusso di lavoro in cui i metadati specifici dei file guidano le decisioni di elaborazione; esattamente il tipo di problema che le meta map risolvono elegantemente.

#### Lista di controllo per la preparazione

Pensate di essere pronti a tuffarvi?

- [ ] Comprendo l'obiettivo di questo corso e i suoi prerequisiti
- [ ] Il mio codespace è attivo e funzionante
- [ ] Ho impostato correttamente la mia directory di lavoro
- [ ] Comprendo il compito assegnato

Se potete spuntare tutte le caselle, siete pronti a partire.

---

## 1. Caricare i metadati da un foglio dati

Aprite il file del flusso di lavoro `main.nf` per esaminare lo stub del flusso di lavoro che vi forniamo come punto di partenza.

```groovy title="main.nf" linenums="1"
#!/usr/bin/env nextflow

workflow  {

    ch_datasheet = channel.fromPath("./data/datasheet.csv")

}
```

Come potete vedere, abbiamo impostato una fabbrica di canali di base per caricare il foglio dati di esempio come file, ma questo non leggerà ancora il contenuto del file.
Iniziamo aggiungendo questa funzionalità.

### 1.1. Leggere il contenuto con `splitCsv`

Dobbiamo scegliere un operatore che analizzi il contenuto del file in modo appropriato con il minimo sforzo da parte nostra.
Poiché il nostro foglio dati è in formato CSV, questo è un lavoro per l'operatore [`splitCsv`](https://www.nextflow.io/docs/latest/reference/operator.html#splitcsv), che carica ogni riga del file come elemento nel canale.

Apportate le seguenti modifiche per aggiungere un'operazione `splitCsv()` al codice di costruzione del canale, più un'operazione `view()` per verificare che il contenuto del file venga caricato correttamente nel canale.

=== "Dopo"

    ```groovy title="main.nf" linenums="3" hl_lines="4-5"
    workflow  {

        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .view()

    }
    ```

=== "Prima"

    ```groovy title="main.nf" linenums="3"
    workflow  {

        ch_datasheet = channel.fromPath("./data/datasheet.csv")

    }
    ```

Notate che stiamo usando l'opzione `header: true` per indicare a Nextflow di leggere la prima riga del file CSV come riga di intestazione.

Vediamo cosa ne esce, vero?
Eseguiamo il flusso di lavoro:

```bash
nextflow run main.nf
```

??? success "Output del comando"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [exotic_albattani] DSL2 - revision: c0d03cec83

    [id:sampleA, character:squirrel, recording:/workspaces/training/side-quests/metadata/data/bonjour.txt]
    [id:sampleB, character:tux, recording:/workspaces/training/side-quests/metadata/data/guten_tag.txt]
    [id:sampleC, character:sheep, recording:/workspaces/training/side-quests/metadata/data/hallo.txt]
    [id:sampleD, character:turkey, recording:/workspaces/training/side-quests/metadata/data/hello.txt]
    [id:sampleE, character:stegosaurus, recording:/workspaces/training/side-quests/metadata/data/hola.txt]
    [id:sampleF, character:moose, recording:/workspaces/training/side-quests/metadata/data/salut.txt]
    [id:sampleG, character:turtle, recording:/workspaces/training/side-quests/metadata/data/ciao.txt]
    ```

Possiamo vedere che l'operatore ha costruito una mappa di coppie chiave-valore per ogni riga del file CSV, con le intestazioni delle colonne come chiavi per i valori corrispondenti.

Ogni voce della mappa corrisponde a una colonna nel nostro foglio dati:

- `id`
- `character`
- `recording`

Ottimo! Questo rende facile accedere a campi specifici da ogni file.
Ad esempio, potremmo accedere all'ID del file con `id` o al percorso del file txt con `recording`.

??? info "(Opzionale) Ulteriori informazioni sulle map"

    In Groovy, il linguaggio di programmazione su cui è costruito Nextflow, una map è una struttura dati chiave-valore simile ai dizionari in Python, agli oggetti in JavaScript o agli hash in Ruby.

    Ecco uno script eseguibile che mostra come definire una map e accedere al suo contenuto in pratica:

    ```groovy title="examples/map_demo.nf"
    #!/usr/bin/env nextflow

    // Crea una semplice map
    def my_map = [id:'sampleA', character:'squirrel']

    // Stampa l'intera map
    println "map: ${my_map}"

    // Accede ai singoli valori usando la notazione con punto
    println "id: ${my_map.id}"
    println "character: ${my_map.character}"
    ```

    Anche se non ha un blocco `workflow` vero e proprio, Nextflow può eseguirlo come se fosse un flusso di lavoro:

    ```bash
    nextflow run examples/map_demo.nf
    ```

    Ed ecco cosa potete aspettarvi di vedere nell'output:

    ```console title="Output"
     N E X T F L O W   ~  version 25.10.2

    Launching `map_demo.nf` [cheesy_plateau] DSL2 - revision: fae5b8496e

    map: [id:sampleA, character:squirrel]
    id: sampleA
    character: squirrel
    ```

### 1.2. Selezionare campi specifici con `map`

Supponiamo di voler accedere alla colonna `character` del foglio dati e stamparla.
Possiamo usare l'operatore `map` di Nextflow per iterare su ogni elemento nel nostro canale e selezionare specificamente la voce `character` dall'oggetto map.

Apportate le seguenti modifiche al flusso di lavoro:

=== "Dopo"

    ```groovy title="main.nf" linenums="3" hl_lines="5-7"
    workflow  {

        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map{ row ->
                row.character
            }
            .view()

    }
    ```

=== "Prima"

    ```groovy title="main.nf" linenums="3"
    workflow  {

        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .view()

    }
    ```

Ora eseguiamo di nuovo il flusso di lavoro:

```bash
nextflow run main.nf
```

??? success "Output del comando"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [exotic_albattani] DSL2 - revision: c0d03cec83

    squirrel
    tux
    sheep
    turkey
    stegosaurus
    moose
    turtle
    ```

Ottimo! Abbiamo sfruttato la struttura map derivata dal nostro foglio dati per accedere ai valori delle singole colonne per ogni riga.

Ora che abbiamo letto con successo il foglio dati e abbiamo accesso ai dati in ogni riga, possiamo iniziare a implementare la logica del nostro pipeline.

### 1.3. Organizzare i metadati in una 'meta map'

Nello stato attuale del flusso di lavoro, i file di input (sotto la chiave `recording`) e i metadati associati (`id`, `character`) sono tutti sullo stesso piano, come se fossero tutti in un unico grande sacchetto.
La conseguenza pratica è che ogni processo che consuma questo canale dovrebbe essere configurato tenendo presente questa struttura:

```groovy
    input:
    tuple val(id), val(character), file(recording)
```

Va bene finché il numero di colonne nel foglio dati non cambia.
Tuttavia, se aggiungete anche solo una colonna al foglio dati, la forma del canale non corrisponderà più a ciò che il processo si aspetta, e il flusso di lavoro produrrà errori.
Rende anche il processo difficile da condividere con altri che potrebbero avere dati di input leggermente diversi, e potreste finire per dover codificare variabili nel processo che non sono necessarie per il blocco script.

Per evitare questo problema, dobbiamo trovare un modo per mantenere la struttura del canale coerente indipendentemente dal numero di colonne che contiene il foglio dati.

Possiamo farlo raccogliendo tutti i metadati in un elemento all'interno della tupla, che chiameremo la mappa dei metadati, o più semplicemente 'meta map'.

Apportate le seguenti modifiche all'operazione `map`:

=== "Dopo"

    ```groovy title="main.nf" linenums="5" hl_lines="4"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                    [ [id: row.id, character: row.character], row.recording ]
            }
            .view()
    ```

=== "Prima"

    ```groovy title="main.nf" linenums="5" hl_lines="4"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map{ row ->
                row.character
            }
            .view()
    ```

Abbiamo ristrutturato gli elementi del canale in una tupla composta da due elementi: la meta map e l'oggetto file corrispondente.

Eseguiamo il flusso di lavoro:

```bash
nextflow run main.nf
```

??? success "Output del comando"

    ```console title="View meta map"
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [lethal_booth] DSL2 - revision: 0d8f844c07

    [[id:sampleA, character:squirrel], /workspaces/training/side-quests/metadata/data/bonjour.txt]
    [[id:sampleB, character:tux], /workspaces/training/side-quests/metadata/data/guten_tag.txt]
    [[id:sampleC, character:sheep], /workspaces/training/side-quests/metadata/data/hallo.txt]
    [[id:sampleD, character:turkey], /workspaces/training/side-quests/metadata/data/hello.txt]
    [[id:sampleE, character:stegosaurus], /workspaces/training/side-quests/metadata/data/hola.txt]
    [[id:sampleF, character:moose], /workspaces/training/side-quests/metadata/data/salut.txt]
    [[id:sampleG, character:turtle], /workspaces/training/side-quests/metadata/data/ciao.txt]
    ```

Ora, ogni elemento nel canale contiene prima la meta map e poi l'oggetto file corrispondente:

```console title="Example output structure"
[
  [id:sampleA, character:squirrel],
  /workspaces/training/side-quests/metadata/data/bonjour.txt
]
```

Di conseguenza, aggiungere più colonne nel foglio dati renderà disponibili più metadati nella mappa `meta`, ma non cambierà la forma del canale.
Questo ci permette di scrivere processi che consumano il canale senza dover codificare gli elementi dei metadati nella specifica di input:

```groovy title="Syntax example"
    input:
    tuple val(meta), file(recording)
```

Questa è una convenzione ampiamente utilizzata per organizzare i metadati nei flussi di lavoro Nextflow.

### Takeaway

In questa sezione, avete imparato:

- **Perché i metadati sono importanti:** Mantenere i metadati insieme ai dati preserva le informazioni importanti sui file durante tutto il flusso di lavoro.
- **Come leggere i fogli dati:** Usare `splitCsv` per leggere file CSV con informazioni di intestazione e trasformare le righe in dati strutturati
- **Come creare una meta map:** Separare i metadati dai dati dei file usando la struttura tupla `[ [id:valore, ...], file ]`

---

## 2. Manipolare i metadati

Ora che abbiamo caricato i nostri metadati, facciamo qualcosa con essi!

Useremo uno strumento chiamato [`langid`](https://github.com/saffsd/langid.py) per identificare la lingua contenuta nel file di registrazione di ogni creatura.
Lo strumento viene fornito pre-addestrato su un insieme di lingue e, dato un frammento di testo, produrrà una previsione della lingua e un punteggio di probabilità associato, entrambi su `stdout`.

### 2.1. Importare il processo ed esaminare il codice

Vi forniamo un modulo di processo pre-scritto chiamato `IDENTIFY_LANGUAGE` che racchiude lo strumento `langid`, quindi dovete solo aggiungere un'istruzione include prima del blocco workflow.

Apportate la seguente modifica al flusso di lavoro:

=== "Dopo"

    ```groovy title="main.nf" linenums="1" hl_lines="3"
    #!/usr/bin/env nextflow

    include { IDENTIFY_LANGUAGE } from './modules/langid.nf'

    workflow {
    ```

=== "Prima"

    ```groovy title="main.nf" linenums="1"
    #!/usr/bin/env nextflow

    workflow {
    ```

Potete aprire il file del modulo per esaminarne il codice:

```groovy title="modules/langid.nf" linenums="1" hl_lines="9 12"
#!/usr/bin/env nextflow

// Usa langid per prevedere la lingua di ogni file di input
process IDENTIFY_LANGUAGE {

    container 'community.wave.seqera.io/library/pip_langid:b2269f456a5629ff'

    input:
    tuple val(meta), path(file)

    output:
    tuple val(meta), path(file), stdout

    script:
    """
    langid < ${file} -l en,de,fr,es,it | sed -E "s/.*\\('([a-z]+)'.*/\\1/" | tr -d '\\n'
    """
}
```

Come potete vedere, la definizione di input usa la stessa struttura `tuple val(meta), path(file)` che abbiamo appena applicato al nostro canale di input.

La definizione di output è strutturata come una tupla con una struttura simile a quella dell'input, tranne per il fatto che contiene anche `stdout` come terzo elemento.
Questo schema `tuple val(meta), path(file), <output>` mantiene i metadati associati sia ai dati di input che agli output mentre scorrono attraverso il pipeline.

Notate che stiamo usando il qualificatore di output [`stdout`](https://www.nextflow.io/docs/latest/process.html#outputs) di Nextflow perché lo strumento stampa il suo output direttamente sulla console invece di scrivere un file; e usiamo `sed` nella riga di comando per rimuovere il punteggio di probabilità, pulire la stringa rimuovendo i caratteri di nuova riga e restituire solo la previsione della lingua.

### 2.2. Aggiungere una chiamata a `IDENTIFY_LANGUAGE`

Ora che il processo è disponibile per il flusso di lavoro, possiamo aggiungere una chiamata al processo `IDENTIFY_LANGUAGE` per eseguirlo sul canale dati.

Apportate le seguenti modifiche al flusso di lavoro:

=== "Dopo"

    ```groovy title="main.nf" linenums="7" hl_lines="7-9"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [[id: row.id, character: row.character], row.recording]
            }

        // Esegue langid per identificare la lingua di ogni saluto
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out.view()
    ```

=== "Prima"

    ```groovy title="main.nf" linenums="7" hl_lines="6"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [[id: row.id, character: row.character], row.recording]
            }
            .view()
    ```

Notate che abbiamo rimosso l'operazione `.view()` originale nella costruzione del canale.

Ora possiamo eseguire il flusso di lavoro:

```bash
nextflow run main.nf
```

??? success "Output del comando"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [voluminous_mcnulty] DSL2 - revision: f9bcfebabb

    executor >  local (7)
    [4e/f722fe] IDENTIFY_LANGUAGE (7) [100%] 7 of 7 ✔
    [[id:sampleA, character:squirrel], /workspaces/training/side-quests/metadata/work/eb/f7148ebdd898fbe1136bec6a714acb/bonjour.txt, fr]
    [[id:sampleB, character:tux], /workspaces/training/side-quests/metadata/work/16/71d72410952c22cd0086d9bca03680/guten_tag.txt, de]
    [[id:sampleD, character:turkey], /workspaces/training/side-quests/metadata/work/c4/b7562adddc1cc0b7d414ec45d436eb/hello.txt, en]
    [[id:sampleC, character:sheep], /workspaces/training/side-quests/metadata/work/ea/04f5d979429e4455e14b9242fb3b45/hallo.txt, de]
    [[id:sampleF, character:moose], /workspaces/training/side-quests/metadata/work/5a/6c2b84bf8fadb98e28e216426be079/salut.txt, fr]
    [[id:sampleE, character:stegosaurus], /workspaces/training/side-quests/metadata/work/af/ee7c69bcab891c40d0529305f6b9e7/hola.txt, es]
    [[id:sampleG, character:turtle], /workspaces/training/side-quests/metadata/work/4e/f722fe47271ba7ebcd69afa42964ca/ciao.txt, it]
    ```

Eccellente! Ora abbiamo una previsione su quale lingua parla ogni personaggio.

E come notato in precedenza, abbiamo anche incluso il file di input e la meta map nell'output, il che significa che entrambi rimangono associati alle nuove informazioni che abbiamo appena prodotto.
Questo si rivelerà utile nel passo successivo.

!!! note "Nota"

    Più in generale, questo schema di mantenere la meta map associata ai risultati rende più facile associare risultati correlati che condividono gli stessi identificatori.

    Come avrete già imparato, non potete fare affidamento sull'ordine degli elementi nei canali per abbinare i risultati tra di essi.
    Dovete invece usare le chiavi per associare correttamente i dati, e le meta map forniscono una struttura ideale per questo scopo.

    Esploriamo questo caso d'uso in dettaglio nella side quest [Splitting & Grouping](../splitting_and_grouping/).

### 2.3. Arricchire i metadati con gli output dei processi

Dato che i risultati che abbiamo appena prodotto sono di per sé una forma di metadati sul contenuto dei file, sarebbe utile aggiungerli alla nostra meta map.

Tuttavia, non vogliamo modificare la meta map esistente sul posto.
Da un punto di vista tecnico, è _possibile_ farlo, ma non è sicuro.

Quindi, invece, creeremo una nuova meta map contenente il contenuto della meta map esistente più una nuova coppia chiave-valore `lang: lang_id` che contiene le nuove informazioni, usando l'operatore `+` (una funzionalità di Groovy).
E combineremo questo con un'operazione [`map`](https://www.nextflow.io/docs/latest/operator.html#map) per sostituire la vecchia mappa con quella nuova.

Ecco le modifiche che dovete apportare al flusso di lavoro:

=== "Dopo"

    ```groovy title="main.nf" linenums="13" hl_lines="3-7"
        // Esegue langid per identificare la lingua di ogni saluto
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out
            .map { meta, file, lang_id ->
                [meta + [lang: lang_id], file]
            }
            .view()
    ```

=== "Prima"

    ```groovy title="main.nf" linenums="13" hl_lines="3"
        // Esegue langid per identificare la lingua di ogni saluto
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out.view()
    ```

Se non avete ancora familiarità con l'operatore `+`, o se questo vi sembra confuso, prendetevi qualche minuto per leggere la spiegazione dettagliata qui sotto.

??? info "Creazione della nuova meta map usando l'operatore `+`"

    **Prima di tutto, dovete sapere che possiamo unire il contenuto di due map usando l'operatore Groovy `+`.**

    Supponiamo di avere le seguenti map:

    ```groovy
    map1 = [id: 'sampleA', character: 'squirrel']
    map2 = [lang: 'fr']
    ```

    Possiamo unirle così:

    ```groovy
    new_map = map1 + map2
    ```

    Il contenuto di `new_map` sarà:

    ```groovy
    [id: 'sampleA', character: 'squirrel', lang: 'fr']
    ```

    Ottimo!

    **Ma cosa succede se dovete aggiungere un campo che non fa già parte di una map?**

    Supponiamo di ricominciare da `map1`, ma la previsione della lingua non è nella sua propria map (non esiste `map2`).
    Invece, è contenuta in una variabile chiamata `lang_id`, e sapete di voler memorizzare il suo valore (`'fr'`) con la chiave `lang`.

    Potete effettivamente fare quanto segue:

    ```groovy
    new_map = [map1 + [lang: lang_id]]
    ```

    Qui, `[lang: new_info]` crea una nuova map senza nome al volo, e `map1 + ` unisce `map1` con la nuova map senza nome, producendo lo stesso contenuto di `new_map` come prima.

    Elegante, vero?

    **Ora trasponiamo questo nel contesto di un'operazione `channel.map()` di Nextflow.**

    Il codice diventa:

    ```groovy
    .map { map1, lang_id ->
        [map1 + [lang: lang_id]]
    }
    ```

    Questo fa quanto segue:

    - `map1, lang_id ->` prende i due elementi nella tupla
    - `[map1 + [lang: lang_id]]` crea la nuova map come descritto sopra

    L'output è una singola map senza nome con lo stesso contenuto di `new_map` nel nostro esempio precedente.
    Quindi abbiamo effettivamente trasformato:

    ```groovy
    [id: 'sampleA', character: 'squirrel'], 'fr'
    ```

    in:

    ```groovy
    [id: 'sampleA', character: 'squirrel', lang: 'fr']
    ```

    Speriamo che possiate vedere che se cambiamo `map1` in `meta`, è fondamentalmente tutto ciò di cui abbiamo bisogno per aggiungere la previsione della lingua alla nostra meta map nel nostro flusso di lavoro.

    Tranne per una cosa!

    Nel caso del nostro flusso di lavoro, **dobbiamo anche tenere conto della presenza dell'oggetto `file` nella tupla**, che è composta da `meta, file, lang_id`.

    Quindi il codice qui diventerebbe:

    ```groovy
    .map { meta, file, lang_id ->
        [meta + [lang: lang_id], file]
    }
    ```

    Se avete difficoltà a capire perché il `file` sembra spostarsi nell'operazione `map`, immaginate che invece di `[meta + [lang: lang_id], file]`, quella riga reciti `[new_map, file]`.
    Questo dovrebbe rendere più chiaro che stiamo semplicemente lasciando il `file` nella sua posizione originale in seconda posizione nella tupla. Abbiamo semplicemente preso il valore `new_info` e lo abbiamo incorporato nella mappa che si trova in prima posizione.

    **E questo ci riporta alla struttura del canale `tuple val(meta), path(file)`!**

Una volta che siete sicuri di capire cosa fa questo codice, eseguite il flusso di lavoro per vedere se ha funzionato:

```bash
nextflow run main.nf -resume
```

??? success "Output del comando"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [cheeky_fermat] DSL2 - revision: d096281ee4

    [4e/f722fe] IDENTIFY_LANGUAGE (7) [100%] 7 of 7, cached: 7 ✔
    [[id:sampleA, character:squirrel, lang:fr], /workspaces/training/side-quests/metadata/work/eb/f7148ebdd898fbe1136bec6a714acb/bonjour.txt]
    [[id:sampleB, character:tux, lang:de], /workspaces/training/side-quests/metadata/work/16/71d72410952c22cd0086d9bca03680/guten_tag.txt]
    [[id:sampleC, character:sheep, lang:de], /workspaces/training/side-quests/metadata/work/ea/04f5d979429e4455e14b9242fb3b45/hallo.txt]
    [[id:sampleD, character:turkey, lang:en], /workspaces/training/side-quests/metadata/work/c4/b7562adddc1cc0b7d414ec45d436eb/hello.txt]
    [[id:sampleF, character:moose, lang:fr], /workspaces/training/side-quests/metadata/work/5a/6c2b84bf8fadb98e28e216426be079/salut.txt]
    [[id:sampleE, character:stegosaurus, lang:es], /workspaces/training/side-quests/metadata/work/af/ee7c69bcab891c40d0529305f6b9e7/hola.txt]
    [[id:sampleG, character:turtle, lang:it], /workspaces/training/side-quests/metadata/work/4e/f722fe47271ba7ebcd69afa42964ca/ciao.txt]
    ```

Sì, funziona!
Abbiamo riorganizzato ordinatamente l'output del processo da `meta, file, lang_id` in modo che `lang_id` sia ora una delle chiavi nella meta map, e le tuple del canale si adattano nuovamente al modello `meta, file`.

<!-- TODO (future) Should we also show how to remove a key using subMap?! Or note where to find that. -->

### 2.4. Assegnare un gruppo linguistico usando i condizionali

Ora che abbiamo le nostre previsioni linguistiche, usiamo le informazioni per assegnare nuovi raggruppamenti.

Nei nostri dati di esempio, le lingue usate dai nostri personaggi possono essere raggruppate in lingue germaniche (inglese, tedesco) e lingue romanze (francese, spagnolo, italiano).
Potrebbe essere utile avere quella classificazione prontamente disponibile in seguito nel pipeline, quindi aggiungiamo queste informazioni nella meta map.

E, buone notizie, questo è un altro caso che si presta perfettamente all'uso dell'operatore `map`!

Nello specifico, definiremo una variabile chiamata `lang_group`, useremo una semplice logica condizionale per determinare quale valore assegnare a `lang_group` per ogni dato.

La sintassi generale sarà simile a questa:

```groovy
.map { meta, file ->

    // la logica condizionale che definisce lang_group va qui

    [meta + [lang_group: lang_group], file]
}
```

Potete vedere che questo è molto simile all'operazione di unione di map al volo che abbiamo usato nel passo precedente.
Dobbiamo solo scrivere le istruzioni condizionali.

Ecco la logica condizionale che vogliamo applicare:

- Definire una variabile chiamata `lang_group` con valore predefinito `'unknown'`.
- Se `lang` è tedesco (`'de'`) o inglese (`'en'`), cambiare `lang_group` in `germanic`.
- Altrimenti se `lang` è incluso in un elenco contenente francese (`'fr'`), spagnolo (`'es'`) e italiano (`'it'`), cambiare `lang_group` in `romance`.

Provate a scriverlo voi stessi se sapete già come scrivere istruzioni condizionali in Nextflow.

!!! tip "Suggerimento"

    Potete accedere al valore di `lang` all'interno dell'operazione map con `meta.lang`.

Dovreste finire per apportare le seguenti modifiche al flusso di lavoro:

=== "Dopo"

    ```groovy title="main.nf" linenums="13" hl_lines="7-19"
        // Esegue langid per identificare la lingua di ogni saluto
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out
            .map { meta, file, lang_id ->
                [meta + [lang: lang_id], file]
            }
            .map { meta, file ->

                def lang_group = "unknown"
                if (meta.lang.equals("de") || meta.lang.equals('en')) {
                    lang_group = "germanic"
                }
                else if (meta.lang in ["fr", "es", "it"]) {
                    lang_group = "romance"
                }

                [meta + [lang_group: lang_group], file]
            }
            .view()
    ```

=== "Prima"

    ```groovy title="main.nf" linenums="13" hl_lines="7"
        // Esegue langid per identificare la lingua di ogni saluto
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out
            .map { meta, file, lang_id ->
                [meta + [lang: lang_id], file]
            }
            .view()
    ```

Ecco i punti chiave:

- Usiamo `def lang_group = "unknown"` per creare la variabile `lang_group` con valore predefinito impostato su `unknown`.
- Usiamo una struttura `if {} else if {}` per la logica condizionale, con test `.equals()` alternativi per le due lingue germaniche, e un test di esistenza in un elenco per le tre lingue romanze.
- Usiamo l'operazione di unione `meta + [lang_group:lang_group]` come in precedenza per generare la meta map aggiornata.

<!-- TODO (future) Add note/links to relevant docs in additional resources section -->

Una volta che tutto ha senso, eseguite di nuovo il flusso di lavoro per vedere il risultato:

```bash
nextflow run main.nf -resume
```

??? success "Output del comando"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [wise_almeida] DSL2 - revision: 46778c3cd0

    [da/652cc6] IDENTIFY_LANGUAGE (7) [100%] 7 of 7, cached: 7 ✔
    [[id:sampleA, character:squirrel, lang:fr, lang_group:romance], /workspaces/training/side-quests/metadata/data/bonjour.txt]
    [[id:sampleB, character:tux, lang:de, lang_group:germanic], /workspaces/training/side-quests/metadata/data/guten_tag.txt]
    [[id:sampleC, character:sheep, lang:de, lang_group:germanic], /workspaces/training/side-quests/metadata/data/hallo.txt]
    [[id:sampleD, character:turkey, lang:en, lang_group:germanic], /workspaces/training/side-quests/metadata/data/hello.txt]
    [[id:sampleE, character:stegosaurus, lang:es, lang_group:romance], /workspaces/training/side-quests/metadata/data/hola.txt]
    [[id:sampleF, character:moose, lang:fr, lang_group:romance], /workspaces/training/side-quests/metadata/data/salut.txt]
    [[id:sampleG, character:turtle, lang:it, lang_group:romance], /workspaces/training/side-quests/metadata/data/ciao.txt]
    ```

Come potete vedere, gli elementi del canale mantengono la loro struttura `[meta, file]`, ma la meta map ora include questa nuova classificazione.

### Takeaway

In questa sezione, avete imparato come:

- **Applicare i metadati di input ai canali di output**: Copiare i metadati in questo modo ci permette di associare i risultati in seguito in base al contenuto dei metadati.
- **Creare chiavi personalizzate**: Avete creato due nuove chiavi nella vostra meta map, unendole con `meta + [nuova_chiave:valore]` alla meta map esistente. Una basata su un valore calcolato da un processo, e una basata su una condizione impostata nell'operatore `map`.

Queste vi permettono di associare metadati nuovi ed esistenti ai file man mano che progredite nel vostro pipeline.
Anche se non state usando i metadati come parte di un processo, mantenere la meta map associata ai dati in questo modo rende facile tenere insieme tutte le informazioni rilevanti.

---

## 3. Utilizzare le informazioni della meta map in un processo

Ora che sapete come creare e aggiornare la meta map, possiamo arrivare alla parte davvero divertente: usare effettivamente i metadati in un processo.

Più specificamente, aggiungeremo un secondo passo al nostro flusso di lavoro per disegnare ogni animale come ASCII art e farlo pronunciare il testo registrato in un fumetto.
Lo faremo usando uno strumento chiamato [`cowpy`](https://github.com/jeffbuttars/cowpy).

??? info "Cosa fa `cowpy`?"

    `cowpy` è uno strumento da riga di comando che genera ASCII art per visualizzare input di testo arbitrari in modo divertente.
    È un'implementazione Python del classico strumento [cowsay](https://en.wikipedia.org/wiki/Cowsay) di Tony Monroe.

    ```console
    cowpy "Hello Nextflow"
    ```

    ```console
    ______________________________________________________
    < Hello Nextflow >
    ------------------------------------------------------
        \   ^__^
          \  (oo)\_______
            (__)\       )\/\
              ||----w |
              ||     ||
    ```

    Facoltativamente, potete selezionare un personaggio (o 'cowacter') da usare al posto della mucca predefinita.

    ```console
    cowpy "Hello Nextflow" -c tux
    ```

    ```console
    __________________
    < Hello Nextflow >
    ------------------
      \
        \
            .--.
          |o_o |
          |:_/ |
          //   \ \
        (|     | )
        /'\_   _/`\
        \___)=(___/
    ```

Se avete seguito il corso Hello Nextflow, avete già visto questo strumento in azione.
In caso contrario, non preoccupatevi; copriremo tutto ciò che dovete sapere man mano che procediamo.

### 3.1. Importare il processo ed esaminare il codice

Vi forniamo un modulo di processo pre-scritto chiamato `COWPY` che racchiude lo strumento `cowpy`, quindi dovete solo aggiungere un'istruzione include prima del blocco workflow.

Apportate la seguente modifica al flusso di lavoro:

=== "Dopo"

    ```groovy title="main.nf" linenums="1" hl_lines="4"
    #!/usr/bin/env nextflow

    include { IDENTIFY_LANGUAGE } from './modules/langid.nf'
    include { COWPY } from './modules/cowpy.nf'

    workflow {
    ```

=== "Prima"

    ```groovy title="main.nf" linenums="1"
    #!/usr/bin/env nextflow

    include { IDENTIFY_LANGUAGE } from './modules/langid.nf'

    workflow {
    ```

Potete aprire il file del modulo per esaminarne il codice:

```groovy title="modules/cowpy.nf" linenums="1"
#!/usr/bin/env nextflow

// Genera ASCII art con cowpy
process COWPY {

    publishDir "results/", mode: 'copy'

    container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

    input:
    path input_file
    val character

    output:
    path "cowpy-${input_file}"

    script:
    """
    cat ${input_file} | cowpy -c ${character} > cowpy-${input_file}
    """
}
```

Come potete vedere, questo processo è attualmente progettato per ricevere un file di input (contenente il testo da visualizzare) e un valore che specifica il personaggio da disegnare in ASCII art, solitamente fornito a livello di flusso di lavoro tramite un parametro da riga di comando.

### 3.2. Passare un campo della meta map come input

Quando abbiamo usato lo strumento `cowpy` nel corso Hello Nextflow, abbiamo usato un parametro da riga di comando per determinare quale personaggio usare per disegnare l'immagine finale.
Aveva senso, perché stavamo generando solo un'immagine per ogni esecuzione del pipeline.

Tuttavia, in questo tutorial, vogliamo generare un'immagine appropriata per ogni soggetto che stiamo elaborando, quindi usare un parametro da riga di comando sarebbe troppo limitante.

Buone notizie: abbiamo una colonna `character` nel nostro foglio dati e quindi nella nostra meta map.
Usiamola per impostare il personaggio che il processo dovrebbe usare per ogni voce.

A tal fine, dovremo fare tre cose:

1. Dare un nome al canale di output proveniente dal processo precedente in modo da poterci operare più comodamente.
2. Determinare come accedere alle informazioni di interesse
3. Aggiungere una chiamata al secondo processo e fornire le informazioni in modo appropriato.

Iniziamo.

#### 3.2.1. Nominare il canale di output precedente

Abbiamo applicato le manipolazioni precedenti direttamente sul canale di output del primo processo, `IDENTIFY_LANGUAGE.out`.
Per poter alimentare il contenuto di quel canale al processo successivo (e farlo in modo chiaro e facile da leggere) vogliamo dargli un nome proprio, `ch_languages`.

Possiamo farlo usando l'operatore [`set`](https://www.nextflow.io/docs/latest/reference/operator.html#set).

Nel flusso di lavoro principale, sostituite l'operatore `.view()` con `.set { ch_languages }`, e aggiungete una riga per verificare che possiamo fare riferimento al canale per nome.

=== "Dopo"

    ```groovy title="main.nf" linenums="14" hl_lines="19 21 22"
        // Esegue langid per identificare la lingua di ogni saluto
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out
            .map { meta, file, lang_id ->
                [meta + [lang: lang_id], file]
            }
            .map { meta, file ->

                def lang_group = "unknown"
                if (meta.lang.equals("de") || meta.lang.equals('en')) {
                    lang_group = "germanic"
                }
                else if (meta.lang in ["fr", "es", "it"]) {
                    lang_group = "romance"
                }

                [meta + [lang_group: lang_group], file]
            }
            .set { ch_languages }

        // Temporaneo: sbirciare in ch_languages
        ch_languages.view()
    ```

=== "Prima"

    ```groovy title="main.nf" linenums="14" hl_lines="19"
        // Esegue langid per identificare la lingua di ogni saluto
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out
            .map { meta, file, lang_id ->
                [meta + [lang: lang_id], file]
            }
            .map { meta, file ->

                def lang_group = "unknown"
                if (meta.lang.equals("de") || meta.lang.equals('en')) {
                    lang_group = "germanic"
                }
                else if (meta.lang in ["fr", "es", "it"]) {
                    lang_group = "romance"
                }

                [meta + [lang_group: lang_group], file]
            }
            .view()
    ```

Eseguiamo:

```bash
nextflow run main.nf
```

??? success "Output del comando"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `./main.nf` [friendly_austin] DSL2 - revision: 3dbe460fd6

    [36/cca6a7] IDENTIFY_LANGUAGE (7) | 7 of 7 ✔
    [[id:sampleB, character:tux, lang:de, lang_group:germanic], /workspaces/training/side-quests/metadata/work/e2/6db2402d83cf72081bcd2d11784714/guten_tag.txt]
    [[id:sampleA, character:squirrel, lang:fr, lang_group:romance], /workspaces/training/side-quests/metadata/work/6c/114c818317d169457d6e7336d5d55b/bonjour.txt]
    [[id:sampleC, character:sheep, lang:de, lang_group:germanic], /workspaces/training/side-quests/metadata/work/55/68c69c5efb527f3604ddb3daab8057/hallo.txt]
    [[id:sampleD, character:turkey, lang:en, lang_group:germanic], /workspaces/training/side-quests/metadata/work/2a/4752055ccb5d1370b0ef9da41d3993/hello.txt]
    [[id:sampleE, character:stegosaurus, lang:es, lang_group:romance], /workspaces/training/side-quests/metadata/work/f4/fcd3186dc666d5d239ffa6c37d125d/hola.txt]
    [[id:sampleF, character:moose, lang:fr, lang_group:romance], /workspaces/training/side-quests/metadata/work/c3/3b2627f733f278a7088332a5806108/salut.txt]
    [[id:sampleG, character:turtle, lang:it, lang_group:romance], /workspaces/training/side-quests/metadata/work/36/cca6a7dbfa26ac24f9329787a32e9d/ciao.txt]
    ```

Questo conferma che ora possiamo fare riferimento al canale per nome.

#### 3.2.2. Accedere ai metadati del file e del personaggio

Sappiamo dal codice del modulo che il processo `COWPY` si aspetta di ricevere un file di testo e un valore `character`.
Per scrivere la chiamata al processo `COWPY`, dobbiamo solo sapere come estrarre l'oggetto file corrispondente e i metadati da ogni elemento nel canale.

Come spesso accade, il modo più semplice per farlo è usare un'operazione `map`.

Il nostro canale contiene tuple strutturate come `[meta, file]`, quindi possiamo accedere direttamente all'oggetto `file`, e possiamo accedere al valore `character` memorizzato all'interno della meta map facendo riferimento ad esso come `meta.character`.

Nel flusso di lavoro principale, apportate le seguenti modifiche al codice:

=== "Dopo"

    ```groovy title="main.nf" linenums="34"
        // Temporaneo: accedere al file e al personaggio
        ch_languages.map { meta, file -> file }.view { file -> "File: " + file }
        ch_languages.map { meta, file -> meta.character }.view { character -> "Character: " + character }
    ```

=== "Prima"

    ```groovy title="main.nf" linenums="34"
        // Temporaneo: sbirciare in ch_languages
        ch_languages.view()
    ```

Notate che stiamo usando closure (come `{ file -> "File: " + file }`) per rendere l'output delle operazioni `.view` più leggibile.

Eseguiamo:

```bash
nextflow run main.nf -resume
```

??? success "Output del comando"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `./main.nf` [cheesy_cantor] DSL2 - revision: 15af9c1ec7

    [43/05df08] IDENTIFY_LANGUAGE (7) [100%] 7 of 7, cached: 7 ✔
    Character: squirrel
    File: /workspaces/training/side-quests/metadata/work/8d/4b9498bbccb7a74f04e41877cdc3e5/bonjour.txt
    File: /workspaces/training/side-quests/metadata/work/d3/604274985406e40d79021dea658e60/guten_tag.txt
    Character: tux
    Character: turkey
    File: /workspaces/training/side-quests/metadata/work/d4/fafcc9415b61d2b0fea872e6a05e8a/hello.txt
    File: /workspaces/training/side-quests/metadata/work/02/468ac9efb27f636715e8144b37e9a7/hallo.txt
    Character: sheep
    Character: moose
    Character: stegosaurus
    File: /workspaces/training/side-quests/metadata/work/d4/61a7e1188b4f2742bc72004e226eca/salut.txt
    File: /workspaces/training/side-quests/metadata/work/ae/68364be238c11149c588bf6fc858b1/hola.txt
    File: /workspaces/training/side-quests/metadata/work/43/05df081af5d879ab52e5828fa0357e/ciao.txt
    Character: turtle
    ```

_I percorsi dei file e i valori dei personaggi potrebbero apparire in un ordine diverso nel vostro output._

Questo conferma che siamo in grado di accedere al file e al personaggio per ogni elemento nel canale.

#### 3.2.3. Chiamare il processo `COWPY`

Ora mettiamo tutto insieme e chiamiamo effettivamente il processo `COWPY` sul canale `ch_languages`.

Nel flusso di lavoro principale, apportate le seguenti modifiche al codice:

=== "Dopo"

    ```groovy title="main.nf" linenums="34"
        // Esegue cowpy per generare ASCII art
        COWPY(
            ch_languages.map { meta, file -> file },
            ch_languages.map { meta, file -> meta.character }
        )
    ```

=== "Prima"

    ```groovy title="main.nf" linenums="34"
        // Temporaneo: accedere al file e al personaggio
        ch_languages.map { meta, file -> [file, meta.character] }
            .view()
    ```

Vedete che copiamo semplicemente le due operazioni map (senza le istruzioni `.view()`) come input alla chiamata del processo.
Assicuratevi solo di non dimenticare la virgola tra di esse!

È un po' macchinoso, ma vedremo come migliorarlo nella prossima sezione.

Eseguiamo:

```bash
nextflow run main.nf -resume
```

??? success "Output del comando"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [suspicious_crick] DSL2 - revision: 25541014c5

    executor >  local (7)
    [43/05df08] IDENTIFY_LANGUAGE (7) [100%] 7 of 7, cached: 7 ✔
    [e7/317c18] COWPY (6)             [100%] 7 of 7 ✔
    ```

Se guardate nella directory dei risultati, dovreste vedere i singoli file contenenti l'ASCII art di ogni saluto pronunciato dal personaggio corrispondente.

??? abstract "Contenuto della directory e del file di esempio"

    ```console
    results/
    ├── cowpy-bonjour.txt
    ├── cowpy-ciao.txt
    ├── cowpy-guten_tag.txt
    ├── cowpy-hallo.txt
    ├── cowpy-hello.txt
    ├── cowpy-hola.txt
    └── cowpy-salut.txt
    ```

    ```text title="results/cowpy-bonjour.txt"
     _________________
    / Bonjour         \
    \ Salut, à demain /
    -----------------
      \
        \
                      _ _
          | \__/|  .~    ~.
          /oo `./      .'
          {o__,   \    {
            / .  . )    \
            `-` '-' \    }
          .(   _(   )_.'
          '---.~_ _ _|
    ```

Questo dimostra che siamo stati in grado di usare le informazioni nella meta map per parametrizzare il comando nel secondo passo del pipeline.

Tuttavia, come notato sopra, parte del codice coinvolto era un po' macchinoso, poiché dovevamo estrarre i metadati ancora nel contesto del corpo del flusso di lavoro.
Questo approccio funziona bene per usare un piccolo numero di campi dalla meta map, ma si adatterebbe male se volessimo usarne molti di più.

Esiste un altro operatore chiamato `multiMap()` che ci permette di semplificare un po' questo, ma anche così non è ideale.

??? info "(Opzionale) Versione alternativa con `multiMap()`"

    Nel caso ve lo stiate chiedendo, non potevamo semplicemente scrivere una singola operazione `map()` che restituisce sia il `file` che il `character`, perché li restituirebbe come tupla.
    Dovevamo scrivere due operazioni `map()` separate per alimentare gli elementi `file` e `character` al processo separatamente.

    Tecnicamente c'è un altro modo per farlo attraverso una singola operazione di mapping, usando l'operatore `multiMap()`, che è in grado di emettere più canali.
    Ad esempio, potreste sostituire la chiamata a `COWPY` sopra con il seguente codice:

    === "Dopo"

        ```groovy title="main.nf" linenums="34"
            // Esegue cowpy per generare ASCII art
            COWPY(
                ch_languages.multiMap { meta, file ->
                    file: file
                    character: meta.character
                }
            )
        ```

    === "Prima"

        ```groovy title="main.nf" linenums="34"
            // Esegue cowpy per generare ASCII art
            COWPY(
                ch_languages.map { meta, file -> file },
                ch_languages.map { meta, file -> meta.character }
            )
        ```

    Questo produce esattamente lo stesso risultato.

In entrambi i casi, è scomodo dover fare un po' di estrazione a livello del flusso di lavoro.

Sarebbe meglio se potessimo passare l'intera meta map al processo e scegliere ciò di cui abbiamo bisogno una volta lì.

### 3.3. Passare e usare l'intera meta map

Il punto della meta map è dopotutto passare tutti i metadati insieme come un pacchetto.
L'unico motivo per cui non potevamo farlo sopra è che il processo non è configurato per accettare una meta map.
Ma poiché controlliamo il codice del processo, possiamo cambiarlo.

Modifichiamo il processo `COWPY` per accettare la struttura tupla `[meta, file]` che abbiamo usato nel primo processo in modo da poter semplificare il flusso di lavoro.

A tal fine, dovremo fare tre cose:

1. Modificare le definizioni di input del modulo del processo `COWPY`
2. Aggiornare il comando del processo per usare la meta map
3. Aggiornare la chiamata al processo nel corpo del flusso di lavoro

Pronti? Andiamo!

#### 3.3.1. Modificare l'input del modulo `COWPY`

Apportate le seguenti modifiche al file del modulo `cowpy.nf`:

=== "Dopo"

    ```groovy title="cowpy.nf" linenums="10" hl_lines="2"
    input:
    tuple val(meta), path(input_file)
    ```

=== "Prima"

    ```groovy title="cowpy.nf" linenums="10" hl_lines="2-3"
    input:
    path(input_file)
    val character
    ```

Questo ci permette di usare la struttura tupla `[meta, file]` che abbiamo trattato in precedenza nel tutorial.

Notate che non abbiamo aggiornato la definizione di output del processo per restituire la meta map, al fine di mantenere il tutorial snello, ma sentitevi liberi di farlo voi stessi come esercizio seguendo il modello del processo `IDENTIFY_LANGUAGE`.

#### 3.3.2. Aggiornare il comando per usare il campo della meta map

L'intera meta map è ora disponibile all'interno del processo, quindi possiamo fare riferimento alle informazioni che contiene direttamente dall'interno del blocco del comando.

Apportate le seguenti modifiche al file del modulo `cowpy.nf`:

=== "Dopo"

    ```groovy title="cowpy.nf" linenums="16" hl_lines="3"
    script:
    """
    cat ${input_file} | cowpy -c ${meta.character} > cowpy-${input_file}
    """
    ```

=== "Prima"

    ```groovy title="cowpy.nf" linenums="16" hl_lines="3"
    script:
    """
    cat ${input_file} | cowpy -c ${character} > cowpy-${input_file}
    """
    ```

Abbiamo sostituito il riferimento al valore `character` precedentemente passato come input autonomo con il valore contenuto nella meta map, a cui facciamo riferimento usando `meta.character`.

Ora aggiorniamo di conseguenza la chiamata al processo.

#### 3.3.3. Aggiornare la chiamata al processo ed eseguirlo

Il processo ora si aspetta che il suo input usi la struttura tupla `[meta, file]`, che è ciò che il processo precedente restituisce, quindi possiamo semplicemente passare l'intero canale `ch_languages` al processo `COWPY`.

Apportate le seguenti modifiche al flusso di lavoro principale:

=== "Dopo"

    ```groovy title="main.nf" linenums="34" hl_lines="2"
    // Esegue cowpy per generare ASCII art
    COWPY(ch_languages)
    ```

=== "Prima"

    ```groovy title="main.nf" linenums="34" hl_lines="3-4"
    // Esegue cowpy per generare ASCII art
    COWPY(
        ch_languages.map { meta, file -> file },
        ch_languages.map { meta, file -> meta.character }
    )
    ```

Questo semplifica notevolmente la chiamata!

Eliminiamo i risultati dell'esecuzione precedente ed eseguiamo:

```bash
rm -r results
nextflow run main.nf
```

??? success "Output del comando"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [wise_sammet] DSL2 - revision: 99797b1e92

    executor >  local (14)
    [5d/dffd4e] process > IDENTIFY_LANGUAGE (7) [100%] 7 of 7 ✔
    [25/9243df] process > COWPY (7)             [100%] 7 of 7 ✔
    ```

Se guardate nella directory dei risultati, dovreste vedere gli stessi output di prima, _cioè_ singoli file contenenti l'ASCII art di ogni saluto pronunciato dal personaggio corrispondente.

??? abstract "Contenuto della directory"

    ```console
    ./results/
    ├── cowpy-bonjour.txt
    ├── cowpy-ciao.txt
    ├── cowpy-guten_tag.txt
    ├── cowpy-hallo.txt
    ├── cowpy-hello.txt
    ├── cowpy-hola.txt
    └── cowpy-salut.txt
    ```

Quindi questo produce gli stessi risultati di prima con un codice più semplice.

Naturalmente, questo presuppone che siate in grado di modificare il codice del processo.
In alcuni casi, potreste dover fare affidamento su processi esistenti che non siete liberi di modificare, il che limita le vostre opzioni.
La buona notizia, se state pianificando di usare moduli dal progetto [nf-core](https://nf-co.re/), è che i moduli nf-core sono tutti configurati per usare la struttura tupla `[meta, file]` come standard.

### 3.4. Risoluzione dei problemi relativi agli input obbligatori mancanti

Il valore `character` è necessario affinché il processo `COWPY` venga eseguito con successo.
Se non impostiamo un valore predefinito per esso in un file di configurazione, DOBBIAMO fornire un valore per esso nel foglio dati.

**Cosa succede se non lo facciamo?**
Dipende da cosa contiene il foglio dati di input e quale versione del flusso di lavoro stiamo eseguendo.

#### 3.4.1. La colonna character esiste ma è vuota

Supponiamo di eliminare il valore del personaggio per una delle voci nel nostro foglio dati per simulare un errore di raccolta dati:

```csv title="datasheet.csv" linenums="1" hl_lines="2"
id,character,recording
sampleA,,/workspaces/training/side-quests/metadata/data/bonjour.txt
sampleB,tux,/workspaces/training/side-quests/metadata/data/guten_tag.txt
sampleC,sheep,/workspaces/training/side-quests/metadata/data/hallo.txt
sampleD,turkey,/workspaces/training/side-quests/metadata/data/hello.txt
sampleE,stegosaurus,/workspaces/training/side-quests/metadata/data/hola.txt
sampleF,moose,/workspaces/training/side-quests/metadata/data/salut.txt
sampleG,turtle,/workspaces/training/side-quests/metadata/data/ciao.txt
```

Per entrambe le versioni del flusso di lavoro che abbiamo usato sopra, la chiave `character` verrà creata per tutte le voci quando il foglio dati viene letto, ma per `sampleA` il valore sarà una stringa vuota.

Questo causerà un errore.

??? failure "Output del comando"

    ```console hl_lines="8 11 16 28"
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [marvelous_hirsch] DSL2 - revision: 0dfeee3cc1

    executor >  local (9)
    [c1/c5dd4f] process > IDENTIFY_LANGUAGE (7) [ 85%] 6 of 7
    [d3/b7c415] process > COWPY (2)             [  0%] 0 of 6
    ERROR ~ Error executing process > 'COWPY (1)'

    Caused by:
      Process `COWPY (1)` terminated with an error exit status (2)


    Command executed:

      cat bonjour.txt | cowpy -c  > cowpy-bonjour.txt

    Command exit status:
      2

    Command output:
      (empty)

    Command error:
      usage: cowpy [-h] [-l] [-L] [-t] [-u] [-e EYES] [-c COWACTER] [-E] [-r] [-x]
                  [-C]
                  [msg ...]
      cowpy: error: argument -c/--cowacter: expected one argument

    Work dir:
      /workspaces/training/side-quests/metadata/work/ca/9d49796612a54dec5ed466063c809b

    Container:
      community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273

    Tip: you can try to figure out what's wrong by changing to the process work dir and showing the script file named `.command.sh`

    -- Check '.nextflow.log' file for details
    ```

Quando Nextflow esegue la riga di comando `cowpy` per quel campione, `${meta.character}` viene riempito con una stringa vuota nella riga di comando `cowpy`, quindi lo strumento `cowpy` genera un errore dicendo che non è stato fornito alcun valore per l'argomento `-c`.

#### 3.4.2. La colonna character non esiste nel foglio dati

Ora supponiamo di eliminare completamente la colonna `character` dal nostro foglio dati:

```csv title="datasheet.csv" linenums="1"
id,recording
sampleA,/workspaces/training/side-quests/metadata/data/bonjour.txt
sampleB,/workspaces/training/side-quests/metadata/data/guten_tag.txt
sampleC,/workspaces/training/side-quests/metadata/data/hallo.txt
sampleD,/workspaces/training/side-quests/metadata/data/hello.txt
sampleE,/workspaces/training/side-quests/metadata/data/hola.txt
sampleF,/workspaces/training/side-quests/metadata/data/salut.txt
sampleG,/workspaces/training/side-quests/metadata/data/ciao.txt
```

In questo caso la chiave `character` non verrà creata affatto quando il foglio dati viene letto.

##### 3.4.2.1. Valore acceduto a livello del flusso di lavoro

Se stiamo usando la versione del codice che abbiamo scritto nella sezione 3.2, Nextflow tenterà di accedere alla chiave `character` nella meta map PRIMA di chiamare il processo `COWPY`.

Non troverà alcun elemento che corrisponda all'istruzione, quindi non eseguirà `COWPY` affatto.

??? success "Output del comando"

    ```console hl_lines="7"
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [desperate_montalcini] DSL2 - revision: 0dfeee3cc1

    executor >  local (7)
    [1a/df2544] process > IDENTIFY_LANGUAGE (7) [100%] 7 of 7 ✔
    [-        ] process > COWPY                 -
    ```

Per quanto riguarda Nextflow, questo flusso di lavoro è stato eseguito con successo!
Tuttavia, nessuno degli output che vogliamo verrà prodotto.

##### 3.4.2.2. Valore acceduto a livello del processo

Se stiamo usando la versione nella sezione 3.3, Nextflow passerà l'intera meta map al processo `COWPY` e tenterà di eseguire il comando.

Questo causerà un errore, ma diverso rispetto al primo caso.

??? failure "Output del comando"

    ```console hl_lines="8 11 16"
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [jovial_bohr] DSL2 - revision: eaaf375827

    executor >  local (9)
    [0d/ada9db] process > IDENTIFY_LANGUAGE (5) [ 85%] 6 of 7
    [06/28065f] process > COWPY (2)             [  0%] 0 of 6
    ERROR ~ Error executing process > 'COWPY (2)'

    Caused by:
      Process `COWPY (2)` terminated with an error exit status (1)


    Command executed:

      cat guten_tag.txt | cowpy -c null > cowpy-guten_tag.txt

    Command exit status:
      1

    Command output:
      (empty)

    Command error:
      Traceback (most recent call last):
        File "/opt/conda/bin/cowpy", line 10, in <module>
          sys.exit(main())
                  ~~~~^^
        File "/opt/conda/lib/python3.13/site-packages/cowpy/cow.py", line 1215, in main
          print(cow(eyes=args.eyes,
                ~~~^^^^^^^^^^^^^^^^
                tongue=args.tongue,
                ^^^^^^^^^^^^^^^^^^^
                thoughts=args.thoughts
                ^^^^^^^^^^^^^^^^^^^^^^
                    ).milk(msg)
                    ^
      TypeError: 'str' object is not callable

    Work dir:
      /workspaces/training/side-quests/metadata/work/06/28065f7d9fd7d22bba084aa941b6d6

    Container:
      community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273

    Tip: you can replicate the issue by changing to the process work dir and entering the command `bash .command.run`

    -- Check '.nextflow.log' file for details
    ```

Questo accade perché `meta.character` non esiste, quindi il nostro tentativo di accedervi restituisce `null`. Di conseguenza, Nextflow inserisce letteralmente `null` nella riga di comando, che ovviamente non viene riconosciuto dallo strumento `cowpy`.

#### 3.4.3. Soluzioni

Oltre a fornire un valore predefinito come parte della configurazione del flusso di lavoro, ci sono due cose che possiamo fare per gestire questo in modo più robusto:

1. Implementare la validazione dell'input nel vostro flusso di lavoro per garantire che il foglio dati contenga tutte le informazioni richieste. Potete trovare un'[introduzione alla validazione dell'input](../hello_nf-core/05_input_validation.md) nel corso di formazione Hello nf-core. <!-- TODO (future) pending a proper Validation side quest -->

2. Se volete assicurarvi che chiunque usi il vostro modulo di processo possa identificare immediatamente gli input obbligatori, potete anche rendere la proprietà dei metadati richiesta un input esplicito.

Ecco un esempio di come funzionerebbe.

Prima, a livello del processo, aggiornate la definizione di input come segue:

=== "Dopo"

    ```groovy title="cowpy.nf" linenums="12" hl_lines="2"
        input:
        tuple val(meta), val(character), path(input_file)
    ```

=== "Prima"

    ```groovy title="cowpy.nf" linenums="12" hl_lines="2"
        input:
        tuple val(meta), path(input_file)
    ```

Poi, a livello del flusso di lavoro, usate un'operazione di mapping per estrarre la proprietà `character` dai metadati e renderla un componente esplicito della tupla di input:

=== "Dopo"

    ```groovy title="main.nf" linenums="37" hl_lines="1"
        COWPY(ch_languages.map{meta, file -> [meta, meta.character, file]})
    ```

=== "Prima"

    ```groovy title="main.nf" linenums="37" hl_lines="1"
        COWPY(ch_languages)
    ```

Questo approccio ha il vantaggio di mostrare esplicitamente che `character` è obbligatorio, e rende il processo più facile da ridistribuire in altri contesti.

Questo evidenzia un importante principio di progettazione:

**Usate la meta map per informazioni opzionali e descrittive, ma estraete i valori obbligatori come input espliciti.**

La meta map è eccellente per mantenere le strutture dei canali pulite e prevenire strutture di canali arbitrarie, ma per gli elementi obbligatori che sono direttamente referenziati in un processo, estrarli come input espliciti crea un codice più robusto e manutenibile.

### Takeaway

In questa sezione, avete imparato come utilizzare i metadati per personalizzare l'esecuzione di un processo, accedendovi sia a livello del flusso di lavoro che a livello del processo.

---

## Esercizio supplementare

Se volete esercitarvi a usare le informazioni della meta map dall'interno di un processo, provate a usare altri elementi della meta map come `lang` e `lang_group` per personalizzare come vengono nominati e/o organizzati gli output.

Ad esempio, provate a modificare il codice per produrre questo risultato:

```console title="Results directory contents"
results/
├── germanic
│   ├── de-guten_tag.txt
│   ├── de-hallo.txt
│   └── en-hello.txt
└── romance
    ├── es-hola.txt
    ├── fr-bonjour.txt
    ├── fr-salut.txt
    └── it-ciao.txt
```

<!-- TODO (future) Provide worked out solution -->
<!-- the renaming should use the meta inside the process -->
<!-- the output org should use the meta in the workflow outputs -->

---

## Riepilogo

In questa side quest, avete esplorato come lavorare efficacemente con i metadati nei flussi di lavoro Nextflow.

Questo schema di mantenere i metadati espliciti e allegati ai dati è una best practice fondamentale in Nextflow, che offre diversi vantaggi rispetto alla codifica fissa delle informazioni sui file:

- I metadati dei file rimangono associati ai file durante tutto il flusso di lavoro
- Il comportamento dei processi può essere personalizzato per file
- L'organizzazione dell'output può riflettere i metadati dei file
- Le informazioni sui file possono essere espanse durante l'esecuzione del pipeline

Applicare questo schema nel vostro lavoro vi permetterà di costruire flussi di lavoro bioinformatici robusti e manutenibili.

### Schemi chiave

1.  **Lettura e strutturazione dei metadati:** Leggere file CSV e creare meta map organizzate che rimangono associate ai vostri file di dati.

    ```groovy
    channel.fromPath('datasheet.csv')
      .splitCsv(header: true)
      .map { row ->
          [ [id:row.id, character:row.character], row.recording ]
      }
    ```

2.  **Espansione dei metadati durante il flusso di lavoro:** Aggiungere nuove informazioni ai vostri metadati man mano che il vostro pipeline progredisce, aggiungendo output dei processi e derivando valori attraverso la logica condizionale.

    - Aggiungere nuove chiavi basate sull'output del processo

    ```groovy
    .map { meta, file, lang ->
      [ meta + [lang:lang], file ]
    }
    ```

    - Aggiungere nuove chiavi usando una clausola condizionale

    ```groovy
    .map{ meta, file ->
        if ( meta.lang.equals("de") || meta.lang.equals('en') ){
            lang_group = "germanic"
        } else if ( meta.lang in ["fr", "es", "it"] ) {
            lang_group = "romance"
        } else {
            lang_group = "unknown"
        }
    }
    ```

3.  **Personalizzazione del comportamento dei processi:** Usare i metadati all'interno del processo.

    ```groovy
    cat $input_file | cowpy -c ${meta.character} > cowpy-${input_file}
    ```

### Risorse aggiuntive

- [map](https://www.nextflow.io/docs/latest/operator.html#map)
- [stdout](https://www.nextflow.io/docs/latest/process.html#outputs)

---

## Cosa c'è dopo?

Tornate al [menu delle Side Quest](../) o cliccate il pulsante in basso a destra della pagina per passare all'argomento successivo nell'elenco.
