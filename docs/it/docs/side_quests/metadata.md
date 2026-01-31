# Metadati e meta map

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduzione assistita da IA - [scopri di più e suggerisci miglioramenti](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

In qualsiasi analisi scientifica, raramente lavoriamo solo con i file di dati grezzi.
Ogni file è accompagnato da informazioni aggiuntive: cosa rappresenta, da dove proviene e cosa lo rende speciale.
Queste informazioni extra sono ciò che chiamiamo metadati.

I metadati sono dati che descrivono altri dati.
I metadati tengono traccia di dettagli importanti sui file e sulle condizioni sperimentali, e aiutano ad adattare le analisi alle caratteristiche uniche di ciascun dataset.

Pensate ad esso come al catalogo di una biblioteca: mentre i libri contengono il contenuto effettivo (dati grezzi), le schede del catalogo forniscono informazioni essenziali su ciascun libro—quando è stato pubblicato, chi lo ha scritto, dove trovarlo (metadati).
Nelle pipeline Nextflow, i metadati possono essere utilizzati per:

- Tracciare informazioni specifiche sui file durante tutto il workflow
- Configurare i processi in base alle caratteristiche dei file
- Raggruppare file correlati per analisi congiunte

### Obiettivi di apprendimento

In questa side quest, esploreremo come gestire i metadati nei workflow.
Partendo da un semplice datasheet (spesso chiamato samplesheet in bioinformatica) contenente informazioni di base sui file, imparerete come:

- Leggere e analizzare i metadati dei file da file CSV
- Creare e manipolare meta map
- Aggiungere nuovi campi di metadati durante l'esecuzione del workflow
- Utilizzare i metadati per personalizzare il comportamento dei processi

Queste competenze vi aiuteranno a costruire pipeline più robuste e flessibili in grado di gestire relazioni tra file complesse e requisiti di elaborazione.

### Prerequisiti

Prima di intraprendere questa side quest, dovrebbe:

- Aver completato il tutorial [Hello Nextflow](../hello_nextflow/README.md) o un corso equivalente per principianti.
- Essere a vostro agio nell'utilizzo di concetti e meccanismi di base di Nextflow (processi, canali, operatori)

---

## 0. Iniziare

#### Aprire il codespace di formazione

Se non lo ha ancora fatto, assicuratevi di aprire l'ambiente di formazione come descritto in [Configurazione dell'Ambiente](../envsetup/index.md).

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Spostarsi nella directory del progetto

Spostiamoci nella directory dove si trovano i file per questo tutorial.

```bash
cd side-quests/metadata
```

Può configurare VSCode per concentrarsi su questa directory:

```bash
code .
```

#### Esaminare i materiali

Troverà un file di workflow principale e una directory `data` contenente un datasheet e una manciata di file di dati.

??? abstract "Contenuti della directory"

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

Il workflow nel file `main.nf` è una bozza che espanderà gradualmente in un workflow completamente funzionante.

Il datasheet elenca i percorsi ai file di dati e alcuni metadati associati, organizzati in 3 colonne:

- `id`: autoesplicativo, un ID assegnato al file
- `character`: un nome di personaggio, che useremo più avanti per disegnare creature diverse
- `data`: percorsi a file `.txt` che contengono saluti in lingue diverse

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

Ogni file di dati contiene del testo di saluto in una di cinque lingue (fr: francese, de: tedesco, es: spagnolo, it: italiano, en: inglese).

Le forniremo anche uno strumento di analisi linguistica containerizzato chiamato `langid`.

#### Esaminare l'assegnazione

La vostra sfida è scrivere un workflow Nextflow che:

1. **Identifichi** automaticamente la lingua in ciascun file
2. **Raggruppi** i file per famiglia linguistica (lingue germaniche vs lingue romanze)
3. **Personalizzi** l'elaborazione per ciascun file in base alla sua lingua e ai suoi metadati
4. **Organizzi** gli output per gruppo linguistico

Questo rappresenta un pattern di workflow tipico in cui i metadati specifici dei file guidano le decisioni di elaborazione; esattamente il tipo di problema che le meta map risolvono in modo elegante.

#### Lista di controllo della preparazione

Pensa di essere pronto per iniziare?

- [ ] Comprendo l'obiettivo di questo corso e i suoi prerequisiti
- [ ] Il mio codespace è attivo e funzionante
- [ ] Ho impostato la mia directory di lavoro in modo appropriato
- [ ] Comprendo l'assegnazione

Se può spuntare tutte le caselle, è pronto per iniziare.

---

## 1. Caricare metadati da un datasheet

Aprite il file di workflow `main.nf` per esaminare la bozza di workflow che Le forniamo come punto di partenza.

```groovy title="main.nf" linenums="1"
#!/usr/bin/env nextflow

workflow  {

    ch_datasheet = channel.fromPath("./data/datasheet.csv")

}
```

Può vedere che abbiamo configurato un channel factory di base per caricare il datasheet di esempio come file, ma questo non leggerà ancora il contenuto del file.
Iniziamo aggiungendo questo.

### 1.1. Leggere il contenuto con `splitCsv`

Dobbiamo scegliere un operatore che analizzerà il contenuto del file in modo appropriato con il minimo sforzo da parte nostra.
Poiché il nostro datasheet è in formato CSV, questo è un lavoro per l'operatore [`splitCsv`](https://www.nextflow.io/docs/latest/reference/operator.html#splitcsv), che carica ogni riga nel file come elemento nel canale.

Apporti le seguenti modifiche per aggiungere un'operazione `splitCsv()` al codice di costruzione del canale, più un'operazione `view()` per verificare che il contenuto del file venga caricato correttamente nel canale.

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

Noti che stiamo usando l'opzione `header: true` per dire a Nextflow di leggere la prima riga del file CSV come riga di intestazione.

Vediamo cosa ne esce, va bene?
Esegua il workflow:

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

Possiamo vedere che l'operatore ha costruito una mappa di coppie chiave-valore per ogni riga nel file CSV, con le intestazioni delle colonne come chiavi per i valori corrispondenti.

Ogni voce della mappa corrisponde a una colonna nel nostro datasheet:

- `id`
- `character`
- `recording`

Questo è ottimo! Rende facile accedere a campi specifici da ciascun file.
Ad esempio, potremmo accedere all'ID del file con `id` o al percorso del file txt con `recording`.

??? info "(Opzionale) Maggiori informazioni sulle map"

    In Groovy, il linguaggio di programmazione su cui è costruito Nextflow, una map è una struttura dati chiave-valore simile ai dizionari in Python, agli oggetti in JavaScript o agli hash in Ruby.

    Ecco uno script eseguibile che mostra come è possibile definire una map e accedere ai suoi contenuti nella pratica:

    ```groovy title="examples/map_demo.nf"
    #!/usr/bin/env nextflow

    // Crea una map semplice
    def my_map = [id:'sampleA', character:'squirrel']

    // Stampa l'intera map
    println "map: ${my_map}"

    // Accede ai singoli valori usando la notazione con punto
    println "id: ${my_map.id}"
    println "character: ${my_map.character}"
    ```

    Anche se non ha un blocco `workflow` appropriato, Nextflow può eseguirlo come se fosse un workflow:

    ```bash
    nextflow run examples/map_demo.nf
    ```

    Ed ecco cosa può aspettarsi di vedere nell'output:

    ```console title="Output"
     N E X T F L O W   ~  version 25.10.2

    Launching `map_demo.nf` [cheesy_plateau] DSL2 - revision: fae5b8496e

    map: [id:sampleA, character:squirrel]
    id: sampleA
    character: squirrel
    ```

### 1.2. Selezionare campi specifici con `map`

Diciamo che vogliamo accedere alla colonna `character` dal datasheet e stamparla.
Possiamo usare l'operatore Nextflow `map` per iterare su ogni elemento nel nostro canale e selezionare specificamente la voce `character` dall'oggetto map.

Apporti le seguenti modifiche al workflow:

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

Ora eseguite nuovamente il workflow:

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

Successo! Abbiamo sfruttato la struttura della map derivata dal nostro datasheet per accedere ai valori delle singole colonne per ogni riga.

Ora che abbiamo letto con successo il datasheet e abbiamo accesso ai dati in ogni riga, possiamo iniziare a implementare la logica della nostra pipeline.

### 1.3. Organizzare i metadati in una 'meta map'

Nello stato attuale del workflow, i file di input (sotto la chiave `recording`) e i metadati associati (`id`, `character`) sono tutti sullo stesso piano, come se fossero tutti in un grande contenitore.
La conseguenza pratica è che ogni processo che consuma questo canale dovrebbe essere configurato con questa struttura in mente:

```groovy
    input:
    tuple val(id), val(character), file(recording)
```

Va bene finché il numero di colonne nel datasheet non cambia.
Tuttavia, se aggiunge anche solo una colonna al datasheet, la forma del canale non corrisponderà più a ciò che il processo si aspetta, e il workflow produrrà errori.
Rende anche difficile condividere il processo con altri che potrebbero avere dati di input leggermente diversi, e potrebbe finire per dover codificare variabili nel processo che non sono necessarie per il blocco script.

Per evitare questo problema, dobbiamo trovare un modo per mantenere la struttura del canale coerente indipendentemente da quante colonne contiene quel datasheet.

Possiamo farlo raccogliendo tutti i metadati in un elemento all'interno della tupla, che chiameremo meta map, o più semplicemente 'meta map'.

Apporti le seguenti modifiche all'operazione `map`:

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

Abbiamo ristrutturato gli elementi del nostro canale in una tupla composta da due elementi, la meta map e l'oggetto file corrispondente.

Eseguiamo il workflow:

```bash
nextflow run main.nf
```

??? success "Output del comando"

    ```console title="Visualizza meta map"
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

```console title="Esempio di struttura di output"
[
  [id:sampleA, character:squirrel],
  /workspaces/training/side-quests/metadata/data/bonjour.txt
]
```

Di conseguenza, l'aggiunta di più colonne nel datasheet renderà disponibili più metadati nella map `meta`, ma non cambierà la forma del canale.
Questo ci consente di scrivere processi che consumano il canale senza dover codificare gli elementi dei metadati nella specifica di input:

```groovy title="Esempio di sintassi"
    input:
    tuple val(meta), file(recording)
```

Questa è una convenzione ampiamente utilizzata per organizzare i metadati nei workflow Nextflow.

### Concetti chiave

In questa sezione, ha imparato:

- **Perché i metadati sono importanti:** Mantenere i metadati con i vostri dati preserva informazioni importanti sui file durante tutto il workflow.
- **Come leggere i datasheet:** Utilizzare `splitCsv` per leggere file CSV con informazioni di intestazione e trasformare le righe in dati strutturati
- **Come creare una meta map:** Separare i metadati dai dati dei file utilizzando la struttura di tupla `[ [id:value, ...], file ]`

---

## 2. Manipolare i metadati

Ora che abbiamo caricato i nostri metadati, facciamo qualcosa con essi!

Useremo uno strumento chiamato [`langid`](https://github.com/saffsd/langid.py) per identificare la lingua contenuta nel file di registrazione di ogni creatura.
Lo strumento è pre-addestrato su un insieme di lingue e, dato un frammento di testo, produrrà una previsione linguistica e un punteggio di probabilità associato, entrambi su `stdout`.

### 2.1. Importare il processo ed esaminare il codice

Le forniamo un modulo di processo pre-scritto chiamato `IDENTIFY_LANGUAGE` che avvolge lo strumento `langid`, quindi deve solo aggiungere un'istruzione include prima del blocco workflow.

Apporti la seguente modifica al workflow:

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

Può aprire il file del modulo per esaminarne il codice:

```groovy title="modules/langid.nf" linenums="1" hl_lines="9 12"
#!/usr/bin/env nextflow

// Usa langid per prevedere la lingua di ciascun file di input
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

Come potete vedere, la definizione di input utilizza la stessa struttura `tuple val(meta), path(file)` che abbiamo appena applicato al nostro canale di input.

La definizione di output è strutturata come una tupla con una struttura simile a quella dell'input, tranne per il fatto che contiene anche `stdout` come terzo elemento.
Questo pattern `tuple val(meta), path(file), <output>` mantiene i metadati associati sia ai dati di input che agli output mentre scorrono attraverso la pipeline.

Noti che stiamo usando il qualificatore di output [`stdout`](https://www.nextflow.io/docs/latest/process.html#outputs) di Nextflow qui perché lo strumento stampa il suo output direttamente sulla console piuttosto che scrivere un file; e usiamo `sed` nella riga di comando per rimuovere il punteggio di probabilità, pulire la stringa rimuovendo i caratteri di nuova riga e restituire solo la previsione linguistica.

### 2.2. Aggiungere una chiamata a `IDENTIFY_LANGUAGE`

Ora che il processo è disponibile per il workflow, possiamo aggiungere una chiamata al processo `IDENTIFY_LANGUAGE` per eseguirlo sul canale dati.

Apporti le seguenti modifiche al workflow:

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

Noti che abbiamo rimosso l'operazione `.view()` originale nella costruzione del canale.

Possiamo ora eseguire il workflow:

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

Eccellente! Ora abbiamo una previsione per quale lingua parla ogni personaggio.

E come notato in precedenza, abbiamo anche incluso il file di input e la meta map nell'output, il che significa che entrambi rimangono associati alle nuove informazioni che abbiamo appena prodotto.
Questo si rivelerà utile nel prossimo passaggio.

!!! note

    Più in generale, questo pattern di mantenere la meta map associata ai risultati rende più facile associare risultati correlati che condividono gli stessi identificatori.

    Come avrà già imparato, non potete fare affidamento sull'ordine degli elementi nei canali per abbinare i risultati tra di essi.
    Invece, deve usare chiavi per associare correttamente i dati, e le meta map forniscono una struttura ideale per questo scopo.

    Esploriamo questo caso d'uso in dettaglio nella side quest [Splitting & Grouping](./splitting_and_grouping.md).

### 2.3. Aumentare i metadati con gli output dei processi

Dato che i risultati che abbiamo appena prodotto sono di per sé una forma di metadati sui contenuti dei file, sarebbe utile aggiungerli alla nostra meta map.

Tuttavia, non vogliamo modificare la meta map esistente sul posto.
Da un punto di vista tecnico, è _possibile_ farlo, ma non è sicuro.

Quindi invece, creeremo una nuova meta map contenente il contenuto della meta map esistente più una nuova coppia chiave-valore `lang: lang_id` che contiene le nuove informazioni, utilizzando l'operatore `+` (una caratteristica di Groovy).
E combineremo questo con un'operazione [`map`](https://www.nextflow.io/docs/latest/operator.html#map) per sostituire la vecchia map con quella nuova.

Ecco le modifiche che deve apportare al workflow:

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

Se non ha ancora familiarità con l'operatore `+`, o se questo sembra confuso, si prenda qualche minuto per esaminare la spiegazione dettagliata qui sotto.

??? info "Creazione della nuova meta map usando l'operatore `+`"

    **Prima di tutto, deve sapere che possiamo unire i contenuti di due map usando l'operatore Groovy `+`.**

    Diciamo che abbiamo le seguenti map:

    ```groovy
    map1 = [id: 'sampleA', character: 'squirrel']
    map2 = [lang: 'fr']
    ```

    Possiamo unirle in questo modo:

    ```groovy
    new_map = map1 + map2
    ```

    Il contenuto di `new_map` sarà:

    ```groovy
    [id: 'sampleA', character: 'squirrel', lang: 'fr']
    ```

    Ottimo!

    **Ma cosa succede se deve aggiungere un campo che non fa già parte di una map?**

    Diciamo che riparte da `map1`, ma la previsione linguistica non è nella sua map (non c'è `map2`).
    Invece, è contenuta in una variabile chiamata `lang_id`, e sa che vuole memorizzare il suo valore (`'fr'`) con la chiave `lang`.

    Può effettivamente fare quanto segue:

    ```groovy
    new_map = [map1 + [lang: lang_id]]
    ```

    Qui, `[lang: new_info]` crea una nuova map senza nome al volo, e `map1 + ` unisce `map1` con la nuova map senza nome, producendo gli stessi contenuti di `new_map` come prima.

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
    - `[map1 + [lang: lang_id]]` crea la nuova map come dettagliato sopra

    L'output è una singola map senza nome con gli stessi contenuti di `new_map` nel nostro esempio sopra.
    Quindi abbiamo effettivamente trasformato:

    ```groovy
    [id: 'sampleA', character: 'squirrel'], 'it'
    ```

    in:

    ```groovy
    [id: 'sampleA', character: 'squirrel', lang: 'fr']
    ```

    Si spera che possa vedere che se cambiamo `map1` in `meta`, è fondamentalmente tutto ciò di cui abbiamo bisogno per aggiungere la previsione linguistica alla nostra meta map nel nostro workflow.

    Tranne una cosa!

    Nel caso del nostro workflow, **dobbiamo anche tenere conto della presenza dell'oggetto `file` nella tupla**, che è composta da `meta, file, lang_id`.

    Quindi il codice qui diventerebbe:

    ```groovy
    .map { meta, file, lang_id ->
        [meta + [lang: lang_id], file]
    }
    ```

    Se avete difficoltà a capire perché il `file` sembra muoversi nell'operazione `map`, immagini che invece di `[meta + [lang: lang_id], file]`, quella riga legga `[new_map, file]`.
    Questo dovrebbe rendere più chiaro che stiamo semplicemente lasciando il `file` nella sua posizione originale in seconda posizione nella tupla. Abbiamo appena preso il valore `new_info` e lo abbiamo incorporato nella map che è in prima posizione.

    **E questo ci riporta alla struttura del canale `tuple val(meta), path(file)`!**

Una volta che siete sicuri di aver capito cosa fa questo codice, eseguite il workflow per vedere se ha funzionato:

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

Sì, è corretto!
Abbiamo riorganizzato ordinatamente l'output del processo da `meta, file, lang_id` in modo che `lang_id` sia ora una delle chiavi nella meta map, e le tuple del canale si adattino nuovamente al modello `meta, file`.

<!-- TODO (future) Should we also show how to remove a key using subMap?! Or note where to find that. -->

### 2.4. Assegnare un gruppo linguistico usando i condizionali

Ora che abbiamo le nostre previsioni linguistiche, usiamo le informazioni per assegnare alcuni nuovi raggruppamenti.

Nei nostri dati di esempio, le lingue utilizzate dai nostri personaggi possono essere raggruppate in lingue germaniche (inglese, tedesco) e lingue romanze (francese, spagnolo, italiano).
Potrebbe essere utile avere quella classificazione facilmente disponibile da qualche parte più avanti nella pipeline, quindi aggiungiamo quelle informazioni nella meta map.

E, buone notizie, questo è ancora un altro caso che si presta perfettamente all'uso dell'operatore `map`!

Nello specifico, definiremo una variabile chiamata `lang_group`, useremo una logica condizionale semplice per determinare quale valore assegnare a `lang_group` per ciascun pezzo di dati.

La sintassi generale sarà così:

```groovy
.map { meta, file ->

    // la logica condizionale che definisce lang_group va qui

    [meta + [lang_group: lang_group], file]
}
```

Può vedere che questo è molto simile all'operazione di fusione di map al volo che abbiamo usato nel passaggio precedente.
Dobbiamo solo scrivere le istruzioni condizionali.

Ecco la logica condizionale che vogliamo applicare:

- Definire una variabile chiamata `lang_group` con valore predefinito `'unknown'`.
- Se `lang` è tedesco (`'de'`) o inglese (`'en'`), cambiare `lang_group` in `germanic`.
- Altrimenti se `lang` è incluso in un elenco contenente francese (`'fr'`), spagnolo (`'es'`) e italiano (`'it'`), cambiare `lang_group` in `romance`.

Provi a scriverlo voi stessi se sa già come scrivere istruzioni condizionali in Nextflow.

!!! tip

    Può accedere al valore di `lang` all'interno dell'operazione map con `meta.lang`.

Dovrebbe finire per apportare le seguenti modifiche al workflow:

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
- Usiamo l'operazione di fusione `meta + [lang_group:lang_group]` come precedentemente per generare la meta map aggiornata.

<!-- TODO (future) Add note/links to relevant docs in additional resources section -->

Una volta che tutto ha senso, eseguite nuovamente il workflow per vedere il risultato:

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

### Concetti chiave

In questa sezione, ha imparato come:

- **Applicare i metadati di input ai canali di output**: Copiare i metadati in questo modo ci consente di associare i risultati successivamente in base al contenuto dei metadati.
- **Creare chiavi personalizzate**: Ha creato due nuove chiavi nella vostra meta map, fondendole con `meta + [new_key:value]` nella meta map esistente. Una basata su un valore calcolato da un processo e una basata su una condizione impostata nell'operatore `map`.

Questi vi consentono di associare metadati nuovi ed esistenti con i file man mano che procede attraverso la vostra pipeline.
Anche se non sta utilizzando i metadati come parte di un processo, mantenere la meta map associata ai dati in questo modo rende facile tenere insieme tutte le informazioni rilevanti.

---

## 3. Utilizzare le informazioni della meta map in un processo

Ora che sa come creare e aggiornare la meta map, possiamo arrivare alla parte davvero divertente: usare effettivamente i metadati in un processo.

Più specificamente, aggiungeremo un secondo passaggio al nostro workflow per disegnare ogni animale come arte ASCII e fargli dire il testo registrato in un fumetto.
Lo faremo utilizzando uno strumento chiamato [`cowpy`](https://github.com/jeffbuttars/cowpy).

??? info "Cosa fa `cowpy`?"

    `cowpy` è uno strumento da riga di comando che genera arte ASCII per visualizzare input di testo arbitrari in modo divertente.
    È un'implementazione python del classico strumento [cowsay](https://en.wikipedia.org/wiki/Cowsay) di Tony Monroe.

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

    Opzionalmente, può selezionare un personaggio (o 'cowacter') da utilizzare al posto della mucca predefinita.

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

Se ha seguito il corso Hello Nextflow, ha già visto questo strumento in azione.
In caso contrario, non si preoccupi; tratteremo tutto ciò che deve sapere mentre procediamo.

### 3.1. Importare il processo ed esaminare il codice

Le forniamo un modulo di processo pre-scritto chiamato `COWPY` che avvolge lo strumento `cowpy`, quindi deve solo aggiungere un'istruzione include prima del blocco workflow.

Apporti la seguente modifica al workflow:

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

Può aprire il file del modulo per esaminarne il codice:

```groovy title="modules/cowpy.nf" linenums="1"
#!/usr/bin/env nextflow

// Genera arte ASCII con cowpy
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

Come potete vedere, questo processo è attualmente progettato per prendere un file di input (contenente il testo da visualizzare) e un valore che specifica il personaggio che dovrebbe essere disegnato in arte ASCII, solitamente fornito a livello di workflow da un parametro da riga di comando.

### 3.2. Passare un campo della meta map come input

Quando abbiamo usato lo strumento `cowpy` nel corso Hello Nextflow, abbiamo usato un parametro da riga di comando per determinare quale personaggio utilizzare per disegnare l'immagine finale.
Aveva senso, perché stavamo generando solo un'immagine per esecuzione della pipeline.

Tuttavia, in questo tutorial, vogliamo generare un'immagine appropriata per ogni soggetto che stiamo elaborando, quindi l'uso di un parametro da riga di comando sarebbe troppo limitante.

Buone notizie: abbiamo una colonna `character` nel nostro datasheet e quindi, nella nostra meta map.
Usiamo quella per impostare il personaggio che il processo dovrebbe usare per ogni voce.

A tal fine, dovremo fare tre cose:

1. Dare un nome al canale di output proveniente dal processo precedente in modo da potervi operare più comodamente.
2. Determinare come accedere alle informazioni di interesse
3. Aggiungere una chiamata al secondo processo e fornire le informazioni in modo appropriato.

Iniziamo.

#### 3.2.1. Nominare il canale di output precedente

Abbiamo applicato le manipolazioni precedenti direttamente sul canale di output del primo processo, `IDENTIFY_LANGUAGE.out`.
Per alimentare il contenuto di quel canale al processo successivo (e farlo in modo chiaro e facile da leggere) vogliamo dargli un nome proprio, `ch_languages`.

Possiamo farlo usando l'operatore [`set`](https://www.nextflow.io/docs/latest/reference/operator.html#set).

Nel workflow principale, sostituisca l'operatore `.view()` con `.set { ch_languages }`, e aggiunga una riga per testare che possiamo riferirci al canale per nome.

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

        // Temporaneo: sbirciare dentro ch_languages
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

Eseguiamo questo:

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

Questo conferma che ora possiamo riferirci al canale per nome.

#### 3.2.2. Accedere al file e ai metadati del personaggio

Sappiamo dall'esame del codice del modulo che il processo `COWPY` si aspetta di ricevere un file di testo e un valore `character`.
Per scrivere la chiamata al processo `COWPY`, dobbiamo solo sapere come estrarre il corrispondente oggetto file e i metadati da ciascun elemento nel canale.

Come spesso accade, il modo più semplice per farlo è usare un'operazione `map`.

Il nostro canale contiene tuple strutturate come `[meta, file]`, quindi possiamo accedere direttamente all'oggetto `file`, e possiamo accedere al valore `character` memorizzato all'interno della meta map riferendoci ad esso come `meta.character`.

Nel workflow principale, apporti le seguenti modifiche al codice:

=== "Dopo"

    ```groovy title="main.nf" linenums="34"
        // Temporaneo: accedere al file e al personaggio
        ch_languages.map { meta, file -> file }.view { file -> "File: " + file }
        ch_languages.map { meta, file -> meta.character }.view { character -> "Character: " + character }
    ```

=== "Prima"

    ```groovy title="main.nf" linenums="34"
        // Temporaneo: sbirciare dentro ch_languages
        ch_languages.view()
    ```

Noti che stiamo usando closure (come `{ file -> "File: " + file }`) per rendere l'output delle operazioni `.view` più leggibile.

Eseguiamo questo:

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

_I percorsi dei file e i valori dei personaggi potrebbero uscire in un ordine diverso nel vostro output._

Questo conferma che siamo in grado di accedere al file e al personaggio per ogni elemento nel canale.

#### 3.2.3. Chiamare il processo `COWPY`

Ora mettiamo tutto insieme e chiamiamo effettivamente il processo `COWPY` sul canale `ch_languages`.

Nel workflow principale, apporti le seguenti modifiche al codice:

=== "Dopo"

    ```groovy title="main.nf" linenums="34"
        // Esegue cowpy per generare arte ASCII
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

Vede che copiamo semplicemente le due operazioni map (meno le istruzioni `.view()`) come input alla chiamata del processo.
Assicuratevi solo di non dimenticare la virgola tra di loro!

È un po' goffo, ma vedremo come migliorarlo nella prossima sezione.

Eseguiamo questo:

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

Se guarda nella directory results, dovrebbe vedere i singoli file contenenti l'arte ASCII di ogni saluto pronunciato dal personaggio corrispondente.

??? abstract "Contenuti della directory e file di esempio"

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

Questo mostra che siamo stati in grado di utilizzare le informazioni nella meta map per parametrizzare il comando nel secondo passaggio della pipeline.

Tuttavia, come notato sopra, parte del codice coinvolto era un po' goffo, poiché abbiamo dovuto decomprimere i metadati mentre eravamo ancora nel contesto del corpo del workflow.
Questo approccio funziona bene per l'utilizzo di un piccolo numero di campi dalla meta map, ma scalrebbe male se volessimo usarne molti di più.

C'è un altro operatore chiamato `multiMap()` che ci consente di semplificare questo un po', ma anche allora non è ideale.

??? info "(Opzionale) Versione alternativa con `multiMap()`"

    Nel caso si stia chiedendo, non potevamo semplicemente scrivere una singola operazione `map()` che restituisse sia il `file` che il `character`, perché ciò li restituirebbe come tupla.
    Abbiamo dovuto scrivere due operazioni `map()` separate per alimentare gli elementi `file` e `character` al processo separatamente.

    Tecnicamente c'è un altro modo per farlo attraverso una singola operazione di mapping, usando l'operatore `multiMap()`, che è in grado di emettere più canali.
    Ad esempio, potrebbe sostituire la chiamata a `COWPY` sopra con il seguente codice:

    === "Dopo"

        ```groovy title="main.nf" linenums="34"
            // Esegue cowpy per generare arte ASCII
            COWPY(
                ch_languages.multiMap { meta, file ->
                    file: file
                    character: meta.character
                }
            )
        ```

    === "Prima"

        ```groovy title="main.nf" linenums="34"
            // Esegue cowpy per generare arte ASCII
            COWPY(
                ch_languages.map { meta, file -> file },
                ch_languages.map { meta, file -> meta.character }
            )
        ```

    Questo produce esattamente lo stesso risultato.

In entrambi i casi, è scomodo dover fare della decompressione a livello di workflow.

Sarebbe meglio se potessimo alimentare l'intera meta map nel processo e scegliere ciò di cui abbiamo bisogno una volta lì.

### 3.3. Passare e utilizzare l'intera meta map

Il punto della meta map è dopo tutto di passare tutti i metadati insieme come un pacchetto.
L'unica ragione per cui non potevamo farlo sopra è che il processo non è configurato per accettare una meta map.
Ma poiché controlliamo il codice del processo, possiamo cambiarlo.

Modifichiamo il processo `COWPY` per accettare la struttura di tupla `[meta, file]` che abbiamo usato nel primo processo in modo da poter semplificare il workflow.

A tal fine, dovremo fare tre cose:

1. Modificare le definizioni di input del modulo del processo `COWPY`
2. Aggiornare il comando del processo per utilizzare la meta map
3. Aggiornare la chiamata del processo nel corpo del workflow

Pronto? Andiamo!

#### 3.3.1. Modificare l'input del modulo `COWPY`

Apporti le seguenti modifiche al file del modulo `cowpy.nf`:

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

Questo ci consente di utilizzare la struttura di tupla `[meta, file]` che abbiamo trattato in precedenza nel tutorial.

Noti che non abbiamo aggiornato la definizione di output del processo per restituire la meta map, al fine di mantenere il tutorial snello, ma si senta libero di farlo voi stessi come esercizio seguendo il modello del processo `IDENTIFY_LANGUAGE`.

#### 3.3.2. Aggiornare il comando per utilizzare il campo della meta map

L'intera meta map è ora disponibile all'interno del processo, quindi possiamo riferirci alle informazioni che contiene direttamente dall'interno del blocco comando.

Apporti le seguenti modifiche al file del modulo `cowpy.nf`:

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

Abbiamo sostituito il riferimento al valore `character` precedentemente passato come input autonomo con il valore contenuto nella meta map, a cui ci riferiamo usando `meta.character`.

Ora aggiorniamo la chiamata del processo di conseguenza.

#### 3.3.3. Aggiornare la chiamata del processo ed eseguirlo

Il processo ora si aspetta che il suo input utilizzi la struttura di tupla `[meta, file]`, che è ciò che il processo precedente restituisce, quindi possiamo semplicemente passare l'intero canale `ch_languages` al processo `COWPY`.

Apporti le seguenti modifiche al workflow principale:

=== "Dopo"

    ```groovy title="main.nf" linenums="34" hl_lines="2"
    // Esegue cowpy per generare arte ASCII
