# Parte 3: Usare un modulo nf-core

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduzione assistita da IA - [ulteriori informazioni e suggerimenti per miglioramenti](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

In questa terza parte del corso di formazione Hello nf-core, vi mostriamo come trovare, installare e utilizzare un modulo nf-core esistente nella vostra pipeline.

Uno dei grandi vantaggi di lavorare con nf-core è la possibilità di sfruttare moduli pre-costruiti e testati dal repository [nf-core/modules](https://github.com/nf-core/modules).
Invece di scrivere ogni processo da zero, potete installare e utilizzare moduli mantenuti dalla comunità che seguono le best practice.

Per dimostrare come funziona, sostituiremo il modulo personalizzato `collectGreetings` con il modulo `cat/cat` da nf-core/modules nella pipeline `core-hello`.

??? info "Come iniziare da questa sezione"

    Questa sezione del corso presuppone che abbiate completato la [Parte 2: Riscrivere Hello per nf-core](./02_rewrite_hello.md) e che abbiate una pipeline `core-hello` funzionante.

    Se non avete completato la Parte 2 o volete ripartire da zero per questa parte, potete usare la soluzione `core-hello-part2` come punto di partenza.
    Eseguite questo comando dalla directory `hello-nf-core/`:

    ```bash
    cp -r solutions/core-hello-part2 core-hello
    cd core-hello
    ```

    Questo vi fornisce una pipeline nf-core completamente funzionale pronta per l'aggiunta di moduli.
    Potete verificare che funzioni correttamente eseguendo il seguente comando:

    ```bash
    nextflow run . --outdir core-hello-results -profile test,docker --validate_params false
    ```

---

## 1. Trovare e installare un modulo nf-core adatto

Per prima cosa, impariamo come trovare un modulo nf-core esistente e installarlo nella nostra pipeline.

Puntiamo a sostituire il processo `collectGreetings`, che usa il comando Unix `cat` per concatenare più file di saluti in uno solo.
Concatenare file è un'operazione molto comune, quindi è ragionevole pensare che possa già esistere un modulo in nf-core progettato per questo scopo.

Addentriamoci.

### 1.1. Esplorare i moduli disponibili sul sito web nf-core

Il progetto nf-core mantiene un catalogo centralizzato di moduli su [https://nf-co.re/modules](https://nf-co.re/modules).

Andate alla pagina dei moduli nel vostro browser web e usate la barra di ricerca per cercare 'concatenate'.

![risultati ricerca moduli](./img/module-search-results.png)

Come potete vedere, ci sono parecchi risultati, molti dei quali moduli progettati per concatenare tipi di file molto specifici.
Tra questi, dovreste vedere uno chiamato `cat_cat` che è generico.

!!! note "Convenzione di denominazione dei moduli"

    Il trattino basso (`_`) viene usato come sostituto del carattere slash (`/`) nei nomi dei moduli.

    I moduli nf-core seguono la convenzione di denominazione `software/comando` quando uno strumento fornisce più comandi, come `samtools/view` (pacchetto samtools, comando view) o `gatk/haplotypecaller` (pacchetto GATK, comando HaplotypeCaller).
    Per strumenti che forniscono un solo comando principale, i moduli usano un singolo livello come `fastqc` o `multiqc`.

Cliccate sul riquadro del modulo `cat_cat` per visualizzare la documentazione del modulo.

La pagina del modulo mostra:

- Una breve descrizione: "A module for concatenation of gzipped or uncompressed files"
- Comando di installazione: `nf-core modules install cat/cat`
- Struttura dei canali di input e output
- Parametri disponibili

### 1.2. Elencare i moduli disponibili dalla riga di comando

In alternativa, potete anche cercare moduli direttamente dalla riga di comando usando gli strumenti nf-core.

```bash
nf-core modules list remote
```

Questo mostrerà un elenco di tutti i moduli disponibili nel repository nf-core/modules, anche se è un po' meno comodo se non conoscete già il nome del modulo che state cercando.
Tuttavia, se lo conoscete, potete inviare l'elenco a `grep` per trovare moduli specifici:

```bash
nf-core modules list remote | grep 'cat/cat'
```

??? success "Output del comando"

    ```console
    │ cat/cat
    ```

Tenete presente che l'approccio con `grep` estrarrà solo risultati con il termine di ricerca nel loro nome, il che non funzionerebbe per `cat_cat`.

### 1.3. Ottenere informazioni dettagliate sul modulo

Per vedere informazioni dettagliate su un modulo specifico dalla riga di comando, usate il comando `info`:

```bash
nf-core modules info cat/cat
```

Questo mostra la documentazione sul modulo, inclusi i suoi input, output e informazioni di base sull'utilizzo.

??? success "Output del comando"

    ```console

                                              ,--./,-.
              ___     __   __   __   ___     /,-._.--~\
        |\ | |__  __ /  ` /  \ |__) |__         }  {
        | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                              `._,._,'

        nf-core/tools version 3.4.1 - https://nf-co.re


    ╭─ Module: cat/cat  ─────────────────────────────────────────────────╮
    │ 🌐 Repository: https://github.com/nf-core/modules.git              │
    │ 🔧 Tools: cat                                                      │
    │ 📖 Description: A module for concatenation of gzipped or           │
    │ uncompressed files                                                 │
    ╰────────────────────────────────────────────────────────────────────╯
                      ╷                                          ╷
    📥 Inputs        │Description                               │Pattern
    ╺━━━━━━━━━━━━━━━━━┿━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┿━━━━━━━╸
    input[0]         │                                          │
    ╶─────────────────┼──────────────────────────────────────────┼───────╴
      meta  (map)     │Groovy Map containing sample information  │
                      │e.g. [ id:'test', single_end:false ]      │
    ╶─────────────────┼──────────────────────────────────────────┼───────╴
      files_in  (file)│List of compressed / uncompressed files   │      *
                      ╵                                          ╵
                          ╷                                 ╷
    📥 Outputs           │Description                      │     Pattern
    ╺━━━━━━━━━━━━━━━━━━━━━┿━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┿━━━━━━━━━━━━╸
    file_out             │                                 │
    ╶─────────────────────┼─────────────────────────────────┼────────────╴
      meta  (map)         │Groovy Map containing sample     │
                          │information                      │
    ╶─────────────────────┼─────────────────────────────────┼────────────╴
      ${prefix}  (file)   │Concatenated file. Will be       │ ${file_out}
                          │gzipped if file_out ends with    │
                          │".gz"                            │
    ╶─────────────────────┼─────────────────────────────────┼────────────╴
    versions             │                                 │
    ╶─────────────────────┼─────────────────────────────────┼────────────╴
      versions.yml  (file)│File containing software versions│versions.yml
                          ╵                                 ╵

    💻  Installation command: nf-core modules install cat/cat

    ```

Queste sono esattamente le stesse informazioni che potete trovare sul sito web.

### 1.4. Installare il modulo cat/cat

Ora che abbiamo trovato il modulo che vogliamo, dobbiamo aggiungerlo al codice sorgente della nostra pipeline.

La buona notizia è che il progetto nf-core include alcuni strumenti per rendere questa parte facile.
Nello specifico, il comando `nf-core modules install` permette di automatizzare il recupero del codice e renderlo disponibile al vostro progetto in un singolo passaggio.

Andate nella directory della vostra pipeline ed eseguite il comando di installazione:

```bash
cd core-hello
nf-core modules install cat/cat
```

Lo strumento potrebbe prima chiedervi di specificare un tipo di repository.
(Se non lo fa, saltate alla sezione "Infine, lo strumento procederà all'installazione del modulo.")

??? success "Output del comando"

    ```console

                                          ,--./,-.
          ___     __   __   __   ___     /,-._.--~\
    |\ | |__  __ /  ` /  \ |__) |__         }  {
    | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                          `._,._,'

    nf-core/tools version 3.4.1 - https://nf-co.re


    WARNING  'repository_type' not defined in .nf-core.yml
    ? Is this repository a pipeline or a modules repository? (Use arrow keys)
    » Pipeline
      Modules repository
    ```

Se succede, premete invio per accettare la risposta predefinita (`Pipeline`) e continuare.

Lo strumento offrirà quindi di modificare la configurazione del vostro progetto per evitare questo prompt in futuro.

??? success "Output del comando"

    ```console
        INFO     To avoid this prompt in the future, add the 'repository_type' key to your .nf-core.yml file.
        ? Would you like me to add this config now? [y/n] (y):
    ```

Tanto vale approfittare di questo comodo strumento!
Premete invio per accettare la risposta predefinita (sì).

Infine, lo strumento procederà all'installazione del modulo.

??? success "Output del comando"

    ```console
    INFO Config added to '.nf-core.yml'
    INFO Reinstalling modules found in 'modules.json' but missing from directory:
    INFO Installing 'cat/cat'
    INFO Use the following statement to include this module:

        include { CAT_CAT } from '../modules/nf-core/cat/cat/main'
    ```

Il comando automaticamente:

- Scarica i file del modulo in `modules/nf-core/cat/cat/`
- Aggiorna `modules.json` per tracciare il modulo installato
- Vi fornisce l'istruzione `include` corretta da usare nel vostro flusso di lavoro

!!! tip

    Assicuratevi sempre che la vostra directory di lavoro corrente sia la radice del vostro progetto pipeline prima di eseguire il comando di installazione del modulo.

Verifichiamo che il modulo sia stato installato correttamente:

```bash
tree -L 4 modules
```

??? abstract "Directory contents"

    ```console
    modules
    ├── local
    │   ├── collectGreetings.nf
    │   ├── convertToUpper.nf
    │   ├── cowpy.nf
    │   └── sayHello.nf
    └── nf-core
        └── cat
            └── cat
                ├── environment.yml
                ├── main.nf
                ├── meta.yml
                └── tests

    5 directories, 7 files
    ```

Potete anche verificare l'installazione chiedendo all'utility nf-core di elencare i moduli installati localmente:

```bash
nf-core modules list local
```

??? success "Output del comando"

    ```console
    INFO     Repository type: pipeline
    INFO     Modules installed in '.':

    ┏━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━┓
    ┃ Module Name ┃ Repository      ┃ Version SHA ┃ Message                                ┃ Date       ┃
    ┡━━━━━━━━━━━━━╇━━━━━━━━━━━━━━━━━╇━━━━━━━━━━━━━╇━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╇━━━━━━━━━━━━┩
    │ cat/cat     │ nf-core/modules │ 41dfa3f     │ update meta.yml of all modules (#8747) │ 2025-07-07 │
    └─────────────┴─────────────────┴─────────────┴────────────────────────────────────────┴────────────┘
    ```

Questo conferma che il modulo `cat/cat` fa ora parte del codice sorgente del vostro progetto.

Tuttavia, per utilizzare effettivamente il nuovo modulo, dobbiamo importarlo nella nostra pipeline.

### 1.5. Aggiornare le importazioni dei moduli

Sostituiamo l'istruzione `include` per il modulo `collectGreetings` con quella per `CAT_CAT` nella sezione delle importazioni del flusso di lavoro `workflows/hello.nf`.

Come promemoria, lo strumento di installazione del modulo ci ha fornito l'istruzione esatta da usare:

```groovy title="Istruzione di import prodotta dal comando install"
include { CAT_CAT } from '../modules/nf-core/cat/cat/main'`
```

Notate che la convenzione nf-core è di usare lettere maiuscole per i nomi dei moduli quando li importiamo.

Aprite [core-hello/workflows/hello.nf](core-hello/workflows/hello.nf) e fate la seguente sostituzione:

=== "Dopo"

    ```groovy title="core-hello/workflows/hello.nf" linenums="1" hl_lines="10"
    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    include { paramsSummaryMap       } from 'plugin/nf-schema'
    include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
    include { sayHello               } from '../modules/local/sayHello.nf'
    include { convertToUpper         } from '../modules/local/convertToUpper.nf'
    include { CAT_CAT                } from '../modules/nf-core/cat/cat/main'
    include { cowpy                  } from '../modules/local/cowpy.nf'
    ```

=== "Prima"

    ```groovy title="core-hello/workflows/hello.nf" linenums="1" hl_lines="10"
    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    include { paramsSummaryMap       } from 'plugin/nf-schema'
    include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
    include { sayHello               } from '../modules/local/sayHello.nf'
    include { convertToUpper         } from '../modules/local/convertToUpper.nf'
    include { collectGreetings       } from '../modules/local/collectGreetings.nf'
    include { cowpy                  } from '../modules/local/cowpy.nf'
    ```

Notate come il percorso per il modulo nf-core differisce dai moduli locali:

- **Modulo nf-core**: `'../modules/nf-core/cat/cat/main'` (fa riferimento a `main.nf`)
- **Modulo locale**: `'../modules/local/collectGreetings.nf'` (riferimento a file singolo)

Il modulo è ora disponibile per il flusso di lavoro, quindi tutto ciò che dobbiamo fare è sostituire la chiamata a `collectGreetings` per usare `CAT_CAT`. Giusto?

Non così in fretta.

A questo punto, potreste essere tentati di tuffarvi e iniziare a modificare il codice, ma vale la pena prendersi un momento per esaminare attentamente cosa si aspetta il nuovo modulo e cosa produce.

Affronteremo questo aspetto come una sezione separata perché coinvolge un nuovo meccanismo che non abbiamo ancora trattato: le mappe di metadati.

!!! note

    Potete opzionalmente eliminare il file `collectGreetings.nf`:

    ```bash
    rm modules/local/collectGreetings.nf
    ```

    Tuttavia, potreste volerlo conservare come riferimento per comprendere le differenze tra moduli locali e moduli nf-core.

### Takeaway

Sapete come trovare un modulo nf-core e renderlo disponibile al vostro progetto.

### Cosa c'è dopo?

Valutare cosa richiede un nuovo modulo e identificare eventuali modifiche importanti necessarie per integrarlo in una pipeline.

---

## 2. Valutare i requisiti del nuovo modulo

Nello specifico, dobbiamo esaminare l'**interfaccia** del modulo, cioè le sue definizioni di input e output, e confrontarla con l'interfaccia del modulo che stiamo cercando di sostituire.
Questo ci permetterà di determinare se possiamo semplicemente trattare il nuovo modulo come una sostituzione diretta o se dovremo adattare parte del collegamento.

Idealmente questo è qualcosa che dovreste fare _prima_ ancora di installare il modulo, ma meglio tardi che mai.
(Per inciso, esiste un comando `uninstall` per eliminare i moduli che decidete di non volere più.)

!!! note

    Il processo CAT_CAT include una gestione piuttosto intelligente di diversi tipi di compressione, estensioni di file e così via che non sono strettamente rilevanti per quello che stiamo cercando di mostrarvi qui, quindi ignoreremo la maggior parte di esso e ci concentreremo solo sulle parti importanti.

### 2.1. Confrontare le interfacce dei due moduli

Come promemoria, ecco come appare l'interfaccia del nostro modulo `collectGreetings`:

```groovy title="modules/local/collectGreetings.nf (estratto)" linenums="1" hl_lines="6-7 10"
process collectGreetings {

    publishDir 'results', mode: 'copy'

    input:
        path input_files
        val batch_name

    output:
        path "COLLECTED-${batch_name}-output.txt" , emit: outfile
```

Il modulo `collectGreetings` accetta due input:

- `input_files` contiene uno o più file di input da elaborare;
- `batch_name` è un valore che usiamo per assegnare un nome specifico dell'esecuzione al file di output, che è una forma di metadati.

Al completamento, `collectGreetings` produce un singolo percorso di file, emesso con il tag `outfile`.

In confronto, l'interfaccia del modulo `cat/cat` è più complessa:

```groovy title="modules/nf-core/cat/cat/main.nf (estratto)" linenums="1" hl_lines="11 14"
process CAT_CAT {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pigz:2.3.4' :
        'biocontainers/pigz:2.3.4' }"

    input:
    tuple val(meta), path(files_in)

    output:
    tuple val(meta), path("${prefix}"), emit: file_out
    path "versions.yml"               , emit: versions
```

Il modulo CAT_CAT accetta un singolo input, ma quell'input è una tupla contenente due cose:

- `meta` è una struttura contenente metadati, chiamata metamap;
- `files_in` contiene uno o più file di input da elaborare, equivalente a `input_files` di `collectGreetings`.

Al completamento, CAT_CAT fornisce i suoi output in due parti:

- Un'altra tupla contenente la metamap e il file di output concatenato, emessa con il tag `file_out`;
- Un file `versions.yml` che cattura informazioni sulla versione del software utilizzata, emesso con il tag `versions`.

Notate anche che per impostazione predefinita, il file di output sarà nominato in base a un identificatore che fa parte dei metadati (codice non mostrato qui).

Questo può sembrare molto da tenere a mente guardando solo il codice, quindi ecco un diagramma per aiutarvi a visualizzare come tutto si incastra.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nf-core/img/module_comparison.svg"
</figure>

Potete vedere che i due moduli hanno requisiti di input simili in termini di contenuto (un insieme di file di input più alcuni metadati) ma aspettative molto diverse su come quel contenuto è confezionato.
Ignorando il file delle versioni per ora, anche il loro output principale è equivalente (un file concatenato), tranne che CAT_CAT emette anche la metamap insieme al file di output.

Le differenze di confezionamento saranno abbastanza facili da gestire, come vedrete tra poco.
Tuttavia, per comprendere la parte della metamap, dobbiamo introdurvi a un contesto aggiuntivo.

### 2.2. Comprendere le metamap

Vi abbiamo appena detto che il modulo CAT_CAT si aspetta una mappa di metadati come parte della sua tupla di input.
Prendiamoci qualche minuto per dare un'occhiata più da vicino a cosa sia.

La **mappa di metadati**, spesso chiamata **metamap** in breve, è una mappa in stile Groovy contenente informazioni su unità di dati.
Nel contesto delle pipeline Nextflow, le unità di dati possono essere qualsiasi cosa vogliate: campioni individuali, batch di campioni o interi dataset.

Per convenzione, una metamap nf-core è chiamata `meta` e contiene il campo obbligatorio `id`, che viene usato per nominare gli output e tracciare le unità di dati.

Ad esempio, una tipica mappa di metadati potrebbe apparire così:

```groovy title="Esempio di metamap a livello di campione"
[id: 'sample1', single_end: false, strandedness: 'forward']
```

O in un caso in cui i metadati sono allegati a livello di batch:

```groovy title="Esempio di metamap a livello di batch"
[id: 'batch1', date: '25.10.01']
```

Ora mettiamo questo nel contesto del processo `CAT_CAT`, che si aspetta che i file di input siano confezionati in una tupla con una metamap, e produce la metamap come parte della tupla di output.

```groovy title="modules/nf-core/cat/cat/main.nf (estratto)" linenums="1" hl_lines="2 5"
input:
tuple val(meta), path(files_in)

output:
tuple val(meta), path("${prefix}"), emit: file_out
```

Di conseguenza, ogni unità di dati viaggia attraverso la pipeline con i metadati rilevanti allegati.
I processi successivi possono quindi accedere facilmente anche a quei metadati.

Ricordate quando vi abbiamo detto che il file prodotto da `CAT_CAT` sarà nominato in base a un identificatore che fa parte dei metadati?
Questo è il codice rilevante:

```groovy title="modules/nf-core/cat/cat/main.nf (estratto)" linenums="35"
prefix   = task.ext.prefix ?: "${meta.id}${getFileSuffix(file_list[0])}"
```

Questo si traduce approssimativamente come segue: se viene fornito un `prefix` tramite il sistema di parametri esterni dell'attività (`task.ext`), usatelo per nominare il file di output; altrimenti createne uno usando `${meta.id}`, che corrisponde al campo `id` nella metamap.

Potete immaginare il canale di input che entra in questo modulo con contenuti come questo:

```groovy title="Esempio di contenuti del canale di input"
ch_input = [[[id: 'batch1', date: '25.10.01'], ['file1A.txt', 'file1B.txt']],
            [[id: 'batch2', date: '25.10.26'], ['file2A.txt', 'file2B.txt']],
            [[id: 'batch3', date: '25.11.14'], ['file3A.txt', 'file3B.txt']]]
```

Quindi i contenuti del canale di output che escono così:

```groovy title="Esempio di contenuti del canale di output"
ch_input = [[[id: 'batch1', date: '25.10.01'], 'batch1.txt'],
            [[id: 'batch2', date: '25.10.26'], 'batch2.txt'],
            [[id: 'batch3', date: '25.11.14'], 'batch3.txt']]
```

Come accennato in precedenza, la configurazione di input `tuple val(meta), path(files_in)` è un pattern standard usato in tutti i moduli nf-core.

Speriamo che possiate iniziare a vedere quanto questo possa essere utile.
Non solo vi permette di nominare gli output in base ai metadati, ma potete anche fare cose come usarli per applicare valori di parametri diversi, e in combinazione con operatori specifici, potete persino raggruppare, ordinare o filtrare i dati mentre fluiscono attraverso la pipeline.

!!! note "Ulteriori informazioni sui metadati"

    Per un'introduzione completa al lavoro con i metadati nei flussi di lavoro Nextflow, incluso come leggere i metadati dai samplesheet e usarli per personalizzare l'elaborazione, consultate la side quest [Metadata in workflows](../side_quests/metadata).

### 2.3. Riepilogare le modifiche da apportare

In base a quanto abbiamo esaminato, queste sono le principali modifiche che dobbiamo apportare alla nostra pipeline per utilizzare il modulo `cat/cat`:

- Creare una metamap contenente il nome del batch;
- Confezionare la metamap in una tupla con l'insieme di file di input da concatenare (provenienti da `convertToUpper`);
- Cambiare la chiamata da `collectGreetings()` a `CAT_CAT`;
- Estrarre il file di output dalla tupla prodotta dal processo `CAT_CAT` prima di passarlo a `cowpy`.

Questo dovrebbe bastare! Ora che abbiamo un piano, siamo pronti a tuffarci.

### Takeaway

Sapete come valutare l'interfaccia di input e output di un nuovo modulo per identificare i suoi requisiti, e avete imparato come le metamap vengono usate dalle pipeline nf-core per mantenere i metadati strettamente associati ai dati mentre fluiscono attraverso una pipeline.

### Cosa c'è dopo?

Integrare il nuovo modulo in un flusso di lavoro.

---

## 3. Integrare CAT_CAT nel flusso di lavoro `hello.nf`

Ora che sapete tutto sulle metamap (o abbastanza per gli scopi di questo corso, comunque), è il momento di implementare effettivamente le modifiche che abbiamo delineato sopra.

Per chiarezza, suddivideremo questo e copriremo ogni passaggio separatamente.

!!! note

    Tutte le modifiche mostrate di seguito vengono apportate alla logica del flusso di lavoro nel blocco `main` nel file del flusso di lavoro `core-hello/workflows/hello.nf`.

### 3.1. Creare una mappa di metadati

Per prima cosa, dobbiamo creare una mappa di metadati per `CAT_CAT`, tenendo presente che i moduli nf-core richiedono che la metamap abbia almeno un campo `id`.

Poiché non abbiamo bisogno di altri metadati, possiamo mantenerla semplice e usare qualcosa del genere:

```groovy title="Esempio di sintassi"
def cat_meta = [id: 'test']
```

Tranne che non vogliamo codificare il valore di `id`; vogliamo usare il valore del parametro `params.batch`.
Quindi il codice diventa:

```groovy title="Esempio di sintassi"
def cat_meta = [id: params.batch]
```

Sì, è letteralmente così semplice creare una metamap di base.

Aggiungiamo queste righe dopo la chiamata a `convertToUpper`, rimuovendo la chiamata a `collectGreetings`:

=== "Dopo"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26" hl_lines="7-8"
        // emette un saluto
        sayHello(ch_samplesheet)

        // converte il saluto in maiuscolo
        convertToUpper(sayHello.out)

        // crea una mappa di metadati con il nome del batch come ID
        def cat_meta = [ id: params.batch ]

        // genera arte ASCII dei saluti con cowpy
        cowpy(collectGreetings.out.outfile, params.character)
    ```

=== "Prima"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26" hl_lines="7-8"
        // emette un saluto
        sayHello(ch_samplesheet)

        // converte il saluto in maiuscolo
        convertToUpper(sayHello.out)

        // raccoglie tutti i saluti in un file
        collectGreetings(convertToUpper.out.collect(), params.batch)

        // genera arte ASCII dei saluti con cowpy
        cowpy(collectGreetings.out.outfile, params.character)
    ```

Questo crea una semplice mappa di metadati dove l'`id` è impostato sul nome del nostro batch (che sarà `test` quando si usa il profilo test).

### 3.2. Creare un canale con tuple di metadati

Successivamente, trasformiamo il canale di file in un canale di tuple contenenti metadati e file:

=== "Dopo"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26" hl_lines="10-11"
        // emette un saluto
        sayHello(ch_samplesheet)

        // converte il saluto in maiuscolo
        convertToUpper(sayHello.out)

        // crea una mappa di metadati con il nome del batch come ID
        def cat_meta = [ id: params.batch ]

        // crea un canale con metadati e file in formato tupla
        ch_for_cat = convertToUpper.out.collect().map { files -> tuple(cat_meta, files) }

        // genera arte ASCII dei saluti con cowpy
        cowpy(collectGreetings.out.outfile, params.character)
    ```

=== "Prima"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26"
        // emette un saluto
        sayHello(ch_samplesheet)

        // converte il saluto in maiuscolo
        convertToUpper(sayHello.out)

        // crea una mappa di metadati con il nome del batch come ID
        def cat_meta = [ id: params.batch ]

        // genera arte ASCII dei saluti con cowpy
        cowpy(collectGreetings.out.outfile, params.character)
    ```

La riga che abbiamo aggiunto realizza due cose:

- `.collect()` raccoglie tutti i file dall'output di `convertToUpper` in un'unica lista
- `.map { files -> tuple(cat_meta, files) }` crea una tupla di `[metadati, file]` nel formato che `CAT_CAT` si aspetta

Questo è tutto ciò che dobbiamo fare per configurare la tupla di input per `CAT_CAT`.

### 3.3. Chiamare il modulo CAT_CAT

Ora chiamiamo `CAT_CAT` sul canale appena creato:

=== "Dopo"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26" hl_lines="13-14"
        // emette un saluto
        sayHello(ch_samplesheet)

        // converte il saluto in maiuscolo
        convertToUpper(sayHello.out)

        // crea una mappa di metadati con il nome del batch come ID
        def cat_meta = [ id: params.batch ]

        // crea un canale con metadati e file in formato tupla
        ch_for_cat = convertToUpper.out.collect().map { files -> tuple(cat_meta, files) }

        // concatena i file usando il modulo nf-core cat/cat
        CAT_CAT(ch_for_cat)

        // genera arte ASCII dei saluti con cowpy
        cowpy(collectGreetings.out.outfile, params.character)
    ```

=== "Prima"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26"
        // emette un saluto
        sayHello(ch_samplesheet)

        // converte il saluto in maiuscolo
        convertToUpper(sayHello.out)

        // crea una mappa di metadati con il nome del batch come ID
        def cat_meta = [ id: params.batch ]

        // crea un canale con metadati e file in formato tupla
        ch_for_cat = convertToUpper.out.collect().map { files -> tuple(cat_meta, files) }

        // genera arte ASCII dei saluti con cowpy
        cowpy(collectGreetings.out.outfile, params.character)
    ```

Questo completa la parte più complicata di questa sostituzione, ma non abbiamo ancora finito: dobbiamo ancora aggiornare come passiamo l'output concatenato al processo `cowpy`.

### 3.4. Estrarre il file di output dalla tupla per `cowpy`

In precedenza, il processo `collectGreetings` produceva semplicemente un file che potevamo passare direttamente a `cowpy`.
Tuttavia, il processo `CAT_CAT` produce una tupla che include la metamap oltre al file di output.

Poiché `cowpy` non accetta ancora tuple di metadati (lo sistemeremo nella prossima parte del corso), dobbiamo estrarre il file di output dalla tupla prodotta da `CAT_CAT` prima di passarlo a `cowpy`:

=== "Dopo"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26" hl_lines="16-17 20"
        // emette un saluto
        sayHello(ch_samplesheet)

        // converte il saluto in maiuscolo
        convertToUpper(sayHello.out)

        // crea una mappa di metadati con il nome del batch come ID
        def cat_meta = [ id: params.batch ]

        // crea un canale con metadati e file in formato tupla
        ch_for_cat = convertToUpper.out.collect().map { files -> tuple(cat_meta, files) }

        // concatena i saluti
        CAT_CAT(ch_for_cat)

        // estrae il file dalla tupla poiché cowpy non usa ancora i metadati
        ch_for_cowpy = CAT_CAT.out.file_out.map{ meta, file -> file }

        // genera arte ASCII dei saluti con cowpy
        cowpy(ch_for_cowpy, params.character)
    ```

=== "Prima"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26" hl_lines="17"
        // emette un saluto
        sayHello(ch_samplesheet)

        // converte il saluto in maiuscolo
        convertToUpper(sayHello.out)

        // crea una mappa di metadati con il nome del batch come ID
        def cat_meta = [ id: params.batch ]

        // crea un canale con metadati e file in formato tupla
        ch_for_cat = convertToUpper.out.collect().map { files -> tuple(cat_meta, files) }

        // concatena i saluti
        CAT_CAT(ch_for_cat)

        // genera arte ASCII dei saluti con cowpy
        cowpy(collectGreetings.out.outfile, params.character)
    ```

L'operazione `.map{ meta, file -> file }` estrae il file dalla tupla `[metadati, file]` prodotta da `CAT_CAT` in un nuovo canale, `ch_for_cowpy`.

Quindi è solo questione di passare `ch_for_cowpy` a `cowpy` invece di `collectGreetings.out.outfile` in quell'ultima riga.

!!! note

    Nella prossima parte del corso, aggiorneremo `cowpy` per lavorare direttamente con le tuple di metadati, quindi questo passaggio di estrazione non sarà più necessario.

### 3.5. Testare il flusso di lavoro

Testiamo che il flusso di lavoro funzioni con il modulo `cat/cat` appena integrato:

```bash
nextflow run . --outdir core-hello-results -profile test,docker --validate_params false
```

Questo dovrebbe essere eseguito ragionevolmente velocemente.

??? success "Output del comando"

    ```console
    N E X T F L O W ~ version 25.04.3

        Launching `./main.nf` [evil_pike] DSL2 - revision: b9e9b3b8de

        Input/output options
          input                     : /workspaces/training/hello-nf-core/core-hello/assets/greetings.csv
          outdir                    : core-hello-results

        Institutional config options
          config_profile_name       : Test profile
          config_profile_description: Minimal test dataset to check pipeline function

        Generic options
          validate_params           : false
          trace_report_suffix       : 2025-10-30_18-50-58

        Core Nextflow options
          runName                   : evil_pike
          containerEngine           : docker
          launchDir                 : /workspaces/training/hello-nf-core/core-hello
          workDir                   : /workspaces/training/hello-nf-core/core-hello/work
          projectDir                : /workspaces/training/hello-nf-core/core-hello
          userName                  : root
          profile                   : test,docker
          configFiles               : /workspaces/training/hello-nf-core/core-hello/nextflow.config

        !! Only displaying parameters that differ from the pipeline defaults !!
        ------------------------------------------------------
        executor >  local (8)
        [b3/f005fd] CORE_HELLO:HELLO:sayHello (3)       [100%] 3 of 3 ✔
        [08/f923d0] CORE_HELLO:HELLO:convertToUpper (3) [100%] 3 of 3 ✔
        [34/3729a9] CORE_HELLO:HELLO:CAT_CAT (test)     [100%] 1 of 1 ✔
        [24/df918a] CORE_HELLO:HELLO:cowpy              [100%] 1 of 1 ✔
        -[core/hello] Pipeline completed successfully-
    ```

Notate che `CAT_CAT` ora appare nell'elenco di esecuzione dei processi invece di `collectGreetings`.

E questo è tutto! Ora stiamo usando un modulo robusto curato dalla comunità invece di codice personalizzato di livello prototipo per quel passaggio nella pipeline.

### Takeaway

Ora sapete come:

- Trovare e installare moduli nf-core
- Valutare i requisiti di un modulo nf-core
- Creare una semplice mappa di metadati da usare con un modulo nf-core
- Integrare un modulo nf-core nel vostro flusso di lavoro

### Cosa c'è dopo?

Imparare ad adattare i vostri moduli locali per seguire le convenzioni nf-core.
Vi mostreremo anche come creare nuovi moduli nf-core da un template usando gli strumenti nf-core.
