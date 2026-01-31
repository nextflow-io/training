# Parte 3: Usare un modulo nf-core

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduzione assistita da IA - [scopri di più e suggerisci miglioramenti](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

In questa terza parte del corso di formazione Hello nf-core, mostreremo come trovare, installare e utilizzare un modulo nf-core esistente nella propria pipeline.

Uno dei grandi vantaggi di lavorare con nf-core è la possibilità di sfruttare moduli pre-costruiti e testati dal repository [nf-core/modules](https://github.com/nf-core/modules).
Invece di scrivere ogni processo da zero, è possibile installare e utilizzare moduli mantenuti dalla comunità che seguono le best practice.

Per dimostrare come funziona, sostituiremo il modulo personalizzato `collectGreetings` con il modulo `cat/cat` da nf-core/modules nella pipeline `core-hello`.

??? info "Come iniziare da questa sezione"

    Questa sezione del corso presuppone che sia stata completata la [Parte 2: Riscrivere Hello per nf-core](./02_rewrite_hello.md) e si disponga di una pipeline `core-hello` funzionante.

    Se non è stata completata la Parte 2 o si desidera iniziare da zero per questa parte, è possibile utilizzare la soluzione `core-hello-part2` come punto di partenza.
    Eseguire questo comando dalla directory `hello-nf-core/`:

    ```bash
    cp -r solutions/core-hello-part2 core-hello
    cd core-hello
    ```

    Questo fornisce una pipeline nf-core completamente funzionale pronta per l'aggiunta di moduli.
    È possibile verificare che funzioni correttamente eseguendo il seguente comando:

    ```bash
    nextflow run . --outdir core-hello-results -profile test,docker --validate_params false
    ```

---

## 1. Trovare e installare un modulo nf-core adatto

Per prima cosa, impariamo come trovare un modulo nf-core esistente e installarlo nella nostra pipeline.

L'obiettivo sarà sostituire il processo `collectGreetings`, che utilizza il comando Unix `cat` per concatenare più file di saluto in uno solo.
La concatenazione di file è un'operazione molto comune, quindi è ragionevole pensare che potrebbe già esistere un modulo in nf-core progettato per questo scopo.

Iniziamo.

### 1.1. Esplorare i moduli disponibili sul sito web nf-core

Il progetto nf-core mantiene un catalogo centralizzato di moduli su [https://nf-co.re/modules](https://nf-co.re/modules).

Navigare alla pagina dei moduli nel browser web e utilizzare la barra di ricerca per cercare 'concatenate'.

![risultati ricerca moduli](./img/module-search-results.png)

Come potete vedere, ci sono diversi risultati, molti dei quali moduli progettati per concatenare tipi di file molto specifici.
Tra questi, dovrebbe essercene uno chiamato `cat_cat` che è generico.

!!! note "Convenzione di denominazione dei moduli"

    Il trattino basso (`_`) viene utilizzato come sostituto del carattere barra (`/`) nei nomi dei moduli.

    I moduli nf-core seguono la convenzione di denominazione `software/comando` quando uno strumento fornisce più comandi, come `samtools/view` (pacchetto samtools, comando view) o `gatk/haplotypecaller` (pacchetto GATK, comando HaplotypeCaller).
    Per gli strumenti che forniscono solo un comando principale, i moduli utilizzano un singolo livello come `fastqc` o `multiqc`.

Fare clic sulla casella del modulo `cat_cat` per visualizzare la documentazione del modulo.

La pagina del modulo mostra:

- Una breve descrizione: "Un modulo per la concatenazione di file compressi o non compressi con gzip"
- Comando di installazione: `nf-core modules install cat/cat`
- Struttura del canale di input e output
- Parametri disponibili

### 1.2. Elencare i moduli disponibili dalla riga di comando

In alternativa, è possibile cercare moduli direttamente dalla riga di comando utilizzando gli strumenti nf-core.

```bash
nf-core modules list remote
```

Questo visualizzerà un elenco di tutti i moduli disponibili nel repository nf-core/modules, anche se è un po' meno conveniente se non si conosce già il nome del modulo che si sta cercando.
Tuttavia, se lo si conosce, è possibile inviare l'elenco a `grep` per trovare moduli specifici:

```bash
nf-core modules list remote | grep 'cat/cat'
```

??? success "Output del comando"

    ```console
    │ cat/cat
    ```

Tenere presente che l'approccio con `grep` estrarrà solo risultati con il termine di ricerca nel loro nome, il che non funzionerebbe per `cat_cat`.

### 1.3. Ottenere informazioni dettagliate sul modulo

Per visualizzare informazioni dettagliate su un modulo specifico dalla riga di comando, utilizzare il comando `info`:

```bash
nf-core modules info cat/cat
```

Questo visualizza la documentazione sul modulo, inclusi input, output e informazioni di base sull'utilizzo.

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
      meta  (map)     │Groovy Map contenente informazioni sul    │
                      │campione es. [ id:'test', single_end:false│
                      │]                                         │
    ╶─────────────────┼──────────────────────────────────────────┼───────╴
      files_in  (file)│Elenco di file compressi / non compressi │      *
                      ╵                                          ╵
                          ╷                                 ╷
    📥 Outputs           │Description                      │     Pattern
    ╺━━━━━━━━━━━━━━━━━━━━━┿━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┿━━━━━━━━━━━━╸
    file_out             │                                 │
    ╶─────────────────────┼─────────────────────────────────┼────────────╴
      meta  (map)         │Groovy Map contenente            │
                          │informazioni sul campione        │
    ╶─────────────────────┼─────────────────────────────────┼────────────╴
      ${prefix}  (file)   │File concatenato. Sarà compresso │ ${file_out}
                          │con gzip se file_out termina con │
                          │".gz"                            │
    ╶─────────────────────┼─────────────────────────────────┼────────────╴
    versions             │                                 │
    ╶─────────────────────┼─────────────────────────────────┼────────────╴
      versions.yml  (file)│File contenente le versioni del  │versions.yml
                          │software                         │
                          ╵                                 ╵

    💻  Comando di installazione: nf-core modules install cat/cat

    ```

Queste sono esattamente le stesse informazioni che si possono trovare sul sito web.

### 1.4. Installare il modulo cat/cat

Ora che abbiamo trovato il modulo che vogliamo, dobbiamo aggiungerlo al codice sorgente della nostra pipeline.

La buona notizia è che il progetto nf-core include degli strumenti per rendere questa parte facile.
Specificamente, il comando `nf-core modules install` permette di automatizzare il recupero del codice e renderlo disponibile al proprio progetto in un singolo passaggio.

Navigare nella directory della pipeline ed eseguire il comando di installazione:

```bash
cd core-hello
nf-core modules install cat/cat
```

Lo strumento potrebbe prima richiedere di specificare un tipo di repository.
(In caso contrario, saltare alla parte "Infine, lo strumento procederà con l'installazione del modulo.")

??? success "Output del comando"

    ```console

                                          ,--./,-.
          ___     __   __   __   ___     /,-._.--~\
    |\ | |__  __ /  ` /  \ |__) |__         }  {
    | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                          `._,._,'

    nf-core/tools version 3.4.1 - https://nf-co.re


    WARNING  'repository_type' non definito in .nf-core.yml
    ? Questo repository è una pipeline o un repository di moduli? (Usa i tasti freccia)
    » Pipeline
      Repository di moduli
    ```

Se appare, premere invio per accettare la risposta predefinita (`Pipeline`) e continuare.

Lo strumento offrirà quindi di modificare la configurazione del progetto per evitare questo prompt in futuro.

??? success "Output del comando"

    ```console
        INFO     Per evitare questo prompt in futuro, aggiungere la chiave 'repository_type' al file .nf-core.yml.
        ? Desideri che aggiunga questa configurazione ora? [y/n] (y):
    ```

Tanto vale approfittare di questo strumento conveniente!
Premere invio per accettare la risposta predefinita (sì).

Infine, lo strumento procederà con l'installazione del modulo.

??? success "Output del comando"

    ```console
    INFO Configurazione aggiunta a '.nf-core.yml'
    INFO Reinstallazione dei moduli trovati in 'modules.json' ma mancanti dalla directory:
    INFO Installazione di 'cat/cat'
    INFO Utilizzare la seguente istruzione per includere questo modulo:

        include { CAT_CAT } from '../modules/nf-core/cat/cat/main'
    ```

Il comando automaticamente:

- Scarica i file del modulo in `modules/nf-core/cat/cat/`
- Aggiorna `modules.json` per tracciare il modulo installato
- Fornisce l'istruzione `include` corretta da utilizzare nel workflow

!!! tip

    Assicurarsi sempre che la directory di lavoro corrente sia la radice del progetto della pipeline prima di eseguire il comando di installazione del modulo.

Verifichiamo che il modulo sia stato installato correttamente:

```bash
tree -L 4 modules
```

??? abstract "Contenuto della directory"

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

È anche possibile verificare l'installazione chiedendo all'utility nf-core di elencare i moduli installati localmente:

```bash
nf-core modules list local
```

??? success "Output del comando"

    ```console
    INFO     Tipo di repository: pipeline
    INFO     Moduli installati in '.':

    ┏━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━┓
    ┃ Nome Modulo ┃ Repository      ┃ Version SHA ┃ Messaggio                              ┃ Data       ┃
    ┡━━━━━━━━━━━━━╇━━━━━━━━━━━━━━━━━╇━━━━━━━━━━━━━╇━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╇━━━━━━━━━━━━┩
    │ cat/cat     │ nf-core/modules │ 41dfa3f     │ update meta.yml of all modules (#8747) │ 2025-07-07 │
    └─────────────┴─────────────────┴─────────────┴────────────────────────────────────────┴────────────┘
    ```

Questo conferma che il modulo `cat/cat` fa ora parte del codice sorgente del progetto.

Tuttavia, per utilizzare effettivamente il nuovo modulo, dobbiamo importarlo nella nostra pipeline.

### 1.5. Aggiornare le importazioni dei moduli

Sostituiamo l'istruzione `include` per il modulo `collectGreetings` con quella per `CAT_CAT` nella sezione delle importazioni del workflow `workflows/hello.nf`.

Come promemoria, lo strumento di installazione del modulo ci ha fornito l'esatta istruzione da utilizzare:

```groovy title="Istruzione di importazione prodotta dal comando di installazione"
include { CAT_CAT } from '../modules/nf-core/cat/cat/main'`
```

Si noti che la convenzione nf-core prevede l'uso di maiuscole per i nomi dei moduli quando li si importa.

Aprire [core-hello/workflows/hello.nf](core-hello/workflows/hello.nf) e apportare la seguente sostituzione:

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

Notare come il percorso per il modulo nf-core differisca dai moduli locali:

- **modulo nf-core**: `'../modules/nf-core/cat/cat/main'` (riferimento a `main.nf`)
- **modulo locale**: `'../modules/local/collectGreetings.nf'` (riferimento a singolo file)

Il modulo è ora disponibile per il workflow, quindi tutto ciò che dobbiamo fare è sostituire la chiamata a `collectGreetings` per usare `CAT_CAT`. Giusto?

Non così in fretta.

A questo punto, si potrebbe essere tentati di immergersi e iniziare a modificare il codice, ma vale la pena prendersi un momento per esaminare attentamente cosa si aspetta il nuovo modulo e cosa produce.

Affronteremo questo aspetto come sezione separata perché coinvolge un nuovo meccanismo che non abbiamo ancora trattato: le mappe di metadati.

!!! note

    È possibile eliminare facoltativamente il file `collectGreetings.nf`:

    ```bash
    rm modules/local/collectGreetings.nf
    ```

    Tuttavia, potrebbe essere utile mantenerlo come riferimento per comprendere le differenze tra moduli locali e nf-core.

### Riepilogo

Ora si sa come trovare un modulo nf-core e renderlo disponibile al proprio progetto.

### Prossimi passi

Valutare cosa richiede un nuovo modulo e identificare eventuali modifiche importanti necessarie per integrarlo in una pipeline.

---

## 2. Valutare i requisiti del nuovo modulo

Nello specifico, dobbiamo esaminare l'**interfaccia** del modulo, cioè le sue definizioni di input e output, e confrontarla con l'interfaccia del modulo che stiamo cercando di sostituire.
Questo ci permetterà di determinare se possiamo semplicemente trattare il nuovo modulo come una sostituzione diretta o se dovremo adattare parte del collegamento.

Idealmente questo è qualcosa che si dovrebbe fare _prima_ di installare il modulo, ma meglio tardi che mai.
(Tanto per la cronaca, esiste un comando `uninstall` per eliminare i moduli che si decide di non volere più.)

!!! note

    Il processo CAT_CAT include una gestione piuttosto intelligente di diversi tipi di compressione, estensioni di file e così via che non sono strettamente rilevanti per ciò che stiamo cercando di mostrarvi qui, quindi ignoreremo la maggior parte di essi e ci concentreremo solo sulle parti che sono importanti.

### 2.1. Confrontare le interfacce dei due moduli

Come promemoria, questa è l'interfaccia del nostro modulo `collectGreetings`:

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
- `batch_name` è un valore che utilizziamo per assegnare un nome specifico all'esecuzione al file di output, che è una forma di metadati.

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

Si noti anche che per impostazione predefinita, il file di output sarà denominato in base a un identificatore che fa parte dei metadati (codice non mostrato qui).

Questo potrebbe sembrare molto da tenere a mente guardando solo il codice, quindi ecco un diagramma per aiutare a visualizzare come tutto si integra.

<figure class="excalidraw">
--8<-- "docs/hello_nf-core/img/module_comparison.svg"
</figure>

Potete vedere che i due moduli hanno requisiti di input simili in termini di contenuto (un insieme di file di input più alcuni metadati) ma aspettative molto diverse per come quel contenuto è confezionato.
Ignorando per ora il file delle versioni, anche il loro output principale è equivalente (un file concatenato), tranne che CAT_CAT emette anche la metamap insieme al file di output.

Le differenze di confezionamento saranno abbastanza facili da gestire, come vedrete tra poco.
Tuttavia, per comprendere la parte della metamap, dobbiamo introdurre qualche contesto aggiuntivo.

### 2.2. Comprendere le metamap

Abbiamo appena detto che il modulo CAT_CAT si aspetta una mappa di metadati come parte della sua tupla di input.
Prendiamoci qualche minuto per dare un'occhiata più da vicino a cosa sia.

La **mappa di metadati**, spesso chiamata **metamap** in breve, è una mappa in stile Groovy contenente informazioni sulle unità di dati.
Nel contesto delle pipeline Nextflow, le unità di dati possono essere qualsiasi cosa si desideri: campioni individuali, batch di campioni o interi dataset.

Per convenzione, una metamap nf-core è denominata `meta` e contiene il campo obbligatorio `id`, che viene utilizzato per denominare gli output e tracciare le unità di dati.

Ad esempio, una tipica mappa di metadati potrebbe apparire così:

```groovy title="Esempio di metamap a livello di campione"
[id: 'sample1', single_end: false, strandedness: 'forward']
```

O nel caso in cui i metadati siano associati a livello di batch:

```groovy title="Esempio di metamap a livello di batch"
[id: 'batch1', date: '25.10.01']
```

Ora mettiamo questo nel contesto del processo `CAT_CAT`, che si aspetta che i file di input siano confezionati in una tupla con una metamap, ed emette anche la metamap come parte della tupla di output.

```groovy title="modules/nf-core/cat/cat/main.nf (estratto)" linenums="1" hl_lines="2 5"
input:
tuple val(meta), path(files_in)

output:
tuple val(meta), path("${prefix}"), emit: file_out
```

Di conseguenza, ogni unità di dati viaggia attraverso la pipeline con i metadati rilevanti allegati.
I processi successivi possono quindi accedere facilmente anche a quei metadati.

Ricordate quando vi abbiamo detto che il file prodotto da `CAT_CAT` sarà denominato in base a un identificatore che fa parte dei metadati?
Questo è il codice rilevante:

```groovy title="modules/nf-core/cat/cat/main.nf (estratto)" linenums="35"
prefix   = task.ext.prefix ?: "${meta.id}${getFileSuffix(file_list[0])}"
```

Questo si traduce approssimativamente come segue: se un `prefix` viene fornito tramite il sistema di parametri esterni delle attività (`task.ext`), usarlo per denominare il file di output; altrimenti crearne uno usando `${meta.id}`, che corrisponde al campo `id` nella metamap.

Potete immaginare il canale di input che arriva a questo modulo con contenuti come questo:

```groovy title="Esempio di contenuto del canale di input"
ch_input = [[[id: 'batch1', date: '25.10.01'], ['file1A.txt', 'file1B.txt']],
            [[id: 'batch2', date: '25.10.26'], ['file2A.txt', 'file2B.txt']],
            [[id: 'batch3', date: '25.11.14'], ['file3A.txt', 'file3B.txt']]]
```

Quindi il contenuto del canale di output che esce così:

```groovy title="Esempio di contenuto del canale di output"
ch_input = [[[id: 'batch1', date: '25.10.01'], 'batch1.txt'],
            [[id: 'batch2', date: '25.10.26'], 'batch2.txt'],
            [[id: 'batch3', date: '25.11.14'], 'batch3.txt']]
```

Come accennato in precedenza, la configurazione di input `tuple val(meta), path(files_in)` è un pattern standard utilizzato in tutti i moduli nf-core.

Speriamo che possiate iniziare a vedere quanto questo possa essere utile.
Non solo permette di denominare gli output in base ai metadati, ma è possibile anche fare cose come usarli per applicare valori di parametri diversi, e in combinazione con operatori specifici, è possibile persino raggruppare, ordinare o filtrare i dati mentre fluiscono attraverso la pipeline.

!!! note "Ulteriori informazioni sui metadati"

    Per un'introduzione completa al lavoro con i metadati nei workflow Nextflow, incluso come leggere i metadati dai samplesheet e utilizzarli per personalizzare l'elaborazione, consultare la missione secondaria [Metadati nei workflow](../side_quests/metadata).

### 2.3. Riepilogare le modifiche da apportare

In base a quanto abbiamo esaminato, queste sono le modifiche principali che dobbiamo apportare alla nostra pipeline per utilizzare il modulo `cat/cat`:

- Creare una metamap contenente il nome del batch;
- Confezionare la metamap in una tupla con l'insieme di file di input da concatenare (provenienti da `convertToUpper`);
- Cambiare la chiamata da `collectGreetings()` a `CAT_CAT`;
- Estrarre il file di output dalla tupla prodotta dal processo `CAT_CAT` prima di passarlo a `cowpy`.

Questo dovrebbe essere sufficiente! Ora che abbiamo un piano, siamo pronti ad immergerci.

### Riepilogo

Si sa come valutare l'interfaccia di input e output di un nuovo modulo per identificare i suoi requisiti, e si è appreso come le metamap vengono utilizzate dalle pipeline nf-core per mantenere i metadati strettamente associati ai dati mentre fluiscono attraverso una pipeline.

### Prossimi passi

Integrare il nuovo modulo in un workflow.

---

## 3. Integrare CAT_CAT nel workflow `hello.nf`

Ora che sapete tutto sulle metamap (o abbastanza per gli scopi di questo corso, comunque), è tempo di implementare effettivamente le modifiche che abbiamo delineato sopra.

Per chiarezza, suddivideremo questo processo e tratteremo ogni passaggio separatamente.

!!! note

    Tutte le modifiche mostrate di seguito sono apportate alla logica del workflow nel blocco `main` nel file del workflow `core-hello/workflows/hello.nf`.

### 3.1. Creare una mappa di metadati

Per prima cosa, dobbiamo creare una mappa di metadati per `CAT_CAT`, tenendo presente che i moduli nf-core richiedono che la metamap abbia almeno un campo `id`.

Poiché non abbiamo bisogno di altri metadati, possiamo mantenerla semplice e usare qualcosa come questo:

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
        // emettere un saluto
        sayHello(ch_samplesheet)

        // convertire il saluto in maiuscolo
        convertToUpper(sayHello.out)

        // creare mappa di metadati con il nome del batch come ID
        def cat_meta = [ id: params.batch ]

        // generare arte ASCII dei saluti con cowpy
        cowpy(collectGreetings.out.outfile, params.character)
    ```

=== "Prima"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26" hl_lines="7-8"
        // emettere un saluto
        sayHello(ch_samplesheet)

        // convertire il saluto in maiuscolo
        convertToUpper(sayHello.out)

        // raccogliere tutti i saluti in un file
        collectGreetings(convertToUpper.out.collect(), params.batch)

        // generare arte ASCII dei saluti con cowpy
        cowpy(collectGreetings.out.outfile, params.character)
    ```

Questo crea una semplice mappa di metadati dove l'`id` è impostato sul nome del nostro batch (che sarà `test` quando si usa il profilo test).

### 3.2. Creare un canale con tuple di metadati

Successivamente, trasformare il canale di file in un canale di tuple contenenti metadati e file:

=== "Dopo"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26" hl_lines="10-11"
        // emettere un saluto
        sayHello(ch_samplesheet)

        // convertire il saluto in maiuscolo
        convertToUpper(sayHello.out)

        // creare mappa di metadati con il nome del batch come ID
        def cat_meta = [ id: params.batch ]

        // creare un canale con metadati e file in formato tupla
        ch_for_cat = convertToUpper.out.collect().map { files -> tuple(cat_meta, files) }

        // generare arte ASCII dei saluti con cowpy
        cowpy(collectGreetings.out.outfile, params.character)
    ```

=== "Prima"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26"
        // emettere un saluto
        sayHello(ch_samplesheet)

        // convertire il saluto in maiuscolo
        convertToUpper(sayHello.out)

        // creare mappa di metadati con il nome del batch come ID
        def cat_meta = [ id: params.batch ]

        // generare arte ASCII dei saluti con cowpy
        cowpy(collectGreetings.out.outfile, params.character)
    ```

La riga che abbiamo aggiunto realizza due cose:

- `.collect()` raccoglie tutti i file dall'output di `convertToUpper` in un singolo elenco
- `.map { files -> tuple(cat_meta, files) }` crea una tupla di `[metadati, file]` nel formato che `CAT_CAT` si aspetta

Questo è tutto ciò che dobbiamo fare per configurare la tupla di input per `CAT_CAT`.

### 3.3. Chiamare il modulo CAT_CAT

Ora chiamare `CAT_CAT` sul canale appena creato:

=== "Dopo"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26" hl_lines="13-14"
        // emettere un saluto
        sayHello(ch_samplesheet)

        // convertire il saluto in maiuscolo
        convertToUpper(sayHello.out)

        // creare mappa di metadati con il nome del batch come ID
        def cat_meta = [ id: params.batch ]

        // creare un canale con metadati e file in formato tupla
        ch_for_cat = convertToUpper.out.collect().map { files -> tuple(cat_meta, files) }

        // concatenare i file utilizzando il modulo nf-core cat/cat
        CAT_CAT(ch_for_cat)

        // generare arte ASCII dei saluti con cowpy
        cowpy(collectGreetings.out.outfile, params.character)
    ```

=== "Prima"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26"
        // emettere un saluto
        sayHello(ch_samplesheet)

        // convertire il saluto in maiuscolo
        convertToUpper(sayHello.out)

        // creare mappa di metadati con il nome del batch come ID
        def cat_meta = [ id: params.batch ]

        // creare un canale con metadati e file in formato tupla
        ch_for_cat = convertToUpper.out.collect().map { files -> tuple(cat_meta, files) }

        // generare arte ASCII dei saluti con cowpy
        cowpy(collectGreetings.out.outfile, params.character)
    ```

Questo completa la parte più complessa di questa sostituzione, ma non abbiamo ancora finito: dobbiamo ancora aggiornare il modo in cui passiamo l'output concatenato al processo `cowpy`.

### 3.4. Estrarre il file di output dalla tupla per `cowpy`

In precedenza, il processo `collectGreetings` produceva solo un file che potevamo passare a `cowpy` direttamente.
Tuttavia, il processo `CAT_CAT` produce una tupla che include la metamap oltre al file di output.

Poiché `cowpy` non accetta ancora tuple di metadati (risolveremo questo problema nella prossima parte del corso), dobbiamo estrarre il file di output dalla tupla prodotta da `CAT_CAT` prima di passarlo a `cowpy`:

=== "Dopo"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26" hl_lines="16-17 20"
        // emettere un saluto
        sayHello(ch_samplesheet)

        // convertire il saluto in maiuscolo
        convertToUpper(sayHello.out)

        // creare mappa di metadati con il nome del batch come ID
        def cat_meta = [ id: params.batch ]

        // creare un canale con metadati e file in formato tupla
        ch_for_cat = convertToUpper.out.collect().map { files -> tuple(cat_meta, files) }

        // concatenare i saluti
        CAT_CAT(ch_for_cat)

        // estrarre il file dalla tupla poiché cowpy non usa ancora metadati
        ch_for_cowpy = CAT_CAT.out.file_out.map{ meta, file -> file }

        // generare arte ASCII dei saluti con cowpy
        cowpy(ch_for_cowpy, params.character)
    ```

=== "Prima"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26" hl_lines="17"
        // emettere un saluto
        sayHello(ch_samplesheet)

        // convertire il saluto in maiuscolo
        convertToUpper(sayHello.out)

        // creare mappa di metadati con il nome del batch come ID
        def cat_meta = [ id: params.batch ]

        // creare un canale con metadati e file in formato tupla
        ch_for_cat = convertToUpper.out.collect().map { files -> tuple(cat_meta, files) }

        // concatenare i saluti
        CAT_CAT(ch_for_cat)

        // generare arte ASCII dei saluti con cowpy
        cowpy(collectGreetings.out.outfile, params.character)
    ```

L'operazione `.map{ meta, file -> file }` estrae il file dalla tupla `[metadati, file]` prodotta da `CAT_CAT` in un nuovo canale, `ch_for_cowpy`.

Quindi è solo questione di passare `ch_for_cowpy` a `cowpy` invece di `collectGreetings.out.outfile` in quest'ultima riga.

!!! note

    Nella prossima parte del corso, aggiorneremo `cowpy` per lavorare direttamente con tuple di metadati, quindi questo passaggio di estrazione non sarà più necessario.

### 3.5. Testare il workflow

Testiamo che il workflow funzioni con il modulo `cat/cat` appena integrato:

```bash
nextflow run . --outdir core-hello-results -profile test,docker --validate_params false
```

Questo dovrebbe eseguirsi ragionevolmente rapidamente.

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

Si noti che `CAT_CAT` appare ora nell'elenco di esecuzione dei processi invece di `collectGreetings`.

E questo è tutto! Stiamo ora utilizzando un modulo robusto curato dalla comunità invece di codice personalizzato di livello prototipale per quel passaggio nella pipeline.

### Riepilogo

Ora si sa come:

- Trovare e installare moduli nf-core
- Valutare i requisiti di un modulo nf-core
- Creare una semplice mappa di metadati da utilizzare con un modulo nf-core
- Integrare un modulo nf-core nel proprio workflow

### Prossimi passi

Imparare ad adattare i propri moduli locali per seguire le convenzioni nf-core.
Mostreremo anche come creare nuovi moduli nf-core da un template utilizzando gli strumenti nf-core.
