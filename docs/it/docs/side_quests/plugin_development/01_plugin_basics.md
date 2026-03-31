# Parte 1: Nozioni di Base sui Plugin

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduzione assistita da IA - [scopri di più e suggerisci miglioramenti](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

In questa sezione, imparerete come i plugin estendono Nextflow, poi proverete tre plugin diversi per vederli in azione.

---

## 1. Come funzionano i plugin

I plugin estendono Nextflow attraverso diversi tipi di estensione:

| Tipo di Estensione    | Cosa fa                                                    | Esempio                      |
| --------------------- | ---------------------------------------------------------- | ---------------------------- |
| Funzioni              | Aggiunge funzioni personalizzate chiamabili dai flussi di lavoro | `samplesheetToList()`        |
| Monitor del flusso di lavoro | Risponde a eventi come il completamento di un'attività | Log personalizzati, alert Slack |
| Executor              | Aggiunge backend per l'esecuzione delle attività           | AWS Batch, Kubernetes        |
| Filesystem            | Aggiunge backend di archiviazione                          | S3, Azure Blob               |

Le funzioni e i monitor del flusso di lavoro (chiamati "trace observer" nell'API di Nextflow) sono i tipi più comuni per gli autori di plugin.
Gli executor e i filesystem sono tipicamente creati dai fornitori di piattaforme.

I prossimi esercizi mostrano plugin di tipo funzione e un plugin observer, così potrete vedere entrambi i tipi in azione.

---

## 2. Usare plugin di tipo funzione

I plugin di tipo funzione aggiungono funzioni chiamabili che importate nei vostri flussi di lavoro.
Ne proverete due: nf-hello (un esempio semplice) e nf-schema (un plugin reale ampiamente utilizzato).
Entrambi gli esercizi modificano la stessa pipeline `hello.nf`, così potrete vedere come i plugin migliorano un flusso di lavoro esistente.

### 2.1. nf-hello: sostituire codice scritto a mano

Il plugin [nf-hello](https://github.com/nextflow-io/nf-hello) fornisce una funzione `randomString` che genera stringhe casuali.
La pipeline definisce già la propria versione inline di questa funzione, che sostituirete con quella del plugin.

#### 2.1.1. Vedere il punto di partenza

Esaminate la pipeline:

```bash
cat hello.nf
```

```groovy title="Output"
#!/usr/bin/env nextflow

params.input = 'greetings.csv'

/**
 * Genera una stringa alfanumerica casuale
 */
def randomString(int length) {
    def chars = ('a'..'z') + ('A'..'Z') + ('0'..'9')
    def random = new Random()
    return (1..length).collect { chars[random.nextInt(chars.size())] }.join()
}

process SAY_HELLO {
    input:
        val greeting
    output:
        stdout
    script:
    """
    echo '$greeting'
    """
}

workflow {
    greeting_ch = channel.fromPath(params.input)
        .splitCsv(header: true)
        .map { row -> "${row.greeting}_${randomString(8)}" }
    SAY_HELLO(greeting_ch)
    SAY_HELLO.out.view { result -> "Output: ${result.trim()}" }
}
```

La pipeline definisce la propria funzione `randomString` inline, poi la usa per aggiungere un ID casuale a ciascun saluto.

Eseguitela:

```bash
nextflow run hello.nf
```

```console title="Output"
Output: Hello_aBcDeFgH
Output: Bonjour_xYzWvUtS
Output: Holà_qRsPdMnK
Output: Ciao_jLhGfEcB
Output: Hallo_tNwOiAuR
```

L'ordine dell'output e le stringhe casuali saranno diversi, e se eseguite di nuovo lo script otterrete un insieme diverso di saluti casuali.

#### 2.1.2. Configurare il plugin

Sostituite la funzione inline con una proveniente dal plugin. Aggiungete questo al vostro `nextflow.config`:

```groovy title="nextflow.config"
// Configurazione per gli esercizi di sviluppo dei plugin
plugins {
    id 'nf-hello@0.5.0'
}
```

I plugin vengono dichiarati in `nextflow.config` usando il blocco `plugins {}`.
Nextflow li scarica automaticamente dal [Nextflow Plugin Registry](https://registry.nextflow.io/), un repository centrale di plugin della community e ufficiali.

#### 2.1.3. Usare la funzione del plugin

Sostituite la funzione `randomString` inline con la versione del plugin:

=== "Dopo"

    ```groovy title="hello.nf" hl_lines="3"
    #!/usr/bin/env nextflow

    include { randomString } from 'plugin/nf-hello'

    params.input = 'greetings.csv'

    process SAY_HELLO {
        input:
            val greeting
        output:
            stdout
        script:
        """
        echo '$greeting'
        """
    }

    workflow {
        greeting_ch = channel.fromPath(params.input)
            .splitCsv(header: true)
            .map { row -> "${row.greeting}_${randomString(8)}" }
        SAY_HELLO(greeting_ch)
        SAY_HELLO.out.view { result -> "Output: ${result.trim()}" }
    }
    ```

=== "Prima"

    ```groovy title="hello.nf" hl_lines="5-12"
    #!/usr/bin/env nextflow

    params.input = 'greetings.csv'

    /**
     * Genera una stringa alfanumerica casuale
     */
    def randomString(int length) {
        def chars = ('a'..'z') + ('A'..'Z') + ('0'..'9')
        def random = new Random()
        return (1..length).collect { chars[random.nextInt(chars.size())] }.join()
    }

    process SAY_HELLO {
        input:
            val greeting
        output:
            stdout
        script:
        """
        echo '$greeting'
        """
    }

    workflow {
        greeting_ch = channel.fromPath(params.input)
            .splitCsv(header: true)
            .map { row -> "${row.greeting}_${randomString(8)}" }
        SAY_HELLO(greeting_ch)
        SAY_HELLO.out.view { result -> "Output: ${result.trim()}" }
    }
    ```

L'istruzione `include` importa `randomString` da una libreria collaudata, testata e mantenuta da un gruppo più ampio di contributori che possono individuare e correggere i bug.
Invece di avere ogni pipeline con la propria copia della funzione, ogni pipeline che usa il plugin ottiene la stessa implementazione verificata.
Questo riduce il codice duplicato e il carico di manutenzione che ne deriva.
La sintassi `#!groovy include { function } from 'plugin/plugin-id'` è la stessa `include` usata per i moduli Nextflow, con il prefisso `plugin/`.
Potete vedere il [codice sorgente di `randomString`](https://github.com/nextflow-io/nf-hello/blob/e67bddebfa589c7ae51f41bf780c92068dc09e93/plugins/nf-hello/src/main/nextflow/hello/HelloExtension.groovy#L110) nel repository nf-hello su GitHub.

#### 2.1.4. Eseguirla

```bash
nextflow run hello.nf
```

```console title="Output"
Pipeline is starting! 🚀
Output: Hello_yqvtclcc
Output: Bonjour_vwwpyzcs
Output: Holà_wrghmgab
Output: Ciao_noniajuy
Output: Hallo_tvrtuxtp
Pipeline complete! 👋
```

(Le vostre stringhe casuali saranno diverse.)

L'output ha ancora suffissi casuali, ma ora `randomString` proviene dal plugin nf-hello invece che dal codice inline.
I messaggi "Pipeline is starting!" e "Pipeline complete!" sono nuovi.
Provengono dal componente observer del plugin, che esplorerete nella Parte 5.

Nextflow scarica i plugin automaticamente la prima volta che vengono utilizzati, quindi qualsiasi pipeline che dichiara `nf-hello@0.5.0` ottiene la stessa funzione `randomString` testata senza dover copiare codice tra i progetti.

Avete ora visto i tre passaggi per usare un plugin di tipo funzione: dichiararlo in `nextflow.config`, importare la funzione con `include`, e chiamarla nel vostro flusso di lavoro.
Il prossimo esercizio applica questi stessi passaggi a un plugin reale.

### 2.2. nf-schema: parsing CSV con validazione

Il plugin [nf-schema](https://github.com/nextflow-io/nf-schema) è uno dei plugin Nextflow più ampiamente utilizzati.
Fornisce `samplesheetToList`, una funzione che analizza file CSV/TSV usando uno schema JSON che definisce le colonne e i tipi attesi.

La pipeline attualmente legge `greetings.csv` usando `splitCsv` e un `map` manuale, ma nf-schema può sostituire questo approccio con un parsing validato e guidato dallo schema.
Un file di schema JSON (`greetings_schema.json`) è già fornito nella directory dell'esercizio.

??? info "Cos'è uno schema?"

    Uno schema è una descrizione formale di come appaiono i dati validi.
    Definisce cose come quali colonne sono attese, che tipo deve avere ciascun valore (stringa, numero, ecc.) e quali campi sono obbligatori.

    Pensatelo come un contratto: se i dati di input non corrispondono allo schema, lo strumento può individuare il problema subito invece di lasciare che causi errori confusi più avanti nella pipeline.

#### 2.2.1. Esaminare lo schema

```bash
cat greetings_schema.json
```

```json title="Output"
{
  "$schema": "https://json-schema.org/draft/2020-12/schema",
  "type": "array",
  "items": {
    "type": "object",
    "properties": {
      "greeting": {
        "type": "string",
        "description": "The greeting text"
      },
      "language": {
        "type": "string",
        "description": "The language of the greeting"
      }
    },
    "required": ["greeting"]
  }
}
```

Lo schema definisce due colonne (`greeting` e `language`) e contrassegna `greeting` come obbligatoria.
Se qualcuno passa un CSV privo della colonna `greeting`, nf-schema intercetta l'errore prima che la pipeline venga eseguita.

#### 2.2.2. Aggiungere nf-schema alla configurazione

Aggiornate `nextflow.config` per includere entrambi i plugin:

=== "Dopo"

    ```groovy title="nextflow.config" hl_lines="3"
    plugins {
        id 'nf-hello@0.5.0'
        id 'nf-schema@2.6.1'
    }
    ```

=== "Prima"

    ```groovy title="nextflow.config"
    plugins {
        id 'nf-hello@0.5.0'
    }
    ```

#### 2.2.3. Aggiornare hello.nf per usare samplesheetToList

Sostituite l'input con `splitCsv` usando `samplesheetToList`:

=== "Dopo"

    ```groovy title="hello.nf" hl_lines="4 20 21 22"
    #!/usr/bin/env nextflow

    include { randomString } from 'plugin/nf-hello'
    include { samplesheetToList } from 'plugin/nf-schema'

    params.input = 'greetings.csv'

    process SAY_HELLO {
        input:
            val greeting
        output:
            stdout
        script:
        """
        echo '$greeting'
        """
    }

    workflow {
        def samplesheet_list = samplesheetToList(params.input, 'greetings_schema.json')
        greeting_ch = Channel.fromList(samplesheet_list)
            .map { row -> "${row[0]}_${randomString(8)}" }
        SAY_HELLO(greeting_ch)
        SAY_HELLO.out.view { result -> "Output: ${result.trim()}" }
    }
    ```

=== "Prima"

    ```groovy title="hello.nf" hl_lines="19 20 21"
    #!/usr/bin/env nextflow

    include { randomString } from 'plugin/nf-hello'

    params.input = 'greetings.csv'

    process SAY_HELLO {
        input:
            val greeting
        output:
            stdout
        script:
        """
        echo '$greeting'
        """
    }

    workflow {
        greeting_ch = channel.fromPath(params.input)
            .splitCsv(header: true)
            .map { row -> "${row.greeting}_${randomString(8)}" }
        SAY_HELLO(greeting_ch)
        SAY_HELLO.out.view { result -> "Output: ${result.trim()}" }
    }
    ```

Il codice personalizzato di parsing con `splitCsv` e `map` viene sostituito da `samplesheetToList`, una funzione collaudata e testata che valida anche il samplesheet rispetto allo schema prima che la pipeline venga eseguita.
Questo riduce il carico di manutenzione della logica di parsing scritta a mano, migliorando al contempo l'esperienza per gli utenti della pipeline, che ricevono messaggi di errore chiari quando il loro input non corrisponde al formato atteso.
Ogni riga diventa una lista di valori nell'ordine delle colonne, quindi `row[0]` è il saluto e `row[1]` è la lingua.

#### 2.2.4. Eseguirla

```bash
nextflow run hello.nf
```

```console title="Output"
Pipeline is starting! 🚀
Output: Hello_diozjdwm
Output: Bonjour_speathmm
Output: Holà_dllxnzap
Output: Ciao_wzueddzc
Output: Hallo_hsxwrjbh
Pipeline complete! 👋
```

(Le vostre stringhe casuali saranno diverse.)

L'output è lo stesso, ma ora lo schema valida la struttura del CSV prima che la pipeline venga eseguita.
Nelle pipeline reali con samplesheet complessi e molte colonne, questo tipo di validazione previene errori che `splitCsv` + `map` manuali non rileverebbero.

#### 2.2.5. Vedere la validazione in azione

Per vedere cosa intercetta la validazione dello schema, provate a introdurre degli errori in `greetings.csv`.

Rinominate la colonna obbligatoria `greeting` in `message`:

```csv title="greetings.csv" hl_lines="1"
message,language
Hello,English
Bonjour,French
Holà,Spanish
Ciao,Italian
Hallo,German
```

Eseguite la pipeline:

```bash
nextflow run hello.nf
```

```console title="Output"
ERROR ~ Validation of samplesheet failed!

The following errors have been detected in greetings.csv:

-> Entry 1: Missing required field(s): greeting
-> Entry 2: Missing required field(s): greeting
-> Entry 3: Missing required field(s): greeting
-> Entry 4: Missing required field(s): greeting
-> Entry 5: Missing required field(s): greeting
```

La pipeline si rifiuta di eseguire perché lo schema richiede una colonna `greeting` e non riesce a trovarla.

Ora ripristinate la colonna obbligatoria ma rinominate la colonna opzionale `language` in `lang`:

```csv title="greetings.csv" hl_lines="1"
greeting,lang
Hello,English
Bonjour,French
Holà,Spanish
Ciao,Italian
Hallo,German
```

```bash
nextflow run hello.nf
```

Questa volta la pipeline viene eseguita, ma stampa un avviso:

```console title="Output (partial)"
WARN: Found the following unidentified headers in greetings.csv:
	- lang
```

Le colonne obbligatorie causano errori bloccanti; le colonne opzionali causano avvisi.
Questo è il tipo di feedback precoce che fa risparmiare tempo nel debugging nelle pipeline reali con decine di colonne.

#### 2.2.6. Configurare il comportamento della validazione

L'avviso riguardante `lang` è utile, ma potete controllarne la severità attraverso la configurazione.
I plugin possono includere i propri scope di configurazione che ne controllano il comportamento.
Il plugin nf-schema include lo scope di configurazione `validation`; modificando le impostazioni qui potete cambiare il comportamento di nf-schema.

Aggiungete un blocco `validation` a `nextflow.config` per fare in modo che le intestazioni non riconosciute causino un errore invece di un avviso:

=== "Dopo"

    ```groovy title="nextflow.config" hl_lines="6-10"
    plugins {
        id 'nf-hello@0.5.0'
        id 'nf-schema@2.6.1'
    }

    validation {
        logging {
            unrecognisedHeaders = "error"
        }
    }
    ```

=== "Prima"

    ```groovy title="nextflow.config"
    plugins {
        id 'nf-hello@0.5.0'
        id 'nf-schema@2.6.1'
    }
    ```

Eseguite di nuovo la pipeline con la colonna `lang` ancora presente:

```bash
nextflow run hello.nf
```

```console title="Output (partial)"
Found the following unidentified headers in greetings.csv:
	- lang
 -- Check script 'hello.nf' at line: 20 or see '.nextflow.log' file for more details
```

La pipeline ora fallisce invece di emettere un avviso.
Il codice della pipeline non è cambiato; è cambiata solo la configurazione.

Ripristinate `greetings.csv` al suo stato originale e rimuovete il blocco `validation` prima di continuare:

```csv title="greetings.csv"
greeting,language
Hello,English
Bonjour,French
Holà,Spanish
Ciao,Italian
Hallo,German
```

```groovy title="nextflow.config"
plugins {
    id 'nf-hello@0.5.0'
    id 'nf-schema@2.6.1'
}
```

Sia nf-hello che nf-schema sono plugin di tipo funzione: forniscono funzioni che importate con `include` e chiamate nel codice del vostro flusso di lavoro.
Il prossimo esercizio mostra un tipo diverso di plugin che funziona senza alcuna istruzione `include`.

---

## 3. Usare un plugin observer: nf-co2footprint

Non tutti i plugin forniscono funzioni da importare.
Il plugin [nf-co2footprint](https://github.com/nextflow-io/nf-co2footprint) usa un **trace observer** per monitorare l'utilizzo delle risorse della vostra pipeline e stimarne l'impronta di carbonio.
Non è necessario modificare alcun codice della pipeline; basta aggiungerlo alla configurazione.

### 3.1. Aggiungere nf-co2footprint alla configurazione

Aggiornate `nextflow.config`:

=== "Dopo"

    ```groovy title="nextflow.config" hl_lines="4"
    plugins {
        id 'nf-hello@0.5.0'
        id 'nf-schema@2.6.1'
        id 'nf-co2footprint@1.2.0'
    }
    ```

=== "Prima"

    ```groovy title="nextflow.config"
    plugins {
        id 'nf-hello@0.5.0'
        id 'nf-schema@2.6.1'
    }
    ```

### 3.2. Eseguire la pipeline

```bash
nextflow run hello.nf
```

Il plugin produce diversi messaggi INFO e WARN durante l'esecuzione.
Questi sono normali per un piccolo esempio eseguito su una macchina locale:

```console title="Output (partial)"
nf-co2footprint plugin  ~  version 1.2.0
WARN - [nf-co2footprint] Target zone null not found. Attempting to retrieve carbon intensity for fallback zone GLOBAL.
INFO - [nf-co2footprint] Using fallback carbon intensity from GLOBAL from CI table: 480.0 gCO₂eq/kWh.
WARN - [nf-co2footprint] Executor 'null' not mapped.
WARN - [nf-co2footprint] Fallback to: `machineType = null`, `pue = 1.0`. ...
...
WARN - [nf-co2footprint] No CPU model detected. Using default CPU power draw value (11.41 W).
WARN - [nf-co2footprint] 🔁 Requested memory is null for task 2. Using maximum consumed memory/`peak_rss` (0 GB) for CO₂e footprint computation.
```

Gli avvisi riguardanti zona, executor, modello CPU e memoria appaiono perché il plugin non riesce a rilevare i dettagli hardware completi di un ambiente di formazione locale.
In un ambiente di produzione (ad es., un cluster HPC o cloud), questi valori sarebbero disponibili e le stime più accurate.

Alla fine, cercate una riga simile a:

```console title="Output (partial)"
🌱 The workflow run used 126.76 uWh of electricity, resulting in the release of 60.84 ug of CO₂ equivalents into the atmosphere.
```

(I vostri numeri saranno diversi.)

### 3.3. Visualizzare il report

Il plugin genera file di output nella vostra directory di lavoro:

```bash
ls co2footprint_*
```

```console title="Output"
co2footprint_report_<timestamp>.html
co2footprint_summary_<timestamp>.txt
co2footprint_trace_<timestamp>.txt
```

Esaminate il riepilogo:

```bash
cat co2footprint_summary_*.txt
```

```console title="Output"
Total CO₂e footprint measures of this workflow run (including cached tasks):
  CO₂e emissions: 60.84 ug
  Energy consumption: 126.76 uWh
  CO₂e emissions (market): -

Which equals:
  - 3.48E-7 km travelled by car
  - It takes one tree 0.17s to sequester the equivalent amount of CO₂ from the atmosphere
  - 1.22E-7 % of a flight from Paris to London
```

(I vostri numeri saranno diversi.)

La prima sezione mostra i dati grezzi di energia ed emissioni.
La sezione "Which equals" mette quei numeri in prospettiva convertendoli in equivalenti familiari.
Il riepilogo include anche una sezione che elenca le opzioni di configurazione del plugin e una citazione all'articolo di ricerca [Green Algorithms](https://doi.org/10.1002/advs.202100707) su cui si basa il metodo di calcolo.

### 3.4. Configurare il plugin

L'avviso "Target zone null" della sezione 3.2 è apparso perché il plugin non aveva nessuna posizione geografica configurata.
Il plugin nf-co2footprint definisce uno scope di configurazione `co2footprint` dove potete impostare la vostra posizione geografica.

Aggiungete un blocco `co2footprint` a `nextflow.config`:

=== "Dopo"

    ```groovy title="nextflow.config" hl_lines="7-9"
    plugins {
        id 'nf-hello@0.5.0'
        id 'nf-schema@2.6.1'
        id 'nf-co2footprint@1.2.0'
    }

    co2footprint {
        location = 'GB'
    }
    ```

=== "Prima"

    ```groovy title="nextflow.config"
    plugins {
        id 'nf-hello@0.5.0'
        id 'nf-schema@2.6.1'
        id 'nf-co2footprint@1.2.0'
    }
    ```

!!! tip "Suggerimento"

    Usate il codice del vostro paese se preferite (ad es., `'US'`, `'DE'`, `'FR'`).

Eseguite la pipeline:

```bash
nextflow run hello.nf
```

```console title="Output (partial)"
INFO - [nf-co2footprint] Using fallback carbon intensity from GB from CI table: 163.92 gCO₂eq/kWh.
```

L'avviso sulla zona è scomparso.
Il plugin ora usa l'intensità di carbonio specifica per GB (163.92 gCO₂eq/kWh) invece del valore globale di fallback (480.0 gCO₂eq/kWh).

!!! note "Nota"

    Potreste anche vedere un messaggio `WARN: Unrecognized config option 'co2footprint.location'`.
    È puramente cosmetico e può essere ignorato tranquillamente; il plugin legge comunque il valore correttamente.

Nella Parte 6, creerete uno scope di configurazione per il vostro plugin.

Questo plugin funziona interamente attraverso il meccanismo observer, agganciandosi agli eventi del ciclo di vita del flusso di lavoro per raccogliere metriche sulle risorse e generare il suo report quando la pipeline è completata.

Avete ora provato plugin di tipo funzione (importati con `include`) e un plugin observer (attivato solo tramite configurazione).
Questi sono i due tipi di estensione più comuni, ma come mostra la tabella nella sezione 1, i plugin possono anche aggiungere executor e filesystem.

---

## 4. Scoprire i plugin

Il [Nextflow Plugin Registry](https://registry.nextflow.io/) è il punto centrale per trovare i plugin disponibili.

![La pagina del plugin nf-hello su registry.nextflow.io](img/plugin-registry-nf-hello.png)

Ogni pagina di plugin mostra la sua descrizione, le versioni disponibili, le istruzioni di installazione e i link alla documentazione.

---

## 5. Prepararsi allo sviluppo di plugin

Le sezioni seguenti (Parti 2-6) usano un file di pipeline separato, `greet.nf`, che dipende da nf-schema ma non da nf-hello o nf-co2footprint.

Aggiornate `nextflow.config` per mantenere solo nf-schema:

```groovy title="nextflow.config"
// Configurazione per gli esercizi di sviluppo dei plugin
plugins {
    id 'nf-schema@2.6.1'
}
```

Rimuovete i file di output di co2footprint:

```bash
rm -f co2footprint_*
```

Il file `hello.nf` conserva il lavoro della Parte 1 come riferimento; d'ora in poi lavorerete con `greet.nf`.

---

## Takeaway

Avete usato tre plugin diversi:

- **nf-hello**: Un plugin di tipo funzione che fornisce `randomString`, importato con `include`
- **nf-schema**: Un plugin di tipo funzione che fornisce `samplesheetToList` per il parsing CSV con validazione dello schema
- **nf-co2footprint**: Un plugin observer che monitora automaticamente l'utilizzo delle risorse, senza necessità di `include`

Pattern chiave:

- I plugin vengono dichiarati in `nextflow.config` con `#!groovy plugins { id 'plugin-name@version' }`
- I plugin di tipo funzione richiedono `#!groovy include { function } from 'plugin/plugin-id'`
- I plugin observer funzionano automaticamente una volta dichiarati nella configurazione
- I plugin possono definire scope di configurazione (ad es., `#!groovy validation {}`, `#!groovy co2footprint {}`) per personalizzare il comportamento
- Il [Nextflow Plugin Registry](https://registry.nextflow.io/) elenca i plugin disponibili

---

## Cosa c'è dopo?

Le sezioni seguenti mostrano come costruire il proprio plugin.
Se non siete interessati allo sviluppo di plugin, potete fermarvi qui o passare direttamente al [Riepilogo](summary.md).

[Continua alla Parte 2 :material-arrow-right:](02_create_project.md){ .md-button .md-button--primary }
