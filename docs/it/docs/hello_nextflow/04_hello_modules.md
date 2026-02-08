# Parte 4: Hello Modules

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduzione assistita da IA - [scopri di più e suggerisci miglioramenti](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/43Ot-f0iOME?si=0AWnXB7xqHAzJdJV&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1&amp;cc_lang_pref=it" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } Guarda [l'intera playlist](https://youtube.com/playlist?list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&si=eF7cLR62goy-lc6n) sul canale YouTube di Nextflow.

:green_book: La trascrizione del video è disponibile [qui](./transcripts/04_hello_modules.md).
///

Questa sezione tratta come organizzare il codice del workflow per rendere lo sviluppo e la manutenzione della pipeline più efficienti e sostenibili.
Nello specifico, dimostreremo come usare i [**moduli**](https://nextflow.io/docs/latest/module.html).

In Nextflow, un **modulo** è un file di codice autonomo, spesso incapsulante una singola definizione di processo.
Per usare un modulo in un workflow, basta aggiungere una singola istruzione `include` al file di codice del workflow; poi potete integrare il processo nel workflow nello stesso modo in cui fareste normalmente.
Questo rende possibile riutilizzare le definizioni dei processi in più workflow senza produrre copie multiple del codice.

Quando abbiamo iniziato a sviluppare il nostro workflow, abbiamo scritto tutto in un singolo file di codice.
Ora sposteremo i processi in moduli individuali.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/modules.svg"
</figure>

Questo renderà il nostro codice più condivisibile, flessibile e manutenibile.

??? info "Come iniziare da questa sezione"

    Questa sezione del corso presuppone che abbiate completato le Parti 1-3 del corso [Hello Nextflow](./index.md), ma se avete familiarità con i concetti base trattati in quelle sezioni, potete iniziare da qui senza fare nulla di speciale.

---

## 0. Riscaldamento: Eseguire `hello-modules.nf`

Useremo lo script del workflow `hello-modules.nf` come punto di partenza.
È equivalente allo script prodotto seguendo la Parte 3 di questo corso di formazione, tranne che abbiamo cambiato le destinazioni dell'output:

```groovy title="hello-modules.nf" linenums="37" hl_lines="3 7 11 15"
output {
    first_output {
        path 'hello_modules'
        mode 'copy'
    }
    uppercased {
        path 'hello_modules'
        mode 'copy'
    }
    collected {
        path 'hello_modules'
        mode 'copy'
    }
    batch_report {
        path 'hello_modules'
        mode 'copy'
    }
}
```

Solo per assicurarci che tutto funzioni, eseguite lo script una volta prima di apportare modifiche:

```bash
nextflow run hello-modules.nf
```

??? success "Output del comando"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-modules.nf` [hopeful_avogadro] DSL2 - revision: b09af1237d

    executor >  local (7)
    [0f/8795c9] sayHello (3)       [100%] 3 of 3 ✔
    [6a/eb2510] convertToUpper (3) [100%] 3 of 3 ✔
    [af/479117] collectGreetings   [100%] 1 of 1 ✔
    ```

Come in precedenza, troverete i file di output nella directory specificata nel blocco `output` (qui, `results/hello_modules/`).

??? abstract "Directory contents"

    ```console
    results/hello_modules/
    ├── Bonjour-output.txt
    ├── COLLECTED-batch-output.txt
    ├── Hello-output.txt
    ├── Holà-output.txt
    ├── batch-report.txt
    ├── UPPER-Bonjour-output.txt
    ├── UPPER-Hello-output.txt
    └── UPPER-Holà-output.txt
    ```

Se tutto ha funzionato, siete pronti a imparare come modularizzare il codice del workflow.

---

## 1. Creare una directory per memorizzare i moduli

È buona pratica memorizzare i moduli in una directory specifica.
Potete chiamare quella directory come volete, ma la convenzione è chiamarla `modules/`.

```bash
mkdir modules
```

---

## 2. Creare un modulo per `sayHello()`

Nella sua forma più semplice, trasformare un processo esistente in un modulo è poco più di un'operazione di copia-incolla.
Creeremo uno stub di file per il modulo, copieremo il codice rilevante e poi lo cancelleremo dal file del workflow principale.

Poi tutto ciò che dovremo fare è aggiungere un'istruzione `include` in modo che Nextflow sappia di importare il codice rilevante a runtime.

### 2.1. Creare uno stub di file per il nuovo modulo

Creiamo un file vuoto per il modulo chiamato `sayHello.nf`.

```bash
touch modules/sayHello.nf
```

Questo ci dà un posto dove mettere il codice del processo.

### 2.2. Spostare il codice del processo `sayHello` nel file del modulo

Copiate l'intera definizione del processo dal file del workflow al file del modulo.

```groovy title="modules/sayHello.nf" linenums="1"
/*
 * Use echo to print 'Hello World!' to a file
 */
process sayHello {

    input:
    val greeting

    output:
    path "${greeting}-output.txt"

    script:
    """
    echo '${greeting}' > '${greeting}-output.txt'
    """
}
```

Una volta fatto, cancellate la definizione del processo dal file del workflow.

### 2.3. Aggiungere una dichiarazione di include prima del blocco workflow

La sintassi per includere un processo da un modulo è abbastanza semplice:

```groovy title="Sintassi: dichiarazione di include"
include { <PROCESS_NAME> } from '<percorso_al_modulo>'
```

Inseriamola sopra il blocco `params` e compiliamola appropriatamente.

=== "Dopo"

    ```groovy title="hello-modules.nf" linenums="44" hl_lines="1 2"
    // Include i moduli
    include { sayHello } from './modules/sayHello.nf'

    /*
    * Pipeline parameters
    */
    params {
        greeting: Path = 'data/greetings.csv'
        batch: String = 'batch'
    }
    ```

=== "Prima"

    ```groovy title="hello-modules.nf" linenums="44"
    /*
    * Pipeline parameters
    */
    params {
        greeting: Path = 'data/greetings.csv'
        batch: String = 'batch'
    }
    ```

Vedete che abbiamo inserito il nome del processo, `sayHello`, e il percorso al file contenente il codice del modulo, `./modules/sayHello.nf`.

### 2.4. Eseguire il workflow

Stiamo eseguendo il workflow essenzialmente con lo stesso codice e input di prima, quindi eseguiamo con il flag `-resume` e vediamo cosa succede.

```bash
nextflow run hello-modules.nf -resume
```

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-modules.nf` [romantic_poisson] DSL2 - revision: 96edfa9ad3

    [f6/cc0107] sayHello (1)       | 3 of 3, cached: 3 ✔
    [3c/4058ba] convertToUpper (2) | 3 of 3, cached: 3 ✔
    [1a/bc5901] collectGreetings   | 1 of 1, cached: 1 ✔
    ```

Questo dovrebbe funzionare molto velocemente perché tutto è in cache.
Sentitevi liberi di controllare gli output pubblicati.

Nextflow ha riconosciuto che è ancora tutto lo stesso lavoro da fare, anche se il codice è diviso in più file.

### Takeaway

Sapete come estrarre un processo in un modulo locale e sapete che fare questo non compromette la ripristinabilità del workflow.

### Cosa c'è dopo?

Fare pratica creando più moduli.
Una volta che ne avete fatto uno, potete farne un milione di più...
Ma per ora ne facciamo solo altri due.

---

## 3. Modularizzare il processo `convertToUpper()`

### 3.1. Creare uno stub di file per il nuovo modulo

Create un file vuoto per il modulo chiamato `convertToUpper.nf`.

```bash
touch modules/convertToUpper.nf
```

### 3.2. Spostare il codice del processo `convertToUpper` nel file del modulo

Copiate l'intera definizione del processo dal file del workflow al file del modulo.

```groovy title="modules/convertToUpper.nf" linenums="1"
/*
 * Use a text replacement tool to convert the greeting to uppercase
 */
process convertToUpper {

    input:
    path input_file

    output:
    path "UPPER-${input_file}"

    script:
    """
    cat '${input_file}' | tr '[a-z]' '[A-Z]' > 'UPPER-${input_file}'
    """
}
```

Una volta fatto, cancellate la definizione del processo dal file del workflow.

### 3.3. Aggiungere una dichiarazione di include prima del blocco `params`

Inserite la dichiarazione di include sopra il blocco `params` e compilatela appropriatamente.

=== "Dopo"

    ```groovy title="hello-modules.nf" linenums="23" hl_lines="3"
    // Include i moduli
    include { sayHello } from './modules/sayHello.nf'
    include { convertToUpper } from './modules/convertToUpper.nf'

    /*
    * Pipeline parameters
    */
    params {
        greeting: Path = 'data/greetings.csv'
        batch: String = 'batch'
    }
    ```

=== "Prima"

    ```groovy title="hello-modules.nf" linenums="23"
    // Include i moduli
    include { sayHello } from './modules/sayHello.nf'

    /*
    * Pipeline parameters
    */
    params {
        greeting: Path = 'data/greetings.csv'
        batch: String = 'batch'
    }
    ```

Questo dovrebbe iniziare a sembrare molto familiare.

### 3.4. Eseguire di nuovo il workflow

Eseguite questo con il flag `-resume`.

```bash
nextflow run hello-modules.nf -resume
```

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-modules.nf` [nauseous_heisenberg] DSL2 - revision: a04a9f2da0

    [c9/763d42] sayHello (3)       | 3 of 3, cached: 3 ✔
    [60/bc6831] convertToUpper (3) | 3 of 3, cached: 3 ✔
    [1a/bc5901] collectGreetings   | 1 of 1, cached: 1 ✔
    ```

Questo dovrebbe ancora produrre lo stesso output di prima.

Due fatti, ne manca ancora uno!

---

## 4. Modularizzare il processo `collectGreetings()`

### 4.1. Creare uno stub di file per il nuovo modulo

Create un file vuoto per il modulo chiamato `collectGreetings.nf`.

```bash
touch modules/collectGreetings.nf
```

### 4.2. Spostare il codice del processo `collectGreetings` nel file del modulo

Copiate l'intera definizione del processo dal file del workflow al file del modulo.

```groovy title="modules/collectGreetings.nf" linenums="1"
/*
 * Collect uppercase greetings into a single output file
 */
process collectGreetings {

    input:
    path input_files
    val batch_name

    output:
    path "COLLECTED-${batch_name}-output.txt", emit: outfile
    path "${batch_name}-report.txt", emit: report

    script:
    count_greetings = input_files.size()
    """
    cat ${input_files} > 'COLLECTED-${batch_name}-output.txt'
    echo 'There were ${count_greetings} greetings in this batch.' > '${batch_name}-report.txt'
    """
}
```

Una volta fatto, cancellate la definizione del processo dal file del workflow.

### 4.3. Aggiungere una dichiarazione di include prima del blocco `params`

Inserite la dichiarazione di include sopra il blocco `params` e compilatela appropriatamente.

=== "Dopo"

    ```groovy title="hello-modules.nf" linenums="3" hl_lines="4"
    // Include i moduli
    include { sayHello } from './modules/sayHello.nf'
    include { convertToUpper } from './modules/convertToUpper.nf'
    include { collectGreetings } from './modules/collectGreetings.nf'

    /*
    * Pipeline parameters
    */
    params {
        greeting: Path = 'data/greetings.csv'
        batch: String = 'batch'
    }
    ```

=== "Prima"

    ```groovy title="hello-modules.nf" linenums="3"
    // Include i moduli
    include { sayHello } from './modules/sayHello.nf'
    include { convertToUpper } from './modules/convertToUpper.nf'

    /*
    * Pipeline parameters
    */
    params {
        greeting: Path = 'data/greetings.csv'
        batch: String = 'batch'
    }
    ```

Ultimo!

### 4.4. Eseguire il workflow

Eseguite questo con il flag `-resume`.

```bash
nextflow run hello-modules.nf -resume
```

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-modules.nf` [friendly_coulomb] DSL2 - revision: 7aa2b9bc0f

    [f6/cc0107] sayHello (1)       | 3 of 3, cached: 3 ✔
    [3c/4058ba] convertToUpper (2) | 3 of 3, cached: 3 ✔
    [1a/bc5901] collectGreetings   | 1 of 1, cached: 1 ✔
    ```

Questo dovrebbe ancora produrre lo stesso output di prima.

### Takeaway

Sapete come modularizzare più processi in un workflow.

Congratulazioni, avete fatto tutto questo lavoro e non è cambiato assolutamente nulla nel funzionamento della pipeline!

A parte gli scherzi, ora il vostro codice è più modulare, e se decidete di scrivere un'altra pipeline che chiama uno di quei processi, vi basta digitare una breve istruzione `include` per usare il modulo rilevante.
Questo è meglio che fare copia-incolla del codice, perché se in seguito decidete di migliorare il modulo, tutte le vostre pipeline erediteranno i miglioramenti.

### Cosa c'è dopo?

Prendetevi una breve pausa se ne avete voglia.

Quando siete pronti, passate alla [**Parte 5: Hello Containers**](./05_hello_containers.md) per imparare come usare i container per gestire le dipendenze software in modo più conveniente e riproducibile.

---

## Quiz

<quiz>
Cos'è un modulo in Nextflow?
- [ ] Un file di configurazione
- [x] Un file autonomo che può contenere definizioni di processi
- [ ] Una definizione di workflow
- [ ] Un operatore di canali

Per approfondire: [2. Creare un modulo per `sayHello()`](#2-creare-un-modulo-per-sayhello)
</quiz>

<quiz>
Quale convenzione viene tipicamente usata per memorizzare i file dei moduli?
- [ ] Nella stessa directory del workflow
- [ ] In una directory `bin/`
- [x] In una directory `modules/`
- [ ] In una directory `lib/`

Per approfondire: [1. Creare una directory per memorizzare i moduli](#1-creare-una-directory-per-memorizzare-i-moduli)
</quiz>

<quiz>
Qual è la sintassi corretta per usare un modulo?

- [ ] `#!groovy import { SAYHELLO } from './modules/sayhello.nf'`
- [ ] `#!groovy require { SAYHELLO } from './modules/sayhello.nf'`
- [x] `#!groovy include { SAYHELLO } from './modules/sayhello.nf'`
- [ ] `#!groovy load { SAYHELLO } from './modules/sayhello.nf'`

Per approfondire: [2.3. Aggiungere una dichiarazione di include](#23-aggiungere-una-dichiarazione-di-include-prima-del-blocco-workflow)
</quiz>

<quiz>
Cosa succede alla funzionalità `-resume` quando si usano i moduli?
- [ ] Non funziona più
- [ ] Richiede configurazione aggiuntiva
- [x] Funziona come prima
- [ ] Funziona solo per i moduli locali
</quiz>

<quiz>
Quali sono i vantaggi dell'uso dei moduli? (Selezionate tutte le risposte applicabili)
- [x] Riutilizzabilità del codice tra workflow
- [x] Manutenzione più facile
- [x] Migliore organizzazione del codice del workflow
- [ ] Velocità di esecuzione più rapida
</quiz>
