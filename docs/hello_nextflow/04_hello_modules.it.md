# Parte 4: Hello Modules

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/Xxp_menS0E8?si=0AWnXB7xqHAzJdJV&amp;list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } Guarda [tutta la playlist](https://www.youtube.com/playlist?list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik) sul canale Youtube Nextflow.

:green_book: La trascrizione di questo video è disponibile [qui](./transcripts/03_hello_workflow.md).
///

Questa sezione spiega come organizzare il codice del workflow per rendere più efficiente e sostenibile lo sviluppo e la manutenzione della pipeline.
In particolare, dimostreremo come utilizzare i **moduli**.

In Nextflow, un **modulo** è una singola definizione di processo incapsulata in un file di codice indipendente.
Per utilizzare un modulo in un workflow, basta aggiungere una singola riga di dichiarazione di importazione al file di codice del workflow; quindi è possibile integrare il processo nel workflow come si farebbe normalmente.

Quando abbiamo iniziato a sviluppare il nostro workflow, abbiamo inserito tutto in un unico file di codice.

La suddivisione dei processi in singoli moduli consente di riutilizzare le definizioni dei processi in più workflow senza produrre più copie del codice.
Questo rende il codice più condivisibile, flessibile e manutenibile.

!!! note

    È anche possibile incapsulare una sezione di un workflow come un "sottoflusso" che può essere importato in una pipeline più ampia, ma ciò esula dallo scopo di questo corso.

---

## 0. Riscaldamento: Eseguire `hello-modules.nf`

Utilizzeremo lo script del workflow `hello-modules.nf` come punto di partenza.
È equivalente allo script prodotto lavorando alla Parte 3 di questo corso di formazione.

Per assicurarsi che tutto funzioni, eseguire lo script una volta prima di apportare qualsiasi modifica:

```bash
nextflow run hello-modules.nf
```

```console title="Output"
 N E X T F L O W   ~  version 24.10.0

Launching `hello-modules.nf` [festering_nobel] DSL2 - revision: eeca64cdb1

executor >  local (7)
[25/648bdd] sayHello (2)       | 3 of 3 ✔
[60/bc6831] convertToUpper (1) | 3 of 3 ✔
[1a/bc5901] collectGreetings   | 1 of 1 ✔
there were 3 greetings in this batch
```

Come in precedenza, i file di output si trovano nella directory `results` (specificata dalla direttiva `publishDir`).

```console title="Directory contents"
results
├── Bonjour-output.txt
├── COLLECTED-output.txt
├── COLLECTED-test-batch-output.txt
├── COLLECTED-trio-output.txt
├── Hello-output.txt
├── Holà-output.txt
├── UPPER-Bonjour-output.txt
├── UPPER-Hello-output.txt
└── UPPER-Holà-output.txt
```

!!! note

    Potrebbe anche essere rimasto un file chiamato `output.txt` se si è lavorato alla Parte 2 nello stesso ambiente.

Se questo ha funzionato, siete pronti per imparare a modularizzare il codice del workflow.

---

## 1. Creare una cartella per memorizzare i moduli

È buona norma memorizzare i moduli in una directory specifica.
Si può chiamare questa cartella come si vuole, ma la convenzione è di chiamarla `modules/`.

```bash
mkdir modules
```

!!! note

    Qui mostriamo come usare i moduli locali, cioè i moduli memorizzati localmente nello stesso repository del resto del codice del workflow, in contrasto con i moduli remoti, che sono memorizzati in altri repository (remoti). Per maggiori informazioni sui moduli remoti, si veda la [documentazione](https://www.nextflow.io/docs/latest/module.html).

---

## 2. Creare un modulo per `sayHello()`

Nella sua forma più semplice, trasformare un processo esistente in un modulo è poco più di un'operazione di copia-incolla.
Creeremo un file stub per il modulo, copieremo il codice pertinente e lo cancelleremo dal file principale del workflow.

A questo punto, basterà aggiungere una dichiarazione di importazione, in modo che Nextflow sappia che deve inserire il codice in questione in fase di esecuzione.

### 2.1. Creare un file stub per il nuovo modulo

Creiamo un file vuoto per il modulo, chiamato `sayHello.nf`.

```bash
touch modules/sayHello.nf
```

Questo ci dà un posto dove mettere il codice del processo.

### 2.2. Spostare il codice del processo `sayHello' nel file del modulo

Copiare l'intera definizione del processo dal file del workflow al file del modulo, assicurandosi di copiare anche lo shebang `#!/usr/bin/env nextflow`.

```groovy title="modules/sayHello.nf" linenums="1"
#!/usr/bin/env nextflow

/*
 * Usa echo per stampare 'Hello World!' A un file.
 */
process sayHello {

    publishDir 'results', mode: 'copy'

    input:
        val greeting

    output:
        path "${greeting}-output.txt"

    script:
    """
    echo '$greeting' > '$greeting-output.txt'
    """
}
```

Una volta fatto ciò, eliminate la definizione del processo dal file del workflow, ma assicuratevi di lasciare lo shebang al suo posto.

### 2.3. Aggiungere una dichiarazione di importazione prima del blocco del workflow

La sintassi per importare un modulo locale è abbastanza semplice:

```groovy title="Syntax: Import declaration"
include { <MODULE_NAME> } from '<path_to_module>'
```

Inseriamo questo blocco sopra il blocco del workflow e compiliamolo in modo appropriato.

=== "Dopo"

    ```groovy title="hello-modules.nf" linenums="50" hl_lines="1 2"
    // Include modules
    include { sayHello } from './modules/sayHello.nf'

    workflow {
    ```

=== "Prima"

    ```groovy title="hello-modules.nf" linenums="50"
    workflow {
    ```

### 2.4. Eseguite il workflow per verificare che faccia la stessa cosa di prima

Stiamo eseguendo il workflow essenzialmente con lo stesso codice e gli stessi input di prima, quindi eseguiamolo con il flag `resume` e vediamo cosa succede.

```bash
nextflow run hello-modules.nf -resume
```

L'esecuzione è molto rapida perché tutto viene memorizzato nella cache.

```console title="Output" linenums="1"
 N E X T F L O W   ~  version 24.10.0

Launching `hello-modules.nf` [romantic_poisson] DSL2 - revision: 96edfa9ad3

[f6/cc0107] sayHello (1)       | 3 of 3, cached: 3 ✔
[3c/4058ba] convertToUpper (2) | 3 of 3, cached: 3 ✔
[1a/bc5901] collectGreetings   | 1 of 1, cached: 1 ✔
there were 3 greetings in this batch
```

Nextflow ha riconosciuto che il lavoro da fare è sempre lo stesso, anche se il codice è suddiviso in più file.

### Conclusione

Sapete come estrarre un processo in un modulo locale e sapete che questo non rompe la riprendibilità del workflow.

### Cosa c'è dopo?

Esercitatevi a creare altri moduli.
Una volta che ne avete fatto uno, potete farne un milione...
Ma per ora facciamone solo altri due.

---

## 3. Modularizzare il processo `convertToUpper()`

### 3.1. Creare un file stub per il nuovo modulo

Creare un file vuoto per il modulo, chiamato `convertToUpper.nf`.

```bash
touch modules/convertToUpper.nf
```

### 3.2. Spostare il codice del processo `convertToUpper` nel file del modulo

Copiare l'intera definizione del processo dal file del workflow al file del modulo, assicurandosi di copiare anche lo shebang `#!/usr/bin/env nextflow`.

```groovy title="modules/convertToUpper.nf" linenums="1"
#!/usr/bin/env nextflow

/*
 * Utilizzare uno strumento di sostituzione del testo per convertire il saluto in maiuscolo.
 */
process convertToUpper {

    publishDir 'results', mode: 'copy'

    input:
        path input_file

    output:
        path "UPPER-${input_file}"

    script:
    """
    cat '$input_file' | tr '[a-z]' '[A-Z]' > 'UPPER-${input_file}'
    """
}
```

Una volta fatto ciò, eliminate la definizione del processo dal file del workflow, ma assicuratevi di lasciare lo shebang al suo posto.

### 3.3. Aggiungere una dichiarazione di importazione prima del blocco del workflow

Inserite la dichiarazione di importazione sopra il blocco del workflow e compilatela in modo appropriato.

=== "Dopo"

    ```groovy title="hello-modules.nf" linenums="31" hl_lines="3"
    // Include modules
    include { sayHello } from './modules/sayHello.nf'
    include { convertToUpper } from './modules/convertToUpper.nf'

    workflow {
    ```

=== "Prima"

    ```groovy title="hello-modules.nf" linenums="31"
    // Include modules
    include { sayHello } from './modules/sayHello.nf'

    workflow {
    ```

### 3.4. Eseguite il workflow per verificare che faccia la stessa cosa di prima

Eseguire questa operazione con il flag `-resume`.

```bash
nextflow run hello-modules.nf -resume
```

Il risultato dovrebbe essere lo stesso di quello precedente.

```console title="Output" linenums="1"
 N E X T F L O W   ~  version 24.10.0

Launching `hello-modules.nf` [nauseous_heisenberg] DSL2 - revision: a04a9f2da0

[c9/763d42] sayHello (3)       | 3 of 3, cached: 3 ✔
[60/bc6831] convertToUpper (3) | 3 of 3, cached: 3 ✔
[1a/bc5901] collectGreetings   | 1 of 1, cached: 1 ✔
There were 3 greetings in this batch
```

Due fatti, uno ancora da fare!

---

## 4. Modularizzare il processo `collectGreetings()`

### 4.1. Creare un file stub per il nuovo modulo

Creare un file vuoto per il modulo, chiamato `collectGreetings.nf`.

```bash
touch modules/collectGreetings.nf
```

### 4.2. Spostare il codice del processo `collectGreetings` nel file del modulo###

Copiare l'intera definizione del processo dal file del workflow al file del modulo, assicurandosi di copiare anche lo shebang `#!/usr/bin/env nextflow`.

```groovy title="modules/collectGreetings.nf" linenums="1"
#!/usr/bin/env nextflow

/*
 * Raccogli i saluti maiuscoli in un singolo file di output
 */
process collectGreetings {

    publishDir 'results', mode: 'copy'

    input:
        path input_files
        val batch_name

    output:
        path "COLLECTED-${batch_name}-output.txt" , emit: outfile
        val count_greetings , emit: count

    script:
        count_greetings = input_files.size()
    """
    cat ${input_files} > 'COLLECTED-${batch_name}-output.txt'
    """
}
```

Una volta fatto ciò, eliminate la definizione del processo dal file del workflow, ma assicuratevi di lasciare lo shebang al suo posto.

### 4.3. Aggiungere una dichiarazione di importazione prima del blocco del workflow###

Inserite la dichiarazione di importazione sopra il blocco del workflow e compilatela in modo appropriato.

=== "Dopo"

    ```groovy title="hello-modules.nf" linenums="9" hl_lines="4"
    // Include modules
    include { sayHello } from './modules/sayHello.nf'
    include { convertToUpper } from './modules/convertToUpper.nf'
    include { collectGreetings } from './modules/collectGreetings.nf'

    workflow {
    ```

=== "Prima"

    ```groovy title="hello-modules.nf" linenums="9"
    // Include modules
    include { sayHello } from './modules/sayHello.nf'
    include { convertToUpper } from './modules/convertToUpper.nf'

    workflow {
    ```

### 4.4. Eseguite il workflow per verificare che faccia la stessa cosa di prima

Eseguire questa operazione con il flag `-resume`.

```bash
nextflow run hello-modules.nf -resume
```

Il risultato dovrebbe essere lo stesso di quello precedente.

```console title="Output" linenums="1"
 N E X T F L O W   ~  version 24.10.0

Launching `hello-modules.nf` [friendly_coulomb] DSL2 - revision: 7aa2b9bc0f

[f6/cc0107] sayHello (1)       | 3 of 3, cached: 3 ✔
[3c/4058ba] convertToUpper (2) | 3 of 3, cached: 3 ✔
[1a/bc5901] collectGreetings   | 1 of 1, cached: 1 ✔
There were 3 greetings in this batch
```

### Conclusione

Sapete come modulare più processi in un workflow.

Congratulazioni, avete fatto tutto questo lavoro e non è cambiato assolutamente nulla nel funzionamento della pipeline!

A parte gli scherzi, ora il codice è più modulare e se si decide di scrivere un'altra pipeline che richiama uno di questi processi, è sufficiente digitare una breve dichiarazione di importazione per utilizzare il modulo corrispondente.
È meglio che copiare il codice, perché se in seguito si decide di migliorare il modulo, tutte le pipeline erediteranno i miglioramenti.

### Cosa c'è dopo?

Fate una breve pausa se ne avete voglia.
Quando siete pronti, passate alla Parte 5 per imparare a usare i container per gestire le dipendenze del software in modo più pratico e riproducibile.
