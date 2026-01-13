# Parte 2: Hello Channels

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/lJ41WMMm44M?si=xCItHLiOQWqoqBB9&amp;list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } Guarda [tutta la playlist](https://www.youtube.com/playlist?list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik) sul canale Youtube Nextflow.

:green_book: La trascrizione di questo video è disponibile [qui](./transcripts/02_hello_channels.md).
///

Nella parte 1 di questo corso (Hello World), ti abbiamo mostrato come fornire un input variabile a un processo fornendo direttamente l'input nella chiamata di processo: `sayHello(params.greet)`.
Era un approccio volutamente semplificato.
In pratica, questo approccio ha grandi limitazioni; vale a dire che funziona solo per casi molto semplici in cui vogliamo eseguire il processo solo una volta, su un singolo valore.
Nella maggior parte dei casi d'uso del workflow realistici, vogliamo elaborare più valori (dati sperimentali per più campioni, ad esempio), quindi abbiamo bisogno di un modo più sofisticato per gestire gli input.

Ecco a cosa servono i **canali** di Nextflow.
I canali sono code progettate per gestire gli input in modo efficiente e spostarli da un passaggio all'altro in workflows in più fasi, fornendo parallelismo integrato e molti vantaggi aggiuntivi.

In questa parte del corso, imparerai come utilizzare un canale per gestire più input da una varietà di fonti diverse.
Imparerai anche a usare **operatori** per trasformare i contenuti del canale secondo necessità.

_Per la formazione sull'uso dei canali per collegare i passaggi in un workflow multi-step, vedere la parte 3 di questo corso. _

---

## 0. Riscaldamento: esegui `hello-channels.nf`

Useremo lo script del workflow `hello-channels.nf` come punto di partenza.
È equivalente allo script prodotto lavorando attraverso la Parte 1 di questo corso di formazione.

Solo per assicurarti che tutto funzioni, esegui lo script una volta prima di apportare qualsiasi modifica:

```bash
nextflow run hello-channels.nf --greeting 'Hello Channels!'
```

```console title="Output"
 N E X T F L O W   ~  version 24.10.0

Launching `hello-channels.nf` [insane_lichterman] DSL2 - revision: c33d41f479

executor >  local (1)
[86/9efa08] sayHello | 1 of 1 ✔
```

Come in precedenza, troverai il file di output denominato `output.txt` nella directory `results` (specificata dalla direttiva `publishDir`).

```console title="results/output.txt" linenums="1"
Hello Channels!
```

Se ha funzionato per te, sei pronto a conoscere i canali.

---

## 1. Fornire input variabili tramite un canale esplicitamente

Creeremo un **canale** per passare l'input variabile al processo `sayHello()` invece di fare affidamento sulla gestione implicita, che ha alcune limitazioni.

### 1.1. Crea un canale di input

Ci sono una varietà di **fabbriche di canali** che possiamo utilizzare per creare un canale.
Per mantenere le cose semplici per ora, useremo la fabbrica di canale più semplice, chiamata `channel.of`, che creerà un canale contenente un singolo valore.
Funzionalmente questo sarà simile a come lo avevamo impostato prima, ma invece di avere Nextflow a creare implicitamente un canale, lo stiamo facendo esplicitamente ora.

Questa è la riga di codice che useremo:

```console title="Syntax"
greeting_ch = channel.of('Hello Channels!')
```

Questo crea un canale chiamato `greeting_ch` usando la fabbrica del canale `channel.of()`, che imposta un semplice canale di coda e carica la stringa `'Hello Channels!' ` da usare come valore di saluto.

!!! note

    Stiamo temporaneamente tornando alle stringhe codificate invece di utilizzare un parametro CLI per motivi di leggibilità. Torneremo a utilizzare i parametri CLI una volta che avremo coperto ciò che sta accadendo a livello del canale.

Nel blocco del workflow, aggiungi il codice di fabbrica del canale:

=== "Dopo"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="3 4"
    workflow {

        // create a channel for inputs
        greeting_ch = channel.of('Hello Channels!')

        // emit a greeting
        sayHello(params.greeting)
    }
    ```

=== "Prima"

    ```groovy title="hello-channels.nf" linenums="27"
    workflow {

        // emit a greeting
        sayHello(params.greeting)
    }
    ```

Questo non è ancora funzionale poiché non abbiamo ancora cambiato l'input alla chiamata di processo.

### 1.2. Aggiungi il canale come input alla chiamata di processo

Ora dobbiamo effettivamente collegare il nostro canale appena creato alla chiamata di processo `sayHello()`, sostituendo il parametro CLI che stavamo fornendo direttamente prima.

Nel blocco del workflow, apportare la seguente modifica del codice:

=== "Dopo"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="7"
    workflow {

        // create a channel for inputs
        greeting_ch = channel.of('Hello Channels!')

        // emit a greeting
        sayHello(greeting_ch)
    }
    ```

=== "Prima"

    ```groovy title="hello-channels.nf" linenums="27"
    workflow {

        // create a channel for inputs
        greeting_ch = channel.of('Hello Channels!')

        // emit a greeting
        sayHello(params.greeting)
    }
    ```

Questo dice a Nextflow di eseguire il processo `sayHello` sui contenuti del canale `greeting_ch`.

Ora il nostro workflow è correttamente funzionale; è l'equivalente esplicito di scrivere `sayHello('Hello Channels!') `.

### 1.3. Esegui di nuovo il comando del workflow

Eseguiamolo!

```bash
nextflow run hello-channels.nf
```

Se hai apportato entrambe le modifiche correttamente, dovresti ottenere un'altra esecuzione riuscita:

```console title="Output" linenums="1"
 N E X T F L O W   ~  version 24.10.0

Launching `hello-channels.nf` [nice_heisenberg] DSL2 - revision: 41b4aeb7e9

executor >  local (1)
[3b/f2b109] sayHello (1) | 1 of 1 ✔
```

Puoi controllare la directory dei risultati per assicurarti che il risultato sia ancora lo stesso di prima.

```console title="results/output.txt" linenums="1"
Hello Channels!
```

Finora stiamo solo modificando progressivamente il codice per aumentare la flessibilità del nostro workflow ottenendo allo stesso tempo lo stesso risultato finale.

!!! note

    Può sembrare che stiamo scrivendo più codice senza alcun beneficio tangibile, ma il valore diventerà chiaro non appena inizieremo a gestire più input.

### Conclusioni

Sai come utilizzare una fabbrica di canale di base per fornire un input a un processo.

### Cosa c'è dopo?

Scopri come utilizzare i canali per far in modo che il workflow iteri su più valori di input.

---

## 2. Modifica il workflow per eseguire su più valori di input

I workflows in genere vengono eseguiti su lotti di input che devono essere elaborati in blocco, quindi vogliamo aggiornare il workflow per accettare più valori di input.

### 2.1. Carica più saluti nel canale di input

Convenientemente, la fabbrica di canali `channel.of()` che abbiamo utilizzato è abbastanza felice di accettare più di un valore, quindi non abbiamo affatto bisogno di modificarlo.

Dobbiamo solo caricare più valori nel canale.

#### 2.1.1. Aggiungi altri saluti

Nel blocco del workflow, apportare la seguente modifica del codice:

=== "Dopo"

    ```groovy title="hello-channels.nf" linenums="29" hl_lines="2"
    // create a channel for inputs
    greeting_ch = channel.of('Hello','Bonjour','Holà')
    ```

=== "Prima"

    ```groovy title="hello-channels.nf" linenums="29"
    // create a channel for inputs
    greeting_ch = channel.of('Hello Channels')
    ```

La documentazione ci dice che questo dovrebbe funzionare. Può davvero essere così semplice?

#### 2.1.2. Esegui il comando e guarda l'output del registro

Proviamolo.

```bash
nextflow run hello-channels.nf
```

Certamente sembra funzionare bene:

```console title="Output" linenums="1"
 N E X T F L O W   ~  version 24.10.0

Launching `hello-channels.nf` [suspicious_lamport] DSL2 - revision: 778deadaea

executor >  local (3)
[cd/77a81f] sayHello (3) | 3 of 3 ✔
```

Tuttavia... Questo sembra indicare che sono state effettuate chiamate "3 di 3" per il processo, il che è incoraggiante, ma questo ci mostra solo una singola esecuzione del processo, con un percorso di sottodirectory (`cd/77a81f`).
Cosa sta succedendo?

Per impostazione predefinita, il sistema di registrazione ANSI scrive la registrazione da più chiamate allo stesso processo sulla stessa linea.
Fortunatamente, possiamo disabilitare quel comportamento per vedere l'elenco completo delle chiamate di processo.

#### 2.1.3. Esegui di nuovo il comando con l'opzione `-ansi-log false`

Per espandere la registrazione per visualizzare una riga per chiamata di processo, aggiungi `-ansi-log false` al comando.

```bash
nextflow run hello-channels.nf -ansi-log false
```

Questa volta vediamo tutte e tre le esecuzioni di processo e le loro sottodirectory di lavoro associate elencate nell'output:

```console title="Output" linenums="1"
N E X T F L O W  ~  version 24.10.0
Launching `hello-channels.nf` [pensive_poitras] DSL2 - revision: 778deadaea
[76/f61695] Submitted process > sayHello (1)
[6e/d12e35] Submitted process > sayHello (3)
[c1/097679] Submitted process > sayHello (2)
```

È molto meglio; almeno per un semplice workflow.
Per un workflow complesso o un gran numero di input, avere l'output completo dell'elenco sul terminale potrebbe diventare un po' travolgente, quindi potresti non scegliere di usare `-ansi-log false` in quei casi.

!!! note

    Il modo in cui viene riportato lo stato è un po' diverso tra le due modalità di registrazione.
    In modalità condensata, Nextflow segnala se le chiamate sono state completate con successo o meno.
    In questa modalità espansa, riporta solo che sono stati inviati.

Detto questo, abbiamo un altro problema. Se guardi nella directory `results`, c'è solo un file: `output.txt`!

```console title="Directory contents"
results
└── output.txt
```

Che succede? Non dovremmo aspettarci un file separato per saluto di input, quindi tre file in tutto?
Tutti e tre i saluti sono andati in un unico file?

Puoi controllare il contenuto di `output.txt`; troverai solo uno dei tre, contenente uno dei tre saluti che abbiamo fornito.

```console title="output.txt" linenums="1"
Bonjour
```

Potresti ricordare che abbiamo codificato il nome del file di output per il processo `sayHello`, quindi tutte e tre le chiamate hanno prodotto un file chiamato `output.txt`.
Puoi controllare le sottodirectory di lavoro per ciascuno dei tre processi; ognuno di essi contiene un file chiamato `output.txt` come previsto.

Finché i file di output rimangono lì, isolati dagli altri processi, va bene.
Ma quando la direttiva `publishDir` copia ciascuno di essi nella stessa directory `results`, qualunque sia stato copiato lì per primo viene sovrascritto da quello successivo, e così via.

### 2.2. Assicurati che i nomi dei file di output siano univoci

Possiamo continuare a pubblicare tutti gli output nella stessa directory dei risultati, ma dobbiamo assicurarci che abbiano nomi univoci.
In particolare, dobbiamo modificare il primo processo per generare un nome di file in modo dinamico in modo che i nomi dei file finali siano unici.

Quindi, come facciamo a rendere unici i nomi dei file?
Un modo comune per farlo è utilizzare un pezzo unico di metadati dagli input (ricevuti dal canale di input) come parte del nome del file di output.
Qui, per comodità, useremo solo il saluto stesso poiché è solo una stringa breve e lo antemeremo al nome del file di output di base.

#### 2.2.1. Costruisci un nome di file di output dinamico

Nel blocco di processo, apporta le seguenti modifiche al codice:

=== "Dopo"

    ```groovy title="hello-channels.nf" linenums="6" hl_lines="9 13"
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

=== "Prima"

    ```groovy title="hello-channels.nf" linenums="6"
    process sayHello {

        publishDir 'results', mode: 'copy'

        input:
            val greeting

        output:
            path 'output.txt'

        script:
        """
        echo '$greeting' > output.txt
        """
    }
    ```

Assicurati di sostituire `output.txt` sia nella definizione di output che nel blocco di comando `script:`.

!!! tip

    Nella definizione di output, DEVI usare virgolette doppie attorno all'espressione del nome del file di output (NON virgolette singole), altrimenti fallirà.

Questo dovrebbe produrre un nome di file di output univoco ogni volta che viene chiamato il processo, in modo che possa essere distinto dagli output di altre iterazioni dello stesso processo nella directory di output.

#### 2.2.2. Esegui il workflow

Eseguiamolo:

```bash
nextflow run hello-channels.nf
```

Tornando alla vista di riepilogo, l'output ha di nuovo questo aspetto:

```console title="Output" linenums="1"
 N E X T F L O W   ~  version 24.10.0

Launching `hello-channels.nf` [astonishing_bell] DSL2 - revision: f57ff44a69

executor >  local (3)
[2d/90a2e2] sayHello (1) | 3 of 3 ✔
```

È importante sottolineare che ora abbiamo tre nuovi file oltre a quello che avevamo già nella directory `risultati`:

```console title="Directory contents"
results
├── Bonjour-output.txt
├── Hello-output.txt
├── Holà-output.txt
└── output.txt
```

Ognuno di loro ha i contenuti previsti:

```console title="Bonjour-output.txt" linenums="1"
Bonjour
```

```console title="Hello-output.txt" linenums="1"
Hello
```

```console title="Holà-output.txt" linenums="1"
Holà
```

Successo! Ora possiamo aggiungere tutti i saluti che vorremmo senza preoccuparci che i file di output vengano sovrascritti.

!!! note

    In pratica, nominare i file in base ai dati di input stessi è quasi sempre impraticabile.
    Il modo migliore per generare nomi di file dinamici è passare i metadati a un processo insieme ai file di input.
    I metadati sono in genere forniti tramite un "foglio campione" o equivalenti.
    Imparerai come farlo più avanti nel tuo corso di formazione Nextflow.

### Conclusioni

Sai come alimentare più elementi di input attraverso un canale.

### Cosa c'è dopo?

Impara a usare un operatore per trasformare i contenuti di un canale.

---

## 3. Usa un operatore per trasformare il contenuto di un canale

In Nextflow, [operatori](https://www.nextflow.io/docs/latest/reference/operator.html) ci consentono di trasformare il contenuto di un canale.

Ti abbiamo appena mostrato come gestire più elementi di input che sono stati codificati direttamente nella fabbrica del canale.
E se volessimo fornire quegli input multipli in una forma diversa?

Ad esempio, immagina di impostare una variabile di input contenente una serie di elementi come questa:

`greetings_array = ['Ciao','Bonjour','Holà']`

Possiamo caricarlo nel nostro canale di output e aspettarci che funzioni? Scopriamolo.

### 3.1. Fornire una serie di valori come input al canale

Il buon senso suggerisce che dovremmo essere in grado di passare semplicemente una serie di valori invece di un singolo valore. Giusto?

#### 3.1.1. Imposta la variabile di input

Prendiamo la variabile `greetings_array` che abbiamo appena immaginato e rendiamola una realtà aggiungendola al blocco del workflow:

=== "Dopo"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="3 4"
    workflow {

        // declare an array of input greetings
        greetings_array = ['Hello','Bonjour','Holà']

        // create a channel for inputs
        greeting_ch = channel.of('Hello','Bonjour','Holà')
    ```

=== "Prima"

    ```groovy title="hello-channels.nf" linenums="27"
    workflow {

        // create a channel for inputs
        greeting_ch = channel.of('Hello','Bonjour','Holà')
    ```

#### 3.1.2. Imposta l'array di saluti come input alla fabbrica del canale

Sostituiremo i valori `'Hello','Bonjour', 'Holà'` attualmente codificati nella fabbrica del canale con il `greetings_array` che abbiamo appena creato.

Nel blocco del workflow, apportare la seguente modifica:

=== "Dopo"

    ```groovy title="hello-channels.nf" linenums="32" hl_lines="2"
        // create a channel for inputs
        greeting_ch = channel.of(greetings_array)
    ```

=== "Prima"

    ```groovy title="hello-channels.nf" linenums="32"
        // create a channel for inputs
        greeting_ch = channel.of('Hello','Bonjour','Holà')
    ```

#### 3.1.3. Esegui il workflow

Proviamo a eseguire questo:

```bash
nextflow run hello-channels.nf
```

Oh no! Nextflow genera un errore che inizia così:

```console title="Output" linenums="1"
 N E X T F L O W   ~  version 24.10.0

Launching `hello-channels.nf` [friendly_koch] DSL2 - revision: 97256837a7

executor >  local (1)
[22/57e015] sayHello (1) | 0 of 1
ERROR ~ Error executing process > 'sayHello (1)'

Caused by:
  Missing output file(s) `[Hello, Bonjour, Holà]-output.txt` expected by process `sayHello (1)`
```

Sembra che Nextflow abbia tentato di eseguire una singola chiamata di processo, usando `[Hello, Bonjour, Holà]` come valore di stringa, invece di utilizzare le tre stringhe nell'array come valori separati.

Come facciamo a ottenere Nextflow per decomprimere l'array e caricare le singole stringhe nel canale?

### 3.2. Usa un operatore per trasformare i contenuti del canale

È qui che entrano in gioco gli **operatori**.

Se sfogli [elenco di operatori](https://www.nextflow.io/docs/latest/reference/operator.html) nella documentazione di Nextflow, troverai [`flatten()`](https://www.nextflow.io/docs/latest/reference/operator.html#flatten), che fa esattamente ciò di cui abbiamo bisogno: decomprimere il contenuto di un array e li emette come singoli elementi.

!!! note

    È tecnicamente possibile ottenere gli stessi risultati utilizzando una fabbrica di canali diversa, [`channel.fromList`](https://nextflow.io/docs/latest/reference/channel.html#fromlist), che include una fase di mappatura implicita nella sua operazione.
    Qui abbiamo scelto di non usarlo per dimostrare l'uso di un operatore su un caso d'uso abbastanza semplice.

#### 3.2.1. Aggiungi l'operatore `flatten()`

Per applicare l'operatore `flatten()` al nostro canale di input, lo apponiamo alla dichiarazione di fabbrica del canale.

Nel blocco del workflow, apportare la seguente modifica del codice:

=== "Dopo"

    ```groovy title="hello-channels.nf" linenums="31" hl_lines="3"
        // create a channel for inputs
        greeting_ch = channel.of(greetings_array)
                             .flatten()
    ```

=== "Prima"

    ```groovy title="hello-channels.nf" linenums="31"
        // create a channel for inputs
        greeting_ch = channel.of(greetings_array)

    ```

Qui abbiamo aggiunto l'operatore sulla riga successiva per la leggibilità, ma puoi aggiungere operatori sulla stessa riga della fabbrica del canale se preferisci, come questo: `greeting_ch = channel.of(greetings_array).flatten()`

#### 3.2.2. Aggiungi `view()` per ispezionare i contenuti del canale

Potremmo eseguirlo subito per verificare se funziona, ma già che ci siamo, aggiungeremo anche un paio di operatori [`view()`](https://www.nextflow.io/docs/latest/reference/operator.html#view), che ci consentono di ispezionare il contenuto di un canale.
Puoi pensare a `view()` come a uno strumento di debug, come un'istruzione `print()` in Python, o il suo equivalente in altre lingue.

Nel blocco del workflow, apportare la seguente modifica del codice:

=== "Dopo"

    ```groovy title="hello-channels.nf" linenums="31" hl_lines="3-5"
        // create a channel for inputs
        greeting_ch = channel.of(greetings_array)
                             .view { greeting -> "Before flatten: $greeting" }
                             .flatten()
                             .view { greeting -> "After flatten: $greeting" }
    ```

=== "Prima"

    ```groovy title="hello-channels.nf" linenums="31"
        // create a channel for inputs
        greeting_ch = channel.of(greetings_array)
                             .flatten()
    ```

Stiamo usando un operatore _closure_ qui - le parentesi graffe.
Questo codice viene eseguito per ogni elemento nel canale.
Definiamo una variabile temporanea per il valore interno, qui chiamata `salunto` (potrebbe essere qualsiasi cosa).
Questa variabile viene utilizzata solo nell'ambito di tale chiusura.

In questo esempio, `$greeting` rappresenta ogni singolo elemento caricato in un canale.

!!! note "Nota su `$it`"

    In alcune pipeline potresti vedere una variabile speciale chiamata `$it` utilizzata all'interno delle chiusure dell'operatore.
    Questa è una variabile _implicita_ che consente un accesso abbreviato alla variabile interna,
    Senza bisogno di definirlo con un `->`.

    Preferiamo essere espliciti per aiutare la chiarezza del codice, in quanto tale la sintassi `$it` è scoraggiata e verrà lentamente eliminata dal linguaggio Nextflow.

#### 3.2.3. Esegui il workflow

Infine, puoi provare a eseguire di nuovo il workflow!

```bash
nextflow run hello-channels.nf
```

Questa volta funziona E ci dà una visione aggiuntiva di come appaiono i contenuti del canale prima e dopo aver eseguito l'operatore `flatten()`:

```console title="Output" linenums="1"
 N E X T F L O W   ~  version 24.10.0

Launching `hello-channels.nf` [tiny_elion] DSL2 - revision: 1d834f23d2

executor >  local (3)
[8e/bb08f3] sayHello (2) | 3 of 3 ✔
Before flatten: [Hello, Bonjour, Holà]
After flatten: Hello
After flatten: Bonjour
After flatten: Holà
```

Vedi che otteniamo una singola istruzione `Before flatten:` perché a quel punto il canale contiene un elemento, l'array originale.
Quindi otteniamo tre istruzioni separate `Dopo appiattire:`, una per ogni saluto, che ora sono singoli elementi nel canale.

È importante sottolineare che questo significa che ogni elemento può ora essere elaborato separatamente dal workflow.

!!! tip

    Dovresti eliminare o commentare le istruzioni `view()` prima di andare avanti.

    ```groovy title="hello-channels.nf" linenums="31"
    // crea un canale per gli input
    greeting_ch = channel.of(greetings_array)
                         .flatten()
    ```

    Li abbiamo lasciati nel file di soluzione `hello-channels-3.nf` a scopo di riferimento.

### Conclusioni

Sai come usare un operatore come `flatten()` per trasformare il contenuto di un canale e come usare l'operatore `view()` per ispezionare il contenuto del canale prima e dopo aver applicato un operatore.

### Cosa c'è dopo?

Scopri come fare in modo che il workflow prenda un file come fonte di valori di input.

---

## 4. Usa un operatore per analizzare i valori di input da un file CSV

Spesso accade che, quando vogliamo eseguire su più ingressi, i valori di input sono contenuti in un file.
Ad esempio, abbiamo preparato un file CSV chiamato `greetings.csv` contenente diversi saluti, uno su ogni riga (come una colonna di dati).

```csv title="greetings.csv" linenums="1"
Hello
Bonjour
Holà
```

Quindi ora dobbiamo modificare il nostro workflow per leggere i valori di un file del genere.

### 4.1. Modifica lo script per aspettarti un file CSV come fonte di saluti

Per iniziare, dovremo apportare due modifiche chiave allo script:

- Cambia il parametro di input per puntare al file CSV
- Passa a una fabbrica di canali progettata per gestire un file

#### 4.1.1. Cambia il parametro di input per puntare al file CSV

Ricordi il parametro `params.greeting` che abbiamo impostato nella parte 1?
Lo aggiorneremo per puntare al file CSV contenente i nostri saluti.

Nel blocco del workflow, apportare la seguente modifica del codice:

=== "Dopo"

    ```groovy title="hello-channels.nf" linenums="25" hl_lines="4"
    /*
     * Pipeline parameters
     */
    params.greeting = 'greetings.csv'
    ```

=== "Prima"

    ```groovy title="hello-channels.nf" linenums="25"
    /*
     * Pipeline parameters
     */
    params.greeting = ['Hello','Bonjour','Holà']
    ```

#### 4.1.2. Passa a una fabbrica di canali progettata per gestire un file

Dal momento che ora vogliamo usare un file invece di semplici stringhe come input, non possiamo usare la fabbrica di canali `channel.of()` di prima.
Dobbiamo passare all'utilizzo di una nuova fabbrica di canali, [`channel.fromPath()`](https://www.nextflow.io/docs/latest/reference/channel.html#channel-path), che ha alcune funzionalità integrate per la gestione dei percorsi dei file.

Nel blocco del flusso di lavoro, apportare la seguente modifica del codice:

=== "Dopo"

    ```groovy title="hello-channels.nf" linenums="31" hl_lines="1 2"
        // create a channel for inputs from a CSV file
        greeting_ch = channel.fromPath(params.greeting)

    ```

=== "Prima"

    ```groovy title="hello-channels.nf" linenums="31"
        // create a channel for inputs
        greeting_ch = channel.of(greetings_array)
                             .flatten()
    ```

#### 4.1.3. Esegui il workflow

Proviamo a eseguire il workflow con la nuova fabbrica di canali e il file di input.

```bash
nextflow run hello-channels.nf
```

Oh no, non funziona. Ecco l'inizio dell'output della console e il messaggio di errore:

```console title="Output" linenums="1"
 N E X T F L O W   ~  version 24.10.0

Launching `hello-channels.nf` [adoring_bhabha] DSL2 - revision: 8ce25edc39

[-        ] sayHello | 0 of 1
ERROR ~ Error executing process > 'sayHello (1)'

Caused by:
  File `/workspaces/training/hello-nextflow/data/greetings.csv-output.txt` is outside the scope of the process work directory: /workspaces/training/hello-nextflow/work/e3/c459b3c8f4029094cc778c89a4393d


Command executed:

  echo '/workspaces/training/hello-nextflow/data/greetings.csv' > '/workspaces/training/hello-nextflow/data/greetings.
```

The `Command executed:` bit (lines 13-15) is especially helpful here.

Questo può sembrare un po' familiare.
Sembra che Nextflow abbia tentato di eseguire una singola chiamata di processo utilizzando il percorso del file stesso come valore stringa.
Quindi ha risolto correttamente il percorso del file, ma in realtà non ha analizzato il suo contenuto, che è quello che volevamo.

Come facciamo a fare in modo che Nextflow apra il file e carichi i suoi contenuti nel canale?

Sembra che abbiamo bisogno di un altro [operatore](https://www.nextflow.io/docs/latest/reference/operator.html)!

### 4.2. Usa l'operatore `splitCsv()` per analizzare il file

Guardando di nuovo l'elenco degli operatori, troviamo [`splitCsv()`](https://www.nextflow.io/docs/latest/reference/operator.html#splitCsv), che è progettato per analizzare e dividere il testo in formato CSV.

#### 4.2.1. Applica `splitCsv()` al canale

Per applicare l'operatore, lo assediamo alla linea di fabbrica del canale come in precedenza.

Nel blocco del workflow, apportare la seguente modifica del codice:

=== "Dopo"

    ```groovy title="hello-channels.nf" linenums="31" hl_lines="3-5"
    // create a channel for inputs from a CSV file
    greeting_ch = channel.fromPath(params.greeting)
                         .view { csv -> "Before splitCsv: $csv" }
                         .splitCsv()
                         .view { csv -> "After splitCsv: $csv" }
    ```

=== "Prima"

    ```groovy title="hello-channels.nf" linenums="31"
    // create a channel for inputs from a CSV file
    greeting_ch = channel.fromPath(params.greeting)

    ```

Come puoi vedere, includiamo anche le dichiarazioni prima/dopo la vista mentre ci siamo.

#### 4.2.2. Esegui di nuovo il workflow

Proviamo a eseguire il workflow con la logica di parsing CSV aggiunta.

```bash
nextflow run hello-channels.nf
```

È interessante notare che anche questo fallisce, ma con un errore diverso. L'output e l'errore della console iniziano così:

```console title="Output" linenums="1"
 N E X T F L O W   ~  version 24.10.0

Launching `hello-channels.nf` [stoic_ride] DSL2 - revision: a0e5de507e

executor >  local (3)
[42/8fea64] sayHello (1) | 0 of 3
Before splitCsv: /workspaces/training/hello-nextflow/greetings.csv
After splitCsv: [Hello]
After splitCsv: [Bonjour]
After splitCsv: [Holà]
ERROR ~ Error executing process > 'sayHello (2)'

Caused by:
  Missing output file(s) `[Bonjour]-output.txt` expected by process `sayHello (2)`


Command executed:

  echo '[Bonjour]' > '[Bonjour]-output.txt'
```

Questa volta Nextflow ha analizzato il contenuto del file (eviva!) Ma sono state aggiunte parentesi intorno ai saluti.

Per farla breve, `splitCsv()` legge ogni riga in un array e ogni valore separato da virgole nella riga diventa un elemento nell'array.
Quindi qui ci dà tre array contenenti un elemento ciascuno.

!!! note

    Anche se questo comportamento sembra scomodo in questo momento, sarà estremamente utile in seguito quando ci occuperemo di file di input con più colonne di dati.

We could solve this by using `flatten()`, which you already know.
However, there's another operator called `map()` that's more appropriate to use here and is really useful to know; it pops up a lot in Nextflow pipelines.

### 4.3. Usa l'operatore `map()` per estrarre i saluti

L'operatore `map()` è un piccolo strumento molto pratico che ci consente di fare tutti i tipi di mappature ai contenuti di un canale.

In questo caso, lo useremo per estrarre quell'elemento che vogliamo da ogni riga del nostro file.
Ecco come appare la sintassi:

```groovy title="Syntax"
.map { item -> item[0] }
```

Ciò significa "per ogni elemento nel canale, prendi il primo di qualsiasi elemento che contiene".

Quindi applichiamolo al nostro parsing CSV.

#### 4.3.1. Applica `map()` al canale

Nel blocco del workflow, apportare la seguente modifica del codice:

=== "Dopo"

    ```groovy title="hello-channels.nf" linenums="31" hl_lines="6-8"
    // create a channel for inputs from a CSV file
    greeting_ch = channel.fromPath(params.greeting)
                         .view { csv -> "Before splitCsv: $csv" }
                         .splitCsv()
                         .view { csv -> "After splitCsv: $csv" }
                         .map { item -> item[0] }
                         .view { csv -> "After map: $csv" }
    ```

=== "Prima"

    ```groovy title="hello-channels.nf" linenums="31"
    // create a channel for inputs from a CSV file
    greeting_ch = channel.fromPath(params.greeting)
                         .view { csv -> "Before splitCsv: $csv" }
                         .splitCsv()
                         .view { csv -> "After splitCsv: $csv" }

    ```

Ancora una volta includiamo un'altra chiamata `view()` per confermare che l'operatore fa ciò che ci aspettiamo.

#### 4.3.2. Esegui il workflow ancora una volta

Eseguiamolo ancora una volta:

```bash
nextflow run hello-channels.nf
```

Questa volta dovrebbe funzionare senza errori.

```console title="Output" linenums="1"
 N E X T F L O W   ~  version 24.10.0

Launching `hello-channels.nf` [tiny_heisenberg] DSL2 - revision: 845b471427

executor >  local (3)
[1a/1d19ab] sayHello (2) | 3 of 3 ✔
Before splitCsv: /workspaces/training/hello-nextflow/greetings.csv
After splitCsv: [Hello]
After splitCsv: [Bonjour]
After splitCsv: [Holà]
After map: Hello
After map: Bonjour
After map: Holà
```

Guardando l'output delle istruzioni `view()`, vediamo quanto segue:

- Una singola istruzione `Before splitCsv:`: a quel punto il canale contiene un elemento, il percorso del file originale.
- Tre istruzioni separate `After splitCsv:`: una per ogni saluto, ma ognuna è contenuta all'interno di un array che corrisponde a quella riga nel file.
- Tre istruzioni separate `After map:`: una per ogni saluto, che ora sono elementi individuali nel canale.

Puoi anche guardare i file di output per verificare che ogni saluto sia stato correttamente estratto ed elaborato attraverso il workflow.

Abbiamo ottenuto lo stesso risultato di prima, ma ora abbiamo molta più flessibilità per aggiungere più elementi al canale di saluti che vogliamo elaborare modificando un file di input, senza modificare alcun codice.

!!! note

    Qui abbiamo avuto tutti i saluti su una riga nel file CSV.
    Puoi provare ad aggiungere più colonne al file CSV e vedere cosa succede; ad esempio, prova quanto segue:

    ```csv title="greetings.csv"
    Hello,English
    Bonjour,French
    Holà,Spanish
    ```

    Puoi anche provare a sostituire `.map { item -> item[0] }` con `.flatten()` e vedere cosa succede a seconda di quante righe e colonne hai nel file di input.

    Imparerai approcci più avanzati per gestire input complessi in una formazione successiva.

### Conclusioni

Sai come usare gli operatori `splitCsv()` e `map()` per leggere in un file di valori di input e gestirli in modo appropriato.

Più in generale, hai una comprensione di base di come Nextflow utilizza **canali** per gestire gli input dei processi e **operatori** per trasformare i loro contenuti.

### Cosa c'è dopo?

Fai una grande pausa, hai lavorato sodo in questo!
Quando sei pronto, passa alla parte 3 per imparare come aggiungere altri passaggi e collegarli insieme in un workflow corretto.
