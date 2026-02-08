# Parte 3: Hello Workflow

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduzione assistita da IA - [scopri di più e suggerisci miglioramenti](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/_aO56V3iXGI?si=Irl9nAQniDyICp2b&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1&amp;cc_lang_pref=it" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } Guarda [l'intera playlist](https://youtube.com/playlist?list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&si=eF7cLR62goy-lc6n) sul canale YouTube di Nextflow.

:green_book: La trascrizione del video è disponibile [qui](./transcripts/03_hello_workflow.md).
///

La maggior parte dei flussi di lavoro reali coinvolge più di un passaggio.
In questo modulo di formazione, imparerete come connettere i processi insieme in un flusso di lavoro multi-step.

Questo vi insegnerà il modo Nextflow per ottenere quanto segue:

1. Far fluire i dati da un processo al successivo
2. Raccogliere gli output da più chiamate di processo in una singola chiamata di processo
3. Passare parametri aggiuntivi a un processo
4. Gestire output multipli provenienti da un processo

Per dimostrare, continueremo a costruire sull'esempio Hello World indipendente dal dominio delle Parti 1 e 2.
Questa volta, apporteremo le seguenti modifiche al nostro flusso di lavoro per riflettere meglio come le persone costruiscono workflow reali:

1. Aggiungere un secondo passaggio che converte il saluto in maiuscolo.
2. Aggiungere un terzo passaggio che raccoglie tutti i saluti trasformati e li scrive in un singolo file.
3. Aggiungere un parametro per nominare il file di output finale e passarlo come input secondario al passaggio di raccolta.
4. Far sì che il passaggio di raccolta riporti anche una semplice statistica su ciò che è stato elaborato.

??? info "Come iniziare da questa sezione"

    Questa sezione del corso presuppone che abbiate completato le Parti 1-2 del corso [Hello Nextflow](./index.md), ma se avete familiarità con i concetti base trattati in quelle sezioni, potete iniziare da qui senza fare nulla di speciale.

---

## 0. Riscaldamento: Eseguire `hello-workflow.nf`

Useremo lo script del flusso di lavoro `hello-workflow.nf` come punto di partenza.
È equivalente allo script prodotto seguendo la Parte 2 di questo corso di formazione, tranne che abbiamo rimosso le istruzioni `view()` e cambiato la destinazione dell'output:

```groovy title="hello-workflow.nf" linenums="37" hl_lines="3"
output {
    first_output {
        path 'hello_workflow'
        mode 'copy'
    }
}
```

Questo diagramma riassume l'operazione attuale del flusso di lavoro.
Dovrebbe sembrarvi familiare, tranne che ora stiamo mostrando esplicitamente che gli output del processo sono impacchettati in un canale, proprio come lo erano gli input.
Useremo quel canale di output a breve.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-workflow-channels.svg"
</figure>

Solo per assicurarci che tutto funzioni, eseguiamo lo script una volta prima di apportare modifiche:

```bash
nextflow run hello-workflow.nf
```

??? success "Output del comando"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-workflow.nf` [admiring_lamarr] DSL2 - revision: 4d4053520d

    executor >  local (3)
    [b1/5826b5] process > sayHello (2) [100%] 3 of 3 ✔
    ```

Come in precedenza, troverete i file di output nella posizione specificata nel blocco `output`.
Per questo capitolo, è sotto `results/hello_workflow/`.

??? abstract "Contenuti della directory"

    ```console
    results/hello_workflow
    ├── Bonjour-output.txt
    ├── Hello-output.txt
    └── Holà-output.txt
    ```

Se tutto ha funzionato, siete pronti a imparare come assemblare un flusso di lavoro multi-step.

---

## 1. Aggiungere un secondo passaggio al flusso di lavoro

Aggiungeremo un passaggio per convertire ogni saluto in maiuscolo.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-multistep.svg"
</figure>

A tal fine, dobbiamo fare tre cose:

- Definire il comando che useremo per fare la conversione in maiuscolo.
- Scrivere un nuovo processo che incapsula il comando di conversione in maiuscolo.
- Chiamare il nuovo processo nel blocco workflow e configurarlo per prendere l'output del processo `sayHello()` come input.

### 1.1. Definire il comando di conversione in maiuscolo e testarlo nel terminale

Per effettuare la conversione dei saluti in maiuscolo, useremo un classico strumento UNIX chiamato `tr` per 'text replacement' (sostituzione di testo), con la seguente sintassi:

```bash title="Sintassi"
tr '[a-z]' '[A-Z]'
```

Questa è una sostituzione di testo molto semplice che non tiene conto delle lettere accentate, quindi per esempio 'Holà' diventerà 'HOLà', ma farà un lavoro abbastanza buono per dimostrare i concetti di Nextflow e questo è ciò che conta.

Per testarlo, possiamo eseguire il comando `echo 'Hello World'` e fare il pipe del suo output al comando `tr`:

```bash
echo 'Hello World' | tr '[a-z]' '[A-Z]' > UPPER-output.txt
```

L'output è un file di testo chiamato `UPPER-output.txt` che contiene la versione maiuscola della stringa `Hello World`.

??? abstract "Contenuti del file"

    ```console title="UPPER-output.txt"
    HELLO WORLD
    ```

Questo è fondamentalmente ciò che cercheremo di fare con il nostro flusso di lavoro.

### 1.2. Scrivere il passaggio di conversione in maiuscolo come processo Nextflow

Possiamo modellare il nostro nuovo processo sul primo, dato che vogliamo usare tutti gli stessi componenti.

Aggiungete la seguente definizione di processo allo script del flusso di lavoro, subito sotto il primo:

```groovy title="hello-workflow.nf" linenums="20"
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

In questo, componiamo il nome del secondo file di output basandoci sul nome del file di input, similmente a quello che abbiamo fatto originariamente per l'output del primo processo.

### 1.3. Aggiungere una chiamata al nuovo processo nel blocco workflow

Ora dobbiamo dire a Nextflow di chiamare effettivamente il processo che abbiamo appena definito.

Nel blocco workflow, effettuate la seguente modifica al codice:

=== "Dopo"

    ```groovy title="hello-workflow.nf" linenums="44" hl_lines="10-11"
    workflow {

        main:
        // crea un canale per gli input da un file CSV
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        // emette un saluto
        sayHello(greeting_ch)
        // converte il saluto in maiuscolo
        convertToUpper()

        publish:
        first_output = sayHello.out
    }
    ```

=== "Prima"

    ```groovy title="hello-workflow.nf" linenums="44"
    workflow {

        main:
        // crea un canale per gli input da un file CSV
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        // emette un saluto
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

Questo non è ancora funzionale perché non abbiamo specificato cosa dovrebbe essere l'input al processo `convertToUpper()`.

### 1.4. Passare l'output del primo processo al secondo processo

Ora dobbiamo far fluire l'output del processo `sayHello()` nel processo `convertToUpper()`.

Convenientemente, Nextflow impacchetta automaticamente l'output di un processo in un canale, come mostrato nel diagramma nella sezione di riscaldamento.
Possiamo riferirci al canale di output di un processo come `<process>.out`.

Quindi l'output del processo `sayHello` è un canale chiamato `sayHello.out`, che possiamo collegare direttamente alla chiamata di `convertToUpper()`.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-multistep-connector.svg"
</figure>

Nel blocco workflow, effettuate la seguente modifica al codice:

=== "Dopo"

    ```groovy title="hello-workflow.nf" linenums="53" hl_lines="2"
        // converte il saluto in maiuscolo
        convertToUpper(sayHello.out)
    ```

=== "Prima"

    ```groovy title="hello-workflow.nf" linenums="53" hl_lines="2"
        // converte il saluto in maiuscolo
        convertToUpper()
    ```

Per un caso semplice come questo (un output verso un input), questo è tutto ciò che dobbiamo fare per connettere due processi!

### 1.5. Configurare la pubblicazione dell'output del flusso di lavoro

Infine, aggiorniamo gli output del flusso di lavoro per pubblicare anche i risultati del secondo processo.

#### 1.5.1. Aggiornare la sezione `publish:` del blocco `workflow`

Nel blocco `workflow`, effettuate la seguente modifica al codice:

=== "Dopo"

    ```groovy title="hello-workflow.nf" linenums="56" hl_lines="3"
        publish:
        first_output = sayHello.out
        uppercased = convertToUpper.out
    }
    ```

=== "Prima"

    ```groovy title="hello-workflow.nf" linenums="56"
        publish:
        first_output = sayHello.out
    }
    ```

La logica è la stessa di prima.

#### 1.5.2. Aggiornare il blocco `output`

Nel blocco `output`, effettuate la seguente modifica al codice:

=== "Dopo"

    ```groovy title="hello-workflow.nf" linenums="61" hl_lines="6-9"
    output {
        first_output {
            path 'hello_workflow'
            mode 'copy'
        }
        uppercased {
            path 'hello_workflow'
            mode 'copy'
        }
    }
    ```

=== "Prima"

    ```groovy title="hello-workflow.nf" linenums="61"
    output {
        first_output {
            path 'hello_workflow'
            mode 'copy'
        }
    }
    ```

Ancora una volta, la logica è la stessa di prima.

Questo vi mostra che potete controllare le impostazioni di output a un livello molto granulare, per ogni singolo output.
Sentitevi liberi di provare a cambiare i percorsi o la modalità di pubblicazione per uno dei processi per vedere cosa succede.

Naturalmente, questo significa che stiamo ripetendo alcune informazioni qui, il che potrebbe diventare scomodo se volessimo aggiornare la posizione per tutti gli output nello stesso modo.
Più avanti nel corso, imparerete come configurare queste impostazioni per output multipli in modo strutturato.

### 1.6. Eseguire il flusso di lavoro con `-resume`

Testiamo questo usando il flag `-resume`, dato che abbiamo già eseguito con successo il primo passaggio del flusso di lavoro.

```bash
nextflow run hello-workflow.nf -resume
```

??? success "Output del comando"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-workflow.nf` [high_cantor] DSL2 - revision: d746983511

    executor >  local (3)
    [ab/816321] process > sayHello (3)       [100%] 3 of 3, cached: 3 ✔
    [e0/ecf81b] process > convertToUpper (3) [100%] 3 of 3 ✔
    ```

Ora c'è una riga extra nell'output della console che corrisponde al nuovo processo che abbiamo appena aggiunto.

Troverete gli output nella directory `results/hello_workflow` come impostato nel blocco `output`.

??? abstract "Contenuti della directory"

    ```console
    results/hello_workflow/
    ├── Bonjour-output.txt
    ├── Hello-output.txt
    ├── Holà-output.txt
    ├── UPPER-Bonjour-output.txt
    ├── UPPER-Hello-output.txt
    └── UPPER-Holà-output.txt
    ```

Comodo! Ma vale comunque la pena dare un'occhiata dentro la directory di lavoro di una delle chiamate al secondo processo.

??? abstract "Contenuti della directory"

    ```console
    work/e0/ecf81b4cacc648b9b994218d5b29d7/
    ├── Holà-output.txt -> /workspaces/training/hello-nextflow/work/ab/81632178cd37e9e815959278808819/Holà-output.txt
    └── UPPER-Holà-output.txt
    ```

Notate che ci sono due file `*-output`: l'output del primo processo così come l'output del secondo.

L'output del primo processo è lì perché Nextflow lo ha **staged** lì per avere tutto il necessario per l'esecuzione nella stessa sottodirectory.

Tuttavia, è in realtà un link simbolico che punta al file originale nella sottodirectory della prima chiamata di processo.
Per impostazione predefinita, quando si esegue su una singola macchina come stiamo facendo qui, Nextflow usa link simbolici piuttosto che copie per fare lo staging dei file di input e intermedi.

Ora, prima di procedere, pensate a come tutto ciò che abbiamo fatto è stato connettere l'output di `sayHello` all'input di `convertToUpper` e i due processi hanno potuto essere eseguiti in serie.
Nextflow ha fatto il lavoro duro di gestire i singoli file di input e output e passarli tra i due comandi per noi.

Questa è una delle ragioni per cui i canali di Nextflow sono così potenti: si occupano del lavoro noioso coinvolto nel connettere insieme i passaggi del flusso di lavoro.

### Takeaway

Sapete come concatenare processi insieme fornendo l'output di un passaggio come input al passaggio successivo.

### Cosa c'è dopo?

Imparare come raccogliere gli output da chiamate di processo in batch e alimentarli in un singolo processo.

---

## 2. Aggiungere un terzo passaggio per raccogliere tutti i saluti

Quando usiamo un processo per applicare una trasformazione a ciascuno degli elementi in un canale, come stiamo facendo qui ai saluti multipli, a volte vogliamo raccogliere elementi dal canale di output di quel processo, e alimentarli in un altro processo che esegue qualche tipo di analisi o somma.

Per dimostrare, aggiungeremo un nuovo passaggio alla nostra pipeline che raccoglie tutti i saluti maiuscoli prodotti dal processo `convertToUpper` e li scrive in un singolo file.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-collect.svg"
</figure>

Senza rovinare la sorpresa, questo coinvolgerà un operatore molto utile.

### 2.1. Definire il comando di raccolta e testarlo nel terminale

Il passaggio di raccolta che vogliamo aggiungere al nostro flusso di lavoro userà il comando `cat` per concatenare più saluti maiuscoli in un singolo file.

Eseguiamo il comando da solo nel terminale per verificare che funzioni come previsto, proprio come abbiamo fatto in precedenza.

Eseguite il seguente nel vostro terminale:

```bash
echo 'Hello' | tr '[a-z]' '[A-Z]' > UPPER-Hello-output.txt
echo 'Bonjour' | tr '[a-z]' '[A-Z]' > UPPER-Bonjour-output.txt
echo 'Holà' | tr '[a-z]' '[A-Z]' > UPPER-Holà-output.txt
cat UPPER-Hello-output.txt UPPER-Bonjour-output.txt UPPER-Holà-output.txt > COLLECTED-output.txt
```

L'output è un file di testo chiamato `COLLECTED-output.txt` che contiene le versioni maiuscole dei saluti originali.

??? abstract "Contenuti del file"

    ```console title="COLLECTED-output.txt"
    HELLO
    BONJOUR
    HOLà
    ```

Questo è il risultato che vogliamo ottenere con il nostro flusso di lavoro.

### 2.2. Creare un nuovo processo per fare il passaggio di raccolta

Creiamo un nuovo processo e chiamiamolo `collectGreetings()`.
Possiamo iniziare a scriverlo basandoci su ciò che abbiamo visto prima.

#### 2.2.1. Scrivere le parti 'ovvie' del processo

Aggiungete la seguente definizione di processo allo script del flusso di lavoro:

```groovy title="hello-workflow.nf" linenums="37"
/*
 * Collect uppercase greetings into a single output file
 */
process collectGreetings {

    input:
    ???

    output:
    path "COLLECTED-output.txt"

    script:
    """
    cat ??? > 'COLLECTED-output.txt'
    """
}
```

Questo è ciò che possiamo scrivere con fiducia basandoci su ciò che avete imparato finora.
Ma questo non è funzionale!
Omette la/le definizione/i di input e la prima metà del comando script perché dobbiamo capire come scriverlo.

#### 2.2.2. Definire gli input di `collectGreetings()`

Dobbiamo raccogliere i saluti da tutte le chiamate al processo `convertToUpper()`.
Cosa sappiamo di poter ottenere dal passaggio precedente nel flusso di lavoro?

Il canale emesso da `convertToUpper()` conterrà i percorsi ai singoli file contenenti i saluti maiuscoli.
Questo equivale a uno slot di input; chiamiamolo `input_files` per semplicità.

Nel blocco process, effettuate la seguente modifica al codice:

=== "Dopo"

    ```groovy title="hello-workflow.nf" linenums="42" hl_lines="2"
          input:
          path input_files
    ```

=== "Prima"

    ```groovy title="hello-workflow.nf" linenums="42" hl_lines="2"
          input:
          ???
    ```

Notate che usiamo il prefisso `path` anche se ci aspettiamo che questo contenga più file.

#### 2.2.3. Comporre il comando di concatenazione

Qui le cose potrebbero diventare un po' complicate, perché dobbiamo essere in grado di gestire un numero arbitrario di file di input.
Nello specifico, non possiamo scrivere il comando in anticipo, quindi dobbiamo dire a Nextflow come comporlo a runtime basandosi su quali input fluiscono nel processo.

In altre parole, se abbiamo un canale di input contenente l'elemento `[file1.txt, file2.txt, file3.txt]`, abbiamo bisogno che Nextflow lo trasformi in `cat file1.txt file2.txt file3.txt`.

Fortunatamente, Nextflow è perfettamente in grado di farlo per noi se scriviamo semplicemente `cat ${input_files}` nel comando script.

Nel blocco process, effettuate la seguente modifica al codice:

=== "Dopo"

    ```groovy title="hello-workflow.nf" linenums="54" hl_lines="3"
        script:
        """
        cat ${input_files} > 'COLLECTED-output.txt'
        """
    ```

=== "Prima"

    ```groovy title="hello-workflow.nf" linenums="54"
        script:
        """
        cat ??? > 'COLLECTED-output.txt'
        """
    ```

In teoria questo dovrebbe gestire qualsiasi numero arbitrario di file di input.

!!! tip "Suggerimento"

    Alcuni strumenti da riga di comando richiedono di fornire un argomento (come `-input`) per ogni file di input.
    In quel caso, dovremmo fare un po' di lavoro extra per comporre il comando.
    Potete vedere un esempio di questo nel corso di formazione [Nextflow for Genomics](../../nf4_science/genomics/).

### 2.3. Aggiungere il passaggio di raccolta al flusso di lavoro

Ora dovremmo solo aver bisogno di chiamare il processo di raccolta sull'output del passaggio di conversione in maiuscolo.
Quello è anche un canale, chiamato `convertToUpper.out`.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-collect-connector.svg"
</figure>

#### 2.3.1. Connettere le chiamate dei processi

Nel blocco workflow, effettuate la seguente modifica al codice:

=== "Dopo"

    ```groovy title="hello-workflow.nf" linenums="75" hl_lines="4 5"
        // converte il saluto in maiuscolo
        convertToUpper(sayHello.out)

        // raccoglie tutti i saluti in un file
        collectGreetings(convertToUpper.out)
    }
    ```

=== "Prima"

    ```groovy title="hello-workflow.nf" linenums="75"
        // converte il saluto in maiuscolo
        convertToUpper(sayHello.out)
    }
    ```

Questo connette l'output di `convertToUpper()` all'input di `collectGreetings()`.

#### 2.3.2. Eseguire il flusso di lavoro con `-resume`

Proviamo.

```bash
nextflow run hello-workflow.nf -resume
```

??? success "Output del comando"

    ```console hl_lines="8"
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-workflow.nf` [mad_gilbert] DSL2 - revision: 6acfd5e28d

    executor >  local (3)
    [79/33b2f0] sayHello (2)         | 3 of 3, cached: 3 ✔
    [99/79394f] convertToUpper (3)   | 3 of 3, cached: 3 ✔
    [47/50fe4a] collectGreetings (1) | 3 of 3 ✔
    ```

Viene eseguito con successo, incluso il terzo passaggio.

Tuttavia, guardate il numero di chiamate per `collectGreetings()` nell'ultima riga.
Ce ne aspettavamo solo una, ma ce ne sono tre.

Ora date un'occhiata ai contenuti del file di output finale.

??? abstract "Contenuti del file"

    ```console title="results/COLLECTED-output.txt"
    Holà
    ```

Oh no. Il passaggio di raccolta è stato eseguito individualmente su ogni saluto, il che NON è quello che volevamo.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-collect-no-operator.svg"
</figure>

Dobbiamo fare qualcosa per dire esplicitamente a Nextflow che vogliamo che quel terzo passaggio venga eseguito su tutti gli elementi nel canale emesso da `convertToUpper()`.

### 2.4. Usare un operatore per raccogliere i saluti in un singolo input

Sì, ancora una volta la risposta al nostro problema è un operatore.

Nello specifico, useremo l'operatore [`collect()`](https://nextflow.io/docs/latest/reference/operator.html#collect) dal nome appropriato.

#### 2.4.1. Aggiungere l'operatore `collect()`

Questa volta apparirà un po' diverso perché non stiamo aggiungendo l'operatore nel contesto di una channel factory; lo stiamo aggiungendo a un canale di output.

Prendiamo il `convertToUpper.out` e aggiungiamo l'operatore `collect()`, che ci dà `convertToUpper.out.collect()`.
Possiamo collegarlo direttamente alla chiamata del processo `collectGreetings()`.

Nel blocco workflow, effettuate la seguente modifica al codice:

=== "Dopo"

    ```groovy title="hello-workflow.nf" linenums="73" hl_lines="2"
        // raccoglie tutti i saluti in un file
        collectGreetings(convertToUpper.out.collect())
    }
    ```

=== "Prima"

    ```groovy title="hello-workflow.nf" linenums="73" hl_lines="2"
        // raccoglie tutti i saluti in un file
        collectGreetings(convertToUpper.out)
    }
    ```

#### 2.4.2. Aggiungere alcune istruzioni `view()`

Includiamo anche un paio di istruzioni `view()` per visualizzare gli stati prima e dopo dei contenuti del canale.

=== "Dopo"

    ```groovy title="hello-workflow.nf" linenums="73" hl_lines="4-6"
        // raccoglie tutti i saluti in un file
        collectGreetings(convertToUpper.out.collect())

        // istruzioni view opzionali
        convertToUpper.out.view { contents -> "Prima di collect: $contents" }
        convertToUpper.out.collect().view { contents -> "Dopo collect: $contents" }
    }
    ```

=== "Prima"

    ```groovy title="hello-workflow.nf" linenums="73"
        // raccoglie tutti i saluti in un file
        collectGreetings(convertToUpper.out.collect())
    }
    ```

Le istruzioni `view()` possono andare dove volete; le abbiamo messe subito dopo la chiamata per leggibilità.

#### 2.4.3. Eseguire di nuovo il flusso di lavoro con `-resume`

Proviamo:

```bash
nextflow run hello-workflow.nf -resume
```

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-workflow.nf` [soggy_franklin] DSL2 - revision: bc8e1b2726

    [d6/cdf466] sayHello (1)       | 3 of 3, cached: 3 ✔
    [99/79394f] convertToUpper (2) | 3 of 3, cached: 3 ✔
    [1e/83586c] collectGreetings   | 1 of 1 ✔
    Prima di collect: /workspaces/training/hello-nextflow/work/b3/d52708edba8b864024589285cb3445/UPPER-Bonjour-output.txt
    Prima di collect: /workspaces/training/hello-nextflow/work/99/79394f549e3040dfc2440f69ede1fc/UPPER-Hello-output.txt
    Prima di collect: /workspaces/training/hello-nextflow/work/aa/56bfe7cf00239dc5badc1d04b60ac4/UPPER-Holà-output.txt
    Dopo collect: [/workspaces/training/hello-nextflow/work/b3/d52708edba8b864024589285cb3445/UPPER-Bonjour-output.txt, /workspaces/training/hello-nextflow/work/99/79394f549e3040dfc2440f69ede1fc/UPPER-Hello-output.txt, /workspaces/training/hello-nextflow/work/aa/56bfe7cf00239dc5badc1d04b60ac4/UPPER-Holà-output.txt]
    ```

Viene eseguito con successo, anche se l'output del log potrebbe apparire un po' più disordinato di questo (lo abbiamo ripulito per leggibilità).

Questa volta il terzo passaggio è stato chiamato solo una volta!
Guardando l'output delle istruzioni `view()`, vediamo quanto segue:

- Tre istruzioni `Prima di collect:`, una per ogni saluto: a quel punto i percorsi dei file sono elementi individuali nel canale.
- Una singola istruzione `Dopo collect:`: i tre percorsi dei file sono ora impacchettati in un singolo elemento.

Possiamo riassumere questo con il seguente diagramma:

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-collect-WITH-operator.svg"
</figure>

Infine, potete dare un'occhiata ai contenuti del file di output per assicurarvi che tutto abbia funzionato correttamente.

??? abstract "Contenuti del file"

    ```console title="results/COLLECTED-output.txt"
    BONJOUR
    HELLO
    HOLà
    ```

Questa volta abbiamo tutti e tre i saluti nel file di output finale. Successo!

!!! note "Nota"

    Se eseguite questo più volte senza `-resume`, vedrete che l'ordine dei saluti cambia da un'esecuzione all'altra.
    Questo vi mostra che l'ordine in cui gli elementi fluiscono attraverso le chiamate dei processi non è garantito essere consistente.

#### 2.4.4. Rimuovere le istruzioni `view()` per leggibilità

Prima di passare alla prossima sezione, raccomandiamo di cancellare le istruzioni `view()` per evitare di ingombrare l'output della console.

=== "Dopo"

    ```groovy title="hello-workflow.nf" linenums="73"
        // raccoglie tutti i saluti in un file
        collectGreetings(convertToUpper.out.collect())
    ```

=== "Prima"

    ```groovy title="hello-workflow.nf" linenums="73" hl_lines="4-6"
        // raccoglie tutti i saluti in un file
        collectGreetings(convertToUpper.out.collect())

        // istruzioni view opzionali
        convertToUpper.out.view { contents -> "Prima di collect: $contents" }
        convertToUpper.out.collect().view { contents -> "Dopo collect: $contents" }
    ```

Questa è fondamentalmente l'operazione inversa dal punto 2.4.2.

### Takeaway

Sapete come raccogliere gli output da un batch di chiamate di processo e alimentarli in un passaggio di analisi o somma congiunta.

Per ricapitolare, questo è ciò che avete costruito finora:

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-collect.svg"
</figure>

### Cosa c'è dopo?

Imparare come passare più di un input a un processo.

---

## 3. Passare parametri aggiuntivi a un processo

Vogliamo essere in grado di nominare il file di output finale con qualcosa di specifico per poter elaborare batch successivi di saluti senza sovrascrivere i risultati finali.

A tal fine, apporteremo i seguenti perfezionamenti al flusso di lavoro:

- Modificare il processo di raccolta per accettare un nome definito dall'utente per il file di output (`batch_name`)
- Aggiungere un parametro da riga di comando al flusso di lavoro (`--batch`) e passarlo al processo di raccolta

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-collect-batch.svg"
</figure>

### 3.1. Modificare il processo di raccolta

Dovremo dichiarare l'input aggiuntivo e integrarlo nel nome del file di output.

#### 3.1.1. Dichiarare l'input aggiuntivo

Buone notizie: possiamo dichiarare quante variabili di input vogliamo nella definizione del processo.
Chiamiamo questa `batch_name`.

Nel blocco process, effettuate la seguente modifica al codice:

=== "Dopo"

    ```groovy title="hello-workflow.nf" linenums="42" hl_lines="3"
        input:
        path input_files
        val batch_name
    ```

=== "Prima"

    ```groovy title="hello-workflow.nf" linenums="42"
        input:
        path input_files
    ```

Potete configurare i vostri processi per aspettarsi quanti input volete.
Al momento, questi sono tutti configurati come input obbligatori; _dovete_ fornire un valore affinché il flusso di lavoro funzioni.

Imparerete come gestire input obbligatori vs. opzionali più avanti nel vostro percorso con Nextflow.

#### 3.1.2. Usare la variabile `batch_name` nel nome del file di output

Possiamo inserire la variabile nel nome del file di output nello stesso modo in cui abbiamo composto nomi di file dinamici prima.

Nel blocco process, effettuate la seguente modifica al codice:

=== "Dopo"

    ```groovy title="hello-workflow.nf" linenums="46" hl_lines="2 6"
        output:
        path "COLLECTED-${batch_name}-output.txt"

        script:
        """
        cat ${input_files} > 'COLLECTED-${batch_name}-output.txt'
        """
    ```

=== "Prima"

    ```groovy title="hello-workflow.nf" linenums="46" hl_lines="2 6"
        output:
        path "COLLECTED-output.txt"

        script:
        """
        cat ${input_files} > 'COLLECTED-output.txt'
        """
    ```

Questo configura il processo per usare il valore `batch_name` per generare un nome di file specifico per l'output finale del flusso di lavoro.

### 3.2. Aggiungere un parametro da riga di comando `batch`

Ora abbiamo bisogno di un modo per fornire il valore per `batch_name` e alimentarlo nella chiamata del processo.

#### 3.2.1. Usare `params` per configurare il parametro

Sapete già come usare il sistema `params` per dichiarare parametri CLI.
Usiamolo per dichiarare un parametro `batch` (con un valore predefinito perché siamo pigri).

Nella sezione dei parametri della pipeline, effettuate le seguenti modifiche al codice:

=== "Dopo"

    ```groovy title="hello-workflow.nf" linenums="55" hl_lines="6"
    /*
     * Pipeline parameters
     */
    params {
        input: Path = 'data/greetings.csv'
        batch: String = 'batch'
    }
    ```

=== "Prima"

    ```groovy title="hello-workflow.nf" linenums="55"
    /*
     * Pipeline parameters
     */
    params {
        input: Path = 'data/greetings.csv'
    }
    ```

Proprio come abbiamo dimostrato per `--input`, potete sovrascrivere quel valore predefinito specificando un valore con `--batch` sulla riga di comando.

#### 3.2.2. Passare il parametro `batch` al processo

Per fornire il valore del parametro al processo, dobbiamo aggiungerlo nella chiamata del processo.

Nel blocco workflow, effettuate la seguente modifica al codice:

=== "Dopo"

    ```groovy title="hello-workflow.nf" linenums="74" hl_lines="2"
        // raccoglie tutti i saluti in un file
        collectGreetings(convertToUpper.out.collect(), params.batch)
    ```

=== "Prima"

    ```groovy title="hello-workflow.nf" linenums="74" hl_lines="2"
        // raccoglie tutti i saluti in un file
        collectGreetings(convertToUpper.out.collect())
    ```

Vedete che per fornire input multipli a un processo, li si elenca semplicemente nelle parentesi della chiamata, separati da virgole.

!!! warning "Avviso"

    DOVETE fornire gli input al processo ESATTAMENTE NELLO STESSO ORDINE in cui sono elencati nel blocco di definizione degli input del processo.

### 3.3. Eseguire il flusso di lavoro

Proviamo a eseguirlo con un nome di batch sulla riga di comando.

```bash
nextflow run hello-workflow.nf -resume --batch trio
```

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-workflow.nf` [confident_rutherford] DSL2 - revision: bc58af409c

    executor >  local (1)
    [79/33b2f0] sayHello (2)       | 3 of 3, cached: 3 ✔
    [99/79394f] convertToUpper (2) | 3 of 3, cached: 3 ✔
    [b5/f19efe] collectGreetings   | 1 of 1 ✔
    ```

Viene eseguito con successo e produce l'output desiderato:

??? abstract "Contenuti del file"

    ```console title="results/COLLECTED-trio-output.txt"
    HELLO
    BONJOUR
    HOLà
    ```

Ora, purché specifichiamo il parametro appropriatamente, le esecuzioni successive su altri batch di input non sovrascriveranno i risultati precedenti.

### Takeaway

Sapete come passare più di un input a un processo.

### Cosa c'è dopo?

Imparare come emettere output multipli e gestirli comodamente.

---

## 4. Aggiungere un output al passaggio di raccolta

Finora abbiamo usato processi che producevano solo un output ciascuno.
Abbiamo potuto accedere ai rispettivi output molto comodamente usando la sintassi `<process>.out`, che abbiamo usato sia nel contesto di passare un output al processo successivo (es. `convertToUpper(sayHello.out)`) che nel contesto della sezione `publish:` (es. `first_output = sayHello.out`).

Cosa succede quando un processo produce più di uno?
Come gestiamo gli output multipli?
Possiamo selezionare e usare un output specifico?

Tutte ottime domande, e la risposta breve è sì, possiamo!

Output multipli saranno impacchettati in canali separati.
Possiamo scegliere di dare nomi a quei canali di output, il che rende facile riferirsi a loro individualmente in seguito, oppure possiamo riferirci a loro per indice.

A scopo dimostrativo, diciamo che vogliamo contare il numero di saluti che vengono raccolti per un dato batch di input e riportarlo in un file.

### 4.1. Modificare il processo per contare e produrre il numero di saluti

Questo richiederà due modifiche chiave alla definizione del processo: abbiamo bisogno di un modo per contare i saluti e scrivere un file di report, poi dobbiamo aggiungere quel file di report al blocco `output` del processo.

#### 4.1.1. Contare il numero di saluti raccolti

Convenientemente, Nextflow ci permette di aggiungere codice arbitrario nel blocco `script:` della definizione del processo, il che è davvero utile per fare cose come questa.

Questo significa che possiamo usare la funzione incorporata `size()` di Nextflow per ottenere il numero di file nell'array `input_files`, e scrivere il risultato su file con un comando `echo`.

Nel blocco process `collectGreetings`, effettuate le seguenti modifiche al codice:

=== "Dopo"

    ```groovy title="hello-workflow.nf" linenums="55" hl_lines="2 5"
        script:
        count_greetings = input_files.size()
        """
        cat ${input_files} > 'COLLECTED-${batch_name}-output.txt'
        echo 'There were ${count_greetings} greetings in this batch.' > '${batch_name}-report.txt'
        """
    ```

=== "Prima"

    ```groovy title="hello-workflow.nf" linenums="55"
        script:
        """
        cat ${input_files} > 'COLLECTED-${batch_name}-output.txt'
        """
    ```

La variabile `count_greetings` sarà calcolata a runtime.

#### 4.1.2. Emettere il file di report e nominare gli output

In principio tutto ciò che dobbiamo fare è aggiungere il file di report al blocco `output:`.

Tuttavia, già che ci siamo, aggiungeremo anche alcuni tag `emit:` alle nostre dichiarazioni di output. Questi ci permetteranno di selezionare gli output per nome invece di dover usare indici posizionali.

Nel blocco process, effettuate la seguente modifica al codice:

=== "Dopo"

    ```groovy title="hello-workflow.nf" linenums="46" hl_lines="2 3"
        output:
        path "COLLECTED-${batch_name}-output.txt", emit: outfile
        path "${batch_name}-report.txt", emit: report
    ```

=== "Prima"

    ```groovy title="hello-workflow.nf" linenums="46"
        output:
        path "COLLECTED-${batch_name}-output.txt"
    ```

I tag `emit:` sono opzionali, e avremmo potuto aggiungere un tag solo a uno degli output.
Ma come dice il detto, perché non entrambi?

!!! tip "Suggerimento"

    Se non nominate gli output di un processo usando `emit:`, potete comunque accedervi individualmente usando il rispettivo indice (a base zero).
    Per esempio, usereste `<process>.out[0]` per ottenere il primo output, `<process>.out[1]` per ottenere il secondo output, e così via.

    Preferiamo nominare gli output perché altrimenti, è troppo facile prendere l'indice sbagliato per errore, specialmente quando il processo produce molti output.

### 4.2. Aggiornare gli output del flusso di lavoro

Ora che abbiamo due output provenienti dal processo `collectGreetings`, l'output `collectGreetings.out` contiene due canali:

- `collectGreetings.out.outfile` contiene il file di output finale
- `collectGreetings.out.report` contiene il file di report

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-collect-report.svg"
</figure>

Dobbiamo aggiornare gli output del flusso di lavoro di conseguenza.

#### 4.2.1. Aggiornare la sezione `publish:`

Nel `blocco workflow`, effettuate la seguente modifica al codice:

=== "Dopo"

    ```groovy title="hello-workflow.nf" linenums="80" hl_lines="4 5"
        publish:
        first_output = sayHello.out
        uppercased = convertToUpper.out
        collected = collectGreetings.out.outfile
        batch_report = collectGreetings.out.report
    ```

=== "Prima"

    ```groovy title="hello-workflow.nf" linenums="80" hl_lines="4"
        publish:
        first_output = sayHello.out
        uppercased = convertToUpper.out
        collected = collectGreetings.out
    ```

Come potete vedere, riferirsi a output specifici dei processi è ora banale.
Quando aggiungeremo un altro passaggio alla nostra pipeline nella Parte 5 (Containers), saremo in grado di riferirci facilmente a `collectGreetings.out.outfile` e passarlo al nuovo processo (spoiler: il nuovo processo si chiama `cowpy`).

Ma per ora, finiamo di aggiornare gli output a livello di flusso di lavoro.

#### 4.2.2. Aggiornare il blocco `output`

Nel blocco `output`, effettuate la seguente modifica al codice:

=== "Dopo"

    ```groovy title="hello-workflow.nf" linenums="86" hl_lines="14-17"
    output {
        first_output {
            path 'hello_workflow'
            mode 'copy'
        }
        uppercased {
            path 'hello_workflow'
            mode 'copy'
        }
        collected {
            path 'hello_workflow'
            mode 'copy'
        }
        batch_report {
            path 'hello_workflow'
            mode 'copy'
        }
    }
    ```

=== "Prima"

    ```groovy title="hello-workflow.nf" linenums="80"
    output {
        first_output {
            path 'hello_workflow'
            mode 'copy'
        }
        uppercased {
            path 'hello_workflow'
            mode 'copy'
        }
        collected {
            path 'hello_workflow'
            mode 'copy'
        }
    }
    ```

Non abbiamo bisogno di aggiornare la definizione dell'output `collected` dato che quel nome non è cambiato.
Dobbiamo solo aggiungere il nuovo output.

### 4.3. Eseguire il flusso di lavoro

Proviamo a eseguirlo con l'attuale batch di saluti.

```bash
nextflow run hello-workflow.nf -resume --batch trio
```

??? success "Output del comando"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-workflow.nf` [ecstatic_wilson] DSL2 - revision: c80285f8c8

    executor >  local (1)
    [c5/4c6ca9] sayHello (3)       [100%] 3 of 3, cached: 3 ✔
    [0e/6cbc59] convertToUpper (3) [100%] 3 of 3, cached: 3 ✔
    [02/61ead2] collectGreetings   [100%] 1 of 1 ✔
    ```

Se guardate nella directory `results/hello_workflow/`, troverete il nuovo file di report, `trio-report.txt`.
Apritelo per verificare che il flusso di lavoro abbia correttamente riportato il conteggio dei saluti che sono stati elaborati.

??? abstract "Contenuti del file"

    ```txt title="trio-report.txt"
    There were 3 greetings in this batch.
    ```

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-collect-4-way.svg"
</figure>

Sentitevi liberi di aggiungere più saluti al CSV e testare cosa succede.

### Takeaway

Sapete come far emettere a un processo output multipli nominati e come gestirli appropriatamente a livello di flusso di lavoro.

Più in generale, capite i principi chiave coinvolti nel connettere processi insieme in modi comuni.

### Cosa c'è dopo?

Prendetevi una pausa extra lunga, ve la siete guadagnata.

Quando siete pronti, passate alla [**Parte 4: Hello Modules**](./04_hello_modules.md) per imparare come modularizzare il vostro codice per una migliore manutenibilità ed efficienza del codice.

---

## Quiz

<quiz>
Come si accede all'output di un processo nel blocco workflow?
- [ ] `process.output`
- [ ] `output.processName`
- [x] `processName.out`
- [ ] `get(processName)`

Per approfondire: [1.4. Passare l'output del primo processo al secondo processo](#14-passare-loutput-del-primo-processo-al-secondo-processo)
</quiz>

<quiz>
Cosa determina l'ordine di esecuzione dei processi in Nextflow?
- [ ] L'ordine in cui i processi sono scritti nel blocco workflow
- [ ] Ordine alfabetico per nome del processo
- [x] Le dipendenze dei dati tra i processi
- [ ] Ordine casuale per esecuzione parallela

Per approfondire: [1.4. Passare l'output del primo processo al secondo processo](#14-passare-loutput-del-primo-processo-al-secondo-processo)
</quiz>

<quiz>
Quale operatore dovrebbe sostituire `???` per raccogliere tutti gli output in una singola lista per il processo a valle?

```groovy hl_lines="4"
workflow {
    greetings_ch = Channel.of('Hello', 'Bonjour', 'Hola')
    SAYHELLO(greetings_ch)
    GATHER_ALL(SAYHELLO.out.???)
}
```

- [ ] `flatten()`
- [x] `collect()`
- [ ] `mix()`
- [ ] `join()`

Per approfondire: [2.4. Usare un operatore per raccogliere i saluti in un singolo input](#24-usare-un-operatore-per-raccogliere-i-saluti-in-un-singolo-input)
</quiz>

<quiz>
Quando si dovrebbe usare l'operatore `collect()`?
- [ ] Quando si vogliono elaborare gli elementi in parallelo
- [ ] Quando si ha bisogno di filtrare i contenuti del canale
- [x] Quando un processo a valle ha bisogno di tutti gli elementi da un processo a monte
- [ ] Quando si vogliono dividere i dati tra più processi

Per approfondire: [2.4. Usare un operatore per raccogliere i saluti in un singolo input](#24-usare-un-operatore-per-raccogliere-i-saluti-in-un-singolo-input)
</quiz>

<quiz>
Come si accede a un output nominato da un processo?
- [ ] `processName.outputName`
- [ ] `processName.get(outputName)`
- [x] `processName.out.outputName`
- [ ] `output.processName.outputName`

Per approfondire: [4.1.2. Emettere il file di report e nominare gli output](#412-emettere-il-file-di-report-e-nominare-gli-output)
</quiz>

<quiz>
Qual è la sintassi corretta per nominare un output in un processo?
- [ ] `name: outputName`
- [ ] `output: outputName`
- [x] `emit: outputName`
- [ ] `label: outputName`

Per approfondire: [4.1.2. Emettere il file di report e nominare gli output](#412-emettere-il-file-di-report-e-nominare-gli-output)
</quiz>

<quiz>
Quando si forniscono input multipli a un processo, cosa deve essere vero?
- [ ] Tutti gli input devono essere dello stesso tipo
- [ ] Gli input devono essere forniti in ordine alfabetico
- [x] L'ordine degli input deve corrispondere all'ordine definito nel blocco input
- [ ] Possono essere forniti solo due input alla volta

Per approfondire: [3. Passare parametri aggiuntivi a un processo](#3-passare-parametri-aggiuntivi-a-un-processo)
</quiz>
