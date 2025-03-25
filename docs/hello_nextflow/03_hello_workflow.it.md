# Parte 3: Hello Workflow

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/zJP7cUYPEbA?si=Irl9nAQniDyICp2b&amp;list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } Guarda [l'intera playlist sul canale Youtube di Nextflowl](https://www.youtube.com/playlist?list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik).
///

La maggior parte dei flussi di lavoro nel mondo reale coinvolge più di un passaggio. 
In questo modulo di formazione, imparerai come connettere i processi insieme in un flusso di lavoro a più fasi.

Questo ti insegnerà il modo Nextflow per ottenere i seguenti obiettivi:

1. Far fluire i dati da un processo al successivo
2. Raccogliere gli output da più chiamate di processo in una singola chiamata
3. Passare più di un input a un processo
4. Gestire più output provenienti da un processo

Per dimostrare ciò, continueremo a costruire sull'esempio del dominio indipendente Hello World delle Parti 1 e 2. 
Questa volta, apporteremo le seguenti modifiche al nostro flusso di lavoro per riflettere meglio su come le persone costruiscono flussi di lavoro reali:

1. Aggiungere un secondo passaggio che converte il saluto in maiuscolo.
2. Aggiungere un terzo passaggio che raccoglie tutti i saluti trasformati e li scrive in un unico file.
3. Aggiungere un parametro per nominare il file di output finale e passarlo come input secondario al passaggio di raccolta.
4. Far sì che il passaggio di raccolta produca anche una statistica semplice su ciò che è stato elaborato.

---

## 0. Warmup: Run 'hello-workflow.nf'

Useremo lo script del workflow 'hello-workflow.nf' come punto di partenza.
Esso è equivalente allo script prodotto lavorando attraverso la Parte 2 di questo corso di formazione.

Per essere sicuri che tutto funzioni, esegui lo script una volta prima di apportare qualsiasi modifica

```bash
nextflow run hello-workflow.nf
```

```console title="Output"
 N E X T F L O W   ~  version 24.10.0

Launching `hello-workflow.nf` [stupefied_sammet] DSL2 - revision: b9e466930b

executor >  local (3)
[2a/324ce6] sayHello (3) | 3 of 3 ✔
```

Come in precedenza, troverai i file di output nella directory 'results' (specificata dalla direttiva 'publishDir').


```console title="Directory contents"
results
├── Bonjour-output.txt
├── Hello-output.txt
└── Holà-output.txt
```

!!! nota

    Potrebbe esserci anche un file chiamato 'output.txt' se hai lavorato attraverso la Parte 2 nello stesso ambiente.

Se tutto ha funzionato, sei pronto per imparare come assemblare un flusso di lavoro a più fasi.

---

## 1. Aggiungi un secondo step al workflow
Aggiungeremo un passaggio per convertire il saluto in maiuscolo.
A tal fine, dobbiamo fare tre cose:

-Definire il comando che useremo per eseguire la conversione in maiuscolo.
-Scrivere un nuovo processo che racchiuda il comando per la conversione in maiuscolo.
-Chiamare il nuovo processo nel blocco del flusso di lavoro e configurarlo per prendere l'output del processo 'sayHello()' come input.

## 1.1 Definire il comando per la conversione in maiuscolo e testarlo nel terminale

Per eseguire la conversione dei saluti in maiuscolo, useremo uno strumento UNIX classico chiamato 'tr' per 'text replacement' (sostituzione del testo), con la seguente sintassi:

```bash title="Syntax"
tr '[a-z]' '[A-Z]'
```

Questa è una sostituzione di testo molto semplice che non tiene conto delle lettere accentate, quindi ad esempio 'Holà' diventerà 'HOLà', ma andrà bene lo stesso per dimostrare i concetti di Nextflow, e questo è ciò che conta.

Per testarlo, possiamo eseguire il comando 'echo 'Hello World'' e passare il suo output al comando 'tr':

```bash
echo 'Hello World' | tr '[a-z]' '[A-Z]' > UPPER-output.txt
```

Il risultato è un file di testo chiamato 'UPPER-output.txt' che contiene la versione in maiuscolo della stringa 'Hello World':

```console title="UPPER-output.txt"
HELLO WORLD
```

Questo è ciò che proveremo a fare con il nostro flusso di lavoro.

### 1.1 Scrivere il passaggio di conversione in maiuscolo come un processo Nextflow

Possiamo modellare il nostro nuovo processo sul primo, poiché vogliamo usare gli stessi componenti.

Aggiungi la seguente definizione del processo allo script del flusso di lavoro:

```groovy title="hello-workflow.nf" linenums="22"
/*
 * Use a text replacement tool to convert the greeting to uppercase
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

Qui, componiamo il secondo nome del file di output in base al nome del file di input, in modo simile a quanto abbiamo fatto inizialmente per l'output del primo processo.

!!!nota 

   Nextflow determinerà l'ordine delle operazioni in base alla concatenazione degli input e degli output, quindi l'ordine delle definizioni dei processi nello script del flusso di lavoro non è importante.
   Tuttavia, ti consigliamo di essere gentile con i tuoi collaboratori e con il futuro te stesso, e cercare di scriverle in un ordine logico per motivi di leggibilità."

### 1.2 Aggiungi una chiamata al nuovo processo nel blocco del flusso di lavoro

Ora dobbiamo dire a Nextflow di chiamare effettivamente il processo che abbiamo appena definito.

Nel blocco del flusso di lavoro, apporta la seguente modifica al codice:

_Prima:_

```groovy title="hello-workflow.nf" linenums="53"
    // emit a greeting
    sayHello(greeting_ch)
}
```

_Dopo:_

```groovy title="hello-workflow.nf" linenums="53"
    // emit a greeting
    sayHello(greeting_ch)

    // convert the greeting to uppercase
    convertToUpper()
}
```

Questo non è ancora funzionante perché non abbiamo specificato cosa deve essere l'input per il processo 'convertToUpper()'.

### 1.3 Passare l'output del primo processo al secondo processo

Ora dobbiamo fare in modo che l'output del processo 'sayHello()' fluisca nel processo 'convertToUpper()'.

Comodamente, Nextflow impacchetta automaticamente l'output di un processo in un canale chiamato '<process>.out'.
Quindi, l'output del processo 'sayHello' è un canale chiamato 'sayHello.out', che possiamo collegare direttamente alla chiamata a 'convertToUpper()'.

Nel blocco del flusso di lavoro, apporta la seguente modifica al codice:

_Prima:_

```groovy title="hello-workflow.nf" linenums="56"
    // convert the greeting to uppercase
    convertToUpper()
}
```

_Dopo:_

```groovy title="hello-workflow.nf" linenums="56"
    // convert the greeting to uppercase
    convertToUpper(sayHello.out)
}
```

Per un caso semplice come questo (un output a un input), è tutto ciò che dobbiamo fare per connettere due processi!

### 1.4 Esegui di nuovo il flusso di lavoro con -resume

Eseguiamo di nuovo il flusso di lavoro utilizzando il flag '-resume', poiché abbiamo già eseguito con successo il primo passaggio del flusso di lavoro.

```bash
nextflow run hello-workflow.nf -resume
```

Dovresti vedere il seguete output:

```console title="Output" linenums="1"
 N E X T F L O W   ~  version 24.10.0

Launching `hello-workflow.nf` [disturbed_darwin] DSL2 - revision: 4e252c048f

executor >  local (3)
[79/33b2f0] sayHello (2)       | 3 of 3, cached: 3 ✔
[b3/d52708] convertToUpper (3) | 3 of 3 ✔
```

Ora c'è una riga extra nell'output della console (riga 7), che corrisponde al nuovo processo che abbiamo appena aggiunto.

Diamo un'occhiata all'interno della directory di lavoro di una delle chiamate al secondo processo.

```console title="Directory contents"
work/b3/d52708edba8b864024589285cb3445/
├── Bonjour-output.txt -> /workspaces/training/hello-nextflow/work/79/33b2f0af8438486258d200045bd9e8/Bonjour-output.txt
└── UPPER-Bonjour-output.txt
```

Troviamo due file di output: l'output del primo processo e l'output del secondo.

L'output del primo processo è lì perché Nextflow lo ha messo in quella directory per avere tutto il necessario per l'esecuzione all'interno della stessa sottodirectory.
Tuttavia, in realtà si tratta di un collegamento simbolico che punta al file originale nella sottodirectory della prima chiamata al processo.
Per impostazione predefinita, quando si esegue su una singola macchina, come stiamo facendo qui, Nextflow utilizza collegamenti simbolici anziché copie per mettere in scena i file di input e i file intermedi.

Troverai anche i file di output finali nella directory 'results', poiché abbiamo usato la direttiva 'publishDir' anche nel secondo processo.

```console title="Directory contents"
results
├── Bonjour-output.txt
├── Hello-output.txt
├── Holà-output.txt
├── UPPER-Bonjour-output.txt
├── UPPER-Hello-output.txt
└── UPPER-Holà-output.txt
```

Pensiamo a come tutto ciò che abbiamo fatto è stato connettere l'output di 'sayHello' all'input di 'convertToUpper' e i due processi potrebbero essere eseguiti in serie.
Nextflow ha fatto il lavoro difficile di gestire i singoli file di input e output e passarli tra i due comandi per noi.

Questa è una delle ragioni per cui i canali di Nextflow sono così potenti: si occupano del lavoro noioso coinvolto nel connettere i passaggi del flusso di lavoro.

### Conclusioni

Ora sai come aggiungere un secondo passaggio che prende l'output del primo passaggio come input.

### Cosa c'è dopo?

Impara come raccogliere gli output da chiamate di processo in batch e passarli a un singolo processo.

---







