# Parte 5: Hello Containers

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/5PyOWjKnNmg?si=QinuAnFwFj-Z8CrO&amp;list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } Guarda [tutta la playlist](https://www.youtube.com/playlist?list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik) sul canale Youtube Nextflow.

:green_book: La trascrizione di questo video è disponibile [qui](./transcripts/04_hello_modules.md).
///

Nella parte 1-4 del corso di training, hai imparato come usare gli elementi di base di Nextflow per assemblare un semplice workflow capace di processare alcuni testi, a fare esecuzioni parallele se ci sono più input e collezionare i risultati per esecuzioni future.

Tuttavia, si era limitati agli strumenti unix di base disponibili nel proprio ambiente.
Le attività del mondo reale spesso richiedono vari strumenti e pacchetti non inclusi default.
In genere, è necessario installare questi strumenti, gestire le loro dipendenze e risolvere eventuali conflitti.

Tutto ciò è molto noiso e fastidioso, quindi vi mostreremo come usare i container per risolvere questo problema in modo molto più comodo.

Un container è un'unità di software leggera, autonoma ed eseguibile creata da un'immagine di container che incluede tutto ciò che serve per eseguitre un'applicazione, compresi codici librerie di sistema e impostazioni.

!!! note

    We'll be teaching this using the technology [Docker](https://www.docker.com/get-started/), but Nextflow supports [several other container technologies](https://www.nextflow.io/docs/latest/container.html#) as well.

---

## 0. Warmup: Eseguire `hello-containers.nf`

Utilizzeremo lo script del workflow hello-containers.nf come punto di paretnza per la seconda sezione.

Just to make sure everything is working, run the script once before making any changes:
È equivalente allo script prodotto lavorando alla Parte 4 di questo corso di formazione.

```bash
nextflow run hello-containers.nf
```

Questo dovrebbe produrre il seguente risultato:

```console title="Output"
 N E X T F L O W   ~  version 24.10.0

Launching `hello-containers.nf` [tender_becquerel] DSL2 - revision: f7cat8e223

executor >  local (7)
[bd/4bb541] sayHello (1)         [100%] 3 of 3 ✔
[85/b627e8] convertToUpper (3)   [100%] 3 of 3 ✔
[7d/f7961c] collectGreetings     [100%] 1 of 1 ✔
```

Come in precendenza, i file output si trovano nella directory `results` (specificata dalla direttiva `publishDir`).

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

    There may also be a file named `output.txt` left over if you worked through Part 2 in the same environment.

Se questo ha funzionato, siete pronti per imparare a usare i container.

---

## 1. Utilizzare un container 'manualmente'

Quello che vogliamo fare è aggiungere un passo al nostro workflow che utilizzerà un container per l'esecuzione.

Tuttavia, prima di iniziare a usarli in Nextflow, esamineremo alcuni concetti e operazioni di base, per consolidare la comprensione di cosa sono i container.

### 1.1. Estrarre l'immagine del container

Per utilizzare un container, di solito si scarica o si 'pull' un'immagine del container da un reggistro dei container e poi si eseguel'immagien del container per creare un'istanza del container.

La sintassi generale è la seguente:

```bash title="Syntax"
docker pull '<container>'
```

La parte `docker pull` è l'istruzione al sistema di container di prelevare un'immagine di container da un repository.

La parte `'<container>'` è l'indirizzo URI dell'imagine del conatiner.

a titolo d'esempio estraiamo un'immagine ocntainer che contenga [cowpy](https://github.com/jeffbuttars/cowpy), un'implementazione di Python di uno strumento chiamato cowsay che genera arte ASCII per visualizzare in modo divertente input di testo arbitrari.

Esistono vari repository dove è possibile trovare i container pubblicati.
Noi abbiamo utilizzato il servizio di [Seqera Containers](https://seqera.io/containers/) per generare quest'immagine Docker container dal pacchetto cowpy Conda: `'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'`.

Eseguire il comando di pull richiesto:

```bash
docker pull 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
```

In questo modo si ottiene il seguente output di console mebntro il sistema scarica le immagini:

```console title="Output"
Unable to find image 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273' locally
131d6a1b707a8e65: Pulling from library/cowpy
dafa2b0c44d2: Pull complete
dec6b097362e: Pull complete
f88da01cff0b: Pull complete
4f4fb700ef54: Pull complete
92dc97a3ef36: Pull complete
403f74b0f85e: Pull complete
10b8c00c10a5: Pull complete
17dc7ea432cc: Pull complete
bb36d6c3110d: Pull complete
0ea1a16bbe82: Pull complete
030a47592a0a: Pull complete
622dd7f15040: Pull complete
895fb5d0f4df: Pull complete
Digest: sha256:fa50498b32534d83e0a89bb21fec0c47cc03933ac95c6b6587df82aaa9d68db3
Status: Downloaded newer image for community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273
community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273
```

Una volta completato il download, si ha una copia locale dell'immagine del container.

### 1.2. Utilizzare il container per eseguire cowpy come comando singolo

Un modo molto comune di utilizzare i container è quello di eseguirli direttamente, cioè in modo non interattivo.
Questo è ottimo per eseguire comandi una tantum.

La sintassi generale è la seguente:

```bash title="Syntax"
docker run --rm '<container>' [tool command]
```

La parte `docker run --rm '<container>'` è l'istruzione al sistema di container di creare un'istanza da un'immagine di container ed eseguire un comando in essa.
Il flag `--rm` indica al sistema di chiudere l'istanza del container al termine del comando.

La sintassi di `[tool command]` dipende dallo strumento che si sta usando e da come è impostato il conatiner.
Cominciamo con `cowpy`.

Completamente assemblato, il comadno di esecuzione del container si presenta come segue:

```bash
docker run --rm 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273' cowpy
```

Esegui per produrre il seguente output:

```console title="Output"
 ______________________________________________________
< Cowacter, eyes:default, tongue:False, thoughts:False >
 ------------------------------------------------------
     \   ^__^
      \  (oo)\_______
         (__)\       )\/\
           ||----w |
           ||     ||
```

Il sistema avvia il container, esegue il comando cowpy con i suoi parametri, invia l'output alla console e infine chiude l'istanza del container.

### 1.3. Ustilizzare il container per eseguire cowpy in modo interattivo

E anche possibile eseguire un conatiner in modo interattivo, in modo da avere un prompt di shell all'interno del container e poter giocare con i comandi.

#### 1.3.1. Avviare il container

Per eseguire il container in modo interattivo, basta aggiungere `-it` al comando `docker run'.
Facoltativamente, possiamo specificare la shell che vogliamo usare all'inetrno del conatiner, aggiungendo ad esempio `/bin/bash` al comando.

```bash
docker run --rm -it 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273' /bin/bash
```

Notate che il prompt cambia in qualcosa di simile a`(base) root@b645838b3314:/tmp#`, che indica che ora siete all'interno del conatienr.

E' possibile verificarlo eseguendo 'ls'per elencare il contenuto della directory:

```bash
ls /
```

```console title="Output"
bin    dev    etc    home   lib    media  mnt    opt    proc   root   run    sbin   srv    sys    tmp    usr    var
```

Si può notare che il filesystem all'interno del conatiner è diverso da quello del sistema host. .

!!! note

    Quando esegui un container, questo è isolato dal sistema host per impostazione predefinita.
    Ciò significa che il container non può accedere ad alcun file sul sistema host a meno che tu non glielo consenta esplicitamente.


    Imparerai come farlo in un minuto.

#### 1.3.2. Eseguire i comandi dello strumento desiderato

Ora che sei all'interno del container, puoi eseguire direttamente il comando `cowpy` e fornirgli alcuni parametri.
Ad esempio, la documentazione dello strumento afferma che possiamo modificare il carattere ('cowacter') con `-c`.

```bash
cowpy "Hello Containers" -c tux
```

Ora l'output mostra il pinguino di Linux, Tux, invece della mucca predefinita, perché abbiamo specificato il parametro `-c tux`.

```console title="Output"
 __________________
< Hello Containers >
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

Poiché ti trovi all'interno del container, puoi eseguire il comando cowpy tutte le volte che vuoi, variando i parametri di input, senza dover utilizzare i comandi Docker.

!!! Tip

    Utilizzare il flag '-c' per selezionare un carattere diverso, incluso:
    `beavis`, `cheese`, `daemon`, `dragonandcow`, `ghostbusters`, `kitty`, `moose`, `milk`, `stegosaurus`, `turkey`, `turtle`, `tux`

Questo è carino. Sarebbe ancora più carino se potessimo alimentare il nostro `greetings.csv` come input in questo.
Ma poiché non abbiamo accesso al file system, non possiamo farlo.

Risolviamolo.

#### 1.3.3. Uscire dal container

Per uscire dal container, puoi digitare `exit` al prompt oppure usare la scorciatoia da tastiera ++ctrl+d++.

```bash
exit
```

Il prompt dovrebbe ora tornare a essere quello che era prima dell'avvio del container.

#### 1.3.4. Montare i dati nel container

Quando si esegue un container, questo è isolato dal sistema host per impostazione predefinita.
Ciò significa che il container non può accedere ad alcun file sul sistema host a meno che non gli venga esplicitamente consentito di farlo.

Un modo per farlo è **montare** un **volume** dal sistema host nel container utilizzando la seguente sintassi:

```bash title="Syntax"
-v <outside_path>:<inside_path>
```

Nel nostro caso `<outside_path>` sarà la directory di lavoro corrente, quindi possiamo semplicemente usare un punto (`.`), e `<inside_path>` è solo un nome che inventiamo; Chiamiamolo `/data`.

Per montare un volume, sostituiamo i percorsi e aggiungiamo l'argomento di montaggio del volume al comando docker run come segue:

```bash
docker run --rm -it -v .:/data 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273' /bin/bash
```

Questo monta la directory di lavoro corrente come un volume che sarà accessibile sotto `/data` all'interno del container.

Puoi verificare che funzioni elencando il contenuto di `/data`:

```bash
ls /data
```

A seconda della parte di questa formazione che hai già svolto in precedenza, il risultato riportato di seguito potrebbe essere leggermente diverso.

```console title="Output"
demo-params.json  hello-channels.nf  hello-workflow.nf  modules          results
greetings.csv     hello-modules.nf   hello-world.nf     nextflow.config  work
```

<!-- ls output may need to be updated -->

Ora puoi vedere il contenuto della directory `data` dall'interno del container, incluso il file `greetings.csv`.

Questo ha effettivamente creato un tunnel attraverso il muro del container che puoi usare per accedere a quella parte del tuo file system.

#### 1.3.5. Utilizza i dati montati

Ora che abbiamo montato la directory `data` nel container, possiamo usare il comando `cowpy` per visualizzare il contenuto del file `greetings.csv`.

Per fare questo, useremo`cat /data/greetings.csv | ` per inviare il contenuto del file CSV al comando `cowpy`.

```bash
cat /data/greetings.csv | cowpy -c turkey
```

Questo produce l'ASCII art desiderata di un tacchino che snocciola i nostri saluti di esempio:

```console title="Output"
 _________
/ Hello   \
| Bonjour |
\ Holà    /
 ---------
  \                                  ,+*^^*+___+++_
   \                           ,*^^^^              )
    \                       _+*                     ^**+_
     \                    +^       _ _++*+_+++_,         )
              _+^^*+_    (     ,+*^ ^          \+_        )
             {       )  (    ,(    ,_+--+--,      ^)      ^\
            { (\@)    } f   ,(  ,+-^ __*_*_  ^^\_   ^\       )
           {:;-/    (_+*-+^^^^^+*+*<_ _++_)_    )    )      /
          ( /  (    (        ,___    ^*+_+* )   <    <      \
           U _/     )    *--<  ) ^\-----++__)   )    )       )
            (      )  _(^)^^))  )  )\^^^^^))^*+/    /       /
          (      /  (_))_^)) )  )  ))^^^^^))^^^)__/     +^^
         (     ,/    (^))^))  )  ) ))^^^^^^^))^^)       _)
          *+__+*       (_))^)  ) ) ))^^^^^^))^^^^^)____*^
          \             \_)^)_)) ))^^^^^^^^^^))^^^^)
           (_             ^\__^^^^^^^^^^^^))^^^^^^^)
             ^\___            ^\__^^^^^^))^^^^^^^^)\\
                  ^^^^^\uuu/^^\uuu/^^^^\^\^\^\^\^\^\^\
                     ___) >____) >___   ^\_\_\_\_\_\_\)
                    ^^^//\\_^^//\\_^       ^(\_\_\_\)
                      ^^^ ^^ ^^^ ^
```

Sentiti libero di giocare con questo comando.
Quando hai finito, esci dal container come in precedenza:

```bash
exit
```

Ti ritroverai nel tuo guscio normale.

### Takeaway

Sai come estrarre uncontainer ed eseguirlo singolarmente o in modo interattivo. Sai anche come rendere accessibili i tuoi dati dall'interno del container, il che ti consente di provare qualsiasi strumento a cui sei interessato su dati reali senza dover installare alcun software sul tuo sistema.

### Prossimo passo?

Scopri come utilizzare i contenitori per l'esecuzione dei processi Nextflow.

---

## 2. Utilizzare i container in Nextflow

Nextflow ha un supporto integrato per l'esecuzione di processi all'interno di contenitori per consentirti di eseguire strumenti che non hai installato nel tuo ambiente di elaborazione.
Ciò significa che puoi utilizzare qualsiasi immagine di container desideri per eseguire i tuoi processi e Nextflow si occuperà di estrarre l'immagine, montare i dati ed eseguire il processo al suo interno.

Per dimostrarlo, aggiungeremo un passaggio `cowpy` alla pipeline che stiamo sviluppando, dopo il passaggio `collectGreetings`.

### 2.1. Scrivi un modulo `cowpy`

#### 2.1.1. Crea uno stub di file per il nuovo modulo

Crea un file vuoto per il modulo chiamato `cowpy.nf`.

```bash
touch modules/cowpy.nf
```

Questo ci fornisce un posto dove mettere il codice del processo.

#### 2.1.2. Copia il codice del processo `cowpy` nel file del modulo

Possiamo modellare il nostro processo `cowpy` sugli altri processi che abbiamo scritto in precedenza.

```groovy title="modules/cowpy.nf" linenums="1"
#!/usr/bin/env nextflow

// Generate ASCII art with cowpy
process cowpy {

    publishDir 'results', mode: 'copy'

    input:
        path input_file
        val character

    output:
        path "cowpy-${input_file}"

    script:
    """
    cat $input_file | cowpy -c "$character" > cowpy-${input_file}
    """

}
```

L'output sarà un nuovo file di testo contenente l'ASCII art generata dallo strumento `cowpy`.

### 2.2. Aggiungi cowpy al flusso di lavoro

Ora dobbiamo importare il modulo e chiamare il processo.

#### 2.2.1. Importa il processo `cowpy` in `hello-containers.nf`

Inserire la dichiarazione di importazione sopra il blocco del workflow e compilarla in modo appropriato.

=== "Dopo"

    ```groovy title="hello-containers.nf" linenums="9" hl_lines="5"
    // Include modules
    include { sayHello } from './modules/sayHello.nf'
    include { convertToUpper } from './modules/convertToUpper.nf'
    include { collectGreetings } from './modules/collectGreetings.nf'
    include { cowpy } from './modules/cowpy.nf'

    workflow {
    ```

=== "Prima"

    ```groovy title="hello-containers.nf" linenums="9"
    // Include modules
    include { sayHello } from './modules/sayHello.nf'
    include { convertToUpper } from './modules/convertToUpper.nf'
    include { collectGreetings } from './modules/collectGreetings.nf'

    workflow {
    ```

#### 2.2.2. Aggiungere una chiamata al processo `cowpy` nel flusso di lavoro

Colleghiamo il processo `cowpy()` all'output del processo `collectGreetings()`, che come ricorderai produce due output:

- `collectGreetings.out.outfile` contiene gli otput dei file
- `collectGreetings.out.count`contiene il conteggio dei saluti per batch

Nel blocco del flusso di lavoro, apportare la seguente modifica al codice:

=== "Dopo"

    ```groovy title="hello-containers.nf" linenums="28" hl_lines="7 8"
        // collect all the greetings into one file
        collectGreetings(convertToUpper.out.collect(), params.batch)

        // emit a message about the size of the batch
        collectGreetings.out.count.view{ num_greetings -> "There were $num_greetings greetings in this batch" }

        // generate ASCII art of the greetings with cowpy
        cowpy(collectGreetings.out.outfile, params.character)
    ```

=== "Prima"

    ```groovy title="hello-containers.nf" linenums="28"
        // collect all the greetings into one file
        collectGreetings(convertToUpper.out.collect(), params.batch)

        // emit a message about the size of the batch
        collectGreetings.out.count.view{ num_greetings -> "There were $num_greetings greetings in this batch" }
    ```

Si noti che includiamo un nuovo parametro CLI, `params.character`, per specificare quale carattere vogliamo che dica i saluti.

#### 2.2.3. Imposta un valore predefinito per `params.character`

Ci piace essere pigri e saltare la digitazione dei parametri nelle nostre righe di comando.

=== "Dopo"

    ```groovy title="hello-containers.nf" linenums="3" hl_lines="6"
    /*
     * Pipeline parameters
     */
    params.greeting = 'greetings.csv'
    params.batch = 'test-batch'
    params.character = 'turkey'
    ```

=== "Prima"

    ```groovy title="hello-containers.nf" linenums="3"
    /*
     * Pipeline parameters
     */
    params.greeting = 'greetings.csv'
    params.batch = 'test-batch'
    ```

Dovrebbe essere tutto ciò di cui abbiamo bisogno per far funzionare tutto.

#### 2.2.4. Eseguite il workflow per verificarne il funzionamento

Eseguire questa operazione con il flag `-resume`.

```bash
nextflow run hello-containers.nf -resume
```

Oh no, c'è un errore!

```console title="Output"
 N E X T F L O W   ~  version 24.10.0

Launching `hello-containers.nf` [special_lovelace] DSL2 - revision: 028a841db1

executor >  local (1)
[f6/cc0107] sayHello (1)       | 3 of 3, cached: 3 ✔
[2c/67a06b] convertToUpper (3) | 3 of 3, cached: 3 ✔
[1a/bc5901] collectGreetings   | 1 of 1, cached: 1 ✔
[b2/488871] cowpy             | 0 of 1
There were 3 greetings in this batch
ERROR ~ Error executing process > 'cowpy'

Caused by:
  Process `cowpy` terminated with an error exit status (127)

Command executed:

  cat COLLECTED-test-batch-output.txt | cowpy -c "turkey" > cowpy-COLLECTED-test-batch-output.txt

Command exit status:
  127

Command output:
  (empty)

Command error:
  .command.sh: line 2: cowpy: command not found

(trimmed output)
```

Questo codice di errore, `error exit status (127)`, significa che l'eseguibile richiesto non è stato trovato.
Naturalmente, dato che stiamo chiamando lo strumento `cowpy` ma non abbiamo ancora specificato un container.

### 2.3. Utilizzare un container per l'esecuzione

Dobbiamo specificare un container e dire a Nextflow di usarlo per il processo `cowpy()`.

#### 2.3.1. Specificare un container per il processo `cowpy` da utilizzare

Modificare il modulo `cowpy.nf` per aggiungere la direttiva `container` alla definizione del processo come segue:

=== "Dopo"

    ```groovy title="modules/cowpy.nf" linenums="4" hl_lines="4"
    process cowpy {

        publishDir 'containers/results', mode: 'copy'
        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
    ```

=== "Prima"

    ```groovy title="modules/cowpy.nf" linenums="4"
    process cowpy {

        publishDir 'containers/results', mode: 'copy'
    ```

indica a Nextflow che, se l'uso di Docker è abilitato, deve usare l'immagine del container specificata qui per eseguire il processo.

#### 2.3.2. Abilitare l'uso di Docker tramite il file `nextflow.config

Qui anticipiamo leggermente l'argomento della prossima e ultima parte di questo corso (Parte 6), che riguarda la configurazione.

Uno dei modi principali che Nextflow offre per configurare l'esecuzione del workflow è l'uso di un file `nextflow.config`. Quando un file di questo tipo è presente nella directory corrente, Nextflow lo caricherà automaticamente e applicherà la configurazione in esso contenuta.
Abbiamo fornito un file `nextflow.config` con una singola riga di codice che disabilita Docker: `docker.enabled = false`.

Ora, passiamo a `true` per abilitare Docker:

=== "Dopo"

    ```console title="nextflow.config" linenums="1" hl_lines="1"
    docker.enabled = true
    ```

=== "Prima"

    ```console title="nextflow.config" linenums="1"
    docker.enabled = false
    ```

!!! note

    It is possible to enable Docker execution from the command-line, on a per-run basis, using the `-with-docker <container>` parameter.
    However, that only allows us to specify one container for the entire workflow, whereas the approach we just showed you allows us to specify a different container per process.
    This is better for modularity, code maintenance and reproducibility.

#### 2.3.3. Eseguire il workflow con Docker abilitato

Eseguire il workflow con il flag `resume`:

```bash
nextflow run hello-containers.nf -resume
```

Questa volta funziona davvero.

```console title="Output" linenums="1"
 N E X T F L O W   ~  version 24.10.0

Launching `hello-containers.nf` [elegant_brattain] DSL2 - revision: 028a841db1

executor >  local (1)
[95/fa0bac] sayHello (3)       | 3 of 3, cached: 3 ✔
[92/32533f] convertToUpper (3) | 3 of 3, cached: 3 ✔
[aa/e697a2] collectGreetings   | 1 of 1, cached: 1 ✔
[7f/caf718] cowpy              | 1 of 1 ✔
There were 3 greetings in this batch
```

È possibile trovare l'output di cowpy nella directory `results`.

```console title="results/cowpy-COLLECTED-test-batch-output.txt"
 _________
/ HOLà    \
| HELLO   |
\ BONJOUR /
 ---------
  \                                  ,+*^^*+___+++_
   \                           ,*^^^^              )
    \                       _+*                     ^**+_
     \                    +^       _ _++*+_+++_,         )
              _+^^*+_    (     ,+*^ ^          \+_        )
             {       )  (    ,(    ,_+--+--,      ^)      ^\
            { (\@)    } f   ,(  ,+-^ __*_*_  ^^\_   ^\       )
           {:;-/    (_+*-+^^^^^+*+*<_ _++_)_    )    )      /
          ( /  (    (        ,___    ^*+_+* )   <    <      \
           U _/     )    *--<  ) ^\-----++__)   )    )       )
            (      )  _(^)^^))  )  )\^^^^^))^*+/    /       /
          (      /  (_))_^)) )  )  ))^^^^^))^^^)__/     +^^
         (     ,/    (^))^))  )  ) ))^^^^^^^))^^)       _)
          *+__+*       (_))^)  ) ) ))^^^^^^))^^^^^)____*^
          \             \_)^)_)) ))^^^^^^^^^^))^^^^)
           (_             ^\__^^^^^^^^^^^^))^^^^^^^)
             ^\___            ^\__^^^^^^))^^^^^^^^)\\
                  ^^^^^\uuu/^^\uuu/^^^^\^\^\^\^\^\^\^\
                     ___) >____) >___   ^\_\_\_\_\_\_\)
                    ^^^//\\_^^//\\_^       ^(\_\_\_\)
                      ^^^ ^^ ^^^ ^
```

Si vede che il personaggio sta pronunciando tutti i saluti, proprio come ha fatto quando abbiamo eseguito il comando `cowpy` sul file `greetings.csv` dall'interno del container.

<!-- considering a side quest where we show how to use a conditional to skip the collect step if we want to emit the cowpy'ed greetings individually, and how to use metadata management to assign a specific character to each greeting, maybe do some cross products etc -->

#### 2.3.4. Controllare come Nextflow ha lanciato il task containerizzato

Diamo un'occhiata alla sottodirectory di lavoro di una delle chiamate al processo cowpy per capire meglio come Nextflow lavora con i container.

Controllare l'output del comando nextflow run per trovare l'ID della chiamata al processo cowpy.
Quindi navigare nella sottodirectory work.
Al suo interno si trova il file .command.run che contiene tutti i comandi eseguiti da Nextflow per conto dell'utente durante l'esecuzione della pipeline.
Aprite il file .command.run e cercate nxf_launch; dovreste vedere qualcosa di simile:

```bash
nxf_launch() {
    docker run -i --cpu-shares 1024 -e "NXF_TASK_WORKDIR" -v /workspaces/training/hello-nextflow/work:/workspaces/training/hello-nextflow/work -w "$NXF_TASK_WORKDIR" --name $NXF_BOXID community.wave.seqera.io/library/pip_cowpy:131d6a1b707a8e65 /bin/bash -ue /workspaces/training/hello-nextflow/work/7f/caf7189fca6c56ba627b75749edcb3/.command.sh
}
```

Come si può vedere, Nextflow utilizza il comando docker run per lanciare la chiamata di processo.
Inoltre, monta la corrispondente sottodirectory di lavoro nel container, imposta di conseguenza la directory di lavoro all'interno del container ed esegue il nostro script bash templato nel file .command.sh.
Tutto il duro lavoro che abbiamo dovuto fare manualmente nella sezione precedente viene svolto da Nextflow!

### Takeaway

Sapete come usare i container in Nextflow per eseguire i processi.

### E ora?

Fate una pausa!
Quando sarete pronti, passate alla Parte 6 per imparare a configurare l'esecuzione della pipeline in base alla vostra infrastruttura e a gestire la configurazione di input e parametri.
È l'ultima parte e il gioco è fatto!
