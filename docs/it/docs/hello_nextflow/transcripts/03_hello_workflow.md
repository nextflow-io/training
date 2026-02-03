# Parte 3: Hello Workflow - Trascrizione

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduzione assistita da IA - [scopri di più e suggerisci miglioramenti](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/zJP7cUYPEbA?si=Irl9nAQniDyICp2b&amp;list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Note importanti"

    Questa pagina mostra solo la trascrizione. Per istruzioni complete passo dopo passo, tornare al [materiale del corso](../03_hello_workflow.md).

    I numeri di sezione mostrati nella trascrizione sono forniti solo a scopo indicativo e potrebbero non includere tutti i numeri di sezione presenti nei materiali.

## Benvenuto

Salve, benvenuti alla parte tre del corso di formazione "Hello Nextflow".

Questo capitolo si intitola "Hello Workflow".

Nel capitolo due, abbiamo costruito un semplice workflow di un processo, ma in realtà, i pipeline sono utili perché possono concatenare più passaggi di analisi insieme.

In questo capitolo, prenderemo quell'esempio iniziale e lo estenderemo per renderlo un po' più realistico.

Aggiungeremo alcuni passaggi aggiuntivi e vedremo come utilizziamo i canali per collegare questi passaggi.

Esamineremo più attività, che possono essere raccolte in un singolo processo, e vedremo i processi che possono avere più input e più output.

Bene, iniziamo.

Quindi, cominciamo. Come prima. Andiamo su training.nextflow.io. Hello Nextflow, capitolo tre. Hello Workflow. E apriamo il nostro spazio di lavoro. Ho pulito tutti i miei file di lavoro dai capitoli precedenti e aprirò Hello Workflow.

Ora questo è lo stesso file su cui abbiamo lavorato finora, quindi dovrebbe sembrare familiare. Abbiamo il nostro processo say hello. Abbiamo il nostro params.greeting con il suo file greetings CSV, e abbiamo il nostro workflow in fondo, che carica quel file CSV, crea il canale e lo passa al nostro processo.

## 0. Riscaldamento: Eseguire hello-workflow.nf

Se volete, possiamo provarlo e verificare che funzioni come previsto. Aprite un terminale e digitate nextflow run hello workflow nf e premete invio.

Bene, ottimo. I nostri tre processi vengono eseguiti. Abbiamo la nostra directory results con i nostri tre output. Bonjour. Hello. Holà. Quindi chiudiamo quei file, chiudiamo il terminale e torniamo allo script.

## 1. Aggiungere un secondo passaggio al workflow

Okay. Per il nostro esempio, rimaniamo semplici e cerchiamo di rimanere agnostici dal dominio. Quindi il nostro secondo processo manipolerà semplicemente queste stringhe, queste parole, in modo semplice. Useremo il comando Unix translate per prendere questi file e renderli tutti maiuscoli. Lo facciamo con il comando "tr".

## 1.1. Definire il comando di conversione in maiuscolo e testarlo nel terminale

Possiamo provare questo direttamente nel terminale bash e vedere se funziona. Quindi digitate echo, Hello World, e poi passatelo con il carattere pipe a tr, e gli diamo un pattern di riconoscimento, a-z e a cosa dovrebbe tradurlo. A-Z in maiuscolo.

Questo è molto semplice perché sta letteralmente convertendo i caratteri A-Z. Quindi non funzionerà su nulla che abbia accenti o cose simili. Ma ai fini dell'esempio, dovreste capire l'idea.

Premo invio e stampa sul terminale HELLO WORLD in maiuscolo. E proprio come prima, potremmo reindirizzare questo a un file se volessimo. Outfile.

Okay. Puliamo questo.

## 1.1. Scrivere il passaggio di conversione in maiuscolo come processo Nextflow

Torniamo al nostro script e scriviamo un nuovo processo per gestire questo comando bash. Copierò il processo precedente, lo incollerò sotto e lo chiamerò convert to upper. Per maiuscolo. Userò lo stesso publishDir results, ma farò alcuni cambiamenti qui. Invece di prendere un val, prenderò un path input file, e avrò un prefisso qui upper, in modo che i nostri file di output non sovrascrivano l'output. E userò il nome della variabile dall'input. E poi cambierò lo script qui sotto, e invece userò cat sul file di input e proprio come abbiamo fatto in Bash TR, a-z, upper input file .txt. Okay, clicchiamo su salva.

## 1.2. Aggiungere una chiamata al nuovo processo nel blocco workflow

Ora se scorro verso il basso, dobbiamo effettivamente chiamare questo processo. Aggiungere semplicemente il processo in uno script non è sufficiente. Dobbiamo dire a Nextflow che dobbiamo eseguire questo processo e dove farlo.

Quindi scriverò qui, convert to upper e

okay, stiamo ricevendo un errore qui che dice che si aspetta un argomento. Certo, dobbiamo passare qualcosa a questo processo in modo che abbia effettivamente qualcosa da fare.

## 1.3. Passare l'output del primo processo al secondo processo

Quello che faremo è prendere l'output da questo processo. Quindi prendo il nome, say hello, e quando faccio dot out.

Per un esempio semplice come questo, dove abbiamo un processo che ha solo un output e lo passiamo a un nuovo processo, quindi ha un input, questo dovrebbe essere tutto ciò di cui abbiamo bisogno. Quindi cliccherò su salva, aprirò il terminale e proviamo a eseguirlo di nuovo.

## 1.4. Eseguire nuovamente il workflow

Ora, non ho ripulito la mia directory work dall'ultima volta che ho eseguito questo workflow. Lo eseguirò di nuovo e userò questa opportunità per mostrare come funziona il caching parziale. Quindi se faccio singolo trattino resume. Si spera che dovrebbe riutilizzare gli output da quel primo processo, che erano esattamente gli stessi dell'ultima volta che l'ho eseguito. Ma ora abbiamo un nuovo processo qui che non è stato eseguito prima, che viene eseguito da zero. E infatti, potete vedere che il primo processo ha utilizzato gli output memorizzati nella cache, e il secondo output ha eseguito tre di tre. Potete anche vedere che abbiamo entrambi i nostri processi qui ora, il nostro primo processo, say hello, eseguito tre volte, e il nostro secondo processo convert to upper eseguito tre volte.

Se eseguo questo di nuovo, come promemoria, con -ansi-log false, dovremmo vedere che vengono eseguite sei diverse attività di processo, tre per ciascuna di esse. Quindi questo sta facendo esattamente quello che speravamo. Il primo processo viene eseguito tre volte, passando quegli output a un secondo processo, che viene poi eseguito tre volte.

Quindi diamo un'occhiata all'interno della directory work e vediamo come Nextflow gestisce questi input di file. Se prendo questa directory hash qui dal secondo processo, possiamo usare di nuovo il comando tree con -a solo per guardare questi file. Potete vedere qui che abbiamo il nostro file di input, che è il file Bonjour-output.txt, e in realtà è un symlink. Questo è ciò che ci mostra questa freccia, e sta puntando al file nella directory work precedente.

Questo ha senso. Nextflow gestisce l'esecuzione di ogni attività nella propria directory incapsulata, quindi è completamente autonoma. Tuttavia, deve fornire i file da passaggi precedenti come input. Piuttosto che raggiungere al di fuori della directory work per ottenere quei file, Nextflow li prepara nella directory work.

Se abbiamo un file system condiviso come qui, lo fa utilizzando un symlink in modo da non utilizzare spazio file aggiuntivo. Se utilizziamo l'archiviazione cloud con bucket in posizioni diverse, recupererebbe quei file e li copierebbe effettivamente nella directory work.

Diamo un'occhiata al file command sh. Se faccio code work, command sh, potete vedere, infatti, sta accedendo a quel file dalla directory locale. Quindi tutto è molto autonomo e pulito.

Possiamo anche controllare la directory results e assicurarci che questi file siano stati emessi correttamente. E infatti, in results, possiamo vedere tutti i file di output dal primo processo e tutti i file di output dal secondo. E sono tutti in maiuscolo come speravamo.

Qui inizia a brillare la potenza di Nextflow. Con un codice molto minimale, Nextflow ha gestito l'esecuzione in parallelo di queste attività con incapsulamento pulito all'interno di directory work separate e preparazione dei file di input e output e pubblicazione dei file, tutto automaticamente per noi, semplicemente pronto all'uso. Quindi potete vedere come, man mano che aumentiamo la complessità dei nostri workflow di analisi, questa funzionalità sia davvero, davvero preziosa.

## 2. Aggiungere un terzo passaggio per raccogliere tutti i saluti

Okay. Questi passaggi erano uno a uno. Avevamo un output dal primo processo che andava a un input per il secondo processo. Successivamente, parleremo di come raccogliere questi diversi output in una singola attività di processo, che è di nuovo una cosa molto comune da fare. Quindi apriamo rapidamente il terminale e facciamo una prova a secco di questo.

## 2.1. Definire il comando di raccolta e testarlo nel terminale

Farò un trucco e copierò il codice bash di esempio dal materiale di formazione e premerò semplicemente invio.

Quello che possiamo vedere qui è che abbiamo eseguito questo comando echo tre volte su tre diversi file di output, che posso vedere qui. E poi abbiamo usato il comando cat per stampare l'output di ciascuno di questi tre diversi file e reindirizzarlo a un singolo file raccolto.

E se faccio "cat COLLECTED-output", potete vedere che ha il contenuto di quei tre diversi file, ora in un singolo file.

## 2.2. Creare un nuovo processo per eseguire il passaggio di raccolta

Quindi vediamo se possiamo replicare la stessa cosa all'interno del nostro pipeline Nextflow.

Scorriamo verso l'alto e creiamo un terzo processo. Copierò questo precedente, e questa volta lo chiamerò Collect Greetings.

Nel terminale bash, l'abbiamo chiamato collected output txt. Quindi dirò lo stesso path output qui. E farò il reindirizzamento qui, quindi viene salvato nello stesso modo.

Okay. Dobbiamo cambiare ciò che accade all'inizio di quel comando, e dobbiamo pensare a quale sia il file di input qui. In effetti, questo processo prenderà più file di input. Manterrò path e cambierò questo in una nuova variabile chiamata input files, plurale.

Poi farò di nuovo cat su di essi come abbiamo fatto nel nostro script bash. E userò qui la variabile.

Ora, potreste pensare che questo non funzionerebbe. Abbiamo visto in precedenza fallimenti in cui un array di stringhe o un array di percorsi è stato passato a un processo e ciò ha causato un errore. Ma in realtà, qui Nextflow gestirà questo automaticamente per noi nel modo giusto. Prenderà diversi file di input, e stamperà semplicemente i diversi percorsi dei file qui.

Naturalmente aiuta che il comando cat possa prendere una serie di nomi di file come questo. Se stessi usando un comando diverso che richiedesse un argomento prima di ogni percorso di file o qualcosa del genere, dovremmo avere un po' più di codice qui e logica per poter gestire l'iterazione di questi percorsi di file. Ma in questo caso, dovrebbe semplicemente funzionare.

## 2.3. Aggiungere il passaggio di raccolta al workflow

Okay, scendiamo al workflow e aggiungiamo il nostro nuovo processo. Collect greetings. E di nuovo, prendiamo l'output da convert to upper out. Salviamo questo.

Proviamolo. nextflow run hello workflow.

Okay, il workflow è stato eseguito, ma qualcosa è un po' strano qui. Abbiamo tre esecuzioni del primo passaggio, che ci aspettiamo. Tre attività per il secondo, ma abbiamo anche tre attività alla fine quando ci aspettavamo di avere solo una singola attività qui che unisce tutti gli output.

Se andiamo nella nostra directory results. Vediamo anche che l'output raccolto ha solo un singolo valore invece di tutti e tre. Questo perché quel file di output è stato sovrascritto tre volte con tre valori diversi.

Questo ha senso perché abbiamo passato un output a un input qui nello stesso modo in cui abbiamo fatto nel passaggio precedente.

## 2.4. Utilizzare un operatore per raccogliere i saluti in un singolo input

Quindi abbiamo bisogno di un operatore qui per prendere questo canale con tre elementi e ridurli a un singolo elemento, in modo che quel processo finale venga eseguito solo una volta.

Per farlo, useremo l'operatore collect. Posso farlo direttamente all'interno del workflow. Posso fare .out e concatenare un operatore qui alla fine .collect.

Clicco salva. E poi ai fini di questa formazione, farò anche alcuni operatori view come abbiamo fatto prima, così possiamo dare un'occhiata a questo canale prima e dopo aver utilizzato l'operatore collect, in modo da poter capire cosa sta succedendo.

Prenderò questo canale, rimuoverò il collect e dot view greetings, e poi duplicherò questa riga, aggiungerò l'operatore collect. E cambierò quello in after.

Questo è separato da dove stiamo chiamando questo, ma va bene perché stiamo usando le stesse chiamate di operatori sullo stesso canale di output.

Okay, clicchiamo su salva e proviamolo nel terminale. Farò lgoing nextflow run. Hello, workflow. Rieseguiamo il nostro script.

Okay. Questo sembra meglio. Come prima possiamo vedere i primi due processi eseguiti tre volte e ora il nostro processo finale è stato eseguito solo una volta.

Se guardiamo cosa è stato stampato dall'operatore view, qui sotto, abbiamo detto before collect, che è questo output qui, e che viene stampato tre volte. E potete vedere che c'è un singolo percorso per ognuno di questi. E poi after collect, potete vedere che abbiamo questo array di tre percorsi. Quindi è come ci aspettavamo.

Okay, controlliamo il file results e vediamo se è quello che ci aspettiamo questa volta. Infatti, ora ci sono tre righe nel file - ha concatenato con successo questi tre output in un singolo file di output. Fantastico.

Okay, farò pulizia e passiamo al prossimo passaggio. E cancellerò queste istruzioni view solo per mantenere le cose pulite.

## 3. Passare più di un input a un processo per nominare il file di output finale in modo univoco

Okay. Finora, tutti i nostri processi hanno preso solo un singolo input. Ora faremo un esercizio in cui aggiungiamo più di un input a un processo per vedere come funziona. Per farlo, useremo questo esempio collect greetings.

Ogni volta che eseguivo il workflow, sovrascriveva quel file nella directory results, che potrebbe non essere quello che vogliamo.

## 3.1. Modificare il processo di raccolta per accettare un nome definito dall'utente per il file di output

Quindi per questo esempio, passeremo un parametro aggiuntivo in modo da poter personalizzare il nome del file di output.

Aggiungere un secondo input a un processo è molto semplice. Aggiungo semplicemente una seconda riga nel blocco input. Questa volta sarà un valore, piuttosto che un percorso, perché vogliamo passare una stringa e lo chiamerò batch underscore name.

Posso ora usare questa variabile nel blocco script, e dirò collected dash dollar batch name.

Sto usando parentesi graffe qui intorno al nome della variabile. Questo è solo per tenerlo separato dal resto di una stringa, e probabilmente non è necessario in questo caso, ma penso che renda più facile la lettura.

Okay. Infine, ricordate di aggiornare il percorso di output perché ora il nome del file è cambiato, quindi farò la stessa cosa e metterò il batch name nell'output del percorso come previsto.

## 3.2. Aggiungere un parametro batch da riga di comando

Ora dobbiamo passare un batch name da qualche parte, e creerò un secondo parametro per farlo in modo da poterlo fare sulla riga di comando quando eseguiamo il workflow.

Quindi farò params batch name, e per impostazione predefinita, chiamiamolo test batch. Ora posso usare questa variabile di parametri speciale in basso, dove chiamiamo il processo.

E infatti VS Code ci sta dicendo che non ci sono abbastanza argomenti per questo processo ora, e che si aspetta un secondo input.

Semplicemente faccio virgola e passo la nostra nuova variabile e l'errore scompare.

Notate che l'ordine degli input qui è davvero importante. Il primo input del processo era il percorso, e il secondo input è il nome. Se cambio l'ordine qui, devo anche cambiare l'ordine quando chiamo il processo. Altrimenti Nextflow passerà il canale sbagliato all'input sbagliato.

## 3.3. Eseguire il workflow

Okay, proviamo e vediamo se funziona. Facciamo "nextflow run hello-workflow. Okay, è stato eseguito come prima. Diamo un'occhiata nella directory results.

Infatti, il nostro nome di file qui ora si chiama "collected test batch output txt". Fantastico.

E ora vediamo se possiamo sovrascriverlo eseguendolo di nuovo. Questa volta farò --batch_name per corrispondere a quel nome di variabile di parametro speciale qui. E lo chiamerò demo output.

Eseguiamo di nuovo il workflow e vedremo se succede qualcosa.

Okay, ora abbiamo un collected demo output .txt. E poiché questo nome di file è diverso da quello, non l'ha sovrascritto. Entrambi sono ora presenti nella directory results.

## 4. Aggiungere un output al passaggio di raccolta

Okay, quindi abbiamo mostrato come dare più input a un processo, ma per quanto riguarda più output? Per questo esempio, calcoleremo il numero di saluti che vengono elaborati e lo emetteremo come output secondario per questo passaggio collect greeting.

## 4.1. Modificare il processo per contare ed emettere il numero di saluti

Faremo un piccolo trucco qui. I processi Nextflow hanno questo blocco script con una stringa multi-linea, e quella viene passata come output bash al dot command dot sh. Ma possiamo effettivamente scrivere qualsiasi codice personalizzato sopra quello, e quello verrà eseguito come parte di un'attività ma non incluso all'interno dello script bash.

Una delle funzioni integrate nella sintassi Nextflow si chiama size. Quindi prenderò l'input del percorso, e dirò count underscore greetings, solo per definire un nome di variabile. Prenderò i file di input e chiamerò "size" su di esso.

Questa funzione conterà la dimensione di questo canale di input e lo assegnerà a una variabile.

Possiamo ora restituire quella variabile come parte del blocco output. Quindi diciamo, val, perché è un valore, non un file. E count greetings.

Ora questo è sufficiente da solo, e potremmo ora accedere a questi diversi output da questo processo. Tuttavia, dovremmo accedervi in modo posizionale. Quindi usando una chiave di indice come zero e uno.

Per rendere un po' più facile ottenere gli output, possiamo nominarli e lo facciamo utilizzando un'istruzione emit.

Quindi facciamo virgola emit out file o come voglio chiamarlo. E qui faccio emit count. Questo è fondamentalmente solo un decoratore, che ci aiuta a scrivere codice leggermente più pulito in modo da poter facilmente fare riferimento agli output specifici in seguito nel blocco workflow.

## 4.2. Riportare l'output alla fine del workflow

Okay. Se scorro verso il basso al blocco workflow, posso ora prendere gli output di collect greetings, fare collect greetings, dot out, e possiamo vedere i nostri due output denominati sono suggeriti qui dall'estensione VS Code. Molto utile.

Quindi farò dot count per ottenere il valore count che abbiamo appena creato, e farò view, in modo che venga stampato nella riga di comando. Così possiamo vederlo quando eseguiamo il workflow.

Scriviamo qualcosa nella closure qui solo per renderlo un po' più carino. num greetings, there were greetings greetings.

E in realtà non ci interessa l'altro output perché non lo stiamo usando come input per nessun altro processo. Ma potete vedere come potremmo facilmente passarlo come input a un altro processo se volessimo, a valle.

## 4.3. Eseguire il workflow

Cliccheremo su salva. Diamo un'occhiata al terminale e proviamolo.

Okay, fantastico. Eccoci qui. There are three greetings. È esattamente giusto.

Okay, ottimo lavoro. Questa è la fine di questo capitolo. Abbiamo finito per essere arrivati fin qui. Ora state iniziando a costruire un workflow piuttosto realistico, dove siamo in grado di gestire input e output e logica all'interno del nostro workflow.

Man mano che questi file di workflow diventano più lunghi, iniziano a diventare un po' ingombranti. Quindi nel prossimo capitolo, esamineremo come possiamo modularizzare il codice Nextflow in file separati in modo che sia più facile trovare e mantenere il codice all'interno del workflow.

Unitevi a noi nel prossimo video per il capitolo quattro. Hello Modules.

[Trascrizione video successivo :octicons-arrow-right-24:](04_hello_modules.md)
