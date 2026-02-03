# Parte 2: Hello Channels - Trascrizione

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduzione assistita da IA - [scopri di più e suggerisci miglioramenti](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/lJ41WMMm44M?si=xCItHLiOQWqoqBB9&amp;list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Note importanti"

    Questa pagina mostra solo la trascrizione. Per le istruzioni complete passo dopo passo, ritorni al [materiale del corso](../02_hello_channels.md).

    I numeri di sezione mostrati nella trascrizione sono forniti solo a scopo indicativo e potrebbero non includere tutti i numeri di sezione presenti nei materiali.

## Benvenuto

Salve, benvenuti alla parte due di Hello Nextflow.

Questo capitolo si chiama Hello Channels. Parleremo di questo aspetto fondamentale di Nextflow.

I channel sono gli elementi che collegano i diversi passaggi nella sua pipeline, il modo in cui i suoi dati e la logica fluiscono attraverso il suo workflow.

Bene, iniziamo.

Iniziamo andando su training.nextflow.io

Hello Nextflow nella barra laterale e facendo clic sulla parte due, Hello Channels.

Tutto il materiale è scritto qui, così può seguire al suo ritmo e recuperare qualsiasi cosa possa aver perso.

Una volta aperto il sito web, può caricare Codespaces e continueremo da dove eravamo rimasti alla fine del capitolo precedente.

## 0. Riscaldamento: Eseguire hello-channels.nf

Per questo capitolo, modificheremo un file diverso. Questo si chiama Hello Channels, quindi può trovarlo nella barra laterale, faccia doppio clic per aprirlo.

Ora, se proviene appena dal capitolo uno, questo file le sembrerà molto familiare. Il punto di partenza qui è sostanzialmente dove finiamo il capitolo uno, con il nostro process chiamato sayHello, il nostro input, output, il nostro publishDir e il nostro params.greeting, e il nostro workflow semplice.

Stiamo iniziando con un nuovo file, quindi è un terreno di gioco livellato per tutti, ma può continuare con il suo file precedente se preferisce.

Nota, ho anche eliminato tutti i file .nextflow\* e le directory work qui, solo per avere un punto di partenza pulito. Non importa se lo fa o no, dipende da Lei.

Bene. Iniziamo verificando che questa pipeline funzioni ancora come ci aspettiamo. Porterò su il terminale qui.

Faccio "nextflow run hello-channels.nf" e premo invio.

Eseguirà quel piccolo workflow, esegue il nostro passaggio sayHello, genera una directory work con quell'hash, ed ecco la nostra cartella results e c'è il nostro file di output, proprio come ci aspettavamo dal nostro params.greeting predefinito.

Quindi è ottimo. Esattamente come nel capitolo uno, funziona come ci aspettiamo.

## 1. Fornire input variabili tramite un channel in modo esplicito

Nel capitolo uno, stava già utilizzando i channel, solo che non se ne rendeva conto. Quando abbiamo specificato una stringa qui, Nextflow ha automaticamente creato un channel attorno a quella stringa per noi, proprio perché sapeva che stavamo chiamando un process, quindi avevamo bisogno di un channel di input.

La prima cosa che faremo è renderla esplicita scrivendo effettivamente il channel stesso.

## 1.1. Creare un channel di input

Quindi andrò al workflow qui in fondo allo script, e dirò greeting_ch. Questa è una convenzione che spesso usiamo nel codice Nextflow: avere un underscore ch alla fine di un nome di variabile quando è un channel, solo per rendere facile identificare che è un channel, ma non è obbligatorio farlo. Uguale channel of Hello Channels.

Quello che abbiamo appena usato è qualcosa chiamato "Channel Factory" nel linguaggio Nextflow. È questa cosa qui, stiamo impostando questa variabile su un nuovo channel, e questa channel factory qui sta creando un channel per noi in un modo particolare.

Ci sono una manciata di diverse channel factory che Nextflow ha, per creare channel da diversi tipi di input. Dot of è la più semplice, e accetta semplicemente qualsiasi stringa che le forniamo.

Noti che quando passo il mouse su queste parole in VS Code, l'estensione Nextflow mi sta dando un popup che spiega cosa fa questa sintassi, e c'è anche un testo "leggi di più" in fondo a quella finestra popup.

Se faccio clic su quello, aprirà la documentazione Nextflow in una nuova scheda e mi porterà direttamente alla documentazione per questa cosa specifica. In questo caso per channel.of.

## 1.2. Aggiungere il channel come input alla chiamata del process

Noti che l'estensione ci sta anche dando un avviso, dicendo che abbiamo creato un nuovo channel qui, ma non viene utilizzato da nulla.

Quindi, sistemiamolo. Prenderò il nuovo nome del channel e sostituirò questo params.greeting con il nostro nuovo channel.

Noti che non stiamo più usando il flag della riga di comando --greeting ora, params.greeting non viene utilizzato, stiamo tornando a codificare questa stringa in modo fisso. Va bene. Sto solo cercando di mantenere le cose semplici. Torneremo più tardi e useremo di nuovo i params.

## 1.3. Eseguire nuovamente il comando workflow

Bene, verifichiamo solo che questo funzioni. Portiamo su il terminale e notiamo di nuovo. Nextflow run hello channels. Verifico output.txt, ed eccolo.

Ottimo, un esempio un po' noioso, che fa esattamente la stessa cosa di prima, ma ora almeno la logica è un po' più chiara. Stiamo essendo espliciti nello scrivere un nuovo channel.

Abbiamo effettivamente appena scritto più codice per fare la stessa cosa. Ma questo inizierà ad avere più senso man mano che diventiamo un po' più complicati nel modo in cui creiamo i nostri channel.

## 2. Modificare il workflow per eseguirlo su valori di input multipli

Bene, rendiamo questo un po' più interessante. È molto raro che si voglia eseguire una pipeline Nextflow su un singolo input, quindi diamogli diversi input.

## 2.1. Caricare più saluti nel channel di input

Dalla documentazione qui. Copierò questi diversi stringhe, tre di esse. Hello, Bonjour, Olà. Oh, ottieni Speranza. Copilot sta suggerendo un paio di altri. Quindi premiamo tab e inseriamo quelli.

La documentazione Nextflow qui ci dice che possiamo dare valori multipli a questo operatore, quindi dovrebbe funzionare, ma proviamolo e vediamo cosa succede.

## 2.1.2. Eseguire il comando e guardare l'output del log

Bene. Sì e no. Vediamo. Dice che cinque di cinque attività sono state eseguite qui, ma ci mostra solo un hash, il che è un po' strano. Va bene. Tutto è come previsto qui. Per impostazione predefinita. Nextflow utilizza un tipo speciale di output per un terminale chiamato codici di controllo ANSI, il che significa che sovrascrive determinate righe per dare una bella visualizzazione compressa di tutti i diversi processi che vengono eseguiti.

Questo ha molto più senso quando si hanno workflow più grandi e si stanno eseguendo centinaia o migliaia di campioni diversi. Si può semplicemente generare così tanto output sul terminale che è impossibile guardarlo, mentre questa visualizzazione di aggiornamento le dà un progresso in tempo reale.

## 2.1.3. Eseguire nuovamente il comando con l'opzione -ansi-log false

Se vuole, può eseguirlo di nuovo, e questa volta userò un argomento core Nextflow aggiuntivo con un singolo trattino che dice, "-ansi-log false". Questo utilizza la versione precedente dell'output del log Nextflow. E qui può vedere tutti i processi individuali che sono stati lanciati.

Sta a Lei decidere se farlo o no. L'output da Nextflow è esattamente lo stesso in entrambi i casi.

## 2.2. Assicurarsi che i nomi dei file di output siano univoci

Bene, diamo un'occhiata ai file di output, poi andremo ai results. Ma c'è solo un singolo file di output. Cosa è successo? Abbiamo visto che il process era stato eseguito molte volte. Possiamo andare nella directory work e vedere tutti i diversi hash, tutte le attività sono state eseguite correttamente. Ma se ricorda nel nostro process qui, stiamo salvando tutto in un file output.txt e poi pubblicando quello in questa directory.

Quindi lo stesso file è stato creato cinque volte, e poi è stato sovrascritto cinque volte. E abbiamo solo qualsiasi attività sia stata eseguita per ultima.

## 2.2.1. Costruire un nome di file di output dinamico

Il modo in cui risolviamo questo è utilizzando un nome di file di output dinamico. Qui abbiamo già una variabile chiamata greeting all'interno del process, quindi possiamo usarla nel nome del file di output. La copio e faccio $greeting-output.txt.

Circonderò questo tra virgolette, solo per evitare che bash si confonda con eventuali spazi che potrebbero insinuarsi qui. E poi prenderò lo stesso nome di file e aggiornerò l'output qui.

È davvero importante che l'output corrisponda a questo, perché altrimenti, questo file non verrà trovato e Nextflow andrà in crash.

Farò un'altra modifica davvero importante, che è cambierò questi apici singoli per apici doppi. Noti che il colore del codice è cambiato quando l'ho fatto. Questa variabile viene espansa solo se usiamo apici doppi. Se uso apici singoli qui, viene usata come valore letterale, e otterrei un singolo file chiamato $greeting-output, che non è quello che voglio.

## 2.2.2. Eseguire il workflow

Quindi rimettiamo gli apici doppi e proviamo.

Pulirò solo la mia directory prima di iniziare, così è facile vedere i nuovi file. Eliminerò tutto ciò che si chiama .nextflow, work e results.

E eseguirò di nuovo quel comando Nextflow e vediamo quali file vengono creati. Quindi esegue i cinque processi là. Se stava guardando molto attentamente, potrebbe aver visto quella riga aggiornarsi mentre era in esecuzione.

E ora possiamo andare nella directory results, e infatti, abbiamo cinque output diversi, e sono tutti prefissati con il saluto diverso.

Se apro ognuno di questi, vedremo che ognuno contiene il saluto corrispondente. Fantastico. È quello che vogliamo.

## 3. Utilizzare un operatore per trasformare il contenuto di un channel

Bene, quindi ora sappiamo cosa sono i channel e sappiamo cosa sono le channel factory. E gli operatori? Questo è un altro termine per parte del linguaggio Nextflow, che è una serie di funzioni che ci permettono di operare sui channel per fare certe cose con loro. Nextflow, viene fornito con una suite di operatori, che ci permettono di manipolare i channel in una varietà di modi diversi.

## 3.1. Fornire un array di valori come input al channel

Lavoriamo attraverso questo con un esempio. Diciamo che vogliamo prendere queste stringhe di input, ma invece di metterle direttamente in una channel factory, vogliamo definirle come un array.

## 3.1.1. Impostare la variabile di input

Quindi prenderò queste e lo farò come una nuova riga sopra e dirò, greetings, array.

Ecco qua. Prenderò quella variabile array e la metterò nel channel.of, e premo salva.

## 3.1.3. Eseguire il workflow

Ora, vediamo cosa succede. Torno al mio terminale. Pulirò di nuovo tutti quei file temporanei. E eseguiamo il workflow.

Non bene. Bene. Si è rotto. Va bene. Me lo aspettavo questa volta. Il debug di cosa va storto quando un workflow Nextflow fallisce è una parte fondamentale dell'essere uno sviluppatore Nextflow. Questo accadrà molto ed è importante capire cosa dice il messaggio di errore e come gestirlo.

I messaggi di errore Nextflow sono in realtà abbastanza strutturati. Ci dice quale process è andato storto. Ci dà un messaggio di errore per una ragione. Dice quale era il comando che ha cercato di eseguire all'interno di quella particolare attività, quale era lo stato di uscita, quale era l'output su dove si trovava quella directory work dell'attività.

Noti che posso option, fare clic su questo in VS Code e lo apre in una barra laterale così posso andare direttamente lì e visualizzare tutti questi file nascosti, di cui abbiamo parlato nel capitolo precedente, incluso il file .command.sh. Questo può vedere che è lo stesso dei comandi che è stato eseguito qui.

Guardando questo file, possiamo avere un'idea di cosa potrebbe essere andato storto qui invece di eseguire una singola attività per ciascun elemento nell'array come ha fatto l'ultima volta, ha semplicemente fornito l'intero array in una volta come una stringa. Quindi dobbiamo decomprimere quell'array in valori individuali prima di passarlo nel channel. Torniamo indietro e vediamo se possiamo farlo usando un operatore.

## 3.2. Utilizzare un operatore per trasformare il contenuto del channel

In questo caso, non cambieremo l'array prima di passarlo nel channel. Adatteremo il channel in modo che si comporti nel modo che ci aspettiamo. Lo faremo usando l'operatore flatten può fare dot iniziare a digitare e possiamo vedere che l'estensione VS Code inizia a suggerire tutti i diversi operatori che abbiamo disponibili.

## 3.2.1. Aggiungere l'operatore flatten()

E selezionerò flatten. Noti che lo spazio bianco non importa in questo contesto per Nextflow. Quindi può mettere questi operatori su una nuova riga se vuole. Quindi posso far cadere questo qui e indentarlo in modo che si trovi sotto ".of" e vedrà che le persone spesso concatenano molti operatori come questo su un channel e lo indentano in questo modo in modo che sia più facile da leggere.

Può anche vedere, come prima posso passare il mouse su questo e leggere cosa sta facendo l'operatore flatten, e seguire anche un link alla documentazione se voglio.

Quindi questo operatore sta prendendo questo channel, che ha un singolo array al suo interno, e separando i valori dell'array.

## 3.2.2. Aggiungere view() per ispezionare il contenuto del channel

Possiamo sbirciare nei channel utilizzando lo speciale operatore view, e ne aggiungerò un paio qui. Questo è un po' come usare istruzioni di stampa in altri linguaggi. Quindi farò dot view e poi userò queste parentesi graffe.

Questo si chiama closure. Questo fondamentalmente fornisce codice aggiuntivo all'operatore view, che eseguirà su ciascun elemento all'interno del channel. In questo caso, dirò greeting before flatten. Greeting.

Sto definendo una variabile qui, che è solo all'interno dell'ambito di questa closure. Quindi questa variabile viene utilizzata solo qui e potrei chiamarla come voglio. Non importa davvero. Sto solo usando greeting per renderla facile da leggere.

In alcune pipeline Nextflow, potrebbe vedere persone usare una variabile implicita speciale chiamata "$it". Così. Questa è una variabile speciale all'interno del codice Nextflow, che è una scorciatoia in modo da non dover fare la piccola definizione di una variabile. Tuttavia, nel tempo stiamo pensando, questo non è molto chiaro per le persone che sono nuove a Nextflow, e scoraggiamo l'uso di "$it" ora.

Quindi mi atterrò al comportamento precedente di greeting e lo userò così perché è più esplicito ed è più chiaro su cosa sta succedendo.

Copierò quindi questa riga e farò esattamente la stessa cosa di nuovo dopo gli argomenti flatten. L'operatore view è un po' speciale perché fa qualcosa sugli elementi, ma continua anche a passarli all'operatore successivo così possiamo concatenarlo nel mezzo di una catena di operazioni come questa, e stamperà lo stato lì e continuerà. Quindi speriamo che questo ci mostri come appare il channel prima e dopo l'operatore flatten.

## 3.2.3. Eseguire il workflow

Proviamolo. Pulisci. Pulisci tutto nello spazio di lavoro. Esegui di nuovo la pipeline.

Bene, quindi possiamo vedere che ha eseguito i nostri cinque processi. Di nuovo, non si è schiantato con un errore, quindi è decisamente buono. E ora abbiamo il before flatten e sicuramente abbiamo il nostro array e abbiamo after flatten, stampato cinque volte una volta per ciascun elemento dell'array. È esattamente quello che speravamo. Quindi è davvero una buona notizia. E questo si adatta esattamente a quello che ci aspetteremmo dal codice.

Non abbiamo più bisogno di queste istruzioni di debug, quindi posso commentarle o eliminarle. Le eliminerò solo per mantenere il mio codice bello e pulito. Bene, ottimo. Questo esempio ora funziona bene e possiamo iniziare a vedere come i channel possono fare una logica un po' più complicata.

## 4. Utilizzare un operatore per analizzare i valori di input da un file CSV

Ora proveremo a farlo usando un file con una serie di input invece. Questo è un modo molto comune per scrivere pipeline Nextflow usando un foglio campione o un CSV di metadati.

## 4.1. Modificare lo script per aspettarsi un file CSV come fonte di saluti

Se vado sulla barra laterale, può vedere greetings.csv nel repository di esempio, e questo è un file CSV molto, molto semplice che contiene solo tre righe con tre saluti diversi. Vediamo se possiamo usare questo file CSV all'interno del nostro workflow.

Ora tornerò a usare params come abbiamo fatto nel capitolo uno, in modo da poter avere un input da riga di comando.

Eliminerò questo array greetings.

## 4.1.1. Cambiare il parametro di input per puntare al file CSV

Imposterò params greeting al nome del file, che è greetings.csv, e userò questa variabile speciale per generare il channel. Metterò quello lì, e gli errori scompaiono. Ricorda che questo sta impostando questa variabile per impostazione predefinita ora. Quindi se eseguo la pipeline senza alcun argomento, userà greetings.csv, ma potrei fare --greeting per sovrascrivere questa variabile se volessi.

## 4.1.2. Passare a una channel factory progettata per gestire un file

Bene, stiamo passando un file ora piuttosto che una stringa o un array di stringhe, quindi probabilmente abbiamo bisogno di una channel factory diversa.

Ci libereremo di "of" che abbiamo usato finora, e useremo invece .fromPath. Questo fa esattamente quello che sembra. Crea un channel con percorsi invece di valori, usando un nome di file stringa o glob. Rimuoverò anche l'operatore flatten poiché non ne abbiamo più bisogno, ora che stiamo passando un file.

## 4.1.3. Eseguire il workflow

Premerò salva, aprirò il terminale, eseguirò il workflow, e poi vedremo cosa succede.

Bene. Si è schiantato di nuovo. Non si preoccupi. Me lo aspettavo anche questo. Diamo un'occhiata al messaggio di errore e vediamo se possiamo capire cosa sta andando storto. Qui possiamo vedere il comando eseguito, e un po' come prima dove avevamo l'intero array stampato. Ora abbiamo il percorso del file che viene fatto echo nel comando, piuttosto che passare attraverso il contenuto del file.

## 4.2. Utilizzare l'operatore splitCsv() per analizzare il file

Quindi per usare il contenuto del file invece, abbiamo bisogno di un altro operatore. L'operatore che useremo per questo è chiamato splitCsv. Ha senso, perché è un file CSV che stiamo caricando.

## 4.2.1. Applicare splitCsv() al channel

Ok, quindi splitCsv. Chiudi parentesi. Non abbiamo bisogno di alcun argomento qui. E di nuovo, userò alcuni operatori view per dare un'idea di cosa sta succedendo qui.

.view csv after splitCsv. Before split Csv.

## 4.2.2. Eseguire nuovamente il workflow

Bene, proviamo a eseguire questo e vedere cosa succede.

Bene, abbiamo un po' più di output questa volta, ma è ancora fallito. Possiamo guardare le istruzioni view, e qui può vedere before split CSV, e abbiamo un percorso di file come abbiamo visto nel messaggio di errore precedente. After split CSV, ora abbiamo tre valori corrispondenti alle tre righe nel file CSV.

Tuttavia, può vedere che ognuno di questi valori è circondato da parentesi quadre. Quindi ognuno di quelli era un array a sé stante, e questo ci ha dato la stessa area che avevamo prima dove sta cercando di fare echo di un array piuttosto che solo una singola stringa.

Se pensiamo a un file CSV, questo ha un certo senso. Tipicamente, un file CSV avrà righe e colonne, quindi split CSV fa un array bidimensionale. La prima dimensione dell'array è ogni riga, e poi c'è una seconda dimensione, che è ogni colonna per ogni riga.

Quindi qui abbiamo solo un singolo valore su ogni riga, quindi abbiamo una singola colonna, quindi abbiamo un array di un elemento per ogni riga del file.

Va bene. Abbiamo solo bisogno di un altro operatore per collassare quell'array per ogni riga del file CSV analizzato. Puliamo questo. Sbarazzarsi di un terminale e vedere cosa possiamo fare.

## 4.3. Utilizzare l'operatore map() per estrarre i saluti

Ora potremmo usare di nuovo l'operatore flatten, che abbiamo usato prima. Abbiamo visto come quello può collassare un array in una serie di valori, che funzionerebbe molto bene qui. Ma userò l'opportunità per dimostrare un altro operatore, che è molto comune all'interno dei workflow chiamato operatore map.

## 4.3.1. Applicare map() al channel

Farò dot map e farò item item[0].

Se scrive molto altro codice in linguaggi, potrebbe già avere familiarità con l'operatore map. Prende un iterabile, come un array o un channel, e fa un'operazione su ciascun valore di quello.

Qui stiamo dicendo che dovremmo definire una variabile chiamata item all'interno dell'ambito di questa closure, e poi vogliamo restituire, solo il primo valore in quell'array. Quindi item indice zero.

Questo sta effettivamente appiattendo l'array. Può vedere come potremmo estendere questo per essere più complesso, però: se il nostro file CSV avesse sei colonne, ma siamo interessati solo alla quarta colonna, potremmo accedere a un indice specifico qui. O fare qualsiasi altro tipo di operazione sul valore prima di passarlo all'elaborazione a valle.

Quindi l'operatore map è estremamente flessibile e molto potente per modificare i channel in volo. Mettiamo un'altra istruzione view solo per poter vedere cosa sta facendo nella nostra esecuzione. Può giudicare quella riga e spostarla verso il basso. E after map.

## 4.3.2. Eseguire il workflow un'ultima volta

Portiamo su il terminale e proviamo a eseguire il workflow.

Bene, nessun errore questa volta. È un buon segno. Ora possiamo passare attraverso tutti questi diversi output dalle istruzioni view. Before split CSV, avevamo un singolo percorso. After split CSV, avevamo gli array di valore singolo, e poi after map, abbiamo solo i valori senza alcuna sintassi di array. Andiamo alla directory results, ed ecco i nostri file di output che si comportano esattamente come volevamo.

C'è un piccolo bonus qui. Può effettivamente vedere che gli operatori view sono leggermente mescolati nell'ordine in cui hanno fatto l'output. Questo perché Nextflow sta facendo la parallelizzazione di queste diverse attività. Quindi dopo aver diviso il CSV, ci sono tre elementi in questo channel, e sta gestendo l'elaborazione di quei tre elementi in parallelo automaticamente. Ciò significa che l'ordine degli output è stocastico e può variare. In questo caso, è semplicemente successo che alcuni degli operatori view sono stati restituiti dopo che il passaggio successivo era stato completato, e quindi è arrivato in questo ordine.

Se eseguo di nuovo lo stesso workflow. Allora infatti, è arrivato in un ordine diverso e questa volta abbiamo gli split CSV e le map nell'ordine che ci aspetteremmo.

Quindi tenga solo presente, non può fare affidamento sull'ordine degli output da un'attività del process perché Nextflow sta gestendo questa parallelizzazione per Lei automaticamente. Nextflow lo fa per Lei con la sua logica di flusso dati, e questo è il vero potere di Nextflow.

Bene, questo è probabilmente uno dei capitoli più importanti di tutta la formazione. Una volta compresi i channel, le channel factory e gli operatori, inizia a entrare nella forza di Nextflow e ciò che lo rende unico come linguaggio di programmazione. Questa funzionalità consente a Nextflow di parallelizzare tutti i suoi workflow per Lei e generare una logica di workflow estremamente complessa con una sintassi molto pulita e un modello di flusso di dati push. Può essere un concetto un po' strano all'inizio, ma una volta che si abitua a scrivere codice come questo, si sentirà rapidamente naturale e prima che se ne renda conto, scriverà fantastici workflow.

Faccia una pausa, una tazza di tè, cammini un po' e passiamo al capitolo tre, dove iniziamo ad estendere questi concetti in workflow più complessi. Ci vediamo nel prossimo video.

[Trascrizione video successivo :octicons-arrow-right-24:](03_hello_workflow.md)
