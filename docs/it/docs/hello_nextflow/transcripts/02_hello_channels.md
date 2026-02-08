# Parte 2: Hello Channels - Trascrizione Video

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduzione assistita da IA - [maggiori informazioni e suggerimenti per migliorare](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/yDR66fzAMOg?si=xCItHLiOQWqoqBB9&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Note importanti"

    Questa pagina mostra solo la trascrizione. Per le istruzioni complete passo dopo passo, tornate al [materiale del corso](../02_hello_channels.md).

    I numeri di sezione mostrati nella trascrizione sono forniti solo a scopo indicativo e potrebbero non includere tutti i numeri di sezione presenti nei materiali.

## Benvenuti

Ciao e bentornati alla Parte 2 di Hello Nextflow. Questo capitolo si chiama Hello Channels.

I canali sono come la colla nella vostra pipeline Nextflow. Sono i pezzi che tengono insieme tutti i diversi processi, che Nextflow usa per passare tutte le informazioni e orchestrare il vostro flusso di lavoro.

C'è un'altra parte relativa ai canali che sono gli operatori. Questi sono fondamentalmente funzioni che possiamo usare sui canali per modificarne i contenuti. Addentriamoci in VS Code e vediamo dove siamo.

Sono molto ingrandito su questo VS Code, quindi per mantenere le cose pulite e ordinate, ho rimosso tutti i file _.nextflow\*_ e la directory _work/_ e _results/_ e tutto dal Capitolo Uno. Sto solo ricominciando da capo qui. Ma non preoccupatevi troppo di questo. Se non volete, potete lasciare quei file lì. Non causeranno problemi.

Inizieremo lavorando su _hello-channels.nf_ per questo capitolo, e se lo apro, dovrebbe sembrare molto simile al file su cui stavamo lavorando in precedenza. Potrebbero esserci parti diverse in sezioni diverse dello script, ma tutto dovrebbe essere fondamentalmente lo stesso.

Una cosa che è diversa è che il path nel blocco output qui è ora _hello_channels_ per questa parte, il che significa che i file dei risultati saranno memorizzati in una sottodirectory diversa nei vostri risultati se avete ancora quella lì. Quindi dovrebbe essere un posto bello e pulito per iniziare senza confondersi con gli output.

Okay, quindi ricordiamo velocemente cosa fa questo script quando eseguiamo questo flusso di lavoro. Facciamo _"nextflow run hello-channels.nf"_. Possiamo fare _"--input myinput"_, e quando eseguiamo questo, userà questo parametro, params.input, che è stato passato come variabile per il processo sayHello qui sopra, che va in greeting e viene salvato in output.txt. E possiamo vedere questo nel file dei risultati. Fantastico.

## 1. Fornire input variabili tramite un canale esplicitamente

Va bene. Ma è, è abbastanza semplicistico. Abbiamo una variabile in questo parametro, che va in un processo che viene eseguito una volta, e non scala davvero. E non possiamo dargli molti file diversi da creare qui. Non possiamo dargli molti saluti diversi. Ne abbiamo solo uno.

In realtà, Nextflow riguarda tutto lo scaling della vostra analisi. Quindi probabilmente volete che faccia più di una cosa. E lo facciamo con i _canali_.

I canali sono un concetto un po' unico per molte persone che iniziano con Nextflow. Viene da questo tipo di concetti di programmazione funzionale, e può richiedere un po' di tempo per capirlo, ma una volta che scatta, sbloccano davvero il potere di Nextflow ed è la chiave per come scrivete i vostri flussi di lavoro.

## 1.1. Creare un canale di input

Iniziamo prendendo questo script e facendolo usare un _canale_ invece di solo un _param_.

Andiamo giù al workflow, che è dove tutta la nostra logica del flusso di lavoro riguarda il collegare le cose insieme. E sto per andare qui dentro e creerò un nuovo canale.

Creare un nuovo canale.

E lo chiamerò "_greeting_ch"_. Questa è una convenzione per fare "_\_ch"_ in questo modo, solo così potete ricordare che questa variabile è un canale. Ma potete chiamarlo come volete.

E poi dirò uguale, e farò _"channel.of"._

Channel è come lo spazio dei nomi per tutto ciò che riguarda i canali. "c" minuscola se avete usato Nextflow prima. E il _".of"_ è qualcosa chiamato fabbrica di canali, che è fondamentalmente un modo per creare un canale.

Ci sono molte fabbriche di canali diverse. Se faccio solo "." qui, potete vedere che VS Code ne suggerisce un sacco, ma _".of"_ è il più semplice e prende solo un input qui.

Quindi posso fare delle parentesi e dirò _"Hello Channels!"_.

Fantastico. Ho un canale. Fantastico. Posso salvare, potrei eseguirlo di nuovo, ma non succederà nulla di interessante. VS Code mi ha dato una linea di avviso arancione qui sotto e mi ha detto che questo è impostato: l'hai creato, ma in realtà non l'hai mai usato per nulla. Questo canale non viene consumato.

Okay, quindi come lo usiamo? Molto semplice. Prenderò questo, lo copierò, ed eliminerò _params.input_ e metterò _"greeting_ch"_ qui invece. Quindi passeremo questo canale come input a sayHello.

Notate che ho codificato questa stringa per ora. Questo è un po' un passo indietro dopo il nostro bel param che abbiamo usato alla fine del capitolo precedente, ma mantiene le cose semplici per iniziare così potete vedere la logica.

Okay, vado nel mio terminale ed eseguirò di nuovo il flusso di lavoro. Senza _"--input"_ questa volta, e verrà eseguito e userà quel canale che abbiamo creato e si spera dovremmo avere un file qui in _results/hello_channels/_ e ora dice "Hello Channels!". Fantastico. Quindi è quello che speriamo, dal nostro canale qui. Fantastico.

## 1.4. Usare view() per ispezionare i contenuti del canale

Un'altra cosa da aggiungere qui, solo una rapida introduzione a un'altra funzione che possiamo usare sui canali chiamata "_.view"_.

Questa è analoga al comando _print_ in Python o altri linguaggi a cui potreste essere abituati, e scarica semplicemente i contenuti di questo canale sul terminale quando lo eseguiamo.

Quindi faccio "_.view"_, e poi se rieseguo di nuovo il flusso di lavoro, dovrebbe stampare sul terminale quali sono i contenuti di quel canale, al momento in cui l'abbiamo creato.

Infatti, potete vedere che è stampato sul terminale qui. _"Hello Channels!"_.

Notate che potete spezzare queste cose su più righe se volete, e infatti, il formattatore automatico di Nextflow cercherà di farlo per voi. Lo spazio bianco non è davvero importante qui, quindi potete concatenare queste cose una dopo l'altra.

## 2. Modificare il flusso di lavoro per eseguire su più valori di input

Okay, quindi il nostro canale ha una cosa dentro che è bella, ma è fondamentalmente la stessa di prima. Quindi rendiamola un po' più complicata. Aggiungiamo qualche altra cosa al nostro canale.

La fabbrica di canali "_.of()"_ può prendere più elementi, quindi scriviamone qualcun altro. Faremo _Hello, Bonjour, Hej_. E poi possiamo eseguire di nuovo questo flusso di lavoro e vedremo cosa succede.

Dovrebbe eseguire di nuovo. E ora abbiamo stampato. _"Hello", "Bonjour"_ e _"Hej"_ sul terminale con la nostra dichiarazione view. Fantastico.

## 2.1.2. Eseguire il comando e guardare l'output del log

Potreste pensare che abbiamo finito a questo punto. Ma in realtà c'è un po' un trucco qui, che ci farà inciampare. Se guardiamo il nostro file di output qui. Potete vedere che ha _"Hello"_ dentro, ma non ha nessuno degli altri output. Infatti, è solo questo.

Se eseguiamo questo flusso di lavoro più volte, potremmo anche vedere che a volte ha _"Bonjour"_, a volte ha _"Hej"_. È un po' casuale.

Se guardiamo il terminale, possiamo vedere che è stato eseguito tre volte e possiamo vedere i diversi output di view. Ma se vado nella directory di lavoro, posso fare _"cat work"_. Metto questo hash ed espando quello e _output.txt_. Potete vedere che questo file nella directory di lavoro è diverso dalla directory dei risultati, e questo è _"Hej"._ Quindi c'è qualcosa che non funziona bene qui.

E la chiave è che, abbiamo tre attività che sono state eseguite. L'output di Nextflow cerca di riassumere questo mentre l'elaborazione procede, così non prende completamente il controllo di tutto il vostro terminale, e quel Logging ANSI usa codici di escape ANSI, ha fondamentalmente sovrascritto le altre attività. Quindi vi mostra solo l'ultima che è stata aggiornata.

## 2.1.3. Eseguire di nuovo il comando con l'opzione -ansi-log false

Ci sono alcune cose che possiamo fare per capire effettivamente questo un po' meglio. Possiamo guardare nella directory di lavoro stessa e potete vedere tutte le diverse directory di lavoro lì, ma è un po' confuso perché saranno mescolate con diverse esecuzioni di Nextflow.

Oppure possiamo dire a Nextflow di non usare i codici di escape ANSI.

Quindi se eseguo di nuovo il comando, ma questa volta dico _"-ansi-log false"_ per disattivarlo, potrei anche usare le variabili d'ambiente _$NO_COLOR_ o _"$NXF_ANSI_LOG=false"_. Quindi usa il tipo di stile più vecchio di logging di Nextflow senza nessuno di questi codici di escape. Stampa direttamente su un terminale senza aggiornamenti intelligenti.

E ora possiamo vedere tutti e tre questi processi che sono stati eseguiti. E ognuno di loro ha il proprio hash di attività. E se andiamo in queste directory di lavoro, vedremo i tre diversi saluti che abbiamo specificato.

Quindi ha un po' più senso ora. Spero che capiate che Nextflow stava facendo questo, stava solo essendo un po' intelligente con quello che vi mostrava nel terminale con quelle directory di lavoro.

Tuttavia, questo è risolto per un problema con le directory di lavoro, ma non ha risolto un problema con il file di output. Abbiamo ancora solo un file di output che dice _"Hello"_.

## 2.2. Assicurarsi che i nomi dei file di output siano unici

Ora per capire questo, dobbiamo tornare al nostro script del flusso di lavoro. Stiamo generando il nostro canale qui, lo stiamo passando al nostro processo, e se guardiamo il processo, stiamo scrivendo il saluto in un file chiamato _"output.txt"_ e passando quel file di output indietro al blocco output qui sotto, pubblicandolo.

Tuttavia, ogni tre volte che questo processo viene eseguito queste tre diverse attività. Generano tutte un file chiamato _"output.txt"_, tutti quei file di output vengono pubblicati nella directory dei risultati, e si sovrascrivono tutti a vicenda. Quindi qualunque file di risultato otteniate lì è solo l'ultimo che è stato generato, ma ha cancellato tutti gli altri. Non è davvero quello che vogliamo.

## 2.2.1. Costruire un nome di file di output dinamico

Ci sono diversi modi per gestire questo, ma il più semplice per ora è solo creare nomi di file univoci diversi. Quindi ogni volta che l'attività viene eseguita con un saluto diverso, genererà un file di output diverso, che non entrerà più in conflitto quando pubblicato. E poi avremo tre file di output univoci.

Facciamo questo esattamente nello stesso modo. Possiamo usare questa variabile ovunque all'interno del blocco script e possiamo usarla più volte.

Quindi posso incollarlo qui, _"$\{greeting\}\_output.txt"_, e poi devo anche incollarlo qui sopra perché non stiamo più creando un file chiamato _output.txt_. Quindi se non aggiorno questo, Nextflow andrà in crash con un errore dicendo che si aspettava un file, che non è mai stato generato.

Quindi devo fare lo stesso lì e devo usare virgolette doppie, non virgolette singole, così questa variabile viene compresa.

Okay, proviamolo e vediamo se ha funzionato. Eseguiremo di nuovo il flusso di lavoro. Si spera ci mostrerà le tre diverse attività all'interno delle tre diverse directory di lavoro. E infatti, potete vedere nella cartella dei risultati qui sopra a sinistra. Ora abbiamo tre file diversi con tre nomi di file diversi e ognuno con i contenuti diversi che ci aspettiamo. Quindi i file non si cancellano più a vicenda, e tutto è lì come ci aspettiamo.

Questa è un po' una configurazione banale che abbiamo attraversato qui, ma sottolinea alcuni dei concetti chiave che dovete capire su come funziona la pubblicazione dei file, e alcune delle cose in cui potreste cadere come trappole. Quindi spero che possiate evitare questo nei vostri flussi di lavoro.

Vale anche la pena notare che quello che abbiamo fatto qui è un po' impratico nelle situazioni di vita reale. Abbiamo preso alcuni dati di input e stiamo usando quei dati, ma stiamo anche nominando il file dopo quei dati, cosa che di solito non potete fare.

Quindi nelle pipeline Nextflow più mature e reali, passerete spesso in giro un oggetto meta con tutti i metadati associati a un dato campione. Potete quindi creare nomi di file dinamici basati su quello, il che è molto più pratico.

Se siete interessati a come fare questo con le best practice, c'è una side quest su _training.nextflow.io_, che riguarda specificamente i metadati e le meta map, quindi potete approfondire lì per maggiori dettagli.

## 3. Fornire più input tramite un array

Okay. Ora esploreremo un po' come sono strutturati i canali e come differiscono da altri tipi di strutture dati nel linguaggio di codifica. E penserò un po' a come potrei potenzialmente usare un array, che potrebbe essere un concetto familiare se venite da altri linguaggi.

Posso usare un array in un canale? Proviamolo. Creerò un array, e l'ho copiato dai documenti, _"greetings_array"_ e _"Hello", "Bonjour"_ e _"Holà"_. E poi lo metterò qui invece delle mie stringhe codificate. Quindi dirò "channel.of" _"greetings_array"_, passando questo array in un canale. Proviamolo.

Tiro su il terminale, ed eseguo la pipeline.

Okay. Potete vedere che la dichiarazione view qui ha stampato il nostro array come previsto, ma poi tutto questo testo rosso, o non sarà rosso se avete ancora _"-ansi-log"_ spento, ma tutto questo testo rosso ci sta dicendo che qualcosa è andato storto.

Non abbiamo più un bel segno di spunta verde qui. Abbiamo una croce rossa, e se rendo solo questo un po' più largo così è più facile da leggere, Nextflow ci sta dicendo cosa è andato storto.

Quindi scomponiamo questo sezione per sezione. Dice che l'errore è stato causato da, e poi la ragione dell'errore, che è file di output mancanti. Quindi fondamentalmente quel blocco output diceva che questo file doveva essere creato e non lo era. Poi dice che questo è il comando che è stato eseguito. Quindi questo è fondamentalmente il contenuto di quel file _.command.sh_. Questo è come appariva dopo che tutte quelle variabili erano state inserite.

E potete vedere qui il nostro comando echo è stato effettivamente eseguito solo una volta e ha usato l'intero array, ma in una rappresentazione stringa, che non è davvero quello che volevamo.

E poi il comando è uscito così, e quella era la directory di lavoro dove possiamo andare e vedere i file per capire un po' di più.

Okay. Quindi quello che è successo è stato. Nextflow ha solo passato questo intero array come un singolo elemento del canale al processo, il che significava che il processo è stato eseguito solo una volta. Aveva un'attività e non ha usato i dati in una struttura che ci aspettavamo.

## 3.2. Usare un operatore per trasformare i contenuti del canale

Quindi dobbiamo fare qualcosa a questo canale prima, prima che possa essere usato. E questo sta preparando il terreno per usare gli operatori, che sono funzioni speciali che possiamo usare sui canali per manipolare i contenuti del canale.

In questo caso, useremo qualcosa chiamato _flatten_. Che passiamo alla fine del canale qui. Quindi creiamo il canale e poi eseguiamo _flatten_. E ancora, se ci passo sopra, mi mostra la documentazione per questo comando direttamente in VS Code, il che è molto utile. Potete anche trovare tutti questi documenti sul sito web di Nextflow, la documentazione.

Potrei semplicemente eseguire questo codice ora e vedere se funziona, ma è anche una bella opportunità per introdurre come fare codice dinamico all'interno degli operatori e all'interno del codice Nextflow, che sono chiamati closure.

Quindi aggiungerò di nuovo un comando view qui prima di eseguire _flatten_. E qui questo ha queste parentesi graffe, che è la closure dinamica. E c'è solo del codice arbitrario qui dentro che verrà eseguito, nel contesto di un operatore view.

Qui, questo sta dicendo prendi il saluto, che è l'input dell'operatore view, e quello è qui. Potrei chiamare questo come voglio, potrei chiamarlo _"foo"_ e devo solo riferirmi ad esso come _"foo"_ dopo. E poi dico con questo, restituisci questo.

E poi impostando restituire una stringa che dice prima del flatten per una variabile. molto semplice.

Ora ne aggiungerò un altro esattamente uguale, ma dirò dopo _flatten_.

Quindi quello che fa questo, perché questo viene eseguito in sequenza, vedrete come appare il canale prima di eseguire _flatten_, e poi di nuovo dopo aver eseguito _flatten_.

E poi questo canale greeting viene ancora creato, quindi verrà ancora passato al processo. E si spera ora il flusso di lavoro verrà eseguito. Proviamolo.

Fantastico. Quindi prima di tutto è che la pipeline non è andata in crash questa volta. Abbiamo avuto tre processi che sono stati eseguiti correttamente e abbiamo un piccolo segno di spunta. E poi possiamo vedere che le nostre dichiarazioni view hanno funzionato.

Abbiamo prima _flatten_, che è quell'array che abbiamo visto prima dal fallimento, e poi abbiamo tre volte il dopo _flatten_ è stato chiamato dove abbiamo _"Hello", "Bonjour"_, e tutti quegli altri tre elementi separati nell'array, che ora sono come speravamo, tre elementi separati nel canale.

E potete vedere che l'operatore _view_ è stato eseguito tre volte. E questo perché questo canale dopo _flatten_ ha ora tre elementi. E quindi l'operatore viene chiamato tre volte.

Molto velocemente, vorrei solo menzionare che quando stavo creando fabbriche di canali prima, ho fatto _"."_, e poi abbiamo visto che c'erano molti modi diversi per creare canali, e uno di questi si chiama "_fromList"_. E quello è effettivamente specificamente progettato per fare questa stessa operazione. Quindi avremmo potuto semplicemente fare from list greetings array, e funzionerà. È una sintassi leggermente più pulita e bella. Ma per gli scopi di questa dimostrazione, volevamo renderlo un po' più passo dopo passo così poteste vedere come il canale viene manipolato e come diversi operatori possono cambiare cosa c'è nel contenuto di un canale.

## 4. Leggere i valori di input da un file CSV

Okay, come possiamo rendere questo un po' più realistico? Probabilmente non vorrete creare un sacco di codice nella vostra pipeline Nextflow con array codificati. Probabilmente vorrete prendere i dati dall'esterno quando lanciate, e quei dati saranno quasi certamente in file.

Quindi la prossima cosa che faremo è replicare questo, ma invece di prendere i dati da un singolo parametro CLI o da una stringa o array codificati, li prenderemo da un file.

Quindi eliminiamo il nostro greetings array. E ora cambieremo di nuovo questa fabbrica di canali. Ho appena detto che ce n'erano un gruppo tra cui scegliere e ce n'è uno chiamato _".fromPath"_. E gli dirò di, in questo caso, prendere _params.input_, che sta tornando al nostro input che stavamo usando prima.

Ora quel parametro non è davvero pronto per essere usato ancora. Stiamo ancora dicendo che è una stringa ed è codificato qui con un default, ma potremmo sovrascrivere quella stringa. Ora vogliamo che questo sia un file invece. Quindi il tipo è diverso. Non è più una _String_. È un _Path_.

E poi possiamo impostare il default se vogliamo, ancora, a un Path. E se guardo in esplora sulla sinistra, potete vedere in questo repository, in questa directory di lavoro, ho una directory chiamata data. Ho un file lì chiamato _"greetings.csv"._

Quindi posso semplicemente impostare il default qui a _"data/greetings.csv"_. Ora, quando eseguo questa pipeline di nuovo senza opzioni da riga di comando, userà questo valore predefinito. Sa che è un path, quindi sa che dovrebbe gestirlo come un path e non una stringa.

E poi passerà quello in una fabbrica di canali da questo _params.input_ e creerà il nostro canale, che poi verrà usato in questo processo chiamato _sayHello_. Proviamolo.

Okay. Fallito. Non preoccupatevi. Questo era previsto. E se state seguendo il materiale di formazione, vedrete che era previsto anche lì. Vediamo cosa sta succedendo qui.

Ha provato a eseguire la pipeline. Ha provato a eseguire il processo, e ha ottenuto un errore abbastanza simile a quello che abbiamo visto prima.

Qui dice: abbiamo provato a eseguire _echo_, ma invece di echeggiare i contenuti di questo file CSV, ha solo echeggiato il path. E potete vedere che è il path assoluto completo qui a questo file CSV.

E poi infatti, perché ha provato a scrivere quello in questo path davvero complicato, non sapeva davvero cosa fare. Ed era fuori dall'ambito della directory di lavoro del processo.

Ho menzionato all'inizio che Nextflow incapsula ogni attività eseguita all'interno di una speciale directory di lavoro. E se provate a scrivere su dati, che sono al di fuori di quella directory di lavoro, Nextflow vi fermerà come precauzione di sicurezza. E questo è quello che è successo qui. Proviamo a scrivere su un path assoluto e Nextflow è fallito e ci ha impedito.

## 4.2. Usare l'operatore splitCsv() per analizzare il file

Okay, diamo un'occhiata a questo canale e vediamo com'è. Possiamo fare _".view"_, e l'ho copiato dal sito web. Quindi _.view_, e abbiamo una closure dinamica qui e diciamo un nome di variabile "_csv"_ come input. Quindi quello è il contenuto del canale, e diciamo prima di splitCsv, e questo è come appare.

Se lo eseguo di nuovo, fallirà ancora, ma ci mostrerà cosa c'è dentro questo canale. Non è particolarmente eccitante. È quella variabile _path_. Quindi potete vedere che è solo una stringa qui perché viene stampata su un terminale, ma è un oggetto _path_, che contiene le informazioni e i metadati su questo file.

Non vogliamo passare i metadati del file all'input. Vogliamo passare i contenuti di quel file. Se guardiamo il file _greetings.csv_, potete vedere qui che ha queste diverse variabili qui. _Hello, Bonjour, Holà_ di nuovo. E queste sono le cose che davvero vogliamo passare al nostro processo, non solo il file stesso come un singolo oggetto.

Quindi dobbiamo analizzare questo file CSV. Dobbiamo disimballarlo, arrivare ai contenuti del file CSV, e poi passare i contenuti all'interno del canale al processo.

Come probabilmente potete dire dal messaggio di log, vogliamo usare _splitCsv_, che è un altro operatore, un altro operatore di canale. Quindi se faccio "_dot" "s"_, e poi potete vedere che è auto-suggerito. Oops, _splitCsv_ e alcune parentesi.

E poi dopo _splitCsv_, metterò un'altra dichiarazione _view_ solo così possiamo vedere come appare dopo. Eseguiamo la pipeline e vediamo cosa abbiamo.

Okay. È ancora fallito, ma in un modo nuovo ed eccitante, che è un progresso.

Questa volta ancora, abbiamo qualche problema con il nostro script, che è stato renderizzato. Ora. Non abbiamo più il path finale, ma abbiamo un array di variabili, che sembra molto simile all'errore che avevamo prima quando stavamo passando un array come input fisso.

Con il nostro logging dall'operatore view, possiamo vedere prima che _splitCsv_ fosse il path. E infatti, dopo _splitCsv_, abbiamo tre output diversi e ognuno di quegli output sembra molto simile a ciascuna delle righe del file _greetings.csv_, il che ha senso.

Quindi quello che è successo qui è che Nextflow ha analizzato questo file CSV dandoci tre oggetti, un array per ogni riga del file CSV. Quindi poi tre volte abbiamo passato un array di variabili al canale invece di un singolo valore stringa.

Okay, quindi l'ultima volta che avevamo questo problema, abbiamo usato _flatten_. Proviamo molto velocemente. Proviamo flatten e vediamo cosa succede.

Posso chiamare queste variabili, qualunque cosa. Quindi lo chiamerò _myarray_ perché non è più davvero un CSV. Proviamo ad eseguirlo di nuovo e vediamo cosa succede con _flatten_.

Quindi questa volta eseguiremo, abbiamo analizzato il CSV in tre oggetti array, e poi l'abbiamo appiattito. E questa volta, è passato. E la pipeline Nextflow è stata eseguita. Tuttavia potete vedere che _flatten_ va davvero in città e appiattisce tutto. E quindi otteniamo tre voci array indipendenti per ogni riga. E quindi ha eseguito il processo tre volte per ogni riga di un CSV. E ora abbiamo un sacco di file di risultati, e 123, 456, e tutti i tipi di cose, non solo quella prima colonna del CSV, che è quello che volevamo davvero.

## 4.3. Usare l'operatore map() per estrarre i saluti

Quindi come arriviamo solo alla prima colonna? Se flatten è troppo semplicistico qui, abbiamo bisogno di un operatore più complesso dove possiamo effettivamente personalizzare e dirgli cosa vogliamo dal CSV.

Per fare questo, useremo _map_. Fondamentalmente _map_ dice solo, esegui del codice, qualche funzione su ogni elemento che mi viene dato e fai qualche tipo di trasformazione su di esso. E poiché è così flessibile, lo vedrete comparire nel codice Nextflow tutto il tempo.

Da solo, non fa nulla. Quindi non vogliamo parentesi normali, vogliamo una closure qui e dobbiamo dirgli cosa fare. Quindi dirò _"row"_, perché gli vengono date righe dal CSV, quindi è un nome di variabile logico. È l'input. E voglio restituire solo il primo elemento di quell'array.

Gli array in Nextflow sono zero based, quindi diremo solo il primo elemento, che è row zero. Se volessimo la seconda colonna, potrei essere uno o la terza colonna essere due, e così via. Possiamo restituire quello che vogliamo qui, ma restituirò solo il primo valore.

E ora, possiamo eseguire la pipeline di nuovo e vedere se fa quello che ci aspettiamo.

Infatti, dopo _splitCsv_ abbiamo i nostri array, e poi dopo la _map,_ abbiamo le nostre belle stringhe pulite, solo _"Hello", "Bonjour"_ e _"Holà"_. E la pipeline sta ora facendo quello che vogliamo. Fantastico.

Quindi possiamo sbarazzarci di tutti questi comandi view ora. Non ne abbiamo più bisogno.

## Riepilogo

Abbiamo finito il nostro tipo di debugging e questo è il codice con cui finiamo. Prendendo il nostro parametro CLI chiamato _input_, che è classificato come un _Path_. Nextflow trova il path, lo carica, e capisce il file CSV. Restituisce tutte le diverse righe. E poi mappiamo solo il primo elemento di quella riga nel canale che ci dà il contenuto del canale, che viene passato al processo.

E il processo viene eseguito su ogni elemento nel canale, che è tre. E esegue il processo tre volte, dandogli tre attività. E quei risultati vengono poi pubblicati dal flusso di lavoro, raccolti dall'output del processo. Pubblicati da un flusso di lavoro e salvati nel blocco output in una sottodirectory chiamata _"hello_channels"_.

Piuttosto bello. Stiamo arrivando ora a qualcosa che assomiglia di più a una pipeline Nextflow di vita reale che potreste eseguire per qualche analisi reale.

## Takeaway

Okay. Si spera che ora stiate acquisendo una sensazione di cosa sono i canali e gli operatori di Nextflow e come gli operatori lavorano sui canali e come potete crearli.

I canali, come ho detto all'inizio di questo video, sono la colla di Nextflow. E potete vedere qui che possiamo prendere input diversi e manipolarli e prendere quei dati e poi passarli alla logica del flusso di lavoro downstream.

E questo blocco workflow qui è davvero dove costruite tutta quella parallelizzazione e tutta la logica intelligente, e spiegate a Nextflow come costruire il vostro DAG del flusso di lavoro, e come orchestrare la vostra pipeline.

I canali non sono il concetto più facile da capire. Quindi fate una pausa, pensateci un po', magari leggete di nuovo il materiale, e assicuratevi davvero di aver assimilato questi concetti perché questo è la chiave per la vostra comprensione di Nextflow e meglio capite i canali e i diversi operatori di canale e le diverse fabbriche di canali. Più vi divertirete scrivendo Nextflow e più potenti saranno le vostre pipeline.

Questo non è lo stesso della programmazione normale in Python o altri linguaggi. Non stiamo usando dichiarazioni _if_ qui, questa è programmazione di flusso funzionale usando canali e operatori. Quindi è un po' diverso, ma è anche super potente.

Questa è la fine di questo capitolo. Andate a fare una piccola pausa e vi vedrò nel prossimo video per la parte tre dove passeremo attraverso Hello Workflow, e parleremo un po' di più sui flussi di lavoro.

Proprio come il capitolo precedente, ci sono alcune domande del quiz in fondo alla pagina web qui, quindi potete farle velocemente e assicurarvi di aver capito tutte le diverse parti del materiale che abbiamo appena fatto. E a parte questo, vi vedrò nel prossimo video. Grazie mille.

Okay.
