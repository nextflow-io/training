# Parte 3: Hello Workflow - Trascrizione Video

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduzione assistita da IA - [maggiori informazioni e suggerimenti per miglioramenti](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/_aO56V3iXGI?si=Irl9nAQniDyICp2b&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Note importanti"

    Questa pagina mostra solo la trascrizione. Per le istruzioni complete passo dopo passo, tornate al [materiale del corso](../03_hello_workflow.md).

    I numeri di sezione mostrati nella trascrizione sono forniti solo a scopo indicativo e potrebbero non includere tutti i numeri di sezione presenti nei materiali.

## Benvenuto e riepilogo

Ciao e bentornati alla parte tre di Hello Nextflow. Questa parte si chiama Hello Workflow, ed è proprio in questa parte del corso dove iniziamo davvero a giustificare il nome pipeline o flusso di lavoro.

Prenderemo il nostro semplice script di pipeline finora con il suo unico processo, e inizieremo ad aggiungere ulteriori processi e vedremo come Nextflow gestisce questa orchestrazione e il flusso di dati attraverso la pipeline.

Torniamo ai nostri code spaces. Vedrete che ho eliminato tutte le directory .nextflow\* e le directory di lavoro e tutto il resto per cercare di mantenerlo pulito. Non preoccupatevi se avete ancora quei file rimasti dalle parti precedenti del corso.

Lavoreremo da un file chiamato hello-workflow.nf. Come prima, questo rappresenta fondamentalmente lo script che abbiamo costruito fino a questo punto, e ci dà un punto di partenza pulito. E ancora, giù nell'output possiamo vedere che il percorso ora è hello_workflow. Quindi i file pubblicati dovrebbero andare in una sottodirectory diversa nella vostra cartella dei risultati.

Per ricapitolare dove siamo finora, abbiamo un singolo processo qui, con un input greeting, un output greeting file. E poi il semplice script Bash, che fa solo un comando echo in un file.

Abbiamo un singolo input del flusso di lavoro, il blocco params qui, dove diciamo che si aspetta un percorso, e il valore predefinito è data/greetings.csv, che è questo file qui sopra.

Poi nel flusso di lavoro stesso, abbiamo un blocco main. Stiamo creando un canale. Stiamo analizzando il CSV in righe e poi prendendo il primo elemento di ogni array, e stiamo passando quel canale a quel processo, che poi genera tre attività, e stiamo pubblicando dal flusso di lavoro, gli output di quel processo.

E poi infine, nel blocco output, stiamo dicendo a Nextflow di pubblicare questi file da questo canale alla directory chiamata hello_workflow. E di copiare quei file piuttosto che collegarli con soft link.

## 1. Aggiungere un secondo passaggio al flusso di lavoro

Okay, in questa parte aggiungeremo un secondo processo al nostro flusso di lavoro. Prenderemo gli output del processo sayHello, e li elaboreremo in un secondo passaggio, che convertirà tutte le lettere all'interno di quei file in convertToUppercase.

Questo è solo un esempio stupido, è solo un po' di semplice elaborazione di stringhe ancora, ma vi mostra come possiamo prendere la logica, all'interno del flusso di lavoro.

Useremo un comando bash chiamato "tr" per questo, che sta per translate. È un comando Unix che esiste da sempre. Se non lo conoscete, non vi biasimo. Non penso di averlo mai usato prima della formazione, ma potete provarlo molto velocemente sul terminale. Se faccio "echo 'hello world'" e poi pipe a 'tr' e poi tra virgolette dite intervallo di caratteri, quindi da A a Z, minuscolo, e poi volete fare da A a Z maiuscolo. E dice semplicemente, traduci queste lettere in queste lettere.

E quando premo invio, potete vedere che ora ha messo tutto in maiuscolo. Molto carino se vi piace urlare alle persone.

Quindi questo è uno stile molto semplice di comando bash che useremo nel nostro secondo processo.

## 1.2. Scrivere il passaggio di conversione in maiuscolo come processo Nextflow

Quindi se torno al mio script, barò un po' e copierò semplicemente il codice dalla, dalla documentazione della formazione. Ma potete vedere esattamente cosa sta succedendo.

Abbiamo un nuovo processo qui. Questo l'abbiamo chiamato convertToUpper, ma potremmo chiamarlo come vogliamo.

Abbiamo un singolo input path, come abbiamo fatto prima. Non è un canale di valori, è un canale di percorsi. E poi un singolo output.

Nel blocco script facciamo "cat" sul file di input. E possiamo metterlo tra parentesi graffe se vogliamo. e che prende quella variabile. E eseguiamo lo stesso comando bash nella pipe e scriviamo i risultati in un file con questo nome di file, e quello viene raccolto dall'output path.

Ora dobbiamo fare qualcosa con questo nuovo processo. Quindi andremo giù al flusso di lavoro dove costruiamo la diversa logica di un flusso di lavoro, e dopo quel primo processo, eseguiremo il nostro secondo processo. Quindi convertToUpper è il nome del processo qui.

Prende un input quindi non possiamo semplicemente chiamarlo da solo. Vogliamo elaborare l'output del primo processo. Quindi proprio come abbiamo fatto con questo, sayHello out dove stiamo pubblicando quei risultati. Vogliamo usare quegli stessi risultati qui come input, quindi possiamo copiarli e metterli lì.

Vogliamo il processo sayHello ".out", e Nextflow sa che questo significa un semplice singolo record di output qui, che è questo file. Quindi quello sarà poi passato come input a un secondo processo.

## 1.5. Configurare la pubblicazione degli output del flusso di lavoro

Okay. E infine, così che salviamo effettivamente i risultati di questo secondo processo, dobbiamo anche pubblicarli dal flusso di lavoro, e poi definirli nel blocco output, stessa sintassi di prima. Quindi possiamo copiare questo e dire second outputs, o come volete chiamarlo.

Prendere il nome del processo che ci interessa, convertToUpper out, e poi qui giù nel blocco output. Aggiungere questo e potremmo fare gli stessi attributi qui. Quindi vogliamo anche questi file nella sottodirectory Hello Workflow, e vogliamo anche copiarli.

Ottimo. Proviamo ad eseguirlo. Quindi se apro il terminale e faccio "nextflow run hello-workflow.nf", e vedremo cosa fa. Vediamo se sembra diverso dalle parti precedenti.

Quindi lancia Nextflow. Nella documentazione, dice di farlo con "-resume", ma ho eliminato tutta la mia directory di lavoro, quindi non avrebbe fatto alcuna differenza qui. Ma se l'avete fatto, allora funzionerà anche quello.

E sembra quasi esattamente lo stesso. Ma potete vedere ora c'è una seconda riga di output qui, dove potete vedere il nome del secondo processo che abbiamo appena aggiunto. E sicuramente, potete vedere che è stato eseguito tre volte con successo.

Brillante. Se avessi avuto le mie directory di lavoro precedenti e avessi fatto questo con "-resume", queste sarebbero state, messe in cache solo il primo passaggio nella pipeline. Perché quegli output erano esattamente gli stessi, quindi Nextflow avrebbe saputo di riutilizzarli di nuovo.

E quindi potete vedere come potete usare -resume per costruire iterativamente il vostro flusso di lavoro, passo dopo passo, se ne avete bisogno.

Okay, diamo un'occhiata alla directory dei risultati qui sopra e vediamo se ha funzionato. Possiamo vedere che abbiamo alcuni altri file qui sopra. Abbiamo i nostri file originali come prima dal primo processo. E sicuramente, abbiamo i nostri file upper e le lettere sono tutte maiuscole, quindi ha funzionato. È davvero bello da vedere.

È anche interessante controllare all'interno di queste directory di lavoro. Come prima l'hash qui corrisponde alle directory di lavoro. Quindi se guardo in "ls work", e poi espando quello, vedremo i diversi file qui.

Vediamo il file di output dal primo processo, che è stato inserito qui come input. E possiamo vedere il nuovo file di output che è stato generato.

Ora se faccio questo con "-la" per elencare e mostrare tutti i file, vedremo alcune cose in più. Prima di tutto, vedrete che questo file è in realtà un soft link al primo processo. Questo è fondamentalmente sempre un soft link se può esserlo, per risparmiare spazio su file. Non stiamo pubblicando i file qui e fa semplicemente riferimento a quel file da una prima attività in una seconda attività in modo che tutto sia incapsulato all'interno di quella directory di lavoro, e sicuro e isolato da tutto il resto.

E questo deve essere lì perché se guardiamo il file .command.sh, quindi se faccio "cat work/b8/56\*", potete vedere che le parti di file qui sono relative, quindi sta facendo cat di quel file di input, che è stato collegato con soft link nella stessa directory di lavoro.

Quindi è così che apparirà ogni directory di lavoro. Quando la guardate in Nextflow, avrete tutti i file di input lì messi in stage in quella directory di lavoro. E poi avrete anche tutti i file di output che sono stati creati. Quindi è ottimo. Sembra come ci aspettavamo.

## 2.1. Definire il comando di raccolta e testarlo nel terminale

Okay, torniamo al nostro flusso di lavoro. Qual è il prossimo passaggio che vogliamo fare?

Abbiamo due processi ora e stanno prendendo questo singolo file CSV, analizzandolo e dividendolo. E poi abbiamo tre attività per ciascuno di questi processi e Nextflow gestisce la parallelizzazione di tutto questo, quindi tutto viene eseguito fianco a fianco dove possibile.

Quel tipo di modo di dividere il lavoro per eseguire cose in parallelo è molto comune. E l'inverso di quello è poi raccogliere tutto insieme. Quindi è quello che faremo con il nostro processo finale nel flusso di lavoro è avremo un terzo qui, che prende questi tre diversi output e li combina tutti in un singolo file.

Possiamo farlo abbastanza semplicemente in un terminale, solo per avere un'idea di come sarà questo.

Se vado alla cartella dei risultati. Quindi, "cd results/hello_workflow/", e abbiamo tutti i file UPPER qui. Posso solo usare "cat", che usiamo per stampare i contenuti di quel file, e potete dare più file a "cat" e li leggerà uno dopo l'altro.

Quindi posso dire "UPPER-\*", che mi dà la stessa lista di tre nomi di file con l'espansione Bash. E posso dire combined.txt. Penso che nella documentazione, elenca i nomi esatti dei file, ma sta facendo la stessa cosa.

Ora, se uso "cat combined.txt", possiamo vedere che abbiamo i contenuti del file di tutti e tre quei file.

Quindi è fondamentalmente tutto ciò che questo processo farà è cercheremo di dargli tutti i diversi file di output da un processo precedente in una singola attività di processo, e poi faremo "cat" insieme e salveremo il file di output.

## 2.2. Creare un nuovo processo per fare il passaggio di raccolta

Okay, quindi aggiungiamo il nostro nuovo processo. Incollerò questo dai materiali di formazione, e potete vedere che ci ha lasciato un po' di esercizio per il lettore qui con questi punti interrogativi. Ma potete vedere il profilo generale del processo è fondamentalmente quello che abbiamo appena fatto nel terminale, dove stiamo facendo "cat" di un gruppo di file di input e scrivendolo in un file di output qui chiamato collected, e poi l'output si aspetta quel singolo percorso di nuovo.

Quindi abbiamo bisogno di qualche tipo di input qui e saranno un insieme di percorsi. Quindi ancora, definiamo un canale input path e chiamiamolo input_files. Ora, questo precedentemente ci ha dato un singolo percorso qui, ma un percorso può anche avere più file qui, anche se è ancora una singola dichiarazione.

Lo copierò qui sotto perché vogliamo fare "cat" di questi file. E potreste pensare che abbiamo alcuni problemi qui con la stampa di un array o cose del genere, ma Nextflow è generalmente abbastanza sensato quando si tratta di questo. E se gli viene dato un canale con più file in esso come questo, li metterà tutti insieme con separatori di spazio. Quindi questo ci darà la sintassi corretta.

È ottimo. Quindi ora colleghiamo il nostro nuovo processo. Vado giù al flusso di lavoro. Dirò combine the outputs, il nuovo nome del processo, e proprio come prima. Prenderò questo processo precedente, convertToUpper e farò ".out".

Ottimo. Proviamolo e vediamo se funziona nel terminale. Se torno semplicemente su un paio di directory e poi rieseguo il comando Nextflow, e vedremo cosa succede.

Quindi il flusso di lavoro è stato lanciato e ora potete vedere che abbiamo tre diversi nomi di processo, il che è ottimo. I primi due sembrano entrambi gli stessi di prima, e il terzo nuovo viene eseguito, il che è buono.

Tuttavia, c'è qualcosa di un po' strano qui. Volevamo combinare quei file di output in un singolo file, eppure questo processo che possiamo vedere è stato eseguito tre volte, non una.

Sicuramente, se andiamo in una di queste directory di lavoro. E facciamo "cat work/" "collected", allora vedremo. C'è solo una singola parola qui, non tre.

E quindi quello che è successo è che Nextflow ha continuato quella parallelizzazione proprio come ha fatto nei passaggi precedenti. E questo processo ci ha dato un canale con tre elementi, e quei tre elementi del canale sono stati passati a quel nostro processo downstream, che ha generato tre attività di processo.

Ha fondamentalmente provato a raccogliere tre volte separate e ogni volta aveva solo un singolo file, quindi ha solo fatto cat singolo file in un output, e infatti, possiamo vedere anche quello nel file .command.sh.

Se faccio .command.sh, possiamo vedere che ha solo un singolo nome di file qui e solo un singolo file è stato messo in stage in quella directory di lavoro.

## 2.3. Aggiungere il passaggio di raccolta al flusso di lavoro

Quindi in qualche modo dobbiamo dire a Nextflow di riunire tutti quegli output da un processo precedente e darli a questo processo downstream come un singolo elemento del canale, piuttosto che tre.

Facciamo questo con un operatore di canale chiamato _collect_.

Questo è un operatore super utile, che vedrete nelle pipeline Nextflow tutto il tempo. Questo è un canale qui, questo canale di output, proprio come quello che abbiamo creato in alto. E quindi possiamo aggiungere operatori di canale ad esso proprio come abbiamo fatto prima. Possiamo solo fare punto, e poi in questo caso, collect, parentesi.

E questo è tutto ciò di cui abbiamo bisogno. Questo manipolerà poi questo canale prima che sia passato in questo processo.

Se volete vedere cosa gli sta succedendo, possiamo anche visualizzarlo qui. Quindi qui, questo non è correlato all'esecuzione di questo processo affatto, quindi potrei metterlo in qualsiasi punto dopo aver eseguito quel processo. Ma prendiamo lo stesso canale di output, e lo stiamo guardando con .view, e poi lo stiamo guardando di nuovo con .collect.view.

E quando eseguiamo questo, ci mostrerà le due diverse strutture di quel canale, prima e dopo collect. Quindi proviamolo ora. Okay, ho appena ridotto un po' perché alcuni degli output sono piuttosto lunghi, ma se eseguo la pipeline, vedremo se funziona.

Spero che un terzo processo venga eseguito solo una volta, perché sta raccogliendo gli output e sicuramente, potete vedere collectGreetings come uno di uno. Quindi è stata eseguita solo un'attività.

E poi se guardiamo le dichiarazioni view, abbiamo tre dichiarazioni view per i tre elementi di prima, con un percorso di file in ciascuno.

E poi dopo quella dichiarazione collect, è stata attivata solo una volta perché c'è un singolo elemento in quel canale. E ora abbiamo questa lista di tre diversi percorsi di file.

È esattamente quello che speravamo. E potete vedere si spera, questo è fondamentalmente l'inverso di quell'operatore "map" che abbiamo fatto per passare dagli array CSV in elementi del canale separati. Ora stiamo prendendo elementi del canale separati e rimettendoli in un singolo array.

Ottimo, possiamo eliminare queste dichiarazioni view. Non ne abbiamo più bisogno. Possiamo passare al prossimo passaggio.

Prima di andare oltre, e prima che mi dimentichi, aggiungerò una nuova dichiarazione publish qui. Third output. Potete chiamarlo qualcosa di più semantico e descrittivo nel vostro flusso di lavoro. E poi lo aggiungerò al blocco output di nuovo e dirò path 'hello_workflow' mode 'copy'. Così che il file di output generato da questo processo sia salvato nella nostra cartella dei risultati qui sopra.

Solo per verificare rapidamente che funzioni. Dovrebbe essere un po' più pulito ora perché non abbiamo quelle dichiarazioni view. E, vedremo se otteniamo il nostro nuovo file di output qui sopra. Una di, un'attività è stata eseguita, ho ottenuto un nuovo file chiamato collected, e ora abbiamo tutte e tre quelle parole. Fantastico. Cosa c'è dopo?

## 3. Passare parametri aggiuntivi a un processo

Okay. Successivamente guarderemo alla gestione di più input in un singolo processo. Finora potete vedere che tutti i nostri processi stanno solo prendendo una cosa come input. Hanno tutti una singola riga sotto il loro input.

Dimostreremo questo permettendo a Nextflow di specificare un diverso identificatore di batch in modo che forse eseguite questo flusso di lavoro più volte e potete dargli un ID batch diverso ogni volta.

Aggiungerò semplicemente una seconda riga nell'input qui per collectGreetings. E lo chiamerò "val", perché questa è una stringa. Ora è un valore, non un percorso, e lo chiamerò "batch_name".

Poi modificherò lo script qui sotto per usare questa variabile, e proverò a metterla nello stesso posto del materiale di formazione. Quindi la metto nel mezzo di questo percorso file COLLECTED-$\{batch_name\}-output.

Non è ancora finita. Ricordate che dobbiamo dire a Nextflow quali saranno i nomi dei file di output. Quindi dobbiamo anche fare la stessa cosa qui sopra: COLLECTED-$\{batch_name\}-output.txt".

Fantastico. Nextflow ora sta ottenendo un secondo input variabile e lo sta interpolando nello script e nell'output.

Un'ultima cosa, ora dobbiamo trovare dove questo viene chiamato, e dobbiamo passare il secondo input al processo. Questo è proprio come qualsiasi altro input in una funzione in qualsiasi altro linguaggio.

Proprio come abbiamo fatto prima nella formazione, userò lo speciale "params" qui, e lo chiameremo "params.batch" così che possiamo avere un'opzione CLI -- batch. E ora potete vedere che il nostro processo qui ha due input separati solo separati da virgola, che vengono passati.

È davvero importante ottenere l'ordine giusto, quindi l'ordine degli argomenti qui per channel e poi il param deve corrispondere. Il canale e il batch name lì. Questo è solo corrispondenza posizionale.

Okay. Posso eseguire questa pipeline ora subito con --batch, ma facciamo prima la cosa giusta e definiamolo nell'input qui in Params. Quindi lo aggiungerò a batch e poi diremo che è una stringa e diamogli un valore predefinito. Quindi chiamiamolo semplicemente batch. Okay? Ora proviamo ad eseguire il flusso di lavoro.

--batch Trio. Penso che dica nel materiale di formazione, ma potremmo usare qualsiasi stringa vogliamo lì. E si spera vedremo quel file di output dei risultati venire qui sopra.

E sicuramente, COLLECTED-trio-output - ha funzionato correttamente. Ha rinominato il nostro file. E potete immaginare ora questo è utile perché se eseguo di nuovo quello con un nome batch diverso, come replicate_two, allora ci darà un nome batch diverso qui sopra.

E e non sovrascriverà poi i file di output in questo caso. Quindi è carino.

## 4. Aggiungere un output al passaggio di raccolta

Okay, quindi ora abbiamo più input al nostro processo qui. Ma cosa succede se vogliamo creare più output? Il nostro esempio qui allora è che creeremo un report per questo processo, solo dicendo questo è quanti file sono stati raccolti.

E lo faremo con un comando echo qui. Quindi possiamo dire echo. There were, copierò questo da un materiale di formazione, così non dovete guardarmi digitarlo.

There were $\{count_greetings\} greetings in this batch, e salvare quello in un nuovo file ora chiamato $\{batch_name\}, quindi stessa variabile, possiamo riutilizzarlo tutte le volte che vogliamo, report.txt.

## 4.1.1. Contare il numero di saluti raccolti

Dobbiamo effettivamente calcolare questo in qualche modo. Potremmo fare quella logica nello script Bash se volessimo, usando la logica Bash. Tuttavia, possiamo anche solo fare scripting direttamente all'interno del codice Nextflow, purché sia all'interno del blocco script nel processo e sopra la sezione citata.

Qualsiasi cosa qui non sarà inclusa nello script finale renderizzato, e sarà solo eseguita da Nextflow quando renderizza un'attività.

Quindi qui stiamo solo facendo un po' di logica. Stiamo creando una nuova variabile chiamata count_greetings. Prendiamo il canale dei file di input qui, e stiamo chiamando .size() su di esso.

Okay, quella funzione mi darà un numero qui in questa variabile, e ora il nostro avviso è sparito perché questa variabile viene definita.

Okay, quindi stiamo creando quel secondo file nella directory di lavoro, ma dobbiamo dire a Nextflow di aspettarselo come output pubblicato di questo processo. Quindi lo facciamo con esattamente la stessa sintassi che abbiamo fatto per il primo file.

Diciamo path perché è, ancora, potremmo pubblicare una variabile qui se volessimo con "val", ma diremo "path". E poi il nome del file previsto. Notate che non è evidenziato qui. È perché ho usato virgolette singole. Devo usare virgolette doppie.

## 4.1.2. Emettere il file di report e nominare gli output

Okay, è ottimo. E potremmo ora iniziare ad accedere a questi output qui sotto proprio come ho fatto qui. Ma ora è un array di oggetti diversi, quindi potrei fare collectGreetings.out[0] per ottenere il primo, o uno per ottenere il secondo, che è il nostro nuovo report.

Ma non mi piace molto farlo perché è abbastanza facile sbagliare il conteggio dell'indice. E state lì a contare righe molto e aggiungete un nuovo output e improvvisamente tutto si rompe. Quindi

è molto più carino riferire tutto per nome invece. E possiamo farlo con una chiave speciale qui chiamata "emit".

Quindi possiamo chiamare questo come vogliamo. Chiamiamolo emit outfile, e emit reports. Se le definite e potete farlo su una o molte, sta a voi. Ora posso andare qui sotto e invece posso andare punto out punto reports e chiamarlo semplicemente per nome, che è molto più facile capire il vostro codice quando lo leggete, ed è più sicuro per le modifiche nel codice.

Ho aggiunto .out.report qui, ma in realtà ho bisogno di avere due diversi output pubblicati. Quindi rinominerò come qualcosa di più interessante come collected e report ed è quello che l'ho chiamato? L'ho chiamato out file, scusa. Quindi quel nome emit qui outfile e report. perché stiamo pubblicando due diversi canali di output e quindi dobbiamo riferirli entrambi nel blocco publish.

Poi dobbiamo anche definire questi nel blocco output. Quindi ho rinominato quello collected, e ancora, per reports, un po' verboso qui, ma è davvero utile quando entrate a leggere un nuovo flusso di lavoro, vedere tutti i diversi output qui, tutti i diversi canali elencati fianco a fianco, e ci sono modi per rendere questo meno verboso, che toccheremo più tardi.

Okay, proviamolo ed eseguiamo il nostro flusso di lavoro e vediamo cosa succede.

Si spera ora dovrebbe essere eseguito fondamentalmente come prima. E otterremo un nuovo file di output qui sopra chiamato replicate_two, report. E eccolo lì. È aperto e dice ci sono tre saluti nel batch, che è quello che ci aspettavamo, quindi è perfetto.

Se vado nella directory di lavoro qui solo per provare a usare è stato eseguito nel codice Nextflow, piuttosto che nello script bash, posso andare a cat work/ command.sh, e vedrete qui che sta solo facendo echo di questa stringa direttamente. There were three greetings in this batch, e quindi quella variabile è stata interpolata da Nextflow. È stata calcolata nel blocco script prima che scrivesse il file .command.sh. Quindi il calcolo della variabile risultante è fondamentalmente hardcoded in questo prima che sia eseguito sul vostro ambiente di calcolo in questo caso.

E quindi potete vedere quella separazione tra lo script. Blocco qui e qualsiasi cosa sopra di esso. Spero che abbia senso.

## Takeaway e quiz

Okay, questa è la fine di questa parte di Hello Nextflow. Quindi come prima, andate a controllare il quiz. Fatelo sulla pagina web o nella CLI, fate alcune domande e controllate semplicemente di aver compreso parte del materiale che abbiamo trattato. Vedete se c'è qualcosa lì che evidenzia qualcosa che potreste non aver capito. Non troppe domande. Carino e facile da fare. Oppure potete farlo sulla pagina web qui sotto pure.

E prendetevi una piccola pausa, fate una piccola passeggiata e tornate e unitevi a noi nella parte quattro di Hello, Nextflow, dove parleremo di moduli. Grazie mille.
