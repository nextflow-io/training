# Parte 1: Hello World - Trascrizione

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduzione assistita da IA - [scopri di più e suggerisci miglioramenti](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/8X2hHI-9vms?si=F0t9LFYYLjAWoyRXj&amp;list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Note importanti"

    Questa pagina mostra solo la trascrizione. Per le istruzioni complete passo dopo passo, tornare al [materiale del corso](../01_hello_world.md).

    I numeri delle sezioni mostrati nella trascrizione sono forniti solo a scopo indicativo e potrebbero non includere tutti i numeri di sezione presenti nei materiali.

## Benvenuti

Salve, benvenuti al Capitolo Uno di Hello Nextflow.

In questa prima parte di un corso di sei parti, esamineremo le basi fondamentali di Nextflow. Inizieremo eseguendo alcuni comandi in un terminale, e poi prenderemo quei comandi Bash e vedremo come integrarli in uno script Nextflow.

Proveremo ad eseguire quella prima pipeline Nextflow, vedremo cosa fa Nextflow, dove viene eseguito, quali file crea e qual è lo scopo di quei file.

Bene, cominciamo.

## training.nextflow.io

Prima di tutto, andate su training.nextflow.io. Proprio come prima, tutto il materiale è scritto qui, e lo seguirò passo dopo passo. Mostrerò il mio schermo mentre eseguo i passaggi della formazione, ma tutto ciò che sto dicendo è nel materiale di formazione quindi potete seguirlo al vostro ritmo, e potete trovare tutto scritto lì.

Questo video ha anche i sottotitoli abilitati, quindi sentitevi liberi di attivarli e seguire esattamente ciò che sto dicendo mentre lo dico.

Okay, andiamo a Hello Nextflow. Questo è il corso che faremo oggi, e abbiamo già fatto l'orientamento nel primo video, quindi andremo direttamente alla parte uno. Hello World.

Okay, lascerò ora questo materiale di formazione e passerò al mio ambiente Code Spaces. Questo è ciò che abbiamo configurato nel primo video. Spero che abbiate qualcosa di molto simile a questo nel vostro sistema. Sto usando VS Code e sto guardando il materiale di formazione e ho cambiato directory nella directory hello Nextflow.

## 0. Riscaldamento: Eseguire Hello World direttamente

Okay. Iniziamo con un paio di nozioni di base, che si spera sembrino familiari a tutti. Inizierò semplicemente scrivendo un comando molto semplice nel terminale. Qui sotto dirò 'echo Hello World!"' premo invio e, senza sorprese, il terminale fa quello che chiedo e restituisce quella stringa. Hello world.

Okay, poi premerò su per ottenere quel comando e lo modificherò un po'. Questa volta reindirizziamo quell'output a un file. Lo scriverò invece in output.txt e premerò invio niente sul terminale questa volta perché l'output non è arrivato al terminale. È andato in quel file.

Posso poi leggere quel file facendo 'cat output.txt' premo tab lì per espandere automaticamente il nome del file ed ecco fatto. Il file è lì.

Posso anche vedere quel file nella barra laterale nell'esploratore file in VS Code. Posso farci doppio clic e aprirlo qui. Se volete aprirlo in VS Code senza cliccare nulla, potete anche fare "code" e poi "output.txt" e fa la stessa cosa.

Ottimo. Questo è il primo passo. Molto semplice.

## 1. Esaminare lo script iniziale del workflow Hello World

Okay. Ora faremo esattamente la stessa cosa, ma in Nextflow, invece che direttamente nel terminale.

Useremo il primo script di esempio per iniziare, questo file si chiama Hello World. Posso fare "ls" per visualizzarlo in un terminale, e sono su Mac, quindi posso fare comando clic per aprire quel file, o avrei potuto semplicemente fare doppio clic nella barra laterale qui.

Ci sono alcune cose che possiamo vedere in questo file. Proprio in alto, c'è un'istruzione hash che dice che questo è un file Nextflow e che potrebbe essere eseguito. Ci sono alcuni commenti qui, solo commenti di codice normali in grigio chiaro, che non influenzano l'esecuzione, e aiutano solo a leggere lo script.

E poi ci sono due strutture principali. C'è un process qui e un workflow.

I processes in Nextflow sono i passaggi della pipeline. Sono le parti che effettivamente eseguono la logica e fanno l'elaborazione.

Il workflow poi in basso unisce questi processes e governa la logica del workflow, come tutto si connette l'uno all'altro.

Inizieremo guardando un process. Torneremo al workflow tra un momento.

## 1.2 La definizione del process

Quindi ogni process inizia con una parola chiave process. Ha un nome e poi ha alcune parentesi graffe e tutto all'interno di quelle parentesi graffe è quel singolo process.

Un process deve avere una sezione script, e qui contenuto c'è uno snippet bash in una stringa multilinea, che è la parte del codice che viene effettivamente eseguita nell'ambiente di calcolo.

Abbiamo anche un'istruzione output qui, che dice a Nextflow quali file sono previsti essere creati dallo script. Notate che l'output qui ha una parola chiave path, che dice a Nextflow che questo è un file, non un valore, o una stringa.

All'interno del blocco script, questa è solo un'istruzione bash normale, ed è esattamente la stessa di quella che abbiamo scritto nel terminale. Stiamo facendo echo di hello world in un file chiamato output.txt. Questo output.txt viene poi raccolto dalla definizione di output. La definizione di output in realtà non sta facendo nulla. Sta solo dicendo a Nextflow cosa aspettarsi, e se questo file non fosse stato creato, Nextflow avrebbe generato un errore.

Notate che questo esempio non è ottimo perché abbiamo codificato rigidamente il nome del file qui, output.txt e output.txt. Se uno di questi fosse stato modificato, ciò avrebbe causato un errore nel nostro workflow.

C'è un modo migliore per farlo con variabili, che tratteremo tra un minuto.

## 1.3 La definizione del workflow

Okay. Scendendo al workflow, possiamo vedere che abbiamo un commento e poi eseguiamo il process chiamato sayHello. Questa è la stessa parola chiave che è qui sopra. Questo è semplice quanto può essere un workflow. Stiamo solo chiamando un singolo process senza input variabili, quindi non lo stiamo collegando a nient'altro. Nella parte successiva di questo corso, parleremo di come rendere questo più potente usando input variabili e collegando le cose con i channels.

## 2. Eseguire il workflow

Okay, questo è tutto ciò di cui abbiamo bisogno. Vediamo se possiamo eseguirlo e vedere cosa succede. Pulirò solo il terminale e poi farò "nextflow run", e chiamerò il nome del file, che è hello-world.nf. Questo è tutto ciò di cui abbiamo bisogno per eseguire una pipeline Nextflow. Questa pipeline non prende alcun input, quindi non abbiamo bisogno di altri argomenti.

Premiamo invio e vediamo cosa succede.

Okay. Si spera che abbiate un output che assomigli a questo. Abbiamo alcune informazioni che ci dicono che Nextflow è stato eseguito e quale versione stava usando. Ci dice quale script è stato lanciato e ci dà un nome generato casualmente per questa particolare esecuzione del workflow. In questo caso, il mio si chiamava "gloomy_crick".

La parte più importante però, è che ci dice quali passaggi sono stati eseguiti nella pipeline. Potete vedere che il nostro process chiamato sayHello è stato eseguito, ed è stato eseguito una volta ed era completo al cento per cento.

Questa parte qui è l'hash per quella particolare attività del workflow. Ogni process viene eseguito una o più volte, e ciascuna di quelle esecuzioni è chiamata attività.

## 2.2. Trovare l'output e i log nella directory di lavoro

Ogni attività ottiene la propria directory isolata dove viene eseguita, quindi è separata dal resto dell'esecuzione del workflow. Questo hash corrisponde alla struttura dei file all'interno della directory di lavoro. Se faccio "tree work", possiamo vedere a0, e poi una versione più lunga di un hash breve, e poi il nostro file output.txt. Potete anche vederlo in una barra laterale.

Potete vedere nella barra laterale che ci sono alcuni file aggiuntivi qui. Il motivo per cui questi non sono apparsi in un terminale è perché sono file nascosti, iniziano con un punto. E infatti, se faccio "tree -a" per tutti, e "work", possiamo vederli qui.

Questi file punto sono presenti in ogni singola directory di lavoro. Che Nextflow crea, e ognuno ha un'attività leggermente diversa. In primo luogo .command.begin include solo alcune istruzioni per Nextflow che configura l'attività prima che venga eseguita. .command.run sono le istruzioni effettive eseguite da Nextflow stesso. Poi .command.sh è probabilmente quello più interessante. Questo è lo script che è stato risolto dal nostro blocco script del process.

Se lo apro, potete vedere che abbiamo il nostro "echo Hello World" nel file output.txt. Questo è esattamente lo stesso del nostro process in questo caso, ma se abbiamo variabili all'interno del nostro codice Nextflow, ogni attività avrà un .command.sh diverso, e potete vedere come quelle variabili sono state risolte.

Gli altri file riguardano come l'attività è stata eseguita. Quindi .command.err, .log e .out sono l'errore standard, l'output standard e i due combinati. E .exitcode dice a Nextflow come questa attività è stata eseguita con quale codice di uscita, se è riuscita o meno.

Infine, abbiamo il nostro file output.txt e certamente, "Hello World" questo è ciò che ci aspettavamo e questo è ciò che è stato creato.

Okay, ottimo. Questa è stata la vostra prima esecuzione Nextflow. Congratulazioni. È davvero così semplice.

Successivamente, vedremo come farlo in modo un po' più conveniente in modo da non dover modificare il codice ogni volta che vogliamo apportare una modifica al modo in cui viene eseguita la pipeline.

## 3. Gestire le esecuzioni del workflow

Questa struttura di directory è ottima per mantenere tutte le attività separate e tutto organizzato, ma ovviamente non è molto conveniente trovare i file di output. Non volete scavare attraverso molte directory nidificate cercando di trovare i risultati della vostra pipeline.

## 3.1. Pubblicare gli output

La buona notizia è che non è previsto che lo facciate. Le directory di lavoro sono davvero solo per Nextflow da usare. Quindi quello che faremo è che useremo una funzione di Nextflow chiamata "publishDir".

Torniamo al nostro workflow, andiamo al process. Possiamo aggiungere una nuova istruzione qui chiamata direttiva. Questo è ciò che Nextflow chiama queste cose in cima ai processes che aumentano il funzionamento, e quella che useremo si chiama publishDir.

Potete vedere che ho iniziato a digitare qui e l'estensione Nextflow per VS Code mi ha suggerito la direttiva, quindi posso semplicemente premere invio.

Okay, seguirò questo con una directory chiamata "results" e le diremo di copiare i file di output lì. Quindi dirò mode copy. Ottimo. premerò salva e rilanciamo il workflow.

nextflow run hello-world.nf

Viene eseguito esattamente allo stesso modo. Anche se notate abbiamo un hash leggermente diverso questa volta. Nextflow utilizzerà un hash diverso ogni volta che eseguite il workflow. E abbiamo un insieme diverso di directory di lavoro di conseguenza. Aree, una chiamata EB invece, ma potete vedere che tutti i file sono gli stessi. Tuttavia, ciò che è nuovo questa volta è che abbiamo anche una directory chiamata "results".

All'interno di "results" qui abbiamo il nostro file di output. Questo è ciò che abbiamo detto a Nextflow di fare. Abbiamo detto, salva i file di risultati in una directory chiamata "results" e copiali lì. E quindi questo è ora molto più facile da trovare. È proprio lì accanto a dove abbiamo lanciato un workflow e tutti i diversi file possono essere organizzati lì come desideriamo, indipendentemente da dove o come Nextflow ha eseguito l'esecuzione effettiva.

Notate che publishDir può gestire symlink, il che è buono se state lavorando su un file system condiviso e volete risparmiare spazio. E inoltre non dovete definire tutti i file che vengono creati da un process come output.

Nextflow copierà solo le cose che sono definite in questo blocco output. Quindi se avete file intermedi creati dal passaggio, che non sono necessari a valle di questo process, semplicemente non li definite in output e non appariranno nel publishDir. Quindi questo è un modo per mantenere i vostri file di output da una pipeline puliti e eliminare facilmente i file intermedi una volta terminato il posto di lavoro.

Una nota rapida qui. C'è una nuova sintassi Nextflow in arrivo chiamata workflow output definitions, che alla fine sostituirà publishDir. Questo ci dà un modo per definire tutti gli output da un workflow a livello di pipeline giù nel blocco workflow. Questo è descritto nei documenti Nextflow se volete provarlo. Ma per ora, publishDir sarà in circolazione per un po', quindi lo manteniamo ancora nella formazione per il 2025.

## 3.2. Rilanciare un workflow con -resume

Okay. Ho menzionato che la directory di lavoro qui ora ha due serie di risultati con un hash diverso da ogni volta che eseguiamo il workflow. Questo è buono. Tuttavia, a volte non vogliamo ricalcolare i passaggi ogni volta se non ne abbiamo bisogno.

Forse state costruendo iterativamente il vostro workflow e state aggiungendo passaggi e volete che i primi passaggi riutilizzino semplicemente le versioni memorizzate nella cache. O forse qualcosa è andato storto sul vostro sistema di calcolo a metà del vostro workflow e volete che continui da dove si era interrotto, ma saltando i passaggi che aveva già completato.

Nextflow ha funzionalità integrate per questo chiamate resume. Proviamole. Quindi prima di tutto, darò solo un'occhiata alla directory di lavoro così possiamo ricordare cosa c'era.

E poi farò "nextflow run hello-world.nf" e aggiungerò un singolo comando qui, "-resume".

Notate, singolo trattino, questo è davvero importante. Lo eseguirò e l'output sembrerà fondamentalmente esattamente lo stesso, con un paio di piccole differenze.

Notate qui dice "cached" in grigio. Ciò significa che Nextflow non ha eseguito l'attività. Questa volta ha trovato qualcosa che corrispondeva ai requisiti e ha riutilizzato direttamente quegli output piuttosto che rieseguire il passaggio.

E certamente, se guardate l'hash qui, potete vedere che questo corrisponde all'hash esistente che avevamo da un'esecuzione precedente.

## 3.3. Eliminare le directory di lavoro precedenti

Okay. Ma se state sviluppando in modo iterativo, accumulerete molti di questi file del workflow. Questo può essere un problema se potreste essere a corto di spazio.

Nextflow può aiutarci a pulire queste directory di lavoro con un paio di comandi di aiuto. Se faccio "nextflow log". Questo mi darà un elenco di tutte le diverse esecuzioni del workflow che ho fatto in questa directory, e hanno i nomi di esecuzione qui. Potete vedere quello gloomy quick, che è stato il primo che abbiamo eseguito, e poi questi due nuovi.

Possiamo ora prendere quel nome e usarlo con il comando "nextflow clean". Posso specificare un singolo nome di esecuzione. O ancora meglio, posso dire a Nextflow di eliminare tutto da prima di un singolo nome del workflow con "-before", e inserirò "stupefied_shaw". Quella è stata la mia esecuzione più recente, "-n".

Il comando "-n" ha detto a Nextflow di farlo come prova senza eliminare effettivamente nulla per davvero, e ci dice quali delle directory hash sarebbero state rimosse. Certamente, è solo quella della prima esecuzione. Entrambe le seconde esecuzioni usano la stessa directory hash.

Lo eseguirò di nuovo, ma ora invece di "-n" per prova, farò "-f" per forzare e ha rimosso quella directory hash. Ora se faccio "tree work", possiamo vedere, abbiamo solo questo file di output rimasto.

Ottimo. Quindi siamo riusciti a pulire un sacco di spazio su disco lì.

Un paio di cose da notare quando si eliminano le directory di lavoro, se fate symlink di cose alla vostra directory di risultati, quelle fonti di symlink saranno ora eliminate e i vostri risultati saranno persi per sempre. Quindi ecco perché usare la modalità copy è una cosa più sicura da fare, e generalmente ciò che raccomandiamo.

In secondo luogo, la funzionalità resume di Nextflow si basa su queste directory di lavoro. Quindi se le eliminate e eseguite di nuovo Nextflow, la funzionalità resume non funzionerà più. Quindi sta a voi tenere traccia di quali cose potreste avere bisogno o meno, ed eliminare le cose solo quando siete sicuri che sia sicuro farlo.

L'altra cosa che possiamo fare è che possiamo semplicemente eliminare l'intera directory di lavoro se abbiamo terminato l'esecuzione del nostro workflow e siamo sicuri di non averne più bisogno.

Quindi posso fare "rm -r work". So che non c'era nulla di importante lì. Ho i miei risultati che mi interessano nella directory results dove li abbiamo copiati. E quindi era sicuro eliminare la directory di lavoro. Sta a voi quale di questi approcci usate.

## 4. Usare un input variabile passato sulla riga di comando

Okay, cosa c'è dopo? Ho menzionato che avevamo codificato rigidamente alcuni dei valori nel nostro script del workflow qui, il file output.txt, e che potrebbe esserci un modo migliore per farlo.

Iniziamo con questo. Quello che faremo sono tre cose. Aggiungeremo un nuovo input al process. Diremo allo script del process come usare quell'input, e poi lo collegheremo nel workflow in modo da poterlo usare dinamicamente con un flag della riga di comando quando eseguiamo Nextflow.

Quindi prima di tutto. Aggiungiamo un blocco input qui. Proprio come output. Questa è una nuova sezione per il process, e dirò, "val greeting".

Notate qui, sto dicendo "val", che dice che questa è una variabile, non un path.

Posso poi scendere nello script e poi posso togliere questo testo codificato rigidamente qui e fare $greeting. Questo funziona proprio come qualsiasi altro linguaggio di programmazione. Stiamo definendo una variabile qui e la stiamo referenziando all'interno di questo blocco script. Quando Nextflow esegue questo process, la variabile sarà interpolata. E quando andiamo a guardare quel file .command.sh, vedremo la stringa effettiva codificata rigidamente qui invece.

## 4.1.3. Configurare un parametro CLI e fornirlo come input alla chiamata del process

Okay, ma dove forniamo la variabile? Successivamente andiamo giù alla sezione workflow, e potete vedere che l'estensione qui sta dicendo, ora ci aspettiamo un input, e mi ha dato un avviso.

Ora, la cosa più semplice che potremmo fare è semplicemente codificarlo rigidamente. Potrei scrivere "Hello World" e fornire quell'input stringa al process. Ma ancora, questo non risolverebbe davvero alcun problema. Dovremmo ancora tornare indietro e modificare il codice della pipeline ogni volta che vogliamo cambiare qualcosa, il che non va bene.

La buona notizia è che Nextflow ha un sistema integrato per gestire gli argomenti della riga di comando chiamati parametri. Quindi invece, posso usare una di queste variabili speciali chiamate params e posso chiamarla come voglio, ma dirò greeting in modo che corrisponda alla logica del workflow.

Premo salva e vediamo cosa possiamo fare con questo.

Quindi se torno al terminale. Quindi facciamo "nextflow run hello-world.nf". Proprio come prima, ma la differenza chiave è che facciamo --greeting

Notate, ci sono due trattini qui perché questo è un parametro. Quando abbiamo ripreso il workflow prima, era un singolo trattino. Questo perché resume è un'opzione core di Nextflow, e questo è un parametro che è specifico per la nostra pipeline.

Non confondete i due. È facile farlo. Se faceste --resume invece di un solo trattino, allora sarebbe "params.resume", che non farebbe nulla. Allo stesso modo, se faceste un singolo trattino qui, Nextflow non lo riconoscerebbe come argomento chiave.

Quindi è --greeting, che corrisponde a parameters greeting.

Posso ora seguire quello con qualsiasi testo voglio. Quindi sono in Svezia al momento, quindi dirò, "Hej världen".

Quindi eseguiamolo, vediamo cosa succede, momento della verità.

Okay, quindi potete vedere che il process è stato eseguito di nuovo, proprio come prima, sayHello con una singola esecuzione.

Questo avrà sovrascritto il file che era nella directory publishDir "results". E quindi fate attenzione quando rieseguite i file perché le cose nella directory pubblicata verranno sovrascritte.

Posso ora fare "code results/output.txt", e certamente, il nostro output è stato aggiornato e ora dice "Hej världen".

## 4.2. Usare valori predefiniti per i parametri della riga di comando

Okay, questo è ottimo. Ma il problema ora è che il nostro workflow si basa sul fatto che definiamo sempre questo parametro, ed è bello avere valori predefiniti sensati in modo che le cose vengano eseguite in modo sensato per il vostro workflow a meno che non sovrascriviate i valori predefiniti.

Quindi il modo in cui lo facciamo è impostando un valore predefinito per il parametro nel nostro script del workflow.

Quindi se torno al mio file hello-world.nf, posso andare nello script appena sopra workflow, digitare "params.greeting" e definirlo come qualsiasi altra variabile. Quindi mettiamo una stringa qui e diciamo "Holà mundo!"

Ora questo parametro ha un valore predefinito definito, che sarà usato qui, o possiamo ancora sovrascriverlo sulla riga di comando con --greeting, proprio come abbiamo fatto prima.

Quindi verifichiamo che funzioni. "nextflow run hello-world.nf"

Nessun argomento della riga di comando questa volta, e verifichiamo se ha fatto la cosa giusta.

"code results/output.txt". Ed eccolo. Abbiamo il nostro valore predefinito.

Okay, proviamo di nuovo, solo per verificare che non vi stia dicendo bugie. Eseguiamolo di nuovo, ma facciamo --greeting, e usiamo l'esempio dal materiale di formazione, diciamo "Konnichiwa!"

Riesegue il workflow, e certamente, il nostro file di output in alto è appena stato aggiornato con il nuovo valore che abbiamo fornito sulla riga di comando.

Ottimo. Questo è un aspetto davvero centrale nella scrittura di qualsiasi workflow Nextflow. Definire valori predefiniti sensati nel codice della pipeline, ma renderlo molto facile da configurare per l'utente finale avendo argomenti della riga di comando sul terminale.

Notate che l'utente finale può sovrascrivere la configurazione in più posti diversi. Potete avere un file di configurazione nella vostra directory home, che viene applicato a ogni singola esecuzione Nextflow che fate. Potete avere un file di configurazione in una directory di lancio. Potete avere un file di configurazione in una directory della pipeline. Tutte queste diverse posizioni di configurazione vengono caricate in un ordine specifico, che è descritto nei documenti Nextflow.

Okay, questa è la fine della sezione uno. Abbiamo avuto il nostro primo script del workflow in Nextflow con un process e un workflow. Abbiamo esaminato input, output, script e pubblicazione, e come collegare parametri e un canale di input nel nostro process.

Congratulazioni, il vostro primo passo verso la scrittura di codice Nextflow è completo.

Fate una piccola pausa e ci vediamo tra qualche minuto per il capitolo due.

[Trascrizione video successiva :octicons-arrow-right-24:](02_hello_channels.md)
