# Parte 1: Hello World - Trascrizione Video

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduzione assistita da IA - [scopri di più e suggerisci miglioramenti](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/tOukLxWCHiA?si=F0t9LFYLjAWoyRXj&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Note importanti"

    Questa pagina mostra solo la trascrizione. Per le istruzioni complete passo dopo passo, tornate al [materiale del corso](../01_hello_world.md).

    I numeri di sezione mostrati nella trascrizione sono forniti solo a scopo indicativo e potrebbero non includere tutti i numeri di sezione presenti nei materiali.

## Benvenuti

Ciao e bentornati.

Siete ora nella Parte Uno del corso "Hello Nextflow" chiamata "Hello World". In questo capitolo, inizieremo a costruire una comprensione delle basi fondamentali di Nextflow.

Quindi, si spera che ora siate configurati in Codespaces o in un ambiente equivalente con VS Code in esecuzione, e abbiate la vostra cartella Hello Nextflow nell'area di lavoro nell'Explorer con tutti questi diversi file qui.

Inizieremo facendo alcune cose molto basilari nel terminale usando Bash, e poi vedremo se riusciamo a fare le stesse cose all'interno di Nextflow così da farvi capire come appare la sintassi.

## 0. Riscaldamento

Quindi iniziamo davvero in modo semplice. Iniziamo semplicemente con "echo", per stampare qualcosa in un terminale. "Hello World". Premo invio e questo va al terminale. Hello World. Si spera che questo non sia una sorpresa per chi sta guardando questo corso.

Ok, facciamo qualcosa con questo. Invece di stamparlo semplicemente nel terminale, scriviamolo in un file. Premerò il tasto freccia su della tastiera, che scorre attraverso la cronologia di Bash, quindi mi dà il mio ultimo comando, e aggiungerò alla fine di esso, un piccolo simbolo maggiore di, che reindirizza l'output da questo comando a un file, e lo chiamerò output.txt.

Invio di nuovo, per eseguire quel comando, niente nel terminale questa volta, ma possiamo vedere sul lato sinistro che è apparso il nuovo file qui, chiamato output.txt.

Possiamo visualizzarlo in un terminale con qualcosa come cat. Quindi cat output.txt e sicuramente dice "Hello World". Possiamo anche fare doppio clic su di esso e si apre nell'editor di codice in VS Code.

## 1.1. Esaminiamo il codice

Va bene. Vi ho detto che era semplice. Cosa c'è dopo? Proviamo a prendere questo processo e rifarlo, ma questa volta, facciamolo all'interno di Nextflow.

Come ho detto, tutti i diversi capitoli di questo corso iniziano con uno script e questo si chiama Hello World. Quindi troverò Hello World. Lo visualizza in anteprima se ci clicco una volta, farò doppio clic per aprirlo nell'editor qui. E rimuoverò velocemente il terminale.

Ora questo è uno script molto semplice, quindi è semplice quanto possibile. È lungo solo 22 righe e fa sostanzialmente la stessa cosa. Infatti, parte di questo dovrebbe sembrare familiare. Possiamo vedere il nostro comando bash che reindirizza a un file lì.

Ok. Cos'altro? Inoltre, in questo file, possiamo iniziare a vedere alcuni dei concetti fondamentali di Nextflow. Abbiamo un processo in rosso qui e un flusso di lavoro. Queste sono le parole chiave speciali e la terminologia speciale in Nextflow.

## 1.1.1. La definizione del processo

Diversi processi all'interno di un flusso di lavoro racchiudono diverse unità logiche del vostro flusso di lavoro. Ogni processo fa una cosa.

Quando lo eseguiamo, genera un'attività o più attività, che sono i passaggi effettivi di esecuzione di una pipeline. Tutti i processi vengono poi orchestrati all'interno di un blocco workflow, che vediamo in basso, e in questo caso esegue solo quel singolo processo.

Il nome del processo segue questa parola chiave qui, e questo può essere praticamente qualsiasi cosa. E poi il contenuto del processo è all'interno di queste parentesi graffe.

C'è davvero solo un requisito per il processo, che è che includa una sorta di blocco script o exec. Questo è nelle triple virgolette qui, e questo è lo script bash che viene scritto nella directory di lavoro quando eseguiamo la pipeline ed è la cosa che effettivamente viene eseguita sul vostro computer o server.

Questo è tipicamente bash, ma potete anche inserire un diverso hash bang qui in alto, e potrebbe essere uno script Python o uno script R. Non importa. Qualunque cosa sia in questo script verrà eseguita.

C'è un'altra cosa che abbiamo aggiunto in questo processo qui, che è la dichiarazione di output. Questo dice a Nextflow che questo processo si aspetta un file di output chiamato output.txt. Dice che è un path, quindi dovrebbe essere gestito come un file, non diciamo, se questo fosse val, direbbe che è come una variabile o un valore.

Notate che questo non sta creando questo file. Non lo sta effettivamente generando. Questo viene fatto dallo script qui sotto. Sta solo dicendo a Nextflow di aspettarsi un file di output con questo nome file.

## 1.1.2. La definizione del flusso di lavoro

Ok. E poi in basso abbiamo un flusso di lavoro qui, e di nuovo, abbiamo una dichiarazione. Questo si chiama Main. Questo è l'equivalente del flusso di lavoro di un blocco script, se volete. È la parte del flusso di lavoro che fa qualcosa. E in questo caso, stiamo dicendo, chiama il processo chiamato sayHello.

Normalmente, ovviamente, la vostra pipeline sembrerà molto più complessa di questa. Probabilmente avrete più di un processo e userete i canali per orchestrare il flusso di dati tra di loro. Arriveremo a questo nelle prossime parti di questo corso, ma per ora, questo è sufficiente. Questa è una pipeline valida, che dovrebbe funzionare.

Posso anche cliccare su preview DAG qui in VS Code. Il DAG è una rappresentazione di una struttura di flusso di dati nella pipeline, e possiamo vederlo renderizzato sul lato come un diagramma mermaid. In questo caso è molto semplice. C'è una casella, che è il flusso di lavoro e un processo, che si chiama sayHello, ma potrebbe sembrare più interessante man mano che andiamo avanti.

## 1.2. Eseguiamo il flusso di lavoro

Ok, proviamo a eseguire questo flusso di lavoro e vediamo cosa succede.

Riporterò su il terminale in basso, cancellerò l'output, e digiterò Nextflow Run. E poi digiterò semplicemente il nome dello script, che è hello-world.nf. E premerò invio.

Ok, ha alcune cose standard in alto, che ci dicono che Nextflow è stato eseguito e quale versione era in esecuzione e qual era il nome dello script e tutto il resto.

E davvero la cosa importante che stiamo cercando qui è _qui_, che è un riepilogo delle diverse attività che sono state eseguite.

Se il vostro assomiglia a questo con un piccolo segno di spunta verde, allora complimenti. Avete appena eseguito la vostra prima pipeline. Fantastico.

Ci dice qui il nome del processo, che è stato eseguito, che si chiamava Say Hello, e ci ha detto che è stato eseguito una volta e che ha avuto successo. Questo si aggiorna man mano che andate avanti, quindi quando state eseguendo una pipeline più grande, vedrete i progressi rappresentati qui. Ma poiché questo è così piccolo, viene eseguito praticamente immediatamente.

## 1.2.2. Troviamo l'output e i registri nella directory di lavoro

Ora quando eseguite una pipeline Nextflow, ognuno di quei processi viene assemblato insieme, e ogni processo, come ho detto prima, può generare attività una o più. Quindi in questo caso, avevamo una singola attività da questo processo. È stata eseguita solo una volta e questo è stato fatto sotto questo _hash_ dell'attività.

Nextflow non gestisce i file nella vostra directory di lavoro direttamente, crea una cartella speciale chiamata work. E se faccio "ls", vedremo che è apparsa qui: _work_, e all'interno ci sono sottodirectory per ogni singola attività che viene eseguita. E questo corrisponde a questo hash. Quindi potete vedere se vado a "ls work/c4", e poi è troncato, ma inizia con 203, e quella è la directory di lavoro, che è stata creata da questo processo quando abbiamo eseguito la pipeline. E potete vederla anche sul lato.

Quando elenco quei file, potete vedere che il file output.txt è stato generato. Potete vederlo anche qui. E ci sono un mucchio di file nascosti, che non vengono mostrati con il mio "ls" normale.

Se clicco su output.txt, sicuramente, abbiamo il nostro output. Fantastico. Quindi la pipeline ha funzionato.

Potrebbe sembrare un sacco di codice standard per eseguire quello che era essenzialmente uno script bash di una riga, ma avrà più senso man mano che i nostri processi diventano più complicati. E questa directory work con Nextflow e questi file, che vengono creati è davvero la spina dorsale di ciò che rende Nextflow così potente.

Ogni attività, ogni elemento di una pipeline è isolato da ogni altra attività. È riproducibile. Non entrano in conflitto l'uno con l'altro, e tutto può essere eseguito in parallelo. È in realtà un modo davvero carino quando ci si abitua a causa di questo isolamento che potete entrare e vedere esattamente cosa è successo per una singola attività e fare debug.

Diamo una rapida occhiata a questi altri file nella directory di lavoro. Dall'alto verso il basso, abbiamo un file chiamato _.command.begin_. Questo è vuoto. È solo quello che viene chiamato un file sentinella, creato da Nextflow dicendo, ok, sto iniziando l'attività. Niente di interessante lì.

Poi c'è _.command.error_, _.command.log_ e _.command.out_. Questi sono tutti output dal comando bash o da questo script che è stato eseguito. Questo è standard error. Questo è standard out, e questo è i due di loro combinati come sono usciti. Quindi ottenete l'ordine logico.

Ok, quelli erano tutti vuoti anche per questo, quindi non molto interessanti, ma le cose diventano più interessanti quando arrivate a _.command.run_.

Questo è tipicamente uno script molto lungo. E questo è ciò che Nextflow effettivamente esegue. Se entrate qui, inizierete a vedere tutta la logica interna di Nextflow e vedere cosa sta facendo e come sta eseguendo il vostro processo. Questo dipenderà da dove state eseguendo, se stiamo eseguendo localmente o inviandolo come un job a SLURM, nel qual caso avremo intestazioni SLURM in alto. Tutte queste diverse configurazioni.

Generalmente, però, non avete davvero bisogno di guardare mai in questo file. È autogenerato da Nextflow e non c'è niente di veramente particolarmente unico per la vostra pipeline, che è in esso. Ma questo è davvero il nucleo di ciò che viene eseguito.

Il prossimo è molto più interessante. _.command.sh_ è lo script generato, che proviene dal vostro processo, e qui potete vedere che Nextflow ha aggiunto l'intestazione Bash, e poi ha eseguito il nostro comando, che era nel nostro blocco script.

E questo è tutto ciò che fa il file _.command.run_, esegue semplicemente questo file _.command.sh_.

Questo è davvero utile, ed è quello che di solito finite per guardare di più quando state cercando di fare debug di qualcosa e verificare che la logica della vostra pipeline Nextflow stia facendo ciò che vi aspettate che faccia.

Infine, abbiamo un file chiamato _.exitcode_, e questo cattura semplicemente il codice di uscita da un'attività, che in questo caso ha avuto successo. Quindi il codice di uscita era zero.

Se qualcosa va storto, finite la memoria o qualcos'altro e fallisce, allora questo è molto utile per capire cosa è andato storto.

## 1.3. Eseguiamo di nuovo il flusso di lavoro

Un'altra cosa da capire sulle directory di lavoro è che se continuo a eseguire questa pipeline ripetutamente, quindi se faccio _"nextflow run hello-world.nf"_, farà esattamente la stessa cosa, ma questa volta avrà un nuovo id di attività. Potete vedere che questo hash qui è diverso, e ora se guardo in work, ci sono due directory hash. E queste sono, di nuovo, separate l'una dall'altra.

Quindi ogni volta che eseguite un flusso di lavoro Nextflow, a meno che non usiate il resume, che usa la cache, toccheremo più tardi, rieseguirà quei processi in nuove directory di lavoro, che sono separate l'una dall'altra. Non avrete collisioni di nomi di file, non avrete problemi del genere. Tutto è isolato e pulito.

E se entriamo in questa directory, potete vedere tutti gli stessi file e lo stesso _output.txt_, che è stato ricreato da zero.

## 2. Pubblichiamo gli output

Ok, questo è ottimo per Nextflow per se stesso, mentre sta eseguendo la vostra pipeline in modo che tutte le cose siano separate l'una dall'altra e pulite e possano essere gestite.

Ma non è super utile se siete una persona che cerca di esplorare i vostri risultati. Non volete davvero scavare attraverso migliaia e migliaia di diverse directory di lavoro cercando di trovare i vostri file di risultato. E non siete davvero destinati a farlo. Le directory di lavoro non sono destinate ad essere lo stato finale di dove vengono creati i vostri file.

Lo facciamo pubblicando i nostri file.

## 2.1.1. Dichiariamo l'output del processo sayHello

Quindi se torno al nostro script, lavoreremo nel nostro blocco workflow qui. Diremo quali file aspettarci, quali file ci interessano, e poi creeremo un nuovo blocco sotto chiamato blocco output.

Questa è la nuova sintassi, che è arrivata con il parser di sintassi ed è predefinita nella versione 26.04 di Nextflow. Quindi se avete usato Nextflow un po' prima, questa è una delle cose che è nuova.

Quindi abbiamo il blocco main, e poi dirò publish e dirò a Nextflow cosa aspettarsi dalla pubblicazione. Lo chiameremo _first_output_, e lo chiameremo, _sayHello.out_.

Ho fatto accidentalmente un errore di battitura lì, ma questa è una buona opportunità per sottolineare anche alcune delle caratteristiche dell'estensione Nextflow VS Code. Potete vedere che subito mi ha dato una piccola linea rossa ondulata sotto questo dicendo che qualcosa non va. E se ci passo sopra, mi dirà che questa variabile non è definita. Non so cosa sia.

È abbastanza ovvio in questo caso, ho fatto un errore di battitura. Intendevo digitare, sayHello, e poi la linea ondulata scompare.

Ora è viola. Il parser di sintassi Nextflow sa che questo è un processo e quando ci passo sopra, mi dà una rappresentazione ridotta di come appare questo processo. Quindi posso vedere molto rapidamente a colpo d'occhio che non prende alcun input e ci dà questo output. Quindi lavorare in VS Code con questa estensione vi dà molte informazioni contestuali mentre state scrivendo codice.

Notate che possiamo riferirci all'output da questo processo con la sintassi _.out_. E al momento possiamo chiamarlo come vogliamo, è solo un nome di variabile arbitrario.

## 2.1.2. Aggiungiamo un blocco output: allo script

Dove diventa importante è quando facciamo il nostro nuovo blocco qui, e questo è sotto il blocco workflow ora, non siamo più all'interno del workflow. Parentesi graffe di nuovo. E questo è dove diciamo semplicemente a Nextflow dove mettere tutti i file, che vengono creati dal flusso di lavoro.

Ora prenderò questo nome di variabile, che ho creato qui, e lo metterò lì e metterò alcune parentesi graffe per questo. E dirò a Nextflow di usare un path. Ops. Path, tra virgolette. E userò il punto. Questo dice semplicemente a Nextflow di mettere il file nella radice della directory results. Quindi nessuna sottodirectory o altro.

Proviamo a eseguire di nuovo il nostro flusso di lavoro. Se faccio _"nextflow run hello-world.nf"_, allora si spera che dovrebbe sembrare praticamente esattamente lo stesso. Niente è davvero cambiato con Nextflow qui. Sta eseguendo le stesse cose. Le sta solo facendo di nuovo nelle directory di lavoro.

Ma ora se faccio _"ls results/"_, vedrete che c'è una nuova directory qui che è stata creata chiamata results, che è la directory di base predefinita per la pubblicazione del flusso di lavoro. E lì c'è un file chiamato _output.txt_.

Se faccio _"ls -l results"_, vedrete che questo è in realtà un soft link alla directory di lavoro. Quindi questo non è un file reale, è collegato alla directory di lavoro e ha raccolto tutti i file lì per noi.

## 2.2. Impostiamo una posizione personalizzata

"Results" è il nome predefinito per questo path. Se eseguo di nuovo il flusso di lavoro, e questa volta faccio _dash_ singolo trattino, questo è, perché è un'opzione core di Nextflow. _"-output-dir **my**results"._ Potrei anche fare solo _"-o"_ per abbreviare. Quindi imposterà una directory di base diversa per dove vengono archiviati i file e ancora una volta, qui in _myresults/_, ora abbiamo un _output.txt_.

È fantastico, ma probabilmente non vogliamo tutti i file solo nella radice. Vogliamo un po' di organizzazione, quindi possiamo anche creare una sottodirectory qui chiamata come vogliamo. Diciamo _"path 'hello_world'"_, e lo eseguo di nuovo. _"nextflow run hello-world.nf"_. Dovrebbe andare nella directory results in una sottodirectory e sicuramente, ora sotto results qui in alto abbiamo _hello_world/_ e abbiamo _output.txt_.

Cosa importante da notare, il vecchio file _output.txt_ è ancora lì. La directory results non viene cancellata quando fate questo. Solo i nuovi file vengono copiati lì. Sovrascriveranno i file che sono già lì se hanno lo stesso nome file, ma non cancelleranno quelli vecchi. Quindi dovete stare un po' attenti quando rieseguite le pipeline. Se non volete che siano sopra i file che sono già lì. Assicuratevi di usare una directory vuota e pulita.

## 2.3. Impostiamo la modalità di pubblicazione su copia

Ok, ho menzionato che questi file sono soft link, quindi se faccio _"ls -l results/hello_world/"_, potete vedere che sta facendo un soft link alla directory di lavoro. Questo è generalmente una buona cosa se state lavorando su qualcosa come HPC, e questi sono file davvero enormi e non volete duplicarli, perché significa che i file sono archiviati solo una volta sul file system.

Tuttavia, significa che se cancellate la directory di lavoro: se faccio _"rm -r work"_ e cancello tutti quei file intermedi che sono stati creati. Ora, se provo a leggere questo file _"results/hello_world/"_. Starà puntando come un soft link a un file che non esiste più e i dati sono andati persi per sempre e sono irrecuperabili, il che forse non è fantastico.

Quindi generalmente, dico che è buona pratica copiare i file invece di fare soft link se potete, perché è più sicuro. Siate solo consapevoli che userà il doppio dello spazio su disco a meno che non cancelliate quelle directory di lavoro.

Per farlo con il blocco output, andrò al primo output qui. Ho impostato il path prima e ora imposterò il mode e potete vedere mentre digito, l'estensione VS code sta, suggerendo cose che sa che è una direttiva di output qui. E dirò copy. Premo salva.

Rieseguiamo il flusso di lavoro. Creerà di nuovo i file, nuova directory di lavoro.

Ora, se vado a _"ls -l results/hello_world/"_ potete vedere che questo è un file reale e non è più un soft link, e Nextflow lo ha copiato. Buono a sapersi. Quindi path e mode sono cose che vi ritroverete a scrivere abbastanza spesso.

Ora, ovviamente, questo è molto semplice. Renderemo questo più complesso e potente man mano che andiamo avanti, e vedrete come rendere queste cose dinamiche e non troppo verbose.

## 2.4. Nota sulle direttive publishDir a livello di processo

Ora, ho detto quando abbiamo iniziato su questo, che questa è una forma di sintassi abbastanza nuova. È disponibile solo nelle ultime versioni di Nextflow mentre registro questo, e si chiama Workflow Outputs.

Se usate questo, è fantastico. Sblocca molte altre funzionalità interessanti all'interno di Nextflow, come, Nextflow Lineage per aiutare a tracciare l'eredità di questi file mentre vengono creati, e presto sarà il predefinito in 26.04. E in una data successiva in futuro, questo sarà l'unico modo per scrivere i vostri flussi di lavoro.

Tuttavia, poiché siamo in questa fase di transizione in questo momento, potreste benissimo vedere pipeline in circolazione, che usate qualcosa chiamato publishDir, che è il vecchio modo di farlo, e questo è definito non a livello di flusso di lavoro e output, ma questo è definito a livello di processo.

E questa dichiarazione dice sostanzialmente la stessa cosa. Dice, pubblica i file di risultato in una directory chiamata results, e usa una modalità copy. Quindi potete vedere che la sintassi è molto simile. Ma quando state scrivendo nuove pipeline ora, cercate di non usare questa direttiva publishDir, anche se la vedete, nei risultati AI o nella documentazione o in altre pipeline, perché questo è il vecchio modo di farlo.

Nel 2026 dovremmo tutti usare workflow outputs.

Questo è tutto documentato, se state facendo questo e avete usato Nextflow prima, potete andare ai documenti Nextflow qui, nextflow.io/docs/. E se scorro verso il basso fino ai tutorial, c'è un tutorial chiamato _Migrating to Workflow Outputs_.

È davvero buono. Passa attraverso tutta la sintassi, come è equivalente alla vecchia sintassi, perché l'abbiamo cambiata, e, ha una timeline e tutto. E passa attraverso tutti i diversi scenari con carichi e carichi di esempi. Quindi potete facilmente convertire il codice Nextflow esistente alla nuova sintassi.

## 3.1. Modifichiamo il processo sayHello per aspettarsi un input variabile

Ok, quindi abbiamo il nostro script semplice, che sta eseguendo un processo, creando un file, dicendo a Nextflow che è un output, e poi stiamo dicendo a Nextflow dove salvare quel file. È un buon inizio.

Ma sarebbe più interessante se non fosse tutto hardcoded. Quindi poi, pensiamo a come dire a Nextflow che questo processo può prendere un input variabile, che è qualcosa che possiamo controllare a runtime quando lanciamo un flusso di lavoro.

Dobbiamo fare alcune cose diverse per far sì che questo accada.

Innanzitutto, dobbiamo dire a questo processo che può accettare una variabile di input e digitiamo _input_ qui come un nuovo blocco di dichiarazione. E lo chiameremo _"val greeting"_.

Il bit val è l'equivalente di un path qui sotto. Dice a Nextflow che questa è una variabile, come una stringa in questo caso. E se ci passate sopra di nuovo, vi dice dall'estensione di cosa significa questo.

Poi diremo a Nextflow cosa fare con questo. Non è sufficiente dire semplicemente che c'è una variabile. Dovete dire nello script come usare quella variabile. E quindi mi libererò di questa stringa hardcoded qui, e metterò una variabile.

Lo farò velocemente senza parentesi graffe solo per mostrarvi che questo è, consentito, e questo è il vecchio modo di farlo. Ma ora con la nuova sintassi, raccomandiamo davvero di metterlo all'interno di parentesi graffe come questa, e rende davvero chiaro che questo viene interpolato da Nextflow qui.

Fantastico. Quindi _"input greeting"_ va in _$\{greeting\}._ L'ultima cosa è che dobbiamo dire a Nextflow a livello di flusso di lavoro che questo processo ora prende un input. E per farlo, gli daremo sostanzialmente una variabile.

## 3.2. Impostiamo un parametro da riga di comando per catturare l'input dell'utente

Potremmo hardcodarlo di nuovo, come Hello World, e funzionerebbe bene, ma ovviamente non ci dà davvero alcun vantaggio. Volevamo essere in grado di configurarlo a runtime, quindi vogliamo essere in grado di farlo sulla CLI, quando lanciate Nextflow.

E il modo in cui lo facciamo è un concetto speciale di Nextflow chiamato _params_. Lo chiameremo _params.input_.

Quello che fa è esporre questa variabile input sulla CLI ed è lì che usiamo un doppio trattino quando lanciamo Nextflow.

Posso chiamarlo come voglio, posso chiamarlo _hello, greeting_. Non importa. Qualunque cosa faccia lì sarà esposta come un'opzione CLI quando lanciamo una pipeline. E questo è un vero trucco magico di Nextflow perché significa che potete costruire il vostro script di flusso di lavoro molto rapidamente con questi parametri, e state essenzialmente costruendo una CLI personalizzata per la vostra pipeline, rendendola davvero facile da personalizzare diverse opzioni al volo quando lanciate.

Quindi. Proviamolo. Torniamo al nostro terminale. Abbiamo il nostro comando _"nextflow run"_ qui. E ora farò _"--input"_, che corrisponde al _"params.input"_ che abbiamo visto prima. Penso che nei documenti sia in francese. A Geraldine piace parlare francese. Lo farò in svedese perché vivo in Svezia. quindi dirò, "_Hej Världen_" e premerò invio.

Posso usare virgolette singole o doppie, influisce solo su come Bash lo interpreta.

Esegue la pipeline Nextflow esattamente nello stesso modo. Potete vedere la directory di lavoro e tutto è lo stesso. Ma ora se vado su _"results/hello_world/output"_. Possiamo vedere il nostro bel svedese qui invece.

Quindi abbiamo passato dinamicamente un input da una CLI a un parametro. Abbiamo passato quello come input al processo e il processo lo ha interpretato e lo ha messo in un blocco script, che ha poi cambiato dinamicamente l'output di quel risultato dello script. Abbastanza bello.

Logica abbastanza complessa con pochissima sintassi qui. E si spera che possiate vedere come questo ora inizia a scalare. E questo è come costruiamo davvero la logica e la personalizzabilità delle nostre pipeline nello script Nextflow.

## 3.4. Usiamo valori predefiniti per i parametri da riga di comando

Ok, è fantastico. Il problema però ora è che, ogni singola volta che eseguo questa pipeline, devo fare dash, input perché venga eseguita.

Se provo a eseguire senza questo parametro, ora Nextflow genererà un errore dicendo che aveva bisogno di questo parametro e non è stato impostato. e quindi non sapeva cosa fare.

Questa è una cosa nuova interessante, tra l'altro. In passato, Nextflow sarebbe stato eseguito semplicemente con una stringa vuota, e avreste avuto tutti i tipi di errori strani, che sarebbero stati difficili da capire. Ma nel nuovo parser di sintassi Nextflow, è un po' più attento e ve lo dice subito.

Quindi non vogliamo sempre specificare ogni singola opzione. È buona pratica specificare valori predefiniti sensati. Quindi come lo facciamo nel nostro script?

Noterete che quando abbiamo scritto questo, abbiamo messo _params.input_ direttamente dove lo stiamo usando. Quindi la soluzione ovvia è che definiamo un valore predefinito, e lo facciamo in alto nello script qui in un blocco params speciale nel flusso di lavoro. Questo è nello script del flusso di lavoro qui.

Di nuovo, un po' di nuova sintassi qui, quindi prestate attenzione. Questa è roba davvero interessante. Abbiamo il nome del parametro, che sarà previsto qui.

E poi dopo questo carattere due punti, stiamo definendo un tipo della variabile. Non dovete farlo, potete semplicemente lasciarlo vuoto, ma è davvero carino. Dice a Nextflow che ci aspettiamo una stringa e trattala come tale.

Se vogliamo un numero invece, per esempio, potremmo scrivere float, e questo direbbe che vogliamo un numero in virgola mobile. E se proviamo a eseguire con quello, allora genererà un errore. Se gli diamo una stringa, che non è un float. E lo passerà anche come tale. Come se facciamo string, allora sa che è una stringa. E anche se ha zeri iniziali ed è tutto numerico, lo passerà comunque come una stringa effettiva.

Quindi quella sicurezza del tipo è una caratteristica molto nuova di Nextflow, ma davvero potente per rendere il vostro codice più sicuro da scrivere e da eseguire.

Poi dopo quello abbiamo un simbolo uguale e poi il valore predefinito qui. Nextflow è stato scritto a Barcellona originariamente, quindi sembra appropriato che abbiamo un po' di, spagnolo qui, _"Holà mundo!"_ come predefinito.

Giusto salverò quello script, tornerò indietro, eseguirò di nuovo lo script senza _--input_. E questa volta dovrebbe essere eseguito e creerà il nostro nuovo file su in _results_. E in questo file ora dice _"Holà mundo!"_.

Questo è solo un valore predefinito però, quindi non significa che non possiamo ancora fare la stessa cosa di prima. Se torno indietro e trovo il mio vecchio script qui, _"Hej Världen"_, perché faccio _--input_ sulla riga di comando, questo sovrascriverà quel valore predefinito e lo userà di nuovo nel file output.txt.

Quindi questo nello script è solo il valore predefinito che sto impostando.

Man mano che costruiamo il nostro flusso di lavoro per essere più complesso e includere più parametri, questo blocco params in alto nello script inizierà a raccoglierli tutti in un unico posto.

E finite con questa simmetria abbastanza bella nel vostro script, dove avete effettivamente tutti i vostri input del flusso di lavoro qui e i vostri output del flusso di lavoro in basso. Ed è molto chiaro quale sia l'interfaccia del vostro flusso di lavoro verso il mondo esterno. Quindi potete prendere una nuova pipeline molto rapidamente con la nuova sintassi e capire come usarla.

Un'ultima cosa interessante. Non dobbiamo impostare un valore predefinito con questo. Se facciamo params input ma non impostiamo un valore predefinito, allora dice a Nextflow che questo parametro è richiesto, e di nuovo, la pipeline non riuscirà a essere eseguita senza di esso, ma vi darà un messaggio di errore più utile piuttosto che qualcosa sul fatto che sia null.

Quindi dice che ci aspettiamo che il suo input sia richiesto, ma non è stato specificato sulla riga di comando. Molto carino.

Ok, quindi si spera che ora sia chiaro come impostare la vostra pipeline Nextflow con input variabili e parametri, come impostare il valore predefinito, impostare, i tipi, potrebbe essere un flag booleano vero falso o un intero o diversi tipi qui. Come passarli nel vostro flusso di lavoro, dove passa, e poi interpola nel vostro processo. E poi sapete anche come personalizzarli sulla riga di comando quando lanciate Nextflow. Questo sta iniziando a sembrare più interessante del nostro semplice comando bash.

## 4. Gestiamo le esecuzioni del flusso di lavoro

Ok. Cosa c'è dopo? Per la parte finale di questo capitolo, parleremo un po' di come gestire tutte le diverse esecuzioni del flusso di lavoro. Se guardate nella mia barra laterale qui e nell'Explorer sotto work, vedrete che ho eseguito un mucchio di diverse pipeline e queste directory di lavoro stanno diventando abbastanza lunghe, ce ne sono molte.

E l'altra cosa è, come ho detto prima, ogni volta che rieseguo questa pipeline, sta creando un nuovo set di directory di lavoro, e sta rieseguendo tutti i processi da zero, il che è una buona cosa. Questo è il comportamento previsto. È riproducibile e sta rigenerando tutto fresco. Ma ovviamente, se state eseguendo processi che richiedono molto tempo, è fastidioso dover sempre iniziare la vostra pipeline dall'inizio se si è bloccata a metà, o se cambiate qualcosa alla fine della pipeline.

## 4.1. Rilanciamo un flusso di lavoro con -resume

Fortunatamente, Nextflow è davvero, bravo a sapere cosa è stato precedentemente eseguito e cosa è disponibile, e riutilizzare quei vecchi risultati è molto, semplice. Aggiungiamo semplicemente un, nuovo flag alla fine del comando _"-resume"_.

Ora, notate che ci sono due trattini su input perché è il parametro. C'è solo un trattino su resume perché è un'opzione core di Nextflow.

Fa inciampare le persone tutto il tempo, anche se avete usato Nextflow per molto tempo. Quindi ricordate sempre uno o due trattini. Dipende se è un'opzione core di Nextflow.

Ok, quindi ora faccio _-resume_ e eseguo esattamente lo stesso flusso di lavoro di nuovo. E questa volta dovrebbe sembrare praticamente esattamente lo stesso con una differenza chiave.

Nell'output qui, potete vedere che i risultati sono stati memorizzati nella cache. E infatti, questo hash dell'attività qui è esattamente lo stesso dell'esecuzione precedente, e ha semplicemente riutilizzato quella directory di lavoro nella sua interezza. Gli input e gli output e lo script erano tutti non modificati. E quindi prende semplicemente quel file da quello e se ci sono passaggi a valle nel processo, li passerebbe al passaggio successivo nella pipeline.

Quindi sta ancora eseguendo l'intera pipeline dall'inizio alla fine, ma sta usando risultati memorizzati nella cache per ognuna di quelle attività, dove può.

Ora, quando fate _-resume_, riprende semplicemente l'ultima esecuzione della pipeline nella vostra directory di lavoro, qualunque essa fosse. Ma potete effettivamente riprendere da qualsiasi esecuzione precedente che avete fatto lì. E ne abbiamo fatte parecchie ora.

## 4.2. Ispezioniamo il registro delle esecuzioni passate

Per guardarle tutte, possiamo fare _"nextflow log"_ invece di _"nextflow run"_, e questo ci darà un bell'output che mostra tutte queste diverse.. Ho bisogno di rendere il mio schermo un po' più piccolo così possiamo vederlo, tutte queste diverse esecuzioni quando le abbiamo fatte, l'id di sessione, il comando e tutto.

E possiamo guardare qui e possiamo prendere il nome di esecuzione di una qualsiasi di queste e poi riprendere una di quelle specifiche. Quindi posso tornare indietro e posso riprendere quella chiamata _hungry_ekeblad_. E la metto semplicemente dopo il _resume_.

Se siete curiosi, tra l'altro, tutti questi aggettivi e nomi di scienziati sono nel codice sorgente di Nextflow. È un modo davvero buono per ottenere la vostra prima pull request a Nextflow andando a trovarlo e aggiungendo il vostro scienziato preferito.

E comunque, quindi ho fatto quello e è tornato indietro e ha guardato i risultati memorizzati nella cache da questa esecuzione del flusso di lavoro, si è reso conto che poteva ancora riutilizzarli, e lo ha fatto. Quindi ho ottenuto di nuovo i risultati memorizzati nella cache.

## 4.3. Cancelliamo le directory di lavoro più vecchie

È fantastico. E se voglio pulire queste directory di lavoro? Ce ne sono carichi qui. Ci sono carichi di file. Forse so per certo che voglio riprendere dalle ultime due esecuzioni della pipeline, ma non mi interessano tutte quelle prima.

Allora posso sceglierne una qui e posso usare un altro comando Nextflow, che è _"nextflow clean"_, e posso fare _"nextflow clean"_, farò _"-before"_, e il particolare nome di esecuzione, che in questo caso era _reverent_pike_ e farò _"-n"_, che dice a Nextflow di fare solo un'esecuzione di prova. Quindi mi dice solo cosa cancellerà. Senza effettivamente fare nulla, quindi rimuoverebbe queste directory di lavoro.

Sembra sensato. Quindi farò lo stesso comando di nuovo, ma invece di _"-n"_ farò _"-f"_ per effettivamente fare la pulizia. E questa volta ha effettivamente rimosso tutte queste directory. E se entro e guardo le directory di lavoro, ora sembra molto più leggero. Fantastico.

Quindi questo è come pulire tutte le vostre directory di lavoro locali in un modo abbastanza sicuro senza distruggere completamente la cache. Quindi potete ancora riprendere se volete.

Se mai dimenticate quali sono questi flag per ogni comando Nextflow potete fare _"nextflow help"_, e poi il nome del comando. Quindi se faccio _"nextflow help clean"_, potete vedere tutte le diverse opzioni: _-after, -before, -but_, tutti i diversi modi per configurare questo comportamento di pulizia. Abbastanza bello.

## Takeaway

Ok, questa è la fine della parte uno di Hello Nextflow. È un inizio abbastanza intenso del corso, ma si spera che ora abbiate una comprensione abbastanza buona di come appare uno script Nextflow; con diverse parti chiave, i processi, i flussi di lavoro, gli output e i parametri. Sapete come configurarli con override di base dalla riga di comando, come creare un blocco di input dinamico con uno script dinamico e sapete come gestire tutte le vostre esecuzioni del carico di lavoro: vedere cosa avete già eseguito, riprendere, pulire. C'è un sacco di roba. Avete fatto molta strada. Quindi se volete fare una pausa e fare una breve passeggiata e una tazza di tè, ora è probabilmente un buon momento. Ve la siete meritata.

Da qui in poi, stiamo sostanzialmente costruendo su questa fondazione. Come possiamo rendere questo più complesso, più potente? Come possiamo renderlo più flessibile? Fare le cose che vogliamo fare la nostra analisi su scala.

## Quiz

Ora se scorrete verso il basso fino alla parte uno, hello world, sulla pagina web vedrete un piccolo quiz e questa è una cosa nuova che abbiamo fatto per questa versione della formazione Nextflow. E potete passare attraverso e mettervi alla prova per verificare di aver capito tutto il materiale che abbiamo fatto in questo capitolo.

Questo non viene inviato a noi o altro, è solo memorizzato nel vostro browser. Quindi non sappiamo quali sono le vostre risposte, ma è solo un piccolo autocontrollo per assicurarvi di non aver perso nulla o frainteso nulla. E potete provarlo quante volte volete.

Se siete come me, forse volete rimanere nel terminale nella vostra istanza VS Code, nel qual caso potete digitare il comando _quiz_ e poi dirgli semplicemente, quale capitolo state facendo. Quindi facciamo _"Hello World"_, e poi potete fare esattamente le stesse, domande del quiz, che sono nel browser web, ma solo nel vostro terminale.

Bello. Ok. Spero che vi piaccia. Divertitevi un po' e, ci vediamo nel prossimo capitolo tra un minuto per parlare tutto sui canali Nextflow.
