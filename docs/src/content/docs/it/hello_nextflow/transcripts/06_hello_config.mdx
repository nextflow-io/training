# Parte 6: Hello Config - Trascrizione Video

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduzione assistita da IA - [scopri di più e suggerisci miglioramenti](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/FcZTiE25TeA?si=tnXTi6mRkITY0zW_&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Note importanti"

    Questa pagina mostra solo la trascrizione. Per le istruzioni complete passo dopo passo, tornate al [materiale del corso](../06_hello_config.md).

    I numeri di sezione mostrati nella trascrizione sono forniti solo a scopo indicativo e potrebbero non includere tutti i numeri di sezione presenti nei materiali.

## Benvenuti

Ciao e bentornati alla Parte Sei di Hello Nextflow. Questa sezione è tutta dedicata ai file di configurazione ed è l'ultima parte di questo corso.

Nextflow è particolarmente efficace in due aspetti: riproducibilità e portabilità. I file di configurazione sono dove vediamo brillare davvero il secondo di questi aspetti. La capacità di configurare una pipeline Nextflow per essere eseguita in modi diversi e funzionare su sistemi differenti, senza dover modificare il codice della pipeline sottostante.

Questo superpotere permette alle pipeline Nextflow di essere riutilizzate da altre persone in luoghi diversi, o su diverse infrastrutture a cui potreste avere accesso voi stessi.

Significa che potete sviluppare il codice della pipeline sul vostro laptop, caricarlo sul cloud, eseguirlo sul vostro HPC, ed è lo stesso codice della pipeline e funziona ovunque.

In questa sezione, affronteremo alcuni argomenti. Inizieremo con come Nextflow gestisce i file di configurazione, da dove li carica, come si scrivono e come si strutturano, e quella separazione tra la pipeline stessa e ciò che dovrebbe andare in un file di configurazione.

Poi passeremo ad alcuni casi d'uso comuni come cambiare dove vengono memorizzati i file di output, e anche come far funzionare la pipeline su diverse infrastrutture, sia utilizzando diversi tipi di packaging software sia inviando job a diverse infrastrutture.

## Gerarchie dei file di configurazione

Okay, iniziamo. Quando si tratta di caricare file di configurazione, Nextflow può attingere da molti posti diversi, il che è una cosa positiva e può anche essere una cosa leggermente rischiosa perché a volte può essere un po' difficile sapere da dove sta ottenendo un file di configurazione e in che ordine carica le cose.

Quindi vi raccomando davvero di cliccare su questo link qui, che ci porta alla documentazione di Nextflow. E in questa pagina di configurazione, elenca i luoghi chiave da cui viene caricata la configurazione, e soprattutto, l'ordine in cui queste cose vengono caricate.

Quindi potete vedere, potete mettere un file di configurazione nella vostra directory home di Nextflow, che tipicamente è ".nextflow" nella vostra directory home. E quel file sarà sempre caricato da ogni esecuzione di Nextflow sul vostro sistema.

Il prossimo posto dove guardare è un file nella radice del vostro repository o directory della pipeline chiamato "nextflow.config".

Poi dopo quello, un altro file chiamato "nextflow.config", ma questa volta nella directory da cui state lanciando Nextflow: la directory di lancio.

Infine, potete fornire percorsi di file di configurazione sulla riga di comando con un argomento "-c", e potete farlo più volte. E vengono applicati nell'ordine in cui li specificate.

Potete fornire file di configurazione in tutti questi luoghi se volete, e saranno caricati iterativamente, ciascuno sovrascrivendo il precedente solo negli scope di configurazione dove si sovrappongono.

Questo è un sistema davvero potente perché significa che potete impostare valori predefiniti sensati e poi diventare gradualmente sempre più specifici man mano che vi concentrate su quella configurazione.

## 0. Riscaldamento: Eseguire hello-config.nf

Okay, chiudiamo questo e saltiamo nei nostri Codespaces e iniziamo. Come prima ho fatto pulizia qui, ho rimosso le mie precedenti directory dei risultati, le mie directory Nextflow e work e così via. Non preoccupatevi se avete ancora quei file in giro. È solo perché sono molto ingrandito e quindi le cose diventano disordinate molto rapidamente altrimenti.

Lavoreremo con hello-config.nf, l'ultimo file nella nostra directory, e questo dovrebbe seguire da dove ci eravamo fermati nella sezione precedente.

Quindi abbiamo i nostri quattro diversi processi, che vengono inclusi dai file modulo. Abbiamo i nostri parametri della pipeline, il nostro blocco workflow dove stiamo chiamando i diversi processi e collegando insieme i canali, pubblicando i canali di output, e poi il blocco output in fondo dove definiamo dove quei file dovrebbero essere memorizzati e come dovrebbero essere copiati.

Abbiamo già anche un file "nextflow.config" dal capitolo precedente, dove abilitiamo Docker, e costruiremo su questo file oggi.

Come prima, abbiamo cambiato il percorso di output in questo script principale in hello config, solo così non si scontra con i risultati precedenti che avete generato.

Okay, controlliamo velocemente che tutto funzioni ancora come ci aspettiamo. Apriamo un terminale e facciamo nextflow run hello-config.nf. Nextflow si carica. Dovrebbe eseguire i nostri quattro diversi processi. Generare della bella arte ASCII usando cowpy e poi salvare i nostri risultati nei file dei risultati in quella directory.

Posso dare una rapida occhiata qui solo per assicurarmi che questi file sembrino come ci aspettiamo, e infatti, c'è il nostro gigantesco tacchino. Ottimo.

## 1.1. Spostare i valori predefiniti in nextflow.config

Ora la prima cosa che faremo è iniziare a spostare alcune cose dal nostro script nel nostro file di configurazione.

E ciò che ci interessa sono principalmente i parametri in questa fase. Vogliamo prendere i valori predefiniti e metterli nel file di configurazione, così è più chiaro quali sono i valori predefiniti ed è più facile per le persone sovrascriverli.

Prenderò questo blocco params qui dallo script e lo metterò nel file di configurazione. E dobbiamo stare un po' attenti qui, perché in questo momento la sintassi è leggermente diversa tra configurazione e script. Il file di configurazione non può accettare dichiarazioni di tipo perché non stiamo davvero definendo questi parametri, li stiamo solo referenziando. Quindi eliminerò quelli.

Ma per il resto è molto simile. Abbiamo un blocco params e poi abbiamo i nostri diversi parametri di input, parametro batch, parametro character.

Posso ora tornare al mio script e non ho più bisogno di definire questi valori predefiniti perché questi valori sono ora nel mio file Nextflow config.

Tuttavia, lascio i nomi dei parametri e i loro tipi, così che Nextflow conosca quelle informazioni e possa ancora fare tutta la sicurezza dei tipi e tutto il resto.

Okay. Salviamo quei file e controlliamo velocemente che tutto funzioni ancora come prima. Non dovrebbero esserci cambiamenti qui. Abbiamo mantenuto gli stessi valori. Abbiamo solo spostato dove sono stati definiti.

Ottimo.

## 1.2. Usare un file di configurazione specifico per l'esecuzione

Ora, finora abbiamo lanciato Nextflow dalla stessa directory dove abbiamo il nostro script della pipeline. Quindi la nostra directory di lancio e la nostra directory della pipeline sono più o meno la stessa cosa.

Per mostrare come possiamo avere diversi file di configurazione con diverse directory di lancio, creeremo ora una nuova sottodirectory.

Quindi dirò mkdir, e la chiameremo tux-run.

E poi farò cd, cambierò directory in tux-run. E notate che ora siamo nella nostra directory di lavoro che non è più nella stessa directory degli script della pipeline.

Okay, creiamo un nuovo file "nextflow.config". Quindi touch nextflow config, e apriamolo semplicemente in VS Code. Potete vedere anche nella barra laterale qui che siamo ora in questa sottodirectory.

Ora possiamo prendere lo stesso blocco params che avevamo nel nextflow.config di livello superiore, copiarlo qui e ora possiamo cambiare questi valori.

Innanzitutto, i dati sono ora un percorso relativo diverso perché siamo in una sottodirectory, quindi dobbiamo aggiornarlo. E poi cambieremo batch in experiment, e cambieremo il character da Turkey a tux.

Ora clicchiamo salva lì, e proviamolo. Proprio come con data, devo ora dire ../ per arrivare allo script. Quindi è Hello config. E premo invio.

Il codice della pipeline non è cambiato affatto, ma ora avremo due set di configurazione che si caricano, e il file di configurazione della directory di lancio dovrebbe sovrascrivere i valori predefiniti, che erano impostati nel nextflow.config della pipeline, e dovremmo ottenere diversi set di risultati.

Infatti, all'interno della nostra directory qui, all'interno di tux-run, potete vedere che abbiamo una directory dot Nextflow e una directory work e questo perché queste vengono create sempre nella vostra directory di lancio. Quindi queste sono diverse dalle directory work e results che avevamo dalle esecuzioni precedenti.

Ora, se guardo in results, possiamo vedere il nostro collected e c'è il nostro piccolo personaggio tux. Quindi potete vedere che quei parametri sono stati interpretati correttamente.

## 1.3. Usare un file di parametri

Okay. Prima quando parlavo dei diversi file di configurazione che potevano essere caricati, ho tralasciato un altro posto da cui possiamo ottenere la configurazione.

Potete ottenerla dalla riga di comando come abbiamo visto con dash dash nomi dei parametri, ma possiamo anche fornire un file YAML o JSON, solo di parametri.

Il file di configurazione può avere tutti i diversi tipi di scope, ma questi file sono solo parametri, ed è un modo user-friendly per fornire molti parametri contemporaneamente, e forse un modo un po' più riproducibile perché li scrivete su file, quindi è facile ottenerli in una fase successiva.

Quindi torniamo al nostro terminale e prima che ce ne dimentichiamo, assicuriamoci di tornare su di una directory, così non sono più nella sottodirectory, e guarderò il file YAML che abbiamo qui chiamato test-params.yaml.

Quindi se faccio semplicemente code test-params.yaml, potete vedere che questo è solo un normale file YAML. Niente di speciale. Con le chiavi che sono i nostri nomi di parametri, con la formattazione YAML quindi due punti qui, e poi un valore.

Notate che questo non è codice Nextflow, quindi non possiamo mettere cose come variabili qui. Questi sono solo valori statici.

Inoltre, poiché JSON in realtà viene analizzato come YAML, possiamo anche avere un file test-params.json, che sembra molto simile. È solo un formato dati diverso.

Quindi abbiamo due diversi file di test qui e abbiamo variabili leggermente diverse.

Okay, quindi come li diamo a Nextflow? È molto semplice. Facciamo Nextflow run hello config, come prima. E invece di "-c" per file di configurazione, o caricare quei nomi di file predefiniti, facciamo -params-file. Singolo trattino perché è un'opzione core di Nextflow.

E poi passiamo il percorso per quel file. Quindi farò "-params-file test-params.yaml", e vedremo se quelli vengono caricati correttamente.

Okay. È stato eseguito. Ricordiamoci cosa c'era in questo file YAML. Quindi il batch era impostato su YAML, quindi è così che dovrebbe essere chiamato, e dovrebbe avere uno stegosauro. Quindi andiamo su e guardiamo in results. E abbiamo COLLECTED-yaml. Quindi vediamo se abbiamo uno Stegosauro. Fantastico, uno Stegosauro che indossa un cappello. È quello che ci piace.

Quindi ha funzionato davvero bene, ed è esattamente lo stesso con il file JSON. Cambiamo semplicemente l'estensione del file qui e Nextflow sa come leggerlo.

E in questo caso, dovremmo avere un batch chiamato JSON e dovremmo avere una tartaruga. Diamo un'occhiata. Meraviglioso. Uno dei miei strumenti CLI preferiti.

## 2.1. Personalizzare la directory di output con -output-dir

Okay, quindi questo è stato principalmente pensare agli input della pipeline e cambiare i parametri. E gli output?

Ora, anche se abbiamo cambiato le sottodirectory usando params, potreste aver notato che tutti i nostri file stanno ancora andando in results.

Possiamo cambiare quella directory base in cui tutti i file vengono pubblicati con un flag della riga di comando chiamato -output-dir. Quindi se faccio Nextflow run hello config, e poi faccio -output-dir, e lo chiameremo "custom-outdir-cli". Non riesco a digitare. Solo così ricordiamo da dove vengono questi file.

Questa è un'opzione core di Nextflow ed è molto nuova. È stata aggiunta solo di recente, e questa è una delle cose che possiamo fare con il nuovo parser del linguaggio e tutto il resto.

È un po' lungo da digitare. Potete anche chiamarlo semplicemente "-o" se volete. Quindi se torno indietro. Posso semplicemente accorciarlo a "-o", che è un po' più semplice.

Okay. Eseguiamo quello. Non abbiamo cambiato nulla nella nostra pipeline o anche nella nostra configurazione a questo punto, e dovrebbe si spera salvare tutti i nostri risultati in una diversa directory di livello superiore. E potete immaginare di poter impostare questo praticamente a qualsiasi percorso vogliate.

È appena arrivato in alto. Abbiamo un custom-outdir-cli, e tutti i file sono organizzati lì esattamente nello stesso modo, con le stesse sottodirectory e nomi di file. Quindi questo è un modo davvero facile per cambiare semplicemente dove la pipeline pubblica i suoi risultati, senza pensare troppo a come quei risultati sono organizzati.

## 2.1.2. Rimuovere i percorsi hardcoded dal blocco output

Se guardo in questa directory, possiamo vedere che abbiamo ancora una sottodirectory chiamata Hello Config, che sembra un po' ridondante ora.

Quindi carichiamo di nuovo il nostro script e possiamo ora rimuovere quella sottodirectory dal blocco output in fondo. Perché non ne abbiamo davvero più bisogno. Quindi possiamo farlo ora, eliminarlo da qui. E poi se è solo questo, potete eliminarlo completamente o lasciarlo come stringa vuota. Lo lascerò come stringa vuota per ora, perché torneremo e metteremo alcune cose diverse al suo posto in futuro. Ma se non vi importa delle sottodirectory, è più pulito rimuovere completamente la dichiarazione path lì.

Okay, salviamo. Proviamolo velocemente di nuovo. In realtà rimuoverò la mia directory "custom-outdir-cli" così non ci confondiamo con file esistenti lì. Perché ricordate, quando pubblicate cose, non rimuove i file che c'erano già. Aggiunge solo nuovi. Eseguiamo di nuovo quel comando, custom-outdir-cli.

E ora se fate "ls custom-outdir-cli", non c'è più directory lì chiamata Hello Config.

## 2.2.1. Impostare outputDir nel file di configurazione

Okay, il flag della riga di comando qui, "-o" o "-output-dir" è buono. Ma che ne dite di impostare i valori predefiniti per questo nella configurazione? Come lo facciamo?

Apro il file "nextflow.config", chiudo tutto il resto e mi sbarazzo di quello. Possiamo aggiungere una nuova opzione di configurazione qui, che ho appena copiato dal sito web del materiale di formazione, e si chiama outputDir.

Non è sotto nessuno scope. Non è sotto params o altro. È di livello superiore, e possiamo impostarlo su una stringa. Ora una cosa semplice da fare è semplicemente cambiarlo in qualsiasi cosa diversa da results come stringa hardcoded. Ma poiché questo è in un file di configurazione Nextflow, possiamo essere un po' più intelligenti qui e includere anche variabili.

E potete vedere qui che abbiamo incluso una variabile params, params.batch, che fa parte di questa stringa. Questo significa che possiamo riutilizzare variabili che provengono da altri posti. E in questo caso, se facciamo --batch, quando eseguiamo Nextflow Pipeline, avremo una sottodirectory nel nostro percorso personalizzato basata su quale era il nome del batch.

Okay, quindi proviamo questo e diamo una rapida occhiata per vedere come appaiono i risultati. Quindi se faccio Nextflow run hello config e --batch my_run. Ricordiamoci come appariva la configurazione. Quindi è custom-outdir-config.

Tree custom-outdir-config. E potete vedere che il batch era chiamato my_run. E poi abbiamo quella sottodirectory chiamata my_run. Quindi quel percorso di file dinamico ha funzionato.

E non solo, non è più andato in una directory results predefinita, e non ho dovuto specificare nulla sulla riga di comando per cambiare la directory base. Quindi abbiamo reimpostato con successo il valore predefinito per l'outputDir predefinito.

## 2.2.2. Sottodirectory con nomi di batch e processi

Okay, portiamo questo un po' più avanti. Questa è una variabile dinamica all'interno del file di configurazione. E lo script stesso? Ora, finora abbiamo avuto questi percorsi qui e anche questi possono essere dinamici. Quindi invece di hardcodare semplicemente qualcosa, possiamo mettere delle parentesi graffe e mettere qualcosa di dinamico.

Quindi per esempio, abbiamo i nostri processi chiamati sayHello. Potremmo fare sayHello.name, che è un attributo del processo, che è un po' noioso perché è solo "sayHello" in questo caso. Ma è variabile.

Quindi questo vi dà un'idea. Quindi possiamo metterlo qui e dire convertToUpper.name, collectGreetings.name, collectGreetings.name di nuovo, e cowpy.

Ora quando eseguiamo, la directory base sarà ancora custom-outdir-config. E sarà in una sottodirectory chiamata params.batch, ma le sottodirectory sotto quella dovrebbero essere organizzate per nome di processo.

Proviamolo e vediamo se funziona. Quindi rimuoverò la directory precedente così non ci confondiamo, e userò esattamente lo stesso comando Nextflow Run.

Dovrebbe essere eseguito nello stesso modo. Potrei usare dash resume su tutti questi per renderlo un po' più veloce e usare i risultati calcolati in precedenza. Ora, se faccio tree custom-outdir-config, potete vedere che non è in results, è nella nostra directory base con il nome del batch. E potete vedere che tutti i risultati sono ora organizzati all'interno di sottodirectory nominate dopo il processo. Quindi abbiamo due posti diversi dove stiamo definendo percorsi di output dinamici qui.

Okay. Ultima cosa, aggiungiamo di nuovo quelle cartelle intermedie, che avevamo prima perché erano piuttosto carine. Intermediates.

E possiamo anche pensare un po' a questo params.batch, forse come sviluppatore della pipeline mi piaceva davvero averlo nella sottodirectory, ma se gli utenti finali della pipeline stanno impostando "-o" o -output-dir sulla CLI, sta completamente sovrascrivendo questa intera dichiarazione, e perdiamo quella sottodirectory.

Quindi quello che possiamo fare è prendere quel percorso dinamico fuori dall'outputDir config, che verrebbe sovrascritto, e metterlo nel percorso di output, che non viene sovrascritto.

Quindi possiamo fare params.batch slash intermediates slash sayHello.name, e fare tutto questo in una stringa tra virgolette doppie, così viene interpolata da Nextflow.

Posso ora copiare, ops. Copiare questi giù agli altri processi. Ricordate di metterli tutti tra virgolette. E rimuovere intermediates da questi output particolari.

Okay? Sembra leggermente più complesso ora, ma potete vedere che stiamo davvero iniziando a costruire una bella struttura di directory di output organizzata nel nostro codice.

E ciò che è davvero bello è che questa complessità extra nel codice non passa attraverso la CLI. Quindi possiamo eseguire il nostro comando con -output-dir e qualsiasi variabile batch, pensando solo a come eseguire la pipeline e non pensando troppo a cosa c'è nel codice. E i nostri file di output saranno costruiti davvero bene in un modo molto ben organizzato, il che è bello per le persone che usano la pipeline fondamentalmente.

Ottimo. Mentre scrivo questo, mi rendo conto di aver fatto un errore. Vediamo se qualcuno mi ha beccato qui. Abbiamo collectGreetings.name, quindi qualcosa è andato un po' storto. E sì, infatti, ho accidentalmente dimenticato di mettere questi tra parentesi graffe.

Quindi ricordate, fate attenzione quando scrivete il vostro codice e assicuratevi di dire a Nextflow cosa è una variabile e cosa è solo una stringa. Perché farà esattamente quello che gli dite di fare. E niente di più. Come tutti i buoni computer. Okay, questo dovrebbe sistemarlo.

## 2.3. Impostare la modalità di pubblicazione a livello di workflow

C'è un pezzo di questo script, che ancora non amo, che è il fatto che stiamo scrivendo mode copy ancora e ancora, e se c'è una cosa che non ci piace, è ripeterci.

Quindi possiamo pulire questo un po' prendendo questo e spostandolo nella configurazione. E infatti, possiamo impostarlo per l'intera pipeline in una volta sola. Quindi non dobbiamo dirlo più volte.

Andiamo al nostro file di configurazione e abbiamo un nuovo scope qui chiamato workflow. E possiamo fare parentesi graffe o possiamo fare la notazione con il punto. Non fa alcuna differenza. Il sito web del materiale di formazione usa la notazione con il punto. Posso dire output e possiamo mescolare e abbinare, quindi mode equals copy. Ottimo.

E ora possiamo tornare qui e eliminare questi. Ora potremmo lasciarli in posizione. La configurazione sta fondamentalmente sovrascrivendo ciò che è scritto qui, ma dato che l'abbiamo nella configurazione a livello di pipeline, e questi due file vengono forniti insieme, non c'è davvero motivo di farlo due volte.

Okay. Solo un controllo di sanità, perché apparentemente facciamo errori. Eseguiamo di nuovo e controlliamo solo che stiamo usando correttamente la modalità copy per pubblicare i file. Quindi eseguiremo di nuovo lo script e questa volta abbiamo messo i risultati in una directory chiamata config-output-mode, vediamo come appaiono i file lì.

E poi se faccio "ls -l" per guardare batch, e possiamo guardare cowpy, per esempio. E dovremmo vedere, sì, che questo è un file vero qui, che non è un soft link, quindi quell'attributo di configurazione è stato applicato correttamente.

## 3. Selezionare una tecnologia di packaging software

Okay. Finora ci siamo concentrati sugli input e sugli output, i file con cui il flusso di lavoro sta lavorando. Ma che ne dite dell'infrastruttura? Ho detto all'inizio che Nextflow vi permette di eseguire la stessa pipeline su diverse configurazioni di calcolo. Quindi come appare?

Per mostrare questo, passeremo dall'uso di Docker per eseguire cowpy, e invece useremo Conda per fare la stessa cosa.

Posso farlo molto semplicemente. Se vado a code, "nextflow.config". Se ricordate in alto, abbiamo definito docker.enabled in precedenza, e nell'ultimo capitolo così che potessimo usare il container con cowpy dentro.

Dirò a Nextflow di non usare Docker. Imposterò quello su false. E dirò Conda enabled equals true. Quindi dire a Nextflow, per favore usa Conda.

Ora abilitare solo Conda non è sufficiente di per sé. Proprio come abbiamo fatto con Docker, dobbiamo dire a Nextflow dove può ottenere il software di cui ha bisogno.

Quindi se saltiamo nei moduli qui. E apriamo lo script cowpy. Possiamo vedere che abbiamo una dichiarazione container in alto. E il container viene usato da Docker, ma anche Singularity, Apptainer, e molti degli altri strumenti software.

Ma non può essere usato per Conda, quindi abbiamo una dichiarazione separata chiamata "conda", e potremmo semplicemente scrivere "cowpy". E questo lascerà alla risoluzione del pacchetto conda capire il modo migliore per risolvere, quello secondo il vostro ambiente conda locale.

Oppure è buona pratica fare quello che dice il sito web del materiale di formazione, che è definire uno specifico canale conda con la sua notazione a doppi due punti, e sicuramente definire una versione specifica del software così che ogni persona che esegue la pipeline otterrà la stessa versione.

Notate che i container sono un po' superiori in questo aspetto, perché quando installate qualcosa con Conda, risolverà comunque tutte le dipendenze per quel pacchetto, e possono cambiare nel tempo. Chiamato dependency drift.

Quindi i container, tuttavia, bloccano l'intero stack di un'intera pila di dipendenze software fino in fondo, quindi potete essere un po' più sicuri che A, funzionerà, e B, sarà riproducibile.

Quindi se siete in grado di usare Docker o Singularity o Apptainer, lo raccomanderei decisamente.

Ora ciò che è bello di questo è che il file modulo, che è scritto dallo sviluppatore della pipeline, ora ha sia Container che Conda, e quindi stiamo dicendo alla persona che sta eseguendo questa pipeline, non ci importa quale soluzione di packaging software usate. Funzionerà sia con Docker che con Conda, e questo è dove ottenere il software in entrambi i casi.

Possiamo aprire il terminale e proviamolo. Quindi Nextflow run hello config --batch conda. E la prima volta che questo viene eseguito con conda, sarà un po' lento quando arriva a quel particolare processo, perché deve eseguire "conda install".

E sta creando un ambiente conda speciale solo per questo processo. Quindi non sta usando il mio ambiente conda globale, che ho sul mio terminale. Sta creandone uno solo per quel processo. Questo è buono perché evita cose come conflitti di dipendenze tra diversi processi nel vostro flusso di lavoro. Se i vostri processi hanno strumenti che necessitano di diverse versioni di Python o cose del genere, va bene perché stanno usando diversi ambienti conda.

Nextflow memorizza nella cache questi ambienti conda localmente, potete vedere che vi dice dove si trova quel percorso, è nella directory work qui. E quindi la prossima volta che eseguo questo script con Conda, sarà molto più veloce perché troverà quell'ambiente conda esistente e lo riutilizzerà semplicemente. Ma la prima volta che lo facciamo, deve andare a recuperarlo, risolverlo, scaricare tutte le dipendenze e configurare tutto.

Okay, ottimo, è stato eseguito. Possiamo semplicemente ricordarci cosa è attualmente configurata la pipeline per usare. Se guardiamo nel file di configurazione, era "custom-outdir-config" in questo momento per me. Vediamo se vado a quella directory base. E ho fatto --batch conda. C'è la nostra sottodirectory conda. Quindi ha funzionato e c'è il nostro output cowpy.

Quindi ha recuperato cowpy, l'ha installato sul mio sistema locale usando conda, e ha eseguito il processo. E ciò che è fantastico è che, come utente finale, non ho dovuto pensare affatto a nessuna gestione del software lì. Nextflow l'ha semplicemente sistemato per me. Ho detto, devo usare conda su questo sistema. Lo sviluppatore della pipeline ha detto quali pacchetti mi servivano. E Nextflow ha fatto il resto. Molto potente.

Notate che potete effettivamente usare una miscela di diverse tecnologie. Quindi posso abilitare Docker per processi specifici, e conda per altri processi, o dire che alcuni processi dovrebbero usare semplicemente qualsiasi software locale avevo installato. Questo è piuttosto insolito, ma è possibile, e in alcuni casi, per esempio, se state usando certi software che potrebbero essere difficili da pacchettizzare in Docker, avete una via d'uscita.

## 4. Selezionare una piattaforma di esecuzione

Quindi questo è il packaging software. L'altra parte della portabilità su altri sistemi è dove i job vengono effettivamente eseguiti. Al momento, sto eseguendo fondamentalmente sul mio laptop o in questi Codespaces, che è un singolo computer. Non c'è niente di speciale. Nextflow sta essendo un po' intelligente nel parallelizzare i job nel miglior modo possibile, ma è tutto su un sistema.

Ora, se state eseguendo su un HPC, probabilmente avete qualche tipo di scheduler di job come SLURM o PBS o qualcosa del genere, e invierete job a quello scheduler e distribuirà tutti i job a diversi nodi di calcolo.

Un altro modo di eseguire è sul cloud. Quindi forse state usando AWS Batch, o Azure Cloud, o Google. E tutti questi funzionano in un sistema simile dove avete uno scheduler e inviate job e vengono inviati a posti diversi per essere calcolati.

Ora nel lontano passato quando ho iniziato a fare bioinformatica, il software di tutti per eseguire analisi era molto legato alla loro infrastruttura computazionale, il che rendeva quasi impossibile replicare.

Ma con questa separazione della configurazione in Nextflow, e con la capacità di Nextflow di interagire con molti diversi backend di infrastruttura di calcolo, è molto semplice prendere la nostra pipeline senza modificare affatto il codice della pipeline e semplicemente cambiarlo.

## 4.1. Puntare a un backend diverso

Quindi se andiamo al nostro file "nextflow.config", e possiamo ora mettere della configurazione a livello di processo. Quindi se metto in alto lo scope process e posso impostare l'executor, e qui è impostato su local, che è il predefinito.

Notate che poiché questo è a livello di processo, possiamo indirizzare cose a diversi processi. E quindi potete effettivamente configurare gli executor per essere specifici del processo e avere un'esecuzione ibrida, dove alcuni job potrebbero essere eseguiti localmente, ovunque il job Nextflow venga eseguito. Alcuni vengono inviati a diversi HPC e alcuni potrebbero essere inviati al cloud. Potete essere intelligenti quanto volete.

Ora, è molto difficile fare una demo in un ambiente di formazione come questo perché non ho un HPC a cui inviare. Ma quello che posso fare è se digito slurm, possiamo barare un po' e potete farvi un'idea di questo.

E questo è davvero più interessante per le persone che sono abituate a eseguire su SLURM e sanno come appaiono gli header SLURM. Ma se faccio Nextflow run, hello config. Fallirà perché proverà a inviare job a un cluster che non esiste. Quindi otterremo qualche tipo di errore su sbatch che non è disponibile.

Sì, scritto. Questo è lo strumento. Questo è lo strumento CLI che usate per inviare job a un cluster slurm. Ma quello che possiamo fare è andare a guardare nella nostra directory work qui facendo command click, aprire quella directory e guardare il .command.run. E potete vedere in cima al file .command.run, abbiamo i nostri header sbatch, che dicono a un cluster SLURM teorico come gestire questa sottomissione di job.

E quindi potete vedere che Nextflow sta essendo intelligente, sta facendo tutte le cose giuste. È solo che non avevamo un cluster a cui inviare.

## 5. Controllare le allocazioni di risorse di calcolo

Cos'altro è diverso tra diverse infrastrutture di calcolo? Un'altra cosa è quante risorse disponibili avete, e infatti, in molti ambienti di calcolo, è un requisito che dobbiate specificare quante CPU e quanta memoria un job necessita.

Di nuovo, Nextflow astrae questo per noi, così che non sia più specifico per un singolo tipo di ambiente di calcolo, e possiamo digitare nello scope a livello di processo qui. CPUs equals one, memory equals two gigabytes. La nostra pipeline non è molto esigente, quindi dovrebbe andare bene.

Ora, ho appena indovinato questi numeri qui, ma come fate a sapere qual è una quantità sensata di risorse da usare? È un lavoro piuttosto difficile andare a scavare attraverso tutti questi diversi processi di una grande pipeline di molti campioni e capire quale fosse l'utilizzo delle risorse.

Quindi un buon approccio per questo è impostare questi valori a numeri alti per iniziare, solo così che la vostra pipeline venga eseguita senza errori, e poi chiedere a Nextflow di generare un report di utilizzo per voi.

Questo è super facile da fare, quindi tornerò a un terminale. Oh, devo ricordarmi di reimpostare quello su local così che la mia pipeline venga effettivamente eseguita. E dirò Nextflow run, e userò un flag della riga di comando -with-report.

E posso lasciarlo vuoto e darà un nome file predefinito, ma gli darò un nome file specifico così che venga salvato in un posto specifico.

Premo Invio, e la pipeline viene eseguita esattamente come normale, ma quando finisce, genererà un bel report HTML per me.

Quindi nella barra laterale qui, ho questo file HTML. Se lo stessi eseguendo localmente, lo aprirei semplicemente. Poiché sono in Codespaces, farò clic destro su quello e clic su download, che lo scaricherà sul mio computer locale. E posso semplicemente aprirlo facilmente nel browser web.

Nextflow può generare un report come questo per qualsiasi pipeline ed ha alcune informazioni davvero belle. Quindi è buona pratica salvare sempre queste cose. Ci dice quando abbiamo eseguito, dove abbiamo eseguito, se ha avuto successo o no, quali parametri sono stati usati, quale era il comando CLI, cose del genere.

E ci sono anche questi grafici sull'utilizzo delle risorse. Quindi ci dice quale percentuale di chiamate CPU sono state usate per ogni processo come box plot qui, perché ci sono molte attività per ogni processo, quindi possiamo vedere la distribuzione.

Potete vedere i nostri processi qui, cowpy e collectGreetings avevano solo una singola attività, quindi è solo una singola linea. E abbiamo sia CPU che memoria e durata del job, ed erano molto veloci.

Se state usando Seqera Platform, tra l'altro, ottenete gli stessi grafici integrati nell'interfaccia Platform senza dover fare nulla. Quindi avete sempre queste informazioni a portata di mano.

Okay, quindi possiamo usare questo report e su un'esecuzione reale, e farci un'idea di quante CPU e quanta memoria viene usata dalla nostra pipeline e tornare e mettere quei valori di nuovo nel nostro file di configurazione, così che la prossima volta forse non richiediamo così tanto. E possiamo essere un po' più parsimoniosi.

Ora potete diventare davvero intelligenti nel configurare i file di configurazione della pipeline. E di nuovo, se state usando Seqera Platform, cercate un piccolo pulsante che sembra una lampadina. Perché se cliccate quello, genererà un file di configurazione altamente ottimizzato, che è adattato specificamente per i vostri dati, la vostra esecuzione e la vostra pipeline. Per eseguirla nel modo più efficiente possibile.

Ma per ora, dirò che in realtà il numero predefinito di CPU che Nextflow stava dando andava bene e ma serve solo un gigabyte di memoria.

## 5.3. Impostare le allocazioni di risorse per un processo specifico

Ora, nella vita reale, è piuttosto insolito che tutti i processi nella vostra pipeline abbiano bisogno degli stessi requisiti. Potreste avere qualcosa come MultiQC come strumento di reporting, che ha bisogno di molto poco in termini di risorse e viene eseguito abbastanza velocemente.

E poi forse avete qualcosa che sta indicizzando un genoma di riferimento o facendo qualche allineamento o facendo qualche altro lavoro. Non importa cosa sia, che richiede molte risorse. E quindi per queste diverse sottomissioni di job a uno scheduler, volete dare diverse quantità di risorse.

Sotto questo scope process, possiamo definire una configurazione, che punta a processi specifici in modi diversi.

Qui stiamo usando withName, possiamo anche usare label, e questi possono usare un pattern per puntare a uno o più processi. Qui stiamo solo dicendo qualsiasi processo che abbia un nome cowpy impostato a due gigabyte di memoria e due CPU, e poiché questo è un selettore più specifico rispetto a quello di processo di livello superiore, questo viene sovrascritto in questi casi, quindi potete costruire un bel file di configurazione qui, che adatta davvero tutti i vostri diversi processi nella vostra pipeline per renderli davvero efficienti.

## 5.5. Aggiungere limiti di risorse

Ora come sviluppatore della pipeline, probabilmente conosco abbastanza bene gli strumenti, e voglio che tutto venga eseguito il più velocemente e il meglio possibile. Quindi potrebbe essere che metta numeri piuttosto alti per alcuni di questi perché so che verrà eseguito molto più velocemente se do a cowpy 20 CPU invece di due.

Va bene finché non andate a eseguire sul vostro laptop o su GitHub Actions Continuous Integration test, o qualche altro sistema, che forse non ha 20 CPU disponibili.

Ora quando provate a eseguire la pipeline, si bloccherà perché Nextflow dirà, non posso inviare questo job da nessuna parte. Non ho le risorse disponibili.

Ora per evitare quel crash duro, possiamo aggiungere un po' più di configurazione, che è specifica per il nostro sistema ora, chiamata limiti di risorse. E appare così. È sotto lo scope process di nuovo.

E limiti di risorse, potete specificare fondamentalmente il tetto di ciò che avete disponibile. È una mappa qui, e potete, all'interno di questa mappa, potete impostare la memoria, le CPU e il tempo.

Ora quello che succede è quando Nextflow invia un'attività da un processo, guarda cosa è richiesto e fondamentalmente fa solo un minimo tra quello e quello. Quindi se abbiamo richiesto 20 CPU, ma solo quattro sono disponibili, richiederà quattro. La pipeline non si blocca e usa il più vicino possibile a ciò che è stato progettato dallo sviluppatore della pipeline.

## 6. Usare i profili per passare tra configurazioni preimpostate

Okay. Ho detto che i limiti di risorse qui potrebbero essere specifici del sistema, e forse ho un file Nextflow config nella mia pipeline, e so che le persone useranno questo in una gamma di posti diversi. Ora, invece di forzare tutti a creare il proprio file Nextflow config ogni singola volta, quello che posso fare è raggruppare diversi preset di configurazione insieme in profili di configurazione.

Scorrerò un po' qui e poi appena dopo params, perché l'ordine del file di configurazione qui è importante, il file di configurazione viene caricato sequenzialmente, quindi metterò questi profili dopo tutto il resto così che sovrascriva i parametri definiti in precedenza. E incollerò questi profili dal materiale di formazione.

Quindi c'è un nuovo scope di livello superiore chiamato profiles. Possiamo avere nomi arbitrari qui. Quindi abbiamo my_laptop e univ_hpc. E qui possiamo vedere che stiamo impostando gli altri stessi parametri di configurazione che avevamo prima. Ora solo all'interno di un profilo. Quindi abbiamo un executor locale per eseguire sul mio laptop e sto inviando a un cluster SLURM sull'HPC.

Sto usando Docker localmente, conda sull'HPC, e il sistema HPC ha limiti di risorse molto più alti.

Ora posso eseguire la pipeline con l'opzione CLI -profile, dire quale profilo voglio usare. Quindi userò my_laptop, e Nextflow applicherà tutta la configurazione all'interno di quello scope del profilo. Quindi posso provarlo ora. È lo stesso comando di prima. Nextflow run hello config, e faccio dash profile, singolo trattino perché è l'opzione core di Nextflow, dash profile my_laptop.

Ora applicherà in batch tutta quell'opzione di configurazione. Oh, e potete vedere, ho detto prima che questo potrebbe succedere che il requisito del processo, ha chiesto quattro CPU e ne ho solo due su questa istanza Codespaces.

Quindi questa è una buona opportunità solo per provare i limiti di risorse del processo, e dire che ho solo due CPU sul mio laptop, o in questi Codespaces. Ora se lo eseguiamo di nuovo, dovrebbe limitare quel requisito a due e si spera che la pipeline venga eseguita. Ottimo.

## 6.2. Creare un profilo di parametri di test

Notate che questi profili non devono avere solo configurazione sulla loro infrastruttura. Potete avere raggruppamenti di qualsiasi configurazione qui, inclusi parametri.

Quindi un'altra cosa che vedrete molto spesso nelle pipeline delle persone è un profilo test, che include parametri, che normalmente inviereste su base per utente. Ma qui abbiamo, fondamentalmente diversi valori predefiniti sensati per quando voglio eseguire casi di test.

E questo è fantastico perché non devo necessariamente andare a specificare tutte queste cose, che potrebbero essere parametri richiesti. Altrimenti posso semplicemente dire dash profile test e verrà eseguito subito.

Ora qualcosa da notare è che i profili possono anche essere combinati più di uno. Quindi posso fare profile my_laptop qui, e poi aggiungere anche test. Non faccio profile due volte. Faccio solo un elenco separato da virgole qui senza spazi. E applicherà questi profili in ordine. Quindi prenderà la configurazione dal profilo my_laptop, e poi applicherà la configurazione test sopra.

Davvero conveniente e potete vedere come potete impostare molti gruppi predefiniti sensati qui per rendere facile eseguire la vostra pipeline.

## 6.3. Usare nextflow config per vedere la configurazione risolta

Si spera, vi ho convinto che la risoluzione della configurazione di Nextflow è potente, ma non vi biasimerei se state diventando un po' strabici a questo punto dopo che ho detto circa 20 modi diversi per fornire configurazione e dare tutti questi diversi strati come una buccia di cipolla.

Quindi se mai vi sentite insicuri su quale sia la configurazione finale risolta per Nextflow, sappiate che c'è un comando chiamato "nextflow config", e possiamo eseguirlo e ci dirà quale è la configurazione risolta nella nostra posizione attuale.

Quindi quando lo eseguo qui, trova il file "nextflow.config" nella directory di lavoro corrente, e processa tutta la diversa configurazione, e mi dà l'output risolto.

Notate che il file di configurazione Nextflow può anche prendere l'opzione CLI profile. Quindi se gli dico di risolvere nei profili my_laptop e test, e potete vedere che ha anche applicato i limiti di risorse qui dall'opzione di configurazione my_laptop e anche impostato i parametri, che erano nel test.

Quindi questo è un modo carino solo per esplorare come funziona la risoluzione della configurazione, se siete per niente insicuri.

## Conclusione

Okay, questo è tutto. Questa è la configurazione Nextflow in poche parole. Potete fare molte cose con la configurazione. È davvero potente. Ma questi sono la maggior parte dei casi d'uso comuni che vi troverete a fare, e questi concetti si applicano a tutte le diverse opzioni.

Datevi una pacca sulla spalla perché questa è la fine del corso di formazione Hello Nextflow. Si spera che ora siate sicuri sia nello scrivere la vostra pipeline Nextflow da zero, configurarla ed eseguirla, e conoscete tutti i dettagli e le cose a cui prestare attenzione.

C'è un altro quiz che potete provare sulla pagina di formazione sulla configurazione. Quindi andate giù e provatelo e assicuratevi di aver capito tutte queste parti sulla configurazione.

E unitevi a noi nell'ultimo video solo per una rapida conclusione su alcuni dei prossimi passi che potrebbero essere buoni da fare dopo questo corso di formazione.

Grazie per essere rimasti con noi. Ben fatto e ci vediamo nel prossimo video.
