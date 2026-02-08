# Parte 5: Hello Containers - Trascrizione Video

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduzione assistita da IA - [scopri di più e suggerisci miglioramenti](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/Xqr--bKEN9U?si=QinuAnFwFj-Z8CrO&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Note importanti"

    Questa pagina mostra solo la trascrizione. Per le istruzioni complete passo dopo passo, tornate al [materiale del corso](../05_hello_containers.md).

    I numeri delle sezioni mostrati nella trascrizione sono forniti solo a scopo indicativo e potrebbero non includere tutti i numeri di sezione presenti nei materiali.

## Benvenuto e contesto

Ciao e bentornati a Hello Nextflow. Questa è la parte cinque chiamata Hello Containers. E in questa parte del corso parleremo di come incapsulare i requisiti software per una pipeline, in modo che le persone che eseguono la pipeline non debbano preoccuparsi di installare il software.

Se lavorate in bioinformatica da tanto tempo quanto me, potreste ricordare quelli che spesso chiamo i brutti vecchi tempi, quando per eseguire la pipeline di qualcun altro o replicare il loro lavoro, passavate ore o giorni cercando di installare tutti i diversi strumenti software che avevano usato, alle stesse versioni, cercando di compilarli sulla vostra macchina, ed era un incubo. Era davvero difficile.

Se lavoravate su un HPC, potreste aver usato i moduli di ambiente dove i sysadmin cercavano di installare il software per voi, il che andava bene, ma era comunque imperfetto.

Ma ora abbiamo modi migliori per fare questo. Nextflow ha supporto integrato per diverse tecnologie di container software. Docker è quella più comune. È quella che useremo oggi. Funziona bene in Codespaces. Funziona bene sul vostro computer locale e funziona bene nel cloud.

Ma anche Singularity o Apptainer, che sono molto comuni sui sistemi HPC e funzionano in modo praticamente identico. O Podman, Shifter, ci sono un sacco di altri che sono tutti molto simili.

L'unica cosa extra, che è un po' simile ma non proprio uguale, che Nextflow supporta è Conda. E Nextflow può gestire ambienti Conda per voi su base per processo, il che è molto meglio che gestire i vostri ambienti Conda. E ancora, può essere distribuito con una pipeline.

Inizieremo questo capitolo parlando un po' delle tecnologie container e Docker e di come funzionano. E faremo la prima metà manualmente in Docker in modo che capiate cosa succede sotto il cofano e come funziona. Perché è davvero importante capire cosa sta facendo Nextflow e come capire cosa sta facendo il vostro flusso di lavoro quando viene eseguito.

Quindi. Saltiamo nei nostri Codespaces. Ora ho ripulito tutto di nuovo, ma se andiamo in Hello Containers, dovreste vedere che tutti i nostri script e tutto il resto sono lì, gli stessi della fine del capitolo sui moduli. Quindi abbiamo i nostri diversi moduli qui, che ho creato nella directory modules.

Sono ancora lì. Devono esserci perché possa funzionare. E il flusso di lavoro e l'output sono tutti uguali, tranne che abbiamo cambiato il percorso di pubblicazione dell'output in Hello Containers, in modo che i vostri file finiscano in quella directory.

Possiamo eseguire questo ora per verificare che funzioni, se volete, oppure possiamo continuare con il terminale.

## 1. Usare un container 'manualmente'

Useremo Docker per gestire i nostri container, e posso verificare che sia installato sui miei Codespaces facendo "docker -v", che mi mostra la versione installata e tutto il resto, e che funzioni correttamente.

Ora i container e Docker hanno due concetti che sono davvero importanti. Uno si chiama image, e uno si chiama container. L'image è lo snapshot, se volete, dell'intero file system che userete, e il container è l'ambiente in esecuzione. Quindi create un container usando un'image.

Una volta che siete in quel container, funziona tipicamente come un intero sistema operativo. È isolato dal mondo esterno. È separato da tutto il resto, e questa è una cosa buona. È così che otteniamo una riproducibilità così buona con Nextflow.

Perché per le attività eseguite all'interno di un container, non sono contaminate da alcun file di configurazione sul vostro sistema locale. Nessun'altra influenza esterna, vengono eseguite nel loro piccolo ambiente isolato. I file vengono quindi prodotti in modo molto, molto riproducibile perché state usando le stesse librerie sottostanti, tutte le stesse dipendenze, esattamente lo stesso software per ogni persona che esegue su ogni diverso ambiente di calcolo. Il che francamente penso sia fantastico e incredibile che funzioni. E ancora, ancora ad oggi mi stupisce che questo sia possibile.

## 1.1. Scaricare l'image del container

Quindi proveremo ad usare alcune image Docker e Docker, quando lo eseguite sul vostro sistema, ha un registro docker sul vostro computer, o in questo caso, nel code space, che tiene traccia di tutte le diverse image che sono state scaricate e usate in passato, e dei diversi livelli di cui sono composte.

Possiamo vedere quali image abbiamo localmente con Docker facendo "docker image ls". E in questo caso potete vedere che ci sono un sacco di image Docker qui, che hanno tutte a che fare con la configurazione di questi Codespaces. Hanno tutte a che fare con i dev container e cose del genere. Quindi non dovete preoccuparvene troppo, ma man mano che aggiungiamo più image e le scarichiamo, man mano che questo corso va avanti, potete controllare quella lista e vedrete che il registro locale tiene traccia di tutte queste cose che abbiamo scaricato.

Ma ne prenderemo una nuova facendo "docker pull". E questo dice a Docker di recuperare una nuova image dal web.

Mettiamo poi l'URI per quel container. Ora questa potrebbe essere un'image docker che avete costruito localmente e poi caricato su internet. Potrebbe essere un'image che qualcun altro ha fatto. Ci sono molti, molti, molti modi diversi per creare image Docker, ma probabilmente uno dei modi più semplici è esternalizzare questo, e far sì che qualcun altro lo faccia per voi.

E quello che useremo in questo tutorial è un servizio di Seqera chiamato Seqera Containers.

Ora, Seqera Containers è totalmente gratuito, e usa un pezzo di software open source che sviluppiamo chiamato Wave, che è stato costruito per gestire i container in modo complementare a Nextflow. E gestisce molti dei casi d'uso comuni che ci troviamo a gestire con Nextflow.

È molto comune che il software di cui abbiamo bisogno sia pacchettizzato in Conda, nei canali Bioconda o conda-forge o altri canali più specifici per dominio. E Wave e Seqera Containers sono davvero bravi a costruire image da questo.

Quindi posso andare su questa interfaccia web e giocheremo con il pacchetto chiamato "cowpy". Quindi digito il nome del pacchetto che voglio. Cerca, l'ha trovato sul Python package index, quindi posso usarlo. Oppure se aspetto un po' più a lungo, sta cercando bioconda e conda-forge. E potete vedere, posso specificare qualsiasi canale con qui. Quindi se volete trovare un canale Nvidia o qualsiasi altra cosa, dovrebbe funzionare anche quello.

E poi posso specificare se voglio che costruisca un'image docker per me o un'image singularity e anche quale architettura CPU voglio usare. Quindi amd64 o arm64.

E una volta che i risultati di bioconda sono elencati, posso ora vedere tutte le diverse versioni disponibili. Lo inserirò. E ora potrei continuare a cercare e ottenere più pacchetti da Conda se voglio e comporre questo container come voglio, ma ne voglio solo quello. Quindi cliccherò su Get Container.

Ora, qualcun altro ha già richiesto lo stesso container prima ed è restituito da un registro, quindi lo otteniamo immediatamente. Ma se nessun altro avesse mai chiesto questo pacchetto software o questa combinazione di pacchetti software, Wave e Seqera Containers lo costruirebbero al volo per noi.

Possiamo copiare questo URL e possiamo anche vedere i dettagli della build. E questo ci mostra cosa ha fatto il servizio nel backend. Ha creato un file di ambiente conda. Un dockerfile, e poi questo è il processo di build docker in esecuzione. Ha anche eseguito una scansione, una scansione di sicurezza, quindi potete vedere eventuali CVE. E vi dice quando è stato creato.

Wave e Seqera Containers possono fare molto di più di questo, ma questo è un caso d'uso semplice, che è il più comune. E dovrei dire che queste image sono ospitate per almeno cinque anni. Quindi potete costruire questi URL nelle vostre pipeline e sapere che non spariranno presto.

Quindi ho il mio URL per la mia image docker per cowpy.

Posso ora fare "docker pull" di quell'URL, e recupererà tutti i diversi livelli e scaricherà questa image in modo che sia disponibile per me localmente.

## 1.2. Usare il container per eseguire cowpy come comando singolo

Okay, ora proviamo ad usarlo davvero. Quindi ora userò un comando "docker run" invece di "docker pull", e userò il flag "--rm", che dice semplicemente a Docker di chiudere questo container una volta che ha finito di fare quello che gli ho chiesto. E poi metto l'identificatore per il container, che è solo un URI.

E poi alla fine, specifico il comando che voglio che Docker esegua all'interno del container generato da questa image. Dirò solo cowpy, che è il nome dello strumento installato da Conda Forge, che è disponibile all'interno dell'image.

Premo invio e voilà. Abbiamo eseguito cowpy su un sistema. Abbiamo una piccola mucca che ci dà alcune informazioni.

Notate che cowpy non è installato sul mio sistema locale. Quindi se lo eseguo senza tutta la roba di Docker, dice, comando non trovato. Quindi questo ha scaricato un'image. Ha creato un container usando Docker, e poi è entrato in quel container e ha eseguito questo comando per noi e ci ha restituito l'output nel nostro terminale. Molto, molto bello.

## 1.3. Usare il container per eseguire cowpy interattivamente

Okay, faremo un passo avanti ora ed eseguiremo questo container interattivamente e daremo un'occhiata in giro, così possiamo vedere cosa sta succedendo all'interno del container.

Quindi se torno indietro e prendo il mio comando run e toglierò cowpy alla fine, perché in realtà non voglio eseguire cowpy. Voglio eseguire un terminale Bash.

E poi tornerò qui e farò "-it", che sta per Interactive e Terminal o TTY, e premerò invio.

E ora potete vedere il prompt, la parte prima che digito, è cambiato. Questo era il prompt di Codespaces dove diceva la directory, e ora dice base e root e tmp. Quindi ora sono all'interno del container, e se faccio "ls", vedrete che i file che vedo in questa directory sono diversi dai file che ho nel mio workspace.

E infatti, non posso vedere nessuno dei file dal mio workspace locale di codespaces o dal mio disco locale all'interno del container Docker. Il runtime del container docker è completamente isolato e non può scrivere o leggere alcun file da un file system host esterno.

Posso, tuttavia, vedere il software installato all'interno del container ed eseguirlo. Quindi posso eseguire cowpy e possiamo vedere un po' di più su come usare cowpy. Qui posso fare "cowpy 'Hello World'" e questo gli dice di mettere effettivamente la mia citazione dentro una piccola bolla di discorso. E potete anche eseguire diversi tipi di mucche, quindi non deve essere una mucca. Potete fare un "-c". E sono in Svezia, quindi sceglierò un alce. Molto carino. Gli ho dato delle corna.

E ce ne sono un sacco di diversi con cui potete giocare, che potete vedere descritti nei documenti di formazione.

## 1.3.4. Montare dati nel container

Okay. Sarebbe bello se potessimo eseguire cowpy sui file nel nostro file system.

Naturalmente, non è molto utile avere solo il container e nessun accesso a nulla. Potrebbe essere sicuro e riproducibile, ma non è molto utile.

Quindi come facciamo? Uscirò da questo container Docker digitando exit, e potete vedere che il prompt ci dice che ora siamo di nuovo nei nostri Codespaces normali.

Ed eseguirò di nuovo lo stesso comando. Ma questa volta aggiungerò alcuni flag aggiuntivi qui. E quello importante è "-v", che sta per montare un volume, che è fondamentalmente una parte di uno spazio disco.

Il "-v" prende due parti: c'è una stringa e poi due punti e una stringa. E la prima parte è il file system locale, che dovrebbe essere montato nel container. E poi la seconda parte è dove dovrebbe finire all'interno del container.

Ora voglio solo caricare tutto il mio file system locale qui. Quindi "." è la directory di lavoro corrente. Quindi farò solo "." e poi ":", e poi metteremo questo in una nuova directory all'interno del container chiamata "my_project". Potrebbe davvero chiamarsi in qualsiasi modo.

E poi eseguirò di nuovo.

Nella directory di lavoro dove sono scaricato, che è /tmp, i file non ci sono. Ma se faccio "ls my_project", eccolo: tutti gli stessi file che avevamo localmente sui nostri Codespaces sono ora disponibili all'interno del container in quel percorso.

Questo è accesso in lettura e scrittura, quindi posso creare nuovi file in questa directory e appariranno sul mio file system host. Questa directory particolare, quindi si comporta esattamente come se fossi fuori dal container, quindi posso ora leggere e scrivere e fare cose.

## 1.3.5. Usare i dati montati

Okay, dimostriamo solo che possiamo farlo. Faccio "cat /my_project/data/greetings.csv". Se ricordate, i contenuti di questo file sembrano così. Posso ora inviare questo a cowpy e la mucca stamperà i diversi output di quel file nella sua piccola bolla di discorso, il che è un po' divertente.

Quindi potete vedere, possiamo ora usare il software nel container per interagire con i file sul nostro sistema host.

Okay, usciamo di nuovo e proseguiremo con il resto del materiale di formazione.

## 2. Usare i container in Nextflow

Quindi è davvero bello usare i container. Spero che abbia senso. E potete vedere il valore di questi container e perché è utile per eseguire software di analisi.

Ma come facciamo tutto questo stesso processo all'interno di Nextflow? Non vogliamo eseguire un sacco di comandi Docker noi stessi. Vogliamo solo lasciare che Nextflow gestisca tutto questo per noi.

Quindi lavoriamo su questo. Aggiungeremo un nuovo processo alla nostra pipeline, per eseguire cowpy. Okay, quindi creiamo un nuovo modulo per il nostro nuovo processo. Quindi andiamo in modules, chiamiamolo cowPy.nf, e poi copierò il codice dal materiale di formazione qui.

Ma potete vedere che il processo è molto semplice. Assomiglia molto a quelli che abbiamo fatto finora, abbiamo un blocco input con un path, che è il nostro file di input, e anche un value qui, quindi questo sarà un carattere, quindi potremmo usare di nuovo un alce se vogliamo.

E poi un output, che è un singolo file qui, un path e poi uno script. E stiamo facendo la stessa cosa che abbiamo fatto interattivamente all'interno del container: stiamo facendo "cat" per leggere il file di input. Stiamo inviando quel contenuto a cowpy. Stiamo scegliendo un carattere specifico basato su quell'input, stiamo scrivendo in un file di output chiamato cowpy input file, che viene poi inviato all'output.

Ottimo. Includiamo quello. Quindi include \{ cowpy \} from "./modules/cowpy.nf", l'ho chiamato cowpy? Sì.

E poi chiamiamo il nostro nuovo processo qui nel blocco main del flusso di lavoro. Quindi eseguiamo cowpy. E prenderemo il nostro nuovo processo cowpy e diremo collectGreetings.out.

E poi se ricordate, c'erano due output per questo modulo. Uno chiamato outfile e uno chiamato report. L'estensione VS Code ci sta suggerendo automaticamente questi e vogliamo .outfile.

Potete sempre saltare in questo processo qui. O passate sopra e dovrebbe mostrarvi rapidamente quali erano gli output. E possiamo anche fare command click su di esso e aprirà il file del modulo se volete vedere più in dettaglio.

Quindi eccoci qui. Quello è l'outfile lì, e quello è il path. Quindi ora sarà il file di input per il nostro processo cowpy. Fantastico.

Ora se ricordate, un processo cowpy ha due input. Avevamo anche il canale value per il carattere. Quindi possiamo aggiungere "params.character" qui. Avrei potuto codificarlo direttamente se avessi voluto, ma rendiamolo un'opzione CLI così possiamo fare dash, dash character.

Giusto. Ora devo definire il parametro di input che abbiamo appena chiamato e dargli un default. Quindi character, String. E mi piace l'alce, quindi lo imposterò su moose per default.

Giusto, proviamo ad eseguirlo. Quindi se faccio Nextflow run hello containers, vedremo cosa succede.

Avrei potuto usare dash resume se avessi le vecchie directory di lavoro in giro. E ancora, questi primi processi sarebbero stati memorizzati nella cache e sarebbe stato un po' più veloce, ma dovrebbe essere sostanzialmente lo stesso.

Ora possiamo vedere subito che ha generato un errore quando è arrivato al nostro nuovo processo, ci sta dicendo qui che c'è stato un errore nell'esecuzione del processo cowpy ed è uscito con uno status di uscita 127. Questo è il comando che ha cercato di eseguire. Sembra giusto, sembra come ci aspettavamo. Sta prendendo quel nome di file di output, che sembra giusto, lo sta eseguendo con un carattere alce e sta cercando di salvarlo.

Ma potete vedere l'errore del comando qui che dice che il comando cowpy non è stato trovato. E ha senso perché non abbiamo effettivamente detto a Nextflow di usare ancora un container. Abbiamo solo dato il comando cowpy. E come ho detto prima, cowpy non è installato sul nostro sistema locale. Quindi quando ha cercato di eseguirlo, è fallito.

## 2.3.1. Specificare un container per cowpy

Dobbiamo dire a Nextflow che c'è un container disponibile e può usarlo. Quindi come facciamo?

Se andiamo nel nostro modulo qui, aggiungeremo una nuova dichiarazione in alto chiamata "container". E la imposteremo su una stringa.

Ora, se ricordate, in Seqera Containers, posso copiare quell'URL e lo inserisco tra virgolette qui.

Ora torniamo indietro e proviamo ad eseguirlo di nuovo.

Vediamo se funziona questa volta.

Sfortunatamente, fallisce esattamente nello stesso modo, anche se ora abbiamo definito un container per il processo da eseguire. Quindi per usare la nostra image docker, dobbiamo dire a Nextflow di abilitare l'uso di Docker quando eseguiamo il flusso di lavoro.

E lo faremo creando un nuovo file di configurazione. Quindi dirò touch nextflow.config.

Questo è un nome di file speciale dove se è nella directory di lavoro mentre lancio la pipeline, verrà caricato automaticamente. Quindi se vado in questo file Nextflow dot config, potete vedere che in realtà esiste già, che avevo dimenticato. E abbiamo docker.enabled qui già, ma è impostato su false, che è il default.

Quindi se cambio quello in equals True invece, docker.enabled. E ci sono documenti di riferimento per tutti questi scope di configurazione nei documenti Nextflow. E inoltre potete vedere che quando passo sopra con un'estensione VS Code, carica i documenti specifici per questo e mi dice cosa significa e come impostarlo.

Quindi ora l'abbiamo impostato su true, e se eseguo Nextflow di nuovo, Nextflow saprà ora di recuperare quell'image docker per noi se non l'abbiamo già localmente, e poi eseguire quel processo con quell'ambiente container.

E quindi possiamo vedere che è stato eseguito con successo e abbiamo un piccolo segno di spunta accanto a cowpy. Fantastico. Se vado su e guardo nella directory dei risultati, il file non c'è ancora. E questo perché dobbiamo ancora pubblicare questo file di output proprio come tutti gli altri.

Quindi andiamo al blocco published all'interno del flusso di lavoro, diciamo mycowpy equals cowpy.out.

E poi qui nel blocco output, mycowpy, parentesi graffe path. Oops. Hello containers. Mode, copy.

Se eseguo di nuovo ora, dovrebbe essere eseguito esattamente nello stesso modo. Avrei potuto usare dash resume e dimentico ogni volta. E poi vado su e ora abbiamo un nuovo file creato chiamato cowpy-COLLECTED, e c'è il mio alce che dice BONJOUR, HELLO, HOLÀ. Fantastico.

Ora naturalmente potrei anche passare ora "--character". Quali sono le diverse opzioni? Penso che ci sia un Turkey? Quindi posso usare character Turkey. Verrà eseguito esattamente nello stesso modo. Ho perso un'altra opportunità di usare dash resume, e ora se carichiamo il nostro file e ora abbiamo un Turkey. Fantastico.

## 2.3.4. Ispezionare come Nextflow ha lanciato l'attività containerizzata

Okay. Ultima piccola cosa. Eseguiamo rapidamente questo comando di nuovo, resume questa volta, e diamo una rapida occhiata nella directory di lavoro per vedere cosa sta facendo Nextflow sotto il cofano per far funzionare tutto questo per noi.

Questa volta è super veloce, andiamo in questa directory di lavoro, cd work/. Ora se ricordate abbiamo un sacco di dot file qui e quello che ci interessa in questo caso è quello che ho detto che quasi mai dobbiamo guardare, chiamato .command.run.

Se faccio code dot command run, lo aprirà nell'editor. E posso cercare in questo file e se scorro verso il basso dovrei vedere Docker run. E potete vedere che Nextflow sta facendo il comando docker run per noi, quando Docker è abilitato in una configurazione. Ha un sacco di diversi flag e cose qui, ma potete vedere il flag "-v" che abbiamo usato noi stessi quando stavamo eseguendo. E potete vedere che sta montando la directory workspace locale nel container, in modo che il container possa accedere ai nostri file di input e salvare gli output. E poi alla fine, sta anche eseguendo .command.sh, che è lo script generato, che ha il comando cowpy dentro.

E quindi potete vedere che Nextflow sta prendendo la logica del flusso di lavoro, che è la roba che effettivamente ci interessa, che è specifica per la nostra analisi, e sta facendo tutte le cose intelligenti dietro le quinte per far funzionare Docker sul nostro sistema.

E lo sta facendo in modo davvero portabile in modo che un utente finale della pipeline possa cambiare la tecnologia che sta usando: Docker, Singularity, Apptainer, Conda. Questo non importa davvero alla logica della pipeline, ma Nextflow gestirà tutte le esigenze infrastrutturali sottostanti, in modo che funzioni ovunque.

E questo è davvero il superpotere di Nextflow. È riproducibilità e portabilità. E con Nextflow potete effettivamente condividere il vostro flusso di lavoro e altre persone possono eseguirlo sui loro sistemi e funzionerà semplicemente.

È una cosa davvero, davvero difficile da fare, e ora sapete come farlo anche con i vostri flussi di lavoro.

Okay, questo è tutto per questo capitolo. Se andate in fondo al corso, troverete di nuovo un quiz sui container. Spero che abbia tutto senso. È un modo davvero bello di lavorare con l'analisi. E se siete nuovi ai container, spero di avervi convinto che è la strada da percorrere, e non tornerete mai indietro.

Ma con questo, fate una piccola pausa forse, e vi unirò tra qualche minuto per passare attraverso la parte finale sei di Hello Nextflow, che è tutta sulla configurazione.

Grazie mille.
