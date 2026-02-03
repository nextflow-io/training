# Parte 5: Hello Containers - Trascrizione

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduzione assistita da IA - [scopri di più e suggerisci miglioramenti](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/5PyOWjKnNmg?si=QinuAnFwFj-Z8CrO&amp;list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Note importanti"

    Questa pagina mostra solo la trascrizione. Per le istruzioni complete passo dopo passo, torni al [materiale del corso](../05_hello_containers.md).

    I numeri di sezione mostrati nella trascrizione sono forniti solo a scopo indicativo e potrebbero non includere tutti i numeri di sezione presenti nei materiali.

## Benvenuto

Salve, benvenuti alla Parte Cinque del corso di formazione Hello Nextflow.

Questo capitolo si intitola Hello Containers. Parleremo di come Nextflow si integra con strumenti come Docker e Singularity per utilizzare container software per fornire software agli utenti del Suo pipeline.

Questo significa che quando le persone eseguono il Suo pipeline, non devono andare ad installare tutti i diversi strumenti da soli. Nextflow lo farà per loro.

I container sono una tecnologia estremamente potente e cruciale per la riproducibilità e la facilità d'uso. Inizieremo con una breve introduzione ai container stessi, eseguendo alcuni comandi docker manualmente, e poi prenderemo gli stessi container e li inseriremo nel nostro pipeline Nextflow.

Bene. Iniziamo.

Quindi, proprio come prima, iniziamo caricando il materiale di formazione. Vada su training.nextflow.io. Hello Nextflow, Capitolo Cinque, Hello Containers.

Salto nel mio ambiente Codespaces e sulla sinistra qui vediamo hello containers punto nf.

Proprio come prima, questo è lo stesso script con cui abbiamo finito il capitolo quattro precedente, quindi dovrebbe essere familiare.

Abbiamo i nostri parametri da linea di comando per specificare il file di input e il nome del batch. Stiamo includendo i nostri tre moduli, e abbiamo il nostro workflow dove eseguiamo i tre processi.

## 0. Riscaldamento: Eseguire hello-containers.nf

Si senta libero di eseguire nuovamente questo workflow e verificare che stia producendo gli output che si aspetta. Per ora, in realtà lo chiuderò e mi immergerò nel terminale.

## 1. Usare un container 'manualmente'

Per iniziare questo capitolo, faremo un po' di riepilogo sulla tecnologia dei container. Se è molto abituato a docker o singularity o altre tecnologie di container, allora lo consideri un ripasso, o si senta libero di saltarlo completamente.

Nextflow supporta molti tipi diversi di tecnologie di container. Ciò include Docker, Singularity, Podman, Shifter, Charliecloud e altri.

In questa formazione, ci concentreremo su Docker. Questo viene preinstallato negli spazi di codice ed è una delle tecnologie di container più popolari, specialmente se sta sviluppando sul Suo computer o sul Suo laptop.

Se sta lavorando in un ambiente accademico su un HPC condiviso, potrebbe scoprire che Singularity è disponibile e Docker no. Va bene. Tutti i concetti sono esattamente gli stessi. Alcuni dei comandi manuali sono diversi, ma se comprende Docker, comprenderà anche singularity.

In effetti, Singularity è installato anche nell'ambiente Code Spaces. Quindi, se vuole, può provare a eseguire le stesse attività utilizzando Singularity invece di Docker.

Bene, cos'è la tecnologia dei container? L'idea dietro Docker è che può recuperare un'immagine da una fonte remota. Scaricarla sulla Sua macchina locale e quindi creare un container basato su quell'immagine.

Questo container in esecuzione è un po' come una macchina virtuale in esecuzione sul Suo computer. È isolato dal Suo ambiente e viene preconfezionato con un sistema operativo e un insieme di software disponibili.

## 1.1. Scaricare l'immagine del container

La sintassi di cui abbiamo bisogno per recuperare un'immagine preesistente è "docker pull". Quindi lo digiterò nel mio terminale, ma ora abbiamo bisogno di un'immagine con cui giocare.

Può creare immagini Lei stesso. Può trovarle su registri pubblici come Docker Hub o quay.io. Ma un modo davvero efficace per ottenere immagini rapidamente è utilizzare Seqera Containers.

Questo è un servizio comunitario gratuito che abbiamo costruito nel 2024, che può utilizzare senza login o altro.

Se va su seqera.io/containers o clicca su containers in alto qui, Le viene presentata un'interfaccia di ricerca e può digitare il nome di qualsiasi strumento disponibile in Conda o sul Python Package Index.

Per impostazione predefinita, cerca i canali Bioconda e Conda Forge, ma può prefissare qualsiasi canale Conda. Sono qui se vuole.

Per un po' di divertimento, usiamo cowpy. Digiterò cowpy. Mi dà risultati da Python Package Index e Conda Forge. Cliccherò su quello per aggiungerlo al mio container. Potrei aggiungere più pacchetti qui se volessi. Seleziono Docker, seleziono linux/amd64 e clicco Get Container.

Questo costruisce l'immagine per me su richiesta se non è già stata creata, e mi dà un URL che posso copiare.

Se è interessato, può cliccare view Build Details, e questo La porta a una pagina, che mostra il file di ambiente conda che è stato utilizzato e il log di build completo per la build, insieme ai risultati della scansione di sicurezza.

Se torno ai miei code spaces, ora posso incollare questo nome di container e premere invio.

Docker ora scarica tutti i diversi layer all'interno di questa immagine del container, e ora ci dice che questa immagine è disponibile per l'uso.

## Scaricare un'immagine Singularity

Se sta utilizzando singularity, il processo è sostanzialmente lo stesso. Selezioniamo i nostri pacchetti di immagini, selezioniamo cowpy. Ora scegliamo Singularity invece di Docker e clicchiamo Get Container. Questo ci dà un URL dell'immagine usando oras://. Oppure, se preferisce, può utilizzare https:// selezionando quella casella. Copi quell'URL. Ora vada a Code Spaces. Abbiamo effettivamente Apptainer installato in questo spazio, che è lo stesso di Singularity, ma sono alias l'uno dell'altro. Quindi farò apptainer pull e poi lo chiamerò cowpy sif, ma può chiamarlo come vuole. Incollo l'URL. E questo scaricherà quell'immagine per me.

Potrei fare ls -lh e vedere cowpy.sif

Singularity è diverso da Docker, in quanto singularity memorizza tutte le immagini in file flat, mentre Docker ha un registro dove mantiene tutti i layer separatamente sulla Sua macchina host, e ha un demone in esecuzione per tenere traccia di tutto questo.

## 1.2. Usare il container per eseguire cowpy come comando singolo

Bene, torniamo a Docker. Ora possiamo provare a eseguire questa immagine che abbiamo creato facendo docker run.

Farò dash dash rm, che esegue semplicemente un'esecuzione singola dell'immagine. E incollerò l'URL dell'immagine. E poi infine, termini questo con un comando che vuole eseguire.

L'immagine che abbiamo generato aveva cowpy installato, quindi proviamo cowpy.

Ecco fatto. Ha eseguito il nostro comando. Non ho cowpy installato localmente. Può vedere se provo a eseguirlo, non esiste. Tuttavia, in questo comando, l'ho eseguito usando Docker e ha generato correttamente questo output.

## 1.3. Usare il container per eseguire cowpy interattivamente

Possiamo andare oltre se vogliamo e avviare un container interattivamente e guardarci dentro. Ancora, faccio "docker run dash dash rm". Ora farò dash it, che dice a Docker che vogliamo un terminale interattivo. Faccio di nuovo l'URL dell'immagine, e questa volta, invece di fare cowpy, farò bin bash perché il comando che vogliamo eseguire è bash.

Questo ci porta dentro questo container in esecuzione e può vedere che il prompt è cambiato ora.

Se faccio LS slash può vedere che le directory qui sono diverse.

Se apro un secondo terminale qui sulla destra, che è solo in esecuzione in GitHub Code Spaces e faccio LS slash, vede che abbiamo directory come workspaces e temp, mentre qui in Docker è diverso.

Quindi questo ambiente è completamente separato all'interno di Docker e isolato dal mio ambiente host. Questa è una cosa buona, perché isola l'esecuzione di questo comando nell'immagine Docker e lo mantiene riproducibile tra persone diverse su sistemi host diversi.

Se vuole utilizzare dati dal Suo sistema host all'interno dell'immagine Docker, deve montarlo esplicitamente nel container.

Lo faremo tra un secondo.

## 1.3.2. Eseguire il/i comando/i dello strumento desiderato/i

Prima però, vediamo se possiamo eseguire cowpy. Ancora una volta, il comando è disponibile ora direttamente sulla linea di comando, e possiamo iniziare a fare cose più complesse e passare argomenti. Hello containers e invece della mucca, facciamo il pinguino tux. Vediamo cos'altro abbiamo.

Facciamo cheese. Meraviglioso. Che ne dice di Dragon e Cow? Abbastanza bene.

## 1.3.3. Uscire dal container

Bene. Non posso fare molto di più perché non ho dati in questo container. Quindi usciamo da questa immagine in esecuzione e vediamo se possiamo montare alcuni dati nel container. Posso farlo facendo control D o digitando exit. Bene, sono tornato nel mio normale spazio di codice GitHub.

## 1.3.4. Montare dati nel container

Per montare alcuni dati nel container Docker, devo usare dash V. Quindi prenderò il mio comando docker precedente, tornerò all'inizio faccio dash v. Farò "." per la directory di lavoro locale corrente, e poi due punti per dire dove dovrebbe essere montata nella directory host e faccio slash data. Quindi questo sta montando questa particolare directory nel container in slash data.

Ora se faccio LS slash possiamo vedere che abbiamo una nuova directory chiamata data, e se faccio LS data, può vedere tutti i file che abbiamo nella barra laterale qui. Fantastico.

## 1.3.5. Usare i dati montati

Ora possiamo iniziare a utilizzare alcuni dei file che sono sul sistema host all'interno dell'immagine Docker. Quindi posso dire cat data greetings csv. Se ricorda, questo è il nostro file CSV con i nostri diversi saluti di prima, e posso passarlo a cowpy. Fantastico. Ora stiamo andando da qualche parte.

Bene. Questo è abbastanza per eseguire Docker interattivamente. Spero che ora abbia un'idea di cosa sia approssimativamente Docker e come usarlo sia per eseguire un comando in modo singolo, sia per usare un'immagine interattivamente. Se sta utilizzando singularity. I comandi sono tutti molto simili tranne che fa cose come apptainer exec o apptainer run, o singularity exec o singularity run.

## 2. Usare i container in Nextflow

Successivamente torneremo al nostro workflow Nextflow e vedremo come utilizzare questa tecnologia all'interno del pipeline Nextflow.

Chiudiamo il terminale e apriamo di nuovo Hello Containers.

## 2.1. Scrivere un modulo cowpy

Per restare con il nostro esempio cowpy, creiamo un nuovo processo nel nostro workflow, che utilizza cowpy. Andiamo su moduli, creiamo un nuovo file e chiamiamolo cowpy nf. Ora barare un po' e copierò il codice dal materiale di formazione e premo salva. E diamo un'occhiata.

Quindi questo è un processo semplice. Spero che ora comprenda come appaiono i blocchi di costruzione di un processo. Abbiamo il nostro publishDir di nuovo, che va in results. Abbiamo due input, un file di input e una stringa chiamata character. Abbiamo un output cowpy input file, e abbiamo uno script che sembra esattamente uguale a quello che abbiamo eseguito manualmente all'interno della nostra immagine docker un secondo fa: cat per stampare un file, passandolo a cowpy, dicendo quale tipo di personaggio cowpy vogliamo utilizzare, e producendo quello nel file di output, che passiamo come output qui.

## 2.2. Aggiungere cowpy al workflow

Bene, torniamo al nostro workflow, importiamo questo nuovo processo. Quindi cowpy da modules cowpy nf. Creiamo un nuovo parametro in modo da poter specificare quale personaggio vogliamo. Diciamo Turkey per impostazione predefinita. E poi chiamiamo questo nuovo processo alla fine del workflow,

cowpy. E usiamo l'output qui da Collect Greetings. Quindi collect greetings out, out file qui. E poi abbiamo bisogno di un secondo argomento, che è il nuovo parametro che abbiamo appena creato. params punto character.

## 2.2.4. Eseguire il workflow per verificare che funzioni

Bene, vediamo se il nostro nuovo processo funziona. Nextflow run hello containers. Questo dovrebbe eseguire quei primi tre processi e poi provare a eseguire cowpy alla fine.

Abbiamo un errore. Quello che sta dicendo qui, cowpy ha avuto un errore e aveva uno stato di uscita 127 e infatti, command sh cowpy command not found.

Non abbiamo detto a Nextflow che abbiamo un'immagine Docker disponibile per cowpy, quindi ha provato a eseguirlo sul nostro sistema host e non abbiamo cowpy installato sul nostro sistema host, quindi ha innescato un errore.

## 2.3. Usare un container per eseguirlo

Quindi quello che dobbiamo fare è che dobbiamo dire a Nextflow che abbiamo un container disponibile. Andiamo al nostro processo cowpy e aggiungeremo una nuova direttiva in cima al processo chiamata container.

Quindi troviamo la nostra immagine, copiamo l'URL e lo mettiamo in una stringa.

Questo non è sufficiente di per sé perché un pipeline Nextflow può avere diversi modi per specificare il software. Potrei anche fare conda conda-forge cowpy, per esempio. E Nextflow deve sapere quale di queste tecnologie vuole utilizzare.

## 2.3.2. Abilitare l'uso di Docker tramite il file nextflow.config

Quindi, per eseguire con Docker abilitato, ci anticiperemo leggermente e useremo il file di configurazione Nextflow, che è qualcosa che tratteremo in modo più dettagliato nel prossimo capitolo. Può vedere in questa directory che abbiamo un file chiamato Nextflow Config, e qui ha già docker.enabled False.

Lo cambieremo in True per abilitare Docker, e poi possiamo provare a eseguire di nuovo il workflow.

## 2.3.3. Eseguire il workflow con Docker abilitato

Nextflow run hello containers nf e questa volta cowpy è stato eseguito con successo. Guardiamo in Results. cowpy collected test e c'è il nostro Turkey. Meraviglioso.

Quindi in background lì, Nextflow sapeva di avere un container disponibile per quel processo.

Ha recuperato l'immagine e ha eseguito i comandi per noi.

## 2.3.4. Ispezionare come Nextflow ha lanciato l'attività containerizzata

Se è curioso, possiamo effettivamente vedere esattamente cosa ha fatto guardando nella directory work. Se faccio code work, e poi l'hash e poi command run, che se ricorda è il file effettivo che viene eseguito per quell'attività, possiamo entrare e possiamo cercare una funzione chiamata NXF launch. E qui può vedere l'esatto comando docker che Nextflow ha utilizzato, che assomiglia molto a quello che stavamo facendo manualmente nel terminale prima. Docker run. Collegamento di questa directory host nel container, e poi specificando l'URL del container.

Quindi non c'è magia qui. È solo che Nextflow sta automaticamente facendo il lavoro pesante per Lei in un modo che significa che può facilmente specificare container nel Suo pipeline, che sono poi prontamente disponibili per chiunque altro utilizzi il Suo workflow. E quelle persone non devono più pensare alla gestione del software per eseguire il Suo pipeline di analisi.

Molto, molto semplice, molto conveniente, e anche davvero riproducibile. Buono a tutto tondo.

Bene, ottimo lavoro. Questa è la fine del Capitolo Cinque. Ci raggiunga nel prossimo video per la parte sei, che è la parte finale di questa formazione Hello Nextflow, dove parleremo della configurazione Nextflow in modo più dettagliato.

Ci vediamo nel prossimo video.

[Trascrizione del prossimo video :octicons-arrow-right-24:](06_hello_config.md)
