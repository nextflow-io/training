# Parte 6: Hello Config - Trascrizione

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduzione assistita da IA - [scopri di più e suggerisci miglioramenti](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/IuDO2HeKvXk?si=tnXTi6mRkITY0zW_&amp;list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Note importanti"

    Questa pagina mostra solo la trascrizione. Per le istruzioni complete passo-passo, tornare al [materiale del corso](../06_hello_config.md).

    I numeri di sezione mostrati nella trascrizione sono forniti solo a scopo indicativo e potrebbero non includere tutti i numeri di sezione presenti nei materiali.

## Benvenuto

Salve, benvenuti alla parte sei del corso di formazione Hello Nextflow.

Questo capitolo si chiama Hello Config, ed è la parte finale del nostro corso di formazione.

In questo capitolo parleremo della configurazione di Nextflow. La configurazione di Nextflow è davvero potente. Ci permette di eseguire la stessa pipeline su diverse infrastrutture di calcolo con diversi sistemi di provisioning del software e diverse opzioni nella pipeline stessa.

Questo significa che è possibile prendere pipeline Nextflow create da altre persone ed eseguirle sul proprio sistema, anche se potrebbero essere state create per un'infrastruttura completamente diversa. Questa capacità di configurare Nextflow rende i flussi di lavoro veramente portabili e condivisibili.

In questo capitolo utilizzeremo il flusso di lavoro che abbiamo costruito nelle parti precedenti, ma non modificheremo affatto il codice del flusso di lavoro. Esamineremo solo il nostro file di configurazione Nextflow e vedremo come la modifica della configurazione altera il modo in cui Nextflow viene eseguito.

Bene, iniziamo.

Come prima, iniziamo andando su training.nextflow.io. Andiamo sulla sinistra su Hello Nextflow e capitolo sei, Hello config. Ora entrerò nel mio ambiente GitHub Codespaces e controllerò lo script che utilizzeremo.

## 0. Preparazione: Verificare che Docker sia abilitato ed eseguire il flusso di lavoro Hello Config

Questo si chiama Hello Config, e parte da dove eravamo prima. Quindi appare esattamente uguale con i nostri tre parametri: greetings per il file CSV, batch per il nome del batch di output e character per il nome cowpy. Abbiamo i nostri quattro import dei diversi processi, e poi abbiamo un flusso di lavoro dove li concateniamo insieme.

Sto per chiudere questo file ora perché non toccheremo affatto il file Nextflow in questo capitolo. Lavoreremo esclusivamente all'interno del file di configurazione. Se guardo il file Nextflow.config che abbiamo brevemente esaminato nel precedente capitolo cinque, possiamo vedere che abbiamo una singola istruzione qui: Docker enabled equals true, che sta dicendo a Nextflow di usare Docker quando esegue questo flusso di lavoro.

Sto usando Nextflow.config nella root della pipeline qui, che viene caricato automaticamente quando eseguo Nextflow. Ma ricordate, Nextflow può caricare file di configurazione da più posizioni.

Se controllo con i documenti Nextflow vado su Configuration, potete vedere un elenco di questi luoghi e una priorità in cui vengono caricati.

Bene. Verifichiamo che il nostro flusso di lavoro si stia eseguendo come previsto. Apro un terminale. Faccio Nextflow. Run. Hello, config. E premo invio. Dovremmo avere quei quattro processi in esecuzione, che terminano con un comando cowpy. Infatti, ha funzionato correttamente. Avevo Docker abilitato, ha scaricato Docker ed ha eseguito cowpy per me, proprio come ha fatto alla fine del capitolo cinque.

## 1. Determinare quale tecnologia di packaging del software utilizzare

Bene. Diciamo che sto eseguendo su un HPC e non ho Docker installato. La cosa migliore da fare in questo scenario sarebbe usare Singularity o Apptainer. Se dovessi farlo, andrei nel modulo cowpy e cambierei questo container per usare l'immagine singularity come ho mostrato nel capitolo precedente, con un oras://, che potete anche ottenere da Seqera Containers.

Andrei poi su Nextflow.config impostando Docker enabled a false e farei singularity enabled equals true. O, se si usa Apptainer, apptainer enabled equals true e funzionerebbe.

Nextflow supporta anche altre tecnologie oltre ai container, qualcosa con cui potreste avere familiarità è conda. Qui possiamo fare conda enabled equals true e impostare Docker a false. conda non usa la stessa direttiva container. Invece, possiamo aggiungerne una nuova qui chiamata conda. Specifichiamo quindi il pacchetto conda che vogliamo usare. È buona prassi essere il più specifici possibile per cercare di rendere la pipeline il più riproducibile possibile. Quindi specificherò il canale conda, conda-forge, e poi cowpy, e la versione esatta, che era 1.1.5.

Potrei anche scrivere semplicemente cowpy se volessi, ma potrebbe risolversi in una versione diversa di cowpy in diverse esecuzioni della pipeline.

La cosa bella di questo è che non ho toccato affatto la direttiva docker. Questa immagine Docker è ancora lì. Sto solo fornendo due alternative ora, e queste possono essere attivate o disattivate usando solo un file di configurazione.

## 1.3. Eseguire il flusso di lavoro per verificare che possa utilizzare Conda

Conda è ora abilitato, quindi proviamolo.

Ottimo. È in esecuzione e potete vedere che c'è un messaggio da Nextflow qui che dice che Nextflow sta creando un ambiente conda per me, e sta usando questa posizione di cache.

In background, Nextflow sta eseguendo comandi "conda create" per me per creare un nuovo ambiente conda isolato con solo i pacchetti che voglio, e poi installa e recupera quei pacchetti conda in modo che possa eseguire il processo.

Potete vedere che ci è voluto un po' di tempo lì perché stava creando l'ambiente e installando il software per la prima volta. Tuttavia, ha memorizzato questo ambiente nella cache, quindi se eseguo di nuovo lo stesso comando Nextflow, dovrebbe essere molto più veloce perché riutilizzerà lo stesso ambiente conda.

Una delle cose belle di questo è che queste direttive possono essere specificate a livello di processo, non solo per l'intero flusso di lavoro. Quindi, se volete, potete mescolare e abbinare quale tecnologia viene utilizzata per diversi processi.

## 2. Allocare risorse di calcolo con le direttive di processo

Il file di configurazione Nextflow può fare molto di più che solo il packaging del software. Possiamo anche dire a Nextflow come eseguire effettivamente i passaggi nella pipeline. Un esempio è dire a un sistema host quali risorse dovrebbero essere rese disponibili per ogni attività in esecuzione.

Per impostazione predefinita, Nextflow non fornisce molto. Fornisce una singola CPU e solo due gigabyte di memoria per ogni processo.

Questo è probabilmente qualcosa che vorremmo cambiare, in modo che i processi che richiedono molto tempo per essere eseguiti possano avere più risorse ed eseguirsi più rapidamente, ma può essere difficile sapere cosa allocare a un processo. Nextflow ha alcuni bei trucchi che possono aiutarvi in questo.

## 2.1. Eseguire il flusso di lavoro per generare un report di utilizzo delle risorse

Eseguiamo di nuovo il flusso di lavoro. Questa volta, aggiungerò un argomento aggiuntivo, che è dash with reports. È un'opzione principale di Nextflow, quindi è un singolo trattino. E poi qualsiasi nome di file mi piaccia. In questo caso, lo chiamerò report config one html.

Eseguirò di nuovo il flusso di lavoro. Verrà eseguito esattamente come prima, ma mi darà un report di supporto aggiuntivo, che potete vedere è ora apparso qui nella barra laterale.

Farò clic destro su questo file, cliccherò download, che lo scarica da GitHub Codespaces al mio sistema locale, in modo da poterlo poi visualizzare facilmente nel browser web qui sopra.

Questo report può essere generato per qualsiasi esecuzione Nextflow, e contiene molte informazioni. Inizia in alto con alcuni metadati sul comando utilizzato, quando il flusso di lavoro è stato eseguito, quanto tempo ha impiegato, ma scorrendo verso il basso, otteniamo informazioni più dettagliate sulle risorse che sono state utilizzate da ogni passaggio della pipeline.

Poiché ogni processo viene eseguito più volte per attività diverse, abbiamo un box plot che mostra la variazione delle risorse che abbiamo utilizzato per ogni processo.

Se scorro un po' più in basso, vedo informazioni simili sulla memoria utilizzata e sulla durata del lavoro. Anche lettura e scrittura del disco.

Potete immaginare per una grande pipeline con attività di lunga durata, questo può essere molto informativo su come mettere a punto la configurazione delle risorse che state richiedendo in modo da non richiederne troppe, ma anche in modo da poterne fornire abbastanza che si esegua rapidamente.

Se continuo a scorrere il report, vediamo anche una tabella delle attività, che ci mostra informazioni dettagliate su ogni singola attività che è stata eseguita nel flusso di lavoro. Questo include informazioni come lo script risolto, che è stato eseguito.

Bene, torniamo al nostro file di configurazione. Abbiamo visto che non avevamo davvero bisogno di molto per il nostro flusso di lavoro, quindi diciamo a Nextflow che abbiamo bisogno solo di un gigabyte di memoria per ogni processo nel flusso di lavoro.

Ora quando lo definiamo così a livello di processo, questo viene applicato a ogni singolo processo nella pipeline.

## 2.3. Impostare le allocazioni di risorse per un singolo processo

Per il gusto della discussione, fingiamo che cowpy stia davvero facendo molto lavoro pesante e abbia bisogno di più risorse rispetto alle altre attività. Possiamo definire un blocco extra di configurazione qui, che si applica solo a quel processo usando, with name cowpy.

Questo è chiamato selettore di configurazione, e possiamo definire diversi pattern qui per corrispondere a diversi processi. Per esempio, potrei fare cow star. Poi seguo con alcune parentesi graffe e diamogli due gigabyte di memoria invece di uno e diciamo due CPU.

Ora Nextflow darà a ogni processo nel flusso di lavoro un gigabyte tranne questa richiesta, che è più specifica. Quindi la sovrascrive. E solo per qualsiasi processo che si chiama cowpy, otterrà due giga di memoria e due CPU.

Notate che Nextflow è intelligente riguardo all'utilizzo delle risorse. Quindi se iniziate a mettere questi numeri a valori più alti, vedrete che Nextflow inizia a mettere in coda le sottomissioni di lavoro una dopo l'altra, piuttosto che eseguirle tutte in parallelo, in modo che non richieda eccessivamente le risorse disponibili.

## 2.4. Eseguire il flusso di lavoro con la configurazione modificata

Proviamo ad eseguire di nuovo il flusso di lavoro e salviamo un nuovo report questa volta.

Bene, possiamo scaricare questo file e dare un'occhiata.

Sì, non sorprendentemente, sembra praticamente esattamente lo stesso perché questo è un flusso di lavoro fittizio, che non sta facendo nulla di reale. Ma potete immaginare come questo approccio iterativo di definizione dei limiti e di esecuzione di flussi di lavoro reali con questo tipo di reportistica vi permetta di fare un approccio basato sull'evidenza per impostare una configurazione appropriata e sfruttare davvero al massimo le risorse computazionali che avete a disposizione.

Potete iniziare ad essere davvero intelligenti riguardo a questo. Nextflow ha una capacità integrata di riprovare in caso di fallimenti, e potete sfruttare nel vostro file di configurazione utilizzando una closure come questa e impostando dinamicamente le risorse che vengono rese disponibili. Quindi qui ho detto a Nextflow di moltiplicare quei due gigabyte per il tentativo di retry. Quindi il secondo retry otterrà quattro giga, il terzo retry otterrà sei giga e così via. Questo va un po' oltre lo scopo di questo corso di formazione, ma se siete interessati, consultate i documenti Nextflow, che hanno una bella sezione sulla logica di retry dinamico.

## 2.5. Aggiungere limiti di risorse

Ora, una cosa che potreste notare su questo è che questo tipo di cosa potrebbe rendere piuttosto facile andare accidentalmente oltre le risorse disponibili sul vostro sistema. Se richiedete più risorse di quelle disponibili Nextflow genererà un errore sulla vostra configurazione e interromperà l'esecuzione. Per evitare ciò, potete usare qualcosa chiamato limiti di risorse.

Sotto l'ambito del processo, nel nostro flusso di lavoro, possiamo definire limiti di risorse come questo, che prende un array, e possiamo specificare la memoria massima, le CPU e il tempo che sono disponibili su questo sistema.

Impostare valori alti qui non aumenta la quantità di risorse richieste. Useremo ancora un gigabyte nelle nostre richieste, ma significa che se una qualsiasi di queste richieste arriva a 750, raggiungerà quel limite massimo e non verrà richiesto nulla di più, il che significa che Nextflow continuerà a funzionare e non si bloccherà a causa di risorse non disponibili.

Quindi questa è una bella salvaguardia da usare, specialmente se state usando una logica dinamica con la vostra allocazione di risorse.

L'altra situazione in cui questo è davvero utile è se state usando pipeline che sono pubbliche e non controllate da voi. Potrebbero venire con valori predefiniti di configurazione, e Nextflow prenderà automaticamente l'approccio giusto di sogliare qualsiasi richiesta di risorse per eseguire sul vostro sistema.

Bene, ottimo. Abbiamo parlato di software. Abbiamo parlato di allocazione di risorse, e abbiamo descritto diversi ambiti di configurazione, sia per tutti i processi che per processi specifici.

## 3. Utilizzare un file di parametri per memorizzare i parametri del flusso di lavoro

Bene, ora rivolgeremo la nostra attenzione ai parametri. Possiamo definire parametri nel file di configurazione proprio come abbiamo fatto prima nello script Nextflow. Quindi params.greeting equals hello o usare l'ambito params e impostare foo equals bar.

E questo è ottimo per impostare i valori predefiniti per il vostro flusso di lavoro. Tuttavia, quando si eseguono pipeline, può essere carino specificare parametri in un file JSON o YAML.

Usare un file come questo è molto meglio che specificare opzioni da riga di comando con dash dash. Poiché quando si esegue un flusso di lavoro, potrebbe essere necessario specificare molti parametri e può essere noioso scriverli tutti in una singola CLI e soggetto a errori. Inoltre, è improbabile che ricordiate tutti i parametri che avete usato, quindi se li codificate in un file, è più facile avviare di nuovo il flusso di lavoro, usando gli stessi parametri in futuro.

Abbiamo un file di esempio qui chiamato test params, e potete vedere che specifica i tre parametri che abbiamo nel nostro flusso di lavoro con tre valori diversi. Personalmente, trovo YAML più facile da scrivere rispetto a JSON. Quindi solo per dimostrare che funziona, creerò un nuovo file chiamato Test yaml e copierò questi dentro, eliminerò le virgolette. E salverò.

Questi file JSON e YAML possono essere più facili da scrivere poiché hanno una sintassi più familiare. Ma notate che questi sono solo per parametri e accettano solo sintassi chiave-valore come questa.

## 3.1. Eseguire il flusso di lavoro utilizzando un file di parametri

Proviamolo. Faccio lo stesso comando di prima. Elimino il report e farò dash params file test params yaml.

No, questa è un'opzione principale di Nextflow, quindi è un singolo trattino.

Bene. Ha eseguito il flusso di lavoro e ha usato i parametri in quel file YAML invece che io specifichi tutti sulla riga di comando. Potrebbe sembrare eccessivo solo per questo semplice esempio, ma potete immaginare se avete 10 o 20 parametri diversi, può essere una seccatura digitare manualmente, e questo è semplicemente molto più facile da modificare in un editor di codice e conservare per motivi di riproducibilità.

## 3. Determinare quale executor dovrebbe essere utilizzato per eseguire il lavoro

Bene. Abbiamo parlato di packaging del software con Docker e conda. Abbiamo parlato dei requisiti di risorse di processo con CPU e memoria. E abbiamo parlato un po' di come specificare parametri quando si eseguono flussi di lavoro.

Le parti finali della configurazione riguardano davvero l'esecuzione, l'infrastruttura di calcolo sottostante stessa, e questo è il vero gioiello della corona di Nextflow: che possiamo eseguire questi stessi flussi di lavoro su più diverse infrastrutture di calcolo.

In realtà passerò al materiale di formazione scritto per un secondo. In questa parte della formazione, possiamo vedere alcuni esempi diversi di come diversi executor, in questo caso, schedulatori HPC, definiscono i requisiti di risorse necessari per sottomettere un lavoro.

Quindi per Slurm, avete queste intestazioni SBATCH, che definiscono dash dash mem e il numero di CPU. Se state usando PBS, avete intestazioni diverse, e se usate Grid Engine, avete ancora intestazioni diverse.

Potete immaginare che sia ancora più diverso se volete eseguire sul cloud, sia esso AWS batch, Google Cloud, Azure o altro.

Ognuna di queste infrastrutture di calcolo sottostanti è chiamata executor e Nextflow sa come comunicare con tutti questi diversi executor per sottomettere lavori con la sintassi corretta.

La buona notizia è che non dovete sapere di questo. Tutto quello che dovete fare è dire a Nextflow, quale executor usare.

## 3.1. Indirizzare un backend diverso

Torniamo al nostro file di configurazione e al processo facciamo executor, e scriverò local.

Local è in realtà quello predefinito, se non specificate nessun altro executor, local è quello che verrà usato, e questo significa semplicemente il vostro sistema host, ovunque abbiate lanciato Nextflow,

Potrei specificare invece, Slurm. E questo sottometterebbe lavori Slurm, oppure potrei dire AWS batch, e questo sottometterebbe lavori ad AWS batch.

Avete bisogno di una configurazione aggiuntiva in alcuni casi, per esempio, eseguire sul cloud avrà bisogno di certe credenziali, ma davvero questo è il nucleo di esso, e può essere semplice come una o due righe di configurazione per eseguire il vostro flusso di lavoro in un ambiente di calcolo completamente diverso.

Anche se stiamo eseguendo su un sistema semplice all'interno di Codespaces, posso ancora giocare un po' con questo e fingere che stiamo eseguendo su Slurm. Se poi lancio di nuovo il flusso di lavoro, Nextflow run, hello config. Fallirà perché non sarà in grado di sottomettere lavori a Slurm. Ma possiamo ancora andare nelle directory di lavoro e vedere cosa ha fatto Nextflow. Quindi se andiamo in questa directory di lavoro e guardiamo Command Run. Potete vedere in cima a questo file, ora abbiamo queste righe di intestazione sbatch, che hanno provato a specificare le risorse necessarie per il lavoro Slurm.

## 4. Utilizzare i profili per selezionare configurazioni preimpostate

Bene, ci siamo quasi. La parte finale di questo capitolo parla dei profili di configurazione. Se state eseguendo la vostra pipeline su diversi sistemi, potrebbe essere fastidioso avere tutti questi diversi file di configurazione Nextflow, che dovete specificare ogni volta.

Invece, potete codificare raggruppamenti di configurazione all'interno del vostro file di configurazione Nextflow, e attivare o disattivare quei gruppi usando un flag di profilo. Vediamo come appare.

## 4.1. Creare profili per passare tra sviluppo locale ed esecuzione su HPC

Creeremo due profili nel nostro esempio qui, uno per il mio laptop e uno per un sistema HPC più pesante. Barò un po' e copierò semplicemente il codice dal materiale di formazione e lo metterò qui.

Abbiamo un nuovo ambito chiamato profiles, e poi abbiamo un nome per ogni profilo, che può essere qualsiasi cosa. E all'interno di quello abbiamo una configurazione, che sembra esattamente la stessa della configurazione di primo livello che abbiamo già scritto. Quindi di nuovo, abbiamo l'ambito del processo. L'ambito di Docker.

Sul profilo chiamato my laptop. Sto dicendo di eseguire usando l'executor local, quindi sul mio sistema host e di usare Docker.

Sul profilo university HPC qui sto dicendo di usare Slurm per sottomettere lavori, di usare conda invece di Docker, e sto specificando diversi limiti di risorse, che potrebbero corrispondere alla dimensione del sistema dei nodi sull'HPC che sto usando.

Per impostazione predefinita, nessuna di questa configurazione verrà utilizzata quando eseguo Nextflow, devo specificare che voglio usare uno di questi profili.

## 4.2. Eseguire il flusso di lavoro con un profilo

Facciamo nextflow run hello config. E farò dash profile, singolo trattino perché è un'opzione principale di Nextflow. E poi il nome che gli ho dato, che è my laptop. Nextflow dovrebbe ora usare il blocco di configurazione che è stato specificato all'interno di quel profilo di configurazione, e applicarlo quando esegue Nextflow. Se volessi usare l'altro blocco di configurazione, devo solo cambiare quel nome di profilo. Molto più facile da ricordare. Molto più facile da usare.

## 4.3. Creare un profilo di test

Notate, i profili possono avere qualsiasi tipo di configurazione, quindi non deve essere correlato al vostro ambiente di esecuzione. Per esempio, creiamo un nuovo profilo qui, che ha un insieme di parametri. Possiamo cambiare questo in tux e cambiare in my profile, e ora quando facciamo profile test, specificherà questi parametri, che sovrascriveranno i parametri specificati al livello superiore del flusso di lavoro.

Quando eseguite Nextflow, potete concatenare più profili e verranno applicati in sequenza.

## 4.4. Eseguire il flusso di lavoro localmente con il profilo di test

Quindi posso prendere il comando precedente e fare virgola test. Questo applicherà la configurazione my laptop per prima, e poi applicherà la configurazione test. Se c'è qualche sovrapposizione, allora il profilo a destra sovrascriverà qualsiasi configurazione nei profili precedenti. Se premo invio, vediamo cosa succede.

Bene, abbiamo un nuovo file di risultati qui. Potete vedere My Profile, che ho specificato come una delle opzioni. E possiamo anche vedere cowpy, my profile, e infatti, c'è tux. Quindi ha funzionato.

## Conclusione

Bene! Fantastico. È tutto. Siete arrivati alla fine del corso. Avete un po' di coriandoli di celebrazione. Ben fatto per aver finito questo capitolo.

[Trascrizione del video successivo :octicons-arrow-right-24:](07_next_steps.md)
