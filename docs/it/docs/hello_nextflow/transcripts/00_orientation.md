# Orientamento - Trascrizione Video

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduzione assistita da IA - [scopri di più e suggerisci miglioramenti](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/G3CV-FcV-rc?si=nyLvwhrSB2m1NPc5&amp;list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Nota importante"

    Questa pagina mostra solo la trascrizione. Per le istruzioni complete passo dopo passo, torni al [materiale del corso](../00_orientation.md).

## Benvenuto

Salve, benvenuto a Hello Nextflow. Mi chiamo Phil Ewels. Sono Product Manager per Open Source in Seqera, e sono lieto di essere qui oggi per accompagnarla attraverso questo primo corso di formazione su Nextflow.

Esamineremo le basi di Nextflow, spiegando come scrivere ed eseguire pipeline e come configurarle.

E costruirà la sua semplice pipeline multi-step. Tratteremo terminologia come operatori e channel factories, e alla fine del corso, sarà pronto per iniziare a costruire le sue pipeline bioinformatiche.

Se ha domande, la preghiamo di contattarci su community.seqera.io. Abbiamo una comunità Nextflow molto attiva con una sezione dedicata alla formazione, quindi ci faccia sapere dove ha difficoltà e qualcuno sarà in grado di aiutare.

Bene. Iniziamo.

## Sito Web di Formazione

Tutto il materiale di formazione per i corsi Nextflow si trova su training.nextflow.io. Può accedervi nel suo browser web. Lo avvii ora e possiamo dare un'occhiata.

Eseguirò questo con la versione 2.1.1. Pubblichiamo piccoli aggiornamenti e correzioni qua e là, quindi non si preoccupi se è leggermente diverso, ma se il materiale è troppo cambiato, può sempre usare questo selettore di versione in alto per scegliere la versione esatta dei materiali che tratterò.

Se preferisce la modalità chiara, può cambiare il tema del sito web qui.

Veda le traduzioni qui, anche se al momento della registrazione, è disponibile solo in inglese, che copre questo nuovo materiale.

E veda anche tutto il codice sorgente per il sito web di formazione e tutto ciò con cui lavoreremo su GitHub.

La homepage qui elenca tutti i diversi corsi di materiale di formazione che abbiamo. Quindi scorrendo verso il basso, vedremo Nextflow per principianti con il corso Hello Nextflow che faremo qui. Può vedere tutti gli altri corsi che abbiamo anche, che funzionano in modo simile.

## Configurazione dell'Ambiente

In realtà inizierò usando questo primo in alto, che è comune per tutti i corsi di formazione, e riguarda specificamente la configurazione del nostro ambiente.

Cliccando, mi porta a questa sezione, e possiamo vedere le istruzioni per sviluppare localmente. Se vuole usare il suo laptop con la sua copia di VS Code e le sue installazioni software, o quello che ci aspettiamo che facciano la maggior parte delle persone, che è usare qualcosa chiamato GitHub Codespaces.

Codespaces è un servizio fornito da GitHub dove eseguono un server web nel cloud, a cui può connettersi. Quel server ha VS Code installato, dove può eseguirlo nel suo browser web, o se preferisce, connetterlo alla sua installazione locale di VS Code. Tutti i calcoli, tutti i file, tutte le modifiche avvengono in remoto, il che significa che tutto il software di cui ha bisogno viene preinstallato ed è lo stesso per tutti.

## Creazione di un GitHub Codespace

Per creare il codespace con tutto ciò di cui abbiamo bisogno, cerchi i pulsanti nei documenti del materiale, che dicono "Open in GitHub Codespaces". Cliccherò ora, aprendolo in una nuova scheda. E mi viene presentata questa pagina web. Ora può vedere che è preconfigurato per essere impostato con nextflow-io training.

Posso semplicemente cliccare crea nuovo codespace. Ma in realtà raccomandiamo di usare una macchina leggermente più grande per la formazione Nextflow con quattro CPU invece di due. Può cambiare quale versione del materiale utilizza. Quindi questo è predefinito su 2.1.1 perché è la versione dei documenti da cui ho seguito il collegamento. Ma potrei anche impostarlo su un branch specifico del repository se voglio.

Ora cliccherò crea codespace. E inizierà a configurare l'ambiente per me.

## Creazione del Codespace

Ora, la prima volta che lo fa, ci vorrà parecchio tempo, quindi ora è un buon momento per andare a prendere una tazza di tè. Si metta comodo, chiacchieri con la persona seduta accanto a lei.

Se è interessato, può cliccare building codespace qui sotto per vedere i log della configurazione. E può vedere qui che sta scaricando un'immagine Docker con tutto ciò di cui ho bisogno e configurando l'ambiente.

Ora, deve aspettare così solo la prima volta che crea un codespace. Se va su github.com/codespaces qui, vedrà tutti i diversi Codespaces che ha aperti. Ecco quello che ho appena creato. La prossima volta che lo fa, può andare qui e può selezionare il codespace precedente e tornare direttamente ad esso. Ed è un processo molto, molto più veloce per avviare quell'ambiente esistente. Questo manterrà anche tutte le modifiche che ha apportato a VS Code e ai file, quindi non perderà i suoi progressi se esce e torna.

Può cliccare i tre punti qui per eseguire altre azioni. Ad esempio, se l'ha configurato con due CPU e ora ne vuole quattro, può cambiare il tipo di macchina. O se vuole ricominciare da capo e da zero, può eliminare il codespace.

## Introduzione a VS Code

Okay, Codespaces ha finito di configurare il mio ambiente e ora mi viene presentato VS Code nel browser web.

Se è abituato a VS Code. Questo sembrerà molto familiare se non l'ha mai usato prima, è abbastanza semplice. Ci sono alcune parti diverse della pagina di cui deve essere consapevole.

Qui a sinistra, abbiamo la barra laterale. Può vedere l'Explorer configurato con tutti i diversi file nel repository GitHub dal repository di formazione.

Su questi pulsanti in basso a sinistra, possono esserci strumenti diversi. Nella barra laterale. Posso cercare tutti i file in tutto il progetto. Posso lavorare con Git, posso lavorare con GitHub, tutte cose diverse come quella.

In alto qui c'è il menu principale. L'esploratore di file è quello che avremo aperto per la maggior parte del tempo qui, e può fare clic destro su uno qualsiasi di questi file e fare le cose normali che si aspetterebbe. Potrebbe dover cliccare attraverso alcuni avvisi come questo dove dice taglia copia e può anche scaricare sulla sua macchina locale.

Quando il codespace si carica, ci dà un'anteprima del file markdown in quest'area principale qui. Questo è proprio lo stesso di quello che viene visualizzato su github.com. Posso chiuderlo e se faccio doppio clic su quel file Readme, vedrà che lo apre come codice nell'editor di codice e proprio come con qualsiasi altro file, possiamo modificare questo codice direttamente.

Infine qui in basso, abbiamo la finestra del terminale. Stavo guardando i log mentre si costruiva, quindi è quello che sta mostrando la cosa corrente. Posso anche premere questo pulsante più per avviare una nuova sessione di terminale. Questo non è in esecuzione sulla mia macchina. Ricordi, questo è in esecuzione nel cloud, e se faccio tree tre alla profondità di due, vedrà tutti gli stessi file qui, che erano sulla sinistra.

## Mostrare solo i file "hello-nextflow"

Questo repository GitHub contiene tutti i diversi set di formazione, non solo quello che stiamo facendo. Quindi se vuole, può concentrarsi solo sulla cartella Hello Nextflow. Un modo per pulire un po' questo è andare al menu file e poi aggiungere cartella all'area di lavoro.

Facciamo clic su quello, andiamo a training. Hello nextflow, e clicchiamo aggiungi. Aggiornerà il suo schermo. E poi nell'Explorer, ora abbiamo due diverse aree di lavoro, quella che avevamo prima per training e una con solo Hello Nextflow.

Se vuole, può fare clic destro su training e fare clic su rimuovi cartella dall'area di lavoro per eliminarla completamente dalla barra laterale.

Ora abbiamo solo i file per questo particolare corso di formazione nella barra laterale. Posso nascondere quell'avviso e ora posso fare la stessa cosa nel terminale qui e fare CD per cambiare directory. Hello, Nextflow. E di nuovo, abbiamo gli stessi file qui, che sono nella barra laterale.

## Hello Nextflow: file

Guardando questi file per il corso Hello Nextflow.

Abbiamo un gruppo di file .nf, che sono per Nextflow, e c'è uno di questi file per ciascuno dei capitoli del corso di formazione. Lavoreremo su questi file e li modificheremo negli esercizi.

Abbiamo anche un file nextflow.config, che ha solo impostazioni di configurazione di base per eseguire Nextflow in questo ambiente, di cui non deve davvero preoccuparsi a questo punto. Un file greetings.csv, che useremo per elaborare i dati, che verrà introdotto nella prossima parte di questo corso, e un file test-params.json, che verrà utilizzato nella parte sei e che può ignorare per ora.

Questi file Nextflow sono solo l'inizio di ogni esercizio. Se vuole vedere come dovrebbero apparire quando sono finiti, può andare in una directory solutions e ci sono le risposte per ogni parte del corso di formazione, quindi può vedere una versione funzionante di ciò che sta cercando di raggiungere.

## Apertura di un terminale

Se in qualsiasi momento chiude il terminale e non riesce a ricordare come tornare indietro, non si preoccupi. Questi pulsanti in alto a destra aprono e chiudono diversi pannelli nell'area di lavoro. Quindi clicchi questo per il pannello inferiore e riapparirà. E si assicuri solo di aver selezionato terminal qui. Può anche cliccare questo pulsante qui, la freccia sul lato destro del terminale per renderlo a schermo intero.

Mi vedrà farlo molto spesso perché ho VS Code ingrandito in modo che possa leggere il testo. A seconda delle dimensioni del suo schermo, potrebbe o meno aver bisogno di farlo. Lo stesso vale per ridurre al minimo il pannello laterale.

Bene. Questo è abbastanza per l'ambiente. Penso che siamo pronti per iniziare. Mi raggiunga nel prossimo video per il capitolo uno.

[Prossima trascrizione video :octicons-arrow-right-24:](01_hello_world.md)
