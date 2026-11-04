# Parte 4: Hello Modules - Trascrizione Video

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduzione assistita da IA - [scopri di più e suggerisci miglioramenti](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/43Ot-f0iOME?si=0AWnXB7xqHAzJdJV&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Note importanti"

    Questa pagina mostra solo la trascrizione. Per le istruzioni complete passo dopo passo, tornate al [materiale del corso](../04_hello_modules.md).

    I numeri di sezione mostrati nella trascrizione sono forniti solo a scopo indicativo e potrebbero non includere tutti i numeri di sezione presenti nei materiali.

## Benvenuti

Ciao e bentornati alla parte quattro di Hello Nextflow. Questa sezione è tutta dedicata ai moduli, ed è una sezione piuttosto breve del corso. Non scriveremo molto codice, si tratta più di come organizziamo il codice nella nostra pipeline.

Fino ad ora, abbiamo messo tutto in un singolo file, il che va bene, ed è in realtà così che costruivamo le pipeline Nextflow ai vecchi tempi.

Ma man mano che la pipeline cresce, lo script diventa sempre più lungo e sempre più difficile da navigare, da mantenere, e significa anche che non possiamo davvero condividere nessuna parte del codice.

I moduli Nextflow ci permettono di estrarre i processi da quello script principale e poi importarli. Questo significa che il codice è più facile da navigare e significa anche che possiamo condividere quel codice dei moduli tra diverse pipeline.

Questo piccolo diagramma sulla pagina principale della documentazione mostra bene il concetto. Invece di un unico enorme script, includeremo questi file di moduli separati, da diversi script di moduli, e tutto verrà inserito nel flusso di lavoro, ma continuerà a funzionare esattamente nello stesso modo.

Quindi saltiamo in GitHub Codespaces e diamo un'occhiata in giro. Come prima, ho ripulito un po' il mio workspace qui. Ho rimosso le vecchie directory Nextflow e la directory di lavoro e così via. Ma non importa se avete ancora quei file in giro.

Comincerò a lavorare nel file hello modules, che è fondamentalmente dove l'abbiamo lasciato alla fine del capitolo precedente. Abbiamo i nostri tre processi qui dentro. Abbiamo un paio di parametri, il blocco workflow, dove eseguiamo quei tre processi e li colleghiamo insieme con i canali. Poi pubblichiamo i canali di output e abbiamo il blocco output che dice come pubblicare quei file.

## 1. Creare una directory per memorizzare i moduli

Ora, come ho detto, non scriveremo o modificheremo molto codice. Sposteremo semplicemente il codice che abbiamo già. I file dei moduli Nextflow tipicamente hanno un singolo processo al loro interno, e per convenzione normalmente li conserviamo in una directory chiamata modules. Ma potete chiamarla come volete. Ma io manterrò una directory modules nel mio repository qui, e poi creerò un file per ogni processo. Quindi dirò nuovo file, sayHello.nf.

## 2. Creare un modulo per sayHello()

Ora prenderò il mio processo e semplicemente selezionerò questo codice, lo taglierò dal file principale hello modules e lo incollerò qui.

Ovviamente questo non fa nulla da solo. Il nostro script principale ha ancora bisogno di quel processo, quindi dobbiamo riportarlo dentro in qualche modo. E lo facciamo con l'istruzione include.

Quindi scrivo include e delle parentesi graffe, e poi prendo il nome del processo. E dico from, e poi gli do un percorso di file relativo. Quindi dice, inizia con ./ perché è relativo da dove questo script è salvato. Quindi è modules sayHello.nf.

Notate che l'estensione VS Code è piuttosto utile qui. Mi dice, se riesce a trovare questo file e se riesce a trovare un processo, che sto nominando. Se metto deliberatamente un errore di battitura qui, mi dà subito un errore e mi dirà che non riesce a trovare questo processo che sto cercando di importare. Quindi tenete d'occhio eventuali errori che trovate.

E questo è davvero tutto. Abbiamo ancora il nostro processo qui. Non sono necessarie modifiche qui sotto. Il processo ha lo stesso nome ed è eseguito esattamente nello stesso modo. È solo che il codice effettivo del processo è ora in un file separato.

Possiamo eseguire di nuovo il flusso di lavoro Nextflow, funzionerà esattamente nello stesso modo. E questo è fondamentalmente il resto di questo capitolo del corso, spostare questi tre processi nei loro file.

Quindi facciamolo ora. Creerò rapidamente un nuovo file di modulo per il secondo processo: convertToUpper.nf. Taglierò quel codice, lo incollerò lì. E poi includerò quello. Perfetto.

E poi creerò un nuovo file per collectGreetings.nf. Taglio quello.

Un sacco di tagliare, tagliare e copiare e incollare.

E ora il nostro script principale del flusso di lavoro sembra improvvisamente molto, molto più corto, molto più accessibile e molto più facile da leggere.

E potete vedere come il progetto ora inizia a costruirsi con i nostri diversi file. Possiamo immergerci nei dettagli nei punti che vogliamo. Navigare per trovare passaggi specifici nella pipeline molto più facilmente, e ottenere rapidamente una panoramica di ciò che sta facendo la pipeline.

## Navigare i moduli con VS Code

Ora, ovviamente, lo svantaggio di fare questo è che se avete una grande pipeline, avrete molti file di moduli e potrebbero essere organizzati in più sottodirectory o ogni genere di cose. Ora, ancora, un piccolo suggerimento qui. L'estensione VS Code è piuttosto brava a navigare la vostra base di codice per voi e anche a dirvi del codice lì.

Potete vedere che VS Code capisce cos'è questo processo e mi dà una piccola panoramica quando passo sopra, così posso vedere senza dover andare a cercare il codice sorgente, quali sono gli input e gli output, che è tipicamente la cosa più importante quando lo uso in un flusso di lavoro.

E anche se tengo premuto command, sono su un Mac, e clicco sul nome del processo, apre il file direttamente subito. Lo carica. Quindi posso saltare direttamente lì senza nemmeno pensare a quali siano i percorsi effettivi dei file. E questo funziona ovunque, posso farlo anche dove i processi vengono chiamati. Quindi questo rende davvero veloce.

## 4.4. Eseguire il flusso di lavoro

Ok, controlliamo solo che la pipeline funzioni ancora come ci aspettiamo. Quindi apro il terminale. Facciamo "nextflow run hello modules", e vediamo se si esegue senza problemi.

Si spera che tutto il punto di questo sia che la pipeline è fondamentalmente invariata, quindi non dovreste davvero vedere cambiamenti rispetto a quando l'abbiamo eseguita prima. L'output qui sembra esattamente lo stesso, e potete vedere la nostra directory results con tutti gli stessi file, quindi è fantastico. Nessun cambiamento è buono.

## Una nota su nf-core/modules

Prima di concludere, voglio toccare rapidamente il potere della collaborazione quando si tratta di moduli. Questi file sono nel mio repository, quindi non è ovvio subito come potremmo collaborare su di essi. E ci sono molti modi diversi in cui potete farlo, ma probabilmente l'esempio più grande e più conosciuto di questo è nf-core.

Se vado sul sito web di nf-core, vado su resources, e modules. Potete vedere che nf-core ha un'enorme libreria di moduli, quasi appena sotto 1700 moduli quando lo visualizzo. E quindi posso digitare il nome di uno qualsiasi dei miei strumenti preferiti, andare a vedere se qualcun altro ha già scritto un modulo per esso, e vedere questo processo di modulo pre-scritto qui con tutti gli input, gli output, i container software, tutte queste informazioni, e potete vedere sul lato qui, quante diverse pipeline nf-core stanno tutte usando questo singolo processo condiviso.

Questo è un esempio un po' estremo, ma potete vedere che questo sta davvero riutilizzando questo codice. E se clicco attraverso al sorgente GitHub per questo, è esattamente lo stesso di quello che stiamo facendo. È solo un grande processo in un file.

Ora dal lato nf-core, facciamo alcuni trucchi per essere in grado di condividere quei file e portarli in diversi repository. E se volete saperne di più su questo, andate a dare un'occhiata al corso che abbiamo sull'uso e la costruzione con nf-core specificamente. Ma volevo solo darvi un'idea di quanto possa essere potente questo concetto di riutilizzo del codice.

## Conclusione

Bene, questo è tutto per i moduli. Vi avevo detto che era una sezione breve del corso. Date un'occhiata al quiz, assicuratevi di capire e assicuratevi che tutto funzioni ancora correttamente. E ci vediamo nel prossimo video, che è tutto sui container software. Grazie mille.
