# Orientamento - Trascrizione Video

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduzione assistita da IA - [scopri di più e suggerisci miglioramenti](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/PIjOdFaYwWaTV4zd?si=y8lAedhEHWaTV4zd&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Nota importante"

    Questa pagina mostra solo la trascrizione. Per le istruzioni complete passo dopo passo, tornate al [materiale del corso](../00_orientation.md).

## Benvenuti

Ciao e benvenuti a Hello Nextflow. Mi chiamo Phil Ewels. Sono Product Manager per il Software Open Source presso Seqera, l'azienda dietro Nextflow.

Questo corso è un'introduzione pratica alla costruzione di flussi di lavoro con Nextflow. È progettato per persone completamente nuove a Nextflow che vogliono sviluppare le proprie pipeline.

Gli esempi sono tutti semplici elaborazioni di testo, così potete concentrarvi sui concetti di Nextflow senza bisogno di competenze specifiche di dominio, solo una certa familiarità con la riga di comando.

Affronteremo le basi di Nextflow: scrivere processi, collegarli in flussi di lavoro multi-step, gestire le dipendenze software con i container e configurare le pipeline per diversi ambienti di calcolo. Alla fine, avrete costruito una pipeline funzionante da zero.

Questo corso si concentra sullo _sviluppo_ di pipeline. Se volete solo _eseguire_ pipeline esistenti senza addentrarvi troppo nel codice, abbiamo un corso più breve "Nextflow Run" che potrebbe essere più adatto a voi.

Una volta acquisite le basi qui, abbiamo anche corsi di approfondimento che applicano questi concetti ad analisi scientifiche reali. Vi insegneremo come usare le pipeline e le best practice della community nf-core.

Se vi bloccate, andate su community.seqera.io. C'è un forum della community attivo con una sezione dedicata proprio alle domande sulla formazione. Potete usarlo in qualsiasi momento, tuttavia organizziamo anche settimane di formazione trimestrali con persone disponibili specificamente per aiutare. Quindi se state seguendo la formazione durante una di queste, non siate timidi e chiedete aiuto.

Potete anche provare a chiedere aiuto a Seqera AI. È ottimo nello spiegare il codice Nextflow e nell'aiutarvi con il debugging.

Quando siete pronti per eseguire Nextflow su larga scala, Seqera Platform è il posto migliore per farlo. Funziona sulla vostra infrastruttura senza alcun vendor lock-in, con tutto dal lancio delle pipeline al monitoraggio in tempo reale, agli ambienti di analisi interattivi. Ma per ora, concentriamoci solo sui fondamentali.

Bene, iniziamo.

## training.nextflow.io

Okay. La prima cosa da notare è che tutti i corsi di formazione su training.nextflow.io sono molto interattivi. L'idea è che seguiate il materiale di formazione e le mie istruzioni, e affrontiamo insieme il materiale di formazione. Quindi avrete bisogno di due cose: il vostro laptop e questo sito web aperto. E questo è praticamente tutto.

Questa è la homepage come appare oggi mentre registro. Potete vedere che c'è una panoramica delle diverse cose, il background e i diversi corsi che abbiamo, e la lista cresce continuamente.

Nextflow for newcomers è dove siamo. Ci sono due corsi qui dentro, Nextflow Run, che è un corso diverso, e Hello Nextflow, che è quello che ci interessa.

E potete anche vedere tutti i diversi corsi nella barra laterale. Posso passare a Hello Nextflow e possiamo vedere tutti i diversi capitoli che affronteremo insieme.

Ci sono un paio di altre cose importanti da notare qui. Innanzitutto, il materiale di formazione è versionato, quindi potete vedere qui in alto. Dice 3.0 latest, che al momento della registrazione è l'ultima versione stabile. Questo cambierà nel tempo. Pubblichiamo nuovi corsi e aggiorniamo il materiale nel tempo. Quindi se è 3.1 o 3.2, non preoccupatevi troppo. Se è 4.0, allora probabilmente c'è un nuovo video, e dovreste forse andare a cercarlo perché probabilmente ci saranno aggiornamenti significativi.

Un altro menu a tendina in alto è questo, per la lingua. Ora questo è completamente nuovo per la versione 3.0. Abbiamo preso il materiale precedentemente tradotto, che era stato fatto da umani, a mano, e l'abbiamo passato in un LLM e abbiamo configurato tutta questa nuova infrastruttura per mantenere diverse traduzioni del materiale di formazione usando la traduzione LLM.

Quindi ora abbiamo tutte queste fantastiche traduzioni qui. Quindi se volete seguire in coreano, potete caricare l'intero sito web in coreano. E seguire lì. Lo stesso per tutte queste altre lingue, hindi e tedesco e così via. Io seguirò in inglese. Questa è la lingua principale in cui scriviamo il materiale.

Un paio di altri pulsanti se vi piace avere la modalità chiara. Invece della modalità scura, potete seguire il sito web in modalità chiara qui in alto.

E poi anche tutto ciò che guardiamo è in un singolo repository GitHub, che è open source, chiamato nextflow-io/training. E se cliccate questo pulsante in qualsiasi momento, vi porterà al repository GitHub. Torneremo su questo tra un minuto.

## Configurazione di GitHub Codespaces

Okay, quindi ora avete questo aperto nella scheda del browser. Andiamo su Hello Nextflow e clicchiamo. Potete vedere nella pagina introduttiva, ci dice alcuni dei requisiti, la panoramica e il piano della lezione di cosa affronteremo approssimativamente, e poi ci addentriamo nell'iniziare.

Ci sono diversi modi per fare questo tutorial interattivo. Se vi va, siete i benvenuti a farlo localmente sul vostro computer con la vostra installazione di Nextflow. Se clicchiamo su Environment Options, potete vedere che ci sono più dettagli su come farlo usando Devcontainer locali o potete anche semplicemente installare tutto il software localmente, con installazione manuale.

Stiamo lavorando per far funzionare bene questo con Seqera Studios, quindi è un'altra opzione. Ma la più comune al momento è usare GitHub Codespaces.

Codespaces configura un ambiente sandbox su un server remoto gestito da GitHub. Ed è gratuito per una certa quantità di utilizzo, che di solito va bene per la formazione. E vi configurerà con un'istanza VS Code, un IDE dove potete accedere a tutti i file dal repository, eseguire Nextflow e tutto il resto. E abbiamo preconfigurato Codespaces per voi. Quindi ha tutto ciò di cui avete bisogno.

La bellezza di questo è che è solo un clic per configurare un Codespace. È lo stesso per tutti, e sappiamo che avete già tutti i prerequisiti installati, quindi è bello e veloce.

Quindi la prima cosa da fare è andare su "Getting Started". Cercate questo pulsante, che dice _Open in Codespaces_. Farò cmd + clic per aprirlo in una nuova scheda, e ci porta a GitHub.

Ecco come appare. Possiamo vedere che abbiamo impostato tutte le opzioni qui per voi. Se volete, potete cliccare change options. Alcune cose che potete fare qui. Potete dare una macchina istanza più grande, per esempio, se scoprite che va in crash perché finisce la memoria o qualcosa del genere. O impostare versioni specifiche del materiale di formazione. Ma di solito potete semplicemente andare con quello che abbiamo configurato qui e potete vederlo. In questo caso sta usando la release 3.0.

Quindi cliccherò create new Codespace. E questo mi porta dentro.

Notate anche, dice no Codespace to resume lì. Se ho precedentemente creato un Codespace, cliccando di nuovo quel pulsante sul materiale di formazione mi porterà alla stessa pagina e elencherà tutti i Codespace che ho già in esecuzione. Quindi potete semplicemente saltare direttamente in essi e continuare da dove avete lasciato. Quindi non importa se avete chiuso il vostro laptop.

Si spengono automaticamente dopo pochi minuti di inattività, ma non è un problema. Potete semplicemente riavviarli.

Una volta avviato un nuovo Codespace, rimarrà su questa pagina così e caricherà per un bel po'. Quindi ora è un buon momento per fare una breve pausa. Forse avete dimenticato di andare in bagno o volete una tazza di tè prima di iniziare? Andate ora mentre aspettate questo, perché girerà lì per un po'.

Velocemente mentre aspettiamo che carichi, andrò anche su github.com/codespaces e mostrerò questa è la pagina panoramica dove potete vedere tutti i diversi Codespace che avete attualmente in esecuzione.

Potete vedere che ne ho uno qui per nextflow-io/training. Nessuna modifica, perché non ho ancora fatto nulla. La quantità di risorse che sta usando, e potete vedere al momento si sta configurando. Posso andare qui, cliccare questo piccolo menu a tendina e cliccare delete. Quindi se accidentalmente configurate più Codespace e non ne state usando alcuni, potete eliminare quelli vecchi e fare pulizia.

Infine, un altro modo per entrare in questo. Se andiamo al repository GitHub. E questo funziona per qualsiasi repository GitHub. Cliccate code. Avete comandi per clonare il repository localmente. E c'è una scheda chiamata Codespaces. E di nuovo, potete crearne uno nuovo, e potete vedere quelli che sono già in esecuzione.

Quindi di nuovo, se dimenticate come avete creato il vostro Codespace, potete sempre tornarci in questo modo.

## L'interfaccia VS Code

Okay, la costruzione è finita e ora sta iniziando a caricare GitHub Codespaces. Non ci vuole sempre così tanto, quindi non preoccupatevi. È solo la prima volta che create il Codespace. Se tornate in uno che esiste già, è molto più veloce.

Non siate troppo impazienti se questa è la prima volta, non ha ancora finito, anche se sta iniziando a darci un'interfaccia.

Ma mentre aspettiamo che le cose finali siano configurate, vi guiderò attraverso l'interfaccia nel caso non abbiate molta familiarità con VS Code.

Innanzitutto, c'è la barra laterale della chat per le cose AI, che non ci serve. Quindi la chiuderò, me ne libererò e libererò un po' di spazio.

Sulla sinistra, abbiamo l'esploratore file che ci mostra tutti i file nel repository Git, che è lo spazio di lavoro che abbiamo creato. Notate, questi non sono file locali. Questo è tutto sul server remoto dove stiamo lavorando. Potete trascinare e rilasciare file locali e cose del genere, ma per la maggior parte, non ci penseremo oggi. Stiamo lavorando puramente da remoto.

Ci sono altri strumenti in questa barra laterale, per esempio, search. Quindi potete cercare tutti i file in un repository in una volta. E se stessimo facendo lavoro di sviluppo sul repository di formazione, potremmo fare integrazione con il controllo sorgente con Git e debugging e altre cose.

Altre cose sono, c'è una finestra principale di editing del codice qui in alto, che ha appena caricato un'anteprima del readme, che è per il materiale di formazione. Quindi in questo caso sta visualizzando markdown, ma normalmente questo sarà un editor di codice.

E poi sotto abbiamo il terminale, che è dove eseguiremo tutti i nostri comandi e interagiremo direttamente con Nextflow.

Tutto nel Codespace è preinstallato, quindi il comando Nextflow è già lì e così via.

Okay. Quando arrivate a questo punto, dovrebbe essere quasi finito. Potete vedere ora ha scaricato il language server Nextflow e ha configurato alcune estensioni per noi in VS code, inclusa l'estensione Nextflow, che sarà utile. Quindi posso chiudere quello e posso chiudere il README.md.

E ora potete vedere che ho qualcosa in più sulla sinistra. Sono un po' ingrandito qui, ma se riduco lo zoom potete vedere che uno dei pulsanti dice Nextflow con l'icona Nextflow. E ha alcune cose carine qui dentro per esplorare il progetto e cose del genere, a cui torneremo più tardi.

Okay. Nel caso perdiate uno di questi pannelli, questi pulsanti in alto a destra sono davvero utili e questi semplicemente mostrano e nascondono le cose. Quindi quello mostra e nasconde l'Explorer, mostra e nasconde il terminale in basso. E così via.

Userò questi parecchio perché sono molto ingrandito, quindi cerco di aiutarvi a vedere tutto il testo sul mio schermo, ed è utile poter andare a schermo intero con il terminale e poi nasconderlo quando guardiamo il codice. Ma la maggior parte del tempo potete semplicemente avere tutte queste cose aperte allo stesso tempo.

Okay, cos'altro guardare? Non molto altro. Notate che Nextflow, come ho detto, è installato. Quindi posso digitare "nextflow -version" e dovrebbe apparire dicendo quale versione abbiamo installato.

C'è anche altra roba installata qui. Alla fine di ogni capitolo, abbiamo una serie di domande quiz, per esempio, sul sito web. E potete anche farle nel terminale se volete digitando quiz.

Ci sono altre scorciatoie da tastiera che userò, giusto nel caso siate curiosi. Per esempio, proprio allora ho premuto cmd+K sul mio Mac e questo ha pulito il terminale, per sbarazzarsi di tutto l'output precedente. Quindi è bello per tenere le cose pulite. Se mi vedete fare così è come lo faccio.

E anche se siete nuovi al terminale, ricordate che potete usare tab per auto completare, che userò molto per auto completare i percorsi.

Quindi posso vedere sulla sinistra qui c'è una cartella chiamata Hello Nextflow, che è quella su cui lavoreremo. Se faccio "ls" per elencare i file, posso fare "hel", premere tab, auto completa. E quindi questo è un modo molto veloce per completare i percorsi.

## Aprire solo la cartella Hello Nextflow

Okay. Questo è fantastico. C'è molta roba in questo repository però.

Ci sono tutti i file per generare il sito web, e ci sono più corsi diversi qui dentro, e potete farlo da questa radice e semplicemente cliccare nella cartella "Hello Nextflow". Ma è bello concentrarsi puramente su questo.

Potete impostare questo come vostro spazio di lavoro con un po' di clic qui intorno e impostando una directory di progetto e roba del genere. Ma il modo più semplice è digitare code, che è il comando CLI per lanciare VS Code, e poi "hello-nextflow".

Questo aprirà una nuova scheda del browser e potete chiudere quella vecchia. E sembra esattamente lo stesso. Ma ora potete vedere che siamo in questa sottodirectory e tutti gli altri file sono invisibili, e abbiamo una configurazione più pulita.

Potete vedere qui che anche la directory di lavoro corrente è ora all'interno della cartella Hello Nextflow. Quindi bello e pulito. Non dobbiamo preoccuparci di essere nel posto sbagliato. Okay.

## Nuova Sintassi Nextflow per il 2026

C'è una cosa speciale che devo menzionare a questo punto. Proprio ora, all'inizio del 2026, stiamo iniziando a introdurre diverse funzionalità in Nextflow, e una delle grandi novità è un nuovo parser di sintassi del linguaggio all'interno di Nextflow.

Fondamentalmente il motore che legge i vostri file Nextflow e li comprende, per il runtime. Ci sono alcuni cambiamenti alla sintassi, ed è davvero importante che usiate Nextflow con il parser di sintassi corretto abilitato.

Avete bisogno di due cose per questo. Avete bisogno di una versione aggiornata di Nextflow e dovete assicurarvi che sia abilitato.

Se faccio "nextflow -version" di nuovo, vedrete che Codespaces sta eseguendo la 25.10.2 e 25.10 è la versione minima per poter usare questa roba.

Se state usando la 26.04, che per me non è ancora uscita, ma lo farà presto. Allora questa eseguirà il nuovo parser di sintassi per impostazione predefinita, e non dovete fare nient'altro.

Ma se state eseguendo la 25.10, dovete abilitare lo strict syntax parser, come viene chiamato, o v2 syntax parser.

Questo viene fatto con una variabile d'ambiente. È già impostata in Codespaces, quindi non dovete fare nulla. Ma se state eseguendo localmente, dovete impostare questo, e posso verificarlo facendo "echo $NXF_SYNTAX_PARSER", e dovrebbe essere impostato su v2.

Quindi se state eseguendo localmente, fate semplicemente "export NXF_SYNTAX_PARSER=v2". Semplice così. Ma ricordatevi di farlo. Perché altrimenti vedrete alcune strane discrepanze ed errori mentre andiamo avanti.

Se avete qualche dubbio su qualsiasi cosa riguardo alla versione Nextflow e al syntax parser, innanzitutto, ricordate, non dovete preoccuparvi se siete in Codespaces. Tutto dovrebbe essere configurato correttamente. Ma in secondo luogo, se andate sul materiale di formazione Nextflow, se andate giù, parlate dei requisiti di versione, c'è un link qui che vi porta alla pagina di aiuto sulle versioni di esplorazione, e questo spiega tutto in dettaglio.

Vale la pena leggere questo se avete un momento. Perché aiuta a chiarire quali sono alcuni dei diversi termini che potreste sentire quando iniziate a usare Nextflow. Cose come DSL1, DSL2, syntax parser one, syntax parser two, e così via. Quindi vale la pena dare un'occhiata e questo ripete parte di quello che ho appena detto.

È anche davvero utile se avete precedentemente scritto codice Nextflow e state tornando per un ripasso. Vi dice alcune delle cose che cambiano e vi collega a parti della documentazione Nextflow, che vi dice come aggiornare il vostro codice Nextflow.

## File del corso

Okay. L'ultima cosa per familiarizzare è solo vedere i file, che sono in questa directory. Potete guardare nella barra laterale o spesso nel materiale di formazione, usiamo il comando tree, -L, che è il numero di livelli in cui guardare. Diremo due, e se faccio questo a schermo intero, vedrete che questo rispecchia esattamente quello che vediamo sulla barra laterale lì, ma esclude i file nascosti, che iniziano con un punto.

Quindi i file \*.nf, sta per Nextflow. Quindi questi sono i file di script Nextflow, e c'è un file starter qui per ciascuno dei diversi capitoli del materiale di formazione, che apriremo ed esploreremo e poi modificheremo.

Cambieremo questi file man mano che andiamo avanti, e quindi alla fine di ogni capitolo, i file dovrebbero apparire praticamente uguali all'inizio del capitolo per quello successivo. Ma vi diamo questi file diversi così potete sempre iniziare da capo e non preoccuparvi troppo di rovinare la sintassi.

Se avete bisogno di confrontare con qualcosa che dovrebbe sicuramente funzionare. Potete controllare nella cartella solutions, e questo è come uno stato finale per ciascuno dei capitoli, così potete confrontare quello che avete scritto con quello che c'è lì.

C'è una directory data. Questa ha solo un file greetings.csv, che useremo come esempio, dati di input in parte del corso, e cose come un file di configurazione e alcuni parametri, che descriveremo più avanti nel corso.

## Conclusione

Okay, quindi ora si spera che tutto funzioni. Il vostro schermo sembra uguale al mio e capite come accedere a tutto e quali sono tutti i diversi file.

Se scorrete fino in fondo alla pagina su getting started, piccola casella di controllo dovreste dire che capisco cosa sto facendo. Il mio ambiente è attivo e funzionante e avete impostato, la vostra directory di lavoro correttamente sulla cartella "Hello Nextflow".

Se avete spuntato tutti quelli e sembrano verdi. Possiamo continuare al prossimo video e al prossimo capitolo, che è la parte uno. Hello World. Ci vediamo tra un momento.
