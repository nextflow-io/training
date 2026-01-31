# Parte 4: Hello Modules - Trascrizione

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduzione assistita da IA - [scopri di più e suggerisci miglioramenti](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/Xxp_menS0E8?si=0AWnXB7xqHAzJdJV&amp;list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Note importanti"

    Questa pagina mostra solo la trascrizione. Per le istruzioni complete passo dopo passo, torni al [materiale del corso](../04_hello_modules.md).

    I numeri delle sezioni mostrati nella trascrizione sono forniti solo a scopo indicativo e potrebbero non includere tutti i numeri di sezione presenti nei materiali.

## Benvenuto

Salve, benvenuto alla Parte Quattro del corso di formazione Hello Nextflow.

Questo capitolo si chiama Hello Modules, e parleremo di come modularizzare il codice Nextflow. Ciò che faremo è prendere il nostro script di workflow unico e dividerlo in file separati.

Questo rende il codice più facile da navigare e mantenere man mano che il workflow diventa più grande, e rende anche possibile condividere moduli tra pipeline in modo che se ha più pipeline che utilizzano lo stesso strumento, debba scrivere quel process solo una volta.

Un esempio classico di questo è il repository dei moduli nf-core, che ha migliaia di diversi strumenti in moduli pronti all'uso, che può installare e utilizzare nel suo workflow.

Nextflow può anche lavorare con sub workflow, che sono come moduli, ma hanno più processi. Questo è al di fuori dello scopo di questa formazione, ma funziona sostanzialmente nello stesso modo.

Bene. Diamo un'occhiata.

Come al solito, inizi andando su training.nextflow.io.

Vada su "Hello Nextflow" nella barra laterale, e stiamo facendo la parte quattro: "Hello Modules".

Ora sto per saltare nel mio ambiente GitHub Code Spaces e dare un'occhiata al file "hello-modules".

Proprio come prima, stiamo partendo dal punto finale del capitolo precedente, quindi questo script dovrebbe risultare familiare. Abbiamo i nostri tre processi, say hello, convert to upper e collect greetings, e in un workflow semplice, che esegue questi tre comandi ed emette un messaggio alla fine. Abbiamo due parametri chiamati greeting e batch, che specifica il nome, che viene utilizzato per il file di output raccolto alla fine.

## 0. Riscaldamento: Eseguire hello-modules.nf

Possiamo verificare che questo workflow funzioni ancora come ci aspettiamo facendo nextflow run hello modules.

Ottimo. Ha eseguito tre attività con ciascuno di questi processi, un'attività collect, e ci ha detto che ci sono tre saluti in questo batch. Se andiamo in results, abbiamo i nostri diversi file di output qui, incluso l'output di test batch raccolto.

## 1. Creare una directory per memorizzare i moduli

Giusto. Facciamo un po' di modularizzazione.

È generalmente una buona idea mettere i moduli in una sottocartella nel repository della pipeline, solo per mantenere le cose ordinate. Può chiamarla come vuole, ma per convenzione di solito la chiamiamo modules.

Quindi procediamo, andiamo in un terminale e facciamo make the modules. Può vederla apparire nella barra laterale e VS Code qui.

## 2. Creare un modulo per sayHello()

Poi creerò un nuovo file per il mio primo modulo. Può fare "touch" o "code" o può farlo nella barra laterale, non importa davvero. Quindi farò code modules e lo chiamerò come il process. Quindi sayHello.nf. NF è un'estensione di file tradizionale per i file Nextflow.

Salverò qui e possiamo vedere che il nostro nuovo file di modulo appare.

## 2.2. Spostare il codice del process sayHello nel file del modulo

Bene, successivamente prenderò il codice del modulo dal workflow. Prenderò anche lo shebang qui e lo copierò per primo in modo che sia chiaramente un file Nextflow. E poi prenderò questo process e lo taglierò. Quindi lo rimuoverò dal mio script di workflow principale e lo incollerò in questo nuovo modulo.

Questo è tutto il contenuto che questo file di modulo conterrà. Solo un singolo process, nessun workflow, nessuna logica, solo un process da solo.

Ora posso chiudere questo file.

## 2.3. Aggiungere una dichiarazione import prima del blocco workflow

Ora al mio workflow manca quel primo process, quindi dobbiamo riportarlo importandolo. La sintassi per questo è molto simile ad altri linguaggi di programmazione, quindi potrebbe risultare familiare. Facciamo include parentesi graffe, il nome del process, say hello, e poi from il percorso del file modules, say hello, nf. Fantastico.

Un paio di trucchi qui. L'estensione VS Code è intelligente riguardo a questo. Riconosce questo percorso del file e può passarci sopra con il mouse e fare follow link. Oppure sono su Mac, posso fare option click e apre questo file. Quindi possiamo saltare rapidamente ad esso.

Questo nome di process è ora utilizzato dal workflow qui sotto, e possiamo fare la stessa cosa qui. Ci mostra un po' di informazioni su quel process, e di nuovo, posso tenere premuto option, cliccarci sopra, e lo aprirà nell'editor.

Quindi è un modo davvero veloce quando ha molti file per i suoi diversi processi per navigare rapidamente nella sua base di codice in VS Code.

Bene. Questo è sostanzialmente tutto per questo capitolo. Ora facciamo semplicemente la stessa cosa di nuovo per gli altri processi.

## 3. Modularizzare il process convertToUpper()

Quindi creiamo un nuovo file qui. Lo chiamiamo Convert to upper nf. Di nuovo, copiamo lo shebang. E poi tagliamo il process.

Copiamo il nome del process lì, includiamo una nuova dichiarazione include con il nuovo nome del process.

## 4. Modularizzare il process collectGreetings()

E poi facciamo lo stesso per il terzo process. Nuovo file, connect Greetings,

facciamo lo shebang. Tagliamo il process, incolliamo il process, e facciamo una nuova dichiarazione include.

Ora può vedere qui che ho una sottolineatura di errore qui che dice invalid include source. E questo è in realtà un errore genuino che ho fatto perché mi stavo muovendo un po' troppo velocemente. Se guarda attentamente, può vedere che ho perso la T in convert to upper

Quindi VS Code molto utilmente mi ha detto che ho fatto un errore lì. Se correggo quel nome file, l'errore scompare. È un buon esempio del motivo per cui il controllo degli errori all'interno di VS Code è così utile per scrivere codice Nextflow. Altrimenti non l'avrei notato e l'avrei scoperto solo molto più tardi quando avessi provato ad eseguire il workflow.

Il nostro script di pipeline principale ora sembra molto più semplice. Non ha alcun process, abbiamo solo tre dichiarazioni include e il nostro workflow. Non abbiamo cambiato nessuna logica del workflow. Non abbiamo cambiato nessun codice del process, quindi si spera che funzioni esattamente nello stesso modo.

## 4.4. Eseguire il workflow per verificare che faccia la stessa cosa di prima

Verifichiamo. Aprirò un terminale e eseguirò esattamente lo stesso comando di prima.

Certo, ha eseguito i nostri processi, say hello, convert to upper collect greetings, e ci ha dato di nuovo tre saluti.

Quindi abbiamo spostato il nostro codice, ma non abbiamo cambiato nulla su come il workflow viene eseguito ed è completamente invariato. L'unica differenza è che ora abbiamo un codice più pulito, più facile da mantenere, e più facile da condividere con altri.

E questo è tutto. È stato un capitolo breve. È un concetto semplice, ma è molto potente e fondamentale per come scriviamo workflow Nextflow più complessi. Quindi è importante che lo comprenda e prenda l'abitudine di utilizzarlo.

Nel prossimo capitolo, avremo un po' di cambiamento di ritmo e smetteremo di pensare così tanto alla sintassi di scrittura del codice Nextflow, e penseremo un po' a come utilizziamo il software nei processi stessi. Ci raggiunga nella parte cinque per Hello Containers.

[Prossima trascrizione video :octicons-arrow-right-24:](05_hello_containers.md)
