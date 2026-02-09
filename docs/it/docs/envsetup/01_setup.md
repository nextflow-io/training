# GitHub Codespaces

GitHub Codespaces è una piattaforma basata sul web che ci permette di fornire un ambiente preconfigurato per la formazione, supportato da macchine virtuali nel cloud.
La piattaforma è gestita da Github (di proprietà di Microsoft) ed è accessibile gratuitamente (con quote di utilizzo) a chiunque abbia un account Github.

!!! warning "Attenzione"

    Gli account associati a organizzazioni potrebbero essere soggetti a determinate restrizioni aggiuntive.
    Se questo è il vostro caso, potreste dover utilizzare un account personale indipendente, oppure utilizzare un'installazione locale.

## Creare un account GitHub

Potete creare un account GitHub gratuito dalla [home page di GitHub](https://github.com/).

## Avviare il vostro GitHub Codespace

Una volta effettuato l'accesso a GitHub, aprite questo link nel vostro browser per aprire l'ambiente di formazione Nextflow: <https://codespaces.new/nextflow-io/training?quickstart=1&ref=master>

In alternativa, potete cliccare sul pulsante mostrato qui sotto, che viene ripetuto in ogni corso di formazione (tipicamente nella pagina di Orientamento).

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

Dovreste visualizzare una pagina dove potete creare un nuovo GitHub Codespace:

![Create a GitHub Codespace](img/codespaces_create.png)

### Configurazione

Per l'uso generale, non dovreste aver bisogno di configurare nulla.
A meno che non sia specificato diversamente nel corso che state iniziando, potete semplicemente cliccare sul pulsante principale per continuare.

Tuttavia, è possibile personalizzare l'ambiente cliccando sul pulsante "Change options".

??? info "Opzioni di configurazione"

    Se cliccate sul pulsante "Change options", avrete la possibilità di personalizzare quanto segue:

    #### Branch

    Questo vi permette di selezionare una versione diversa dei materiali di formazione.
    Il branch `master` generalmente contiene correzioni di bug e materiali che sono stati sviluppati e approvati di recente ma non sono ancora stati pubblicati sul sito web.
    Altri branch contengono lavori in corso che potrebbero non essere completamente funzionali.

    #### Machine type

    Questo vi permette di personalizzare la macchina virtuale che utilizzerete per seguire la formazione.

    Utilizzare una macchina con più core vi permette di sfruttare maggiormente la capacità di Nextflow di parallelizzare l'esecuzione del flusso di lavoro.
    Tuttavia, consumerà la vostra quota gratuita più velocemente, quindi non raccomandiamo di modificare questa impostazione a meno che non sia consigliato nelle istruzioni del corso che state pianificando di seguire.

    Vedere 'Quote di GitHub Codespaces' qui sotto per maggiori dettagli sulle quote.

### Tempo di avvio

Aprire un nuovo ambiente GitHub Codespaces per la prima volta può richiedere diversi minuti, perché il sistema deve configurare la vostra macchina virtuale, quindi non preoccupatevi se c'è un tempo di attesa.
Tuttavia, non dovrebbe richiedere più di cinque minuti.

## Navigare nell'interfaccia di formazione

Una volta caricato il vostro GitHub Codespaces, dovreste vedere qualcosa di simile a quanto segue (che potrebbe aprirsi in modalità chiara a seconda delle preferenze del vostro account):

![GitHub Codespaces welcome](img/codespaces_welcome.png)

Questa è l'interfaccia dell'IDE VSCode, una popolare applicazione di sviluppo di codice che raccomandiamo di utilizzare per lo sviluppo con Nextflow.

- **L'editor principale** è dove si apriranno il codice Nextflow e altri file di testo. Qui è dove modificherete il codice. Quando aprite il codespace, questo vi mostrerà un'anteprima del file `README.md`.
- **Il terminale** sotto l'editor principale vi permette di eseguire comandi. Qui è dove eseguirete tutte le righe di comando fornite nelle istruzioni del corso.
- **La barra laterale** vi permette di personalizzare il vostro ambiente ed eseguire attività di base (copiare, incollare, aprire file, cercare, git, ecc.). Per impostazione predefinita è aperta sull'esploratore di file, che vi permette di sfogliare i contenuti del repository. Cliccando su un file nell'esploratore lo aprirete nella finestra dell'editor principale.

Potete regolare le proporzioni relative dei pannelli della finestra come preferite.

<!-- TODO (future) Link to development best practices side quest? -->

## Altre note sull'utilizzo di GitHub Codespaces

### Riprendere una sessione

Una volta creato un ambiente, potete facilmente riprenderlo o riavviarlo e continuare da dove avevate lasciato.
Il vostro ambiente andrà in timeout dopo 30 minuti di inattività e salverà le vostre modifiche per un massimo di 2 settimane.

Potete riaprire un ambiente da <https://github.com/codespaces/>.
Gli ambienti precedenti saranno elencati.
Cliccate su una sessione per riprenderla.

![List GitHub Codespace sessions](img/codespaces_list.png)

Se avete salvato l'URL del vostro precedente ambiente GitHub Codespaces, potete semplicemente aprirlo nel vostro browser.
In alternativa, cliccate sullo stesso pulsante che avete utilizzato per crearlo in primo luogo:

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

Dovreste vedere la sessione precedente, l'opzione predefinita è riprenderla:

![Resume a GitHub Codespace](img/codespaces_resume.png)

### Salvare file sulla vostra macchina locale

Per salvare qualsiasi file dal pannello dell'esploratore, fate clic destro sul file e selezionate `Download`.

### Gestire le quote di GitHub Codespaces

GitHub Codespaces vi offre fino a 15 GB-mese di storage al mese, e 120 ore-core al mese.
Questo equivale a circa 60 ore di runtime dell'ambiente predefinito utilizzando lo spazio di lavoro standard (2 core, 8 GB RAM e 32 GB di storage).

Potete crearli con più risorse (vedere la spiegazione sopra), ma questo consumerà il vostro utilizzo gratuito più velocemente e avrete meno ore di accesso a questo spazio.
Per esempio, se selezionate una macchina a 4 core invece dei 2 core predefiniti, la vostra quota si esaurirà nella metà del tempo.

Opzionalmente, potete acquistare l'accesso a più risorse.

Per maggiori informazioni, consultate la documentazione di GitHub:
[About billing for GitHub Codespaces](https://docs.github.com/en/billing/managing-billing-for-your-products/managing-billing-for-github-codespaces/about-billing-for-github-codespaces)
