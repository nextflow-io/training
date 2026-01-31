# Installazione manuale

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduzione assistita da IA - [scopri di più e suggerisci miglioramenti](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

È possibile installare tutto il necessario per eseguire la formazione nel proprio ambiente locale manualmente.

Qui abbiamo documentato come farlo su sistemi standard compatibili con POSIX (assumendo un computer personale come un laptop).
Tenga presente che alcuni dettagli potrebbero essere diversi a seconda del vostro sistema specifico.

!!! tip "Suggerimento"

    Prima di procedere, ha considerato l'utilizzo dell'[approccio Devcontainers](03_devcontainer.md)?
    Fornisce tutti gli strumenti e le dipendenze necessarie senza richiedere l'installazione manuale.

## Requisiti software generali

Nextflow può essere utilizzato su qualsiasi sistema compatibile con POSIX (Linux, macOS, Windows Subsystem for Linux, ecc.) con Java installato.
I nostri corsi di formazione hanno alcuni requisiti aggiuntivi.

In totale, dovrà avere installato il seguente software:

- Bash o shell equivalente
- [Java 11 (o successivo, fino a 21)](https://www.oracle.com/technetwork/java/javase/downloads/index.html)
- [Git](https://git-scm.com/)
- [Docker](https://docs.docker.com/get-docker/)
- [Conda](https://conda.io/) 4.5 (o successivo)
- [VSCode](https://code.visualstudio.com) con l'[estensione Nextflow](https://www.nextflow.io/docs/latest/developer-env.html#devenv-nextflow)

L'applicazione VSCode è tecnicamente opzionale ma raccomandiamo fortemente di utilizzarla per lavorare attraverso i corsi così come per il vostro lavoro di sviluppo Nextflow in generale.

Il manuale della documentazione Nextflow fornisce istruzioni per l'installazione di queste dipendenze in [Configurazione dell'ambiente](https://www.nextflow.io/docs/latest/developer-env.html).

## Nextflow e strumenti nf-core

Dovrà installare Nextflow stesso, oltre agli strumenti nf-core, come descritto negli articoli collegati di seguito:

- [Installazione Nextflow](https://www.nextflow.io/docs/latest/install.html)
- [Strumenti nf-core](https://nf-co.re/docs/nf-core-tools/installation)

Raccomandiamo di utilizzare l'opzione di auto-installazione per Nextflow e l'opzione PyPI per gli strumenti nf-core.

!!! warning "Compatibilità delle versioni"

    <!-- Any update to this content needs to be copied to the home page -->
    **A partire da gennaio 2026, tutti i nostri corsi di formazione Nextflow richiedono Nextflow versione 25.10.2 o successiva, con la sintassi strict v2 attivata, salvo diversa indicazione.**

    Per maggiori informazioni sui requisiti di versione e sulla sintassi strict v2, consultare la guida [Versioni di Nextflow](../info/nxf_versions.md).

    Le versioni precedenti del materiale formativo corrispondenti alla sintassi precedente sono disponibili tramite il selettore di versione nella barra del menu di questa pagina web.

## Materiali formativi

Il modo più semplice per scaricare i materiali formativi è clonare l'intero repository utilizzando questo comando:

```bash
git clone https://github.com/nextflow-io/training.git
```

Ogni corso ha la propria directory.
Per lavorare attraverso un corso, apra una finestra del terminale (idealmente, dall'interno dell'applicazione VSCode) e faccia `cd` nella directory pertinente.

Può quindi seguire le istruzioni del corso fornite sul sito web.
