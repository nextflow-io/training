# Installazione manuale

È possibile installare tutto il necessario per seguire la formazione nel proprio ambiente locale manualmente.

Qui abbiamo documentato come farlo su sistemi standard compatibili con POSIX (assumendo una macchina personale come un laptop).
Tenete presente che alcuni dettagli potrebbero essere diversi a seconda del vostro sistema specifico.

!!! tip

    Prima di procedere, avete considerato l'approccio [Devcontainers](03_devcontainer.md)?
    Fornisce tutti gli strumenti e le dipendenze necessarie senza richiedere un'installazione manuale.

## Requisiti software generali

Nextflow può essere utilizzato su qualsiasi sistema compatibile con POSIX (Linux, macOS, Windows Subsystem for Linux, ecc.) con Java installato.
I nostri corsi di formazione hanno alcuni requisiti aggiuntivi.

In totale, dovrete avere installato il software seguente:

- Bash o shell equivalente
- [Java 11 (o successivo, fino alla 21)](https://www.oracle.com/technetwork/java/javase/downloads/index.html)
- [Git](https://git-scm.com/)
- [Docker](https://docs.docker.com/get-docker/)
- [Conda](https://conda.io/) 4.5 (o successivo)
- [VSCode](https://code.visualstudio.com) con l'[estensione Nextflow](https://www.nextflow.io/docs/latest/developer-env.html#devenv-nextflow)

L'applicazione VSCode è tecnicamente opzionale ma raccomandiamo vivamente di utilizzarla per seguire i corsi così come per il vostro lavoro di sviluppo con Nextflow in generale.

Il manuale della documentazione di Nextflow fornisce istruzioni per installare queste dipendenze nella sezione [Environment setup](https://www.nextflow.io/docs/latest/developer-env.html).

## Nextflow e strumenti nf-core

Dovrete installare Nextflow stesso, più gli strumenti nf-core, come dettagliato negli articoli linkati di seguito:

- [Installazione di Nextflow](https://www.nextflow.io/docs/latest/install.html)
- [Strumenti nf-core](https://nf-co.re/docs/nf-core-tools/installation)

Raccomandiamo di utilizzare l'opzione di auto-installazione per Nextflow e l'opzione PyPI per gli strumenti nf-core.

!!! warning "Compatibilità delle versioni"

    <!-- Any update to this content needs to be copied to the home page -->
    **A partire da gennaio 2026, tutti i nostri corsi di formazione su Nextflow richiedono Nextflow versione 25.10.2 o successiva, con la sintassi strict v2 attivata, salvo diversa indicazione.**

    Per maggiori informazioni sui requisiti di versione e sulla sintassi strict v2, consultate la guida [Versioni di Nextflow](../info/nxf_versions.md).

    Le versioni precedenti del materiale di formazione corrispondenti alla sintassi precedente sono disponibili tramite il selettore di versione nella barra del menu di questa pagina web.

## Materiali di formazione

Il modo più semplice per scaricare i materiali di formazione è clonare l'intero repository utilizzando questo comando:

```bash
git clone https://github.com/nextflow-io/training.git
```

Ogni corso ha la propria directory.
Per seguire un corso, aprite una finestra del terminale (idealmente, dall'interno dell'applicazione VSCode) ed eseguite `cd` nella directory pertinente.

Potete quindi seguire le istruzioni del corso fornite sul sito web.
