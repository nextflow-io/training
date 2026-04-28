---
title: Home
description: Benvenuti nel portale di formazione della community Nextflow!
hide:
  - toc
  - footer
---

# Nextflow Training

<div class="grid cards" markdown>

-   :material-book-open-variant:{ .lg .middle } __Corsi self-service__

    ---

    **Benvenuti nel portale di formazione della community Nextflow!**

    I corsi di formazione elencati di seguito sono progettati per essere utilizzati come risorsa self-service.
    Potete seguirli autonomamente in qualsiasi momento, sia nell'ambiente web che mettiamo a disposizione tramite Github Codespaces, sia nel vostro ambiente locale.

    [Esplora i corsi :material-arrow-right:](#catalog-of-nextflow-training-courses){ .md-button .md-button--primary .mt-1 }

-   :material-information-outline:{ .lg .middle } __Informazioni aggiuntive__

    ---

    ??? warning "Compatibilità delle versioni"

        <!-- Qualsiasi aggiornamento a questo contenuto deve essere copiato nella pagina di installazione locale -->
        **A partire da gennaio 2026, tutti i nostri corsi di formazione Nextflow richiedono Nextflow versione 25.10.2 o successiva, con la sintassi strict attivata, salvo diversa indicazione.**

        Per ulteriori informazioni sui requisiti di versione e sulla sintassi strict, consultate la [guida alla migrazione della documentazione Nextflow](https://nextflow.io/docs/latest/strict-syntax.html).

        Le versioni precedenti del materiale di formazione corrispondenti alla sintassi precedente sono disponibili tramite il selettore di versione nella barra dei menu di questa pagina web.

    ??? terminal "Opzioni di ambiente"

        Mettiamo a disposizione un ambiente di formazione web-based in cui tutto il necessario per seguire la formazione è preinstallato, disponibile tramite Github Codespaces (richiede un account GitHub gratuito).

        [![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

        Se questa opzione non soddisfa le vostre esigenze, consultate le altre [opzioni di ambiente](./envsetup/index.md).

    ??? learning "Eventi di formazione"

        Se preferite seguire la formazione Nextflow nell'ambito di un evento strutturato, ci sono molte opportunità per farlo. Vi consigliamo di consultare le seguenti opzioni:

        - **[Training Weeks]()** organizzate trimestralmente dal team della Community
        - **[Seqera Events](https://seqera.io/events/)** includono eventi di formazione in presenza organizzati da Seqera (cercate 'Seqera Sessions' e 'Nextflow Summit')
        - **[Nextflow Ambassadors]()** organizzano eventi per la loro community locale
        - **[nf-core events](https://nf-co.re/events)** includono hackathon della community

    ??? people "Informazioni per i formatori"

        Se siete un istruttore che organizza le proprie sessioni di formazione, siete i benvenuti a utilizzare i nostri materiali direttamente dal portale di formazione, a condizione di attribuire il giusto credito. Consultate la sezione 'Crediti e contributi' qui sotto per i dettagli.

        Inoltre, ci farebbe molto piacere sentire la vostra opinione su come potremmo supportare meglio le vostre attività di formazione! Contattateci all'indirizzo [community@seqera.io](mailto:community@seqera.io) o sul forum della community (consultate la pagina [Aiuto](help.md)).

    ??? licensing "Licenza open-source e politica di contribuzione"

        [![Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International (CC BY-NC-SA 4.0)](assets/img/cc_by-nc-sa.svg){ align=right }](https://creativecommons.org/licenses/by-nc-sa/4.0/)

        Questo materiale di formazione è sviluppato e mantenuto da [Seqera](https://seqera.io) e rilasciato con una licenza open-source ([CC BY-NC-SA](https://creativecommons.org/licenses/by-nc-sa/4.0/)) a beneficio della community. Se desiderate utilizzare questo materiale in un modo che esuli dall'ambito della licenza (tenete presente le limitazioni sull'uso commerciale e la ridistribuzione), contattateci all'indirizzo [community@seqera.io](mailto:community@seqera.io) per discutere la vostra richiesta.

        Accogliamo con piacere miglioramenti, correzioni e segnalazioni di bug dalla community. Ogni pagina ha un'icona :material-file-edit-outline: in alto a destra che rimanda al repository del codice, dove potete segnalare problemi o proporre modifiche al materiale di formazione tramite una pull request. Consultate il file `README.md` nel repository per ulteriori dettagli.

</div>

!!! note "Traduzione assistita da IA"

    Questa traduzione è stata creata utilizzando l'intelligenza artificiale e revisionata da traduttori umani.
    Apprezziamo il vostro feedback e i suggerimenti per miglioramenti.
    Consultate la nostra [guida alla traduzione](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md) per maggiori informazioni.

## Catalog of Nextflow training courses

<div class="grid cards" markdown>

-   :material-walk:{ .lg .middle } __Percorso introduttivo__

    ---

    ### :material-compass:{.nextflow-primary} Nextflow per i Principianti {.mt-1}

    Corsi indipendenti dal dominio, pensati per chi è completamente nuovo a Nextflow. Ogni corso è composto da una serie di moduli di formazione progettati per aiutare i partecipanti a sviluppare progressivamente le proprie competenze.

    ??? courses "**Hello Nextflow:** Impara a sviluppare le tue pipeline"

        Questo corso copre i componenti fondamentali del linguaggio Nextflow con un livello di dettaglio sufficiente per sviluppare pipeline semplici ma pienamente funzionali, oltre agli elementi chiave di progettazione, sviluppo e configurazione delle pipeline.

        [Inizia la formazione Hello Nextflow :material-arrow-right:](hello_nextflow/index.md){ .md-button .md-button--secondary }

    ??? courses "**Nextflow Run:** Impara a eseguire pipeline esistenti"

        Un'introduzione concisa all'esecuzione e alla configurazione di pipeline Nextflow, basata sul corso per sviluppatori Hello Nextflow ma con meno enfasi sul codice. Copre l'esecuzione, gli output, la struttura di base del codice e la configurazione per diversi ambienti di calcolo.

        [Inizia la formazione Nextflow Run :material-arrow-right:](nextflow_run/index.md){ .md-button .md-button--secondary }

    ---

    ### :material-microscope:{.nextflow-primary} Nextflow per la Scienza {.mt-1}

    Impara ad applicare i concetti e i componenti presentati in 'Hello Nextflow' a casi d'uso scientifici specifici.

    ??? courses "**Nextflow for Genomics** (variant calling)"

        Per i ricercatori che desiderano imparare a sviluppare le proprie pipeline di genomica. Il corso utilizza un caso d'uso di variant calling per dimostrare come sviluppare una pipeline di genomica semplice ma funzionale.

        [Inizia la formazione Nextflow for Genomics :material-arrow-right:](nf4_science/genomics/index.md){ .md-button .md-button--secondary }

    ??? courses "**Nextflow for RNAseq** (bulk RNAseq)"

        Per i ricercatori che desiderano imparare a sviluppare le proprie pipeline RNAseq. Il corso utilizza un caso d'uso di elaborazione bulk RNAseq per dimostrare come sviluppare una pipeline RNAseq semplice ma funzionale.

        [Inizia la formazione Nextflow for RNAseq :material-arrow-right:](nf4_science/rnaseq/index.md){ .md-button .md-button--secondary }

    ??? courses "**Nextflow for Imaging** (spatial omics)"

        Per i ricercatori nel campo dell'imaging e della spatial omics che desiderano imparare a eseguire e personalizzare pipeline di analisi. Il corso utilizza la pipeline nf-core/molkart per fornire una pipeline biologicamente rilevante e dimostrare come eseguire, configurare e gestire gli input per i flussi di lavoro Nextflow.

        [Inizia la formazione Nextflow for Imaging :material-arrow-right:](nf4_science/imaging/index.md){ .md-button .md-button--secondary }

-   :material-run:{ .lg .middle } __Percorso avanzato__

    ---

    ### :material-bridge:{.nextflow-primary} Da Nextflow a nf-core {.mt-1}

    Impara a utilizzare il codice e le best practice del progetto community [nf-core](https://nf-co.re/).

    Questi corsi vi aiutano a passare dai fondamenti di Nextflow alle best practice di nf-core.
    Capite come e perché la community nf-core sviluppa le pipeline, e come potete contribuire e riutilizzare queste tecniche.

    ??? courses "**Hello nf-core:** Inizia con nf-core"

        Per gli sviluppatori che desiderano imparare a eseguire e sviluppare pipeline conformi a [nf-core](https://nf-co.re/). Il corso copre la struttura delle pipeline nf-core con un livello di dettaglio sufficiente per sviluppare pipeline semplici ma pienamente funzionali che seguono il template nf-core e le best practice di sviluppo, nonché per utilizzare i moduli nf-core esistenti.

        [Inizia la formazione Hello nf-core :material-arrow-right:](hello_nf-core/index.md){ .md-button .md-button--secondary }

    ---

    ### :material-rocket-launch:{.nextflow-primary} Formazione Avanzata su Nextflow {.mt-1}

    Apprendi concetti e meccanismi avanzati per sviluppare e distribuire pipeline Nextflow per affrontare casi d'uso del mondo reale.

    ??? courses "**Side Quests:** Approfondimenti su argomenti specifici"

        Mini-corsi autonomi pensati per gli sviluppatori Nextflow che desiderano ampliare il proprio repertorio e/o approfondire le proprie competenze su argomenti particolari. Sono presentati in sequenza ma possono essere seguiti in qualsiasi ordine (consultate le dipendenze nella panoramica di ogni mini-corso).

        [Esplora i Side Quests :material-arrow-right:](side_quests/index.md){ .md-button .md-button--secondary }

    ??? courses "**Training Collections:** Percorsi di apprendimento consigliati attraverso i Side Quests"

        Le Training Collections combinano più Side Quests per offrire un'esperienza di apprendimento completa attorno a un tema o caso d'uso specifico.

        [Esplora le Training Collections :material-arrow-right:](training_collections/index.md){ .md-button .md-button--secondary }

</div>

!!! info "Cerchi materiali di formazione archiviati?"

    I materiali di formazione precedenti (Fundamentals Training, Advanced Training e altri corsi sperimentali) sono stati rimossi dal portale di formazione in quanto incompatibili con la sintassi strict di Nextflow 3.0.
    Se avete bisogno di accedere a questi materiali, sono disponibili nella [cronologia git](https://github.com/nextflow-io/training) precedente a gennaio 2026.

---

<div markdown class="homepage_logos">

![Seqera](assets/img/seqera_logo.png#only-light)

![Seqera](assets/img/seqera_logo_dark.png#only-dark)

</div>
