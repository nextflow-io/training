---
title: Home
description: Benvenuti nel portale di formazione della community Nextflow!
hide:
  - toc
  - footer
---

# Formazione Nextflow

<div class="grid cards" markdown>

-   :material-book-open-variant:{ .lg .middle } __Corsi self-service__

    ---

    **Benvenuti nel portale di formazione della community Nextflow!**

    I corsi di formazione elencati di seguito sono progettati per essere utilizzati come risorse self-service.
    È possibile seguirli autonomamente in qualsiasi momento, sia nell'ambiente web che forniamo tramite Github Codespaces, sia nel proprio ambiente locale.

    [Esplora i corsi :material-arrow-right:](#catalogo-dei-corsi-di-formazione-nextflow){ .md-button .md-button--primary .mt-1 }

-   :material-information-outline:{ .lg .middle } __Informazioni aggiuntive__

    ---

    ??? warning "Compatibilità delle versioni"

        <!-- Any update to this content needs to be copied to the local installation page -->
        **A partire da gennaio 2026, tutti i nostri corsi di formazione Nextflow richiedono Nextflow versione 25.10.2 o successiva, con la sintassi strict attivata, salvo diversa indicazione.**

        Per maggiori informazioni sui requisiti di versione e sulla sintassi strict, consultare la [guida alla migrazione nella documentazione Nextflow](https://nextflow.io/docs/latest/strict-syntax.html).

        Le versioni precedenti del materiale formativo corrispondenti alla sintassi precedente sono disponibili tramite il selettore di versione nella barra del menu di questa pagina web.

    ??? terminal "Opzioni per l'ambiente"

        Forniamo un ambiente di formazione basato sul web dove tutto il necessario per seguire la formazione è preinstallato, disponibile tramite Github Codespaces (richiede un account GitHub gratuito).

        [![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

        Se questa opzione non soddisfa le vostre esigenze, consulti le altre [Opzioni per l'ambiente](./envsetup/index.md).

    ??? learning "Eventi formativi"

        Se preferisce seguire la formazione Nextflow come parte di un evento strutturato, ci sono molte opportunità per farlo. vi consigliamo di esplorare le seguenti opzioni:

        - **[Training Weeks]()** organizzate trimestralmente dal team Community
        - **[Seqera Events](https://seqera.io/events/)** includono eventi formativi in presenza organizzati da Seqera (cercare 'Seqera Sessions' e 'Nextflow Summit')
        - **[Nextflow Ambassadors]()** organizzano eventi per le loro community locali
        - **[nf-core events](https://nf-co.re/events)** includono hackathon della community

    ??? people "Informazioni per i formatori"

        Se siete un formatore che organizza le proprie sessioni di formazione, siete liberi di utilizzare i nostri materiali direttamente dal portale di formazione, a condizione di attribuire il corretto credito. Vedere 'Crediti e contributi' di seguito per i dettagli.

        Inoltre, ci piacerebbe ricevere il vostro feedback su come potremmo supportare meglio i vostri sforzi formativi! Vi preghiamo di contattarci all'indirizzo [community@seqera.io](mailto:community@seqera.io) o sul forum della community (vedere la pagina [Aiuto](help.md)).

    ??? licensing "Licenza open-source e politica sui contributi"

        [![Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International (CC BY-NC-SA 4.0)](assets/img/cc_by-nc-sa.svg){ align=right }](https://creativecommons.org/licenses/by-nc-sa/4.0/)

        Questo materiale formativo è sviluppato e mantenuto da [Seqera](https://seqera.io) e rilasciato con una licenza open-source ([CC BY-NC-SA](https://creativecommons.org/licenses/by-nc-sa/4.0/)) a beneficio della community. Se desiderate utilizzare questo materiale in un modo che non rientra nell'ambito della licenza (si noti le limitazioni sull'uso commerciale e la redistribuzione), Vi preghiamo di contattarci all'indirizzo [community@seqera.io](mailto:community@seqera.io) per discutere la vostra richiesta.

        Accogliamo con piacere miglioramenti, correzioni e segnalazioni di bug dalla community. Ogni pagina ha un'icona :material-file-edit-outline: in alto a destra che collega al repository del codice, dove è possibile segnalare problemi o proporre modifiche al materiale sorgente della formazione tramite una pull request. Consultare il `README.md` nel repository per maggiori dettagli.

</div>

!!! note "Traduzione assistita da IA"

    Questa traduzione è stata creata utilizzando l'intelligenza artificiale e revisionata da traduttori umani.
    Apprezziamo il tuo feedback e i tuoi suggerimenti per miglioramenti.
    Consulta la nostra [guida alla traduzione](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md) per maggiori informazioni.

## Catalogo dei corsi di formazione Nextflow

<div class="grid cards" markdown>

-   :material-walk:{ .lg .middle } __Percorso introduttivo__

    ---

    ### :material-compass:{.nextflow-primary} Nextflow per principianti {.mt-1}

    Corsi indipendenti dal dominio applicativo, pensati per chi è completamente nuovo a Nextflow. Ogni corso consiste in una serie di moduli formativi progettati per aiutare gli studenti a sviluppare le proprie competenze progressivamente.

    ??? courses "**Hello Nextflow:** Impara a sviluppare le tue pipeline"

        Questo corso copre i componenti fondamentali del linguaggio Nextflow con sufficiente dettaglio per consentire lo sviluppo di pipeline semplici ma pienamente funzionanti, oltre agli elementi chiave della progettazione, sviluppo e configurazione delle pipeline.

        [Inizia la formazione Hello Nextflow :material-arrow-right:](hello_nextflow/index.md){ .md-button .md-button--secondary }

    ??? courses "**Nextflow Run:** Impara a eseguire pipeline esistenti"

        Un'introduzione concisa all'esecuzione e alla configurazione di pipeline Nextflow, basata sul corso Hello Nextflow per sviluppatori ma con meno enfasi sul codice. Copre esecuzione, output, struttura di base del codice e configurazione per diversi ambienti di calcolo.

        [Inizia la formazione Nextflow Run :material-arrow-right:](nextflow_run/index.md){ .md-button .md-button--secondary }

    ---

    ### :material-microscope:{.nextflow-primary} Nextflow per la Scienza {.mt-1}

    Impara ad applicare i concetti e i componenti presentati in 'Hello Nextflow' a casi d'uso scientifici specifici.

    ??? courses "**Nextflow per la Genomica** (variant calling)"

        Per i ricercatori che desiderano imparare a sviluppare le proprie pipeline genomiche. Il corso utilizza un caso d'uso di variant calling per dimostrare come sviluppare una pipeline genomica semplice ma funzionale.

        [Inizia la formazione Nextflow per la Genomica :material-arrow-right:](nf4_science/genomics/){ .md-button .md-button--secondary }

    ??? courses "**Nextflow per RNAseq** (bulk RNAseq)"

        Per i ricercatori che desiderano imparare a sviluppare le proprie pipeline RNAseq. Il corso utilizza un caso d'uso di elaborazione bulk RNAseq per dimostrare come sviluppare una pipeline RNAseq semplice ma funzionale.

        [Inizia la formazione Nextflow per RNAseq :material-arrow-right:](nf4_science/rnaseq/){ .md-button .md-button--secondary }

    ??? courses "**Nextflow per l'Imaging** (spatial omics)"

        Per i ricercatori in imaging e spatial omics che desiderano imparare a eseguire e personalizzare pipeline di analisi. Il corso utilizza la pipeline nf-core/molkart per fornire una pipeline biologicamente rilevante e dimostrare come eseguire, configurare e gestire gli input per i workflow Nextflow.

        [Inizia la formazione Nextflow per l'Imaging :material-arrow-right:](nf4_science/imaging/){ .md-button .md-button--secondary }

-   :material-run:{ .lg .middle } __Percorso avanzato__

    ---

    ### :material-bridge:{.nextflow-primary} Da Nextflow a nf-core {.mt-1}

    Impara a utilizzare il codice e le best practice del progetto community [nf-core](https://nf-co.re/).

    Questi corsi aiutano a passare dai fondamenti di Nextflow alle best practice di nf-core.
    Comprenda come e perché la community nf-core costruisce le pipeline, e come contribuire e riutilizzare queste tecniche.

    ??? courses "**Hello nf-core:** Inizia con nf-core"

        Per gli sviluppatori che desiderano imparare a eseguire e sviluppare pipeline conformi a [nf-core](https://nf-co.re/). Il corso copre la struttura delle pipeline nf-core con sufficiente dettaglio per consentire lo sviluppo di pipeline semplici ma pienamente funzionali che seguono il template nf-core e le best practice di sviluppo, oltre all'utilizzo dei moduli nf-core esistenti.

        [Inizia la formazione Hello nf-core :material-arrow-right:](hello_nf-core/index.md){ .md-button .md-button--secondary }

    ---

    ### :material-rocket-launch:{.nextflow-primary} Formazione Nextflow Avanzata {.mt-1}

    Impara concetti e meccanismi avanzati per sviluppare e distribuire pipeline Nextflow per affrontare casi d'uso reali.

    ??? courses "**Side Quests:** Approfondimenti su argomenti specifici"

        Mini-corsi indipendenti pensati per sviluppatori Nextflow che desiderano ampliare il proprio raggio d'azione e/o approfondire le proprie competenze su argomenti particolari. Sono presentati in modo lineare ma possono essere seguiti in qualsiasi ordine (vedere le dipendenze nella panoramica di ogni mini-corso).

        [Esplora le Side Quests :material-arrow-right:](side_quests/){ .md-button .md-button--secondary }

    ??? courses "**Training Collections:** Percorsi di apprendimento consigliati attraverso le Side Quests"

        Le Training Collections combinano più Side Quests per fornire un'esperienza di apprendimento completa su un tema o caso d'uso particolare.

        [Esplora le Training Collections :material-arrow-right:](training_collections/){ .md-button .md-button--secondary }

</div>

!!! info "Cerchi materiali formativi archiviati?"

    I materiali formativi più vecchi (Fundamentals Training, Advanced Training e altri corsi sperimentali) sono stati rimossi dal portale di formazione in quanto incompatibili con la sintassi strict di Nextflow 3.0.
    Se hai bisogno di accedere a questi materiali, sono disponibili nella [cronologia git](https://github.com/nextflow-io/training) prima di gennaio 2026.

---

<div markdown class="homepage_logos">

![Seqera](assets/img/seqera_logo.png#only-light)

![Seqera](assets/img/seqera_logo_dark.png#only-dark)

</div>
