---
title: Home
description: Benvenuti nel portale di formazione della community Nextflow!
hide:
  - toc
  - footer
---

# Nextflow Training

<div class="grid cards" markdown>

-   :material-book-open-variant:{ .lg .middle } __Corsi in autonomia__

    ---

    **Benvenuti nel portale di formazione della community Nextflow!**

    Seguite i corsi qui sotto al vostro ritmo, nel nostro ambiente web o nel vostro.
    Ogni corso è pratico, con esercizi orientati agli obiettivi che potete completare in modo indipendente.

    [Esplora i corsi :material-arrow-down:](#catalog-of-nextflow-training-courses){ .md-button .md-button--primary .mt-1 }

-   :material-account-group-outline:{ .lg .middle } __Eventi di formazione__

    ---

    **Cercate qualcosa oltre ai corsi in autonomia?**

    Trovate eventi di formazione strutturati, indicazioni per organizzare le vostre sessioni di formazione e la nostra licenza open-source e la politica di contribuzione.

    [Vedi gli eventi di formazione :material-arrow-right:](training_events.md){ .md-button .md-button--secondary .mt-1 }

</div>

!!! note "Traduzione assistita da IA"

    Questa traduzione è stata creata utilizzando l'intelligenza artificiale e revisionata da traduttori umani.
    Apprezziamo il vostro feedback e i suggerimenti per miglioramenti.
    Consultate la nostra [guida alla traduzione](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md) per maggiori informazioni.

## Catalog of Nextflow training courses

<div class="grid cards" markdown>

-   :material-account:{ .lg .middle } __Per gli utenti__

    ---

    ### :material-play-circle:{.nextflow-primary} Eseguire pipeline {.mt-1}

    Imparate a eseguire pipeline esistenti senza scrivere codice.

    ??? courses "**Nextflow Run:** Eseguire pipeline con Nextflow"

        Un'introduzione rapida all'esecuzione di pipeline Nextflow che non richiede la comprensione del codice. Copre il lancio delle pipeline, il recupero degli output, l'utilizzo dei container e la configurazione dell'esecuzione a livello base.

        [Vedi la formazione :material-arrow-right:](nextflow_run/index.md){ .md-button .md-button--secondary }

    ??? courses "**Use nf-core:** Trovare ed eseguire pipeline curate dalla community"

        Un'introduzione rapida alla ricerca, all'esecuzione e alla configurazione di pipeline del progetto community nf-core, partendo da una pipeline demo minimale per poi scalare fino a una pipeline di analisi a livello di produzione.

        [Vedi la formazione :material-arrow-right:](nfcore_use/index.md){ .md-button .md-button--secondary }

    ??? courses "**Scale with Seqera:** Lanciare e monitorare pipeline su larga scala"

        Un'introduzione pratica al lancio e al monitoraggio di pipeline Nextflow con Seqera Platform, sia dall'interfaccia web che dalla riga di comando.

        [Vedi la formazione :material-arrow-right:](seqera_scale/index.md){ .md-button .md-button--secondary }

    ---

    ### :material-tune:{.nextflow-primary} Gestire l'esecuzione {.mt-1}

    Imparate a gestire l'esecuzione delle pipeline in modo efficace.

    ??? courses "**Execution Config:** Configurare le pipeline come un professionista"

        Un'introduzione pratica alla configurazione dell'esecuzione di pipeline Nextflow: adattamento a diversi ambienti di calcolo, controllo delle allocazioni di risorse e dei tentativi di riesecuzione, e passaggio tra profili di configurazione predefiniti.

        [Vedi la formazione :material-arrow-right:](execution_config/index.md){ .md-button .md-button--secondary }

    !!! info compact "Altri argomenti in arrivo"

        L'ottimizzazione delle prestazioni, l'esecuzione su HPC/cloud e altro ancora sono previsti per questa sezione.
        Votate su cosa trattare nella nostra [breve sondaggio di interesse](https://seqera.typeform.com/to/JCs91e8v).

-   :material-code-tags:{ .lg .middle } __Per gli sviluppatori__

    ---

    ### :material-wrench:{.nextflow-primary} Scrivere pipeline {.mt-1}

    Imparate a sviluppare le vostre pipeline Nextflow.

    ??? courses "**Hello Nextflow:** Sviluppare le proprie pipeline da zero"

        Questo corso copre i componenti fondamentali del linguaggio Nextflow con un livello di dettaglio sufficiente per sviluppare pipeline semplici ma pienamente funzionali, oltre agli elementi chiave di progettazione, sviluppo e configurazione delle pipeline.

        [Vedi la formazione :material-arrow-right:](hello_nextflow/index.md){ .md-button .md-button--secondary }

    ??? courses "**Build with nf-core:** Usare gli strumenti e le regole di nf-core"

        Per gli sviluppatori Nextflow che desiderano imparare a sviluppare pipeline conformi a [nf-core](https://nf-co.re/).
        Il corso copre la struttura delle pipeline nf-core con un livello di dettaglio sufficiente per sviluppare pipeline semplici ma pienamente funzionali che sfruttano il template nf-core e le best practice di sviluppo, nonché l'utilizzo di moduli nf-core esistenti.

        [Vedi la formazione :material-arrow-right:](nfcore_build/index.md){ .md-button .md-button--secondary }

    ??? catalog "**Side Quests:** Approfondire argomenti avanzati di Nextflow"

        Una raccolta di mini-corsi autonomi pensati per gli sviluppatori Nextflow che desiderano ampliare le proprie competenze e/o approfondire argomenti specifici.
        Sono presentati in modo lineare ma possono essere seguiti in qualsiasi ordine (vedere le dipendenze nella panoramica di ogni mini-corso).

        [Sfoglia i Side Quests :material-arrow-right:](side_quests/index.md){ .md-button .md-button--secondary }

    ---

    ### :material-microscope:{.nextflow-primary} Nextflow per la Scienza {.mt-1}

    Imparate a sviluppare pipeline Nextflow per applicazioni scientifiche specifiche.

    ??? courses "**Genomics:** Sviluppare una pipeline per la variant calling"

        Un corso per ricercatori che desiderano imparare a sviluppare le proprie pipeline di genomica, utilizzando un caso d'uso di variant calling per illustrare i pattern di sviluppo Nextflow essenziali.

        [Vedi la formazione :material-arrow-right:](nf4_science/genomics/index.md){ .md-button .md-button--secondary }

    ??? courses "**RNAseq:** Sviluppare una pipeline per l'elaborazione di RNAseq bulk"

        Un corso per ricercatori che desiderano imparare a sviluppare le proprie pipeline RNAseq, utilizzando un caso d'uso di elaborazione RNAseq bulk per illustrare i pattern di sviluppo Nextflow essenziali.

        [Vedi la formazione :material-arrow-right:](nf4_science/rnaseq/index.md){ .md-button .md-button--secondary }

    ??? courses "**Bioimaging:** Eseguire e configurare pipeline di imaging"

        Un corso per ricercatori che desiderano imparare a eseguire e configurare pipeline di bioimaging, utilizzando nf-core/molkart per illustrare i pattern d'uso Nextflow essenziali.

        [Vedi la formazione :material-arrow-right:](nf4_science/imaging/index.md){ .md-button .md-button--secondary }

</div>

## Setup e Supporto

<div class="grid cards mb-4" markdown>

-   :material-cog-outline:{ .lg .middle } __Ambiente di formazione__

    ---

    Opzioni per configurare il vostro ambiente per le formazioni Nextflow.

    [Vedi gli ambienti di formazione :material-arrow-right:](envsetup/index.md){ .md-button .md-button--secondary }

-   :material-tag-outline:{ .lg .middle } __Versioni di Nextflow__

    ---

    Comprendere e gestire l'evoluzione delle versioni della sintassi di Nextflow.

    [Verifica i requisiti di versione :material-arrow-right:](info/nxf_versions.md){ .md-button .md-button--secondary }

-   :material-file-code-outline:{ .lg .middle } __La pipeline Hello__

    ---

    Riepilogo di cosa fa la pipeline Hello e come è strutturata.

    [Leggi il riepilogo :material-arrow-right:](info/hello_pipeline.md){ .md-button .md-button--secondary }

-   :material-lifebuoy:{ .lg .middle } __Ottenere aiuto__

    ---

    Risorse utili quando avete un problema con la formazione Nextflow.

    [Trova aiuto :material-arrow-right:](help.md){ .md-button .md-button--secondary }

</div>
