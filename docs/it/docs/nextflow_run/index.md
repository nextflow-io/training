---
title: Nextflow Run
hide:
  - toc
page_type: index_page
index_type: course
additional_information:
  technical_requirements: true
  learning_objectives:
    - Avviare e gestire pipeline Nextflow dalla linea di comando
    - Capire come i canali e gli operatori abilitano flussi di lavoro multi-input e multi-step efficienti
    - Usare i container per gestire le dipendenze software e garantire la riproducibilità
    - Configurare l'esecuzione della pipeline e gli output
    - Generare report di esecuzione, ispezionare la cronologia delle esecuzioni passate e pulire le vecchie directory di lavoro
    - Eseguire pipeline direttamente da repository remoti come GitHub
  audience_prerequisites:
    - "**Pubblico:** Questo corso è progettato per chi è completamente nuovo a Nextflow e vuole eseguire pipeline esistenti."
    - "**Competenze:** Si assume una certa familiarità con la linea di comando, concetti base di scripting e formati di file comuni."
    - "**Dominio:** Gli esercizi sono tutti indipendenti dal dominio, quindi non è richiesta alcuna conoscenza scientifica pregressa."
---

# Nextflow Run

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduzione assistita da IA - [scopri di più e suggerisci miglioramenti](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

**Nextflow Run è un'introduzione pratica all'esecuzione di workflow di analisi dati riproducibili e scalabili.**

Lavorando attraverso una serie di esercizi orientati agli obiettivi, imparerai i fondamenti dell'avvio e della gestione di pipeline Nextflow, capirai come i canali e gli operatori abilitano l'elaborazione parallela di input multipli, e utilizzerai i container per gestire le dipendenze software.

Acquisirai le competenze e la sicurezza per iniziare a eseguire workflow con Nextflow.

<!-- additional_information -->

## Panoramica del corso

Questo corso è pratico, con esercizi orientati agli obiettivi strutturati per introdurre le informazioni gradualmente.

Eseguirai diverse versioni di una pipeline Nextflow che elabora input di testo, partendo da una semplice versione a singolo step e progredendo verso una versione multi-step che prende un file CSV di input, esegue alcuni step di trasformazione, e produce un singolo file di testo contenente ASCII art generata da uno strumento containerizzato.

Questo corso si concentra sull'esecuzione di pipeline (dal nome del comando principale `nextflow run`).
Se cerchi un'introduzione allo sviluppo di pipeline Nextflow, consulta [Hello Nextflow](../hello_nextflow/index.md).

!!! note "Nota"

    Stai cercando la versione precedente di questo corso? È stata sostituita dalla versione in questa pagina, ma è ancora consultabile nella [release 3.6.1](https://training.nextflow.io/3.6.1/nextflow_run/) del sito di formazione.

### Piano delle lezioni

| Capitolo del corso                                                       | Riepilogo                                                                                                                 | Durata stimata |
| ------------------------------------------------------------------------ | ------------------------------------------------------------------------------------------------------------------------- | -------------- |
| [Parte 1: Eseguire Nextflow](./01_run_nextflow.md)                       | Avviare e gestire pipeline Nextflow, e comprendere i meccanismi essenziali del flusso di lavoro                           | 25 min         |
| [Parte 2: Configurare la pipeline](./02_configure_pipeline.md)           | Configurare l'esecuzione della pipeline e gli output usando `nextflow.config`                                             | 20 min         |
| [Parte 3: Gestire le esecuzioni del workflow](./03_manage_executions.md) | Generare report di esecuzione, ispezionare la cronologia delle esecuzioni passate e pulire le vecchie directory di lavoro | 10 min         |
| [Parte 4: Eseguire pipeline remote](./04_remote_repositories.md)         | Eseguire una pipeline direttamente da GitHub e bloccarla a una revisione specifica                                        | 10 min         |

Al termine di questo corso, sarai ben preparato per affrontare i prossimi passi nel tuo percorso verso l'esecuzione di workflow riproducibili per le tue esigenze di calcolo scientifico.

Pronto per iniziare il corso?

[Inizia a imparare :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }
