---
title: Nextflow Run
hide:
  - toc
page_type: index_page
index_type: course
additional_information:
  technical_requirements: true
  learning_objectives:
    - Lanciare e gestire l'esecuzione di flussi di lavoro Nextflow
    - Trovare e interpretare gli output (risultati) e i file di log
    - Riconoscere i componenti principali di Nextflow in un semplice flusso di lavoro multi-step
    - Configurare l'esecuzione della pipeline per eseguirla su piattaforme di calcolo comuni, inclusi HPC e cloud
    - Riassumere le best practice per riproducibilità, portabilità e riutilizzo del codice che rendono le pipeline FAIR, inclusa la modularità del codice e i container software
  audience_prerequisites:
    - "**Pubblico:** Questo corso è progettato per chi è completamente nuovo a Nextflow e vuole eseguire pipeline esistenti."
    - "**Competenze:** Si presume una certa familiarità con la riga di comando, concetti di scripting di base e formati di file comuni."
    - "**Dominio:** Gli esercizi sono tutti indipendenti dal dominio, quindi non è richiesta alcuna conoscenza scientifica preliminare."
---

# Nextflow Run

**Nextflow Run è un'introduzione pratica all'esecuzione di flussi di lavoro di analisi dati riproducibili e scalabili.**

Lavorando attraverso esempi pratici ed esercizi guidati, imparerete i fondamenti dell'utilizzo di Nextflow, incluso come eseguire pipeline, gestire file e dipendenze software, parallelizzare l'esecuzione senza sforzo ed eseguire flussi di lavoro in diversi ambienti di calcolo.

Acquisirete le competenze e la fiducia per iniziare a eseguire flussi di lavoro con Nextflow.

<!-- additional_information -->

## Panoramica del corso

### Cosa farete

Questo corso è pratico, con esercizi orientati agli obiettivi strutturati per introdurre le informazioni gradualmente.

Eseguirete diverse versioni di una pipeline Nextflow che elabora input di testo.
Inizierete con una versione semplice che consiste in un singolo passaggio e alla fine progredirete verso una versione multi-step che prende un file CSV di input di testo tabulari, esegue alcuni passaggi di trasformazione e produce un singolo file di testo contenente un'immagine ASCII di un personaggio che dice il testo trasformato.

Questo corso si concentra sull'esecuzione di pipeline (dal nome del comando principale `nextflow run`).
Se cercate un'introduzione allo sviluppo di pipeline Nextflow, consultate [Hello Nextflow](../hello_nextflow/index.md).

### Piano delle lezioni

Abbiamo suddiviso questo in tre parti che si concentreranno ciascuna su aspetti specifici dell'esecuzione e della gestione di pipeline scritte in Nextflow.

| Capitolo del corso                                | Riepilogo                                                                                                                                  | Durata stimata |
| ------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------ | -------------- |
| [Parte 1: Operazioni di base](./01_basics.md)     | Lanciare e gestire l'esecuzione di un semplice flusso di lavoro                                                                            | 30 minuti      |
| [Parte 2: Eseguire pipeline reali](./02_pipeline.md) | Elaborare input complessi, eseguire flussi di lavoro multi-step, utilizzare container e parallelizzare l'esecuzione senza sforzo           | 60 minuti      |
| [Parte 3: Configurazione dell'esecuzione](./03_config.md) | Personalizzare il comportamento della pipeline e ottimizzare l'utilizzo in diversi ambienti computazionali                                 | 60 minuti      |

Alla fine di questo corso, sarete ben preparati per affrontare i prossimi passi nel vostro percorso per eseguire flussi di lavoro riproducibili per le vostre esigenze di calcolo scientifico.

Pronti a seguire il corso?

[Inizia a imparare :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }
