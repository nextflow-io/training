---
title: Nextflow Run
hide:
  - toc
page_type: index_page
index_type: course
additional_information:
  technical_requirements: true
  learning_objectives:
    - Avviare e gestire l'esecuzione di workflow Nextflow
    - Trovare e interpretare output (risultati) e file di log
    - Riconoscere i componenti principali di Nextflow in un semplice workflow multi-step
    - Configurare l'esecuzione di pipeline per piattaforme di calcolo comuni inclusi HPC e cloud
    - Riassumere le best practice per riproducibilità, portabilità e riuso del codice che rendono le pipeline FAIR, inclusa la modularità del codice e i container software
  audience_prerequisites:
    - "**Pubblico:** Questo corso è progettato per chi è completamente nuovo a Nextflow e vuole eseguire pipeline esistenti."
    - "**Competenze:** Si assume una certa familiarità con la linea di comando, concetti base di scripting e formati di file comuni."
    - "**Dominio:** Gli esercizi sono tutti indipendenti dal dominio, quindi non è richiesta alcuna conoscenza scientifica pregressa."
---

# Nextflow Run

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduzione assistita da IA - [scopri di più e suggerisci miglioramenti](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

**Nextflow Run è un'introduzione pratica all'esecuzione di workflow di analisi dati riproducibili e scalabili.**

Lavorando attraverso esempi pratici ed esercizi guidati, imparerai i fondamenti dell'utilizzo di Nextflow, incluso come eseguire pipeline, gestire file e dipendenze software, parallelizzare l'esecuzione senza sforzo, e far girare workflow su diversi ambienti di calcolo.

Acquisirai le competenze e la sicurezza per iniziare a eseguire workflow con Nextflow.

<!-- additional_information -->

## Panoramica del corso

### Cosa farai

Questo corso è pratico, con esercizi orientati agli obiettivi strutturati per introdurre le informazioni gradualmente.

Eseguirai diverse versioni di una pipeline Nextflow che elabora input di testo.
Inizierai con una versione semplice che consiste in un singolo step, e progredirai eventualmente verso una versione multi-step che prende un file CSV di input di testo tabulari, esegue alcuni step di trasformazione, e produce un singolo file di testo contenente un'immagine ASCII di un personaggio che pronuncia il testo trasformato.

Questo corso si concentra sull'esecuzione di pipeline (dal nome del comando principale `nextflow run`).
Se cerchi un'introduzione allo sviluppo di pipeline Nextflow, consulta [Hello Nextflow](../hello_nextflow/index.md).

### Piano delle lezioni

Abbiamo suddiviso questo corso in tre parti che si concentreranno ciascuna su aspetti specifici dell'esecuzione e gestione di pipeline scritte in Nextflow.

| Capitolo del corso                                   | Riepilogo                                                                                                              | Durata stimata |
| ---------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------- | -------------- |
| [Parte 1: Operazioni base](./01_basics.md)           | Avvio e gestione dell'esecuzione di un semplice workflow                                                               | 30 min         |
| [Parte 2: Eseguire pipeline reali](./02_pipeline.md) | Elaborazione di input complessi, esecuzione di workflow multi-step, utilizzo di container e parallelizzazione semplice | 60 min         |
| [Parte 3: Configurazione](./03_config.md)            | Personalizzazione del comportamento della pipeline e ottimizzazione dell'utilizzo in diversi ambienti di calcolo       | 60 min         |

Al termine di questo corso, sarai ben preparato per affrontare i prossimi passi nel tuo percorso verso l'esecuzione di workflow riproducibili per le tue esigenze di calcolo scientifico.

Pronto per iniziare il corso?

[Inizia a imparare :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }
