---
title: Hello Nextflow
hide:
  - toc
page_type: index_page
index_type: course
additional_information:
  technical_requirements: true
  learning_objectives:
    - Avviare e gestire l'esecuzione di workflow Nextflow
    - Trovare e interpretare gli output (risultati) e i file di log generati da Nextflow
    - Risolvere problemi di base
    - Costruire un workflow semplice multi-step dai componenti fondamentali di Nextflow
    - Distinguere tra tipi essenziali di channel factory e operatori e utilizzarli efficacemente in un workflow semplice
    - Configurare l'esecuzione della pipeline per funzionare su piattaforme di calcolo comuni inclusi HPC e cloud
    - Applicare le best practice per riproducibilità, portabilità e riutilizzo del codice che rendono le pipeline FAIR, inclusa la modularità del codice e i container software
  audience_prerequisites:
    - "**Pubblico:** Questo corso è progettato per chi è completamente nuovo a Nextflow e desidera sviluppare le proprie pipeline."
    - "**Competenze:** Si presume una certa familiarità con la riga di comando, concetti di base di scripting e formati di file comuni."
    - "**Dominio:** Gli esercizi sono tutti indipendenti dal dominio applicativo, quindi non è richiesta alcuna conoscenza scientifica preliminare."
  videos_playlist: https://youtube.com/playlist?list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&si=eF7cLR62goy-lc6n
---

# Hello Nextflow

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduzione assistita da IA - [scopri di più e suggerisci miglioramenti](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

**Hello Nextflow è un'introduzione pratica alla costruzione di workflow di analisi dati riproducibili e scalabili.**

Attraverso esempi pratici ed esercizi guidati, imparerete i fondamenti dello sviluppo di pipeline con Nextflow, incluso come definire processi, connetterli in pipeline, gestire file e dipendenze software, parallelizzare l'esecuzione senza sforzo ed eseguire workflow in diversi ambienti di calcolo.

Acquisirete le competenze e la sicurezza per iniziare a sviluppare e eseguire i vostri workflow con Nextflow.

<!-- additional_information -->

## Panoramica del corso

Questo corso è progettato per essere pratico, con esercizi orientati agli obiettivi strutturati per introdurre le informazioni gradualmente.

Svilupperete una semplice pipeline Nextflow che prende alcuni input di testo, esegue alcuni passaggi di trasformazione e produce un singolo file di testo contenente un'immagine ASCII di un personaggio che dice il testo trasformato.

### Piano delle lezioni

Per evitare di sovraccaricarvi con concetti e codice, abbiamo suddiviso questo in sei parti che si concentreranno ciascuna su aspetti specifici dello sviluppo di pipeline con Nextflow.

| Capitolo del corso                                    | Riepilogo                                                                                                                                        | Durata stimata |
| ----------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------ | -------------- |
| [Parte 1: Hello World](./01_hello_world.md)           | Componenti e principi di base coinvolti nell'assemblaggio e nell'esecuzione di un workflow Nextflow                                              | 30 min         |
| [Parte 2: Hello Channels](./02_hello_channels.md)     | Utilizzo di canali e operatori per elaborare input e parallelizzare l'esecuzione senza sforzo                                                    | 45 min         |
| [Parte 3: Hello Workflow](./03_hello_workflow.md)     | Utilizzo dei canali per concatenare più step e gestire il trasferimento di dati tra gli step                                                     | 60 min         |
| [Parte 4: Hello Modules](./04_hello_modules.md)       | Applicazione dei principi di modularità del codice per aumentare la riusabilità e ridurre l'onere di manutenzione                                | 20 min         |
| [Parte 5: Hello Containers](./05_hello_containers.md) | Utilizzo dei container come meccanismo per gestire le dipendenze software e aumentare la riproducibilità                                         | 60 min         |
| [Parte 6: Hello Config](./06_hello_config.md)         | Personalizzazione del comportamento della pipeline e ottimizzazione dell'utilizzo in diversi ambienti computazionali | 60 min         |

Al termine di questo corso, sarete ben preparati per affrontare i prossimi passi nel vostro percorso per sviluppare workflow riproducibili per le vostre esigenze di calcolo scientifico.

Pronti a seguire il corso?

[Inizia :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }

<!-- Clearfix for float -->
<div style="content: ''; clear: both; display: table;"></div>
