---
title: Execution Config
hide:
  - toc
page_type: index_page
index_type: course
additional_information:
  technical_requirements: true
  learning_objectives:
    - Passare da una tecnologia di packaging software a un'altra tra Docker e Conda
    - Selezionare una piattaforma di esecuzione e capire come Nextflow adatta l'esecuzione delle attività
    - Controllare l'allocazione delle risorse di calcolo e riprovare automaticamente le attività che falliscono
    - Definire e combinare profili per passare tra configurazioni predefinite
  audience_prerequisites:
    - "**Pubblico:** Questo corso è pensato per chi sa già come lanciare pipeline Nextflow in locale e vuole configurare l'esecuzione in modo più approfondito."
    - "**Competenze:** È richiesta una certa familiarità con la riga di comando."
    - "**Corsi:** È necessario aver completato [Nextflow Run](../nextflow_run/index.md) o comunque essere a proprio agio nell'eseguire una pipeline locale con `nextflow run`."
---

# Execution Config

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduzione assistita da IA - [scopri di più e suggerisci miglioramenti](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

**Execution Config è un'introduzione pratica all'adattamento dell'esecuzione di pipeline Nextflow a diversi ambienti di calcolo.**

Attraverso esercizi orientati agli obiettivi, imparerete a cambiare tecnologia di packaging software, selezionare una piattaforma di esecuzione, controllare l'allocazione delle risorse di calcolo e i tentativi di ripetizione, e raggruppare la configurazione in profili intercambiabili.

Acquisirete le competenze e la sicurezza per configurare l'esecuzione di pipeline Nextflow come dei professionisti.

<!-- additional_information -->

## Panoramica del corso

Questo corso è pratico e si basa sulle competenze acquisite in [Nextflow Run](../nextflow_run/index.md).

Prenderete la stessa pipeline multi-step di quel corso e ne adatterete progressivamente la configurazione a diversi ambienti di calcolo, per poi raggruppare tutto in profili tra cui potrete passare in fase di esecuzione.

### Piano delle lezioni

| Capitolo del corso                                                             | Sommario                                                                                    | Durata stimata |
| ------------------------------------------------------------------------------ | ------------------------------------------------------------------------------------------- | -------------- |
| [Parte 1: Adattarsi all'ambiente di calcolo](./01_packaging_and_execution.md)  | Cambiare tecnologia di packaging software e selezionare una piattaforma di esecuzione       | 20 min         |
| [Parte 2: Gestire le risorse di calcolo e i fallimenti](./02_resources_and_retries.md) | Controllare l'allocazione delle risorse e riprovare automaticamente le attività che falliscono | 15 min         |
| [Parte 3: Usare i profili per cambiare configurazione](./03_profiles.md)       | Definire e combinare profili, e ispezionare la configurazione completamente risolta          | 15 min         |

Al termine di questo corso, sarete a vostro agio nel configurare pipeline Nextflow per una serie di ambienti di calcolo e nel passare da uno all'altro con il minimo sforzo.

Pronti a iniziare il corso?

[Inizia a imparare :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }
