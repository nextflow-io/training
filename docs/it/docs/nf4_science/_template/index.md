---
title: Nextflow per {DOMAIN}
hide:
  - toc
page_type: index_page
index_type: course
additional_information:
  technical_requirements: true
  learning_objectives:
    - Scrivere un flusso di lavoro lineare per applicare {METHOD} a un singolo campione
    - Gestire file accessori come {ACCESSORY_FILES} in modo appropriato
    - Sfruttare il paradigma dataflow di Nextflow per parallelizzare l'elaborazione per campione
    - Implementare l'aggregazione multi-campione utilizzando gli operatori di canale rilevanti
  audience_prerequisites:
    - "**Pubblico:** Questo corso è progettato per ricercatori in {DOMAIN} e campi correlati che desiderano sviluppare o personalizzare pipeline di analisi dati."
    - "**Competenze:** Si presume una certa familiarità con la riga di comando, concetti di scripting di base e formati di file comuni in {DOMAIN}."
    - "**Prerequisiti:** Concetti fondamentali di Nextflow e strumenti trattati in [Hello Nextflow](../../hello_nextflow/)."
---

# Nextflow per {DOMAIN}

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduzione assistita da IA - [scopri di più e suggerisci miglioramenti](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

**Un corso pratico che applica Nextflow a un caso d'uso reale in {DOMAIN}: {METHOD_SHORT_DESCRIPTION}.**

Questo corso si basa sulla formazione per principianti [Hello Nextflow](../../hello_nextflow/) e dimostra come utilizzare Nextflow nel contesto specifico del dominio {DOMAIN}.
Implementerete una pipeline {METHOD} con [{TOOL_A}]({TOOL_A_URL}) e [{TOOL_B}]({TOOL_B_URL}).

<!-- additional_information -->

## Panoramica del corso

Questo corso è pratico, con esercizi orientati agli obiettivi strutturati per introdurre le informazioni gradualmente.

Inizierete eseguendo gli strumenti di analisi manualmente nel terminale per comprendere la metodologia, quindi costruirete progressivamente una pipeline Nextflow che automatizza e scala l'analisi.

### Piano delle lezioni

Abbiamo suddiviso il corso in tre parti che si concentrano ciascuna su aspetti specifici dell'applicazione di Nextflow a un caso d'uso in {DOMAIN}.

| Capitolo del corso                                        | Riepilogo                                                                                                      | Durata stimata |
| --------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------- | -------------- |
| [Parte 1: Panoramica del metodo](./01_method.md)          | Comprensione della metodologia {METHOD} ed esecuzione manuale degli strumenti                                  | 30 min         |
| [Parte 2: Elaborazione singolo campione](./02_single_sample.md) | Costruzione di una pipeline che {PART2_SUMMARY}, quindi scalatura a più campioni                               | 60 min         |
| [Parte 3: Aggregazione multi-campione](./03_multi_sample.md)   | Aggiunta di {AGGREGATION_SUMMARY} multi-campione utilizzando operatori di canale per aggregare output per campione | 45 min         |

Al termine di questo corso, sarete in grado di applicare concetti fondamentali di Nextflow e strumenti a un caso d'uso tipico in {DOMAIN}.

Pronti per iniziare il corso?

[Iniziamo :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }

<!-- Clearfix for float -->
<div style="content: ''; clear: both; display: table;"></div>
