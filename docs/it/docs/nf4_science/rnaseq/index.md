---
title: RNAseq
hide:
  - toc
page_type: index_page
index_type: course
additional_information:
  technical_requirements: true
  learning_objectives:
    - Scrivere un flusso di lavoro lineare per applicare metodi di elaborazione e QC di base per RNAseq
    - Gestire appropriatamente file specifici del dominio come FASTQ e risorse di genomi di riferimento
    - Gestire dati di sequenziamento single-end e paired-end
    - Sfruttare il paradigma dataflow di Nextflow per parallelizzare l'elaborazione RNAseq per campione
    - Aggregare report di QC attraverso più fasi e campioni utilizzando operatori di canale rilevanti
  audience_prerequisites:
    - "**Pubblico:** Questo corso è progettato per ricercatori in trascrittomica e campi correlati che desiderano sviluppare o personalizzare pipeline di analisi dati."
    - "**Competenze:** Si presuppone una certa familiarità con la riga di comando, concetti di scripting di base e formati di file RNAseq comuni."
    - "**Prerequisiti:** Concetti e strumenti fondamentali di Nextflow trattati in [Hello Nextflow](../../hello_nextflow/)."
---

# Nextflow per RNAseq

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduzione assistita da IA - [scopri di più e suggerisci miglioramenti](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

**Un corso pratico che applica Nextflow a un caso d'uso reale di trascrittomica: elaborazione di RNAseq bulk con Trim Galore, HISAT2 e FastQC.**

Questo corso si basa sulla formazione per principianti [Hello Nextflow](../../hello_nextflow/) e dimostra come utilizzare Nextflow nel contesto specifico dell'analisi RNAseq bulk.
Implementerete una pipeline di elaborazione che rimuove le sequenze adattatrici, allinea le letture a un genoma di riferimento ed esegue il controllo qualità (QC) in diverse fasi.

<!-- additional_information -->

## Panoramica del corso

Questo corso è pratico, con esercizi orientati agli obiettivi strutturati per introdurre le informazioni gradualmente.

Inizierete eseguendo manualmente gli strumenti di elaborazione nel terminale per comprendere la metodologia, quindi costruirete progressivamente una pipeline Nextflow che automatizza e scala l'analisi.

### Piano delle lezioni

Abbiamo suddiviso il corso in tre parti che si concentrano ciascuna su aspetti specifici dell'applicazione di Nextflow a un caso d'uso RNAseq.

| Capitolo del corso                                                         | Riepilogo                                                                                                                   | Durata stimata |
| -------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------- | -------------- |
| [Parte 1: Panoramica del metodo](./01_method.md)                           | Comprensione della metodologia di elaborazione RNAseq ed esecuzione manuale degli strumenti                                 | 30 minuti      |
| [Parte 2: Implementazione per singolo campione](./02_single-sample.md)     | Costruzione di una pipeline che elabora, allinea e controlla la qualità di un singolo campione, quindi scala a più campioni | 60 minuti      |
| [Parte 3: Implementazione multi-campione paired-end](./03_multi-sample.md) | Estensione della pipeline per gestire dati paired-end e aggregare report di QC attraverso i campioni                        | 45 minuti      |

Al termine di questo corso, sarete in grado di applicare concetti e strumenti fondamentali di Nextflow a un caso d'uso tipico di RNAseq.

Pronti per iniziare il corso?

[Inizia :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }

<!-- Clearfix for float -->
<div style="content: ''; clear: both; display: table;"></div>
