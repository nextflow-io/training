---
title: Nextflow per la Genomica
hide:
  - toc
page_type: index_page
index_type: course
additional_information:
  technical_requirements: true
  learning_objectives:
    - Scrivere un flusso di lavoro lineare per applicare il variant calling a un singolo campione
    - Gestire appropriatamente i file accessori come i file di indice e le risorse del genoma di riferimento
    - Sfruttare il paradigma dataflow di Nextflow per parallelizzare il variant calling per campione
    - Implementare il joint calling multi-campione utilizzando gli operatori di canale rilevanti
  audience_prerequisites:
    - "**Pubblico:** Questo corso è progettato per ricercatori in genomica e campi correlati che desiderano sviluppare o personalizzare pipeline di analisi dati."
    - "**Competenze:** Si presume una certa familiarità con la riga di comando, concetti di scripting di base e formati di file genomici comuni."
    - "**Prerequisiti:** Concetti fondamentali di Nextflow e strumenti trattati in [Hello Nextflow](../../hello_nextflow/)."
---

# Nextflow per la Genomica

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduzione assistita da IA - [scopri di più e suggerisci miglioramenti](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

**Un corso pratico che applica Nextflow a un caso d'uso reale di genomica: variant calling con GATK.**

Questo corso si basa sulla formazione per principianti [Hello Nextflow](../../hello_nextflow/) e dimostra come utilizzare Nextflow nel contesto specifico del dominio della genomica.
Implementerete una pipeline di variant calling con [GATK](https://gatk.broadinstitute.org/) (Genome Analysis Toolkit), un pacchetto software ampiamente utilizzato per l'analisi di dati di sequenziamento ad alto rendimento.

<!-- additional_information -->

## Panoramica del corso

Questo corso è pratico, con esercizi orientati agli obiettivi strutturati per introdurre le informazioni gradualmente.

Inizierete eseguendo manualmente gli strumenti di variant calling nel terminale per comprendere la metodologia, quindi costruirete progressivamente una pipeline Nextflow che automatizza e scala l'analisi.

### Piano delle lezioni

Abbiamo suddiviso il corso in tre parti che si concentrano ciascuna su aspetti specifici dell'applicazione di Nextflow a un caso d'uso di genomica.

| Capitolo del corso                                                          | Riepilogo                                                                                                          | Durata stimata |
| --------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------ | -------------- |
| [Parte 1: Panoramica del metodo](./01_method.md)                            | Comprensione della metodologia di variant calling ed esecuzione manuale degli strumenti                            | 30 min         |
| [Parte 2: Variant calling per campione](./02_per_sample_variant_calling.md) | Costruzione di una pipeline che indicizza i file BAM e chiama le varianti, quindi scala a più campioni             | 60 min         |
| [Parte 3: Joint calling su una coorte](./03_joint_calling.md)               | Aggiunta del joint genotyping multi-campione utilizzando operatori di canale per aggregare gli output per campione | 45 min         |

Al termine di questo corso, sarete in grado di applicare concetti fondamentali di Nextflow e strumenti a un tipico caso d'uso di genomica.

Pronti per iniziare il corso?

[Inizia :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }

<!-- Clearfix for float -->
<div style="content: ''; clear: both; display: table;"></div>
