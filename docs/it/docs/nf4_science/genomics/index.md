---
title: Nextflow per la Genomica
hide:
  - toc
page_type: index_page
index_type: course
additional_information:
  technical_requirements: true
  learning_objectives:
    - Scrivere un workflow lineare per applicare il variant calling a un singolo campione
    - Gestire appropriatamente file accessori come file indice e risorse del genoma di riferimento
    - Sfruttare il paradigma di dataflow di Nextflow per parallelizzare il variant calling per campione
    - Implementare il joint calling multi-campione utilizzando operatori di canale rilevanti
  audience_prerequisites:
    - "**Pubblico:** Questo corso è progettato per ricercatori nel campo della genomica e settori correlati che desiderano sviluppare o personalizzare pipeline di analisi dei dati."
    - "**Competenze:** Si presume una certa familiarità con la riga di comando, concetti di scripting di base e formati di file comunemente utilizzati nella genomica."
    - "**Prerequisiti:** Concetti e strumenti fondamentali di Nextflow trattati in [Hello Nextflow](../../hello_nextflow/)."
---

# Nextflow per la Genomica

**Un corso pratico che applica Nextflow a un caso d'uso reale della genomica: variant calling con GATK.**

Questo corso si basa sulla formazione per principianti [Hello Nextflow](../../hello_nextflow/) e dimostra come utilizzare Nextflow nel contesto specifico del dominio della genomica.
Implementerete una pipeline di variant calling con [GATK](https://gatk.broadinstitute.org/) (Genome Analysis Toolkit), un pacchetto software ampiamente utilizzato per analizzare dati di sequenziamento ad alto rendimento.

<!-- additional_information -->

## Panoramica del corso

Questo corso è pratico, con esercizi orientati agli obiettivi strutturati per introdurre le informazioni in modo graduale.

Inizierete eseguendo manualmente gli strumenti di variant calling nel terminale per comprendere la metodologia, quindi costruirete progressivamente una pipeline Nextflow che automatizza e scala l'analisi.

### Piano delle lezioni

Abbiamo suddiviso questo corso in tre parti che si concentrano ciascuna su aspetti specifici dell'applicazione di Nextflow a un caso d'uso della genomica.

| Capitolo del corso                                                          | Sommario                                                                                                       | Durata stimata |
| --------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------- | -------------- |
| [Parte 1: Panoramica del metodo](./01_method.md)                            | Comprensione della metodologia di variant calling ed esecuzione manuale degli strumenti                        | 30 minuti      |
| [Parte 2: Variant calling per campione](./02_per_sample_variant_calling.md) | Costruzione di una pipeline che indicizza file BAM e chiama varianti, poi scala a più campioni                 | 60 minuti      |
| [Parte 3: Joint calling su una coorte](./03_joint_calling.md)               | Aggiunta del joint genotyping multi-campione utilizzando operatori di canale per aggregare output per campione | 45 minuti      |

Al termine di questo corso, sarete in grado di applicare concetti e strumenti fondamentali di Nextflow a un caso d'uso tipico della genomica.

Pronti a seguire il corso?

[Iniziamo :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }

<!-- Clearfix for float -->
<div style="content: ''; clear: both; display: table;"></div>
