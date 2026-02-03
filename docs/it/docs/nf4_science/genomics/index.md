---
title: Nextflow per la Genomica
hide:
  - toc
---

# Nextflow per la Genomica

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduzione assistita da IA - [scopri di più e suggerisci miglioramenti](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Questo corso di formazione è destinato ai ricercatori nel campo della genomica e settori correlati che sono interessati a sviluppare o personalizzare pipeline di analisi dei dati.
Si basa sulla formazione per principianti [Hello Nextflow](../../hello_nextflow/) e dimostra come utilizzare Nextflow nel contesto specifico del dominio della genomica.

In particolare, questo corso dimostra come implementare una semplice pipeline di variant calling con [GATK](https://gatk.broadinstitute.org/) (Genome Analysis Toolkit), un pacchetto software ampiamente utilizzato per analizzare dati di sequenziamento ad alto rendimento.

Cominciamo! Clicchi sul pulsante "Open in GitHub Codespaces" qui sotto per avviare l'ambiente di formazione (preferibilmente in una scheda separata), quindi continui a leggere durante il caricamento.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

## Obiettivi di apprendimento

Seguendo questo corso, imparerà ad applicare concetti e strumenti fondamentali di Nextflow a un caso d'uso tipico della genomica.

Al termine di questo workshop sarà in grado di:

- Scrivere un workflow lineare per applicare il variant calling a un singolo campione
- Gestire appropriatamente file accessori come file indice e risorse del genoma di riferimento
- Sfruttare il paradigma di dataflow di Nextflow per parallelizzare il variant calling per campione
- Implementare il variant calling multi-campione utilizzando operatori di canale rilevanti
- Implementare test di pipeline per singolo passaggio e end-to-end che gestiscono appropriatamente le idiosincrasie specifiche della genomica

<!-- TODO for future expansion: add metadata/samplesheet handling -->

## Prerequisiti

Il corso presuppone una familiarità minima con quanto segue:

- Strumenti e formati di file comunemente utilizzati in questo dominio scientifico
- Esperienza con la riga di comando
- Concetti e strumenti fondamentali di Nextflow trattati nella formazione per principianti [Hello Nextflow](../../hello_nextflow/)

Per i requisiti tecnici e la configurazione dell'ambiente, consulti il mini-corso [Configurazione dell'Ambiente](../../envsetup/).
