---
title: Nextflow per RNAseq
hide:
  - toc
---

# Nextflow per RNAseq

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduzione assistita da IA - [maggiori informazioni e suggerimenti per miglioramenti](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Questo corso di formazione è destinato a ricercatori nel campo della trascrittomica e discipline correlate che sono interessati a sviluppare o personalizzare pipeline di analisi dati.
Si basa sul corso per principianti [Hello Nextflow](../../hello_nextflow/) e dimostra come utilizzare Nextflow nel contesto specifico dell'analisi RNAseq bulk.

In particolare, questo corso dimostra come implementare una semplice pipeline di elaborazione RNAseq bulk per rimuovere le sequenze adattatrici, allineare le letture a un genoma di riferimento ed eseguire il controllo qualità (QC) in diverse fasi.

Iniziamo! Cliccate sul pulsante "Open in GitHub Codespaces" qui sotto per avviare l'ambiente di formazione (preferibilmente in una scheda separata), poi continuate a leggere mentre si carica.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

## Obiettivi di apprendimento

Lavorando attraverso questo corso, imparerete come applicare i concetti e gli strumenti fondamentali di Nextflow a un caso d'uso tipico di RNAseq.

Al termine di questo workshop sarete in grado di:

- Scrivere un flusso di lavoro lineare per applicare metodi di base di elaborazione e QC per RNAseq
- Gestire appropriatamente file specifici del dominio come FASTQ e risorse del genoma di riferimento
- Gestire dati di sequenziamento single-end e paired-end
- Sfruttare il paradigma dataflow di Nextflow per parallelizzare l'elaborazione RNAseq per campione
- Aggregare report QC attraverso più fasi e campioni utilizzando operatori di canale rilevanti

<!-- TODO
- Configure pipeline execution and manage and optimize resource allocations
- Implement per-step and end-to-end pipeline tests that handle RNAseq-specific idiosyncrasies appropriately
-->
<!-- TODO for future expansion: add metadata/samplesheet handling -->

## Prerequisiti

Il corso presuppone una minima familiarità con quanto segue:

- Strumenti e formati di file comunemente utilizzati in questo dominio scientifico
- Esperienza con la riga di comando
- Concetti e strumenti fondamentali di Nextflow trattati nel corso per principianti [Hello Nextflow](../../hello_nextflow/)

Per i requisiti tecnici e la configurazione dell'ambiente, consultate il mini-corso [Configurazione dell'Ambiente](../../envsetup/).
