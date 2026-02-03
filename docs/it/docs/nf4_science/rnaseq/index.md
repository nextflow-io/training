---
title: Nextflow per RNAseq
hide:
  - toc
---

# Nextflow per RNAseq

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduzione assistita da IA - [scopri di più e suggerisci miglioramenti](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Questo corso di formazione è destinato a ricercatori in trascrittomica e campi correlati interessati allo sviluppo o alla personalizzazione di pipeline di analisi dati.
Si basa sulla formazione per principianti [Hello Nextflow](../../hello_nextflow/) e dimostra come utilizzare Nextflow nel contesto specifico dell'analisi RNAseq bulk.

In particolare, questo corso dimostra come implementare una semplice pipeline di elaborazione RNAseq bulk per rimuovere le sequenze adattatrici, allineare le letture a un genoma di riferimento ed eseguire il controllo qualità (QC) in diverse fasi.

Cominciamo! Clicchi sul pulsante "Open in GitHub Codespaces" qui sotto per avviare l'ambiente di formazione (preferibilmente in una scheda separata), quindi continui a leggere mentre si carica.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

## Obiettivi di apprendimento

Seguendo questo corso, imparerà ad applicare i concetti e gli strumenti fondamentali di Nextflow a un caso d'uso tipico di RNAseq.

Al termine di questo workshop sarà in grado di:

- Scrivere un workflow lineare per applicare metodi di elaborazione e QC di base per RNAseq
- Gestire appropriatamente file specifici del dominio come FASTQ e risorse di genomi di riferimento
- Gestire dati di sequenziamento single-end e paired-end
- Sfruttare il paradigma dataflow di Nextflow per parallelizzare l'elaborazione RNAseq per campione
- Aggregare report di QC attraverso più fasi e campioni utilizzando operatori di canale appropriati

<!-- TODO
- Configure pipeline execution and manage and optimize resource allocations
- Implement per-step and end-to-end pipeline tests that handle RNAseq-specific idiosyncrasies appropriately
-->
<!-- TODO for future expansion: add metadata/samplesheet handling -->

## Prerequisiti

Il corso presuppone una minima familiarità con quanto segue:

- Strumenti e formati di file comunemente utilizzati in questo dominio scientifico
- Esperienza con la riga di comando
- Concetti e strumenti fondamentali di Nextflow trattati nella formazione per principianti [Hello Nextflow](../../hello_nextflow/)

Per i requisiti tecnici e la configurazione dell'ambiente, consulti il mini-corso [Configurazione dell'Ambiente](../../envsetup/).
