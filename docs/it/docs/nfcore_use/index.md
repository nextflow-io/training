---
title: Usare nf-core
hide:
  - toc
page_type: index_page
index_type: course
additional_information:
  technical_requirements: true
  learning_objectives:
    - Trovare, recuperare ed eseguire pipeline della community nf-core
    - Configurare l'esecuzione della pipeline usando parametri e file di configurazione
    - Capire come le pipeline nf-core validano i parametri e i dati di input
    - Eseguire una pipeline su scala produttiva (nf-core/rnaseq) e sovrascrivere le sue allocazioni di risorse predefinite
  audience_prerequisites:
    - "**Pubblico:** Questo corso è pensato per chi sa già eseguire pipeline Nextflow in locale, è alle prime armi con nf-core e vuole eseguire pipeline della community già esistenti."
    - "**Competenze:** Si presuppone una certa familiarità con la riga di comando, concetti base di scripting e i formati di file più comuni."
    - "**Corsi:** È necessario aver completato [Nextflow Run](../nextflow_run/index.md) o comunque essere a proprio agio nell'eseguire una pipeline locale con `nextflow run`."
    - "**Dominio:** Gli esercizi utilizzano pipeline bioinformatiche, ma non è richiesta alcuna conoscenza scientifica pregressa."
---

# Usare nf-core

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduzione assistita da IA - [scopri di più e suggerisci miglioramenti](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

**Usare nf-core è un'introduzione pratica alla ricerca, all'esecuzione e alla configurazione delle pipeline della community nf-core.**

Attraverso esempi pratici ed esercizi guidati, imparerete a trovare e recuperare le pipeline nf-core, a eseguirle usando i loro profili di test integrati e a personalizzarne l'esecuzione tramite parametri e file di configurazione.

Acquisirete le competenze e la sicurezza necessarie per iniziare a eseguire pipeline nf-core nelle vostre analisi.

<!-- additional_information -->

## Panoramica del corso

Questo corso è pratico, con esercizi orientati agli obiettivi e strutturati per introdurre le informazioni in modo graduale.

Inizierete con `nf-core/demo`, una pipeline minimale mantenuta dal progetto nf-core a scopo formativo, per poi applicare quanto appreso a `nf-core/rnaseq`, una pipeline di produzione ampiamente utilizzata per l'analisi del sequenziamento RNA bulk.

Questo corso si concentra sull'esecuzione delle pipeline.
Se cercate un'introduzione allo sviluppo di pipeline compatibili con nf-core, consultate [Build with nf-core](../nfcore_build/index.md).

### Piano delle lezioni

| Capitolo del corso                                                   | Sommario                                                                                                          | Durata stimata |
| -------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------- | -------------- |
| [Parte 1: Eseguire una pipeline demo](./01_run_demo.md)              | Trovare e recuperare una pipeline nf-core ed eseguirla usando il suo profilo di test                              | 20 min         |
| [Parte 2: Configurare l'esecuzione della pipeline](./02_configure_execution.md) | Impostare i parametri, capire la validazione e personalizzare l'allocazione delle risorse e gli argomenti degli strumenti | 20 min |
| [Parte 3: Eseguire una pipeline di produzione](./03_run_production_pipeline.md) | Scaricare ed eseguire nf-core/rnaseq e sovrascrivere le sue allocazioni di risorse predefinite              | 20 min         |

Al termine di questo corso, sarete in grado di sfruttare la ricchezza di pipeline della community offerte dal progetto nf-core.

Pronti a iniziare il corso?

[Inizia a imparare :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }
