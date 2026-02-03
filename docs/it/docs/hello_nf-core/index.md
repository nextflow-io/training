---
title: Hello nf-core
hide:
  - toc
page_type: index_page
index_type: course
additional_information:
  technical_requirements: true
  learning_objectives:
    - Recuperare, lanciare e gestire l'esecuzione di pipeline nf-core
    - Descrivere la struttura del codice e l'organizzazione dei progetti delle pipeline nf-core
    - Creare una pipeline di base compatibile con nf-core da un template
    - Aggiornare un workflow Nextflow semplice per adattarlo agli standard nf-core
    - Aggiungere moduli nf-core a una pipeline compatibile con nf-core
    - Contribuire con i propri moduli a nf-core
    - Validare input e parametri utilizzando gli strumenti nf-core
  audience_prerequisites:
    - "**Pubblico:** Questo corso è progettato per studenti che hanno già familiarità con Nextflow di base e vogliono imparare a utilizzare le risorse e le best practice di nf-core."
    - "**Competenze:** Si presume familiarità con la riga di comando, concetti di scripting di base e formati di file comuni."
    - "**Corsi:** È necessario aver completato il corso [Hello Nextflow](../hello_nextflow/index.md) o equivalente."
    - "**Dominio:** Gli esercizi sono tutti indipendenti dal dominio, quindi non è richiesta alcuna conoscenza scientifica preliminare."
---

# Hello nf-core

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduzione assistita da IA - [scopri di più e suggerisci miglioramenti](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

**Hello nf-core è un'introduzione pratica all'uso delle risorse e delle best practice di nf-core.**

![nf-core logo](./img/nf-core-logo.png)

Attraverso esempi pratici ed esercizi guidati, imparerete a utilizzare e sviluppare moduli e pipeline compatibili con nf-core, e a utilizzare efficacemente gli strumenti nf-core.

Acquisirà le competenze e la sicurezza per iniziare a sviluppare pipeline secondo le best practice di nf-core.

<!-- additional_information -->

## Panoramica del corso

Questo corso è progettato per essere pratico, con esercizi orientati agli obiettivi strutturati per introdurre le informazioni gradualmente.

Verrà introdotto a [**nf-core**](https://nf-co.re/), uno sforzo comunitario per sviluppare e mantenere un insieme curato di pipeline scientifiche costruite utilizzando Nextflow, oltre a strumenti e linee guida rilevanti che promuovono lo sviluppo aperto, i test e la revisione tra pari ([Nat Biotechnol 38, 276–278 (2020)](https://www.nature.com/articles/s41587-020-0439-x), [Genome Biol 26, 228 (2025)](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-025-03673-9)).

Le pipeline sviluppate dalla comunità nf-core sono progettate per essere modulari, scalabili e portabili, consentendo ai ricercatori di adattarle ed eseguirle facilmente utilizzando i propri dati e risorse di calcolo.
Le linee guida sulle best practice applicate dal progetto garantiscono inoltre che le pipeline siano robuste, ben documentate e validate su dataset reali.
Questo aiuta ad aumentare l'affidabilità e la riproducibilità delle analisi scientifiche e, in definitiva, consente ai ricercatori di accelerare le proprie scoperte scientifiche.

Non copriremo tutto ciò che c'è da sapere sulle pipeline nf-core in questo corso, perché nf-core comprende molte funzionalità e convenzioni sviluppate dalla comunità nel corso degli anni.
Ci concentreremo invece sui concetti essenziali che La aiuteranno a iniziare e a capire come funziona nf-core.

### Piano delle lezioni

Abbiamo suddiviso il corso in cinque parti, ciascuna focalizzata su aspetti specifici dell'utilizzo delle risorse nf-core.

| Capitolo del corso                                                        | Riepilogo                                                                                                                                                         | Durata stimata |
| ------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------- | -------------- |
| [Parte 1: Eseguire una pipeline demo](./01_run_demo.md)                   | Eseguire una pipeline nf-core esistente ed esaminare la sua struttura di codice per comprendere cosa rende queste pipeline diverse dai workflow Nextflow di base  | 30 minuti      |
| [Parte 2: Riscrivere Hello per nf-core](./02_rewrite_hello.md)            | Adattare un workflow esistente allo scaffold del template nf-core, partendo dal workflow semplice prodotto nel corso [Hello Nextflow](../hello_nextflow/index.md) | 60 minuti      |
| [Parte 3: Usare un modulo nf-core](./03_use_module.md)                    | Esplorare la libreria di moduli della comunità e imparare a integrare moduli pre-costruiti e testati che incapsulano strumenti bioinformatici comuni              | 30 minuti      |
| [Parte 4: Creare un modulo nf-core](./04_make_module.md)                  | Creare il proprio modulo in stile nf-core utilizzando la struttura specifica, le convenzioni di denominazione e i requisiti di metadati stabiliti da nf-core      | 30 minuti      |
| [Parte 5: Aggiungere la validazione dell'input](./05_input_validation.md) | Implementare la validazione dell'input sia per i parametri da riga di comando che per i file di dati di input utilizzando nf-schema                               | 30 minuti      |

Al termine di questo corso, sarà in grado di sfruttare l'enorme ricchezza di risorse offerte dal progetto nf-core.

Pronto a iniziare il corso?

[Inizia a imparare :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }
