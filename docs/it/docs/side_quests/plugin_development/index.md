---
title: Sviluppo di Plugin
hide:
  - toc
---

# Sviluppo di Plugin

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduzione assistita da IA - [scopri di più e suggerisci miglioramenti](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Il sistema di plugin di Nextflow consente di estendere il linguaggio con funzioni personalizzate, hook di monitoraggio, backend di esecuzione e molto altro.
I plugin permettono alla community di aggiungere funzionalità a Nextflow senza modificarne il nucleo, rendendoli ideali per condividere funzionalità riutilizzabili tra pipeline.

Durante questa formazione, imparerete come utilizzare plugin esistenti e, facoltativamente, come crearne di propri.

## Pubblico e prerequisiti

La Parte 1 riguarda l'utilizzo di plugin esistenti ed è rilevante per tutti gli utenti di Nextflow.
Le Parti 2-6 riguardano la creazione di plugin personalizzati e coinvolgono codice Groovy e strumenti di build.
Non è richiesta alcuna esperienza precedente con Java o Groovy.

**Prerequisiti**

- Un account GitHub OPPURE un'installazione locale come descritto [qui](../../envsetup/02_local).
- Aver completato il corso [Hello Nextflow](../../hello_nextflow/index.md) o equivalente.
- Java 21 o versione successiva (incluso nell'ambiente di formazione; necessario solo per le Parti 2-6).

**Directory di lavoro:** `side-quests/plugin_development`

## Obiettivi di apprendimento

Al termine di questa formazione, sarete in grado di:

**Utilizzo dei plugin (Parte 1):**

- Installare e configurare plugin esistenti nei vostri flussi di lavoro
- Importare e utilizzare le funzioni dei plugin

**Sviluppo di plugin (Parti 2-6):**

- Creare un nuovo progetto plugin utilizzando il generatore di progetti integrato di Nextflow
- Implementare funzioni personalizzate richiamabili dai flussi di lavoro
- Compilare, testare e installare il plugin in locale
- Monitorare gli eventi del flusso di lavoro (ad esempio, completamento di un'attività, avvio/fine della pipeline) per log personalizzati o notifiche
- Aggiungere opzioni di configurazione per rendere i plugin personalizzabili
- Distribuire il proprio plugin

## Piano delle lezioni

#### Parte 1: Nozioni di base sui plugin

Utilizzare plugin esistenti in un flusso di lavoro Nextflow e configurarne il comportamento.

#### Parte 2: Creare un progetto plugin

Generare un nuovo progetto plugin ed esaminarne la struttura.

#### Parte 3: Funzioni personalizzate

Implementare funzioni personalizzate, compilare il plugin ed eseguirlo in un flusso di lavoro.

#### Parte 4: Testing

Scrivere ed eseguire unit test utilizzando il framework Spock.

#### Parte 5: Monitoraggio del flusso di lavoro

Rispondere a eventi come il completamento di un'attività per costruire un contatore di attività.

#### Parte 6: Configurazione e distribuzione

Leggere le impostazioni da `nextflow.config` per rendere il plugin personalizzabile, quindi imparare come condividerlo.

Pronti a iniziare il corso?

[Inizia a imparare :material-arrow-right:](01_plugin_basics.md){ .md-button .md-button--primary }
