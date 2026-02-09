# Nextflow run for Imaging

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduzione assistita da IA - [maggiori informazioni e suggerimenti per miglioramenti](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Questo corso di formazione è destinato a ricercatori nel campo dell'imaging e della biologia spaziale interessati a eseguire e personalizzare pipeline di analisi dati.
Insegna i concetti fondamentali di Nextflow relativi all'esecuzione, organizzazione e configurazione dei flussi di lavoro utilizzando [nf-core/molkart](https://nf-co.re/molkart), una pipeline per l'elaborazione di dati di trascrittomica spaziale Molecular Cartography.
Le competenze che acquisirete qui sono trasferibili a qualsiasi pipeline Nextflow o nf-core.

Iniziamo! Cliccate sul pulsante "Open in GitHub Codespaces" qui sotto per avviare l'ambiente di formazione (preferibilmente in una scheda separata), poi continuate a leggere mentre si carica.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

## Obiettivi di apprendimento

Lavorando attraverso questo corso, imparerete ad applicare i concetti e gli strumenti fondamentali di Nextflow all'esecuzione di pipeline di analisi di imaging.

Al termine di questo workshop sarete in grado di:

- Lanciare un flusso di lavoro Nextflow localmente e monitorarne l'esecuzione
- Trovare e interpretare gli output (risultati) e i file di log generati da Nextflow
- Eseguire una pipeline nf-core con dati di test e input personalizzati
- Configurare l'esecuzione della pipeline utilizzando profili e file di parametri
- Gestire gli input utilizzando samplesheet e parametri da riga di comando

## Pubblico e prerequisiti

Questo corso presuppone una minima familiarità con quanto segue:

- Esperienza con la riga di comando
- Familiarità di base con i formati di file di imaging (immagini TIFF, dati tabulari)

Per i requisiti tecnici e la configurazione dell'ambiente, consultate il mini-corso [Environment Setup](../../envsetup/).
