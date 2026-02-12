# Riepilogo del corso

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduzione assistita da IA - [scopri di più e suggerisci miglioramenti](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Congratulazioni per aver completato il corso di formazione Nextflow per RNAseq!

## Il vostro percorso

Avete iniziato eseguendo manualmente gli strumenti di elaborazione RNAseq nel terminale per comprendere la metodologia.
Poi avete costruito una pipeline Nextflow per un singolo campione per automatizzare il processo, l'avete scalata per gestire più campioni in parallelo e l'avete estesa per gestire dati paired-end e aggregare report QC tra i campioni.

### Cosa avete costruito

- Una pipeline di elaborazione RNAseq che prende file FASTQ come input e produce read trimmate, allineamenti e report QC aggregati come output.
- Processi per il trimming (Trim Galore), l'allineamento (HISAT2), il controllo qualità (FastQC) e l'aggregazione dei report (MultiQC) memorizzati in file modulo separati.
- La pipeline parallelizza automaticamente l'elaborazione dei campioni di input utilizzando il paradigma dataflow di Nextflow.
- La pipeline finale gestisce dati di sequenziamento paired-end.

### Competenze acquisite

Attraverso questo corso pratico, avete imparato come:

- Scrivere un flusso di lavoro lineare per applicare metodi di base di elaborazione e QC RNAseq
- Gestire appropriatamente file specifici del dominio come FASTQ e risorse del genoma di riferimento
- Gestire dati di sequenziamento single-end e paired-end
- Sfruttare il paradigma dataflow di Nextflow per parallelizzare l'elaborazione RNAseq per campione
- Aggregare report QC attraverso più passaggi e campioni utilizzando operatori di canale pertinenti

Ora siete pronti per iniziare ad applicare Nextflow ai flussi di lavoro di analisi RNAseq nel vostro lavoro.

## Prossimi passi per sviluppare le vostre competenze

Ecco i nostri suggerimenti principali su cosa fare dopo:

- Applicare Nextflow ad altri casi d'uso di analisi scientifica con [Nextflow for Science](../index.md)
- Iniziare con nf-core con [Hello nf-core](../../hello_nf-core/index.md)
- Esplorare funzionalità più avanzate di Nextflow con le [Side Quests](../../side_quests/index.md)

Infine, vi consigliamo di dare un'occhiata a [**Seqera Platform**](https://seqera.io/), una piattaforma basata su cloud sviluppata dai creatori di Nextflow che rende ancora più semplice lanciare e gestire i vostri flussi di lavoro, nonché gestire i vostri dati ed eseguire analisi in modo interattivo in qualsiasi ambiente.

## Ottenere aiuto

Per risorse di aiuto e supporto della community, consultate la [pagina Aiuto](../../help.md).

## Questionario di feedback

Prima di proseguire, dedicate un minuto a completare il questionario del corso! Il vostro feedback ci aiuta a migliorare i nostri materiali di formazione per tutti.

[Compila il questionario :material-arrow-right:](survey.md){ .md-button .md-button--primary }
