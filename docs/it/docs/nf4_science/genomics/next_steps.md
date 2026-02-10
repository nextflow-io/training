# Riepilogo del corso

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduzione assistita da IA - [scopri di più e suggerisci miglioramenti](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Congratulazioni per aver completato il corso di formazione Nextflow for Genomics! 🎉

## Il vostro percorso

Avete iniziato eseguendo manualmente gli strumenti di variant calling nel terminale per comprendere la metodologia.
Poi avete costruito una pipeline Nextflow per un singolo campione per automatizzare il processo, l'avete scalata per gestire più campioni in parallelo e avete aggiunto il joint genotyping multi-campione utilizzando gli operatori di canale.

### Cosa avete costruito

- Una pipeline di variant calling che prende file BAM come input e produce VCF joint-called come output.
- Tre processi (`SAMTOOLS_INDEX`, `GATK_HAPLOTYPECALLER` e `GATK_JOINTGENOTYPING`) memorizzati in file modulo separati.
- La pipeline scala automaticamente per qualsiasi numero di campioni di input utilizzando il paradigma dataflow di Nextflow.
- I risultati vengono pubblicati in una directory chiamata `results/`.

### Competenze acquisite

Attraverso questo corso pratico, avete imparato a:

- Scrivere un flusso di lavoro lineare per applicare il variant calling a un singolo campione
- Gestire appropriatamente i file accessori come i file indice e le risorse del genoma di riferimento
- Sfruttare il paradigma dataflow di Nextflow per parallelizzare il variant calling per campione
- Implementare il joint calling multi-campione utilizzando gli operatori di canale rilevanti

Ora siete pronti per iniziare ad applicare Nextflow ai flussi di lavoro di analisi genomica nel vostro lavoro.

## Prossimi passi per sviluppare le vostre competenze

Ecco i nostri principali suggerimenti su cosa fare dopo:

- Applicate Nextflow ad altri casi d'uso di analisi scientifica con [Nextflow for Science](../index.md)
- Iniziate con nf-core con [Hello nf-core](../../hello_nf-core/index.md)
- Esplorate funzionalità più avanzate di Nextflow con le [Side Quests](../../side_quests/index.md)

Infine, vi consigliamo di dare un'occhiata a [**Seqera Platform**](https://seqera.io/), una piattaforma basata su cloud sviluppata dai creatori di Nextflow che rende ancora più facile lanciare e gestire i vostri flussi di lavoro, oltre a gestire i vostri dati ed eseguire analisi in modo interattivo in qualsiasi ambiente.

## Ottenere aiuto

Per risorse di aiuto e supporto della comunità, consultate la [pagina di aiuto](../../help.md).

## Sondaggio di feedback

Prima di proseguire, dedicate un minuto a completare il sondaggio del corso! Il vostro feedback ci aiuta a migliorare i nostri materiali di formazione per tutti.

[Partecipa al sondaggio :material-arrow-right:](survey.md){ .md-button .md-button--primary }
