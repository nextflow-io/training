# Riepilogo del corso

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduzione assistita da IA - [scopri di più e suggerisci miglioramenti](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Congratulazioni per aver completato il corso di formazione Use nf-core! 🎉

<!-- placeholder for video -->

## Il vostro percorso

Avete iniziato trovando e recuperando la pipeline `nf-core/demo`, poi avete imparato a eseguirla usando il suo profilo di test ed esaminarne gli output.
Successivamente, avete configurato la sua esecuzione tramite parametri della pipeline e file di configurazione, e avete visto come le pipeline nf-core validano i parametri e i dati di input.
Infine, avete applicato le stesse competenze a `nf-core/rnaseq`, una pipeline su scala produttiva, e avete imparato come sovrascrivere le allocazioni di risorse predefinite per adattarle all'hardware a vostra disposizione.

### Cosa avete imparato

Siete ora in grado di trovare, recuperare, eseguire e configurare le pipeline nf-core.

- Le pipeline nf-core vengono recuperate con `nextflow pull` e seguono un'organizzazione del codice standard.
- Ogni pipeline nf-core include un profilo `test` per una rapida validazione su un piccolo dataset.
- I parametri della pipeline (impostati tramite `--param_name` o `-params-file`) e la configurazione (impostata tramite `-c`) hanno scopi diversi: input e opzioni di analisi da un lato, logistica di esecuzione come l'allocazione delle risorse dall'altro.
- Le pipeline nf-core validano automaticamente i parametri e i file di input, intercettando gli errori prima che venga eseguita qualsiasi operazione.
- I valori predefiniti delle risorse vengono assegnati tramite label (`process_low`, `process_medium`, `process_high`) definiti in `conf/base.config`, che potete sovrascrivere con un file di configurazione personalizzato.

### Competenze acquisite

Attraverso questo corso pratico, avete imparato come:

- Trovare una pipeline nf-core sul sito nf-co.re e recuperarne il codice sorgente
- Eseguire una pipeline usando il suo profilo di test integrato ed esaminarne gli output
- Ottenere aiuto, impostare parametri e comprendere la validazione dei parametri e dell'input
- Personalizzare l'allocazione delle risorse e gli argomenti degli strumenti tramite file di configurazione
- Recuperare ed eseguire una pipeline su scala produttiva, e sovrascrivere i suoi label di risorse predefiniti

Siete ora dotati delle conoscenze fondamentali per iniziare a eseguire le pipeline nf-core per le vostre analisi.

## Prossimi passi per sviluppare le vostre competenze

Ecco i nostri principali suggerimenti su cosa fare dopo:

- Lanciate e monitorate queste pipeline su larga scala con [Scale with Seqera](../seqera_scale/index.md)
- Non limitatevi a eseguire le pipeline nf-core, sviluppatele! Imparate le best practice di nf-core con [Build with nf-core](../nfcore_build/index.md)
- Siete nuovi a Nextflow? Iniziate con [Nextflow Run](../nextflow_run/index.md)
- Applicate Nextflow a un caso d'uso di analisi scientifica con [Nextflow for Science](../nf4_science/index.md)
- Esplorate funzionalità più avanzate di Nextflow con le [Side Quests](../side_quests/index.md)

## Ottenere aiuto

Per risorse di supporto e assistenza dalla community, consultate la [pagina di aiuto](../help.md).

## Sondaggio di feedback

Prima di proseguire, dedicate un minuto a completare il sondaggio del corso! Il vostro feedback ci aiuta a migliorare i materiali di formazione per tutti.

[Partecipa al sondaggio :material-arrow-right:](survey.md){ .md-button .md-button--primary }
