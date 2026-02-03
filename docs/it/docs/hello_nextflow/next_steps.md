# Riepilogo del corso

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduzione assistita da IA - [scopri di più e suggerisci miglioramenti](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Congratulazioni per aver completato il corso di formazione Hello Nextflow! 🎉

<!--
<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/xHOcx_4Ancg?si=Lp8hS8RdaMwbp5j5&amp;list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } Guarda [l'intera playlist sul canale YouTube di Nextflow](https://www.youtube.com/playlist?list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik).

:green_book: Potete leggere la [trascrizione del video](./transcripts/07_next_steps.md) insieme al video.
///
-->

## Il vostro percorso

Avete iniziato con un workflow molto semplice che eseguiva un comando codificato in modo fisso.
Nel corso delle sei parti, avete trasformato quel workflow di base in una pipeline modulare multi-step che esercita le funzionalità chiave di Nextflow inclusi channel, operatori, supporto integrato per container e opzioni di configurazione.

### Cosa avete costruito

- La forma finale del workflow Hello prende come input un file CSV contenente saluti testuali.
- I quattro step sono implementati come processi Nextflow (`sayHello`, `convertToUpper`, `collectGreetings` e `cowpy`) memorizzati in file di modulo separati.
- I risultati vengono pubblicati in una directory chiamata `results/`.
- L'output finale della pipeline è un file di testo semplice contenente arte ASCII di un personaggio che dice i saluti in maiuscolo.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_pipeline_complete.svg"
</figure>

1. **`sayHello`:** Scrive ogni saluto nel proprio file di output (ad es. "Hello-output.txt")
2. **`convertToUpper`:** Converte ogni saluto in maiuscolo (ad es. "HELLO")
3. **`collectGreetings`:** Raccoglie tutti i saluti in maiuscolo in un singolo file batch
4. **`cowpy`:** Genera arte ASCII utilizzando lo strumento `cowpy`

La configurazione del workflow supporta la fornitura di input e parametri in modo flessibile e riproducibile.

### Competenze acquisite

Attraverso questo corso pratico, avete imparato a:

- Descrivere e utilizzare i componenti fondamentali di Nextflow sufficienti per costruire un workflow multi-step semplice
- Descrivere concetti successivi come operatori e channel factory
- Avviare un workflow Nextflow localmente
- Trovare e interpretare output (risultati) e file di log generati da Nextflow
- Risolvere problemi di base

Ora siete equipaggiati con le conoscenze fondamentali per iniziare a sviluppare le vostre pipeline in Nextflow.

## Prossimi passi per sviluppare le vostre competenze

Ecco i nostri 3 principali suggerimenti su cosa fare dopo:

- Applicare Nextflow a un caso d'uso di analisi scientifica con [Nextflow for Science](../nf4_science/index.md)
- Iniziare con nf-core con [Hello nf-core](../../hello_nf-core/index.md)
- Esplorare funzionalità Nextflow più avanzate con le [Side Quests](../side_quests/index.md)

Infine, vi consigliamo di dare un'occhiata a [**Seqera Platform**](https://seqera.io/), una piattaforma basata su cloud sviluppata dai creatori di Nextflow che rende ancora più facile avviare e gestire i workflow, oltre a gestire i dati ed eseguire analisi in modo interattivo in qualsiasi ambiente.

## Sondaggio di feedback

Prima di procedere, dedicate un minuto per completare il sondaggio del corso! Il vostro feedback ci aiuta a migliorare i nostri materiali formativi per tutti.

[Vai al sondaggio :material-arrow-right:](survey.md){ .md-button .md-button--primary }
