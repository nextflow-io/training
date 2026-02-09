# Riepilogo del corso

Congratulazioni per aver completato il corso di formazione Hello Nextflow! 🎉

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/wnmUAfDL01E?si=y8lAedhEHWaTV4zd&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } Guardate l'[intera playlist sul canale YouTube di Nextflow](https://youtube.com/playlist?list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&si=eF7cLR62goy-lc6n).

:green_book: Potete leggere la [trascrizione del video](./transcripts/07_next_steps.md) insieme al video.
///

## Il vostro percorso

Siete partiti da un flusso di lavoro molto semplice che eseguiva un comando hardcoded.
Nel corso di sei parti, avete trasformato quel flusso di lavoro di base in una pipeline modulare multi-step che utilizza funzionalità chiave di Nextflow tra cui canali, operatori, supporto integrato per i container e opzioni di configurazione.

### Cosa avete costruito

- La forma finale del flusso di lavoro Hello prende come input un file CSV contenente saluti testuali.
- I quattro step sono implementati come processi Nextflow (`sayHello`, `convertToUpper`, `collectGreetings` e `cowpy`) memorizzati in file di modulo separati.
- I risultati vengono pubblicati in una directory chiamata `results/`.
- L'output finale della pipeline è un file di testo semplice contenente ASCII art di un personaggio che pronuncia i saluti in maiuscolo.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_pipeline_complete.svg"
</figure>

1. **`sayHello`:** Scrive ogni saluto nel proprio file di output (_es._ "Hello-output.txt")
2. **`convertToUpper`:** Converte ogni saluto in maiuscolo (_es._ "HELLO")
3. **`collectGreetings`:** Raccoglie tutti i saluti in maiuscolo in un singolo file batch
4. **`cowpy`:** Genera ASCII art utilizzando lo strumento `cowpy`

La configurazione del flusso di lavoro supporta la fornitura di input e parametri in modo flessibile e riproducibile.

### Competenze acquisite

Attraverso questo corso pratico, avete imparato a:

- Descrivere e utilizzare i componenti principali di Nextflow sufficienti per costruire un semplice flusso di lavoro multi-step
- Descrivere concetti di livello successivo come operatori e fabbriche di canali
- Lanciare un flusso di lavoro Nextflow localmente
- Trovare e interpretare gli output (risultati) e i file di log generati da Nextflow
- Risolvere problemi di base

Ora siete equipaggiati con le conoscenze fondamentali per iniziare a sviluppare le vostre pipeline in Nextflow.

## Prossimi passi per sviluppare le vostre competenze

Ecco i nostri 3 suggerimenti principali su cosa fare dopo:

- Applicate Nextflow a un caso d'uso di analisi scientifica con [Nextflow for Science](../nf4_science/index.md)
- Iniziate con nf-core con [Hello nf-core](../hello_nf-core/index.md)
- Esplorate funzionalità più avanzate di Nextflow con le [Side Quests](../side_quests/index.md)

Infine, vi consigliamo di dare un'occhiata a [**Seqera Platform**](https://seqera.io/), una piattaforma basata su cloud sviluppata dai creatori di Nextflow che rende ancora più facile lanciare e gestire i vostri flussi di lavoro, oltre a gestire i vostri dati ed eseguire analisi in modo interattivo in qualsiasi ambiente.

## Questionario di feedback

Prima di proseguire, prendetevi un minuto per completare il questionario del corso! Il vostro feedback ci aiuta a migliorare i nostri materiali di formazione per tutti.

[Compila il questionario :material-arrow-right:](survey.md){ .md-button .md-button--primary }
