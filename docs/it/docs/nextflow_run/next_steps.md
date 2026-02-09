# Riepilogo del corso

Congratulazioni per aver completato il corso di formazione Nextflow Run! 🎉

<!-- placeholder for video -->

## Il vostro percorso

Siete partiti da un flusso di lavoro molto semplice e avete imparato a eseguirlo, trovare gli output e gestirne l'esecuzione.
Poi, avete lavorato attraverso versioni sempre più complesse di quel flusso di lavoro e avete imparato a riconoscere i concetti e i meccanismi essenziali che alimentano le pipeline Nextflow, inclusi canali e operatori, modularizzazione del codice e container.
Infine, avete imparato come personalizzare la configurazione di una pipeline per adattarla alle vostre preferenze e alla vostra infrastruttura computazionale.

### Cosa avete imparato

Ora siete in grado di gestire l'esecuzione della pipeline Hello, descrivere come è strutturata e identificare i principali pezzi di codice coinvolti.

- La forma finale del flusso di lavoro Hello prende come input un file CSV contenente saluti testuali.
- I quattro passaggi sono implementati come processi Nextflow (`sayHello`, `convertToUpper`, `collectGreetings` e `cowpy`) memorizzati in file modulo separati.
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

Attraverso questo corso pratico, avete imparato come:

- Lanciare un flusso di lavoro Nextflow localmente
- Trovare e interpretare gli output (risultati) e i file di log generati da Nextflow
- Riconoscere i componenti core di Nextflow che costituiscono un semplice flusso di lavoro multi-step
- Descrivere concetti avanzati come operatori e fabbriche di canali
- Configurare pipeline per diversi ambienti computazionali

Ora siete equipaggiati con le conoscenze fondamentali per iniziare a integrare pipeline Nextflow esistenti nel vostro lavoro.

## Prossimi passi per sviluppare le vostre competenze

Ecco i nostri migliori suggerimenti su cosa fare dopo:

- Non limitatevi a eseguire Nextflow, scrivetelo! Diventate sviluppatori Nextflow con [Hello Nextflow](../hello_nextflow/index.md)
- Applicate Nextflow a un caso d'uso di analisi scientifica con [Nextflow for Science](../nf4_science/index.md)
- Iniziate con nf-core con [Hello nf-core](../hello_nf-core/index.md)
- Imparate tecniche di risoluzione dei problemi con la [Debugging Side Quest](../side_quests/debugging.md)

Infine, vi consigliamo di dare un'occhiata a [**Seqera Platform**](https://seqera.io/), una piattaforma basata su cloud sviluppata dai creatori di Nextflow che rende ancora più facile lanciare e gestire i vostri flussi di lavoro, oltre a gestire i vostri dati ed eseguire analisi in modo interattivo in qualsiasi ambiente.

## Ottenere aiuto

Per risorse di aiuto e supporto della comunità, consultate la [pagina di Aiuto](../help.md).

## Sondaggio di feedback

Prima di proseguire, dedicate un minuto a completare il sondaggio del corso! Il vostro feedback ci aiuta a migliorare i nostri materiali di formazione per tutti.

[Partecipa al sondaggio :material-arrow-right:](survey.md){ .md-button .md-button--primary }
