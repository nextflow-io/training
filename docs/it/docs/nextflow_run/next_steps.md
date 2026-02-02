# Riepilogo del corso

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduzione assistita da IA - [scopri di più e suggerisci miglioramenti](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Congratulazioni per aver completato il corso di formazione Nextflow Run!

<!-- placeholder for video -->

## Il tuo percorso

Hai iniziato con un workflow molto basilare, e hai imparato a eseguirlo, trovare gli output e gestire la sua esecuzione.
Poi, hai lavorato attraverso versioni sempre più complesse di quel workflow e hai imparato a riconoscere i concetti e i meccanismi essenziali che alimentano le pipeline Nextflow, inclusi channel e operatori, modularizzazione del codice e container.
Infine, hai imparato come personalizzare la configurazione di una pipeline per adattarla alle tue preferenze e alla tua infrastruttura di calcolo.

### Cosa hai imparato

Ora sei in grado di gestire l'esecuzione della pipeline Hello, descrivere come è strutturata e identificare i principali pezzi di codice coinvolti.

- La forma finale del workflow Hello prende come input un file CSV contenente saluti di testo.
- I quattro step sono implementati come process Nextflow (`sayHello`, `convertToUpper`, `collectGreetings` e `cowpy`) memorizzati in file modulo separati.
- I risultati vengono pubblicati in una directory chiamata `results/`.
- L'output finale della pipeline è un file di testo semplice contenente arte ASCII di un personaggio che pronuncia i saluti in maiuscolo.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_pipeline_complete.svg"
</figure>

1. **`sayHello`:** Scrive ogni saluto nel suo file di output (es. "Hello-output.txt")
2. **`convertToUpper`:** Converte ogni saluto in maiuscolo (es. "HELLO")
3. **`collectGreetings`:** Raccoglie tutti i saluti maiuscoli in un singolo file batch
4. **`cowpy`:** Genera arte ASCII usando lo strumento `cowpy`

La configurazione del workflow supporta la fornitura di input e parametri in modo flessibile e riproducibile.

### Competenze acquisite

Attraverso questo corso pratico, hai imparato come:

- Lanciare un workflow Nextflow localmente
- Trovare e interpretare output (risultati) e file di log generati da Nextflow
- Riconoscere i componenti principali di Nextflow che costituiscono un semplice workflow multi-step
- Descrivere concetti avanzati come operatori e channel factory
- Configurare pipeline per diversi ambienti di calcolo

Ora sei equipaggiato con le conoscenze fondamentali per iniziare a integrare le pipeline Nextflow esistenti nel tuo lavoro.

## Prossimi passi per sviluppare le tue competenze

Ecco i nostri migliori suggerimenti su cosa fare dopo:

- Non limitarti a eseguire Nextflow, scrivilo! Diventa uno sviluppatore Nextflow con [Hello Nextflow](../hello_nextflow/index.md)
- Applica Nextflow a un caso d'uso di analisi scientifica con [Nextflow for Science](../nf4_science/index.md)
- Inizia con nf-core con [Hello nf-core](../hello_nf-core/index.md)
- Impara tecniche di troubleshooting con la [Debugging Side Quest](../side_quests/debugging.md)

Infine, ti raccomandiamo di dare un'occhiata a [**Seqera Platform**](https://seqera.io/), una piattaforma cloud-based sviluppata dai creatori di Nextflow che rende ancora più facile lanciare e gestire i tuoi workflow, oltre a gestire i tuoi dati e eseguire analisi interattivamente in qualsiasi ambiente.

## Ottenere aiuto

Per risorse di aiuto e supporto della community, consulta la [pagina di Aiuto](../help.md).

## Sondaggio di feedback

Prima di proseguire, prenditi un minuto per completare il sondaggio del corso! Il tuo feedback ci aiuta a migliorare i nostri materiali di formazione per tutti.

[Compila il sondaggio :material-arrow-right:](survey.md){ .md-button .md-button--primary }
