# Per iniziare

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduzione assistita da IA - [scopri di più e suggerisci miglioramenti](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

## Avviare un ambiente di formazione

Per utilizzare l'ambiente preconfigurato che forniamo su GitHub Codespaces, clicca sul pulsante "Open in GitHub Codespaces" qui sotto. Per altre opzioni, consulta [Opzioni di ambiente](../envsetup/index.md).

Consigliamo di aprire l'ambiente di formazione in una nuova scheda o finestra del browser (usa il tasto destro del mouse, ctrl-click o cmd-click a seconda del dispositivo) in modo da poter continuare a leggere mentre l'ambiente si carica.
Sarà necessario tenere queste istruzioni aperte in parallelo per seguire il corso.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

### Nozioni di base sull'ambiente

Questo ambiente di formazione contiene tutto il software, il codice e i dati necessari per seguire il corso, quindi non è necessario installare nulla.

Il codespace è configurato con un'interfaccia VSCode, che include un esplora file, un editor di codice e un terminale.
Tutte le istruzioni fornite durante il corso (ad es. 'apri il file', 'modifica il codice' o 'esegui questo comando') si riferiscono a queste tre parti dell'interfaccia VSCode, salvo diversa indicazione.

Se stai seguendo questo corso in autonomia, ti invitiamo a familiarizzare con le [nozioni di base sull'ambiente](../envsetup/01_setup.md) per ulteriori dettagli.

## Prepararsi al lavoro

Una volta avviato il codespace, ci sono due cose da fare prima di iniziare: impostare la directory di lavoro e dare un'occhiata ai materiali forniti.

### Impostare la directory di lavoro

Per impostazione predefinita, il codespace si apre alla radice di tutti i corsi di formazione.
Per questo corso, spostatevi nella directory `seqera-scale/`:

```bash
cd seqera-scale/
```

Poi impostate VSCode in modo che si concentri su questa directory, così solo i file pertinenti appariranno nella barra laterale dell'esplora file:

```bash
code .
```

!!! tip "Suggerimento"

    Se per qualsiasi motivo vi spostate da questa directory (ad es. il codespace va in sospensione), potete sempre usare il percorso completo per tornarci, supponendo che stiate eseguendo il corso nell'ambiente di formazione GitHub Codespaces:

    ```bash
    cd /workspaces/training/seqera-scale
    ```

### Esplorare i materiali forniti

Potete esplorare i materiali del corso usando l'esplora file sulla sinistra, oppure con il comando `tree`.
Eseguite il seguente comando dal terminale per vedere la struttura completa:

```bash
tree -a .
```

??? abstract "Contenuto della directory"

    ```console
    .
    └── .seqera_config
    ```

Il file **`.seqera_config`** è uno stub che compilerete durante la sezione 3 per configurare la CLI `tw` con il vostro token di accesso Seqera e il workspace.

## Lista di controllo per la preparazione

Pensate di essere pronti a iniziare?

- [ ] Comprendo l'obiettivo di questo corso e i suoi prerequisiti
- [ ] Il mio ambiente è attivo e funzionante
- [ ] Ho impostato correttamente la mia directory di lavoro

Se riuscite a spuntare tutte le caselle, siete pronti per partire.

**Per continuare con [Parte 1: Lanciare pipeline dall'interfaccia web](./01_run_with_seqera.md), clicca sulla freccia nell'angolo in basso a destra di questa pagina.**
