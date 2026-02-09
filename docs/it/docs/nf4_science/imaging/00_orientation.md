# Orientamento

Questo orientamento presuppone che abbiate già aperto l'ambiente di formazione cliccando sul pulsante "Open in GitHub Codespaces".
Se non l'avete ancora fatto, fatelo ora, idealmente in una seconda finestra o scheda del browser in modo da poter consultare queste istruzioni.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://github.com/codespaces/new?hide_repo_select=true&ref=master&repo=290790519&skip_quickstart=true&machine=premiumLinux&devcontainer_path=.devcontainer%2Fdevcontainer.json)

!!!warning "Requisito dimensione macchina"

    Assicuratevi di selezionare una **macchina a 8 core** quando create il vostro Codespace per questo corso di formazione. I flussi di lavoro di bioimaging richiedono risorse di calcolo aggiuntive.

## GitHub Codespaces

L'ambiente GitHub Codespaces contiene tutto il software, il codice e i dati necessari per seguire questo corso di formazione, quindi non dovete installare nulla voi stessi.
Tuttavia, avete bisogno di un account GitHub (gratuito) per accedere, e se non avete familiarità con l'interfaccia dovreste prendervi qualche minuto per familiarizzare completandone il mini-corso [GitHub Codespaces Orientation](../../envsetup/index.md).

## Pre-scaricare le immagini Docker

Una volta aperto il vostro Codespace, scarichiamo in anticipo tutte le immagini Docker di cui avremo bisogno per questo corso di formazione.
Questo farà risparmiare tempo in seguito e garantirà un'esecuzione fluida dei flussi di lavoro.

Aprite una nuova scheda del terminale ed eseguite il seguente comando:

```bash
nextflow run nf-core/molkart -profile docker,test -stub -resume --outdir results
```

Questo comando scaricherà tutte le immagini Docker necessarie in background.
Potete continuare con il resto dell'orientamento mentre questo è in esecuzione.

!!!tip

    Il flag `-stub` consente alla pipeline di essere eseguita rapidamente senza elaborare dati reali, il che è perfetto per scaricare le immagini. Potete monitorare il progresso nella scheda del terminale.

## Directory di lavoro

Durante questo corso di formazione, lavoreremo nella directory `nf4-science/imaging/`.

Cambiate directory ora eseguendo questo comando nel terminale:

```bash
cd nf4-science/imaging/
```

!!!tip

    Se per qualsiasi motivo vi spostate fuori da questa directory, potete sempre usare il percorso completo per ritornarvi, assumendo che stiate eseguendo questo nell'ambiente di formazione GitHub Codespaces:

    ```bash
    cd /workspaces/training/nf4-science/imaging
    ```

**Ora, per iniziare il corso, cliccate sulla freccia nell'angolo in basso a destra di questa pagina.**
