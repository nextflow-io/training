# Parte 2: Lanciare pipeline dalla riga di comando

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduzione assistita da IA - [scopri di più e suggerisci miglioramenti](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Nella [Parte 1](./01_run_with_seqera.md), abbiamo lanciato nf-core/rnaseq dall'interfaccia web di Seqera.
Ora facciamo lo stesso dalla riga di comando usando il CLI `tw`, e aggiungiamo una nuova pipeline al nostro workspace.

---

## 1. Lanciare pipeline dalla riga di comando

Nella vista delle esecuzioni, clicca sulla scheda **Command line**.
Vedrai il comando `nextflow run` esatto che Platform ha costruito e inviato per conto tuo — lo stesso tipo di comando che hai eseguito manualmente nel corso Use nf-core.

Platform non sostituisce Nextflow; lo orchestra.
Tutto ciò che puoi fare tramite l'interfaccia web, puoi farlo anche da un terminale usando il CLI `tw`, lo strumento a riga di comando per interagire con le API di Platform.
Questo è utile per automatizzare i lanci da script o pipeline CI/CD.

Lo faremo ora dallo stesso codespace che hai usato per i corsi precedenti.

### 1.1. Installare il CLI tw

Esegui i seguenti comandi nel terminale del tuo Codespace per scaricare e installare il binario `tw`:

```bash
curl -fsSL https://github.com/seqeralabs/tower-cli/releases/latest/download/tw-linux-x86_64 -o tw
chmod +x tw
sudo mv tw /usr/local/bin/
```

Verifica l'installazione:

```bash
tw --version
```

??? success "Output del comando"

    ```console
    Tower CLI version 0.40.0 (build 26db579)
    ```

Il CLI `tw` è installato e pronto per essere configurato.

### 1.2. Ottenere un token di accesso

Il CLI `tw` si autentica con Seqera usando un token di accesso personale.

1. Nell'interfaccia web di Seqera, clicca sul tuo avatar nell'angolo in alto a destra e seleziona **Your tokens**.
2. Clicca su **Add token**, assegnagli un nome (ad es. `training`), e clicca su **Add**.
3. Copia il valore del token — verrà mostrato una sola volta.
   Se non lo salvi subito da qualche parte, dovrai generarne un altro.

### 1.3. Configurare il CLI

Per comodità, configureremo un file di configurazione contenente il token di accesso
appena generato e l'identificatore del workspace.

Apri il file `.seqera_config` in questa directory nell'editor e imposta le due variabili:

- **`TOWER_ACCESS_TOKEN`**: il token generato nella sezione 1.2
- **`TOWER_WORKSPACE_ID`**: l'ID numerico del tuo workspace (la colonna `ID` in `tw workspaces list`, che eseguirai nella sezione 1.4)

Una volta inseriti i valori, carica la configurazione:

```bash
source .seqera_config
```

Verifica la connessione:

```bash
tw info
```

??? success "Output del comando"

    ```console
        Details
    -------------------------+-----------------------------
     Tower API endpoint      | https://api.cloud.seqera.io
     Tower API version       | 1.150.0
     Tower version           | 26.1.0-cycle54
     CLI version             | 0.30.0 (fde9dec)
     CLI minimum API version | 1.148.0
     Authenticated user      | <your-name>

    System health status
    ---------------------------------------+----
     Remote API server connection check    | OK
     Tower API version check               | OK
     Authentication API credential's token | OK
    ```

Il CLI `tw` è ora autenticato e connesso al tuo account Seqera.
Esegui `source .seqera_config` all'inizio di ogni sessione Codespace per ricaricare la configurazione.

!!! tip "Suggerimento"

    Se il tuo workspace non ha un ambiente di calcolo primario impostato, puoi aggiungere `export TOWER_COMPUTE_ENV=<compute-env-name>` al tuo file di configurazione per impostarne uno predefinito.
    Qualsiasi valore di configurazione può essere sovrascritto dalla riga di comando passando il flag esplicitamente (ad es. `--compute-env other-env`).
    Consulta il [riferimento del CLI tw](https://docs.seqera.io/platform/latest/cli/reference) per l'elenco completo delle opzioni e delle variabili d'ambiente.

### 1.4. Esplorare il workspace dal CLI

Elenca i workspace a cui hai accesso:

```bash
tw workspaces list
```

??? success "Output del comando"

    ```console
    Available workspaces:
    ID              | Name            | Full Name                | Visibility
    --------------- | --------------- | ------------------------ | ----------
    <workspace-id>  | my-workspace    | my-org/my-workspace      | PRIVATE
    ```

Visualizza le esecuzioni nel tuo workspace, inclusa quella di nf-core/rnaseq che hai appena lanciato:

```bash
tw runs list
```

??? success "Output del comando"

    ```console
    Pipeline runs for my-org/my-workspace workspace:
    ID        | Status   | Name          | Pipeline               | Run name
    --------- | -------- | ------------- | ---------------------- | --------
    <run-id>  | RUNNING  | ...           | nf-core/rnaseq         | happy_curie
    ```

La stessa esecuzione che stai monitorando nell'interfaccia web è visibile qui.

!!! note "Nota"

    Poiché `TOWER_WORKSPACE_ID` è impostato in `.seqera_config`, puoi omettere `--workspace` da tutti i comandi `tw`.
    Senza la configurazione, dovresti passarlo esplicitamente:

    ```bash
    tw runs list --workspace <org>/<workspace>
    ```

Tutto ciò che è visibile nell'interfaccia web è accessibile dal CLI.

### 1.5. Lanciare nf-core/rnaseq dal CLI

La pipeline che hai aggiunto al tuo workspace nella [Parte 1](./01_run_with_seqera.md) è disponibile per nome nel CLI.
Lanciala con il profilo `test`:

```bash
tw launch nf-core-rnaseq --profile test
```

??? success "Output del comando"

    ```console
    Launching pipeline nf-core-rnaseq
    Run name: focused_einstein
    https://cloud.seqera.io/orgs/my-org/workspaces/my-workspace/watch/<run-id>
    ```

Apri il link nel browser e conferma che l'esecuzione appare nel pannello **Runs**.

Una volta che riesci a vederla in esecuzione, hai confermato che il CLI e l'interfaccia web sono due viste sullo stesso workspace.

!!! note "Nota"

    Puoi anche passare un URL GitHub completo direttamente a `tw launch` senza aggiungere prima la pipeline a un workspace.
    Tuttavia, aggiungere la pipeline esplicitamente prima di lanciarla è generalmente preferibile: salva la configurazione della pipeline per le esecuzioni future, la rende disponibile per nome e la rende visibile a tutti i membri del workspace nel Launchpad.

    È possibile aggiungere una pipeline a un workspace direttamente dalla riga di comando usando `tw`.
    La sezione successiva mostra come farlo con la pipeline nf-core/demo.

### Takeaway

Sai come autenticare il CLI `tw`, ispezionare il tuo workspace e lanciare una pipeline salvata dal terminale.

### Cosa c'è dopo?

Aggiungere una nuova pipeline al workspace dalla riga di comando e lanciarla.

---

## 2. Aggiungere una nuova pipeline ed eseguirla

Qualsiasi pipeline Nextflow su GitHub può essere aggiunta al tuo workspace con `tw pipelines add`, purché abbia un punto di ingresso `main.nf` e un `nextflow.config` nella sua radice.
nf-core/demo è un buon esempio con cui esercitarsi: l'hai già eseguita nel corso Use nf-core, quindi sai cosa fa e cosa aspettarti.

### 2.1. Aggiungere nf-core/demo al workspace

Esegui il seguente comando per registrare la pipeline nel tuo workspace:

```bash
tw pipelines add \
    --name nf-core-demo \
    https://github.com/nf-core/demo
```

??? success "Output del comando"

    ```console
    New pipeline 'nf-core-demo' added at my-org/my-workspace workspace
    ```

La pipeline è ora registrata e apparirà nel Launchpad.

### 2.2. Verificare che appaia nel Launchpad

Elenca le pipeline nel tuo workspace per confermare che sia stata aggiunta:

```bash
tw pipelines list
```

??? success "Output del comando"

    ```console
    Pipelines at my-org/my-workspace:
    ID   | Name            | Repository
    ---- | --------------- | -----------------------------------------
    ...  | nf-core-demo    | https://github.com/nf-core/demo
    ...  | nf-core-rnaseq  | ...
    ```

Apri il tuo workspace nel browser e clicca su **Launchpad** per confermare che nf-core/demo appaia ora accanto a nf-core/rnaseq.

!!! tip "Suggerimento"

    Puoi anche aggiungere pipeline tramite l'interfaccia web: nella barra laterale sinistra, clicca su **Launchpad**, poi su **Add pipeline**, e compila il modulo di conseguenza.

Clicca sul pulsante **Launch** nella voce nf-core/demo per aprire il suo modulo di lancio.
Vedrai che i parametri `input` e `outdir` sono evidenziati in rosso — sono campi obbligatori senza valori predefiniti, perché `tw pipelines add` registra solo il sorgente della pipeline senza pre-configurare alcun parametro.
Le prossime due sezioni illustrano come fornire questi valori: prima tramite il modulo web, poi dalla riga di comando.

### 2.3. Lanciare nf-core/demo dall'interfaccia web

Con il modulo di lancio aperto, compila i due parametri obbligatori.

Per `input`, inserisci l'URL del samplesheet di test dal profilo test di nf-core/demo.
Puoi trovarlo in `conf/test.config` all'interno del repository della pipeline, che hai esaminato nel corso Use nf-core:

```
https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv
```

Per `outdir`, inserisci un percorso di archiviazione cloud dove la pipeline può scrivere i suoi risultati.
Usa il bucket configurato per il tuo workspace, con una sottodirectory per mantenere le esecuzioni organizzate:

```
s3://my-bucket/demo-results
```

Una volta compilati entrambi i campi, clicca sul pulsante blu **Launch**.

L'esecuzione appare nel pannello **Runs** e dovrebbe completarsi in pochi minuti sul dataset di test.
Clicca sull'esecuzione per esplorare la tabella delle attività e gli eventuali report di esecuzione.

### 2.4. Lanciare nf-core/demo dal CLI

A differenza di `nextflow run`, il comando `tw launch` non accetta flag di parametri individuali come `--input` o `--outdir`.
I parametri devono essere forniti tramite un file in formato YAML o JSON, passato con `--params-file`.
Questo favorisce la riproducibilità: un file di parametri salvato documenta esattamente quali valori sono stati usati per un'esecuzione, rendendo facile ripetere o condividere una configurazione di esecuzione.

Crea un file di parametri nella tua directory di lavoro:

```bash
touch params.yaml
```

Aprilo nell'editor e aggiungi il percorso di output:

```yaml title="params.yaml"
outdir: "s3://my-bucket/demo-results-cli"
```

Ora puoi lanciare la pipeline usando il profilo `test` (che fornisce il samplesheet `input`) e il file di parametri (che fornisce `outdir`):

```bash
tw launch nf-core-demo --profile test --params-file params.yaml
```

??? success "Output del comando"

    ```console
    Launching pipeline nf-core-demo
    Run name: focused_feynman
    https://cloud.seqera.io/orgs/my-org/workspaces/my-workspace/watch/<run-id>
    ```

Apri il link per confermare che l'esecuzione appaia nel pannello **Runs**.

!!! tip "Suggerimento"

    Puoi includere il file di parametri durante la fase di configurazione iniziale se desideri impostare alcuni valori predefiniti, insieme ad alcune proprietà aggiuntive per riprodurre ciò che abbiamo fatto in precedenza tramite il modulo web:

    ```bash
    tw pipelines add \
      --name nf-core-demo-gg \
      --description "Demo pipeline with defaults for testing" \
      --params-file params.yaml \
      --profile test \
      https://github.com/nf-core/demo
    ```

### Takeaway

Sai come aggiungere qualsiasi pipeline Nextflow ospitata su GitHub al tuo workspace e lanciarla, sia dall'interfaccia web compilando i parametri manualmente, sia dal CLI `tw` combinando un profilo con un file di parametri.

---

## Riepilogo

In questa parte hai imparato a:

- Autenticare il CLI `tw` e lanciare una pipeline salvata dal terminale
- Aggiungere una nuova pipeline da GitHub usando il CLI e verificare che appaia nel Launchpad
- Lanciare una pipeline dall'interfaccia web di Seqera compilando manualmente i parametri obbligatori
- Lanciare una pipeline dal CLI usando un profilo Nextflow e un file di parametri
