# Parte 2: Configurare la pipeline

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduzione assistita da IA - [scopri di più e suggerisci miglioramenti](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Nella [Parte 1](./01_run_nextflow.md), abbiamo eseguito una pipeline completa a più fasi che elabora più input in parallelo usando i container.
Ora vedremo come configurare il comportamento della pipeline usando `nextflow.config`: prima esaminando il file di configurazione che abbiamo già fornito, poi esplorando un paio di altri modi per fornire la configurazione, e infine controllando come e dove vengono pubblicati gli output.

---

## 1. Esaminare il file di configurazione principale

Nextflow rileva automaticamente `nextflow.config` dalla directory di lavoro e applica le sue impostazioni a ogni esecuzione.

Vi forniamo un file di configurazione che copre quattro aree: packaging del software, impostazioni dei processi, parametri della pipeline e profili di esecuzione.

??? full-code "nextflow.config"

    ```groovy title="nextflow.config" linenums="1"
    /*
     * Packaging del software
     */
    docker.enabled = true

    /*
     * Impostazioni dei processi
     */
    process {
        cpus = 1
        memory = 1.GB
    }

    /*
     * Parametri della pipeline
     */
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }

    /*
     * Profili
     */
    profiles {
        test {
            params.input = 'data/greetings.csv'
            params.batch = 'test'
            params.character = 'tux'
        }
        conda {
            docker.enabled = false
            conda.enabled = true
        }
    }
    ```

Esaminiamo ciascuna sezione, poi mettiamo i profili in pratica eseguendo la pipeline con uno di essi.

!!! note "Nota"

    Questa configurazione copre l'esecuzione locale su una singola macchina.
    Nextflow supporta anche gli scheduler HPC (SLURM, PBS, LSF) e gli executor cloud (AWS Batch, Google Cloud Batch, Azure Batch), tutti configurati attraverso lo stesso meccanismo `nextflow.config`.
    Consultate la [Parte 1: Adattarsi all'ambiente di calcolo](../execution_config/01_packaging_and_execution.md) nel corso [Execution Config](../execution_config/index.md) per una guida completa di queste opzioni.

### 1.1. Packaging del software

Il packaging del software è il modo in cui Nextflow fornisce gli strumenti effettivi di cui i vostri processi hanno bisogno, che si tratti di un'immagine container, un ambiente Conda, o altro.

```groovy title="nextflow.config" linenums="1"
/*
 * Packaging del software
 */
docker.enabled = true
```

Questa riga abilita Docker per ogni processo.
Qualsiasi processo che dichiara una direttiva `container` viene eseguito all'interno dell'immagine specificata.

### 1.2. Impostazioni dei processi

Ricordate che un processo è un singolo passaggio nella vostra pipeline, come `sayHello` o `cowpy`.
Nextflow vi permette di configurare diversi aspetti di come ciascuno viene effettivamente eseguito: quanta CPU e memoria ottiene, quale container o ambiente Conda utilizza, e altro ancora.

```groovy title="nextflow.config" linenums="6"
/*
 * Impostazioni dei processi
 */
process {
    cpus = 1
    memory = 1.GB
}
```

Questo limita ogni processo a una singola CPU e 1 GB di memoria.

Nextflow vi permette anche di impostare valori diversi per singoli processi con nome o gruppi di processi; imparerete come farlo nella [Parte 2: Gestire le risorse di calcolo e i fallimenti](../execution_config/02_resources_and_retries.md#12-set-resource-allocations-for-a-specific-process) del corso [Execution Config](../execution_config/index.md).

### 1.3. Parametri della pipeline

I parametri sono gli input da riga di comando della pipeline, gli stessi flag `--input`, `--batch` e `--character` che avete già impostato direttamente dalla riga di comando.
Impostare qui i valori predefiniti significa che non dovete digitarli ogni volta, anche se, come vedrete più avanti in questa parte, ci sono un paio di altri modi per fornirli.

```groovy title="nextflow.config" linenums="14"
/*
 * Parametri della pipeline
 */
params {
    input = 'data/greetings.csv'
    batch = 'batch'
    character = 'turkey'
}
```

Questi valori predefiniti entrano in gioco ogni volta che un parametro non viene fornito dalla riga di comando, quindi eseguire `nextflow run main.nf` senza flag funziona comunque.

### 1.4. Profili

I profili vi permettono di raggruppare un insieme di impostazioni sotto un unico nome, così potete passare da una configurazione completa all'altra con un solo flag invece di modificare i valori manualmente ogni volta.

```groovy title="nextflow.config" linenums="23"
/*
 * Profili
 */
profiles {
    test {
        params.input = 'data/greetings.csv'
        params.batch = 'test'
        params.character = 'tux'
    }
    conda {
        docker.enabled = false
        conda.enabled = true
    }
}
```

Il profilo `test` sovrascrive tre parametri per eseguire la pipeline con un insieme di input piccolo e ben definito; ogni pipeline nf-core ne include uno per la validazione rapida, ed è una convenzione che vale la pena seguire anche nelle vostre pipeline.

Il profilo `conda` cambia il packaging del software da Docker a Conda.

Si attiva un profilo passando `-profile <nome>` dalla riga di comando.

Mettiamo in pratica il profilo `test`.

```bash
nextflow run main.nf -profile test
```

??? success "Output del comando"

    ```console hl_lines="6 7 8 9"
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [reverent_heisenberg] revision: ce74f81996

    executor >  local (8)
    [3d/8a12c7] sayHello (3)       | 3 of 3 ✔
    [e7/e0934f] convertToUpper (2) | 3 of 3 ✔
    [1d/616569] collectGreetings   | 1 of 1 ✔
    [44/7d46cf] cowpy              | 1 of 1 ✔

    Outputs:

      /workspaces/training/nextflow-run/results

      first_output:
        - test/intermediates/Bonjour-output.txt
        - test/intermediates/Hello-output.txt
        - test/intermediates/Hola-output.txt

      uppercased:
        - test/intermediates/UPPER-Bonjour-output.txt
        - test/intermediates/UPPER-Hello-output.txt
        - test/intermediates/UPPER-Hola-output.txt

      collected: test/intermediates/COLLECTED-test-output.txt

      batch_report: test/test-report.txt

      cowpy_art: test/cowpy-COLLECTED-test-output.txt
    ```

La pipeline viene eseguita con `batch = 'test'` e `character = 'tux'`.
Controllate `results/test/`: il nome del batch fa ora parte del percorso della directory, e l'arte ASCII mostra il pinguino tux invece del tacchino.

!!! note "Nota"

    Potete attivare più profili contemporaneamente, e usare `nextflow config -profile <nome>,<nome>` per vedere il risultato completamente risolto prima di eseguire qualsiasi cosa.
    La combinazione di profili e il modo in cui Nextflow risolve i conflitti tra di essi è trattata in dettaglio nella [Parte 3: Usare i profili per cambiare configurazione](../execution_config/03_profiles.md) del corso [Execution Config](../execution_config/index.md).

### Takeaway

Sapete cosa fanno gli elementi più comuni di un file `nextflow.config` e come attivare un profilo.

### Cosa c'è dopo?

Impareremo un paio di altri modi per fornire valori di configurazione senza modificare il file `nextflow.config` principale, utili per configurare singole esecuzioni e per condividere un insieme preciso di impostazioni con qualcun altro.

---

## 2. Fornire la configurazione tramite file supplementari

Impostare i valori predefiniti in `nextflow.config` funziona bene per i valori che cambiano raramente.
Nextflow vi offre anche due meccanismi più mirati: un file di configurazione specifico per l'esecuzione, per adattare l'esecuzione a un ambiente particolare, e un file di parametri per condividere un insieme preciso di valori di input con un collaboratore.

### 2.1. Usare un file di configurazione specifico per l'esecuzione

Supponiamo che stiate spostando la pipeline su una macchina che non ha Docker e vogliate dare a ogni processo più risorse con cui lavorare.
Create un nuovo file di configurazione con solo le sovrascritture di cui avete bisogno:

```groovy title="custom.config" linenums="1"
process {
    cpus = 2
    memory = 2.GB
}

docker.enabled = false
conda.enabled = true
```

Passatelo insieme alla vostra pipeline principale con `-c`:

```bash
nextflow run main.nf -c custom.config
```

??? success "Output del comando"

    ```console hl_lines="6 7 8 9"
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [exotic_cray] revision: ce74f81996

    executor >  local (8)
    [77/e02315] sayHello (1)       | 3 of 3 ✔
    [a6/ccf44b] convertToUpper (3) | 3 of 3 ✔
    [56/fd1296] collectGreetings   | 1 of 1 ✔
    [2c/205a94] cowpy              | 1 of 1 ✔

    Outputs:

      /workspaces/training/nextflow-run/results

      first_output:
        - batch/intermediates/Bonjour-output.txt
        - batch/intermediates/Hello-output.txt
        - batch/intermediates/Hola-output.txt

      uppercased:
        - batch/intermediates/UPPER-Bonjour-output.txt
        - batch/intermediates/UPPER-Hello-output.txt
        - batch/intermediates/UPPER-Hola-output.txt

      collected: batch/intermediates/COLLECTED-batch-output.txt

      batch_report: batch/batch-report.txt

      cowpy_art: batch/cowpy-COLLECTED-batch-output.txt
    ```

Nextflow unisce `custom.config` sopra il `nextflow.config` della pipeline, quindi ogni processo ora ottiene 2 CPU e 2 GB di memoria invece dei valori predefiniti, e viene eseguito tramite Conda invece di Docker.
`cowpy` è l'unico processo con un pacchetto Conda dichiarato insieme al suo container, quindi è quello per cui vedrete Nextflow costruire effettivamente un ambiente:

```console
Creating env using conda: conda-forge::cowpy==1.1.5 [cache /path/to/work/conda/env-898314d566668b6587ad714ae06b8520]
```

Un file piccolo che sovrascrive solo l'allocazione delle risorse e il packaging, senza toccare i parametri della pipeline, è esattamente il pattern che le pipeline nf-core si aspettano dalle configurazioni istituzionali.
Consultate il repository [nf-core/configs](https://github.com/nf-core/configs) per esempi reali.

Questo vi offre un modo usa-e-getta per adattare una pipeline a un nuovo ambiente senza toccare la configurazione normale.

### 2.2. Usare un file di parametri

Supponiamo invece che dobbiate condividere un insieme preciso di parametri di esecuzione con un collaboratore, o registrarli per una pubblicazione.

Nextflow vi permette di fornire [file di parametri](https://nextflow.io/docs/latest/config.html#parameter-file) in formato YAML o JSON, che sono un modo più semplice per distribuire un insieme esatto e riproducibile di valori.

Un file di parametri chiamato `test-params.yaml` è già fornito nella vostra directory di lavoro:

```yaml title="test-params.yaml" linenums="1"
input: "data/greetings.csv"
batch: "yaml"
character: "stegosaurus"
```

La sintassi usa i due punti (`:`) invece dei segni di uguale (`=`) usati in `nextflow.config`, poiché questo file è YAML semplice anziché Groovy.

!!! info "Info"

    È fornita anche una versione JSON, `test-params.json`. Sentitevi liberi di provarla da soli; la sintassi per passarla è identica.

Passate il file con `-params-file`:

```bash
nextflow run main.nf -params-file test-params.yaml
```

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [sharp_faraday] revision: ce74f81996

    executor >  local (8)
    [1c/9ff63e] sayHello (1)       | 3 of 3 ✔
    [3b/bb5691] convertToUpper (2) | 3 of 3 ✔
    [cd/2c1f6e] collectGreetings   | 1 of 1 ✔
    [89/c333bc] cowpy              | 1 of 1 ✔

    Outputs:

      /workspaces/training/nextflow-run/results

      first_output:
        - yaml/intermediates/Bonjour-output.txt
        - yaml/intermediates/Hello-output.txt
        - yaml/intermediates/Hola-output.txt

      uppercased:
        - yaml/intermediates/UPPER-Bonjour-output.txt
        - yaml/intermediates/UPPER-Hello-output.txt
        - yaml/intermediates/UPPER-Hola-output.txt

      collected: yaml/intermediates/COLLECTED-yaml-output.txt

      batch_report: yaml/yaml-report.txt

      cowpy_art: yaml/cowpy-COLLECTED-yaml-output.txt
    ```

??? abstract "Contenuto del file"

    ```console title="results/yaml/cowpy-COLLECTED-yaml-output.txt"
     _________
    / BONJOUR \
    | HOLA    |
    \ HELLO   /
     ---------
    \                             .       .
     \                           / `.   .' "
      \                  .---.  <    > <    >  .---.
       \                 |    \  \ - ~ ~ - /  /    |
             _____          ..-~             ~-..-~
            |     |   \~~~\.'                    `./~~~/
           ---------   \__/                        \__/
          .'  O    \     /               /       \  "
         (_____,    `._.'               |         }  \/~~~/
          `----.          /       }     |        /    \__/
                `-.      |       /      |       /      `. ,~~|
                    ~-.__|      /_ - ~ ^|      /- _      `..-'
                        |     /        |     /     ~-.     `-. _  _  _
                        |_____|        |_____|         ~ - . _ _ _ _ _>
    ```

Un file di parametri è particolarmente utile quando una pipeline ha più di una manciata di parametri: vi permette di fornirli tutti in una volta, senza una riga di comando lunghissima o alcuna modifica allo script del flusso di lavoro, ed è facile da distribuire insieme ai vostri risultati.

### Takeaway

Conoscete altri due modi per fornire la configurazione: un file di configurazione specifico per l'esecuzione per adattare l'esecuzione a un nuovo ambiente, e un file di parametri per condividere valori di input esatti e riproducibili.

### Cosa c'è dopo?

Impareremo come controllare come e dove vengono pubblicati gli output della pipeline.

---

## 3. Gestire gli output della pipeline

L'autore di una pipeline decide come gli output sono organizzati nel codice, ma non è necessario toccare quel codice per controllare dove finiscono o come ci arrivano.
Nextflow vi offre modi a livello di configurazione per farlo: impostare una directory di output di base e scegliere se i file vengono copiati o collegati tramite symlink.

### 3.1. Personalizzare la directory di output

Per impostazione predefinita, Nextflow pubblica gli output sotto `results/`.
Puntate altrove con `-output-dir` (o la sua forma abbreviata, `-o`):

```bash
nextflow run main.nf -output-dir outputs
```

??? success "Output del comando"

    ```console hl_lines="6 7 8 9"
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [serene_kimura] revision: ce74f81996

    executor >  local (8)
    [31/df5c15] sayHello (3)       | 3 of 3 ✔
    [08/8bb2e5] convertToUpper (1) | 3 of 3 ✔
    [e5/5814da] collectGreetings   | 1 of 1 ✔
    [cf/8ab8c6] cowpy              | 1 of 1 ✔

    Outputs:

      /workspaces/training/nextflow-run/outputs

      first_output:
        - batch/intermediates/Bonjour-output.txt
        - batch/intermediates/Hello-output.txt
        - batch/intermediates/Hola-output.txt

      uppercased:
        - batch/intermediates/UPPER-Bonjour-output.txt
        - batch/intermediates/UPPER-Hello-output.txt
        - batch/intermediates/UPPER-Hola-output.txt

      collected: batch/intermediates/COLLECTED-batch-output.txt

      batch_report: batch/batch-report.txt

      cowpy_art: batch/cowpy-COLLECTED-batch-output.txt
    ```

??? abstract "Contenuto della directory"

    ```console
    outputs/batch
    ├── batch-report.txt
    ├── cowpy-COLLECTED-batch-output.txt
    └── intermediates
        ├── Bonjour-output.txt
        ├── COLLECTED-batch-output.txt
        ├── Hello-output.txt
        ├── Hola-output.txt
        ├── UPPER-Bonjour-output.txt
        ├── UPPER-Hello-output.txt
        └── UPPER-Hola-output.txt
    ```

Gli output ora finiscono sotto `outputs/batch/` invece del valore predefinito `results/batch/`.
Il codice della pipeline decide ancora la struttura all'interno di quella directory di base, come le sottodirectory `batch/` e `intermediates/`; `-output-dir` controlla solo dove inizia quella struttura.

`-output-dir` è in realtà solo una scorciatoia da riga di comando per l'opzione di configurazione `outputDir`, quindi può essere inserita ovunque sia possibile inserire la configurazione: direttamente in `nextflow.config`, all'interno di un profilo, o in un file overlay `-c` come quello usato in precedenza in questa parte.
Ad esempio, questo frammento mostra la stessa impostazione inserita direttamente in `nextflow.config` invece di essere passata dalla riga di comando:

```groovy title="nextflow.config"
outputDir = 'outputs'
```

Consultate [Configuration file](https://nextflow.io/docs/latest/config.html) nella documentazione di riferimento di Nextflow per l'elenco completo dei posti in cui un'opzione di configurazione come questa può essere inserita.

### 3.2. Scegliere come vengono pubblicati gli output

Per impostazione predefinita, Nextflow pubblica gli output come symlink che puntano alle posizioni degli output sotto `work/`, non come copie reali:

```console
$ ls -l results/batch/intermediates/Hello-output.txt
lrwxr-xr-x  ...  Hello-output.txt -> /workspaces/training/nextflow-run/work/b7/b8c4d1.../Hello-output.txt
```

Gli autori della pipeline possono impostare la 'modalità di pubblicazione' su `'copy'` o `'move'` per ogni singolo processo nel codice del flusso di lavoro.
Di solito lo fanno per gli output finali della pipeline, lasciando il comportamento predefinito `'symlink'` per i file intermedi che possono essere eliminati una volta completata l'intera pipeline.

Questo evita di duplicare i dati su disco, ma significa che non potete eliminare le directory delle attività sotto `work/` senza rompere il collegamento, perdendo la possibilità di usare `-resume`.
Se volete che tutti i file di output vengano copiati correttamente, impostate [`workflow.output.mode`](https://nextflow.io/docs/latest/reference/config.html#workflow) su `'copy'` nella configurazione della pipeline. (A differenza di `-output-dir`, non esiste un flag da riga di comando per questo; è solo configurazione.)

Provate a impostarlo in `nextflow.config`:

=== "Dopo"

    ```groovy title="nextflow.config" hl_lines="7"
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }

    workflow.output.mode = 'copy'
    ```

=== "Prima"

    ```groovy title="nextflow.config"
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }
    ```

Poi eseguite la pipeline, cambiando il nome del batch in modo da poter vedere la differenza negli output:

```bash
nextflow run main.nf --batch withmode
```

??? success "Output del comando"

    ```console hl_lines="6 7 8 9"
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [angry_noether] revision: ce74f81996

    executor >  local (8)
    [41/2b478d] sayHello (3)       | 3 of 3 ✔
    [bf/dd2840] convertToUpper (2) | 3 of 3 ✔
    [ea/364e97] collectGreetings   | 1 of 1 ✔
    [10/76fe7b] cowpy              | 1 of 1 ✔

    Outputs:

      /workspaces/training/nextflow-run/results

      first_output:
        - withmode/intermediates/Bonjour-output.txt
        - withmode/intermediates/Hello-output.txt
        - withmode/intermediates/Hola-output.txt

      uppercased:
        - withmode/intermediates/UPPER-Bonjour-output.txt
        - withmode/intermediates/UPPER-Hello-output.txt
        - withmode/intermediates/UPPER-Hola-output.txt

      collected: withmode/intermediates/COLLECTED-withmode-output.txt

      batch_report: withmode/withmode-report.txt

      cowpy_art: withmode/cowpy-COLLECTED-withmode-output.txt
    ```

Date un'occhiata a uno dei file di output come prima:

```console
$ ls -l results/withmode/intermediates/Hello-output.txt
-rw-r--r--  ...  Hello-output.txt
```

Ora è un file reale e indipendente che rimarrà disponibile anche se `work/` viene ripulita.

!!! warning "Avviso"

    L'impostazione `workflow.output.mode` riempie solo un valore predefinito per gli output che non hanno già una modalità impostata nel codice della pipeline.
    Non può sovrascrivere una modalità che l'autore ha codificato direttamente, indipendentemente da ciò che impostate.

### Takeaway

Sapete come personalizzare la directory di output di base e scegliere tra output copiati e collegati tramite symlink, il tutto senza toccare il codice della pipeline.

### Cosa c'è dopo?

Passate alla [Parte 3](./03_manage_executions.md), dove imparerete come ispezionare la cronologia delle esecuzioni passate, generare report di esecuzione e ripulire le vecchie directory di lavoro.

---

## Riepilogo

In questa parte avete imparato a:

- Configurare il comportamento della pipeline usando `nextflow.config` e i profili
- Fornire la configurazione tramite un file di configurazione specifico per l'esecuzione o un file di parametri
- Personalizzare la directory di output e scegliere tra output copiati e collegati tramite symlink
