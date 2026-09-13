# Parte 2: Gestire le risorse di calcolo e i fallimenti

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduzione assistita da IA - [scopri di più e suggerisci miglioramenti](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>


Nella [Parte 1](./01_packaging_and_execution.md), hai adattato dove e come vengono eseguite le attività di una pipeline.
Qui adatteremo quante risorse di calcolo riceve ogni attività e cosa succede quando un'attività fallisce nonostante la nostra migliore stima sull'allocazione.

---

## 1. Controllare le allocazioni delle risorse di calcolo

Per impostazione predefinita, Nextflow alloca una singola CPU a ogni processo tramite la direttiva `cpus`, e non impone un limite di memoria a meno che non ne venga impostato uno:

```groovy title="Built-in configuration"
process {
    cpus = 1
}
```

Sai già da [Nextflow Run](../nextflow_run/index.md) che la configurazione di questa pipeline imposta `memory` a 1 GB per tutti i processi.
Ma come fai a sapere quali valori usare effettivamente per le tue pipeline?

### 1.1. Generare un report sull'utilizzo delle risorse

Hai già generato un report di esecuzione con `-with-report` in [Nextflow Run](../nextflow_run/02_configure_pipeline.md).
Quello stesso report è il modo in cui scopri quanta CPU e memoria i tuoi processi hanno effettivamente bisogno: esegui il flusso di lavoro con alcune allocazioni predefinite, registra l'utilizzo effettivo, poi regola di conseguenza.

```bash
nextflow run main.nf -with-report report-config-1.html
```

Il report è un file HTML che puoi aprire in un browser.
Suddivide il tempo di esecuzione e l'utilizzo delle risorse per processo, inclusa la percentuale delle risorse allocate effettivamente utilizzata.
Ecco cosa mostra per `cowpy` con i valori predefiniti attuali (1 CPU, 1 GB di memoria):

| Metrica               | Valore |
| --------------------- | ------ |
| Utilizzo CPU          | 116%   |
| Memoria massima usata | 6.4 MB |
| Memoria allocata      | 1 GB   |

`cowpy` usa ben meno dell'1% della sua allocazione di 1 GB; il `%cpu` superiore al 100% significa semplicemente che usa brevemente più di una CPU di elaborazione all'interno del container, in brevi raffiche.

Consulta [Reports](https://nextflow.io/docs/latest/reports.html) per l'elenco completo delle funzionalità disponibili.

### 1.2. Impostare le allocazioni di risorse per un processo specifico

Il report sopra mostra `cowpy` comodamente all'interno della sua allocazione attuale, ma supponiamo che tu voglia dargli più margine comunque, ad esempio perché ti aspetti input più grandi in produzione.
Puoi sovrascrivere i valori predefiniti per un singolo processo con `withName`.

=== "Dopo"

    ```groovy title="nextflow.config" linenums="6" hl_lines="5-6"
    process {
        memory = 1.GB
        withName: 'cowpy' {
            conda = 'conda-forge::cowpy==1.1.5'
            memory = 2.GB
            cpus = 2
        }
    }
    ```

=== "Prima"

    ```groovy title="nextflow.config" linenums="6"
    process {
        memory = 1.GB
        withName: 'cowpy' {
            conda = 'conda-forge::cowpy==1.1.5'
        }
    }
    ```

Con questa configurazione, ogni processo richiede 1 GB di memoria e una singola CPU, tranne `cowpy`, che richiede 2 GB e 2 CPU (oltre all'impostazione `conda` della [Parte 1](./01_packaging_and_execution.md)).

!!! info "Info"

    Se la tua macchina ha poche CPU e ne allochi un numero elevato per processo, le chiamate alle attività potrebbero mettersi in coda l'una dietro l'altra, poiché Nextflow non richiederà più CPU di quelle disponibili.

Eseguiamolo di nuovo con un nome di report diverso, così possiamo confrontare prima e dopo.

```bash
nextflow run main.nf -with-report report-config-2.html
```

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [voluminous_venter] revision: c3c85dec78

    executor >  local (8)
    [a1/0e96d4] sayHello (1)       | 3 of 3 ✔
    [a3/7173a3] convertToUpper (2) | 3 of 3 ✔
    [4f/a8ae3d] collectGreetings   | 1 of 1 ✔
    [91/3724f8] cowpy              | 1 of 1 ✔

    Outputs:

      /workspaces/training/config-exec/results

      first_output:
        - full_pipeline/intermediates/Bonjour-output.txt
        - full_pipeline/intermediates/Hola-output.txt
        - full_pipeline/intermediates/Hello-output.txt

      uppercased:
        - full_pipeline/intermediates/UPPER-Bonjour-output.txt
        - full_pipeline/intermediates/UPPER-Hello-output.txt
        - full_pipeline/intermediates/UPPER-Hola-output.txt

      collected: full_pipeline/intermediates/COLLECTED-batch-output.txt

      batch_report: full_pipeline/batch-report.txt

      cowpy_art: full_pipeline/cowpy-COLLECTED-batch-output.txt
    ```

Confrontando i due report per `cowpy`:

| Metrica               | Prima (1 CPU, 1 GB) | Dopo (2 CPU, 2 GB) |
| --------------------- | ------------------- | ------------------ |
| Memoria massima usata | 6.4 MB              | 6.4 MB             |
| Utilizzo CPU          | 116%                | 118%               |

Raddoppiare l'allocazione non ha cambiato affatto l'utilizzo effettivo, il che ci dice che il valore originale di 1 GB / 1 CPU era già generoso per questo carico di lavoro di prova.
Su una pipeline reale che elabora dati non banali, ci si aspetterebbe che i numeri differiscano in modo significativo tra i processi, ed è esattamente per questo che si fa il profiling prima di decidere cosa allocare, invece di tirare a indovinare.

### 1.3. Aggiungere limiti alle risorse

A seconda della tua infrastruttura di calcolo, potrebbero esserci vincoli rigidi su ciò che puoi richiedere, ad esempio un limite a livello di cluster.
La direttiva `resourceLimits` ti permette di impostare questi limiti:

```groovy title="Syntax example"
process {
    resourceLimits = [
        memory: 750.GB,
        cpus: 200,
        time: 30.d
    ]
}
```

Nextflow traduce questi valori in qualsiasi formato si aspetti l'executor di destinazione.
Se un processo richiede più del limite, la richiesta viene ridotta al limite anziché rifiutata.

!!! warning "Avviso"

    Questo non è qualcosa che puoi eseguire nell'ambiente di formazione, poiché richiede un'infrastruttura HPC per avere effetto.

??? info "Configurazioni di riferimento istituzionali"

    Il progetto nf-core mantiene una [raccolta di file di configurazione](https://nf-co.re/configs/) condivisi da istituzioni di tutto il mondo, che coprono un'ampia gamma di executor HPC e cloud.
    Sono un utile punto di partenza indipendentemente dal fatto che la tua istituzione sia tra quelle elencate.

### Takeaway

Sai come generare un report di profiling per valutare l'utilizzo delle risorse, sovrascrivere le allocazioni di risorse per un processo specifico e limitare le allocazioni con `resourceLimits`.

### Cosa c'è dopo?

Scopriamo come fare in modo che una pipeline si recuperi automaticamente quando un'attività fallisce, indipendentemente dal fatto che la nostra stima sull'allocazione delle risorse fosse corretta.

---

## 2. Gestire i fallimenti delle attività con i retry

Il profiling ti dice di cosa ha bisogno un processo nella maggior parte dei casi, ma i carichi di lavoro reali variano: un'allocazione comoda per la maggior parte degli input può comunque essere troppo ridotta per uno insolitamente grande, e le stime possono semplicemente essere sbagliate.
Invece di lasciare che una singola attività fallita mandi in crash l'intera esecuzione, Nextflow può riprovare automaticamente un'attività fallita, opzionalmente assegnandole più risorse a ogni tentativo.

### 2.1. Riprovare automaticamente un'attività fallita

Per vedere questo in azione, impostiamo deliberatamente l'allocazione di memoria di `cowpy` al di sotto di ciò di cui ha effettivamente bisogno: ricordiamo dalla [1.1](#11-generate-a-resource-utilization-report) che raggiunge un picco di circa 6.4 MB, quindi 6 MB dovrebbero essere appena insufficienti.

=== "Dopo"

    ```groovy title="nextflow.config" linenums="6" hl_lines="5-7"
    process {
        memory = 1.GB
        withName: 'cowpy' {
            conda = 'conda-forge::cowpy==1.1.5'
            memory = 6.MB
            errorStrategy = 'retry'
            maxRetries = 2
        }
    }
    ```

=== "Prima"

    ```groovy title="nextflow.config" linenums="6"
    process {
        memory = 1.GB
        withName: 'cowpy' {
            conda = 'conda-forge::cowpy==1.1.5'
            memory = 2.GB
            cpus = 2
        }
    }
    ```

`errorStrategy` dice a Nextflow cosa fare quando un'attività fallisce: `'retry'` risottomette l'attività invece di fermare l'intera pipeline.
`maxRetries` limita il numero di tentativi aggiuntivi prima che Nextflow si arrenda.

```bash
nextflow run main.nf
```

??? failure "Output del comando (abbreviato)"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [desperate_brazil] revision: c3c85dec78

    executor >  local (10)
    [67/fe1f49] sayHello (1)       | 3 of 3 ✔
    [8a/f13335] convertToUpper (1) | 3 of 3 ✔
    [39/1b24ed] collectGreetings   | 1 of 1 ✔
    [7a/d5eb6f] cowpy              | 0 of 1, retries: 2 ✘
    [9d/b79eb3] NOTE: Process `cowpy` terminated with an error exit status (137) -- Execution is retried (1)
    [6d/1d9d84] NOTE: Process `cowpy` terminated with an error exit status (137) -- Execution is retried (2)
    ERROR ~ Error executing process > 'cowpy'

    Caused by:
      Process `cowpy` terminated with an error exit status (137)

    Command executed:
      cat COLLECTED-batch-output.txt | cowpy -c "turkey" > cowpy-COLLECTED-batch-output.txt

    Command exit status:
      137

    Command output:
      (empty)

    Command error:
      /usr/local/bin/_activate_current_env.sh: line 35:    14 Killed                  micromamba activate "${ENV_NAME:-base}"

    Work dir:
      /workspaces/training/config-exec/work/7a/d5eb6feeac0eed18d95d3da7a7aeb4

    Tip: when you have fixed the problem you can continue the execution adding the option `-resume` to the run command line

    -- Check '.nextflow.log' file for details
    ```

Il codice di uscita 137 è il segnale standard per un'interruzione per mancanza di memoria: il container non aveva abbastanza memoria per eseguire `cowpy`.
Nextflow ha riprovato l'attività due volte, tre tentativi in totale, corrispondenti a `maxRetries = 2`.
Poiché l'allocazione di memoria non è mai cambiata tra i tentativi, ogni tentativo ha incontrato lo stesso ostacolo; una volta esauriti i retry, Nextflow riporta il fallimento per intero e ferma la pipeline, uscendo con uno stato non zero.

Riprovare da solo non risolve nulla se la causa sottostante non cambia tra i tentativi.

### 2.2. Aumentare le risorse a ogni retry

All'interno di una direttiva di processo, `task.attempt` contiene il numero del tentativo corrente, a partire da 1.
Puoi usarlo in una closure per scalare un'allocazione di risorse verso l'alto a ogni retry.

=== "Dopo"

    ```groovy title="nextflow.config" linenums="6" hl_lines="5 7"
    process {
        memory = 1.GB
        withName: 'cowpy' {
            conda = 'conda-forge::cowpy==1.1.5'
            memory = { 6.MB * task.attempt }
            errorStrategy = 'retry'
            maxRetries = 3
        }
    }
    ```

=== "Prima"

    ```groovy title="nextflow.config" linenums="6"
    process {
        memory = 1.GB
        withName: 'cowpy' {
            conda = 'conda-forge::cowpy==1.1.5'
            memory = 6.MB
            errorStrategy = 'retry'
            maxRetries = 2
        }
    }
    ```

Eseguiamo di nuovo il flusso di lavoro:

```bash
nextflow run main.nf
```

??? success "Output del comando (abbreviato)"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [grave_joliot] revision: c3c85dec78

    executor >  local (9)
    [0f/211b8a] sayHello (2)       | 3 of 3 ✔
    [26/301ea2] convertToUpper (3) | 3 of 3 ✔
    [22/5a895b] collectGreetings   | 1 of 1 ✔
    [e1/beee86] cowpy              | 1 of 1, retries: 1 ✔
    [b6/7aed6a] NOTE: Process `cowpy` terminated with an error exit status (137) -- Execution is retried (1)

    Outputs:

      /workspaces/training/config-exec/results

      first_output:
        - full_pipeline/intermediates/Hola-output.txt
        - full_pipeline/intermediates/Hello-output.txt
        - full_pipeline/intermediates/Bonjour-output.txt

      uppercased:
        - full_pipeline/intermediates/UPPER-Hola-output.txt
        - full_pipeline/intermediates/UPPER-Hello-output.txt
        - full_pipeline/intermediates/UPPER-Bonjour-output.txt

      collected: full_pipeline/intermediates/COLLECTED-batch-output.txt

      batch_report: full_pipeline/batch-report.txt

      cowpy_art: full_pipeline/cowpy-COLLECTED-batch-output.txt
    ```

Il primo tentativo fallisce ancora a 6 MB, ma il retry viene eseguito con 12 MB (`6.MB * 2`) e ha successo, e la pipeline si completa con tutti gli output pubblicati.

!!! warning "Avviso"

    L'output della console include ancora una riga `NOTE:` che riporta il primo tentativo fallito, anche se la pipeline nel suo complesso ha avuto successo: Nextflow registra ogni retry individualmente, ma un fallimento riprovato non influisce sull'esito complessivo.
    Controlla il riepilogo `Outputs:`, o lo stato di uscita del comando, per confermare se l'esecuzione ha effettivamente avuto successo.

Consulta [Dynamic computing resources](https://nextflow.io/docs/latest/process.html#dynamic-task-resources) nella documentazione di Nextflow per pattern di retry più avanzati, incluso il ridimensionamento basato sul tipo specifico di errore verificatosi.

### Takeaway

Sai come fare in modo che una pipeline riprovi automaticamente le attività fallite e come scalare le allocazioni di risorse a ogni retry usando `task.attempt`.

### Cosa c'è dopo?

Passa alla [Parte 3](./03_profiles.md), dove imparerai come raggruppare configurazioni come questa in profili intercambiabili.

---

## Riepilogo

In questa parte hai imparato a:

- Generare un report di profiling delle risorse e impostare allocazioni di risorse per processo
- Limitare le richieste di risorse con `resourceLimits`
- Riprovare automaticamente un'attività fallita con `errorStrategy` e `maxRetries`
- Scalare un'allocazione di risorse verso l'alto a ogni retry usando `task.attempt`
