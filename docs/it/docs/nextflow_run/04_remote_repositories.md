# Parte 4: Eseguire pipeline remote

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduzione assistita da IA - [scopri di più e suggerisci miglioramenti](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Finora abbiamo eseguito script di flusso di lavoro memorizzati localmente.
In pratica, spesso vorremo eseguire pipeline pubblicate in repository remoti, come GitHub, senza doverle scaricare manualmente.

Nextflow rende tutto questo semplice: è possibile eseguire qualsiasi pipeline direttamente dall'URL di un repository Git.

---

## 1. Eseguire una pipeline da GitHub

La sintassi di base per eseguire una pipeline remota è `nextflow run <repository>`, dove `<repository>` può essere il percorso di un repository GitHub come `nextflow-io/hello`, un URL completo, oppure un percorso verso GitLab, Bitbucket o un altro servizio di hosting Git.

### 1.1. Avviare la pipeline

Eseguiamo la pipeline demo ufficiale "hello" di Nextflow.
Si tratta di una pipeline diversa, molto più semplice di quella che abbiamo usato in questo corso: è precedente alla pipeline "Hello" utilizzata in questa formazione e si limita a stampare un saluto per ciascuna di alcune lingue predefinite, quindi non aspettiamoci l'input CSV o l'arte ASCII a cui siamo abituati.

```bash
nextflow run nextflow-io/hello
```

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Pulling nextflow-io/hello ...
     downloaded from https://github.com/nextflow-io/hello.git
    Launching `https://github.com/nextflow-io/hello` [sleepy_swanson] revision: 3c2cdc9823 [master]

    executor >  local (4)
    [ba/08236d] sayHello (4) | 4 of 4 ✔
    Ciao world!

    Hello world!

    Bonjour world!

    Hola world!
    ```

### 1.2. Trovare dove la pipeline viene memorizzata nella cache

La prima volta che si esegue una pipeline remota, Nextflow la scarica e la memorizza nella cache localmente.
Le esecuzioni successive riutilizzano la versione in cache, a meno che non si richieda esplicitamente un aggiornamento.

Per impostazione predefinita, Nextflow salva le pipeline scaricate in `$NXF_HOME/assets`.
Per trovare dove è stata salvata una pipeline specifica e quali revisioni sono disponibili, possiamo chiedere direttamente a Nextflow:

```bash
nextflow info nextflow-io/hello
```

??? success "Output del comando"

    ```console
     project name: nextflow-io/hello
     repository  : https://github.com/nextflow-io/hello
     local path  : /workspaces/.nextflow/assets/.repos/nextflow-io/hello
     main script : main.nf
     revisions   :
     > master (default)
       mybranch
       testing
       v1.1 [t]
       v1.2 [t]
       v1.3 [t]
    ```

    Nextflow contrassegna con `>` ogni revisione già estratta localmente; le altre sono disponibili ma non ancora recuperate in una copia di lavoro.

È anche possibile elencare tutte le pipeline scaricate finora con `nextflow list`:

```bash
nextflow list
```

??? success "Output del comando"

    ```console
    nextflow-io/hello
    ```

Il corso [Use nf-core](../nfcore_use/01_run_demo.md#12-retrieve-the-pipeline-code) approfondisce questo meccanismo di cache, incluso come esplorare il codice sorgente di una pipeline scaricata.

### Takeaway

Sappiamo come eseguire una pipeline direttamente da un repository GitHub senza doverla scaricare manualmente, e dove trovarla localmente in seguito.

### Cosa c'è dopo?

Vediamo come fissare una versione specifica di una pipeline remota per garantire la riproducibilità.

---

## 2. Specificare una versione per la riproducibilità

Per impostazione predefinita, Nextflow esegue l'ultima revisione del branch predefinito.
È possibile fissare una versione (tag), un branch o un commit specifico usando il flag `-r`.

### 2.1. Fissare una revisione specifica

```bash
nextflow run nextflow-io/hello -r v1.3
```

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Pulling nextflow-io/hello:v1.3 ...
     downloaded from https://github.com/nextflow-io/hello.git
    Launching `https://github.com/nextflow-io/hello` [sick_carson] revision: 2ce0b0e294 [v1.3]

    executor >  local (4)
    [61/e11f77] sayHello (4) | 4 of 4 ✔
    Ciao world!

    Bonjour world!

    Hello world!

    Hola world!
    ```

Nextflow recupera questa revisione la prima volta che viene richiesta, da cui le righe `Pulling` e `downloaded from`; le richieste successive della stessa revisione passano direttamente a `Launching`.
Fissare una revisione esatta è fondamentale per la riproducibilità.
Garantisce che noi e i nostri collaboratori eseguiamo esattamente lo stesso codice della pipeline, indipendentemente da cosa è cambiato nel repository nel frattempo.

### 2.2. Le revisioni si applicano solo alla singola invocazione

Fissare una revisione con `-r` influisce solo sull'esecuzione in cui viene specificata: non cambia ciò che un successivo `nextflow run` senza flag utilizzerà.
Proviamo a eseguire nuovamente la pipeline senza `-r`:

```bash
nextflow run nextflow-io/hello
```

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `https://github.com/nextflow-io/hello` [hungry_maxwell] revision: 3c2cdc9823 [master]

    executor >  local (4)
    [ba/08236d] sayHello (4) | 4 of 4 ✔
    Ciao world!

    Hello world!

    Bonjour world!

    Hola world!
    ```

Anche se l'esecuzione precedente aveva fissato esplicitamente `v1.3`, questa esecuzione torna direttamente al branch predefinito (`master`).
Nextflow mantiene una copia di lavoro locale separata per ogni revisione utilizzata, come mostrano i marcatori `>` in `nextflow info`, ma non ricorda mai quale sia stata l'ultima eseguita.
È possibile trovare il nome del branch predefinito di una pipeline eseguendo `nextflow info <pipeline>`; è quello contrassegnato con `(default)`.
La riproducibilità è interamente nostra responsabilità: passiamo sempre `-r` esplicitamente quando è importante, invece di assumere che una revisione fissata in un'esecuzione precedente sia ancora valida.

### Takeaway

Sappiamo come fissare una pipeline remota a una versione, un branch o un commit specifico per un'esecuzione riproducibile, e che il pin si applica solo a quella singola invocazione, non alle esecuzioni successive.

### Cosa c'è dopo?

Abbiamo coperto i fondamentali dell'esecuzione e della gestione delle pipeline Nextflow.
Consulta il [Riepilogo del corso](next_steps.md) per sapere come proseguire.

---

## Riepilogo

In questa parte abbiamo imparato a:

- Eseguire una pipeline direttamente da un repository GitHub senza doverla scaricare
- Fissare una pipeline remota a una revisione specifica per la riproducibilità, e capire che il pin si applica solo a quella singola invocazione
