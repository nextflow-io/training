# Parte 3: Gestire le esecuzioni del flusso di lavoro

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduzione assistita da IA - [scopri di più e suggerisci miglioramenti](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Man mano che eseguite e rieseguite le pipeline, accumulate la cronologia delle esecuzioni e le vecchie directory `work/`.
Nella [Parte 1](./01_run_nextflow.md#23-use-resume-to-skip-completed-work) avete già usato `-resume` per saltare il lavoro già completato.
Qui imparerete come generare report su un'esecuzione, ispezionare la cronologia delle esecuzioni passate con [`nextflow log`](https://nextflow.io/docs/latest/reference/cli.html#log), ed eliminare le vecchie directory di lavoro non più necessarie con [`nextflow clean`](https://nextflow.io/docs/latest/reference/cli.html#clean).

---

## 1. Generare report della pipeline

Nextflow può generare diversi tipi di report su un'esecuzione, ciascuno aggiunto con il proprio flag `-with-*`: un report di esecuzione (`-with-report`), una timeline di esecuzione (`-with-timeline`), un file di traccia delle attività (`-with-trace`) e un diagramma del flusso di lavoro (`-with-dag`).
Qui genereremo i primi due; per gli altri, consultate la sezione [Execution reports](https://nextflow.io/docs/latest/reports.html) nella documentazione di riferimento di Nextflow.

### 1.1. Generare un report di esecuzione

Aggiungete `-with-report` a qualsiasi comando `nextflow run` per generare un report HTML al termine della pipeline:

```bash
nextflow run main.nf --input data/greetings.csv --character turkey -with-report
```

??? success "Output del comando"

    ```console hl_lines="6 7 8 9"
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [intergalactic_dalembert] revision: ce74f81996

    executor >  local (8)
    [34/23f10f] sayHello (2)       | 3 of 3 ✔
    [af/cab69d] convertToUpper (3) | 3 of 3 ✔
    [9e/d73afb] collectGreetings   | 1 of 1 ✔
    [3c/392db0] cowpy              | 1 of 1 ✔

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

Nextflow scrive il report in un file chiamato `report-<timestamp>.html` nella directory di lavoro.
Apritelo in un browser per vedere un riepilogo dell'esecuzione, una tabella di ogni attività con il suo stato e il tempo di esecuzione, e grafici sull'utilizzo delle risorse suddivisi per processo.

La scheda **Tasks** elenca ogni attività eseguita dalla pipeline, con il nome del processo, lo stato e l'utilizzo delle risorse:

![Tabella delle attività nel report di esecuzione](img/execution_report_tasks.png)

Il report è particolarmente utile quando una pipeline impiega più tempo del previsto o un'attività fallisce: la tabella delle attività mostra esattamente dove è stato impiegato il tempo e quali attività hanno avuto successo o sono fallite.

### 1.2. Generare una timeline di esecuzione

Aggiungete `-with-timeline` a un'esecuzione per ottenere una visualizzazione in stile diagramma di Gantt di quando è stata eseguita ogni attività:

```bash
nextflow run main.nf --input data/greetings.csv --character turkey -with-timeline
```

??? success "Output del comando"

    ```console hl_lines="6 7 8 9"
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [jolly_noyce] revision: ce74f81996

    executor >  local (8)
    [ad/e92ef3] sayHello (3)       | 3 of 3 ✔
    [2a/df8a8d] convertToUpper (2) | 3 of 3 ✔
    [be/7fb72a] collectGreetings   | 1 of 1 ✔
    [63/dc9bd6] cowpy              | 1 of 1 ✔

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

Nextflow scrive la timeline in un file chiamato `timeline-<timestamp>.html`.
Apritelo in un browser per vedere una barra per ogni attività, posizionata e dimensionata in base a quando è stata eseguita e quanto tempo ha impiegato:

![Timeline di esecuzione](img/execution_timeline.png)

La timeline rende immediatamente visibile la forma "fan-out poi fan-in" della [Parte 1](./01_run_nextflow.md#31-run-the-workflow): le tre attività `sayHello` vengono eseguite in parallelo, poi le tre attività `convertToUpper`, quindi `collectGreetings` e `cowpy` vengono eseguite una dopo l'altra poiché ciascuna dipende da tutto ciò che la precede.

### Takeaway

Sapete come generare un report HTML di esecuzione con `-with-report` e una timeline di esecuzione con `-with-timeline`, e dove cercare gli altri tipi di report supportati da Nextflow.

### Cosa c'è dopo?

Imparate come ispezionare la cronologia delle esecuzioni passate.

---

## 2. Ispezionare il log delle esecuzioni passate

Che stiate sviluppando una pipeline o eseguendola in produzione, prima o poi avrete bisogno di cercare informazioni sulle esecuzioni passate.

### 2.1. Il file di cronologia

Ogni volta che lanciate un flusso di lavoro Nextflow, viene scritta una riga in un file di log chiamato `history`, all'interno di una directory nascosta chiamata `.nextflow` nella directory di lavoro corrente.

??? abstract "Contenuto del file"

    ```txt title=".nextflow/history" linenums="1"
    2026-09-13 01:47:58	2.8s	determined_ramanujan	OK	ce74f81996b8ce999853fdbb25ba7969	03078950-5011-414f-992b-8fe1bd219ac2	nextflow run main.nf --input data/greetings.csv --character turkey
    2026-09-13 01:48:03	3s	prickly_cuvier	OK	ce74f81996b8ce999853fdbb25ba7969	5a9e8bf7-a7db-4c2c-b1d8-3bd4c7266528	nextflow run main.nf --input data/greetings.csv --character tux
    2026-09-13 01:48:09	2.3s	elegant_panini	OK	ce74f81996b8ce999853fdbb25ba7969	517ccc29-9db7-4d8d-8ea0-c8db6e9ac653	nextflow run main.nf --input data/greetings.csv --character stegosaurus
    2026-09-13 01:48:14	1.5s	deadly_lamport	OK	ce74f81996b8ce999853fdbb25ba7969	517ccc29-9db7-4d8d-8ea0-c8db6e9ac653	nextflow run main.nf --input data/greetings.csv --character stegosaurus -resume
    ```

Ogni riga fornisce il timestamp, la durata, il nome dell'esecuzione, lo stato, l'ID di revisione, l'ID di sessione e la riga di comando completa per un'esecuzione lanciata da questa directory.

Osservate le ultime due righe: sono due invocazioni separate (una normale, una con `-resume`) dello stesso identico comando, e condividono lo stesso ID di sessione.
L'ID di sessione cambia solo quando si lancia un'esecuzione genuinamente nuova; l'uso di `-resume` lo mantiene, ed è così che Nextflow sa quale cache riutilizzare.

### 2.2. Usare `nextflow log` per una visualizzazione più leggibile

Leggere il file di cronologia grezzo funziona, ma `nextflow log` formatta le stesse informazioni con un'intestazione:

```bash
nextflow log
```

??? success "Output del comando"

    ```console linenums="1"
    TIMESTAMP          	DURATION	RUN NAME            	STATUS	REVISION ID	SESSION ID                          	COMMAND
    2026-09-13 01:47:58	2.8s    	determined_ramanujan	OK    	ce74f81996 	03078950-5011-414f-992b-8fe1bd219ac2	nextflow run main.nf --input data/greetings.csv --character turkey
    2026-09-13 01:48:03	3s      	prickly_cuvier      	OK    	ce74f81996 	5a9e8bf7-a7db-4c2c-b1d8-3bd4c7266528	nextflow run main.nf --input data/greetings.csv --character tux
    2026-09-13 01:48:09	2.3s    	elegant_panini      	OK    	ce74f81996 	517ccc29-9db7-4d8d-8ea0-c8db6e9ac653	nextflow run main.nf --input data/greetings.csv --character stegosaurus
    2026-09-13 01:48:14	1.5s    	deadly_lamport      	OK    	ce74f81996 	517ccc29-9db7-4d8d-8ea0-c8db6e9ac653	nextflow run main.nf --input data/greetings.csv --character stegosaurus -resume
    ```

Nextflow raggruppa le informazioni di caching utilizzate per `-resume` in `.nextflow/cache`, indicizzate per ID di sessione.
Ecco perché cercare il nome dell'esecuzione o l'ID di sessione corretto è il primo passo ogni volta che è necessario investigare o ripulire un'esecuzione passata.

### Takeaway

Sapete dove Nextflow registra la cronologia delle esecuzioni passate e come ispezionarla con `nextflow log`.

### Cosa c'è dopo?

Imparate come rimuovere le vecchie directory di lavoro non più necessarie.

---

## 3. Eliminare le vecchie directory di lavoro

Ogni esecuzione lascia le proprie directory delle attività in `work/`, anche dopo aver copiato gli output che vi interessano in `results/`.
Eseguite abbastanza pipeline durante lo sviluppo e quelle sottodirectory si accumulano, quindi Nextflow mette a disposizione `nextflow clean` per rimuovere quelle non più necessarie.

### 3.1. Determinare i criteri di eliminazione

`nextflow clean` supporta diversi modi per selezionare cosa rimuovere; consultate la [documentazione di riferimento](https://www.nextflow.io/docs/latest/reference/cli.html#clean) per l'elenco completo.
Qui elimineremo tutto dalle esecuzioni precedenti a una determinata esecuzione, usando il suo nome.

Cercate l'esecuzione più recente che volete mantenere usando `nextflow log`; nell'[esempio della sezione 2.2](#22-use-nextflow-log-for-a-friendlier-view) si tratta di `elegant_panini`, l'ultima esecuzione normale prima di quella con `-resume`.
Il nome dell'esecuzione è la stringa in due parti generata automaticamente mostrata nella riga della console `Launching (...)`, oppure nella colonna `RUN NAME` di `nextflow log`.

### 3.2. Eseguire una prova a secco

Aggiungete prima `-n` per verificare cosa eliminerebbe un determinato comando senza eliminare effettivamente nulla:

```bash
nextflow clean -before elegant_panini -n
```

??? success "Output del comando"

    ```console
    Would remove /workspaces/training/nextflow-run/work/e5/9ec3fe24deb5d1ffa478dff5e5e1d4
    Would remove /workspaces/training/nextflow-run/work/75/36b1ece86b213c331852d0a3dc659b
    Would remove /workspaces/training/nextflow-run/work/be/8157808f28699259b5ff37e4bf9f72
    Would remove /workspaces/training/nextflow-run/work/f3/d828b7ba6b315eb08e0f96168119f5
    Would remove /workspaces/training/nextflow-run/work/94/7a2d7ef67d81b6e66430f300816698
    Would remove /workspaces/training/nextflow-run/work/eb/e37c5c8638c45916de3d9794d6d22f
    Would remove /workspaces/training/nextflow-run/work/3f/f56a705a6370f6c615aff22d3f240a
    Would remove /workspaces/training/nextflow-run/work/25/4922fa35ba2c6087777ddce20e9f2f
    Would remove /workspaces/training/nextflow-run/work/a4/2de7bcc2f8fc974261f7280c1480f7
    Would remove /workspaces/training/nextflow-run/work/48/8c49bded446d4fcb78c8125318c5e3
    Would remove /workspaces/training/nextflow-run/work/8e/c9834a9bdd7c21de2991d959201a9d
    Would remove /workspaces/training/nextflow-run/work/d5/a9055dc8c817c6c1f436494685955c
    Would remove /workspaces/training/nextflow-run/work/27/9473df9e55e35c8869d9571b480926
    Would remove /workspaces/training/nextflow-run/work/d3/ebc066e1d844232b9f5ee0a7d87464
    Would remove /workspaces/training/nextflow-run/work/91/17ef1815ed06b8f177d0135d8ef0dc
    Would remove /workspaces/training/nextflow-run/work/12/0443f9fcff27a552d4bc73aaedf24a
    ```

Sono 16 directory di attività: le 8 attività dell'esecuzione `turkey` più le 8 dell'esecuzione `tux`, esattamente quante ci si aspetterebbe per due esecuzioni complete di questa pipeline a quattro processi.
L'esecuzione `elegant_panini` stessa, e le attività in cache riutilizzate dall'esecuzione con `-resume`, vengono lasciate intatte.

Il vostro output elencherà nomi di directory diversi, e il numero di righe dipende da quante esecuzioni avete effettuato; se non vedete nessuna riga, o il nome dell'esecuzione non corrisponde a nessuno nel vostro log, oppure non c'è nulla da eliminare prima di essa.

### 3.3. Procedere con l'eliminazione

Una volta che la prova a secco sembra corretta, rieseguite lo stesso comando con `-f` al posto di `-n`:

```bash
nextflow clean -before elegant_panini -f
```

??? success "Output del comando"

    ```console
    Removed /workspaces/training/nextflow-run/work/e5/9ec3fe24deb5d1ffa478dff5e5e1d4
    Removed /workspaces/training/nextflow-run/work/75/36b1ece86b213c331852d0a3dc659b
    Removed /workspaces/training/nextflow-run/work/be/8157808f28699259b5ff37e4bf9f72
    Removed /workspaces/training/nextflow-run/work/f3/d828b7ba6b315eb08e0f96168119f5
    Removed /workspaces/training/nextflow-run/work/94/7a2d7ef67d81b6e66430f300816698
    Removed /workspaces/training/nextflow-run/work/eb/e37c5c8638c45916de3d9794d6d22f
    Removed /workspaces/training/nextflow-run/work/3f/f56a705a6370f6c615aff22d3f240a
    Removed /workspaces/training/nextflow-run/work/25/4922fa35ba2c6087777ddce20e9f2f
    Removed /workspaces/training/nextflow-run/work/a4/2de7bcc2f8fc974261f7280c1480f7
    Removed /workspaces/training/nextflow-run/work/48/8c49bded446d4fcb78c8125318c5e3
    Removed /workspaces/training/nextflow-run/work/8e/c9834a9bdd7c21de2991d959201a9d
    Removed /workspaces/training/nextflow-run/work/d5/a9055dc8c817c6c1f436494685955c
    Removed /workspaces/training/nextflow-run/work/27/9473df9e55e35c8869d9571b480926
    Removed /workspaces/training/nextflow-run/work/d3/ebc066e1d844232b9f5ee0a7d87464
    Removed /workspaces/training/nextflow-run/work/91/17ef1815ed06b8f177d0135d8ef0dc
    Removed /workspaces/training/nextflow-run/work/12/0443f9fcff27a552d4bc73aaedf24a
    ```

`nextflow clean` svuota le directory delle attività ma lascia al loro posto le directory padre a due caratteri (come `e5/`).

!!! warning "Avviso"

    Eliminare le directory di lavoro delle esecuzioni passate le rimuove dalla cache di Nextflow e cancella qualsiasi output conservato solo lì.
    Questo compromette la capacità di Nextflow di riprendere l'esecuzione senza rieseguire i processi corrispondenti, quindi ripulite solo le esecuzioni da cui siete certi di non dover riprendere.
    Questo è anche il motivo per cui vale la pena pubblicare tutto ciò che vi interessa in `results/` con `mode 'copy'` piuttosto che fare affidamento sulla directory `work/` o su una modalità di pubblicazione `symlink`.

### Takeaway

Sapete come rimuovere le vecchie directory di lavoro con `nextflow clean`, e perché farlo comporta la rinuncia alla possibilità di riprendere da quelle esecuzioni.

### Cosa c'è dopo?

Imparate come eseguire pipeline direttamente da repository remoti come GitHub nella [Parte 4](./04_remote_repositories.md).

---

## Riepilogo

In questa parte avete imparato a:

- Generare un report HTML di esecuzione con `-with-report` e una timeline di esecuzione con `-with-timeline`
- Ispezionare la cronologia delle esecuzioni passate con `nextflow log`
- Rimuovere le vecchie directory di lavoro con `nextflow clean`, e comprendere il compromesso sul resume che ne deriva
