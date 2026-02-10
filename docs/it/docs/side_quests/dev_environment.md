# Ambiente di Sviluppo

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduzione assistita da IA - [scopri di più e suggerisci miglioramenti](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Gli Ambienti di Sviluppo Integrati (IDE) moderni possono trasformare radicalmente la vostra esperienza di sviluppo in Nextflow. Questa quest secondaria si concentra specificamente sull'utilizzo di VS Code e della sua estensione Nextflow per scrivere codice più velocemente, individuare errori tempestivamente e navigare workflow complessi in modo efficiente.

!!! note "Questo non è un tutorial tradizionale"

    A differenza di altri moduli di formazione, questa guida è organizzata come una raccolta di suggerimenti rapidi, consigli ed esempi pratici piuttosto che come un tutorial passo-passo. Ogni sezione può essere esplorata indipendentemente in base ai vostri interessi e alle attuali esigenze di sviluppo. Sentitevi liberi di spostarvi tra le sezioni e concentrarvi sulle funzionalità che saranno immediatamente più utili per lo sviluppo del vostro workflow.

## Cosa dovreste sapere prima

Questa guida presuppone che abbiate completato il corso di formazione [Hello Nextflow](../hello_nextflow/) e che abbiate familiarità con i concetti fondamentali di Nextflow, tra cui:

- **Struttura base del workflow**: Comprensione dei process, workflow e come si collegano tra loro
- **Operazioni sui channel**: Creazione di channel, passaggio di dati tra process e utilizzo di operatori base
- **Moduli e organizzazione**: Creazione di moduli riutilizzabili e utilizzo delle istruzioni include
- **Nozioni di base sulla configurazione**: Utilizzo di `nextflow.config` per parametri, direttive dei process e profili

## Cosa imparerete qui

Questa guida si concentra sulle **funzionalità di produttività dell'IDE** che vi renderanno sviluppatori Nextflow più efficienti:

- **Evidenziazione avanzata della sintassi**: Comprendere cosa VS Code vi sta mostrando sulla struttura del codice
- **Auto-completamento intelligente**: Sfruttare suggerimenti contestuali per una scrittura del codice più rapida
- **Rilevamento degli errori e diagnostica**: Individuare errori di sintassi prima di eseguire il workflow
- **Navigazione del codice**: Spostarvi rapidamente tra process, moduli e definizioni
- **Formattazione e organizzazione**: Mantenere uno stile di codice coerente e leggibile
- **Sviluppo assistito da AI** (facoltativo): Utilizzare strumenti AI moderni integrati con il vostro IDE

!!! info "Perché le funzionalità dell'IDE adesso?"

    Probabilmente avete già utilizzato VS Code durante il corso [Hello Nextflow](../hello_nextflow/), ma abbiamo mantenuto il focus sull'apprendimento dei fondamenti di Nextflow piuttosto che sulle funzionalità dell'IDE. Ora che avete familiarità con i concetti base di Nextflow come process, workflow, channel e moduli, siete pronti per sfruttare le sofisticate funzionalità dell'IDE che vi renderanno sviluppatori più efficienti.

    Pensate a questo come a un "potenziamento" del vostro ambiente di sviluppo - lo stesso editor che avete utilizzato ha capacità molto più potenti che diventano veramente preziose una volta compreso ciò con cui vi stanno aiutando.

---

## 0. Configurazione e Preparazione

Configuriamo uno spazio di lavoro specificamente per esplorare le funzionalità dell'IDE:

```bash title="Navigare alla directory delle funzionalità IDE"
cd side-quests/ide_features
```

Aprite questa directory in VS Code:

```bash title="Aprire VS Code nella directory corrente"
code .
```

La directory `ide_features` contiene workflow di esempio che dimostrano varie funzionalità dell'IDE:

```bash title="Mostrare la struttura della directory"
tree .
```

```console title="Struttura del progetto"
tree .
.
├── basic_workflow.nf
├── complex_workflow.nf
├── data
│   ├── sample_001.fastq.gz
│   ├── sample_002.fastq.gz
│   ├── sample_003.fastq.gz
│   ├── sample_004.fastq.gz
│   ├── sample_005.fastq.gz
│   └── sample_data.csv
├── modules
│   ├── fastqc.nf
│   ├── star.nf
│   └── utils.nf
└── nextflow.config

3 directories, 12 files
```

!!! note "Informazioni sui File di Esempio"

    - `basic_workflow.nf` è un workflow base funzionante che potete eseguire e modificare
    - `complex_workflow.nf` è progettato solo per scopi illustrativi per dimostrare le funzionalità di navigazione - potrebbe non essere eseguito con successo ma mostra una struttura realistica di workflow multi-file

### Scorciatoie da Tastiera

Alcune delle funzionalità in questa guida utilizzeranno scorciatoie da tastiera facoltative. Potreste accedere a questo materiale tramite GitHub Codespaces nel browser, e in questo caso a volte le scorciatoie potrebbero non funzionare come previsto perché sono utilizzate per altre funzioni nel vostro sistema.

Se state eseguendo VS Code localmente, come probabilmente farete quando effettivamente scriverete workflow, le scorciatoie funzioneranno come descritto.

Se state utilizzando un Mac, alcune (non tutte) scorciatoie da tastiera utilizzeranno "cmd" invece di "ctrl", e lo indicheremo nel testo come `Ctrl/Cmd`.

### 0.1. Installazione dell'Estensione Nextflow

!!! note "State già utilizzando Devcontainers?"

    Se state lavorando in **GitHub Codespaces** o utilizzando un **devcontainer locale**, l'estensione Nextflow è probabilmente già installata e configurata. Potete saltare i passaggi di installazione manuale qui sotto e procedere direttamente all'esplorazione delle funzionalità dell'estensione.

Per installare manualmente l'estensione:

1. Aprite VS Code
2. Andate alla vista Estensioni cliccando sull'icona delle estensioni a sinistra: ![icona estensioni](img/extensions_icon.png) (scorciatoia `Ctrl/Cmd+Shift+X` se state eseguendo VSCode localmente)
3. Cercate "Nextflow"
4. Installate l'estensione ufficiale Nextflow

![Installare l'Estensione Nextflow](img/install_extension.png)

### 0.2. Layout dello Spazio di Lavoro

Poiché avete utilizzato VS Code durante tutto Hello Nextflow, avete già familiarità con le nozioni di base. Ecco come organizzare il vostro spazio di lavoro in modo efficiente per questa sessione:

- **Area Editor**: Per visualizzare e modificare file. Potete dividere questa area in più pannelli per confrontare file affiancati.
- **Esplora File** click (![icona esplora file](img/files_icon.png)) (`Ctrl/Cmd+Shift+E`): I file e le cartelle locali sul vostro sistema. Tenetelo aperto a sinistra per navigare tra i file
- **Terminale Integrato** (`Ctrl+Shift+` backtick sia per Windows che per MacOS): Un terminale per interagire con il computer in basso. Utilizzatelo per eseguire Nextflow o altri comandi.
- **Pannello Problemi** (`Ctrl+Shift+M`): VS Code mostrerà qui eventuali errori e problemi rilevati. Questo è utile per evidenziare problemi a colpo d'occhio.

Potete trascinare i pannelli o nasconderli (`Ctrl/Cmd+B` per attivare/disattivare la barra laterale) per personalizzare il vostro layout mentre lavoriamo sugli esempi.

### Takeaway

Avete configurato VS Code con l'estensione Nextflow e comprendete il layout dello spazio di lavoro per uno sviluppo efficiente.

### Prossimo passo

Imparate come l'evidenziazione della sintassi vi aiuta a comprendere la struttura del codice Nextflow a colpo d'occhio.

---

## 1. Evidenziazione della Sintassi e Struttura del Codice

Ora che il vostro spazio di lavoro è configurato, esploriamo come l'evidenziazione della sintassi di VS Code vi aiuta a leggere e scrivere codice Nextflow in modo più efficace.

### 1.1. Elementi di Sintassi Nextflow

Aprite `basic_workflow.nf` per vedere l'evidenziazione della sintassi in azione:

![Dimostrazione Sintassi](img/syntax_showcase.png)

Notate come VS Code evidenzia:

- **Parole chiave** (`process`, `workflow`, `input`, `output`, `script`) con colori distinti
- **Stringhe letterali** e **parametri** con stili diversi
- **Commenti** con un colore tenue
- **Variabili** e **chiamate a funzione** con enfasi appropriata
- **Blocchi di codice** con guide di indentazione appropriate

!!! note "Colori Dipendenti dal Tema"

    I colori specifici che vedete dipenderanno dal vostro tema VS Code (modalità scura/chiara), impostazioni dei colori ed eventuali personalizzazioni effettuate. L'aspetto importante è che diversi elementi sintattici siano visivamente distinti l'uno dall'altro, rendendo la struttura del codice più facile da comprendere indipendentemente dalla combinazione di colori scelta.

### 1.2. Comprensione della Struttura del Codice

L'evidenziazione della sintassi vi aiuta a identificare rapidamente:

- **Confini dei process**: Distinzione chiara tra diversi process
- **Blocchi input/output**: Facile individuazione delle definizioni di flusso dati
- **Blocchi script**: I comandi effettivi che vengono eseguiti
- **Operazioni sui channel**: Passaggi di trasformazione dei dati
- **Direttive di configurazione**: Impostazioni specifiche dei process

Questa organizzazione visiva diventa inestimabile quando si lavora con workflow complessi contenenti molteplici process e flussi di dati intricati.

### Takeaway

Comprendete come l'evidenziazione della sintassi di VS Code vi aiuta a leggere la struttura del codice Nextflow e identificare diversi elementi del linguaggio per uno sviluppo più rapido.

### Prossimo passo

Imparate come l'auto-completamento intelligente accelera la scrittura del codice con suggerimenti contestuali.

---

## 2. Auto-completamento Intelligente

Le funzionalità di auto-completamento di VS Code vi aiutano a scrivere codice più velocemente e con meno errori suggerendo opzioni appropriate in base al contesto.

### 2.1. Suggerimenti Contestuali

Le opzioni di auto-completamento variano a seconda di dove vi trovate nel vostro codice:

#### Operazioni sui Channel

Aprite nuovamente `basic_workflow.nf` e provate a digitare `channel.` nel blocco workflow:

![Auto-completamento channel](img/autocomplete_channel.png)

Vedrà suggerimenti per:

- `fromPath()` - Creare channel da percorsi di file
- `fromFilePairs()` - Creare channel da file accoppiati
- `of()` - Creare channel da valori
- `fromSRA()` - Creare channel da accessi SRA
- E molti altri...

Questo vi aiuta a trovare rapidamente il factory channel giusto da utilizzare senza dover ricordare i nomi esatti dei metodi.

Potete anche scoprire gli operatori disponibili da applicare ai channel. Per esempio, digitate `FASTQC.out.html.` per vedere le operazioni disponibili:

![Auto-completamento operazioni channel](img/autocomplete_operators.png)

#### Direttive dei Process

All'interno di un blocco script di un process, digiti `task.` per vedere le proprietà runtime disponibili:

![Auto-completamento proprietà task](img/autocomplete_task.png)

#### Configurazione

Apra nextflow.config e digiti `process.` in qualsiasi punto per vedere le direttive dei process disponibili:

![Auto-completamento config](img/autocomplete_config.png)

Vedrà suggerimenti per:

- `executor`
- `memory`
- `cpus`

Questo fa risparmiare tempo durante la configurazione dei process e funziona attraverso diversi ambiti di configurazione. Per esempio, provi a digitare `docker.` per vedere le opzioni di configurazione specifiche di Docker.

### Takeaway

Può utilizzare l'auto-completamento intelligente di VS Code per scoprire operazioni sui channel, direttive dei process e opzioni di configurazione disponibili senza memorizzare la sintassi.

### Prossimo passo

Imparate come il rilevamento degli errori in tempo reale La aiuta a individuare problemi prima di eseguire il vostro workflow, semplicemente leggendo il codice.

## 3. Rilevamento degli Errori e Diagnostica

Il rilevamento degli errori in tempo reale di VS Code La aiuta a individuare problemi prima di eseguire il vostro workflow.

### 3.1. Rilevamento degli Errori di Sintassi

Creiamo un errore deliberato per vedere il rilevamento in azione. Apra `basic_workflow.nf` e cambi il nome del process da `FASTQC` a `FASTQ` (o qualsiasi altro nome non valido). VS Code evidenzierà immediatamente l'errore nel blocco workflow con una sottolineatura ondulata rossa:

![Sottolineatura errore](img/error_underline.png)

### 3.2. Pannello Problemi

Oltre all'evidenziazione dei singoli errori, VS Code fornisce un Pannello Problemi centralizzato che aggrega tutti gli errori, avvisi e messaggi informativi nel vostro spazio di lavoro. Lo apra con `Ctrl/Cmd+Shift+M` e utilizzi l'icona filtro per mostrare solo gli errori rilevanti per il file corrente:

![Filtrare il pannello problemi](img/active_file.png)

Clicchi su qualsiasi problema per saltare direttamente alla riga problematica

![Pannello Problemi](img/problems_panel.png)

Corregga l'errore cambiando il nome del process di nuovo in `FASTQC`.

### 3.3. Pattern di Errori Comuni

Gli errori comuni nella sintassi Nextflow includono:

- **Parentesi mancanti**: `{` o `}` non corrispondenti
- **Blocchi incompleti**: Sezioni richieste mancanti nei process
- **Sintassi non valida**: DSL Nextflow malformato
- **Errori di battitura nelle parole chiave**: Direttive dei processes scritte in modo errato
- **Mancata corrispondenza dei channel**: Incompatibilità di tipo

Il language server Nextflow evidenzia questi problemi nel pannello Problemi. Può controllarli preventivamente per evitare errori di sintassi durante l'esecuzione di una pipeline.

### Takeaway

Può utilizzare il rilevamento degli errori e il Pannello Problemi di VS Code per individuare errori di sintassi e problemi prima di eseguire il vostro workflow, risparmiando tempo e prevenendo frustrazioni.

### Prossimo passo

Imparate come navigare in modo efficiente tra process, moduli e definizioni in workflow complessi.

---

## 4. Navigazione del Codice e Gestione dei Simboli

Una navigazione efficiente è cruciale quando si lavora con workflow complessi che si estendono su più file. Per comprendere questo, sostituisca la definizione del process in `basic_workflow.nf` con un import per il modulo che Le abbiamo fornito:

=== "Dopo"

    ```groovy title="basic_workflow.nf" linenums="3"
    include { FASTQC } from './modules/fastqc.nf'
    ```

=== "Prima"

    ```groovy title="basic_workflow.nf" linenums="3"
    process FASTQC {
        tag "${sample_id}"
        publishDir "${params.output_dir}/fastqc", mode: 'copy'

        input:
        tuple val(sample_id), path(reads)

        output:
        tuple val(sample_id), path("*.html"), emit: html
        tuple val(sample_id), path("*.zip"), emit: zip

        script:
        def args = task.ext.args ?: ''
        """
        fastqc \\
            ${args} \\
            --threads ${task.cpus} \\
            ${reads}
        """
    }
    ```

### 4.1. Vai alla Definizione

Se passa il mouse sopra un nome di process come `FASTQC`, vedrà un popup con l'interfaccia del modulo (input e output):

![Vai alla definizione](img/syntax.png)

Questa funzionalità è particolarmente preziosa durante la scrittura di workflow, poiché Le consente di comprendere l'interfaccia del modulo senza aprire direttamente il file del modulo.

Può navigare rapidamente a qualsiasi definizione di process, modulo o variabile usando **Ctrl/Cmd-click**. Passi il mouse sul link al file del modulo in cima allo script e segua il link come suggerito:

![Seguire il link](img/follow_link.png)

La stessa cosa funziona per i nomi dei process. Torni a `basic_workflow.nf` e provi questo sul nome del process `FASTQC` nel blocco workflow. Questo La collega direttamente al nome del process (che è lo stesso del file del modulo in questo esempio, ma potrebbe essere a metà di un file molto più grande).

Per tornare dove era, usi **Alt+←** (o **Ctrl+-** su Mac). Questo è un modo potente per esplorare il codice senza perdere la vostra posizione.

Ora esploriamo la navigazione in un workflow più complesso usando `complex_workflow.nf` (il file solo illustrativo menzionato in precedenza). Questo workflow contiene molteplici process definiti in file di moduli separati, oltre ad alcuni inline. Sebbene strutture multi-file complesse possano essere difficili da navigare manualmente, la capacità di saltare alle definizioni rende l'esplorazione molto più gestibile.

1. Apra `complex_workflow.nf`
2. Navighi alle definizioni dei moduli
3. Usi **Alt+←** (o **Ctrl+-**) per navigare indietro
4. Navighi al nome del process `FASTQC` nel blocco workflow. Questo La collega direttamente al nome del process (che è lo stesso del file del modulo in questo esempio, ma potrebbe essere a metà di un file molto più grande).
5. Navighi nuovamente indietro
6. Navighi al process `TRIM_GALORE` nel blocco workflow. Questo è definito inline, quindi non La porterà a un file separato, ma Le mostrerà comunque la definizione del process, e potrà comunque navigare indietro a dove era.

### 4.2. Navigazione dei Simboli

Con `complex_workflow.nf` ancora aperto, può ottenere una panoramica di tutti i simboli nel file digitando `@` nella barra di ricerca in cima a VSCode (la scorciatoia da tastiera è `Ctrl/Cmd+Shift+O`, ma potrebbe non funzionare in Codespaces). Questo apre il pannello di navigazione dei simboli, che elenca tutti i simboli nel file corrente:

![Navigazione simboli](img/symbols.png)

Questo mostra:

- Tutte le definizioni dei process
- Definizioni dei workflows (ci sono due workflows definiti in questo file)
- Definizioni delle funzioni

Inizi a digitare per filtrare i risultati.

### 4.3. Trova Tutti i Riferimenti

Comprendere dove un process o una variabile viene utilizzata in tutto il vostro codebase può essere molto utile. Per esempio, se volete trovare tutti i riferimenti al process `FASTQC`, iniziate navigando alla sua definizione. Potete farlo aprendo direttamente `modules/fastqc.nf`, o utilizzando la funzione di navigazione rapida di VS Code con `Ctrl/Cmd-click` come abbiamo fatto sopra. Una volta alla definizione del process, cliccate con il tasto destro sul nome del process `FASTQC` e selezionate "Find All References" dal menu contestuale per vedere tutte le istanze dove viene utilizzato.

![Trova riferimenti](img/references.png)

Questa funzionalità visualizza tutte le istanze dove `FASTQC` è referenziato all'interno del vostro spazio di lavoro, incluso il suo utilizzo nei due workflows distinti. Questa comprensione è cruciale per valutare il potenziale impatto delle modifiche al process `FASTQC`.

### 4.4. Pannello Outline

Il pannello Outline, situato nella barra laterale Explorer (cliccate ![Icona Explorer](img/files_icon.png)), fornisce una panoramica conveniente di tutti i simboli nel vostro file corrente. Questa funzionalità vi consente di navigare e gestire rapidamente la struttura del vostro codice visualizzando funzioni, variabili e altri elementi chiave in una vista gerarchica.

![Pannello outline](img/outline.png)

Utilizzi il pannello Outline per navigare rapidamente a diverse parti del vostro codice senza usare il browser dei file.

### 4.5. Visualizzazione DAG

L'estensione Nextflow di VS Code potete visualizzare il vostro workflow come un Grafo Aciclico Diretto (DAG). Questo vi aiuta a comprendere il flusso dei dati e le dipendenze tra i processes. Aprite `complex_workflow.nf` e cliccate il pulsante "Preview DAG" sopra `workflow {` (il secondo blocco `workflow` in questo file):

![Anteprima DAG](img/dag_preview.png)

Questo è solo il workflow 'principale', ma può anche visualizzare in anteprima il DAG per i workflows interni cliccando il pulsante "Preview DAG" sopra il workflow `RNASEQ_PIPELINE {` più in alto:

![Anteprima DAG workflow interno](img/dag_preview_inner.png)

Per questo workflow, potete utilizzare i nodi nel DAG per navigare alle corrispondenti definizioni dei processes nel codice. Clicchi su un nodo, e La porterà alla definizione del process rilevante nell'editor. Particolarmente quando un workflow cresce a dimensioni considerevoli, questo può davvero aiutarvi a navigare nel codice e comprendere come i processes sono collegati.

### Takeaway

Può navigare workflow complessi in modo efficiente usando vai-alla-definizione, ricerca dei simboli, trova riferimenti e visualizzazione DAG per comprendere la struttura del codice e le dipendenze.

### Prossimo passo

Imparate come lavorare efficacemente su più file interconnessi in progetti Nextflow più grandi.

## 5. Lavorare su Più File

Lo sviluppo Nextflow reale implica lavorare con molteplici file interconnessi. Esploriamo come VS Code La aiuta a gestire progetti complessi in modo efficiente.

### 5.1. Navigazione Rapida dei File

Con `complex_workflow.nf` aperto, noterà che importa diversi moduli. Pratichiamo la navigazione rapida tra di essi.

Prema **Ctrl+P** (o **Cmd+P**) e inizi a digitare "fast":

VS Code Le mostrerà i file corrispondenti. Selezioni `modules/fastqc.nf` per saltarci istantaneamente. Questo è molto più veloce che cliccare attraverso l'esplora file quando sa approssimativamente quale file sta cercando.

Provi questo con altri pattern:

- Digiti "star" per trovare il file del modulo di allineamento STAR (`star.nf`)
- Digiti "utils" per trovare il file delle funzioni di utilità (`utils.nf`)
- Digiti "config" per saltare ai file di configurazione (`nextflow.config`)

### 5.2. Editor Diviso per Sviluppo Multi-file

Quando si lavora con i moduli, spesso è necessario vedere sia il workflow principale che le definizioni dei moduli simultaneamente. Configuriamo questo:

1. Apra `complex_workflow.nf`
2. Apra `modules/fastqc.nf` in una nuova scheda
3. Clicchi con il tasto destro sulla scheda `modules/fastqc.nf` e selezioni "Split Right"
4. Ora potete vedere entrambi i file affiancati

![Editor diviso](img/split_editor.png)

Questo è inestimabile quando:

- Si controllano le interfacce dei moduli durante la scrittura delle chiamate del workflow, e l'anteprima non è sufficiente
- Si confrontano processes simili tra diversi moduli
- Si effettua il debug del flusso dati tra workflow e moduli

### 5.3. Ricerca a Livello di Progetto

A volte è necessario trovare dove pattern specifici sono utilizzati nell'intero progetto. Prema `Ctrl/Cmd+Shift+F` per aprire il pannello di ricerca.

Provi a cercare `publishDir` nell'intero spazio di lavoro:

![Ricerca progetto](img/project_search.png)

Questo Le mostra ogni file che utilizza directory di pubblicazione, aiutandoLa a:

- Comprendere i pattern di organizzazione dell'output
- Trovare esempi di direttive specifiche
- Garantire coerenza tra i moduli

### Takeaway

Può gestire progetti multi-file complessi usando la navigazione rapida dei file, editor divisi e ricerca a livello di progetto per lavorare in modo efficiente su workflows e moduli.

### Prossimo passo

Imparate come le funzionalità di formattazione e manutenzione del codice mantengono i vostri workflows organizzati e leggibili.

---

## 6. Formattazione e Manutenzione del Codice

Una corretta formattazione del codice è essenziale non solo per l'estetica ma anche per migliorare la leggibilità, la comprensione e la facilità di aggiornamento di workflows complessi.

### 6.1. Formattazione Automatica in Azione

Apra `basic_workflow.nf` e rovini deliberatamente la formattazione:

- Rimuova alcune indentazioni: Evidenzi l'intero documento e prema `shift+tab` molte volte per rimuovere quante più indentazioni possibile.
- Aggiunga spazi extra in punti casuali: nell'istruzione `channel.fromPath`, aggiunga 30 spazi dopo `(`.
- Interrompa alcune righe in modo scomodo: Aggiunga una nuova riga tra l'operatore `.view {` e la stringa `Processing sample:` ma non aggiunga una nuova riga corrispondente prima della parentesi di chiusura `}`.

Ora prema `Shift+Alt+F` (o `Shift+Option+F` su MacOS) per la formattazione automatica:

VS Code immediatamente:

- Corregge l'indentazione per mostrare chiaramente la struttura del process
- Allinea elementi simili in modo coerente
- Rimuove spazi bianchi non necessari
- Mantiene interruzioni di riga leggibili

Noti che la formattazione automatica potrebbe non risolvere ogni problema di stile del codice. Il language server Nextflow mira a mantenere il vostro codice ordinato, ma rispetta anche le vostre preferenze personali in certe aree. Per esempio, se rimuove l'indentazione all'interno del blocco `script` di un process, il formattatore la lascerà così com'è, poiché potrebbe preferire intenzionalmente quello stile.

Attualmente, non esiste un'applicazione rigorosa dello stile per Nextflow, quindi il language server offre una certa flessibilità. Tuttavia, applicherà in modo coerente le regole di formattazione attorno alle definizioni di metodi e funzioni per mantenere la chiarezza.

### 6.2. Funzionalità di Organizzazione del Codice

#### Commento Rapido

Selezionate un blocco di codice nel vostro workflow e prema **Ctrl+/** (o **Cmd+/**) per commentarlo:

```groovy
// workflow {
//     ch_input = channel.fromPath(params.input)
//         .splitCsv(header: true)
//         .map { row -> [row.sample_id, file(row.fastq_path)] }
//
//     FASTQC(ch_input)
// }
```

Questo è perfetto per:

- Disabilitare temporaneamente parti di workflows durante lo sviluppo
- Aggiungere commenti esplicativi a operazioni sui channel complesse
- Documentare sezioni del workflow

Usi nuovamente **Ctrl+/** (o **Cmd+/**) per decommentare il codice.

#### Ripiegamento del Codice per Panoramica

In `complex_workflow.nf`, notate le piccole frecce accanto alle definizioni dei processes. Cliccatele per ripiegare (comprimere) i processes:

![Ripiegamento codice](img/code_folding.png)

Questo Le dà una panoramica di alto livello della struttura del vostro workflow senza perdersi nei dettagli di implementazione.

#### Corrispondenza delle Parentesi

Posizioni il cursore accanto a qualsiasi parentesi `{` o `}` e VS Code evidenzia la parentesi corrispondente. Usi **Ctrl+Shift+\\** (o **Cmd+Shift+\\**) per saltare tra parentesi corrispondenti.

Questo è cruciale per:

- Comprendere i confini dei processes
- Trovare parentesi mancanti o extra
- Navigare strutture di workflow annidate

#### Selezione e Modifica Multi-linea

Per modificare molteplici righe simultaneamente, VS Code offre potenti capacità multi-cursore:

- **Selezione multi-linea**: Tenga premuto **Ctrl+Alt** (o **Cmd+Option** per MacOS) e usi i tasti freccia per selezionare molteplici righe
- **Indentazione multi-linea**: Selezioni molteplici righe e usi **Tab** per indentare o **Shift+Tab** per rimuovere l'indentazione di interi blocchi

Questo è particolarmente utile per:

- Indentare interi blocchi di processes in modo coerente
- Aggiungere commenti a molteplici righe contemporaneamente
- Modificare definizioni di parametri simili attraverso molteplici processes

### Takeaway

Può mantenere un codice pulito e leggibile usando formattazione automatica, funzionalità di commento, ripiegamento del codice, corrispondenza delle parentesi e modifica multi-linea per organizzare workflow complessi in modo efficiente.

### Prossimo passo

Imparate come VS Code si integra con il vostro workflow di sviluppo più ampio oltre alla semplice modifica del codice.

---

## 7. Integrazione del Workflow di Sviluppo

VS Code si integra bene con il vostro workflow di sviluppo oltre alla semplice modifica del codice.

### 7.1. Integrazione del Controllo Versione

!!! note "Codespaces e Integrazione Git"

    Se sta lavorando in **GitHub Codespaces**, alcune funzionalità di integrazione Git potrebbero non funzionare come previsto, in particolare le scorciatoie da tastiera per Source Control. Potrebbe anche aver rifiutato di aprire la directory come repository Git durante la configurazione iniziale, il che va bene per scopi di formazione.

Se il vostro progetto è un repository git (come lo è questo), VS Code mostra:

- File modificati con indicatori colorati
- Stato Git nella barra di stato
- Viste diff inline
- Capacità di commit e push

Aprite il pannello Source Control usando il pulsante source control (![Icona Source control](img/source_control_icon.png)) (`Ctrl+Shift+G` o `Cmd+Shift+G` se state lavorando con VSCode localmente) per vedere le modifiche git e effettuare commit direttamente nell'editor.

![Pannello Source Control](img/source_control.png)

### 7.2. Esecuzione e Ispezione dei Workflows

Eseguiamo un workflow e quindi ispezioniamo i risultati. Nel terminale integrato (`Ctrl+Shift+` backtick sia in Windows che in MacOS), eseguite il workflow base:

```bash title="Eseguire il workflow base"
nextflow run basic_workflow.nf --input data/sample_data.csv --output_dir results
```

Mentre il workflow è in esecuzione, vedrà l'output in tempo reale nel terminale. Dopo il completamento, potete utilizzare VS Code per ispezionare i risultati senza lasciare il vostro editor:

1. **Navigare alle directory di lavoro**: Utilizzi l'esplora file o il terminale per navigare `.nextflow/work`
2. **Aprire file di log**: Clicchi sui percorsi dei file di log nell'output del terminale per aprirli direttamente in VS Code
3. **Ispezionare gli output**: Navighi le directory dei risultati pubblicati nell'esplora file
4. **Visualizzare report di esecuzione**: Apra report HTML direttamente in VS Code o nel vostro browser

Questo mantiene tutto in un unico posto piuttosto che cambiare tra molteplici applicazioni.

### Takeaway

Può integrare VS Code con il controllo versione e l'esecuzione del workflow per gestire l'intero processo di sviluppo da un'unica interfaccia.

### Prossimo passo

Vedete come tutte queste funzionalità dell'IDE lavorano insieme nel vostro workflow di sviluppo quotidiano.

---

## 8. Riepilogo e Note Rapide

Ecco alcune note rapide su ciascuna delle funzionalità dell'IDE discusse sopra:

### 8.1. Iniziare una Nuova Funzionalità

1. **Apertura rapida file** (`Ctrl+P` o `Cmd+P`) per trovare moduli esistenti rilevanti
2. **Editor diviso** per visualizzare processes simili affiancati
3. **Navigazione simboli** (`Ctrl+Shift+O` o `Cmd+Shift+O`) per comprendere la struttura del file
4. **Auto-completamento** per scrivere nuovo codice rapidamente

### 8.2. Debug dei Problemi

1. **Pannello problemi** (`Ctrl+Shift+M` o `Cmd+Shift+M`) per vedere tutti gli errori contemporaneamente
2. **Vai alla definizione** (`Ctrl-click` o `Cmd-click`) per comprendere le interfacce dei processes
3. **Trova tutti i riferimenti** per vedere come i processes vengono utilizzati
4. **Ricerca a livello di progetto** per trovare pattern o problemi simili

### 8.3. Refactoring e Miglioramento

1. **Ricerca a livello di progetto** (`Ctrl+Shift+F` o `Cmd+Shift+F`) per trovare pattern
2. **Formattazione automatica** (`Shift+Alt+F` o `Shift+Option+F`) per mantenere la coerenza
3. **Ripiegamento del codice** per concentrarsi sulla struttura
4. **Integrazione Git** per tracciare le modifiche

---

## Riepilogo

Ha ora avuto un tour rapido delle funzionalità dell'IDE di VS Code per lo sviluppo Nextflow. Questi strumenti La renderanno significativamente più produttivo:

- **Riducendo gli errori** attraverso il controllo della sintassi in tempo reale
- **Accelerando lo sviluppo** con auto-completamento intelligente
- **Migliorando la navigazione** in workflows multi-file complessi
- **Mantenendo la qualità** attraverso formattazione coerente
- **Migliorando la comprensione** attraverso evidenziazione avanzata e visualizzazione della struttura

Non ci aspettiamo che ricordi tutto, ma ora sa che queste funzionalità esistono e sarà in grado di trovarle quando ne avrà bisogno. Man mano che continua a sviluppare workflows Nextflow, queste funzionalità dell'IDE diventeranno una seconda natura, permettendoLe di concentrarsi sulla scrittura di codice di alta qualità piuttosto che lottare con sintassi e struttura.

### Prossimo passo

Applichi queste competenze dell'IDE mentre lavora attraverso altri moduli di formazione, per esempio:

- **[nf-test](nf-test.md)**: Creare suite di test complete per i vostri workflows
- **[Hello nf-core](../../hello_nf-core/)**: Costruire pipeline di qualità per la produzione con standard della comunità

Il vero potere di queste funzionalità dell'IDE emerge mentre lavora su progetti più grandi e complessi. Inizi a incorporarle gradualmente nel vostro workflow - nel giro di poche sessioni, diventeranno una seconda natura e trasformeranno il vostro approccio allo sviluppo Nextflow.

Dall'individuare errori prima che La rallentino alla navigazione di codebase complessi con facilità, questi strumenti La renderanno uno sviluppatore più sicuro ed efficiente.

Buona programmazione!
