# Ambiente di Sviluppo

Gli Ambienti di Sviluppo Integrati (IDE) moderni possono trasformare radicalmente la vostra esperienza di sviluppo con Nextflow. Questa side quest si concentra specificamente sull'utilizzo di VS Code e della sua estensione Nextflow per scrivere codice più velocemente, individuare errori in anticipo e navigare flussi di lavoro complessi in modo efficiente.

!!! note "Questo non è un tutorial tradizionale"

    A differenza di altri moduli di formazione, questa guida è organizzata come una raccolta di suggerimenti rapidi, consigli ed esempi pratici piuttosto che come un tutorial passo-passo. Ogni sezione può essere esplorata indipendentemente in base ai vostri interessi e alle esigenze attuali di sviluppo. Sentitevi liberi di saltare da una sezione all'altra e concentrarvi sulle funzionalità che saranno più immediatamente utili per lo sviluppo dei vostri flussi di lavoro.

## Cosa dovreste sapere prima

Questa guida presuppone che abbiate completato il corso di formazione [Hello Nextflow](../hello_nextflow/) e che abbiate familiarità con i concetti fondamentali di Nextflow, tra cui:

- **Struttura base del flusso di lavoro**: Comprensione dei processi, dei flussi di lavoro e di come si connettono tra loro
- **Operazioni sui canali**: Creazione di canali, passaggio di dati tra processi e utilizzo di operatori di base
- **Moduli e organizzazione**: Creazione di moduli riutilizzabili e utilizzo delle istruzioni include
- **Nozioni di base sulla configurazione**: Utilizzo di `nextflow.config` per parametri, direttive di processo e profili

## Cosa imparerete qui

Questa guida si concentra sulle **funzionalità di produttività dell'IDE** che vi renderanno sviluppatori Nextflow più efficienti:

- **Evidenziazione avanzata della sintassi**: Comprendere cosa VS Code vi mostra sulla struttura del vostro codice
- **Auto-completamento intelligente**: Sfruttare suggerimenti contestuali per scrivere codice più velocemente
- **Rilevamento errori e diagnostica**: Individuare errori di sintassi prima di eseguire il flusso di lavoro
- **Navigazione del codice**: Spostarsi rapidamente tra processi, moduli e definizioni
- **Formattazione e organizzazione**: Mantenere uno stile di codice coerente e leggibile
- **Sviluppo assistito da IA** (opzionale): Utilizzare strumenti IA moderni integrati con il vostro IDE

!!! info "Perché le funzionalità dell'IDE ora?"

    Probabilmente avete già utilizzato VS Code durante il corso [Hello Nextflow](../hello_nextflow/), ma abbiamo mantenuto il focus sull'apprendimento dei fondamenti di Nextflow piuttosto che sulle funzionalità dell'IDE. Ora che avete familiarità con i concetti base di Nextflow come processi, flussi di lavoro, canali e moduli, siete pronti a sfruttare le sofisticate funzionalità dell'IDE che vi renderanno sviluppatori più efficienti.

    Pensate a questo come a un "salto di livello" del vostro ambiente di sviluppo - lo stesso editor che avete utilizzato ha capacità molto più potenti che diventano veramente preziose una volta che comprendete cosa vi stanno aiutando a fare.

---

## 0. Setup e Riscaldamento

Impostiamo uno spazio di lavoro specifico per esplorare le funzionalità dell'IDE:

```bash title="Navigate to the IDE features directory"
cd side-quests/ide_features
```

Aprite questa directory in VS Code:

```bash title="Open VS Code in current directory"
code .
```

La directory `ide_features` contiene flussi di lavoro di esempio che dimostrano varie funzionalità dell'IDE:

```bash title="Show directory structure"
tree .
```

```console title="Project structure"
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

    - `basic_workflow.nf` è un flusso di lavoro base funzionante che potete eseguire e modificare
    - `complex_workflow.nf` è progettato solo a scopo illustrativo per dimostrare le funzionalità di navigazione - potrebbe non essere eseguito con successo ma mostra una struttura realistica di flusso di lavoro multi-file

### Scorciatoie da Tastiera

Alcune delle funzionalità in questa guida utilizzeranno scorciatoie da tastiera opzionali. Potreste accedere a questo materiale tramite GitHub Codespaces nel browser, e in questo caso a volte le scorciatoie potrebbero non funzionare come previsto perché sono utilizzate per altre cose nel vostro sistema.

Se state eseguendo VS Code localmente, come probabilmente farete quando scriverete effettivamente i flussi di lavoro, le scorciatoie funzioneranno come descritto.

Se state usando un Mac, alcune (non tutte) scorciatoie da tastiera useranno "cmd" invece di "ctrl", e lo indicheremo nel testo come `Ctrl/Cmd`.

### 0.1. Installazione dell'Estensione Nextflow

!!! note "State Già Usando Devcontainer?"

    Se state lavorando in **GitHub Codespaces** o utilizzando un **devcontainer locale**, l'estensione Nextflow è probabilmente già installata e configurata per voi. Potete saltare i passaggi di installazione manuale qui sotto e procedere direttamente all'esplorazione delle funzionalità dell'estensione.

Per installare l'estensione manualmente:

1. Aprite VS Code
2. Andate alla vista Estensioni cliccando l'icona delle estensioni a sinistra: ![icona estensioni](img/extensions_icon.png) (scorciatoia `Ctrl/Cmd+Shift+X` se state eseguendo VSCode localmente)
3. Cercate "Nextflow"
4. Installate l'estensione ufficiale Nextflow

![Installa Estensione Nextflow](img/install_extension.png)

### 0.2. Layout dello Spazio di Lavoro

Dato che avete utilizzato VS Code durante Hello Nextflow, avete già familiarità con le basi. Ecco come organizzare il vostro spazio di lavoro in modo efficiente per questa sessione:

- **Area Editor**: Per visualizzare e modificare file. Potete dividere questa area in più pannelli per confrontare file affiancati.
- **Esplora File** clic (![icona esplora file](img/files_icon.png)) (`Ctrl/Cmd+Shift+E`): I file e le cartelle locali sul vostro sistema. Tenetelo aperto a sinistra per navigare tra i file
- **Terminale Integrato** (`Ctrl+Shift+` backtick sia per Windows che MacOS): Un terminale per interagire con il computer in basso. Usatelo per eseguire Nextflow o altri comandi.
- **Pannello Problemi** (`Ctrl+Shift+M`): VS Code mostrerà qui eventuali errori e problemi rilevati. Questo è utile per evidenziare problemi a colpo d'occhio.

Potete trascinare i pannelli o nasconderli (`Ctrl/Cmd+B` per attivare/disattivare la barra laterale) per personalizzare il vostro layout mentre lavoriamo sugli esempi.

### Takeaway

Avete configurato VS Code con l'estensione Nextflow e comprendete il layout dello spazio di lavoro per uno sviluppo efficiente.

### Cosa c'è dopo?

Imparate come l'evidenziazione della sintassi vi aiuta a comprendere la struttura del codice Nextflow a colpo d'occhio.

---

## 1. Evidenziazione della Sintassi e Struttura del Codice

Ora che il vostro spazio di lavoro è configurato, esploriamo come l'evidenziazione della sintassi di VS Code vi aiuta a leggere e scrivere codice Nextflow in modo più efficace.

### 1.1. Elementi della Sintassi Nextflow

Aprite `basic_workflow.nf` per vedere l'evidenziazione della sintassi in azione:

![Showcase Sintassi](img/syntax_showcase.png)

Notate come VS Code evidenzia:

- **Parole chiave** (`process`, `workflow`, `input`, `output`, `script`) con colori distinti
- **Stringhe letterali** e **parametri** con stili diversi
- **Commenti** in un colore attenuato
- **Variabili** e **chiamate a funzioni** con enfasi appropriata
- **Blocchi di codice** con guide di indentazione appropriate

!!! note "Colori Dipendenti dal Tema"

    I colori specifici che vedete dipenderanno dal vostro tema VS Code (modalità scura/chiara), dalle impostazioni dei colori e da eventuali personalizzazioni che avete fatto. L'importante è che i diversi elementi sintattici siano visivamente distinti l'uno dall'altro, rendendo la struttura del codice più facile da comprendere indipendentemente dallo schema di colori scelto.

### 1.2. Comprensione della Struttura del Codice

L'evidenziazione della sintassi vi aiuta a identificare rapidamente:

- **Confini dei processi**: Distinzione chiara tra diversi processi
- **Blocchi input/output**: Facile individuazione delle definizioni del flusso di dati
- **Blocchi script**: I comandi effettivi che vengono eseguiti
- **Operazioni sui canali**: Passaggi di trasformazione dei dati
- **Direttive di configurazione**: Impostazioni specifiche del processo

Questa organizzazione visiva diventa inestimabile quando si lavora con flussi di lavoro complessi contenenti più processi e flussi di dati intricati.

### Takeaway

Comprendete come l'evidenziazione della sintassi di VS Code vi aiuta a leggere la struttura del codice Nextflow e identificare diversi elementi del linguaggio per uno sviluppo più rapido.

### Cosa c'è dopo?

Imparate come l'auto-completamento intelligente accelera la scrittura del codice con suggerimenti contestuali.

---

## 2. Auto-completamento Intelligente

Le funzionalità di auto-completamento di VS Code vi aiutano a scrivere codice più velocemente e con meno errori suggerendo opzioni appropriate in base al contesto.

### 2.1. Suggerimenti Contestuali

Le opzioni di auto-completamento variano a seconda di dove vi trovate nel vostro codice:

#### Operazioni sui Canali

Aprite di nuovo `basic_workflow.nf` e provate a digitare `channel.` nel blocco workflow:

![Auto-completamento canali](img/autocomplete_channel.png)

Vedrete suggerimenti per:

- `fromPath()` - Crea canale da percorsi di file
- `fromFilePairs()` - Crea canale da file accoppiati
- `of()` - Crea canale da valori
- `fromSRA()` - Crea canale da accessi SRA
- E molti altri...

Questo vi aiuta a trovare rapidamente la fabbrica di canali giusta da usare senza dover ricordare i nomi esatti dei metodi.

Potete anche scoprire gli operatori disponibili da applicare ai canali. Per esempio, digitate `FASTQC.out.html.` per vedere le operazioni disponibili:

![Auto-completamento operatori canali](img/autocomplete_operators.png)

#### Direttive di Processo

All'interno di un blocco script di processo, digitate `task.` per vedere le proprietà runtime disponibili:

![Auto-completamento proprietà task](img/autocomplete_task.png)

#### Configurazione

Aprite nextflow.config e digitate `process.` ovunque per vedere le direttive di processo disponibili:

![Auto-completamento config](img/autocomplete_config.png)

Vedrete suggerimenti per:

- `executor`
- `memory`
- `cpus`

Questo fa risparmiare tempo quando si configurano i processi e funziona attraverso diversi ambiti di configurazione. Per esempio, provate a digitare `docker.` per vedere le opzioni di configurazione specifiche di Docker.

### Takeaway

Potete usare l'auto-completamento intelligente di VS Code per scoprire operazioni sui canali disponibili, direttive di processo e opzioni di configurazione senza memorizzare la sintassi.

### Cosa c'è dopo?

Imparate come il rilevamento errori in tempo reale vi aiuta a individuare problemi prima di eseguire il vostro flusso di lavoro, semplicemente leggendo il codice.

## 3. Rilevamento Errori e Diagnostica

Il rilevamento errori in tempo reale di VS Code vi aiuta a individuare problemi prima di eseguire il vostro flusso di lavoro.

### 3.1. Rilevamento Errori di Sintassi

Creiamo un errore deliberato per vedere il rilevamento in azione. Aprite `basic_workflow.nf` e cambiate il nome del processo da `FASTQC` a `FASTQ` (o qualsiasi altro nome non valido). VS Code evidenzierà immediatamente l'errore nel blocco workflow con una sottolineatura ondulata rossa:

![Sottolineatura errore](img/error_underline.png)

### 3.2. Pannello Problemi

Oltre all'evidenziazione dei singoli errori, VS Code fornisce un pannello Problemi centralizzato che aggrega tutti gli errori, avvisi e messaggi informativi nel vostro spazio di lavoro. Apritelo con `Ctrl/Cmd+Shift+M` e usate l'icona filtro per mostrare solo gli errori rilevanti per il file corrente:

![Filtra il pannello problemi](img/active_file.png)

Cliccate su qualsiasi problema per saltare direttamente alla riga problematica

![Pannello Problemi](img/problems_panel.png)

Correggete l'errore cambiando il nome del processo di nuovo in `FASTQC`.

### 3.3. Pattern di Errori Comuni

Gli errori comuni nella sintassi Nextflow includono:

- **Parentesi mancanti**: `{` o `}` non corrispondenti
- **Blocchi incompleti**: Sezioni richieste mancanti nei processi
- **Sintassi non valida**: DSL Nextflow malformato
- **Errori di battitura nelle parole chiave**: Direttive di processo scritte male
- **Mancata corrispondenza dei canali**: Incompatibilità di tipo

Il language server Nextflow evidenzia questi problemi nel pannello Problemi. Potete controllarli in anticipo per evitare errori di sintassi durante l'esecuzione di una pipeline.

### Takeaway

Potete usare il rilevamento errori di VS Code e il pannello Problemi per individuare errori di sintassi e problemi prima di eseguire il vostro flusso di lavoro, risparmiando tempo e prevenendo frustrazioni.

### Cosa c'è dopo?

Imparate come navigare in modo efficiente tra processi, moduli e definizioni in flussi di lavoro complessi.

---

## 4. Navigazione del Codice e Gestione dei Simboli

Una navigazione efficiente è cruciale quando si lavora con flussi di lavoro complessi che si estendono su più file. Per comprendere questo, sostituite la definizione del processo in `basic_workflow.nf` con un import per il modulo che vi abbiamo fornito:

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

Se passate il mouse sopra un nome di processo come `FASTQC`, vedrete un popup con l'interfaccia del modulo (input e output):

![Vai alla definizione](img/syntax.png)

Questa funzionalità è particolarmente preziosa quando si scrivono flussi di lavoro, poiché vi permette di comprendere l'interfaccia del modulo senza aprire direttamente il file del modulo.

Potete navigare rapidamente a qualsiasi definizione di processo, modulo o variabile usando **Ctrl/Cmd-clic**. Passate il mouse sopra il link al file del modulo in cima allo script e seguite il link come suggerito:

![Segui link](img/follow_link.png)

La stessa cosa funziona per i nomi dei processi. Tornate a `basic_workflow.nf` e provate questo sul nome del processo `FASTQC` nel blocco workflow. Questo vi collega direttamente al nome del processo (che è lo stesso del file del modulo in questo esempio, ma potrebbe essere a metà di un file molto più grande).

Per tornare dove eravate, usate **Alt+←** (o **Ctrl+-** su Mac). Questo è un modo potente per esplorare il codice senza perdere la vostra posizione.

Ora esploriamo la navigazione in un flusso di lavoro più complesso usando `complex_workflow.nf` (il file solo illustrativo menzionato prima). Questo flusso di lavoro contiene più processi definiti in file di modulo separati, così come alcuni inline. Mentre le strutture multi-file complesse possono essere difficili da navigare manualmente, la capacità di saltare alle definizioni rende l'esplorazione molto più gestibile.

1. Aprite `complex_workflow.nf`
2. Navigate alle definizioni dei moduli
3. Usate **Alt+←** (o **Ctrl+-**) per navigare indietro
4. Navigate al nome del processo `FASTQC` nel blocco workflow. Questo vi collega direttamente al nome del processo (che è lo stesso del file del modulo in questo esempio, ma potrebbe essere a metà di un file molto più grande).
5. Navigate di nuovo indietro
6. Navigate al processo `TRIM_GALORE` nel blocco workflow. Questo è definito inline, quindi non vi porterà a un file separato, ma vi mostrerà comunque la definizione del processo, e potrete comunque navigare indietro a dove eravate.

### 4.2. Navigazione dei Simboli

Con `complex_workflow.nf` ancora aperto, potete ottenere una panoramica di tutti i simboli nel file digitando `@` nella barra di ricerca in cima a VSCode (la scorciatoia da tastiera è `Ctrl/Cmd+Shift+O`, ma potrebbe non funzionare in Codespaces). Questo apre il pannello di navigazione dei simboli, che elenca tutti i simboli nel file corrente:

![Navigazione simboli](img/symbols.png)

Questo mostra:

- Tutte le definizioni di processo
- Definizioni di flusso di lavoro (ci sono due flussi di lavoro definiti in questo file)
- Definizioni di funzioni

Iniziate a digitare per filtrare i risultati.

### 4.3. Trova Tutti i Riferimenti

Comprendere dove un processo o una variabile viene utilizzato in tutto il vostro codebase può essere molto utile. Per esempio, se volete trovare tutti i riferimenti al processo `FASTQC`, iniziate navigando alla sua definizione. Potete farlo aprendo direttamente `modules/fastqc.nf`, o usando la funzionalità di navigazione rapida di VS Code con `Ctrl/Cmd-clic` come abbiamo fatto sopra. Una volta alla definizione del processo, fate clic destro sul nome del processo `FASTQC` e selezionate "Find All References" dal menu contestuale per vedere tutte le istanze in cui viene utilizzato.

![Trova riferimenti](img/references.png)

Questa funzionalità mostra tutte le istanze in cui `FASTQC` è referenziato nel vostro spazio di lavoro, incluso il suo utilizzo nei due flussi di lavoro distinti. Questa intuizione è cruciale per valutare il potenziale impatto delle modifiche al processo `FASTQC`.

### 4.4. Pannello Outline

Il pannello Outline, situato nella barra laterale Explorer (clic ![icona Explorer](img/files_icon.png)), fornisce una panoramica conveniente di tutti i simboli nel vostro file corrente. Questa funzionalità vi permette di navigare rapidamente e gestire la struttura del vostro codice visualizzando funzioni, variabili e altri elementi chiave in una vista gerarchica.

![Pannello outline](img/outline.png)

Usate il pannello Outline per navigare rapidamente a diverse parti del vostro codice senza usare il browser dei file.

### 4.5. Visualizzazione DAG

L'estensione Nextflow di VS Code può visualizzare il vostro flusso di lavoro come un Grafo Aciclico Diretto (DAG). Questo vi aiuta a comprendere il flusso di dati e le dipendenze tra i processi. Aprite `complex_workflow.nf` e cliccate il pulsante "Preview DAG" sopra `workflow {` (il secondo blocco `workflow` in questo file):

![Anteprima DAG](img/dag_preview.png)

Questo è solo il flusso di lavoro 'entry', ma potete anche visualizzare in anteprima il DAG per i flussi di lavoro interni cliccando il pulsante "Preview DAG" sopra il workflow `RNASEQ_PIPELINE {` più in alto:

![Anteprima DAG flusso di lavoro interno](img/dag_preview_inner.png)

Per questo flusso di lavoro, potete usare i nodi nel DAG per navigare alle corrispondenti definizioni di processo nel codice. Cliccate su un nodo, e vi porterà alla definizione del processo rilevante nell'editor. Particolarmente quando un flusso di lavoro cresce fino a una grande dimensione, questo può davvero aiutarvi a navigare nel codice e comprendere come i processi sono connessi.

### Takeaway

Potete navigare flussi di lavoro complessi in modo efficiente usando vai-alla-definizione, ricerca simboli, trova riferimenti e visualizzazione DAG per comprendere la struttura del codice e le dipendenze.

### Cosa c'è dopo?

Imparate come lavorare efficacemente attraverso più file interconnessi in progetti Nextflow più grandi.

## 5. Lavorare su Più File

Lo sviluppo Nextflow reale coinvolge il lavoro con più file interconnessi. Esploriamo come VS Code vi aiuta a gestire progetti complessi in modo efficiente.

### 5.1. Navigazione Rapida dei File

Con `complex_workflow.nf` aperto, noterete che importa diversi moduli. Pratichiamo la navigazione rapida tra di essi.

Premete **Ctrl+P** (o **Cmd+P**) e iniziate a digitare "fast":

VS Code vi mostrerà i file corrispondenti. Selezionate `modules/fastqc.nf` per saltarci istantaneamente. Questo è molto più veloce che cliccare attraverso l'esplora file quando sapete approssimativamente quale file state cercando.

Provate questo con altri pattern:

- Digitate "star" per trovare il file del modulo di allineamento STAR (`star.nf`)
- Digitate "utils" per trovare il file delle funzioni di utilità (`utils.nf`)
- Digitate "config" per saltare ai file di configurazione (`nextflow.config`)

### 5.2. Editor Diviso per Sviluppo Multi-file

Quando si lavora con i moduli, spesso è necessario vedere sia il flusso di lavoro principale che le definizioni dei moduli simultaneamente. Impostiamo questo:

1. Aprite `complex_workflow.nf`
2. Aprite `modules/fastqc.nf` in una nuova scheda
3. Fate clic destro sulla scheda `modules/fastqc.nf` e selezionate "Split Right"
4. Ora potete vedere entrambi i file affiancati

![Editor diviso](img/split_editor.png)

Questo è inestimabile quando:

- Si controllano le interfacce dei moduli mentre si scrivono chiamate al flusso di lavoro, e l'anteprima non è sufficiente
- Si confrontano processi simili attraverso diversi moduli
- Si esegue il debug del flusso di dati tra flusso di lavoro e moduli

### 5.3. Ricerca a Livello di Progetto

A volte è necessario trovare dove pattern specifici sono utilizzati in tutto il vostro progetto. Premete `Ctrl/Cmd+Shift+F` per aprire il pannello di ricerca.

Provate a cercare `publishDir` nello spazio di lavoro:

![Ricerca progetto](img/project_search.png)

Questo vi mostra ogni file che usa directory di pubblicazione, aiutandovi a:

- Comprendere i pattern di organizzazione dell'output
- Trovare esempi di direttive specifiche
- Garantire coerenza tra i moduli

### Takeaway

Potete gestire progetti multi-file complessi usando navigazione rapida dei file, editor divisi e ricerca a livello di progetto per lavorare in modo efficiente attraverso flussi di lavoro e moduli.

### Cosa c'è dopo?

Imparate come le funzionalità di formattazione e manutenzione del codice mantengono i vostri flussi di lavoro organizzati e leggibili.

---

## 6. Formattazione e Manutenzione del Codice

Una corretta formattazione del codice è essenziale non solo per l'estetica ma anche per migliorare la leggibilità, la comprensione e la facilità di aggiornamento di flussi di lavoro complessi.

### 6.1. Formattazione Automatica in Azione

Aprite `basic_workflow.nf` e rovinate deliberatamente la formattazione:

- Rimuovete alcune indentazioni: Evidenziate l'intero documento e premete `shift+tab` molte volte per rimuovere quante più indentazioni possibile.
- Aggiungete spazi extra in posti casuali: nell'istruzione `channel.fromPath`, aggiungete 30 spazi dopo la `(`.
- Spezzate alcune righe in modo scomodo: Aggiungete una nuova riga tra l'operatore `.view {` e la stringa `Processing sample:` ma non aggiungete una nuova riga corrispondente prima della parentesi di chiusura `}`.

Ora premete `Shift+Alt+F` (o `Shift+Option+F` su MacOS) per formattare automaticamente:

VS Code immediatamente:

- Corregge l'indentazione per mostrare chiaramente la struttura del processo
- Allinea elementi simili in modo coerente
- Rimuove spazi bianchi non necessari
- Mantiene interruzioni di riga leggibili

Notate che la formattazione automatica potrebbe non risolvere ogni problema di stile del codice. Il language server Nextflow mira a mantenere il vostro codice ordinato, ma rispetta anche le vostre preferenze personali in certe aree. Per esempio, se rimuovete l'indentazione all'interno del blocco `script` di un processo, il formattatore la lascerà così com'è, poiché potreste preferire intenzionalmente quello stile.

Attualmente, non c'è un'applicazione rigorosa dello stile per Nextflow, quindi il language server offre una certa flessibilità. Tuttavia, applicherà costantemente regole di formattazione intorno alle definizioni di metodi e funzioni per mantenere la chiarezza.

### 6.2. Funzionalità di Organizzazione del Codice

#### Commento Rapido

Selezionate un blocco di codice nel vostro flusso di lavoro e premete **Ctrl+/** (o **Cmd+/**) per commentarlo:

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

- Disabilitare temporaneamente parti di flussi di lavoro durante lo sviluppo
- Aggiungere commenti esplicativi a operazioni sui canali complesse
- Documentare sezioni del flusso di lavoro

Usate **Ctrl+/** (o **Cmd+/**) di nuovo per decommentare il codice.

#### Ripiegamento del Codice per Panoramica

In `complex_workflow.nf`, notate le piccole frecce accanto alle definizioni dei processi. Cliccatele per ripiegare (collassare) i processi:

![Ripiegamento codice](img/code_folding.png)

Questo vi dà una panoramica ad alto livello della struttura del vostro flusso di lavoro senza perdervi nei dettagli di implementazione.

#### Corrispondenza Parentesi

Posizionate il cursore accanto a qualsiasi parentesi `{` o `}` e VS Code evidenzia la parentesi corrispondente. Usate **Ctrl+Shift+\\** (o **Cmd+Shift+\\**) per saltare tra parentesi corrispondenti.

Questo è cruciale per:

- Comprendere i confini dei processi
- Trovare parentesi mancanti o extra
- Navigare strutture di flusso di lavoro annidate

#### Selezione e Modifica Multi-riga

Per modificare più righe simultaneamente, VS Code offre potenti capacità multi-cursore:

- **Selezione multi-riga**: Tenete premuto **Ctrl+Alt** (o **Cmd+Option** per MacOS) e usate i tasti freccia per selezionare più righe
- **Indentazione multi-riga**: Selezionate più righe e usate **Tab** per indentare o **Shift+Tab** per de-indentare interi blocchi

Questo è particolarmente utile per:

- Indentare interi blocchi di processo in modo coerente
- Aggiungere commenti a più righe contemporaneamente
- Modificare definizioni di parametri simili attraverso più processi

### Takeaway

Potete mantenere codice pulito e leggibile usando formattazione automatica, funzionalità di commento, ripiegamento del codice, corrispondenza parentesi e modifica multi-riga per organizzare flussi di lavoro complessi in modo efficiente.

### Cosa c'è dopo?

Imparate come VS Code si integra con il vostro flusso di lavoro di sviluppo più ampio oltre alla semplice modifica del codice.

---

## 7. Integrazione del Flusso di Lavoro di Sviluppo

VS Code si integra bene con il vostro flusso di lavoro di sviluppo oltre alla semplice modifica del codice.

### 7.1. Integrazione Controllo Versione

!!! note "Codespaces e Integrazione Git"

    Se state lavorando in **GitHub Codespaces**, alcune funzionalità di integrazione Git potrebbero non funzionare come previsto, in particolare le scorciatoie da tastiera per il Controllo Sorgente. Potreste anche aver rifiutato di aprire la directory come repository Git durante la configurazione iniziale, il che va bene per scopi di formazione.

Se il vostro progetto è un repository git (come questo), VS Code mostra:

- File modificati con indicatori colorati
- Stato Git nella barra di stato
- Viste diff inline
- Capacità di commit e push

Aprite il pannello Controllo Sorgente usando il pulsante controllo sorgente (![icona Controllo sorgente](img/source_control_icon.png)) (`Ctrl+Shift+G` o `Cmd+Shift+G` se state lavorando con VSCode localmente) per vedere le modifiche git e preparare commit direttamente nell'editor.

![Pannello Controllo Sorgente](img/source_control.png)

### 7.2. Esecuzione e Ispezione dei Flussi di Lavoro

Eseguiamo un flusso di lavoro e poi ispezioniamo i risultati. Nel terminale integrato (`Ctrl+Shift+` backtick sia in Windows che MacOS), eseguite il flusso di lavoro base:

```bash title="Run the basic workflow"
nextflow run basic_workflow.nf --input data/sample_data.csv --output_dir results
```

Mentre il flusso di lavoro è in esecuzione, vedrete l'output in tempo reale nel terminale. Dopo il completamento, potete usare VS Code per ispezionare i risultati senza lasciare il vostro editor:

1. **Navigate alle directory di lavoro**: Usate l'esplora file o il terminale per navigare `.nextflow/work`
2. **Aprite file di log**: Cliccate sui percorsi dei file di log nell'output del terminale per aprirli direttamente in VS Code
3. **Ispezionate gli output**: Navigate le directory dei risultati pubblicati nell'esplora file
4. **Visualizzate i report di esecuzione**: Aprite i report HTML direttamente in VS Code o nel vostro browser

Questo mantiene tutto in un unico posto piuttosto che passare tra più applicazioni.

### Takeaway

Potete integrare VS Code con il controllo versione e l'esecuzione del flusso di lavoro per gestire l'intero processo di sviluppo da un'unica interfaccia.

### Cosa c'è dopo?

Vedete come tutte queste funzionalità dell'IDE lavorano insieme nel vostro flusso di lavoro di sviluppo quotidiano.

---

## 8. Riepilogo e Note Rapide

Ecco alcune note rapide su ciascuna delle funzionalità dell'IDE discusse sopra:

### 8.1. Iniziare una Nuova Funzionalità

1. **Apertura rapida file** (`Ctrl+P` o `Cmd+P`) per trovare moduli esistenti rilevanti
2. **Editor diviso** per visualizzare processi simili affiancati
3. **Navigazione simboli** (`Ctrl+Shift+O` o `Cmd+Shift+O`) per comprendere la struttura del file
4. **Auto-completamento** per scrivere nuovo codice rapidamente

### 8.2. Debug dei Problemi

1. **Pannello problemi** (`Ctrl+Shift+M` o `Cmd+Shift+M`) per vedere tutti gli errori contemporaneamente
2. **Vai alla definizione** (`Ctrl-clic` o `Cmd-clic`) per comprendere le interfacce dei processi
3. **Trova tutti i riferimenti** per vedere come i processi sono utilizzati
4. **Ricerca a livello di progetto** per trovare pattern o problemi simili

### 8.3. Refactoring e Miglioramento

1. **Ricerca a livello di progetto** (`Ctrl+Shift+F` o `Cmd+Shift+F`) per trovare pattern
2. **Formattazione automatica** (`Shift+Alt+F` o `Shift+Option+F`) per mantenere la coerenza
3. **Ripiegamento del codice** per concentrarsi sulla struttura
4. **Integrazione Git** per tracciare le modifiche

---

## Riepilogo

Avete ora avuto un tour veloce delle funzionalità dell'IDE di VS Code per lo sviluppo Nextflow. Questi strumenti vi renderanno significativamente più produttivi:

- **Riducendo gli errori** attraverso il controllo della sintassi in tempo reale
- **Accelerando lo sviluppo** con auto-completamento intelligente
- **Migliorando la navigazione** in flussi di lavoro multi-file complessi
- **Mantenendo la qualità** attraverso formattazione coerente
- **Migliorando la comprensione** attraverso evidenziazione avanzata e visualizzazione della struttura

Non ci aspettiamo che ricordiate tutto, ma ora sapete che queste funzionalità esistono e sarete in grado di trovarle quando ne avrete bisogno. Man mano che continuerete a sviluppare flussi di lavoro Nextflow, queste funzionalità dell'IDE diventeranno una seconda natura, permettendovi di concentrarvi sulla scrittura di codice di alta qualità piuttosto che lottare con sintassi e struttura.

### Cosa c'è dopo?

Applicate queste competenze dell'IDE mentre lavorate attraverso altri moduli di formazione, per esempio:

- **[nf-test](nf-test.md)**: Creare suite di test complete per i vostri flussi di lavoro
- **[Hello nf-core](../../hello_nf-core/)**: Costruire pipeline di qualità produttiva con standard della comunità

Il vero potere di queste funzionalità dell'IDE emerge mentre lavorate su progetti più grandi e complessi. Iniziate a incorporarle nel vostro flusso di lavoro gradualmente—entro poche sessioni, diventeranno una seconda natura e trasformeranno il modo in cui affrontate lo sviluppo Nextflow.

Dall'individuare errori prima che vi rallentino alla navigazione di codebase complessi con facilità, questi strumenti vi renderanno sviluppatori più sicuri ed efficienti.

Buona programmazione!
