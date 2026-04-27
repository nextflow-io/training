# Ambiente di Sviluppo

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduzione assistita da IA - [scopri di più e suggerisci miglioramenti](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

I moderni Ambienti di Sviluppo Integrati (IDE) possono trasformare radicalmente la vostra esperienza di sviluppo con Nextflow. Questa missione secondaria si concentra specificamente sull'utilizzo di VS Code e della sua estensione Nextflow per scrivere codice più velocemente, individuare gli errori in anticipo e navigare in flussi di lavoro complessi in modo efficiente.

!!! note "Questo non è un tutorial tradizionale"

    A differenza degli altri moduli di formazione, questa guida è organizzata come una raccolta di suggerimenti rapidi, consigli pratici ed esempi concreti, piuttosto che come un tutorial passo dopo passo. Ogni sezione può essere esplorata in modo indipendente in base ai vostri interessi e alle vostre esigenze di sviluppo attuali. Sentitevi liberi di saltare da una sezione all'altra e di concentrarvi sulle funzionalità che saranno più immediatamente utili per il vostro sviluppo di flussi di lavoro.

## Cosa dovreste già sapere

Questa guida presuppone che abbiate completato il corso di formazione [Hello Nextflow](../hello_nextflow/) e che abbiate familiarità con i concetti fondamentali di Nextflow, tra cui:

- **Struttura di base del flusso di lavoro**: Comprensione dei processi, dei flussi di lavoro e di come si collegano tra loro
- **Operazioni sui canali**: Creazione di canali, passaggio di dati tra processi e utilizzo degli operatori di base
- **Moduli e organizzazione**: Creazione di moduli riutilizzabili e utilizzo delle istruzioni include
- **Nozioni di configurazione**: Utilizzo di `nextflow.config` per parametri, direttive di processo e profili

## Cosa imparerete qui

Questa guida si concentra sulle **funzionalità di produttività dell'IDE** che vi renderanno sviluppatori Nextflow più efficienti:

- **Evidenziazione avanzata della sintassi**: Capire cosa VS Code vi mostra sulla struttura del vostro codice
- **Auto-completamento intelligente**: Sfruttare i suggerimenti contestuali per scrivere codice più velocemente
- **Rilevamento degli errori e diagnostica**: Individuare gli errori di sintassi prima di eseguire il vostro flusso di lavoro
- **Navigazione nel codice**: Spostarsi rapidamente tra processi, moduli e definizioni
- **Formattazione e organizzazione**: Mantenere uno stile di codice coerente e leggibile
- **Sviluppo assistito dall'IA** (opzionale): Utilizzo di moderni strumenti di IA integrati nel vostro IDE

!!! info "Perché le funzionalità dell'IDE adesso?"

    Probabilmente avete già utilizzato VS Code durante il corso [Hello Nextflow](../hello_nextflow/), ma abbiamo mantenuto il focus sull'apprendimento dei fondamentali di Nextflow piuttosto che sulle funzionalità dell'IDE. Ora che avete familiarità con i concetti di base di Nextflow come processi, flussi di lavoro, canali e moduli, siete pronti a sfruttare le sofisticate funzionalità dell'IDE che vi renderanno sviluppatori più efficienti.

    Pensate a questo come a un "salto di livello" nel vostro ambiente di sviluppo: lo stesso editor che avete già utilizzato ha capacità molto più potenti che diventano davvero preziose una volta che capite cosa vi stanno aiutando a fare.

---

## 0. Setup e Riscaldamento

Configuriamo un workspace specifico per esplorare le funzionalità dell'IDE:

```bash title="Navigate to the IDE features directory"
cd side-quests/ide_features
```

Aprite questa directory in VS Code:

```bash title="Open VS Code in current directory"
code .
```

La directory `ide_features` contiene flussi di lavoro di esempio che illustrano varie funzionalità dell'IDE:

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

!!! note "Informazioni sui file di esempio"

    - `basic_workflow.nf` è un flusso di lavoro di base funzionante che potete eseguire e modificare
    - `complex_workflow.nf` è progettato solo a scopo illustrativo per dimostrare le funzionalità di navigazione - potrebbe non essere eseguito con successo, ma mostra una struttura realistica di flusso di lavoro multi-file

### Scorciatoie da tastiera

Alcune delle funzionalità descritte in questa guida utilizzano scorciatoie da tastiera opzionali. Se accedete a questo materiale tramite GitHub Codespaces nel browser, in alcuni casi le scorciatoie potrebbero non funzionare come previsto perché vengono utilizzate per altre funzioni del sistema.

Se eseguite VS Code localmente, come probabilmente farete quando scriverete effettivamente i flussi di lavoro, le scorciatoie funzioneranno come descritto.

Se utilizzate un Mac, alcune (non tutte) le scorciatoie da tastiera useranno "cmd" invece di "ctrl", e lo indicheremo nel testo come `Ctrl/Cmd`.

### 0.1. Installazione dell'estensione Nextflow

!!! note "Già usando i Devcontainer?"

    Se state lavorando in **GitHub Codespaces** o utilizzando un **devcontainer locale**, l'estensione Nextflow è probabilmente già installata e configurata per voi. Potete saltare i passaggi di installazione manuale qui sotto e procedere direttamente all'esplorazione delle funzionalità dell'estensione.

Per installare l'estensione manualmente:

1. Aprite VS Code
2. Andate alla vista Estensioni cliccando sull'icona delle estensioni a sinistra: ![icona estensioni](img/extensions_icon.png) (scorciatoia `Ctrl/Cmd+Shift+X` se eseguite VSCode localmente)
3. Cercate "Nextflow"
4. Installate l'estensione ufficiale Nextflow

![Installa l'estensione Nextflow](img/install_extension.png)

### 0.2. Layout del Workspace

Poiché avete già utilizzato VS Code durante Hello Nextflow, avete già familiarità con le basi. Ecco come organizzare il vostro workspace in modo efficiente per questa sessione:

- **Area Editor**: Per visualizzare e modificare i file. Potete dividerla in più pannelli per confrontare i file affiancati.
- **File Explorer** (clic su ![icona file explorer](img/files_icon.png)) (`Ctrl/Cmd+Shift+E`): I file e le cartelle locali sul vostro sistema. Tenetelo aperto a sinistra per navigare tra i file
- **Terminale Integrato** (`Ctrl+Shift+` backtick sia per Windows che per MacOS): Un terminale per interagire con il computer in basso. Usatelo per eseguire Nextflow o altri comandi.
- **Pannello Problemi** (`Ctrl+Shift+M`): VS Code mostrerà qui tutti gli errori e i problemi rilevati. È utile per evidenziare i problemi a colpo d'occhio.

Potete trascinare i pannelli o nasconderli (`Ctrl/Cmd+B` per attivare/disattivare la barra laterale) per personalizzare il layout mentre lavoriamo sugli esempi.

### Takeaway

Avete VS Code configurato con l'estensione Nextflow e comprendete il layout del workspace per uno sviluppo efficiente.

### Cosa c'è dopo?

Scoprite come l'evidenziazione della sintassi vi aiuta a comprendere la struttura del codice Nextflow a colpo d'occhio.

---

## 1. Evidenziazione della Sintassi e Struttura del Codice

Ora che il vostro workspace è configurato, esploriamo come l'evidenziazione della sintassi di VS Code vi aiuta a leggere e scrivere codice Nextflow in modo più efficace.

### 1.1. Elementi di Sintassi Nextflow

Aprite `basic_workflow.nf` per vedere l'evidenziazione della sintassi in azione:

![Presentazione della sintassi](img/syntax_showcase.png)

Notate come VS Code evidenzia:

- **Parole chiave** (`process`, `workflow`, `input`, `output`, `script`) con colori distinti
- **Stringhe letterali** e **parametri** con stili diversi
- **Commenti** con un colore attenuato
- **Variabili** e **chiamate di funzione** con un'enfasi appropriata
- **Blocchi di codice** con guide di indentazione corrette

!!! note "Colori dipendenti dal tema"

    I colori specifici che vedrete dipenderanno dal tema di VS Code (modalità scura/chiara), dalle impostazioni dei colori e da eventuali personalizzazioni apportate. L'importante è che i diversi elementi sintattici siano visivamente distinti tra loro, rendendo la struttura del codice più facile da comprendere indipendentemente dallo schema di colori scelto.

### 1.2. Comprensione della Struttura del Codice

L'evidenziazione della sintassi vi aiuta a identificare rapidamente:

- **Confini dei processi**: Distinzione chiara tra i diversi processi
- **Blocchi input/output**: Facile individuazione delle definizioni del flusso di dati
- **Blocchi script**: I comandi effettivi in esecuzione
- **Operazioni sui canali**: Passaggi di trasformazione dei dati
- **Direttive di configurazione**: Impostazioni specifiche del processo

Questa organizzazione visiva diventa preziosa quando si lavora con flussi di lavoro complessi contenenti più processi e flussi di dati intricati.

### Takeaway

Capite come l'evidenziazione della sintassi di VS Code vi aiuta a leggere la struttura del codice Nextflow e a identificare i diversi elementi del linguaggio per uno sviluppo più rapido.

### Cosa c'è dopo?

Scoprite come l'auto-completamento intelligente accelera la scrittura del codice con suggerimenti contestuali.

---

## 2. Auto-completamento Intelligente

Le funzionalità di auto-completamento di VS Code vi aiutano a scrivere codice più velocemente e con meno errori, suggerendo opzioni appropriate in base al contesto.

### 2.1. Suggerimenti Contestuali

Le opzioni di auto-completamento variano a seconda di dove vi trovate nel codice:

#### Operazioni sui Canali

Aprite di nuovo `basic_workflow.nf` e provate a digitare `channel.` nel blocco workflow:

![Auto-completamento dei canali](img/autocomplete_channel.png)

Vedrete suggerimenti per:

- `fromPath()` - Crea un canale da percorsi di file
- `fromFilePairs()` - Crea un canale da file accoppiati
- `of()` - Crea un canale da valori
- `fromSRA()` - Crea un canale da accession SRA
- E molti altri...

Questo vi aiuta a trovare rapidamente la fabbrica di canali giusta da usare senza dover ricordare i nomi esatti dei metodi.

Potete anche scoprire gli operatori disponibili da applicare ai canali. Ad esempio, digitate `FASTQC.out.html.` per vedere le operazioni disponibili:

![Auto-completamento degli operatori sui canali](img/autocomplete_operators.png)

#### Direttive di Processo

All'interno di un blocco script di un processo, digitate `task.` per vedere le proprietà di runtime disponibili:

![Auto-completamento delle proprietà task](img/autocomplete_task.png)

#### Configurazione

Aprite nextflow.config e digitate `process.` in qualsiasi punto per vedere le direttive di processo disponibili:

![Auto-completamento della configurazione](img/autocomplete_config.png)

Vedrete suggerimenti per:

- `executor`
- `memory`
- `cpus`

Questo fa risparmiare tempo nella configurazione dei processi e funziona in diversi ambiti di configurazione. Ad esempio, provate a digitare `docker.` per vedere le opzioni di configurazione specifiche di Docker.

### Takeaway

Potete usare l'auto-completamento intelligente di VS Code per scoprire le operazioni sui canali, le direttive di processo e le opzioni di configurazione disponibili senza dover memorizzare la sintassi.

### Cosa c'è dopo?

Scoprite come il rilevamento degli errori in tempo reale vi aiuta a individuare i problemi prima di eseguire il vostro flusso di lavoro, semplicemente leggendo il codice.

## 3. Rilevamento degli Errori e Diagnostica

Il rilevamento degli errori in tempo reale di VS Code vi aiuta a individuare i problemi prima di eseguire il vostro flusso di lavoro.

### 3.1. Rilevamento degli Errori di Sintassi

Creiamo un errore deliberato per vedere il rilevamento in azione. Aprite `basic_workflow.nf` e cambiate il nome del processo da `FASTQC` a `FASTQ` (o qualsiasi altro nome non valido). VS Code evidenzierà immediatamente l'errore nel blocco workflow con una sottolineatura ondulata rossa:

![Sottolineatura dell'errore](img/error_underline.png)

### 3.2. Pannello Problemi

Oltre all'evidenziazione dei singoli errori, VS Code fornisce un pannello Problemi centralizzato che aggrega tutti gli errori, gli avvisi e i messaggi informativi nel vostro workspace. Apritelo con `Ctrl/Cmd+Shift+M` e usate l'icona del filtro per mostrare solo gli errori relativi al file corrente:

![Filtra il pannello problemi](img/active_file.png)

Cliccate su qualsiasi problema per passare direttamente alla riga problematica

![Pannello Problemi](img/problems_panel.png)

Correggete l'errore riportando il nome del processo a `FASTQC`.

### 3.3. Pattern di Errori Comuni

Gli errori comuni nella sintassi Nextflow includono:

- **Parentesi mancanti**: `{` o `}` non corrispondenti
- **Blocchi incompleti**: Sezioni obbligatorie mancanti nei processi
- **Sintassi non valida**: DSL Nextflow malformato
- **Errori di battitura nelle parole chiave**: Direttive di processo scritte in modo errato
- **Mancata corrispondenza dei canali**: Incompatibilità di tipo

Il language server Nextflow evidenzia questi problemi nel pannello Problemi. Potete verificarli in anticipo per evitare errori di sintassi durante l'esecuzione di una pipeline.

### Takeaway

Potete usare il rilevamento degli errori di VS Code e il pannello Problemi per individuare errori di sintassi e problemi prima di eseguire il vostro flusso di lavoro, risparmiando tempo e prevenendo frustrazioni.

### Cosa c'è dopo?

Scoprite come navigare in modo efficiente tra processi, moduli e definizioni in flussi di lavoro complessi.

---

## 4. Navigazione nel Codice e Gestione dei Simboli

Una navigazione efficiente è fondamentale quando si lavora con flussi di lavoro complessi che si estendono su più file. Per capirlo, sostituite la definizione del processo in `basic_workflow.nf` con un'importazione per il modulo che vi abbiamo fornito:

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

Se passate il mouse sul nome di un processo come `FASTQC`, vedrete un popup con l'interfaccia del modulo (input e output):

![Vai alla definizione](img/syntax.png)

Questa funzionalità è particolarmente preziosa quando si creano flussi di lavoro, poiché consente di comprendere l'interfaccia del modulo senza aprire direttamente il file del modulo.

Potete navigare rapidamente verso qualsiasi definizione di processo, modulo o variabile usando **Ctrl/Cmd-click**. Passate il mouse sul link al file del modulo nella parte superiore dello script e seguite il link come suggerito:

![Segui il link](img/follow_link.png)

Lo stesso vale per i nomi dei processi. Tornate a `basic_workflow.nf` e provate questo sul nome del processo `FASTQC` nel blocco workflow. Questo vi porta direttamente al nome del processo (che in questo esempio coincide con il file del modulo, ma potrebbe trovarsi a metà di un file molto più grande).

Per tornare al punto precedente, usate **Alt+←** (o **Ctrl+-** su Mac). Questo è un modo potente per esplorare il codice senza perdere il proprio posto.

Ora esploriamo la navigazione in un flusso di lavoro più complesso usando `complex_workflow.nf` (il file solo illustrativo menzionato in precedenza). Questo flusso di lavoro contiene più processi definiti in file di modulo separati, oltre ad alcuni inline. Sebbene le strutture multi-file complesse possano essere difficili da navigare manualmente, la possibilità di saltare alle definizioni rende l'esplorazione molto più gestibile.

1. Aprite `complex_workflow.nf`
2. Navigate verso le definizioni dei moduli
3. Usate **Alt+←** (o **Ctrl+-**) per tornare indietro
4. Navigate verso il nome del processo `FASTQC` nel blocco workflow. Questo vi porta direttamente al nome del processo (che in questo esempio coincide con il file del modulo, ma potrebbe trovarsi a metà di un file molto più grande).
5. Tornate indietro di nuovo
6. Navigate verso il processo `TRIM_GALORE` nel blocco workflow. Questo è definito inline, quindi non vi porterà a un file separato, ma vi mostrerà comunque la definizione del processo, e potrete comunque tornare al punto precedente.

### 4.2. Navigazione tra Simboli

Con `complex_workflow.nf` ancora aperto, potete ottenere una panoramica di tutti i simboli nel file digitando `@` nella barra di ricerca in cima a VSCode (la scorciatoia da tastiera è `Ctrl/Cmd+Shift+O`, ma potrebbe non funzionare in Codespaces). Questo apre il pannello di navigazione dei simboli, che elenca tutti i simboli nel file corrente:

![Navigazione tra simboli](img/symbols.png)

Questo mostra:

- Tutte le definizioni dei processi
- Le definizioni dei flussi di lavoro (in questo file sono definiti due flussi di lavoro)
- Le definizioni delle funzioni

Iniziate a digitare per filtrare i risultati.

### 4.3. Trova Tutti i Riferimenti

Capire dove un processo o una variabile viene utilizzato in tutto il codebase può essere molto utile. Ad esempio, se volete trovare tutti i riferimenti al processo `FASTQC`, iniziate navigando verso la sua definizione. Potete farlo aprendo direttamente `modules/fastqc.nf`, oppure usando la funzionalità di navigazione rapida di VS Code con `Ctrl/Cmd-click` come abbiamo fatto sopra. Una volta alla definizione del processo, fate clic destro sul nome del processo `FASTQC` e selezionate "Find All References" dal menu contestuale per vedere tutte le istanze in cui viene utilizzato.

![Trova riferimenti](img/references.png)

Questa funzionalità mostra tutte le istanze in cui `FASTQC` è referenziato nel vostro workspace, incluso il suo utilizzo nei due flussi di lavoro distinti. Questa informazione è fondamentale per valutare il potenziale impatto delle modifiche al processo `FASTQC`.

### 4.4. Pannello Outline

Il pannello Outline, situato nella barra laterale Explorer (clic su ![icona Explorer](img/files_icon.png)), fornisce una comoda panoramica di tutti i simboli nel file corrente. Questa funzionalità consente di navigare rapidamente e gestire la struttura del codice visualizzando funzioni, variabili e altri elementi chiave in una vista gerarchica.

![Pannello Outline](img/outline.png)

Usate il pannello Outline per navigare rapidamente verso diverse parti del codice senza usare il browser dei file.

### 4.5. Visualizzazione DAG

L'estensione Nextflow di VS Code può visualizzare il vostro flusso di lavoro come un Grafo Aciclico Diretto (DAG). Questo vi aiuta a comprendere il flusso di dati e le dipendenze tra i processi. Aprite `complex_workflow.nf` e cliccate sul pulsante "Preview DAG" sopra `workflow {` (il secondo blocco `workflow` in questo file):

![Anteprima DAG](img/dag_preview.png)

Questo è solo il flusso di lavoro 'entry', ma potete anche visualizzare in anteprima il DAG per i flussi di lavoro interni cliccando sul pulsante "Preview DAG" sopra il flusso di lavoro `RNASEQ_PIPELINE {` più in alto:

![Anteprima DAG flusso di lavoro interno](img/dag_preview_inner.png)

Per questo flusso di lavoro, potete usare i nodi nel DAG per navigare verso le corrispondenti definizioni di processo nel codice. Cliccate su un nodo e vi porterà alla definizione del processo rilevante nell'editor. Soprattutto quando un flusso di lavoro cresce fino a raggiungere grandi dimensioni, questo può davvero aiutarvi a navigare nel codice e a capire come i processi sono collegati.

### Takeaway

Potete navigare in flussi di lavoro complessi in modo efficiente usando vai-alla-definizione, ricerca di simboli, trova riferimenti e visualizzazione DAG per comprendere la struttura del codice e le dipendenze.

### Cosa c'è dopo?

Scoprite come lavorare efficacemente su più file interconnessi in progetti Nextflow più grandi.

## 5. Lavorare su Più File

Lo sviluppo reale con Nextflow implica lavorare con più file interconnessi. Esploriamo come VS Code vi aiuta a gestire progetti complessi in modo efficiente.

### 5.1. Navigazione Rapida tra File

Con `complex_workflow.nf` aperto, noterete che importa diversi moduli. Pratichiamo la navigazione rapida tra di essi.

Premete **Ctrl+P** (o **Cmd+P**) e iniziate a digitare "fast":

VS Code vi mostrerà i file corrispondenti. Selezionate `modules/fastqc.nf` per saltarvi istantaneamente. Questo è molto più veloce che cliccare nel file explorer quando sapete approssimativamente quale file state cercando.

Provate con altri pattern:

- Digitate "star" per trovare il file del modulo di allineamento STAR (`star.nf`)
- Digitate "utils" per trovare il file delle funzioni di utilità (`utils.nf`)
- Digitate "config" per saltare ai file di configurazione (`nextflow.config`)

### 5.2. Editor Diviso per lo Sviluppo Multi-file

Quando si lavora con i moduli, spesso è necessario vedere sia il flusso di lavoro principale che le definizioni dei moduli contemporaneamente. Configuriamo questo:

1. Aprite `complex_workflow.nf`
2. Aprite `modules/fastqc.nf` in una nuova scheda
3. Fate clic destro sulla scheda `modules/fastqc.nf` e selezionate "Split Right"
4. Ora potete vedere entrambi i file affiancati

![Editor diviso](img/split_editor.png)

Questo è prezioso quando:

- Si controllano le interfacce dei moduli mentre si scrivono le chiamate al flusso di lavoro, e l'anteprima non è sufficiente
- Si confrontano processi simili tra moduli diversi
- Si esegue il debug del flusso di dati tra flusso di lavoro e moduli

### 5.3. Ricerca nell'Intero Progetto

A volte è necessario trovare dove vengono utilizzati pattern specifici in tutto il progetto. Premete `Ctrl/Cmd+Shift+F` per aprire il pannello di ricerca.

Provate a cercare `publishDir` nell'intero workspace:

![Ricerca nel progetto](img/project_search.png)

Questo vi mostra ogni file che utilizza le directory di pubblicazione, aiutandovi a:

- Comprendere i pattern di organizzazione dell'output
- Trovare esempi di direttive specifiche
- Garantire la coerenza tra i moduli

### Takeaway

Potete gestire progetti multi-file complessi usando la navigazione rapida tra file, gli editor divisi e la ricerca nell'intero progetto per lavorare in modo efficiente tra flussi di lavoro e moduli.

### Cosa c'è dopo?

Scoprite come le funzionalità di formattazione e manutenzione del codice mantengono i vostri flussi di lavoro organizzati e leggibili.

---

## 6. Formattazione e Manutenzione del Codice

Una corretta formattazione del codice è essenziale non solo per l'estetica, ma anche per migliorare la leggibilità, la comprensione e la facilità di aggiornamento di flussi di lavoro complessi.

### 6.1. Formattazione Automatica in Azione

Aprite `basic_workflow.nf` e rovinate deliberatamente la formattazione:

- Rimuovete alcune indentazioni: Evidenziate l'intero documento e premete `shift+tab` molte volte per rimuovere quante più indentazioni possibile.
- Aggiungete spazi extra in punti casuali: nell'istruzione `channel.fromPath`, aggiungete 30 spazi dopo il `(`.
- Spezzate alcune righe in modo scomodo: Aggiungete una nuova riga tra l'operatore `.view {` e la stringa `Processing sample:` ma non aggiungete una corrispondente nuova riga prima della parentesi di chiusura `}`.

Ora premete `Shift+Alt+F` (o `Shift+Option+F` su MacOS) per la formattazione automatica:

VS Code immediatamente:

- Corregge l'indentazione per mostrare chiaramente la struttura del processo
- Allinea elementi simili in modo coerente
- Rimuove gli spazi bianchi non necessari
- Mantiene interruzioni di riga leggibili

Notate che la formattazione automatica potrebbe non risolvere ogni problema di stile del codice. Il language server Nextflow mira a mantenere il codice ordinato, ma rispetta anche le vostre preferenze personali in alcune aree. Ad esempio, se rimuovete l'indentazione all'interno del blocco `script` di un processo, il formattatore la lascerà così com'è, poiché potreste preferire intenzionalmente quello stile.

Attualmente non esiste un'applicazione rigorosa dello stile per Nextflow, quindi il language server offre una certa flessibilità. Tuttavia, applicherà in modo coerente le regole di formattazione attorno alle definizioni di metodi e funzioni per mantenere la chiarezza.

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

- Disabilitare temporaneamente parti dei flussi di lavoro durante lo sviluppo
- Aggiungere commenti esplicativi a operazioni complesse sui canali
- Documentare le sezioni del flusso di lavoro

Usate di nuovo **Ctrl+/** (o **Cmd+/**) per decommentare il codice.

#### Piegatura del Codice per una Panoramica

In `complex_workflow.nf`, notate le piccole frecce accanto alle definizioni dei processi. Cliccatele per piegare (comprimere) i processi:

![Piegatura del codice](img/code_folding.png)

Questo vi dà una panoramica di alto livello della struttura del vostro flusso di lavoro senza perdervi nei dettagli implementativi.

#### Corrispondenza delle Parentesi

Posizionate il cursore accanto a qualsiasi parentesi `{` o `}` e VS Code evidenzia la parentesi corrispondente. Usate **Ctrl+Shift+\\** (o **Cmd+Shift+\\**) per saltare tra le parentesi corrispondenti.

Questo è fondamentale per:

- Comprendere i confini dei processi
- Trovare parentesi mancanti o in eccesso
- Navigare in strutture di flusso di lavoro annidate

#### Selezione e Modifica Multi-riga

Per modificare più righe contemporaneamente, VS Code offre potenti capacità multi-cursore:

- **Selezione multi-riga**: Tenete premuto **Ctrl+Alt** (o **Cmd+Option** per MacOS) e usate i tasti freccia per selezionare più righe
- **Indentazione multi-riga**: Selezionate più righe e usate **Tab** per indentare o **Shift+Tab** per ridurre l'indentazione di interi blocchi

Questo è particolarmente utile per:

- Indentare interi blocchi di processo in modo coerente
- Aggiungere commenti a più righe contemporaneamente
- Modificare definizioni di parametri simili in più processi

### Takeaway

Potete mantenere un codice pulito e leggibile usando la formattazione automatica, le funzionalità di commento, la piegatura del codice, la corrispondenza delle parentesi e la modifica multi-riga per organizzare flussi di lavoro complessi in modo efficiente.

### Cosa c'è dopo?

Scoprite come VS Code si integra con il vostro flusso di sviluppo più ampio, al di là della semplice modifica del codice.

---

## 7. Integrazione nel Flusso di Sviluppo

VS Code si integra bene con il vostro flusso di sviluppo al di là della semplice modifica del codice.

### 7.1. Integrazione con il Controllo di Versione

!!! note "Codespaces e integrazione Git"

    Se state lavorando in **GitHub Codespaces**, alcune funzionalità di integrazione Git potrebbero non funzionare come previsto, in particolare le scorciatoie da tastiera per il Source Control. Potreste anche aver rifiutato di aprire la directory come repository Git durante la configurazione iniziale, il che va bene per scopi formativi.

Se il vostro progetto è un repository git (come questo), VS Code mostra:

- File modificati con indicatori colorati
- Stato Git nella barra di stato
- Viste diff inline
- Capacità di commit e push

Aprite il pannello Source Control usando il pulsante source control (![icona Source control](img/source_control_icon.png)) (`Ctrl+Shift+G` o `Cmd+Shift+G` se lavorate con VSCode localmente) per vedere le modifiche git e preparare i commit direttamente nell'editor.

![Pannello Source Control](img/source_control.png)

### 7.2. Esecuzione e Ispezione dei Flussi di Lavoro

Eseguiamo un flusso di lavoro e poi ispezioniamo i risultati. Nel terminale integrato (`Ctrl+Shift+` backtick sia per Windows che per MacOS), eseguite il flusso di lavoro di base:

```bash title="Run the basic workflow"
nextflow run basic_workflow.nf --input data/sample_data.csv --output_dir results
```

Mentre il flusso di lavoro è in esecuzione, vedrete l'output in tempo reale nel terminale. Al termine, potete usare VS Code per ispezionare i risultati senza lasciare l'editor:

1. **Navigate nelle directory di lavoro**: Usate il file explorer o il terminale per sfogliare `.nextflow/work`
2. **Aprite i file di log**: Cliccate sui percorsi dei file di log nell'output del terminale per aprirli direttamente in VS Code
3. **Ispezionate gli output**: Sfogliate le directory dei risultati pubblicati nel file explorer
4. **Visualizzate i report di esecuzione**: Aprite i report HTML direttamente in VS Code o nel browser

Questo mantiene tutto in un unico posto invece di passare tra più applicazioni.

### Takeaway

Potete integrare VS Code con il controllo di versione e l'esecuzione dei flussi di lavoro per gestire l'intero processo di sviluppo da un'unica interfaccia.

### Cosa c'è dopo?

Vedete come tutte queste funzionalità dell'IDE lavorano insieme nel vostro flusso di sviluppo quotidiano.

---

## 8. Riepilogo e note rapide

Ecco alcune note rapide su ciascuna delle funzionalità dell'IDE discusse sopra:

### 8.1. Iniziare una Nuova Funzionalità

1. **Apertura rapida dei file** (`Ctrl+P` o `Cmd+P`) per trovare i moduli esistenti rilevanti
2. **Editor diviso** per visualizzare processi simili affiancati
3. **Navigazione tra simboli** (`Ctrl+Shift+O` o `Cmd+Shift+O`) per comprendere la struttura del file
4. **Auto-completamento** per scrivere nuovo codice rapidamente

### 8.2. Debug dei Problemi

1. **Pannello Problemi** (`Ctrl+Shift+M` o `Cmd+Shift+M`) per vedere tutti gli errori contemporaneamente
2. **Vai alla definizione** (`Ctrl-click` o `Cmd-click`) per comprendere le interfacce dei processi
3. **Trova tutti i riferimenti** per vedere come vengono utilizzati i processi
4. **Ricerca nell'intero progetto** per trovare pattern o problemi simili

### 8.3. Refactoring e Miglioramento

1. **Ricerca nell'intero progetto** (`Ctrl+Shift+F` o `Cmd+Shift+F`) per trovare pattern
2. **Formattazione automatica** (`Shift+Alt+F` o `Shift+Option+F`) per mantenere la coerenza
3. **Piegatura del codice** per concentrarsi sulla struttura
4. **Integrazione Git** per tracciare le modifiche

---

## Riepilogo

Avete ora fatto un tour rapido delle funzionalità IDE di VS Code per lo sviluppo con Nextflow. Questi strumenti vi renderanno significativamente più produttivi grazie a:

- **Riduzione degli errori** attraverso il controllo della sintassi in tempo reale
- **Accelerazione dello sviluppo** con l'auto-completamento intelligente
- **Miglioramento della navigazione** in flussi di lavoro complessi multi-file
- **Mantenimento della qualità** attraverso una formattazione coerente
- **Miglioramento della comprensione** attraverso l'evidenziazione avanzata e la visualizzazione della struttura

Non ci aspettiamo che ricordiate tutto, ma ora che sapete che queste funzionalità esistono, sarete in grado di trovarle quando ne avrete bisogno. Man mano che continuate a sviluppare flussi di lavoro Nextflow, queste funzionalità dell'IDE diventeranno una seconda natura, permettendovi di concentrarvi sulla scrittura di codice di alta qualità piuttosto che lottare con la sintassi e la struttura.

### Cosa c'è dopo?

Applicate queste competenze IDE mentre lavorate su altri moduli di formazione, ad esempio:

- **[nf-test](nf-test.md)**: Create suite di test complete per i vostri flussi di lavoro
- **[Hello nf-core](../../hello_nf-core/)**: Costruite pipeline di qualità produttiva con gli standard della community

Il vero potere di queste funzionalità IDE emerge quando si lavora su progetti più grandi e complessi. Iniziate a incorporarle gradualmente nel vostro flusso di lavoro: nel giro di poche sessioni, diventeranno una seconda natura e trasformeranno il vostro approccio allo sviluppo con Nextflow.

Dal rilevamento degli errori prima che vi rallentino alla navigazione in codebase complessi con facilità, questi strumenti vi renderanno sviluppatori più sicuri ed efficienti.

Buona programmazione!
