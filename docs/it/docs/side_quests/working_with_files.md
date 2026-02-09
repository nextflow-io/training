# Elaborazione di file di input

I flussi di lavoro di analisi scientifica spesso comportano l'elaborazione di un gran numero di file.
Nextflow fornisce strumenti potenti per gestire i file in modo efficiente, aiutandovi a organizzare ed elaborare i vostri dati con un codice minimo.

### Obiettivi di apprendimento

In questa side quest, esploreremo come Nextflow gestisce i file, dalle operazioni di base alle tecniche più avanzate per lavorare con collezioni di file.
Imparerete come estrarre metadati dai nomi dei file, un requisito comune nelle pipeline di analisi scientifica.

Al termine di questa side quest, sarete in grado di:

- Creare oggetti Path da stringhe di percorsi di file utilizzando il metodo `file()` di Nextflow
- Accedere agli attributi dei file come nome, estensione e directory padre
- Gestire in modo trasparente sia file locali che remoti utilizzando URI
- Utilizzare i canali per automatizzare la gestione dei file con `channel.fromPath()` e `channel.fromFilePairs()`
- Estrarre e strutturare metadati dai nomi dei file utilizzando la manipolazione di stringhe
- Raggruppare file correlati utilizzando pattern matching ed espressioni glob
- Integrare le operazioni sui file nei processi Nextflow con una corretta gestione degli input
- Organizzare gli output dei processi utilizzando strutture di directory basate sui metadati

Queste competenze vi aiuteranno a costruire flussi di lavoro in grado di gestire diversi tipi di input di file con grande flessibilità.

### Prerequisiti

Prima di affrontare questa side quest, dovreste:

- Aver completato il tutorial [Hello Nextflow](../../hello_nextflow/) o un corso equivalente per principianti.
- Essere a vostro agio nell'uso di concetti e meccanismi di base di Nextflow (processi, canali, operatori)

---

## 0. Per iniziare

#### Aprire il codespace di formazione

Se non l'avete ancora fatto, assicuratevi di aprire l'ambiente di formazione come descritto in [Configurazione dell'Ambiente](../envsetup/index.md).

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Spostarsi nella directory del progetto

Spostiamoci nella directory dove si trovano i file per questo tutorial.

```bash
cd side-quests/working_with_files
```

Potete impostare VSCode per concentrarsi su questa directory:

```bash
code .
```

#### Esaminare i materiali

Troverete un semplice file di flusso di lavoro chiamato `main.nf`, una directory `modules` contenente due file di moduli, e una directory `data` contenente alcuni file di dati di esempio.

??? abstract "Directory contents"

    ```console
    .
    ├── data
    │   ├── patientA_rep1_normal_R1_001.fastq.gz
    │   ├── patientA_rep1_normal_R2_001.fastq.gz
    │   ├── patientA_rep1_tumor_R1_001.fastq.gz
    │   ├── patientA_rep1_tumor_R2_001.fastq.gz
    │   ├── patientA_rep2_normal_R1_001.fastq.gz
    │   ├── patientA_rep2_normal_R2_001.fastq.gz
    │   ├── patientA_rep2_tumor_R1_001.fastq.gz
    │   ├── patientA_rep2_tumor_R2_001.fastq.gz
    │   ├── patientB_rep1_normal_R1_001.fastq.gz
    │   ├── patientB_rep1_normal_R2_001.fastq.gz
    │   ├── patientB_rep1_tumor_R1_001.fastq.gz
    │   ├── patientB_rep1_tumor_R2_001.fastq.gz
    │   ├── patientC_rep1_normal_R1_001.fastq.gz
    │   ├── patientC_rep1_normal_R2_001.fastq.gz
    │   ├── patientC_rep1_tumor_R1_001.fastq.gz
    │   └── patientC_rep1_tumor_R2_001.fastq.gz
    ├── main.nf
    └── modules
        ├── analyze_reads.nf
        └── count_lines.nf
    ```

Questa directory contiene dati di sequenziamento paired-end da tre pazienti (A, B, C).

Per ogni paziente, abbiamo campioni di tipo `tumor` (tipicamente provenienti da biopsie tumorali) o `normal` (prelevati da tessuto sano o sangue).
Se non avete familiarità con l'analisi del cancro, sappiate semplicemente che questo corrisponde a un modello sperimentale che utilizza campioni tumore/normale appaiati per eseguire analisi contrastive.

Per il paziente A in particolare, abbiamo due set di replicati tecnici (ripetizioni).

I file di dati di sequenziamento sono denominati con una convenzione tipica `_R1_` e `_R2_` per quelle che sono note come 'forward reads' e 'reverse reads'.

_Non preoccupatevi se non avete familiarità con questo disegno sperimentale, non è fondamentale per comprendere questo tutorial._

#### Esaminare l'assegnazione

La vostra sfida è scrivere un flusso di lavoro Nextflow che:

1. **Carichi** i file di input utilizzando i metodi di gestione dei file di Nextflow
2. **Estragga** i metadati (ID paziente, replicato, tipo di campione) dalla struttura del nome del file
3. **Raggruppi** i file appaiati (R1/R2) insieme utilizzando `channel.fromFilePairs()`
4. **Elabori** i file con un modulo di analisi fornito
5. **Organizzi** gli output in una struttura di directory basata sui metadati estratti

#### Checklist di preparazione

Pensate di essere pronti per iniziare?

- [ ] Comprendo l'obiettivo di questo corso e i suoi prerequisiti
- [ ] Il mio codespace è attivo e funzionante
- [ ] Ho impostato la mia directory di lavoro in modo appropriato
- [ ] Comprendo l'assegnazione

Se potete spuntare tutte le caselle, siete pronti per partire.

---

## 1. Operazioni di base sui file

### 1.1. Identificare il tipo di un oggetto con `.class`

Date un'occhiata al file di flusso di lavoro `main.nf`:

```groovy title="main.nf" linenums="1"
#!/usr/bin/env nextflow

workflow {

    // Create a Path object from a string path
    myFile = 'data/patientA_rep1_normal_R1_001.fastq.gz'

    println "${myFile} is of class ${myFile.class}"
}
```

Questo è un mini-flusso di lavoro (senza processi) che fa riferimento a un singolo percorso di file nel suo flusso di lavoro, quindi lo stampa sulla console, insieme alla sua classe.

??? info "Cos'è `.class`?"

    In Nextflow, `.class` ci dice con che tipo di oggetto stiamo lavorando. È come chiedere "che tipo di cosa è questa?" per scoprire se è una stringa, un numero, un file o qualcos'altro.
    Questo ci aiuterà a illustrare la differenza tra una semplice stringa e un oggetto Path nelle prossime sezioni.

Eseguiamo il flusso di lavoro:

```bash
nextflow run main.nf
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [romantic_chandrasekhar] DSL2 - revision: 5a4a89bc3a

    data/patientA_rep1_normal_R1_001.fastq.gz is of class java.lang.String
    ```

Come potete vedere, Nextflow ha stampato la stringa di percorso esattamente come l'abbiamo scritta.

Questo è solo output di testo; Nextflow non ha ancora fatto nulla di speciale con esso.
Abbiamo anche confermato che per quanto riguarda Nextflow, questa è solo una stringa (di classe `java.lang.String`).
Ha senso, dato che non abbiamo ancora detto a Nextflow che corrisponde a un file.

### 1.2. Creare un oggetto Path con file()

Possiamo dire a Nextflow come gestire i file creando [oggetti Path](https://www.nextflow.io/docs/latest/reference/stdlib-types.html#path) da stringhe di percorso.

Nel nostro flusso di lavoro, possiamo convertire la stringa di percorso `data/patientA_rep1_normal_R1_001.fastq.gz` in un oggetto Path utilizzando il metodo `file()`, che fornisce accesso alle proprietà e alle operazioni sui file.

Modificate il `main.nf` per racchiudere la stringa con `file()` come segue:

=== "Dopo"

    ```groovy title="main.nf" linenums="5" hl_lines="2"
        // Create a Path object from a string path
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        println "${myFile} is of class ${myFile.class}"
    ```

=== "Prima"

    ```groovy title="main.nf" linenums="5" hl_lines="2"
        // Create a Path object from a string path
        myFile = 'data/patientA_rep1_normal_R1_001.fastq.gz'

        println "${myFile} is of class ${myFile.class}"
    ```

Ora eseguite nuovamente il flusso di lavoro:

```bash
nextflow run main.nf
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [kickass_coulomb] DSL2 - revision: 5af44b1b59

    /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz is of class class sun.nio.fs.UnixPath
    ```

Questa volta, vedete il percorso assoluto completo invece del percorso relativo che abbiamo fornito come input.

Nextflow ha convertito la nostra stringa in un oggetto Path e l'ha risolto nella posizione effettiva del file sul sistema.
Il percorso del file sarà ora assoluto, come in `/workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz`.

Notate anche che la classe dell'oggetto Path è `sun.nio.fs.UnixPath`: questo è il modo di Nextflow di rappresentare i file locali.
Come vedremo più avanti, i file remoti avranno nomi di classe diversi (come `nextflow.file.http.XPath` per i file HTTP), ma funzionano tutti esattamente allo stesso modo e possono essere utilizzati in modo identico nei vostri flussi di lavoro.

!!! tip

    **La differenza chiave:**

    - **Stringa di percorso**: Solo testo che Nextflow tratta come caratteri
    - **Oggetto Path**: Un riferimento intelligente al file con cui Nextflow può lavorare

    Pensatela così: una stringa di percorso è come scrivere un indirizzo su carta, mentre un oggetto Path è come avere l'indirizzo caricato in un dispositivo GPS che sa come navigare fino a lì e può dirvi dettagli sul viaggio.

### 1.3. Accedere agli attributi del file

Perché questo è utile? Beh, ora che Nextflow comprende che `myFile` è un oggetto Path e non solo una stringa, possiamo accedere ai vari attributi dell'oggetto Path.

Aggiorniamo il nostro flusso di lavoro per stampare gli attributi del file integrati:

=== "Dopo"

    ```groovy title="main.nf" linenums="5" hl_lines="4-9"
        // Create a Path object from a string path
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Print file attributes
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    ```

=== "Prima"

    ```groovy title="main.nf" linenums="5" hl_lines="4"
        // Create a Path object from a string path
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        println "${myFile} is of class ${myFile.class}"
    ```

Eseguite il flusso di lavoro:

```bash
nextflow run main.nf
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [ecstatic_ampere] DSL2 - revision: f3fa3dcb48

    File object class: sun.nio.fs.UnixPath
    File name: patientA_rep1_normal_R1_001.fastq.gz
    Simple name: patientA_rep1_normal_R1_001
    Extension: gz
    Parent directory: /workspaces/training/side-quests/working_with_files/data
    ```

Vedete i vari attributi del file stampati sulla console sopra.

### 1.4. Fornire il file a un processo

La differenza tra stringhe e oggetti Path diventa critica quando iniziate a costruire flussi di lavoro effettivi con processi.
Finora abbiamo verificato che Nextflow sta ora trattando il nostro file di input come un file, ma vediamo se possiamo effettivamente eseguire qualcosa su quel file in un processo.

#### 1.4.1. Importare il processo ed esaminare il codice

Vi forniamo un modulo di processo pre-scritto chiamato `COUNT_LINES` che prende un file di input e conta quante righe contiene.

Per utilizzare il processo nel flusso di lavoro, dovete solo aggiungere un'istruzione include prima del blocco workflow:

=== "Dopo"

    ```groovy title="main.nf" linenums="1" hl_lines="3"
    #!/usr/bin/env nextflow

    include { COUNT_LINES } from './modules/count_lines.nf'

    workflow {
    ```

=== "Prima"

    ```groovy title="main.nf" linenums="1"
    #!/usr/bin/env nextflow

    workflow {
    ```

Potete aprire il file del modulo per esaminarne il codice:

```groovy title="modules/count_lines.nf" linenums="1"
#!/usr/bin/env nextflow

process COUNT_LINES {
    debug true

    input:
    path input_file

    script:
    """
    set -o pipefail
    echo "Processing file: $input_file"
    gzip -dc $input_file | wc -l
    """
}
```

Come potete vedere, è uno script abbastanza semplice che decomprime il file e conta quante righe contiene.

??? info "Cosa fa `debug true`?"

    La direttiva `debug true` nella definizione del processo fa sì che Nextflow stampi l'output dal vostro script (come il conteggio delle righe "40") direttamente nel log di esecuzione.
    Senza questo, vedreste solo lo stato di esecuzione del processo ma non l'output effettivo dal vostro script.

    Per maggiori informazioni sul debug dei flussi di lavoro Nextflow, consultate la side quest [Debugging Nextflow Workflows](debugging.md).

#### 1.4.2. Aggiungere una chiamata a `COUNT_LINES`

Ora che il processo è disponibile per il flusso di lavoro, possiamo aggiungere una chiamata al processo `COUNT_LINES` per eseguirlo sul file di input.

Apportate le seguenti modifiche al flusso di lavoro:

=== "Dopo"

    ```groovy title="main.nf" linenums="7" hl_lines="11-12"
        // Create a Path object from a string path
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Print file attributes
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"

        // Count the lines in the file
        COUNT_LINES(myFile)
    ```

=== "Prima"

    ```groovy title="main.nf" linenums="7" hl_lines="4-9"
        // Create a Path object from a string path
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Print file attributes
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    ```

E ora eseguite il flusso di lavoro:

```bash
nextflow run main.nf
```

??? success "Command output"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [cheeky_hypatia] DSL2 - revision: 281d13c414

    File object class: class sun.nio.fs.UnixPath
    File name: patientA_rep1_normal_R1_001.fastq.gz
    Simple name: patientA_rep1_normal_R1_001
    Extension: gz
    Parent directory: /workspaces/training/side-quests/working_with_files/data
    executor >  local (1)
    [e9/341c05] COUNT_LINES [100%] 1 of 1 ✔
    Processing file: /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz
    40
    ```

Questo mostra che siamo in grado di operare sul file in modo appropriato all'interno di un processo.

Nello specifico, Nextflow ha eseguito con successo le seguenti operazioni:

- Ha messo in staging il file nella directory di lavoro
- Ha decompresso il file .gz
- Ha contato le righe (40 righe in questo caso)
- È stato completato senza errori

La chiave di questa operazione fluida è che stiamo dicendo esplicitamente a Nextflow che il nostro input è un file e dovrebbe essere trattato come tale.

### 1.5. Risolvere errori di base nell'input dei file

Questo spesso confonde i nuovi arrivati a Nextflow, quindi prendiamoci qualche minuto per vedere cosa succede quando lo fate in modo sbagliato.

Ci sono due punti principali in cui potete sbagliare la gestione dei file: a livello del flusso di lavoro e a livello del processo.

#### 1.5.1. Errore a livello di flusso di lavoro

Vediamo cosa succede se torniamo a trattare il file come una stringa quando specifichiamo l'input nel blocco workflow.

Apportate le seguenti modifiche al flusso di lavoro, assicurandovi di commentare le istruzioni di stampa specifiche del percorso:

=== "Dopo"

    ```groovy title="main.nf" linenums="7" hl_lines="2 6-11"
        // Create a Path object from a string path
        myFile = 'data/patientA_rep1_normal_R1_001.fastq.gz'

        // Print file attributes
        println "File object class: ${myFile.class}"
        /*
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
        */

        // Count the lines in the file
        COUNT_LINES(myFile)
    ```

=== "Prima"

    ```groovy title="main.nf" linenums="7" hl_lines="4-9"
        // Create a Path object from a string path
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Print file attributes
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"

        // Count the lines in the file
        COUNT_LINES(myFile)
    ```

E ora eseguite il flusso di lavoro:

```bash
nextflow run main.nf
```

??? failure "Command output"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [friendly_goodall] DSL2 - revision: ae50609b20

    [-        ] COUNT_LINES -
    ERROR ~ Error executing process > 'COUNT_LINES'

    Caused by:
      Not a valid path value: 'data/patientA_rep1_normal_R1_001.fastq.gz'



    Tip: view the complete command output by changing to the process work dir and entering the command `cat .command.out`

    -- Check '.nextflow.log' file for details
    ```

Questa è la parte importante:

```console
Not a valid path value: 'data/patientA_rep1_normal_R1_001.fastq.gz'
```

Quando specificate un input `path`, Nextflow verifica che stiate passando riferimenti a file effettivi, non solo stringhe.
Questo errore vi sta dicendo che `'data/patientA_rep1_normal_R1_001.fastq.gz'` non è un valore di percorso valido perché è una stringa, non un oggetto Path.

Nextflow ha rilevato immediatamente il problema e si è fermato prima ancora di avviare il processo.

#### 1.5.2. Errore a livello di processo

L'altro posto in cui potremmo dimenticare di specificare che vogliamo che Nextflow tratti l'input come un file è nella definizione del processo.

!!! warning "Mantenere l'errore del flusso di lavoro dalla sezione 1.5.1"

    Perché questo test funzioni correttamente, mantenete il flusso di lavoro nel suo stato errato (usando una semplice stringa invece di `file()`).
    Quando combinato con `val` nel processo, questo produce l'errore mostrato di seguito.

Apportate la seguente modifica al modulo:

=== "Dopo"

    ```groovy title="modules/count_lines.nf" linenums="3" hl_lines="5"
    process COUNT_LINES {
        debug true

        input:
        val input_file
    ```

=== "Prima"

    ```groovy title="modules/count_lines.nf" linenums="3" hl_lines="5"
    process COUNT_LINES {
        debug true

        input:
        path input_file
    ```

E ora eseguite nuovamente il flusso di lavoro:

```bash
nextflow run main.nf
```

??? failure "Command output"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [soggy_golick] DSL2 - revision: ae50609b20

    executor >  local (1)
    [b3/b3023c] COUNT_LINES [  0%] 0 of 1 ✘
    ERROR ~ Error executing process > 'COUNT_LINES'

    Caused by:
      Process `COUNT_LINES` terminated with an error exit status (1)


    Command executed:

      set -o pipefail
      echo "Processing file: data/patientA_rep1_normal_R1_001.fastq.gz"
      gzip -dc data/patientA_rep1_normal_R1_001.fastq.gz | wc -l

    Command exit status:
      1

    Command output:
      Processing file: data/patientA_rep1_normal_R1_001.fastq.gz
      0

    Command error:
      Processing file: data/patientA_rep1_normal_R1_001.fastq.gz
      gzip: data/patientA_rep1_normal_R1_001.fastq.gz: No such file or directory
      0

    Work dir:
      /workspaces/training/side-quests/working_with_files/work/b3/b3023cb2ccb986851301d8e369e79f

    Tip: when you have fixed the problem you can continue the execution adding the option `-resume` to the run command line

    -- Check '.nextflow.log' file for details
    ```

Questo mostra molti dettagli sull'errore perché il processo è impostato per produrre informazioni di debug, come notato sopra.

Queste sono le sezioni più rilevanti:

```console
Command executed:

  set -o pipefail
  echo "Processing file: data/patientA_rep1_normal_R1_001.fastq.gz"
  gzip -dc data/patientA_rep1_normal_R1_001.fastq.gz | wc -l
```

```console
Command error:
  Processing file: data/patientA_rep1_normal_R1_001.fastq.gz
  gzip: data/patientA_rep1_normal_R1_001.fastq.gz: No such file or directory
  0
```

Questo dice che il sistema non è riuscito a trovare il file; tuttavia se cercate il percorso, c'è un file con quel nome in quella posizione.

Quando abbiamo eseguito questo, Nextflow ha passato il valore stringa allo script, ma non ha _messo in staging_ il file effettivo nella directory di lavoro.
Quindi il processo ha cercato di utilizzare la stringa relativa, `data/patientA_rep1_normal_R1_001.fastq.gz`, ma quel file non esiste all'interno della directory di lavoro del processo.

Presi insieme, questi due esempi vi mostrano quanto sia importante dire a Nextflow se un input deve essere gestito come un file.

!!! note

    Assicuratevi di tornare indietro e correggere entrambi gli errori intenzionali prima di continuare con la prossima sezione.

### Takeaway

- Stringhe di percorso vs oggetti Path: Le stringhe sono solo testo, gli oggetti Path sono riferimenti intelligenti ai file
- Il metodo `file()` converte una stringa di percorso in un oggetto Path con cui Nextflow può lavorare
- Potete accedere alle proprietà del file come `name`, `simpleName`, `extension` e `parent` [utilizzando gli attributi del file](https://www.nextflow.io/docs/latest/working-with-files.html#getting-file-attributes)
- L'uso di oggetti Path invece di stringhe consente a Nextflow di gestire correttamente i file nel vostro flusso di lavoro
- Risultati dell'Input del Processo: Una corretta gestione dei file richiede oggetti Path, non stringhe, per garantire che i file siano correttamente messi in staging e accessibili per l'uso da parte dei processi.

---

## 2. Utilizzo di file remoti

Una delle caratteristiche chiave di Nextflow è la capacità di passare senza problemi tra file locali (sulla stessa macchina) e file remoti accessibili via internet.

Se lo fate correttamente, non dovreste mai aver bisogno di cambiare la logica del vostro flusso di lavoro per accogliere file provenienti da posizioni diverse.
Tutto ciò che dovete fare per utilizzare un file remoto è specificare il prefisso appropriato nel percorso del file quando lo fornite al flusso di lavoro.

Ad esempio, `/path/to/data` non ha prefisso, indicando che è un percorso di file locale 'normale', mentre `s3://path/to/data` include il prefisso `s3://`, indicando che si trova nell'object storage S3 di Amazon.

Sono supportati molti protocolli diversi:

- HTTP(S)/FTP (http://, https://, ftp://)
- Amazon S3 (s3://)
- Azure Blob Storage (az://)
- Google Cloud Storage (gs://)

Per utilizzare uno qualsiasi di questi, specificate semplicemente il prefisso rilevante nella stringa, che viene quindi tecnicamente chiamata Uniform Resource Identifier (URI) invece di percorso di file.
Nextflow gestirà l'autenticazione e il posizionamento dei file nel posto giusto, scaricando o caricando e tutte le altre operazioni sui file che vi aspettereste.

Il punto di forza chiave di questo sistema è che ci consente di passare tra ambienti senza cambiare alcuna logica della pipeline.
Ad esempio, potete sviluppare con un piccolo set di test locale prima di passare a un set di test su larga scala situato in storage remoto semplicemente cambiando l'URI.

### 2.1. Utilizzare un file da internet

Testiamo questo sostituendo il percorso locale che stiamo fornendo al nostro flusso di lavoro con un percorso HTTPS che punta a una copia degli stessi dati memorizzati in Github.

!!! warning

    Questo funzionerà solo se avete una connessione internet attiva.

Aprite nuovamente `main.nf` e cambiate il percorso di input come segue:

=== "Dopo"

    ```groovy title="main.nf" linenums="2" hl_lines="2"
        // Using a remote file from the internet
        myFile = file('https://raw.github.com/nextflow-io/training/master/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz')

        // Print file attributes
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    ```

=== "Prima"

    ```groovy title="main.nf" linenums="2" hl_lines="2"
        // Create a Path object from a string path
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Print file attributes
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    ```

Eseguiamo il flusso di lavoro:

```bash
nextflow run main.nf
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [insane_swartz] DSL2 - revision: fff18abe6d

    File object class: class nextflow.file.http.XPath
    File name: patientA_rep1_normal_R1_001.fastq.gz
    Simple name: patientA_rep1_normal_R1_001
    Extension: gz
    Parent directory: /nextflow-io/training/master/side-quests/working_with_files/data
    executor >  local (1)
    [8a/2ab7ca] COUNT_LINES [100%] 1 of 1 ✔
    Processing file: patientA_rep1_normal_R1_001.fastq.gz
    40
    ```

Funziona! Potete vedere che è cambiato molto poco.

L'unica differenza nell'output della console è che la classe dell'oggetto path è ora `nextflow.file.http.XPath`, mentre per il percorso locale la classe era `sun.nio.fs.UnixPath`.
Non dovete ricordare queste classi; lo menzioniamo solo per dimostrare che Nextflow identifica e gestisce le diverse posizioni in modo appropriato.

Dietro le quinte, Nextflow ha scaricato il file in una directory di staging situata all'interno della directory di lavoro.
Quel file in staging può quindi essere trattato come un file locale e collegato simbolicamente nella directory del processo rilevante.

Potete verificare che ciò sia accaduto qui guardando il contenuto della directory di lavoro situata al valore hash del processo.

??? abstract "Work directory contents"

    Se l'hash del processo era `8a/2ab7ca`, potreste esplorare la directory di lavoro:

    ```console
    $ ls -la work/8a/2ab7ca*/
    total 16
    drwxr-xr-x  6 user  staff   192 Jan 28 10:00 .
    drwxr-xr-x  3 user  staff    96 Jan 28 10:00 ..
    -rw-r--r--  1 user  staff     0 Jan 28 10:00 .command.begin
    -rw-r--r--  1 user  staff   127 Jan 28 10:00 .command.sh
    lrwxr-xr-x  1 user  staff    89 Jan 28 10:00 patientA_rep1_normal_R1_001.fastq.gz -> /path/to/work/stage/.../patientA_rep1_normal_R1_001.fastq.gz
    ```

    Il collegamento simbolico punta a una copia in staging del file remoto che Nextflow ha scaricato automaticamente.

Notate che per file più grandi, il passaggio di download richiederà tempo extra rispetto all'esecuzione su file locali.
Tuttavia, Nextflow verifica se ha già una copia in staging per evitare download non necessari.
Quindi se eseguite nuovamente sullo stesso file e non avete eliminato il file in staging, Nextflow utilizzerà la copia in staging.

Questo mostra quanto sia facile passare tra dati locali e remoti utilizzando Nextflow, che è una caratteristica chiave di Nextflow.

!!! note

    L'unica eccezione importante a questo principio è che non potete utilizzare pattern glob o percorsi di directory con HTTPS perché HTTPS non può elencare più file, quindi dovete specificare URL di file esatti.
    Tuttavia, altri protocolli di storage come blob storage (`s3://`, `az://`, `gs://`) possono utilizzare sia glob che percorsi di directory.

    Ecco come potreste utilizzare pattern glob con cloud storage:

    ```groovy title="Cloud storage examples (not runnable in this environment)"
    // S3 with glob patterns - would match multiple files
    ch_s3_files = channel.fromPath('s3://my-bucket/data/*.fastq.gz')

    // Azure Blob Storage with glob patterns
    ch_azure_files = channel.fromPath('az://container/data/patient*_R{1,2}.fastq.gz')

    // Google Cloud Storage with glob patterns
    ch_gcs_files = channel.fromPath('gs://bucket/data/sample_*.fastq.gz')
    ```

    Vi mostreremo come lavorare con i glob in pratica nella prossima sezione.

### 2.2. Tornare al file locale

Torneremo a utilizzare i nostri file di esempio locali per il resto di questa side quest, quindi riportiamo l'input del flusso di lavoro al file originale:

=== "Dopo"

    ```groovy title="main.nf" linenums="2" hl_lines="2"
        // Create a Path object from a string path
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Print file attributes
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    ```

=== "Prima"

    ```groovy title="main.nf" linenums="2" hl_lines="2"
        // Create a Path object from a string path
        myFile = file('https://raw.github.com/nextflow-io/training/master/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz')

        // Print file attributes
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    ```

### Takeaway

- I dati remoti sono accessibili utilizzando un URI (HTTP, FTP, S3, Azure, Google Cloud)
- Nextflow scaricherà e metterà automaticamente in staging i dati nel posto giusto, purché questi percorsi vengano forniti ai processi
- Non scrivete logica per scaricare o caricare file remoti!
- I file locali e remoti producono tipi di oggetti diversi ma funzionano in modo identico
- **Importante**: HTTP/HTTPS funzionano solo con singoli file (nessun pattern glob)
- Il cloud storage (S3, Azure, GCS) supporta sia singoli file che pattern glob
- Potete passare senza problemi tra sorgenti di dati locali e remote senza cambiare la logica del codice (purché il protocollo supporti le operazioni richieste)

---

## 3. Utilizzo della fabbrica di canali `fromPath()`

Finora abbiamo lavorato con un singolo file alla volta, ma in Nextflow, tipicamente vorremo creare un canale di input con più file di input da elaborare.

Un modo ingenuo per farlo sarebbe combinare il metodo `file()` con [`channel.of()`](https://www.nextflow.io/docs/latest/reference/channel.html#of) in questo modo:

```groovy title="Syntax example"
ch_files = channel.of([file('data/patientA_rep1_normal_R1_001.fastq.gz')],
                      [file('data/patientA_rep1_normal_R1_001.fastq.gz')])
```

Funziona, ma è goffo.

!!! tip "Quando usare `file()` vs `channel.fromPath()`"

    - Usate `file()` quando avete bisogno di un singolo oggetto Path per manipolazione diretta (controllare se un file esiste, leggere i suoi attributi, o passarlo a una singola invocazione di processo)
    - Usate `channel.fromPath()` quando avete bisogno di un canale che può contenere più file, specialmente con pattern glob, o quando i file fluiranno attraverso più processi

È qui che entra in gioco [`channel.fromPath()`](https://www.nextflow.io/docs/latest/reference/channel.html#frompath): una comoda fabbrica di canali che raggruppa tutte le funzionalità di cui abbiamo bisogno per generare un canale da una o più stringhe di file statici così come pattern glob.

### 3.1. Aggiungere la fabbrica di canali

Aggiorniamo il nostro flusso di lavoro per utilizzare `channel.fromPath`.

=== "Dopo"

    ```groovy title="main.nf" linenums="7" hl_lines="1-3"
        // Load files with channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R1_001.fastq.gz')
        ch_files.view { myFile -> "Found file: $myFile" }

        // Print file attributes
        /* Commentateli per ora, ci torneremo!
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
        */

        // Count the lines in the file
        // COUNT_LINES(myFile)
    ```

=== "Prima"

    ```groovy title="main.nf" linenums="7" hl_lines="1-2"
        // Create a Path object from a string path
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Print file attributes
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"

        // Count the lines in the file
        COUNT_LINES(myFile)
    ```

Abbiamo anche commentato il codice che stampa gli attributi per ora, e aggiunto un'istruzione `.view` per stampare solo il nome del file invece.

Eseguite il flusso di lavoro:

```bash
nextflow run main.nf
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [grave_meucci] DSL2 - revision: b09964a583

    Found file: /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz
    ```

Come potete vedere, il percorso del file viene caricato come oggetto di tipo `Path` nel canale.
Questo è simile a ciò che avrebbe fatto `file()`, tranne che ora abbiamo un canale in cui possiamo caricare più file se vogliamo.

L'uso di `channel.fromPath()` è un modo conveniente per creare un nuovo canale popolato da un elenco di file.

### 3.2. Visualizzare gli attributi dei file nel canale

Nel nostro primo tentativo di utilizzare la fabbrica di canali, abbiamo semplificato il codice e stampato solo il nome del file.

Torniamo a stampare gli attributi completi del file:

=== "Dopo"

    ```groovy title="main.nf" linenums="7" hl_lines="3-9 12"
        // Load files with channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R1_001.fastq.gz')
        ch_files.view { myFile ->
            println "File object class: ${myFile.class}"
            println "File name: ${myFile.name}"
            println "Simple name: ${myFile.simpleName}"
            println "Extension: ${myFile.extension}"
            println "Parent directory: ${myFile.parent}"
        }

        // Count the lines in the file
        COUNT_LINES(ch_files)
    ```

=== "Prima"

    ```groovy title="main.nf" linenums="7" hl_lines="3"
        // Load files with channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R1_001.fastq.gz')
        ch_files.view { myFile -> "Found file: $myFile" }

        // Count the lines in the file
        // COUNT_LINES(ch_files)
    ```

Stiamo anche riabilitando la chiamata al processo `COUNT_LINES` per verificare che l'elaborazione dei file funzioni ancora correttamente con il nostro approccio basato sui canali.

Eseguite il flusso di lavoro:

```bash
nextflow run main.nf
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [furious_swanson] DSL2 - revision: c35c34950d

    executor >  local (1)
    [9d/6701a6] COUNT_LINES (1) [100%] 1 of 1 ✔
    File object class: sun.nio.fs.UnixPath
    File name: patientA_rep1_normal_R1_001.fastq.gz
    Simple name: patientA_rep1_normal_R1_001
    Extension: gz
    Parent directory: /workspaces/training/side-quests/working_with_files/data
    Processing file: patientA_rep1_normal_R1_001.fastq.gz
    40
    ```

Ed eccolo, stessi risultati di prima ma ora abbiamo il file in un canale, quindi possiamo aggiungerne altri.

### 3.3. Utilizzo di un glob per corrispondere a più file

Ci sono diversi modi in cui potremmo caricare più file nel canale.
Qui vi mostreremo come utilizzare i pattern glob, che sono un modo conveniente per corrispondere e recuperare nomi di file e directory basati su caratteri jolly.
Il processo di corrispondenza di questi pattern è chiamato "globbing" o "espansione del nome file".

!!! note

    Come notato in precedenza, Nextflow supporta il globbing per gestire file di input e output nella maggior parte dei casi, tranne con i percorsi di file HTTPS perché HTTPS non può elencare più file.

Diciamo che vogliamo recuperare entrambi i file in una coppia di file associati a un dato paziente, `patientA`:

```console
patientA_rep1_normal_R1_001.fastq.gz
patientA_rep1_normal_R2_001.fastq.gz
```

Poiché l'unica differenza tra i nomi dei file è il numero di replicato, _cioè_ il numero dopo `R`, possiamo utilizzare il carattere jolly `*` per sostituire il numero come segue:

```console
patientA_rep1_normal_R*_001.fastq.gz
```

Questo è il pattern glob di cui abbiamo bisogno.

Ora tutto ciò che dobbiamo fare è aggiornare il percorso del file nella fabbrica di canali per utilizzare quel pattern glob come segue:

=== "Dopo"

    ```groovy title="main.nf" linenums="7"
      // Load files with channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
    ```

=== "Prima"

    ```groovy title="main.nf" linenums="7"
      // Load files with channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R1_001.fastq.gz')
    ```

Nextflow riconoscerà automaticamente che questo è un pattern glob e lo gestirà in modo appropriato.

Eseguite il flusso di lavoro per testarlo:

```bash
nextflow run main.nf
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [boring_sammet] DSL2 - revision: d2aa789c9a

    executor >  local (2)
    [3c/a65de5] COUNT_LINES (2) [100%] 2 of 2 ✔
    File object class: class sun.nio.fs.UnixPath
    File name: patientA_rep1_normal_R1_001.fastq.gz
    Simple name: patientA_rep1_normal_R1_001
    Extension: gz
    Parent directory: /workspaces/training/side-quests/working_with_files/data
    File object class: class sun.nio.fs.UnixPath
    File name: patientA_rep1_normal_R2_001.fastq.gz
    Simple name: patientA_rep1_normal_R2_001
    Extension: gz
    Parent directory: /workspaces/training/side-quests/working_with_files/data
    Processing file: patientA_rep1_normal_R1_001.fastq.gz
    40

    Processing file: patientA_rep1_normal_R2_001.fastq.gz
    40
    ```

Come potete vedere, ora abbiamo due oggetti Path nel nostro canale, il che mostra che Nextflow ha eseguito correttamente l'espansione del nome file, e ha caricato ed elaborato entrambi i file come previsto.

Utilizzando questo metodo, possiamo recuperare quanti file vogliamo semplicemente cambiando il pattern glob. Se lo rendessimo più generoso, ad esempio sostituendo tutte le parti variabili dei nomi dei file con `*` (_ad es._ `data/patient*_rep*_*_R*_001.fastq.gz`) potremmo prendere tutti i file di esempio nella directory `data`.

### Takeaway

- `channel.fromPath()` crea un canale con file che corrispondono a un pattern
- Ogni file viene emesso come elemento separato nel canale
- Possiamo utilizzare un pattern glob per corrispondere a più file
- I file vengono automaticamente convertiti in oggetti Path con attributi completi
- Il metodo `.view()` consente l'ispezione del contenuto del canale

---

## 4. Estrazione di metadati di base dai nomi dei file

Nella maggior parte dei domini scientifici, è molto comune avere metadati codificati nei nomi dei file che contengono i dati.
Ad esempio, in bioinformatica, i file contenenti dati di sequenziamento sono spesso denominati in modo da codificare informazioni sul campione, condizione, replicato e numero di lettura.

Se i nomi dei file sono costruiti secondo una convenzione coerente, potete estrarre quei metadati in modo standardizzato e utilizzarli nel corso della vostra analisi.
Questo è un grande 'se', ovviamente, e dovreste essere molto cauti ogni volta che vi affidate alla struttura del nome del file; ma la realtà è che questo approccio è molto ampiamente utilizzato, quindi diamo un'occhiata a come si fa in Nextflow.

Nel caso dei nostri dati di esempio, sappiamo che i nomi dei file includono metadati strutturati in modo coerente.
Ad esempio, il nome del file `patientA_rep1_normal_R2_001` codifica quanto segue:

- ID paziente: `patientA`
- ID replicato: `rep1`
- tipo di campione: `normal` (in contrapposizione a `tumor`)
- set di letture: `R1` (in contrapposizione a `R2`)

Modificheremo il nostro flusso di lavoro per recuperare queste informazioni in tre passaggi:

1. Recuperare il `simpleName` del file, che include i metadati
2. Separare i metadati utilizzando un metodo chiamato `tokenize()`
3. Utilizzare una map per organizzare i metadati

!!! warning

    Non dovreste mai codificare informazioni sensibili nei nomi dei file, come nomi di pazienti o altre caratteristiche identificative, poiché ciò può compromettere la privacy dei pazienti o altre restrizioni di sicurezza rilevanti.

### 4.1. Recuperare il `simpleName`

Il `simpleName` è un attributo del file che corrisponde al nome del file privato del suo percorso e dell'estensione.

Apportate le seguenti modifiche al flusso di lavoro:

=== "Dopo"

    ```groovy title="main.nf" linenums="7" hl_lines="3-6"
        // Load files with channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
        ch_files.map { myFile ->
            [ myFile.simpleName, myFile ]
        }
        .view()
    ```

=== "Prima"

    ```groovy title="main.nf" linenums="7" hl_lines="3-9"
        // Load files with channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
        ch_files.view { myFile ->
            println "File object class: ${myFile.class}"
            println "File name: ${myFile.name}"
            println "Simple name: ${myFile.simpleName}"
            println "Extension: ${myFile.extension}"
            println "Parent directory: ${myFile.parent}"
        }
    ```

Questo recupera il `simpleName` e lo associa all'oggetto file completo utilizzando un'operazione `map()`.

Eseguite il flusso di lavoro per testare che funzioni:

```bash
nextflow run main.nf
```

??? success "Command output"

    ```console hl_lines="7-8"
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [suspicious_mahavira] DSL2 - revision: ae8edc4e48

    executor >  local (2)
    [e9/55774b] COUNT_LINES (2) [100%] 2 of 2 ✔
    [patientA_rep1_normal_R2_001, /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R2_001.fastq.gz]
    [patientA_rep1_normal_R1_001, /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz]
    Processing file: patientA_rep1_normal_R1_001.fastq.gz
    40

    Processing file: patientA_rep1_normal_R2_001.fastq.gz
    40
    ```

Ogni elemento nel canale è ora una tupla contenente il `simpleName` e l'oggetto file originale.

### 4.2. Estrarre i metadati dal `simplename`

A questo punto, i metadati che vogliamo sono incorporati nel `simplename`, ma non possiamo accedere direttamente ai singoli elementi.
Quindi dobbiamo dividere il `simplename` nei suoi componenti.
Fortunatamente, quei componenti sono semplicemente separati da underscore nel nome del file originale, quindi possiamo applicare un metodo comune di Nextflow chiamato `tokenize()` che è perfetto per questo compito.

Apportate le seguenti modifiche al flusso di lavoro:

=== "Dopo"

    ```groovy title="main.nf" linenums="7" hl_lines="4"
        // Load files with channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
        ch_files.map { myFile ->
            [ myFile.simpleName.tokenize('_'), myFile ]
        }
    ```

=== "Prima"

    ```groovy title="main.nf" linenums="7" hl_lines="4"
        // Load files with channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
        ch_files.map { myFile ->
            [ myFile.simpleName, myFile ]
        }
    ```

Il metodo `tokenize()` dividerà la stringa `simpleName` ovunque trovi underscore, e restituirà un elenco contenente le sottostringhe.

Eseguite il flusso di lavoro:

```bash
nextflow run main.nf
```

??? success "Command output"

    ```console hl_lines="7-8"
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [gigantic_gauss] DSL2 - revision: a39baabb57

    executor >  local (2)
    [e7/da2f4b] COUNT_LINES (2) [100%] 2 of 2 ✔
    [[patientA, rep1, normal, R2, 001], /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R2_001.fastq.gz]
    [[patientA, rep1, normal, R1, 001], /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz]
    Processing file: patientA_rep1_normal_R2_001.fastq.gz
    40

    Processing file: patientA_rep1_normal_R1_001.fastq.gz
    40
    ```

Ora la tupla per ogni elemento nel nostro canale contiene l'elenco di metadati (_ad es._ `[patientA, rep1, normal, R1, 001]`) e l'oggetto file originale.

È fantastico!
Abbiamo scomposto le informazioni del nostro paziente da una singola stringa in un elenco di stringhe.
Ora possiamo gestire ogni parte delle informazioni del paziente separatamente.

### 4.3. Utilizzare una map per organizzare i metadati

I nostri metadati sono solo un elenco piatto al momento.
È abbastanza facile da usare ma difficile da leggere.

```console
[patientA, rep1, normal, R1, 001]
```

Qual è l'elemento all'indice 3? Potete dirlo senza fare riferimento alla spiegazione originale della struttura dei metadati?

Questa è una grande opportunità per utilizzare un archivio chiave-valore, dove ogni elemento ha un set di chiavi e i loro valori associati, quindi potete facilmente fare riferimento a ogni chiave per ottenere il valore corrispondente.

Nel nostro esempio, ciò significa passare da questa organizzazione:

```groovy
data = [patientA, 1, normal, R1]

println data[3]
```

A questa:

```groovy
data = [id: patientA, replicate: 1, type: normal, readNum: 1]

println data.readNum
```

In Nextflow, questo si chiama [map](https://nextflow.io/docs/latest/script.html#maps).

Convertiamo ora il nostro elenco piatto in una map.
Apportate le seguenti modifiche al flusso di lavoro:

=== "Dopo"

    ```groovy title="main.nf" linenums="7" hl_lines="4-13"
        // Load files with channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
        ch_files.map { myFile ->
            def (patient, replicate, type, readNum) = myFile.simpleName.tokenize('_')
            [
              [
                id: patient,
                replicate: replicate.replace('rep', ''),
                type: type,
                readNum: readNum.replace('R', ''),
              ],
              myFile
            ]
        }
    ```

=== "Prima"

    ```groovy title="main.nf" linenums="7" hl_lines="4"
        // Load files with channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
        ch_files.map { myFile ->
            [ myFile.simpleName.tokenize('_'), myFile ]
        }
    ```

Le modifiche chiave qui sono:

- **Assegnazione destrutturante**: `def (patient, replicate, type, readNum) = ...` estrae i valori tokenizzati in variabili nominate in una riga
- **Sintassi letterale della map**: `[id: patient, replicate: ...]` crea una map dove ogni chiave (come `id`) è associata a un valore (come `patient`)
- **Struttura annidata**: L'elenco esterno `[..., myFile]` accoppia la map dei metadati con l'oggetto file originale

Abbiamo anche semplificato un paio di stringhe di metadati utilizzando un metodo di sostituzione di stringhe chiamato `replace()` per rimuovere alcuni caratteri non necessari (_ad es._ `replicate.replace('rep', '')` per mantenere solo il numero dagli ID dei replicati).

Eseguiamo nuovamente il flusso di lavoro:

```bash
nextflow run main.nf
```

??? success "Command output"

    ```console hl_lines="7-8"
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [infallible_swartz] DSL2 - revision: 7f4e68c0cb

    executor >  local (2)
    [1b/e7fb27] COUNT_LINES (1) [100%] 2 of 2 ✔
    [[id:patientA, replicate:1, type:normal, readNum:2], /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R2_001.fastq.gz]
    [[id:patientA, replicate:1, type:normal, readNum:1], /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz]
    Processing file: patientA_rep1_normal_R2_001.fastq.gz
    40

    Processing file: patientA_rep1_normal_R1_001.fastq.gz
    40
    ```

Ora i metadati sono etichettati in modo ordinato (_ad es._ `[id:patientA, replicate:1, type:normal, readNum:2]`) quindi è molto più facile capire cosa è cosa.

Sarà anche molto più facile utilizzare effettivamente elementi di metadati nel flusso di lavoro, e renderà il nostro codice più facile da leggere e più manutenibile.

### Takeaway

- Possiamo gestire i nomi dei file in Nextflow con la potenza di un linguaggio di programmazione completo
- Possiamo trattare i nomi dei file come stringhe per estrarre informazioni rilevanti
- L'uso di metodi come `tokenize()` e `replace()` ci consente di manipolare stringhe nel nome del file
- L'operazione `.map()` trasforma gli elementi del canale preservando la struttura
- I metadati strutturati (map) rendono il codice più leggibile e manutenibile rispetto agli elenchi posizionali

Prossimamente, vedremo come gestire file di dati appaiati.

---

## 5. Gestione di file di dati appaiati

Molti disegni sperimentali producono file di dati appaiati che beneficiano di essere gestiti in modo esplicitamente appaiato.
Ad esempio, in bioinformatica, i dati di sequenziamento sono spesso generati sotto forma di letture appaiate, ovvero stringhe di sequenza che provengono dallo stesso frammento di DNA (spesso chiamate 'forward' e 'reverse' perché vengono lette da estremità opposte).

Questo è il caso dei nostri dati di esempio, dove R1 e R2 si riferiscono ai due set di letture.

```console
data/patientA_rep1_normal_R1_001.fastq.gz
data/patientA_rep1_normal_R2_001.fastq.gz
```

Nextflow fornisce una fabbrica di canali specializzata per lavorare con file appaiati come questo chiamata `channel.fromFilePairs()`, che raggruppa automaticamente i file in base a un pattern di denominazione condiviso. Ciò consente di associare i file appaiati più strettamente con meno sforzo.

Modificheremo il nostro flusso di lavoro per sfruttare questo.
Ci vorranno due passaggi:

1. Cambiare la fabbrica di canali in `channel.fromFilePairs()`
2. Estrarre e mappare i metadati

### 5.1. Cambiare la fabbrica di canali in `channel.fromFilePairs()`

Per utilizzare `channel.fromFilePairs`, dobbiamo specificare il pattern che Nextflow dovrebbe utilizzare per identificare i due membri in una coppia.

Tornando ai nostri dati di esempio, possiamo formalizzare il pattern di denominazione come segue:

```console
data/patientA_rep1_normal_R{1,2}_001.fastq.gz
```

Questo è simile al pattern glob che abbiamo usato prima, tranne che questo enumera specificamente le sottostringhe (sia `1` che `2` che vengono subito dopo la R) che identificano i due membri della coppia.

Aggiorniamo il flusso di lavoro `main.nf` di conseguenza:

=== "Dopo"

    ```groovy title="main.nf" linenums="7" hl_lines="1-2"
        // Load files with channel.fromFilePairs
        ch_files = channel.fromFilePairs('data/patientA_rep1_normal_R{1,2}_001.fastq.gz')
        /* Commentate il mapping per ora, ci torneremo!
        ch_files.map { myFile ->
            def (sample, replicate, type, readNum) = myFile.simpleName.tokenize('_')
            [
                [
                    id: sample,
                    replicate: replicate.replace('rep', ''),
                    type: type,
                    readNum: readNum,
                ],
                myFile
            ]
        }
        */
        .view()
    ```

=== "Prima"

    ```groovy title="main.nf" linenums="7" hl_lines="1-2"
        // Load files with channel.fromFilePairs
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
        ch_files.map { myFile ->
            def (sample, replicate, type, readNum) = myFile.simpleName.tokenize('_')
            [
                [
                    id: sample,
                    replicate: replicate.replace('rep', ''),
                    type: type,
                    readNum: readNum,
                ],
                myFile
            ]
        }
        .view()
    ```

Abbiamo cambiato la fabbrica di canali e adattato il pattern di corrispondenza dei file, e mentre ci eravamo, abbiamo commentato l'operazione map.
La aggiungeremo di nuovo più tardi, con alcune modifiche.

Eseguite il flusso di lavoro per testarlo:

```bash
nextflow run main.nf
```

??? failure "Command output"

    ```console hl_lines="7-8"
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [angry_koch] DSL2 - revision: 44fdf66105

    [-        ] COUNT_LINES -
    [-        ] COUNT_LINES -
    [patientA_rep1_normal_R, [/workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz, /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R2_001.fastq.gz]]
    ERROR ~ Error executing process > 'COUNT_LINES (1)'

    Caused by:
      Not a valid path value: 'patientA_rep1_normal_R'



    Tip: when you have fixed the problem you can continue the execution adding the option `-resume` to the run command line

    -- Check '.nextflow.log' file for details
    ```

Uh-oh, questa volta l'esecuzione è fallita!

La parte rilevante del messaggio di errore è qui:

```console
Not a valid path value: 'patientA_rep1_normal_R'
```

Questo perché abbiamo cambiato la fabbrica di canali.
Fino ad ora, il canale di input originale conteneva solo i percorsi dei file.
Tutta la manipolazione dei metadati che abbiamo fatto non ha effettivamente influenzato il contenuto del canale.

Ora che stiamo usando la fabbrica di canali `.fromFilePairs`, il contenuto del canale risultante è diverso.
Vediamo solo un elemento del canale, composto da una tupla contenente due elementi: la parte del `simpleName` condivisa dai due file, che funge da identificatore, e una tupla contenente i due oggetti file, nel formato `id, [ file1, file2 ]`.

È fantastico, perché Nextflow ha fatto il lavoro difficile di estrarre il nome del paziente esaminando il prefisso condiviso e usandolo come identificatore del paziente.

Tuttavia, questo rompe il nostro flusso di lavoro attuale.
Se volessimo ancora eseguire `COUNT_LINES` allo stesso modo senza cambiare il processo, dovremmo applicare un'operazione di mapping per estrarre i percorsi dei file.
Ma non lo faremo, perché il nostro obiettivo finale è utilizzare un processo diverso, `ANALYZE_READS`, che gestisce le coppie di file in modo appropriato.

Quindi commentiamo semplicemente (o eliminiamo) la chiamata a `COUNT_LINES` e andiamo avanti.

=== "Dopo"

    ```groovy title="main.nf" linenums="26" hl_lines="2"
        // Count the lines in the file
        // COUNT_LINES(ch_files)
    ```

=== "Prima"

    ```groovy title="main.nf" linenums="26" hl_lines="2"
        // Count the lines in the file
        COUNT_LINES(ch_files)
    ```

Potete anche commentare o eliminare l'istruzione include di `COUNT_LINES`, ma ciò non avrà alcun effetto funzionale.

Ora eseguiamo nuovamente il flusso di lavoro:

```bash
nextflow run main.nf
```

??? success "Command output"

    ```console hl_lines="5"
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [fabulous_davinci] DSL2 - revision: 22b53268dc

    [patientA_rep1_normal_R, [/workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz, /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R2_001.fastq.gz]]
    ```

Evviva, questa volta il flusso di lavoro ha successo!

Tuttavia, dobbiamo ancora estrarre il resto dei metadati dal campo `id`.

### 5.2. Estrarre e organizzare i metadati dalle coppie di file

La nostra operazione `map` di prima non funzionerà perché non corrisponde alla struttura dei dati, ma possiamo modificarla per farla funzionare.

Abbiamo già accesso all'identificatore effettivo del paziente nella stringa che `fromFilePairs()` ha usato come identificatore, quindi possiamo usarlo per estrarre i metadati senza ottenere il `simpleName` dall'oggetto Path come abbiamo fatto prima.

Decommentate l'operazione map nel flusso di lavoro e apportate le seguenti modifiche:

=== "Dopo"

    ```groovy title="main.nf" linenums="7" hl_lines="3-4 9 11 13"
        // Load files with channel.fromFilePairs
        ch_files = channel.fromFilePairs('data/patientA_rep1_normal_R{1,2}_001.fastq.gz')
        ch_files.map { id, files ->
            def (sample, replicate, type) = id.tokenize('_')
            [
                [
                    id: sample,
                    replicate: replicate.replace('rep', ''),
                    type: type
                ],
                files
            ]
        }
        .view()
    ```

=== "Prima"

    ```groovy title="main.nf" linenums="7" hl_lines="3-5 11 13"
        // Load files with channel.fromFilePairs
        ch_files = channel.fromFilePairs('data/patientA_rep1_normal_R{1,2}_001.fastq.gz')
        /* Commentate il mapping per ora, ci torneremo!
        ch_files.map { myFile ->
            def (sample, replicate, type, readNum) = myFile.simpleName.tokenize('_')
            [
                [
                    id: sample,
                    replicate: replicate.replace('rep', ''),
                    type: type,
                    readNum: readNum,
                ],
                myFile
            ]
        }
        */
        .view()
    ```

Questa volta la map parte da `id, files` invece che solo da `myFile`, e `tokenize()` viene applicato a `id` invece che a `myFile.simpleName`.

Notate anche che abbiamo eliminato `readNum` dalla riga `tokenize()`; qualsiasi sottostringa che non nominiamo specificamente (partendo da sinistra) verrà silenziosamente eliminata.
Possiamo farlo perché i file appaiati sono ora strettamente associati, quindi non abbiamo più bisogno di `readNum` nella map dei metadati.

Eseguiamo il flusso di lavoro:

```bash
nextflow run main.nf
```

??? success "Command output"

    ```console

    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [prickly_stonebraker] DSL2 - revision: f62ab10a3f

    [[id:patientA, replicate:1, type:normal], [/workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz, /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R2_001.fastq.gz]]
    ```

Ed eccolo: abbiamo la map dei metadati (`[id:patientA, replicate:1, type:normal]`) nella prima posizione della tupla di output, seguita dalla tupla di file appaiati, come previsto.

Naturalmente, questo raccoglierà ed elaborerà solo quella specifica coppia di file.
Se volete sperimentare con l'elaborazione di più coppie, potete provare ad aggiungere caratteri jolly nel pattern di input e vedere cosa succede.
Ad esempio, provate a usare `data/patientA_rep1_*_R{1,2}_001.fastq.gz`

### Takeaway

- [`channel.fromFilePairs()` trova e accoppia automaticamente file correlati](https://www.nextflow.io/docs/latest/reference/channel.html#fromfilepairs)
- Questo semplifica la gestione delle letture paired-end nella vostra pipeline
- I file appaiati possono essere raggruppati come tuple `[id, [file1, file2]]`
- L'estrazione dei metadati può essere fatta dall'ID del file appaiato piuttosto che dai singoli file

---

## 6. Utilizzo delle operazioni sui file nei processi

Ora mettiamo tutto insieme in un semplice processo per rafforzare come utilizzare le operazioni sui file all'interno di un processo Nextflow.

Vi forniamo un modulo di processo pre-scritto chiamato `ANALYZE_READS` che prende una tupla di metadati e una coppia di file di input e li analizza.
Potremmo immaginare che questo stia facendo allineamento di sequenze, o chiamata di varianti o qualsiasi altro passaggio che abbia senso per questo tipo di dati.

Iniziamo.

### 6.1. Importare il processo ed esaminare il codice

Per utilizzare questo processo nel flusso di lavoro, dobbiamo solo aggiungere un'istruzione include del modulo prima del blocco workflow.

Apportate la seguente modifica al flusso di lavoro:

=== "Dopo"

    ```groovy title="main.nf" linenums="1" hl_lines="3"
    #!/usr/bin/env nextflow

    include { ANALYZE_READS } from './modules/analyze_reads.nf'

    workflow {
    ```

=== "Prima"

    ```groovy title="main.nf" linenums="1"
    #!/usr/bin/env nextflow

    workflow {
    ```

Potete aprire il file del modulo per esaminarne il codice:

```groovy title="modules/analyze_reads.nf - process example" linenums="1"
#!/usr/bin/env nextflow

process ANALYZE_READS {
    tag { meta.id }

    publishDir { "results/${meta.id}" }, mode: 'copy'

    input:
    tuple val(meta), path(files)

    output:
    tuple val(meta.id), path("${meta.id}_stats.txt")

    script:
    """
    echo "Sample metadata: ${meta.id}" > ${meta.id}_stats.txt
    echo "Replicate: ${meta.replicate}" >> ${meta.id}_stats.txt
    echo "Type: ${meta.type}" >> ${meta.id}_stats.txt
    echo "Read 1: ${files[0]}" >> ${meta.id}_stats.txt
    echo "Read 2: ${files[1]}" >> ${meta.id}_stats.txt
    echo "Read 1 size: \$(gunzip -dc ${files[0]} | wc -l | awk '{print \$1/4}') reads" >> ${meta.id}_stats.txt
    echo "Read 2 size: \$(gunzip -dc ${files[1]} | wc -l | awk '{print \$1/4}') reads" >> ${meta.id}_stats.txt
    """
}
```

!!! note

    Le direttive `tag` e `publishDir` utilizzano la sintassi closure (`{ ... }`) invece dell'interpolazione di stringhe (`"${...}"`).
    Questo perché queste direttive fanno riferimento a variabili di input (`meta`) che non sono disponibili fino al runtime.
    La sintassi closure rinvia la valutazione fino a quando il processo non viene effettivamente eseguito.

!!! note

    Stiamo chiamando la nostra map di metadati `meta` per convenzione.
    Per un approfondimento sulle meta map, consultate la side quest [Metadata and meta maps](./metadata.md).

### 6.2. Chiamare il processo nel flusso di lavoro

Ora che il processo è disponibile per il flusso di lavoro, possiamo aggiungere una chiamata al processo `ANALYZE_READS` per eseguirlo.

Per eseguirlo sui nostri dati di esempio, dovremo fare due cose:

1. Dare un nome al canale rimappato
2. Aggiungere una chiamata al processo

#### 6.2.1. Nominare il canale di input rimappato

In precedenza abbiamo applicato le manipolazioni di mapping direttamente al canale di input.
Per fornire il contenuto rimappato al processo `ANALYZE_READS` (e farlo in modo chiaro e facile da leggere) vogliamo creare un nuovo canale chiamato `ch_samples`.

Possiamo farlo usando l'operatore [`set`](https://www.nextflow.io/docs/latest/reference/operator.html#set).

Nel flusso di lavoro principale, sostituite l'operatore `.view()` con `.set { ch_samples }`, e aggiungete una riga per testare che possiamo fare riferimento al canale per nome.

=== "Dopo"

    ```groovy title="main.nf" linenums="7" hl_lines="14 16-17"
        // Load files with channel.fromFilePairs
        ch_files = channel.fromFilePairs('data/patientA_rep1_normal_R{1,2}_001.fastq.gz')
        ch_files.map { id,  files ->
           def (sample, replicate, type, readNum) = id.tokenize('_')
           [
               [
                   id: sample,
                   replicate: replicate.replace('rep', ''),
                   type: type
               ],
               files
           ]
        }
            .set { ch_samples }

        // Temporary: peek into ch_samples
        ch_samples.view()
    ```

=== "Prima"

    ```groovy title="main.nf" linenums="7" hl_lines="14"
        // Load files with channel.fromFilePairs
        ch_files = channel.fromFilePairs('data/patientA_rep1_normal_R{1,2}_001.fastq.gz')
        ch_files.map { id,  files ->
           def (sample, replicate, type, readNum) = id.tokenize('_')
           [
               [
                   id: sample,
                   replicate: replicate.replace('rep', ''),
                   type: type
               ],
               files
           ]
        }
        .view()
    }
    ```

Eseguiamo questo:

```bash
nextflow run main.nf
```

??? success "Command output"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [goofy_kirch] DSL2 - revision: 3313283e42

    [[id:patientA, replicate:1, type:normal], [/workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz, /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R2_001.fastq.gz]]
    ```

Questo conferma che ora possiamo fare riferimento al canale per nome.

#### 6.2.2. Chiamare il processo sui dati

Ora chiamiamo effettivamente il processo `ANALYZE_READS` sul canale `ch_samples`.

Nel flusso di lavoro principale, apportate le seguenti modifiche al codice:

=== "Dopo"

    ```groovy title="main.nf" linenums="23"
        // Run the analysis
        ANALYZE_READS(ch_samples)
    ```

=== "Prima"

    ```groovy title="main.nf" linenums="23"
        // Temporary: peek into ch_samples
        ch_samples.view()
    ```

Eseguiamo questo:

```bash
nextflow run main.nf
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `./main.nf` [shrivelled_cori] DSL2 - revision: b546a31769

    executor >  local (1)
    [b5/110360] process > ANALYZE_READS (patientA) [100%] 1 of 1 ✔
    ```

Questo processo è impostato per pubblicare i suoi output in una directory `results`, quindi date un'occhiata lì.

??? abstract "Directory and file contents"

    ```console
    results
    └── patientA
        └── patientA_stats.txt
    ```

    ```txt title="patientA_stats.txt"
    Sample metadata: patientA
    Replicate: 1
    Type: normal
    Read 1: patientA_rep1_normal_R1_001.fastq.gz
    Read 2: patientA_rep1_normal_R2_001.fastq.gz
    Read 1 size: 10 reads
    Read 2 size: 10 reads
    ```

Il processo ha preso i nostri input e creato un nuovo file contenente i metadati del paziente, come progettato.
Splendido!

### 6.3. Includere molti più pazienti

Naturalmente, questo sta elaborando solo una singola coppia di file per un singolo paziente, che non è esattamente il tipo di alto throughput che sperate di ottenere con Nextflow.
Probabilmente vorrete elaborare molti più dati alla volta.

Ricordate che `channel.fromPath()` accetta un _glob_ come input, il che significa che può accettare qualsiasi numero di file che corrispondono al pattern.
Pertanto, se vogliamo includere tutti i pazienti, possiamo semplicemente modificare la stringa di input per includere più pazienti, come notato di sfuggita in precedenza.

Facciamo finta di voler essere il più avidi possibile.
Apportate le seguenti modifiche al flusso di lavoro:

=== "Dopo"

    ```groovy title="main.nf" linenums="7" hl_lines="2"
        // Load files with channel.fromFilePairs
        ch_files = channel.fromFilePairs('data/*_R{1,2}_001.fastq.gz')
    ```

=== "Prima"

    ```groovy title="main.nf" linenums="7" hl_lines="2"
        // Load files with channel.fromFilePairs
        ch_files = channel.fromFilePairs('data/patientA_rep1_normal_R{1,2}_001.fastq.gz')
    ```

Eseguite nuovamente la pipeline:

```bash
nextflow run main.nf
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `./main.nf` [big_stonebraker] DSL2 - revision: f7f9b8a76c

    executor >  local (8)
    [d5/441891] process > ANALYZE_READS (patientC) [100%] 8 of 8 ✔
    ```

La directory results dovrebbe ora contenere risultati per tutti i dati disponibili.

??? abstract "Directory contents"

    ```console
    results
    ├── patientA
    │   └── patientA_stats.txt
    ├── patientB
    │   └── patientB_stats.txt
    └── patientC
        └── patientC_stats.txt
    ```

Successo! Abbiamo analizzato tutti i pazienti in una volta sola! Giusto?

Forse no.
Se guardate più attentamente, abbiamo un problema: abbiamo due replicati per patientA, ma solo un file di output!
Stiamo sovrascrivendo il file di output ogni volta.

### 6.4. Rendere i file pubblicati unici

Poiché abbiamo accesso ai metadati del paziente, possiamo usarli per rendere i file pubblicati unici includendo metadati differenzianti, sia nella struttura della directory che nei nomi dei file stessi.

Apportate la seguente modifica al flusso di lavoro:

=== "Dopo"

    ```groovy title="modules/analyze_reads.nf" linenums="6"
        publishDir { "results/${meta.type}/${meta.id}/${meta.replicate}" }, mode: 'copy'
    ```

=== "Prima"

    ```groovy title="modules/analyze_reads.nf" linenums="6"
        publishDir { "results/${meta.id}" }, mode: 'copy'
    ```

Qui mostriamo l'opzione di utilizzare livelli di directory aggiuntivi per tenere conto dei tipi di campione e dei replicati, ma potreste sperimentare facendolo a livello di nome file.

Ora eseguite la pipeline un'ultima volta, ma assicuratevi di rimuovere prima la directory results per darvi uno spazio di lavoro pulito:

```bash
rm -r results
nextflow run main.nf
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `./main.nf` [insane_swartz] DSL2 - revision: fff18abe6d

    executor >  local (8)
    [e3/449081] process > ANALYZE_READS (patientC) [100%] 8 of 8 ✔
    ```

Controllate ora la directory results:

??? abstract "Directory contents"

    ```console
    results/
    ├── normal
    │   ├── patientA
    │   │   ├── 1
    │   │   │   └── patientA_stats.txt
    │   │   └── 2
    │   │       └── patientA_stats.txt
    │   ├── patientB
    │   │   └── 1
    │   │       └── patientB_stats.txt
    │   └── patientC
    │       └── 1
    │           └── patientC_stats.txt
    └── tumor
        ├── patientA
        │   ├── 1
        │   │   └── patientA_stats.txt
        │   └── 2
        │       └── patientA_stats.txt
        ├── patientB
        │   └── 1
        │       └── patientB_stats.txt
        └── patientC
            └── 1
                └── patientC_stats.txt
    ```

Ed eccolo, tutti i nostri metadati, ordinatamente organizzati. Questo è un successo!

C'è molto di più che potete fare una volta che avete i vostri metadati caricati in una map come questa:

1. Creare directory di output organizzate in base agli attributi del paziente
2. Prendere decisioni nei processi in base alle proprietà del paziente
3. Dividere, unire e ricombinare i dati in base ai valori dei metadati

Questo pattern di mantenere i metadati espliciti e collegati ai dati (piuttosto che codificati nei nomi dei file) è una best practice fondamentale in Nextflow che consente di costruire flussi di lavoro di analisi robusti e manutenibili.
Potete saperne di più su questo nella side quest [Metadata and meta maps](./metadata.md).

### Takeaway

- La direttiva `publishDir` può organizzare gli output in base ai valori dei metadati
- I metadati nelle tuple consentono un'organizzazione strutturata dei risultati
- Questo approccio crea flussi di lavoro manutenibili con una chiara provenienza dei dati
- I processi possono prendere tuple di metadati e file come input
- La direttiva `tag` fornisce l'identificazione del processo nei log di esecuzione
- La struttura del flusso di lavoro separa la creazione del canale dall'esecuzione del processo

---

## Riepilogo

In questa side quest, avete imparato come lavorare con i file in Nextflow, dalle operazioni di base alle tecniche più avanzate per gestire collezioni di file.

L'applicazione di queste tecniche nel vostro lavoro vi consentirà di costruire flussi di lavoro più efficienti e manutenibili, specialmente quando lavorate con un gran numero di file con convenzioni di denominazione complesse.

### Pattern chiave

1.  **Operazioni di Base sui File:** Abbiamo creato oggetti Path con `file()` e acceduto agli attributi del file come nome, estensione e directory padre, imparando la differenza tra stringhe e oggetti Path.

    - Creare un oggetto Path con `file()`

    ```groovy
    myFile = file('path/to/file.txt')
    ```

    - Ottenere gli attributi del file

    ```groovy
    println myFile.name       // file.txt
    println myFile.baseName   // file
    println myFile.extension  // txt
    println myFile.parent     // path/to
    ```

2.  **Utilizzo di File Remoti**: Abbiamo imparato come passare in modo trasparente tra file locali e remoti utilizzando URI, dimostrando la capacità di Nextflow di gestire file da varie sorgenti senza cambiare la logica del flusso di lavoro.

    - File locale

    ```groovy
    myFile = file('path/to/file.txt')
    ```

    - FTP

    ```groovy
    myFile = file('ftp://path/to/file.txt')
    ```

    - HTTPS

    ```groovy
    myFile = file('https://path/to/file.txt')
    ```

    - Amazon S3

    ```groovy
    myFile = file('s3://path/to/file.txt')
    ```

    - Azure Blob Storage

    ```groovy
    myFile = file('az://path/to/file.txt')
    ```

    - Google Cloud Storage

    ```groovy
    myFile = file('gs://path/to/file.txt')
    ```

3.  **Caricamento di file utilizzando la fabbrica di canali `fromPath()`:** Abbiamo creato canali da pattern di file con `channel.fromPath()` e visualizzato i loro attributi di file, inclusi i tipi di oggetto.

    - Creare un canale da un pattern di file

    ```groovy
     ch_files = channel.fromPath('data/*.fastq.gz')
    ```

    - Ottenere gli attributi del file

    ```groovy
     ch_files.view { myFile ->
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    }
    ```

4.  **Estrazione dei Metadati del Paziente dai Nomi dei File:** Abbiamo utilizzato `tokenize()` e `replace()` per estrarre e strutturare i metadati dai nomi dei file, convertendoli in map organizzate.

    ```groovy
    def name = file.name.tokenize('_')
    def patientId = name[0]
    def replicate = name[1].replace('rep', '')
    def type = name[2]
    def readNum = name[3].replace('R', '')
    ```

5.  **Semplificazione con channel.fromFilePairs:** Abbiamo utilizzato `channel.fromFilePairs()` per accoppiare automaticamente file correlati ed estrarre metadati dagli ID dei file appaiati.

    ```groovy
    ch_pairs = channel.fromFilePairs('data/*_R{1,2}_001.fastq.gz')
    ```

6.  **Utilizzo delle Operazioni sui File nei Processi:** Abbiamo integrato le operazioni sui file nei processi Nextflow con una corretta gestione degli input, utilizzando `publishDir` per organizzare gli output in base ai metadati.

    - Associare una meta map con gli input del processo

    ```groovy
    ch_files = channel.fromFilePairs('data/patientA_rep1_normal_R{1,2}_001.fastq.gz')
    ch_files.map { id,  files ->
        def (sample, replicate, type, readNum) = id.tokenize('_')
        [
            [
                id: sample,
                replicate: replicate.replace('rep', ''),
                type: type
            ],
             files
        ]
    }
        .set { ch_samples }

    ANALYZE_READS(ch_samples)
    ```

    - Organizzare gli output in base ai metadati

    ```groovy
    publishDir { "results/${meta.type}/${meta.id}/${meta.replicate}" }, mode: 'copy'
    ```

### Risorse aggiuntive

- [Documentazione Nextflow: Working with Files](https://www.nextflow.io/docs/latest/working-with-files.html)
- [channel.fromPath](https://www.nextflow.io/docs/latest/reference/channel.html#frompath)
- [channel.fromFilePairs](https://www.nextflow.io/docs/latest/reference/channel.html#fromfilepairs)

---

## Cosa c'è dopo?

Tornate al [menu delle Side Quest](./index.md) o cliccate il pulsante in basso a destra della pagina per passare al prossimo argomento nell'elenco.
