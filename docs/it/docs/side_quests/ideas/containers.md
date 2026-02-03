# Parte 1: Maggiori Informazioni sui Container

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduzione assistita da IA - [scopri di più e suggerisci miglioramenti](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

[TODO]

---

## 1. Come trovare o creare immagini di container

Alcuni sviluppatori di software forniscono immagini di container per i loro software disponibili su registri di container come Docker Hub, ma molti non lo fanno.
In questa sezione opzionale, le mostreremo due modi per ottenere un'immagine di container per gli strumenti che desidera utilizzare nelle sue pipeline Nextflow: utilizzando Seqera Containers e costruendo l'immagine del container da sé.

Otterrà/costruirà un'immagine di container per il pacchetto pip `quote`, che sarà utilizzato nell'esercizio alla fine di questa sezione.

### 1.1. Ottenere un'immagine di container da Seqera Containers

Seqera Containers è un servizio gratuito che costruisce immagini di container per strumenti installabili tramite pip e conda (incluso bioconda).
Navighi su [Seqera Containers](https://www.seqera.io/containers/) e cerchi il pacchetto pip `quote`.

![Seqera Containers](img/seqera-containers-1.png)

Clicchi su "+Add" e poi su "Get Container" per richiedere un'immagine di container per il pacchetto pip `quote`.

![Seqera Containers](img/seqera-containers-2.png)

Se è la prima volta che viene costruito un container di comunità per questa versione del pacchetto, potrebbero essere necessari alcuni minuti per completare l'operazione.
Clicchi per copiare l'URI (ad es. `community.wave.seqera.io/library/pip_quote:ae07804021465ee9`) dell'immagine di container che è stata creata per voi.

Ora potete utilizzare l'immagine di container per eseguire il comando `quote` e ottenere un detto casuale di Grace Hopper.

```bash
docker run --rm community.wave.seqera.io/library/pip_quote:ae07804021465ee9 quote "Grace Hopper"
```

Output:

```console title="Output"
Humans are allergic to change. They love to say, 'We've always done it
this way.' I try to fight that. That's why I have a clock on my wall
that runs counter-clockwise.
```

### 1.2. Costruire l'immagine del container da sé

Utilizziamo alcuni dettagli di costruzione dal sito web di Seqera Containers per costruire noi stessi l'immagine di container per il pacchetto pip `quote`.
Ritornate al sito web di Seqera Containers e cliccate sul pulsante "Build Details".

Il primo elemento che esamineremo è il `Dockerfile`, un tipo di file di script che contiene tutti i comandi necessari per costruire l'immagine di container.
Abbiamo aggiunto alcuni commenti esplicativi al Dockerfile qui sotto per aiutarla a comprendere cosa fa ciascuna parte.

```Dockerfile title="Dockerfile"
# Inizia dall'immagine docker di base micromamba
FROM mambaorg/micromamba:1.5.10-noble
# Copia il file conda.yml nel container
COPY --chown=$MAMBA_USER:$MAMBA_USER conda.yml /tmp/conda.yml
# Installa varie utilità per Nextflow da utilizzare e i pacchetti nel file conda.yml
RUN micromamba install -y -n base -f /tmp/conda.yml \
    && micromamba install -y -n base conda-forge::procps-ng \
    && micromamba env export --name base --explicit > environment.lock \
    && echo ">> CONDA_LOCK_START" \
    && cat environment.lock \
    && echo "<< CONDA_LOCK_END" \
    && micromamba clean -a -y
# Esegui il container come utente root
USER root
# Imposta la variabile d'ambiente PATH per includere la directory di installazione di micromamba
ENV PATH="$MAMBA_ROOT_PREFIX/bin:$PATH"
```

Il secondo elemento che esamineremo è il file `conda.yml`, che contiene l'elenco dei pacchetti che devono essere installati nell'immagine di container.

```conda.yml title="conda.yml"
channels:
- conda-forge
- bioconda
dependencies:
- pip
- pip:
  - quote==3.0.0 #
```

Copiate il contenuto di questi file negli stub situati nella directory `containers/build`, quindi eseguite il seguente comando per costruire l'immagine di container da sé.

!!! Note "Nota"

    Utilizziamo il flag `-t quote:latest` per etichettare l'immagine di container con il nome `quote` e il tag `latest`.
    Saremo in grado di utilizzare questo tag per riferirci all'immagine di container quando la eseguiremo su questo sistema.

```bash
docker build -t quote:latest containers/build
```

Dopo che la costruzione è terminata, potete eseguire l'immagine di container che avete appena costruito.

```bash
docker run --rm quote:latest quote "Margaret Oakley Dayhoff"
```

### Riepilogo

Ha imparato due modi diversi per ottenere un'immagine di container per uno strumento che desidera utilizzare nelle sue pipeline Nextflow: utilizzando Seqera Containers e costruendo l'immagine di container da sé.

### Qual è il prossimo passo?

Ha tutto ciò di cui ha bisogno per continuare al [prossimo capitolo](./04_hello_genomics.md) di questa serie di formazione.
Può anche continuare con un esercizio opzionale per recuperare citazioni di pionieri dell'informatica/biologia utilizzando il container `quote` e visualizzarle utilizzando il container `cowsay`.

---

## 2. Fare in modo che la mucca citi scienziati famosi

Questa sezione contiene alcuni esercizi avanzati, per praticare ciò che ha imparato finora.
Completare questi esercizi _non è richiesto_ per comprendere le parti successive della formazione, ma fornisce un modo divertente per rafforzare le sue conoscenze capendo come fare in modo che la mucca citi scienziati famosi.

```console title="cowsay-output-Grace-Hopper.txt"
  _________________________________________________
 /                                                 \
| Humans are allergic to change. They love to       |
| say, 'We've always done it this way.' I try to fi |
| ght that. That's why I have a clock on my wall th |
| at runs counter-clockwise.                        |
| -Grace Hopper                                     |
 \                                                 /
  =================================================
                                                 \
                                                  \
                                                    ^__^
                                                    (oo)\_______
                                                    (__)\       )\/\
                                                        ||----w |
                                                        ||     ||
```

### 2.1. Modificare lo script `hello-containers.nf` per utilizzare un processo getQuote

Abbiamo un elenco di pionieri dell'informatica e della biologia nel file `containers/data/pioneers.csv`.
Ad alto livello, per completare questo esercizio dovrà:

- Modificare il `params.input_file` predefinito per puntare al file `pioneers.csv`.
- Creare un processo `getQuote` che utilizza il container `quote` per recuperare una citazione per ogni input.
- Connettere l'output del processo `getQuote` al processo `cowsay` per visualizzare la citazione.

Per l'immagine di container `quote`, potete utilizzare quella che avete costruito da sé nell'esercizio avanzato precedente o utilizzare quella ottenuta da Seqera Containers.

!!! Hint "Suggerimento"

    Una buona scelta per il blocco `script` del suo processo getQuote potrebbe essere:
        ```groovy
        script:
            def safe_author = author.tokenize(' ').join('-')
            """
            quote "$author" > quote-${safe_author}.txt
            echo "-${author}" >> quote-${safe_author}.txt
            """
        ```

Può trovare una soluzione a questo esercizio in `containers/solutions/hello-containers-4.1.nf`.

### 2.2. Modificare la sua pipeline Nextflow per permetterle di eseguire nelle modalità `quote` e `sayHello`.

Aggiunga della logica di ramificazione alla sua pipeline per permetterle di accettare input destinati sia a `quote` che a `sayHello`.
Ecco un esempio di come utilizzare un'istruzione `if` in un workflow Nextflow:

```groovy title="hello-containers.nf"
workflow {
    if (params.quote) {
        ...
    }
    else {
        ...
    }
    cowSay(text_ch)
}
```

!!! Hint "Suggerimento"

    Può utilizzare `new_ch = processName.out` per assegnare un nome al canale di output di un processo.

Può trovare una soluzione a questo esercizio in `containers/solutions/hello-containers-4.2.nf`.

### Riepilogo

Sa come utilizzare i container in Nextflow per eseguire processi e come costruire della logica di ramificazione nelle sue pipeline!

### Qual è il prossimo passo?

Festeggi, faccia una pausa e beva dell'acqua!

Quando sarà pronto, passi alla Parte 3 di questa serie di formazione per imparare come applicare ciò che ha imparato finora a un caso d'uso di analisi dati più realistico.
