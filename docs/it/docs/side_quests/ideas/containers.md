# Parte 1: Altri Container

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduzione assistita da IA - [scopri di più e suggerisci miglioramenti](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

[TODO]

---

## 1. Come trovare o creare immagini container

Alcuni sviluppatori di software forniscono immagini container per i loro software disponibili su registri di container come Docker Hub, ma molti non lo fanno.
In questa sezione opzionale, vi mostreremo due modi per ottenere un'immagine container per gli strumenti che volete utilizzare nelle vostre pipeline Nextflow: utilizzando Seqera Containers e costruendo voi stessi l'immagine container.

Otterrete/costruirete un'immagine container per il pacchetto pip `quote`, che verrà utilizzato nell'esercizio alla fine di questa sezione.

### 1.1. Ottenere un'immagine container da Seqera Containers

Seqera Containers è un servizio gratuito che costruisce immagini container per strumenti installabili tramite pip e conda (incluso bioconda).
Navigate su [Seqera Containers](https://www.seqera.io/containers/) e cercate il pacchetto pip `quote`.

![Seqera Containers](img/seqera-containers-1.png)

Cliccate su "+Add" e poi su "Get Container" per richiedere un'immagine container per il pacchetto pip `quote`.

![Seqera Containers](img/seqera-containers-2.png)

Se questa è la prima volta che viene costruito un container della community per questa versione del pacchetto, potrebbero volerci alcuni minuti per completare l'operazione.
Cliccate per copiare l'URI (ad es. `community.wave.seqera.io/library/pip_quote:ae07804021465ee9`) dell'immagine container che è stata creata per voi.

Ora potete utilizzare l'immagine container per eseguire il comando `quote` e ottenere una citazione casuale di Grace Hopper.

```bash
docker run --rm community.wave.seqera.io/library/pip_quote:ae07804021465ee9 quote "Grace Hopper"
```

Output:

```console title="Output"
Humans are allergic to change. They love to say, 'We've always done it
this way.' I try to fight that. That's why I have a clock on my wall
that runs counter-clockwise.
```

### 1.2. Costruire l'immagine container da soli

Utilizziamo alcuni dettagli di costruzione dal sito web di Seqera Containers per costruire noi stessi l'immagine container per il pacchetto pip `quote`.
Tornate al sito web di Seqera Containers e cliccate sul pulsante "Build Details".

Il primo elemento che esamineremo è il `Dockerfile`, un tipo di file di script che contiene tutti i comandi necessari per costruire l'immagine container.
Abbiamo aggiunto alcuni commenti esplicativi al Dockerfile qui sotto per aiutarvi a capire cosa fa ogni parte.

```Dockerfile title="Dockerfile"
# Parte dall'immagine docker base micromamba
FROM mambaorg/micromamba:1.5.10-noble
# Copia il file conda.yml nel container
COPY --chown=$MAMBA_USER:$MAMBA_USER conda.yml /tmp/conda.yml
# Installa varie utility per l'uso di Nextflow e i pacchetti nel file conda.yml
RUN micromamba install -y -n base -f /tmp/conda.yml \
    && micromamba install -y -n base conda-forge::procps-ng \
    && micromamba env export --name base --explicit > environment.lock \
    && echo ">> CONDA_LOCK_START" \
    && cat environment.lock \
    && echo "<< CONDA_LOCK_END" \
    && micromamba clean -a -y
# Esegue il container come utente root
USER root
# Imposta la variabile d'ambiente PATH per includere la directory di installazione di micromamba
ENV PATH="$MAMBA_ROOT_PREFIX/bin:$PATH"
```

Il secondo elemento che esamineremo è il file `conda.yml`, che contiene l'elenco dei pacchetti che devono essere installati nell'immagine container.

```conda.yml title="conda.yml"
channels:
- conda-forge
- bioconda
dependencies:
- pip
- pip:
  - quote==3.0.0 #
```

Copiate il contenuto di questi file negli stub situati nella directory `containers/build`, quindi eseguite il seguente comando per costruire voi stessi l'immagine container.

!!! Note "Nota"

    Utilizziamo il flag `-t quote:latest` per etichettare l'immagine container con il nome `quote` e il tag `latest`.
    Potremo utilizzare questo tag per fare riferimento all'immagine container quando la eseguiremo su questo sistema.

```bash
docker build -t quote:latest containers/build
```

Dopo che la costruzione è terminata, potete eseguire l'immagine container che avete appena costruito.

```bash
docker run --rm quote:latest quote "Margaret Oakley Dayhoff"
```

### Takeaway

Avete imparato due modi diversi per ottenere un'immagine container per uno strumento che volete utilizzare nelle vostre pipeline Nextflow: utilizzando Seqera Containers e costruendo voi stessi l'immagine container.

### Cosa c'è dopo?

Avete tutto ciò che vi serve per continuare al [prossimo capitolo](./04_hello_genomics.md) di questa serie di formazione.
Potete anche continuare con un esercizio opzionale per recuperare citazioni di pionieri dell'informatica/biologia utilizzando il container `quote` e visualizzarle utilizzando il container `cowsay`.

---

## 2. Fare citare alla mucca scienziati famosi

Questa sezione contiene alcuni esercizi avanzati, per mettere in pratica ciò che avete imparato finora.
Svolgere questi esercizi _non è necessario_ per comprendere le parti successive della formazione, ma fornisce un modo divertente per consolidare i vostri apprendimenti capendo come fare citare alla mucca scienziati famosi.

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
Ad alto livello, per completare questo esercizio dovrete:

- Modificare il `params.input_file` predefinito per puntare al file `pioneers.csv`.
- Creare un processo `getQuote` che utilizza il container `quote` per recuperare una citazione per ogni input.
- Collegare l'output del processo `getQuote` al processo `cowsay` per visualizzare la citazione.

Per l'immagine container `quote`, potete utilizzare quella che avete costruito voi stessi nell'esercizio avanzato precedente o utilizzare quella che avete ottenuto da Seqera Containers.

!!! Hint "Suggerimento"

    Una buona scelta per il blocco `script` del vostro processo getQuote potrebbe essere:
        ```groovy
        script:
            def safe_author = author.tokenize(' ').join('-')
            """
            quote "$author" > quote-${safe_author}.txt
            echo "-${author}" >> quote-${safe_author}.txt
            """
        ```

Potete trovare una soluzione a questo esercizio in `containers/solutions/hello-containers-4.1.nf`.

### 2.2. Modificare la vostra pipeline Nextflow per consentirle di eseguire nelle modalità `quote` e `sayHello`.

Aggiungete della logica di ramificazione alla vostra pipeline per consentirle di accettare input destinati sia a `quote` che a `sayHello`.
Ecco un esempio di come utilizzare un'istruzione `if` in un flusso di lavoro Nextflow:

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

    Potete utilizzare `new_ch = processName.out` per assegnare un nome al canale di output di un processo.

Potete trovare una soluzione a questo esercizio in `containers/solutions/hello-containers-4.2.nf`.

### Takeaway

Sapete come utilizzare i container in Nextflow per eseguire processi e come costruire della logica di ramificazione nelle vostre pipeline!

### Cosa c'è dopo?

Festeggiate, prendetevi una pausa per fare stretching e bevete dell'acqua!

Quando siete pronti, passate alla Parte 3 di questa serie di formazione per imparare come applicare ciò che avete imparato finora a un caso d'uso di analisi dati più realistico.
