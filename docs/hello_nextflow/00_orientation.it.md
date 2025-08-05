# Orientation

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/G3CV-FcV-rc?si=nyLvwhrSB2m1NPc5&amp;list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } Guarda [tutta la playlist](https://www.youtube.com/playlist?list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik) sul canale Youtube Nextflow.

:green_book: La trascrizione di questo video è disponibile [qui](./transcripts/00_orientation.md).
///

## GitHub Codespaces

L'ambiente GitHub Codespaces contiene tutto il software, il codice e i dati necessari per lavorare attraverso questo corso di formazione, quindi non devi installare nulla da solo.
Tuttavia, hai bisogno di un account GitHub (gratuito) per effettuare l'accesso e dovresti prenderti qualche minuto per familiarizzare con l'interfaccia.

Se non l'hai ancora fatto, segui il mini-corso [Environment Setup](../../envsetup/) prima di proseguire.

## Working directory

Durante questo corso di formazione, lavoreremo nella directory `hello-nextflow/`.

Cambia directory ora eseguendo questo comando nel terminale:

```bash
cd hello-nextflow/
```

!!!tip

    Se per qualsiasi motivo esci da questa directory, puoi sempre utilizzare il percorso completo per tornarci, supponendo che tu stia eseguendo questo comando nell'ambiente di formazione Github Codespaces:

    ```bash
    cd /workspaces/training/hello-nextflow
    ```

Diamo ora un'occhiata al contenuto di questa directory.

## Materiali forniti

Puoi esplorare il contenuto di questa directory utilizzando l'esploratore di file sul lato sinistro dell'area di lavoro di training.
In alternativa, puoi utilizzare il comando `tree`.

Durante il corso, utilizziamo l'output `tree` per rappresentare la struttura e il contenuto delle directory in un formato leggibile, a volte con piccole modifiche per chiarezza.

Qui generiamo un indice fino al secondo livello in basso:

```bash
tree . -L 2
```

Se esegui questo comando all'interno di `hello-nextflow`, dovresti vedere il seguente output:

```console title="Directory contents"
.
├── greetings.csv
├── hello-channels.nf
├── hello-config.nf
├── hello-containers.nf
├── hello-modules.nf
├── hello-workflow.nf
├── hello-world.nf
├── nextflow.config
├── solutions
│   ├── 1-hello-world
│   ├── 2-hello-channels
│   ├── 3-hello-workflow
│   ├── 4-hello-modules
│   ├── 5-hello-containers
│   └── 6-hello-config
└── test-params.json

7 directories, 9 files
```

**Ecco un riepilogo di ciò che dovresti sapere per iniziare:**

- **I file `.nf`** sono script di workflow denominati in base alla parte del corso in cui vengono utilizzati.

- **Il file `nextflow.config`** è un file di configurazione che imposta le proprietà minime dell'ambiente.
  Per ora puoi ignorarlo.

- **Il file `greetings.csv`** contiene i dati di input che utilizzeremo nella maggior parte del corso. È descritto nella Parte 1, quando lo introduciamo per la prima volta.

- **Il file `test-params.json`** è un file che utilizzeremo nella Parte 6. Per ora puoi ignorarlo.

- **La directory `solutions`** contiene gli script del workflow completati che risultano da ogni passaggio del corso.
  Sono pensati per essere utilizzati come riferimento per controllare il tuo lavoro e risolvere eventuali problemi.
  Il nome e il numero nel nome del file corrispondono al passaggio della parte pertinente del corso.
  Ad esempio, il file `hello-world-4.nf` è il risultato previsto del completamento dei passaggi da 1 a 4 della Parte 1: Hello World.

**Ora, per iniziare il corso, clicca sulla freccia nell'angolo in basso a destra di questa pagina.**
