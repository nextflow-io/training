# Orientamento

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduzione assistita da IA - [scopri di più e suggerisci miglioramenti](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

L'ambiente di formazione contiene tutto il software, il codice e i dati necessari per seguire questo corso di formazione, quindi non è necessario installare nulla autonomamente.
Tuttavia, è necessario un account (gratuito) per accedere, e si dovrebbe dedicare qualche minuto per familiarizzare con l'interfaccia.

Se non lo ha ancora fatto, segua [questo link](../../../envsetup/) prima di procedere oltre.

## Materiali forniti

Durante questo corso di formazione, lavoreremo nella directory `nf4-science/genomics/`, nella quale è necessario spostarsi quando si apre lo spazio di lavoro della formazione.
Questa directory contiene tutti i file di codice, i dati di test e i file accessori di cui avrà bisogno.

Si senta libero di esplorare i contenuti di questa directory; il modo più semplice per farlo è utilizzare l'esploratore di file sul lato sinistro dello spazio di lavoro della formazione nell'interfaccia VSCode.
In alternativa, potete utilizzare il comando `tree`.
Durante il corso, utilizziamo l'output di `tree` per rappresentare la struttura e i contenuti della directory in forma leggibile, talvolta con modifiche minori per chiarezza.

Qui generiamo un indice dei contenuti fino al secondo livello:

```bash
tree . -L 2
```

Se esegue questo comando all'interno di `nf4-science/genomics`, dovrebbe vedere il seguente output:

```console title="Contenuti della directory"

.
├── data
│   ├── bam
│   ├── ref
│   ├── sample_bams.txt
│   └── samplesheet.csv
├── genomics-1.nf
├── genomics-2.nf
├── genomics-3.nf
├── genomics-4.nf
├── nextflow.config
└── solutions
    ├── modules
    ├── nf-test.config
    └── tests

6 directories, 8 files

```

!!!note "Nota"

    Non si preoccupi se questo sembra molto; esamineremo le parti rilevanti in ogni fase del corso.
    Questo è solo per fornirLe una panoramica.

**Ecco un riepilogo di ciò che dovrebbe sapere per iniziare:**

- **I file `.nf`** sono script di workflow denominati in base alla parte del corso in cui vengono utilizzati.

- **Il file `nextflow.config`** è un file di configurazione che imposta proprietà minime dell'ambiente.
  Può ignorarlo per ora.

- **La directory `data`** contiene dati di input e risorse correlate, descritte più avanti nel corso.

- **La directory `solutions`** contiene file di moduli e configurazioni di test che risultano dalle Parti 3 e 4 del corso.
  Sono destinati ad essere utilizzati come riferimento per verificare il vostro lavoro e risolvere eventuali problemi.

!!!tip "Suggerimento"

    Se per qualsiasi motivo vi spostate fuori da questa directory, potete sempre eseguire questo comando per ritornarvi:

    ```bash
    cd /workspaces/training/nf4-science/genomics
    ```

Ora, per iniziare il corso, cliccate sulla freccia nell'angolo in basso a destra di questa pagina.
