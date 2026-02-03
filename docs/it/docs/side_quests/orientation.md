# Orientamento

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduzione assistita da IA - [scopri di più e suggerisci miglioramenti](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

L'ambiente GitHub Codespaces contiene tutto il software, il codice e i dati necessari per seguire questo corso di formazione, quindi non è necessario installare nulla autonomamente.
Tuttavia, è necessario un account (gratuito) per effettuare l'accesso, e dovrebbe dedicare alcuni minuti per familiarizzare con l'interfaccia.

Se non lo ha ancora fatto, segua [questo link](../../envsetup/) prima di proseguire.

## Materiali forniti

Durante questo corso di formazione, lavoreremo nella directory `side-quests/`.
Questa directory contiene tutti i file di codice, i dati di test e i file accessori di cui avrà bisogno.

Si senta libero di esplorare il contenuto di questa directory; il modo più semplice per farlo è utilizzare l'explorer dei file sul lato sinistro dell'area di lavoro di GitHub Codespaces.
In alternativa, potete utilizzare il comando `tree`.
Durante il corso, utilizziamo l'output di `tree` per rappresentare la struttura e il contenuto delle directory in forma leggibile, talvolta con modifiche minori per maggiore chiarezza.

Qui generiamo un indice dei contenuti fino al secondo livello:

```bash
tree . -L 2
```

Se esegue questo comando all'interno di `side-quests`, dovrebbe vedere il seguente output:

```console title="Contenuto della directory"
.
├── metadata
├── nf-core
├── nf-test
├── solutions
├── splitting_and_grouping
└── workflows_of_workflows
```

**Ecco un riepilogo di ciò che dovrebbe sapere per iniziare:**

- **Ogni directory corrisponde a una side quest individuale.**
  I loro contenuti sono dettagliati nella pagina della corrispondente side quest.

- **La directory `solutions`** contiene gli script di workflow e/o modulo completati che risultano dall'esecuzione dei vari passaggi di ciascuna side quest.
  Sono destinati ad essere utilizzati come riferimento per verificare il suo lavoro e risolvere eventuali problemi.

!!!tip "Suggerimento"

    Se per qualsiasi motivo doveste uscire da questa directory, potete sempre eseguire questo comando per tornarvi:

    ```bash
    cd /workspaces/training/side-quests
    ```

Ora, per iniziare il corso, cliccate sulla freccia nell'angolo in basso a destra di questa pagina.
