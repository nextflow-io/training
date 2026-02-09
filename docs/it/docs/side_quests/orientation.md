# Orientamento

L'ambiente GitHub Codespaces contiene tutto il software, il codice e i dati necessari per seguire questo corso di formazione, quindi non è necessario installare nulla.
Tuttavia, è necessario un account (gratuito) per effettuare il login, e dovreste prendervi qualche minuto per familiarizzare con l'interfaccia.

Se non l'avete ancora fatto, seguite [questo link](../../envsetup/) prima di procedere.

## Materiali forniti

Durante questo corso di formazione, lavoreremo nella directory `side-quests/`.
Questa directory contiene tutti i file di codice, i dati di test e i file accessori di cui avrete bisogno.

Sentitevi liberi di esplorare i contenuti di questa directory; il modo più semplice per farlo è utilizzare l'esploratore di file sul lato sinistro dell'area di lavoro di GitHub Codespaces.
In alternativa, potete utilizzare il comando `tree`.
Durante il corso, utilizziamo l'output di `tree` per rappresentare la struttura e i contenuti della directory in forma leggibile, a volte con piccole modifiche per maggiore chiarezza.

Qui generiamo un indice dei contenuti fino al secondo livello:

```bash
tree . -L 2
```

Se eseguite questo comando all'interno di `side-quests`, dovreste vedere il seguente output:

```console title="Directory contents"
.
├── metadata
├── nf-core
├── nf-test
├── solutions
├── splitting_and_grouping
└── workflows_of_workflows
```

**Ecco un riepilogo di ciò che dovreste sapere per iniziare:**

- **Ogni directory corrisponde a una singola side quest.**
  I loro contenuti sono dettagliati nella pagina della corrispondente side quest.

- **La directory `solutions`** contiene gli script di flusso di lavoro e/o modulo completati che risultano dall'esecuzione dei vari passaggi di ciascuna side quest.
  Sono destinati ad essere utilizzati come riferimento per verificare il vostro lavoro e risolvere eventuali problemi.

!!!tip

    Se per qualsiasi motivo vi spostate da questa directory, potete sempre eseguire questo comando per ritornarvi:

    ```bash
    cd /workspaces/training/side-quests
    ```

Ora, per iniziare il corso, cliccate sulla freccia nell'angolo in basso a destra di questa pagina.
