---
title: Versioni di Nextflow
description: Comprendere e gestire l'evoluzione delle versioni della sintassi di Nextflow
hide:
  - toc
  - footer
---

## Versione attuale della sintassi Nextflow supportata e requisiti

A partire dalla versione 3.0 del portale di formazione, tutti i nostri corsi di formazione sono basati sulla versione 25.10.2 di Nextflow, salvo diversa indicazione nella pagina indice del corso (ad eccezione dei materiali deprecati o archiviati che potrebbero non includere un avviso di versione).

Poiché i corsi ora utilizzano input tipizzati a livello di workflow così come direttive di output a livello di workflow, richiedono l'uso del parser di sintassi V2.
Se prevede di utilizzare l'ambiente che forniamo attraverso [Github Codespaces](../envsetup/01_setup.md) o [devcontainer locali](../envsetup/03_devcontainer.md), non è necessario fare nulla a meno che non sia specificatamente indicato nelle istruzioni del corso.
Tuttavia, se prevede di seguire le formazioni nel proprio ambiente ([Installazione manuale](../envsetup/02_local.md)), dovrà assicurarsi di utilizzare Nextflow versione 25.10.2 o successiva con il parser di sintassi v2 abilitato.

## Versioni precedenti dei materiali di formazione

I nostri materiali di formazione sono versionati da febbraio 2025.

Può accedere alle versioni precedenti dei materiali di formazione che funzionano con versioni di Nextflow **precedenti alla 25.10.2** tramite il menu a discesa in cima a ogni pagina che mostra la versione numerata dei materiali di formazione.
Quando seleziona una versione precedente dei materiali di formazione, i link all'ambiente di formazione specificheranno automaticamente la versione corrispondente dell'ambiente.

## Altre informazioni sulle versioni della sintassi Nextflow

Nextflow ha due concetti di versionamento distinti che vengono talvolta confusi: **versioni DSL** e **versioni del parser di sintassi**.

**DSL1 vs DSL2** si riferisce a modi fondamentalmente diversi di scrivere pipeline Nextflow.
DSL1 era la sintassi originale dove i process erano implicitamente connessi attraverso i channel.
DSL2, introdotto in Nextflow 20.07, ha aggiunto funzionalità di modularità: la capacità di importare process e workflow da altri file, blocchi `workflow` espliciti e output dei process nominati.
DSL1 è stato deprecato in Nextflow 22.03 e rimosso in 22.12.
Tutto il codice Nextflow moderno utilizza DSL2.

**Parser di sintassi v1 vs v2** si riferisce a parser diversi che funzionano entrambi con codice DSL2.
Il parser v1 è quello originale, più permissivo.
Il parser v2 è più rigoroso e abilita nuove funzionalità del linguaggio come la tipizzazione statica (input e output tipizzati) e direttive di output a livello di workflow.
Il parser v2 fornisce anche messaggi di errore migliori e rileva più errori in fase di parsing piuttosto che a runtime.
Il parser v2 diventerà predefinito in Nextflow 26.04.

In sintesi: DSL2 è il linguaggio che si scrive; la versione del parser di sintassi determina con quanta rigidità viene interpretato quel linguaggio e quali funzionalità avanzate sono disponibili.

### Controllare e impostare la versione di Nextflow

Può verificare quale versione di Nextflow è installata sul sistema utilizzando il comando `nextflow --version`.

Per maggiori informazioni su come aggiornare la versione di Nextflow, consulti la documentazione di riferimento su [Aggiornamento di Nextflow](https://www.nextflow.io/docs/latest/updating-nextflow.html).

### Abilitare il parser di sintassi v2

Per **abilitare** il parser di sintassi v2 per la sessione corrente, eseguite il seguente comando nel terminale:

```bash
export NXF_SYNTAX_PARSER=v2
```

Per rendere questa impostazione permanente (in attesa che v2 diventi predefinito in Nextflow 26.04), aggiunga il comando export al profilo della shell (`~/.bashrc`, `~/.zshrc`, ecc.):

```bash
echo 'export NXF_SYNTAX_PARSER=v2' >> ~/.bashrc
source ~/.bashrc
```

Si noti che la variabile d'ambiente `NXF_SYNTAX_PARSER=v2` è un requisito temporaneo.
Da Nextflow 26.04 in poi, il parser v2 diventerà predefinito e questa impostazione non sarà più necessaria.

### Disabilitare il parser di sintassi v2

Per **disabilitare** il parser di sintassi v2 per la sessione corrente, eseguite il seguente comando nel terminale:

```bash
export NXF_SYNTAX_PARSER=v1
```

<!-- Will it be possible to disable it in versions after 26.04? -->

### Migrazione del codice esistente

Per una guida sulla migrazione del codice esistente per conformarsi alle versioni più recenti di Nextflow, consulti le [Note di migrazione](https://www.nextflow.io/docs/latest/migrations/index.html) nella documentazione di riferimento.

Questi due articoli sono particolarmente utili per la migrazione alla release più recente:

- [Migrazione agli output del workflow](https://www.nextflow.io/docs/latest/tutorials/workflow-outputs.html)
- [Migrazione ai tipi statici](https://www.nextflow.io/docs/latest/tutorials/static-types.html)

Entrambe queste funzionalità sono trattate come parte della formazione per principianti a partire dalla versione 3.0 dei materiali di formazione.

A seconda della generazione di codice Nextflow che intende migrare, potrebbe essere in grado di completare la maggior parte del lavoro tramite il linter di Nextflow utilizzando il comando `nextflow lint -format`.
Consulti il riferimento CLI per [`lint`](https://www.nextflow.io/docs/latest/reference/cli.html#lint) per maggiori dettagli.

Speriamo che questo sia utile.
Se avete bisogno di aiuto, ci contatti su Slack o sul forum.

---

<div markdown class="homepage_logos">

![Seqera](../assets/img/seqera_logo.png#only-light)

![Seqera](../assets/img/seqera_logo_dark.png#only-dark)

</div>
