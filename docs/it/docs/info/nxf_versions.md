---
title: Versioni di Nextflow
description: Comprendere e gestire l'evoluzione delle versioni della sintassi di Nextflow
hide:
  - toc
  - footer
---

## Versione della sintassi di Nextflow attualmente supportata e requisiti

A partire dalla versione 3.0 del portale di formazione, tutti i nostri corsi di formazione si basano sulla versione 25.10.2 di Nextflow, salvo diversa indicazione nella pagina indice del corso (ad eccezione dei materiali deprecati o archiviati che potrebbero non includere un avviso sulla versione).

Poiché i corsi ora utilizzano input tipizzati a livello di flusso di lavoro e direttive di output a livello di flusso di lavoro, richiedono l'uso del parser di sintassi V2.
Se prevedete di utilizzare l'ambiente che forniamo tramite [Github Codespaces](../envsetup/01_setup.md) o [devcontainer locali](../envsetup/03_devcontainer.md), non dovete fare nulla a meno che non sia specificato nelle istruzioni del corso.
Tuttavia, se prevedete di seguire le formazioni nel vostro ambiente ([Installazione manuale](../envsetup/02_local.md)), dovrete assicurarvi di utilizzare Nextflow versione 25.10.2 o successiva con il parser di sintassi v2 abilitato.

## Versioni precedenti dei materiali di formazione

I nostri materiali di formazione sono versionati da febbraio 2025.

Potete accedere alle versioni precedenti dei materiali di formazione che funzionano con versioni di Nextflow **precedenti alla 25.10.2** tramite il menu a tendina nella parte superiore di ogni pagina che mostra la versione numerata dei materiali di formazione.
Quando selezionate una versione precedente dei materiali di formazione, i link all'ambiente di formazione specificheranno automaticamente la versione corrispondente dell'ambiente.

## Altre informazioni sulle versioni della sintassi di Nextflow

Nextflow ha due concetti di versionamento distinti che a volte vengono confusi: **versioni DSL** e **versioni del parser di sintassi**.

**DSL1 vs DSL2** si riferisce a modi fondamentalmente diversi di scrivere pipeline Nextflow.
DSL1 era la sintassi originale in cui i processi erano collegati implicitamente attraverso i canali.
DSL2, introdotto in Nextflow 20.07, ha aggiunto funzionalità di modularità: la possibilità di importare processi e flussi di lavoro da altri file, blocchi `workflow` espliciti e output di processo nominati.
DSL1 è stato deprecato in Nextflow 22.03 e rimosso nella 22.12.
Tutto il codice Nextflow moderno utilizza DSL2.

**Parser di sintassi v1 vs v2** si riferisce a parser diversi che entrambi funzionano con il codice DSL2.
Il parser v1 è il parser originale, più permissivo.
Il parser v2 è più rigoroso e abilita nuove funzionalità del linguaggio come la tipizzazione statica (input e output tipizzati) e direttive di output a livello di flusso di lavoro.
Il parser v2 fornisce anche messaggi di errore migliori e rileva più errori al momento del parsing piuttosto che durante l'esecuzione.
Il parser v2 diventerà il predefinito in Nextflow 26.04.

In sintesi: DSL2 è il linguaggio che scrivete; la versione del parser di sintassi determina quanto rigorosamente quel linguaggio viene interpretato e quali funzionalità avanzate sono disponibili.

### Verificare e impostare la versione di Nextflow

Potete verificare quale versione di Nextflow è installata sul vostro sistema utilizzando il comando `nextflow --version`.

Per maggiori informazioni su come aggiornare la vostra versione di Nextflow, consultate la documentazione di riferimento su [Aggiornamento di Nextflow](https://www.nextflow.io/docs/latest/updating-nextflow.html).

### Abilitare il parser di sintassi v2

Per **abilitare** il parser di sintassi v2 per la sessione corrente, eseguite il seguente comando nel vostro terminale:

```bash
export NXF_SYNTAX_PARSER=v2
```

Per rendere questa impostazione permanente (in attesa che v2 diventi il predefinito in Nextflow 26.04), aggiungete il comando export al vostro profilo shell (`~/.bashrc`, `~/.zshrc`, ecc.):

```bash
echo 'export NXF_SYNTAX_PARSER=v2' >> ~/.bashrc
source ~/.bashrc
```

Notate che la variabile d'ambiente `NXF_SYNTAX_PARSER=v2` è un requisito temporaneo.
Da Nextflow 26.04 in poi, il parser v2 diventerà il predefinito e questa impostazione non sarà più necessaria.

### Disabilitare il parser di sintassi v2

Per **disabilitare** il parser di sintassi v2 per la sessione corrente, eseguite il seguente comando nel vostro terminale:

```bash
export NXF_SYNTAX_PARSER=v1
```

<!-- Sarà possibile disabilitarlo nelle versioni successive alla 26.04? -->

### Migrazione del codice esistente

Per indicazioni sulla migrazione del codice esistente per conformarsi alle versioni più recenti di Nextflow, consultate le [Note di migrazione](https://www.nextflow.io/docs/latest/migrations/index.html) nella documentazione di riferimento.

Questi due articoli sono particolarmente utili per la migrazione alla versione più recente:

- [Migrazione agli output del flusso di lavoro](https://www.nextflow.io/docs/latest/tutorials/workflow-outputs.html)
- [Migrazione ai tipi statici](https://www.nextflow.io/docs/latest/tutorials/static-types.html)

Entrambe queste funzionalità sono trattate come parte della formazione per principianti a partire dalla versione 3.0 dei materiali di formazione.

A seconda della generazione del codice Nextflow che intendete migrare, potreste essere in grado di completare la maggior parte del lavoro utilizzando il linter di Nextflow con il comando `nextflow lint -format`.
Consultate il riferimento CLI per [`lint`](https://www.nextflow.io/docs/latest/reference/cli.html#lint) per maggiori dettagli.

Speriamo che questo sia utile.
Se avete bisogno di aiuto, contattateci su Slack o sul forum.

---

<div markdown class="homepage_logos">

![Seqera](../assets/img/seqera_logo.png#only-light)

![Seqera](../assets/img/seqera_logo_dark.png#only-dark)

</div>
