# Premiers pas

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduction assistée par IA - [en savoir plus et suggérer des améliorations](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>


## Démarrer un environnement de formation

Pour utiliser l'environnement préconstruit que nous fournissons sur GitHub Codespaces, cliquez sur le bouton "Open in GitHub Codespaces" ci-dessous. Pour d'autres options, consultez [Options d'environnement](../envsetup/index.md).

Nous vous recommandons d'ouvrir l'environnement de formation dans un nouvel onglet ou une nouvelle fenêtre de navigateur (utilisez le clic droit, ctrl+clic ou cmd+clic selon votre équipement) afin de pouvoir continuer à lire pendant que l'environnement se charge.
Vous devrez garder ces instructions ouvertes en parallèle pour suivre le cours.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

### Bases de l'environnement

Cet environnement de formation contient tous les logiciels, le code et les données nécessaires pour suivre le cours de formation, vous n'avez donc rien à installer vous-même.

Le codespace est configuré avec une interface VSCode, qui comprend un explorateur de fichiers, un éditeur de code et un terminal.
Toutes les instructions données pendant le cours (par exemple, "ouvrir le fichier", "modifier le code" ou "exécuter cette commande") font référence à ces trois parties de l'interface VSCode, sauf indication contraire.

Si vous suivez ce cours par vous-même, veuillez vous familiariser avec les [bases de l'environnement](../envsetup/01_setup.md) pour plus de détails.

### Prérequis de version

Ce cours nécessite Nextflow 25.10.2 ou une version ultérieure, avec le parseur de syntaxe v2 activé (valeur par défaut dans la version 25.10+).
Si vous utilisez un environnement local ou personnalisé, veuillez vous assurer que vous utilisez les paramètres corrects tels que documentés [ici](../info/nxf_versions.md).

## Se préparer à travailler

Une fois votre codespace en cours d'exécution, il y a deux choses à faire avant de commencer : définir votre répertoire de travail et jeter un coup d'œil aux matériaux fournis.

### Définir le répertoire de travail

Par défaut, le codespace s'ouvre à la racine de tous les cours de formation.
Pour ce cours, changez vers le répertoire `config-exec/` :

```bash
cd config-exec/
```

Ensuite, configurez VSCode pour qu'il se concentre sur ce répertoire, afin que seuls les fichiers pertinents apparaissent dans la barre latérale de l'explorateur de fichiers :

```bash
code .
```

!!! tip "Astuce"

    Si pour une raison quelconque vous quittez ce répertoire (par exemple, si votre codespace se met en veille), vous pouvez toujours utiliser le chemin complet pour y revenir, en supposant que vous exécutez ceci dans l'environnement de formation GitHub Codespaces :

    ```bash
    cd /workspaces/training/config-exec
    ```

### Explorer les matériaux fournis

Vous pouvez explorer les matériaux du cours à l'aide de l'explorateur de fichiers sur la gauche, ou avec la commande `tree`.
Exécutez la commande suivante depuis le terminal pour voir la structure complète :

```bash
tree . -L 2
```

??? abstract "Contenu du répertoire"

    ```console
    .
    ├── data
    │   └── greetings.csv
    ├── main.nf
    ├── modules
    │   ├── collectGreetings.nf
    │   ├── convertToUpper.nf
    │   ├── cowpy.nf
    │   └── sayHello.nf
    └── nextflow.config
    ```

Les fichiers **`main.nf`** et **`modules/`** correspondent au même pipeline multi-étapes issu de [Nextflow Run](../nextflow_run/index.md), et le fichier **`nextflow.config`** est la même configuration que vous avez déjà vue là-bas.
Vous les étendrez tous les deux au fil des exercices.

Le répertoire **`data/`** contient le fichier d'entrée CSV que le pipeline lit.

## Liste de vérification

Vous pensez être prêt·e à vous lancer ?

- [ ] Je comprends l'objectif de ce cours et ses prérequis
- [ ] Mon environnement est opérationnel
- [ ] J'ai défini mon répertoire de travail de manière appropriée

Si vous pouvez cocher toutes les cases, vous êtes prêt·e à commencer.

**Pour continuer vers [Partie 1 : Adapter à votre environnement de calcul](./01_packaging_and_execution.md), cliquez sur la flèche dans le coin inférieur droit de cette page.**
