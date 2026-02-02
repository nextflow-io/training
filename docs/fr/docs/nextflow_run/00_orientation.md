# Démarrage

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduction assistée par IA - [en savoir plus et suggérer des améliorations](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

## Lancer un environnement de formation

Pour utiliser l'environnement pré-construit que nous fournissons sur GitHub Codespaces, cliquez sur le bouton « Open in GitHub Codespaces » ci-dessous. Pour d'autres options, consultez [Options d'environnement](../envsetup/index.md).

Nous vous recommandons d'ouvrir l'environnement de formation dans un nouvel onglet ou une nouvelle fenêtre de navigateur (utilisez le clic droit, ctrl-clic ou cmd-clic selon votre équipement) afin de pouvoir continuer à lire pendant que l'environnement se charge.
Vous devrez garder ces instructions ouvertes en parallèle pour suivre la formation.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

### Les bases de l'environnement

Cet environnement de formation contient tous les logiciels, le code et les données nécessaires pour suivre la formation, vous n'avez donc pas besoin d'installer quoi que ce soit vous-même.

Le codespace est configuré avec une interface VSCode, qui comprend un explorateur de fichiers, un éditeur de code et un terminal shell.
Toutes les instructions données pendant la formation (par exemple « ouvrir le fichier », « modifier le code » ou « exécuter cette commande ») se réfèrent à ces trois parties de l'interface VSCode, sauf indication contraire.

Si vous suivez cette formation par vous-même, veuillez vous familiariser avec les [bases de l'environnement](../envsetup/01_setup.md) pour plus de détails.

### Exigences de version

Cette formation est conçue pour Nextflow 25.10.2 ou ultérieur **avec le parseur de syntaxe v2 ACTIVÉ**.
Si vous utilisez un environnement local ou personnalisé, veuillez vous assurer que vous utilisez les paramètres corrects comme documenté [ici](../info/nxf_versions.md).

## Se préparer à travailler

Une fois votre codespace en cours d'exécution, il y a deux choses que vous devez faire avant de plonger dans la formation : définir votre répertoire de travail pour cette formation spécifique, et examiner les matériaux fournis.

### Définir le répertoire de travail

Par défaut, le codespace s'ouvre avec le répertoire de travail défini à la racine de toutes les formations, mais pour cette formation, nous allons travailler dans le répertoire `nextflow-run/`.

Changez de répertoire maintenant en exécutant cette commande dans le terminal :

```bash
cd nextflow-run/
```

Vous pouvez configurer VSCode pour se concentrer sur ce répertoire, de sorte que seuls les fichiers pertinents apparaissent dans la barre latérale de l'explorateur de fichiers :

```bash
code .
```

!!! tip "Astuce"

    Si pour une raison quelconque vous sortez de ce répertoire (par exemple, votre codespace s'endort), vous pouvez toujours utiliser le chemin complet pour y revenir, en supposant que vous exécutez cela dans l'environnement de formation GitHub Codespaces :

    ```bash
    cd /workspaces/training/nextflow-run
    ```

Maintenant, examinons le contenu.

### Explorer les matériaux fournis

Vous pouvez explorer le contenu de ce répertoire en utilisant l'explorateur de fichiers sur le côté gauche de l'espace de travail de formation.
Alternativement, vous pouvez utiliser la commande `tree`.

Tout au long de la formation, nous utilisons la sortie de `tree` pour représenter la structure et le contenu des répertoires sous une forme lisible, parfois avec des modifications mineures pour plus de clarté.

Ici, nous générons une table des matières jusqu'au deuxième niveau :

```bash
tree . -L 2
```

??? abstract "Contenu du répertoire"

    ```console
    .
    ├── 1-hello.nf
    ├── 2a-inputs.nf
    ├── 2b-multistep.nf
    ├── 2c-modules.nf
    ├── 2d-container.nf
    ├── 3-main.nf
    ├── data
    │   └── greetings.csv
    ├── modules
    │   ├── collectGreetings.nf
    │   ├── convertToUpper.nf
    │   ├── cowpy.nf
    │   └── sayHello.nf
    ├── nextflow.config
    ├── solutions
    │   ├── 3-main.nf
    │   ├── modules
    │   └── nextflow.config
    ├── test-params.json
    └── test-params.yaml
    ```

Cliquez sur la boîte colorée pour développer la section et voir son contenu.
Nous utilisons des sections repliables comme celle-ci pour afficher la sortie attendue des commandes ainsi que le contenu des répertoires et des fichiers de manière concise.

- **Les fichiers `.nf`** sont des scripts de workflow numérotés en fonction de la partie de la formation où ils sont utilisés.

- **Le fichier `nextflow.config`** est un fichier de configuration qui définit les propriétés minimales de l'environnement.
  Vous pouvez l'ignorer pour l'instant.

- **Le fichier `greetings.csv`** sous `data/` contient les données d'entrée que nous utiliserons dans la majeure partie de la formation. Il est décrit dans la Partie 2 (Exécuter des pipelines), lorsque nous l'introduisons pour la première fois.

- **Les fichiers `test-params.*`** sont des fichiers de configuration que nous utiliserons dans la Partie 3 (Configuration). Vous pouvez les ignorer pour l'instant.

- **Le répertoire `solutions`** contient l'état final du workflow et de ses fichiers accessoires (config et modules) qui résultent de l'achèvement de la formation.
  Ils sont destinés à être utilisés comme référence pour vérifier votre travail et résoudre les problèmes.

## Liste de vérification de préparation

Vous pensez être prêt·e à plonger ?

- [ ] Je comprends l'objectif de cette formation et ses prérequis
- [ ] Mon environnement est opérationnel
- [ ] J'ai défini mon répertoire de travail de manière appropriée

Si vous pouvez cocher toutes les cases, vous êtes prêt·e à commencer.

**Pour continuer vers [Partie 1 : Exécuter les opérations de base](./01_basics.md), cliquez sur la flèche en bas à droite de cette page.**
