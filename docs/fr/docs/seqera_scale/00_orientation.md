# Premiers pas

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduction assistée par IA - [en savoir plus et suggérer des améliorations](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

## Démarrer un environnement de formation

Pour utiliser l'environnement préconstruit que nous fournissons sur GitHub Codespaces, cliquez sur le bouton "Open in GitHub Codespaces" ci-dessous. Pour d'autres options, consultez [Options d'environnement](../envsetup/index.md).

Nous vous recommandons d'ouvrir l'environnement de formation dans un nouvel onglet ou une nouvelle fenêtre de navigateur (utilisez le clic droit, ctrl+clic ou cmd+clic selon votre équipement) afin de pouvoir continuer à lire pendant que l'environnement se charge.
Vous devrez garder ces instructions ouvertes en parallèle pour suivre le cours.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

### Bases de l'environnement

Cet environnement de formation contient tous les logiciels, le code et les données nécessaires pour suivre le cours, vous n'avez donc rien à installer vous-même.

Le codespace est configuré avec une interface VSCode, qui comprend un explorateur de fichiers, un éditeur de code et un terminal.
Toutes les instructions données pendant le cours (par exemple, « ouvrir le fichier », « modifier le code » ou « exécuter cette commande ») font référence à ces trois parties de l'interface VSCode, sauf indication contraire.

Si vous suivez ce cours de manière autonome, veuillez vous familiariser avec les [bases de l'environnement](../envsetup/01_setup.md) pour plus de détails.

## Se préparer à travailler

Une fois votre codespace en cours d'exécution, il y a deux choses à faire avant de commencer : définir votre répertoire de travail et examiner les ressources fournies.

### Définir le répertoire de travail

Par défaut, le codespace s'ouvre à la racine de tous les cours de formation.
Pour ce cours, accédez au répertoire `seqera-scale/` :

```bash
cd seqera-scale/
```

Configurez ensuite VSCode pour qu'il se concentre sur ce répertoire, afin que seuls les fichiers pertinents apparaissent dans la barre latérale de l'explorateur de fichiers :

```bash
code .
```

!!! tip "Astuce"

    Si pour une raison quelconque vous quittez ce répertoire (par exemple, si votre codespace se met en veille), vous pouvez toujours utiliser le chemin complet pour y revenir, en supposant que vous exécutez ceci dans l'environnement de formation GitHub Codespaces :

    ```bash
    cd /workspaces/training/seqera-scale
    ```

### Explorer les ressources fournies

Vous pouvez explorer les ressources du cours à l'aide de l'explorateur de fichiers sur la gauche, ou avec la commande `tree`.
Exécutez la commande suivante depuis le terminal pour voir la structure complète :

```bash
tree -a .
```

??? abstract "Contenu du répertoire"

    ```console
    .
    └── .seqera_config
    ```

Le fichier **`.seqera_config`** est un modèle que vous remplirez lors de la section 3 pour configurer le CLI `tw` avec votre token d'accès Seqera et votre workspace.

## Liste de vérification

Vous pensez être prêt·e à vous lancer ?

- [ ] Je comprends l'objectif de ce cours et ses prérequis
- [ ] Mon environnement est opérationnel
- [ ] J'ai correctement défini mon répertoire de travail

Si vous pouvez cocher toutes les cases, vous êtes prêt·e à commencer.

**Pour continuer vers [Partie 1 : Lancer des pipelines depuis l'interface web](./01_run_with_seqera.md), cliquez sur la flèche dans le coin inférieur droit de cette page.**
