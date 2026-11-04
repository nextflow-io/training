# Premiers pas

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduction assistée par IA - [en savoir plus et suggérer des améliorations](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

## Démarrer un environnement de formation

Pour utiliser l'environnement pré-configuré que nous fournissons sur GitHub Codespaces, cliquez sur le bouton « Open in GitHub Codespaces » ci-dessous. Pour d'autres options, consultez [Options d'environnement](../../envsetup/index.md).

Nous vous recommandons d'ouvrir l'environnement de formation dans un nouvel onglet ou une nouvelle fenêtre de navigateur (utilisez le clic droit, ctrl-clic ou cmd-clic selon votre équipement) afin de pouvoir continuer à lire pendant le chargement de l'environnement.
Vous devrez garder ces instructions ouvertes en parallèle pour suivre le cours.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

### Bases de l'environnement

Cet environnement de formation contient tous les logiciels, le code et les données nécessaires pour suivre le cours de formation, vous n'avez donc rien à installer vous-même.

Le codespace est configuré avec une interface VSCode, qui comprend un explorateur de système de fichiers, un éditeur de code et un terminal shell.
Toutes les instructions données pendant le cours (par exemple « ouvrir le fichier », « modifier le code » ou « exécuter cette commande ») font référence à ces trois parties de l'interface VSCode sauf indication contraire.

Si vous suivez ce cours par vous-même, veuillez vous familiariser avec les [bases de l'environnement](../../envsetup/01_setup.md) pour plus de détails.

### Exigences de version

Cette formation est conçue pour Nextflow 25.10.2 ou ultérieur **avec l'analyseur de syntaxe v2 ACTIVÉ**.
Si vous utilisez un environnement local ou personnalisé, assurez-vous d'utiliser les paramètres corrects comme documenté [ici](../../info/nxf_versions.md).

## Préparez-vous à travailler

Une fois votre codespace en cours d'exécution, vous devez faire deux choses avant de plonger dans la formation : définir votre répertoire de travail pour ce cours spécifique et examiner les ressources fournies.

### Définir le répertoire de travail

Par défaut, le codespace s'ouvre avec le répertoire de travail défini à la racine de tous les cours de formation, mais pour ce cours, nous travaillerons dans le répertoire `nf4-science/{DOMAIN_DIR}/`.

Changez de répertoire maintenant en exécutant cette commande dans le terminal :

```bash
cd nf4-science/{DOMAIN_DIR}/
```

Vous pouvez configurer VSCode pour qu'il se concentre sur ce répertoire, afin que seuls les fichiers pertinents s'affichent dans la barre latérale de l'explorateur de fichiers :

```bash
code .
```

!!! tip "Astuce"

    Si pour une raison quelconque vous sortez de ce répertoire (par exemple, votre codespace se met en veille), vous pouvez toujours utiliser le chemin complet pour y revenir, en supposant que vous exécutez ceci dans l'environnement de formation GitHub Codespaces :

    ```bash
    cd /workspaces/training/nf4-science/{DOMAIN_DIR}
    ```

Maintenant, examinons le contenu.

### Explorer les ressources fournies

Vous pouvez explorer le contenu de ce répertoire en utilisant l'explorateur de fichiers sur le côté gauche de l'espace de travail de formation.
Alternativement, vous pouvez utiliser la commande `tree`.

Tout au long du cours, nous utilisons la sortie de `tree` pour représenter la structure et le contenu des répertoires sous une forme lisible, parfois avec des modifications mineures pour plus de clarté.

Ici, nous générons une table des matières jusqu'au deuxième niveau :

```bash
tree . -L 2
```

??? abstract "Contenu du répertoire"

    ```console
    .
    ├── data
    │   ├── {DATA_SUBDIRS}
    │   └── samplesheet.csv
    ├── {DOMAIN_DIR}.nf
    ├── modules
    │   ├── {TOOL_A_MODULE}.nf
    │   └── {TOOL_B_MODULE}.nf
    ├── nextflow.config
    └── solutions
        ├── modules
        ├── part2
        └── part3

    N directories, N files
    ```

Cliquez sur la boîte colorée pour développer la section et afficher son contenu.
Nous utilisons des sections repliables comme celle-ci pour afficher la sortie de commande attendue ainsi que le contenu des répertoires et des fichiers de manière concise.

- **Le fichier `{DOMAIN_DIR}.nf`** est un script de workflow que vous construirez au fil du cours.

- **Le répertoire `modules`** contient des fichiers de modules squelettes que vous remplirez pendant le cours.

- **Le fichier `nextflow.config`** est un fichier de configuration qui définit les propriétés minimales de l'environnement.
  Vous pouvez l'ignorer pour l'instant.

- **Le répertoire `data`** contient les données d'entrée et les ressources associées, décrites plus tard dans le cours.

- **Le répertoire `solutions`** contient des fichiers de modules complétés et des solutions spécifiques à chaque partie qui peuvent servir de point de départ pour la partie suivante.
  Ils sont destinés à être utilisés comme référence pour vérifier votre travail et résoudre tout problème.

## Liste de vérification de préparation

Pensez-vous être prêt·e à vous lancer ?

- [ ] Je comprends l'objectif de ce cours et ses prérequis
- [ ] Mon environnement est opérationnel
- [ ] J'ai défini mon répertoire de travail de manière appropriée

Si vous pouvez cocher toutes les cases, vous êtes prêt·e à commencer.

**Pour continuer vers la [Partie 1 : Aperçu de la méthode](./01_method.md), cliquez sur la flèche dans le coin inférieur droit de cette page.**
