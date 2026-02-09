# Premiers pas

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduction assistée par IA - [en savoir plus et suggérer des améliorations](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

## Démarrer un environnement de formation

Pour utiliser l'environnement pré-construit que nous fournissons sur GitHub Codespaces, cliquez sur le bouton « Open in GitHub Codespaces » ci-dessous. Pour d'autres options, consultez [Options d'environnement](../../envsetup/index.md).

Nous recommandons d'ouvrir l'environnement de formation dans un nouvel onglet ou une nouvelle fenêtre de navigateur (utilisez clic droit, ctrl-clic ou cmd-clic selon votre équipement) afin de pouvoir continuer à lire pendant le chargement de l'environnement.
Vous devrez garder ces instructions ouvertes en parallèle pour suivre le cours.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

### Bases de l'environnement

Cet environnement de formation contient tous les logiciels, le code et les données nécessaires pour suivre cette formation, vous n'avez donc rien à installer vous-même.

Le codespace est configuré avec une interface VSCode, qui comprend un explorateur de système de fichiers, un éditeur de code et un terminal shell.
Toutes les instructions données pendant le cours (par exemple « ouvrez le fichier », « modifiez le code » ou « exécutez cette commande ») font référence à ces trois parties de l'interface VScode sauf indication contraire.

Si vous suivez ce cours par vous-même, veuillez vous familiariser avec les [bases de l'environnement](../../envsetup/01_setup.md) pour plus de détails.

### Exigences de version

Cette formation est conçue pour Nextflow 25.10.2 ou ultérieur **avec le parseur de syntaxe v2 ACTIVÉ**.
Si vous utilisez un environnement local ou personnalisé, assurez-vous d'utiliser les paramètres corrects comme documenté [ici](../../info/nxf_versions.md).

## Préparez-vous à travailler

Une fois votre codespace en cours d'exécution, vous devez faire deux choses avant de plonger dans la formation : définir votre répertoire de travail pour ce cours spécifique et jeter un œil au matériel fourni.

### Définir le répertoire de travail

Par défaut, le codespace s'ouvre avec le répertoire de travail défini à la racine de tous les cours de formation, mais pour ce cours, nous travaillerons dans le répertoire `nf4-science/genomics/`.

Changez de répertoire maintenant en exécutant cette commande dans le terminal :

```bash
cd nf4-science/genomics/
```

Vous pouvez configurer VSCode pour qu'il se concentre sur ce répertoire, de sorte que seuls les fichiers pertinents s'affichent dans la barre latérale de l'explorateur de fichiers :

```bash
code .
```

!!! tip

    Si pour une raison quelconque vous sortez de ce répertoire (par exemple si votre codespace se met en veille), vous pouvez toujours utiliser le chemin complet pour y revenir, en supposant que vous travaillez dans l'environnement de formation Github Codespaces :

    ```bash
    cd /workspaces/training/nf4-science/genomics
    ```

Maintenant, jetons un œil au contenu.

### Explorer le matériel fourni

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
    │   ├── bam
    │   ├── ref
    │   ├── sample_bams.txt
    │   └── samplesheet.csv
    ├── genomics.nf
    ├── modules
    │   ├── gatk_haplotypecaller.nf
    │   └── samtools_index.nf
    ├── nextflow.config
    └── solutions
        ├── modules
        ├── nf-test.config
        ├── part2
        └── tests

    8 directories, 8 files
    ```

Cliquez sur la boîte colorée pour développer la section et voir son contenu.
Nous utilisons des sections repliables comme celle-ci pour afficher la sortie de commande attendue ainsi que le contenu des répertoires et des fichiers de manière concise.

- **Le fichier `genomics.nf`** est un script de workflow que vous construirez au fil du cours.

- **Le répertoire `modules`** contient des fichiers de modules squelettes que vous remplirez pendant le cours.

- **Le fichier `nextflow.config`** est un fichier de configuration qui définit des propriétés minimales de l'environnement.
  Vous pouvez l'ignorer pour l'instant.

- **Le répertoire `data`** contient les données d'entrée et les ressources associées, décrites plus tard dans le cours.

- **Le répertoire `solutions`** contient des fichiers de modules complétés et une solution pour la Partie 2 qui peut servir de point de départ pour la Partie 3.
  Ils sont destinés à être utilisés comme référence pour vérifier votre travail et résoudre tout problème.

## Liste de vérification de préparation

Pensez-vous être prêt·e à vous lancer ?

- [ ] Je comprends l'objectif de ce cours et ses prérequis
- [ ] Mon environnement est opérationnel
- [ ] J'ai défini mon répertoire de travail de manière appropriée

Si vous pouvez cocher toutes les cases, vous êtes prêt·e à commencer.

**Pour continuer vers [Partie 1 : Aperçu de la méthode et tests manuels](./01_method.md), cliquez sur la flèche dans le coin inférieur droit de cette page.**
