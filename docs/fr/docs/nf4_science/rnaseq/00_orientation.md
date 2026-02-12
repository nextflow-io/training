# Premiers pas

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduction assistée par IA - [en savoir plus et suggérer des améliorations](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

## Démarrer un environnement de formation

Pour utiliser l'environnement pré-configuré que nous fournissons sur GitHub Codespaces, cliquez sur le bouton « Open in GitHub Codespaces » ci-dessous. Pour d'autres options, consultez [Options d'environnement](../../envsetup/index.md).

Nous vous recommandons d'ouvrir l'environnement de formation dans un nouvel onglet ou une nouvelle fenêtre de navigateur (utilisez le clic droit, ctrl-clic ou cmd-clic selon votre équipement) afin de pouvoir continuer à lire pendant le chargement de l'environnement.
Vous devrez garder ces instructions ouvertes en parallèle pour suivre le cours.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

### Bases de l'environnement

Cet environnement de formation contient tous les logiciels, le code et les données nécessaires pour suivre ce cours de formation, vous n'avez donc rien à installer vous-même.

Le codespace est configuré avec une interface VSCode, qui comprend un explorateur de système de fichiers, un éditeur de code et un terminal shell.
Toutes les instructions données pendant le cours (par exemple « ouvrir le fichier », « modifier le code » ou « exécuter cette commande ») font référence à ces trois parties de l'interface VSCode sauf indication contraire.

Si vous suivez ce cours par vous-même, veuillez vous familiariser avec les [bases de l'environnement](../../envsetup/01_setup.md) pour plus de détails.

### Exigences de version

Cette formation est conçue pour Nextflow 25.10.2 ou ultérieur **avec l'analyseur de syntaxe v2 ACTIVÉ**.
Si vous utilisez un environnement local ou personnalisé, veuillez vous assurer d'utiliser les paramètres corrects comme documenté [ici](../../info/nxf_versions.md).

## Préparez-vous à travailler

Une fois votre codespace en cours d'exécution, vous devez faire deux choses avant de plonger dans la formation : définir votre répertoire de travail pour ce cours spécifique et examiner le matériel fourni.

### Définir le répertoire de travail

Par défaut, le codespace s'ouvre avec le répertoire de travail défini à la racine de tous les cours de formation, mais pour ce cours, nous travaillerons dans le répertoire `nf4-science/rnaseq/`.

Changez de répertoire maintenant en exécutant cette commande dans le terminal :

```bash
cd nf4-science/rnaseq/
```

Vous pouvez configurer VSCode pour se concentrer sur ce répertoire, afin que seuls les fichiers pertinents s'affichent dans la barre latérale de l'explorateur de fichiers :

```bash
code .
```

!!! tip "Astuce"

    Si pour une raison quelconque vous sortez de ce répertoire (par exemple, votre codespace se met en veille), vous pouvez toujours utiliser le chemin complet pour y retourner, en supposant que vous travaillez dans l'environnement de formation GitHub Codespaces :

    ```bash
    cd /workspaces/training/nf4-science/rnaseq
    ```

Maintenant, examinons le contenu.

### Explorer le matériel fourni

Vous pouvez explorer le contenu de ce répertoire en utilisant l'explorateur de fichiers sur le côté gauche de l'espace de travail de formation.
Alternativement, vous pouvez utiliser la commande `tree`.

Tout au long du cours, nous utilisons la sortie de `tree` pour représenter la structure et le contenu des répertoires sous une forme lisible, parfois avec des modifications mineures pour plus de clarté.

Ici, nous générons une table des matières jusqu'au troisième niveau :

```bash
tree . -L 3
```

??? abstract "Contenu du répertoire"

    ```console
    .
    ├── data
    │   ├── genome.fa
    │   ├── paired-end.csv
    │   ├── reads
    │   │   ├── ENCSR000COQ1_1.fastq.gz
    │   │   ├── ENCSR000COQ1_2.fastq.gz
    │   │   ├── ENCSR000COQ2_1.fastq.gz
    │   │   ├── ENCSR000COQ2_2.fastq.gz
    │   │   ├── ENCSR000COR1_1.fastq.gz
    │   │   ├── ENCSR000COR1_2.fastq.gz
    │   │   ├── ENCSR000COR2_1.fastq.gz
    │   │   ├── ENCSR000COR2_2.fastq.gz
    │   │   ├── ENCSR000CPO1_1.fastq.gz
    │   │   ├── ENCSR000CPO1_2.fastq.gz
    │   │   ├── ENCSR000CPO2_1.fastq.gz
    │   │   └── ENCSR000CPO2_2.fastq.gz
    │   └── single-end.csv
    ├── nextflow.config
    ├── rnaseq.nf
    └── solutions
        ├── modules
        │   ├── fastqc.nf
        │   ├── fastqc_pe.nf
        │   ├── hisat2_align.nf
        │   ├── hisat2_align_pe.nf
        │   ├── multiqc.nf
        │   ├── trim_galore.nf
        │   └── trim_galore_pe.nf
        ├── rnaseq-2.1.nf
        ├── rnaseq-2.2.nf
        ├── rnaseq-2.3.nf
        ├── rnaseq-3.1.nf
        ├── rnaseq-3.2.nf
        └── rnaseq_pe-3.3.nf
    ```

Cliquez sur la boîte colorée pour développer la section et afficher son contenu.
Nous utilisons des sections repliables comme celle-ci pour afficher la sortie de commande attendue ainsi que le contenu des répertoires et des fichiers de manière concise.

- **Le fichier `rnaseq.nf`** est une ébauche de script de workflow que vous développerez au fur et à mesure du cours.

- **Le répertoire `modules`** contient des ébauches de modules de processus que vous compléterez pendant le cours.

- **Le fichier `nextflow.config`** est un fichier de configuration qui définit les propriétés minimales de l'environnement.
  Vous pouvez l'ignorer pour le moment.

- **Le répertoire `data`** contient les données d'entrée et les ressources associées, décrites plus tard dans le cours.

- **Le répertoire `solutions`** contient les scripts de workflow et les modules complets qui résultent de chaque étape du cours.
  Ils sont destinés à être utilisés comme référence pour vérifier votre travail et résoudre tout problème.
  La solution de la Partie 2 peut être utilisée comme point de départ pour la Partie 3.

## Liste de vérification de préparation

Pensez-vous être prêt·e à vous lancer ?

- [ ] Je comprends l'objectif de ce cours et ses prérequis
- [ ] Mon environnement est opérationnel
- [ ] J'ai défini mon répertoire de travail de manière appropriée

Si vous pouvez cocher toutes les cases, vous êtes prêt·e à commencer.

**Pour continuer vers la [Partie 1 : Aperçu de la méthode](./01_method.md), cliquez sur la flèche dans le coin inférieur droit de cette page.**
