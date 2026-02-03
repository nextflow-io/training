# Premiers pas

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduction assistée par IA - [en savoir plus et suggérer des améliorations](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

## Démarrer un environnement de formation

Pour utiliser l'environnement préconstruit que nous fournissons sur GitHub Codespaces, cliquez sur le bouton « Open in GitHub Codespaces » ci-dessous. Pour d'autres options, consultez [Options d'environnement](../envsetup/index.md).

Nous vous recommandons d'ouvrir l'environnement de formation dans un nouvel onglet ou une nouvelle fenêtre de navigateur (utilisez le clic droit, ctrl+clic ou cmd+clic selon votre équipement) afin de pouvoir continuer à lire pendant le chargement de l'environnement.
Vous devrez garder ces instructions ouvertes en parallèle pour suivre le cours.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

### Notions de base de l'environnement

Cet environnement de formation contient tous les logiciels, le code et les données nécessaires pour suivre le cours de formation, vous n'avez donc rien à installer vous-même.

Le codespace est configuré avec une interface VSCode, qui comprend un explorateur de système de fichiers, un éditeur de code et un terminal shell.
Toutes les instructions données pendant le cours (par exemple « ouvrir le fichier », « modifier le code » ou « exécuter cette commande ») font référence à ces trois parties de l'interface VSCode, sauf indication contraire.

Si vous suivez ce cours par vous-même, veuillez vous familiariser avec les [notions de base de l'environnement](../envsetup/01_setup.md) pour plus de détails.

### Exigences de version

Cette formation est conçue pour **Nextflow 25.10.2** ou version ultérieure **avec l'analyseur de syntaxe v2 DÉSACTIVÉ**.

#### Si vous utilisez notre environnement de formation :

Vous DEVEZ exécuter la commande suivante avant d'aller plus loin :

```bash
export NXF_SYNTAX_PARSER=v1
```

#### Si vous utilisez un environnement local ou personnalisé :

Veuillez vous assurer d'utiliser les paramètres corrects comme documenté [ici](../info/nxf_versions.md).

La formation nécessite en outre **nf-core tools 3.4.1**.
Si vous utilisez une version différente des outils nf-core, vous pourriez rencontrer des difficultés à suivre.

Vous pouvez vérifier quelle version est installée dans votre environnement à l'aide de la commande `nf-core --version`.

## Préparez-vous à travailler

Une fois votre codespace en cours d'exécution, vous devez faire deux choses avant de vous plonger dans la formation : définir votre répertoire de travail pour ce cours spécifique et examiner les ressources fournies.

### Définir le répertoire de travail

Par défaut, le codespace s'ouvre avec le répertoire de travail défini à la racine de tous les cours de formation, mais pour ce cours, nous travaillerons dans le répertoire `hello-nf-core/`.

Changez de répertoire maintenant en exécutant cette commande dans le terminal :

```bash
cd hello-nf-core/
```

!!! tip "Astuce"

    Si pour une raison quelconque vous sortez de ce répertoire (par exemple, votre codespace se met en veille), vous pouvez toujours utiliser le chemin complet pour y revenir, en supposant que vous exécutez ceci dans l'environnement de formation Github Codespaces :

    ```bash
    cd /workspaces/training/hello-nf-core
    ```

Maintenant, jetons un coup d'œil au contenu de ce répertoire.

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
    ├── greetings.csv
    ├── original-hello
    │   ├── hello.nf
    │   ├── modules
    │   └── nextflow.config
    └── solutions
        ├── composable-hello
        ├── core-hello-part2
        ├── core-hello-part3
        ├── core-hello-part4
        ├── core-hello-part5
        └── core-hello-start
    ```

Cliquez sur la case colorée pour développer la section et afficher son contenu.
Nous utilisons des sections repliables comme celle-ci pour inclure la sortie de commande attendue de manière concise.

- **Le fichier `greetings.csv`** est un CSV contenant des données colonnes minimales que nous utilisons à des fins de test.

- **Le répertoire `original-hello`** contient une copie du code source produit en suivant la série de formation complète Hello Nextflow (avec Docker activé).

- **Le répertoire `solutions`** contient les scripts de workflow complétés qui résultent de chaque étape du cours.
  Ils sont destinés à être utilisés comme référence pour vérifier votre travail et résoudre les problèmes éventuels.

## Liste de vérification de préparation

Pensez-vous être prêt·e à vous lancer ?

- [ ] Je comprends l'objectif de ce cours et ses prérequis
- [ ] Mon environnement est opérationnel
- [ ] Je me suis assuré que l'analyseur de syntaxe est défini sur **v1**
- [ ] J'ai défini mon répertoire de travail de manière appropriée

Si vous pouvez cocher toutes les cases, vous êtes prêt·e à commencer.

**Pour continuer vers la Partie 1, cliquez sur la flèche dans le coin inférieur droit de cette page.**
