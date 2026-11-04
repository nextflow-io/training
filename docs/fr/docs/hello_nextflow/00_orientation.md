# Premiers pas

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduction assistée par IA - [en savoir plus et suggérer des améliorations](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/PIjOdFaYwWA?si=gZxlXgkVxuLEzOsC" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } Voir [la playlist complète](https://youtube.com/playlist?list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&si=eF7cLR62goy-lc6n) sur la chaîne YouTube Nextflow.

:green_book: La transcription de la vidéo est disponible [ici](./transcripts/00_orientation.md).
///

!!! tip "Astuce"

    Les vidéos YouTube ont des super pouvoirs !

    - :fontawesome-solid-closed-captioning: Sous-titres de haute qualité (révisés manuellement). Activez-les avec l'icône :material-subtitles:
    - :material-bookmark: Chapitres vidéo dans la timeline qui correspondent aux titres de la page.

## Démarrer un environnement de formation

Pour utiliser l'environnement pré-construit que nous fournissons sur GitHub Codespaces, cliquez sur le bouton « Open in GitHub Codespaces » ci-dessous. Pour d'autres options, voir [Options d'environnement](../envsetup/index.md).

Nous vous recommandons d'ouvrir l'environnement de formation dans un nouvel onglet ou une nouvelle fenêtre de navigateur (utilisez le clic droit, ctrl-clic ou cmd-clic selon votre équipement) afin de pouvoir continuer à lire pendant le chargement de l'environnement.
Vous devrez garder ces instructions ouvertes en parallèle pour suivre le cours.

[![Ouvrir dans GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

### Bases de l'environnement

Cet environnement de formation contient tous les logiciels, le code et les données nécessaires pour suivre le cours de formation, vous n'avez donc rien à installer vous-même.

Le codespace est configuré avec une interface VSCode, qui comprend un explorateur de système de fichiers, un éditeur de code et un terminal shell.
Toutes les instructions données pendant le cours (par exemple « ouvrir le fichier », « modifier le code » ou « exécuter cette commande ») font référence à ces trois parties de l'interface VSCode, sauf indication contraire.

Si vous suivez ce cours par vous-même, veuillez vous familiariser avec les [bases de l'environnement](../envsetup/01_setup.md) pour plus de détails.

### Exigences de version

Cette formation est conçue pour Nextflow 25.10.2 ou ultérieur **avec le parseur de syntaxe v2 ACTIVÉ**.
Si vous utilisez un environnement local ou personnalisé, veuillez vous assurer que vous utilisez les bons paramètres comme documenté [ici](../info/nxf_versions.md).

## Se préparer à travailler

Une fois votre codespace en cours d'exécution, il y a deux choses que vous devez faire avant de plonger dans la formation : définir votre répertoire de travail pour ce cours spécifique et examiner les supports fournis.

### Définir le répertoire de travail

Par défaut, le codespace s'ouvre avec le répertoire de travail défini à la racine de tous les cours de formation, mais pour ce cours, nous travaillerons dans le répertoire `hello-nextflow/`.

Changez de répertoire maintenant en exécutant cette commande dans le terminal :

```bash
cd hello-nextflow/
```

Vous pouvez configurer VSCode pour se concentrer sur ce répertoire, de sorte que seuls les fichiers pertinents s'affichent dans la barre latérale de l'explorateur de fichiers :

```bash
code .
```

!!! tip "Astuce"

    Si pour une raison quelconque vous sortez de ce répertoire (par exemple, votre codespace se met en veille), vous pouvez toujours utiliser le chemin complet pour y revenir, en supposant que vous l'exécutez dans l'environnement de formation GitHub Codespaces :

    ```bash
    cd /workspaces/training/hello-nextflow
    ```

Examinons maintenant le contenu.

### Explorer les supports fournis

Vous pouvez explorer le contenu de ce répertoire en utilisant l'explorateur de fichiers sur le côté gauche de l'espace de travail de formation.
Sinon, vous pouvez utiliser la commande `tree`.

Tout au long du cours, nous utilisons la sortie de `tree` pour représenter la structure des répertoires et le contenu sous une forme lisible, parfois avec des modifications mineures pour plus de clarté.

Ici, nous générons une table des matières jusqu'au deuxième niveau :

```bash
tree . -L 2
```

??? abstract "Contenu du répertoire"

    ```console
    .
    ├── data
    │   └── greetings.csv
    ├── hello-channels.nf
    ├── hello-config.nf
    ├── hello-containers.nf
    ├── hello-modules.nf
    ├── hello-workflow.nf
    ├── hello-world.nf
    ├── nextflow.config
    ├── solutions
    │   ├── 1-hello-world
    │   ├── 2-hello-channels
    │   ├── 3-hello-workflow
    │   ├── 4-hello-modules
    │   ├── 5-hello-containers
    │   └── 6-hello-config
    ├── test-params.json
    └── test-params.yaml
    ```

Cliquez sur la boîte colorée pour développer la section et voir son contenu.
Nous utilisons des sections pliables comme celle-ci pour inclure la sortie de commande attendue de manière concise.

- **Les fichiers `.nf`** sont des scripts de workflow nommés en fonction de la partie du cours où ils sont utilisés.

- **Le fichier `nextflow.config`** est un fichier de configuration qui définit des propriétés d'environnement minimales.
  Vous pouvez l'ignorer pour l'instant.

- **Le fichier `greetings.csv`** sous `data/` contient les données d'entrée que nous utiliserons dans la majeure partie du cours. Il est décrit dans la Partie 2 (Channels), lorsque nous l'introduisons pour la première fois.

- **Les fichiers `test-params.*`** sont des fichiers de configuration que nous utiliserons dans la Partie 6 (Configuration). Vous pouvez les ignorer pour l'instant.

- **Le répertoire `solutions`** contient les scripts de workflow complétés qui résultent de chaque étape du cours.
  Ils sont destinés à être utilisés comme référence pour vérifier votre travail et résoudre les éventuels problèmes.

## Liste de vérification de préparation

Vous pensez être prêt·e à vous lancer ?

- [ ] Je comprends l'objectif de ce cours et ses prérequis
- [ ] Mon environnement est opérationnel
- [ ] J'ai défini mon répertoire de travail de manière appropriée

Si vous pouvez cocher toutes les cases, vous êtes prêt·e à commencer.

**Pour continuer vers la [Partie 1 : Hello World](./01_hello_world.md), cliquez sur la flèche dans le coin inférieur droit de cette page.**
