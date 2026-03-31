# Environnement de développement

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduction assistée par IA - [en savoir plus et suggérer des améliorations](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Les environnements de développement intégrés (IDE) modernes peuvent transformer radicalement votre expérience de développement Nextflow. Cette quête secondaire se concentre spécifiquement sur l'utilisation de VS Code et de son extension Nextflow pour écrire du code plus rapidement, détecter les erreurs tôt et naviguer efficacement dans des workflows complexes.

!!! note "Ce n'est pas un tutoriel traditionnel"

    Contrairement aux autres modules de formation, ce guide est organisé comme une collection d'astuces rapides, de conseils et d'exemples pratiques plutôt que comme un tutoriel pas à pas. Chaque section peut être explorée indépendamment selon vos intérêts et vos besoins de développement actuels. N'hésitez pas à naviguer librement et à vous concentrer sur les fonctionnalités qui seront les plus immédiatement utiles à votre développement de workflows.

## Ce que vous devez savoir au préalable

Ce guide suppose que vous avez suivi le cours de formation [Hello Nextflow](../hello_nextflow/) et que vous êtes à l'aise avec les concepts fondamentaux de Nextflow, notamment :

- **Structure de base d'un workflow** : Comprendre les processus, les workflows et la façon dont ils s'articulent
- **Opérations sur les canaux** : Créer des canaux, transmettre des données entre les processus et utiliser les opérateurs de base
- **Modules et organisation** : Créer des modules réutilisables et utiliser les instructions `include`
- **Bases de la configuration** : Utiliser `nextflow.config` pour les paramètres, les directives de processus et les profils

## Ce que vous apprendrez ici

Ce guide se concentre sur les **fonctionnalités de productivité de l'IDE** qui feront de vous un·e développeur·se Nextflow plus efficace :

- **Coloration syntaxique avancée** : Comprendre ce que VS Code vous indique sur la structure de votre code
- **Auto-complétion intelligente** : Exploiter les suggestions contextuelles pour écrire du code plus rapidement
- **Détection d'erreurs et diagnostics** : Détecter les erreurs de syntaxe avant d'exécuter votre workflow
- **Navigation dans le code** : Se déplacer rapidement entre les processus, les modules et les définitions
- **Formatage et organisation** : Maintenir un style de code cohérent et lisible
- **Développement assisté par IA** (optionnel) : Utiliser les outils d'IA modernes intégrés à votre IDE

!!! info "Pourquoi les fonctionnalités IDE maintenant ?"

    Vous avez probablement déjà utilisé VS Code pendant le cours [Hello Nextflow](../hello_nextflow/), mais nous avons gardé le focus sur l'apprentissage des fondamentaux de Nextflow plutôt que sur les fonctionnalités de l'IDE. Maintenant que vous êtes à l'aise avec les concepts de base de Nextflow comme les processus, les workflows, les canaux et les modules, vous êtes prêt·e à exploiter les fonctionnalités sophistiquées de l'IDE qui feront de vous un·e développeur·se plus efficace.

    Considérez cela comme une montée en niveau de votre environnement de développement — le même éditeur que vous utilisez depuis le début dispose de capacités bien plus puissantes qui deviennent vraiment précieuses une fois que vous comprenez ce qu'elles vous aident à accomplir.

---

## 0. Configuration et échauffement

Mettons en place un espace de travail spécifiquement pour explorer les fonctionnalités de l'IDE :

```bash title="Navigate to the IDE features directory"
cd side-quests/ide_features
```

Ouvrez ce répertoire dans VS Code :

```bash title="Open VS Code in current directory"
code .
```

Le répertoire `ide_features` contient des exemples de workflows qui illustrent diverses fonctionnalités de l'IDE :

```bash title="Show directory structure"
tree .
```

```console title="Project structure"
tree .
.
├── basic_workflow.nf
├── complex_workflow.nf
├── data
│   ├── sample_001.fastq.gz
│   ├── sample_002.fastq.gz
│   ├── sample_003.fastq.gz
│   ├── sample_004.fastq.gz
│   ├── sample_005.fastq.gz
│   └── sample_data.csv
├── modules
│   ├── fastqc.nf
│   ├── star.nf
│   └── utils.nf
└── nextflow.config

3 directories, 12 files
```

!!! note "À propos des fichiers d'exemple"

    - `basic_workflow.nf` est un workflow de base fonctionnel que vous pouvez exécuter et modifier
    - `complex_workflow.nf` est conçu uniquement à des fins d'illustration pour démontrer les fonctionnalités de navigation — il peut ne pas s'exécuter avec succès, mais il montre une structure de workflow multi-fichiers réaliste

### Raccourcis clavier

Certaines fonctionnalités de ce guide utilisent des raccourcis clavier optionnels. Si vous accédez à ce contenu via GitHub Codespaces dans un navigateur, il est possible que certains raccourcis ne fonctionnent pas comme prévu, car ils sont utilisés à d'autres fins dans votre système.

Si vous exécutez VS Code localement, comme vous le ferez probablement lorsque vous écrirez réellement des workflows, les raccourcis fonctionneront comme décrit.

Si vous utilisez un Mac, certains raccourcis clavier (pas tous) utiliseront "cmd" à la place de "ctrl", et nous l'indiquerons dans le texte sous la forme `Ctrl/Cmd`.

### 0.1. Installation de l'extension Nextflow

!!! note "Vous utilisez déjà des Devcontainers ?"

    Si vous travaillez dans **GitHub Codespaces** ou si vous utilisez un **devcontainer local**, l'extension Nextflow est probablement déjà installée et configurée pour vous. Vous pouvez ignorer les étapes d'installation manuelle ci-dessous et passer directement à l'exploration des fonctionnalités de l'extension.

Pour installer l'extension manuellement :

1. Ouvrez VS Code
2. Accédez à la vue Extensions en cliquant sur l'icône des extensions à gauche : ![icône des extensions](img/extensions_icon.png) (raccourci `Ctrl/Cmd+Shift+X` si vous exécutez VSCode localement)
3. Recherchez "Nextflow"
4. Installez l'extension Nextflow officielle

![Installer l'extension Nextflow](img/install_extension.png)

### 0.2. Organisation de l'espace de travail

Puisque vous avez utilisé VS Code tout au long de Hello Nextflow, vous êtes déjà familier·ère avec les bases. Voici comment organiser efficacement votre espace de travail pour cette session :

- **Zone d'édition** : Pour visualiser et modifier les fichiers. Vous pouvez la diviser en plusieurs volets pour comparer des fichiers côte à côte.
- **Explorateur de fichiers** (cliquez sur ![icône de l'explorateur de fichiers](img/files_icon.png)) (`Ctrl/Cmd+Shift+E`) : Les fichiers et dossiers locaux de votre système. Gardez-le ouvert à gauche pour naviguer entre les fichiers.
- **Terminal intégré** (`Ctrl+Shift+` backtick pour Windows et MacOS) : Un terminal pour interagir avec l'ordinateur en bas de l'écran. Utilisez-le pour exécuter Nextflow ou d'autres commandes.
- **Panneau Problèmes** (`Ctrl+Shift+M`) : VS Code affichera ici toutes les erreurs et problèmes détectés. Cela est utile pour repérer les problèmes d'un coup d'œil.

Vous pouvez déplacer les panneaux ou les masquer (`Ctrl/Cmd+B` pour basculer la barre latérale) afin de personnaliser votre disposition au fil des exemples.

### À retenir

Vous avez configuré VS Code avec l'extension Nextflow et vous comprenez l'organisation de l'espace de travail pour un développement efficace.

### Et ensuite ?

Découvrez comment la coloration syntaxique vous aide à comprendre la structure du code Nextflow en un coup d'œil.

---

## 1. Coloration syntaxique et structure du code

Maintenant que votre espace de travail est configuré, explorons comment la coloration syntaxique de VS Code vous aide à lire et à écrire du code Nextflow plus efficacement.

### 1.1. Éléments de syntaxe Nextflow

Ouvrez `basic_workflow.nf` pour voir la coloration syntaxique en action :

![Aperçu de la syntaxe](img/syntax_showcase.png)

Remarquez comment VS Code met en évidence :

- Les **mots-clés** (`process`, `workflow`, `input`, `output`, `script`) avec des couleurs distinctes
- Les **chaînes de caractères** et les **paramètres** avec un style différent
- Les **commentaires** dans une couleur atténuée
- Les **variables** et les **appels de fonctions** avec une emphase appropriée
- Les **blocs de code** avec des guides d'indentation appropriés

!!! note "Couleurs dépendantes du thème"

    Les couleurs spécifiques que vous verrez dépendront de votre thème VS Code (mode sombre/clair), de vos paramètres de couleur et des personnalisations que vous avez effectuées. L'important est que les différents éléments syntaxiques soient visuellement distingués les uns des autres, ce qui facilite la compréhension de la structure du code quel que soit le schéma de couleurs choisi.

### 1.2. Comprendre la structure du code

La coloration syntaxique vous aide à identifier rapidement :

- **Les limites des processus** : Distinction claire entre les différents processus
- **Les blocs d'entrée/sortie** : Faciles à repérer pour les définitions de flux de données
- **Les blocs de script** : Les commandes réellement exécutées
- **Les opérations sur les canaux** : Les étapes de transformation des données
- **Les directives de configuration** : Les paramètres spécifiques aux processus

Cette organisation visuelle devient indispensable lorsque vous travaillez avec des workflows complexes contenant plusieurs processus et des flux de données élaborés.

### À retenir

Vous comprenez comment la coloration syntaxique de VS Code vous aide à lire la structure du code Nextflow et à identifier les différents éléments du langage pour un développement plus rapide.

### Et ensuite ?

Découvrez comment l'auto-complétion intelligente accélère l'écriture du code grâce à des suggestions contextuelles.

---

## 2. Auto-complétion intelligente

Les fonctionnalités d'auto-complétion de VS Code vous aident à écrire du code plus rapidement et avec moins d'erreurs en suggérant des options appropriées selon le contexte.

### 2.1. Suggestions contextuelles

Les options d'auto-complétion varient selon l'endroit où vous vous trouvez dans votre code :

#### Opérations sur les canaux

Ouvrez à nouveau `basic_workflow.nf` et essayez de taper `channel.` dans le bloc workflow :

![Auto-complétion des canaux](img/autocomplete_channel.png)

Vous verrez des suggestions pour :

- `fromPath()` — Créer un canal à partir de chemins de fichiers
- `fromFilePairs()` — Créer un canal à partir de fichiers appariés
- `of()` — Créer un canal à partir de valeurs
- `fromSRA()` — Créer un canal à partir d'accessions SRA
- Et bien d'autres...

Cela vous aide à trouver rapidement la bonne fabrique de canal à utiliser sans avoir à mémoriser les noms de méthodes exacts.

Vous pouvez également découvrir les opérateurs disponibles pour les canaux. Par exemple, tapez `FASTQC.out.html.` pour voir les opérations disponibles :

![Auto-complétion des opérateurs de canaux](img/autocomplete_operators.png)

#### Directives de processus

Dans un bloc script d'un processus, tapez `task.` pour voir les propriétés d'exécution disponibles :

![Auto-complétion des propriétés de tâche](img/autocomplete_task.png)

#### Configuration

Ouvrez `nextflow.config` et tapez `process.` n'importe où pour voir les directives de processus disponibles :

![Auto-complétion de la configuration](img/autocomplete_config.png)

Vous verrez des suggestions pour :

- `executor`
- `memory`
- `cpus`

Cela fait gagner du temps lors de la configuration des processus et fonctionne dans différentes portées de configuration. Par exemple, essayez de taper `docker.` pour voir les options de configuration spécifiques à Docker.

### À retenir

Vous pouvez utiliser l'auto-complétion intelligente de VS Code pour découvrir les opérations disponibles sur les canaux, les directives de processus et les options de configuration sans mémoriser la syntaxe.

### Et ensuite ?

Découvrez comment la détection d'erreurs en temps réel vous aide à repérer les problèmes avant d'exécuter votre workflow, simplement en lisant le code.

## 3. Détection d'erreurs et diagnostics

La détection d'erreurs en temps réel de VS Code vous aide à repérer les problèmes avant d'exécuter votre workflow.

### 3.1. Détection des erreurs de syntaxe

Créons une erreur délibérée pour voir la détection en action. Ouvrez `basic_workflow.nf` et changez le nom du processus de `FASTQC` en `FASTQ` (ou tout autre nom invalide). VS Code mettra immédiatement en évidence l'erreur dans le bloc workflow avec un soulignement rouge ondulé :

![Soulignement d'erreur](img/error_underline.png)

### 3.2. Panneau Problèmes

Au-delà de la mise en évidence individuelle des erreurs, VS Code fournit un panneau Problèmes centralisé qui regroupe toutes les erreurs, avertissements et messages d'information dans votre espace de travail. Ouvrez-le avec `Ctrl/Cmd+Shift+M` et utilisez l'icône de filtre pour n'afficher que les erreurs relatives au fichier courant :

![Filtrer le panneau Problèmes](img/active_file.png)

Cliquez sur n'importe quel problème pour accéder directement à la ligne concernée.

![Panneau Problèmes](img/problems_panel.png)

Corrigez l'erreur en remettant le nom du processus à `FASTQC`.

### 3.3. Erreurs courantes

Les erreurs courantes dans la syntaxe Nextflow incluent :

- **Accolades manquantes** : `{` ou `}` non appariées
- **Blocs incomplets** : Sections requises manquantes dans les processus
- **Syntaxe invalide** : DSL Nextflow malformé
- **Fautes de frappe dans les mots-clés** : Directives de processus mal orthographiées
- **Incompatibilités de canaux** : Incompatibilités de types

Le serveur de langage Nextflow met en évidence ces problèmes dans le panneau Problèmes. Vous pouvez les vérifier tôt pour éviter les erreurs de syntaxe lors de l'exécution d'un pipeline.

### À retenir

Vous pouvez utiliser la détection d'erreurs de VS Code et le panneau Problèmes pour repérer les erreurs de syntaxe et les problèmes avant d'exécuter votre workflow, ce qui vous fait gagner du temps et évite les frustrations.

### Et ensuite ?

Découvrez comment naviguer efficacement entre les processus, les modules et les définitions dans des workflows complexes.

---

## 4. Navigation dans le code et gestion des symboles

Une navigation efficace est essentielle lorsque vous travaillez avec des workflows complexes répartis sur plusieurs fichiers. Pour comprendre cela, remplacez la définition du processus dans `basic_workflow.nf` par un import du module que nous vous avons fourni :

=== "Après"

    ```groovy title="basic_workflow.nf" linenums="3"
    include { FASTQC } from './modules/fastqc.nf'
    ```

=== "Avant"

    ```groovy title="basic_workflow.nf" linenums="3"
    process FASTQC {
        tag "${sample_id}"
        publishDir "${params.output_dir}/fastqc", mode: 'copy'

        input:
        tuple val(sample_id), path(reads)

        output:
        tuple val(sample_id), path("*.html"), emit: html
        tuple val(sample_id), path("*.zip"), emit: zip

        script:
        def args = task.ext.args ?: ''
        """
        fastqc \\
            ${args} \\
            --threads ${task.cpus} \\
            ${reads}
        """
    }
    ```

### 4.1. Aller à la définition

Si vous survolez un nom de processus comme `FASTQC`, vous verrez une fenêtre contextuelle avec l'interface du module (entrées et sorties) :

![Aller à la définition](img/syntax.png)

Cette fonctionnalité est particulièrement précieuse lors de la création de workflows, car elle vous permet de comprendre l'interface du module sans ouvrir directement le fichier du module.

Vous pouvez naviguer rapidement vers n'importe quelle définition de processus, de module ou de variable en utilisant **Ctrl/Cmd+clic**. Survolez le lien vers le fichier du module en haut du script et suivez le lien comme suggéré :

![Suivre le lien](img/follow_link.png)

La même chose fonctionne pour les noms de processus. Revenez à `basic_workflow.nf` et essayez cela sur le nom du processus `FASTQC` dans le bloc workflow. Cela vous amène directement au nom du processus (qui est le même que le fichier du module dans cet exemple, mais pourrait se trouver au milieu d'un fichier beaucoup plus grand).

Pour revenir à l'endroit où vous étiez, utilisez **Alt+←** (ou **Ctrl+-** sur Mac). C'est un moyen puissant d'explorer le code sans perdre votre position.

Explorons maintenant la navigation dans un workflow plus complexe en utilisant `complex_workflow.nf` (le fichier d'illustration mentionné précédemment). Ce workflow contient plusieurs processus définis dans des fichiers de modules séparés, ainsi que quelques processus en ligne. Bien que les structures multi-fichiers complexes puissent être difficiles à naviguer manuellement, la possibilité de sauter aux définitions rend l'exploration beaucoup plus gérable.

1. Ouvrez `complex_workflow.nf`
2. Naviguez vers les définitions de modules
3. Utilisez **Alt+←** (ou **Ctrl+-**) pour revenir en arrière
4. Naviguez vers le nom du processus `FASTQC` dans le bloc workflow. Cela vous amène directement au nom du processus (qui est le même que le fichier du module dans cet exemple, mais pourrait se trouver au milieu d'un fichier beaucoup plus grand).
5. Revenez en arrière
6. Naviguez vers le processus `TRIM_GALORE` dans le bloc workflow. Celui-ci est défini en ligne, donc il ne vous amènera pas à un fichier séparé, mais il vous montrera quand même la définition du processus, et vous pourrez toujours revenir à l'endroit où vous étiez.

### 4.2. Navigation par symboles

Avec `complex_workflow.nf` toujours ouvert, vous pouvez obtenir un aperçu de tous les symboles du fichier en tapant `@` dans la barre de recherche en haut de VSCode (le raccourci clavier est `Ctrl/Cmd+Shift+O`, mais il peut ne pas fonctionner dans Codespaces). Cela ouvre le panneau de navigation par symboles, qui liste tous les symboles du fichier courant :

![Navigation par symboles](img/symbols.png)

Cela affiche :

- Toutes les définitions de processus
- Les définitions de workflow (il y a deux workflows définis dans ce fichier)
- Les définitions de fonctions

Commencez à taper pour filtrer les résultats.

### 4.3. Trouver toutes les références

Comprendre où un processus ou une variable est utilisé dans votre base de code peut être très utile. Par exemple, si vous souhaitez trouver toutes les références au processus `FASTQC`, commencez par naviguer vers sa définition. Vous pouvez le faire en ouvrant directement `modules/fastqc.nf`, ou en utilisant la fonctionnalité de navigation rapide de VS Code avec `Ctrl/Cmd+clic` comme nous l'avons fait précédemment. Une fois à la définition du processus, faites un clic droit sur le nom du processus `FASTQC` et sélectionnez "Find All References" dans le menu contextuel pour voir toutes les instances où il est utilisé.

![Trouver les références](img/references.png)

Cette fonctionnalité affiche toutes les instances où `FASTQC` est référencé dans votre espace de travail, y compris son utilisation dans les deux workflows distincts. Cette information est cruciale pour évaluer l'impact potentiel des modifications apportées au processus `FASTQC`.

### 4.4. Panneau Plan

Le panneau Plan, situé dans la barre latérale de l'Explorateur (cliquez sur ![icône de l'Explorateur](img/files_icon.png)), fournit un aperçu pratique de tous les symboles de votre fichier courant. Cette fonctionnalité vous permet de naviguer rapidement et de gérer la structure de votre code en affichant les fonctions, les variables et d'autres éléments clés dans une vue hiérarchique.

![Panneau Plan](img/outline.png)

Utilisez le panneau Plan pour naviguer rapidement vers différentes parties de votre code sans utiliser le navigateur de fichiers.

### 4.5. Visualisation du DAG

L'extension Nextflow de VS Code peut visualiser votre workflow sous forme de graphe acyclique dirigé (DAG). Cela vous aide à comprendre le flux de données et les dépendances entre les processus. Ouvrez `complex_workflow.nf` et cliquez sur le bouton "Preview DAG" au-dessus de `workflow {` (le deuxième bloc `workflow` dans ce fichier) :

![Aperçu du DAG](img/dag_preview.png)

Il s'agit uniquement du workflow d'entrée, mais vous pouvez également prévisualiser le DAG pour les workflows internes en cliquant sur le bouton "Preview DAG" au-dessus du workflow `RNASEQ_PIPELINE {` plus haut :

![Aperçu du DAG pour le workflow interne](img/dag_preview_inner.png)

Pour ce workflow, vous pouvez utiliser les nœuds du DAG pour naviguer vers les définitions de processus correspondantes dans le code. Cliquez sur un nœud et il vous amènera à la définition du processus correspondant dans l'éditeur. En particulier lorsqu'un workflow devient très grand, cela peut vraiment vous aider à naviguer dans le code et à comprendre comment les processus sont connectés.

### À retenir

Vous pouvez naviguer efficacement dans des workflows complexes en utilisant aller à la définition, la recherche de symboles, trouver les références et la visualisation du DAG pour comprendre la structure du code et les dépendances.

### Et ensuite ?

Découvrez comment travailler efficacement sur plusieurs fichiers interconnectés dans des projets Nextflow plus importants.

## 5. Travailler sur plusieurs fichiers

Le développement Nextflow réel implique de travailler avec plusieurs fichiers interconnectés. Explorons comment VS Code vous aide à gérer efficacement des projets complexes.

### 5.1. Navigation rapide entre les fichiers

Avec `complex_workflow.nf` ouvert, vous remarquerez qu'il importe plusieurs modules. Pratiquons la navigation rapide entre eux.

Appuyez sur **Ctrl+P** (ou **Cmd+P**) et commencez à taper "fast" :

VS Code vous montrera les fichiers correspondants. Sélectionnez `modules/fastqc.nf` pour y accéder instantanément. C'est bien plus rapide que de cliquer dans l'explorateur de fichiers lorsque vous savez approximativement quel fichier vous cherchez.

Essayez cela avec d'autres motifs :

- Tapez "star" pour trouver le fichier du module d'alignement STAR (`star.nf`)
- Tapez "utils" pour trouver le fichier des fonctions utilitaires (`utils.nf`)
- Tapez "config" pour accéder aux fichiers de configuration (`nextflow.config`)

### 5.2. Éditeur divisé pour le développement multi-fichiers

Lorsque vous travaillez avec des modules, vous avez souvent besoin de voir simultanément le workflow principal et les définitions des modules. Configurons cela :

1. Ouvrez `complex_workflow.nf`
2. Ouvrez `modules/fastqc.nf` dans un nouvel onglet
3. Faites un clic droit sur l'onglet `modules/fastqc.nf` et sélectionnez "Split Right"
4. Vous pouvez maintenant voir les deux fichiers côte à côte

![Éditeur divisé](img/split_editor.png)

Cela est indispensable lorsque :

- Vous vérifiez les interfaces des modules lors de l'écriture des appels de workflow, et l'aperçu ne suffit pas
- Vous comparez des processus similaires dans différents modules
- Vous déboguez le flux de données entre le workflow et les modules

### 5.3. Recherche dans tout le projet

Parfois, vous devez trouver où des motifs spécifiques sont utilisés dans l'ensemble de votre projet. Appuyez sur `Ctrl/Cmd+Shift+F` pour ouvrir le panneau de recherche.

Essayez de rechercher `publishDir` dans tout l'espace de travail :

![Recherche dans le projet](img/project_search.png)

Cela vous montre chaque fichier qui utilise des répertoires de publication, vous aidant à :

- Comprendre les schémas d'organisation des sorties
- Trouver des exemples de directives spécifiques
- Assurer la cohérence entre les modules

### À retenir

Vous pouvez gérer des projets multi-fichiers complexes en utilisant la navigation rapide entre les fichiers, les éditeurs divisés et la recherche dans tout le projet pour travailler efficacement sur les workflows et les modules.

### Et ensuite ?

Découvrez comment les fonctionnalités de formatage et de maintenance du code maintiennent vos workflows organisés et lisibles.

---

## 6. Formatage et maintenance du code

Un formatage correct du code est essentiel non seulement pour l'esthétique, mais aussi pour améliorer la lisibilité, la compréhension et la facilité de mise à jour des workflows complexes.

### 6.1. Formatage automatique en action

Ouvrez `basic_workflow.nf` et dégradez délibérément le formatage :

- Supprimez une partie de l'indentation : Sélectionnez l'intégralité du document et appuyez plusieurs fois sur `Shift+Tab` pour supprimer autant d'indentations que possible.
- Ajoutez des espaces supplémentaires à des endroits aléatoires : dans l'instruction `channel.fromPath`, ajoutez 30 espaces après le `(`.
- Cassez certaines lignes de manière maladroite : Ajoutez une nouvelle ligne entre l'opérateur `.view {` et la chaîne `Processing sample:`, mais n'ajoutez pas de nouvelle ligne correspondante avant la parenthèse fermante `}`.

Appuyez maintenant sur `Shift+Alt+F` (ou `Shift+Option+F` sur MacOS) pour formater automatiquement :

VS Code effectue immédiatement :

- La correction de l'indentation pour montrer clairement la structure des processus
- L'alignement cohérent des éléments similaires
- La suppression des espaces blancs inutiles
- Le maintien de sauts de ligne lisibles

Notez que le formatage automatique peut ne pas résoudre tous les problèmes de style de code. Le serveur de langage Nextflow vise à garder votre code ordonné, mais il respecte également vos préférences personnelles dans certains domaines. Par exemple, si vous supprimez l'indentation à l'intérieur du bloc `script` d'un processus, le formateur le laissera tel quel, car vous pourriez intentionnellement préférer ce style.

Actuellement, il n'existe pas d'application stricte du style pour Nextflow, donc le serveur de langage offre une certaine flexibilité. Cependant, il appliquera systématiquement des règles de formatage autour des définitions de méthodes et de fonctions pour maintenir la clarté.

### 6.2. Fonctionnalités d'organisation du code

#### Commentaire rapide

Sélectionnez un bloc de code dans votre workflow et appuyez sur **Ctrl+/** (ou **Cmd+/**) pour le commenter :

```groovy
// workflow {
//     ch_input = channel.fromPath(params.input)
//         .splitCsv(header: true)
//         .map { row -> [row.sample_id, file(row.fastq_path)] }
//
//     FASTQC(ch_input)
// }
```

C'est parfait pour :

- Désactiver temporairement des parties de workflows pendant le développement
- Ajouter des commentaires explicatifs à des opérations de canaux complexes
- Documenter des sections de workflow

Utilisez à nouveau **Ctrl+/** (ou **Cmd+/**) pour décommenter le code.

#### Repliement de code pour une vue d'ensemble

Dans `complex_workflow.nf`, remarquez les petites flèches à côté des définitions de processus. Cliquez dessus pour replier (réduire) les processus :

![Repliement de code](img/code_folding.png)

Cela vous donne une vue d'ensemble de haut niveau de la structure de votre workflow sans vous perdre dans les détails d'implémentation.

#### Correspondance des accolades

Placez votre curseur à côté d'une accolade `{` ou `}` et VS Code met en évidence l'accolade correspondante. Utilisez **Ctrl+Shift+\\** (ou **Cmd+Shift+\\**) pour sauter entre les accolades correspondantes.

C'est essentiel pour :

- Comprendre les limites des processus
- Trouver les accolades manquantes ou en trop
- Naviguer dans les structures de workflow imbriquées

#### Sélection et édition multi-lignes

Pour modifier plusieurs lignes simultanément, VS Code offre de puissantes capacités multi-curseurs :

- **Sélection multi-lignes** : Maintenez **Ctrl+Alt** (ou **Cmd+Option** pour MacOS) et utilisez les touches fléchées pour sélectionner plusieurs lignes
- **Indentation multi-lignes** : Sélectionnez plusieurs lignes et utilisez **Tab** pour indenter ou **Shift+Tab** pour désindenter des blocs entiers

C'est particulièrement utile pour :

- Indenter des blocs de processus entiers de manière cohérente
- Ajouter des commentaires à plusieurs lignes à la fois
- Modifier des définitions de paramètres similaires dans plusieurs processus

### À retenir

Vous pouvez maintenir un code propre et lisible en utilisant le formatage automatique, les fonctionnalités de commentaire, le repliement de code, la correspondance des accolades et l'édition multi-lignes pour organiser efficacement des workflows complexes.

### Et ensuite ?

Découvrez comment VS Code s'intègre à votre workflow de développement plus large, au-delà de la simple édition de code.

---

## 7. Intégration dans le workflow de développement

VS Code s'intègre bien à votre workflow de développement au-delà de la simple édition de code.

### 7.1. Intégration du contrôle de version

!!! note "Codespaces et intégration Git"

    Si vous travaillez dans **GitHub Codespaces**, certaines fonctionnalités d'intégration Git peuvent ne pas fonctionner comme prévu, notamment les raccourcis clavier pour le contrôle de source. Vous avez peut-être également refusé d'ouvrir le répertoire en tant que dépôt Git lors de la configuration initiale, ce qui est acceptable à des fins de formation.

Si votre projet est un dépôt git (comme c'est le cas ici), VS Code affiche :

- Les fichiers modifiés avec des indicateurs colorés
- Le statut Git dans la barre d'état
- Les vues de différences en ligne
- Les capacités de commit et de push

Ouvrez le panneau Contrôle de source en utilisant le bouton de contrôle de source (![icône du contrôle de source](img/source_control_icon.png)) (`Ctrl+Shift+G` ou `Cmd+Shift+G` si vous travaillez avec VSCode localement) pour voir les modifications git et effectuer des commits directement dans l'éditeur.

![Panneau Contrôle de source](img/source_control.png)

### 7.2. Exécution et inspection des workflows

Exécutons un workflow puis inspectons les résultats. Dans le terminal intégré (`Ctrl+Shift+` backtick pour Windows et MacOS), exécutez le workflow de base :

```bash title="Run the basic workflow"
nextflow run basic_workflow.nf --input data/sample_data.csv --output_dir results
```

Pendant l'exécution du workflow, vous verrez la sortie en temps réel dans le terminal. Après la fin de l'exécution, vous pouvez utiliser VS Code pour inspecter les résultats sans quitter votre éditeur :

1. **Naviguer vers les répertoires de travail** : Utilisez l'explorateur de fichiers ou le terminal pour parcourir `.nextflow/work`
2. **Ouvrir les fichiers journaux** : Cliquez sur les chemins des fichiers journaux dans la sortie du terminal pour les ouvrir directement dans VS Code
3. **Inspecter les sorties** : Parcourez les répertoires de résultats publiés dans l'explorateur de fichiers
4. **Visualiser les rapports d'exécution** : Ouvrez les rapports HTML directement dans VS Code ou dans votre navigateur

Cela permet de tout garder au même endroit plutôt que de basculer entre plusieurs applications.

### À retenir

Vous pouvez intégrer VS Code au contrôle de version et à l'exécution des workflows pour gérer l'ensemble de votre processus de développement depuis une seule interface.

### Et ensuite ?

Voyez comment toutes ces fonctionnalités IDE fonctionnent ensemble dans votre workflow de développement quotidien.

---

## 8. Récapitulatif et notes rapides

Voici quelques notes rapides sur chacune des fonctionnalités IDE abordées ci-dessus :

### 8.1. Démarrer une nouvelle fonctionnalité

1. **Ouverture rapide de fichiers** (`Ctrl+P` ou `Cmd+P`) pour trouver les modules existants pertinents
2. **Éditeur divisé** pour visualiser des processus similaires côte à côte
3. **Navigation par symboles** (`Ctrl+Shift+O` ou `Cmd+Shift+O`) pour comprendre la structure du fichier
4. **Auto-complétion** pour écrire rapidement du nouveau code

### 8.2. Déboguer des problèmes

1. **Panneau Problèmes** (`Ctrl+Shift+M` ou `Cmd+Shift+M`) pour voir toutes les erreurs en même temps
2. **Aller à la définition** (`Ctrl+clic` ou `Cmd+clic`) pour comprendre les interfaces des processus
3. **Trouver toutes les références** pour voir comment les processus sont utilisés
4. **Recherche dans tout le projet** pour trouver des motifs ou des problèmes similaires

### 8.3. Refactorisation et amélioration

1. **Recherche dans tout le projet** (`Ctrl+Shift+F` ou `Cmd+Shift+F`) pour trouver des motifs
2. **Formatage automatique** (`Shift+Alt+F` ou `Shift+Option+F`) pour maintenir la cohérence
3. **Repliement de code** pour se concentrer sur la structure
4. **Intégration Git** pour suivre les modifications

---

## Résumé

Vous venez de faire un tour rapide des fonctionnalités IDE de VS Code pour le développement Nextflow. Ces outils vous rendront considérablement plus productif·ve en :

- **Réduisant les erreurs** grâce à la vérification syntaxique en temps réel
- **Accélérant le développement** avec l'auto-complétion intelligente
- **Améliorant la navigation** dans des workflows complexes multi-fichiers
- **Maintenant la qualité** grâce à un formatage cohérent
- **Améliorant la compréhension** grâce à la coloration avancée et à la visualisation de la structure

Nous ne nous attendons pas à ce que vous vous souveniez de tout, mais maintenant que vous savez que ces fonctionnalités existent, vous serez en mesure de les retrouver quand vous en aurez besoin. Au fur et à mesure que vous continuerez à développer des workflows Nextflow, ces fonctionnalités IDE deviendront une seconde nature, vous permettant de vous concentrer sur l'écriture de code de qualité plutôt que de vous battre avec la syntaxe et la structure.

### Et ensuite ?

Appliquez ces compétences IDE en travaillant sur d'autres modules de formation, par exemple :

- **[nf-test](nf-test.md)** : Créez des suites de tests complètes pour vos workflows
- **[Hello nf-core](../../hello_nf-core/)** : Construisez des pipelines de qualité production avec les standards de la communauté

La véritable puissance de ces fonctionnalités IDE se révèle lorsque vous travaillez sur des projets plus grands et plus complexes. Commencez à les intégrer progressivement dans votre workflow — en quelques sessions, elles deviendront une seconde nature et transformeront votre approche du développement Nextflow.

De la détection des erreurs avant qu'elles ne vous ralentissent à la navigation aisée dans des bases de code complexes, ces outils feront de vous un·e développeur·se plus confiant·e et plus efficace.

Bon codage !
